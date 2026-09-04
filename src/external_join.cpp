/*
  Copyright (c) 2026 Jouni Siren and GCSA2 contributors

  External-memory prefix-doubling join. The implementation favors a strict
  resident-memory bound and restart-friendly immutable files over throughput.
*/

#include <gcsa/path_graph.h>
#include <gcsa/compressed_block.h>
#include <gcsa/resources.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <csignal>
#include <cstdint>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <queue>
#include <spawn.h>
#include <stdexcept>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>

extern char** environ;

namespace gcsa
{

namespace
{

constexpr std::uint64_t JOIN_HEADER_MAGIC = 0x314e494f4a534347ULL; // "GCSJOIN1"
constexpr std::uint64_t JOIN_FOOTER_MAGIC = 0x31444e454a534347ULL; // "GCSJEND1"
constexpr std::uint32_t JOIN_FORMAT_VERSION = 1;
// Every sorted join run has an immutable companion that carries exact key
// groups and a one-byte per-record aggregate for range tiling.  Planning can
// consequently stream compact summaries instead of reopening every wide
// JoinRecord after the bounded radix sample.
constexpr std::uint64_t JOIN_GROUP_HEADER_MAGIC = 0x3150555247534347ULL; // "GCSGRUP1"
constexpr std::uint64_t JOIN_GROUP_FOOTER_MAGIC = 0x31444e4547534347ULL; // "GCSGEND1"
constexpr std::uint64_t JOIN_DETAIL_HEADER_MAGIC = 0x314c544447534347ULL; // "GCSGDTL1"
constexpr std::uint64_t JOIN_DETAIL_FOOTER_MAGIC = 0x3154454447534347ULL; // "GCSGDET1"
constexpr std::uint32_t JOIN_GROUP_FORMAT_VERSION = 1;
constexpr size_type JOIN_BASE_FIXED_BYTES = 32 * KILOBYTE;
constexpr size_type JOIN_LABEL_COUNT = PathLabel::LABEL_LENGTH + 1;
constexpr size_type JOIN_HEADER_BYTES = 8 + 4 + 4 + 4 + 8;
constexpr size_type JOIN_FOOTER_BYTES = 8 + 8 + 8;
constexpr size_type JOIN_RECORD_BYTES =
  4 + 3 * sizeof(node_type) + 4 + 8 + JOIN_LABEL_COUNT * sizeof(PathNode::rank_type);
constexpr size_type JOIN_GROUP_HEADER_BYTES = 8 + 4 + 4 + 4 + 4 + 3 * sizeof(std::uint64_t);
constexpr size_type JOIN_GROUP_FOOTER_BYTES = 8 + 4 * sizeof(std::uint64_t);
constexpr size_type JOIN_DETAIL_HEADER_BYTES = 8 + 4 + 4 + 4 + 4 + 2 * sizeof(std::uint64_t);
constexpr size_type JOIN_DETAIL_FOOTER_BYTES = 8 + 3 * sizeof(std::uint64_t);
constexpr size_type JOIN_GROUP_SUMMARY_BYTES = sizeof(node_type) + 5 * sizeof(std::uint64_t);
constexpr std::uint8_t JOIN_DETAIL_ORDER_MASK = 0x1F;
constexpr std::uint8_t JOIN_DETAIL_SORTED = 0x80;
constexpr std::uint8_t JOIN_DETAIL_RESERVED = 0x60;
// Join records are compact on disk, so issuing one system call per record can
// turn a sequential spill into syscall-bound I/O. Keep all phase-local I/O
// caches byte-sized and account for them in the same join reservation.
constexpr size_type JOIN_IO_BUFFER_BYTES = 64 * KILOBYTE;
constexpr size_type JOIN_IO_BUFFER_RECORDS = JOIN_IO_BUFFER_BYTES / JOIN_RECORD_BYTES;
constexpr size_type JOIN_PARALLEL_SORT_MIN_RECORDS = 64 * 1024;
// Run generation may simultaneously have a join writer, path/rank source
// caches, and path/rank output caches. Merge reader caches are accounted per
// input when selecting fan-in below.
// The output run has a payload cache plus summary and compact-detail sidecar
// caches. Account for all of them in the same reservation as merge readers.
constexpr size_type JOIN_FIXED_BYTES = JOIN_BASE_FIXED_BYTES + 7 * JOIN_IO_BUFFER_BYTES;
static_assert(JOIN_IO_BUFFER_RECORDS > 0, "join record exceeds the I/O buffer");
static_assert(PathLabel::LABEL_LENGTH <= JOIN_DETAIL_ORDER_MASK,
  "join group detail cannot encode PathNode order");
// Clean page cache is charged to a cgroup's MemoryMax. Keep a small rolling
// tail and evict completed sequential prefixes so disk, not cache, remains the
// authoritative working set. The larger flush interval amortizes fdatasync().
constexpr off_t JOIN_CACHE_TAIL_BYTES = 64 * MEGABYTE;
constexpr off_t JOIN_CACHE_FLUSH_BYTES = 512 * MEGABYTE;

enum JoinKeyKind : std::uint32_t
{
  LEFT_BY_TO = 1,
  RIGHT_BY_FROM = 2
};

struct JoinRecord
{
  logical_file_id_t logical;
  node_type key;
  PathNode node;
  std::uint64_t ordinal;
  PathNode::rank_type labels[JOIN_LABEL_COUNT];
};

struct JoinGroupSummary
{
  node_type key;
  size_type begin, end, count, order_sum;
  size_type bypass_paths, bypass_ranks;

  JoinGroupSummary() : key(0), begin(0), end(0), count(0), order_sum(0),
    bypass_paths(0), bypass_ranks(0) { }
};

struct JoinGroupSidecar
{
  std::string summary_name, detail_name;
  size_type groups;
  std::uint64_t summary_checksum, detail_checksum;

  JoinGroupSidecar(const std::string& summaries = std::string(),
    const std::string& details = std::string(), size_type group_count = 0,
    std::uint64_t summaries_checksum = 0, std::uint64_t details_checksum = 0) :
    summary_name(summaries), detail_name(details), groups(group_count),
    summary_checksum(summaries_checksum), detail_checksum(details_checksum) { }

  bool complete() const
  {
    return !this->summary_name.empty() && !this->detail_name.empty();
  }
};

struct JoinRun
{
  std::string name;
  size_type records;
  std::uint64_t checksum;
  JoinGroupSidecar sidecar;

  JoinRun(const std::string& path = std::string(), size_type count = 0,
    std::uint64_t payload_checksum = 0,
    const JoinGroupSidecar& group_sidecar = JoinGroupSidecar()) :
    name(path), records(count), checksum(payload_checksum), sidecar(group_sidecar) { }
};

void
removeJoinRun(JoinRun& run)
{
  TempFile::remove(run.name);
  TempFile::remove(run.sidecar.summary_name);
  TempFile::remove(run.sidecar.detail_name);
}

std::runtime_error
joinError(const std::string& message, const std::string& path = std::string())
{
  return std::runtime_error("externalPathGraphExtend(): " + message +
    (path.empty() ? std::string() : ": " + path));
}

size_type checkedJoinAdd(size_type left, size_type right, const char* description);
size_type checkedJoinMultiply(size_type left, size_type right, const char* description);
size_type storedJoinBytes(const std::string& path);

void
adviseSequential(int descriptor)
{
#if defined(POSIX_FADV_SEQUENTIAL)
  // Advice is an optimization. Correctness and durability do not depend on a
  // filesystem implementing it, so unsupported advice is deliberately ignored.
  static_cast<void>(::posix_fadvise(descriptor, 0, 0, POSIX_FADV_SEQUENTIAL));
#else
  static_cast<void>(descriptor);
#endif
}

void
discardCachedRange(int descriptor, off_t offset, off_t bytes)
{
#if defined(POSIX_FADV_DONTNEED)
  if(bytes >= 0)
  {
    static_cast<void>(::posix_fadvise(descriptor, offset, bytes, POSIX_FADV_DONTNEED));
  }
#else
  static_cast<void>(descriptor); static_cast<void>(offset); static_cast<void>(bytes);
#endif
}

void
trimWrittenCache(int descriptor, off_t& released, bool complete,
  const std::string& path)
{
  off_t end = ::lseek(descriptor, 0, SEEK_CUR);
  if(end < 0) { throw joinError("cannot determine output position", path); }
  if(!complete && end - released < JOIN_CACHE_FLUSH_BYTES) { return; }
  if(::fdatasync(descriptor) != 0) { throw joinError("cannot sync sequential output", path); }
  off_t discard_end = (complete ? end : std::max(released, end - JOIN_CACHE_TAIL_BYTES));
  if(discard_end > released)
  {
    discardCachedRange(descriptor, released, discard_end - released);
    released = discard_end;
  }
}

void
trimReadCache(int descriptor, off_t consumed, off_t& released)
{
  // A pathological key may replay a range. Resetting the watermark on a large
  // backwards seek ensures each replay still has a bounded cache footprint.
  if(consumed + JOIN_CACHE_FLUSH_BYTES < released) { released = 0; }
  if(consumed - released < JOIN_CACHE_FLUSH_BYTES) { return; }
  off_t discard_end = std::max(released, consumed - JOIN_CACHE_TAIL_BYTES);
  discardCachedRange(descriptor, released, discard_end - released);
  released = discard_end;
}

void
writeAll(int descriptor, const void* source, size_type bytes, const std::string& path)
{
  const std::uint8_t* data = static_cast<const std::uint8_t*>(source);
  while(bytes > 0)
  {
    ssize_t written = ::write(descriptor, data, bytes);
    if(written < 0 && errno == EINTR) { continue; }
    if(written <= 0) { throw joinError("write failed", path); }
    DiskIO::write_volume += static_cast<size_type>(written);
    data += written; bytes -= written;
  }
}

void
pwriteAll(int descriptor, const void* source, size_type bytes, off_t offset,
  const std::string& path)
{
  const std::uint8_t* data = static_cast<const std::uint8_t*>(source);
  while(bytes > 0)
  {
    ssize_t written = ::pwrite(descriptor, data, bytes, offset);
    if(written < 0 && errno == EINTR) { continue; }
    if(written <= 0) { throw joinError("write failed", path); }
    DiskIO::write_volume += static_cast<size_type>(written);
    data += written; bytes -= written; offset += written;
  }
}

void
preadAll(int descriptor, void* target, size_type bytes, off_t offset,
  const std::string& path)
{
  std::uint8_t* data = static_cast<std::uint8_t*>(target);
  while(bytes > 0)
  {
    ssize_t read_bytes = ::pread(descriptor, data, bytes, offset);
    if(read_bytes < 0 && errno == EINTR) { continue; }
    if(read_bytes == 0) { throw joinError("unexpected end of file", path); }
    if(read_bytes < 0) { throw joinError("read failed", path); }
    DiskIO::read_volume += static_cast<size_type>(read_bytes);
    data += read_bytes; bytes -= read_bytes; offset += read_bytes;
  }
}

template<class Value>
void
encodeLittle(std::uint8_t*& target, Value value)
{
  for(size_type i = 0; i < sizeof(Value); i++)
  {
    target[i] = static_cast<std::uint8_t>(value >> (8 * i));
  }
  target += sizeof(Value);
}

template<class Value>
Value
decodeLittle(const std::uint8_t*& source)
{
  Value result = 0;
  for(size_type i = 0; i < sizeof(Value); i++)
  {
    result |= static_cast<Value>(source[i]) << (8 * i);
  }
  source += sizeof(Value);
  return result;
}

std::uint64_t
joinChecksum(const void* source, size_type bytes,
  std::uint64_t result = 1469598103934665603ULL)
{
  const std::uint8_t* data = static_cast<const std::uint8_t*>(source);
  for(size_type i = 0; i < bytes; i++)
  {
    result ^= data[i]; result *= 1099511628211ULL;
  }
  return result;
}

void
encodeJoinRecord(const JoinRecord& source,
  std::array<std::uint8_t, JOIN_RECORD_BYTES>& buffer)
{
  std::uint8_t* out = buffer.data();
  encodeLittle<std::uint32_t>(out, source.logical.value);
  encodeLittle<node_type>(out, source.key);
  encodeLittle<node_type>(out, source.node.from);
  encodeLittle<node_type>(out, source.node.to);
  encodeLittle<std::uint8_t>(out, source.node.predecessors());
  encodeLittle<std::uint8_t>(out, static_cast<std::uint8_t>(source.node.order()));
  encodeLittle<std::uint8_t>(out, static_cast<std::uint8_t>(source.node.lcp()));
  encodeLittle<std::uint8_t>(out, 0); // Reserved for format-compatible flags.
  encodeLittle<std::uint64_t>(out, source.ordinal);
  for(size_type i = 0; i < JOIN_LABEL_COUNT; i++)
  {
    encodeLittle<PathNode::rank_type>(out, source.labels[i]);
  }
}

void
decodeJoinRecord(const std::array<std::uint8_t, JOIN_RECORD_BYTES>& buffer,
  JoinRecord& target)
{
  const std::uint8_t* in = buffer.data();
  target.logical = logical_file_id_t(decodeLittle<std::uint32_t>(in));
  target.key = decodeLittle<node_type>(in);
  target.node.from = decodeLittle<node_type>(in);
  target.node.to = decodeLittle<node_type>(in);
  std::uint8_t predecessors = decodeLittle<std::uint8_t>(in);
  std::uint8_t order = decodeLittle<std::uint8_t>(in);
  std::uint8_t lcp = decodeLittle<std::uint8_t>(in);
  static_cast<void>(decodeLittle<std::uint8_t>(in));
  target.ordinal = decodeLittle<std::uint64_t>(in);
  target.node.fields = 0;
  target.node.setPredecessors(predecessors);
  target.node.setOrder(order); target.node.setLCP(lcp); target.node.setPointer(0);
  // Order 0 is valid for a pruned sorted range and still carries one rank.
  if(order > PathLabel::LABEL_LENGTH || lcp > order)
  {
    throw joinError("invalid join record path metadata");
  }
  for(size_type i = 0; i < JOIN_LABEL_COUNT; i++)
  {
    target.labels[i] = decodeLittle<PathNode::rank_type>(in);
  }
}

void
encodeJoinGroupSummary(const JoinGroupSummary& source,
  std::array<std::uint8_t, JOIN_GROUP_SUMMARY_BYTES>& buffer)
{
  if(source.end < source.begin || source.count != source.end - source.begin)
  {
    throw joinError("invalid join-key group half-open range");
  }
  std::uint8_t* out = buffer.data();
  encodeLittle<node_type>(out, source.key);
  encodeLittle<std::uint64_t>(out, source.begin);
  encodeLittle<std::uint64_t>(out, source.end);
  encodeLittle<std::uint64_t>(out, source.order_sum);
  encodeLittle<std::uint64_t>(out, source.bypass_paths);
  encodeLittle<std::uint64_t>(out, source.bypass_ranks);
}

void
decodeJoinGroupSummary(const std::array<std::uint8_t, JOIN_GROUP_SUMMARY_BYTES>& buffer,
  JoinGroupSummary& target)
{
  const std::uint8_t* in = buffer.data();
  target.key = decodeLittle<node_type>(in);
  target.begin = decodeLittle<std::uint64_t>(in);
  target.end = decodeLittle<std::uint64_t>(in);
  target.order_sum = decodeLittle<std::uint64_t>(in);
  target.bypass_paths = decodeLittle<std::uint64_t>(in);
  target.bypass_ranks = decodeLittle<std::uint64_t>(in);
  if(target.end < target.begin)
  {
    throw joinError("invalid join-key group half-open range");
  }
  target.count = target.end - target.begin;
  if(target.bypass_paths > target.count || target.bypass_ranks < target.bypass_paths)
  {
    throw joinError("invalid join-key group aggregate");
  }
}

class JoinGroupSidecarWriter
{
public:
  JoinGroupSidecarWriter(logical_file_id_t logical, JoinKeyKind kind) :
    logical_id(logical), key_kind(kind), summary_name(TempFile::getName("gcsa_join_groups")),
    detail_name(TempFile::getName("gcsa_join_group_detail")),
    summary_descriptor(-1), detail_descriptor(-1), group_count(0), detail_records(0),
    summary_checksum(1469598103934665603ULL),
    detail_checksum(1469598103934665603ULL), summary_released(0), detail_released(0),
    active(), have_active(false), finished(false), summary_buffer(), detail_buffer()
  {
    this->summary_descriptor = ::open(this->summary_name.c_str(),
      O_CREAT | O_EXCL | O_WRONLY, 0644);
    if(this->summary_descriptor < 0)
    {
      throw joinError("cannot create join group sidecar", this->summary_name);
    }
    this->detail_descriptor = ::open(this->detail_name.c_str(),
      O_CREAT | O_EXCL | O_WRONLY, 0644);
    if(this->detail_descriptor < 0)
    {
      int saved_errno = errno; ::close(this->summary_descriptor);
      this->summary_descriptor = -1;
      errno = saved_errno;
      throw joinError("cannot create join group detail sidecar", this->detail_name);
    }
    adviseSequential(this->summary_descriptor); adviseSequential(this->detail_descriptor);
    this->writeHeaders();
    this->summary_buffer.reserve(std::max(static_cast<size_type>(1),
      JOIN_IO_BUFFER_BYTES / JOIN_GROUP_SUMMARY_BYTES));
    this->detail_buffer.reserve(JOIN_IO_BUFFER_BYTES);
  }

  ~JoinGroupSidecarWriter()
  {
    if(this->summary_descriptor >= 0) { ::close(this->summary_descriptor); }
    if(this->detail_descriptor >= 0) { ::close(this->detail_descriptor); }
    if(!this->finished)
    {
      TempFile::remove(this->summary_name); TempFile::remove(this->detail_name);
    }
  }

  void writeRecord(const JoinRecord& record, size_type offset)
  {
    if(offset != this->detail_records)
    {
      throw joinError("join group sidecar record offset mismatch", this->summary_name);
    }
    if(!this->have_active)
    {
      this->startGroup(record.key, offset);
    }
    else if(record.key != this->active.key)
    {
      if(record.key < this->active.key)
      {
        throw joinError("join run is not ordered by join key", this->summary_name);
      }
      this->finishGroup(offset); this->startGroup(record.key, offset);
    }

    size_type order = record.node.order();
    if(order > PathLabel::LABEL_LENGTH || order > JOIN_DETAIL_ORDER_MASK)
    {
      throw joinError("join group detail has invalid path order", this->detail_name);
    }
    std::uint8_t detail = static_cast<std::uint8_t>(order);
    if(record.node.sorted()) { detail |= JOIN_DETAIL_SORTED; }
    this->detail_buffer.push_back(detail);
    this->detail_checksum = joinChecksum(&detail, sizeof(detail), this->detail_checksum);
    this->detail_records++;
    this->active.count++;
    this->active.order_sum = checkedJoinAdd(this->active.order_sum, order,
      "group order sum");
    if(record.node.sorted())
    {
      this->active.bypass_paths++;
      this->active.bypass_ranks = checkedJoinAdd(this->active.bypass_ranks,
        order + 1, "group bypass ranks");
    }
    if(this->detail_buffer.size() >= JOIN_IO_BUFFER_BYTES) { this->flushDetails(); }
  }

  JoinGroupSidecar finish(size_type records, std::uint64_t run_checksum)
  {
    if(this->finished)
    {
      throw joinError("join group sidecar already finished", this->summary_name);
    }
    if(this->detail_records != records)
    {
      throw joinError("join group sidecar record count mismatch", this->summary_name);
    }
    if(this->have_active) { this->finishGroup(records); }
    this->flushSummaries(); this->flushDetails();
    this->writeFinalMetadata(records, run_checksum);
    this->writeFooters(records, run_checksum);
    trimWrittenCache(this->summary_descriptor, this->summary_released, true,
      this->summary_name);
    trimWrittenCache(this->detail_descriptor, this->detail_released, true,
      this->detail_name);
    if(::close(this->summary_descriptor) != 0)
    {
      throw joinError("cannot close join group sidecar", this->summary_name);
    }
    this->summary_descriptor = -1;
    if(::close(this->detail_descriptor) != 0)
    {
      throw joinError("cannot close join group detail sidecar", this->detail_name);
    }
    this->detail_descriptor = -1; this->finished = true;
    return JoinGroupSidecar(this->summary_name, this->detail_name, this->group_count,
      this->summary_checksum, this->detail_checksum);
  }

  void discard()
  {
    if(this->summary_descriptor >= 0)
    {
      ::close(this->summary_descriptor); this->summary_descriptor = -1;
    }
    if(this->detail_descriptor >= 0)
    {
      ::close(this->detail_descriptor); this->detail_descriptor = -1;
    }
    TempFile::remove(this->summary_name); TempFile::remove(this->detail_name);
    this->finished = true;
  }

private:
  void writeHeaders()
  {
    std::array<std::uint8_t, JOIN_GROUP_HEADER_BYTES> summary_header;
    std::uint8_t* summary = summary_header.data();
    encodeLittle<std::uint64_t>(summary, JOIN_GROUP_HEADER_MAGIC);
    encodeLittle<std::uint32_t>(summary, JOIN_GROUP_FORMAT_VERSION);
    encodeLittle<std::uint32_t>(summary, static_cast<std::uint32_t>(this->key_kind));
    encodeLittle<std::uint32_t>(summary, this->logical_id.value);
    encodeLittle<std::uint32_t>(summary, 0);
    encodeLittle<std::uint64_t>(summary, 0);
    encodeLittle<std::uint64_t>(summary, 0);
    encodeLittle<std::uint64_t>(summary, 0);
    writeAll(this->summary_descriptor, summary_header.data(), summary_header.size(),
      this->summary_name);

    std::array<std::uint8_t, JOIN_DETAIL_HEADER_BYTES> detail_header;
    std::uint8_t* details = detail_header.data();
    encodeLittle<std::uint64_t>(details, JOIN_DETAIL_HEADER_MAGIC);
    encodeLittle<std::uint32_t>(details, JOIN_GROUP_FORMAT_VERSION);
    encodeLittle<std::uint32_t>(details, static_cast<std::uint32_t>(this->key_kind));
    encodeLittle<std::uint32_t>(details, this->logical_id.value);
    encodeLittle<std::uint32_t>(details, 0);
    encodeLittle<std::uint64_t>(details, 0);
    encodeLittle<std::uint64_t>(details, 0);
    writeAll(this->detail_descriptor, detail_header.data(), detail_header.size(),
      this->detail_name);
  }

  void startGroup(node_type key, size_type begin)
  {
    this->active = JoinGroupSummary();
    this->active.key = key; this->active.begin = begin; this->active.end = begin;
    this->have_active = true;
  }

  void finishGroup(size_type end)
  {
    if(!this->have_active || end < this->active.begin ||
       this->active.count != end - this->active.begin)
    {
      throw joinError("invalid join group sidecar range", this->summary_name);
    }
    this->active.end = end;
    std::array<std::uint8_t, JOIN_GROUP_SUMMARY_BYTES> encoded;
    encodeJoinGroupSummary(this->active, encoded);
    this->summary_checksum = joinChecksum(encoded.data(), encoded.size(),
      this->summary_checksum);
    this->summary_buffer.push_back(encoded); this->group_count++;
    this->have_active = false;
    if(this->summary_buffer.size() * JOIN_GROUP_SUMMARY_BYTES >= JOIN_IO_BUFFER_BYTES)
    {
      this->flushSummaries();
    }
  }

  void flushSummaries()
  {
    if(this->summary_buffer.empty()) { return; }
    writeAll(this->summary_descriptor, this->summary_buffer.data(),
      this->summary_buffer.size() * JOIN_GROUP_SUMMARY_BYTES, this->summary_name);
    this->summary_buffer.clear();
    if(::lseek(this->summary_descriptor, 0, SEEK_CUR) - this->summary_released >=
       JOIN_CACHE_FLUSH_BYTES)
    {
      trimWrittenCache(this->summary_descriptor, this->summary_released, false,
        this->summary_name);
    }
  }

  void flushDetails()
  {
    if(this->detail_buffer.empty()) { return; }
    writeAll(this->detail_descriptor, this->detail_buffer.data(),
      this->detail_buffer.size(), this->detail_name);
    this->detail_buffer.clear();
    if(::lseek(this->detail_descriptor, 0, SEEK_CUR) - this->detail_released >=
       JOIN_CACHE_FLUSH_BYTES)
    {
      trimWrittenCache(this->detail_descriptor, this->detail_released, false,
        this->detail_name);
    }
  }

  void writeFinalMetadata(size_type records, std::uint64_t run_checksum)
  {
    std::array<std::uint8_t, 3 * sizeof(std::uint64_t)> summary_metadata;
    std::uint8_t* summary = summary_metadata.data();
    encodeLittle<std::uint64_t>(summary, records);
    encodeLittle<std::uint64_t>(summary, run_checksum);
    encodeLittle<std::uint64_t>(summary, this->group_count);
    pwriteAll(this->summary_descriptor, summary_metadata.data(), summary_metadata.size(),
      24, this->summary_name);

    std::array<std::uint8_t, 2 * sizeof(std::uint64_t)> detail_metadata;
    std::uint8_t* details = detail_metadata.data();
    encodeLittle<std::uint64_t>(details, records);
    encodeLittle<std::uint64_t>(details, run_checksum);
    pwriteAll(this->detail_descriptor, detail_metadata.data(), detail_metadata.size(),
      24, this->detail_name);
  }

  void writeFooters(size_type records, std::uint64_t run_checksum)
  {
    std::array<std::uint8_t, JOIN_GROUP_FOOTER_BYTES> summary_footer;
    std::uint8_t* summary = summary_footer.data();
    encodeLittle<std::uint64_t>(summary, JOIN_GROUP_FOOTER_MAGIC);
    encodeLittle<std::uint64_t>(summary, records);
    encodeLittle<std::uint64_t>(summary, run_checksum);
    encodeLittle<std::uint64_t>(summary, this->group_count);
    encodeLittle<std::uint64_t>(summary, this->summary_checksum);
    writeAll(this->summary_descriptor, summary_footer.data(), summary_footer.size(),
      this->summary_name);

    std::array<std::uint8_t, JOIN_DETAIL_FOOTER_BYTES> detail_footer;
    std::uint8_t* details = detail_footer.data();
    encodeLittle<std::uint64_t>(details, JOIN_DETAIL_FOOTER_MAGIC);
    encodeLittle<std::uint64_t>(details, records);
    encodeLittle<std::uint64_t>(details, run_checksum);
    encodeLittle<std::uint64_t>(details, this->detail_checksum);
    writeAll(this->detail_descriptor, detail_footer.data(), detail_footer.size(),
      this->detail_name);
  }

  JoinGroupSidecarWriter(const JoinGroupSidecarWriter&);
  JoinGroupSidecarWriter& operator=(const JoinGroupSidecarWriter&);

  logical_file_id_t logical_id;
  JoinKeyKind key_kind;
  std::string summary_name, detail_name;
  int summary_descriptor, detail_descriptor;
  size_type group_count, detail_records;
  std::uint64_t summary_checksum, detail_checksum;
  off_t summary_released, detail_released;
  JoinGroupSummary active;
  bool have_active, finished;
  std::vector<std::array<std::uint8_t, JOIN_GROUP_SUMMARY_BYTES>> summary_buffer;
  std::vector<std::uint8_t> detail_buffer;
};

bool
joinRecordLess(const JoinRecord& left, const JoinRecord& right)
{
  if(left.logical != right.logical) { return left.logical < right.logical; }
  if(left.key != right.key) { return left.key < right.key; }
  size_type order = std::min(left.node.order(), right.node.order());
  for(size_type i = 0; i < order; i++)
  {
    if(left.labels[i] != right.labels[i]) { return left.labels[i] < right.labels[i]; }
  }
  if(left.node.order() != right.node.order()) { return left.node.order() < right.node.order(); }
  if(left.node.from != right.node.from) { return left.node.from < right.node.from; }
  if(left.node.to != right.node.to) { return left.node.to < right.node.to; }
  if(left.node.predecessors() != right.node.predecessors())
  {
    return left.node.predecessors() < right.node.predecessors();
  }
  if(left.node.lcp() != right.node.lcp()) { return left.node.lcp() < right.node.lcp(); }
  for(size_type i = order; i < left.node.ranks(); i++)
  {
    if(left.labels[i] != right.labels[i]) { return left.labels[i] < right.labels[i]; }
  }
  return left.ordinal < right.ordinal;
}

class JoinFileWriter
{
public:
  JoinFileWriter(const std::string& path, logical_file_id_t logical, JoinKeyKind kind,
    size_type expected_records, const TempFileCodecParameters& codec,
    ExternalPathJoinStats* stats) :
    name(path), logical_id(logical), key_kind(kind), descriptor(-1),
    expected_count(expected_records), count(0),
    checksum(1469598103934665603ULL), cache_released(0), finished(false),
    statistics(stats), sidecar(logical, kind), compressed(),
    batch_records(JOIN_IO_BUFFER_RECORDS), buffer()
  {
    try
    {
      if(codec.enabled())
      {
        const size_type minimum_block = std::max(JOIN_RECORD_BYTES,
          std::max(JOIN_HEADER_BYTES, JOIN_FOOTER_BYTES));
        if(codec.block_size < minimum_block)
        {
          throw joinError("compressed join-run block is smaller than one logical record",
            this->name);
        }
        this->compressed.reset(new CompressedBlockWriter(this->name,
          codec.block_size, CompressedBlockWriter::ZSTD, codec.level,
          codec.workers));
        this->batch_records = std::min(JOIN_IO_BUFFER_RECORDS,
          codec.block_size / JOIN_RECORD_BYTES);
      }
      else
      {
        this->descriptor = ::open(this->name.c_str(),
          O_CREAT | O_EXCL | O_WRONLY, 0644);
        if(this->descriptor < 0)
        {
          throw joinError("cannot create join run", this->name);
        }
        adviseSequential(this->descriptor);
      }
      this->buffer.reserve(this->batch_records);
      std::array<std::uint8_t, JOIN_HEADER_BYTES> header;
      std::uint8_t* out = header.data();
      encodeLittle<std::uint64_t>(out, JOIN_HEADER_MAGIC);
      encodeLittle<std::uint32_t>(out, JOIN_FORMAT_VERSION);
      encodeLittle<std::uint32_t>(out, static_cast<std::uint32_t>(this->key_kind));
      encodeLittle<std::uint32_t>(out, this->logical_id.value);
      encodeLittle<std::uint64_t>(out, this->expected_count);
      this->writeLogicalRecord(header.data(), header.size());
    }
    catch(...)
    {
      if(this->descriptor >= 0) { ::close(this->descriptor); this->descriptor = -1; }
      this->compressed.reset();
      ::unlink(this->name.c_str());
      throw;
    }
  }

  ~JoinFileWriter()
  {
    if(this->descriptor >= 0) { ::close(this->descriptor); }
    if(!this->finished) { ::unlink(this->name.c_str()); }
  }

  void writeRecord(const JoinRecord& record)
  {
    this->sidecar.writeRecord(record, this->count);
    std::array<std::uint8_t, JOIN_RECORD_BYTES> encoded;
    encodeJoinRecord(record, encoded);
    this->checksum = joinChecksum(encoded.data(), encoded.size(), this->checksum);
    this->buffer.push_back(encoded);
    this->count++;
    if(this->buffer.size() >= this->batch_records)
    {
      this->flushRecords();
    }
  }

  JoinRun finish()
  {
    if(this->finished) { throw joinError("join run already finished", this->name); }
    try
    {
      if(this->count != this->expected_count)
      {
        throw joinError("join run record count differs from its declared count",
          this->name);
      }
      this->flushRecords();
      std::array<std::uint8_t, JOIN_FOOTER_BYTES> footer;
      std::uint8_t* out = footer.data();
      encodeLittle<std::uint64_t>(out, JOIN_FOOTER_MAGIC);
      encodeLittle<std::uint64_t>(out, this->count);
      encodeLittle<std::uint64_t>(out, this->checksum);
      this->writeLogicalRecord(footer.data(), footer.size());
      if(this->compressed)
      {
        this->compressed->finish();
      }
      else
      {
        trimWrittenCache(this->descriptor, this->cache_released, true, this->name);
        if(::close(this->descriptor) != 0)
        {
          this->descriptor = -1;
          throw joinError("cannot close join run", this->name);
        }
        this->descriptor = -1;
      }
      JoinGroupSidecar group_sidecar = this->sidecar.finish(this->count,
        this->checksum);
      if(this->statistics != nullptr)
      {
        size_type logical_bytes = checkedJoinAdd(JOIN_HEADER_BYTES,
          checkedJoinAdd(checkedJoinMultiply(this->count, JOIN_RECORD_BYTES,
            "join-run logical payload bytes"), JOIN_FOOTER_BYTES,
            "join-run logical footer bytes"), "join-run logical bytes");
        struct stat info;
        if(::stat(this->name.c_str(), &info) != 0 || info.st_size < 0 ||
           static_cast<std::uintmax_t>(info.st_size) >
             std::numeric_limits<size_type>::max())
        {
          throw joinError("cannot determine stored join-run bytes", this->name);
        }
        this->statistics->join_run_logical_bytes = checkedJoinAdd(
          this->statistics->join_run_logical_bytes, logical_bytes,
          "aggregate join-run logical bytes");
        this->statistics->join_run_stored_bytes = checkedJoinAdd(
          this->statistics->join_run_stored_bytes,
          static_cast<size_type>(info.st_size), "aggregate stored join-run bytes");
        if(this->compressed)
        {
          this->statistics->compressed_join_runs = checkedJoinAdd(
            this->statistics->compressed_join_runs, 1,
            "compressed join-run count");
        }
      }
      this->finished = true;
      return JoinRun(this->name, this->count, this->checksum, group_sidecar);
    }
    catch(...)
    {
      this->sidecar.discard(); throw;
    }
  }

private:
  void writeLogicalRecord(const void* data, size_type bytes)
  {
    if(this->compressed) { this->compressed->writeRecord(data, bytes); }
    else { writeAll(this->descriptor, data, bytes, this->name); }
  }

  void flushRecords()
  {
    if(this->buffer.empty()) { return; }
    size_type bytes = this->buffer.size() * JOIN_RECORD_BYTES;
    if(this->compressed)
    {
      // Preserve the raw writer's aligned byte staging so large runs do not
      // pay one vector insertion and one API call per 108-byte join record.
      this->compressed->writeRecord(this->buffer.data(), bytes);
    }
    else { writeAll(this->descriptor, this->buffer.data(), bytes, this->name); }
    this->buffer.clear();
    if(!this->compressed)
    {
      off_t written = JOIN_HEADER_BYTES + this->count * JOIN_RECORD_BYTES;
      if(written - this->cache_released >= JOIN_CACHE_FLUSH_BYTES)
      {
        trimWrittenCache(this->descriptor, this->cache_released, false, this->name);
      }
    }
  }

  JoinFileWriter(const JoinFileWriter&);
  JoinFileWriter& operator=(const JoinFileWriter&);

  std::string name;
  logical_file_id_t logical_id;
  JoinKeyKind key_kind;
  int descriptor;
  size_type expected_count, count;
  std::uint64_t checksum;
  off_t cache_released;
  bool finished;
  ExternalPathJoinStats* statistics;
  JoinGroupSidecarWriter sidecar;
  std::unique_ptr<CompressedBlockWriter> compressed;
  size_type batch_records;
  std::vector<std::array<std::uint8_t, JOIN_RECORD_BYTES>> buffer;
};

class JoinFileReader
{
public:
  JoinFileReader(const JoinRun& run, logical_file_id_t logical, JoinKeyKind kind) :
    name(run.name), descriptor(-1), record_count(0), expected_checksum(0),
    cache_released(0), buffer(JOIN_IO_BUFFER_RECORDS), buffer_first(0),
    buffer_records(0), compressed()
  {
    try
    {
      if(CompressedBlockReader::isFramed(this->name))
      {
        this->compressed.reset(new CompressedBlockReader(this->name));
      }
      else
      {
        this->descriptor = ::open(this->name.c_str(), O_RDONLY);
        if(this->descriptor < 0)
        {
          throw joinError("cannot open join run", this->name);
        }
        adviseSequential(this->descriptor);
      }
      std::array<std::uint8_t, JOIN_HEADER_BYTES> header;
      this->readLogical(header.data(), header.size(), 0);
      const std::uint8_t* in = header.data();
      std::uint64_t magic = decodeLittle<std::uint64_t>(in);
      std::uint32_t version = decodeLittle<std::uint32_t>(in);
      std::uint32_t actual_kind = decodeLittle<std::uint32_t>(in);
      std::uint32_t actual_logical = decodeLittle<std::uint32_t>(in);
      this->record_count = decodeLittle<std::uint64_t>(in);
      if(magic != JOIN_HEADER_MAGIC || version != JOIN_FORMAT_VERSION ||
         actual_kind != static_cast<std::uint32_t>(kind) ||
         actual_logical != logical.value || this->record_count != run.records)
      {
        throw joinError("join run header mismatch", this->name);
      }
      if(this->record_count > (std::numeric_limits<size_type>::max() -
         JOIN_HEADER_BYTES - JOIN_FOOTER_BYTES) / JOIN_RECORD_BYTES)
      {
        throw joinError("join run is too large for this build", this->name);
      }
      size_type logical_bytes = JOIN_HEADER_BYTES +
        this->record_count * JOIN_RECORD_BYTES + JOIN_FOOTER_BYTES;
      if(this->compressed)
      {
        // Outer records are operational batches and may vary between writers;
        // the byte-identical inner JOIN v1 stream is the compatibility contract.
        if(this->compressed->logicalSize() != logical_bytes)
        {
          throw joinError("framed join run length mismatch", this->name);
        }
      }
      else
      {
        struct stat info;
        if(::fstat(this->descriptor, &info) != 0 || info.st_size < 0 ||
           static_cast<size_type>(info.st_size) != logical_bytes)
        {
          throw joinError("join run length mismatch", this->name);
        }
      }
      std::array<std::uint8_t, JOIN_FOOTER_BYTES> footer;
      this->readLogical(footer.data(), footer.size(),
        JOIN_HEADER_BYTES + this->record_count * JOIN_RECORD_BYTES);
      in = footer.data();
      if(decodeLittle<std::uint64_t>(in) != JOIN_FOOTER_MAGIC ||
         decodeLittle<std::uint64_t>(in) != this->record_count)
      {
        throw joinError("join run footer mismatch", this->name);
      }
      this->expected_checksum = decodeLittle<std::uint64_t>(in);
      if(run.checksum != 0 && this->expected_checksum != run.checksum)
      {
        throw joinError("join run checksum identity mismatch", this->name);
      }
    }
    catch(...)
    {
      if(this->descriptor >= 0) { ::close(this->descriptor); this->descriptor = -1; }
      this->compressed.reset();
      throw;
    }
  }

  ~JoinFileReader()
  {
    if(this->descriptor >= 0)
    {
      discardCachedRange(this->descriptor, 0, 0); ::close(this->descriptor);
    }
  }

  size_type size() const { return this->record_count; }
  // The checksum is part of the immutable run identity. Partition plans bind
  // to it so a same-length replacement run can never inherit stale ranges.
  std::uint64_t checksum() const { return this->expected_checksum; }

  void read(size_type index, JoinRecord& record) const
  {
    if(index >= this->record_count) { throw joinError("join record index out of range", this->name); }
    if(this->buffer_records == 0 || index < this->buffer_first ||
       index >= this->buffer_first + this->buffer_records)
    {
      this->refill(index);
    }
    decodeJoinRecord(this->buffer[index - this->buffer_first], record);
  }

  void validate() const
  {
    std::uint64_t checksum = 1469598103934665603ULL;
    if(this->compressed)
    {
      // Rewind and consume the complete logical stream in order. This checks
      // the framed whole-stream checksum in addition to the join payload
      // checksum used as the immutable run identity.
      this->compressed->seekBlock(0);
      std::array<std::uint8_t, JOIN_HEADER_BYTES> header;
      if(this->compressed->read(header.data(), header.size()) != header.size())
      {
        throw joinError("truncated framed join run header", this->name);
      }
      for(size_type first = 0; first < this->record_count;
          first += this->buffer.size())
      {
        size_type records = std::min(this->buffer.size(),
          this->record_count - first);
        size_type bytes = records * JOIN_RECORD_BYTES;
        if(this->compressed->read(this->buffer.data(), bytes) != bytes)
        {
          throw joinError("truncated framed join run payload", this->name);
        }
        checksum = joinChecksum(this->buffer.data(), bytes, checksum);
      }
      std::array<std::uint8_t, JOIN_FOOTER_BYTES> footer;
      if(this->compressed->read(footer.data(), footer.size()) != footer.size())
      {
        throw joinError("truncated framed join run footer", this->name);
      }
    }
    else
    {
      for(size_type first = 0; first < this->record_count;
          first += this->buffer.size())
      {
        size_type records = std::min(this->buffer.size(),
          this->record_count - first);
        preadAll(this->descriptor, this->buffer.data(),
          records * JOIN_RECORD_BYTES,
          JOIN_HEADER_BYTES + first * JOIN_RECORD_BYTES, this->name);
        checksum = joinChecksum(this->buffer.data(),
          records * JOIN_RECORD_BYTES, checksum);
        trimReadCache(this->descriptor,
          JOIN_HEADER_BYTES + (first + records) * JOIN_RECORD_BYTES,
          this->cache_released);
      }
    }
    if(checksum != this->expected_checksum) { throw joinError("join run checksum mismatch", this->name); }
    if(this->descriptor >= 0) { discardCachedRange(this->descriptor, 0, 0); }
    this->cache_released = 0;
    this->buffer_records = 0;
  }

private:
  void readLogical(void* target, size_type bytes, size_type offset) const
  {
    if(this->compressed)
    {
      if(this->compressed->readAt(offset, target, bytes) != bytes)
      {
        throw joinError("unexpected end of framed join run", this->name);
      }
    }
    else
    {
      preadAll(this->descriptor, target, bytes, static_cast<off_t>(offset),
        this->name);
    }
  }

  void refill(size_type first) const
  {
    this->buffer_first = first;
    this->buffer_records = std::min(this->buffer.size(), this->record_count - first);
    this->readLogical(this->buffer.data(),
      this->buffer_records * JOIN_RECORD_BYTES,
      JOIN_HEADER_BYTES + first * JOIN_RECORD_BYTES);
    if(this->descriptor >= 0)
    {
      trimReadCache(this->descriptor,
        JOIN_HEADER_BYTES + (first + this->buffer_records) * JOIN_RECORD_BYTES,
        this->cache_released);
    }
  }

  JoinFileReader(const JoinFileReader&);
  JoinFileReader& operator=(const JoinFileReader&);

  std::string name;
  int descriptor;
  size_type record_count;
  std::uint64_t expected_checksum;
  mutable off_t cache_released;
  mutable std::vector<std::array<std::uint8_t, JOIN_RECORD_BYTES>> buffer;
  mutable size_type buffer_first, buffer_records;
  mutable std::unique_ptr<CompressedBlockReader> compressed;
};

size_type
joinRawReaderReservation()
{
  return sizeof(JoinRecord) + sizeof(size_type) +
    sizeof(std::unique_ptr<JoinFileReader>) + JOIN_IO_BUFFER_BYTES + 128;
}

size_type
joinRawScanBaseBudget()
{
  return checkedJoinAdd(JOIN_FIXED_BYTES,
    checkedJoinMultiply(2, joinRawReaderReservation(),
      "minimum join readers"), "raw join scan base");
}

size_type
joinRawMinimumBudget()
{
  // Even the exact-minimum join must have one materialization record beyond
  // its two readers and fixed stream state.
  return checkedJoinAdd(joinRawScanBaseBudget(), sizeof(JoinRecord),
    "minimum raw join budget");
}

size_type
joinRunWriterMemory(const TempFileCodecParameters& codec)
{
  if(!codec.enabled()) { return 0; }
  return CompressedBlockWriter::workingMemoryEstimate(codec.block_size,
    CompressedBlockWriter::ZSTD, codec.level, codec.workers);
}

size_type
joinRunWriterExtra(const TempFileCodecParameters& codec)
{
  // JOIN_FIXED_BYTES owns the aligned JOIN staging buffer in both formats.
  // Framing adds its block, compressed output, and codec context on top.
  return joinRunWriterMemory(codec);
}

size_type
joinRunReaderMemory(const TempFileCodecParameters& codec)
{
  if(!codec.enabled()) { return 0; }
  return CompressedBlockReader::workingMemoryEstimate(codec.block_size);
}

size_type
joinRunReaderMemory(const JoinRun& run)
{
  if(!CompressedBlockReader::isFramed(run.name)) { return 0; }
  std::uint64_t block = CompressedBlockReader::declaredBlockSize(run.name);
  if(block > std::numeric_limits<size_type>::max())
  {
    throw joinError("compressed join-run block exceeds this build", run.name);
  }
  return CompressedBlockReader::workingMemoryEstimate(
    static_cast<size_type>(block));
}

size_type
joinRunReaderPairMemory(const JoinRun& left, const JoinRun& right)
{
  return checkedJoinAdd(joinRunReaderMemory(left), joinRunReaderMemory(right),
    "join-run reader codec workspaces");
}

size_type
joinRunReadMinimumBudget(const TempFileCodecParameters& codec)
{
  return checkedJoinAdd(joinRawMinimumBudget(),
    checkedJoinMultiply(2, joinRunReaderMemory(codec),
      "minimum compressed join readers"),
    "minimum compressed join scan budget");
}

size_type
joinRunSortMinimumBudget(const TempFileCodecParameters& codec)
{
  return checkedJoinAdd(joinRunReadMinimumBudget(codec),
    joinRunWriterExtra(codec), "minimum compressed join sorter budget");
}

TempFileCodecParameters
boundedJoinRunCodec(const TempFileCodecParameters& requested,
  size_type join_sort_budget)
{
  if(!requested.enabled()) { return requested; }
  TempFileCodecParameters result = requested;
  const size_type minimum_block = std::max(JOIN_RECORD_BYTES,
    std::max(JOIN_HEADER_BYTES, JOIN_FOOTER_BYTES));
  result.block_size = std::max(minimum_block,
    std::min(result.block_size, join_sort_budget));
  auto fits = [&](size_type block_bytes) -> bool
  {
    try
    {
      result.block_size = block_bytes;
      // This is the join sorter's own phase share. The concurrently active
      // label sorter has an independent reservation outside this ceiling.
      return joinRunSortMinimumBudget(result) <= join_sort_budget;
    }
    catch(const std::exception&)
    {
      return false;
    }
  };
  while(result.block_size > minimum_block && !fits(result.block_size))
  {
    result.block_size = std::max(minimum_block, result.block_size / 2);
  }
  if(!fits(result.block_size))
  {
    if(requested.compression == TempCompression::AUTO)
    {
      result.compression = TempCompression::NONE;
      return result;
    }
    throw joinError(
      "join-sort budget cannot admit compressed join-run writer and readers");
  }
  return result;
}

class JoinGroupSidecarReader
{
public:
  JoinGroupSidecarReader(const JoinRun& run, logical_file_id_t logical,
    JoinKeyKind kind) :
    summary_name(run.sidecar.summary_name), detail_name(run.sidecar.detail_name),
    summary_descriptor(-1), detail_descriptor(-1), record_count(run.records),
    group_count(0), expected_summary_checksum(0), expected_detail_checksum(0),
    summary_cache_released(0), detail_cache_released(0),
    summary_buffer(std::max(static_cast<size_type>(1),
      JOIN_IO_BUFFER_BYTES / JOIN_GROUP_SUMMARY_BYTES)),
    summary_buffer_first(0), summary_buffer_records(0),
    detail_buffer(JOIN_IO_BUFFER_BYTES)
  {
    if(!run.sidecar.complete())
    {
      throw joinError("join run has no group sidecar", run.name);
    }
    this->summary_descriptor = ::open(this->summary_name.c_str(), O_RDONLY);
    this->detail_descriptor = ::open(this->detail_name.c_str(), O_RDONLY);
    if(this->summary_descriptor < 0 || this->detail_descriptor < 0)
    {
      throw joinError("cannot open join group sidecar", this->summary_name);
    }
    adviseSequential(this->summary_descriptor); adviseSequential(this->detail_descriptor);

    std::array<std::uint8_t, JOIN_GROUP_HEADER_BYTES> summary_header;
    preadAll(this->summary_descriptor, summary_header.data(), summary_header.size(), 0,
      this->summary_name);
    const std::uint8_t* input = summary_header.data();
    std::uint64_t summary_magic = decodeLittle<std::uint64_t>(input);
    std::uint32_t summary_version = decodeLittle<std::uint32_t>(input);
    std::uint32_t summary_kind = decodeLittle<std::uint32_t>(input);
    std::uint32_t summary_logical = decodeLittle<std::uint32_t>(input);
    std::uint32_t summary_reserved = decodeLittle<std::uint32_t>(input);
    size_type summary_records = decodeLittle<std::uint64_t>(input);
    std::uint64_t summary_run_checksum = decodeLittle<std::uint64_t>(input);
    this->group_count = decodeLittle<std::uint64_t>(input);
    if(summary_magic != JOIN_GROUP_HEADER_MAGIC ||
       summary_version != JOIN_GROUP_FORMAT_VERSION ||
       summary_kind != static_cast<std::uint32_t>(kind) ||
       summary_logical != logical.value || summary_reserved != 0 ||
       summary_records != run.records || summary_run_checksum != run.checksum ||
       this->group_count != run.sidecar.groups ||
       this->group_count > (std::numeric_limits<size_type>::max() -
         JOIN_GROUP_HEADER_BYTES - JOIN_GROUP_FOOTER_BYTES) / JOIN_GROUP_SUMMARY_BYTES)
    {
      throw joinError("join group sidecar header mismatch", this->summary_name);
    }
    struct stat summary_info;
    size_type summary_bytes = JOIN_GROUP_HEADER_BYTES +
      this->group_count * JOIN_GROUP_SUMMARY_BYTES + JOIN_GROUP_FOOTER_BYTES;
    if(::fstat(this->summary_descriptor, &summary_info) != 0 ||
       summary_info.st_size < 0 || static_cast<size_type>(summary_info.st_size) != summary_bytes)
    {
      throw joinError("join group sidecar length mismatch", this->summary_name);
    }
    std::array<std::uint8_t, JOIN_GROUP_FOOTER_BYTES> summary_footer;
    preadAll(this->summary_descriptor, summary_footer.data(), summary_footer.size(),
      JOIN_GROUP_HEADER_BYTES + this->group_count * JOIN_GROUP_SUMMARY_BYTES,
      this->summary_name);
    input = summary_footer.data();
    if(decodeLittle<std::uint64_t>(input) != JOIN_GROUP_FOOTER_MAGIC ||
       decodeLittle<std::uint64_t>(input) != run.records ||
       decodeLittle<std::uint64_t>(input) != run.checksum ||
       decodeLittle<std::uint64_t>(input) != this->group_count)
    {
      throw joinError("join group sidecar footer mismatch", this->summary_name);
    }
    this->expected_summary_checksum = decodeLittle<std::uint64_t>(input);
    if(this->expected_summary_checksum != run.sidecar.summary_checksum)
    {
      throw joinError("join group sidecar checksum identity mismatch", this->summary_name);
    }

    std::array<std::uint8_t, JOIN_DETAIL_HEADER_BYTES> detail_header;
    preadAll(this->detail_descriptor, detail_header.data(), detail_header.size(), 0,
      this->detail_name);
    input = detail_header.data();
    std::uint64_t detail_magic = decodeLittle<std::uint64_t>(input);
    std::uint32_t detail_version = decodeLittle<std::uint32_t>(input);
    std::uint32_t detail_kind = decodeLittle<std::uint32_t>(input);
    std::uint32_t detail_logical = decodeLittle<std::uint32_t>(input);
    std::uint32_t detail_reserved = decodeLittle<std::uint32_t>(input);
    size_type detail_records = decodeLittle<std::uint64_t>(input);
    std::uint64_t detail_run_checksum = decodeLittle<std::uint64_t>(input);
    if(detail_magic != JOIN_DETAIL_HEADER_MAGIC ||
       detail_version != JOIN_GROUP_FORMAT_VERSION ||
       detail_kind != static_cast<std::uint32_t>(kind) ||
       detail_logical != logical.value || detail_reserved != 0 ||
       detail_records != run.records || detail_run_checksum != run.checksum ||
       detail_records > std::numeric_limits<size_type>::max() -
         JOIN_DETAIL_HEADER_BYTES - JOIN_DETAIL_FOOTER_BYTES)
    {
      throw joinError("join group detail header mismatch", this->detail_name);
    }
    struct stat detail_info;
    size_type detail_bytes = JOIN_DETAIL_HEADER_BYTES + detail_records +
      JOIN_DETAIL_FOOTER_BYTES;
    if(::fstat(this->detail_descriptor, &detail_info) != 0 ||
       detail_info.st_size < 0 || static_cast<size_type>(detail_info.st_size) != detail_bytes)
    {
      throw joinError("join group detail length mismatch", this->detail_name);
    }
    std::array<std::uint8_t, JOIN_DETAIL_FOOTER_BYTES> detail_footer;
    preadAll(this->detail_descriptor, detail_footer.data(), detail_footer.size(),
      JOIN_DETAIL_HEADER_BYTES + detail_records, this->detail_name);
    input = detail_footer.data();
    if(decodeLittle<std::uint64_t>(input) != JOIN_DETAIL_FOOTER_MAGIC ||
       decodeLittle<std::uint64_t>(input) != run.records ||
       decodeLittle<std::uint64_t>(input) != run.checksum)
    {
      throw joinError("join group detail footer mismatch", this->detail_name);
    }
    this->expected_detail_checksum = decodeLittle<std::uint64_t>(input);
    if(this->expected_detail_checksum != run.sidecar.detail_checksum)
    {
      throw joinError("join group detail checksum identity mismatch", this->detail_name);
    }
  }

  ~JoinGroupSidecarReader()
  {
    if(this->summary_descriptor >= 0)
    {
      discardCachedRange(this->summary_descriptor, 0, 0);
      ::close(this->summary_descriptor);
    }
    if(this->detail_descriptor >= 0)
    {
      discardCachedRange(this->detail_descriptor, 0, 0);
      ::close(this->detail_descriptor);
    }
  }

  size_type size() const { return this->group_count; }
  size_type records() const { return this->record_count; }
  std::uint64_t summaryChecksum() const { return this->expected_summary_checksum; }
  std::uint64_t detailChecksum() const { return this->expected_detail_checksum; }

  void read(size_type index, JoinGroupSummary& summary) const
  {
    if(index >= this->group_count)
    {
      throw joinError("join group sidecar index out of range", this->summary_name);
    }
    if(this->summary_buffer_records == 0 || index < this->summary_buffer_first ||
       index >= this->summary_buffer_first + this->summary_buffer_records)
    {
      this->refillSummaries(index);
    }
    decodeJoinGroupSummary(this->summary_buffer[index - this->summary_buffer_first],
      summary);
    if(summary.end > this->record_count)
    {
      throw joinError("join group sidecar range exceeds join run", this->summary_name);
    }
  }

  JoinGroupSummary summarize(const JoinGroupSummary& parent, size_type begin,
    size_type end) const
  {
    if(begin >= end || begin < parent.begin || end > parent.end)
    {
      throw joinError("invalid join-key subgroup sidecar range", this->detail_name);
    }
    if(begin == parent.begin && end == parent.end) { return parent; }
    JoinGroupSummary result;
    result.key = parent.key; result.begin = begin; result.end = end;
    result.count = end - begin;
    for(size_type offset = begin; offset < end; )
    {
      size_type bytes = std::min(this->detail_buffer.size(), end - offset);
      preadAll(this->detail_descriptor, this->detail_buffer.data(), bytes,
        JOIN_DETAIL_HEADER_BYTES + offset, this->detail_name);
      trimReadCache(this->detail_descriptor, JOIN_DETAIL_HEADER_BYTES + offset + bytes,
        this->detail_cache_released);
      for(size_type i = 0; i < bytes; i++)
      {
        std::uint8_t detail = this->detail_buffer[i];
        size_type order = detail & JOIN_DETAIL_ORDER_MASK;
        if((detail & JOIN_DETAIL_RESERVED) != 0 || order > PathLabel::LABEL_LENGTH)
        {
          throw joinError("invalid join group detail record", this->detail_name);
        }
        result.order_sum = checkedJoinAdd(result.order_sum, order,
          "subgroup order sum");
        if((detail & JOIN_DETAIL_SORTED) != 0)
        {
          result.bypass_paths++;
          result.bypass_ranks = checkedJoinAdd(result.bypass_ranks, order + 1,
            "subgroup bypass ranks");
        }
      }
      offset += bytes;
    }
    return result;
  }

  void validate() const
  {
    std::uint64_t summary_checksum = 1469598103934665603ULL;
    size_type cursor = 0;
    node_type previous_key = 0;
    bool have_previous = false;
    std::array<std::uint8_t, JOIN_GROUP_SUMMARY_BYTES> encoded;
    for(size_type group = 0; group < this->group_count; group++)
    {
      preadAll(this->summary_descriptor, encoded.data(), encoded.size(),
        JOIN_GROUP_HEADER_BYTES + group * JOIN_GROUP_SUMMARY_BYTES,
        this->summary_name);
      summary_checksum = joinChecksum(encoded.data(), encoded.size(), summary_checksum);
      trimReadCache(this->summary_descriptor, JOIN_GROUP_HEADER_BYTES +
        (group + 1) * JOIN_GROUP_SUMMARY_BYTES, this->summary_cache_released);
      JoinGroupSummary summary;
      decodeJoinGroupSummary(encoded, summary);
      if(summary.begin != cursor || summary.end > this->record_count ||
         (have_previous && summary.key <= previous_key))
      {
        throw joinError("join group sidecar is not canonical", this->summary_name);
      }
      cursor = summary.end; previous_key = summary.key; have_previous = true;
    }
    if(cursor != this->record_count || summary_checksum != this->expected_summary_checksum)
    {
      throw joinError("join group sidecar checksum mismatch", this->summary_name);
    }
    std::uint64_t detail_checksum = 1469598103934665603ULL;
    for(size_type offset = 0; offset < this->record_count; )
    {
      size_type bytes = std::min(this->detail_buffer.size(), this->record_count - offset);
      preadAll(this->detail_descriptor, this->detail_buffer.data(), bytes,
        JOIN_DETAIL_HEADER_BYTES + offset, this->detail_name);
      detail_checksum = joinChecksum(this->detail_buffer.data(), bytes, detail_checksum);
      trimReadCache(this->detail_descriptor, JOIN_DETAIL_HEADER_BYTES + offset + bytes,
        this->detail_cache_released);
      for(size_type i = 0; i < bytes; i++)
      {
        size_type order = this->detail_buffer[i] & JOIN_DETAIL_ORDER_MASK;
        if((this->detail_buffer[i] & JOIN_DETAIL_RESERVED) != 0 ||
           order > PathLabel::LABEL_LENGTH)
        {
          throw joinError("invalid join group detail record", this->detail_name);
        }
      }
      offset += bytes;
    }
    if(detail_checksum != this->expected_detail_checksum)
    {
      throw joinError("join group detail checksum mismatch", this->detail_name);
    }
    discardCachedRange(this->summary_descriptor, 0, 0);
    discardCachedRange(this->detail_descriptor, 0, 0);
    this->summary_cache_released = 0; this->detail_cache_released = 0;
    this->summary_buffer_records = 0;
  }

private:
  void refillSummaries(size_type first) const
  {
    this->summary_buffer_first = first;
    this->summary_buffer_records = std::min(this->summary_buffer.size(),
      this->group_count - first);
    preadAll(this->summary_descriptor, this->summary_buffer.data(),
      this->summary_buffer_records * JOIN_GROUP_SUMMARY_BYTES,
      JOIN_GROUP_HEADER_BYTES + first * JOIN_GROUP_SUMMARY_BYTES,
      this->summary_name);
    trimReadCache(this->summary_descriptor, JOIN_GROUP_HEADER_BYTES +
      (first + this->summary_buffer_records) * JOIN_GROUP_SUMMARY_BYTES,
      this->summary_cache_released);
  }

  JoinGroupSidecarReader(const JoinGroupSidecarReader&);
  JoinGroupSidecarReader& operator=(const JoinGroupSidecarReader&);

  std::string summary_name, detail_name;
  int summary_descriptor, detail_descriptor;
  size_type record_count, group_count;
  std::uint64_t expected_summary_checksum, expected_detail_checksum;
  mutable off_t summary_cache_released, detail_cache_released;
  mutable std::vector<std::array<std::uint8_t, JOIN_GROUP_SUMMARY_BYTES>> summary_buffer;
  mutable size_type summary_buffer_first, summary_buffer_records;
  mutable std::vector<std::uint8_t> detail_buffer;
};

struct JoinHeapComparator
{
  const std::vector<JoinRecord>* records;

  explicit JoinHeapComparator(const std::vector<JoinRecord>* source = nullptr) : records(source) { }

  bool operator()(size_type left, size_type right) const
  {
    if(joinRecordLess((*this->records)[left], (*this->records)[right])) { return false; }
    if(joinRecordLess((*this->records)[right], (*this->records)[left])) { return true; }
    return left > right;
  }
};

JoinRun
mergeJoinRuns(const std::vector<JoinRun>& inputs, logical_file_id_t logical,
  JoinKeyKind kind, size_type byte_budget,
  const TempFileCodecParameters& codec,
  bool verify_payloads, ExternalPathJoinStats* stats)
{
  if(inputs.empty()) { throw joinError("cannot merge an empty run set"); }
  size_type total_records = 0, reader_codec_bytes = 0;
  for(size_type i = 0; i < inputs.size(); i++)
  {
    total_records = checkedJoinAdd(total_records, inputs[i].records,
      "merged join-run records");
    reader_codec_bytes = checkedJoinAdd(reader_codec_bytes,
      joinRunReaderMemory(inputs[i]), "merged join-run reader workspaces");
  }
  size_type merge_reservation = checkedJoinAdd(
    checkedJoinAdd(JOIN_FIXED_BYTES, joinRunWriterExtra(codec),
      "join merge writer workspace"),
    checkedJoinAdd(reader_codec_bytes,
      checkedJoinMultiply(inputs.size(), joinRawReaderReservation(),
        "join merge raw readers"),
      "join merge readers"), "join merge reservation");
  if(merge_reservation > byte_budget)
  {
    throw joinError("join merge exceeds its byte budget");
  }
  std::vector<std::unique_ptr<JoinFileReader>> readers;
  std::vector<JoinRecord> current(inputs.size());
  std::vector<size_type> offsets(inputs.size(), 0);
  readers.reserve(inputs.size());
  for(size_type i = 0; i < inputs.size(); i++)
  {
    readers.emplace_back(new JoinFileReader(inputs[i], logical, kind));
    if(verify_payloads) { readers.back()->validate(); }
    if(readers.back()->size() > 0) { readers.back()->read(0, current[i]); }
  }
  std::priority_queue<size_type, std::vector<size_type>, JoinHeapComparator>
    queue{JoinHeapComparator(&current)};
  for(size_type i = 0; i < readers.size(); i++)
  {
    if(readers[i]->size() > 0) { queue.push(i); }
  }
  std::string output_name = TempFile::getName("gcsa_join_run");
  JoinFileWriter output(output_name, logical, kind, total_records, codec, stats);
  while(!queue.empty())
  {
    size_type best = queue.top(); queue.pop();
    output.writeRecord(current[best]); offsets[best]++;
    if(offsets[best] < readers[best]->size())
    {
      readers[best]->read(offsets[best], current[best]); queue.push(best);
    }
  }
  if(stats != nullptr)
  {
    stats->merge_operations++;
    stats->max_records_resident = std::max(stats->max_records_resident,
      2 * static_cast<size_type>(inputs.size()));
    stats->max_bytes_resident = std::max(stats->max_bytes_resident,
      merge_reservation);
  }
  return output.finish();
}

class ExternalJoinSorter
{
public:
  ExternalJoinSorter(logical_file_id_t logical, JoinKeyKind kind,
    size_type byte_budget, size_type requested_fan_in, bool verify_payloads,
    ExternalPathJoinStats* stats, const TempFileCodecParameters& run_codec) :
    logical_id(logical), key_kind(kind), budget(byte_budget), fan_in(0),
    run_records(0), verify_runs(verify_payloads), statistics(stats),
    codec(run_codec)
  {
    if(this->budget < joinRunSortMinimumBudget(this->codec))
    {
      throw joinError("memory budget is too small for the configured join-run codec");
    }
    size_type writer_extra = joinRunWriterExtra(this->codec);
    size_type per_input = checkedJoinAdd(joinRawReaderReservation(),
      joinRunReaderMemory(this->codec), "compressed join merge reader");
    size_type merge_fixed = checkedJoinAdd(JOIN_FIXED_BYTES, writer_extra,
      "compressed join merge writer");
    size_type maximum_fan_in = (this->budget - merge_fixed) / per_input;
    this->fan_in = std::min(std::max(static_cast<size_type>(2), requested_fan_in),
      maximum_fan_in);
    if(this->fan_in < 2) { throw joinError("memory budget cannot support a two-way join merge"); }
    // Sorting and encoding are separate phases. Preserve the current in-place
    // sort capacity whenever the codec fits its existing slack, while also
    // proving that the sorted vector and framed writer can coexist during
    // emission.
    size_type sort_capacity = (this->budget - JOIN_FIXED_BYTES) /
      (2 * sizeof(JoinRecord));
    size_type emit_capacity = (this->budget - merge_fixed) /
      sizeof(JoinRecord);
    this->run_records = std::max(static_cast<size_type>(1),
      std::min(sort_capacity, emit_capacity));
    this->buffer.reserve(this->run_records);
    if(this->statistics != nullptr)
    {
      this->statistics->max_records_resident = std::max(
        this->statistics->max_records_resident, this->run_records);
      this->statistics->max_bytes_resident = std::max(
        this->statistics->max_bytes_resident,
        std::max(JOIN_FIXED_BYTES + 2 * this->run_records * sizeof(JoinRecord),
          merge_fixed + this->run_records * sizeof(JoinRecord)));
    }
  }

  void add(const JoinRecord& record)
  {
    this->buffer.push_back(record);
    if(this->buffer.size() == this->run_records) { this->flush(); }
  }

  JoinRun finish()
  {
    this->flush();
    std::vector<JoinRecord>().swap(this->buffer);
    std::vector<JoinRun> remaining;
    for(size_type level = 0; level < this->levels.size(); level++)
    {
      remaining.insert(remaining.end(), this->levels[level].begin(), this->levels[level].end());
      this->levels[level].clear();
    }
    if(remaining.empty())
    {
      std::string empty_name = TempFile::getName("gcsa_join_run");
      JoinFileWriter empty(empty_name, this->logical_id, this->key_kind, 0,
        this->codec, this->statistics);
      return empty.finish();
    }
    while(remaining.size() > 1)
    {
      std::vector<JoinRun> next;
      for(size_type first = 0; first < remaining.size(); first += this->fan_in)
      {
        size_type last = std::min(remaining.size(), first + this->fan_in);
        if(last - first == 1) { next.push_back(remaining[first]); continue; }
        std::vector<JoinRun> group(remaining.begin() + first, remaining.begin() + last);
        JoinRun merged = mergeJoinRuns(group, this->logical_id,
          this->key_kind, this->budget, this->codec,
          this->verify_runs, this->statistics);
        for(size_type i = 0; i < group.size(); i++) { removeJoinRun(group[i]); }
        next.push_back(merged);
      }
      remaining.swap(next);
    }
    return remaining.front();
  }

private:
  void flush()
  {
    if(this->buffer.empty()) { return; }
    // Sorting is in-place inside the run's existing byte reservation. Large
    // runs use the construction thread pool; tiny runs stay sequential to
    // avoid making forced-spill tests and small indexes pay team startup cost.
    if(this->buffer.size() >= JOIN_PARALLEL_SORT_MIN_RECORDS && omp_get_max_threads() > 1)
    {
      parallelQuickSort(this->buffer.begin(), this->buffer.end(), joinRecordLess);
      if(this->statistics != nullptr) { this->statistics->join_parallel_sorts++; }
    }
    else
    {
      sequentialSort(this->buffer.begin(), this->buffer.end(), joinRecordLess);
    }
    std::string name = TempFile::getName("gcsa_join_run");
    JoinFileWriter writer(name, this->logical_id, this->key_kind,
      this->buffer.size(), this->codec, this->statistics);
    for(size_type i = 0; i < this->buffer.size(); i++) { writer.writeRecord(this->buffer[i]); }
    JoinRun run = writer.finish();
    // Release the sort allocation before a leveled merge reserves its reader
    // state. Keeping vector capacity here would make the two phase budgets add.
    std::vector<JoinRecord>().swap(this->buffer);
    if(this->statistics != nullptr) { this->statistics->initial_runs++; }
    this->addRun(run, 0);
    this->buffer.reserve(this->run_records);
  }

  void addRun(const JoinRun& run, size_type level)
  {
    if(this->levels.size() <= level) { this->levels.resize(level + 1); }
    this->levels[level].push_back(run);
    if(this->levels[level].size() < this->fan_in) { return; }
    std::vector<JoinRun> inputs; inputs.swap(this->levels[level]);
    JoinRun merged = mergeJoinRuns(inputs, this->logical_id,
      this->key_kind, this->budget, this->codec,
      this->verify_runs, this->statistics);
    for(size_type i = 0; i < inputs.size(); i++) { removeJoinRun(inputs[i]); }
    this->addRun(merged, level + 1);
  }

  logical_file_id_t logical_id;
  JoinKeyKind key_kind;
  size_type budget, fan_in, run_records;
  bool verify_runs;
  ExternalPathJoinStats* statistics;
  TempFileCodecParameters codec;
  std::vector<JoinRecord> buffer;
  std::vector<std::vector<JoinRun>> levels;
};

class PathShardReader
{
public:
  PathShardReader(const PathGraph& graph, size_type file) :
    path_name(graph.path_names.at(file)), rank_name(graph.rank_names.at(file)),
    path_descriptor(-1), rank_descriptor(-1), path_count(graph.path_counts.at(file)),
    rank_count(graph.rank_counts.at(file)), path_offset(0), rank_offset(0),
    path_cache_released(0), rank_cache_released(0),
    path_buffer(std::max(static_cast<size_type>(1), JOIN_IO_BUFFER_BYTES / sizeof(PathNode))),
    rank_buffer(std::max(static_cast<size_type>(1),
      JOIN_IO_BUFFER_BYTES / sizeof(PathNode::rank_type))),
    path_buffer_first(0), path_buffer_records(0),
    rank_buffer_first(0), rank_buffer_records(0), compressed_path(), compressed_rank()
  {
    const bool path_framed = CompressedBlockReader::isFramed(this->path_name);
    const bool rank_framed = CompressedBlockReader::isFramed(this->rank_name);
    if(path_framed != rank_framed)
    {
      throw joinError("path/rank shard uses inconsistent storage formats",
        this->path_name);
    }
    if(path_framed)
    {
      this->compressed_path.reset(new CompressedBlockReader(this->path_name));
      this->compressed_rank.reset(new CompressedBlockReader(this->rank_name));
      if(this->compressed_path->logicalSize() !=
           checkedJoinMultiply(this->path_count, sizeof(PathNode), "path shard logical bytes") ||
         this->compressed_rank->logicalSize() !=
           checkedJoinMultiply(this->rank_count, sizeof(PathNode::rank_type), "rank shard logical bytes"))
      {
        throw joinError("compressed path shard length does not match metadata",
          this->path_name);
      }
      return;
    }

    this->path_descriptor = ::open(this->path_name.c_str(), O_RDONLY);
    this->rank_descriptor = ::open(this->rank_name.c_str(), O_RDONLY);
    if(this->path_descriptor < 0 || this->rank_descriptor < 0)
    {
      throw joinError("cannot open path shard", this->path_name);
    }
    adviseSequential(this->path_descriptor); adviseSequential(this->rank_descriptor);
    struct stat paths, ranks;
    if(::fstat(this->path_descriptor, &paths) != 0 || ::fstat(this->rank_descriptor, &ranks) != 0 ||
       static_cast<size_type>(paths.st_size) != this->path_count * sizeof(PathNode) ||
       static_cast<size_type>(ranks.st_size) != this->rank_count * sizeof(PathNode::rank_type))
    {
      throw joinError("path shard length does not match metadata", this->path_name);
    }
  }

  ~PathShardReader()
  {
    if(this->path_descriptor >= 0)
    {
      discardCachedRange(this->path_descriptor, 0, 0); ::close(this->path_descriptor);
    }
    if(this->rank_descriptor >= 0)
    {
      discardCachedRange(this->rank_descriptor, 0, 0); ::close(this->rank_descriptor);
    }
  }

  bool read(JoinRecord& record)
  {
    if(this->path_offset >= this->path_count)
    {
      if(this->rank_offset != this->rank_count)
      {
        throw joinError("path shard rank count mismatch", this->rank_name);
      }
      return false;
    }
    if(this->path_buffer_records == 0 || this->path_offset < this->path_buffer_first ||
       this->path_offset >= this->path_buffer_first + this->path_buffer_records)
    {
      this->refillPaths();
    }
    record.node = this->path_buffer[this->path_offset - this->path_buffer_first];
    // Pruning may collapse a unique equal-label range to order 0. Such a
    // sorted path still owns one boundary rank (ranks() == order() + 1), can
    // bypass left expansion, and remains a valid right-side join target.
    if(record.node.pointer() != this->rank_offset ||
       record.node.order() > PathLabel::LABEL_LENGTH || record.node.lcp() > record.node.order() ||
       record.node.ranks() > this->rank_count - this->rank_offset)
    {
      throw joinError("invalid path or rank pointer", this->path_name);
    }
    this->readRanks(record.labels, record.node.ranks());
    for(size_type i = record.node.ranks(); i < JOIN_LABEL_COUNT; i++)
    {
      record.labels[i] = PathLabel::NO_RANK;
    }
    this->path_offset++;
    record.node.setPointer(0);
    return true;
  }

private:
  void refillPaths()
  {
    this->path_buffer_first = this->path_offset;
    this->path_buffer_records = std::min(this->path_buffer.size(),
      this->path_count - this->path_buffer_first);
    const size_type bytes = this->path_buffer_records * sizeof(PathNode);
    if(this->compressed_path)
    {
      if(this->compressed_path->read(this->path_buffer.data(), bytes) != bytes)
      {
        throw joinError("short compressed path shard", this->path_name);
      }
    }
    else
    {
      preadAll(this->path_descriptor, this->path_buffer.data(), bytes,
        this->path_buffer_first * sizeof(PathNode), this->path_name);
    }
    if(this->path_descriptor >= 0)
    {
      trimReadCache(this->path_descriptor,
        (this->path_buffer_first + this->path_buffer_records) * sizeof(PathNode),
        this->path_cache_released);
    }
  }

  void refillRanks()
  {
    this->rank_buffer_first = this->rank_offset;
    this->rank_buffer_records = std::min(this->rank_buffer.size(),
      this->rank_count - this->rank_buffer_first);
    const size_type bytes = this->rank_buffer_records * sizeof(PathNode::rank_type);
    if(this->compressed_rank)
    {
      if(this->compressed_rank->read(this->rank_buffer.data(), bytes) != bytes)
      {
        throw joinError("short compressed rank shard", this->rank_name);
      }
    }
    else
    {
      preadAll(this->rank_descriptor, this->rank_buffer.data(), bytes,
        this->rank_buffer_first * sizeof(PathNode::rank_type), this->rank_name);
    }
    if(this->rank_descriptor >= 0)
    {
      trimReadCache(this->rank_descriptor,
        (this->rank_buffer_first + this->rank_buffer_records) * sizeof(PathNode::rank_type),
        this->rank_cache_released);
    }
  }

  void readRanks(PathNode::rank_type* target, size_type records)
  {
    while(records > 0)
    {
      if(this->rank_buffer_records == 0 || this->rank_offset < this->rank_buffer_first ||
         this->rank_offset >= this->rank_buffer_first + this->rank_buffer_records)
      {
        this->refillRanks();
      }
      size_type offset = this->rank_offset - this->rank_buffer_first;
      size_type copied = std::min(records, this->rank_buffer_records - offset);
      std::copy(this->rank_buffer.begin() + offset,
        this->rank_buffer.begin() + offset + copied, target);
      target += copied; records -= copied; this->rank_offset += copied;
    }
  }

  PathShardReader(const PathShardReader&);
  PathShardReader& operator=(const PathShardReader&);

  std::string path_name, rank_name;
  int path_descriptor, rank_descriptor;
  size_type path_count, rank_count, path_offset, rank_offset;
  off_t path_cache_released, rank_cache_released;
  std::vector<PathNode> path_buffer;
  std::vector<PathNode::rank_type> rank_buffer;
  size_type path_buffer_first, path_buffer_records;
  size_type rank_buffer_first, rank_buffer_records;
  std::unique_ptr<CompressedBlockReader> compressed_path, compressed_rank;
};

void
scanJoinSide(const PathGraph& graph, const std::vector<size_type>& shards,
  logical_file_id_t logical, JoinKeyKind kind, ExternalJoinSorter& sorter,
  ExternalPathJoinStats* stats)
{
  std::uint64_t ordinal = 0;
  for(size_type shard : shards)
  {
    PathShardReader reader(graph, shard);
    JoinRecord record;
    while(reader.read(record))
    {
      record.logical = logical; record.ordinal = ordinal++;
      if(kind == RIGHT_BY_FROM)
      {
        record.key = record.node.from; sorter.add(record);
        if(stats != nullptr) { stats->right_records++; }
      }
      else if(!record.node.sorted())
      {
        record.key = record.node.to; sorter.add(record);
        if(stats != nullptr) { stats->left_records++; }
      }
    }
  }
}

size_type
maximumCompressedShardReaderBytes(const PathGraph& graph)
{
  size_type maximum = 0;
  for(size_type file = 0; file < graph.files(); file++)
  {
    const bool path_framed = CompressedBlockReader::isFramed(
      graph.path_names[file]);
    const bool rank_framed = CompressedBlockReader::isFramed(
      graph.rank_names[file]);
    if(path_framed != rank_framed)
    {
      throw joinError("path/rank shard uses inconsistent storage formats",
        graph.path_names[file]);
    }
    if(!path_framed) { continue; }
    size_type pair = checkedJoinAdd(
      CompressedBlockReader::workingMemoryEstimate(
        CompressedBlockReader::declaredBlockSize(graph.path_names[file])),
      CompressedBlockReader::workingMemoryEstimate(
        CompressedBlockReader::declaredBlockSize(graph.rank_names[file])),
      "compressed source shard reader bytes");
    maximum = std::max(maximum, pair);
  }
  return maximum;
}

JoinRecord
extendRecord(const JoinRecord& left, const JoinRecord& right)
{
  JoinRecord result;
  result.logical = left.logical; result.key = 0; result.ordinal = 0;
  result.node.from = left.node.from; result.node.to = right.node.to;
  result.node.fields = 0;
  if(right.node.sorted()) { result.node.makeSorted(); }
  result.node.setPredecessors(left.node.predecessors());
  size_type order = left.node.order() + right.node.order();
  if(order > PathLabel::LABEL_LENGTH) { throw joinError("generated path label is too long"); }
  result.node.setOrder(order);
  result.node.setLCP(left.node.order() + right.node.lcp());
  result.node.setPointer(0);
  size_type offset = 0;
  for(size_type i = 0; i < left.node.order(); i++) { result.labels[offset++] = left.labels[i]; }
  for(size_type i = 0; i < right.node.ranks(); i++) { result.labels[offset++] = right.labels[i]; }
  for(; offset < JOIN_LABEL_COUNT; offset++) { result.labels[offset] = PathLabel::NO_RANK; }
  return result;
}

void
writeLabelRecord(const JoinRecord& record, ExternalPathSortSink& output,
  bool bypass, ExternalPathJoinStats* stats)
{
  output.write(record.node, record.labels);
  if(stats != nullptr)
  {
    stats->direct_label_records++;
    stats->intermediate_path_bytes_avoided += sizeof(PathNode) +
      record.node.ranks() * sizeof(PathNode::rank_type);
    if(bypass) { stats->sorted_bypass++; }
    else { stats->generated_records++; }
  }
}

void
joinSortedRuns(const JoinRun& left_run, const JoinRun& right_run,
  logical_file_id_t logical, size_type byte_budget, ExternalPathSortSink& output,
  bool verify_payloads, ExternalPathJoinStats* stats,
  size_type left_begin = 0, size_type left_end = std::numeric_limits<size_type>::max(),
  size_type right_begin = 0, size_type right_end = std::numeric_limits<size_type>::max(),
  bool emit_bypass = true)
{
  size_type reader_codec_bytes = joinRunReaderPairMemory(left_run, right_run);
  size_type resident_base = checkedJoinAdd(joinRawScanBaseBudget(),
    reader_codec_bytes, "compressed join scan base");
  size_type required_budget = checkedJoinAdd(resident_base,
    sizeof(JoinRecord), "compressed join scan minimum");
  if(byte_budget < required_budget)
  {
    throw joinError("join scan budget cannot admit both compressed run readers");
  }
  JoinFileReader left(left_run, logical, LEFT_BY_TO);
  JoinFileReader right(right_run, logical, RIGHT_BY_FROM);
  if(verify_payloads) { left.validate(); right.validate(); }
  left_end = std::min(left_end, left.size());
  right_end = std::min(right_end, right.size());
  if(left_begin > left_end || right_begin > right_end)
  {
    throw joinError("invalid worker join range");
  }
  if(stats != nullptr)
  {
    stats->max_bytes_resident = std::max(stats->max_bytes_resident, resident_base);
  }
  size_type left_offset = left_begin, right_offset = right_begin;
  JoinRecord left_record, right_record;
  while(left_offset < left_end && right_offset < right_end)
  {
    left.read(left_offset, left_record); right.read(right_offset, right_record);
    if(left_record.key < right_record.key) { left_offset++; continue; }
    if(right_record.key < left_record.key)
    {
      if(emit_bypass && right_record.node.sorted())
      {
        writeLabelRecord(right_record, output, true, stats);
      }
      right_offset++;
      continue;
    }
    node_type key = left_record.key;
    size_type left_limit = left_offset;
    do
    {
      left_limit++;
      if(left_limit >= left_end) { break; }
      left.read(left_limit, left_record);
    } while(left_record.key == key);
    size_type right_limit = right_offset;
    do
    {
      right_limit++;
      if(right_limit >= right_end) { break; }
      right.read(right_limit, right_record);
    } while(right_record.key == key);
    // Sorted paths no longer extend on the left, but they remain part of the
    // next generation and may still be right-hand extension targets. Emit the
    // bypass while this right-key group is already resident in the join scan.
    for(size_type j = right_offset; j < right_limit; j++)
    {
      right.read(j, right_record);
      if(emit_bypass && right_record.node.sorted())
      {
        writeLabelRecord(right_record, output, true, stats);
      }
    }
    size_type group_records = (left_limit - left_offset) + (right_limit - right_offset);
    if(stats != nullptr && group_records >
       std::max(static_cast<size_type>(1),
         (byte_budget - resident_base) / sizeof(JoinRecord)))
    {
      stats->blocked_key_groups++;
    }
    // Do not materialize either side of a pathological key. Keep a bounded
    // block from the smaller side and replay the other disk-backed range for
    // each block. This reduces the former per-record replay to per-block I/O;
    // the compact JoinRecord block retains all label data needed for extension.
    size_type left_count = left_limit - left_offset;
    size_type right_count = right_limit - right_offset;
    size_type block_capacity = std::max(static_cast<size_type>(1),
      (byte_budget - resident_base) / sizeof(JoinRecord));
    // A pathological key must never become a single in-memory side just
    // because one side happens to fit the reservation. Keep at least two
    // deterministic blocks whenever that side has multiple records.
    size_type smaller_side = std::min(left_count, right_count);
    if(smaller_side > 1) { block_capacity = std::min(block_capacity, smaller_side / 2); }
    std::vector<JoinRecord> block;
    block.reserve(block_capacity);
    if(stats != nullptr)
    {
      stats->max_records_resident = std::max(stats->max_records_resident, block_capacity);
      stats->max_bytes_resident = std::max(stats->max_bytes_resident,
        checkedJoinAdd(resident_base,
          checkedJoinMultiply(block_capacity, sizeof(JoinRecord), "join block bytes"),
          "join block reservation"));
    }
    if(left_count <= right_count)
    {
      for(size_type begin = left_offset; begin < left_limit; begin += block_capacity)
      {
        if(stats != nullptr) { stats->blocked_key_blocks++; }
        size_type end = std::min(left_limit, begin + block_capacity);
        block.clear();
        for(size_type i = begin; i < end; i++) { left.read(i, left_record); block.push_back(left_record); }
        for(const JoinRecord& cached_left : block)
        {
          for(size_type j = right_offset; j < right_limit; j++)
          {
            right.read(j, right_record);
            writeLabelRecord(extendRecord(cached_left, right_record), output, false, stats);
          }
        }
      }
    }
    else
    {
      for(size_type i = left_offset; i < left_limit; i++)
      {
        left.read(i, left_record);
        for(size_type begin = right_offset; begin < right_limit; begin += block_capacity)
        {
          if(stats != nullptr) { stats->blocked_key_blocks++; }
          size_type end = std::min(right_limit, begin + block_capacity);
          block.clear();
          for(size_type j = begin; j < end; j++) { right.read(j, right_record); block.push_back(right_record); }
          for(const JoinRecord& cached_right : block)
          {
            writeLabelRecord(extendRecord(left_record, cached_right), output, false, stats);
          }
        }
      }
    }
    left_offset = left_limit; right_offset = right_limit;
  }
  // The left stream can end before the right stream. Such right records cannot
  // participate in joins, but sorted records must still bypass expansion.
  while(right_offset < right_end)
  {
    right.read(right_offset++, right_record);
    if(emit_bypass && right_record.node.sorted())
    {
      writeLabelRecord(right_record, output, true, stats);
    }
  }
}

struct JoinPartition
{
  size_type left_begin, left_end, right_begin, right_end;
  size_type expected_paths, expected_ranks, expected_bytes;
  size_type generated_paths, bypass_paths;
  bool emit_bypass;

  JoinPartition() :
    left_begin(0), left_end(0), right_begin(0), right_end(0),
    expected_paths(0), expected_ranks(0), expected_bytes(0),
    generated_paths(0), bypass_paths(0),
    emit_bypass(true) { }
};

// Planning must not retain one entry per distinct key: a realistic graph can
// have many more keys than fit in the join reservation. This fixed-size,
// deterministic systematic sample builds a bounded MSD radix range-pack tree.
// The subsequent exact group pass remains the authority for ranges and output
// counts, so a missed or underestimated key can only affect performance.
constexpr size_type JOIN_PLAN_SAMPLE_RECORDS = 4096;
constexpr size_type JOIN_PLAN_RADIX_BITS = 4;
constexpr size_type JOIN_PLAN_MAX_BYTES = 4 * MEGABYTE;
constexpr size_type JOIN_PLAN_PACK_BYTES = 4 * sizeof(std::uint64_t);
constexpr size_type JOIN_PLAN_FIXED_BYTES = 8 + 4 + 4 + 18 * sizeof(std::uint64_t);
constexpr size_type JOIN_PLAN_RUNTIME_OVERHEAD = 6 * JOIN_IO_BUFFER_BYTES +
  2 * JOIN_PLAN_SAMPLE_RECORDS * sizeof(node_type) + 64 * KILOBYTE;
constexpr std::uint64_t JOIN_PLAN_MAGIC = 0x334e4c5044534347ULL; // "GCSDPLN3"
constexpr std::uint32_t JOIN_PLAN_VERSION = 3;

struct JoinKeySample
{
  std::vector<node_type> keys;
};

// A pack is one leaf in a deterministic MSD radix tree. Prefix contains the
// high `bits` bits of the join key; estimates are scheduling hints only. The
// leaves tile the complete uint64_t key space in increasing order, including
// empty sampled ranges, so an unsampled key still has an unambiguous pack.
struct JoinRadixPack
{
  node_type prefix;
  size_type bits, estimated_input_bytes, estimated_output_bytes;

  JoinRadixPack(node_type prefix_value = 0, size_type prefix_bits = 0,
    size_type input_bytes = 0, size_type output_bytes = 0) :
    prefix(prefix_value), bits(prefix_bits),
    estimated_input_bytes(input_bytes), estimated_output_bytes(output_bytes) { }
};

struct JoinRadixPlan
{
  std::vector<JoinRadixPack> packs;
  std::vector<node_type> boundaries;
  std::vector<node_type> heavy_keys;
  size_type sampled_records, splits, max_bits, serialized_limit;
  size_type left_groups, right_groups;
  std::uint64_t left_summary_checksum, right_summary_checksum;
  std::uint64_t left_detail_checksum, right_detail_checksum;
  bool capped;

  JoinRadixPlan() : sampled_records(0), splits(0), max_bits(0),
    serialized_limit(0), left_groups(0), right_groups(0),
    left_summary_checksum(0), right_summary_checksum(0),
    left_detail_checksum(0), right_detail_checksum(0), capped(false) { }
};

void
bindJoinRadixPlanSidecars(JoinRadixPlan& plan, const JoinRun& left,
  const JoinRun& right)
{
  if(!left.sidecar.complete() || !right.sidecar.complete())
  {
    throw joinError("cannot bind a radix plan without join group sidecars");
  }
  plan.left_groups = left.sidecar.groups; plan.right_groups = right.sidecar.groups;
  plan.left_summary_checksum = left.sidecar.summary_checksum;
  plan.right_summary_checksum = right.sidecar.summary_checksum;
  plan.left_detail_checksum = left.sidecar.detail_checksum;
  plan.right_detail_checksum = right.sidecar.detail_checksum;
}

bool
joinRadixPlanSidecarsMatch(const JoinRadixPlan& plan, const JoinRun& left,
  const JoinRun& right)
{
  return left.sidecar.complete() && right.sidecar.complete() &&
    plan.left_groups == left.sidecar.groups &&
    plan.right_groups == right.sidecar.groups &&
    plan.left_summary_checksum == left.sidecar.summary_checksum &&
    plan.right_summary_checksum == right.sidecar.summary_checksum &&
    plan.left_detail_checksum == left.sidecar.detail_checksum &&
    plan.right_detail_checksum == right.sidecar.detail_checksum;
}

JoinKeySample
sampleJoinKeys(const JoinFileReader& reader)
{
  JoinKeySample result;
  size_type samples = std::min(reader.size(), JOIN_PLAN_SAMPLE_RECORDS);
  if(samples == 0) { return result; }
  result.keys.reserve(samples);
  JoinRecord record;
  for(size_type sample = 0; sample < samples; sample++)
  {
    // Avoid overflowing sample * reader.size() for a multi-terabyte run.
    size_type quotient = reader.size() / samples;
    size_type remainder = reader.size() % samples;
    size_type offset = sample * quotient + (sample * remainder) / samples;
    reader.read(offset, record);
    result.keys.push_back(record.key);
  }
  return result;
}

size_type
sampleKeyCount(const JoinKeySample& sample, node_type key)
{
  auto range = std::equal_range(sample.keys.begin(), sample.keys.end(), key);
  return static_cast<size_type>(range.second - range.first);
}

size_type
scaleSampleCount(size_type count, size_type total, size_type samples)
{
  if(count == 0 || total == 0) { return 0; }
  // ceil(count * total / samples), without overflow.
  size_type quotient = total / samples, remainder = total % samples;
  size_type remainder_product = checkedJoinMultiply(count, remainder,
    "sampled key remainder");
  return checkedJoinAdd(checkedJoinMultiply(count, quotient,
    "sampled key multiplicity"),
    remainder_product / samples + (remainder_product % samples != 0),
    "sampled key estimate");
}

size_type
saturatingPlanAdd(size_type left, size_type right)
{
  if(left > std::numeric_limits<size_type>::max() - right)
  {
    return std::numeric_limits<size_type>::max();
  }
  return left + right;
}

size_type
saturatingPlanMultiply(size_type left, size_type right)
{
  if(left != 0 && right > std::numeric_limits<size_type>::max() / left)
  {
    return std::numeric_limits<size_type>::max();
  }
  return left * right;
}

bool
keyHasPrefix(node_type key, node_type prefix, size_type bits)
{
  return (bits == 0 || key >> (64 - bits) == prefix);
}

size_type
estimatedPrefixOutput(const JoinKeySample& left_sample,
  const JoinKeySample& right_sample, size_type left_total,
  size_type right_total, node_type prefix, size_type bits)
{
  const size_type minimum_output_bytes = sizeof(PathNode) +
    3 * sizeof(PathNode::rank_type);
  size_type result = 0;
  node_type lower = (bits == 0 ? node_type(0) : prefix << (64 - bits));
  auto left = std::lower_bound(left_sample.keys.begin(), left_sample.keys.end(), lower);
  while(left != left_sample.keys.end() && keyHasPrefix(*left, prefix, bits))
  {
    node_type key = *left;
    auto left_end = std::upper_bound(left, left_sample.keys.end(), key);
    size_type right_count = sampleKeyCount(right_sample, key);
    if(right_count == 0) { left = left_end; continue; }
    size_type estimated_left = scaleSampleCount(
      static_cast<size_type>(left_end - left), left_total,
      left_sample.keys.size());
    size_type estimated_right = scaleSampleCount(right_count, right_total,
      right_sample.keys.size());
    size_type paths = saturatingPlanMultiply(estimated_left, estimated_right);
    result = saturatingPlanAdd(result,
      saturatingPlanMultiply(paths, minimum_output_bytes));
    left = left_end;
  }
  return result;
}

size_type
joinPlanBytes(size_type packs, size_type heavy_keys)
{
  size_type pack_bytes = saturatingPlanMultiply(packs, JOIN_PLAN_PACK_BYTES);
  size_type heavy_bytes = saturatingPlanMultiply(heavy_keys, sizeof(node_type));
  return saturatingPlanAdd(JOIN_PLAN_FIXED_BYTES,
    saturatingPlanAdd(pack_bytes, heavy_bytes));
}

bool
joinPlanCanFit(const JoinRadixPlan& plan, size_type extra_packs,
  size_type extra_heavy_keys)
{
  return (joinPlanBytes(saturatingPlanAdd(plan.packs.size(), extra_packs),
    saturatingPlanAdd(plan.heavy_keys.size(), extra_heavy_keys)) <=
    plan.serialized_limit);
}

node_type
radixPackLower(const JoinRadixPack& pack)
{
  return (pack.bits == 0 ? node_type(0) :
    (pack.bits == 64 ? pack.prefix : pack.prefix << (64 - pack.bits)));
}

void
finalizeJoinRadixPlan(JoinRadixPlan& plan)
{
  plan.boundaries.clear();
  plan.boundaries.reserve(plan.packs.empty() ? 0 : plan.packs.size() - 1);
  for(size_type i = 1; i < plan.packs.size(); i++)
  {
    plan.boundaries.push_back(radixPackLower(plan.packs[i]));
  }
  std::sort(plan.heavy_keys.begin(), plan.heavy_keys.end());
  plan.heavy_keys.erase(std::unique(plan.heavy_keys.begin(), plan.heavy_keys.end()),
    plan.heavy_keys.end());
}

void
planRadixPrefix(const JoinKeySample& left_sample,
  const JoinKeySample& right_sample, size_type left_total,
  size_type right_total, node_type prefix, size_type bits,
  size_type target_bytes, size_type reserved_packs, JoinRadixPlan& result)
{
  auto left_begin = std::lower_bound(left_sample.keys.begin(),
    left_sample.keys.end(), (bits == 0 ? node_type(0) : prefix << (64 - bits)));
  auto right_begin = std::lower_bound(right_sample.keys.begin(),
    right_sample.keys.end(), (bits == 0 ? node_type(0) : prefix << (64 - bits)));
  auto prefix_end = [prefix, bits](node_type key)
  {
    return !keyHasPrefix(key, prefix, bits);
  };
  auto left_end = std::find_if(left_begin, left_sample.keys.end(), prefix_end);
  auto right_end = std::find_if(right_begin, right_sample.keys.end(), prefix_end);
  size_type left_hits = static_cast<size_type>(left_end - left_begin);
  size_type right_hits = static_cast<size_type>(right_end - right_begin);

  size_type estimated_left = scaleSampleCount(left_hits, left_total,
    left_sample.keys.size());
  size_type estimated_right = scaleSampleCount(right_hits, right_total,
    right_sample.keys.size());
  size_type input_bytes = saturatingPlanMultiply(
    saturatingPlanAdd(estimated_left, estimated_right), JOIN_RECORD_BYTES);
  size_type output_bytes = estimatedPrefixOutput(left_sample, right_sample,
    left_total, right_total, prefix, bits);
  size_type estimated_bytes = std::max(input_bytes, output_bytes);

  bool empty = (left_hits == 0 && right_hits == 0);
  node_type minimum_key = std::numeric_limits<node_type>::max();
  node_type maximum_key = 0;
  if(left_hits > 0)
  {
    minimum_key = std::min(minimum_key, *left_begin);
    maximum_key = std::max(maximum_key, *(left_end - 1));
  }
  if(right_hits > 0)
  {
    minimum_key = std::min(minimum_key, *right_begin);
    maximum_key = std::max(maximum_key, *(right_end - 1));
  }
  bool terminal = (empty || estimated_bytes <= target_bytes || bits == 64 ||
    minimum_key == maximum_key);
  // A split replaces this prefix with at least 16 child packs. Keep enough
  // space for every later sibling before allowing deeper recursion. When the
  // byte cap is reached, the exact pass simply partitions the coarser pack.
  if(!terminal && !joinPlanCanFit(result,
      reserved_packs + (static_cast<size_type>(1) << JOIN_PLAN_RADIX_BITS), 0))
  {
    result.capped = true; terminal = true;
  }
  if(terminal)
  {
    if(!joinPlanCanFit(result, reserved_packs + 1, 0))
    {
      throw joinError("MSD range-pack planner violated its byte reservation");
    }
    result.packs.emplace_back(prefix, bits, input_bytes, output_bytes);
    if(!empty && output_bytes > target_bytes &&
       sampleKeyCount(left_sample, minimum_key) > 0 &&
       sampleKeyCount(right_sample, minimum_key) > 0)
    {
      if(joinPlanCanFit(result, reserved_packs, 1))
      {
        result.heavy_keys.push_back(minimum_key);
      }
      else { result.capped = true; }
    }
    result.max_bits = std::max(result.max_bits, bits);
    return;
  }

  size_type next_bits = bits + JOIN_PLAN_RADIX_BITS;
  result.splits++;
  result.max_bits = std::max(result.max_bits, next_bits);
  for(size_type digit = 0; digit < (static_cast<size_type>(1) << JOIN_PLAN_RADIX_BITS); digit++)
  {
    node_type child = (prefix << JOIN_PLAN_RADIX_BITS) | digit;
    planRadixPrefix(left_sample, right_sample, left_total, right_total,
      child, next_bits, target_bytes,
      reserved_packs + ((static_cast<size_type>(1) << JOIN_PLAN_RADIX_BITS) -
        digit - 1), result);
  }
}

JoinRadixPlan
sampleJoinRadixPlan(const JoinFileReader& left, const JoinFileReader& right,
  size_type target_bytes, size_type serialized_limit)
{
  JoinKeySample left_sample = sampleJoinKeys(left);
  JoinKeySample right_sample = sampleJoinKeys(right);
  JoinRadixPlan result;
  result.serialized_limit = serialized_limit;
  result.sampled_records = left_sample.keys.size() + right_sample.keys.size();
  planRadixPrefix(left_sample, right_sample, left.size(), right.size(),
    0, 0, target_bytes, 0, result);
  finalizeJoinRadixPlan(result);
  return result;
}

size_type
joinPlanSerializedLimit(size_type memory_budget, size_type reader_codec_bytes)
{
  const size_type minimum = JOIN_PLAN_FIXED_BYTES + JOIN_PLAN_PACK_BYTES;
  size_type runtime_overhead = checkedJoinAdd(JOIN_PLAN_RUNTIME_OVERHEAD,
    reader_codec_bytes, "compressed join planner overhead");
  if(memory_budget <= runtime_overhead + 4 * minimum)
  {
    return minimum;
  }
  // On fresh planning, vector capacity can approach twice the final pack
  // count; on restore, the payload, decoded packs, and derived boundaries
  // overlap. One quarter covers the worst of those representations.
  return std::min(JOIN_PLAN_MAX_BYTES,
    (memory_budget - runtime_overhead) / 4);
}

template<class Value>
void
appendPlanValue(std::vector<std::uint8_t>& target, Value value)
{
  size_type offset = target.size();
  if(offset > JOIN_PLAN_MAX_BYTES || sizeof(Value) > JOIN_PLAN_MAX_BYTES - offset)
  {
    throw joinError("MSD range-pack plan exceeds the format limit");
  }
  target.resize(offset + sizeof(Value));
  std::uint8_t* output = target.data() + offset;
  encodeLittle<Value>(output, value);
}

template<class Value>
Value
readPlanValue(const std::vector<std::uint8_t>& source, size_type& offset)
{
  if(offset > source.size() || sizeof(Value) > source.size() - offset)
  {
    throw joinError("truncated MSD range-pack plan");
  }
  const std::uint8_t* input = source.data() + offset;
  Value result = decodeLittle<Value>(input); offset += sizeof(Value);
  return result;
}

std::string
joinRadixPlanTaskName(const std::string& generation,
  logical_file_id_t logical, size_type target_bytes, size_type serialized_limit,
  std::uint64_t left_checksum, std::uint64_t right_checksum)
{
  return generation + "-msd-plan-l" + std::to_string(logical.value) +
    "-v3" +
    "-t" + std::to_string(target_bytes) +
    "-m" + std::to_string(serialized_limit) +
    "-lc" + std::to_string(left_checksum) +
    "-rc" + std::to_string(right_checksum);
}

ArtifactIdentity
joinRadixPlanArtifact(const std::string& task)
{
  return ArtifactIdentity(task, "join-plan", "msd-range-packs",
    "join-msd-plan-v3");
}

void
validateJoinRadixPacks(const JoinRadixPlan& plan)
{
  if(plan.packs.empty()) { throw joinError("MSD range-pack plan is empty"); }
  typedef unsigned __int128 wide_type;
  wide_type cursor = 0;
  const wide_type universe = (static_cast<wide_type>(1) << 64);
  for(const JoinRadixPack& pack : plan.packs)
  {
    if(pack.bits > 64 || (pack.bits == 0 && pack.prefix != 0) ||
       (pack.bits < 64 && pack.bits > 0 &&
        pack.prefix >= (static_cast<node_type>(1) << pack.bits)))
    {
      throw joinError("invalid MSD range-pack prefix");
    }
    wide_type lower = static_cast<wide_type>(radixPackLower(pack));
    wide_type width = (static_cast<wide_type>(1) << (64 - pack.bits));
    if(lower != cursor || lower + width > universe)
    {
      throw joinError("MSD range packs do not tile the join-key space");
    }
    cursor = lower + width;
  }
  if(cursor != universe)
  {
    throw joinError("MSD range packs do not cover the join-key space");
  }
}

std::vector<std::uint8_t>
encodeJoinRadixPlan(const JoinRadixPlan& plan, logical_file_id_t logical,
  size_type target_bytes, size_type left_records, size_type right_records,
  std::uint64_t left_checksum, std::uint64_t right_checksum)
{
  validateJoinRadixPacks(plan);
  std::vector<std::uint8_t> result;
  size_type serialized_bytes = joinPlanBytes(plan.packs.size(),
    plan.heavy_keys.size());
  if(serialized_bytes > plan.serialized_limit ||
     serialized_bytes > JOIN_PLAN_MAX_BYTES)
  {
    throw joinError("MSD range-pack plan exceeds its memory reservation");
  }
  result.reserve(serialized_bytes);
  appendPlanValue<std::uint64_t>(result, JOIN_PLAN_MAGIC);
  appendPlanValue<std::uint32_t>(result, JOIN_PLAN_VERSION);
  appendPlanValue<std::uint32_t>(result, logical.value);
  appendPlanValue<std::uint64_t>(result, target_bytes);
  appendPlanValue<std::uint64_t>(result, left_records);
  appendPlanValue<std::uint64_t>(result, right_records);
  appendPlanValue<std::uint64_t>(result, left_checksum);
  appendPlanValue<std::uint64_t>(result, right_checksum);
  appendPlanValue<std::uint64_t>(result, plan.left_groups);
  appendPlanValue<std::uint64_t>(result, plan.right_groups);
  appendPlanValue<std::uint64_t>(result, plan.left_summary_checksum);
  appendPlanValue<std::uint64_t>(result, plan.right_summary_checksum);
  appendPlanValue<std::uint64_t>(result, plan.left_detail_checksum);
  appendPlanValue<std::uint64_t>(result, plan.right_detail_checksum);
  appendPlanValue<std::uint64_t>(result, plan.serialized_limit);
  appendPlanValue<std::uint64_t>(result, plan.sampled_records);
  appendPlanValue<std::uint64_t>(result, plan.splits);
  appendPlanValue<std::uint64_t>(result, plan.max_bits);
  appendPlanValue<std::uint64_t>(result, plan.capped ? 1 : 0);
  appendPlanValue<std::uint64_t>(result, plan.packs.size());
  appendPlanValue<std::uint64_t>(result, plan.heavy_keys.size());
  for(const JoinRadixPack& pack : plan.packs)
  {
    appendPlanValue<node_type>(result, pack.prefix);
    appendPlanValue<std::uint64_t>(result, pack.bits);
    appendPlanValue<std::uint64_t>(result, pack.estimated_input_bytes);
    appendPlanValue<std::uint64_t>(result, pack.estimated_output_bytes);
  }
  for(node_type key : plan.heavy_keys)
  {
    appendPlanValue<node_type>(result, key);
  }
  if(result.size() != serialized_bytes)
  {
    throw joinError("MSD range-pack serialized length mismatch");
  }
  return result;
}

JoinRadixPlan
decodeJoinRadixPlan(const std::vector<std::uint8_t>& payload,
  logical_file_id_t logical, size_type target_bytes,
  size_type left_records, size_type right_records,
  std::uint64_t left_checksum, std::uint64_t right_checksum,
  size_type serialized_limit, const JoinRun& left_run, const JoinRun& right_run)
{
  size_type offset = 0;
  if(readPlanValue<std::uint64_t>(payload, offset) != JOIN_PLAN_MAGIC ||
     readPlanValue<std::uint32_t>(payload, offset) != JOIN_PLAN_VERSION ||
     readPlanValue<std::uint32_t>(payload, offset) != logical.value ||
     readPlanValue<std::uint64_t>(payload, offset) != target_bytes ||
     readPlanValue<std::uint64_t>(payload, offset) != left_records ||
     readPlanValue<std::uint64_t>(payload, offset) != right_records ||
     readPlanValue<std::uint64_t>(payload, offset) != left_checksum ||
     readPlanValue<std::uint64_t>(payload, offset) != right_checksum)
  {
    throw joinError("MSD range-pack plan does not match its join inputs");
  }
  JoinRadixPlan result;
  result.left_groups = readPlanValue<std::uint64_t>(payload, offset);
  result.right_groups = readPlanValue<std::uint64_t>(payload, offset);
  result.left_summary_checksum = readPlanValue<std::uint64_t>(payload, offset);
  result.right_summary_checksum = readPlanValue<std::uint64_t>(payload, offset);
  result.left_detail_checksum = readPlanValue<std::uint64_t>(payload, offset);
  result.right_detail_checksum = readPlanValue<std::uint64_t>(payload, offset);
  if(!joinRadixPlanSidecarsMatch(result, left_run, right_run) ||
     readPlanValue<std::uint64_t>(payload, offset) != serialized_limit)
  {
    throw joinError("MSD range-pack plan does not match its join sidecars");
  }
  result.sampled_records = readPlanValue<std::uint64_t>(payload, offset);
  result.splits = readPlanValue<std::uint64_t>(payload, offset);
  result.max_bits = readPlanValue<std::uint64_t>(payload, offset);
  result.capped = (readPlanValue<std::uint64_t>(payload, offset) != 0);
  size_type packs = readPlanValue<std::uint64_t>(payload, offset);
  size_type heavy_keys = readPlanValue<std::uint64_t>(payload, offset);
  result.serialized_limit = serialized_limit;
  if(packs == 0 || packs > (payload.size() - offset) / JOIN_PLAN_PACK_BYTES)
  {
    throw joinError("invalid MSD range-pack count");
  }
  result.packs.reserve(packs);
  for(size_type i = 0; i < packs; i++)
  {
    node_type prefix = readPlanValue<node_type>(payload, offset);
    size_type bits = readPlanValue<std::uint64_t>(payload, offset);
    size_type input_bytes = readPlanValue<std::uint64_t>(payload, offset);
    size_type output_bytes = readPlanValue<std::uint64_t>(payload, offset);
    result.packs.emplace_back(prefix, bits, input_bytes, output_bytes);
  }
  if(heavy_keys > (payload.size() - offset) / sizeof(node_type) ||
     offset + heavy_keys * sizeof(node_type) != payload.size())
  {
    throw joinError("invalid MSD range-pack heavy-key count");
  }
  result.heavy_keys.reserve(heavy_keys);
  for(size_type i = 0; i < heavy_keys; i++)
  {
    result.heavy_keys.push_back(readPlanValue<node_type>(payload, offset));
  }
  if(joinPlanBytes(result.packs.size(), result.heavy_keys.size()) != payload.size() ||
     payload.size() > serialized_limit ||
     !std::is_sorted(result.heavy_keys.begin(), result.heavy_keys.end()) ||
     std::adjacent_find(result.heavy_keys.begin(), result.heavy_keys.end()) !=
       result.heavy_keys.end())
  {
    throw joinError("MSD range-pack plan is not canonical");
  }
  validateJoinRadixPacks(result);
  finalizeJoinRadixPlan(result);
  return result;
}

JoinRadixPlan
loadOrCreateJoinRadixPlan(const JoinFileReader& left,
  const JoinFileReader& right, const JoinRun& left_run, const JoinRun& right_run,
  logical_file_id_t logical,
  size_type target_bytes, size_type memory_budget, BuildWorkspace* workspace,
  const std::string& checkpoint_task, size_type reader_codec_bytes,
  ExternalPathJoinStats* stats)
{
  JoinRadixPlan result;
  size_type serialized_limit = joinPlanSerializedLimit(memory_budget,
    reader_codec_bytes);
  bool restored = false;
  std::string task;
  if(workspace != nullptr && !checkpoint_task.empty())
  {
    task = joinRadixPlanTaskName(checkpoint_task, logical, target_bytes,
      serialized_limit, left.checksum(), right.checksum());
    if(workspace->task_completed(task, "join-plan"))
    {
      result = decodeJoinRadixPlan(workspace->read_artifact_payload(
        joinRadixPlanArtifact(task), logical, physical_shard_id_t(0),
        serialized_limit), logical, target_bytes, left.size(), right.size(),
        left.checksum(), right.checksum(), serialized_limit, left_run, right_run);
      restored = true;
    }
  }
  if(!restored)
  {
    result = sampleJoinRadixPlan(left, right, target_bytes, serialized_limit);
    bindJoinRadixPlanSidecars(result, left_run, right_run);
    if(workspace != nullptr && !task.empty())
    {
      std::vector<std::uint8_t> payload = encodeJoinRadixPlan(result, logical,
        target_bytes, left.size(), right.size(), left.checksum(), right.checksum());
      ArtifactIdentity identity = joinRadixPlanArtifact(task);
      BuildWorkspace::ArtifactWriter writer = workspace->open_artifact(identity,
        logical, physical_shard_id_t(0), "msd-key-prefix", "all");
      writer.write(payload.data(), payload.size());
      std::vector<BuildWorkspace::ArtifactRef> artifacts;
      artifacts.push_back(writer.finish(result.packs.size() +
        result.heavy_keys.size()));
      workspace->commit_task(task, "join-plan", artifacts);
    }
  }
  if(stats != nullptr)
  {
    stats->sampled_plan_records += result.sampled_records;
    stats->radix_plan_bins += result.packs.size();
    stats->radix_plan_splits += result.splits;
    stats->radix_plan_max_bits = std::max(stats->radix_plan_max_bits,
      result.max_bits);
    stats->radix_plan_capped += (result.capped ? 1 : 0);
    size_type planner_bound = saturatingPlanAdd(
      saturatingPlanAdd(JOIN_PLAN_RUNTIME_OVERHEAD, reader_codec_bytes),
      saturatingPlanMultiply(4, serialized_limit));
    stats->max_bytes_resident = std::max(stats->max_bytes_resident,
      std::min(memory_budget, planner_bound));
    if(restored)
    {
      stats->restored_radix_plans++;
      stats->restored_sidecar_metadata += 2;
    }
  }
  return result;
}

size_type
checkedJoinAdd(size_type left, size_type right, const char* description)
{
  if(left > std::numeric_limits<size_type>::max() - right)
  {
    throw joinError(std::string("partition estimate overflow: ") + description);
  }
  return left + right;
}

size_type
checkedJoinMultiply(size_type left, size_type right, const char* description)
{
  if(left != 0 && right > std::numeric_limits<size_type>::max() / left)
  {
    throw joinError(std::string("partition estimate overflow: ") + description);
  }
  return left * right;
}

JoinGroupSummary
readJoinGroup(const JoinGroupSidecarReader& reader, size_type group,
  ExternalPathJoinStats* stats)
{
  JoinGroupSummary result;
  reader.read(group, result);
  if(result.count == 0)
  {
    throw joinError("join group sidecar contains an empty group");
  }
  if(stats != nullptr) { stats->sidecar_plan_groups++; }
  return result;
}

void
joinGroupOutput(const JoinGroupSummary* left, const JoinGroupSummary* right,
  bool emit_bypass, size_type& paths, size_type& ranks, size_type& bytes)
{
  paths = 0; ranks = 0;
  if(right != nullptr && emit_bypass)
  {
    paths = right->bypass_paths;
    ranks = right->bypass_ranks;
  }
  if(left != nullptr && right != nullptr && left->key == right->key)
  {
    size_type generated = checkedJoinMultiply(left->count, right->count,
      "generated paths");
    paths = checkedJoinAdd(paths, generated, "partition paths");
    size_type left_orders = checkedJoinMultiply(right->count, left->order_sum,
      "left order contribution");
    size_type right_orders = checkedJoinMultiply(left->count, right->order_sum,
      "right order contribution");
    ranks = checkedJoinAdd(ranks, generated, "generated terminal ranks");
    ranks = checkedJoinAdd(ranks, left_orders, "generated left ranks");
    ranks = checkedJoinAdd(ranks, right_orders, "generated right ranks");
  }
  bytes = checkedJoinAdd(checkedJoinMultiply(paths, sizeof(PathNode),
    "partition path bytes"), checkedJoinMultiply(ranks,
    sizeof(PathNode::rank_type), "partition rank bytes"),
    "partition payload bytes");
}

JoinGroupSummary
readJoinSubgroup(const JoinGroupSidecarReader& reader,
  const JoinGroupSummary& parent, size_type begin, size_type end,
  ExternalPathJoinStats* stats)
{
  if(stats != nullptr)
  {
    stats->sidecar_plan_detail_records = checkedJoinAdd(
      stats->sidecar_plan_detail_records, end - begin,
      "sidecar detail planning records");
  }
  return reader.summarize(parent, begin, end);
}

JoinGroupSummary
remainingJoinSubgroup(const JoinGroupSummary& parent,
  const JoinGroupSummary& prefix)
{
  if(prefix.begin != parent.begin || prefix.end > parent.end ||
     prefix.count > parent.count || prefix.order_sum > parent.order_sum ||
     prefix.bypass_paths > parent.bypass_paths ||
     prefix.bypass_ranks > parent.bypass_ranks)
  {
    throw joinError("invalid join-key subgroup summary");
  }
  JoinGroupSummary result;
  result.key = parent.key; result.begin = prefix.end; result.end = parent.end;
  result.count = parent.count - prefix.count;
  result.order_sum = parent.order_sum - prefix.order_sum;
  result.bypass_paths = parent.bypass_paths - prefix.bypass_paths;
  result.bypass_ranks = parent.bypass_ranks - prefix.bypass_ranks;
  return result;
}

void
planJoinGroupRecursive(const JoinGroupSidecarReader& left_reader,
  const JoinGroupSidecarReader& right_reader, const JoinGroupSummary& left,
  const JoinGroupSummary& right, size_type target_bytes, bool emit_bypass,
  std::vector<JoinPartition>& result, ExternalPathJoinStats* stats)
{
  size_type paths = 0, ranks = 0, bytes = 0;
  const JoinGroupSummary* left_ptr = (left.count == 0 ? nullptr : &left);
  joinGroupOutput(left_ptr, &right, emit_bypass, paths, ranks, bytes);
  size_type inputs = checkedJoinAdd(left.count, right.count,
    "partition input records");
  size_type input_bytes = checkedJoinMultiply(inputs, JOIN_RECORD_BYTES,
    "partition input bytes");
  const size_type rank_limit = (static_cast<size_type>(1) << 40);

  bool fits = (bytes <= target_bytes && input_bytes <= target_bytes &&
    ranks <= rank_limit);
  if(fits || (left.count <= 1 && right.count <= 1))
  {
    if(ranks > rank_limit)
    {
      throw joinError("one generated path exceeds the 40-bit rank pointer range");
    }
    if(paths == 0) { return; }
    JoinPartition partition;
    partition.left_begin = left.begin; partition.left_end = left.end;
    partition.right_begin = right.begin; partition.right_end = right.end;
    partition.expected_paths = paths; partition.expected_ranks = ranks;
    partition.expected_bytes = bytes; partition.emit_bypass = emit_bypass;
    partition.bypass_paths = (emit_bypass ? right.bypass_paths : 0);
    partition.generated_paths = paths - partition.bypass_paths;
    result.push_back(partition);
    return;
  }

  // Split by record ranges instead of materializing a high-multiplicity key.
  // Splitting the right side gives each child a disjoint bypass range. When the
  // left side is split, only its first child emits the shared right bypass.
  bool split_left = (left.count > 1 &&
    (right.count <= 1 || left.count >= right.count));
  if(stats != nullptr)
  {
    stats->recursive_splits++;
    if(split_left) { stats->left_range_splits++; }
    else { stats->right_range_splits++; }
  }
  if(split_left)
  {
    size_type middle = left.begin + left.count / 2;
    JoinGroupSummary first = readJoinSubgroup(left_reader, left,
      left.begin, middle, stats);
    JoinGroupSummary second = remainingJoinSubgroup(left, first);
    planJoinGroupRecursive(left_reader, right_reader, first, right,
      target_bytes, emit_bypass, result, stats);
    planJoinGroupRecursive(left_reader, right_reader, second, right,
      target_bytes, false, result, stats);
  }
  else
  {
    if(right.count <= 1)
    {
      throw joinError("cannot split an oversized join-key partition");
    }
    size_type middle = right.begin + right.count / 2;
    JoinGroupSummary first = readJoinSubgroup(right_reader, right,
      right.begin, middle, stats);
    JoinGroupSummary second = remainingJoinSubgroup(right, first);
    planJoinGroupRecursive(left_reader, right_reader, left, first,
      target_bytes, emit_bypass, result, stats);
    planJoinGroupRecursive(left_reader, right_reader, left, second,
      target_bytes, emit_bypass, result, stats);
  }
}

std::vector<JoinPartition>
planJoinPartitions(const JoinRun& left_run, const JoinRun& right_run,
  logical_file_id_t logical, size_type target_bytes, size_type memory_budget,
  bool verify_payloads,
  ExternalPathJoinStats* stats, BuildWorkspace* workspace = nullptr,
  const std::string& checkpoint_task = std::string())
{
  size_type reader_codec_bytes = joinRunReaderPairMemory(left_run, right_run);
  size_type minimum_planner = checkedJoinAdd(JOIN_PLAN_RUNTIME_OVERHEAD,
    reader_codec_bytes, "compressed join planner minimum");
  minimum_planner = checkedJoinAdd(minimum_planner,
    4 * (JOIN_PLAN_FIXED_BYTES + JOIN_PLAN_PACK_BYTES),
    "minimum join plan representations");
  if(memory_budget < minimum_planner)
  {
    throw joinError("join planner budget cannot admit both compressed run readers");
  }
  JoinFileReader left(left_run, logical, LEFT_BY_TO);
  JoinFileReader right(right_run, logical, RIGHT_BY_FROM);
  JoinGroupSidecarReader left_groups(left_run, logical, LEFT_BY_TO);
  JoinGroupSidecarReader right_groups(right_run, logical, RIGHT_BY_FROM);
  if(verify_payloads)
  {
    left.validate(); right.validate(); left_groups.validate(); right_groups.validate();
  }
  target_bytes = std::max(static_cast<size_type>(1), target_bytes);
  JoinRadixPlan radix_plan = loadOrCreateJoinRadixPlan(left, right, left_run,
    right_run, logical, target_bytes, memory_budget, workspace,
    checkpoint_task, reader_codec_bytes, stats);
  size_type radix_boundary = 0;

  std::vector<JoinPartition> result;
  size_type left_offset = 0, right_offset = 0;
  size_type left_group_index = 0, right_group_index = 0;
  JoinPartition current;
  current.left_begin = 0; current.right_begin = 0;
  size_type current_input_bytes = 0;
  auto flush_current = [&]()
  {
    current.left_end = left_offset; current.right_end = right_offset;
    if(current.expected_paths > 0) { result.push_back(current); }
    current = JoinPartition();
    current.left_begin = left_offset; current.right_begin = right_offset;
    current_input_bytes = 0;
  };
  JoinGroupSummary left_group, right_group;
  bool have_left = false, have_right = false;
  while(have_left || have_right || left_group_index < left_groups.size() ||
    right_group_index < right_groups.size())
  {
    if(!have_left && left_group_index < left_groups.size())
    {
      left_group = readJoinGroup(left_groups, left_group_index++, stats);
      if(left_group.begin != left_offset)
      {
        throw joinError("left join group sidecar is not contiguous");
      }
      have_left = true;
    }
    if(!have_right && right_group_index < right_groups.size())
    {
      right_group = readJoinGroup(right_groups, right_group_index++, stats);
      if(right_group.begin != right_offset)
      {
        throw joinError("right join group sidecar is not contiguous");
      }
      have_right = true;
    }
    const JoinGroupSummary* selected_left = nullptr;
    const JoinGroupSummary* selected_right = nullptr;
    if(have_left && (!have_right || left_group.key < right_group.key))
    {
      selected_left = &left_group;
    }
    else if(have_right && (!have_left || right_group.key < left_group.key))
    {
      selected_right = &right_group;
    }
    else
    {
      selected_left = &left_group; selected_right = &right_group;
    }

    node_type selected_key = (selected_left == nullptr ?
      selected_right->key : selected_left->key);
    bool crossed_radix_boundary = false;
    while(radix_boundary < radix_plan.boundaries.size() &&
          radix_plan.boundaries[radix_boundary] <= selected_key)
    {
      crossed_radix_boundary = true; radix_boundary++;
    }
    if(crossed_radix_boundary &&
       (current.expected_paths > 0 || current_input_bytes > 0))
    {
      flush_current();
      if(stats != nullptr) { stats->radix_boundary_flushes++; }
    }

    size_type group_paths = 0, group_ranks = 0, group_bytes = 0;
    joinGroupOutput(selected_left, selected_right, true,
      group_paths, group_ranks, group_bytes);

    // A skewed key is a legal partition boundary but not an excuse to require
    // its Cartesian product in RAM. Tile the larger side while replaying the
    // smaller side from its bounded join-run reader. This also keeps each
    // physical rank sidecar below PathNode's 40-bit pointer limit.
    size_type group_inputs = checkedJoinAdd(
      (selected_left == nullptr ? 0 : selected_left->count),
      (selected_right == nullptr ? 0 : selected_right->count),
      "group input records");
    size_type group_input_bytes = checkedJoinMultiply(group_inputs,
      JOIN_RECORD_BYTES, "group input bytes");
    bool oversized_group = (group_bytes > target_bytes ||
      group_input_bytes > target_bytes ||
      group_ranks > (static_cast<size_type>(1) << 40));
    // The bounded histogram lets a repeated semantic key begin range tiling
    // before it monopolizes an otherwise balanced worker task. Exact counts
    // above still decide the output contract and retain logical-file identity.
    bool sampled_heavy_group = (selected_left != nullptr && selected_right != nullptr &&
      std::binary_search(radix_plan.heavy_keys.begin(),
        radix_plan.heavy_keys.end(), selected_left->key));
    if((oversized_group || sampled_heavy_group) && selected_right != nullptr &&
       (selected_left != nullptr || selected_right->bypass_paths > 0))
    {
      flush_current();
      JoinGroupSummary empty_left;
      empty_left.key = selected_right->key;
      empty_left.begin = left_offset; empty_left.end = left_offset;
      planJoinGroupRecursive(left_groups, right_groups,
        (selected_left == nullptr ? empty_left : *selected_left),
        *selected_right, target_bytes, true, result, stats);
      if(selected_left != nullptr) { left_offset = selected_left->end; }
      right_offset = selected_right->end;
      have_left = false; have_right = false;
      current.left_begin = left_offset; current.right_begin = right_offset;
      continue;
    }

    // Keep ordinary partitions close to the requested checkpoint size. Flush
    // before adding the next complete key group, because key boundaries are the
    // cheapest recovery boundaries and never duplicate Cartesian products.
    const size_type rank_limit = (static_cast<size_type>(1) << 40);
    if((current.expected_paths > 0 || current_input_bytes > 0) &&
       (current.expected_bytes > target_bytes - std::min(target_bytes, group_bytes) ||
        current_input_bytes > target_bytes - std::min(target_bytes, group_input_bytes) ||
        current.expected_ranks > rank_limit - std::min(rank_limit, group_ranks)))
    {
      flush_current();
    }

    current.expected_paths = checkedJoinAdd(current.expected_paths,
      group_paths, "task paths");
    size_type group_bypass = (selected_right == nullptr ? 0 :
      selected_right->bypass_paths);
    current.bypass_paths = checkedJoinAdd(current.bypass_paths,
      group_bypass, "task bypass paths");
    current.generated_paths = checkedJoinAdd(current.generated_paths,
      group_paths - group_bypass, "task generated paths");
    current.expected_ranks = checkedJoinAdd(current.expected_ranks,
      group_ranks, "task ranks");
    current.expected_bytes = checkedJoinAdd(current.expected_bytes,
      group_bytes, "task bytes");
    if(selected_left != nullptr)
    {
      left_offset = selected_left->end;
      have_left = false;
    }
    if(selected_right != nullptr)
    {
      right_offset = selected_right->end;
      have_right = false;
    }
    current_input_bytes = checkedJoinAdd(current_input_bytes,
      group_input_bytes,
      "task input bytes");
    current.left_end = left_offset; current.right_end = right_offset;

    if(current.expected_bytes >= target_bytes || current_input_bytes >= target_bytes)
    {
      flush_current();
    }
  }
  flush_current();
  if(left_offset != left.size() || right_offset != right.size())
  {
    throw joinError("join group sidecar did not cover its join run");
  }
  if(result.empty())
  {
    current.left_begin = 0; current.left_end = left.size();
    current.right_begin = 0; current.right_end = right.size();
    result.push_back(current);
  }
  return result;
}

constexpr std::uint64_t WORKER_TASK_MAGIC = 0x314b535441534347ULL;   // "GCSATSK1"
constexpr std::uint64_t WORKER_RESULT_MAGIC = 0x3153455241534347ULL; // "GCSARES1"
constexpr std::uint32_t WORKER_FORMAT_VERSION = 5;
constexpr size_type WORKER_CONTROL_LIMIT = MEGABYTE;

struct ExternalJoinWorkerTask
{
  logical_file_id_t logical;
  JoinRun left, right;
  JoinPartition partition;
  std::string output_path, output_rank, temp_directory;
  size_type sort_budget, join_block_budget, fan_in, threads;
  TempFileCodecParameters codec;
  bool verify_payloads;
};

struct ExternalJoinWorkerResult
{
  size_type paths, ranks, bytes, stored_bytes;
  size_type generated, bypassed, label_runs, label_merge_passes, label_parallel_sorts;
  size_type grouped_records, group_headers, context_bytes_saved;
  size_type max_records, max_bytes, blocks, bytes_read, bytes_written;

  ExternalJoinWorkerResult() : paths(0), ranks(0), bytes(0), stored_bytes(0), generated(0),
    bypassed(0), label_runs(0), label_merge_passes(0), label_parallel_sorts(0),
    grouped_records(0), group_headers(0), context_bytes_saved(0),
    max_records(0), max_bytes(0), blocks(0), bytes_read(0), bytes_written(0) { }
};

template<class Value>
void
appendWorkerValue(std::vector<std::uint8_t>& target, Value value)
{
  size_type old_size = target.size();
  target.resize(old_size + sizeof(Value));
  std::uint8_t* output = target.data() + old_size;
  encodeLittle<Value>(output, value);
}

void
appendWorkerString(std::vector<std::uint8_t>& target, const std::string& value)
{
  appendWorkerValue<std::uint64_t>(target, value.size());
  target.insert(target.end(), value.begin(), value.end());
}

template<class Value>
Value
readWorkerValue(const std::vector<std::uint8_t>& source, size_type& offset)
{
  if(offset > source.size() || sizeof(Value) > source.size() - offset)
  {
    throw joinError("truncated worker control file");
  }
  const std::uint8_t* input = source.data() + offset;
  Value result = decodeLittle<Value>(input);
  offset += sizeof(Value); return result;
}

std::string
readWorkerString(const std::vector<std::uint8_t>& source, size_type& offset)
{
  size_type length = readWorkerValue<std::uint64_t>(source, offset);
  if(offset > source.size() || length > source.size() - offset)
  {
    throw joinError("truncated worker control string");
  }
  std::string result(reinterpret_cast<const char*>(source.data() + offset), length);
  offset += length; return result;
}

std::vector<std::uint8_t>
readWorkerControl(const std::string& path)
{
  std::ifstream input(path.c_str(), std::ios_base::binary);
  if(!input) { throw joinError("cannot open worker control file", path); }
  input.seekg(0, std::ios_base::end);
  std::streamoff end = input.tellg();
  if(end < 0 || static_cast<size_type>(end) > WORKER_CONTROL_LIMIT)
  {
    throw joinError("invalid worker control file size", path);
  }
  input.seekg(0, std::ios_base::beg);
  std::vector<std::uint8_t> result(static_cast<size_type>(end));
  if(!result.empty() && !DiskIO::read(input, result.data(), result.size()))
  {
    throw joinError("short worker control file", path);
  }
  return result;
}

void
writeWorkerControl(const std::string& path, const std::vector<std::uint8_t>& data)
{
  std::string partial = path + "." +
    std::to_string(static_cast<unsigned long long>(getpid())) + ".partial";
  int descriptor = ::open(partial.c_str(), O_CREAT | O_EXCL | O_WRONLY, 0644);
  if(descriptor < 0) { throw joinError("cannot create worker control file", partial); }
  try
  {
    writeAll(descriptor, data.data(), data.size(), partial);
    if(::fdatasync(descriptor) != 0) { throw joinError("cannot sync worker control file", partial); }
    if(::close(descriptor) != 0) { descriptor = -1; throw joinError("cannot close worker control file", partial); }
    descriptor = -1;
    if(::rename(partial.c_str(), path.c_str()) != 0)
    {
      throw joinError("cannot publish worker control file", path);
    }
  }
  catch(...)
  {
    if(descriptor >= 0) { ::close(descriptor); }
    ::unlink(partial.c_str()); throw;
  }
}

std::vector<std::uint8_t>
encodeWorkerTask(const ExternalJoinWorkerTask& task)
{
  std::vector<std::uint8_t> data;
  appendWorkerValue<std::uint64_t>(data, WORKER_TASK_MAGIC);
  appendWorkerValue<std::uint32_t>(data, WORKER_FORMAT_VERSION);
  appendWorkerValue<std::uint32_t>(data, task.logical.value);
  appendWorkerValue<std::uint64_t>(data, task.left.records);
  appendWorkerValue<std::uint64_t>(data, task.right.records);
  appendWorkerValue<std::uint64_t>(data, task.left.checksum);
  appendWorkerValue<std::uint64_t>(data, task.right.checksum);
  appendWorkerValue<std::uint64_t>(data, task.partition.left_begin);
  appendWorkerValue<std::uint64_t>(data, task.partition.left_end);
  appendWorkerValue<std::uint64_t>(data, task.partition.right_begin);
  appendWorkerValue<std::uint64_t>(data, task.partition.right_end);
  appendWorkerValue<std::uint64_t>(data, task.partition.expected_paths);
  appendWorkerValue<std::uint64_t>(data, task.partition.expected_ranks);
  appendWorkerValue<std::uint64_t>(data, task.partition.expected_bytes);
  appendWorkerValue<std::uint64_t>(data, task.partition.generated_paths);
  appendWorkerValue<std::uint64_t>(data, task.partition.bypass_paths);
  appendWorkerValue<std::uint64_t>(data, task.sort_budget);
  appendWorkerValue<std::uint64_t>(data, task.join_block_budget);
  appendWorkerValue<std::uint64_t>(data, task.fan_in);
  appendWorkerValue<std::uint64_t>(data, task.threads);
  appendWorkerValue<std::uint8_t>(data,
    static_cast<std::uint8_t>(task.codec.compression));
  appendWorkerValue<std::uint64_t>(data, task.codec.block_size);
  appendWorkerValue<std::uint64_t>(data, task.codec.workers);
  appendWorkerValue<std::uint32_t>(data,
    static_cast<std::uint32_t>(static_cast<std::int32_t>(task.codec.level)));
  appendWorkerValue<std::uint8_t>(data, task.partition.emit_bypass ? 1 : 0);
  appendWorkerValue<std::uint8_t>(data, task.verify_payloads ? 1 : 0);
  appendWorkerString(data, task.left.name);
  appendWorkerString(data, task.right.name);
  appendWorkerString(data, task.output_path);
  appendWorkerString(data, task.output_rank);
  appendWorkerString(data, task.temp_directory);
  return data;
}

ExternalJoinWorkerTask
decodeWorkerTask(const std::vector<std::uint8_t>& data)
{
  size_type offset = 0;
  if(readWorkerValue<std::uint64_t>(data, offset) != WORKER_TASK_MAGIC ||
     readWorkerValue<std::uint32_t>(data, offset) != WORKER_FORMAT_VERSION)
  {
    throw joinError("incompatible worker task format");
  }
  ExternalJoinWorkerTask task;
  task.logical = logical_file_id_t(readWorkerValue<std::uint32_t>(data, offset));
  task.left.records = readWorkerValue<std::uint64_t>(data, offset);
  task.right.records = readWorkerValue<std::uint64_t>(data, offset);
  task.left.checksum = readWorkerValue<std::uint64_t>(data, offset);
  task.right.checksum = readWorkerValue<std::uint64_t>(data, offset);
  task.partition.left_begin = readWorkerValue<std::uint64_t>(data, offset);
  task.partition.left_end = readWorkerValue<std::uint64_t>(data, offset);
  task.partition.right_begin = readWorkerValue<std::uint64_t>(data, offset);
  task.partition.right_end = readWorkerValue<std::uint64_t>(data, offset);
  task.partition.expected_paths = readWorkerValue<std::uint64_t>(data, offset);
  task.partition.expected_ranks = readWorkerValue<std::uint64_t>(data, offset);
  task.partition.expected_bytes = readWorkerValue<std::uint64_t>(data, offset);
  task.partition.generated_paths = readWorkerValue<std::uint64_t>(data, offset);
  task.partition.bypass_paths = readWorkerValue<std::uint64_t>(data, offset);
  task.sort_budget = readWorkerValue<std::uint64_t>(data, offset);
  task.join_block_budget = readWorkerValue<std::uint64_t>(data, offset);
  task.fan_in = readWorkerValue<std::uint64_t>(data, offset);
  task.threads = readWorkerValue<std::uint64_t>(data, offset);
  const std::uint8_t compression = readWorkerValue<std::uint8_t>(data, offset);
  if(compression > static_cast<std::uint8_t>(TempCompression::ZSTD))
  {
    throw joinError("invalid temporary compression mode in worker task");
  }
  task.codec.compression = static_cast<TempCompression>(compression);
  task.codec.block_size = readWorkerValue<std::uint64_t>(data, offset);
  task.codec.workers = readWorkerValue<std::uint64_t>(data, offset);
  task.codec.level = static_cast<std::int32_t>(
    readWorkerValue<std::uint32_t>(data, offset));
  if(task.codec.block_size == 0 || task.codec.workers == 0 ||
     task.codec.level < -5 || task.codec.level > 22)
  {
    throw joinError("invalid temporary compression settings in worker task");
  }
  task.partition.emit_bypass = (readWorkerValue<std::uint8_t>(data, offset) != 0);
  task.verify_payloads = (readWorkerValue<std::uint8_t>(data, offset) != 0);
  task.left.name = readWorkerString(data, offset);
  task.right.name = readWorkerString(data, offset);
  task.output_path = readWorkerString(data, offset);
  task.output_rank = readWorkerString(data, offset);
  task.temp_directory = readWorkerString(data, offset);
  if(offset != data.size()) { throw joinError("worker task has trailing data"); }
  return task;
}

std::vector<std::uint8_t>
encodeWorkerResult(const ExternalJoinWorkerResult& result)
{
  std::vector<std::uint8_t> data;
  appendWorkerValue<std::uint64_t>(data, WORKER_RESULT_MAGIC);
  appendWorkerValue<std::uint32_t>(data, WORKER_FORMAT_VERSION);
  appendWorkerValue<std::uint64_t>(data, result.paths);
  appendWorkerValue<std::uint64_t>(data, result.ranks);
  appendWorkerValue<std::uint64_t>(data, result.bytes);
  appendWorkerValue<std::uint64_t>(data, result.stored_bytes);
  appendWorkerValue<std::uint64_t>(data, result.generated);
  appendWorkerValue<std::uint64_t>(data, result.bypassed);
  appendWorkerValue<std::uint64_t>(data, result.label_runs);
  appendWorkerValue<std::uint64_t>(data, result.label_merge_passes);
  appendWorkerValue<std::uint64_t>(data, result.label_parallel_sorts);
  appendWorkerValue<std::uint64_t>(data, result.grouped_records);
  appendWorkerValue<std::uint64_t>(data, result.group_headers);
  appendWorkerValue<std::uint64_t>(data, result.context_bytes_saved);
  appendWorkerValue<std::uint64_t>(data, result.max_records);
  appendWorkerValue<std::uint64_t>(data, result.max_bytes);
  appendWorkerValue<std::uint64_t>(data, result.blocks);
  appendWorkerValue<std::uint64_t>(data, result.bytes_read);
  appendWorkerValue<std::uint64_t>(data, result.bytes_written);
  return data;
}

ExternalJoinWorkerResult
decodeWorkerResult(const std::vector<std::uint8_t>& data)
{
  size_type offset = 0;
  if(readWorkerValue<std::uint64_t>(data, offset) != WORKER_RESULT_MAGIC ||
     readWorkerValue<std::uint32_t>(data, offset) != WORKER_FORMAT_VERSION)
  {
    throw joinError("incompatible worker result format");
  }
  ExternalJoinWorkerResult result;
  result.paths = readWorkerValue<std::uint64_t>(data, offset);
  result.ranks = readWorkerValue<std::uint64_t>(data, offset);
  result.bytes = readWorkerValue<std::uint64_t>(data, offset);
  result.stored_bytes = readWorkerValue<std::uint64_t>(data, offset);
  result.generated = readWorkerValue<std::uint64_t>(data, offset);
  result.bypassed = readWorkerValue<std::uint64_t>(data, offset);
  result.label_runs = readWorkerValue<std::uint64_t>(data, offset);
  result.label_merge_passes = readWorkerValue<std::uint64_t>(data, offset);
  result.label_parallel_sorts = readWorkerValue<std::uint64_t>(data, offset);
  result.grouped_records = readWorkerValue<std::uint64_t>(data, offset);
  result.group_headers = readWorkerValue<std::uint64_t>(data, offset);
  result.context_bytes_saved = readWorkerValue<std::uint64_t>(data, offset);
  result.max_records = readWorkerValue<std::uint64_t>(data, offset);
  result.max_bytes = readWorkerValue<std::uint64_t>(data, offset);
  result.blocks = readWorkerValue<std::uint64_t>(data, offset);
  result.bytes_read = readWorkerValue<std::uint64_t>(data, offset);
  result.bytes_written = readWorkerValue<std::uint64_t>(data, offset);
  if(offset != data.size()) { throw joinError("worker result has trailing data"); }
  return result;
}

} // namespace

int
externalPathJoinWorker(const std::string& task_file)
{
  try
  {
    ExternalJoinWorkerTask task = decodeWorkerTask(readWorkerControl(task_file));
    if(task.partition.left_end > task.left.records ||
       task.partition.right_end > task.right.records ||
       task.partition.left_begin > task.partition.left_end ||
       task.partition.right_begin > task.partition.right_end)
    {
      throw joinError("worker task range exceeds its immutable input run");
    }
    if(task.partition.expected_ranks > (static_cast<size_type>(1) << 40))
    {
      throw joinError("worker partition exceeds the 40-bit rank pointer range");
    }
    TempFile::setDirectory(task.temp_directory);
    omp_set_num_threads(std::max(static_cast<size_type>(1), task.threads));
    ::unlink(task.output_path.c_str()); ::unlink(task.output_rank.c_str());
    ::unlink((task.output_path + ".partial").c_str());
    ::unlink((task.output_rank + ".partial").c_str());
    ::unlink((task_file + ".result").c_str());

    size_type read_start = DiskIO::read_volume.load();
    size_type write_start = DiskIO::write_volume.load();
    PathGraph output(1, 0, 0);
    TempFile::remove(output.path_names[0]); TempFile::remove(output.rank_names[0]);
    output.path_names[0] = task.output_path;
    output.rank_names[0] = task.output_rank;
    output.logical_file_ids[0] = task.logical;
    size_type committed_bytes = 0;
    ExternalPathSortStats sort_stats;
    ExternalPathJoinStats join_stats;
    const size_type peak_output_bytes = externalPathGraphShardPeakBytes(
      task.partition.expected_paths, task.partition.expected_ranks,
      task.sort_budget, task.codec);
    ExternalPathSortSink sink(output, 0, task.sort_budget, task.fan_in,
      peak_output_bytes, committed_bytes, &sort_stats, task.codec);
    joinSortedRuns(task.left, task.right, task.logical, task.join_block_budget,
      sink, task.verify_payloads, &join_stats,
      task.partition.left_begin, task.partition.left_end,
      task.partition.right_begin, task.partition.right_end,
      task.partition.emit_bypass);
    sink.finish();
    if(output.path_counts[0] != task.partition.expected_paths ||
       output.rank_counts[0] != task.partition.expected_ranks ||
       sink.bytes() != task.partition.expected_bytes ||
       join_stats.generated_records != task.partition.generated_paths ||
       join_stats.sorted_bypass != task.partition.bypass_paths)
    {
      throw joinError("worker output does not match its deterministic partition plan");
    }

    ExternalJoinWorkerResult result;
    result.paths = output.path_counts[0]; result.ranks = output.rank_counts[0];
    result.bytes = sink.bytes(); result.stored_bytes = sink.storedBytes();
    result.generated = join_stats.generated_records;
    result.bypassed = join_stats.sorted_bypass; result.label_runs = sort_stats.runs;
    result.label_merge_passes = sort_stats.merge_passes;
    result.label_parallel_sorts = sort_stats.parallel_sorts;
    result.grouped_records = sort_stats.grouped_records;
    result.group_headers = sort_stats.group_headers;
    result.context_bytes_saved = sort_stats.context_bytes_saved;
    result.max_records = sort_stats.max_records_resident;
    result.max_bytes = checkedJoinAdd(sort_stats.max_bytes_resident,
      join_stats.max_bytes_resident, "worker combined resident bytes");
    result.blocks = join_stats.blocked_key_blocks;
    result.bytes_read = DiskIO::read_volume.load() - read_start;
    result.bytes_written = DiskIO::write_volume.load() - write_start;
    // The parent owns these completed files and will either admit them into the
    // next PathGraph or remove them if another worker fails.
    output.delete_files = false;
    writeWorkerControl(task_file + ".result", encodeWorkerResult(result));
    return EXIT_SUCCESS;
  }
  catch(const std::exception& error)
  {
    std::cerr << "externalPathJoinWorker(): " << error.what() << std::endl;
    return EXIT_FAILURE;
  }
}

namespace
{

size_type
appendOutputShard(PathGraph& graph, logical_file_id_t logical,
  physical_shard_id_t physical)
{
  size_type file = graph.files();
  graph.path_names.push_back(TempFile::getName(PathGraph::PREFIX));
  graph.rank_names.push_back(TempFile::getName(PathGraph::PREFIX));
  graph.path_counts.push_back(0); graph.rank_counts.push_back(0);
  graph.logical_file_ids.push_back(logical);
  graph.physical_shard_ids.push_back(physical);
  return file;
}

size_type
pathPairPeakForCodec(size_type paths, size_type ranks,
  size_type output_pairs, const TempFileCodecParameters& codec)
{
  const size_type path_bytes = checkedJoinMultiply(paths, sizeof(PathNode),
    "peak path output bytes");
  const size_type rank_bytes = checkedJoinMultiply(ranks,
    sizeof(PathNode::rank_type), "peak rank output bytes");
  if(!codec.enabled())
  {
    return checkedJoinAdd(path_bytes, rank_bytes, "peak raw path/rank bytes");
  }
  std::uint64_t path_peak = CompressedBlockWriter::maximumTemporaryBytes(
    path_bytes, sizeof(PathNode), codec.block_size);
  std::uint64_t rank_peak = CompressedBlockWriter::maximumTemporaryBytes(
    rank_bytes, (PathLabel::LABEL_LENGTH + 1) * sizeof(PathNode::rank_type),
    codec.block_size);
  if(path_peak > std::numeric_limits<size_type>::max() ||
     rank_peak > std::numeric_limits<size_type>::max() - path_peak)
  {
    throw joinError("framed path/rank output size overflow");
  }
  size_type result = static_cast<size_type>(path_peak + rank_peak);
  if(output_pairs > 1)
  {
    // Splitting a stream can add at most one partial block per new file. Each
    // framed pair also repeats two headers/footers. maximumTemporaryBytes()
    // charges 72 bytes per peak block and 88 fixed bytes per stream.
    const size_type split_overhead = checkedJoinMultiply(output_pairs - 1,
      2 * (72 + 88), "split framed path/rank overhead");
    result = checkedJoinAdd(result, split_overhead,
      "split framed path/rank peak bytes");
  }
  return result;
}

// Workers deliberately produce independent physical shards. Before exposing a
// generation, collapse those shards per logical input so downstream mergers do
// not allocate one pair of buffered readers (and FDs) for every partition.
void
compactLogicalJoinShards(PathGraph& source, size_type size_limit,
  const ConstructionParameters& parameters, size_type label_fan_in,
  MemoryBudget& memory, size_type& committed_bytes,
  ExternalPathJoinStats* stats)
{
  std::map<logical_file_id_t, std::vector<size_type>> groups;
  for(size_type file = 0; file < source.files(); file++)
  {
    groups[source.logicalFile(file)].push_back(file);
  }
  bool needs_compaction = false;
  for(auto& group : groups)
  {
    std::sort(group.second.begin(), group.second.end(), [&source](size_type left, size_type right)
    {
      if(source.physicalShard(left) != source.physicalShard(right))
      {
        return source.physicalShard(left) < source.physicalShard(right);
      }
      return left < right;
    });
    needs_compaction = needs_compaction || group.second.size() > 1;
  }
  if(!needs_compaction) { return; }

  // Each input is already in the complete path-label order used by
  // pathSortLess().  Re-sorting the partition output would therefore turn a
  // merge-only operation into a full extra external sort.  Keep two file
  // descriptors per reader, two for the output pair, and two descriptors of
  // reserve for the surrounding construction.
  const size_type available_memory = static_cast<size_type>(memory.available());
  const size_type base_reader_bytes = 2 * JOIN_IO_BUFFER_BYTES + sizeof(JoinRecord) +
    4 * sizeof(size_type);
  size_type maximum_reader_codec = 0;
  for(size_type file = 0; file < source.files(); file++)
  {
    if(CompressedBlockReader::isFramed(source.path_names[file]))
    {
      size_type codec_bytes = checkedJoinAdd(
        CompressedBlockReader::workingMemoryEstimate(
          CompressedBlockReader::declaredBlockSize(source.path_names[file])),
        CompressedBlockReader::workingMemoryEstimate(
          CompressedBlockReader::declaredBlockSize(source.rank_names[file])),
        "logical merge compressed reader bytes");
      maximum_reader_codec = std::max(maximum_reader_codec, codec_bytes);
    }
  }
  const size_type reader_bytes = checkedJoinAdd(base_reader_bytes,
    maximum_reader_codec, "logical merge reader reservation");

  TempFileCodecParameters output_codec = parameters.getTempFileCodecParameters();
  const size_type minimum_block = std::max(sizeof(PathNode),
    (PathLabel::LABEL_LENGTH + 1) * sizeof(PathNode::rank_type));
  size_type writer_bytes = 2 * JOIN_IO_BUFFER_BYTES;
  if(output_codec.enabled())
  {
    output_codec.block_size = std::max(minimum_block,
      std::min(output_codec.block_size, available_memory));
    auto writer_memory = [&]() -> size_type
    {
      return checkedJoinMultiply(2,
        CompressedBlockWriter::workingMemoryEstimate(output_codec.block_size,
          CompressedBlockWriter::ZSTD, output_codec.level,
          output_codec.workers),
        "logical merge compressed writer bytes");
    };
    const size_type two_readers = checkedJoinMultiply(2, reader_bytes,
      "minimum logical merge readers");
    while(output_codec.block_size > minimum_block &&
          (writer_memory() > available_memory ||
           two_readers > available_memory - writer_memory()))
    {
      output_codec.block_size = std::max(minimum_block,
        output_codec.block_size / 2);
    }
    writer_bytes = writer_memory();
    if(writer_bytes > available_memory ||
       two_readers > available_memory - writer_bytes)
    {
      if(output_codec.compression == TempCompression::AUTO)
      {
        output_codec.compression = TempCompression::NONE;
        writer_bytes = 2 * JOIN_IO_BUFFER_BYTES;
      }
      else
      {
        return; // Optional compaction cannot fit the explicit codec workspace.
      }
    }
  }
  size_type merge_fan_in = std::min(label_fan_in,
    (parameters.getMaxOpenFiles() - 4) / 2);
  if(merge_fan_in > 0)
  {
    merge_fan_in = std::min(merge_fan_in,
      (available_memory > writer_bytes ? (available_memory - writer_bytes) / reader_bytes : 0));
  }
  if(merge_fan_in < 2)
  {
    return; // Compaction is optional; the existing LABEL streams are valid.
  }
  bool needs_merge = false;
  for(const auto& group : groups)
  {
    size_type target_pairs = std::max(static_cast<size_type>(1),
      (parameters.getMaxOpenFiles() - 4) / (2 * groups.size()));
    if(group.second.size() > target_pairs) { needs_merge = true; break; }
  }
  if(!needs_merge) { return; }
  // Account for all reader state, heap entries, and both byte-bounded output
  // buffers. The reservation is held for the entire compacting pass.
  size_type merge_reservation = checkedJoinAdd(writer_bytes,
    checkedJoinMultiply(merge_fan_in, reader_bytes, "logical merge reader bytes"),
    "logical merge reservation");
  MemoryBudget::Reservation merge_memory = memory.reserve(merge_reservation,
    "logical-shard-merge");
  size_type staged_committed = committed_bytes;
  std::vector<std::string> new_outputs;
  std::map<std::string, size_type> new_output_charges;
  struct NewOutputCleanup
  {
    std::vector<std::string>& names; bool keep;
    NewOutputCleanup(std::vector<std::string>& output_names) : names(output_names), keep(false) { }
    ~NewOutputCleanup()
    {
      if(!this->keep) { for(std::string& name : this->names) { TempFile::remove(name); } }
    }
  } cleanup(new_outputs);

  struct HeapLess
  {
    const std::vector<std::unique_ptr<PathShardReader>>* readers;
    const std::vector<JoinRecord>* records;
    bool operator()(size_type left, size_type right) const
    {
      const JoinRecord& a = records->at(left), &b = records->at(right);
      size_type order = std::min(a.node.order(), b.node.order());
      for(size_type i = 0; i < order; i++)
      {
        if(a.labels[i] != b.labels[i]) { return a.labels[i] > b.labels[i]; }
      }
      if(a.node.order() != b.node.order()) { return a.node.order() > b.node.order(); }
      if(a.node.from != b.node.from) { return a.node.from > b.node.from; }
      if(a.node.to != b.node.to) { return a.node.to > b.node.to; }
      if(a.node.predecessors() != b.node.predecessors())
      {
        return a.node.predecessors() > b.node.predecessors();
      }
      if(a.node.lcp() != b.node.lcp()) { return a.node.lcp() > b.node.lcp(); }
      for(size_type i = 0; i < a.node.ranks(); i++)
      {
        if(a.labels[i] != b.labels[i]) { return a.labels[i] > b.labels[i]; }
      }
      return false; // Equal records are byte-equivalent after pointer rewriting.
    }
  };

  auto merge_batch = [&](const PathGraph& input, const std::vector<size_type>& shards,
    PathGraph& output, logical_file_id_t logical, size_type& next_physical) -> bool
  {
    size_type paths = 0, ranks = 0;
    for(size_type shard : shards)
    {
      paths = checkedJoinAdd(paths, input.path_counts[shard], "logical merge path count");
      ranks = checkedJoinAdd(ranks, input.rank_counts[shard], "logical merge rank count");
    }
    const size_type pointer_limit = (static_cast<size_type>(1) << 40);
    const size_type rank_capacity = pointer_limit - PathLabel::LABEL_LENGTH;
    const size_type output_pairs = (ranks == 0 ? 1 :
      1 + (ranks - 1) / rank_capacity);
    const size_type peak_output_bytes = pathPairPeakForCodec(paths, ranks,
      output_pairs, output_codec);
    if(peak_output_bytes > size_limit ||
       staged_committed > size_limit - peak_output_bytes)
    {
      return false; // Preserve the original bounded streams when duplication will not fit.
    }
    int path_fd = -1, rank_fd = -1;
    std::unique_ptr<CompressedBlockWriter> compressed_path, compressed_rank;
    std::string partial_path, partial_rank;
    size_type file = 0, written_paths = 0, written_ranks = 0;
    size_type emitted_paths = 0, emitted_ranks = 0;
    off_t path_released = 0, rank_released = 0;
    std::vector<std::uint8_t> path_buffer, rank_buffer;
    auto close_output = [&]()
    {
      if(compressed_path)
      {
        compressed_path->finish(); compressed_rank->finish();
        compressed_path.reset(); compressed_rank.reset();
      }
      if(path_fd >= 0) { trimWrittenCache(path_fd, path_released, true, partial_path); }
      if(rank_fd >= 0) { trimWrittenCache(rank_fd, rank_released, true, partial_rank); }
      int path_error = (path_fd >= 0 ? ::close(path_fd) : 0); path_fd = -1;
      int rank_error = (rank_fd >= 0 ? ::close(rank_fd) : 0); rank_fd = -1;
      if(path_error != 0 || rank_error != 0) { throw joinError("cannot close logical merge output"); }
      // Register both final names before either rename. If the second install
      // fails after the first succeeds, generation rollback still owns the
      // partially installed pair.
      new_outputs.push_back(output.path_names[file]); new_outputs.push_back(output.rank_names[file]);
      const size_type stored_path_bytes = storedJoinBytes(partial_path);
      const size_type stored_rank_bytes = storedJoinBytes(partial_rank);
      new_output_charges[output.path_names[file]] = stored_path_bytes;
      new_output_charges[output.rank_names[file]] = stored_rank_bytes;
      if(::rename(partial_path.c_str(), output.path_names[file].c_str()) != 0 ||
         ::rename(partial_rank.c_str(), output.rank_names[file].c_str()) != 0)
      {
        throw joinError("cannot install logical merge output");
      }
      output.path_counts[file] = written_paths; output.rank_counts[file] = written_ranks;
      output.path_count += written_paths; output.rank_count += written_ranks;
      const size_type output_bytes = checkedJoinAdd(stored_path_bytes,
        stored_rank_bytes, "stored logical-merge output bytes");
      staged_committed = checkedJoinAdd(staged_committed, output_bytes,
        "committed stored output bytes");
    };
    auto open_output = [&]()
    {
      file = appendOutputShard(output, logical, physical_shard_id_t(next_physical++));
      partial_path = output.path_names[file] + ".partial";
      partial_rank = output.rank_names[file] + ".partial";
      std::remove(partial_path.c_str()); std::remove(partial_rank.c_str());
      if(output_codec.enabled())
      {
        compressed_path.reset(new CompressedBlockWriter(partial_path,
          output_codec.block_size, CompressedBlockWriter::ZSTD,
          output_codec.level, output_codec.workers));
        compressed_rank.reset(new CompressedBlockWriter(partial_rank,
          output_codec.block_size, CompressedBlockWriter::ZSTD,
          output_codec.level, output_codec.workers));
        written_paths = 0; written_ranks = 0;
        return;
      }
      path_fd = ::open(partial_path.c_str(), O_CREAT | O_EXCL | O_WRONLY, 0600);
      rank_fd = ::open(partial_rank.c_str(), O_CREAT | O_EXCL | O_WRONLY, 0600);
      if(path_fd < 0 || rank_fd < 0)
      {
        throw joinError("cannot create logical merge output: " + std::string(std::strerror(errno)),
          (path_fd < 0 ? partial_path : partial_rank));
      }
      adviseSequential(path_fd); adviseSequential(rank_fd);
      written_paths = 0; written_ranks = 0; path_released = 0; rank_released = 0;
      path_buffer.clear(); rank_buffer.clear();
      path_buffer.reserve(JOIN_IO_BUFFER_BYTES); rank_buffer.reserve(JOIN_IO_BUFFER_BYTES);
    };
    auto flush = [&](std::vector<std::uint8_t>& buffer, int descriptor, off_t& released,
      const std::string& name)
    {
      if(!buffer.empty()) { writeAll(descriptor, buffer.data(), buffer.size(), name); buffer.clear(); }
      trimWrittenCache(descriptor, released, false, name);
    };
    try
    {
      std::vector<std::unique_ptr<PathShardReader>> readers;
      std::vector<JoinRecord> records(shards.size()); readers.reserve(shards.size());
      std::priority_queue<size_type, std::vector<size_type>, HeapLess> queue((HeapLess{ &readers, &records }));
      for(size_type i = 0; i < shards.size(); i++)
      {
        readers.emplace_back(new PathShardReader(input, shards[i]));
        if(readers.back()->read(records[i])) { queue.push(i); }
      }
      open_output();
      while(!queue.empty())
      {
        size_type best = queue.top(); queue.pop();
        PathNode node = records[best].node;
        if(node.ranks() > pointer_limit)
        {
          throw joinError("one path record exceeds the 40-bit rank pointer range");
        }
        if(written_ranks > pointer_limit - node.ranks())
        {
          if(!compressed_path)
          {
            flush(path_buffer, path_fd, path_released, partial_path);
            flush(rank_buffer, rank_fd, rank_released, partial_rank);
          }
          close_output(); open_output();
        }
        node.setPointer(written_ranks);
        const std::uint8_t* path_bytes = reinterpret_cast<const std::uint8_t*>(&node);
        const std::uint8_t* rank_bytes = reinterpret_cast<const std::uint8_t*>(records[best].labels);
        if(compressed_path)
        {
          compressed_path->writeRecord(path_bytes, sizeof(node));
          if(node.ranks() > 0)
          {
            compressed_rank->writeRecord(rank_bytes,
              node.ranks() * sizeof(PathNode::rank_type));
          }
          written_paths++; written_ranks += node.ranks();
          emitted_paths++; emitted_ranks += node.ranks();
          if(readers[best]->read(records[best])) { queue.push(best); }
          continue;
        }
        if(path_buffer.size() + sizeof(node) > JOIN_IO_BUFFER_BYTES)
        {
          flush(path_buffer, path_fd, path_released, partial_path);
        }
        if(rank_buffer.size() + node.ranks() * sizeof(PathNode::rank_type) > JOIN_IO_BUFFER_BYTES)
        {
          flush(rank_buffer, rank_fd, rank_released, partial_rank);
        }
        path_buffer.insert(path_buffer.end(), path_bytes, path_bytes + sizeof(node));
        rank_buffer.insert(rank_buffer.end(), rank_bytes,
          rank_bytes + node.ranks() * sizeof(PathNode::rank_type));
        if(path_buffer.size() >= JOIN_IO_BUFFER_BYTES) { flush(path_buffer, path_fd, path_released, partial_path); }
        if(rank_buffer.size() >= JOIN_IO_BUFFER_BYTES) { flush(rank_buffer, rank_fd, rank_released, partial_rank); }
        written_paths++; written_ranks += node.ranks();
        emitted_paths++; emitted_ranks += node.ranks();
        if(readers[best]->read(records[best])) { queue.push(best); }
      }
      if(!compressed_path)
      {
        flush(path_buffer, path_fd, path_released, partial_path);
        flush(rank_buffer, rank_fd, rank_released, partial_rank);
      }
      close_output();
      if(emitted_paths != paths || emitted_ranks != ranks)
      {
        throw joinError("logical shard merge output count mismatch");
      }
      return true;
    }
    catch(...)
    {
      if(path_fd >= 0) { ::close(path_fd); }
      if(rank_fd >= 0) { ::close(rank_fd); }
      compressed_path.reset(); compressed_rank.reset();
      std::remove(partial_path.c_str()); std::remove(partial_rank.c_str());
      throw;
    }
  };

  PathGraph compacted(0, source.k(), source.step());
  compacted.delete_files = false; // It borrows retained source aliases until commit.
  std::vector<size_type> obsolete;
  size_type next_physical = 0;
  for(const auto& group : groups)
  {
    logical_file_id_t logical = group.first;
    std::vector<size_type> shards = group.second;
    if(shards.size() == 1)
    {
      size_type file = shards.front();
      compacted.path_names.push_back(source.path_names[file]); compacted.rank_names.push_back(source.rank_names[file]);
      compacted.path_counts.push_back(source.path_counts[file]); compacted.rank_counts.push_back(source.rank_counts[file]);
      compacted.logical_file_ids.push_back(logical); compacted.physical_shard_ids.push_back(physical_shard_id_t(next_physical++));
      compacted.path_count += source.path_counts[file]; compacted.rank_count += source.rank_counts[file];
      continue;
    }
    // Keep a bounded downstream frontier: one path/rank pair per retained
    // stream, divided fairly among logical groups. The batch width is then
    // derived from the group size and that target, but never exceeds the
    // reader/descriptor limit.
    const size_type target_pairs = std::max(static_cast<size_type>(1),
      (parameters.getMaxOpenFiles() - 4) / (2 * groups.size()));
    PathGraph current(0, source.k(), source.step()); current.delete_files = false;
    std::vector<size_type> origin;
    for(size_type source_file : shards)
    {
      current.path_names.push_back(source.path_names[source_file]); current.rank_names.push_back(source.rank_names[source_file]);
      current.path_counts.push_back(source.path_counts[source_file]); current.rank_counts.push_back(source.rank_counts[source_file]);
      current.logical_file_ids.push_back(logical); current.physical_shard_ids.push_back(source.physical_shard_ids[source_file]);
      current.path_count += source.path_counts[source_file]; current.rank_count += source.rank_counts[source_file];
      origin.push_back(source_file);
    }
    while(current.files() > target_pairs)
    {
      size_type batch_width = (current.files() + target_pairs - 1) / target_pairs;
      batch_width = std::min(batch_width, merge_fan_in);
      PathGraph pass(0, source.k(), source.step()); pass.delete_files = false;
      std::vector<size_type> next_origin;
      for(size_type first = 0; first < current.files(); first += batch_width)
      {
        size_type last = std::min(current.files(), first + batch_width);
        if(last - first == 1)
        {
          size_type file = first;
          pass.path_names.push_back(current.path_names[file]); pass.rank_names.push_back(current.rank_names[file]);
          pass.path_counts.push_back(current.path_counts[file]); pass.rank_counts.push_back(current.rank_counts[file]);
          pass.logical_file_ids.push_back(logical); pass.physical_shard_ids.push_back(current.physical_shard_ids[file]);
          pass.path_count += current.path_counts[file]; pass.rank_count += current.rank_counts[file];
          next_origin.push_back(origin[file]); continue;
        }
        std::vector<size_type> batch;
        for(size_type file = first; file < last; file++) { batch.push_back(file); }
        if(!merge_batch(current, batch, pass, logical, next_physical)) { return; }
        for(size_type file = first; file < last; file++)
        {
          if(origin[file] != std::numeric_limits<size_type>::max()) { obsolete.push_back(origin[file]); }
        }
        next_origin.push_back(std::numeric_limits<size_type>::max());
      }
      current.swap(pass); current.delete_files = false;
      origin.swap(next_origin);
    }
    for(size_type file = 0; file < current.files(); file++)
    {
      compacted.path_names.push_back(current.path_names[file]); compacted.rank_names.push_back(current.rank_names[file]);
      compacted.path_counts.push_back(current.path_counts[file]); compacted.rank_counts.push_back(current.rank_counts[file]);
      compacted.logical_file_ids.push_back(logical); compacted.physical_shard_ids.push_back(current.physical_shard_ids[file]);
      compacted.path_count += current.path_counts[file]; compacted.rank_count += current.rank_counts[file];
    }
  }
  // Intermediate generations stayed available until every logical group
  // succeeded. Retire only those not retained by the final frontier, charging
  // exactly the bytes whose unlink succeeded.
  for(std::string& name : new_outputs)
  {
    bool final_name = (std::find(compacted.path_names.begin(), compacted.path_names.end(), name) != compacted.path_names.end() ||
      std::find(compacted.rank_names.begin(), compacted.rank_names.end(), name) != compacted.rank_names.end());
    if(!final_name)
    {
      if(std::remove(name.c_str()) == 0)
      {
        staged_committed -= new_output_charges.at(name);
      }
    }
  }
  // Physical ids are artifact locators, not semantic graph ids. Intermediate
  // passes may mix retained source ids with newly allocated ids, so normalize
  // the published frontier to a deterministic collision-free sequence.
  for(size_type file = 0; file < compacted.files(); file++)
  {
    compacted.physical_shard_ids[file] = physical_shard_id_t(file);
  }
  // Install the complete frontier before retiring any predecessor. On every
  // earlier return/exception NewOutputCleanup removes only new outputs.
  source.delete_files = false; source.swap(compacted); source.delete_files = true;
  cleanup.keep = true; // The new frontier is now authoritative, even if cleanup leaks.
  for(size_type source_file : obsolete)
  {
    const size_type path_bytes = storedJoinBytes(compacted.path_names[source_file]);
    const size_type rank_bytes = storedJoinBytes(compacted.rank_names[source_file]);
    // Failed cleanup leaves an unreachable predecessor on disk. Keep those
    // bytes charged rather than risking deletion of the installed frontier.
    if(std::remove(compacted.path_names[source_file].c_str()) == 0) { staged_committed -= path_bytes; }
    if(std::remove(compacted.rank_names[source_file].c_str()) == 0) { staged_committed -= rank_bytes; }
  }
  committed_bytes = staged_committed;
  if(stats != nullptr)
  {
    stats->max_records_resident = std::max(stats->max_records_resident, merge_fan_in);
    stats->max_bytes_resident = std::max(stats->max_bytes_resident,
      merge_reservation);
  }
}

std::string
joinPartitionTaskName(const std::string& generation,
  logical_file_id_t logical, const JoinPartition& partition,
  std::uint64_t left_checksum, std::uint64_t right_checksum)
{
  return generation + "-join-l" + std::to_string(logical.value) +
    "-lb" + std::to_string(partition.left_begin) +
    "-le" + std::to_string(partition.left_end) +
    "-rb" + std::to_string(partition.right_begin) +
    "-re" + std::to_string(partition.right_end) +
    "-b" + std::to_string(partition.emit_bypass ? 1 : 0) +
    "-lc" + std::to_string(left_checksum) +
    "-rc" + std::to_string(right_checksum);
}

ArtifactIdentity
joinPartitionArtifact(const std::string& task, const std::string& name,
  const std::string& kind)
{
  return ArtifactIdentity(task, "join-partition", name, kind);
}

BuildWorkspace::ArtifactRef
checkpointJoinPartitionArtifact(BuildWorkspace& workspace,
  const ArtifactIdentity& identity, logical_file_id_t logical,
  const std::string& source, size_type records, size_type expected_bytes,
  size_type buffer_bytes, const std::string& sort_order)
{
  static_cast<void>(sort_order);
  // Worker outputs are immutable once their result marker is read. Adopt the
  // raw payload so a colocated workspace keeps the same inode; the workspace
  // primitive retains its bounded-copy fallback for other filesystems.
  return workspace.adopt_raw_payload(identity, logical, physical_shard_id_t(0),
    source, records, expected_bytes, buffer_bytes);
}

size_type
storedJoinBytes(const std::string& path)
{
  struct stat info;
  if(::stat(path.c_str(), &info) != 0 || info.st_size < 0 ||
     static_cast<std::uintmax_t>(info.st_size) >
       std::numeric_limits<size_type>::max())
  {
    throw joinError("cannot determine stored shard bytes", path);
  }
  return static_cast<size_type>(info.st_size);
}

void
checkpointJoinPartition(BuildWorkspace& workspace, const std::string& task,
  logical_file_id_t logical, const JoinPartition& partition,
  const std::string& path_name, const std::string& rank_name,
  size_type buffer_bytes)
{
  std::vector<BuildWorkspace::ArtifactRef> artifacts;
  artifacts.push_back(checkpointJoinPartitionArtifact(workspace,
    joinPartitionArtifact(task, "paths", "path-nodes-v1"), logical,
    path_name, partition.expected_paths,
    storedJoinBytes(path_name), buffer_bytes, "label"));
  artifacts.push_back(checkpointJoinPartitionArtifact(workspace,
    joinPartitionArtifact(task, "ranks", "path-ranks-v1"), logical,
    rank_name, partition.expected_ranks,
    storedJoinBytes(rank_name), buffer_bytes, "path-order"));
  workspace.commit_task(task, "join-partition", artifacts);
}

void
restoreJoinPartition(const BuildWorkspace& workspace, const std::string& task,
  logical_file_id_t logical, const JoinPartition& partition,
  const std::string& path_name,
  const std::string& rank_name, size_type buffer_bytes)
{
  workspace.restore_adopted_payload(
    joinPartitionArtifact(task, "paths", "path-nodes-v1"),
    logical, physical_shard_id_t(0), path_name, partition.expected_paths,
    storedJoinBytes(workspace.artifact_path(
      joinPartitionArtifact(task, "paths", "path-nodes-v1"), logical,
      physical_shard_id_t(0))), buffer_bytes, true);
  workspace.restore_adopted_payload(
    joinPartitionArtifact(task, "ranks", "path-ranks-v1"),
    logical, physical_shard_id_t(0), rank_name, partition.expected_ranks,
    storedJoinBytes(workspace.artifact_path(
      joinPartitionArtifact(task, "ranks", "path-ranks-v1"), logical,
      physical_shard_id_t(0))), buffer_bytes, true);
}

struct ActiveJoinWorker
{
  pid_t pid;
  size_type shard, maximum_stored_bytes;
  JoinPartition partition;
  std::string task_file, checkpoint_name;
  MemoryBudget::Reservation reservation;

  ActiveJoinWorker(pid_t process, size_type output_shard,
    const JoinPartition& range, const std::string& task,
    const std::string& checkpoint, size_type stored_limit,
    MemoryBudget::Reservation memory) :
    pid(process), shard(output_shard), maximum_stored_bytes(stored_limit),
    partition(range), task_file(task),
    checkpoint_name(checkpoint),
    reservation(std::move(memory)) { }

  ActiveJoinWorker(ActiveJoinWorker&&) = default;
  ActiveJoinWorker& operator=(ActiveJoinWorker&&) = default;

private:
  ActiveJoinWorker(const ActiveJoinWorker&);
  ActiveJoinWorker& operator=(const ActiveJoinWorker&);
};

pid_t
spawnJoinWorker(const std::string& executable, const std::string& task_file)
{
  std::vector<std::string> arguments;
  arguments.push_back(executable);
  arguments.push_back("gcsa-worker-task");
  arguments.push_back(task_file);
  std::vector<char*> argv;
  for(std::string& argument : arguments) { argv.push_back(&argument[0]); }
  argv.push_back(nullptr);
  pid_t child = -1;
  int result = ::posix_spawnp(&child, executable.c_str(), nullptr, nullptr,
    argv.data(), environ);
  if(result != 0)
  {
    throw joinError("cannot spawn external join worker: " +
      std::string(std::strerror(result)), executable);
  }
  return child;
}

void
validateWorkerOutput(const std::string& path, size_type expected_bytes)
{
  if(CompressedBlockReader::isFramed(path))
  {
    if(CompressedBlockReader::declaredLogicalSize(path) != expected_bytes)
    {
      throw joinError("compressed worker output length mismatch", path);
    }
    return;
  }
  struct stat info;
  if(::stat(path.c_str(), &info) != 0 || info.st_size < 0 ||
     static_cast<size_type>(info.st_size) != expected_bytes)
  {
    throw joinError("worker output length mismatch", path);
  }
}

struct ConcurrentPathBudgets
{
  size_type sort, join;

  ConcurrentPathBudgets(size_type sort_bytes, size_type join_bytes) :
    sort(sort_bytes), join(join_bytes) { }
};

ConcurrentPathBudgets
allocateConcurrentPathBudgets(const ConstructionParameters& parameters,
  size_type available_bytes, size_type minimum_join,
  const std::string& context)
{
  const size_type minimum_sort = externalPathGraphSortMinimumBudget();
  const size_type minimum_total = checkedJoinAdd(minimum_sort, minimum_join,
    "minimum concurrent path workspace");
  if(available_bytes < minimum_total)
  {
    throw joinError(context + " cannot admit concurrent label sorting and join blocking");
  }

  // Expert values are caps, while automatic values are preferred shares. An
  // automatic 75/25 split may place one side below its irreducible stream
  // buffers under tiny test budgets even though both phases fit together.
  // Rebalance that case instead of rejecting a feasible aggregate goal.
  const bool automatic_sort = parameters.sortRunSizeIsAutomatic();
  const bool automatic_join = parameters.joinPartitionSizeIsAutomatic();
  const size_type sort_cap = (automatic_sort ? available_bytes :
    parameters.getSortRunSize(available_bytes));
  const size_type join_cap = (automatic_join ? available_bytes :
    parameters.getJoinPartitionSize(available_bytes));
  if(sort_cap < minimum_sort)
  {
    throw joinError(context + " label-sort cap is below the implementation minimum");
  }
  if(join_cap < minimum_join)
  {
    throw joinError(context + " join cap is below the implementation minimum");
  }

  // Caps constrain allocations; they are not requests to reserve that many
  // bytes before considering the other concurrent phase. In particular,
  // setting both caps equal to the global memory ceiling must not hand almost
  // everything to the first phase and leave the join at its minimum. Begin at
  // the normal 75/25 operating point, then give capacity rejected by one cap
  // to the other phase. If both caps are restrictive, unused bytes remain
  // outside the working set.
  size_type join_bytes = std::max(minimum_join, available_bytes / 4);
  join_bytes = std::min(join_bytes, available_bytes - minimum_sort);
  size_type sort_bytes = available_bytes - join_bytes;
  sort_bytes = std::min(sort_bytes, sort_cap);
  join_bytes = std::min(join_bytes, join_cap);

  size_type remaining = available_bytes - sort_bytes - join_bytes;
  size_type add_sort = std::min(remaining, sort_cap - sort_bytes);
  sort_bytes += add_sort; remaining -= add_sort;
  size_type add_join = std::min(remaining, join_cap - join_bytes);
  join_bytes += add_join;
  return ConcurrentPathBudgets(sort_bytes, join_bytes);
}

void
runJoinWorkerPartitions(const JoinRun& left, const JoinRun& right,
  logical_file_id_t logical, const std::vector<JoinPartition>& partitions,
  size_type size_limit, const ConstructionParameters& parameters,
  size_type label_fan_in, PathGraph& next, size_type& next_physical,
  size_type& committed_bytes, ExternalPathJoinStats* stats,
  BuildWorkspace* workspace, const std::string& checkpoint_task)
{
  if(parameters.getWorkerExecutable().empty())
  {
    throw joinError("process workers require a worker executable");
  }
  size_type requested = std::min(parameters.getProcessWorkers(), partitions.size());
  size_type minimum_join = checkedJoinAdd(joinRawMinimumBudget(),
    joinRunReaderPairMemory(left, right),
    "minimum compressed worker join scan");
  size_type minimum_worker = checkedJoinAdd(externalPathGraphSortMinimumBudget(),
    minimum_join, "minimum worker reservation");
  // Leave room for the coordinator, allocator metadata, dynamic-library state,
  // and the checkpoint copy buffer. Tiny test budgets retain the same fraction
  // instead of relying on a fixed margin that would make them unusable.
  size_type safety_margin = parameters.getMemoryLimitBytes() / 16;
  size_type usable_memory = parameters.getMemoryLimitBytes() - safety_margin;
  size_type maximum_by_memory = usable_memory / minimum_worker;
  size_type concurrency = std::min(requested, maximum_by_memory);
  if(concurrency < 1)
  {
    throw joinError("memory limit cannot admit one external join worker process");
  }
  size_type worker_reservation = usable_memory / concurrency;
  ConcurrentPathBudgets worker_budgets = allocateConcurrentPathBudgets(parameters,
    worker_reservation, minimum_join, "per-worker workspace");
  size_type sort_budget = worker_budgets.sort;
  size_type join_block_budget = worker_budgets.join;
  size_type worker_threads = std::max(static_cast<size_type>(1),
    static_cast<size_type>(omp_get_max_threads()) / concurrency);
  TempFileCodecParameters worker_codec = parameters.getTempFileCodecParameters();
  // Process workers are the outer parallelism layer. Bound zstd's inner pool
  // before calculating the same framed-output limit the child will enforce.
  worker_codec.workers = std::min(worker_codec.workers, worker_threads);
  size_type checkpoint_buffer = std::min(parameters.getIOBufferSize(),
    std::max(static_cast<size_type>(KILOBYTE), worker_reservation / 8));

  size_type planned_bytes = 0;
  for(const JoinPartition& partition : partitions)
  {
    const size_type peak_bytes = externalPathGraphShardPeakBytes(
      partition.expected_paths, partition.expected_ranks,
      sort_budget, worker_codec);
    planned_bytes = checkedJoinAdd(planned_bytes, peak_bytes,
      "all worker peak output bytes");
    if(partition.expected_ranks > (static_cast<size_type>(1) << 40))
    {
      throw joinError("one join partition exceeds the 40-bit rank pointer range");
    }
  }
  if(planned_bytes > size_limit || committed_bytes > size_limit - planned_bytes)
  {
    throw joinError("configured disk limit exceeded by deterministic partition plan");
  }

  MemoryBudget memory(parameters.getMemoryLimitBytes(), safety_margin);
  std::vector<ActiveJoinWorker> active;
  active.reserve(concurrency);
  auto collect_one = [&]()
  {
    if(active.empty()) { return; }
    ActiveJoinWorker& worker = active.front();
    int status = 0;
    pid_t waited;
    do { waited = ::waitpid(worker.pid, &status, 0); }
    while(waited < 0 && errno == EINTR);
    worker.pid = -1;
    if(waited < 0 || !WIFEXITED(status) || WEXITSTATUS(status) != EXIT_SUCCESS)
    {
      throw joinError("external join worker failed", worker.task_file);
    }
    // The completed child no longer owns a working set. Return its tokens
    // before validating and checkpointing with the parent's bounded buffer.
    worker.reservation = MemoryBudget::Reservation();
    std::string result_file = worker.task_file + ".result";
    ExternalJoinWorkerResult result = decodeWorkerResult(readWorkerControl(result_file));
    if(result.max_bytes > worker_reservation)
    {
      throw joinError("worker exceeded its combined join/sort reservation", result_file);
    }
    if(result.paths != worker.partition.expected_paths ||
       result.ranks != worker.partition.expected_ranks ||
       result.bytes != worker.partition.expected_bytes)
    {
      throw joinError("worker result disagrees with partition plan", result_file);
    }
    validateWorkerOutput(next.path_names[worker.shard],
      result.paths * sizeof(PathNode));
    validateWorkerOutput(next.rank_names[worker.shard],
      result.ranks * sizeof(PathNode::rank_type));
    const size_type actual_stored_bytes = checkedJoinAdd(
      storedJoinBytes(next.path_names[worker.shard]),
      storedJoinBytes(next.rank_names[worker.shard]),
      "worker stored output bytes");
    if(result.stored_bytes != actual_stored_bytes ||
       result.stored_bytes > worker.maximum_stored_bytes)
    {
      throw joinError("worker stored output exceeds its admitted peak", result_file);
    }
    if(workspace != nullptr && !worker.checkpoint_name.empty())
    {
      checkpointJoinPartition(*workspace, worker.checkpoint_name, logical,
        worker.partition, next.path_names[worker.shard],
        next.rank_names[worker.shard], checkpoint_buffer);
    }
    next.path_counts[worker.shard] = result.paths;
    next.rank_counts[worker.shard] = result.ranks;
    next.path_count += result.paths; next.rank_count += result.ranks;
    committed_bytes = checkedJoinAdd(committed_bytes, result.stored_bytes,
      "committed worker output bytes");
    DiskIO::read_volume += result.bytes_read;
    DiskIO::write_volume += result.bytes_written;
    if(stats != nullptr)
    {
      stats->generated_records += result.generated;
      stats->sorted_bypass += result.bypassed;
      stats->direct_label_records += result.paths;
      stats->intermediate_path_bytes_avoided += result.bytes;
      stats->label_sort_runs += result.label_runs;
      stats->label_merge_passes = std::max(stats->label_merge_passes,
        result.label_merge_passes);
      stats->label_parallel_sorts = checkedJoinAdd(stats->label_parallel_sorts,
        result.label_parallel_sorts, "worker parallel label sorts");
      stats->grouped_expansion_records = checkedJoinAdd(
        stats->grouped_expansion_records, result.grouped_records,
        "worker grouped expansion records");
      stats->expansion_context_bytes_saved = checkedJoinAdd(
        stats->expansion_context_bytes_saved, result.context_bytes_saved,
        "worker expansion context bytes saved");
      stats->max_records_resident = std::max(stats->max_records_resident,
        checkedJoinMultiply(result.max_records, concurrency,
          "concurrent worker resident records"));
      stats->max_bytes_resident = std::max(stats->max_bytes_resident,
        checkedJoinMultiply(result.max_bytes, concurrency,
          "concurrent worker resident bytes"));
      stats->blocked_key_blocks = checkedJoinAdd(stats->blocked_key_blocks,
        result.blocks, "worker join blocks");
      stats->worker_processes++;
    }
    std::remove(result_file.c_str());
    TempFile::remove(worker.task_file);
    active.erase(active.begin());
  };

  try
  {
    for(const JoinPartition& partition : partitions)
    {
      const size_type peak_output_bytes = externalPathGraphShardPeakBytes(
        partition.expected_paths, partition.expected_ranks,
        sort_budget, worker_codec);
      size_type shard = appendOutputShard(next, logical,
        physical_shard_id_t(next_physical++));
      std::string partition_checkpoint = (checkpoint_task.empty() ? std::string() :
        joinPartitionTaskName(checkpoint_task, logical, partition,
          left.checksum, right.checksum));
      if(workspace != nullptr && !partition_checkpoint.empty() &&
         workspace->task_completed(partition_checkpoint, "join-partition"))
      {
        // Do not overlap restore buffers with fully admitted child working
        // sets. Completed tasks are cheap to restore and never respawn.
        while(!active.empty()) { collect_one(); }
        restoreJoinPartition(*workspace, partition_checkpoint, logical, partition,
          next.path_names[shard], next.rank_names[shard], checkpoint_buffer);
        validateWorkerOutput(next.path_names[shard],
          partition.expected_paths * sizeof(PathNode));
        validateWorkerOutput(next.rank_names[shard],
          partition.expected_ranks * sizeof(PathNode::rank_type));
        next.path_counts[shard] = partition.expected_paths;
        next.rank_counts[shard] = partition.expected_ranks;
        next.path_count += partition.expected_paths;
        next.rank_count += partition.expected_ranks;
        const size_type restored_bytes = checkedJoinAdd(
          storedJoinBytes(next.path_names[shard]),
          storedJoinBytes(next.rank_names[shard]),
          "restored worker output bytes");
        if(restored_bytes > peak_output_bytes)
        {
          throw joinError("restored worker output exceeds its admitted peak");
        }
        committed_bytes = checkedJoinAdd(committed_bytes, restored_bytes,
          "committed restored output bytes");
        if(stats != nullptr)
        {
          stats->restored_partitions++;
          stats->generated_records += partition.generated_paths;
          stats->sorted_bypass += partition.bypass_paths;
          stats->direct_label_records += partition.expected_paths;
          stats->intermediate_path_bytes_avoided += partition.expected_bytes;
        }
        continue;
      }
      while(active.size() >= concurrency) { collect_one(); }
      ExternalJoinWorkerTask task;
      task.logical = logical; task.left = left; task.right = right;
      task.partition = partition;
      task.output_path = next.path_names[shard];
      task.output_rank = next.rank_names[shard];
      task.temp_directory = TempFile::temp_dir;
      task.sort_budget = sort_budget; task.join_block_budget = join_block_budget;
      task.fan_in = label_fan_in;
      task.threads = worker_threads;
      task.codec = worker_codec;
      task.verify_payloads = parameters.getVerifyWorkspace();
      std::string task_file = TempFile::getName("gcsa_join_worker_task");
      writeWorkerControl(task_file, encodeWorkerTask(task));
      MemoryBudget::Reservation reservation = memory.reserve(worker_reservation,
        "join-worker-" + std::to_string(shard));
      pid_t child = -1;
      try
      {
        child = spawnJoinWorker(parameters.getWorkerExecutable(), task_file);
      }
      catch(...)
      {
        TempFile::remove(task_file);
        throw;
      }
      active.emplace_back(child, shard, partition, task_file,
        partition_checkpoint, peak_output_bytes,
        std::move(reservation));
    }
    while(!active.empty()) { collect_one(); }
  }
  catch(...)
  {
    for(ActiveJoinWorker& worker : active)
    {
      if(worker.pid > 0)
      {
        ::kill(worker.pid, SIGTERM);
        while(::waitpid(worker.pid, nullptr, 0) < 0 && errno == EINTR) { }
      }
      std::string result_file = worker.task_file + ".result";
      std::remove(result_file.c_str());
      TempFile::remove(worker.task_file);
    }
    throw;
  }
  if(stats != nullptr)
  {
    stats->join_partitions += partitions.size();
    stats->max_bytes_resident = std::max(stats->max_bytes_resident,
      static_cast<size_type>(memory.stats().maximum));
  }
}

} // namespace

size_type
externalPathJoinMinimumBudget()
{
  return joinRawMinimumBudget();
}

void
externalPathGraphExtend(PathGraph& graph, size_type size_limit,
  const ConstructionParameters& parameters, ExternalPathJoinStats* stats,
  BuildWorkspace* workspace, const std::string& checkpoint_task)
{
  if(stats != nullptr) { *stats = ExternalPathJoinStats(); }
  if(graph.logical_file_ids.size() != graph.files() ||
     graph.physical_shard_ids.size() != graph.files())
  {
    throw joinError("path graph shard identity metadata is incomplete");
  }
  size_type memory_budget = parameters.getMemoryLimitBytes();
  // Distribution and join-to-label expansion are different lifetimes. The
  // former shares memory only with one source shard reader, so forcing it into
  // the later concurrent join share can create hundreds of thousands of tiny
  // runs even when the complete distribution sort fits the global ceiling.
  // Give distribution its own phase budget, bounded by the explicit join cap
  // when present, and account the source decoder and codec inside it. The
  // label-sort cap belongs to the later output sorter; using it here would
  // break configurations where the two algorithms have different minima.
  size_type source_codec_bytes = maximumCompressedShardReaderBytes(graph);
  if(source_codec_bytes >= memory_budget)
  {
    throw joinError(
      "memory limit cannot admit one compressed source shard reader");
  }
  size_type maximum_distribution_sort = memory_budget - source_codec_bytes;
  if(!(parameters.joinPartitionSizeIsAutomatic()))
  {
    maximum_distribution_sort = std::min(maximum_distribution_sort,
      parameters.getJoinPartitionSize(memory_budget));
  }
  TempFileCodecParameters join_codec = boundedJoinRunCodec(
    parameters.getTempFileCodecParameters(), maximum_distribution_sort);
  if(maximum_distribution_sort < joinRunSortMinimumBudget(join_codec))
  {
    throw joinError("distribution sort cap is below the join-run minimum");
  }

  // Expansion reads two completed join runs while the label sorter is live;
  // only this lifetime needs the concurrent 75/25 split.
  size_type minimum_join = joinRunReadMinimumBudget(join_codec);
  ConcurrentPathBudgets direct_budgets = allocateConcurrentPathBudgets(parameters,
    memory_budget, minimum_join, "configured memory limit");
  size_type direct_sort_budget = direct_budgets.sort;
  size_type direct_join_budget = direct_budgets.join;
  size_type distribution_sort_budget = maximum_distribution_sort;
  if(stats != nullptr)
  {
    stats->distribution_sort_budget = distribution_sort_budget;
    stats->label_sort_budget = direct_sort_budget;
    stats->join_block_budget = direct_join_budget;
  }
  if(parameters.getMaxOpenFiles() < 6)
  {
    throw joinError("max-open-files must be at least 6 for the external join");
  }
  // During join-run compaction, two generated-output descriptors and one
  // merge-output descriptor are also open. Label sorting happens after those
  // writers close and can use one additional input descriptor.
  size_type join_fan_in = std::min(parameters.getMergeFanIn(),
    parameters.getMaxOpenFiles() - 3);
  size_type label_fan_in = std::min(parameters.getMergeFanIn(),
    parameters.getMaxOpenFiles() - 1);

  std::map<logical_file_id_t, std::vector<size_type>> groups;
  for(size_type file = 0; file < graph.files(); file++)
  {
    groups[graph.logicalFile(file)].push_back(file);
  }
  for(auto& group : groups)
  {
    std::sort(group.second.begin(), group.second.end(), [&graph](size_type left, size_type right)
    {
      if(graph.physicalShard(left) != graph.physicalShard(right))
      {
        return graph.physicalShard(left) < graph.physicalShard(right);
      }
      return left < right;
    });
  }

  PathGraph next(0, 2 * graph.k(), graph.step() + 1);
  size_type next_physical = 0, committed_bytes = 0;
  for(const auto& group : groups)
  {
    logical_file_id_t logical = group.first;

    ExternalJoinSorter right_sorter(logical, RIGHT_BY_FROM,
      distribution_sort_budget,
      join_fan_in, parameters.getVerifyWorkspace(), stats, join_codec);
    scanJoinSide(graph, group.second, logical, RIGHT_BY_FROM, right_sorter, stats);
    JoinRun right = right_sorter.finish();

    ExternalJoinSorter left_sorter(logical, LEFT_BY_TO,
      distribution_sort_budget,
      join_fan_in, parameters.getVerifyWorkspace(), stats, join_codec);
    scanJoinSide(graph, group.second, logical, LEFT_BY_TO, left_sorter, stats);
    JoinRun left = left_sorter.finish();

    if(stats != nullptr && source_codec_bytes > 0)
    {
      stats->max_bytes_resident = std::max(stats->max_bytes_resident,
        checkedJoinAdd(distribution_sort_budget, source_codec_bytes,
          "compressed distribution scan peak"));
    }

    if(parameters.getProcessWorkers() > 1)
    {
      std::vector<JoinPartition> partitions = planJoinPartitions(left, right,
        logical, parameters.getCheckpointBytes(), memory_budget,
        parameters.getVerifyWorkspace(), stats, workspace, checkpoint_task);
      runJoinWorkerPartitions(left, right, logical, partitions, size_limit,
        parameters, label_fan_in, next, next_physical, committed_bytes, stats,
        workspace, checkpoint_task);
    }
    else
    {
      size_type output_file = appendOutputShard(next, logical,
        physical_shard_id_t(next_physical++));
      // Both join sorters have released their run-sized buffers. Reuse that
      // working set for label run generation and feed expansion records
      // directly into it, avoiding a complete unsorted pair and reread.
      ExternalPathSortStats sort_stats;
      ExternalPathJoinStats direct_join_stats;
      ExternalPathSortSink output(next, output_file, direct_sort_budget,
        label_fan_in, size_limit, committed_bytes,
        (stats == nullptr ? nullptr : &sort_stats),
        parameters.getTempFileCodecParameters());
      joinSortedRuns(left, right, logical, direct_join_budget, output,
        parameters.getVerifyWorkspace(), (stats == nullptr ? nullptr : &direct_join_stats));
      output.finish();
      if(stats != nullptr)
      {
        stats->join_partitions++;
        stats->label_sort_runs += sort_stats.runs;
        stats->label_merge_passes = std::max(stats->label_merge_passes,
          sort_stats.merge_passes);
        stats->label_parallel_sorts += sort_stats.parallel_sorts;
        stats->grouped_expansion_records = checkedJoinAdd(
          stats->grouped_expansion_records, sort_stats.grouped_records,
          "direct grouped expansion records");
        stats->expansion_context_bytes_saved = checkedJoinAdd(
          stats->expansion_context_bytes_saved, sort_stats.context_bytes_saved,
          "direct expansion context bytes saved");
        stats->generated_records += direct_join_stats.generated_records;
        stats->sorted_bypass += direct_join_stats.sorted_bypass;
        stats->direct_label_records += direct_join_stats.direct_label_records;
        stats->intermediate_path_bytes_avoided += direct_join_stats.intermediate_path_bytes_avoided;
        stats->blocked_key_groups += direct_join_stats.blocked_key_groups;
        stats->blocked_key_blocks += direct_join_stats.blocked_key_blocks;
        stats->max_records_resident = std::max(stats->max_records_resident,
          checkedJoinAdd(sort_stats.max_records_resident,
            direct_join_stats.max_records_resident, "direct join resident records"));
        stats->max_bytes_resident = std::max(stats->max_bytes_resident,
          checkedJoinAdd(sort_stats.max_bytes_resident,
            direct_join_stats.max_bytes_resident, "direct join resident bytes"));
        // Label run generation overlaps with exactly two bounded join readers.
        // Include them in the reported peak instead of treating the sequential
        // sort and join budgets as independent maxima.
      }
    }
    removeJoinRun(left); removeJoinRun(right);
  }

  MemoryBudget compaction_memory(memory_budget, memory_budget / 16);
  compactLogicalJoinShards(next, size_limit, parameters, label_fan_in,
    compaction_memory, committed_bytes, stats);
  graph.clear(); graph.swap(next);
}

} // namespace gcsa
