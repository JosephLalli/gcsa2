/*
  Copyright (c) 2026 Jouni Siren and GCSA2 contributors

  External-memory prefix-doubling join. The implementation favors a strict
  resident-memory bound and restart-friendly immutable files over throughput.
*/

#include <gcsa/path_graph.h>
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
constexpr size_type JOIN_BASE_FIXED_BYTES = 32 * KILOBYTE;
constexpr size_type JOIN_LABEL_COUNT = PathLabel::LABEL_LENGTH + 1;
constexpr size_type JOIN_HEADER_BYTES = 8 + 4 + 4 + 4 + 8;
constexpr size_type JOIN_FOOTER_BYTES = 8 + 8 + 8;
constexpr size_type JOIN_RECORD_BYTES =
  4 + 3 * sizeof(node_type) + 4 + 8 + JOIN_LABEL_COUNT * sizeof(PathNode::rank_type);
// Join records are compact on disk, so issuing one system call per record can
// turn a sequential spill into syscall-bound I/O. Keep all phase-local I/O
// caches byte-sized and account for them in the same join reservation.
constexpr size_type JOIN_IO_BUFFER_BYTES = 64 * KILOBYTE;
constexpr size_type JOIN_IO_BUFFER_RECORDS = JOIN_IO_BUFFER_BYTES / JOIN_RECORD_BYTES;
constexpr size_type JOIN_PARALLEL_SORT_MIN_RECORDS = 64 * 1024;
// Run generation may simultaneously have a join writer, path/rank source
// caches, and path/rank output caches. Merge reader caches are accounted per
// input when selecting fan-in below.
constexpr size_type JOIN_FIXED_BYTES = JOIN_BASE_FIXED_BYTES + 5 * JOIN_IO_BUFFER_BYTES;
static_assert(JOIN_IO_BUFFER_RECORDS > 0, "join record exceeds the I/O buffer");
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

struct JoinRun
{
  std::string name;
  size_type records;

  JoinRun(const std::string& path = std::string(), size_type count = 0) :
    name(path), records(count) { }
};

std::runtime_error
joinError(const std::string& message, const std::string& path = std::string())
{
  return std::runtime_error("externalPathGraphExtend(): " + message +
    (path.empty() ? std::string() : ": " + path));
}

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
  if(order == 0 || order > PathLabel::LABEL_LENGTH || lcp > order)
  {
    throw joinError("invalid join record path metadata");
  }
  for(size_type i = 0; i < JOIN_LABEL_COUNT; i++)
  {
    target.labels[i] = decodeLittle<PathNode::rank_type>(in);
  }
}

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
  JoinFileWriter(const std::string& path, logical_file_id_t logical, JoinKeyKind kind) :
    name(path), logical_id(logical), key_kind(kind), descriptor(-1), count(0),
    checksum(1469598103934665603ULL), cache_released(0), finished(false), buffer()
  {
    this->descriptor = ::open(this->name.c_str(), O_CREAT | O_EXCL | O_WRONLY, 0644);
    if(this->descriptor < 0) { throw joinError("cannot create join run", this->name); }
    adviseSequential(this->descriptor);
    std::array<std::uint8_t, JOIN_HEADER_BYTES> header;
    std::uint8_t* out = header.data();
    encodeLittle<std::uint64_t>(out, JOIN_HEADER_MAGIC);
    encodeLittle<std::uint32_t>(out, JOIN_FORMAT_VERSION);
    encodeLittle<std::uint32_t>(out, static_cast<std::uint32_t>(this->key_kind));
    encodeLittle<std::uint32_t>(out, this->logical_id.value);
    encodeLittle<std::uint64_t>(out, 0);
    writeAll(this->descriptor, header.data(), header.size(), this->name);
    this->buffer.reserve(JOIN_IO_BUFFER_RECORDS);
  }

  ~JoinFileWriter()
  {
    if(this->descriptor >= 0) { ::close(this->descriptor); }
    if(!this->finished) { ::unlink(this->name.c_str()); }
  }

  void writeRecord(const JoinRecord& record)
  {
    std::array<std::uint8_t, JOIN_RECORD_BYTES> encoded;
    encodeJoinRecord(record, encoded);
    this->checksum = joinChecksum(encoded.data(), encoded.size(), this->checksum);
    this->buffer.push_back(encoded);
    this->count++;
    if(this->buffer.size() >= JOIN_IO_BUFFER_RECORDS) { this->flushRecords(); }
  }

  JoinRun finish()
  {
    if(this->finished) { throw joinError("join run already finished", this->name); }
    this->flushRecords();
    std::array<std::uint8_t, 8> encoded_count;
    std::uint8_t* count_out = encoded_count.data();
    encodeLittle<std::uint64_t>(count_out, this->count);
    pwriteAll(this->descriptor, encoded_count.data(), encoded_count.size(), 20, this->name);
    std::array<std::uint8_t, JOIN_FOOTER_BYTES> footer;
    std::uint8_t* out = footer.data();
    encodeLittle<std::uint64_t>(out, JOIN_FOOTER_MAGIC);
    encodeLittle<std::uint64_t>(out, this->count);
    encodeLittle<std::uint64_t>(out, this->checksum);
    writeAll(this->descriptor, footer.data(), footer.size(), this->name);
    trimWrittenCache(this->descriptor, this->cache_released, true, this->name);
    if(::close(this->descriptor) != 0) { throw joinError("cannot close join run", this->name); }
    this->descriptor = -1; this->finished = true;
    return JoinRun(this->name, this->count);
  }

private:
  void flushRecords()
  {
    if(this->buffer.empty()) { return; }
    writeAll(this->descriptor, this->buffer.data(),
      this->buffer.size() * JOIN_RECORD_BYTES, this->name);
    this->buffer.clear();
    off_t written = JOIN_HEADER_BYTES + this->count * JOIN_RECORD_BYTES;
    if(written - this->cache_released >= JOIN_CACHE_FLUSH_BYTES)
    {
      trimWrittenCache(this->descriptor, this->cache_released, false, this->name);
    }
  }

  JoinFileWriter(const JoinFileWriter&);
  JoinFileWriter& operator=(const JoinFileWriter&);

  std::string name;
  logical_file_id_t logical_id;
  JoinKeyKind key_kind;
  int descriptor;
  size_type count;
  std::uint64_t checksum;
  off_t cache_released;
  bool finished;
  std::vector<std::array<std::uint8_t, JOIN_RECORD_BYTES>> buffer;
};

class JoinFileReader
{
public:
  JoinFileReader(const JoinRun& run, logical_file_id_t logical, JoinKeyKind kind) :
    name(run.name), descriptor(-1), record_count(0), expected_checksum(0),
    cache_released(0), buffer(JOIN_IO_BUFFER_RECORDS), buffer_first(0),
    buffer_records(0)
  {
    this->descriptor = ::open(this->name.c_str(), O_RDONLY);
    if(this->descriptor < 0) { throw joinError("cannot open join run", this->name); }
    adviseSequential(this->descriptor);
    std::array<std::uint8_t, JOIN_HEADER_BYTES> header;
    preadAll(this->descriptor, header.data(), header.size(), 0, this->name);
    const std::uint8_t* in = header.data();
    std::uint64_t magic = decodeLittle<std::uint64_t>(in);
    std::uint32_t version = decodeLittle<std::uint32_t>(in);
    std::uint32_t actual_kind = decodeLittle<std::uint32_t>(in);
    std::uint32_t actual_logical = decodeLittle<std::uint32_t>(in);
    this->record_count = decodeLittle<std::uint64_t>(in);
    if(magic != JOIN_HEADER_MAGIC || version != JOIN_FORMAT_VERSION ||
       actual_kind != static_cast<std::uint32_t>(kind) || actual_logical != logical.value ||
       this->record_count != run.records)
    {
      throw joinError("join run header mismatch", this->name);
    }
    if(this->record_count > (std::numeric_limits<size_type>::max() -
       JOIN_HEADER_BYTES - JOIN_FOOTER_BYTES) / JOIN_RECORD_BYTES)
    {
      throw joinError("join run is too large for this build", this->name);
    }
    struct stat info;
    if(::fstat(this->descriptor, &info) != 0 || static_cast<size_type>(info.st_size) !=
       JOIN_HEADER_BYTES + this->record_count * JOIN_RECORD_BYTES + JOIN_FOOTER_BYTES)
    {
      throw joinError("join run length mismatch", this->name);
    }
    std::array<std::uint8_t, JOIN_FOOTER_BYTES> footer;
    preadAll(this->descriptor, footer.data(), footer.size(),
      JOIN_HEADER_BYTES + this->record_count * JOIN_RECORD_BYTES, this->name);
    in = footer.data();
    if(decodeLittle<std::uint64_t>(in) != JOIN_FOOTER_MAGIC ||
       decodeLittle<std::uint64_t>(in) != this->record_count)
    {
      throw joinError("join run footer mismatch", this->name);
    }
    this->expected_checksum = decodeLittle<std::uint64_t>(in);
  }

  ~JoinFileReader()
  {
    if(this->descriptor >= 0)
    {
      discardCachedRange(this->descriptor, 0, 0); ::close(this->descriptor);
    }
  }

  size_type size() const { return this->record_count; }

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
    for(size_type first = 0; first < this->record_count; first += this->buffer.size())
    {
      size_type records = std::min(this->buffer.size(), this->record_count - first);
      preadAll(this->descriptor, this->buffer.data(), records * JOIN_RECORD_BYTES,
        JOIN_HEADER_BYTES + first * JOIN_RECORD_BYTES, this->name);
      checksum = joinChecksum(this->buffer.data(), records * JOIN_RECORD_BYTES, checksum);
      trimReadCache(this->descriptor,
        JOIN_HEADER_BYTES + (first + records) * JOIN_RECORD_BYTES, this->cache_released);
    }
    if(checksum != this->expected_checksum) { throw joinError("join run checksum mismatch", this->name); }
    discardCachedRange(this->descriptor, 0, 0);
    this->cache_released = 0;
    this->buffer_records = 0;
  }

private:
  void refill(size_type first) const
  {
    this->buffer_first = first;
    this->buffer_records = std::min(this->buffer.size(), this->record_count - first);
    preadAll(this->descriptor, this->buffer.data(),
      this->buffer_records * JOIN_RECORD_BYTES,
      JOIN_HEADER_BYTES + first * JOIN_RECORD_BYTES, this->name);
    trimReadCache(this->descriptor,
      JOIN_HEADER_BYTES + (first + this->buffer_records) * JOIN_RECORD_BYTES,
      this->cache_released);
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
  JoinKeyKind kind, bool verify_payloads, ExternalPathJoinStats* stats)
{
  if(inputs.empty()) { throw joinError("cannot merge an empty run set"); }
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
  JoinFileWriter output(output_name, logical, kind);
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
      JOIN_FIXED_BYTES + inputs.size() *
      (JOIN_IO_BUFFER_BYTES + sizeof(JoinRecord) + sizeof(size_type) +
       sizeof(std::unique_ptr<JoinFileReader>)));
  }
  return output.finish();
}

class ExternalJoinSorter
{
public:
  ExternalJoinSorter(logical_file_id_t logical, JoinKeyKind kind,
    size_type byte_budget, size_type requested_fan_in, bool verify_payloads,
    ExternalPathJoinStats* stats) :
    logical_id(logical), key_kind(kind), budget(byte_budget), fan_in(0),
    run_records(0), verify_runs(verify_payloads), statistics(stats)
  {
    if(this->budget < externalPathJoinMinimumBudget())
    {
      throw joinError("memory budget is too small for an external join");
    }
    size_type per_input = sizeof(JoinRecord) + sizeof(size_type) +
      sizeof(std::unique_ptr<JoinFileReader>) + JOIN_IO_BUFFER_BYTES + 128;
    size_type maximum_fan_in = (this->budget - JOIN_FIXED_BYTES) / per_input;
    this->fan_in = std::min(std::max(static_cast<size_type>(2), requested_fan_in),
      maximum_fan_in);
    if(this->fan_in < 2) { throw joinError("memory budget cannot support a two-way join merge"); }
    this->run_records = std::max(static_cast<size_type>(1),
      (this->budget - JOIN_FIXED_BYTES) / (2 * sizeof(JoinRecord)));
    this->buffer.reserve(this->run_records);
    if(this->statistics != nullptr)
    {
      this->statistics->max_records_resident = std::max(
        this->statistics->max_records_resident, this->run_records);
      this->statistics->max_bytes_resident = std::max(
        this->statistics->max_bytes_resident,
        JOIN_FIXED_BYTES + this->run_records * sizeof(JoinRecord));
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
      JoinFileWriter empty(empty_name, this->logical_id, this->key_kind);
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
        JoinRun merged = mergeJoinRuns(group, this->logical_id, this->key_kind,
          this->verify_runs, this->statistics);
        for(size_type i = 0; i < group.size(); i++) { TempFile::remove(group[i].name); }
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
    JoinFileWriter writer(name, this->logical_id, this->key_kind);
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
    JoinRun merged = mergeJoinRuns(inputs, this->logical_id, this->key_kind,
      this->verify_runs, this->statistics);
    for(size_type i = 0; i < inputs.size(); i++) { TempFile::remove(inputs[i].name); }
    this->addRun(merged, level + 1);
  }

  logical_file_id_t logical_id;
  JoinKeyKind key_kind;
  size_type budget, fan_in, run_records;
  bool verify_runs;
  ExternalPathJoinStats* statistics;
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
    rank_buffer_first(0), rank_buffer_records(0)
  {
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
    if(record.node.pointer() != this->rank_offset || record.node.order() == 0 ||
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
    preadAll(this->path_descriptor, this->path_buffer.data(),
      this->path_buffer_records * sizeof(PathNode),
      this->path_buffer_first * sizeof(PathNode), this->path_name);
    trimReadCache(this->path_descriptor,
      (this->path_buffer_first + this->path_buffer_records) * sizeof(PathNode),
      this->path_cache_released);
  }

  void refillRanks()
  {
    this->rank_buffer_first = this->rank_offset;
    this->rank_buffer_records = std::min(this->rank_buffer.size(),
      this->rank_count - this->rank_buffer_first);
    preadAll(this->rank_descriptor, this->rank_buffer.data(),
      this->rank_buffer_records * sizeof(PathNode::rank_type),
      this->rank_buffer_first * sizeof(PathNode::rank_type), this->rank_name);
    trimReadCache(this->rank_descriptor,
      (this->rank_buffer_first + this->rank_buffer_records) * sizeof(PathNode::rank_type),
      this->rank_cache_released);
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
  JoinFileReader left(left_run, logical, LEFT_BY_TO);
  JoinFileReader right(right_run, logical, RIGHT_BY_FROM);
  if(verify_payloads) { left.validate(); right.validate(); }
  left_end = std::min(left_end, left.size());
  right_end = std::min(right_end, right.size());
  if(left_begin > left_end || right_begin > right_end)
  {
    throw joinError("invalid worker join range");
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
         (byte_budget - JOIN_FIXED_BYTES) / sizeof(JoinRecord)))
    {
      stats->blocked_key_groups++;
    }
    // This is the unsplittable heavy-key fallback. Replaying the right range
    // for each left record is I/O-expensive, but resident path data is constant.
    for(size_type i = left_offset; i < left_limit; i++)
    {
      left.read(i, left_record);
      for(size_type j = right_offset; j < right_limit; j++)
      {
        right.read(j, right_record);
        JoinRecord generated = extendRecord(left_record, right_record);
        writeLabelRecord(generated, output, false, stats);
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

struct JoinGroupSummary
{
  node_type key;
  size_type begin, end, count, order_sum;
  size_type bypass_paths, bypass_ranks;

  JoinGroupSummary() : key(0), begin(0), end(0), count(0), order_sum(0),
    bypass_paths(0), bypass_ranks(0) { }
};

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
readJoinGroup(const JoinFileReader& reader, size_type begin, size_type limit)
{
  if(begin >= limit) { throw joinError("cannot inspect an empty join-key group"); }
  JoinGroupSummary result;
  result.begin = begin;
  JoinRecord record;
  reader.read(begin, record); result.key = record.key;
  size_type offset = begin;
  while(offset < limit)
  {
    reader.read(offset, record);
    if(record.key != result.key) { break; }
    result.count++;
    result.order_sum = checkedJoinAdd(result.order_sum, record.node.order(),
      "group order sum");
    if(record.node.sorted())
    {
      result.bypass_paths++;
      result.bypass_ranks = checkedJoinAdd(result.bypass_ranks,
        record.node.ranks(), "bypass ranks");
    }
    offset++;
  }
  result.end = offset;
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
readJoinSubgroup(const JoinFileReader& reader, const JoinGroupSummary& parent,
  size_type begin, size_type end)
{
  if(begin >= end || begin < parent.begin || end > parent.end)
  {
    throw joinError("invalid join-key subgroup range");
  }
  JoinGroupSummary result;
  result.key = parent.key; result.begin = begin; result.end = end;
  JoinRecord record;
  for(size_type offset = begin; offset < end; offset++)
  {
    reader.read(offset, record);
    if(record.key != parent.key)
    {
      throw joinError("join-key subgroup crosses a key boundary");
    }
    result.count++;
    result.order_sum = checkedJoinAdd(result.order_sum, record.node.order(),
      "subgroup order sum");
    if(record.node.sorted())
    {
      result.bypass_paths++;
      result.bypass_ranks = checkedJoinAdd(result.bypass_ranks,
        record.node.ranks(), "subgroup bypass ranks");
    }
  }
  return result;
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
planJoinGroupRecursive(const JoinFileReader& left_reader,
  const JoinFileReader& right_reader, const JoinGroupSummary& left,
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
      left.begin, middle);
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
      right.begin, middle);
    JoinGroupSummary second = remainingJoinSubgroup(right, first);
    planJoinGroupRecursive(left_reader, right_reader, left, first,
      target_bytes, emit_bypass, result, stats);
    planJoinGroupRecursive(left_reader, right_reader, left, second,
      target_bytes, emit_bypass, result, stats);
  }
}

std::vector<JoinPartition>
planJoinPartitions(const JoinRun& left_run, const JoinRun& right_run,
  logical_file_id_t logical, size_type target_bytes, bool verify_payloads,
  ExternalPathJoinStats* stats)
{
  JoinFileReader left(left_run, logical, LEFT_BY_TO);
  JoinFileReader right(right_run, logical, RIGHT_BY_FROM);
  if(verify_payloads) { left.validate(); right.validate(); }
  target_bytes = std::max(static_cast<size_type>(1), target_bytes);

  std::vector<JoinPartition> result;
  size_type left_offset = 0, right_offset = 0;
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
  while(left_offset < left.size() || right_offset < right.size())
  {
    if(!have_left && left_offset < left.size())
    {
      left_group = readJoinGroup(left, left_offset, left.size()); have_left = true;
    }
    if(!have_right && right_offset < right.size())
    {
      right_group = readJoinGroup(right, right_offset, right.size()); have_right = true;
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
    if(oversized_group && selected_right != nullptr &&
       (selected_left != nullptr || selected_right->bypass_paths > 0))
    {
      flush_current();
      JoinGroupSummary empty_left;
      empty_left.key = selected_right->key;
      empty_left.begin = left_offset; empty_left.end = left_offset;
      planJoinGroupRecursive(left, right,
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
constexpr std::uint32_t WORKER_FORMAT_VERSION = 1;
constexpr size_type WORKER_CONTROL_LIMIT = MEGABYTE;

struct ExternalJoinWorkerTask
{
  logical_file_id_t logical;
  JoinRun left, right;
  JoinPartition partition;
  std::string output_path, output_rank, temp_directory;
  size_type sort_budget, fan_in, threads;
  bool verify_payloads;
};

struct ExternalJoinWorkerResult
{
  size_type paths, ranks, bytes;
  size_type generated, bypassed, label_runs, label_merge_passes, label_parallel_sorts;
  size_type max_records, max_bytes, bytes_read, bytes_written;

  ExternalJoinWorkerResult() : paths(0), ranks(0), bytes(0), generated(0),
    bypassed(0), label_runs(0), label_merge_passes(0), label_parallel_sorts(0),
    max_records(0), max_bytes(0), bytes_read(0), bytes_written(0) { }
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
  appendWorkerValue<std::uint64_t>(data, task.fan_in);
  appendWorkerValue<std::uint64_t>(data, task.threads);
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
  task.fan_in = readWorkerValue<std::uint64_t>(data, offset);
  task.threads = readWorkerValue<std::uint64_t>(data, offset);
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
  appendWorkerValue<std::uint64_t>(data, result.generated);
  appendWorkerValue<std::uint64_t>(data, result.bypassed);
  appendWorkerValue<std::uint64_t>(data, result.label_runs);
  appendWorkerValue<std::uint64_t>(data, result.label_merge_passes);
  appendWorkerValue<std::uint64_t>(data, result.label_parallel_sorts);
  appendWorkerValue<std::uint64_t>(data, result.max_records);
  appendWorkerValue<std::uint64_t>(data, result.max_bytes);
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
  result.generated = readWorkerValue<std::uint64_t>(data, offset);
  result.bypassed = readWorkerValue<std::uint64_t>(data, offset);
  result.label_runs = readWorkerValue<std::uint64_t>(data, offset);
  result.label_merge_passes = readWorkerValue<std::uint64_t>(data, offset);
  result.label_parallel_sorts = readWorkerValue<std::uint64_t>(data, offset);
  result.max_records = readWorkerValue<std::uint64_t>(data, offset);
  result.max_bytes = readWorkerValue<std::uint64_t>(data, offset);
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
    ExternalPathSortSink sink(output, 0, task.sort_budget, task.fan_in,
      task.partition.expected_bytes, committed_bytes, &sort_stats);
    joinSortedRuns(task.left, task.right, task.logical, task.sort_budget,
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
    result.bytes = sink.bytes(); result.generated = join_stats.generated_records;
    result.bypassed = join_stats.sorted_bypass; result.label_runs = sort_stats.runs;
    result.label_merge_passes = sort_stats.merge_passes;
    result.label_parallel_sorts = sort_stats.parallel_sorts;
    result.max_records = sort_stats.max_records_resident;
    result.max_bytes = sort_stats.max_bytes_resident +
      2 * (JOIN_IO_BUFFER_BYTES + sizeof(JoinRecord) + 256);
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

std::string
joinPartitionTaskName(const std::string& generation,
  logical_file_id_t logical, const JoinPartition& partition)
{
  return generation + "-join-l" + std::to_string(logical.value) +
    "-lb" + std::to_string(partition.left_begin) +
    "-le" + std::to_string(partition.left_end) +
    "-rb" + std::to_string(partition.right_begin) +
    "-re" + std::to_string(partition.right_end) +
    "-b" + std::to_string(partition.emit_bypass ? 1 : 0);
}

ArtifactIdentity
joinPartitionArtifact(const std::string& task, const std::string& name,
  const std::string& kind)
{
  return ArtifactIdentity(task, "join-partition", name, kind);
}

BuildWorkspace::ArtifactRef
copyJoinPartitionArtifact(BuildWorkspace& workspace,
  const ArtifactIdentity& identity, logical_file_id_t logical,
  const std::string& source, size_type records, size_type expected_bytes,
  size_type buffer_bytes, const std::string& sort_order)
{
  std::ifstream input;
  input.rdbuf()->pubsetbuf(nullptr, 0);
  input.open(source.c_str(), std::ios_base::binary);
  if(!input) { throw joinError("cannot checkpoint worker output", source); }
  input.seekg(0, std::ios_base::end);
  std::streamoff end = input.tellg();
  input.seekg(0, std::ios_base::beg);
  if(end < 0 || static_cast<size_type>(end) != expected_bytes)
  {
    throw joinError("worker output changed before checkpoint", source);
  }
  BuildWorkspace::ArtifactWriter writer = workspace.open_artifact(identity,
    logical, physical_shard_id_t(0), sort_order, "all");
  std::vector<std::uint8_t> buffer(std::max(static_cast<size_type>(1),
    buffer_bytes));
  size_type remaining = expected_bytes;
  while(remaining > 0)
  {
    size_type bytes = std::min(remaining, buffer.size());
    input.read(reinterpret_cast<char*>(buffer.data()), bytes);
    if(input.gcount() != static_cast<std::streamsize>(bytes))
    {
      throw joinError("short read while checkpointing worker output", source);
    }
    writer.write(buffer.data(), bytes); remaining -= bytes;
  }
  return writer.finish(records);
}

void
checkpointJoinPartition(BuildWorkspace& workspace, const std::string& task,
  logical_file_id_t logical, const JoinPartition& partition,
  const std::string& path_name, const std::string& rank_name,
  size_type buffer_bytes)
{
  std::vector<BuildWorkspace::ArtifactRef> artifacts;
  artifacts.push_back(copyJoinPartitionArtifact(workspace,
    joinPartitionArtifact(task, "paths", "path-nodes-v1"), logical,
    path_name, partition.expected_paths,
    checkedJoinMultiply(partition.expected_paths, sizeof(PathNode),
      "checkpoint path bytes"), buffer_bytes, "label"));
  artifacts.push_back(copyJoinPartitionArtifact(workspace,
    joinPartitionArtifact(task, "ranks", "path-ranks-v1"), logical,
    rank_name, partition.expected_ranks,
    checkedJoinMultiply(partition.expected_ranks, sizeof(PathNode::rank_type),
      "checkpoint rank bytes"), buffer_bytes, "path-order"));
  workspace.commit_task(task, "join-partition", artifacts);
}

void
restoreJoinPartition(const BuildWorkspace& workspace, const std::string& task,
  logical_file_id_t logical, const std::string& path_name,
  const std::string& rank_name, size_type buffer_bytes)
{
  workspace.restore_artifact(joinPartitionArtifact(task, "paths", "path-nodes-v1"),
    logical, physical_shard_id_t(0), path_name, buffer_bytes);
  workspace.restore_artifact(joinPartitionArtifact(task, "ranks", "path-ranks-v1"),
    logical, physical_shard_id_t(0), rank_name, buffer_bytes);
}

struct ActiveJoinWorker
{
  pid_t pid;
  size_type shard;
  JoinPartition partition;
  std::string task_file, checkpoint_name;
  MemoryBudget::Reservation reservation;

  ActiveJoinWorker(pid_t process, size_type output_shard,
    const JoinPartition& range, const std::string& task,
    const std::string& checkpoint,
    MemoryBudget::Reservation memory) :
    pid(process), shard(output_shard), partition(range), task_file(task),
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
  struct stat info;
  if(::stat(path.c_str(), &info) != 0 || info.st_size < 0 ||
     static_cast<size_type>(info.st_size) != expected_bytes)
  {
    throw joinError("worker output length mismatch", path);
  }
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
  size_type minimum_worker = externalPathGraphSortMinimumBudget() +
    2 * (JOIN_IO_BUFFER_BYTES + sizeof(JoinRecord) + 256);
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
  size_type sort_budget = std::min(parameters.getSortRunSize(),
    worker_reservation - 2 * (JOIN_IO_BUFFER_BYTES + sizeof(JoinRecord) + 256));
  if(sort_budget < externalPathGraphSortMinimumBudget())
  {
    throw joinError("per-worker label-sort budget is below the implementation minimum");
  }
  size_type worker_threads = std::max(static_cast<size_type>(1),
    static_cast<size_type>(omp_get_max_threads()) / concurrency);
  size_type checkpoint_buffer = std::min(parameters.getIOBufferSize(),
    std::max(static_cast<size_type>(KILOBYTE), worker_reservation / 8));

  size_type planned_bytes = 0;
  for(const JoinPartition& partition : partitions)
  {
    planned_bytes = checkedJoinAdd(planned_bytes, partition.expected_bytes,
      "all worker output bytes");
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
    if(workspace != nullptr && !worker.checkpoint_name.empty())
    {
      checkpointJoinPartition(*workspace, worker.checkpoint_name, logical,
        worker.partition, next.path_names[worker.shard],
        next.rank_names[worker.shard], checkpoint_buffer);
    }
    next.path_counts[worker.shard] = result.paths;
    next.rank_counts[worker.shard] = result.ranks;
    next.path_count += result.paths; next.rank_count += result.ranks;
    committed_bytes += result.bytes;
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
      stats->max_records_resident = std::max(stats->max_records_resident,
        checkedJoinMultiply(result.max_records, concurrency,
          "concurrent worker resident records"));
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
      size_type shard = appendOutputShard(next, logical,
        physical_shard_id_t(next_physical++));
      std::string partition_checkpoint = (checkpoint_task.empty() ? std::string() :
        joinPartitionTaskName(checkpoint_task, logical, partition));
      if(workspace != nullptr && !partition_checkpoint.empty() &&
         workspace->task_completed(partition_checkpoint, "join-partition"))
      {
        // Do not overlap restore buffers with fully admitted child working
        // sets. Completed tasks are cheap to restore and never respawn.
        while(!active.empty()) { collect_one(); }
        restoreJoinPartition(*workspace, partition_checkpoint, logical,
          next.path_names[shard], next.rank_names[shard], checkpoint_buffer);
        validateWorkerOutput(next.path_names[shard],
          partition.expected_paths * sizeof(PathNode));
        validateWorkerOutput(next.rank_names[shard],
          partition.expected_ranks * sizeof(PathNode::rank_type));
        next.path_counts[shard] = partition.expected_paths;
        next.rank_counts[shard] = partition.expected_ranks;
        next.path_count += partition.expected_paths;
        next.rank_count += partition.expected_ranks;
        committed_bytes += partition.expected_bytes;
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
      task.sort_budget = sort_budget; task.fan_in = label_fan_in;
      task.threads = worker_threads;
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
        partition_checkpoint,
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
  // The smallest useful task must support two buffered merge inputs and a
  // bounded output writer, not merely a handful of in-memory records.
  size_type per_input = JOIN_IO_BUFFER_BYTES + sizeof(JoinRecord) +
    sizeof(size_type) + sizeof(std::unique_ptr<JoinFileReader>) + 128;
  return JOIN_FIXED_BYTES + 2 * per_input;
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
  if(memory_budget < externalPathJoinMinimumBudget() ||
     memory_budget < externalPathGraphSortMinimumBudget())
  {
    throw joinError("configured memory limit is below the external path minimum");
  }
  size_type join_budget = std::min(memory_budget, parameters.getJoinPartitionSize());
  size_type sort_budget = std::min(memory_budget, parameters.getSortRunSize());
  if(join_budget < externalPathJoinMinimumBudget() ||
     sort_budget < externalPathGraphSortMinimumBudget())
  {
    throw joinError("join partition or sort run budget is below the implementation minimum");
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

    ExternalJoinSorter right_sorter(logical, RIGHT_BY_FROM, join_budget,
      join_fan_in, parameters.getVerifyWorkspace(), stats);
    scanJoinSide(graph, group.second, logical, RIGHT_BY_FROM, right_sorter, stats);
    JoinRun right = right_sorter.finish();

    ExternalJoinSorter left_sorter(logical, LEFT_BY_TO, join_budget,
      join_fan_in, parameters.getVerifyWorkspace(), stats);
    scanJoinSide(graph, group.second, logical, LEFT_BY_TO, left_sorter, stats);
    JoinRun left = left_sorter.finish();

    if(parameters.getProcessWorkers() > 1)
    {
      std::vector<JoinPartition> partitions = planJoinPartitions(left, right,
        logical, parameters.getCheckpointBytes(), parameters.getVerifyWorkspace(), stats);
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
      ExternalPathSortSink output(next, output_file, sort_budget,
        label_fan_in, size_limit, committed_bytes,
        (stats == nullptr ? nullptr : &sort_stats));
      joinSortedRuns(left, right, logical, join_budget, output,
        parameters.getVerifyWorkspace(), stats);
      output.finish();
      if(stats != nullptr)
      {
        stats->join_partitions++;
        stats->label_sort_runs += sort_stats.runs;
        stats->label_merge_passes = std::max(stats->label_merge_passes,
          sort_stats.merge_passes);
        stats->label_parallel_sorts += sort_stats.parallel_sorts;
        stats->max_records_resident = std::max(stats->max_records_resident,
          sort_stats.max_records_resident);
        stats->max_bytes_resident = std::max(stats->max_bytes_resident,
          sort_stats.max_bytes_resident);
        // Label run generation overlaps with exactly two bounded join readers.
        // Include them in the reported peak instead of treating the sequential
        // sort and join budgets as independent maxima.
        size_type direct_pipeline_bytes = sort_stats.max_bytes_resident +
          2 * (JOIN_IO_BUFFER_BYTES + sizeof(JoinRecord) + 256);
        stats->max_bytes_resident = std::max(stats->max_bytes_resident,
          direct_pipeline_bytes);
      }
    }
    TempFile::remove(left.name); TempFile::remove(right.name);
  }

  graph.clear(); graph.swap(next);
}

} // namespace gcsa
