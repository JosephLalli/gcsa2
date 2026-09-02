/*
  Copyright (c) 2026 Jouni Siren and GCSA2 contributors

  External-memory prefix-doubling join. The implementation favors a strict
  resident-memory bound and restart-friendly immutable files over throughput.
*/

#include <gcsa/path_graph.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <cstdint>
#include <fcntl.h>
#include <limits>
#include <map>
#include <memory>
#include <queue>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>

namespace gcsa
{

namespace
{

constexpr std::uint64_t JOIN_HEADER_MAGIC = 0x314e494f4a534347ULL; // "GCSJOIN1"
constexpr std::uint64_t JOIN_FOOTER_MAGIC = 0x31444e454a534347ULL; // "GCSJEND1"
constexpr std::uint32_t JOIN_FORMAT_VERSION = 1;
constexpr size_type JOIN_FIXED_BYTES = 32 * KILOBYTE;
constexpr size_type JOIN_LABEL_COUNT = PathLabel::LABEL_LENGTH + 1;
constexpr size_type JOIN_HEADER_BYTES = 8 + 4 + 4 + 4 + 8;
constexpr size_type JOIN_FOOTER_BYTES = 8 + 8 + 8;
constexpr size_type JOIN_RECORD_BYTES =
  4 + 3 * sizeof(node_type) + 4 + 8 + JOIN_LABEL_COUNT * sizeof(PathNode::rank_type);
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
    checksum(1469598103934665603ULL), cache_released(0), finished(false)
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
  }

  ~JoinFileWriter()
  {
    if(this->descriptor >= 0) { ::close(this->descriptor); }
    if(!this->finished) { ::unlink(this->name.c_str()); }
  }

  void writeRecord(const JoinRecord& record)
  {
    std::array<std::uint8_t, JOIN_RECORD_BYTES> buffer;
    encodeJoinRecord(record, buffer);
    writeAll(this->descriptor, buffer.data(), buffer.size(), this->name);
    this->checksum = joinChecksum(buffer.data(), buffer.size(), this->checksum);
    this->count++;
    off_t written = JOIN_HEADER_BYTES + this->count * JOIN_RECORD_BYTES;
    if(written - this->cache_released >= JOIN_CACHE_FLUSH_BYTES)
    {
      trimWrittenCache(this->descriptor, this->cache_released, false, this->name);
    }
  }

  JoinRun finish()
  {
    if(this->finished) { throw joinError("join run already finished", this->name); }
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
};

class JoinFileReader
{
public:
  JoinFileReader(const JoinRun& run, logical_file_id_t logical, JoinKeyKind kind) :
    name(run.name), descriptor(-1), record_count(0), expected_checksum(0),
    cache_released(0)
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
    std::array<std::uint8_t, JOIN_RECORD_BYTES> buffer;
    preadAll(this->descriptor, buffer.data(), buffer.size(),
      JOIN_HEADER_BYTES + index * JOIN_RECORD_BYTES, this->name);
    decodeJoinRecord(buffer, record);
    trimReadCache(this->descriptor,
      JOIN_HEADER_BYTES + (index + 1) * JOIN_RECORD_BYTES, this->cache_released);
  }

  void validate() const
  {
    std::uint64_t checksum = 1469598103934665603ULL;
    std::array<std::uint8_t, JOIN_RECORD_BYTES> buffer;
    for(size_type i = 0; i < this->record_count; i++)
    {
      preadAll(this->descriptor, buffer.data(), buffer.size(),
        JOIN_HEADER_BYTES + i * JOIN_RECORD_BYTES, this->name);
      checksum = joinChecksum(buffer.data(), buffer.size(), checksum);
      trimReadCache(this->descriptor,
        JOIN_HEADER_BYTES + (i + 1) * JOIN_RECORD_BYTES, this->cache_released);
    }
    if(checksum != this->expected_checksum) { throw joinError("join run checksum mismatch", this->name); }
    discardCachedRange(this->descriptor, 0, 0);
    this->cache_released = 0;
  }

private:
  JoinFileReader(const JoinFileReader&);
  JoinFileReader& operator=(const JoinFileReader&);

  std::string name;
  int descriptor;
  size_type record_count;
  std::uint64_t expected_checksum;
  mutable off_t cache_released;
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
      (sizeof(JoinRecord) + sizeof(size_type) + sizeof(std::unique_ptr<JoinFileReader>)));
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
      sizeof(std::unique_ptr<JoinFileReader>) + 128;
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
    std::sort(this->buffer.begin(), this->buffer.end(), joinRecordLess);
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
    path_cache_released(0), rank_cache_released(0)
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
    preadAll(this->path_descriptor, &record.node, sizeof(PathNode),
      this->path_offset * sizeof(PathNode), this->path_name);
    if(record.node.pointer() != this->rank_offset || record.node.order() == 0 ||
       record.node.order() > PathLabel::LABEL_LENGTH || record.node.lcp() > record.node.order() ||
       record.node.ranks() > this->rank_count - this->rank_offset)
    {
      throw joinError("invalid path or rank pointer", this->path_name);
    }
    preadAll(this->rank_descriptor, record.labels,
      record.node.ranks() * sizeof(PathNode::rank_type),
      this->rank_offset * sizeof(PathNode::rank_type), this->rank_name);
    for(size_type i = record.node.ranks(); i < JOIN_LABEL_COUNT; i++)
    {
      record.labels[i] = PathLabel::NO_RANK;
    }
    this->rank_offset += record.node.ranks(); this->path_offset++;
    record.node.setPointer(0);
    trimReadCache(this->path_descriptor,
      this->path_offset * sizeof(PathNode), this->path_cache_released);
    trimReadCache(this->rank_descriptor,
      this->rank_offset * sizeof(PathNode::rank_type), this->rank_cache_released);
    return true;
  }

private:
  PathShardReader(const PathShardReader&);
  PathShardReader& operator=(const PathShardReader&);

  std::string path_name, rank_name;
  int path_descriptor, rank_descriptor;
  size_type path_count, rank_count, path_offset, rank_offset;
  off_t path_cache_released, rank_cache_released;
};

class PathPairWriter
{
public:
  PathPairWriter(const PathGraph& graph, size_type file, size_type size_limit,
    size_type& committed_bytes) :
    path_name(graph.path_names.at(file)), rank_name(graph.rank_names.at(file)),
    path_descriptor(-1), rank_descriptor(-1), limit(size_limit),
    predecessor_bytes(committed_bytes), path_count(0), rank_count(0), bytes(0), closed(false)
  {
    this->path_descriptor = ::open(this->path_name.c_str(), O_CREAT | O_TRUNC | O_WRONLY, 0644);
    this->rank_descriptor = ::open(this->rank_name.c_str(), O_CREAT | O_TRUNC | O_WRONLY, 0644);
    if(this->path_descriptor < 0 || this->rank_descriptor < 0)
    {
      throw joinError("cannot create generated path pair", this->path_name);
    }
    adviseSequential(this->path_descriptor); adviseSequential(this->rank_descriptor);
  }

  ~PathPairWriter()
  {
    if(this->path_descriptor >= 0) { ::close(this->path_descriptor); }
    if(this->rank_descriptor >= 0) { ::close(this->rank_descriptor); }
  }

  void writeRecord(const JoinRecord& record)
  {
    size_type record_bytes = sizeof(PathNode) + record.node.ranks() * sizeof(PathNode::rank_type);
    if(record_bytes > this->limit || this->predecessor_bytes > this->limit - record_bytes ||
       this->bytes > this->limit - this->predecessor_bytes - record_bytes)
    {
      throw joinError("configured disk limit exceeded while generating paths");
    }
    if(this->rank_count >= (static_cast<size_type>(1) << 40))
    {
      throw joinError("rank sidecar exceeds the 40-bit PathNode pointer format");
    }
    PathNode node = record.node; node.setPointer(this->rank_count);
    writeAll(this->path_descriptor, &node, sizeof(node), this->path_name);
    writeAll(this->rank_descriptor, record.labels,
      node.ranks() * sizeof(PathNode::rank_type), this->rank_name);
    this->path_count++; this->rank_count += node.ranks(); this->bytes += record_bytes;
    if(static_cast<off_t>(this->path_count * sizeof(PathNode)) -
       this->path_cache_released >= JOIN_CACHE_FLUSH_BYTES)
    {
      trimWrittenCache(this->path_descriptor, this->path_cache_released, false,
        this->path_name);
    }
    if(static_cast<off_t>(this->rank_count * sizeof(PathNode::rank_type)) -
       this->rank_cache_released >= JOIN_CACHE_FLUSH_BYTES)
    {
      trimWrittenCache(this->rank_descriptor, this->rank_cache_released, false,
        this->rank_name);
    }
  }

  void finish(size_type& committed_bytes)
  {
    if(this->closed) { return; }
    trimWrittenCache(this->path_descriptor, this->path_cache_released, true,
      this->path_name);
    trimWrittenCache(this->rank_descriptor, this->rank_cache_released, true,
      this->rank_name);
    if(::close(this->path_descriptor) != 0 || ::close(this->rank_descriptor) != 0)
    {
      throw joinError("cannot close generated path pair", this->path_name);
    }
    this->path_descriptor = -1; this->rank_descriptor = -1;
    committed_bytes += this->bytes; this->closed = true;
  }

  size_type paths() const { return this->path_count; }
  size_type ranks() const { return this->rank_count; }

private:
  PathPairWriter(const PathPairWriter&);
  PathPairWriter& operator=(const PathPairWriter&);

  std::string path_name, rank_name;
  int path_descriptor, rank_descriptor;
  size_type limit, predecessor_bytes, path_count, rank_count, bytes;
  off_t path_cache_released = 0, rank_cache_released = 0;
  bool closed;
};

void
scanJoinSide(const PathGraph& graph, const std::vector<size_type>& shards,
  logical_file_id_t logical, JoinKeyKind kind, ExternalJoinSorter& sorter,
  PathPairWriter* bypass, ExternalPathJoinStats* stats)
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
        if(record.node.sorted())
        {
          if(bypass == nullptr) { throw joinError("missing sorted-path bypass writer"); }
          bypass->writeRecord(record);
          if(stats != nullptr) { stats->sorted_bypass++; }
        }
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
joinSortedRuns(const JoinRun& left_run, const JoinRun& right_run,
  logical_file_id_t logical, size_type byte_budget, PathPairWriter& output,
  bool verify_payloads, ExternalPathJoinStats* stats)
{
  JoinFileReader left(left_run, logical, LEFT_BY_TO);
  JoinFileReader right(right_run, logical, RIGHT_BY_FROM);
  if(verify_payloads) { left.validate(); right.validate(); }
  size_type left_offset = 0, right_offset = 0;
  JoinRecord left_record, right_record;
  while(left_offset < left.size() && right_offset < right.size())
  {
    left.read(left_offset, left_record); right.read(right_offset, right_record);
    if(left_record.key < right_record.key) { left_offset++; continue; }
    if(right_record.key < left_record.key) { right_offset++; continue; }
    node_type key = left_record.key;
    size_type left_limit = left_offset;
    do
    {
      left_limit++;
      if(left_limit >= left.size()) { break; }
      left.read(left_limit, left_record);
    } while(left_record.key == key);
    size_type right_limit = right_offset;
    do
    {
      right_limit++;
      if(right_limit >= right.size()) { break; }
      right.read(right_limit, right_record);
    } while(right_record.key == key);
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
        output.writeRecord(extendRecord(left_record, right_record));
        if(stats != nullptr) { stats->generated_records++; }
      }
    }
    left_offset = left_limit; right_offset = right_limit;
  }
}

} // namespace

size_type
externalPathJoinMinimumBudget()
{
  return JOIN_FIXED_BYTES + 8 * sizeof(JoinRecord);
}

void
externalPathGraphExtend(PathGraph& graph, size_type size_limit,
  const ConstructionParameters& parameters, ExternalPathJoinStats* stats)
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

  PathGraph next(groups.size(), 2 * graph.k(), graph.step() + 1);
  size_type output_file = 0, committed_bytes = 0;
  for(const auto& group : groups)
  {
    logical_file_id_t logical = group.first;
    next.logical_file_ids[output_file] = logical;
    next.physical_shard_ids[output_file] = physical_shard_id_t(output_file);
    PathPairWriter output(next, output_file, size_limit, committed_bytes);

    ExternalJoinSorter right_sorter(logical, RIGHT_BY_FROM, join_budget,
      join_fan_in, parameters.getVerifyWorkspace(), stats);
    scanJoinSide(graph, group.second, logical, RIGHT_BY_FROM, right_sorter, &output, stats);
    JoinRun right = right_sorter.finish();

    ExternalJoinSorter left_sorter(logical, LEFT_BY_TO, join_budget,
      join_fan_in, parameters.getVerifyWorkspace(), stats);
    scanJoinSide(graph, group.second, logical, LEFT_BY_TO, left_sorter, nullptr, stats);
    JoinRun left = left_sorter.finish();

    joinSortedRuns(left, right, logical, join_budget, output,
      parameters.getVerifyWorkspace(), stats);
    TempFile::remove(left.name); TempFile::remove(right.name);
    output.finish(committed_bytes);
    next.path_counts[output_file] = output.paths();
    next.rank_counts[output_file] = output.ranks();
    next.path_count += output.paths(); next.rank_count += output.ranks();
    ExternalPathSortStats sort_stats;
    externalPathGraphSort(next, output_file, sort_budget,
      label_fan_in, (stats == nullptr ? nullptr : &sort_stats));
    if(stats != nullptr)
    {
      stats->label_sort_runs += sort_stats.runs;
      stats->label_merge_passes = std::max(stats->label_merge_passes,
        sort_stats.merge_passes);
      stats->max_records_resident = std::max(stats->max_records_resident,
        sort_stats.max_records_resident);
      stats->max_bytes_resident = std::max(stats->max_bytes_resident,
        sort_stats.max_bytes_resident);
    }
    output_file++;
  }

  graph.clear(); graph.swap(next);
}

} // namespace gcsa
