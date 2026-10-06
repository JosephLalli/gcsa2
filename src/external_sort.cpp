/*
  Copyright (c) 2026 GCSA2 contributors

  Bounded fixed-record external sorting.
*/

#include <gcsa/external_sort.h>
#include <gcsa/internal.h>

#include <algorithm>
#include <array>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <limits>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <unistd.h>
#include <vector>

namespace gcsa
{

namespace
{

// Stream buffers are disabled below, so this covers container bookkeeping,
// heap state, and a small amount of allocator slack that every sort/merge
// phase needs in addition to byte-accounted record storage.
constexpr size_type FIXED_BYTES = 16 * KILOBYTE;
constexpr size_type OUTPUT_BUFFER_BYTES = 64 * KILOBYTE;
// Clean and dirty file-backed pages count against Linux cgroup memory.max.
// Flush completed output prefixes before asking the kernel to evict them, and
// keep only a small sequential tail warm. These are cache bounds, not record
// buffers, and therefore do not change the sort's userspace memory budget.
constexpr off_t CACHE_TAIL_BYTES = 64 * MEGABYTE;
constexpr off_t CACHE_FLUSH_BYTES = 512 * MEGABYTE;
constexpr size_type CACHE_CHECK_RECORDS = 64 * 1024;

void
fail(const std::string& message)
{
  throw std::runtime_error("ExternalFixedRecordSorter: " + message);
}

struct CacheDescriptor
{
  int value;

  CacheDescriptor(const std::string& name, int flags) : value(::open(name.c_str(), flags))
  {
    if(this->value < 0) { fail("cannot open cache window for " + name); }
#if defined(POSIX_FADV_SEQUENTIAL)
    // Advice is non-semantic; unsupported filesystems may safely ignore it.
    static_cast<void>(::posix_fadvise(this->value, 0, 0, POSIX_FADV_SEQUENTIAL));
#endif
  }

  CacheDescriptor(CacheDescriptor&& another) noexcept : value(another.value)
  {
    another.value = -1;
  }

  CacheDescriptor& operator=(CacheDescriptor&& another) noexcept
  {
    if(this != &another)
    {
      if(this->value >= 0) { ::close(this->value); }
      this->value = another.value; another.value = -1;
    }
    return *this;
  }

  ~CacheDescriptor() { if(this->value >= 0) { ::close(this->value); } }

private:
  CacheDescriptor(const CacheDescriptor&);
  CacheDescriptor& operator=(const CacheDescriptor&);
};

void
discardCache(int descriptor, off_t offset, off_t bytes)
{
#if defined(POSIX_FADV_DONTNEED)
  if(bytes > 0)
  {
    static_cast<void>(::posix_fadvise(descriptor, offset, bytes,
      POSIX_FADV_DONTNEED));
  }
#else
  static_cast<void>(descriptor); static_cast<void>(offset);
  static_cast<void>(bytes);
#endif
}

void
trimReadCache(int descriptor, off_t consumed, off_t& released,
  bool complete = false)
{
  if(!complete && consumed - released < CACHE_FLUSH_BYTES) { return; }
  off_t discard_end = (complete ? consumed :
    std::max(released, consumed - CACHE_TAIL_BYTES));
  if(discard_end > released)
  {
    discardCache(descriptor, released, discard_end - released);
    released = discard_end;
  }
}

void
trimWrittenCache(std::ofstream& output, int descriptor, off_t written,
  off_t& released, bool complete, const std::string& name)
{
  if(!complete && written - released < CACHE_FLUSH_BYTES) { return; }
  output.flush();
  if(!output || ::fdatasync(descriptor) != 0)
  {
    fail("cannot sync sequential output " + name);
  }
  off_t discard_end = (complete ? written :
    std::max(released, written - CACHE_TAIL_BYTES));
  if(discard_end > released)
  {
    discardCache(descriptor, released, discard_end - released);
    released = discard_end;
  }
}

off_t
checkedByteOffset(size_type records, size_type record_bytes)
{
  if(record_bytes != 0 &&
     records > static_cast<size_type>(std::numeric_limits<off_t>::max()) / record_bytes)
  {
    fail("file offset exceeds the platform limit");
  }
  return static_cast<off_t>(records * record_bytes);
}

void
updateStats(ExternalFixedRecordSortStats* stats, size_type records,
  size_type bytes)
{
  if(stats == nullptr) { return; }
  stats->max_records_resident = std::max(stats->max_records_resident, records);
  stats->max_bytes_resident = std::max(stats->max_bytes_resident, bytes);
}

size_type
fileRecords(const std::string& name, size_type record_bytes)
{
  std::ifstream input;
  input.rdbuf()->pubsetbuf(nullptr, 0);
  input.open(name.c_str(), std::ios_base::binary);
  if(!input) { fail("cannot open " + name); }
  input.seekg(0, std::ios_base::end);
  std::streamoff bytes = input.tellg();
  if(bytes < 0 || static_cast<std::uint64_t>(bytes) % record_bytes != 0)
  {
    fail("truncated fixed-record file " + name);
  }
  return static_cast<size_type>(static_cast<std::uint64_t>(bytes) / record_bytes);
}

void
writeBytes(std::ofstream& output, const std::uint8_t* data, size_type bytes,
  const std::string& name)
{
  if(bytes == 0) { return; }
  output.write(reinterpret_cast<const char*>(data), bytes);
  if(!output) { fail("cannot write " + name); }
  DiskIO::write_volume += bytes;
}

void
writeEmpty(const std::string& name)
{
  std::ofstream output;
  output.rdbuf()->pubsetbuf(nullptr, 0);
  output.open(name.c_str(), std::ios_base::binary | std::ios_base::trunc);
  if(!output) { fail("cannot create " + name); }
}

size_type
outputBufferBytes(size_type available, size_type record_bytes)
{
  size_type result = std::min(OUTPUT_BUFFER_BYTES, available / 4);
  return std::max(record_bytes, result - result % record_bytes);
}

struct SortPlan
{
  size_type record_bytes, budget, fan_in, run_records, merge_records;
  size_type run_output_bytes;
  bool total_order, ascending_u64;
};

SortPlan
makePlan(size_type record_bytes, size_type byte_budget, size_type requested_fan_in,
  bool total_order = false,
  ExternalFixedRecordSorter::RecordOrder order =
    ExternalFixedRecordSorter::RecordOrder::COMPARATOR)
{
  if(record_bytes == 0) { fail("record width must be nonzero"); }
  constexpr size_type max_in_place_words = 4;
  total_order = total_order && record_bytes % sizeof(std::uint64_t) == 0 &&
    record_bytes / sizeof(std::uint64_t) >= 1 &&
    record_bytes / sizeof(std::uint64_t) <= max_in_place_words;
  bool ascending_u64 = false;
#if defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
  ascending_u64 = total_order &&
    order == ExternalFixedRecordSorter::RecordOrder::ASCENDING_U64 &&
    record_bytes == sizeof(std::uint64_t);
#else
  static_cast<void>(order);
#endif
  if(byte_budget < ExternalFixedRecordSorter::minimumBudget(record_bytes))
  {
    fail("byte budget is too small");
  }

  size_type available = byte_budget - FIXED_BYTES;
  size_type output_bytes = outputBufferBytes(available, record_bytes);
  if(output_bytes > available) { fail("byte budget cannot hold one output record"); }

  // A merge owns one record buffer per input and one output buffer. Restrict
  // fan-in before choosing the per-stream buffer so no phase can exceed the
  // declared working-set budget.
  size_type merge_space = available;
  size_type maximum_fan_in = merge_space / (2 * record_bytes);
  size_type fan_in = std::min(std::max(static_cast<size_type>(2), requested_fan_in),
    maximum_fan_in);
  if(fan_in < 2) { fail("byte budget cannot support a two-way merge"); }
  size_type merge_records = merge_space / ((fan_in + 1) * record_bytes);
  if(merge_records == 0) { fail("byte budget cannot buffer a merge"); }

  // Run formation stores record bytes and, when comparator ties are observable,
  // an offset permutation. The output staging buffer coexists with the sorted
  // run until it has been emitted.
  size_type run_space = available - output_bytes;
  size_type run_records = run_space /
    (total_order ? record_bytes : record_bytes + sizeof(size_type));
  if(run_records == 0) { fail("byte budget cannot hold one sortable record"); }

  return { record_bytes, byte_budget, fan_in, run_records, merge_records,
    output_bytes, total_order, ascending_u64 };
}

struct RunReader
{
  std::ifstream input;
  std::vector<std::uint8_t> buffer;
  size_type record_bytes, capacity, records, offset;
  off_t total_bytes, bytes_read, cache_released;
  bool at_end;
  CacheDescriptor cache;

  RunReader(const std::string& name, size_type width, size_type requested) :
    input(), buffer(), record_bytes(width), capacity(0), records(0), offset(0),
    total_bytes(0), bytes_read(0), cache_released(0), at_end(false),
    cache(name, O_RDONLY)
  {
    input.rdbuf()->pubsetbuf(nullptr, 0);
    input.open(name.c_str(), std::ios_base::binary);
    if(!input) { fail("cannot open run " + name); }
    size_type total = fileRecords(name, width);
    this->total_bytes = checkedByteOffset(total, width);
    capacity = std::max(static_cast<size_type>(1), std::min(requested, total));
    if(total == 0) { at_end = true; return; }
    buffer.resize(capacity * record_bytes);
    refill();
  }

  void refill()
  {
    trimReadCache(this->cache.value, this->bytes_read, this->cache_released,
      this->bytes_read == this->total_bytes);
    input.read(reinterpret_cast<char*>(buffer.data()), buffer.size());
    std::streamsize bytes = input.gcount();
    if(bytes < 0 || static_cast<size_type>(bytes) % record_bytes != 0)
    {
      fail("truncated run");
    }
    DiskIO::read_volume += static_cast<size_type>(bytes);
    this->bytes_read += bytes;
    records = static_cast<size_type>(bytes) / record_bytes;
    offset = 0; at_end = (records == 0);
  }

  const void* current() const
  {
    return buffer.data() + offset * record_bytes;
  }

  void advance()
  {
    offset++;
    if(offset == records) { refill(); }
  }

  RunReader(const RunReader&) = delete;
  RunReader& operator=(const RunReader&) = delete;
  RunReader(RunReader&& another) noexcept :
    input(std::move(another.input)), buffer(std::move(another.buffer)),
    record_bytes(another.record_bytes), capacity(another.capacity),
    records(another.records), offset(another.offset),
    total_bytes(another.total_bytes), bytes_read(another.bytes_read),
    cache_released(another.cache_released), at_end(another.at_end),
    cache(std::move(another.cache))
  {
  }

  RunReader& operator=(RunReader&&) = delete;
};

struct HeapCompare
{
  const std::vector<RunReader>* readers;
  const ExternalFixedRecordSorter::Comparator* compare;

  bool operator()(size_type left, size_type right) const
  {
    int result = (*compare)((*readers)[left].current(), (*readers)[right].current());
    if(result != 0) { return result > 0; }
    // Equal records are intentionally deterministic across the input runs.
    return left > right;
  }
};

template<size_type WORDS>
void
sortAsWords(std::uint8_t* records, size_type count,
  const ExternalFixedRecordSorter::Comparator& compare)
{
  typedef std::array<std::uint64_t, WORDS> word_record;
  word_record* values = reinterpret_cast<word_record*>(records);
  std::sort(values, values + count,
    [&compare](const word_record& left, const word_record& right)
    {
      return compare(left.data(), right.data()) < 0;
    });
}

void
sortRecordsInPlace(std::uint8_t* records, size_type count, const SortPlan& plan,
  const ExternalFixedRecordSorter::Comparator& compare)
{
  if(count < 2) { return; }
  if(plan.ascending_u64)
  {
    std::uint64_t* values = reinterpret_cast<std::uint64_t*>(records);
    std::sort(values, values + count); return;
  }
  switch(plan.record_bytes / sizeof(std::uint64_t))
  {
    case 1: sortAsWords<1>(records, count, compare); return;
    case 2: sortAsWords<2>(records, count, compare); return;
    case 3: sortAsWords<3>(records, count, compare); return;
    case 4: sortAsWords<4>(records, count, compare); return;
    default: fail("in-place sort was planned for an unsupported record width");
  }
}

std::string
makeRunName()
{
  return TempFile::getName("gcsa_fixed_sort_run");
}

std::string
writeSortedRun(const std::uint8_t* records, size_type count,
  const SortPlan& plan, const std::string& name)
{
  std::ofstream output;
  output.rdbuf()->pubsetbuf(nullptr, 0);
  output.open(name.c_str(), std::ios_base::binary | std::ios_base::trunc);
  if(!output) { fail("cannot create run " + name); }
  CacheDescriptor cache(name, O_RDWR);
  off_t written = 0, cache_released = 0;
  size_type remaining = checkedByteOffset(count, plan.record_bytes);
  const std::uint8_t* cursor = records;
  while(remaining > 0)
  {
    size_type chunk = std::min(remaining, plan.run_output_bytes);
    writeBytes(output, cursor, chunk, name);
    written += static_cast<off_t>(chunk);
    cursor += chunk; remaining -= chunk;
    trimWrittenCache(output, cache.value, written, cache_released, false, name);
  }
  trimWrittenCache(output, cache.value, written, cache_released, true, name);
  return name;
}

std::string
writeSortedRun(const std::uint8_t* records, const std::vector<size_type>& order,
  size_type count, const SortPlan& plan, const std::string& name)
{
  std::ofstream output;
  output.rdbuf()->pubsetbuf(nullptr, 0);
  output.open(name.c_str(), std::ios_base::binary | std::ios_base::trunc);
  if(!output) { fail("cannot create run " + name); }
  CacheDescriptor cache(name, O_RDWR);
  off_t written = 0, cache_released = 0;
  std::vector<std::uint8_t> staging(plan.run_output_bytes);
  size_type used = 0;
  for(size_type i = 0; i < count; i++)
  {
    size_type index = order[i];
    if(used + plan.record_bytes > staging.size())
    {
      writeBytes(output, staging.data(), used, name);
      written += static_cast<off_t>(used);
      trimWrittenCache(output, cache.value, written, cache_released, false, name);
      used = 0;
    }
    std::memcpy(staging.data() + used,
      records + index * plan.record_bytes, plan.record_bytes);
    used += plan.record_bytes;
  }
  writeBytes(output, staging.data(), used, name);
  written += static_cast<off_t>(used);
  trimWrittenCache(output, cache.value, written, cache_released, true, name);
  return name;
}

std::string
mergeRuns(const std::vector<std::string>& names, const SortPlan& plan,
  const ExternalFixedRecordSorter::Comparator& compare,
  ExternalFixedRecordSortStats* stats)
{
  if(names.empty()) { fail("cannot merge an empty run set"); }
  if(names.size() == 1) { return names.front(); }
  if(stats != nullptr) { stats->merge_operations++; }

  std::vector<RunReader> readers;
  readers.reserve(names.size());
  size_type buffered_records = 0;
  for(const std::string& name : names)
  {
    readers.emplace_back(name, plan.record_bytes, plan.merge_records);
    buffered_records += readers.back().capacity;
  }
  std::vector<std::uint8_t> output_buffer(plan.merge_records * plan.record_bytes);
  updateStats(stats, buffered_records + plan.merge_records,
    FIXED_BYTES + (buffered_records + plan.merge_records) * plan.record_bytes);

  std::string output_name = makeRunName();
  std::ofstream output;
  output.rdbuf()->pubsetbuf(nullptr, 0);
  output.open(output_name.c_str(), std::ios_base::binary | std::ios_base::trunc);
  if(!output) { fail("cannot create merged run " + output_name); }
  CacheDescriptor output_cache(output_name, O_RDWR);
  off_t output_bytes = 0, output_cache_released = 0;

  std::priority_queue<size_type, std::vector<size_type>, HeapCompare> queue(
    HeapCompare{ &readers, &compare });
  for(size_type i = 0; i < readers.size(); i++)
  {
    if(!readers[i].at_end) { queue.push(i); }
  }
  size_type used = 0;
  while(!queue.empty())
  {
    size_type best = queue.top(); queue.pop();
    std::memcpy(output_buffer.data() + used * plan.record_bytes,
      readers[best].current(), plan.record_bytes);
    used++; readers[best].advance();
    if(!readers[best].at_end) { queue.push(best); }
    if(used == plan.merge_records)
    {
      writeBytes(output, output_buffer.data(), used * plan.record_bytes, output_name);
      output_bytes += checkedByteOffset(used, plan.record_bytes);
      trimWrittenCache(output, output_cache.value, output_bytes,
        output_cache_released, false, output_name);
      used = 0;
    }
  }
  writeBytes(output, output_buffer.data(), used * plan.record_bytes, output_name);
  output_bytes += checkedByteOffset(used, plan.record_bytes);
  trimWrittenCache(output, output_cache.value, output_bytes,
    output_cache_released, true, output_name);
  return output_name;
}

void
removeRuns(std::vector<std::string>& names)
{
  for(std::string& name : names) { TempFile::remove(name); }
  names.clear();
}

void
addRun(std::vector<std::vector<std::string>>& levels, std::string name,
  size_type level, const SortPlan& plan,
  const ExternalFixedRecordSorter::Comparator& compare,
  ExternalFixedRecordSortStats* stats)
{
  if(levels.size() <= level) { levels.resize(level + 1); }
  levels[level].push_back(name);
  if(levels[level].size() < plan.fan_in) { return; }
  std::vector<std::string> inputs; inputs.swap(levels[level]);
  std::string merged = mergeRuns(inputs, plan, compare, stats);
  removeRuns(inputs);
  if(stats != nullptr) { stats->merge_passes = std::max(stats->merge_passes, level + 1); }
  addRun(levels, merged, level + 1, plan, compare, stats);
}

std::string
finalRun(std::vector<std::vector<std::string>>& levels, const SortPlan& plan,
  const ExternalFixedRecordSorter::Comparator& compare,
  ExternalFixedRecordSortStats* stats)
{
  std::vector<std::string> runs;
  for(std::vector<std::string>& level : levels)
  {
    runs.insert(runs.end(), level.begin(), level.end()); level.clear();
  }
  while(runs.size() > plan.fan_in)
  {
    std::vector<std::string> compacted;
    for(size_type first = 0; first < runs.size(); first += plan.fan_in)
    {
      size_type last = std::min(runs.size(), first + plan.fan_in);
      if(last - first == 1)
      {
        compacted.push_back(runs[first]); continue;
      }
      std::vector<std::string> group(runs.begin() + first, runs.begin() + last);
      compacted.push_back(mergeRuns(group, plan, compare, stats));
      removeRuns(group);
    }
    runs.swap(compacted);
    if(stats != nullptr) { stats->merge_passes++; }
  }
  if(runs.empty()) { return std::string(); }
  if(runs.size() == 1) { return runs.front(); }
  std::string merged = mergeRuns(runs, plan, compare, stats);
  removeRuns(runs);
  return merged;
}

std::string
sortToRun(const std::string& input_name, const SortPlan& plan,
  const ExternalFixedRecordSorter::Comparator& compare,
  ExternalFixedRecordSortStats* stats)
{
  size_type total = fileRecords(input_name, plan.record_bytes);
  if(total == 0) { return std::string(); }
  std::ifstream input;
  input.rdbuf()->pubsetbuf(nullptr, 0);
  input.open(input_name.c_str(), std::ios_base::binary);
  if(!input) { fail("cannot open " + input_name); }
  CacheDescriptor input_cache(input_name, O_RDONLY);
  off_t input_bytes = 0, input_cache_released = 0;

  std::vector<std::vector<std::string>> levels;
  for(size_type remaining = total; remaining > 0; )
  {
    // Run formation and merge are separate budget phases. In particular,
    // addRun() can synchronously merge a full level, so release both the raw
    // run and its offset permutation before it is allowed to allocate merge
    // readers/output buffers. Reallocation for the next run is intentional.
    // RETIRED, NOT REMOVED: allocating at the planned capacity.
    //
    // count was computed one line after the allocation that should have used
    // it, so a run sized for the memory limit rather than for the input was
    // fully value-initialised regardless of how few records existed. With
    // --gcsa-sort-run-size 23G the plan plans 1,350,560,518 records while the
    // chr21 redundancy stream holds 563,772,793, so 12.589 GB was memset and
    // never read, on the phase whose page cache then collapsed from 12.67 to
    // 1.00 GiB. Sizing to count changes no run boundary and no output byte:
    // count itself is unchanged.
    //
    // Set SORT_RUNS_AT_PLANNED_CAPACITY to true to restore the old sizing.
    constexpr bool SORT_RUNS_AT_PLANNED_CAPACITY = false;
    size_type count = std::min(remaining, plan.run_records);
    size_type buffer_records = (SORT_RUNS_AT_PLANNED_CAPACITY ?
      plan.run_records : count);
    std::vector<std::uint8_t> records(buffer_records * plan.record_bytes);
    std::vector<size_type> order(plan.total_order ? 0 : buffer_records);
    input.read(reinterpret_cast<char*>(records.data()), count * plan.record_bytes);
    if(input.gcount() != static_cast<std::streamsize>(count * plan.record_bytes))
    {
      fail("short read from " + input_name);
    }
    DiskIO::read_volume += count * plan.record_bytes;
    input_bytes += checkedByteOffset(count, plan.record_bytes);
    trimReadCache(input_cache.value, input_bytes, input_cache_released,
      count == remaining);
    std::string run;
    if(plan.total_order)
    {
      sortRecordsInPlace(records.data(), count, plan, compare);
      updateStats(stats, buffer_records,
        FIXED_BYTES + buffer_records * plan.record_bytes + plan.run_output_bytes);
      run = writeSortedRun(records.data(), count, plan, makeRunName());
    }
    else
    {
      std::iota(order.begin(), order.begin() + count, 0);
      std::sort(order.begin(), order.begin() + count,
        [&records, &plan, &compare](size_type left, size_type right)
        {
          int result = compare(records.data() + left * plan.record_bytes,
            records.data() + right * plan.record_bytes);
          return (result != 0 ? result < 0 : left < right);
        });
      updateStats(stats, buffer_records,
        FIXED_BYTES + buffer_records * (plan.record_bytes + sizeof(size_type)) +
        plan.run_output_bytes);
      run = writeSortedRun(records.data(), order, count, plan, makeRunName());
    }
    if(stats != nullptr) { stats->runs++; }
    std::vector<std::uint8_t>().swap(records);
    std::vector<size_type>().swap(order);
    addRun(levels, run, 0, plan, compare, stats);
    remaining -= count;
  }
  return finalRun(levels, plan, compare, stats);
}

void
installRun(std::string run, const std::string& output_name)
{
  if(std::rename(run.c_str(), output_name.c_str()) != 0)
  {
    // output_name is caller-owned temporary state in every construction route;
    // replacing it here avoids exposing a stale previous result as a success.
    std::remove(output_name.c_str());
    if(std::rename(run.c_str(), output_name.c_str()) != 0)
    {
      TempFile::remove(run); fail("cannot install sorted output " + output_name);
    }
  }
  TempFile::remove(run); // Forget the renamed path from TempFile's registry.
}

void
consumeRun(const std::string& run, const SortPlan& plan,
  const ExternalFixedRecordSorter::Comparator& compare,
  const ExternalFixedRecordSorter::Reducer& reducer,
  const std::string& output_name, ExternalFixedRecordSortStats* stats)
{
  // Reducers often emit one compact result per group. Supply an explicit,
  // budget-accounted output buffer so key/start deduplication does not turn
  // into one syscall per result while still avoiding a hidden stream buffer.
  std::vector<char> output_buffer(plan.merge_records * plan.record_bytes);
  std::ofstream output;
  output.rdbuf()->pubsetbuf(output_buffer.data(), output_buffer.size());
  output.open(output_name.c_str(), std::ios_base::binary | std::ios_base::trunc);
  if(!output) { fail("cannot create reduced output " + output_name); }
  CacheDescriptor output_cache(output_name, O_RDWR);
  off_t output_cache_released = 0;
  RunReader reader(run, plan.record_bytes, plan.merge_records);
  if(reader.at_end)
  {
    trimWrittenCache(output, output_cache.value, 0, output_cache_released,
      true, output_name);
    return;
  }
  std::vector<std::uint8_t> previous(plan.record_bytes);
  std::memcpy(previous.data(), reader.current(), plan.record_bytes);
  reader.advance();
  bool first = true;
  size_type processed = 0;
  updateStats(stats, reader.capacity + 1,
    FIXED_BYTES + (reader.capacity + 1) * plan.record_bytes + output_buffer.size());
  while(!reader.at_end)
  {
    bool last = (compare(previous.data(), reader.current()) != 0);
    reducer(previous.data(), first, last, output);
    first = last;
    std::memcpy(previous.data(), reader.current(), plan.record_bytes);
    reader.advance();
    processed++;
    if(processed % CACHE_CHECK_RECORDS == 0)
    {
      output.flush();
      std::streamoff position = output.tellp();
      if(position < 0) { fail("cannot inspect reduced output " + output_name); }
      trimWrittenCache(output, output_cache.value, static_cast<off_t>(position),
        output_cache_released, false, output_name);
    }
  }
  reducer(previous.data(), first, true, output);
  output.flush();
  std::streamoff position = output.tellp();
  if(position < 0) { fail("cannot inspect reduced output " + output_name); }
  trimWrittenCache(output, output_cache.value, static_cast<off_t>(position),
    output_cache_released, true, output_name);
  DiskIO::write_volume += static_cast<size_type>(position);
  output.close();
  if(!output) { fail("cannot write reduced output " + output_name); }
}

} // namespace

size_type
ExternalFixedRecordSorter::minimumBudget(size_type record_bytes)
{
  if(record_bytes == 0) { return std::numeric_limits<size_type>::max(); }
  if(record_bytes > (std::numeric_limits<size_type>::max() - FIXED_BYTES) / 8)
  {
    return std::numeric_limits<size_type>::max();
  }
  // Enough for a one-record sort buffer, its offset, a one-record output
  // buffer, and the three buffers of a two-way merge.
  return FIXED_BYTES + 8 * (record_bytes + sizeof(size_type));
}

void
ExternalFixedRecordSorter::sort(const std::string& input_name,
  const std::string& output_name, size_type record_bytes, size_type byte_budget,
  size_type requested_fan_in, const Comparator& compare,
  ExternalFixedRecordSortStats* stats, bool total_order, RecordOrder order)
{
  if(!compare) { fail("missing comparator"); }
  if(stats != nullptr) { *stats = ExternalFixedRecordSortStats(); }
  SortPlan plan = makePlan(record_bytes, byte_budget, requested_fan_in,
    total_order, order);
  std::string run = sortToRun(input_name, plan, compare, stats);
  if(run.empty()) { writeEmpty(output_name); return; }
  installRun(run, output_name);
}

void
ExternalFixedRecordSorter::sortAndReduce(const std::string& input_name,
  const std::string& output_name, size_type record_bytes, size_type byte_budget,
  size_type requested_fan_in, const Comparator& compare, const Reducer& reducer,
  ExternalFixedRecordSortStats* stats, bool total_order, RecordOrder order)
{
  if(!compare || !reducer) { fail("missing comparator or reducer"); }
  if(stats != nullptr) { *stats = ExternalFixedRecordSortStats(); }
  SortPlan plan = makePlan(record_bytes, byte_budget, requested_fan_in,
    total_order, order);
  std::string run = sortToRun(input_name, plan, compare, stats);
  if(run.empty()) { writeEmpty(output_name); return; }
  try
  {
    consumeRun(run, plan, compare, reducer, output_name, stats);
  }
  catch(...)
  {
    TempFile::remove(run); throw;
  }
  TempFile::remove(run);
}

} // namespace gcsa
