#include <gcsa/path_graph.h>
#include <gcsa/external_sort.h>
#include <gcsa/path_sort_run.h>

#include <sdsl/wt_algorithm.hpp>

#include <algorithm>
#include <array>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <limits>
#include <memory>
#include <queue>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>

namespace gcsa
{

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr size_type PathLabel::LABEL_LENGTH;
constexpr PathLabel::rank_type PathLabel::NO_RANK;

constexpr size_type PathNode::LABEL_LENGTH;

constexpr size_type PathGraph::UNKNOWN;

constexpr size_type MergedGraph::UNKNOWN;

//------------------------------------------------------------------------------

// Other class variables.

const std::string PathGraph::PREFIX = "gcsa";

const std::string MergedGraph::PREFIX = "gcsa";

//------------------------------------------------------------------------------

PathNode::PathNode(const KMer& kmer, WriteBuffer<PathNode::rank_type>& labels)
{
  this->from = kmer.from; this->to = kmer.to;
  this->fields = 0;

  if(kmer.sorted()) { this->makeSorted(); }
  this->setPredecessors(Key::predecessors(kmer.key));
  this->setOrder(1); this->setLCP(1);

  this->setPointer(labels.size());
  labels.push_back(Key::label(kmer.key));
  labels.push_back(0);  // Dummy value; the last label is not in use.
}

PathNode::PathNode(const PathNode& source,
    const std::vector<PathNode::rank_type>& old_labels, std::vector<PathNode::rank_type>& new_labels)
{
  this->from = source.from; this->to = source.to;
  this->fields = source.fields;

  this->setPointer(new_labels.size());
  for(size_type i = 0, j = source.pointer(); i < source.ranks(); i++, j++)
  {
    new_labels.push_back(old_labels[j]);
  }
}

PathNode::PathNode(const PathNode& left, const PathNode& right,
    const std::vector<PathNode::rank_type>& old_labels, std::vector<PathNode::rank_type>& new_labels)
{
  this->from = left.from; this->to = right.to;
  if(right.sorted()) { this->makeSorted(); }

  this->fields = 0;
  this->setPredecessors(left.predecessors());

  size_type left_order = left.order();
  size_type new_order = left_order + right.order();
  this->setOrder(new_order); this->setLCP(left_order + right.lcp());

  this->setPointer(new_labels.size());
  for(size_type i = 0, j = left.pointer(); i < left_order; i++, j++)
  {
    new_labels.push_back(old_labels[j]);
  }
  for(size_type i = 0, j = right.pointer(); i < right.ranks(); i++, j++)
  {
    new_labels.push_back(old_labels[j]);
  }
}

PathNode::PathNode()
{
  this->from = 0; this->to = 0;
  this->fields = 0;
}

PathNode::PathNode(const PathNode& source)
{
  this->copy(source);
}

PathNode::PathNode(PathNode&& source) noexcept
{
  *this = std::move(source);
}

PathNode::~PathNode()
{
}

PathNode&
PathNode::operator= (const PathNode& source)
{
  if(&source != this)
  {
    this->copy(source);
  }
  return *this;
}

void
PathNode::copy(const PathNode& source)
{
  this->from = source.from; this->to = source.to;
  this->fields = source.fields;
}

PathNode&
PathNode::operator= (PathNode&& source) noexcept
{
  if(&source != this)
  {
    this->from = std::move(source.from);
    this->to = std::move(source.to);
    this->fields = std::move(source.fields);
  }
  return *this;
}

void
PathNode::print(std::ostream& out, const std::vector<PathNode::rank_type>& labels) const
{
  this->print(out, labels.data());
}

void
PathNode::print(std::ostream& out, const rank_type* labels) const
{
  out << "(" << Node::decode(this->from) << " -> " << Node::decode(this->to);
  out << "; o" << this->order();
  out << "; l" << this->lcp();
  for(size_type i = 0; i < this->order(); i++)
  {
    out << (i == 0 ? "; [" : ", ") << this->firstLabel(i, labels);
  }
  for(size_type i = 0; i < this->order(); i++)
  {
    out << (i == 0 ? " to " : ", ") << this->lastLabel(i, labels);
  }
  out << "])";
}

//------------------------------------------------------------------------------

LCP::LCP()
{
}

LCP::LCP(const std::vector<key_type>& keys, size_type _kmer_length)
{
  this->kmer_length = _kmer_length;
  this->total_keys = keys.size();

  sdsl::int_vector<8> buffer(keys.size(), 0);
  for(size_type i = 1; i < keys.size(); i++)
  {
    buffer[i] = Key::lcp(keys[i - 1], keys[i], this->kmer_length);
  }
  directConstruct(this->kmer_lcp, buffer);
}

LCP::LCP(const std::string& key_name, size_type key_count,
  size_type _kmer_length, size_type buffer_bytes)
{
  this->kmer_length = _kmer_length;
  this->total_keys = key_count;
  sdsl::int_vector<8> buffer(key_count, 0);
  if(key_count == 0)
  {
    directConstruct(this->kmer_lcp, buffer); return;
  }
  if(buffer_bytes < sizeof(key_type)) { buffer_bytes = sizeof(key_type); }
  size_type records_per_block = std::max(static_cast<size_type>(1),
    buffer_bytes / sizeof(key_type));
  std::ifstream input;
  input.rdbuf()->pubsetbuf(nullptr, 0);
  input.open(key_name.c_str(), std::ios_base::binary);
  if(!input) { throw std::runtime_error("LCP: cannot open key stream " + key_name); }
  std::vector<key_type> keys(records_per_block);
  key_type previous = 0;
  size_type seen = 0;
  while(seen < key_count)
  {
    size_type count = std::min(records_per_block, key_count - seen);
    input.read(reinterpret_cast<char*>(keys.data()), count * sizeof(key_type));
    if(input.gcount() != static_cast<std::streamsize>(count * sizeof(key_type)))
    {
      throw std::runtime_error("LCP: truncated key stream " + key_name);
    }
    for(size_type i = 0; i < count; i++)
    {
      if(seen + i > 0) { buffer[seen + i] = Key::lcp(previous, keys[i], this->kmer_length); }
      previous = keys[i];
    }
    seen += count;
  }
  key_type extra;
  input.read(reinterpret_cast<char*>(&extra), sizeof(extra));
  if(input.gcount() != 0)
  {
    throw std::runtime_error("LCP: key stream has trailing records " + key_name);
  }
  directConstruct(this->kmer_lcp, buffer);
}

range_type
LCP::min_lcp(const PathNode& a, const PathNode& b, const std::vector<LCP::rank_type>& labels) const
{
  return this->min_lcp(a, b, labels.data(), labels.data());
}

range_type
LCP::max_lcp(const PathNode& a, const PathNode& b, const std::vector<LCP::rank_type>& labels) const
{
  return this->max_lcp(a, b, labels.data(), labels.data());
}

range_type
LCP::min_lcp(const PathNode& a, const PathNode& b,
  const LCP::rank_type* a_labels, const LCP::rank_type* b_labels) const
{
  size_type order = std::min(a.order(), b.order());
  range_type lcp(0, 0);
  while(lcp.first < order && a.firstLabel(lcp.first, a_labels) == b.lastLabel(lcp.first, b_labels))
  {
    lcp.first++;
  }
  if(lcp.first < order)
  {
    size_type left = a.firstLabel(lcp.first, a_labels) + 1;
    size_type right = std::min((size_type)(b.lastLabel(lcp.first, b_labels)), this->total_keys - 1);
    lcp.second = sdsl::quantile_freq(this->kmer_lcp, left, right, 0).first;
  }
  return lcp;
}

range_type
LCP::max_lcp(const PathNode& a, const PathNode& b,
  const LCP::rank_type* a_labels, const LCP::rank_type* b_labels) const
{
  size_type order = std::min(a.order(), b.order());
  range_type lcp(0, 0);
  while(lcp.first < order && a.lastLabel(lcp.first, a_labels) == b.firstLabel(lcp.first, b_labels))
  {
    lcp.first++;
  }
  if(lcp.first < order)
  {
    size_type left = a.lastLabel(lcp.first, a_labels) + 1;
    size_type right = b.firstLabel(lcp.first, b_labels);
    lcp.second = sdsl::quantile_freq(this->kmer_lcp, left, right, 0).first;
  }
  return lcp;
}

void
LCP::swap(LCP& another) noexcept
{
  std::swap(this->kmer_length, another.kmer_length);
  std::swap(this->total_keys, another.total_keys);
  this->kmer_lcp.swap(another.kmer_lcp);
}

//------------------------------------------------------------------------------

/*
  This structure combines a PathNode and its label. It also stores the identifier of
  its source file.
*/

struct PriorityNode
{
  typedef PathLabel::rank_type rank_type;

  constexpr static size_type LABEL_LENGTH = PathLabel::LABEL_LENGTH;
  constexpr static rank_type NO_RANK = PathLabel::NO_RANK;

  rank_type file;
  rank_type label[LABEL_LENGTH + 1];
  PathNode  node;

  inline bool eof() const { return (this->label[0] == NO_RANK); }

  // Compare by first labels.
  inline bool operator< (const PriorityNode& another) const
  {
    size_type order = std::min(this->node.order(), another.node.order());
    for(size_type i = 0; i < order; i++)
    {
      if(this->label[i] != another.label[i]) { return (this->label[i] < another.label[i]); }
    }
    return (this->node.order() < another.node.order());
  }

  inline rank_type firstLabel(size_type i) const { return this->label[i]; }

  inline size_type bytes() const { return this->node.bytes(); }
};

constexpr size_type PriorityNode::LABEL_LENGTH;
constexpr PriorityNode::rank_type PriorityNode::NO_RANK;

//------------------------------------------------------------------------------

/*
  This structure builds a PathGraph. It expects a stream of PriorityNodes. When all paths
  for a certain file have been written, call sort() for that file.
*/

struct PathGraphBuilder
{
  PathGraph graph;
  std::vector<WriteBuffer<PathNode>> path_files;
  std::vector<WriteBuffer<PathNode::rank_type>> rank_files;
  size_type limit;  // Bytes of disk space.

  constexpr static size_type WRITE_BUFFER_SIZE = MEGABYTE;  // PathNodes per thread.

  PathGraphBuilder(size_type file_count, size_type path_order, size_type step, size_type size_limit);
  void close();

  /*
    The file number is assumed to be valid.
    The first call is not thread safe, while the bulk write() is.
  */
  void write(PriorityNode& path);
  void write(std::vector<PathNode>& paths, std::vector<PathNode::rank_type>& labels, size_type memory_limit, size_type file);

  void sort(size_type file, size_type byte_budget, size_type fan_in);
};

constexpr size_type PathGraphBuilder::WRITE_BUFFER_SIZE;

PathGraphBuilder::PathGraphBuilder(size_type file_count, size_type path_order, size_type step, size_type size_limit) :
  graph(file_count, path_order, step),
  path_files(file_count), rank_files(file_count),
  limit(size_limit)
{
  for(size_type file = 0; file < file_count; file++)
  {
    this->path_files[file].open(this->graph.path_names[file]);
    this->rank_files[file].open(this->graph.rank_names[file]);
  }
}

void
PathGraphBuilder::close()
{
  for(size_type file = 0; file < this->path_files.size(); file++)
  {
    this->path_files[file].close();
    this->rank_files[file].close();
  }
}

inline void
writePath(PathNode& path, const PathNode::rank_type* labels,
  WriteBuffer<PathNode>& path_file, WriteBuffer<PathNode::rank_type>& rank_file)
{
  size_type old_ptr = path.pointer();
  size_type limit = old_ptr + path.ranks();

  path.setPointer(rank_file.size());
  path_file.push_back(path);

  for(size_type i = old_ptr; i < limit; i++) { rank_file.push_back(labels[i]); }
}

void
PathGraphBuilder::write(PriorityNode& path)
{
  if(this->graph.bytes() + path.bytes() > this->limit)
  {
    std::cerr << "PathGraphBuilder::write(): Size limit exceeded, construction aborted" << std::endl;
    std::exit(EXIT_SIZE_LIMIT_EXCEEDED);
  }
  writePath(path.node, path.label, this->path_files[path.file], this->rank_files[path.file]);

  this->graph.path_counts[path.file]++; this->graph.path_count++;
  this->graph.rank_counts[path.file] += path.node.ranks();
  this->graph.rank_count += path.node.ranks();
}

void
PathGraphBuilder::write(std::vector<PathNode>& paths, std::vector<PathNode::rank_type>& labels, size_type memory_limit, size_type file)
{
  size_type bytes_required = paths.size() * sizeof(PathNode) + labels.size() * sizeof(PathNode::rank_type);
  #pragma omp critical
  {
    if(bytes_required + this->graph.bytes() > this->limit)
    {
      std::cerr << "PathGraphBuilder::write(): Size limit exceeded, construction aborted" << std::endl;
      std::exit(EXIT_SIZE_LIMIT_EXCEEDED);
    }
    for(size_type i = 0; i < paths.size(); i++)
    {
      writePath(paths[i], labels.data(), this->path_files[file], this->rank_files[file]);
    }
    this->graph.path_counts[file] += paths.size(); this->graph.path_count += paths.size();
    this->graph.rank_counts[file] += labels.size(); this->graph.rank_count += labels.size();
  }
  
  // Check if this file will violate the RAM limit parameter when it is loaded
  size_type future_memory = this->graph.path_counts[file] * sizeof(PathNode) + this->graph.rank_counts[file] * sizeof(PathNode::rank_type);
  if (future_memory > memory_limit) {
    std::cerr << "PathGraphBuilder::write(): Memory use of file " << file << " of kmer paths ("
      << inGigabytes(future_memory) << " GB) exceeds memory limit (" << inGigabytes(memory_limit) << " GB)" << std::endl;
    std::exit(EXIT_SIZE_LIMIT_EXCEEDED);
  }
  
  paths.clear(); labels.clear();
}

void
PathGraphBuilder::sort(size_type file, size_type byte_budget, size_type fan_in)
{
  this->path_files[file].close();
  this->rank_files[file].close();
  externalPathGraphSort(this->graph, file, byte_budget, fan_in);

  if(Verbosity::level >= Verbosity::FULL)
  {
    std::cerr << "PathGraphBuilder::sort(): File " << file << ": Sorted "
              << this->graph.path_counts[file] << " paths" << std::endl;
  }
}

//------------------------------------------------------------------------------

namespace
{

// A run record is deliberately self-contained: the rank pointer in node is
// meaningful only while reading the source/final PathGraph pair. The transient
// disk representation is versioned and prefix-compressed by path_sort_run.cpp.
using PathSortRecord = PathSortRunRecord;

// The source PathGraph uses record-at-a-time decoding because each PathNode
// points to a variable number of ranks. Explicit stream buffers avoid one
// syscall per decoded record while remaining part of the declared budget.
constexpr size_type PATH_SORT_STREAM_BUFFER_BYTES = 64 * KILOBYTE;
constexpr size_type PATH_SORT_RUN_BUFFER_BYTES = 64 * KILOBYTE;
constexpr size_type PATH_SORT_PARALLEL_MIN_RECORDS = 64 * 1024;
// Covers the two source stream buffers, stream state, vector control blocks,
// allocator slack, and the heap.
constexpr size_type PATH_SORT_FIXED_BYTES =
  8192 + 2 * PATH_SORT_STREAM_BUFFER_BYTES;
// Linux charges clean filesystem cache to a cgroup's memory ceiling. The
// sorter therefore retains only a small sequential tail and periodically
// syncs completed output prefixes before making them reclaimable. These
// values intentionally match the external join policy.
constexpr off_t PATH_SORT_CACHE_TAIL_BYTES = 64 * MEGABYTE;
constexpr off_t PATH_SORT_CACHE_FLUSH_BYTES = 512 * MEGABYTE;
constexpr size_type PATH_SORT_SOURCE_CACHE_CHECK_RECORDS = 1024 * 1024;

inline size_type
pathSortReaderReservation()
{
  return PATH_SORT_RUN_BUFFER_BYTES + sizeof(PathSortRunReader) +
    2 * sizeof(size_type);
}

inline size_type
pathSortOutputRecordReservation()
{
  return sizeof(PathNode) +
    (PathLabel::LABEL_LENGTH + 1) * sizeof(PathNode::rank_type);
}

inline void
updatePathSortStats(ExternalPathSortStats* stats, size_type records, size_type bytes)
{
  if(stats == nullptr) { return; }
  stats->max_records_resident = std::max(stats->max_records_resident, records);
  stats->max_bytes_resident = std::max(stats->max_bytes_resident, bytes);
}

inline void
externalSortFailure(const std::string& message)
{
  std::cerr << "externalPathGraphSort(): " << message << std::endl;
  std::exit(EXIT_FAILURE);
}

inline void
recordPathSortRun(ExternalPathSortStats* stats, const PathSortRunWriter& writer)
{
  if(stats == nullptr) { return; }
  if(stats->run_uncompressed_bytes > std::numeric_limits<size_type>::max() - writer.uncompressedBytes() ||
     stats->run_compressed_bytes > std::numeric_limits<size_type>::max() - writer.bytes())
  {
    externalSortFailure("path-sort run byte counter overflow");
  }
  stats->run_uncompressed_bytes += writer.uncompressedBytes();
  stats->run_compressed_bytes += writer.bytes();
}

inline off_t
pathSortByteOffset(size_type records, size_type record_bytes)
{
  if(record_bytes != 0 && records > static_cast<size_type>(std::numeric_limits<off_t>::max()) / record_bytes)
  {
    externalSortFailure("file offset exceeds the platform limit");
  }
  return static_cast<off_t>(records * record_bytes);
}

inline void
addPathSortBytes(off_t& total, size_type records, size_type record_bytes)
{
  off_t bytes = pathSortByteOffset(records, record_bytes);
  if(total > std::numeric_limits<off_t>::max() - bytes)
  {
    externalSortFailure("file offset exceeds the platform limit");
  }
  total += bytes;
}

inline int
openPathSortCacheDescriptor(const std::string& name, int flags)
{
  int descriptor = ::open(name.c_str(), flags);
  if(descriptor < 0) { externalSortFailure("cannot open cache window for " + name); }
#if defined(POSIX_FADV_SEQUENTIAL)
  // This is an optimization only. Durability depends on fdatasync(), and
  // correctness does not depend on whether the filesystem honors the advice.
  static_cast<void>(::posix_fadvise(descriptor, 0, 0, POSIX_FADV_SEQUENTIAL));
#endif
  return descriptor;
}

inline void
discardPathSortCache(int descriptor, off_t offset, off_t bytes)
{
#if defined(POSIX_FADV_DONTNEED)
  if(bytes > 0)
  {
    static_cast<void>(::posix_fadvise(descriptor, offset, bytes, POSIX_FADV_DONTNEED));
  }
#else
  static_cast<void>(descriptor); static_cast<void>(offset); static_cast<void>(bytes);
#endif
}

inline void
trimPathSortReadCache(int descriptor, off_t consumed, off_t& released, bool complete = false)
{
  if(!complete && consumed - released < PATH_SORT_CACHE_FLUSH_BYTES) { return; }
  off_t discard_end = (complete ? consumed : std::max(released, consumed - PATH_SORT_CACHE_TAIL_BYTES));
  if(discard_end > released)
  {
    discardPathSortCache(descriptor, released, discard_end - released);
    released = discard_end;
  }
}

inline void
trimPathSortWrittenCache(std::ofstream& output, int descriptor, off_t written,
  off_t& released, bool complete, const std::string& name)
{
  if(!complete && written - released < PATH_SORT_CACHE_FLUSH_BYTES) { return; }
  output.flush();
  if(!output || ::fdatasync(descriptor) != 0)
  {
    externalSortFailure("cannot sync sequential output " + name);
  }
  off_t discard_end = (complete ? written : std::max(released, written - PATH_SORT_CACHE_TAIL_BYTES));
  if(discard_end > released)
  {
    discardPathSortCache(descriptor, released, discard_end - released);
    released = discard_end;
  }
}

inline bool
pathSortLess(const PathSortRecord& a, const PathSortRecord& b)
{
  size_type order = std::min(a.node.order(), b.node.order());
  for(size_type i = 0; i < order; i++)
  {
    if(a.labels[i] != b.labels[i]) { return (a.labels[i] < b.labels[i]); }
  }
  if(a.node.order() != b.node.order()) { return (a.node.order() < b.node.order()); }
  if(a.node.from != b.node.from) { return (a.node.from < b.node.from); }
  if(a.node.to != b.node.to) { return (a.node.to < b.node.to); }
  if(a.node.predecessors() != b.node.predecessors()) { return (a.node.predecessors() < b.node.predecessors()); }
  if(a.node.lcp() != b.node.lcp()) { return (a.node.lcp() < b.node.lcp()); }
  for(size_type i = 0; i < a.node.ranks(); i++)
  {
    if(a.labels[i] != b.labels[i]) { return (a.labels[i] < b.labels[i]); }
  }
  return false; // Equal records are byte-equivalent after rank pointers are rewritten.
}

inline void
sortPathRecords(std::vector<PathSortRecord>& records, ExternalPathSortStats* stats)
{
  // The balanced parallel quicksort works in place, so all record storage is
  // still covered by the current run reservation. Equivalent records become
  // byte-identical when writeSortedPathPair() replaces their rank pointers.
  if(records.size() >= PATH_SORT_PARALLEL_MIN_RECORDS && omp_get_max_threads() > 1)
  {
    parallelQuickSort(records.begin(), records.end(), pathSortLess);
    if(stats != nullptr) { stats->parallel_sorts++; }
  }
  else
  {
    sequentialSort(records.begin(), records.end(), pathSortLess);
  }
}

inline void
validatePathSortRecord(const PathSortRecord& record, size_type expected_pointer,
  size_type rank_count)
{
  if(record.node.ranks() > PathLabel::LABEL_LENGTH + 1 || record.node.lcp() > record.node.order())
  {
    externalSortFailure("invalid path order or lcp");
  }
  if(expected_pointer > rank_count || record.node.pointer() != expected_pointer ||
     record.node.ranks() > rank_count - expected_pointer)
  {
    externalSortFailure("invalid rank pointer or rank count");
  }
}

inline bool
readPathSortRecord(std::ifstream& paths, std::ifstream& ranks, PathSortRecord& record,
  size_type& rank_offset, size_type rank_count)
{
  paths.read(reinterpret_cast<char*>(&record.node), sizeof(PathNode));
  DiskIO::read_volume += paths.gcount();
  if(paths.eof() && paths.gcount() == 0) { return false; }
  if(paths.gcount() != sizeof(PathNode)) { externalSortFailure("truncated path file"); }
  validatePathSortRecord(record, rank_offset, rank_count);
  size_type count = record.node.ranks();
  ranks.read(reinterpret_cast<char*>(record.labels), count * sizeof(PathNode::rank_type));
  DiskIO::read_volume += ranks.gcount();
  if(ranks.gcount() != (std::streamsize)(count * sizeof(PathNode::rank_type))) { externalSortFailure("truncated rank file"); }
  rank_offset += count;
  return true;
}

inline void
syncPathSortFile(const std::string& name)
{
  int descriptor = ::open(name.c_str(), O_RDONLY);
  if(descriptor < 0 || ::fsync(descriptor) != 0)
  {
    if(descriptor >= 0) { ::close(descriptor); }
    externalSortFailure("cannot sync " + name);
  }
  ::close(descriptor);
}

struct PathSortHeapComparator
{
  const std::vector<std::unique_ptr<PathSortRunReader>>* readers;

  explicit PathSortHeapComparator(
    const std::vector<std::unique_ptr<PathSortRunReader>>* _readers) : readers(_readers) { }

  bool operator() (size_type a, size_type b) const
  {
    const PathSortRecord& left = (*this->readers)[a]->current();
    const PathSortRecord& right = (*this->readers)[b]->current();
    if(pathSortLess(left, right)) { return false; }
    if(pathSortLess(right, left)) { return true; }
    return (a > b);
  }
};

void
writePathSortRun(const std::string& name, const std::vector<PathSortRecord>& records,
  ExternalPathSortStats* stats)
{
  PathSortRunWriter output(name, records.size(), PATH_SORT_RUN_BUFFER_BYTES);
  for(const PathSortRecord& record : records) { output.write(record); }
  output.finish(); recordPathSortRun(stats, output);
}

std::string
mergePathSortRuns(const std::vector<std::string>& inputs, size_type,
  ExternalPathSortStats* stats)
{
  if(stats != nullptr) { stats->merge_operations++; }
  std::string output_name = TempFile::getName("gcsa_path_sort_run");
  std::vector<std::unique_ptr<PathSortRunReader>> readers;
  readers.reserve(inputs.size());
  size_type total_records = 0;
  for(size_type i = 0; i < inputs.size(); i++)
  {
    readers.emplace_back(new PathSortRunReader(inputs[i], PATH_SORT_RUN_BUFFER_BYTES));
    if(total_records > std::numeric_limits<size_type>::max() - readers.back()->records())
    {
      externalSortFailure("merged path-sort record count overflow");
    }
    total_records += readers.back()->records();
  }
  PathSortRunWriter output(output_name, total_records, PATH_SORT_RUN_BUFFER_BYTES);
  updatePathSortStats(stats, readers.size() + 1,
    PATH_SORT_FIXED_BYTES + PATH_SORT_RUN_BUFFER_BYTES + readers.size() *
      (PATH_SORT_RUN_BUFFER_BYTES + sizeof(PathSortRunReader) + 2 * sizeof(size_type)));
  std::priority_queue<size_type, std::vector<size_type>, PathSortHeapComparator>
    queue{PathSortHeapComparator(&readers)};
  for(size_type i = 0; i < readers.size(); i++)
  {
    if(!(readers[i]->atEnd())) { queue.push(i); }
  }
  while(!queue.empty())
  {
    size_type best = queue.top(); queue.pop();
    output.write(readers[best]->current()); readers[best]->advance();
    if(!(readers[best]->atEnd())) { queue.push(best); }
  }
  output.finish(); recordPathSortRun(stats, output);
  return output_name;
}

void
writeSortedPathPair(const std::string& run_name, const std::string& path_name,
  const std::string& rank_name, size_type expected_paths, size_type expected_ranks,
  size_type buffer_records, ExternalPathSortStats* stats)
{
  PathSortRunReader reader(run_name, PATH_SORT_RUN_BUFFER_BYTES);
  std::ofstream paths, ranks;
  paths.rdbuf()->pubsetbuf(nullptr, 0); ranks.rdbuf()->pubsetbuf(nullptr, 0);
  paths.open(path_name.c_str(), std::ios_base::binary);
  ranks.open(rank_name.c_str(), std::ios_base::binary);
  if(!paths || !ranks) { externalSortFailure("cannot create sorted path pair"); }
  int path_cache_descriptor = openPathSortCacheDescriptor(path_name, O_RDWR);
  int rank_cache_descriptor = openPathSortCacheDescriptor(rank_name, O_RDWR);
  off_t path_bytes = 0, rank_bytes = 0;
  off_t path_cache_released = 0, rank_cache_released = 0;
  size_type path_buffer_records = std::max((size_type)1,
    std::min(buffer_records, expected_paths));
  size_type maximum_rank_records =
    (buffer_records > std::numeric_limits<size_type>::max() / (PathLabel::LABEL_LENGTH + 1)
      ? std::numeric_limits<size_type>::max()
      : buffer_records * (PathLabel::LABEL_LENGTH + 1));
  size_type rank_buffer_records = std::max((size_type)1,
    std::min(maximum_rank_records, expected_ranks));
  std::vector<PathNode> path_buffer; path_buffer.reserve(path_buffer_records);
  std::vector<PathNode::rank_type> rank_buffer;
  rank_buffer.reserve(rank_buffer_records);
  updatePathSortStats(stats, 1 + path_buffer_records,
    PATH_SORT_FIXED_BYTES + reader.bufferBytes() + sizeof(PathSortRunReader) +
    path_buffer_records * sizeof(PathNode) +
    rank_buffer_records * sizeof(PathNode::rank_type));
  size_type path_count = 0, rank_count = 0;
  while(!(reader.atEnd()))
  {
    PathNode node = reader.current().node;
    if(rank_buffer.size() + node.ranks() > rank_buffer.capacity())
    {
      DiskIO::write(ranks, rank_buffer.data(), rank_buffer.size());
      addPathSortBytes(rank_bytes, rank_buffer.size(), sizeof(PathNode::rank_type));
      trimPathSortWrittenCache(ranks, rank_cache_descriptor, rank_bytes,
        rank_cache_released, false, rank_name);
      rank_buffer.clear();
    }
    node.setPointer(rank_count);
    path_buffer.push_back(node);
    for(size_type i = 0; i < node.ranks(); i++) { rank_buffer.push_back(reader.current().labels[i]); }
    path_count++; rank_count += node.ranks(); reader.advance();
    if(rank_buffer.size() >= rank_buffer.capacity())
    {
      DiskIO::write(ranks, rank_buffer.data(), rank_buffer.size());
      addPathSortBytes(rank_bytes, rank_buffer.size(), sizeof(PathNode::rank_type));
      trimPathSortWrittenCache(ranks, rank_cache_descriptor, rank_bytes,
        rank_cache_released, false, rank_name);
      rank_buffer.clear();
    }
    if(path_buffer.size() == path_buffer.capacity())
    {
      DiskIO::write(paths, path_buffer.data(), path_buffer.size());
      addPathSortBytes(path_bytes, path_buffer.size(), sizeof(PathNode));
      trimPathSortWrittenCache(paths, path_cache_descriptor, path_bytes,
        path_cache_released, false, path_name);
      path_buffer.clear();
    }
  }
  if(!path_buffer.empty())
  {
    DiskIO::write(paths, path_buffer.data(), path_buffer.size());
    addPathSortBytes(path_bytes, path_buffer.size(), sizeof(PathNode));
  }
  if(!rank_buffer.empty())
  {
    DiskIO::write(ranks, rank_buffer.data(), rank_buffer.size());
    addPathSortBytes(rank_bytes, rank_buffer.size(), sizeof(PathNode::rank_type));
  }
  trimPathSortWrittenCache(paths, path_cache_descriptor, path_bytes,
    path_cache_released, true, path_name);
  trimPathSortWrittenCache(ranks, rank_cache_descriptor, rank_bytes,
    rank_cache_released, true, rank_name);
  paths.close(); ranks.close();
  ::close(path_cache_descriptor); ::close(rank_cache_descriptor);
  if(path_count != expected_paths || rank_count != expected_ranks) { externalSortFailure("sorted path pair has incorrect counts"); }
}

void
addPathSortRun(std::vector<std::vector<std::string>>& levels, const std::string& run,
  size_type level, size_type fan_in, size_type buffer_records,
  ExternalPathSortStats* stats)
{
  if(levels.size() <= level) { levels.resize(level + 1); }
  levels[level].push_back(run);
  if(levels[level].size() < fan_in) { return; }

  std::vector<std::string> inputs; inputs.swap(levels[level]);
  std::string merged = mergePathSortRuns(inputs, buffer_records, stats);
  for(size_type i = 0; i < inputs.size(); i++) { TempFile::remove(inputs[i]); }
  if(stats != nullptr) { stats->merge_passes = std::max(stats->merge_passes, level + 1); }
  addPathSortRun(levels, merged, level + 1, fan_in, buffer_records, stats);
}

} // namespace

struct ExternalPathSortSink::Impl
{
  Impl(PathGraph& target, size_type target_file, size_type byte_budget,
    size_type requested_fan_in, size_type size_limit,
    size_type& already_committed, ExternalPathSortStats* statistics) :
    graph(target), file(target_file), limit(size_limit),
    committed_bytes(already_committed), stats(statistics), fan_in(0),
    merge_records(0), run_records(0), path_count(0), rank_count(0),
    payload_bytes(0), complete(false)
  {
    if(this->file >= this->graph.files()) { externalSortFailure("invalid sink file number"); }
    if(this->graph.path_counts[this->file] != 0 || this->graph.rank_counts[this->file] != 0)
    {
      externalSortFailure("streaming sort sink target is not empty");
    }
    if(byte_budget < externalPathGraphSortMinimumBudget())
    {
      externalSortFailure("streaming sort sink byte budget is too small");
    }
    if(this->stats != nullptr) { *(this->stats) = ExternalPathSortStats(); }

    size_type record_bytes = sizeof(PathSortRecord);
    size_type available = byte_budget - PATH_SORT_FIXED_BYTES;
    size_type reader_bytes = pathSortReaderReservation();
    size_type merge_available = available - PATH_SORT_RUN_BUFFER_BYTES;
    size_type maximum_fan_in = merge_available / reader_bytes;
    this->fan_in = std::min(std::max(static_cast<size_type>(2), requested_fan_in),
      maximum_fan_in);
    if(this->fan_in < 2)
    {
      externalSortFailure("streaming sort sink cannot support a two-way merge");
    }
    size_type materialize_available = available - reader_bytes;
    this->merge_records = materialize_available /
      pathSortOutputRecordReservation();
    if(this->merge_records == 0)
    {
      externalSortFailure("streaming sort sink cannot buffer a merge");
    }

    // The producer and this sink overlap only during the final join scan. The
    // quarter-budget reserve covers its two bounded join readers, allocator
    // metadata, OpenMP worker stacks, and the rolling dirty-cache tail.
    size_type run_available = available - PATH_SORT_RUN_BUFFER_BYTES;
    size_type run_bytes = run_available - run_available / 4;
    this->run_records = std::max(static_cast<size_type>(1),
      run_bytes / record_bytes);
    this->records.reserve(this->run_records);
    updatePathSortStats(this->stats, this->run_records,
      PATH_SORT_FIXED_BYTES + PATH_SORT_RUN_BUFFER_BYTES +
      this->run_records * record_bytes);
  }

  ~Impl()
  {
    for(size_type level = 0; level < this->levels.size(); level++)
    {
      for(std::string& run : this->levels[level]) { TempFile::remove(run); }
    }
  }

  void write(const PathNode& source, const PathNode::rank_type* labels)
  {
    if(this->complete) { externalSortFailure("cannot write to a finished streaming sort sink"); }
    if(source.ranks() > PathLabel::LABEL_LENGTH + 1 || source.lcp() > source.order())
    {
      externalSortFailure("invalid path supplied to streaming sort sink");
    }
    if(source.ranks() > 0 && labels == nullptr)
    {
      externalSortFailure("missing labels for streaming sort sink record");
    }
    size_type record_bytes = sizeof(PathNode) +
      source.ranks() * sizeof(PathNode::rank_type);
    if(record_bytes > this->limit ||
       this->committed_bytes > this->limit - record_bytes ||
       this->payload_bytes > this->limit - this->committed_bytes - record_bytes)
    {
      externalSortFailure("configured disk limit exceeded while generating sorted paths");
    }
    if(this->rank_count >= (static_cast<size_type>(1) << 40))
    {
      externalSortFailure("rank sidecar exceeds the 40-bit PathNode pointer format");
    }

    PathSortRecord record;
    record.node = source;
    // Run records are self-contained; the pointer is rewritten only when the
    // final rank sidecar is emitted.
    record.node.setPointer(0);
    size_type i = 0;
    for(; i < source.ranks(); i++) { record.labels[i] = labels[i]; }
    for(; i < PathLabel::LABEL_LENGTH + 1; i++)
    {
      record.labels[i] = PathLabel::NO_RANK;
    }
    this->records.push_back(record);
    this->path_count++;
    this->rank_count += source.ranks();
    this->payload_bytes += record_bytes;
    if(this->records.size() == this->run_records) { this->flush(); }
  }

  void flush()
  {
    if(this->records.empty()) { return; }
    sortPathRecords(this->records, this->stats);
    std::string name = TempFile::getName("gcsa_path_sort_run");
    writePathSortRun(name, this->records, this->stats);

    // A leveled merge may start immediately in addPathSortRun(). Release the
    // run buffer first so sort and merge reservations never coexist.
    std::vector<PathSortRecord>().swap(this->records);
    addPathSortRun(this->levels, name, 0, this->fan_in,
      this->merge_records, this->stats);
    if(this->stats != nullptr) { this->stats->runs++; }
    this->records.reserve(this->run_records);
  }

  std::string finishRuns()
  {
    this->flush();
    std::vector<PathSortRecord>().swap(this->records);
    std::vector<std::string> runs;
    for(size_type level = 0; level < this->levels.size(); level++)
    {
      runs.insert(runs.end(), this->levels[level].begin(), this->levels[level].end());
      this->levels[level].clear();
    }
    while(runs.size() > this->fan_in)
    {
      std::vector<std::string> compacted;
      for(size_type first = 0; first < runs.size(); first += this->fan_in)
      {
        size_type last = std::min(runs.size(), first + this->fan_in);
        std::vector<std::string> group(runs.begin() + first, runs.begin() + last);
        std::string merged = mergePathSortRuns(group, this->merge_records, this->stats);
        for(std::string& run : group) { TempFile::remove(run); }
        compacted.push_back(merged);
      }
      runs.swap(compacted);
      if(this->stats != nullptr) { this->stats->merge_passes++; }
    }
    if(runs.empty())
    {
      std::string name = TempFile::getName("gcsa_path_sort_run");
      writePathSortRun(name, this->records, this->stats);
      runs.push_back(name);
    }
    if(runs.size() == 1) { return runs.front(); }
    std::string merged = mergePathSortRuns(runs, this->merge_records, this->stats);
    for(std::string& run : runs) { TempFile::remove(run); }
    return merged;
  }

  void finish()
  {
    if(this->complete) { return; }
    std::string merged = this->finishRuns();
    const std::string& final_path = this->graph.path_names[this->file];
    const std::string& final_rank = this->graph.rank_names[this->file];
    std::string partial_path = final_path + ".partial";
    std::string partial_rank = final_rank + ".partial";
    std::remove(partial_path.c_str()); std::remove(partial_rank.c_str());
    writeSortedPathPair(merged, partial_path, partial_rank,
      this->path_count, this->rank_count, this->merge_records, this->stats);
    TempFile::remove(merged);
    syncPathSortFile(partial_path); syncPathSortFile(partial_rank);
    if(std::rename(partial_path.c_str(), final_path.c_str()) != 0 ||
       std::rename(partial_rank.c_str(), final_rank.c_str()) != 0)
    {
      externalSortFailure("cannot atomically install streaming sorted path pair");
    }

    this->graph.path_counts[this->file] = this->path_count;
    this->graph.rank_counts[this->file] = this->rank_count;
    this->graph.path_count += this->path_count;
    this->graph.rank_count += this->rank_count;
    this->committed_bytes += this->payload_bytes;
    this->complete = true;
  }

  PathGraph& graph;
  size_type file, limit;
  size_type& committed_bytes;
  ExternalPathSortStats* stats;
  size_type fan_in, merge_records, run_records;
  size_type path_count, rank_count, payload_bytes;
  bool complete;
  std::vector<PathSortRecord> records;
  std::vector<std::vector<std::string>> levels;
};

ExternalPathSortSink::ExternalPathSortSink(PathGraph& graph, size_type file,
  size_type byte_budget, size_type fan_in, size_type size_limit,
  size_type& committed_bytes, ExternalPathSortStats* stats) :
  impl(new Impl(graph, file, byte_budget, fan_in, size_limit,
    committed_bytes, stats))
{
}

ExternalPathSortSink::~ExternalPathSortSink()
{
  delete this->impl;
}

void
ExternalPathSortSink::write(const PathNode& node, const PathNode::rank_type* labels)
{
  this->impl->write(node, labels);
}

void
ExternalPathSortSink::finish()
{
  this->impl->finish();
}

size_type
ExternalPathSortSink::paths() const
{
  return this->impl->path_count;
}

size_type
ExternalPathSortSink::ranks() const
{
  return this->impl->rank_count;
}

size_type
ExternalPathSortSink::bytes() const
{
  return this->impl->payload_bytes;
}

void
externalPathGraphSort(PathGraph& graph, size_type file, size_type byte_budget, size_type fan_in,
  ExternalPathSortStats* stats)
{
  if(file >= graph.files()) { externalSortFailure("invalid file number"); }
  if(byte_budget < externalPathGraphSortMinimumBudget()) { externalSortFailure("byte budget is too small"); }
  if(stats != nullptr) { *stats = ExternalPathSortStats(); }
  size_type record_bytes = sizeof(PathSortRecord);
  size_type available = byte_budget - PATH_SORT_FIXED_BYTES;
  size_type reader_bytes = pathSortReaderReservation();
  size_type merge_available = available - PATH_SORT_RUN_BUFFER_BYTES;
  size_type max_fan_in = merge_available / reader_bytes;
  fan_in = std::min(std::max((size_type)2, fan_in), max_fan_in);
  if(fan_in < 2) { externalSortFailure("byte budget cannot support a two-way merge"); }
  size_type materialize_available = available - reader_bytes;
  size_type merge_records = materialize_available /
    pathSortOutputRecordReservation();
  if(merge_records == 0) { externalSortFailure("byte budget cannot buffer a merge"); }
  // Use three quarters for the in-place record sort. The remaining quarter is
  // deliberately much larger than the byte-counted stream buffers, OpenMP
  // task stacks, allocator metadata, and rolling output-cache window. Keeping
  // that safety reserve while crossing common one-run thresholds avoids a
  // complete read/write merge pass. Compute the slack first to avoid overflow
  // on a large byte limit.
  size_type run_available = available - PATH_SORT_RUN_BUFFER_BYTES;
  size_type run_bytes = run_available - run_available / 4;
  size_type run_records = std::max((size_type)1, run_bytes / record_bytes);

  std::array<char, PATH_SORT_STREAM_BUFFER_BYTES> path_stream_buffer;
  std::array<char, PATH_SORT_STREAM_BUFFER_BYTES> rank_stream_buffer;
  std::ifstream paths, ranks;
  paths.rdbuf()->pubsetbuf(path_stream_buffer.data(), path_stream_buffer.size());
  ranks.rdbuf()->pubsetbuf(rank_stream_buffer.data(), rank_stream_buffer.size());
  graph.open(paths, ranks, file);
  int source_path_cache_descriptor = openPathSortCacheDescriptor(graph.path_names[file], O_RDONLY);
  int source_rank_cache_descriptor = openPathSortCacheDescriptor(graph.rank_names[file], O_RDONLY);
  off_t source_path_cache_released = 0, source_rank_cache_released = 0;
  if(graph.path_counts[file] > std::numeric_limits<size_type>::max() / sizeof(PathNode) ||
     graph.rank_counts[file] > std::numeric_limits<size_type>::max() / sizeof(PathNode::rank_type) ||
     fileSize(paths) != graph.path_counts[file] * sizeof(PathNode) ||
     fileSize(ranks) != graph.rank_counts[file] * sizeof(PathNode::rank_type))
  {
    externalSortFailure("source path/rank file length does not match metadata");
  }
  // A budget is an upper bound, not an allocation target. Reserving the full
  // budget for a tiny run caused small builds to touch many gigabytes.
  run_records = std::max((size_type)1,
    std::min(run_records, graph.path_counts[file]));
  // At most fan_in - 1 runs remain at each level. Run metadata therefore uses
  // O(fan_in * log(number_of_runs)) memory rather than one string per run.
  std::vector<std::vector<std::string>> levels;
  std::vector<PathSortRecord> records; records.reserve(run_records);
  updatePathSortStats(stats, run_records, PATH_SORT_FIXED_BYTES +
    PATH_SORT_RUN_BUFFER_BYTES + run_records * record_bytes);
  size_type rank_offset = 0, path_count = 0;
  PathSortRecord record;
  while(readPathSortRecord(paths, ranks, record, rank_offset, graph.rank_counts[file]))
  {
    records.push_back(record); path_count++;
    if(path_count % PATH_SORT_SOURCE_CACHE_CHECK_RECORDS == 0)
    {
      trimPathSortReadCache(source_path_cache_descriptor,
        pathSortByteOffset(path_count, sizeof(PathNode)), source_path_cache_released);
      trimPathSortReadCache(source_rank_cache_descriptor,
        pathSortByteOffset(rank_offset, sizeof(PathNode::rank_type)), source_rank_cache_released);
    }
    if(records.size() == run_records)
    {
      sortPathRecords(records, stats);
      std::string name = TempFile::getName("gcsa_path_sort_run");
      writePathSortRun(name, records, stats); records.clear();
      addPathSortRun(levels, name, 0, fan_in, merge_records, stats);
      if(stats != nullptr) { stats->runs++; }
    }
  }
  if(rank_offset != graph.rank_counts[file] || path_count != graph.path_counts[file])
  {
    externalSortFailure("path/rank sidecar counts do not match metadata");
  }
  trimPathSortReadCache(source_path_cache_descriptor,
    pathSortByteOffset(path_count, sizeof(PathNode)), source_path_cache_released, true);
  trimPathSortReadCache(source_rank_cache_descriptor,
    pathSortByteOffset(rank_offset, sizeof(PathNode::rank_type)), source_rank_cache_released, true);
  ::close(source_path_cache_descriptor); ::close(source_rank_cache_descriptor);
  if(!records.empty())
  {
    sortPathRecords(records, stats);
    std::string name = TempFile::getName("gcsa_path_sort_run");
    writePathSortRun(name, records, stats);
    addPathSortRun(levels, name, 0, fan_in, merge_records, stats);
    if(stats != nullptr) { stats->runs++; }
  }
  paths.close(); ranks.close();
  // Run construction and run merging are sequential phases. Drop the run
  // buffer capacity before allocating per-reader merge buffers so the byte
  // budget remains a ceiling instead of becoming a per-phase allowance.
  std::vector<PathSortRecord>().swap(records);

  std::vector<std::string> runs;
  for(size_type level = 0; level < levels.size(); level++)
  {
    runs.insert(runs.end(), levels[level].begin(), levels[level].end());
  }
  while(runs.size() > fan_in)
  {
    std::vector<std::string> compacted;
    for(size_type first = 0; first < runs.size(); first += fan_in)
    {
      size_type last = std::min(runs.size(), first + fan_in);
      std::vector<std::string> group(runs.begin() + first, runs.begin() + last);
      std::string merged = mergePathSortRuns(group, merge_records, stats);
      for(size_type i = 0; i < group.size(); i++) { TempFile::remove(group[i]); }
      compacted.push_back(merged);
    }
    runs.swap(compacted);
    if(stats != nullptr) { stats->merge_passes++; }
  }
  if(runs.empty())
  {
    std::string name = TempFile::getName("gcsa_path_sort_run");
    writePathSortRun(name, records, stats); runs.push_back(name);
  }
  // A run is already in final label order. Sending a single run through the
  // k-way merger only copies the complete payload once before the normal
  // path/rank materialization pass, which is especially costly for the common
  // production case where a generous run budget creates exactly one run.
  std::string merged;
  if(runs.size() == 1)
  {
    merged = runs.front();
  }
  else
  {
    merged = mergePathSortRuns(runs, merge_records, stats);
    for(size_type i = 0; i < runs.size(); i++) { TempFile::remove(runs[i]); }
  }

  std::string final_path = TempFile::getName(PathGraph::PREFIX), final_rank = TempFile::getName(PathGraph::PREFIX);
  std::string partial_path = final_path + ".partial", partial_rank = final_rank + ".partial";
  writeSortedPathPair(merged, partial_path, partial_rank, graph.path_counts[file], graph.rank_counts[file], merge_records, stats);
  TempFile::remove(merged);
  syncPathSortFile(partial_path); syncPathSortFile(partial_rank);
  if(std::rename(partial_path.c_str(), final_path.c_str()) != 0 || std::rename(partial_rank.c_str(), final_rank.c_str()) != 0)
  {
    externalSortFailure("cannot atomically install sorted path pair");
  }
  std::string old_path = graph.path_names[file], old_rank = graph.rank_names[file];
  graph.path_names[file] = final_path; graph.rank_names[file] = final_rank;
  if(graph.delete_files)
  {
    TempFile::remove(old_path); TempFile::remove(old_rank);
  }
}

size_type
externalPathGraphSortMinimumBudget()
{
  // The largest minimum phase is a two-way merge: two independently decoded
  // input byte streams plus one prefix-compressed output stream. All three are
  // fixed-size; no decoded input vector scales with the run length.
  return PATH_SORT_FIXED_BYTES + PATH_SORT_RUN_BUFFER_BYTES +
    2 * pathSortReaderReservation();
}

//------------------------------------------------------------------------------

struct PathGraphMerger;

/*
  PathGraphMerger normally only needs a sliding window. An equal-label range
  can, however, be arbitrarily long before prune() decides how to emit it.
  Keep that range addressable on disk once its resident budget is exhausted.
*/
template<class Element>
struct SpillableGroup
{
  std::vector<Element> memory;
  std::string          filename;
  int                  file;
  size_type            elements, disk_elements, offset, byte_limit, write_bytes, read_bytes, advised_write;
  mutable size_type    cache_offset;
  mutable std::vector<Element> read_cache;
  size_type*           spill_counter;

  SpillableGroup(size_type limit, size_type* counter = nullptr) : file(-1), elements(0), disk_elements(0), offset(0),
    byte_limit(std::max(2 * static_cast<size_type>(sizeof(Element)), limit)),
    write_bytes(this->byte_limit / 2), read_bytes(this->byte_limit - this->write_bytes), advised_write(0), cache_offset(0),
    spill_counter(counter)
  {
    this->memory.reserve(this->byte_limit / sizeof(Element));
  }
  ~SpillableGroup() { this->clear(); }

  inline bool spilled() const { return this->file >= 0; }
  inline size_type size() const { return this->elements - this->offset; }
  inline bool buffered(size_type i) const
  {
    return (i >= this->offset && i < this->elements);
  }

  void push_back(const Element& value)
  {
    if(!this->spilled() && (this->memory.size() + 1) * sizeof(Element) > this->byte_limit)
    {
      this->filename = TempFile::getName("gcsa_prune_group");
      this->file = ::open(this->filename.c_str(), O_CREAT | O_TRUNC | O_RDWR, 0600);
      if(this->file < 0) { throw std::runtime_error("PathGraph::prune(): cannot create spill file"); }
#if defined(POSIX_FADV_SEQUENTIAL)
      static_cast<void>(::posix_fadvise(this->file, 0, 0, POSIX_FADV_SEQUENTIAL));
#endif
      if(this->spill_counter != nullptr) { (*this->spill_counter)++; }
      this->flush();
    }
    this->memory.push_back(value);
    this->elements++;
    if(this->spilled() && this->memory.size() * sizeof(Element) >= this->write_bytes) { this->flush(); }
  }

  Element get(size_type i) const
  {
    if(i < this->offset || i >= this->elements) { throw std::out_of_range("PathGraph::prune(): spill index"); }
    if(!this->spilled()) { return this->memory[i - this->offset]; }
    if(i >= this->disk_elements) { return this->memory[i - this->disk_elements]; }
    if(i < this->cache_offset || i >= this->cache_offset + this->read_cache.size())
    {
      this->discard(this->cache_offset, this->read_cache.size());
      this->cache_offset = i;
      size_type count = std::min(this->read_bytes / sizeof(Element), this->disk_elements - i);
      this->read_cache.resize(count); this->read(this->read_cache.data(), count, i);
    }
    return this->read_cache[i - this->cache_offset];
  }

  void set(size_type i, const Element& value)
  {
    if(i < this->offset || i >= this->elements) { throw std::out_of_range("PathGraph::prune(): spill index"); }
    if(!this->spilled()) { this->memory[i - this->offset] = value; return; }
    if(i >= this->disk_elements) { this->memory[i - this->disk_elements] = value; return; }
    this->flush(); this->write(&value, 1, i);
    if(::fdatasync(this->file) != 0) { throw std::runtime_error("PathGraph::prune(): spill sync failed"); }
    this->discard(i, 1); sdsl::util::clear(this->read_cache);
  }

  void seek(size_type i)
  {
    if(i > this->elements) { throw std::out_of_range("PathGraph::prune(): spill seek"); }
    if(!this->spilled())
    {
      this->memory.erase(this->memory.begin(), this->memory.begin() + (i - this->offset));
    }
    this->offset = i;
  }

  void clear()
  {
    if(this->file >= 0) { ::close(this->file); this->file = -1; }
    if(!this->filename.empty()) { TempFile::remove(this->filename); this->filename.clear(); }
    sdsl::util::clear(this->memory); sdsl::util::clear(this->read_cache);
    this->elements = 0; this->disk_elements = 0; this->offset = 0; this->advised_write = 0; this->cache_offset = 0;
  }

private:
  void flush()
  {
    if(this->memory.empty()) { return; }
    this->write(this->memory.data(), this->memory.size(), this->disk_elements);
    this->disk_elements += this->memory.size(); sdsl::util::clear(this->memory);
    this->memory.reserve(std::max(static_cast<size_type>(1), this->write_bytes / sizeof(Element)));
    if(::fdatasync(this->file) != 0) { throw std::runtime_error("PathGraph::prune(): spill sync failed"); }
    this->discard(this->advised_write, this->disk_elements - this->advised_write);
    this->advised_write = this->disk_elements;
  }

  static size_type bytesFor(size_type count)
  {
    if(count > std::numeric_limits<size_type>::max() / sizeof(Element))
    {
      throw std::overflow_error("PathGraph::prune(): spill byte count overflow");
    }
    return count * sizeof(Element);
  }

  static off_t fileOffset(size_type start)
  {
    size_type offset = bytesFor(start);
    if(offset > static_cast<size_type>(std::numeric_limits<off_t>::max()))
    {
      throw std::overflow_error("PathGraph::prune(): spill offset overflow");
    }
    return static_cast<off_t>(offset);
  }

  void discard(size_type start, size_type count) const
  {
#if defined(POSIX_FADV_DONTNEED)
    if(count > 0 && this->file >= 0)
    {
      size_type bytes = bytesFor(count); off_t offset = fileOffset(start);
      static_cast<void>(::posix_fadvise(this->file, offset,
        static_cast<off_t>(bytes), POSIX_FADV_DONTNEED));
    }
#else
    static_cast<void>(start); static_cast<void>(count);
#endif
  }

  void write(const Element* values, size_type count, size_type start) const
  {
    if(values == nullptr && count != 0) { throw std::invalid_argument("PathGraph::prune(): null spill write"); }
    const char* data = reinterpret_cast<const char*>(values);
    size_type bytes = bytesFor(count), done = 0; off_t offset = fileOffset(start);
    while(done < bytes)
    {
      ssize_t result = ::pwrite(this->file, data + done, bytes - done, offset + static_cast<off_t>(done));
      if(result < 0 && errno == EINTR) { continue; }
      if(result <= 0) { throw std::runtime_error("PathGraph::prune(): spill write failed"); }
      DiskIO::write_volume += static_cast<size_type>(result);
      done += static_cast<size_type>(result);
    }
  }

  void read(Element* values, size_type count, size_type start) const
  {
    if(values == nullptr && count != 0) { throw std::invalid_argument("PathGraph::prune(): null spill read"); }
    char* data = reinterpret_cast<char*>(values);
    size_type bytes = bytesFor(count), done = 0; off_t offset = fileOffset(start);
    while(done < bytes)
    {
      ssize_t result = ::pread(this->file, data + done, bytes - done, offset + static_cast<off_t>(done));
      if(result < 0 && errno == EINTR) { continue; }
      if(result <= 0) { throw std::runtime_error("PathGraph::prune(): spill read failed"); }
      DiskIO::read_volume += static_cast<size_type>(result);
      done += static_cast<size_type>(result);
    }
  }
};

struct PathRange
{
  size_type  from, to;
  range_type left_lcp, range_lcp, right_lcp;

  inline range_type range() const { return range_type(this->from, this->to); }
  inline size_type length() const { return this->to + 1 - this->from; }

  PathRange() : from(0), to(0), left_lcp(0, 0), range_lcp(0, 0), right_lcp(0, 0) { }
  PathRange(size_type start, size_type stop, range_type _left_lcp, PathGraphMerger& merger);
};

struct PathRangeRecord
{
  size_type from, to;
  size_type left_first, left_second;
  size_type range_first, range_second;
  size_type right_first, right_second;
};

template<class Element> struct SpillableDequeCodec;

template<>
struct SpillableDequeCodec<PathRange>
{
  typedef PathRangeRecord record_type;

  static record_type encode(const PathRange& value)
  {
    return { value.from, value.to,
      value.left_lcp.first, value.left_lcp.second,
      value.range_lcp.first, value.range_lcp.second,
      value.right_lcp.first, value.right_lcp.second };
  }

  static PathRange decode(const record_type& value)
  {
    PathRange result;
    result.from = value.from; result.to = value.to;
    result.left_lcp = range_type(value.left_first, value.left_second);
    result.range_lcp = range_type(value.range_first, value.range_second);
    result.right_lcp = range_type(value.right_first, value.right_second);
    return result;
  }
};

/*
  An append-only spill file with deque-relative access. The only mutation the
  merger needs is replacing the current front after a batched prefix removal;
  keep that replacement as a one-record overlay instead of dirtying a random
  disk page. All accessors return copies, so spilling or popping never leaves
  an invalid reference in the caller.
*/
template<class Element>
struct SpillableDeque
{
  typedef SpillableDequeCodec<Element> codec_type;
  typedef typename codec_type::record_type record_type;

  std::vector<record_type> memory;
  std::string          filename;
  int                  file;
  size_type            begin_index, end_index, disk_end;
  size_type            resident_records, write_records, read_records;
  mutable size_type    cache_offset;
  mutable std::vector<record_type> read_cache;
  bool                 replaced_front;
  Element              front_value;
  size_type*           spill_counter;

  SpillableDeque(size_type byte_limit, size_type* counter = nullptr) :
    file(-1), begin_index(0), end_index(0), disk_end(0),
    resident_records(std::max(static_cast<size_type>(2), byte_limit / sizeof(record_type))),
    write_records(std::max(static_cast<size_type>(1), this->resident_records / 2)),
    read_records(std::max(static_cast<size_type>(1), this->resident_records - this->write_records)),
    cache_offset(0), replaced_front(false), front_value(), spill_counter(counter)
  {
    this->memory.reserve(this->resident_records);
  }

  ~SpillableDeque() { this->release(); }

  inline bool spilled() const { return (this->file >= 0); }
  inline bool empty() const { return (this->size() == 0); }
  inline size_type size() const { return (this->end_index - this->begin_index); }

  Element front() const
  {
    if(this->empty()) { throw std::out_of_range("PathGraphMerger: empty range deque"); }
    return (this->replaced_front ? this->front_value : this->at(0));
  }

  Element back() const
  {
    if(this->empty()) { throw std::out_of_range("PathGraphMerger: empty range deque"); }
    return (this->size() == 1 && this->replaced_front ? this->front_value : this->at(this->size() - 1));
  }

  Element operator[](size_type i) const
  {
    if(i == 0 && this->replaced_front)
    {
      if(this->empty()) { throw std::out_of_range("PathGraphMerger: empty range deque"); }
      return this->front_value;
    }
    return this->at(i);
  }

  void replaceFront(const Element& value)
  {
    if(this->empty()) { throw std::out_of_range("PathGraphMerger: empty range deque"); }
    this->front_value = value; this->replaced_front = true;
  }

  void pushBack(const Element& value)
  {
    if(!this->spilled() && this->memory.capacity() < this->resident_records)
    {
      this->memory.reserve(this->resident_records);
    }
    if(!this->spilled() && this->memory.size() == this->resident_records)
    {
      this->openSpill(); this->flush();
    }
    if(this->spilled() && this->memory.size() == this->write_records) { this->flush(); }
    this->memory.push_back(codec_type::encode(value)); this->end_index++;
  }

  void popFront(size_type count)
  {
    if(count > this->size()) { throw std::out_of_range("PathGraphMerger: range deque pop"); }
    if(count == 0) { return; }
    this->replaced_front = false;

    if(!this->spilled())
    {
      this->memory.erase(this->memory.begin(), this->memory.begin() + count);
      this->begin_index = 0; this->end_index = this->memory.size();
      return;
    }

    size_type old_begin = this->begin_index;
    this->begin_index += count;
    size_type disk_stop = std::min(this->begin_index, this->disk_end);
    if(disk_stop > old_begin) { this->discard(old_begin, disk_stop - old_begin); }
    sdsl::util::clear(this->read_cache); this->cache_offset = 0;

    if(this->begin_index > this->disk_end)
    {
      size_type consumed = this->begin_index - this->disk_end;
      this->memory.erase(this->memory.begin(), this->memory.begin() + consumed);
      this->disk_end = this->begin_index;
    }
    if(this->empty()) { this->resetEmpty(); }
  }

  void clear()
  {
    this->release();
  }

private:
  Element at(size_type i) const
  {
    if(i >= this->size()) { throw std::out_of_range("PathGraphMerger: range deque index"); }
    size_type absolute = this->begin_index + i;
    if(!this->spilled()) { return codec_type::decode(this->memory[absolute]); }
    if(absolute >= this->disk_end) { return codec_type::decode(this->memory[absolute - this->disk_end]); }
    if(absolute < this->cache_offset || absolute >= this->cache_offset + this->read_cache.size())
    {
      this->discard(this->cache_offset, this->read_cache.size());
      this->cache_offset = absolute;
      size_type count = std::min(this->read_records, this->disk_end - absolute);
      this->read_cache.resize(count); this->read(this->read_cache.data(), count, absolute);
    }
    return codec_type::decode(this->read_cache[absolute - this->cache_offset]);
  }

  void openSpill()
  {
    this->filename = TempFile::getName("gcsa_prune_ranges");
    this->file = ::open(this->filename.c_str(), O_CREAT | O_TRUNC | O_RDWR, 0600);
    if(this->file < 0) { throw std::runtime_error("PathGraphMerger: cannot create range spill file"); }
#if defined(POSIX_FADV_SEQUENTIAL)
    static_cast<void>(::posix_fadvise(this->file, 0, 0, POSIX_FADV_SEQUENTIAL));
#endif
    if(this->spill_counter != nullptr) { (*this->spill_counter)++; }
  }

  void flush()
  {
    if(this->memory.empty()) { return; }
    this->write(this->memory.data(), this->memory.size(), this->disk_end);
    size_type start = this->disk_end; this->disk_end += this->memory.size();
    sdsl::util::clear(this->memory); this->memory.reserve(this->write_records);
    if(::fdatasync(this->file) != 0) { throw std::runtime_error("PathGraphMerger: range spill sync failed"); }
    this->discard(start, this->disk_end - start);
  }

  void resetEmpty()
  {
    sdsl::util::clear(this->memory); this->memory.reserve(this->write_records);
    sdsl::util::clear(this->read_cache);
    if(::ftruncate(this->file, 0) != 0) { throw std::runtime_error("PathGraphMerger: range spill truncate failed"); }
    this->begin_index = 0; this->end_index = 0; this->disk_end = 0;
    this->cache_offset = 0;
  }

  void release()
  {
    if(this->file >= 0) { ::close(this->file); this->file = -1; }
    if(!this->filename.empty()) { TempFile::remove(this->filename); this->filename.clear(); }
    sdsl::util::clear(this->memory); sdsl::util::clear(this->read_cache);
    this->begin_index = 0; this->end_index = 0; this->disk_end = 0;
    this->cache_offset = 0; this->replaced_front = false;
  }

  static size_type bytesFor(size_type count)
  {
    if(count > std::numeric_limits<size_type>::max() / sizeof(record_type))
    {
      throw std::overflow_error("PathGraphMerger: range spill byte count overflow");
    }
    return count * sizeof(record_type);
  }

  static off_t fileOffset(size_type start)
  {
    size_type offset = bytesFor(start);
    if(offset > static_cast<size_type>(std::numeric_limits<off_t>::max()))
    {
      throw std::overflow_error("PathGraphMerger: range spill offset overflow");
    }
    return static_cast<off_t>(offset);
  }

  void discard(size_type start, size_type count) const
  {
#if defined(POSIX_FADV_DONTNEED)
    if(count > 0 && this->file >= 0)
    {
      size_type bytes = bytesFor(count);
      static_cast<void>(::posix_fadvise(this->file, fileOffset(start),
        static_cast<off_t>(bytes), POSIX_FADV_DONTNEED));
    }
#else
    static_cast<void>(start); static_cast<void>(count);
#endif
  }

  void write(const record_type* values, size_type count, size_type start) const
  {
    const char* data = reinterpret_cast<const char*>(values);
    size_type bytes = bytesFor(count), done = 0; off_t offset = fileOffset(start);
    while(done < bytes)
    {
      ssize_t result = ::pwrite(this->file, data + done, bytes - done,
        offset + static_cast<off_t>(done));
      if(result < 0 && errno == EINTR) { continue; }
      if(result <= 0) { throw std::runtime_error("PathGraphMerger: range spill write failed"); }
      DiskIO::write_volume += static_cast<size_type>(result);
      done += static_cast<size_type>(result);
    }
  }

  void read(record_type* values, size_type count, size_type start) const
  {
    char* data = reinterpret_cast<char*>(values);
    size_type bytes = bytesFor(count), done = 0; off_t offset = fileOffset(start);
    while(done < bytes)
    {
      ssize_t result = ::pread(this->file, data + done, bytes - done,
        offset + static_cast<off_t>(done));
      if(result < 0 && errno == EINTR) { continue; }
      if(result <= 0) { throw std::runtime_error("PathGraphMerger: range spill read failed"); }
      DiskIO::read_volume += static_cast<size_type>(result);
      done += static_cast<size_type>(result);
    }
  }
};

/*
  Keep a fixed number of path/rank descriptor pairs and read one current head
  synchronously. This removes the former two 1 MiB buffers, two descriptors,
  and two reader threads per physical shard. The heap still owns one compact
  PriorityNode per shard because exact k-way merge ordering requires a head.
*/
struct PathGraphInputCache
{
  struct Entry
  {
    size_type file, stamp;
    int path, rank;
    off_t path_high, rank_high, path_released, rank_released;

    Entry() : file(PathGraph::UNKNOWN), stamp(0), path(-1), rank(-1),
      path_high(0), rank_high(0), path_released(0), rank_released(0) { }
  };

  const PathGraph& graph;
  std::vector<Entry> entries;
  size_type clock;
  PathGraphMergeStats* stats;

  constexpr static size_type MAX_PAIRS = 32;
  constexpr static off_t CACHE_TAIL = 64 * KILOBYTE;

  PathGraphInputCache(const PathGraph& source, PathGraphMergeStats* merge_stats) :
    graph(source), entries(), clock(0), stats(merge_stats)
  {
    this->entries.reserve(MAX_PAIRS);
  }

  ~PathGraphInputCache() { this->close(); }

  void read(size_type file, size_type offset, PathNode& node,
    PathNode::rank_type* labels)
  {
    Entry& entry = this->get(file);
    off_t path_offset = this->checkedOffset(offset, sizeof(PathNode));
    this->preadAll(entry.path, &node, sizeof(node), path_offset);
    entry.path_high = std::max(entry.path_high,
      path_offset + static_cast<off_t>(sizeof(node)));
    this->trim(entry.path, entry.path_high, entry.path_released);

    if(node.ranks() > PathLabel::LABEL_LENGTH + 1)
    {
      throw std::runtime_error("PathGraphMerger: invalid path rank count");
    }
    size_type rank_bytes = node.ranks() * sizeof(PathNode::rank_type);
    off_t rank_offset = this->checkedOffset(node.pointer(), sizeof(PathNode::rank_type));
    this->preadAll(entry.rank, labels, rank_bytes, rank_offset);
    entry.rank_high = std::max(entry.rank_high,
      rank_offset + static_cast<off_t>(rank_bytes));
    this->trim(entry.rank, entry.rank_high, entry.rank_released);
  }

  void close()
  {
    for(Entry& entry : this->entries) { this->close(entry); }
    this->entries.clear();
  }

private:
  Entry& get(size_type file)
  {
    this->clock++;
    for(Entry& entry : this->entries)
    {
      if(entry.file == file) { entry.stamp = this->clock; return entry; }
    }

    Entry* target = nullptr;
    if(this->entries.size() < MAX_PAIRS)
    {
      this->entries.push_back(Entry()); target = &(this->entries.back());
      if(this->stats != nullptr)
      {
        this->stats->max_open_input_pairs = std::max(
          this->stats->max_open_input_pairs, static_cast<size_type>(this->entries.size()));
      }
    }
    else
    {
      target = &(*std::min_element(this->entries.begin(), this->entries.end(),
        [](const Entry& left, const Entry& right) { return left.stamp < right.stamp; }));
      this->close(*target); *target = Entry();
    }

    target->file = file; target->stamp = this->clock;
    target->path = ::open(this->graph.path_names[file].c_str(), O_RDONLY);
    if(target->path < 0) { throw std::runtime_error("PathGraphMerger: cannot open path input"); }
    target->rank = ::open(this->graph.rank_names[file].c_str(), O_RDONLY);
    if(target->rank < 0)
    {
      ::close(target->path); target->path = -1;
      throw std::runtime_error("PathGraphMerger: cannot open rank input");
    }
#if defined(POSIX_FADV_SEQUENTIAL)
    static_cast<void>(::posix_fadvise(target->path, 0, 0, POSIX_FADV_SEQUENTIAL));
    static_cast<void>(::posix_fadvise(target->rank, 0, 0, POSIX_FADV_SEQUENTIAL));
#endif
    return *target;
  }

  void close(Entry& entry)
  {
    if(entry.path >= 0)
    {
      this->discard(entry.path, entry.path_released, entry.path_high);
      ::close(entry.path); entry.path = -1;
    }
    if(entry.rank >= 0)
    {
      this->discard(entry.rank, entry.rank_released, entry.rank_high);
      ::close(entry.rank); entry.rank = -1;
    }
  }

  static off_t checkedOffset(size_type records, size_type width)
  {
    if(width != 0 && records > static_cast<size_type>(std::numeric_limits<off_t>::max()) / width)
    {
      throw std::overflow_error("PathGraphMerger: input offset overflow");
    }
    return static_cast<off_t>(records * width);
  }

  static void preadAll(int descriptor, void* target, size_type bytes, off_t offset)
  {
    char* data = reinterpret_cast<char*>(target); size_type done = 0;
    while(done < bytes)
    {
      ssize_t result = ::pread(descriptor, data + done, bytes - done,
        offset + static_cast<off_t>(done));
      if(result < 0 && errno == EINTR) { continue; }
      if(result <= 0) { throw std::runtime_error("PathGraphMerger: truncated input"); }
      DiskIO::read_volume += static_cast<size_type>(result);
      done += static_cast<size_type>(result);
    }
  }

  static void trim(int descriptor, off_t high, off_t& released)
  {
    if(high <= released + 2 * CACHE_TAIL) { return; }
    off_t keep_from = high - CACHE_TAIL;
    discard(descriptor, released, keep_from); released = keep_from;
  }

  static void discard(int descriptor, off_t from, off_t to)
  {
#if defined(POSIX_FADV_DONTNEED)
    if(to > from)
    {
      static_cast<void>(::posix_fadvise(descriptor, from, to - from,
        POSIX_FADV_DONTNEED));
    }
#else
    static_cast<void>(descriptor); static_cast<void>(from); static_cast<void>(to);
#endif
  }
};

constexpr size_type PathGraphInputCache::MAX_PAIRS;
constexpr off_t PathGraphInputCache::CACHE_TAIL;

/*
  This structure reads a buffered stream of PriorityNodes in sorted order and outputs a
  stream of ranges of PriorityNodes with the same label.
*/

struct PathGraphMerger
{
  const PathGraph&                              graph;
  const LCP&                                    lcp;

  // Buffers.
  SpillableDeque<PathRange>                     ranges;
  SpillableGroup<PriorityNode>                  buffer;

  // Input cache and priority queue.
  PathGraphInputCache                           input_files;
  std::vector<size_type>                        offsets;
  PriorityQueue<PriorityNode>                   inputs;

  PathGraphMerger(const PathGraph& path_graph, const LCP& kmer_lcp,
    size_type group_buffer_bytes = MEGABYTE,
    PathGraphMergeStats* stats = nullptr);
  void close();

  inline size_type size() const { return this->graph.size(); }

  /*
    Iterates through ranges of paths with the same label.
  */
  range_type first();
  range_type next();
  inline bool atEnd(range_type range) const { return (range.first >= this->size()); }

  inline range_type range_lcp(size_type i, size_type j) const // i < j
  {
    PriorityNode left = this->buffer.get(i), right = this->buffer.get(j);
    return this->lcp.min_lcp(left.node, right.node, left.label, right.label);
  }

  inline range_type border_lcp(size_type i, size_type j) const // i < j
  {
    PriorityNode left = this->buffer.get(i), right = this->buffer.get(j);
    return this->lcp.max_lcp(left.node, right.node, left.label, right.label);
  }

  /*
    Extends the front range forward into a maximal range of paths such that

    (a) comp(next_range) returns true for each additional next_range; and
    (b) the paths share a common prefix no other path has.
  */
  template<class FromComparator>
  range_type extendRange(FromComparator& comp);

  /*
    Merges the path nodes into paths[range.second] of the front range.
  */
  void mergePathNodes();

  // Find the rightmost path with the same label.
  size_type rangeEnd(size_type start);

  // Add the next PriorityNode to buffer.
  void bufferNext();

  // Read the next PriorityNode from the file.
  void read(PriorityNode& path);
};

PathGraphMerger::PathGraphMerger(const PathGraph& path_graph, const LCP& kmer_lcp,
  size_type group_buffer_bytes, PathGraphMergeStats* stats) :
  graph(path_graph), lcp(kmer_lcp),
  ranges(group_buffer_bytes, (stats == nullptr ? nullptr : &(stats->range_spills))),
  buffer(group_buffer_bytes, (stats == nullptr ? nullptr : &(stats->priority_spills))),
  input_files(path_graph, stats),
  offsets(path_graph.files()), inputs(path_graph.files())
{
  if(stats != nullptr) { *stats = PathGraphMergeStats(); }
  for(size_type file = 0; file < path_graph.files(); file++)
  {
    this->offsets[file] = 0;
    this->inputs[file].file = file; this->read(this->inputs[file]);
  }
  this->inputs.heapify();
}

void
PathGraphMerger::close()
{
  this->ranges.clear();
  this->buffer.clear();
  this->input_files.close();
  this->offsets.clear();
  this->inputs.clear();
}

range_type
PathGraphMerger::first()
{
  if(this->buffer.offset > 0)
  {
    std::cerr << "PathGraphMerger::first(): Cannot seek back in the buffer" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  this->ranges.clear();

  this->ranges.pushBack(PathRange(0, this->rangeEnd(0), range_type(0, 0), *this));
  return this->ranges.front().range();
}

range_type
PathGraphMerger::next()
{
  PathRange temp = this->ranges.front(); this->ranges.popFront(1);
  if(this->ranges.empty())
  {
    this->ranges.pushBack(
      PathRange(temp.to + 1, this->rangeEnd(temp.to + 1), temp.right_lcp, *this));
  }
  this->buffer.seek(this->ranges.front().from);
  return this->ranges.front().range();
}

template<class FromComparator>
range_type
PathGraphMerger::extendRange(FromComparator& comp)
{
  PathRange range = this->ranges.front();
  size_type curr = 1;
  range_type parent_lcp = range.range_lcp;

  while(true)
  {
    // Find the next range.
    if(curr >= this->ranges.size())
    {
      PathRange temp = this->ranges.back();
      this->ranges.pushBack(PathRange(temp.to + 1, this->rangeEnd(temp.to + 1), temp.right_lcp, *this));
    }
    PathRange next_range = this->ranges[curr];
    if(next_range.from >= this->size()) { break; }

    // Is this a suffix tree node with the same start nodes as range and lcp > range.left_lcp?
    if(!(comp(next_range.range()))) { break; }
    parent_lcp = std::min(parent_lcp, next_range.left_lcp);
    if(parent_lcp <= range.left_lcp) { break; }
    if(next_range.right_lcp >= parent_lcp) { curr++; continue; }

    // Replace range with [range.from, next_range.to].
    range.to = next_range.to;
    range.range_lcp = parent_lcp; range.right_lcp = next_range.right_lcp;
    this->ranges.popFront(curr);
    this->ranges.replaceFront(range); curr = 1;
  }

  return range.range();
}

void
PathGraphMerger::mergePathNodes()
{
  PathRange range = this->ranges.front();
  if(range.length() == 1)
  {
    PriorityNode node = this->buffer.get(range.to);
    node.node.makeSorted(); this->buffer.set(range.to, node);
    return;
  }

  PriorityNode first = this->buffer.get(range.from);
  PriorityNode merged = this->buffer.get(range.to);

  size_type order = range.range_lcp.first;
  if(range.range_lcp.second > 0)
  {
    range_type ranks(first.node.firstLabel(order, first.label),
      merged.node.lastLabel(order, merged.label));
    merged.label[order] = ranks.first;
    merged.label[order + 1] = ranks.second;
    order++;
  }

  merged.node.makeSorted();
  merged.node.setOrder(order);
  merged.node.setLCP(range.range_lcp.first);
  for(size_type i = range.from; i < range.to; i++)
  {
    merged.node.addPredecessors(this->buffer.get(i).node);
  }
  this->buffer.set(range.to, merged);
}

size_type
PathGraphMerger::rangeEnd(size_type start)
{
  if(!(this->buffer.buffered(start))) { this->bufferNext(); }

  size_type stop = start;
  while(stop + 1 < this->size())
  {
    if(!(this->buffer.buffered(stop + 1))) { this->bufferNext(); }
    if(this->buffer.get(start) < this->buffer.get(stop + 1)) { break; }
    stop++;
  }
  return stop;
}

void
PathGraphMerger::bufferNext()
{
  this->buffer.push_back(this->inputs[0]);  // Add to the buffer.
  this->read(this->inputs[0]);  // Read the next.
  this->inputs.down(0); // Restore heap order.
}

void
PathGraphMerger::read(PriorityNode& path)
{
  if(this->offsets[path.file] >= this->graph.path_counts[path.file])
  {
    path.node.setOrder(1);
    path.node.setLCP(1);
    path.label[0] = PriorityNode::NO_RANK;
  }
  else
  {
    this->input_files.read(path.file, this->offsets[path.file], path.node, path.label);
    path.node.setPointer(0);  // Label is now stored in the PriorityNode.
    this->offsets[path.file]++;
  }
}

PathRange::PathRange(size_type start, size_type stop, range_type _left_lcp, PathGraphMerger& merger) :
  from(start), to(stop),
  left_lcp(_left_lcp),
  range_lcp(0, 0),
  right_lcp(0, 0)
{
  if(stop < merger.size()) { this->range_lcp = merger.range_lcp(start, stop); }
  if(stop + 1 < merger.size()) { this->right_lcp = merger.border_lcp(stop, stop + 1); }
}

//------------------------------------------------------------------------------

PathGraph::PathGraph(const InputGraph& source, sdsl::int_vector<0>& distinct_labels)
{
  this->path_count = 0; this->rank_count = 0; this->range_count = 0;
  this->order = source.k(); this->doubling_steps = 0;
  this->unique = UNKNOWN; this->redundant = UNKNOWN;
  this->unsorted = UNKNOWN; this->nondeterministic = UNKNOWN;
  this->delete_files = true;

  if(source.files() > std::numeric_limits<uint32_t>::max())
  {
    std::cerr << "PathGraph::PathGraph(): Too many logical input files" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  for(size_type file = 0; file < source.files(); file++)
  {
    this->logical_file_ids.push_back(logical_file_id_t(static_cast<uint32_t>(file)));
    this->physical_shard_ids.push_back(physical_shard_id_t(file));
    std::string path_name = TempFile::getName(PREFIX);
    this->path_names.push_back(path_name);
    std::string rank_name = TempFile::getName(PREFIX);
    this->rank_names.push_back(rank_name);

    // Read KMers, sort them, and convert the keys labels to the ranks of those labels.
    std::vector<KMer> kmers;
    source.read(kmers, file);
    parallelQuickSort(kmers.begin(), kmers.end());
    size_type current_rank = 0;
    for(size_type i = 0; i < kmers.size(); i++)
    {
      while(Key::label(kmers[i].key) > distinct_labels[current_rank]) { current_rank++; }
      kmers[i].key = Key::replace(kmers[i].key, current_rank);
    }

    // Convert the KMers to PathNodes.
    WriteBuffer<PathNode> path_buffer(path_name);
    WriteBuffer<PathNode::rank_type> rank_buffer(rank_name);
    for(size_type i = 0; i < kmers.size(); i++)
    {
      path_buffer.push_back(PathNode(kmers[i], rank_buffer));
    }
    this->path_counts.push_back(path_buffer.size()); this->path_count += path_buffer.size();
    this->rank_counts.push_back(rank_buffer.size()); this->rank_count += rank_buffer.size();
    path_buffer.close(); rank_buffer.close();
  }

  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    std::cerr << "PathGraph::PathGraph(): " << this->size() << " paths with "
              << this->ranks() << " ranks" << std::endl;
    std::cerr << "PathGraph::PathGraph(): " << inGigabytes(this->bytes()) << " GB in "
              << this->files() << " file(s)" << std::endl;
  }
}

PathGraph::PathGraph(size_type file_count, size_type path_order, size_type steps) :
  path_names(file_count), rank_names(file_count), path_counts(file_count, 0), rank_counts(file_count, 0),
  logical_file_ids(), physical_shard_ids(),
  path_count(0), rank_count(0), range_count(0), order(path_order), doubling_steps(steps),
  unique(0), redundant(0), unsorted(0), nondeterministic(0),
  delete_files(true)
{
  if(file_count > std::numeric_limits<uint32_t>::max())
  {
    std::cerr << "PathGraph::PathGraph(): Too many logical input files" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  this->logical_file_ids.reserve(file_count);
  this->physical_shard_ids.reserve(file_count);
  for(size_type file = 0; file < this->files(); file++)
  {
    this->logical_file_ids.push_back(logical_file_id_t(static_cast<uint32_t>(file)));
    this->physical_shard_ids.push_back(physical_shard_id_t(file));
    this->path_names[file] = TempFile::getName(PREFIX);
    this->rank_names[file] = TempFile::getName(PREFIX);
  }
}

PathGraph::PathGraph(const std::string& path_name, const std::string& rank_name)
{
  this->path_count = 0; this->rank_count = 0; this->range_count = 0;
  this->order = 0; this->doubling_steps = 0;
  this->unique = 0; this->redundant = 0;
  this->unsorted = 0; this->nondeterministic = 0;
  this->delete_files = false;

  this->path_names.push_back(path_name);
  this->logical_file_ids.push_back(logical_file_id_t(0));
  this->physical_shard_ids.push_back(physical_shard_id_t(0));
  std::ifstream path_file(path_name, std::ios_base::binary);
  if(!path_file)
  {
    std::cerr << "PathGraph::PathGraph(): Cannot open path file " << path_name << std::endl;
    std::exit(EXIT_FAILURE);
  }
  this->path_count = fileSize(path_file) / sizeof(PathNode);
  this->path_counts.push_back(this->path_count);
  path_file.close();

  this->rank_names.push_back(rank_name);
  std::ifstream rank_file(rank_name, std::ios_base::binary);
  if(!rank_file)
  {
    std::cerr << "PathGraph::PathGraph(): Cannot open rank file " << rank_name << std::endl;
    std::exit(EXIT_FAILURE);
  }
  this->rank_count = fileSize(rank_file) / sizeof(PathNode::rank_type);
  this->rank_counts.push_back(this->rank_count);
  rank_file.close();
}

PathGraph::~PathGraph()
{
  this->clear();
}

void
PathGraph::clear()
{
  if(this->delete_files)
  {
    for(size_type file = 0; file < this->files(); file++)
    {
      TempFile::remove(this->path_names[file]);
      TempFile::remove(this->rank_names[file]);
    }
  }
  this->path_names.clear();
  this->rank_names.clear();
  this->path_counts.clear();
  this->rank_counts.clear();
  this->logical_file_ids.clear();
  this->physical_shard_ids.clear();

  this->path_count = 0; this->rank_count = 0;
  this->order = 0;
  this->unique = UNKNOWN; this->redundant = UNKNOWN;
  this->unsorted = UNKNOWN; this->nondeterministic = UNKNOWN;
}

void
PathGraph::swap(PathGraph& another) noexcept
{
  this->path_names.swap(another.path_names);
  this->rank_names.swap(another.rank_names);
  this->path_counts.swap(another.path_counts);
  this->rank_counts.swap(another.rank_counts);
  this->logical_file_ids.swap(another.logical_file_ids);
  this->physical_shard_ids.swap(another.physical_shard_ids);

  std::swap(this->path_count, another.path_count);
  std::swap(this->rank_count, another.rank_count);
  std::swap(this->range_count, another.range_count);
  std::swap(this->order, another.order);
  std::swap(this->doubling_steps, another.doubling_steps);

  std::swap(this->unique, another.unique);
  std::swap(this->redundant, another.redundant);
  std::swap(this->unsorted, another.unsorted);
  std::swap(this->nondeterministic, another.nondeterministic);
  std::swap(this->delete_files, another.delete_files);
}

void
PathGraph::open(std::ifstream& path_file, std::ifstream& rank_file, size_type file) const
{
  if(file >= this->files())
  {
    std::cerr << "PathGraph::open(): Invalid file number: " << file << std::endl;
    std::exit(EXIT_FAILURE);
  }

  path_file.open(this->path_names[file].c_str(), std::ios_base::binary);
  if(!path_file)
  {
    std::cerr << "PathGraph::open(): Cannot open path file " << this->path_names[file] << std::endl;
    std::exit(EXIT_FAILURE);
  }

  rank_file.open(this->rank_names[file].c_str(), std::ios_base::binary);
  if(!rank_file)
  {
    std::cerr << "PathGraph::open(): Cannot open rank file " << this->rank_names[file] << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

//------------------------------------------------------------------------------

struct SameFromLogicalFile
{
  const PathGraphMerger& merger;
  node_type              from;
  logical_file_id_t      logical_file;
  bool                   same_from, same_file;

  SameFromLogicalFile(const PathGraphMerger& source, range_type range) :
    merger(source), from(source.buffer.get(range.first).node.from),
    logical_file(source.graph.logicalFile(source.buffer.get(range.first).file)),
    same_from(true), same_file(true)
  {
    for(size_type i = range.first + 1; i <= range.second; i++)
    {
      PriorityNode node = this->merger.buffer.get(i);
      if(node.node.from != this->from) { this->same_from = false; }
      if(this->merger.graph.logicalFile(node.file) != this->logical_file)
      {
        this->same_file = false;
      }
    }
  }

  inline bool operator() (range_type range) const
  {
    for(size_type i = range.first; i <= range.second; i++)
    {
      PriorityNode node = this->merger.buffer.get(i);
      if(node.node.from != this->from ||
         this->merger.graph.logicalFile(node.file) != this->logical_file)
      {
        return false;
      }
    }
    return true;
  }
};

void
PathGraph::prune(const LCP& lcp, size_type size_limit,
  size_type group_buffer_bytes, PathGraphMergeStats* stats)
{
  size_type old_path_count = this->size();

  PathGraphMerger merger(*this, lcp, group_buffer_bytes, stats);
  PathGraphBuilder builder(this->files(), this->k(), this->step(), size_limit);
  builder.graph.logical_file_ids = this->logical_file_ids;
  builder.graph.physical_shard_ids = this->physical_shard_ids;
  for(range_type range = merger.first(); !(merger.atEnd(range)); range = merger.next())
  {
    SameFromLogicalFile same_from(merger, range);
    if(same_from.same_from)
    {
      if(same_from.same_file)
      {
        range = merger.extendRange(same_from);
        merger.mergePathNodes();
        PriorityNode node = merger.buffer.get(range.second);
        builder.write(node);
        builder.graph.unique++;
      }
      else
      {
        // FIXME Later: Write just one path per file.
        for(size_type i = range.first; i <= range.second; i++)
        {
          PriorityNode node = merger.buffer.get(i);
          node.node.makeSorted(); builder.write(node);
        }
        builder.graph.redundant += Range::length(range);
      }
    }
    else
    {
      for(size_type i = range.first; i <= range.second; i++)
      {
        PriorityNode node = merger.buffer.get(i);
        if(node.node.sorted()) { builder.graph.nondeterministic++; }
        else { builder.graph.unsorted++; }
        builder.write(node);
      }
    }
    builder.graph.range_count++;
  }
  merger.close(); builder.close();
  this->clear(); this->swap(builder.graph);

  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    std::cerr << "PathGraph::prune(): " << old_path_count << " -> " << this->size() << " paths ("
              << this->ranges() << " ranges)" << std::endl;
    std::cerr << "PathGraph::prune(): "
              << this->unique << " unique, "
              << this->redundant << " redundant, "
              << this->unsorted << " unsorted, "
              << this->nondeterministic << " nondeterministic paths" << std::endl;
    std::cerr << "PathGraph::prune(): " << inGigabytes(this->bytes()) << " GB in "
              << this->files() << " file(s)" << std::endl;
  }
}

//------------------------------------------------------------------------------

void
PathGraph::extend(size_type size_limit, size_type memory_limit)
{
  size_type old_path_count = this->size();

  PathGraphBuilder builder(this->files(), 2 * this->k(), this->step() + 1, size_limit);
  builder.graph.logical_file_ids = this->logical_file_ids;
  builder.graph.physical_shard_ids = this->physical_shard_ids;
  for(size_type file = 0; file < this->files(); file++)
  {
    // Read the current file.
    std::vector<PathNode> paths;
    std::vector<PathNode::rank_type> labels;
    this->read(paths, labels, file);

    // Initialization.
    PathFromComparator from_c;  // Sort the paths by from.
    parallelQuickSort(paths.begin(), paths.end(), from_c);
    ValueIndex<PathNode, FromGetter> from_index(paths);
    size_type threads = omp_get_max_threads();

    // Create thread-specific buffers.
    std::vector<std::vector<PathNode>> temp_nodes(threads);
    std::vector<std::vector<PathNode::rank_type>> temp_labels(threads);
    for(size_type thread = 0; thread < threads; thread++)
    {
      temp_nodes[thread].reserve(PathGraphBuilder::WRITE_BUFFER_SIZE);
      temp_labels[thread].reserve(((1 << builder.graph.step()) + 1) * PathGraphBuilder::WRITE_BUFFER_SIZE);
    }

    // Create the next generation.
    size_type chunk_size = getChunkSize(paths.size(), MEGABYTE);
    #pragma omp parallel for schedule(dynamic, chunk_size)
    for(size_type i = 0; i < paths.size(); i++)
    {
      size_type thread = omp_get_thread_num();
      if(paths[i].sorted())
      {
        temp_nodes[thread].push_back(PathNode(paths[i], labels, temp_labels[thread]));
        if(temp_nodes[thread].size() >= PathGraphBuilder::WRITE_BUFFER_SIZE)
        {
          builder.write(temp_nodes[thread], temp_labels[thread], memory_limit, file);
        }
      }
      else
      {
        size_type first = from_index.find(paths[i].to);
        for(size_type j = first; j < paths.size() && paths[j].from == paths[i].to; j++)
        {
          temp_nodes[thread].push_back(PathNode(paths[i], paths[j], labels, temp_labels[thread]));
          if(temp_nodes[thread].size() >= PathGraphBuilder::WRITE_BUFFER_SIZE)
          {
            builder.write(temp_nodes[thread], temp_labels[thread], memory_limit, file);
          }
        }
      }
    }
    for(size_type thread = 0; thread < threads; thread++)
    {
      builder.write(temp_nodes[thread], temp_labels[thread], memory_limit, file);
    }
    sdsl::util::clear(paths); sdsl::util::clear(labels);

    if(Verbosity::level >= Verbosity::FULL)
    {
      std::cerr << "PathGraph::extend(): File " << file << ": Created " << builder.graph.path_counts[file]
                << " order-" << builder.graph.k() << " paths" << std::endl;
    }

    // Sort the next generation.
    // The legacy extension route still delegates label ordering to the external
    // sorter. Give it a conservative fraction of the extension allowance, but
    // never less than the sorter's explicit stream/buffer minimum.
    if(memory_limit < externalPathGraphSortMinimumBudget())
    {
      externalSortFailure("extension memory limit is below the label-sort minimum");
    }
    size_type sort_budget = std::max(externalPathGraphSortMinimumBudget(), memory_limit / 8);
    builder.sort(file, sort_budget, 8);
  }
  builder.close();
  this->clear(); this->swap(builder.graph);

  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    std::cerr << "PathGraph::extend(): " << old_path_count << " -> " << this->size() << " paths ("
              << this->ranks() << " ranks)" << std::endl;
    std::cerr << "PathGraph::extend(): " << inGigabytes(this->bytes()) << " GB in "
              << this->files() << " file(s)" << std::endl;
  }
}

void
PathGraph::debugExtend()
{
  for(size_type file = 0; file < this->files(); file++)
  {
    // Read the current file.
    std::vector<PathNode> paths;
    std::vector<PathNode::rank_type> labels;
    this->read(paths, labels, file);

    // Initialization.
    PathFromComparator from_c;  // Sort the paths by from.
    parallelQuickSort(paths.begin(), paths.end(), from_c);
    ValueIndex<PathNode, FromGetter> from_index(paths);

    // Count the paths in the next generation.
    size_type sorted_paths = 0, unsorted_paths = 0;
    for(size_type i = 0; i < paths.size(); i++)
    {
      if(paths[i].sorted())
      {
        sorted_paths++;
      }
      else
      {
        size_type start = from_index.find(paths[i].to);
        size_type limit = start;
        while(limit < paths.size() && paths[limit].from == paths[i].to) { limit++; }
        unsorted_paths += limit - start;
        if(limit - start >= 100 && Verbosity::level >= Verbosity::FULL)
        {
          std::cerr << "PathGraph::debugExtend(): File " << file << ", path " << i
                    << " has " << (limit - start) << " extensions" << std::endl;
          std::cerr << "PathGraph::debugExtend(): The path is: ";
          paths[i].print(std::cerr, labels); std::cerr << std::endl;
        }
      }
    }

    if(Verbosity::level >= Verbosity::EXTENDED)
    {
      std::cerr << "PathGraph::debugExtend(): File " << file << ": " << this->size() << " -> "
                << (sorted_paths + unsorted_paths) << " paths ("
                << sorted_paths << " sorted, " << unsorted_paths << " unsorted)" << std::endl;
    }
  }
}

//------------------------------------------------------------------------------

void
PathGraph::read(std::vector<PathNode>& paths, std::vector<PathNode::rank_type>& labels, size_type file) const
{
  paths.resize(this->path_counts[file]);
  labels.resize(this->rank_counts[file]);

  std::ifstream path_file, rank_file;
  this->open(path_file, rank_file, file);
  if(!DiskIO::read(path_file, paths.data(), this->path_counts[file]))
  {
    std::cerr << "PathGraph::read(): Unexpected EOF in " << this->path_names[file] << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if(!DiskIO::read(rank_file, labels.data(), this->rank_counts[file]))
  {
    std::cerr << "PathGraph::read(): Unexpected EOF in " << this->rank_names[file] << std::endl;
    std::exit(EXIT_FAILURE);
  }
  path_file.close(); rank_file.close();

  if(Verbosity::level >= Verbosity::FULL)
  {
    std::cerr << "PathGraph::read(): File " << file << ": Read " << paths.size()
              << " order-" << this->k()<< " paths" << std::endl;
  }
}

//------------------------------------------------------------------------------

template<class Element>
struct SequentialRecordWriter
{
  std::string name;
  int descriptor;
  std::vector<Element> buffer;
  off_t written, released;

  constexpr static off_t CACHE_TAIL = 8 * MEGABYTE;
  constexpr static off_t CACHE_FLUSH = 64 * MEGABYTE;

  SequentialRecordWriter(const std::string& filename, size_type buffer_bytes) :
    name(filename), descriptor(-1), buffer(), written(0), released(0)
  {
    size_type records = std::max(static_cast<size_type>(1), buffer_bytes / sizeof(Element));
    this->buffer.reserve(records);
    this->descriptor = ::open(this->name.c_str(), O_CREAT | O_TRUNC | O_WRONLY, 0600);
    if(this->descriptor < 0) { throw std::runtime_error("SameFromSet: cannot create raw set file"); }
#if defined(POSIX_FADV_SEQUENTIAL)
    static_cast<void>(::posix_fadvise(this->descriptor, 0, 0, POSIX_FADV_SEQUENTIAL));
#endif
  }

  ~SequentialRecordWriter() { if(this->descriptor >= 0) { ::close(this->descriptor); } }

  void pushBack(const Element& value)
  {
    if(this->buffer.size() == this->buffer.capacity()) { this->flush(); }
    this->buffer.push_back(value);
  }

  void close()
  {
    if(this->descriptor < 0) { return; }
    this->flush(); this->syncAndTrim(true);
    if(::close(this->descriptor) != 0)
    {
      this->descriptor = -1;
      throw std::runtime_error("SameFromSet: cannot close raw set file");
    }
    this->descriptor = -1;
  }

private:
  void flush()
  {
    if(this->buffer.empty()) { return; }
    const char* data = reinterpret_cast<const char*>(this->buffer.data());
    size_type bytes = this->buffer.size() * sizeof(Element), done = 0;
    while(done < bytes)
    {
      ssize_t result = ::write(this->descriptor, data + done, bytes - done);
      if(result < 0 && errno == EINTR) { continue; }
      if(result <= 0) { throw std::runtime_error("SameFromSet: raw set write failed"); }
      DiskIO::write_volume += static_cast<size_type>(result);
      done += static_cast<size_type>(result);
    }
    this->written += static_cast<off_t>(bytes); this->buffer.clear();
    this->syncAndTrim(false);
  }

  void syncAndTrim(bool complete)
  {
    if(!complete && this->written - this->released < CACHE_FLUSH) { return; }
    if(::fdatasync(this->descriptor) != 0)
    {
      throw std::runtime_error("SameFromSet: raw set sync failed");
    }
    off_t discard_end = (complete ? this->written : this->written - CACHE_TAIL);
#if defined(POSIX_FADV_DONTNEED)
    if(discard_end > this->released)
    {
      static_cast<void>(::posix_fadvise(this->descriptor, this->released,
        discard_end - this->released, POSIX_FADV_DONTNEED));
    }
#endif
    this->released = discard_end;
  }
};

template<class Element>
constexpr off_t SequentialRecordWriter<Element>::CACHE_TAIL;
template<class Element>
constexpr off_t SequentialRecordWriter<Element>::CACHE_FLUSH;

template<class Element>
struct SequentialRecordReader
{
  int descriptor;
  std::vector<Element> buffer;
  size_type total, loaded, offset;
  off_t released;

  SequentialRecordReader(const std::string& filename, size_type buffer_bytes) :
    descriptor(::open(filename.c_str(), O_RDONLY)), buffer(),
    total(0), loaded(0), offset(0), released(0)
  {
    if(this->descriptor < 0) { throw std::runtime_error("SameFromSet: cannot open set file"); }
    struct stat info;
    if(::fstat(this->descriptor, &info) != 0 || info.st_size < 0 ||
       info.st_size % static_cast<off_t>(sizeof(Element)) != 0)
    {
      ::close(this->descriptor); this->descriptor = -1;
      throw std::runtime_error("SameFromSet: invalid set file");
    }
    this->total = static_cast<size_type>(info.st_size / sizeof(Element));
    size_type records = std::max(static_cast<size_type>(1), buffer_bytes / sizeof(Element));
    this->buffer.reserve(records);
#if defined(POSIX_FADV_SEQUENTIAL)
    static_cast<void>(::posix_fadvise(this->descriptor, 0, 0, POSIX_FADV_SEQUENTIAL));
#endif
  }

  ~SequentialRecordReader()
  {
    if(this->descriptor >= 0)
    {
      this->discard(static_cast<off_t>(this->loaded * sizeof(Element)));
      ::close(this->descriptor);
    }
  }

  bool next(Element& value)
  {
    if(this->offset == this->buffer.size())
    {
      if(this->loaded == this->total) { return false; }
      this->refill();
    }
    value = this->buffer[this->offset++]; return true;
  }

private:
  void refill()
  {
    this->discard(static_cast<off_t>(this->loaded * sizeof(Element)));
    size_type count = std::min(static_cast<size_type>(this->buffer.capacity()),
      this->total - this->loaded);
    this->buffer.resize(count); size_type bytes = count * sizeof(Element), done = 0;
    off_t file_offset = static_cast<off_t>(this->loaded * sizeof(Element));
    while(done < bytes)
    {
      ssize_t result = ::pread(this->descriptor,
        reinterpret_cast<char*>(this->buffer.data()) + done, bytes - done,
        file_offset + static_cast<off_t>(done));
      if(result < 0 && errno == EINTR) { continue; }
      if(result <= 0) { throw std::runtime_error("SameFromSet: truncated set file"); }
      DiskIO::read_volume += static_cast<size_type>(result);
      done += static_cast<size_type>(result);
    }
    this->loaded += count; this->offset = 0;
  }

  void discard(off_t through)
  {
#if defined(POSIX_FADV_DONTNEED)
    if(through > this->released)
    {
      static_cast<void>(::posix_fadvise(this->descriptor, this->released,
        through - this->released, POSIX_FADV_DONTNEED));
    }
#endif
    this->released = through;
  }
};

struct SameFromSet
{
  const PathGraphMerger& merger;
  std::string            selected;
  size_type              selected_nodes, budget, stream_buffer;
  PathGraphMergeStats*   stats;

  SameFromSet(const PathGraphMerger& source, size_type group_buffer_bytes,
    PathGraphMergeStats* merge_stats) :
    merger(source), selected(), selected_nodes(0),
    budget(std::max(ExternalFixedRecordSorter::minimumBudget(sizeof(node_type)),
      group_buffer_bytes)),
    stream_buffer(std::max(static_cast<size_type>(sizeof(node_type)),
      std::min(static_cast<size_type>(64 * KILOBYTE), this->budget / 4))),
    stats(merge_stats)
  {
  }

  ~SameFromSet() { if(!this->selected.empty()) { TempFile::remove(this->selected); } }

  std::pair<std::string, size_type> fromNodes(range_type range) const
  {
    std::string raw = TempFile::getName("gcsa_same_from_raw");
    std::string reduced = TempFile::getName("gcsa_same_from_set");
    size_type count = 0;
    try
    {
      {
        SequentialRecordWriter<node_type> output(raw, this->stream_buffer);
        for(size_type i = range.first; i <= range.second; i++)
        {
          output.pushBack(this->merger.buffer.get(i).node.from);
        }
        output.close();
      }
      if(this->stats != nullptr) { this->stats->from_set_sorts++; }
      ExternalFixedRecordSorter::sortAndReduce(raw, reduced, sizeof(node_type), this->budget, 2,
        [](const void* left, const void* right) {
          node_type a, b; std::memcpy(&a, left, sizeof(a)); std::memcpy(&b, right, sizeof(b));
          return (a < b ? -1 : (a > b ? 1 : 0));
        },
        [&count](const void* record, bool first, bool, std::ostream& stream) {
          if(first) { stream.write(reinterpret_cast<const char*>(record), sizeof(node_type)); count++; }
        });
    }
    catch(...)
    {
      TempFile::remove(raw); TempFile::remove(reduced); throw;
    }
    TempFile::remove(raw); return std::make_pair(reduced, count);
  }

  bool operator() (range_type range)
  {
    std::pair<std::string, size_type> next_set = this->fromNodes(range);
    bool equal = (next_set.second == this->selected_nodes);
    try
    {
      if(equal)
      {
        SequentialRecordReader<node_type> left(this->selected, this->stream_buffer);
        SequentialRecordReader<node_type> right(next_set.first, this->stream_buffer);
        node_type a, b;
        for(size_type i = 0; equal && i < this->selected_nodes; i++)
        {
          equal = (left.next(a) && right.next(b) && a == b);
        }
        if(equal) { equal = (!left.next(a) && !right.next(b)); }
      }
    }
    catch(...)
    {
      TempFile::remove(next_set.first); throw;
    }
    TempFile::remove(next_set.first); return equal;
  }

  void select(range_type range)
  {
    std::pair<std::string, size_type> result = this->fromNodes(range);
    if(result.second == 0)
    {
      TempFile::remove(result.first);
      throw std::runtime_error("SameFromSet: empty selected set");
    }
    if(!this->selected.empty()) { TempFile::remove(this->selected); }
    this->selected = result.first; this->selected_nodes = result.second;
  }

  template<class Callback>
  node_type streamAfterFirst(Callback callback) const
  {
    SequentialRecordReader<node_type> input(this->selected, this->stream_buffer);
    node_type first_value;
    if(!input.next(first_value)) { throw std::runtime_error("SameFromSet: empty selected set"); }
    node_type value;
    size_type seen = 1;
    while(input.next(value)) { callback(value); seen++; }
    if(seen != this->selected_nodes)
    {
      throw std::runtime_error("SameFromSet: selected set size changed");
    }
    return first_value;
  }
};

MergedGraph::MergedGraph(const PathGraph& source, const DeBruijnGraph& mapper,
  const LCP& kmer_lcp, size_type size_limit, size_type group_buffer_bytes,
  PathGraphMergeStats* stats) :
  path_name(TempFile::getName(PREFIX)), rank_name(TempFile::getName(PREFIX)),
  from_name(TempFile::getName(PREFIX)), lcp_name(TempFile::getName(PREFIX)),
  path_count(0), rank_count(0), from_count(0),
  order(source.k()),
  next(mapper.alpha.sigma + 1, 0), next_from(mapper.alpha.sigma + 1, 0)
{
  WriteBuffer<PathNode>            path_file(this->path_name);
  WriteBuffer<PathNode::rank_type> rank_file(this->rank_name);
  WriteBuffer<range_type>          from_file(this->from_name);
  WriteBuffer<uint8_t>             lcp_file(this->lcp_name);

  /*
     Initialize next[comp] to be the the rank of the first kmer starting with
     the corresponding character. Later, next[comp] is transformed into the rank
     of the first path with firstLabel(0) >= next[comp].
  */
  for(size_type comp = 0; comp < mapper.alpha.sigma; comp++)
  {
    this->next[comp] = mapper.charRange(comp).first;
  }
  this->next[mapper.alpha.sigma] = ~(size_type)0;
  this->next_from[mapper.alpha.sigma] = ~(size_type)0;

  PathGraphMerger merger(source, kmer_lcp, group_buffer_bytes, stats);
  SameFromSet same_from_set(merger, group_buffer_bytes, stats);
  size_type curr_comp = 0;  // Used to transform next.

  size_type bytes = 0;
  for(range_type range = merger.first(); !(merger.atEnd(range)); range = merger.next())
  {
    range_type path_lcp = merger.ranges.front().left_lcp; // Write this to the LCP array.
    same_from_set.select(range);
    range = merger.extendRange(same_from_set);
    merger.mergePathNodes();
    PriorityNode curr = merger.buffer.get(range.second);

    // Write the actual data
    bytes += curr.node.bytes() + (same_from_set.selected_nodes - 1) * sizeof(range_type) + 1;
    if(bytes > size_limit)
    {
      std::cerr << "MergedGraph::MergedGraph(): Size limit exceeded, construction aborted" << std::endl;
      std::exit(EXIT_SIZE_LIMIT_EXCEEDED);
    }
    curr.node.from = same_from_set.streamAfterFirst(
      [&](node_type from) { from_file.push_back(range_type(this->path_count, from)); });
    writePath(curr.node, curr.label, path_file, rank_file);
    lcp_file.push_back(path_lcp.first * mapper.order() + path_lcp.second);

    // Update the counts and the pointers to paths starting with each comp value.
    while(curr.firstLabel(0) >= this->next[curr_comp])
    {
      this->next[curr_comp] = this->path_count;
      this->next_from[curr_comp] = this->from_count;
      curr_comp++;
    }
    this->path_count++;
    this->rank_count += curr.node.ranks();
    this->from_count += same_from_set.selected_nodes - 1;
  }
  merger.close();
  path_file.close(); rank_file.close(); from_file.close(); lcp_file.close();

  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    std::cerr << "MergedGraph::MergedGraph(): " << this->size() << " paths with "
              << this->ranks() << " ranks and " << this->extra() << " additional start nodes" << std::endl;
    std::cerr << "MergedGraph::MergedGraph(): " << inGigabytes(this->bytes()) << " GB" << std::endl;
  }
}

MergedGraph::~MergedGraph()
{
  this->clear();
}

void
MergedGraph::clear()
{
  TempFile::remove(this->path_name);
  TempFile::remove(this->rank_name);
  TempFile::remove(this->from_name);
  TempFile::remove(this->lcp_name);

  this->path_count = 0; this->rank_count = 0; this->from_count = 0;
  this->order = 0;

  for(size_type i = 0; i < this->next.size(); i++) { this->next[i] = 0; }
  for(size_type i = 0; i < this->next_from.size(); i++) { this->next_from[i] = 0; }
}

//------------------------------------------------------------------------------

} // namespace gcsa
