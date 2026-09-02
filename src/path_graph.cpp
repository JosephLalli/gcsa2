#include <gcsa/path_graph.h>

#include <sdsl/wt_algorithm.hpp>

#include <algorithm>
#include <cstdio>
#include <deque>
#include <fcntl.h>
#include <limits>
#include <queue>
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

/* A run record is deliberately self-contained: the rank pointer in node is
   meaningful only while reading the source/final PathGraph pair. */
struct PathSortRecord
{
  PathNode node;
  PathNode::rank_type labels[PathLabel::LABEL_LENGTH + 1];
};

// Covers stream state, vector control blocks, allocator slack, and the heap.
constexpr size_type PATH_SORT_FIXED_BYTES = 8192;
// Linux charges clean filesystem cache to a cgroup's memory ceiling. The
// sorter therefore retains only a small sequential tail and periodically
// syncs completed output prefixes before making them reclaimable. These
// values intentionally match the external join policy.
constexpr off_t PATH_SORT_CACHE_TAIL_BYTES = 64 * MEGABYTE;
constexpr off_t PATH_SORT_CACHE_FLUSH_BYTES = 512 * MEGABYTE;
constexpr size_type PATH_SORT_SOURCE_CACHE_CHECK_RECORDS = 1024 * 1024;

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

struct PathSortRunReader
{
  std::ifstream file;
  std::vector<PathSortRecord> buffer;
  size_type offset, buffer_records, total_records;
  bool at_end;
  int cache_descriptor;
  off_t bytes_read, cache_released;

  PathSortRunReader(const std::string& name, size_type requested_records) :
    file(), buffer(), offset(0), buffer_records(1), total_records(0), at_end(false),
    cache_descriptor(-1), bytes_read(0), cache_released(0)
  {
    // The explicit byte-counted vector is the only input buffer. Otherwise
    // libstdc++ adds an unaccounted buffer for every merge input stream.
    this->file.rdbuf()->pubsetbuf(nullptr, 0);
    this->file.open(name.c_str(), std::ios_base::binary);
    if(!this->file) { externalSortFailure("cannot open run " + name); }
    size_type bytes = fileSize(this->file);
    if(bytes % sizeof(PathSortRecord) != 0) { externalSortFailure("truncated run " + name); }
    this->total_records = bytes / sizeof(PathSortRecord);
    this->cache_descriptor = openPathSortCacheDescriptor(name, O_RDONLY);
    this->buffer_records = std::max((size_type)1,
      std::min(requested_records, this->total_records));
    this->buffer.resize(this->buffer_records);
    this->refill();
  }

  PathSortRunReader(PathSortRunReader&& another) noexcept :
    file(std::move(another.file)), buffer(std::move(another.buffer)),
    offset(another.offset), buffer_records(another.buffer_records),
    total_records(another.total_records), at_end(another.at_end),
    cache_descriptor(another.cache_descriptor), bytes_read(another.bytes_read),
    cache_released(another.cache_released)
  {
    another.cache_descriptor = -1;
  }

  ~PathSortRunReader()
  {
    if(this->cache_descriptor >= 0) { ::close(this->cache_descriptor); }
  }

  void refill()
  {
    bool complete = (this->bytes_read == pathSortByteOffset(this->total_records, sizeof(PathSortRecord)));
    trimPathSortReadCache(this->cache_descriptor, this->bytes_read, this->cache_released, complete);
    this->file.read(reinterpret_cast<char*>(this->buffer.data()), this->buffer.size() * sizeof(PathSortRecord));
    std::streamsize bytes = this->file.gcount();
    DiskIO::read_volume += bytes;
    if(bytes % (std::streamsize)sizeof(PathSortRecord) != 0) { externalSortFailure("truncated run"); }
    this->bytes_read += bytes;
    this->buffer.resize(bytes / sizeof(PathSortRecord));
    this->offset = 0;
    this->at_end = this->buffer.empty();
  }

  const PathSortRecord& current() const { return this->buffer[this->offset]; }

  void advance()
  {
    this->offset++;
    if(this->offset == this->buffer.size())
    {
      this->buffer.resize(this->buffer_records);
      this->refill();
    }
  }

  PathSortRunReader(const PathSortRunReader&) = delete;
  PathSortRunReader& operator=(const PathSortRunReader&) = delete;
  PathSortRunReader& operator=(PathSortRunReader&&) = delete;
};

struct PathSortHeapComparator
{
  const std::vector<PathSortRunReader>* readers;

  explicit PathSortHeapComparator(const std::vector<PathSortRunReader>* _readers) : readers(_readers) { }

  bool operator() (size_type a, size_type b) const
  {
    const PathSortRecord& left = (*this->readers)[a].current();
    const PathSortRecord& right = (*this->readers)[b].current();
    if(pathSortLess(left, right)) { return false; }
    if(pathSortLess(right, left)) { return true; }
    return (a > b);
  }
};

void
writePathSortRun(const std::string& name, const std::vector<PathSortRecord>& records)
{
  std::ofstream output;
  output.rdbuf()->pubsetbuf(nullptr, 0);
  output.open(name.c_str(), std::ios_base::binary);
  if(!output) { externalSortFailure("cannot create run " + name); }
  int cache_descriptor = openPathSortCacheDescriptor(name, O_RDWR);
  off_t written = 0, cache_released = 0;
  if(!records.empty()) { DiskIO::write(output, records.data(), records.size()); }
  written = pathSortByteOffset(records.size(), sizeof(PathSortRecord));
  trimPathSortWrittenCache(output, cache_descriptor, written, cache_released, true, name);
  output.close();
  ::close(cache_descriptor);
}

std::string
mergePathSortRuns(const std::vector<std::string>& inputs, size_type buffer_records,
  ExternalPathSortStats* stats)
{
  std::string output_name = TempFile::getName("gcsa_path_sort_run");
  std::ofstream output;
  output.rdbuf()->pubsetbuf(nullptr, 0);
  output.open(output_name.c_str(), std::ios_base::binary);
  if(!output) { externalSortFailure("cannot create merged run"); }
  int output_cache_descriptor = openPathSortCacheDescriptor(output_name, O_RDWR);
  off_t output_bytes = 0, output_cache_released = 0;
  std::vector<PathSortRunReader> readers;
  readers.reserve(inputs.size());
  size_type input_buffer_records = 0, total_records = 0;
  for(size_type i = 0; i < inputs.size(); i++)
  {
    readers.emplace_back(inputs[i], buffer_records);
    input_buffer_records += readers.back().buffer_records;
    total_records += readers.back().total_records;
  }
  size_type output_records = std::max((size_type)1,
    std::min(buffer_records, total_records));
  std::vector<PathSortRecord> output_buffer;
  output_buffer.reserve(output_records);
  updatePathSortStats(stats, input_buffer_records + output_records,
    PATH_SORT_FIXED_BYTES + (input_buffer_records + output_records) * sizeof(PathSortRecord) +
    inputs.size() * (sizeof(PathSortRunReader) + 2 * sizeof(size_type)));
  std::priority_queue<size_type, std::vector<size_type>, PathSortHeapComparator>
    queue{PathSortHeapComparator(&readers)};
  for(size_type i = 0; i < readers.size(); i++)
  {
    if(!readers[i].at_end) { queue.push(i); }
  }
  while(!queue.empty())
  {
    size_type best = queue.top(); queue.pop();
    output_buffer.push_back(readers[best].current()); readers[best].advance();
    if(!readers[best].at_end) { queue.push(best); }
    if(output_buffer.size() >= output_buffer.capacity())
    {
      DiskIO::write(output, output_buffer.data(), output_buffer.size());
      addPathSortBytes(output_bytes, output_buffer.size(), sizeof(PathSortRecord));
      trimPathSortWrittenCache(output, output_cache_descriptor, output_bytes,
        output_cache_released, false, output_name);
      output_buffer.clear();
    }
  }
  if(!output_buffer.empty())
  {
    DiskIO::write(output, output_buffer.data(), output_buffer.size());
    addPathSortBytes(output_bytes, output_buffer.size(), sizeof(PathSortRecord));
  }
  trimPathSortWrittenCache(output, output_cache_descriptor, output_bytes,
    output_cache_released, true, output_name);
  output.close();
  ::close(output_cache_descriptor);
  return output_name;
}

void
writeSortedPathPair(const std::string& run_name, const std::string& path_name,
  const std::string& rank_name, size_type expected_paths, size_type expected_ranks,
  size_type buffer_records, ExternalPathSortStats* stats)
{
  PathSortRunReader reader(run_name, buffer_records);
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
  updatePathSortStats(stats, reader.buffer_records + path_buffer_records,
    PATH_SORT_FIXED_BYTES + reader.buffer_records * sizeof(PathSortRecord) +
    path_buffer_records * sizeof(PathNode) +
    rank_buffer_records * sizeof(PathNode::rank_type));
  size_type path_count = 0, rank_count = 0;
  while(!reader.at_end)
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

void
externalPathGraphSort(PathGraph& graph, size_type file, size_type byte_budget, size_type fan_in,
  ExternalPathSortStats* stats)
{
  if(file >= graph.files()) { externalSortFailure("invalid file number"); }
  if(byte_budget < externalPathGraphSortMinimumBudget()) { externalSortFailure("byte budget is too small"); }
  if(stats != nullptr) { *stats = ExternalPathSortStats(); }
  size_type record_bytes = sizeof(PathSortRecord);
  size_type available = byte_budget - PATH_SORT_FIXED_BYTES;
  size_type reader_overhead = sizeof(PathSortRunReader) + 2 * sizeof(size_type);
  size_type max_fan_in = (available - record_bytes) / (record_bytes + reader_overhead);
  fan_in = std::min(std::max((size_type)2, fan_in), max_fan_in);
  if(fan_in < 2) { externalSortFailure("byte budget cannot support a two-way merge"); }
  size_type merge_available = available - fan_in * reader_overhead;
  size_type merge_records = merge_available / ((fan_in + 1) * record_bytes);
  if(merge_records == 0) { externalSortFailure("byte budget cannot buffer a merge"); }
  // Reserve two thirds for sorting records and one third for allocator / I/O slack.
  size_type run_records = std::max((size_type)1, available / (3 * record_bytes));

  std::ifstream paths, ranks;
  paths.rdbuf()->pubsetbuf(nullptr, 0); ranks.rdbuf()->pubsetbuf(nullptr, 0);
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
  updatePathSortStats(stats, run_records, PATH_SORT_FIXED_BYTES + run_records * record_bytes);
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
      std::sort(records.begin(), records.end(), pathSortLess);
      std::string name = TempFile::getName("gcsa_path_sort_run");
      writePathSortRun(name, records); records.clear();
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
    std::sort(records.begin(), records.end(), pathSortLess);
    std::string name = TempFile::getName("gcsa_path_sort_run");
    writePathSortRun(name, records);
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
    writePathSortRun(name, records); runs.push_back(name);
  }
  std::string merged = mergePathSortRuns(runs, merge_records, stats);
  for(size_type i = 0; i < runs.size(); i++) { TempFile::remove(runs[i]); }

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
  // One record for each side of a two-way merge, one output record, and fixed state.
  return PATH_SORT_FIXED_BYTES + 3 * sizeof(PathSortRecord) +
    2 * (sizeof(PathSortRunReader) + 2 * sizeof(size_type));
}

//------------------------------------------------------------------------------

struct PathGraphMerger;

struct PathRange
{
  size_type  from, to;
  range_type left_lcp, range_lcp, right_lcp;

  inline range_type range() const { return range_type(this->from, this->to); }
  inline size_type length() const { return this->to + 1 - this->from; }

  PathRange(size_type start, size_type stop, range_type _left_lcp, PathGraphMerger& merger);
};

/*
  This structure reads a buffered stream of PriorityNodes in sorted order and outputs a
  stream of ranges of PriorityNodes with the same label.
*/

struct PathGraphMerger
{
  const PathGraph&                              graph;
  const LCP&                                    lcp;

  // Buffers.
  std::deque<PathRange>                         ranges;
  BufferWindow<PriorityNode>                    buffer;

  // Priority queue.
  std::vector<ReadBuffer<PathNode>>             path_files;
  std::vector<ReadBuffer<PathNode::rank_type>>  rank_files;
  std::vector<size_type>                        offsets;
  PriorityQueue<PriorityNode>                   inputs;

  PathGraphMerger(const PathGraph& path_graph, const LCP& kmer_lcp);
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
    return this->lcp.min_lcp(this->buffer[i].node, this->buffer[j].node,
      this->buffer[i].label, this->buffer[j].label);
  }

  inline range_type border_lcp(size_type i, size_type j) const // i < j
  {
    return this->lcp.max_lcp(this->buffer[i].node, this->buffer[j].node,
      this->buffer[i].label, this->buffer[j].label);
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

PathGraphMerger::PathGraphMerger(const PathGraph& path_graph, const LCP& kmer_lcp) :
  graph(path_graph), lcp(kmer_lcp),
  path_files(path_graph.files()), rank_files(path_graph.files()),
  offsets(path_graph.files()), inputs(path_graph.files())
{
  for(size_type file = 0; file < path_graph.files(); file++)
  {
    this->path_files[file].open(path_graph.path_names[file]);
    this->rank_files[file].open(path_graph.rank_names[file]);
    this->offsets[file] = 0;
    this->inputs[file].file = file; this->read(this->inputs[file]);
  }
  this->inputs.heapify();
}

void
PathGraphMerger::close()
{
  sdsl::util::clear(this->ranges);
  sdsl::util::clear(this->buffer);

  for(size_type file = 0; file < this->graph.files(); file++)
  {
    this->path_files[file].close();
    this->rank_files[file].close();
  }
  this->path_files.clear();
  this->rank_files.clear();
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

  this->ranges.push_back(PathRange(0, this->rangeEnd(0), range_type(0, 0), *this));
  return this->ranges.front().range();
}

range_type
PathGraphMerger::next()
{
  PathRange temp = this->ranges.front(); this->ranges.pop_front();
  if(this->ranges.empty())
  {
    this->ranges.push_back(
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
      const PathRange& temp = this->ranges.back();
      this->ranges.push_back(PathRange(temp.to + 1, this->rangeEnd(temp.to + 1), temp.right_lcp, *this));
    }
    const PathRange& next_range = this->ranges[curr];
    if(next_range.from >= this->size()) { break; }

    // Is this a suffix tree node with the same start nodes as range and lcp > range.left_lcp?
    if(!(comp(next_range.range()))) { break; }
    parent_lcp = std::min(parent_lcp, next_range.left_lcp);
    if(parent_lcp <= range.left_lcp) { break; }
    if(next_range.right_lcp >= parent_lcp) { curr++; continue; }

    // Replace range with [range.from, next_range.to].
    range.to = next_range.to;
    range.range_lcp = parent_lcp; range.right_lcp = next_range.right_lcp;
    for(size_type i = 0; i <= curr; i++) { this->ranges.pop_front(); }
    this->ranges.push_front(range); curr = 1;
  }

  return range.range();
}

void
PathGraphMerger::mergePathNodes()
{
  PathRange& range = this->ranges.front();
  if(range.length() == 1)
  {
    this->buffer[range.to].node.makeSorted();
    return;
  }

  size_type order = range.range_lcp.first;
  if(range.range_lcp.second > 0)
  {
    range_type ranks(this->buffer[range.from].node.firstLabel(order, this->buffer[range.from].label),
      this->buffer[range.to].node.lastLabel(order, this->buffer[range.to].label));
    this->buffer[range.to].label[order] = ranks.first;
    this->buffer[range.to].label[order + 1] = ranks.second;
    order++;
  }

  this->buffer[range.to].node.makeSorted();
  this->buffer[range.to].node.setOrder(order);
  this->buffer[range.to].node.setLCP(range.range_lcp.first);
  for(size_type i = range.from; i < range.to; i++)
  {
    this->buffer[range.to].node.addPredecessors(this->buffer[i].node);
  }
}

size_type
PathGraphMerger::rangeEnd(size_type start)
{
  if(!(this->buffer.buffered(start))) { this->bufferNext(); }

  size_type stop = start;
  while(stop + 1 < this->size())
  {
    if(!(this->buffer.buffered(stop + 1))) { this->bufferNext(); }
    if(this->buffer[start] < this->buffer[stop + 1]) { break; }
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
    this->path_files[path.file].seek(this->offsets[path.file]);
    path.node = this->path_files[path.file][this->offsets[path.file]];
    this->rank_files[path.file].seek(path.node.pointer());
    for(size_type i = 0; i < path.node.ranks(); i++)
    {
      path.label[i] = this->rank_files[path.file][path.node.pointer() + i];
    }
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
    merger(source), from(source.buffer[range.first].node.from),
    logical_file(source.graph.logicalFile(source.buffer[range.first].file)),
    same_from(true), same_file(true)
  {
    for(size_type i = range.first + 1; i <= range.second; i++)
    {
      if(this->merger.buffer[i].node.from != this->from) { this->same_from = false; }
      if(this->merger.graph.logicalFile(this->merger.buffer[i].file) != this->logical_file)
      {
        this->same_file = false;
      }
    }
  }

  inline bool operator() (range_type range) const
  {
    for(size_type i = range.first; i <= range.second; i++)
    {
      if(this->merger.buffer[i].node.from != this->from ||
         this->merger.graph.logicalFile(this->merger.buffer[i].file) != this->logical_file)
      {
        return false;
      }
    }
    return true;
  }
};

void
PathGraph::prune(const LCP& lcp, size_type size_limit)
{
  size_type old_path_count = this->size();

  PathGraphMerger merger(*this, lcp);
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
        builder.write(merger.buffer[range.second]);
        builder.graph.unique++;
      }
      else
      {
        // FIXME Later: Write just one path per file.
        for(size_type i = range.first; i <= range.second; i++)
        {
          merger.buffer[i].node.makeSorted();
          builder.write(merger.buffer[i]);
        }
        builder.graph.redundant += Range::length(range);
      }
    }
    else
    {
      for(size_type i = range.first; i <= range.second; i++)
      {
        if(merger.buffer[i].node.sorted()) { builder.graph.nondeterministic++; }
        else { builder.graph.unsorted++; }
        builder.write(merger.buffer[i]);
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
    // Sorting is deliberately single-threaded and uses only a conservative
    // fraction of the extension memory allowance for run and merge buffers.
    size_type sort_budget = std::max((size_type)sizeof(PathSortRecord), memory_limit / 8);
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

struct SameFromSet
{
  const PathGraphMerger& merger;
  std::vector<node_type> nodes, buffer;

  SameFromSet(const PathGraphMerger& source) :
    merger(source)
  {
  }

  inline void fromNodes(range_type range, std::vector<node_type>& to) const
  {
    to.clear();
    node_type prev = ~(node_type)0;
    for(size_type i = range.first; i <= range.second; i++)
    {
      node_type curr = this->merger.buffer[i].node.from;
      if(curr != prev) { to.push_back(curr); prev = curr; }
    }
    if(to.size() > 1) { removeDuplicates(to, false); }
  }

  inline bool operator() (range_type range)
  {
    this->fromNodes(range, this->buffer);

    // Manual comparison guarantees using a single thread.
    if(this->buffer.size() != this->nodes.size()) { return false; }
    for(size_type i = 0; i < this->buffer.size(); i++)
    {
      if(this->buffer[i] != this->nodes[i]) { return false; }
    }
    return true;
  }

  inline void select(range_type range)
  {
    this->fromNodes(range, this->nodes);
  }
};

MergedGraph::MergedGraph(const PathGraph& source, const DeBruijnGraph& mapper, const LCP& kmer_lcp, size_type size_limit) :
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

  PathGraphMerger merger(source, kmer_lcp);
  SameFromSet same_from_set(merger);
  size_type curr_comp = 0;  // Used to transform next.

  size_type bytes = 0;
  for(range_type range = merger.first(); !(merger.atEnd(range)); range = merger.next())
  {
    range_type path_lcp = merger.ranges.front().left_lcp; // Write this to the LCP array.
    same_from_set.select(range);
    range = merger.extendRange(same_from_set);
    merger.mergePathNodes();
    PriorityNode& curr = merger.buffer[range.second];
    curr.node.from = same_from_set.nodes[0];

    // Write the actual data
    bytes += curr.node.bytes() + (same_from_set.nodes.size() - 1) * sizeof(range_type) + 1;
    if(bytes > size_limit)
    {
      std::cerr << "MergedGraph::MergedGraph(): Size limit exceeded, construction aborted" << std::endl;
      std::exit(EXIT_SIZE_LIMIT_EXCEEDED);
    }
    writePath(curr.node, curr.label, path_file, rank_file);
    for(size_type i = 1; i < same_from_set.nodes.size(); i++)
    {
      from_file.push_back(range_type(this->path_count, same_from_set.nodes[i]));
    }
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
    this->from_count += same_from_set.nodes.size() - 1;
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
