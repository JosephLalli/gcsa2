#include <gcsa/path_graph.h>
#include <gcsa/path_graph_external.h>
#include <gcsa/compressed_block.h>
#include <gcsa/external_sort.h>
#include <gcsa/path_sort_run.h>

#include <sdsl/wt_algorithm.hpp>

#include <algorithm>
#include <array>
#include <cerrno>
#include <condition_variable>
#include <deque>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fcntl.h>
#include <functional>
#include <limits>
#include <memory>
#include <queue>
#include <stdexcept>
#include <thread>
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

std::atomic<size_type> LCP::range_minimum_queries(0);
std::atomic<size_type> LCP::range_minimum_span(0);

namespace
{

/*
  The diagnostic range-minimum counters are summed per thread and folded into
  the shared atomics in batches. Parallel prune and merge workers each issue
  hundreds of millions of queries, and two fetch_adds per query on one shared
  cache line made the workers queue for the line rather than descend the
  wavelet tree. A thread folds its remainder in when it exits; a total read
  while threads run can miss at most one batch per thread.
*/
struct RangeMinimumTally
{
  constexpr static size_type BATCH = 4096;
  size_type queries = 0, span = 0;

  ~RangeMinimumTally() { this->flush(); }

  void flush()
  {
    if(this->queries == 0) { return; }
    LCP::range_minimum_queries.fetch_add(this->queries, std::memory_order_relaxed);
    LCP::range_minimum_span.fetch_add(this->span, std::memory_order_relaxed);
    this->queries = 0; this->span = 0;
  }
};

thread_local RangeMinimumTally range_minimum_tally;

inline void
tallyRangeMinimum(size_type left, size_type right)
{
  RangeMinimumTally& tally = range_minimum_tally;
  tally.queries++;
  tally.span += (right >= left ? right - left + 1 : 0);
  if(tally.queries >= RangeMinimumTally::BATCH) { tally.flush(); }
}

} // anonymous namespace

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
    // Diagnostic only, and on the library's hottest loop: 753 M descents in a
    // chr21 merge, where two relaxed fetch_adds measured 3.46 s against plain
    // increments. Gate them on the verbosity that prints them, so a production
    // run pays one well-predicted branch and the measurement is still there
    // when it is asked for.
    if(Verbosity::level >= Verbosity::EXTENDED) { tallyRangeMinimum(left, right); }
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
    // Diagnostic only, and on the library's hottest loop: 753 M descents in a
    // chr21 merge, where two relaxed fetch_adds measured 3.46 s against plain
    // increments. Gate them on the verbosity that prints them, so a production
    // run pays one well-predicted branch and the measurement is still there
    // when it is asked for.
    if(Verbosity::level >= Verbosity::EXTENDED) { tallyRangeMinimum(left, right); }
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

/*
  Pruning may preserve thousands of physical shards. Keeping a WriteBuffer pair
  open for every shard makes both resident buffers and descriptors proportional
  to that physical layout. This cache bounds both resources while retaining one
  append-only byte stream per shard. Counts in PathGraph remain authoritative,
  so reopening a shard cannot change rank pointers or semantic file identity.
*/
struct PathGraphOutputCache
{
  struct Entry
  {
    size_type file, stamp;
    int path, rank;
    off_t path_high, rank_high, path_released, rank_released;
    std::vector<PathNode> paths;
    std::vector<PathNode::rank_type> ranks;

    Entry() : file(PathGraph::UNKNOWN), stamp(0), path(-1), rank(-1),
      path_high(0), rank_high(0), path_released(0), rank_released(0),
      paths(), ranks() { }
  };

  PathGraph& graph;
  std::vector<Entry> entries;
  std::vector<std::uint64_t> path_digests, rank_digests;
  // Files rewritten by a sort after they were written; the sort records their
  // checksums, and the digests accumulated here describe the unsorted bytes.
  std::vector<bool> rewritten;
  bool closed;
  size_type clock, max_pairs, pair_buffer_bytes;
  PathGraphMergeStats* stats;

  constexpr static off_t CACHE_FLUSH_BYTES = 4 * MEGABYTE;

  PathGraphOutputCache(PathGraph& target, size_type total_buffer_bytes,
    size_type requested_pairs, PathGraphMergeStats* merge_stats) :
    graph(target), entries(),
    path_digests(target.files(), 1469598103934665603ULL),
    rank_digests(target.files(), 1469598103934665603ULL),
    rewritten(target.files(), false), closed(false), clock(0),
    max_pairs(std::max(static_cast<size_type>(1), std::min(requested_pairs,
      total_buffer_bytes / minimumPairBytes()))),
    // Cap the per-pair staging buffer. With output_pairs == 1 -- prune's normal
    // case -- total_buffer_bytes / max_pairs is the entire group buffer, and
    // these pages are genuinely written, so the whole share was resident. Once
    // the framed-codec path makes pathMergeInputBudget reach memory_limit/16
    // (1.44 GiB at the chr21 configuration) that one number is handed to four
    // independent structures which each size themselves from all of it.
    // Sequential writes saturate well below the flush watermark, past which the
    // pages are released anyway, so the watermark is the natural ceiling.
    pair_buffer_bytes(std::max(minimumPairBytes(),
      std::min(static_cast<size_type>(CACHE_FLUSH_BYTES),
        total_buffer_bytes / this->max_pairs))), stats(merge_stats)
  {
    this->entries.reserve(this->max_pairs);
  }

  ~PathGraphOutputCache()
  {
    // Construction code calls close() on the success path. During exception
    // unwinding, best-effort closure leaves the incomplete PathGraph temporary
    // files owned by its normal cleanup machinery.
    try { this->close(); } catch(...) { }
  }

  void append(PriorityNode& source)
  {
    const size_type file = source.file;
    if(file >= this->graph.files())
    {
      throw std::out_of_range("PathGraphBuilder: invalid output shard");
    }
    Entry& entry = this->get(file);
    this->makeRoom(entry, source.node.ranks());
    source.node.setPointer(this->graph.rank_counts[file]);
    entry.paths.push_back(source.node);
    for(size_type i = 0; i < source.node.ranks(); i++)
    {
      entry.ranks.push_back(source.label[i]);
    }
    this->advanceCounts(file, source.node.ranks());
    if(this->bufferedBytes(entry) >= this->pair_buffer_bytes) { this->flush(entry); }
  }

  void append(size_type file, PathNode& node,
    const PathNode::rank_type* labels)
  {
    if(file >= this->graph.files() || (node.ranks() > 0 && labels == nullptr))
    {
      throw std::invalid_argument("PathGraphBuilder: invalid output record");
    }
    Entry& entry = this->get(file);
    this->makeRoom(entry, node.ranks());
    node.setPointer(this->graph.rank_counts[file]);
    entry.paths.push_back(node);
    if(node.ranks() > 0)
    {
      entry.ranks.insert(entry.ranks.end(), labels, labels + node.ranks());
    }
    this->advanceCounts(file, node.ranks());
    if(this->bufferedBytes(entry) >= this->pair_buffer_bytes) { this->flush(entry); }
  }

  void closeFile(size_type file)
  {
    for(Entry& entry : this->entries)
    {
      if(entry.file == file) { this->close(entry); return; }
    }
  }

  void markRewritten(size_type file) { this->rewritten.at(file) = true; }

  void close()
  {
    if(this->closed) { return; }
    // A failed write/sync must not strand every descriptor after the first
    // failing entry. Close all entries, then rethrow the first I/O error.
    std::exception_ptr first_error;
    for(Entry& entry : this->entries)
    {
      try { this->close(entry); }
      catch(...) { if(first_error == nullptr) { first_error = std::current_exception(); } }
    }
    if(first_error != nullptr) { std::rethrow_exception(first_error); }
    for(size_type file = 0; file < this->graph.files(); file++)
    {
      if(this->rewritten[file]) { continue; }
      this->graph.path_checksums.at(file).record(this->graph.path_names[file],
        this->path_digests[file]);
      this->graph.rank_checksums.at(file).record(this->graph.rank_names[file],
        this->rank_digests[file]);
    }
    this->closed = true;
  }

private:
  Entry& get(size_type file)
  {
    if(this->closed) { throw std::logic_error("cannot append to a closed path writer"); }
    this->clock++;
    Entry* target = nullptr;
    for(Entry& entry : this->entries)
    {
      if(entry.file == file) { entry.stamp = this->clock; return entry; }
      if(entry.file == PathGraph::UNKNOWN) { target = &entry; }
    }
    if(target == nullptr && this->entries.size() < this->max_pairs)
    {
      this->entries.push_back(Entry()); target = &(this->entries.back());
      if(this->stats != nullptr)
      {
        this->stats->max_open_output_pairs = std::max(
          this->stats->max_open_output_pairs, static_cast<size_type>(this->entries.size()));
      }
    }
    if(target == nullptr)
    {
      target = &(*std::min_element(this->entries.begin(), this->entries.end(),
        [](const Entry& left, const Entry& right) { return left.stamp < right.stamp; }));
      this->close(*target);
    }
    this->open(*target, file); target->stamp = this->clock;
    return *target;
  }

  void open(Entry& entry, size_type file)
  {
    entry.file = file;
    entry.path = ::open(this->graph.path_names[file].c_str(),
      O_CREAT | O_WRONLY | O_APPEND, 0600);
    if(entry.path < 0)
    {
      entry.file = PathGraph::UNKNOWN;
      throw std::runtime_error("PathGraphBuilder: cannot open path output");
    }
    entry.rank = ::open(this->graph.rank_names[file].c_str(),
      O_CREAT | O_WRONLY | O_APPEND, 0600);
    if(entry.rank < 0)
    {
      ::close(entry.path); entry.path = -1; entry.file = PathGraph::UNKNOWN;
      throw std::runtime_error("PathGraphBuilder: cannot open rank output");
    }
    entry.path_high = ::lseek(entry.path, 0, SEEK_END);
    entry.rank_high = ::lseek(entry.rank, 0, SEEK_END);
    if(entry.path_high < 0 || entry.rank_high < 0)
    {
      this->close(entry);
      throw std::runtime_error("PathGraphBuilder: cannot seek output shard");
    }
    entry.path_released = entry.path_high; entry.rank_released = entry.rank_high;
#if defined(POSIX_FADV_SEQUENTIAL)
    static_cast<void>(::posix_fadvise(entry.path, 0, 0, POSIX_FADV_SEQUENTIAL));
    static_cast<void>(::posix_fadvise(entry.rank, 0, 0, POSIX_FADV_SEQUENTIAL));
#endif
    // Both arrays are preallocated for the largest individual PathNode before
    // sharing the remaining allowance. makeRoom() flushes before either vector
    // would grow, avoiding allocator doubling beyond the declared pool.
    size_type extra_bytes = this->pair_buffer_bytes - minimumPairBytes();
    size_type path_bytes = sizeof(PathNode) + extra_bytes / 2;
    size_type rank_bytes = maximumRankBytes() + (extra_bytes - extra_bytes / 2);
    entry.paths.reserve(std::max(static_cast<size_type>(1), path_bytes / sizeof(PathNode)));
    entry.ranks.reserve(std::max(static_cast<size_type>(1), rank_bytes / sizeof(PathNode::rank_type)));
  }

  void makeRoom(Entry& entry, size_type ranks)
  {
    if(ranks > PathLabel::LABEL_LENGTH + 1)
    {
      throw std::runtime_error("PathGraphBuilder: path label exceeds disk-record limit");
    }
    if(entry.paths.size() + 1 > entry.paths.capacity() ||
       entry.ranks.size() + ranks > entry.ranks.capacity())
    {
      this->flush(entry);
    }
    if(entry.paths.size() + 1 > entry.paths.capacity() ||
       entry.ranks.size() + ranks > entry.ranks.capacity())
    {
      throw std::runtime_error("PathGraphBuilder: output buffer reservation is too small");
    }
  }

  void close(Entry& entry)
  {
    if(entry.file == PathGraph::UNKNOWN) { return; }
    std::exception_ptr first_error;
    try
    {
      this->flush(entry);
      this->syncAndTrim(entry.path, entry.path_high, entry.path_released, true);
      this->syncAndTrim(entry.rank, entry.rank_high, entry.rank_released, true);
    }
    catch(...) { first_error = std::current_exception(); }
    int path_result = (entry.path >= 0 ? ::close(entry.path) : 0);
    int rank_result = (entry.rank >= 0 ? ::close(entry.rank) : 0);
    entry.path = entry.rank = -1; entry.file = PathGraph::UNKNOWN; entry.stamp = 0;
    if(first_error == nullptr && (path_result != 0 || rank_result != 0))
    {
      first_error = std::make_exception_ptr(
        std::runtime_error("PathGraphBuilder: cannot close output shard"));
    }
    if(first_error != nullptr) { std::rethrow_exception(first_error); }
  }

  void flush(Entry& entry)
  {
    if(!(entry.paths.empty()))
    {
      this->writeAll(entry.path, entry.paths.data(),
        entry.paths.size() * sizeof(PathNode));
      this->path_digests[entry.file] = BuildWorkspace::checksum(entry.paths.data(),
        entry.paths.size() * sizeof(PathNode), this->path_digests[entry.file]);
      entry.path_high += static_cast<off_t>(entry.paths.size() * sizeof(PathNode));
      entry.paths.clear();
    }
    if(!(entry.ranks.empty()))
    {
      this->writeAll(entry.rank, entry.ranks.data(),
        entry.ranks.size() * sizeof(PathNode::rank_type));
      this->rank_digests[entry.file] = BuildWorkspace::checksum(entry.ranks.data(),
        entry.ranks.size() * sizeof(PathNode::rank_type), this->rank_digests[entry.file]);
      entry.rank_high += static_cast<off_t>(entry.ranks.size() * sizeof(PathNode::rank_type));
      entry.ranks.clear();
    }
    this->syncAndTrim(entry.path, entry.path_high, entry.path_released, false);
    this->syncAndTrim(entry.rank, entry.rank_high, entry.rank_released, false);
  }

  void advanceCounts(size_type file, size_type ranks)
  {
    this->graph.path_counts[file]++; this->graph.path_count++;
    this->graph.rank_counts[file] += ranks; this->graph.rank_count += ranks;
  }

  static size_type bufferedBytes(const Entry& entry)
  {
    return entry.paths.size() * sizeof(PathNode) +
      entry.ranks.size() * sizeof(PathNode::rank_type);
  }

  constexpr static size_type maximumRankBytes()
  {
    return (PathLabel::LABEL_LENGTH + 1) * sizeof(PathNode::rank_type);
  }

  constexpr static size_type minimumPairBytes()
  {
    return sizeof(PathNode) + maximumRankBytes();
  }

  static void writeAll(int descriptor, const void* source, size_type bytes)
  {
    const char* data = reinterpret_cast<const char*>(source); size_type done = 0;
    while(done < bytes)
    {
      ssize_t result = ::write(descriptor, data + done, bytes - done);
      if(result < 0 && errno == EINTR) { continue; }
      if(result <= 0) { throw std::runtime_error("PathGraphBuilder: output write failed"); }
      DiskIO::write_volume += static_cast<size_type>(result);
      done += static_cast<size_type>(result);
    }
  }

  void syncAndTrim(int descriptor, off_t high, off_t& released, bool complete)
  {
    if(!complete && high - released < CACHE_FLUSH_BYTES) { return; }
    if(::fdatasync(descriptor) != 0)
    {
      throw std::runtime_error("PathGraphBuilder: output sync failed");
    }
    off_t tail = static_cast<off_t>(std::min(this->pair_buffer_bytes,
      static_cast<size_type>(std::numeric_limits<off_t>::max())));
    off_t discard_end = (complete ? high : std::max(released, high - tail));
#if defined(POSIX_FADV_DONTNEED)
    if(discard_end > released)
    {
      static_cast<void>(::posix_fadvise(descriptor, released,
        discard_end - released, POSIX_FADV_DONTNEED));
    }
#endif
    released = discard_end;
  }
};

constexpr off_t PathGraphOutputCache::CACHE_FLUSH_BYTES;

struct PathGraphBuilder
{
  PathGraph graph;
  PathGraphOutputCache output_files;
  size_type limit;  // Bytes of disk space.

  constexpr static size_type WRITE_BUFFER_SIZE = MEGABYTE;  // PathNodes per thread.

  PathGraphBuilder(size_type file_count, size_type path_order, size_type step,
    size_type size_limit, size_type writer_buffer_bytes = MEGABYTE,
    size_type max_writer_pairs = 16, PathGraphMergeStats* stats = nullptr);
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

PathGraphBuilder::PathGraphBuilder(size_type file_count, size_type path_order,
  size_type step, size_type size_limit, size_type writer_buffer_bytes,
  size_type max_writer_pairs, PathGraphMergeStats* stats) :
  graph(file_count, path_order, step),
  output_files(this->graph, writer_buffer_bytes, max_writer_pairs, stats),
  limit(size_limit)
{
  // Materialize every empty shard without retaining a descriptor. Some input
  // shards legitimately produce no paths, and later checkpoint/sort code must
  // still observe a valid empty path/rank pair.
  for(size_type file = 0; file < this->graph.files(); file++)
  {
    int path = ::open(this->graph.path_names[file].c_str(),
      O_CREAT | O_TRUNC | O_WRONLY, 0600);
    int rank = ::open(this->graph.rank_names[file].c_str(),
      O_CREAT | O_TRUNC | O_WRONLY, 0600);
    if(path < 0 || rank < 0)
    {
      if(path >= 0) { ::close(path); }
      if(rank >= 0) { ::close(rank); }
      throw std::runtime_error("PathGraphBuilder: cannot create output shard");
    }
    int path_result = ::close(path), rank_result = ::close(rank);
    if(path_result != 0 || rank_result != 0)
    {
      throw std::runtime_error("PathGraphBuilder: cannot close empty output shard");
    }
  }
}

void
PathGraphBuilder::close()
{
  this->output_files.close();
}

void
PathGraphBuilder::write(PriorityNode& path)
{
  if(this->graph.bytes() + path.bytes() > this->limit)
  {
    std::cerr << "PathGraphBuilder::write(): Size limit exceeded, construction aborted" << std::endl;
    std::exit(EXIT_SIZE_LIMIT_EXCEEDED);
  }
  this->output_files.append(path);
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
      size_type pointer = paths[i].pointer();
      if(pointer > labels.size() || paths[i].ranks() > labels.size() - pointer)
      {
        throw std::runtime_error("PathGraphBuilder::write(): invalid rank pointer");
      }
      this->output_files.append(file, paths[i], labels.data() + pointer);
    }
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
  this->output_files.closeFile(file);
  externalPathGraphSort(this->graph, file, byte_budget, fan_in);
  this->output_files.markRewritten(file);

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
     stats->run_compressed_bytes > std::numeric_limits<size_type>::max() - writer.bytes() ||
     stats->grouped_records > std::numeric_limits<size_type>::max() - writer.groupReferences() ||
     stats->group_headers > std::numeric_limits<size_type>::max() - writer.groupHeaders() ||
     stats->context_bytes_saved > std::numeric_limits<size_type>::max() - writer.contextBytesSaved())
  {
    externalSortFailure("path-sort run statistics overflow");
  }
  stats->run_uncompressed_bytes += writer.uncompressedBytes();
  stats->run_compressed_bytes += writer.bytes();
  stats->grouped_records += writer.groupReferences();
  stats->group_headers += writer.groupHeaders();
  stats->context_bytes_saved += writer.contextBytesSaved();
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
  size_type buffer_records, ExternalPathSortStats* stats,
  const TempFileCodecParameters& codec, std::uint64_t& path_digest,
  std::uint64_t& rank_digest)
{
  path_digest = rank_digest = 1469598103934665603ULL;
  PathSortRunReader reader(run_name, PATH_SORT_RUN_BUFFER_BYTES);
  if(codec.enabled())
  {
    const CompressedBlockWriter::Mode mode = CompressedBlockWriter::ZSTD;
    CompressedBlockWriter paths(path_name, codec.block_size, mode, codec.level,
      codec.workers);
    CompressedBlockWriter ranks(rank_name, codec.block_size, mode, codec.level,
      codec.workers);
    paths.trackStoredChecksum(); ranks.trackStoredChecksum();
    updatePathSortStats(stats, 1,
      PATH_SORT_FIXED_BYTES + reader.bufferBytes() + sizeof(PathSortRunReader) +
      2 * CompressedBlockWriter::workingMemoryEstimate(
        codec.block_size, mode, codec.level, codec.workers));
    size_type path_count = 0, rank_count = 0;
    while(!(reader.atEnd()))
    {
      PathNode node = reader.current().node;
      node.setPointer(rank_count);
      paths.writeRecord(&node, sizeof(node));
      if(node.ranks() > 0)
      {
        ranks.writeRecord(reader.current().labels,
          node.ranks() * sizeof(PathNode::rank_type));
      }
      path_count++; rank_count += node.ranks(); reader.advance();
    }
    paths.finish(); ranks.finish();
    path_digest = paths.storedChecksum(); rank_digest = ranks.storedChecksum();
    if(path_count != expected_paths || rank_count != expected_ranks)
    {
      externalSortFailure("compressed sorted path pair has incorrect counts");
    }
    return;
  }

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
      rank_digest = BuildWorkspace::checksum(rank_buffer.data(),
        rank_buffer.size() * sizeof(PathNode::rank_type), rank_digest);
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
      rank_digest = BuildWorkspace::checksum(rank_buffer.data(),
        rank_buffer.size() * sizeof(PathNode::rank_type), rank_digest);
      addPathSortBytes(rank_bytes, rank_buffer.size(), sizeof(PathNode::rank_type));
      trimPathSortWrittenCache(ranks, rank_cache_descriptor, rank_bytes,
        rank_cache_released, false, rank_name);
      rank_buffer.clear();
    }
    if(path_buffer.size() == path_buffer.capacity())
    {
      DiskIO::write(paths, path_buffer.data(), path_buffer.size());
      path_digest = BuildWorkspace::checksum(path_buffer.data(),
        path_buffer.size() * sizeof(PathNode), path_digest);
      addPathSortBytes(path_bytes, path_buffer.size(), sizeof(PathNode));
      trimPathSortWrittenCache(paths, path_cache_descriptor, path_bytes,
        path_cache_released, false, path_name);
      path_buffer.clear();
    }
  }
  if(!path_buffer.empty())
  {
    DiskIO::write(paths, path_buffer.data(), path_buffer.size());
    path_digest = BuildWorkspace::checksum(path_buffer.data(),
      path_buffer.size() * sizeof(PathNode), path_digest);
    addPathSortBytes(path_bytes, path_buffer.size(), sizeof(PathNode));
  }
  if(!rank_buffer.empty())
  {
    DiskIO::write(ranks, rank_buffer.data(), rank_buffer.size());
    rank_digest = BuildWorkspace::checksum(rank_buffer.data(),
      rank_buffer.size() * sizeof(PathNode::rank_type), rank_digest);
    addPathSortBytes(rank_bytes, rank_buffer.size(), sizeof(PathNode::rank_type));
  }
  trimPathSortWrittenCache(paths, path_cache_descriptor, path_bytes,
    path_cache_released, true, path_name);
  trimPathSortWrittenCache(ranks, rank_cache_descriptor, rank_bytes,
    rank_cache_released, true, rank_name);
  paths.close(); ranks.close();
  ::close(path_cache_descriptor); ::close(rank_cache_descriptor);
  if(paths.fail() || ranks.fail()) { externalSortFailure("cannot close sorted path pair"); }
  if(path_count != expected_paths || rank_count != expected_ranks) { externalSortFailure("sorted path pair has incorrect counts"); }
}

TempFileCodecParameters
boundedPathCodec(const TempFileCodecParameters& requested, size_type byte_budget)
{
  if(!requested.enabled()) { return requested; }
  TempFileCodecParameters result = requested;
  const size_type minimum_block = std::max(sizeof(PathNode),
    (PathLabel::LABEL_LENGTH + 1) * sizeof(PathNode::rank_type));
  result.block_size = std::max(minimum_block,
    std::min(result.block_size, byte_budget));
  auto fits = [&](size_type block_bytes) -> bool
  {
    const size_type one = CompressedBlockWriter::workingMemoryEstimate(
      block_bytes, CompressedBlockWriter::ZSTD, result.level, result.workers);
    return (one <= byte_budget / 2);
  };
  while(result.block_size > minimum_block && !fits(result.block_size))
  {
    result.block_size = std::max(minimum_block, result.block_size / 2);
  }
  if(!fits(result.block_size))
  {
    // AUTO is a feasibility policy: under a deliberately tiny test budget,
    // retain the raw writer rather than rejecting a build that otherwise fits.
    // Explicit ZSTD remains strict and reports that its codec workspace cannot
    // be admitted under the configured ceiling.
    if(requested.compression == TempCompression::AUTO)
    {
      result.compression = TempCompression::NONE;
      return result;
    }
    externalSortFailure("memory budget cannot admit two compressed path writers");
  }
  return result;
}

} // anonymous namespace

size_type
pathMergeInputBudget(const ConstructionParameters& parameters)
{
  return std::max(static_cast<size_type>(1),
    std::min(parameters.getIOBufferSize(),
      parameters.getMemoryLimitBytes() / 16));
}

size_type
pathMergeOutputPairs(size_type max_open_files, size_type logical_inputs)
{
  // Two descriptors for the possible spill files, then two per cache entry.
  // The output cache can never occupy more entries than there are logical
  // outputs, and prune() emits exactly one per logical input.
  const size_type cache_pairs = (max_open_files < 6 ? 2 : (max_open_files - 2) / 2);
  return std::max(static_cast<size_type>(1),
    std::min(logical_inputs, cache_pairs / 2));
}

size_type
pathMergeCeilingBudget(const ConstructionParameters& parameters)
{
  return std::max(static_cast<size_type>(1),
    parameters.getMemoryLimitBytes() / 16);
}

size_type
pathMergeInputPairs(size_type max_open_files, size_type logical_inputs)
{
  // Whatever the output cache cannot occupy belongs to the input side, capped
  // by the narrower final merge so that one committed block size serves both.
  // Both merges keep every admitted pair open for a whole pass and visit the
  // shards in round-robin label order, the access pattern an LRU degrades to a
  // 100% miss rate on the moment it holds fewer pairs than there are shards.
  // Splitting the allowance evenly instead left prune() 31 pairs for the 62
  // shards compaction had retained, and every 24-byte PathNode then cost a
  // whole block decode. MergedGraph spends fourteen descriptors on its four
  // sequential outputs, two spill files, and the two-way SameFromSet sorter.
  if(max_open_files < 16) { return 1; }
  const size_type cache_pairs = (max_open_files - 2) / 2;
  return std::max(static_cast<size_type>(1),
    std::min(cache_pairs - pathMergeOutputPairs(max_open_files, logical_inputs),
      (max_open_files - 14) / 2));
}

size_type
pathMergeInputPairs(const ConstructionParameters& parameters)
{
  // One logical input is the widest case: every further one claims an output
  // entry. Producers need that upper bound, because a block small enough for
  // it is small enough for every narrower merge.
  return pathMergeInputPairs(parameters.getMaxOpenFiles(),
    static_cast<size_type>(1));
}


size_type
pathGraphFramedPairBytes(const PathGraph& source)
{
  size_type result = 0;
  for(size_type file = 0; file < source.files(); file++)
  {
    size_type pair_bytes = 0;
    if(CompressedBlockReader::isFramed(source.path_names[file]))
    {
      pair_bytes += CompressedBlockReader::workingMemoryEstimate(
        CompressedBlockReader::declaredBlockSize(source.path_names[file]));
    }
    if(CompressedBlockReader::isFramed(source.rank_names[file]))
    {
      pair_bytes += CompressedBlockReader::workingMemoryEstimate(
        CompressedBlockReader::declaredBlockSize(source.rank_names[file]));
    }
    result = std::max(result, pair_bytes);
  }
  return result;
}

size_type
pathMergeInputBudget(const ConstructionParameters& parameters,
  const PathGraph& source)
{
  const size_type requested = pathMergeInputBudget(parameters);
  const size_type pair_bytes = pathGraphFramedPairBytes(source);
  const size_type ceiling = pathMergeCeilingBudget(parameters);
  if(pair_bytes == 0)
  {
    // Raw shards hold no decoded blocks, but this budget also sizes the
    // merger's equal-label range and priority groups, which spill to disk once
    // they outgrow it. (Raw read windows are sized separately and do not grow
    // with it.) Capping it at --io-buffer-size (64 MiB) left each of step 1's
    // seven prune workers about 9 MiB: on the joint chr2+chr18 graph that
    // prune wrote 97.7 GiB of spill and output for 23.9 GB of output, at 1.23
    // cores, while the framed steps after it, funded from the same sixteenth
    // of --memory-limit as below, did not spill. Group sizes change only where
    // records are stored, never the merge order or the output.
    return std::max(requested, ceiling);
  }

  // --io-buffer-size sizes stream buffers; it was never a memory ceiling, and
  // no entry point above this library exposes it, so a generation whose blocks
  // were committed by an earlier configuration could not ask for more. Derive
  // the workspace from the global memory goal instead: reserve enough for one
  // legal unit of progress -- one decoded path block plus one decoded rank
  // block, which is the irreducible cost of opening a committed pair and
  // cannot be subdivided the way an oversized request normally is -- and then
  // grow toward holding every shard the merge may open, because a cache one
  // entry short of the shard count misses on every record. A sixteenth of
  // --memory-limit remains the bound, since the merger keeps a range deque, a
  // priority group and an output cache of the same size beside this one.
  const size_type pairs = std::max(static_cast<size_type>(1),
    std::min(pathMergeInputPairs(parameters), source.files()));
  size_type wanted = (pair_bytes > ceiling / pairs ? ceiling : pairs * pair_bytes);
  size_type budget = std::max(requested, std::min(wanted, ceiling));
  if(budget < pair_bytes)
  {
    // One pair still fits the ceiling: take it, so the merge can make progress
    // at reduced concurrency instead of refusing.
    budget = std::min(pair_bytes, ceiling);
  }
  return budget;
}

size_type
pathMergeInputCacheBudget(const ConstructionParameters& parameters,
  const PathGraph& source)
{
  const size_type pair_bytes = pathGraphFramedPairBytes(source);
  if(pair_bytes == 0) { return 0; }

  // Ask for every resident shard pair first. Anything less is not a smaller
  // cache but a broken one: the merge visits the shards round-robin, so one
  // entry short of the count misses on every record. The remainder of this
  // independent quarter-memory ceiling is a capacity, not an eager
  // allocation: PathGraphInputCache gives it to a shared block prefetch pool,
  // which uses only what its one-block-per-reader queue can keep productive.
  // If the pairs themselves consume the ceiling, prefetch simply receives no
  // bytes and the previous synchronous behavior remains intact.
  return std::max(pair_bytes, parameters.getMemoryLimitBytes() / 4);
}

size_type
mergeAdmissibleBlockSize(size_type byte_budget, size_type pairs,
  size_type requested_block)
{
  if(pairs == 0) { pairs = 1; }
  // Two decoded blocks and two decoder contexts per open pair.
  const size_type per_pair_share = byte_budget / (2 * pairs);
  // ConstructionParameters::setCompressionBlockSize() will not store anything
  // smaller, so returning less would be silently rounded back up.
  const size_type minimum_block = 64 * KILOBYTE;
  size_type block = std::max(minimum_block, requested_block);
  while(block > minimum_block &&
        CompressedBlockReader::workingMemoryEstimate(block) > per_pair_share)
  {
    block /= 2;
  }
  return std::max(minimum_block, std::min(block, requested_block));
}

namespace
{

size_type
pathSortLogicalBytes(size_type records, size_type record_bytes,
  const char* description)
{
  if(record_bytes != 0 &&
     records > std::numeric_limits<size_type>::max() / record_bytes)
  {
    externalSortFailure(std::string(description) + " overflow");
  }
  return records * record_bytes;
}

size_type
pathSortAddBytes(size_type left, size_type right, const char* description)
{
  if(left > std::numeric_limits<size_type>::max() - right)
  {
    externalSortFailure(std::string(description) + " overflow");
  }
  return left + right;
}

size_type
pathPairPeakBytes(size_type paths, size_type ranks,
  const TempFileCodecParameters& codec)
{
  const size_type path_bytes = pathSortLogicalBytes(paths, sizeof(PathNode),
    "path shard bytes");
  const size_type rank_bytes = pathSortLogicalBytes(ranks,
    sizeof(PathNode::rank_type), "rank shard bytes");
  if(!codec.enabled())
  {
    return pathSortAddBytes(path_bytes, rank_bytes, "path/rank shard bytes");
  }

  const std::uint64_t framed_paths =
    CompressedBlockWriter::maximumTemporaryBytes(path_bytes,
      sizeof(PathNode), codec.block_size);
  const std::uint64_t framed_ranks =
    CompressedBlockWriter::maximumTemporaryBytes(rank_bytes,
      (PathLabel::LABEL_LENGTH + 1) * sizeof(PathNode::rank_type),
      codec.block_size);
  if(framed_paths > std::numeric_limits<size_type>::max() ||
     framed_ranks > std::numeric_limits<size_type>::max() - framed_paths)
  {
    externalSortFailure("framed path/rank shard size overflow");
  }
  return static_cast<size_type>(framed_paths + framed_ranks);
}

size_type
pathSortStoredBytes(const std::string& name)
{
  struct stat info;
  if(::stat(name.c_str(), &info) != 0 || info.st_size < 0 ||
     static_cast<std::uintmax_t>(info.st_size) >
       std::numeric_limits<size_type>::max())
  {
    externalSortFailure("cannot determine stored shard bytes for " + name);
  }
  return static_cast<size_type>(info.st_size);
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

size_type
externalPathGraphShardPeakBytes(size_type paths, size_type ranks,
  size_type sort_byte_budget, const TempFileCodecParameters& requested_codec)
{
  if(sort_byte_budget < externalPathGraphSortMinimumBudget())
  {
    externalSortFailure("streaming sort byte budget is too small");
  }
  const size_type available = sort_byte_budget - PATH_SORT_FIXED_BYTES;
  const size_type materialize_available = available - pathSortReaderReservation();
  return pathPairPeakBytes(paths, ranks,
    boundedPathCodec(requested_codec, materialize_available));
}

struct ExternalPathSortSink::Impl
{
  Impl(PathGraph& target, size_type target_file, size_type byte_budget,
    size_type requested_fan_in, size_type size_limit,
    size_type& already_committed, ExternalPathSortStats* statistics,
    const TempFileCodecParameters& requested_codec) :
    graph(target), file(target_file), limit(size_limit),
    committed_bytes(already_committed), stats(statistics), fan_in(0),
    merge_records(0), run_records(0), path_count(0), rank_count(0),
    payload_bytes(0), stored_bytes(0), codec(), complete(false)
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
    this->codec = boundedPathCodec(requested_codec, materialize_available);
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
    const size_type record_bytes = sizeof(PathNode) +
      source.ranks() * sizeof(PathNode::rank_type);
    const size_type next_paths = this->path_count + 1;
    const size_type next_ranks = this->rank_count + source.ranks();
    const size_type peak_bytes = pathPairPeakBytes(next_paths, next_ranks,
      this->codec);
    if(peak_bytes > this->limit ||
       this->committed_bytes > this->limit - peak_bytes)
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
    std::uint64_t path_digest, rank_digest;
    std::string merged = this->finishRuns();
    const std::string& final_path = this->graph.path_names[this->file];
    const std::string& final_rank = this->graph.rank_names[this->file];
    std::string partial_path = final_path + ".partial";
    std::string partial_rank = final_rank + ".partial";
    std::remove(partial_path.c_str()); std::remove(partial_rank.c_str());
    writeSortedPathPair(merged, partial_path, partial_rank,
      this->path_count, this->rank_count, this->merge_records, this->stats,
      this->codec, path_digest, rank_digest);
    TempFile::remove(merged);
    syncPathSortFile(partial_path); syncPathSortFile(partial_rank);
    if(std::rename(partial_path.c_str(), final_path.c_str()) != 0 ||
       std::rename(partial_rank.c_str(), final_rank.c_str()) != 0)
    {
      externalSortFailure("cannot atomically install streaming sorted path pair");
    }

    this->graph.path_counts[this->file] = this->path_count;
    this->graph.path_checksums.at(this->file).record(final_path, path_digest);
    this->graph.rank_checksums.at(this->file).record(final_rank, rank_digest);
    this->graph.rank_counts[this->file] = this->rank_count;
    this->graph.path_count += this->path_count;
    this->graph.rank_count += this->rank_count;
    this->stored_bytes = pathSortAddBytes(pathSortStoredBytes(final_path),
      pathSortStoredBytes(final_rank), "stored path/rank shard bytes");
    if(this->stored_bytes > this->limit ||
       this->committed_bytes > this->limit - this->stored_bytes)
    {
      externalSortFailure("installed path/rank shard exceeds the configured disk limit");
    }
    this->committed_bytes += this->stored_bytes;
    this->complete = true;
  }

  PathGraph& graph;
  size_type file, limit;
  size_type& committed_bytes;
  ExternalPathSortStats* stats;
  size_type fan_in, merge_records, run_records;
  size_type path_count, rank_count, payload_bytes, stored_bytes;
  TempFileCodecParameters codec;
  bool complete;
  std::vector<PathSortRecord> records;
  std::vector<std::vector<std::string>> levels;
};

ExternalPathSortSink::ExternalPathSortSink(PathGraph& graph, size_type file,
  size_type byte_budget, size_type fan_in, size_type size_limit,
  size_type& committed_bytes, ExternalPathSortStats* stats,
  const TempFileCodecParameters& codec) :
  impl(new Impl(graph, file, byte_budget, fan_in, size_limit,
    committed_bytes, stats, codec))
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

size_type
ExternalPathSortSink::storedBytes() const
{
  return this->impl->stored_bytes;
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
  std::uint64_t path_digest, rank_digest;
  writeSortedPathPair(merged, partial_path, partial_rank, graph.path_counts[file],
    graph.rank_counts[file], merge_records, stats, TempFileCodecParameters(), path_digest, rank_digest);
  TempFile::remove(merged);
  syncPathSortFile(partial_path); syncPathSortFile(partial_rank);
  if(std::rename(partial_path.c_str(), final_path.c_str()) != 0 || std::rename(partial_rank.c_str(), final_rank.c_str()) != 0)
  {
    externalSortFailure("cannot atomically install sorted path pair");
  }
  std::string old_path = graph.path_names[file], old_rank = graph.rank_names[file];
  graph.path_names[file] = final_path; graph.rank_names[file] = final_rank;
  graph.path_checksums.at(file).record(final_path, path_digest);
  graph.rank_checksums.at(file).record(final_rank, rank_digest);
  // The rewritten pair is raw (writeSortedPathPair above passes a default,
  // disabled codec), so any total recorded for the previous shards is stale.
  // Fall back to the logical payload rather than carry a wrong charge.
  graph.stored_bytes = PathGraph::UNKNOWN;
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
  Spill files are scratch: nothing reads them after a crash, so they need no
  durability boundary. They are written back only so that the
  POSIX_FADV_DONTNEED which follows can drop their pages, which a cgroup
  charges like any other memory and which stay resident while dirty.
  fdatasync() achieved that, but an appended file also has a new size to
  commit, so every call waited for an ext4 journal commit: step 1 of the
  joint chr2+chr18 prune flushed every ~4.5 MiB and its blocked thread time
  was mostly journal commits or the dirty-page throttle. sync_file_range()
  writes back and waits for exactly the written range and leaves the
  metadata alone. A file system that does not support it gets fdatasync().
*/
static void
writeBackSpillRange(int file, off_t offset, size_type bytes, const char* failure)
{
  if(bytes == 0) { return; }
#if defined(SYNC_FILE_RANGE_WRITE)
  if(::sync_file_range(file, offset, static_cast<off_t>(bytes),
       SYNC_FILE_RANGE_WAIT_BEFORE | SYNC_FILE_RANGE_WRITE |
       SYNC_FILE_RANGE_WAIT_AFTER) == 0)
  {
    return;
  }
  if(errno != ENOSYS && errno != EINVAL && errno != ESPIPE)
  {
    throw std::runtime_error(failure);
  }
#else
  static_cast<void>(offset);
#endif
  if(::fdatasync(file) != 0) { throw std::runtime_error(failure); }
}

/*
  PathGraphMerger normally only needs a sliding window. An equal-label range
  can, however, be arbitrarily long before prune() decides how to emit it.
  Keep that range addressable on disk once its resident budget is exhausted.
*/
template<class Element>
struct SpillableGroup
{
  std::vector<Element> memory;
  /*
    Index of the first live element in `memory`. seek() used to erase the dead
    prefix, and because Element (PriorityNode, 96 bytes) is not trivially
    copyable, vector::erase moves element by element rather than issuing one
    memmove -- once per output range, of order 400 million times in MergedGraph
    alone. Release popped from a deque's front; a deque is unavailable here
    because flush() hands memory.data() to a bulk write. Advancing a head and
    compacting only when the dead prefix passes half the buffer makes it
    amortised O(1) per element.

    Invariant: `head` is zero whenever the group has spilled, because flush()
    clears the buffer and resets it. Only the resident branches add it.
  */
  size_type head;
  std::string          filename;
  int                  file;
  size_type            elements, disk_elements, offset, file_begin;
  size_type            byte_limit, write_bytes, read_bytes, advised_write;
  mutable size_type    cache_offset;
  mutable std::vector<Element> read_cache;
  size_type*           spill_counter;

  SpillableGroup(size_type limit, size_type* counter = nullptr) : head(0), file(-1), elements(0), disk_elements(0), offset(0), file_begin(0),
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
    if(!this->spilled() &&
       (this->liveSize() + 1) * sizeof(Element) > this->byte_limit)
    {
      this->filename = TempFile::getName("gcsa_prune_group");
      this->file = ::open(this->filename.c_str(), O_CREAT | O_TRUNC | O_RDWR, 0600);
      if(this->file < 0) { throw std::runtime_error("PathGraph::prune(): cannot create spill file"); }
#if defined(POSIX_FADV_SEQUENTIAL)
      static_cast<void>(::posix_fadvise(this->file, 0, 0, POSIX_FADV_SEQUENTIAL));
#endif
      if(this->spill_counter != nullptr) { (*this->spill_counter)++; }
      // The merger may have already discarded many small label ranges before
      // the first oversized range appears. Spill indexes are absolute merger
      // positions, while the new file starts at the current live window.
      // Establish that coordinate origin before flush(); otherwise get(i)
      // treats an absolute index as an index into the unflushed tail and can
      // return a corrupted PriorityNode (including a bogus source shard).
      this->file_begin = this->offset;
      this->disk_elements = this->offset;
      this->advised_write = this->offset;
      this->cache_offset = this->offset;
      this->flush();
    }
    this->memory.push_back(value);
    this->elements++;
    if(this->spilled() && this->liveSize() * sizeof(Element) >= this->write_bytes) { this->flush(); }
  }

  Element get(size_type i) const
  {
    if(i < this->offset || i >= this->elements) { throw std::out_of_range("PathGraph::prune(): spill index"); }
    if(!this->spilled()) { return this->memory[this->head + i - this->offset]; }
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
    if(!this->spilled()) { this->memory[this->head + i - this->offset] = value; return; }
    if(i >= this->disk_elements) { this->memory[i - this->disk_elements] = value; return; }
    this->flush(); this->write(&value, 1, i);
    writeBackSpillRange(this->file, this->fileOffset(i), bytesFor(1),
      "PathGraph::prune(): spill sync failed");
    this->discard(i, 1); sdsl::util::clear(this->read_cache);
  }

  size_type liveSize() const { return this->memory.size() - this->head; }

  void seek(size_type i)
  {
    if(i > this->elements) { throw std::out_of_range("PathGraph::prune(): spill seek"); }
    if(!this->spilled())
    {
      this->head += (i - this->offset);
      if(this->head > this->memory.size() / 2)
      {
        this->memory.erase(this->memory.begin(),
          this->memory.begin() + this->head);
        this->head = 0;
      }
    }
    this->offset = i;
  }

  void clear()
  {
    if(this->file >= 0) { ::close(this->file); this->file = -1; }
    if(!this->filename.empty()) { TempFile::remove(this->filename); this->filename.clear(); }
    sdsl::util::clear(this->memory); this->head = 0;
    sdsl::util::clear(this->read_cache);
    this->elements = 0; this->disk_elements = 0; this->offset = 0; this->file_begin = 0;
    this->advised_write = 0; this->cache_offset = 0;
  }

private:
  void flush()
  {
    if(this->liveSize() == 0)
    {
      sdsl::util::clear(this->memory); this->head = 0; return;
    }
    const size_type start = this->disk_elements, count = this->liveSize();
    this->write(this->memory.data() + this->head, count, start);
    this->disk_elements += count;
    sdsl::util::clear(this->memory); this->head = 0;
    this->memory.reserve(std::max(static_cast<size_type>(1), this->write_bytes / sizeof(Element)));
    writeBackSpillRange(this->file, this->fileOffset(start), bytesFor(count),
      "PathGraph::prune(): spill sync failed");
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

  off_t fileOffset(size_type start) const
  {
    if(start < this->file_begin)
    {
      throw std::out_of_range("PathGraph::prune(): spill offset before live window");
    }
    size_type offset = bytesFor(start - this->file_begin);
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
    writeBackSpillRange(this->file, fileOffset(start), bytesFor(this->disk_end - start),
      "PathGraphMerger: range spill sync failed");
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
  Keep a fixed number of path/rank descriptor pairs. Framed readers share one
  bounded prefetch/decode pool; raw readers retain the synchronous windowed
  path. The heap still owns one compact PriorityNode per shard because exact
  k-way merge ordering requires a head.
*/
struct PathGraphInputCache
{
  struct Window
  {
    std::vector<char> data;
    off_t start;
    size_type valid;

    explicit Window(size_type bytes = 0) : data(bytes), start(0), valid(0) { }

    void reset() { this->start = 0; this->valid = 0; }

    bool contains(off_t offset, size_type bytes) const
    {
      if(offset < this->start) { return false; }
      size_type relative = static_cast<size_type>(offset - this->start);
      return (relative <= this->valid && bytes <= this->valid - relative);
    }
  };

  struct Entry
  {
    size_type file, stamp;
    int path, rank;
    off_t path_high, rank_high, path_released, rank_released;
    Window path_window, rank_window;
    std::unique_ptr<CompressedBlockReader> compressed_path, compressed_rank;

    explicit Entry(size_type window_bytes = 0) :
      file(PathGraph::UNKNOWN), stamp(0), path(-1), rank(-1),
      path_high(0), rank_high(0), path_released(0), rank_released(0),
      path_window(window_bytes), rank_window(window_bytes),
      compressed_path(), compressed_rank() { }

    Entry(Entry&&) = default;
    Entry& operator=(Entry&&) = default;

    Entry(const Entry&) = delete;
    Entry& operator=(const Entry&) = delete;

    void reset(size_type new_file, size_type new_stamp)
    {
      this->file = new_file; this->stamp = new_stamp;
      this->path = -1; this->rank = -1;
      this->path_high = 0; this->rank_high = 0;
      this->path_released = 0; this->rank_released = 0;
      this->path_window.reset(); this->rank_window.reset();
      this->compressed_path.reset(); this->compressed_rank.reset();
    }
  };

  const PathGraph& graph;
  std::vector<Entry> entries;
  // Logical shard -> index into `entries`, or NO_ENTRY. get() ran a linear scan
  // of the open shards on every record; with ~30 resident pairs and 1.2e9 paths
  // through the merge that is of order 3e10 comparisons, and it grows with the
  // shard count, so it degrades exactly as the input gets larger. Indices rather
  // than pointers because emplace_back may reallocate `entries`.
  std::vector<size_type> file_to_entry;
  size_type clock, max_pairs, window_bytes, compressed_pair_bytes;
  size_type minimum_prefetch_job_bytes;
  size_type oversized_pair_bytes;
  bool transient_descriptors;
  PathGraphMergeStats* stats;
  std::shared_ptr<CompressedBlockPrefetchPool> prefetch_pool;

  constexpr static size_type NO_ENTRY = ~static_cast<size_type>(0);

  constexpr static off_t CACHE_TAIL = 64 * KILOBYTE;
  constexpr static size_type MIN_WINDOW = 4 * KILOBYTE;
  constexpr static size_type MAX_WINDOW = 64 * KILOBYTE;
  constexpr static size_type MAX_PREFETCH_WORKERS = 8;

  PathGraphInputCache(const PathGraph& source, PathGraphMergeStats* merge_stats,
    size_type requested_pairs, size_type byte_budget,
    size_type max_prefetch_workers = MAX_PREFETCH_WORKERS) :
    graph(source), entries(), file_to_entry(source.files(), NO_ENTRY), clock(0),
    max_pairs(std::max(static_cast<size_type>(1),
      std::min(requested_pairs, std::max(static_cast<size_type>(1), source.files())))),
    window_bytes(0), compressed_pair_bytes(0),
    minimum_prefetch_job_bytes(PathGraph::UNKNOWN), oversized_pair_bytes(0),
    transient_descriptors(false),
    stats(merge_stats), prefetch_pool()
  {
    // A framed reader retains one decoded block. Bound the LRU cardinality by
    // bytes as well as descriptors, so a large --max-open-files value cannot
    // multiply the compression block into an unaccounted resident peak.
    bool all_framed = true;
    for(size_type file = 0; file < source.files(); file++)
    {
      const bool path_framed = CompressedBlockReader::isFramed(source.path_names[file]);
      const bool rank_framed = CompressedBlockReader::isFramed(source.rank_names[file]);
      if(path_framed != rank_framed)
      {
        throw std::runtime_error(
          "PathGraphMerger: path/rank shard storage formats differ");
      }
      if(!path_framed) { all_framed = false; }
      if(path_framed)
      {
        size_type path_block = CompressedBlockReader::declaredBlockSize(
          source.path_names[file]);
        size_type rank_block = CompressedBlockReader::declaredBlockSize(
          source.rank_names[file]);
        size_type pair_bytes = CompressedBlockReader::workingMemoryEstimate(
          path_block) + CompressedBlockReader::workingMemoryEstimate(rank_block);
        this->compressed_pair_bytes = std::max(this->compressed_pair_bytes,
          pair_bytes);
        this->minimum_prefetch_job_bytes = std::min(
          this->minimum_prefetch_job_bytes,
          std::min(CompressedBlockReader::prefetchWorkingMemoryEstimate(path_block),
            CompressedBlockReader::prefetchWorkingMemoryEstimate(rank_block)));
      }
    }
    if(this->compressed_pair_bytes > 0)
    {
      // Parallel mergers must retain decoded blocks from all active shards.
      // A tiny descriptor LRU otherwise repeatedly decompresses the same
      // blocks while the merge heap alternates between shards. Transient
      // synchronous readers separate decoded-cache capacity from open FDs.
      this->transient_descriptors = all_framed && max_prefetch_workers == 0;
      if(this->transient_descriptors)
      {
        this->max_pairs = std::max(static_cast<size_type>(1), source.files());
      }
      // One open pair is the irreducible minimum for a merge, and a committed
      // block size cannot be renegotiated: the shards on disk declare it, and
      // no operational flag rewrites them. Refusing here made a workspace
      // written with a larger block permanently unresumable, so admit the
      // single pair and report the overshoot instead of aborting. Fresh shards
      // do not reach this branch, because mergeAdmissibleBlockSize() bounds the
      // block the writer commits by what this cache can hold.
      // PathGraphMerger resets *stats after this cache is constructed, so the
      // value is kept here and copied across once that reset has happened.
      if(this->compressed_pair_bytes > byte_budget)
      {
        this->oversized_pair_bytes = this->compressed_pair_bytes;
      }
      this->max_pairs = std::min(this->max_pairs,
        std::max(static_cast<size_type>(1),
          byte_budget / this->compressed_pair_bytes));
      // A generation does not normally mix raw and framed shards. Avoid
      // allocating raw pread windows in every framed LRU entry.
      this->window_bytes = 0;
      this->entries.reserve(this->max_pairs);
      size_type resident_bytes = byte_budget;
      if(this->compressed_pair_bytes <= byte_budget / this->max_pairs)
      {
        resident_bytes = this->compressed_pair_bytes * this->max_pairs;
      }
      size_type prefetch_bytes = (byte_budget > resident_bytes ?
        byte_budget - resident_bytes : 0);
      size_type available_threads = static_cast<size_type>(omp_get_max_threads());
      size_type requested_workers = (available_threads > 1 ?
        std::min(max_prefetch_workers, available_threads - 1) : 0);
      if(requested_workers > 0 &&
         this->minimum_prefetch_job_bytes != PathGraph::UNKNOWN &&
         prefetch_bytes >= this->minimum_prefetch_job_bytes)
      {
        size_type workers = std::min(requested_workers,
          prefetch_bytes / this->minimum_prefetch_job_bytes);
        this->prefetch_pool.reset(new CompressedBlockPrefetchPool(
          prefetch_bytes, workers));
      }
      return;
    }
    // At most 2 * max_pairs windows can coexist. Reserve at most one quarter
    // of the caller's merge workspace for them, leaving the rest for equal-
    // label groups, range metadata, and phase-specific state.
    size_type candidate = byte_budget / this->max_pairs / 8;
    if(candidate >= MIN_WINDOW)
    {
      this->window_bytes = std::min(MAX_WINDOW,
        candidate - (candidate % MIN_WINDOW));
    }
    this->entries.reserve(this->max_pairs);
  }

  ~PathGraphInputCache() { this->close(); }

  // When every shard pair fits, open all pooled readers before consuming the
  // first heap head. Their deferred block-0 jobs then form large, independent
  // preads instead of a constructor-time sequence of read/decode stalls.
  // A cache shared by consecutive mergers reports into the current one.
  void setStats(PathGraphMergeStats* merge_stats) { this->stats = merge_stats; }

  void prime()
  {
    if(!(this->prefetch_pool) || this->max_pairs < this->graph.files()) { return; }
    if(this->stats != nullptr)
    {
      this->stats->prefetch_workers = this->prefetch_pool->workers();
    }
    for(size_type file = 0; file < this->graph.files(); ++file)
    {
      static_cast<void>(this->get(file));
    }
  }

  void read(size_type file, size_type offset, PathNode& node,
    PathNode::rank_type* labels)
  {
    Entry& entry = this->get(file);
    off_t path_offset = this->checkedOffset(offset, sizeof(PathNode));
    off_t path_limit = this->checkedOffset(this->graph.path_counts[file],
      sizeof(PathNode));
    this->readWindow(entry.path, entry.compressed_path.get(),
      entry.path_window, &node, sizeof(node),
      path_offset, path_limit, true);
    if(this->stats != nullptr) { this->stats->path_input_reads++; }
    entry.path_high = std::max(entry.path_high,
      path_offset + static_cast<off_t>(sizeof(node)));
    if(entry.path >= 0) { this->trim(entry.path, entry.path_high, entry.path_released); }

    if(node.ranks() > PathLabel::LABEL_LENGTH + 1)
    {
      throw std::runtime_error("PathGraphMerger: invalid path rank count");
    }
    size_type rank_bytes = node.ranks() * sizeof(PathNode::rank_type);
    if(node.pointer() > this->graph.rank_counts[file] ||
       node.ranks() > this->graph.rank_counts[file] - node.pointer())
    {
      throw std::runtime_error("PathGraphMerger: path rank range is outside its sidecar");
    }
    off_t rank_offset = this->checkedOffset(node.pointer(), sizeof(PathNode::rank_type));
    off_t rank_limit = this->checkedOffset(this->graph.rank_counts[file],
      sizeof(PathNode::rank_type));
    this->readWindow(entry.rank, entry.compressed_rank.get(),
      entry.rank_window, labels, rank_bytes,
      rank_offset, rank_limit, false);
    if(this->stats != nullptr) { this->stats->rank_input_reads++; }
    entry.rank_high = std::max(entry.rank_high,
      rank_offset + static_cast<off_t>(rank_bytes));
    if(entry.rank >= 0) { this->trim(entry.rank, entry.rank_high, entry.rank_released); }
  }

  void close()
  {
    for(Entry& entry : this->entries) { this->close(entry); }
    this->refreshPrefetchStats();
    this->entries.clear();
    // The index map must die with the entries it points into.
    std::fill(this->file_to_entry.begin(), this->file_to_entry.end(), NO_ENTRY);
    this->prefetch_pool.reset();
  }

private:
  Entry& get(size_type file)
  {
    this->clock++;
    if(file < this->file_to_entry.size() &&
       this->file_to_entry[file] != NO_ENTRY)
    {
      Entry& entry = this->entries[this->file_to_entry[file]];
      entry.stamp = this->clock;
      return entry;
    }

    Entry* target = nullptr;
    size_type target_index = 0;
    if(this->entries.size() < this->max_pairs)
    {
      target_index = this->entries.size();
      this->entries.emplace_back(this->window_bytes); target = &(this->entries.back());
      if(this->stats != nullptr)
      {
        this->stats->max_open_input_pairs = std::max(
          this->stats->max_open_input_pairs,
          (this->transient_descriptors ? static_cast<size_type>(1) :
            static_cast<size_type>(this->entries.size())));
        this->stats->max_input_buffer_bytes = std::max(
          this->stats->max_input_buffer_bytes,
          (this->compressed_pair_bytes > 0 ? this->compressed_pair_bytes :
            2 * this->window_bytes) * static_cast<size_type>(this->entries.size()));
      }
    }
    else
    {
      auto victim = std::min_element(this->entries.begin(), this->entries.end(),
        [](const Entry& left, const Entry& right) { return left.stamp < right.stamp; });
      target_index = static_cast<size_type>(victim - this->entries.begin());
      target = &(*victim);
      // Drop the evicted shard's mapping before the entry is repurposed.
      if(target->file != PathGraph::UNKNOWN &&
         target->file < this->file_to_entry.size())
      {
        this->file_to_entry[target->file] = NO_ENTRY;
      }
      this->close(*target);
    }

    target->reset(file, this->clock);
    if(file < this->file_to_entry.size()) { this->file_to_entry[file] = target_index; }
    const bool path_framed = CompressedBlockReader::isFramed(
      this->graph.path_names[file]);
    const bool rank_framed = CompressedBlockReader::isFramed(
      this->graph.rank_names[file]);
    if(path_framed != rank_framed)
    {
      throw std::runtime_error(
        "PathGraphMerger: path/rank shard storage formats differ");
    }
    if(path_framed)
    {
      if(this->prefetch_pool)
      {
        target->compressed_path.reset(new CompressedBlockReader(
          this->graph.path_names[file], this->prefetch_pool));
        target->compressed_rank.reset(new CompressedBlockReader(
          this->graph.rank_names[file], this->prefetch_pool));
      }
      else
      {
        CompressedBlockReader::FileAccess access = (this->transient_descriptors ?
          CompressedBlockReader::FileAccess::TRANSIENT :
          CompressedBlockReader::FileAccess::PERSISTENT);
        target->compressed_path.reset(new CompressedBlockReader(
          this->graph.path_names[file], access));
        target->compressed_rank.reset(new CompressedBlockReader(
          this->graph.rank_names[file], access));
      }
      if(target->compressed_path->logicalSize() !=
           this->graph.path_counts[file] * sizeof(PathNode) ||
         target->compressed_rank->logicalSize() !=
           this->graph.rank_counts[file] * sizeof(PathNode::rank_type))
      {
        throw std::runtime_error(
          "PathGraphMerger: compressed shard length does not match metadata");
      }
      return *target;
    }

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
    entry.compressed_path.reset(); entry.compressed_rank.reset();
  }

  void refreshPrefetchStats()
  {
    if(this->stats == nullptr || !(this->prefetch_pool)) { return; }
    CompressedBlockPrefetchPool::Stats pool = this->prefetch_pool->stats();
    this->stats->prefetch_workers = pool.worker_threads;
    this->stats->prefetch_submitted = pool.submitted;
    this->stats->prefetch_completed = pool.completed;
    this->stats->prefetch_consumed = pool.consumed;
    this->stats->prefetch_ready_hits = pool.ready_hits;
    this->stats->prefetch_waits = pool.waits;
    this->stats->prefetch_synchronous_blocks = pool.synchronous_blocks;
    this->stats->prefetch_cancelled = pool.cancelled;
    this->stats->prefetch_errors = pool.errors;
    this->stats->prefetch_physical_bytes = pool.physical_bytes;
    this->stats->prefetch_decoded_bytes = pool.decoded_bytes;
    this->stats->prefetch_wait_nanoseconds = pool.wait_nanoseconds;
    this->stats->max_prefetch_bytes = pool.peak_bytes;
    size_type resident = this->compressed_pair_bytes;
    if(resident > 0 && this->stats->max_open_input_pairs <=
       std::numeric_limits<size_type>::max() / resident)
    {
      resident *= this->stats->max_open_input_pairs;
      if(pool.peak_bytes <= std::numeric_limits<size_type>::max() - resident)
      {
        this->stats->max_input_buffer_bytes = std::max(
          this->stats->max_input_buffer_bytes, resident + pool.peak_bytes);
      }
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

  void readWindow(int descriptor, CompressedBlockReader* compressed,
    Window& window, void* target, size_type bytes, off_t offset, off_t limit,
    bool path_stream)
  {
    if(offset < 0 || limit < offset ||
       bytes > static_cast<size_type>(limit - offset))
    {
      throw std::runtime_error("PathGraphMerger: input range is outside its file");
    }
    if(bytes == 0) { return; }

    if(compressed != nullptr)
    {
      if(compressed->readAt(static_cast<std::uint64_t>(offset), target, bytes) != bytes)
      {
        throw std::runtime_error("PathGraphMerger: truncated compressed input");
      }
      if(this->stats != nullptr) { this->stats->direct_input_reads++; }
      return;
    }

    if(window.data.empty() || bytes > window.data.size())
    {
      this->preadAll(descriptor, target, bytes, offset);
      if(this->stats != nullptr) { this->stats->direct_input_reads++; }
      return;
    }
    if(!(window.contains(offset, bytes)))
    {
      size_type available = static_cast<size_type>(limit - offset);
      size_type refill = std::min(static_cast<size_type>(window.data.size()), available);
      this->preadAll(descriptor, window.data.data(), refill, offset);
      window.start = offset; window.valid = refill;
      if(this->stats != nullptr)
      {
        if(path_stream) { this->stats->path_input_refills++; }
        else { this->stats->rank_input_refills++; }
      }
    }
    size_type relative = static_cast<size_type>(offset - window.start);
    std::memcpy(target, window.data.data() + relative, bytes);
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

constexpr off_t PathGraphInputCache::CACHE_TAIL;
constexpr size_type PathGraphInputCache::MIN_WINDOW;
constexpr size_type PathGraphInputCache::MAX_WINDOW;
constexpr size_type PathGraphInputCache::MAX_PREFETCH_WORKERS;

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
  // A merger owns its cache unless the caller lends one that outlives it, so
  // consecutive partitions on one worker reuse decoded blocks.
  std::unique_ptr<PathGraphInputCache>          owned_input;
  PathGraphInputCache&                          input_files;
  std::vector<size_type>                        offsets, end_offsets;
  size_type                                     path_count;
  PriorityQueue<PriorityNode>                   inputs;

  PathGraphMerger(const PathGraph& path_graph, const LCP& kmer_lcp,
    size_type group_buffer_bytes = MEGABYTE,
    PathGraphMergeStats* stats = nullptr,
    size_type max_input_pairs = 32,
    size_type input_cache_bytes = 0,
    PathNode::rank_type lower_rank = 0,
    PathNode::rank_type upper_rank = PathLabel::NO_RANK,
    size_type max_prefetch_workers = PathGraphInputCache::MAX_PREFETCH_WORKERS,
    const size_type* start_offsets = nullptr, const size_type* stop_offsets = nullptr,
    PathGraphInputCache* shared_input = nullptr);
  void close();

  inline size_type size() const { return this->path_count; }

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
  size_type group_buffer_bytes, PathGraphMergeStats* stats,
  size_type max_input_pairs, size_type input_cache_bytes,
  PathNode::rank_type lower_rank, PathNode::rank_type upper_rank,
  size_type max_prefetch_workers,
  const size_type* start_offsets, const size_type* stop_offsets,
  PathGraphInputCache* shared_input) :
  graph(path_graph), lcp(kmer_lcp),
  ranges(group_buffer_bytes, (stats == nullptr ? nullptr : &(stats->range_spills))),
  buffer(group_buffer_bytes, (stats == nullptr ? nullptr : &(stats->priority_spills))),
  owned_input(shared_input != nullptr ? nullptr : new PathGraphInputCache(path_graph, stats,
    max_input_pairs, (input_cache_bytes == 0 ? group_buffer_bytes : input_cache_bytes),
    max_prefetch_workers)),
  input_files(shared_input != nullptr ? *shared_input : *owned_input),
  offsets(path_graph.files()), end_offsets(path_graph.files()), path_count(0),
  inputs(path_graph.files())
{
  if(stats != nullptr)
  {
    *stats = PathGraphMergeStats();
    // The input cache is a member, so it was constructed before this reset.
    stats->oversized_input_pair_bytes = this->input_files.oversized_pair_bytes;
  }
  if(shared_input != nullptr) { shared_input->setStats(stats); }
  this->input_files.prime();
  auto lower_bound = [&](size_type file, PathNode::rank_type key) -> size_type
  {
    size_type low = 0, high = path_graph.path_counts[file];
    PriorityNode probe; probe.file = file;
    while(low < high)
    {
      size_type mid = low + (high - low) / 2;
      this->input_files.read(file, mid, probe.node, probe.label);
      if(probe.firstLabel(0) < key) { low = mid + 1; }
      else { high = mid; }
    }
    return low;
  };
  for(size_type file = 0; file < path_graph.files(); file++)
  {
    // A caller that located the bounds for many partitions at once passes
    // them in; searching again here costs block decodes per partition.
    this->offsets[file] = (start_offsets != nullptr ? start_offsets[file] :
      (lower_rank == 0 ? 0 : lower_bound(file, lower_rank)));
    this->end_offsets[file] = (stop_offsets != nullptr ? stop_offsets[file] :
      (upper_rank == PathLabel::NO_RANK ?
        path_graph.path_counts[file] : lower_bound(file, upper_rank)));
    if(this->end_offsets[file] < this->offsets[file])
    {
      throw std::runtime_error("PathGraphMerger: invalid label-rank partition");
    }
    if(this->end_offsets[file] - this->offsets[file] >
       std::numeric_limits<size_type>::max() - this->path_count)
    {
      throw std::overflow_error("PathGraphMerger: partition path count overflow");
    }
    this->path_count += this->end_offsets[file] - this->offsets[file];
    this->inputs[file].file = file; this->read(this->inputs[file]);
  }
  this->inputs.heapify();
}

void
PathGraphMerger::close()
{
  this->ranges.clear();
  this->buffer.clear();
  if(this->owned_input) { this->input_files.close(); }
  this->offsets.clear();
  this->end_offsets.clear(); this->path_count = 0;
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
  if(this->size() == 0) { return range_type(0, 0); }

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
  if(this->offsets[path.file] >= this->end_offsets[path.file])
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
  this->path_checksums.resize(source.files());
  this->rank_checksums.resize(source.files());
  this->path_count = 0; this->rank_count = 0; this->range_count = 0;
  this->order = source.k(); this->doubling_steps = 0;
  this->unique = UNKNOWN; this->redundant = UNKNOWN;
  this->unsorted = UNKNOWN; this->nondeterministic = UNKNOWN;
  this->stored_bytes = UNKNOWN;
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
  path_checksums(file_count), rank_checksums(file_count),
  logical_file_ids(), physical_shard_ids(),
  path_count(0), rank_count(0), range_count(0), order(path_order), doubling_steps(steps),
  unique(0), redundant(0), unsorted(0), nondeterministic(0),
  stored_bytes(UNKNOWN),
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
  this->path_checksums.resize(1); this->rank_checksums.resize(1);
  this->path_count = 0; this->rank_count = 0; this->range_count = 0;
  this->order = 0; this->doubling_steps = 0;
  this->unique = 0; this->redundant = 0;
  this->unsorted = 0; this->nondeterministic = 0;
  this->stored_bytes = UNKNOWN;
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
  this->path_checksums.clear();
  this->rank_checksums.clear();
  this->logical_file_ids.clear();
  this->physical_shard_ids.clear();

  this->path_count = 0; this->rank_count = 0;
  this->order = 0;
  this->unique = UNKNOWN; this->redundant = UNKNOWN;
  this->unsorted = UNKNOWN; this->nondeterministic = UNKNOWN;
  this->stored_bytes = UNKNOWN;
}

void
PathGraph::swap(PathGraph& another) noexcept
{
  this->path_names.swap(another.path_names);
  this->rank_names.swap(another.rank_names);
  this->path_counts.swap(another.path_counts);
  this->rank_counts.swap(another.rank_counts);
  this->path_checksums.swap(another.path_checksums);
  this->rank_checksums.swap(another.rank_checksums);
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
  std::swap(this->stored_bytes, another.stored_bytes);
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

namespace
{

struct PruneOutputLayout
{
  std::vector<size_type> source_to_output;
  std::vector<logical_file_id_t> logical;
  std::vector<physical_shard_id_t> physical;
};

PruneOutputLayout
pruneOutputLayout(const PathGraph& source)
{
  PruneOutputLayout result;
  result.source_to_output.resize(source.files(), PathGraph::UNKNOWN);
  std::map<logical_file_id_t, size_type> outputs;
  for(size_type file = 0; file < source.files(); file++)
  {
    logical_file_id_t logical = source.logicalFile(file);
    auto inserted = outputs.emplace(logical, result.logical.size());
    if(inserted.second)
    {
      result.logical.push_back(logical);
      // The serial route preserves the first physical identity. Parallel
      // publication normalizes all retained shards after partition assembly.
      result.physical.push_back(source.physicalShard(file));
    }
    result.source_to_output[file] = inserted.first->second;
  }
  return result;
}

struct PruneSizeLimit { };
struct PruneCancelled { };

bool
reservePruneOutput(std::atomic<size_type>& used, size_type bytes,
  size_type limit)
{
  size_type current = used.load(std::memory_order_relaxed);
  while(true)
  {
    if(current > limit || bytes > limit - current) { return false; }
    if(used.compare_exchange_weak(current, current + bytes,
      std::memory_order_relaxed, std::memory_order_relaxed))
    {
      return true;
    }
  }
}

/*
  The only safe independent prune ranges are separated by a zero k-mer LCP.
  Use wavelet-tree select rather than scanning all keys: this setup cost is
  proportional to the number of root components, not the graph or key count.
*/
bool
pruneRootPartitions(const LCP& lcp,
  std::vector<std::pair<PathNode::rank_type, PathNode::rank_type>>& ranges,
  size_type depth)
{
  ranges.clear();
  if(lcp.total_keys < 2 || lcp.total_keys >= PathLabel::NO_RANK ||
     lcp.kmer_lcp.size() != lcp.total_keys)
  {
    return false;
  }

  // Split wherever adjacent keys share fewer than depth leading characters.
  std::vector<size_type> boundaries; boundaries.push_back(0);
  for(size_type value = 0; value < depth; value++)
  {
    const LCP::rank_type shared = value;
    size_type count = lcp.kmer_lcp.rank(lcp.total_keys, shared);
    for(size_type occurrence = 1; occurrence <= count; occurrence++)
    {
      size_type boundary = lcp.kmer_lcp.select(occurrence, shared);
      if(boundary > 0 && boundary < lcp.total_keys) { boundaries.push_back(boundary); }
    }
  }
  boundaries.push_back(lcp.total_keys);
  std::sort(boundaries.begin(), boundaries.end());
  boundaries.erase(std::unique(boundaries.begin(), boundaries.end()),
    boundaries.end());
  if(boundaries.size() < 3) { return false; }

  for(size_type i = 1; i < boundaries.size(); i++)
  {
    if(boundaries[i] <= boundaries[i - 1] ||
       (i + 1 < boundaries.size() && lcp.kmer_lcp[boundaries[i]] >= depth))
    {
      ranges.clear(); return false;
    }
    ranges.push_back(std::make_pair(
      static_cast<PathNode::rank_type>(boundaries[i - 1]),
      static_cast<PathNode::rank_type>(boundaries[i])));
  }
  return (ranges.size() > 1);
}

// Experimental: split prune partitions where adjacent keys share fewer than
// this many leading characters. 1, the default, splits only where the first
// character changes, the zero-LCP boundaries no group can cross.
size_type
pruneSplitDepth()
{
  const char* value = std::getenv("GCSA_EXPERIMENTAL_PRUNE_SPLIT_DEPTH");
  if(value == nullptr || *value == '\0') { return 1; }
  const unsigned long long depth = std::strtoull(value, nullptr, 10);
  return (depth == 0 ? 1 : static_cast<size_type>(std::min<unsigned long long>(depth, 16)));
}

/*
  Runs visit(file, input) once for every shard, on up to `threads` threads that
  each read through their own single-pair input cache, and rethrows the first
  error. The shard passes before a parallel prune or merge are independent per
  shard, and each decodes whole framed blocks, so they belong on every thread
  the memory allows rather than on the workers alone.
*/
template<class Visit>
void
forEachShard(const PathGraph& source, size_type pair_buffer_bytes,
  size_type threads, Visit visit)
{
  const size_type files = source.files();
  if(files == 0) { return; }
  threads = std::max(static_cast<size_type>(1), std::min(threads, files));
  std::atomic<size_type> next_file(0);
  std::vector<std::exception_ptr> errors(threads);
  auto run = [&](size_type thread)
  {
    try
    {
      PathGraphInputCache input(source, nullptr, 1, pair_buffer_bytes, 0);
      for(size_type file = next_file++; file < files; file = next_file++)
      {
        visit(file, input);
      }
      input.close();
    }
    catch(...) { errors[thread] = std::current_exception(); }
  };
  std::vector<std::thread> pool;
  for(size_type thread = 1; thread < threads; thread++)
  {
    try { pool.emplace_back(run, thread); }
    catch(const std::system_error&) { break; } // The started threads and this one finish the queue.
  }
  run(0);
  for(std::thread& worker : pool) { worker.join(); }
  for(const std::exception_ptr& error : errors) { if(error) { std::rethrow_exception(error); } }
}

// Threads for the shard passes: every OpenMP thread, as many as shards, and
// no more single-pair caches than the input cache budget holds at once.
size_type
shardPassThreads(const PathGraph& source, size_type pair_bytes, size_type cache_budget)
{
  size_type threads = static_cast<size_type>(std::max(1, omp_get_max_threads()));
  threads = std::min(threads, std::max(static_cast<size_type>(1), source.files()));
  if(pair_bytes > 0) { threads = std::min(threads, std::max(static_cast<size_type>(1), cache_budget / pair_bytes)); }
  return threads;
}

// Merge inputs are sorted within each physical shard, so the last record has
// its largest first rank. Read exactly that record through the normal random-
// access path before workers create output; otherwise lower_bound(total_keys)
// could silently exclude a corrupt out-of-domain tail from every partition.
bool
shardTailsInDomain(const PathGraph& source, const LCP& lcp,
  size_type pair_buffer_bytes, size_type threads)
{
  std::atomic<bool> in_domain(true);
  forEachShard(source, pair_buffer_bytes, threads,
    [&](size_type file, PathGraphInputCache& input)
    {
      if(source.path_counts[file] == 0 || !in_domain.load(std::memory_order_relaxed)) { return; }
      PriorityNode tail; tail.file = file;
      input.read(file, source.path_counts[file] - 1, tail.node, tail.label);
      if(tail.firstLabel(0) >= lcp.total_keys) { in_domain.store(false, std::memory_order_relaxed); }
    });
  return in_domain.load();
}

void
validatePruneShardTails(const PathGraph& source, const LCP& lcp,
  size_type pair_buffer_bytes, size_type threads)
{
  if(!shardTailsInDomain(source, lcp, pair_buffer_bytes, threads))
  {
    throw std::runtime_error(
      "PathGraph::prune(): input label rank is outside the LCP key domain");
  }
}

template<bool report_progress>
void
prunePathRange(const PathGraph& source, const LCP& lcp,
  PathNode::rank_type lower_rank, PathNode::rank_type upper_rank,
  size_type group_buffer_bytes, size_type input_cache_bytes,
  size_type input_pairs, const PruneOutputLayout& layout,
  PathGraphBuilder& builder, PathGraphMergeStats* stats,
  std::atomic<size_type>* output_bytes, size_type size_limit,
  std::atomic<bool>* cancelled, PriorityNode* first_output,
  PriorityNode* last_output, bool* has_output, ProgressReporter* progress,
  bool* first_range_open = nullptr,
  const size_type* start_offsets = nullptr, const size_type* stop_offsets = nullptr,
  PathGraphInputCache* shared_input = nullptr)
{
  PathGraphMerger merger(source, lcp, group_buffer_bytes, stats, input_pairs,
    input_cache_bytes, lower_rank, upper_rank,
    // Worker-level merge concurrency replaces the decoder pool. Keeping both
    // would oversubscribe CPUs and duplicate cache reservations.
    (cancelled == nullptr ? PathGraphInputCache::MAX_PREFETCH_WORKERS : 0),
    start_offsets, stop_offsets, shared_input);

  auto write_output = [&](PriorityNode node)
  {
    if(cancelled != nullptr && cancelled->load(std::memory_order_relaxed))
    {
      throw PruneCancelled();
    }
    if(node.file >= layout.source_to_output.size() ||
       layout.source_to_output[node.file] == PathGraph::UNKNOWN)
    {
      throw std::runtime_error("PathGraph::prune(): invalid source shard identity");
    }
    if(output_bytes != nullptr &&
       !reservePruneOutput(*output_bytes, node.bytes(), size_limit))
    {
      if(cancelled != nullptr) { cancelled->store(true, std::memory_order_relaxed); }
      throw PruneSizeLimit();
    }
    if(has_output != nullptr)
    {
      if(!(*has_output) && first_output != nullptr) { *first_output = node; }
      if(last_output != nullptr) { *last_output = node; }
      *has_output = true;
    }
    node.file = layout.source_to_output[node.file];
    builder.write(node);
  };

  for(range_type range = merger.first(); !(merger.atEnd(range));
      range = merger.next())
  {
    if(cancelled != nullptr && cancelled->load(std::memory_order_relaxed))
    {
      throw PruneCancelled();
    }
    SameFromLogicalFile same_from(merger, range);
    if(same_from.same_from)
    {
      if(same_from.same_file)
      {
        range = merger.extendRange(same_from);
        // Only a first range that merged its whole partition could have
        // continued into the next one; the caller decides whether it might.
        if(first_range_open != nullptr && range.first == 0 &&
           range.second + 1 == merger.size())
        {
          *first_range_open = true;
        }
        merger.mergePathNodes();
        write_output(merger.buffer.get(range.second));
        builder.graph.unique++;
      }
      else
      {
        for(size_type i = range.first; i <= range.second; i++)
        {
          PriorityNode node = merger.buffer.get(i);
          node.node.makeSorted(); write_output(node);
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
        write_output(node);
      }
    }
    builder.graph.range_count++;
    if(report_progress && progress != nullptr)
    {
      progress->advance(range.second - range.first + 1);
    }
  }
  merger.close();
}

void
addPruneStats(PathGraphMergeStats& total, const PathGraphMergeStats& part)
{
  total.priority_spills += part.priority_spills;
  total.range_spills += part.range_spills;
  total.from_set_sorts += part.from_set_sorts;
  total.path_input_reads += part.path_input_reads;
  total.rank_input_reads += part.rank_input_reads;
  total.path_input_refills += part.path_input_refills;
  total.rank_input_refills += part.rank_input_refills;
  total.direct_input_reads += part.direct_input_reads;
  total.prefetch_submitted += part.prefetch_submitted;
  total.prefetch_completed += part.prefetch_completed;
  total.prefetch_consumed += part.prefetch_consumed;
  total.prefetch_ready_hits += part.prefetch_ready_hits;
  total.prefetch_waits += part.prefetch_waits;
  total.prefetch_synchronous_blocks += part.prefetch_synchronous_blocks;
  total.prefetch_cancelled += part.prefetch_cancelled;
  total.prefetch_errors += part.prefetch_errors;
  total.prefetch_physical_bytes += part.prefetch_physical_bytes;
  total.prefetch_decoded_bytes += part.prefetch_decoded_bytes;
  total.prefetch_wait_nanoseconds += part.prefetch_wait_nanoseconds;
}

template<class Getter>
size_type
pruneConcurrentPeak(const std::vector<PathGraphMergeStats>& parts,
  size_type workers, Getter get)
{
  std::vector<size_type> peaks;
  peaks.reserve(parts.size());
  for(const PathGraphMergeStats& part : parts) { peaks.push_back(get(part)); }
  std::sort(peaks.begin(), peaks.end(), std::greater<size_type>());
  size_type result = 0;
  for(size_type i = 0; i < std::min(workers, static_cast<size_type>(peaks.size())); i++)
  {
    if(peaks[i] > std::numeric_limits<size_type>::max() - result)
    {
      return std::numeric_limits<size_type>::max();
    }
    result += peaks[i];
  }
  return result;
}

struct PrunePartition
{
  PathNode::rank_type lower_rank, upper_rank;
  std::unique_ptr<PathGraph> graph;
  PathGraphMergeStats stats;
  PriorityNode first_output, last_output;
  bool has_output, first_range_open;
  // Per-shard record offsets of lower_rank and upper_rank.
  std::vector<size_type> start_offsets, stop_offsets;

  PrunePartition() : lower_rank(0), upper_rank(0), graph(), stats(),
    first_output(), last_output(), has_output(false), first_range_open(false) { }
};

void
buildPrunePartition(PrunePartition& result, const PathGraph& source,
  const LCP& lcp, size_type size_limit, size_type group_buffer_bytes,
  size_type input_cache_bytes, size_type input_pairs, size_type output_pairs,
  const PruneOutputLayout& layout, std::atomic<size_type>& output_bytes,
  std::atomic<bool>& cancelled, PathGraphInputCache* shared_input = nullptr)
{
  PathGraphBuilder builder(layout.logical.size(), source.k(), source.step(),
    size_limit, group_buffer_bytes, output_pairs, &(result.stats));
  builder.graph.logical_file_ids = layout.logical;
  builder.graph.physical_shard_ids = layout.physical;
  prunePathRange<false>(source, lcp, result.lower_rank, result.upper_rank,
    group_buffer_bytes, input_cache_bytes, input_pairs, layout, builder,
    &(result.stats), &output_bytes, size_limit, &cancelled,
    &(result.first_output), &(result.last_output), &(result.has_output), nullptr,
    &(result.first_range_open),
    (result.start_offsets.empty() ? nullptr : result.start_offsets.data()),
    (result.stop_offsets.empty() ? nullptr : result.stop_offsets.data()),
    shared_input);
  builder.close();

  // Retain only the lightweight graph metadata between tasks. In particular,
  // do not retain every completed builder's output-cache vector capacity; that
  // would turn a worker bound into a partition-count bound.
  result.graph.reset(new PathGraph(0, source.k(), source.step()));
  result.graph->swap(builder.graph);
}

/*
  Locate every partition bound in every shard once, before the workers start.
  Each worker's merger would otherwise binary-search every shard for its own
  two bounds, and on framed shards every probe that lands in a new block
  decodes the whole block: at 131 partitions over 33 chr18 shards that setup
  was 95% of the prune. Bounds ascend, so each shard is walked once from the
  previous bound. Shards are independent, so the pass runs on every thread the
  caller allows (it ran on the prune's seven workers and took 4.7 to 5.9 s per
  chr18 step). Partition is any type with lower_rank, upper_rank, start_offsets
  and stop_offsets; the prune and the final merge share this pass.

  The walk only moves forward, so each reader's single decoded block is
  decoded once. A galloping search overshoots and then binary-searches back
  across block boundaries, decoding the same 16 MiB blocks again on each
  crossing, which matters once bounds are a fraction of a block apart (the
  final merge's 123 to 441 key ranges over 51 chr18 shards). Here the stride
  doubles up to a cap of about a sixteenth of the expected gap between
  bounds (at least 4,096 records), and the bound is then found by a forward
  scan of at most one stride from the last record below it.
*/
template<class Partition>
void
locatePartitionBounds(const PathGraph& source, std::vector<Partition>& partitions,
  size_type pair_buffer_bytes, size_type threads)
{
  const size_type files = source.files(), parts = partitions.size();
  if(parts == 0) { return; }
  for(Partition& part : partitions)
  {
    part.start_offsets.assign(files, 0); part.stop_offsets.assign(files, 0);
  }
  forEachShard(source, pair_buffer_bytes, threads,
    [&](size_type file, PathGraphInputCache& input)
    {
      PriorityNode probe; probe.file = file;
      const size_type count = source.path_counts[file];
      auto below = [&](size_type offset, PathNode::rank_type key) -> bool
      {
        input.read(file, offset, probe.node, probe.label);
        return (probe.firstLabel(0) < key);
      };
      const size_type max_step = std::max(static_cast<size_type>(4096), count / ((parts + 1) * 16));
      size_type low = 0; // Every record before low is below the current key.
      for(size_type i = 0; i <= parts; i++)
      {
        const PathNode::rank_type key = (i < parts ?
          partitions[i].lower_rank : partitions[parts - 1].upper_rank);
        if(key != 0)
        {
          for(size_type step = 1; low < count; step = std::min(2 * step, max_step))
          {
            const size_type probe = std::min(low + step, count) - 1;
            if(!below(probe, key))
            {
              while(low < probe && below(low, key)) { low++; }
              break;
            }
            low = probe + 1;
          }
        }
        if(i < parts) { partitions[i].start_offsets[file] = low; }
        if(i > 0) { partitions[i - 1].stop_offsets[file] = low; }
      }
    });
}

// Declines a requested parallel prune or merge. The serial route that follows
// resets the stats, so the reason travels back to the caller to record.
bool
declineParallel(const char** reason_out, const char* reason)
{
  if(reason_out != nullptr) { *reason_out = reason; }
  return false;
}

bool
buildPrunedGraphParallel(PathGraph& source, const LCP& lcp,
  size_type size_limit, size_type group_buffer_bytes,
  PathGraphMergeStats* stats, size_type max_open_files,
  size_type input_cache_bytes, size_type requested_workers,
  const PruneOutputLayout& layout, const char** fallback_reason,
  size_type concurrent_open_files)
{
  std::vector<std::pair<PathNode::rank_type, PathNode::rank_type>> ranges;
  if(requested_workers < 2 || source.size() == 0 || layout.logical.empty()) { return false; }
  const size_type split_depth = pruneSplitDepth();
  if(!pruneRootPartitions(lcp, ranges, split_depth))
  {
    return declineParallel(fallback_reason, "no zero-LCP root partitions");
  }

  const size_type minimum_group = std::max(
    2 * static_cast<size_type>(sizeof(PriorityNode)),
    static_cast<size_type>(sizeof(PathNode)) +
      (PathLabel::LABEL_LENGTH + 1) * sizeof(PathNode::rank_type));
  size_type workers = std::min(requested_workers,
    static_cast<size_type>(std::max(1, omp_get_max_threads())));
  workers = std::min(workers, static_cast<size_type>(ranges.size()));
  workers = std::min(workers, group_buffer_bytes / minimum_group);

  bool all_framed = (source.files() > 0);
  for(size_type file = 0; file < source.files(); file++)
  {
    const bool path_framed = CompressedBlockReader::isFramed(source.path_names[file]);
    const bool rank_framed = CompressedBlockReader::isFramed(source.rank_names[file]);
    if(path_framed != rank_framed)
    {
      throw std::runtime_error(
        "PathGraph::prune(): path/rank shard storage formats differ");
    }
    all_framed = all_framed && path_framed;
  }

  // A worker merge alternates among every physical shard. Parallelism is only
  // useful if every worker retains the entire decoder/window set. Framed
  // readers use transient descriptors; raw/mixed readers also need all pairs
  // inside the global descriptor allowance.
  const size_type pair_bytes = pathGraphFramedPairBytes(source);
  const size_type raw_pair_budget = 8 * static_cast<size_type>(4 * KILOBYTE);
  const size_type resident_pair = std::max(pair_bytes,
    (all_framed ? static_cast<size_type>(0) : raw_pair_budget));
  if(source.files() > 0 &&
     resident_pair > std::numeric_limits<size_type>::max() / source.files())
  {
    return declineParallel(fallback_reason, "input window size overflows");
  }
  const size_type resident_per_worker = resident_pair * source.files();
  const size_type cache_budget = (input_cache_bytes == 0 ?
    group_buffer_bytes : input_cache_bytes);
  if(resident_per_worker > 0)
  {
    workers = std::min(workers, cache_budget / resident_per_worker);
  }
  // Raw readers hold their descriptors for the whole pass, so the workers'
  // descriptors are real and come from the separate concurrent budget. Framed
  // readers open theirs transiently and stay under the merge ceiling, which
  // also sizes each framed worker's input share below.
  const size_type admission_files = (all_framed ? max_open_files :
    std::max(max_open_files, concurrent_open_files));
  workers = std::min(workers, admission_files / 6);
  if(!all_framed)
  {
    if(source.files() > (std::numeric_limits<size_type>::max() - 4) / 2)
    {
      return declineParallel(fallback_reason, "too many input shards for the open-file limit");
    }
    workers = std::min(workers,
      admission_files / (2 + 2 * source.files() + 2));
  }
  if(workers < 2)
  {
    return declineParallel(fallback_reason, "memory, cache or open-file limits admit fewer than two workers");
  }

  // The exact output-pair allocation depends on the per-worker FD share. If a
  // rounded share cannot retain the full raw/mixed input set, reduce
  // concurrency rather than admit a known cache-thrashing execution.
  size_type worker_open_files = 0, output_pairs = 0, input_pairs = 0;
  while(workers >= 2)
  {
    worker_open_files = admission_files / workers;
    output_pairs = pathMergeOutputPairs(worker_open_files, layout.logical.size());
    input_pairs = (all_framed ? pathMergeInputPairs(worker_open_files,
      layout.logical.size()) : source.files());
    size_type required_fds = 2 + 2 * output_pairs +
      2 * (all_framed ? static_cast<size_type>(1) : input_pairs);
    if(required_fds <= worker_open_files) { break; }
    workers--;
  }
  if(workers < 2)
  {
    return declineParallel(fallback_reason, "open-file share cannot hold each worker's shards");
  }

  const size_type worker_group = group_buffer_bytes / workers;
  const size_type worker_cache = (input_cache_bytes == 0 ? worker_group :
    input_cache_bytes / workers);
  if(resident_per_worker > worker_cache)
  {
    return declineParallel(fallback_reason, "input cache cannot hold each worker's windows");
  }

  // The shard passes run before any worker and hold one input pair (two
  // descriptors) per thread, so they may use more threads than the workers.
  const size_type pass_pair = (pair_bytes > 0 ? pair_bytes : raw_pair_budget);
  const size_type pass_threads = std::max(workers, std::min(admission_files / 2,
    shardPassThreads(source, pair_bytes, cache_budget)));
  validatePruneShardTails(source, lcp, pass_pair, pass_threads);

  std::vector<PrunePartition> partitions(ranges.size());
  for(size_type i = 0; i < ranges.size(); i++)
  {
    partitions[i].lower_rank = ranges[i].first;
    partitions[i].upper_rank = ranges[i].second;
  }
  const auto locate_start = std::chrono::steady_clock::now();
  locatePartitionBounds(source, partitions, pass_pair, pass_threads);
  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    std::cerr << "PathGraph::prune(): located " << partitions.size() << " partition bounds in "
              << source.files() << " shard(s) in "
              << std::chrono::duration<double>(std::chrono::steady_clock::now() - locate_start).count()
              << " seconds" << std::endl;
  }

  std::atomic<size_type> output_bytes(0), next_partition(0);
  std::atomic<bool> cancelled(false);
  std::vector<std::exception_ptr> errors(partitions.size());
  // Workers take runs of adjacent partitions and keep one input cache across
  // them. Adjacent partitions start in the same compressed blocks, and each
  // shard's offsets only ascend within a worker, so a block one partition
  // decoded serves the next instead of being decoded again. Four runs per
  // worker keep the tail balanced.
  const size_type run_length = std::max(static_cast<size_type>(1),
    partitions.size() / (4 * workers));
  auto run = [&]()
  {
    std::unique_ptr<PathGraphInputCache> shared_input;
    for(size_type start = next_partition.fetch_add(run_length, std::memory_order_relaxed);
        start < partitions.size();
        start = next_partition.fetch_add(run_length, std::memory_order_relaxed))
    {
      const size_type stop = std::min(start + run_length, partitions.size());
      for(size_type i = start; i < stop; i++)
      {
        if(cancelled.load(std::memory_order_relaxed)) { return; }
        try
        {
          if(!shared_input)
          {
            shared_input.reset(new PathGraphInputCache(source, nullptr, input_pairs,
              worker_cache, 0));
          }
          buildPrunePartition(partitions[i], source, lcp, size_limit,
            worker_group, worker_cache, input_pairs, output_pairs, layout,
            output_bytes, cancelled, shared_input.get());
        }
        catch(...)
        {
          errors[i] = std::current_exception();
          cancelled.store(true, std::memory_order_relaxed);
        }
      }
    }
    if(shared_input) { shared_input->close(); }
  };

  std::vector<std::thread> threads;
  threads.reserve(workers - 1);
  std::exception_ptr startup_error;
  try
  {
    for(size_type i = 1; i < workers; i++) { threads.emplace_back(run); }
  }
  catch(const std::system_error&)
  {
    // The shared queue lets the threads that did start, plus this thread,
    // finish safely under the smaller actual concurrency.
  }
  catch(...)
  {
    startup_error = std::current_exception();
    cancelled.store(true, std::memory_order_relaxed);
  }
  const size_type actual_workers = threads.size() + 1;
  if(!startup_error) { run(); }
  for(std::thread& thread : threads) { thread.join(); }
  if(startup_error) { std::rethrow_exception(startup_error); }

  std::exception_ptr first_error;
  bool size_limit_exceeded = false;
  for(const std::exception_ptr& error : errors)
  {
    if(!error) { continue; }
    try { std::rethrow_exception(error); }
    catch(const PruneCancelled&) { }
    catch(const PruneSizeLimit&) { size_limit_exceeded = true; }
    catch(...) { if(!first_error) { first_error = error; } }
  }
  if(first_error) { std::rethrow_exception(first_error); }
  if(size_limit_exceeded)
  {
    std::cerr << "PathGraphBuilder::write(): Size limit exceeded, construction aborted" << std::endl;
    std::exit(EXIT_SIZE_LIMIT_EXCEEDED);
  }

  // Validate the actual emitted neighbors. A malformed or unusual label range
  // that spans a planned root split cannot be published as independent shards;
  // discard the attempt and let the exact serial route decide the groups.
  if(split_depth > 1)
  {
    // Every comparison inside a partition shares at least split_depth leading
    // characters, so only a partition whose first range merged all of it can
    // reach a split. The serial merger carries that range across only if
    // (a) the LCP across the right split exceeds the LCP across the left one,
    // which is the range's left_lcp; and (b) every partition below that right
    // LCP is also one merged range from the same start node and logical file,
    // because extendRange() stops at the first group that fails either test.
    // Such a span is one group in the serial pass, so prune it again as one
    // partition (stitch it) and repeat until no merged range crosses a split.
    // The keys below '#A', '#C', '#G' and '#T' all start at the end node's
    // last offset, so a split inside one of those subtrees is always stitched.
    auto split_lcp = [&](size_type left, size_type right) -> range_type
    {
      return lcp.max_lcp(partitions[left].last_output.node, partitions[right].first_output.node,
        partitions[left].last_output.label, partitions[right].first_output.label);
    };
    auto same_start = [&](size_type a, size_type b) -> bool
    {
      return partitions[a].first_output.node.from == partitions[b].first_output.node.from &&
        source.logicalFile(partitions[a].first_output.file) ==
        source.logicalFile(partitions[b].first_output.file);
    };
    size_type stitches = 0, stitched_partitions = 0, logged = 0;
    while(true)
    {
      // A node merged in an earlier step carries a label interval. If one
      // reaches past its partition's split, its neighbors' LCPs depend on keys
      // in later partitions: stitch through the partition holding its end.
      size_type stitch_first = 0, stitch_last = 0;
      bool found = false;
      std::vector<size_type> nonempty;
      for(size_type i = 0; !found && i < partitions.size(); i++)
      {
        const PrunePartition& part = partitions[i];
        if(!part.has_output) { continue; }
        const PathNode::rank_type last = part.last_output.node.lastLabel(0, part.last_output.label);
        if(last >= part.upper_rank)
        {
          size_type j = i + 1;
          while(j < partitions.size() && partitions[j].upper_rank <= last) { j++; }
          if(j >= partitions.size())
          {
            std::cerr << "PathGraph::prune(): warning: a label range reaches past the last partition; "
                      << "pruning this step serially" << std::endl;
            return declineParallel(fallback_reason, "a label range reaches past a prefix split");
          }
          found = true; stitch_first = i; stitch_last = j;
          if(Verbosity::level >= Verbosity::EXTENDED && logged < 32)
          {
            logged++;
            std::cerr << "PathGraph::prune(): partition " << i << " of " << partitions.size()
                      << " (keys " << part.lower_rank << " to " << part.upper_rank
                      << "): a label interval reaches key " << last << std::endl;
          }
          break;
        }
        nonempty.push_back(i);
      }

      for(size_type n = 0; !found && n + 1 < nonempty.size(); n++)
      {
        const size_type i = nonempty[n];
        if(!(partitions[i].first_range_open)) { continue; }
        const range_type left = (n == 0 ? range_type(0, 0) : split_lcp(nonempty[n - 1], i));
        const range_type right = split_lcp(i, nonempty[n + 1]);
        bool crosses = (right > left);
        size_type checked = n + 1;
        if(crosses)
        {
          for(; checked < nonempty.size(); checked++)
          {
            const size_type j = nonempty[checked];
            if(split_lcp(nonempty[checked - 1], j) < right) { break; }
            if(!(partitions[j].first_range_open) || !same_start(i, j)) { crosses = false; break; }
          }
        }
        if(Verbosity::level >= Verbosity::EXTENDED && logged < 32)
        {
          logged++;
          std::cerr << "PathGraph::prune(): merged partition " << i << " of " << partitions.size()
                    << " (keys " << partitions[i].lower_rank << " to " << partitions[i].upper_rank
                    << "): left split LCP (" << left.first << ", " << left.second
                    << "), right split LCP (" << right.first << ", " << right.second << "), "
                    << (crosses ? "continues past the split" :
                        (right <= left ? "closed: right split is not deeper" :
                         "closed: a later partition differs")) << std::endl;
        }
        if(crosses)
        {
          found = true; stitch_first = i; stitch_last = nonempty[checked - 1];
        }
      }
      if(!found) { break; }

      // Release the replaced partitions' output reservations, then prune the
      // span on this thread; the other workers have finished.
      PrunePartition stitched;
      stitched.lower_rank = partitions[stitch_first].lower_rank;
      stitched.upper_rank = partitions[stitch_last].upper_rank;
      stitched.start_offsets = partitions[stitch_first].start_offsets;
      stitched.stop_offsets = partitions[stitch_last].stop_offsets;
      for(size_type i = stitch_first; i <= stitch_last; i++)
      {
        if(!(partitions[i].graph)) { continue; }
        size_type released = partitions[i].graph->bytes();
        output_bytes.fetch_sub(std::min(released, output_bytes.load(std::memory_order_relaxed)),
          std::memory_order_relaxed);
      }
      try
      {
        buildPrunePartition(stitched, source, lcp, size_limit, worker_group, worker_cache,
          input_pairs, output_pairs, layout, output_bytes, cancelled);
      }
      catch(const PruneSizeLimit&)
      {
        std::cerr << "PathGraphBuilder::write(): Size limit exceeded, construction aborted" << std::endl;
        std::exit(EXIT_SIZE_LIMIT_EXCEEDED);
      }
      stitches++; stitched_partitions += stitch_last + 1 - stitch_first;
      if(Verbosity::level >= Verbosity::EXTENDED)
      {
        std::cerr << "PathGraph::prune(): stitched partitions " << stitch_first << " to " << stitch_last
                  << " (keys " << stitched.lower_rank << " to " << stitched.upper_rank << ")" << std::endl;
      }
      partitions[stitch_first] = std::move(stitched);
      partitions.erase(partitions.begin() + stitch_first + 1, partitions.begin() + stitch_last + 1);
    }
    if(stitches > 0)
    {
      std::cerr << "PathGraph::prune(): " << stitches << " stitch(es) re-pruned "
                << stitched_partitions << " partition(s) into wider spans; " << partitions.size()
                << " partition(s) remain" << std::endl;
    }
  }
  const PriorityNode* previous = nullptr;
  for(size_type i = 0; split_depth == 1 && i < partitions.size(); i++)
  {
    const PrunePartition& part = partitions[i];
    if(!part.has_output) { continue; }
    if(previous != nullptr &&
       lcp.max_lcp(previous->node, part.first_output.node,
         previous->label, part.first_output.label) != range_type(0, 0))
    {
      std::cerr << "PathGraph::prune(): warning: root partitions share a nonzero LCP; "
                << "pruning this step serially" << std::endl;
      return declineParallel(fallback_reason, "root partitions share a nonzero LCP");
    }
    previous = &(part.last_output);
  }

  PathGraph combined(0, source.k(), source.step());
  combined.delete_files = false; // It borrows partition files until commit.
  std::vector<std::vector<bool>> retain(partitions.size(),
    std::vector<bool>(layout.logical.size(), false));
  std::vector<bool> logical_retained(layout.logical.size(), false);
  for(size_type part = 0; part < partitions.size(); part++)
  {
    if(!(partitions[part].graph)) { continue; }
    for(size_type logical = 0; logical < layout.logical.size(); logical++)
    {
      if(partitions[part].graph->path_counts[logical] > 0)
      {
        retain[part][logical] = true; logical_retained[logical] = true;
      }
    }
  }
  // Preserve every semantic input even if pruning emitted no records for it.
  for(size_type logical = 0; logical < layout.logical.size(); logical++)
  {
    if(!logical_retained[logical]) { retain[0][logical] = true; }
  }

  size_type next_physical = 0;
  for(size_type part = 0; part < partitions.size(); part++)
  {
    PathGraph& graph = *(partitions[part].graph);
    combined.range_count += graph.range_count;
    combined.unique += graph.unique; combined.redundant += graph.redundant;
    combined.unsorted += graph.unsorted;
    combined.nondeterministic += graph.nondeterministic;
    for(size_type logical = 0; logical < layout.logical.size(); logical++)
    {
      if(!retain[part][logical]) { continue; }
      combined.path_names.push_back(graph.path_names[logical]);
      combined.rank_names.push_back(graph.rank_names[logical]);
      combined.path_counts.push_back(graph.path_counts[logical]);
      combined.rank_counts.push_back(graph.rank_counts[logical]);
      combined.path_checksums.push_back(graph.path_checksums[logical]);
      combined.rank_checksums.push_back(graph.rank_checksums[logical]);
      combined.logical_file_ids.push_back(layout.logical[logical]);
      combined.physical_shard_ids.push_back(
        physical_shard_id_t(next_physical++));
      combined.path_count += graph.path_counts[logical];
      combined.rank_count += graph.rank_counts[logical];
    }
  }

  if(stats != nullptr)
  {
    *stats = PathGraphMergeStats();
    std::vector<PathGraphMergeStats> part_stats;
    part_stats.reserve(partitions.size());
    for(const PrunePartition& part : partitions)
    {
      addPruneStats(*stats, part.stats); part_stats.push_back(part.stats);
    }
    stats->max_open_input_pairs = pruneConcurrentPeak(part_stats, actual_workers,
      [](const PathGraphMergeStats& value) { return value.max_open_input_pairs; });
    stats->max_open_output_pairs = pruneConcurrentPeak(part_stats, actual_workers,
      [](const PathGraphMergeStats& value) { return value.max_open_output_pairs; });
    stats->max_input_buffer_bytes = pruneConcurrentPeak(part_stats, actual_workers,
      [](const PathGraphMergeStats& value) { return value.max_input_buffer_bytes; });
    stats->prefetch_workers = pruneConcurrentPeak(part_stats, actual_workers,
      [](const PathGraphMergeStats& value) { return value.prefetch_workers; });
    stats->max_prefetch_bytes = pruneConcurrentPeak(part_stats, actual_workers,
      [](const PathGraphMergeStats& value) { return value.max_prefetch_bytes; });
    stats->oversized_input_pair_bytes = pruneConcurrentPeak(part_stats, actual_workers,
      [](const PathGraphMergeStats& value) { return value.oversized_input_pair_bytes; });
    stats->max_from_set_nodes = pruneConcurrentPeak(part_stats, actual_workers,
      [](const PathGraphMergeStats& value) { return value.max_from_set_nodes; });
    stats->prune_requested_workers = requested_workers;
    stats->prune_workers = actual_workers;
    stats->prune_partitions = partitions.size();
  }

  // Remove unretained empty pairs before transferring ownership. No payload is
  // copied or linked, so the writer's closed-file checksum identity remains
  // valid for every retained shard. All potentially allocating work above is
  // complete while the partition graphs still own their files.
  for(size_type part = 0; part < partitions.size(); part++)
  {
    PathGraph& graph = *(partitions[part].graph);
    for(size_type logical = 0; logical < layout.logical.size(); logical++)
    {
      if(!retain[part][logical])
      {
        TempFile::remove(graph.path_names[logical]);
        TempFile::remove(graph.rank_names[logical]);
      }
    }
    graph.delete_files = false;
  }

  combined.delete_files = true;
  source.clear(); source.swap(combined);
  return true;
}

} // anonymous namespace

void
PathGraph::prune(const LCP& lcp, size_type size_limit,
  size_type group_buffer_bytes, PathGraphMergeStats* stats,
  size_type max_open_files, size_type input_cache_bytes,
  size_type prune_workers, size_type concurrent_open_files)
{
  size_type old_path_count = this->size();

  // Each live input or output pair consumes two descriptors. Neither side may
  // scale with shard count, but the two sides do not need equal shares.
  if(max_open_files < 6)
  {
    throw std::runtime_error("PathGraph::prune(): max-open-files must be at least 6");
  }

  const size_type requested_workers = std::max(static_cast<size_type>(1),
    prune_workers);
  PruneOutputLayout layout = pruneOutputLayout(*this);
  const char* fallback_reason = nullptr;
  if(buildPrunedGraphParallel(*this, lcp, size_limit, group_buffer_bytes,
       stats, max_open_files, input_cache_bytes, requested_workers, layout,
       &fallback_reason, concurrent_open_files))
  {
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
    return;
  }

  // The priority group and range deque may each own one spill descriptor; the
  // rest is shared by two caches with two descriptors per entry. An even split
  // spends half of it on output entries nothing will occupy, because there is
  // one output per logical input, and starves the side that does scale with
  // shard count. That is not a fairness question: the merger visits the input
  // shards in round-robin label order, so an input cache one entry short of
  // the shard count misses on every record, and a framed miss decodes whole
  // blocks to deliver one 24-byte PathNode.
  size_type output_pairs = pathMergeOutputPairs(max_open_files,
    layout.logical.size());
  size_type input_pairs = pathMergeInputPairs(max_open_files,
    layout.logical.size());

  PathGraphBuilder builder(layout.logical.size(), this->k(), this->step(), size_limit,
    group_buffer_bytes, output_pairs, stats);
  builder.graph.logical_file_ids = layout.logical;
  builder.graph.physical_shard_ids = layout.physical;
  ProgressReporter progress("PathGraph::prune()", this->size(), "paths");
  prunePathRange<true>(*this, lcp, 0, PathLabel::NO_RANK,
    group_buffer_bytes, input_cache_bytes, input_pairs, layout, builder, stats,
    nullptr, size_limit, nullptr, nullptr, nullptr, nullptr, &progress);
  progress.finish();
  builder.close();
  if(stats != nullptr)
  {
    stats->prune_requested_workers = requested_workers;
    stats->prune_workers = 1;
    stats->prune_partitions = (this->size() == 0 ? 0 : 1);
    stats->prune_parallel_fallbacks = (requested_workers > 1 ? 1 : 0);
    stats->prune_fallback_reason = fallback_reason;
  }
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

  const bool path_framed = CompressedBlockReader::isFramed(this->path_names[file]);
  const bool rank_framed = CompressedBlockReader::isFramed(this->rank_names[file]);
  if(path_framed != rank_framed)
  {
    throw std::runtime_error("PathGraph::read(): path/rank storage formats differ");
  }
  if(path_framed)
  {
    CompressedBlockReader path_file(this->path_names[file]);
    CompressedBlockReader rank_file(this->rank_names[file]);
    size_type path_bytes = paths.size() * sizeof(PathNode);
    size_type rank_bytes = labels.size() * sizeof(PathNode::rank_type);
    if(path_file.logicalSize() != path_bytes || rank_file.logicalSize() != rank_bytes ||
       path_file.read(paths.data(), path_bytes) != path_bytes ||
       rank_file.read(labels.data(), rank_bytes) != rank_bytes)
    {
      throw std::runtime_error(
        "PathGraph::read(): compressed shard length does not match metadata");
    }
    return;
  }

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
  std::unique_ptr<OutputStream> direct; // Set when temp writes bypass the page cache.
  std::vector<Element> buffer;
  off_t written, released;

  constexpr static off_t CACHE_TAIL = 8 * MEGABYTE;
  constexpr static off_t CACHE_FLUSH = 64 * MEGABYTE;

  SequentialRecordWriter(const std::string& filename, size_type buffer_bytes) :
    name(filename), descriptor(-1), direct(), buffer(), written(0), released(0)
  {
    size_type records = std::max(static_cast<size_type>(1), buffer_bytes / sizeof(Element));
    this->buffer.reserve(records);
    if(OutputStream::directWrites())
    {
      // The writer pool keeps blocks in flight while the merge continues,
      // instead of this thread waiting in the dirty-page throttle and in an
      // fdatasync every CACHE_FLUSH bytes.
      this->direct.reset(new OutputStream(this->name,
        std::ios_base::binary | std::ios_base::out | std::ios_base::trunc));
      if(!this->direct->is_open()) { throw std::runtime_error("SameFromSet: cannot create raw set file"); }
      return;
    }
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
    if(this->direct)
    {
      this->flush();
      this->direct->close();
      const bool failed = this->direct->fail();
      this->direct.reset();
      if(failed) { throw std::runtime_error("SameFromSet: cannot close raw set file"); }
      // The same durability at close as the buffered route.
      int synced = ::open(this->name.c_str(), O_RDONLY);
      const bool sync_failed = (synced < 0 || ::fdatasync(synced) != 0);
      if(synced >= 0) { ::close(synced); }
      if(sync_failed) { throw std::runtime_error("SameFromSet: raw set sync failed"); }
      return;
    }
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
    if(this->direct)
    {
      this->direct->write(data, static_cast<std::streamsize>(bytes));
      if(!*this->direct) { throw std::runtime_error("SameFromSet: raw set write failed"); }
      DiskIO::write_volume += bytes;
      this->written += static_cast<off_t>(bytes); this->buffer.clear();
      return;
    }
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
  struct Result
  {
    std::string name;
    size_type nodes;
    bool in_memory;

    Result() : name(), nodes(0), in_memory(false) { }
    Result(const std::string& file_name, size_type count) :
      name(file_name), nodes(count), in_memory(false) { }
    explicit Result(size_type count) : name(), nodes(count), in_memory(true) { }
  };

  const PathGraphMerger& merger;
  std::string            selected;
  size_type              selected_nodes, budget, stream_buffer;
  PathGraphMergeStats*   stats;
  // Keep the selected set and the next candidate within one byte budget.
  // The vectors retain their allocation across groups, avoiding per-group I/O
  // and allocation for the common small-set case.
  std::vector<node_type> selected_memory, scratch;
  size_type              memory_records;
  bool                   selected_in_memory;

  SameFromSet(const PathGraphMerger& source, size_type group_buffer_bytes,
    PathGraphMergeStats* merge_stats) :
    merger(source), selected(), selected_nodes(0),
    budget(std::max(ExternalFixedRecordSorter::minimumBudget(sizeof(node_type)),
      group_buffer_bytes)),
    stream_buffer(std::max(static_cast<size_type>(sizeof(node_type)),
      std::min(static_cast<size_type>(64 * KILOBYTE), this->budget / 4))),
    stats(merge_stats), selected_memory(), scratch(),
    memory_records(group_buffer_bytes / (2 * sizeof(node_type))),
    selected_in_memory(false)
  {
    this->selected_memory.reserve(this->memory_records);
    this->scratch.reserve(this->memory_records);
  }

  ~SameFromSet() { if(!this->selected.empty()) { TempFile::remove(this->selected); } }

  Result fromNodes(range_type range)
  {
    size_type records = Range::length(range);
    if(records <= this->memory_records)
    {
      // Release (origin/master:src/path_graph.cpp:1100-1109) collapsed adjacent
      // runs while gathering and sorted only if more than one distinct value
      // survived. The fork dropped both. With 266,602,272 additional start
      // nodes over 399,778,113 final paths, a large share of ranges hold a
      // single distinct start node and take the zero-sort path. The result is
      // the same sorted unique set either way: adjacent collapsing cannot
      // remove a distinct value, and the sort/unique below still runs whenever
      // more than one survives.
      this->scratch.clear();
      node_type prev = ~static_cast<node_type>(0);
      for(size_type i = range.first; i <= range.second; i++)
      {
        node_type curr = this->merger.buffer.get(i).node.from;
        if(curr != prev) { this->scratch.push_back(curr); prev = curr; }
      }
      if(this->scratch.size() > 1)
      {
        std::sort(this->scratch.begin(), this->scratch.end());
        this->scratch.erase(std::unique(this->scratch.begin(), this->scratch.end()),
          this->scratch.end());
      }
      if(this->stats != nullptr)
      {
        this->stats->max_from_set_nodes = std::max(this->stats->max_from_set_nodes,
          this->scratch.size());
      }
      return Result(this->scratch.size());
    }

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
        }, nullptr, true,
        ExternalFixedRecordSorter::RecordOrder::ASCENDING_U64);
    }
    catch(...)
    {
      TempFile::remove(raw); TempFile::remove(reduced); throw;
    }
    TempFile::remove(raw); return Result(reduced, count);
  }

  bool operator() (range_type range)
  {
    Result next_set = this->fromNodes(range);
    bool equal = (next_set.nodes == this->selected_nodes);
    try
    {
      if(equal)
      {
        if(this->selected_in_memory && next_set.in_memory)
        {
          equal = (this->selected_memory == this->scratch);
        }
        else if(this->selected_in_memory)
        {
          SequentialRecordReader<node_type> right(next_set.name, this->stream_buffer);
          node_type value;
          for(size_type i = 0; equal && i < this->selected_nodes; i++)
          {
            equal = (right.next(value) && value == this->selected_memory[i]);
          }
          if(equal) { equal = !right.next(value); }
        }
        else if(next_set.in_memory)
        {
          SequentialRecordReader<node_type> left(this->selected, this->stream_buffer);
          node_type value;
          for(size_type i = 0; equal && i < this->selected_nodes; i++)
          {
            equal = (left.next(value) && value == this->scratch[i]);
          }
          if(equal) { equal = !left.next(value); }
        }
        else
        {
          SequentialRecordReader<node_type> left(this->selected, this->stream_buffer);
          SequentialRecordReader<node_type> right(next_set.name, this->stream_buffer);
          node_type a, b;
          for(size_type i = 0; equal && i < this->selected_nodes; i++)
          {
            equal = (left.next(a) && right.next(b) && a == b);
          }
          if(equal) { equal = (!left.next(a) && !right.next(b)); }
        }
      }
    }
    catch(...)
    {
      if(!next_set.in_memory) { TempFile::remove(next_set.name); }
      throw;
    }
    if(!next_set.in_memory) { TempFile::remove(next_set.name); }
    return equal;
  }

  void select(range_type range)
  {
    Result result = this->fromNodes(range);
    if(result.nodes == 0)
    {
      if(!result.in_memory) { TempFile::remove(result.name); }
      throw std::runtime_error("SameFromSet: empty selected set");
    }
    if(!this->selected.empty()) { TempFile::remove(this->selected); }
    this->selected.clear(); this->selected_nodes = result.nodes;
    this->selected_in_memory = result.in_memory;
    if(result.in_memory)
    {
      this->selected_memory.swap(this->scratch); this->scratch.clear();
    }
    else
    {
      this->selected_memory.clear(); this->selected = result.name;
    }
  }

  template<class Callback>
  node_type streamAfterFirst(Callback callback) const
  {
    if(this->selected_in_memory)
    {
      if(this->selected_memory.empty()) { throw std::runtime_error("SameFromSet: empty selected set"); }
      for(size_type i = 1; i < this->selected_memory.size(); i++) { callback(this->selected_memory[i]); }
      if(this->selected_memory.size() != this->selected_nodes)
      {
        throw std::runtime_error("SameFromSet: selected set size changed");
      }
      return this->selected_memory.front();
    }
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

struct MergedGraphSizeLimit { };

struct MergedGraphPartition
{
  PathNode::rank_type lower_rank, upper_rank;
  std::string path_name, rank_name, from_name, lcp_name;
  size_type path_count, rank_count, from_count;
  PriorityNode first_output, last_output;
  bool has_output;
  PathGraphMergeStats stats;

  MergedGraphPartition() :
    lower_rank(0), upper_rank(0), path_name(), rank_name(), from_name(), lcp_name(),
    path_count(0), rank_count(0), from_count(0), first_output(), last_output(),
    has_output(false), stats() { }

  ~MergedGraphPartition() { this->clear(); }

  MergedGraphPartition(const MergedGraphPartition&) = delete;
  MergedGraphPartition& operator=(const MergedGraphPartition&) = delete;

  void clear()
  {
    TempFile::remove(this->path_name); TempFile::remove(this->rank_name);
    TempFile::remove(this->from_name); TempFile::remove(this->lcp_name);
    this->path_count = 0; this->rank_count = 0; this->from_count = 0;
    this->has_output = false;
  }
};

static void
initializeMergedGraphNext(MergedGraph& graph, const DeBruijnGraph& mapper)
{
  for(size_type comp = 0; comp < mapper.alpha.sigma; comp++)
  {
    graph.next[comp] = mapper.charRange(comp).first;
    graph.next_from[comp] = 0;
  }
  graph.next[mapper.alpha.sigma] = ~(size_type)0;
  graph.next_from[mapper.alpha.sigma] = ~(size_type)0;
}

size_type
mergedGraphOutputBufferBytes(size_type group_buffer_bytes)
{
  return std::max(static_cast<size_type>(1), group_buffer_bytes / 8);
}

/*
  Writes records to a framed (compressed) file on its own thread. The serial
  path merge only copies records into a batch and hands full batches over a
  short queue; compression and the write happen on the stream's thread, so the
  merge thread does not wait on either unless the queue is full.
*/
template<class Element>
struct CompressedRecordStream
{
  constexpr static size_type QUEUE_DEPTH = 4;

  std::unique_ptr<CompressedBlockWriter> writer;
  std::vector<Element> batch;
  size_type batch_records;
  std::deque<std::vector<Element>> queue;
  std::mutex mtx;
  std::condition_variable ready, room;
  bool closing;
  std::exception_ptr error;
  std::thread worker;

  CompressedRecordStream(const std::string& filename, size_type buffer_bytes,
    const TempFileCodecParameters& codec) :
    writer(new CompressedBlockWriter(filename, codec.block_size,
      CompressedBlockWriter::ZSTD, codec.level, 1)),
    batch(), batch_records(std::max(static_cast<size_type>(1), buffer_bytes / sizeof(Element))),
    queue(), closing(false), error()
  {
    this->batch.reserve(this->batch_records);
    this->worker = std::thread([this]() { this->run(); });
  }

  ~CompressedRecordStream() { this->stop(); }

  void pushBack(const Element& value)
  {
    this->batch.push_back(value);
    if(this->batch.size() >= this->batch_records) { this->handOff(); }
  }

  void close()
  {
    if(!(this->batch.empty())) { this->handOff(); }
    this->stop();
    if(this->error) { std::rethrow_exception(this->error); }
    this->writer->finish();
  }

private:
  void handOff()
  {
    std::unique_lock<std::mutex> lock(this->mtx);
    this->room.wait(lock, [this]() { return this->queue.size() < QUEUE_DEPTH || this->error; });
    if(this->error) { std::rethrow_exception(this->error); }
    this->queue.push_back(std::move(this->batch));
    this->batch = std::vector<Element>();
    this->batch.reserve(this->batch_records);
    this->ready.notify_one();
  }

  void stop()
  {
    {
      std::lock_guard<std::mutex> lock(this->mtx);
      this->closing = true;
    }
    this->ready.notify_one();
    if(this->worker.joinable()) { this->worker.join(); }
  }

  void run()
  {
    while(true)
    {
      std::vector<Element> records;
      {
        std::unique_lock<std::mutex> lock(this->mtx);
        this->ready.wait(lock, [this]() { return !(this->queue.empty()) || this->closing; });
        if(this->queue.empty()) { return; }
        records = std::move(this->queue.front()); this->queue.pop_front();
      }
      this->room.notify_one();
      try
      {
        for(const Element& record : records) { this->writer->writeRecord(&record, sizeof(Element)); }
      }
      catch(...)
      {
        std::lock_guard<std::mutex> lock(this->mtx);
        this->error = std::current_exception();
        this->room.notify_all();
        return;
      }
    }
  }
};

template<class PathWriter, class RankWriter, class FromWriter>
static void
buildMergedGraphSerial(MergedGraph& result, const PathGraph& source,
  const DeBruijnGraph& mapper, const LCP& kmer_lcp, size_type size_limit,
  size_type group_buffer_bytes, PathGraphMergeStats* stats,
  size_type max_open_files, size_type input_cache_bytes,
  PathWriter& path_file, RankWriter& rank_file, FromWriter& from_file)
{
  constexpr size_type fixed_descriptors = 14;
  size_type input_pairs = std::max(static_cast<size_type>(1),
    (max_open_files - fixed_descriptors) / 2);
  size_type output_buffer_bytes = mergedGraphOutputBufferBytes(group_buffer_bytes);
  // The LCP stream stays raw: its consumer, the LCP array builder, reads it raw.
  SequentialRecordWriter<uint8_t> lcp_file(result.lcp_name, output_buffer_bytes);

  PathGraphMerger merger(source, kmer_lcp, group_buffer_bytes, stats, input_pairs,
    input_cache_bytes);
  SameFromSet same_from_set(merger, group_buffer_bytes, stats);
  size_type curr_comp = 0;
  size_type bytes = 0;
  ProgressReporter progress("MergedGraph()", source.size(), "paths");
  for(range_type range = merger.first(); !(merger.atEnd(range)); range = merger.next())
  {
    range_type path_lcp = merger.ranges.front().left_lcp;
    same_from_set.select(range);
    range = merger.extendRange(same_from_set);
    merger.mergePathNodes();
    PriorityNode curr = merger.buffer.get(range.second);

    bytes += curr.node.bytes() + (same_from_set.selected_nodes - 1) * sizeof(range_type) + 1;
    if(bytes > size_limit)
    {
      std::cerr << "MergedGraph::MergedGraph(): Size limit exceeded, construction aborted" << std::endl;
      std::exit(EXIT_SIZE_LIMIT_EXCEEDED);
    }
    curr.node.from = same_from_set.streamAfterFirst(
      [&](node_type from) { from_file.pushBack(range_type(result.path_count, from)); });
    curr.node.setPointer(result.rank_count);
    path_file.pushBack(curr.node);
    for(size_type i = 0; i < curr.node.ranks(); i++) { rank_file.pushBack(curr.label[i]); }
    lcp_file.pushBack(path_lcp.first * mapper.order() + path_lcp.second);

    while(curr.firstLabel(0) >= result.next[curr_comp])
    {
      result.next[curr_comp] = result.path_count;
      result.next_from[curr_comp] = result.from_count;
      curr_comp++;
    }
    result.path_count++;
    result.rank_count += curr.node.ranks();
    result.from_count += same_from_set.selected_nodes - 1;
    progress.advance(range.second - range.first + 1);
  }
  progress.finish();
  merger.close();
  path_file.close(); rank_file.close(); from_file.close(); lcp_file.close();
  if(stats != nullptr)
  {
    stats->merge_workers = 1;
    stats->merge_partitions = (source.size() == 0 ? 0 : 1);
  }
}

static bool
mergedGraphRootPartitions(const DeBruijnGraph& mapper, const LCP& kmer_lcp,
  std::vector<std::pair<PathNode::rank_type, PathNode::rank_type>>& result)
{
  result.clear();
  if(kmer_lcp.total_keys != mapper.size() ||
     kmer_lcp.total_keys >= PathLabel::NO_RANK)
  {
    return false;
  }

  std::vector<size_type> boundaries;
  boundaries.reserve(mapper.alpha.sigma + 1);
  for(size_type comp = 0; comp < mapper.alpha.sigma; comp++)
  {
    size_type boundary = mapper.charRange(comp).first;
    if(boundary > kmer_lcp.total_keys) { return false; }
    boundaries.push_back(boundary);
  }
  boundaries.push_back(kmer_lcp.total_keys);
  if(boundaries.size() < 3 || boundaries.front() != 0 ||
     !std::is_sorted(boundaries.begin(), boundaries.end()))
  {
    return false;
  }

  // A zero k-mer LCP is a root edge in the suffix tree. Once extendRange()
  // reaches it, parent_lcp becomes zero and cannot remain greater than the
  // current range's nonnegative left LCP. Hence no merge group can cross it.
  for(size_type i = 1; i + 1 < boundaries.size(); i++)
  {
    size_type boundary = boundaries[i];
    if(boundary > 0 && boundary < kmer_lcp.total_keys &&
       kmer_lcp.kmer_lcp[boundary] != 0)
    {
      return false;
    }
  }
  for(size_type i = 1; i < boundaries.size(); i++)
  {
    result.push_back(std::make_pair(
      static_cast<PathNode::rank_type>(boundaries[i - 1]),
      static_cast<PathNode::rank_type>(boundaries[i])));
  }
  return (result.size() > 1);
}

static bool
reserveMergedGraphBytes(std::atomic<size_type>& used, size_type amount,
  size_type limit)
{
  size_type current = used.load(std::memory_order_relaxed);
  while(true)
  {
    if(current > limit || amount > limit - current) { return false; }
    if(used.compare_exchange_weak(current, current + amount,
      std::memory_order_relaxed, std::memory_order_relaxed))
    {
      return true;
    }
  }
}

static void
buildMergedGraphPartition(MergedGraphPartition& result, const PathGraph& source,
  const DeBruijnGraph& mapper, const LCP& kmer_lcp,
  size_type group_buffer_bytes, size_type input_cache_bytes,
  size_type input_pairs, std::atomic<size_type>& output_bytes,
  size_type size_limit, std::atomic<bool>& cancelled)
{
  result.path_name = TempFile::getName("gcsa_merged_part_path");
  result.rank_name = TempFile::getName("gcsa_merged_part_rank");
  result.from_name = TempFile::getName("gcsa_merged_part_from");
  result.lcp_name = TempFile::getName("gcsa_merged_part_lcp");
  size_type output_buffer_bytes = std::max(static_cast<size_type>(1),
    group_buffer_bytes / 8);
  SequentialRecordWriter<PathNode> path_file(result.path_name, output_buffer_bytes);
  SequentialRecordWriter<PathNode::rank_type> rank_file(result.rank_name, output_buffer_bytes);
  SequentialRecordWriter<range_type> from_file(result.from_name, output_buffer_bytes);
  SequentialRecordWriter<uint8_t> lcp_file(result.lcp_name, output_buffer_bytes);
  PathGraphMerger merger(source, kmer_lcp, group_buffer_bytes, &(result.stats),
    input_pairs, input_cache_bytes, result.lower_rank, result.upper_rank, 0);
  SameFromSet same_from_set(merger, group_buffer_bytes, &(result.stats));

  for(range_type range = merger.first(); !(merger.atEnd(range)); range = merger.next())
  {
    if(cancelled.load(std::memory_order_relaxed)) { break; }
    range_type path_lcp = merger.ranges.front().left_lcp;
    same_from_set.select(range);
    range = merger.extendRange(same_from_set);
    merger.mergePathNodes();
    PriorityNode curr = merger.buffer.get(range.second);
    if(same_from_set.selected_nodes == 0)
    {
      throw std::runtime_error("MergedGraph: empty parallel from-set");
    }
    size_type extra_nodes = same_from_set.selected_nodes - 1;
    if(extra_nodes > std::numeric_limits<size_type>::max() / sizeof(range_type))
    {
      throw std::overflow_error("MergedGraph: parallel output byte count overflow");
    }
    size_type record_bytes = curr.node.bytes();
    size_type from_bytes = extra_nodes * sizeof(range_type);
    if(record_bytes > std::numeric_limits<size_type>::max() - from_bytes - 1)
    {
      throw std::overflow_error("MergedGraph: parallel output byte count overflow");
    }
    if(!reserveMergedGraphBytes(output_bytes, record_bytes + from_bytes + 1,
      size_limit))
    {
      throw MergedGraphSizeLimit();
    }

    PriorityNode boundary = curr; boundary.node.setPointer(0);
    if(!result.has_output) { result.first_output = boundary; result.has_output = true; }
    result.last_output = boundary;
    curr.node.from = same_from_set.streamAfterFirst(
      [&](node_type from) { from_file.pushBack(range_type(result.path_count, from)); });
    curr.node.setPointer(result.rank_count);
    path_file.pushBack(curr.node);
    for(size_type i = 0; i < curr.node.ranks(); i++) { rank_file.pushBack(curr.label[i]); }
    lcp_file.pushBack(path_lcp.first * mapper.order() + path_lcp.second);
    result.path_count++;
    result.rank_count += curr.node.ranks();
    result.from_count += extra_nodes;
  }
  merger.close();
  path_file.close(); rank_file.close(); from_file.close(); lcp_file.close();
}

static size_type
mergedRecordBytes(size_type records, size_type width, const char* context)
{
  if(width != 0 && records > std::numeric_limits<size_type>::max() / width)
  {
    throw std::overflow_error(std::string("MergedGraph: ") + context + " byte count overflow");
  }
  size_type bytes = records * width;
  if(bytes > static_cast<size_type>(std::numeric_limits<off_t>::max()))
  {
    throw std::overflow_error(std::string("MergedGraph: ") + context + " exceeds file offsets");
  }
  return bytes;
}

static void
checkMergedPartSize(const std::string& name, size_type bytes)
{
  struct stat info;
  if(::stat(name.c_str(), &info) != 0 || info.st_size < 0 ||
     static_cast<size_type>(info.st_size) != bytes)
  {
    throw std::runtime_error("MergedGraph: partition stream length mismatch");
  }
}

static void
readMergedPart(int descriptor, void* target, size_type bytes, size_type offset)
{
  char* data = reinterpret_cast<char*>(target); size_type done = 0;
  while(done < bytes)
  {
    ssize_t got = ::pread(descriptor, data + done, bytes - done,
      static_cast<off_t>(offset + done));
    if(got < 0 && errno == EINTR) { continue; }
    if(got <= 0) { throw std::runtime_error("MergedGraph: truncated partition stream"); }
    DiskIO::read_volume += static_cast<size_type>(got);
    done += static_cast<size_type>(got);
  }
}

static void
writeMergedPart(int descriptor, const void* source, size_type bytes,
  size_type offset)
{
  const char* data = reinterpret_cast<const char*>(source); size_type done = 0;
  while(done < bytes)
  {
    ssize_t put = ::pwrite(descriptor, data + done, bytes - done,
      static_cast<off_t>(offset + done));
    if(put < 0 && errno == EINTR) { continue; }
    if(put <= 0) { throw std::runtime_error("MergedGraph: partition stream write failed"); }
    DiskIO::write_volume += static_cast<size_type>(put);
    done += static_cast<size_type>(put);
  }
}

static int
openMergedPart(const std::string& name, size_type expected_bytes)
{
  checkMergedPartSize(name, expected_bytes);
  int descriptor = ::open(name.c_str(), O_RDONLY);
  if(descriptor < 0) { throw std::runtime_error("MergedGraph: cannot open partition stream"); }
#if defined(POSIX_FADV_SEQUENTIAL)
  static_cast<void>(::posix_fadvise(descriptor, 0, 0, POSIX_FADV_SEQUENTIAL));
#endif
  return descriptor;
}

constexpr size_type MERGED_COPY_CACHE_TAIL_BYTES = 8 * MEGABYTE;
constexpr size_type MERGED_COPY_CACHE_FLUSH_BYTES = 64 * MEGABYTE;

static void
discardMergedCache(int descriptor, size_type begin, size_type end)
{
#if defined(POSIX_FADV_DONTNEED)
  if(end > begin)
  {
    static_cast<void>(::posix_fadvise(descriptor, static_cast<off_t>(begin),
      static_cast<off_t>(end - begin), POSIX_FADV_DONTNEED));
  }
#else
  static_cast<void>(descriptor); static_cast<void>(begin); static_cast<void>(end);
#endif
}

static void
trimMergedReadCache(int descriptor, size_type consumed, size_type& released,
  bool complete = false)
{
  if(consumed < released)
  {
    throw std::runtime_error("MergedGraph: non-monotone partition read");
  }
  if(!complete && consumed - released < MERGED_COPY_CACHE_FLUSH_BYTES) { return; }
  size_type discard_end = (complete ? consumed :
    std::max(released, consumed - MERGED_COPY_CACHE_TAIL_BYTES));
  discardMergedCache(descriptor, released, discard_end); released = discard_end;
}

static void
trimMergedWrittenCache(int descriptor, size_type written, size_type& released,
  bool complete = false)
{
  if(written < released)
  {
    throw std::runtime_error("MergedGraph: non-monotone partition write");
  }
  if(written == released ||
     (!complete && written - released < MERGED_COPY_CACHE_FLUSH_BYTES))
  {
    return;
  }
  // The output ranges of different partitions are disjoint. fdatasync() may
  // flush concurrent ranges on this shared descriptor as well, which is safe;
  // each worker advises away only its own completed prefix.
  if(::fdatasync(descriptor) != 0)
  {
    throw std::runtime_error("MergedGraph: cannot sync parallel output range");
  }
  size_type discard_end = (complete ? written :
    std::max(released, written - MERGED_COPY_CACHE_TAIL_BYTES));
  discardMergedCache(descriptor, released, discard_end); released = discard_end;
}

static void
copyMergedRaw(const std::string& input_name, size_type input_bytes,
  int output, size_type output_offset, size_type buffer_bytes)
{
  int input = openMergedPart(input_name, input_bytes);
  try
  {
    std::vector<char> buffer(std::max(static_cast<size_type>(1), buffer_bytes));
    size_type input_released = 0, output_released = output_offset;
    for(size_type offset = 0; offset < input_bytes; )
    {
      size_type bytes = std::min(static_cast<size_type>(buffer.size()), input_bytes - offset);
      readMergedPart(input, buffer.data(), bytes, offset);
      writeMergedPart(output, buffer.data(), bytes, output_offset + offset);
      offset += bytes;
      trimMergedReadCache(input, offset, input_released);
      trimMergedWrittenCache(output, output_offset + offset, output_released);
    }
    trimMergedReadCache(input, input_bytes, input_released, true);
    trimMergedWrittenCache(output, output_offset + input_bytes,
      output_released, true);
    ::close(input);
  }
  catch(...)
  {
    ::close(input); throw;
  }
}

static void
copyMergedPartition(const MergedGraphPartition& part,
  size_type path_base, size_type rank_base, size_type from_base,
  int path_output, int rank_output, int from_output, int lcp_output,
  size_type buffer_bytes)
{
  const size_type path_bytes = mergedRecordBytes(part.path_count, sizeof(PathNode), "path partition");
  const size_type rank_bytes = mergedRecordBytes(part.rank_count, sizeof(PathNode::rank_type), "rank partition");
  const size_type from_bytes = mergedRecordBytes(part.from_count, sizeof(range_type), "from partition");
  const size_type lcp_bytes = mergedRecordBytes(part.path_count, sizeof(uint8_t), "LCP partition");

  int input = openMergedPart(part.path_name, path_bytes);
  try
  {
    size_type records_per_block = std::max(static_cast<size_type>(1),
      buffer_bytes / sizeof(PathNode));
    std::vector<PathNode> buffer(records_per_block);
    size_type input_released = 0;
    size_type output_released = path_base * sizeof(PathNode);
    for(size_type offset = 0; offset < part.path_count; )
    {
      size_type records = std::min(records_per_block, part.path_count - offset);
      size_type bytes = records * sizeof(PathNode);
      readMergedPart(input, buffer.data(), bytes, offset * sizeof(PathNode));
      for(size_type i = 0; i < records; i++)
      {
        if(buffer[i].pointer() > std::numeric_limits<size_type>::max() - rank_base ||
           buffer[i].pointer() + rank_base > ((static_cast<size_type>(1) << 40) - 1))
        {
          throw std::overflow_error("MergedGraph: rebased rank pointer overflow");
        }
        buffer[i].setPointer(buffer[i].pointer() + rank_base);
      }
      writeMergedPart(path_output, buffer.data(), bytes,
        (path_base + offset) * sizeof(PathNode));
      offset += records;
      trimMergedReadCache(input, offset * sizeof(PathNode), input_released);
      trimMergedWrittenCache(path_output,
        (path_base + offset) * sizeof(PathNode), output_released);
    }
    trimMergedReadCache(input, path_bytes, input_released, true);
    trimMergedWrittenCache(path_output,
      (path_base + part.path_count) * sizeof(PathNode), output_released, true);
    ::close(input);
  }
  catch(...)
  {
    ::close(input); throw;
  }

  copyMergedRaw(part.rank_name, rank_bytes, rank_output,
    rank_base * sizeof(PathNode::rank_type), buffer_bytes);

  input = openMergedPart(part.from_name, from_bytes);
  try
  {
    size_type records_per_block = std::max(static_cast<size_type>(1),
      buffer_bytes / sizeof(range_type));
    std::vector<range_type> buffer(records_per_block);
    size_type input_released = 0;
    size_type output_released = from_base * sizeof(range_type);
    for(size_type offset = 0; offset < part.from_count; )
    {
      size_type records = std::min(records_per_block, part.from_count - offset);
      size_type bytes = records * sizeof(range_type);
      readMergedPart(input, buffer.data(), bytes, offset * sizeof(range_type));
      for(size_type i = 0; i < records; i++)
      {
        if(buffer[i].first > std::numeric_limits<size_type>::max() - path_base)
        {
          throw std::overflow_error("MergedGraph: rebased from pointer overflow");
        }
        buffer[i].first += path_base;
      }
      writeMergedPart(from_output, buffer.data(), bytes,
        (from_base + offset) * sizeof(range_type));
      offset += records;
      trimMergedReadCache(input, offset * sizeof(range_type), input_released);
      trimMergedWrittenCache(from_output,
        (from_base + offset) * sizeof(range_type), output_released);
    }
    trimMergedReadCache(input, from_bytes, input_released, true);
    trimMergedWrittenCache(from_output,
      (from_base + part.from_count) * sizeof(range_type), output_released, true);
    ::close(input);
  }
  catch(...)
  {
    ::close(input); throw;
  }

  copyMergedRaw(part.lcp_name, lcp_bytes, lcp_output,
    path_base * sizeof(uint8_t), buffer_bytes);
}

static void
addMergedGraphStats(PathGraphMergeStats& total, const PathGraphMergeStats& part)
{
  total.priority_spills += part.priority_spills;
  total.range_spills += part.range_spills;
  total.from_set_sorts += part.from_set_sorts;
  total.path_input_reads += part.path_input_reads;
  total.rank_input_reads += part.rank_input_reads;
  total.path_input_refills += part.path_input_refills;
  total.rank_input_refills += part.rank_input_refills;
  total.direct_input_reads += part.direct_input_reads;
  total.prefetch_submitted += part.prefetch_submitted;
  total.prefetch_completed += part.prefetch_completed;
  total.prefetch_consumed += part.prefetch_consumed;
  total.prefetch_ready_hits += part.prefetch_ready_hits;
  total.prefetch_waits += part.prefetch_waits;
  total.prefetch_synchronous_blocks += part.prefetch_synchronous_blocks;
  total.prefetch_cancelled += part.prefetch_cancelled;
  total.prefetch_errors += part.prefetch_errors;
  total.prefetch_physical_bytes += part.prefetch_physical_bytes;
  total.prefetch_decoded_bytes += part.prefetch_decoded_bytes;
  total.prefetch_wait_nanoseconds += part.prefetch_wait_nanoseconds;
}

template<class Getter>
static size_type
mergedConcurrentPeak(const std::vector<MergedGraphPartition>& partitions,
  size_type workers, Getter get)
{
  std::vector<size_type> peaks;
  peaks.reserve(partitions.size());
  for(const MergedGraphPartition& part : partitions) { peaks.push_back(get(part.stats)); }
  std::sort(peaks.begin(), peaks.end(), std::greater<size_type>());
  size_type result = 0;
  for(size_type i = 0; i < std::min(workers, static_cast<size_type>(peaks.size())); i++)
  {
    if(peaks[i] > std::numeric_limits<size_type>::max() - result)
    {
      return std::numeric_limits<size_type>::max();
    }
    result += peaks[i];
  }
  return result;
}

static bool
buildMergedGraphParallel(MergedGraph& result, const PathGraph& source,
  const DeBruijnGraph& mapper, const LCP& kmer_lcp, size_type size_limit,
  size_type group_buffer_bytes, PathGraphMergeStats* stats,
  size_type max_open_files, size_type input_cache_bytes,
  size_type requested_workers, const char** fallback_reason)
{
  std::vector<std::pair<PathNode::rank_type, PathNode::rank_type>> ranges;
  if(requested_workers < 2) { return false; }
  if(!mergedGraphRootPartitions(mapper, kmer_lcp, ranges))
  {
    return declineParallel(fallback_reason, "no zero-LCP root partitions");
  }

  constexpr size_type fixed_descriptors = 14;
  constexpr size_type minimum_descriptors = fixed_descriptors + 2;
  size_type active_ranges = 0;
  for(const auto& range : ranges)
  {
    if(range.first < range.second) { active_ranges++; }
  }
  size_type workers = std::min(requested_workers, active_ranges);
  workers = std::min(workers,
    static_cast<size_type>(std::max(1, omp_get_max_threads())));
  workers = std::min(workers, max_open_files / minimum_descriptors);
  bool all_framed = (source.files() > 0);
  for(size_type file = 0; all_framed && file < source.files(); file++)
  {
    all_framed = CompressedBlockReader::isFramed(source.path_names[file]) &&
      CompressedBlockReader::isFramed(source.rank_names[file]);
  }
  if(!all_framed && source.files() > 0)
  {
    // Raw/mixed generations use persistent descriptors. Keep every shard's
    // input window resident or fall back instead of creating an eviction loop.
    if(source.files() >
       (std::numeric_limits<size_type>::max() - fixed_descriptors) / 2)
    {
      return declineParallel(fallback_reason, "too many input shards for the open-file limit");
    }
    workers = std::min(workers,
      max_open_files / (fixed_descriptors + 2 * source.files()));
  }
  const size_type minimum_group = ExternalFixedRecordSorter::minimumBudget(sizeof(node_type));
  if(group_buffer_bytes < minimum_group)
  {
    return declineParallel(fallback_reason, "group buffer below the sorter minimum");
  }
  workers = std::min(workers, group_buffer_bytes / minimum_group);
  size_type pair_bytes = pathGraphFramedPairBytes(source);
  size_type cache_budget = (input_cache_bytes == 0 ? group_buffer_bytes : input_cache_bytes);
  if(pair_bytes > 0)
  {
    size_type shards = std::max(static_cast<size_type>(1), source.files());
    if(pair_bytes > std::numeric_limits<size_type>::max() / shards)
    {
      return declineParallel(fallback_reason, "input cache size overflows");
    }
    size_type resident_bytes = pair_bytes * shards;
    if(cache_budget < resident_bytes)
    {
      return declineParallel(fallback_reason, "input cache cannot hold every shard's window");
    }
    workers = std::min(workers, cache_budget / resident_bytes);
  }
  if(workers < 2)
  {
    return declineParallel(fallback_reason, "memory, cache or open-file limits admit fewer than two workers");
  }

  size_type worker_group = std::max(static_cast<size_type>(1), group_buffer_bytes / workers);
  size_type worker_cache = (input_cache_bytes == 0 ? 0 :
    std::max(static_cast<size_type>(1), input_cache_bytes / workers));
  size_type worker_open_files = max_open_files / workers;
  size_type input_pairs = std::max(static_cast<size_type>(1),
    (worker_open_files - fixed_descriptors) / 2);
  std::vector<MergedGraphPartition> partitions(ranges.size());
  for(size_type i = 0; i < ranges.size(); i++)
  {
    partitions[i].lower_rank = ranges[i].first;
    partitions[i].upper_rank = ranges[i].second;
  }

  std::atomic<size_type> output_bytes(0);
  std::atomic<bool> cancelled(false);
  std::vector<std::exception_ptr> errors(partitions.size());
  #pragma omp parallel for schedule(dynamic, 1) num_threads(static_cast<int>(workers))
  for(std::int64_t i = 0; i < static_cast<std::int64_t>(partitions.size()); i++)
  {
    if(cancelled.load(std::memory_order_relaxed)) { continue; }
    if(partitions[i].lower_rank == partitions[i].upper_rank) { continue; }
    try
    {
      buildMergedGraphPartition(partitions[i], source, mapper, kmer_lcp,
        worker_group, worker_cache, input_pairs, output_bytes, size_limit,
        cancelled);
    }
    catch(...)
    {
      errors[i] = std::current_exception();
      cancelled.store(true, std::memory_order_relaxed);
    }
  }

  std::exception_ptr first_error;
  bool size_limit_exceeded = false;
  for(const std::exception_ptr& error : errors)
  {
    if(!error) { continue; }
    try { std::rethrow_exception(error); }
    catch(const MergedGraphSizeLimit&) { size_limit_exceeded = true; }
    catch(...) { if(!first_error) { first_error = error; } }
  }
  if(first_error || size_limit_exceeded)
  {
    for(MergedGraphPartition& part : partitions) { part.clear(); }
    if(first_error) { std::rethrow_exception(first_error); }
    std::cerr << "MergedGraph::MergedGraph(): Size limit exceeded, construction aborted" << std::endl;
    std::exit(EXIT_SIZE_LIMIT_EXCEEDED);
  }

  // Validate the actual emitted neighbors as well as the planned root split.
  // If an unusual alphabet violates the zero-LCP premise, discard all local
  // streams and let the byte-for-byte serial implementation decide the groups.
  const PriorityNode* previous = nullptr;
  for(const MergedGraphPartition& part : partitions)
  {
    if(!part.has_output) { continue; }
    if(previous != nullptr &&
       kmer_lcp.max_lcp(previous->node, part.first_output.node,
         previous->label, part.first_output.label) != range_type(0, 0))
    {
      for(MergedGraphPartition& cleanup : partitions) { cleanup.clear(); }
      std::cerr << "MergedGraph: warning: root partitions share a nonzero LCP; "
                << "merging serially" << std::endl;
      return declineParallel(fallback_reason, "root partitions share a nonzero LCP");
    }
    previous = &(part.last_output);
  }

  std::vector<size_type> path_base(partitions.size() + 1, 0);
  std::vector<size_type> rank_base(partitions.size() + 1, 0);
  std::vector<size_type> from_base(partitions.size() + 1, 0);
  for(size_type i = 0; i < partitions.size(); i++)
  {
    if(partitions[i].path_count > std::numeric_limits<size_type>::max() - path_base[i] ||
       partitions[i].rank_count > std::numeric_limits<size_type>::max() - rank_base[i] ||
       partitions[i].from_count > std::numeric_limits<size_type>::max() - from_base[i])
    {
      for(MergedGraphPartition& cleanup : partitions) { cleanup.clear(); }
      throw std::overflow_error("MergedGraph: parallel prefix count overflow");
    }
    path_base[i + 1] = path_base[i] + partitions[i].path_count;
    rank_base[i + 1] = rank_base[i] + partitions[i].rank_count;
    from_base[i + 1] = from_base[i] + partitions[i].from_count;
  }
  result.path_count = path_base.back();
  result.rank_count = rank_base.back();
  result.from_count = from_base.back();
  if(result.bytes() != output_bytes.load(std::memory_order_relaxed))
  {
    for(MergedGraphPartition& cleanup : partitions) { cleanup.clear(); }
    throw std::runtime_error("MergedGraph: parallel output byte accounting mismatch");
  }

  // The partition streams remain live while the four final streams are
  // materialized. Charge both complete logical copies before ftruncate() can
  // grow the final files. If the remaining construction budget admits the
  // output but not this parallel-copy peak, remove the partitions and use the
  // one-copy serial route instead.
  const size_type merged_bytes = result.bytes();
  if(merged_bytes > size_limit - merged_bytes)
  {
    for(MergedGraphPartition& cleanup : partitions) { cleanup.clear(); }
    result.path_count = 0; result.rank_count = 0; result.from_count = 0;
    return declineParallel(fallback_reason, "disk limit cannot hold two copies of the merged graph");
  }

  bool has_output = false;
  PathNode::rank_type max_label = 0;
  for(const MergedGraphPartition& part : partitions)
  {
    if(part.has_output)
    {
      max_label = part.last_output.firstLabel(0); has_output = true;
    }
  }
  for(size_type comp = 0; comp < mapper.alpha.sigma &&
      has_output && result.next[comp] <= max_label; comp++)
  {
    PathNode::rank_type threshold = static_cast<PathNode::rank_type>(result.next[comp]);
    auto found = std::lower_bound(ranges.begin(), ranges.end(), threshold,
      [](const std::pair<PathNode::rank_type, PathNode::rank_type>& range,
         PathNode::rank_type value) { return range.first < value; });
    if(found == ranges.end() || found->first != threshold)
    {
      for(MergedGraphPartition& cleanup : partitions) { cleanup.clear(); }
      throw std::runtime_error("MergedGraph: component boundary is not a merge partition");
    }
    size_type partition = static_cast<size_type>(found - ranges.begin());
    result.next[comp] = path_base[partition];
    result.next_from[comp] = from_base[partition];
  }

  std::array<int, 4> outputs = {{ -1, -1, -1, -1 }};
  auto close_outputs = [&]()
  {
    for(int& descriptor : outputs)
    {
      if(descriptor >= 0) { ::close(descriptor); descriptor = -1; }
    }
  };
  auto remove_outputs = [&]()
  {
    TempFile::remove(result.path_name); TempFile::remove(result.rank_name);
    TempFile::remove(result.from_name); TempFile::remove(result.lcp_name);
  };
  try
  {
    const std::array<std::string*, 4> names = {{
      &(result.path_name), &(result.rank_name), &(result.from_name), &(result.lcp_name) }};
    const std::array<size_type, 4> sizes = {{
      mergedRecordBytes(result.path_count, sizeof(PathNode), "path output"),
      mergedRecordBytes(result.rank_count, sizeof(PathNode::rank_type), "rank output"),
      mergedRecordBytes(result.from_count, sizeof(range_type), "from output"),
      mergedRecordBytes(result.path_count, sizeof(uint8_t), "LCP output") }};
    for(size_type i = 0; i < outputs.size(); i++)
    {
      outputs[i] = ::open(names[i]->c_str(), O_CREAT | O_TRUNC | O_RDWR, 0600);
      if(outputs[i] < 0 || ::ftruncate(outputs[i], static_cast<off_t>(sizes[i])) != 0)
      {
        throw std::runtime_error("MergedGraph: cannot create parallel output stream");
      }
    }

    cancelled.store(false, std::memory_order_relaxed);
    std::fill(errors.begin(), errors.end(), std::exception_ptr());
    size_type copy_buffer = std::max(static_cast<size_type>(sizeof(PathNode)),
      std::min(static_cast<size_type>(64 * KILOBYTE), worker_group / 4));
    #pragma omp parallel for schedule(dynamic, 1) num_threads(static_cast<int>(workers))
    for(std::int64_t i = 0; i < static_cast<std::int64_t>(partitions.size()); i++)
    {
      if(cancelled.load(std::memory_order_relaxed)) { continue; }
      if(partitions[i].lower_rank == partitions[i].upper_rank) { continue; }
      try
      {
        copyMergedPartition(partitions[i], path_base[i], rank_base[i], from_base[i],
          outputs[0], outputs[1], outputs[2], outputs[3], copy_buffer);
      }
      catch(...)
      {
        errors[i] = std::current_exception();
        cancelled.store(true, std::memory_order_relaxed);
      }
    }
    for(const std::exception_ptr& error : errors)
    {
      if(error) { std::rethrow_exception(error); }
    }
    for(int descriptor : outputs)
    {
      if(::fdatasync(descriptor) != 0)
      {
        throw std::runtime_error("MergedGraph: cannot sync parallel output stream");
      }
    }
    for(size_type i = 0; i < outputs.size(); i++)
    {
      discardMergedCache(outputs[i], 0, sizes[i]);
    }
    close_outputs();
  }
  catch(...)
  {
    close_outputs(); remove_outputs();
    for(MergedGraphPartition& part : partitions) { part.clear(); }
    throw;
  }

  if(stats != nullptr)
  {
    *stats = PathGraphMergeStats();
    for(const MergedGraphPartition& part : partitions)
    {
      addMergedGraphStats(*stats, part.stats);
    }
    // Tasks beyond the worker count run later. Sum only the largest `workers`
    // local peaks so the diagnostics remain a conservative simultaneous bound
    // instead of growing with the number of root-prefix partitions.
    stats->max_open_input_pairs = mergedConcurrentPeak(partitions, workers,
      [](const PathGraphMergeStats& value) { return value.max_open_input_pairs; });
    stats->max_open_output_pairs = mergedConcurrentPeak(partitions, workers,
      [](const PathGraphMergeStats& value) { return value.max_open_output_pairs; });
    stats->max_input_buffer_bytes = mergedConcurrentPeak(partitions, workers,
      [](const PathGraphMergeStats& value) { return value.max_input_buffer_bytes; });
    stats->prefetch_workers = mergedConcurrentPeak(partitions, workers,
      [](const PathGraphMergeStats& value) { return value.prefetch_workers; });
    stats->max_prefetch_bytes = mergedConcurrentPeak(partitions, workers,
      [](const PathGraphMergeStats& value) { return value.max_prefetch_bytes; });
    stats->oversized_input_pair_bytes = mergedConcurrentPeak(partitions, workers,
      [](const PathGraphMergeStats& value) { return value.oversized_input_pair_bytes; });
    stats->max_from_set_nodes = mergedConcurrentPeak(partitions, workers,
      [](const PathGraphMergeStats& value) { return value.max_from_set_nodes; });
    stats->merge_workers = workers;
    stats->merge_partitions = partitions.size();
  }
  for(MergedGraphPartition& part : partitions) { part.clear(); }
  return true;
}

MergedGraph::MergedGraph(const PathGraph& source, const DeBruijnGraph& mapper,
  const LCP& kmer_lcp, size_type size_limit, size_type group_buffer_bytes,
  PathGraphMergeStats* stats, size_type max_open_files,
  size_type input_cache_bytes, size_type merge_workers,
  const TempFileCodecParameters* output_codec) :
  path_name(TempFile::getName(PREFIX)), rank_name(TempFile::getName(PREFIX)),
  from_name(TempFile::getName(PREFIX)), lcp_name(TempFile::getName(PREFIX)),
  path_count(0), rank_count(0), from_count(0),
  order(source.k()),
  next(mapper.alpha.sigma + 1, 0), next_from(mapper.alpha.sigma + 1, 0)
{
  // Four sequential outputs coexist with at most eight descriptors in the
  // two-way SameFromSet sorter (its source remains open during run merging)
  // and two optional prune spill descriptors. A
  // merger input pair needs two more. Account those fixed descriptors before
  // sizing the LRU input cache instead of silently exceeding --max-open-files.
  constexpr size_type fixed_descriptors = 14;
  if(max_open_files < fixed_descriptors + 2)
  {
    throw std::runtime_error("MergedGraph: max-open-files must be at least 16");
  }
  static_cast<void>(fixed_descriptors);
  initializeMergedGraphNext(*this, mapper);
  const char* fallback_reason = nullptr;
  if(!buildMergedGraphParallel(*this, source, mapper, kmer_lcp, size_limit,
      group_buffer_bytes, stats, max_open_files, input_cache_bytes,
      merge_workers, &fallback_reason))
  {
    this->path_count = 0; this->rank_count = 0; this->from_count = 0;
    initializeMergedGraphNext(*this, mapper);
    const size_type output_buffer_bytes = mergedGraphOutputBufferBytes(group_buffer_bytes);
    if(output_codec != nullptr && output_codec->enabled())
    {
      // Paths, ranks and start nodes are framed, each compressed on its own
      // thread; every reader goes through ReadBuffer, which reads both forms.
      CompressedRecordStream<PathNode> path_file(this->path_name, output_buffer_bytes, *output_codec);
      CompressedRecordStream<PathNode::rank_type> rank_file(this->rank_name, output_buffer_bytes, *output_codec);
      CompressedRecordStream<range_type> from_file(this->from_name, output_buffer_bytes, *output_codec);
      buildMergedGraphSerial(*this, source, mapper, kmer_lcp, size_limit,
        group_buffer_bytes, stats, max_open_files, input_cache_bytes,
        path_file, rank_file, from_file);
    }
    else
    {
      SequentialRecordWriter<PathNode> path_file(this->path_name, output_buffer_bytes);
      SequentialRecordWriter<PathNode::rank_type> rank_file(this->rank_name, output_buffer_bytes);
      SequentialRecordWriter<range_type> from_file(this->from_name, output_buffer_bytes);
      buildMergedGraphSerial(*this, source, mapper, kmer_lcp, size_limit,
        group_buffer_bytes, stats, max_open_files, input_cache_bytes,
        path_file, rank_file, from_file);
    }
    if(stats != nullptr) { stats->merge_fallback_reason = fallback_reason; }
  }

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
