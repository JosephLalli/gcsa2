#ifndef GCSA_PATH_GRAPH_H
#define GCSA_PATH_GRAPH_H

#include <gcsa/dbg.h>
#include <gcsa/files.h>
#include <gcsa/internal.h>
#include <gcsa/workspace.h>

#include <sdsl/wt_blcd.hpp>

namespace gcsa
{

/*
  path_graph.h: Internal graph manipulation methods.
*/

//------------------------------------------------------------------------------

struct PathLabel
{
  typedef std::uint32_t rank_type;

  // This should be at least 1 << ConstructionParameters::MAX_STEPS.
  constexpr static size_type LABEL_LENGTH = 16;

  // Labels starting with NO_RANK will be after real labels in lexicographic order.
  // We also use NO_RANK for padding last labels.
  constexpr static rank_type NO_RANK = ~(rank_type)0;

  rank_type label[LABEL_LENGTH];
  bool      first;

  inline bool operator< (const PathLabel& another) const
  {
    for(size_type i = 0; i < LABEL_LENGTH; i++)
    {
      if(this->label[i] != another.label[i]) { return (this->label[i] < another.label[i]); }
    }
    return (this->first && !(another.first));
  }

  inline bool operator<= (const PathLabel& another) const
  {
    for(size_type i = 0; i < LABEL_LENGTH; i++)
    {
      if(this->label[i] != another.label[i]) { return (this->label[i] < another.label[i]); }
    }
    return (this->first || !(another.first));
  }
};

//------------------------------------------------------------------------------

/*
  The node type used during doubling. As in the original GCSA, from and to are nodes
  in the original graph, denoting a path as a semiopen range [from, to). If
  to == -1, the path will not be extended, because it already has a unique label.
  rank_type is the integer type used to store ranks of the original kmers.
  During edge generation, 'to' node will be used to store indegree and the outdegree.

  The rank sequences are stored in an std::vector at position 'pointer()'. The stored
  sequence consists of the first label ('order()' ranks) followed by one rank for the
  diverging last rank of the last label. If the first and the last ranks are identical,
  the last rank is a dummy value.

  There are also alternative versions of most operations using raw pointers to the
  label array.
*/

struct PathNode
{
  typedef PathLabel::rank_type rank_type;

  constexpr static size_type LABEL_LENGTH = PathLabel::LABEL_LENGTH;

  node_type from, to;

  inline bool sorted() const { return (this->to == ~(node_type)0); }
  inline void makeSorted() { this->to = ~(node_type)0; }

//------------------------------------------------------------------------------

  /*
    From low-order to high-order bits:

    8 bits   which predecessor comp values exist
    8 bits   length of the kmer rank sequences representing the path label range
    8 bits   lcp of the above sequences
    40 bits  pointer to the label data
  */
  size_type fields;

  inline byte_type predecessors() const { return (this->fields & 0xFF); }
  inline void setPredecessors(byte_type preds)
  {
    this->fields &= ~(size_type)0xFF;
    this->fields |= (size_type)preds;
  }
  inline bool hasPredecessor(comp_type comp) const
  {
    return (this->fields & (1 << comp));
  }
  inline void addPredecessors(const PathNode& another)
  {
    this->fields |= another.predecessors();
  }

  // Order is the length of the kmer rank sequences representing the path label range.
  inline size_type order() const { return ((this->fields >> 8) & 0xFF); }
  inline void setOrder(size_type new_order)
  {
    this->fields &= ~(size_type)0xFF00;
    this->fields |= new_order << 8;
  }

  // LCP is the length of the common prefix of kmer rank sequences.
  inline size_type lcp() const { return ((this->fields >> 16) & 0xFF); }
  inline void setLCP(size_type new_lcp)
  {
    this->fields &= ~(size_type)0xFF0000;
    this->fields |= new_lcp << 16;
  }

  inline size_type pointer() const { return (this->fields >> 24); }
  inline void setPointer(size_type new_pointer)
  {
    this->fields &= 0xFFFFFF;
    this->fields |= new_pointer << 24;
  }

//------------------------------------------------------------------------------

  inline size_type ranks() const { return this->order() + 1; }
  inline size_type bytes() const { return sizeof(*this) + this->ranks() * sizeof(rank_type); }

  template<class ArrayType>
  inline rank_type firstLabel(size_type i, ArrayType& labels) const
  {
    return labels[this->pointer() + i];
  }

  template<class ArrayType>
  inline rank_type lastLabel(size_type i, ArrayType& labels) const
  {
    if(i < this->lcp()) { return labels[this->pointer() + i]; }
    else { return labels[this->pointer() + this->order()]; }
  }

//------------------------------------------------------------------------------

  PathNode(const KMer& kmer, WriteBuffer<rank_type>& labels);

  PathNode(const PathNode& source,
    const std::vector<rank_type>& old_labels, std::vector<rank_type>& new_labels);

  PathNode(const PathNode& left, const PathNode& right,
    const std::vector<rank_type>& old_labels, std::vector<rank_type>& new_labels);

  void print(std::ostream& out, const std::vector<rank_type>& labels) const;
  void print(std::ostream& out, const rank_type* labels) const;

  PathNode();
  explicit PathNode(std::vector<rank_type>& labels);
  PathNode(PathNode&& source) noexcept;
  ~PathNode();

  inline void swap(PathNode& another) noexcept
  {
    if(&another != this)
    {
      std::swap(this->from, another.from); std::swap(this->to, another.to);
      std::swap(this->fields, another.fields);
    }
  }

  PathNode& operator= (PathNode&& source) noexcept;

  /*
    These are dangerous, because the nodes will share the same label. Changing one will
    change the other as well.
  */
  PathNode(const PathNode& source);
  PathNode& operator= (const PathNode& source);
  void copy(const PathNode& source);
};

// Compares the first labels.
struct PathFirstComparator
{
  const std::vector<PathNode::rank_type>& labels;

  explicit PathFirstComparator(const std::vector<PathNode::rank_type>& _labels) : labels(_labels) { }

  inline bool operator() (const PathNode& a, const PathNode& b) const
  {
    size_type order = std::min(a.order(), b.order());
    for(size_type i = 0, a_ptr = a.pointer(), b_ptr = b.pointer(); i < order; i++, a_ptr++, b_ptr++)
    {
      if(labels[a_ptr] != labels[b_ptr]) { return (labels[a_ptr] < labels[b_ptr]); }
    }
    return (a.order() < b.order());
  }
};

// Compares the 'from' nodes.
struct PathFromComparator
{
  inline bool operator() (const PathNode& a, const PathNode& b) const
  {
    return (a.from < b.from);
  }
};

struct FromGetter
{
  inline static size_type get(const PathNode& path)
  {
    return path.from;
  }
};

//------------------------------------------------------------------------------

struct LCP
{
  typedef PathNode::rank_type rank_type;

  size_type       kmer_length, total_keys;
  sdsl::wt_blcd<> kmer_lcp; // Faster than proper RMQ for small values.

  LCP();
  LCP(const std::vector<key_type>& keys, size_type _kmer_length);

  // Disk-first equivalent of the vector constructor. key_name is a raw,
  // label-sorted unique key stream; only the final LCP construction buffer is
  // materialized.
  LCP(const std::string& key_name, size_type key_count,
    size_type _kmer_length, size_type buffer_bytes);

  /*
    Computes the minimal/maximal lcp of the path labels corresponding to path nodes a and b.
    a must be before b in lexicographic order, and the ranges must not overlap.
    The returned lcp value is a pair (x,y), where x is the lcp of the PathNode labels
    and y is the lcp of the first diverging kmers.
  */
  range_type min_lcp(const PathNode& a, const PathNode& b, const std::vector<rank_type>& labels) const;
  range_type max_lcp(const PathNode& a, const PathNode& b, const std::vector<rank_type>& labels) const;

  range_type min_lcp(const PathNode& a, const PathNode& b,
    const rank_type* a_labels, const rank_type* b_labels) const;
  range_type max_lcp(const PathNode& a, const PathNode& b,
    const rank_type* a_labels, const rank_type* b_labels) const;

  // Increments the lcp by 1.
  inline range_type increment(range_type lcp) const
  {
    lcp.second++;
    lcp.first += lcp.second / this->kmer_length;
    lcp.second %= this->kmer_length;
    return lcp;
  }

  void swap(LCP& another) noexcept;
};

//------------------------------------------------------------------------------

/*
  Diagnostics for the bounded merge used by prune() and MergedGraph. The
  counters are optional and do not affect construction semantics.
*/
struct PathGraphMergeStats
{
  size_type priority_spills, range_spills, from_set_sorts;
  size_type max_open_input_pairs;

  PathGraphMergeStats() :
    priority_spills(0), range_spills(0), from_set_sorts(0),
    max_open_input_pairs(0) { }
};

//------------------------------------------------------------------------------

/*
  A path graph is a set of files. Each file in the input graph becomes two temporary
  files: one for the paths and another for the rank sequences corresponding to path
  labels. The PathNodes in each file are sorted by their labels, and the read() member
  functions will also return the PathNodes in sorted order. The labels are stored in
  the same order as the PathNodes.
*/

struct PathGraph
{
  std::vector<std::string> path_names, rank_names;
  std::vector<size_type>   path_counts, rank_counts;

  // `file` in the legacy code is a physical stream index. It must never be
  // used as source-graph identity: one logical graph may have many spill
  // shards. These strong, parallel identifiers make that distinction explicit.
  std::vector<logical_file_id_t> logical_file_ids;
  std::vector<physical_shard_id_t> physical_shard_ids;

  size_type path_count, rank_count, range_count;
  size_type order, doubling_steps;

  size_type unique, redundant, unsorted, nondeterministic;

  bool delete_files;

  constexpr static size_type UNKNOWN = ~(size_type)0;
  const static std::string PREFIX;  // gcsa

  PathGraph(const InputGraph& source, sdsl::int_vector<0>& distinct_labels);
  PathGraph(size_type file_count, size_type path_order, size_type steps);
  PathGraph(const std::string& path_name, const std::string& rank_name);  // For debugging.
  ~PathGraph();

  void clear();
  void swap(PathGraph& another) noexcept;

  void open(std::ifstream& path_file, std::ifstream& rank_file, size_type file) const;

  inline size_type size() const { return this->path_count; }
  inline size_type ranks() const { return this->rank_count; }
  inline size_type ranges() const { return this->range_count; } // Only works after prune().
  inline size_type k() const { return this->order; }
  inline size_type step() const { return this->doubling_steps; }
  inline size_type files() const { return this->path_names.size(); }
  inline logical_file_id_t logicalFile(size_type file) const { return this->logical_file_ids.at(file); }
  inline physical_shard_id_t physicalShard(size_type file) const { return this->physical_shard_ids.at(file); }

  inline size_type bytes() const
  {
    return this->size() * sizeof(PathNode) + this->ranks() * sizeof(PathNode::rank_type);
  }

  /*
    The size limit (in bytes) is the space available for the new graph. The current graph is
    not taken into account, so you may want to use something like:

      path_graph.prune(lcp, total_size_limit - path_graph.bytes())
  */
  // Keep a bounded resident window while pruning one potentially very large
  // equal-label range. The default is deliberately small enough that the
  // external route cannot retain an adversarial range in RAM.
  void prune(const LCP& lcp, size_type size_limit,
    size_type group_buffer_bytes = MEGABYTE,
    PathGraphMergeStats* stats = nullptr);
  void extend(size_type size_limit, size_type memory_limit);

  void debugExtend();

  void read(std::vector<PathNode>& paths, std::vector<PathNode::rank_type>& labels, size_type file) const;

  PathGraph(const PathGraph&) = delete;
  PathGraph& operator= (const PathGraph&) = delete;
};

struct ExternalPathSortStats
{
  size_type runs, merge_operations, merge_passes, parallel_sorts;
  size_type max_records_resident, max_bytes_resident;
  size_type run_uncompressed_bytes, run_compressed_bytes;

  ExternalPathSortStats() :
    runs(0), merge_operations(0), merge_passes(0), parallel_sorts(0),
    max_records_resident(0), max_bytes_resident(0),
    run_uncompressed_bytes(0), run_compressed_bytes(0) { }
};

// Sort one physical PathGraph file without materializing it in memory.
size_type externalPathGraphSortMinimumBudget();
void externalPathGraphSort(PathGraph& graph, size_type file,
  size_type byte_budget, size_type fan_in, ExternalPathSortStats* stats = nullptr);

/*
  A bounded sink for records that are produced incrementally but must become a
  label-sorted PathGraph shard. Full buffers are sorted into immutable runs;
  finish() performs bounded multi-pass merging and installs the final
  path/rank pair. The producer never has to materialize an unsorted pair and
  then read it back through externalPathGraphSort().

  The target shard must be empty. committed_bytes accounts for other target
  shards under the same disk limit and is advanced only after finish().
*/
class ExternalPathSortSink
{
public:
  ExternalPathSortSink(PathGraph& graph, size_type file,
    size_type byte_budget, size_type fan_in, size_type size_limit,
    size_type& committed_bytes, ExternalPathSortStats* stats = nullptr);
  ~ExternalPathSortSink();

  void write(const PathNode& node, const PathNode::rank_type* labels);
  void finish();

  size_type paths() const;
  size_type ranks() const;
  size_type bytes() const;

private:
  struct Impl;
  Impl* impl;

  ExternalPathSortSink(const ExternalPathSortSink&);
  ExternalPathSortSink& operator=(const ExternalPathSortSink&);
};

struct ExternalPathJoinStats
{
  size_type left_records, right_records, sorted_bypass, generated_records;
  size_type initial_runs, merge_operations, blocked_key_groups, blocked_key_blocks;
  size_type join_parallel_sorts, label_sort_runs, label_merge_passes, label_parallel_sorts;
  size_type direct_label_records, intermediate_path_bytes_avoided;
  size_type join_partitions, worker_processes, restored_partitions;
  size_type recursive_splits, left_range_splits, right_range_splits;
  size_type max_records_resident, max_bytes_resident;

  ExternalPathJoinStats() :
    left_records(0), right_records(0), sorted_bypass(0), generated_records(0),
    initial_runs(0), merge_operations(0), blocked_key_groups(0), blocked_key_blocks(0),
    join_parallel_sorts(0), label_sort_runs(0), label_merge_passes(0), label_parallel_sorts(0),
    direct_label_records(0), intermediate_path_bytes_avoided(0),
    join_partitions(0), worker_processes(0), restored_partitions(0),
    recursive_splits(0), left_range_splits(0), right_range_splits(0),
    max_records_resident(0), max_bytes_resident(0) { }
};

// Prefix doubling with a bounded external sort-merge join. Input shards with
// the same logical ID are joined together, regardless of physical layout.
size_type externalPathJoinMinimumBudget();
void externalPathGraphExtend(PathGraph& graph, size_type size_limit,
  const ConstructionParameters& parameters, ExternalPathJoinStats* stats = nullptr,
  BuildWorkspace* workspace = nullptr, const std::string& checkpoint_task = std::string());

// Hidden child-process entry point used by build_gcsa and vg. The task file is
// versioned and contains immutable run ranges plus unique output paths.
int externalPathJoinWorker(const std::string& task_file);

//------------------------------------------------------------------------------

/*
  A merged graph is a path graph with the path nodes in one file, the labels in another
  file, and the additional start nodes in a third file. The PathNodes are sorted by their
  labels and their label pointers are set correctly.
*/

struct MergedGraph
{
  std::string path_name, rank_name, from_name, lcp_name;

  size_type path_count, rank_count, from_count;
  size_type order;

  std::vector<size_type> next;      // paths[next[comp]] is the first path starting with comp.
  std::vector<size_type> next_from; // Where to find the corresponding additional start nodes.

  constexpr static size_type UNKNOWN = ~(size_type)0;
  const static std::string PREFIX;  // gcsa

  /*
    The size limit (in bytes) is the space available for the new graph. The current graph is
    not taken into account, so you may want to use something like:

      MergedGraph merged_graph(source, mapper, kmer_lcp, total_size_limit - source.bytes())
  */
  MergedGraph(const PathGraph& source, const DeBruijnGraph& mapper,
    const LCP& kmer_lcp, size_type size_limit,
    size_type group_buffer_bytes = MEGABYTE,
    PathGraphMergeStats* stats = nullptr);
  ~MergedGraph();

  void clear();

  inline size_type size() const { return this->path_count; }
  inline size_type ranks() const { return this->rank_count; }
  inline size_type extra() const { return this->from_count; }
  inline size_type k() const { return this->order; }

  inline size_type path_bytes() const { return this->size() * sizeof(PathNode); }
  inline size_type rank_bytes() const { return this->ranks() * sizeof(PathNode::rank_type); }
  inline size_type from_bytes() const { return this->extra() * sizeof(range_type); }
  inline size_type lcp_bytes() const { return this->size(); }

  inline size_type bytes() const
  {
    return this->path_bytes() + this->rank_bytes() + this->from_bytes() + this->lcp_bytes();
  }

  MergedGraph(const MergedGraph&) = delete;
  MergedGraph& operator= (const MergedGraph&) = delete;
};

//------------------------------------------------------------------------------

} // namespace gcsa

#endif // GCSA_UTILS_H
