/*
  Copyright (c) 2026 GCSA2 contributors

  Fork-owned declarations for the external-memory path-graph route.

  These are separated from path_graph.h so that upstream's header stays
  recognisable: a reader opening path_graph.h to understand PathNode, PathGraph
  or MergedGraph meets those types, not the fork's budget arithmetic. Nothing
  here is referenced by the legacy in-memory route, and none of it needs access
  to PathNode's internals -- it is free functions, their statistics structs, and
  one sink class, all of which take PathGraph by reference.

  Include this header explicitly from the translation units that build or drive
  the external route. path_graph.h deliberately does not include it.
*/

#ifndef GCSA_PATH_GRAPH_EXTERNAL_H
#define GCSA_PATH_GRAPH_EXTERNAL_H

#include <gcsa/path_graph.h>
#include <gcsa/resources.h>

namespace gcsa
{

//------------------------------------------------------------------------------

struct ExternalPathSortStats
{
  size_type runs, merge_operations, merge_passes, parallel_sorts;
  size_type max_records_resident, max_bytes_resident;
  size_type run_uncompressed_bytes, run_compressed_bytes;
  size_type grouped_records, group_headers, context_bytes_saved;
  // Codec counters cover every immutable run written, including merge output;
  // they measure transient I/O avoided, not distinct final graph records.

  ExternalPathSortStats() :
    runs(0), merge_operations(0), merge_passes(0), parallel_sorts(0),
    max_records_resident(0), max_bytes_resident(0),
    run_uncompressed_bytes(0), run_compressed_bytes(0),
    grouped_records(0), group_headers(0), context_bytes_saved(0) { }
};

// Sort one physical PathGraph file without materializing it in memory.
size_type externalPathGraphSortMinimumBudget();
void externalPathGraphSort(PathGraph& graph, size_type file,
  size_type byte_budget, size_type fan_in, ExternalPathSortStats* stats = nullptr);

// Worst-case live filesystem bytes for the final path/rank pair produced by a
// streaming label sorter. The bound includes framed commit metadata and its
// temporary index copy; raw output is exactly its logical payload size.
size_type externalPathGraphShardPeakBytes(size_type paths, size_type ranks,
  size_type sort_byte_budget, const TempFileCodecParameters& codec);

// Largest framed block whose decoded working set still lets `pairs` shard pairs
// stay open at once inside `byte_budget`.
//
// The label merge in prune() and MergedGraph reads every open shard
// sequentially, so a block within this bound is decoded exactly once per pass.
// A larger one forces the input cache to evict and reconstruct readers, and
// each miss then decodes whole blocks to deliver one 24-byte PathNode. Choosing
// a block against the label-sort budget alone is what committed shards the
// merger could not open at all.
size_type mergeAdmissibleBlockSize(size_type byte_budget, size_type pairs,
  size_type requested_block);

// Merge workspace, and the framed pair counts prune() and MergedGraph divide
// their --max-open-files allowance into for `logical_inputs` logical inputs.
// Producers of generated shards bound their compression block against these,
// and the post-join compactor retains no more shards than the input count, so
// what one phase commits the next can open: both merges visit their shards in
// round-robin label order, and a cache one entry short of the shard count
// misses on every record.
size_type pathMergeInputBudget(const ConstructionParameters& parameters);
size_type pathMergeOutputPairs(size_type max_open_files, size_type logical_inputs);
size_type pathMergeInputPairs(size_type max_open_files, size_type logical_inputs);
size_type pathMergeInputPairs(const ConstructionParameters& parameters);

// The ceiling that workspace may grow to: a sixteenth of --memory-limit, the
// share this buffer has always been bounded by. Producers size their block
// against this rather than against --io-buffer-size, which sizes stream
// buffers and which no entry point above this library exposes.
size_type pathMergeCeilingBudget(const ConstructionParameters& parameters);

// Reservation for decoded framed input blocks while merging `source`. Every
// shard a merge opens stays open for the whole pass and they are visited in
// round-robin label order, so a cache one entry short of the shard count
// misses on every record: resident shard pairs are funded first. Any room
// left under the independent quarter-memory ceiling is available to the
// shared read/decode prefetch pool; prefetch is disabled rather than evicting
// a resident pair when there is no surplus. The cache is deliberately
// separate from the equal-label group buffer, which scales with one label
// range rather than with the number of shards.
size_type pathMergeInputCacheBudget(const ConstructionParameters& parameters,
  const PathGraph& source);

// The same workspace for a generation that already exists. A committed shard
// declares its own block size and nothing rewrites it, so this raises the
// buffer to hold one decoded path/rank pair of source. --memory-limit remains
// the ceiling: a workspace whose block needs more than a sixteenth of it is
// reported with the limit it does need, not silently over-allocated. A raw
// (unframed) generation gets that sixteenth as well, because the same budget
// sizes the merger's spillable equal-label groups and its raw read windows.
size_type pathMergeInputBudget(const ConstructionParameters& parameters,
  const PathGraph& source);

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
    size_type& committed_bytes, ExternalPathSortStats* stats = nullptr,
    const TempFileCodecParameters& codec = TempFileCodecParameters());
  ~ExternalPathSortSink();

  void write(const PathNode& node, const PathNode::rank_type* labels);
  void finish();

  size_type paths() const;
  size_type ranks() const;
  size_type bytes() const;
  size_type storedBytes() const;

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
  size_type sampled_plan_records, radix_plan_bins, radix_plan_splits;
  size_type radix_plan_max_bits, radix_plan_capped;
  size_type radix_boundary_flushes, restored_radix_plans;
  // Exact process-worker planning streams these compact sidecars. The full
  // JoinRecord counter must remain zero: the 4096-record sample is the only
  // planner read of the wide run payload.
  size_type sidecar_plan_groups, sidecar_plan_detail_records;
  size_type full_record_plan_rescans, restored_sidecar_metadata;
  // These aggregate PathSortRun codec references across all label-run passes.
  // A record rewritten by a merge may therefore contribute more than once.
  size_type grouped_expansion_records, expansion_context_bytes_saved;
  // Join-run byte totals include every immutable primary run generation,
  // including merge rewrites. Logical bytes describe the unchanged
  // fixed-record stream. Sidecar totals include every immutable group/detail
  // pair written alongside those runs, in their logical and installed forms.
  size_type join_run_logical_bytes, join_run_stored_bytes, compressed_join_runs;
  size_type join_sidecar_logical_bytes, join_sidecar_stored_bytes;
  size_type compressed_join_sidecars;
  // Top-level operational allocations. Distribution has its own lifetime;
  // label sorting and join blocking are concurrent and must sum to the limit.
  size_type distribution_sort_budget, label_sort_budget, join_block_budget;
  size_type max_records_resident, max_bytes_resident;
  // Post-join compaction: merged batches, and the most that ran at once.
  size_type compaction_batches, compaction_concurrency;
  // Key-range distribution: the most ranges one logical input was split into,
  // the most range workers that ran at once, and range plans restored from
  // the workspace.
  size_type distribution_ranges, distribution_concurrency, restored_range_plans;

  ExternalPathJoinStats() :
    left_records(0), right_records(0), sorted_bypass(0), generated_records(0),
    initial_runs(0), merge_operations(0), blocked_key_groups(0), blocked_key_blocks(0),
    join_parallel_sorts(0), label_sort_runs(0), label_merge_passes(0), label_parallel_sorts(0),
    direct_label_records(0), intermediate_path_bytes_avoided(0),
    join_partitions(0), worker_processes(0), restored_partitions(0),
    recursive_splits(0), left_range_splits(0), right_range_splits(0),
    sampled_plan_records(0), radix_plan_bins(0), radix_plan_splits(0),
    radix_plan_max_bits(0), radix_plan_capped(0),
    radix_boundary_flushes(0), restored_radix_plans(0),
    sidecar_plan_groups(0), sidecar_plan_detail_records(0),
    full_record_plan_rescans(0), restored_sidecar_metadata(0),
    grouped_expansion_records(0), expansion_context_bytes_saved(0),
    join_run_logical_bytes(0), join_run_stored_bytes(0), compressed_join_runs(0),
    join_sidecar_logical_bytes(0), join_sidecar_stored_bytes(0),
    compressed_join_sidecars(0),
    distribution_sort_budget(0), label_sort_budget(0), join_block_budget(0),
    max_records_resident(0), max_bytes_resident(0),
    compaction_batches(0), compaction_concurrency(0),
    distribution_ranges(0), distribution_concurrency(0), restored_range_plans(0) { }
};

// Prefix doubling with a bounded external sort-merge join. Input shards with
// the same logical ID are joined together, regardless of physical layout.
size_type externalPathJoinMinimumBudget();
void externalPathGraphExtend(PathGraph& graph, size_type size_limit,
  const ConstructionParameters& parameters, ExternalPathJoinStats* stats = nullptr,
  BuildWorkspace* workspace = nullptr, const std::string& checkpoint_task = std::string());

// Collapses the join workers' physical shards per logical input to at most the
// pairs the next merge holds open. externalPathGraphExtend() calls it after the
// joins; it is declared for the tests of its concurrent batches.
void compactLogicalJoinShards(PathGraph& source, size_type size_limit,
  const ConstructionParameters& parameters, size_type label_fan_in,
  MemoryBudget& memory, size_type& committed_bytes,
  ExternalPathJoinStats* stats);

// Hidden child-process entry point used by build_gcsa and vg. The task file is
// versioned and contains immutable run ranges plus unique output paths.
int externalPathJoinWorker(const std::string& task_file);

//------------------------------------------------------------------------------

} // namespace gcsa

#endif // GCSA_PATH_GRAPH_EXTERNAL_H
