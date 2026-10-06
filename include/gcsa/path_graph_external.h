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

namespace gcsa
{

//------------------------------------------------------------------------------

struct ExternalPathSortStats
{
  size_type runs, merge_passes, max_records_resident, max_bytes_resident;

  ExternalPathSortStats() :
    runs(0), merge_passes(0), max_records_resident(0), max_bytes_resident(0) { }
};

// Sort one physical PathGraph file without materializing it in memory.
size_type externalPathGraphSortMinimumBudget();
void externalPathGraphSort(PathGraph& graph, size_type file,
  size_type byte_budget, size_type fan_in, ExternalPathSortStats* stats = nullptr);

class ExternalPathSortSink
{
public:
  ExternalPathSortSink(PathGraph& graph, size_type file,
    size_type byte_budget, size_type fan_in, size_type size_limit,
    size_type& committed_bytes, ExternalPathSortStats* stats = nullptr);
  ~ExternalPathSortSink();

  void write(const PathNode& node, const PathNode::rank_type* labels);
  void finish();

private:
  struct Impl;
  Impl* impl;

  ExternalPathSortSink(const ExternalPathSortSink&) = delete;
  ExternalPathSortSink& operator=(const ExternalPathSortSink&) = delete;
};

// Merge workspace and pair counts used by prune() and MergedGraph.
size_type pathMergeInputBudget(const ConstructionParameters& parameters);
size_type pathMergeOutputPairs(size_type max_open_files, size_type logical_inputs);
size_type pathMergeInputPairs(size_type max_open_files, size_type logical_inputs);

struct ExternalPathJoinStats
{
  size_type left_records, right_records, sorted_bypass, generated_records;
  size_type initial_runs, merge_operations, blocked_key_groups;
  size_type max_records_resident, max_bytes_resident;

  ExternalPathJoinStats() :
    left_records(0), right_records(0), sorted_bypass(0), generated_records(0),
    initial_runs(0), merge_operations(0), blocked_key_groups(0),
    max_records_resident(0), max_bytes_resident(0) { }
};

// Prune independent zero-LCP roots using a privately admitted thread limit.
// Falls back to the serial route when complete worker caches do not fit.
void externalPathGraphPrune(PathGraph& graph, const LCP& lcp, size_type size_limit,
  size_type group_buffer_bytes, size_type max_open_files,
  size_type requested_workers, PathGraphMergeStats* stats = nullptr);

// Prefix doubling with a bounded external sort-merge join. Input shards with
// the same logical ID are joined together, regardless of physical layout.
size_type externalPathJoinMinimumBudget();
void externalPathGraphExtend(PathGraph& graph, size_type size_limit,
  const ConstructionParameters& parameters, ExternalPathJoinStats* stats = nullptr);
size_type externalPathGraphExtendWithThreads(PathGraph& graph,
  size_type size_limit, const ConstructionParameters& parameters,
  size_type requested_threads, ExternalPathJoinStats* stats = nullptr);

//------------------------------------------------------------------------------

} // namespace gcsa

#endif // GCSA_PATH_GRAPH_EXTERNAL_H
