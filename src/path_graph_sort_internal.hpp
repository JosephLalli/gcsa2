#ifndef GCSA_PATH_GRAPH_SORT_INTERNAL_HPP
#define GCSA_PATH_GRAPH_SORT_INTERNAL_HPP

#include <gcsa/path_graph_external.h>

namespace gcsa
{

// Internal construction seam. The caller supplies the thread count captured
// before external construction temporarily serializes ambient OpenMP work.
void externalPathGraphSortWithThreads(PathGraph& graph, size_type file,
  size_type byte_budget, size_type fan_in, ExternalPathSortStats* stats,
  size_type available_threads);

} // namespace gcsa

#endif // GCSA_PATH_GRAPH_SORT_INTERNAL_HPP
