#ifndef GCSA_CHECKPOINT_H
#define GCSA_CHECKPOINT_H

#include <gcsa/path_graph.h>
#include <gcsa/workspace.h>

namespace gcsa
{

// Durable phase checkpoint for a complete PathGraph run set. Payloads are
// copied through a bounded byte buffer and retain semantic logical IDs.
bool pathGraphCheckpointExists(const BuildWorkspace& workspace,
  const std::string& task, const std::string& phase);
void checkpointPathGraph(BuildWorkspace& workspace, const PathGraph& graph,
  const std::string& task, const std::string& phase, size_type buffer_bytes);
void restorePathGraph(const BuildWorkspace& workspace, PathGraph& graph,
  const std::string& task, const std::string& phase, size_type buffer_bytes);

} // namespace gcsa

#endif // GCSA_CHECKPOINT_H
