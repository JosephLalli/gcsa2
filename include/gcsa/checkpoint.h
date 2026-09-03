#ifndef GCSA_CHECKPOINT_H
#define GCSA_CHECKPOINT_H

#include <gcsa/path_graph.h>
#include <gcsa/workspace.h>

namespace gcsa
{

// Durable phase checkpoint for a complete PathGraph run set. Same-filesystem
// payloads are adopted by immutable hardlink; other filesystems use a bounded
// copy. Physical storage never changes semantic logical-file identity.
bool pathGraphCheckpointExists(const BuildWorkspace& workspace,
  const std::string& task, const std::string& phase);
void checkpointPathGraph(BuildWorkspace& workspace, const PathGraph& graph,
  const std::string& task, const std::string& phase, size_type buffer_bytes);
void restorePathGraph(const BuildWorkspace& workspace, PathGraph& graph,
  const std::string& task, const std::string& phase, size_type buffer_bytes);
void restorePathGraph(const BuildWorkspace& workspace, PathGraph& graph,
  const std::string& task, const std::string& phase, size_type buffer_bytes,
  bool verify_checksum);

} // namespace gcsa

#endif // GCSA_CHECKPOINT_H
