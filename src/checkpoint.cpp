/*
  Copyright (c) 2026 Jouni Siren and GCSA2 contributors

  Durable PathGraph phase checkpoints.
*/

#include <gcsa/checkpoint.h>

#include <algorithm>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace gcsa
{

namespace
{

constexpr std::uint64_t PATH_GRAPH_CHECKPOINT_MAGIC = 0x3154504b434750ULL;
constexpr std::uint32_t PATH_GRAPH_CHECKPOINT_VERSION = 1;
constexpr size_type PATH_GRAPH_METADATA_FIXED = 8 + 4 + 10 * 8;
constexpr size_type PATH_GRAPH_METADATA_PER_FILE = 4 + 8 + 8 + 8;
constexpr size_type PATH_GRAPH_METADATA_LIMIT = 256 * MEGABYTE;

template<class Value>
void
appendLittle(std::vector<std::uint8_t>& target, Value value)
{
  for(size_type i = 0; i < sizeof(Value); i++)
  {
    target.push_back(static_cast<std::uint8_t>(value >> (8 * i)));
  }
}

template<class Value>
Value
readLittle(const std::vector<std::uint8_t>& source, size_type& offset)
{
  if(offset > source.size() || sizeof(Value) > source.size() - offset)
  {
    throw std::runtime_error("restorePathGraph(): truncated checkpoint metadata");
  }
  Value result = 0;
  for(size_type i = 0; i < sizeof(Value); i++)
  {
    result |= static_cast<Value>(source[offset + i]) << (8 * i);
  }
  offset += sizeof(Value); return result;
}

ArtifactIdentity
metadataIdentity(const std::string& task, const std::string& phase)
{
  return ArtifactIdentity(task, phase, "path-graph-metadata", "path-graph-v1");
}

ArtifactIdentity
pathIdentity(const std::string& task, const std::string& phase, size_type file)
{
  return ArtifactIdentity(task, phase, "paths-" + std::to_string(file), "path-nodes-v1");
}

ArtifactIdentity
rankIdentity(const std::string& task, const std::string& phase, size_type file)
{
  return ArtifactIdentity(task, phase, "ranks-" + std::to_string(file), "path-ranks-v1");
}

BuildWorkspace::ArtifactRef
copyIntoArtifact(BuildWorkspace& workspace, const ArtifactIdentity& identity,
  logical_file_id_t logical, physical_shard_id_t shard, const std::string& source,
  size_type records, size_type expected_bytes, size_type buffer_bytes,
  const std::string& sort_order)
{
  std::ifstream input;
  input.rdbuf()->pubsetbuf(nullptr, 0);
  input.open(source.c_str(), std::ios_base::binary);
  if(!input) { throw std::runtime_error("checkpointPathGraph(): cannot open " + source); }
  input.seekg(0, std::ios_base::end);
  size_type actual_bytes = static_cast<size_type>(input.tellg());
  input.seekg(0, std::ios_base::beg);
  if(actual_bytes != expected_bytes)
  {
    throw std::runtime_error("checkpointPathGraph(): file length does not match metadata: " + source);
  }
  BuildWorkspace::ArtifactWriter writer = workspace.open_artifact(
    identity, logical, shard, sort_order, "all");
  std::vector<std::uint8_t> buffer(std::max(static_cast<size_type>(1), buffer_bytes));
  size_type remaining = actual_bytes;
  while(remaining > 0)
  {
    size_type want = std::min(remaining, buffer.size());
    input.read(reinterpret_cast<char*>(buffer.data()), want);
    if(input.gcount() != static_cast<std::streamsize>(want))
    {
      throw std::runtime_error("checkpointPathGraph(): short read from " + source);
    }
    writer.write(buffer.data(), want); remaining -= want;
  }
  return writer.finish(records);
}

size_type
checkedPayloadBytes(size_type records, size_type record_bytes,
  const std::string& source)
{
  if(record_bytes != 0 && records > std::numeric_limits<size_type>::max() / record_bytes)
  {
    throw std::runtime_error(
      "checkpointPathGraph(): record count overflows payload size: " + source);
  }
  return records * record_bytes;
}

} // namespace

bool
pathGraphCheckpointExists(const BuildWorkspace& workspace,
  const std::string& task, const std::string& phase)
{
  return workspace.task_completed(task, phase);
}

void
checkpointPathGraph(BuildWorkspace& workspace, const PathGraph& graph,
  const std::string& task, const std::string& phase, size_type buffer_bytes)
{
  if(buffer_bytes == 0) { throw std::invalid_argument("checkpoint buffer must be nonzero"); }
  if(graph.files() > (PATH_GRAPH_METADATA_LIMIT - PATH_GRAPH_METADATA_FIXED) /
     PATH_GRAPH_METADATA_PER_FILE)
  {
    throw std::runtime_error("checkpointPathGraph(): shard manifest exceeds format limit");
  }
  std::vector<std::uint8_t> metadata;
  metadata.reserve(PATH_GRAPH_METADATA_FIXED + graph.files() * PATH_GRAPH_METADATA_PER_FILE);
  appendLittle<std::uint64_t>(metadata, PATH_GRAPH_CHECKPOINT_MAGIC);
  appendLittle<std::uint32_t>(metadata, PATH_GRAPH_CHECKPOINT_VERSION);
  appendLittle<std::uint64_t>(metadata, graph.files());
  appendLittle<std::uint64_t>(metadata, graph.k());
  appendLittle<std::uint64_t>(metadata, graph.step());
  appendLittle<std::uint64_t>(metadata, graph.size());
  appendLittle<std::uint64_t>(metadata, graph.ranks());
  appendLittle<std::uint64_t>(metadata, graph.ranges());
  appendLittle<std::uint64_t>(metadata, graph.unique);
  appendLittle<std::uint64_t>(metadata, graph.redundant);
  appendLittle<std::uint64_t>(metadata, graph.unsorted);
  appendLittle<std::uint64_t>(metadata, graph.nondeterministic);
  for(size_type file = 0; file < graph.files(); file++)
  {
    appendLittle<std::uint32_t>(metadata, graph.logicalFile(file).value);
    appendLittle<std::uint64_t>(metadata, graph.physicalShard(file).value);
    appendLittle<std::uint64_t>(metadata, graph.path_counts[file]);
    appendLittle<std::uint64_t>(metadata, graph.rank_counts[file]);
  }

  std::vector<BuildWorkspace::ArtifactRef> artifacts;
  artifacts.reserve(1 + 2 * graph.files());
  ArtifactIdentity metadata_id = metadataIdentity(task, phase);
  BuildWorkspace::ArtifactWriter metadata_writer = workspace.open_artifact(
    metadata_id, logical_file_id_t(0), physical_shard_id_t(0), "manifest", "all");
  metadata_writer.write(metadata.data(), metadata.size());
  artifacts.push_back(metadata_writer.finish(graph.files()));
  for(size_type file = 0; file < graph.files(); file++)
  {
    artifacts.push_back(copyIntoArtifact(workspace, pathIdentity(task, phase, file),
      graph.logicalFile(file), graph.physicalShard(file), graph.path_names[file],
      graph.path_counts[file], checkedPayloadBytes(graph.path_counts[file],
        sizeof(PathNode), graph.path_names[file]),
      buffer_bytes, "label"));
    artifacts.push_back(copyIntoArtifact(workspace, rankIdentity(task, phase, file),
      graph.logicalFile(file), graph.physicalShard(file), graph.rank_names[file],
      graph.rank_counts[file], checkedPayloadBytes(graph.rank_counts[file],
        sizeof(PathNode::rank_type), graph.rank_names[file]),
      buffer_bytes, "path-order"));
  }
  workspace.commit_task(task, phase, artifacts);
}

void
restorePathGraph(const BuildWorkspace& workspace, PathGraph& graph,
  const std::string& task, const std::string& phase, size_type buffer_bytes)
{
  if(!workspace.task_completed(task, phase))
  {
    throw std::runtime_error("restorePathGraph(): phase is not committed: " + task + "/" + phase);
  }
  std::vector<std::uint8_t> metadata = workspace.read_artifact_payload(
    metadataIdentity(task, phase), logical_file_id_t(0), physical_shard_id_t(0),
    PATH_GRAPH_METADATA_LIMIT);
  size_type offset = 0;
  if(readLittle<std::uint64_t>(metadata, offset) != PATH_GRAPH_CHECKPOINT_MAGIC ||
     readLittle<std::uint32_t>(metadata, offset) != PATH_GRAPH_CHECKPOINT_VERSION)
  {
    throw std::runtime_error("restorePathGraph(): incompatible checkpoint metadata");
  }
  size_type files = readLittle<std::uint64_t>(metadata, offset);
  size_type order = readLittle<std::uint64_t>(metadata, offset);
  size_type step = readLittle<std::uint64_t>(metadata, offset);
  size_type total_paths = readLittle<std::uint64_t>(metadata, offset);
  size_type total_ranks = readLittle<std::uint64_t>(metadata, offset);
  size_type ranges = readLittle<std::uint64_t>(metadata, offset);
  size_type unique = readLittle<std::uint64_t>(metadata, offset);
  size_type redundant = readLittle<std::uint64_t>(metadata, offset);
  size_type unsorted = readLittle<std::uint64_t>(metadata, offset);
  size_type nondeterministic = readLittle<std::uint64_t>(metadata, offset);
  if(files > (metadata.size() - offset) / PATH_GRAPH_METADATA_PER_FILE ||
     offset + files * PATH_GRAPH_METADATA_PER_FILE != metadata.size())
  {
    throw std::runtime_error("restorePathGraph(): invalid checkpoint shard count");
  }

  PathGraph restored(files, order, step);
  restored.path_count = 0; restored.rank_count = 0;
  restored.range_count = ranges; restored.unique = unique; restored.redundant = redundant;
  restored.unsorted = unsorted; restored.nondeterministic = nondeterministic;
  for(size_type file = 0; file < files; file++)
  {
    restored.logical_file_ids[file] = logical_file_id_t(
      readLittle<std::uint32_t>(metadata, offset));
    restored.physical_shard_ids[file] = physical_shard_id_t(
      readLittle<std::uint64_t>(metadata, offset));
    restored.path_counts[file] = readLittle<std::uint64_t>(metadata, offset);
    restored.rank_counts[file] = readLittle<std::uint64_t>(metadata, offset);
    restored.path_count += restored.path_counts[file];
    restored.rank_count += restored.rank_counts[file];
    workspace.restore_artifact(pathIdentity(task, phase, file),
      restored.logicalFile(file), restored.physicalShard(file),
      restored.path_names[file], buffer_bytes);
    workspace.restore_artifact(rankIdentity(task, phase, file),
      restored.logicalFile(file), restored.physicalShard(file),
      restored.rank_names[file], buffer_bytes);
  }
  if(restored.path_count != total_paths || restored.rank_count != total_ranks)
  {
    throw std::runtime_error("restorePathGraph(): checkpoint totals do not match shards");
  }
  graph.clear(); graph.swap(restored);
}

} // namespace gcsa
