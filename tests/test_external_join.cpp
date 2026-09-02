#include <gcsa/path_graph.h>
#include <gcsa/checkpoint.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <filesystem>
#include <limits>
#include <string>
#include <unistd.h>
#include <vector>

using namespace gcsa;

namespace
{

struct TestRecord
{
  node_type from, to;
  byte_type predecessors;
  PathNode::rank_type label;
  bool sorted;
};

struct SemanticRecord
{
  node_type from, to;
  byte_type predecessors;
  size_type order, lcp;
  std::vector<PathNode::rank_type> labels;

  bool operator<(const SemanticRecord& another) const
  {
    if(labels != another.labels) { return labels < another.labels; }
    if(from != another.from) { return from < another.from; }
    if(to != another.to) { return to < another.to; }
    if(predecessors != another.predecessors) { return predecessors < another.predecessors; }
    if(order != another.order) { return order < another.order; }
    return lcp < another.lcp;
  }

  bool operator==(const SemanticRecord& another) const
  {
    return from == another.from && to == another.to &&
      predecessors == another.predecessors && order == another.order &&
      lcp == another.lcp && labels == another.labels;
  }
};

void require(bool condition)
{
  if(!condition) { std::abort(); }
}

void
writePathPair(const std::string& path_name, const std::string& rank_name,
  std::vector<TestRecord> records)
{
  std::sort(records.begin(), records.end(), [](const TestRecord& left, const TestRecord& right)
  {
    return left.label < right.label;
  });
  std::ofstream paths(path_name.c_str(), std::ios_base::binary);
  std::ofstream ranks(rank_name.c_str(), std::ios_base::binary);
  require(paths && ranks);
  size_type pointer = 0;
  for(const TestRecord& record : records)
  {
    PathNode node;
    node.from = record.from; node.to = record.to; node.fields = 0;
    if(record.sorted) { node.makeSorted(); }
    node.setPredecessors(record.predecessors);
    node.setOrder(1); node.setLCP(1); node.setPointer(pointer);
    PathNode::rank_type label[2] = { record.label, 0 };
    paths.write(reinterpret_cast<const char*>(&node), sizeof(node));
    ranks.write(reinterpret_cast<const char*>(label), sizeof(label));
    pointer += 2;
  }
}

void
appendShard(PathGraph& graph, const std::string& path_name,
  const std::string& rank_name, logical_file_id_t logical,
  physical_shard_id_t physical)
{
  std::ifstream paths(path_name.c_str(), std::ios_base::binary);
  std::ifstream ranks(rank_name.c_str(), std::ios_base::binary);
  require(paths && ranks);
  paths.seekg(0, std::ios_base::end); ranks.seekg(0, std::ios_base::end);
  size_type path_count = static_cast<size_type>(paths.tellg()) / sizeof(PathNode);
  size_type rank_count = static_cast<size_type>(ranks.tellg()) / sizeof(PathNode::rank_type);
  graph.path_names.push_back(path_name); graph.rank_names.push_back(rank_name);
  graph.path_counts.push_back(path_count); graph.rank_counts.push_back(rank_count);
  graph.logical_file_ids.push_back(logical); graph.physical_shard_ids.push_back(physical);
  graph.path_count += path_count; graph.rank_count += rank_count;
}

std::vector<SemanticRecord>
readGraph(const PathGraph& graph)
{
  std::vector<SemanticRecord> result;
  for(size_type file = 0; file < graph.files(); file++)
  {
    std::vector<PathNode> paths;
    std::vector<PathNode::rank_type> ranks;
    graph.read(paths, ranks, file);
    for(const PathNode& path : paths)
    {
      SemanticRecord record;
      record.from = path.from; record.to = path.to;
      record.predecessors = path.predecessors();
      record.order = path.order(); record.lcp = path.lcp();
      for(size_type i = 0; i < path.ranks(); i++)
      {
        record.labels.push_back(ranks[path.pointer() + i]);
      }
      result.push_back(record);
    }
  }
  std::sort(result.begin(), result.end());
  return result;
}

} // namespace

int main()
{
  omp_set_num_threads(1);
  std::vector<TestRecord> left =
  {
    { 1, 100, 1, 10, false },
    { 2, 100, 2, 20, false },
    { 3, 100, 4, 30, false }
  };
  std::vector<TestRecord> right;
  // Exceed both the minimum-budget run buffer and the two-reader merge cache
  // so this remains a forced multi-run, pathological-key spill test even when
  // phase I/O is buffered in bytes.
  for(size_type i = 0; i < 1200; i++)
  {
    right.push_back({ 100, 1000 + i, static_cast<byte_type>(1),
      static_cast<PathNode::rank_type>(1000 + i), false });
  }
  // A sorted path bypasses left expansion, but remains a valid right target.
  right.push_back({ 100, std::numeric_limits<node_type>::max(), 8, 5000, true });

  std::vector<TestRecord> combined = left;
  combined.insert(combined.end(), right.begin(), right.end());
  std::string base = "test_external_join_" +
    std::to_string(static_cast<unsigned long long>(getpid()));
  std::string combined_path = base + ".combined.path";
  std::string combined_rank = base + ".combined.rank";
  std::string left_path = base + ".left.path";
  std::string left_rank = base + ".left.rank";
  std::string right_path = base + ".right.path";
  std::string right_rank = base + ".right.rank";
  writePathPair(combined_path, combined_rank, combined);
  writePathPair(left_path, left_rank, left);
  writePathPair(right_path, right_rank, right);

  PathGraph legacy(combined_path, combined_rank);
  legacy.order = 1;
  legacy.extend(GIGABYTE, 64 * MEGABYTE);

  PathGraph external(left_path, left_rank);
  external.order = 1;
  external.logical_file_ids[0] = logical_file_id_t(7);
  external.physical_shard_ids[0] = physical_shard_id_t(101);
  appendShard(external, right_path, right_rank,
    logical_file_id_t(7), physical_shard_id_t(202));

  ConstructionParameters parameters;
  size_type memory_budget = std::max(externalPathJoinMinimumBudget(),
    externalPathGraphSortMinimumBudget());
  parameters.setMemoryLimitBytes(memory_budget);
  parameters.setJoinPartitionSize(externalPathJoinMinimumBudget());
  parameters.setSortRunSize(externalPathGraphSortMinimumBudget());
  parameters.setMergeFanIn(2);
  parameters.setWorkDirectory(".");
  ExternalPathJoinStats stats;
  externalPathGraphExtend(external, GIGABYTE, parameters, &stats);

  require(readGraph(legacy) == readGraph(external));
  require(external.files() == 1);
  require(external.logicalFile(0) == logical_file_id_t(7));
  require(stats.generated_records == left.size() * right.size());
  require(stats.sorted_bypass == 1);
  require(stats.initial_runs > 2 && stats.merge_operations > 0);
  require(stats.blocked_key_groups > 0);
  require(stats.max_bytes_resident <= memory_budget);

  char workspace_root[] = "/tmp/gcsa-path-checkpoint-XXXXXX";
  require(mkdtemp(workspace_root) != nullptr);
  BuildWorkspace::Settings semantic;
  semantic["fixture"] = "external-join";
  BuildWorkspace workspace(workspace_root, semantic, BuildWorkspace::Settings(),
    BuildWorkspace::NEW_WORKSPACE);
  checkpointPathGraph(workspace, external, "step-01", "extend", 4096);
  require(pathGraphCheckpointExists(workspace, "step-01", "extend"));
  PathGraph restored(0, 1, 0);
  restorePathGraph(workspace, restored, "step-01", "extend", 4096);
  require(readGraph(restored) == readGraph(external));
  require(restored.logicalFile(0) == logical_file_id_t(7));

  // Physical shards with different logical IDs must never join.
  PathGraph separated(left_path, left_rank);
  separated.order = 1;
  separated.logical_file_ids[0] = logical_file_id_t(1);
  separated.physical_shard_ids[0] = physical_shard_id_t(101);
  appendShard(separated, right_path, right_rank,
    logical_file_id_t(2), physical_shard_id_t(202));
  // Exercise the optional full payload-checksum pass as well as the default
  // header/footer/length validation used above.
  parameters.setVerifyWorkspace();
  externalPathGraphExtend(separated, GIGABYTE, parameters);
  require(separated.files() == 2);
  require(separated.size() == 1); // Only the already-sorted bypass path remains.

  std::remove(combined_path.c_str()); std::remove(combined_rank.c_str());
  std::remove(left_path.c_str()); std::remove(left_rank.c_str());
  std::remove(right_path.c_str()); std::remove(right_rank.c_str());
  std::filesystem::remove_all(workspace_root);
  return 0;
}
