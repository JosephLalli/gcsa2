#include <gcsa/path_graph.h>
#include <gcsa/path_graph_external.h>
#include <gcsa/checkpoint.h>
#include <gcsa/compressed_block.h>

#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iterator>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <fcntl.h>
#include <map>
#include <sys/resource.h>
#include <sys/wait.h>
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
  size_type order = 1;
  size_type lcp = 1;
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

void
requireAt(bool condition, const char* expression, int line)
{
  if(!condition)
  {
    std::fprintf(stderr, "test_external_join:%d: requirement failed: %s\n",
      line, expression);
    std::abort();
  }
}

// Keep assertions useful under the production -O3 test build, where adjacent
// abort sites can otherwise collapse onto one misleading source line.
#define require(value) requireAt((value), #value, __LINE__)

void
requireRecordedChecksum(const ClosedPayloadChecksum& provenance, const std::string& path)
{
  std::ifstream input(path.c_str(), std::ios::binary);
  require(bool(input));
  std::vector<char> payload((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
  require(provenance.matches(path, payload.size()));
  require(provenance.value == BuildWorkspace::checksum(payload.data(), payload.size()));
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
    node.setOrder(record.order); node.setLCP(record.lcp); node.setPointer(pointer);
    PathNode::rank_type label[2] = { record.label, 0 };
    paths.write(reinterpret_cast<const char*>(&node), sizeof(node));
    ranks.write(reinterpret_cast<const char*>(label),
      node.ranks() * sizeof(PathNode::rank_type));
    pointer += node.ranks();
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

std::vector<std::uint8_t>
readBytes(const std::string& name)
{
  std::ifstream input(name.c_str(), std::ios_base::binary);
  require(static_cast<bool>(input));
  input.seekg(0, std::ios_base::end);
  std::streamoff end = input.tellg();
  require(end >= 0);
  input.seekg(0, std::ios_base::beg);
  std::vector<std::uint8_t> result(static_cast<size_type>(end));
  if(!result.empty())
  {
    input.read(reinterpret_cast<char*>(result.data()), result.size());
    require(input.gcount() == static_cast<std::streamsize>(result.size()));
  }
  return result;
}

void
requireByteIdenticalGraph(const PathGraph& expected, const PathGraph& actual)
{
  require(expected.files() == actual.files());
  for(size_type file = 0; file < expected.files(); file++)
  {
    require(expected.logicalFile(file) == actual.logicalFile(file));
    require(expected.physicalShard(file) == actual.physicalShard(file));
    require(expected.path_counts[file] == actual.path_counts[file]);
    require(expected.rank_counts[file] == actual.rank_counts[file]);
    require(readBytes(expected.path_names[file]) == readBytes(actual.path_names[file]));
    require(readBytes(expected.rank_names[file]) == readBytes(actual.rank_names[file]));
  }
}

void
requireSameGraph(const PathGraph& expected_graph, const PathGraph& actual_graph,
  const char* context)
{
  std::vector<SemanticRecord> expected = readGraph(expected_graph);
  std::vector<SemanticRecord> actual = readGraph(actual_graph);
  if(expected != actual)
  {
    size_type mismatch = 0;
    while(mismatch < expected.size() && mismatch < actual.size() &&
          expected[mismatch] == actual[mismatch]) { mismatch++; }
    std::cerr << context << ": semantic graph mismatch at " << mismatch
              << " (expected " << expected.size() << " records, got "
              << actual.size() << ")" << std::endl;
    if(mismatch < expected.size() && mismatch < actual.size())
    {
      std::cerr << "expected from/to/order/lcp/labels: " << expected[mismatch].from
                << "/" << expected[mismatch].to << "/" << expected[mismatch].order
                << "/" << expected[mismatch].lcp << "/" << expected[mismatch].labels.size()
                << "; got " << actual[mismatch].from << "/" << actual[mismatch].to
                << "/" << actual[mismatch].order << "/" << actual[mismatch].lcp
                << "/" << actual[mismatch].labels.size() << std::endl;
    }
    std::abort();
  }
}

// Every workspace entry whose name starts with `prefix`: task completion
// records by name, artifacts by name and bytes.
std::map<std::string, std::vector<std::uint8_t>>
workspaceEntries(const std::string& root, const std::string& prefix)
{
  std::map<std::string, std::vector<std::uint8_t>> result;
  for(const auto& entry : std::filesystem::directory_iterator(root))
  {
    const std::string name = entry.path().filename().string();
    if(name.compare(0, prefix.size(), prefix) != 0) { continue; }
    result[name] = (entry.path().extension() == ".bin" ?
      readBytes(entry.path().string()) : std::vector<std::uint8_t>());
  }
  return result;
}

// Test-only fault injection for the re-executed worker: when
// GCSA_TEST_JOIN_WORKER_FAIL_DIR names a directory, each worker claims the next
// ordinal by creating a file there, and the one whose ordinal equals
// GCSA_TEST_JOIN_WORKER_FAIL_AT exits before joining, as a crashed worker would.
bool
injectedWorkerFailure()
{
  const char* directory = std::getenv("GCSA_TEST_JOIN_WORKER_FAIL_DIR");
  const char* at = std::getenv("GCSA_TEST_JOIN_WORKER_FAIL_AT");
  if(directory == nullptr || at == nullptr) { return false; }
  for(size_type ordinal = 0; ; ordinal++)
  {
    const std::string claim = std::string(directory) + "/" + std::to_string(ordinal);
    int descriptor = ::open(claim.c_str(), O_CREAT | O_EXCL | O_WRONLY, 0644);
    if(descriptor < 0)
    {
      if(errno == EEXIST) { continue; }
      return false;
    }
    ::close(descriptor);
    return ordinal == std::strtoull(at, nullptr, 10);
  }
}

bool
labelSorted(const PathGraph& graph)
{
  for(size_type file = 0; file < graph.files(); file++)
  {
    std::vector<PathNode> paths;
    std::vector<PathNode::rank_type> ranks;
    graph.read(paths, ranks, file);
    for(size_type i = 1; i < paths.size(); i++)
    {
      const PathNode& left = paths[i - 1]; const PathNode& right = paths[i];
      auto less = [&](const PathNode& a, const PathNode& b)
      {
        size_type order = std::min(a.order(), b.order());
        for(size_type j = 0; j < order; j++)
        {
          if(ranks[a.pointer() + j] != ranks[b.pointer() + j])
          {
            return ranks[a.pointer() + j] < ranks[b.pointer() + j];
          }
        }
        if(a.order() != b.order()) { return a.order() < b.order(); }
        if(a.from != b.from) { return a.from < b.from; }
        if(a.to != b.to) { return a.to < b.to; }
        if(a.predecessors() != b.predecessors()) { return a.predecessors() < b.predecessors(); }
        if(a.lcp() != b.lcp()) { return a.lcp() < b.lcp(); }
        for(size_type j = 0; j < a.ranks(); j++)
        {
          if(ranks[a.pointer() + j] != ranks[b.pointer() + j])
          {
            return ranks[a.pointer() + j] < ranks[b.pointer() + j];
          }
        }
        return false;
      };
      if(less(right, left)) { return false; }
    }
  }
  return true;
}

} // namespace

int main(int argc, char** argv)
{
  if(argc == 3 && std::string(argv[1]) == "gcsa-worker-task")
  {
    if(injectedWorkerFailure()) { return EXIT_FAILURE; }
    return externalPathJoinWorker(argv[2]);
  }
  omp_set_num_threads(1);

  // The disk ceiling applies to the peak physical framed pair, not merely its
  // logical PathNode/rank payload. Verify rejection in a child because the
  // legacy sorter error path deliberately terminates the current process.
  TempFileCodecParameters limit_codec(TempCompression::ZSTD,
    64 * KILOBYTE, 1, 1);
  const size_type limit_sort_budget = 8 * MEGABYTE;
  PathNode limit_node;
  limit_node.from = 1; limit_node.to = 2; limit_node.fields = 0;
  limit_node.setPredecessors(1); limit_node.setOrder(1);
  limit_node.setLCP(1); limit_node.setPointer(0);
  PathNode::rank_type limit_ranks[PathLabel::LABEL_LENGTH + 1] = { 17, 19 };
  const size_type logical_pair_bytes = sizeof(PathNode) +
    limit_node.ranks() * sizeof(PathNode::rank_type);
  const size_type peak_pair_bytes = externalPathGraphShardPeakBytes(1,
    limit_node.ranks(), limit_sort_budget, limit_codec);
  require(peak_pair_bytes > logical_pair_bytes);
  {
    PathGraph rejected_output(1, 0, 0);
    pid_t child = ::fork();
    require(child >= 0);
    if(child == 0)
    {
      size_type committed = 0;
      ExternalPathSortSink sink(rejected_output, 0, limit_sort_budget, 2,
        logical_pair_bytes, committed, nullptr, limit_codec);
      sink.write(limit_node, limit_ranks);
      std::_Exit(EXIT_SUCCESS);
    }
    int status = 0;
    require(::waitpid(child, &status, 0) == child);
    require(WIFEXITED(status) && WEXITSTATUS(status) == EXIT_FAILURE);
  }
  {
    PathGraph admitted_output(1, 0, 0);
    size_type committed = 0;
    ExternalPathSortSink sink(admitted_output, 0, limit_sort_budget, 2,
      peak_pair_bytes, committed, nullptr, limit_codec);
    sink.write(limit_node, limit_ranks);
    sink.finish();
    require(sink.storedBytes() == committed);
    require(committed > logical_pair_bytes && committed <= peak_pair_bytes);
  }

  std::vector<TestRecord> left;
  for(size_type i = 0; i < 40; i++)
  {
    left.push_back({ i + 1, 100, static_cast<byte_type>(1 << (i % 4)),
      static_cast<PathNode::rank_type>(10 + i), false });
  }
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
  // Pruning can collapse a unique range to order 0. It still has one boundary
  // rank, bypasses left expansion, and can be joined on its from node.
  right.push_back({ 100, std::numeric_limits<node_type>::max(), 8, 6000, true, 0, 0 });

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
  // extend() sorts each file after writing it; the recorded checksum must
  // describe the sorted bytes, not the unsorted ones written first.
  for(size_type file = 0; file < legacy.files(); file++)
  {
    requireRecordedChecksum(legacy.path_checksums[file], legacy.path_names[file]);
    requireRecordedChecksum(legacy.rank_checksums[file], legacy.rank_names[file]);
  }

  PathGraph external(left_path, left_rank);
  external.order = 1;
  external.logical_file_ids[0] = logical_file_id_t(7);
  external.physical_shard_ids[0] = physical_shard_id_t(101);
  appendShard(external, right_path, right_rank,
    logical_file_id_t(7), physical_shard_id_t(202));

  ConstructionParameters parameters;
  size_type memory_budget = externalPathJoinMinimumBudget() +
    externalPathGraphSortMinimumBudget();
  parameters.setMemoryLimitBytes(memory_budget);
  parameters.setJoinPartitionSize(externalPathJoinMinimumBudget());
  parameters.setSortRunSize(externalPathGraphSortMinimumBudget());
  parameters.setMergeFanIn(2);
  parameters.setWorkDirectory(".");
  parameters.setTempCompression("none");
  ExternalPathJoinStats stats;
  externalPathGraphExtend(external, GIGABYTE, parameters, &stats);

  requireSameGraph(legacy, external, "single-process external join");
  require(external.files() == 1);
  require(external.logicalFile(0) == logical_file_id_t(7));
  require(stats.generated_records == left.size() * right.size());
  require(stats.sorted_bypass == 2);
  require(stats.direct_label_records == stats.generated_records + stats.sorted_bypass);
  require(stats.intermediate_path_bytes_avoided > 0);
  require(stats.initial_runs > 2 && stats.merge_operations > 0);
  require(stats.label_sort_runs > 2 && stats.label_merge_passes > 0);
  require(stats.blocked_key_groups > 0);
  require(stats.compressed_join_runs == 0);
  require(stats.join_run_stored_bytes == stats.join_run_logical_bytes);
  require(stats.compressed_join_sidecars == 0);
  require(stats.join_sidecar_stored_bytes == stats.join_sidecar_logical_bytes);
  require(stats.max_bytes_resident <= memory_budget);

  // A tiny single-process budget forces the same high-fanout key through
  // several bounded disk-backed blocks (1201 records cannot fit in the
  // join-block reservation), while retaining the legacy result exactly.
  PathGraph blocked_graph(left_path, left_rank);
  blocked_graph.order = 1;
  blocked_graph.logical_file_ids[0] = logical_file_id_t(7);
  blocked_graph.physical_shard_ids[0] = physical_shard_id_t(101);
  appendShard(blocked_graph, right_path, right_rank,
    logical_file_id_t(7), physical_shard_id_t(202));
  ConstructionParameters blocked_parameters = parameters;
  ExternalPathJoinStats blocked_stats;
  externalPathGraphExtend(blocked_graph, GIGABYTE, blocked_parameters, &blocked_stats);
  require(readGraph(legacy) == readGraph(blocked_graph));
  require(blocked_stats.blocked_key_groups > 0);
  require(blocked_stats.blocked_key_blocks > 1);
  require(blocked_stats.max_bytes_resident <= blocked_parameters.getMemoryLimitBytes());

  // These are caps, not two simultaneous reservations. Distribution sorting
  // precedes label generation and may use the complete phase-local sort cap;
  // setting both expert caps to the global ceiling must therefore produce one
  // right run and one left run instead of starving distribution down to the
  // join minimum and creating a deep merge cascade.
  PathGraph roomy_graph(left_path, left_rank);
  roomy_graph.order = 1;
  roomy_graph.logical_file_ids[0] = logical_file_id_t(7);
  roomy_graph.physical_shard_ids[0] = physical_shard_id_t(101);
  appendShard(roomy_graph, right_path, right_rank,
    logical_file_id_t(7), physical_shard_id_t(202));
  ConstructionParameters roomy_parameters;
  roomy_parameters.setMemoryLimitBytes(8 * MEGABYTE);
  roomy_parameters.setSortRunSize(8 * MEGABYTE);
  roomy_parameters.setJoinPartitionSize(8 * MEGABYTE);
  roomy_parameters.setMergeFanIn(2);
  roomy_parameters.setWorkDirectory(".");
  roomy_parameters.setTempCompression("none");
  ExternalPathJoinStats roomy_stats;
  externalPathGraphExtend(roomy_graph, GIGABYTE, roomy_parameters,
    &roomy_stats);
  requireSameGraph(legacy, roomy_graph, "phase-local distribution budget");
  require(roomy_stats.initial_runs == 2);
  require(roomy_stats.distribution_sort_budget == 8 * MEGABYTE);
  require(roomy_stats.label_sort_budget + roomy_stats.join_block_budget <=
    roomy_parameters.getMemoryLimitBytes());
  require(roomy_stats.max_bytes_resident <=
    roomy_parameters.getMemoryLimitBytes());

  // The aggregate goal normally derives a 75/25 sort/join split. At 1 MiB,
  // however, the nominal join quarter is smaller than its fixed stream
  // buffers. The allocator must rebalance the still-feasible aggregate budget
  // instead of treating the preferred ratio as an infeasible hard cap.
  require(MEGABYTE >= externalPathJoinMinimumBudget() +
    externalPathGraphSortMinimumBudget());
  PathGraph automatic_graph(left_path, left_rank);
  automatic_graph.order = 1;
  automatic_graph.logical_file_ids[0] = logical_file_id_t(7);
  automatic_graph.physical_shard_ids[0] = physical_shard_id_t(101);
  appendShard(automatic_graph, right_path, right_rank,
    logical_file_id_t(7), physical_shard_id_t(202));
  ConstructionParameters automatic_parameters;
  automatic_parameters.setMemoryLimitBytes(MEGABYTE);
  automatic_parameters.setMergeFanIn(2);
  automatic_parameters.setWorkDirectory(".");
  // Make even the minimum framed writer workspace infeasible. AUTO must keep
  // the byte-identical raw JOIN stream instead of rejecting an otherwise
  // feasible aggregate memory goal.
  automatic_parameters.setCompressionWorkers(256);
  ExternalPathJoinStats automatic_stats;
  externalPathGraphExtend(automatic_graph, GIGABYTE,
    automatic_parameters, &automatic_stats);
  requireSameGraph(legacy, automatic_graph, "automatic aggregate budget");
  require(automatic_stats.compressed_join_runs == 0);
  require(automatic_stats.max_bytes_resident <=
    automatic_parameters.getMemoryLimitBytes());

  PathGraph rejected_graph(left_path, left_rank);
  rejected_graph.order = 1;
  rejected_graph.logical_file_ids[0] = logical_file_id_t(7);
  rejected_graph.physical_shard_ids[0] = physical_shard_id_t(101);
  appendShard(rejected_graph, right_path, right_rank,
    logical_file_id_t(7), physical_shard_id_t(202));
  ConstructionParameters rejected_parameters = automatic_parameters;
  rejected_parameters.setTempCompression("zstd");
  bool rejected_zstd = false;
  try
  {
    externalPathGraphExtend(rejected_graph, GIGABYTE, rejected_parameters);
  }
  catch(const std::runtime_error& error)
  {
    rejected_zstd = (std::string(error.what()).find(
      "join-sort budget cannot admit compressed join-run") != std::string::npos);
  }
  require(rejected_zstd);

  char workspace_root[] = "/tmp/gcsa-path-checkpoint-XXXXXX";
  require(mkdtemp(workspace_root) != nullptr);
  BuildWorkspace::Settings semantic;
  semantic["fixture"] = "external-join";
  BuildWorkspace workspace(workspace_root, semantic, BuildWorkspace::Settings(),
    BuildWorkspace::NEW_WORKSPACE);

  // The same skewed key must be splittable across fork-free worker processes.
  // Every physical result shard retains logical ID 7, and the PathGraph merger
  // must observe exactly the same semantic record multiset as the legacy path.
  PathGraph process_graph(left_path, left_rank);
  process_graph.order = 1;
  process_graph.logical_file_ids[0] = logical_file_id_t(7);
  process_graph.physical_shard_ids[0] = physical_shard_id_t(101);
  appendShard(process_graph, right_path, right_rank,
    logical_file_id_t(7), physical_shard_id_t(202));
  ConstructionParameters process_parameters = parameters;
  process_parameters.setMemoryLimitBytes(8 * MEGABYTE);
  // Four MiB admits the framed run writer plus two framed run readers and,
  // on the next doubling step, the compressed source-shard decoder pair.
  process_parameters.setJoinPartitionSize(4 * MEGABYTE);
  process_parameters.setSortRunSize(2 * MEGABYTE);
  process_parameters.setCheckpointBytes(32 * KILOBYTE);
  process_parameters.setProcessWorkers(2);
  process_parameters.setTempCompression("zstd");
  // 1 KiB frames force both the 40-summary stream and the recursively
  // split 1201-detail stream across several blocks. Planning consequently
  // exercises framed readAt() for nonzero subgroup offsets, rather than only
  // a sequential first-block read.
  process_parameters.setCompressionBlockSize(KILOBYTE);
  process_parameters.setCompressionWorkers(1);
  // The setter clamps to MIN_OPEN_FILES (64). Compaction retains a bounded
  // multi-run frontier rather than pretending this request limits it to 8.
  process_parameters.setMaxOpenFiles(8);
  process_parameters.setWorkerExecutable(argv[0]);
  omp_set_num_threads(4);
  ExternalPathJoinStats process_stats;
  externalPathGraphExtend(process_graph, GIGABYTE,
    process_parameters, &process_stats, &workspace, "step-process");
  omp_set_num_threads(1);
  require(readGraph(legacy) == readGraph(process_graph));
  // Many worker partitions remain a bounded set of individually LABEL-sorted
  // streams for the downstream merger.
  require(process_graph.files() > 1);
  require(process_graph.files() <= (process_parameters.getMaxOpenFiles() - 4) / 2);
  require(labelSorted(process_graph));
  for(size_type file = 0; file < process_graph.files(); file++)
  {
    require(process_graph.logicalFile(file) == logical_file_id_t(7));
    require(process_graph.physicalShard(file) == physical_shard_id_t(file));
    require(CompressedBlockReader::isFramed(process_graph.path_names[file]));
    require(CompressedBlockReader::isFramed(process_graph.rank_names[file]));
  }
  require(process_stats.join_partitions > 1);
  // The sampled planner must turn this single high-fanout semantic key into
  // several deterministic range tasks rather than one worker monopoly.
  require(process_stats.join_partitions >= 4);
  require(process_stats.worker_processes == process_stats.join_partitions);
  require(process_stats.recursive_splits > 0);
  require(process_stats.left_range_splits > 0);
  require(process_stats.right_range_splits > 0);
  require(process_stats.generated_records == left.size() * right.size());
  require(process_stats.sorted_bypass == 2);
  require(process_stats.sampled_plan_records > 0);
  require(process_stats.radix_plan_bins >= 1);
  require(process_stats.sidecar_plan_groups > 0);
  require(process_stats.sidecar_plan_detail_records > 0);
  require(process_stats.sidecar_plan_groups * 48 > KILOBYTE);
  require(process_stats.sidecar_plan_detail_records > KILOBYTE);
  require(process_stats.full_record_plan_rescans == 0);
  require(process_stats.compressed_join_runs > 0);
  require(process_stats.join_run_stored_bytes <
    process_stats.join_run_logical_bytes);
  require(process_stats.compressed_join_sidecars > 0);
  require(process_stats.join_sidecar_stored_bytes <
    process_stats.join_sidecar_logical_bytes);
  require(process_stats.grouped_expansion_records > 0);
  require(process_stats.expansion_context_bytes_saved > 0);
  require(process_stats.max_bytes_resident <= process_parameters.getMemoryLimitBytes());

  // Rebuild the cross-file, high-fanout worker plan without a workspace. The
  // compact sidecars must produce the same physical streams, not merely the
  // same logical record multiset.
  PathGraph fresh_process(left_path, left_rank);
  fresh_process.order = 1;
  fresh_process.logical_file_ids[0] = logical_file_id_t(7);
  fresh_process.physical_shard_ids[0] = physical_shard_id_t(101);
  appendShard(fresh_process, right_path, right_rank,
    logical_file_id_t(7), physical_shard_id_t(202));
  ConstructionParameters fresh_parameters = process_parameters;
  fresh_parameters.setVerifyWorkspace();
  ExternalPathJoinStats fresh_stats;
  externalPathGraphExtend(fresh_process, GIGABYTE,
    fresh_parameters, &fresh_stats);
  requireSameGraph(process_graph, fresh_process, "fresh sidecar worker plan");
  requireByteIdenticalGraph(process_graph, fresh_process);
  require(fresh_stats.full_record_plan_rescans == 0);

  PathGraph resumed_process(left_path, left_rank);
  resumed_process.order = 1;
  resumed_process.logical_file_ids[0] = logical_file_id_t(7);
  resumed_process.physical_shard_ids[0] = physical_shard_id_t(101);
  appendShard(resumed_process, right_path, right_rank,
    logical_file_id_t(7), physical_shard_id_t(202));
  ExternalPathJoinStats resumed_stats;
  externalPathGraphExtend(resumed_process, GIGABYTE,
    process_parameters, &resumed_stats, &workspace, "step-process");
  require(readGraph(process_graph) == readGraph(resumed_process));
  requireByteIdenticalGraph(process_graph, resumed_process);
  require(resumed_process.files() == process_graph.files());
  require(labelSorted(resumed_process));
  for(size_type file = 0; file < resumed_process.files(); file++)
  {
    require(resumed_process.physicalShard(file) == physical_shard_id_t(file));
  }
  require(resumed_stats.worker_processes == 0);
  require(resumed_stats.restored_partitions == resumed_stats.join_partitions);
  require(resumed_stats.restored_radix_plans == 1);
  require(resumed_stats.restored_sidecar_metadata == 2);
  require(resumed_stats.full_record_plan_rescans == 0);

  // Consume the restored compressed shards in another doubling step. This is
  // the important mixed-format boundary: the distribution sorter and two
  // framed decode workspaces must share the same join budget rather than
  // adding an untracked block cache on top of it.
  PathGraph second_legacy(combined_path, combined_rank);
  second_legacy.order = 1;
  second_legacy.extend(GIGABYTE, 64 * MEGABYTE);
  second_legacy.extend(GIGABYTE, 64 * MEGABYTE);
  ConstructionParameters second_parameters = process_parameters;
  second_parameters.setProcessWorkers(1);
  ExternalPathJoinStats second_stats;
  externalPathGraphExtend(resumed_process, GIGABYTE, second_parameters,
    &second_stats);
  requireSameGraph(second_legacy, resumed_process,
    "compressed source doubling step");
  require(second_stats.max_bytes_resident <=
    second_parameters.getMemoryLimitBytes());

  // Spread exact join keys across the most-significant nibble of node_type.
  // The 4 KiB target cannot admit their combined output, but each key fits by
  // itself. This forces a sampled MSD split and proves that the exact scan
  // flushes at the persisted range-pack boundaries without changing results.
  std::vector<TestRecord> radix_left, radix_right, radix_combined;
  for(size_type digit = 1; digit <= 8; digit++)
  {
    node_type key = (static_cast<node_type>(digit) << 60) | 0x1234;
    for(size_type i = 0; i < 8; i++)
    {
      radix_left.push_back({ 10000 + 100 * digit + i, key,
        static_cast<byte_type>(1U << (i % 4)),
        static_cast<PathNode::rank_type>(100 * digit + i), false });
      radix_right.push_back({ key, 20000 + 100 * digit + i,
        static_cast<byte_type>(1),
        static_cast<PathNode::rank_type>(1000 + 100 * digit + i), false });
    }
  }
  radix_combined = radix_left;
  radix_combined.insert(radix_combined.end(), radix_right.begin(), radix_right.end());
  std::string radix_combined_path = base + ".radix-combined.path";
  std::string radix_combined_rank = base + ".radix-combined.rank";
  std::string radix_left_path = base + ".radix-left.path";
  std::string radix_left_rank = base + ".radix-left.rank";
  std::string radix_right_path = base + ".radix-right.path";
  std::string radix_right_rank = base + ".radix-right.rank";
  writePathPair(radix_combined_path, radix_combined_rank, radix_combined);
  writePathPair(radix_left_path, radix_left_rank, radix_left);
  writePathPair(radix_right_path, radix_right_rank, radix_right);

  PathGraph radix_legacy(radix_combined_path, radix_combined_rank);
  radix_legacy.order = 1;
  radix_legacy.extend(GIGABYTE, 64 * MEGABYTE);
  PathGraph radix_external(radix_left_path, radix_left_rank);
  radix_external.order = 1;
  radix_external.logical_file_ids[0] = logical_file_id_t(11);
  radix_external.physical_shard_ids[0] = physical_shard_id_t(301);
  appendShard(radix_external, radix_right_path, radix_right_rank,
    logical_file_id_t(11), physical_shard_id_t(302));
  ConstructionParameters radix_parameters = process_parameters;
  radix_parameters.setCheckpointBytes(4 * KILOBYTE);
  ExternalPathJoinStats radix_stats;
  externalPathGraphExtend(radix_external, GIGABYTE, radix_parameters,
    &radix_stats, &workspace, "step-msd");
  require(readGraph(radix_legacy) == readGraph(radix_external));
  require(radix_stats.radix_plan_splits > 0);
  require(radix_stats.radix_plan_bins > 1);
  require(radix_stats.radix_plan_max_bits >= 4);
  require(radix_stats.radix_boundary_flushes > 0);
  require(radix_stats.full_record_plan_rescans == 0);

  PathGraph radix_resumed(radix_left_path, radix_left_rank);
  radix_resumed.order = 1;
  radix_resumed.logical_file_ids[0] = logical_file_id_t(11);
  radix_resumed.physical_shard_ids[0] = physical_shard_id_t(301);
  appendShard(radix_resumed, radix_right_path, radix_right_rank,
    logical_file_id_t(11), physical_shard_id_t(302));
  ExternalPathJoinStats radix_resumed_stats;
  externalPathGraphExtend(radix_resumed, GIGABYTE, radix_parameters,
    &radix_resumed_stats, &workspace, "step-msd");
  require(readGraph(radix_external) == readGraph(radix_resumed));
  requireByteIdenticalGraph(radix_external, radix_resumed);
  require(radix_resumed_stats.restored_radix_plans == 1);
  require(radix_resumed_stats.restored_sidecar_metadata == 2);
  require(radix_resumed_stats.worker_processes == 0);

  // A same-length run with different payload must not inherit either the plan
  // or completed worker ranges. Run checksums are part of both task identities.
  std::vector<TestRecord> changed_left = radix_left;
  changed_left.front().label += 50000;
  std::vector<TestRecord> changed_combined = changed_left;
  changed_combined.insert(changed_combined.end(), radix_right.begin(), radix_right.end());
  std::string changed_combined_path = base + ".changed-combined.path";
  std::string changed_combined_rank = base + ".changed-combined.rank";
  std::string changed_left_path = base + ".changed-left.path";
  std::string changed_left_rank = base + ".changed-left.rank";
  writePathPair(changed_combined_path, changed_combined_rank, changed_combined);
  writePathPair(changed_left_path, changed_left_rank, changed_left);
  PathGraph changed_legacy(changed_combined_path, changed_combined_rank);
  changed_legacy.order = 1;
  changed_legacy.extend(GIGABYTE, 64 * MEGABYTE);
  PathGraph changed_external(changed_left_path, changed_left_rank);
  changed_external.order = 1;
  changed_external.logical_file_ids[0] = logical_file_id_t(11);
  changed_external.physical_shard_ids[0] = physical_shard_id_t(301);
  appendShard(changed_external, radix_right_path, radix_right_rank,
    logical_file_id_t(11), physical_shard_id_t(302));
  ExternalPathJoinStats changed_stats;
  externalPathGraphExtend(changed_external, GIGABYTE, radix_parameters,
    &changed_stats, &workspace, "step-msd");
  require(readGraph(changed_legacy) == readGraph(changed_external));
  require(changed_stats.restored_radix_plans == 0);
  require(changed_stats.restored_partitions == 0);
  require(changed_stats.worker_processes == changed_stats.join_partitions);

  // Requesting multiple workers must remain valid when the deterministic plan
  // contains only one task. This is common for small logical inputs.
  PathGraph one_partition(left_path, left_rank);
  one_partition.order = 1;
  one_partition.logical_file_ids[0] = logical_file_id_t(7);
  one_partition.physical_shard_ids[0] = physical_shard_id_t(101);
  appendShard(one_partition, right_path, right_rank,
    logical_file_id_t(7), physical_shard_id_t(202));
  ConstructionParameters one_parameters = process_parameters;
  one_parameters.setCheckpointBytes(GIGABYTE);
  ExternalPathJoinStats one_stats;
  externalPathGraphExtend(one_partition, GIGABYTE, one_parameters, &one_stats);
  require(readGraph(legacy) == readGraph(one_partition));
  require(one_partition.files() == 1);
  require(one_stats.join_partitions == 1);
  require(one_stats.worker_processes == 1);

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
  require(separated.size() == 2); // Only the already-sorted bypass paths remain.

  // Join workers leave many physical shards per logical input, and compaction
  // merges them in batches. With threads to spare, the batches of a pass run
  // at once: the published frontier must be the one the serial loop publishes,
  // stream for stream, and the step must report that the batches overlapped.
  // Stored bytes are not compared, because a concurrent writer compresses in
  // its own thread rather than with the configured zstd workers.
  {
    ConstructionParameters compaction_parameters;
    compaction_parameters.setMemoryLimitBytes(2 * GIGABYTE);
    compaction_parameters.setTempCompression("zstd");
    compaction_parameters.setCompressionBlockSize(KILOBYTE);
    compaction_parameters.setCompressionWorkers(2);
    const size_type logical_inputs = 2, records = 50;
    const size_type target = pathMergeInputPairs(compaction_parameters.getMaxOpenFiles(),
      logical_inputs) / logical_inputs;
    // Several multi-shard batches per logical input, and a trailing singleton.
    const size_type shards = 4 * target + 1;
    const size_type label_fan_in = std::min(compaction_parameters.getMergeFanIn(),
      compaction_parameters.getMaxOpenFiles() - 1);
    auto build_shards = [&](PathGraph& graph, const std::string& tag)
    {
      for(size_type logical = 0; logical < logical_inputs; logical++)
      {
        for(size_type shard = 0; shard < shards; shard++)
        {
          const std::string name = base + ".compaction." + tag + "." +
            std::to_string(logical) + "." + std::to_string(shard);
          std::vector<TestRecord> shard_records;
          for(size_type i = 0; i < records; i++)
          {
            shard_records.push_back({ i + 1, i + 2, static_cast<byte_type>(1),
              static_cast<PathNode::rank_type>((i * shards + shard) * logical_inputs + logical + 1),
              false });
          }
          writePathPair(name + ".path", name + ".rank", shard_records);
          appendShard(graph, name + ".path", name + ".rank",
            logical_file_id_t(logical), physical_shard_id_t(logical * shards + shard));
        }
      }
    };
    PathGraph serial_graph(0, 2, 1), concurrent_graph(0, 2, 1);
    build_shards(serial_graph, "serial");
    build_shards(concurrent_graph, "concurrent");
    ExternalPathJoinStats serial_stats, concurrent_stats;
    size_type serial_committed = 0, concurrent_committed = 0;
    MemoryBudget serial_memory(compaction_parameters.getMemoryLimitBytes(),
      compaction_parameters.getMemoryLimitBytes() / 16);
    MemoryBudget concurrent_memory(compaction_parameters.getMemoryLimitBytes(),
      compaction_parameters.getMemoryLimitBytes() / 16);
    omp_set_num_threads(1);
    compactLogicalJoinShards(serial_graph, GIGABYTE, compaction_parameters,
      label_fan_in, serial_memory, serial_committed, &serial_stats);
    omp_set_num_threads(4);
    compactLogicalJoinShards(concurrent_graph, GIGABYTE, compaction_parameters,
      label_fan_in, concurrent_memory, concurrent_committed, &concurrent_stats);
    omp_set_num_threads(1);
    require(serial_stats.compaction_batches > logical_inputs);
    require(serial_stats.compaction_concurrency == 1);
    require(concurrent_stats.compaction_batches == serial_stats.compaction_batches);
    require(concurrent_stats.compaction_concurrency > 1);
    require(concurrent_stats.max_bytes_resident <= compaction_parameters.getMemoryLimitBytes());
    require(serial_graph.files() == concurrent_graph.files());
    require(serial_graph.files() <= target * logical_inputs);
    require(serial_graph.path_count == concurrent_graph.path_count);
    require(serial_graph.rank_count == concurrent_graph.rank_count);
    require(serial_graph.path_count == logical_inputs * shards * records);
    require(labelSorted(concurrent_graph));
    for(size_type file = 0; file < serial_graph.files(); file++)
    {
      require(serial_graph.logicalFile(file) == concurrent_graph.logicalFile(file));
      require(serial_graph.physicalShard(file) == concurrent_graph.physicalShard(file));
      require(serial_graph.path_counts[file] == concurrent_graph.path_counts[file]);
      require(serial_graph.rank_counts[file] == concurrent_graph.rank_counts[file]);
      std::vector<PathNode> serial_paths, concurrent_paths;
      std::vector<PathNode::rank_type> serial_ranks, concurrent_ranks;
      serial_graph.read(serial_paths, serial_ranks, file);
      concurrent_graph.read(concurrent_paths, concurrent_ranks, file);
      require(serial_paths.size() == concurrent_paths.size());
      require(std::memcmp(serial_paths.data(), concurrent_paths.data(),
        serial_paths.size() * sizeof(PathNode)) == 0);
      require(serial_ranks == concurrent_ranks);
    }
  }

  // Concurrent batches share the descriptor budget instead of each claiming
  // all of it. A transient user service runs with a soft limit of 1,024
  // descriptors, which unbounded 64-thread compaction of the joint chr2+chr18
  // build exhausted; here the budget equals the ceiling, and the process limit
  // sits exactly that far above the descriptors already open.
  {
    ConstructionParameters limited_parameters;
    limited_parameters.setMemoryLimitBytes(2 * GIGABYTE);
    limited_parameters.setTempCompression("zstd");
    limited_parameters.setCompressionBlockSize(KILOBYTE);
    limited_parameters.setMaxOpenFiles(ConstructionParameters::MIN_OPEN_FILES);
    const size_type max_open = limited_parameters.getMaxOpenFiles();
    limited_parameters.setConcurrentOpenFiles(max_open);
    const size_type logical_inputs = 2, records = 50;
    const size_type target = pathMergeInputPairs(max_open, logical_inputs) / logical_inputs;
    const size_type shards = 4 * target + 1, width = (shards + target - 1) / target;
    auto build_limited = [&](PathGraph& target_graph, const std::string& tag)
    {
      for(size_type logical = 0; logical < logical_inputs; logical++)
      {
        for(size_type shard = 0; shard < shards; shard++)
        {
          const std::string name = base + ".limited." + tag + "." + std::to_string(logical) + "." + std::to_string(shard);
          std::vector<TestRecord> shard_records;
          for(size_type i = 0; i < records; i++)
          {
            shard_records.push_back({ i + 1, i + 2, static_cast<byte_type>(1),
              static_cast<PathNode::rank_type>((i * shards + shard) * logical_inputs + logical + 1), false });
          }
          writePathPair(name + ".path", name + ".rank", shard_records);
          appendShard(target_graph, name + ".path", name + ".rank",
            logical_file_id_t(logical), physical_shard_id_t(logical * shards + shard));
        }
      }
    };
    PathGraph graph(0, 2, 1);
    build_limited(graph, "ceiling");
    ExternalPathJoinStats stats;
    size_type committed = 0;
    MemoryBudget memory(limited_parameters.getMemoryLimitBytes(),
      limited_parameters.getMemoryLimitBytes() / 16);
    size_type open_now = 0; // Includes the listing's own descriptor, closed again below.
    for(const auto& entry : std::filesystem::directory_iterator("/proc/self/fd")) { static_cast<void>(entry); open_now++; }
    struct rlimit original, limited;
    require(::getrlimit(RLIMIT_NOFILE, &original) == 0);
    limited = original; limited.rlim_cur = open_now - 1 + max_open;
    require(::setrlimit(RLIMIT_NOFILE, &limited) == 0);
    omp_set_num_threads(8);
    bool completed = true;
    try
    {
      compactLogicalJoinShards(graph, GIGABYTE, limited_parameters,
        std::min(limited_parameters.getMergeFanIn(), max_open - 1), memory, committed, &stats);
    }
    catch(const std::exception&) { completed = false; }
    omp_set_num_threads(1);
    require(::setrlimit(RLIMIT_NOFILE, &original) == 0);
    require(completed);
    require(stats.compaction_concurrency > 1);
    require(stats.compaction_concurrency * (2 * width + 4) <= max_open);
    require(graph.files() <= target * logical_inputs);
    require(graph.path_count == logical_inputs * shards * records);
    require(labelSorted(graph));

    // A larger concurrent budget runs more of the same batches at once and
    // leaves their width, and so the published frontier, unchanged.
    limited_parameters.setConcurrentOpenFiles(1024);
    PathGraph budgeted(0, 2, 1);
    build_limited(budgeted, "budgeted");
    ExternalPathJoinStats budgeted_stats;
    size_type budgeted_committed = 0;
    MemoryBudget budgeted_memory(limited_parameters.getMemoryLimitBytes(),
      limited_parameters.getMemoryLimitBytes() / 16);
    omp_set_num_threads(8);
    compactLogicalJoinShards(budgeted, GIGABYTE, limited_parameters,
      std::min(limited_parameters.getMergeFanIn(), max_open - 1), budgeted_memory,
      budgeted_committed, &budgeted_stats);
    omp_set_num_threads(1);
    require(budgeted_stats.compaction_batches == stats.compaction_batches);
    require(budgeted_stats.compaction_concurrency > stats.compaction_concurrency);
    require(budgeted.files() == graph.files());
    require(budgeted.path_count == graph.path_count);
    for(size_type file = 0; file < graph.files(); file++)
    {
      require(budgeted.path_counts[file] == graph.path_counts[file]);
      require(budgeted.rank_counts[file] == graph.rank_counts[file]);
    }
  }

  // Key-range distribution. A logical input of at least 2 x 65,536 paths is
  // split into key ranges, one per thread, whose distributions run at once.
  // Every route must produce the serial route's records; the ranged process
  // worker route must also restore its range plan and every partition when
  // the same step runs again on its workspace.
  {
    const size_type chain = 4 * 65536;
    std::vector<TestRecord> chain_records;
    chain_records.reserve(chain);
    for(size_type i = 0; i < chain; i++)
    {
      // Path i runs from node i + 1 to node i + 2, so it extends path i + 1.
      chain_records.push_back({ i + 1, i + 2, static_cast<byte_type>(1 << (i % 4)),
        static_cast<PathNode::rank_type>(i), false });
    }
    const std::string chain_path = base + ".chain.path";
    const std::string chain_rank = base + ".chain.rank";
    writePathPair(chain_path, chain_rank, chain_records);
    auto chain_graph = [&](PathGraph& graph)
    {
      graph.order = 1;
      graph.logical_file_ids[0] = logical_file_id_t(3);
      graph.physical_shard_ids[0] = physical_shard_id_t(30);
    };
    ConstructionParameters range_parameters;
    range_parameters.setMemoryLimitBytes(64 * MEGABYTE);
    range_parameters.setCheckpointBytes(512 * KILOBYTE);
    range_parameters.setWorkDirectory(".");
    range_parameters.setTempCompression("none");

    PathGraph serial_chain(chain_path, chain_rank);
    chain_graph(serial_chain);
    ExternalPathJoinStats serial_chain_stats;
    omp_set_num_threads(1);
    externalPathGraphExtend(serial_chain, GIGABYTE, range_parameters, &serial_chain_stats);
    require(serial_chain_stats.distribution_ranges == 1);
    require(serial_chain_stats.distribution_concurrency == 1);
    require(serial_chain_stats.generated_records == chain - 1);
    require(serial_chain.size() == chain - 1);

    // Direct joins, one output shard per range.
    PathGraph direct_chain(chain_path, chain_rank);
    chain_graph(direct_chain);
    ExternalPathJoinStats direct_chain_stats;
    omp_set_num_threads(4);
    externalPathGraphExtend(direct_chain, GIGABYTE, range_parameters, &direct_chain_stats);
    omp_set_num_threads(1);
    require(direct_chain_stats.distribution_ranges == 4);
    // The single shard is read once, in chunks, by one reader per thread.
    require(direct_chain_stats.distribution_concurrency == 4);
    require(direct_chain_stats.left_records == serial_chain_stats.left_records);
    require(direct_chain_stats.right_records == serial_chain_stats.right_records);
    require(direct_chain_stats.max_bytes_resident <= range_parameters.getMemoryLimitBytes());
    require(direct_chain.files() == 4);
    requireSameGraph(serial_chain, direct_chain, "key-range distribution, direct joins");

    // Process workers over every range's partitions, on a workspace.
    char range_root[] = "/tmp/gcsa-range-plan-XXXXXX";
    require(mkdtemp(range_root) != nullptr);
    BuildWorkspace::Settings range_semantic;
    range_semantic["fixture"] = "external-join-ranges";
    BuildWorkspace range_workspace(range_root, range_semantic, BuildWorkspace::Settings(),
      BuildWorkspace::NEW_WORKSPACE);
    ConstructionParameters worker_range_parameters = range_parameters;
    worker_range_parameters.setProcessWorkers(2);
    worker_range_parameters.setWorkerExecutable(argv[0]);
    PathGraph worker_chain(chain_path, chain_rank);
    chain_graph(worker_chain);
    ExternalPathJoinStats worker_chain_stats;
    omp_set_num_threads(4);
    externalPathGraphExtend(worker_chain, GIGABYTE, worker_range_parameters,
      &worker_chain_stats, &range_workspace, "step-ranges");
    omp_set_num_threads(1);
    require(worker_chain_stats.distribution_ranges == 4);
    require(worker_chain_stats.distribution_concurrency == 4); // one reader per thread, in chunks
    require(worker_chain_stats.restored_range_plans == 0);
    require(worker_chain_stats.join_partitions >= 4);
    require(worker_chain_stats.worker_processes == worker_chain_stats.join_partitions);
    requireSameGraph(serial_chain, worker_chain, "key-range distribution, process workers");

    // The same step again: the committed range plan gives the same ranges
    // whatever the thread count, so the runs, their checksums and every
    // partition task name match, and nothing is joined again.
    PathGraph resumed_chain(chain_path, chain_rank);
    chain_graph(resumed_chain);
    ExternalPathJoinStats resumed_chain_stats;
    omp_set_num_threads(2);
    externalPathGraphExtend(resumed_chain, GIGABYTE, worker_range_parameters,
      &resumed_chain_stats, &range_workspace, "step-ranges");
    omp_set_num_threads(1);
    require(resumed_chain_stats.restored_range_plans == 1);
    require(resumed_chain_stats.distribution_ranges == 4);
    require(resumed_chain_stats.distribution_concurrency == 2); // one reader per thread, in chunks
    require(resumed_chain_stats.restored_radix_plans > 0);
    require(resumed_chain_stats.restored_partitions == worker_chain_stats.join_partitions);
    require(resumed_chain_stats.worker_processes == 0);
    requireByteIdenticalGraph(worker_chain, resumed_chain);

    std::remove(chain_path.c_str()); std::remove(chain_rank.c_str());
    std::filesystem::remove_all(range_root);
  }

  // One worker pool across logical inputs. Inputs 7 and 9 each plan several
  // partitions. Joined in one step, their partitions share a single pool,
  // launched largest first across both inputs, yet each input's output shards,
  // their numbering, and their partition task names and artifacts must be
  // exactly those of the input joined by itself (which is what a per-input pool
  // produced), and a step whose pool fails midway must resume to the same bytes.
  {
    std::vector<TestRecord> second_left, second_right;
    for(size_type i = 0; i < 30; i++)
    {
      second_left.push_back({ 5000 + i, 200, static_cast<byte_type>(1 << (i % 4)),
        static_cast<PathNode::rank_type>(20000 + i), false });
    }
    for(size_type i = 0; i < 900; i++)
    {
      second_right.push_back({ 200, 6000 + i, static_cast<byte_type>(1),
        static_cast<PathNode::rank_type>(21000 + i), false });
    }
    const std::string second_left_path = base + ".pool-left.path";
    const std::string second_left_rank = base + ".pool-left.rank";
    const std::string second_right_path = base + ".pool-right.path";
    const std::string second_right_rank = base + ".pool-right.rank";
    writePathPair(second_left_path, second_left_rank, second_left);
    writePathPair(second_right_path, second_right_rank, second_right);
    auto first_input = [&](PathGraph& graph)
    {
      graph.logical_file_ids[0] = logical_file_id_t(7);
      graph.physical_shard_ids[0] = physical_shard_id_t(101);
      appendShard(graph, right_path, right_rank, logical_file_id_t(7), physical_shard_id_t(202));
    };
    auto second_shards = [&](PathGraph& graph)
    {
      appendShard(graph, second_right_path, second_right_rank,
        logical_file_id_t(9), physical_shard_id_t(302));
    };
    BuildWorkspace::Settings pool_semantic;
    pool_semantic["fixture"] = "external-join-pool";
    // Few enough partitions per input that the step does not compact them, so
    // its output graph is the workers' shards themselves.
    ConstructionParameters pool_parameters = process_parameters;
    pool_parameters.setCheckpointBytes(320 * KILOBYTE);
    auto pool_step = [&](PathGraph& graph, const std::string& root,
      BuildWorkspace::OpenMode mode, ExternalPathJoinStats& pool_stats)
    {
      BuildWorkspace pool_workspace(root, pool_semantic, BuildWorkspace::Settings(), mode);
      graph.order = 1;
      omp_set_num_threads(4);
      try
      {
        externalPathGraphExtend(graph, GIGABYTE, pool_parameters, &pool_stats,
          &pool_workspace, "step-pool");
      }
      catch(...) { omp_set_num_threads(1); throw; }
      omp_set_num_threads(1);
    };
    auto make_root = [](const char* stem)
    {
      std::string pattern = std::string("/tmp/") + stem + "-XXXXXX";
      std::vector<char> name(pattern.begin(), pattern.end()); name.push_back(0);
      require(mkdtemp(name.data()) != nullptr);
      return std::string(name.data());
    };
    const std::string alone_first_root = make_root("gcsa-pool-first");
    const std::string alone_second_root = make_root("gcsa-pool-second");
    const std::string pooled_root = make_root("gcsa-pool-both");
    const std::string interrupted_root = make_root("gcsa-pool-interrupted");
    const std::string claims_root = make_root("gcsa-pool-claims");
    // The interrupted step below kills a running worker, whose partially
    // written output shard (<name>.partial and its index spool) the
    // coordinator does not name; keep the step's temporaries in a directory
    // this block removes.
    const std::string pool_temp_root = make_root("gcsa-pool-temp");
    TempFile::setDirectory(pool_temp_root);

    PathGraph alone_first(left_path, left_rank);
    first_input(alone_first);
    ExternalPathJoinStats alone_first_stats;
    pool_step(alone_first, alone_first_root, BuildWorkspace::NEW_WORKSPACE, alone_first_stats);
    PathGraph alone_second(second_left_path, second_left_rank);
    alone_second.logical_file_ids[0] = logical_file_id_t(9);
    alone_second.physical_shard_ids[0] = physical_shard_id_t(301);
    second_shards(alone_second);
    ExternalPathJoinStats alone_second_stats;
    pool_step(alone_second, alone_second_root, BuildWorkspace::NEW_WORKSPACE, alone_second_stats);
    const size_type first_partitions = alone_first_stats.join_partitions;
    const size_type second_partitions = alone_second_stats.join_partitions;
    // Two workers run in every step below, so each worker's reservation, and
    // with it the framed bytes it writes, is the same alone and pooled. Three
    // partitions per input keep one committed for each input when the pool
    // fails on its last launch below.
    require(first_partitions >= 3 && second_partitions >= 3);
    require(alone_first_stats.join_pool_inputs == 1);
    require(alone_first_stats.join_pool_concurrent_inputs == 0);
    // Too few shards to compact: the step's output is the workers' shards.
    require(alone_first.files() == first_partitions);
    require(alone_second.files() == second_partitions);

    auto both_inputs = [&](PathGraph& graph)
    {
      first_input(graph);
      appendShard(graph, second_left_path, second_left_rank,
        logical_file_id_t(9), physical_shard_id_t(301));
      second_shards(graph);
    };
    PathGraph pooled(left_path, left_rank);
    both_inputs(pooled);
    ExternalPathJoinStats pooled_stats;
    pool_step(pooled, pooled_root, BuildWorkspace::NEW_WORKSPACE, pooled_stats);
    std::cerr << "test_external_join: one pool for " << pooled_stats.join_partitions
              << " partitions of " << pooled_stats.join_pool_inputs
              << " logical inputs (" << first_partitions << " + " << second_partitions
              << "), at most " << pooled_stats.join_pool_concurrent_inputs
              << " inputs with a worker running at once" << std::endl;
    require(pooled_stats.join_pool_inputs == 2);
    require(pooled_stats.join_pool_concurrent_inputs == 2);
    require(pooled_stats.join_partitions == first_partitions + second_partitions);
    require(pooled_stats.worker_processes == pooled_stats.join_partitions);
    require(pooled_stats.generated_records ==
      alone_first_stats.generated_records + alone_second_stats.generated_records);
    // Input 7's shards come first, numbered 0..n-1, then input 9's, each
    // byte-identical to the same input's shard when joined alone.
    require(pooled.files() == first_partitions + second_partitions);
    for(size_type file = 0; file < pooled.files(); file++)
    {
      const bool first = (file < first_partitions);
      const PathGraph& alone = (first ? alone_first : alone_second);
      const size_type alone_file = (first ? file : file - first_partitions);
      require(pooled.physicalShard(file) == physical_shard_id_t(file));
      require(pooled.logicalFile(file) == alone.logicalFile(alone_file));
      require(pooled.path_counts[file] == alone.path_counts[alone_file]);
      require(pooled.rank_counts[file] == alone.rank_counts[alone_file]);
      require(readBytes(pooled.path_names[file]) == readBytes(alone.path_names[alone_file]));
      require(readBytes(pooled.rank_names[file]) == readBytes(alone.rank_names[alone_file]));
    }
    // The same partition tasks, plans and artifacts as the two inputs alone.
    std::map<std::string, std::vector<std::uint8_t>> expected_tasks =
      workspaceEntries(alone_first_root, "step-pool-");
    const std::map<std::string, std::vector<std::uint8_t>> second_tasks =
      workspaceEntries(alone_second_root, "step-pool-");
    expected_tasks.insert(second_tasks.begin(), second_tasks.end());
    require(expected_tasks.size() > first_partitions + second_partitions);
    require(workspaceEntries(pooled_root, "step-pool-") == expected_tasks);

    // Disk admission. Each input's plan must fit the limit by itself, as it had
    // to in its own pool, and the pool needs no more: at the larger input's
    // planned peaks, eight workers' running peaks exceed the limit, so launches
    // wait for running workers, whose stored pairs are far below their peaks.
    // One byte less, and that input's own plan is refused before any launch.
    {
      ConstructionParameters wide_parameters = pool_parameters;
      wide_parameters.setMemoryLimitBytes(32 * MEGABYTE);
      wide_parameters.setProcessWorkers(8);
      auto wide_step = [&](PathGraph& graph, size_type limit, ExternalPathJoinStats& wide_stats)
      {
        both_inputs(graph);
        graph.order = 1;
        omp_set_num_threads(4);
        try { externalPathGraphExtend(graph, limit, wide_parameters, &wide_stats); }
        catch(...) { omp_set_num_threads(1); throw; }
        omp_set_num_threads(1);
      };
      PathGraph wide(left_path, left_rank);
      ExternalPathJoinStats wide_stats;
      wide_step(wide, GIGABYTE, wide_stats);
      const size_type input_limit = wide_stats.join_largest_input_plan_bytes;
      require(wide_stats.join_pool_disk_waits == 0);
      PathGraph limited(left_path, left_rank);
      ExternalPathJoinStats limited_stats;
      wide_step(limited, input_limit, limited_stats);
      std::cerr << "test_external_join: at the larger input's planned "
                << input_limit << " bytes, the pool waited for "
                << limited_stats.join_pool_disk_waits << " workers to free disk" << std::endl;
      require(limited_stats.join_pool_disk_waits > 0);
      require(limited.storedBytes() <= input_limit);
      requireByteIdenticalGraph(wide, limited);
      PathGraph refused(left_path, left_rank);
      ExternalPathJoinStats refused_stats;
      bool refused_plan = false;
      try { wide_step(refused, input_limit - 1, refused_stats); }
      catch(const std::runtime_error& error)
      {
        refused_plan = (std::string(error.what()).find(
          "exceeded by deterministic partition plan: logical input ") != std::string::npos);
      }
      require(refused_plan);
      require(refused_stats.worker_processes == 0);
    }

    // The same step again on its workspace restores every partition.
    PathGraph restored_pool(left_path, left_rank);
    both_inputs(restored_pool);
    ExternalPathJoinStats restored_pool_stats;
    pool_step(restored_pool, pooled_root, BuildWorkspace::RESUME, restored_pool_stats);
    require(restored_pool_stats.restored_partitions == pooled_stats.join_partitions);
    require(restored_pool_stats.worker_processes == 0);
    requireByteIdenticalGraph(pooled, restored_pool);

    // A worker fails on the pool's last launch. The coordinator stops its
    // sibling and the step throws with the partitions collected so far
    // committed, from both inputs. Resuming launches only the rest.
    PathGraph failing(left_path, left_rank);
    both_inputs(failing);
    ExternalPathJoinStats failing_stats;
    const size_type last_launch = first_partitions + second_partitions - 1;
    ::setenv("GCSA_TEST_JOIN_WORKER_FAIL_DIR", claims_root.c_str(), 1);
    ::setenv("GCSA_TEST_JOIN_WORKER_FAIL_AT", std::to_string(last_launch).c_str(), 1);
    bool interrupted = false;
    try
    {
      pool_step(failing, interrupted_root, BuildWorkspace::NEW_WORKSPACE, failing_stats);
    }
    catch(const std::runtime_error& error)
    {
      interrupted = (std::string(error.what()).find("external join worker failed") != std::string::npos);
    }
    ::unsetenv("GCSA_TEST_JOIN_WORKER_FAIL_DIR");
    ::unsetenv("GCSA_TEST_JOIN_WORKER_FAIL_AT");
    require(interrupted);
    size_type committed_first = 0, committed_second = 0;
    for(const auto& entry : workspaceEntries(interrupted_root, "step-pool-join-l"))
    {
      if(entry.first.find("--join-partition.complete") == std::string::npos) { continue; }
      if(entry.first.compare(0, 18, "step-pool-join-l7-") == 0) { committed_first++; }
      if(entry.first.compare(0, 18, "step-pool-join-l9-") == 0) { committed_second++; }
    }
    std::cerr << "test_external_join: interrupted pool committed " << committed_first
              << " + " << committed_second << " partitions" << std::endl;
    require(committed_first >= 1 && committed_second >= 1);
    require(committed_first + committed_second < pooled_stats.join_partitions);
    PathGraph resumed_pool(left_path, left_rank);
    both_inputs(resumed_pool);
    ExternalPathJoinStats resumed_pool_stats;
    pool_step(resumed_pool, interrupted_root, BuildWorkspace::RESUME, resumed_pool_stats);
    require(resumed_pool_stats.restored_partitions == committed_first + committed_second);
    require(resumed_pool_stats.worker_processes ==
      pooled_stats.join_partitions - resumed_pool_stats.restored_partitions);
    requireByteIdenticalGraph(pooled, resumed_pool);
    require(workspaceEntries(interrupted_root, "step-pool-") == expected_tasks);

    std::remove(second_left_path.c_str()); std::remove(second_left_rank.c_str());
    std::remove(second_right_path.c_str()); std::remove(second_right_rank.c_str());
    TempFile::setDirectory(std::string());
    for(const std::string& root : { alone_first_root, alone_second_root, pooled_root,
      interrupted_root, claims_root, pool_temp_root })
    {
      std::filesystem::remove_all(root);
    }
  }

  std::remove(combined_path.c_str()); std::remove(combined_rank.c_str());
  std::remove(left_path.c_str()); std::remove(left_rank.c_str());
  std::remove(right_path.c_str()); std::remove(right_rank.c_str());
  std::remove(radix_combined_path.c_str()); std::remove(radix_combined_rank.c_str());
  std::remove(radix_left_path.c_str()); std::remove(radix_left_rank.c_str());
  std::remove(radix_right_path.c_str()); std::remove(radix_right_rank.c_str());
  std::remove(changed_combined_path.c_str()); std::remove(changed_combined_rank.c_str());
  std::remove(changed_left_path.c_str()); std::remove(changed_left_rank.c_str());
  std::filesystem::remove_all(workspace_root);
  return 0;
}
