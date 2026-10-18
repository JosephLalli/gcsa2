#include <gcsa/path_graph.h>
#include <gcsa/path_graph_external.h>
#include <gcsa/checkpoint.h>
#include <gcsa/compressed_block.h>

#include <algorithm>
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
#include <sys/resource.h>
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

int main()
{
  omp_set_num_threads(1);

  // The disk ceiling applies to the peak physical framed pair, not merely its
  // logical PathNode/rank payload. The sink must report rejection without
  // terminating the process, so threaded task cleanup can unwind normally.
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
    bool rejected = false;
    try
    {
      size_type committed = 0;
      ExternalPathSortSink sink(rejected_output, 0, limit_sort_budget, 2,
        logical_pair_bytes, committed, nullptr, limit_codec);
      sink.write(limit_node, limit_ranks);
    }
    catch(const std::runtime_error& error)
    {
      rejected = (std::string(error.what()).find(
        "configured disk limit exceeded") != std::string::npos);
    }
    require(rejected);
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

  // The same skewed key must be splittable across fork-free join threads.
  // Every physical result shard retains logical ID 7, and the PathGraph merger
  // must observe exactly the same semantic record multiset as the legacy path.
  PathGraph threaded_graph(left_path, left_rank);
  threaded_graph.order = 1;
  threaded_graph.logical_file_ids[0] = logical_file_id_t(7);
  threaded_graph.physical_shard_ids[0] = physical_shard_id_t(101);
  appendShard(threaded_graph, right_path, right_rank,
    logical_file_id_t(7), physical_shard_id_t(202));
  ConstructionParameters threaded_parameters = parameters;
  threaded_parameters.setMemoryLimitBytes(8 * MEGABYTE);
  // Four MiB admits the framed run writer plus two framed run readers and,
  // on the next doubling step, the compressed source-shard decoder pair.
  threaded_parameters.setJoinPartitionSize(4 * MEGABYTE);
  threaded_parameters.setSortRunSize(2 * MEGABYTE);
  threaded_parameters.setCheckpointBytes(32 * KILOBYTE);
  threaded_parameters.setTempCompression("zstd");
  // 1 KiB frames force both the 40-summary stream and the recursively
  // split 1201-detail stream across several blocks. Planning consequently
  // exercises framed readAt() for nonzero subgroup offsets, rather than only
  // a sequential first-block read.
  threaded_parameters.setCompressionBlockSize(KILOBYTE);
  threaded_parameters.setCompressionWorkers(1);
  // The setter clamps to MIN_OPEN_FILES (64). Compaction retains a bounded
  // multi-run frontier rather than pretending this request limits it to 8.
  threaded_parameters.setMaxOpenFiles(8);
  omp_set_num_threads(4);
  ExternalPathJoinStats threaded_stats;
  externalPathGraphExtend(threaded_graph, GIGABYTE,
    threaded_parameters, &threaded_stats, &workspace, "step-threaded");
  omp_set_num_threads(1);
  require(readGraph(legacy) == readGraph(threaded_graph));
  // Many task partitions remain a bounded set of individually LABEL-sorted
  // streams for the downstream merger.
  require(threaded_graph.files() > 1);
  require(threaded_graph.files() <= (threaded_parameters.getMaxOpenFiles() - 4) / 2);
  require(labelSorted(threaded_graph));
  for(size_type file = 0; file < threaded_graph.files(); file++)
  {
    require(threaded_graph.logicalFile(file) == logical_file_id_t(7));
    require(threaded_graph.physicalShard(file) == physical_shard_id_t(file));
    require(CompressedBlockReader::isFramed(threaded_graph.path_names[file]));
    require(CompressedBlockReader::isFramed(threaded_graph.rank_names[file]));
  }
  require(threaded_stats.join_partitions > 1);
  // The sampled planner must turn this single high-fanout semantic key into
  // several deterministic range tasks rather than one task monopoly.
  require(threaded_stats.join_partitions >= 4);
  require(threaded_stats.join_tasks == threaded_stats.join_partitions);
  require(threaded_stats.join_concurrency > 1);
  require(threaded_stats.recursive_splits > 0);
  require(threaded_stats.left_range_splits > 0);
  require(threaded_stats.right_range_splits > 0);
  require(threaded_stats.generated_records == left.size() * right.size());
  require(threaded_stats.sorted_bypass == 2);
  require(threaded_stats.sampled_plan_records > 0);
  require(threaded_stats.radix_plan_bins >= 1);
  require(threaded_stats.sidecar_plan_groups > 0);
  require(threaded_stats.sidecar_plan_detail_records > 0);
  require(threaded_stats.sidecar_plan_groups * 48 > KILOBYTE);
  require(threaded_stats.sidecar_plan_detail_records > KILOBYTE);
  require(threaded_stats.full_record_plan_rescans == 0);
  require(threaded_stats.compressed_join_runs > 0);
  require(threaded_stats.join_run_stored_bytes <
    threaded_stats.join_run_logical_bytes);
  require(threaded_stats.compressed_join_sidecars > 0);
  require(threaded_stats.join_sidecar_stored_bytes <
    threaded_stats.join_sidecar_logical_bytes);
  require(threaded_stats.grouped_expansion_records > 0);
  require(threaded_stats.expansion_context_bytes_saved > 0);
  require(threaded_stats.max_bytes_resident <= threaded_parameters.getMemoryLimitBytes());

  // Rebuild the cross-file, high-fanout threaded plan without a workspace. The
  // compact sidecars must produce the same physical streams, not merely the
  // same logical record multiset.
  PathGraph fresh_threaded(left_path, left_rank);
  fresh_threaded.order = 1;
  fresh_threaded.logical_file_ids[0] = logical_file_id_t(7);
  fresh_threaded.physical_shard_ids[0] = physical_shard_id_t(101);
  appendShard(fresh_threaded, right_path, right_rank,
    logical_file_id_t(7), physical_shard_id_t(202));
  ConstructionParameters fresh_parameters = threaded_parameters;
  fresh_parameters.setVerifyWorkspace();
  ExternalPathJoinStats fresh_stats;
  omp_set_num_threads(4);
  externalPathGraphExtend(fresh_threaded, GIGABYTE,
    fresh_parameters, &fresh_stats);
  omp_set_num_threads(1);
  requireSameGraph(threaded_graph, fresh_threaded, "fresh sidecar threaded plan");
  requireByteIdenticalGraph(threaded_graph, fresh_threaded);
  require(fresh_stats.full_record_plan_rescans == 0);

  PathGraph resumed_threaded(left_path, left_rank);
  resumed_threaded.order = 1;
  resumed_threaded.logical_file_ids[0] = logical_file_id_t(7);
  resumed_threaded.physical_shard_ids[0] = physical_shard_id_t(101);
  appendShard(resumed_threaded, right_path, right_rank,
    logical_file_id_t(7), physical_shard_id_t(202));
  ExternalPathJoinStats resumed_stats;
  omp_set_num_threads(4);
  externalPathGraphExtend(resumed_threaded, GIGABYTE,
    threaded_parameters, &resumed_stats, &workspace, "step-threaded");
  omp_set_num_threads(1);
  require(readGraph(threaded_graph) == readGraph(resumed_threaded));
  requireByteIdenticalGraph(threaded_graph, resumed_threaded);
  require(resumed_threaded.files() == threaded_graph.files());
  require(labelSorted(resumed_threaded));
  for(size_type file = 0; file < resumed_threaded.files(); file++)
  {
    require(resumed_threaded.physicalShard(file) == physical_shard_id_t(file));
  }
  require(resumed_stats.join_tasks == 0);
  require(resumed_stats.join_concurrency == 0);
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
  ConstructionParameters second_parameters = threaded_parameters;
  ExternalPathJoinStats second_stats;
  externalPathGraphExtend(resumed_threaded, GIGABYTE, second_parameters,
    &second_stats);
  requireSameGraph(second_legacy, resumed_threaded,
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
  ConstructionParameters radix_parameters = threaded_parameters;
  radix_parameters.setCheckpointBytes(4 * KILOBYTE);
  ExternalPathJoinStats radix_stats;
  omp_set_num_threads(4);
  externalPathGraphExtend(radix_external, GIGABYTE, radix_parameters,
    &radix_stats, &workspace, "step-msd");
  omp_set_num_threads(1);
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
  omp_set_num_threads(4);
  externalPathGraphExtend(radix_resumed, GIGABYTE, radix_parameters,
    &radix_resumed_stats, &workspace, "step-msd");
  omp_set_num_threads(1);
  require(readGraph(radix_external) == readGraph(radix_resumed));
  requireByteIdenticalGraph(radix_external, radix_resumed);
  require(radix_resumed_stats.restored_radix_plans == 1);
  require(radix_resumed_stats.restored_sidecar_metadata == 2);
  require(radix_resumed_stats.join_tasks == 0);

  // A same-length run with different payload must not inherit either the plan
  // or completed task ranges. Run checksums are part of both task identities.
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
  omp_set_num_threads(4);
  externalPathGraphExtend(changed_external, GIGABYTE, radix_parameters,
    &changed_stats, &workspace, "step-msd");
  omp_set_num_threads(1);
  require(readGraph(changed_legacy) == readGraph(changed_external));
  require(changed_stats.restored_radix_plans == 0);
  require(changed_stats.restored_partitions == 0);
  require(changed_stats.join_tasks == changed_stats.join_partitions);
  require(changed_stats.join_concurrency > 1);

  // Requesting multiple threads must remain valid when the deterministic plan
  // contains only one task. This is common for small logical inputs.
  PathGraph one_partition(left_path, left_rank);
  one_partition.order = 1;
  one_partition.logical_file_ids[0] = logical_file_id_t(7);
  one_partition.physical_shard_ids[0] = physical_shard_id_t(101);
  appendShard(one_partition, right_path, right_rank,
    logical_file_id_t(7), physical_shard_id_t(202));
  ConstructionParameters one_parameters = threaded_parameters;
  one_parameters.setCheckpointBytes(GIGABYTE);
  ExternalPathJoinStats one_stats;
  omp_set_num_threads(4);
  externalPathGraphExtend(one_partition, GIGABYTE, one_parameters, &one_stats);
  omp_set_num_threads(1);
  require(readGraph(legacy) == readGraph(one_partition));
  require(one_partition.files() == 1);
  require(one_stats.join_partitions == 1);
  require(one_stats.join_tasks == 1);
  require(one_stats.join_concurrency == 1);

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

  // Join tasks leave many physical shards per logical input, and compaction
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
    MemoryBudget serial_memory(compaction_parameters.getMemoryLimitBytes());
    MemoryBudget concurrent_memory(compaction_parameters.getMemoryLimitBytes());
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
    MemoryBudget memory(limited_parameters.getMemoryLimitBytes());
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
    MemoryBudget budgeted_memory(limited_parameters.getMemoryLimitBytes());
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
  // Every route must produce the serial route's records; the ranged threaded
  // task route must also restore its range plan and every partition when
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

    // Threaded joins without a workspace. The checkpoint-sized task shards
    // are compacted to the bounded frontier the downstream merger can open.
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
    require(direct_chain_stats.join_tasks == direct_chain_stats.join_partitions);
    require(direct_chain_stats.join_concurrency > 1);
    require(direct_chain_stats.join_partitions > direct_chain.files());
    require(direct_chain_stats.compaction_batches > 0);
    require(direct_chain.files() <= pathMergeInputPairs(range_parameters));
    for(size_type file = 0; file < direct_chain.files(); file++)
    {
      require(direct_chain.logicalFile(file) == logical_file_id_t(3));
      require(direct_chain.physicalShard(file) == physical_shard_id_t(file));
    }
    requireSameGraph(serial_chain, direct_chain, "key-range distribution, threaded joins");

    // Threaded tasks over every range's partitions, on a workspace.
    char range_root[] = "/tmp/gcsa-range-plan-XXXXXX";
    require(mkdtemp(range_root) != nullptr);
    BuildWorkspace::Settings range_semantic;
    range_semantic["fixture"] = "external-join-ranges";
    BuildWorkspace range_workspace(range_root, range_semantic, BuildWorkspace::Settings(),
      BuildWorkspace::NEW_WORKSPACE);
    ConstructionParameters workspace_range_parameters = range_parameters;
    PathGraph workspace_chain(chain_path, chain_rank);
    chain_graph(workspace_chain);
    ExternalPathJoinStats workspace_chain_stats;
    omp_set_num_threads(4);
    externalPathGraphExtend(workspace_chain, GIGABYTE, workspace_range_parameters,
      &workspace_chain_stats, &range_workspace, "step-ranges");
    omp_set_num_threads(1);
    require(workspace_chain_stats.distribution_ranges == 4);
    require(workspace_chain_stats.distribution_concurrency == 4); // one reader per thread, in chunks
    require(workspace_chain_stats.restored_range_plans == 0);
    require(workspace_chain_stats.join_partitions >= 4);
    require(workspace_chain_stats.join_tasks == workspace_chain_stats.join_partitions);
    require(workspace_chain_stats.join_concurrency > 1);
    requireSameGraph(serial_chain, workspace_chain, "key-range distribution, workspace tasks");
    requireByteIdenticalGraph(direct_chain, workspace_chain);

    // The same step again: the committed range plan gives the same ranges
    // whatever the thread count, so the runs, their checksums and every
    // partition task name match, and nothing is joined again.
    PathGraph resumed_chain(chain_path, chain_rank);
    chain_graph(resumed_chain);
    ExternalPathJoinStats resumed_chain_stats;
    omp_set_num_threads(2);
    externalPathGraphExtend(resumed_chain, GIGABYTE, workspace_range_parameters,
      &resumed_chain_stats, &range_workspace, "step-ranges");
    omp_set_num_threads(1);
    require(resumed_chain_stats.restored_range_plans == 1);
    require(resumed_chain_stats.distribution_ranges == 4);
    require(resumed_chain_stats.distribution_concurrency == 2); // one reader per thread, in chunks
    require(resumed_chain_stats.restored_radix_plans > 0);
    require(resumed_chain_stats.restored_partitions == workspace_chain_stats.join_partitions);
    require(resumed_chain_stats.join_tasks == 0);
    require(resumed_chain_stats.join_concurrency == 0);
    requireByteIdenticalGraph(workspace_chain, resumed_chain);

    std::remove(chain_path.c_str()); std::remove(chain_rank.c_str());
    std::filesystem::remove_all(range_root);
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
