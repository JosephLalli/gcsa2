#include <gcsa/path_graph.h>

#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <sys/wait.h>
#include <type_traits>
#include <unistd.h>
#include <vector>

using namespace gcsa;

// Storage topology must not be usable as graph semantics by accident.
static_assert(!std::is_convertible<physical_shard_id_t, logical_file_id_t>::value,
  "physical shard IDs must not convert to logical input IDs");
static_assert(!std::is_convertible<logical_file_id_t, physical_shard_id_t>::value,
  "logical input IDs must not convert to physical shard IDs");

static void require(bool condition)
{
  if(!condition) { std::abort(); }
}

struct TestRecord
{
  node_type from, to;
  byte_type predecessors;
  size_type order, lcp;
  std::vector<PathNode::rank_type> labels;
};

static void write_input(const std::string& paths, const std::string& ranks,
  const std::vector<TestRecord>& records)
{
  std::ofstream path_out(paths.c_str(), std::ios_base::binary);
  std::ofstream rank_out(ranks.c_str(), std::ios_base::binary);
  size_type pointer = 0;
  for(size_type i = 0; i < records.size(); i++)
  {
    const TestRecord& source = records[i];
    require(source.labels.size() == source.order + 1);
    PathNode node;
    node.from = source.from; node.to = source.to; node.fields = 0;
    node.setPredecessors(source.predecessors); node.setOrder(source.order);
    node.setLCP(source.lcp); node.setPointer(pointer);
    path_out.write(reinterpret_cast<const char*>(&node), sizeof(node));
    rank_out.write(reinterpret_cast<const char*>(source.labels.data()),
      source.labels.size() * sizeof(PathNode::rank_type));
    pointer += source.labels.size();
  }
}

static std::vector<PathNode::rank_type> read_ranks(const std::string& name)
{
  std::ifstream input(name.c_str(), std::ios_base::binary);
  std::vector<PathNode::rank_type> result;
  PathNode::rank_type value;
  while(input.read(reinterpret_cast<char*>(&value), sizeof(value))) { result.push_back(value); }
  return result;
}

static bool test_record_less(const TestRecord& a, const TestRecord& b)
{
  size_type order = std::min(a.order, b.order);
  for(size_type i = 0; i < order; i++)
  {
    if(a.labels[i] != b.labels[i]) { return (a.labels[i] < b.labels[i]); }
  }
  if(a.order != b.order) { return (a.order < b.order); }
  if(a.from != b.from) { return (a.from < b.from); }
  if(a.to != b.to) { return (a.to < b.to); }
  if(a.predecessors != b.predecessors) { return (a.predecessors < b.predecessors); }
  if(a.lcp != b.lcp) { return (a.lcp < b.lcp); }
  return a.labels < b.labels;
}

static void check_sorted(PathGraph& graph, std::vector<TestRecord> expected)
{
  std::sort(expected.begin(), expected.end(), test_record_less);
  require(graph.files() == 1);
  std::vector<PathNode> paths;
  std::vector<PathNode::rank_type> ranks;
  graph.read(paths, ranks, 0);
  require(paths.size() == expected.size());
  size_type pointer = 0;
  for(size_type i = 0; i < paths.size(); i++)
  {
    require(paths[i].from == expected[i].from);
    require(paths[i].to == expected[i].to);
    require(paths[i].predecessors() == expected[i].predecessors);
    require(paths[i].order() == expected[i].order);
    require(paths[i].lcp() == expected[i].lcp);
    require(paths[i].pointer() == pointer);
    for(size_type j = 0; j < paths[i].ranks(); j++)
    {
      require(ranks[pointer + j] == expected[i].labels[j]);
    }
    pointer += paths[i].ranks();
  }
  require(pointer == ranks.size());
}

static void expect_sort_failure(const std::string& path, const std::string& rank,
  size_type budget, bool invalid_pointer, bool trailing_rank)
{
  pid_t child = fork();
  require(child >= 0);
  if(child == 0)
  {
    PathGraph graph(path, rank);
    if(invalid_pointer)
    {
      std::fstream input(path.c_str(), std::ios_base::in | std::ios_base::out | std::ios_base::binary);
      PathNode node; input.read(reinterpret_cast<char*>(&node), sizeof(node));
      node.setPointer(1); input.seekp(0); input.write(reinterpret_cast<const char*>(&node), sizeof(node));
      input.close();
    }
    if(trailing_rank)
    {
      std::ofstream output(rank.c_str(), std::ios_base::binary | std::ios_base::app);
      PathNode::rank_type extra = 999; output.write(reinterpret_cast<const char*>(&extra), sizeof(extra));
    }
    externalPathGraphSort(graph, 0, budget, 2);
    _exit(0);
  }
  int status = 0; require(waitpid(child, &status, 0) == child);
  require(!WIFEXITED(status) || WEXITSTATUS(status) != 0);
}

int main()
{
  const std::vector<TestRecord> pattern = {
    { 8, 20, 3, 1, 1, { 2, 90 } },
    { 7, 21, 2, 2, 1, { 1, 5, 80 } },
    { 9, 22, 1, 1, 1, { 1, 70 } },
    { 2, 23, 1, 1, 1, { 2, 60 } },
    { 6, 24, 2, 2, 1, { 1, 5, 70 } },
    { 1, 25, 1, 2, 1, { 1, 5, 70 } }
  };
  std::vector<TestRecord> records;
  for(size_type copy = 0; copy < 6; copy++)
  {
    for(size_type i = 0; i < pattern.size(); i++)
    {
      TestRecord record = pattern[i];
      record.from += 100 * copy; record.to += 100 * copy;
      records.push_back(record);
    }
  }
  const std::string base = "test_external_path_sort_" + std::to_string((unsigned long long)getpid());
  const std::string path_a = base + ".a.path", rank_a = base + ".a.rank";
  const std::string path_b = base + ".b.path", rank_b = base + ".b.rank";
  const std::string path_c = base + ".c.path", rank_c = base + ".c.rank";
  const std::string path_d = base + ".d.path", rank_d = base + ".d.rank";
  write_input(path_a, rank_a, records); write_input(path_b, rank_b, records);
  write_input(path_c, rank_c, records); write_input(path_d, rank_d, records);

  PathGraph first(path_a, rank_a);
  first.logical_file_ids[0] = logical_file_id_t(17);
  first.physical_shard_ids[0] = physical_shard_id_t(9001);
  ExternalPathSortStats first_stats;
  size_type budget = externalPathGraphSortMinimumBudget();
  externalPathGraphSort(first, 0, budget, 2, &first_stats);
  check_sorted(first, records);
  require(first_stats.runs > 2);
  require(first_stats.merge_passes >= 2);
  require(first_stats.max_records_resident < records.size());
  require(first_stats.max_bytes_resident <= budget);
  require(first.logicalFile(0) == logical_file_id_t(17));
  require(first.physicalShard(0) == physical_shard_id_t(9001));
  std::vector<PathNode::rank_type> first_ranks = read_ranks(first.rank_names[0]);

  PathGraph second(path_b, rank_b);
  externalPathGraphSort(second, 0, budget, 2);
  check_sorted(second, records);
  require(first_ranks == read_ranks(second.rank_names[0]));
  require(first.files() == 1 && second.files() == 1);

  expect_sort_failure(path_c, rank_c, budget - 1, false, false);
  expect_sort_failure(path_c, rank_c, budget, true, false);
  expect_sort_failure(path_d, rank_d, budget, false, true);

  std::remove(first.path_names[0].c_str()); std::remove(first.rank_names[0].c_str());
  std::remove(second.path_names[0].c_str()); std::remove(second.rank_names[0].c_str());
  std::remove(path_c.c_str()); std::remove(rank_c.c_str());
  std::remove(path_d.c_str()); std::remove(rank_d.c_str());
  return 0;
}
