#include <gcsa/path_graph.h>
#include <gcsa/path_graph_external.h>

#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <dirent.h>
#include <fstream>
#include <map>
#include <stdexcept>
#include <sys/resource.h>
#include <limits>
#include <set>
#include <string>
#include <unistd.h>
#include <vector>

using namespace gcsa;

static void require(bool value) { if(!value) { std::abort(); } }

static void write_input(const std::string& path_name, const std::string& rank_name,
  size_type records, node_type from)
{
  std::ofstream paths(path_name.c_str(), std::ios_base::binary);
  std::ofstream ranks(rank_name.c_str(), std::ios_base::binary);
  for(size_type i = 0; i < records; i++)
  {
    PathNode node; node.from = from; node.to = i + 1; node.fields = 0;
    node.setPredecessors(1); node.setOrder(1); node.setLCP(1); node.setPointer(2 * i);
    PathNode::rank_type label[] = { 1, 9 };
    paths.write(reinterpret_cast<const char*>(&node), sizeof(node));
    ranks.write(reinterpret_cast<const char*>(label), sizeof(label));
  }
}

// A large equal-label range need not be the first range in the merge. This is
// the important late-spill shape: the in-memory window has a nonzero absolute
// offset when it is first converted into a disk-backed window.
static void write_delayed_spill_input(const std::string& path_name,
  const std::string& rank_name, size_type records)
{
  std::ofstream paths(path_name.c_str(), std::ios_base::binary);
  std::ofstream ranks(rank_name.c_str(), std::ios_base::binary);
  for(size_type i = 0; i < records; i++)
  {
    PathNode node; node.from = (i == 0 ? 7 : 42); node.to = i + 1; node.fields = 0;
    node.setPredecessors(1); node.setOrder(1); node.setLCP(1); node.setPointer(2 * i);
    PathNode::rank_type label[] = {
      static_cast<PathNode::rank_type>(i == 0 ? 0 : 1), 9
    };
    paths.write(reinterpret_cast<const char*>(&node), sizeof(node));
    ranks.write(reinterpret_cast<const char*>(label), sizeof(label));
  }
}

static std::vector<char> contents(const std::string& name)
{
  std::ifstream input(name.c_str(), std::ios_base::binary);
  return std::vector<char>((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
}

static void initialize_graph(PathGraph& result, const std::string& base,
  const std::vector<logical_file_id_t>& logical, size_type records = 80)
{
  for(size_type i = 0; i < logical.size(); i++)
  {
    write_input(base + "." + std::to_string(i) + ".path", base + "." + std::to_string(i) + ".rank", records, 42);
  }
  result.order = 1;
  result.logical_file_ids[0] = logical[0];
  result.physical_shard_ids[0] = physical_shard_id_t(100);
  for(size_type i = 1; i < logical.size(); i++)
  {
    result.path_names.push_back(base + "." + std::to_string(i) + ".path");
    result.rank_names.push_back(base + "." + std::to_string(i) + ".rank");
    result.path_counts.push_back(records); result.rank_counts.push_back(2 * records);
    result.path_count += records; result.rank_count += 2 * records;
    result.logical_file_ids.push_back(logical[i]);
    result.physical_shard_ids.push_back(physical_shard_id_t(100 + i));
  }
}

static void remove_inputs(const std::string& base, size_type files)
{
  for(size_type i = 0; i < files; i++)
  {
    std::remove((base + "." + std::to_string(i) + ".path").c_str());
    std::remove((base + "." + std::to_string(i) + ".rank").c_str());
  }
}

static PathGraphMergeStats compare_prune(const std::string& base,
  const std::vector<logical_file_id_t>& logical, const LCP& lcp,
  size_type records = 80)
{
  const std::string reference_base = base + ".reference", spilled_base = base + ".spilled";
  write_input(reference_base + ".0.path", reference_base + ".0.rank", records, 42);
  write_input(spilled_base + ".0.path", spilled_base + ".0.rank", records, 42);
  PathGraph reference(reference_base + ".0.path", reference_base + ".0.rank");
  PathGraph spilled(spilled_base + ".0.path", spilled_base + ".0.rank");
  initialize_graph(reference, reference_base, logical, records);
  initialize_graph(spilled, spilled_base, logical, records);
  reference.prune(lcp, MEGABYTE);
  PathGraphMergeStats stats;
  const size_type one_record = sizeof(PathNode) +
    (PathLabel::LABEL_LENGTH + 1) * sizeof(PathNode::rank_type);
  spilled.prune(lcp, MEGABYTE, one_record, &stats, 8);
  require(stats.priority_spills > 0);
  require(reference.files() == spilled.files());
  std::set<logical_file_id_t> distinct_logical(logical.begin(), logical.end());
  require(reference.files() == distinct_logical.size());
  for(size_type file = 0; file < reference.files(); file++)
  {
    require(reference.logicalFile(file) == spilled.logicalFile(file));
    require(reference.physicalShard(file) == spilled.physicalShard(file));
    require(contents(reference.path_names[file]) == contents(spilled.path_names[file]));
    require(contents(reference.rank_names[file]) == contents(spilled.rank_names[file]));
  }
  require(reference.size() == spilled.size());
  remove_inputs(base + ".reference", logical.size());
  remove_inputs(base + ".spilled", logical.size());
  return stats;
}

static void compare_delayed_spill(const std::string& base, const LCP& lcp)
{
  const std::string reference_path = base + ".reference.path";
  const std::string reference_rank = base + ".reference.rank";
  const std::string spilled_path = base + ".spilled.path";
  const std::string spilled_rank = base + ".spilled.rank";
  write_delayed_spill_input(reference_path, reference_rank, 80);
  write_delayed_spill_input(spilled_path, spilled_rank, 80);
  PathGraph reference(reference_path, reference_rank);
  PathGraph spilled(spilled_path, spilled_rank);
  reference.order = 1; spilled.order = 1;
  reference.prune(lcp, MEGABYTE);
  PathGraphMergeStats stats;
  const size_type one_record = sizeof(PathNode) +
    (PathLabel::LABEL_LENGTH + 1) * sizeof(PathNode::rank_type);
  spilled.prune(lcp, MEGABYTE, one_record, &stats, 8);
  require(stats.priority_spills > 0);
  require(reference.size() == spilled.size());
  require(contents(reference.path_names[0]) == contents(spilled.path_names[0]));
  require(contents(reference.rank_names[0]) == contents(spilled.rank_names[0]));
  std::remove(reference_path.c_str()); std::remove(reference_rank.c_str());
  std::remove(spilled_path.c_str()); std::remove(spilled_rank.c_str());
}

static void write_merged_fixture(const std::string& name)
{
  Alphabet alpha;
  byte_type edge = static_cast<byte_type>(1 << alpha.char2comp['A']);
  std::vector<std::string> labels = { "AAA", "AAC", "AAG", "AAT" };
  std::vector<KMer> kmers;
  for(const std::string& label : labels)
  {
    key_type key = Key::encode(alpha, label, edge, edge);
    for(size_type i = 0; i < 16; i++)
    {
      kmers.push_back(KMer(key, Node::encode(100 + i, 0), Node::encode(200 + i, 0)));
    }
  }
  std::ofstream output(name.c_str(), std::ios_base::binary);
  writeBinary(output, kmers, 3); output.close(); require(static_cast<bool>(output));
}

static void compare_merged_graph(const std::string& base)
{
  const std::string input_name = base + ".graph";
  write_merged_fixture(input_name);
  ConstructionParameters parameters;
  Alphabet alpha;
  {
    InputGraph input(std::vector<std::string>{ input_name }, true, parameters, alpha);
    std::vector<key_type> keys; input.readKeys(keys);
    require(keys.size() == 4);
    DeBruijnGraph mapper(keys, input.k(), input.alpha);
    LCP lcp(keys, input.k());
    sdsl::int_vector<0> distinct(keys.size(), 0,
      bit_length(Key::label(keys.back())));
    for(size_type i = 0; i < keys.size(); i++) { distinct[i] = Key::label(keys[i]); }
    PathGraph paths(input, distinct);

    PathGraphMergeStats reference_stats, spilled_stats;
    MergedGraph reference(paths, mapper, lcp, GIGABYTE, MEGABYTE,
      &reference_stats);
    MergedGraph spilled(paths, mapper, lcp, GIGABYTE, 1, &spilled_stats, 16);

    require(reference.size() == spilled.size());
    require(reference.ranks() == spilled.ranks());
    require(reference.extra() == spilled.extra());
    require(reference.next == spilled.next);
    require(reference.next_from == spilled.next_from);
    require(contents(reference.path_name) == contents(spilled.path_name));
    require(contents(reference.rank_name) == contents(spilled.rank_name));
    require(contents(reference.from_name) == contents(spilled.from_name));
    require(contents(reference.lcp_name) == contents(spilled.lcp_name));
    require(reference.size() == 1); require(reference.extra() == 15);
    require(reference_stats.priority_spills == 0);
    require(reference_stats.range_spills == 0);
    require(reference_stats.from_set_sorts == 0);
    require(spilled_stats.priority_spills > 0);
    require(spilled_stats.range_spills > 0);
    require(spilled_stats.from_set_sorts >= 4);
    require(spilled_stats.max_open_input_pairs <= 1);
    require(2 * spilled_stats.max_open_input_pairs + 14 <= 16);

    std::ifstream path_input(reference.path_name.c_str(), std::ios_base::binary);
    PathNode node; path_input.read(reinterpret_cast<char*>(&node), sizeof(node));
    require(path_input.gcount() == sizeof(node));
    require(node.from == Node::encode(100, 0));
    std::ifstream from_input(reference.from_name.c_str(), std::ios_base::binary);
    for(size_type i = 1; i < 16; i++)
    {
      range_type extra;
      from_input.read(reinterpret_cast<char*>(&extra), sizeof(extra));
      require(from_input.gcount() == sizeof(extra));
      require(extra.first == 0); require(extra.second == Node::encode(100 + i, 0));
    }
  }
  std::remove(input_name.c_str());
}

struct PruneFixtureRecord
{
  PathNode::rank_type rank;
  node_type from;
  bool sorted;
};

static std::vector<key_type> parallel_prune_keys(Alphabet& alpha)
{
  const std::vector<std::string> labels = {
    "AAA", "AAC", "CAA", "CAC", "GAA", "GAC", "TAA", "TAC"
  };
  std::vector<key_type> result;
  for(const std::string& label : labels)
  {
    comp_type comp = alpha.char2comp[label.front()];
    byte_type edge = static_cast<byte_type>(1 << comp);
    result.push_back(Key::encode(alpha, label, edge, edge));
  }
  std::sort(result.begin(), result.end());
  return result;
}

static std::vector<std::vector<PruneFixtureRecord>> parallel_prune_records()
{
  std::vector<std::vector<PruneFixtureRecord>> result(3);
  auto add = [&](size_type file, PathNode::rank_type rank, size_type count,
    node_type from, bool sorted = false)
  {
    for(size_type i = 0; i < count; i++)
    {
      result[file].push_back({ rank, from, sorted && i == 0 });
    }
  };

  // Ranks 0/1 share a nonzero-LCP root component, one semantic input, and a
  // start node, so extendRange() must join their groups across two physical
  // shards. Rank 2 is duplicated across logical inputs. Rank 3 mixes sorted
  // and unsorted paths. Ranks 4/5 extend, while rank 6 has the same start and
  // logical identity but sits beyond an exact zero-LCP boundary and must not.
  add(0, 0, 5000, Node::encode(100, 0));
  add(1, 0, 4000, Node::encode(100, 0));
  add(0, 1, 6, Node::encode(100, 0));
  add(1, 1, 4, Node::encode(100, 0));
  add(0, 2, 3, Node::encode(200, 0));
  add(2, 2, 3, Node::encode(200, 0));
  add(0, 3, 2, Node::encode(301, 0), true);
  add(1, 3, 2, Node::encode(302, 0));
  add(1, 4, 8, Node::encode(500, 0));
  add(1, 5, 5, Node::encode(500, 0));
  add(1, 6, 5, Node::encode(500, 0));
  add(2, 6, 2, Node::encode(600, 0));
  add(1, 7, 3, Node::encode(500, 0));
  add(2, 7, 2, Node::encode(601, 0));
  for(auto& records : result)
  {
    std::stable_sort(records.begin(), records.end(),
      [](const PruneFixtureRecord& left, const PruneFixtureRecord& right)
      {
        return left.rank < right.rank;
      });
  }
  return result;
}

static void write_parallel_prune_pair(const std::string& path_name,
  const std::string& rank_name, const std::vector<PruneFixtureRecord>& records,
  bool corrupt = false, bool invalid_tail_rank = false)
{
  std::ofstream paths(path_name, std::ios::binary);
  std::ofstream ranks(rank_name, std::ios::binary);
  for(size_type i = 0; i < records.size(); i++)
  {
    PathNode node;
    node.from = records[i].from; node.to = Node::encode(10000 + i, 0);
    node.fields = 0; node.setPredecessors(1);
    node.setOrder(1); node.setLCP(1); node.setPointer(2 * i);
    if(records[i].sorted) { node.makeSorted(); }
    if(corrupt && i == 0) { node.setPointer(2 * records.size() + 17); }
    PathNode::rank_type label[] = {
      (invalid_tail_rank && i + 1 == records.size() ? 12345 : records[i].rank), 0
    };
    paths.write(reinterpret_cast<const char*>(&node), sizeof(node));
    ranks.write(reinterpret_cast<const char*>(label), sizeof(label));
  }
  paths.close(); ranks.close();
  require(static_cast<bool>(paths) && static_cast<bool>(ranks));
}

static void initialize_parallel_prune_graph(PathGraph& graph,
  const std::string& base, bool corrupt = false, bool invalid_tail_rank = false)
{
  auto records = parallel_prune_records(); records.emplace_back();
  const logical_file_id_t logical[] = {
    logical_file_id_t(7), logical_file_id_t(7), logical_file_id_t(8), logical_file_id_t(9)
  };
  for(size_type file = 0; file < records.size(); file++)
  {
    std::string path = base + "." + std::to_string(file) + ".path";
    std::string rank = base + "." + std::to_string(file) + ".rank";
    write_parallel_prune_pair(path, rank, records[file], corrupt && file == 0,
      invalid_tail_rank && file == 0);
    graph.path_names.push_back(path); graph.rank_names.push_back(rank);
    graph.path_counts.push_back(records[file].size());
    graph.rank_counts.push_back(2 * records[file].size());
    graph.logical_file_ids.push_back(logical[file]);
    graph.physical_shard_ids.push_back(physical_shard_id_t(100 + file));
    graph.path_count += records[file].size(); graph.rank_count += 2 * records[file].size();
  }
}

struct LogicalPruneRecord
{
  node_type from, to;
  size_type fields;
  std::vector<PathNode::rank_type> label;

  bool operator==(const LogicalPruneRecord& another) const
  {
    return (this->from == another.from && this->to == another.to &&
      this->fields == another.fields && this->label == another.label);
  }
};

static std::map<logical_file_id_t, std::vector<LogicalPruneRecord>>
logical_prune_records(const PathGraph& graph)
{
  std::map<logical_file_id_t, std::vector<LogicalPruneRecord>> result;
  for(size_type file = 0; file < graph.files(); file++)
  {
    std::vector<PathNode> paths;
    std::vector<PathNode::rank_type> labels;
    graph.read(paths, labels, file);
    for(PathNode node : paths)
    {
      size_type pointer = node.pointer();
      require(pointer <= labels.size());
      require(node.ranks() <= labels.size() - pointer);
      LogicalPruneRecord record;
      record.from = node.from; record.to = node.to;
      record.label.assign(labels.begin() + pointer,
        labels.begin() + pointer + node.ranks());
      node.setPointer(0); record.fields = node.fields;
      result[graph.logicalFile(file)].push_back(record);
    }
  }
  return result;
}

static size_type directory_entries(const std::string& directory)
{
  DIR* stream = ::opendir(directory.c_str()); require(stream != nullptr);
  size_type result = 0;
  while(struct dirent* entry = ::readdir(stream))
  {
    std::string name(entry->d_name);
    if(name != "." && name != "..") { result++; }
  }
  require(::closedir(stream) == 0);
  return result;
}

static void compare_parallel_prune(const std::string& root)
{
  Alphabet alpha;
  auto keys = parallel_prune_keys(alpha);
  LCP lcp(keys, 3); DeBruijnGraph mapper(keys, 3, alpha);
  const size_type budget = 512 * KILOBYTE;
  PathGraph serial(0, 3, 0);
  initialize_parallel_prune_graph(serial, root + "/serial");
  serial.prune(lcp, GIGABYTE, budget);
  MergedGraph serial_final(serial, mapper, lcp, GIGABYTE, MEGABYTE);
  for(size_type workers : { 1, 2, 4 })
  {
    for(size_type repetition = 0; repetition < 3; repetition++)
    {
      PathGraph candidate(0, 3, 0);
      initialize_parallel_prune_graph(candidate, root + "/candidate");
      PathGraphMergeStats stats;
      externalPathGraphPrune(candidate, lcp, GIGABYTE, budget, 128, workers, &stats);
      require(logical_prune_records(serial) == logical_prune_records(candidate));
      require(candidate.size() == serial.size()); require(candidate.ranks() == serial.ranks());
      require(candidate.ranges() == serial.ranges()); require(candidate.unique == serial.unique);
      require(candidate.redundant == serial.redundant);
      require(candidate.unsorted == serial.unsorted);
      require(candidate.nondeterministic == serial.nondeterministic);
      std::set<logical_file_id_t> logical(candidate.logical_file_ids.begin(), candidate.logical_file_ids.end());
      require(logical == std::set<logical_file_id_t>{ logical_file_id_t(7), logical_file_id_t(8), logical_file_id_t(9) });
      if(workers > 1) { require(candidate.files() > serial.files()); require(stats.priority_spills > 0); }
      else { require(candidate.files() == serial.files()); }
      MergedGraph merged(candidate, mapper, lcp, GIGABYTE, MEGABYTE);
      require(contents(serial_final.path_name) == contents(merged.path_name));
      require(contents(serial_final.rank_name) == contents(merged.rank_name));
      require(contents(serial_final.from_name) == contents(merged.from_name));
      require(contents(serial_final.lcp_name) == contents(merged.lcp_name));
    }
  }
  // Both a small buffer allowance and a small descriptor allowance select
  // the exact serial route before any worker output can be adopted.
  for(bool small_memory : { false, true })
  {
    PathGraph fallback(0, 3, 0); initialize_parallel_prune_graph(fallback, root + "/fallback");
    externalPathGraphPrune(fallback, lcp, GIGABYTE,
      (small_memory ? 16 * KILOBYTE : budget), (small_memory ? 128 : 10), 4);
    require(logical_prune_records(serial) == logical_prune_records(fallback));
    require(fallback.files() == serial.files());
  }
  for(bool invalid_tail : { false, true })
  {
    PathGraph broken(0, 3, 0);
    initialize_parallel_prune_graph(broken, root + "/broken", !invalid_tail, invalid_tail);
    const auto paths = broken.path_names; const size_type count = broken.size();
    const size_type before = directory_entries(root);
    bool failed = false;
    try { externalPathGraphPrune(broken, lcp, GIGABYTE, budget, 128, 4); }
    catch(const std::runtime_error&) { failed = true; }
    require(failed); require(broken.path_names == paths); require(broken.size() == count);
    require(directory_entries(root) == before);
    for(const auto& path : paths) { require(::access(path.c_str(), F_OK) == 0); }
  }
  PathGraph empty(1, 3, 0); LCP empty_lcp;
  { std::ofstream paths(empty.path_names[0]), ranks(empty.rank_names[0]); }
  externalPathGraphPrune(empty, empty_lcp, GIGABYTE, budget, 128, 4);
  require(empty.size() == 0 && empty.files() == 1);
}

static void compare_unsafe_prune_boundary(const std::string& root)
{
  Alphabet alpha; auto keys = parallel_prune_keys(alpha); LCP lcp(keys, 3);
  PathGraph serial(1, 3, 0), candidate(1, 3, 0);
  for(PathGraph* graph : { &serial, &candidate })
  {
    std::ofstream paths(graph->path_names[0], std::ios::binary);
    std::ofstream ranks(graph->rank_names[0], std::ios::binary);
    for(size_type i = 0; i < 2; i++)
    {
      PathNode node; node.from = Node::encode(10 + i, 0); node.to = Node::encode(20 + i, 0);
      node.fields = 0; node.setOrder(1); node.setLCP(i == 0 ? 0 : 1);
      node.setPredecessors(1); node.setPointer(2 * i);
      PathNode::rank_type label[] = { static_cast<PathNode::rank_type>(2 * i), 2 };
      paths.write(reinterpret_cast<const char*>(&node), sizeof(node));
      ranks.write(reinterpret_cast<const char*>(label), sizeof(label));
    }
    graph->path_count = graph->path_counts[0] = 2;
    graph->rank_count = graph->rank_counts[0] = 4;
  }
  serial.prune(lcp, MEGABYTE, MEGABYTE);
  externalPathGraphPrune(candidate, lcp, MEGABYTE, MEGABYTE, 128, 4);
  require(candidate.files() == serial.files());
  require(logical_prune_records(serial) == logical_prune_records(candidate));
}

int main()
{
  const char* configured_tmp = std::getenv("TMPDIR");
  const std::string temp_root =
    ((configured_tmp != nullptr && configured_tmp[0] != 0) ? configured_tmp : "/tmp");
  const std::string root_template = temp_root + "/gcsa-prune-spill-XXXXXX";
  std::vector<char> root_buffer(root_template.begin(), root_template.end());
  root_buffer.push_back(0);
  require(mkdtemp(root_buffer.data()) != nullptr);
  const std::string root(root_buffer.data()); TempFile::setDirectory(root);
  std::vector<key_type> keys(16); for(size_type i = 0; i < keys.size(); i++) { keys[i] = i; }
  LCP lcp(keys, 1);

  // One large equal-label group exercises the spill path against the former in-memory result.
  compare_prune(std::string(root) + "/one", { logical_file_id_t(7) }, lcp);
  // First consume a small range, then spill. Absolute merger offsets must not
  // be confused with offsets relative to a newly created spill file.
  compare_delayed_spill(std::string(root) + "/delayed", lcp);
  // Physical shards 0/1 represent logical input 7; shard 2 is logical input 8.
  // The spill path must retain exactly the same multi-logical-file output semantics.
  compare_prune(std::string(root) + "/many", { logical_file_id_t(7), logical_file_id_t(7), logical_file_id_t(8) }, lcp);
  // Reader threads, MiB buffers, and descriptors must not scale with the
  // number of physical shards. Current-head metadata remains one per shard.
  std::vector<logical_file_id_t> shards(40, logical_file_id_t(7));
  PathGraphMergeStats shard_stats = compare_prune(std::string(root) + "/shards",
    shards, lcp, 4);
  require(shard_stats.max_open_input_pairs <= 2);
  require(shard_stats.max_open_input_pairs < shards.size());
  require(shard_stats.max_open_output_pairs <= 2);
  require(shard_stats.max_open_output_pairs < shards.size());
  require(2 * shard_stats.max_open_input_pairs +
    2 * shard_stats.max_open_output_pairs + 2 <= 8);

  // Use the real construction route so mapper ranks and LCP invariants are
  // validated rather than synthesized by the test.
  compare_merged_graph(std::string(root) + "/merged");
  compare_parallel_prune(root);
  compare_unsafe_prune_boundary(root);
  require(directory_entries(root) == 0);
  require(rmdir(root.c_str()) == 0);
  return 0;
}
