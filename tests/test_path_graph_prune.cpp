#include <gcsa/path_graph.h>

#include <cstdlib>
#include <cstdio>
#include <fstream>
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
    // In-memory SameFromSet groups must not invoke the external sorter.
    require(reference_stats.from_set_sorts == 0);
    require(spilled_stats.priority_spills > 0);
    require(spilled_stats.range_spills > 0);
    // A tiny group budget cannot retain even one selected/candidate pair, so
    // it takes the external path while producing byte-identical output above.
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

int main()
{
  char root[] = "/tmp/gcsa-prune-spill-XXXXXX";
  require(mkdtemp(root) != nullptr); TempFile::setDirectory(root);
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
  rmdir(root);
  return 0;
}
