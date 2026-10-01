#include <gcsa/path_graph.h>
#include <gcsa/path_graph_external.h>
#include <gcsa/compressed_block.h>
#include <gcsa/support.h>

#include <cstdlib>
#include <cstdio>
#include <dirent.h>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <sys/resource.h>
#include <unistd.h>
#include <vector>

using namespace gcsa;

static void require(bool value) { if(!value) { std::abort(); } }

// A framed reader must consume one descriptor, not an istream descriptor plus
// a second pread/cache-advice descriptor. Lowering the process ceiling makes
// that resource contract executable: the 40-pair fixture below fits with one
// descriptor per stream and deterministically fails with two.
class ScopedFileLimit
{
public:
  explicit ScopedFileLimit(rlim_t ceiling) : original(), changed(false)
  {
    require(::getrlimit(RLIMIT_NOFILE, &(this->original)) == 0);
    if(this->original.rlim_cur > ceiling)
    {
      struct rlimit limited = this->original;
      limited.rlim_cur = ceiling;
      require(::setrlimit(RLIMIT_NOFILE, &limited) == 0);
      this->changed = true;
    }
  }

  ~ScopedFileLimit()
  {
    if(this->changed)
    {
      require(::setrlimit(RLIMIT_NOFILE, &(this->original)) == 0);
    }
  }

private:
  struct rlimit original;
  bool changed;
};

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

enum PruneFixtureStorage { FIXTURE_RAW, FIXTURE_FRAMED, FIXTURE_MIXED };

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
  add(0, 0, 32, Node::encode(100, 0));
  add(1, 0, 24, Node::encode(100, 0));
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
  bool framed, bool corrupt = false, bool invalid_tail_rank = false)
{
  std::ofstream raw_paths, raw_ranks;
  std::unique_ptr<CompressedBlockWriter> framed_paths, framed_ranks;
  if(framed)
  {
    framed_paths.reset(new CompressedBlockWriter(path_name, 4 * KILOBYTE,
      CompressedBlockWriter::ZSTD));
    framed_ranks.reset(new CompressedBlockWriter(rank_name, 4 * KILOBYTE,
      CompressedBlockWriter::ZSTD));
  }
  else
  {
    raw_paths.open(path_name.c_str(), std::ios_base::binary);
    raw_ranks.open(rank_name.c_str(), std::ios_base::binary);
  }

  for(size_type i = 0; i < records.size(); i++)
  {
    PathNode node;
    node.from = records[i].from;
    node.to = Node::encode(10000 + i, 0);
    node.fields = 0; node.setPredecessors(1);
    node.setOrder(1); node.setLCP(1); node.setPointer(2 * i);
    if(records[i].sorted) { node.makeSorted(); }
    if(corrupt && i == 0) { node.setPointer(2 * records.size() + 17); }
    PathNode::rank_type label[2] = {
      (invalid_tail_rank && i + 1 == records.size() ?
        static_cast<PathNode::rank_type>(12345) : records[i].rank), 0 };
    if(framed)
    {
      framed_paths->writeRecord(&node, sizeof(node));
      framed_ranks->writeRecord(label, sizeof(label));
    }
    else
    {
      raw_paths.write(reinterpret_cast<const char*>(&node), sizeof(node));
      raw_ranks.write(reinterpret_cast<const char*>(label), sizeof(label));
    }
  }
  if(framed) { framed_paths->finish(); framed_ranks->finish(); }
}

static void initialize_parallel_prune_graph(PathGraph& graph,
  const std::string& base, PruneFixtureStorage storage, bool corrupt = false,
  bool invalid_tail_rank = false,
  const std::vector<std::vector<PruneFixtureRecord>>* fixture = nullptr)
{
  const std::vector<std::vector<PruneFixtureRecord>> records =
    (fixture != nullptr ? *fixture : parallel_prune_records());
  const logical_file_id_t logical[3] = {
    logical_file_id_t(7), logical_file_id_t(7), logical_file_id_t(8)
  };
  for(size_type file = 0; file < records.size(); file++)
  {
    const bool framed = (storage == FIXTURE_FRAMED ||
      (storage == FIXTURE_MIXED && file == 1));
    std::string path_name = base + "." + std::to_string(file) + ".path";
    std::string rank_name = base + "." + std::to_string(file) + ".rank";
    write_parallel_prune_pair(path_name, rank_name, records[file], framed,
      corrupt && file == 0, invalid_tail_rank && file == 0);
    graph.path_names.push_back(path_name); graph.rank_names.push_back(rank_name);
    graph.path_counts.push_back(records[file].size());
    graph.rank_counts.push_back(2 * records[file].size());
    graph.path_checksums.push_back(ClosedPayloadChecksum());
    graph.rank_checksums.push_back(ClosedPayloadChecksum());
    graph.logical_file_ids.push_back(logical[file]);
    graph.physical_shard_ids.push_back(physical_shard_id_t(100 + file));
    graph.path_count += records[file].size();
    graph.rank_count += 2 * records[file].size();
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

// Shard s holds record i with first label shards * i + s, so a label merge
// alternates shards on every record. With fewer admitted pairs than shards,
// every read evicts and the framed branch re-decodes a whole block per record.
static void write_interleaved_shard(const std::string& path_name,
  const std::string& rank_name, size_type shards, size_type shard,
  size_type records, size_type block_bytes)
{
  std::ofstream raw_paths, raw_ranks;
  std::unique_ptr<CompressedBlockWriter> framed_paths, framed_ranks;
  if(block_bytes == 0)
  {
    raw_paths.open(path_name.c_str(), std::ios_base::binary);
    raw_ranks.open(rank_name.c_str(), std::ios_base::binary);
  }
  else
  {
    framed_paths.reset(new CompressedBlockWriter(path_name, block_bytes,
      CompressedBlockWriter::ZSTD));
    framed_ranks.reset(new CompressedBlockWriter(rank_name, block_bytes,
      CompressedBlockWriter::ZSTD));
  }
  for(size_type i = 0; i < records; i++)
  {
    PathNode node;
    node.from = 1 + shards * i + shard; node.to = i + 1; node.fields = 0;
    node.setPredecessors(1); node.setOrder(1); node.setLCP(1);
    node.setPointer(2 * i);
    PathNode::rank_type label[2] = {
      static_cast<PathNode::rank_type>(shards * i + shard), 9 };
    if(block_bytes == 0)
    {
      raw_paths.write(reinterpret_cast<const char*>(&node), sizeof(node));
      raw_ranks.write(reinterpret_cast<const char*>(label), sizeof(label));
    }
    else
    {
      framed_paths->writeRecord(&node, sizeof(node));
      framed_ranks->writeRecord(label, sizeof(label));
    }
  }
  if(block_bytes != 0) { framed_paths->finish(); framed_ranks->finish(); }
}

// single_logical is the shape compactLogicalJoinShards() leaves behind: many
// physical shards of one named source graph, and therefore one prune() output.
static void build_interleaved_graph(PathGraph& graph,
  const std::vector<std::string>& path_names,
  const std::vector<std::string>& rank_names, size_type records,
  bool single_logical = false)
{
  graph.order = 1;
  graph.logical_file_ids[0] = logical_file_id_t(0);
  graph.physical_shard_ids[0] = physical_shard_id_t(100);
  graph.path_counts[0] = records; graph.rank_counts[0] = 2 * records;
  for(size_type s = 1; s < path_names.size(); s++)
  {
    graph.path_names.push_back(path_names[s]);
    graph.rank_names.push_back(rank_names[s]);
    graph.path_counts.push_back(records);
    graph.rank_counts.push_back(2 * records);
    graph.logical_file_ids.push_back(logical_file_id_t(single_logical ? 0 : s));
    graph.physical_shard_ids.push_back(physical_shard_id_t(100 + s));
  }
  graph.path_count = path_names.size() * records;
  graph.rank_count = 2 * path_names.size() * records;
}

// Framed shards must survive the prune that reads them, and must cost what
// their payload costs.
//
// The block a generation commits is chosen by the producer, but it is admitted
// by this consumer: PathGraphInputCache holds two decoded blocks and two
// decoder contexts per open pair inside pathMergeInputBudget(). Sizing the
// block against the label-sort budget alone committed shards the merger could
// not open at all, and even when one pair did fit, admitting fewer pairs than
// shards made every read evict and re-decode.
static void compare_framed_prune(const std::string& base, const LCP& lcp,
  size_type shards, size_type records)
{
  ConstructionParameters parameters;
  parameters.setMemoryLimitBytes(size_type(2) << 30);
  const size_type budget = pathMergeInputBudget(parameters);
  const size_type pairs = pathMergeInputPairs(parameters);
  const size_type requested = parameters.getCompressionBlockSize();
  const size_type admitted = mergeAdmissibleBlockSize(budget, pairs, requested);

  // The clamp must leave room for every pair the merger may hold open, and
  // must not silently disable compression to get there.
  require(admitted <= requested);
  require(pairs * 2 * CompressedBlockReader::workingMemoryEstimate(admitted) <= budget);

  std::vector<std::string> raw_paths, raw_ranks, framed_paths, framed_ranks;
  std::vector<std::string> legacy_paths, legacy_ranks;
  for(size_type s = 0; s < shards; s++)
  {
    const std::string tag = std::to_string(s);
    raw_paths.push_back(base + ".raw." + tag + ".path");
    raw_ranks.push_back(base + ".raw." + tag + ".rank");
    framed_paths.push_back(base + ".framed." + tag + ".path");
    framed_ranks.push_back(base + ".framed." + tag + ".rank");
    legacy_paths.push_back(base + ".legacy." + tag + ".path");
    legacy_ranks.push_back(base + ".legacy." + tag + ".rank");
    write_interleaved_shard(raw_paths[s], raw_ranks[s], shards, s, records, 0);
    write_interleaved_shard(framed_paths[s], framed_ranks[s], shards, s, records,
      admitted);
    // A workspace committed before the block was bounded keeps its own block,
    // and no operational flag can rewrite it.
    write_interleaved_shard(legacy_paths[s], legacy_ranks[s], shards, s, records,
      requested);
  }

  PathGraph reference(raw_paths[0], raw_ranks[0]);
  PathGraph framed(framed_paths[0], framed_ranks[0]);
  PathGraph legacy(legacy_paths[0], legacy_ranks[0]);
  build_interleaved_graph(reference, raw_paths, raw_ranks, records);
  build_interleaved_graph(framed, framed_paths, framed_ranks, records);
  build_interleaved_graph(legacy, legacy_paths, legacy_ranks, records);

  PathGraphMergeStats raw_stats, framed_stats, legacy_stats;
  reference.prune(lcp, GIGABYTE, budget, &raw_stats, 128);
  framed.prune(lcp, GIGABYTE, budget, &framed_stats, 128);
  // Refusing here would make an already-committed workspace unresumable, so
  // the oversized pair is admitted and reported instead.
  legacy.prune(lcp, GIGABYTE, budget, &legacy_stats, 128);

  require(framed.files() == reference.files());
  require(framed.size() == reference.size());
  require(legacy.size() == reference.size());
  for(size_type file = 0; file < reference.files(); file++)
  {
    require(contents(reference.path_names[file]) == contents(framed.path_names[file]));
    require(contents(reference.rank_names[file]) == contents(framed.rank_names[file]));
    require(contents(reference.path_names[file]) == contents(legacy.path_names[file]));
  }

  // Every shard stays open, so each block is decoded once instead of once per
  // record, and the framed merge is no worse than the raw one it replaced.
  require(framed_stats.max_open_input_pairs == shards);
  require(framed_stats.max_open_input_pairs == raw_stats.max_open_input_pairs);
  require(framed_stats.oversized_input_pair_bytes == 0);
  require(legacy_stats.oversized_input_pair_bytes > budget);

  for(size_type s = 0; s < shards; s++)
  {
    std::remove(raw_paths[s].c_str()); std::remove(raw_ranks[s].c_str());
    std::remove(framed_paths[s].c_str()); std::remove(framed_ranks[s].c_str());
    std::remove(legacy_paths[s].c_str()); std::remove(legacy_ranks[s].c_str());
  }
}

// One logical input over many physical shards, which is what the post-join
// compactor retains and what prune() then has to merge.
//
// prune() emits one output shard per logical input, so an even split of the
// descriptor allowance reserved half of it for output entries nothing occupies
// and left the input cache (max_open_files - 2) / 4 pairs: 31 at the shipped
// 128, for a shard set compaction sizes at up to (max_open_files - 4) / 2.
// The merger visits the shards in round-robin label order, exactly the access
// pattern LRU degrades to a 100% miss rate on, so one entry short of the shard
// count made every read evict. A framed miss then decodes a whole block to
// deliver one 24-byte PathNode, and a raw miss rereads the window.
static void compare_wide_logical_input(const std::string& base, const LCP& lcp,
  size_type shards, size_type records)
{
  ConstructionParameters parameters;
  parameters.setMemoryLimitBytes(size_type(2) << 30);
  const size_type budget = pathMergeInputBudget(parameters);
  const size_type block = mergeAdmissibleBlockSize(budget,
    pathMergeInputPairs(parameters), parameters.getCompressionBlockSize());

  std::vector<std::string> raw_paths, raw_ranks, framed_paths, framed_ranks;
  for(size_type s = 0; s < shards; s++)
  {
    const std::string tag = std::to_string(s);
    raw_paths.push_back(base + ".raw." + tag + ".path");
    raw_ranks.push_back(base + ".raw." + tag + ".rank");
    framed_paths.push_back(base + ".framed." + tag + ".path");
    framed_ranks.push_back(base + ".framed." + tag + ".rank");
    write_interleaved_shard(raw_paths[s], raw_ranks[s], shards, s, records, 0);
    write_interleaved_shard(framed_paths[s], framed_ranks[s], shards, s, records,
      block);
  }

  PathGraph reference(raw_paths[0], raw_ranks[0]);
  PathGraph framed(framed_paths[0], framed_ranks[0]);
  build_interleaved_graph(reference, raw_paths, raw_ranks, records, true);
  build_interleaved_graph(framed, framed_paths, framed_ranks, records, true);
  const size_type logical_bytes = reference.bytes();

  PathGraphMergeStats raw_stats, framed_stats;
  size_type raw_read = DiskIO::read_volume;
  reference.prune(lcp, GIGABYTE, budget, &raw_stats, 128);
  raw_read = DiskIO::read_volume - raw_read;
  size_type framed_read = DiskIO::read_volume;
  framed.prune(lcp, GIGABYTE, budget, &framed_stats, 128);
  framed_read = DiskIO::read_volume - framed_read;

  // One logical input means one output shard, and compression must not change
  // what pruning produces.
  require(reference.files() == 1); require(framed.files() == 1);
  require(framed.size() == reference.size());
  require(contents(reference.path_names[0]) == contents(framed.path_names[0]));
  require(contents(reference.rank_names[0]) == contents(framed.rank_names[0]));

  // Every shard stays open, on one output entry plus two spill descriptors.
  require(raw_stats.max_open_output_pairs == 1);
  require(framed_stats.max_open_output_pairs == 1);
  require(raw_stats.max_open_input_pairs == shards);
  require(framed_stats.max_open_input_pairs == shards);
  require(2 * framed_stats.max_open_input_pairs +
    2 * framed_stats.max_open_output_pairs + 2 <= 128);
  require(framed_stats.oversized_input_pair_bytes == 0);
  require(framed_stats.max_input_buffer_bytes <= budget);

  // The merge therefore costs its own logical bytes, not one block per record.
  require(raw_read <= 2 * logical_bytes);
  require(framed_read <= 2 * logical_bytes);

  for(size_type s = 0; s < shards; s++)
  {
    std::remove(raw_paths[s].c_str()); std::remove(raw_ranks[s].c_str());
    std::remove(framed_paths[s].c_str()); std::remove(framed_ranks[s].c_str());
  }
}

// Exercise the path merger across several blocks per stream. The pooled
// constructor queues block 0 for every resident path/rank reader, while exact
// cursor reads queue each following block without changing merge order.
static void compare_prefetched_multiblock(const std::string& base,
  const LCP& lcp, size_type shards, size_type records)
{
  const size_type block = 4 * KILOBYTE;
  const size_type job =
    CompressedBlockReader::prefetchWorkingMemoryEstimate(block);
  const size_type pair_bytes = 2 *
    CompressedBlockReader::workingMemoryEstimate(block);
  const size_type cache = shards * pair_bytes + 8 * job;

  std::vector<std::string> raw_paths, raw_ranks, framed_paths, framed_ranks;
  size_type expected_blocks = 0;
  for(size_type s = 0; s < shards; ++s)
  {
    const std::string tag = std::to_string(s);
    raw_paths.push_back(base + ".raw." + tag + ".path");
    raw_ranks.push_back(base + ".raw." + tag + ".rank");
    framed_paths.push_back(base + ".framed." + tag + ".path");
    framed_ranks.push_back(base + ".framed." + tag + ".rank");
    write_interleaved_shard(raw_paths[s], raw_ranks[s], shards, s, records, 0);
    write_interleaved_shard(framed_paths[s], framed_ranks[s], shards, s,
      records, block);
    CompressedBlockReader paths(framed_paths[s]);
    CompressedBlockReader ranks(framed_ranks[s]);
    expected_blocks += paths.blocks() + ranks.blocks();
  }

  PathGraph reference(raw_paths[0], raw_ranks[0]);
  PathGraph framed(framed_paths[0], framed_ranks[0]);
  // Keep distinct logical inputs and physical shard identities here. The pool
  // must be observational even when prune projects the global merge into
  // several semantic outputs, not only for the common one-logical-input case.
  build_interleaved_graph(reference, raw_paths, raw_ranks, records);
  build_interleaved_graph(framed, framed_paths, framed_ranks, records);
  const size_type logical_input_bytes = framed.bytes();

  reference.prune(lcp, GIGABYTE, MEGABYTE, nullptr, 128);
  PathGraphMergeStats stats;
  framed.prune(lcp, GIGABYTE, MEGABYTE, &stats, 128, cache);

  require(framed.files() == reference.files());
  require(framed.size() == reference.size());
  for(size_type file = 0; file < reference.files(); ++file)
  {
    require(reference.logicalFile(file) == framed.logicalFile(file));
    require(reference.physicalShard(file) == framed.physicalShard(file));
    require(contents(reference.path_names[file]) == contents(framed.path_names[file]));
    require(contents(reference.rank_names[file]) == contents(framed.rank_names[file]));
  }
  require(stats.max_open_input_pairs == shards);
  require(stats.prefetch_workers > 0 && stats.prefetch_workers <= 8);
  require(stats.prefetch_submitted > 2 * shards);
  require(stats.prefetch_consumed + stats.prefetch_synchronous_blocks ==
    expected_blocks);
  require(stats.prefetch_ready_hits + stats.prefetch_waits ==
    stats.prefetch_consumed);
  require(stats.prefetch_physical_bytes > 0);
  require(stats.prefetch_decoded_bytes == logical_input_bytes);
  require(stats.max_prefetch_bytes <= cache - shards * pair_bytes);
  require(stats.max_input_buffer_bytes <= cache);

  for(size_type s = 0; s < shards; ++s)
  {
    std::remove(raw_paths[s].c_str()); std::remove(raw_ranks[s].c_str());
    std::remove(framed_paths[s].c_str()); std::remove(framed_ranks[s].c_str());
  }
}

// The v0.10 chr21 failure, reduced: shards committed with 16 MiB blocks, the
// default 64 MiB --io-buffer-size, several shards, and a global memory limit
// small enough that the two decode workspaces do not fit the I/O buffer.
//
// The merge workspace must come from the global budget rather than being
// capped by --io-buffer-size, so that one legal unit of progress -- one
// decoded path block plus one decoded rank block -- is always reservable.
// Concurrency is what gives way when the budget is tight, never the build.
static void compare_committed_block_resume(const std::string& base,
  const LCP& lcp, size_type shards, size_type records)
{
  ConstructionParameters parameters;
  // Large enough that a quarter of it funds every shard's decode workspace,
  // which is the regime the chr21 resume runs in (23 GiB against 43 shards),
  // and still small enough that the 64 MiB I/O buffer cannot hold one pair.
  parameters.setMemoryLimitBytes(size_type(8) << 30);
  parameters.setIOBufferSize(64 * MEGABYTE);
  parameters.setCompressionBlockSize(16 * MEGABYTE);
  const size_type committed_block = parameters.getCompressionBlockSize();
  const size_type pair_bytes =
    2 * CompressedBlockReader::workingMemoryEstimate(committed_block);
  const size_type ceiling = pathMergeCeilingBudget(parameters);

  // The shape of the failure: the I/O buffer cannot hold one pair, the global
  // budget can. Nothing below depends on the arithmetic staying exactly here,
  // but the test is pointless if the fixture stops reproducing it.
  require(pair_bytes > pathMergeInputBudget(parameters));
  require(pair_bytes <= ceiling);

  std::vector<std::string> raw_paths, raw_ranks, framed_paths, framed_ranks;
  for(size_type s = 0; s < shards; s++)
  {
    const std::string tag = std::to_string(s);
    raw_paths.push_back(base + ".raw." + tag + ".path");
    raw_ranks.push_back(base + ".raw." + tag + ".rank");
    framed_paths.push_back(base + ".framed." + tag + ".path");
    framed_ranks.push_back(base + ".framed." + tag + ".rank");
    write_interleaved_shard(raw_paths[s], raw_ranks[s], shards, s, records, 0);
    write_interleaved_shard(framed_paths[s], framed_ranks[s], shards, s, records,
      committed_block);
  }

  PathGraph reference(raw_paths[0], raw_ranks[0]);
  PathGraph framed(framed_paths[0], framed_ranks[0]);
  build_interleaved_graph(reference, raw_paths, raw_ranks, records, true);
  build_interleaved_graph(framed, framed_paths, framed_ranks, records, true);

  // The two workspaces the production caller derives for this generation. The
  // group buffer scales with one label range; the input cache scales with the
  // shard count, and must hold every shard: sharing one number between them is
  // what left the chr21 resume holding 22 of its 43 committed shards, which
  // under round-robin label order misses on every record.
  const size_type budget = pathMergeInputBudget(parameters, framed);
  const size_type cache = pathMergeInputCacheBudget(parameters, framed);
  require(budget >= pair_bytes);
  require(budget <= ceiling);
  require(cache >= shards * pair_bytes);
  require(cache <= std::max(pair_bytes, parameters.getMemoryLimitBytes() / 4));

  PathGraphMergeStats stats;
  const size_type one_record = sizeof(PathNode) +
    (PathLabel::LABEL_LENGTH + 1) * sizeof(PathNode::rank_type);
  reference.prune(lcp, GIGABYTE, one_record, nullptr, 128);
  framed.prune(lcp, GIGABYTE, budget, &stats, 128, cache);

  // It finishes, and it finishes with the same index the raw shards produce.
  require(framed.files() == 1);
  require(framed.size() == reference.size());
  require(contents(reference.path_names[0]) == contents(framed.path_names[0]));
  require(contents(reference.rank_names[0]) == contents(framed.rank_names[0]));

  // Every reservation stayed inside the ceiling, and no pair was admitted
  // above the budget it was granted.
  require(stats.oversized_input_pair_bytes == 0);
  require(stats.max_open_input_pairs == shards);
  require(stats.max_open_input_pairs * pair_bytes <= cache);
  require(stats.prefetch_workers > 0 && stats.prefetch_workers <= 8);
  require(stats.prefetch_submitted == 2 * shards);
  require(stats.prefetch_consumed + stats.prefetch_synchronous_blocks ==
    2 * shards);
  require(stats.prefetch_ready_hits + stats.prefetch_waits ==
    stats.prefetch_consumed);
  require(stats.prefetch_physical_bytes > 0);
  require(stats.max_prefetch_bytes <= cache - shards * pair_bytes);
  require(stats.max_input_buffer_bytes <= cache);
  require(2 * stats.max_open_input_pairs + 2 * stats.max_open_output_pairs + 2 <= 128);

  for(size_type s = 0; s < shards; s++)
  {
    std::remove(raw_paths[s].c_str()); std::remove(raw_ranks[s].c_str());
    std::remove(framed_paths[s].c_str()); std::remove(framed_ranks[s].c_str());
  }
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
    // The merger reads each shard monotonically. Bounded 64 KiB windows must
    // therefore turn many logical records into a small number of preads.
    require(reference_stats.path_input_reads == paths.size());
    require(reference_stats.rank_input_reads == paths.size());
    require(reference_stats.path_input_refills > 0);
    require(reference_stats.rank_input_refills > 0);
    require(reference_stats.path_input_refills < reference_stats.path_input_reads);
    require(reference_stats.rank_input_refills < reference_stats.rank_input_reads);
    require(reference_stats.direct_input_reads == 0);
    require(reference_stats.max_input_buffer_bytes <= MEGABYTE / 4);
    // In-memory SameFromSet groups must not invoke the external sorter.
    require(reference_stats.from_set_sorts == 0);
    require(spilled_stats.priority_spills > 0);
    require(spilled_stats.range_spills > 0);
    // A tiny group budget cannot retain even one selected/candidate pair, so
    // it takes the external path while producing byte-identical output above.
    require(spilled_stats.from_set_sorts >= 4);
    require(spilled_stats.max_open_input_pairs <= 1);
    require(spilled_stats.max_input_buffer_bytes == 0);
    require(spilled_stats.direct_input_reads ==
      spilled_stats.path_input_reads + spilled_stats.rank_input_reads);
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

static void write_parallel_merged_fixture(const std::string& name)
{
  Alphabet alpha;
  const std::vector<std::pair<std::string, size_type>> groups = {
    { "AAA", 5000 }, { "AAC", 5000 },
    { "CAA", 300 },  { "CAC", 300 },
    { "GAA", 17 },   { "GAC", 17 },
    { "TAA", 3 },    { "TAC", 3 }
  };
  std::vector<KMer> kmers;
  for(const auto& group : groups)
  {
    comp_type comp = alpha.char2comp[group.first.front()];
    byte_type edge = static_cast<byte_type>(1 << comp);
    key_type key = Key::encode(alpha, group.first, edge, edge);
    size_type component_base = static_cast<size_type>(comp) * 100000;
    for(size_type i = 0; i < group.second; i++)
    {
      // The two labels in each root component have exactly the same start set.
      // extendRange() must merge across that nonzero-LCP label boundary, while
      // the zero-LCP A/C/G/T boundaries remain independent worker tasks.
      kmers.push_back(KMer(key, Node::encode(component_base + i + 1, 0),
        Node::encode(component_base + i + 50001, 0)));
    }
  }
  std::ofstream output(name.c_str(), std::ios_base::binary);
  writeBinary(output, kmers, 3); output.close(); require(static_cast<bool>(output));
}

static void require_same_merged_graph(const MergedGraph& expected,
  const MergedGraph& actual)
{
  require(expected.size() == actual.size());
  require(expected.ranks() == actual.ranks());
  require(expected.extra() == actual.extra());
  require(expected.next == actual.next);
  require(expected.next_from == actual.next_from);
  require(contents(expected.path_name) == contents(actual.path_name));
  require(contents(expected.rank_name) == contents(actual.rank_name));
  require(contents(expected.from_name) == contents(actual.from_name));
  require(contents(expected.lcp_name) == contents(actual.lcp_name));
}

static void compare_parallel_merged_graph(const std::string& base)
{
  const std::string input_name = base + ".graph";
  write_parallel_merged_fixture(input_name);
  ConstructionParameters parameters;
  Alphabet alpha;
  int previous_threads = omp_get_max_threads();
  int previous_dynamic = omp_get_dynamic();
  omp_set_dynamic(0); omp_set_num_threads(8);
  {
    InputGraph input(std::vector<std::string>{ input_name }, true, parameters, alpha);
    std::vector<key_type> keys; input.readKeys(keys);
    require(keys.size() == 8);
    DeBruijnGraph mapper(keys, input.k(), input.alpha);
    LCP lcp(keys, input.k());
    sdsl::int_vector<0> distinct(keys.size(), 0,
      bit_length(Key::label(keys.back())));
    for(size_type i = 0; i < keys.size(); i++) { distinct[i] = Key::label(keys[i]); }
    PathGraph paths(input, distinct);

    const size_type parallel_budget = 256 * KILOBYTE;
    PathGraphMergeStats serial_stats, two_stats, four_stats, fd_stats,
      scratch_limited_stats;
    MergedGraph serial(paths, mapper, lcp, GIGABYTE, MEGABYTE,
      &serial_stats, 128, 0, 1);
    MergedGraph two(paths, mapper, lcp, GIGABYTE, parallel_budget,
      &two_stats, 128, 0, 2);
    MergedGraph four(paths, mapper, lcp, GIGABYTE, parallel_budget,
      &four_stats, 128, 0, 4);
    // The 32-descriptor ceiling admits only two complete worker footprints.
    MergedGraph fd_limited(paths, mapper, lcp, GIGABYTE, parallel_budget,
      &fd_stats, 32, 0, 8);
    // Parallel copy would need two logical copies of the merged streams. The
    // output itself fits this limit exactly, so the bounded implementation
    // must discard its local partitions and complete through the serial route.
    MergedGraph scratch_limited(paths, mapper, lcp, serial.bytes(),
      parallel_budget, &scratch_limited_stats, 128, 0, 4);

    require_same_merged_graph(serial, two);
    require_same_merged_graph(serial, four);
    require_same_merged_graph(serial, fd_limited);
    require_same_merged_graph(serial, scratch_limited);
    require(serial.size() == 4); // Each within-component adjacent pair merged.
    require(serial_stats.merge_workers == 1);
    require(two_stats.merge_workers == 2);
    require(four_stats.merge_workers == 4);
    require(fd_stats.merge_workers == 2);
    require(scratch_limited_stats.merge_workers == 1);
    // One requested worker is not a fallback; a disk-limited request is, and
    // the stats name why.
    require(serial_stats.merge_fallback_reason == nullptr);
    require(scratch_limited_stats.merge_fallback_reason != nullptr &&
      std::string(scratch_limited_stats.merge_fallback_reason).find("disk") != std::string::npos);
    // Seven alphabet components are scheduled, including the empty $, N, and
    // # intervals; A is deliberately much larger than C/G/T.
    require(four_stats.merge_partitions == mapper.alpha.sigma);
    require(four_stats.merge_partitions > four.size());
    require(four_stats.priority_spills > 0);
    require(four_stats.from_set_sorts > 0);
    require(four_stats.prefetch_workers == 0);
    require(2 * four_stats.max_open_input_pairs +
      14 * four_stats.merge_workers <= 128);
    require(four_stats.max_input_buffer_bytes <= parallel_budget);

    // Fifty-seven interleaved framed shards reproduce the chr21 compactor's
    // shape. Dividing 128 FDs among four workers must not reduce the decoded
    // cache to a few pairs and decompress one whole block for every record.
    {
      constexpr size_type shard_count = 57;
      constexpr size_type block_bytes = 64 * KILOBYTE;
      PathGraph wide(shard_count, paths.order, paths.doubling_steps);
      std::vector<PathNode> nodes;
      std::vector<PathNode::rank_type> labels;
      paths.read(nodes, labels, 0);
      size_type stored_bytes = 0;
      for(size_type shard = 0; shard < shard_count; shard++)
      {
        CompressedBlockWriter path_output(wide.path_names[shard], block_bytes,
          CompressedBlockWriter::ZSTD);
        CompressedBlockWriter rank_output(wide.rank_names[shard], block_bytes,
          CompressedBlockWriter::ZSTD);
        wide.logical_file_ids[shard] = logical_file_id_t(0);
        for(size_type i = shard; i < nodes.size(); i += shard_count)
        {
          PathNode node = nodes[i];
          size_type pointer = node.pointer();
          node.setPointer(wide.rank_counts[shard]);
          path_output.writeRecord(&node, sizeof(node));
          rank_output.writeRecord(labels.data() + pointer,
            node.ranks() * sizeof(PathNode::rank_type));
          wide.path_counts[shard]++;
          wide.rank_counts[shard] += node.ranks();
        }
        path_output.finish(); rank_output.finish();
        wide.path_count += wide.path_counts[shard];
        wide.rank_count += wide.rank_counts[shard];
        stored_bytes += CompressedBlockReader(wide.path_names[shard]).physicalSize();
        stored_bytes += CompressedBlockReader(wide.rank_names[shard]).physicalSize();
      }
      PathGraphMergeStats wide_stats;
      size_type before = DiskIO::read_volume;
      ScopedFileLimit limit(128);
      MergedGraph merged(wide, mapper, lcp, GIGABYTE, MEGABYTE,
        &wide_stats, 128, 256 * MEGABYTE, 4);
      size_type read_bytes = DiskIO::read_volume - before;
      require_same_merged_graph(serial, merged);
      require(wide_stats.merge_workers == 4);
      require(wide_stats.max_open_input_pairs <= 4);
      require(wide_stats.max_input_buffer_bytes <= 256 * MEGABYTE);
      // Allow all workers' boundary probes plus ordered-output assembly. This
      // bound is far below repeatedly decoding a shard per merged record.
      require(read_bytes < 32 * stored_bytes + 8 * serial.bytes());
      // A mixed old/new workspace cannot use the all-framed descriptor mode.
      // Its 57 persistent pairs only fit one worker under the same FD cap.
      {
        std::ofstream path_output(wide.path_names[0], std::ios::binary);
        std::ofstream rank_output(wide.rank_names[0], std::ios::binary);
        size_type rank_offset = 0;
        for(size_type i = 0; i < nodes.size(); i += shard_count)
        {
          PathNode node = nodes[i];
          size_type pointer = node.pointer(); node.setPointer(rank_offset);
          path_output.write(reinterpret_cast<const char*>(&node), sizeof(node));
          rank_output.write(reinterpret_cast<const char*>(labels.data() + pointer),
            node.ranks() * sizeof(PathNode::rank_type));
          rank_offset += node.ranks();
        }
      }
      PathGraphMergeStats mixed_stats;
      MergedGraph mixed(wide, mapper, lcp, GIGABYTE, MEGABYTE,
        &mixed_stats, 128, 256 * MEGABYTE, 4);
      require_same_merged_graph(serial, mixed);
      require(mixed_stats.merge_workers == 1);
      require(mixed_stats.merge_fallback_reason != nullptr);
    }

    // A completely empty frontier still takes the parallel assembly route,
    // creates four valid empty streams, and preserves the serial next tables.
    PathGraph empty_paths(2, input.k(), 0);
    PathGraphMergeStats empty_serial_stats, empty_parallel_stats;
    MergedGraph empty_serial(empty_paths, mapper, lcp, GIGABYTE,
      parallel_budget, &empty_serial_stats, 128, 0, 1);
    MergedGraph empty_parallel(empty_paths, mapper, lcp, GIGABYTE,
      parallel_budget, &empty_parallel_stats, 128, 0, 4);
    require_same_merged_graph(empty_serial, empty_parallel);
    require(empty_parallel.size() == 0 && empty_parallel.ranks() == 0 &&
      empty_parallel.extra() == 0);
    require(empty_parallel_stats.merge_workers == 4);
    require(empty_parallel_stats.merge_partitions == mapper.alpha.sigma);
  }
  omp_set_num_threads(previous_threads); omp_set_dynamic(previous_dynamic);
  std::remove(input_name.c_str());
}

static std::vector<std::vector<PruneFixtureRecord>> distinct_start_prune_records()
{
  std::vector<std::vector<PruneFixtureRecord>> result(3);
  for(PathNode::rank_type rank = 0; rank < 8; rank++)
  {
    for(size_type i = 0; i < 2; i++)
    {
      result[0].push_back({ rank, Node::encode(1000 + 2 * rank, 0), false });
      result[1].push_back({ rank, Node::encode(1001 + 2 * rank, 0), false });
    }
  }
  return result;
}

static const char* prune_storage_name(PruneFixtureStorage storage)
{
  if(storage == FIXTURE_RAW) { return "raw"; }
  if(storage == FIXTURE_FRAMED) { return "framed"; }
  return "mixed";
}

static void compare_parallel_prune(const std::string& root,
  PruneFixtureStorage storage)
{
  Alphabet alpha;
  std::vector<key_type> keys = parallel_prune_keys(alpha);
  LCP lcp(keys, 3);
  DeBruijnGraph mapper(keys, 3, alpha);
  const std::string tag = prune_storage_name(storage);
  const size_type group_budget = 4 * KILOBYTE;
  const size_type cache_budget = 64 * MEGABYTE;

  int previous_threads = omp_get_max_threads();
  int previous_dynamic = omp_get_dynamic();
  omp_set_dynamic(0); omp_set_num_threads(std::max(4, previous_threads));
  {
    PathGraph serial(0, 3, 0), parallel(0, 3, 0);
    initialize_parallel_prune_graph(serial,
      root + "/parallel-prune-" + tag + "-serial", storage);
    initialize_parallel_prune_graph(parallel,
      root + "/parallel-prune-" + tag + "-parallel", storage);

    PathGraphMergeStats serial_stats, parallel_stats;
    serial.prune(lcp, GIGABYTE, group_budget, &serial_stats, 128,
      cache_budget, 1);
    parallel.prune(lcp, GIGABYTE, group_budget, &parallel_stats, 128,
      cache_budget, 4);

    require(logical_prune_records(serial) == logical_prune_records(parallel));
    require(serial.size() == parallel.size());
    require(serial.ranks() == parallel.ranks());
    require(serial.ranges() == parallel.ranges());
    require(serial.unique == parallel.unique);
    require(serial.redundant == parallel.redundant);
    require(serial.unsorted == parallel.unsorted);
    require(serial.nondeterministic == parallel.nondeterministic);
    require(parallel.files() > serial.files());
    require(parallel_stats.prune_requested_workers == 4);
    require(parallel_stats.prune_workers == 4);
    require(parallel_stats.prune_partitions == 4);
    require(parallel_stats.prune_parallel_fallbacks == 0);
    require(parallel_stats.priority_spills > 0);
    require(parallel_stats.prefetch_workers == 0);
    require(parallel_stats.max_input_buffer_bytes <= cache_budget);
    require(2 * parallel_stats.max_open_input_pairs +
      2 * parallel_stats.max_open_output_pairs +
      2 * parallel_stats.prune_workers <= 128);
    for(size_type file = 0; file < serial.files(); file++)
    {
      require(serial.path_checksums[file].matches(serial.path_names[file],
        serial.path_counts[file] * sizeof(PathNode)));
      require(serial.rank_checksums[file].matches(serial.rank_names[file],
        serial.rank_counts[file] * sizeof(PathNode::rank_type)));
    }
    for(size_type file = 0; file < parallel.files(); file++)
    {
      require(parallel.physicalShard(file) == physical_shard_id_t(file));
      require(parallel.path_checksums[file].matches(parallel.path_names[file],
        parallel.path_counts[file] * sizeof(PathNode)));
      require(parallel.rank_checksums[file].matches(parallel.rank_names[file],
        parallel.rank_counts[file] * sizeof(PathNode::rank_type)));
    }

    // Retaining root-partition shards must be transparent to the final merge,
    // including its pointer and LCP streams, not merely to another prune.
    PathGraphMergeStats final_serial_stats, final_parallel_stats;
    MergedGraph serial_final(serial, mapper, lcp, GIGABYTE, MEGABYTE,
      &final_serial_stats, 128, cache_budget, 1);
    MergedGraph parallel_final(parallel, mapper, lcp, GIGABYTE, MEGABYTE,
      &final_parallel_stats, 128, cache_budget, 1);
    require_same_merged_graph(serial_final, parallel_final);

    if(storage == FIXTURE_RAW)
    {
      // Six descriptors can make serial progress but cannot retain one full
      // raw three-shard input set in two workers. The request must fall back
      // before opening worker outputs and retain the exact serial result.
      PathGraph limited(0, 3, 0);
      initialize_parallel_prune_graph(limited,
        root + "/parallel-prune-low-fd", storage);
      PathGraphMergeStats limited_stats;
      limited.prune(lcp, GIGABYTE, group_budget, &limited_stats, 6,
        cache_budget, 4);
      require(logical_prune_records(serial) == logical_prune_records(limited));
      require(limited_stats.prune_requested_workers == 4);
      require(limited_stats.prune_workers == 1);
      require(limited_stats.prune_partitions == 1);
      require(limited_stats.prune_parallel_fallbacks == 1);

      // A separate budget for the workers admits them while the six-descriptor
      // ceiling keeps sizing the merge: each worker holds the three raw input
      // pairs, two output pairs and two spill descriptors, twelve in all.
      PathGraph budgeted(0, 3, 0);
      initialize_parallel_prune_graph(budgeted,
        root + "/parallel-prune-worker-budget", storage);
      PathGraphMergeStats budgeted_stats;
      budgeted.prune(lcp, GIGABYTE, group_budget, &budgeted_stats, 6,
        cache_budget, 4, 64);
      require(logical_prune_records(serial) == logical_prune_records(budgeted));
      require(budgeted.unique == serial.unique);
      require(budgeted.unsorted == serial.unsorted);
      require(budgeted_stats.prune_workers == 4);
      require(budgeted_stats.prune_parallel_fallbacks == 0);
      require(2 * budgeted_stats.max_open_input_pairs +
        2 * budgeted_stats.max_open_output_pairs +
        2 * budgeted_stats.prune_workers <= 64);
    }
  }
  omp_set_num_threads(previous_threads); omp_set_dynamic(previous_dynamic);
}

static void compare_parallel_prune_failure_cleanup(const std::string& root)
{
  Alphabet alpha;
  std::vector<key_type> keys = parallel_prune_keys(alpha);
  LCP lcp(keys, 3);
  int previous_threads = omp_get_max_threads();
  int previous_dynamic = omp_get_dynamic();
  omp_set_dynamic(0); omp_set_num_threads(std::max(4, previous_threads));
  {
    auto require_clean_failure = [&](const std::string& name,
      bool corrupt_pointer, bool invalid_tail_rank)
    {
      PathGraph broken(0, 3, 0);
      initialize_parallel_prune_graph(broken,
        root + "/parallel-prune-" + name, FIXTURE_RAW,
        corrupt_pointer, invalid_tail_rank);
      const size_type original_size = broken.size();
      const std::vector<std::string> original_paths = broken.path_names;
      const size_type before = directory_entries(root);
      bool failed = false;
      try
      {
        PathGraphMergeStats stats;
        broken.prune(lcp, GIGABYTE, 4 * KILOBYTE, &stats, 128,
          64 * MEGABYTE, 4);
      }
      catch(const std::runtime_error&) { failed = true; }
      require(failed);
      require(broken.size() == original_size);
      require(broken.path_names == original_paths);
      for(const std::string& path : original_paths)
      {
        require(::access(path.c_str(), F_OK) == 0);
      }
      require(directory_entries(root) == before);
    };
    require_clean_failure("broken-pointer", true, false);
    require_clean_failure("invalid-tail-rank", false, true);
  }
  omp_set_num_threads(previous_threads); omp_set_dynamic(previous_dynamic);
}

// Splitting at three leading characters gives the fixture's eight keys a
// partition each. Ranks 0/1 share a start node and merge across AAA|AAC, so the
// parallel attempt must notice and fall back to the exact serial result; with
// two start nodes per key nothing merges across keys, and eight partitions run
// and publish the serial records.
// Partitions start with a fresh merge heap, so records whose labels are equal
// can leave in a different order than the serial pass emits them; compare those
// as multisets.
static std::map<logical_file_id_t, std::vector<LogicalPruneRecord>>
label_ordered_records(const PathGraph& graph)
{
  std::map<logical_file_id_t, std::vector<LogicalPruneRecord>> result = logical_prune_records(graph);
  for(auto& entry : result)
  {
    std::stable_sort(entry.second.begin(), entry.second.end(),
      [](const LogicalPruneRecord& left, const LogicalPruneRecord& right)
      {
        if(left.label != right.label) { return left.label < right.label; }
        if(left.from != right.from) { return left.from < right.from; }
        if(left.to != right.to) { return left.to < right.to; }
        return left.fields < right.fields;
      });
  }
  return result;
}

static void compare_split_depth_prune(const std::string& root)
{
  Alphabet alpha;
  std::vector<key_type> keys = parallel_prune_keys(alpha);
  LCP lcp(keys, 3);
  const size_type group_budget = 4 * KILOBYTE;
  const size_type cache_budget = 64 * MEGABYTE;
  int previous_threads = omp_get_max_threads();
  int previous_dynamic = omp_get_dynamic();
  omp_set_dynamic(0); omp_set_num_threads(std::max(4, previous_threads));
  require(::setenv("GCSA_EXPERIMENTAL_PRUNE_SPLIT_DEPTH", "3", 1) == 0);
  {
    PathGraph serial(0, 3, 0), parallel(0, 3, 0);
    initialize_parallel_prune_graph(serial, root + "/split-depth-merging-serial", FIXTURE_RAW);
    initialize_parallel_prune_graph(parallel, root + "/split-depth-merging-parallel", FIXTURE_RAW);
    PathGraphMergeStats serial_stats, parallel_stats;
    serial.prune(lcp, GIGABYTE, group_budget, &serial_stats, 128, cache_budget, 1);
    parallel.prune(lcp, GIGABYTE, group_budget, &parallel_stats, 128, cache_budget, 4);
    require(logical_prune_records(serial) == logical_prune_records(parallel));
    require(parallel.unique == serial.unique && parallel.unsorted == serial.unsorted);
    require(parallel_stats.prune_parallel_fallbacks == 1);
    require(parallel_stats.prune_fallback_reason != nullptr &&
      std::string(parallel_stats.prune_fallback_reason).find("prefix split") != std::string::npos);
  }
  {
    const std::vector<std::vector<PruneFixtureRecord>> records = distinct_start_prune_records();
    PathGraph serial(0, 3, 0), parallel(0, 3, 0);
    initialize_parallel_prune_graph(serial, root + "/split-depth-distinct-serial", FIXTURE_RAW,
      false, false, &records);
    initialize_parallel_prune_graph(parallel, root + "/split-depth-distinct-parallel", FIXTURE_RAW,
      false, false, &records);
    PathGraphMergeStats serial_stats, parallel_stats;
    serial.prune(lcp, GIGABYTE, group_budget, &serial_stats, 128, cache_budget, 1);
    parallel.prune(lcp, GIGABYTE, group_budget, &parallel_stats, 128, cache_budget, 4);
    require(label_ordered_records(serial) == label_ordered_records(parallel));
    require(parallel.unique == serial.unique && parallel.unsorted == serial.unsorted);
    require(parallel_stats.prune_parallel_fallbacks == 0);
    require(parallel_stats.prune_partitions == 8);
    require(parallel_stats.prune_workers == 4);
  }
  require(::unsetenv("GCSA_EXPERIMENTAL_PRUNE_SPLIT_DEPTH") == 0);
  omp_set_num_threads(previous_threads); omp_set_dynamic(previous_dynamic);
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

  // A generated generation is framed by default, and the next prune has to be
  // able to open it.
  std::vector<key_type> framed_keys(3 * 200 + 2);
  for(size_type i = 0; i < framed_keys.size(); i++) { framed_keys[i] = i; }
  LCP framed_lcp(framed_keys, 1);
  compare_framed_prune(std::string(root) + "/framed", framed_lcp, 3, 200);

  // Compaction leaves one logical input spread over many shards, and the merge
  // must hold all of them open rather than a quarter of the allowance.
  std::vector<key_type> wide_keys(40 * 200 + 2);
  for(size_type i = 0; i < wide_keys.size(); i++) { wide_keys[i] = i; }
  LCP wide_lcp(wide_keys, 1);
  {
    ScopedFileLimit descriptor_limit(128);
    compare_wide_logical_input(std::string(root) + "/wide", wide_lcp, 40, 200);
  }

  // The shared pool must also hand off later blocks, not just overlap the
  // initial shard heads. Make the concurrency requirement explicit instead
  // of inheriting a potentially single-threaded test environment.
  int previous_threads = omp_get_max_threads();
  omp_set_num_threads(std::max(9, previous_threads));
  std::vector<key_type> prefetch_keys(6 * 1000 + 2);
  for(size_type i = 0; i < prefetch_keys.size(); i++) { prefetch_keys[i] = i; }
  LCP prefetch_lcp(prefetch_keys, 1);
  compare_prefetched_multiblock(std::string(root) + "/prefetch",
    prefetch_lcp, 6, 1000);

  // The configuration that aborted the chr21 v0.10 run: 16 MiB committed
  // blocks and a 64 MiB I/O buffer under a global limit that can afford them.
  std::vector<key_type> resume_keys(24 * 60 + 2);
  for(size_type i = 0; i < resume_keys.size(); i++) { resume_keys[i] = i; }
  LCP resume_lcp(resume_keys, 1);
  compare_committed_block_resume(std::string(root) + "/committed", resume_lcp,
    24, 60);
  omp_set_num_threads(previous_threads);

  // Use the real construction route so mapper ranks and LCP invariants are
  // validated rather than synthesized by the test.
  compare_merged_graph(std::string(root) + "/merged");
  // Force the bounded parallel route at several worker counts. The fixture has
  // empty and highly skewed root components, spills its dominant group, and
  // contains adjacent distinct-label groups that must not be cut apart.
  compare_parallel_merged_graph(std::string(root) + "/parallel-merged");
  // Repeated-round pruning uses exact zero-LCP root partitions. Exercise raw,
  // framed, and mixed input generations, retained physical output shards,
  // logical-file semantics, forced spills, bounded-resource fallback, final
  // merge equality, checksum provenance, and exception cleanup.
  compare_parallel_prune(root, FIXTURE_RAW);
  compare_parallel_prune(root, FIXTURE_FRAMED);
  compare_parallel_prune(root, FIXTURE_MIXED);
  compare_split_depth_prune(root);
  compare_parallel_prune_failure_cleanup(root);
  rmdir(root.c_str());
  return 0;
}
