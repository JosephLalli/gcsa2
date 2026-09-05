#include <gcsa/path_graph.h>
#include <gcsa/compressed_block.h>
#include <gcsa/support.h>

#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <limits>
#include <memory>
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
  compare_wide_logical_input(std::string(root) + "/wide", wide_lcp, 40, 200);

  // The configuration that aborted the chr21 v0.10 run: 16 MiB committed
  // blocks and a 64 MiB I/O buffer under a global limit that can afford them.
  std::vector<key_type> resume_keys(24 * 60 + 2);
  for(size_type i = 0; i < resume_keys.size(); i++) { resume_keys[i] = i; }
  LCP resume_lcp(resume_keys, 1);
  compare_committed_block_resume(std::string(root) + "/committed", resume_lcp,
    24, 60);

  // Use the real construction route so mapper ranks and LCP invariants are
  // validated rather than synthesized by the test.
  compare_merged_graph(std::string(root) + "/merged");
  rmdir(root);
  return 0;
}
