#include <gcsa/final_events.h>
#include <gcsa/internal.h>

#include <array>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>

using namespace gcsa;

namespace
{

void require(bool value) { if(!value) { std::abort(); } }

void writeByte(const std::string& filename, std::streamoff offset,
  std::uint8_t value)
{
  std::fstream file(filename.c_str(), std::ios_base::binary |
    std::ios_base::in | std::ios_base::out);
  require(static_cast<bool>(file));
  file.seekp(offset); file.put(static_cast<char>(value)); file.flush();
  require(static_cast<bool>(file));
}

void writeInteger(const std::string& filename, std::streamoff offset,
  std::uint64_t value)
{
  std::array<std::uint8_t, 8> encoded = {};
  for(size_type i = 0; i < encoded.size(); i++)
  {
    encoded[i] = static_cast<std::uint8_t>(value >> (8 * i));
  }
  std::fstream file(filename.c_str(), std::ios_base::binary |
    std::ios_base::in | std::ios_base::out);
  require(static_cast<bool>(file));
  file.seekp(offset);
  file.write(reinterpret_cast<const char*>(encoded.data()), encoded.size());
  file.flush(); require(static_cast<bool>(file));
}

void initSupports(GCSA& index)
{
  for(size_type comp = 0; comp < index.alpha.sigma; comp++)
  {
    sdsl::util::init_support(index.fast_rank[comp], &(index.fast_bwt[comp]));
    sdsl::util::init_support(index.sparse_rank[comp], &(index.sparse_bwt[comp]));
  }
  sdsl::util::init_support(index.edge_rank, &(index.edges));
  sdsl::util::init_support(index.sampled_path_rank, &(index.sampled_paths));
  sdsl::util::init_support(index.sample_select, &(index.samples));
}

BuildWorkspace::ArtifactRef checkpointWrapped(BuildWorkspace& workspace,
  const ArtifactIdentity& identity, physical_shard_id_t shard,
  const std::string& source, size_type records, const std::string& order)
{
  std::ifstream input(source.c_str(), std::ios_base::binary);
  require(static_cast<bool>(input));
  BuildWorkspace::ArtifactWriter writer = workspace.open_artifact(identity,
    logical_file_id_t(0), shard, order, "all");
  std::array<char, 64> buffer = {};
  while(input)
  {
    input.read(buffer.data(), buffer.size());
    if(input.gcount() > 0)
    {
      writer.write(buffer.data(), static_cast<size_type>(input.gcount()));
    }
  }
  require(input.eof());
  return writer.finish(records);
}

void checkpointLegacyWrappedFinalEvents(BuildWorkspace& workspace,
  const FinalEventFiles& files, const FinalEventMetadata& metadata)
{
  std::vector<BuildWorkspace::ArtifactRef> artifacts;
  artifacts.push_back(checkpointWrapped(workspace,
    ArtifactIdentity("final", "events", "metadata", "final-events-meta-v1"),
    physical_shard_id_t(0), files.metadata, 1, "metadata"));
  artifacts.push_back(checkpointWrapped(workspace,
    ArtifactIdentity("final", "events", "bwt-masks", "bwt-mask-u8-v1"),
    physical_shard_id_t(1), files.bwt_masks, metadata.paths, "path-rank"));
  for(size_type comp = 0; comp < metadata.sigma; comp++)
  {
    artifacts.push_back(checkpointWrapped(workspace,
      ArtifactIdentity("final", "events", "edge-destinations-" +
        std::to_string(comp), "path-rank-u64le-v1"),
      physical_shard_id_t(100 + comp), files.edge_destinations[comp],
      metadata.bwt_counts[comp], "source-path-rank"));
  }
  artifacts.push_back(checkpointWrapped(workspace,
    ArtifactIdentity("final", "events", "sampled-paths", "path-rank-u64le-v1"),
    physical_shard_id_t(2), files.sample_positions, metadata.sampled_paths,
    "path-rank"));
  artifacts.push_back(checkpointWrapped(workspace,
    ArtifactIdentity("final", "events", "sample-ids", "node-id-u64le-v1"),
    physical_shard_id_t(3), files.sample_ids, metadata.sample_ids,
    "sample-order"));
  artifacts.push_back(checkpointWrapped(workspace,
    ArtifactIdentity("final", "events", "sample-ends", "sample-rank-u64le-v1"),
    physical_shard_id_t(4), files.sample_ends, metadata.sampled_paths,
    "sample-order"));
  artifacts.push_back(checkpointWrapped(workspace,
    ArtifactIdentity("final", "events", "occurrences", "path-value-u64le-v1"),
    physical_shard_id_t(5), files.occurrences, metadata.occurrence_items,
    "path-rank"));
  artifacts.push_back(checkpointWrapped(workspace,
    ArtifactIdentity("final", "events", "redundancy", "path-rank-u64le-v1"),
    physical_shard_id_t(6), files.redundant, metadata.redundant,
    "path-rank"));
  workspace.commit_task("final", "events", artifacts);
}

GCSA legacyEquivalent(const Alphabet& alphabet)
{
  GCSA expected;
  expected.header.path_nodes = 4; expected.header.edges = 6; expected.header.order = 8;
  sdsl::int_vector<64> counts(alphabet.sigma, 0);
  counts[0] = 2; counts[1] = 2; counts[2] = 1; counts[5] = 1;
  expected.alpha = Alphabet(counts, alphabet.char2comp, alphabet.comp2char);
  expected.fast_bwt.resize(alphabet.sigma); expected.fast_rank.resize(alphabet.sigma);
  expected.sparse_bwt.resize(alphabet.sigma); expected.sparse_rank.resize(alphabet.sigma);
  std::vector<GCSA::bit_vector> raw(alphabet.sigma, GCSA::bit_vector(4, 0));
  raw[0][0] = raw[1][0] = raw[1][1] = 1;
  raw[2][2] = raw[5][2] = raw[0][3] = 1;
  expected.sparse_bwt[0] = raw[0];
  for(size_type comp = 1; comp <= alphabet.fast_chars; comp++)
  {
    expected.fast_bwt[comp] = raw[comp];
  }
  for(size_type comp = alphabet.fast_chars + 1; comp < alphabet.sigma; comp++)
  {
    expected.sparse_bwt[comp] = raw[comp];
  }

  GCSA::bit_vector edges(6, 0);
  edges[1] = edges[2] = edges[4] = edges[5] = 1;
  expected.edges = edges;
  GCSA::bit_vector sampled(4, 0); sampled[0] = sampled[3] = 1;
  expected.sampled_paths = sampled;
  expected.stored_samples = sdsl::int_vector<0>(3, 0, bit_length(100));
  expected.stored_samples[0] = 5; expected.stored_samples[1] = 7;
  expected.stored_samples[2] = 100;
  expected.samples = GCSA::bit_vector(3, 0);
  expected.samples[1] = expected.samples[2] = 1;

  CounterArray occurrences(4, 4), redundant(3, 4);
  occurrences.increment(0, 2); occurrences.increment(2, 20);
  redundant.increment(0); redundant.increment(2, 2);
  expected.extra_pointers = SadaSparse(occurrences);
  expected.redundant_pointers = SadaCount(redundant);
  initSupports(expected);
  return expected;
}

} // namespace

int main()
{
  char root[] = "/tmp/gcsa-final-events-XXXXXX";
  require(::mkdtemp(root) != nullptr);
  TempFile::setDirectory(root);

  ConstructionParameters parameters;
  parameters.setMemoryLimitBytes(64 * KILOBYTE);
  parameters.setSortRunSize(64 * KILOBYTE);
  parameters.setIOBufferSize(128);
  parameters.setMergeFanIn(2);

  // A path may represent more mapped start nodes than fit in RAM. Exercise
  // both reusable in-memory sets and the forced external sort/dedup path with
  // the minimum legal byte reservation.
  const size_type node_set_budget = SpillableNodeSet::minimumBudget();
  MemoryBudget node_memory(2 * node_set_budget);
  {
    SpillableNodeSet current(node_set_budget, 2, node_memory);
    SpillableNodeSet predecessor(node_set_budget, 2, node_memory);
    for(node_type value = 2048; value > 0; value--)
    {
      current.push_back(value + 100);
      current.push_back(value + 100); // Deduplicate across spilled runs.
      predecessor.push_back(value + 99);
    }
    current.finish(); predecessor.finish();
    require(current.spilled() && predecessor.spilled());
    require(current.size() == 2048 && predecessor.size() == 2048);

    current.rewind(); predecessor.rewind();
    node_type curr = 0, pred = 0;
    for(node_type expected = 101; expected <= 2148; expected++)
    {
      require(current.next(curr) && predecessor.next(pred));
      require(curr == expected && curr == pred + 1);
    }
    require(!current.next(curr) && !predecessor.next(pred));

    // Reuse releases spill artifacts but retains the same admitted workspace.
    current.clear();
    current.push_back(9); current.push_back(3); current.push_back(9);
    current.finish();
    require(!current.spilled() && current.size() == 2);
    current.rewind();
    require(current.next(curr) && curr == 3);
    require(current.next(curr) && curr == 9);
    require(!current.next(curr));

    // A rewind must reproduce the complete sorted stream; the final scan may
    // make three current-set passes before emitting samples.
    current.rewind();
    require(current.next(curr) && curr == 3);
    require(current.next(curr) && curr == 9);
    require(!current.next(curr));

    MemoryBudget::Stats node_stats = node_memory.stats();
    require(node_stats.current == 2 * node_set_budget);
    require(node_stats.maximum == 2 * node_set_budget);
  }
  require(node_memory.stats().current == 0);

  bool tiny_node_budget_rejected = false;
  try
  {
    MemoryBudget tiny_memory(node_set_budget);
    SpillableNodeSet invalid(node_set_budget - 1, 2, tiny_memory);
  }
  catch(const std::runtime_error&) { tiny_node_budget_rejected = true; }
  require(tiny_node_budget_rejected);

  // Components 1..4 use the fast representation in the default alphabet;
  // components 0, 5, and 6 exercise the sparse sides of the split.
  Alphabet alphabet;
  FinalEventFiles files(alphabet.sigma);
  FinalEventMetadata metadata;
  {
    MemoryBudget budget(4096);
    FinalEventWriter writer(files, alphabet.sigma, 32, budget);
    writer.path((1U << 0) | (1U << 1));
    writer.edge(0, 0); writer.edge(1, 0);
    writer.sampledPath(0); writer.sample(5); writer.sample(7); writer.sampleEnd();
    writer.occurrence(0, 2);

    writer.path(1U << 1); writer.edge(1, 1);

    writer.path((1U << 2) | (1U << 5));
    writer.edge(2, 2); writer.edge(5, 3);
    writer.occurrence(2, 20);

    writer.path(1U << 0); writer.edge(0, 2);
    writer.sampledPath(3); writer.sample(100); writer.sampleEnd();

    // Deliberately unsorted with a duplicate; external sorting must produce
    // counts [1, 0, 2] for the three suffix-tree slots.
    writer.redundancy(2); writer.redundancy(0); writer.redundancy(2);
    metadata = writer.finish();
  }
  metadata.fast_chars = alphabet.fast_chars;
  sortFinalRedundancy(files, parameters);
  writeFinalEventMetadata(files, metadata);

  BuildWorkspace::Settings semantic;
  semantic["fixture"] = "final-events";
  BuildWorkspace workspace(root, semantic, BuildWorkspace::Settings(),
    BuildWorkspace::NEW_WORKSPACE);
  checkpointFinalEvents(workspace, files, metadata, 64);
  require(workspace.task_completed("final", "events"));

  // Same-filesystem checkpointing adopts the immutable payload inode. This is
  // the property that removes a complete event-stream rewrite.
  struct stat source_status, artifact_status;
  require(::stat(files.bwt_masks.c_str(), &source_status) == 0);
  const std::string mask_artifact = workspace.artifact_path(
    ArtifactIdentity("final", "events", "bwt-masks", "bwt-mask-u8-v1"),
    logical_file_id_t(0), physical_shard_id_t(1));
  require(::stat(mask_artifact.c_str(), &artifact_status) == 0);
  require(source_status.st_dev == artifact_status.st_dev);
  require(source_status.st_ino == artifact_status.st_ino);
  require(artifact_status.st_nlink >= 2);

  FinalEventFiles restored(alphabet.sigma);
  FinalEventMetadata restored_metadata;
  require(restoreFinalEvents(workspace, restored, restored_metadata, 4,
    alphabet.sigma, 64));
  require(restored_metadata.total_edges == 6);
  require(restored_metadata.occurrence_extra == 22);
  require(restored_metadata.redundant == 3);

  struct stat restored_status;
  require(::stat(restored.bwt_masks.c_str(), &restored_status) == 0);
  require(restored_status.st_dev == artifact_status.st_dev);
  require(restored_status.st_ino == artifact_status.st_ino);

  // Workspaces produced before raw-payload adoption remain resumable. Their
  // generic artifact header/footer makes the committed file larger than the
  // expected raw stream, selecting the checked compatibility restore path.
  const std::string legacy_root = std::string(root) + "/legacy-workspace";
  BuildWorkspace legacy_workspace(legacy_root, semantic,
    BuildWorkspace::Settings(), BuildWorkspace::NEW_WORKSPACE);
  checkpointLegacyWrappedFinalEvents(legacy_workspace, files, metadata);
  FinalEventFiles legacy_restored(alphabet.sigma);
  FinalEventMetadata legacy_metadata;
  require(restoreFinalEvents(legacy_workspace, legacy_restored, legacy_metadata,
    4, alphabet.sigma, 64));
  require(legacy_metadata.total_edges == metadata.total_edges);
  require(legacy_metadata.sample_ids == metadata.sample_ids);
  legacy_restored.clear();

  GCSA observed;
  observed.header.path_nodes = 4; observed.header.edges = 6; observed.header.order = 8;
  buildFinalComponents(observed, alphabet, restored, restored_metadata, parameters);
  require(observed.extra_pointers.count(0, 0) == 2);
  require(observed.extra_pointers.count(1, 1) == 0);
  require(observed.extra_pointers.count(2, 2) == 20);
  require(observed.redundant_pointers.count(0, 0) == 1);
  require(observed.redundant_pointers.count(1, 1) == 0);
  require(observed.redundant_pointers.count(2, 2) == 2);
  require(observed.sampled(0) && observed.sampled(3));
  require(observed.sampleCount() == 3);
  require(observed.sample(0) == 5 && observed.sample(1) == 7 && observed.sample(2) == 100);

  GCSA expected = legacyEquivalent(alphabet);
  std::ostringstream observed_bytes, expected_bytes;
  observed.serialize(observed_bytes); expected.serialize(expected_bytes);
  require(observed_bytes.str() == expected_bytes.str());

  // Direct packing preserves the public serialization exactly while keeping
  // component construction local to the writer. Loading proves that the
  // stream order still matches GCSA::load().
  const std::string packed_name = std::string(root) + "/packed.gcsa";
  storeFinalComponents(observed.header, alphabet, restored, restored_metadata,
    parameters, packed_name);
  std::ifstream packed_input(packed_name.c_str(), std::ios_base::binary);
  std::ostringstream packed_bytes; packed_bytes << packed_input.rdbuf();
  require(packed_bytes.str() == expected_bytes.str());
  std::istringstream packed_stream(packed_bytes.str());
  GCSA packed; packed.load(packed_stream);
  require(packed.extra_pointers.count(2, 2) == 20);
  require(packed.redundant_pointers.count(2, 2) == 2);
  require(packed.sampleCount() == 3 && packed.sample(2) == 100);

  const auto rejects = [&](const FinalEventMetadata& candidate_metadata)
  {
    GCSA candidate;
    try
    {
      buildFinalComponents(candidate, alphabet, restored, candidate_metadata,
        parameters);
    }
    catch(const std::runtime_error&) { return true; }
    return false;
  };

  // Metadata and streams are validated before unchecked succinct-builder
  // assumptions can become out-of-bounds writes.
  FinalEventMetadata invalid_metadata = restored_metadata;
  invalid_metadata.sample_bits = 65;
  require(rejects(invalid_metadata));

  // Path 1 originally has only component 1. Adding a second component-5 bit
  // exceeds the declared sparse-BWT one-count and must fail closed.
  writeByte(restored.bwt_masks, 1, static_cast<std::uint8_t>((1U << 1) | (1U << 5)));
  require(rejects(restored_metadata));
  writeByte(restored.bwt_masks, 1, static_cast<std::uint8_t>(1U << 1));

  // A corrupt occurrence value larger than the declared total previously
  // underflowed unsigned subtraction before set_unsafe().
  writeInteger(restored.occurrences, 8,
    std::numeric_limits<std::uint64_t>::max());
  require(rejects(restored_metadata));
  writeInteger(restored.occurrences, 8, 2);

  invalid_metadata = restored_metadata;
  invalid_metadata.sample_bits = 1;
  require(rejects(invalid_metadata));

  // Raw truncation is rejected before SDSL sees an inconsistent universe.
  require(::truncate(restored.occurrences.c_str(), 16) == 0);
  bool rejected = false;
  try { buildFinalComponents(observed, alphabet, restored, restored_metadata, parameters); }
  catch(const std::runtime_error&) { rejected = true; }
  require(rejected);

  restored.clear(); files.clear();
  std::filesystem::remove_all(root);
  return 0;
}
