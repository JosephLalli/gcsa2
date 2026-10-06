#include <gcsa/final_events.h>
#include <gcsa/internal.h>
#include "final_events_internal.hpp"

#include <array>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

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

void appendInteger(std::ostream& out, std::uint64_t value)
{
  std::array<std::uint8_t, 8> encoded = {};
  for(size_type i = 0; i < encoded.size(); i++)
  {
    encoded[i] = static_cast<std::uint8_t>(value >> (8 * i));
  }
  out.write(reinterpret_cast<const char*>(encoded.data()), encoded.size());
}

std::string readFile(const std::string& filename)
{
  std::ifstream input(filename.c_str(), std::ios_base::binary);
  return std::string(std::istreambuf_iterator<char>(input),
    std::istreambuf_iterator<char>());
}

bool hasComponentSpool(const std::string& directory)
{
  for(const auto& entry : std::filesystem::directory_iterator(directory))
  {
    if(entry.path().filename().string().find("gcsa_final_component") !=
       std::string::npos)
    {
      return true;
    }
  }
  return false;
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

// Build through the path production uses -- storeFinalComponents publishes a
// .gcsa and GCSA::load() reads it -- then hand back the loaded index. The
// resident assembler these tests used to call was a second implementation of
// this against a different sink, reachable from nothing that ships.
GCSA
buildViaStore(const Alphabet& alphabet, const FinalEventFiles& files,
  const FinalEventMetadata& metadata, const ConstructionParameters& parameters,
  size_type path_nodes, size_type edges, size_type order)
{
  GCSAHeader header;
  header.path_nodes = path_nodes; header.edges = edges; header.order = order;
  std::string published = TempFile::getName("gcsa_final_components_test");
  GCSA result;
  try
  {
    storeFinalComponents(header, alphabet, files, metadata, parameters, published);
    if(!sdsl::load_from_file(result, published))
    {
      throw std::runtime_error("buildViaStore(): cannot load the published index");
    }
  }
  catch(...) { TempFile::remove(published); throw; }
  TempFile::remove(published);
  return result;
}

} // namespace

int main()
{
  const char* temp_directory = std::getenv("TMPDIR");
  std::string root_storage = std::string(
    temp_directory != nullptr && temp_directory[0] != '\0' ? temp_directory : "/tmp") +
    "/gcsa-final-events-XXXXXX";
  char* root = root_storage.data();
  require(::mkdtemp(root) != nullptr);
  TempFile::setDirectory(root);

  ConstructionParameters parameters;
  parameters.setMemoryLimitBytes(64 * KILOBYTE);

  // Cross the bit_vector_il threshold where SDSL stores its auxiliary
  // breadth-first rank samples. Exact byte equality here covers both the
  // interleaved block layout and the bounded 1024-sample construction.
  const size_type large_paths = 4 * MEGABYTE;
  const std::string large_masks = std::string(root) + "/large-masks";
  GCSA::bit_vector expected_dense(large_paths, 0);
  size_type expected_ones = 0;
  {
    std::ofstream output(large_masks.c_str(), std::ios_base::binary);
    require(static_cast<bool>(output));
    std::array<std::uint8_t, 4096> block = {};
    for(size_type first = 0; first < large_paths; first += block.size())
    {
      size_type count = std::min<size_type>(block.size(), large_paths - first);
      for(size_type offset = 0; offset < count; offset++)
      {
        size_type path = first + offset;
        bool set = (path % 7 == 0 || path % 1021 == 0);
        block[offset] = (set ? static_cast<std::uint8_t>(1U << 1) : 0);
        if(set) { expected_dense[path] = 1; expected_ones++; }
      }
      output.write(reinterpret_cast<const char*>(block.data()), count);
      require(static_cast<bool>(output));
    }
  }
  std::ostringstream streamed_fast, expected_fast;
  serializeFastBWTComponent(streamed_fast, large_masks, large_paths,
    expected_ones, 1, parameters);
  GCSA::fast_vector expected_interleaved(expected_dense);
  expected_interleaved.serialize(expected_fast);
  require(streamed_fast.str() == expected_fast.str());

  // The same immutable mask stream can be replayed to encode all Elias--Fano
  // members directly. Compare bytes above the SDSL fast-select threshold so
  // the test covers multiple select blocks and both high-bit supports.
  std::ostringstream streamed_sparse_bwt, expected_sparse_bwt;
  serializeSparseBWTComponent(streamed_sparse_bwt, large_masks, large_paths,
    expected_ones, 1, parameters);
  GCSA::sparse_vector expected_sparse(expected_dense);
  expected_sparse.serialize(expected_sparse_bwt);
  require(streamed_sparse_bwt.str() == expected_sparse_bwt.str());
  require(std::filesystem::remove(large_masks));

  // Exercise packed sample values crossing 64-bit word boundaries.
  const std::array<std::uint64_t, 5> packed_values = {
    1, 131071, 17, 65535, 42
  };
  const size_type packed_width = 17;
  const std::string packed_ids = std::string(root) + "/packed-ids";
  { std::ofstream create(packed_ids.c_str(), std::ios_base::binary); }
  sdsl::int_vector<0> expected_ids(packed_values.size(), 0, packed_width);
  for(size_type i = 0; i < packed_values.size(); i++)
  {
    writeInteger(packed_ids, i * 8, packed_values[i]);
    expected_ids[i] = packed_values[i];
  }
  std::ostringstream streamed_ids, expected_id_bytes;
  serializeSampleIds(streamed_ids, packed_ids, packed_values.size(),
    packed_width, parameters);
  expected_ids.serialize(expected_id_bytes);
  require(streamed_ids.str() == expected_id_bytes.str());
  require(std::filesystem::remove(packed_ids));

  // Sample boundaries are a monotone event stream. Direct serialization must
  // match the ordinary bit_vector followed by its select_support_mcl without
  // allocating a bit for every sample ID.
  const size_type boundary_count = 5 * 4096 + 17;
  const size_type boundary_universe = boundary_count * 7 + 11;
  const std::string boundary_file = std::string(root) + "/sample-boundaries";
  GCSA::bit_vector expected_boundaries(boundary_universe, 0);
  {
    std::ofstream output(boundary_file.c_str(), std::ios_base::binary);
    require(static_cast<bool>(output));
    for(size_type i = 0; i < boundary_count; i++)
    {
      const std::uint64_t endpoint = (i + 1 == boundary_count ?
        boundary_universe - 1 : i * 7);
      output.write(reinterpret_cast<const char*>(&endpoint), sizeof(endpoint));
      expected_boundaries[endpoint] = 1;
    }
    require(static_cast<bool>(output));
  }
  GCSA::bit_vector::select_1_type expected_boundary_select;
  sdsl::util::init_support(expected_boundary_select, &expected_boundaries);
  std::ostringstream streamed_boundaries, expected_boundary_bytes;
  serializeSampleBoundaries(streamed_boundaries, boundary_file,
    boundary_universe, boundary_count, parameters);
  expected_boundaries.serialize(expected_boundary_bytes);
  expected_boundary_select.serialize(expected_boundary_bytes);
  require(streamed_boundaries.str() == expected_boundary_bytes.str());
  require(std::filesystem::remove(boundary_file));

  // Exercise a complete long select block as well as the fast initializer's
  // final partial block. The complete block uses its actual last position for
  // packed width, while the partial block uses the bitvector universe.
  const size_type long_boundary_count = 4097;
  const size_type long_boundary_universe = 1000003;
  const std::string long_boundary_file = std::string(root) + "/long-sample-boundaries";
  GCSA::bit_vector expected_long_boundaries(long_boundary_universe, 0);
  {
    std::ofstream output(long_boundary_file.c_str(), std::ios_base::binary);
    require(static_cast<bool>(output));
    for(size_type i = 0; i < long_boundary_count; i++)
    {
      const std::uint64_t endpoint = (i + 1 == long_boundary_count ?
        long_boundary_universe - 1 : i * 200);
      appendInteger(output, endpoint);
      expected_long_boundaries[endpoint] = 1;
    }
    require(static_cast<bool>(output));
  }
  GCSA::bit_vector::select_1_type expected_long_boundary_select;
  sdsl::util::init_support(expected_long_boundary_select,
    &expected_long_boundaries);
  std::ostringstream streamed_long_boundaries, expected_long_boundary_bytes;
  serializeSampleBoundaries(streamed_long_boundaries, long_boundary_file,
    long_boundary_universe, long_boundary_count, parameters);
  expected_long_boundaries.serialize(expected_long_boundary_bytes);
  expected_long_boundary_select.serialize(expected_long_boundary_bytes);
  require(streamed_long_boundaries.str() == expected_long_boundary_bytes.str());
  require(std::filesystem::remove(long_boundary_file));

  // The fast initializer's cross-block lookahead can also change a full block
  // from compact to long. Pin that private SDSL quirk because final GCSA files
  // must remain byte-compatible with the existing loader and serializer.
  const size_type lookahead_boundary_count = 4097;
  const size_type lookahead_boundary_universe = 1000000;
  const std::string lookahead_boundary_file =
    std::string(root) + "/lookahead-sample-boundaries";
  GCSA::bit_vector expected_lookahead_boundaries(
    lookahead_boundary_universe, 0);
  {
    std::ofstream output(lookahead_boundary_file.c_str(), std::ios_base::binary);
    require(static_cast<bool>(output));
    for(size_type i = 0; i < lookahead_boundary_count; i++)
    {
      const std::uint64_t endpoint = (i + 1 == lookahead_boundary_count ?
        lookahead_boundary_universe - 1 : i);
      appendInteger(output, endpoint);
      expected_lookahead_boundaries[endpoint] = 1;
    }
    require(static_cast<bool>(output));
  }
  GCSA::bit_vector::select_1_type expected_lookahead_boundary_select;
  sdsl::util::init_support(expected_lookahead_boundary_select,
    &expected_lookahead_boundaries);
  std::ostringstream streamed_lookahead_boundaries,
    expected_lookahead_boundary_bytes;
  serializeSampleBoundaries(streamed_lookahead_boundaries,
    lookahead_boundary_file, lookahead_boundary_universe,
    lookahead_boundary_count, parameters);
  expected_lookahead_boundaries.serialize(expected_lookahead_boundary_bytes);
  expected_lookahead_boundary_select.serialize(
    expected_lookahead_boundary_bytes);
  require(streamed_lookahead_boundaries.str() ==
    expected_lookahead_boundary_bytes.str());
  require(std::filesystem::remove(lookahead_boundary_file));

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

    // Reuse releases spill artifacts but retains the same memory reservation.
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

  GCSA observed = buildViaStore(alphabet, files, metadata,
    parameters, 4, 6, 8);
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

  // The direct redundancy writer covers duplicate positions, a zero run, and
  // the final suffix-tree slot in the same native SadaCount byte layout.
  std::ostringstream streamed_redundancy, expected_redundancy;
  serializeRedundantPointers(streamed_redundancy, files.redundant, 4, 3,
    parameters);
  expected.redundant_pointers.serialize(expected_redundancy);
  require(streamed_redundancy.str() == expected_redundancy.str());

  // Cross select_support_mcl's 4096-one superblock boundary and its final
  // partial-block path. This is a direct byte comparison against SDSL.
  const size_type large_slots = 5 * 4096 + 17;
  const std::string large_redundancy = std::string(root) + "/large-redundancy";
  std::vector<size_type> large_counts(large_slots, 0);
  size_type large_events = 0;
  {
    std::ofstream output(large_redundancy.c_str(), std::ios_base::binary);
    require(static_cast<bool>(output));
    for(size_type slot = 0; slot < large_slots; slot++)
    {
      large_counts[slot] = (slot % 5 == 0 ? 2 : (slot % 11 == 0 ? 1 : 0));
      for(size_type i = 0; i < large_counts[slot]; i++)
      {
        std::uint64_t value = slot;
        output.write(reinterpret_cast<const char*>(&value), sizeof(value));
        large_events++;
      }
    }
    require(static_cast<bool>(output));
  }
  GCSA::bit_vector large_data(large_slots + large_events, 0);
  size_type large_tail = 0;
  for(size_type slot = 0; slot < large_slots; slot++)
  {
    large_tail += large_counts[slot];
    large_data[slot + large_tail] = 1;
  }
  SadaCount large_native;
  large_native.data = large_data;
  sdsl::util::init_support(large_native.select, &(large_native.data));
  std::ostringstream streamed_large_redundancy, expected_large_redundancy;
  serializeRedundantPointers(streamed_large_redundancy, large_redundancy,
    large_slots + 1, large_events, parameters);
  large_native.serialize(expected_large_redundancy);
  require(streamed_large_redundancy.str() == expected_large_redundancy.str());
  require(std::filesystem::remove(large_redundancy));

  // Both Elias--Fano vectors cross the select_support_mcl fast-construction
  // threshold. Exact equality proves the direct writer preserves SadaSparse's
  // low/high vectors and both high-bit select supports without resident
  // sd_vector_builder instances.
  const size_type occurrence_paths = 1000003;
  const std::string large_occurrence_file = std::string(root) + "/large-occurrences";
  CounterArray large_occurrences(occurrence_paths, 4);
  size_type occurrence_items = 0, occurrence_extra = 0;
  {
    std::ofstream output(large_occurrence_file.c_str(), std::ios_base::binary);
    require(static_cast<bool>(output));
    for(size_type path = 3; path < occurrence_paths; path += 19)
    {
      const std::uint64_t extra = 1 + (occurrence_items % 13);
      const std::uint64_t encoded_path = path;
      appendInteger(output, encoded_path);
      appendInteger(output, extra);
      large_occurrences.increment(path, extra);
      occurrence_items++; occurrence_extra += extra;
    }
    require(static_cast<bool>(output));
  }
  SadaSparse native_occurrences(large_occurrences);
  std::ostringstream streamed_occurrences, expected_occurrences;
  serializeOccurrencePointers(streamed_occurrences, large_occurrence_file,
    occurrence_paths, occurrence_items, occurrence_extra, parameters);
  native_occurrences.serialize(expected_occurrences);
  require(streamed_occurrences.str() == expected_occurrences.str());
  require(std::filesystem::remove(large_occurrence_file));

  // Regression for select_support_mcl::init_fast(): for each nonterminal full
  // 4096-one block, SDSL includes the first position of the next block when it
  // chooses the miniblock width. This dense fixture has exactly 24 blocks and
  // exposed a byte mismatch in the first streaming implementation.
  const size_type dense_occurrence_paths = 196608;
  const std::string dense_occurrence_file =
    std::string(root) + "/dense-occurrences";
  CounterArray dense_occurrences(dense_occurrence_paths, 2);
  size_type dense_occurrence_items = 0;
  {
    std::ofstream output(dense_occurrence_file.c_str(), std::ios_base::binary);
    require(static_cast<bool>(output));
    for(size_type path = 0; path < dense_occurrence_paths; path += 2)
    {
      const std::uint64_t encoded_path = path, extra = 1;
      appendInteger(output, encoded_path);
      appendInteger(output, extra);
      dense_occurrences.increment(path, extra);
      dense_occurrence_items++;
    }
    require(static_cast<bool>(output));
  }
  SadaSparse native_dense_occurrences(dense_occurrences);
  std::ostringstream streamed_dense_occurrences, expected_dense_occurrences;
  serializeOccurrencePointers(streamed_dense_occurrences,
    dense_occurrence_file, dense_occurrence_paths, dense_occurrence_items,
    dense_occurrence_items, parameters);
  native_dense_occurrences.serialize(expected_dense_occurrences);
  require(streamed_dense_occurrences.str() == expected_dense_occurrences.str());
  require(std::filesystem::remove(dense_occurrence_file));

  // The same fast-initializer lookahead can change the representation, not
  // merely its packed width. Keep the first 4096 paths dense and move the next
  // path far enough away that SDSL classifies the first select block as long.
  // Build the native reference with sparse builders so the large path universe
  // does not require a correspondingly large CounterArray fixture.
  const size_type lookahead_occurrence_paths = 835000000;
  const size_type lookahead_occurrence_items = 50000;
  const size_type select_block_size = 4096;
  const size_type lookahead_second_block = lookahead_occurrence_paths -
    (lookahead_occurrence_items - select_block_size);
  const auto lookahead_params = SadaSparse::sd_vector::get_params(
    lookahead_occurrence_paths, lookahead_occurrence_items);
  const size_type lookahead_high_bits = lookahead_params.second;
  const size_type lookahead_logn = sdsl::bits::hi(
    ((lookahead_high_bits + 63) / 64) * 64) + 1;
  const size_type lookahead_logn4 = lookahead_logn * lookahead_logn *
    lookahead_logn * lookahead_logn;
  require(lookahead_high_bits >= 100000);
  require(((select_block_size - 1) >> lookahead_params.first) +
    (select_block_size - 1) <= lookahead_logn4);
  require((lookahead_second_block >> lookahead_params.first) +
    select_block_size > lookahead_logn4);
  const std::string lookahead_occurrence_file =
    std::string(root) + "/lookahead-occurrences";
  sdsl::sd_vector_builder lookahead_filter_builder(
    lookahead_occurrence_paths, lookahead_occurrence_items);
  sdsl::sd_vector_builder lookahead_values_builder(
    lookahead_occurrence_items, lookahead_occurrence_items);
  {
    std::ofstream output(lookahead_occurrence_file.c_str(),
      std::ios_base::binary);
    require(static_cast<bool>(output));
    for(size_type i = 0; i < lookahead_occurrence_items; i++)
    {
      const size_type path = (i < select_block_size ? i :
        lookahead_second_block + (i - select_block_size));
      appendInteger(output, path); appendInteger(output, 1);
      lookahead_filter_builder.set(path); lookahead_values_builder.set(i);
    }
    require(static_cast<bool>(output));
  }
  SadaSparse native_lookahead_occurrences;
  native_lookahead_occurrences.filter =
    SadaSparse::sd_vector(lookahead_filter_builder);
  native_lookahead_occurrences.values =
    SadaSparse::sd_vector(lookahead_values_builder);
  sdsl::util::init_support(native_lookahead_occurrences.filter_rank,
    &(native_lookahead_occurrences.filter));
  sdsl::util::init_support(native_lookahead_occurrences.value_select,
    &(native_lookahead_occurrences.values));
  std::ostringstream streamed_lookahead_occurrences,
    expected_lookahead_occurrences;
  serializeOccurrencePointers(streamed_lookahead_occurrences,
    lookahead_occurrence_file, lookahead_occurrence_paths,
    lookahead_occurrence_items, lookahead_occurrence_items, parameters);
  native_lookahead_occurrences.serialize(expected_lookahead_occurrences);
  require(streamed_lookahead_occurrences.str() ==
    expected_lookahead_occurrences.str());
  require(std::filesystem::remove(lookahead_occurrence_file));

  const std::string empty_occurrence_file = std::string(root) + "/empty-occurrences";
  { std::ofstream create(empty_occurrence_file.c_str(), std::ios_base::binary); }
  // A large empty filter still has a large Elias--Fano high vector and a
  // nonempty select_0 support; keep it above init_fast()'s threshold.
  const size_type empty_occurrence_paths = 200000;
  CounterArray empty_occurrences(empty_occurrence_paths, 2);
  SadaSparse native_empty_occurrences(empty_occurrences);
  std::ostringstream streamed_empty_occurrences, expected_empty_occurrences;
  serializeOccurrencePointers(streamed_empty_occurrences, empty_occurrence_file,
    empty_occurrence_paths, 0, 0, parameters);
  native_empty_occurrences.serialize(expected_empty_occurrences);
  require(streamed_empty_occurrences.str() == expected_empty_occurrences.str());
  require(std::filesystem::remove(empty_occurrence_file));

  const std::string malformed_occurrence_file =
    std::string(root) + "/malformed-occurrences";
  {
    std::ofstream output(malformed_occurrence_file.c_str(),
      std::ios_base::binary);
    appendInteger(output, 7); appendInteger(output, 1);
    appendInteger(output, 7); appendInteger(output, 1);
    require(static_cast<bool>(output));
  }
  bool rejected_malformed_occurrences = false;
  try
  {
    std::ostringstream ignored;
    serializeOccurrencePointers(ignored, malformed_occurrence_file,
      10, 2, 2, parameters);
  }
  catch(const std::runtime_error&) { rejected_malformed_occurrences = true; }
  require(rejected_malformed_occurrences);
  require(std::filesystem::remove(malformed_occurrence_file));

  // Direct packing preserves the public serialization exactly while keeping
  // component construction local to the writer. Loading proves that the
  // stream order still matches GCSA::load().
  const std::string packed_name = std::string(root) + "/packed.gcsa";
  storeFinalComponents(observed.header, alphabet, files, metadata,
    parameters, packed_name);
  std::ifstream packed_input(packed_name.c_str(), std::ios_base::binary);
  std::ostringstream packed_bytes; packed_bytes << packed_input.rdbuf();
  require(packed_bytes.str() == expected_bytes.str());
  std::istringstream packed_stream(packed_bytes.str());
  GCSA packed; packed.load(packed_stream);
  require(packed.extra_pointers.count(2, 2) == 20);
  require(packed.redundant_pointers.count(2, 2) == 2);
  require(packed.sampleCount() == 3 && packed.sample(2) == 100);

  // The internal construction caller may admit independent component
  // encoders. Their spools must still publish the exact ordinary byte stream.
  ConstructionParameters parallel_parameters = parameters;
  parallel_parameters.setMemoryLimitBytes(192 * KILOBYTE);
  const std::string parallel_name = std::string(root) + "/parallel.gcsa";
  storeFinalComponentsConcurrent(observed.header, alphabet, files, metadata,
    parallel_parameters, parallel_name, 4);
  require(readFile(parallel_name) == expected_bytes.str());
  GCSA parallel;
  require(sdsl::load_from_file(parallel, parallel_name));
  require(parallel.extra_pointers.count(2, 2) == 20);
  require(parallel.redundant_pointers.count(2, 2) == 2);
  require(parallel.sampleCount() == 3 && parallel.sample(2) == 100);
  require(std::filesystem::remove(parallel_name));

  const std::string memory_fallback_name =
    std::string(root) + "/memory-fallback.gcsa";
  storeFinalComponentsConcurrent(observed.header, alphabet, files, metadata,
    parameters, memory_fallback_name, 4);
  require(readFile(memory_fallback_name) == expected_bytes.str());
  require(std::filesystem::remove(memory_fallback_name));

  // A low process descriptor limit admits the unchanged serial route before
  // starting workers. The resulting bytes remain identical.
  const std::string fallback_name = std::string(root) + "/fallback.gcsa";
  pid_t fallback_child = ::fork(); require(fallback_child >= 0);
  if(fallback_child == 0)
  {
    struct rlimit limit = { 16, 16 };
    if(::setrlimit(RLIMIT_NOFILE, &limit) != 0) { _exit(2); }
    try
    {
      storeFinalComponentsConcurrent(observed.header, alphabet, files,
        metadata, parallel_parameters, fallback_name, 4);
    }
    catch(...) { _exit(3); }
    _exit(0);
  }
  int fallback_status = 0;
  require(::waitpid(fallback_child, &fallback_status, 0) == fallback_child);
  require(WIFEXITED(fallback_status) && WEXITSTATUS(fallback_status) == 0);
  require(readFile(fallback_name) == expected_bytes.str());
  require(std::filesystem::remove(fallback_name));

  // A semantic failure reached by a component worker must cancel and join
  // its peers, remove every spool, and leave no partial publication.
  writeInteger(files.occurrences, 8,
    std::numeric_limits<std::uint64_t>::max());
  const std::string failed_name = std::string(root) + "/failed.gcsa";
  bool parallel_failure_rejected = false;
  try
  {
    storeFinalComponentsConcurrent(observed.header, alphabet, files, metadata,
      parallel_parameters, failed_name, 4);
  }
  catch(const std::runtime_error&) { parallel_failure_rejected = true; }
  require(parallel_failure_rejected);
  require(!std::filesystem::exists(failed_name));
  require(!std::filesystem::exists(failed_name + "." +
    std::to_string(static_cast<unsigned long long>(::getpid())) + ".partial"));
  require(!hasComponentSpool(root));
  writeInteger(files.occurrences, 8, 2);



  const auto rejects = [&](const FinalEventMetadata& candidate_metadata)
  {
    try
    {
      buildViaStore(alphabet, files, candidate_metadata, parameters, 4, 6, 8);
    }
    catch(const std::runtime_error&) { return true; }
    return false;
  };

  // Metadata and streams are validated before unchecked succinct-builder
  // assumptions can become out-of-bounds writes.
  FinalEventMetadata invalid_metadata = metadata;
  invalid_metadata.sample_bits = 65;
  require(rejects(invalid_metadata));

  // Path 1 originally has only component 1. Adding a second component-5 bit
  // exceeds the declared sparse-BWT one-count and must fail closed.
  writeByte(files.bwt_masks, 1, static_cast<std::uint8_t>((1U << 1) | (1U << 5)));
  require(rejects(metadata));
  writeByte(files.bwt_masks, 1, static_cast<std::uint8_t>(1U << 1));

  // A corrupt occurrence value larger than the declared total previously
  // underflowed unsigned subtraction before set_unsafe().
  writeInteger(files.occurrences, 8,
    std::numeric_limits<std::uint64_t>::max());
  require(rejects(metadata));
  writeInteger(files.occurrences, 8, 2);

  invalid_metadata = metadata;
  invalid_metadata.sample_bits = 1;
  require(rejects(invalid_metadata));

  // Raw truncation is rejected before SDSL sees an inconsistent universe.
  require(::truncate(files.occurrences.c_str(), 16) == 0);
  bool rejected = false;
  try { buildViaStore(alphabet, files, metadata, parameters, 4, 6, 8); }
  catch(const std::runtime_error&) { rejected = true; }
  require(rejected);

  files.clear();
  std::filesystem::remove_all(root);
  return 0;
}
