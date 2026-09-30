#include <gcsa/algorithms.h>
#include <gcsa/compressed_block.h>
#include <gcsa/path_graph_external.h>
#include <gcsa/files.h>
#include <gcsa/gcsa.h>
#include <gcsa/lcp.h>
#include <gcsa/path_graph.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>

using namespace gcsa;

namespace
{

void
requireAt(bool condition, size_type line)
{
  if(!condition)
  {
    std::cerr << "test_construction_resume: requirement failed at line "
              << line << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

#define require(condition) requireAt((condition), __LINE__)

std::string
makeTempRoot(const std::string& stem)
{
  const char* configured = std::getenv("TMPDIR");
  std::filesystem::path parent =
    (configured != nullptr && configured[0] != '\0' ? configured : "/tmp");
  std::string pattern = (parent / (stem + "-XXXXXX")).string();
  std::vector<char> mutable_pattern(pattern.begin(), pattern.end());
  mutable_pattern.push_back('\0');
  require(mkdtemp(mutable_pattern.data()) != nullptr);
  return std::string(mutable_pattern.data());
}

std::vector<char>
readFile(const std::string& filename)
{
  std::ifstream input(filename.c_str(), std::ios_base::binary);
  require(static_cast<bool>(input));
  return std::vector<char>(std::istreambuf_iterator<char>(input),
    std::istreambuf_iterator<char>());
}

void
store(const GCSA& index, const LCPArray& lcp, const std::string& prefix)
{
  require(sdsl::store_to_file(index, prefix + GCSA::EXTENSION));
  require(sdsl::store_to_file(lcp, prefix + LCPArray::EXTENSION));
}

ConstructionParameters
externalParameters(const std::string& workspace, size_type memory_limit)
{
  ConstructionParameters parameters;
  parameters.setSteps(2);
  parameters.setWorkDirectory(workspace);
  parameters.setMemoryLimitBytes(memory_limit);
  parameters.setIOBufferSize(KILOBYTE);
  parameters.setSortRunSize(externalPathGraphSortMinimumBudget());
  parameters.setJoinPartitionSize(externalPathJoinMinimumBudget());
  parameters.setMergeFanIn(2);
  parameters.setMaxOpenFiles(ConstructionParameters::MIN_OPEN_FILES);
  return parameters;
}

bool
hasWorkspaceArtifact(const std::string& workspace, const std::string& prefix)
{
  for(const std::filesystem::directory_entry& entry :
      std::filesystem::directory_iterator(workspace))
  {
    const std::string name = entry.path().filename().string();
    if(entry.is_regular_file() && name.compare(0, prefix.size(), prefix) == 0 &&
       entry.path().extension() == ".bin")
    {
      return true;
    }
  }
  return false;
}

std::map<std::string, std::vector<char>>
workspaceArtifacts(const std::string& workspace, const std::string& prefix)
{
  std::map<std::string, std::vector<char>> result;
  for(const std::filesystem::directory_entry& entry :
      std::filesystem::directory_iterator(workspace))
  {
    const std::string name = entry.path().filename().string();
    if(entry.is_regular_file() && name.compare(0, prefix.size(), prefix) == 0 &&
       entry.path().extension() == ".bin")
    {
      result[name] = readFile(entry.path().string());
    }
  }
  return result;
}

std::vector<char>
readLogicalFile(const std::string& filename)
{
  if(!CompressedBlockReader::isFramed(filename)) { return readFile(filename); }
  CompressedBlockReader input(filename);
  require(input.logicalSize() <= std::numeric_limits<size_type>::max());
  std::vector<char> result(static_cast<size_type>(input.logicalSize()));
  if(!result.empty())
  {
    require(input.read(result.data(), result.size()) == result.size());
  }
  char extra = 0;
  require(input.read(&extra, 1) == 0);
  return result;
}

std::map<std::string, std::vector<char>>
workspaceLogicalArtifacts(const std::string& workspace,
  const std::string& prefix)
{
  std::map<std::string, std::vector<char>> result;
  for(const std::filesystem::directory_entry& entry :
      std::filesystem::directory_iterator(workspace))
  {
    const std::string name = entry.path().filename().string();
    if(entry.is_regular_file() && name.compare(0, prefix.size(), prefix) == 0 &&
       entry.path().extension() == ".bin")
    {
      result[name] = readLogicalFile(entry.path().string());
    }
  }
  return result;
}

void
enableSmallFramedCodec(ConstructionParameters& parameters)
{
  parameters.setTempCompression("zstd");
  parameters.setCompressionBlockSize(64 * KILOBYTE);
  parameters.setCompressionWorkers(1);
  // The raw minima used by externalParameters() do not include zstd's
  // writer and reader workspaces.  These established framed-join caps admit
  // the minimum supported codec block while remaining small enough to
  // exercise external runs.
  parameters.setJoinPartitionSize(4 * MEGABYTE);
  parameters.setSortRunSize(2 * MEGABYTE);
}

void
enableSmallFramedStreams(ConstructionParameters& parameters)
{
  parameters.setIOBufferSize(16 * KILOBYTE);
  enableSmallFramedCodec(parameters);
}

void
appendDeBruijn(size_type t, size_type period, size_type order,
  size_type alphabet_size, std::vector<size_type>& state,
  std::vector<size_type>& sequence)
{
  if(t > order)
  {
    if(order % period == 0)
    {
      for(size_type i = 1; i <= period; i++) { sequence.push_back(state[i]); }
    }
    return;
  }
  state[t] = state[t - period];
  appendDeBruijn(t + 1, period, order, alphabet_size, state, sequence);
  for(size_type value = state[t - period] + 1; value < alphabet_size; value++)
  {
    state[t] = value;
    appendDeBruijn(t + 1, t, order, alphabet_size, state, sequence);
  }
}

void
writeSparseSkewCycle(const std::string& filename)
{
  constexpr size_type order = 6;
  const std::string alphabet = "ACGT";
  std::vector<size_type> state(alphabet.size() * order + 1, 0), sequence;
  appendDeBruijn(1, 1, order, alphabet.size(), state, sequence);
  require(sequence.size() == 4096);

  std::ofstream output(filename.c_str());
  require(static_cast<bool>(output));
  // Reuse 257 graph positions around the 4096-label cycle. The same from node
  // therefore recurs in many nonadjacent final paths and across every bounded
  // preparation batch, while the de Bruijn order supplies changing LCP-stack
  // boundaries. This makes byte equality exercise ordered prev-occurrence and
  // suffix-stack state rather than merely disjoint per-path sets.
  constexpr size_type graph_positions = 257;
  for(size_type i = 0; i < sequence.size(); i++)
  {
    std::string label;
    for(size_type j = 0; j < order; j++)
    {
      label.push_back(alphabet[sequence[(i + j) % sequence.size()]]);
    }
    char predecessor = alphabet[sequence[(i + sequence.size() - 1) % sequence.size()]];
    char successor = alphabet[sequence[(i + order) % sequence.size()]];
    size_type current_node = i % graph_positions;
    size_type next = ((i + 1) % sequence.size()) % graph_positions;
    output << label << '\t' << (10 + current_node / 256) << ':'
           << (current_node % 256 + 1) << '\t'
           << predecessor << '\t' << successor << '\t'
           << (10 + next / 256) << ':' << (next % 256 + 1) << '\n';
  }

  // A separate N self-loop component carries enough distinct starts to exceed
  // the minimum SpillableNodeSet collection capacity. Each duplicate has a
  // real incoming and outgoing loop, keeping the input graph consistent while
  // its one-character predecessor mask remains sparse. The 4096 DNA paths force
  // many batches and mostly exercise the unsampled continuation comparison.
  const std::string skew_label(order, 'N');
  for(size_type id = 1; id <= 3; id++)
  {
    output << skew_label << '\t' << id << ":1\tN\tN\t"
           << id << ":1\n";
  }
  for(size_type copy = 0; copy < 2048; copy++)
  {
    size_type id = 1000 + copy;
    output << skew_label << '\t' << id << ":1\tN\tN\t"
           << id << ":1\n";
  }
  output.close(); require(static_cast<bool>(output));
}

} // namespace

int
main()
{
  unsetenv("GCSA_EXPERIMENTAL_FINAL_STATE");
  unsetenv("GCSA_EXPERIMENTAL_STATE_WORKERS");
  omp_set_num_threads(1);
  Verbosity::set(Verbosity::SILENT);

  const std::string legacy_root = makeTempRoot("gcsa-resume-legacy");
  const std::string legacy_parallel_root =
    makeTempRoot("gcsa-resume-legacy-parallel");
  const std::string workspace_root = makeTempRoot("gcsa-resume-work");
  const std::string parallel_root = makeTempRoot("gcsa-resume-parallel");
  const std::string state_one_root = makeTempRoot("gcsa-resume-state-one");
  const std::string state_parallel_root =
    makeTempRoot("gcsa-resume-state-parallel");
  const std::string fallback_root = makeTempRoot("gcsa-resume-fallback");
  const std::string empty_root = makeTempRoot("gcsa-resume-empty");
  const std::string input_name = legacy_root + "/sparse-skew.gcsa2";
  writeSparseSkewCycle(input_name);
  const std::string legacy_prefix = std::string(legacy_root) + "/index";
  const std::string legacy_parallel_prefix =
    std::string(legacy_parallel_root) + "/index";
  const std::string external_prefix = std::string(workspace_root) + "/index";
  const std::string parallel_prefix = std::string(parallel_root) + "/index";
  const std::string state_one_prefix = std::string(state_one_root) + "/index";
  const std::string state_parallel_prefix =
    std::string(state_parallel_root) + "/index";
  const std::string fallback_prefix = std::string(fallback_root) + "/index";
  const std::string staged_prefix = std::string(workspace_root) + "/staged";
  const std::string mapping_name = std::string(legacy_root) + "/mapping";
  // Exercise mapping-before-deduplication in the final start-node set: the
  // mapping is deliberately non-monotone and collapses two input node ids.
  NodeMapping mapping(1);
  mapping.insert(2); mapping.insert(2); mapping.insert(1);
  require(sdsl::store_to_file(mapping, mapping_name));
  // Prefix doubling streams join results directly into label sorting, so the
  // configured ceiling must admit both independently byte-bounded stages.
  const size_type external_minimum = externalPathJoinMinimumBudget() +
    externalPathGraphSortMinimumBudget();

  // Establish the byte-for-byte reference through the original in-memory
  // extend route. The common external label sorter is still exercised.
  TempFile::setDirectory(legacy_root);
  {
    ConstructionParameters parameters;
    parameters.setSteps(2);
    parameters.setMemoryLimitBytes(MEGABYTE);
    InputGraph graph({ input_name }, false, parameters, Alphabet(), mapping_name);
    GCSA index(graph, parameters);
    LCPArray lcp(graph, parameters);
    require(verifyIndex(index, &lcp, graph));
    store(index, lcp, legacy_prefix);
  }

  // The fixture has thousands of distinct labels joined across two doubling
  // rounds. A caller may still configure OpenMP threads for the legacy route,
  // but its extender requires one physical shard per logical input. Verify the
  // constructor keeps prune serial and preserves the public oracle bytes.
  omp_set_num_threads(4);
  TempFile::setDirectory(legacy_parallel_root);
  {
    ConstructionParameters parameters;
    parameters.setSteps(2);
    parameters.setMemoryLimitBytes(MEGABYTE);
    InputGraph graph({ input_name }, false, parameters, Alphabet(), mapping_name);
    GCSA index(graph, parameters);
    LCPArray lcp(graph, parameters);
    require(verifyIndex(index, &lcp, graph));
    store(index, lcp, legacy_parallel_prefix);
  }
  require(readFile(legacy_prefix + GCSA::EXTENSION) ==
    readFile(legacy_parallel_prefix + GCSA::EXTENSION));
  require(readFile(legacy_prefix + LCPArray::EXTENSION) ==
    readFile(legacy_parallel_prefix + LCPArray::EXTENSION));
  omp_set_num_threads(1);

  // Stop after a durable, fine-grained prefix-doubling checkpoint.
  TempFile::setDirectory(workspace_root);
  {
    ConstructionParameters parameters = externalParameters(workspace_root,
      external_minimum + 256 * KILOBYTE);
    parameters.setStopAfter("step-01-prune");
    bool stopped = false;
    try
    {
      InputGraph graph({ input_name }, false, parameters, Alphabet(), mapping_name);
      GCSA index(graph, parameters);
    }
    catch(const ConstructionStopped& event)
    {
      stopped = (event.completed_phase == "step-01-prune");
    }
    require(stopped);
    require(std::filesystem::exists(
      std::string(workspace_root) + "/step-01--prune.complete"));
  }

  // Resume through the final ordered scan, but stop before component assembly.
  // This verifies that a completed event set is independently durable.
  std::ostringstream serial_log;
  {
    ConstructionParameters parameters = externalParameters(workspace_root,
      external_minimum + 64 * MEGABYTE);
    enableSmallFramedStreams(parameters);
    parameters.setResume();
    // Cleanup is operational and may be enabled only after an older frontier
    // already exists. Resumption must catch up without discarding this step's
    // reusable join tasks.
    parameters.setCleanObsolete();
    parameters.setStopAfter("final-events");
    bool stopped = false;
    std::streambuf* old_stderr = std::cerr.rdbuf(serial_log.rdbuf());
    Verbosity::set(Verbosity::EXTENDED);
    try
    {
      InputGraph graph({ input_name }, false, parameters, Alphabet(), mapping_name);
      GCSA index(graph, parameters);
    }
    catch(const ConstructionStopped& event)
    {
      stopped = (event.completed_phase == "final-events");
    }
    catch(...)
    {
      std::cerr.rdbuf(old_stderr); Verbosity::set(Verbosity::SILENT); throw;
    }
    std::cerr.rdbuf(old_stderr); Verbosity::set(Verbosity::SILENT);
    require(stopped);
    require(std::filesystem::exists(
      std::string(workspace_root) + "/final--events.complete"));
    require(!hasWorkspaceArtifact(workspace_root, "initial--paths--"));
  }

  // Rebuild the same event frontier with component-wise predecessor workers.
  // A 16 KiB event buffer forces framed event payloads, while 64 descriptors
  // is the supported minimum and the small sort/join caps retain forced
  // spilling.
  // The final event artifacts themselves must be byte-identical before either
  // workspace is allowed to assemble the public index.
  omp_set_num_threads(4);
  TempFile::setDirectory(parallel_root);
  std::ostringstream parallel_log;
  {
    ConstructionParameters parameters = externalParameters(parallel_root,
      external_minimum + 64 * MEGABYTE);
    enableSmallFramedStreams(parameters);
    InputGraph graph({ input_name }, false, parameters, Alphabet(), mapping_name);
    std::streambuf* old_stderr = std::cerr.rdbuf(parallel_log.rdbuf());
    Verbosity::set(Verbosity::EXTENDED);
    try
    {
      GCSA::buildAndStore(graph, parameters,
        parallel_prefix + GCSA::EXTENSION);
      LCPArray::buildAndStore(graph, parameters,
        parallel_prefix + LCPArray::EXTENSION);
    }
    catch(...)
    {
      std::cerr.rdbuf(old_stderr); Verbosity::set(Verbosity::SILENT); throw;
    }
    std::cerr.rdbuf(old_stderr); Verbosity::set(Verbosity::SILENT);
    require(std::filesystem::exists(
      std::string(parallel_root) + "/final--events.complete"));
    GCSA index;
    LCPArray lcp;
    require(sdsl::load_from_file(index,
      parallel_prefix + GCSA::EXTENSION));
    require(sdsl::load_from_file(lcp,
      parallel_prefix + LCPArray::EXTENSION));
    require(verifyIndex(index, &lcp, graph));
  }
  const std::string predecessor_prefix =
    "final predecessor lookup: 4 worker(s), ";
  size_type predecessor_line = parallel_log.str().find(predecessor_prefix);
  require(predecessor_line != std::string::npos);
  size_type predecessor_batches = std::stoull(parallel_log.str().substr(
    predecessor_line + predecessor_prefix.size()));
  require(predecessor_batches > 1);
  const std::string preparation_prefix =
    "final ordered-state preparation: 4 worker(s), ";
  size_type preparation_line = parallel_log.str().find(preparation_prefix);
  require(preparation_line != std::string::npos);
  size_type prepared_paths = std::stoull(parallel_log.str().substr(
    preparation_line + preparation_prefix.size()));
  require(prepared_paths >= 4096);
  require(parallel_log.str().find("1 spill fallback(s)", preparation_line) !=
    std::string::npos);
  const auto scanReadBytes = [](const std::string& log) -> size_type {
    const std::string prefix = "final scan read bytes: ";
    size_type at = log.find(prefix); require(at != std::string::npos);
    return std::stoull(log.substr(at + prefix.size()));
  };
  const size_type serial_scan_reads = scanReadBytes(serial_log.str());
  const size_type parallel_scan_reads = scanReadBytes(parallel_log.str());
  std::cerr << "final scan regression fixture read bytes: serial "
            << serial_scan_reads << ", parallel " << parallel_scan_reads << '\n';
  // Allow initial prefetch alignment, but not one full-buffer reread per
  // component and batch. Sampling and lookup must remain monotone streams.
  require(parallel_scan_reads <= serial_scan_reads + 1024 * KILOBYTE);
  const std::string spill_prefix = "final start-node sets: ";
  size_type spill_line = parallel_log.str().find(spill_prefix);
  require(spill_line != std::string::npos);
  size_type spills = std::stoull(parallel_log.str().substr(
    spill_line + spill_prefix.size()));
  require(spills > 0);
  std::map<std::string, std::vector<char>> serial_events =
    workspaceArtifacts(workspace_root, "final--events--");
  std::map<std::string, std::vector<char>> parallel_events =
    workspaceArtifacts(parallel_root, "final--events--");
  require(!serial_events.empty());
  require(serial_events == parallel_events);

  require(readFile(legacy_prefix + GCSA::EXTENSION) ==
    readFile(parallel_prefix + GCSA::EXTENSION));
  require(readFile(legacy_prefix + LCPArray::EXTENSION) ==
    readFile(parallel_prefix + LCPArray::EXTENSION));

  // Exercise the opt-in final-state planner through the same bounded batches
  // with one and four state workers. The 257 recurring ranks skip many batch
  // boundaries, the mapping collapses duplicate nodes before rank dispatch,
  // and the oversized N path forces the serial spill path to update the same
  // owner shard. Sorted event checkpoints and both public products must remain
  // byte-identical to the legacy and default external routes.
  const auto runExperimentalState = [&](const std::string& root,
    const std::string& prefix, const char* workers) -> std::string
  {
    require(setenv("GCSA_EXPERIMENTAL_FINAL_STATE", "1", 1) == 0);
    require(setenv("GCSA_EXPERIMENTAL_STATE_WORKERS", workers, 1) == 0);
    TempFile::setDirectory(root);
    std::ostringstream log;
    ConstructionParameters parameters = externalParameters(root,
      external_minimum + 64 * MEGABYTE);
    enableSmallFramedStreams(parameters);
    InputGraph graph({ input_name }, false, parameters, Alphabet(), mapping_name);
    std::streambuf* old_stderr = std::cerr.rdbuf(log.rdbuf());
    Verbosity::set(Verbosity::EXTENDED);
    try
    {
      GCSA::buildAndStore(graph, parameters, prefix + GCSA::EXTENSION);
      LCPArray::buildAndStore(graph, parameters, prefix + LCPArray::EXTENSION);
    }
    catch(...)
    {
      std::cerr.rdbuf(old_stderr); Verbosity::set(Verbosity::SILENT);
      unsetenv("GCSA_EXPERIMENTAL_FINAL_STATE");
      unsetenv("GCSA_EXPERIMENTAL_STATE_WORKERS");
      throw;
    }
    std::cerr.rdbuf(old_stderr); Verbosity::set(Verbosity::SILENT);
    unsetenv("GCSA_EXPERIMENTAL_FINAL_STATE");
    unsetenv("GCSA_EXPERIMENTAL_STATE_WORKERS");
    GCSA index;
    LCPArray lcp;
    require(sdsl::load_from_file(index, prefix + GCSA::EXTENSION));
    require(sdsl::load_from_file(lcp, prefix + LCPArray::EXTENSION));
    require(verifyIndex(index, &lcp, graph));
    return log.str();
  };

  const std::string state_one_log = runExperimentalState(
    state_one_root, state_one_prefix, "1");
  const std::string state_parallel_log = runExperimentalState(
    state_parallel_root, state_parallel_prefix, "4");
  require(state_one_log.find(
    "experimental final state: requested 1 worker(s), active 1 worker(s)") !=
    std::string::npos);
  require(state_parallel_log.find(
    "experimental final state: requested 4 worker(s), active 4 worker(s)") !=
    std::string::npos);
  require(state_one_log.find("experimental final-state timings:") !=
    std::string::npos);
  require(state_parallel_log.find("experimental final-state LCP replay:") !=
    std::string::npos);
  require(workspaceArtifacts(state_one_root, "final--events--") ==
    serial_events);
  require(workspaceArtifacts(state_parallel_root, "final--events--") ==
    serial_events);
  require(readFile(state_one_prefix + GCSA::EXTENSION) ==
    readFile(legacy_prefix + GCSA::EXTENSION));
  require(readFile(state_parallel_prefix + GCSA::EXTENSION) ==
    readFile(legacy_prefix + GCSA::EXTENSION));
  require(readFile(state_one_prefix + LCPArray::EXTENSION) ==
    readFile(legacy_prefix + LCPArray::EXTENSION));
  require(readFile(state_parallel_prefix + LCPArray::EXTENSION) ==
    readFile(legacy_prefix + LCPArray::EXTENSION));

  // Request workers under the supported descriptor floor, but keep the 1 KiB
  // I/O buffer below the minimum useful batch. The final scan must retain the
  // serial reader layout and produce the same framed event bytes and indexes.
  require(setenv("GCSA_EXPERIMENTAL_FINAL_STATE", "1", 1) == 0);
  require(setenv("GCSA_EXPERIMENTAL_STATE_WORKERS", "4", 1) == 0);
  TempFile::setDirectory(fallback_root);
  std::ostringstream fallback_log;
  {
    ConstructionParameters parameters = externalParameters(fallback_root,
      external_minimum + 64 * MEGABYTE);
    enableSmallFramedCodec(parameters);
    InputGraph graph({ input_name }, false, parameters, Alphabet(), mapping_name);
    std::streambuf* old_stderr = std::cerr.rdbuf(fallback_log.rdbuf());
    Verbosity::set(Verbosity::EXTENDED);
    try
    {
      GCSA::buildAndStore(graph, parameters,
        fallback_prefix + GCSA::EXTENSION);
      LCPArray::buildAndStore(graph, parameters,
        fallback_prefix + LCPArray::EXTENSION);
    }
    catch(...)
    {
      std::cerr.rdbuf(old_stderr); Verbosity::set(Verbosity::SILENT); throw;
    }
    std::cerr.rdbuf(old_stderr); Verbosity::set(Verbosity::SILENT);
    GCSA index;
    LCPArray lcp;
    require(sdsl::load_from_file(index,
      fallback_prefix + GCSA::EXTENSION));
    require(sdsl::load_from_file(lcp,
      fallback_prefix + LCPArray::EXTENSION));
    require(verifyIndex(index, &lcp, graph));
  }
  unsetenv("GCSA_EXPERIMENTAL_FINAL_STATE");
  unsetenv("GCSA_EXPERIMENTAL_STATE_WORKERS");
  require(fallback_log.str().find(
    "final predecessor lookup: 1 worker(s), 0 bounded batch(es)") !=
    std::string::npos);
  require(fallback_log.str().find(
    "serial fallback: batch admission too small") != std::string::npos);
  require(fallback_log.str().find(
    "final ordered-state preparation: 1 worker(s), 0 paths / 0 from-ranks prepared") !=
    std::string::npos);
  require(fallback_log.str().find(
    "experimental final state: requested 4 worker(s), active 1 worker(s)") !=
    std::string::npos);
  require(fallback_log.str().find(
    "serial fallback: predecessor batch admission") != std::string::npos);
  // The final-event writer bounds its frame size by the I/O buffer, so this
  // deliberate 1 KiB fallback uses different physical frame boundaries than
  // the 16 KiB serial run. Compare decoded streams here; the matched-buffer
  // serial/parallel comparison above remains byte-for-byte.
  require(workspaceLogicalArtifacts(fallback_root, "final--events--") ==
    workspaceLogicalArtifacts(workspace_root, "final--events--"));
  require(readFile(legacy_prefix + GCSA::EXTENSION) ==
    readFile(fallback_prefix + GCSA::EXTENSION));
  require(readFile(legacy_prefix + LCPArray::EXTENSION) ==
    readFile(fallback_prefix + LCPArray::EXTENSION));

  omp_set_num_threads(1);
  TempFile::setDirectory(workspace_root);

  // A second resume changes the operational RAM ceiling again, restores the
  // immutable events, and assembles byte-identical final components.
  {
    ConstructionParameters parameters = externalParameters(workspace_root,
      external_minimum + MEGABYTE);
    parameters.setResume();
    InputGraph graph({ input_name }, false, parameters, Alphabet(), mapping_name);
    // Production builds both external products through buildAndStore, not
    // through the constructors: build_gcsa and both vg entry points branch at
    // the call site. Exercise the paths that ship, then load what they
    // published.
    const std::string external_gcsa = external_prefix + GCSA::EXTENSION;
    const std::string external_lcp = external_prefix + LCPArray::EXTENSION;
    GCSA::buildAndStore(graph, parameters, external_gcsa);
    LCPArray::buildAndStore(graph, parameters, external_lcp);
    GCSA index; require(sdsl::load_from_file(index, external_gcsa));
    LCPArray lcp; require(sdsl::load_from_file(lcp, external_lcp));
    require(verifyIndex(index, &lcp, graph));
  }

  require(readFile(legacy_prefix + GCSA::EXTENSION) ==
    readFile(external_prefix + GCSA::EXTENSION));
  require(readFile(legacy_prefix + LCPArray::EXTENSION) ==
    readFile(external_prefix + LCPArray::EXTENSION));

  // Reuse the same completed event frontier through the production direct
  // packer. This exercises the constructor-to-file lifecycle, not just the
  // lower-level component helper, and proves ordinary GCSA::load()/queries
  // remain compatible.
  {
    ConstructionParameters parameters = externalParameters(workspace_root,
      external_minimum + MEGABYTE);
    parameters.setResume();
    InputGraph graph({ input_name }, false, parameters, Alphabet(), mapping_name);
    GCSA::buildAndStore(graph, parameters, staged_prefix + GCSA::EXTENSION);
    LCPArray::buildAndStore(graph, parameters, staged_prefix + LCPArray::EXTENSION);
    GCSA staged;
    LCPArray lcp;
    require(sdsl::load_from_file(staged, staged_prefix + GCSA::EXTENSION));
    require(sdsl::load_from_file(lcp, staged_prefix + LCPArray::EXTENSION));
    require(verifyIndex(staged, &lcp, graph));
  }
  require(readFile(legacy_prefix + GCSA::EXTENSION) ==
    readFile(staged_prefix + GCSA::EXTENSION));
  require(readFile(legacy_prefix + LCPArray::EXTENSION) ==
    readFile(staged_prefix + LCPArray::EXTENSION));

  // A semantic change must refuse reuse even though operational parameters are
  // deliberately allowed to change between invocations.
  bool semantic_change_refused = false;
  try
  {
    ConstructionParameters parameters = externalParameters(workspace_root,
      external_minimum + 128 * KILOBYTE);
    parameters.setSteps(1);
    parameters.setResume();
    InputGraph graph({ input_name }, false, parameters, Alphabet(), mapping_name);
    GCSA index(graph, parameters);
  }
  catch(const std::runtime_error&)
  {
    semantic_change_refused = true;
  }
  require(semantic_change_refused);

  // Direct construction must publish a valid empty index as well. Previously
  // the constructor returned before the staged packer ran, leaving no output
  // (or, worse, leaving a stale destination from an earlier build).
  {
    const std::string empty_input = std::string(empty_root) + "/empty.graph";
    const std::string legacy_empty = std::string(empty_root) + "/legacy.gcsa";
    const std::string staged_empty = std::string(empty_root) + "/staged.gcsa";
    std::vector<KMer> kmers;
    std::ofstream output(empty_input.c_str(), std::ios_base::binary);
    writeBinary(output, kmers, 4); output.close();
    require(static_cast<bool>(output));

    ConstructionParameters legacy_parameters;
    InputGraph legacy_graph({ empty_input }, true, legacy_parameters);
    GCSA legacy_empty_index(legacy_graph, legacy_parameters);
    require(sdsl::store_to_file(legacy_empty_index, legacy_empty));

    ConstructionParameters external_parameters = externalParameters(empty_root,
      external_minimum + 64 * KILOBYTE);
    InputGraph external_graph({ empty_input }, true, external_parameters);
    GCSA::buildAndStore(external_graph, external_parameters, staged_empty);
    LCPArray::buildAndStore(external_graph, external_parameters,
      staged_empty + LCPArray::EXTENSION);
    GCSA loaded;
    LCPArray loaded_lcp;
    require(sdsl::load_from_file(loaded, staged_empty));
    require(sdsl::load_from_file(loaded_lcp, staged_empty + LCPArray::EXTENSION));
    require(loaded.size() == 0);
    require(loaded_lcp.size() == 0);
    require(readFile(legacy_empty) == readFile(staged_empty));
  }

  std::filesystem::remove_all(legacy_root);
  std::filesystem::remove_all(legacy_parallel_root);
  std::filesystem::remove_all(workspace_root);
  std::filesystem::remove_all(parallel_root);
  std::filesystem::remove_all(state_one_root);
  std::filesystem::remove_all(state_parallel_root);
  std::filesystem::remove_all(fallback_root);
  std::filesystem::remove_all(empty_root);
  return 0;
}
