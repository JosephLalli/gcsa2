#include <gcsa/algorithms.h>
#include <gcsa/files.h>
#include <gcsa/gcsa.h>
#include <gcsa/lcp.h>
#include <gcsa/path_graph.h>

#include <filesystem>
#include <fstream>
#include <iterator>
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

} // namespace

int
main()
{
  omp_set_num_threads(1);
  Verbosity::set(Verbosity::SILENT);

  char legacy_root[] = "/tmp/gcsa-resume-legacy-XXXXXX";
  char workspace_root[] = "/tmp/gcsa-resume-work-XXXXXX";
  char empty_root[] = "/tmp/gcsa-resume-empty-XXXXXX";
  require(mkdtemp(legacy_root) != nullptr);
  require(mkdtemp(workspace_root) != nullptr);
  require(mkdtemp(empty_root) != nullptr);
  const std::string input_name = "tests/cycle.gcsa2";
  const std::string legacy_prefix = std::string(legacy_root) + "/index";
  const std::string external_prefix = std::string(workspace_root) + "/index";
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
  {
    ConstructionParameters parameters = externalParameters(workspace_root,
      external_minimum + 128 * KILOBYTE);
    parameters.setResume();
    // Cleanup is operational and may be enabled only after an older frontier
    // already exists. Resumption must catch up without discarding this step's
    // reusable join tasks.
    parameters.setCleanObsolete();
    parameters.setStopAfter("final-events");
    bool stopped = false;
    try
    {
      InputGraph graph({ input_name }, false, parameters, Alphabet(), mapping_name);
      GCSA index(graph, parameters);
    }
    catch(const ConstructionStopped& event)
    {
      stopped = (event.completed_phase == "final-events");
    }
    require(stopped);
    require(std::filesystem::exists(
      std::string(workspace_root) + "/final--events.complete"));
    require(!hasWorkspaceArtifact(workspace_root, "initial--paths--"));
  }

  // A second resume changes the operational RAM ceiling again, restores the
  // immutable events, and assembles byte-identical final components.
  {
    ConstructionParameters parameters = externalParameters(workspace_root,
      external_minimum + 64 * KILOBYTE);
    parameters.setResume();
    InputGraph graph({ input_name }, false, parameters, Alphabet(), mapping_name);
    GCSA index(graph, parameters);
    LCPArray lcp(graph, parameters);
    require(verifyIndex(index, &lcp, graph));
    store(index, lcp, external_prefix);
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
      external_minimum + 64 * KILOBYTE);
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
  std::filesystem::remove_all(workspace_root);
  std::filesystem::remove_all(empty_root);
  return 0;
}
