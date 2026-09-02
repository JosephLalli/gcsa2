#include <gcsa/algorithms.h>
#include <gcsa/files.h>
#include <gcsa/gcsa.h>
#include <gcsa/lcp.h>

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
require(bool condition)
{
  if(!condition) { std::exit(EXIT_FAILURE); }
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
  parameters.setSortRunSize(32 * KILOBYTE);
  parameters.setJoinPartitionSize(64 * KILOBYTE);
  parameters.setMergeFanIn(2);
  parameters.setMaxOpenFiles(8);
  return parameters;
}

} // namespace

int
main()
{
  omp_set_num_threads(1);
  Verbosity::set(Verbosity::SILENT);

  char legacy_root[] = "/tmp/gcsa-resume-legacy-XXXXXX";
  char workspace_root[] = "/tmp/gcsa-resume-work-XXXXXX";
  require(mkdtemp(legacy_root) != nullptr);
  require(mkdtemp(workspace_root) != nullptr);
  const std::string input_name = "tests/cycle.gcsa2";
  const std::string legacy_prefix = std::string(legacy_root) + "/index";
  const std::string external_prefix = std::string(workspace_root) + "/index";

  // Establish the byte-for-byte reference through the original in-memory
  // extend route. The common external label sorter is still exercised.
  TempFile::setDirectory(legacy_root);
  {
    ConstructionParameters parameters;
    parameters.setSteps(2);
    parameters.setMemoryLimitBytes(MEGABYTE);
    InputGraph graph({ input_name }, false, parameters);
    GCSA index(graph, parameters);
    LCPArray lcp(graph, parameters);
    require(verifyIndex(index, &lcp, graph));
    store(index, lcp, legacy_prefix);
  }

  // Stop after a durable, fine-grained prefix-doubling checkpoint.
  TempFile::setDirectory(workspace_root);
  {
    ConstructionParameters parameters = externalParameters(workspace_root, 128 * KILOBYTE);
    parameters.setStopAfter("step-01-prune");
    bool stopped = false;
    try
    {
      InputGraph graph({ input_name }, false, parameters);
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

  // Resume with a different operational RAM ceiling. Semantic parameters and
  // input checksums stay fixed, while the committed prune is restored.
  {
    ConstructionParameters parameters = externalParameters(workspace_root, 96 * KILOBYTE);
    parameters.setResume();
    InputGraph graph({ input_name }, false, parameters);
    GCSA index(graph, parameters);
    LCPArray lcp(graph, parameters);
    require(verifyIndex(index, &lcp, graph));
    store(index, lcp, external_prefix);
  }

  require(readFile(legacy_prefix + GCSA::EXTENSION) ==
    readFile(external_prefix + GCSA::EXTENSION));
  require(readFile(legacy_prefix + LCPArray::EXTENSION) ==
    readFile(external_prefix + LCPArray::EXTENSION));

  // A semantic change must refuse reuse even though operational parameters are
  // deliberately allowed to change between invocations.
  bool semantic_change_refused = false;
  try
  {
    ConstructionParameters parameters = externalParameters(workspace_root, 96 * KILOBYTE);
    parameters.setSteps(1);
    parameters.setResume();
    InputGraph graph({ input_name }, false, parameters);
    GCSA index(graph, parameters);
  }
  catch(const std::runtime_error&)
  {
    semantic_change_refused = true;
  }
  require(semantic_change_refused);

  std::filesystem::remove_all(legacy_root);
  std::filesystem::remove_all(workspace_root);
  return 0;
}
