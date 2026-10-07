#include <gcsa/algorithms.h>
#include <gcsa/gcsa.h>
#include <gcsa/lcp.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

using namespace gcsa;

namespace
{

void
requireAt(bool condition, size_type line)
{
  if(!condition)
  {
    std::cerr << "test_lcp_overlap: requirement failed at line " << line << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

#define require(condition) requireAt((condition), __LINE__)

std::string
makeRoot(const std::string& stem)
{
  const char* configured = std::getenv("TMPDIR");
  std::filesystem::path parent =
    (configured != nullptr && configured[0] != '\0' ? configured : "/tmp");
  std::string pattern = (parent / (stem + "-XXXXXX")).string();
  std::vector<char> mutable_pattern(pattern.begin(), pattern.end());
  mutable_pattern.push_back('\0');
  require(::mkdtemp(mutable_pattern.data()) != nullptr);
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

std::string
writeCycle(const std::string& root)
{
  const std::string input = root + "/cycle.gcsa2";
  std::ofstream out(input.c_str());
  require(static_cast<bool>(out));
  out << "AC\t1:0\tG\tG\t2:0\n"
      << "CG\t2:0\tA\tA\t3:0\n"
      << "GA\t3:0\tC\tC\t1:0\n";
  out.close(); require(static_cast<bool>(out));
  return input;
}

ConstructionParameters
parameters(const std::string& root, size_type memory)
{
  ConstructionParameters result;
  result.setSteps(1);
  result.setWorkDirectory(root);
  result.setMemoryLimitBytes(memory);
  return result;
}

void
requireNoTransient(const std::string& root)
{
  for(const auto& entry : std::filesystem::recursive_directory_iterator(root))
  {
    if(!entry.is_regular_file()) { continue; }
    const std::string name = entry.path().filename().string();
    require(name.find(".staged") == std::string::npos);
    require(name.find(".partial") == std::string::npos);
  }
}

struct Outputs
{
  std::vector<char> gcsa, lcp;
};

Outputs
buildSuccess(const std::string& root, int threads, size_type memory,
  const std::string& expected_log)
{
  TempFile::setDirectory(root);
  omp_set_num_threads(threads);
  ConstructionParameters build_parameters = parameters(root, memory);
  const std::string input = writeCycle(root);
  const std::string gcsa_name = root + "/index.gcsa";
  const std::string lcp_name = root + "/index.lcp";
  std::ostringstream log;
  std::streambuf* old_stderr = std::cerr.rdbuf(log.rdbuf());
  const size_type old_verbosity = Verbosity::level;
  Verbosity::set(Verbosity::EXTENDED);
  try
  {
    InputGraph graph({ input }, false, build_parameters);
    GCSA::buildAndStore(graph, build_parameters, gcsa_name, lcp_name);
    GCSA index; LCPArray lcp;
    require(sdsl::load_from_file(index, gcsa_name));
    require(sdsl::load_from_file(lcp, lcp_name));
    require(verifyIndex(index, &lcp, graph));
  }
  catch(...)
  {
    std::cerr.rdbuf(old_stderr); Verbosity::set(old_verbosity); throw;
  }
  std::cerr.rdbuf(old_stderr); Verbosity::set(old_verbosity);
  require(log.str().find(expected_log) != std::string::npos);
  requireNoTransient(root);
  return { readFile(gcsa_name), readFile(lcp_name) };
}

} // namespace

int
main()
{
  const size_type overlap_memory = 512 * MEGABYTE;
  const std::string serial_root = makeRoot("gcsa-lcp-serial");
  const Outputs serial = buildSuccess(serial_root, 1, overlap_memory,
    "serial LCP fallback: only one construction thread");

  const std::string overlap_root = makeRoot("gcsa-lcp-overlap");
  const Outputs overlap = buildSuccess(overlap_root, 2, overlap_memory,
    "overlapping LCP storage with final component storage");
  require(serial.gcsa == overlap.gcsa); require(serial.lcp == overlap.lcp);

  const std::string fallback_root = makeRoot("gcsa-lcp-fallback");
  const Outputs fallback = buildSuccess(fallback_root, 2, 64 * MEGABYTE,
    "serial LCP fallback: aggregate memory target");
  require(serial.gcsa == fallback.gcsa); require(serial.lcp == fallback.lcp);

  // A failed GCSA target joins the LCP worker and leaves the prior LCP target
  // untouched; an early-finishing worker never publishes ahead of the GCSA.
  const std::string gcsa_failure_root = makeRoot("gcsa-lcp-gcsa-failure");
  TempFile::setDirectory(gcsa_failure_root); omp_set_num_threads(2);
  {
    ConstructionParameters build_parameters =
      parameters(gcsa_failure_root, overlap_memory);
    const std::string input = writeCycle(gcsa_failure_root);
    const std::string lcp_name = gcsa_failure_root + "/index.lcp";
    {
      std::ofstream old(lcp_name.c_str(), std::ios_base::binary);
      old << "prior-lcp"; old.close(); require(static_cast<bool>(old));
    }
    bool failed = false;
    try
    {
      InputGraph graph({ input }, false, build_parameters);
      GCSA::buildAndStore(graph, build_parameters,
        gcsa_failure_root + "/missing/index.gcsa", lcp_name);
    }
    catch(const std::exception&) { failed = true; }
    require(failed);
    require(readFile(lcp_name) == std::vector<char>(
      { 'p', 'r', 'i', 'o', 'r', '-', 'l', 'c', 'p' }));
  }
  requireNoTransient(gcsa_failure_root);

  // A failed LCP target is reported only after the GCSA is durable. The leaf
  // remains owned by InputGraph, allowing the ordinary constructor to retry.
  const std::string lcp_failure_root = makeRoot("gcsa-lcp-lcp-failure");
  TempFile::setDirectory(lcp_failure_root); omp_set_num_threads(2);
  {
    ConstructionParameters build_parameters =
      parameters(lcp_failure_root, overlap_memory);
    const std::string input = writeCycle(lcp_failure_root);
    const std::string gcsa_name = lcp_failure_root + "/index.gcsa";
    InputGraph graph({ input }, false, build_parameters);
    bool failed = false;
    try
    {
      GCSA::buildAndStore(graph, build_parameters, gcsa_name,
        lcp_failure_root + "/missing/index.lcp");
    }
    catch(const std::exception&) { failed = true; }
    require(failed); require(std::filesystem::exists(gcsa_name));
    require(readFile(gcsa_name) == serial.gcsa);
    require(!(graph.lcp_name.empty()));
    require(std::filesystem::exists(graph.lcp_name));
    GCSA index; require(sdsl::load_from_file(index, gcsa_name));
    LCPArray recovered(graph, build_parameters);
    require(verifyIndex(index, &recovered, graph));
  }
  requireNoTransient(lcp_failure_root);

  std::filesystem::remove_all(serial_root);
  std::filesystem::remove_all(overlap_root);
  std::filesystem::remove_all(fallback_root);
  std::filesystem::remove_all(gcsa_failure_root);
  std::filesystem::remove_all(lcp_failure_root);
  return EXIT_SUCCESS;
}
