#include <gcsa/algorithms.h>
#include <gcsa/files.h>
#include <gcsa/gcsa.h>
#include <gcsa/lcp.h>

#include <cstdio>
#include <cstdlib>
#include <dirent.h>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

using namespace gcsa;

namespace
{

void require(bool value) { if(!value) { std::abort(); } }

void
writeFixture(const std::string& name, size_type copies, bool corrupt = false)
{
  std::ofstream output(name.c_str(), std::ios_base::trunc); require(static_cast<bool>(output));
  for(size_type i = 0; i < copies; i++)
  {
    output << "AC\t" << (corrupt ? "4:0" : "1:0") << "\tG\tG\t2:0\n";
    output << "CG\t2:0\tA\tA\t3:0\n";
    output << "GA\t3:0\tC\tC\t1:0\n";
  }
  output.close(); require(static_cast<bool>(output));
}

size_type
directoryEntries(const std::string& name)
{
  DIR* directory = opendir(name.c_str()); require(directory != nullptr);
  size_type result = 0;
  while(dirent* entry = readdir(directory))
  {
    std::string value(entry->d_name);
    if(value != "." && value != "..") { result++; }
  }
  closedir(directory); return result;
}

struct VerifyResult
{
  bool success;
  std::string output, errors;
};

template<class Verify>
VerifyResult
captureVerification(const Verify& verify)
{
  std::ostringstream output, errors;
  std::streambuf* old_output = std::cout.rdbuf(output.rdbuf());
  std::streambuf* old_errors = std::cerr.rdbuf(errors.rdbuf());
  try
  {
    bool success = verify();
    std::cout.rdbuf(old_output); std::cerr.rdbuf(old_errors);
    return { success, output.str(), errors.str() };
  }
  catch(...)
  {
    std::cout.rdbuf(old_output); std::cerr.rdbuf(old_errors); throw;
  }
}

bool
contains(const std::string& text, const std::string& pattern)
{
  return text.find(pattern) != std::string::npos;
}

} // namespace

int
main()
{
  omp_set_num_threads(1);
  Verbosity::set(Verbosity::SILENT);

  const char* temporary = std::getenv("TMPDIR");
  std::string root_template = (temporary == nullptr ? "/tmp" : temporary);
  root_template += "/gcsa-external-verify-XXXXXX";
  std::vector<char> root_buffer(root_template.begin(), root_template.end());
  root_buffer.push_back(0);
  char* root = mkdtemp(root_buffer.data());
  require(root != nullptr); TempFile::setDirectory(root);
  const std::string input_name = std::string(root) + "/repeated.gcsa2";
  const size_type byte_budget = 20 * KILOBYTE;

  {
    writeFixture(input_name, 200);
    ConstructionParameters parameters;
    parameters.setSteps(2); parameters.setMemoryLimitBytes(MEGABYTE);
    InputGraph graph({ input_name }, false, parameters);
    GCSA index(graph, parameters); LCPArray lcp(graph, parameters);
    const size_type retained_files = directoryEntries(root);

    VerifyResult ordinary = captureVerification([&]()
    {
      return verifyIndex(index, &lcp, graph);
    });
    VerifyResult bounded = captureVerification([&]()
    {
      return verifyIndex(index, &lcp, graph, byte_budget);
    });
    require(ordinary.success && bounded.success);
    require(ordinary.errors.empty() && bounded.errors.empty());
    require(contains(ordinary.output, "Queried the index with 3 patterns"));
    require(contains(bounded.output, "Queried the index with 3 patterns"));
    require(contains(ordinary.output, "Index verification complete"));
    require(contains(bounded.output, "Index verification complete"));
    require(directoryEntries(root) == retained_files);

    // Preserve the release verifier's failure and detailed occurrence diagnostics.
    writeFixture(input_name, 200, true);
    ordinary = captureVerification([&]()
    {
      return verifyIndex(index, &lcp, graph);
    });
    bounded = captureVerification([&]()
    {
      return verifyIndex(index, &lcp, graph, byte_budget);
    });
    require(!ordinary.success && !bounded.success);
    require(ordinary.errors == bounded.errors);
    require(contains(bounded.errors, "Expected 4:0, got 1:0"));
    require(contains(bounded.errors, "verifyIndex(): Expected { 4:0 }"));
    require(contains(bounded.errors, "verifyIndex(): Got { 1:0 }"));
    require(directoryEntries(root) == retained_files);

    bool small_budget = false;
    try { verifyIndex(index, &lcp, graph, 1); }
    catch(const std::invalid_argument&) { small_budget = true; }
    require(small_budget);
    require(directoryEntries(root) == retained_files);

    // The release verifier accepts this whole-label subset. The bounded route
    // additionally detects that an InputGraph source changed after validation.
    writeFixture(input_name, 1);
    ordinary = captureVerification([&]()
    {
      return verifyIndex(index, &lcp, graph);
    });
    require(ordinary.success);
    bool truncated = false;
    try { verifyIndex(index, &lcp, graph, byte_budget); }
    catch(const std::runtime_error& error)
    {
      truncated = contains(error.what(), "input record count changed");
    }
    require(truncated);
    require(directoryEntries(root) == retained_files);
  }

  std::remove(input_name.c_str());
  TempFile::setDirectory(TempFile::DEFAULT_TEMP_DIR);
  require(rmdir(root) == 0);
  return 0;
}
