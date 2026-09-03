#include <gcsa/algorithms.h>
#include <gcsa/files.h>
#include <gcsa/gcsa.h>
#include <gcsa/lcp.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <unistd.h>
#include <vector>

using namespace gcsa;

namespace
{

void require(bool value) { if(!value) { std::abort(); } }

} // namespace

int
main()
{
  omp_set_num_threads(1);
  Verbosity::set(Verbosity::SILENT);

  char root[] = "/tmp/gcsa-external-verify-XXXXXX";
  require(mkdtemp(root) != nullptr); TempFile::setDirectory(root);
  const std::string input_name = std::string(root) + "/repeated.gcsa2";
  {
    std::ofstream output(input_name.c_str()); require(static_cast<bool>(output));
    // Each label has a deliberately much larger input group than the forced
    // external sort run, while preserving the small, valid cycle graph.
    for(size_type i = 0; i < 200; i++)
    {
      output << "AC\t1:0\tG\tG\t2:0\n";
      output << "CG\t2:0\tA\tA\t3:0\n";
      output << "GA\t3:0\tC\tC\t1:0\n";
    }
  }

  ConstructionParameters parameters;
  parameters.setSteps(2); parameters.setMemoryLimitBytes(MEGABYTE);
  InputGraph graph({ input_name }, false, parameters);
  GCSA index(graph, parameters); LCPArray lcp(graph, parameters);

  // Preserve the established in-memory verifier as the equivalence reference.
  std::vector<KMer> kmers; graph.read(kmers);
  require(verifyIndex(index, &lcp, kmers, graph.k(), graph.mapping));

  const size_type budget = verifyIndexMinimumBudget();
  ExternalFixedRecordSortStats statistics;
  require(verifyIndex(index, &lcp, graph, budget, 2, &statistics));
  require(statistics.runs > 2);
  require(statistics.max_bytes_resident <= budget);

  std::remove(input_name.c_str()); rmdir(root);
  return 0;
}
