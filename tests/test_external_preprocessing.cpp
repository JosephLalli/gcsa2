#include <gcsa/external_preprocessing.h>

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <set>
#include <string>
#include <unistd.h>
#include <vector>

using namespace gcsa;

namespace
{

void require(bool value) { if(!value) { std::abort(); } }

KMer makeKMer(size_type label, size_type from, size_type to)
{
  // k=4 consumes the low 12 label bits. All records have real predecessor
  // and successor edges, while repeated labels exercise key merging.
  key_type key = (label << 16) | (static_cast<key_type>(1) << 8) | 2;
  return KMer(key, Node::encode(from, 0), Node::encode(to, 0));
}

void writeGraph(const std::string& name, const std::vector<KMer>& kmers)
{
  std::vector<KMer> copy = kmers;
  std::ofstream output(name.c_str(), std::ios_base::binary);
  require(static_cast<bool>(output));
  writeBinary(output, copy, 4); output.close();
}

ConstructionParameters parameters(const std::string& root)
{
  ConstructionParameters result;
  result.setWorkDirectory(root);
  result.setMemoryLimitBytes(128 * KILOBYTE);
  return result;
}

} // namespace

int main()
{
  omp_set_num_threads(1);
  Verbosity::set(Verbosity::SILENT);
  const char* configured = std::getenv("TMPDIR");
  std::string pattern = std::string(
    configured != nullptr && configured[0] != '\0' ? configured : "/tmp") +
    "/gcsa-external-preprocess-XXXXXX";
  std::vector<char> root(pattern.begin(), pattern.end()); root.push_back('\0');
  require(mkdtemp(root.data()) != nullptr);
  TempFile::setDirectory(root.data());
  const std::string first = std::string(root.data()) + "/first.graph";
  const std::string second = std::string(root.data()) + "/second.graph";
  std::vector<KMer> left, right;
  std::set<size_type> labels;
  for(size_type i = 0; i < 100; i++)
  {
    size_type label = 100 + (97 * i) % 61;
    left.push_back(makeKMer(label, 10 + i, 1000 + i)); labels.insert(label);
    right.push_back(makeKMer(label, 10000 + i, 20000 + i)); labels.insert(label);
  }
  writeGraph(first, left); writeGraph(second, right);

  ConstructionParameters build_parameters = parameters(root.data());
  InputGraph graph({ first, second }, true, build_parameters);
  {
    ExternalInputPreprocessor preprocessor(graph, build_parameters);
    DeBruijnGraph mapper;
    LCP lcp;
    sdsl::int_vector<0> last_char;
    // Production builds these supports at separate lifecycle boundaries. Keep
    // the unit fixture on that route so delaying them cannot silently drift
    // from the compatibility wrapper.
    preprocessor.buildLCP(lcp);
    require(lcp.total_keys == labels.size());
    require(mapper.size() == 0);
    preprocessor.buildMapper(mapper);
    preprocessor.buildLastCharacters(last_char);
    require(mapper.size() == labels.size());
    require(last_char.size() == labels.size());
    sdsl::sd_vector<> starts;
    preprocessor.buildStartNodes(starts);
    sdsl::sd_vector<>::rank_1_type start_rank;
    sdsl::util::init_support(start_rank, &starts);
    require(start_rank(starts.size()) == left.size() + right.size());

    PathGraph initial(0, graph.k(), 0);
    preprocessor.buildInitialPathGraph(initial);
    require(initial.size() == left.size() + right.size());
    require(initial.ranks() == 2 * initial.size());
    require(initial.files() >= graph.files());
    bool saw_first = false, saw_second = false;
    for(size_type file = 0; file < initial.files(); file++)
    {
      logical_file_id_t logical = initial.logicalFile(file);
      require(logical == logical_file_id_t(0) || logical == logical_file_id_t(1));
      require(initial.physicalShard(file).value >> 32 == logical.value);
      saw_first = saw_first || logical == logical_file_id_t(0);
      saw_second = saw_second || logical == logical_file_id_t(1);
      std::vector<PathNode> paths;
      std::vector<PathNode::rank_type> ranks;
      initial.read(paths, ranks, file);
      for(size_type i = 1; i < paths.size(); i++)
      {
        require(ranks[paths[i - 1].pointer()] <= ranks[paths[i].pointer()]);
      }
    }
    require(saw_first && saw_second);
  }

  // Text scanning uses the same bounded callback contract, including a line
  // whose two destinations must be split across one-record blocks.
  const std::string text = std::string(root.data()) + "/input.gcsa2";
  {
    std::ofstream output(text.c_str());
    output << "ACGT\t1:0\tA\tC\t2:0,3:0\n";
    output << "CGTA\t2:0\tC\tG\t4:0\n";
  }
  InputGraph text_graph({ text }, false, build_parameters);
  size_type text_records = 0, text_blocks = 0;
  text_graph.scanKMerBlocks(0, sizeof(KMer),
    [&text_records, &text_blocks](size_type, const std::vector<KMer>& block)
    {
      require(block.size() == 1); text_records += block.size(); text_blocks++;
    });
  require(text_records == 3 && text_blocks == 3);

  std::filesystem::remove_all(root.data());
  return 0;
}
