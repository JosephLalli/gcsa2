#include <gcsa/lcp.h>
#include <gcsa/workspace.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

using namespace gcsa;

namespace
{

void requireAt(bool value, size_type line)
{
  if(!value)
  {
    std::fprintf(stderr, "test_lcp_streaming: requirement failed at line %llu\n",
      static_cast<unsigned long long>(line));
    std::abort();
  }
}

#define require(value) requireAt((value), __LINE__)

std::vector<char>
readFile(const std::string& filename)
{
  std::ifstream input(filename.c_str(), std::ios_base::binary);
  require(static_cast<bool>(input));
  return std::vector<char>(std::istreambuf_iterator<char>(input),
    std::istreambuf_iterator<char>());
}

void store(const LCPArray& lcp, const std::string& filename)
{
  require(sdsl::store_to_file(lcp, filename));
}

void writeLeaves(const std::string& filename, const std::vector<std::uint8_t>& leaves)
{
  std::ofstream output(filename.c_str(), std::ios_base::binary);
  require(static_cast<bool>(output));
  output.write(reinterpret_cast<const char*>(leaves.data()), leaves.size());
  require(static_cast<bool>(output));
}

std::vector<std::vector<std::uint8_t>>
expectedLevels(const std::vector<std::uint8_t>& leaves, size_type branching)
{
  std::vector<std::vector<std::uint8_t>> result(1, leaves);
  while(result.back().size() > 1)
  {
    const std::vector<std::uint8_t>& source = result.back();
    std::vector<std::uint8_t> target;
    for(size_type begin = 0; begin < source.size(); begin += branching)
    {
      std::uint8_t minimum = source[begin];
      for(size_type i = begin + 1; i < std::min(begin + branching,
          static_cast<size_type>(source.size())); i++)
      {
        minimum = std::min(minimum, source[i]);
      }
      target.push_back(minimum);
    }
    result.push_back(target);
  }
  return result;
}

LCPArray
legacyLCP(const std::string& leaf_filename, size_type branching)
{
  ConstructionParameters parameters;
  parameters.setLCPBranching(branching);
  InputGraph graph({ "tests/cycle.gcsa2" }, false, parameters);
  graph.lcp_name = leaf_filename;
  LCPArray result(graph, parameters);
  // The test owns the plain leaf file; do not let InputGraph remove it.
  graph.lcp_name.clear();
  return result;
}

LCPArray
externalLCP(const std::string& leaf_filename, const std::string& work_directory,
  size_type branching, bool resume = false)
{
  ConstructionParameters parameters;
  parameters.setLCPBranching(branching);
  parameters.setWorkDirectory(work_directory);
  parameters.setMemoryLimitBytes(2);
  parameters.setIOBufferSize(2);
  if(resume) { parameters.setResume(); }
  InputGraph graph({ "tests/cycle.gcsa2" }, false, parameters);
  graph.lcp_name = leaf_filename;
  LCPArray result(graph, parameters);
  graph.lcp_name.clear();
  return result;
}

} // namespace

int
main()
{
  omp_set_num_threads(1);
  Verbosity::set(Verbosity::SILENT);
  char root[] = "/tmp/gcsa-lcp-streaming-XXXXXX";
  require(mkdtemp(root) != nullptr);
  TempFile::setDirectory(root);
  const std::string leaf = std::string(root) + "/merged-leaf.lcp";
  const size_type branching = 4, byte_budget = 2;

  // 65 leaves force five levels: 65, 17, 5, 2, 1. The final one-leaf group
  // makes the tail minimum observable at every level above it.
  std::vector<std::uint8_t> leaves(65);
  for(size_type i = 0; i < leaves.size(); i++)
  {
    leaves[i] = static_cast<std::uint8_t>((19 * i + 7) % 31);
  }
  leaves[0] = 31; leaves.back() = 0;
  writeLeaves(leaf, leaves);

  LCPArray legacy = legacyLCP(leaf, branching);
  std::vector<std::vector<std::uint8_t>> expected = expectedLevels(leaves, branching);
  require(expected.size() == 5 && expected[1].back() == leaves.back());

  BuildWorkspace::Settings semantic;
  semantic["fixture"] = "lcp-streaming";
  LCPStreamingStats streaming_stats;
  LCPArray streamed;
  {
    BuildWorkspace workspace(root, semantic, BuildWorkspace::Settings(),
      BuildWorkspace::NEW_WORKSPACE);
    streamed = LCPArray(leaf, branching, byte_budget, &workspace, &streaming_stats);

    require(streaming_stats.levels == expected.size());
    require(streaming_stats.generated_levels == expected.size() - 1);
    require(streaming_stats.checkpointed_levels == expected.size() - 1);
    require(streaming_stats.restored_levels == 0);
    require(streaming_stats.max_bytes_resident <= byte_budget);
    for(size_type level = 1; level < expected.size(); level++)
    {
      require(workspace.task_completed("lcp-streaming", "level-" +
        std::to_string(level)));
    }
  }

  require(streamed.header.size == leaves.size() && streamed.branching() == branching);
  require(streamed.levels() == expected.size());
  size_type offset = 0;
  for(size_type level = 0; level < expected.size(); level++)
  {
    require(streamed.offsets[level] == offset);
    for(size_type i = 0; i < expected[level].size(); i++)
    {
      require(streamed[offset + i] == expected[level][i]);
    }
    offset += expected[level].size();
  }
  require(streamed.offsets[streamed.levels()] == offset);

  const std::string legacy_file = std::string(root) + "/legacy.lcp";
  const std::string streamed_file = std::string(root) + "/streamed.lcp";
  store(legacy, legacy_file); store(streamed, streamed_file);
  require(readFile(legacy_file) == readFile(streamed_file));

  // Existing checkpoints must be restored rather than regenerated, while the
  // serialized result and all query data remain byte-identical.
  BuildWorkspace resumed_workspace(root, semantic, BuildWorkspace::Settings(),
    BuildWorkspace::RESUME);
  LCPStreamingStats resumed_stats;
  LCPArray resumed(leaf, branching, byte_budget, &resumed_workspace, &resumed_stats);
  require(resumed_stats.generated_levels == 0 && resumed_stats.checkpointed_levels == 0);
  require(resumed_stats.restored_levels == expected.size() - 1);
  require(resumed_stats.max_bytes_resident <= byte_budget);
  const std::string resumed_file = std::string(root) + "/resumed.lcp";
  store(resumed, resumed_file);
  require(readFile(streamed_file) == readFile(resumed_file));

  LCPArray loaded;
  {
    std::ifstream input(streamed_file.c_str(), std::ios_base::binary);
    require(static_cast<bool>(input)); loaded.load(input); require(static_cast<bool>(input));
  }
  require(readFile(streamed_file) == readFile(resumed_file));
  require(loaded.size() == streamed.size() && loaded.values() == streamed.values());
  for(size_type i = 0; i < streamed.values(); i++) { require(loaded[i] == streamed[i]); }

  // The same bounded stream route works without durable artifacts.
  LCPStreamingStats transient_stats;
  LCPArray transient(leaf, branching, byte_budget, nullptr, &transient_stats);
  require(transient_stats.generated_levels == expected.size() - 1);
  require(transient_stats.checkpointed_levels == 0 && transient_stats.max_bytes_resident <= byte_budget);
  const std::string transient_file = std::string(root) + "/transient.lcp";
  store(transient, transient_file);
  require(readFile(streamed_file) == readFile(transient_file));

  // The existing InputGraph constructor selects the same stream route under
  // external-memory parameters and owns a dedicated child workspace.
  LCPArray automatic = externalLCP(leaf, root, branching);
  const std::string automatic_file = std::string(root) + "/automatic.lcp";
  store(automatic, automatic_file);
  require(readFile(legacy_file) == readFile(automatic_file));
  require(std::filesystem::exists(std::string(root) + "/lcp-levels/build.json"));
  LCPArray automatic_resumed = externalLCP(leaf, root, branching, true);
  const std::string automatic_resumed_file = std::string(root) + "/automatic-resumed.lcp";
  store(automatic_resumed, automatic_resumed_file);
  require(readFile(automatic_file) == readFile(automatic_resumed_file));

  // --resume must initialize the child workspace if a previous external
  // construction finished before this LCP stage existed.
  char missing_root[] = "/tmp/gcsa-lcp-streaming-resume-XXXXXX";
  require(mkdtemp(missing_root) != nullptr);
  LCPArray fresh_resume = externalLCP(leaf, missing_root, branching, true);
  const std::string fresh_resume_file = std::string(missing_root) + "/fresh-resume.lcp";
  store(fresh_resume, fresh_resume_file);
  require(readFile(legacy_file) == readFile(fresh_resume_file));
  require(std::filesystem::exists(std::string(missing_root) + "/lcp-levels/build.json"));
  std::filesystem::remove_all(missing_root);

  // The all-zero width and the single-leaf hierarchy match the legacy
  // bit-compressed representation, including its one-bit value width.
  const std::string zero_leaf = std::string(root) + "/zero-leaf.lcp";
  std::vector<std::uint8_t> zero_values(17, 0);
  writeLeaves(zero_leaf, zero_values);
  LCPArray zero_legacy = legacyLCP(zero_leaf, branching);
  LCPStreamingStats zero_stats;
  LCPArray zero_streamed(zero_leaf, branching, byte_budget, nullptr, &zero_stats);
  const std::string zero_legacy_file = std::string(root) + "/zero-legacy.lcp";
  const std::string zero_streamed_file = std::string(root) + "/zero-streamed.lcp";
  store(zero_legacy, zero_legacy_file); store(zero_streamed, zero_streamed_file);
  require(zero_streamed.data.width() == zero_legacy.data.width());
  require(readFile(zero_legacy_file) == readFile(zero_streamed_file));
  require(zero_stats.levels > 1 && zero_stats.max_bytes_resident <= byte_budget);

  const std::string one_leaf = std::string(root) + "/one-leaf.lcp";
  writeLeaves(one_leaf, std::vector<std::uint8_t>(1, 0));
  LCPArray one_legacy = legacyLCP(one_leaf, branching);
  LCPStreamingStats one_stats;
  LCPArray one_streamed(one_leaf, branching, byte_budget, nullptr, &one_stats);
  const std::string one_legacy_file = std::string(root) + "/one-legacy.lcp";
  const std::string one_streamed_file = std::string(root) + "/one-streamed.lcp";
  store(one_legacy, one_legacy_file); store(one_streamed, one_streamed_file);
  require(one_streamed.levels() == 1 && one_streamed.offsets[0] == 0 &&
    one_streamed.offsets[1] == 1 && one_stats.generated_levels == 0);
  require(one_streamed.data.width() == one_legacy.data.width());
  require(readFile(one_legacy_file) == readFile(one_streamed_file));

  bool rejected = false;
  try { LCPArray invalid(leaf, branching, byte_budget - 1); }
  catch(const std::invalid_argument&) { rejected = true; }
  require(rejected);

  std::filesystem::remove_all(root);
  return 0;
}
