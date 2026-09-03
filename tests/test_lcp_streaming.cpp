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

LCPArray load(const std::string& filename)
{
  LCPArray result;
  std::ifstream input(filename.c_str(), std::ios_base::binary);
  require(static_cast<bool>(input)); result.load(input); require(static_cast<bool>(input));
  return result;
}

void compareArrays(const LCPArray& left, const LCPArray& right)
{
  require(left.size() == right.size());
  require(left.values() == right.values());
  require(left.levels() == right.levels());
  for(size_type i = 0; i < left.values(); i++) { require(left[i] == right[i]); }
  for(size_type i = 0; i <= left.levels(); i++)
  {
    require(left.offsets[i] == right.offsets[i]);
  }
}

void compareQueries(const LCPArray& left, const LCPArray& right)
{
  compareArrays(left, right);
  for(size_type i = 0; i < left.size(); i++)
  {
    require(left.psv(i) == right.psv(i));
    require(left.psev(i) == right.psev(i));
    require(left.nsv(i) == right.nsv(i));
    require(left.nsev(i) == right.nsev(i));
  }
  for(size_type begin = 0; begin < left.size(); begin++)
  {
    for(size_type end = begin; end < left.size(); end++)
    {
      require(left.rmq(begin, end) == right.rmq(begin, end));
    }
  }
}

void writeLeaves(const std::string& filename, const std::vector<std::uint8_t>& leaves)
{
  std::ofstream output(filename.c_str(), std::ios_base::binary);
  require(static_cast<bool>(output));
  if(!leaves.empty())
  {
    output.write(reinterpret_cast<const char*>(leaves.data()), leaves.size());
  }
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

void
externalStore(const std::string& leaf_filename, const std::string& work_directory,
  const std::string& output_filename, size_type branching, bool resume = false)
{
  ConstructionParameters parameters;
  parameters.setLCPBranching(branching);
  parameters.setWorkDirectory(work_directory);
  parameters.setMemoryLimitBytes(2);
  parameters.setIOBufferSize(2);
  if(resume) { parameters.setResume(); }
  InputGraph graph({ "tests/cycle.gcsa2" }, false, parameters);
  graph.lcp_name = leaf_filename;
  LCPArray::buildAndStore(graph, parameters, output_filename);
  graph.lcp_name.clear();
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
  const std::string direct_file = std::string(root) + "/direct.lcp";
  store(legacy, legacy_file); store(streamed, streamed_file);
  require(readFile(legacy_file) == readFile(streamed_file));
  LCPStreamingStats direct_stats;
  LCPArray::buildAndStore(leaf, branching, byte_budget, direct_file,
    nullptr, &direct_stats);
  require(readFile(legacy_file) == readFile(direct_file));
  require(direct_stats.levels == expected.size());
  require(direct_stats.generated_levels == expected.size() - 1);
  require(direct_stats.max_bytes_resident <= byte_budget);
  LCPArray direct = load(direct_file);
  compareQueries(legacy, direct);

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
  const std::string direct_resumed_file = std::string(root) + "/direct-resumed.lcp";
  LCPStreamingStats direct_resumed_stats;
  LCPArray::buildAndStore(leaf, branching, byte_budget, direct_resumed_file,
    &resumed_workspace, &direct_resumed_stats);
  require(direct_resumed_stats.generated_levels == 0);
  require(direct_resumed_stats.restored_levels == expected.size() - 1);
  require(direct_resumed_stats.max_bytes_resident <= byte_budget);
  require(readFile(streamed_file) == readFile(direct_resumed_file));

  LCPArray loaded = load(streamed_file);
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

  // The public InputGraph overload owns the same resumable child workspace
  // while publishing directly, rather than returning a resident hierarchy.
  char direct_root[] = "/tmp/gcsa-lcp-direct-XXXXXX";
  require(mkdtemp(direct_root) != nullptr);
  const std::string automatic_direct_file = std::string(direct_root) + "/automatic.lcp";
  externalStore(leaf, direct_root, automatic_direct_file, branching);
  require(readFile(legacy_file) == readFile(automatic_direct_file));
  require(std::filesystem::exists(std::string(direct_root) + "/lcp-levels/build.json"));
  const std::string automatic_direct_resumed = std::string(direct_root) + "/resumed.lcp";
  externalStore(leaf, direct_root, automatic_direct_resumed, branching, true);
  require(readFile(automatic_direct_file) == readFile(automatic_direct_resumed));
  std::filesystem::remove_all(direct_root);

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
  const std::string zero_direct_file = std::string(root) + "/zero-direct.lcp";
  LCPStreamingStats zero_direct_stats;
  LCPArray::buildAndStore(zero_leaf, branching, byte_budget, zero_direct_file,
    nullptr, &zero_direct_stats);
  require(readFile(zero_legacy_file) == readFile(zero_direct_file));
  compareQueries(zero_legacy, load(zero_direct_file));
  require(zero_direct_stats.max_bytes_resident <= byte_budget);

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
  const std::string one_direct_file = std::string(root) + "/one-direct.lcp";
  LCPStreamingStats one_direct_stats;
  LCPArray::buildAndStore(one_leaf, branching, byte_budget, one_direct_file,
    nullptr, &one_direct_stats);
  require(readFile(one_legacy_file) == readFile(one_direct_file));
  compareQueries(one_legacy, load(one_direct_file));
  require(one_direct_stats.max_bytes_resident <= byte_budget);

  // The full eight-bit width starts each output word with one bits, matching
  // the resident constructor's untouched final-word padding.
  const std::string full_leaf = std::string(root) + "/full-width-leaf.lcp";
  std::vector<std::uint8_t> full_values = { 255, 9, 7, 201, 3, 19, 41, 2, 88 };
  writeLeaves(full_leaf, full_values);
  LCPArray full_legacy = legacyLCP(full_leaf, branching);
  const std::string full_legacy_file = std::string(root) + "/full-width-legacy.lcp";
  const std::string full_direct_file = std::string(root) + "/full-width-direct.lcp";
  store(full_legacy, full_legacy_file);
  LCPStreamingStats full_direct_stats;
  LCPArray::buildAndStore(full_leaf, branching, byte_budget, full_direct_file,
    nullptr, &full_direct_stats);
  require(readFile(full_legacy_file) == readFile(full_direct_file));
  compareArrays(full_legacy, load(full_direct_file));
  require(full_direct_stats.max_bytes_resident <= byte_budget);

  // With a tiny seven-bit hierarchy, SDSL compresses in place without
  // reallocating and retains historical source-word padding. Direct packing
  // intentionally delegates this bounded sub-cacheline case back to SDSL.
  const std::string padding_leaf = std::string(root) + "/legacy-padding-leaf.lcp";
  std::vector<std::uint8_t> padding_values = { 127, 65, 9, 3, 88, 17, 4 };
  writeLeaves(padding_leaf, padding_values);
  LCPArray padding_legacy = legacyLCP(padding_leaf, 64);
  const std::string padding_legacy_file = std::string(root) + "/legacy-padding.lcp";
  const std::string padding_direct_file = std::string(root) + "/direct-padding.lcp";
  store(padding_legacy, padding_legacy_file);
  LCPStreamingStats padding_direct_stats;
  LCPArray::buildAndStore(padding_leaf, 64, byte_budget, padding_direct_file,
    nullptr, &padding_direct_stats);
  require(readFile(padding_legacy_file) == readFile(padding_direct_file));
  compareArrays(padding_legacy, load(padding_direct_file));
  require(padding_direct_stats.max_bytes_resident <= byte_budget);

  // Empty raw hierarchies still carry a single empty leaf level and must use
  // the caller's branching factor in the serialized header.
  const std::string empty_leaf = std::string(root) + "/empty-leaf.lcp";
  writeLeaves(empty_leaf, std::vector<std::uint8_t>());
  LCPArray empty_legacy = legacyLCP(empty_leaf, branching);
  const std::string empty_legacy_file = std::string(root) + "/empty-legacy.lcp";
  const std::string empty_direct_file = std::string(root) + "/empty-direct.lcp";
  store(empty_legacy, empty_legacy_file);
  LCPStreamingStats empty_direct_stats;
  LCPArray::buildAndStore(empty_leaf, branching, byte_budget, empty_direct_file,
    nullptr, &empty_direct_stats);
  require(readFile(empty_legacy_file) == readFile(empty_direct_file));
  compareArrays(empty_legacy, load(empty_direct_file));
  require(empty_direct_stats.levels == 1 && empty_direct_stats.max_bytes_resident <= byte_budget);

  bool rejected = false;
  try { LCPArray invalid(leaf, branching, byte_budget - 1); }
  catch(const std::invalid_argument&) { rejected = true; }
  require(rejected);

  std::filesystem::remove_all(root);
  return 0;
}
