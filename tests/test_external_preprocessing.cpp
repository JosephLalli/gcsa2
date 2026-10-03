#include <gcsa/external_preprocessing.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <set>
#include <string>
#include <sys/stat.h>
#include <sys/wait.h>
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
  result.setSortRunSize(ExternalFixedRecordSorter::minimumBudget(
    sizeof(KMer) + sizeof(std::uint64_t)));
  result.setIOBufferSize(KILOBYTE);
  result.setMergeFanIn(2);
  // This deliberately small physical-shard interval is independent of logical
  // input identity and forces a many-shard initial PathGraph.
  result.setCheckpointRecords(5);
  return result;
}

std::string readBytes(const std::string& name)
{
  std::ifstream input(name.c_str(), std::ios_base::binary);
  require(static_cast<bool>(input));
  return std::string(std::istreambuf_iterator<char>(input),
    std::istreambuf_iterator<char>());
}

// The reduced stream artifacts, named as ExternalInputPreprocessor names them.
const ArtifactIdentity KEY_ARTIFACT("preprocess", "key-start-streams",
  "unique-keys", "key-stream-v1");
const ArtifactIdentity START_ARTIFACT("preprocess", "key-start-streams",
  "unique-start-nodes", "node-stream-v1");
const std::string PREPROCESS_MARKER = "preprocess--key-start-streams.complete";

} // namespace

int main()
{
  omp_set_num_threads(1);
  Verbosity::set(Verbosity::SILENT);
  char root[] = "/tmp/gcsa-external-preprocess-XXXXXX";
  require(mkdtemp(root) != nullptr);
  TempFile::setDirectory(root);
  const std::string first = std::string(root) + "/first.graph";
  const std::string second = std::string(root) + "/second.graph";
  std::vector<KMer> left, right;
  std::set<size_type> labels;
  for(size_type i = 0; i < 100; i++)
  {
    size_type label = 100 + (97 * i) % 61;
    left.push_back(makeKMer(label, 10 + i, 1000 + i)); labels.insert(label);
    right.push_back(makeKMer(label, 10000 + i, 20000 + i)); labels.insert(label);
  }
  writeGraph(first, left); writeGraph(second, right);

  ConstructionParameters build_parameters = parameters(root);
  InputGraph graph({ first, second }, true, build_parameters);
  BuildWorkspace::Settings semantic;
  semantic["fixture"] = "external-preprocessing";

  // The reduced streams are copied into their artifacts while the sorter
  // writes them, and published only after both reductions. Crash right after
  // the key artifact's rename, before the task marker exists: recovery must
  // drop that renamed artifact and the start artifact's partial file, and a
  // rerun must publish the same artifacts under the same marker as the
  // uninterrupted build below (compared at the end). The child is forked
  // before this process starts any helper thread.
  const std::string crash_root = std::string(root) + "/crash";
  require(mkdir(crash_root.c_str(), 0755) == 0);
  {
    BuildWorkspace created(crash_root, semantic, BuildWorkspace::Settings(),
      BuildWorkspace::NEW_WORKSPACE);
  }
  pid_t child = fork();
  require(child >= 0);
  if(child == 0)
  {
    setenv("GCSA_WORKSPACE_CRASH_POINT", "artifact-after-rename", 1);
    try
    {
      BuildWorkspace crashing(crash_root, semantic);
      ExternalInputPreprocessor preprocessor(graph, build_parameters, &crashing);
      preprocessor.keyCount();
    }
    catch(...) { ::_exit(88); }
    ::_exit(87); // The requested boundary was not reached.
  }
  int status = 0;
  require(waitpid(child, &status, 0) == child);
  require(WIFEXITED(status) && WEXITSTATUS(status) == 86);
  BuildWorkspace workspace(root, semantic, BuildWorkspace::Settings(),
    BuildWorkspace::NEW_WORKSPACE);

  {
    ExternalInputPreprocessor preprocessor(graph, build_parameters, &workspace);
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
    require(initial.files() > graph.files());
    bool saw_first = false, saw_second = false;
    for(size_type file = 0; file < initial.files(); file++)
    {
      logical_file_id_t logical = initial.logicalFile(file);
      require(logical == logical_file_id_t(0) || logical == logical_file_id_t(1));
      require(initial.physicalShard(file).value >> 32 == logical.value);
      saw_first = saw_first || logical == logical_file_id_t(0);
      saw_second = saw_second || logical == logical_file_id_t(1);
      // Checkpoint adoption trusts the digest each shard's writer computed on
      // its helper thread, so it must be the checksum of the file as closed.
      const std::string path_bytes = readBytes(initial.path_names[file]);
      const std::string rank_bytes = readBytes(initial.rank_names[file]);
      require(initial.path_checksums[file].matches(initial.path_names[file], path_bytes.size()));
      require(initial.rank_checksums[file].matches(initial.rank_names[file], rank_bytes.size()));
      require(initial.path_checksums[file].value ==
        BuildWorkspace::checksum(path_bytes.data(), path_bytes.size()));
      require(initial.rank_checksums[file].value ==
        BuildWorkspace::checksum(rank_bytes.data(), rank_bytes.size()));
      std::vector<PathNode> paths;
      std::vector<PathNode::rank_type> ranks;
      initial.read(paths, ranks, file);
      for(size_type i = 1; i < paths.size(); i++)
      {
        require(ranks[paths[i - 1].pointer()] <= ranks[paths[i].pointer()]);
      }
    }
    require(saw_first && saw_second);
    require(preprocessor.stats().key_sort.runs > 2);
    require(preprocessor.stats().start_sort.runs > 2);
    require(preprocessor.stats().key_sort.merge_operations > 0);
    require(preprocessor.stats().key_sort.max_bytes_resident <=
      build_parameters.getSortRunSize());
    require(preprocessor.stats().physical_shards == initial.files());
  }
  require(workspace.task_completed("preprocess", "key-start-streams"));

  // Both artifacts pass full validation (header, checksum over the payload,
  // footer, and the task marker's record), and their payloads are exactly the
  // reduced streams: one merged key per label in label order, and the sorted
  // distinct start nodes.
  workspace.validate_artifact(KEY_ARTIFACT, logical_file_id_t(0), physical_shard_id_t(0));
  workspace.validate_artifact(START_ARTIFACT, logical_file_id_t(0), physical_shard_id_t(1));
  {
    std::vector<key_type> keys;
    for(size_type label : labels) { keys.push_back(makeKMer(label, 0, 0).key); }
    std::vector<std::uint8_t> payload = workspace.read_artifact_payload(KEY_ARTIFACT,
      logical_file_id_t(0), physical_shard_id_t(0), MEGABYTE);
    require(payload.size() == keys.size() * sizeof(key_type));
    require(std::memcmp(payload.data(), keys.data(), payload.size()) == 0);
    std::set<node_type> nodes;
    for(const KMer& kmer : left) { nodes.insert(kmer.from); }
    for(const KMer& kmer : right) { nodes.insert(kmer.from); }
    std::vector<node_type> starts(nodes.begin(), nodes.end());
    payload = workspace.read_artifact_payload(START_ARTIFACT,
      logical_file_id_t(0), physical_shard_id_t(1), MEGABYTE);
    require(payload.size() == starts.size() * sizeof(node_type));
    require(std::memcmp(payload.data(), starts.data(), payload.size()) == 0);
  }

  // Resume must restore the two reduced streams instead of rerunning their
  // forced spill/merge phases. (Initial KMer shard materialization is covered
  // separately by checkpointPathGraph in GCSA construction.)
  {
    ExternalInputPreprocessor resumed(graph, build_parameters, &workspace);
    require(resumed.keyCount() == labels.size());
    require(resumed.startNodeCount() == left.size() + right.size());
    require(resumed.stats().key_sort.runs == 0);
    require(resumed.stats().start_sort.runs == 0);
  }

  // Text scanning uses the same bounded callback contract, including a line
  // whose two destinations must be split across one-record blocks.
  const std::string text = std::string(root) + "/input.gcsa2";
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

  // Recovery from the crash above, then the rerun.
  {
    BuildWorkspace recovered(crash_root, semantic);
    require(!recovered.task_completed("preprocess", "key-start-streams"));
    require(access(recovered.artifact_path(KEY_ARTIFACT, logical_file_id_t(0),
      physical_shard_id_t(0)).c_str(), F_OK) != 0);
    for(const auto& entry : std::filesystem::directory_iterator(crash_root))
    {
      require(entry.path().filename().string().find(".partial") == std::string::npos);
    }
    ExternalInputPreprocessor rerun(graph, build_parameters, &recovered);
    require(rerun.keyCount() == labels.size());
    require(rerun.startNodeCount() == left.size() + right.size());
    require(recovered.task_completed("preprocess", "key-start-streams"));
  }
  require(readBytes(crash_root + "/" + PREPROCESS_MARKER) ==
    readBytes(std::string(root) + "/" + PREPROCESS_MARKER));

  std::filesystem::remove_all(root);
  return 0;
}
