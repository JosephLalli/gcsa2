#include <gcsa/checkpoint.h>
#include <gcsa/compressed_block.h>
#include <gcsa/path_graph_external.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <fcntl.h>
#include <iterator>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>

using namespace gcsa;

static void require(bool condition)
{
  if(!condition) { std::fprintf(stderr, "checkpoint checksum test failure\n"); std::abort(); }
}

static std::vector<char> bytes(const std::string& path)
{
  std::ifstream input(path.c_str(), std::ios::binary);
  require(bool(input));
  return std::vector<char>(std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>());
}

static void verify(const ClosedPayloadChecksum& provenance, const std::string& path)
{
  const auto payload = bytes(path);
  require(provenance.matches(path, payload.size()));
  require(provenance.value == BuildWorkspace::checksum(payload.data(), payload.size()));
}

static void writer_test(const std::string& root, CompressedBlockWriter::Mode mode,
  size_type workers, size_type count)
{
  const std::string path = root + "/writer";
  CompressedBlockWriter writer(path, 4 * 1024 * 1024, mode, 1, workers);
  bool rejected = false;
  try { writer.storedChecksum(); } catch(const std::logic_error&) { rejected = true; }
  require(rejected);
  writer.trackStoredChecksum();
  std::vector<uint64_t> record(16384);
  for(size_type n = 0; n < count; ++n)
  {
    for(size_type j = 0; j < record.size(); ++j) { record[j] = n * 197 + j % 127; }
    writer.writeRecord(record.data(), record.size() * sizeof(uint64_t));
  }
  writer.finish();
  const auto payload = bytes(path);
  require(writer.storedChecksum() == BuildWorkspace::checksum(payload.data(), payload.size()));
  require(unlink(path.c_str()) == 0);
}

static void checkpoint_test(const std::string& root, TempCompression compression,
  size_type count, const std::string& task)
{
  PathGraph graph(1, 16, 0);
  size_type committed = 0;
  ExternalPathSortSink sink(graph, 0, 16 * MEGABYTE, 4, 256 * MEGABYTE,
    committed, nullptr, TempFileCodecParameters(compression, 64 * 1024, 2, 1));
  for(size_type i = 0; i < count; ++i)
  {
    PathNode node;
    node.from = count - i; node.to = count + i; node.fields = 0;
    node.setOrder(1); node.setLCP(1); node.setPredecessors(2); node.setPointer(0);
    PathNode::rank_type ranks[2] = { static_cast<PathNode::rank_type>(i % 97), 0 };
    sink.write(node, ranks);
  }
  sink.finish();
  verify(graph.path_checksums[0], graph.path_names[0]);
  verify(graph.rank_checksums[0], graph.rank_names[0]);
  const auto paths = bytes(graph.path_names[0]), ranks = bytes(graph.rank_names[0]);
  const uint64_t total = paths.size() + ranks.size();
  BuildWorkspace workspace(root + "/" + task, BuildWorkspace::Settings(),
    BuildWorkspace::Settings(), BuildWorkspace::NEW_WORKSPACE);
  uint64_t scanned = BuildWorkspace::adoption_checksum_scan_bytes.load();
  uint64_t reused = BuildWorkspace::adoption_checksum_reused_bytes.load();
  checkpointPathGraph(workspace, graph, "first", "write", 7777);
  require(BuildWorkspace::adoption_checksum_scan_bytes.load() == scanned);
  require(BuildWorkspace::adoption_checksum_reused_bytes.load() - reused == total);
  PathGraph restored(0, 0, 0);
  restorePathGraph(workspace, restored, "first", "write", 7777, true);
  require(bytes(restored.path_names[0]) == paths && bytes(restored.rank_names[0]) == ranks);
  require(!restored.path_checksums[0].valid && !restored.rank_checksums[0].valid);
  scanned = BuildWorkspace::adoption_checksum_scan_bytes.load();
  reused = BuildWorkspace::adoption_checksum_reused_bytes.load();
  checkpointPathGraph(workspace, restored, "second", "write", 7777);
  require(BuildWorkspace::adoption_checksum_scan_bytes.load() - scanned == total);
  require(BuildWorkspace::adoption_checksum_reused_bytes.load() == reused);
}

static void stale_test(const std::string& root)
{
  PathGraph graph(1, 16, 0);
  size_type committed = 0;
  ExternalPathSortSink sink(graph, 0, 4 * MEGABYTE, 2, 64 * MEGABYTE, committed);
  PathNode node; node.from = 1; node.to = 2; node.fields = 0;
  node.setOrder(1); node.setLCP(1); node.setPredecessors(1); node.setPointer(0);
  PathNode::rank_type ranks[2] = { 1, 0 };
  sink.write(node, ranks); sink.finish();
  const auto old = bytes(graph.path_names[0]);
  const auto signature = graph.path_checksums[0];
  const std::string replacement = root + "/replacement";
  { std::ofstream output(replacement.c_str(), std::ios::binary);
    node.from = 3; output.write(reinterpret_cast<const char*>(&node), sizeof(node)); }
  require(rename(replacement.c_str(), graph.path_names[0].c_str()) == 0);
  // Even a same-size replacement with restored mtime cannot reuse the digest.
  struct timespec times[2] = { signature.identity.st_atim, signature.identity.st_mtim };
  require(utimensat(AT_FDCWD, graph.path_names[0].c_str(), times, 0) == 0);
  require(!signature.matches(graph.path_names[0], old.size()));
  const uint64_t scanned = BuildWorkspace::adoption_checksum_scan_bytes.load();
  BuildWorkspace workspace(root + "/stale", BuildWorkspace::Settings(),
    BuildWorkspace::Settings(), BuildWorkspace::NEW_WORKSPACE);
  checkpointPathGraph(workspace, graph, "stale", "write", 4096);
  require(BuildWorkspace::adoption_checksum_scan_bytes.load() - scanned == old.size());
  PathGraph restored(0, 0, 0);
  restorePathGraph(workspace, restored, "stale", "write", 4096, true);
  require(bytes(restored.path_names[0]) == bytes(graph.path_names[0]));
  require(bytes(restored.path_names[0]) != old);
}

int main()
{
  const char* tmp = std::getenv("TMPDIR");
  std::string pattern = std::string(tmp ? tmp : "/tmp") + "/gcsa-checksums-XXXXXX";
  std::vector<char> buffer(pattern.begin(), pattern.end()); buffer.push_back(0);
  require(mkdtemp(buffer.data()) != nullptr);
  const std::string root(buffer.data()); TempFile::setDirectory(root);
  writer_test(root, CompressedBlockWriter::RAW, 1, 0);
  writer_test(root, CompressedBlockWriter::RAW, 1, 67);
  writer_test(root, CompressedBlockWriter::ZSTD, 1, 67);
  writer_test(root, CompressedBlockWriter::ZSTD, 4, 67);
  checkpoint_test(root, TempCompression::NONE, 0, "empty-raw");
  checkpoint_test(root, TempCompression::NONE, 4097, "raw");
  checkpoint_test(root, TempCompression::AUTO, 0, "empty-framed");
  checkpoint_test(root, TempCompression::AUTO, 4097, "framed");
  stale_test(root);
  std::fprintf(stderr, "checkpoint checksums: raw/framed reuse, verified resume, unknown/stale fallback passed\n");
}
