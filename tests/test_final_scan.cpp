/*
  The pipelined final-event scan must write exactly the serial scan's events.
  Each case builds the same randomized input twice up to the committed event
  frontier, once with GCSA_SERIAL_FINAL_SCAN=1 and once pipelined, and
  requires every final-event artifact to be byte-identical. One case then
  completes both indexes and compares the public files.
*/

#include <gcsa/algorithms.h>
#include <gcsa/files.h>
#include <gcsa/gcsa.h>
#include <gcsa/lcp.h>
#include <gcsa/path_graph_external.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <omp.h>
#include <unistd.h>

using namespace gcsa;

namespace
{

void
requireAt(bool condition, size_type line)
{
  if(!condition)
  {
    std::cerr << "test_final_scan: requirement failed at line " << line << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

#define require(condition) requireAt((condition), __LINE__)

std::string
makeTempRoot(const std::string& stem)
{
  const char* configured = std::getenv("TMPDIR");
  std::filesystem::path parent =
    (configured != nullptr && configured[0] != '\0' ? configured : "/tmp");
  std::string pattern = (parent / (stem + "-XXXXXX")).string();
  std::vector<char> mutable_pattern(pattern.begin(), pattern.end());
  mutable_pattern.push_back('\0');
  require(mkdtemp(mutable_pattern.data()) != nullptr);
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

std::map<std::string, std::vector<char>>
workspaceArtifacts(const std::string& workspace, const std::string& prefix)
{
  std::map<std::string, std::vector<char>> result;
  for(const std::filesystem::directory_entry& entry :
      std::filesystem::directory_iterator(workspace))
  {
    const std::string name = entry.path().filename().string();
    if(entry.is_regular_file() && name.compare(0, prefix.size(), prefix) == 0 &&
       entry.path().extension() == ".bin")
    {
      result[name] = readFile(entry.path().string());
    }
  }
  return result;
}

/*
  A randomized input. DNA k-mers are every occurrence in one cyclic sequence
  over `chars`, so the de Bruijn structure is consistent; characters left out
  of `chars` give alphabet components without paths. A random sequence of
  `length` repeats k-mers, so paths double, merge and gather several start
  nodes; a de Bruijn sequence (debruijn, length ignored) has every k-mer once.
  Graph positions follow the sequence in runs of consecutive offsets, which
  makes chains of unsampled paths whose start nodes continue their
  predecessor's, except for shared_percent of positions drawn from a pool of
  61 shared nodes, which break chains and make multi-node paths and ranks
  recurring across batches. The N self-loop path gets heavy_copies start
  nodes; past a slot's arena of four nodes per batch path it is a serial path
  in the spillable set.

  A linear input (linear) is '#', the random sequence and '$', as vg's source
  and sink overlay reads it: k-mers running past the end are padded with
  '$', the source's predecessor is '$' and the sink's successor '#'. Paths
  following the source then have the sink as predecessor and are sampled.

  A random cycle is a graph without tips, which vg breaks with its source
  and sink before GCSA2 sees it, and its external construction stops at the
  edge component on every scan route ("invalid edge stream"); so do some of
  these single-character-sink linear inputs. Their events are still the
  scan's complete output, so such cases compare events only; the de Bruijn
  case, like test_construction_resume's fixture, builds completely.
*/
struct ScanCase
{
  std::uint64_t seed;
  std::string chars;
  bool debruijn, linear;
  size_type length, shared_percent, heavy_copies;
  size_type io_buffer;      // batch size is about io_buffer / 248 paths
  int threads;
  size_type sample_period;
  bool mapping, framed, previous_on_disk;
};

constexpr size_type K = 5;
constexpr size_type SHARED_POOL = 61, SHARED_FIRST_ID = 5;

// Appends a de Bruijn sequence of the given order over alphabet_size values.
void
appendDeBruijn(size_type t, size_type period, size_type order,
  size_type alphabet_size, std::vector<size_type>& state,
  std::vector<size_type>& sequence)
{
  if(t > order)
  {
    if(order % period == 0)
    {
      for(size_type i = 1; i <= period; i++) { sequence.push_back(state[i]); }
    }
    return;
  }
  state[t] = state[t - period];
  appendDeBruijn(t + 1, period, order, alphabet_size, state, sequence);
  for(size_type value = state[t - period] + 1; value < alphabet_size; value++)
  {
    state[t] = value;
    appendDeBruijn(t + 1, t, order, alphabet_size, state, sequence);
  }
}

void
writeInput(const std::string& filename, const ScanCase& scan)
{
  std::mt19937_64 rng(scan.seed);
  std::ofstream output(filename.c_str());
  require(static_cast<bool>(output));
  std::string sequence;
  if(scan.debruijn)
  {
    std::vector<size_type> state(scan.chars.size() * K + 1, 0), values;
    appendDeBruijn(1, 1, K, scan.chars.size(), state, values);
    for(size_type value : values) { sequence.push_back(scan.chars[value]); }
  }
  else
  {
    sequence.resize(scan.length);
    for(char& c : sequence) { c = scan.chars[rng() % scan.chars.size()]; }
    if(scan.linear) { sequence = "#" + sequence + "$"; }
  }
  const size_type length = sequence.size();
  std::vector<std::string> nodes(length);
  for(size_type i = 0; i < length; i++)
  {
    if(scan.linear && (i == 0 || i + 1 == length))
    {
      nodes[i] = (i == 0 ? "1:0" : "2:0");
    }
    else if(rng() % 100 < scan.shared_percent)
    {
      nodes[i] = std::to_string(SHARED_FIRST_ID + rng() % SHARED_POOL) + ":" +
        std::to_string(rng() % 3);
    }
    else
    {
      nodes[i] = std::to_string(100 + i / 1000) + ":" + std::to_string(i % 1000);
    }
  }
  for(size_type i = 0; i < length; i++)
  {
    std::string label;
    char successor = sequence[(i + K) % length];
    if(scan.linear)
    {
      for(size_type j = 0; j < K; j++)
      {
        label.push_back(i + j < length ? sequence[i + j] : '$');
      }
      if(i + K >= length) { successor = (i + 1 == length ? '#' : '$'); }
    }
    else
    {
      for(size_type j = 0; j < K; j++) { label.push_back(sequence[(i + j) % length]); }
    }
    output << label << '\t' << nodes[i] << '\t'
           << sequence[(i + length - 1) % length] << '\t'
           << successor << '\t' << nodes[(i + 1) % length] << '\n';
  }
  const std::string heavy_label(K, 'N');
  for(size_type copy = 0; copy < scan.heavy_copies; copy++)
  {
    const size_type id = 50000 + copy;
    output << heavy_label << '\t' << id << ":1\tN\tN\t" << id << ":1\n";
  }
  output.close(); require(static_cast<bool>(output));
}

// A non-monotone mapping of the shared pool that collapses pairs of ids, so
// start nodes are mapped before they are deduplicated.
void
writeMapping(const std::string& filename)
{
  NodeMapping mapping(SHARED_FIRST_ID);
  for(size_type j = 0; j < SHARED_POOL; j++)
  {
    mapping.insert(SHARED_FIRST_ID + ((j * 37) % SHARED_POOL) / 2);
  }
  require(sdsl::store_to_file(mapping, filename));
}

ConstructionParameters
caseParameters(const ScanCase& scan, const std::string& workspace)
{
  ConstructionParameters parameters;
  parameters.setSteps(2);
  parameters.setWorkDirectory(workspace);
  parameters.setMemoryLimitBytes(externalPathJoinMinimumBudget() +
    externalPathGraphSortMinimumBudget() + 64 * MEGABYTE);
  parameters.setSortRunSize(externalPathGraphSortMinimumBudget());
  parameters.setJoinPartitionSize(externalPathJoinMinimumBudget());
  parameters.setMergeFanIn(2);
  parameters.setMaxOpenFiles(ConstructionParameters::MIN_OPEN_FILES);
  parameters.setSamplePeriod(scan.sample_period);
  if(scan.framed)
  {
    parameters.setTempCompression("zstd");
    parameters.setCompressionBlockSize(64 * KILOBYTE);
    parameters.setCompressionWorkers(1);
    parameters.setJoinPartitionSize(4 * MEGABYTE);
    parameters.setSortRunSize(2 * MEGABYTE);
  }
  parameters.setIOBufferSize(scan.io_buffer);
  return parameters;
}

void
setSwitch(const char* name, bool on)
{
  if(on) { require(::setenv(name, "1", 1) == 0); }
  else { require(::unsetenv(name) == 0); }
}

/*
  Builds through the final events, or the whole index when `prefix` is set,
  and returns the EXTENDED log.
*/
std::string
build(const ScanCase& scan, const std::string& input, const std::string& mapping,
  const std::string& workspace, bool pipelined, const std::string& prefix = std::string())
{
  setSwitch("GCSA_SERIAL_FINAL_SCAN", !pipelined);
  setSwitch("GCSA_FINAL_SCAN_PREVIOUS_ON_DISK", scan.previous_on_disk);
  omp_set_num_threads(scan.threads);
  TempFile::setDirectory(workspace);
  ConstructionParameters parameters = caseParameters(scan, workspace);
  if(prefix.empty()) { parameters.setStopAfter("final-events"); }
  std::ostringstream log;
  std::streambuf* old_stderr = std::cerr.rdbuf(log.rdbuf());
  Verbosity::set(Verbosity::EXTENDED);
  bool stopped = false;
  try
  {
    InputGraph graph({ input }, false, parameters, Alphabet(), mapping);
    if(prefix.empty()) { GCSA index(graph, parameters); }
    else
    {
      GCSA::buildAndStore(graph, parameters, prefix + GCSA::EXTENSION);
      LCPArray::buildAndStore(graph, parameters, prefix + LCPArray::EXTENSION);
    }
  }
  catch(const ConstructionStopped& event)
  {
    stopped = (event.completed_phase == "final-events");
  }
  catch(...)
  {
    std::cerr.rdbuf(old_stderr); Verbosity::set(Verbosity::SILENT); throw;
  }
  std::cerr.rdbuf(old_stderr); Verbosity::set(Verbosity::SILENT);
  require(stopped == prefix.empty());
  setSwitch("GCSA_SERIAL_FINAL_SCAN", false);
  setSwitch("GCSA_FINAL_SCAN_PREVIOUS_ON_DISK", false);
  return log.str();
}

// The number before `suffix` on the pipeline log line.
size_type
pipelineCount(const std::string& log, const std::string& suffix)
{
  const std::string prefix = "final scan pipeline: ";
  size_type line = log.find(prefix); require(line != std::string::npos);
  const std::string text = log.substr(line, log.find('\n', line) - line);
  size_type at = text.find(suffix); require(at != std::string::npos);
  size_type start = text.rfind(' ', at - 1) + 1;
  return std::stoull(text.substr(start, at - start));
}

} // namespace

int
main()
{
  Verbosity::set(Verbosity::SILENT);

  const std::vector<ScanCase> cases = {
    // seed, chars, de Bruijn, linear, length, shared %, N copies, I/O
    // buffer, threads, sample period, mapping, framed, previous on disk
    { 1, "ACGT", true, false, 0, 30, 2051, 16 * KILOBYTE, 4, 64, true, true, false },
    { 2, "ACGT", false, false, 6000, 20, 2051, 16 * KILOBYTE, 4, 64, false, true, false },
    { 3, "ACG", false, false, 5000, 10, 0, 20 * KILOBYTE, 2, 1, true, false, false },
    { 4, "AT", false, false, 5000, 30, 400, 32 * KILOBYTE, 8, 7, false, true, true },
    { 5, "ACGT", false, false, 6000, 5, 1000, 17 * KILOBYTE, 3, 1024, true, true, false },
    { 6, "CGT", false, false, 5000, 50, 300, 40 * KILOBYTE, 4, 2, true, false, true },
    { 7, "ACGT", false, false, 8000, 0, 2051, 24 * KILOBYTE, 5, 1000000, false, false, false },
    { 8, "ACGT", false, true, 6000, 10, 1000, 16 * KILOBYTE, 4, 64, true, true, false },
    { 9, "AT", false, true, 5000, 0, 0, 21 * KILOBYTE, 6, 1000000, false, false, true },
  };

  for(size_type number = 0; number < cases.size(); number++)
  {
    const ScanCase& scan = cases[number];
    const std::string root = makeTempRoot("gcsa-final-scan");
    const std::string input = root + "/input.gcsa2";
    writeInput(input, scan);
    std::string mapping;
    if(scan.mapping) { mapping = root + "/mapping"; writeMapping(mapping); }
    const std::string serial_root = root + "/serial", pipelined_root = root + "/pipelined";
    std::filesystem::create_directory(serial_root);
    std::filesystem::create_directory(pipelined_root);

    const std::string serial_log = build(scan, input, mapping, serial_root, false);
    const std::string pipelined_log = build(scan, input, mapping, pipelined_root, true);
    require(serial_log.find("final scan pipeline:") == std::string::npos);
    require(pipelined_log.find("final scan pipeline:") != std::string::npos);
    const size_type batches = pipelineCount(pipelined_log, " batch(es)");
    const size_type serial_paths = pipelineCount(pipelined_log, " serial path(s)");
    require(batches > 1);
    // At 1,000 start nodes or more the N path outgrows every case's arena;
    // a random path with a large set may become serial as well.
    require(scan.heavy_copies < 1000 || serial_paths >= 1);
    require((pipelined_log.find("previous occurrences on disk") != std::string::npos) ==
      scan.previous_on_disk);

    const std::map<std::string, std::vector<char>> serial_events =
      workspaceArtifacts(serial_root, "final--events--");
    const std::map<std::string, std::vector<char>> pipelined_events =
      workspaceArtifacts(pipelined_root, "final--events--");
    require(serial_events.size() >= 7);
    require(serial_events == pipelined_events);
    std::cerr << "test_final_scan: case " << (number + 1) << ": " << batches
              << " batches, " << serial_paths << " serial path(s), "
              << serial_events.size() << " identical event artifacts" << std::endl;

    // Complete the de Bruijn case from scratch on both routes: the public
    // files must match, and the pipelined index must answer every query.
    if(scan.debruijn)
    {
      const std::string serial_full = root + "/serial-full", pipelined_full = root + "/pipelined-full";
      std::filesystem::create_directory(serial_full);
      std::filesystem::create_directory(pipelined_full);
      build(scan, input, mapping, serial_full, false, serial_full + "/index");
      build(scan, input, mapping, pipelined_full, true, pipelined_full + "/index");
      require(readFile(serial_full + "/index" + GCSA::EXTENSION) ==
        readFile(pipelined_full + "/index" + GCSA::EXTENSION));
      require(readFile(serial_full + "/index" + LCPArray::EXTENSION) ==
        readFile(pipelined_full + "/index" + LCPArray::EXTENSION));
      ConstructionParameters parameters;
      InputGraph graph({ input }, false, parameters, Alphabet(), mapping);
      GCSA index; LCPArray lcp;
      require(sdsl::load_from_file(index, pipelined_full + "/index" + GCSA::EXTENSION));
      require(sdsl::load_from_file(lcp, pipelined_full + "/index" + LCPArray::EXTENSION));
      require(verifyIndex(index, &lcp, graph));
    }
    std::filesystem::remove_all(root);
  }
  return 0;
}
