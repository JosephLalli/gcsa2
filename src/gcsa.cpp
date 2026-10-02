#include <gcsa/algorithms.h>
#include <gcsa/path_graph_external.h>
#include <gcsa/checkpoint.h>
#include <gcsa/disk_array.h>
#include <gcsa/external_preprocessing.h>
#include <gcsa/final_events.h>
#include <gcsa/internal.h>
#include <gcsa/path_graph.h>
#include <gcsa/workspace.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <filesystem>
#include <fcntl.h>
#include <iomanip>
#include <limits>
#include <memory>
#include <random>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>
#include <unordered_set>

namespace gcsa
{

//------------------------------------------------------------------------------

// Other class variables

const std::string GCSA::EXTENSION = ".gcsa";

namespace
{

std::string
doublingTask(size_type step)
{
  std::ostringstream result;
  result << "step-" << std::setw(2) << std::setfill('0') << step;
  return result.str();
}

// Catch up journaled cleanup when it is enabled only on a resumed run. A
// restored prune still needs this step's completed join ranges, while a
// restored extend supersedes them. Every deletion is anchored directly to the
// validated current frontier rather than depending on an intermediate
// predecessor that may itself already have been retired.
void
retireBeforeFrontier(BuildWorkspace& workspace,
  const std::string& successor_task, const std::string& successor_phase,
  size_type completed_step, bool extend_complete)
{
  if(workspace.task_completed("initial", "paths"))
  {
    workspace.retire_obsolete("initial", "paths",
      successor_task, successor_phase);
  }
  for(size_type step = 1; step <= completed_step; step++)
  {
    const std::string task = doublingTask(step);
    if(!(task == successor_task && successor_phase == "prune") &&
       workspace.task_completed(task, "prune"))
    {
      workspace.retire_obsolete(task, "prune",
        successor_task, successor_phase);
    }
    if(!(task == successor_task && successor_phase == "extend") &&
       workspace.task_completed(task, "extend"))
    {
      workspace.retire_obsolete(task, "extend",
        successor_task, successor_phase);
    }
    if(step < completed_step || extend_complete)
    {
      workspace.retire_obsolete_family(task + "-join-", "join-partition",
        successor_task, successor_phase);
      workspace.retire_obsolete_family(task + "-msd-plan-", "join-plan",
        successor_task, successor_phase);
      workspace.retire_obsolete_family(task + "-range-plan-", "range-plan",
        successor_task, successor_phase);
    }
  }
}

std::uint64_t
constructionFileChecksum(const std::string& filename, size_type buffer_bytes)
{
  std::ifstream input;
  input.rdbuf()->pubsetbuf(nullptr, 0);
  input.open(filename.c_str(), std::ios_base::binary);
  if(!input) { throw std::runtime_error("GCSA::GCSA(): cannot checksum " + filename); }
  std::vector<std::uint8_t> buffer(std::max(static_cast<size_type>(1), buffer_bytes));
  std::uint64_t checksum = 1469598103934665603ULL;
  while(input)
  {
    input.read(reinterpret_cast<char*>(buffer.data()), buffer.size());
    std::streamsize bytes = input.gcount();
    if(bytes > 0) { checksum = BuildWorkspace::checksum(buffer.data(), bytes, checksum); }
  }
  if(!input.eof()) { throw std::runtime_error("GCSA::GCSA(): failed to checksum " + filename); }
  return checksum;
}

BuildWorkspace::Settings
constructionSemanticSettings(const InputGraph& graph,
  const ConstructionParameters& parameters, size_type checksum_buffer)
{
  BuildWorkspace::Settings settings;
  settings["gcsa_format"] = std::to_string(Version::GCSA_VERSION);
  settings["lcp_format"] = std::to_string(Version::LCP_VERSION);
  settings["binary_input"] = (graph.binary ? "true" : "false");
  // A workspace generated before disk-first preprocessing cannot safely mix
  // old initial/path checkpoints with the v1 reduced key/start artifacts.
  settings["external_preprocessing"] = "v1";
  settings["external_final_events"] = "v1";
  // Keep the v2 compatibility tag: v3 metadata and the framed reader both
  // deliberately accept the older raw payloads, so changing the semantic
  // fingerprint would reject workspaces this build can safely resume.
  settings["path_graph_checkpoint"] = "v2-adopted-raw";
  // Both formats affect deterministic transient ordering and therefore the
  // semantic range names used by resumable join-partition tasks.
  settings["external_join_planner"] = "v2-msd-range-pack";
  settings["path_sort_run_codec"] = "v2-left-context-groups";
  settings["kmer_length"] = std::to_string(graph.k());
  settings["doubling_steps"] = std::to_string(parameters.getSteps());
  settings["sample_period"] = std::to_string(parameters.getSamplePeriod());
  settings["lcp_branching"] = std::to_string(parameters.getLCPBranching());
  settings["input_count"] = std::to_string(graph.files());
  for(size_type file = 0; file < graph.files(); file++)
  {
    struct stat info;
    if(::stat(graph.filenames[file].c_str(), &info) != 0)
    {
      throw std::runtime_error("GCSA::GCSA(): cannot stat " + graph.filenames[file]);
    }
    std::string prefix = "input_" + std::to_string(file) + "_";
    settings[prefix + "path"] = graph.filenames[file];
    settings[prefix + "bytes"] = std::to_string(static_cast<std::uint64_t>(info.st_size));
    settings[prefix + "checksum"] = std::to_string(
      constructionFileChecksum(graph.filenames[file], checksum_buffer));
    settings[prefix + "logical_id"] = std::to_string(file);
  }
  settings["mapping_path"] = graph.mapping_name;
  if(!graph.mapping_name.empty())
  {
    struct stat info;
    if(::stat(graph.mapping_name.c_str(), &info) != 0)
    {
      throw std::runtime_error("GCSA::GCSA(): cannot stat " + graph.mapping_name);
    }
    settings["mapping_bytes"] = std::to_string(static_cast<std::uint64_t>(info.st_size));
    settings["mapping_checksum"] = std::to_string(
      constructionFileChecksum(graph.mapping_name, checksum_buffer));
  }
  return settings;
}

BuildWorkspace::Settings
constructionOperationalSettings(const ConstructionParameters& parameters)
{
  BuildWorkspace::Settings settings;
  settings["memory_limit"] = std::to_string(parameters.getMemoryLimitBytes());
  settings["disk_limit"] = std::to_string(parameters.getLimitBytes());
  settings["io_buffer_size"] = std::to_string(parameters.getIOBufferSize());
  settings["sort_run_size"] = std::to_string(parameters.getSortRunSize());
  settings["join_partition_size"] = std::to_string(parameters.getJoinPartitionSize());
  settings["sort_run_size_mode"] = (parameters.sortRunSizeIsAutomatic() ? "auto" : "explicit");
  settings["join_partition_size_mode"] = (parameters.joinPartitionSizeIsAutomatic() ? "auto" : "explicit");
  settings["merge_fan_in"] = std::to_string(parameters.getMergeFanIn());
  settings["max_open_files"] = std::to_string(parameters.getMaxOpenFiles());
  settings["process_workers"] = std::to_string(parameters.getProcessWorkers());
  settings["temp_compression"] = tempCompressionName(parameters.getTempCompression());
  settings["compression_block_size"] = std::to_string(parameters.getCompressionBlockSize());
  settings["compression_workers"] = std::to_string(parameters.getCompressionWorkers());
  settings["compression_level"] = std::to_string(parameters.getCompressionLevel());
  settings["clean_obsolete"] = (parameters.getCleanObsolete() ? "true" : "false");
  settings["threads"] = std::to_string(omp_get_max_threads());
  return settings;
}

// A sub-phase timer and I/O probe.
//
// The construction path has four timers, and each spans several distinct
// operations: "Merging" covers the mapper build and the merged-graph merge,
// "Construction" covers the final-event scan and the component build. A phase total therefore cannot say which operation costs
// what, which is exactly the question an optimization pass has to answer.
// These probes print under EXTENDED verbosity only and serialize nothing.
struct SubPhaseProbe
{
  std::string name;
  double start;
  size_type read_start, write_start;

  explicit SubPhaseProbe(const std::string& phase_name) :
    name(phase_name), start(readTimer()),
    read_start(DiskIO::read_volume), write_start(DiskIO::write_volume) { }

  void report() const
  {
    if(Verbosity::level < Verbosity::EXTENDED) { return; }
    std::cerr << "GCSA::GCSA(): subphase " << this->name << ": "
              << (readTimer() - this->start) << " seconds, "
              << inGigabytes(DiskIO::read_volume - this->read_start) << " GB read, "
              << inGigabytes(DiskIO::write_volume - this->write_start) << " GB written"
              << std::endl;
  }
};

void
stopAfterCommittedPhase(const ConstructionParameters& parameters,
  const std::string& completed_phase)
{
  if(parameters.getStopAfter() == completed_phase)
  {
    throw ConstructionStopped(completed_phase);
  }
}

// The direct final-component path normally publishes through
// storeFinalComponents(). An empty InputGraph has no event frontier, but it
// must still replace the destination with the ordinary serialized empty GCSA.
// Keep the same write/sync/rename/sync protocol used for nonempty indexes.
void
storeEmptyIndexAtomically(const GCSA& index, const std::string& filename)
{
  const std::string partial = filename + "." +
    std::to_string(static_cast<unsigned long long>(::getpid())) + ".partial";
  bool published = false;
  try
  {
    std::ofstream output(partial.c_str(), std::ios_base::binary | std::ios_base::trunc);
    if(!output) { throw std::runtime_error("cannot open partial empty index " + partial); }
    index.serialize(output);
    output.flush(); output.close();
    if(!output) { throw std::runtime_error("cannot write partial empty index " + partial); }

    int descriptor = ::open(partial.c_str(), O_RDONLY);
    if(descriptor < 0 || ::fdatasync(descriptor) != 0)
    {
      if(descriptor >= 0) { ::close(descriptor); }
      throw std::runtime_error("cannot sync partial empty index " + partial);
    }
    if(::close(descriptor) != 0 || ::rename(partial.c_str(), filename.c_str()) != 0)
    {
      throw std::runtime_error("cannot publish empty index " + filename);
    }
    published = true;

    std::filesystem::path parent = std::filesystem::path(filename).parent_path();
    if(parent.empty()) { parent = "."; }
    descriptor = ::open(parent.c_str(), O_RDONLY | O_DIRECTORY);
    if(descriptor < 0 || ::fsync(descriptor) != 0)
    {
      if(descriptor >= 0) { ::close(descriptor); }
      throw std::runtime_error("cannot sync empty index directory " + filename);
    }
    if(::close(descriptor) != 0)
    {
      throw std::runtime_error("cannot close empty index directory " + filename);
    }
  }
  catch(...)
  {
    if(!published) { ::unlink(partial.c_str()); }
    throw;
  }
}

} // namespace

//------------------------------------------------------------------------------

GCSA::GCSA() :
  header(), alpha(),
  fast_bwt(this->alpha.sigma), fast_rank(this->alpha.sigma),
  sparse_bwt(this->alpha.sigma), sparse_rank(this->alpha.sigma),
  edges(), edge_rank(),
  sampled_paths(), sampled_path_rank(),
  stored_samples(), samples(), sample_select(),
  extra_pointers(), redundant_pointers()
{
}

GCSA::GCSA(const GCSA& source)
{
  this->copy(source);
}

GCSA::GCSA(GCSA&& source) noexcept
{
  *this = std::move(source);
}

GCSA::~GCSA()
{
}

void
GCSA::swap(GCSA& another) noexcept
{
  if(this != &another)
  {
    this->header.swap(another.header);
    this->alpha.swap(another.alpha);

    this->fast_bwt.swap(another.fast_bwt);
    this->fast_rank.swap(another.fast_rank);

    this->sparse_bwt.swap(another.sparse_bwt);
    this->sparse_rank.swap(another.sparse_rank);

    this->edges.swap(another.edges);
    sdsl::util::swap_support(this->edge_rank, another.edge_rank, &(this->edges), &(another.edges));

    this->sampled_paths.swap(another.sampled_paths);
    sdsl::util::swap_support(this->sampled_path_rank, another.sampled_path_rank, &(this->sampled_paths), &(another.sampled_paths));

    this->stored_samples.swap(another.stored_samples);
    this->samples.swap(another.samples);
    sdsl::util::swap_support(this->sample_select, another.sample_select, &(this->samples), &(another.samples));

    this->extra_pointers.swap(another.extra_pointers);
    this->redundant_pointers.swap(another.redundant_pointers);

    this->setVectors();
  }
}

GCSA&
GCSA::operator=(const GCSA& source)
{
  if(this != &source) { this->copy(source); }
  return *this;
}

GCSA&
GCSA::operator=(GCSA&& source) noexcept
{
  if(this != &source)
  {
    this->header = std::move(source.header);
    this->alpha = std::move(source.alpha);

    this->fast_bwt = std::move(source.fast_bwt);
    this->fast_rank = std::move(source.fast_rank);

    this->sparse_bwt = std::move(source.sparse_bwt);
    this->sparse_rank = std::move(source.sparse_rank);

    this->edges = std::move(source.edges);
    this->edge_rank = std::move(source.edge_rank);

    this->sampled_paths = std::move(source.sampled_paths);
    this->sampled_path_rank = std::move(source.sampled_path_rank);

    this->stored_samples = std::move(source.stored_samples);
    this->samples = std::move(source.samples);
    this->sample_select = std::move(source.sample_select);

    this->extra_pointers = std::move(source.extra_pointers);
    this->redundant_pointers = std::move(source.redundant_pointers);

    this->setVectors();
  }
  return *this;
}

GCSA::size_type
GCSA::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += this->header.serialize(out, child, "header");
  written_bytes += this->alpha.serialize(out, child, "alpha");

  for(size_type comp = 0; comp < this->alpha.sigma; comp++)
  {
    written_bytes += this->fast_bwt[comp].serialize(out, child, "fast_bwt");
  }
  for(size_type comp = 0; comp < this->alpha.sigma; comp++)
  {
    written_bytes += this->fast_rank[comp].serialize(out, child, "fast_rank");
  }

  for(size_type comp = 0; comp < this->alpha.sigma; comp++)
  {
    written_bytes += this->sparse_bwt[comp].serialize(out, child, "sparse_bwt");
  }
  for(size_type comp = 0; comp < this->alpha.sigma; comp++)
  {
    written_bytes += this->sparse_rank[comp].serialize(out, child, "sparse_rank");
  }

  written_bytes += this->edges.serialize(out, child, "edges");
  written_bytes += this->edge_rank.serialize(out, child, "edge_rank");

  written_bytes += this->sampled_paths.serialize(out, child, "sampled_paths");
  written_bytes += this->sampled_path_rank.serialize(out, child, "sampled_path_rank");

  written_bytes += this->stored_samples.serialize(out, child, "stored_samples");
  written_bytes += this->samples.serialize(out, child, "samples");
  written_bytes += this->sample_select.serialize(out, child, "sample_select");

  written_bytes += this->extra_pointers.serialize(out, child, "extra_pointers");
  written_bytes += this->redundant_pointers.serialize(out, child, "redundant_pointers");

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
GCSA::load(std::istream& in)
{
  this->header.load(in);
  if(!(this->header.check()))
  {
    std::stringstream ss;
    ss << "GCSA::load(): Invalid header: " << this->header;
    throw std::runtime_error(ss.str());
  }
  this->alpha.load(in);

  this->fast_bwt.resize(this->alpha.sigma); this->fast_rank.resize(this->alpha.sigma);
  for(size_type comp = 0; comp < this->alpha.sigma; comp++) { this->fast_bwt[comp].load(in); }
  for(size_type comp = 0; comp < this->alpha.sigma; comp++) { this->fast_rank[comp].load(in, &(this->fast_bwt[comp])); }

  this->sparse_bwt.resize(this->alpha.sigma); this->sparse_rank.resize(this->alpha.sigma);
  for(size_type comp = 0; comp < this->alpha.sigma; comp++) { this->sparse_bwt[comp].load(in); }
  for(size_type comp = 0; comp < this->alpha.sigma; comp++) { this->sparse_rank[comp].load(in, &(this->sparse_bwt[comp])); }

  this->edges.load(in);
  this->edge_rank.load(in, &(this->edges));

  this->sampled_paths.load(in);
  this->sampled_path_rank.load(in, &(this->sampled_paths));

  this->stored_samples.load(in);
  this->samples.load(in);
  this->sample_select.load(in, &(this->samples));

  this->extra_pointers.load(in);
  this->redundant_pointers.load(in);
}

void
GCSA::copy(const GCSA& source)
{
  this->header = source.header;
  this->alpha = source.alpha;

  this->fast_bwt = source.fast_bwt;
  this->fast_rank = source.fast_rank;

  this->sparse_bwt = source.sparse_bwt;
  this->sparse_rank = source.sparse_rank;

  this->edges = source.edges;
  this->edge_rank = source.edge_rank;

  this->sampled_paths = source.sampled_paths;
  this->sampled_path_rank = source.sampled_path_rank;

  this->stored_samples = source.stored_samples;
  this->samples = source.samples;
  this->sample_select = source.sample_select;

  this->extra_pointers = source.extra_pointers;
  this->redundant_pointers = source.redundant_pointers;

  this->setVectors();
}

void
GCSA::setVectors()
{
  for(size_type comp = 0; comp < this->alpha.sigma; comp++)
  {
    this->fast_rank[comp].set_vector(&(this->fast_bwt[comp]));
    this->sparse_rank[comp].set_vector(&(this->sparse_bwt[comp]));
  }

  this->edge_rank.set_vector(&(this->edges));

  this->sampled_path_rank.set_vector(&(this->sampled_paths));

  this->sample_select.set_vector(&(this->samples));
}

//------------------------------------------------------------------------------

struct MergedGraphReader
{
  ReadBuffer<PathNode>            paths;
  ReadBuffer<PathNode::rank_type> labels;
  ReadBuffer<range_type>          from_nodes;

  size_type path, rank, from, path_count;

  const DeBruijnGraph*            mapper;
  const sdsl::int_vector<0>*      last_char;

  void init(const MergedGraph& graph, const DeBruijnGraph* _mapper,
    const sdsl::int_vector<0>* _last_char,
    size_type buffer_bytes = ReadBuffer<PathNode>::DEFAULT_BUFFER_BYTES,
    bool release_cache = false);
  void init(const MergedGraph& graph, size_type comp,
    size_type buffer_bytes = ReadBuffer<PathNode>::DEFAULT_BUFFER_BYTES,
    bool release_cache = false, bool read_from_nodes = true);
  void initFromOnly(const MergedGraph& graph,
    size_type buffer_bytes = ReadBuffer<PathNode>::DEFAULT_BUFFER_BYTES,
    bool release_cache = false, size_type comp = MergedGraph::UNKNOWN);
  void close();

  void seek(bool seek_labels = true);

  inline void advance()
  {
    if(this->path >= this->path_count || this->path + 1 >= this->path_count)
    {
      return;
    }
    this->path++;
    this->seek();
  }

  void predecessor(comp_type comp, PathLabel& first, PathLabel& last);

  /*
    Does paths[path + offset] intersect with the given range of labels?
  */
  bool intersect(const PathLabel& first, const PathLabel& last, size_type offset);

  void fromNodes(std::vector<node_type>& results, const NodeMapping& mapping);
  void fromNodes(SpillableNodeSet& results, const NodeMapping& mapping);
  void fromNodes(SpillableNodeSet& results, const NodeMapping& mapping,
    node_type primary);
};

/*
  Source-label reader for bounded predecessor batches. It deliberately omits
  the from-node stream: the ordered reader below remains responsible for that
  stream, so preparing a batch neither advances nor rereads global scan state.
*/
struct MergedGraphLabelReader
{
  ReadBuffer<PathNode> paths;
  ReadBuffer<PathNode::rank_type> labels;
  size_type path;

  MergedGraphLabelReader() : paths(), labels(), path(0) { }

  void init(const MergedGraph& graph, size_type buffer_bytes)
  {
    this->paths.open(graph.path_name, buffer_bytes, true);
    this->labels.open(graph.rank_name, buffer_bytes, true);
    this->path = 0;
    if(this->paths.size() > 0)
    {
      this->paths.seek(0);
      this->labels.seek(this->paths[0].pointer());
    }
  }

  void advance()
  {
    if(this->path + 1 >= this->paths.size()) { return; }
    this->path++;
    this->paths.seek(this->path);
    this->labels.seek(this->paths[this->path].pointer());
  }

  void close()
  {
    this->paths.close(); this->labels.close(); this->path = 0;
  }

  MergedGraphLabelReader(const MergedGraphLabelReader&) = delete;
  MergedGraphLabelReader& operator=(const MergedGraphLabelReader&) = delete;
};

void
MergedGraphReader::init(const MergedGraph& graph,
  const DeBruijnGraph* _mapper, const sdsl::int_vector<0>* _last_char,
  size_type buffer_bytes, bool release_cache)
{
  this->paths.open(graph.path_name, buffer_bytes, release_cache);
  this->labels.open(graph.rank_name, buffer_bytes, release_cache);
  this->from_nodes.open(graph.from_name, buffer_bytes, release_cache);

  this->path = this->rank = this->from = 0; this->path_count = graph.size();
  this->seek();

  this->mapper = _mapper;
  this->last_char = _last_char;
}

void
MergedGraphReader::init(const MergedGraph& graph, size_type comp,
  size_type buffer_bytes, bool release_cache, bool read_from_nodes)
{
  this->paths.open(graph.path_name, buffer_bytes, release_cache);
  this->labels.open(graph.rank_name, buffer_bytes, release_cache);
  if(read_from_nodes)
  {
    this->from_nodes.open(graph.from_name, buffer_bytes, release_cache);
  }

  this->path = graph.next[comp]; this->path_count = graph.size();
  this->from = graph.next_from[comp];
  this->seek();

  this->mapper = nullptr;
  this->last_char = nullptr;
}

void
MergedGraphReader::initFromOnly(const MergedGraph& graph,
  size_type buffer_bytes, bool release_cache, size_type comp)
{
  this->from_nodes.open(graph.from_name, buffer_bytes, release_cache);
  this->path = (comp == MergedGraph::UNKNOWN ? 0 : graph.next[comp]);
  this->from = (comp == MergedGraph::UNKNOWN ? 0 : graph.next_from[comp]);
  this->rank = 0; this->path_count = graph.size();
  this->seek(false);
  this->mapper = nullptr; this->last_char = nullptr;
}

void
MergedGraphReader::close()
{
  this->paths.close(),
  this->labels.close();
  this->from_nodes.close();

  this->path = this->rank = this->from = this->path_count = 0;
}

void
MergedGraphReader::seek(bool seek_labels)
{
  // A sparse alphabet component may have no destination paths. Keep its
  // reader at a valid end frontier; callers that never observe the component
  // need no dummy path, while an inconsistent edge still fails the explicit
  // destination bounds check in the external scan.
  if(this->path >= this->path_count)
  {
    this->rank = this->labels.size();
    this->from = this->from_nodes.size();
    return;
  }
  if(this->paths.isOpen())
  {
    this->paths.seek(this->path);
    this->rank = this->paths[this->path].pointer();
    if(seek_labels) { this->labels.seek(this->rank); }
  }
  // Move the window before walking. init() starts `from` at the component's
  // first start node, deep in the stream for every component but the first,
  // and indexing a ReadBuffer beyond its window first reads every element in
  // between into it. Without this seek the final scan's per-component readers
  // read each stream prefix from offset 0 into memory during setup; on the
  // joint chr2+chr18 graph that setup took about 88 s while RSS rose from
  // 15.5 to 32.1 GiB and stayed there until the scan ended, consistent with
  // an estimated 47 GiB of prefix reads. The walk itself is unchanged.
  this->from_nodes.seek(this->from);
  while(this->from < this->from_nodes.size() && this->from_nodes[this->from].first < this->path)
  {
    this->from++;
  }
  this->from_nodes.seek(this->from);
}

template<class LabelReader>
void
predecessorRange(const PathNode& curr, comp_type comp,
  const DeBruijnGraph& mapper, const sdsl::int_vector<0>& last_char,
  const LabelReader& label,
  PathLabel& first, PathLabel& last)
{
  size_type i = 0;
  first.first = true; last.first = false;

  // Handle the common prefix of the labels.
  while(i < curr.lcp())
  {
    PathNode::rank_type rank = label(i);
    first.label[i] = last.label[i] = mapper.LF(rank, comp);
    comp = last_char[rank];
    i++;
  }

  // Handle the diverging suffixes of the labels.
  comp_type first_comp = comp, last_comp = comp;
  if(i < curr.order())
  {
    PathNode::rank_type first_rank = label(i);
    PathNode::rank_type last_rank = label(i + 1);
    first.label[i] = mapper.LF(first_rank, first_comp);
    first_comp = last_char[first_rank];
    last.label[i] = mapper.LF(last_rank, last_comp);
    last_comp = last_char[last_rank];
    i++;
  }
  if(i < PathLabel::LABEL_LENGTH)
  {
    first.label[i] = mapper.node_rank(mapper.alpha.C[first_comp]);
    last.label[i] = mapper.node_rank(mapper.alpha.C[last_comp + 1]) - 1;
    i++;
  }

  // Pad the labels.
  while(i < PathLabel::LABEL_LENGTH)
  {
    first.label[i] = 0; last.label[i] = PathLabel::NO_RANK; i++;
  }
}

void
MergedGraphReader::predecessor(comp_type comp, PathLabel& first, PathLabel& last)
{
  const PathNode& curr = this->paths[this->path];
  predecessorRange(curr, comp, *(this->mapper), *(this->last_char),
    [this, &curr](size_type i) { return this->labels[curr.pointer() + i]; },
    first, last);
}

inline PathLabel
firstLabel(const PathNode& path, ReadBuffer<PathNode::rank_type>& labels)
{
  PathLabel res; res.first = true;
  size_type limit = std::min(path.order(), PathLabel::LABEL_LENGTH);
  for(size_type i = 0; i < limit; i++) { res.label[i] = path.firstLabel(i, labels); }
  for(size_type i = limit; i < PathLabel::LABEL_LENGTH; i++) { res.label[i] = 0; }
  return res;
}

inline PathLabel
lastLabel(const PathNode& path, ReadBuffer<PathNode::rank_type>& labels)
{
  PathLabel res; res.first = false;
  size_type limit = std::min(path.order(), PathLabel::LABEL_LENGTH);
  for(size_type i = 0; i < limit; i++) { res.label[i] = path.lastLabel(i, labels); }
  for(size_type i = limit; i < PathLabel::LABEL_LENGTH; i++) { res.label[i] = PathLabel::NO_RANK; }
  return res;
}

/*
  Does the path node intersect with the given range of labels?
*/
bool
MergedGraphReader::intersect(const PathLabel& first, const PathLabel& last, size_type offset)
{
  PathLabel my_first = firstLabel(this->paths[this->path + offset], this->labels);
  if(my_first <= first)
  {
    PathLabel my_last = lastLabel(this->paths[this->path + offset], this->labels);
    return (first <= my_last);
  }
  else
  {
    return my_first <= last;
  }
}

void
MergedGraphReader::fromNodes(std::vector<node_type>& results, const NodeMapping& mapping)
{
  results.clear();
  results.push_back(this->paths[this->path].from);

  size_type old_pointer = this->from;
  while(this->from < this->from_nodes.size() && this->from_nodes[this->from].first == this->path)
  {
    results.push_back(this->from_nodes[this->from].second); this->from++;
  }
  this->from = old_pointer;

  Node::map(results, mapping);
  removeDuplicates(results, false);
}

void
MergedGraphReader::fromNodes(SpillableNodeSet& results, const NodeMapping& mapping)
{
  this->fromNodes(results, mapping, this->paths[this->path].from);
}

void
MergedGraphReader::fromNodes(SpillableNodeSet& results,
  const NodeMapping& mapping, node_type primary)
{
  results.clear();
  const auto mapped = [&mapping](node_type node) {
    return (mapping.empty() ? node :
      Node::encode(mapping(Node::id(node)), Node::offset(node), Node::rc(node)));
  };
  results.push_back(mapped(primary));

  size_type old_pointer = this->from;
  while(this->from < this->from_nodes.size() && this->from_nodes[this->from].first == this->path)
  {
    results.push_back(mapped(this->from_nodes[this->from].second)); this->from++;
  }
  this->from = old_pointer;
  results.finish();
}

//------------------------------------------------------------------------------

namespace
{

struct ExternalFinalScanStats
{
  DiskBackedArray64::Stats previous_occurrences, suffix_tree_stack;
  size_type from_node_spills, maximum_from_nodes;
  size_type predecessor_workers, predecessor_batches;
  size_type maximum_predecessor_batch, predecessor_buffer_bytes;
  size_type prepared_from_paths, prepared_from_ranks;
  size_type prepared_from_spill_fallbacks, from_preparation_buffer_bytes;
  bool predecessor_parallel_fallback;
  bool restored;

  ExternalFinalScanStats() : previous_occurrences(), suffix_tree_stack(),
    from_node_spills(0), maximum_from_nodes(0), predecessor_workers(1),
    predecessor_batches(0), maximum_predecessor_batch(0),
    predecessor_buffer_bytes(0), prepared_from_paths(0),
    prepared_from_ranks(0), prepared_from_spill_fallbacks(0),
    from_preparation_buffer_bytes(0), predecessor_parallel_fallback(false),
    restored(false) { }
};

/*
  Immutable work for the unordered half of one final-scan batch.

  Each component worker follows exactly the same monotone destination-reader
  trajectory as the serial scan, but it owns a distinct reader and destination
  slot. Path workers independently sort/deduplicate the mapped current-node
  slice and resolve each node's immutable from-rank. The globally ordered
  consumer later emits those destinations and alone mutates prev_occ, the
  suffix-tree stack, and the event streams.

  Current-node payload is held in separately admitted arenas. A path that does
  not fit is the last path in its batch and uses the existing spillable set, so
  no variable-size allocation can inflate this fixed record.
*/
struct FinalPredecessorWork
{
  PathNode path;
  std::array<PathNode::rank_type, PathLabel::LABEL_LENGTH + 1> labels;
  std::array<size_type, FinalEventMetadata::MAX_SIGMA> destinations;
  std::array<node_type, FinalEventMetadata::MAX_SIGMA> destination_primary;
  size_type from_offset, from_count;
  bool from_spilled, sample_by_period;
};

constexpr size_type MIN_FINAL_PREDECESSOR_BATCH = 64;

size_type
checkedProduct(size_type first, size_type second, const std::string& label)
{
  if(second != 0 && first > std::numeric_limits<size_type>::max() / second)
  {
    throw std::runtime_error("GCSA::GCSA(): " + label + " size overflows");
  }
  return first * second;
}

node_type
mappedFinalNode(node_type node, const NodeMapping& mapping)
{
  return (mapping.empty() ? node :
    Node::encode(mapping(Node::id(node)), Node::offset(node), Node::rc(node)));
}

size_type
finalFromRank(node_type node, const sdsl::sd_vector<>& from_nodes)
{
  // successor() performs the membership test and rank lookup in one succinct
  // descent. This helper is pure over the immutable sd_vector and is therefore
  // safe to run in the unordered half of a batch.
  if(node >= from_nodes.size())
  {
    throw std::runtime_error(
      "GCSA::GCSA(): start node " + Node::decode(node) +
      " (encoded " + std::to_string(node) +
      ") is absent from the external start-node index of universe " +
      std::to_string(from_nodes.size()));
  }
  auto occurrence = from_nodes.successor(node);
  if(occurrence->second != node)
  {
    throw std::runtime_error(
      "GCSA::GCSA(): start node " + Node::decode(node) +
      " (encoded " + std::to_string(node) +
      ") is absent from the external start-node index of universe " +
      std::to_string(from_nodes.size()));
  }
  return occurrence->first;
}

template<class Callback>
void
forEachFinalFromNode(const FinalPredecessorWork& work,
  const std::vector<node_type>& prepared_nodes, SpillableNodeSet& spilled_nodes,
  Callback callback)
{
  if(work.from_spilled)
  {
    spilled_nodes.rewind();
    node_type node;
    while(spilled_nodes.next(node)) { callback(node); }
    return;
  }
  for(size_type j = 0; j < work.from_count; j++)
  {
    callback(prepared_nodes[work.from_offset + j]);
  }
}

/*
  FORK: sizing a final-scan stream buffer.

  All three of the final scan's buffer sizes -- readers, writers and the
  start-node sets -- are getIOBufferSize() unless the budget is too small to
  afford it. That was sixty lines of nested max/min a reader had to evaluate by
  hand to discover it is a constant. The share test states it instead: above a
  comfortable budget it passes and the answer is the target; only a tight budget
  falls through to the arithmetic.

  Substituting the chr21 configuration, the fall-through binds below roughly
  6.4 GiB available for readers, 1.63 GiB for writers and 1 GiB for the node
  sets -- so on every run this subsystem has actually been used for, the answer
  is the constant. The arithmetic is not dead, though: it is the only thing that
  lets a deliberately tight --gcsa-memory-limit run proceed instead of failing
  its reservation, which is the case the whole subsystem exists to serve. Hence
  a branch rather than a deletion.

  The tight branch is the original expression verbatim, and the ample branch is
  what that expression returns whenever `share >= target` and the target clears
  the floor, so this cannot change a size on any budget.
*/
size_type
finalScanBuffer(size_type target, size_type minimum, size_type available,
  size_type streams_sharing)
{
  size_type share = available /
    std::max(static_cast<size_type>(1), streams_sharing);
  if(share >= target && target >= minimum) { return target; }
  return std::max(minimum, std::min(target, share));
}

/*
  Scan the final MergedGraph once and publish compact immutable events. This is
  intentionally one atomic task: prev_occ and the suffix-tree traversal stack
  are mutable disk arrays, and a mid-scan checkpoint would require a committed,
  idempotent assignment-log protocol. A crash before the task marker therefore
  discards the event artifacts and restarts this scan, while a completed scan
  is reusable across operational memory/thread changes.
*/
FinalEventMetadata
produceExternalFinalEvents(const MergedGraph& merged_graph,
  const DeBruijnGraph& mapper, const sdsl::int_vector<0>& last_char,
  const sdsl::sd_vector<>& from_nodes,
  size_type unique_from_nodes, const InputGraph& graph,
  const ConstructionParameters& parameters, BuildWorkspace& workspace,
  FinalEventFiles& files, size_type checkpoint_buffer,
  ExternalFinalScanStats* stats)
{
  FinalEventMetadata metadata;
  FinalEventChecksums event_checksums(graph.alpha.sigma);
  if(restoreFinalEvents(workspace, files, metadata, merged_graph.size(),
    graph.alpha.sigma, checkpoint_buffer, parameters.getVerifyWorkspace()))
  {
    if(stats != nullptr) { stats->restored = true; }
    return metadata;
  }

  if(graph.alpha.sigma == 0 || graph.alpha.sigma > FinalEventMetadata::MAX_SIGMA)
  {
    throw std::runtime_error("GCSA::GCSA(): external final events require an alphabet of at most 8 components");
  }
  const size_type memory_limit = parameters.getMemoryLimitBytes();
  const size_type safety_margin = memory_limit / 8;
  MemoryBudget memory(memory_limit, safety_margin);

  const size_type requested_predecessor_workers = std::min(graph.alpha.sigma,
    static_cast<size_type>(std::max(1, omp_get_max_threads())));
  // Parallel batches move the ordered path/label streams to a source reader.
  // Each destination keeps lookup path/label streams and a separate sampling
  // from-node stream. All are monotone; the stream count matches the serial
  // route and no buffer is rewound after preparing a batch.
  const size_type reader_streams = 3 * (graph.alpha.sigma + 1) + 1;
  const size_type minimum_reader = std::max(sizeof(PathNode),
    std::max(sizeof(PathNode::rank_type), sizeof(range_type)));
  const size_type writer_streams = graph.alpha.sigma + 6;
  const size_type reader_minimum = checkedProduct(2 * reader_streams,
    minimum_reader, "minimum final reader reservation");
  const size_type writer_minimum = checkedProduct(writer_streams,
    static_cast<size_type>(16), "minimum final writer reservation");
  const size_type array_minimum =
    (unique_from_nodes > 0 ? DiskBackedArray64::minimumCacheBytes() : 0) +
    (merged_graph.size() > 0 ? DiskBackedArray64::minimumCacheBytes() : 0);

  // One merged path can represent an arbitrarily large set of original start
  // nodes. Reserve two reusable external sets before sizing the reader/writer
  // caches, because current and predecessor sets must coexist for the
  // continuation test. Their in-memory collection buffers are small for the
  // common case; an oversized set becomes a sorted, deduplicated disk stream.
  const size_type minimum_node_set = SpillableNodeSet::minimumBudget();
  const size_type non_set_minimum = reader_minimum + writer_minimum + array_minimum;
  if(memory.available() < non_set_minimum + 2 * minimum_node_set)
  {
    throw std::runtime_error("GCSA::GCSA(): memory limit cannot hold the minimum final-scan workspace");
  }
  size_type desired_node_set = std::max(minimum_node_set,
    std::min(parameters.getIOBufferSize(), parameters.getSortRunSize()));
  size_type node_set_cap = (memory.available() - non_set_minimum) / 2;
  size_type node_set_budget = std::min(node_set_cap,
    finalScanBuffer(desired_node_set, minimum_node_set, memory.available(), 16));

  // ExternalFixedRecordSorter owns two descriptors per merge input and four
  // fixed descriptors. The predecessor set may be sorted while the current
  // set has one read descriptor open. Account those on top of every final scan
  // reader, writer, and mutable disk array; this turns --max-open-files into a
  // real ceiling rather than an advisory fan-in value.
  const size_type fixed_descriptors = reader_streams + writer_streams + 2;
  const size_type sorter_fixed_descriptors = 5;
  if(parameters.getMaxOpenFiles() < fixed_descriptors +
     sorter_fixed_descriptors + 2 * 2)
  {
    throw std::runtime_error("GCSA::GCSA(): max-open-files cannot hold the final-scan streams");
  }
  size_type node_set_fan_in = std::min(parameters.getMergeFanIn(),
    (parameters.getMaxOpenFiles() - fixed_descriptors -
      sorter_fixed_descriptors) / 2);
  node_set_fan_in = std::max(static_cast<size_type>(2), node_set_fan_in);

  std::string previous_name = TempFile::getName("gcsa_final_prev_occ");
  std::string stack_name = TempFile::getName("gcsa_final_lcp_stack");
  try
  {
    {
      SpillableNodeSet curr_from(node_set_budget, node_set_fan_in, memory);
      SpillableNodeSet pred_from(node_set_budget, node_set_fan_in, memory);

      // Every ReadBuffer may simultaneously own its foreground window and its
      // asynchronous refill vector. Admit both windows only after the two
      // persistent start-node workspaces, so later writer reservations cannot
      // wait on memory that this single-threaded phase will never release.
      size_type reader_buffer = finalScanBuffer(parameters.getIOBufferSize(),
        minimum_reader, memory.available(), 4 * reader_streams);
      size_type reader_reservation_bytes = checkedProduct(2 * reader_streams,
        reader_buffer, "final reader reservation");
      MemoryBudget::Reservation reader_reservation = memory.reserve(
        reader_reservation_bytes, "final-merged-graph-readers");

      size_type writer_buffer = finalScanBuffer(parameters.getIOBufferSize(),
        static_cast<size_type>(16), memory.available(), 2 * writer_streams);
      if(checkedProduct(writer_streams, writer_buffer,
          "final writer reservation") > memory.available())
      {
        throw std::runtime_error("GCSA::GCSA(): memory limit cannot hold final event buffers");
      }
      FinalEventWriter output(files, graph.alpha.sigma, writer_buffer, memory,
        parameters.getTempFileCodecParameters());

      // Predecessor ranges are independent by alphabet component, while mapped
      // current-node sets and their immutable from-ranks are independent by
      // path. Admit both fixed records and node/rank arenas before handing the
      // remaining bytes to the mutable arrays. Lookup and sampling own separate
      // monotone streams; the source-label reader replaces two streams omitted
      // from the ordered reader, keeping reader_streams unchanged. A tight
      // budget or small configured I/O buffer keeps the original serial scan as
      // an explicit fallback.
      size_type predecessor_batch_capacity = 0;
      size_type prepared_from_capacity = 0;
      MemoryBudget::Reservation predecessor_work_reservation;
      std::vector<FinalPredecessorWork> predecessor_work;
      std::vector<node_type> prepared_from_nodes;
      std::vector<size_type> prepared_from_ranks;
      if(requested_predecessor_workers > 1 && memory.available() > array_minimum)
      {
        const size_type staging_ceiling =
          (memory.available() - array_minimum) / 2;
        const size_type desired_batch_capacity =
          parameters.getIOBufferSize() / sizeof(FinalPredecessorWork);
        const size_type desired_work_bytes = checkedProduct(
          desired_batch_capacity, sizeof(FinalPredecessorWork),
          "final predecessor work");
        const size_type desired_node_bytes = checkedProduct(
          parameters.getIOBufferSize(), 2, "final from-node arenas");
        if(desired_work_bytes >
           std::numeric_limits<size_type>::max() - desired_node_bytes)
        {
          throw std::runtime_error(
            "GCSA::GCSA(): final preparation staging size overflows");
        }
        const size_type staging_bytes = std::min(staging_ceiling,
          desired_work_bytes + desired_node_bytes);
        const size_type bytes_per_minimum_path = sizeof(FinalPredecessorWork) +
          sizeof(node_type) + sizeof(size_type);
        predecessor_batch_capacity = std::min(desired_batch_capacity,
          staging_bytes / bytes_per_minimum_path);
        if(predecessor_batch_capacity >= MIN_FINAL_PREDECESSOR_BATCH)
        {
          const size_type work_bytes = checkedProduct(predecessor_batch_capacity,
            sizeof(FinalPredecessorWork), "final predecessor batch");
          prepared_from_capacity = (staging_bytes - work_bytes) /
            (sizeof(node_type) + sizeof(size_type));
          if(prepared_from_capacity < predecessor_batch_capacity)
          {
            throw std::runtime_error(
              "GCSA::GCSA(): final preparation arena cannot hold one node per path");
          }
          const size_type from_node_bytes = checkedProduct(prepared_from_capacity,
            sizeof(node_type), "final prepared from nodes");
          const size_type from_rank_bytes = checkedProduct(prepared_from_capacity,
            sizeof(size_type), "final prepared from ranks");
          const size_type admitted_bytes = work_bytes + from_node_bytes +
            from_rank_bytes;
          predecessor_work_reservation = memory.reserve(admitted_bytes,
            "final-unordered-preparation");
          predecessor_work.resize(predecessor_batch_capacity);
          prepared_from_nodes.resize(prepared_from_capacity);
          prepared_from_ranks.resize(prepared_from_capacity);
          if(stats != nullptr)
          {
            stats->predecessor_workers = requested_predecessor_workers;
            stats->predecessor_buffer_bytes = work_bytes;
            stats->from_preparation_buffer_bytes =
              from_node_bytes + from_rank_bytes;
          }
        }
        else
        {
          predecessor_batch_capacity = prepared_from_capacity = 0;
        }
      }
      if(stats != nullptr && requested_predecessor_workers > 1 &&
         predecessor_batch_capacity == 0)
      {
        stats->predecessor_parallel_fallback = true;
      }

      size_type cache_bytes = memory.available();
      size_type previous_cache = cache_bytes / 2;
      size_type stack_cache = cache_bytes - previous_cache;
      if(unique_from_nodes > 0 && previous_cache < DiskBackedArray64::minimumCacheBytes())
      {
        throw std::runtime_error("GCSA::GCSA(): memory limit cannot hold previous-occurrence cache");
      }
      if(merged_graph.size() > 0 && stack_cache < DiskBackedArray64::minimumCacheBytes())
      {
        throw std::runtime_error("GCSA::GCSA(): memory limit cannot hold suffix-tree stack cache");
      }
      DiskBackedArray64 previous(previous_name, unique_from_nodes,
        previous_cache, memory, true, 64 * KILOBYTE);
      // RETIRED, NOT REMOVED: the whole-graph capacity below.
      //
      // The array was sized at 3 * merged_graph.size() elements -- 9.59 GB for
      // chr21 -- and DiskBackedArray64 zero-fills its cache, so it reserved
      // 8.94 GiB against the budget during the one sub-step that already runs
      // pinned at the cgroup limit. The traversal cannot use it: it pushes
      // frames with strictly increasing curr_lcp popped against lcp_array,
      // which is a ReadBuffer<uint8_t>, so the stack holds at most one frame
      // per distinct LCP value and the highest index it can touch is
      // 3 * 255 + 2. The chr21 log confirms it empirically -- "ST stack 1
      // reads / 1 writes", one 64 KiB block across 399,778,113 paths.
      //
      // Set WHOLE_GRAPH_SUFFIX_TREE_STACK to true to restore the old sizing if
      // a future LCP representation widens past a byte and this bound is ever
      // suspected; the bound below throws rather than overruns, so the failure
      // would be loud either way.
      constexpr bool WHOLE_GRAPH_SUFFIX_TREE_STACK = false;
      // Derive the bound from the LCP stream's own element type rather than a
      // literal, so widening that stream widens this with it. MergedGraph
      // writes the array as SequentialRecordWriter<merged_lcp_type> and the
      // scan below reads it as ReadBuffer<merged_lcp_type>; curr_lcp is one
      // such value plus one, so the distinct values are [0, max + 1] and the
      // stack holds one frame of three elements for each.
      typedef std::uint8_t merged_lcp_type;
      constexpr size_type distinct_lcp_values =
        static_cast<size_type>(std::numeric_limits<merged_lcp_type>::max()) + 2;
      size_type stack_capacity = (WHOLE_GRAPH_SUFFIX_TREE_STACK ?
        checkedProduct(merged_graph.size(), 3, "suffix-tree stack") :
        checkedProduct(distinct_lcp_values, 3, "suffix-tree stack"));
      DiskBackedArray64 stack(stack_name, stack_capacity,
        stack_cache, memory, true, 64 * KILOBYTE);

      const size_type scan_read_start = DiskIO::read_volume;
      std::vector<MergedGraphReader> reader(graph.alpha.sigma + 1);
      std::vector<MergedGraphReader> sampling_reader(
        predecessor_batch_capacity > 0 ? graph.alpha.sigma : 0);
      if(predecessor_batch_capacity > 0)
      {
        reader[0].initFromOnly(merged_graph, reader_buffer, true);
      }
      else
      {
        reader[0].init(merged_graph, &mapper, &last_char, reader_buffer, true);
      }
      for(size_type comp = 0; comp < graph.alpha.sigma; comp++)
      {
        reader[comp + 1].init(merged_graph, comp, reader_buffer, true,
          predecessor_batch_capacity == 0);
        if(predecessor_batch_capacity > 0)
        {
          sampling_reader[comp].initFromOnly(merged_graph, reader_buffer, true, comp);
        }
      }
      MergedGraphLabelReader predecessor_source;
      if(predecessor_batch_capacity > 0)
      {
        predecessor_source.init(merged_graph, reader_buffer);
      }
      ReadBuffer<merged_lcp_type> lcp_array;
      lcp_array.open(merged_graph.lcp_name, reader_buffer, true);

      PathLabel first, last;
      size_type stack_size = 0;
      ProgressReporter scan_progress("final event scan",
        merged_graph.size(), "paths");
      std::array<std::exception_ptr, FinalEventMetadata::MAX_SIGMA>
        predecessor_errors;
      for(size_type i = 0; i < merged_graph.size();)
      {
        size_type batch_paths = 1;
        if(predecessor_batch_capacity > 0)
        {
          const size_type batch_limit = std::min(predecessor_batch_capacity,
            merged_graph.size() - i);
          size_type prepared_from_used = 0;
          batch_paths = 0;
          curr_from.clear();
          for(size_type offset = 0; offset < batch_limit; offset++)
          {
            FinalPredecessorWork& work = predecessor_work[offset];
            work.path = predecessor_source.paths[predecessor_source.path];
            if(work.path.ranks() > work.labels.size())
            {
              throw std::runtime_error(
                "GCSA::GCSA(): final predecessor label exceeds the fixed bound");
            }
            for(size_type rank = 0; rank < work.path.ranks(); rank++)
            {
              work.labels[rank] =
                predecessor_source.labels[work.path.pointer() + rank];
            }
            work.destinations.fill(MergedGraph::UNKNOWN);
            work.destination_primary.fill(MergedGraph::UNKNOWN);
            work.from_offset = prepared_from_used;
            work.from_count = 0;
            work.from_spilled = false;
            work.sample_by_period = false;

            if(reader[0].path != i + offset)
            {
              throw std::runtime_error(
                "GCSA::GCSA(): current from-node reader lost batch alignment");
            }
            const auto append_from = [&](node_type unmapped) {
              node_type node = mappedFinalNode(unmapped, graph.mapping);
              if(!(work.from_spilled) &&
                 prepared_from_used < prepared_from_capacity)
              {
                prepared_from_nodes[prepared_from_used++] = node;
                work.from_count++;
                return;
              }
              if(!(work.from_spilled))
              {
                work.from_spilled = true;
                curr_from.clear();
                for(size_type j = 0; j < work.from_count; j++)
                {
                  node_type buffered =
                    prepared_from_nodes[work.from_offset + j];
                  curr_from.push_back(buffered);
                  if(buffered % parameters.getSamplePeriod() == 0)
                  {
                    work.sample_by_period = true;
                  }
                }
              }
              curr_from.push_back(node);
              if(node % parameters.getSamplePeriod() == 0)
              {
                work.sample_by_period = true;
              }
            };

            append_from(work.path.from);
            while(reader[0].from < reader[0].from_nodes.size() &&
                  reader[0].from_nodes[reader[0].from].first == reader[0].path)
            {
              append_from(reader[0].from_nodes[reader[0].from].second);
              reader[0].from++;
            }
            reader[0].from_nodes.seek(reader[0].from);
            if(work.from_spilled)
            {
              curr_from.finish();
              work.from_count = curr_from.size();
              if(stats != nullptr) { stats->prepared_from_spill_fallbacks++; }
            }

            predecessor_source.advance();
            reader[0].advance();
            batch_paths++;
            // curr_from is the one persistent overflow slot. Ending the batch
            // here lets the ordered consumer reach this path before the slot is
            // reused, without rereading any prefix of the from-node stream.
            if(work.from_spilled) { break; }
          }

          predecessor_errors.fill(std::exception_ptr());
          std::exception_ptr from_preparation_error;
          std::atomic<bool> batch_failed(false);
#pragma omp parallel num_threads(requested_predecessor_workers)
          {
#pragma omp for schedule(dynamic, 64) nowait
            for(std::int64_t item = 0;
                item < static_cast<std::int64_t>(batch_paths); item++)
            {
              size_type offset = static_cast<size_type>(item);
              FinalPredecessorWork& work = predecessor_work[offset];
              if(work.from_spilled) { continue; }
              try
              {
                auto begin = prepared_from_nodes.begin() + work.from_offset;
                auto end = begin + work.from_count;
                std::sort(begin, end);
                end = std::unique(begin, end);
                work.from_count = static_cast<size_type>(end - begin);
                for(size_type j = 0; j < work.from_count; j++)
                {
                  node_type node = prepared_from_nodes[work.from_offset + j];
                  prepared_from_ranks[work.from_offset + j] =
                    finalFromRank(node, from_nodes);
                  if(node % parameters.getSamplePeriod() == 0)
                  {
                    work.sample_by_period = true;
                  }
                }
              }
              catch(...)
              {
#pragma omp critical(gcsa_final_from_preparation_error)
                {
                  if(!(from_preparation_error))
                  {
                    from_preparation_error = std::current_exception();
                  }
                }
                batch_failed.store(true, std::memory_order_relaxed);
              }
            }

#pragma omp for schedule(static, 1) nowait
            for(std::int64_t component = 0;
                component < static_cast<std::int64_t>(graph.alpha.sigma);
                component++)
            {
              size_type comp = static_cast<size_type>(component);
              try
              {
                MergedGraphReader& destination_reader = reader[comp + 1];
                PathLabel worker_first, worker_last;
                for(size_type offset = 0; offset < batch_paths; offset++)
                {
                  if(batch_failed.load(std::memory_order_relaxed)) { break; }
                  FinalPredecessorWork& work = predecessor_work[offset];
                  if(!(work.path.hasPredecessor(comp))) { continue; }
                  predecessorRange(work.path, comp, mapper, last_char,
                    [&work](size_type rank) { return work.labels[rank]; },
                    worker_first, worker_last);
                  if(!(destination_reader.intersect(worker_first, worker_last, 0)))
                  {
                    destination_reader.advance();
                  }
                  if(destination_reader.path >= merged_graph.size())
                  {
                    throw std::runtime_error(
                      "GCSA::GCSA(): final edge destination is outside the merged graph");
                  }
                  work.destinations[comp] = destination_reader.path;
                  work.destination_primary[comp] =
                    destination_reader.paths[destination_reader.path].from;
                }
              }
              catch(...)
              {
                predecessor_errors[comp] = std::current_exception();
                batch_failed.store(true, std::memory_order_relaxed);
              }
            }
          }
          if(from_preparation_error)
          {
            std::rethrow_exception(from_preparation_error);
          }
          for(size_type comp = 0; comp < graph.alpha.sigma; comp++)
          {
            if(predecessor_errors[comp])
            {
              std::rethrow_exception(predecessor_errors[comp]);
            }
          }
          if(stats != nullptr)
          {
            stats->predecessor_batches++;
            stats->maximum_predecessor_batch = std::max(
              stats->maximum_predecessor_batch, batch_paths);
          }
        }

        for(size_type batch_offset = 0; batch_offset < batch_paths;
            batch_offset++, i++)
        {
          scan_progress.advance();
          // Close any predecessor spill reader retained by the preceding path
          // before current-set collection can invoke the external sorter.
          pred_from.clear();
          size_type indegree = 0, pred_comp = 0;
          byte_type predecessor_mask = 0;
          bool sample_this = false;
          const PathNode& source_path = (predecessor_batch_capacity > 0 ?
            predecessor_work[batch_offset].path :
            reader[0].paths[reader[0].path]);
          for(size_type comp = 0; comp < graph.alpha.sigma; comp++)
          {
            if(!(source_path.hasPredecessor(comp))) { continue; }
            size_type destination = MergedGraph::UNKNOWN;
            if(predecessor_batch_capacity > 0)
            {
              destination = predecessor_work[batch_offset].destinations[comp];
            }
            else
            {
              reader[0].predecessor(comp, first, last);
              if(!(reader[comp + 1].intersect(first, last, 0)))
              {
                reader[comp + 1].advance();
              }
              destination = reader[comp + 1].path;
            }
            if(destination >= merged_graph.size())
            {
              throw std::runtime_error("GCSA::GCSA(): final edge destination is outside the merged graph");
            }
            predecessor_mask |= static_cast<byte_type>(static_cast<size_type>(1) << comp);
            output.edge(comp, destination);
            indegree++; pred_comp = comp;
          }
          output.path(predecessor_mask);

          size_type curr_from_size = 0;
          if(predecessor_batch_capacity == 0)
          {
            reader[0].fromNodes(curr_from, graph.mapping);
            curr_from_size = curr_from.size();
          }
          else
          {
            curr_from_size = predecessor_work[batch_offset].from_count;
          }
          if(curr_from_size == 0)
          {
            throw std::runtime_error("GCSA::GCSA(): merged path has no start node");
          }
          if(stats != nullptr)
          {
            bool current_spilled = (predecessor_batch_capacity == 0 ?
              curr_from.spilled() :
              predecessor_work[batch_offset].from_spilled);
            stats->from_node_spills += current_spilled;
            stats->maximum_from_nodes = std::max(stats->maximum_from_nodes,
              curr_from_size);
            if(predecessor_batch_capacity > 0 && !current_spilled)
            {
              stats->prepared_from_paths++;
              stats->prepared_from_ranks += curr_from_size;
            }
          }
          output.occurrence(i, curr_from_size - 1);

          lcp_array.seek(i);
          size_type curr_lcp = lcp_array[i] + (i > 0 ? 1 : 0);
          while(stack_size > 0 && stack.get(3 * (stack_size - 1)) > curr_lcp)
          {
            stack_size--;
          }
          if(stack_size > 0 && stack.get(3 * (stack_size - 1)) == curr_lcp)
          {
            stack.set(3 * (stack_size - 1) + 2, i);
          }
          else
          {
            // The bounded capacity above assumes one frame per distinct LCP
            // value. Fail loudly rather than write past the array if a future
            // LCP representation ever breaks that assumption.
            if(3 * stack_size + 2 >= stack_capacity)
            {
              throw std::runtime_error(
                "GCSA::GCSA(): suffix-tree stack exceeded its bounded capacity");
            }
            stack.set(3 * stack_size, curr_lcp);
            stack.set(3 * stack_size + 1, i);
            stack.set(3 * stack_size + 2, i);
            stack_size++;
          }

          // Rank resolution is immutable and was prepared in parallel for
          // bounded in-memory slices. The spill/serial route resolves it here.
          // Only this callback touches previous occurrence state or emits a
          // redundancy event, so path order remains the exact serial order.
          const auto update_occurrence = [&](size_type rank) {
            size_type prior = previous.get(rank);
            if(prior > 0)
            {
              size_type low = 0, high = stack_size;
              while(low < high)
              {
                size_type middle = low + (high - low) / 2;
                if(stack.get(3 * middle + 2) < prior) { low = middle + 1; }
                else { high = middle; }
              }
              if(low >= stack_size)
              {
                throw std::runtime_error("GCSA::GCSA(): invalid previous-occurrence suffix-tree state");
              }
              size_type first_time = stack.get(3 * low + 1);
              if(first_time == 0)
              {
                throw std::runtime_error("GCSA::GCSA(): redundancy position underflows");
              }
              output.redundancy(first_time - 1);
            }
            previous.set(rank, i + 1);
          };

          node_type node;
          if(predecessor_batch_capacity > 0 &&
             !(predecessor_work[batch_offset].from_spilled))
          {
            const FinalPredecessorWork& work =
              predecessor_work[batch_offset];
            for(size_type j = 0; j < work.from_count; j++)
            {
              update_occurrence(prepared_from_ranks[work.from_offset + j]);
            }
          }
          else
          {
            curr_from.rewind();
            while(curr_from.next(node))
            {
              update_occurrence(finalFromRank(node, from_nodes));
            }
          }

          if(indegree > 1) { sample_this = true; }
          if(source_path.hasPredecessor(Alphabet::SINK_COMP))
          {
            sample_this = true;
          }
          if(predecessor_batch_capacity > 0)
          {
            sample_this = sample_this ||
              predecessor_work[batch_offset].sample_by_period;
          }
          else
          {
            curr_from.rewind();
            while(curr_from.next(node))
            {
              if(node % parameters.getSamplePeriod() == 0)
              {
                sample_this = true; break;
              }
            }
          }

          if(!sample_this)
          {
            if(indegree == 0)
            {
              throw std::runtime_error("GCSA::GCSA(): unsampled path has no predecessor");
            }
            if(predecessor_batch_capacity > 0)
            {
              MergedGraphReader& sample_reader = sampling_reader[pred_comp];
              size_type destination =
                predecessor_work[batch_offset].destinations[pred_comp];
              if(destination < sample_reader.path)
              {
                throw std::runtime_error("GCSA::GCSA(): sampling destination moved backwards");
              }
              sample_reader.path = destination;
              sample_reader.seek(false);
              sample_reader.fromNodes(pred_from, graph.mapping,
                predecessor_work[batch_offset].destination_primary[pred_comp]);
            }
            else { reader[pred_comp + 1].fromNodes(pred_from, graph.mapping); }
            if(stats != nullptr)
            {
              stats->from_node_spills += pred_from.spilled();
              stats->maximum_from_nodes = std::max(stats->maximum_from_nodes,
                pred_from.size());
            }
            if(pred_from.size() != curr_from_size) { sample_this = true; }
            else
            {
              pred_from.rewind();
              if(predecessor_batch_capacity > 0)
              {
                const FinalPredecessorWork& work =
                  predecessor_work[batch_offset];
                bool continuation = true;
                forEachFinalFromNode(work, prepared_from_nodes, curr_from,
                  [&pred_from, &continuation](node_type curr_node) {
                    node_type pred_node;
                    if(continuation &&
                       (!pred_from.next(pred_node) || curr_node != pred_node + 1))
                    {
                      continuation = false;
                    }
                  });
                if(!continuation) { sample_this = true; }
              }
              else
              {
                node_type curr_node, pred_node;
                curr_from.rewind();
                while(curr_from.next(curr_node))
                {
                  if(!pred_from.next(pred_node) || curr_node != pred_node + 1)
                  {
                    sample_this = true; break;
                  }
                }
              }
            }
          }

          if(sample_this)
          {
            output.sampledPath(i);
            if(predecessor_batch_capacity > 0)
            {
              forEachFinalFromNode(predecessor_work[batch_offset],
                prepared_from_nodes, curr_from,
                [&output](node_type sample) { output.sample(sample); });
            }
            else
            {
              curr_from.rewind();
              while(curr_from.next(node)) { output.sample(node); }
            }
            output.sampleEnd();
          }
          if(predecessor_batch_capacity == 0) { reader[0].advance(); }
        }
      }
      predecessor_source.close();
      for(MergedGraphReader& current : reader) { current.close(); }
      for(MergedGraphReader& current : sampling_reader) { current.close(); }
      lcp_array.close();
      scan_progress.finish();
      if(Verbosity::level >= Verbosity::EXTENDED)
      {
        std::cerr << "GCSA::GCSA(): final scan read bytes: "
                  << (DiskIO::read_volume - scan_read_start) << std::endl;
      }
      previous.flush(false); stack.flush(false);
      if(stats != nullptr)
      {
        stats->previous_occurrences = previous.stats();
        stats->suffix_tree_stack = stack.stats();
      }
      metadata = output.finish();
      event_checksums = output.checksums();
    }

    TempFile::remove(previous_name); TempFile::remove(stack_name);
    metadata.fast_chars = graph.alpha.fast_chars;
    if(metadata.paths != merged_graph.size())
    {
      throw std::runtime_error("GCSA::GCSA(): final event path count mismatch");
    }
    // The redundancy stream is committed unsorted: component construction
    // needs only its per-slot counts.
    writeFinalEventMetadata(files, metadata);
    checkpointFinalEvents(workspace, files, metadata, checkpoint_buffer,
      &event_checksums);
  }
  catch(...)
  {
    TempFile::remove(previous_name); TempFile::remove(stack_name); throw;
  }
  return metadata;
}

} // namespace

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

// " (reason)" for a declined parallel prune or merge, or nothing.
std::string
fallbackReason(const char* reason)
{
  return (reason == nullptr ? std::string() : std::string(" (") + reason + ")");
}

/*
  FORK: reporting for the external route.

  These were 120 lines inline in GCSA::GCSA, where they were 22% of a 672-line
  constructor and hid its control flow. Each one only reads a statistics struct
  and writes to std::cerr, so lifting them changes nothing and lets the
  constructor read as a sequence of phases again. Each keeps its own verbosity
  guard so the call site is a single unconditional line.
*/

void
reportPruneMergeStats(const PathGraphMergeStats& merge_stats)
{
  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    std::cerr << "PathGraph::prune(): "
              << merge_stats.prune_workers << " of "
              << merge_stats.prune_requested_workers << " requested worker(s) over "
              << merge_stats.prune_partitions << " root partition(s), "
              << merge_stats.prune_parallel_fallbacks << " serial fallback(s)"
              << fallbackReason(merge_stats.prune_fallback_reason) << "; "
              << merge_stats.priority_spills << " path-group spills, "
              << merge_stats.range_spills << " range spills, "
              << merge_stats.max_open_input_pairs << " input pairs and "
              << merge_stats.max_open_output_pairs << " output pairs open at peak; "
              << merge_stats.path_input_refills << "/"
              << merge_stats.rank_input_refills << " path/rank window refills for "
              << merge_stats.path_input_reads << " paths using "
              << formatBytes(merge_stats.max_input_buffer_bytes)
              << " of bounded input windows"
              << std::endl;
    if(merge_stats.prefetch_workers > 0)
    {
      std::cerr << "PathGraph::prune(): framed prefetch: "
                << merge_stats.prefetch_workers << " workers, "
                << merge_stats.prefetch_consumed << "/"
                << merge_stats.prefetch_submitted << " blocks consumed/submitted ("
                << merge_stats.prefetch_completed << " completed, "
                << merge_stats.prefetch_ready_hits << " ready, "
                << merge_stats.prefetch_waits << " waited, "
                << merge_stats.prefetch_synchronous_blocks << " synchronous, "
                << merge_stats.prefetch_cancelled << " cancelled, "
                << merge_stats.prefetch_errors << " errors), "
                << formatBytes(merge_stats.prefetch_physical_bytes)
                << " stored bytes requested -> "
                << formatBytes(merge_stats.prefetch_decoded_bytes)
                << " successfully decoded work, "
                << (static_cast<double>(merge_stats.prefetch_wait_nanoseconds) / 1.0e9)
                << " seconds waiting, peak "
                << formatBytes(merge_stats.max_prefetch_bytes)
                << std::endl;
    }
  }
}


void
reportJoinStats(const ExternalPathJoinStats& join_stats)
{
  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    std::cerr << "externalPathGraphExtend(): "
              << join_stats.left_records << " left, "
              << join_stats.right_records << " right, "
              << join_stats.generated_records << " generated, "
              << join_stats.sorted_bypass << " bypassed" << std::endl;
    std::cerr << "externalPathGraphExtend(): "
              << join_stats.initial_runs << " join runs, "
              << join_stats.merge_operations << " join merges, "
              << join_stats.join_partitions << " join partitions, "
              << join_stats.worker_processes << " worker processes, "
              << join_stats.restored_partitions << " restored partitions, "
              << join_stats.recursive_splits << " recursive splits, "
              << join_stats.label_sort_runs << " label runs, "
              << join_stats.label_merge_passes << " label merge passes, "
              << join_stats.blocked_key_groups << " blocked key groups" << std::endl;
    std::cerr << "externalPathGraphExtend(): "
              << join_stats.join_parallel_sorts << " parallel join sorts, "
              << join_stats.label_parallel_sorts << " parallel label sorts, "
              << join_stats.compaction_batches << " compaction batches ("
              << join_stats.compaction_concurrency << " at once)" << std::endl;
    std::cerr << "externalPathGraphExtend(): distribution in "
              << join_stats.distribution_ranges << " key ranges ("
              << join_stats.distribution_concurrency << " at once, "
              << join_stats.restored_range_plans << " range plans restored)" << std::endl;
    std::cerr << "externalPathGraphExtend(): phase budgets: distribution "
              << formatBytes(join_stats.distribution_sort_budget)
              << "; concurrent label sort "
              << formatBytes(join_stats.label_sort_budget) << " + join block "
              << formatBytes(join_stats.join_block_budget) << std::endl;
    std::cerr << "externalPathGraphExtend(): sampled "
              << join_stats.sampled_plan_records << " join keys into "
              << join_stats.radix_plan_bins << " MSD range packs ("
              << join_stats.radix_plan_splits << " radix splits, "
              << join_stats.radix_plan_capped << " budget-capped plans, "
              << join_stats.radix_boundary_flushes << " exact boundary flushes, "
              << join_stats.restored_radix_plans << " restored plans)" << std::endl;
    std::cerr << "externalPathGraphExtend(): streamed "
              << join_stats.direct_label_records << " records directly to label runs; avoided "
              << formatBytes(join_stats.intermediate_path_bytes_avoided)
              << " of intermediate path/rank I/O per direction" << std::endl;
    std::cerr << "externalPathGraphExtend(): emitted "
              << join_stats.grouped_expansion_records
              << " compact left-context references across transient run passes; avoided "
              << formatBytes(join_stats.expansion_context_bytes_saved)
              << " of transient run payload" << std::endl;
    if(join_stats.join_run_logical_bytes > 0)
    {
      double ratio = static_cast<double>(join_stats.join_run_stored_bytes) /
        static_cast<double>(join_stats.join_run_logical_bytes);
      std::cerr << "externalPathGraphExtend(): stored "
                << formatBytes(join_stats.join_run_stored_bytes) << " for "
                << formatBytes(join_stats.join_run_logical_bytes)
                << " of logical fixed-record join runs ("
                << join_stats.compressed_join_runs << " compressed runs, "
                << ratio << " stored/logical)" << std::endl;
    }
    if(join_stats.join_sidecar_logical_bytes > 0)
    {
      double ratio = static_cast<double>(join_stats.join_sidecar_stored_bytes) /
        static_cast<double>(join_stats.join_sidecar_logical_bytes);
      std::cerr << "externalPathGraphExtend(): stored "
                << formatBytes(join_stats.join_sidecar_stored_bytes) << " for "
                << formatBytes(join_stats.join_sidecar_logical_bytes)
                << " of logical join group/detail sidecars ("
                << join_stats.compressed_join_sidecars << " compressed sidecars, "
                << ratio << " stored/logical)" << std::endl;
    }
    std::cerr << "externalPathGraphExtend(): maximum bounded workspace "
              << formatBytes(join_stats.max_bytes_resident) << " ("
              << join_stats.max_records_resident << " records)" << std::endl;
  }
}


void
reportFinalMergeStats(const PathGraphMergeStats& final_merge_stats)
{
  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    std::cerr << "GCSA::GCSA(): LCP range minima: "
              << LCP::range_minimum_queries.load() << " wavelet descents over "
              << LCP::range_minimum_span.load() << " total key positions"
              << std::endl;
    std::cerr << "MergedGraph: "
              << final_merge_stats.merge_workers << " worker(s)"
              << fallbackReason(final_merge_stats.merge_fallback_reason) << " over "
              << final_merge_stats.merge_partitions << " root partition(s); "
              << final_merge_stats.priority_spills << " path-group spills, "
              << final_merge_stats.range_spills << " range spills, "
              << final_merge_stats.from_set_sorts << " external from-set sorts (peak "
              << final_merge_stats.max_from_set_nodes << " resident nodes) and "
              << final_merge_stats.max_open_input_pairs << " input pairs open at peak; "
              << final_merge_stats.path_input_refills << "/"
              << final_merge_stats.rank_input_refills
              << " path/rank window refills for "
              << final_merge_stats.path_input_reads << " paths using "
              << formatBytes(final_merge_stats.max_input_buffer_bytes)
              << " of bounded input windows"
              << std::endl;
    if(final_merge_stats.prefetch_workers > 0)
    {
      std::cerr << "MergedGraph: framed prefetch: "
                << final_merge_stats.prefetch_workers << " workers, "
                << final_merge_stats.prefetch_consumed << "/"
                << final_merge_stats.prefetch_submitted
                << " blocks consumed/submitted ("
                << final_merge_stats.prefetch_completed << " completed, "
                << final_merge_stats.prefetch_ready_hits << " ready, "
                << final_merge_stats.prefetch_waits << " waited, "
                << final_merge_stats.prefetch_synchronous_blocks
                << " synchronous, "
                << final_merge_stats.prefetch_cancelled << " cancelled, "
                << final_merge_stats.prefetch_errors << " errors), "
                << formatBytes(final_merge_stats.prefetch_physical_bytes)
                << " stored bytes requested -> "
                << formatBytes(final_merge_stats.prefetch_decoded_bytes)
                << " successfully decoded work, "
                << (static_cast<double>(final_merge_stats.prefetch_wait_nanoseconds) /
                    1.0e9)
                << " seconds waiting, peak "
                << formatBytes(final_merge_stats.max_prefetch_bytes)
                << std::endl;
    }
  }
}


void
reportFinalEventStats(const FinalEventMetadata& event_metadata,
  const ExternalFinalScanStats& event_stats)
{
  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    std::cerr << "GCSA::GCSA(): final events "
              << (event_stats.restored ? "restored" : "produced") << ", "
              << event_metadata.total_edges << " edges, "
              << event_metadata.sample_ids << " sample ids, "
              << event_metadata.occurrence_items << " nonzero occurrence events, "
              << event_metadata.redundant << " redundancy events" << std::endl;
    if(!event_stats.restored)
    {
      std::cerr << "GCSA::GCSA(): disk caches: prev-occ "
                << event_stats.previous_occurrences.block_reads << " reads / "
                << event_stats.previous_occurrences.block_writes << " writes; ST stack "
                << event_stats.suffix_tree_stack.block_reads << " reads / "
                << event_stats.suffix_tree_stack.block_writes << " writes" << std::endl;
      std::cerr << "GCSA::GCSA(): final start-node sets: "
                << event_stats.from_node_spills << " spills, maximum "
                << event_stats.maximum_from_nodes << " unique nodes" << std::endl;
      std::cerr << "GCSA::GCSA(): final predecessor lookup: "
                << event_stats.predecessor_workers << " worker(s), "
                << event_stats.predecessor_batches << " bounded batch(es), maximum "
                << event_stats.maximum_predecessor_batch << " paths / "
                << formatBytes(event_stats.predecessor_buffer_bytes)
                << " admitted";
      if(event_stats.predecessor_parallel_fallback)
      {
        std::cerr << " (serial fallback: batch admission too small)";
      }
      std::cerr << std::endl;
      std::cerr << "GCSA::GCSA(): final ordered-state preparation: "
                << event_stats.predecessor_workers << " worker(s), "
                << event_stats.prepared_from_paths << " paths / "
                << event_stats.prepared_from_ranks << " from-ranks prepared, "
                << event_stats.prepared_from_spill_fallbacks
                << " spill fallback(s), "
                << formatBytes(event_stats.from_preparation_buffer_bytes)
                << " admitted" << std::endl;
    }
  }
}

//------------------------------------------------------------------------------

GCSA::GCSA(InputGraph& graph, const ConstructionParameters& parameters) :
  GCSA(graph, parameters, nullptr)
{
}

void
GCSA::buildAndStore(InputGraph& graph, const ConstructionParameters& parameters,
  const std::string& filename)
{
  if(!parameters.externalMemory())
  {
    throw std::invalid_argument("GCSA::buildAndStore() requires a durable external-memory workspace");
  }
  if(graph.size() == 0)
  {
    GCSA empty;
    storeEmptyIndexAtomically(empty, filename);
    return;
  }
  GCSA builder(graph, parameters, &filename);
}

GCSA::GCSA(InputGraph& graph, const ConstructionParameters& parameters,
  const std::string* direct_output) :
  GCSA()
{
  double start = readTimer();

  if(graph.size() == 0) { return; }
  size_type bytes_required = graph.size() * (sizeof(PathNode) + 2 * sizeof(PathNode::rank_type));
  if(bytes_required > parameters.getLimitBytes())
  {
    std::cerr << "GCSA::GCSA(): The input is too large: " << (bytes_required / GIGABYTE_DOUBLE) << " GB" << std::endl;
    std::cerr << "GCSA::GCSA(): Construction aborted" << std::endl;
    std::exit(EXIT_SIZE_LIMIT_EXCEEDED);
  }

  size_type checkpoint_buffer = std::max(static_cast<size_type>(1),
    std::min(parameters.getIOBufferSize(), parameters.getMemoryLimitBytes() / 16));
  std::unique_ptr<BuildWorkspace> workspace;
  if(parameters.externalMemory())
  {
    BuildWorkspace::Settings semantic = constructionSemanticSettings(
      graph, parameters, checkpoint_buffer);
    BuildWorkspace::Settings operational = constructionOperationalSettings(parameters);
    workspace.reset(new BuildWorkspace(parameters.getWorkDirectory(), semantic, operational,
      parameters.getResume() ? BuildWorkspace::RESUME : BuildWorkspace::NEW_WORKSPACE));
    if(Verbosity::level >= Verbosity::BASIC)
    {
      // This is the aggregate allocator/scheduler goal used to choose spill
      // runs and admit join workers. Library allocations and the embedding
      // application's graph remain outside it, so deployments should report
      // their hard cgroup cap separately.
      std::cerr << "GCSA::GCSA(): External working-set target "
                << formatBytes(parameters.getMemoryLimitBytes())
                << ", disk limit " << formatBytes(parameters.getLimitBytes())
                << std::endl;
    }
  }

  // Extract key and start-node facts. The legacy branch retains its historical
  // vectors exactly. With a workspace, ExternalInputPreprocessor instead
  // streams immutable record files through bounded sort/merge passes and
  // checkpoints the reduced key/start streams for resume. Only the LCP support
  // is needed during prefix doubling; mapper, last-character, and start-node
  // supports are deliberately delayed so they do not inflate the longest and
  // largest spill phase's resident set.
  DeBruijnGraph mapper;
  LCP lcp;
  sdsl::int_vector<0> last_char;
  sdsl::sd_vector<> from_nodes;
  size_type unique_from_nodes = 0;
  std::unique_ptr<ExternalInputPreprocessor> external_preprocessor;
  sdsl::int_vector<0> distinct_labels;
  if(parameters.externalMemory())
  {
    external_preprocessor.reset(new ExternalInputPreprocessor(graph, parameters,
      workspace.get()));
    external_preprocessor->buildLCP(lcp);
  }
  else
  {
    std::vector<key_type> keys;
    graph.readKeys(keys);
    mapper = DeBruijnGraph(keys, graph.k(), graph.alpha);
    lcp = LCP(keys, graph.k());
    Key::lastChars(keys, last_char);
    distinct_labels = sdsl::int_vector<0>(keys.size(), 0,
      bit_length(Key::label(keys.back())));
    for(size_type i = 0; i < keys.size(); i++) { distinct_labels[i] = Key::label(keys[i]); }
    sdsl::util::clear(keys);

    // Determine the existing start nodes. Because the information is only used for index
    // construction, we use NodeMapping to map the node ids used for construction to the
    // node ids reported by locate().
    std::vector<node_type> from_node_buffer;
    graph.readFrom(from_node_buffer, true);
    from_nodes = sdsl::sd_vector<>(from_node_buffer.begin(), from_node_buffer.end());
    unique_from_nodes = from_node_buffer.size();
    sdsl::util::clear(from_node_buffer);
  }
  sdsl::sd_vector<>::rank_1_type from_rank;

  // Create or restore the initial PathGraph. Resumption still rebuilds the
  // compact key/LCP support above, but never repeats a committed doubling phase.
  PathGraph path_graph(0, graph.k(), 0);
  const auto checkpoint_path_graph = [&](const std::string& task,
      const std::string& phase) {
    const std::uint64_t scan_before =
      BuildWorkspace::adoption_checksum_scan_bytes.load(
        std::memory_order_relaxed);
    const std::uint64_t reused_before =
      BuildWorkspace::adoption_checksum_reused_bytes.load(
        std::memory_order_relaxed);
    const double checkpoint_start = readTimer();
    checkpointPathGraph(*workspace, path_graph, task, phase,
      checkpoint_buffer);
    const double checkpoint_stop = readTimer();
    const std::uint64_t scan_after =
      BuildWorkspace::adoption_checksum_scan_bytes.load(
        std::memory_order_relaxed);
    const std::uint64_t reused_after =
      BuildWorkspace::adoption_checksum_reused_bytes.load(
        std::memory_order_relaxed);
    if(Verbosity::level >= Verbosity::BASIC)
    {
      std::cerr << "GCSA::GCSA(): Checkpoint " << task << "/" << phase
                << ": " << (checkpoint_stop - checkpoint_start)
                << " seconds, adoption checksum scanned "
                << formatBytes(scan_after - scan_before) << ", reused "
                << formatBytes(reused_after - reused_before)
                << " of writer provenance" << std::endl;
    }
  };
  size_type first_step = 1;
  bool restored_prune = false, restored_graph = false;
  // Tracks the durable checkpoint backing the current in-memory frontier.
  // --clean-obsolete may retire it only after the next frontier has its own
  // synced task marker. Operational cleanup therefore cannot invalidate a
  // graph still needed after a crash.
  std::string frontier_task, frontier_phase;
  if(workspace)
  {
    for(size_type step = parameters.getSteps(); step > 0; step--)
    {
      std::string task = doublingTask(step);
      if(pathGraphCheckpointExists(*workspace, task, "extend"))
      {
        if(Verbosity::level >= Verbosity::BASIC)
        {
          std::cerr << "GCSA::GCSA(): Restoring checkpoint "
                    << task << "/extend" << std::endl;
        }
        restorePathGraph(*workspace, path_graph, task, "extend", checkpoint_buffer,
          parameters.getVerifyWorkspace());
        frontier_task = task; frontier_phase = "extend";
        if(parameters.getCleanObsolete())
        {
          retireBeforeFrontier(*workspace, task, "extend", step, true);
        }
        first_step = step + 1; restored_graph = true; break;
      }
      if(pathGraphCheckpointExists(*workspace, task, "prune"))
      {
        if(Verbosity::level >= Verbosity::BASIC)
        {
          std::cerr << "GCSA::GCSA(): Restoring checkpoint "
                    << task << "/prune" << std::endl;
        }
        restorePathGraph(*workspace, path_graph, task, "prune", checkpoint_buffer,
          parameters.getVerifyWorkspace());
        frontier_task = task; frontier_phase = "prune";
        if(parameters.getCleanObsolete())
        {
          retireBeforeFrontier(*workspace, task, "prune", step, false);
        }
        first_step = step; restored_prune = true; restored_graph = true; break;
      }
    }
    if(!restored_graph && pathGraphCheckpointExists(*workspace, "initial", "paths"))
    {
      if(Verbosity::level >= Verbosity::BASIC)
      {
        std::cerr << "GCSA::GCSA(): Restoring checkpoint initial/paths" << std::endl;
      }
      restorePathGraph(*workspace, path_graph, "initial", "paths", checkpoint_buffer,
        parameters.getVerifyWorkspace());
      frontier_task = "initial"; frontier_phase = "paths";
      restored_graph = true;
    }
  }
  if(!restored_graph)
  {
    PathGraph initial_graph(0, graph.k(), 0);
    if(parameters.externalMemory())
    {
      external_preprocessor->buildInitialPathGraph(initial_graph);
    }
    else
    {
      PathGraph legacy_initial(graph, distinct_labels);
      initial_graph.swap(legacy_initial);
    }
    path_graph.swap(initial_graph);
    if(workspace)
    {
      checkpoint_path_graph("initial", "paths");
      frontier_task = "initial"; frontier_phase = "paths";
      stopAfterCommittedPhase(parameters, "initial");
    }
  }
  if(!parameters.externalMemory()) { sdsl::util::clear(distinct_labels); }
  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    double stop = readTimer();
    std::cerr << "GCSA::GCSA(): Preprocessing: " << (stop - start) << " seconds, "
              << inGigabytes(memoryUsage()) << " GB" << std::endl;
    start = stop;
  }

  // Prefix-doubling.
  if(Verbosity::level >= Verbosity::BASIC)
  {
    std::cerr << "GCSA::GCSA(): Prefix-doubling from path length " << path_graph.k() << std::endl;
  }
  for(size_type step = first_step; step <= parameters.getSteps(); step++)
  {
    if(Verbosity::level >= Verbosity::BASIC)
    {
      std::cerr << "GCSA::GCSA(): Step " << step << " (path length " << path_graph.k() << " -> "
                << (2 * path_graph.k()) << ")" << std::endl;
    }
    std::string task = doublingTask(step);
    if(!(restored_prune && step == first_step))
    {
      size_type prune_buffer = pathMergeInputBudget(parameters, path_graph);
      PathGraphMergeStats merge_stats;
      const size_type merge_inputs = path_graph.files();
      // Parallel prune may retain several physical shards for one logical
      // input. The external extender joins across those shards by logical id;
      // legacy PathGraph::extend() still joins only within each physical file.
      // Keep that historical route serial until it has the same shard-aware
      // contract, regardless of the caller's OpenMP thread count.
      const size_type prune_workers = (parameters.externalMemory() ?
        static_cast<size_type>(std::max(1, omp_get_max_threads())) : 1);
      const double prune_start = readTimer();
      path_graph.prune(lcp, path_graph.remainingLimit(parameters.getLimitBytes()),
        prune_buffer, &merge_stats, parameters.getMaxOpenFiles(),
        pathMergeInputCacheBudget(parameters, path_graph),
        prune_workers, parameters.getConcurrentOpenFiles());
      const double prune_stop = readTimer();
      if(Verbosity::level >= Verbosity::BASIC)
      {
        std::cerr << "GCSA::GCSA(): Prune " << task << ": "
                  << (prune_stop - prune_start) << " seconds; "
                  << merge_stats.prune_workers << " of "
                  << merge_stats.prune_requested_workers
                  << " requested worker(s), "
                  << merge_stats.prune_partitions << " partition(s), "
                  << merge_stats.prune_parallel_fallbacks
                  << " serial fallback(s)"
                  << fallbackReason(merge_stats.prune_fallback_reason) << std::endl;
      }
      // The merger visits its shards in round-robin label order, so a cache one
      // entry short of the shard count misses on every record rather than on a
      // fraction of them. That cliff is worth a word at the default verbosity:
      // it is invisible in the output and the remedy is one flag.
      if(merge_stats.prune_workers == 1 &&
         merge_stats.max_open_input_pairs < merge_inputs &&
         Verbosity::level >= Verbosity::BASIC)
      {
        std::cerr << "PathGraph::prune(): warning: only "
                  << merge_stats.max_open_input_pairs << " of " << merge_inputs
                  << " input shards fit --max-open-files "
                  << parameters.getMaxOpenFiles()
                  << "; the merge will reread the rest. Raise --max-open-files"
                  << " to about " << (4 * merge_inputs + 2) << "." << std::endl;
      }
      reportPruneMergeStats(merge_stats);
      if(workspace)
      {
        checkpoint_path_graph(task, "prune");
        if(parameters.getCleanObsolete() && !frontier_task.empty())
        {
          workspace->retire_obsolete(frontier_task, frontier_phase, task, "prune");
        }
        frontier_task = task; frontier_phase = "prune";
        stopAfterCommittedPhase(parameters, task + "-prune");
      }
    }
    restored_prune = false;
    if(parameters.externalMemory())
    {
      ExternalPathJoinStats join_stats;
      externalPathGraphExtend(path_graph,
        path_graph.remainingLimit(parameters.getLimitBytes()), parameters, &join_stats,
        workspace.get(), task);
      reportJoinStats(join_stats);
    }
    else
    {
      path_graph.extend(path_graph.remainingLimit(parameters.getLimitBytes()),
        parameters.getMemoryLimitBytes());
    }
    if(workspace)
    {
      checkpoint_path_graph(task, "extend");
      if(parameters.getCleanObsolete() && !frontier_task.empty())
      {
        workspace->retire_obsolete(frontier_task, frontier_phase, task, "extend");
      }
      if(parameters.getCleanObsolete())
      {
        // The full generation checkpoint now supersedes fine-grained join
        // outputs and its sampled MSD plan. A resume repeats this family scan,
        // so a crash between individual retirement markers only delays space
        // reclamation; it cannot lose the committed frontier.
        workspace->retire_obsolete_family(task + "-join-", "join-partition",
          task, "extend");
        workspace->retire_obsolete_family(task + "-msd-plan-", "join-plan",
          task, "extend");
        workspace->retire_obsolete_family(task + "-range-plan-", "range-plan",
          task, "extend");
      }
      frontier_task = task; frontier_phase = "extend";
      stopAfterCommittedPhase(parameters, task + "-extend");
    }
  }
  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    double stop = readTimer();
    std::cerr << "GCSA::GCSA(): Prefix-doubling: " << (stop - start) << " seconds, "
              << inGigabytes(memoryUsage()) << " GB" << std::endl;
    start = stop;
  }

  // Merge the paths into the nodes of a maximally pruned de Bruijn graph.
  if(Verbosity::level >= Verbosity::BASIC)
  {
    std::cerr << "GCSA::GCSA(): Merging the paths" << std::endl;
  }
  if(parameters.externalMemory())
  {
    // The mapper is read-only and only participates in the final merge and
    // event scan. Build it after doubling rather than carrying it through all
    // generated path generations.
    SubPhaseProbe mapper_probe("merge/mapper");
    external_preprocessor->buildMapper(mapper);
    mapper_probe.report();
  }
  size_type merge_buffer = pathMergeInputBudget(parameters, path_graph);
  size_type merge_cache = pathMergeInputCacheBudget(parameters, path_graph);
  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    // The two numbers that decide the merge's entire input working set.
    std::cerr << "GCSA::GCSA(): merge group buffer " << inMegabytes(merge_buffer)
              << " MB, framed input cache " << inMegabytes(merge_cache)
              << " MB over " << path_graph.files() << " shard(s)" << std::endl;
  }
  // The partitioned merge copies every partition into place, doubling the
  // merge's scratch writes, and its chr21 gain did not clear run-to-run
  // spread; it runs only when requested.
  const char* parallel_merge = std::getenv("GCSA_EXPERIMENTAL_PARALLEL_MERGE");
  size_type merge_workers = 1;
  if(parameters.externalMemory() && parallel_merge != nullptr && std::string(parallel_merge) == "1")
  {
    merge_workers = static_cast<size_type>(std::max(1, omp_get_max_threads()));
  }
  PathGraphMergeStats final_merge_stats;
  SubPhaseProbe merge_probe("merge/merged-graph");
  // Opt-in (GCSA_IO_COMPRESS_MERGE=1, with temporary compression on): the
  // serial merge writes its paths, ranks and start nodes framed, compressed on
  // one thread per stream. On chr18 that cut the merge's writes from 25.8 to
  // 4.1 GiB and the final scan's device reads from 72.7 to 15.3 GiB, but the
  // CPU-bound scan took 0:01:27 longer, so it stays off by default.
  const TempFileCodecParameters merge_codec = parameters.getTempFileCodecParameters();
  const char* compress_merge = std::getenv("GCSA_IO_COMPRESS_MERGE");
  const bool framed_merge = parameters.externalMemory() && compress_merge != nullptr &&
    std::string(compress_merge) == "1";
  MergedGraph merged_graph(path_graph, mapper, lcp,
    path_graph.remainingLimit(parameters.getLimitBytes()), merge_buffer,
    &final_merge_stats, parameters.getMaxOpenFiles(), merge_cache, merge_workers,
    (framed_merge ? &merge_codec : nullptr));
  merge_probe.report();
  reportFinalMergeStats(final_merge_stats);
  this->header.path_nodes = merged_graph.size();
  this->header.order = merged_graph.k();
  path_graph.clear();
  sdsl::util::clear(lcp);
  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    double stop = readTimer();
    std::cerr << "GCSA::GCSA(): Merging: " << (stop - start) << " seconds, "
              << inGigabytes(memoryUsage()) << " GB" << std::endl;
    start = stop;
  }

  // Structures used for building GCSA.
  if(Verbosity::level >= Verbosity::BASIC)
  {
    std::cerr << "GCSA::GCSA(): Building the index" << std::endl;
  }
  size_type occ_count = 0, red_count = 0;
  size_type final_sample_count = 0, final_sampled_positions = 0;
  if(parameters.externalMemory())
  {
    // These supports are only used by the final event scan. Delaying them until
    // after PathGraph/LCP release avoids retaining four key/start-proportional
    // structures during prefix doubling. The external scan takes its start-node
    // ranks from sd_vector::successor, so it needs no rank support at all --
    // only the legacy route below builds one.
    external_preprocessor->buildLastCharacters(last_char);
    external_preprocessor->buildStartNodes(from_nodes);
    unique_from_nodes = external_preprocessor->startNodeCount();

    FinalEventFiles event_files(graph.alpha.sigma);
    ExternalFinalScanStats event_stats;
    SubPhaseProbe events_probe("construct/final-event-scan");
    FinalEventMetadata event_metadata = produceExternalFinalEvents(merged_graph,
      mapper, last_char, from_nodes, unique_from_nodes, graph,
      parameters, *workspace, event_files, checkpoint_buffer, &event_stats);
    events_probe.report();
    stopAfterCommittedPhase(parameters, "final-events");

    // No final component needs the construction mapper or mutable scan state.
    // Releasing them before component assembly is an important part of the
    // working-set reduction.
    sdsl::util::clear(last_char); sdsl::util::clear(from_nodes);
    sdsl::util::clear(mapper);
    this->header.edges = event_metadata.total_edges;
    if(direct_output != nullptr)
    {
      SubPhaseProbe store_probe("construct/store-components");
      storeFinalComponents(this->header, graph.alpha, event_files,
        event_metadata, parameters, *direct_output);
      store_probe.report();
    }
    else
    {
      // The external route publishes its components straight to the output
      // file; it has never had a resident assembler that ships. Every
      // production caller reaches it through GCSA::buildAndStore, so a null
      // target here is a programming error rather than a configuration.
      throw std::invalid_argument(
        "GCSA::GCSA(): external-memory construction requires a publication"
        " target; use GCSA::buildAndStore()");
    }
    if(event_metadata.occurrence_extra >
       std::numeric_limits<size_type>::max() - event_metadata.paths)
    {
      throw std::runtime_error("GCSA::GCSA(): occurrence count overflows");
    }
    occ_count = event_metadata.paths + event_metadata.occurrence_extra;
    red_count = event_metadata.redundant;
    final_sample_count = event_metadata.sample_ids;
    final_sampled_positions = event_metadata.sampled_paths;

    reportFinalEventStats(event_metadata, event_stats);
  }
  else
  {
  sdsl::util::init_support(from_rank, &(from_nodes));
  sdsl::int_vector<64> counts(graph.alpha.sigma, 0); // alpha
  std::vector<bit_vector> bwt(graph.alpha.sigma); // fast_bwt, sparse_bwt
  for(size_type comp = 0; comp < bwt.size(); comp++) { bwt[comp] = bit_vector(merged_graph.size(), 0); }
  CounterArray outdegrees(merged_graph.size(), 4); // edges
  bit_vector sampled_positions(merged_graph.size(), 0); // sampled_paths
  std::vector<node_type> sample_buffer; // stored_samples
  this->samples = bit_vector(merged_graph.size() + merged_graph.extra(), 0);

  // Structures used for building counting support.
  // Invariant: The previous occurrence of from node x was at path prev_occ[from_rank(x)] - 1.
  CounterArray occurrences(merged_graph.size(), 4), redundant(merged_graph.size() - 1, 4);
  sdsl::int_vector<0> prev_occ(unique_from_nodes, 0, bit_length(merged_graph.size()));
  std::vector<size_type> node_lcp, first_time, last_time;

  // Read pointers to the MergedGraph files.
  std::vector<MergedGraphReader> reader(graph.alpha.sigma + 1);
  reader[0].init(merged_graph, &mapper, &last_char);
  for(size_type comp = 0; comp < graph.alpha.sigma; comp++)
  {
    reader[comp + 1].init(merged_graph, comp);
  }
  ReadBuffer<uint8_t> lcp_array; lcp_array.open(merged_graph.lcp_name);

  // The actual construction.
  PathLabel first, last;
  size_type total_edges = 0, sample_bits = 0;
  std::vector<node_type> pred_from, curr_from;
  for(size_type i = 0; i < merged_graph.size(); i++, reader[0].advance())
  {
    // Find the predecessors.
    size_type indegree = 0, pred_comp = 0;
    bool sample_this = false;
    for(size_type comp = 0; comp < graph.alpha.sigma; comp++)
    {
      if(!(reader[0].paths[reader[0].path].hasPredecessor(comp))) { continue; }

      // Find the predecessor of paths[i] with comp and the path intersecting it.
      reader[0].predecessor(comp, first, last);
      if(!(reader[comp + 1].intersect(first, last, 0)))
      {
        reader[comp + 1].advance();
      }

      // Add the edge.
      bwt[comp][i] = 1; counts[comp]++;
      indegree++; outdegrees.increment(reader[comp + 1].path); total_edges++;
      pred_comp = comp; // For sampling.
    }

    /*
      Get the start nodes and update the occurrences/redundant arrays.

      We traverse the ST in inorder using the LCP array. For each internal node, we
      record the LCP value and the first and the last times (positions) we have
      encountered that value within the subtree. If we have encountered the current
      from node before, the LCA of the previous and current occurrences is the
      highest ST node we have encountered after the previous occurrence. We then
      increment the redundant array at the first encounter with that node.
    */
    reader[0].fromNodes(curr_from, graph.mapping);
    occurrences.increment(i, curr_from.size() - 1);
    lcp_array.seek(i);
    size_type curr_lcp = lcp_array[i] + (i > 0 ? 1 : 0); // Handle LCP[0] as -1.
    while(!(node_lcp.empty()) && node_lcp.back() > curr_lcp)
    {
      node_lcp.pop_back(); first_time.pop_back(); last_time.pop_back();
    }
    if(!(node_lcp.empty()) && node_lcp.back() == curr_lcp) { last_time.back() = i; }
    else { node_lcp.push_back(curr_lcp); first_time.push_back(i); last_time.push_back(i); }
    for(size_type j = 0; j < curr_from.size(); j++)
    {
      size_type temp = from_rank(curr_from[j]);
      if(prev_occ[temp] > 0)
      {
        size_type pos = std::lower_bound(last_time.begin(), last_time.end(), prev_occ[temp]) - last_time.begin();
        redundant.increment(first_time[pos] - 1);
      }
      prev_occ[temp] = i + 1;
    }

    /*
      Simple cases for sampling the node:
      - multiple predecessors
      - at the beginning of the source node with no real predecessors
      - if a from node is divisible by the sample period
    */
    if(indegree > 1) { sample_this = true; }
    if(reader[0].paths[reader[0].path].hasPredecessor(Alphabet::SINK_COMP)) { sample_this = true; }
    for(size_type k = 0; k < curr_from.size(); k++)
    {
      if(curr_from[k] % parameters.getSamplePeriod() == 0) { sample_this = true; break; }
    }

    // Sample if the start nodes cannot be derived from the only predecessor.
    if(!sample_this)
    {
      reader[pred_comp + 1].fromNodes(pred_from, graph.mapping);
      if(pred_from.size() != curr_from.size()) { sample_this = true; }
      else
      {
        for(size_type k = 0; k < curr_from.size(); k++)
        {
          if(curr_from[k] != pred_from[k] + 1) { sample_this = true; break; }
        }
      }
    }

    // Store the samples.
    if(sample_this)
    {
      sampled_positions[i] = 1;
      for(size_type k = 0; k < curr_from.size(); k++)
      {
        sample_bits = std::max(sample_bits, bit_length(curr_from[k]));
        sample_buffer.push_back(curr_from[k]);
      }
      this->samples[sample_buffer.size() - 1] = 1;
    }
  }
  for(size_type i = 0; i < reader.size(); i++) { reader[i].close(); }
  lcp_array.close();
  sdsl::util::clear(last_char); sdsl::util::clear(from_nodes); sdsl::util::clear(prev_occ);
  this->header.edges = total_edges;

  // Initialize alpha.
  this->alpha = Alphabet(counts, graph.alpha.char2comp, graph.alpha.comp2char);
  sdsl::util::clear(mapper);

  // Initialize extra_pointers and redundant_pointers.
  occ_count = occurrences.sum() + occurrences.size(); red_count = redundant.sum();
  this->extra_pointers = SadaSparse(occurrences);
  this->redundant_pointers = SadaCount(redundant);
  sdsl::util::clear(occurrences); sdsl::util::clear(redundant);

  // Initialize bwt.
  this->fast_bwt.resize(this->alpha.sigma); this->fast_rank.resize(this->alpha.sigma);
  this->sparse_bwt.resize(this->alpha.sigma); this->sparse_rank.resize(this->alpha.sigma);
  this->sparse_bwt[0] = bwt[0]; sdsl::util::clear(bwt[0]);
  for(size_type comp = 1; comp <= graph.alpha.fast_chars; comp++)
  {
    this->fast_bwt[comp] = bwt[comp]; sdsl::util::clear(bwt[comp]);
  }
  for(size_type comp = graph.alpha.fast_chars + 1; comp < graph.alpha.sigma; comp++)
  {
    this->sparse_bwt[comp] = bwt[comp]; sdsl::util::clear(bwt[comp]);
  }

  // Initialize bitvectors (edges, sampled_positions, samples).
  bit_vector edge_buffer(total_edges, 0); total_edges = 0;
  for(size_type i = 0; i < merged_graph.size(); i++)
  {
    total_edges += outdegrees[i];
    edge_buffer[total_edges - 1] = 1;
  }
  outdegrees.clear();
  this->edges = edge_buffer; sdsl::util::clear(edge_buffer);
  this->sampled_paths = sampled_positions; sdsl::util::clear(sampled_positions);
  this->samples.resize(sample_buffer.size());
  this->initSupport();

  // Initialize stored_samples.
  this->stored_samples = sdsl::int_vector<0>(sample_buffer.size(), 0, sample_bits);
  for(size_type i = 0; i < sample_buffer.size(); i++) { this->stored_samples[i] = sample_buffer[i]; }
  sdsl::util::clear(sample_buffer);
  final_sample_count = this->sampleCount();
  final_sampled_positions = this->sampledPositions();
  }

  // Transfer the LCP array from MergedGraph to InputGraph.
  TempFile::remove(graph.lcp_name);
  graph.lcp_name = merged_graph.lcp_name;
  merged_graph.lcp_name.clear();

  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    double stop = readTimer();
    std::cerr << "GCSA::GCSA(): Construction: " << (stop - start) << " seconds, "
              << inGigabytes(memoryUsage()) << " GB" << std::endl;
  }
  if(Verbosity::level >= Verbosity::BASIC)
  {
    std::cerr << "GCSA::GCSA(): " << this->size() << " paths, " << this->edgeCount() << " edges" << std::endl;
    std::cerr << "GCSA::GCSA(): " << occ_count << " pointers (" << red_count << " redundant)" << std::endl;
    std::cerr << "GCSA::GCSA(): " << final_sample_count << " samples at "
              << final_sampled_positions << " positions" << std::endl;
  }
}

void
GCSA::initSupport()
{
  for(size_type comp = 0; comp < this->alpha.sigma; comp++)
  {
    sdsl::util::init_support(this->fast_rank[comp], &(this->fast_bwt[comp]));
    sdsl::util::init_support(this->sparse_rank[comp], &(this->sparse_bwt[comp]));
  }

  sdsl::util::init_support(this->edge_rank, &(this->edges));
  sdsl::util::init_support(this->sampled_path_rank, &(this->sampled_paths));
  sdsl::util::init_support(this->sample_select, &(this->samples));
}

//------------------------------------------------------------------------------

void
GCSA::LF_fast(range_type range, std::vector<range_type>& results) const
{
  for(size_type comp = 1; comp <= this->alpha.fast_chars; comp++) { results[comp] = Range::empty_range(); }
  if(Range::empty(range)) { return; }

  if(range.first == range.second) // Single path node.
  {
    for(size_type comp = 1; comp <= this->alpha.fast_chars; comp++)
    {
      if(this->fast_bwt[comp][range.first])
      {
        results[comp].first = results[comp].second = this->edge_rank(this->LF(this->fast_rank, range.first, comp));
      }
    }
  }
  else  // General case.
  {
    for(size_type comp = 1; comp <= this->alpha.fast_chars; comp++)
    {
      results[comp] = this->LF(this->fast_rank, range, comp);
      if(!Range::empty(results[comp])) { results[comp] = this->pathNodeRange(results[comp]); }
    }
  }
}

void
GCSA::LF_all(range_type range, std::vector<range_type>& results) const
{
  for(size_type comp = 1; comp + 1 < this->alpha.sigma; comp++) { results[comp] = Range::empty_range(); }
  if(Range::empty(range)) { return; }

  if(range.first == range.second) // Single path node.
  {
    for(size_type comp = 1; comp <= this->alpha.fast_chars; comp++)
    {
      if(this->fast_bwt[comp][range.first])
      {
        results[comp].first = results[comp].second = this->edge_rank(this->LF(this->fast_rank, range.first, comp));
      }
    }
    for(size_type comp = this->alpha.fast_chars + 1; comp + 1 < this->alpha.sigma; comp++)
    {
      if(this->sparse_bwt[comp][range.first])
      {
        results[comp].first = results[comp].second = this->edge_rank(this->LF(this->sparse_rank, range.first, comp));
      }
    }
  }
  else  // General case.
  {
    for(size_type comp = 1; comp + 1 < this->alpha.sigma; comp++)
    {
      results[comp] = this->LF(range, comp);
    }
  }
}

//------------------------------------------------------------------------------

size_type
GCSA::count(range_type range) const
{
  if(Range::empty(range) || range.second >= this->size()) { return 0; }
  size_type res = this->extra_pointers.count(range.first, range.second) + Range::length(range);
  if(range.second > range.first) { res -= this->redundant_pointers.count(range.first, range.second - 1); }
  return res;
}

//------------------------------------------------------------------------------

void
GCSA::locate(size_type path_node, std::vector<node_type>& results, bool append, bool sort) const
{
  if(!append) { results.clear(); }
  if(path_node >= this->size())
  {
    if(sort) { removeDuplicates(results, false); }
    return;
  }

  this->locateInternal(path_node, results);
  if(sort) { removeDuplicates(results, false); }
}

void
GCSA::locate(range_type range, std::vector<node_type>& results, bool append, bool sort) const
{
  if(!append) { results.clear(); }
  if(Range::empty(range) || range.second >= this->size())
  {
    if(sort) { removeDuplicates(results, false); }
    return;
  }

  for(size_type i = range.first; i <= range.second; i++)
  {
    this->locateInternal(i, results);
  }
  if(sort) { removeDuplicates(results, false); }
}

void
GCSA::locate(range_type range, const OccurrenceCallback& report) const
{
  if(!report) { throw std::invalid_argument("GCSA::locate(): missing occurrence callback"); }
  if(Range::empty(range) || range.second >= this->size()) { return; }
  for(size_type path_node = range.first; path_node <= range.second; path_node++)
  {
    size_type steps = 0, sample_path = path_node;
    while(!(this->sampled(sample_path)))
    {
      sample_path = this->LF(sample_path); steps++;
    }
    size_type sample = this->firstSample(sample_path);
    do
    {
      report(this->sample(sample) + steps); sample++;
    }
    while(!(this->lastSample(sample - 1)));
  }
}

void
GCSA::locate(range_type range, size_type max_positions, std::vector<node_type>& results) const
{
  results.clear();

  size_type total_positions = this->count(range);
  if(total_positions <= 0) { return; }
  max_positions = std::min(max_positions, total_positions);

  std::mt19937_64 rng(range.first ^ range.second);
  if(max_positions >= total_positions / 2)  // Just locate everything.
  {
    this->locate(range, results);
  }
  else  // Locate at random positions until we have enough distinct occurrences.
  {
    std::unordered_set<node_type, size_type(*)(size_type)> found(16, wang_hash_64);
    while(found.size() < max_positions)
    {
      size_type pos = range.first + rng() % Range::length(range);
      this->locateInternal(pos, results);
      for(node_type node : results) { found.insert(node); }
      results.clear();
    }
    for(node_type node : found) { results.push_back(node); }
  }

  // If there are too many results, select a random subset of them.
  if(results.size() > max_positions)
  {
    deterministicShuffle(results, rng, false);
    results.resize(max_positions);
  }
  sequentialSort(results.begin(), results.end());
}

void
GCSA::locateInternal(size_type path_node, std::vector<node_type>& results) const
{
  size_type steps = 0;
  while(!(this->sampled(path_node)))
  {
    path_node = this->LF(path_node);
    steps++;
  }

  size_type sample = this->firstSample(path_node);
  do
  {
    results.push_back(this->sample(sample) + steps); sample++;
  }
  while(!(this->lastSample(sample - 1)));
}

//------------------------------------------------------------------------------

} // namespace gcsa
