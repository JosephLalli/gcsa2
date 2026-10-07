#include <gcsa/algorithms.h>
#include <external_configuration.hpp>
#include <gcsa/path_graph_external.h>
#include <gcsa/disk_array.h>
#include <gcsa/external_preprocessing.h>
#include <gcsa/final_events.h>
#include <gcsa/internal.h>
#include <gcsa/lcp.h>
#include <gcsa/path_graph.h>

#include <algorithm>
#include <cstdlib>
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
formatBytes(size_type bytes)
{
  const char* suffixes[] = { "B", "KiB", "MiB", "GiB", "TiB", "PiB" };
  long double value = bytes;
  size_type suffix = 0;
  while(value >= 1024.0L && suffix + 1 < sizeof(suffixes) / sizeof(suffixes[0]))
  {
    value /= 1024.0L; suffix++;
  }
  std::ostringstream out;
  out << std::fixed << std::setprecision(value < 10.0L && suffix > 0 ? 2 : 1)
      << static_cast<double>(value) << ' ' << suffixes[suffix];
  return out.str();
}

std::string
doublingTask(size_type step)
{
  std::ostringstream result;
  result << "step-" << std::setw(2) << std::setfill('0') << step;
  return result.str();
}

class SerialExternalConstruction
{
public:
  SerialExternalConstruction() : thread_count(omp_get_max_threads())
  {
    omp_set_num_threads(1);
  }

  ~SerialExternalConstruction()
  {
    omp_set_num_threads(this->thread_count);
  }

private:
  int thread_count;
};

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
    bool release_cache = false);
  void close();

  void seek();

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
  size_type buffer_bytes, bool release_cache)
{
  this->paths.open(graph.path_name, buffer_bytes, release_cache);
  this->labels.open(graph.rank_name, buffer_bytes, release_cache);
  this->from_nodes.open(graph.from_name, buffer_bytes, release_cache);

  this->path = graph.next[comp]; this->path_count = graph.size();
  this->from = graph.next_from[comp];
  this->seek();

  this->mapper = nullptr;
  this->last_char = nullptr;
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
MergedGraphReader::seek()
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
  if(this->paths.descriptor >= 0)
  {
    this->paths.seek(this->path);
    this->rank = this->paths[this->path].pointer();
    this->labels.seek(this->rank);
  }
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
  results.clear();
  const auto mapped = [&mapping](node_type node) {
    return (mapping.empty() ? node :
      Node::encode(mapping(Node::id(node)), Node::offset(node), Node::rc(node)));
  };
  results.push_back(mapped(this->paths[this->path].from));

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

size_type
checkedProduct(size_type first, size_type second, const std::string& label)
{
  if(second != 0 && first > std::numeric_limits<size_type>::max() / second)
  {
    throw std::runtime_error("GCSA::GCSA(): " + label + " size overflows");
  }
  return first * second;
}

// Use the target when a fair stream share can hold it. Under a tight budget,
// keep the same share calculation while respecting the caller's minimum.
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
  Scan the final MergedGraph once and publish compact immutable event streams.
  Prev-occ and the suffix-tree traversal stack remain bounded disk arrays for
  the lifetime of this fresh scan.
*/
FinalEventMetadata
produceExternalFinalEvents(const MergedGraph& merged_graph,
  const DeBruijnGraph& mapper, const sdsl::int_vector<0>& last_char,
  const sdsl::sd_vector<>& from_nodes,
  size_type unique_from_nodes, const InputGraph& graph,
  const ConstructionParameters& parameters, FinalEventFiles& files)
{
  FinalEventMetadata metadata;

  if(graph.alpha.sigma == 0 || graph.alpha.sigma > FinalEventMetadata::MAX_SIGMA)
  {
    throw std::runtime_error("GCSA::GCSA(): external final events require an alphabet of at most 8 components");
  }
  const size_type memory_limit = parameters.getMemoryLimitBytes();
  const size_type safety_margin = memory_limit / 8;
  MemoryBudget memory(memory_limit, safety_margin);

  const size_type reader_streams = 3 * (graph.alpha.sigma + 1);
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
    std::min(externalIOBufferSize(parameters), externalSortRunSize(parameters)));
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
  if(externalMaxOpenFiles() < fixed_descriptors +
     sorter_fixed_descriptors + 2 * 2)
  {
    throw std::runtime_error("GCSA::GCSA(): max-open-files cannot hold the final-scan streams");
  }
  size_type node_set_fan_in = std::min(externalMergeFanIn(),
    (externalMaxOpenFiles() - fixed_descriptors -
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
      size_type reader_buffer = finalScanBuffer(externalIOBufferSize(parameters),
        minimum_reader, memory.available(), 4 * reader_streams);
      size_type reader_reservation_bytes = checkedProduct(2 * reader_streams,
        reader_buffer, "final reader reservation");
      MemoryBudget::Reservation reader_reservation = memory.reserve(
        reader_reservation_bytes);

      size_type writer_buffer = finalScanBuffer(externalIOBufferSize(parameters),
        static_cast<size_type>(16), memory.available(), 2 * writer_streams);
      if(checkedProduct(writer_streams, writer_buffer,
          "final writer reservation") > memory.available())
      {
        throw std::runtime_error("GCSA::GCSA(): memory limit cannot hold final event buffers");
      }
      FinalEventWriter output(files, graph.alpha.sigma, writer_buffer, memory);

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
      // Derive the bound from the LCP stream's own element type rather than a
      // literal, so widening that stream widens this with it. MergedGraph
      // writes the array as SequentialRecordWriter<merged_lcp_type> and the
      // scan below reads it as ReadBuffer<merged_lcp_type>; curr_lcp is one
      // such value plus one, so the distinct values are [0, max + 1] and the
      // stack holds one frame of three elements for each.
      typedef std::uint8_t merged_lcp_type;
      constexpr size_type distinct_lcp_values =
        static_cast<size_type>(std::numeric_limits<merged_lcp_type>::max()) + 2;
      size_type stack_capacity = checkedProduct(distinct_lcp_values, 3,
        "suffix-tree stack");
      DiskBackedArray64 stack(stack_name, stack_capacity,
        stack_cache, memory, true, 64 * KILOBYTE);

      std::vector<MergedGraphReader> reader(graph.alpha.sigma + 1);
      reader[0].init(merged_graph, &mapper, &last_char, reader_buffer, true);
      for(size_type comp = 0; comp < graph.alpha.sigma; comp++)
      {
        reader[comp + 1].init(merged_graph, comp, reader_buffer, true);
      }
      ReadBuffer<merged_lcp_type> lcp_array;
      lcp_array.open(merged_graph.lcp_name, reader_buffer, true);

      PathLabel first, last;
      size_type stack_size = 0;
      for(size_type i = 0; i < merged_graph.size(); i++, reader[0].advance())
      {
        // Close any predecessor spill reader retained by the preceding path
        // before current-set collection can invoke the external sorter.
        pred_from.clear();
        size_type indegree = 0, pred_comp = 0;
        byte_type predecessor_mask = 0;
        bool sample_this = false;
        for(size_type comp = 0; comp < graph.alpha.sigma; comp++)
        {
          if(!(reader[0].paths[reader[0].path].hasPredecessor(comp))) { continue; }
          reader[0].predecessor(comp, first, last);
          if(!(reader[comp + 1].intersect(first, last, 0)))
          {
            reader[comp + 1].advance();
          }
          if(reader[comp + 1].path >= merged_graph.size())
          {
            throw std::runtime_error("GCSA::GCSA(): final edge destination is outside the merged graph");
          }
          predecessor_mask |= static_cast<byte_type>(static_cast<size_type>(1) << comp);
          output.edge(comp, reader[comp + 1].path);
          indegree++; pred_comp = comp;
        }
        output.path(predecessor_mask);

        reader[0].fromNodes(curr_from, graph.mapping);
        if(curr_from.size() == 0)
        {
          throw std::runtime_error("GCSA::GCSA(): merged path has no start node");
        }
        output.occurrence(i, curr_from.size() - 1);

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

        node_type node;
        curr_from.rewind();
        while(curr_from.next(node))
        {
          // sd_vector::operator[] performs its own select_0 and backward scan,
          // and rank_support_sd::rank performs another, so testing membership
          // and then taking the rank costs two descents for one lookup.
          // successor(i) yields both: its one_iterator dereferences to
          // (rank, position). ValueIndex::find uses the same idiom.
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
          size_type rank = occurrence->first;
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
        }

        if(indegree > 1) { sample_this = true; }
        if(reader[0].paths[reader[0].path].hasPredecessor(Alphabet::SINK_COMP))
        {
          sample_this = true;
        }
        curr_from.rewind();
        while(curr_from.next(node))
        {
          if(node % parameters.getSamplePeriod() == 0)
          {
            sample_this = true; break;
          }
        }

        if(!sample_this)
        {
          if(indegree == 0)
          {
            throw std::runtime_error("GCSA::GCSA(): unsampled path has no predecessor");
          }
          reader[pred_comp + 1].fromNodes(pred_from, graph.mapping);
          if(pred_from.size() != curr_from.size()) { sample_this = true; }
          else
          {
            node_type curr_node, pred_node;
            curr_from.rewind(); pred_from.rewind();
            while(curr_from.next(curr_node))
            {
              if(!pred_from.next(pred_node) || curr_node != pred_node + 1)
              {
                sample_this = true; break;
              }
            }
          }
        }

        if(sample_this)
        {
          output.sampledPath(i);
          curr_from.rewind();
          while(curr_from.next(node)) { output.sample(node); }
          output.sampleEnd();
        }
      }
      for(MergedGraphReader& current : reader) { current.close(); }
      lcp_array.close();
      previous.flush(false); stack.flush(false);
      metadata = output.finish();
    }

    TempFile::remove(previous_name); TempFile::remove(stack_name);
    metadata.fast_chars = graph.alpha.fast_chars;
    if(metadata.paths != merged_graph.size())
    {
      throw std::runtime_error("GCSA::GCSA(): final event path count mismatch");
    }
    sortFinalRedundancy(files, parameters);
  }
  catch(...)
  {
    TempFile::remove(previous_name); TempFile::remove(stack_name); throw;
  }
  return metadata;
}

} // namespace

//------------------------------------------------------------------------------

GCSA::GCSA(InputGraph& graph, const ConstructionParameters& parameters) :
  GCSA(graph, parameters, nullptr)
{
}

void
GCSA::buildAndStore(InputGraph& graph, const ConstructionParameters& parameters,
  const std::string& gcsa_filename, const std::string& lcp_filename)
{
  if(!parameters.externalMemory())
  {
    throw std::invalid_argument("GCSA::buildAndStore() requires external-memory construction");
  }
  // TempFile is process-global. This binds the configured directory at the
  // construction boundary, but does not isolate simultaneous constructions
  // that select different directories.
  TempFile::setDirectory(parameters.work_directory);
  SerialExternalConstruction serial_construction;
  if(graph.size() == 0)
  {
    GCSA empty;
    storeEmptyIndexAtomically(empty, gcsa_filename);
  }
  else { GCSA builder(graph, parameters, &gcsa_filename); }
  LCPArray::buildAndStore(graph, parameters, lcp_filename);
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

  if(parameters.externalMemory())
  {
    if(Verbosity::level >= Verbosity::BASIC)
    {
      // This is the aggregate allocator/scheduler goal used to choose spill
      // runs. Library allocations and the embedding application's graph remain
      // outside it, so deployments should report their hard cgroup cap separately.
      std::cerr << "GCSA::GCSA(): External working-set target "
                << formatBytes(parameters.getMemoryLimitBytes())
                << ", disk limit " << formatBytes(parameters.getLimitBytes())
                << std::endl;
    }
  }

  // Extract key and start-node facts. The legacy branch retains its historical
  // vectors exactly. The external preprocessor instead streams immutable
  // record files through bounded sort/merge passes. Only the LCP support is
  // needed during prefix doubling; mapper, last-character, and start-node
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
    external_preprocessor.reset(new ExternalInputPreprocessor(graph, parameters));
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

  // Create the initial PathGraph for this fresh construction.
  PathGraph path_graph(0, graph.k(), 0);
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
  for(size_type step = 1; step <= parameters.getSteps(); step++)
  {
    if(Verbosity::level >= Verbosity::BASIC)
    {
      std::cerr << "GCSA::GCSA(): Step " << step << " (path length " << path_graph.k() << " -> "
                << (2 * path_graph.k()) << ")" << std::endl;
    }
    std::string task = doublingTask(step);
    size_type prune_buffer = pathMergeInputBudget(parameters);
    const double prune_start = readTimer();
    path_graph.prune(lcp, path_graph.remainingLimit(parameters.getLimitBytes()),
      prune_buffer, externalMaxOpenFiles());
    const double prune_stop = readTimer();
    if(Verbosity::level >= Verbosity::BASIC)
    {
      std::cerr << "GCSA::GCSA(): Prune " << task << ": "
                << (prune_stop - prune_start) << " seconds" << std::endl;
    }
    if(parameters.externalMemory())
    {
      externalPathGraphExtend(path_graph,
        path_graph.remainingLimit(parameters.getLimitBytes()), parameters);
    }
    else
    {
      path_graph.extend(path_graph.remainingLimit(parameters.getLimitBytes()),
        parameters.getMemoryLimitBytes());
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
    external_preprocessor->buildMapper(mapper);
  }
  size_type merge_buffer = pathMergeInputBudget(parameters);
  MergedGraph merged_graph(path_graph, mapper, lcp,
    path_graph.remainingLimit(parameters.getLimitBytes()), merge_buffer,
    externalMaxOpenFiles());
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
    FinalEventMetadata event_metadata = produceExternalFinalEvents(merged_graph,
      mapper, last_char, from_nodes, unique_from_nodes, graph,
      parameters, event_files);

    // No final component needs the construction mapper or mutable scan state.
    // Releasing them before component assembly is an important part of the
    // working-set reduction.
    sdsl::util::clear(last_char); sdsl::util::clear(from_nodes);
    sdsl::util::clear(mapper);
    this->header.edges = event_metadata.total_edges;
    if(direct_output != nullptr)
    {
      storeFinalComponents(this->header, graph.alpha, event_files,
        event_metadata, parameters, *direct_output);
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
