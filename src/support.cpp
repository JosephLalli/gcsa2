#include <gcsa/support.h>
#include <gcsa/internal.h>

namespace gcsa
{

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr size_type ConstructionParameters::DOUBLING_STEPS;
constexpr size_type ConstructionParameters::MAX_STEPS;
constexpr size_type ConstructionParameters::SIZE_LIMIT;
constexpr size_type ConstructionParameters::ABSOLUTE_LIMIT;
constexpr size_type ConstructionParameters::MEMORY_LIMIT;
constexpr size_type ConstructionParameters::ABSOLUTE_MEMORY_LIMIT;
constexpr size_type ConstructionParameters::SAMPLE_PERIOD;
constexpr size_type ConstructionParameters::LCP_BRANCHING;
constexpr size_type ConstructionParameters::IO_BUFFER_SIZE;
constexpr size_type ConstructionParameters::SORT_RUN_SIZE;
constexpr size_type ConstructionParameters::JOIN_PARTITION_SIZE;
constexpr size_type ConstructionParameters::MERGE_FAN_IN;
constexpr size_type ConstructionParameters::MAX_OPEN_FILES;
constexpr size_type ConstructionParameters::MIN_OPEN_FILES;
constexpr size_type ConstructionParameters::CHECKPOINT_RECORDS;
constexpr size_type ConstructionParameters::CHECKPOINT_BYTES;

constexpr Alphabet::size_type Alphabet::MAX_SIGMA;
constexpr Alphabet::size_type Alphabet::SOURCE_COMP;
constexpr Alphabet::size_type Alphabet::SINK_COMP;
constexpr Alphabet::size_type Alphabet::FAST_CHARS;

constexpr size_type Key::GCSA_CHAR_WIDTH;
constexpr key_type Key::CHAR_MASK;
constexpr size_type Key::MAX_LENGTH;
constexpr key_type Key::PRED_SUCC_MASK;

constexpr size_type Node::OFFSET_BITS;
constexpr size_type Node::ID_OFFSET;
constexpr size_type Node::ORIENTATION_MASK;
constexpr size_type Node::OFFSET_MASK;

//------------------------------------------------------------------------------

// Other class variables.

/*
  The default alphabet interprets \0 and $ as endmarkers, ACGT and acgt as ACGT,
  # as a the label of the source node, and the and the remaining characters as N.
*/

const sdsl::int_vector<8> Alphabet::DEFAULT_CHAR2COMP =
{
  0, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 6,  0, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,

  5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,

  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,

  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};

const sdsl::int_vector<8> Alphabet::DEFAULT_COMP2CHAR = { '$', 'A', 'C', 'G', 'T', 'N', '#' };

//------------------------------------------------------------------------------

ConstructionParameters::ConstructionParameters() :
  doubling_steps(DOUBLING_STEPS), size_limit(SIZE_LIMIT * GIGABYTE),
  memory_limit(MEMORY_LIMIT * GIGABYTE), sample_period(SAMPLE_PERIOD),
  lcp_branching(LCP_BRANCHING),
  resume(false), keep_work(false), verify_workspace(false), clean_obsolete(false),
  allow_path_explosion(false), io_buffer_size(IO_BUFFER_SIZE),
  sort_run_size(SORT_RUN_SIZE), join_partition_size(JOIN_PARTITION_SIZE),
  merge_fan_in(MERGE_FAN_IN), max_open_files(MAX_OPEN_FILES),
  process_workers(PROCESS_WORKERS),
  checkpoint_records(CHECKPOINT_RECORDS), checkpoint_bytes(CHECKPOINT_BYTES)
{
}

void
ConstructionParameters::setSteps(size_type steps)
{
  this->doubling_steps = Range::bound(steps, 1, MAX_STEPS);
}

void
ConstructionParameters::setLimit(size_type gigabytes)
{
  this->size_limit = Range::bound(gigabytes, 1, ABSOLUTE_LIMIT) * GIGABYTE;
}

void
ConstructionParameters::setLimitBytes(size_type bytes)
{
  this->size_limit = Range::bound(bytes, 1, ABSOLUTE_LIMIT * GIGABYTE);
}

void
ConstructionParameters::reduceLimit(size_type bytes)
{
  if(bytes > this->size_limit) { this->size_limit = 0; }
  else { this->size_limit -= bytes; }
}

void
ConstructionParameters::setMemoryLimit(size_type gigabytes)
{
  this->memory_limit = Range::bound(gigabytes, 1, ABSOLUTE_MEMORY_LIMIT) * GIGABYTE;
}

void
ConstructionParameters::setMemoryLimitBytes(size_type bytes)
{
  this->memory_limit = Range::bound(bytes, 1, ABSOLUTE_MEMORY_LIMIT * GIGABYTE);
}

void
ConstructionParameters::setSamplePeriod(size_type period)
{
  this->sample_period = std::max((size_type)1, period);
}

void
ConstructionParameters::setLCPBranching(size_type factor)
{
  this->lcp_branching = std::max((size_type)2, factor);
}

void
ConstructionParameters::setWorkDirectory(const std::string& directory)
{
  this->work_directory = directory;
  while(this->work_directory.length() > 1 && this->work_directory.back() == '/')
  {
    this->work_directory.pop_back();
  }
}

void ConstructionParameters::setResume(bool value) { this->resume = value; }
void ConstructionParameters::setKeepWork(bool value) { this->keep_work = value; }
void ConstructionParameters::setVerifyWorkspace(bool value) { this->verify_workspace = value; }
void ConstructionParameters::setCleanObsolete(bool value) { this->clean_obsolete = value; }
void ConstructionParameters::setAllowPathExplosion(bool value) { this->allow_path_explosion = value; }

void
ConstructionParameters::setIOBufferSize(size_type bytes)
{
  this->io_buffer_size = std::max((size_type)KILOBYTE, bytes);
}

void
ConstructionParameters::setSortRunSize(size_type bytes)
{
  this->sort_run_size = std::max((size_type)KILOBYTE, bytes);
}

void
ConstructionParameters::setJoinPartitionSize(size_type bytes)
{
  this->join_partition_size = std::max((size_type)KILOBYTE, bytes);
}

size_type
ConstructionParameters::getSortRunSize() const
{
  return this->getSortRunSize(this->memory_limit);
}

size_type
ConstructionParameters::getSortRunSize(size_type available_bytes) const
{
  if(!(this->sortRunSizeIsAutomatic()))
  {
    return std::min(this->sort_run_size, available_bytes);
  }
  // Subtract rather than multiply to avoid overflow near the absolute limit.
  size_type join_share = available_bytes / 4;
  return available_bytes - join_share;
}

size_type
ConstructionParameters::getJoinPartitionSize() const
{
  return this->getJoinPartitionSize(this->memory_limit);
}

size_type
ConstructionParameters::getJoinPartitionSize(size_type available_bytes) const
{
  if(!(this->joinPartitionSizeIsAutomatic()))
  {
    return std::min(this->join_partition_size, available_bytes);
  }
  return available_bytes / 4;
}

void
ConstructionParameters::setMergeFanIn(size_type fan_in)
{
  this->merge_fan_in = std::max((size_type)2, fan_in);
}

void
ConstructionParameters::setMaxOpenFiles(size_type files)
{
  // The final event scan has more fixed streams than PathGraph pruning and
  // merging. Keep enough headroom for its bounded two-way sorter as well.
  this->max_open_files = std::max(MIN_OPEN_FILES, files);
}

void
ConstructionParameters::setProcessWorkers(size_type workers)
{
  this->process_workers = std::max((size_type)1, workers);
}

void
ConstructionParameters::setWorkerExecutable(const std::string& executable)
{
  this->worker_executable = executable;
}

void
ConstructionParameters::setCheckpointRecords(size_type records)
{
  this->checkpoint_records = std::max((size_type)1, records);
}

void
ConstructionParameters::setCheckpointBytes(size_type bytes)
{
  this->checkpoint_bytes = std::max((size_type)KILOBYTE, bytes);
}

void
ConstructionParameters::setStopAfter(const std::string& phase)
{
  this->stop_after = phase;
}

//------------------------------------------------------------------------------

Alphabet::Alphabet() :
  char2comp(DEFAULT_CHAR2COMP), comp2char(DEFAULT_COMP2CHAR),
  C(DEFAULT_COMP2CHAR.size() + 1, 0),
  sigma(DEFAULT_COMP2CHAR.size()), fast_chars(FAST_CHARS)
{
}

Alphabet::Alphabet(const Alphabet& a)
{
  this->copy(a);
}

Alphabet::Alphabet(Alphabet&& a) noexcept
{
  *this = std::move(a);
}

Alphabet::Alphabet(const sdsl::int_vector<64>& counts,
  const sdsl::int_vector<8>& _char2comp, const sdsl::int_vector<8>& _comp2char) :
  char2comp(_char2comp), comp2char(_comp2char),
  C(_comp2char.size() + 1, 0),
  sigma(_comp2char.size()), fast_chars(FAST_CHARS)
{
  for(size_type i = 0; i < counts.size(); i++) { this->C[i + 1] = this->C[i] + counts[i]; }
}

Alphabet::~Alphabet()
{
}

void
Alphabet::copy(const Alphabet& a)
{
  this->char2comp = a.char2comp;
  this->comp2char = a.comp2char;
  this->C = a.C;
  this->sigma = a.sigma;
  this->fast_chars = a.fast_chars;
}

void
Alphabet::swap(Alphabet& a) noexcept
{
  if(this != &a)
  {
    this->char2comp.swap(a.char2comp);
    this->comp2char.swap(a.comp2char);
    this->C.swap(a.C);
    std::swap(this->sigma, a.sigma);
    std::swap(this->fast_chars, a.fast_chars);
  }
}

Alphabet&
Alphabet::operator=(const Alphabet& a)
{
  if(this != &a) { this->copy(a); }
  return *this;
}

Alphabet&
Alphabet::operator=(Alphabet&& a) noexcept
{
  if(this != &a)
  {
    this->char2comp = std::move(a.char2comp);
    this->comp2char = std::move(a.comp2char);
    this->C = std::move(a.C);
    this->sigma = a.sigma;
    this->fast_chars = a.fast_chars;
  }
  return *this;
}

Alphabet::size_type
Alphabet::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;
  written_bytes += this->char2comp.serialize(out, child, "char2comp");
  written_bytes += this->comp2char.serialize(out, child, "comp2char");
  written_bytes += this->C.serialize(out, child, "C");
  written_bytes += sdsl::write_member(this->sigma, out, child, "sigma");
  written_bytes += sdsl::write_member(this->fast_chars, out, child, "fast_chars");
  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
Alphabet::load(std::istream& in)
{
  this->char2comp.load(in);
  this->comp2char.load(in);
  this->C.load(in);
  sdsl::read_member(this->sigma, in);
  sdsl::read_member(this->fast_chars, in);
}

//------------------------------------------------------------------------------

NodeMapping::NodeMapping() :
  first_node(~(size_type)0), next_node(~(size_type)0)
{
}

NodeMapping::NodeMapping(const NodeMapping& source)
{
  this->copy(source);
}

NodeMapping::NodeMapping(NodeMapping&& source) noexcept
{
  *this = std::move(source);
}

NodeMapping::~NodeMapping()
{
}

NodeMapping::NodeMapping(size_type first_node_id) :
  first_node(first_node_id), next_node(first_node_id)
{
}

void
NodeMapping::copy(const NodeMapping& source)
{
  this->first_node = source.first_node;
  this->next_node = source.next_node;
  this->mapping = source.mapping;
}

void
NodeMapping::swap(NodeMapping& source) noexcept
{
  if(this != &source)
  {
    std::swap(this->first_node, source.first_node);
    std::swap(this->next_node, source.next_node);
    this->mapping.swap(source.mapping);
  }
}

NodeMapping&
NodeMapping::operator=(const NodeMapping& source)
{
  if(this != &source) { this->copy(source); }
  return *this;
}

NodeMapping&
NodeMapping::operator=(NodeMapping&& source) noexcept
{
  if(this != &source)
  {
    this->first_node = source.first_node;
    this->next_node = source.next_node;
    this->mapping = std::move(source.mapping);
  }
  return *this;
}

NodeMapping::size_type
NodeMapping::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;
  written_bytes += sdsl::write_member(this->first_node, out, child, "first_node");
  written_bytes += sdsl::write_member(this->next_node, out, child, "next_node");

  // Serialize the data.
  size_type data_bytes = this->size() * sizeof(size_type);
  sdsl::structure_tree_node* data_node = sdsl::structure_tree::add_child(child, "mapping", "std::vector<gcsa::size_type>");
  if(this->size() > 0) { DiskIO::write(out, this->mapping.data(), this->size(), false); }
  sdsl::structure_tree::add_size(data_node, data_bytes);
  written_bytes += data_bytes;

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
NodeMapping::load(std::istream& in)
{
  sdsl::read_member(this->first_node, in);
  sdsl::read_member(this->next_node, in);
  this->mapping.resize(this->next_node - this->first_node);
  if(this->size() > 0) { DiskIO::read(in, this->mapping.data(), this->size(), false); }
}

size_type
NodeMapping::insert(size_type node_id)
{
  this->mapping.push_back(node_id);
  return this->next_node++;
}

//------------------------------------------------------------------------------

SadaCount::SadaCount()
{
}

SadaCount::SadaCount(const SadaCount& source)
{
  this->copy(source);
}

SadaCount::SadaCount(SadaCount&& source) noexcept
{
  *this = std::move(source);
}

SadaCount::~SadaCount()
{
}

void
SadaCount::swap(SadaCount& another) noexcept
{
  if(this != &another)
  {
    this->data.swap(another.data);
    sdsl::util::swap_support(this->select, another.select, &(this->data), &(another.data));
  }
}

SadaCount&
SadaCount::operator=(const SadaCount& source)
{
  if(this != &source) { this->copy(source); }
  return *this;
}

SadaCount&
SadaCount::operator=(SadaCount&& source) noexcept
{
  if(this != &source)
  {
    this->data = std::move(source.data);
    this->select = std::move(source.select);
    this->setVectors();
  }
  return *this;
}

SadaCount::size_type
SadaCount::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += this->data.serialize(out, child, "data");
  written_bytes += this->select.serialize(out, child, "select");

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
SadaCount::load(std::istream& in)
{
  this->data.load(in);
  this->select.load(in, &(this->data));
}

void
SadaCount::copy(const SadaCount& source)
{
  this->data = source.data;
  this->select = source.select;
  this->setVectors();
}

void
SadaCount::setVectors()
{
  this->select.set_vector(&(this->data));
}

//------------------------------------------------------------------------------

SadaSparse::SadaSparse()
{
}

SadaSparse::SadaSparse(const SadaSparse& source)
{
  this->copy(source);
}

SadaSparse::SadaSparse(SadaSparse&& source) noexcept
{
  *this = std::move(source);
}

SadaSparse::~SadaSparse()
{
}

void
SadaSparse::swap(SadaSparse& another) noexcept
{
  if(this != &another)
  {
    this->filter.swap(another.filter);
    sdsl::util::swap_support(this->filter_rank, another.filter_rank,
      &(this->filter), &(another.filter));

    this->values.swap(another.values);
    sdsl::util::swap_support(this->value_select, another.value_select,
      &(this->values), &(another.values));
  }
}

SadaSparse&
SadaSparse::operator=(const SadaSparse& source)
{
  if(this != &source) { this->copy(source); }
  return *this;
}

SadaSparse&
SadaSparse::operator=(SadaSparse&& source) noexcept
{
  if(this != &source)
  {
    this->filter = std::move(source.filter);
    this->filter_rank = std::move(source.filter_rank);

    this->values = std::move(source.values);
    this->value_select = std::move(source.value_select);

    this->setVectors();
  }
  return *this;
}

SadaSparse::size_type
SadaSparse::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += this->filter.serialize(out, child, "filter");
  written_bytes += this->filter_rank.serialize(out, child, "filter_rank");

  written_bytes += this->values.serialize(out, child, "values");
  written_bytes += this->value_select.serialize(out, child, "value_select");

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
SadaSparse::load(std::istream& in)
{
  this->filter.load(in);
  this->filter_rank.load(in, &(this->filter));

  this->values.load(in);
  this->value_select.load(in, &(this->values));
}

void
SadaSparse::copy(const SadaSparse& source)
{
  this->filter = source.filter;
  this->filter_rank = source.filter_rank;

  this->values = source.values;
  this->value_select = source.value_select;

  this->setVectors();
}

void
SadaSparse::setVectors()
{
  this->filter_rank.set_vector(&(this->filter));
  this->value_select.set_vector(&(this->values));
}

//------------------------------------------------------------------------------

std::string
Key::decode(key_type key, size_type kmer_length, const Alphabet& alpha)
{
  key = label(key);
  kmer_length = std::min(kmer_length, Key::MAX_LENGTH);

  std::string res(kmer_length, '\0');
  for(size_type i = 1; i <= kmer_length; i++)
  {
    res[kmer_length - i] = alpha.comp2char[key & CHAR_MASK];
    key >>= GCSA_CHAR_WIDTH;
  }

  return res;
}

void
Key::lastChars(const std::vector<key_type>& keys, sdsl::int_vector<0>& last_char)
{
  sdsl::util::clear(last_char);
  last_char = sdsl::int_vector<0>(keys.size(), 0, GCSA_CHAR_WIDTH);
  for(size_type i = 0; i < keys.size(); i++) { last_char[i] = Key::last(keys[i]); }
}

//------------------------------------------------------------------------------

node_type
Node::encode(const std::string& token)
{
  size_t separator = 0;
  size_type node = std::stoul(token, &separator);
  if(separator + 1 >= token.length())
  {
    std::cerr << "Node::encode(): Invalid position token " << token << std::endl;
    return 0;
  }

  bool reverse_complement = false;
  if(token[separator + 1] == '-')
  {
    reverse_complement = true;
    separator++;
  }

  std::string temp = token.substr(separator + 1);
  size_type offset = std::stoul(temp);
  if(offset > OFFSET_MASK)
  {
    std::cerr << "Node::encode(): Offset " << offset << " too large" << std::endl;
    return 0;
  }

  return encode(node, offset, reverse_complement);
}

std::string
Node::decode(node_type node)
{
  std::ostringstream ss;
  ss << id(node) << ':';
  if(rc(node)) { ss << '-'; }
  ss << offset(node);
  return ss.str();
}

void
Node::map(std::vector<node_type>& nodes, const NodeMapping& mapping)
{
  if(mapping.empty()) { return; }
  for(node_type& node : nodes)
  {
    node = encode(mapping(id(node)), offset(node), rc(node));
  }
}

//------------------------------------------------------------------------------

KMer::KMer()
{
}

KMer::KMer(const std::vector<std::string>& tokens, const Alphabet& alpha, size_type successor)
{
  byte_type predecessors = chars(tokens[2], alpha);
  byte_type successors = chars(tokens[3], alpha);
  this->key = Key::encode(alpha, tokens[0], predecessors, successors);
  this->from = Node::encode(tokens[1]);
  this->to = Node::encode(tokens[successor]);
}

byte_type
KMer::chars(const std::string& token, const Alphabet& alpha)
{
  byte_type val = 0;
  for(size_type i = 0; i < token.length(); i += 2) { val |= 1 << alpha.char2comp[token[i]]; }
  return val;
}

std::ostream&
operator<< (std::ostream& out, const KMer& kmer)
{
  out << "(key " << Key::label(kmer.key)
      << ", in " << (size_type)(Key::predecessors(kmer.key))
      << ", out " << (size_type)(Key::successors(kmer.key))
      << ", from " << Node::decode(kmer.from)
      << ", to " << Node::decode(kmer.to) << ")";
  return out;
}

//------------------------------------------------------------------------------

} // namespace gcsa
