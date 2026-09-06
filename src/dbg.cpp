#include <gcsa/dbg.h>
#include <gcsa/internal.h>

#include <fstream>
#include <stdexcept>
#include <vector>

namespace gcsa
{

namespace
{

template<class Callback>
void
scanKeys(const std::string& name, size_type expected_records,
  size_type buffer_bytes, const Callback& callback)
{
  if(buffer_bytes < sizeof(key_type)) { buffer_bytes = sizeof(key_type); }
  // Clamp to the records that exist. expected_records is already in hand and is
  // validated against below, so a 64 MiB buffer_bytes was value-initialising
  // 64 MiB twice per DeBruijnGraph for a stream that may hold far less.
  size_type records_per_block = std::max(static_cast<size_type>(1),
    std::min(buffer_bytes / sizeof(key_type),
      std::max(static_cast<size_type>(1), expected_records)));
  std::ifstream input;
  input.rdbuf()->pubsetbuf(nullptr, 0);
  input.open(name.c_str(), std::ios_base::binary);
  if(!input) { throw std::runtime_error("DeBruijnGraph: cannot open key stream " + name); }
  std::vector<key_type> buffer(records_per_block);
  size_type seen = 0;
  while(seen < expected_records)
  {
    size_type count = std::min(records_per_block, expected_records - seen);
    if(!DiskIO::read(input, buffer.data(), count))
    {
      throw std::runtime_error("DeBruijnGraph: truncated key stream " + name);
    }
    for(size_type i = 0; i < count; i++) { callback(buffer[i], seen + i); }
    seen += count;
  }
  key_type extra;
  input.read(reinterpret_cast<char*>(&extra), sizeof(extra));
  if(input.gcount() != 0)
  {
    throw std::runtime_error("DeBruijnGraph: key stream has trailing records " + name);
  }
}

} // namespace

//------------------------------------------------------------------------------

// Other class variables.

const std::string DeBruijnGraph::EXTENSION = ".dbg";

//------------------------------------------------------------------------------

DeBruijnGraph::DeBruijnGraph()
{
  this->node_count = 0;
  this->graph_order = 0;
}

DeBruijnGraph::DeBruijnGraph(const DeBruijnGraph& g)
{
  this->copy(g);
}

DeBruijnGraph::DeBruijnGraph(DeBruijnGraph&& g) noexcept
{
  *this = std::move(g);
}

DeBruijnGraph::~DeBruijnGraph()
{
}

void
DeBruijnGraph::copy(const DeBruijnGraph& g)
{
  this->node_count = g.node_count;
  this->graph_order = g.graph_order;

  this->alpha = g.alpha;

  this->bwt = g.bwt;
  this->bwt_rank = g.bwt_rank;

  this->nodes = g.nodes;
  this->node_rank = g.node_rank;

  this->setVectors();
}

void
DeBruijnGraph::swap(DeBruijnGraph& g) noexcept
{
  if(this != &g)
  {
    std::swap(this->node_count, g.node_count);
    std::swap(this->graph_order, g.graph_order);

    this->alpha.swap(g.alpha);

    this->bwt.swap(g.bwt);
    sdsl::util::swap_support(this->bwt_rank, g.bwt_rank, &(this->bwt), &(g.bwt));

    this->nodes.swap(g.nodes);
    sdsl::util::swap_support(this->node_rank, g.node_rank, &(this->nodes), &(g.nodes));
  }
}

DeBruijnGraph&
DeBruijnGraph::operator=(const DeBruijnGraph& g)
{
  if(this != &g) { this->copy(g); }
  return *this;
}

DeBruijnGraph&
DeBruijnGraph::operator=(DeBruijnGraph&& g) noexcept
{
  if(this != &g)
  {
    this->node_count = std::move(g.node_count);
    this->graph_order = std::move(g.graph_order);

    this->alpha = std::move(g.alpha);

    this->bwt = std::move(g.bwt);
    this->bwt_rank = std::move(g.bwt_rank);

    this->nodes = std::move(g.nodes);
    this->node_rank = std::move(g.node_rank);

    this->setVectors();
  }
  return *this;
}

void
DeBruijnGraph::setVectors()
{
  this->bwt_rank.set_vector(&(this->bwt));
  this->node_rank.set_vector(&(this->nodes));
}

DeBruijnGraph::size_type
DeBruijnGraph::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += sdsl::write_member(this->node_count, out, child, "node_count");
  written_bytes += sdsl::write_member(this->graph_order, out, child, "graph_order");

  written_bytes += this->alpha.serialize(out, child, "alpha");

  written_bytes += this->bwt.serialize(out, child, "bwt");
  written_bytes += this->bwt_rank.serialize(out, child, "bwt_rank");

  written_bytes += this->nodes.serialize(out, child, "nodes");
  written_bytes += this->node_rank.serialize(out, child, "node_rank");

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
DeBruijnGraph::load(std::istream& in)
{
  sdsl::read_member(this->node_count, in);
  sdsl::read_member(this->graph_order, in);

  this->alpha.load(in);

  this->bwt.load(in);
  this->bwt_rank.load(in, &(this->bwt));

  this->nodes.load(in);
  this->node_rank.load(in, &(this->nodes));
}

//------------------------------------------------------------------------------

DeBruijnGraph::DeBruijnGraph(const std::vector<key_type>& keys, size_type kmer_length, const Alphabet& alphabet)
{
  this->node_count = keys.size();
  this->graph_order = kmer_length;

  size_type total_edges = 0;
  for(size_type i = 0; i < keys.size(); i++) { total_edges += sdsl::bits::lt_cnt[Key::predecessors(keys[i])]; }

  sdsl::int_vector<64> counts(alphabet.sigma, 0);
  bit_vector bwt_buffer(alphabet.sigma * total_edges, 0);
  bit_vector node_buffer(total_edges, 0);
  for(size_type i = 0, edge_pos = 0; i < keys.size(); i++)
  {
    size_type pred = Key::predecessors(keys[i]);
    for(size_type j = 0; j < alphabet.sigma; j++)
    {
      if(pred & (((size_type)1) << j))
      {
        bwt_buffer[j * this->size() + i] = 1;
        counts[j]++;
      }
    }
    edge_pos += sdsl::bits::lt_cnt[Key::successors(keys[i])];
    node_buffer[edge_pos - 1] = 1;
  }
  this->alpha = Alphabet(counts, alphabet.char2comp, alphabet.comp2char);
  this->bwt = bwt_buffer; sdsl::util::clear(bwt_buffer);
  this->nodes = node_buffer; sdsl::util::clear(node_buffer);

  this->initSupport();
}

DeBruijnGraph::DeBruijnGraph(const std::string& key_name, size_type key_count,
  size_type kmer_length, const Alphabet& alphabet, size_type buffer_bytes)
{
  this->node_count = key_count;
  this->graph_order = kmer_length;

  size_type total_edges = 0;
  sdsl::int_vector<64> counts(alphabet.sigma, 0);
  scanKeys(key_name, key_count, buffer_bytes, [&total_edges, &counts](key_type key, size_type)
  {
    total_edges += sdsl::bits::lt_cnt[Key::predecessors(key)];
    for(size_type comp = 0; comp < counts.size(); comp++)
    {
      if(Key::predecessors(key) & (static_cast<size_type>(1) << comp)) { counts[comp]++; }
    }
  });

  bit_vector bwt_buffer(alphabet.sigma * this->size(), 0);
  bit_vector node_buffer(total_edges, 0);
  size_type edge_pos = 0;
  scanKeys(key_name, key_count, buffer_bytes,
    [&bwt_buffer, &node_buffer, &edge_pos, this, &alphabet](key_type key, size_type rank)
  {
    size_type predecessors = Key::predecessors(key);
    for(size_type comp = 0; comp < alphabet.sigma; comp++)
    {
      if(predecessors & (static_cast<size_type>(1) << comp))
      {
        bwt_buffer[comp * this->size() + rank] = 1;
      }
    }
    size_type outdegree = sdsl::bits::lt_cnt[Key::successors(key)];
    if(outdegree == 0)
    {
      throw std::runtime_error("DeBruijnGraph: key stream contains a node without successors");
    }
    edge_pos += outdegree; node_buffer[edge_pos - 1] = 1;
  });
  if(edge_pos != total_edges)
  {
    throw std::runtime_error("DeBruijnGraph: inconsistent key stream edge count");
  }
  this->alpha = Alphabet(counts, alphabet.char2comp, alphabet.comp2char);
  this->bwt = bwt_buffer; sdsl::util::clear(bwt_buffer);
  this->nodes = node_buffer; sdsl::util::clear(node_buffer);
  this->initSupport();
}

//------------------------------------------------------------------------------

void
DeBruijnGraph::initSupport()
{
  sdsl::util::init_support(this->bwt_rank, &(this->bwt));
  sdsl::util::init_support(this->node_rank, &(this->nodes));
}

//------------------------------------------------------------------------------

} // namespace gcsa
