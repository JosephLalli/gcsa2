#ifndef GCSA_DBG_H
#define GCSA_DBG_H

#include <gcsa/support.h>

namespace gcsa
{

/*
  dbg.h: Specialization of GCSA for de Bruijn graphs.
*/

//------------------------------------------------------------------------------

/*
  This is a specialization of GCSA for de Bruijn graphs. Because the graph is
  reverse deterministic, we can use indicator bitvectors for encoding the BWT.
  This simplifies LF() to two rank() operations.
*/

class DeBruijnGraph
{
public:
  typedef gcsa::size_type       size_type;
  typedef sdsl::bit_vector      bit_vector;
  typedef sdsl::bit_vector_il<> fast_bit_vector;

//------------------------------------------------------------------------------

  DeBruijnGraph();
  DeBruijnGraph(const DeBruijnGraph& g);
  DeBruijnGraph(DeBruijnGraph&& g) noexcept;
  ~DeBruijnGraph();

  void swap(DeBruijnGraph& g) noexcept;
  DeBruijnGraph& operator=(const DeBruijnGraph& g);
  DeBruijnGraph& operator=(DeBruijnGraph&& g) noexcept;

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  const static std::string EXTENSION; // .dbg

//------------------------------------------------------------------------------

  /*
    The input vector must be sorted and contain only unique kmers of length 16 or less.
    Each kmer must have at least one predecessor and one successor.
  */
  DeBruijnGraph(const std::vector<key_type>& keys, size_type kmer_length, const Alphabet& alphabet);

//------------------------------------------------------------------------------

  inline size_type size() const { return this->node_count; }
  inline size_type order() const { return this->graph_order; }
  inline size_type edgeCount() const { return this->nodes.size(); }

  inline size_type LF(size_type node, comp_type comp) const
  {
    return this->node_rank(this->bwt_rank(comp * this->size() + node));
  }

  inline range_type LF(range_type range, comp_type comp) const
  {
    return range_type(this->LF(range.first, comp), this->LF(range.second + 1, comp) - 1);
  }

  inline range_type charRange(comp_type comp) const
  {
    return range_type(this->node_rank(this->alpha.C[comp]), this->node_rank(this->alpha.C[comp + 1]) - 1);
  }

  template<class Iterator>
  range_type find(Iterator begin, Iterator end) const
  {
    if(begin == end) { return range_type(0, this->size() - 1); }
    if((size_type)(end - begin) > this->order())
    {
      std::cerr << "DeBruijnGraph::find(): Query length exceeds " << this->order() << std::endl;
      return range_type(1, 0);
    }

    --end;
    range_type range = this->charRange(this->alpha.char2comp[*end]);
    while(!Range::empty(range) && end != begin)
    {
      --end;
      range = this->LF(range, this->alpha.char2comp[*end]);
    }

    return range;
  }

  template<class Container>
  range_type find(const Container& pattern) const
  {
    return this->find(pattern.begin(), pattern.end());
  }

  template<class Element>
  range_type find(const Element* pattern, size_type length) const
  {
    return this->find(pattern, pattern + length);
  }

//------------------------------------------------------------------------------

  size_type                    node_count;
  size_type                    graph_order;

  Alphabet                     alpha;

  // If node i has predecessor comp, bit comp * size + i is set.
  fast_bit_vector              bwt;
  fast_bit_vector::rank_1_type bwt_rank;

  // The last outgoing edge from each path is marked with an 1-bit.
  fast_bit_vector              nodes;
  fast_bit_vector::rank_1_type node_rank;

//------------------------------------------------------------------------------

private:
  void copy(const DeBruijnGraph& g);
  void setVectors();
  void initSupport();
};  // class DeBruijnGraph

//------------------------------------------------------------------------------

} // namespace gcsa

#endif // GCSA_DBG_H
