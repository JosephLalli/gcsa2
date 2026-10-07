/*
  Copyright (c) 2026 GCSA2 contributors

  Disk-first input preprocessing for GCSA construction.
*/

#ifndef GCSA_EXTERNAL_PREPROCESSING_H
#define GCSA_EXTERNAL_PREPROCESSING_H

#include <gcsa/dbg.h>
#include <gcsa/external_sort.h>
#include <gcsa/path_graph.h>

namespace gcsa
{

/*
  Produces the immutable input facts needed by construction without materializing
  global KMer/key/start-node vectors. Key and start streams use bounded scratch
  artifacts; physical initial PathGraph shards retain their ordinary temporary
  file lifecycle.
*/
class ExternalInputPreprocessor
{
public:
  ExternalInputPreprocessor(const InputGraph& graph,
    const ConstructionParameters& parameters);
  ~ExternalInputPreprocessor();

  // Construction does not need every key support at the same time. Splitting
  // the scans lets GCSA retain only the LCP support during prefix doubling,
  // then create the mapper and last-character vector close to their use.
  void buildMapper(DeBruijnGraph& mapper);
  void buildLCP(LCP& lcp);
  void buildLastCharacters(sdsl::int_vector<0>& last_char);

  // Build the mapped, unique start-node set in two bounded scans.
  void buildStartNodes(sdsl::sd_vector<>& from_nodes);

  // Sort each logical KMer stream deterministically and split it into label-
  // sorted physical PathGraph shards. Every output shard retains its logical ID.
  void buildInitialPathGraph(PathGraph& result);

  size_type startNodeCount() const;

  ExternalInputPreprocessor(const ExternalInputPreprocessor&) = delete;
  ExternalInputPreprocessor& operator=(const ExternalInputPreprocessor&) = delete;

private:
  const InputGraph& graph_;
  const ConstructionParameters& parameters_;
  std::string key_name_, start_name_;
  size_type key_count_, start_count_;
  bool prepared_;

  void prepare();
  size_type sortBudget() const;
  size_type ioBufferBytes() const;
};

} // namespace gcsa

#endif // GCSA_EXTERNAL_PREPROCESSING_H
