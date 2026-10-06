/*
  Copyright (c) 2026 GCSA2 contributors

  External event streams for final GCSA construction.
*/

#ifndef GCSA_FINAL_EVENTS_H
#define GCSA_FINAL_EVENTS_H

#include <gcsa/external_sort.h>
#include <gcsa/gcsa.h>
#include <gcsa/resources.h>
#include <external_configuration.hpp>

#include <array>
#include <cstdint>
#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

namespace gcsa
{

/*
  A mapped start-node set used by the final merged-graph scan.

  Most paths have only a few start nodes, so the set remains an in-memory
  sorted vector. If a single path has more nodes than the reserved workspace,
  records are streamed to disk, externally sorted and deduplicated, and then
  exposed through the same rewind()/next() interface. The reservation covers
  mutually exclusive collection, external-sort, and readback phases.

  This class is public only so the forced-spill behavior can be unit tested;
  it is a construction helper and is not part of the query interface.
*/
class SpillableNodeSet
{
public:
  SpillableNodeSet(size_type byte_budget, size_type merge_fan_in,
    MemoryBudget& memory);
  ~SpillableNodeSet();

  void clear();
  void push_back(node_type value);
  void finish();

  size_type size() const;
  bool spilled() const;
  void rewind();
  bool next(node_type& value);

  static size_type minimumBudget();

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;

  SpillableNodeSet(const SpillableNodeSet&) = delete;
  SpillableNodeSet& operator=(const SpillableNodeSet&) = delete;
};

/*
  The final merged-graph scan is globally ordered and therefore deliberately
  remains a single atomic task. It writes these immutable streams and commits
  them together. Component construction can then restart without repeating the
  scan, and no raw array proportional to all paths has to coexist with the
  other raw arrays.

  All integer records are explicitly little-endian uint64 values. Occurrence
  records are pairs (path, extra_occurrences); BWT masks are single bytes.
*/
struct FinalEventMetadata
{
  constexpr static size_type MAX_SIGMA = 8;

  size_type paths, sigma, fast_chars, total_edges;
  size_type sampled_paths, sample_ids, sample_bits;
  size_type occurrence_items, occurrence_extra, redundant;
  std::array<size_type, MAX_SIGMA> bwt_counts;

  FinalEventMetadata();
};

struct FinalEventFiles
{
  std::string metadata, bwt_masks;
  std::vector<std::string> edge_destinations;
  std::string sample_positions, sample_ids, sample_ends;
  std::string occurrences, redundant;
  bool delete_files;

  explicit FinalEventFiles(size_type sigma = 0);
  FinalEventFiles(FinalEventFiles&& source);
  FinalEventFiles& operator=(FinalEventFiles&& source);
  ~FinalEventFiles();

  void clear();

  FinalEventFiles(const FinalEventFiles&) = delete;
  FinalEventFiles& operator=(const FinalEventFiles&) = delete;
};

/*
  A byte-budgeted writer set. buffer_bytes is the budget for each individual
  stream, and every buffer reserves from the shared MemoryBudget. The caller
  chooses buffer_bytes after accounting for the number of streams.
*/
class FinalEventWriter
{
public:
  FinalEventWriter(const FinalEventFiles& files, size_type sigma,
    size_type buffer_bytes, MemoryBudget& budget);
  ~FinalEventWriter();

  void path(byte_type predecessor_mask);
  void edge(comp_type comp, size_type source_path);
  void sampledPath(size_type path);
  void sample(node_type node);
  void sampleEnd();
  void occurrence(size_type path, size_type extra);
  void redundancy(size_type path);

  // Closes and fdatasyncs every stream. The returned counts describe the
  // closed payloads; semantic fields such as paths/fast_chars are completed
  // by the merged-graph scan.
  FinalEventMetadata finish();

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;

  FinalEventWriter(const FinalEventWriter&) = delete;
  FinalEventWriter& operator=(const FinalEventWriter&) = delete;
};

// Sort the (potentially nonmonotone) redundancy positions externally.
// stats is optional and diagnostic: run count, merge passes and peak resident
// bytes are what decide whether this sort is worth replacing.
void sortFinalRedundancy(FinalEventFiles& files,
  const ConstructionParameters& parameters,
  ExternalFixedRecordSortStats* stats = nullptr);

// Write/read and validate the small versioned metadata record.
void writeFinalEventMetadata(const FinalEventFiles& files,
  const FinalEventMetadata& metadata);
FinalEventMetadata readFinalEventMetadata(const FinalEventFiles& files);

/*
  Stream one fast BWT component directly in the existing bit_vector_il<512>
  serialization. This avoids materializing both a dense bitvector and its
  interleaved replacement. Public only for exact-format regression tests.
*/
void serializeFastBWTComponent(std::ostream& out,
  const std::string& mask_file, size_type paths, size_type expected_ones,
  comp_type comp, const ConstructionParameters& parameters);

/*
  Stream one sparse BWT component directly in the existing sd_vector<>
  serialization. The mask file is replayed for each Elias--Fano member so no
  builder or completed component is proportional to the path count.
*/
void serializeSparseBWTComponent(std::ostream& out,
  const std::string& mask_file, size_type paths, size_type expected_ones,
  comp_type comp, const ConstructionParameters& parameters);

// Stream the stored-sample int_vector<0> in its existing packed format.
void serializeSampleIds(std::ostream& out, const std::string& sample_file,
  size_type samples, size_type sample_bits,
  const ConstructionParameters& parameters);

// Stream the ordinary sample-boundary bit_vector and select_support_mcl from
// the monotone sample-end event stream. Public for exact-format tests.
void serializeSampleBoundaries(std::ostream& out,
  const std::string& sample_end_file, size_type sample_ids,
  size_type sampled_paths, const ConstructionParameters& parameters);

// Stream the existing SadaSparse occurrence-pointer representation from the
// sorted (path, extra-count) event file. Public for exact-format tests.
void serializeOccurrencePointers(std::ostream& out,
  const std::string& occurrence_file, size_type paths,
  size_type occurrence_items, size_type occurrence_extra,
  const ConstructionParameters& parameters);

// Stream SadaCount's ordinary bit_vector and select_support_mcl payload from
// sorted redundancy events without materializing the dense unary vector.
// Public only for exact-format regression tests.
void serializeRedundantPointers(std::ostream& out,
  const std::string& redundancy_file, size_type paths, size_type redundant,
  const ConstructionParameters& parameters);

/*
  Serialize final components directly from immutable event streams in the
  normal GCSA::load() order. The result has the normal public .gcsa format and
  is published with a synced atomic rename.
*/
void storeFinalComponents(const GCSAHeader& header,
  const Alphabet& source_alphabet, const FinalEventFiles& files,
  const FinalEventMetadata& metadata,
  const ConstructionParameters& parameters, const std::string& filename);

} // namespace gcsa

#endif // GCSA_FINAL_EVENTS_H
