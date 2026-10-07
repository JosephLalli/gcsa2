#ifndef GCSA_FINAL_EVENTS_INTERNAL_HPP
#define GCSA_FINAL_EVENTS_INTERNAL_HPP

#include <gcsa/final_events.h>

namespace gcsa
{

// Internal construction seam. The caller supplies the thread count captured
// before external construction temporarily serializes ambient OpenMP work.
void storeFinalComponentsConcurrent(const GCSAHeader& header,
  const Alphabet& source_alphabet, const FinalEventFiles& files,
  const FinalEventMetadata& metadata,
  const ConstructionParameters& parameters, const std::string& filename,
  size_type available_threads, size_type reserved_descriptors = 0);

} // namespace gcsa

#endif // GCSA_FINAL_EVENTS_INTERNAL_HPP
