#ifndef GCSA_EXTERNAL_CONFIGURATION_H
#define GCSA_EXTERNAL_CONFIGURATION_H

#include <gcsa/support.h>

namespace gcsa
{

inline size_type externalIOBufferSize(const ConstructionParameters&)
{
  return 64 * MEGABYTE;
}

inline size_type externalSortRunSize(const ConstructionParameters&, size_type available_bytes)
{
  return available_bytes - available_bytes / 4;
}

inline size_type externalSortRunSize(const ConstructionParameters& parameters)
{
  return externalSortRunSize(parameters, parameters.getMemoryLimitBytes());
}

inline size_type externalMergeFanIn()
{
  return 64;
}

inline size_type externalMaxOpenFiles()
{
  return 128;
}

inline size_type externalPreprocessingShardRecords()
{
  return 16 * MILLION;
}

} // namespace gcsa

#endif // GCSA_EXTERNAL_CONFIGURATION_H
