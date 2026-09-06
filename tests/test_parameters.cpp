#include <gcsa/support.h>

#include <cstdlib>
#include <stdexcept>

using namespace gcsa;

static void require(bool value)
{
  if(!value) { std::abort(); }
}

int main()
{
  require(parseBytes("64M") == 64 * MEGABYTE);
  require(parseBytes("1.5GiB") == GIGABYTE + GIGABYTE / 2);
  require(parseBytes("40T") == 40 * KILOBYTE * GIGABYTE);
  bool invalid = false;
  try { parseBytes("12XB"); }
  catch(const std::invalid_argument&) { invalid = true; }
  require(invalid);

  ConstructionParameters parameters;
  parameters.setMemoryLimitBytes(parseBytes("96G"));
  require(parameters.sortRunSizeIsAutomatic());
  require(parameters.joinPartitionSizeIsAutomatic());
  require(parameters.getSortRunSize() == 72 * GIGABYTE);
  require(parameters.getJoinPartitionSize() == 24 * GIGABYTE);
  require(parameters.getSortRunSize(8 * GIGABYTE) == 6 * GIGABYTE);
  require(parameters.getJoinPartitionSize(8 * GIGABYTE) == 2 * GIGABYTE);
  parameters.setMaxOpenFiles(1);
  require(parameters.getMaxOpenFiles() == ConstructionParameters::MIN_OPEN_FILES);
  parameters.setWorkDirectory("/tmp/example.gcsa-work/");
  parameters.setResume(); parameters.setKeepWork();
  parameters.setLimitBytes(parseBytes("40T"));
  parameters.setIOBufferSize(parseBytes("64M"));
  parameters.setSortRunSize(parseBytes("4G"));
  parameters.setJoinPartitionSize(parseBytes("16G"));
  require(!(parameters.sortRunSizeIsAutomatic()));
  require(!(parameters.joinPartitionSizeIsAutomatic()));
  require(parameters.getSortRunSize() == 4 * GIGABYTE);
  require(parameters.getJoinPartitionSize() == 16 * GIGABYTE);
  parameters.setMergeFanIn(32); parameters.setMaxOpenFiles(96);
  parameters.setProcessWorkers(4);
  parameters.setTempCompression("zstd");
  parameters.setCompressionBlockSize(parseBytes("8M"));
  parameters.setCompressionWorkers(3);
  parameters.setCompressionLevel(2);
  parameters.setWorkerExecutable("/tmp/build_gcsa");
  require(parameters.externalMemory());
  require(parameters.getWorkDirectory() == "/tmp/example.gcsa-work");
  require(parameters.getResume() && parameters.getKeepWork());
  require(parameters.getMemoryLimitBytes() == 96 * GIGABYTE);
  require(parameters.getLimitBytes() == 40 * KILOBYTE * GIGABYTE);
  require(parameters.getMergeFanIn() == 32 && parameters.getMaxOpenFiles() == 96);
  require(parameters.getProcessWorkers() == 4);
  require(parameters.getTempCompression() == TempCompression::ZSTD);
  require(std::string(tempCompressionName(parameters.getTempCompression())) == "zstd");
  require(parameters.getCompressionBlockSize() == 8 * MEGABYTE);
  require(parameters.getCompressionWorkers() == 3);
  require(parameters.getCompressionLevel() == 2);
  require(parameters.getWorkerExecutable() == "/tmp/build_gcsa");

  invalid = false;
  try { parameters.setTempCompression("gzip"); }
  catch(const std::invalid_argument&) { invalid = true; }
  require(invalid);
  // The disk ceiling depends on the route, and the two setters may be called
  // in either order, so re-clamping must make them commute.
  {
    ConstructionParameters external;
    external.setWorkDirectory("/tmp/gcsa-limit-probe");
    external.setLimit(65536);                      // 64 TiB, above the release cap
    require(external.getLimitBytes() == 65536 * GIGABYTE);

    ConstructionParameters reordered;
    reordered.setLimit(65536);                     // clamped to the release cap
    require(reordered.getLimitBytes() == ConstructionParameters::ABSOLUTE_LIMIT * GIGABYTE);
    reordered.setWorkDirectory("/tmp/gcsa-limit-probe");
    reordered.setLimit(65536);                     // now the external cap applies
    require(reordered.getLimitBytes() == 65536 * GIGABYTE);

    ConstructionParameters legacy;
    legacy.setWorkDirectory("/tmp/gcsa-limit-probe");
    legacy.setLimit(65536);
    legacy.setWorkDirectory("");                   // leaving the external route
    require(legacy.getLimitBytes() == ConstructionParameters::ABSOLUTE_LIMIT * GIGABYTE);
  }

  return 0;
}
