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
  parameters.setWorkDirectory("/tmp/example.gcsa-work/");
  parameters.setResume(); parameters.setKeepWork();
  parameters.setMemoryLimitBytes(parseBytes("96G"));
  parameters.setLimitBytes(parseBytes("40T"));
  parameters.setIOBufferSize(parseBytes("64M"));
  parameters.setSortRunSize(parseBytes("4G"));
  parameters.setJoinPartitionSize(parseBytes("16G"));
  parameters.setMergeFanIn(32); parameters.setMaxOpenFiles(96);
  require(parameters.externalMemory());
  require(parameters.getWorkDirectory() == "/tmp/example.gcsa-work");
  require(parameters.getResume() && parameters.getKeepWork());
  require(parameters.getMemoryLimitBytes() == 96 * GIGABYTE);
  require(parameters.getLimitBytes() == 40 * KILOBYTE * GIGABYTE);
  require(parameters.getMergeFanIn() == 32 && parameters.getMaxOpenFiles() == 96);
  return 0;
}
