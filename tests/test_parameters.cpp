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
  const std::string original_temp_dir = TempFile::temp_dir;
  require(parseBytes("64M") == 64 * MEGABYTE);
  require(parseBytes("1.5GiB") == GIGABYTE + GIGABYTE / 2);
  require(parseBytes("40T") == 40 * KILOBYTE * GIGABYTE);
  bool invalid = false;
  try { parseBytes("12XB"); }
  catch(const std::invalid_argument&) { invalid = true; }
  require(invalid);

  ConstructionParameters parameters;
  parameters.setMemoryLimitBytes(parseBytes("96G"));
  parameters.setWorkDirectory("/tmp/example.gcsa-work/");
  parameters.setLimitBytes(parseBytes("40T"));
  require(parameters.externalMemory());
  require(TempFile::temp_dir == original_temp_dir);
  require(parameters.getMemoryLimitBytes() == 96 * GIGABYTE);
  require(parameters.getLimitBytes() == 40 * KILOBYTE * GIGABYTE);
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
    require(TempFile::temp_dir == original_temp_dir);
  }

  return 0;
}
