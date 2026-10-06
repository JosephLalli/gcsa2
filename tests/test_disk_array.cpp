#include <gcsa/disk_array.h>

#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <unistd.h>

using namespace gcsa;

namespace
{

void require(bool value) { if(!value) { std::abort(); } }

} // namespace

int
main()
{
  char root[] = "/tmp/gcsa-disk-array-XXXXXX";
  require(mkdtemp(root) != nullptr);
  std::string filename = std::string(root) + "/array.bin";

  // Two 64-byte cache blocks force repeated direct-map replacement while the
  // logical array remains much larger than the configured resident payload.
  MemoryBudget budget(256);
  {
    DiskBackedArray64 array(filename, 1024, 160, budget, true, 64);
    for(size_type i = 0; i < array.size(); i += 17)
    {
      require(array.get(i) == 0);
      array.set(i, 3 * i + 7);
    }
    array.set(1, 99); array.set(1023, 12345);
    array.flush();
    require(array.stats().replacements > 0);
    require(budget.stats().maximum <= 160);
  }
  require(budget.stats().current == 0);

  {
    DiskBackedArray64 array(filename, 1024, 80, budget, false, 32);
    for(size_type i = 0; i < array.size(); i++)
    {
      std::uint64_t expected = (i % 17 == 0 ? 3 * i + 7 : 0);
      if(i == 1) { expected = 99; }
      if(i == 1023) { expected = 12345; }
      require(array.get(i) == expected);
    }
  }

  bool rejected = false;
  try
  {
    DiskBackedArray64 array(filename, 1024,
      DiskBackedArray64::minimumCacheBytes() - 1, budget, false, 8);
  }
  catch(const std::runtime_error&) { rejected = true; }
  require(rejected);

  require(::unlink(filename.c_str()) == 0);
  require(::rmdir(root) == 0);
  return 0;
}
