#include <gcsa/internal.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

using namespace gcsa;

namespace
{

void require(bool value)
{
  if(!value) { std::abort(); }
}

bool exists(const std::string& filename)
{
  struct stat status;
  return (::stat(filename.c_str(), &status) == 0);
}

} // namespace

int
main()
{
  char root[] = "/tmp/gcsa-adaptive-record-store-XXXXXX";
  require(mkdtemp(root) != nullptr);
  TempFile::setDirectory(root);

  // Pre-admitted resident data does not create a raw temporary, and growing
  // from two to four records charges the old and new capacities together.
  {
    MemoryBudget budget(8 * sizeof(size_type));
    std::shared_ptr<AdaptiveRecordStore<size_type>> store(
      new AdaptiveRecordStore<size_type>(budget, "adaptive-test-resident"));
    require(store->tryReserveRecords(2));
    store->append(7); store->append(11);
    require(store->resident() && store->count() == 2);
    require(store->residentData()[0] == 7 && store->residentData()[1] == 11);

    ReadBuffer<size_type> reader;
    reader.open(store, 3 * sizeof(size_type));
    require(reader.isOpen());
    require(reader[0] == 7 && reader[1] == 11);
    reader.seek(1); require(reader[1] == 11);
    reader.close();

    store->append(13);
    require(store->resident() && store->capacity() == 4);
    require(budget.stats().maximum == 6 * sizeof(size_type));
    size_type values[3] = { 0, 0, 0 };
    store->readAt(0, values, 3);
    require(values[0] == 7 && values[1] == 11 && values[2] == 13);
    store.reset();
    require(budget.stats().current == 0);
  }

  // A failed growth admission spills lazily to the historic raw stream. The
  // reader keeps the store alive while it uses the regular asynchronous file
  // path. A bulk append remains one contiguous raw store write.
  {
    MemoryBudget budget(3 * sizeof(size_type));
    std::shared_ptr<AdaptiveRecordStore<size_type>> store(
      new AdaptiveRecordStore<size_type>(budget, "adaptive-test-spill"));
    require(store->tryReserveRecords(2));
    store->append(3); store->append(5);
    size_type tail[3] = { 8, 13, 21 };
    store->append(tail, 3);
    require(store->spilled() && store->count() == 5);
    const std::string spill_name = store->filename();
    require(exists(spill_name));
    require(store->at(0) == 3 && store->at(1) == 5 && store->at(2) == 8 &&
      store->at(3) == 13 && store->at(4) == 21);

    ReadBuffer<size_type> reader;
    reader.open(store, 2 * sizeof(size_type));
    require(reader[4] == 21);
    reader.seek(0); require(reader[0] == 3);
    reader.close();

    // A named checkpoint must work after overflow: it copies the complete
    // raw payload out of the owned temporary and leaves its lifecycle with
    // the checkpoint owner.
    const std::string checkpoint = std::string(root) + "/overflow-checkpoint.raw";
    store->materialize(checkpoint);
    require(store->filename() == checkpoint);
    require(store->fileOwnership() ==
      AdaptiveRecordStore<size_type>::FileOwnership::BORROWED);
    require(!exists(spill_name) && exists(checkpoint));
    require(store->releaseMaterializedFile() == checkpoint);
    store.reset();
    require(exists(checkpoint));
    require(std::remove(checkpoint.c_str()) == 0);
    require(budget.stats().current == 0);
  }

  // An explicit named materialization is borrowed: destroying the store must
  // not delete a path whose lifetime belongs to a caller or checkpoint.
  {
    const std::string filename = std::string(root) + "/borrowed.raw";
    MemoryBudget budget(4 * sizeof(size_type));
    {
      AdaptiveRecordStore<size_type> store(budget, "adaptive-test-borrowed");
      require(store.tryReserveRecords(3));
      store.append(17); store.append(19); store.append(23);
      store.materialize(filename);
      require(store.spilled());
      require(store.fileOwnership() ==
        AdaptiveRecordStore<size_type>::FileOwnership::BORROWED);
    }
    require(exists(filename));
    std::ifstream input(filename.c_str(), std::ios_base::binary);
    size_type values[3] = { 0, 0, 0 };
    input.read(reinterpret_cast<char*>(values), sizeof(values));
    require(static_cast<bool>(input));
    require(values[0] == 17 && values[1] == 19 && values[2] == 23);
    require(std::remove(filename.c_str()) == 0);
    require(budget.stats().current == 0);
  }

  // Failed pre-admission itself remains file-free; only append overflow is
  // allowed to choose the temporary spill route.
  {
    MemoryBudget budget(sizeof(size_type));
    AdaptiveRecordStore<size_type> store(budget, "adaptive-test-decline");
    require(!store.tryReserveRecords(2));
    require(store.resident() && store.count() == 0);
  }

  // Releasing peak scratch leaves its retained portion accurately charged.
  {
    MemoryBudget budget(10 * sizeof(size_type));
    MemoryBudget::Reservation peak = budget.tryReserve(8 * sizeof(size_type),
      "adaptive-test-peak");
    require(static_cast<bool>(peak));
    peak.shrinkTo(3 * sizeof(size_type));
    require(peak.bytes() == 3 * sizeof(size_type));
    require(budget.stats().current == 3 * sizeof(size_type));
    peak.shrinkTo(0);
    require(!static_cast<bool>(peak));
    require(budget.stats().current == 0);
  }

  TempFile::setDirectory(TempFile::DEFAULT_TEMP_DIR);
  require(::rmdir(root) == 0);
  return 0;
}
