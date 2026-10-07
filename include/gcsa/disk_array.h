/*
  Copyright (c) 2026 GCSA2 contributors

  A bounded block cache for mutable construction arrays.
*/

#ifndef GCSA_DISK_ARRAY_H
#define GCSA_DISK_ARRAY_H

#include <gcsa/resources.h>
#include <gcsa/utils.h>

#include <cstdint>
#include <string>
#include <vector>

namespace gcsa
{

/*
  DiskBackedArray64 is deliberately not a memory mapping. Only cache_bytes of
  array payload can become resident through this object, and every cache byte
  is covered by a MemoryBudget reservation. The direct-mapped replacement
  policy may perform extra I/O on an adversarial access pattern, but it cannot
  grow metadata or resident state with the logical array size.

  A newly created file is logically zero-filled with ftruncate(). Call flush()
  before a later phase consumes the array.
*/
class DiskBackedArray64
{
public:
  struct Stats
  {
    size_type block_reads, block_writes, replacements;

    Stats() : block_reads(0), block_writes(0), replacements(0) { }
  };

  constexpr static size_type DEFAULT_BLOCK_BYTES = MEGABYTE;

  // cache_bytes includes payload, tags, and dirty flags. The file is created
  // and truncated when create_file is true; otherwise its exact length is
  // validated before use.
  DiskBackedArray64(const std::string& filename, size_type elements,
    size_type cache_bytes, MemoryBudget& budget, bool create_file = true,
    size_type requested_block_bytes = DEFAULT_BLOCK_BYTES);
  ~DiskBackedArray64();

  size_type size() const { return this->elements_; }
  std::uint64_t get(size_type index);
  void set(size_type index, std::uint64_t value);

  // Flush dirty blocks. sync_file adds an fdatasync durability boundary.
  void flush(bool sync_file = true);
  const Stats& stats() const { return this->stats_; }

  static size_type minimumCacheBytes();

private:
  constexpr static std::uint64_t EMPTY_TAG = ~std::uint64_t(0);

  std::string filename_;
  int descriptor_;
  size_type elements_, array_bytes_, block_bytes_, block_elements_, slots_;
  std::vector<std::uint64_t> payload_, tags_;
  std::vector<std::uint8_t> dirty_;
  MemoryBudget::Reservation reservation_;
  Stats stats_;

  size_type slotFor(size_type block) const { return block % this->slots_; }
  std::uint64_t* load(size_type block);
  void writeSlot(size_type slot);

  DiskBackedArray64(const DiskBackedArray64&);
  DiskBackedArray64& operator=(const DiskBackedArray64&);
};

} // namespace gcsa

#endif // GCSA_DISK_ARRAY_H
