#ifndef GCSA_COMPRESSED_BLOCK_H
#define GCSA_COMPRESSED_BLOCK_H

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

namespace gcsa
{

// A portable temporary-file stream. Records are never split between blocks.
class CompressedBlockWriter
{
public:
  enum Mode { RAW, ZSTD };

  // The estimate includes one uncompressed block, one worst-case zstd output
  // buffer, and the workspace the encoder will create: a single compression
  // context when the block holds one zstd job, and otherwise one context per
  // engaged worker plus the multi-threaded round buffer and job-output pool.
  // A block shorter than one job engages one worker whatever is requested, so
  // the writer configures and this function charges only the usable workers.
  // It can be used to reserve memory before construction.
  static std::size_t workingMemoryEstimate(std::size_t block_bytes, Mode mode,
    int zstd_level = 3, std::size_t compression_workers = 1);

  // Conservative peak filesystem bytes while committing one framed stream.
  // The bound includes the construction-only index sidecar while the same
  // index is being appended to the final artifact. Every logical record must
  // be at most maximum_record_bytes and is kept within one block.
  static std::uint64_t maximumTemporaryBytes(std::uint64_t logical_bytes,
    std::size_t maximum_record_bytes, std::size_t block_bytes);

  CompressedBlockWriter(const std::string& filename, std::size_t block_bytes,
    Mode mode = ZSTD, int zstd_level = 3,
    std::size_t compression_workers = 1);
  ~CompressedBlockWriter();
  CompressedBlockWriter(const CompressedBlockWriter&) = delete;
  CompressedBlockWriter& operator=(const CompressedBlockWriter&) = delete;

  // Zero-byte records are rejected: every record belongs to one stored block.
  void writeRecord(const void* data, std::size_t bytes);
  void finish();
  std::uint64_t records() const { return record_count; }
  std::uint64_t bytes() const { return byte_count; }

private:
  void flushBlock();
  std::string final_name, temporary_name, index_name;
  std::ofstream output, index_output;
  std::vector<std::uint8_t> buffer;
  std::size_t block_size;
  Mode requested_mode;
  int compression_level;
  std::size_t compression_workers;
  void* compression_context;
  std::uint64_t record_count, buffer_records, byte_count, whole_checksum;
  bool completed;
};

class CompressedBlockReader
{
public:
  explicit CompressedBlockReader(const std::string& filename);
  ~CompressedBlockReader();
  static bool isFramed(const std::string& filename);
  // Reads and validates only the fixed header; it does not allocate or decode
  // a data block. Schedulers use this before admitting a reader task.
  static std::uint64_t declaredBlockSize(const std::string& filename);
  // Validates the fixed header, footer, and index extent without allocating a
  // data block. This is safe to call before acquiring a reader reservation.
  static std::uint64_t declaredLogicalSize(const std::string& filename);
  static std::size_t workingMemoryEstimate(std::size_t block_bytes);
  std::size_t read(void* data, std::size_t bytes);
  std::size_t readAt(std::uint64_t offset, void* data, std::size_t bytes);
  void seekUncompressedByte(std::uint64_t offset);
  void seekBlock(std::uint64_t block);
  std::uint64_t blocks() const { return block_count; }
  std::uint64_t bytes() const { return logical_bytes; }
  std::uint64_t logicalSize() const { return logical_bytes; }
  std::uint64_t physicalSize() const { return physical_bytes; }
  std::uint64_t records() const { return record_count; }
  std::uint64_t blockSize() const { return block_size; }
  std::size_t indexMemoryBytes() const { return 0; }

private:
  void loadBlock(std::size_t block);
  void indexEntry(std::uint64_t block, std::uint64_t& physical,
    std::uint64_t& logical);
  void checkSequentialChecksum();
  void releaseConsumedCache(std::uint64_t physical_offset);
  std::ifstream input;
  std::vector<std::uint8_t> current;
  std::size_t current_block, current_offset;
  std::uint64_t index_offset, block_count, block_size, current_logical;
  std::uint64_t logical_bytes, physical_bytes, record_count, whole_checksum;
  std::uint64_t sequential_bytes, sequential_checksum;
  bool sequential;
  // Framed streams were the only readers in the codebase without the
  // release-behind discipline that workspace, external_sort, path_sort_run and
  // internal.h all share, so their clean pages stayed charged to the cgroup for
  // the life of the phase. The descriptor exists only to issue the advice; all
  // reading still goes through `input`.
  int cache_descriptor;
  std::uint64_t cache_released;
};

} // namespace gcsa

#endif // GCSA_COMPRESSED_BLOCK_H
