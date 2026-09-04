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

  // The estimate includes one uncompressed block and a conservative zstd output
  // buffer for every active compression worker, and can be used to reserve
  // memory before construction.
  static std::size_t workingMemoryEstimate(std::size_t block_bytes, Mode mode,
    int zstd_level = 3, std::size_t compression_workers = 1);

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
  std::ifstream input;
  std::vector<std::uint8_t> current;
  std::size_t current_block, current_offset;
  std::uint64_t index_offset, block_count, block_size, current_logical;
  std::uint64_t logical_bytes, physical_bytes, record_count, whole_checksum;
  std::uint64_t sequential_bytes, sequential_checksum;
  bool sequential;
};

} // namespace gcsa

#endif // GCSA_COMPRESSED_BLOCK_H
