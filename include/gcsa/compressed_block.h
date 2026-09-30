#ifndef GCSA_COMPRESSED_BLOCK_H
#define GCSA_COMPRESSED_BLOCK_H

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

namespace gcsa
{

class CompressedBlockReader;

// A shared pool with a hard bound on pool-owned read/decode/cache allocations.
// The queue contains metadata only and each reader may have at most one queued,
// running, or cached block. Decode workspace is reserved before a worker
// starts and shrinks to the decoded block's actual allocation once the result
// enters the cache. Native thread stacks and small control objects are outside
// this data-buffer bound; callers separately reserve the reader's current block.
class CompressedBlockPrefetchPool
{
public:
  struct Stats
  {
    std::uint64_t submitted, completed, consumed, ready_hits;
    std::uint64_t waits, wait_nanoseconds, synchronous_blocks;
    std::uint64_t cancelled, errors, physical_bytes, decoded_bytes;
    std::size_t current_bytes, peak_bytes, worker_threads;

    Stats();
  };

  CompressedBlockPrefetchPool(std::size_t maximum_bytes,
    std::size_t worker_threads);
  ~CompressedBlockPrefetchPool();
  CompressedBlockPrefetchPool(const CompressedBlockPrefetchPool&) = delete;
  CompressedBlockPrefetchPool& operator=(
    const CompressedBlockPrefetchPool&) = delete;

  std::size_t memoryLimit() const;
  std::size_t usedBytes() const;
  std::size_t workers() const;
  Stats stats() const;

private:
  struct State;
  std::unique_ptr<State> state;

  std::uint64_t registerReader();
  bool submit(std::uint64_t owner, std::uint64_t block, int descriptor,
    std::uint64_t block_size, std::uint64_t physical,
    std::uint64_t logical, std::uint64_t next_physical,
    std::uint64_t next_logical);
  bool take(std::uint64_t owner, std::uint64_t block,
    std::vector<std::uint8_t>& result);
  void cancel(std::uint64_t owner);
  void recordSynchronous(std::uint64_t physical_bytes,
    std::uint64_t decoded_bytes);

  friend class CompressedBlockReader;
};

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
  // Opt-in before the first record, for payloads that will be checkpointed.
  // The stored-byte digest includes framing, compressed bytes, index and footer.
  void trackStoredChecksum();
  std::uint64_t storedChecksum() const;
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
  std::uint64_t stored_checksum;
  bool track_stored_checksum;
  bool completed;
};

class CompressedBlockReader
{
public:
  enum class FileAccess { PERSISTENT, TRANSIENT };
  explicit CompressedBlockReader(const std::string& filename);
  // Transient readers retain decoded blocks but open their immutable backing
  // file only for an actual block/index read. They do not use async prefetch.
  CompressedBlockReader(const std::string& filename, FileAccess access);
  // This form validates the stream metadata but queues block 0 instead of
  // decoding it synchronously. Readers sharing `pool` are serviced by the
  // same fixed workers and byte bound.
  CompressedBlockReader(const std::string& filename,
    const std::shared_ptr<CompressedBlockPrefetchPool>& pool);
  ~CompressedBlockReader();
  static bool isFramed(const std::string& filename);
  // Reads and validates only the fixed header; it does not allocate or decode
  // a data block. Schedulers use this before admitting a reader task.
  static std::uint64_t declaredBlockSize(const std::string& filename);
  // Validates the fixed header, footer, and index extent without allocating a
  // data block. This is safe to call before acquiring a reader reservation.
  static std::uint64_t declaredLogicalSize(const std::string& filename);
  static std::size_t workingMemoryEstimate(std::size_t block_bytes);
  // Peak reservation for one asynchronous job. The decoded result is charged
  // at its actual allocation after the read/decode temporaries are gone.
  static std::size_t prefetchWorkingMemoryEstimate(std::size_t block_bytes);
  std::size_t read(void* data, std::size_t bytes);
  std::size_t readAt(std::uint64_t offset, void* data, std::size_t bytes);
  // A sequential reader normally advises the kernel to drop the pages of
  // blocks it has consumed. Several readers scanning one stream at once must
  // not: the leader would evict what the others are about to read.
  void retainPageCache() { this->retain_cache = true; }
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
  CompressedBlockReader(const std::string& filename,
    const std::shared_ptr<CompressedBlockPrefetchPool>& pool,
    bool defer_first_block, bool transient_descriptor = false);
  struct DescriptorLease
  {
    explicit DescriptorLease(CompressedBlockReader& reader);
    ~DescriptorLease();
    DescriptorLease(const DescriptorLease&) = delete;
    DescriptorLease& operator=(const DescriptorLease&) = delete;
    CompressedBlockReader& reader;
    bool acquired;
  };
  void loadBlock(std::size_t block);
  void prefetchBlock(std::size_t block);
  void resetPrefetch();
  void blockExtent(std::size_t block, std::uint64_t& physical,
    std::uint64_t& logical, std::uint64_t& next_physical,
    std::uint64_t& next_logical);
  void indexEntry(std::uint64_t block, std::uint64_t& physical,
    std::uint64_t& logical);
  void checkSequentialChecksum();
  void releaseConsumedCache(std::uint64_t physical_offset);
  std::ifstream input;
  std::string source_name;
  bool transient_descriptor;
  std::vector<std::uint8_t> current;
  std::size_t current_block, current_offset;
  std::uint64_t index_offset, block_count, block_size, current_logical;
  std::uint64_t logical_bytes, physical_bytes, record_count, whole_checksum;
  std::uint64_t sequential_bytes, sequential_checksum;
  bool sequential;
  // One persistent descriptor serves synchronous pread(), pooled pread(), and
  // release-behind advice. Metadata is initially validated through `input`,
  // which is closed before this descriptor is opened, so a resident framed
  // stream still costs exactly one descriptor.
  int cache_descriptor;
  std::uint64_t cache_released;
  bool retain_cache;
  std::shared_ptr<CompressedBlockPrefetchPool> prefetch_pool;
  std::uint64_t prefetch_owner;
  bool extent_cached;
  std::size_t cached_extent_block;
  std::uint64_t cached_physical, cached_logical;
  std::uint64_t cached_next_physical, cached_next_logical;
};

} // namespace gcsa

#endif // GCSA_COMPRESSED_BLOCK_H
