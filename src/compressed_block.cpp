#include <gcsa/compressed_block.h>
#include <gcsa/internal.h>

#define ZSTD_STATIC_LINKING_ONLY
#include <zstd.h>

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <condition_variable>
#include <cstring>
#include <deque>
#include <exception>
#include <fcntl.h>
#include <limits>
#include <map>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <unistd.h>

namespace gcsa
{

namespace
{

constexpr std::uint64_t FILE_MAGIC = 0x314b4c4241534347ULL;
constexpr std::uint64_t FOOTER_MAGIC = 0x3152544f4f464347ULL;
constexpr std::uint64_t FNV_OFFSET = 1469598103934665603ULL;
constexpr std::uint32_t BLOCK_MAGIC = 0x314b4c42U;
constexpr std::uint32_t FORMAT_VERSION = 2;
// Same release-behind thresholds the other bounded streams use
// (path_sort_run.cpp:32-33, external_sort.cpp, internal.h).
constexpr std::uint64_t BLOCK_CACHE_TAIL_BYTES = 64 * MEGABYTE;
constexpr std::uint64_t BLOCK_CACHE_FLUSH_BYTES = 512 * MEGABYTE;

/*
  RETIRED, NOT REMOVED: re-hashing every delivered byte for a whole-stream
  digest.

  loadBlock already verifies each decoded block against its own stored checksum
  and cross-checks block extents against the on-disk index, so every byte handed
  to a caller has been verified once. read() then hashed the same bytes a second
  time to reproduce the footer's whole-stream digest. checksum() is byte-serial
  FNV-1a with the multiply on the loop-carried dependency chain, roughly 4-5
  cycles per byte or ~0.7 GB/s, and the chr21 merge alone reads 104.53 GB
  through this path -- about 300 s of hashing of which half is the duplicate.
  Every other framed reader in the build pays the same rate.

  What the second pass detected and the first does not: a reordering of whole
  blocks that preserved every individual block checksum AND the index. Nothing
  in the fork produces such a file; it would take a deliberate edit.

  Set VERIFY_WHOLE_STREAM_CHECKSUM to true to restore it. Do NOT substitute a
  faster digest: the value is stored in the footer, so changing it changes the
  artifact format and breaks resume for committed workspaces.
*/
constexpr bool VERIFY_WHOLE_STREAM_CHECKSUM = false;
constexpr std::uint32_t HEADER_SIZE = 32;
constexpr std::uint32_t FOOTER_SIZE = 56;
constexpr std::uint32_t RAW_BLOCK = 0;
constexpr std::uint32_t ZSTD_BLOCK = 1;
constexpr std::uint64_t INDEX_ENTRY_SIZE = 16;
constexpr std::uint64_t BLOCK_HEADER_SIZE = 40;

std::uint64_t
checksum(const std::uint8_t* data, std::size_t bytes,
  std::uint64_t value = FNV_OFFSET)
{
  for(std::size_t i = 0; i < bytes; ++i)
  {
    value = (value ^ data[i]) * 1099511628211ULL;
  }
  return value;
}

template<class Value>
void
writeLittle(std::ostream& output, Value value)
{
  std::uint8_t bytes[sizeof(Value)];
  for(std::size_t i = 0; i < sizeof(Value); ++i)
  {
    bytes[i] = static_cast<std::uint8_t>(value >> (8 * i));
  }
  output.write(reinterpret_cast<const char*>(bytes), sizeof(bytes));
  DiskIO::write_volume += sizeof(bytes);
  if(!output) { throw std::runtime_error("compressed block: write failure"); }
}

template<class Value>
Value
readLittle(std::istream& input, const char* context)
{
  std::uint8_t bytes[sizeof(Value)];
  input.read(reinterpret_cast<char*>(bytes), sizeof(bytes));
  DiskIO::read_volume += static_cast<std::size_t>(input.gcount());
  if(input.gcount() != static_cast<std::streamsize>(sizeof(bytes)))
  {
    throw std::runtime_error(std::string("compressed block: truncated ") + context);
  }
  Value value = 0;
  for(std::size_t i = 0; i < sizeof(Value); ++i)
  {
    value |= static_cast<Value>(bytes[i]) << (8 * i);
  }
  return value;
}

std::uint64_t
streamPosition(std::ostream& output)
{
  std::streampos position = output.tellp();
  if(position < 0) { throw std::runtime_error("compressed block: tell failure"); }
  return static_cast<std::uint64_t>(position);
}

void
preadExactCounted(int descriptor, void* data, std::size_t bytes,
  std::uint64_t offset, const char* context)
{
  if(descriptor < 0 ||
     offset > static_cast<std::uint64_t>(std::numeric_limits<off_t>::max()) ||
     bytes > static_cast<std::uint64_t>(std::numeric_limits<off_t>::max()) - offset)
  {
    throw std::runtime_error(std::string("compressed block: invalid ") + context +
      " offset");
  }
  std::uint8_t* output = static_cast<std::uint8_t*>(data);
  std::size_t done = 0;
  while(done < bytes)
  {
    ssize_t result = ::pread(descriptor, output + done, bytes - done,
      static_cast<off_t>(offset + done));
    if(result < 0 && errno == EINTR) { continue; }
    if(result <= 0)
    {
      throw std::runtime_error(std::string("compressed block: truncated ") + context);
    }
    done += static_cast<std::size_t>(result);
    DiskIO::read_volume += static_cast<size_type>(result);
  }
}

template<class Value>
Value
readLittleMemory(const std::uint8_t* bytes)
{
  Value value = 0;
  for(std::size_t i = 0; i < sizeof(Value); ++i)
  {
    value |= static_cast<Value>(bytes[i]) << (8 * i);
  }
  return value;
}

std::string
uniqueTemporaryName(const std::string& filename, const char* suffix)
{
  std::string pattern = filename + suffix + ".XXXXXX";
  std::vector<char> writable(pattern.begin(), pattern.end());
  writable.push_back(0);
  int descriptor = mkstemp(writable.data());
  if(descriptor < 0)
  {
    throw std::runtime_error("compressed block: temp create failure");
  }
  if(::close(descriptor) != 0)
  {
    ::unlink(writable.data());
    throw std::runtime_error("compressed block: temp close failure");
  }
  return writable.data();
}

void
syncFile(const std::string& filename)
{
  int descriptor = ::open(filename.c_str(), O_RDONLY);
  if(descriptor < 0)
  {
    throw std::runtime_error("compressed block: file sync open failure");
  }
  int result = ::fdatasync(descriptor);
  int close_result = ::close(descriptor);
  if(result != 0 || close_result != 0)
  {
    throw std::runtime_error("compressed block: fdatasync failure");
  }
}

void
syncDirectory(const std::string& filename)
{
  std::string::size_type slash = filename.find_last_of('/');
  std::string directory = (slash == std::string::npos ? "." :
    (slash == 0 ? "/" : filename.substr(0, slash)));
  int descriptor = ::open(directory.c_str(), O_RDONLY | O_DIRECTORY);
  if(descriptor < 0)
  {
    throw std::runtime_error("compressed block: directory sync open failure");
  }
  int result = ::fsync(descriptor);
  int close_result = ::close(descriptor);
  if(result != 0 || close_result != 0)
  {
    throw std::runtime_error("compressed block: directory sync failure");
  }
}

std::size_t
estimateSum(std::size_t left, std::size_t right)
{
  if(left > std::numeric_limits<std::size_t>::max() - right)
  {
    throw std::overflow_error("compressed block worker memory estimate overflow");
  }
  return left + right;
}

std::size_t
estimateProduct(std::size_t count, std::size_t bytes)
{
  if(bytes != 0 && count > std::numeric_limits<std::size_t>::max() / bytes)
  {
    throw std::overflow_error("compressed block worker memory estimate overflow");
  }
  return count * bytes;
}

// A multi-threaded encoder cuts one compression call into jobs of
// 1 << max(20, windowLog + 2) bytes, bounded by ZSTD_c_jobSize. The window is
// the one ZSTD_getCParams() derives for this level and this block, so the job
// size is known before any context exists.
std::size_t
compressionJobBytes(std::size_t block_bytes, int zstd_level)
{
  ZSTD_compressionParameters parameters = ZSTD_getCParams(zstd_level,
    static_cast<unsigned long long>(block_bytes), 0);
  unsigned job_log = std::max(20u, parameters.windowLog + 2);
  if(job_log >= 8 * sizeof(std::size_t))
  {
    return std::numeric_limits<std::size_t>::max();
  }
  std::size_t job_bytes = static_cast<std::size_t>(1) << job_log;
  ZSTD_bounds bounds = ZSTD_cParam_getBounds(ZSTD_c_jobSize);
  if(!ZSTD_isError(bounds.error) && bounds.upperBound > 0 &&
     job_bytes > static_cast<std::size_t>(bounds.upperBound))
  {
    job_bytes = static_cast<std::size_t>(bounds.upperBound);
  }
  return job_bytes;
}

// No worker starts before a whole job exists, so a block holding one job keeps
// one worker busy however many were requested. The writer configures this
// count rather than the request: the surplus workers would only enlarge the
// job pool. This also subsumes ZSTDMT_JOBSIZE_MIN, below which zstd disables
// multi-threading itself, because a job is never smaller than one mebibyte.
std::size_t
effectiveCompressionWorkers(std::size_t block_bytes, int zstd_level,
  std::size_t workers)
{
  if(workers <= 1 || block_bytes == 0) { return 1; }
  std::size_t job_bytes = compressionJobBytes(block_bytes, zstd_level);
  return std::min(workers, 1 + (block_bytes - 1) / job_bytes);
}

} // namespace

std::size_t
CompressedBlockWriter::workingMemoryEstimate(std::size_t block_bytes,
  Mode mode, int zstd_level, std::size_t workers)
{
  if(block_bytes == 0 || workers == 0)
  {
    throw std::invalid_argument(
      "compressed block size and worker count must be nonzero");
  }
  if(mode != RAW && mode != ZSTD)
  {
    throw std::invalid_argument("invalid compressed block mode");
  }

  std::size_t output_bytes = (mode == ZSTD ?
    ZSTD_compressBound(block_bytes) : block_bytes);
  std::size_t context_bytes = (mode == ZSTD ?
    ZSTD_estimateCCtxSize(zstd_level) : 0);
  if(ZSTD_isError(context_bytes) ||
     output_bytes > std::numeric_limits<std::size_t>::max() - block_bytes ||
     context_bytes > std::numeric_limits<std::size_t>::max() -
       block_bytes - output_bytes)
  {
    throw std::overflow_error("compressed block memory estimate overflow");
  }

  // The writer holds one uncompressed block and one worst-case output buffer
  // whatever the worker count: zstd replicates neither per worker.
  std::size_t total = block_bytes + output_bytes + context_bytes;
  std::size_t active = (mode == ZSTD ?
    effectiveCompressionWorkers(block_bytes, zstd_level, workers) : 1);
  if(active > 1)
  {
    // A multi-threaded context adds exactly what ZSTDMT_initCStream_internal()
    // allocates: one compression context per worker, a round buffer of
    // workers + 3 jobs, and a pool of at most BUF_POOL_MAX_NB_BUFFERS,
    // 2 * workers + 3, job output buffers. All three follow the job size
    // rather than the block size, so they are charged from the level and the
    // block the encoder will actually be given.
    std::size_t job_bytes = compressionJobBytes(block_bytes, zstd_level);
    std::size_t job_output = ZSTD_compressBound(job_bytes);
    if(job_output < job_bytes)
    {
      throw std::overflow_error("compressed block worker memory estimate overflow");
    }
    total = estimateSum(total, estimateProduct(active - 1, context_bytes));
    total = estimateSum(total, estimateProduct(active + 3, job_bytes));
    total = estimateSum(total, estimateProduct(2 * active + 3, job_output));
  }
  return total;
}

std::uint64_t
CompressedBlockWriter::maximumTemporaryBytes(std::uint64_t logical_bytes,
  std::size_t maximum_record_bytes, std::size_t block_bytes)
{
  if(block_bytes == 0 || maximum_record_bytes == 0 ||
     maximum_record_bytes > block_bytes)
  {
    throw std::invalid_argument("invalid compressed block record bound");
  }

  // Before a non-final block is flushed, the next record did not fit. As that
  // record is at most maximum_record_bytes, every completed block contains at
  // least block_bytes - maximum_record_bytes + 1 logical bytes. This is much
  // tighter than assuming one block per record while remaining valid for the
  // variable-width rank records.
  const std::uint64_t minimum_fill =
    static_cast<std::uint64_t>(block_bytes - maximum_record_bytes) + 1;
  const std::uint64_t blocks = (logical_bytes == 0 ? 0 :
    1 + (logical_bytes - 1) / minimum_fill);
  const std::uint64_t per_block_peak = BLOCK_HEADER_SIZE + 2 * INDEX_ENTRY_SIZE;
  const std::uint64_t fixed = HEADER_SIZE + FOOTER_SIZE;
  if(blocks > (std::numeric_limits<std::uint64_t>::max() - fixed) /
       per_block_peak)
  {
    throw std::overflow_error("compressed block temporary size overflow");
  }
  const std::uint64_t metadata = fixed + blocks * per_block_peak;
  if(logical_bytes > std::numeric_limits<std::uint64_t>::max() - metadata)
  {
    throw std::overflow_error("compressed block temporary size overflow");
  }
  return logical_bytes + metadata;
}

CompressedBlockWriter::CompressedBlockWriter(const std::string& filename,
  std::size_t block_bytes, Mode mode, int zstd_level, std::size_t workers) :
  final_name(filename), temporary_name(), index_name(), output(), index_output(),
  buffer(), block_size(block_bytes), requested_mode(mode),
  compression_level(zstd_level), compression_workers(workers),
  compression_context(nullptr), record_count(0), buffer_records(0),
  byte_count(0), whole_checksum(FNV_OFFSET), completed(false)
{
  if(block_bytes == 0 || workers == 0 || (mode != RAW && mode != ZSTD) ||
     zstd_level < ZSTD_minCLevel() || zstd_level > ZSTD_maxCLevel())
  {
    throw std::invalid_argument("invalid compressed block configuration");
  }
  // Requesting more workers than this block can occupy allocates a larger job
  // pool without adding a job. Configure what workingMemoryEstimate() charges,
  // so the reservation and the encoder cannot drift apart.
  if(mode == ZSTD)
  {
    this->compression_workers =
      effectiveCompressionWorkers(block_bytes, zstd_level, workers);
  }

  try
  {
    this->temporary_name = uniqueTemporaryName(filename, ".partial");
    this->index_name = uniqueTemporaryName(filename, ".index");
    this->output.open(this->temporary_name,
      std::ios::binary | std::ios::trunc);
    this->index_output.open(this->index_name,
      std::ios::binary | std::ios::trunc);
    if(!this->output || !this->index_output)
    {
      throw std::runtime_error("compressed block: temporary open failure");
    }

    if(mode == ZSTD)
    {
      ZSTD_CCtx* context = ZSTD_createCCtx();
      if(context == nullptr)
      {
        throw std::runtime_error("compressed block: cannot create zstd context");
      }
      this->compression_context = context;
      std::size_t result = ZSTD_CCtx_setParameter(context,
        ZSTD_c_compressionLevel, zstd_level);
      if(ZSTD_isError(result))
      {
        throw std::runtime_error("compressed block: invalid zstd compression level");
      }
      if(this->compression_workers > 1)
      {
        result = ZSTD_CCtx_setParameter(context, ZSTD_c_nbWorkers,
          this->compression_workers);
        if(ZSTD_isError(result))
        {
          throw std::runtime_error(
            "compressed block: zstd library has no multithread support");
        }
      }
    }

    this->buffer.reserve(block_bytes);
    writeLittle(this->output, FILE_MAGIC);
    writeLittle(this->output, FORMAT_VERSION);
    writeLittle(this->output, HEADER_SIZE);
    writeLittle(this->output, static_cast<std::uint64_t>(block_bytes));
    writeLittle(this->output, static_cast<std::uint32_t>(mode));
    writeLittle(this->output,
      static_cast<std::uint32_t>(this->compression_workers));
  }
  catch(...)
  {
    if(this->compression_context != nullptr)
    {
      ZSTD_freeCCtx(static_cast<ZSTD_CCtx*>(this->compression_context));
      this->compression_context = nullptr;
    }
    if(!this->temporary_name.empty()) { ::unlink(this->temporary_name.c_str()); }
    if(!this->index_name.empty()) { ::unlink(this->index_name.c_str()); }
    throw;
  }
}

CompressedBlockWriter::~CompressedBlockWriter()
{
  this->output.close();
  this->index_output.close();
  if(this->compression_context != nullptr)
  {
    ZSTD_freeCCtx(static_cast<ZSTD_CCtx*>(this->compression_context));
  }
  if(!this->completed && !this->temporary_name.empty())
  {
    ::unlink(this->temporary_name.c_str());
  }
  if(!this->index_name.empty()) { ::unlink(this->index_name.c_str()); }
}

void
CompressedBlockWriter::writeRecord(const void* data, std::size_t bytes)
{
  if(this->completed)
  {
    throw std::logic_error("compressed block writer already finished");
  }
  if(data == nullptr || bytes == 0 || bytes > this->block_size)
  {
    throw std::invalid_argument("invalid compressed block record");
  }
  if(bytes > this->block_size - this->buffer.size()) { this->flushBlock(); }
  if(this->record_count == std::numeric_limits<std::uint64_t>::max() ||
     bytes > std::numeric_limits<std::uint64_t>::max() - this->byte_count)
  {
    throw std::overflow_error("compressed block record count overflow");
  }

  const std::uint8_t* begin = static_cast<const std::uint8_t*>(data);
  this->buffer.insert(this->buffer.end(), begin, begin + bytes);
  ++this->record_count;
  ++this->buffer_records;
  this->byte_count += bytes;
  this->whole_checksum = checksum(begin, bytes, this->whole_checksum);
}

void
CompressedBlockWriter::flushBlock()
{
  if(this->buffer.empty()) { return; }
  writeLittle(this->index_output, streamPosition(this->output));
  writeLittle(this->index_output,
    this->byte_count - static_cast<std::uint64_t>(this->buffer.size()));

  std::vector<std::uint8_t> compressed;
  std::uint32_t codec = RAW_BLOCK;
  if(this->requested_mode == ZSTD)
  {
    compressed.resize(ZSTD_compressBound(this->buffer.size()));
    std::size_t result = ZSTD_compress2(
      static_cast<ZSTD_CCtx*>(this->compression_context),
      compressed.data(), compressed.size(),
      this->buffer.data(), this->buffer.size());
    if(ZSTD_isError(result))
    {
      throw std::runtime_error("compressed block: zstd compression failed");
    }
    compressed.resize(result);
    if(compressed.size() < this->buffer.size()) { codec = ZSTD_BLOCK; }
    else { compressed.clear(); }
  }
  const std::vector<std::uint8_t>& payload =
    (codec == ZSTD_BLOCK ? compressed : this->buffer);

  writeLittle(this->output, BLOCK_MAGIC);
  writeLittle(this->output, codec);
  writeLittle(this->output, this->buffer_records);
  writeLittle(this->output,
    static_cast<std::uint64_t>(this->buffer.size()));
  writeLittle(this->output,
    static_cast<std::uint64_t>(payload.size()));
  writeLittle(this->output,
    checksum(this->buffer.data(), this->buffer.size()));
  this->output.write(reinterpret_cast<const char*>(payload.data()),
    payload.size());
  DiskIO::write_volume += payload.size();
  if(!this->output)
  {
    throw std::runtime_error("compressed block: payload write failure");
  }
  this->buffer.clear();
  this->buffer_records = 0;
}

void
CompressedBlockWriter::finish()
{
  if(this->completed) { return; }
  this->flushBlock();
  this->index_output.flush();
  if(!this->index_output)
  {
    throw std::runtime_error("compressed block: index flush failure");
  }
  syncFile(this->index_name);

  std::ifstream index_input(this->index_name, std::ios::binary);
  if(!index_input)
  {
    throw std::runtime_error("compressed block: index reopen failure");
  }
  std::uint64_t index_offset = streamPosition(this->output);
  std::uint64_t blocks = 0;
  std::uint8_t entry[INDEX_ENTRY_SIZE];
  while(true)
  {
    index_input.read(reinterpret_cast<char*>(entry), sizeof(entry));
    std::streamsize bytes = index_input.gcount();
    if(bytes == 0)
    {
      if(index_input.eof()) { break; }
      throw std::runtime_error("compressed block: index read failure");
    }
    if(bytes != static_cast<std::streamsize>(sizeof(entry)))
    {
      throw std::runtime_error("compressed block: truncated temporary index");
    }
    this->output.write(reinterpret_cast<const char*>(entry), sizeof(entry));
    DiskIO::read_volume += static_cast<std::size_t>(bytes);
    DiskIO::write_volume += sizeof(entry);
    if(!this->output)
    {
      throw std::runtime_error("compressed block: index write failure");
    }
    ++blocks;
  }
  index_input.close();
  this->index_output.close();

  writeLittle(this->output, FOOTER_MAGIC);
  writeLittle(this->output, FORMAT_VERSION);
  writeLittle(this->output, FOOTER_SIZE);
  writeLittle(this->output, blocks);
  writeLittle(this->output, this->record_count);
  writeLittle(this->output, this->byte_count);
  writeLittle(this->output, this->whole_checksum);
  writeLittle(this->output, index_offset);
  this->output.flush();
  this->output.close();
  if(this->output.fail())
  {
    throw std::runtime_error("compressed block: output close failure");
  }

  // The sidecar is construction-only. Removing it before publication means a
  // successfully committed artifact never leaves a permanent companion file.
  if(::unlink(this->index_name.c_str()) != 0)
  {
    throw std::runtime_error("compressed block: index cleanup failure");
  }
  this->index_name.clear();
  syncFile(this->temporary_name);
  if(::rename(this->temporary_name.c_str(), this->final_name.c_str()) != 0)
  {
    throw std::runtime_error("compressed block: rename failure");
  }
  syncDirectory(this->final_name);
  this->completed = true;
}

bool
CompressedBlockReader::isFramed(const std::string& filename)
{
  std::ifstream input(filename, std::ios::binary);
  try
  {
    return input && readLittle<std::uint64_t>(input, "header") == FILE_MAGIC &&
      readLittle<std::uint32_t>(input, "header") == FORMAT_VERSION &&
      readLittle<std::uint32_t>(input, "header") == HEADER_SIZE;
  }
  catch(...)
  {
    return false;
  }
}

std::uint64_t
CompressedBlockReader::declaredBlockSize(const std::string& filename)
{
  std::ifstream input(filename, std::ios::binary);
  if(!input || readLittle<std::uint64_t>(input, "header") != FILE_MAGIC ||
     readLittle<std::uint32_t>(input, "header") != FORMAT_VERSION ||
     readLittle<std::uint32_t>(input, "header") != HEADER_SIZE)
  {
    throw std::runtime_error("compressed block: invalid header");
  }
  std::uint64_t result = readLittle<std::uint64_t>(input, "header");
  if(result == 0)
  {
    throw std::runtime_error("compressed block: invalid block size");
  }
  return result;
}

std::uint64_t
CompressedBlockReader::declaredLogicalSize(const std::string& filename)
{
  std::ifstream input(filename, std::ios::binary);
  if(!input || readLittle<std::uint64_t>(input, "header") != FILE_MAGIC ||
     readLittle<std::uint32_t>(input, "header") != FORMAT_VERSION ||
     readLittle<std::uint32_t>(input, "header") != HEADER_SIZE)
  {
    throw std::runtime_error("compressed block: invalid header");
  }
  std::uint64_t block_size = readLittle<std::uint64_t>(input, "header");
  std::uint32_t mode = readLittle<std::uint32_t>(input, "header");
  readLittle<std::uint32_t>(input, "header");
  if(block_size == 0 ||
     (mode != CompressedBlockWriter::RAW && mode != CompressedBlockWriter::ZSTD))
  {
    throw std::runtime_error("compressed block: invalid header values");
  }

  input.seekg(0, std::ios::end);
  std::streamoff end = input.tellg();
  if(end < 0 || end < static_cast<std::streamoff>(HEADER_SIZE + FOOTER_SIZE))
  {
    throw std::runtime_error("compressed block: truncated footer");
  }
  input.seekg(end - FOOTER_SIZE);
  if(readLittle<std::uint64_t>(input, "footer") != FOOTER_MAGIC ||
     readLittle<std::uint32_t>(input, "footer") != FORMAT_VERSION ||
     readLittle<std::uint32_t>(input, "footer") != FOOTER_SIZE)
  {
    throw std::runtime_error("compressed block: invalid footer");
  }
  std::uint64_t blocks = readLittle<std::uint64_t>(input, "footer");
  readLittle<std::uint64_t>(input, "footer"); // Record count.
  std::uint64_t logical_bytes = readLittle<std::uint64_t>(input, "footer");
  readLittle<std::uint64_t>(input, "footer"); // Whole-stream checksum.
  std::uint64_t index_offset = readLittle<std::uint64_t>(input, "footer");
  if(blocks > std::numeric_limits<std::uint64_t>::max() / INDEX_ENTRY_SIZE)
  {
    throw std::runtime_error("compressed block: invalid offset index size");
  }
  std::uint64_t index_bytes = blocks * INDEX_ENTRY_SIZE;
  if(index_offset < HEADER_SIZE ||
     index_offset > std::numeric_limits<std::uint64_t>::max() - index_bytes ||
     index_offset + index_bytes != static_cast<std::uint64_t>(end - FOOTER_SIZE))
  {
    throw std::runtime_error("compressed block: invalid offset index");
  }
  return logical_bytes;
}

std::size_t
CompressedBlockReader::workingMemoryEstimate(std::size_t block_bytes)
{
  std::size_t context_bytes = ZSTD_estimateDCtxSize();
  if(block_bytes == 0 || ZSTD_isError(context_bytes) ||
     block_bytes > std::numeric_limits<std::size_t>::max() / 2 ||
     context_bytes > std::numeric_limits<std::size_t>::max() - 2 * block_bytes)
  {
    throw std::invalid_argument("invalid compressed block size");
  }
  return 2 * block_bytes + context_bytes;
}

std::size_t
CompressedBlockReader::prefetchWorkingMemoryEstimate(std::size_t block_bytes)
{
  // The fixed header is read into a stack buffer. The reservation covers the
  // encoded payload, decoded payload, and zstd's one-shot decoder context.
  return CompressedBlockReader::workingMemoryEstimate(block_bytes);
}

CompressedBlockPrefetchPool::Stats::Stats() :
  submitted(0), completed(0), consumed(0), ready_hits(0),
  waits(0), wait_nanoseconds(0), synchronous_blocks(0),
  cancelled(0), errors(0), physical_bytes(0), decoded_bytes(0),
  current_bytes(0), peak_bytes(0), worker_threads(0)
{
}

struct CompressedBlockPrefetchPool::State
{
  enum Status { QUEUED, RUNNING, READY, FAILED, CANCELLED };

  struct Entry
  {
    std::uint64_t owner, block, block_size;
    int descriptor;
    std::uint64_t physical, logical, next_physical, next_logical;
    std::size_t reservation;
    Status status;
    bool cancel_requested, cancellation_counted, physical_accounted;
    std::uint64_t physical_read;
    std::vector<std::uint8_t> data;
    std::exception_ptr error;

    Entry(std::uint64_t reader, std::uint64_t block_id, int input,
      std::uint64_t declared_block_size, std::uint64_t physical_offset,
      std::uint64_t logical_offset, std::uint64_t next_physical_offset,
      std::uint64_t next_logical_offset, std::size_t bytes) :
      owner(reader), block(block_id), block_size(declared_block_size),
      descriptor(input), physical(physical_offset), logical(logical_offset),
      next_physical(next_physical_offset), next_logical(next_logical_offset),
      reservation(bytes), status(QUEUED), cancel_requested(false),
      cancellation_counted(false), physical_accounted(false),
      physical_read(0), data(), error()
    {
    }
  };

  std::size_t limit, used, peak;
  std::uint64_t next_owner;
  bool stopping;
  Stats counters;
  std::deque<std::shared_ptr<Entry>> queue;
  std::map<std::uint64_t, std::shared_ptr<Entry>> owners;
  std::vector<std::thread> worker_threads;
  mutable std::mutex mutex;
  std::condition_variable work, changed;

  State(std::size_t maximum_bytes, std::size_t workers) :
    limit(maximum_bytes), used(0), peak(0), next_owner(1), stopping(false),
    counters(), queue(), owners(), worker_threads(), mutex(), work(), changed()
  {
    this->counters.worker_threads = workers;
    try
    {
      for(std::size_t i = 0; i < workers; ++i)
      {
        this->worker_threads.push_back(std::thread(&State::run, this));
      }
    }
    catch(...)
    {
      {
        std::lock_guard<std::mutex> lock(this->mutex);
        this->stopping = true; this->work.notify_all();
      }
      for(std::thread& worker : this->worker_threads)
      {
        if(worker.joinable()) { worker.join(); }
      }
      throw;
    }
  }

  static void preadExact(int descriptor, void* target, std::size_t bytes,
    std::uint64_t offset, std::uint64_t& physical_read)
  {
    std::uint8_t* output = static_cast<std::uint8_t*>(target);
    std::size_t done = 0;
    while(done < bytes)
    {
      ssize_t result = ::pread(descriptor, output + done, bytes - done,
        static_cast<off_t>(offset + done));
      if(result < 0 && errno == EINTR) { continue; }
      if(result <= 0)
      {
        throw std::runtime_error("compressed block prefetch: truncated block");
      }
      done += static_cast<std::size_t>(result);
      physical_read += static_cast<std::uint64_t>(result);
    }
  }

  template<class Value>
  static Value readMemory(const std::uint8_t*& cursor,
    const std::uint8_t* limit)
  {
    if(static_cast<std::size_t>(limit - cursor) < sizeof(Value))
    {
      throw std::runtime_error("compressed block prefetch: truncated header");
    }
    Value result = 0;
    for(std::size_t i = 0; i < sizeof(Value); ++i)
    {
      result |= static_cast<Value>(cursor[i]) << (8 * i);
    }
    cursor += sizeof(Value);
    return result;
  }

  static std::vector<std::uint8_t> decode(const Entry& entry,
    std::uint64_t& physical_read)
  {
    if(entry.descriptor < 0 || entry.next_physical <= entry.physical ||
       entry.next_logical <= entry.logical ||
       entry.next_physical - entry.physical < BLOCK_HEADER_SIZE)
    {
      throw std::runtime_error("compressed block prefetch: invalid extent");
    }

    std::uint8_t header[BLOCK_HEADER_SIZE];
    preadExact(entry.descriptor, header, sizeof(header), entry.physical,
      physical_read);
    const std::uint8_t* cursor = header;
    const std::uint8_t* header_limit = header + sizeof(header);
    std::uint32_t magic = readMemory<std::uint32_t>(cursor, header_limit);
    std::uint32_t codec = readMemory<std::uint32_t>(cursor, header_limit);
    std::uint64_t records = readMemory<std::uint64_t>(cursor, header_limit);
    std::uint64_t raw_bytes = readMemory<std::uint64_t>(cursor, header_limit);
    std::uint64_t stored_bytes = readMemory<std::uint64_t>(cursor, header_limit);
    std::uint64_t expected_checksum =
      readMemory<std::uint64_t>(cursor, header_limit);

    bool valid_codec_extent =
      (codec == RAW_BLOCK && stored_bytes == raw_bytes) ||
      (codec == ZSTD_BLOCK && stored_bytes < raw_bytes);
    if(magic != BLOCK_MAGIC || !valid_codec_extent || records == 0 ||
       raw_bytes == 0 || raw_bytes > entry.block_size || records > raw_bytes ||
       raw_bytes > std::numeric_limits<std::size_t>::max() ||
       stored_bytes > std::numeric_limits<std::size_t>::max() ||
       stored_bytes != entry.next_physical - entry.physical - BLOCK_HEADER_SIZE ||
       raw_bytes != entry.next_logical - entry.logical)
    {
      throw std::runtime_error("compressed block prefetch: invalid block header");
    }

    std::vector<std::uint8_t> encoded(static_cast<std::size_t>(stored_bytes));
    preadExact(entry.descriptor, encoded.data(), encoded.size(),
      entry.physical + BLOCK_HEADER_SIZE, physical_read);
    std::vector<std::uint8_t> decoded(static_cast<std::size_t>(raw_bytes));
    if(codec == RAW_BLOCK)
    {
      decoded.swap(encoded);
    }
    else
    {
      std::size_t result = ZSTD_decompress(decoded.data(), decoded.size(),
        encoded.data(), encoded.size());
      if(ZSTD_isError(result) || result != decoded.size())
      {
        throw std::runtime_error("compressed block prefetch: zstd decompression failed");
      }
    }
    if(checksum(decoded.data(), decoded.size()) != expected_checksum)
    {
      throw std::runtime_error("compressed block prefetch: block checksum failure");
    }
    return decoded;
  }

  bool runnable() const
  {
    for(const std::shared_ptr<Entry>& entry : this->queue)
    {
      if(!(entry->cancel_requested) && entry->status == QUEUED &&
         entry->reservation <= this->limit - this->used)
      {
        return true;
      }
    }
    return false;
  }

  std::deque<std::shared_ptr<Entry>>::iterator nextRunnable()
  {
    return std::find_if(this->queue.begin(), this->queue.end(),
      [&](const std::shared_ptr<Entry>& entry)
      {
        return !(entry->cancel_requested) && entry->status == QUEUED &&
          entry->reservation <= this->limit - this->used;
      });
  }

  void countCancellation(const std::shared_ptr<Entry>& entry)
  {
    if(!(entry->cancellation_counted))
    {
      entry->cancellation_counted = true;
      ++this->counters.cancelled;
    }
  }

  std::uint64_t claimPhysical(const std::shared_ptr<Entry>& entry)
  {
    if(entry->physical_accounted) { return 0; }
    entry->physical_accounted = true;
    return entry->physical_read;
  }

  void run()
  {
    for(;;)
    {
      std::shared_ptr<Entry> entry;
      {
        std::unique_lock<std::mutex> lock(this->mutex);
        this->work.wait(lock,
          [&]() { return this->stopping || this->runnable(); });
        if(this->stopping) { return; }
        auto next = this->nextRunnable();
        if(next == this->queue.end()) { continue; }
        entry = *next; this->queue.erase(next);
        entry->status = RUNNING;
        this->used += entry->reservation;
        this->peak = std::max(this->peak, this->used);
      }

      std::vector<std::uint8_t> decoded;
      std::exception_ptr error;
      std::uint64_t physical_read = 0;
      try
      {
        decoded = State::decode(*entry, physical_read);
      }
      catch(...)
      {
        error = std::current_exception();
      }

      std::unique_lock<std::mutex> lock(this->mutex);
      entry->physical_read = physical_read;
      this->counters.physical_bytes += physical_read;
      if(!error)
      {
        ++this->counters.completed;
        this->counters.decoded_bytes += decoded.size();
      }

      if(entry->cancel_requested || this->stopping)
      {
        // Free the decoded allocation before releasing its reservation and
        // waking another worker. This keeps the byte limit true even during
        // cancellation.
        lock.unlock();
        std::vector<std::uint8_t>().swap(decoded);
        lock.lock();
        this->used -= entry->reservation;
        entry->reservation = 0; entry->status = CANCELLED;
        this->countCancellation(entry);
        auto found = this->owners.find(entry->owner);
        if(!this->stopping && found != this->owners.end() &&
           found->second == entry)
        {
          this->owners.erase(found);
        }
      }
      else if(error)
      {
        this->used -= entry->reservation;
        entry->reservation = 0; entry->status = FAILED;
        entry->error = error; ++this->counters.errors;
      }
      else
      {
        std::size_t cached_bytes = decoded.capacity();
        if(cached_bytes > entry->reservation)
        {
          this->used -= entry->reservation;
          entry->reservation = 0; entry->status = FAILED;
          entry->error = std::make_exception_ptr(std::runtime_error(
            "compressed block prefetch: decoded allocation exceeds reservation"));
          ++this->counters.errors;
        }
        else
        {
          this->used -= entry->reservation - cached_bytes;
          entry->reservation = cached_bytes;
          entry->data.swap(decoded); entry->status = READY;
        }
      }
      this->changed.notify_all(); this->work.notify_all();
      // Drop the worker's ownership while the state mutex is still held. The
      // queue/owner map retains live entries, and cancelled entries can be
      // destroyed here without racing a consumer's final shared_ptr release.
      entry.reset();
    }
  }

  void shutdown()
  {
    {
      std::lock_guard<std::mutex> lock(this->mutex);
      this->stopping = true;
      for(auto& item : this->owners)
      {
        item.second->cancel_requested = true;
      }
      this->work.notify_all(); this->changed.notify_all();
    }
    for(std::thread& worker : this->worker_threads)
    {
      if(worker.joinable()) { worker.join(); }
    }

    std::uint64_t physical = 0;
    {
      std::lock_guard<std::mutex> lock(this->mutex);
      for(auto& item : this->owners)
      {
        std::shared_ptr<Entry> entry = item.second;
        this->countCancellation(entry);
        physical += this->claimPhysical(entry);
        std::vector<std::uint8_t>().swap(entry->data);
      }
      this->owners.clear(); this->queue.clear(); this->used = 0;
    }
    DiskIO::read_volume += static_cast<size_type>(physical);
  }
};

CompressedBlockPrefetchPool::CompressedBlockPrefetchPool(
  std::size_t maximum_bytes, std::size_t worker_threads) : state()
{
  if(maximum_bytes == 0 || worker_threads == 0)
  {
    throw std::invalid_argument("invalid compressed block prefetch pool");
  }
  this->state.reset(new State(maximum_bytes, worker_threads));
}

CompressedBlockPrefetchPool::~CompressedBlockPrefetchPool()
{
  if(this->state) { this->state->shutdown(); }
}

std::size_t
CompressedBlockPrefetchPool::memoryLimit() const
{
  return this->state->limit;
}

std::size_t
CompressedBlockPrefetchPool::usedBytes() const
{
  std::lock_guard<std::mutex> lock(this->state->mutex);
  return this->state->used;
}

std::size_t
CompressedBlockPrefetchPool::workers() const
{
  return this->state->worker_threads.size();
}

CompressedBlockPrefetchPool::Stats
CompressedBlockPrefetchPool::stats() const
{
  std::lock_guard<std::mutex> lock(this->state->mutex);
  Stats result = this->state->counters;
  result.current_bytes = this->state->used;
  result.peak_bytes = this->state->peak;
  return result;
}

std::uint64_t
CompressedBlockPrefetchPool::registerReader()
{
  std::lock_guard<std::mutex> lock(this->state->mutex);
  if(this->state->stopping || this->state->next_owner == 0)
  {
    throw std::runtime_error("compressed block prefetch: owner space exhausted");
  }
  return this->state->next_owner++;
}

bool
CompressedBlockPrefetchPool::submit(std::uint64_t owner, std::uint64_t block,
  int descriptor, std::uint64_t block_size, std::uint64_t physical,
  std::uint64_t logical, std::uint64_t next_physical,
  std::uint64_t next_logical)
{
  if(block_size > std::numeric_limits<std::size_t>::max()) { return false; }
  std::size_t reservation =
    CompressedBlockReader::prefetchWorkingMemoryEstimate(
      static_cast<std::size_t>(block_size));
  std::lock_guard<std::mutex> lock(this->state->mutex);
  if(this->state->stopping || descriptor < 0 || reservation > this->state->limit)
  {
    return false;
  }
  auto found = this->state->owners.find(owner);
  if(found != this->state->owners.end())
  {
    if(found->second->block == block) { return true; }
    throw std::logic_error("compressed block prefetch: reader queued two blocks");
  }
  std::shared_ptr<State::Entry> entry(new State::Entry(owner, block,
    descriptor, block_size, physical, logical, next_physical, next_logical,
    reservation));
  this->state->owners[owner] = entry;
  this->state->queue.push_back(entry);
  ++this->state->counters.submitted;
  this->state->work.notify_one();
  return true;
}

bool
CompressedBlockPrefetchPool::take(std::uint64_t owner, std::uint64_t block,
  std::vector<std::uint8_t>& result)
{
  std::uint64_t physical = 0;
  std::exception_ptr error;
  bool waited = false, ready = false;
  std::chrono::steady_clock::time_point wait_start;
  {
    std::unique_lock<std::mutex> lock(this->state->mutex);
    auto found = this->state->owners.find(owner);
    if(found == this->state->owners.end()) { return false; }
    std::shared_ptr<State::Entry> entry = found->second;
    if(entry->block != block)
    {
      throw std::logic_error("compressed block prefetch: unexpected block");
    }
    if(entry->status == State::QUEUED)
    {
      // A speculative job that has not started must never block the serial
      // merge behind cached work for another shard. Remove it and use the
      // reader's already-reserved synchronous workspace.
      auto queued = std::find(this->state->queue.begin(),
        this->state->queue.end(), entry);
      if(queued != this->state->queue.end()) { this->state->queue.erase(queued); }
      entry->cancel_requested = true; entry->status = State::CANCELLED;
      this->state->countCancellation(entry);
      this->state->owners.erase(found);
      this->state->changed.notify_all();
      return false;
    }
    ready = (entry->status == State::READY);
    if(entry->status == State::RUNNING)
    {
      waited = true; wait_start = std::chrono::steady_clock::now();
      this->state->changed.wait(lock,
        [&]() { return entry->status != State::RUNNING; });
    }
    if(waited)
    {
      std::uint64_t elapsed = static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(
          std::chrono::steady_clock::now() - wait_start).count());
      ++this->state->counters.waits;
      this->state->counters.wait_nanoseconds += elapsed;
    }
    if(entry->status == State::READY)
    {
      if(ready) { ++this->state->counters.ready_hits; }
      result.swap(entry->data);
      this->state->used -= entry->reservation;
      entry->reservation = 0;
      ++this->state->counters.consumed;
    }
    else if(entry->status == State::FAILED)
    {
      error = entry->error;
    }
    else
    {
      return false;
    }
    physical = this->state->claimPhysical(entry);
    this->state->owners.erase(owner);
    this->state->work.notify_all();
  }
  DiskIO::read_volume += static_cast<size_type>(physical);
  if(error) { std::rethrow_exception(error); }
  return true;
}

void
CompressedBlockPrefetchPool::cancel(std::uint64_t owner)
{
  if(owner == 0) { return; }
  std::shared_ptr<State::Entry> entry;
  std::uint64_t physical = 0;
  {
    std::unique_lock<std::mutex> lock(this->state->mutex);
    auto found = this->state->owners.find(owner);
    if(found == this->state->owners.end()) { return; }
    entry = found->second; entry->cancel_requested = true;
    if(entry->status == State::QUEUED)
    {
      auto queued = std::find(this->state->queue.begin(),
        this->state->queue.end(), entry);
      if(queued != this->state->queue.end()) { this->state->queue.erase(queued); }
      entry->status = State::CANCELLED;
      this->state->countCancellation(entry);
      this->state->owners.erase(found);
    }
    else if(entry->status == State::RUNNING)
    {
      this->state->changed.wait(lock, [&]()
      {
        return this->state->owners.find(owner) == this->state->owners.end();
      });
    }
    else
    {
      std::vector<std::uint8_t>().swap(entry->data);
      this->state->used -= entry->reservation;
      entry->reservation = 0; entry->status = State::CANCELLED;
      this->state->countCancellation(entry);
      this->state->owners.erase(found);
    }
    physical = this->state->claimPhysical(entry);
    this->state->work.notify_all(); this->state->changed.notify_all();
  }
  DiskIO::read_volume += static_cast<size_type>(physical);
}

void
CompressedBlockPrefetchPool::recordSynchronous(
  std::uint64_t physical_bytes, std::uint64_t decoded_bytes)
{
  std::lock_guard<std::mutex> lock(this->state->mutex);
  ++this->state->counters.synchronous_blocks;
  this->state->counters.physical_bytes += physical_bytes;
  this->state->counters.decoded_bytes += decoded_bytes;
}

CompressedBlockReader::CompressedBlockReader(const std::string& filename) :
  CompressedBlockReader(filename,
    std::shared_ptr<CompressedBlockPrefetchPool>(), false)
{
}

CompressedBlockReader::CompressedBlockReader(const std::string& filename,
  FileAccess access) :
  CompressedBlockReader(filename, std::shared_ptr<CompressedBlockPrefetchPool>(),
    access == FileAccess::TRANSIENT, access == FileAccess::TRANSIENT)
{
}

CompressedBlockReader::CompressedBlockReader(const std::string& filename,
  const std::shared_ptr<CompressedBlockPrefetchPool>& pool) :
  CompressedBlockReader(filename, pool, static_cast<bool>(pool))
{
}

CompressedBlockReader::CompressedBlockReader(const std::string& filename,
  const std::shared_ptr<CompressedBlockPrefetchPool>& pool,
  bool defer_first_block, bool transient) :
  input(filename, std::ios::binary), source_name(filename),
  transient_descriptor(transient), current(), current_block(0),
  current_offset(0), index_offset(0), block_count(0), block_size(0),
  current_logical(0), logical_bytes(0), physical_bytes(0), record_count(0),
  whole_checksum(0), sequential_bytes(0), sequential_checksum(FNV_OFFSET),
  sequential(true), cache_descriptor(-1), cache_released(0), retain_cache(false),
  prefetch_pool(pool), prefetch_owner(0), extent_cached(false),
  cached_extent_block(0), cached_physical(0), cached_logical(0),
  cached_next_physical(0), cached_next_logical(0)
{
  if(!this->input)
  {
    throw std::runtime_error("compressed block: open failure");
  }
  if(readLittle<std::uint64_t>(this->input, "header") != FILE_MAGIC ||
     readLittle<std::uint32_t>(this->input, "header") != FORMAT_VERSION ||
     readLittle<std::uint32_t>(this->input, "header") != HEADER_SIZE)
  {
    throw std::runtime_error("compressed block: invalid header");
  }
  this->block_size = readLittle<std::uint64_t>(this->input, "header");
  if(this->block_size == 0)
  {
    throw std::runtime_error("compressed block: invalid block size");
  }
  std::uint32_t mode = readLittle<std::uint32_t>(this->input, "header");
  if(mode != CompressedBlockWriter::RAW && mode != CompressedBlockWriter::ZSTD)
  {
    throw std::runtime_error("compressed block: invalid requested codec");
  }
  // Format-v2 artifacts written by early builds stored zero in this reserved
  // field. New writers store the operational worker count, but readers do not
  // need it and deliberately retain compatibility with both encodings.
  readLittle<std::uint32_t>(this->input, "header");

  this->input.seekg(0, std::ios::end);
  std::streamoff end = this->input.tellg();
  if(end < 0 || end < static_cast<std::streamoff>(HEADER_SIZE + FOOTER_SIZE))
  {
    throw std::runtime_error("compressed block: truncated footer");
  }
  this->physical_bytes = static_cast<std::uint64_t>(end);
  this->input.seekg(end - FOOTER_SIZE);
  if(readLittle<std::uint64_t>(this->input, "footer") != FOOTER_MAGIC ||
     readLittle<std::uint32_t>(this->input, "footer") != FORMAT_VERSION ||
     readLittle<std::uint32_t>(this->input, "footer") != FOOTER_SIZE)
  {
    throw std::runtime_error("compressed block: invalid footer");
  }
  this->block_count = readLittle<std::uint64_t>(this->input, "footer");
  this->record_count = readLittle<std::uint64_t>(this->input, "footer");
  this->logical_bytes = readLittle<std::uint64_t>(this->input, "footer");
  this->whole_checksum = readLittle<std::uint64_t>(this->input, "footer");
  this->index_offset = readLittle<std::uint64_t>(this->input, "footer");

  if(this->block_count > std::numeric_limits<std::size_t>::max() ||
     this->block_count > std::numeric_limits<std::uint64_t>::max() /
       INDEX_ENTRY_SIZE)
  {
    throw std::runtime_error("compressed block: invalid offset index size");
  }
  std::uint64_t index_bytes = this->block_count * INDEX_ENTRY_SIZE;
  if(this->index_offset < HEADER_SIZE ||
     this->index_offset > std::numeric_limits<std::uint64_t>::max() - index_bytes ||
     this->index_offset + index_bytes !=
       static_cast<std::uint64_t>(end - FOOTER_SIZE))
  {
    throw std::runtime_error("compressed block: invalid offset index");
  }

  if(this->block_count == 0)
  {
    if(this->logical_bytes != 0 || this->record_count != 0 ||
       this->whole_checksum != FNV_OFFSET)
    {
      throw std::runtime_error("compressed block: invalid empty stream");
    }
  }
  else
  {
    std::uint64_t physical, logical;
    this->indexEntry(0, physical, logical);
    if(physical != HEADER_SIZE || logical != 0)
    {
      throw std::runtime_error("compressed block: invalid first index");
    }
    this->indexEntry(this->block_count - 1, physical, logical);
    if(physical >= this->index_offset || logical >= this->logical_bytes)
    {
      throw std::runtime_error("compressed block: invalid last index");
    }
  }

  // Do not retain a second descriptor merely to support positional reads and
  // cache advice. Close the metadata stream first, then use one raw descriptor
  // for every operation after construction.
  this->input.close();
  if(!this->transient_descriptor)
  {
    this->cache_descriptor = ::open(filename.c_str(), O_RDONLY);
    if(this->cache_descriptor < 0)
    {
      throw std::runtime_error("compressed block: descriptor open failure");
    }
#if defined(POSIX_FADV_SEQUENTIAL)
    static_cast<void>(::posix_fadvise(this->cache_descriptor, 0, 0,
      POSIX_FADV_SEQUENTIAL));
#endif
  }

  try
  {
    if(this->block_count > 0)
    {
      if(defer_first_block && this->prefetch_pool)
      {
        this->prefetch_owner = this->prefetch_pool->registerReader();
        this->prefetchBlock(0);
      }
      else if(!defer_first_block)
      {
        this->loadBlock(0);
      }
    }
  }
  catch(...)
  {
    if(this->prefetch_pool && this->prefetch_owner != 0)
    {
      this->prefetch_pool->cancel(this->prefetch_owner);
      this->prefetch_owner = 0;
    }
    ::close(this->cache_descriptor); this->cache_descriptor = -1;
    throw;
  }
}

CompressedBlockReader::DescriptorLease::DescriptorLease(
  CompressedBlockReader& owner) : reader(owner), acquired(false)
{
  // Metadata validation still uses the constructor's ifstream. Nested leases
  // reuse the outer descriptor and never close it prematurely.
  if(this->reader.cache_descriptor < 0 && !this->reader.input.is_open())
  {
    this->reader.cache_descriptor = ::open(this->reader.source_name.c_str(),
      O_RDONLY);
    if(this->reader.cache_descriptor < 0)
    {
      throw std::runtime_error("compressed block: descriptor reopen failure");
    }
    this->acquired = this->reader.transient_descriptor;
  }
}

CompressedBlockReader::DescriptorLease::~DescriptorLease()
{
  if(this->acquired)
  {
    ::close(this->reader.cache_descriptor);
    this->reader.cache_descriptor = -1;
  }
}

void
CompressedBlockReader::prefetchBlock(std::size_t block)
{
  if(!(this->prefetch_pool) || this->prefetch_owner == 0 ||
     block >= this->block_count)
  {
    return;
  }
  std::uint64_t physical, logical;
  std::uint64_t next_physical, next_logical;
  this->blockExtent(block, physical, logical, next_physical, next_logical);
  this->prefetch_pool->submit(this->prefetch_owner, block,
    this->cache_descriptor, this->block_size, physical, logical,
    next_physical, next_logical);
}

void
CompressedBlockReader::resetPrefetch()
{
  if(!(this->prefetch_pool)) { return; }
  if(this->prefetch_owner != 0)
  {
    this->prefetch_pool->cancel(this->prefetch_owner);
  }
  this->extent_cached = false;
  this->prefetch_owner = this->prefetch_pool->registerReader();
}

void
CompressedBlockReader::blockExtent(std::size_t block,
  std::uint64_t& physical, std::uint64_t& logical,
  std::uint64_t& next_physical, std::uint64_t& next_logical)
{
  if(block >= this->block_count)
  {
    throw std::out_of_range("compressed block extent past end");
  }
  if(this->extent_cached && this->cached_extent_block == block)
  {
    physical = this->cached_physical; logical = this->cached_logical;
    next_physical = this->cached_next_physical;
    next_logical = this->cached_next_logical;
    return;
  }

  // Fetch enough consecutive entries to validate both boundaries and their
  // predecessor relation in one syscall. The hot sequential path caches this
  // extent from prefetch submission until the block is consumed.
  const std::uint64_t first = (block == 0 ? block : block - 1);
  const std::uint64_t last = std::min<std::uint64_t>(
    this->block_count - 1, block + 1);
  const std::size_t entries = static_cast<std::size_t>(last - first + 1);
  std::uint8_t encoded[3 * INDEX_ENTRY_SIZE];
  DescriptorLease descriptor(*this);
  preadExactCounted(this->cache_descriptor, encoded,
    entries * INDEX_ENTRY_SIZE, this->index_offset + first * INDEX_ENTRY_SIZE,
    "offset index extent");

  std::uint64_t entry_physical[3], entry_logical[3];
  for(std::size_t i = 0; i < entries; ++i)
  {
    const std::uint8_t* entry = encoded + i * INDEX_ENTRY_SIZE;
    entry_physical[i] = readLittleMemory<std::uint64_t>(entry);
    entry_logical[i] = readLittleMemory<std::uint64_t>(
      entry + sizeof(std::uint64_t));
    if(entry_physical[i] < HEADER_SIZE ||
       entry_physical[i] >= this->index_offset ||
       entry_logical[i] > this->logical_bytes ||
       (i > 0 && (entry_physical[i] <= entry_physical[i - 1] ||
         entry_logical[i] <= entry_logical[i - 1])))
    {
      throw std::runtime_error("compressed block: invalid offset index extent");
    }
  }

  const std::size_t current = static_cast<std::size_t>(block - first);
  physical = entry_physical[current]; logical = entry_logical[current];
  if(block + 1 < this->block_count)
  {
    next_physical = entry_physical[current + 1];
    next_logical = entry_logical[current + 1];
  }
  else
  {
    next_physical = this->index_offset; next_logical = this->logical_bytes;
  }
  this->extent_cached = true; this->cached_extent_block = block;
  this->cached_physical = physical; this->cached_logical = logical;
  this->cached_next_physical = next_physical;
  this->cached_next_logical = next_logical;
}

void
CompressedBlockReader::indexEntry(std::uint64_t block,
  std::uint64_t& physical, std::uint64_t& logical)
{
  if(block >= this->block_count)
  {
    throw std::out_of_range("compressed block index past end");
  }
  const std::uint64_t entry_offset =
    this->index_offset + block * INDEX_ENTRY_SIZE;
  DescriptorLease descriptor(*this);
  if(this->cache_descriptor >= 0)
  {
    std::uint8_t entry[INDEX_ENTRY_SIZE];
    preadExactCounted(this->cache_descriptor, entry, sizeof(entry),
      entry_offset, "offset index");
    physical = readLittleMemory<std::uint64_t>(entry);
    logical = readLittleMemory<std::uint64_t>(entry + sizeof(std::uint64_t));
  }
  else
  {
    this->input.clear();
    this->input.seekg(entry_offset);
    physical = readLittle<std::uint64_t>(this->input, "offset index");
    logical = readLittle<std::uint64_t>(this->input, "offset index");
  }
  if(physical < HEADER_SIZE || physical >= this->index_offset ||
     logical > this->logical_bytes)
  {
    throw std::runtime_error("compressed block: invalid index entry");
  }
  if(block > 0)
  {
    std::uint64_t previous_offset =
      this->index_offset + (block - 1) * INDEX_ENTRY_SIZE;
    std::uint64_t previous_physical, previous_logical;
    if(this->cache_descriptor >= 0)
    {
      std::uint8_t entry[INDEX_ENTRY_SIZE];
      preadExactCounted(this->cache_descriptor, entry, sizeof(entry),
        previous_offset, "offset index");
      previous_physical = readLittleMemory<std::uint64_t>(entry);
      previous_logical = readLittleMemory<std::uint64_t>(
        entry + sizeof(std::uint64_t));
    }
    else
    {
      this->input.seekg(previous_offset);
      previous_physical = readLittle<std::uint64_t>(this->input, "offset index");
      previous_logical = readLittle<std::uint64_t>(this->input, "offset index");
    }
    if(physical <= previous_physical || logical <= previous_logical)
    {
      throw std::runtime_error("compressed block: unordered index");
    }
  }
}

CompressedBlockReader::~CompressedBlockReader()
{
  if(this->prefetch_pool && this->prefetch_owner != 0)
  {
    this->prefetch_pool->cancel(this->prefetch_owner);
    this->prefetch_owner = 0;
  }
  if(this->cache_descriptor >= 0) { ::close(this->cache_descriptor); }
}

/*
  Release the pages of blocks already consumed, keeping a tail so a reader that
  steps back within the current window still hits cache. Only a stream that has
  read strictly forward is eligible: `sequential` is cleared by any backward or
  random seek, and the final event scan does seek backward.
*/
void
CompressedBlockReader::releaseConsumedCache(std::uint64_t physical_offset)
{
#if defined(POSIX_FADV_DONTNEED)
  if(this->cache_descriptor < 0 || !this->sequential || this->retain_cache) { return; }
  if(physical_offset < this->cache_released + BLOCK_CACHE_FLUSH_BYTES) { return; }
  std::uint64_t discard_end = physical_offset - BLOCK_CACHE_TAIL_BYTES;
  if(discard_end <= this->cache_released) { return; }
  static_cast<void>(::posix_fadvise(this->cache_descriptor,
    static_cast<off_t>(this->cache_released),
    static_cast<off_t>(discard_end - this->cache_released), POSIX_FADV_DONTNEED));
  this->cache_released = discard_end;
#else
  static_cast<void>(physical_offset);
#endif
}

void
CompressedBlockReader::loadBlock(std::size_t block)
{
  if(block >= this->block_count)
  {
    this->current.clear();
    this->current_block = block;
    this->current_offset = 0;
    this->current_logical = this->logical_bytes;
    return;
  }

  DescriptorLease descriptor(*this);
  std::uint64_t physical, logical;
  std::uint64_t next_physical, next_logical;
  this->blockExtent(block, physical, logical, next_physical, next_logical);

  if(this->prefetch_pool && this->prefetch_owner != 0)
  {
    // The old block has been consumed. Release it before taking the new block,
    // then transfer the cached allocation directly into the reader while the
    // pool lock is held. When take() releases pool capacity, the allocation is
    // already covered by this reader's separately reserved current-block
    // workspace; there is no unaccounted pool-to-reader handoff interval.
    std::vector<std::uint8_t>().swap(this->current);
    if(this->prefetch_pool->take(this->prefetch_owner, block, this->current))
    {
      if(this->current.empty() || this->current.size() != next_logical - logical)
      {
        throw std::runtime_error("compressed block: invalid prefetched block");
      }
      this->current_block = block;
      this->current_offset = 0;
      this->current_logical = logical;
      this->extent_cached = false;
      this->releaseConsumedCache(next_physical);
      if(this->sequential && block + 1 < this->block_count)
      {
        this->prefetchBlock(block + 1);
      }
      return;
    }
  }

  std::uint8_t header[BLOCK_HEADER_SIZE];
  preadExactCounted(this->cache_descriptor, header, sizeof(header), physical,
    "block header");
  if(readLittleMemory<std::uint32_t>(header) != BLOCK_MAGIC)
  {
    throw std::runtime_error("compressed block: invalid block header");
  }
  std::uint32_t codec = readLittleMemory<std::uint32_t>(header + 4);
  std::uint64_t records = readLittleMemory<std::uint64_t>(header + 8);
  std::uint64_t raw_bytes = readLittleMemory<std::uint64_t>(header + 16);
  std::uint64_t stored_bytes = readLittleMemory<std::uint64_t>(header + 24);
  std::uint64_t expected_checksum =
    readLittleMemory<std::uint64_t>(header + 32);
  bool valid_codec_extent =
    (codec == RAW_BLOCK && stored_bytes == raw_bytes) ||
    (codec == ZSTD_BLOCK && stored_bytes < raw_bytes);
  if(!valid_codec_extent || records == 0 || raw_bytes == 0 ||
     raw_bytes > this->block_size ||
     records > raw_bytes ||
     raw_bytes > std::numeric_limits<std::size_t>::max() ||
     stored_bytes > std::numeric_limits<std::size_t>::max() ||
     physical > next_physical || next_physical - physical < BLOCK_HEADER_SIZE ||
     stored_bytes != next_physical - physical - BLOCK_HEADER_SIZE ||
     logical > next_logical || raw_bytes != next_logical - logical)
  {
    throw std::runtime_error("compressed block: invalid block extent");
  }

  std::vector<std::uint8_t> encoded(static_cast<std::size_t>(stored_bytes));
  preadExactCounted(this->cache_descriptor, encoded.data(), encoded.size(),
    physical + BLOCK_HEADER_SIZE, "block payload");
  this->current.resize(static_cast<std::size_t>(raw_bytes));
  if(codec == RAW_BLOCK)
  {
    this->current.swap(encoded);
  }
  else if(codec == ZSTD_BLOCK)
  {
    std::size_t result = ZSTD_decompress(this->current.data(),
      this->current.size(), encoded.data(), encoded.size());
    if(ZSTD_isError(result) || result != this->current.size())
    {
      throw std::runtime_error("compressed block: zstd decompression failed");
    }
  }
  else
  {
    throw std::runtime_error("compressed block: invalid block codec");
  }
  if(checksum(this->current.data(), this->current.size()) != expected_checksum)
  {
    throw std::runtime_error("compressed block: block checksum failure");
  }
  if(this->prefetch_pool)
  {
    this->prefetch_pool->recordSynchronous(
      next_physical - physical, raw_bytes);
  }
  this->current_block = block;
  this->current_offset = 0;
  this->current_logical = logical;
  this->extent_cached = false;
  // Everything before the next block's offset has been consumed by a strictly
  // forward reader.
  this->releaseConsumedCache(next_physical);
  if(this->prefetch_pool && this->sequential &&
     block + 1 < this->block_count)
  {
    this->prefetchBlock(block + 1);
  }
}

void
CompressedBlockReader::checkSequentialChecksum()
{
  if(VERIFY_WHOLE_STREAM_CHECKSUM &&
     this->sequential && this->sequential_bytes == this->logical_bytes &&
     this->sequential_checksum != this->whole_checksum)
  {
    throw std::runtime_error("compressed block: stream checksum failure");
  }
}

std::size_t
CompressedBlockReader::read(void* data, std::size_t bytes)
{
  if(bytes > 0 && data == nullptr)
  {
    throw std::invalid_argument("null compressed block read buffer");
  }
  std::uint8_t* output = static_cast<std::uint8_t*>(data);
  std::size_t copied = 0;
  if(bytes > 0 && this->current.empty() &&
     this->current_block < this->block_count)
  {
    this->loadBlock(this->current_block);
  }
  while(copied < bytes && this->current_block < this->block_count)
  {
    if(this->current_offset == this->current.size())
    {
      this->loadBlock(this->current_block + 1);
      continue;
    }
    std::size_t chunk = std::min(bytes - copied,
      this->current.size() - this->current_offset);
    std::memcpy(output + copied,
      this->current.data() + this->current_offset, chunk);
    if(VERIFY_WHOLE_STREAM_CHECKSUM && this->sequential)
    {
      this->sequential_checksum = checksum(
        this->current.data() + this->current_offset, chunk,
        this->sequential_checksum);
      this->sequential_bytes += chunk;
    }
    this->current_offset += chunk;
    copied += chunk;
  }
  this->checkSequentialChecksum();
  return copied;
}

std::size_t
CompressedBlockReader::readAt(std::uint64_t offset, void* data,
  std::size_t bytes)
{
  std::uint64_t cursor = (this->current_block < this->block_count ?
    this->current_logical + this->current_offset : this->logical_bytes);
  if(offset == cursor) { return this->read(data, bytes); }
  this->seekUncompressedByte(offset);
  return this->read(data, bytes);
}

void
CompressedBlockReader::seekBlock(std::uint64_t block)
{
  if(block > this->block_count)
  {
    throw std::out_of_range("compressed block seek past end");
  }
  this->resetPrefetch();
  this->sequential = (block == 0);
  this->sequential_bytes = 0;
  this->sequential_checksum = FNV_OFFSET;
  this->loadBlock(block);
}

void
CompressedBlockReader::seekUncompressedByte(std::uint64_t offset)
{
  if(offset > this->logical_bytes)
  {
    throw std::out_of_range("compressed block seek past end");
  }
  std::uint64_t cursor = (this->current_block < this->block_count ?
    this->current_logical + this->current_offset : this->logical_bytes);
  if(offset == cursor) { return; }
  if(offset == this->logical_bytes)
  {
    this->resetPrefetch();
    this->sequential = false;
    this->loadBlock(this->block_count);
    return;
  }
  if(this->current_block < this->block_count &&
     offset >= this->current_logical &&
     offset < this->current_logical + this->current.size())
  {
    if(this->sequential) { this->resetPrefetch(); }
    this->sequential = this->sequential && (offset == this->sequential_bytes);
    this->current_offset = offset - this->current_logical;
    return;
  }

  DescriptorLease descriptor(*this);
  std::uint64_t low = 0, high = this->block_count;
  while(low + 1 < high)
  {
    std::uint64_t middle = low + (high - low) / 2;
    std::uint64_t physical, logical;
    this->indexEntry(middle, physical, logical);
    if(logical <= offset) { low = middle; }
    else { high = middle; }
  }
  std::uint64_t physical, logical;
  this->indexEntry(low, physical, logical);
  this->resetPrefetch();
  this->sequential = (offset == 0);
  this->sequential_bytes = 0;
  this->sequential_checksum = FNV_OFFSET;
  this->loadBlock(low);
  this->current_offset = offset - logical;
}

} // namespace gcsa
