#include <gcsa/compressed_block.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <limits.h>
#include <sys/resource.h>
#include <unistd.h>

using namespace gcsa;

namespace
{

void
require(bool condition)
{
  if(!condition) { std::abort(); }
}

std::vector<std::uint8_t>
readAll(const std::string& filename)
{
  std::ifstream input(filename.c_str(), std::ios::binary);
  return std::vector<std::uint8_t>(
    std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>());
}

bool
rejected(const std::string& filename)
{
  try
  {
    CompressedBlockReader input(filename);
    std::uint8_t byte;
    while(input.read(&byte, 1)) { }
  }
  catch(const std::runtime_error&)
  {
    return true;
  }
  return false;
}

// OutputStream writes the same bytes as std::ofstream in whichever mode the
// environment selects: around the direct writer's 4 KiB alignment and 4 MiB
// block boundaries, through a flush in mid-stream, with tellp() positions.
void
outputStreamTest(const std::string& root)
{
  const std::size_t block = std::size_t(4) << 20;
  const std::size_t sizes[] = { 0, 1, 4095, 4096, 4097, block - 1, block, block + 1, 10 * (std::size_t(1) << 20) + 123 };
  std::vector<std::uint8_t> data(10 * (std::size_t(1) << 20) + 200);
  for(std::size_t i = 0; i < data.size(); i++) { data[i] = static_cast<std::uint8_t>((i * 2654435761u) >> 13); }
  const std::string name = root + "/output-stream";
  for(std::size_t size : sizes)
  {
    for(bool flush_midway : { false, true })
    {
      const std::size_t half = size / 2;
      {
        OutputStream output(name, std::ios::binary | std::ios::trunc);
        require(output.is_open() && output.good());
        output.write(reinterpret_cast<const char*>(data.data()), half);
        if(flush_midway)
        {
          output.flush();
          require(output.good());
          require(readAll(name) == std::vector<std::uint8_t>(data.begin(), data.begin() + half));
        }
        require(static_cast<std::size_t>(output.tellp()) == half);
        output.write(reinterpret_cast<const char*>(data.data() + half), size - half);
        require(static_cast<std::size_t>(output.tellp()) == size);
        output.close();
        require(!output.fail() && !output.is_open());
      }
      require(readAll(name) == std::vector<std::uint8_t>(data.begin(), data.begin() + size));
    }
  }
  require(::unlink(name.c_str()) == 0);
}

} // namespace

int
main()
{
  const char* configured_tmp = std::getenv("TMPDIR");
  std::string pattern = std::string(configured_tmp != nullptr ? configured_tmp : "/tmp") +
    "/gcsa-compressed-block-XXXXXX";
  std::vector<char> root_buffer(pattern.begin(), pattern.end());
  root_buffer.push_back(0); char* root = root_buffer.data();
  require(mkdtemp(root) != nullptr);
  std::string first = std::string(root) + "/first";
  std::string second = std::string(root) + "/second";
  std::string parallel = std::string(root) + "/parallel";
  std::string parallel_copy = std::string(root) + "/parallel-copy";

  std::vector<std::uint8_t> data(83);
  for(std::size_t i = 0; i < data.size(); ++i)
  {
    data[i] = static_cast<std::uint8_t>(i * 37);
  }
  {
    CompressedBlockWriter output(first, 16, CompressedBlockWriter::ZSTD);
    for(std::size_t i = 0; i < data.size(); i += 7)
    {
      output.writeRecord(data.data() + i,
        std::min<std::size_t>(7, data.size() - i));
    }
    output.finish();
  }
  {
    CompressedBlockReader input(first);
    require(CompressedBlockReader::isFramed(first));
    require(CompressedBlockReader::declaredBlockSize(first) == 16);
    require(input.blocks() > 2);
    require(input.logicalSize() == data.size());
    require(input.physicalSize() > input.logicalSize());
    require(input.blockSize() == 16);
    require(input.indexMemoryBytes() == 0);
    require(CompressedBlockReader::workingMemoryEstimate(input.blockSize()) >=
      input.blockSize());
    std::vector<std::uint8_t> observed(data.size());
    require(input.read(observed.data(), observed.size()) == observed.size());
    require(observed == data);
    input.seekUncompressedByte(31);
    std::uint8_t byte;
    require(input.read(&byte, 1) == 1 && byte == data[31]);
    input.seekBlock(2);
    require(input.read(&byte, 1) == 1 && byte == data[28]);
    require(input.readAt(4, &byte, 1) == 1 && byte == data[4]);
  }
  // Decoded blocks can outlive their backing descriptors. More resident
  // readers than the process FD ceiling must still support random/sequential
  // reads, and a cache hit must not reopen the file.
  {
    struct rlimit original;
    require(getrlimit(RLIMIT_NOFILE, &original) == 0);
    struct rlimit limited = original;
    limited.rlim_cur = std::min<rlim_t>(original.rlim_cur, 64);
    require(setrlimit(RLIMIT_NOFILE, &limited) == 0);
    std::vector<std::unique_ptr<CompressedBlockReader>> readers;
    for(std::size_t i = 0; i < 128; i++)
    {
      readers.emplace_back(new CompressedBlockReader(first,
        CompressedBlockReader::FileAccess::TRANSIENT));
      std::uint8_t byte;
      require(readers.back()->readAt(2, &byte, 1) == 1 && byte == data[2]);
    }
    const std::string hidden = first + ".hidden";
    require(std::rename(first.c_str(), hidden.c_str()) == 0);
    for(auto& reader : readers)
    {
      std::uint8_t byte;
      require(reader->readAt(3, &byte, 1) == 1 && byte == data[3]);
      bool failed = false;
      try { reader->readAt(60, &byte, 1); }
      catch(const std::runtime_error&) { failed = true; }
      require(failed);
    }
    require(std::rename(hidden.c_str(), first.c_str()) == 0);
    for(auto& reader : readers)
    {
      std::uint8_t byte;
      require(reader->readAt(60, &byte, 1) == 1 && byte == data[60]);
      reader->seekBlock(0);
      std::vector<std::uint8_t> observed(data.size());
      require(reader->read(observed.data(), observed.size()) == observed.size());
      require(observed == data);
    }
    readers.clear();
    require(setrlimit(RLIMIT_NOFILE, &original) == 0);
  }
  // Multiple readers share fixed workers and a single byte bound. Block 0 is
  // queued by the pooled constructor, allowing a caller such as the path
  // merger to open all shard streams before it consumes their heads.
  {
    const std::size_t job =
      CompressedBlockReader::prefetchWorkingMemoryEstimate(16);
    std::shared_ptr<CompressedBlockPrefetchPool> pool(
      new CompressedBlockPrefetchPool(2 * job, 2));
    {
      CompressedBlockReader first_input(first, pool);
      CompressedBlockReader second_input(first, pool);
      for(std::size_t attempt = 0; attempt < 1000 &&
          pool->stats().completed < 2; ++attempt)
      {
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
      }
      require(pool->stats().completed >= 2);

      std::vector<std::uint8_t> first_observed(data.size());
      for(std::size_t offset = 0; offset < data.size(); offset += 7)
      {
        std::size_t bytes = std::min<std::size_t>(7, data.size() - offset);
        require(first_input.readAt(offset, first_observed.data() + offset,
          bytes) == bytes);
      }
      std::vector<std::uint8_t> second_observed(data.size());
      require(second_input.read(second_observed.data(), second_observed.size()) ==
        second_observed.size());
      require(first_observed == data); require(second_observed == data);

      // Loading block 0 queues block 1. A nonsequential seek cancels that one
      // speculative block and falls back to the ordinary seekable reader.
      CompressedBlockReader random_input(first, pool);
      std::uint8_t byte;
      require(random_input.read(&byte, 1) == 1 && byte == data[0]);
      require(random_input.readAt(50, &byte, 1) == 1 && byte == data[50]);

      CompressedBlockPrefetchPool::Stats stats = pool->stats();
      require(stats.worker_threads == 2);
      require(stats.submitted >= 3); require(stats.completed >= 2);
      require(stats.consumed >= 2); require(stats.ready_hits >= 2);
      require(stats.cancelled >= 1); require(stats.synchronous_blocks >= 1);
      require(stats.physical_bytes > 0); require(stats.decoded_bytes > 0);
      require(stats.peak_bytes <= pool->memoryLimit());
      require(pool->usedBytes() <= pool->memoryLimit());
    }
    require(pool->usedBytes() == 0);
  }
  require(readAll(first).size() <=
    CompressedBlockWriter::maximumTemporaryBytes(data.size(), 7, 16));
  // Even an empty framed stream has a fixed header/footer. The peak bound also
  // models the temporary index copy, which is empty in this case.
  require(CompressedBlockWriter::maximumTemporaryBytes(0, 7, 16) == 88);
  {
    CompressedBlockWriter output(second, 16, CompressedBlockWriter::ZSTD);
    for(std::size_t i = 0; i < data.size(); i += 7)
    {
      output.writeRecord(data.data() + i,
        std::min<std::size_t>(7, data.size() - i));
    }
    output.finish();
  }
  require(readAll(first) == readAll(second));

  {
    CompressedBlockWriter output(second, 128, CompressedBlockWriter::RAW);
    output.writeRecord(data.data(), data.size());
    output.finish();
    CompressedBlockReader input(second);
    std::vector<std::uint8_t> observed(data.size());
    require(input.read(observed.data(), observed.size()) == observed.size());
    require(observed == data);
  }
  {
    std::vector<std::uint8_t> zeros(64, 0);
    CompressedBlockWriter output(second, 64, CompressedBlockWriter::ZSTD);
    output.writeRecord(zeros.data(), zeros.size());
    output.finish();
    std::vector<std::uint8_t> encoded = readAll(second);
    require(encoded.size() > 72 && encoded[32 + 4] == 1);
  }
  {
    CompressedBlockWriter output(second, 8, CompressedBlockWriter::RAW);
    output.finish();
    CompressedBlockReader input(second);
    require(input.blocks() == 0 && input.bytes() == 0);
  }

  // Internal zstd workers are explicitly represented in the scheduler's
  // reservation and produce a normal, seekable framed stream.
  {
    std::vector<std::uint8_t> zeros(1024 * 1024, 0);
    CompressedBlockWriter output(parallel, 128 * 1024,
      CompressedBlockWriter::ZSTD, 1, 2);
    for(std::size_t offset = 0; offset < zeros.size(); offset += 16 * 1024)
    {
      output.writeRecord(zeros.data() + offset, 16 * 1024);
    }
    output.finish();
    CompressedBlockReader input(parallel);
    std::vector<std::uint8_t> observed(zeros.size());
    require(input.read(observed.data(), observed.size()) == observed.size());
    require(observed == zeros);
    // zstd starts no worker before a whole job exists, and a job is never
    // smaller than one mebibyte, so a 128 KiB block runs single-threaded
    // whatever is requested. Charging the request would reserve a job pool the
    // encoder never allocates.
    require(CompressedBlockWriter::workingMemoryEstimate(128 * 1024,
      CompressedBlockWriter::ZSTD, 1, 2) ==
      CompressedBlockWriter::workingMemoryEstimate(128 * 1024,
        CompressedBlockWriter::ZSTD, 1, 1));
    // A block that does hold several jobs engages them, and is charged for the
    // per-worker contexts, round buffer and job-output pool they allocate.
    require(CompressedBlockWriter::workingMemoryEstimate(16 * 1024 * 1024,
      CompressedBlockWriter::ZSTD, 1, 4) >
      CompressedBlockWriter::workingMemoryEstimate(16 * 1024 * 1024,
        CompressedBlockWriter::ZSTD, 1, 1));
    CompressedBlockWriter copy(parallel_copy, 128 * 1024,
      CompressedBlockWriter::ZSTD, 1, 2);
    for(std::size_t offset = 0; offset < zeros.size(); offset += 16 * 1024)
    {
      copy.writeRecord(zeros.data() + offset, 16 * 1024);
    }
    copy.finish();
    require(readAll(parallel) == readAll(parallel_copy));
  }

  // Atomic publication must sync the containing directory even when the
  // caller supplies a relative output name.
  {
    char previous[PATH_MAX];
    require(getcwd(previous, sizeof(previous)) != nullptr);
    require(chdir(root) == 0);
    CompressedBlockWriter output("relative", 128, CompressedBlockWriter::RAW);
    output.writeRecord(data.data(), data.size());
    output.finish();
    require(CompressedBlockReader::isFramed("relative"));
    require(::unlink("relative") == 0);
    require(chdir(previous) == 0);
  }

  std::vector<std::uint8_t> corrupt = readAll(first);
  corrupt.resize(corrupt.size() - 1);
  {
    std::ofstream output(second.c_str(), std::ios::binary | std::ios::trunc);
    output.write(reinterpret_cast<const char*>(corrupt.data()), corrupt.size());
  }
  require(rejected(second));

  corrupt = readAll(first);
  corrupt[72] ^= 1;
  {
    std::ofstream output(second.c_str(), std::ios::binary | std::ios::trunc);
    output.write(reinterpret_cast<const char*>(corrupt.data()), corrupt.size());
  }
  require(rejected(second));
  // The deferred constructor must surface corruption found by a worker, retain
  // no failed cache allocation, and close its sole descriptor during unwind.
  {
    const std::size_t job =
      CompressedBlockReader::prefetchWorkingMemoryEstimate(16);
    std::shared_ptr<CompressedBlockPrefetchPool> pool(
      new CompressedBlockPrefetchPool(job, 1));
    bool did_reject = false;
    try
    {
      CompressedBlockReader input(second, pool);
      for(std::size_t attempt = 0; attempt < 1000 &&
          pool->stats().errors == 0; ++attempt)
      {
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
      }
      std::uint8_t byte;
      static_cast<void>(input.read(&byte, 1));
    }
    catch(const std::runtime_error&)
    {
      did_reject = true;
    }
    require(did_reject); require(pool->stats().errors == 1);
    require(pool->usedBytes() == 0);
  }

  // The writer only labels a block as zstd when it is smaller than the raw
  // payload. Reject an inconsistent extent before allocating the payload.
  corrupt = readAll(first);
  require(corrupt[36] == 0);
  corrupt[36] = 1;
  {
    std::ofstream output(second.c_str(), std::ios::binary | std::ios::trunc);
    output.write(reinterpret_cast<const char*>(corrupt.data()), corrupt.size());
  }
  require(rejected(second));

  outputStreamTest(root);

  require(::unlink(first.c_str()) == 0);
  require(::unlink(second.c_str()) == 0);
  require(::unlink(parallel.c_str()) == 0);
  require(::unlink(parallel_copy.c_str()) == 0);
  require(::rmdir(root) == 0);
  return 0;
}
