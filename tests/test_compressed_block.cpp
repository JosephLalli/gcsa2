#include <compressed_block.hpp>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <filesystem>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

#include <limits.h>
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
    require(input.logicalSize() == data.size());
    require(CompressedBlockReader::workingMemoryEstimate(16) >= 16);
    std::vector<std::uint8_t> observed(data.size());
    require(input.read(observed.data(), observed.size()) == observed.size());
    require(observed == data);
    input.seekUncompressedByte(31);
    std::uint8_t byte;
    require(input.read(&byte, 1) == 1 && byte == data[31]);
    require(input.readAt(4, &byte, 1) == 1 && byte == data[4]);
  }
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
    require(input.logicalSize() == 0);
  }

  // Abandoning a writer must remove both construction files and leave the
  // intended destination unpublished.
  const std::string abandoned = std::string(root) + "/abandoned";
  {
    CompressedBlockWriter output(abandoned, 16, CompressedBlockWriter::ZSTD);
    output.writeRecord(data.data(), 7);
  }
  require(!std::filesystem::exists(abandoned));
  for(const auto& entry : std::filesystem::directory_iterator(root))
  {
    require(entry.path().filename().string().find("abandoned.") != 0);
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
  require(CompressedBlockReader::isFramed(second));
  require(rejected(second));

  corrupt = readAll(first);
  corrupt[8]++;
  {
    std::ofstream output(second.c_str(), std::ios::binary | std::ios::trunc);
    output.write(reinterpret_cast<const char*>(corrupt.data()), corrupt.size());
  }
  require(CompressedBlockReader::isFramed(second));
  require(rejected(second));

  corrupt = readAll(first);
  corrupt[72] ^= 1;
  {
    std::ofstream output(second.c_str(), std::ios::binary | std::ios::trunc);
    output.write(reinterpret_cast<const char*>(corrupt.data()), corrupt.size());
  }
  require(rejected(second));
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

  require(::unlink(first.c_str()) == 0);
  require(::unlink(second.c_str()) == 0);
  require(::rmdir(root) == 0);
  return 0;
}
