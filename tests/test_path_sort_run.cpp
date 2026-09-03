#include <gcsa/path_sort_run.h>

#include <algorithm>
#include <cstdlib>
#include <fcntl.h>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

using namespace gcsa;

namespace
{

void require(bool value) { if(!value) { std::abort(); } }

size_type bytes(const std::string& filename)
{
  struct stat info; require(::stat(filename.c_str(), &info) == 0);
  require(info.st_size >= 0); return static_cast<size_type>(info.st_size);
}

bool equal(const PathSortRunRecord& left, const PathSortRunRecord& right)
{
  if(left.node.from != right.node.from || left.node.to != right.node.to ||
     left.node.predecessors() != right.node.predecessors() ||
     left.node.order() != right.node.order() || left.node.lcp() != right.node.lcp())
  {
    return false;
  }
  for(size_type i = 0; i < left.node.ranks(); i++)
  {
    if(left.labels[i] != right.labels[i]) { return false; }
  }
  return true;
}

void consume(const std::string& filename, const std::vector<PathSortRunRecord>& expected)
{
  PathSortRunReader reader(filename, pathSortRunMinimumBuffer());
  require(reader.records() == expected.size());
  size_type offset = 0;
  while(!(reader.atEnd()))
  {
    require(offset < expected.size());
    require(equal(reader.current(), expected[offset]));
    reader.advance(); offset++;
  }
  require(offset == expected.size());
}

} // namespace

int main()
{
  char root[] = "/tmp/gcsa-path-sort-codec-XXXXXX";
  require(mkdtemp(root) != nullptr);
  TempFile::setDirectory(root);
  const std::string run = std::string(root) + "/labels.run";
  const std::string corrupt = std::string(root) + "/corrupt.run";
  const std::string truncated = std::string(root) + "/truncated.run";

  std::vector<PathSortRunRecord> records(1000);
  for(size_type i = 0; i < records.size(); i++)
  {
    PathSortRunRecord& record = records[i];
    record.node.from = i + 1; record.node.to = i + 2; record.node.fields = 0;
    record.node.setPredecessors(static_cast<byte_type>(1U << (i % 4)));
    record.node.setOrder(16); record.node.setLCP(15);
    record.node.setPointer(i * 17); // The codec deliberately normalizes this.
    for(size_type rank = 0; rank < 17; rank++)
    {
      record.labels[rank] = (rank < 15 ? static_cast<PathNode::rank_type>(rank + 7) :
        static_cast<PathNode::rank_type>(i * 2 + rank));
    }
  }

  PathSortRunWriter writer(run, records.size(), pathSortRunMinimumBuffer());
  for(const PathSortRunRecord& record : records) { writer.write(record); }
  writer.finish();
  require(writer.records() == records.size());
  require(writer.bytes() == bytes(run));
  require(writer.bytes() < writer.uncompressedBytes() / 2);
  consume(run, records);

  // Corruption is detected by the payload checksum even when every encoded
  // field remains structurally readable.
  {
    std::ifstream input(run.c_str(), std::ios_base::binary);
    std::ofstream output(corrupt.c_str(), std::ios_base::binary);
    output << input.rdbuf(); output.close();
    int descriptor = ::open(corrupt.c_str(), O_RDWR); require(descriptor >= 0);
    std::uint8_t value = 0;
    off_t checksum_byte = static_cast<off_t>(bytes(corrupt) - 1);
    require(::pread(descriptor, &value, 1, checksum_byte) == 1);
    value ^= 1;
    require(::pwrite(descriptor, &value, 1, checksum_byte) == 1);
    require(::close(descriptor) == 0);
  }
  bool rejected = false;
  try { consume(corrupt, records); }
  catch(const std::runtime_error&) { rejected = true; }
  require(rejected);

  // A file whose footer is incomplete is rejected from header length alone.
  {
    std::ifstream input(run.c_str(), std::ios_base::binary);
    std::ofstream output(truncated.c_str(), std::ios_base::binary);
    std::vector<char> data(bytes(run) - 1);
    input.read(data.data(), data.size()); require(input.gcount() == static_cast<std::streamsize>(data.size()));
    output.write(data.data(), data.size()); output.close();
  }
  rejected = false;
  try { PathSortRunReader reader(truncated, pathSortRunMinimumBuffer()); }
  catch(const std::runtime_error&) { rejected = true; }
  require(rejected);

  std::remove(run.c_str()); std::remove(corrupt.c_str());
  std::remove(truncated.c_str()); rmdir(root);
  return 0;
}
