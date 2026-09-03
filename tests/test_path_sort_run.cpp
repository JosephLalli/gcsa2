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

void drain(const std::string& filename)
{
  PathSortRunReader reader(filename, pathSortRunMinimumBuffer());
  while(!(reader.atEnd())) { reader.advance(); }
}

void copyAndMutate(const std::string& source, const std::string& target,
  off_t offset, std::uint8_t mask)
{
  std::ifstream input(source.c_str(), std::ios_base::binary);
  std::ofstream output(target.c_str(), std::ios_base::binary);
  output << input.rdbuf(); output.close();
  int descriptor = ::open(target.c_str(), O_RDWR); require(descriptor >= 0);
  std::uint8_t value = 0;
  require(::pread(descriptor, &value, 1, offset) == 1);
  value ^= mask;
  require(::pwrite(descriptor, &value, 1, offset) == 1);
  require(::close(descriptor) == 0);
}

} // namespace

int main()
{
  char root[] = "/tmp/gcsa-path-sort-codec-XXXXXX";
  require(mkdtemp(root) != nullptr);
  TempFile::setDirectory(root);
  const std::string run = std::string(root) + "/labels.run";
  const std::string corrupt = std::string(root) + "/corrupt.run";
  const std::string bad_flags = std::string(root) + "/bad-flags.run";
  const std::string bad_reference = std::string(root) + "/bad-reference.run";
  const std::string bad_prefix = std::string(root) + "/bad-prefix.run";
  const std::string bad_reserved = std::string(root) + "/bad-reserved.run";
  const std::string truncated = std::string(root) + "/truncated.run";
  const off_t header_bytes = 40, footer_bytes = 40;

  std::vector<PathSortRunRecord> records(1000);
  for(size_type i = 0; i < records.size(); i++)
  {
    PathSortRunRecord& record = records[i];
    size_type group = i / 100;
    // Each group models one left path expanding to many right paths. The run
    // stores the immutable left context once and uses compact references for
    // the remaining 99 records, while labels remain independently decodable.
    record.node.from = group + 1; record.node.to = i + 2; record.node.fields = 0;
    record.node.setPredecessors(static_cast<byte_type>(1U << (group % 4)));
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
  require(writer.groupHeaders() == 10);
  require(writer.groupReferences() == records.size() - writer.groupHeaders());
  require(writer.contextBytesSaved() ==
    writer.groupHeaders() * 4 + writer.groupReferences() * 14);
  require(writer.bytes() < writer.uncompressedBytes() / 2);
  consume(run, records);

  // Corruption is detected by the payload checksum even when the modified
  // from-node byte remains structurally readable.
  copyAndMutate(run, corrupt, header_bytes + 1, 1);
  bool rejected = false;
  try { drain(corrupt); }
  catch(const std::runtime_error&) { rejected = true; }
  require(rejected);

  // Structural validation fails before a corrupt record can be admitted to a
  // merge. These offsets are stable parts of the explicit v2 wire format.
  copyAndMutate(run, bad_flags, header_bytes, 0x80);
  rejected = false;
  try { drain(bad_flags); }
  catch(const std::runtime_error&) { rejected = true; }
  require(rejected);

  copyAndMutate(run, bad_reference, header_bytes, 0x01);
  rejected = false;
  try { drain(bad_reference); }
  catch(const std::runtime_error&) { rejected = true; }
  require(rejected);

  // flags + from + predecessors + order + to + lcp precede label prefix.
  copyAndMutate(run, bad_prefix, header_bytes + 1 + 8 + 1 + 1 + 8 + 1, 0xFF);
  rejected = false;
  try { drain(bad_prefix); }
  catch(const std::runtime_error&) { rejected = true; }
  require(rejected);

  // The footer's reserved field must stay zero for forward-compatible format
  // negotiation; silently accepting it would make corruption ambiguous.
  copyAndMutate(run, bad_reserved,
    static_cast<off_t>(bytes(run)) - footer_bytes + 12, 0x01);
  rejected = false;
  try { drain(bad_reserved); }
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
  std::remove(bad_flags.c_str()); std::remove(bad_reference.c_str());
  std::remove(bad_prefix.c_str()); std::remove(bad_reserved.c_str());
  std::remove(truncated.c_str()); rmdir(root);
  return 0;
}
