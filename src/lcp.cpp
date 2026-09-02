#include <gcsa/lcp.h>

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <limits>
#include <stack>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>

#include <gcsa/internal.h>
#include <gcsa/workspace.h>

namespace gcsa
{

namespace
{

const std::string LCP_STREAMING_TASK = "lcp-streaming";

std::string
levelPhase(size_type level)
{
  return "level-" + std::to_string(level);
}

ArtifactIdentity
levelArtifact(size_type level)
{
  return ArtifactIdentity(LCP_STREAMING_TASK, levelPhase(level),
    "range-minimum-level", "lcp-level-u8-v1");
}

[[noreturn]] void
systemError(const std::string& action, const std::string& filename)
{
  int error = errno;
  throw std::runtime_error(action + " " + filename + ": " + std::strerror(error));
}

size_type
rawFileSize(const std::string& filename)
{
  struct stat status;
  if(::stat(filename.c_str(), &status) != 0) { systemError("cannot stat", filename); }
  if(status.st_size < 0)
  {
    throw std::runtime_error("negative file size for " + filename);
  }
  if(static_cast<std::uintmax_t>(status.st_size) > std::numeric_limits<size_type>::max())
  {
    throw std::runtime_error("file is too large for LCP construction: " + filename);
  }
  return static_cast<size_type>(status.st_size);
}

int
openForRead(const std::string& filename)
{
  int fd = ::open(filename.c_str(), O_RDONLY);
  if(fd < 0) { systemError("cannot open", filename); }
  return fd;
}

int
openForWrite(const std::string& filename)
{
  int fd = ::open(filename.c_str(), O_WRONLY | O_CREAT | O_EXCL, 0644);
  if(fd < 0) { systemError("cannot create", filename); }
  return fd;
}

void
closeChecked(int fd, const std::string& filename)
{
  if(::close(fd) != 0) { systemError("cannot close", filename); }
}

void
readExact(int fd, std::uint8_t* buffer, size_type bytes, const std::string& filename)
{
  size_type offset = 0;
  while(offset < bytes)
  {
    size_type remaining = bytes - offset;
    size_type request = std::min(remaining,
      static_cast<size_type>(std::numeric_limits<ssize_t>::max()));
    ssize_t got = ::read(fd, buffer + offset, static_cast<size_t>(request));
    if(got < 0)
    {
      if(errno == EINTR) { continue; }
      systemError("cannot read", filename);
    }
    if(got == 0) { throw std::runtime_error("unexpected EOF in " + filename); }
    offset += static_cast<size_type>(got);
  }
}

void
requireEnd(int fd, const std::string& filename)
{
  std::uint8_t extra = 0;
  while(true)
  {
    ssize_t got = ::read(fd, &extra, 1);
    if(got < 0)
    {
      if(errno == EINTR) { continue; }
      systemError("cannot read", filename);
    }
    if(got > 0) { throw std::runtime_error("unexpected trailing data in " + filename); }
    return;
  }
}

void
writeExact(int fd, const std::uint8_t* buffer, size_type bytes, const std::string& filename)
{
  size_type offset = 0;
  while(offset < bytes)
  {
    size_type remaining = bytes - offset;
    size_type request = std::min(remaining,
      static_cast<size_type>(std::numeric_limits<ssize_t>::max()));
    ssize_t written = ::write(fd, buffer + offset, static_cast<size_t>(request));
    if(written < 0)
    {
      if(errno == EINTR) { continue; }
      systemError("cannot write", filename);
    }
    if(written == 0) { throw std::runtime_error("short write to " + filename); }
    offset += static_cast<size_type>(written);
  }
}

void
noteBuffers(LCPStreamingStats& stats, size_type bytes)
{
  stats.max_bytes_resident = std::max(stats.max_bytes_resident, bytes);
}

size_type
roundedDivision(size_type value, size_type divisor)
{
  return value / divisor + (value % divisor != 0);
}

std::vector<size_type>
levelSizes(size_type leaves, size_type branching)
{
  std::vector<size_type> result(1, leaves);
  while(result.back() > 1)
  {
    result.push_back(roundedDivision(result.back(), branching));
  }
  return result;
}

std::uint8_t
maximumValue(const std::string& filename, size_type values, size_type byte_budget,
  LCPStreamingStats& stats)
{
  std::vector<std::uint8_t> buffer(byte_budget);
  int input = openForRead(filename);
  noteBuffers(stats, buffer.size());
  std::uint8_t result = 0;
  try
  {
    size_type remaining = values;
    while(remaining > 0)
    {
      size_type bytes = std::min(remaining, static_cast<size_type>(buffer.size()));
      readExact(input, buffer.data(), bytes, filename);
      for(size_type i = 0; i < bytes; i++) { result = std::max(result, buffer[i]); }
      remaining -= bytes;
    }
    requireEnd(input, filename);
    closeChecked(input, filename); input = -1;
  }
  catch(...)
  {
    if(input >= 0) { ::close(input); }
    throw;
  }
  return result;
}

void
copyIntoData(LCPArray& target, size_type offset, const std::string& filename,
  size_type values, size_type byte_budget, LCPStreamingStats& stats)
{
  std::vector<std::uint8_t> buffer(byte_budget);
  int input = openForRead(filename);
  noteBuffers(stats, buffer.size());
  try
  {
    size_type copied = 0;
    while(copied < values)
    {
      size_type bytes = std::min(values - copied, static_cast<size_type>(buffer.size()));
      readExact(input, buffer.data(), bytes, filename);
      for(size_type i = 0; i < bytes; i++) { target.data[offset + copied + i] = buffer[i]; }
      copied += bytes;
    }
    requireEnd(input, filename);
    closeChecked(input, filename); input = -1;
  }
  catch(...)
  {
    if(input >= 0) { ::close(input); }
    throw;
  }
}

void
flushLevelOutput(int output, std::vector<std::uint8_t>& buffer, size_type& used,
  const std::string& filename, BuildWorkspace::ArtifactWriter* checkpoint)
{
  if(used == 0) { return; }
  writeExact(output, buffer.data(), used, filename);
  if(checkpoint != nullptr) { checkpoint->write(buffer.data(), used); }
  used = 0;
}

void
buildLevel(const std::string& source, size_type source_values,
  const std::string& target, size_type target_values, size_type branching,
  size_type byte_budget, BuildWorkspace::ArtifactWriter* checkpoint,
  LCPStreamingStats& stats)
{
  const size_type input_bytes = byte_budget / 2;
  const size_type output_bytes = byte_budget - input_bytes;
  std::vector<std::uint8_t> input_buffer(input_bytes), output_buffer(output_bytes);
  int input = openForRead(source), output = -1;
  noteBuffers(stats, input_buffer.size() + output_buffer.size());
  try
  {
    output = openForWrite(target);
    size_type remaining = source_values, group_size = 0, produced = 0, used = 0;
    std::uint8_t minimum = std::numeric_limits<std::uint8_t>::max();
    while(remaining > 0)
    {
      size_type bytes = std::min(remaining, static_cast<size_type>(input_buffer.size()));
      readExact(input, input_buffer.data(), bytes, source);
      for(size_type i = 0; i < bytes; i++)
      {
        minimum = std::min(minimum, input_buffer[i]); group_size++;
        if(group_size == branching)
        {
          output_buffer[used++] = minimum; produced++;
          if(used == output_buffer.size())
          {
            flushLevelOutput(output, output_buffer, used, target, checkpoint);
          }
          group_size = 0; minimum = std::numeric_limits<std::uint8_t>::max();
        }
      }
      remaining -= bytes;
    }
    // The final partial branch has a parent too. Omitting it changes both the
    // hierarchy offsets and the range minima above it.
    if(group_size > 0)
    {
      output_buffer[used++] = minimum; produced++;
    }
    flushLevelOutput(output, output_buffer, used, target, checkpoint);
    requireEnd(input, source);
    if(produced != target_values)
    {
      throw std::runtime_error("wrong number of values while building LCP level");
    }
    closeChecked(input, source); input = -1;
    closeChecked(output, target); output = -1;
  }
  catch(...)
  {
    if(input >= 0) { ::close(input); }
    if(output >= 0) { ::close(output); }
    throw;
  }
}

void
checkLevelLength(const std::string& filename, size_type expected)
{
  if(rawFileSize(filename) != expected)
  {
    throw std::runtime_error("LCP level length does not match its expected value count: " + filename);
  }
}

void
matchLegacyPadding(LCPArray& target, std::uint8_t maximum)
{
  const size_type values = target.values();
  const size_type width = bit_length(maximum);
  if(width == 8 || values > std::numeric_limits<size_type>::max() / BYTE_BITS)
  {
    return;
  }
  const size_type raw_bits = values * BYTE_BITS;
  const size_type raw_words = roundedDivision(raw_bits, WORD_BITS);
  if(raw_words > std::numeric_limits<size_type>::max() / WORD_BITS) { return; }
  const size_type raw_capacity = raw_words * WORD_BITS;
  if(raw_capacity != target.data.capacity()) { return; }

  // sdsl::util::bit_compress() preserves unused bits when reducing the width
  // without crossing an allocation boundary. Recreate that historical padding
  // exactly. This condition is limited to a sub-cacheline hierarchy, while
  // large construction never retains a raw multi-level copy.
  sdsl::int_vector<0> legacy(values, ~(std::uint8_t)0, 8);
  for(size_type i = 0; i < values; i++) { legacy[i] = target.data[i]; }
  sdsl::util::bit_compress(legacy);
  target.data.swap(legacy);
}

size_type
streamingBudget(const ConstructionParameters& parameters)
{
  return std::min(parameters.getMemoryLimitBytes(), parameters.getIOBufferSize());
}

std::uint64_t
rawFileChecksum(const std::string& filename, size_type bytes, size_type byte_budget)
{
  std::vector<std::uint8_t> buffer(byte_budget);
  int input = openForRead(filename);
  std::uint64_t checksum = 1469598103934665603ULL;
  try
  {
    size_type remaining = bytes;
    while(remaining > 0)
    {
      size_type current = std::min(remaining, static_cast<size_type>(buffer.size()));
      readExact(input, buffer.data(), current, filename);
      checksum = BuildWorkspace::checksum(buffer.data(), current, checksum);
      remaining -= current;
    }
    requireEnd(input, filename);
    closeChecked(input, filename); input = -1;
  }
  catch(...)
  {
    if(input >= 0) { ::close(input); }
    throw;
  }
  return checksum;
}

bool
hasWorkspaceManifest(const std::string& directory)
{
  const std::string manifest = directory + "/build.json";
  if(::access(manifest.c_str(), F_OK) == 0) { return true; }
  if(errno == ENOENT) { return false; }
  systemError("cannot access", manifest);
}

BuildWorkspace::Settings
lcpSemanticSettings(size_type branching, size_type leaf_bytes, std::uint64_t checksum)
{
  BuildWorkspace::Settings settings;
  settings["lcp_streaming"] = "v1";
  settings["lcp_branching"] = std::to_string(branching);
  settings["leaf_bytes"] = std::to_string(leaf_bytes);
  settings["leaf_checksum"] = std::to_string(checksum);
  return settings;
}

BuildWorkspace::Settings
lcpOperationalSettings(size_type byte_budget)
{
  BuildWorkspace::Settings settings;
  settings["stream_buffer_bytes"] = std::to_string(byte_budget);
  return settings;
}

} // namespace

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr size_type STNode::UNKNOWN;

//------------------------------------------------------------------------------

// Other class variables.

const std::string LCPArray::EXTENSION = ".lcp";

//------------------------------------------------------------------------------

std::ostream&
operator<< (std::ostream& out, const STNode& node)
{
  out << "(" << node.left_lcp << ", " << node.range() << " at depth " << node.lcp()
      << ", " << node.right_lcp << ")";
  return out;
}

//------------------------------------------------------------------------------

LCPArray::LCPArray()
{
}

LCPArray::LCPArray(const LCPArray& source)
{
  this->copy(source);
}

LCPArray::LCPArray(LCPArray&& source) noexcept
{
  *this = std::move(source);
}

LCPArray::~LCPArray()
{
}

void
LCPArray::copy(const LCPArray& source)
{
  this->header = source.header;
  this->data = source.data;
  this->offsets = source.offsets;
}

void
LCPArray::swap(LCPArray& another) noexcept
{
  if(this != &another)
  {
    this->header.swap(another.header);
    this->data.swap(another.data);
    this->offsets.swap(another.offsets);
  }
}

LCPArray&
LCPArray::operator=(const LCPArray& source)
{
  if(this != &source) { this->copy(source); }
  return *this;
}

LCPArray&
LCPArray::operator=(LCPArray&& source) noexcept
{
  if(this != &source)
  {
    this->header = std::move(source.header);
    this->data = std::move(source.data);
    this->offsets = std::move(source.offsets);
  }
  return *this;
}

LCPArray::size_type
LCPArray::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += this->header.serialize(out, child, "header");
  written_bytes += this->data.serialize(out, child, "data");
  written_bytes += this->offsets.serialize(out, child, "offsets");

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
LCPArray::load(std::istream& in)
{
  this->header.load(in);
  if(!(this->header.check()))
  {
    std::stringstream ss;
    ss << "LCP::load(): Invalid header: " << this->header;
    throw std::runtime_error(ss.str());
  }

  this->data.load(in);
  this->offsets.load(in);
}

//------------------------------------------------------------------------------

/*
  Tree operations on the range minimum tree. If level is required, it refers to the
  level of the parameter node.
*/

inline size_type
rmtRoot(const LCPArray& lcp)
{
  return lcp.values() - 1;
}

inline size_type
rmtParent(const LCPArray& lcp, size_type node, size_type level)
{
  return lcp.offsets[level + 1] + (node - lcp.offsets[level]) / lcp.branching();
}

inline bool
rmtIsFirst(const LCPArray& lcp, size_type node, size_type level)
{
  return ((node - lcp.offsets[level]) % lcp.branching() == 0);
}

inline size_type
rmtFirstSibling(const LCPArray& lcp, size_type node, size_type level)
{
  return node - (node - lcp.offsets[level]) % lcp.branching();
}

inline size_type
rmtLastSibling(const LCPArray& lcp, size_type first_child, size_type level)
{
  return std::min(lcp.offsets[level + 1], first_child + lcp.branching()) - 1;
}

inline size_type
rmtFirstChild(const LCPArray& lcp, size_type node, size_type level)
{
  return lcp.offsets[level - 1] + (node - lcp.offsets[level]) * lcp.branching();
}

inline size_type
rmtLastChild(const LCPArray& lcp, size_type node, size_type level)
{
  return rmtLastSibling(lcp, rmtFirstChild(lcp, node, level), level - 1);
}

inline size_type
rmtLevel(const LCPArray& lcp, size_type node)
{
  size_type level = 0;
  while(lcp.offsets[level + 1] <= node) { level++; }
  return level;
}

//------------------------------------------------------------------------------

LCPArray::LCPArray(const InputGraph& graph, const ConstructionParameters& parameters)
{
  double start = readTimer();

  if(graph.size() == 0) { return; }
  if(graph.lcp_name.empty())
  {
    std::cerr << "LCPArray::LCPArray(): The input graph does not contain the LCP file" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if(parameters.externalMemory())
  {
    const size_type byte_budget = streamingBudget(parameters);
    if(byte_budget < 2)
    {
      throw std::runtime_error("external LCP construction requires at least two bytes of stream budget");
    }
    const size_type leaf_bytes = rawFileSize(graph.lcp_name);
    const std::uint64_t leaf_checksum = rawFileChecksum(graph.lcp_name,
      leaf_bytes, byte_budget);
    const std::string workspace_directory = parameters.getWorkDirectory() + "/lcp-levels";
    BuildWorkspace::Settings semantic = lcpSemanticSettings(parameters.getLCPBranching(),
      leaf_bytes, leaf_checksum);
    BuildWorkspace::OpenMode mode = (parameters.getResume() &&
      hasWorkspaceManifest(workspace_directory) ? BuildWorkspace::RESUME :
      BuildWorkspace::NEW_WORKSPACE);
    BuildWorkspace workspace(workspace_directory, semantic,
      lcpOperationalSettings(byte_budget), mode);
    LCPArray streamed(graph.lcp_name, parameters.getLCPBranching(), byte_budget,
      &workspace);
    this->swap(streamed);
    return;
  }

  std::ifstream in(graph.lcp_name, std::ios_base::binary);
  if(!in)
  {
    std::cerr << "LCPArray::LCPArray(): Cannot open LCP file " << graph.lcp_name << std::endl;
    std::exit(EXIT_FAILURE);
  }

  this->header.branching = parameters.getLCPBranching();

  // Determine the number of levels.
  this->header.size = fileSize(in);
  size_type level_count = 1, level_size = this->size();
  while(level_size > 1)
  {
    level_count++; level_size = (level_size + this->branching() - 1) / this->branching();
  }

  // Initialize offsets.
  this->offsets = sdsl::int_vector<64>(level_count + 1, 0);
  level_size = this->size();
  size_type total_size = 0;
  for(size_type level = 0; level < this->levels(); level++)
  {
    total_size += level_size;
    this->offsets[level + 1] = total_size;
    level_size = (level_size + this->branching() - 1) / this->branching();
  }

  // Initialize data.
  this->data = sdsl::int_vector<0>(total_size, ~(uint8_t)0, 8);
  if(!DiskIO::read(in, reinterpret_cast<uint8_t*>(this->data.data()), this->size()))
  {
    std::cerr << "LCPArray::LCPArray(): Unexpected EOF in " << graph.lcp_name << std::endl;
    std::exit(EXIT_FAILURE);
  }
  in.close();
  for(size_type level = 0; level + 1 < this->levels(); level++)
  {
    for(size_type i = this->offsets[level]; i < this->offsets[level + 1]; i++)
    {
      size_type par = rmtParent(*this, i, level);
      if(this->data[i] < this->data[par]) { this->data[par] = this->data[i]; }
    }
  }
  sdsl::util::bit_compress(this->data);

  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    double seconds = readTimer() - start;
    std::cerr << "LCPArray::LCPArray(): Construction: " << seconds << " seconds, "
              << inGigabytes(memoryUsage()) << " GB" << std::endl;
  }
  if(Verbosity::level >= Verbosity::BASIC)
  {
    std::cerr << "LCPArray::LCPArray(): " << this->values() << " values at " << this->levels()
              << " levels (branching factor " << this->branching() << ")" << std::endl;
  }
}

//------------------------------------------------------------------------------

LCPArray::LCPArray(const std::string& leaf_filename, size_type branching,
  size_type byte_budget, BuildWorkspace* workspace, LCPStreamingStats* stats)
{
  LCPStreamingStats local_stats;
  LCPStreamingStats& telemetry = (stats == nullptr ? local_stats : *stats);
  telemetry = LCPStreamingStats();
  if(branching < 2)
  {
    throw std::invalid_argument("streaming LCP construction requires branching >= 2");
  }
  if(byte_budget < 2)
  {
    throw std::invalid_argument("streaming LCP construction requires a byte budget of at least 2");
  }

  this->header.branching = branching;
  this->header.size = rawFileSize(leaf_filename);
  std::vector<size_type> sizes = levelSizes(this->size(), branching);
  telemetry.levels = sizes.size();

  this->offsets = sdsl::int_vector<64>(sizes.size() + 1, 0);
  size_type total_size = 0;
  for(size_type level = 0; level < sizes.size(); level++)
  {
    if(total_size > std::numeric_limits<size_type>::max() - sizes[level])
    {
      throw std::runtime_error("LCP hierarchy size overflows the address space");
    }
    total_size += sizes[level]; this->offsets[level + 1] = total_size;
  }

  // All values in a range-minimum level are leaf values. Determining their
  // maximum first lets data be allocated directly in its final packed width,
  // avoiding an all-level uncompressed intermediate.
  std::uint8_t maximum = maximumValue(leaf_filename, this->size(), byte_budget, telemetry);
  this->data = sdsl::int_vector<0>(total_size,
    (maximum == std::numeric_limits<std::uint8_t>::max() ? ~(std::uint8_t)0 : 0),
    bit_length(maximum));
  copyIntoData(*this, this->offsets[0], leaf_filename, this->size(), byte_budget, telemetry);

  std::string current = leaf_filename;
  bool current_is_temporary = false;
  try
  {
    for(size_type level = 1; level < sizes.size(); level++)
    {
      const size_type target_values = sizes[level];
      const ArtifactIdentity identity = levelArtifact(level);
      std::string next = TempFile::getName("gcsa_lcp_level");
      try
      {
        if(workspace != nullptr && workspace->task_completed(identity.task, identity.phase))
        {
          // restore_artifact uses its own bounded copy buffer while no raw
          // level input/output buffers are live.
          workspace->restore_artifact(identity, logical_file_id_t(0),
            physical_shard_id_t(level), next, byte_budget);
          checkLevelLength(next, target_values); telemetry.restored_levels++;
        }
        else
        {
          if(workspace != nullptr)
          {
            BuildWorkspace::ArtifactWriter writer = workspace->open_artifact(identity,
              logical_file_id_t(0), physical_shard_id_t(level), "position", "all");
            buildLevel(current, sizes[level - 1], next, target_values, branching,
              byte_budget, &writer, telemetry);
            BuildWorkspace::ArtifactRef artifact = writer.finish(target_values);
            std::vector<BuildWorkspace::ArtifactRef> artifacts(1, artifact);
            workspace->commit_task(identity.task, identity.phase, artifacts);
            telemetry.checkpointed_levels++;
          }
          else
          {
            buildLevel(current, sizes[level - 1], next, target_values, branching,
              byte_budget, nullptr, telemetry);
          }
          checkLevelLength(next, target_values); telemetry.generated_levels++;
        }

        // Leaf-to-root offsets are fixed before the stream begins. Once this
        // level has populated final packed data, its raw predecessor is no
        // longer needed to generate a later level.
        copyIntoData(*this, this->offsets[level], next, target_values,
          byte_budget, telemetry);
      }
      catch(...)
      {
        TempFile::remove(next); throw;
      }
      if(current_is_temporary) { TempFile::remove(current); }
      current = next; current_is_temporary = true;
    }
    if(current_is_temporary) { TempFile::remove(current); }
  }
  catch(...)
  {
    if(current_is_temporary) { TempFile::remove(current); }
    throw;
  }
  matchLegacyPadding(*this, maximum);
}

//------------------------------------------------------------------------------

LCPArray::node_type
LCPArray::parent(const LCPArray::node_type& node) const
{
  if(node == this->root()) { return this->root(); }

  size_type node_lcp = std::max(node.left_lcp, node.right_lcp);
  range_type left(node.sp, node.left_lcp), right(node.ep + 1, node.right_lcp);
  if(node.left_lcp == node_lcp)
  {
    left = this->psv(node.sp);
    if(left == this->notFound()) { left = range_type(0, 0); }
  }
  if(node.right_lcp == node_lcp)
  {
    right = this->nsv(node.ep + 1);
    if(right == this->notFound()) { right = range_type(this->size(), 0); }
  }

  return node_type(left.first, right.first - 1, left.second, right.second, node_lcp);
}

LCPArray::node_type
LCPArray::parent(range_type range) const
{
  return this->parent(this->nodeFor(range));
}

//------------------------------------------------------------------------------

size_type
LCPArray::depth(const LCPArray::node_type& node) const
{
  if(node.lcp() != node_type::UNKNOWN) { return node.lcp(); }
  return this->depth(node.range());
}

size_type
LCPArray::depth(LCPArray::node_type& node) const
{
  if(node.lcp() == node_type::UNKNOWN) { node.node_lcp = this->depth(node.range()); }
  return node.lcp();
}

size_type
LCPArray::depth(range_type range) const
{
  if(Range::length(range) <= 1) { return node_type::UNKNOWN; }
  range_type res = this->rmq(range.first + 1, range.second);
  return (res == this->notFound() ? node_type::UNKNOWN : res.second);
}

//------------------------------------------------------------------------------

/*
  Find the last value less than 'val' between 'from' (inclusive) and 'to' (exclusive).
  Return value will be lcp.notFound() if no such value exists.
*/
template<class Comparator>
range_type
psv(const LCPArray& lcp, size_type from, size_type to, size_type val, const Comparator& comp)
{
  while(to > from)
  {
    to--;
    if(comp(lcp[to], val)) { return range_type(to, lcp[to]); }
  }
  return lcp.notFound();
}

template<class Comparator>
range_type
psv(const LCPArray& lcp, size_type to, const Comparator& comp)
{
  if(to == 0 || to >= lcp.size()) { return lcp.notFound(); }

  // Find the children of the lowest common ancestor of psv(to) and 'to'.
  size_type level = 0, val = lcp[to];
  range_type res = lcp.notFound();
  while(to != rmtRoot(lcp))
  {
    res = gcsa::psv(lcp, rmtFirstSibling(lcp, to, level), to, val, comp);
    if(res.first < lcp.values()) { break; }
    to = rmtParent(lcp, to, level); level++;
  }
  if(res.first >= lcp.values()) { return res; } // Not found.

  // Go to the leaf containing psv(to).
  while(level > 0)
  {
    size_type from = rmtFirstChild(lcp, res.first, level); level--;
    res = gcsa::psv(lcp, from, rmtLastSibling(lcp, from, level) + 1, val, comp);
  }

  return res;
}

range_type
LCPArray::psv(size_type pos) const
{
  return gcsa::psv(*this, pos, std::less<size_type>());
}

range_type
LCPArray::psev(size_type pos) const
{
  return gcsa::psv(*this, pos, std::less_equal<size_type>());
}

//------------------------------------------------------------------------------

/*
  Find the first value less than 'val' between 'from' and 'to' (inclusive).
  Return value will be lcp.notFound() if no such value exists.
*/
template<class Comparator>
range_type
nsv(const LCPArray& lcp, size_type from, size_type to, size_type val, const Comparator& comp)
{
  for(size_type i = from; i <= to; i++)
  {
    if(comp(lcp[i], val)) { return range_type(i, lcp[i]); }
  }
  return lcp.notFound();
}

template<class Comparator>
range_type
nsv(const LCPArray& lcp, size_type from, const Comparator& comp)
{
  if(from + 1 >= lcp.size()) { return lcp.notFound(); }

  // Find the children of the lowest common ancestor for 'from' and nsv(from).
  size_type level = 0, val = lcp[from];
  range_type res = lcp.notFound();
  while(from != rmtRoot(lcp))
  {
    res = gcsa::nsv(lcp, from + 1, rmtLastSibling(lcp, from, level), val, comp);
    if(res.first < lcp.values()) { break; }
    from = rmtParent(lcp, from, level); level++;
  }
  if(res.first >= lcp.values()) { return res; }

  // Go to the leaf containing nsv(to).
  while(level > 0)
  {
    from = rmtFirstChild(lcp, res.first, level); level--;
    res = gcsa::nsv(lcp, from, rmtLastSibling(lcp, from, level), val, comp);
  }

  return res;
}

range_type
LCPArray::nsv(size_type pos) const
{
  return gcsa::nsv(*this, pos, std::less<size_type>());
}

range_type
LCPArray::nsev(size_type pos) const
{
  return gcsa::nsv(*this, pos, std::less_equal<size_type>());
}

//------------------------------------------------------------------------------

inline void
updateRes(const LCPArray& lcp, range_type& res, size_type i)
{
  if(lcp.data[i] < res.second) { res.first = i; res.second = lcp.data[i]; }
}

range_type
LCPArray::rmq(size_type sp, size_type ep) const
{
  if(sp > ep || ep >= this->size()) { return this->notFound(); }
  if(sp == ep) { return range_type(sp, this->data[sp]); }

  /*
    Search for a subtree containing the rmq, maintaining the following invariants:
      - left < right
      - nodes before subtree(left) are processed
      - nodes after subtree(right) are in tail
      - res contains (i, lcp[i]) for the rmq in the processed range
  */
  range_type res(this->values(), this->size());
  size_type level = 0, left = sp, right = ep;
  std::stack<range_type> tail;
  while(true)
  {
    size_type left_par = rmtParent(*this, left, level), right_par = rmtParent(*this, right, level);
    if(left_par == right_par)
    {
      for(size_type i = left; i <= right; i++) { updateRes(*this, res, i); }
      break;
    }

    size_type left_child = rmtFirstChild(*this, left_par, level + 1);
    if(left != left_child)
    {
      size_type last_child = rmtLastSibling(*this, left_child, level);
      for(size_type i = left; i <= last_child; i++) { updateRes(*this, res, i); }
      left_par++;
    }

    size_type right_child = rmtLastChild(*this, right_par, level + 1);
    if(right != right_child)
    {
      size_type first_child = rmtFirstSibling(*this, right_child, level);
      for(size_type i = right; i >= first_child; i--) { tail.push(range_type(i, this->data[i])); }
      right_par--;
    }

    if(left_par >= right_par)
    {
      if(left_par == right_par) { updateRes(*this, res, left_par); }
      break;
    }
    left = left_par; right = right_par; level++;
  }

  // Check the tail.
  while(!(tail.empty()))
  {
    range_type temp = tail.top(); tail.pop();
    if(temp.second < res.second) { res = temp; }
  }

  // Find the leftmost leaf in subtree(res.first) containing LCP value res.second.
  level = rmtLevel(*this, res.first);
  while(level > 0)
  {
    res.first = rmtFirstChild(*this, res.first, level); level--;
    while(this->data[res.first] != res.second) { res.first++; }
  }

  return res;
}

range_type
LCPArray::rmq(range_type range) const
{
  return this->rmq(range.first, range.second);
}

//------------------------------------------------------------------------------

} // namespace gcsa
