#include <gcsa/lcp.h>
#include <external_configuration.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <cerrno>
#include <cstring>
#include <filesystem>
#include <fcntl.h>
#include <fstream>
#include <limits>
#include <stack>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>

#include <gcsa/internal.h>

namespace gcsa
{

static void buildLCPAndStore(const std::string& leaf_filename,
  size_type branching, size_type byte_budget, const std::string& filename);

namespace
{

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
    DiskIO::read_volume += static_cast<size_type>(got);
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
    if(got > 0)
    {
      DiskIO::read_volume += static_cast<size_type>(got);
      throw std::runtime_error("unexpected trailing data in " + filename);
    }
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
    DiskIO::write_volume += static_cast<size_type>(written);
    offset += static_cast<size_type>(written);
  }
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
maximumValue(const std::string& filename, size_type values, size_type byte_budget)
{
  std::vector<std::uint8_t> buffer(byte_budget);
  int input = openForRead(filename);
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
flushLevelOutput(int output, std::vector<std::uint8_t>& buffer, size_type& used,
  const std::string& filename)
{
  if(used == 0) { return; }
  writeExact(output, buffer.data(), used, filename);
  used = 0;
}

void
buildLevel(const std::string& source, size_type source_values,
  const std::string& target, size_type target_values, size_type branching,
  size_type byte_budget)
{
  const size_type input_bytes = byte_budget / 2;
  const size_type output_bytes = byte_budget - input_bytes;
  std::vector<std::uint8_t> input_buffer(input_bytes), output_buffer(output_bytes);
  int input = openForRead(source), output = -1;
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
            flushLevelOutput(output, output_buffer, used, target);
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
    flushLevelOutput(output, output_buffer, used, target);
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

/*
  Visit raw levels in final leaf-to-root serialization order while retaining
  only the current and next raw files. The consumer must finish with a level
  before returning; a later iteration may remove that temporary file.
*/
template<class LevelConsumer>
void
forEachRawLevel(const std::string& leaf_filename,
  const std::vector<size_type>& sizes, size_type branching,
  size_type generation_budget, LevelConsumer consume)
{
  std::string current = leaf_filename;
  bool current_is_temporary = false;
  try
  {
    consume(0, current, sizes[0]);
    for(size_type level = 1; level < sizes.size(); level++)
    {
      const size_type target_values = sizes[level];
      std::string next = TempFile::getName("gcsa_lcp_level");
      try
      {
        buildLevel(current, sizes[level - 1], next, target_values, branching,
          generation_budget);
        checkLevelLength(next, target_values);
        consume(level, next, target_values);
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
}

size_type
streamingBudget(const ConstructionParameters& parameters)
{
  return std::min(parameters.getMemoryLimitBytes(), externalIOBufferSize(parameters));
}

bool
usesLegacyPadding(size_type values, size_type width)
{
  if(width == BYTE_BITS || values > std::numeric_limits<size_type>::max() / BYTE_BITS ||
     values > std::numeric_limits<size_type>::max() / width)
  {
    return false;
  }
  return roundedDivision(values * BYTE_BITS, WORD_BITS) ==
    roundedDivision(values * width, WORD_BITS);
}

/*
  SDSL's legacy constructor first stored the hierarchy at width eight in an
  all-ones allocation and then compressed it in place. If compression keeps
  the same allocation-sized word count, the unused tail bits retain pieces of
  that old byte representation. The equality above implies

      (8 - width) * values < 64,

  so the complete compatibility case contains fewer than 64 values and at
  most eight words. Reproduce it in fixed stack storage instead of allocating
  an uncharged int_vector that could exceed a tiny configured stream budget.
*/
class LegacyPackedLCPWriter
{
public:
  LegacyPackedLCPWriter(std::ostream& output, size_type expected_values,
    size_type value_width) :
    out(output), expected(expected_values), width(value_width), values(0),
    words()
  {
    if(this->width == 0 || this->width >= BYTE_BITS ||
       !usesLegacyPadding(this->expected, this->width))
    {
      throw std::runtime_error("invalid legacy packed LCP dimensions");
    }
    this->words.fill(std::numeric_limits<std::uint64_t>::max());
    sdsl::int_vector<0>::write_header(this->expected * this->width,
      static_cast<std::uint8_t>(this->width), this->out);
    if(!(this->out)) { throw std::runtime_error("cannot write packed LCP header"); }
  }

  void write(std::uint8_t value)
  {
    if(this->values >= this->expected ||
       value >= (static_cast<size_type>(1) << this->width))
    {
      throw std::runtime_error("LCP value exceeds its legacy packed representation");
    }
    const size_type bit_position = this->values * BYTE_BITS;
    const size_type word = bit_position / WORD_BITS;
    const size_type offset = bit_position % WORD_BITS;
    const std::uint64_t mask = std::numeric_limits<std::uint8_t>::max();
    this->words[word] &= ~(mask << offset);
    this->words[word] |= static_cast<std::uint64_t>(value) << offset;
    this->values++;
  }

  void finish()
  {
    if(this->values != this->expected)
    {
      throw std::runtime_error("legacy packed LCP value count mismatch");
    }

    // This is the loop used by sdsl::util::bit_compress(). Reading at width
    // eight and writing at the smaller width in place is safe because the
    // write cursor never overtakes the read cursor.
    const std::uint64_t* read_data = this->words.data();
    std::uint64_t* write_data = this->words.data();
    std::uint8_t read_offset = 0, write_offset = 0;
    for(size_type i = 0; i < this->values; i++)
    {
      std::uint64_t value = sdsl::bits::read_int_and_move(read_data,
        read_offset, BYTE_BITS);
      sdsl::bits::write_int_and_move(write_data, value,
        write_offset, static_cast<std::uint8_t>(this->width));
    }

    const size_type output_words = roundedDivision(this->expected * this->width,
      WORD_BITS);
    this->out.write(reinterpret_cast<const char*>(this->words.data()),
      output_words * sizeof(std::uint64_t));
    if(!(this->out)) { throw std::runtime_error("cannot write legacy packed LCP data"); }
  }

private:
  std::ostream& out;
  size_type expected, width, values;
  std::array<std::uint64_t, 8> words;
};

/*
  Stream the payload of an sdsl::int_vector<0>. Unused bits are initialized in
  the same way as the resident constructor: zero normally and one when the
  maximum value requires the full eight-bit width.
*/
class PackedLCPWriter
{
public:
  PackedLCPWriter(std::ostream& output, size_type expected_values,
    size_type value_width, bool one_padding, size_type buffer_bytes) :
    out(output), expected(expected_values), width(value_width), values(0),
    word(one_padding ? std::numeric_limits<std::uint64_t>::max() : 0),
    initial_word(word), bit_offset(0), words_written(0),
    buffer(buffer_bytes / sizeof(std::uint64_t)), buffered(0)
  {
    if(this->width == 0 || this->width > BYTE_BITS ||
       this->expected > std::numeric_limits<size_type>::max() / this->width)
    {
      throw std::runtime_error("invalid packed LCP dimensions");
    }
    sdsl::int_vector<0>::write_header(this->expected * this->width,
      static_cast<std::uint8_t>(this->width), this->out);
    if(!(this->out)) { throw std::runtime_error("cannot write packed LCP header"); }
  }

  size_type bufferBytes() const
  {
    return this->buffer.size() * sizeof(std::uint64_t);
  }

  void write(std::uint8_t value)
  {
    if(this->values >= this->expected ||
       (this->width < BYTE_BITS && value >= (static_cast<size_type>(1) << this->width)))
    {
      throw std::runtime_error("LCP value exceeds its packed representation");
    }
    std::uint64_t remaining_value = value;
    size_type remaining_bits = this->width;
    while(remaining_bits > 0)
    {
      size_type available = WORD_BITS - this->bit_offset;
      size_type take = std::min(remaining_bits, available);
      std::uint64_t mask = (static_cast<std::uint64_t>(1) << take) - 1;
      this->word &= ~(mask << this->bit_offset);
      this->word |= (remaining_value & mask) << this->bit_offset;
      remaining_value >>= take;
      remaining_bits -= take;
      this->bit_offset += take;
      if(this->bit_offset == WORD_BITS)
      {
        this->writeWord(this->word);
        this->word = this->initial_word; this->bit_offset = 0;
      }
    }
    this->values++;
  }

  void finish()
  {
    if(this->values != this->expected)
    {
      throw std::runtime_error("packed LCP value count mismatch");
    }
    if(this->bit_offset > 0) { this->writeWord(this->word); }
    this->flush();
    size_type expected_words = roundedDivision(this->expected * this->width,
      WORD_BITS);
    if(this->words_written != expected_words)
    {
      throw std::runtime_error("packed LCP word count mismatch");
    }
  }

private:
  void writeWord(std::uint64_t value)
  {
    if(this->buffer.empty())
    {
      this->out.write(reinterpret_cast<const char*>(&value), sizeof(value));
      if(!(this->out)) { throw std::runtime_error("cannot write packed LCP data"); }
    }
    else
    {
      this->buffer[this->buffered++] = value;
      if(this->buffered == this->buffer.size()) { this->flush(); }
    }
    this->words_written++;
  }

  void flush()
  {
    if(this->buffered == 0) { return; }
    this->out.write(reinterpret_cast<const char*>(this->buffer.data()),
      this->buffered * sizeof(std::uint64_t));
    if(!(this->out)) { throw std::runtime_error("cannot write packed LCP data"); }
    this->buffered = 0;
  }

  std::ostream& out;
  size_type expected, width, values;
  std::uint64_t word, initial_word;
  size_type bit_offset, words_written;
  std::vector<std::uint64_t> buffer;
  size_type buffered;
};

template<class ValueConsumer>
void
consumeRawLevel(const std::string& filename, size_type values,
  size_type byte_budget, ValueConsumer consume)
{
  if(byte_budget == 0) { throw std::runtime_error("zero LCP input buffer budget"); }
  std::vector<std::uint8_t> buffer(byte_budget);
  int input = openForRead(filename);
  try
  {
    size_type remaining = values;
    while(remaining > 0)
    {
      size_type bytes = std::min(remaining, static_cast<size_type>(buffer.size()));
      readExact(input, buffer.data(), bytes, filename);
      for(size_type i = 0; i < bytes; i++) { consume(buffer[i]); }
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
}

std::atomic<std::uint64_t> lcp_store_counter(0);

template<class Writer>
void
publishLCP(const std::string& filename, Writer writer)
{
  std::string partial = filename + "." +
    std::to_string(static_cast<unsigned long long>(::getpid())) + "." +
    std::to_string(lcp_store_counter.fetch_add(1)) + ".partial";
  std::ofstream out(partial.c_str(), std::ios_base::binary | std::ios_base::trunc);
  if(!out) { throw std::runtime_error("cannot create partial LCP array " + partial); }
  try
  {
    writer(out);
    const std::streampos end_position = out.tellp();
    const std::streamoff output_size = end_position - std::streampos(0);
    if(output_size < 0 || static_cast<std::uintmax_t>(output_size) >
       std::numeric_limits<size_type>::max())
    {
      throw std::runtime_error("cannot determine partial LCP array size " + partial);
    }
    const size_type output_bytes = static_cast<size_type>(output_size);
    out.flush(); out.close();
    if(!out) { throw std::runtime_error("cannot finish partial LCP array " + partial); }
    int descriptor = ::open(partial.c_str(), O_RDONLY);
    if(descriptor < 0 || ::fdatasync(descriptor) != 0)
    {
      if(descriptor >= 0) { ::close(descriptor); }
      throw std::runtime_error("cannot sync partial LCP array " + partial);
    }
    if(::close(descriptor) != 0 || ::rename(partial.c_str(), filename.c_str()) != 0)
    {
      throw std::runtime_error("cannot publish LCP array " + filename);
    }
    std::filesystem::path parent = std::filesystem::path(filename).parent_path();
    if(parent.empty()) { parent = "."; }
    descriptor = ::open(parent.c_str(), O_RDONLY | O_DIRECTORY);
    if(descriptor < 0 || ::fsync(descriptor) != 0)
    {
      if(descriptor >= 0) { ::close(descriptor); }
      throw std::runtime_error("cannot sync LCP output directory " + filename);
    }
    if(::close(descriptor) != 0)
    {
      throw std::runtime_error("cannot close LCP output directory " + filename);
    }
    DiskIO::write_volume += output_bytes;
  }
  catch(...)
  {
    out.close(); ::unlink(partial.c_str()); throw;
  }
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

void
LCPArray::buildAndStore(const InputGraph& graph,
  const ConstructionParameters& parameters, const std::string& filename)
{
  if(!(parameters.externalMemory()))
  {
    throw std::invalid_argument("LCPArray::buildAndStore() requires external-memory construction");
  }
  if(graph.size() == 0)
  {
    LCPArray empty;
    publishLCP(filename, [&](std::ostream& out) { empty.serialize(out); });
    return;
  }
  if(graph.lcp_name.empty())
  {
    throw std::runtime_error("LCPArray::buildAndStore(): input graph has no LCP leaf file");
  }

  const size_type byte_budget = streamingBudget(parameters);
  if(byte_budget < 2)
  {
    throw std::runtime_error("external LCP construction requires at least two bytes of stream budget");
  }
  buildLCPAndStore(graph.lcp_name, parameters.getLCPBranching(), byte_budget,
    filename);
}

static void
buildLCPAndStore(const std::string& leaf_filename,
  size_type branching, size_type byte_budget, const std::string& filename)
{
  if(branching < 2)
  {
    throw std::invalid_argument("streaming LCP construction requires branching >= 2");
  }
  if(byte_budget < 2)
  {
    throw std::invalid_argument("streaming LCP construction requires a byte budget of at least 2");
  }

  const size_type leaves = rawFileSize(leaf_filename);
  const std::vector<size_type> sizes = levelSizes(leaves, branching);
  sdsl::int_vector<64> offsets(sizes.size() + 1, 0);
  size_type total_values = 0;
  for(size_type level = 0; level < sizes.size(); level++)
  {
    if(total_values > std::numeric_limits<size_type>::max() - sizes[level])
    {
      throw std::runtime_error("LCP hierarchy size overflows the address space");
    }
    total_values += sizes[level]; offsets[level + 1] = total_values;
  }

  const std::uint8_t maximum = maximumValue(leaf_filename, leaves, byte_budget);
  const size_type width = bit_length(maximum);
  const bool legacy_padding = usesLegacyPadding(total_values, width);
  LCPHeader header;
  header.size = leaves; header.branching = branching;

  publishLCP(filename, [&](std::ostream& out)
  {
    header.serialize(out);
    if(!out) { throw std::runtime_error("cannot write LCP file header"); }

    if(legacy_padding)
    {
      // When the old eight-bit and final packed forms occupy the same words,
      // bit_compress() retains padding from the all-ones source. The fixed-size
      // writer preserves those bits without a hidden heap allocation.
      LegacyPackedLCPWriter data(out, total_values, width);
      forEachRawLevel(leaf_filename, sizes, branching, byte_budget,
        [&](size_type, const std::string& level_file, size_type values)
      {
        consumeRawLevel(level_file, values, byte_budget,
          [&](std::uint8_t value) { data.write(value); });
      });
      data.finish();
    }
    else
    {
      size_type output_buffer_bytes = 0;
      if(byte_budget >= 2 * sizeof(std::uint64_t))
      {
        output_buffer_bytes = (byte_budget / 2 / sizeof(std::uint64_t)) *
          sizeof(std::uint64_t);
      }
      PackedLCPWriter data(out, total_values, width,
        maximum == std::numeric_limits<std::uint8_t>::max(),
        output_buffer_bytes);
      const size_type input_budget = byte_budget - data.bufferBytes();
      if(input_budget < 2)
      {
        throw std::runtime_error("LCP stream budget cannot admit level generation");
      }
      forEachRawLevel(leaf_filename, sizes, branching, input_budget,
        [&](size_type, const std::string& level_file, size_type values)
      {
        consumeRawLevel(level_file, values, input_budget,
          [&](std::uint8_t value) { data.write(value); });
      });
      data.finish();
    }
    offsets.serialize(out);
    if(!out) { throw std::runtime_error("cannot write LCP offsets"); }
  });
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
