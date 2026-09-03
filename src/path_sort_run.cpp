/*
  Copyright (c) 2026 GCSA2 contributors

  Prefix-compressed transient path-label sort runs.
*/

#include <gcsa/path_sort_run.h>

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <limits>
#include <stdexcept>
#include <sys/stat.h>
#include <type_traits>
#include <unistd.h>

namespace gcsa
{

namespace
{

constexpr std::uint64_t RUN_HEADER_MAGIC = 0x32524f5350534347ULL; // "GCSPSOR2"
constexpr std::uint64_t RUN_FOOTER_MAGIC = 0x32444e4550534347ULL; // "GCSPEND2"
constexpr std::uint32_t RUN_FORMAT_VERSION = 1;
constexpr std::uint32_t RUN_ENDIAN_MARKER = 0x01020304U;
constexpr size_type RUN_HEADER_BYTES = 8 + 4 + 4 + 8 + 8 + 8;
constexpr size_type RUN_FOOTER_BYTES = 8 + 4 + 4 + 8 + 8 + 8;
constexpr size_type RUN_MINIMUM_BUFFER_BYTES = 256;
constexpr off_t RUN_CACHE_TAIL_BYTES = 64 * MEGABYTE;
constexpr off_t RUN_CACHE_FLUSH_BYTES = 512 * MEGABYTE;
constexpr std::uint64_t FNV_OFFSET = 1469598103934665603ULL;
constexpr std::uint64_t FNV_PRIME = 1099511628211ULL;
static_assert(sizeof(node_type) == sizeof(std::uint64_t),
  "path-sort run format requires 64-bit node identifiers");
static_assert(sizeof(size_type) == sizeof(std::uint64_t),
  "path-sort run format requires 64-bit PathNode fields");

std::runtime_error
runError(const std::string& message, const std::string& filename)
{
  return std::runtime_error("PathSortRun: " + message +
    (filename.empty() ? std::string() : ": " + filename));
}

template<class Value>
void
encodeLittle(std::uint8_t* output, Value value)
{
  typedef typename std::make_unsigned<Value>::type unsigned_type;
  unsigned_type source = static_cast<unsigned_type>(value);
  for(size_type i = 0; i < sizeof(Value); i++)
  {
    output[i] = static_cast<std::uint8_t>(source & 0xFFU); source >>= 8;
  }
}

template<class Value>
Value
decodeLittle(const std::uint8_t* input)
{
  typedef typename std::make_unsigned<Value>::type unsigned_type;
  unsigned_type result = 0;
  for(size_type i = 0; i < sizeof(Value); i++)
  {
    result |= static_cast<unsigned_type>(input[i]) << (8 * i);
  }
  return static_cast<Value>(result);
}

void
updateChecksum(std::uint64_t& checksum, const void* data, size_type bytes)
{
  const std::uint8_t* input = static_cast<const std::uint8_t*>(data);
  for(size_type i = 0; i < bytes; i++)
  {
    checksum ^= input[i]; checksum *= FNV_PRIME;
  }
}

void
writeAll(int descriptor, const void* data, size_type bytes,
  const std::string& filename)
{
  const std::uint8_t* input = static_cast<const std::uint8_t*>(data);
  size_type written = 0;
  while(written < bytes)
  {
    ssize_t result = ::write(descriptor, input + written, bytes - written);
    if(result < 0 && errno == EINTR) { continue; }
    if(result <= 0) { throw runError("short write", filename); }
    written += static_cast<size_type>(result);
    DiskIO::write_volume += static_cast<size_type>(result);
  }
}

void
pwriteAll(int descriptor, const void* data, size_type bytes, off_t offset,
  const std::string& filename)
{
  const std::uint8_t* input = static_cast<const std::uint8_t*>(data);
  size_type written = 0;
  while(written < bytes)
  {
    ssize_t result = ::pwrite(descriptor, input + written, bytes - written,
      offset + static_cast<off_t>(written));
    if(result < 0 && errno == EINTR) { continue; }
    if(result <= 0) { throw runError("short header write", filename); }
    written += static_cast<size_type>(result);
    DiskIO::write_volume += static_cast<size_type>(result);
  }
}

void
adviseSequential(int descriptor)
{
#if defined(POSIX_FADV_SEQUENTIAL)
  static_cast<void>(::posix_fadvise(descriptor, 0, 0, POSIX_FADV_SEQUENTIAL));
#else
  static_cast<void>(descriptor);
#endif
}

void
discardCache(int descriptor, off_t offset, off_t bytes)
{
#if defined(POSIX_FADV_DONTNEED)
  if(bytes > 0)
  {
    // Cache advice is an optimization and must never alter construction
    // correctness on filesystems that do not implement it.
    static_cast<void>(::posix_fadvise(descriptor, offset, bytes,
      POSIX_FADV_DONTNEED));
  }
#else
  static_cast<void>(descriptor); static_cast<void>(offset);
  static_cast<void>(bytes);
#endif
}

void
appendEncoded(std::vector<std::uint8_t>& target, size_type& offset,
  std::uint64_t value, size_type bytes)
{
  for(size_type i = 0; i < bytes; i++)
  {
    target[offset++] = static_cast<std::uint8_t>(value & 0xFFU); value >>= 8;
  }
}

std::vector<std::uint8_t>
runHeader(size_type records, size_type payload_bytes, std::uint64_t checksum)
{
  std::vector<std::uint8_t> result(RUN_HEADER_BYTES); size_type offset = 0;
  appendEncoded(result, offset, RUN_HEADER_MAGIC, 8);
  appendEncoded(result, offset, RUN_FORMAT_VERSION, 4);
  appendEncoded(result, offset, RUN_ENDIAN_MARKER, 4);
  appendEncoded(result, offset, records, 8);
  appendEncoded(result, offset, payload_bytes, 8);
  appendEncoded(result, offset, checksum, 8);
  return result;
}

std::vector<std::uint8_t>
runFooter(size_type records, size_type payload_bytes, std::uint64_t checksum)
{
  std::vector<std::uint8_t> result(RUN_FOOTER_BYTES); size_type offset = 0;
  appendEncoded(result, offset, RUN_FOOTER_MAGIC, 8);
  appendEncoded(result, offset, RUN_FORMAT_VERSION, 4);
  appendEncoded(result, offset, 0, 4);
  appendEncoded(result, offset, records, 8);
  appendEncoded(result, offset, payload_bytes, 8);
  appendEncoded(result, offset, checksum, 8);
  return result;
}

size_type
checkedAdd(size_type left, size_type right, const std::string& filename)
{
  if(left > std::numeric_limits<size_type>::max() - right)
  {
    throw runError("byte count overflow", filename);
  }
  return left + right;
}

size_type
checkedBufferBytes(size_type requested, const std::string& filename)
{
  if(requested < RUN_MINIMUM_BUFFER_BYTES)
  {
    throw runError("I/O buffer is below the format minimum", filename);
  }
  return requested;
}

} // namespace

PathSortRunWriter::PathSortRunWriter(const std::string& filename,
  size_type expected_records, size_type buffer_bytes) :
  filename_(filename), file_(-1),
  buffer_(checkedBufferBytes(buffer_bytes, filename)), used_(0),
  expected_records_(expected_records), records_(0), payload_bytes_(0),
  total_bytes_(0), checksum_(FNV_OFFSET), cache_released_(0),
  previous_ranks_(0), finished_(false)
{
  this->file_ = ::open(this->filename_.c_str(),
    O_CREAT | O_TRUNC | O_RDWR, 0600);
  if(this->file_ < 0) { throw runError("cannot create run", this->filename_); }
  try
  {
    adviseSequential(this->file_);
    std::fill(this->previous_, this->previous_ + PathLabel::LABEL_LENGTH + 1,
      PathLabel::NO_RANK);
    std::vector<std::uint8_t> header = runHeader(this->expected_records_, 0, 0);
    this->append(header.data(), header.size(), false);
  }
  catch(...)
  {
    ::close(this->file_); this->file_ = -1;
    ::unlink(this->filename_.c_str()); throw;
  }
}

PathSortRunWriter::~PathSortRunWriter()
{
  if(this->file_ >= 0) { ::close(this->file_); }
  if(!(this->finished_)) { ::unlink(this->filename_.c_str()); }
}

void
PathSortRunWriter::append(const void* data, size_type bytes, bool payload)
{
  if(bytes == 0) { return; }
  if(this->finished_) { throw runError("write after finish", this->filename_); }
  if(payload)
  {
    this->payload_bytes_ = checkedAdd(this->payload_bytes_, bytes, this->filename_);
    updateChecksum(this->checksum_, data, bytes);
  }
  const std::uint8_t* input = static_cast<const std::uint8_t*>(data);
  while(bytes > 0)
  {
    if(this->used_ == this->buffer_.size()) { this->flush(false); }
    size_type copied = std::min(bytes, this->buffer_.size() - this->used_);
    std::memcpy(this->buffer_.data() + this->used_, input, copied);
    this->used_ += copied; input += copied; bytes -= copied;
  }
}

void
PathSortRunWriter::flush(bool durable)
{
  if(this->used_ > 0)
  {
    writeAll(this->file_, this->buffer_.data(), this->used_, this->filename_);
    this->total_bytes_ = checkedAdd(this->total_bytes_, this->used_, this->filename_);
    this->used_ = 0;
  }
  if(!durable &&
     static_cast<off_t>(this->total_bytes_) - this->cache_released_ < RUN_CACHE_FLUSH_BYTES)
  {
    return;
  }
  if(::fdatasync(this->file_) != 0)
  {
    throw runError("cannot sync run", this->filename_);
  }
  off_t written = static_cast<off_t>(this->total_bytes_);
  off_t discard_end = (durable ? written :
    std::max(this->cache_released_, written - RUN_CACHE_TAIL_BYTES));
  if(discard_end > this->cache_released_)
  {
    discardCache(this->file_, this->cache_released_,
      discard_end - this->cache_released_);
    this->cache_released_ = discard_end;
  }
}

void
PathSortRunWriter::write(const PathSortRunRecord& source)
{
  if(this->records_ >= this->expected_records_)
  {
    throw runError("too many records", this->filename_);
  }
  if(source.node.ranks() > PathLabel::LABEL_LENGTH + 1 ||
     source.node.lcp() > source.node.order())
  {
    throw runError("invalid path order or lcp", this->filename_);
  }

  PathNode node = source.node; node.setPointer(0);
  std::uint8_t encoded[sizeof(std::uint64_t)];
  encodeLittle<std::uint64_t>(encoded, node.from);
  this->append(encoded, sizeof(node.from), true);
  encodeLittle<std::uint64_t>(encoded, node.to);
  this->append(encoded, sizeof(node.to), true);
  encodeLittle<std::uint64_t>(encoded, node.fields);
  this->append(encoded, sizeof(node.fields), true);

  size_type ranks = node.ranks(), common = 0;
  while(common < ranks && common < this->previous_ranks_ &&
        source.labels[common] == this->previous_[common]) { common++; }
  std::uint8_t prefix = static_cast<std::uint8_t>(common);
  this->append(&prefix, sizeof(prefix), true);
  std::uint8_t rank_bytes[sizeof(PathNode::rank_type)];
  for(size_type i = common; i < ranks; i++)
  {
    encodeLittle<PathNode::rank_type>(rank_bytes, source.labels[i]);
    this->append(rank_bytes, sizeof(rank_bytes), true);
  }
  std::copy(source.labels, source.labels + ranks, this->previous_);
  this->previous_ranks_ = ranks; this->records_++;
}

void
PathSortRunWriter::writeHeader()
{
  std::vector<std::uint8_t> header = runHeader(this->records_,
    this->payload_bytes_, this->checksum_);
  pwriteAll(this->file_, header.data(), header.size(), 0, this->filename_);
}

void
PathSortRunWriter::finish()
{
  if(this->finished_) { return; }
  if(this->records_ != this->expected_records_)
  {
    throw runError("record count mismatch", this->filename_);
  }
  std::vector<std::uint8_t> footer = runFooter(this->records_,
    this->payload_bytes_, this->checksum_);
  this->append(footer.data(), footer.size(), false);
  this->flush(true);
  this->writeHeader();
  if(::fdatasync(this->file_) != 0)
  {
    throw runError("cannot sync final header", this->filename_);
  }
  discardCache(this->file_, 0, static_cast<off_t>(this->total_bytes_));
  this->finished_ = true;
}

size_type
PathSortRunWriter::uncompressedBytes() const
{
  if(this->records_ > std::numeric_limits<size_type>::max() / sizeof(PathSortRunRecord))
  {
    throw runError("uncompressed byte count overflow", this->filename_);
  }
  return this->records_ * sizeof(PathSortRunRecord);
}

PathSortRunReader::PathSortRunReader(const std::string& filename,
  size_type buffer_bytes) :
  filename_(filename), file_(-1),
  buffer_(checkedBufferBytes(buffer_bytes, filename)),
  begin_(0), end_(0), file_offset_(0), total_records_(0), records_read_(0),
  payload_bytes_(0), payload_read_(0), total_bytes_(0), expected_checksum_(0),
  checksum_(FNV_OFFSET), cache_released_(0), previous_ranks_(0),
  at_end_(false), footer_checked_(false)
{
  this->file_ = ::open(this->filename_.c_str(), O_RDONLY);
  if(this->file_ < 0) { throw runError("cannot open run", this->filename_); }
  try
  {
    adviseSequential(this->file_);
    std::fill(this->previous_, this->previous_ + PathLabel::LABEL_LENGTH + 1,
      PathLabel::NO_RANK);
    this->readHeader();
    if(this->total_records_ == 0)
    {
      this->readFooter(); this->at_end_ = true;
    }
    else { this->readRecord(); }
  }
  catch(...)
  {
    ::close(this->file_); this->file_ = -1; throw;
  }
}

PathSortRunReader::~PathSortRunReader()
{
  if(this->file_ >= 0) { ::close(this->file_); }
}

void
PathSortRunReader::refill()
{
  if(this->begin_ != this->end_)
  {
    throw runError("internal read buffer is not empty", this->filename_);
  }
  if(static_cast<off_t>(this->file_offset_) - this->cache_released_ >=
     RUN_CACHE_FLUSH_BYTES)
  {
    off_t discard_end = std::max(this->cache_released_,
      static_cast<off_t>(this->file_offset_) - RUN_CACHE_TAIL_BYTES);
    discardCache(this->file_, this->cache_released_,
      discard_end - this->cache_released_);
    this->cache_released_ = discard_end;
  }
  ssize_t result;
  do { result = ::read(this->file_, this->buffer_.data(), this->buffer_.size()); }
  while(result < 0 && errno == EINTR);
  if(result <= 0) { throw runError("unexpected end of run", this->filename_); }
  DiskIO::read_volume += static_cast<size_type>(result);
  this->file_offset_ = checkedAdd(this->file_offset_,
    static_cast<size_type>(result), this->filename_);
  this->begin_ = 0; this->end_ = static_cast<size_type>(result);
}

void
PathSortRunReader::readExact(void* target, size_type bytes, bool payload)
{
  std::uint8_t* output = static_cast<std::uint8_t*>(target);
  size_type copied = 0;
  while(copied < bytes)
  {
    if(this->begin_ == this->end_) { this->refill(); }
    size_type chunk = std::min(bytes - copied, this->end_ - this->begin_);
    std::memcpy(output + copied, this->buffer_.data() + this->begin_, chunk);
    this->begin_ += chunk; copied += chunk;
  }
  if(payload)
  {
    if(this->payload_read_ > this->payload_bytes_ ||
       bytes > this->payload_bytes_ - this->payload_read_)
    {
      throw runError("payload length exceeded", this->filename_);
    }
    this->payload_read_ += bytes;
    updateChecksum(this->checksum_, target, bytes);
  }
}

void
PathSortRunReader::readHeader()
{
  std::uint8_t header[RUN_HEADER_BYTES]; this->readExact(header, sizeof(header), false);
  size_type offset = 0;
  std::uint64_t magic = decodeLittle<std::uint64_t>(header + offset); offset += 8;
  std::uint32_t version = decodeLittle<std::uint32_t>(header + offset); offset += 4;
  std::uint32_t endian = decodeLittle<std::uint32_t>(header + offset); offset += 4;
  this->total_records_ = decodeLittle<std::uint64_t>(header + offset); offset += 8;
  this->payload_bytes_ = decodeLittle<std::uint64_t>(header + offset); offset += 8;
  this->expected_checksum_ = decodeLittle<std::uint64_t>(header + offset);
  if(magic != RUN_HEADER_MAGIC || version != RUN_FORMAT_VERSION ||
     endian != RUN_ENDIAN_MARKER)
  {
    throw runError("incompatible run header", this->filename_);
  }
  this->total_bytes_ = checkedAdd(checkedAdd(RUN_HEADER_BYTES,
    this->payload_bytes_, this->filename_), RUN_FOOTER_BYTES, this->filename_);
  struct stat info;
  if(::fstat(this->file_, &info) != 0 || info.st_size < 0 ||
     static_cast<size_type>(info.st_size) != this->total_bytes_)
  {
    throw runError("run length does not match header", this->filename_);
  }
}

void
PathSortRunReader::readRecord()
{
  if(this->records_read_ >= this->total_records_)
  {
    throw runError("too many run records", this->filename_);
  }
  std::uint8_t encoded[sizeof(std::uint64_t)];
  this->readExact(encoded, sizeof(node_type), true);
  this->current_.node.from = decodeLittle<std::uint64_t>(encoded);
  this->readExact(encoded, sizeof(node_type), true);
  this->current_.node.to = decodeLittle<std::uint64_t>(encoded);
  this->readExact(encoded, sizeof(size_type), true);
  this->current_.node.fields = decodeLittle<std::uint64_t>(encoded);
  this->current_.node.setPointer(0);
  if(this->current_.node.ranks() > PathLabel::LABEL_LENGTH + 1 ||
     this->current_.node.lcp() > this->current_.node.order())
  {
    throw runError("invalid path order or lcp", this->filename_);
  }
  std::uint8_t prefix = 0; this->readExact(&prefix, sizeof(prefix), true);
  size_type ranks = this->current_.node.ranks();
  if(prefix > ranks || prefix > this->previous_ranks_)
  {
    throw runError("invalid label prefix", this->filename_);
  }
  std::copy(this->previous_, this->previous_ + prefix, this->current_.labels);
  std::uint8_t rank_bytes[sizeof(PathNode::rank_type)];
  for(size_type i = prefix; i < ranks; i++)
  {
    this->readExact(rank_bytes, sizeof(rank_bytes), true);
    this->current_.labels[i] = decodeLittle<PathNode::rank_type>(rank_bytes);
  }
  std::fill(this->current_.labels + ranks,
    this->current_.labels + PathLabel::LABEL_LENGTH + 1, PathLabel::NO_RANK);
  std::copy(this->current_.labels, this->current_.labels + ranks, this->previous_);
  this->previous_ranks_ = ranks; this->records_read_++;
}

void
PathSortRunReader::readFooter()
{
  if(this->footer_checked_) { return; }
  if(this->payload_read_ != this->payload_bytes_ ||
     this->checksum_ != this->expected_checksum_)
  {
    throw runError("payload checksum or length mismatch", this->filename_);
  }
  std::uint8_t footer[RUN_FOOTER_BYTES]; this->readExact(footer, sizeof(footer), false);
  size_type offset = 0;
  std::uint64_t magic = decodeLittle<std::uint64_t>(footer + offset); offset += 8;
  std::uint32_t version = decodeLittle<std::uint32_t>(footer + offset); offset += 4;
  offset += 4;
  size_type records = decodeLittle<std::uint64_t>(footer + offset); offset += 8;
  size_type payload = decodeLittle<std::uint64_t>(footer + offset); offset += 8;
  std::uint64_t checksum = decodeLittle<std::uint64_t>(footer + offset);
  if(magic != RUN_FOOTER_MAGIC || version != RUN_FORMAT_VERSION ||
     records != this->total_records_ || payload != this->payload_bytes_ ||
     checksum != this->expected_checksum_)
  {
    throw runError("invalid run footer", this->filename_);
  }
  discardCache(this->file_, 0, static_cast<off_t>(this->total_bytes_));
  this->footer_checked_ = true;
}

const PathSortRunRecord&
PathSortRunReader::current() const
{
  if(this->at_end_) { throw runError("read past end", this->filename_); }
  return this->current_;
}

void
PathSortRunReader::advance()
{
  if(this->at_end_) { throw runError("advance past end", this->filename_); }
  if(this->records_read_ == this->total_records_)
  {
    this->readFooter(); this->at_end_ = true;
  }
  else { this->readRecord(); }
}

size_type
pathSortRunMinimumBuffer()
{
  return RUN_MINIMUM_BUFFER_BYTES;
}

} // namespace gcsa
