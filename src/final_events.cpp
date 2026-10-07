/*
  Copyright (c) 2026 GCSA2 contributors

  External event streams for final GCSA construction.
*/

#include <gcsa/final_events.h>
#include <external_configuration.hpp>

#include <compressed_block.hpp>
#include <gcsa/external_sort.h>
#include <gcsa/internal.h>

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <queue>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>

namespace gcsa
{

namespace
{

constexpr off_t CACHE_FLUSH_BYTES = 512 * 1024 * 1024;
constexpr off_t CACHE_TAIL_BYTES = 64 * 1024 * 1024;

std::runtime_error
eventError(const std::string& message, const std::string& path = std::string())
{
  return std::runtime_error("final GCSA events: " + message +
    (path.empty() ? std::string() : ": " + path));
}

void
put64(std::uint8_t* target, std::uint64_t value)
{
  for(size_type i = 0; i < 8; i++)
  {
    target[i] = static_cast<std::uint8_t>(value >> (8 * i));
  }
}

std::uint64_t
get64(const std::uint8_t* source)
{
  std::uint64_t result = 0;
  for(size_type i = 0; i < 8; i++)
  {
    result |= static_cast<std::uint64_t>(source[i]) << (8 * i);
  }
  return result;
}

void
writeAll(int descriptor, const void* source, size_type bytes,
  const std::string& path)
{
  const std::uint8_t* data = static_cast<const std::uint8_t*>(source);
  size_type done = 0;
  while(done < bytes)
  {
    size_type request = std::min(bytes - done,
      static_cast<size_type>(std::numeric_limits<ssize_t>::max()));
    ssize_t written = ::write(descriptor, data + done, request);
    if(written < 0 && errno == EINTR) { continue; }
    if(written <= 0) { throw eventError("write failed", path); }
    done += static_cast<size_type>(written);
  }
  DiskIO::write_volume += bytes;
}

void
readAll(int descriptor, void* target, size_type bytes, const std::string& path)
{
  std::uint8_t* data = static_cast<std::uint8_t*>(target);
  size_type done = 0;
  while(done < bytes)
  {
    size_type request = std::min(bytes - done,
      static_cast<size_type>(std::numeric_limits<ssize_t>::max()));
    ssize_t got = ::read(descriptor, data + done, request);
    if(got < 0 && errno == EINTR) { continue; }
    if(got <= 0) { throw eventError("short read", path); }
    done += static_cast<size_type>(got);
  }
  DiskIO::read_volume += bytes;
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
discardCache(int descriptor, off_t first, off_t bytes)
{
#if defined(POSIX_FADV_DONTNEED)
  if(bytes > 0)
  {
    static_cast<void>(::posix_fadvise(descriptor, first, bytes,
      POSIX_FADV_DONTNEED));
  }
#else
  static_cast<void>(descriptor); static_cast<void>(first); static_cast<void>(bytes);
#endif
}

size_type
checkedBytes(size_type records, size_type width, const std::string& path)
{
  if(width != 0 && records > std::numeric_limits<size_type>::max() / width)
  {
    throw eventError("record stream length overflows", path);
  }
  return records * width;
}

size_type
physicalFileBytes(const std::string& path)
{
  struct stat status;
  if(::stat(path.c_str(), &status) != 0 || status.st_size < 0)
  {
    throw eventError("cannot stat event stream", path);
  }
  if(static_cast<std::uintmax_t>(status.st_size) >
     std::numeric_limits<size_type>::max())
  {
    throw eventError("event stream is too large for this build", path);
  }
  return static_cast<size_type>(status.st_size);
}

size_type
logicalFileBytes(const std::string& path)
{
  if(!CompressedBlockReader::isFramed(path)) { return physicalFileBytes(path); }
  std::uint64_t logical = 0;
  try
  {
    logical = CompressedBlockReader::declaredLogicalSize(path);
  }
  catch(const std::runtime_error& error)
  {
    throw eventError(error.what(), path);
  }
  if(logical > std::numeric_limits<size_type>::max())
  {
    throw eventError("framed logical size overflows", path);
  }
  return static_cast<size_type>(logical);
}

void
requireFileSize(const std::string& path, size_type records, size_type width)
{
  size_type expected = checkedBytes(records, width, path);
  size_type observed = logicalFileBytes(path);
  if(observed != expected)
  {
    throw eventError("event stream has " + std::to_string(observed) +
      " bytes; expected " + std::to_string(expected), path);
  }
}

class BufferedEventWriter
{
public:
  BufferedEventWriter(const std::string& path, size_type buffer_bytes,
    MemoryBudget& budget, bool compressed) :
    path_(path), descriptor_(-1), buffer_(), used_(0), written_(0),
    cache_released_(0), closed_(false),
    reservation_(), framed_()
  {
    if(compressed)
    {
      const CompressedBlockWriter::Mode mode = CompressedBlockWriter::ZSTD;
      reservation_ = budget.reserve(CompressedBlockWriter::workingMemoryEstimate(
        buffer_bytes, mode, 1));
      framed_.reset(new CompressedBlockWriter(path, buffer_bytes, mode,
        1));
      return;
    }
    buffer_bytes = std::max(static_cast<size_type>(16), buffer_bytes);
    reservation_ = budget.reserve(buffer_bytes);
    buffer_.resize(buffer_bytes);
    descriptor_ = ::open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if(descriptor_ < 0) { throw eventError("cannot create event stream", path); }
    adviseSequential(descriptor_);
  }

  ~BufferedEventWriter()
  {
    if(!closed_)
    {
      try { close(); } catch(...) { }
    }
  }

  void byte(std::uint8_t value) { append(&value, 1); }

  void integer(std::uint64_t value)
  {
    std::uint8_t record[8]; put64(record, value); append(record, sizeof(record));
  }

  void pair(std::uint64_t first, std::uint64_t second)
  {
    std::uint8_t record[16]; put64(record, first); put64(record + 8, second);
    append(record, sizeof(record));
  }

  void close()
  {
    if(closed_) { return; }
    if(framed_)
    {
      framed_->finish(); framed_.reset(); closed_ = true;
      reservation_ = MemoryBudget::Reservation(); return;
    }
    flush();
    if(::fdatasync(descriptor_) != 0)
    {
      throw eventError("fdatasync failed", path_);
    }
    discardCache(descriptor_, 0, static_cast<off_t>(written_));
    if(::close(descriptor_) != 0)
    {
      descriptor_ = -1; throw eventError("close failed", path_);
    }
    descriptor_ = -1; closed_ = true;
    std::vector<std::uint8_t>().swap(buffer_);
    reservation_ = MemoryBudget::Reservation();
  }

private:
  std::string path_;
  int descriptor_;
  std::vector<std::uint8_t> buffer_;
  size_type used_, written_;
  off_t cache_released_;
  bool closed_;
  MemoryBudget::Reservation reservation_;
  std::unique_ptr<CompressedBlockWriter> framed_;

  void append(const void* source, size_type bytes)
  {
    if(closed_) { throw eventError("write after close", path_); }
    const std::uint8_t* data = static_cast<const std::uint8_t*>(source);
    if(framed_)
    {
      framed_->writeRecord(data, bytes);
      return;
    }
    while(bytes > 0)
    {
      size_type count = std::min(bytes, buffer_.size() - used_);
      std::memcpy(buffer_.data() + used_, data, count);
      used_ += count; data += count; bytes -= count;
      if(used_ == buffer_.size()) { flush(); }
    }
  }

  void flush()
  {
    if(used_ == 0) { return; }
    writeAll(descriptor_, buffer_.data(), used_, path_);
    written_ += used_; used_ = 0;
    if(static_cast<off_t>(written_) - cache_released_ >= CACHE_FLUSH_BYTES)
    {
      if(::fdatasync(descriptor_) != 0)
      {
        throw eventError("periodic fdatasync failed", path_);
      }
      off_t discard_end = std::max(cache_released_,
        static_cast<off_t>(written_) - CACHE_TAIL_BYTES);
      discardCache(descriptor_, cache_released_, discard_end - cache_released_);
      cache_released_ = discard_end;
    }
  }
};

class BufferedEventReader
{
public:
  BufferedEventReader(const std::string& path, size_type width,
    size_type records, size_type buffer_bytes, MemoryBudget& budget) :
    path_(path), descriptor_(-1), width_(width), records_(records), seen_(0),
    buffer_(), buffered_(0), offset_(0), consumed_bytes_(0), cache_released_(0),
    reservation_(), framed_()
  {
    if(width_ == 0) { throw eventError("zero-width event stream", path); }
    requireFileSize(path_, records_, width_);
    size_type capacity = std::max(width_, buffer_bytes - buffer_bytes % width_);
    if(CompressedBlockReader::isFramed(path_))
    {
      const size_type block = static_cast<size_type>(
        CompressedBlockReader::declaredBlockSize(path_));
      const size_type codec_bytes = CompressedBlockReader::workingMemoryEstimate(block);
      if(codec_bytes > std::numeric_limits<size_type>::max() - capacity)
      {
        throw eventError("framed reader memory estimate overflows", path_);
      }
      // An encoder owns this budget and all of its readers. No other worker
      // can release an existing reader reservation while this constructor
      // waits, so an infeasible simultaneous reader set must fail promptly.
      if(codec_bytes + capacity > budget.available())
      {
        throw eventError("memory limit cannot hold simultaneous final-event readers", path_);
      }
      reservation_ = budget.reserve(codec_bytes + capacity);
      framed_.reset(new CompressedBlockReader(path_));
    }
    else
    {
      if(capacity > budget.available())
      {
        throw eventError("memory limit cannot hold simultaneous final-event readers", path_);
      }
      reservation_ = budget.reserve(capacity);
    }
    buffer_.resize(capacity);
    if(!framed_)
    {
      descriptor_ = ::open(path.c_str(), O_RDONLY);
      if(descriptor_ < 0) { throw eventError("cannot open event stream", path); }
      adviseSequential(descriptor_);
    }
  }

  ~BufferedEventReader()
  {
    if(descriptor_ >= 0)
    {
      discardCache(descriptor_, 0, static_cast<off_t>(consumed_bytes_));
      ::close(descriptor_);
    }
  }

  bool nextByte(std::uint8_t& value)
  {
    if(width_ != 1) { throw eventError("wrong byte reader width", path_); }
    return next(&value);
  }

  bool nextInteger(std::uint64_t& value)
  {
    if(width_ != 8) { throw eventError("wrong integer reader width", path_); }
    std::uint8_t record[8];
    if(!next(record)) { return false; }
    value = get64(record); return true;
  }

  bool nextPair(std::uint64_t& first, std::uint64_t& second)
  {
    if(width_ != 16) { throw eventError("wrong pair reader width", path_); }
    std::uint8_t record[16];
    if(!next(record)) { return false; }
    first = get64(record); second = get64(record + 8); return true;
  }

  void finish()
  {
    if(seen_ != records_)
    {
      throw eventError("event stream was not fully consumed", path_);
    }
  }

private:
  std::string path_;
  int descriptor_;
  size_type width_, records_, seen_;
  std::vector<std::uint8_t> buffer_;
  size_type buffered_, offset_, consumed_bytes_;
  off_t cache_released_;
  MemoryBudget::Reservation reservation_;
  std::unique_ptr<CompressedBlockReader> framed_;

  bool next(void* target)
  {
    if(seen_ == records_) { return false; }
    if(offset_ == buffered_) { refill(); }
    std::memcpy(target, buffer_.data() + offset_, width_);
    offset_ += width_; seen_++; return true;
  }

  void refill()
  {
    size_type remaining = records_ - seen_;
    size_type bytes = std::min(buffer_.size(), checkedBytes(remaining, width_, path_));
    if(framed_)
    {
      if(framed_->read(buffer_.data(), bytes) != bytes)
      {
        throw eventError("short framed event read", path_);
      }
    }
    else { readAll(descriptor_, buffer_.data(), bytes, path_); }
    buffered_ = bytes; offset_ = 0; consumed_bytes_ += bytes;
    if(!framed_ &&
       static_cast<off_t>(consumed_bytes_) - cache_released_ >= CACHE_FLUSH_BYTES)
    {
      off_t discard_end = std::max(cache_released_,
        static_cast<off_t>(consumed_bytes_) - CACHE_TAIL_BYTES);
      discardCache(descriptor_, cache_released_, discard_end - cache_released_);
      cache_released_ = discard_end;
    }
  }
};

void
validateMetadata(const FinalEventMetadata& metadata)
{
  if(metadata.sigma == 0 || metadata.sigma > FinalEventMetadata::MAX_SIGMA)
  {
    throw eventError("invalid alphabet size in metadata");
  }
  if(metadata.fast_chars >= metadata.sigma)
  {
    throw eventError("invalid fast-character count in metadata");
  }
  size_type edge_sum = 0;
  for(size_type comp = 0; comp < metadata.sigma; comp++)
  {
    if(metadata.bwt_counts[comp] > metadata.paths)
    {
      throw eventError("BWT component count exceeds path count");
    }
    if(edge_sum > std::numeric_limits<size_type>::max() - metadata.bwt_counts[comp])
    {
      throw eventError("edge count overflows");
    }
    edge_sum += metadata.bwt_counts[comp];
  }
  if(edge_sum != metadata.total_edges)
  {
    throw eventError("component edge counts do not match total");
  }
  if(metadata.sampled_paths > metadata.paths ||
     (metadata.sampled_paths == 0) != (metadata.sample_ids == 0))
  {
    throw eventError("inconsistent sample counts");
  }
  if(metadata.sample_bits > 64 ||
     (metadata.sample_ids > 0 && metadata.sample_bits == 0) ||
     (metadata.sample_ids == 0 && metadata.sample_bits != 0))
  {
    throw eventError("invalid sample identifier width");
  }
  if(metadata.occurrence_items > metadata.paths ||
     metadata.occurrence_items > metadata.occurrence_extra)
  {
    throw eventError("inconsistent occurrence counts");
  }
  if(metadata.paths == 0 && (metadata.total_edges != 0 || metadata.redundant != 0))
  {
    throw eventError("events exist for an empty graph");
  }
}

/*
  Serialize sdsl::bit_vector_il<512> from a sequential bit source. The SDSL
  constructor interleaves one cumulative count with each eight 64-bit data
  words and retains at most 1024 cumulative samples in breadth-first midpoint
  order. Reproducing that layout here avoids its dense source and destination
  allocations while preserving the public index format byte for byte.
*/
template<class NextBit>
size_type
serializeFastVector(std::ostream& out, size_type bits, NextBit next_bit)
{
  constexpr size_type BLOCK_BITS = 512;
  constexpr size_type WORDS_PER_BLOCK = BLOCK_BITS / 64;
  if(bits > std::numeric_limits<size_type>::max() - BLOCK_BITS)
  {
    throw eventError("fast bitvector length overflows");
  }
  const size_type superblocks = (bits + BLOCK_BITS) / BLOCK_BITS;
  const size_type data_words = (bits + 64) / 64;
  if(data_words > std::numeric_limits<size_type>::max() - superblocks - 1)
  {
    throw eventError("fast bitvector size overflows");
  }
  const size_type stored_words = data_words + superblocks + 1;

  sdsl::write_member(bits, out);
  sdsl::write_member(stored_words, out);
  sdsl::write_member(superblocks, out);
  sdsl::write_member(static_cast<size_type>(9), out);
  sdsl::int_vector<64>::write_header(checkedBytes(stored_words, 64,
    "fast bitvector"), 64, out);

  size_type sample_count = 0;
  if(stored_words > 1024 * 64)
  {
    sample_count = std::min<size_type>(1024,
      static_cast<size_type>(1) << sdsl::bits::hi(superblocks));
  }
  std::vector<std::pair<size_type, size_type>> targets;
  targets.reserve(sample_count);
  std::queue<std::pair<size_type, size_type>> ranges;
  ranges.push(std::make_pair(0, superblocks));
  for(size_type order = 0; order < sample_count; order++)
  {
    const std::pair<size_type, size_type> range = ranges.front();
    ranges.pop();
    const size_type middle = range.first + (range.second - range.first) / 2;
    targets.push_back(std::make_pair(middle, order));
    ranges.push(std::make_pair(range.first, middle));
    ranges.push(std::make_pair(middle + 1, range.second));
  }
  std::sort(targets.begin(), targets.end());
  sdsl::int_vector<64> rank_samples(sample_count, 0);

  const auto write_word = [&out](std::uint64_t word)
  {
    out.write(reinterpret_cast<const char*>(&word), sizeof(word));
    if(!out) { throw eventError("cannot serialize fast bitvector"); }
  };
  size_type target = 0, cumulative = 0, position = 0;
  for(size_type word_index = 0; word_index < data_words; word_index++)
  {
    if(word_index % WORDS_PER_BLOCK == 0)
    {
      const size_type block = word_index / WORDS_PER_BLOCK;
      write_word(cumulative);
      if(target < targets.size() && targets[target].first == block)
      {
        rank_samples[targets[target].second] = cumulative;
        target++;
      }
    }
    std::uint64_t word = 0;
    for(size_type offset = 0; offset < 64 && position < bits;
        offset++, position++)
    {
      if(next_bit()) { word |= static_cast<std::uint64_t>(1) << offset; }
    }
    write_word(word); cumulative += sdsl::bits::cnt(word);
  }
  write_word(cumulative);
  if(position != bits || target != targets.size())
  {
    throw eventError("fast bitvector stream length mismatch");
  }
  rank_samples.serialize(out);
  if(!out) { throw eventError("cannot finish fast bitvector serialization"); }
  return cumulative;
}

template<class NextValue>
void
serializePackedVector(std::ostream& out, size_type values, size_type width,
  NextValue next_value)
{
  const size_type stored_width = (width == 0 ? 64 : width);
  if(stored_width > 64 ||
     values > std::numeric_limits<size_type>::max() / stored_width)
  {
    throw eventError("packed integer vector size overflows");
  }
  const size_type bits = values * stored_width;
  sdsl::int_vector<0>::write_header(bits,
    static_cast<std::uint8_t>(stored_width), out);

  const auto write_word = [&out](std::uint64_t word)
  {
    out.write(reinterpret_cast<const char*>(&word), sizeof(word));
    if(!out) { throw eventError("cannot serialize packed integer vector"); }
  };
  std::uint64_t word = 0;
  size_type used = 0, written_words = 0;
  for(size_type i = 0; i < values; i++)
  {
    const std::uint64_t value = next_value();
    if(stored_width < 64 &&
       value >= (static_cast<std::uint64_t>(1) << stored_width))
    {
      throw eventError("packed integer exceeds its declared width");
    }
    if(stored_width == 64)
    {
      write_word(value); written_words++; continue;
    }

    word |= value << used;
    const size_type next_used = used + stored_width;
    if(next_used >= 64)
    {
      write_word(word); written_words++;
      used = next_used - 64;
      word = (used == 0 ? 0 : value >> (stored_width - used));
    }
    else
    {
      used = next_used;
    }
  }
  if(used > 0) { write_word(word); written_words++; }
  const size_type expected_words = bits / 64 + (bits % 64 != 0);
  if(written_words != expected_words)
  {
    throw eventError("packed integer vector word count mismatch");
  }
}

template<class NextPosition>
void
serializePlainBitVector(std::ostream& out, size_type bits, NextPosition next_position)
{
  if(bits > std::numeric_limits<size_type>::max() - 63)
  {
    throw eventError("bitvector universe overflows");
  }
  sdsl::int_vector<1>::write_header(bits, 1, out);
  const size_type words = (bits + 63) / 64;
  size_type next = 0;
  bool available = next_position(next);
  for(size_type word_index = 0; word_index < words; word_index++)
  {
    std::uint64_t word = 0;
    const size_type first = word_index * 64;
    while(available && next < first + 64)
    {
      if(next < first) { throw eventError("bit positions are not increasing"); }
      word |= static_cast<std::uint64_t>(1) << (next - first);
      available = next_position(next);
    }
    out.write(reinterpret_cast<const char*>(&word), sizeof(word));
    if(!out) { throw eventError("cannot serialize bitvector"); }
  }
  if(available) { throw eventError("bit position exceeds bitvector universe"); }
}

template<class NextPosition>
void
serializeSelectMCL(std::ostream& out, size_type bits, size_type ones,
  NextPosition next_position)
{
  sdsl::write_member(ones, out);
  if(ones == 0) { return; }
  if(ones > bits)
  {
    throw eventError("select one-count exceeds its universe");
  }

  const size_type block_size = 4096;
  if(ones > std::numeric_limits<size_type>::max() - (block_size - 1))
  {
    throw eventError("select one-count overflows");
  }
  const size_type blocks = (ones + block_size - 1) / block_size;
  if(bits > std::numeric_limits<size_type>::max() - 63)
  {
    throw eventError("select universe overflows");
  }
  const size_type logn = sdsl::bits::hi(((bits + 63) / 64) * 64) + 1;
  const size_type logn4 = logn * logn * logn * logn;
  bool has_long_block = false;
  std::string metadata_name = TempFile::getName("gcsa_select_metadata");
  std::string payload_name = TempFile::getName("gcsa_select_payload");
  std::FILE* metadata = std::fopen(metadata_name.c_str(), "w+b");
  if(metadata == nullptr)
  {
    const std::string failed_name = metadata_name;
    TempFile::remove(metadata_name); TempFile::remove(payload_name);
    throw eventError("cannot create select metadata spool", failed_name);
  }
  std::FILE* payload = std::fopen(payload_name.c_str(), "w+b");
  if(payload == nullptr)
  {
    const std::string failed_name = payload_name;
    std::fclose(metadata);
    metadata = nullptr;
    TempFile::remove(metadata_name); TempFile::remove(payload_name);
    throw eventError("cannot create select payload spool", failed_name);
  }

  try
  {
    std::array<size_type, block_size> positions;
    size_type next = 0, seen = 0, previous = 0;
    bool have_previous = false;
    bool available = next_position(next);
    for(size_type block = 0; block < blocks; block++)
    {
      const size_type count = std::min(block_size, ones - seen);
      for(size_type i = 0; i < count; i++)
      {
        if(!available || next >= bits ||
           (have_previous && next <= previous))
        {
          throw eventError("select positions are invalid");
        }
        positions[i] = next; previous = next; have_previous = true;
        seen++; available = next_position(next);
      }
      const bool fast = (bits >= 100000);
      if(fast && count == block_size && available &&
         (next >= bits || next <= positions[count - 1]))
      {
        throw eventError("select lookahead position is invalid");
      }
      // SDSL's fast initializer samples positions 0,64,...,4032, then scans
      // 64 more arguments to choose the block width. For every nonterminal
      // full block this scan includes position 4096: the first one in the next
      // block. Although that off-by-one argument is not stored in this block,
      // reproducing it here is required for byte-identical public indexes.
      const size_type width_position =
        (fast && count == block_size && available ? next : positions[count - 1]);
      const size_type span = width_position - positions[0];
      // init_fast() encodes the final partial block as long and, unlike the
      // slow initializer, leaves its otherwise-unused superblock entry zero.
      const bool partial_fast_block = (fast && count < block_size);
      const bool is_long = (span > logn4 || partial_fast_block);
      has_long_block = has_long_block || is_long;
      // SDSL calls the second field mini_or_long, but its bit convention is
      // slightly surprising: 1 denotes a compact miniblock and 0 denotes a
      // long block. Keep this directory on disk as an explicit 9-byte record
      // per select block instead of an unbudgeted proportional vector.
      std::array<std::uint8_t, 9> block_metadata = {};
      put64(block_metadata.data(), partial_fast_block ? 0 : positions[0]);
      block_metadata[8] = static_cast<std::uint8_t>(!is_long);
      if(std::fwrite(block_metadata.data(), 1, block_metadata.size(), metadata) !=
         block_metadata.size())
      {
        throw eventError("cannot write select metadata spool");
      }
      DiskIO::write_volume += block_metadata.size();
      sdsl::int_vector<0> encoded(is_long ? block_size : 64, 0,
        is_long ? sdsl::bits::hi(partial_fast_block ? bits - 1 :
          width_position) + 1 :
        sdsl::bits::hi(span) + 1);
      if(is_long)
      {
        for(size_type i = 0; i < count; i++) { encoded[i] = positions[i]; }
      }
      else
      {
        for(size_type i = 0; i < count; i += 64)
        {
          encoded[i / 64] = positions[i] - positions[0];
        }
      }
      std::ostringstream serialized;
      encoded.serialize(serialized);
      const std::string bytes = serialized.str();
      if(std::fwrite(bytes.data(), 1, bytes.size(), payload) != bytes.size())
      {
        throw eventError("cannot write select payload spool");
      }
      DiskIO::write_volume += bytes.size();
    }
    if(available || seen != ones)
    {
      throw eventError("select position count does not match metadata");
    }

    const auto rewind_metadata = [&]()
    {
      if(std::fflush(metadata) != 0 || std::fseek(metadata, 0, SEEK_SET) != 0)
      {
        throw eventError("cannot rewind select metadata spool");
      }
    };
    const auto read_metadata = [&](std::array<std::uint8_t, 9>& record)
    {
      if(std::fread(record.data(), 1, record.size(), metadata) != record.size())
      {
        throw eventError("cannot read select metadata spool");
      }
      DiskIO::read_volume += record.size();
    };

    rewind_metadata();
    serializePackedVector(out, blocks, logn, [&]()
    {
      std::array<std::uint8_t, 9> record = {};
      read_metadata(record); return get64(record.data());
    });

    // The indicator vector is omitted entirely when every block is compact.
    if(has_long_block)
    {
      rewind_metadata();
      size_type block = 0;
      serializePlainBitVector(out, blocks, [&](size_type& mini_block)
      {
        while(block < blocks)
        {
          std::array<std::uint8_t, 9> record = {};
          read_metadata(record);
          const size_type current = block++;
          if(record[8] != 0) { mini_block = current; return true; }
        }
        return false;
      });
    }
    else { GCSA::bit_vector empty; empty.serialize(out); }

    if(std::fclose(metadata) != 0)
    {
      metadata = nullptr;
      TempFile::remove(metadata_name);
      throw eventError("cannot close select metadata spool");
    }
    metadata = nullptr;
    TempFile::remove(metadata_name);

    if(std::fflush(payload) != 0 || std::fseek(payload, 0, SEEK_SET) != 0)
    {
      throw eventError("cannot rewind select payload spool");
    }
    std::array<char, 64 * 1024> buffer;
    while(true)
    {
      size_type count = std::fread(buffer.data(), 1, buffer.size(), payload);
      if(count > 0)
      {
        DiskIO::read_volume += count;
        out.write(buffer.data(), count);
        if(!out) { throw eventError("cannot serialize select payload"); }
      }
      if(count < buffer.size())
      {
        if(std::ferror(payload)) { throw eventError("cannot read select payload spool"); }
        break;
      }
    }
    if(std::fclose(payload) != 0)
    {
      payload = nullptr;
      TempFile::remove(payload_name);
      throw eventError("cannot close select payload spool");
    }
    payload = nullptr;
    TempFile::remove(payload_name);
  }
  catch(...)
  {
    if(metadata != nullptr) { std::fclose(metadata); }
    if(payload != nullptr) { std::fclose(payload); }
    TempFile::remove(metadata_name); TempFile::remove(payload_name);
    throw;
  }
}

/*
  Serialize the default sdsl::sd_vector from a replayable monotone position
  stream. The ordinary builder retains both Elias--Fano vectors and their
  select supports. Replaying an immutable event file lets us write the exact
  same low/high layout with only a bounded input buffer and one 4096-position
  select block resident at a time.

  replay(consumer) must create a fresh source and call consumer(next), where
  next(position) yields exactly `ones` strictly increasing positions in the
  half-open universe [0, universe).
*/
template<class Replay>
void
serializeSparseVector(std::ostream& out, size_type universe, size_type ones,
  const Replay& replay)
{
  if(ones > universe)
  {
    throw eventError("sparse-vector one-count exceeds its universe");
  }
  const std::pair<size_type, size_type> params =
    SadaSparse::sd_vector::get_params(universe, ones);
  const size_type low_width = params.first, high_bits = params.second;
  if(low_width == 0 || low_width > 64 || high_bits < ones)
  {
    throw eventError("invalid sparse-vector parameters");
  }

  sdsl::write_member(universe, out);
  sdsl::write_member(static_cast<std::uint8_t>(low_width), out);
  const std::uint64_t low_mask = (low_width == 64 ?
    std::numeric_limits<std::uint64_t>::max() : sdsl::bits::lo_set[low_width]);

  replay([&](const auto& next_position)
  {
    serializePackedVector(out, ones, low_width, [&]()
    {
      size_type position = 0;
      if(!next_position(position))
      {
        throw eventError("sparse-vector position stream ended early");
      }
      return static_cast<std::uint64_t>(position) & low_mask;
    });
  });

  const auto encode_high_positions = [&](const auto& serialize)
  {
    replay([&](const auto& next_position)
    {
      size_type ordinal = 0;
      const auto next_high = [&](size_type& high_position)
      {
        size_type position = 0;
        if(!next_position(position)) { return false; }
        if(position >= universe || ordinal >= ones)
        {
          throw eventError("sparse-vector position exceeds metadata");
        }
        // get_params() currently never returns 64 for a valid size_type
        // universe, but keep the serializer free of an undefined 64-bit shift
        // if that implementation detail changes.
        const size_type high_part = (low_width == 64 ? 0 : position >> low_width);
        if(high_part > std::numeric_limits<size_type>::max() - ordinal)
        {
          throw eventError("sparse-vector high position overflows");
        }
        high_position = high_part + ordinal;
        if(high_position >= high_bits)
        {
          throw eventError("sparse-vector high position exceeds its universe");
        }
        ordinal++;
        return true;
      };
      serialize(next_high);
      if(ordinal != ones)
      {
        throw eventError("sparse-vector position count does not match metadata");
      }
    });
  };

  encode_high_positions([&](const auto& next_high)
  {
    serializePlainBitVector(out, high_bits, next_high);
  });
  encode_high_positions([&](const auto& next_high)
  {
    serializeSelectMCL(out, high_bits, ones, next_high);
  });

  // The high vector has one zero per Elias--Fano bucket. Derive those zero
  // positions from the same monotone high-one stream instead of materializing
  // the high vector for select_0 construction.
  const size_type zeros = high_bits - ones;
  encode_high_positions([&](const auto& next_high)
  {
    size_type cursor = 0, one = 0;
    bool one_available = next_high(one);
    const auto next_zero = [&](size_type& zero)
    {
      while(cursor < high_bits)
      {
        if(one_available && one < cursor)
        {
          throw eventError("sparse-vector high positions are not increasing");
        }
        if(one_available && one == cursor)
        {
          cursor++;
          one_available = next_high(one);
          continue;
        }
        zero = cursor++; return true;
      }
      if(one_available)
      {
        throw eventError("sparse-vector high position exceeds its universe");
      }
      return false;
    };
    serializeSelectMCL(out, high_bits, zeros, next_zero);
  });
  if(!out) { throw eventError("cannot serialize sparse vector"); }
}

void
validatePayloads(const FinalEventFiles& files, const FinalEventMetadata& metadata)
{
  if(files.edge_destinations.size() != metadata.sigma)
  {
    throw eventError("edge stream count does not match alphabet");
  }
  requireFileSize(files.bwt_masks, metadata.paths, 1);
  for(size_type comp = 0; comp < metadata.sigma; comp++)
  {
    requireFileSize(files.edge_destinations[comp], metadata.bwt_counts[comp], 8);
  }
  requireFileSize(files.sample_positions, metadata.sampled_paths, 8);
  requireFileSize(files.sample_ids, metadata.sample_ids, 8);
  requireFileSize(files.sample_ends, metadata.sampled_paths, 8);
  requireFileSize(files.occurrences, metadata.occurrence_items, 16);
  requireFileSize(files.redundant, metadata.redundant, 8);
}

size_type
componentBufferBytes(const ConstructionParameters& parameters, size_type readers)
{
  readers = std::max(static_cast<size_type>(1), readers);
  size_type share = parameters.getMemoryLimitBytes() / (4 * readers);
  return std::max(static_cast<size_type>(16),
    std::min(externalIOBufferSize(parameters), share));
}

class RedundancyPositionReader
{
public:
  RedundancyPositionReader(const std::string& filename, size_type paths,
    size_type redundant, const ConstructionParameters& parameters) :
    slots_(paths > 0 ? paths - 1 : 0), redundant_(redundant), slot_(0),
    cumulative_(0), budget_(parameters.getMemoryLimitBytes(),
      parameters.getMemoryLimitBytes() / 8), input_(new BufferedEventReader(filename, 8, redundant,
      componentBufferBytes(parameters, 1), budget_)),
    current_(0), available_(input_->nextInteger(current_)), finished_(false)
  {
  }

  ~RedundancyPositionReader()
  {
    if(!finished_)
    {
      try { finish(); } catch(...) { }
    }
  }

  bool next(size_type& position)
  {
    if(slot_ == slots_)
    {
      finish(); return false;
    }
    if(available_ && current_ < slot_)
    {
      throw eventError("redundancy stream is not nondecreasing");
    }
    while(available_ && current_ == slot_)
    {
      cumulative_++; available_ = input_->nextInteger(current_);
    }
    position = slot_ + cumulative_;
    slot_++; return true;
  }

private:
  size_type slots_, redundant_, slot_, cumulative_;
  MemoryBudget budget_;
  std::unique_ptr<BufferedEventReader> input_;
  std::uint64_t current_;
  bool available_, finished_;

  void finish()
  {
    if(finished_) { return; }
    if(slot_ != slots_ || available_ || cumulative_ != redundant_)
    {
      throw eventError("redundancy rank exceeds suffix-tree slot universe");
    }
    input_->finish(); finished_ = true;
  }
};

int
compareEncoded64(const void* left, const void* right)
{
  std::uint64_t a = get64(static_cast<const std::uint8_t*>(left));
  std::uint64_t b = get64(static_cast<const std::uint8_t*>(right));
  return (a < b ? -1 : (a > b ? 1 : 0));
}

} // namespace

//------------------------------------------------------------------------------

struct SpillableNodeSet::Impl
{
  size_type budget, fan_in, collection_capacity, count;
  std::vector<node_type> values, read_buffer;
  std::string raw_name, sorted_name;
  int raw_descriptor, read_descriptor;
  size_type raw_bytes, read_values, read_offset, values_seen;
  off_t raw_cache_released, read_cache_released;
  bool finished, on_disk;
  MemoryBudget::Reservation reservation;

  Impl(size_type byte_budget, size_type merge_fan_in, MemoryBudget& memory) :
    budget(byte_budget), fan_in(std::max(static_cast<size_type>(2), merge_fan_in)),
    collection_capacity(0), count(0), values(), read_buffer(), raw_name(),
    sorted_name(), raw_descriptor(-1), read_descriptor(-1), raw_bytes(0),
    read_values(0), read_offset(0), values_seen(0), raw_cache_released(0),
    read_cache_released(0), finished(false), on_disk(false), reservation()
  {
    if(this->budget < SpillableNodeSet::minimumBudget())
    {
      throw eventError("start-node set budget is too small");
    }
    this->reservation = memory.reserve(this->budget);

    // Leave half the reservation unused during collection. It covers vector
    // allocator slack and the small file/stream objects used when the vector
    // is flushed. The external sorter receives the entire reservation only
    // after this vector has released its storage.
    this->collection_capacity = std::max(static_cast<size_type>(1),
      this->budget / (2 * sizeof(node_type)));
    // collection_capacity is a spill threshold, not a size hint. Let the
    // vector grow to actual demand while retaining the full admission
    // reservation for a possible sorter phase.
  }

  ~Impl()
  {
    this->closeRaw(false);
    this->closeRead();
    if(!this->raw_name.empty()) { TempFile::remove(this->raw_name); }
    if(!this->sorted_name.empty()) { TempFile::remove(this->sorted_name); }
  }

  void closeRaw(bool complete)
  {
    if(this->raw_descriptor < 0) { return; }
    if(complete)
    {
      if(::fdatasync(this->raw_descriptor) != 0)
      {
        throw eventError("cannot sync spilled start-node set", this->raw_name);
      }
      discardCache(this->raw_descriptor, this->raw_cache_released,
        static_cast<off_t>(this->raw_bytes) - this->raw_cache_released);
    }
    if(::close(this->raw_descriptor) != 0 && complete)
    {
      this->raw_descriptor = -1;
      throw eventError("cannot close spilled start-node set", this->raw_name);
    }
    this->raw_descriptor = -1;
  }

  void closeRead()
  {
    if(this->read_descriptor >= 0)
    {
      discardCache(this->read_descriptor, this->read_cache_released,
        static_cast<off_t>(this->values_seen * sizeof(node_type)) -
        this->read_cache_released);
      ::close(this->read_descriptor); this->read_descriptor = -1;
    }
    this->read_values = this->read_offset = this->values_seen = 0;
    this->read_cache_released = 0;
  }

  void openRaw()
  {
    if(this->raw_descriptor >= 0) { return; }
    this->raw_name = TempFile::getName("gcsa_final_from_nodes_raw");
    this->raw_descriptor = ::open(this->raw_name.c_str(),
      O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if(this->raw_descriptor < 0)
    {
      throw eventError("cannot create spilled start-node set", this->raw_name);
    }
    adviseSequential(this->raw_descriptor);
  }

  void flushValues()
  {
    if(this->values.empty()) { return; }
    this->openRaw();
    size_type bytes = checkedBytes(this->values.size(), sizeof(node_type),
      this->raw_name);
    writeAll(this->raw_descriptor, this->values.data(), bytes, this->raw_name);
    this->raw_bytes += bytes; this->values.clear();
    if(static_cast<off_t>(this->raw_bytes) - this->raw_cache_released >=
       CACHE_FLUSH_BYTES)
    {
      if(::fdatasync(this->raw_descriptor) != 0)
      {
        throw eventError("cannot periodically sync spilled start-node set",
          this->raw_name);
      }
      off_t discard_end = std::max(this->raw_cache_released,
        static_cast<off_t>(this->raw_bytes) - CACHE_TAIL_BYTES);
      discardCache(this->raw_descriptor, this->raw_cache_released,
        discard_end - this->raw_cache_released);
      this->raw_cache_released = discard_end;
    }
  }

  void reset()
  {
    this->closeRaw(false); this->closeRead();
    if(!this->raw_name.empty()) { TempFile::remove(this->raw_name); this->raw_name.clear(); }
    if(!this->sorted_name.empty()) { TempFile::remove(this->sorted_name); this->sorted_name.clear(); }
    std::vector<node_type>().swap(this->read_buffer);
    // Reuse whatever capacity the previous set grew into; do not re-reserve to
    // the spill threshold for the same reason the constructor does not.
    this->values.clear();
    this->count = this->raw_bytes = 0;
    this->raw_cache_released = 0;
    this->finished = this->on_disk = false;
  }

  void push(node_type value)
  {
    if(this->finished) { throw eventError("append to completed start-node set"); }
    if(this->values.size() == this->collection_capacity) { this->flushValues(); }
    this->values.push_back(value);
  }

  void complete()
  {
    if(this->finished) { throw eventError("start-node set completed twice"); }
    if(this->raw_descriptor < 0)
    {
      sequentialSort(this->values.begin(), this->values.end());
      this->values.resize(std::unique(this->values.begin(), this->values.end()) -
        this->values.begin());
      this->count = this->values.size(); this->finished = true; return;
    }

    this->flushValues(); this->closeRaw(true);
    std::vector<node_type>().swap(this->values);
    this->sorted_name = TempFile::getName("gcsa_final_from_nodes_sorted");
    size_type unique_values = 0;
    try
    {
      ExternalFixedRecordSorter::sortAndReduce(this->raw_name,
        this->sorted_name, sizeof(node_type), this->budget, this->fan_in,
        [](const void* left, const void* right) {
          node_type a, b;
          std::memcpy(&a, left, sizeof(a)); std::memcpy(&b, right, sizeof(b));
          return (a < b ? -1 : (a > b ? 1 : 0));
        },
        [&unique_values](const void* record, bool first, bool, std::ostream& output) {
          if(first)
          {
            output.write(static_cast<const char*>(record), sizeof(node_type));
            unique_values++;
          }
        }, nullptr);
    }
    catch(...)
    {
      TempFile::remove(this->sorted_name); this->sorted_name.clear(); throw;
    }
    TempFile::remove(this->raw_name); this->raw_name.clear();
    this->count = unique_values; this->on_disk = true; this->finished = true;

    size_type read_bytes = std::max(sizeof(node_type),
      std::min(static_cast<size_type>(64 * KILOBYTE), this->budget / 4));
    this->read_buffer.reserve(std::max(static_cast<size_type>(1),
      read_bytes / sizeof(node_type)));
  }

  void startRead()
  {
    if(!this->finished) { throw eventError("read from incomplete start-node set"); }
    if(!this->on_disk) { this->read_offset = 0; return; }
    this->closeRead();
    this->read_descriptor = ::open(this->sorted_name.c_str(), O_RDONLY);
    if(this->read_descriptor < 0)
    {
      throw eventError("cannot open spilled start-node set", this->sorted_name);
    }
    adviseSequential(this->read_descriptor);
  }

  bool readNext(node_type& value)
  {
    if(!this->finished) { throw eventError("read from incomplete start-node set"); }
    if(!this->on_disk)
    {
      if(this->read_offset == this->values.size()) { return false; }
      value = this->values[this->read_offset++]; return true;
    }
    if(this->values_seen == this->count) { return false; }
    if(this->read_offset == this->read_values)
    {
      size_type records = std::min(static_cast<size_type>(this->read_buffer.capacity()),
        this->count - this->values_seen);
      this->read_buffer.resize(records);
      readAll(this->read_descriptor, this->read_buffer.data(),
        checkedBytes(records, sizeof(node_type), this->sorted_name),
        this->sorted_name);
      this->read_values = records; this->read_offset = 0;
    }
    value = this->read_buffer[this->read_offset++]; this->values_seen++;
    if(static_cast<off_t>(this->values_seen * sizeof(node_type)) -
       this->read_cache_released >= CACHE_FLUSH_BYTES)
    {
      off_t discard_end = std::max(this->read_cache_released,
        static_cast<off_t>(this->values_seen * sizeof(node_type)) - CACHE_TAIL_BYTES);
      discardCache(this->read_descriptor, this->read_cache_released,
        discard_end - this->read_cache_released);
      this->read_cache_released = discard_end;
    }
    return true;
  }
};

SpillableNodeSet::SpillableNodeSet(size_type byte_budget,
  size_type merge_fan_in, MemoryBudget& memory) :
  impl_(new Impl(byte_budget, merge_fan_in, memory))
{
}

SpillableNodeSet::~SpillableNodeSet() = default;

void SpillableNodeSet::clear() { this->impl_->reset(); }
void SpillableNodeSet::push_back(node_type value) { this->impl_->push(value); }
void SpillableNodeSet::finish() { this->impl_->complete(); }
size_type SpillableNodeSet::size() const { return this->impl_->count; }
bool SpillableNodeSet::spilled() const { return this->impl_->on_disk; }
void SpillableNodeSet::rewind() { this->impl_->startRead(); }
bool SpillableNodeSet::next(node_type& value) { return this->impl_->readNext(value); }

size_type
SpillableNodeSet::minimumBudget()
{
  return ExternalFixedRecordSorter::minimumBudget(sizeof(node_type));
}

FinalEventMetadata::FinalEventMetadata() :
  paths(0), sigma(0), fast_chars(0), total_edges(0), sampled_paths(0),
  sample_ids(0), sample_bits(0), occurrence_items(0), occurrence_extra(0),
  redundant(0), bwt_counts()
{
  bwt_counts.fill(0);
}

void
serializeFastBWTComponent(std::ostream& out, const std::string& mask_file,
  size_type paths, size_type expected_ones, comp_type comp,
  const ConstructionParameters& parameters)
{
  if(comp >= FinalEventMetadata::MAX_SIGMA)
  {
    throw eventError("fast BWT component is outside the mask width");
  }
  MemoryBudget budget(parameters.getMemoryLimitBytes(),
    parameters.getMemoryLimitBytes() / 8);
  BufferedEventReader masks(mask_file, 1, paths,
    componentBufferBytes(parameters, 1), budget);
  size_type observed = serializeFastVector(out, paths, [&]()
  {
    std::uint8_t mask = 0;
    masks.nextByte(mask);
    return (mask & (static_cast<size_type>(1) << comp)) != 0;
  });
  masks.finish();
  if(observed != expected_ones)
  {
    throw eventError("BWT mask count does not match metadata");
  }
}

void
serializeSparseBWTComponent(std::ostream& out, const std::string& mask_file,
  size_type paths, size_type expected_ones, comp_type comp,
  const ConstructionParameters& parameters)
{
  if(comp >= FinalEventMetadata::MAX_SIGMA || expected_ones > paths)
  {
    throw eventError("invalid sparse BWT component metadata");
  }

  const auto replay = [&](const auto& consume)
  {
    MemoryBudget budget(parameters.getMemoryLimitBytes(),
      parameters.getMemoryLimitBytes() / 8);
    BufferedEventReader masks(mask_file, 1, paths,
      componentBufferBytes(parameters, 1), budget);
    size_type path = 0, observed = 0;
    const auto next = [&](size_type& position)
    {
      while(path < paths)
      {
        std::uint8_t mask = 0;
        masks.nextByte(mask);
        const size_type current = path++;
        if(mask & (static_cast<size_type>(1) << comp))
        {
          position = current; observed++; return true;
        }
      }
      return false;
    };

    consume(next);
    // Packed low parts request exactly expected_ones positions instead of
    // reading to EOF. Drain the remaining mask bytes so extra set bits and
    // truncated streams are still detected in every replay.
    size_type ignored = 0;
    while(next(ignored)) { }
    masks.finish();
    if(observed != expected_ones)
    {
      throw eventError("BWT mask count does not match metadata");
    }
  };

  serializeSparseVector(out, paths, expected_ones, replay);
}

void
serializeSampleIds(std::ostream& out, const std::string& sample_file,
  size_type samples, size_type sample_bits,
  const ConstructionParameters& parameters)
{
  if((samples == 0) != (sample_bits == 0) || sample_bits > 64)
  {
    throw eventError("invalid stored-sample width");
  }
  MemoryBudget budget(parameters.getMemoryLimitBytes(),
    parameters.getMemoryLimitBytes() / 8);
  BufferedEventReader ids(sample_file, 8, samples,
    componentBufferBytes(parameters, 1), budget);
  serializePackedVector(out, samples, sample_bits, [&]()
  {
    std::uint64_t node = 0;
    ids.nextInteger(node); return node;
  });
  ids.finish();
}

void
serializeSampleBoundaries(std::ostream& out,
  const std::string& sample_end_file, size_type sample_ids,
  size_type sampled_paths, const ConstructionParameters& parameters)
{
  if(sampled_paths > sample_ids ||
     ((sampled_paths == 0) != (sample_ids == 0)))
  {
    throw eventError("inconsistent sample-boundary counts");
  }

  // SDSL serializes the supported bitvector before select_support_mcl. Replay
  // the immutable monotone event stream for each component, retaining only one
  // configured input buffer instead of a bit for every stored sample ID.
  const auto pass = [&](const auto& serialize)
  {
    MemoryBudget budget(parameters.getMemoryLimitBytes(),
      parameters.getMemoryLimitBytes() / 8);
    BufferedEventReader ends(sample_end_file, 8, sampled_paths,
      componentBufferBytes(parameters, 1), budget);
    size_type consumed = 0, previous = 0;
    const auto next = [&](size_type& position)
    {
      if(consumed == sampled_paths) { return false; }
      std::uint64_t end = 0;
      ends.nextInteger(end);
      if(end >= sample_ids || (consumed > 0 && end <= previous))
      {
        throw eventError("sample endpoint stream is not strictly increasing");
      }
      position = static_cast<size_type>(end);
      previous = position; consumed++;
      return true;
    };
    serialize(next);
    ends.finish();
    if(consumed != sampled_paths ||
       (sampled_paths > 0 && previous + 1 != sample_ids))
    {
      throw eventError("final sample endpoint does not cover all sample IDs");
    }
  };

  pass([&](const auto& next)
  {
    serializePlainBitVector(out, sample_ids, next);
  });
  pass([&](const auto& next)
  {
    serializeSelectMCL(out, sample_ids, sampled_paths, next);
  });
}

void
serializeOccurrencePointers(std::ostream& out,
  const std::string& occurrence_file, size_type paths,
  size_type occurrence_items, size_type occurrence_extra,
  const ConstructionParameters& parameters)
{
  if(occurrence_items > paths || occurrence_items > occurrence_extra ||
     (paths == 0 && (occurrence_items != 0 || occurrence_extra != 0)))
  {
    throw eventError("inconsistent occurrence counts");
  }

  // Every sparse-vector component replays the same immutable event stream.
  // This deliberately trades sequential reads for a construction whose RSS
  // does not grow with the number of paths or nonzero occurrence values.
  const auto replay = [&](bool cumulative_values, const auto& consume)
  {
    MemoryBudget budget(parameters.getMemoryLimitBytes(),
      parameters.getMemoryLimitBytes() / 8);
    BufferedEventReader input(occurrence_file, 16, occurrence_items,
      componentBufferBytes(parameters, 1), budget);
    size_type seen = 0, tail = 0;
    std::uint64_t previous = 0;
    const auto next = [&](size_type& position)
    {
      if(seen == occurrence_items) { return false; }
      std::uint64_t path = 0, extra = 0;
      input.nextPair(path, extra);
      if(path >= paths || extra == 0 || (seen > 0 && path <= previous) ||
         tail > occurrence_extra || extra > occurrence_extra - tail)
      {
        throw eventError("invalid occurrence event");
      }
      tail += static_cast<size_type>(extra);
      position = (cumulative_values ? tail - 1 : static_cast<size_type>(path));
      previous = path; seen++;
      return true;
    };
    consume(next);
    input.finish();
    if(seen != occurrence_items || tail != occurrence_extra)
    {
      throw eventError("occurrence value total mismatch");
    }
  };

  serializeSparseVector(out, paths, occurrence_items, [&](const auto& consume)
  {
    replay(false, consume);
  });
  SadaSparse::sd_vector::rank_1_type empty_rank;
  empty_rank.serialize(out);

  serializeSparseVector(out, occurrence_extra, occurrence_items,
    [&](const auto& consume)
    {
      replay(true, consume);
    });
  SadaSparse::sd_vector::select_1_type empty_select;
  empty_select.serialize(out);
  if(!out) { throw eventError("cannot serialize occurrence pointers"); }
}

void
serializeRedundantPointers(std::ostream& out, const std::string& redundancy_file,
  size_type paths, size_type redundant, const ConstructionParameters& parameters)
{
  if(paths == 0 && redundant != 0)
  {
    throw eventError("redundancy events exist for an empty graph");
  }
  const size_type slots = (paths > 0 ? paths - 1 : 0);
  if(redundant > std::numeric_limits<size_type>::max() - slots)
  {
    throw eventError("redundancy bitvector length overflows");
  }
  const size_type bits = slots + redundant;
  {
    RedundancyPositionReader input(redundancy_file, paths, redundant, parameters);
    serializePlainBitVector(out, bits, [&](size_type& position)
    {
      return input.next(position);
    });
  }
  {
    RedundancyPositionReader input(redundancy_file, paths, redundant, parameters);
    serializeSelectMCL(out, bits, slots, [&](size_type& position)
    {
      return input.next(position);
    });
  }
}

FinalEventFiles::FinalEventFiles(size_type sigma) :
  bwt_masks(TempFile::getName("gcsa_final_bwt_masks")),
  edge_destinations(),
  sample_positions(TempFile::getName("gcsa_final_sample_paths")),
  sample_ids(TempFile::getName("gcsa_final_sample_ids")),
  sample_ends(TempFile::getName("gcsa_final_sample_ends")),
  occurrences(TempFile::getName("gcsa_final_occurrences")),
  redundant(TempFile::getName("gcsa_final_redundant"))
{
  edge_destinations.reserve(sigma);
  for(size_type comp = 0; comp < sigma; comp++)
  {
    edge_destinations.push_back(TempFile::getName("gcsa_final_edges"));
  }
}

FinalEventFiles::~FinalEventFiles()
{
  this->clear();
}

void
FinalEventFiles::clear()
{
  TempFile::remove(this->bwt_masks);
  for(std::string& path : this->edge_destinations) { TempFile::remove(path); }
  TempFile::remove(this->sample_positions); TempFile::remove(this->sample_ids);
  TempFile::remove(this->sample_ends); TempFile::remove(this->occurrences);
  TempFile::remove(this->redundant);
  this->bwt_masks.clear(); this->edge_destinations.clear();
  this->sample_positions.clear(); this->sample_ids.clear(); this->sample_ends.clear();
  this->occurrences.clear(); this->redundant.clear();
}

struct FinalEventWriter::Impl
{
  size_type sigma;
  std::unique_ptr<BufferedEventWriter> masks;
  std::vector<std::unique_ptr<BufferedEventWriter>> edges;
  std::unique_ptr<BufferedEventWriter> sample_positions, sample_ids, sample_ends;
  std::unique_ptr<BufferedEventWriter> occurrences, redundant;
  FinalEventMetadata metadata;
  bool finished;

  Impl(const FinalEventFiles& files, size_type alphabet_size,
    size_type buffer_bytes, MemoryBudget& budget) :
    sigma(alphabet_size), masks(), edges(), sample_positions(), sample_ids(),
    sample_ends(), occurrences(), redundant(), metadata(), finished(false)
  {
    if(sigma == 0 || sigma > FinalEventMetadata::MAX_SIGMA ||
       files.edge_destinations.size() != sigma)
    {
      throw eventError("invalid writer alphabet size");
    }
    metadata.sigma = sigma;
    const size_type streams = sigma + 6;
    size_type block_size = std::max(static_cast<size_type>(16),
      std::min(static_cast<size_type>(MEGABYTE), buffer_bytes));
    bool compressed = true;
    const auto fits = [&](size_type block) -> bool
    {
      const size_type one = CompressedBlockWriter::workingMemoryEstimate(
        block, CompressedBlockWriter::ZSTD, 1);
      // Preserve the final scan's existing half-budget frontier for its two
      // mutable disk arrays after all simultaneous codec workspaces enter.
      return one <= budget.available() / (2 * streams);
    };
    while(block_size > 16 && !fits(block_size))
    {
      block_size = std::max(static_cast<size_type>(16), block_size / 2);
    }
    if(!fits(block_size)) { compressed = false; }
    masks.reset(new BufferedEventWriter(files.bwt_masks, compressed ?
      block_size : buffer_bytes, budget, compressed));
    edges.reserve(sigma);
    for(size_type comp = 0; comp < sigma; comp++)
    {
      edges.emplace_back(new BufferedEventWriter(files.edge_destinations[comp],
        compressed ? block_size : buffer_bytes, budget, compressed));
    }
    sample_positions.reset(new BufferedEventWriter(files.sample_positions,
      compressed ? block_size : buffer_bytes, budget, compressed));
    sample_ids.reset(new BufferedEventWriter(files.sample_ids,
      compressed ? block_size : buffer_bytes, budget, compressed));
    sample_ends.reset(new BufferedEventWriter(files.sample_ends,
      compressed ? block_size : buffer_bytes, budget, compressed));
    occurrences.reset(new BufferedEventWriter(files.occurrences,
      compressed ? block_size : buffer_bytes, budget, compressed));
    // The external sorter consumes this stream before publication. Keep it raw
    // so the sorted output remains compatible with its fixed-record interface.
    redundant.reset(new BufferedEventWriter(files.redundant,
      buffer_bytes, budget, false));
  }

  void closeAll()
  {
    masks->close();
    for(auto& writer : edges) { writer->close(); }
    sample_positions->close(); sample_ids->close(); sample_ends->close();
    occurrences->close(); redundant->close();
  }
};

FinalEventWriter::FinalEventWriter(const FinalEventFiles& files, size_type sigma,
  size_type buffer_bytes, MemoryBudget& budget) :
  impl_(new Impl(files, sigma, buffer_bytes, budget))
{
}

FinalEventWriter::~FinalEventWriter()
{
}

void
FinalEventWriter::path(byte_type predecessor_mask)
{
  if(this->impl_->finished) { throw eventError("path event after finish"); }
  if(this->impl_->sigma < 8 &&
     (predecessor_mask >> this->impl_->sigma) != 0)
  {
    throw eventError("BWT mask has a component outside the alphabet");
  }
  this->impl_->masks->byte(predecessor_mask); this->impl_->metadata.paths++;
}

void
FinalEventWriter::edge(comp_type comp, size_type source_path)
{
  if(this->impl_->finished || comp >= this->impl_->sigma)
  {
    throw eventError("invalid edge event");
  }
  this->impl_->edges[comp]->integer(source_path);
  this->impl_->metadata.bwt_counts[comp]++;
  this->impl_->metadata.total_edges++;
}

void
FinalEventWriter::sampledPath(size_type path)
{
  this->impl_->sample_positions->integer(path);
  this->impl_->metadata.sampled_paths++;
}

void
FinalEventWriter::sample(node_type node)
{
  this->impl_->sample_ids->integer(node);
  this->impl_->metadata.sample_ids++;
  this->impl_->metadata.sample_bits = std::max(this->impl_->metadata.sample_bits,
    bit_length(node));
}

void
FinalEventWriter::sampleEnd()
{
  if(this->impl_->metadata.sample_ids == 0)
  {
    throw eventError("sample end without a sample ID");
  }
  this->impl_->sample_ends->integer(this->impl_->metadata.sample_ids - 1);
}

void
FinalEventWriter::occurrence(size_type path, size_type extra)
{
  if(extra == 0) { return; }
  this->impl_->occurrences->pair(path, extra);
  this->impl_->metadata.occurrence_items++;
  if(this->impl_->metadata.occurrence_extra >
     std::numeric_limits<size_type>::max() - extra)
  {
    throw eventError("occurrence total overflows");
  }
  this->impl_->metadata.occurrence_extra += extra;
}

void
FinalEventWriter::redundancy(size_type path)
{
  this->impl_->redundant->integer(path); this->impl_->metadata.redundant++;
}

FinalEventMetadata
FinalEventWriter::finish()
{
  if(this->impl_->finished) { throw eventError("event writer already finished"); }
  this->impl_->closeAll(); this->impl_->finished = true;
  return this->impl_->metadata;
}

void
sortFinalRedundancy(FinalEventFiles& files,
  const ConstructionParameters& parameters)
{
  std::string sorted = TempFile::getName("gcsa_final_redundant_sorted");
  try
  {
    size_type usable_memory = parameters.getMemoryLimitBytes() -
      parameters.getMemoryLimitBytes() / 8;
    size_type budget = std::min(usable_memory, externalSortRunSize(parameters));
    const size_type sorter_descriptors = 4;
    if(externalMaxOpenFiles() < sorter_descriptors + 2 * 2)
    {
      throw eventError("max-open-files cannot hold redundancy sort streams");
    }
    size_type fan_in = std::min(externalMergeFanIn(),
      (externalMaxOpenFiles() - sorter_descriptors) / 2);
    ExternalFixedRecordSorter::sort(files.redundant, sorted, 8, budget,
      fan_in, compareEncoded64);
    TempFile::remove(files.redundant); files.redundant = sorted;
  }
  catch(...)
  {
    TempFile::remove(sorted); throw;
  }
}

void
storeFinalComponents(const GCSAHeader& header,
  const Alphabet& source_alphabet, const FinalEventFiles& files,
  const FinalEventMetadata& metadata,
  const ConstructionParameters& parameters, const std::string& filename)
{
  validateMetadata(metadata); validatePayloads(files, metadata);
  if(!header.check() || header.path_nodes != metadata.paths ||
     header.edges != metadata.total_edges || source_alphabet.sigma != metadata.sigma ||
     source_alphabet.fast_chars != metadata.fast_chars)
  {
    throw eventError("header or alphabet does not match final events");
  }

  std::string partial = filename + "." +
    std::to_string(static_cast<unsigned long long>(::getpid())) + ".partial";
  std::ofstream out(partial.c_str(), std::ios_base::binary | std::ios_base::trunc);
  if(!out) { throw eventError("cannot open partial final index", partial); }
  const auto write = [&](const auto& value)
  {
    value.serialize(out);
    if(!out) { throw eventError("cannot write partial final index", partial); }
  };
  try
  {
    sdsl::int_vector<64> counts(metadata.sigma, 0);
    for(size_type comp = 0; comp < metadata.sigma; comp++)
    {
      counts[comp] = metadata.bwt_counts[comp];
    }
    Alphabet alphabet(counts, source_alphabet.char2comp, source_alphabet.comp2char);
    write(header);
    write(alphabet);

    const size_type component_tasks = 2 * metadata.sigma + 6;
    const auto encode = [&](size_type task, std::ostream& out,
      const ConstructionParameters& parameters)
    {
      const auto write = [&](const auto& value)
      {
        value.serialize(out);
        if(!out) { throw eventError("cannot encode final component", partial); }
      };
      if(task < metadata.sigma)
      {
        const size_type comp = task;
        if(comp > 0 && comp <= metadata.fast_chars)
        {
          serializeFastBWTComponent(out, files.bwt_masks, metadata.paths,
            metadata.bwt_counts[comp], static_cast<comp_type>(comp), parameters);
        }
        else { GCSA::fast_vector empty; write(empty); }
      }
      else if(task == metadata.sigma)
      {
        for(size_type comp = 0; comp < metadata.sigma; comp++)
        {
          GCSA::fast_vector::rank_1_type empty; write(empty);
        }
      }
      else if(task <= 2 * metadata.sigma)
      {
        const size_type comp = task - metadata.sigma - 1;
        if(comp > 0 && comp <= metadata.fast_chars)
        {
          GCSA::sparse_vector empty; write(empty);
        }
        else
        {
          serializeSparseBWTComponent(out, files.bwt_masks, metadata.paths,
            metadata.bwt_counts[comp], static_cast<comp_type>(comp), parameters);
        }
      }
      else if(task == 2 * metadata.sigma + 1)
      {
        for(size_type comp = 0; comp < metadata.sigma; comp++)
        {
          GCSA::sparse_vector::rank_1_type empty; write(empty);
        }
      }
      else if(task == 2 * metadata.sigma + 2)
      {
        MemoryBudget budget(parameters.getMemoryLimitBytes(),
          parameters.getMemoryLimitBytes() / 8);
        size_type buffer_bytes = componentBufferBytes(parameters, metadata.sigma);
        std::vector<std::unique_ptr<BufferedEventReader>> readers;
        std::vector<std::uint64_t> current(metadata.sigma);
        std::vector<bool> available(metadata.sigma, false);
        for(size_type comp = 0; comp < metadata.sigma; comp++)
        {
          readers.emplace_back(new BufferedEventReader(files.edge_destinations[comp],
            8, metadata.bwt_counts[comp], buffer_bytes, budget));
          available[comp] = readers.back()->nextInteger(current[comp]);
        }
        size_type path = 0, remaining_degree = 0, consumed_edges = 0;
        size_type observed_paths = serializeFastVector(out, metadata.total_edges, [&]()
        {
          if(remaining_degree == 0)
          {
            if(path >= metadata.paths)
            {
              throw eventError("edge stream exceeds path universe");
            }
            size_type degree = 0;
            for(size_type comp = 0; comp < metadata.sigma; comp++)
            {
              if(available[comp] && current[comp] < path)
              {
                throw eventError("edge stream is not nondecreasing");
              }
              while(available[comp] && current[comp] == path)
              {
                degree++;
                available[comp] = readers[comp]->nextInteger(current[comp]);
              }
            }
            if(degree == 0 || degree > metadata.total_edges ||
               consumed_edges > metadata.total_edges - degree)
            {
              throw eventError("invalid edge stream");
            }
            consumed_edges += degree;
            remaining_degree = degree;
            path++;
          }
          remaining_degree--;
          return (remaining_degree == 0);
        });
        for(size_type comp = 0; comp < metadata.sigma; comp++)
        {
          if(available[comp])
          {
            throw eventError("edge rank exceeds path universe");
          }
          readers[comp]->finish();
        }
        if(path != metadata.paths || remaining_degree != 0 ||
           consumed_edges != metadata.total_edges || observed_paths != metadata.paths)
        {
          throw eventError("edge total mismatch");
        }
        GCSA::fast_vector::rank_1_type empty_rank;
        write(empty_rank);
      }

      else if(task == 2 * metadata.sigma + 3)
      {
        MemoryBudget budget(parameters.getMemoryLimitBytes(),
          parameters.getMemoryLimitBytes() / 8);
        size_type buffer_bytes = componentBufferBytes(parameters, 1);
        {
          BufferedEventReader positions(files.sample_positions, 8,
            metadata.sampled_paths, buffer_bytes, budget);
          std::uint64_t next_sample = 0, previous = 0;
          bool sample_available = positions.nextInteger(next_sample);
          size_type path = 0;
          size_type observed_samples = serializeFastVector(out, metadata.paths, [&]()
          {
            if(sample_available && next_sample < path)
            {
              throw eventError("sampled path stream is not strictly increasing");
            }
            bool result = (sample_available && next_sample == path);
            if(result)
            {
              previous = next_sample;
              sample_available = positions.nextInteger(next_sample);
              if(sample_available && next_sample <= previous)
              {
                throw eventError("sampled path stream is not strictly increasing");
              }
            }
            path++;
            return result;
          });
          positions.finish();
          if(sample_available || observed_samples != metadata.sampled_paths)
          {
            throw eventError("sampled path stream does not match metadata");
          }
        }
        GCSA::fast_vector::rank_1_type empty_rank;
        write(empty_rank);

        serializeSampleIds(out, files.sample_ids, metadata.sample_ids,
          metadata.sample_bits, parameters);
        if(!out) { throw eventError("cannot write partial final index", partial); }
        serializeSampleBoundaries(out, files.sample_ends, metadata.sample_ids,
          metadata.sampled_paths, parameters);
        if(!out) { throw eventError("cannot write partial final index", partial); }
      }

      else if(task == 2 * metadata.sigma + 4)
      {
        serializeOccurrencePointers(out, files.occurrences, metadata.paths,
          metadata.occurrence_items, metadata.occurrence_extra, parameters);
        if(!out) { throw eventError("cannot write partial final index", partial); }
      }
      else
      {
        if(metadata.paths == 0 && metadata.redundant != 0)
        {
          throw eventError("redundancy events exist for an empty graph");
        }
        size_type slots = (metadata.paths > 0 ? metadata.paths - 1 : 0);
        if(metadata.redundant > std::numeric_limits<size_type>::max() - slots)
        {
          throw eventError("redundancy bitvector length overflows");
        }
        serializeRedundantPointers(out, files.redundant, metadata.paths,
          metadata.redundant, parameters);
        if(!out) { throw eventError("cannot write partial final index", partial); }
      }
      if(!out) { throw eventError("cannot encode final component", partial); }
    };
    for(size_type task = 0; task < component_tasks; task++)
    {
      encode(task, out, parameters);
    }

    out.flush();
    out.close();
    if(!out) { throw eventError("cannot finish partial final index", partial); }
    int descriptor = ::open(partial.c_str(), O_RDONLY);
    if(descriptor < 0 || ::fdatasync(descriptor) != 0)
    {
      if(descriptor >= 0) { ::close(descriptor); }
      throw eventError("cannot sync partial final index", partial);
    }
    if(::close(descriptor) != 0 ||
       ::rename(partial.c_str(), filename.c_str()) != 0)
    {
      throw eventError("cannot publish final index", filename);
    }
    std::filesystem::path parent = std::filesystem::path(filename).parent_path();
    if(parent.empty()) { parent = "."; }
    descriptor = ::open(parent.c_str(), O_RDONLY | O_DIRECTORY);
    if(descriptor < 0 || ::fsync(descriptor) != 0)
    {
      if(descriptor >= 0) { ::close(descriptor); }
      throw eventError("cannot sync final index directory", filename);
    }
    if(::close(descriptor) != 0)
    {
      throw eventError("cannot close final index directory", filename);
    }

  }
  catch(...)
  {
    out.close();
    ::unlink(partial.c_str());
    throw;
  }
}

} // namespace gcsa
