/*
  Copyright (c) 2026 GCSA2 contributors

  External event streams for final GCSA construction.
*/

#include <gcsa/final_events.h>

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

constexpr char EVENT_MAGIC[8] = { 'G', 'C', 'S', 'A', 'E', 'V', '1', '\0' };
constexpr std::uint64_t EVENT_VERSION = 1;
constexpr size_type EVENT_METADATA_WORDS = 11 + FinalEventMetadata::MAX_SIGMA;
constexpr off_t CACHE_FLUSH_BYTES = 512 * 1024 * 1024;
constexpr off_t CACHE_TAIL_BYTES = 64 * 1024 * 1024;

const std::string FINAL_TASK = "final";
const std::string FINAL_PHASE = "events";

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
rawFileBytes(const std::string& path)
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

void
requireFileSize(const std::string& path, size_type records, size_type width)
{
  size_type expected = checkedBytes(records, width, path);
  size_type observed = rawFileBytes(path);
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
    MemoryBudget& budget, const std::string& task) :
    path_(path), descriptor_(-1), buffer_(), used_(0), written_(0),
    cache_released_(0), closed_(false), reservation_()
  {
    buffer_bytes = std::max(static_cast<size_type>(16), buffer_bytes);
    reservation_ = budget.reserve(buffer_bytes, task);
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

  void append(const void* source, size_type bytes)
  {
    if(closed_) { throw eventError("write after close", path_); }
    const std::uint8_t* data = static_cast<const std::uint8_t*>(source);
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
    size_type records, size_type buffer_bytes, MemoryBudget& budget,
    const std::string& task) :
    path_(path), descriptor_(-1), width_(width), records_(records), seen_(0),
    buffer_(), buffered_(0), offset_(0), consumed_bytes_(0), cache_released_(0),
    reservation_()
  {
    if(width_ == 0) { throw eventError("zero-width event stream", path); }
    requireFileSize(path_, records_, width_);
    size_type capacity = std::max(width_, buffer_bytes - buffer_bytes % width_);
    reservation_ = budget.reserve(capacity, task);
    buffer_.resize(capacity);
    descriptor_ = ::open(path.c_str(), O_RDONLY);
    if(descriptor_ < 0) { throw eventError("cannot open event stream", path); }
    adviseSequential(descriptor_);
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
    readAll(descriptor_, buffer_.data(), bytes, path_);
    buffered_ = bytes; offset_ = 0; consumed_bytes_ += bytes;
    if(static_cast<off_t>(consumed_bytes_) - cache_released_ >= CACHE_FLUSH_BYTES)
    {
      off_t discard_end = std::max(cache_released_,
        static_cast<off_t>(consumed_bytes_) - CACHE_TAIL_BYTES);
      discardCache(descriptor_, cache_released_, discard_end - cache_released_);
      cache_released_ = discard_end;
    }
  }
};

ArtifactIdentity metadataArtifact()
{
  return ArtifactIdentity(FINAL_TASK, FINAL_PHASE, "metadata", "final-events-meta-v1");
}

ArtifactIdentity maskArtifact()
{
  return ArtifactIdentity(FINAL_TASK, FINAL_PHASE, "bwt-masks", "bwt-mask-u8-v1");
}

ArtifactIdentity edgeArtifact(size_type comp)
{
  return ArtifactIdentity(FINAL_TASK, FINAL_PHASE,
    "edge-destinations-" + std::to_string(comp), "path-rank-u64le-v1");
}

ArtifactIdentity streamArtifact(const std::string& name, const std::string& kind)
{
  return ArtifactIdentity(FINAL_TASK, FINAL_PHASE, name, kind);
}

physical_shard_id_t edgeShard(size_type comp)
{
  return physical_shard_id_t(100 + comp);
}

BuildWorkspace::ArtifactRef
checkpointFile(BuildWorkspace& workspace, const ArtifactIdentity& identity,
  physical_shard_id_t shard, const std::string& path, size_type records,
  size_type buffer_bytes, const std::string& sort_order)
{
  // Event streams are closed and immutable at this point. Give the raw inode
  // its semantic workspace name instead of writing a second full copy. The
  // workspace helper retains a bounded-copy fallback for another filesystem
  // and computes the committed checksum in either case.
  static_cast<void>(sort_order);
  return workspace.adopt_raw_payload(identity, logical_file_id_t(0), shard,
    path, records, rawFileBytes(path),
    std::max(static_cast<size_type>(1), buffer_bytes));
}

void
restoreFile(const BuildWorkspace& workspace, const ArtifactIdentity& identity,
  physical_shard_id_t shard, const std::string& path, size_type records,
  size_type width, size_type buffer_bytes, bool verify_checksum)
{
  const size_type expected_bytes = checkedBytes(records, width, path);
  const std::string source = workspace.artifact_path(identity,
    logical_file_id_t(0), shard);
  struct stat status;
  if(::stat(source.c_str(), &status) != 0 || status.st_size < 0)
  {
    throw eventError("cannot stat committed event stream", source);
  }

  if(static_cast<std::uintmax_t>(status.st_size) == expected_bytes)
  {
    // New workspaces store the raw immutable payload. A normal resume validates
    // task identity and exact length in O(metadata); --verify-workspace asks
    // restore_adopted_payload() for the full checksum scan.
    workspace.restore_adopted_payload(identity, logical_file_id_t(0), shard,
      path, records, expected_bytes,
      std::max(static_cast<size_type>(1), buffer_bytes), verify_checksum);
  }
  else
  {
    // Backward compatibility for v1 workspaces whose event payloads were
    // wrapped in a generic artifact header/footer. This also lets an active
    // construction upgrade without discarding a completed final-event scan.
    workspace.restore_artifact(identity, logical_file_id_t(0), shard, path,
      std::max(static_cast<size_type>(1), buffer_bytes));
  }
}

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

  const size_type block_size = 4096;
  const size_type blocks = (ones + block_size - 1) / block_size;
  const size_type logn = sdsl::bits::hi(((bits + 63) / 64) * 64) + 1;
  const size_type logn4 = logn * logn * logn * logn;
  std::vector<size_type> superblocks; superblocks.reserve(blocks);
  GCSA::bit_vector long_blocks(blocks, 0);
  std::FILE* payload = std::tmpfile();
  if(payload == nullptr) { throw eventError("cannot create select payload spool"); }

  try
  {
    std::array<size_type, block_size> positions;
    size_type next = 0, seen = 0;
    bool available = next_position(next);
    for(size_type block = 0; block < blocks; block++)
    {
      const size_type count = std::min(block_size, ones - seen);
      for(size_type i = 0; i < count; i++)
      {
        if(!available || (i > 0 && next <= positions[i - 1]))
        {
          throw eventError("select positions are not strictly increasing");
        }
        positions[i] = next; seen++; available = next_position(next);
      }
      superblocks.push_back(positions[0]);
      const size_type span = positions[count - 1] - positions[0];
      // init_fast() encodes the final partial block as long.
      const bool fast = (bits >= 100000);
      const bool is_long = (span > logn4 || (fast && count < block_size));
      long_blocks[block] = is_long;
      sdsl::int_vector<0> encoded(is_long ? block_size : 64, 0,
        is_long ? sdsl::bits::hi(fast ? bits - 1 : positions[count - 1]) + 1 :
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
    }
    if(available || seen != ones)
    {
      throw eventError("select position count does not match metadata");
    }
    sdsl::int_vector<0> superblock_vector(blocks, 0, logn);
    for(size_type i = 0; i < blocks; i++) { superblock_vector[i] = superblocks[i]; }
    superblock_vector.serialize(out);
    bool any_long = false;
    for(size_type i = 0; i < blocks; i++) { any_long = any_long || long_blocks[i]; }
    if(any_long) { long_blocks.serialize(out); }
    else { GCSA::bit_vector empty; empty.serialize(out); }
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
        out.write(buffer.data(), count);
        if(!out) { throw eventError("cannot serialize select payload"); }
      }
      if(count < buffer.size())
      {
        if(std::ferror(payload)) { throw eventError("cannot read select payload spool"); }
        break;
      }
    }
    std::fclose(payload);
  }
  catch(...)
  {
    std::fclose(payload); throw;
  }
}

void
validatePayloads(const FinalEventFiles& files, const FinalEventMetadata& metadata)
{
  if(files.edge_destinations.size() != metadata.sigma)
  {
    throw eventError("edge stream count does not match alphabet");
  }
  requireFileSize(files.metadata, 8 + EVENT_METADATA_WORDS * 8, 1);
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
    std::min(parameters.getIOBufferSize(), share));
}

class RedundancyPositionReader
{
public:
  RedundancyPositionReader(const std::string& filename, size_type paths,
    size_type redundant, const ConstructionParameters& parameters) :
    slots_(paths > 0 ? paths - 1 : 0), redundant_(redundant), slot_(0),
    cumulative_(0), budget_(parameters.getMemoryLimitBytes(),
      parameters.getMemoryLimitBytes() / 8), input_(new BufferedEventReader(filename, 8, redundant,
      componentBufferBytes(parameters, 1), budget_, "final-redundancy-reader")),
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

size_type
residentMemoryUsage()
{
#ifdef __linux__
  std::ifstream input("/proc/self/statm");
  size_type total = 0, resident = 0;
  if(input >> total >> resident)
  {
    long page_size = ::sysconf(_SC_PAGESIZE);
    if(page_size > 0 && resident <= std::numeric_limits<size_type>::max() /
      static_cast<size_type>(page_size))
    {
      return resident * static_cast<size_type>(page_size);
    }
  }
#endif
  return 0;
}

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
    this->reservation = memory.reserve(this->budget, "final-start-node-set");

    // Leave half the reservation unused during collection. It covers vector
    // allocator slack and the small file/stream objects used when the vector
    // is flushed. The external sorter receives the entire reservation only
    // after this vector has released its storage.
    this->collection_capacity = std::max(static_cast<size_type>(1),
      this->budget / (2 * sizeof(node_type)));
    this->values.reserve(this->collection_capacity);
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
    this->values.clear();
    if(this->values.capacity() < this->collection_capacity)
    {
      this->values.reserve(this->collection_capacity);
    }
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
        });
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
    componentBufferBytes(parameters, 1), budget, "final-bwt-mask-reader");
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
    componentBufferBytes(parameters, 1), budget, "final-sample-id-reader");
  serializePackedVector(out, samples, sample_bits, [&]()
  {
    std::uint64_t node = 0;
    ids.nextInteger(node); return node;
  });
  ids.finish();
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
  metadata(TempFile::getName("gcsa_final_event_meta")),
  bwt_masks(TempFile::getName("gcsa_final_bwt_masks")),
  edge_destinations(),
  sample_positions(TempFile::getName("gcsa_final_sample_paths")),
  sample_ids(TempFile::getName("gcsa_final_sample_ids")),
  sample_ends(TempFile::getName("gcsa_final_sample_ends")),
  occurrences(TempFile::getName("gcsa_final_occurrences")),
  redundant(TempFile::getName("gcsa_final_redundant")), delete_files(true)
{
  edge_destinations.reserve(sigma);
  for(size_type comp = 0; comp < sigma; comp++)
  {
    edge_destinations.push_back(TempFile::getName("gcsa_final_edges"));
  }
}

FinalEventFiles::FinalEventFiles(FinalEventFiles&& source) :
  metadata(std::move(source.metadata)), bwt_masks(std::move(source.bwt_masks)),
  edge_destinations(std::move(source.edge_destinations)),
  sample_positions(std::move(source.sample_positions)),
  sample_ids(std::move(source.sample_ids)), sample_ends(std::move(source.sample_ends)),
  occurrences(std::move(source.occurrences)), redundant(std::move(source.redundant)),
  delete_files(source.delete_files)
{
  source.delete_files = false;
}

FinalEventFiles&
FinalEventFiles::operator=(FinalEventFiles&& source)
{
  if(this != &source)
  {
    this->clear();
    this->metadata = std::move(source.metadata);
    this->bwt_masks = std::move(source.bwt_masks);
    this->edge_destinations = std::move(source.edge_destinations);
    this->sample_positions = std::move(source.sample_positions);
    this->sample_ids = std::move(source.sample_ids);
    this->sample_ends = std::move(source.sample_ends);
    this->occurrences = std::move(source.occurrences);
    this->redundant = std::move(source.redundant);
    this->delete_files = source.delete_files; source.delete_files = false;
  }
  return *this;
}

FinalEventFiles::~FinalEventFiles()
{
  this->clear();
}

void
FinalEventFiles::clear()
{
  if(this->delete_files)
  {
    TempFile::remove(this->metadata); TempFile::remove(this->bwt_masks);
    for(std::string& path : this->edge_destinations) { TempFile::remove(path); }
    TempFile::remove(this->sample_positions); TempFile::remove(this->sample_ids);
    TempFile::remove(this->sample_ends); TempFile::remove(this->occurrences);
    TempFile::remove(this->redundant);
  }
  this->metadata.clear(); this->bwt_masks.clear(); this->edge_destinations.clear();
  this->sample_positions.clear(); this->sample_ids.clear(); this->sample_ends.clear();
  this->occurrences.clear(); this->redundant.clear(); this->delete_files = false;
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
    masks.reset(new BufferedEventWriter(files.bwt_masks, buffer_bytes, budget,
      "final-event-bwt-mask"));
    edges.reserve(sigma);
    for(size_type comp = 0; comp < sigma; comp++)
    {
      edges.emplace_back(new BufferedEventWriter(files.edge_destinations[comp],
        buffer_bytes, budget, "final-event-edge"));
    }
    sample_positions.reset(new BufferedEventWriter(files.sample_positions,
      buffer_bytes, budget, "final-event-sampled-path"));
    sample_ids.reset(new BufferedEventWriter(files.sample_ids, buffer_bytes,
      budget, "final-event-sample-id"));
    sample_ends.reset(new BufferedEventWriter(files.sample_ends, buffer_bytes,
      budget, "final-event-sample-end"));
    occurrences.reset(new BufferedEventWriter(files.occurrences, buffer_bytes,
      budget, "final-event-occurrence"));
    redundant.reset(new BufferedEventWriter(files.redundant, buffer_bytes,
      budget, "final-event-redundancy"));
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
    size_type budget = std::min(usable_memory, parameters.getSortRunSize());
    const size_type sorter_descriptors = 4;
    if(parameters.getMaxOpenFiles() < sorter_descriptors + 2 * 2)
    {
      throw eventError("max-open-files cannot hold redundancy sort streams");
    }
    size_type fan_in = std::min(parameters.getMergeFanIn(),
      (parameters.getMaxOpenFiles() - sorter_descriptors) / 2);
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
writeFinalEventMetadata(const FinalEventFiles& files,
  const FinalEventMetadata& metadata)
{
  validateMetadata(metadata);
  std::vector<std::uint8_t> bytes(8 + EVENT_METADATA_WORDS * 8, 0);
  std::memcpy(bytes.data(), EVENT_MAGIC, 8);
  size_type word = 0;
  const auto append = [&bytes, &word](size_type value)
  {
    put64(bytes.data() + 8 + 8 * word, value); word++;
  };
  append(EVENT_VERSION); append(metadata.paths); append(metadata.sigma);
  append(metadata.fast_chars); append(metadata.total_edges);
  append(metadata.sampled_paths); append(metadata.sample_ids);
  append(metadata.sample_bits); append(metadata.occurrence_items);
  append(metadata.occurrence_extra); append(metadata.redundant);
  for(size_type comp = 0; comp < FinalEventMetadata::MAX_SIGMA; comp++)
  {
    append(metadata.bwt_counts[comp]);
  }
  if(word != EVENT_METADATA_WORDS) { throw eventError("internal metadata width mismatch"); }

  int descriptor = ::open(files.metadata.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
  if(descriptor < 0) { throw eventError("cannot create metadata", files.metadata); }
  try
  {
    writeAll(descriptor, bytes.data(), bytes.size(), files.metadata);
    if(::fdatasync(descriptor) != 0) { throw eventError("metadata fdatasync failed", files.metadata); }
    if(::close(descriptor) != 0) { descriptor = -1; throw eventError("metadata close failed", files.metadata); }
    descriptor = -1;
  }
  catch(...)
  {
    if(descriptor >= 0) { ::close(descriptor); }
    throw;
  }
}

FinalEventMetadata
readFinalEventMetadata(const FinalEventFiles& files)
{
  requireFileSize(files.metadata, 8 + EVENT_METADATA_WORDS * 8, 1);
  std::vector<std::uint8_t> bytes(8 + EVENT_METADATA_WORDS * 8, 0);
  int descriptor = ::open(files.metadata.c_str(), O_RDONLY);
  if(descriptor < 0) { throw eventError("cannot open metadata", files.metadata); }
  try
  {
    readAll(descriptor, bytes.data(), bytes.size(), files.metadata);
    if(::close(descriptor) != 0) { descriptor = -1; throw eventError("metadata close failed", files.metadata); }
    descriptor = -1;
  }
  catch(...)
  {
    if(descriptor >= 0) { ::close(descriptor); }
    throw;
  }
  if(std::memcmp(bytes.data(), EVENT_MAGIC, 8) != 0)
  {
    throw eventError("invalid metadata magic", files.metadata);
  }
  size_type word = 0;
  const auto next = [&bytes, &word]() -> size_type
  {
    return static_cast<size_type>(get64(bytes.data() + 8 + 8 * word++));
  };
  if(next() != EVENT_VERSION) { throw eventError("unsupported metadata version", files.metadata); }
  FinalEventMetadata metadata;
  metadata.paths = next(); metadata.sigma = next(); metadata.fast_chars = next();
  metadata.total_edges = next(); metadata.sampled_paths = next();
  metadata.sample_ids = next(); metadata.sample_bits = next();
  metadata.occurrence_items = next(); metadata.occurrence_extra = next();
  metadata.redundant = next();
  for(size_type comp = 0; comp < FinalEventMetadata::MAX_SIGMA; comp++)
  {
    metadata.bwt_counts[comp] = next();
  }
  validateMetadata(metadata); return metadata;
}

void
checkpointFinalEvents(BuildWorkspace& workspace,
  const FinalEventFiles& files, const FinalEventMetadata& metadata,
  size_type buffer_bytes)
{
  validatePayloads(files, metadata);
  std::vector<BuildWorkspace::ArtifactRef> artifacts;
  artifacts.push_back(checkpointFile(workspace, metadataArtifact(),
    physical_shard_id_t(0), files.metadata, 1, buffer_bytes, "metadata"));
  artifacts.push_back(checkpointFile(workspace, maskArtifact(),
    physical_shard_id_t(1), files.bwt_masks, metadata.paths, buffer_bytes,
    "path-rank"));
  for(size_type comp = 0; comp < metadata.sigma; comp++)
  {
    artifacts.push_back(checkpointFile(workspace, edgeArtifact(comp), edgeShard(comp),
      files.edge_destinations[comp], metadata.bwt_counts[comp], buffer_bytes,
      "source-path-rank"));
  }
  artifacts.push_back(checkpointFile(workspace,
    streamArtifact("sampled-paths", "path-rank-u64le-v1"), physical_shard_id_t(2),
    files.sample_positions, metadata.sampled_paths, buffer_bytes, "path-rank"));
  artifacts.push_back(checkpointFile(workspace,
    streamArtifact("sample-ids", "node-id-u64le-v1"), physical_shard_id_t(3),
    files.sample_ids, metadata.sample_ids, buffer_bytes, "sample-order"));
  artifacts.push_back(checkpointFile(workspace,
    streamArtifact("sample-ends", "sample-rank-u64le-v1"), physical_shard_id_t(4),
    files.sample_ends, metadata.sampled_paths, buffer_bytes, "sample-order"));
  artifacts.push_back(checkpointFile(workspace,
    streamArtifact("occurrences", "path-value-u64le-v1"), physical_shard_id_t(5),
    files.occurrences, metadata.occurrence_items, buffer_bytes, "path-rank"));
  artifacts.push_back(checkpointFile(workspace,
    streamArtifact("redundancy", "path-rank-u64le-v1"), physical_shard_id_t(6),
    files.redundant, metadata.redundant, buffer_bytes, "path-rank"));
  workspace.commit_task(FINAL_TASK, FINAL_PHASE, artifacts);
}

bool
restoreFinalEvents(const BuildWorkspace& workspace,
  FinalEventFiles& files, FinalEventMetadata& metadata,
  size_type expected_paths, size_type expected_sigma, size_type buffer_bytes,
  bool verify_checksum)
{
  if(!workspace.task_completed(FINAL_TASK, FINAL_PHASE)) { return false; }
  restoreFile(workspace, metadataArtifact(), physical_shard_id_t(0),
    files.metadata, 1, 8 + EVENT_METADATA_WORDS * 8, buffer_bytes,
    verify_checksum);
  metadata = readFinalEventMetadata(files);
  if(metadata.paths != expected_paths || metadata.sigma != expected_sigma ||
     files.edge_destinations.size() != metadata.sigma)
  {
    throw eventError("committed events do not match the merged graph");
  }
  restoreFile(workspace, maskArtifact(), physical_shard_id_t(1),
    files.bwt_masks, metadata.paths, 1, buffer_bytes, verify_checksum);
  for(size_type comp = 0; comp < metadata.sigma; comp++)
  {
    restoreFile(workspace, edgeArtifact(comp), edgeShard(comp),
      files.edge_destinations[comp], metadata.bwt_counts[comp], 8,
      buffer_bytes, verify_checksum);
  }
  restoreFile(workspace, streamArtifact("sampled-paths", "path-rank-u64le-v1"),
    physical_shard_id_t(2), files.sample_positions, metadata.sampled_paths, 8,
    buffer_bytes, verify_checksum);
  restoreFile(workspace, streamArtifact("sample-ids", "node-id-u64le-v1"),
    physical_shard_id_t(3), files.sample_ids, metadata.sample_ids, 8,
    buffer_bytes, verify_checksum);
  restoreFile(workspace, streamArtifact("sample-ends", "sample-rank-u64le-v1"),
    physical_shard_id_t(4), files.sample_ends, metadata.sampled_paths, 8,
    buffer_bytes, verify_checksum);
  restoreFile(workspace, streamArtifact("occurrences", "path-value-u64le-v1"),
    physical_shard_id_t(5), files.occurrences, metadata.occurrence_items, 16,
    buffer_bytes, verify_checksum);
  restoreFile(workspace, streamArtifact("redundancy", "path-rank-u64le-v1"),
    physical_shard_id_t(6), files.redundant, metadata.redundant, 8,
    buffer_bytes, verify_checksum);
  validatePayloads(files, metadata); return true;
}

void
buildFinalComponents(GCSA& index, const Alphabet& source_alphabet,
  const FinalEventFiles& files, const FinalEventMetadata& metadata,
  const ConstructionParameters& parameters)
{
  validateMetadata(metadata);
  validatePayloads(files, metadata);
  if(source_alphabet.sigma != metadata.sigma ||
     source_alphabet.fast_chars != metadata.fast_chars)
  {
    throw eventError("event alphabet does not match the input alphabet");
  }

  sdsl::int_vector<64> counts(metadata.sigma, 0);
  for(size_type comp = 0; comp < metadata.sigma; comp++)
  {
    counts[comp] = metadata.bwt_counts[comp];
  }
  index.alpha = Alphabet(counts, source_alphabet.char2comp,
    source_alphabet.comp2char);

  index.fast_bwt.resize(metadata.sigma); index.fast_rank.resize(metadata.sigma);
  index.sparse_bwt.resize(metadata.sigma); index.sparse_rank.resize(metadata.sigma);
  for(size_type comp = 0; comp < metadata.sigma; comp++)
  {
    MemoryBudget budget(parameters.getMemoryLimitBytes(),
      parameters.getMemoryLimitBytes() / 8);
    size_type buffer_bytes = componentBufferBytes(parameters, 1);
    BufferedEventReader masks(files.bwt_masks, 1, metadata.paths,
      buffer_bytes, budget, "final-bwt-mask-reader");
    if(comp > 0 && comp <= metadata.fast_chars)
    {
      GCSA::bit_vector dense(metadata.paths, 0);
      for(size_type path = 0; path < metadata.paths; path++)
      {
        std::uint8_t mask = 0; masks.nextByte(mask);
        if(mask & (static_cast<size_type>(1) << comp)) { dense[path] = 1; }
      }
      masks.finish(); index.fast_bwt[comp] = dense; sdsl::util::clear(dense);
    }
    else
    {
      sdsl::sd_vector_builder builder(metadata.paths, metadata.bwt_counts[comp]);
      size_type observed = 0;
      for(size_type path = 0; path < metadata.paths; path++)
      {
        std::uint8_t mask = 0; masks.nextByte(mask);
        if(mask & (static_cast<size_type>(1) << comp))
        {
          // Metadata is durable input here. Use the checked insertion so a
          // corrupt one-count cannot overrun sd_vector_builder before the
          // aggregate count check below reports the mismatch.
          builder.set(path); observed++;
        }
      }
      masks.finish();
      if(observed != metadata.bwt_counts[comp])
      {
        throw eventError("BWT mask count does not match metadata");
      }
      index.sparse_bwt[comp] = GCSA::sparse_vector(builder);
    }
  }

  // Each component's source ranks are nondecreasing. A bounded k-way merge
  // reconstructs outdegrees without a full CounterArray.
  {
    MemoryBudget budget(parameters.getMemoryLimitBytes(),
      parameters.getMemoryLimitBytes() / 8);
    size_type buffer_bytes = componentBufferBytes(parameters, metadata.sigma);
    std::vector<std::unique_ptr<BufferedEventReader>> readers;
    std::vector<std::uint64_t> current(metadata.sigma, 0);
    std::vector<bool> available(metadata.sigma, false);
    readers.reserve(metadata.sigma);
    for(size_type comp = 0; comp < metadata.sigma; comp++)
    {
      readers.emplace_back(new BufferedEventReader(files.edge_destinations[comp],
        8, metadata.bwt_counts[comp], buffer_bytes, budget, "final-edge-reader"));
      available[comp] = readers.back()->nextInteger(current[comp]);
    }
    GCSA::bit_vector edge_buffer(metadata.total_edges, 0);
    size_type tail = 0;
    for(size_type path = 0; path < metadata.paths; path++)
    {
      size_type degree = 0;
      for(size_type comp = 0; comp < metadata.sigma; comp++)
      {
        if(available[comp] && current[comp] < path)
        {
          throw eventError("edge stream is not nondecreasing");
        }
        while(available[comp] && current[comp] == path)
        {
          degree++; available[comp] = readers[comp]->nextInteger(current[comp]);
        }
      }
      if(degree == 0) { throw eventError("path has no outgoing edge"); }
      if(degree > metadata.total_edges || tail > metadata.total_edges - degree)
      {
        throw eventError("edge stream exceeds declared total");
      }
      tail += degree; edge_buffer[tail - 1] = 1;
    }
    for(size_type comp = 0; comp < metadata.sigma; comp++)
    {
      if(available[comp])
      {
        throw eventError("edge rank exceeds path universe");
      }
      readers[comp]->finish();
    }
    if(tail != metadata.total_edges) { throw eventError("edge total mismatch"); }
    index.edges = edge_buffer; sdsl::util::clear(edge_buffer);
  }

  // Sampled path positions, node identifiers, and path boundaries are separate
  // monotone streams. This avoids the historical unbounded sample_buffer.
  {
    MemoryBudget budget(parameters.getMemoryLimitBytes(),
      parameters.getMemoryLimitBytes() / 8);
    size_type buffer_bytes = componentBufferBytes(parameters, 1);
    BufferedEventReader positions(files.sample_positions, 8,
      metadata.sampled_paths, buffer_bytes, budget, "final-sample-path-reader");
    GCSA::bit_vector sampled(metadata.paths, 0);
    std::uint64_t previous = 0;
    for(size_type i = 0; i < metadata.sampled_paths; i++)
    {
      std::uint64_t path = 0; positions.nextInteger(path);
      if(path >= metadata.paths || (i > 0 && path <= previous))
      {
        throw eventError("sampled path stream is not strictly increasing");
      }
      sampled[path] = 1; previous = path;
    }
    positions.finish(); index.sampled_paths = sampled; sdsl::util::clear(sampled);

    index.stored_samples = sdsl::int_vector<0>(metadata.sample_ids, 0,
      metadata.sample_bits);
    BufferedEventReader ids(files.sample_ids, 8, metadata.sample_ids,
      buffer_bytes, budget, "final-sample-id-reader");
    for(size_type i = 0; i < metadata.sample_ids; i++)
    {
      std::uint64_t node = 0; ids.nextInteger(node);
      if(metadata.sample_bits < 64 &&
         node >= (static_cast<std::uint64_t>(1) << metadata.sample_bits))
      {
        throw eventError("sample identifier exceeds its declared width");
      }
      index.stored_samples[i] = node;
    }
    ids.finish();

    index.samples = GCSA::bit_vector(metadata.sample_ids, 0);
    BufferedEventReader ends(files.sample_ends, 8, metadata.sampled_paths,
      buffer_bytes, budget, "final-sample-end-reader");
    previous = 0;
    for(size_type i = 0; i < metadata.sampled_paths; i++)
    {
      std::uint64_t end = 0; ends.nextInteger(end);
      if(end >= metadata.sample_ids || (i > 0 && end <= previous))
      {
        throw eventError("sample endpoint stream is not strictly increasing");
      }
      index.samples[end] = 1; previous = end;
    }
    ends.finish();
    if(metadata.sampled_paths > 0 && previous + 1 != metadata.sample_ids)
    {
      throw eventError("final sample endpoint does not cover all sample IDs");
    }
  }

  // Construct Sada-S directly from nonzero occurrence events.
  {
    MemoryBudget budget(parameters.getMemoryLimitBytes(),
      parameters.getMemoryLimitBytes() / 8);
    BufferedEventReader input(files.occurrences, 16, metadata.occurrence_items,
      componentBufferBytes(parameters, 1), budget, "final-occurrence-reader");
    sdsl::sd_vector_builder filter(metadata.paths, metadata.occurrence_items);
    sdsl::sd_vector_builder values(metadata.occurrence_extra,
      metadata.occurrence_items);
    size_type tail = 0; std::uint64_t previous = 0;
    for(size_type i = 0; i < metadata.occurrence_items; i++)
    {
      std::uint64_t path = 0, extra = 0; input.nextPair(path, extra);
      if(path >= metadata.paths || extra == 0 || (i > 0 && path <= previous))
      {
        throw eventError("invalid occurrence event order or value");
      }
      if(tail > metadata.occurrence_extra ||
         extra > metadata.occurrence_extra - tail)
      {
        throw eventError("occurrence value total overflows metadata");
      }
      // Checked insertions turn malformed order, universe, or metadata counts
      // into an exception instead of undefined writes inside SDSL.
      filter.set(path);
      tail += extra; values.set(tail - 1); previous = path;
    }
    input.finish();
    if(tail != metadata.occurrence_extra)
    {
      throw eventError("occurrence value total mismatch");
    }
    index.extra_pointers.filter = SadaSparse::sd_vector(filter);
    index.extra_pointers.values = SadaSparse::sd_vector(values);
    sdsl::util::init_support(index.extra_pointers.filter_rank,
      &(index.extra_pointers.filter));
    sdsl::util::init_support(index.extra_pointers.value_select,
      &(index.extra_pointers.values));
  }

  // Redundancy events were externally sorted after the scan. Duplicate ranks
  // become the unary zero run for that suffix-tree slot.
  {
    if(metadata.paths == 0 && metadata.redundant != 0)
    {
      throw eventError("redundancy events exist for an empty graph");
    }
    if(metadata.redundant > std::numeric_limits<size_type>::max() -
       (metadata.paths > 0 ? metadata.paths - 1 : 0))
    {
      throw eventError("redundancy bitvector length overflows");
    }
    size_type slots = (metadata.paths > 0 ? metadata.paths - 1 : 0);
    index.redundant_pointers.data = GCSA::bit_vector(slots + metadata.redundant, 0);
    MemoryBudget budget(parameters.getMemoryLimitBytes(),
      parameters.getMemoryLimitBytes() / 8);
    BufferedEventReader input(files.redundant, 8, metadata.redundant,
      componentBufferBytes(parameters, 1), budget, "final-redundancy-reader");
    std::uint64_t current = 0;
    bool available = input.nextInteger(current);
    size_type cumulative = 0;
    for(size_type slot = 0; slot < slots; slot++)
    {
      if(available && current < slot)
      {
        throw eventError("redundancy stream is not nondecreasing");
      }
      while(available && current == slot)
      {
        cumulative++; available = input.nextInteger(current);
      }
      index.redundant_pointers.data[slot + cumulative] = 1;
    }
    if(available || cumulative != metadata.redundant)
    {
      throw eventError("redundancy rank exceeds suffix-tree slot universe");
    }
    input.finish();
    sdsl::util::init_support(index.redundant_pointers.select,
      &(index.redundant_pointers.data));
  }

  for(size_type comp = 0; comp < metadata.sigma; comp++)
  {
    sdsl::util::init_support(index.fast_rank[comp], &(index.fast_bwt[comp]));
    sdsl::util::init_support(index.sparse_rank[comp], &(index.sparse_bwt[comp]));
  }
  sdsl::util::init_support(index.edge_rank, &(index.edges));
  sdsl::util::init_support(index.sampled_path_rank, &(index.sampled_paths));
  sdsl::util::init_support(index.sample_select, &(index.samples));
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
  size_type peak_before = memoryUsage();
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

    // The on-disk order is all BWTs followed by all ranks. Reconstructing a
    // component for its rank costs one bounded stream pass but avoids retaining
    // all BWT members merely to reach that later section.
    for(size_type comp = 0; comp < metadata.sigma; comp++)
    {
      if(comp > 0 && comp <= metadata.fast_chars)
      {
        serializeFastBWTComponent(out, files.bwt_masks, metadata.paths,
          metadata.bwt_counts[comp], static_cast<comp_type>(comp), parameters);
        if(!out) { throw eventError("cannot write partial final index", partial); }
      }
      else
      {
        GCSA::fast_vector empty;
        write(empty);
      }
    }
    // rank_support_il serializes no payload; load() simply binds it to the BWT
    // already read above. Rebuilding each fast BWT solely to initialize an
    // empty serialization performed a full mask scan and allocation for no
    // output bytes.
    for(size_type comp = 0; comp < metadata.sigma; comp++)
    {
      GCSA::fast_vector::rank_1_type empty;
      write(empty);
    }
    for(size_type comp = 0; comp < metadata.sigma; comp++)
    {
      if(comp > 0 && comp <= metadata.fast_chars)
      {
        GCSA::sparse_vector empty;
        write(empty);
        continue;
      }
      MemoryBudget budget(parameters.getMemoryLimitBytes(),
        parameters.getMemoryLimitBytes() / 8);
      BufferedEventReader masks(files.bwt_masks, 1, metadata.paths,
        componentBufferBytes(parameters, 1), budget, "final-bwt-mask-reader");
      sdsl::sd_vector_builder builder(metadata.paths, metadata.bwt_counts[comp]);
      size_type observed = 0;
      for(size_type path = 0; path < metadata.paths; path++)
      {
        std::uint8_t mask = 0;
        masks.nextByte(mask);
        if(mask & (static_cast<size_type>(1) << comp))
        {
          builder.set(path);
          observed++;
        }
      }
      masks.finish();
      if(observed != metadata.bwt_counts[comp])
      {
        throw eventError("BWT mask count does not match metadata");
      }
      GCSA::sparse_vector bwt(builder);
      write(bwt);
    }
    for(size_type comp = 0; comp < metadata.sigma; comp++)
    {
      if(comp > 0 && comp <= metadata.fast_chars)
      {
        GCSA::sparse_vector::rank_1_type empty;
        write(empty);
        continue;
      }
      MemoryBudget budget(parameters.getMemoryLimitBytes(),
        parameters.getMemoryLimitBytes() / 8);
      BufferedEventReader masks(files.bwt_masks, 1, metadata.paths,
        componentBufferBytes(parameters, 1), budget, "final-bwt-mask-reader");
      sdsl::sd_vector_builder builder(metadata.paths, metadata.bwt_counts[comp]);
      size_type observed = 0;
      for(size_type path = 0; path < metadata.paths; path++)
      {
        std::uint8_t mask = 0;
        masks.nextByte(mask);
        if(mask & (static_cast<size_type>(1) << comp))
        {
          builder.set(path);
          observed++;
        }
      }
      masks.finish();
      if(observed != metadata.bwt_counts[comp])
      {
        throw eventError("BWT mask count does not match metadata");
      }
      GCSA::sparse_vector bwt(builder);
      GCSA::sparse_vector::rank_1_type rank;
      sdsl::util::init_support(rank, &bwt);
      write(rank);
    }

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
          8, metadata.bwt_counts[comp], buffer_bytes, budget, "final-edge-reader"));
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

    {
      MemoryBudget budget(parameters.getMemoryLimitBytes(),
        parameters.getMemoryLimitBytes() / 8);
      size_type buffer_bytes = componentBufferBytes(parameters, 1);
      {
        BufferedEventReader positions(files.sample_positions, 8,
          metadata.sampled_paths, buffer_bytes, budget,
          "final-sample-path-reader");
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

      GCSA::bit_vector ends_vector(metadata.sample_ids, 0);
      BufferedEventReader ends(files.sample_ends, 8, metadata.sampled_paths,
        buffer_bytes, budget, "final-sample-end-reader");
      std::uint64_t previous = 0;
      for(size_type i = 0; i < metadata.sampled_paths; i++)
      {
        std::uint64_t end = 0;
        ends.nextInteger(end);
        if(end >= metadata.sample_ids || (i > 0 && end <= previous))
        {
          throw eventError("sample endpoint stream is not strictly increasing");
        }
        ends_vector[end] = 1;
        previous = end;
      }
      ends.finish();
      if(metadata.sampled_paths > 0 && previous + 1 != metadata.sample_ids)
      {
        throw eventError("final sample endpoint does not cover all sample IDs");
      }
      write(ends_vector);
      GCSA::bit_vector::select_1_type select;
      sdsl::util::init_support(select, &ends_vector);
      write(select);
    }

    {
      MemoryBudget budget(parameters.getMemoryLimitBytes(),
        parameters.getMemoryLimitBytes() / 8);
      BufferedEventReader input(files.occurrences, 16, metadata.occurrence_items,
        componentBufferBytes(parameters, 1), budget, "final-occurrence-reader");
      sdsl::sd_vector_builder filter(metadata.paths, metadata.occurrence_items);
      sdsl::sd_vector_builder values(metadata.occurrence_extra,
        metadata.occurrence_items);
      size_type tail = 0;
      std::uint64_t previous = 0;
      for(size_type i = 0; i < metadata.occurrence_items; i++)
      {
        std::uint64_t path = 0, extra = 0;
        input.nextPair(path, extra);
        if(path >= metadata.paths || extra == 0 ||
           (i > 0 && path <= previous) || tail > metadata.occurrence_extra ||
           extra > metadata.occurrence_extra - tail)
        {
          throw eventError("invalid occurrence event");
        }
        filter.set(path);
        tail += extra;
        values.set(tail - 1);
        previous = path;
      }
      input.finish();
      if(tail != metadata.occurrence_extra)
      {
        throw eventError("occurrence value total mismatch");
      }
      SadaSparse pointers;
      pointers.filter = SadaSparse::sd_vector(filter);
      pointers.values = SadaSparse::sd_vector(values);
      sdsl::util::init_support(pointers.filter_rank, &(pointers.filter));
      sdsl::util::init_support(pointers.value_select, &(pointers.values));
      write(pointers);
    }
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
    if(Verbosity::level >= Verbosity::EXTENDED)
    {
      size_type peak_after = memoryUsage();
      std::cerr << "GCSA::storeFinalComponents(): RSS current "
                << inGigabytes(residentMemoryUsage()) << " GB, peak increment "
                << inGigabytes(peak_after >= peak_before ? peak_after - peak_before : 0)
                << " GB, goal " << inGigabytes(parameters.getMemoryLimitBytes())
                << " GB (SDSL allocations are outside the byte budget)" << std::endl;
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
