/*
  Copyright (c) 2026 GCSA2 contributors
*/

#include <gcsa/disk_array.h>

#include <gcsa/internal.h>

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <limits>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>

namespace gcsa
{

namespace
{

std::runtime_error
diskArrayError(const std::string& message, const std::string& filename)
{
  return std::runtime_error("DiskBackedArray64: " + message + ": " + filename);
}

void
preadAll(int descriptor, void* target, size_type bytes, off_t offset,
  const std::string& filename)
{
  std::uint8_t* data = static_cast<std::uint8_t*>(target);
  size_type done = 0;
  while(done < bytes)
  {
    ssize_t got = ::pread(descriptor, data + done, bytes - done, offset + done);
    if(got < 0 && errno == EINTR) { continue; }
    if(got <= 0) { throw diskArrayError("short block read", filename); }
    done += static_cast<size_type>(got);
  }
  DiskIO::read_volume += bytes;
}

void
pwriteAll(int descriptor, const void* source, size_type bytes, off_t offset,
  const std::string& filename)
{
  const std::uint8_t* data = static_cast<const std::uint8_t*>(source);
  size_type done = 0;
  while(done < bytes)
  {
    ssize_t written = ::pwrite(descriptor, data + done, bytes - done, offset + done);
    if(written < 0 && errno == EINTR) { continue; }
    if(written <= 0) { throw diskArrayError("block write failed", filename); }
    done += static_cast<size_type>(written);
  }
  DiskIO::write_volume += bytes;
}

} // namespace

constexpr size_type DiskBackedArray64::DEFAULT_BLOCK_BYTES;
constexpr std::uint64_t DiskBackedArray64::EMPTY_TAG;

size_type
DiskBackedArray64::minimumCacheBytes()
{
  return sizeof(std::uint64_t) + sizeof(std::uint64_t) + sizeof(std::uint8_t);
}

DiskBackedArray64::DiskBackedArray64(const std::string& filename,
  size_type elements, size_type cache_bytes, MemoryBudget& budget,
  bool create_file, size_type requested_block_bytes) :
  filename_(filename), descriptor_(-1), elements_(elements), array_bytes_(0),
  block_bytes_(0), block_elements_(0), slots_(0)
{
  if(elements > std::numeric_limits<size_type>::max() / sizeof(std::uint64_t))
  {
    throw diskArrayError("logical length overflows bytes", filename);
  }
  this->array_bytes_ = elements * sizeof(std::uint64_t);
  if(this->array_bytes_ > static_cast<size_type>(std::numeric_limits<off_t>::max()))
  {
    throw diskArrayError("logical length exceeds filesystem offsets", filename);
  }

  int flags = O_RDWR;
  if(create_file) { flags |= O_CREAT | O_TRUNC; }
  this->descriptor_ = ::open(filename.c_str(), flags, 0644);
  if(this->descriptor_ < 0) { throw diskArrayError("cannot open", filename); }

  try
  {
    if(create_file)
    {
      if(::ftruncate(this->descriptor_, static_cast<off_t>(this->array_bytes_)) != 0)
      {
        throw diskArrayError("cannot size", filename);
      }
    }
    else
    {
      struct stat status;
      if(::fstat(this->descriptor_, &status) != 0 || status.st_size < 0 ||
         static_cast<size_type>(status.st_size) != this->array_bytes_)
      {
        throw diskArrayError("file length does not match logical array", filename);
      }
    }

    if(this->elements_ == 0) { return; }
    if(cache_bytes < minimumCacheBytes())
    {
      throw diskArrayError("cache budget is too small", filename);
    }
    requested_block_bytes = std::max(sizeof(std::uint64_t), requested_block_bytes);
    requested_block_bytes -= requested_block_bytes % sizeof(std::uint64_t);
    size_type half_cache = std::max(sizeof(std::uint64_t), cache_bytes / 2);
    half_cache -= half_cache % sizeof(std::uint64_t);
    this->block_bytes_ = std::min(this->array_bytes_,
      std::min(requested_block_bytes, half_cache));
    this->block_bytes_ = std::max(sizeof(std::uint64_t), this->block_bytes_);
    this->block_elements_ = this->block_bytes_ / sizeof(std::uint64_t);

    size_type slot_bytes = this->block_bytes_ + sizeof(std::uint64_t) +
      sizeof(std::uint8_t);
    size_type logical_blocks = (this->array_bytes_ + this->block_bytes_ - 1) /
      this->block_bytes_;
    this->slots_ = std::min(logical_blocks, cache_bytes / slot_bytes);
    if(this->slots_ == 0)
    {
      // Shrink to the largest aligned block that leaves room for its tag.
      this->block_bytes_ = (cache_bytes - sizeof(std::uint64_t) -
        sizeof(std::uint8_t)) / sizeof(std::uint64_t) * sizeof(std::uint64_t);
      if(this->block_bytes_ == 0)
      {
        throw diskArrayError("cache budget cannot hold one value and tag", filename);
      }
      this->block_elements_ = this->block_bytes_ / sizeof(std::uint64_t);
      this->slots_ = 1;
    }

    size_type reserved = this->slots_ *
      (this->block_bytes_ + sizeof(std::uint64_t) + sizeof(std::uint8_t));
    this->reservation_ = budget.reserve(reserved, "disk-backed-array-cache");
    this->payload_.resize(this->slots_ * this->block_elements_, 0);
    this->tags_.resize(this->slots_, EMPTY_TAG);
    this->dirty_.resize(this->slots_, 0);
  }
  catch(...)
  {
    ::close(this->descriptor_); this->descriptor_ = -1;
    throw;
  }
}

DiskBackedArray64::~DiskBackedArray64()
{
  if(this->descriptor_ >= 0)
  {
    try { this->flush(false); } catch(...) { }
    ::close(this->descriptor_); this->descriptor_ = -1;
  }
}

void
DiskBackedArray64::writeSlot(size_type slot)
{
  if(slot >= this->slots_ || this->tags_[slot] == EMPTY_TAG || !this->dirty_[slot])
  {
    return;
  }
  size_type offset = this->tags_[slot] * this->block_bytes_;
  size_type bytes = std::min(this->block_bytes_, this->array_bytes_ - offset);
  pwriteAll(this->descriptor_, this->payload_.data() + slot * this->block_elements_,
    bytes, static_cast<off_t>(offset), this->filename_);
  this->dirty_[slot] = 0; this->stats_.block_writes++;
}

std::uint64_t*
DiskBackedArray64::load(size_type block)
{
  size_type slot = this->slotFor(block);
  if(this->tags_[slot] == block)
  {
    return this->payload_.data() + slot * this->block_elements_;
  }
  if(this->tags_[slot] != EMPTY_TAG)
  {
    this->writeSlot(slot); this->stats_.replacements++;
  }

  size_type offset = block * this->block_bytes_;
  size_type bytes = std::min(this->block_bytes_, this->array_bytes_ - offset);
  std::uint64_t* target = this->payload_.data() + slot * this->block_elements_;
  std::fill(target, target + this->block_elements_, 0);
  preadAll(this->descriptor_, target, bytes, static_cast<off_t>(offset),
    this->filename_);
  this->tags_[slot] = block; this->dirty_[slot] = 0;
  this->stats_.block_reads++;
  return target;
}

std::uint64_t
DiskBackedArray64::get(size_type index)
{
  if(index >= this->elements_) { throw std::out_of_range("DiskBackedArray64::get"); }
  size_type byte_offset = index * sizeof(std::uint64_t);
  size_type block = byte_offset / this->block_bytes_;
  size_type within = (byte_offset % this->block_bytes_) / sizeof(std::uint64_t);
  return this->load(block)[within];
}

void
DiskBackedArray64::set(size_type index, std::uint64_t value)
{
  if(index >= this->elements_) { throw std::out_of_range("DiskBackedArray64::set"); }
  size_type byte_offset = index * sizeof(std::uint64_t);
  size_type block = byte_offset / this->block_bytes_;
  size_type within = (byte_offset % this->block_bytes_) / sizeof(std::uint64_t);
  size_type slot = this->slotFor(block);
  this->load(block)[within] = value; this->dirty_[slot] = 1;
}

void
DiskBackedArray64::flush(bool sync_file)
{
  for(size_type slot = 0; slot < this->slots_; slot++) { this->writeSlot(slot); }
  if(sync_file && ::fdatasync(this->descriptor_) != 0)
  {
    throw diskArrayError("fdatasync failed", this->filename_);
  }
}

} // namespace gcsa
