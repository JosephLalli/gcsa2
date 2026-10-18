#ifndef GCSA_ADAPTIVE_RECORD_STORE_H
#define GCSA_ADAPTIVE_RECORD_STORE_H

#include <gcsa/resources.h>
#include <gcsa/utils.h>

#include <algorithm>
#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

namespace gcsa
{

namespace adaptive_record_store_detail
{

int openRawRecordFile(const std::string& filename);
void closeRawRecordFile(int descriptor) noexcept;
void readRawRecords(int descriptor, void* target, size_type bytes,
  size_type offset, const std::string& filename);
void writeRawRecords(int descriptor, const void* source, size_type bytes,
  size_type offset, const std::string& filename);
void copyRawRecords(int source, int destination, size_type bytes,
  const std::string& source_filename, const std::string& destination_filename);

} // namespace adaptive_record_store_detail

/*
  A raw-record sequence that stays resident while its exact vector capacity is
  admitted by a shared MemoryBudget. If a normal append cannot admit its next
  capacity, it writes the existing native records to a raw temporary file and
  continues there. There is no quota per store: every growth asks the shared
  budget, while the old allocation remains charged until the replacement is
  constructed.

  The record representation is intentionally the raw sizeof(Element) layout
  used by the existing temporary path/rank streams. It therefore does not add
  a serialization format or a trivially-copyable restriction that the legacy
  streams do not have. Structural mutation is single-threaded. After callers
  pre-admit with tryReserveRecords() and resize(), they may fill disjoint
  residentData() elements in parallel, but must not append, resize, or
  materialize until those writes and any readers have finished.
*/
template<class Element>
class AdaptiveRecordStore
{
public:
  enum class FileOwnership { NONE, OWNED_TEMPORARY, BORROWED };

  AdaptiveRecordStore(MemoryBudget& memory, const std::string& task) :
    memory_(&memory), task_(task), data_(nullptr), count_(0), capacity_(0),
    descriptor_(-1), filename_(), ownership_(FileOwnership::NONE),
    resident_memory_()
  {
  }

  ~AdaptiveRecordStore()
  {
    this->clear();
  }

  AdaptiveRecordStore(const AdaptiveRecordStore&) = delete;
  AdaptiveRecordStore& operator=(const AdaptiveRecordStore&) = delete;

  bool resident() const { return (this->descriptor_ < 0); }
  bool spilled() const { return !(this->resident()); }
  size_type count() const { return this->count_; }
  size_type capacity() const { return this->capacity_; }
  FileOwnership fileOwnership() const { return this->ownership_; }

  // This only admits resident capacity. A failed admission leaves the store
  // resident and does not create a temporary file.
  bool tryReserveRecords(size_type records)
  {
    if(this->spilled()) { return false; }
    if(records <= this->capacity_) { return true; }
    return this->growResident(records);
  }

  // The caller must have admitted capacity with tryReserveRecords().
  void resize(size_type records)
  {
    if(this->spilled())
    {
      throw std::logic_error("AdaptiveRecordStore::resize(): store has spilled");
    }
    if(records > this->capacity_)
    {
      throw std::logic_error("AdaptiveRecordStore::resize(): capacity was not admitted");
    }
    if(records < this->count_)
    {
      for(size_type i = this->count_; i > records; i--)
      {
        std::allocator_traits<allocator_type>::destroy(this->allocator_, this->data_ + i - 1);
      }
      this->count_ = records;
      return;
    }

    const size_type old_count = this->count_;
    try
    {
      for(; this->count_ < records; this->count_++)
      {
        std::allocator_traits<allocator_type>::construct(this->allocator_,
          this->data_ + this->count_);
      }
    }
    catch(...)
    {
      while(this->count_ > old_count)
      {
        this->count_--;
        std::allocator_traits<allocator_type>::destroy(this->allocator_,
          this->data_ + this->count_);
      }
      throw;
    }
  }

  Element* residentData()
  {
    this->requireResident("residentData");
    return this->data_;
  }

  const Element* residentData() const
  {
    this->requireResident("residentData");
    return this->data_;
  }

  void append(const Element& value)
  {
    this->append(&value, 1);
  }

  void append(Element&& value)
  {
    this->appendValue(std::move(value));
  }

  // A spilled range is issued as one contiguous raw write. This preserves the
  // append-only byte-stream behavior used by the old temporary files.
  void append(const Element* values, size_type records)
  {
    if(records == 0) { return; }
    if(values == nullptr)
    {
      throw std::invalid_argument("AdaptiveRecordStore::append(): null source");
    }
    if(records > std::numeric_limits<size_type>::max() - this->count_)
    {
      throw std::length_error("AdaptiveRecordStore::append(): record count overflows");
    }

    const size_type required = this->count_ + records;
    this->ensureAppendCapacity(required);
    if(this->spilled())
    {
      adaptive_record_store_detail::writeRawRecords(this->descriptor_, values,
        this->bytesFor(records), this->bytesFor(this->count_), this->filename_);
      this->count_ = required;
      return;
    }

    const size_type old_count = this->count_;
    try
    {
      for(; this->count_ < required; this->count_++)
      {
        std::allocator_traits<allocator_type>::construct(this->allocator_,
          this->data_ + this->count_, values[this->count_ - old_count]);
      }
    }
    catch(...)
    {
      while(this->count_ > old_count)
      {
        this->count_--;
        std::allocator_traits<allocator_type>::destroy(this->allocator_,
          this->data_ + this->count_);
      }
      throw;
    }
  }

  // Reads are expressed in records, rather than bytes, at both resident and
  // spilled storage. The raw copy matches existing temporary-file behavior.
  void readAt(size_type offset, Element* target, size_type records) const
  {
    if(offset > this->count_ || records > this->count_ - offset)
    {
      throw std::out_of_range("AdaptiveRecordStore::readAt(): record range is out of bounds");
    }
    if(records == 0) { return; }
    if(target == nullptr)
    {
      throw std::invalid_argument("AdaptiveRecordStore::readAt(): null target");
    }

    const size_type bytes = this->bytesFor(records);
    if(this->resident())
    {
      std::copy_n(this->data_ + offset, records, target);
    }
    else
    {
      adaptive_record_store_detail::readRawRecords(this->descriptor_, target,
        bytes, this->bytesFor(offset), this->filename_);
    }
  }

  Element at(size_type offset) const
  {
    Element value;
    this->readAt(offset, &value, 1);
    return value;
  }

  // Materialize to a private temporary only when an append must spill or a
  // caller explicitly asks for it. The resulting file is deleted by clear()
  // or destruction unless materialize(filename) was used instead.
  const std::string& materialize()
  {
    if(this->spilled()) { return this->filename_; }
    std::string filename = TempFile::getName("gcsa_adaptive_records");
    try
    {
      this->materialize(filename, FileOwnership::OWNED_TEMPORARY);
    }
    catch(...)
    {
      TempFile::remove(filename); throw;
    }
    return this->filename_;
  }

  // The named path remains under the caller's lifecycle control. This is the
  // explicit route for a process consumer or an intentional checkpoint.
  void materialize(const std::string& filename)
  {
    this->materialize(filename, FileOwnership::BORROWED);
  }

  // Transfers deletion responsibility for the materialized filename to the
  // caller while keeping the store readable through that path.
  std::string releaseMaterializedFile()
  {
    if(this->resident())
    {
      throw std::logic_error("AdaptiveRecordStore::releaseMaterializedFile(): store is resident");
    }
    this->ownership_ = FileOwnership::BORROWED;
    return this->filename_;
  }

  const std::string& filename() const
  {
    if(this->resident())
    {
      throw std::logic_error("AdaptiveRecordStore::filename(): store is resident");
    }
    return this->filename_;
  }

  void clear() noexcept
  {
    if(this->descriptor_ >= 0)
    {
      adaptive_record_store_detail::closeRawRecordFile(this->descriptor_);
      this->descriptor_ = -1;
    }
    if(this->ownership_ == FileOwnership::OWNED_TEMPORARY)
    {
      TempFile::remove(this->filename_);
    }
    this->filename_.clear();
    this->ownership_ = FileOwnership::NONE;
    this->destroyResident();
  }

private:
  typedef std::allocator<Element> allocator_type;

  MemoryBudget* memory_;
  std::string task_;
  allocator_type allocator_;
  Element* data_;
  size_type count_, capacity_;
  int descriptor_;
  std::string filename_;
  FileOwnership ownership_;
  MemoryBudget::Reservation resident_memory_;

  static size_type bytesFor(size_type records)
  {
    if(records > std::numeric_limits<size_type>::max() / sizeof(Element))
    {
      throw std::length_error("AdaptiveRecordStore: record bytes overflow");
    }
    return records * sizeof(Element);
  }

  void requireResident(const char* operation) const
  {
    if(this->spilled())
    {
      throw std::logic_error(std::string("AdaptiveRecordStore::") + operation +
        "(): store has spilled");
    }
  }

  bool growResident(size_type records)
  {
    const size_type bytes = this->bytesFor(records);
    MemoryBudget::Reservation next = this->memory_->tryReserve(bytes, this->task_);
    if(!next) { return false; }

    Element* replacement = nullptr;
    size_type constructed = 0;
    try
    {
      replacement = std::allocator_traits<allocator_type>::allocate(this->allocator_, records);
      for(; constructed < this->count_; constructed++)
      {
        std::allocator_traits<allocator_type>::construct(this->allocator_, replacement + constructed,
          std::move_if_noexcept(this->data_[constructed]));
      }
    }
    catch(...)
    {
      while(constructed > 0)
      {
        constructed--;
        std::allocator_traits<allocator_type>::destroy(this->allocator_, replacement + constructed);
      }
      if(replacement != nullptr)
      {
        std::allocator_traits<allocator_type>::deallocate(this->allocator_, replacement, records);
      }
      throw;
    }

    const size_type old_count = this->count_;
    this->destroyResident();
    this->data_ = replacement;
    this->count_ = old_count;
    this->capacity_ = records;
    this->resident_memory_ = std::move(next);
    return true;
  }

  template<class Value>
  void appendValue(Value&& value)
  {
    if(this->count_ == std::numeric_limits<size_type>::max())
    {
      throw std::length_error("AdaptiveRecordStore::append(): record count overflows");
    }
    this->ensureAppendCapacity(this->count_ + 1);

    if(this->spilled())
    {
      adaptive_record_store_detail::writeRawRecords(this->descriptor_, &value,
        sizeof(Element), this->bytesFor(this->count_), this->filename_);
    }
    else
    {
      std::allocator_traits<allocator_type>::construct(this->allocator_,
        this->data_ + this->count_, std::forward<Value>(value));
    }
    this->count_++;
  }

  void ensureAppendCapacity(size_type required)
  {
    if(!this->resident() || required <= this->capacity_) { return; }
    size_type next = (this->capacity_ == 0 ? 1 : this->capacity_);
    while(next < required)
    {
      if(next > std::numeric_limits<size_type>::max() / 2)
      {
        next = required;
        break;
      }
      next *= 2;
    }
    if(!(this->tryReserveRecords(next))) { this->materialize(); }
  }

  void materialize(const std::string& filename, FileOwnership ownership)
  {
    if(filename.empty())
    {
      throw std::invalid_argument("AdaptiveRecordStore::materialize(): empty filename");
    }
    if(this->spilled())
    {
      if(this->filename_ == filename)
      {
        this->ownership_ = ownership;
        return;
      }

      int descriptor = adaptive_record_store_detail::openRawRecordFile(filename);
      try
      {
        adaptive_record_store_detail::copyRawRecords(this->descriptor_, descriptor,
          this->bytesFor(this->count_), this->filename_, filename);
      }
      catch(...)
      {
        adaptive_record_store_detail::closeRawRecordFile(descriptor); throw;
      }

      std::string old_filename = this->filename_;
      const FileOwnership old_ownership = this->ownership_;
      adaptive_record_store_detail::closeRawRecordFile(this->descriptor_);
      if(old_ownership == FileOwnership::OWNED_TEMPORARY)
      {
        TempFile::remove(old_filename);
      }
      this->descriptor_ = descriptor;
      this->filename_ = filename;
      this->ownership_ = ownership;
      return;
    }

    int descriptor = adaptive_record_store_detail::openRawRecordFile(filename);
    try
    {
      if(this->count_ > 0)
      {
        adaptive_record_store_detail::writeRawRecords(descriptor, this->data_,
          this->bytesFor(this->count_), 0, filename);
      }
    }
    catch(...)
    {
      adaptive_record_store_detail::closeRawRecordFile(descriptor); throw;
    }

    const size_type old_count = this->count_;
    this->destroyResident();
    this->descriptor_ = descriptor;
    this->count_ = old_count;
    this->filename_ = filename;
    this->ownership_ = ownership;
  }

  void destroyResident() noexcept
  {
    if(this->data_ != nullptr)
    {
      for(size_type i = this->count_; i > 0; i--)
      {
        std::allocator_traits<allocator_type>::destroy(this->allocator_, this->data_ + i - 1);
      }
      std::allocator_traits<allocator_type>::deallocate(this->allocator_, this->data_, this->capacity_);
    }
    this->data_ = nullptr;
    this->count_ = 0;
    this->capacity_ = 0;
    this->resident_memory_ = MemoryBudget::Reservation();
  }
};

} // namespace gcsa

#endif // GCSA_ADAPTIVE_RECORD_STORE_H
