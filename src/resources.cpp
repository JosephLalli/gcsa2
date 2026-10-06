#include <gcsa/resources.h>
#include <stdexcept>

namespace gcsa
{

MemoryBudget::Reservation::Reservation() : owner(nullptr), amount(0)
{
}

MemoryBudget::Reservation::Reservation(MemoryBudget* budget, uint64_t bytes) :
  owner(budget), amount(bytes)
{
}

MemoryBudget::Reservation::Reservation(Reservation&& other) :
  owner(other.owner), amount(other.amount)
{
  other.owner = nullptr;
  other.amount = 0;
}

MemoryBudget::Reservation&
MemoryBudget::Reservation::operator=(Reservation&& other)
{
  if(this != &other)
  {
    if(this->owner != nullptr) { this->owner->release(this->amount); }
    this->owner = other.owner;
    this->amount = other.amount;
    other.owner = nullptr;
    other.amount = 0;
  }
  return *this;
}

MemoryBudget::Reservation::~Reservation()
{
  if(this->owner != nullptr) { this->owner->release(this->amount); }
}

MemoryBudget::MemoryBudget(uint64_t ceiling, uint64_t safety_margin) :
  ceiling_(ceiling), margin_(safety_margin), current_(0), maximum_(0)
{
  if(safety_margin > ceiling)
  {
    throw std::invalid_argument("MemoryBudget safety margin exceeds ceiling");
  }
}

MemoryBudget::Reservation
MemoryBudget::reserve(uint64_t bytes)
{
  if(bytes > this->ceiling_ - this->margin_)
  {
    throw std::runtime_error("memory reservation exceeds budget after safety margin");
  }
  std::unique_lock<std::mutex> lock(this->mutex_);
  this->changed_.wait(lock, [this, bytes]
  {
    return bytes <= this->ceiling_ - this->margin_ - this->current_;
  });
  this->current_ += bytes;
  if(this->current_ > this->maximum_) { this->maximum_ = this->current_; }
  return Reservation(this, bytes);
}

void
MemoryBudget::release(uint64_t bytes)
{
  std::lock_guard<std::mutex> lock(this->mutex_);
  this->current_ -= bytes;
  this->changed_.notify_all();
}

MemoryBudget::Stats
MemoryBudget::stats() const
{
  std::lock_guard<std::mutex> lock(this->mutex_);
  return { this->current_, this->maximum_ };
}

uint64_t
MemoryBudget::available() const
{
  std::lock_guard<std::mutex> lock(this->mutex_);
  return this->ceiling_ - this->margin_ - this->current_;
}

} // namespace gcsa
