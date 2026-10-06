#ifndef GCSA_RESOURCES_H
#define GCSA_RESOURCES_H

#include <cstdint>
#include <condition_variable>
#include <mutex>

namespace gcsa
{

// A process-local, deliberately conservative memory admission controller.
class MemoryBudget
{
public:
  struct Stats
  {
    uint64_t current, maximum;
  };
  class Reservation
  {
  public:
    Reservation();
    Reservation(Reservation&& other);
    Reservation& operator=(Reservation&& other);
    ~Reservation();
    uint64_t bytes() const { return amount; }
    explicit operator bool() const { return owner != nullptr; }
  private:
    friend class MemoryBudget;
    Reservation(MemoryBudget* budget, uint64_t bytes);
    Reservation(const Reservation&);
    Reservation& operator=(const Reservation&);
    MemoryBudget* owner;
    uint64_t amount;
  };
  // ceiling includes the margin: no more than ceiling - safety_margin is reservable.
  MemoryBudget(uint64_t ceiling, uint64_t safety_margin = 0);
  Reservation reserve(uint64_t bytes);
  Stats stats() const;
  uint64_t available() const;
private:
  friend class Reservation;
  void release(uint64_t bytes);
  uint64_t ceiling_, margin_, current_, maximum_;
  mutable std::mutex mutex_;
  std::condition_variable changed_;
};

} // namespace gcsa
#endif
