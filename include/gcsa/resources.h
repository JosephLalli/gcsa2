#ifndef GCSA_RESOURCES_H
#define GCSA_RESOURCES_H

#include <cstdint>
#include <condition_variable>
#include <map>
#include <mutex>
#include <string>

namespace gcsa {

// A process-local, deliberately conservative memory admission controller.
class MemoryBudget {
public:
  struct Stats
  {
    uint64_t current, maximum, largest;
    std::string largest_task;
    std::map<std::string, uint64_t> active;
  };
  class Reservation {
  public:
    Reservation();
    Reservation(Reservation&& other);
    Reservation& operator=(Reservation&& other);
    ~Reservation();
    uint64_t bytes() const { return amount; }
    explicit operator bool() const { return owner != 0; }
  private:
    friend class MemoryBudget;
    Reservation(MemoryBudget* budget, const std::string& name, uint64_t bytes);
    Reservation(const Reservation&);
    Reservation& operator=(const Reservation&);
    MemoryBudget* owner; std::string task; uint64_t amount;
  };
  // ceiling includes the margin: no more than ceiling - safety_margin is reservable.
  MemoryBudget(uint64_t ceiling, uint64_t safety_margin = 0);
  Reservation reserve(uint64_t bytes, const std::string& task);
  Stats stats() const;
  uint64_t available() const;
private:
  friend class Reservation;
  void release(const std::string& task, uint64_t bytes);
  uint64_t ceiling_, margin_, current_, maximum_, largest_;
  std::string largest_task_;
  mutable std::mutex mutex_; std::condition_variable changed_;
  std::map<std::string, uint64_t> active_;
};

class DiskBudget {
public:
  DiskBudget(const std::string& path, uint64_t configured_limit, uint64_t free_safety_margin = 0);
  bool can_reserve(uint64_t bytes, std::string* reason = 0) const;
  void require(uint64_t bytes) const;
  uint64_t available() const;
private:
  std::string path_; uint64_t limit_, margin_;
};

} // namespace gcsa
#endif
