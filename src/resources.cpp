#include <gcsa/resources.h>
#include <algorithm>
#include <sys/statvfs.h>
#include <sys/stat.h>
#include <dirent.h>
#include <stdexcept>
#include <sstream>
#include <utility>

namespace gcsa {
MemoryBudget::Reservation::Reservation() : owner(0), amount(0) {}
MemoryBudget::Reservation::Reservation(MemoryBudget* b, const std::string& n, uint64_t a) : owner(b), task(n), amount(a) {}
MemoryBudget::Reservation::Reservation(Reservation&& o) : owner(o.owner), task(o.task), amount(o.amount) { o.owner = 0; o.amount = 0; }
MemoryBudget::Reservation& MemoryBudget::Reservation::operator=(Reservation&& o) { if(this != &o) { if(owner) owner->release(task, amount); owner=o.owner; task=o.task; amount=o.amount; o.owner=0; o.amount=0; } return *this; }
MemoryBudget::Reservation::~Reservation() { if(owner) owner->release(task, amount); }
MemoryBudget::MemoryBudget(uint64_t c, uint64_t m) : ceiling_(c), margin_(m), current_(0), maximum_(0), largest_(0) { if(m > c) throw std::invalid_argument("MemoryBudget safety margin exceeds ceiling"); }
MemoryBudget::Reservation MemoryBudget::reserve(uint64_t bytes, const std::string& task) { if(bytes > ceiling_ - margin_) throw std::runtime_error("memory reservation exceeds budget after safety margin"); std::unique_lock<std::mutex> lock(mutex_); changed_.wait(lock, [this, bytes] { return bytes <= ceiling_ - margin_ - current_; }); current_ += bytes; if(current_ > maximum_) maximum_=current_; if(bytes > largest_) { largest_=bytes; largest_task_=task; } active_[task] += bytes; return Reservation(this, task, bytes); }
void MemoryBudget::release(const std::string& task, uint64_t bytes) { std::lock_guard<std::mutex> lock(mutex_); current_ -= bytes; std::map<std::string,uint64_t>::iterator i=active_.find(task); if(i != active_.end() && (i->second -= bytes) == 0) active_.erase(i); changed_.notify_all(); }
MemoryBudget::Stats MemoryBudget::stats() const { std::lock_guard<std::mutex> lock(mutex_); Stats s={current_,maximum_,largest_,largest_task_,active_}; return s; }
uint64_t MemoryBudget::available() const { std::lock_guard<std::mutex> lock(mutex_); return ceiling_ - margin_ - current_; }
DiskBudget::DiskBudget(const std::string& p,uint64_t l,uint64_t m):path_(p),limit_(l),margin_(m) {}
namespace {
uint64_t live_bytes(const std::string& path)
{
  DIR* directory = opendir(path.c_str());
  if(directory == nullptr)
  {
    throw std::runtime_error("cannot open workspace for disk accounting: " + path);
  }
  uint64_t total = 0;
  for(dirent* entry; (entry = readdir(directory)) != nullptr;)
  {
    if(entry->d_name[0] == '.') { continue; }
    const std::string child = path + "/" + entry->d_name;
    struct stat info;
    if(lstat(child.c_str(), &info) != 0)
    {
      closedir(directory);
      throw std::runtime_error("cannot stat workspace entry: " + child);
    }
    if(S_ISREG(info.st_mode))
    {
      // Adopted checkpoints and restored PathGraph shards are hard links to
      // the same immutable inode. Charge each directory entry only its share
      // of allocated blocks, so link-based checkpointing does not appear to
      // consume the payload again under the configured disk limit.
      const uint64_t links = std::max<uint64_t>(1, info.st_nlink);
      total += (static_cast<uint64_t>(info.st_blocks) * 512) / links;
    }
    else if(S_ISDIR(info.st_mode)) { total += live_bytes(child); }
  }
  closedir(directory);
  return total;
}
}
uint64_t DiskBudget::available() const { struct statvfs v; if(statvfs(path_.c_str(), &v) != 0) throw std::runtime_error("cannot inspect free space for " + path_); uint64_t free_bytes=uint64_t(v.f_bavail)*v.f_frsize, live=live_bytes(path_); if(free_bytes <= margin_ || live >= limit_) return 0; uint64_t safe=free_bytes-margin_, remaining=limit_-live; return remaining < safe ? remaining : safe; }
bool DiskBudget::can_reserve(uint64_t bytes,std::string* reason) const { uint64_t a=available(); if(bytes <= a) return true; if(reason) { std::ostringstream out; out << "disk reservation of " << bytes << " bytes exceeds available safe budget " << a; *reason=out.str(); } return false; }
void DiskBudget::require(uint64_t bytes) const { std::string why; if(!can_reserve(bytes,&why)) throw std::runtime_error(why); }
} // namespace gcsa
