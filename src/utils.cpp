/*
  Copyright (c) 2018, 2019 Jouni Siren
  Copyright (c) 2015, 2016, 2017 Genome Research Ltd.

  Author: Jouni Siren <jouni.siren@iki.fi>

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/

#include <gcsa/utils.h>

#include <atomic>
#include <cerrno>
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cmath>
#include <deque>
#include <iomanip>
#include <limits>
#include <new>
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <thread>

#include <fcntl.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <unistd.h>

#include <gcsa/internal.h>

namespace gcsa
{

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr size_type Verbosity::SILENT;
constexpr size_type Verbosity::BASIC;
constexpr size_type Verbosity::EXTENDED;
constexpr size_type Verbosity::DEFAULT;
constexpr size_type Verbosity::FULL;

constexpr size_type Version::MAJOR_VERSION;
constexpr size_type Version::MINOR_VERSION;
constexpr size_type Version::PATCH_VERSION;
constexpr size_type Version::GCSA_VERSION;
constexpr size_type Version::LCP_VERSION;

//------------------------------------------------------------------------------

// Other class variables.

size_type Verbosity::level = Verbosity::DEFAULT;

//------------------------------------------------------------------------------

void
Verbosity::set(size_type new_level)
{
  level = Range::bound(new_level, SILENT, FULL);
}

std::string
Verbosity::levelName()
{
  switch(level)
  {
    case SILENT:
      return "silent"; break;
    case BASIC:
      return "basic"; break;
    case EXTENDED:
      return "extended"; break;
    case FULL:
      return "full"; break;
  }
  return "unknown";
}

//------------------------------------------------------------------------------

std::string
Version::str(bool verbose)
{
  std::ostringstream ss;
  if(verbose) { ss << "GCSA2 version "; }
  else { ss << "v"; }
  ss << MAJOR_VERSION << "." << MINOR_VERSION << "." << PATCH_VERSION;
  if(verbose) { ss << " (GCSA version " << GCSA_VERSION << ", LCP version " << LCP_VERSION << ")"; }
  return ss.str();
}

void
Version::print(std::ostream& out, const std::string& tool_name, bool verbose, size_type new_lines)
{
  out << tool_name;
  if(verbose) { out << std::endl; }
  else { out << " "; }
  out << str(verbose);
  for(size_type i = 0; i < new_lines; i++) { out << std::endl; }
}

//------------------------------------------------------------------------------

void
printHeader(const std::string& header, size_type indent)
{
  std::string padding;
  if(header.length() + 1 < indent) { padding = std::string(indent - 1 - header.length(), ' '); }
  std::cout << header << ":" << padding;
}

void
printTime(const std::string& header, size_type queries, double seconds, size_type indent)
{
  printHeader(header, indent);
  std::cout << queries << " queries in " << seconds << " seconds ("
            << inMicroseconds(seconds / queries) << " µs/query)" << std::endl;
}

//------------------------------------------------------------------------------

double
readTimer()
{
  return omp_get_wtime();
}

size_type
memoryUsage()
{
  rusage usage;
  getrusage(RUSAGE_SELF, &usage);
#if defined(__APPLE__) && defined(__MACH__)
  return usage.ru_maxrss;
#else
  return KILOBYTE * usage.ru_maxrss;
#endif
}

size_type
readVolume()
{
  return DiskIO::read_volume;
}

size_type
writeVolume()
{
  return DiskIO::write_volume;
}

size_type
parseBytes(const std::string& value)
{
  if(value.empty() || value[0] == '-')
  {
    throw std::invalid_argument("empty or negative byte count: " + value);
  }

  size_t parsed = 0;
  long double number = std::stold(value, &parsed);
  if(!std::isfinite(number) || number < 0.0L)
  {
    throw std::invalid_argument("non-finite or negative byte count: " + value);
  }

  std::string suffix = value.substr(parsed);
  std::transform(suffix.begin(), suffix.end(), suffix.begin(),
    [](unsigned char c) { return std::toupper(c); });
  if(suffix.length() >= 2 && suffix.substr(suffix.length() - 2) == "IB")
  {
    suffix.resize(suffix.length() - 2);
  }
  else if(!(suffix.empty()) && suffix.back() == 'B')
  {
    suffix.pop_back();
  }

  size_type multiplier = 1;
  if(suffix == "K") { multiplier = KILOBYTE; }
  else if(suffix == "M") { multiplier = MEGABYTE; }
  else if(suffix == "G") { multiplier = GIGABYTE; }
  else if(suffix == "T") { multiplier = GIGABYTE * KILOBYTE; }
  else if(suffix == "P") { multiplier = GIGABYTE * MEGABYTE; }
  else if(!suffix.empty()) { throw std::invalid_argument("invalid byte suffix: " + value); }

  long double result = number * static_cast<long double>(multiplier);
  if(result > static_cast<long double>(std::numeric_limits<size_type>::max()))
  {
    throw std::out_of_range("byte count is too large: " + value);
  }
  return static_cast<size_type>(result);
}

std::string
formatBytes(size_type bytes)
{
  const char* suffixes[] = { "B", "KiB", "MiB", "GiB", "TiB", "PiB" };
  long double value = bytes;
  size_type suffix = 0;
  while(value >= 1024.0L && suffix + 1 < sizeof(suffixes) / sizeof(suffixes[0]))
  {
    value /= 1024.0L; suffix++;
  }
  std::ostringstream out;
  out << std::fixed << std::setprecision(value < 10.0L && suffix > 0 ? 2 : 1)
      << static_cast<double>(value) << ' ' << suffixes[suffix];
  return out.str();
}

//------------------------------------------------------------------------------

namespace TempFile
{
  size_type counter = 0;

  const std::string DEFAULT_TEMP_DIR = ".";
  std::string temp_dir = DEFAULT_TEMP_DIR;

  // Guards counter, temp_dir and the handler's set: the external join names
  // temporaries from several threads at once. Defined before the handler so
  // it is still alive when the handler's destructor runs at exit.
  std::mutex mutex;

  // By storing the filenames in a static object, we can delete the remaining
  // temporary files when std::exit() is called.
  struct Handler
  {
    std::set<std::string> filenames;
    ~Handler()
    {
      std::lock_guard<std::mutex> lock(mutex);
      for(auto& filename : this->filenames)
      {
        std::remove(filename.c_str());
      }
    }
  } handler;

  void
  setDirectory(const std::string& directory)
  {
    std::lock_guard<std::mutex> lock(mutex);
    if(directory.empty()) { temp_dir = DEFAULT_TEMP_DIR; }
    else if(directory[directory.length() - 1] != '/') { temp_dir = directory; }
    else { temp_dir = directory.substr(0, directory.length() - 1); }
  }

  std::string
  getName(const std::string& name_part)
  {
    char hostname[32];
    gethostname(hostname, 32); hostname[31] = 0;

    std::lock_guard<std::mutex> lock(mutex);
    std::string filename = temp_dir + '/' + name_part + '_'
      + std::string(hostname) + '_'
      + sdsl::util::to_string(sdsl::util::pid()) + '_'
      + sdsl::util::to_string(counter);
    handler.filenames.insert(filename);
    counter++;

    return filename;
  }

  void
  remove(std::string& filename)
  {
    if(!(filename.empty()))
    {
      std::remove(filename.c_str());
      std::lock_guard<std::mutex> lock(mutex);
      handler.filenames.erase(filename);
      filename.clear();
    }
  }

  void
  forget() {
    std::lock_guard<std::mutex> lock(mutex);
    handler.filenames.clear();
    counter = 0;
  }
} // namespace TempFile

size_type
readRows(const std::string& filename, std::vector<std::string>& rows, bool skip_empty_rows)
{
  std::ifstream input(filename.c_str(), std::ios_base::binary);
  if(!input)
  {
    std::cerr << "readRows(): Cannot open input file " << filename << std::endl;
    return 0;
  }

  size_type chars = 0;
  while(input)
  {
    std::string buf;
    std::getline(input, buf);
    if(skip_empty_rows && buf.empty()) { continue; }
    rows.push_back(buf);
    chars += buf.length();
  }

  input.close();
  return chars;
}

size_type
fileSize(std::ifstream& file)
{
  std::streamoff curr = file.tellg();

  file.seekg(0, std::ios::end);
  std::streamoff size = file.tellg();
  file.seekg(0, std::ios::beg);
  size -= file.tellg();

  file.seekg(curr, std::ios::beg);
  return size;
}

size_type
fileSize(std::ofstream& file)
{
  std::streamoff curr = file.tellp();

  file.seekp(0, std::ios::end);
  std::streamoff size = file.tellp();
  file.seekp(0, std::ios::beg);
  size -= file.tellp();

  file.seekp(curr, std::ios::beg);
  return size;
}

//------------------------------------------------------------------------------

size_type
getChunkSize(size_type n, size_type min_size)
{
  size_type threads = omp_get_max_threads();
  size_type chunks = std::max(threads, (size_type)8) * threads;
  size_type chunk_size = n / chunks;
  return std::max(chunk_size, min_size);
}

//------------------------------------------------------------------------------

namespace
{

// How many advance() calls to take before consulting the clock again. Small
// enough that a slow phase still reports on time, large enough that a fast one
// does not call readTimer() per record.
constexpr size_type PROGRESS_CHECK_STRIDE = 4096;

std::string
formatDuration(double seconds)
{
  if(seconds < 0.0 || !(seconds < 1.0e9)) { return "unknown"; }
  size_type total = static_cast<size_type>(seconds + 0.5);
  size_type hours = total / 3600, minutes = (total % 3600) / 60, secs = total % 60;
  std::ostringstream result;
  if(hours > 0) { result << hours << "h"; }
  if(hours > 0 || minutes > 0)
  {
    if(hours > 0 && minutes < 10) { result << "0"; }
    result << minutes << "m";
    if(secs < 10) { result << "0"; }
  }
  result << secs << "s";
  return result.str();
}

} // anonymous namespace

ProgressReporter::ProgressReporter(const std::string& phase_name,
  size_type total_units, const std::string& unit_name, double interval_seconds) :
  name(phase_name), unit(unit_name), total(total_units), processed(0),
  next_check(PROGRESS_CHECK_STRIDE), start_time(readTimer()),
  last_report(start_time), interval(interval_seconds), done(false)
{
}

void
ProgressReporter::advance(size_type units)
{
  this->processed += units;
  if(this->processed < this->next_check) { return; }
  this->next_check = this->processed + PROGRESS_CHECK_STRIDE;
  if(this->interval <= 0.0 || Verbosity::level < Verbosity::BASIC) { return; }
  if(readTimer() - this->last_report < this->interval) { return; }
  this->report(false);
}

void
ProgressReporter::finish()
{
  if(this->done) { return; }
  this->done = true;
  if(this->interval <= 0.0 || Verbosity::level < Verbosity::BASIC) { return; }
  this->report(true);
}

void
ProgressReporter::report(bool final_report)
{
  double now = readTimer(), elapsed = now - this->start_time;
  this->last_report = now;
  double rate = (elapsed > 0.0 ? this->processed / elapsed : 0.0);

  std::cerr << "  " << this->name << ": " << this->processed;
  if(this->total > 0)
  {
    std::cerr << " / " << this->total << " " << this->unit << " ("
              << std::fixed << std::setprecision(1)
              << (100.0 * this->processed / this->total) << "%)";
  }
  else { std::cerr << " " << this->unit; }
  std::cerr << ", " << formatDuration(elapsed) << " elapsed";
  if(rate > 0.0)
  {
    std::cerr << ", " << static_cast<size_type>(rate) << " " << this->unit << "/s";
    // An estimate from the mean rate so far. A phase whose cost per unit is
    // not uniform will drift, so this is a progress indication and not a
    // schedule.
    if(!final_report && this->total > this->processed)
    {
      std::cerr << ", ~" << formatDuration((this->total - this->processed) / rate)
                << " left";
    }
  }
  std::cerr << std::endl;
}

//------------------------------------------------------------------------------

struct ReadAhead::Source
{
  int descriptor = -1;
  bool owns = false;
  std::uint64_t size = 0;
  std::mutex mtx;
  std::condition_variable idle;
  std::size_t in_flight = 0;
  bool cancelled = false;

  ~Source() { if(this->owns && this->descriptor >= 0) { ::close(this->descriptor); } }
};

namespace
{

// A non-negative integer from the environment, with an optional binary K/M/G suffix.
std::uint64_t
environmentSize(const char* name, std::uint64_t fallback)
{
  const char* value = std::getenv(name);
  if(value == nullptr || *value == '\0') { return fallback; }
  char* end = nullptr;
  errno = 0;
  unsigned long long parsed = std::strtoull(value, &end, 10);
  std::uint64_t scale = 1;
  if(end != value && *end != '\0' && end[1] == '\0')
  {
    switch(std::toupper(static_cast<unsigned char>(*end)))
    {
      case 'K': scale = std::uint64_t(1) << 10; end++; break;
      case 'M': scale = std::uint64_t(1) << 20; end++; break;
      case 'G': scale = std::uint64_t(1) << 30; end++; break;
    }
  }
  if(end == value || *end != '\0' || errno == ERANGE ||
     parsed > std::numeric_limits<std::uint64_t>::max() / scale)
  {
    throw std::invalid_argument(std::string(name) + " must be a size such as 64M, not '" + value + "'");
  }
  return parsed * scale;
}

// Reads queued pieces into a scratch buffer, which leaves them in the page cache.
class ReadAheadPool
{
public:
  static ReadAheadPool& instance()
  {
    // Never destroyed: its detached threads wait on the queue for the life of the process.
    static ReadAheadPool* pool = new ReadAheadPool();
    return *pool;
  }

  void submit(const std::shared_ptr<ReadAhead::Source>& source, std::uint64_t from, std::uint64_t to)
  {
    {
      std::lock_guard<std::mutex> lock(this->mtx);
      for(std::uint64_t offset = from; offset < to; offset += ReadAhead::PIECE_BYTES)
      {
        this->pieces.emplace_back(source, offset);
      }
    }
    this->ready.notify_all();
  }

private:
  ReadAheadPool()
  {
    const std::uint64_t threads = environmentSize("GCSA_IO_READAHEAD_THREADS", 16);
    if(threads == 0) { throw std::invalid_argument("GCSA_IO_READAHEAD_THREADS must be positive"); }
    for(std::uint64_t i = 0; i < threads; i++) { std::thread(&ReadAheadPool::work, this).detach(); }
  }

  void work()
  {
    std::vector<char> scratch(ReadAhead::PIECE_BYTES);
    while(true)
    {
      std::pair<std::shared_ptr<ReadAhead::Source>, std::uint64_t> piece;
      {
        std::unique_lock<std::mutex> lock(this->mtx);
        this->ready.wait(lock, [this]() { return !this->pieces.empty(); });
        piece = std::move(this->pieces.front());
        this->pieces.pop_front();
      }
      ReadAhead::Source& source = *(piece.first);
      {
        std::lock_guard<std::mutex> lock(source.mtx);
        if(source.cancelled) { continue; }
        source.in_flight++;
      }
      // Advisory: after an error or a short read, the reader reads the range itself.
      static_cast<void>(::pread(source.descriptor, scratch.data(), scratch.size(),
        static_cast<off_t>(piece.second)));
      {
        std::lock_guard<std::mutex> lock(source.mtx);
        source.in_flight--;
      }
      source.idle.notify_all();
    }
  }

  std::mutex mtx;
  std::condition_variable ready;
  std::deque<std::pair<std::shared_ptr<ReadAhead::Source>, std::uint64_t>> pieces;
};

} // anonymous namespace

std::uint64_t
ReadAhead::windowBytes()
{
  static const std::uint64_t window = environmentSize("GCSA_IO_READAHEAD_BYTES", 0);
  return window;
}

ReadAhead::ReadAhead(const std::string& filename) :
  source(), requested(0)
{
  if(windowBytes() == 0) { return; }
  // Advisory: if this descriptor cannot be opened the reader still reads the
  // file itself, without read-ahead.
  int descriptor = ::open(filename.c_str(), O_RDONLY | O_CLOEXEC);
  if(descriptor < 0) { return; }
  struct stat status;
  if(::fstat(descriptor, &status) != 0 || status.st_size < 0) { ::close(descriptor); return; }
  this->source = std::make_shared<Source>();
  this->source->descriptor = descriptor;
  this->source->owns = true;
  this->source->size = static_cast<std::uint64_t>(status.st_size);
}

ReadAhead::ReadAhead(int descriptor) :
  source(), requested(0)
{
  if(windowBytes() == 0 || descriptor < 0) { return; }
  struct stat status;
  if(::fstat(descriptor, &status) != 0 || status.st_size < 0) { return; }
  this->source = std::make_shared<Source>();
  this->source->descriptor = descriptor;
  this->source->size = static_cast<std::uint64_t>(status.st_size);
}

ReadAhead::~ReadAhead()
{
  this->release();
}

ReadAhead::ReadAhead(ReadAhead&& another) noexcept :
  source(std::move(another.source)), requested(another.requested)
{
}

ReadAhead&
ReadAhead::operator=(ReadAhead&& another) noexcept
{
  if(this != &another)
  {
    this->release();
    this->source = std::move(another.source);
    this->requested = another.requested;
  }
  return *this;
}

void
ReadAhead::release()
{
  if(!(this->source)) { return; }
  {
    std::unique_lock<std::mutex> lock(this->source->mtx);
    this->source->cancelled = true;
    this->source->idle.wait(lock, [this]() { return this->source->in_flight == 0; });
  }
  this->source.reset();
}

void
ReadAhead::advance(std::uint64_t offset)
{
  if(!(this->source)) { return; }
  const std::uint64_t target = std::min(this->source->size, offset + windowBytes());
  if(this->requested < offset) { this->requested = offset - offset % PIECE_BYTES; }
  // Queue whole pieces, so a reader advancing a record at a time does not queue
  // a job per call; the file's last piece may be partial.
  if(target <= this->requested ||
     (target - this->requested < PIECE_BYTES && target < this->source->size)) { return; }
  ReadAheadPool::instance().submit(this->source, this->requested, target);
  this->requested = target;
}

//------------------------------------------------------------------------------

namespace
{

constexpr std::size_t DIRECT_ALIGNMENT = 4096;
constexpr std::size_t DIRECT_BLOCK_BYTES = std::size_t(4) << 20;

struct AlignedBlock
{
  char* data;

  AlignedBlock() : data(nullptr)
  {
    void* memory = nullptr;
    if(::posix_memalign(&memory, DIRECT_ALIGNMENT, DIRECT_BLOCK_BYTES) != 0) { throw std::bad_alloc(); }
    this->data = static_cast<char*>(memory);
  }
  ~AlignedBlock() { std::free(this->data); }

  AlignedBlock(const AlignedBlock&) = delete;
  AlignedBlock& operator=(const AlignedBlock&) = delete;
};

// Writes all of `bytes` at `offset`; returns 0 or an errno value.
int
pwriteAll(int descriptor, const char* data, std::size_t bytes, std::uint64_t offset)
{
  std::size_t done = 0;
  while(done < bytes)
  {
    ssize_t result = ::pwrite(descriptor, data + done, bytes - done,
      static_cast<off_t>(offset + done));
    if(result < 0 && errno == EINTR) { continue; }
    if(result <= 0) { return (result < 0 ? errno : EIO); }
    done += static_cast<std::size_t>(result);
  }
  return 0;
}

// What a direct stream shares with its blocks in flight.
struct DirectWriteState
{
  int descriptor = -1;
  std::mutex mtx;
  std::condition_variable returned;
  std::vector<std::unique_ptr<AlignedBlock>> free_blocks;
  std::size_t in_flight = 0;
  int error = 0;  // first failed write's errno
};

struct DirectWriteJob
{
  std::shared_ptr<DirectWriteState> state;
  std::unique_ptr<AlignedBlock> block;
  std::uint64_t offset;
};

class DirectWritePool
{
public:
  static DirectWritePool& instance()
  {
    // Never destroyed: its detached threads wait on the queue for the life of the process.
    static DirectWritePool* pool = new DirectWritePool();
    return *pool;
  }

  void submit(DirectWriteJob job)
  {
    {
      std::lock_guard<std::mutex> lock(this->mtx);
      this->jobs.push_back(std::move(job));
    }
    this->ready.notify_one();
  }

private:
  DirectWritePool()
  {
    const std::uint64_t threads = environmentSize("GCSA_IO_DIRECT_WRITE_THREADS", 16);
    if(threads == 0) { throw std::invalid_argument("GCSA_IO_DIRECT_WRITE_THREADS must be positive"); }
    for(std::uint64_t i = 0; i < threads; i++) { std::thread(&DirectWritePool::work, this).detach(); }
  }

  void work()
  {
    while(true)
    {
      DirectWriteJob job;
      {
        std::unique_lock<std::mutex> lock(this->mtx);
        this->ready.wait(lock, [this]() { return !this->jobs.empty(); });
        job = std::move(this->jobs.front());
        this->jobs.pop_front();
      }
      int error = pwriteAll(job.state->descriptor, job.block->data, DIRECT_BLOCK_BYTES, job.offset);
      {
        std::lock_guard<std::mutex> lock(job.state->mtx);
        if(error != 0 && job.state->error == 0) { job.state->error = error; }
        job.state->free_blocks.push_back(std::move(job.block));
        job.state->in_flight--;
      }
      job.state->returned.notify_all();
    }
  }

  std::mutex mtx;
  std::condition_variable ready;
  std::deque<DirectWriteJob> jobs;
};

/*
  Put area-free stream buffer: every write goes through xsputn(), which fills
  the current aligned block and hands full blocks to the pool. The current block
  always starts at an aligned file offset.
*/
class DirectOutputBuffer : public std::streambuf
{
public:
  DirectOutputBuffer() : state(), current(), used(0), block_offset(0), open_(false) { }
  ~DirectOutputBuffer() override { if(this->open_) { this->close(); } }

  // Returns 0, or the errno of a failed open (EINVAL when O_DIRECT is unsupported).
  int open(const std::string& filename)
  {
    int descriptor = ::open(filename.c_str(),
      O_WRONLY | O_CREAT | O_TRUNC | O_DIRECT | O_CLOEXEC, 0666);
    if(descriptor < 0) { return errno; }
    this->state = std::make_shared<DirectWriteState>();
    this->state->descriptor = descriptor;
    const std::uint64_t depth = std::max(std::uint64_t(2),
      environmentSize("GCSA_IO_DIRECT_WRITE_DEPTH", 4));
    for(std::uint64_t i = 1; i < depth; i++)
    {
      this->state->free_blocks.emplace_back(new AlignedBlock());
    }
    this->current.reset(new AlignedBlock());
    this->used = 0; this->block_offset = 0; this->open_ = true;
    return 0;
  }

  bool is_open() const { return this->open_; }

  bool close()
  {
    if(!this->open_) { return false; }
    bool ok = (this->sync() == 0);
    ok = (::close(this->state->descriptor) == 0) && ok;
    this->open_ = false;
    return ok;
  }

protected:
  std::streamsize xsputn(const char* data, std::streamsize count) override
  {
    // After a failed write the current block is gone; the stream reports the short write.
    if(!this->open_ || !this->current || count < 0) { return 0; }
    std::streamsize done = 0;
    while(done < count)
    {
      const std::size_t take = std::min(DIRECT_BLOCK_BYTES - this->used,
        static_cast<std::size_t>(count - done));
      std::memcpy(this->current->data + this->used, data + done, take);
      this->used += take; done += static_cast<std::streamsize>(take);
      if(this->used == DIRECT_BLOCK_BYTES && !this->submitCurrent()) { return done; }
    }
    return done;
  }

  int_type overflow(int_type ch) override
  {
    if(traits_type::eq_int_type(ch, traits_type::eof())) { return traits_type::not_eof(ch); }
    char c = traits_type::to_char_type(ch);
    return (this->xsputn(&c, 1) == 1 ? ch : traits_type::eof());
  }

  // Writes every byte: waits for the blocks in flight, writes the current block
  // padded to the alignment, and truncates the file to its logical size.
  int sync() override
  {
    if(!this->open_) { return -1; }
    std::unique_lock<std::mutex> lock(this->state->mtx);
    this->state->returned.wait(lock, [this]() { return this->state->in_flight == 0; });
    if(this->state->error != 0) { return -1; }
    lock.unlock();
    if(this->used > 0)
    {
      const std::size_t padded = (this->used + DIRECT_ALIGNMENT - 1) / DIRECT_ALIGNMENT * DIRECT_ALIGNMENT;
      std::memset(this->current->data + this->used, 0, padded - this->used);
      if(pwriteAll(this->state->descriptor, this->current->data, padded, this->block_offset) != 0) { return -1; }
    }
    if(::ftruncate(this->state->descriptor, static_cast<off_t>(this->block_offset + this->used)) != 0) { return -1; }
    return 0;
  }

  pos_type seekoff(off_type offset, std::ios_base::seekdir direction, std::ios_base::openmode which) override
  {
    if(offset != 0 || direction != std::ios_base::cur || !(which & std::ios_base::out)) { return pos_type(off_type(-1)); }
    return pos_type(static_cast<off_type>(this->block_offset + this->used));
  }

  pos_type seekpos(pos_type, std::ios_base::openmode) override { return pos_type(off_type(-1)); }

private:
  // Hands the full current block to the pool and takes a free one.
  bool submitCurrent()
  {
    {
      std::lock_guard<std::mutex> lock(this->state->mtx);
      if(this->state->error != 0) { return false; }
      this->state->in_flight++;
    }
    DirectWritePool::instance().submit({ this->state, std::move(this->current), this->block_offset });
    this->block_offset += DIRECT_BLOCK_BYTES; this->used = 0;
    std::unique_lock<std::mutex> lock(this->state->mtx);
    this->state->returned.wait(lock, [this]()
      { return !this->state->free_blocks.empty() || this->state->error != 0; });
    if(this->state->error != 0) { return false; }
    this->current = std::move(this->state->free_blocks.back());
    this->state->free_blocks.pop_back();
    return true;
  }

  std::shared_ptr<DirectWriteState> state;
  std::unique_ptr<AlignedBlock> current;
  std::size_t used;
  std::uint64_t block_offset;
  bool open_;
};

} // anonymous namespace

bool
OutputStream::directWrites()
{
  static const bool enabled = (environmentSize("GCSA_IO_DIRECT_WRITES", 0) != 0);
  return enabled;
}

OutputStream::OutputStream() :
  std::ostream(nullptr), buffer(), direct_io(false)
{
  this->setstate(std::ios_base::badbit);
}

OutputStream::OutputStream(const std::string& filename, std::ios_base::openmode mode) :
  OutputStream()
{
  this->open(filename, mode);
}

OutputStream::~OutputStream()
{
  if(this->is_open()) { this->close(); }
}

void
OutputStream::open(const std::string& filename, std::ios_base::openmode mode)
{
  if(this->is_open()) { this->setstate(std::ios_base::failbit); return; }
  const bool truncating = !(mode & (std::ios_base::app | std::ios_base::in | std::ios_base::ate));
  if(directWrites() && truncating)
  {
    std::unique_ptr<DirectOutputBuffer> direct(new DirectOutputBuffer());
    int error = direct->open(filename);
    if(error == 0)
    {
      this->buffer = std::move(direct); this->direct_io = true;
      this->rdbuf(this->buffer.get());
      return;
    }
    if(error != EINVAL) { this->setstate(std::ios_base::failbit); return; }
    static std::once_flag warned;
    std::call_once(warned, [&filename]()
    {
      std::cerr << "OutputStream: O_DIRECT is unsupported for " << filename
                << "; writing through the page cache" << std::endl;
    });
  }
  std::unique_ptr<std::filebuf> file(new std::filebuf());
  if(this->file_buffer_size >= 0) { file->pubsetbuf(this->file_buffer, this->file_buffer_size); }
  if(file->open(filename, mode | std::ios_base::out) == nullptr)
  {
    this->setstate(std::ios_base::failbit); return;
  }
  this->buffer = std::move(file); this->direct_io = false;
  this->rdbuf(this->buffer.get());
}

bool
OutputStream::is_open() const
{
  if(!(this->buffer)) { return false; }
  return (this->direct_io ?
    static_cast<const DirectOutputBuffer*>(this->buffer.get())->is_open() :
    static_cast<const std::filebuf*>(this->buffer.get())->is_open());
}

void
OutputStream::close()
{
  bool ok = false;
  if(this->buffer)
  {
    ok = (this->direct_io ?
      static_cast<DirectOutputBuffer*>(this->buffer.get())->close() :
      static_cast<std::filebuf*>(this->buffer.get())->close() != nullptr);
  }
  if(!ok) { this->setstate(std::ios_base::failbit); }
}

//------------------------------------------------------------------------------

} // namespace gcsa
