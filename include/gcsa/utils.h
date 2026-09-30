#ifndef GCSA_UTILS_H
#define GCSA_UTILS_H

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <sdsl/bit_vector_il.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/sd_vector.hpp>

// TODO Later: Get rid of OpenMP.
#include <omp.h>

// Parallel sorting is only available with libstdc++ parallel mode.
#ifdef __GLIBCXX__
#include <parallel/algorithm>
#endif

namespace gcsa
{

/*
  utils.h: Common utility methods.
*/

//------------------------------------------------------------------------------

typedef std::uint64_t size_type;
typedef std::uint8_t  char_type;
typedef std::uint8_t  comp_type;
typedef std::uint8_t  byte_type;

constexpr size_type BYTE_BITS = 8;
constexpr size_type WORD_BITS = 64;

constexpr size_type KILOBYTE     = 1024;
constexpr size_type MILLION      = 1000000;
constexpr size_type MEGABYTE     = KILOBYTE * KILOBYTE;
constexpr size_type GIGABYTE     = KILOBYTE * MEGABYTE;

constexpr double KILOBYTE_DOUBLE = 1024.0;
constexpr double MILLION_DOUBLE  = 1000000.0;
constexpr double MEGABYTE_DOUBLE = KILOBYTE_DOUBLE * KILOBYTE_DOUBLE;
constexpr double GIGABYTE_DOUBLE = KILOBYTE_DOUBLE * MEGABYTE_DOUBLE;

// Exit code for exceeding parameter-defined size limit.
constexpr int EXIT_SIZE_LIMIT_EXCEEDED = 0x2BAD512E;

//------------------------------------------------------------------------------

/*
  range_type stores a closed range [first, second]. Empty ranges are indicated by
  first > second. The emptiness check uses +1 to handle the common special case
  [0, -1].
*/

typedef std::pair<size_type, size_type> range_type;

struct Range
{
  inline static size_type length(range_type range)
  {
    return range.second + 1 - range.first;
  }

  inline static bool empty(range_type range)
  {
    return (range.first + 1 > range.second + 1);
  }

  inline static bool empty(size_type sp, size_type ep)
  {
    return (sp + 1 > ep + 1);
  }

  inline static size_type bound(size_type value, range_type bounds)
  {
    return bound(value, bounds.first, bounds.second);
  }

  inline static size_type bound(size_type value, size_type low, size_type high)
  {
    return std::max(std::min(value, high), low);
  }

  inline static range_type empty_range()
  {
    return range_type(1, 0);
  }
};

template<class A, class B>
std::ostream& operator<<(std::ostream& stream, const std::pair<A, B>& data)
{
  return stream << "(" << data.first << ", " << data.second << ")";
}

//------------------------------------------------------------------------------

/*
  Global verbosity setting for index construction. Used in conditions of type
  if(Verbosity::level >= Verbosity::THRESHOLD). While the level can be set directly,
  Verbosity::set() does a few sanity checks.

  SILENT    no status information
  BASIC     basic progress information and statistics on the input and the final index
  EXTENDED  add intermediate statistics
  FULL      add technical information on processing individual files
*/

struct Verbosity
{
  static size_type level;

  static void set(size_type new_level);
  static std::string levelName();

  constexpr static size_type SILENT   = 0;
  constexpr static size_type BASIC    = 1;
  constexpr static size_type EXTENDED = 2;
  constexpr static size_type FULL     = 3;

  // The default is now SILENT, as it's more convenient for users.
  constexpr static size_type DEFAULT = SILENT;
};

//------------------------------------------------------------------------------

struct Version
{
  static std::string str(bool verbose = false);
  static void print(std::ostream& out, const std::string& tool_name, bool verbose = false, size_type new_lines = 2);

  constexpr static size_type MAJOR_VERSION = 1;
  constexpr static size_type MINOR_VERSION = 4;
  constexpr static size_type PATCH_VERSION = 0;

  constexpr static size_type GCSA_VERSION  = 3;
  constexpr static size_type LCP_VERSION   = 1;
};

//------------------------------------------------------------------------------

template<class IntegerType>
inline size_type
bit_length(IntegerType val)
{
  return sdsl::bits::hi(val) + 1;
}

//------------------------------------------------------------------------------

/*
  Thomas Wang's integer hash function. In many implementations, std::hash
  is identity function for integers, which leads to performance issues.
*/

inline size_type
wang_hash_64(size_type key)
{
  key = (~key) + (key << 21); // key = (key << 21) - key - 1;
  key = key ^ (key >> 24);
  key = (key + (key << 3)) + (key << 8); // key * 265
  key = key ^ (key >> 14);
  key = (key + (key << 2)) + (key << 4); // key * 21
  key = key ^ (key >> 28);
  key = key + (key << 31);
  return key;
}

//------------------------------------------------------------------------------

inline double
inMegabytes(size_type bytes)
{
  return bytes / MEGABYTE_DOUBLE;
}

inline double
inGigabytes(size_type bytes)
{
  return bytes / GIGABYTE_DOUBLE;
}

inline double
inBPC(size_type bytes, size_type size)
{
  return (8.0 * bytes) / size;
}

inline double
inMicroseconds(double seconds)
{
  return seconds * MILLION_DOUBLE;
}

const size_type DEFAULT_INDENT = 18;

void printHeader(const std::string& header, size_type indent = DEFAULT_INDENT);
void printTime(const std::string& header, size_type queries, double seconds, size_type indent = DEFAULT_INDENT);

//------------------------------------------------------------------------------

double readTimer();       // Seconds from an arbitrary time point.
size_type memoryUsage();  // Peak memory usage in bytes.

size_type readVolume();   // Only for GCSA construction.
size_type writeVolume();  // Only for GCSA construction.

// Parse a non-negative byte count with an optional binary suffix (K, M, G, T,
// P, optionally followed by B or iB). Throws std::invalid_argument or
// std::out_of_range rather than silently wrapping.
size_type parseBytes(const std::string& value);
std::string formatBytes(size_type bytes);

//------------------------------------------------------------------------------

/*
  Temporary file names have the pattern "prefix_hostname_pid_counter", where
  - prefix is given as an argument to getName();
  - hostname is the name of the host;
  - pid is the process id; and
  - counter is a running counter starting from 0.

  The generated names are stored until the file is deleted with remove(). All
  remaining temporary files are deleted when the program exits (normally or
  with std::exit()).

  getName(), remove(), forget() and setDirectory() hold one mutex, so threads
  may name and remove temporaries concurrently. temp_dir may be read while no
  thread is setting it.
*/

namespace TempFile
{
  extern const std::string DEFAULT_TEMP_DIR;
  extern std::string temp_dir;

  void setDirectory(const std::string& directory);
  std::string getName(const std::string& name_part);
  void remove(std::string& filename);  // Also clears the filename.
  // Forget about current temporary files so that they aren't deleted.
  void forget();
}

// Returns the total length of the rows, excluding line ends.
size_type readRows(const std::string& filename, std::vector<std::string>& rows, bool skip_empty_rows);

size_type fileSize(std::ifstream& file);
size_type fileSize(std::ofstream& file);

//------------------------------------------------------------------------------

/*
  parallelQuickSort() uses less working space than parallelMergeSort(). Calling omp_set_nested(1)
  improves the speed of parallelQuickSort().
*/

template<class Iterator, class Comparator>
void
parallelQuickSort(Iterator first, Iterator last, const Comparator& comp)
{
#ifdef __GLIBCXX__
  int nested = omp_get_nested();
  omp_set_nested(1);
  __gnu_parallel::sort(first, last, comp, __gnu_parallel::balanced_quicksort_tag());
  omp_set_nested(nested);
#else
  std::sort(first, last, comp);
#endif
}

template<class Iterator>
void
parallelQuickSort(Iterator first, Iterator last)
{
#ifdef __GLIBCXX__
  int nested = omp_get_nested();
  omp_set_nested(1);
  __gnu_parallel::sort(first, last, __gnu_parallel::balanced_quicksort_tag());
  omp_set_nested(nested);
#else
  std::sort(first, last);
#endif
}

template<class Iterator, class Comparator>
void
parallelMergeSort(Iterator first, Iterator last, const Comparator& comp)
{
#ifdef __GLIBCXX__
  __gnu_parallel::sort(first, last, comp, __gnu_parallel::multiway_mergesort_tag());
#else
  std::sort(first, last, comp);
#endif
}

template<class Iterator>
void
parallelMergeSort(Iterator first, Iterator last)
{
#ifdef __GLIBCXX__
  __gnu_parallel::sort(first, last, __gnu_parallel::multiway_mergesort_tag());
#else
  std::sort(first, last);
#endif
}

template<class Iterator, class Comparator>
void
sequentialSort(Iterator first, Iterator last, const Comparator& comp)
{
#ifdef __GLIBCXX__
  __gnu_parallel::sort(first, last, comp, __gnu_parallel::sequential_tag());
#else
  std::sort(first, last, comp);
#endif
}

template<class Iterator>
void
sequentialSort(Iterator first, Iterator last)
{
#ifdef __GLIBCXX__
  __gnu_parallel::sort(first, last, __gnu_parallel::sequential_tag());
#else
  std::sort(first, last);
#endif
}

template<class Element>
void
removeDuplicates(std::vector<Element>& vec, bool parallel)
{
  if(parallel) { parallelQuickSort(vec.begin(), vec.end()); }
  else         { sequentialSort(vec.begin(), vec.end()); }
  vec.resize(std::unique(vec.begin(), vec.end()) - vec.begin());
}

template<class Element, class Random>
void
deterministicShuffle(std::vector<Element>& vec, Random& rng, bool parallel)
{
  // We sort to get a deterministic shuffle of the elements instead of the offsets.
  if(parallel) { parallelQuickSort(vec.begin(), vec.end()); }
  else         { sequentialSort(vec.begin(), vec.end()); }
  for(size_type i = vec.size(); i > 0; i--)
  {
    std::swap(vec[i - 1], vec[rng() % i]);
  }
}

//------------------------------------------------------------------------------

/*
  Split the vector approximately evenly between the given number of threads. The comparator
  should return true when it is safe to split between the first argument and the second argument.
*/

template<class Element, class Comparator>
std::vector<range_type>
getBounds(const std::vector<Element>& vec, size_type threads, const Comparator& comp)
{
  std::vector<range_type> bounds(threads);
  for(size_type thread = 0, start = 0; thread < threads; thread++)
  {
    bounds[thread].first = start;
    if(start < vec.size())
    {
      start += std::max((size_type)1, (vec.size() - start) / (threads - thread));
      while(start < vec.size() && !comp(vec[start - 1], vec[start])) { start++; }
    }
    bounds[thread].second = start - 1;
  }
  return bounds;
}

// Get a good chunk size for dynamic scheduling.
size_type getChunkSize(size_type n, size_type min_size);

//------------------------------------------------------------------------------

/*
  GCSA uses a contiguous byte alphabet [0, sigma - 1] internally. Array C is based on the
  number of occurrences of each character in the BWT.
*/

template<class AlphabetType>
inline bool
hasChar(const AlphabetType& alpha, comp_type comp)
{
  return (alpha.C[comp + 1] > alpha.C[comp]);
}

template<class AlphabetType>
inline range_type
charRange(const AlphabetType& alpha, comp_type comp)
{
  return range_type(alpha.C[comp], alpha.C[comp + 1] - 1);
}

template<class AlphabetType>
comp_type
findChar(const AlphabetType& alpha, size_type bwt_pos)
{
  comp_type comp = 0;
  while(alpha.C[comp + 1] <= bwt_pos) { comp++; }
  return comp;
}

template<class BWTType, class AlphabetType>
inline size_type
LF(const BWTType& bwt, const AlphabetType& alpha, size_type i, comp_type comp)
{
  return alpha.C[comp] + bwt.rank(i, comp);
}

template<class BWTType, class AlphabetType>
inline range_type
LF(const BWTType& bwt, const AlphabetType& alpha, range_type range, comp_type comp)
{
  return range_type(LF(bwt, alpha, range.first, comp), LF(bwt, alpha, range.second + 1, comp) - 1);
}

//------------------------------------------------------------------------------

/*
  Progress for a phase whose total work is known in advance.

  A chromosome-scale doubling step can run for hours while printing nothing, so
  the only way to tell a working build from a stuck one has been to watch the
  workspace grow and guess. This reports a completed fraction, a rate and a
  remaining estimate, throttled so that the cost is one clock read per update
  call and one line of output per interval.

  advance() is safe to call in the innermost loop: it compares a counter
  against a threshold and returns, and only consults the clock when that
  threshold is crossed. It reports nothing below BASIC verbosity, and the
  totals it prints are the phase's own units (paths, records, partitions), not
  bytes, so a reader can check them against the counts the phase reports when
  it finishes.
*/
class ProgressReporter
{
public:
  // interval_seconds <= 0 disables reporting entirely.
  ProgressReporter(const std::string& phase_name, size_type total_units,
    const std::string& unit_name = "records", double interval_seconds = 60.0);

  // Advance by `units` and report if the interval has elapsed.
  void advance(size_type units = 1);

  // Report the final state regardless of the interval, once.
  void finish();

private:
  void report(bool final_report);

  std::string name, unit;
  size_type total, processed, next_check;
  double start_time, last_report, interval;
  bool done;
};

//------------------------------------------------------------------------------

/*
  FORK: read-ahead for sequential readers.

  On a shared device with high queue latency, a stream's throughput is bounded
  by the bytes it has in flight: one 1 MiB request per 30 ms is about 33 MB/s.
  A ReadAhead keeps a window ahead of its reader in the page cache by reading
  it in PIECE_BYTES pieces on a shared pool of threads, so the stream has many
  requests in flight while the reader's own reads hit the cache.

  GCSA_IO_READAHEAD_BYTES sets the window (0 or unset disables read-ahead) and
  GCSA_IO_READAHEAD_THREADS the pool size (default 16). A ReadAhead borrows its
  reader's descriptor where the reader has one, so it adds nothing to the
  --max-open-files budget; destroying it cancels queued pieces and waits for
  those in flight, after which the owner may close the descriptor. A reader
  without a descriptor gives a filename, and the ReadAhead opens its own.
  Read-ahead is advisory: its reads only populate the page cache, and the
  reader's own reads report any error.
*/
class ReadAhead
{
public:
  explicit ReadAhead(const std::string& filename);
  explicit ReadAhead(int descriptor);
  ~ReadAhead();

  // The reader's next read starts at byte `offset`.
  void advance(std::uint64_t offset);

  // Window size from GCSA_IO_READAHEAD_BYTES; 0 when read-ahead is disabled.
  static std::uint64_t windowBytes();

  constexpr static std::uint64_t PIECE_BYTES = std::uint64_t(1) << 20;

  struct Source;

  ReadAhead(const ReadAhead&) = delete;
  ReadAhead& operator=(const ReadAhead&) = delete;
  ReadAhead(ReadAhead&& another) noexcept;
  ReadAhead& operator=(ReadAhead&& another) noexcept;

private:
  // Cancels queued pieces and waits for those in flight.
  void release();

  std::shared_ptr<Source> source;
  std::uint64_t requested;
};

//------------------------------------------------------------------------------

/*
  FORK: sequential output that can bypass the page cache.

  Buffered writes on a shared host are paced by the kernel's global dirty-page
  limit, which other processes fill. With GCSA_IO_DIRECT_WRITES=1, an
  OutputStream opened for truncating output ("w" semantics) writes with
  O_DIRECT from aligned 4 MiB blocks, with up to GCSA_IO_DIRECT_WRITE_DEPTH
  blocks per stream (default 4) in flight on a shared pool of
  GCSA_IO_DIRECT_WRITE_THREADS threads (default 16), so its data never enters
  the page cache. flush() writes everything: it pads the final partial block,
  writes it, and truncates the file to its logical size; later writes rewrite
  that block. A filesystem without O_DIRECT falls back to buffered output with
  a warning. The blocks are outside GCSA2's memory budget.

  Otherwise an OutputStream behaves like std::ofstream, which it replaces.
*/
class OutputStream : public std::ostream
{
public:
  OutputStream();
  explicit OutputStream(const std::string& filename,
    std::ios_base::openmode mode = std::ios_base::out);
  ~OutputStream() override;

  void open(const std::string& filename,
    std::ios_base::openmode mode = std::ios_base::out);
  bool is_open() const;
  void close();

  // True when this stream writes with O_DIRECT.
  bool direct() const { return this->direct_io; }

  // What std::ofstream::rdbuf()->pubsetbuf() before open() did: applied to a
  // buffered stream at open(); a direct stream stages through its own blocks.
  void setFileBuffer(char* data, std::streamsize size) { this->file_buffer = data; this->file_buffer_size = size; }

  static bool directWrites();

  OutputStream(const OutputStream&) = delete;
  OutputStream& operator=(const OutputStream&) = delete;

private:
  std::unique_ptr<std::streambuf> buffer;
  bool direct_io;
  char* file_buffer = nullptr;
  std::streamsize file_buffer_size = -1;  // negative: the filebuf's default buffering
};

//------------------------------------------------------------------------------

} // namespace gcsa

#endif // GCSA_UTILS_H
