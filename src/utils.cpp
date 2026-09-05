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

#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

#include <sys/resource.h>
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

  // By storing the filenames in a static object, we can delete the remaining
  // temporary files when std::exit() is called.
  struct Handler
  {
    std::set<std::string> filenames;
    ~Handler()
    {
      for(auto& filename : this->filenames)
      {
        std::remove(filename.c_str());
      }
    }
  } handler;

  void
  setDirectory(const std::string& directory)
  {
    if(directory.empty()) { temp_dir = DEFAULT_TEMP_DIR; }
    else if(directory[directory.length() - 1] != '/') { temp_dir = directory; }
    else { temp_dir = directory.substr(0, directory.length() - 1); }
  }

  std::string
  getName(const std::string& name_part)
  {
    char hostname[32];
    gethostname(hostname, 32); hostname[31] = 0;

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
      handler.filenames.erase(filename);
      filename.clear();
    }
  }

  void
  forget() {
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

} // namespace gcsa
