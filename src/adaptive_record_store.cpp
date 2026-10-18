#include <gcsa/adaptive_record_store.h>

#include <gcsa/internal.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <fcntl.h>
#include <limits>
#include <stdexcept>
#include <unistd.h>

namespace gcsa
{

namespace adaptive_record_store_detail
{

namespace
{

std::runtime_error
recordStoreError(const std::string& operation, const std::string& filename)
{
  return std::runtime_error("AdaptiveRecordStore: " + operation + ": " +
    filename + ": " + std::strerror(errno));
}

off_t
checkedOffset(size_type offset, size_type bytes, const std::string& filename)
{
  if(offset > static_cast<size_type>(std::numeric_limits<off_t>::max()) ||
     bytes > static_cast<size_type>(std::numeric_limits<off_t>::max()) - offset)
  {
    throw std::runtime_error("AdaptiveRecordStore: raw-file offset is too large: " + filename);
  }
  return static_cast<off_t>(offset);
}

} // namespace

int
openRawRecordFile(const std::string& filename)
{
  int flags = O_RDWR | O_CREAT | O_TRUNC;
#ifdef O_CLOEXEC
  flags |= O_CLOEXEC;
#endif
  int descriptor;
  do
  {
    descriptor = ::open(filename.c_str(), flags, 0600);
  }
  while(descriptor < 0 && errno == EINTR);
  if(descriptor < 0) { throw recordStoreError("cannot create raw store", filename); }
  return descriptor;
}

void
closeRawRecordFile(int descriptor) noexcept
{
  if(descriptor < 0) { return; }
  while(::close(descriptor) < 0 && errno == EINTR)
  {
  }
}

void
readRawRecords(int descriptor, void* target, size_type bytes, size_type offset,
  const std::string& filename)
{
  std::uint8_t* data = static_cast<std::uint8_t*>(target);
  size_type done = 0;
  const off_t base = checkedOffset(offset, bytes, filename);
  while(done < bytes)
  {
    const size_type request = std::min(bytes - done,
      static_cast<size_type>(std::numeric_limits<ssize_t>::max()));
    ssize_t read_bytes = ::pread(descriptor, data + done, request,
      base + static_cast<off_t>(done));
    if(read_bytes < 0 && errno == EINTR) { continue; }
    if(read_bytes <= 0) { throw recordStoreError("short raw-store read", filename); }
    done += static_cast<size_type>(read_bytes);
  }
  DiskIO::read_volume += bytes;
}

void
writeRawRecords(int descriptor, const void* source, size_type bytes, size_type offset,
  const std::string& filename)
{
  const std::uint8_t* data = static_cast<const std::uint8_t*>(source);
  size_type done = 0;
  const off_t base = checkedOffset(offset, bytes, filename);
  while(done < bytes)
  {
    const size_type request = std::min(bytes - done,
      static_cast<size_type>(std::numeric_limits<ssize_t>::max()));
    ssize_t written = ::pwrite(descriptor, data + done, request,
      base + static_cast<off_t>(done));
    if(written < 0 && errno == EINTR) { continue; }
    if(written <= 0) { throw recordStoreError("raw-store write failed", filename); }
    done += static_cast<size_type>(written);
  }
  DiskIO::write_volume += bytes;
}

void
copyRawRecords(int source, int destination, size_type bytes,
  const std::string& source_filename, const std::string& destination_filename)
{
  std::array<std::uint8_t, 64 * KILOBYTE> buffer;
  for(size_type offset = 0; offset < bytes; )
  {
    size_type chunk = std::min(bytes - offset, static_cast<size_type>(buffer.size()));
    readRawRecords(source, buffer.data(), chunk, offset, source_filename);
    writeRawRecords(destination, buffer.data(), chunk, offset, destination_filename);
    offset += chunk;
  }
}

} // namespace adaptive_record_store_detail

} // namespace gcsa
