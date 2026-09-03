/*
  Copyright (c) 2026 GCSA2 contributors

  Prefix-compressed transient path-label sort runs.
*/

#ifndef GCSA_PATH_SORT_RUN_H
#define GCSA_PATH_SORT_RUN_H

#include <gcsa/path_graph.h>

#include <cstdint>
#include <string>
#include <sys/types.h>
#include <vector>

namespace gcsa
{

// Self-contained in-memory record used by external path-label sorting. The
// rank pointer in PathNode is normalized on disk, because labels are embedded.
struct PathSortRunRecord
{
  PathNode node;
  PathNode::rank_type labels[PathLabel::LABEL_LENGTH + 1];
};

/*
  Records are stored in a versioned little-endian stream. Each record contains
  the PathNode fields, a one-byte common-prefix length against the previous
  record, and only the remaining ranks. Both classes own fixed-size byte
  buffers and periodically evict completed file-cache prefixes.
*/
class PathSortRunWriter
{
public:
  PathSortRunWriter(const std::string& filename, size_type expected_records,
    size_type buffer_bytes);
  ~PathSortRunWriter();

  void write(const PathSortRunRecord& record);
  void finish();

  size_type records() const { return this->records_; }
  size_type bytes() const { return this->total_bytes_; }
  size_type payloadBytes() const { return this->payload_bytes_; }
  size_type uncompressedBytes() const;

private:
  std::string filename_;
  int file_;
  std::vector<std::uint8_t> buffer_;
  size_type used_, expected_records_, records_, payload_bytes_, total_bytes_;
  std::uint64_t checksum_;
  off_t cache_released_;
  PathNode::rank_type previous_[PathLabel::LABEL_LENGTH + 1];
  size_type previous_ranks_;
  bool finished_;

  void append(const void* data, size_type bytes, bool payload);
  void flush(bool durable);
  void writeHeader();

  PathSortRunWriter(const PathSortRunWriter&) = delete;
  PathSortRunWriter& operator=(const PathSortRunWriter&) = delete;
};

class PathSortRunReader
{
public:
  PathSortRunReader(const std::string& filename, size_type buffer_bytes);
  ~PathSortRunReader();

  const PathSortRunRecord& current() const;
  void advance();
  bool atEnd() const { return this->at_end_; }

  size_type records() const { return this->total_records_; }
  size_type bytes() const { return this->total_bytes_; }
  size_type bufferBytes() const { return this->buffer_.size(); }

private:
  std::string filename_;
  int file_;
  std::vector<std::uint8_t> buffer_;
  size_type begin_, end_, file_offset_, total_records_, records_read_;
  size_type payload_bytes_, payload_read_, total_bytes_;
  std::uint64_t expected_checksum_, checksum_;
  off_t cache_released_;
  PathSortRunRecord current_;
  PathNode::rank_type previous_[PathLabel::LABEL_LENGTH + 1];
  size_type previous_ranks_;
  bool at_end_, footer_checked_;

  void readExact(void* data, size_type bytes, bool payload);
  void refill();
  void readHeader();
  void readRecord();
  void readFooter();

  PathSortRunReader(const PathSortRunReader&) = delete;
  PathSortRunReader& operator=(const PathSortRunReader&) = delete;
  PathSortRunReader(PathSortRunReader&&) = delete;
  PathSortRunReader& operator=(PathSortRunReader&&) = delete;
};

size_type pathSortRunMinimumBuffer();

} // namespace gcsa

#endif // GCSA_PATH_SORT_RUN_H
