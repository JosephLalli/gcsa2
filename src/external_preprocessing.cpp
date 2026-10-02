/*
  Copyright (c) 2026 GCSA2 contributors

  Disk-first input preprocessing for GCSA construction.
*/

#include <gcsa/external_preprocessing.h>
#include <gcsa/internal.h>

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <unistd.h>
#include <vector>

namespace gcsa
{

namespace
{

const std::string PREPROCESS_TASK = "preprocess";
const std::string PREPROCESS_PHASE = "key-start-streams";

ArtifactIdentity
keyArtifact()
{
  return ArtifactIdentity(PREPROCESS_TASK, PREPROCESS_PHASE,
    "unique-keys", "key-stream-v1");
}

ArtifactIdentity
startArtifact()
{
  return ArtifactIdentity(PREPROCESS_TASK, PREPROCESS_PHASE,
    "unique-start-nodes", "node-stream-v1");
}

size_type
recordCount(const std::string& name, size_type record_bytes)
{
  std::ifstream input;
  input.rdbuf()->pubsetbuf(nullptr, 0);
  input.open(name.c_str(), std::ios_base::binary);
  if(!input) { throw std::runtime_error("external preprocessing: cannot open " + name); }
  input.seekg(0, std::ios_base::end);
  std::streamoff bytes = input.tellg();
  if(bytes < 0 || static_cast<std::uint64_t>(bytes) % record_bytes != 0)
  {
    throw std::runtime_error("external preprocessing: truncated record stream " + name);
  }
  return static_cast<size_type>(static_cast<std::uint64_t>(bytes) / record_bytes);
}

class RawWriter
{
public:
  RawWriter(const std::string& name, size_type buffer_bytes, bool track_digest = false) :
    name_(name), output_(), buffer_(std::max(static_cast<size_type>(1), buffer_bytes)), used_(0),
    track_digest_(track_digest), digest_(1469598103934665603ULL)
  {
    output_.setFileBuffer(nullptr, 0);
    output_.open(name.c_str(), std::ios_base::binary | std::ios_base::trunc);
    if(!output_) { throw std::runtime_error("external preprocessing: cannot create " + name); }
  }

  ~RawWriter()
  {
    try { this->close(); } catch(...) { }
  }

  void write(const void* record, size_type bytes)
  {
    if(bytes > buffer_.size())
    {
      this->flush();
      output_.write(reinterpret_cast<const char*>(record), bytes);
      if(!output_) { throw std::runtime_error("external preprocessing: cannot write " + name_); }
      if(track_digest_) { digest_ = BuildWorkspace::checksum(record, bytes, digest_); }
      DiskIO::write_volume += bytes;
      return;
    }
    if(used_ + bytes > buffer_.size()) { this->flush(); }
    std::memcpy(buffer_.data() + used_, record, bytes); used_ += bytes;
  }

  template<class Record>
  void write(const Record& record)
  {
    this->write(&record, sizeof(Record));
  }

  void flush()
  {
    if(used_ == 0) { return; }
    output_.write(reinterpret_cast<const char*>(buffer_.data()), used_);
    if(!output_) { throw std::runtime_error("external preprocessing: cannot write " + name_); }
    if(track_digest_) { digest_ = BuildWorkspace::checksum(buffer_.data(), used_, digest_); }
    DiskIO::write_volume += used_;
    used_ = 0;
  }

  void close()
  {
    if(!output_.is_open()) { return; }
    this->flush(); output_.close();
    if(!output_) { throw std::runtime_error("external preprocessing: cannot close " + name_); }
  }

  std::uint64_t storedChecksum() const
  {
    if(output_.is_open() || !track_digest_)
    { throw std::logic_error("raw checksum requires a closed tracked writer"); }
    return digest_;
  }

private:
  std::string name_;
  OutputStream output_;
  std::vector<std::uint8_t> buffer_;
  size_type used_;
  bool track_digest_;
  std::uint64_t digest_;
};

// A read-only descriptor that RawReader's ReadAhead can borrow. It is a member
// declared before the ReadAhead, so the read-ahead is destroyed first and has
// waited for its pieces in flight before the descriptor closes.
struct InputDescriptor
{
  int value;

  explicit InputDescriptor(const std::string& name) :
    value(::open(name.c_str(), O_RDONLY | O_CLOEXEC))
  {
    if(this->value < 0) { throw std::runtime_error("external preprocessing: cannot open " + name); }
#if defined(POSIX_FADV_SEQUENTIAL)
    // Advice is non-semantic; unsupported filesystems may safely ignore it.
    static_cast<void>(::posix_fadvise(this->value, 0, 0, POSIX_FADV_SEQUENTIAL));
#endif
  }

  ~InputDescriptor() { if(this->value >= 0) { ::close(this->value); } }

  InputDescriptor(const InputDescriptor&) = delete;
  InputDescriptor& operator=(const InputDescriptor&) = delete;
};

/*
  Reads the sorted KMer stream and the unique key and start streams for the
  mapping pass and the start-node and last-character scans. One synchronous
  read per buffer with no read-ahead leaves these passes at whatever the
  kernel's own read-ahead delivers on a busy shared device, so it reads in
  slices and advances a ReadAhead before each, as the sorter's readSliced
  does. The ReadAhead borrows this reader's descriptor and adds none to the
  descriptor budget.
*/
class RawReader
{
public:
  RawReader(const std::string& name, size_type record_bytes,
    size_type expected_records, size_type buffer_bytes) :
    name_(name), record_bytes_(checkedWidth(record_bytes)), input_(name),
    read_ahead_(input_.value), buffer_(), records_(0), offset_(0),
    remaining_(expected_records), position_(0)
  {
    size_type records_per_buffer = std::max(static_cast<size_type>(1),
      buffer_bytes / record_bytes_);
    buffer_.resize(records_per_buffer * record_bytes_);
  }

  bool next(void* record)
  {
    if(offset_ == records_)
    {
      if(remaining_ == 0) { return false; }
      size_type capacity = buffer_.size() / record_bytes_;
      size_type count = std::min(capacity, remaining_);
      if(this->readSliced(buffer_.data(), count * record_bytes_) != count * record_bytes_)
      {
        throw std::runtime_error("external preprocessing: truncated stream " + name_);
      }
      DiskIO::read_volume += count * record_bytes_;
      records_ = count; offset_ = 0; remaining_ -= count;
    }
    std::memcpy(record, buffer_.data() + offset_ * record_bytes_, record_bytes_);
    offset_++; return true;
  }

  void finish()
  {
    if(remaining_ != 0 || offset_ != records_)
    {
      throw std::runtime_error("external preprocessing: stream was not fully consumed " + name_);
    }
    std::vector<std::uint8_t> extra(record_bytes_);
    if(this->readSliced(extra.data(), record_bytes_) != 0)
    {
      throw std::runtime_error("external preprocessing: trailing records in " + name_);
    }
  }

private:
  constexpr static size_type SLICE_BYTES = 16 * MEGABYTE;

  static size_type checkedWidth(size_type record_bytes)
  {
    if(record_bytes == 0) { throw std::invalid_argument("external preprocessing: zero record width"); }
    return record_bytes;
  }

  // Reads up to `bytes`, stopping early only at the end of the file.
  size_type readSliced(std::uint8_t* data, size_type bytes)
  {
    size_type done = 0;
    while(done < bytes)
    {
      read_ahead_.advance(position_);
      ssize_t got = ::read(input_.value, data + done, std::min(SLICE_BYTES, bytes - done));
      if(got < 0 && errno == EINTR) { continue; }
      if(got < 0) { throw std::runtime_error("external preprocessing: cannot read " + name_); }
      if(got == 0) { break; }
      done += static_cast<size_type>(got); position_ += static_cast<std::uint64_t>(got);
    }
    return done;
  }

  std::string name_;
  size_type record_bytes_;
  InputDescriptor input_;
  ReadAhead read_ahead_;
  std::vector<std::uint8_t> buffer_;
  size_type records_, offset_, remaining_;
  std::uint64_t position_;
};

struct InitialKMerRecord
{
  KMer kmer;
  std::uint64_t ordinal;
};

static_assert(std::is_trivially_copyable<InitialKMerRecord>::value,
  "external initial KMer records must be fixed-width binary values");

int
compareKeyLabel(const void* left, const void* right)
{
  key_type a = *static_cast<const key_type*>(left);
  key_type b = *static_cast<const key_type*>(right);
  size_type a_label = Key::label(a), b_label = Key::label(b);
  return (a_label < b_label ? -1 : (a_label > b_label ? 1 : 0));
}

int
compareNode(const void* left, const void* right)
{
  node_type a = *static_cast<const node_type*>(left);
  node_type b = *static_cast<const node_type*>(right);
  return (a < b ? -1 : (a > b ? 1 : 0));
}

int
compareInitialKMer(const void* left, const void* right)
{
  const InitialKMerRecord& a = *static_cast<const InitialKMerRecord*>(left);
  const InitialKMerRecord& b = *static_cast<const InitialKMerRecord*>(right);
  size_type a_label = Key::label(a.kmer.key), b_label = Key::label(b.kmer.key);
  if(a_label != b_label) { return (a_label < b_label ? -1 : 1); }
  // Ordinal is assigned during one logical-file scan, making every normal
  // record distinct and preserving source order within a label deterministically.
  if(a.ordinal != b.ordinal) { return (a.ordinal < b.ordinal ? -1 : 1); }
  if(a.kmer.key != b.kmer.key) { return (a.kmer.key < b.kmer.key ? -1 : 1); }
  if(a.kmer.from != b.kmer.from) { return (a.kmer.from < b.kmer.from ? -1 : 1); }
  return (a.kmer.to < b.kmer.to ? -1 : (a.kmer.to > b.kmer.to ? 1 : 0));
}

node_type
mapStart(node_type node, const NodeMapping& mapping)
{
  if(mapping.empty()) { return node; }
  return Node::encode(mapping(Node::id(node)), Node::offset(node), Node::rc(node));
}

BuildWorkspace::ArtifactRef
checkpointStream(BuildWorkspace& workspace, const ArtifactIdentity& identity,
  physical_shard_id_t shard, const std::string& name, size_type records,
  size_type buffer_bytes)
{
  std::ifstream input;
  input.rdbuf()->pubsetbuf(nullptr, 0);
  input.open(name.c_str(), std::ios_base::binary);
  if(!input) { throw std::runtime_error("external preprocessing: cannot checkpoint " + name); }
  BuildWorkspace::ArtifactWriter writer = workspace.open_artifact(identity,
    logical_file_id_t(0), shard, "value", "all");
  std::vector<std::uint8_t> buffer(std::max(static_cast<size_type>(1), buffer_bytes));
  while(input)
  {
    input.read(reinterpret_cast<char*>(buffer.data()), buffer.size());
    std::streamsize bytes = input.gcount();
    if(bytes > 0)
    {
      DiskIO::read_volume += static_cast<size_type>(bytes);
      writer.write(buffer.data(), static_cast<size_type>(bytes));
    }
  }
  if(!input.eof()) { throw std::runtime_error("external preprocessing: failed while checkpointing " + name); }
  return writer.finish(records);
}

physical_shard_id_t
physicalId(logical_file_id_t logical, std::uint64_t local_shard)
{
  if(local_shard > std::numeric_limits<std::uint32_t>::max())
  {
    throw std::runtime_error("external preprocessing: too many physical shards for one logical input");
  }
  return physical_shard_id_t((static_cast<std::uint64_t>(logical.value) << 32) | local_shard);
}

} // namespace

ExternalInputPreprocessor::ExternalInputPreprocessor(const InputGraph& graph,
  const ConstructionParameters& parameters, BuildWorkspace* workspace) :
  graph_(graph), parameters_(parameters), workspace_(workspace),
  key_name_(), start_name_(), key_count_(0), start_count_(0), prepared_(false), stats_()
{
}

ExternalInputPreprocessor::~ExternalInputPreprocessor()
{
  TempFile::remove(this->key_name_);
  TempFile::remove(this->start_name_);
}

size_type
ExternalInputPreprocessor::sortBudget() const
{
  size_type minimum = std::max(
    ExternalFixedRecordSorter::minimumBudget(sizeof(key_type)),
    std::max(ExternalFixedRecordSorter::minimumBudget(sizeof(node_type)),
      ExternalFixedRecordSorter::minimumBudget(sizeof(InitialKMerRecord))));
  size_type budget = std::min(this->parameters_.getMemoryLimitBytes(),
    this->parameters_.getSortRunSize());
  if(budget < minimum)
  {
    throw std::runtime_error("external preprocessing: memory/sort-run budget is below the fixed-record sort minimum");
  }
  return budget;
}

size_type
ExternalInputPreprocessor::ioBufferBytes() const
{
  size_type budget = this->sortBudget();
  size_type result = std::min(this->parameters_.getIOBufferSize(), budget / 8);
  return std::max(static_cast<size_type>(sizeof(InitialKMerRecord)), result);
}

void
ExternalInputPreprocessor::prepare()
{
  if(this->prepared_) { return; }
  this->key_name_ = TempFile::getName("gcsa_unique_keys");
  this->start_name_ = TempFile::getName("gcsa_unique_starts");
  const size_type io_bytes = this->ioBufferBytes();
  if(this->workspace_ != nullptr &&
     this->workspace_->task_completed(PREPROCESS_TASK, PREPROCESS_PHASE))
  {
    this->workspace_->restore_artifact(keyArtifact(), logical_file_id_t(0),
      physical_shard_id_t(0), this->key_name_, io_bytes);
    this->workspace_->restore_artifact(startArtifact(), logical_file_id_t(0),
      physical_shard_id_t(1), this->start_name_, io_bytes);
    this->key_count_ = recordCount(this->key_name_, sizeof(key_type));
    this->start_count_ = recordCount(this->start_name_, sizeof(node_type));
    this->stats_.unique_keys = this->key_count_;
    this->stats_.unique_start_nodes = this->start_count_;
    this->prepared_ = true; return;
  }

  std::string key_source = TempFile::getName("gcsa_key_source");
  std::string start_source = TempFile::getName("gcsa_start_source");
  try
  {
    {
      // These source writers must be destroyed before either sort begins:
      // retaining their explicit buffers alongside a full sort reservation
      // would violate the external working-set contract.
      RawWriter keys(key_source, io_bytes), starts(start_source, io_bytes);
      const size_type block_bytes = std::max(static_cast<size_type>(sizeof(KMer)),
        std::min(this->sortBudget() / 8, this->parameters_.getMemoryLimitBytes() / 8));
      for(size_type file = 0; file < this->graph_.files(); file++)
      {
        this->graph_.scanKMerBlocks(file, block_bytes,
          [&keys, &starts, this](size_type, const std::vector<KMer>& kmers)
          {
            for(const KMer& kmer : kmers)
            {
              keys.write(kmer.key);
              node_type from = mapStart(kmer.from, this->graph_.mapping);
              starts.write(from);
            }
          });
      }
      keys.close(); starts.close();
    }

    // compareKeyLabel orders on the label alone, so keys that differ in their
    // predecessor/successor bits compare equal, and the run sorter would keep
    // their input order with an offset permutation and one indirect call per
    // comparison: 0:13:01 on one core for the joint chr2+chr18 keys, against
    // 0:00:34 for the start sort over the same 809 million eight-byte records.
    // That order cannot be observed. The label is the key's top bits
    // (Key::label() is key >> 16), so ordering by the whole key is also an
    // order by label, and the merge and the grouping below still use
    // compareKeyLabel, so every label group stays contiguous. Within a group
    // the reducer keeps the label and ORs the low 16 bits (Key::merge()),
    // which commutes, so each group's output is the same in any order.
    // Declaring the whole-key order lets run formation sort the keys as
    // scalars in place, the order the legacy readKeys() uses as well.
    key_type merged_key = 0;
    ExternalFixedRecordSorter::sortAndReduce(key_source, this->key_name_,
      sizeof(key_type), this->sortBudget(), this->parameters_.getMergeFanIn(),
      compareKeyLabel,
      [&merged_key](const void* value, bool first, bool last, std::ostream& output)
      {
        key_type key = *static_cast<const key_type*>(value);
        if(first) { merged_key = key; }
        else { merged_key = Key::merge(merged_key, key); }
        if(last) { output.write(reinterpret_cast<const char*>(&merged_key), sizeof(merged_key)); }
      }, &this->stats_.key_sort, true,
      ExternalFixedRecordSorter::RecordOrder::ASCENDING_U64);
    ExternalFixedRecordSorter::sortAndReduce(start_source, this->start_name_,
      sizeof(node_type), this->sortBudget(), this->parameters_.getMergeFanIn(),
      compareNode,
      [](const void* value, bool first, bool, std::ostream& output)
      {
        if(first) { output.write(reinterpret_cast<const char*>(value), sizeof(node_type)); }
      }, &this->stats_.start_sort, true,
      ExternalFixedRecordSorter::RecordOrder::ASCENDING_U64);
    TempFile::remove(key_source);
    TempFile::remove(start_source);
  }
  catch(...)
  {
    TempFile::remove(key_source); TempFile::remove(start_source);
    throw;
  }

  this->key_count_ = recordCount(this->key_name_, sizeof(key_type));
  this->start_count_ = recordCount(this->start_name_, sizeof(node_type));
  if(this->key_count_ == 0)
  {
    throw std::runtime_error("external preprocessing: input produced no valid keys");
  }
  if(this->workspace_ != nullptr)
  {
    std::vector<BuildWorkspace::ArtifactRef> artifacts;
    artifacts.push_back(checkpointStream(*this->workspace_, keyArtifact(),
      physical_shard_id_t(0), this->key_name_, this->key_count_, io_bytes));
    artifacts.push_back(checkpointStream(*this->workspace_, startArtifact(),
      physical_shard_id_t(1), this->start_name_, this->start_count_, io_bytes));
    this->workspace_->commit_task(PREPROCESS_TASK, PREPROCESS_PHASE, artifacts);
  }
  this->stats_.unique_keys = this->key_count_;
  this->stats_.unique_start_nodes = this->start_count_;
  this->prepared_ = true;
}

size_type
ExternalInputPreprocessor::keyCount() const
{
  const_cast<ExternalInputPreprocessor*>(this)->prepare();
  return this->key_count_;
}

size_type
ExternalInputPreprocessor::startNodeCount() const
{
  const_cast<ExternalInputPreprocessor*>(this)->prepare();
  return this->start_count_;
}

void
ExternalInputPreprocessor::buildKeySupport(DeBruijnGraph& mapper, LCP& lcp,
  sdsl::int_vector<0>& last_char)
{
  this->buildMapper(mapper);
  this->buildLCP(lcp);
  this->buildLastCharacters(last_char);
}

void
ExternalInputPreprocessor::buildMapper(DeBruijnGraph& mapper)
{
  this->prepare();
  const size_type io_bytes = this->ioBufferBytes();
  mapper = DeBruijnGraph(this->key_name_, this->key_count_, this->graph_.k(),
    this->graph_.alpha, io_bytes);
}

void
ExternalInputPreprocessor::buildLCP(LCP& lcp)
{
  this->prepare();
  const size_type io_bytes = this->ioBufferBytes();
  lcp = LCP(this->key_name_, this->key_count_, this->graph_.k(), io_bytes);
}

void
ExternalInputPreprocessor::buildLastCharacters(sdsl::int_vector<0>& last_char)
{
  this->prepare();
  const size_type io_bytes = this->ioBufferBytes();
  last_char = sdsl::int_vector<0>(this->key_count_, 0, Key::GCSA_CHAR_WIDTH);
  RawReader keys(this->key_name_, sizeof(key_type), this->key_count_, io_bytes);
  key_type key = 0;
  for(size_type i = 0; i < this->key_count_; i++)
  {
    if(!keys.next(&key)) { throw std::runtime_error("external preprocessing: missing key"); }
    last_char[i] = Key::last(key);
  }
  keys.finish();
}

void
ExternalInputPreprocessor::buildStartNodes(sdsl::sd_vector<>& from_nodes)
{
  this->prepare();
  if(this->start_count_ == 0)
  {
    from_nodes = sdsl::sd_vector<>(); return;
  }
  const size_type io_bytes = this->ioBufferBytes();
  RawReader first_pass(this->start_name_, sizeof(node_type), this->start_count_, io_bytes);
  node_type previous = 0, maximum = 0;
  for(size_type i = 0; i < this->start_count_; i++)
  {
    node_type node = 0;
    if(!first_pass.next(&node)) { throw std::runtime_error("external preprocessing: missing start node"); }
    if(i > 0 && node <= previous)
    {
      throw std::runtime_error("external preprocessing: unique start stream is not strictly sorted");
    }
    previous = maximum = node;
  }
  first_pass.finish();
  if(maximum == std::numeric_limits<node_type>::max())
  {
    throw std::runtime_error("external preprocessing: start-node universe overflows sd_vector");
  }
  sdsl::sd_vector_builder builder(maximum + 1, this->start_count_);
  RawReader second_pass(this->start_name_, sizeof(node_type), this->start_count_, io_bytes);
  for(size_type i = 0; i < this->start_count_; i++)
  {
    node_type node = 0;
    if(!second_pass.next(&node)) { throw std::runtime_error("external preprocessing: missing start node"); }
    builder.set_unsafe(node);
  }
  second_pass.finish();
  from_nodes = sdsl::sd_vector<>(builder);
}

void
ExternalInputPreprocessor::buildInitialPathGraph(PathGraph& result)
{
  this->prepare();
  result.clear();
  PathGraph initial(0, this->graph_.k(), 0);
  if(this->graph_.files() > std::numeric_limits<std::uint32_t>::max())
  {
    throw std::runtime_error("external preprocessing: too many logical input files");
  }
  const size_type sort_budget = this->sortBudget();
  const size_type io_bytes = this->ioBufferBytes();
  const size_type path_bytes = sizeof(PathNode) + 2 * sizeof(PathNode::rank_type);
  const size_type target_shard_bytes = std::max(path_bytes, sort_budget / 2);
  const size_type paths_per_shard = std::max(static_cast<size_type>(1),
    std::min(this->parameters_.getCheckpointRecords(), target_shard_bytes / path_bytes));

  for(size_type file = 0; file < this->graph_.files(); file++)
  {
    std::string raw_name = TempFile::getName("gcsa_initial_kmers");
    std::string sorted_name = TempFile::getName("gcsa_initial_kmers_sorted");
    try
    {
      {
        // As above, the source writer's bounded buffer is released before the
        // sorter takes its complete sort_budget reservation.
        RawWriter raw(raw_name, io_bytes);
        std::uint64_t ordinal = 0;
        this->graph_.scanKMerBlocks(file,
          std::max(static_cast<size_type>(sizeof(KMer)), sort_budget / 8),
          [&raw, &ordinal](size_type, const std::vector<KMer>& kmers)
          {
            for(const KMer& kmer : kmers)
            {
              InitialKMerRecord record = {};
              record.kmer = kmer; record.ordinal = ordinal++;
              raw.write(record);
            }
          });
        raw.close();
      }
      ExternalFixedRecordSortStats kmer_stats;
      ExternalFixedRecordSorter::sort(raw_name, sorted_name,
        sizeof(InitialKMerRecord), sort_budget, this->parameters_.getMergeFanIn(),
        compareInitialKMer, &kmer_stats, true);
      this->stats_.logical_kmer_sorts++;

      // Mapping uses two bounded readers and two bounded path/rank writers.
      // ioBufferBytes() is at most sortBudget()/8, so their four explicit
      // buffers stay below half the sort budget. The subsequent succinct
      // vectors (mapper/LCP/start set) are final index support, not temporary
      // preprocessing working state, and are intentionally outside this
      // spill-stage accounting.
      RawReader kmers(sorted_name, sizeof(InitialKMerRecord), this->graph_.sizes[file], io_bytes);
      RawReader keys(this->key_name_, sizeof(key_type), this->key_count_, io_bytes);
      key_type key = 0;
      if(!keys.next(&key)) { throw std::runtime_error("external preprocessing: empty key stream"); }
      size_type key_rank = 0, path_count = 0, rank_count = 0;
      std::uint64_t local_shard = 0;
      std::string path_name, rank_name;
      std::unique_ptr<RawWriter> paths, ranks;
      const auto open_shard = [&]()
      {
        path_name = TempFile::getName(PathGraph::PREFIX);
        rank_name = TempFile::getName(PathGraph::PREFIX);
        paths.reset(new RawWriter(path_name, io_bytes, true));
        ranks.reset(new RawWriter(rank_name, io_bytes, true));
        path_count = rank_count = 0;
      };
      const auto close_shard = [&]()
      {
        if(!paths) { return; }
        paths->close(); ranks->close();
        const std::uint64_t path_digest = paths->storedChecksum(), rank_digest = ranks->storedChecksum();
        paths.reset(); ranks.reset();
        if(path_count == 0)
        {
          TempFile::remove(path_name); TempFile::remove(rank_name); return;
        }
        initial.path_names.push_back(path_name); initial.rank_names.push_back(rank_name);
        initial.path_checksums.emplace_back(); initial.rank_checksums.emplace_back();
        initial.path_checksums.back().record(path_name, path_digest);
        initial.rank_checksums.back().record(rank_name, rank_digest);
        initial.path_counts.push_back(path_count); initial.rank_counts.push_back(rank_count);
        initial.logical_file_ids.push_back(logical_file_id_t(static_cast<std::uint32_t>(file)));
        initial.physical_shard_ids.push_back(physicalId(logical_file_id_t(static_cast<std::uint32_t>(file)), local_shard));
        initial.path_count += path_count; initial.rank_count += rank_count;
        local_shard++; this->stats_.physical_shards++;
      };
      open_shard();
      for(size_type i = 0; i < this->graph_.sizes[file]; i++)
      {
        InitialKMerRecord record = {};
        if(!kmers.next(&record)) { throw std::runtime_error("external preprocessing: missing sorted KMer"); }
        const size_type label = Key::label(record.kmer.key);
        while(label > Key::label(key))
        {
          key_rank++;
          if(!keys.next(&key))
          {
            throw std::runtime_error("external preprocessing: KMer label absent from key stream");
          }
        }
        if(label != Key::label(key))
        {
          throw std::runtime_error("external preprocessing: KMer label absent from key stream");
        }
        PathNode node;
        node.from = record.kmer.from; node.to = record.kmer.to; node.fields = 0;
        if(record.kmer.sorted()) { node.makeSorted(); }
        node.setPredecessors(Key::predecessors(record.kmer.key));
        node.setOrder(1); node.setLCP(1); node.setPointer(rank_count);
        PathNode::rank_type labels[2] = {
          static_cast<PathNode::rank_type>(key_rank), static_cast<PathNode::rank_type>(0)
        };
        paths->write(node); ranks->write(labels, sizeof(labels));
        path_count++; rank_count += 2;
        if(path_count == paths_per_shard)
        {
          close_shard(); open_shard();
        }
      }
      // The KMer stream must be fully consumed. The unique key stream usually
      // has labels absent from this particular logical file, so it is valid to
      // stop after the final needed rank and let RawReader close normally.
      kmers.finish(); close_shard();
      TempFile::remove(raw_name);
      TempFile::remove(sorted_name);
    }
    catch(...)
    {
      TempFile::remove(raw_name); TempFile::remove(sorted_name);
      throw;
    }
  }
  result.swap(initial);
}

} // namespace gcsa
