/*
  Copyright (c) 2026 GCSA2 contributors

  Disk-first input preprocessing for GCSA construction.
*/

#include <gcsa/external_preprocessing.h>
#include <gcsa/internal.h>

#include <algorithm>
#include <cerrno>
#include <condition_variable>
#include <cstdio>
#include <cstring>
#include <exception>
#include <fcntl.h>
#include <fstream>
#include <functional>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>
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

/*
  Passes a byte stream to a sink on one helper thread, in order: the producer
  fills one buffer while the helper hands the other to the sink, so the sink's
  writes and checksums overlap the producer's work instead of following it.
  Two buffers of buffer_bytes each; write() blocks only while the helper is
  still busy with the previous buffer. finish() hands over the last bytes,
  waits for the sink, and rethrows the first exception the sink threw; the
  destructor stops the helper without finishing, for unwinding.
*/
class BackgroundWriter
{
public:
  typedef std::function<void(const std::uint8_t* data, size_type bytes)> Sink;

  BackgroundWriter(size_type buffer_bytes, Sink sink) :
    sink_(std::move(sink)), filling_(std::max(static_cast<size_type>(1), buffer_bytes)),
    draining_(filling_.size()), used_(0), pending_(0), busy_(false), stopping_(false),
    error_(), helper_()
  {
    this->helper_ = std::thread(&BackgroundWriter::drain, this);
  }

  ~BackgroundWriter() { this->stop(); }

  void write(const void* data, size_type bytes)
  {
    if(this->used_ + bytes <= this->filling_.size())
    {
      std::memcpy(this->filling_.data() + this->used_, data, bytes);
      this->used_ += bytes; return;
    }
    this->writeSplit(static_cast<const std::uint8_t*>(data), bytes);
  }

  void finish()
  {
    if(this->used_ > 0) { this->handOver(); }
    this->stop();
    if(this->error_) { std::rethrow_exception(this->error_); }
  }

  BackgroundWriter(const BackgroundWriter&) = delete;
  BackgroundWriter& operator=(const BackgroundWriter&) = delete;

private:
  void writeSplit(const std::uint8_t* data, size_type bytes)
  {
    while(bytes > 0)
    {
      if(this->used_ == this->filling_.size()) { this->handOver(); }
      const size_type take = std::min(bytes, this->filling_.size() - this->used_);
      std::memcpy(this->filling_.data() + this->used_, data, take);
      this->used_ += take; data += take; bytes -= take;
    }
  }

  // Waits for the helper to finish the previous buffer, then gives it this one.
  void handOver()
  {
    std::unique_lock<std::mutex> lock(this->mutex_);
    this->idle_.wait(lock, [this]() { return !this->busy_; });
    if(this->error_) { std::rethrow_exception(this->error_); }
    this->filling_.swap(this->draining_);
    this->pending_ = this->used_; this->used_ = 0; this->busy_ = true;
    lock.unlock();
    this->ready_.notify_one();
  }

  // Waits for the buffer in flight and joins the helper.
  void stop()
  {
    if(!this->helper_.joinable()) { return; }
    {
      std::unique_lock<std::mutex> lock(this->mutex_);
      this->idle_.wait(lock, [this]() { return !this->busy_; });
      this->stopping_ = true;
    }
    this->ready_.notify_one();
    this->helper_.join();
  }

  void drain()
  {
    std::unique_lock<std::mutex> lock(this->mutex_);
    while(true)
    {
      this->ready_.wait(lock, [this]() { return this->busy_ || this->stopping_; });
      if(!this->busy_) { return; }
      lock.unlock();
      std::exception_ptr error;
      try { this->sink_(this->draining_.data(), this->pending_); }
      catch(...) { error = std::current_exception(); }
      lock.lock();
      if(error && !this->error_) { this->error_ = error; }
      this->busy_ = false;
      this->idle_.notify_all();
    }
  }

  Sink sink_;
  std::vector<std::uint8_t> filling_, draining_;
  size_type used_, pending_;
  bool busy_, stopping_;
  std::exception_ptr error_;
  std::mutex mutex_;
  std::condition_variable ready_, idle_;
  std::thread helper_;
};

/*
  Copies one reduced key or start stream into its workspace artifact while the
  sorter writes the stream's working file, through the sorter's output tap.
  The working file is still needed, because the mapper and LCP builders and
  the later scans read a raw stream from offset zero; what this removes is the
  separate checkpoint pass that read the file back and wrote and checksummed
  the artifact afterwards (0:00:36-0:02:40 for the joint chr2+chr18 streams).
  The artifact receives exactly the working file's bytes, in order, so its
  payload, length, checksum and header are those the copy pass produced.

  The artifact stays a .partial file until finish(), which preprocessing calls
  where it used to copy, after both reductions and the empty-key check, so
  the artifact renames, their crash points and the task commit keep their
  order. A crash or an exception before then leaves only a partial file, which
  the ArtifactWriter removes when unwinding and recovery removes after a crash.
*/
class ArtifactCopy
{
public:
  ArtifactCopy() : writer_(), copy_(), bytes_(0) { }

  void open(BuildWorkspace& workspace, const ArtifactIdentity& identity,
    physical_shard_id_t shard, size_type buffer_bytes)
  {
    this->writer_ = workspace.open_artifact(identity, logical_file_id_t(0), shard,
      "value", "all");
    this->copy_.reset(new BackgroundWriter(buffer_bytes,
      [this](const std::uint8_t* data, size_type bytes) { this->writer_.write(data, bytes); }));
  }

  // The sorter's tap; empty, so the sorter copies nothing, when not opened.
  ExternalFixedRecordSorter::OutputTap tap()
  {
    if(!this->copy_) { return ExternalFixedRecordSorter::OutputTap(); }
    return [this](const void* data, size_type bytes)
    {
      this->copy_->write(data, bytes); this->bytes_ += bytes;
    };
  }

  // Waits until the artifact has every byte and releases the helper and its
  // buffers before the next phase reserves memory.
  void drain()
  {
    if(!this->copy_) { return; }
    this->copy_->finish(); this->copy_.reset();
  }

  BuildWorkspace::ArtifactRef finish(size_type records, size_type bytes)
  {
    this->drain();
    if(bytes != this->bytes_)
    {
      throw std::runtime_error("external preprocessing: reduced stream and its artifact copy differ in length");
    }
    return this->writer_.finish(records);
  }

  ArtifactCopy(const ArtifactCopy&) = delete;
  ArtifactCopy& operator=(const ArtifactCopy&) = delete;

private:
  BuildWorkspace::ArtifactWriter writer_;
  // Declared after the writer it feeds, so it stops before the writer closes.
  std::unique_ptr<BackgroundWriter> copy_;
  size_type bytes_;
};

/*
  Writes a raw record stream, optionally with the FNV-1a digest of every byte
  that checkpoint adoption reuses instead of rereading the file. The writes
  and the digest run on a BackgroundWriter's helper thread, so the producer
  only copies records into a buffer. On the main thread the digest alone was
  the largest cost of the mapping pass: the initial path and rank shards are
  32 bytes per k-mer (24.1 GiB for the joint chr2+chr18 build) at about
  0.77 GB/s. The bytes, their order and the digest are unchanged; the stream
  is complete, and the digest final, once close() returns.
*/
class RawWriter
{
public:
  RawWriter(const std::string& name, size_type buffer_bytes, bool track_digest = false) :
    name_(name), output_(), track_digest_(track_digest), digest_(1469598103934665603ULL),
    background_()
  {
    output_.setFileBuffer(nullptr, 0);
    output_.open(name.c_str(), std::ios_base::binary | std::ios_base::trunc);
    if(!output_) { throw std::runtime_error("external preprocessing: cannot create " + name); }
    background_.reset(new BackgroundWriter(buffer_bytes,
      [this](const std::uint8_t* data, size_type bytes) { this->sink(data, bytes); }));
  }

  ~RawWriter()
  {
    try { this->close(); } catch(...) { }
  }

  void write(const void* record, size_type bytes)
  {
    background_->write(record, bytes);
  }

  template<class Record>
  void write(const Record& record)
  {
    this->write(&record, sizeof(Record));
  }

  void close()
  {
    if(!output_.is_open()) { return; }
    // The helper writes to output_, so it must stop before output_ closes,
    // also when it failed.
    std::unique_ptr<BackgroundWriter> background(std::move(background_));
    try { if(background) { background->finish(); } }
    catch(...) { background.reset(); output_.close(); throw; }
    background.reset();
    output_.close();
    if(!output_) { throw std::runtime_error("external preprocessing: cannot close " + name_); }
  }

  std::uint64_t storedChecksum() const
  {
    if(output_.is_open() || !track_digest_)
    { throw std::logic_error("raw checksum requires a closed tracked writer"); }
    return digest_;
  }

private:
  // Runs on the helper thread, which alone touches output_ and digest_ until
  // close() has joined it.
  void sink(const std::uint8_t* data, size_type bytes)
  {
    output_.write(reinterpret_cast<const char*>(data), bytes);
    if(!output_) { throw std::runtime_error("external preprocessing: cannot write " + name_); }
    if(track_digest_) { digest_ = BuildWorkspace::checksum(data, bytes, digest_); }
    DiskIO::write_volume += bytes;
  }

  std::string name_;
  OutputStream output_;
  bool track_digest_;
  std::uint64_t digest_;
  // Declared last, so a writer destroyed without close() stops its helper
  // before the stream it writes to.
  std::unique_ptr<BackgroundWriter> background_;
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
  const ConstructionParameters& parameters, BuildWorkspace* workspace,
  MemoryBudget* memory) :
  graph_(graph), parameters_(parameters), workspace_(workspace),
  key_name_(), start_name_(), key_count_(0), start_count_(0), prepared_(false), stats_(),
  memory_(memory), keys_(), starts_(), facts_memory_(), resident_(false),
  predecessor_edges_(0)
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
  if(this->memory_ != nullptr) { budget = std::min(budget, this->memory_->available()); }
  if(budget < minimum)
  {
    throw std::runtime_error("external preprocessing: memory/sort-run budget is below the fixed-record sort minimum");
  }
  return budget;
}

size_type
ExternalInputPreprocessor::ioBufferBytes() const
{
  size_type budget = std::min(this->parameters_.getMemoryLimitBytes(),
    this->parameters_.getSortRunSize());
  size_type result = std::min(this->parameters_.getIOBufferSize(), budget / 8);
  return std::max(static_cast<size_type>(sizeof(InitialKMerRecord)), result);
}

void
ExternalInputPreprocessor::prepare()
{
  if(this->prepared_) { return; }
  // A durable old preprocessing checkpoint remains readable. New small inputs
  // need no key/start files: count-sized vectors and the decode block share
  // the same lifetime budget as the path generations that follow.
  if((this->workspace_ == nullptr ||
      !this->workspace_->task_completed(PREPROCESS_TASK, PREPROCESS_PHASE)) &&
     this->tryPrepareResident())
  {
    this->prepared_ = true; return;
  }
  this->key_name_ = TempFile::getName("gcsa_unique_keys");
  this->start_name_ = TempFile::getName("gcsa_unique_starts");
  const size_type sort_budget = this->sortBudget();
  const size_type io_bytes = std::min(this->ioBufferBytes(), sort_budget / 8);
  MemoryBudget::Reservation phase_memory;
  if(this->memory_ != nullptr)
  {
    phase_memory = this->memory_->tryReserve(sort_budget, "external key/start preprocessing");
    if(!phase_memory) { throw std::runtime_error("cannot admit external preprocessing sort"); }
  }
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
  // With a workspace, each reduced stream is copied into its artifact while
  // it is written (see ArtifactCopy). During a reduction the sorter holds its
  // run reader and one output buffer, two of the fan-in + 1 buffers its merges
  // reserve, so a copy's two io-sized buffers (each at most sortBudget()/8)
  // fit in the working set the sort already reserved.
  ArtifactCopy key_copy, start_copy;
  try
  {
    {
      // These source writers must be destroyed before either sort begins:
      // retaining their explicit buffers alongside a full sort reservation
      // would violate the external working-set contract.
      RawWriter keys(key_source, io_bytes), starts(start_source, io_bytes);
      const size_type block_bytes = std::max(static_cast<size_type>(sizeof(KMer)),
        std::min(sort_budget / 8, this->parameters_.getMemoryLimitBytes() / 8));
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

    if(this->workspace_ != nullptr)
    {
      key_copy.open(*this->workspace_, keyArtifact(), physical_shard_id_t(0), io_bytes);
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
      sizeof(key_type), sort_budget, this->parameters_.getMergeFanIn(),
      compareKeyLabel,
      [&merged_key, this](const void* value, bool first, bool last, std::ostream& output)
      {
        key_type key = *static_cast<const key_type*>(value);
        if(first) { merged_key = key; }
        else { merged_key = Key::merge(merged_key, key); }
        if(last)
        {
          this->predecessor_edges_ += sdsl::bits::lt_cnt[Key::predecessors(merged_key)];
          output.write(reinterpret_cast<const char*>(&merged_key), sizeof(merged_key));
        }
      }, &this->stats_.key_sort, true,
      ExternalFixedRecordSorter::RecordOrder::ASCENDING_U64, key_copy.tap());
    key_copy.drain();
    if(this->workspace_ != nullptr)
    {
      start_copy.open(*this->workspace_, startArtifact(), physical_shard_id_t(1), io_bytes);
    }
    ExternalFixedRecordSorter::sortAndReduce(start_source, this->start_name_,
      sizeof(node_type), sort_budget, this->parameters_.getMergeFanIn(),
      compareNode,
      [](const void* value, bool first, bool, std::ostream& output)
      {
        if(first) { output.write(reinterpret_cast<const char*>(value), sizeof(node_type)); }
      }, &this->stats_.start_sort, true,
      ExternalFixedRecordSorter::RecordOrder::ASCENDING_U64, start_copy.tap());
    start_copy.drain();
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
    artifacts.push_back(key_copy.finish(this->key_count_,
      this->key_count_ * sizeof(key_type)));
    artifacts.push_back(start_copy.finish(this->start_count_,
      this->start_count_ * sizeof(node_type)));
    this->workspace_->commit_task(PREPROCESS_TASK, PREPROCESS_PHASE, artifacts);
  }
  this->stats_.unique_keys = this->key_count_;
  this->stats_.unique_start_nodes = this->start_count_;
  this->prepared_ = true;
}

bool
ExternalInputPreprocessor::tryPrepareResident()
{
  if(this->memory_ == nullptr || this->graph_.size() == 0) { return false; }
  const size_type record_bytes = sizeof(key_type) + sizeof(node_type);
  if(this->graph_.size() > std::numeric_limits<size_type>::max() / record_bytes)
  {
    return false;
  }
  MemoryBudget::Reservation facts = this->memory_->tryReserve(
    this->graph_.size() * record_bytes, "preprocessing keys and starts");
  if(!facts) { return false; }
  const size_type block_records = std::max(static_cast<size_type>(1),
    std::min(this->graph_.size(), this->parameters_.getIOBufferSize() / sizeof(KMer)));
  MemoryBudget::Reservation scratch = this->memory_->tryReserve(
    block_records * sizeof(KMer), "preprocessing decode block");
  if(!scratch) { return false; }

  std::vector<key_type> keys;
  std::vector<node_type> starts;
  keys.reserve(this->graph_.size()); starts.reserve(this->graph_.size());
  // reserve() is permitted to allocate more than requested. Do not leave that
  // excess outside the pool if a different standard library chooses to do so.
  MemoryBudget::Reservation excess = this->memory_->tryReserve(
    (keys.capacity() - this->graph_.size()) * sizeof(key_type) +
    (starts.capacity() - this->graph_.size()) * sizeof(node_type),
    "preprocessing vector excess");
  if(!excess) { return false; }
  if(excess.bytes() != 0)
  {
    // The exact-capacity allocator used by the supported toolchain takes the
    // usual route. Other implementations can use the bounded external path.
    return false;
  }
  for(size_type file = 0; file < this->graph_.files(); file++)
  {
    this->graph_.scanKMerBlocks(file, block_records * sizeof(KMer),
      [&keys, &starts, this](size_type, const std::vector<KMer>& kmers)
      {
        for(const KMer& kmer : kmers)
        {
          if(keys.size() == this->graph_.size())
          {
            throw std::runtime_error("resident preprocessing: input record count changed");
          }
          keys.push_back(kmer.key);
          starts.push_back(mapStart(kmer.from, this->graph_.mapping));
        }
      });
  }
  if(keys.size() != this->graph_.size())
  {
    throw std::runtime_error("resident preprocessing: input record count changed");
  }
  parallelQuickSort(keys.begin(), keys.end());
  size_type distinct = 0;
  for(size_type i = 0; i < keys.size(); i++)
  {
    if(distinct != 0 && Key::label(keys[distinct - 1]) == Key::label(keys[i]))
    {
      keys[distinct - 1] = Key::merge(keys[distinct - 1], keys[i]);
    }
    else { keys[distinct++] = keys[i]; }
  }
  keys.resize(distinct);
  for(key_type key : keys)
  {
    this->predecessor_edges_ += sdsl::bits::lt_cnt[Key::predecessors(key)];
  }
  removeDuplicates(starts, true);
  this->keys_.swap(keys); this->starts_.swap(starts);
  this->facts_memory_ = std::move(facts);
  this->key_count_ = this->keys_.size(); this->start_count_ = this->starts_.size();
  this->stats_.unique_keys = this->key_count_;
  this->stats_.unique_start_nodes = this->start_count_;
  this->resident_ = true;
  return true;
}

void
ExternalInputPreprocessor::materializeFacts()
{
  if(!this->resident_) { return; }
  this->key_name_ = TempFile::getName("gcsa_unique_keys");
  this->start_name_ = TempFile::getName("gcsa_unique_starts");
  for(size_type file = 0; file < 2; file++)
  {
    const std::string& name = (file == 0 ? this->key_name_ : this->start_name_);
    std::ofstream output;
    output.rdbuf()->pubsetbuf(nullptr, 0);
    output.open(name, std::ios_base::binary);
    if(!output) { throw std::runtime_error("resident preprocessing: cannot write " + name); }
    if(file == 0) { DiskIO::write(output, this->keys_.data(), this->keys_.size()); }
    else { DiskIO::write(output, this->starts_.data(), this->starts_.size()); }
    output.close();
    if(!output) { throw std::runtime_error("resident preprocessing: cannot close " + name); }
  }
  sdsl::util::clear(this->keys_); sdsl::util::clear(this->starts_);
  this->facts_memory_ = MemoryBudget::Reservation();
  this->resident_ = false;
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

namespace
{

size_type supportSum(size_type left, size_type right)
{
  if(left > std::numeric_limits<size_type>::max() - right)
  {
    throw std::runtime_error("preprocessing support size overflows");
  }
  return left + right;
}

size_type supportBytes(size_type count, size_type width)
{
  if(count > (std::numeric_limits<size_type>::max() - 16 * KILOBYTE) / width)
  {
    throw std::runtime_error("preprocessing support size overflows");
  }
  return count * width + 16 * KILOBYTE;
}

} // namespace

size_type
ExternalInputPreprocessor::lcpMemoryBytes()
{
  this->prepare();
  // The byte vector, its ram-file copy, wavelet construction and retained
  // bitmap/rank storage coexist. The constant covers the 256-symbol tables.
  return supportSum(supportBytes(this->key_count_, 8), this->ioBufferBytes());
}

size_type
ExternalInputPreprocessor::mapperMemoryBytes()
{
  this->prepare();
  if(this->predecessor_edges_ == 0 && this->key_count_ != 0)
  {
    // Old checkpoints did not record edge counts. Only that explicit resume
    // route needs this extra scan; ordinary reduction counts them as it writes.
    RawReader keys(this->key_name_, sizeof(key_type), this->key_count_, this->ioBufferBytes());
    key_type key = 0;
    while(keys.next(&key))
    {
      this->predecessor_edges_ += sdsl::bits::lt_cnt[Key::predecessors(key)];
    }
    keys.finish();
  }
  const size_type bwt_positions = (this->resident_ ? this->predecessor_edges_ : this->key_count_);
  const size_type bwt = supportBytes(bwt_positions, this->graph_.alpha.sigma);
  const size_type nodes = supportBytes(this->predecessor_edges_, 1);
  // Source bits and interleaved-rank copies overlap during assignment. Three
  // raw-sized copies conservatively cover the pinned bit_vector_il layout.
  return supportSum(supportBytes(supportSum(bwt, nodes) / 8 + 16, 3),
    this->ioBufferBytes());
}

size_type
ExternalInputPreprocessor::lastCharacterMemoryBytes()
{
  this->prepare();
  return supportSum(supportBytes(this->key_count_, 1), this->ioBufferBytes());
}

size_type
ExternalInputPreprocessor::startNodeMemoryBytes()
{
  this->prepare();
  // Elias--Fano low bits are at most 64 per start; high bits, select tables
  // and builder-to-vector overlap are included in the remaining allowance.
  return supportSum(supportBytes(this->start_count_, 24), this->ioBufferBytes());
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
  if(this->resident_)
  {
    mapper = DeBruijnGraph(this->keys_, this->graph_.k(), this->graph_.alpha);
    return;
  }
  const size_type io_bytes = this->ioBufferBytes();
  mapper = DeBruijnGraph(this->key_name_, this->key_count_, this->graph_.k(),
    this->graph_.alpha, io_bytes);
}

void
ExternalInputPreprocessor::buildLCP(LCP& lcp)
{
  this->prepare();
  if(this->resident_) { lcp = LCP(this->keys_, this->graph_.k()); return; }
  const size_type io_bytes = this->ioBufferBytes();
  lcp = LCP(this->key_name_, this->key_count_, this->graph_.k(), io_bytes);
}

void
ExternalInputPreprocessor::buildLastCharacters(sdsl::int_vector<0>& last_char)
{
  this->prepare();
  if(this->resident_) { Key::lastChars(this->keys_, last_char); return; }
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
  if(this->resident_)
  {
    from_nodes = sdsl::sd_vector<>(this->starts_.begin(), this->starts_.end());
    return;
  }
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

bool
ExternalInputPreprocessor::tryBuildResidentInitial(PathGraph& result)
{
  size_type largest_input = 0;
  for(size_type records : this->graph_.sizes) { largest_input = std::max(largest_input, records); }
  if(largest_input > std::numeric_limits<size_type>::max() / sizeof(InitialKMerRecord))
  {
    return false;
  }
  const size_type block_records = std::max(static_cast<size_type>(1),
    std::min(largest_input, this->parameters_.getIOBufferSize() / sizeof(KMer)));
  MemoryBudget::Reservation sorting = this->memory_->tryReserve(
    largest_input * sizeof(InitialKMerRecord), "initial logical KMer sort");
  MemoryBudget::Reservation decoding = this->memory_->tryReserve(
    block_records * sizeof(KMer), "initial decode block");
  if(!sorting || !decoding) { return false; }

  const size_type path_bytes = sizeof(PathNode) + 2 * sizeof(PathNode::rank_type);
  // Storage placement must not change the configured physical shard layout.
  const size_type configured_sort = std::min(this->parameters_.getMemoryLimitBytes(),
    this->parameters_.getSortRunSize());
  const size_type paths_per_shard = std::max(static_cast<size_type>(1),
    std::min(this->parameters_.getCheckpointRecords(),
      std::max(path_bytes, configured_sort / 2) / path_bytes));
  PathGraph initial(0, this->graph_.k(), 0, this->memory_);
  // Admit the complete initial generation alongside its largest logical sort
  // before doing work. A failed admission leaves no partial files or frontier.
  for(size_type file = 0; file < this->graph_.files(); file++)
  {
    if(file > std::numeric_limits<std::uint32_t>::max())
    {
      throw std::runtime_error("resident preprocessing: too many logical input files");
    }
    for(size_type offset = 0, shard = 0; offset < this->graph_.sizes[file]; shard++)
    {
      const size_type count = std::min(paths_per_shard, this->graph_.sizes[file] - offset);
      auto paths = std::make_shared<AdaptiveRecordStore<PathNode>>(
        *this->memory_, "initial paths");
      auto ranks = std::make_shared<AdaptiveRecordStore<PathNode::rank_type>>(
        *this->memory_, "initial ranks");
      if(!paths->tryReserveRecords(count) || !ranks->tryReserveRecords(2 * count))
      {
        return false;
      }
      paths->resize(count); ranks->resize(2 * count);
      initial.path_names.emplace_back(); initial.rank_names.emplace_back();
      initial.path_stores.push_back(std::move(paths));
      initial.rank_stores.push_back(std::move(ranks));
      initial.path_counts.push_back(count); initial.rank_counts.push_back(2 * count);
      initial.path_checksums.emplace_back(); initial.rank_checksums.emplace_back();
      logical_file_id_t logical(static_cast<std::uint32_t>(file));
      initial.logical_file_ids.push_back(logical);
      initial.physical_shard_ids.push_back(physicalId(logical, shard));
      initial.path_count += count; initial.rank_count += 2 * count;
      offset += count;
    }
  }

  size_type output_shard = 0;
  for(size_type file = 0; file < this->graph_.files(); file++)
  {
    std::vector<InitialKMerRecord> records(this->graph_.sizes[file]);
    size_type seen = 0;
    this->graph_.scanKMerBlocks(file, block_records * sizeof(KMer),
      [&records, &seen](size_type, const std::vector<KMer>& kmers)
      {
        for(const KMer& kmer : kmers)
        {
          if(seen == records.size())
          {
            throw std::runtime_error("resident preprocessing: logical input grew");
          }
          records[seen].kmer = kmer; records[seen].ordinal = seen; seen++;
        }
      });
    if(seen != records.size())
    {
      throw std::runtime_error("resident preprocessing: logical input shrank");
    }
    parallelQuickSort(records.begin(), records.end(),
      [](const InitialKMerRecord& a, const InitialKMerRecord& b)
      {
        return (compareInitialKMer(&a, &b) < 0);
      });
    size_type key_rank = 0, local = 0;
    for(const InitialKMerRecord& record : records)
    {
      const size_type label = Key::label(record.kmer.key);
      while(key_rank < this->keys_.size() && Key::label(this->keys_[key_rank]) < label)
      {
        key_rank++;
      }
      if(key_rank == this->keys_.size() || Key::label(this->keys_[key_rank]) != label)
      {
        throw std::runtime_error("resident preprocessing: KMer label absent from keys");
      }
      PathNode node;
      node.from = record.kmer.from; node.to = record.kmer.to; node.fields = 0;
      if(record.kmer.sorted()) { node.makeSorted(); }
      node.setPredecessors(Key::predecessors(record.kmer.key));
      node.setOrder(1); node.setLCP(1); node.setPointer(2 * local);
      initial.path_stores[output_shard]->residentData()[local] = node;
      auto* ranks = initial.rank_stores[output_shard]->residentData();
      ranks[2 * local] = static_cast<PathNode::rank_type>(key_rank);
      ranks[2 * local + 1] = 0;
      if(++local == initial.path_counts[output_shard]) { output_shard++; local = 0; }
    }
    this->stats_.logical_kmer_sorts++;
  }
  this->stats_.physical_shards = initial.files();
  initial.stored_bytes = 0;
  result.swap(initial);
  return true;
}

void
ExternalInputPreprocessor::buildInitialPathGraph(PathGraph& result)
{
  this->prepare();
  if(this->resident_ && this->tryBuildResidentInitial(result)) { return; }
  this->materializeFacts();
  result.clear();
  PathGraph initial(0, this->graph_.k(), 0);
  if(this->graph_.files() > std::numeric_limits<std::uint32_t>::max())
  {
    throw std::runtime_error("external preprocessing: too many logical input files");
  }
  const size_type sort_budget = this->sortBudget();
  const size_type io_bytes = std::min(this->ioBufferBytes(), sort_budget / 8);
  MemoryBudget::Reservation phase_memory;
  if(this->memory_ != nullptr)
  {
    phase_memory = this->memory_->tryReserve(sort_budget, "external initial preprocessing");
    if(!phase_memory) { throw std::runtime_error("cannot admit external initial sort"); }
  }
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

      // Mapping uses two bounded readers and two bounded path/rank writers,
      // each writer double-buffered for its helper thread. ioBufferBytes() is
      // at most sortBudget()/8, so their six explicit buffers stay within
      // three quarters of the sort budget. The subsequent succinct
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
