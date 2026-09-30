#ifndef GCSA_WORKSPACE_H
#define GCSA_WORKSPACE_H

#include <cstddef>
#include <atomic>
#include <cstdint>
#include <map>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>
namespace gcsa {
// In-process provenance from an immutable, successfully closed writer. This is
// not persisted in workspace manifests: reopening an old workspace uses its
// recorded checksums and normal restore validation, not a remembered pathname.
// A live child writer can transmit it in its ephemeral worker-result record.
struct ClosedPayloadChecksum
{
  uint64_t value = 0;
  bool valid = false;
  struct stat identity = {};

  void record(const std::string& path, uint64_t digest);
  bool matches(const std::string& path, uint64_t bytes) const;
};

struct logical_file_id_t
{
  explicit logical_file_id_t(uint32_t x = 0) : value(x) { }
  bool operator==(const logical_file_id_t& another) const { return value == another.value; }
  bool operator!=(const logical_file_id_t& another) const { return value != another.value; }
  bool operator<(const logical_file_id_t& another) const { return value < another.value; }
  uint32_t value;
};

struct physical_shard_id_t
{
  explicit physical_shard_id_t(uint64_t x = 0) : value(x) { }
  bool operator==(const physical_shard_id_t& another) const { return value == another.value; }
  bool operator!=(const physical_shard_id_t& another) const { return value != another.value; }
  bool operator<(const physical_shard_id_t& another) const { return value < another.value; }
  uint64_t value;
};
struct ArtifactIdentity { std::string task, phase, relative_path, kind; ArtifactIdentity(const std::string& t="default",const std::string& p="default",const std::string& r="artifact",const std::string& k="binary"):task(t),phase(p),relative_path(r),kind(k){} };
struct ArtifactHeader { static const uint32_t VERSION=2; std::string kind,sort_order,key_range; logical_file_id_t logical; physical_shard_id_t shard; uint64_t records,bytes,checksum; };
struct ArtifactFooter { static const uint32_t VERSION=2; std::string kind,sort_order,key_range; logical_file_id_t logical; physical_shard_id_t shard; uint64_t records,bytes,checksum; };
class BuildWorkspace {
public:
  // Adoption diagnostics only; no effect on persisted metadata or semantics.
  static std::atomic<uint64_t> adoption_checksum_scan_bytes;
  static std::atomic<uint64_t> adoption_checksum_reused_bytes;
  typedef std::map<std::string,std::string> Settings;
  enum OpenMode { RESUME, NEW_WORKSPACE };
  struct ArtifactRef
  {
    ArtifactIdentity identity;
    logical_file_id_t logical;
    physical_shard_id_t shard;
    uint64_t records, bytes, checksum;

    ArtifactRef(const ArtifactIdentity& i, logical_file_id_t l,
      physical_shard_id_t s, uint64_t r, uint64_t b, uint64_t c) :
      identity(i), logical(l), shard(s), records(r), bytes(b), checksum(c) { }
  };
  class ArtifactWriter {
  public:
    ArtifactWriter(); ArtifactWriter(ArtifactWriter&&); ArtifactWriter& operator=(ArtifactWriter&&); ~ArtifactWriter();
    void write(const void*,size_t); ArtifactRef finish(uint64_t records);
  private:
    friend class BuildWorkspace; ArtifactWriter(BuildWorkspace*,const ArtifactIdentity&,logical_file_id_t,physical_shard_id_t,const std::string&,const std::string&);
    ArtifactWriter(const ArtifactWriter&); ArtifactWriter& operator=(const ArtifactWriter&);
    BuildWorkspace* ws_; ArtifactIdentity id_; logical_file_id_t logical_; physical_shard_id_t shard_; std::string final_,partial_,sort_,range_; int fd_; uint64_t bytes_,sum_; off_t cache_released_; bool done_;
  };
  BuildWorkspace(const std::string&,const Settings&,const Settings& operational=Settings(),OpenMode mode=RESUME);
  const std::string& fingerprint() const { return fingerprint_; } std::string manifest() const;
  std::string artifact_path(const ArtifactIdentity&,logical_file_id_t,physical_shard_id_t) const;
  std::string artifact_path(logical_file_id_t,physical_shard_id_t) const;
  ArtifactWriter open_artifact(const ArtifactIdentity&,logical_file_id_t,physical_shard_id_t,const std::string& sort_order="",const std::string& key_range="");
  void commit_artifact(logical_file_id_t,physical_shard_id_t,uint64_t,const std::vector<uint8_t>&);
  void commit_task(const std::string&,const std::string&,const std::vector<ArtifactRef>&,const std::vector<std::string>& dependencies=std::vector<std::string>());
  // Durably retire artifacts published by predecessor only after successor's
  // completion marker is present and structurally valid. Safe to repeat after
  // a crash; the retirement marker is retained as the recovery journal.
  void retire_obsolete(const std::string& predecessor_task,
    const std::string& predecessor_phase,const std::string& successor_task,
    const std::string& successor_phase);
  // Retire every committed task in one deterministic task-name family. This
  // is used after a complete extend checkpoint supersedes its per-partition
  // join outputs and sampled partition plan.
  void retire_obsolete_family(const std::string& predecessor_task_prefix,
    const std::string& predecessor_phase,const std::string& successor_task,
    const std::string& successor_phase);
  void validate_artifact(const ArtifactIdentity&,logical_file_id_t,physical_shard_id_t) const;
  void validate_artifact(logical_file_id_t,physical_shard_id_t) const;
  bool task_completed(const std::string& task,const std::string& phase) const;
  std::vector<uint8_t> read_artifact_payload(const ArtifactIdentity&,
    logical_file_id_t,physical_shard_id_t,size_t maximum_bytes) const;
  void restore_artifact(const ArtifactIdentity&,logical_file_id_t,
    physical_shard_id_t,const std::string& output_path,size_t buffer_bytes) const;
  // known_checksum may be supplied only by a closed immutable writer that
  // hashed exactly the raw payload while writing it. Same-filesystem adoption
  // can then avoid rereading the payload; copy fallback still verifies it.
  ArtifactRef adopt_raw_payload(const ArtifactIdentity&,logical_file_id_t,
    physical_shard_id_t,const std::string& source_path,uint64_t records,
    uint64_t expected_bytes,size_t buffer_bytes,
    const uint64_t* known_checksum=nullptr);
  void restore_adopted_payload(const ArtifactIdentity&,logical_file_id_t,
    physical_shard_id_t,const std::string& output_path,uint64_t records,
    uint64_t expected_bytes,size_t buffer_bytes,bool verify_checksum=false) const;
  void recover(); static uint64_t checksum(const void*,size_t,uint64_t seed=1469598103934665603ULL);
private:
  std::string directory_,fingerprint_; Settings semantic_,operational_;
  std::string completion_path(const std::string&,const std::string&) const;
  std::string retirement_path(const std::string&,const std::string&,const std::string&,const std::string&) const;
  void retire_marked(const std::string&,const std::string&);
  static std::string safe(const std::string&);
  void ensure_completed(const ArtifactIdentity&,const std::string&,uint64_t) const;
  uint64_t committed_artifact_checksum(const ArtifactIdentity&,
    const std::string&,uint64_t records,uint64_t bytes) const;
};
} // namespace gcsa
#endif
