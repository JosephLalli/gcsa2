#ifndef GCSA_WORKSPACE_H
#define GCSA_WORKSPACE_H

#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <vector>
namespace gcsa {
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
    BuildWorkspace* ws_; ArtifactIdentity id_; logical_file_id_t logical_; physical_shard_id_t shard_; std::string final_,partial_,sort_,range_; int fd_; uint64_t bytes_,sum_; bool done_;
  };
  BuildWorkspace(const std::string&,const Settings&,const Settings& operational=Settings(),OpenMode mode=RESUME);
  const std::string& fingerprint() const { return fingerprint_; } std::string manifest() const;
  std::string artifact_path(const ArtifactIdentity&,logical_file_id_t,physical_shard_id_t) const;
  std::string artifact_path(logical_file_id_t,physical_shard_id_t) const;
  ArtifactWriter open_artifact(const ArtifactIdentity&,logical_file_id_t,physical_shard_id_t,const std::string& sort_order="",const std::string& key_range="");
  void commit_artifact(logical_file_id_t,physical_shard_id_t,uint64_t,const std::vector<uint8_t>&);
  void commit_task(const std::string&,const std::string&,const std::vector<ArtifactRef>&,const std::vector<std::string>& dependencies=std::vector<std::string>());
  void validate_artifact(const ArtifactIdentity&,logical_file_id_t,physical_shard_id_t) const;
  void validate_artifact(logical_file_id_t,physical_shard_id_t) const;
  bool task_completed(const std::string& task,const std::string& phase) const;
  std::vector<uint8_t> read_artifact_payload(const ArtifactIdentity&,
    logical_file_id_t,physical_shard_id_t,size_t maximum_bytes) const;
  void restore_artifact(const ArtifactIdentity&,logical_file_id_t,
    physical_shard_id_t,const std::string& output_path,size_t buffer_bytes) const;
  void recover(); static uint64_t checksum(const void*,size_t,uint64_t seed=1469598103934665603ULL);
private:
  std::string directory_,fingerprint_; Settings semantic_,operational_;
  std::string completion_path(const std::string&,const std::string&) const; static std::string safe(const std::string&);
  void ensure_completed(const ArtifactIdentity&,const std::string&,uint64_t) const;
};
} // namespace gcsa
#endif
