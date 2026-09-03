#include <gcsa/resources.h>
#include <gcsa/workspace.h>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <sys/stat.h>
#include <sys/wait.h>
#include <thread>
#include <unistd.h>
using namespace gcsa;
static void require(bool v) { if(!v) { std::fprintf(stderr,"workspace test failure\n"); std::abort(); } }
static bool invalid(BuildWorkspace& w,const ArtifactIdentity& i) { try { w.validate_artifact(i,logical_file_id_t(7),physical_shard_id_t(9)); } catch(const std::runtime_error&) { return true; } return false; }
static bool restore_invalid(BuildWorkspace& w,const ArtifactIdentity& i,const std::string& output) { try { w.restore_artifact(i,logical_file_id_t(7),physical_shard_id_t(9),output,7777); } catch(const std::runtime_error&) { return access(output.c_str(),F_OK)!=0; } return false; }

static ArtifactIdentity crash_identity(const std::string& name)
{
  return ArtifactIdentity("crash-"+name,"commit","payload","crash-test");
}

static void crash_at(const std::string& root,const BuildWorkspace::Settings& semantic,
  const std::string& point,bool commit_task)
{
  pid_t child=fork();require(child>=0);
  if(child==0)
  {
    setenv("GCSA_WORKSPACE_CRASH_POINT",point.c_str(),1);
    BuildWorkspace workspace(root,semantic);
    ArtifactIdentity identity=crash_identity(point);
    BuildWorkspace::ArtifactWriter writer=workspace.open_artifact(
      identity,logical_file_id_t(11),physical_shard_id_t(13));
    const char payload[]="crash boundary";
    writer.write(payload,sizeof(payload));
    BuildWorkspace::ArtifactRef artifact=writer.finish(1);
    if(commit_task)
    {
      workspace.commit_task(identity.task,identity.phase,
        std::vector<BuildWorkspace::ArtifactRef>(1,artifact));
    }
    ::_exit(87); // The requested boundary was not reached.
  }
  int status=0;require(waitpid(child,&status,0)==child);
  require(WIFEXITED(status)&&WEXITSTATUS(status)==86);
}

static void crash_during_recovery(const std::string& root,
  const BuildWorkspace::Settings& semantic)
{
  pid_t child=fork();require(child>=0);
  if(child==0)
  {
    setenv("GCSA_WORKSPACE_CRASH_POINT","cleanup-after-remove",1);
    BuildWorkspace workspace(root,semantic);
    ::_exit(87);
  }
  int status=0;require(waitpid(child,&status,0)==child);
  require(WIFEXITED(status)&&WEXITSTATUS(status)==86);
}

int main() {
  MemoryBudget b(100,10); MemoryBudget::Reservation first=b.reserve(90,"first"); std::atomic<bool> started(false),got(false);
  std::thread waiter([&]{started=true;MemoryBudget::Reservation second=b.reserve(1,"second");got=true;}); while(!started) std::this_thread::yield(); std::this_thread::sleep_for(std::chrono::milliseconds(25)); require(!got); first=MemoryBudget::Reservation(); waiter.join(); require(got && b.stats().current==0 && b.stats().maximum==90); bool over=false; try { b.reserve(91,"too-large"); } catch(const std::runtime_error&) { over=true; } require(over);
  char root[]="/tmp/gcsa-workspace-XXXXXX"; require(mkdtemp(root)); BuildWorkspace::Settings sem,op; sem["k"]="64"; op["threads"]="2"; BuildWorkspace w(root,sem,op,BuildWorkspace::NEW_WORKSPACE); require(access((std::string(root)+"/build.json").c_str(),F_OK)==0);
  bool new_refused=false, semantic_refused=false; try { BuildWorkspace again(root,sem,op,BuildWorkspace::NEW_WORKSPACE); } catch(const std::runtime_error&) { new_refused=true; } BuildWorkspace::Settings changed=sem; changed["k"]="65"; try { BuildWorkspace bad(root,changed); } catch(const std::runtime_error&) { semantic_refused=true; } require(new_refused && semantic_refused); BuildWorkspace resumed(root,sem,BuildWorkspace::Settings());
  ArtifactIdentity id("phase/task","sort stage","relative/path","payload-kind"); require(!resumed.task_completed(id.task,id.phase)); std::vector<uint8_t> payload(2*1024*1024+17); for(size_t n=0;n<payload.size();++n) payload[n]=uint8_t(n); BuildWorkspace::ArtifactWriter writer=resumed.open_artifact(id,logical_file_id_t(7),physical_shard_id_t(9),"key","a-z"); writer.write(&payload[0],1024*1024); writer.write(&payload[1024*1024],payload.size()-1024*1024); BuildWorkspace::ArtifactRef ref=writer.finish(3); std::vector<BuildWorkspace::ArtifactRef> refs(1,ref); resumed.commit_task(id.task,id.phase,refs,std::vector<std::string>(1,"input-checksum")); resumed.validate_artifact(id,logical_file_id_t(7),physical_shard_id_t(9)); require(resumed.task_completed(id.task,id.phase)); require(resumed.read_artifact_payload(id,logical_file_id_t(7),physical_shard_id_t(9),payload.size())==payload); std::string restored=std::string(root)+"/restored.payload"; resumed.restore_artifact(id,logical_file_id_t(7),physical_shard_id_t(9),restored,7777); { std::ifstream input(restored.c_str(),std::ios::binary); std::vector<uint8_t> copy((std::istreambuf_iterator<char>(input)),std::istreambuf_iterator<char>()); require(copy==payload); } require(unlink(restored.c_str())==0);
  // Same-filesystem raw checkpoints are published and restored by hard link:
  // no payload copy, while the completion marker still governs publication.
  ArtifactIdentity raw_id("raw", "checkpoint", "payload", "raw-v1");
  std::string raw_source = std::string(root) + "/raw-source";
  {
    std::ofstream out(raw_source.c_str(), std::ios::binary);
    out.write(reinterpret_cast<const char*>(payload.data()), payload.size());
  }
  BuildWorkspace::ArtifactRef raw_ref = resumed.adopt_raw_payload(raw_id,
    logical_file_id_t(8), physical_shard_id_t(10), raw_source, 5,
    payload.size(), 7777);
  std::vector<BuildWorkspace::ArtifactRef> raw_refs(1, raw_ref);
  resumed.commit_task(raw_id.task, raw_id.phase, raw_refs);
  std::string raw_stored = resumed.artifact_path(raw_id,
    logical_file_id_t(8), physical_shard_id_t(10));
  std::string raw_restored = std::string(root) + "/raw-restored";
  struct stat source_stat, stored_stat, restored_stat;
  require(stat(raw_source.c_str(), &source_stat) == 0 &&
    stat(raw_stored.c_str(), &stored_stat) == 0 &&
    source_stat.st_ino == stored_stat.st_ino && stored_stat.st_nlink >= 2);
  resumed.restore_adopted_payload(raw_id, logical_file_id_t(8),
    physical_shard_id_t(10), raw_restored, 5, payload.size(), 7777);
  require(stat(raw_restored.c_str(), &restored_stat) == 0 &&
    restored_stat.st_ino == stored_stat.st_ino);
  require(unlink(raw_source.c_str()) == 0);
  BuildWorkspace resumed_raw(root, sem);
  resumed_raw.restore_adopted_payload(raw_id, logical_file_id_t(8),
    physical_shard_id_t(10), raw_source, 5, payload.size(), 7777);
  require(stat(raw_source.c_str(), &source_stat) == 0 &&
    source_stat.st_ino == stored_stat.st_ino);
  {
    std::fstream file(raw_stored.c_str(), std::ios::in | std::ios::out | std::ios::binary);
    file.seekp(0); char bad = 127; file.write(&bad, 1);
  }
  bool raw_bad = false;
  try
  {
    // Ordinary resume trusts the committed checksum and validates identity,
    // record count, and byte length. Explicit verification rereads payloads.
    resumed_raw.restore_adopted_payload(raw_id, logical_file_id_t(8),
      physical_shard_id_t(10), std::string(root) + "/bad-raw", 5,
      payload.size(), 7777, true);
  }
  catch(const std::runtime_error&) { raw_bad = true; }
  require(raw_bad);
  // Disk accounting follows allocated inode blocks rather than summing every
  // hardlink name as another full payload.
  std::string budget_source = std::string(root) + "/budget-source";
  std::string budget_link = std::string(root) + "/budget-link";
  {
    std::ofstream file(budget_source.c_str(), std::ios::binary);
    file.write(reinterpret_cast<const char*>(payload.data()), payload.size());
  }
  DiskBudget linked_budget(root, 100ULL * 1024 * 1024, 0);
  uint64_t before_link = linked_budget.available();
  require(link(budget_source.c_str(), budget_link.c_str()) == 0);
  require(linked_budget.available() == before_link);
  require(unlink(budget_link.c_str()) == 0);
  require(unlink(budget_source.c_str()) == 0);
  std::string good=resumed.artifact_path(id,logical_file_id_t(7),physical_shard_id_t(9)); std::string unmarked=resumed.artifact_path(ArtifactIdentity("other","phase","uncommitted","bin"),logical_file_id_t(1),physical_shard_id_t(2)); BuildWorkspace::ArtifactWriter orphan=resumed.open_artifact(ArtifactIdentity("other","phase","uncommitted","bin"),logical_file_id_t(1),physical_shard_id_t(2)); orphan.write(payload.data(),7); orphan.finish(1); std::string renamed=std::string(root)+"/renamed-successor.bin"; require(rename(unmarked.c_str(),renamed.c_str())==0); { std::ofstream stale((std::string(root)+"/stale.partial").c_str()); stale<<"partial"; } resumed.recover(); require(access(renamed.c_str(),F_OK)!=0 && access((std::string(root)+"/stale.partial").c_str(),F_OK)!=0); resumed.validate_artifact(id,logical_file_id_t(7),physical_shard_id_t(9));
  { std::fstream f(good.c_str(),std::ios::in|std::ios::out|std::ios::binary); f.seekp(54); char x=127; f.write(&x,1); } require(invalid(resumed,id)); require(restore_invalid(resumed,id,restored)); BuildWorkspace::ArtifactWriter repair=resumed.open_artifact(id,logical_file_id_t(7),physical_shard_id_t(9),"key","a-z"); repair.write(payload.data(),payload.size()); ref=repair.finish(3); refs[0]=ref; resumed.commit_task(id.task,id.phase,refs); { std::fstream f(good.c_str(),std::ios::in|std::ios::out|std::ios::binary); f.seekp(0,std::ios::end); std::streamoff n=f.tellp(); f.close(); require(truncate(good.c_str(),n-1)==0); } require(invalid(resumed,id)); require(restore_invalid(resumed,id,restored));
  // Exercise the actual durability boundaries with abrupt child exits. A
  // synced partial, a renamed-but-unmarked artifact, and a partial task marker
  // are all discarded. Once the task marker is renamed, the artifact survives.
  crash_at(root,sem,"artifact-before-rename",false); { BuildWorkspace recovered(root,sem); ArtifactIdentity crash=crash_identity("artifact-before-rename"); require(!recovered.task_completed(crash.task,crash.phase)&&access(recovered.artifact_path(crash,logical_file_id_t(11),physical_shard_id_t(13)).c_str(),F_OK)!=0); }
  crash_at(root,sem,"artifact-after-rename",false); { BuildWorkspace recovered(root,sem); ArtifactIdentity crash=crash_identity("artifact-after-rename"); require(!recovered.task_completed(crash.task,crash.phase)&&access(recovered.artifact_path(crash,logical_file_id_t(11),physical_shard_id_t(13)).c_str(),F_OK)!=0); }
  crash_at(root,sem,"task-before-rename",true); { BuildWorkspace recovered(root,sem); ArtifactIdentity crash=crash_identity("task-before-rename"); require(!recovered.task_completed(crash.task,crash.phase)&&access(recovered.artifact_path(crash,logical_file_id_t(11),physical_shard_id_t(13)).c_str(),F_OK)!=0); }
  crash_at(root,sem,"task-after-rename",true); { BuildWorkspace recovered(root,sem); ArtifactIdentity crash=crash_identity("task-after-rename"); require(recovered.task_completed(crash.task,crash.phase)); recovered.validate_artifact(crash,logical_file_id_t(11),physical_shard_id_t(13)); }

  // Cleanup itself is idempotent. Leave two renamed artifacts without task
  // markers, crash after removing the first, and let the next resume finish.
  ArtifactIdentity cleanup_a=crash_identity("cleanup-a"),cleanup_b=crash_identity("cleanup-b");
  BuildWorkspace::ArtifactWriter cleanup_writer_a=resumed.open_artifact(cleanup_a,logical_file_id_t(21),physical_shard_id_t(1)); cleanup_writer_a.write(payload.data(),19); cleanup_writer_a.finish(1);
  BuildWorkspace::ArtifactWriter cleanup_writer_b=resumed.open_artifact(cleanup_b,logical_file_id_t(21),physical_shard_id_t(2)); cleanup_writer_b.write(payload.data(),23); cleanup_writer_b.finish(1);
  crash_during_recovery(root,sem); { BuildWorkspace recovered(root,sem); require(access(recovered.artifact_path(cleanup_a,logical_file_id_t(21),physical_shard_id_t(1)).c_str(),F_OK)!=0&&access(recovered.artifact_path(cleanup_b,logical_file_id_t(21),physical_shard_id_t(2)).c_str(),F_OK)!=0); }
  DiskBudget disk(root,1,0); std::string why; require(!disk.can_reserve(2,&why) && !why.empty()); return 0;
}
