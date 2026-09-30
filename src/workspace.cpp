#include <gcsa/workspace.h>
#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>
namespace gcsa {
std::atomic<uint64_t> BuildWorkspace::adoption_checksum_scan_bytes(0);
std::atomic<uint64_t> BuildWorkspace::adoption_checksum_reused_bytes(0);
void
ClosedPayloadChecksum::record(const std::string& path, uint64_t digest)
{
  this->valid = false;
  if(::stat(path.c_str(), &this->identity) != 0 ||
     !S_ISREG(this->identity.st_mode) || this->identity.st_size < 0)
  {
    throw std::runtime_error("cannot seal closed payload checksum: " + path);
  }
  this->value = digest;
  this->valid = true;
}

bool
ClosedPayloadChecksum::matches(const std::string& path, uint64_t bytes) const
{
  struct stat current;
  if(!this->valid || ::stat(path.c_str(), &current) != 0 ||
     !S_ISREG(current.st_mode) || current.st_size < 0 ||
     static_cast<uint64_t>(current.st_size) != bytes) { return false; }
  return current.st_dev == this->identity.st_dev &&
    current.st_ino == this->identity.st_ino &&
    current.st_size == this->identity.st_size &&
    current.st_mtim.tv_sec == this->identity.st_mtim.tv_sec &&
    current.st_mtim.tv_nsec == this->identity.st_mtim.tv_nsec &&
    current.st_ctim.tv_sec == this->identity.st_ctim.tv_sec &&
    current.st_ctim.tv_nsec == this->identity.st_ctim.tv_nsec;
}

namespace {
const uint32_t HEAD=0x47534148U, FOOT=0x47534146U; const size_t BUFFER=1024*1024;
const off_t CACHE_TAIL_BYTES=64*1024*1024, CACHE_FLUSH_BYTES=512*1024*1024;
std::atomic<uint64_t> partial_counter(0);
std::string partial_suffix() { return "."+std::to_string((uint64_t)getpid())+"."+std::to_string(partial_counter.fetch_add(1))+".partial"; }
std::runtime_error err(const std::string& what,const std::string& path) { return std::runtime_error(what+": "+path); }
void write_all(int fd,const void* data,size_t n,const std::string& path) { const char* p=(const char*)data; while(n) { ssize_t r=write(fd,p,n); if(r<0 && errno==EINTR) continue; if(r<=0) throw err("write failed",path); p+=r; n-=r; } }
void pwrite_all(int fd,const void* data,size_t n,off_t off,const std::string& path) { const char* p=(const char*)data; while(n) { ssize_t r=pwrite(fd,p,n,off); if(r<0 && errno==EINTR) continue; if(r<=0) throw err("pwrite failed",path); p+=r; n-=r; off+=r; } }
void put(int fd,uint64_t v,size_t n,const std::string& p) { unsigned char b[8]; for(size_t i=0;i<n;i++) b[i]=unsigned(v>>(8*i)); write_all(fd,b,n,p); }
uint64_t readn(int fd,size_t n,const std::string& p) { unsigned char b[8]; size_t at=0; while(at<n) { ssize_t r=read(fd,b+at,n-at); if(r<0 && errno==EINTR) continue; if(r==0) throw err("truncated artifact",p); if(r<0) throw err("read failed",p); at+=r; } uint64_t v=0; for(size_t i=0;i<n;i++) v|=uint64_t(b[i])<<(8*i); return v; }
void put_string(int fd,const std::string& s,const std::string& p) { if(s.size()>1024*1024) throw std::runtime_error("artifact metadata exceeds 1 MiB"); put(fd,s.size(),4,p); write_all(fd,s.data(),s.size(),p); }
std::string get_string(int fd,const std::string& p) { uint64_t n=readn(fd,4,p); if(n>1024*1024) throw err("invalid oversized artifact metadata",p); std::string s(n,'\0'); size_t at=0; while(at<n) { ssize_t r=read(fd,&s[at],n-at); if(r<0 && errno==EINTR) continue; if(r<=0) throw err("truncated artifact metadata",p); at+=r; } return s; }
void sync_file(int fd,const std::string& p) { if(fdatasync(fd)!=0) throw err("fdatasync failed",p); if(close(fd)!=0) throw err("close failed",p); }
void sync_dir(const std::string& d) { int fd=open(d.c_str(),O_RDONLY|O_DIRECTORY); if(fd<0) throw err("open directory for fsync failed",d); if(fsync(fd)!=0) { close(fd); throw err("fsync directory failed",d); } if(close(fd)!=0) throw err("close directory failed",d); }
void advise_sequential(int fd) {
#if defined(POSIX_FADV_SEQUENTIAL)
  static_cast<void>(posix_fadvise(fd,0,0,POSIX_FADV_SEQUENTIAL));
#else
  static_cast<void>(fd);
#endif
}
void discard_cache(int fd,off_t begin,off_t bytes) {
#if defined(POSIX_FADV_DONTNEED)
  if(bytes>0){static_cast<void>(posix_fadvise(fd,begin,bytes,POSIX_FADV_DONTNEED));}
#else
  static_cast<void>(fd);static_cast<void>(begin);static_cast<void>(bytes);
#endif
}
void trim_read_cache(int fd,off_t consumed,off_t& released,bool complete=false) {
  if(!complete&&consumed-released<CACHE_FLUSH_BYTES)return;
  off_t end=(complete?consumed:std::max(released,consumed-CACHE_TAIL_BYTES));
  if(end>released){discard_cache(fd,released,end-released);released=end;}
}
void trim_written_cache(int fd,off_t& released,bool complete,const std::string& p) {
  off_t end=lseek(fd,0,SEEK_CUR);if(end<0)throw err("cannot determine artifact output position",p);
  if(!complete&&end-released<CACHE_FLUSH_BYTES)return;
  if(fdatasync(fd)!=0)throw err("fdatasync failed",p);
  off_t discard_end=(complete?end:std::max(released,end-CACHE_TAIL_BYTES));
  if(discard_end>released){discard_cache(fd,released,discard_end-released);released=discard_end;}
}
uint64_t checksum_raw_file(int fd,uint64_t bytes,std::vector<uint8_t>& buffer,
  const std::string& path) {
  advise_sequential(fd);off_t released=0;uint64_t result=1469598103934665603ULL;
  while(bytes>0) {
    size_t want=std::min<uint64_t>(bytes,buffer.size()),offset=0;
    while(offset<want) {
      ssize_t got=read(fd,buffer.data()+offset,want-offset);
      if(got<0&&errno==EINTR)continue;
      if(got<=0)throw err("truncated raw payload",path);
      offset+=got;
    }
    result=BuildWorkspace::checksum(buffer.data(),want,result);bytes-=want;
    off_t consumed=lseek(fd,0,SEEK_CUR);
    if(consumed<0)throw err("cannot determine raw payload input position",path);
    trim_read_cache(fd,consumed,released);
  }
  char extra;ssize_t got;
  do { got=read(fd,&extra,1); } while(got<0&&errno==EINTR);
  if(got<0)throw err("read failed",path);
  if(got>0)throw err("raw payload has trailing data",path);
  off_t consumed=lseek(fd,0,SEEK_CUR);
  if(consumed<0)throw err("cannot determine raw payload input position",path);
  trim_read_cache(fd,consumed,released,true);
  return result;
}
bool exists(const std::string& p) { return access(p.c_str(),F_OK)==0; }
std::string base(const std::string& p) { size_t x=p.find_last_of('/'); return x==std::string::npos?p:p.substr(x+1); }
std::string parent(const std::string& p) { size_t x=p.find_last_of('/'); return x==std::string::npos?".":(x==0?"/":p.substr(0,x)); }
std::string read_text_file(const std::string& path) {
  int fd=open(path.c_str(),O_RDONLY);if(fd<0)throw err("cannot read workspace record",path);
  std::string text;char buffer[4096];
  for(;;) { ssize_t n=read(fd,buffer,sizeof(buffer));if(n<0&&errno==EINTR)continue;if(n<0){close(fd);throw err("cannot read workspace record",path);}if(n==0)break;text.append(buffer,n); }
  if(close(fd)!=0)throw err("cannot close workspace record",path);
  return text;
}
std::set<std::string> committed_artifact_names(const std::string& path,
  const std::string& fingerprint,bool require_artifacts=true) {
  std::istringstream input(read_text_file(path));std::string line;bool version=false,match=false;std::set<std::string> names;
  while(std::getline(input,line)) {
    if(line=="version=1"){version=true;continue;}
    if(line=="fingerprint="+fingerprint){match=true;continue;}
    if(line.compare(0,9,"artifact\t")!=0)continue;
    std::istringstream fields(line);std::string tag,name;uint64_t records=0,bytes=0,sum=0;
    if(!(fields>>tag>>name>>records>>bytes>>sum) || tag!="artifact" || (fields>>std::ws && !fields.eof()) ||
      name!=base(name) || name.size()<4 || name.substr(name.size()-4)!=".bin")
    {throw err("corrupt task completion record",path);}
    if(require_artifacts) { struct stat st;if(stat((parent(path)+"/"+name).c_str(),&st)!=0||!S_ISREG(st.st_mode))throw err("completion references missing artifact",parent(path)+"/"+name); }
    names.insert(name);
  }
  if(!version||!match||names.empty())throw err("missing or corrupt successor task completion record",path);
  return names;
}
std::set<std::string> referenced_by_other_tasks(const std::string& directory,
  const std::string& exclude) {
  std::set<std::string> result;DIR* d=opendir(directory.c_str());if(d==0)throw err("cannot open workspace",directory);
  for(dirent* entry;(entry=readdir(d));) { std::string name=entry->d_name,path=directory+"/"+name;
    if(path==exclude||name.size()<=9||name.substr(name.size()-9)!=".complete")continue;
    std::istringstream input(read_text_file(path));std::string line;
    while(std::getline(input,line)) { if(line.compare(0,9,"artifact\t")!=0)continue;size_t end=line.find('\t',9);if(end!=std::string::npos)result.insert(line.substr(9,end-9)); }
  }
  closedir(d);return result;
}
struct RetirementRecord
{
  std::string predecessor,successor;
  bool complete;
};
RetirementRecord read_retirement_record(const std::string& path,
  const std::string& expected_fingerprint) {
  std::istringstream input(read_text_file(path));std::string line;
  bool version=false,fingerprint=false,complete=false;RetirementRecord result;
  while(std::getline(input,line)) {
    if(line=="version=1")version=true;
    else if(line=="fingerprint="+expected_fingerprint)fingerprint=true;
    else if(line.compare(0,12,"predecessor=")==0)result.predecessor=line.substr(12);
    else if(line.compare(0,10,"successor=")==0)result.successor=line.substr(10);
    else if(line=="complete=1")complete=true;
  }
  if(!version||!fingerprint||result.predecessor.empty()||result.successor.empty()||
    result.predecessor!=base(result.predecessor)||result.successor!=base(result.successor)||
    result.predecessor.size()<9||result.successor.size()<9||
    result.predecessor.substr(result.predecessor.size()-9)!=".complete"||
    result.successor.substr(result.successor.size()-9)!=".complete")
  {throw err("corrupt retirement record",path);}
  result.complete=complete;return result;
}
void complete_retirement_record(const std::string& path) {
  std::string text=read_text_file(path);
  if(text.find("complete=1\n")!=std::string::npos)return;
  if(text.empty()||text.back()!='\n')text+='\n';
  text+="complete=1\n";
  std::string temporary=path+partial_suffix();
  int fd=open(temporary.c_str(),O_CREAT|O_EXCL|O_WRONLY,0644);
  if(fd<0)throw err("cannot create completed retirement record",temporary);
  try {
    write_all(fd,text.data(),text.size(),temporary);
    if(fdatasync(fd)!=0)throw err("fdatasync failed",temporary);
    int completed_fd=fd;fd=-1;
    if(close(completed_fd)!=0)throw err("close failed",temporary);
    if(rename(temporary.c_str(),path.c_str())!=0)
      throw err("cannot complete retirement record",path);
    sync_dir(parent(path));
  } catch(...) { if(fd>=0)close(fd);unlink(temporary.c_str());throw; }
}
std::string json(const std::string& s) { std::string r; for(size_t i=0;i<s.size();++i) { unsigned char c=s[i]; if(c=='"'||c=='\\') { r+='\\'; r+=c; } else if(c=='\n') r+="\\n"; else if(c=='\r') r+="\\r"; else if(c=='\t') r+="\\t"; else if(c<32) throw std::runtime_error("control character in manifest setting"); else r+=c; } return r; }

// Deterministic crash injection for recovery tests. _exit() intentionally
// skips C++ destructors: the on-disk state must look like an abrupt process or
// machine failure, rather than an exception that lets ArtifactWriter clean its
// partial file. The variable is undocumented user interface and is only set by
// the workspace test executable.
void crash_if_requested(const char* point)
{
  const char* requested=std::getenv("GCSA_WORKSPACE_CRASH_POINT");
  if(requested!=0 && std::string(requested)==point){::_exit(86);}
}
}
uint64_t BuildWorkspace::checksum(const void* d,size_t n,uint64_t h) { const uint8_t* p=(const uint8_t*)d; for(size_t i=0;i<n;i++) { h^=p[i]; h*=1099511628211ULL; } return h; }
std::string BuildWorkspace::safe(const std::string& v) { if(v.empty()) throw std::runtime_error("empty artifact identity component"); std::string out; for(size_t i=0;i<v.size();++i) { char c=v[i]; if((c>='a'&&c<='z')||(c>='A'&&c<='Z')||(c>='0'&&c<='9')||c=='-'||c=='_') out+=c; else out+='_'; } return out; }
BuildWorkspace::BuildWorkspace(const std::string& d,const Settings& s,const Settings& o,OpenMode mode):directory_(d),semantic_(s),operational_(o) { if(mkdir(d.c_str(),0755)!=0 && errno!=EEXIST) throw err("cannot create workspace",d); std::ostringstream sm; sm<<"gcsa-workspace-semantic-v2\n"; for(Settings::const_iterator i=s.begin();i!=s.end();++i) sm<<i->first<<"="<<i->second<<"\n"; std::string semantic=sm.str(); fingerprint_=std::to_string(checksum(semantic.data(),semantic.size())); std::string mp=d+"/build.json"; if(mode==NEW_WORKSPACE && exists(mp)) throw err("refusing NEW_WORKSPACE over existing manifest",mp); if(mode==RESUME && !exists(mp)) throw err("cannot resume workspace without a manifest",mp); if(exists(mp)) { int fd=open(mp.c_str(),O_RDONLY); if(fd<0) throw err("cannot read manifest",mp); std::string data; char b[4096]; for(;;) { ssize_t n=read(fd,b,sizeof(b)); if(n<0&&errno==EINTR)continue; if(n<0){close(fd);throw err("cannot read manifest",mp);} if(n==0)break; data.append(b,n); } close(fd); std::string key="\"semantic_fingerprint\":\""+fingerprint_+"\""; if(data.find(key)==std::string::npos) throw err("semantic fingerprint mismatch while resuming",mp); } std::string tmp=mp+partial_suffix(); int fd=open(tmp.c_str(),O_WRONLY|O_CREAT|O_EXCL,0644); if(fd<0) throw err("cannot create manifest",tmp); try { std::string m=manifest(); write_all(fd,m.data(),m.size(),tmp); sync_file(fd,tmp); if(rename(tmp.c_str(),mp.c_str())!=0) throw err("cannot publish manifest",mp); sync_dir(d); } catch(...) { close(fd); unlink(tmp.c_str()); throw; } recover(); }
std::string BuildWorkspace::manifest() const { std::ostringstream x; x<<"{\"version\":2,\"semantic_fingerprint\":\""<<fingerprint_<<"\",\"semantic\":{"; bool first=true; for(Settings::const_iterator i=semantic_.begin();i!=semantic_.end();++i){if(!first)x<<",";first=false;x<<"\""<<json(i->first)<<"\":\""<<json(i->second)<<"\"";} x<<"},\"operational\":{";first=true;for(Settings::const_iterator i=operational_.begin();i!=operational_.end();++i){if(!first)x<<",";first=false;x<<"\""<<json(i->first)<<"\":\""<<json(i->second)<<"\"";} x<<"}}\n"; return x.str(); }
std::string BuildWorkspace::artifact_path(const ArtifactIdentity& i,logical_file_id_t l,physical_shard_id_t s) const { std::ostringstream x; x<<directory_<<"/"<<safe(i.task)<<"--"<<safe(i.phase)<<"--"<<safe(i.relative_path)<<"--"<<safe(i.kind)<<"--"<<l.value<<"--"<<s.value<<".bin"; return x.str(); }
std::string BuildWorkspace::artifact_path(logical_file_id_t l,physical_shard_id_t s) const { return artifact_path(ArtifactIdentity(),l,s); }
std::string BuildWorkspace::completion_path(const std::string& t,const std::string& p) const { return directory_+"/"+safe(t)+"--"+safe(p)+".complete"; }
std::string BuildWorkspace::retirement_path(const std::string& pt,const std::string& pp,const std::string& st,const std::string& sp) const { return directory_+"/retire--"+safe(pt)+"--"+safe(pp)+"--after--"+safe(st)+"--"+safe(sp)+".retired"; }
BuildWorkspace::ArtifactWriter::ArtifactWriter():ws_(0),fd_(-1),bytes_(0),sum_(1469598103934665603ULL),cache_released_(0),done_(true) {}
BuildWorkspace::ArtifactWriter::ArtifactWriter(BuildWorkspace* w,const ArtifactIdentity& i,logical_file_id_t l,physical_shard_id_t s,const std::string& so,const std::string& kr):ws_(w),id_(i),logical_(l),shard_(s),final_(w->artifact_path(i,l,s)),sort_(so),range_(kr),fd_(-1),bytes_(0),sum_(1469598103934665603ULL),cache_released_(0),done_(false) { partial_=final_+partial_suffix(); fd_=open(partial_.c_str(),O_CREAT|O_EXCL|O_WRONLY,0644); if(fd_<0) throw err("cannot create unique partial artifact",partial_); advise_sequential(fd_);put(fd_,HEAD,4,partial_);put(fd_,ArtifactHeader::VERSION,4,partial_);put_string(fd_,id_.kind,partial_);put(fd_,logical_.value,4,partial_);put(fd_,shard_.value,8,partial_);put(fd_,0,8,partial_);put(fd_,0,8,partial_);put(fd_,0,8,partial_);put_string(fd_,sort_,partial_);put_string(fd_,range_,partial_); }
BuildWorkspace::ArtifactWriter::ArtifactWriter(ArtifactWriter&& o):ws_(o.ws_),id_(o.id_),logical_(o.logical_),shard_(o.shard_),final_(o.final_),partial_(o.partial_),sort_(o.sort_),range_(o.range_),fd_(o.fd_),bytes_(o.bytes_),sum_(o.sum_),cache_released_(o.cache_released_),done_(o.done_){o.fd_=-1;o.done_=true;}
BuildWorkspace::ArtifactWriter& BuildWorkspace::ArtifactWriter::operator=(ArtifactWriter&& o) { if(this!=&o){if(!done_){if(fd_>=0)close(fd_);if(!partial_.empty())unlink(partial_.c_str());}ws_=o.ws_;id_=o.id_;logical_=o.logical_;shard_=o.shard_;final_=o.final_;partial_=o.partial_;sort_=o.sort_;range_=o.range_;fd_=o.fd_;bytes_=o.bytes_;sum_=o.sum_;cache_released_=o.cache_released_;done_=o.done_;o.fd_=-1;o.done_=true;}return *this;}
BuildWorkspace::ArtifactWriter::~ArtifactWriter(){if(!done_){if(fd_>=0)close(fd_);if(!partial_.empty())unlink(partial_.c_str());}}
void BuildWorkspace::ArtifactWriter::write(const void* d,size_t n){if(done_)throw std::runtime_error("write on finished ArtifactWriter");write_all(fd_,d,n,partial_);sum_=BuildWorkspace::checksum(d,n,sum_);bytes_+=n;trim_written_cache(fd_,cache_released_,false,partial_);}
BuildWorkspace::ArtifactRef BuildWorkspace::ArtifactWriter::finish(uint64_t records){if(done_)throw std::runtime_error("ArtifactWriter already finished"); off_t count_off=off_t(24+id_.kind.size()); unsigned char b[8]; for(size_t i=0;i<8;i++)b[i]=unsigned(records>>(8*i));pwrite_all(fd_,b,8,count_off,partial_);for(size_t i=0;i<8;i++)b[i]=unsigned(bytes_>>(8*i));pwrite_all(fd_,b,8,count_off+8,partial_);for(size_t i=0;i<8;i++)b[i]=unsigned(sum_>>(8*i));pwrite_all(fd_,b,8,count_off+16,partial_);put(fd_,FOOT,4,partial_);put(fd_,ArtifactFooter::VERSION,4,partial_);put_string(fd_,id_.kind,partial_);put(fd_,logical_.value,4,partial_);put(fd_,shard_.value,8,partial_);put(fd_,records,8,partial_);put(fd_,bytes_,8,partial_);put(fd_,sum_,8,partial_);put_string(fd_,sort_,partial_);put_string(fd_,range_,partial_);trim_written_cache(fd_,cache_released_,true,partial_);if(close(fd_)!=0)throw err("close failed",partial_);fd_=-1;crash_if_requested("artifact-before-rename");if(rename(partial_.c_str(),final_.c_str())!=0)throw err("cannot atomically publish artifact",final_);sync_dir(ws_->directory_);crash_if_requested("artifact-after-rename");done_=true;return ArtifactRef(id_,logical_,shard_,records,bytes_,sum_);}
BuildWorkspace::ArtifactWriter BuildWorkspace::open_artifact(const ArtifactIdentity& i,logical_file_id_t l,physical_shard_id_t s,const std::string& so,const std::string& kr){return ArtifactWriter(this,i,l,s,so,kr);}
void BuildWorkspace::commit_artifact(logical_file_id_t l,physical_shard_id_t s,uint64_t r,const std::vector<uint8_t>& p){ArtifactIdentity i;ArtifactWriter w=open_artifact(i,l,s);w.write(p.data(),p.size());ArtifactRef ref=w.finish(r);std::vector<ArtifactRef> a(1,ref);commit_task(i.task,i.phase,a);}
void BuildWorkspace::commit_task(const std::string& task,const std::string& phase,const std::vector<ArtifactRef>& a,const std::vector<std::string>& deps){if(a.empty())throw std::runtime_error("cannot commit empty task completion record");std::string final=completion_path(task,phase),tmp=final+partial_suffix();int fd=open(tmp.c_str(),O_CREAT|O_EXCL|O_WRONLY,0644);if(fd<0)throw err("cannot create completion record",tmp);try{std::ostringstream x;x<<"version=1\nfingerprint="<<fingerprint_<<"\n";for(size_t n=0;n<a.size();++n){std::string p=artifact_path(a[n].identity,a[n].logical,a[n].shard);if(!exists(p))throw err("completion references missing artifact",p);x<<"artifact\t"<<base(p)<<"\t"<<a[n].records<<"\t"<<a[n].bytes<<"\t"<<a[n].checksum<<"\n";}for(size_t n=0;n<deps.size();++n)x<<"dependency\t"<<deps[n]<<"\n";std::string text=x.str();write_all(fd,text.data(),text.size(),tmp);sync_file(fd,tmp);crash_if_requested("task-before-rename");if(rename(tmp.c_str(),final.c_str())!=0)throw err("cannot publish completion record",final);sync_dir(directory_);crash_if_requested("task-after-rename");}catch(...){close(fd);unlink(tmp.c_str());throw;}}
void BuildWorkspace::retire_marked(const std::string& predecessor,const std::string& successor)
{
  // Validate the successor before every unlink, including resume after a crash.
  committed_artifact_names(successor,fingerprint_);
  std::set<std::string> obsolete=committed_artifact_names(predecessor,fingerprint_,false);
  std::set<std::string> shared=referenced_by_other_tasks(directory_,predecessor);
  for(std::set<std::string>::const_iterator i=obsolete.begin();i!=obsolete.end();++i)
  {
    if(shared.find(*i)!=shared.end())continue;
    std::string path=directory_+"/"+*i;
    if(unlink(path.c_str())==0)crash_if_requested("retire-after-remove");
    else if(errno!=ENOENT)throw err("cannot retire obsolete artifact",path);
  }
  sync_dir(directory_);
}
void BuildWorkspace::retire_obsolete(const std::string& predecessor_task,
  const std::string& predecessor_phase,const std::string& successor_task,
  const std::string& successor_phase)
{
  const std::string predecessor=completion_path(predecessor_task,predecessor_phase);
  const std::string successor=completion_path(successor_task,successor_phase);
  const std::string final=retirement_path(predecessor_task,predecessor_phase,successor_task,successor_phase);
  if(exists(final))
  {
    RetirementRecord record=read_retirement_record(final,fingerprint_);
    if(record.predecessor!=base(predecessor)||record.successor!=base(successor))
      throw err("retirement record identity mismatch",final);
    if(record.complete)return;
    // The marker is the durable deletion intent. Its replay path accepts
    // already-missing predecessor artifacts, making direct retries idempotent.
    retire_marked(predecessor,successor);
    complete_retirement_record(final);
    return;
  }
  // Do not publish a recovery journal until both task markers are trustworthy.
  committed_artifact_names(successor,fingerprint_);
  // The successor is the authoritative frontier. Accept an already-partially
  // retired predecessor so cleanup can recover even if an older deletion
  // journal was lost after one or more unlinks.
  std::set<std::string> predecessor_artifacts=
    committed_artifact_names(predecessor,fingerprint_,false);
  bool remains=false;
  for(const std::string& artifact:predecessor_artifacts)
  {
    if(exists(directory_+"/"+artifact)){remains=true;break;}
  }
  if(!remains)return;
  const std::string temporary=final+partial_suffix();int fd=open(temporary.c_str(),O_CREAT|O_EXCL|O_WRONLY,0644);
  if(fd<0)throw err("cannot create retirement record",temporary);
  try {
    std::string text="version=1\nfingerprint="+fingerprint_+"\npredecessor="+base(predecessor)+"\nsuccessor="+base(successor)+"\n";
    write_all(fd,text.data(),text.size(),temporary);sync_file(fd,temporary);fd=-1;
    crash_if_requested("retire-before-marker-rename");
    if(rename(temporary.c_str(),final.c_str())!=0)throw err("cannot publish retirement record",final);
    sync_dir(directory_);
  } catch(...) { if(fd>=0)close(fd);unlink(temporary.c_str());throw; }
  retire_marked(predecessor,successor);
  complete_retirement_record(final);
}

void
BuildWorkspace::retire_obsolete_family(
  const std::string& predecessor_task_prefix,
  const std::string& predecessor_phase,
  const std::string& successor_task,
  const std::string& successor_phase)
{
  const std::string prefix=safe(predecessor_task_prefix);
  const std::string suffix="--"+safe(predecessor_phase)+".complete";
  std::vector<std::string> tasks;
  DIR* directory=opendir(directory_.c_str());
  if(directory==0)throw err("cannot open workspace",directory_);
  for(dirent* entry;(entry=readdir(directory));)
  {
    std::string name=entry->d_name;
    if(name.size()<prefix.size()+suffix.size() ||
       name.compare(0,prefix.size(),prefix)!=0 ||
       name.compare(name.size()-suffix.size(),suffix.size(),suffix)!=0)
    {
      continue;
    }
    tasks.push_back(name.substr(0,name.size()-suffix.size()));
  }
  closedir(directory);
  std::sort(tasks.begin(),tasks.end());
  tasks.erase(std::unique(tasks.begin(),tasks.end()),tasks.end());
  for(const std::string& task:tasks)
  {
    if(task==successor_task && predecessor_phase==successor_phase)continue;
    const std::string predecessor=completion_path(task,predecessor_phase);
    std::set<std::string> artifacts=committed_artifact_names(
      predecessor,fingerprint_,false);
    bool remains=false;
    for(const std::string& artifact:artifacts)
    {
      if(exists(directory_+"/"+artifact)){remains=true;break;}
    }
    // A prior retirement journal may already have removed this task. Avoid
    // publishing a redundant marker for a predecessor with nothing left.
    if(!remains)continue;
    this->retire_obsolete(task,predecessor_phase,
      successor_task,successor_phase);
  }
}
uint64_t
BuildWorkspace::committed_artifact_checksum(const ArtifactIdentity& identity,
  const std::string& basename, uint64_t records, uint64_t bytes) const
{
  const std::string path = completion_path(identity.task, identity.phase);
  int descriptor = open(path.c_str(), O_RDONLY);
  if(descriptor < 0) { throw err("missing task completion record", path); }
  std::string text;
  char buffer[4096];
  for(;;)
  {
    ssize_t got = read(descriptor, buffer, sizeof(buffer));
    if(got < 0 && errno == EINTR) { continue; }
    if(got < 0)
    {
      close(descriptor); throw err("cannot read task completion record", path);
    }
    if(got == 0) { break; }
    text.append(buffer, got);
  }
  if(close(descriptor) != 0) { throw err("cannot close task completion record", path); }
  if(text.find("fingerprint=" + fingerprint_ + "\n") == std::string::npos)
  {
    throw err("task completion fingerprint mismatch", path);
  }

  std::istringstream input(text);
  std::string line;
  while(std::getline(input, line))
  {
    std::istringstream fields(line);
    std::string tag, name;
    uint64_t committed_records = 0, committed_bytes = 0, committed_checksum = 0;
    if(fields >> tag >> name >> committed_records >> committed_bytes >> committed_checksum &&
       tag == "artifact" && name == basename)
    {
      const uint64_t any = std::numeric_limits<uint64_t>::max();
      if((records != any && records != committed_records) ||
         (bytes != any && bytes != committed_bytes))
      {
        throw err("artifact count or length disagrees with completion record", path);
      }
      return committed_checksum;
    }
  }
  throw err("artifact not committed by matching completion record", path);
}

void
BuildWorkspace::ensure_completed(const ArtifactIdentity& identity,
  const std::string& basename, uint64_t checksum) const
{
  const uint64_t any = std::numeric_limits<uint64_t>::max();
  if(committed_artifact_checksum(identity, basename, any, any) != checksum)
  {
    throw err("artifact checksum disagrees with completion record",
      completion_path(identity.task, identity.phase));
  }
}
void
BuildWorkspace::validate_artifact(const ArtifactIdentity& i,logical_file_id_t l,
  physical_shard_id_t s) const
{
  std::string p=artifact_path(i,l,s);int fd=open(p.c_str(),O_RDONLY);
  if(fd<0)throw err("missing artifact",p);
  advise_sequential(fd);off_t cache_released=0;
  try
  {
    if(readn(fd,4,p)!=HEAD||readn(fd,4,p)!=ArtifactHeader::VERSION||
      get_string(fd,p)!=i.kind||readn(fd,4,p)!=l.value||readn(fd,8,p)!=s.value)
    {throw err("artifact header mismatch",p);}
    uint64_t rec=readn(fd,8,p),bytes=readn(fd,8,p),sum=readn(fd,8,p);
    std::string sort=get_string(fd,p),range=get_string(fd,p);
    std::vector<char>b(BUFFER);uint64_t left=bytes,actual=1469598103934665603ULL;
    while(left)
    {
      size_t want=left<b.size()?size_t(left):b.size(),at=0;
      while(at<want)
      {
        ssize_t n=read(fd,&b[at],want-at);if(n<0&&errno==EINTR)continue;
        if(n<=0){throw err("truncated artifact payload",p);}
        at+=n;
      }
      actual=checksum(&b[0],want,actual);left-=want;
      off_t consumed=lseek(fd,0,SEEK_CUR);
      if(consumed<0)throw err("cannot determine artifact input position",p);
      trim_read_cache(fd,consumed,cache_released);
    }
    if(actual!=sum||readn(fd,4,p)!=FOOT||readn(fd,4,p)!=ArtifactFooter::VERSION||
      get_string(fd,p)!=i.kind||readn(fd,4,p)!=l.value||readn(fd,8,p)!=s.value||
      readn(fd,8,p)!=rec||readn(fd,8,p)!=bytes||readn(fd,8,p)!=sum||
      get_string(fd,p)!=sort||get_string(fd,p)!=range)
    {throw err("artifact checksum or footer mismatch",p);}
    char extra;ssize_t n=read(fd,&extra,1);
    if(n!=0)throw err("artifact has trailing data",p);
    off_t consumed=lseek(fd,0,SEEK_CUR);
    if(consumed<0)throw err("cannot determine artifact input position",p);
    trim_read_cache(fd,consumed,cache_released,true);
    if(close(fd)!=0){throw err("close failed",p);}
    fd=-1;
    ensure_completed(i,base(p),sum);
  }
  catch(...){if(fd>=0)close(fd);throw;}
}
void BuildWorkspace::validate_artifact(logical_file_id_t l,physical_shard_id_t s)const{validate_artifact(ArtifactIdentity(),l,s);}

bool
BuildWorkspace::task_completed(const std::string& task,const std::string& phase) const
{
  std::string path=this->completion_path(task,phase);
  int fd=open(path.c_str(),O_RDONLY);
  if(fd<0)
  {
    if(errno==ENOENT){return false;}
    throw err("cannot read task completion record",path);
  }
  std::string text;char buffer[4096];
  for(;;)
  {
    ssize_t bytes=read(fd,buffer,sizeof(buffer));
    if(bytes<0&&errno==EINTR){continue;}
    if(bytes<0){close(fd);throw err("cannot read task completion record",path);}
    if(bytes==0){break;}
    text.append(buffer,bytes);
  }
  close(fd);
  return text.find("fingerprint="+this->fingerprint_+"\n")!=std::string::npos;
}

std::vector<uint8_t>
BuildWorkspace::read_artifact_payload(const ArtifactIdentity& identity,
  logical_file_id_t logical,physical_shard_id_t shard,size_t maximum_bytes) const
{
  this->validate_artifact(identity,logical,shard);
  std::string path=this->artifact_path(identity,logical,shard);
  int fd=open(path.c_str(),O_RDONLY);
  if(fd<0){throw err("cannot read artifact",path);}
  advise_sequential(fd);off_t cache_released=0;
  try
  {
    if(readn(fd,4,path)!=HEAD||readn(fd,4,path)!=ArtifactHeader::VERSION)
    {throw err("artifact header mismatch",path);}
    get_string(fd,path);readn(fd,4,path);readn(fd,8,path);readn(fd,8,path);
    uint64_t bytes=readn(fd,8,path);readn(fd,8,path);get_string(fd,path);get_string(fd,path);
    if(bytes>maximum_bytes){throw err("artifact payload exceeds bounded read limit",path);}
    std::vector<uint8_t> result(bytes);
    size_t offset=0;
    while(offset<result.size())
    {
      ssize_t got=read(fd,result.data()+offset,result.size()-offset);
      if(got<0&&errno==EINTR){continue;}
      if(got<=0){throw err("truncated artifact payload",path);}
      offset+=got;
      off_t consumed=lseek(fd,0,SEEK_CUR);
      if(consumed<0){throw err("cannot determine artifact input position",path);}
      trim_read_cache(fd,consumed,cache_released);
    }
    off_t consumed=lseek(fd,0,SEEK_CUR);
    if(consumed<0){throw err("cannot determine artifact input position",path);}
    trim_read_cache(fd,consumed,cache_released,true);
    if(close(fd)!=0){throw err("close failed",path);}fd=-1;return result;
  }
  catch(...){if(fd>=0){close(fd);}throw;}
}

void
BuildWorkspace::restore_artifact(const ArtifactIdentity& identity,
  logical_file_id_t logical,physical_shard_id_t shard,
  const std::string& output_path,size_t buffer_bytes) const
{
  if(buffer_bytes==0){throw std::invalid_argument("artifact restore buffer must be nonzero");}
  std::string source=this->artifact_path(identity,logical,shard);
  int input=open(source.c_str(),O_RDONLY);
  if(input<0){throw err("cannot read artifact",source);}
  std::string partial=output_path+partial_suffix();
  int output=open(partial.c_str(),O_WRONLY|O_CREAT|O_EXCL,0644);
  if(output<0){close(input);throw err("cannot create restored artifact",partial);}
  advise_sequential(input);advise_sequential(output);
  off_t input_cache_released=0,output_cache_released=0;
  try
  {
    if(readn(input,4,source)!=HEAD||readn(input,4,source)!=ArtifactHeader::VERSION||
      get_string(input,source)!=identity.kind||readn(input,4,source)!=logical.value||
      readn(input,8,source)!=shard.value)
    {throw err("artifact header mismatch",source);}
    uint64_t records=readn(input,8,source);
    uint64_t bytes=readn(input,8,source);
    uint64_t expected_checksum=readn(input,8,source);
    std::string sort_order=get_string(input,source);
    std::string key_range=get_string(input,source);
    uint64_t remaining=bytes;
    uint64_t actual_checksum=1469598103934665603ULL;
    std::vector<uint8_t> buffer(buffer_bytes);
    while(remaining>0)
    {
      size_t want=std::min<uint64_t>(remaining,buffer.size()),offset=0;
      while(offset<want)
      {
        ssize_t got=read(input,buffer.data()+offset,want-offset);
        if(got<0&&errno==EINTR){continue;}
        if(got<=0){throw err("truncated artifact payload",source);}
        offset+=got;
      }
      actual_checksum=checksum(buffer.data(),want,actual_checksum);
      write_all(output,buffer.data(),want,partial);remaining-=want;
      off_t consumed=lseek(input,0,SEEK_CUR);
      if(consumed<0){throw err("cannot determine artifact input position",source);}
      trim_read_cache(input,consumed,input_cache_released);
      trim_written_cache(output,output_cache_released,false,partial);
    }
    if(actual_checksum!=expected_checksum||
      readn(input,4,source)!=FOOT||readn(input,4,source)!=ArtifactFooter::VERSION||
      get_string(input,source)!=identity.kind||readn(input,4,source)!=logical.value||
      readn(input,8,source)!=shard.value||readn(input,8,source)!=records||
      readn(input,8,source)!=bytes||readn(input,8,source)!=expected_checksum||
      get_string(input,source)!=sort_order||get_string(input,source)!=key_range)
    {throw err("artifact checksum or footer mismatch",source);}
    char extra;
    ssize_t extra_bytes;
    do { extra_bytes=read(input,&extra,1); } while(extra_bytes<0&&errno==EINTR);
    if(extra_bytes<0){throw err("read failed",source);}
    if(extra_bytes>0){throw err("artifact has trailing data",source);}
    off_t consumed=lseek(input,0,SEEK_CUR);
    if(consumed<0){throw err("cannot determine artifact input position",source);}
    trim_read_cache(input,consumed,input_cache_released,true);
    if(close(input)!=0){throw err("close failed",source);}input=-1;
    // A valid header/footer is not sufficient: only an artifact named by the
    // matching durable task marker may be restored and published.
    this->ensure_completed(identity,base(source),expected_checksum);
    trim_written_cache(output,output_cache_released,true,partial);
    if(close(output)!=0){throw err("close failed",partial);}output=-1;
    if(rename(partial.c_str(),output_path.c_str())!=0)
    {throw err("cannot publish restored artifact",output_path);}
    sync_dir(parent(output_path));
  }
  catch(...)
  {
    if(input>=0){close(input);}if(output>=0){close(output);}unlink(partial.c_str());throw;
  }
}

BuildWorkspace::ArtifactRef
BuildWorkspace::adopt_raw_payload(const ArtifactIdentity& identity,
  logical_file_id_t logical, physical_shard_id_t shard,
  const std::string& source, uint64_t records, uint64_t expected_bytes,
  size_t buffer_bytes, const uint64_t* known_checksum)
{
  if(buffer_bytes == 0)
  {
    throw std::invalid_argument("raw payload buffer must be nonzero");
  }
  struct stat st;
  if(stat(source.c_str(), &st) != 0 || !S_ISREG(st.st_mode) ||
     st.st_size < 0 || static_cast<uint64_t>(st.st_size) != expected_bytes)
  {
    throw err("raw payload length mismatch", source);
  }

  const std::string target = artifact_path(identity, logical, shard);
  const std::string partial = target + partial_suffix();
  std::vector<uint8_t> buffer;
  const auto allocate_buffer = [&]()
  {
    if(buffer.empty()) { buffer.resize(std::max<size_t>(1, buffer_bytes)); }
  };
  int input = -1, output = -1;
  uint64_t sum = 1469598103934665603ULL;
  try
  {
    if(link(source.c_str(), partial.c_str()) == 0)
    {
      // The writer that produced `source` is closed before checkpointing.
      // Syncing through the new link makes its data durable without copying.
      input = open(partial.c_str(), O_RDONLY);
      if(input < 0) { throw err("cannot read linked raw payload", partial); }
      if(fdatasync(input) != 0) { throw err("fdatasync failed", partial); }
      // A writer that incrementally checksummed exactly the closed payload can
      // make adoption metadata-only. Callers without such provenance retain
      // the full validation scan used by the original API.
      if(known_checksum == nullptr)
      {
        allocate_buffer();
        sum = checksum_raw_file(input, expected_bytes, buffer, partial);
        adoption_checksum_scan_bytes += expected_bytes;
      }
      else
      {
        sum = *known_checksum;
        adoption_checksum_reused_bytes += expected_bytes;
      }
      if(close(input) != 0) { throw err("close failed", partial); }
      input = -1;
    }
    else
    {
      // Cross-filesystem and filesystems without hardlinks retain a bounded
      // copy fallback. Compute the checksum during that one unavoidable pass.
      input = open(source.c_str(), O_RDONLY);
      output = open(partial.c_str(), O_WRONLY | O_CREAT | O_EXCL, 0644);
      if(input < 0 || output < 0)
      {
        throw err("cannot create raw payload fallback", partial);
      }
      advise_sequential(input); advise_sequential(output);
      allocate_buffer();
      off_t input_released = 0, output_released = 0;
      uint64_t left = expected_bytes;
      while(left > 0)
      {
        size_t want = std::min<uint64_t>(left, buffer.size()), offset = 0;
        while(offset < want)
        {
          ssize_t got = read(input, buffer.data() + offset, want - offset);
          if(got < 0 && errno == EINTR) { continue; }
          if(got <= 0) { throw err("truncated raw payload", source); }
          offset += got;
        }
        sum = checksum(buffer.data(), want, sum);
        adoption_checksum_scan_bytes += want;
        write_all(output, buffer.data(), want, partial); left -= want;
        off_t consumed = lseek(input, 0, SEEK_CUR);
        if(consumed < 0) { throw err("cannot determine raw payload input position", source); }
        trim_read_cache(input, consumed, input_released);
        trim_written_cache(output, output_released, false, partial);
      }
      char extra; ssize_t got;
      do { got = read(input, &extra, 1); } while(got < 0 && errno == EINTR);
      if(got < 0) { throw err("read failed", source); }
      if(got > 0) { throw err("raw payload has trailing data", source); }
      if(known_checksum != nullptr && sum != *known_checksum)
      {
        throw err("raw payload checksum differs from writer checksum", source);
      }
      off_t consumed = lseek(input, 0, SEEK_CUR);
      if(consumed < 0) { throw err("cannot determine raw payload input position", source); }
      trim_read_cache(input, consumed, input_released, true);
      if(close(input) != 0) { throw err("close failed", source); }
      input = -1;
      trim_written_cache(output, output_released, true, partial);
      if(close(output) != 0) { throw err("close failed", partial); }
      output = -1;
    }

    crash_if_requested("artifact-before-rename");
    if(rename(partial.c_str(), target.c_str()) != 0)
    {
      throw err("cannot atomically publish raw payload", target);
    }
    sync_dir(directory_);
    crash_if_requested("artifact-after-rename");
  }
  catch(...)
  {
    if(input >= 0) { close(input); }
    if(output >= 0) { close(output); }
    unlink(partial.c_str()); throw;
  }
  return ArtifactRef(identity, logical, shard, records, expected_bytes, sum);
}

void
BuildWorkspace::restore_adopted_payload(const ArtifactIdentity& identity,
  logical_file_id_t logical, physical_shard_id_t shard,
  const std::string& output, uint64_t records, uint64_t expected_bytes,
  size_t buffer_bytes, bool verify_checksum) const
{
  if(buffer_bytes == 0)
  {
    throw std::invalid_argument("raw payload buffer must be nonzero");
  }
  const std::string source = artifact_path(identity, logical, shard);
  struct stat st;
  if(stat(source.c_str(), &st) != 0 || !S_ISREG(st.st_mode) || st.st_size < 0 ||
     static_cast<uint64_t>(st.st_size) != expected_bytes)
  {
    throw err("raw payload length mismatch", source);
  }
  const uint64_t expected_checksum = committed_artifact_checksum(
    identity, base(source), records, expected_bytes);
  const std::string partial = output + partial_suffix();
  std::vector<uint8_t> buffer(std::max<size_t>(1, buffer_bytes));
  int input = -1, destination = -1;
  try
  {
    if(link(source.c_str(), partial.c_str()) == 0)
    {
      // Normal resume is O(metadata): committed length and task identity are
      // enough. --verify-workspace explicitly requests the full checksum pass.
      if(verify_checksum)
      {
        input = open(partial.c_str(), O_RDONLY);
        if(input < 0) { throw err("cannot verify linked raw payload", partial); }
        uint64_t actual = checksum_raw_file(input, expected_bytes, buffer, partial);
        if(close(input) != 0) { throw err("close failed", partial); }
        input = -1;
        if(actual != expected_checksum) { throw err("raw payload checksum mismatch", source); }
      }
    }
    else
    {
      input = open(source.c_str(), O_RDONLY);
      destination = open(partial.c_str(), O_WRONLY | O_CREAT | O_EXCL, 0644);
      if(input < 0 || destination < 0)
      {
        throw err("cannot create restored raw payload", partial);
      }
      advise_sequential(input); advise_sequential(destination);
      uint64_t actual = 1469598103934665603ULL, left = expected_bytes;
      off_t input_released = 0, output_released = 0;
      while(left > 0)
      {
        size_t want = std::min<uint64_t>(left, buffer.size()), offset = 0;
        while(offset < want)
        {
          ssize_t got = read(input, buffer.data() + offset, want - offset);
          if(got < 0 && errno == EINTR) { continue; }
          if(got <= 0) { throw err("truncated raw payload", source); }
          offset += got;
        }
        actual = checksum(buffer.data(), want, actual);
        write_all(destination, buffer.data(), want, partial); left -= want;
        off_t consumed = lseek(input, 0, SEEK_CUR);
        if(consumed < 0) { throw err("cannot determine raw payload input position", source); }
        trim_read_cache(input, consumed, input_released);
        trim_written_cache(destination, output_released, false, partial);
      }
      if(actual != expected_checksum) { throw err("raw payload checksum mismatch", source); }
      if(close(input) != 0) { throw err("close failed", source); }
      input = -1;
      trim_written_cache(destination, output_released, true, partial);
      if(close(destination) != 0) { throw err("close failed", partial); }
      destination = -1;
    }
    if(rename(partial.c_str(), output.c_str()) != 0)
    {
      throw err("cannot publish restored raw payload", output);
    }
    sync_dir(parent(output));
  }
  catch(...)
  {
    if(input >= 0) { close(input); }
    if(destination >= 0) { close(destination); }
    unlink(partial.c_str()); throw;
  }
}

void
BuildWorkspace::recover()
{
  DIR* directory=opendir(directory_.c_str());
  if(directory==0){throw err("cannot open workspace",directory_);}

  // A renamed retirement record is a durable promise to finish deleting the
  // predecessor names. Revalidate its successor before doing so on every
  // resume; a missing or damaged successor leaves the predecessor untouched.
  for(dirent* entry;(entry=readdir(directory));)
  {
    std::string name=entry->d_name,path=directory_+"/"+name;
    if(name.size()<=8||name.substr(name.size()-8)!=".retired")continue;
    RetirementRecord record;
    try { record=read_retirement_record(path,fingerprint_); }
    catch(...) { closedir(directory);throw; }
    if(record.complete)continue;
    retire_marked(directory_+"/"+record.predecessor,
      directory_+"/"+record.successor);
    complete_retirement_record(path);
  }

  // A completion record is the commit point. Renamed artifacts without one
  // are intentionally excluded, because a crash may have happened after an
  // artifact rename but before all outputs of its task were durable.
  std::set<std::string> keep;
  rewinddir(directory);
  for(dirent* entry;(entry=readdir(directory));)
  {
    std::string name=entry->d_name;
    if(name.size()<=9||name.substr(name.size()-9)!=".complete"){continue;}
    int descriptor=open((directory_+"/"+name).c_str(),O_RDONLY);
    if(descriptor<0){continue;}
    std::string text;char buffer[4096];ssize_t bytes;
    while((bytes=read(descriptor,buffer,sizeof(buffer)))>0){text.append(buffer,bytes);}
    close(descriptor);
    std::istringstream input(text);std::string line;
    while(std::getline(input,line))
    {
      if(line.compare(0,9,"artifact\t")!=0){continue;}
      size_t separator=line.find('\t',9);
      if(separator!=std::string::npos){keep.insert(line.substr(9,separator-9));}
    }
  }

  rewinddir(directory);
  for(dirent* entry;(entry=readdir(directory));)
  {
    std::string name=entry->d_name,path=directory_+"/"+name;
    // The framed writer uses a bounded on-disk footer-index sidecar. Neither
    // it nor a unique partial is reusable after a crash without a task marker.
    bool partial=(name.find(".partial")!=std::string::npos ||
      name.find(".index.")!=std::string::npos);
    bool uncommitted=(name.size()>4&&name.substr(name.size()-4)==".bin"&&
      keep.find(name)==keep.end());
    if(partial||uncommitted)
    {
      if(unlink(path.c_str())!=0&&errno!=ENOENT)
      {
        closedir(directory);throw err("cannot remove incomplete workspace artifact",path);
      }
      crash_if_requested("cleanup-after-remove");
    }
  }
  closedir(directory);sync_dir(directory_);
}
} // namespace gcsa
