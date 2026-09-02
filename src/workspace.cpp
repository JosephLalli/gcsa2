#include <gcsa/workspace.h>
#include <algorithm>
#include <atomic>
#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <set>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>
namespace gcsa {
namespace {
const uint32_t HEAD=0x47534148U, FOOT=0x47534146U; const size_t BUFFER=1024*1024;
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
bool exists(const std::string& p) { return access(p.c_str(),F_OK)==0; }
std::string base(const std::string& p) { size_t x=p.find_last_of('/'); return x==std::string::npos?p:p.substr(x+1); }
std::string parent(const std::string& p) { size_t x=p.find_last_of('/'); return x==std::string::npos?".":(x==0?"/":p.substr(0,x)); }
std::string json(const std::string& s) { std::string r; for(size_t i=0;i<s.size();++i) { unsigned char c=s[i]; if(c=='"'||c=='\\') { r+='\\'; r+=c; } else if(c=='\n') r+="\\n"; else if(c=='\r') r+="\\r"; else if(c=='\t') r+="\\t"; else if(c<32) throw std::runtime_error("control character in manifest setting"); else r+=c; } return r; }
}
uint64_t BuildWorkspace::checksum(const void* d,size_t n,uint64_t h) { const uint8_t* p=(const uint8_t*)d; for(size_t i=0;i<n;i++) { h^=p[i]; h*=1099511628211ULL; } return h; }
std::string BuildWorkspace::safe(const std::string& v) { if(v.empty()) throw std::runtime_error("empty artifact identity component"); std::string out; for(size_t i=0;i<v.size();++i) { char c=v[i]; if((c>='a'&&c<='z')||(c>='A'&&c<='Z')||(c>='0'&&c<='9')||c=='-'||c=='_') out+=c; else out+='_'; } return out; }
BuildWorkspace::BuildWorkspace(const std::string& d,const Settings& s,const Settings& o,OpenMode mode):directory_(d),semantic_(s),operational_(o) { if(mkdir(d.c_str(),0755)!=0 && errno!=EEXIST) throw err("cannot create workspace",d); std::ostringstream sm; sm<<"gcsa-workspace-semantic-v2\n"; for(Settings::const_iterator i=s.begin();i!=s.end();++i) sm<<i->first<<"="<<i->second<<"\n"; std::string semantic=sm.str(); fingerprint_=std::to_string(checksum(semantic.data(),semantic.size())); std::string mp=d+"/build.json"; if(mode==NEW_WORKSPACE && exists(mp)) throw err("refusing NEW_WORKSPACE over existing manifest",mp); if(mode==RESUME && !exists(mp)) throw err("cannot resume workspace without a manifest",mp); if(exists(mp)) { int fd=open(mp.c_str(),O_RDONLY); if(fd<0) throw err("cannot read manifest",mp); std::string data; char b[4096]; for(;;) { ssize_t n=read(fd,b,sizeof(b)); if(n<0&&errno==EINTR)continue; if(n<0){close(fd);throw err("cannot read manifest",mp);} if(n==0)break; data.append(b,n); } close(fd); std::string key="\"semantic_fingerprint\":\""+fingerprint_+"\""; if(data.find(key)==std::string::npos) throw err("semantic fingerprint mismatch while resuming",mp); } std::string tmp=mp+partial_suffix(); int fd=open(tmp.c_str(),O_WRONLY|O_CREAT|O_EXCL,0644); if(fd<0) throw err("cannot create manifest",tmp); try { std::string m=manifest(); write_all(fd,m.data(),m.size(),tmp); sync_file(fd,tmp); if(rename(tmp.c_str(),mp.c_str())!=0) throw err("cannot publish manifest",mp); sync_dir(d); } catch(...) { close(fd); unlink(tmp.c_str()); throw; } recover(); }
std::string BuildWorkspace::manifest() const { std::ostringstream x; x<<"{\"version\":2,\"semantic_fingerprint\":\""<<fingerprint_<<"\",\"semantic\":{"; bool first=true; for(Settings::const_iterator i=semantic_.begin();i!=semantic_.end();++i){if(!first)x<<",";first=false;x<<"\""<<json(i->first)<<"\":\""<<json(i->second)<<"\"";} x<<"},\"operational\":{";first=true;for(Settings::const_iterator i=operational_.begin();i!=operational_.end();++i){if(!first)x<<",";first=false;x<<"\""<<json(i->first)<<"\":\""<<json(i->second)<<"\"";} x<<"}}\n"; return x.str(); }
std::string BuildWorkspace::artifact_path(const ArtifactIdentity& i,logical_file_id_t l,physical_shard_id_t s) const { std::ostringstream x; x<<directory_<<"/"<<safe(i.task)<<"--"<<safe(i.phase)<<"--"<<safe(i.relative_path)<<"--"<<safe(i.kind)<<"--"<<l.value<<"--"<<s.value<<".bin"; return x.str(); }
std::string BuildWorkspace::artifact_path(logical_file_id_t l,physical_shard_id_t s) const { return artifact_path(ArtifactIdentity(),l,s); }
std::string BuildWorkspace::completion_path(const std::string& t,const std::string& p) const { return directory_+"/"+safe(t)+"--"+safe(p)+".complete"; }
BuildWorkspace::ArtifactWriter::ArtifactWriter():ws_(0),fd_(-1),bytes_(0),sum_(1469598103934665603ULL),done_(true) {}
BuildWorkspace::ArtifactWriter::ArtifactWriter(BuildWorkspace* w,const ArtifactIdentity& i,logical_file_id_t l,physical_shard_id_t s,const std::string& so,const std::string& kr):ws_(w),id_(i),logical_(l),shard_(s),final_(w->artifact_path(i,l,s)),sort_(so),range_(kr),fd_(-1),bytes_(0),sum_(1469598103934665603ULL),done_(false) { partial_=final_+partial_suffix(); fd_=open(partial_.c_str(),O_CREAT|O_EXCL|O_WRONLY,0644); if(fd_<0) throw err("cannot create unique partial artifact",partial_); put(fd_,HEAD,4,partial_);put(fd_,ArtifactHeader::VERSION,4,partial_);put_string(fd_,id_.kind,partial_);put(fd_,logical_.value,4,partial_);put(fd_,shard_.value,8,partial_);put(fd_,0,8,partial_);put(fd_,0,8,partial_);put(fd_,0,8,partial_);put_string(fd_,sort_,partial_);put_string(fd_,range_,partial_); }
BuildWorkspace::ArtifactWriter::ArtifactWriter(ArtifactWriter&& o):ws_(o.ws_),id_(o.id_),logical_(o.logical_),shard_(o.shard_),final_(o.final_),partial_(o.partial_),sort_(o.sort_),range_(o.range_),fd_(o.fd_),bytes_(o.bytes_),sum_(o.sum_),done_(o.done_){o.fd_=-1;o.done_=true;}
BuildWorkspace::ArtifactWriter& BuildWorkspace::ArtifactWriter::operator=(ArtifactWriter&& o) { if(this!=&o){if(!done_){if(fd_>=0)close(fd_);if(!partial_.empty())unlink(partial_.c_str());}ws_=o.ws_;id_=o.id_;logical_=o.logical_;shard_=o.shard_;final_=o.final_;partial_=o.partial_;sort_=o.sort_;range_=o.range_;fd_=o.fd_;bytes_=o.bytes_;sum_=o.sum_;done_=o.done_;o.fd_=-1;o.done_=true;}return *this;}
BuildWorkspace::ArtifactWriter::~ArtifactWriter(){if(!done_){if(fd_>=0)close(fd_);if(!partial_.empty())unlink(partial_.c_str());}}
void BuildWorkspace::ArtifactWriter::write(const void* d,size_t n){if(done_)throw std::runtime_error("write on finished ArtifactWriter");write_all(fd_,d,n,partial_);sum_=BuildWorkspace::checksum(d,n,sum_);bytes_+=n;}
BuildWorkspace::ArtifactRef BuildWorkspace::ArtifactWriter::finish(uint64_t records){if(done_)throw std::runtime_error("ArtifactWriter already finished"); off_t count_off=off_t(24+id_.kind.size()); unsigned char b[8]; for(size_t i=0;i<8;i++)b[i]=unsigned(records>>(8*i));pwrite_all(fd_,b,8,count_off,partial_);for(size_t i=0;i<8;i++)b[i]=unsigned(bytes_>>(8*i));pwrite_all(fd_,b,8,count_off+8,partial_);for(size_t i=0;i<8;i++)b[i]=unsigned(sum_>>(8*i));pwrite_all(fd_,b,8,count_off+16,partial_);put(fd_,FOOT,4,partial_);put(fd_,ArtifactFooter::VERSION,4,partial_);put_string(fd_,id_.kind,partial_);put(fd_,logical_.value,4,partial_);put(fd_,shard_.value,8,partial_);put(fd_,records,8,partial_);put(fd_,bytes_,8,partial_);put(fd_,sum_,8,partial_);put_string(fd_,sort_,partial_);put_string(fd_,range_,partial_);sync_file(fd_,partial_);fd_=-1;if(rename(partial_.c_str(),final_.c_str())!=0)throw err("cannot atomically publish artifact",final_);sync_dir(ws_->directory_);done_=true;return ArtifactRef(id_,logical_,shard_,records,bytes_,sum_);}
BuildWorkspace::ArtifactWriter BuildWorkspace::open_artifact(const ArtifactIdentity& i,logical_file_id_t l,physical_shard_id_t s,const std::string& so,const std::string& kr){return ArtifactWriter(this,i,l,s,so,kr);}
void BuildWorkspace::commit_artifact(logical_file_id_t l,physical_shard_id_t s,uint64_t r,const std::vector<uint8_t>& p){ArtifactIdentity i;ArtifactWriter w=open_artifact(i,l,s);w.write(p.data(),p.size());ArtifactRef ref=w.finish(r);std::vector<ArtifactRef> a(1,ref);commit_task(i.task,i.phase,a);}
void BuildWorkspace::commit_task(const std::string& task,const std::string& phase,const std::vector<ArtifactRef>& a,const std::vector<std::string>& deps){if(a.empty())throw std::runtime_error("cannot commit empty task completion record");std::string final=completion_path(task,phase),tmp=final+partial_suffix();int fd=open(tmp.c_str(),O_CREAT|O_EXCL|O_WRONLY,0644);if(fd<0)throw err("cannot create completion record",tmp);try{std::ostringstream x;x<<"version=1\nfingerprint="<<fingerprint_<<"\n";for(size_t n=0;n<a.size();++n){std::string p=artifact_path(a[n].identity,a[n].logical,a[n].shard);if(!exists(p))throw err("completion references missing artifact",p);x<<"artifact\t"<<base(p)<<"\t"<<a[n].records<<"\t"<<a[n].bytes<<"\t"<<a[n].checksum<<"\n";}for(size_t n=0;n<deps.size();++n)x<<"dependency\t"<<deps[n]<<"\n";std::string text=x.str();write_all(fd,text.data(),text.size(),tmp);sync_file(fd,tmp);if(rename(tmp.c_str(),final.c_str())!=0)throw err("cannot publish completion record",final);sync_dir(directory_);}catch(...){close(fd);unlink(tmp.c_str());throw;}}
void BuildWorkspace::ensure_completed(const ArtifactIdentity& i,const std::string& b,uint64_t sum)const{std::string p=completion_path(i.task,i.phase);int fd=open(p.c_str(),O_RDONLY);if(fd<0)throw err("missing task completion record",p);std::string text;char buf[4096];for(;;){ssize_t n=read(fd,buf,sizeof(buf));if(n<0&&errno==EINTR)continue;if(n<0){close(fd);throw err("cannot read task completion record",p);}if(n==0)break;text.append(buf,n);}close(fd);std::string prefix="artifact\t"+b+"\t";size_t found=text.find(prefix);if(text.find("fingerprint="+fingerprint_+"\n")==std::string::npos||found==std::string::npos||text.find("\t"+std::to_string(sum)+"\n",found+prefix.size())==std::string::npos)throw err("artifact not committed by matching completion record",p);}
void BuildWorkspace::validate_artifact(const ArtifactIdentity& i,logical_file_id_t l,physical_shard_id_t s)const{std::string p=artifact_path(i,l,s);int fd=open(p.c_str(),O_RDONLY);if(fd<0)throw err("missing artifact",p);try{if(readn(fd,4,p)!=HEAD||readn(fd,4,p)!=ArtifactHeader::VERSION||get_string(fd,p)!=i.kind||readn(fd,4,p)!=l.value||readn(fd,8,p)!=s.value)throw err("artifact header mismatch",p);uint64_t rec=readn(fd,8,p),bytes=readn(fd,8,p),sum=readn(fd,8,p);std::string sort=get_string(fd,p),range=get_string(fd,p);std::vector<char>b(BUFFER);uint64_t left=bytes,actual=1469598103934665603ULL;while(left){size_t want=left<b.size()?size_t(left):b.size();size_t at=0;while(at<want){ssize_t n=read(fd,&b[at],want-at);if(n<0&&errno==EINTR)continue;if(n<=0)throw err("truncated artifact payload",p);at+=n;}actual=checksum(&b[0],want,actual);left-=want;}if(actual!=sum||readn(fd,4,p)!=FOOT||readn(fd,4,p)!=ArtifactFooter::VERSION||get_string(fd,p)!=i.kind||readn(fd,4,p)!=l.value||readn(fd,8,p)!=s.value||readn(fd,8,p)!=rec||readn(fd,8,p)!=bytes||readn(fd,8,p)!=sum||get_string(fd,p)!=sort||get_string(fd,p)!=range)throw err("artifact checksum or footer mismatch",p);char extra;ssize_t n=read(fd,&extra,1);if(n!=0)throw err("artifact has trailing data",p);close(fd);ensure_completed(i,base(p),sum);}catch(...){close(fd);throw;}}
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
    }
    close(fd);return result;
  }
  catch(...){close(fd);throw;}
}

void
BuildWorkspace::restore_artifact(const ArtifactIdentity& identity,
  logical_file_id_t logical,physical_shard_id_t shard,
  const std::string& output_path,size_t buffer_bytes) const
{
  this->validate_artifact(identity,logical,shard);
  if(buffer_bytes==0){throw std::invalid_argument("artifact restore buffer must be nonzero");}
  std::string source=this->artifact_path(identity,logical,shard);
  int input=open(source.c_str(),O_RDONLY);
  if(input<0){throw err("cannot read artifact",source);}
  std::string partial=output_path+partial_suffix();
  int output=open(partial.c_str(),O_WRONLY|O_CREAT|O_EXCL,0644);
  if(output<0){close(input);throw err("cannot create restored artifact",partial);}
  try
  {
    if(readn(input,4,source)!=HEAD||readn(input,4,source)!=ArtifactHeader::VERSION)
    {throw err("artifact header mismatch",source);}
    get_string(input,source);readn(input,4,source);readn(input,8,source);readn(input,8,source);
    uint64_t remaining=readn(input,8,source);readn(input,8,source);
    get_string(input,source);get_string(input,source);
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
      write_all(output,buffer.data(),want,partial);remaining-=want;
    }
    close(input);input=-1;sync_file(output,partial);output=-1;
    if(rename(partial.c_str(),output_path.c_str())!=0)
    {throw err("cannot publish restored artifact",output_path);}
    sync_dir(parent(output_path));
  }
  catch(...)
  {
    if(input>=0){close(input);}if(output>=0){close(output);}unlink(partial.c_str());throw;
  }
}

void BuildWorkspace::recover(){DIR*d=opendir(directory_.c_str());if(!d)throw err("cannot open workspace",directory_);std::set<std::string> keep;for(dirent*e;(e=readdir(d));){std::string n=e->d_name;if(n.size()>9&&n.substr(n.size()-9)==".complete"){int fd=open((directory_+"/"+n).c_str(),O_RDONLY);if(fd<0)continue;std::string t;char b[4096];ssize_t z;while((z=read(fd,b,sizeof(b)))>0)t.append(b,z);close(fd);std::istringstream in(t);std::string line;while(std::getline(in,line)){if(line.compare(0,9,"artifact\t")==0){size_t x=line.find('\t',9);if(x!=std::string::npos)keep.insert(line.substr(9,x-9));}}}}rewinddir(d);for(dirent*e;(e=readdir(d));){std::string n=e->d_name,p=directory_+"/"+n;if(n.find(".partial")!=std::string::npos)unlink(p.c_str());else if(n.size()>4&&n.substr(n.size()-4)==".bin"&&keep.find(n)==keep.end())unlink(p.c_str());}closedir(d);sync_dir(directory_);}
} // namespace gcsa
