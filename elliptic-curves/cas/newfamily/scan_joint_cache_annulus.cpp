// Complete combined-score annulus enumeration before any candidate retention.
// Two immutable, disjoint R17XS001 caches; signed projective residues.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>
struct Table { int p; std::vector<std::int64_t> units; std::vector<unsigned char> good; };
struct Row { int n,d; std::int64_t units; int good; };
static void exact(std::istream& in,char* out,std::size_t n) {
 if(!in.read(out,static_cast<std::streamsize>(n))) throw std::runtime_error("short cache frame");
}
static std::uint32_t u32(std::istream& in) {
 unsigned char b[4];exact(in,reinterpret_cast<char*>(b),4);
 return std::uint32_t(b[0])|(std::uint32_t(b[1])<<8)|(std::uint32_t(b[2])<<16)|(std::uint32_t(b[3])<<24);
}
static bool prime(int p) { for(int q=2;q*q<=p;++q)if(p%q==0)return false;return p>=2; }
static int integer(const char* arg,int lo,int hi) {
 std::size_t used=0;std::string s=arg;auto n=std::stoll(s,&used);
 if(used!=s.size() || n<lo || n>hi)throw std::runtime_error("invalid integer argument");return static_cast<int>(n);
}
static void read_cache(const char* path,std::vector<Table>& tables) {
 std::ifstream in(path,std::ios::binary);char marker[8];exact(in,marker,8);
 if(std::string(marker,8)!="R17XS001")throw std::runtime_error("invalid cache magic");
 auto count=u32(in);if(count<1 || count>8192 || tables.size()+count>8192)throw std::runtime_error("invalid prime count");
 int previous=tables.empty()?0:tables.back().p;
 for(std::uint32_t j=0;j<count;++j) {
  int p=static_cast<int>(u32(in));auto length=u32(in);
  if(p<=previous || p<5 || p>32749 || !prime(p) || length!=static_cast<std::uint32_t>(p+1))throw std::runtime_error("invalid or overlapping prime frame");
  previous=p;Table t{p,std::vector<std::int64_t>(length),std::vector<unsigned char>(length)};
  std::vector<unsigned char> bytes(std::size_t(length)*8);exact(in,reinterpret_cast<char*>(bytes.data()),bytes.size());exact(in,reinterpret_cast<char*>(t.good.data()),length);
  for(std::size_t i=0;i<length;++i) {
   std::uint64_t v=0;for(int k=0;k<8;++k)v|=std::uint64_t(bytes[8*i+k])<<(8*k);
   t.units[i]=(v>>63)?-1-static_cast<std::int64_t>(~v):static_cast<std::int64_t>(v);
   if(t.good[i]>1 || (!t.good[i] && t.units[i]) || t.units[i]<-10000000000000LL || t.units[i]>10000000000000LL)throw std::runtime_error("invalid score symbol");
  }
  tables.push_back(std::move(t));
 }
 exact(in,marker,8);if(std::string(marker,8)!="ENDXSC01" || in.peek()!=std::char_traits<char>::eof())throw std::runtime_error("invalid cache end");
}
static bool better(const Row& a,const Row& b) {
 if(a.units!=b.units)return a.units>b.units;
 if(a.good!=b.good)return a.good>b.good;
 if(a.d!=b.d)return a.d<b.d;
 return a.n<b.n; // n is absolute; signed slice order is explicit in the protocol.
}
struct Better { bool operator()(const Row&a,const Row&b)const{return better(a,b);} };
static int power(int a,int n,int p) {
 int r=1;while(n){if(n&1)r=std::int64_t(r)*a%p;a=std::int64_t(a)*a%p;n>>=1;}return r;
}
int main(int argc,char** argv) {
 try {
  if(argc!=11)throw std::runtime_error("usage: joint_scan SHORT_CACHE EXT_CACHE SIGN NUM DEN KEEP SHARD SHARDS INNER");
  int sign=integer(argv[3],-1,1);if(!sign)throw std::runtime_error("nonzero sign required");
  const int cap=16777216;
  int N=integer(argv[4],1,cap),M=integer(argv[5],1,cap),K=integer(argv[6],1,1048576),shards=integer(argv[8],1,cap),shard=integer(argv[7],0,shards-1),inner=integer(argv[9],0,std::max(N,M)-1);
  // argv[10] is a required declared prime count, binding both cache bands.
  int expected_primes=integer(argv[10],1,8192);
  std::vector<Table> tables;read_cache(argv[1],tables);read_cache(argv[2],tables);
  if(static_cast<int>(tables.size())!=expected_primes)throw std::runtime_error("declared prime count differs");
  std::priority_queue<Row,std::vector<Row>,Better> heap;std::uint64_t population=0;
  for(int d=shard+1;d<=M;d+=shards) {
   std::size_t size=std::size_t(N)+1;std::vector<std::int64_t> scores(size);std::vector<int> goods(size);std::int64_t constant=0;int constant_good=0;
   for(const auto&t:tables) {
    int p=t.p;
    if(d%p==0){constant+=t.units[p];constant_good+=t.good[p];continue;}
    int inverse=power(d%p,p-2,p);if(sign<0)inverse=p-inverse;
    std::vector<std::int64_t> cycle(p);std::vector<int> flags(p);int residue=0;
    for(int j=0;j<p;++j){cycle[j]=t.units[residue];flags[j]=t.good[residue];residue+=inverse;if(residue>=p)residue-=p;}
    for(std::size_t offset=0;offset<size;offset+=p) {
     std::size_t count=std::min(std::size_t(p),size-offset);
     for(std::size_t j=0;j<count;++j){scores[offset+j]+=cycle[j];goods[offset+j]+=flags[j];}
    }
   }
   for(int n=1;n<=N;++n) {
    if(std::max(n,d)<=inner || std::gcd(n,d)!=1)continue;
    ++population;Row r{n,d,scores[n]+constant,goods[n]+constant_good};
    if(static_cast<int>(heap.size())<K)heap.push(r);else if(better(r,heap.top())){heap.pop();heap.push(r);}
   }
  }
  std::vector<Row> rows;while(!heap.empty()){rows.push_back(heap.top());heap.pop();}std::sort(rows.begin(),rows.end(),better);
  std::cout<<"JOINT_NAGAO_ANNULUS_V1\nP "<<sign<<' '<<N<<' '<<M<<' '<<K<<' '<<shard<<' '<<shards<<' '<<inner<<' '<<tables.size()<<'\n';
  for(const auto&r:rows)std::cout<<"C "<<sign*r.n<<' '<<r.d<<' '<<r.units<<' '<<r.good<<'\n';
  std::cout<<"S "<<population<<' '<<rows.size()<<'\n';return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
