// Separate prospective scanner. The frozen displayed-score engine remains unchanged.
// Reuse its checked cache framing, integer parser and exact heap ordering.
#define main displayed_score_reference_main
#include "scan_joint_cache_annulus.cpp"
#undef main

struct LocalTable {
 int p,m;
 std::vector<std::int64_t> units;
 std::vector<unsigned char> good;
};
static std::vector<LocalTable> read_local(const char* path) {
 std::ifstream in(path,std::ios::binary);char marker[8];exact(in,marker,8);
 if(std::string(marker,8)!="MW16LC01" || u32(in)!=2)throw std::runtime_error("invalid local cache frame");
 std::vector<LocalTable> out;
 for(int k=0;k<2;++k) {
  int p=static_cast<int>(u32(in)),m=static_cast<int>(u32(in));auto length=u32(in);
  if(p!=(k==0?5:13) || m!=(k==0?125:169) || length!=static_cast<unsigned>(m+m/p))throw std::runtime_error("invalid local projective ring");
  LocalTable t{p,m,std::vector<std::int64_t>(length),std::vector<unsigned char>(length)};
  std::vector<unsigned char> bytes(std::size_t(length)*8);exact(in,reinterpret_cast<char*>(bytes.data()),bytes.size());exact(in,reinterpret_cast<char*>(t.good.data()),length);
  for(std::size_t i=0;i<length;++i) {
   std::uint64_t v=0;for(int j=0;j<8;++j)v|=std::uint64_t(bytes[8*i+j])<<(8*j);
   t.units[i]=(v>>63)?-1-static_cast<std::int64_t>(~v):static_cast<std::int64_t>(v);
   if(t.good[i]>1 || (!t.good[i] && t.units[i]) || t.units[i]<-10000000000000LL || t.units[i]>10000000000000LL)throw std::runtime_error("invalid local score symbol");
  }
  out.push_back(std::move(t));
 }
 exact(in,marker,8);if(std::string(marker,8)!="ENDLC001" || in.peek()!=std::char_traits<char>::eof())throw std::runtime_error("invalid local cache end");
 return out;
}
int main(int argc,char**argv) {
 try {
  if(argc!=12)throw std::runtime_error("usage: corrected_scan SHORT_CACHE EXT_CACHE SIGN NUM DEN KEEP SHARD SHARDS INNER PRIME_COUNT LOCAL_CACHE");
  int sign=integer(argv[3],-1,1);if(!sign)throw std::runtime_error("nonzero sign required");
  const int cap=16777216;
  int N=integer(argv[4],1,cap),M=integer(argv[5],1,cap),K=integer(argv[6],1,1048576),shards=integer(argv[8],1,cap),shard=integer(argv[7],0,shards-1),inner=integer(argv[9],0,std::max(N,M)-1);
  int expected_primes=integer(argv[10],1,8192);
  std::vector<Table> tables;read_cache(argv[1],tables);read_cache(argv[2],tables);
  if(static_cast<int>(tables.size())!=expected_primes)throw std::runtime_error("declared prime count differs");
  auto local=read_local(argv[11]);
  for(int q:{5,13})if(std::none_of(tables.begin(),tables.end(),[q](const Table&t){return t.p==q;}))throw std::runtime_error("corrected prime absent from base score");
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
   for(const auto&t:local) {
    // Projective coordinates over Z/p^k: use n/d if d is a unit,
    // otherwise d/n with n a unit. Nonunits in both coordinates fail gcd below.
    int exponent=t.m-t.m/t.p-1;
    int inverse=d%t.p?power(d%t.m,exponent,t.m):0;
    std::vector<std::int64_t> cycle(t.m);std::vector<int> flags(t.m);
    for(int j=0;j<t.m;++j) {
     int n=sign>0?j:(t.m-j)%t.m,index;
     if(inverse)index=std::int64_t(n)*inverse%t.m;
     else {
      if(n%t.p==0)continue;
      int s=std::int64_t(d%t.m)*power(n,exponent,t.m)%t.m;
      if(s%t.p)throw std::runtime_error("invalid infinity coordinate");
      index=t.m+s/t.p;
     }
     cycle[j]=t.units[index];flags[j]=t.good[index];
    }
    for(std::size_t offset=0;offset<size;offset+=t.m) {
     std::size_t count=std::min(std::size_t(t.m),size-offset);
     for(std::size_t j=0;j<count;++j){scores[offset+j]+=cycle[j];goods[offset+j]+=flags[j];}
    }
   }
   for(int n=1;n<=N;++n) {
    if(std::max(n,d)<=inner || std::gcd(n,d)!=1)continue;
    ++population;Row r{n,d,scores[n]+constant,goods[n]+constant_good};
    if(r.good>expected_primes)throw std::runtime_error("local correction double-counted good primes");
    if(static_cast<int>(heap.size())<K)heap.push(r);else if(better(r,heap.top())){heap.pop();heap.push(r);}
   }
  }
  std::vector<Row> rows;while(!heap.empty()){rows.push_back(heap.top());heap.pop();}std::sort(rows.begin(),rows.end(),better);
  std::cout<<"CORRECTED_MW16_ANNULUS_V1\nP "<<sign<<' '<<N<<' '<<M<<' '<<K<<' '<<shard<<' '<<shards<<' '<<inner<<' '<<tables.size()<<'\n';
  for(const auto&r:rows)std::cout<<"C "<<sign*r.n<<' '<<r.d<<' '<<r.units<<' '<<r.good<<'\n';
  std::cout<<"S "<<population<<' '<<rows.size()<<'\n';return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
