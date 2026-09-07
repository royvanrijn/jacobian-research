// Score an immutable rational-candidate list from quantized projective tables.
// No point counting, parameter enumeration, ranking or candidate refill occurs.
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
struct Candidate { int n, d; std::int64_t score=0; int good=0; };
static void exact(std::istream& in, char* out, std::size_t n) {
  if (!in.read(out, static_cast<std::streamsize>(n))) throw std::runtime_error("short cache frame");
}
static std::uint32_t u32(std::istream& in) {
  unsigned char b[4]; exact(in,reinterpret_cast<char*>(b),4);
  return std::uint32_t(b[0])|(std::uint32_t(b[1])<<8)|(std::uint32_t(b[2])<<16)|(std::uint32_t(b[3])<<24);
}
static bool prime(int p) { for (int q=2;q*q<=p;++q) if (p%q==0) return false; return p>=2; }
int main(int argc,char** argv) {
 try {
  if (argc!=3 && argc!=4) throw std::runtime_error("usage: score_retained_projective_cache CACHE CANDIDATES [MAX_HEIGHT]");
  std::int64_t height=131072;
  if (argc==4) {
   std::size_t used=0;std::string value=argv[3];height=std::stoll(value,&used);
   if (used!=value.size() || height<1 || height>1073741824) throw std::runtime_error("invalid explicit height bound");
  }
  std::ifstream input(argv[2]); std::string magic; std::size_t count;
  if (!(input>>magic>>count) || magic!="R17-CANDIDATES-V1" || count<1 || count>1048576) throw std::runtime_error("invalid candidate header");
  std::vector<Candidate> rows(count);
  for (auto& r:rows) {
   if (!(input>>r.n>>r.d) || r.n==0 || r.n < -height || r.n>height || r.d<1 || r.d>height || std::gcd(r.n,r.d)!=1) throw std::runtime_error("invalid primitive candidate");
  }
  std::string extra; if (input>>extra) throw std::runtime_error("trailing candidate input");
  std::ifstream cache(argv[1],std::ios::binary);char header[8];exact(cache,header,8);
  if (std::string(header,8)!="R17XS001") throw std::runtime_error("invalid cache magic");
  const auto primes=u32(cache);if (primes<1 || primes>8192) throw std::runtime_error("invalid prime count");
  int previous=0;
  for (std::uint32_t j=0;j<primes;++j) {
   int p=static_cast<int>(u32(cache));auto length=u32(cache);
   if (p<=previous || p<5 || p>32749 || !prime(p) || length!=static_cast<std::uint32_t>(p+1)) throw std::runtime_error("invalid prime frame");
   previous=p;std::vector<unsigned char> bytes(std::size_t(length)*8),good(length);std::vector<std::int64_t> units(length);
   exact(cache,reinterpret_cast<char*>(bytes.data()),bytes.size());exact(cache,reinterpret_cast<char*>(good.data()),good.size());
   for (std::size_t i=0;i<length;++i) {
    std::uint64_t v=0;for (int k=0;k<8;++k) v|=std::uint64_t(bytes[8*i+k])<<(8*k);
    units[i]=(v>>63) ? -1-static_cast<std::int64_t>(~v) : static_cast<std::int64_t>(v);
    if (good[i]>1 || (!good[i] && units[i]!=0) || units[i]<-10000000000000LL || units[i]>10000000000000LL) throw std::runtime_error("invalid score symbol");
   }
   std::vector<int> inverse(p);inverse[1]=1;
   for (int i=2;i<p;++i) inverse[i]=static_cast<int>(p-std::int64_t(p/i)*inverse[p%i]%p);
   for (auto& r:rows) {
    int d=r.d%p;int n=r.n%p;if(n<0)n+=p;
    int t=d ? static_cast<int>(std::int64_t(n)*inverse[d]%p) : p;
    r.score+=units[t];r.good+=good[t];
   }
  }
  exact(cache,header,8);
  if (std::string(header,8)!="ENDXSC01" || cache.peek()!=std::char_traits<char>::eof()) throw std::runtime_error("invalid cache end");
  for (std::size_t i=0;i<count;++i) std::cout<<"R "<<i<<" "<<rows[i].score<<" "<<rows[i].good<<"\n";
  std::cout<<"S "<<count<<" "<<primes<<"\n";
 } catch(const std::exception& e) { std::cerr<<e.what()<<"\n";return 1; }
 return 0;
}
