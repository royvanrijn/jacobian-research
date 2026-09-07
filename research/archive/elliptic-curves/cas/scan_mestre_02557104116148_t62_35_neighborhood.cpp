// Fresh local-score/trace neighborhood scan around exact-rank-16 T=62/35.
//
// The search exhausts reduced positive rationals with denominator 1001..5000
// and 11/7 <= T <= 69/35.  Thus it is disjoint from the preceding complete
// denominator<=1000 rectangle.  Two discovery tails are closed independently:
// Nagao local rank score and squared Frobenius-trace distance from T=62/35.
// Only their fixed union is evaluated at the disjoint held band.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

struct LocalSymbol { bool good=false; int trace=0; std::int64_t score=0; };
struct PrimeTable { int prime=0; std::vector<LocalSymbol> symbols; int anchor_trace=0; };
struct Candidate {
  int numerator=0, denominator=1;
  std::int64_t discovery_score=0, held_score=0;
  std::int64_t discovery_similarity=0, held_similarity=0;
  int discovery_good=0, held_good=0;
  int discovery_matches=0, held_matches=0;
};

static constexpr std::array<int,17> DISCOVERY_PRIMES{
  811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919};
static constexpr std::array<int,15> HELD_PRIMES{
  929,937,941,947,953,967,971,977,983,991,997,1009,1013,1019,1021};
static const std::array<std::string,9> A_COEFFICIENTS{{
  "-27546462334108146267","0","-1530834958134000","0",
  "3106634557536","0","-359251200","0","-34992"}};
static const std::array<std::string,13> B_COEFFICIENTS{{
  "55624396621360883431446459126","0","5508337118494541464141200","0",
  "-13510101631102695979884","0","4296246530145998400","0",
  "-634269830133888","0","38799129600","0","2519424"}};

static int decimal_mod(const std::string& s,int p){
  bool neg=!s.empty()&&s[0]=='-'; int x=0;
  for(std::size_t i=neg?1:0;i<s.size();++i)x=(x*10+s[i]-'0')%p;
  return neg&&x?p-x:x;
}
static int mul(std::int64_t a,std::int64_t b,int p){return int((a*b)%p);}
static int power(int a,int n,int p){int r=1;while(n){if(n&1)r=mul(r,a,p);a=mul(a,a,p);n>>=1;}return r;}
template<std::size_t N> static int eval(const std::array<std::string,N>& c,int x,int p){
  int r=0;for(std::size_t o=0;o<N;++o){std::size_t i=N-1-o;r=(mul(r,x,p)+decimal_mod(c[i],p))%p;}return r;
}
static int trace(int a,int b,int p){
  int s=0;for(int x=0;x<p;++x){int rhs=(mul(mul(x,x,p),x,p)+mul(a,x,p)+b)%p;if(!rhs)continue;int q=power(rhs,(p-1)/2,p);s+=q==1?1:-1;}return -s;
}
static LocalSymbol symbol(int residue,bool infinity,int p){
  int a=infinity?decimal_mod(A_COEFFICIENTS.back(),p):eval(A_COEFFICIENTS,residue,p);
  int b=infinity?decimal_mod(B_COEFFICIENTS.back(),p):eval(B_COEFFICIENTS,residue,p);
  if((4*mul(mul(a,a,p),a,p)+27*mul(b,b,p))%p==0)return {};
  int t=trace(a,b,p), order=p+1-t;
  double score=(double(2-t)/double(order))*std::log(double(p));
  return {true,t,std::int64_t(std::llround(score*1e12))};
}
template<std::size_t N> static std::vector<PrimeTable> tables(const std::array<int,N>& ps){
  std::vector<PrimeTable> out;out.reserve(N);
  for(int p:ps){PrimeTable q;q.prime=p;q.symbols.reserve(p+1);for(int r=0;r<p;++r)q.symbols.push_back(symbol(r,false,p));q.symbols.push_back(symbol(0,true,p));
    int anchor=mul(62%p,power(35%p,p-2,p),p);q.anchor_trace=q.symbols[anchor].trace;out.push_back(std::move(q));}
  return out;
}
static void score_band(Candidate& c,const std::vector<PrimeTable>& ts,bool held){
  std::int64_t rank=0,sim=0;int good=0,matches=0;
  for(const auto& q:ts){int p=q.prime,index=p;if(c.denominator%p){index=mul(c.numerator%p,power(c.denominator%p,p-2,p),p);}const auto& s=q.symbols[index];if(!s.good)continue;
    rank+=s.score;++good;int d=s.trace-q.anchor_trace;sim+=std::int64_t(std::llround(-(double(d)*d/(4.0*p))*1e12));if(!d)++matches;}
  if(held){c.held_score=rank;c.held_similarity=sim;c.held_good=good;c.held_matches=matches;}
  else{c.discovery_score=rank;c.discovery_similarity=sim;c.discovery_good=good;c.discovery_matches=matches;}
}
static bool better_rank(const Candidate&a,const Candidate&b){
  if(a.discovery_score!=b.discovery_score)return a.discovery_score>b.discovery_score;
  if(a.discovery_good!=b.discovery_good)return a.discovery_good>b.discovery_good;
  if(a.discovery_similarity!=b.discovery_similarity)return a.discovery_similarity>b.discovery_similarity;
  if(a.denominator!=b.denominator)return a.denominator<b.denominator;return a.numerator<b.numerator;
}
static bool better_trace(const Candidate&a,const Candidate&b){
  if(a.discovery_similarity!=b.discovery_similarity)return a.discovery_similarity>b.discovery_similarity;
  if(a.discovery_matches!=b.discovery_matches)return a.discovery_matches>b.discovery_matches;
  if(a.discovery_score!=b.discovery_score)return a.discovery_score>b.discovery_score;
  if(a.denominator!=b.denominator)return a.denominator<b.denominator;return a.numerator<b.numerator;
}
struct RankCmp{bool operator()(const Candidate&a,const Candidate&b)const{return better_rank(a,b);}};
struct TraceCmp{bool operator()(const Candidate&a,const Candidate&b)const{return better_trace(a,b);}};
static std::uint64_t key(int a,int b){return(std::uint64_t(std::uint32_t(a))<<32)|std::uint32_t(b);}
static std::unordered_set<std::uint64_t> exclusions(const std::string& path){
  std::ifstream in(path);if(!in)throw std::runtime_error("cannot open exclusions");std::unordered_set<std::uint64_t> s;int a,b;while(in>>a>>b){if(a<=0||b<=0||std::gcd(a,b)!=1)throw std::runtime_error("bad exclusion");s.insert(key(a,b));}if(!in.eof())throw std::runtime_error("bad exclusion stream");return s;
}
static std::uint64_t digest(const std::vector<PrimeTable>& ts){
  std::uint64_t h=1469598103934665603ULL;auto mix=[&](std::uint64_t x){for(int i=0;i<8;++i){h^=(x>>(8*i))&255ULL;h*=1099511628211ULL;}};
  for(const auto&q:ts){mix(q.prime);for(const auto&s:q.symbols){mix(s.good);mix(std::uint64_t(std::int64_t(s.trace)));}mix(std::uint64_t(std::int64_t(q.anchor_trace)));}return h;
}
static std::string decimal(std::int64_t u){bool n=u<0;std::uint64_t a=n?-u:u;std::string t=std::to_string(a%1000000000000ULL);return(n?"-":"")+std::to_string(a/1000000000000ULL)+"."+std::string(12-t.size(),'0')+t;}

int main(int argc,char**argv){
  if(argc!=6){std::cerr<<"usage: scanner DEN_MIN DEN_MAX SCORE_KEEP TRACE_KEEP EXCLUSIONS\n";return 2;}
  int dmin=std::atoi(argv[1]),dmax=std::atoi(argv[2]),ks=std::atoi(argv[3]),kt=std::atoi(argv[4]);
  if(dmin<1||dmax<dmin||dmax>20000||ks<1||kt<1||ks>20000||kt>20000)return 2;
  auto excluded=exclusions(argv[5]);auto discovery=tables(DISCOVERY_PRIMES),held=tables(HELD_PRIMES);
  std::priority_queue<Candidate,std::vector<Candidate>,RankCmp> hs;
  std::priority_queue<Candidate,std::vector<Candidate>,TraceCmp> ht;
  std::uint64_t primitive=0,prior=0,evaluated=0;
  for(int b=dmin;b<=dmax;++b){int lo=(11*b+6)/7,hi=(69*b)/35;for(int a=lo;a<=hi;++a){if(std::gcd(a,b)!=1)continue;++primitive;if(excluded.count(key(a,b))){++prior;continue;}++evaluated;Candidate c;c.numerator=a;c.denominator=b;score_band(c,discovery,false);
      if(int(hs.size())<ks)hs.push(c);else if(better_rank(c,hs.top())){hs.pop();hs.push(c);}if(int(ht.size())<kt)ht.push(c);else if(better_trace(c,ht.top())){ht.pop();ht.push(c);}}}
  std::unordered_map<std::uint64_t,Candidate> unioned;while(!hs.empty()){auto c=hs.top();hs.pop();unioned[key(c.numerator,c.denominator)]=c;}while(!ht.empty()){auto c=ht.top();ht.pop();unioned[key(c.numerator,c.denominator)]=c;}
  std::vector<Candidate> out;out.reserve(unioned.size());for(auto&item:unioned){score_band(item.second,held,true);out.push_back(item.second);}std::sort(out.begin(),out.end(),[](const Candidate&a,const Candidate&b){return a.denominator!=b.denominator?a.denominator<b.denominator:a.numerator<b.numerator;});
  Candidate anchor;anchor.numerator=62;anchor.denominator=35;score_band(anchor,discovery,false);score_band(anchor,held,true);
  std::cout<<"MESTRE_02557104116148_T62_35_NEIGHBORHOOD_V1\nD";for(int p:DISCOVERY_PRIMES)std::cout<<' '<<p;std::cout<<"\nH";for(int p:HELD_PRIMES)std::cout<<' '<<p;
  std::cout<<"\nL "<<digest(discovery)<<' '<<digest(held)<<"\nK 62 35 "<<decimal(anchor.discovery_score)<<' '<<decimal(anchor.held_score)<<' '<<anchor.discovery_good<<' '<<anchor.held_good<<' '<<decimal(anchor.discovery_similarity)<<' '<<decimal(anchor.held_similarity)<<' '<<anchor.discovery_matches<<' '<<anchor.held_matches<<'\n';
  for(const auto&c:out)std::cout<<"C "<<c.numerator<<' '<<c.denominator<<' '<<decimal(c.discovery_score)<<' '<<decimal(c.held_score)<<' '<<c.discovery_good<<' '<<c.held_good<<' '<<decimal(c.discovery_similarity)<<' '<<decimal(c.held_similarity)<<' '<<c.discovery_matches<<' '<<c.held_matches<<'\n';
  std::cout<<"S "<<dmin<<' '<<dmax<<' '<<ks<<' '<<kt<<' '<<primitive<<' '<<prior<<' '<<evaluated<<' '<<out.size()<<' '<<excluded.size()<<'\n';
}
