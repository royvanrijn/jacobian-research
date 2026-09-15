// Complete cubic-plus-rational signature gate for the Q80 genus-one k0 chart.
#include <array>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <ctime>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/resource.h>
using namespace std;
constexpr int p=131,Q=p*p*p;
struct F {int a,b,c;};
inline bool eq(F x,F y){return x.a==y.a&&x.b==y.b&&x.c==y.c;}
inline bool zero(F x){return !(x.a|x.b|x.c);}
inline F add(F x,F y){return {(x.a+y.a)%p,(x.b+y.b)%p,(x.c+y.c)%p};}
inline F neg(F x){return {(p-x.a)%p,(p-x.b)%p,(p-x.c)%p};}
inline F sub(F x,F y){return {(p+x.a-y.a)%p,(p+x.b-y.b)%p,(p+x.c-y.c)%p};}
inline F scale(F x,int n){return {x.a*n%p,x.b*n%p,x.c*n%p};}
inline F mul(F x,F y){
 int z=x.b*y.c+x.c*y.b,w=x.c*y.c;
 return {(8*p*p+x.a*y.a-3*z)%p,(8*p*p+x.a*y.b+x.b*y.a-z-3*w)%p,
         (8*p*p+x.a*y.c+x.b*y.b+x.c*y.a-w)%p};
}
F power(F a,int n){F b{1,0,0};while(n){if(n&1)b=mul(b,a);a=mul(a,a);n>>=1;}return b;}
int packed(F x){return x.a+p*x.b+p*p*x.c;}
F unpack(int n){return {n%p,n/p%p,n/(p*p)};}
F fr1,fr2;
F frob(F x){return add(F{x.a,0,0},add(scale(fr1,x.b),scale(fr2,x.c)));}
int norm(F x){
 auto y=frob(x),z=frob(y);auto n=mul(mul(x,y),z);assert(!n.b&&!n.c);return n.a;
}
array<int,p> invp,legendre;
F inverse(F x){assert(!zero(x));F y=frob(x),z=frob(y);F yz=mul(y,z);F n=mul(x,yz);assert(!n.b&&!n.c&&n.a);return scale(yz,invp[n.a]);}
using Poly=array<F,3>;
Poly one(){return {{{1,0,0},{0,0,0},{0,0,0}}};}
Poly pmul(const Poly& u,const Poly& v,F A,F B){
 F o0=mul(u[0],v[0]),o1=add(mul(u[0],v[1]),mul(u[1],v[0]));
 F o2=add(add(mul(u[0],v[2]),mul(u[1],v[1])),mul(u[2],v[0]));
 F o3=add(mul(u[1],v[2]),mul(u[2],v[1])),o4=mul(u[2],v[2]);
 return {sub(o0,mul(B,o3)),sub(sub(o1,mul(A,o3)),mul(B,o4)),sub(o2,mul(A,o4))};
}
Poly ppow(Poly a,int n,F A,F B){Poly b=one();while(n){if(n&1)b=pmul(b,a,A,B);a=pmul(a,a,A,B);n>>=1;}return b;}
bool equalp(const Poly& a,const Poly& b){return eq(a[0],b[0])&&eq(a[1],b[1])&&eq(a[2],b[2]);}
vector<F> roots(F A,F B,int& trials){
 F disc=add(scale(power(A,3),4),scale(power(B,2),27));
 if(zero(disc)){assert(!zero(A));F node=scale(mul(B,inverse(scale(A,2))),p-3);return {node,node,scale(node,p-2)};}
 Poly x{{{0,0,0},{1,0,0},{0,0,0}}};
 // Direct Q-power; the independent implementation will iterate the p-Frobenius.
 if(!equalp(ppow(x,Q,A,B),x))return {};
 for(int k=0;k<64;k++){
   Poly base=x;base[0]=unpack(k);Poly h=ppow(base,(Q-1)/2,A,B);h[0]=sub(h[0],F{1,0,0});F r;
   if(zero(h[2])){
     if(zero(h[1]))continue;
     r=neg(mul(h[0],inverse(h[1])));
   }else{
     F h22=mul(h[2],h[2]);
     F rx=add(sub(mul(h[1],h[1]),mul(h[0],h[2])),mul(A,h22));
     F r0=add(mul(h[1],h[0]),mul(B,h22));
     if(zero(rx)){assert(zero(r0));r=mul(h[1],inverse(h[2]));}
     else r=neg(mul(r0,inverse(rx)));
   }
   assert(zero(add(add(power(r,3),mul(A,r)),B)));
   F dq=neg(add(scale(A,4),scale(mul(r,r),3)));
   F s=power(dq,(Q+1)/4);assert(eq(mul(s,s),dq)&&!zero(s));
   vector<F> out{r,scale(add(neg(r),s),66),scale(sub(neg(r),s),66)};
   sort(out.begin(),out.end(),[](F a,F b){return packed(a)<packed(b);});
   assert(!eq(out[0],out[1])&&!eq(out[1],out[2]));trials=max(trials,k+1);return out;
 }
 throw runtime_error("split cubic exceeded the declared 64-trial factorization bound");
}
F eval(const vector<int>& c,F x){F y{0,0,0};for(int k=(int)c.size()-1;k>=0;k--)y=add(mul(y,x),F{c[k],0,0});return y;}
vector<int> readpoly(){int n;cin>>n;vector<int>v(n);for(int&i:v)cin>>i;return v;}
struct Row {int t;vector<int>roots,codes;bool smooth;int agreement=-1;};
void printrow(const Row&r){
 cout<<"{\"t\":"<<r.t<<",\"smooth\":"<<(r.smooth?"true":"false")<<",\"roots\":[";
 for(size_t i=0;i<r.roots.size();i++)cout<<(i?",":"")<<r.roots[i];
 cout<<"],\"norm_codes\":[";for(size_t i=0;i<r.codes.size();i++)cout<<(i?",":"")<<r.codes[i];
 cout<<"],\"agreement_code\":"<<r.agreement<<"}";
}
int main(int argc,char**argv){
 rlimit cpu{90,95},mem{512UL*1024*1024,512UL*1024*1024};setrlimit(RLIMIT_CPU,&cpu);setrlimit(RLIMIT_AS,&mem);
 clock_t start=clock();int limit=argc>1?stoi(argv[1]):749320;
 fr1=power(F{0,1,0},p);fr2=mul(fr1,fr1);
 for(int i=1;i<p;i++){invp[i]=power(F{i,0,0},p-2).a;legendre[i]=power(F{i,0,0},(p-1)/2).a;}
 auto A=readpoly(),B=readpoly();int n;cin>>n;assert(n==17);
 vector<pair<vector<int>,vector<int>>> sections;
 for(int j=0;j<n;j++){auto num=readpoly(),den=readpoly();sections.push_back({num,den});}
 int nc;cin>>nc;vector<int>agree(1<<17,-1);for(int j=0;j<nc;j++){int c,t;cin>>c>>t;agree[c]=t;}assert(nc==110);
 int ns;cin>>ns;assert(ns==16);vector<array<int,3>> signatures(ns);for(auto& sig:signatures){for(int& c:sig)cin>>c;sort(sig.begin(),sig.end());}
 int visited=0,split=0,nodal=0,maxtrials=0;vector<Row> zeros,candidates;
 for(int id=p;id<Q&&visited<limit;id++){
   F t=unpack(id),ft=frob(t);if(id>packed(ft)||id>packed(frob(ft)))continue;
   visited++;F a=eval(A,t),b=eval(B,t);auto rr=roots(a,b,maxtrials);if(rr.empty())continue;
   bool smooth=!zero(add(scale(power(a,3),4),scale(power(b,2),27)));if(smooth)split++;else nodal++;
   vector<int>codes(3,0);
   for(int j=0;j<17;j++){
     F num=eval(sections[j].first,t),den=eval(sections[j].second,t);
     if(zero(den)){assert(!zero(num));continue;}
     int nd=norm(den);assert(nd);
     for(int k=0;k<3;k++){
       F value=sub(num,mul(rr[k],den));int v;
       if(zero(value))v=norm(add(scale(mul(rr[k],rr[k]),3),a));
       else v=norm(value)*nd%p;
       assert(v);if(legendre[v]==p-1)codes[k]|=1<<j;
     }
   }
   assert((codes[0]^codes[1]^codes[2])==0);
   array<int,3> sig{codes[0],codes[1],codes[2]};sort(sig.begin(),sig.end());if(find(signatures.begin(),signatures.end(),sig)==signatures.end())continue;
   Row row{id,{},codes,smooth};for(F e:rr)row.roots.push_back(packed(e));
   for(int j=0;j<3;j++)if(codes[j]==0&&agree[codes[(j+1)%3]]!=-1)row.agreement=codes[(j+1)%3];
   zeros.push_back(row);if(row.agreement!=-1)candidates.push_back(row);
 }
 cout<<"{\"status\":\""<<(visited==749320?"COMPLETE":"BOUNDED_PREVIEW")<<"\",\"orbits\":"<<visited<<",\"smooth_split\":"<<split<<",\"nodal\":"<<nodal<<",\"max_factor_trials\":"<<maxtrials<<",\"signature_rows\":[";
 for(size_t i=0;i<zeros.size();i++){if(i)cout<<",";printrow(zeros[i]);}
 cout<<"],\"candidate_rows\":[";for(size_t i=0;i<candidates.size();i++){if(i)cout<<",";printrow(candidates[i]);}
 cout<<"],\"cpu_seconds\":"<<double(clock()-start)/CLOCKS_PER_SEC<<"}\n";
}
