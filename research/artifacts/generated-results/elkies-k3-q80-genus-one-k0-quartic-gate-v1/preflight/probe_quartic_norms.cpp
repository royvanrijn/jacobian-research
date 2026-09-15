// Bounded Q80 quartic norm gate; F_131^4 as F_131[i,j], i^2=-1,j^2=1+i.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <ctime>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <sys/resource.h>
using namespace std;
constexpr int p=131,p2=p*p,q=p2*p2;
struct K{int a,b;}; struct F{K a,b;};
bool eq(K x,K y){return x.a==y.a&&x.b==y.b;}
bool eq(F x,F y){return eq(x.a,y.a)&&eq(x.b,y.b);}
bool zero(K x){return !(x.a|x.b);} bool zero(F x){return zero(x.a)&&zero(x.b);}
K add(K x,K y){return {(x.a+y.a)%p,(x.b+y.b)%p};}
K sub(K x,K y){return {(p+x.a-y.a)%p,(p+x.b-y.b)%p};}
K neg(K x){return {(p-x.a)%p,(p-x.b)%p};}
K scale(K x,int c){return {x.a*c%p,x.b*c%p};}
K mul(K x,K y){return {(p*p+x.a*y.a-x.b*y.b)%p,(x.a*y.b+x.b*y.a)%p};}
K sq(K x){return {(p*p+x.a*x.a-x.b*x.b)%p,2*x.a*x.b%p};}
K u(K x){return {(p+x.a-x.b)%p,(x.a+x.b)%p};}
int pack(K x){return x.a+p*x.b;} K unpack2(int x){return {x%p,x/p};}
int pack(F x){return pack(x.a)+p2*pack(x.b);}
F unpack4(int x){return {unpack2(x%p2),unpack2(x/p2)};}
F add(F x,F y){return {add(x.a,y.a),add(x.b,y.b)};}
F sub(F x,F y){return {sub(x.a,y.a),sub(x.b,y.b)};}
F neg(F x){return {neg(x.a),neg(x.b)};}
F scale(F x,int c){return {scale(x.a,c),scale(x.b,c)};}
F mul(F x,F y){K ac=mul(x.a,y.a),bd=mul(x.b,y.b);return {add(ac,u(bd)),sub(sub(mul(add(x.a,x.b),add(y.a,y.b)),ac),bd)};}
F sq(F x){return {add(sq(x.a),u(sq(x.b))),scale(mul(x.a,x.b),2)};}
F power(F x,int n){F y{{1,0},{0,0}};while(n){if(n&1)y=mul(y,x);x=sq(x);n>>=1;}return y;}
array<int,p> invp,leg;
array<int,p2> sqrt2;
K inverse(K x){int n=(x.a*x.a+x.b*x.b)%p;assert(n);return scale(K{x.a,(p-x.b)%p},invp[n]);}
F inverse(F x){K v=inverse(sub(sq(x.a),u(sq(x.b))));return {mul(x.a,v),neg(mul(x.b,v))};}
int norm(F x){K n=sub(sq(x.a),u(sq(x.b)));return (n.a*n.a+n.b*n.b)%p;}
int character(F x){int n=norm(x);assert(n);return leg[n];}
bool root(K x,K& y){int n=sqrt2[pack(x)];if(n<0)return false;y=unpack2(n);return true;}
bool root(F x,F& y){
 if(zero(x.b)){
   K r;if(root(x.a,r)){y={r,{0,0}};return true;}
   if(!root(mul(x.a,K{66,65}),r))return false;
   y={{0,0},r};return true;
 }
 K t;if(!root(sub(sq(x.a),u(sq(x.b))),t))return false;
 K r;if(!root(scale(add(x.a,t),66),r)){
   bool ok=root(scale(sub(x.a,t),66),r);assert(ok);
 }
 assert(!zero(r));y={r,scale(mul(x.b,inverse(r)),66)};assert(eq(sq(y),x));return true;
}
F eval(const vector<int>& c,F t){F y{{0,0},{0,0}};for(int i=int(c.size())-1;i>=0;i--){y=mul(y,t);y.a.a=(y.a.a+c[i])%p;}return y;}
vector<int> readpoly(){int n;cin>>n;assert(cin&&n>=0&&n<1000);vector<int> v(n);for(int& x:v){cin>>x;assert(cin&&x>=0&&x<p);}return v;}
vector<F> roots(F A,F B){
 F halfB=scale(B,66),thirdA=scale(A,44),s;
 F delta=add(sq(halfB),mul(sq(thirdA),thirdA));
 assert(!zero(delta));if(!root(delta,s))return {};
 F U=add(neg(halfB),s);if(zero(U))U=sub(neg(halfB),s);
 assert(!zero(U));F a=power(U,65444427);if(!eq(mul(sq(a),a),U))return {};
 F b=neg(mul(thirdA,inverse(a))),z{{65,19},{0,0}},zz=sq(z);
 vector<F> r{add(a,b),add(mul(z,a),mul(zz,b)),add(mul(zz,a),mul(z,b))};
 for(F e:r)assert(zero(add(add(mul(sq(e),e),mul(A,e)),B)));
 assert(!eq(r[0],r[1])&&!eq(r[0],r[2])&&!eq(r[1],r[2]));
 sort(r.begin(),r.end(),[](F x,F y){return pack(x)<pack(y);});return r;
}
K L(K x){return {117*(x.a+x.b)%p,117*(p+x.a-x.b)%p};}
vector<K> canonical_pairs(){vector<K> v;for(int n=1;n<p2;n++){K b=unpack2(n),c=L(b);assert(eq(L(c),neg(b)));if(n<pack(c)&&n<pack(neg(b))&&n<pack(neg(c)))v.push_back(b);}assert(v.size()==4290);return v;}
uint64_t mix(uint64_t x){x^=x>>30;x*=0xbf58476d1ce4e5b9ULL;x^=x>>27;x*=0x94d049bb133111ebULL;return x^(x>>31);}
int main(int argc,char**argv){
 rlimit cpu{90,95},mem{512UL*1024*1024,512UL*1024*1024};assert(!setrlimit(RLIMIT_CPU,&cpu));assert(!setrlimit(RLIMIT_AS,&mem));clock_t start=clock();
 assert(argc==3);int begin=stoi(argv[1]),end=stoi(argv[2]);assert(0<=begin&&begin<end&&end<=4290);
 sqrt2.fill(-1);for(int i=0;i<p2;i++)sqrt2[pack(sq(unpack2(i)))]=i;
 for(int i=1;i<p;i++){invp[i]=power(F{{i,0},{0,0}},p-2).a.a;leg[i]=power(F{{i,0},{0,0}},65).a.a;}
 assert(sqrt2[pack(K{1,1})]<0);assert(eq(power(F{{0,0},{1,0}},p),F{{0,0},{117,117}}));
 auto A=readpoly(),B=readpoly();int n;cin>>n;assert(n==17);vector<pair<vector<int>,vector<int>>> sections;
 for(int i=0;i<n;i++){auto num=readpoly(),den=readpoly();sections.push_back({num,den});}
 auto pairs=canonical_pairs();array<uint64_t,18> hist{};uint64_t visited=0,split=0,hash=0;vector<array<int,4>> candidates;
 for(int z=begin;z<end;z++)for(int a=0;a<p2;a++){
   F t{unpack2(a),pairs[z]};visited++;F aa=eval(A,t),bb=eval(B,t);auto rr=roots(aa,bb);if(rr.empty())continue;split++;
   int j;
   for(j=0;j<17;j++){
     F num=eval(sections[j].first,t),den=eval(sections[j].second,t);
     if(zero(den)){assert(!zero(num));continue;}
     int nd=character(den);bool bad=false;
     for(int k=0;k<3;k++){
       F val=sub(num,mul(rr[k],den));int c=zero(val)?character(add(scale(sq(rr[k]),3),aa)):character(val)*nd%p;
       if(c==p-1){bad=true;break;}
     }
     if(bad)break;
   }
   hist[j]++;hash+=mix(uint64_t(pack(t))*32+j);
   if(j==17)candidates.push_back({pack(t),pack(rr[0]),pack(rr[1]),pack(rr[2])});
 }
 cout<<"{\"status\":\"COMPLETE_RANGE\",\"pair_begin\":"<<begin<<",\"pair_end\":"<<end<<",\"orbits\":"<<visited<<",\"smooth_split\":"<<split<<",\"first_failed_character\":[";
 for(int i=0;i<18;i++)cout<<(i?",":"")<<hist[i];
 cout<<"],\"split_character_checksum\":\""<<hash<<"\",\"candidate_rows\":[";
 for(size_t i=0;i<candidates.size();i++){if(i)cout<<",";cout<<"[";for(int j=0;j<4;j++)cout<<(j?",":"")<<candidates[i][j];cout<<"]";}
 cout<<"],\"cpu_seconds\":"<<double(clock()-start)/CLOCKS_PER_SEC<<"}\n";
}
