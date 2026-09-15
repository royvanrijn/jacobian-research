// Independent monomial arithmetic modulo z^4-2*z^2+2; Tonelli-Shanks roots.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <ctime>
#include <iostream>
#include <vector>
#include <sys/resource.h>
using namespace std;
constexpr int p=131,p2=p*p,q=p2*p2;
using F=array<int,4>;
const F zero{0,0,0,0},one{1,0,0,0},z{0,1,0,0};
int mod(long long x){x%=p;return int(x<0?x+p:x);}
F add(F x,F y){for(int i=0;i<4;i++)x[i]=mod(x[i]+y[i]);return x;}
F sub(F x,F y){for(int i=0;i<4;i++)x[i]=mod(x[i]-y[i]);return x;}
F scale(F x,int n){for(int& a:x)a=mod(a*n);return x;}
F mul(F x,F y){
 array<long long,7> c{};for(int i=0;i<4;i++)for(int j=0;j<4;j++)c[i+j]+=x[i]*y[j];
 for(int i=6;i>=4;i--){c[i-2]+=2*c[i];c[i-4]-=2*c[i];}
 return {mod(c[0]),mod(c[1]),mod(c[2]),mod(c[3])};
}
F power(F x,int n){F y=one;while(n){if(n&1)y=mul(y,x);x=mul(x,x);n>>=1;}return y;}
// Convert the externally specified tower encoding to the independent monomial basis.
F unpack(int n){int a=n%p,b=n/p%p,c=n/p2%p,d=n/(p*p2);return {mod(a-b),mod(c-d),b,d};}
int pack(F x){return mod(x[0]+x[2])+p*x[2]+p2*mod(x[1]+x[3])+p*p2*x[3];}
array<F,4> frob_columns;
F frob(F x){F y=zero;for(int i=0;i<4;i++)y=add(y,scale(frob_columns[i],x[i]));return y;}
array<int,p> invp,leg;
int norm(F x){F a=frob(x),b=frob(a),c=frob(b),n=mul(mul(x,a),mul(b,c));assert(n[1]==0&&n[2]==0&&n[3]==0);return n[0];}
F inverse(F x){F a=frob(x),b=frob(a),c=frob(b),prod=mul(mul(a,b),c),n=mul(x,prod);assert(n[0]&&n[1]==0&&n[2]==0&&n[3]==0);return scale(prod,invp[n[0]]);}
F ts_constant;
bool sqrtf(F x,F& y){
 if(x==zero){y=zero;return true;}
 int m=4;F c=ts_constant,t=power(x,(q-1)/16);y=power(x,((q-1)/16+1)/2);
 while(t!=one){
   F u=t;int i=0;while(u!=one&&i<m){u=mul(u,u);i++;}
   if(i==m)return false;
   F b=power(c,1<<(m-i-1));y=mul(y,b);c=mul(b,b);t=mul(t,c);m=i;
 }
 assert(mul(y,y)==x);return true;
}
vector<F> roots(F A,F B){
 F s,third=scale(A,44),half=scale(B,66),delta=add(mul(half,half),mul(mul(third,third),third));
 assert(delta!=zero);if(!sqrtf(delta,s))return {};
 // Use the opposite resolvent from the tower producer.
 F U=sub(scale(half,p-1),s);if(U==zero)U=add(scale(half,p-1),s);
 assert(U!=zero);F a=power(U,65444427);if(mul(mul(a,a),a)!=U)return {};
 F b=scale(mul(third,inverse(a)),p-1),w=unpack(65+19*p),w2=mul(w,w);
 vector<F> r{add(a,b),add(mul(w,a),mul(w2,b)),add(mul(w2,a),mul(w,b))};
 for(F e:r)assert(add(add(mul(mul(e,e),e),mul(A,e)),B)==zero);
 sort(r.begin(),r.end(),[](F a,F b){return pack(a)<pack(b);});assert(r[0]!=r[1]&&r[1]!=r[2]);return r;
}
vector<int> readpoly(){int n;cin>>n;assert(cin&&n>=0&&n<1000);vector<int> c(n);for(int& a:c){cin>>a;assert(cin&&a>=0&&a<p);}return c;}
F eval(const vector<int>& a,F t){F y=zero;for(int i=int(a.size())-1;i>=0;i--){y=mul(y,t);y[0]=mod(y[0]+a[i]);}return y;}
uint64_t mix(uint64_t x){x^=x>>30;x*=0xbf58476d1ce4e5b9ULL;x^=x>>27;x*=0x94d049bb133111ebULL;return x^(x>>31);}
int main(int argc,char**argv){
 rlimit cpu{90,95},mem{512UL*1024*1024,512UL*1024*1024};assert(!setrlimit(RLIMIT_CPU,&cpu));assert(!setrlimit(RLIMIT_AS,&mem));clock_t start=clock();
 assert(argc==3);int begin=stoi(argv[1]),end=stoi(argv[2]);assert(0<=begin&&begin<end&&end<=4290);
 for(int i=0;i<4;i++){F e=zero;e[i]=1;frob_columns[i]=power(e,p);}
 for(int i=1;i<p;i++){invp[i]=power(F{i,0,0,0},p-2)[0];leg[i]=power(F{i,0,0,0},65)[0];}
 assert(power(z,(q-1)/2)==scale(one,p-1));ts_constant=power(z,(q-1)/16);
 vector<int> pairs;for(int c=1;c<p2;c++){
   F x=unpack(p2*c);int least=c;F y=x;
   for(int k=1;k<=4;k++){y=frob(y);assert(pack(y)%p2==0);least=min(least,pack(y)/p2);}
   assert(y==x);if(c==least)pairs.push_back(c);
 }
 assert(pairs.size()==4290);
 auto A=readpoly(),B=readpoly();int n;cin>>n;assert(n==17);vector<pair<vector<int>,vector<int>>> basis;
 for(int i=0;i<n;i++){auto num=readpoly(),den=readpoly();basis.push_back({num,den});}
 string extra;assert(!(cin>>extra));
 array<uint64_t,18> hist{};uint64_t visited=0,split=0,hash=0;vector<array<int,4>> candidates;
 for(int k=begin;k<end;k++)for(int a=0;a<p2;a++){
   int id=a+p2*pairs[k];F t=unpack(id);visited++;F aa=eval(A,t),bb=eval(B,t);auto rr=roots(aa,bb);if(rr.empty())continue;split++;
   int j;for(j=0;j<17;j++){
     F num=eval(basis[j].first,t),den=eval(basis[j].second,t);if(den==zero){assert(num!=zero);continue;}
     F value=mul(num,inverse(den));bool bad=false;
     for(int i=0;i<3;i++){
       F c=sub(value,rr[i]);if(c==zero)c=add(scale(mul(rr[i],rr[i]),3),aa);
       int nn=norm(c);assert(nn);if(leg[nn]==p-1){bad=true;break;}
     }
     if(bad)break;
   }
   hist[j]++;hash+=mix(uint64_t(id)*32+j);if(j==17)candidates.push_back({id,pack(rr[0]),pack(rr[1]),pack(rr[2])});
 }
 cout<<"{\"status\":\"COMPLETE_RANGE\",\"pair_begin\":"<<begin<<",\"pair_end\":"<<end<<",\"orbits\":"<<visited<<",\"smooth_split\":"<<split<<",\"first_failed_character\":[";
 for(int i=0;i<18;i++)cout<<(i?",":"")<<hist[i];cout<<"],\"split_character_checksum\":\""<<hash<<"\",\"candidate_rows\":[";
 for(size_t i=0;i<candidates.size();i++){if(i)cout<<",";cout<<"[";for(int j=0;j<4;j++)cout<<(j?",":"")<<candidates[i][j];cout<<"]";}
 cout<<"],\"cpu_seconds\":"<<double(clock()-start)/CLOCKS_PER_SEC<<"}\n";
}
