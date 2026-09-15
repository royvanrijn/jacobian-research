// Independent generic-polynomial cubic-plus-rational signature replay for Q80 k0.
#include <array>
#include <algorithm>
#include <cassert>
#include <ctime>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>
#include <sys/resource.h>
using namespace std;
constexpr int p=131,q=p*p*p;
using F=array<int,3>;
const F z{0,0,0},one{1,0,0},theta{0,1,0};
int mod(long long n){n%=p;return int(n<0?n+p:n);}
F plusf(F a,F b){for(int i=0;i<3;i++)a[i]=mod(a[i]+b[i]);return a;}
F minusf(F a,F b){for(int i=0;i<3;i++)a[i]=mod(a[i]-b[i]);return a;}
F scalef(F a,int c){for(int&i:a)i=mod(i*c);return a;}
F times(F a,F b){
 array<long long,5>c{};for(int i=0;i<3;i++)for(int j=0;j<3;j++)c[i+j]+=a[i]*b[j];
 for(int i=4;i>=3;i--){c[i-3]-=3*c[i];c[i-2]-=c[i];}
 return {mod(c[0]),mod(c[1]),mod(c[2])};
}
F powf(F a,int n){F r=one;for(;n;n>>=1,a=times(a,a))if(n&1)r=times(r,a);return r;}
F inv(F a){assert(a!=z);return a==one?one:powf(a,q-2);}
int pack(F a){return a[0]+p*a[1]+p*p*a[2];}
F unpack(int x){return {x%p,x/p%p,x/(p*p)};}
int normf(F a){
 array<F,3>cols{a,times(a,theta),times(a,times(theta,theta))};
 long long n=0;array<int,3>v{0,1,2};
 do{int sign=1;for(int i=0;i<3;i++)for(int j=i+1;j<3;j++)if(v[i]>v[j])sign=-sign;
    n+=sign*cols[0][v[0]]*cols[1][v[1]]*cols[2][v[2]];
 }while(next_permutation(v.begin(),v.end()));
 return mod(n);
}
using Poly=vector<F>;
void trim(Poly&v){while(v.size()&&v.back()==z)v.pop_back();}
Poly sum(Poly a,const Poly&b){a.resize(max(a.size(),b.size()),z);for(size_t i=0;i<b.size();i++)a[i]=plusf(a[i],b[i]);trim(a);return a;}
Poly scaled(Poly a,F n){for(F&v:a)v=times(v,n);trim(a);return a;}
Poly difference(Poly a,const Poly&b){return sum(a,scaled(b,F{p-1,0,0}));}
Poly product(const Poly&a,const Poly&b){
 if(a.empty()||b.empty())return {};Poly r(a.size()+b.size()-1,z);
 for(size_t i=0;i<a.size();i++)for(size_t j=0;j<b.size();j++)r[i+j]=plusf(r[i+j],times(a[i],b[j]));
 trim(r);return r;
}
pair<Poly,Poly> division(Poly a,Poly b){
 trim(a);trim(b);assert(!b.empty());Poly quotient(a.size()>=b.size()?a.size()-b.size()+1:0,z);F bi=inv(b.back());
 while(a.size()>=b.size()){
   size_t k=a.size()-b.size();F n=times(a.back(),bi);quotient[k]=n;
   for(size_t j=0;j<b.size();j++)a[k+j]=minusf(a[k+j],times(n,b[j]));
   trim(a);
 }
 trim(quotient);return {quotient,a};
}
Poly powerp(Poly a,int n,const Poly&f){Poly r{one};while(n){if(n&1)r=division(product(r,a),f).second;a=division(product(a,a),f).second;n>>=1;}return r;}
Poly gcdp(Poly a,Poly b){while(!b.empty()){Poly r=division(a,b).second;a=b;b=r;}return a.empty()?a:scaled(a,inv(a.back()));}
F evaluate(const Poly&a,F t){F r=z;for(int i=int(a.size())-1;i>=0;i--)r=plusf(times(r,t),a[i]);return r;}
vector<F> allroots(F A,F B,int&maximum){
 Poly f{B,A,z,one},x{z,one};
 F delta=plusf(scalef(powf(A,3),4),scalef(powf(B,2),27));
 if(delta==z){assert(A!=z);F n=scalef(times(B,inv(scalef(A,2))),p-3);return {n,n,scalef(n,p-2)};}
 Poly xp=powerp(x,p,f),xp2=division(product(xp,xp),f).second,u=x;
 // Apply the coefficient p-Frobenius and substitute X^p three times.
 for(int j=0;j<3;j++){
   u.resize(3,z);Poly next{powf(u[0],p)};
   next=sum(next,scaled(xp,powf(u[1],p)));
   next=sum(next,scaled(xp2,powf(u[2],p)));u=next;
 }
 if(u!=x)return {};
 for(int k=1;k<=64;k++){
   // Independent splitter uses theta multiples, rather than prime-field constants.
   Poly h=powerp(Poly{scalef(theta,k),one},(q-1)/2,f);h=difference(h,Poly{one});Poly g=gcdp(f,h);
   if(g.size()!=2&&g.size()!=3)continue;
   Poly linear=g.size()==2?g:division(f,g).first;assert(linear.size()==2);
   F r=scalef(times(linear[0],inv(linear[1])),p-1);assert(evaluate(f,r)==z);
   Poly factor=division(f,Poly{scalef(r,p-1),one}).first;assert(factor.size()==3&&factor[2]==one);
   F d=minusf(times(factor[1],factor[1]),scalef(factor[0],4)),s=powf(d,(q+1)/4);
   assert(times(s,s)==d&&s!=z);
   vector<F> roots{r,scalef(plusf(scalef(factor[1],p-1),s),66),scalef(minusf(scalef(factor[1],p-1),s),66)};
   sort(roots.begin(),roots.end(),[](F a,F b){return pack(a)<pack(b);});
   assert(roots[0]!=roots[1]&&roots[1]!=roots[2]);maximum=max(maximum,k);return roots;
 }
 throw runtime_error("independent splitter reached its declared bound");
}
Poly readpoly(){int n;cin>>n;Poly v(n);for(F&a:v){cin>>a[0];a[1]=a[2]=0;}return v;}
int main(int argc,char**argv){
 rlimit cpu{90,95},mem{512UL*1024*1024,512UL*1024*1024};setrlimit(RLIMIT_CPU,&cpu);setrlimit(RLIMIT_AS,&mem);
 clock_t started=clock();int limit=argc>1?stoi(argv[1]):749320;
 Poly A=readpoly(),B=readpoly();int n;cin>>n;assert(n==17);vector<pair<Poly,Poly>>basis;
 for(int j=0;j<n;j++){Poly num=readpoly(),den=readpoly();basis.push_back({num,den});}
 int nc;cin>>nc;vector<int>agreement(1<<17,-1),leg(p);
 for(int j=0;j<nc;j++){int code,t;cin>>code>>t;agreement[code]=t;}assert(nc==110);
 int ns;cin>>ns;assert(ns==16);vector<array<int,3>> signatures(ns);for(auto& sig:signatures){for(int& c:sig)cin>>c;sort(sig.begin(),sig.end());}
 for(int j=1;j<p;j++)leg[j]=powf(F{j,0,0},(p-1)/2)[0];
 F thp=powf(theta,p),th2p=times(thp,thp);
 auto frob=[&](F x){return plusf(F{x[0],0,0},plusf(scalef(thp,x[1]),scalef(th2p,x[2])));};
 int visited=0,split=0,nodal=0,maxtrials=0,candidates=0;bool first=true;
 cout<<"{\"signature_rows\":[";
 for(int id=p;id<q&&visited<limit;id++){
   F t=unpack(id),fp=frob(t);if(id>pack(fp)||id>pack(frob(fp)))continue;
   visited++;F a=evaluate(A,t),b=evaluate(B,t);auto rr=allroots(a,b,maxtrials);if(rr.empty())continue;
   bool smooth=rr[0]!=rr[1];if(smooth)split++;else nodal++;
   vector<int>codes(3,0);
   for(int j=0;j<17;j++){
     F num=evaluate(basis[j].first,t),den=evaluate(basis[j].second,t);
     if(den==z){assert(num!=z);continue;}
     // Quotient evaluation with a field inverse, independently of norm numerator times denominator.
     F value=times(num,inv(den));
     for(int k=0;k<3;k++){
       F c=minusf(value,rr[k]);if(c==z)c=plusf(scalef(times(rr[k],rr[k]),3),a);
       int norm=normf(c);assert(norm);if(leg[norm]==p-1)codes[k]|=1<<j;
     }
   }
   assert((codes[0]^codes[1]^codes[2])==0);
   array<int,3> sig{codes[0],codes[1],codes[2]};sort(sig.begin(),sig.end());if(find(signatures.begin(),signatures.end(),sig)==signatures.end())continue;
   int agree=-1;for(int k=0;k<3;k++)if(codes[k]==0&&agreement[codes[(k+1)%3]]!=-1)agree=codes[(k+1)%3];
   if(agree!=-1)candidates++;
   if(!first)cout<<",";first=false;
   cout<<"{\"t\":"<<id<<",\"smooth\":"<<(smooth?"true":"false")<<",\"roots\":[";
   for(int k=0;k<3;k++)cout<<(k?",":"")<<pack(rr[k]);
   cout<<"],\"norm_codes\":[";for(int k=0;k<3;k++)cout<<(k?",":"")<<codes[k];
   cout<<"],\"agreement_code\":"<<agree<<"}";
 }
 cout<<"],\"status\":\""<<(visited==749320?"COMPLETE":"BOUNDED_PREVIEW")<<"\",\"orbits\":"<<visited<<",\"smooth_split\":"<<split<<",\"nodal\":"<<nodal
     <<",\"max_factor_trials\":"<<maxtrials<<",\"candidates\":"<<candidates<<",\"cpu_seconds\":"<<double(clock()-started)/CLOCKS_PER_SEC<<"}\n";
}
