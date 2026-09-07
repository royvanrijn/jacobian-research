// Exact fibre character sums at one pinned prime and its quadratic extension.
#include <array>
#include <iostream>
#include <vector>
constexpr int p=131;
int mod(int x) { x%=p; return x<0?x+p:x; }
struct F { int a,b; };
int D;
F add(F x,F y) {return {mod(x.a+y.a),mod(x.b+y.b)};}
F mul(F x,F y) {return {mod(x.a*y.a+D*x.b*y.b),mod(x.a*y.b+x.b*y.a)};}
F eval(const std::vector<int>& c,F t) {
  F v{0,0}; for(auto i=c.rbegin();i!=c.rend();++i) v=add(mul(v,t),{*i,0}); return v;
}
int main() {
  std::vector<int>A(9),B(13),q(3); for(auto *v:{&A,&B,&q}) for(int& c:*v) std::cin>>c;
  if(!std::cin) return 2;
  std::array<int,p> chi{};chi.fill(-1);chi[0]=0;for(int i=1;i<p;++i)chi[i*i%p]=1;
  D=2;while(chi[D]!=-1)++D;
  std::cout<<"D "<<D<<'\n';
  for(int n:{1,2}) {
    std::vector<F> xs,cubes;
    for(int b=0;b<(n==1?1:p);++b)for(int a=0;a<p;++a) {
      F x{a,b};xs.push_back(x);cubes.push_back(mul(mul(x,x),x));
    }
    auto character=[&](F x){return n==1?chi[x.a]:chi[mod(x.a*x.a-D*x.b*x.b)];};
    auto count=[&](F a,F b){
      int sum=0;
      for(size_t i=0;i<xs.size();++i) {
        // Keep the inner loop separate from polynomial evaluation.
        F z=add(add(cubes[i],mul(a,xs[i])),b);sum+=character(z);
      }
      return -sum;
    };
    for(int b=0;b<(n==1?1:(p+1)/2);++b)for(int a=0;a<p;++a) {
      F t{a,b};int trace=count(eval(A,t),eval(B,t));int sign=character(eval(q,t));
      std::cout<<n<<' '<<a<<' '<<b<<' '<<(b?2:1)<<' '<<trace<<' '<<sign*trace<<'\n';
    }
    int trace=count({A.back(),0},{B.back(),0});int sign=character({q.back(),0});
    std::cout<<n<<" -1 0 1 "<<trace<<' '<<sign*trace<<'\n';
  }
}
