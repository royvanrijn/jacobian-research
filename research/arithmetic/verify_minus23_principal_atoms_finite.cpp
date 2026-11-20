#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>

// Complete finite existence check only; the infinite onset is a separate proof.
// Every membership bit is set by a fully bounded primitive integer Q-coordinate.
uint64_t isqrt_exact(uint64_t n){
    uint64_t r=std::sqrt((long double)n);
    while((r+1)*(r+1)<=n)++r;
    while(r*r>n)--r;
    return r;
}
int64_t floorhalf(int64_t n){return n>=0?n/2:-((-n+1)/2);}
uint64_t norm(int64_t x,int64_t y){return x*x+x*y+6*y*y;}
struct Bits{
    std::vector<uint64_t> a;
    explicit Bits(uint64_t n):a((n+64)/64,0){}
    void set(uint64_t n){a[n>>6]|=uint64_t(1)<<(n&63);}
    bool get(uint64_t n)const{return(a[n>>6]>>(n&63))&1;}
};
void write32(std::ostream&f,uint32_t a){for(int i=0;i<4;i++)f.put(char((a>>(8*i))&255));}

int main(int argc,char**argv){
    if(argc!=4){std::cerr<<"usage: checker inclusive_bound output_prefix anchor_limit\n";return 2;}
    uint64_t bound=std::stoull(argv[1]);std::string prefix=argv[2];
    uint64_t anchor_limit=std::stoull(argv[3]);
    if(bound<1000||bound>=(uint64_t(1)<<32)||anchor_limit>bound)return 2;
    auto start=std::chrono::steady_clock::now();
    Bits primitive(bound);primitive.set(1);
    uint32_t ymax=isqrt_exact(4*bound/23);
    std::vector<uint32_t> spf(ymax+1);
    std::iota(spf.begin(),spf.end(),0);
    for(uint32_t p=2;p*p<=ymax;p++)if(spf[p]==p)
        for(uint32_t j=p*p;j<=ymax;j+=p)if(spf[j]==j)spf[j]=p;
    std::vector<uint8_t> coprime(2*isqrt_exact(bound)+3);
    uint64_t row_candidates=0,primitive_positive_y_points=0;
    for(uint32_t y=1;y<=ymax;y++){
        uint64_t s=isqrt_exact(4*bound-23*uint64_t(y)*y);
        int64_t lo=-int64_t((s+y)/2),hi=floorhalf(int64_t(s)-y);
        uint64_t len=uint64_t(hi-lo+1);row_candidates+=len;
        if(len>coprime.size())throw std::runtime_error("row bound");
        std::fill(coprime.begin(),coprime.begin()+len,1);
        uint32_t rest=y;
        while(rest>1){uint32_t p=spf[rest];
            int64_t first=((-lo)%p+p)%p;
            for(uint64_t j=first;j<len;j+=p)coprime[j]=0;
            while(rest%p==0)rest/=p;
        }
        for(uint64_t j=0;j<len;j++)if(coprime[j]){
            int64_t x=lo+int64_t(j);uint64_t n=norm(x,y);
            if(n==0||n>bound)throw std::runtime_error("ellipse endpoint");
            primitive.set(n);primitive_positive_y_points++;
        }
        if(y%2000==0)std::cerr<<"ellipse y="<<y<<"/"<<ymax<<"\n";
    }
    // Independent square enumeration challenges all endpoint and row-gcd code.
    const uint64_t small=std::min<uint64_t>(bound,100000);
    Bits square(small);uint64_t ycap=isqrt_exact(4*small/23)+1;
    int64_t xcap=isqrt_exact(small)+ycap+1;
    for(int64_t x=-xcap;x<=xcap;x++)for(int64_t y=-int64_t(ycap);y<=int64_t(ycap);y++){
        uint64_t n=norm(x,y);
        if(n>0&&n<=small&&std::gcd(x,y)==1)square.set(n);
    }
    for(uint64_t n=0;n<=small;n++)if(square.get(n)!=primitive.get(n))throw std::runtime_error("independent square flags");
    if(primitive.get(4)||primitive.get(2)||primitive.get(3)||!primitive.get(1)||!primitive.get(6))throw std::runtime_error("ordinary/primitive controls");
    std::vector<uint32_t> anchors;
    for(uint32_t a=1;a<=anchor_limit;a++)if(a%23&&primitive.get(a))anchors.push_back(a);
    uint64_t prime_cap=bound/5;
    Bits composite(prime_cap);
    composite.set(0);composite.set(1);
    for(uint64_t p=2;p*p<=prime_cap;p++)if(!composite.get(p))
        for(uint64_t j=p*p;j<=prime_cap;j+=p)composite.set(j);
    std::ofstream witness(prefix+"_assignments.bin",std::ios::binary);
    uint64_t all_prime_count=0,principal_prime_count=0,covered[3]={0,0,0},max_trials[3]={0,0,0};
    const int cs[3]={5,10,15};
    std::vector<std::pair<uint64_t,int>> failures;
    for(uint64_t p=2;p<=prime_cap;p++)if(!composite.get(p)){
        all_prime_count++;
        if(p==23||!primitive.get(p))continue;
        principal_prime_count++;
        for(int i=0;i<3;i++){
            uint64_t n=cs[i]*p;if(n>bound)continue;
            uint64_t trials=0;bool found=false;
            for(uint32_t a:anchors){if(a>=n)break;trials++;
                if(primitive.get(n-a)&&((n-a)%23)){
                    covered[i]++;max_trials[i]=std::max(max_trials[i],trials);
                    write32(witness,p);write32(witness,cs[i]);write32(witness,a);
                    found=true;break;
                }
            }
            if(!found)failures.push_back({p,cs[i]});
        }
    }
    witness.close();
    uint64_t flag_count=0;for(uint64_t a:primitive.a)flag_count+=__builtin_popcountll(a);
    std::ofstream flags(prefix+"_primitive_flags_le64.bin",std::ios::binary);
    for(uint64_t a:primitive.a)for(int i=0;i<8;i++)flags.put(char((a>>(8*i))&255));
    flags.close();
    double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::ofstream out(prefix+"_counts.json");
    out<<"{\n\"inclusive_norm_bound\":"<<bound<<",\n\"prime_cap\":"<<prime_cap
       <<",\n\"all_primes_through_cap\":"<<all_prime_count<<",\n\"principal_split_primes_through_cap\":"<<principal_prime_count
       <<",\n\"primitive_norm_value_count\":"<<flag_count<<",\n\"positive_y_ellipse_candidates\":"<<row_candidates
       <<",\n\"primitive_positive_y_points\":"<<primitive_positive_y_points
       <<",\n\"anchor_limit\":"<<anchor_limit<<",\n\"anchor_count\":"<<anchors.size()
       <<",\n\"coefficient_domain\":[5,10,15],\n\"covered_by_coefficient\":["<<covered[0]<<","<<covered[1]<<","<<covered[2]
       <<"],\n\"maximum_anchor_trials\":["<<max_trials[0]<<","<<max_trials[1]<<","<<max_trials[2]
       <<"],\n\"independent_square_regression_bound\":"<<small<<",\n\"failures\":[";
    for(size_t i=0;i<failures.size();i++){if(i)out<<",";out<<"["<<failures[i].first<<","<<failures[i].second<<"]";}
    out<<"],\n\"elapsed_seconds\":"<<seconds<<"\n}\n";out.close();
    std::cout<<"COMPLETE norm_bound="<<bound<<" anchors="<<anchors.size()<<" covered="<<covered[0]<<","<<covered[1]<<","<<covered[2]<<" failures="<<failures.size()<<" seconds="<<seconds<<"\n";
    return failures.empty()?0:1;
}
