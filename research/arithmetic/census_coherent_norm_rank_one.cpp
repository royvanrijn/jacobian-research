// Exact positive odd-coefficient census for the coherent three-prime norm records.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <vector>
using integer = std::int64_t;
using wide = __int128_t;

integer mod(integer a, integer p) {
    const integer r = a % p;
    return r >= 0 ? r : r+p;
}
integer inv(integer a, integer p) {
    if (p == 1) return 0;
    integer b = p, x = 1, y = 0;
    while (b != 0) {
        integer q = a / b;
        integer r = a % b;
        a = b;
        b = r;
        integer c = x-q*y;
        x = y;
        y = c;
    }
    return mod(x,p);
}
integer isqrt(integer a) {
    integer r = std::sqrt(static_cast<long double>(a));
    while (wide(r+1)*(r+1) <= a) ++r;
    while (wide(r)*r > a) --r;
    return r;
}
std::vector<int> primes_up_to(int n) {
    std::vector<bool> composite(n+1);
    std::vector<int> primes;
    for (int p=2; p<=n; ++p) {
        if (composite[p]) continue;
        primes.push_back(p);
        if (p<=n/p) for (int m=p*p; m<=n; m+=p) composite[m]=true;
    }
    return primes;
}
bool rank_one(int Z, std::array<int,3> p, int x, int y, int z) {
    const integer rad2 = 12LL*Z*Z;
    const integer radius = isqrt(rad2);
    const integer c1 = mod(p[0]*inv(p[1],x),x);
    const integer c3 = mod(p[2]*inv(p[1],y),y);
    for (integer j=-radius; j<=radius; ++j)
        for (integer i=-radius+mod(c1*j+radius,x); i<=radius; i+=x)
            for (integer k=-radius+mod(c3*j+radius,y); k<=radius; k+=y) {
                if (i*i+j*j+k*k > rad2 || (p[2]*i-p[0]*k)%z != 0) continue;
                if (p[1]*i != p[0]*j || p[2]*j != p[1]*k) return false;
            }
    return true;
}
bool squarefree(integer d, const std::vector<int>& primes) {
    for (integer p:primes) {
        if (p*p>d) break;
        if (d%(p*p)==0) return false;
    }
    return true;
}
bool in_field_window(integer E, integer X, const std::vector<int>& primes) {
    for (integer g=1; g<=isqrt(E/X); ++g)
        if (E%(g*g)==0 && E/(g*g)<2*X && squarefree(E/(g*g),primes))
            return true;
    return false;
}
void run(int Z, bool annulus) {
    auto ps=primes_up_to(2*Z-1);
    ps.erase(ps.begin(),std::lower_bound(ps.begin(),ps.end(),Z));
    const int P=ps.size();
    const integer zz=1LL*Z*Z*Z;
    const integer X=isqrt(zz*zz*Z);
    const integer V=isqrt(isqrt(zz*Z*Z));
    const int lo=annulus?V:1, hi=2*V;
    const auto sq_primes=primes_up_to(isqrt(2*X));
    std::vector<integer> powers;
    for (integer p:ps) powers.push_back(p*p*p);
    std::vector<std::vector<std::vector<int>>> roots(P*P, std::vector<std::vector<int>>(hi+1));
    for (int i=0;i<P;++i) for (int j=0;j<P;++j) if (i!=j)
        for (int m=lo|1;m<=hi;m+=2) {
            if (std::gcd(m,ps[i])!=1 || std::gcd(m,ps[j])!=1) continue;
            const integer target=mod(mod(powers[i],m)*inv(mod(powers[j],m),m),m);
            for (integer r=0;r<m;++r) if (r*r%m==target) roots[i*P+j][m].push_back(r);
        }
    struct Pair { int x, y, inverse; };
    std::vector<Pair> pairs;
    for (int x=lo|1;x<=hi;x+=2) for (int y=lo|1;y<=hi;y+=2)
        if (std::gcd(x,y)==1) pairs.push_back({x,y,static_cast<int>(inv(x,y))});
    integer total=0, rank=0, pos=0, unit=0, win=0, rankwin=0;
    for (int i=0;i<P;++i) for (int j=i+1;j<P;++j) for (int k=j+1;k<P;++k)
        for (const auto [x,y,xinv]:pairs) {
            // z/y is a square root of p_i^3/p_j^3 modulo x;
            // z/x is a square root of p_k^3/p_j^3 modulo y.
            const auto& u=roots[i*P+j][x];
            const auto& v=roots[k*P+j][y];
            if (u.empty() || v.empty()) continue;
            const integer period=1LL*x*y;
            const wide A=wide(powers[i])*y*y, C=wide(powers[k])*x*x;
            if (A==C) continue;
            for (integer r:u) {
                const integer v1=y*r%x;
                for (integer s:v) {
                    const integer z0=v1+x*mod((x*s-v1)*xinv,y);
                    // All z in this CRT class are coprime to x and y.
                    for (integer z=lo+mod(z0-lo,period);z<=hi;z+=period) {
                        if (z%2==0 || z%ps[i]==0 || z%ps[k]==0 || (A-C)%z!=0) continue;
                        ++total;
                        const bool one=rank_one(Z,{ps[i],ps[j],ps[k]},x,y,z);
                        if (one) ++rank;
                        const wide B=wide(powers[j])*z*z;
                        const wide num=4*A*C-(A+C-B)*(A+C-B);
                        const wide t=wide(x)*y*z;
                        if (num%(t*t)!=0) std::abort();
                        const integer E=num/(t*t);
                        if (E<=0) continue;
                        if (one) ++pos;
                        if (std::gcd(E,1LL*ps[i]*ps[j]*ps[k])!=1) continue;
                        if (one) ++unit;
                        if (in_field_window(E,X,sq_primes)) {
                            ++win;
                            if (one) ++rankwin;
                        }
                    }
                }
            }
        }
    std::cout<<Z<<" "<<P<<" "<<lo<<" "<<hi<<" "<<X<<" "<<total<<" "<<rank<<" "
             <<pos<<" "<<unit<<" "<<win<<" "<<rankwin<<std::endl;
}
int main() {
    constexpr std::array<int, 8> scales{16, 24, 32, 48, 64, 96, 128, 192};
    std::cout << "Z prime_count lo hi X all rank1 positive_rank1 unit_rank1 window_all window_rank1\n";
    for (int Z : scales) {
        run(Z, true);
        if (Z <= 128) run(Z, false);
    }
}

