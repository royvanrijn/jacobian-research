#include <algorithm>
#include <cstdint>
#include <vector>
#include <array>

using Int=__int128_t;
constexpr Int one=Int(1)<<48;
constexpr int max_coefficients=14;
struct Interval { Int lo,hi; };
struct Words { std::uint64_t low; std::int64_t high; };

Int integer(Words w) { return Int(w.high)*(Int(1)<<64) + Int(w.low); }
Int floor_grid(Int x) { return x>>48; }
Int ceil_grid(Int x) { return -floor_grid(-x); }

// The input endpoint x represents x/2^48. For a product, keep the
// quotient and remainder upon dividing its raw value by 2^48 separately.
// Splitting both operands before multiplication avoids any 256-bit value.
struct Product {
    Int integer,remainder;
    static Product of(Int x,Int y) {
        Int a=floor_grid(x), b=floor_grid(y);
        Int h=x-a*one, k=y-b*one;
        Int l=h*k;
        return {a*b*one + a*k + b*h + floor_grid(l),l&(one-1)};
    }
    bool operator<(Product other) const {
        return integer<other.integer || (integer==other.integer && remainder<other.remainder);
    }
};

struct ProductSum {
    Int integer=0,remainder=0;
    void add(Product p) { integer+=p.integer; remainder+=p.remainder; }
    Int floor() const { return integer+floor_grid(remainder); }
    Int ceil() const { return integer+ceil_grid(remainder); }
    bool positive() const { return floor()>0 || (floor()==0 && (remainder&(one-1))); }
    bool negative() const { return floor()<0; }
};

void accumulate(ProductSum &lo, ProductSum &hi, Interval a, Interval b) {
    if ((a.lo==0 && a.hi==0) || (b.lo==0 && b.hi==0)) return;
    Product ll=Product::of(a.lo,b.lo), lh=Product::of(a.lo,b.hi);
    Product hl=Product::of(a.hi,b.lo), hh=Product::of(a.hi,b.hi);
    lo.add(std::min(std::min(ll,lh),std::min(hl,hh)));
    hi.add(std::max(std::max(ll,lh),std::max(hl,hh)));
}

class Projection {
public:
    int n,g;
    std::vector<Interval> parameters;
    std::array<int,max_coefficients> weights{};
    bool has_sign(const Interval *pivot) const {
        bool pos=true,neg=true;
        for(int k=0;k<n;++k) {
            ProductSum lo,hi;
            for(int j=0;j<g;++j) accumulate(lo,hi,parameters[j],pivot[j*n+k]);
            pos=pos&&lo.positive();
            neg=neg&&hi.negative();
            if (!pos&&!neg) return false;
        }
        return true;
    }
    bool midpoint_sign(const Interval *pivot) const {
        ProductSum lo,hi;
        for(int j=0;j<g;++j) {
            Interval mid{0,0};
            for(int k=0;k<n;++k) {
                mid.lo+=pivot[j*n+k].lo*weights[k];
                mid.hi+=pivot[j*n+k].hi*weights[k];
            }
            accumulate(lo,hi,parameters[j],mid);
        }
        return lo.positive()||hi.negative();
    }
    int verify(const std::vector<Interval> &data,int depth) const {
        for(int r=0;r<n;++r) if(has_sign(&data[r*g*n])) return 1;
        if(depth>=14) return -1;
        bool possible=false;
        for(int r=0;r<n;++r) if(midpoint_sign(&data[r*g*n])) { possible=true; break; }
        if(!possible) return -1;
        std::vector<Interval> left(data.size()),right(data.size());
        for(int i=0;i<n*g;++i) {
            std::array<Interval,max_coefficients> work;
            for(int j=0;j<n;++j) work[j]=data[i*n+j];
            for(int l=n;l>0;--l) {
                left[i*n+n-l]=work[0];
                right[i*n+l-1]=work[l-1];
                for(int j=0;j+1<l;++j) {
                    work[j]={(work[j].lo+work[j+1].lo)>>1,-((-(work[j].hi+work[j+1].hi))>>1)};
                }
            }
        }
        int a=verify(right,depth+1);
        if(a<0) return -1;
        int b=verify(left,depth+1);
        return b<0 ? -1 : a+b;
    }
};

int sfc_internal(int n,int groups,const std::int64_t *sine,
        const Words *transform,const std::int64_t *model,
        const std::int64_t *center, const std::int64_t *derivatives) {
    if(n<2 || n>max_coefficients || groups<1 || groups>256 || (Int(-2)>>1)!=-1) return -2;
    auto read=[&](const std::int64_t *a,int i){return Interval{a[2*i],a[2*i+1]};};
    auto limit=[](Interval a,Int mx) {return -mx<=a.lo && a.lo<=a.hi && a.hi<=mx;};
    Projection p{n,groups,{}};
    for(int k=0,s=1;k<n;++k) {p.weights[k]=s; s=s*(n-1-k)/(k+1);}
    for(int j=0;j<groups;++j) {
        auto a=read(model,j*(n+1));
        if(!limit(a,2*one)) return -2;
        p.parameters.push_back(a);
        for(int k=0;k<n;++k) if(!limit(read(model,j*(n+1)+k+1),512*one)) return -2;
    }
    std::array<Interval,max_coefficients*max_coefficients> tf;
    const int transform_bits = n==14 ? 83 : 84;
    for(int t=0;t<n*n;++t) {
        tf[t]={integer(transform[2*t]),integer(transform[2*t+1])};
        if(n>=12) {
            // A common positive scale leaves every projected sign
            // unchanged. Round its Bernstein entries outwards.
            const int shift = n==12 ? 4 : (n==13 ? 8 : 13);
            if(!limit(tf[t],Int(1)<<(transform_bits+shift))) return -2;
            tf[t]={tf[t].lo>>shift,-((-tf[t].hi)>>shift)};
        }
        if(!limit(tf[t],Int(1)<<transform_bits)) return -2;
    }
    for(int k=0;k<n;++k) {
        if(!limit(read(sine,k),32*one)) return -2;
        if(center && (!limit(read(center,k),512*one)
                    ||!limit(read(derivatives,k),512*one))) return -2;
    }
    auto product=[&](Interval c,const std::int64_t *directions,int k,int r) {
        auto s=read(directions,k>=r?k-r:r-k);
        if(k<r) s={-s.hi,-s.lo};
        auto ll=s.lo*c.lo,lh=s.lo*c.hi,hl=s.hi*c.lo,hh=s.hi*c.hi;
        return Interval{floor_grid(std::min(std::min(ll,lh),std::min(hl,hh))),
                         ceil_grid(std::max(std::max(ll,lh),std::max(hl,hh)))};
    };
    // After division by the grid, the sine and linear products
    // together are bounded by 2^62+2^66. The transform bound is
    // 2^84, so each resulting control is at most
    // For n<=13 the controls are bounded by 13*(2^62+2^66)*2^36
    // < 2^106. At n=14 the tighter transform bound instead gives
    // 14*(2^62+2^66)*2^35 < 2^105. The midpoint weights sum to
    // 2^(n-1). In either case, multiplying by parameters bounded
    // by 2 and summing at most 256 groups stays below 2^127.
    std::vector<Interval> result(n*n*groups);
    for(int r=0;r<n;++r) {
        for(int j=0;j<groups;++j) {
            std::array<Interval,max_coefficients> products;
            for(int k=0;k<n;++k) {
                auto c=read(model,j*(n+1)+k+1);
                products[k]=product(c,(center && j==0)?center:sine,k,r);
                if(center && j==1) {
                    auto extra=product(read(model,k+1),derivatives,k,r);
                    products[k].lo+=extra.lo;
                    products[k].hi+=extra.hi;
                }
            }
            for(int i=0;i<n;++i) {
                ProductSum lo,hi;
                for(int k=0;k<n;++k) {
                    auto a=tf[i*n+k],b=products[k];
                    if(a.lo>=0) {
                        lo.add(Product::of(b.lo>=0?a.lo:a.hi,b.lo));
                        hi.add(Product::of(b.hi>=0?a.hi:a.lo,b.hi));
                    } else {
                        accumulate(lo,hi,a,b);
                    }
                }
                result[(r*groups+j)*n+i]={lo.floor(),hi.ceil()};
            }
        }
        if(p.has_sign(&result[r*groups*n])) return 1;
    }
    return p.verify(result,0);
}


extern "C" int sfc_fused_radial(int n,int groups,const std::int64_t *sine,
        const Words *transform,const std::int64_t *model,
        const std::int64_t *center,const std::int64_t *derivatives) {
    return sfc_internal(n,groups,sine,transform,model,center,derivatives);
}

extern "C" int sfc_radial(int n,int groups,const std::int64_t *sine,
        const Words *transform,const std::int64_t *model) {
    return sfc_internal(n,groups,sine,transform,model,nullptr,nullptr);
}
