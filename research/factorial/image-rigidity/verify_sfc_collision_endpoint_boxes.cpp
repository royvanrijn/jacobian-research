// Search the order 2, 3, 4 window for ordinary three-column factorial supports.
// Usage: binary 0 radius | binary 1 radius | binary split radius | binary filter
// Scan commands print only the modular survivors to stdout, one support per line.
// Filter reads this format from stdin and rescans the supports at its own prime.
#include <algorithm>
#include <array>
#include <cassert>
#include <charconv>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string_view>
#include <vector>

#ifndef PRIME
#define PRIME 65521
#endif
constexpr int p = PRIME;
constexpr bool prime() {
    for (int divisor = 2; divisor * divisor <= p; ++divisor)
        if (p % divisor == 0) return false;
    return true;
}
static_assert(p > 4 && p < 65536 && prime());

int multiply(int a, int b) { return (static_cast<uint32_t>(a) * b) % p; }
int subtract(int a, int b) { return a >= b ? a-b : p+a-b; }
std::array<int,p> inverse;
std::array<int,p> factorial;

// Coefficients of a binary homogeneous form, indexed by the exponent of y.
// GCD treats them as a univariate polynomial after setting z=1.
struct Polynomial {
    int degree;
    std::array<int,9> c{};
};
Polynomial operator+(const Polynomial& a, const Polynomial& b) {
    assert(a.degree == b.degree);
    Polynomial sum{a.degree};
    for (int k=0; k<=sum.degree; ++k) {
        int v=a.c[k]+b.c[k];
        sum.c[k]=v < p ? v : v-p;
    }
    return sum;
}
Polynomial operator-(const Polynomial& a, const Polynomial& b) {
    assert(a.degree == b.degree);
    Polynomial difference{a.degree};
    for (int k=0; k<=difference.degree; ++k) difference.c[k]=subtract(a.c[k],b.c[k]);
    return difference;
}
Polynomial operator*(const Polynomial& a, const Polynomial& b) {
    Polynomial product{a.degree+b.degree};
    assert(product.degree <= 8);
    std::array<uint64_t,9> accumulator{};
    for (int i=0; i<=a.degree; ++i)
        for (int j=0; j<=b.degree; ++j)
            accumulator[i+j] += static_cast<uint32_t>(a.c[i])*b.c[j];
    for (int k=0; k<=product.degree; ++k) product.c[k]=accumulator[k]%p;
    return product;
}
Polynomial operator*(int scalar, const Polynomial& a) {
    Polynomial product{a.degree};
    for (int k=0; k<=a.degree; ++k) product.c[k]=multiply(scalar,a.c[k]);
    return product;
}
void trim(Polynomial& a) {
    while (a.degree>=0 && a.c[a.degree]==0) --a.degree;
}
Polynomial gcd(Polynomial a, Polynomial b) {
    trim(a);
    trim(b);
    while (b.degree>=0) {
        if (b.degree==0) return b;
        const int inv_lead=inverse[b.c[b.degree]];
        for (int shift=a.degree-b.degree; shift>=0; --shift) {
            const int leading=multiply(a.c[b.degree+shift],inv_lead);
            for (int j=0; j<b.degree; ++j) {
                a.c[j+shift]=subtract(a.c[j+shift],multiply(leading,b.c[j]));
            }
        }
        a.degree=std::min(a.degree,b.degree-1);
        trim(a);
        std::swap(a,b);
    }
    return a;
}

using Coefficients = std::array<int,31>;
constexpr int index(int n,int i,int j) {
    const int prefix=n==2 ? 0 : n==3 ? 6 : 16;
    return prefix+i*(n+1)-i*(i-1)/2+j;
}
Polynomial slice(const Coefficients& coefficients,int n,int i) {
    Polynomial value{n-i};
    for (int j=0; j<=n-i; ++j) value.c[j]=coefficients[index(n,i,j)];
    return value;
}

Polynomial quadratic_norm(const Polynomial& B,const Polynomial& C,
                          const Polynomial& linear,const Polynomial& constant) {
    return constant*constant-B*(linear*constant)+C*(linear*linear);
}

// Normalize F2 to x^2 + Bx + C, then reduce F3 and F4 modulo F2.
// Their remainders are Dx+E and Gx+H. Since F2(1,0,0) is nonzero,
// every projective candidate determines [y:z]. The two linear equations
// must be consistent, and each must share a root with the quadratic.
// These conditions are sufficient over the algebraic closure.
bool has_common_root(const Coefficients& coefficients) {
    assert(coefficients[index(2,2,0)]!=0);
    const int normalize=inverse[coefficients[index(2,2,0)]];
    const auto B=normalize*slice(coefficients,2,1);
    const auto C=normalize*slice(coefficients,2,0);
    const auto B2minusC=B*B-C;
    const auto BC=B*C;
    const auto D=slice(coefficients,3,1)-B*slice(coefficients,3,2)
        +coefficients[index(3,3,0)]*B2minusC;
    const auto E=slice(coefficients,3,0)-C*slice(coefficients,3,2)
        +coefficients[index(3,3,0)]*BC;
    const auto G=slice(coefficients,4,1)-B*slice(coefficients,4,2)
        +B2minusC*slice(coefficients,4,3)
        -coefficients[index(4,4,0)]*(B*B2minusC-BC);
    const auto H=slice(coefficients,4,0)-C*slice(coefficients,4,2)
        +BC*slice(coefficients,4,3)
        -coefficients[index(4,4,0)]*(B2minusC*C);
    const auto determinant=D*H-E*G;
    const auto norm3=quadratic_norm(B,C,D,E);
    // The coefficients of y^6 and y^8 test the missing point [y:z]=[1:0].
    const bool possible_infinity=determinant.c[6]==0 && norm3.c[6]==0;
    const auto common=gcd(determinant,norm3);
    if (possible_infinity || common.degree!=0) {
        const auto norm4=quadratic_norm(B,C,G,H);
        return (possible_infinity && norm4.c[8]==0) || gcd(common,norm4).degree!=0;
    }
    return false;
}

using Row=std::array<int,3>;
using Support=std::vector<Row>;
struct Option {
    Row row;
    Coefficients factors;
    unsigned equal_columns;
};
int equal_columns(Row row) {
    return (row[0]==row[1]) + 2*(row[0]==row[2]) + 4*(row[1]==row[2]);
}
Coefficients row_factors(Row row) {
    Coefficients factors;
    for (int n=2; n<=4; ++n)
        for (int i=0; i<=n; ++i)
            for (int j=0; j<=n-i; ++j) {
                const int exponent=i*row[0]+j*row[1]+(n-i-j)*row[2];
                if (exponent<0 || exponent>=p) throw std::out_of_range("factorial argument");
                factors[index(n,i,j)]=factorial[exponent];
            }
    return factors;
}
Coefficients multinomials() {
    Coefficients result;
    for (int n=2; n<=4; ++n)
        for (int i=0; i<=n; ++i)
            for (int j=0; j<=n-i; ++j)
                result[index(n,i,j)]=factorial[n]/(factorial[i]*factorial[j]*factorial[n-i-j]);
    return result;
}
Coefficients apply_row(const Coefficients& coefficients,const Coefficients& factors) {
    Coefficients product;
    for (size_t i=0; i<product.size(); ++i) product[i]=multiply(coefficients[i],factors[i]);
    return product;
}
void print_support(const Support& support) {
    std::printf("%zu",support.size());
    for (const auto& row:support) std::printf(" %d %d %d",row[0],row[1],row[2]);
    std::puts("");
}
struct Scan {
    std::vector<std::vector<Option>> options;
    Support support;
    uint64_t count=0,duplicates=0,survivors=0;

    Scan(const Support& center,int radius): options(center.size()),support(center.size()) {
        for (size_t i=0; i<center.size(); ++i) {
            const auto [a,b,c]=center[i];
            for (int x=std::max(0,a-radius); x<=a+radius; ++x)
                for (int y=std::max(0,b-radius); y<=b+radius; ++y)
                    for (int z=std::max(0,c-radius); z<=c+radius; ++z) {
                        const Row row{x,y,z};
                        options[i].push_back({row,row_factors(row),static_cast<unsigned>(equal_columns(row))});
                    }
        }
    }

    void walk(size_t index,const Coefficients& coefficients,unsigned equal_mask=7) {
        if (index==support.size()) {
            ++count;
            if (equal_mask) {
                ++duplicates;
            } else if (has_common_root(coefficients)) {
                ++survivors;
                print_support(support);
            }
            return;
        }
        for (const auto& option:options[index]) {
            support[index]=option.row;
            walk(index+1,apply_row(coefficients,option.factors),equal_mask&option.equal_columns);
        }
    }
};
int read_radius(std::string_view input) {
    int value = -1;
    const auto result = std::from_chars(input.begin(), input.end(), value);
    return result.ec == std::errc{} && result.ptr == input.end() && value >= 0 ? value : -1;
}

int main(int argc,char** argv) {
    factorial[0]=1;
    inverse[1]=1;
    for (int i=1; i<p; ++i) factorial[i]=multiply(factorial[i-1],i);
    for (int i=2; i<p; ++i) inverse[i]=p-multiply(p/i,inverse[p%i]);
    const auto weights=multinomials();

    if (argc==2 && std::string_view(argv[1])=="filter") {
        size_t rows;
        unsigned long long count=0,survivors=0;
        while (std::scanf("%zu",&rows)==1) {
            Support support(rows);
            auto coefficients=weights;
            for (auto& row:support) {
                if (std::scanf("%d%d%d",&row[0],&row[1],&row[2])!=3) return 2;
                coefficients=apply_row(coefficients,row_factors(row));
            }
            ++count;
            if (has_common_root(coefficients)) {
                ++survivors;
                print_support(support);
            }
        }
        std::fprintf(stderr,"p=%d filtered=%llu survivors=%llu\n",p,count,survivors);
        return 0;
    }

    // Matrix at q=0 and the matrix obtained when the first two rows merge at q=1.
    const Support endpoint0{{5,2,3},{2,5,3},{5,4,4},{4,5,4},{5,5,4}};
    const Support endpoint1{{7,7,6},{5,4,4},{4,5,4},{5,5,4}};

    if (argc==3 && std::string_view(argv[1])=="split") {
        const int radius=read_radius(argv[2]);
        if (radius<0 || radius>6) return 2;
        Scan scan(endpoint0,1);
        const auto [a,b,c]=endpoint1.front();
        uint64_t splits=0;
        for (int x=a-radius; x<=a+radius; ++x)
            for (int y=b-radius; y<=b+radius; ++y)
                for (int z=c-radius; z<=c+radius; ++z)
                    for (int i=0; i<=x; ++i)
                        for (int j=0; j<=y; ++j)
                            for (int k=0; k<=z; ++k) {
                                const Row row1{i,j,k};
                                const Row row2{x-i,y-j,z-k};
                                // Ordinary moments are invariant under swapping these rows.
                                if (row2<row1) continue;
                                ++splits;
                                scan.support[0]=row1;
                                scan.support[1]=row2;
                                const auto coefficients=apply_row(
                                    apply_row(weights,row_factors(row1)),row_factors(row2));
                                scan.walk(2,coefficients,equal_columns(row1)&equal_columns(row2));
                            }
        std::fprintf(stderr,
            "p=%d split radius=%d unordered_splits=%llu tuples=%llu duplicates=%llu survivors=%llu\n",
            p,radius,static_cast<unsigned long long>(splits),
            static_cast<unsigned long long>(scan.count),
            static_cast<unsigned long long>(scan.duplicates),
            static_cast<unsigned long long>(scan.survivors));
        return 0;
    }
    if (argc!=3 || std::string_view(argv[1]).size()!=1
        || (argv[1][0]!='0' && argv[1][0]!='1')) return 2;
    const int radius=read_radius(argv[2]);
    if (radius<0 || radius>(p-1)/4-7) return 2;
    Scan scan(argv[1][0]=='0' ? endpoint0 : endpoint1,radius);
    scan.walk(0,weights);
    std::fprintf(stderr,"p=%d endpoint=%s radius=%d tuples=%llu duplicates=%llu survivors=%llu\n",
        p,argv[1],radius,static_cast<unsigned long long>(scan.count),
        static_cast<unsigned long long>(scan.duplicates),static_cast<unsigned long long>(scan.survivors));
}

