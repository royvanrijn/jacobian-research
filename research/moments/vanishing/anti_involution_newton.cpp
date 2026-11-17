// Bounded deterministic finite-field candidate search; not a certificate.
#include <array>
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

constexpr int P = 41;
using Point = std::array<int,6>;
struct Term { int c; std::array<int,6> e; };
int mod(int x) { x %= P; return x < 0 ? x+P : x; }
int power(int x,int n) { int y=1; while(n) {if(n&1)y=y*x%P;x=x*x%P;n>>=1;}return y; }

void evaluate(const std::vector<std::vector<Term>>& polynomials,const Point& x,
              Point& f,std::array<Point,6>& J) {
    int powers[6][13], inverse[6];
    for(int j=0;j<6;++j) {
        powers[j][0]=1;for(int e=1;e<=12;++e)powers[j][e]=powers[j][e-1]*x[j]%P;
        inverse[j]=x[j]?power(x[j],P-2):0;
    }
    for(int row=0;row<6;++row) {
        f[row]=0;J[row].fill(0);
        for(const auto& term:polynomials[row]) {
            int value=term.c;
            for(int j=0;j<6;++j)value=value*powers[j][term.e[j]]%P;
            f[row]=(f[row]+value)%P;
            for(int j=0;j<6;++j)if(term.e[j]) {
                int d=0;
                if(x[j])d=value*term.e[j]%P*inverse[j]%P;
                else if(term.e[j]==1) {
                    d=term.c;for(int k=0;k<6;++k)if(k!=j)d=d*powers[k][term.e[k]]%P;
                }
                J[row][j]=(J[row][j]+d)%P;
            }
        }
    }
}

bool solve(std::array<Point,6> J,Point rhs,Point& delta,int& determinant) {
    determinant=1;
    for(int col=0;col<6;++col) {
        int pivot=col;while(pivot<6&&!J[pivot][col])++pivot;
        if(pivot==6)return false;
        if(pivot!=col) {std::swap(J[pivot],J[col]);std::swap(rhs[pivot],rhs[col]);determinant=mod(-determinant);}
        int d=J[col][col], inv=power(d,P-2);determinant=determinant*d%P;
        for(int k=col;k<6;++k)J[col][k]=J[col][k]*inv%P;rhs[col]=rhs[col]*inv%P;
        for(int row=0;row<6;++row)if(row!=col&&J[row][col]) {
            int factor=J[row][col];
            for(int k=col;k<6;++k)J[row][k]=mod(J[row][k]-factor*J[col][k]);
            rhs[row]=mod(rhs[row]-factor*rhs[col]);
        }
    }
    delta=rhs;return true;
}

int main(int argc,char** argv) {
    if(argc!=3)return 2;
    std::ifstream input(argv[1]);
    std::vector<std::vector<Term>> polynomials(6);
    for(int row=0;row<6;++row) {
        int size;input>>size;polynomials[row].resize(size);
        for(auto& term:polynomials[row]) {input>>term.c;for(int& e:term.e)input>>e;}
    }
    if(!input)return 3;
    std::mt19937 random(20261116);
    auto start=std::chrono::steady_clock::now();
    int steps=0,done=0;
    for(int restart=0;restart<10000;++restart) {
        if(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()>40)break;
        Point x;for(int& v:x)v=random()%P;
        std::vector<Point> seen;
        ++done;
        for(int iteration=0;iteration<12;++iteration) {
            ++steps;Point f,delta;std::array<Point,6> J;
            evaluate(polynomials,x,f,J);int determinant;
            bool zero=true;for(int v:f)zero&=v==0;
            if(!solve(J,f,delta,determinant))break;
            if(zero) {
                int cyclic=mod((x[4]-x[5])*mod((x[4]+x[5])*(x[4]+x[5])))*power(8,P-2)%P;
                int scale=0;
                for(int l=1;l<P;++l)if(power(l,6)*8%P*cyclic%P==1) {scale=l;break;}
                if(cyclic) {
                    std::ofstream out(argv[2]);
                    out<<"{\"status\":\"candidate_requires_full_verification\",\"prime\":41,\"point\":[";
                    for(int j=0;j<6;++j)out<<(j?",":"")<<x[j];
                    out<<"],\"phase_scale\":"<<scale<<",\"cyclic_determinant\":"<<cyclic
                       <<",\"anti_jacobian_det\":"<<determinant<<",\"restarts\":"<<done
                       <<",\"steps\":"<<steps<<",\"seed\":20261116}\n";
                    std::cout<<"CANDIDATE restart "<<done<<" steps "<<steps<<"\n";return 0;
                }
                break;
            }
            seen.push_back(x);for(int j=0;j<6;++j)x[j]=mod(x[j]-delta[j]);
            bool cycle=false;for(const auto& old:seen)if(old==x)cycle=true;if(cycle)break;
        }
    }
    double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::ofstream out(argv[2]);
    out<<"{\"status\":\"no_candidate_in_bounded_search\",\"prime\":41,\"restarts\":"<<done
       <<",\"steps\":"<<steps<<",\"seconds\":"<<seconds<<",\"seed\":20261116,"
       <<"\"scope\":\"Not an exclusion; full Jacobian required for any candidate\"}\n";
    std::cout<<"No candidate; "<<done<<" restarts, "<<steps<<" steps, "<<seconds<<" seconds\n";
}
