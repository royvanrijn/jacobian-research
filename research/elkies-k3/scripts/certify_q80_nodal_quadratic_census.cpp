// Complete 131^3 polynomial census with the simple root on the nodal quadratic.
#include <array>
#include <ctime>
#include <iostream>
#include <sys/resource.h>
using P = std::array<int,13>;
constexpr int p=131;
int mod(int n) { n%=p; return n<0?n+p:n; }
P mul(const P&a,const P&b) {
    P c{};
    for(int i=0;i<13;i++) if(a[i])
        for(int j=0;j+i<13;j++) if(b[j]) c[i+j]+=a[i]*b[j];
    for(int&i:c)i%=p;
    return c;
}
void print_poly(const P&a) {
    int n=12;while(n>0&&!a[n])--n;
    std::cout<<'[';
    for(int i=0;i<=n;i++){if(i)std::cout<<',';std::cout<<a[i];}
    std::cout<<']';
}
int main() {
    rlimit cpu{30,35},mem{1024UL*1024*1024,1024UL*1024*1024};
    setrlimit(RLIMIT_CPU,&cpu);setrlimit(RLIMIT_AS,&mem);
    auto start=std::clock();P A{},B{};int r0,r1;
    for(int i=0;i<=8;i++)std::cin>>A[i];
    for(int i=0;i<=12;i++)std::cin>>B[i];
    std::cin>>r0>>r1;if(!std::cin)return 2;
    int inv[p]{};for(int a=1;a<p;a++)for(int b=1;b<p;b++)if(a*b%p==1)inv[a]=b;
    bool first=true;long tried=0,accepted=0;
    std::cout<<"{\"status\":\"PASS\",\"records\":[";
    for(int z2=0;z2<p;z2++)for(int z1=0;z1<p;z1++)for(int z0=0;z0<p;z0++) {
        ++tried;P X{},S{},Y{};
        X[0]=mod(r0+88*z0);X[1]=mod(r1+62*z0+88*z1);
        X[2]=mod(z0+62*z1+88*z2);X[3]=mod(z1+62*z2);X[4]=z2;
        auto cube=mul(mul(X,X),X),ax=mul(A,X);P f{};
        for(int i=0;i<13;i++)f[i]=mod(cube[i]+ax[i]+B[i]);
        for(int i=12;i>=2;i--){int c=f[i];S[i-2]=c;f[i]=0;f[i-1]=mod(f[i-1]-62*c);f[i-2]=mod(f[i-2]-88*c);}
        if(f[0]||f[1])return 3;
        int degree=10;while(degree>=0&&!S[degree])--degree;
        if(degree<0)return 4;
        if(degree%2)continue;
        int scalar=S[degree],n=degree/2;Y[n]=1;
        for(int k=n-1;k>=0;k--) {
            int coefficient=0;
            for(int i=k+1;i<=n;i++){int j=n+k-i;if(j>k&&j<=n)coefficient+=Y[i]*Y[j];}
            Y[k]=mod(mod(S[n+k]*inv[scalar])-coefficient)*66%p;
        }
        auto square=mul(Y,Y);bool good=true;
        for(int k=0;k<=degree;k++)if(square[k]*scalar%p!=S[k]){good=false;break;}
        if(!good)continue;
        if(!first)std::cout<<',';first=false;++accepted;
        std::cout<<"{\"z\":["<<z0<<','<<z1<<','<<z2<<"],\"x\":";print_poly(X);
        std::cout<<",\"scalar\":"<<scalar<<",\"s\":";print_poly(Y);std::cout<<'}';
    }
    std::cout<<"],\"tried\":"<<tried<<",\"accepted\":"<<accepted
             <<",\"cpu_seconds\":"<<double(std::clock()-start)/CLOCKS_PER_SEC<<"}\n";
}
