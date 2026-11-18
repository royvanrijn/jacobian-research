#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
constexpr int p=37;
using Poly=array<int,14>;
using Tuple=array<int,7>;
int inv[p];
Poly C5{},C6{},D5{},D6{};
int md(int64_t x){int y=int(x%p);return y<0?y+p:y;}
Poly multiply(const Poly&A,int da,const Poly&B,int db){
 Poly out{};for(int i=0;i<=da;i++)for(int j=0;j<=db;j++)out[i+j]=md(out[i+j]+int64_t(A[i])*B[j]);return out;
}
Poly inverseOperator(const Poly&R,int sign){
 Poly out{};for(int j=12;j>=0;j--)out[j]=md((R[j]-int64_t(sign)*4*(j+1)*(j<12?out[j+1]:0))*inv[(2*j+1)%p]);return out;
}
Poly remainder(Poly T,const Poly&M,int degree){
 for(int j=12;j>=degree;j--){int coefficient=T[j];for(int k=0;k<degree;k++)T[j-degree+k]=md(T[j-degree+k]-int64_t(coefficient)*M[k]);T[j]=0;}return T;
}
void initialize(){
 for(int x=1;x<p;x++)for(int y=1;y<p;y++)if(x*y%p==1){inv[x]=y;break;}
 vector<Poly>C(7),D(7);C[0][0]=D[0][0]=1;C[1][0]=p-1;C[1][1]=1;D[1][0]=D[1][1]=1;
 for(int j=1;j<6;j++)for(int k=0;k<=j+1;k++){C[j+1][k]=md((k?C[j][k-1]:0)-C[j-1][k]);D[j+1][k]=md((k?D[j][k-1]:0)-D[j-1][k]);}
 C5=C[5];C6=C[6];D5=D[5];D6=D[6];
}
Poly numerator(const Poly&P,const Poly&Q,int r,int sign){
 Poly S=multiply(P,3,Q,3),factor{};
 for(int j=0;j<=6;j++)factor[j]=md(25*((sign==1?D6[j]:C6[j])-int64_t(r)*(sign==1?D5[j]:C5[j])));
 return inverseOperator(multiply(factor,6,S,6),sign);
}
int zeroConstraint(const Poly&P,const Poly&Q,int r){
 static const int binom[7][7]={{1},{1,1},{1,2,1},{1,3,3,1},{1,4,6,4,1},{1,5,10,10,5,1},{1,6,15,20,15,6,1}};
 Poly S=multiply(P,3,Q,3),Sigma{};
 for(int j=0;j<=6;j++)for(int k=0;k<=j;k++)Sigma[6+j-2*k]=md(Sigma[6+j-2*k]+S[j]*binom[j][k]);
 int sum=25,power=1;
 for(int j=1;j<=13;j++){int term=md(25*(Sigma[j-1]-int64_t(r)*(j<=12?Sigma[j]:0)));term=md(int64_t(term)*inv[(2*j-1)%p]);sum=md(sum+int64_t(term)*power);power=md(int64_t(power)*r);}
 return md(int64_t(sum)*inv[2]);
}
bool squareDivisible(const Poly&P,const Poly&Q,int r,int sign){
 Poly T=numerator(P,Q,r,sign),M=multiply(sign==1?P:Q,3,sign==1?P:Q,3);
 T=remainder(T,M,6);for(int j=0;j<6;j++)if(T[j])return false;return true;
}
struct Rows{int value[6][4];};
bool gauss(Rows&M,vector<int>&pivots,vector<int>&free){
 int rank=0;for(int c=0;c<3;c++){
  int k=rank;while(k<6&&!M.value[k][c])k++;
  if(k==6){free.push_back(c);continue;}
  if(k!=rank)for(int j=0;j<4;j++)swap(M.value[k][j],M.value[rank][j]);
  int inverse=inv[M.value[rank][c]];for(int j=c;j<4;j++)M.value[rank][j]=md(int64_t(M.value[rank][j])*inverse);
  for(int i=0;i<6;i++)if(i!=rank&&M.value[i][c]){int factor=M.value[i][c];for(int j=c;j<4;j++)M.value[i][j]=md(M.value[i][j]-int64_t(factor)*M.value[rank][j]);}
  pivots.push_back(c);rank++;
 }
 for(int i=rank;i<6;i++)if(M.value[i][3])return false;return true;
}
int determinant(int A[3][3]){
 return md(int64_t(A[0][0])*(A[1][1]*A[2][2]-A[1][2]*A[2][1])-int64_t(A[0][1])*(A[1][0]*A[2][2]-A[1][2]*A[2][0])+int64_t(A[0][2])*(A[1][0]*A[2][1]-A[1][1]*A[2][0]));
}
vector<array<int,3>> cramerOrEnumerate(const Rows&M){
 int A[3][3];for(int i=0;i<3;i++)for(int j=0;j<3;j++)A[i][j]=M.value[i][j];
 int det=determinant(A);vector<array<int,3>>out;
 if(det){array<int,3>q;for(int c=0;c<3;c++){int B[3][3];for(int i=0;i<3;i++)for(int j=0;j<3;j++)B[i][j]=(j==c?M.value[i][3]:A[i][j]);q[c]=md(int64_t(determinant(B))*inv[det]);}out.push_back(q);return out;}
 int pivotRow=-1,pivotColumn=-1;
 for(int i=0;i<3&&pivotRow<0;i++)for(int j=0;j<3;j++)if(M.value[i][j]){pivotRow=i;pivotColumn=j;break;}
 if(pivotRow<0){for(int i=0;i<3;i++)if(M.value[i][3])return out;for(int x=0;x<p;x++)for(int y=0;y<p;y++)for(int z=0;z<p;z++)out.push_back({x,y,z});return out;}
 vector<int>free;for(int j=0;j<3;j++)if(j!=pivotColumn)free.push_back(j);
 for(int x=0;x<p;x++)for(int y=0;y<p;y++){array<int,3>q{};q[free[0]]=x;q[free[1]]=y;q[pivotColumn]=md(int64_t(inv[M.value[pivotRow][pivotColumn]])*(M.value[pivotRow][3]-int64_t(M.value[pivotRow][free[0]])*x-int64_t(M.value[pivotRow][free[1]])*y));bool okay=true;for(int i=0;i<3;i++){int sum=-M.value[i][3];for(int j=0;j<3;j++)sum=md(sum+M.value[i][j]*q[j]);if(sum){okay=false;break;}}if(okay)out.push_back(q);}
 return out;
}
void vectorJson(const Poly&T,int size){cout<<"[";for(int j=0;j<size;j++){if(j)cout<<",";cout<<T[j];}cout<<"]";}
void controls(bool basis=false){
 vector<Tuple>tests={Tuple{1,2,3,4,5,6,7},Tuple{0,0,0,0,0,0,1},Tuple{6,12,8,0,0,0,0},Tuple{28,7,14,28,16,5,17}};
 if(basis){tests.clear();for(int pi=-1;pi<3;pi++)for(int qi=-1;qi<3;qi++)for(int rv=0;rv<=12;rv++){Tuple t{};if(pi>=0)t[pi]=1;if(qi>=0)t[3+qi]=1;t[6]=rv;tests.push_back(t);}}
 for(auto t:tests){Poly P{},Q{};P[0]=t[2];P[1]=t[1];P[2]=t[0];P[3]=1;Q[0]=t[5];Q[1]=t[4];Q[2]=t[3];Q[3]=1;cout<<"{\"tuple\":[";for(int j=0;j<7;j++){if(j)cout<<",";cout<<t[j];}cout<<"],\"X\":";vectorJson(numerator(P,Q,t[6],1),13);cout<<",\"Y\":";vectorJson(numerator(P,Q,t[6],-1),13);cout<<",\"G\":"<<zeroConstraint(P,Q,t[6])<<"}"<<endl;}
}
void census(bool weak){
 int64_t outer=0,candidates=0,ranks[4]={0,0,0,0};vector<Tuple>points;
 for(int a2=0;a2<p;a2++)for(int a1=0;a1<p;a1++)for(int a0=0;a0<p;a0++){
  Poly P{};P[0]=a0;P[1]=a1;P[2]=a2;P[3]=1;
  Poly divisor=weak?P:multiply(P,3,P,3);int degree=weak?3:6;
  Poly A=multiply(D6,6,P,3),B=multiply(D5,5,P,3);int basis[2][6][4]={};
  for(int mode=0;mode<2;mode++)for(int j=0;j<4;j++){Poly rhs{};for(int k=0;k<10;k++)rhs[k+j]=md((mode?-25:25)*(mode?B[k]:A[k]));Poly X=remainder(inverseOperator(rhs,1),divisor,degree);for(int k=0;k<degree;k++)basis[mode][k][j]=X[k];}
  // G at r=0 is25, so that omitted parameter cannot solve the full system.
  for(int rv=1;rv<p;rv++){
   outer++;Rows M{};for(int k=0;k<degree;k++)for(int j=0;j<4;j++){int co=md(basis[0][k][j]+int64_t(rv)*basis[1][k][j]);M.value[k][j]=j==3?md(-co):co;}
   vector<array<int,3>>solutions;
   if(weak)solutions=cramerOrEnumerate(M);
   else{vector<int>pivots,free;if(!gauss(M,pivots,free))continue;ranks[pivots.size()]++;int size=1;for(size_t j=0;j<free.size();j++)size*=p;for(int code=0;code<size;code++){array<int,3>q{};int remain=code;for(int column:free){q[column]=remain%p;remain/=p;}for(size_t k=0;k<pivots.size();k++){int value=M.value[k][3];for(int column:free)value=md(value-int64_t(M.value[k][column])*q[column]);q[pivots[k]]=value;}solutions.push_back(q);}}
   for(auto q:solutions){candidates++;Poly Q{};Q[0]=q[0];Q[1]=q[1];Q[2]=q[2];Q[3]=1;if(!squareDivisible(P,Q,rv,1)||!squareDivisible(P,Q,rv,-1)||zeroConstraint(P,Q,rv))continue;points.push_back({a2,a1,a0,q[2],q[1],q[0],rv});}
  }
 }
 cout<<"{\"mode\":\""<<(weak?"weak_cramer":"strong_gauss")<<"\",\"outer\":"<<outer<<",\"candidates\":"<<candidates<<",\"points\":[";
 for(size_t i=0;i<points.size();i++){if(i)cout<<",";cout<<"[";for(int j=0;j<7;j++){if(j)cout<<",";cout<<points[i][j];}cout<<"]";}cout<<"],\"consistent_ranks\":";if(weak)cout<<"null";else{cout<<"[";for(int j=0;j<4;j++){if(j)cout<<",";cout<<ranks[j];}cout<<"]";}cout<<"}"<<endl;
}
int main(int argc,char**argv){initialize();if(argc>1&&string(argv[1])=="--controls"){controls();return 0;}if(argc>1&&string(argv[1])=="--basis-controls"){controls(true);return 0;}census(false);census(true);}
