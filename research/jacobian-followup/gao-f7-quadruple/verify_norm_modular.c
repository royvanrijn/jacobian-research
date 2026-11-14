/* GPL-2.0-or-later. Exact modular interpolation, no multivariate resultant.
 * Every determinant value is checked by FLINT and elementary C elimination. */
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_mat.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define N 15
#define CHECK(c) do { if (!(c)) { fprintf(stderr,"failed line %d\n",__LINE__); exit(2); } } while (0)
static uint64_t inv_plain(uint64_t a,uint64_t p){
 int64_t t=0,nt=1,r=(int64_t)p,nr=(int64_t)a;
 while(nr){int64_t q=r/nr,z=t-q*nt;t=nt;nt=z;z=r-q*nr;r=nr;nr=z;}
 CHECK(r==1);if(t<0)t+=(int64_t)p;return (uint64_t)t;
}
static uint64_t det_plain(uint64_t W[N][N],uint64_t p){
 uint64_t det=1;
 for(int k=0;k<N;k++){
  int pivot=k;while(pivot<N && W[pivot][k]==0)pivot++;
  if(pivot==N)return 0;
  if(pivot!=k){for(int j=k;j<N;j++){uint64_t t=W[k][j];W[k][j]=W[pivot][j];W[pivot][j]=t;}det=det ? p-det:0;}
  uint64_t v=W[k][k];det=det*v%p;uint64_t vi=inv_plain(v,p);
  for(int i=k+1;i<N;i++){
   uint64_t m=W[i][k]*vi%p;if(!m)continue;W[i][k]=0;
   for(int j=k+1;j<N;j++){uint64_t sub=m*W[k][j]%p;W[i][j]=(W[i][j]+p-sub)%p;}
  }
 }
 return det;
}
static int prime_trial(uint64_t p){if(p<2 || p%2==0)return p==2;for(uint64_t d=3;d*d<=p;d+=2)if(p%d==0)return 0;return 1;}
int main(int argc,char **argv){
 CHECK(argc==2 || argc==3);FILE *fd=fopen(argv[1],"r");CHECK(fd);
 long points,np;CHECK(fscanf(fd,"%ld %ld",&points,&np)==2);CHECK(points==11670 && np==345);
 ulong *primes=malloc((size_t)np*sizeof(ulong));CHECK(primes);
 for(long i=0;i<np;i++)CHECK(fscanf(fd,"%lu",&primes[i])==1);
 fmpz_poly_t polys[18];fmpz_t x;fmpz_init(x);
 for(int i=0;i<18;i++){
  fmpz_poly_init(polys[i]);long len;CHECK(fscanf(fd,"%ld",&len)==1 && len>=0 && len<=11670);
  for(long j=0;j<len;j++){CHECK(fmpz_fread(fd,x)>0);fmpz_poly_set_coeff_fmpz(polys[i],j,x);}
 }
 fclose(fd);fmpz_clear(x);long use=np;if(argc==3){use=strtol(argv[2],NULL,10);CHECK(use>0 && use<=np);}
 ulong *xs=malloc((size_t)points*sizeof(ulong));CHECK(xs);for(long j=0;j<points;j++)xs[j]=(ulong)j;
 ulong *values[18];for(int i=0;i<18;i++){values[i]=malloc((size_t)points*sizeof(ulong));CHECK(values[i]);}
 clock_t start=clock();unsigned long long comparisons=0;
 for(long k=0;k<use;k++){
  ulong p=primes[k];CHECK(p<2147483648UL && p>(ulong)points && prime_trial(p));for(long j=0;j<k;j++)CHECK(primes[j]!=p);
  nmod_poly_t reduced;nmod_poly_init(reduced,p);
  for(int i=0;i<18;i++){fmpz_poly_get_nmod_poly(reduced,polys[i]);nmod_poly_evaluate_nmod_vec_fast(values[i],reduced,xs,points);}
  nmod_poly_clear(reduced);nmod_mat_t M;nmod_mat_init(M,N,N,p);
  for(long j=0;j<points;j++){
   uint64_t W[N][N]={{0}};
   for(int i=0;i<7;i++)for(int t=0;t<9;t++)W[i][i+t]=values[8-t][j];
   for(int i=0;i<8;i++)for(int t=0;t<8;t++)W[7+i][i+t]=values[16-t][j];
   for(int u=0;u<N;u++)for(int v=0;v<N;v++)nmod_mat_entry(M,u,v)=(ulong)W[u][v];
   ulong expected=values[17][j];CHECK(nmod_mat_det(M)==expected);CHECK(det_plain(W,p)==expected);comparisons+=2;
  }
  nmod_mat_clear(M);
  if((k+1)%16==0 || k+1==use){printf("{\"stage\":\"norm_grid\",\"primes\":%ld,\"points_per_prime\":%ld,\"comparisons\":%llu,\"cpu_seconds\":%.3f}\n",k+1,points,comparisons,(double)(clock()-start)/CLOCKS_PER_SEC);fflush(stdout);}
 }
 printf("{\"result\":\"PASS\",\"primes\":%ld,\"complete\":%s,\"comparisons\":%llu,\"flint_version\":\"%s\",\"cpu_seconds\":%.3f}\n",use,use==np?"true":"false",comparisons,FLINT_VERSION,(double)(clock()-start)/CLOCKS_PER_SEC);
 for(int i=0;i<18;i++){fmpz_poly_clear(polys[i]);free(values[i]);}free(xs);free(primes);flint_cleanup();return 0;
}
