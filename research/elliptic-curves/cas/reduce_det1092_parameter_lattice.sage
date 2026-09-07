#!/usr/bin/env sage-python
"""Exact strict local contractions and bounded2D lattice reduction of a K3 base chart."""
import json,hashlib,sys
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,matrix,gcd
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
from research_runtime.store import checkpoint,digest
D=ROOT/'artifacts/local/elliptic-curves/det1092-parameter-lattice-reduction-v1';BASE=ROOT/'artifacts/local/elliptic-curves/det1092-integral-parameter-model-v1/result.json';SUPPORT=ROOT/'artifacts/local/elliptic-curves/det1092-parameter-hessian-support-v1/result.json'
def hashed(p):return hashlib.sha256(p.read_bytes()).hexdigest()
R=PolynomialRing(ZZ,'s');s=R.gen()
def transform(f,n,M):
 a,b,c,d=M.list();return sum(f[i]*(a*s+b)**i*(c*s+d)**(n-i) for i in range(n+1))
def valuation(f,p):return min(c.valuation(p) for c in f if c)
def divide(f,d):
 assert all(c%d==0 for c in f);return R([c//d for c in f])
def resultant(f,g):
 rows=[]
 for i in range(12):rows.append([0]*i+[f[k] for k in range(8,-1,-1)]+[0]*(11-i))
 for i in range(8):rows.append([0]*i+[g[k] for k in range(12,-1,-1)]+[0]*(7-i))
 return matrix(ZZ,rows).det()
def serialize(f,g,M,u):return dict(A_coefficients=list(map(str,f)),B_coefficients=list(map(str,g)),parameter_matrix=list(map(str,M.list())),weierstrass_u=str(u))
def main():
 p=json.loads((D/'protocol.json').read_bytes());assert p['script_sha256']==hashed(Path(__file__)) and all(hashed(ROOT/n)==h for n,h in p['inputs'].items());base=json.loads(BASE.read_bytes());support=json.loads(SUPPORT.read_bytes());F=R(base['A_coefficients']);G=R(base['B_coefficients']);F0,G0=F,G;M=matrix(ZZ,2,[1,0,0,1]);u0=u=QQ(1)/ZZ(base['scale_L']);res0=resultant(F,G);assert res0
 data=dict(status='RUNNING_LOCAL_CONTRACTIONS',protocol_hash=digest(p),steps=[],initial=serialize(F,G,M,u),initial_resultant_sha256=hashlib.sha256(str(res0).encode()).hexdigest(),initial_resultant_valuations={str(q):int(res0.valuation(q)) for q in support['candidate_primes']});checkpoint(D/'result.json',data)
 def save(kind,q,k,N):
  data['steps'].append(dict(kind=kind,prime=int(q),scale_exponent=int(k),base_matrix=list(map(str,N.list())),state=serialize(F,G,M,u)));assert len(data['steps'])<=p['maximum_steps'];checkpoint(D/'result.json',data);print(kind,q,k,'bits',max(abs(c).nbits() for c in F),max(abs(c).nbits() for c in G),flush=True)
 for q in support['candidate_primes']:
  q=ZZ(q);k=min(valuation(F,q)//4,valuation(G,q)//6)
  if k:
   F=divide(F,q**(4*k));G=divide(G,q**(6*k));u*=q**k;save('global_weierstrass_scale',q,k,matrix(ZZ,2,[1,0,0,1]))
  while True:
   vf=valuation(F,q);vg=valuation(G,q);Rp=PolynomialRing(GF(q),'x');polys=[]
   if vf<8:polys.append(Rp(divide(F,q**vf)))
   if vg<12:polys.append(Rp(divide(G,q**vg)))
   assert polys;common=polys[0]
   for h in polys[1:]:common=common.gcd(h)
   roots=sorted(int(r) for r in common.roots(multiplicities=False));candidates=[matrix(ZZ,2,[q,r,0,1]) for r in roots]+[matrix(ZZ,2,[0,1,q,0])];good=[]
   for N in candidates:
    f=transform(F,8,N);g=transform(G,12,N);k=min(valuation(f,q)//4,valuation(g,q)//6)
    if k>=2:good.append((N,f,g,k))
   assert len(good)<=1
   if not good:break
   N,f,g,k=good[0];F=divide(f,q**(4*k));G=divide(g,q**(6*k));M=M*N;u*=q**k
   content=gcd(M.list());M=M/content;u/=content**2;M=matrix(ZZ,M)
   save('strict_base_contraction',q,k,N)
 res=resultant(F,G);assert res==QQ(res0)*M.det()**96*(u0/u)**96 and res
 data.update(status='RUNNING_ARCHIMEDEAN_CHARTS',contracted=serialize(F,G,M,u),contracted_resultant_sha256=hashlib.sha256(str(res).encode()).hexdigest(),charts=[]);checkpoint(D/'result.json',data)
 for exponent in p['lattice_weight_exponents']:
  w=matrix(ZZ,2,[2**max(exponent,0),0,0,2**max(-exponent,0)]);B=M.transpose()*w;reduced=B.LLL();U=reduced*B.inverse();U=matrix(ZZ,U);assert abs(U.det())==1;N=U.transpose();f=transform(F,8,N);g=transform(G,12,N);newM=M*N
  height=max(max(abs(c) for c in f)**3,max(abs(c) for c in g)**2)
  data['charts'].append(dict(weight_exponent=exponent,height=str(height),height_bits=int(height.nbits()),state=serialize(f,g,newM,u)));checkpoint(D/'result.json',data)
 chosen=min(data['charts'],key=lambda r:(ZZ(r['height']),r['state']['parameter_matrix']));data.update(status='PASS',selected=chosen,boundary='Strict determinant-decreasing contractions at the frozen necessary-prime support, followed by17 weighted2D LLL bases. Exact source-family identity and resultant scaling checked. No global arithmetic-model minimality or coefficient optimum claim; equal-resultant intermediate routes and other fibrations are outside this protocol. No point search or parameter-population scan.');checkpoint(D/'result.json',data)
 print('SELECTED bits',chosen['height_bits'],'M',chosen['state']['parameter_matrix'],'u',chosen['state']['weierstrass_u'],flush=True)
if __name__=='__main__':main()
