"""Exact cofactor recovery for the degree87 P1=Q1=0 chart."""
from pathlib import Path
import json,hashlib,time,sympy as S
from flint import fmpz_poly
HERE=Path(__file__).resolve().parent;START=time.monotonic();s,a,b=S.symbols('s a b')
old=json.loads((HERE.parent/'gao-f7-linear-cube/high-root-bivariate-input.json').read_text());L=S.Poly(S.sympify(old['linear_equations'][0]),b);P,Q=map(lambda x:S.Poly(x,a,s),[L.nth(1),L.nth(0)])
def rows(F):
 return [fmpz_poly([int(F.coeff_monomial(a**i*s**j)) for j in range(F.degree(s)+1)]) for i in range(F.degree(a)+1)]
pp,qq=rows(P),rows(Q);M=[[fmpz_poly([]) for j in range(8)] for i in range(8)]
for i in range(5):
 for j,c in enumerate(reversed(pp)):M[i][i+j]=c
for i in range(3):
 for j,c in enumerate(reversed(qq)):M[5+i][i+j]=c
def det(W):
 W=[r[:] for r in W];prev=fmpz_poly([1]);sign=1
 for k in range(len(W)-1):
  z=next((i for i in range(k,len(W)) if W[i][k]),None)
  if z is None:return fmpz_poly([])
  if z!=k:W[z],W[k]=W[k],W[z];sign=-sign
  pivot=W[k][k]
  for i in range(k+1,len(W)):
   left=W[i][k]
   for j in range(k+1,len(W)):
    value=W[i][j]*pivot-left*W[k][j];u,r=divmod(value,prev);assert not r;W[i][j]=u
   W[i][k]=fmpz_poly([])
  prev=pivot
 return sign*W[-1][-1]
U=det(M);co=[]
for col in [6,7]:co.append((-1)**(5+col)*det([[c for j,c in enumerate(r) if j!=col] for i,r in enumerate(M) if i!=5]))
# Remove their common polynomial gcd to keep finite-ring arithmetic small.
g=co[0].gcd(co[1]);num,den=[divmod(v,g)[0] for v in co]
pole=json.loads((HERE/'coordinate-pole-resultant-candidate.json').read_text());H=fmpz_poly(list(map(int,next(r['coefficients'] for r in pole['factors'] if len(r['coefficients'])==88))));assert U==fmpz_poly(list(map(int,pole['coefficients'])))
rec={'status':'CANDIDATE_REQUIRES_REPLAY','source_sha256':hashlib.sha256((HERE.parent/'gao-f7-linear-cube/high-root-bivariate-input.json').read_bytes()).hexdigest(),'variable':'s','resultant_coefficients':list(map(str,U.coeffs())),'projection_coefficients':list(map(str,H.coeffs())),'cofactor_row':5,'numerator_column':6,'denominator_column':7,'raw_numerator':list(map(str,co[0].coeffs())),'raw_denominator':list(map(str,co[1].coeffs())),'common_factor':list(map(str,g.coeffs())),'a_numerator':list(map(str,num.coeffs())),'a_denominator':list(map(str,den.coeffs())),'resultant_degree_bound':118,'cofactor_degree_bound':102}
(HERE/'coordinate-poles-candidate.json').write_text(json.dumps(rec,separators=(',',':'))+'\n');print(json.dumps({'resultant_degree':U.degree(),'projection_degree':H.degree(),'raw_cofactor_degrees':[v.degree() for v in co],'reduced_degrees':[num.degree(),den.degree()],'seconds':time.monotonic()-START}))
