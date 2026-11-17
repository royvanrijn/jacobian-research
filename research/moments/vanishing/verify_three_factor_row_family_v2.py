#!/usr/bin/env python3
"""Exact PID Smith certificates for the full b=c=1 normalized row family."""
from pathlib import Path
import json,time,hashlib,resource
import sympy as S
from sympy.matrices.normalforms import smith_normal_form
from verify_nodal_product_reduction import moments
resource.setrlimit(resource.RLIMIT_AS,(768*1024*1024,768*1024*1024))
resource.setrlimit(resource.RLIMIT_CPU,(35,35))
p=Path(__file__).resolve().parent;start=time.monotonic()
a,b,c,d,e,f,z=S.symbols('a b c d e f z');t=1+2*z*z
raw=moments([[a,d,e],[d,b,f],[e,f,c]],10)
ar=-2*(d*d+e*e+4*z*d*e)/t
row=[S.expand(v.subs({b:1,c:1,f:z})) for v in raw[1::2]]
domain=S.QQ.poly_ring(z);data=[]
for k in range(10):
 base=[(i,j,k-2*i-j) for i in range(k//2+1) for j in range(k-2*i+1)]
 cols=[]
 for n,g in enumerate(row,1):
  deg=2*n
  for i in range((k-deg)//2+1):
   for j in range(k-deg-2*i+1):
    h=S.Poly(S.expand(g*a**i*d**j*e**(k-deg-2*i-j)),a,d,e)
    cols.append([h.coeff_monomial(a**u*d**v*e**w) for u,v,w in base])
 if cols:
  mat=S.Matrix.hstack(*(S.Matrix(v) for v in cols));N=smith_normal_form(mat,domain=domain)
  diag=[S.factor(N[i,i]) for i in range(min(N.shape)) if N[i,i]!=0]
 else:mat=S.zeros(len(base),0);diag=[]

 item={'weighted_degree':k,'relation_matrix_shape':list(mat.shape),'free_rank':len(base)-len(diag),'Smith_nonzero_diagonal':list(map(str,diag))}
 data.append(item);print(json.dumps(item),flush=True)
 (p/'three_factor_row_family_v2_partial.json').write_text(json.dumps(data,indent=2)+'\n')
 assert all(S.degree(v,z)==0 for v in diag)
assert [x['free_rank'] for x in data]==[1,2,3,4,4,4,3,2,1,0]
h=S.Poly(row[1],a,d,e).coeff_monomial(a**2)
u,v,gcd=S.gcdex(t,h,z)
assert gcd==1 and S.expand(u*t+v*h)==1
cert={'degree10_Bezout':{'M2_a_coefficient':str(t),'M4_a2_coefficient':str(h),'u':str(u),'v':str(v),'identity':'u*t+v*h=1; degree9 zero kills every degree10 monomial containing d or e, so M2*a^4 and M4*a^3 kill a^5'},'scope':'QQ[z]-free normalized row algebra, weighted degrees 0..8 (a=2,d=e=1); degree9 and degree10 zero imply all later degrees zero','data':data,'length':24,'elapsed_seconds':time.monotonic()-start,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'terminal':True,'completed':True,'input_checker_sha256':hashlib.sha256((p/'verify_nodal_product_reduction.py').read_bytes()).hexdigest(),'assurance':'Internal symbolic calculation; no independent replay or formal verification'}
(p/'three_factor_row_family_v2_certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
print('END',flush=True)
