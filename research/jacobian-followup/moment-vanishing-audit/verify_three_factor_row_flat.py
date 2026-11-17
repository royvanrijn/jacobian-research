#!/usr/bin/env python3
"""Exact PID Smith certificates for the normalized row open set."""
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
row=[S.fraction(S.cancel(v.subs({b:1,c:1,f:z}).subs(a,ar)))[0] for v in raw[3::2]]
domain=S.QQ.poly_ring(z);data=[]
for k in range(10):
 cols=[]
 for g in row:
  deg=S.Poly(g,d,e).total_degree()
  for i in range(k-deg+1):
   h=S.Poly(S.expand(g*d**i*e**(k-deg-i)),d,e)
   cols.append([h.coeff_monomial(d**j*e**(k-j)) for j in range(k+1)])
 if cols:
  mat=S.Matrix.hstack(*(S.Matrix(v) for v in cols));N=smith_normal_form(mat,domain=domain)
  diag=[S.factor(N[i,i]) for i in range(min(N.shape)) if N[i,i]!=0]
 else:mat=S.zeros(k+1,0);diag=[]
 stripped=[]
 for v in diag:
  poly=S.Poly(v,z)
  while S.rem(poly,S.Poly(t,z)).is_zero:poly=S.quo(poly,S.Poly(t,z))
  stripped.append(str(S.factor(poly.as_expr())))
 item={'binary_degree':k,'relation_matrix_shape':list(mat.shape),'free_rank':k+1-len(diag),'Smith_nonzero_diagonal':list(map(str,diag)),'after_removing_t_powers':stripped}
 data.append(item);print(json.dumps(item),flush=True)
 (p/'three_factor_row_flat_partial.json').write_text(json.dumps(data,indent=2)+'\n')
 assert all(S.degree(S.sympify(v),z)==0 for v in stripped)
assert [x['free_rank'] for x in data]==[1,2,3,4,4,4,3,2,1,0]
cert={'scope':'QQ[z,1/(1+2z^2)]-free normalized row algebra, binary degrees 0..8; degree9 zero implies all later degrees zero','data':data,'length':24,'elapsed_seconds':time.monotonic()-start,'input_checker_sha256':hashlib.sha256((p/'verify_nodal_product_reduction.py').read_bytes()).hexdigest(),'assurance':'Internal symbolic calculation; no independent replay or formal verification'}
(p/'three_factor_row_flat_certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
print('END',flush=True)
