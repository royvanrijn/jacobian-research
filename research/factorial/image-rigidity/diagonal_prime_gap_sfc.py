"""Exact controls for DIAGONAL_PRIME_GAP_SFC.tex.

All-offset and all-prime/all-index claims rely on the written argument.
These controls check its quotient shifts, integral basis, normalizations,
prime unit cancellation and finite-field/literal identities.
"""
from pathlib import Path
from hashlib import sha256
from itertools import combinations
import json
import platform
import sympy as s

E,C,M,d,n,g,u,eps=s.symbols("E C M d n gamma u epsilon")
checks=[]
def check(name,statement):
 if not bool(statement):raise AssertionError(name)
 checks.append(name)
def rising_basis(poly,degree):
 out=[s.Integer(0)]*(degree+1)
 rem=s.expand(poly)
 for k in range(degree,-1,-1):
  out[k]=s.expand(rem).coeff(E,k)
  rem=s.expand(rem-out[k]*s.rf(E,k))
 assert rem==0
 return out

for h in range(1,10):
 V,rem=s.div(s.Poly(s.rf(E,h)-s.rf(C,h),E),s.Poly(E-C,E))
 V=V.as_expr()
 check(f"exact_divided_difference_h{h}",rem.is_zero)
 vs=rising_basis(V,h-1)
 check(f"monic_rising_basis_h{h}",vs[-1]==1)
 for k,v in enumerate(vs):
  check(f"integral_rising_coefficient_h{h}_k{k}",s.Poly(v,C).domain==s.ZZ)
  check(f"closed_telescoping_coefficient_h{h}_k{k}",
        s.expand(v-s.rf(C+k+1,h-1-k))==0)
 # This checks the exact G(E-d)=V(E) commutator shift universally in d,M.
 shifted=V.subs({E:E+d,C:M+d},simultaneous=True)
 check(f"Ore_quotient_shift_h{h}",
       s.expand((E-M)*shifted-(s.rf(E+d,h)-s.rf(M+d,h)))==0)
 vv=V.subs(C,g+d*n+d)
 aa=s.rf(g+d*n+d,h)
 raw=u*s.rf(E,h)+eps*aa*s.Symbol("K")-u*vv*(E-g)
 target=eps*aa*s.Symbol("K")+u*(aa-d*(n+1)*vv)
 check(f"universal_upper_reduction_h{h}",s.expand(raw-target)==0)
 # A deliberately changed quotient must fail the same complete identity.
 wrong=target-u*d*(n+1)
 check(f"corrupted_quotient_rejection_h{h}",s.expand(raw-wrong)!=0)

v2=s.cancel((s.rf(E,2)-s.rf(C,2))/(E-C))
check("offset_two_normalization",s.expand(v2-E-C-1)==0)
v3=s.cancel((s.rf(E,3)-s.rf(M+3,3))/(E-M-3))
check("cubic_rising_basis_normalization",
      s.expand(v3-s.rf(E,2)-(M+5)*E-(M+4)*(M+5))==0)
B0,B1,B2,K=s.symbols("B_n B_next B_next2 kappa_u")
col0=B1-K*B0
col1=(M+5)*B1-B2-K*(M+5)*B0
check("cubic_two_column_elimination",s.expand(col1-(M+5)*col0+B2)==0)
kap=s.Rational(3)*(n+1)/((6*n+4)*(6*n+5)*(6*n+6))
check("cubic_scalar_unit_normalization",
      s.cancel(kap-1/(4*(3*n+2)*(6*n+5)))==0)

x=s.symbols("x")
for p in [3,5,7,11,13]:
 cc=2*p*n+p+1
 aa=s.rf(cc,p)
 kap=p*(n+1)/aa
 unit=1/(2*s.prod(cc+j for j in range(p-1)))
 check(f"exact_unit_cancellation_p{p}",s.cancel(kap-unit)==0)
 check(f"all_other_denominator_factors_are_units_p{p}",
       all((j+1)%p for j in range(p-1)))
for p in [2,3,5,7]:
 mat=s.zeros(p)
 for j in range(p):
  mat[j,j]=-j
  if j+1<p:mat[j+1,j]=1
 check(f"complete_characteristic_polynomial_p{p}",
       s.Poly(mat.charpoly(x).as_expr()-(x**p-x),x,modulus=p).is_zero)
 seed=s.zeros(p,1);seed[0]=1
 basis=s.Matrix.hstack(*[mat**j*seed for j in range(p)])
 check(f"monic_cyclic_basis_p{p}",int(basis.det())%p==1)
 for N in range(1,2*p):
  cols=s.Matrix.hstack(*[mat**k*seed for k in range(N+1,N+p)])
  # The fixed rows1..p-1 give a unit; no search or rank heuristic.
  check(f"fixed_positive_Krylov_minor_p{p}_n{N}",
        int(cols.extract(range(1,p),range(p-1)).det())%p!=0)
  L=[]
  current=mat**(N+1)*seed
  for j in range(p-1):
   L.append(current)
   current=(-mat-j*s.eye(p))*current
  lj=s.Matrix.hstack(*L)
  check(f"fixed_ladder_minor_p{p}_n{N}",
        int(lj.extract(range(1,p),range(p-1)).det())%p!=0)

t=s.symbols("t")
def T(N,G,D):
 return s.Poly(sum(s.binomial(N,k)*s.rf(G,D*k)*t**k
                   for k in range(N+1)),t)
for p,N in [(2,1),(2,3),(3,1),(3,2),(3,4),(5,1),(5,2),(7,1)]:
 lo=T(N,p*N+1,p);hi=T(N+1,p*(N+1)+1,p)
 check(f"literal_ordinary_diagonal_gcd_p{p}_n{N}",s.gcd(lo,hi).degree()==0)

here=Path(__file__).resolve().parent
inputs=[here/"DIAGONAL_PRIME_GAP_SFC.tex",Path(__file__),
        here/"LOWER_TWO_PRIME_GAP_SFC.tex"]
receipt={
 "passed":True,"checks":len(checks),"check_names":checks,
 "scope":"Exact internal quotient/Ore/basis/prime-residue controls; unbounded claims have written proofs",
 "general_offset_scope":"all d>=2,n>=1,gamma>0,1<=h<=d; exact mixed-rank equivalence",
 "diagonal_scope":"ordinary X^p(A+BX^p), every prime and every consecutive positive moment pair",
 "cubic_support":"ordinary support{3,6}, all indices",
 "all_real_gamma_shape_diagonal_claim":False,
 "all_binomial_sfc":False,
 "exact_imports":{
   "general_concomitant":"LOWER_TWO_PRIME_GAP_SFC.tex thm:mixed proof equation(3)",
   "prime_two_branch":"LOWER_TWO_PRIME_GAP_SFC.tex thm:prime at beta1"},
 "source_sha256":{p.name:sha256(p.read_bytes()).hexdigest() for p in inputs},
 "python":platform.python_version(),"sympy":s.__version__,
 "command":"python research/factorial/image-rigidity/diagonal_prime_gap_sfc.py",
 "independent_replay":False,"formal_verification":False,"external_review":False}
(here/"diagonal_prime_gap_sfc_verification.json").write_text(json.dumps(receipt,indent=2)+"\n")
print("PASS",len(checks),"exact controls")
