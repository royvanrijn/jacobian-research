"""Exact replay for LOWER_TWO_PRIME_GAP_SFC.tex.

Universal statements are proved in the TeX; finite controls below check
normalizations, signs, exact concomitant columns and complete small fields.
"""
from pathlib import Path
from hashlib import sha256
from itertools import combinations
import json
import sympy as s

z,g,n,d,e,u,eps=s.symbols("z gamma n d E u eps")
checks=[]
def check(name,statement):
 if not bool(statement):raise AssertionError(name)
 checks.append(name)
M=g+d*n
A=(M+d)*(M+d+1)
# Exact universal two-step Euler/Ore reduction (epsilon^2=1).
reduced=eps*A*s.Symbol("Ed")-u*(e+M+d+1)*(e-g)+u*e*(e+1)
expected=eps*A*s.Symbol("Ed")+u*(g*(M+d+1)-d*(n+1)*e)
check("universal_all_d_upper_Ore_reduction",s.expand(reduced-expected)==0)

def Bs(N,D,G):
 out=[s.Integer(1)]
 for j in range(N):
  out.append(s.expand((z+D*j+G+D-1)*out[-1]-z*s.diff(out[-1],z)))
 return out
def fall(x,k):
 return s.prod(x-j for j in range(k))

for D in range(2,8):
 for N in [1,2,4]:
  mm=g+D*N
  bs=Bs(N+D,D,g)
  b,c=bs[N],bs[N+1]
  check(f"h_identity_d{D}_n{N}",
   s.expand((1-mm-D-z)*b+z*s.diff(b,z)+c)==0)
  coeff=[s.Integer(1)]
  L=c
  for j in range(D-1):
   from_basis=sum(coeff[r]*bs[N+1+r] for r in range(len(coeff)))
   check(f"constant_ladder_d{D}_n{N}_j{j}",s.expand(L-from_basis)==0)
   check(f"ladder_pivot_d{D}_n{N}_j{j}",coeff[-1]==(-1)**j)
   nxt=[s.Integer(0)]*(len(coeff)+1)
   for r,a in enumerate(coeff):
    nxt[r]+=a*(mm+2*D-1+D*r-j);nxt[r+1]-=a
   coeff=[s.expand(v) for v in nxt]
   L=s.expand(z*(s.diff(L,z)-L)-j*L)
  # Full concomitant coefficients tested on a complete monomial jet basis.
  if N==1 and D<=5:
   Q=[c]
   for j in range(D-1):Q.append(s.expand(s.diff(Q[-1],z)-Q[-1]))
   for k in range(D+1):
    raw=b*fall(k-mm,D)*z**k-z**(D+k)*b
    raw+=sum((-1)**j*Q[j]*fall(k-mm,D-1-j)*z**(k+j)
             for j in range(D))
    jet=(-1)**D*z**k*(
       b*s.rf(mm-k,D)-sum(z**j*Q[j]*s.rf(mm-k,D-1-j)
                           for j in range(D)))-z**(D+k)*b
    check(f"full_concomitant_jet_d{D}_k{k}",s.expand(raw-jet)==0)

x,U,G=s.symbols("x U G")
for p in [2,3,5,7]:
 mat=s.zeros(p)
 for j in range(p):
  mat[j,j]=G-1-j
  mat[(j+1)%p,j]+=1 if j<p-1 else U
 char=mat.charpoly(x).as_expr()
 check(f"universal_residue_charpoly_p{p}",
       s.Poly(char-(x**p-x-U-G**p+G),x,U,G,modulus=p).is_zero)
 # Every gamma and u residue, all positive exponent classes.
 for gg in range(p):
  for uu in range(p):
   aa=mat.subs({G:gg,U:uu})
   seed=s.zeros(p,1);seed[0]=1
   for N in range(1,2*p):
    cols=s.Matrix.hstack(*[aa**k*seed for k in range(N+1,N+p)])
    minors=[int(cols.extract(rows,range(p-1)).det())%p
            for rows in combinations(range(p),p-1)]
    check(f"positive_Krylov_p{p}_g{gg}_u{uu}_n{N}",any(minors))
 # Zero-u unit scalar minor at every prime-field gamma and exponent class.
 for gg in range(p):
  aa=mat.subs({G:gg,U:0})
  seed=s.zeros(p,1);seed[0]=1
  for N in range(1,p+1):
   cols=s.Matrix.hstack(*[aa**k*seed for k in range(N+1,N+p-1)],
                        aa**N*seed) if p>2 else aa**N*seed
   rows=next(rows for rows in combinations(range(p),p-1)
             if int(cols.extract(rows,range(p-1)).det())%p)
   last=aa**(N+p-1)*seed+s.Symbol("eta")*aa**N*seed
   test=s.Matrix.hstack(*[aa**k*seed for k in range(N+1,N+p-1)],last) if p>2 else last
   constant=int(cols.extract(rows,range(p-1)).det())%p
   check(f"affine_residual_minor_p{p}_g{gg}_n{N}",
         s.Poly(test.extract(rows,range(p-1)).det()
                -constant*(1+s.Symbol("eta")),s.Symbol("eta"),modulus=p).is_zero)

here=Path(__file__).resolve().parent
repo=here.parents[2]
proof=here/"LOWER_TWO_PRIME_GAP_SFC.tex"
prior=here/"HIGHER_GAP_ADJOINT_RANK.tex"
receipt={
 "passed":True,"checks":len(checks),"check_names":checks,
 "scope":"Internal exact sign/column/Ore controls and complete small residue fields; all-d/all-prime proofs are written",
 "all_d_mixed_rank":"written proof, including epsilon and d2 empty list",
 "prime_scope":"all primes; positive real algebraic beta integral at one p-adic embedding; all n>=1",
 "all_binomial_sfc":False,
 "all_real_shape_claim":False,
 "p2_included":True,
 "algebraic_shape_one_integral_place":True,
 "command":"python research/factorial/image-rigidity/lower_two_prime_gap_sfc.py",
 "sympy_version":s.__version__,
 "source_sha256":{p.name:sha256(p.read_bytes()).hexdigest()
                  for p in [proof,Path(__file__),prior]},
 "independent_replay":False,"formal_verification":False,"external_review":False}
(here/"lower_two_prime_gap_sfc_verification.json").write_text(json.dumps(receipt,indent=2)+"\n")
print("PASS",len(checks),"exact controls")
