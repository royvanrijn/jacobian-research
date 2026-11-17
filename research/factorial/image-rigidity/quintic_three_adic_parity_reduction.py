from pathlib import Path
from math import comb,prod
import hashlib,json
import sympy as s
from flint import fmpz_poly
HERE=Path(__file__).resolve().parent
IMPORTED=HERE/"lower_exponent_initial_and_boundary.tex"
IMPORTED_SHA256="e66e02d0354de5cd92d33eef54efea15c387b9f6a0c1b62ed379295e13f07d45"
assert hashlib.sha256(IMPORTED.read_bytes()).hexdigest()==IMPORTED_SHA256
assert all(label in IMPORTED.read_text() for label in (
 "thm:gap-three-offset-two-rank","cor:gap-three-lower-two-eight-ninths",
 "cor:ordinary-quintic-ray-reduction","thm:p3-index-interpolation-obstruction"))
def vp(n,p=3):
 assert n
 out=0
 while n%p==0:n//=p;out+=1
 return out
def vf(n,p=3):
 out=0
 while n:n//=p;out+=n
 return out
def coeffs(n):
 c=[1]
 for k in range(1,n+1):
  c.append(c[-1]*(n-k+1)*prod(range(2*n+3*k-2,2*n+3*k+1))//k)
 return c
counts={}
# Complete formal operator calculation, with independent parameters n and M.
z,n,M,xi=s.symbols("z n M xi")
def add_op(P,Q):
 out=dict(P)
 for j,a in Q.items():out[j]=s.expand(out.get(j,0)+a)
 return {j:s.factor(a) for j,a in out.items() if a!=0}
def compose_op(P,Q):
 out={}
 for i,a in P.items():
  for j,b in Q.items():
   for k in range(i+1):
    deg=i-k+j
    out[deg]=out.get(deg,0)+s.binomial(i,k)*a*s.diff(b,z,k)
 return {j:s.factor(a) for j,a in out.items() if a!=0}
def rising_E(order):
 out={0:s.Integer(1)}
 for j in range(order):out=compose_op({0:M+j,1:-z},out)
 return out
L={0:3*n/z,1:-(M*(M-1)*(M+1)/z**3+1),
   2:3*M*(M-1)/z**2,3:-3*(M-1)/z,4:s.Integer(1)}
O=add_op(rising_E(5),{j:z**3*a for j,a in rising_E(2).items()})
while max(O)>=4:
 j=max(O);a=O[j]
 correction=compose_op({j-4:a},L)
 O=add_op(O,{k:-v for k,v in correction.items()})
A=(M+3)*(M+4)
expected={0:s.rf(M,5)+(M*(M+1)-6*n*(M+2))*z**3,
 1:3*z*(-M*(M+1)*A+(n+1)*z**3),
 2:3*M*z**2*A,3:-z**3*A}
assert all(s.cancel(O.get(j,0)-expected[j])==0 for j in range(4))
counts["complete_symbolic_Ore_remainder_coefficients"]=4
# Formal B,B',B'' jet at x=xi*z; no sampled auxiliary polynomial.
B0,B1,B2=s.symbols("B0 B1 B2")
logg=(1-M)/z-xi
Fp=xi*B1+logg*B0
Fpp=xi**2*B2+2*logg*xi*B1+(s.diff(logg,z)+logg**2)*B0
a3=-3*(M-1)/z;a2=3*M*(M-1)/z**2
b,c,d=expected[1],expected[2],expected[3]
U=d*(a3*B0-Fp)-c*B0
V=d*(Fpp-s.diff(a3,z)*B0-a3*Fp+a2*B0)-b*B0
x=xi*z
C=(x+M+2)*B0-x*B1
Cx=B0+(x+M+1)*B1-x*B2
D=(x+M+5)*C-x*Cx
H=D+3*(n+1)/A*z**3*B0
for expr in [U+z**2*A*C,V+(3*M+5)/z*U+z*A*H]:
 numerator=s.cancel(expr).as_numer_denom()[0]
 assert s.rem(s.Poly(numerator,xi),s.Poly(xi**3-1,xi)).is_zero
counts["complete_symbolic_gauge_stripped_column_identities"]=2
# Reduction has only3 scalar gamma residues and2 positive-power parity classes.
z=s.symbols("z")
for gamma in range(3):
 def A(f):return s.Poly(z*f-z*s.diff(f,z)+(gamma-1)*f,z, modulus=3).rem(s.Poly(z**3,z,modulus=3)).as_expr()
 powers=[s.Integer(1)]
 for j in range(6):powers.append(A(powers[-1]))
 assert all(s.Poly(powers[j+3]-powers[j+1],z,modulus=3).is_zero for j in range(4))
 def col(f):return s.Matrix([s.Poly(f,z).nth(j)%3 for j in range(3)])
 K=s.Matrix.hstack(*map(col,powers[:3]))
 assert int(K.det())%3==1
 for n in (1,2):
  B,C,D=map(col,powers[n:n+3])
  assert any(int(s.Matrix.hstack(C,B).extract(rows,[0,1]).det())%3 for rows in [(0,1),(0,2),(1,2)])
  for eta in range(3):
   H=D+eta*B
   minors=[int(s.Matrix.hstack(C,H).extract(rows,[0,1]).det())%3 for rows in [(0,1),(0,2),(1,2)]]
   assert (not any(minors))==(eta==2)
counts["complete_F3_gamma_and_positive_power_parity_cases"]=6
counts["complete_F3_eta_rank_cases"]=18
indices=[1,19,28,55]
literal=0;gcd_cases=0
for n in indices:
 c=coeffs(n)
 assert vp(c[1])==1
 assert all(vp(c[k])>k for k in range(2,n+1))
 g=(2*n+1)//3
 for k,a in enumerate(c):
  direct=comb(n,k)*prod(range(2*n+1,2*n+3*k+1))
  assert a==direct
  h=vf(n)-vf(k)-vf(n-k)+vf(2*n+3*k)-vf(2*n)
  assert vp(a)==h
  assert h==k+vp(comb(n,k)*prod(range(g,g+k)))
  literal+=1
 if vp(5*n+4)==2:
  assert fmpz_poly(c).gcd(fmpz_poly(coeffs(n+1))).degree()==0
  gcd_cases+=1
counts["literal_coefficient_valuation_and_block_identities"]=literal
counts["literal_gcd_controls_in_proved_congruence_classes"]=gcd_cases
assert [pow(2,12*k+6,27) for k in range(3)]==[10,1,19]
counts["full_sparse_family_exponent_period_cases"]=3
# The stronger width-one conclusion requires r2; n10 is outside it.
assert vp(5*10+4)==3
counts["outside_width_one_balance_hypothesis_control"]=1
out={"passed":True,
 "scope":"Exact internal symbolic/reduced-field replay and bounded literal controls; universal claims rely on the retained proof",
 "universal_statement":"All ordinary a2d3 pairs at lower index n not10mod27 are coprime; all degree<=5 ordinary binomials reduce to terminal N11mod27",
 "interpolation_consequence":"The retained tau_star and m_star are not nonnegative ordinary integers; no irrationality or transcendence assertion",
 "open_terminal_ray":"N11mod27 remains open",
 "controls":counts,
 "source_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [
  HERE/"QUINTIC_THREE_ADIC_PARITY_REDUCTION.tex",Path(__file__),IMPORTED]},
 "exact_imports":{
  "mixed_adjoint":"lower_exponent_initial_and_boundary.tex thm:gap-three-offset-two-rank and algebraic transfer in cor:gap-three-lower-two-eight-ninths",
  "degree_five_classification":"lower_exponent_initial_and_boundary.tex cor:ordinary-quintic-ray-reduction",
  "interpolated_common_point":"lower_exponent_initial_and_boundary.tex thm:p3-index-interpolation-obstruction"},
 "independent_replay":False,"formal_verification":False,"external_review":False}
(HERE/"quintic_three_adic_parity_reduction_verification.json").write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps(out,indent=2))
