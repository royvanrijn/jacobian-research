"""Exact internal replay of the all-positive-Gamma lower-five sector proof.

The written spectral, GWS, Pick and three-region argument proves the
unbounded claims. These controls certify formulas and normalizations.
"""
from pathlib import Path
from itertools import combinations
import hashlib,json,platform
import sympy as S
HERE=Path(__file__).resolve().parent
PROOF=HERE/"LOWER_FIVE_QUADRATIC_SECTOR_SFC.tex"
IMPORTS={
 "LOWER_TWO_PRIME_GAP_SFC.tex":
 "8453b5b0ccd220573ef106fa26790e1f691885970b8c2116f036ed893ab2c200",
 "CUBIC_ADJOINT_STABILITY_SFC.tex":
 "ea2510e511545b9dd7a49a0796d75c45575687e60770c3359d38c5bc99cb4125",
 "sources/wagner-multivariate-stable-polynomials-2011.pdf":
 "ab72551b3407480e08185ce2b0a9527c8cea385b95ef93dcf904ec90a46170c1",
}
for name,h in IMPORTS.items():
 assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==h,(name,"changed exact input")
counts={}
def count(k,n=1):counts[k]=counts.get(k,0)+n
M,n,u=S.symbols('M n u')
gamma=M-2*n
def E(v):
 r,e,k=v
 return [-2*u*S.diff(r,u)-gamma*u*k,
         r-e-2*u*S.diff(e,u)+u*k,
         e+M*k-2*u*S.diff(k,u)]
v=[S.Integer(1),S.Integer(0),S.Integer(0)];jets={0:v}
for j in range(7):
 ev=E(v);v=[S.factor(ev[k]+j*v[k]) for k in range(3)];jets[j+1]=v
upper=[S.expand(u*jets[5][k]-jets[7][k]) for k in range(3)]
ss=2*(n+1);AA=S.rf(M+2,5);bb=3*ss*(M+4)**2
kk=(M+4)*(M+5)*(M+6);cc=M+2*n+8
UU=AA+bb*u+ss*u*u;WW=ss*u*(kk+cc*u)
assert S.expand(upper[1]+WW)==0
assert S.expand(upper[2]+UU)==0
count("generic_shift_five_rising_remainder_coefficients",2)
assert S.expand(upper[2]+(AA+bb*u+2*ss*u*u))!=0
count("corrupted_shift_coefficient_rejections")
# The generic concomitant/upper two-by-two elimination.
B,C,U,W=S.symbols('B C U W')
assert S.expand(S.Matrix([[W,U],[-C,B]]).det()-(U*C+W*B))==0
count("generic_fulljet_elimination_identities")
# Exact pairing and finite polarization normalization.
g,k,z,a=S.symbols('gamma k z a')
assert S.expand((g+2*k)*(g+2*k+1)-4*(g/2+k)*((g+1)/2+k))==0
count("generic_rising_factorial_pairing_step")
def elementary(values,j):
 return sum((S.prod(c) for c in combinations(values,j)),S.Integer(0))
for index in range(1,6):
 ds=S.symbols('d0:'+str(index));prod=S.prod(ds)
 composed=sum(S.rf(a,j)*z**(index-j)*elementary(ds,j) for j in range(index+1))
 polarized=prod*sum(S.rf(a,j)*elementary([z/d for d in ds],index-j) for j in range(index+1))
 assert S.cancel(polarized-composed)==0
 count("literal_complete_polarization_composition_identities")
 for j in range(index+1):
  assert elementary([z]*index,index-j)==S.binomial(index,j)*z**(index-j)
  count("literal_polarization_diagonal_coefficients")
 Cn=sum(S.binomial(index,j)*S.rf(g/2,j)*S.rf((g+1)/2,j)*z**(index-j) for j in range(index+1))
 t=S.symbols('t')
 Tn=sum(S.binomial(index,j)*S.rf(g,2*j)*t**j for j in range(index+1))
 assert S.expand((4*t)**index*Cn.subs(z,1/(4*t))-Tn)==0
 count("literal_moment_reversal_identities")
# A rising-coefficient nonmonic orientation would fail already at n1.
assert S.expand((1+a*z)-(z+a))!=0
count("wrong_Gamma_orientation_rejections")
# The actual tridiagonal characteristic recurrence, without numerical roots.
alpha=S.symbols('alpha')
for index in range(1,7):
 D0,D1=S.Integer(1),z+alpha+index-1
 for j in range(2,index+1):
  diag=alpha+index+1-2*j
  D0,D1=D1,S.expand((z+diag)*D1+(j-1)*(alpha+index-j+1)*D0)
 expected=sum(S.binomial(index,j)*S.rf(alpha,j)*z**(index-j) for j in range(index+1))
 d=S.Poly(S.expand(D1-expected),z)
 for j in range(index+1):
  assert d.nth(j)==0;count("literal_Gamma_matrix_characteristic_coefficients")
j=S.symbols('j')
assert S.expand((n-1)*(alpha+1)-j*(alpha+n-j)-(n-1-j)*(alpha+1-j))==0
count("generic_largest_offdiagonal_weight_identity")
R=lambda q:4*(n-1)*(q+1)/(q-n+1)**2
assert S.simplify(S.diff(R(alpha),alpha)+4*(n-1)*(alpha+n+1)/(alpha-n+1)**3)==0
count("generic_spectral_ratio_derivative")
margins=[
 (S.Rational(40,9)-R(5*n/2),(696*n+304)/(9*(3*n+2)**2)),
 (3-R(3*n),(20*n+7)/(2*n+1)**2),
 (1-R(6*n),(n*n+30*n+5)/(5*n+1)**2)]
for lhs,rhs in margins:
 assert S.cancel(lhs-rhs)==0;count("generic_three_region_spectral_margins")
 for coeff in S.Poly(S.fraction(S.factor(rhs))[0],n).all_coeffs():
  assert coeff>0;count("positive_spectral_margin_coefficients")
# Generic Pick term and imaginary rational numerator.
aa,rr,uu=S.symbols('a r u')
assert S.cancel(-rr/(uu+aa)+rr/aa-2*uu/(uu+aa)-(rr-2*aa)*uu/(aa*(uu+aa)))==0
count("generic_Pick_partial_fraction_identity")
x,y,A,b,s,c,K=S.symbols('x y A b s c K',real=True)
v=x+S.I*y;vc=x-S.I*y
Uv=A+b*v+s*v*v;Uvc=A+b*vc+s*vc*vc
Wv=s*v*(K+c*v);Wvc=s*vc*(K+c*vc)
imag=S.expand((Wv*Uvc-Wvc*Uv)/(2*S.I))
target=s*y*(K*A+2*c*A*x+(c*b-s*K)*(x*x+y*y))
assert S.expand(imag-target)==0
count("generic_imaginary_numerator_identity")
assert S.cancel((K*A+2*c*A*(-x/(x*x+y*y))+(c*b-s*K)/(x*x+y*y))*(x*x+y*y)
                -(K*A*(x*x+y*y)-2*c*A*x+c*b-s*K))==0
count("generic_negative_reciprocal_disk_identity")
F=2*M*M+(6*n+25)*M+24*n+66
delta=S.expand(cc*bb-ss*kk)
assert S.expand(delta-ss*(M+4)*F)==0
count("generic_positive_delta_factorization")
qn=ss*(M+4)*F;qd=(M+2)*(M+3)*cc**2
assert S.cancel(kk*delta/(cc**2*AA)-qn/qd)==0
count("generic_angle_ratio_cancellation")
assert S.expand(F.subs(M,-3)-(6*n+9))==0
assert S.expand(F.subs(M,-(2*n+8))-(-4*n*n-10*n-6))==0
assert S.expand((6*n+25)**2-8*(24*n+66)-(36*n*n+108*n+97))==0
count("generic_factor_root_order_controls",3)
endpoint=[
 (7*qn.subs(M,8*n)-3*qd.subs(M,8*n),
  8*(64*n**4+1492*n**3+2699*n*n+1634*n+318)),
 (4*qn.subs(M,14*n)-qd.subs(M,14*n),
  32*(98*n**4+1323*n**3+1390*n*n+483*n+54))]
for lhs,rhs in endpoint:
 assert S.expand(lhs-rhs)==0;count("generic_angle_endpoint_identities")
 for coefficient in S.Poly(rhs,n).all_coeffs():
  assert coefficient>0;count("positive_angle_endpoint_coefficients")
assert S.Rational(3,7)-S.Rational(31,49)**2==S.Rational(68,2401)
assert (S.Rational(40,9)-1)/(S.Rational(40,9)+1)==S.Rational(31,49)
count("exact_low_region_angular_comparisons",2)
# Literal Gamma parameters cover all regions/boundaries and nonalgebraic
# hypotheses are supplied by the written proof, not these rational controls.
profiles=[]
for index,beta in [(1,S.Rational(1,2)),(2,S.Rational(1,7)),(3,3),(4,28),(5,100),(6,S.Rational(5,2))]:
 gamma=5*index+beta;MM=gamma+2*index
 for order in [index,index+1]:
  for j in range(order+1):
   actual=S.binomial(order,j)*S.rf(beta,5*order+2*j)
   normalized=S.rf(beta,5*order)*S.binomial(order,j)*S.rf(5*order+beta,2*j)
   assert actual==normalized
   count("literal_Gamma_phase_normalization_coefficients")
 qvalue=S.cancel((qn/qd).subs({n:index,M:MM}))
 if beta<=index:
  assert qvalue>S.Rational(3,7)
 elif beta<=7*index:
  assert qvalue>S.Rational(1,4)
 t=S.symbols('t')
 lower=S.Poly(sum(S.binomial(index,j)*S.rf(gamma,2*j)*t**j for j in range(index+1)),t)
 upper=S.Poly(sum(S.binomial(index+1,j)*S.rf(gamma+5,2*j)*t**j for j in range(index+2)),t)
 assert S.gcd(lower,upper).degree()==0
 count("literal_positive_rational_Gamma_pair_gcd_controls")
 profiles.append({"n":index,"beta":str(beta),"exact_phase_coefficients":True,"exact_pair_coprime":True,
                  "angle_ratio":str(qvalue)})
receipt={
 "passed":True,"total_checks":sum(counts.values()),"controls":counts,
 "scope":{
  "all_positive_real_Gamma_beta":True,
  "lower_five_quadratic_gap_all_indices":True,
  "ordinary_support_five_seven_all_indices":True,
  "lower_moment_uniform_Hurwitz_assumption":False,
  "algebraic_or_p_integral_shape_restriction":False,
  "all_degree_seven_binomial_SFC":False,"unrestricted_FC_or_SFC":False,
 },
 "literal_Gamma_profiles":profiles,
 "source_sha256":{
  PROOF.name:hashlib.sha256(PROOF.read_bytes()).hexdigest(),
  Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),**IMPORTS,
 },
 "exact_imports":{
  "general_fulljet_concomitant":"LOWER_TWO_PRIME_GAP_SFC.tex thm:mixed proof Eq3 atd2",
  "Gamma_matrix_and_auxiliary_bound":"CUBIC_ADJOINT_STABILITY_SFC.tex lem:gamma; Eqs3-4",
  "GWS":"Wagner Theorem4.1, PDFp11 Section4; symmetric multiaffine polarization on an open convex half-plane",
 },
 "python":platform.python_version(),"sympy":S.__version__,
 "independent_replay":False,"formal_verification":False,"external_review":False,
 "novelty_claim":False,"priority_claim":False,
 "command":".venv/bin/python research/factorial/image-rigidity/lower_five_quadratic_sector_sfc.py",
}
out=HERE/"lower_five_quadratic_sector_sfc_verification.json"
out.write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"passed":True,"total_checks":receipt["total_checks"],"controls":counts,"receipt":str(out)},indent=2))
