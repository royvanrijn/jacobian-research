"""Exact internal controls for ordinary support {3,7} residue reduction.
The unbounded claims are proved in the retained TeX, not by these profiles.
Run with the repository virtual environment for SymPy.
"""
from pathlib import Path
from math import comb, factorial
import hashlib,json,platform
import sympy as sp
HERE=Path(__file__).resolve().parent
PROOF=HERE/"LOWER_THREE_FOURTH_POWER_RESIDUE_REDUCTION.tex"
IMPORTS={
 "DIAGONAL_PRIME_GAP_SFC.tex":"db034b5fb29ece9c3a0cb8cbce4c82f692b933e3679b27ee823db5ca5e7f98eb",
 "LOWER_TWO_FOURTH_POWER_GAP_SFC.tex":"e238569c21a0ccc7fb507d910608a734461fefaab4d96d5df0124afe73407631",
 "ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex":"9721978bcb795af9169c695c5fb8348974069c5148e5e0caa8bd3c4d65650e07",
}
for name,digest in IMPORTS.items():
 assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,(name,"changed exact import")
counts={}
def count(k,n=1):counts[k]=counts.get(k,0)+n
M,x,y,kap,u,n,z=sp.symbols("M x y kappa u n z")
# Complete generic rising quotient and triangular column transformation.
C=M+4
V=sp.rf(x,2)+(M+6)*x+(M+5)*(M+6)
assert sp.expand(sp.rf(x,3)-sp.rf(C,3)-(x-C)*V)==0
count("generic_rising_quotient_identities")
Bv=sp.Matrix([1,0,0,0]);Cv=sp.Matrix([0,1,0,0])
Dv=sp.Matrix([0,0,1,0]);Ev=sp.Matrix([0,0,0,1])
L0=Cv;L1=(M+7)*Cv-Dv
L2=(M+7)*(M+6)*Cv-(2*M+17)*Dv+Ev
old=[L0,L1+kap*u*Bv,L2+kap*(M+6)*u*Bv]
new=[old[0],(M+7)*old[0]-old[1]]
new.append(old[2]-(M+7)*(M+6)*new[0]+(2*M+17)*new[1])
wanted=[Cv,Dv-kap*u*Bv,Ev-kap*(M+11)*u*Bv]
for a,b in zip(new,wanted):
 for c in a-b:
  assert sp.expand(c)==0;count("generic_cyclic_column_coefficients")
assert sp.expand(new[2][0]+kap*(M+10)*u)!=0
count("corrupted_column_coefficient_rejections")
ordinary_kappa=4*(n+1)/((7*n+5)*(7*n+6)*(7*n+7))
assert sp.cancel(ordinary_kappa-4/(7*(7*n+5)*(7*n+6)))==0
count("generic_ordinary_kappa_cancellations")
R=[sp.prod(x+4*(n+j) for j in range(i)) for i in range(1,4)]
for i,r in enumerate(R,1):
 assert all(int(c)%2==0 for c in sp.Poly(sp.expand(r-x**i),x,n).coeffs())
 count("generic_reduced_monic_basis_identities")
 assert all(int(c)%4==0 for c in sp.Poly(sp.expand(r.subs(x,1)-1),n).coeffs())
 count("generic_value_congruences_mod4")
 quotient=sp.cancel(r.subs(x,0)/(4*n))
 assert sp.Poly(quotient,n).domain==sp.ZZ
 count("generic_constant_factor_divisibilities")
for i,residue in [(0,1),(1,0),(2,1)]:
 val=sp.diff(R[i],x).subs(x,1)
 assert all(int(c)%2==0 for c in sp.Poly(sp.expand(val-residue),n).coeffs())
 count("generic_derivative_congruences_mod2")
assert any(int(c)%2 for c in sp.Poly(sp.diff(R[2],x).subs(x,1),n).coeffs())
count("corrupted_derivative_congruence_rejections")

def aux(index,gamma):
 p=sp.Integer(1)
 for j in range(index):
  p=sp.expand((z+4*j+gamma+3)*p-z*sp.diff(p,z))
 return p
def cyclic(p):
 poly=sp.Poly(p,z)
 return [sp.expand(sum(poly.nth(4*k+j)*u**k for k in range(poly.degree()//4+1))) for j in range(4)]
literal=[]
for index in range(1,11):
 gamma=3*index+1;MM=7*index+1
 B,C,D,F=[aux(j,gamma) for j in [index,index+1,index+2,index+3]]
 LL1=z*(sp.diff(C,z)-C)
 LL2=z*z*(sp.diff(C,z,2)-2*sp.diff(C,z)+C)
 assert sp.expand(LL1-((MM+7)*C-D))==0
 assert sp.expand(LL2-((MM+7)*(MM+6)*C-(2*MM+17)*D+F))==0
 count("literal_auxiliary_triangular_identities",2)
 kk=sp.Rational(4,7*(7*index+5)*(7*index+6))
 cols=[cyclic(C),cyclic(LL1+kk*u*B),cyclic(LL2+kk*(MM+6)*u*B)]
 second=[sp.expand((MM+7)*a-b) for a,b in zip(cols[0],cols[1])]
 third=[sp.expand(c-(MM+7)*(MM+6)*a+(2*MM+17)*b) for a,b,c in zip(cols[0],second,cols[2])]
 target=[cyclic(C),cyclic(D-kk*u*B),cyclic(F-kk*(MM+11)*u*B)]
 for col,want in zip([cols[0],second,third],target):
  for a,b in zip(col,want):
   assert sp.expand(a-b)==0;count("literal_full_cyclic_entries")
 literal.append({"n":index,"all_twelve_entries_checked":True})
# Complete prime-stage linear algebra used in the exact boundary necessity.
for c in [0,1]:
 A=sp.Matrix([[c,0],[1,(c-1)%2]])
 A2=A*A
 assert all(int(a-b)%2==0 for a,b in zip(A2,A))
 count("complete_prime_stage_idempotent_cases")
 seed=sp.Matrix([1,0]);v=A*seed
 assert any(int(a)%2 for a in v)
 assert all(int(a)%2==0 for a in A*(A+sp.eye(2))*seed)
 count("complete_prime_stage_constant_and_linear_cases",2)
# Full residue periods are finite complete congruence certificates.
bad16=[j for j in range(16) if (7*j+5)*(7*j+6)%16==0]
bad32=[j for j in range(32) if (7*j+5)*(7*j+6)%32==0]
bad64=[j for j in range(64) if
       (j%2==1 and (7*j+5)*(7*j+6)%32==0) or
       (j%2==0 and (7*j+5)*(7*j+6)%64==0)]
assert bad16==[6,13] and bad32==[13,22] and bad64==[13,45,54]
count("complete_residue_states",16+32+64)
def v2(a):
 assert a!=0
 a=abs(int(a));return (a&-a).bit_length()-1
def vf2(a):
 total=0
 while a:
  a//=2;total+=a
 return total
strict=0
for index in range(1,65):
 gamma=3*index+1
 for k in range(1,index+1):
  poch=factorial(3*index+4*k)//factorial(3*index)
  coeff=comb(index,k)*poch
  literal_val=v2(coeff)
  legendre=(vf2(index)-vf2(k)-vf2(index-k)+vf2(3*index+4*k)-vf2(3*index))
  assert literal_val==legendre
  count("direct_integer_vs_Legendre_coefficient_valuations")
  assert literal_val>=3*k
  count("literal_base_radius_coefficients")
  if k>=2:
   assert v2(poch)>=3*k+k//2
   count("literal_interval_2_4_8_bounds")
  if index%2==0 or index%16==13:
   assert literal_val>3*k
   strict+=1
count("literal_strict_radius_coefficients",strict)
# Outside the strict-radius hypothesis the bound need not be strict.
assert v2(factorial(7)//factorial(3))==3
count("outside_strict_radius_hypothesis_controls")
receipt={
 "kind":"internal exact replay; unbounded scope follows from the written proof",
 "passed":True,
 "scope":{
  "ordinary_support_three_seven":True,"lower_indices":"n>=1",
  "proved_terminal_classes":"N not14,46,55 mod64",
  "remaining_lower_classes_mod64":bad64,
  "remaining_terminal_classes_mod64":[14,46,55],
  "all_index_support_three_seven_closed":False,
  "arbitrary_positive_Gamma_shape_extension":False,
  "all_ordinary_univariate_binomials":False,
  "unrestricted_FC_or_SFC":False,
 },
 "controls":counts,"literal_cyclic_profiles":literal,
 "complete_residue_certificates":{"mod16":bad16,"mod32":bad32,"mod64":bad64},
 "exact_imports":{
  "mixed_rank":"DIAGONAL_PRIME_GAP_SFC.tex thm:offset-rank",
  "boundary":"LOWER_TWO_FOURTH_POWER_GAP_SFC.tex thm:boundary and nonresonant unit-basis proof",
  "stage_change":"ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex lem:artin-hasse-closure",
  "nonintegral_constant_normalization":"ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex lem:frobenius-noncollision",
 },
 "source_sha256":{
  PROOF.name:hashlib.sha256(PROOF.read_bytes()).hexdigest(),
  Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),**IMPORTS,
 },
 "python":platform.python_version(),"sympy":sp.__version__,
 "independent_replay":False,"formal_verification":False,"external_review":False,
 "novelty_claim":False,"priority_claim":False,
 "command":".venv/bin/python research/factorial/image-rigidity/lower_three_fourth_power_residue_reduction.py",
}
out=HERE/"lower_three_fourth_power_residue_reduction_verification.json"
out.write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"passed":True,"controls":counts,"receipt":str(out)},indent=2))
