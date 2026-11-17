#!/usr/bin/env python3
"""Exact finite verifier; the written digit and formal-series proof supplies every r."""
from __future__ import annotations
import copy,hashlib,json,math,platform
from fractions import Fraction
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
HERE=Path(__file__).resolve().parent
S,v,x,y=s.symbols("S v x y")
COUNTS={}
def check(key,ok=True,n=1):
 if not bool(ok):raise AssertionError(key)
 COUNTS[key]=COUNTS.get(key,0)+n
def red(f):return s.Poly(s.expand(f),S,v,modulus=7).as_expr()
def same(key,a,b):check(key,red(a-b)==0)
def deriv(f):return red(v*s.diff(f,v)+(S-1)*s.diff(f,S))
def sparse(expr):return {m:int(c)%7 for m,c in s.Poly(expr,v,S,modulus=7).terms() if int(c)%7}
def decode(rows):
 out={}
 for a,b,c in rows:
  assert isinstance(a,int) and a>=0 and isinstance(b,int) and b>=0 and isinstance(c,int) and 0<=c<7
  assert (a,b) not in out
  if c:out[(a,b)]=c
 return out
def dadd(a,b,scale=1):
 out=a.copy()
 for m,c in b.items():
  out[m]=(out.get(m,0)+scale*c)%7
  if not out[m]:out.pop(m)
 return out
def dmul(a,b):
 out={}
 for m,c in a.items():
  for n,d in b.items():
   k=(m[0]+n[0],m[1]+n[1]);out[k]=(out.get(k,0)+c*d)%7
 return {m:c for m,c in out.items() if c}
# These separate integer-list operations verify all finite-algebra orbit data.
def trim(a):
 a=[int(c)%7 for c in a]
 while len(a)>1 and not a[-1]:a.pop()
 return a or [0]
def padd(a,b,scale=1):
 out=[0]*max(len(a),len(b))
 for i,c in enumerate(a):out[i]+=c
 for i,c in enumerate(b):out[i]+=scale*c
 return trim(out)
def pmul(a,b):
 out=[0]*(len(a)+len(b)-1)
 for i,c in enumerate(a):
  for j,d in enumerate(b):out[i+j]+=c*d
 return trim(out)
def pdiv(a,b):
 a=trim(a);b=trim(b);assert b!=[0]
 q=[0]*max(1,len(a)-len(b)+1);iv=pow(b[-1],-1,7)
 while a!=[0] and len(a)>=len(b):
  k=len(a)-len(b);c=a[-1]*iv%7;q[k]=c
  for j,d in enumerate(b):a[k+j]=(a[k+j]-c*d)%7
  a=trim(a)
 return trim(q),a
def pmod(a,b):return pdiv(a,b)[1]
def pgcd(a,b):
 while b!=[0]:a,b=b,pmod(a,b)
 return trim([c*pow(a[-1],-1,7)%7 for c in a]) if a!=[0] else [0]
def ppow(a,k,q):
 out=[1]
 while k:
  if k&1:out=pmod(pmul(out,a),q)
  a=pmod(pmul(a,a),q);k//=2
 return out
def frob(a,q):
 out=[0]*(7*(len(a)-1)+1)
 for j,c in enumerate(a):out[7*j]=c
 return pmod(out,q)
def eval_sp(p,SS):
 out=[0]
 for (ev,eS),c in p.items():
  a=pmul([c],ppow_unreduced(SS,eS))
  out=padd(out,[0]*ev+a)
 return out
def ppow_unreduced(a,k):
 out=[1]
 while k:
  if k&1:out=pmul(out,a)
  a=pmul(a,a);k//=2
 return out
SEEDS={
 3:[6,4-S,2+5*S-S*S],
 4:[1,4+S,6+6*S+S*S,S+3*S*S+S**3],
 5:[6+v*S**5,1-S,2+6*S+6*S*S,3*S+5*S*S+6*S**3,5*S*S+3*S**3+6*S**4],
 6:[1-v*S**5*(S-1),S-2+v*S**6]+[S**(a-2)*(1+5*S+S*S) for a in range(2,6)]}
P=(x+y)+3*(x+y)**2+5*(x+y)**3+5*(x+y)**4+3*(x+y)**5+(x+y)**6
AA=-x-S*y
for b in range(3,7):
 for a in range(1,b):
  expr=P*AA**(a-1)*(AA-a*(1+x+y))
  c=s.Poly(s.expand(expr),x,y).coeff_monomial(x*y**b)
  if b==6 and a==1:c+=-1+v*S**6
  same("generic_P_seed_coefficient_identities",c,SEEDS[b][a])
  check("selected_degree_excludes_all_carries",a+6<b+7)
 if b<=4:c=(b+1)*[0,1,3,5,5,3,1][b+1]
 elif b==5:c=6+6*v*s.Poly(AA**6,x,y).coeff_monomial(x*y**5)
 else:c=1+6*v*s.Poly((x+y)*AA**6,x,y).coeff_monomial(x*y**6)
 same("generic_zero_row_seed_identities",c,SEEDS[b][0])
same("generic_b6_carry_T",s.binomial(7,1)/7,1)
same("generic_b6_carry_A",s.binomial(7,1)/7*(-1)*(-S)**6,-S**6)
# Exact low Artin-Hasse coefficients, including the nontrivial cancellation at degree7.
def vp_frac(a):
 a=Fraction(a)
 if not a:return math.inf
 n,d=a.numerator,a.denominator;r=0
 while n%7==0:n//=7;r+=1
 while d%7==0:d//=7;r-=1
 return r
def residue_frac(a):
 a=Fraction(a);assert a.denominator%7
 return a.numerator%7*pow(a.denominator,-1,7)%7
for a in range(6):
 j=a+6
 c=Fraction(math.factorial(a))*sum((Fraction(1,7**k*math.factorial(k)*math.factorial(j-7*k)) for k in range(j//7+1)),Fraction())
 check("Artin_Hasse_low_section_integrality",vp_frac(c)>=0)
 if a==1:check("Artin_Hasse_degree7_unit",c==Fraction(103,720) and residue_frac(c)==2)
R,t=s.symbols("R t")
cleared=6*R-(R-7*t)*(7*R**6-1)
check("generic_kappa_mod49_identity",s.expand(cleared-(7*(R-R**7-t)+49*t*R**6))==0)
certificate_path=HERE/"gap_six_scaled_power_plus_one_certificate.json"
data=json.loads(certificate_path.read_text())
check("certificate_modulus_and_variable_order",data["modulus"]==7 and data["monomial_exponents"]==["v","S"])
check("complete_four_cases",[c["b"] for c in data["cases"]]==[3,4,5,6])
matrices={};case_summary=[]
for case in data["cases"]:
 b=case["b"];original=[]
 for column in range(5):
  f=(x+1)*(6-x)**column
  dif=[s.expand(sum((-1)**j*s.binomial(a,j)*f.subs(x,j) for j in range(a+1))) for a in range(6)]
  original.append([red(s.binomial(a,b)*sum(s.binomial(b,j)*(-S)**j*dif[a-j] for j in range(b+1)))
                   for a in range(b,6)])
 orig=s.Matrix(6-b,5,lambda i,j:original[j][i])
 rows=[[red(z)] for z in SEEDS[b]]
 for column in range(4):
  old=[row[-1] for row in rows]
  rows[0].append(red(6*old[0]+deriv(old[0])+v*(orig[5-b,column] if b<6 else 0)))
  for a in range(1,b):rows[a].append(red((6-a)*old[a]+deriv(old[a])+a*old[a-1]))
 M=s.Matrix(rows+orig.tolist()).applyfunc(red);matrices[b]=M
 for a in range(6):
  for column in range(5):
   check("complete_universal_matrix_entries",sparse(M[a,column])==decode(case["matrix"][a][column]))
 minors=[];domain=s.GF(7).poly_ring(S,v)
 for omit in range(6):
  A=M.extract([i for i in range(6) if i!=omit],range(5))
  det=red(DomainMatrix.from_Matrix(A).convert_to(domain).det().as_expr())
  minors.append(sparse(det))
  check("complete_maximal_minor_identities",minors[-1]==decode(case["maximal_minors"][omit]))
 q=trim(case["q"]);W=trim(case["V"])
 check("quotient_degree_and_monic",len(q)-1==5*b and q[-1]==1)
 targets={"q":{(0,i):c for i,c in enumerate(q) if c},
          "v_minus_V":dadd({(1,0):1},{(0,i):c for i,c in enumerate(W) if c},-1)}
 monomials=0
 for name in ("q","v_minus_V"):
  witness=case["identities"][name]
  check("certificate_target_matches_orbit_input",decode(witness["polynomial"])==targets[name])
  total={}
  check("complete_six_minor_multipliers",len(witness["multipliers"])==6)
  for multiplier,minor in zip(witness["multipliers"],minors):
   mm=decode(multiplier);monomials+=len(mm);total=dadd(total,dmul(mm,minor))
  check("necessary_relation_sparse_identities",total==targets[name])
 check("quotient_parameter_unit",pgcd(q,W)==[1])
 derivative=trim([(i*c)%7 for i,c in enumerate(q)][1:])
 check("quotient_squarefree_control",pgcd(q,derivative)==[1])
 state=[1];period=case["period"];states=case["orbit_states"]
 check("complete_orbit_state_count",len(states)==period)
 for j,recorded in enumerate(states):
  check("every_orbit_state",state==trim(recorded))
  check("every_orbit_target_unit",pgcd(padd(state,[0,1],-1),q)==[1])
  state=pmod(padd([1],pmul(W,frob(state,q)),-1),q)
  check("every_orbit_state_degree",len(state)<=len(q)-1)
 check("orbit_return_to_initial_state",state==[1])
 case_summary.append({"b":b,"q_degree":len(q)-1,"period":period,
                      "all_target_gcds_are_one":True,"all_transitions_and_return_verified":True,
                      "multiplier_monomials":monomials})
 # A corrupted coefficient in a nonzero multiplier must be rejected.
 if b==3:
  witness=case["identities"]["q"];total={}
  altered=False
  for rr,minor in zip(witness["multipliers"],minors):
   mm=decode(rr)
   if mm and not altered:
    mon=next(iter(mm));mm[mon]=(mm[mon]+1)%7;mm={k:c for k,c in mm.items() if c};altered=True
   total=dadd(total,dmul(mm,minor))
  check("corrupted_membership_witness_rejected",total!=targets["q"])
  bad=trim(states[1]);bad=padd(bad,[1])
  expected=pmod(padd([1],pmul(W,frob([1],q)),-1),q)
  check("corrupted_orbit_state_rejected",bad!=expected)
# Direct characteristic-zero integer auxiliary rows are built separately from the universal formulas.
def act(f,c):
 out=[0]*(len(f)+1)
 for j,a in enumerate(f):out[j]+=(c-j)*a;out[j+1]+=a
 return out
def vp_integer(n):
 v0=0
 while n and n%7==0:n//=7;v0+=1
 return v0
profiles=[]
for b in range(3,7):
 for r in (1,2):
  N=b*7**r+1;c=N+5;B=[1]
  for ell in range(N):
   new=act(B,c)
   for j,a in enumerate(B):new[j]+=6*ell*a
   B=new
  columns=[]
  for j in range(5):
   if j:B=act(B,c)
   columns.append(B)
  Sr=[1]
  for i in range(1,r+1):
   term=[0]*((7**i-1)//6)+[(-1)**i];Sr=padd(Sr,term)
  den=6**N*math.factorial(N);nu=vp_integer(den);iv=pow(den//7**nu%7,-1,7)
  coeff_count=0
  for a in range(6):
   for column in range(5):
    predicted=eval_sp(sparse(matrices[b][a,column]),Sr);actual=[]
    divisor=7**(nu+int(a<b))
    for ell in range(max(0,(len(columns[column])-1-a)//6+1)):
     num=math.factorial(a)*7**ell*columns[column][a+6*ell]
     check("literal_global_division_integrality",num%divisor==0)
     actual.append(num//divisor*iv%7)
    actual=trim(actual)
    check("complete_actual_boundary_entry_polynomials",actual==predicted)
    coeff_count+=max(len(actual),len(predicted))
  profiles.append({"b":b,"r":r,"terminal":N,"full_coefficient_controls":coeff_count})
same("correct_b4_seed_is_constant",matrices[4][0,0],1)
check("incorrect_b4_seed_model_rejected",red(matrices[4][0,0]-(1-v*S**5*(S-1)))!=0)
# Generic moment-bridge coefficient identities (not finite terminal coverage).
g=s.symbols("g")
for k in range(5):
 check("Gamma_shift_bridge_coefficients",s.expand(g*s.rf(g+1,6*k)-(g+6*k)*s.rf(g,6*k))==0)
for N0 in range(2,8):
 for k in range(N0+1):
  check("lower_order_bridge_coefficients",s.binomial(N0-1,k)==(1-s.Rational(k,N0))*s.binomial(N0,k))
imports=["GAP_SIX_TWO_POWER_PLUS_ONE_SFC.tex","ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex",
 "PHYSICAL_PRIME_TOTAL_DEGREE_REDUCTION.tex","gap_three_hurwitz_structure.tex","GAP_SIX_THIN_DIGIT_SFC.tex"]
inputs=["GAP_SIX_SCALED_POWER_PLUS_ONE_SFC.tex",Path(__file__).name,certificate_path.name]+imports
out={"passed":True,"kind":"internal exact finite-certificate replay; unrestricted digit depth comes from written proofs and periodicity",
 "scope":{"new_leading_digits":[3,4,5,6],"new_terminal_family":"N=b*7^r+1,r>=0",
  "unified_family_by_imports":"1<=b<=6,N=b*7^r+1,r>=0",
  "every_digit_depth":True,"ordinary_factorial":True,"all_Gamma_shapes":False,
  "all_index_gap_six_SFC":False,"all_degree_seven_binomial_SFC":False,"unrestricted_FC_or_SFC":False},
 "controls":COUNTS,"total_controls":sum(COUNTS.values()),"finite_algebra_cases":case_summary,
 "literal_actual_profiles":profiles,
 "source_sha256":{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in inputs},
 "python":platform.python_version(),"sympy":s.__version__,"orbit_verifier":"pure Python integer coefficient lists; no factorization or Groebner input",
 "independent_replay":False,"formal_verification":False,"external_review":False,"novelty_claim":False,"priority_claim":False,
 "command":".venv/bin/python research/factorial/image-rigidity/gap_six_scaled_power_plus_one_sfc.py"}
(HERE/"gap_six_scaled_power_plus_one_sfc_verification.json").write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps({k:out[k] for k in ("passed","total_controls","controls","finite_algebra_cases","literal_actual_profiles")},indent=2))
