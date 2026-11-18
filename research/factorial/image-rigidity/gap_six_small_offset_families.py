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
 3:[6,1+2*S,S-2],
 4:[1,3+2*S,2*S-3,-S-1],
 5:[6+v*S,3-3*S,2-2*S,2*S-2,3*S-3],
 6:[1+v*(S-1),S*(v-1),0,0,0,0]}
P=(x+y)+3*(x+y)**2+5*(x+y)**3+5*(x+y)**4+3*(x+y)**5+(x+y)**6
AA=-x-S*y
for c in range(3,7):
 for a in range(1,c):
  z=s.Poly(s.expand(P*AA**(a-1)*(AA-a*(1+x+y))),x,y).coeff_monomial(x**c*y)
  if c==6 and a==1:z+=-1+v*S
  same("generic_divided_seed_coefficient_identities",z,SEEDS[c][a])
  check("all_carry_degree_bound",a+6<c+7)
 if c<=4:z=(c+1)*[0,1,3,5,5,3,1][c+1]
 elif c==5:z=6+6*v*s.Poly(AA**6,x,y).coeff_monomial(x**5*y)
 else:z=1+6*v*s.Poly((x+y)*AA**6,x,y).coeff_monomial(x**6*y)
 same("generic_zero_row_seed_identities",z,SEEDS[c][0])
same("generic_unit6_carry_A",s.binomial(7,1)/7*(-1)**6*(-S),-S)
same("generic_unit6_carry_T",s.binomial(7,1)/7,1)
path=HERE/"gap_six_small_offset_families_certificate.json"
data=json.loads(path.read_text())
check("certificate_ring",data["modulus"]==7 and data["monomial_exponents"]==["v","S"])
check("complete_cases",[a["c"] for a in data["cases"]]==[3,4,5,6])
matrices={};case_summaries=[]
for case in data["cases"]:
 c=case["c"];C=c+5
 ordinary=[[s.Integer(0)] for _ in range(6)]
 if c<6:ordinary[c][0]=red((-1)**c*(1+c*S))
 if c+1<6:ordinary[c+1][0]=red((c+1)*(-1)**(c+1)*S)
 for a in range(6):
  baseline=s.Poly((1+x+y)*AA**a,x,y).coeff_monomial(x**c*y)
  same("generic_ordinary_seed_coefficients",baseline,ordinary[a][0])
 for j in range(4):
  old=[row[-1] for row in ordinary]
  for a in range(6):ordinary[a].append(red((C-a)*old[a]+deriv(old[a])+(a*old[a-1] if a else 0)))
 for a in range(c):
  for j in range(5):same("zero_ordinary_rows_before_division",ordinary[a][j],0)
 rows=[[red(a)] for a in SEEDS[c]]
 for j in range(4):
  old=[row[-1] for row in rows]
  rows[0].append(red(C*old[0]+deriv(old[0])+v*ordinary[5][j]))
  for a in range(1,c):rows[a].append(red((C-a)*old[a]+deriv(old[a])+a*old[a-1]))
 M=s.Matrix(rows+ordinary[c:]).applyfunc(red);matrices[c]=M
 for a in range(6):
  for j in range(5):check("all_universal_matrix_entries",sparse(M[a,j])==decode(case["matrix"][a][j]))
 minors=[];domain=s.GF(7).poly_ring(S,v)
 for omit in range(6):
  B=M.extract([a for a in range(6) if a!=omit],range(5))
  det=red(DomainMatrix.from_Matrix(B).convert_to(domain).det().as_expr())
  minors.append(sparse(det));check("all_recomputed_maximal_minors",minors[-1]==decode(case["maximal_minors"][omit]))
 targets=[]
 for ident in case["identities"]:
  total={}
  check("all_six_minor_multipliers",len(ident["multipliers"])==6)
  for multiplier,minor in zip(ident["multipliers"],minors):total=dadd(total,dmul(decode(multiplier),minor))
  target=decode(ident["polynomial"])
  check("all_exact_necessary_relation_identities",total==target);targets.append(target)
 if c<6:
  q,W=trim(case["q"]),trim(case["W"])
  targetq={(0,i):z for i,z in enumerate(q) if z}
  targetw=dadd({(1,0):1},{(0,i):z for i,z in enumerate(W) if z},-1)
  check("relation_targets_match_orbit_inputs",targetq in targets and targetw in targets)
  check("quotient_monic_degree",q[-1]==1 and len(q)-1==5*c)
  check("quotient_parameter_is_unit",pgcd(q,W)==[1])
  state=[1];period=case["period"];states=case["orbit_states"]
  check("complete_period_length",len(states)==period)
  branch_depths=[]
  for j,stored in enumerate(states):
   check("every_orbit_state",state==trim(stored))
   g=pgcd(padd(state,[0,1],-1),q)
   expected=[3,1] if c==5 and j%6==1 else [1]
   check("every_target_gcd",g==expected)
   if g!=[1]:branch_depths.append(j)
   state=pmod(padd([1],pmul(W,frob(state,q)),-1),q)
  check("orbit_return_to_one",state==[1])
  if c==5:
   ev=lambda p,a:sum(z*pow(a,i,7) for i,z in enumerate(p))%7
   qp=trim([i*z for i,z in enumerate(q)][1:])
   check("simple_linear_residue_branch",ev(q,4)==0 and ev(qp,4)==5 and ev(W,4)==4)
   check("local_disk_divided_difference_unit",5*((4-1)*pow(4,-1,7))%7==2)
   check("residual_point_is_compatible",all(sum(z*pow(4,ev0+es,7) for (ev0,es),z in f.items())%7==0 for f in minors))
  case_summaries.append({"offset":c,"q_degree":len(q)-1,"period":period,
   "exceptional_depths_mod_period":branch_depths,"all_other_target_gcds_are_one":True})
 else:
  check("necessary_S_six",{(0,6):1} in targets)
  f=next(a for a in targets if a!={(0,6):1})
  at0={ev0:z for (ev0,es),z in f.items() if es==0}
  expected=s.Poly((v-1)**5,v,modulus=7)
  check("necessary_v_one_at_S_zero",at0=={j:int(expected.nth(j))%7 for j in range(6)})
  check("residue_point_is_compatible",all(sum(z for (ev0,es),z in f.items() if es==0)%7==0 for f in minors))
  check("complete_even_odd_S_at_v_one_cycle",[(a,(1-a**7)%7) for a in (1,0)]==[(1,0),(0,1)])
  case_summaries.append({"offset":6,"necessary_residue_point":{"S":0,"v":1},
   "first_matrix_even_depth_exclusion":True,"first_matrix_odd_depth_exclusion":False,
   "full_trace_exclusion_all_depths":True})
 if c==3:
  total={};changed=False
  for mm,minor in zip(case["identities"][0]["multipliers"],minors):
   d=decode(mm)
   if d and not changed:
    mon=next(iter(d));d[mon]=(d[mon]+1)%7;d={a:z for a,z in d.items() if z};changed=True
   total=dadd(total,dmul(d,minor))
  check("corrupt_membership_certificate_rejected",total!=targets[0])
check("incorrect_constant_six_raising_rejected",red(matrices[3][0,1]-1)!=0)
# Actual integer auxiliaries: all columns and all polynomial coefficients.
def act(f,C):
 out=[0]*(len(f)+1)
 for j,z in enumerate(f):out[j]+=(C-j)*z;out[j+1]+=z
 return out
def vp_integer(n):
 val=0
 while n and n%7==0:n//=7;val+=1
 return val
profiles=[]
for c in range(3,7):
 for r in (1,2):
  N=7**r+c;C=N+5;B=[1]
  for ell in range(N):
   new=act(B,C)
   for j,z in enumerate(B):new[j]+=6*ell*z
   B=new
  columns=[]
  for j in range(5):
   if j:B=act(B,C)
   columns.append(B)
  Sr=[1]
  for i in range(1,r+1):Sr=padd(Sr,[0]*((7**i-1)//6)+[(-1)**i])
  den=6**N*math.factorial(N);nu=vp_integer(den);iv=pow(den//7**nu%7,-1,7);coefficient_count=0
  for a in range(6):
   for j in range(5):
    predicted=eval_sp(sparse(matrices[c][a,j]),Sr);actual=[];divisor=7**(nu+int(a<c))
    for ell in range(max(0,(len(columns[j])-1-a)//6+1)):
     num=math.factorial(a)*7**ell*columns[j][a+6*ell]
     check("literal_global_division_integrality",num%divisor==0)
     actual.append(num//divisor*iv%7)
    actual=trim(actual)
    check("all_actual_boundary_entry_polynomials",actual==predicted)
    coefficient_count+=max(len(actual),len(predicted))
  profiles.append({"offset":c,"r":r,"terminal":N,"full_coefficient_controls":coefficient_count})
# Trace parity is checked over Z/49, not merely over F7.
MOD49=49
def trim49(a):
 a=[int(z)%49 for z in a]
 while len(a)>1 and not a[-1]:a.pop()
 return a or [0]
def add49(a,b,scale=1):
 out=[0]*max(len(a),len(b))
 for i,z in enumerate(a):out[i]+=z
 for i,z in enumerate(b):out[i]+=scale*z
 return trim49(out)
def mul49(a,b):
 out=[0]*(len(a)+len(b)-1)
 for i,z in enumerate(a):
  for j,w in enumerate(b):out[i+j]+=z*w
 return trim49(out)
q49=[48,7,28,35,14,21,27,1]
def rem49(a):
 a=trim49(a)
 while len(a)>=len(q49):
  d=len(a)-len(q49);z=a[-1]
  for j,w in enumerate(q49):a[d+j]=(a[d+j]-z*w)%49
  a=trim49(a)
 return a
def power49(a,n):
 out=[1]
 while n:
  if n&1:out=rem49(mul49(out,a))
  a=rem49(mul49(a,a));n//=2
 return out
lam=s.symbols("lam")
chi=s.Poly((lam-1)**7-15*lam**6,lam)
check("exact_trace_characteristic_polynomial",
      q49==[int(chi.nth(i))%49 for i in range(8)])
P343,P7=power49([0,1],343),power49([0,1],7)
check("full_mod49_quotient_parity_identity",P343==P7)
check("residue_F49_quotient_identity",
      all(z%7==0 for z in add49(power49([0,1],49),[0,1],-1)))
check("mod49_precision_corruption_rejected",rem49(add49(add49(P343,P7,-1),[7]))!=[0])
def trace49(n):
 vals=[7]
 for k in range(1,n+1):
  if k<=7:z=-sum(q49[7-j]*vals[k-j] for j in range(1,k))-k*q49[7-k]
  else:z=-sum(q49[7-j]*vals[k-j] for j in range(1,8))
  vals.append(z%49)
 return vals[n]
def choose(n,k):return math.comb(n,k) if 0<=k<=n else 0
def binomial_trace49(n):
 return (7+sum((7*choose(n+6*k-1,7*k)+choose(n+6*k-1,7*k-1))*pow(15,k,49)
        for k in range(1,n+1)))%49
def actual_point(n):
 total=7;leading=None
 for k in range(1,n+1):
  a=math.comb(n,k)*math.prod(range(n,n+6*k))
  check("actual_H_coefficient_integrality",a%7**(k-1)==0)
  a//=7**(k-1);total+=(-1)**k*a
  if k==n:leading=a
 check("actual_H_leading_unit",leading%7!=0)
 return total
trace_profiles=[]
check("minus_factorial_trace_parameter",(-math.factorial(6))%49==15)
for n,want in ((13,7),(55,43)):
 newton=trace49(n);binomial=binomial_trace49(n);actual=actual_point(n)
 check("three_algorithm_exact_trace_point_value",newton==binomial==actual%49==want)
 trace_profiles.append({"terminal":n,"trace_mod49":want,
  "actual_H_point_mod49":actual%49,"point_valuation":1 if want%7==0 else 0})
check("odd_trace_Eisenstein_constant",trace_profiles[0]["trace_mod49"]==7)
check("even_trace_disk_is_a_unit",trace_profiles[1]["trace_mod49"]%7!=0)

# Bridge normalization and even-moment proof are characteristic-zero identities.
g=s.symbols("g")
for k in range(5):check("generic_upper_moment_bridge",s.expand(g*s.rf(g+1,6*k)-(g+6*k)*s.rf(g,6*k))==0)
for N in range(2,8):
 for k in range(N+1):check("lower_moment_bridge",s.binomial(N-1,k)==(1-s.Rational(k,N))*s.binomial(N,k))
imports=["GAP_SIX_TWO_POWER_PLUS_ONE_SFC.tex","PHYSICAL_PRIME_TOTAL_DEGREE_REDUCTION.tex",
 "gap_three_hurwitz_structure.tex","GAP_SIX_POWER_PLUS_TWO_SFC.tex","GAP_SIX_SCALED_POWER_PLUS_ONE_SFC.tex"]
names=["GAP_SIX_SMALL_OFFSET_FAMILIES.tex",Path(__file__).name,path.name]+imports
out={"passed":True,"kind":"internal exact certificate replay; written arguments give every depth",
 "scope":{"ordinary_offsets_all_depths":[3,4,5],"unified_offsets_all_depths":[2,3,4,5,6],
 "offset_six_all_depths":True,"generic_odd_prime_trace_parity_proved_in_written_source":True,
 "all_index_gap_six_SFC":False,"all_positive_Gamma_shapes":False,"all_degree_seven_SFC":False,"unrestricted_FC_or_SFC":False},
 "controls":COUNTS,"total_controls":sum(COUNTS.values()),"finite_algebra_cases":case_summaries,
 "actual_boundary_profiles":profiles,"exact_trace_profiles":trace_profiles,"c5_injective_domain":"every Cp point reducing to4; arbitrarily ramified extensions allowed",
 "source_sha256":{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names},
 "python":platform.python_version(),"sympy":s.__version__,"orbit_verifier":"pure Python coefficient lists; no Groebner or factorization input",
 "independent_replay":False,"formal_verification":False,"external_review":False,"novelty_claim":False,"priority_claim":False,
 "command":".venv/bin/python research/factorial/image-rigidity/gap_six_small_offset_families.py"}
(HERE/"gap_six_small_offset_families_verification.json").write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps({k:out[k] for k in ("passed","total_controls","controls","finite_algebra_cases","actual_boundary_profiles")},indent=2))
