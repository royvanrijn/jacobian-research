#!/usr/bin/env python3
"""Exact proof identities and independent verification of the complete30-state certificate."""
from __future__ import annotations
from pathlib import Path
from fractions import Fraction as Q
import copy,hashlib,json,math,platform
import sympy as s
HERE=Path(__file__).resolve().parent;COUNTS={};p=7
def check(k,n=1):COUNTS[k]=COUNTS.get(k,0)+n
def tr(a):
 while len(a)>1 and not a[-1]:a.pop()
 return a
def zero(a):return len(a)==1 and not a[0]
def add(a,b,m=p):
 c=[0]*max(len(a),len(b))
 for i,x in enumerate(a):c[i]+=x
 for i,x in enumerate(b):c[i]+=x
 return tr([x%m for x in c])
def neg(a,m=p):return [(-x)%m for x in a]
def sub(a,b,m=p):return add(a,neg(b,m),m)
def mul(a,b,m=p):
 if zero(a) or zero(b):return [0]
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  if x:
   for j,y in enumerate(b):
    if y:c[i+j]+=x*y
 return tr([x%m for x in c])
def div(a,b):
 a=a[:];q=[0]*max(1,len(a)-len(b)+1);ib=pow(b[-1],-1,p)
 while not zero(a) and len(a)>=len(b):
  k=len(a)-len(b);v=a[-1]*ib%p;q[k]=v
  a=sub(a,[0]*k+[x*v%p for x in b])
 return tr(q),a
def gcd(a,b):
 while not zero(b):a,b=b,div(a,b)[1]
 if zero(a):return a
 z=pow(a[-1],-1,p)
 return [x*z%p for x in a]

def rem(a,q):return div(a,q)[1]
def powmod(a,n,q):
 out=[1]
 while n:
  if n%2:out=rem(mul(out,a),q)
  a=rem(mul(a,a),q);n//=2
 return out
def derivative(a):return tr([(i*a[i])%7 for i in range(1,len(a))] or [0])
def validate_certificate(cert):
 assert cert['prime']==7 and cert['period']==30
 q=cert['modulus'];vv=cert['v_representative'];states=cert['states'];X=[0,1]
 assert q==[2,6,6,3,5,0,0,6,3,4,1] and len(states)==30 and states[0]==[1]
 g=[1]*6;gp=derivative(g);h=add(g,mul(X,gp));den=mul([0]*5+[1],h)
 assert gcd(den,q)==[1]
 assert rem(add(mul(vv,den),gp),q)==[0]
 check('orbit_denominator_and_rational_representative_controls',2)
 state=[1]
 for j in range(30):
  assert states[j]==state,(j,states[j],state)
  assert gcd(sub(state,X),q)==[1]
  state=rem(sub([1],mul(vv,powmod(state,7,q))),q)
  check('complete_orbit_transition_and_unit_states',2)
 assert state==[1]
 check('complete_orbit_period_return')
cert_path=HERE/'gap_six_two_power_plus_one_orbit.json'
cert=json.loads(cert_path.read_text());validate_certificate(cert)
bad=copy.deepcopy(cert);bad['states'][0]=[2]
try:validate_certificate(bad)
except AssertionError:check('corrupted_orbit_certificate_rejected')
else:raise AssertionError('corrupted state accepted')

S,v=s.symbols('S v')
Ds=[]
for b in range(5):
 f=[(j+1)*pow(6-j,b,7)%7 for j in range(7)]
 D=[sum((-1)**j*math.comb(a,j)*f[j] for j in range(a+1))%7 for a in range(7)]
 assert D[6]==0 and D[5]==(-int(b==4))%7 and sum(D[:6])%7==0
 Ds.append(D);check('universal_difference_hyperplane_controls',3)
R0=[s.expand(D[0]+v*(2*S*D[5]-S*S*D[4])) for D in Ds]
R1=[1+S]
for a in range(1,5):
 z=R1[-1];R1.append(s.expand(5*z+(S-1)*s.diff(z,S)+v*s.diff(z,v)+R0[a-1]))
for b in range(5):
 target=S*pow(6,b,7)+pow(5,b,7)+(3*v*S*S if b==4 else 0)
 assert s.Poly(R1[b]-target,S,v,modulus=7).is_zero
 check('universal_divided_row_raising_identities')
g=sum(S**i for i in range(6));gp=s.diff(g,S);h=g+S*gp
qq=s.Poly(s.expand(-gp*h*(1+5*S)+6*S*S*gp*gp-h*h),S,modulus=7)
qcoeff=[int(z)%7 for z in reversed(qq.all_coeffs())]
qcoeff=[z*pow(qcoeff[-1],-1,7)%7 for z in qcoeff]
assert qcoeff==cert['modulus']
check('universal_rank_eliminant_identity')
# Formal all-r carry identities use these exact small universal polynomials.
T,C=s.symbols('T C');P=T+3*T**2+5*T**3+5*T**4+3*T**5+T**6
R=1+T+7*(C+P);t_rel=T-T**7+7*C
difference=s.Poly(s.expand(R-R**7-t_rel),T,C)
assert all(int(z)%49==0 for z in difference.coeffs())
check('universal_mod49_R_carry_identity')
a=s.Symbol('a')
assert s.Poly(-27+6*a-(1-a),a,modulus=7).is_zero
check('universal_final_coefficient_reduction')
assert Q(1,math.factorial(7))+Q(1,7)==Q(103,720)
assert (103%7)*pow(720,-1,7)%7==2
check('exact_Artin_Hasse_seventh_coefficient',2)

def modq(a):
 a=Q(a);assert a.denominator%7
 return a.numerator%7*pow(a.denominator,-1,7)%7
def act(f,c):
 out=[0]*(len(f)+1)
 for j,x in enumerate(f):out[j]+=(c-j)*x;out[j+1]+=x
 return out
def rawcols(n):
 B=[1];c=n+5
 for ell in range(n):
  out=act(B,c)
  for j,x in enumerate(B):out[j]+=6*ell*x
  B=out
 cols=[]
 for a in range(5):
  if a:B=act(B,c)
  cols.append(B)
 return cols
profiles=[]
for r in (1,2,3):
 n=2*7**r+1;den=6**n*math.factorial(n);cols=rawcols(n)
 trace=1+sum((-1)**i*v**((7**i-1)//6) for i in range(1,r+1))
 for b,B in enumerate(cols):
  D=Ds[b]
  for a in range(2,6):
   actual=sum(modq(Q(math.factorial(a)*7**ell*B[a+6*ell],den))*v**ell
       for ell in range(max(0,(len(B)-1-a)//6+1)))
   target=math.comb(a,2)*(D[a]-2*trace*D[a-1]+trace**2*D[a-2])
   assert s.Poly(actual-target,v,modulus=7).is_zero
   check('literal_complete_ordinary_boundary_rows')
  for a in (0,1):
   actual=sum(modq(Q(math.factorial(a)*7**ell*B[a+6*ell],7*den))*v**ell
       for ell in range(max(0,(len(B)-1-a)//6+1)))
   target=(D[0]+v*(2*trace*D[5]-trace**2*D[4]) if not a
           else trace*pow(6,b,7)+pow(5,b,7)+(3*v*trace**2 if b==4 else 0))
   assert s.Poly(actual-target,v,modulus=7).is_zero
   check('literal_complete_divided_boundary_rows')
 profiles.append({'r':r,'terminal':n,'all_auxiliary_column_coefficients':sum(len(B) for B in cols)})
# Low-order diagonal normalization, as an independent fixed-shape recurrence.
t,z=s.symbols('t z');R=s.Integer(1)
for n in range(1,5):
 a=s.Symbol('aa')
 test=R+a*t**n
 coefficient=s.expand(test-test**7-t).coeff(t,n)
 R=s.expand(R+s.solve(coefficient,a)[0]*t**n)
H=1-R
exponential=sum(z**j*H**j/math.factorial(j) for j in range(5))
for eta in (-1,0,1):
 pref=s.series(6*R**(1-eta)/(7*R**6-1),t,0,5).removeO()
 gf=s.expand(pref*exponential)
 for n in range(5):
  B=[1];c=n+eta+5
  for ell in range(n):
   out=act(B,c)
   for j,x in enumerate(B):out[j]+=6*ell*x
   B=out
  target=sum(Q(x,6**n*math.factorial(n))*z**j for j,x in enumerate(B))
  assert s.expand(gf.coeff(t,n)-target)==0
  check('literal_diagonal_GF_normalization_profiles')
imports=['ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex',
         'PHYSICAL_PRIME_TOTAL_DEGREE_REDUCTION.tex','gap_three_hurwitz_structure.tex']
inputs=[HERE/'GAP_SIX_TWO_POWER_PLUS_ONE_SFC.tex',Path(__file__).resolve(),cert_path]+[HERE/x for x in imports]
out={'passed':True,'scope':'All-r formal/digit proof with a complete30-state quotient-ring certificate; literal orders only check identities.',
 'uniform_terminal_family':'N=2*7^r+1 for everyinteger r>=1','all_index_gap_six_SFC_claimed':False,
 'controls':COUNTS,'literal_profiles':profiles,
 'orbit_certificate':cert,'orbit_certificate_sha256':hashlib.sha256(cert_path.read_bytes()).hexdigest(),
 'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
 'python':platform.python_version(),'sympy':s.__version__,
 'independent_replay':False,'formal_verification':False,'external_review':False,
 'novelty_claim':False,'priority_claim':False,
 'command':'.venv/bin/python research/factorial/image-rigidity/gap_six_two_power_plus_one_sfc.py'}
(HERE/'gap_six_two_power_plus_one_sfc_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'passed':True,'controls':COUNTS,'total_controls':sum(COUNTS.values())},indent=2))
