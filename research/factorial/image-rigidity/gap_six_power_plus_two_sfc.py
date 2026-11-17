#!/usr/bin/env python3
"""Exact necessary minor identities, complete6-state orbit, and full literal row controls."""
from pathlib import Path
from fractions import Fraction as Q
import copy,hashlib,itertools,json,math,platform
import sympy as s
HERE=Path(__file__).resolve().parent
p=7;S,v=s.symbols("S v")
COUNTS={}
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

certificate_path=HERE/'gap_six_power_plus_two_certificate.json'
certificate=json.loads(certificate_path.read_text())
assert certificate['schema']=='gap-six-power-plus-two-exact-certificate-v1'
assert certificate['modulus']==7 and certificate['monomial_exponents']==['v','S']
data={'matrix':certificate['matrix'],'minors':certificate['maximal_minors']}
witness=certificate['minor_membership'];orbit=certificate['orbit']
assert witness['modulus']==7 and len(witness['witnesses'])==2
assert all(len(row['multipliers'])==6 for row in witness['witnesses'])
check('mandatory_certificate_schema')
def bp(expr):
 return {m:int(z)%7 for m,z in s.Poly(s.sympify(expr,locals={'S':S,'v':v}),v,S,modulus=7).terms() if int(z)%7}
def ba(a,b):
 out=a.copy()
 for m,z in b.items():
  out[m]=(out.get(m,0)+z)%7
  if not out[m]:out.pop(m)
 return out
def bm(a,b):
 out={}
 for (i,j),x in a.items():
  for (k,l),y in b.items():
   m=(i+k,j+l);out[m]=(out.get(m,0)+x*y)%7
 return {m:x for m,x in out.items() if x}
def bs(a,c):return {m:x*c%7 for m,x in a.items() if x*c%7}
def decode(a):return {(i,j):z for i,j,z in a}
def det(A):
 out={}
 for perm in itertools.permutations(range(len(A))):
  term={(0,0):1}
  for i,j in enumerate(perm):term=bm(term,A[i][j])
  sign=(-1)**sum(perm[i]>perm[j] for i in range(len(A)) for j in range(i+1,len(A)))
  out=ba(out,bs(term,sign))
 return out
M=[[bp(x) for x in row] for row in data['matrix']]
minors=[]
for omit in range(6):
 f=det([row for i,row in enumerate(M) if i!=omit]);minors.append(f)
 assert f==bp(data['minors'][omit])==decode(witness['input_minors'][omit])
 check('exact_maximal_minors')
# Regenerate ALL30 entries solely from the uniform seeds and exact raising.
def r(f):return s.Poly(s.expand(f),S,v,modulus=7).as_expr()
def D(f):return r(v*s.diff(f,v)+(S-1)*s.diff(f,S))
rows=[[s.Integer(1)],[5+4*S],[1+2*S],[4*S],[s.Integer(0)],[s.Integer(0)]]
for j in range(4):
 old=[row[-1] for row in rows]
 for a in range(6):
  cross=v*old[5] if a==0 else old[0] if a==1 else 0 if a==2 else a*old[a-1]
  rows[a].append(r(-a*old[a]+D(old[a])+cross))
for a in range(6):
 for j in range(5):
  assert bp(rows[a][j])==M[a][j];check('all_universal_matrix_entries_from_seeds')
def validate_witnesses(w):
 for row in w['witnesses']:
  out={}
  for mult,f in zip(row['multipliers'],minors):out=ba(out,bm(decode(mult),f))
  assert out==decode(row['polynomial'])
validate_witnesses(witness);check('explicit_necessary_ideal_membership_identities',2)
bad=copy.deepcopy(witness);bad['witnesses'][0]['multipliers'][0][0][2]^=1
try:validate_witnesses(bad)
except AssertionError:check('corrupted_ideal_witness_rejected')
else:raise AssertionError('bad witness accepted')
q=orbit['modulus'];vv=orbit['v_representative'];X=[0,1]
assert decode(witness['witnesses'][0]['polynomial'])==ba({(1,0):1},bs({(0,i):z for i,z in enumerate(vv) if z},-1))
assert decode(witness['witnesses'][1]['polynomial'])=={(0,i):z for i,z in enumerate(q) if z}
def validate_orbit(d):
 assert d['prime']==7 and d['preperiod']==0 and d['period']==6 and len(d['states'])==6
 assert gcd(q,vv)==gcd(q,derivative(q))==[1]
 assert rem(add(add(powmod(vv,8,q),vv),[1]),q)==[0]
 expected=[[1],sub([1],vv),neg(mul([2],vv)),sub(neg([1]),mul([2],vv)),sub(neg([1]),vv),[0]]
 assert d['states']==expected
 state=[1]
 for j in range(6):
  assert d['states'][j]==state
  assert gcd(sub(state,X),q)==[1]
  state=rem(sub([1],mul(vv,powmod(state,7,q))),q)
 assert state==[1]
validate_orbit(orbit);check('complete_orbit_states_transitions_units',14)
bad=copy.deepcopy(orbit);bad['states'][0]=[2]
try:validate_orbit(bad)
except AssertionError:check('corrupted_complete_orbit_rejected')
else:raise AssertionError('bad orbit accepted')
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
for rr in (1,2,3):
 n=7**rr+2;den=6**n*math.factorial(n);cols=rawcols(n)
 trace=1+sum((-1)**i*v**((7**i-1)//6) for i in range(1,rr+1))
 for j,B in enumerate(cols):
  for a in range(6):
   actual=sum(modq(Q(math.factorial(a)*7**ell*B[a+6*ell],den*(7 if a<2 else 1)))*v**ell
     for ell in range(max(0,(len(B)-1-a)//6+1)))
   target=rows[a][j].subs(S,trace)
   assert s.Poly(actual-target,v,modulus=7).is_zero,(n,a,j)
   check('literal_full_entry_polynomial_profiles')
 profiles.append({'r':rr,'terminal':n,'entries':30,'raw_auxiliary_coefficients':sum(len(B) for B in cols)})
# Exact coefficient combinatorics underlying the uniform GF seed proof.
a=s.symbols('a')
assert s.Poly(3*a-30-(5-4*a),a,modulus=7).is_zero
assert 5*3%7==1
assert s.Poly(1-2*a-(1+2*(-a)),a,modulus=7).is_zero
check('uniform_seed_coefficient_identities',3)
imports={
 'GAP_SIX_TWO_POWER_PLUS_ONE_SFC.tex':'7c8b7be9cd416803f0d43e0f76248cb5006f94ffa9feafc5295115d9cf61ffc6',
 'ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex':'9721978bcb795af9169c695c5fb8348974069c5148e5e0caa8bd3c4d65650e07',
 'PHYSICAL_PRIME_TOTAL_DEGREE_REDUCTION.tex':'42c5fc615400a49dd3ae4f19a9a832381607e69ce4f6ff6d255b57ca4cfd156f',
 'gap_three_hurwitz_structure.tex':'22c690a45976c9b72b210ff76bbcd3a7084ddba0aae4af0c7f40ce2b109f705c',
}
for name,digest in imports.items():
 assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,('Imported proof changed',name)
 check('frozen_mathematical_imports')
inputs=[HERE/'GAP_SIX_POWER_PLUS_TWO_SFC.tex',Path(__file__).resolve(),certificate_path]+[HERE/name for name in imports]
out={'passed':True,'kind':'internal exact finite-certificate replay; unrestricted depth is proved in the source',
 'scope':{'ordinary_terminal_family':'N=7^r+2,r>=1','every_digit_depth':True,
  'positive_real_Gamma_extension':False,'all_index_gap_six_SFC':False,
  'ordinary_degree_seven_SFC':False,'all_binomial_SFC':False,'unrestricted_FC_or_SFC':False},
 'controls':COUNTS,'total_controls':sum(COUNTS.values()),'literal_profiles':profiles,'complete_orbit_period':6,
 'source_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs},
 'python':platform.python_version(),'sympy':s.__version__,
 'orbit_verifier':'pure Python integer coefficient lists; no factorization or Groebner input',
 'independent_replay':False,'formal_verification':False,'external_review':False,
 'novelty_claim':False,'priority_claim':False,
 'command':'.venv/bin/python research/factorial/image-rigidity/gap_six_power_plus_two_sfc.py'}
(HERE/'gap_six_power_plus_two_sfc_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'passed':True,'controls':COUNTS,'total_controls':sum(COUNTS.values())},indent=2))
