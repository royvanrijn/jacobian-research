#!/usr/bin/env python3
"""Replay exact uniform matrix identities, complete finite certificates and fixed trace packets."""
from pathlib import Path
from fractions import Fraction as Q
import copy,hashlib,itertools,json,math,platform
import sympy as s
HERE=Path(__file__).resolve().parent;p=7;S,v=s.symbols("S v");COUNTS={}
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
CERT_PATH=HERE/'gap_six_two_digit_leading_two_certificate.json'
cert=json.loads(CERT_PATH.read_text())
assert cert['prime']==7 and set(cert['cases'])==set(map(str,range(2,7)))
EXPECTED={2:(20,114),3:(30,1254),4:(40,1092),5:(30,72),6:(35,728)}
def validate_orbit(data,c):
 deg,period=EXPECTED[c]
 q=data['modulus'];vv=data['v_representative'];states=data['states'];X=[0,1]
 assert len(q)==deg+1 and q[-1]==1 and len(vv)<=deg
 assert data['period']==period and len(states)==period and states[0]==[1]
 assert all(isinstance(z,int) and 0<=z<7 for a in [q,vv]+states for z in a)
 assert all(a and a==tr(a[:]) and len(a)<=deg for a in states)
 assert gcd(q,derivative(q))==gcd(q,vv)==[1]
 if c>=5:assert gcd(q,X)==[1]
 state=[1]
 for j,st in enumerate(states):
  assert st==state
  assert gcd(sub(state,X),q)==([6,1] if c==6 and j%7==0 else [1])
  state=rem(sub([1],mul(vv,powmod(state,7,q))),q)
  check('complete_orbit_states_transitions_and_target_gcds',3)
 assert state==[1]
 if c==6:assert sum(vv)%7==6
 check('complete_orbit_schema_and_returns',4)
def validate_witnesses(case,gens):
 for w in case['necessary_witnesses']:
  assert len(w['multipliers'])==len(gens)
  z={}
  for a,g in zip(w['multipliers'],gens):z=ba(z,bm(decode(a),g))
  assert z==decode(w['polynomial'])
def r(f):return s.Poly(s.expand(f),S,v,modulus=7).as_expr()
def D(f):return r(v*s.diff(f,v)+(S-1)*s.diff(f,S))
x,y=s.symbols('x y');Uxy=x+y;Bxy=-x-S*y
Pxy=Uxy+3*Uxy**2+5*Uxy**3+5*Uxy**4+3*Uxy**5+Uxy**6
CT=s.expand((Uxy**7-x**7-y**7)/7)
CA=s.expand((Bxy**7-(-x)**7-(-S*y)**7)/7)
U=S-1;E=v*S**7+U
Qtrace=1-v**7+2*v*(U+v)**6+v**2*U**7*(U+v)**5
profiles=[];matrices={}
for c in range(2,7):
 case=cert['cases'][str(c)];h=min(c+1,6)
 assert case['offset']==c and case['divided_rows']==list(range(h))
 assert len(case['matrix'])==6 and all(len(row)==5 for row in case['matrix'])
 assert len(case['maximal_minors'])==6
 seeds=[]
 for a in range(6):
  if a<h:
   f=CT+Pxy+6*v*(1+Uxy)*Bxy**6 if a==0 else (
      (CT+Pxy)*Bxy**(a-1)*(Bxy-a*(1+Uxy))
       -a*v*(1+Uxy)*Bxy**(a-1)*CA)
   z=s.expand(f).coeff(x,c).coeff(y,2)
  else:z=s.expand((1+Uxy)*Bxy**a).coeff(x,c).coeff(y,2)
  seeds.append(r(z))
 assert all(bp(z)==decode(a) for z,a in zip(seeds,case['seed_polynomials']))
 check('uniform_all_depth_seed_extractions',6)
 rows=[[z] for z in seeds]
 for j in range(4):
  old=[row[-1] for row in rows]
  for a in range(6):
   cross=(v*old[5] if h<6 else 0) if a==0 else (0 if a==h else a*old[a-1])
   rows[a].append(r((c+5-a)*old[a]+D(old[a])+cross))
 M=[[bp(z) for z in row] for row in rows];matrices[c]=rows
 assert M==[[decode(z) for z in row] for row in case['matrix']]
 check('all_universal_matrix_entries',30)
 minors=[]
 for omit in range(6):
  z=det([row for a,row in enumerate(M) if a!=omit]);minors.append(z)
  assert z==decode(case['maximal_minors'][omit]);check('all_exact_maximal_minors')
 gens=[decode(z) for z in case['necessary_generators']]
 if c<6:assert gens==minors
 else:
  assert all(not z for z in minors[1:])
  assert minors[0]==bm({(0,5):1},gens[0])
  assert gens[1]==bp(Qtrace)
  check('actual_moment_generator_and_localized_minor',2)
 validate_witnesses(case,gens)
 check('explicit_necessary_ideal_witnesses',len(case['necessary_witnesses']))
 q=case['modulus'];vv=case['v_representative']
 ts=[decode(w['polynomial']) for w in case['necessary_witnesses']]
 if c<5:
  assert len(ts)==2 and ts[1]=={(0,i):z for i,z in enumerate(q) if z}
  assert ts[0]==ba({(1,0):1},bs({(0,i):z for i,z in enumerate(vv) if z},-1))
 else:
  exponent=5 if c==5 else 7
  linear=ts[1 if c==5 else 2];qq=ts[2 if c==5 else 3]
  assert qq=={(0,i+exponent):z for i,z in enumerate(q) if z}
  rest=linear.copy();assert rest.pop((1,exponent))==1
  assert all(i==0 for i,j in rest)
  rp=[0]*(max(j for i,j in rest)+1)
  for (i,j),z in rest.items():rp[j]=z
  assert rem(add(mul([0]*exponent+[1],vv),rp),q)==[0]
  at0={i:z for (i,j),z in (ts[0] if c==5 else gens[1]).items() if j==0}
  assert at0==({0:1,1:3,2:6,3:3,4:1} if c==5 else {0:1,1:2,2:3,3:4,4:5,5:6})
  check('localized_and_exceptional_rank_targets',3)
 validate_orbit(case,c)
 for rr in (1,2):
  n=2*7**rr+c;den=6**n*math.factorial(n);cols=rawcols(n)
  ss=1+sum((-1)**i*v**((7**i-1)//6) for i in range(1,rr+1))
  for j,raw in enumerate(cols):
   for a in range(6):
    actual=sum(modq(Q(math.factorial(a)*7**ell*raw[a+6*ell],den*(7 if a<h else 1)))*v**ell
       for ell in range(max(0,(len(raw)-1-a)//6+1)))
    assert s.Poly(actual-rows[a][j].subs(S,ss),v,modulus=7).is_zero
    check('literal_full_entry_polynomial_profiles')
  profiles.append({'offset':c,'r':rr,'terminal':n,'entries':30})
# Corrupted necessary and complete-cycle records must fail.
bad=copy.deepcopy(cert['cases']['2'])
bad['necessary_witnesses'][0]['multipliers'][0][0][2]=(bad['necessary_witnesses'][0]['multipliers'][0][0][2]+1)%7
try:validate_witnesses(bad,[decode(z) for z in bad['necessary_generators']])
except AssertionError:check('corrupted_necessary_witness_rejected')
else:raise AssertionError('bad witness accepted')
bad=copy.deepcopy(cert['cases']['6']);bad['states'][0]=[2]
try:validate_orbit(bad,6)
except AssertionError:check('corrupted_complete_orbit_rejected')
else:raise AssertionError('bad orbit accepted')
# Actual moment-root residue: independently reduce in its seven-dimensional algebra.
tau=s.symbols('tau')
G=v*sum((-U)**j*tau**(6-j) for j in range(7))+1
num=(v*tau**6+1)**6*G**2
rr=s.rem(num,tau**7+tau/v-1,tau,domain=s.GF(7).frac_field(S,v))
assert s.Poly(s.cancel(s.Poly(rr,tau).nth(6)/v-Qtrace),S,v,modulus=7).is_zero
assert s.Poly(E**2*(U+v)**5-v**7*S**7-U**7*Qtrace,S,v,modulus=7).is_zero
check('uniform_actual_reciprocal_trace_identities',2)
for rr in (1,2,3):
 n=2*7**rr+6
 ss=1+sum((-1)**i*v**((7**i-1)//6) for i in range(1,rr+1))
 actual=sum((math.comb(n+6*k-1,7*k-1)%7)*v**(n-k) for k in range(1,n+1))
 assert s.Poly(actual-Qtrace.subs(S,ss),v,modulus=7).is_zero
 check('literal_full_moment_residue_polynomials')
# Trace packets: exact arithmetic in (Z/49)[X]/P_Z, without p-adic approximations.
def zadd(a,b):
 c=[0]*max(len(a),len(b))
 for i,x in enumerate(a):c[i]+=x
 for i,x in enumerate(b):c[i]+=x
 return tr([x%49 for x in c])
def zsub(a,b):return zadd(a,[-x%49 for x in b])
def zmul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return tr([x%49 for x in c])
def zrem(a,q49):
 a=a[:]
 while len(a)>=len(q49):
  j=len(a)-len(q49);x=a[-1]
  for i,y in enumerate(q49):a[i+j]=(a[i+j]-x*y)%49
  tr(a)
 return a
def zpow(a,n,q49):
 out=[1]
 while n:
  if n%2:out=zrem(zmul(out,a),q49)
  a=zrem(zmul(a,a),q49);n//=2
 return out
def smalltraces(q49,n):
 a=[7]
 for k in range(1,n+1):
  z=(-sum(q49[7-j]*a[k-j] for j in range(1,k))-k*q49[7-k] if k<=7
    else -sum(q49[7-j]*a[k-j] for j in range(1,8)))
  a.append(z%49)
 return a
X=[0,1]
minus=cert['minus_one_packets'];q15=minus['trace_polynomial']
assert q15==[48,7,28,35,14,21,27,1]
assert zpow(X,343,q15)==zpow(X,7,q15)
trace15=smalltraces(q15,104)
for c in (5,6):
 expected=minus['point_values'][str(c)]
 assert trace15[14+c]==expected['odd'] and trace15[98+c]==expected['even']
 # Separate binomial formula for the same four exact representatives.
 for n,want in ((14+c,expected['odd']),(98+c,expected['even'])):
  value=(7+sum((7*math.comb(n+6*k-1,7*k)+math.comb(n+6*k-1,7*k-1))*pow(15,k,49)
     for k in range(1,n+1)))%49
  assert value==want
  check('two_algorithm_fixed_minus_one_packet_values',2)
plus=cert['plus_one_packet'];q34=plus['trace_polynomial']
assert q34==[48,7,28,35,14,21,8,1]
assert plus['period_identity_exponents']==[7**8,7]
assert zpow(X,7**8,q34)==zpow(X,7,q34)
remainder=zpow(X,plus['reference_exponent'],q34)
assert plus['reference_exponent']==2*7**7+6
assert remainder==plus['reference_remainder']==[20,29,21,7,28,7,15]
small=smalltraces(q34,6)
assert small==plus['small_power_traces']==[7,41,22,48,8,34,15]
assert sum(a*b for a,b in zip(remainder,small))%49==plus['trace_value']==21
check('uniform_packet_period_and_reference_identities',6)
imports=['GAP_SIX_TWO_POWER_PLUS_ONE_SFC.tex',
 'ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex','PHYSICAL_PRIME_TOTAL_DEGREE_REDUCTION.tex',
 'gap_three_hurwitz_structure.tex','GAP_SIX_SCALED_POWER_PLUS_ONE_SFC.tex',
 'GAP_SIX_SMALL_OFFSET_FAMILIES.tex']
files=[HERE/'GAP_SIX_TWO_DIGIT_LEADING_TWO_SFC.tex',Path(__file__).resolve(),CERT_PATH]+[HERE/z for z in imports]
out={'passed':True,
 'scope':'ordinary support{1,7}, all terminal N=2*7^r+c,c0..6,r>=0; uniform formal proof plus complete finite certificates',
 'all_index_gap_six_SFC_claimed':False,'Gamma_shape_extension_claimed':False,
 'controls':COUNTS,'literal_profiles':profiles,
 'source_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in files},
 'python':platform.python_version(),'sympy':s.__version__,
 'independent_replay':False,'formal_verification':False,'external_review':False,
 'novelty_claim':False,'priority_claim':False,
 'command':'.venv/bin/python research/factorial/image-rigidity/gap_six_two_digit_leading_two_sfc.py'}
(HERE/'gap_six_two_digit_leading_two_sfc_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'passed':True,'controls':COUNTS,'total_controls':sum(COUNTS.values())},indent=2))
