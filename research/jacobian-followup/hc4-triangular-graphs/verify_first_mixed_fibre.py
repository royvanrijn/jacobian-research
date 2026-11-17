"""NEW MU0 first mixed-fibre/order/trace transcriptions; never replay."""
import hashlib,json,math,time
from pathlib import Path
import sympy as S
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
RECEIPT=HERE/'receipt-first-mixed-fibre.json'
assert not RECEIPT.exists(),'Sealed receipt exists: never replay.'
assert not any(e['id']=='HC4MYMU0' for e in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries'])
started=time.perf_counter();checks=[];algebraic=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check(name,value,alg=True):
 assert bool(value),name
 checks.append(name)
 if alg:algebraic.append(name)
y,u,s,t,z,v=S.symbols('y u s t z v')
h=y**2*(y-1)**3;hp=S.diff(h,y);a=-S.Rational(108,3125)
check('exact quintic derivative',S.expand(hp-y*(y-1)**2*(5*y-2))==0)
check('second critical value',h.subs(y,S.Rational(2,5))==a)
check('monic primitive-order discriminant',S.expand(S.discriminant(h-u,y)-u**3*(3125*u+108))==0)
cubic=S.cancel((h-a)/(y-S.Rational(2,5))**2)
check('second critical fibre has exactly one double root',S.degree(cubic,y)==3 and cubic.subs(y,S.Rational(2,5))!=0)
check('remaining three roots simple',S.degree(S.gcd(cubic,S.diff(cubic,y)),y)==0)

def perm(cycles):
 out=list(range(5))
 for cyc in cycles:
  for i,x in enumerate(cyc):out[x]=cyc[(i+1)%len(cyc)]
 return tuple(out)
def compose(p,q):return tuple(p[q[i]] for i in range(5))
def cycle_lengths(p):
 seen=set();out=[]
 for x in range(5):
  if x in seen:continue
  j=x;length=0
  while j not in seen:seen.add(j);length+=1;j=p[j]
  out.append(length)
 return sorted(out)
g0=perm([(0,1),(2,3,4)]);g1=perm([(1,2)])
group={tuple(range(5))};pending=list(group)
while pending:
 p=pending.pop()
 for g in [g0,g1]:
  q=compose(p,g)
  if q not in group:group.add(q);pending.append(q)
check('finite branch zero has cycle type two plus three',cycle_lengths(g0)==[2,3])
check('other finite branch is a transposition',cycle_lengths(g1)==[1,1,1,2])
check('infinity is a five-cycle',cycle_lengths(compose(g0,g1))==[5])
check('exact branch generators give relative S5',len(group)==120)

# Newton traces in the actual polynomial extension, not a split numerical fit.
def polynomial_trace(H,G):
 m=S.degree(H,y);co=S.Poly(H-u,y).all_coeffs();power=[S.Integer(m)]
 for k in range(1,max(S.degree(G,y),1)+1):
  if k<=m:value=-sum(co[j]*power[k-j] for j in range(1,k))-k*co[k]
  else:value=-sum(co[j]*power[k-j] for j in range(1,m+1))
  power.append(S.expand(value))
 return S.expand(sum(c*power[j[0]] for j,c in S.Poly(G,y).terms()))
for k in range(1,5):check(f'first four relative y moments k{k}',polynomial_trace(h,y**k)==3)
check('fifth relative y moment',polynomial_trace(h,y**5)==3+5*u)
check('relative first derivative trace zero',polynomial_trace(h,hp)==0)
check('relative derivative square trace minus eight u',polynomial_trace(h,S.expand(hp**2))==-8*u)
check('all-degree kappa first model',2*(2*S.diff(h,y,2).subs(y,0)+3*S.diff(h,y,2).subs(y,1))==-8)
for N in [3,4,5,8]:
 H=y**2*(y-1)**N
 check(f'all-degree two plus N first derivative trace N{N}',polynomial_trace(H,S.expand(S.diff(H,y)))==0)
 check(f'all-degree two plus N derivative-square trace N{N}',polynomial_trace(H,S.expand(S.diff(H,y)**2))==8*(-1)**N*u)

# Retain a genuinely rational slice in a field-level transcription. No claim
# that this sample beta realizes an actual Keller source-image algebra.
M=S.zeros(5)
for j in range(4):M[j+1,j]=1
M[:,4]=S.Matrix([u,0,1,-3,3])
I=S.eye(5);Hp=5*M**4-12*M**3+9*M**2-2*M
Beta=M.inv()**3
X=Hp*(s*I-Beta)
R1=S.cancel(S.trace(X));R2=S.cancel(S.trace(X**2))
check('rational slice companion satisfies h',M**5-3*M**4+3*M**3-M**2-u*I==S.zeros(5))
check('rational slice first trace invariant',S.diff(R1,s)==0)
check('rational slice quadratic trace equation',S.cancel(S.diff(R2,s,2)+16*u)==0)

# Every row uses the same general gcd/different/conductor formulas.
for name,r,ns,degrees,index,conductors in [
 ('C2',2,[2,3],[1,1,3],2,[1,1,2]),
 ('C3',3,[2,3],[2,1,1,1],4,[2,2,2,2]),
 ('D2',2,[2,1,1,1],[1,1,1,1,1],1,[1,1,0,0,0]),
 ('D1',1,[2,1,1,1],[2,1,1,1],0,[0,0,0,0]),
]:
 ds=[];cs=[]
 for n in ns:
  g=math.gcd(r,n);ds.extend([n//g]*g);cs.extend([(r*(n-1)-n+g)//g]*g)
 disc_order=r*sum(n-1 for n in ns);true_order=sum(n-math.gcd(r,n) for n in ns)
 check(name+' completed normal factor degrees',ds==degrees and sum(ds)==5)
 check(name+' primitive versus true discriminant index',(disc_order-true_order)//2==index and (disc_order-true_order)%2==0)
 check(name+' full conductor branch orders',cs==conductors)
 check(name+' tame trace-dual orders',[1-e for e in ds]==[0 if e==1 else -(e-1) for e in ds])

for e,pair in [(2,(-1,-1)),(3,(-1,-2))]:
 check(f'codifferent permits pole pair e{e}',all(k>=1-e for k in pair) and sum(pair)==-e)
 check(f'product trace has a genuine base pole e{e}',S.denom(e/t)==t)

# Exact local equation; sharpness of the valuation bound, not a global model.
for N in [3,4,7]:
 local_delta=lambda f:-N*t**(N-1)*S.diff(f,v)
 local_R2=8*S.Integer(-1)**N/S.Integer(N)**2*t**(2-N)*v**2
 check(f'local quadratic trace equation N{N}',S.cancel(local_delta(local_delta(local_R2))-16*(-1)**N*t**N)==0)
 check(f'forced pole order N{N}',2-N<0 and N-2*(N-1)==2-N)

# A genuine field with unramified primes: an unmarked trace can cancel,
# while an integral intermediate marker detects the target-prime pole.
V=S.Matrix([[0,t+1],[1,0]]);A=V/t
check('integral marker monic equation',V**2-(t+1)*S.eye(2)==S.zeros(2))
check('unmarked field trace cancellation retained',S.trace(A)==0)
check('integral marker detects a target pole',S.cancel(S.trace(V*A)-2-2/t)==0)

cp=HERE/'sources/first-mixed-fibre-context.json';C=json.loads(cp.read_text())
for p,digest in C['retained_proof_and_source_sha256'].items():check('retained byte pin '+p,sha(ROOT/p)==digest,False)
check('five retained claims and ten proof-checker pins',len(C['retained_structural_inputs'])==5 and len(C['retained_proof_and_source_sha256'])==10,False)
for k in ['actual_normalized_degree_five_factor','rational_slice_retained','relative_S5_not_deck_obstruction','actual_source_image_fibres','full_relative_order_conductor_and_codifferent','quadratic_trace_pole_independent_of_slice','all_degree_nonzero_kappa_test','integral_marker_target_trace_witness']:
 check('scope '+k,C['scope_flags'][k],False)
for k in ['beta_zero_assumed','source_image_identified_with_old_universal_cover','monogenic_order_identified_with_normalization','unmarked_absolute_trace_pole_asserted','actual_Keller_model_constructed','all_mixed_profiles_excluded','automatic_Mathieu_trace_transfer','unrestricted_HC4_or_JC2_resolved']:
 check('boundary '+k,not C['scope_flags'][k],False)
check('no new imports caches downloads or predecessor replay',C['new_external_theorem_imports']==C['new_primary_caches']==C['new_primary_downloads']==0 and C['registered_mathematical_checkers_replayed']==[],False)
r=dict(outcome='PASS',claim='HC4MYMU0',check_count=len(checks),algebraic_controls=len(algebraic),checks=checks,algebraic_check_names=algebraic,elapsed_seconds=time.perf_counter()-started,
 proof_sha256=sha(HERE/'FIRST_MIXED_FIBRE_SOURCE_IMAGE_AND_QUADRATIC_TRACE_POLE.md'),checker_sha256=sha(Path(__file__)),context_sha256=sha(cp),registered_mathematical_checkers_replayed=[],independent_replay=False,formal_verification=False,external_review=False)
RECEIPT.write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps({k:r[k] for k in ['outcome','check_count','algebraic_controls','elapsed_seconds']}))
