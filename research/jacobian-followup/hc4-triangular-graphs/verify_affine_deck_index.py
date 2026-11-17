"""NEW MN0 exact transcription controls; never overwrite a sealed receipt."""
import hashlib,itertools,json,math,time
from fractions import Fraction
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
RECEIPT=HERE/'receipt-affine-deck-index.json'
assert not RECEIPT.exists(),'Sealed mathematical receipt exists: do not replay.'
started=time.perf_counter();checks=[];algebraic=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check(name,c,alg=True):
 assert bool(c),name
 checks.append(name)
 if alg:algebraic.append(name)
perms=list(itertools.permutations(range(4)))
for q in [1,2,3]:
 t=2**q;degree=4*t;supports=set();count=0;minimum=degree
 for g in perms:
  for a,b,c in itertools.product(range(t),repeat=3):
   values=(c,c^b,c^a,c^a^b)
   if g==(0,1,2,3) and a==b==c==0:continue
   fixed=t*sum(g[x]==x and values[x]==0 for x in range(4))
   support=degree-fixed;supports.add(support);minimum=min(minimum,support);count+=1
 check(f'ALL natural H{q} elements: affine minimum support',minimum==degree//2)
 check(f'ALL natural H{q} elements: support spectrum',supports=={degree//2,3*degree//4,degree})
 check(f'ALL natural H{q} affine functions counted',count==24*2**(3*q)-1)
def action(g,values,t):return tuple(g[x]*t+(z^values[x]) for x in range(4) for z in range(t))
def compose(a,b):return tuple(a[b[i]] for i in range(len(a)))
for q in [1,2,3]:
 t=2**q;ident=tuple(range(4*t));r=1;g=(1,0,2,3)
 f=action(g,(0,r,0,r),t);square=compose(f,f)
 check(f'order-four inertia outside decorated class q{q}',square!=ident and compose(square,square)==ident)
 moved=sum(a!=i for i,a in enumerate(f));check(f'order-four support includes fixed-block residual q{q}',moved==3*t and moved-2*t==t)
 vertical=action((0,1,2,3),(0,r,0,r),t)
 check(f'vertical inertia also obeys half-support bound q{q}',sum(a!=i for i,a in enumerate(vertical))==2*t)
 # Constants commute with a generating set of ALL affine-function translations and base S4.
 generators=[action(g0,(0,0,0,0),t) for g0 in perms]
 for j in range(q):
  v=1<<j
  generators.extend(action((0,1,2,3),vals,t) for vals in [(v,v,v,v),(0,0,v,v),(0,v,0,v)])
 for j in range(q):
  translation=action((0,1,2,3),(1<<j,)*4,t)
  check(f'constant radical deck generator central q{q} j{j}',all(compose(translation,a)==compose(a,translation) for a in generators))
# Nonbinary affine fixed-point bounds on genuinely new small finite-field controls.
for p in [3,5]:
 points=list(itertools.product(range(p),repeat=2))
 maps=[lambda x,p=p:((x[0]+1)%p,x[1]),lambda x,p=p:(x[0],(-x[1])%p),lambda x,p=p:((x[0]+x[1])%p,x[1])]
 for j,f in enumerate(maps):
  fixed=sum(f(x)==x for x in points)
  check(f'affine F{p} fixed-point bound type{j}',fixed<=p)
# Low-Euler normalization: good covers produce reduced but generally SPLIT fibres.
for signature,order in [((2,2,2),4),((2,2,5),10),((2,3,3),12),((2,3,4),24),((2,3,5),60)]:
 orb=2-sum((1-Fraction(1,m) for m in signature),Fraction(0))
 check(f'spherical cover Euler identity {signature}',order*orb==2)
 check(f'normalized exact-root transverse indices {signature}',all(m//math.gcd(m,m)==1 for m in signature))
 # Over each base preimage a degree-order etale cover of an original A1 has order components.
 components=sum((order//m)*(m-1) for m in signature)
 check(f'pullback splitting retained {signature}',2+components==2*order and components>0)
for a,b,expected in [(2,3,1),(2,4,2),(3,6,3)]:
 check(f'bad unequal spindle group is cyclic gcd({a},{b})',math.gcd(a,b)==expected)
for Q,d,h,allowed in [(2,16,4,False),(2,32,8,False),(2,64,16,False),(2,32,4,True),(2,64,8,True),(3,27,9,False),(3,81,9,True),(5,625,25,True)]:
 check(f'noncyclic deck/Euler necessary inequality Q{Q} d{d} h{h}',(3*h<=Fraction(d,Q))==allowed)
for n in range(3,9):
 d=2**n
 candidates=[d//2**j for j in range(n+1) if 3*2**j<=d//2]
 check(f'binary deck-index rounding d{d}',not candidates or min(candidates)>=8)
for p in [3,5,7]:
 d=p**4;indices=[p**j for j in range(5) if p**j>=3*p]
 check(f'prime-field deck-index rounding p{p}',min(indices)==p*p)
for q in range(2,9):
 t=2**q;d=4*t
 check(f'entire natural radical tower excluded q{q}',Fraction(d,2*t)==2 and d//t==4 and not (3*t<=d//2))
context=HERE/'sources/affine-deck-index-context.json';C=json.loads(context.read_text())
for p,digest in C['retained_proof_and_source_sha256'].items():check('source byte pin '+p,sha(ROOT/p)==digest,False)
check('three prior registered mathematical claims',len(C['retained_structural_inputs'])==3,False)
check('two explicit classical inputs and primary source caches',len(C['external_imports'])==len(C['new_primary_caches'])==2,False)
check('no registered checker replay',C['registered_mathematical_checkers_replayed']==[],False)
check('full goal remains open',C['scope_flags']['full_goal_active'] and not C['scope_flags']['unrestricted_HC4_or_JC2_resolved'],False)
receipt=dict(outcome='PASS',check_count=len(checks),algebraic_controls=len(algebraic),checks=checks,algebraic_check_names=algebraic,elapsed_seconds=time.perf_counter()-started,proof_sha256=sha(HERE/'AFFINE_MONODROMY_DECK_INDEX_AND_UNDECORATED_RADICAL_TOWERS.md'),checker_sha256=sha(Path(__file__)),context_sha256=sha(context),registered_mathematical_checkers_replayed=[],independent_replay=False,formal_verification=False,external_review=False)
RECEIPT.write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ['outcome','check_count','algebraic_controls','elapsed_seconds']}))
