"""NEW preregistration controls for MM0; refuse to overwrite a sealed receipt."""
import hashlib,json,time
from pathlib import Path
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
RECEIPT=HERE/'receipt-abhyankar-gurjar-framing.json'
assert not RECEIPT.exists(),'Sealed receipt exists: registered mathematical checker must not be replayed.'
started=time.perf_counter();checks=[];algebraic=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check(name,condition,alg=True):
 assert bool(condition),name
 checks.append(name)
 if alg:algebraic.append(name)
x,y,z1,z2,w1,w2,c1,c2,t=s.symbols('x y z1 z2 w1 w2 c1 c2 t')
def expand(p):return s.expand(p)
def top(p):
 p=s.Poly(expand(p),x,y);d=p.total_degree()
 return expand(sum(v*x**a*y**b for (a,b),v in p.terms() if a+b==d))
def trunc(p,d):return expand(sum(v*x**a*y**b for (a,b),v in s.Poly(expand(p),x,y).terms() if a+b<=d))
def T(a,H,m):
 return expand(sum(s.diff(a*H[0]**i*H[1]**(m-i),x,i,y,m-i)/(s.factorial(i)*s.factorial(m-i)) for i in range(m+1)))
def wick(p):
 out=0
 for (a,b,i,j),coef in s.Poly(expand(p),w1,w2,z1,z2).terms():
  if (a,b)==(i,j):out+=coef*s.factorial(a)*s.factorial(b)
 return expand(out)
F=(x+y**2,y+(x+y**2)**2);H=(x-F[0],y-F[1]);phi=y**4;D=4
Inv=(x-(y-x**2)**2,y-x**2)
for name,A,C in [('forward-inverse',F,Inv),('inverse-forward',Inv,F)]:
 check(name,all(expand(a.subs({x:C[0],y:C[1]},simultaneous=True)-b)==0 for a,b in zip(A,(x,y))))
check('unit Jacobian actual plane control',s.det(s.Matrix(F).jacobian([x,y]))==1)
check('identity first jet actual plane control',s.Matrix(F).jacobian([x,y]).subs({x:0,y:0})==s.eye(2))
for a,label in [(s.Integer(1),'one'),(x,'first coordinate vanishing at e2'),(y,'second coordinate'),(x**2+3*x*y+y**3,'mixed amplitude')]:
 for m in range(1,4):
  b=T(a,H,m);expected=expand((-1)**m*s.diff(top(a)*phi**m,y,m)/s.factorial(m))
  check(f'all-amplitude top transcription {label} m{m}',top(b)==expected and s.Poly(b,x,y).total_degree()==s.Poly(a,x,y).total_degree()+(D-1)*m)
for m in range(1,6):
 expected=(-1)**m*s.factorial(D*m)/(s.factorial(m)*s.factorial((D-1)*m))
 check(f'nonzero scalar top coefficient m{m}',top(T(s.Integer(1),H,m)).subs({x:0,y:1})==expected and expected!=0)
cut=5
for a,label in [(s.Integer(1),'one'),(x,'x inverse'),(y,'y inverse'),(F[0],'first target pullback'),(F[1],'second target pullback')]:
 lhs=sum(T(a,H,m) for m in range(cut+1))
 rhs=a.subs({x:Inv[0],y:Inv[1]},simultaneous=True)
 check(f'formal inverse sum through source degree {cut}: {label}',trunc(lhs-rhs,cut)==0)
P=expand(w1*H[0].subs({x:c1+z1,y:c2+z2})+w2*H[1].subs({x:c1+z1,y:c2+z2}))
for a,label in [(s.Integer(1),'one'),(x,'x'),(y,'y')]:
 for m in [1,2]:
  lhs=wick(a.subs({x:c1+z1,y:c2+z2})*P**m)
  rhs=s.factorial(m)*T(a,H,m).subs({x:c1,y:c2})
  check(f'translated Wick contraction {label} m{m}',expand(lhs-rhs)==0)
for a,b in [(0,0),(1,0),(0,2),(2,1)]:
 check(f'circular Wick normalization {a},{b}',wick(w1**a*w2**b*z1**a*z2**b)==s.factorial(a)*s.factorial(b))
check('unbalanced Wick monomial vanishes',wick(w1*z2)==0)
P0=w1*z1**2+w2*z2**3
check('unshifted non-Keller phase pure moments zero',all(wick(P0**m)==0 for m in [1,2,3]))
check('unshifted non-Keller phase fixed markers safe',all(wick(a*P0**3)==0 for a in [1,z1,w1,w1**2,w1*w2,z1*w2]))
check('degree-imbalance control really not Keller',s.det(s.Matrix([x-x**2,y-y**3]).jacobian([x,y]))!=1)
tr=s.diff(H[0],x)+s.diff(H[1],y)
check('mixed-degree correction divergence nonzero',tr!=0 and expand(tr-T(s.Integer(1),H,1))==0)
check('linear tau homotopy not identically Keller',expand(s.det((s.Matrix([x,y])-t*s.Matrix(H)).jacobian([x,y]))-1-t*(t-1)*tr)==0 and tr!=0)
def diagonal(p):
 out=0
 for (a,b,i,j),coef in s.Poly(expand(p),w1,w2,z1,z2).terms():
  if (a,b)==(i,j):out+=coef*x**a*y**b
 return expand(out)
p=w1*z2+w2*z1
check('diagonal projection is not multiplicative',diagonal(p)==0 and diagonal(p**2)==2*x*y)
context=HERE/'sources/abhyankar-gurjar-framing-context.json';C=json.loads(context.read_text())
for p,digest in C['retained_proof_and_source_sha256'].items():check('retained byte pin '+p,sha(ROOT/p)==digest,False)
check('three retained structural inputs',len(C['retained_structural_inputs'])==3,False)
check('no new imports or source caches',C['external_imports']==C['new_primary_caches']==[],False)
check('full original target retained',C['scope_flags']['full_goal_active'] and not C['scope_flags']['unrestricted_HC4_or_JC2_resolved'],False)
check('no old mathematical checker invoked',C['registered_mathematical_checkers_replayed']==[],False)
receipt=dict(outcome='PASS',check_count=len(checks),algebraic_controls=len(algebraic),checks=checks,algebraic_check_names=algebraic,elapsed_seconds=time.perf_counter()-started,proof_sha256=sha(HERE/'UNIVERSAL_ABHYANKAR_GURJAR_NONTERMINATION_AND_MOMENT_TRANSFER.md'),checker_sha256=sha(Path(__file__)),context_sha256=sha(context),sympy_version=s.__version__,registered_mathematical_checkers_replayed=[],independent_replay=False,formal_verification=False,external_review=False)
RECEIPT.write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ['outcome','check_count','algebraic_controls','elapsed_seconds']}))
