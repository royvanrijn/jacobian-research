"""NEW preregistration controls for MH0; refuse registered replay."""
from pathlib import Path
import hashlib,json,time
import sympy as s
R=Path(__file__).resolve().parents[2];H=Path(__file__).resolve().parent
receipt=H/'receipt-iterated-deficit.json';cp=H/'sources/iterated-deficit-context.json'
known={e['id']:e for e in json.loads((R/'MATH_STATUS.json').read_text())['entries']}
assert 'HC4MYMH0' not in known,'registered checker replay forbidden'
assert not receipt.exists(),'receipt already exists;do not replay'
c=json.loads(cp.read_text());start=time.perf_counter();checks=[];alg=0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check(name,v,a=False):
 global alg
 assert bool(v),name
 checks.append(name);alg+=int(a)
for p,d in c['retained_proof_and_source_sha256'].items():check('sealed '+p,sha(R/p)==d)
# Finite-function controls test the universal algebra;not plane Keller examples.
maps=[([0,0,1],2),([1,2,0,0],2),([0,1,1,1,3],3),([1,1,3,4,4,5],3)]
for k,(f,d) in enumerate(maps):
 N=len(f)
 def push(v):
  out=[0]*N
  for x,y in enumerate(f):out[y]+=v[x]
  return out
 one=[1]*N;r=one[:];delta1=[d-z for z in push(one)]
 deficit=[0]*N;images=set(range(N));support=set();powerdelta=delta1[:]
 for n in range(1,6):
  old=deficit[:];r=push(r);deficit=[d**n-z for z in r]
  rec=[d**(n-1)*a+b for a,b in zip(delta1,push(old))]
  total=[0]*N;term=delta1[:]
  for j in range(n):
   total=[a+d**(n-1-j)*b for a,b in zip(total,term)];term=push(term)
  check('pointwise recurrence '+str((k,n)),deficit==rec,True)
  check('whole convolution '+str((k,n)),deficit==total,True)
  check('Euler/counting Fubini '+str((k,n)),sum(deficit)==N*(d**n-1),True)
  check('positivity and support '+str((k,n)),min(deficit)>=0 and {i for i,z in enumerate(deficit) if z>0}==support.union({i for i,z in enumerate(powerdelta) if z>0}),True)
  support={i for i,z in enumerate(deficit) if z>0};powerdelta=push(powerdelta)
t,x,y=s.symbols('t x y')
for m in [1,2,3,5,8]:
 xx=t*t;yy=t**(2*m+1)-t;eq=y*y-x*(x**m-1)**2
 check('nodal normalization equation '+str(m),s.expand(eq.subs({x:xx,y:yy},simultaneous=True))==0,True)
 check('no affine parameter critical point '+str(m),s.gcd(s.diff(xx,t),s.diff(yy,t))==1,True)
 check('node tangent determinant '+str(m),s.det(s.Matrix([[2*t,2*m],[-2*t,2*m]]))==8*m*t,True)
 check('generic normalization parameter recovered '+str(m),s.cancel(yy/(xx**m-1)-t)==0,True)
u,v,w=s.symbols('u v w')
e=s.symbols('e',integer=True,positive=True)
m,n=s.symbols('m n',integer=True,positive=True)
A=u+v**m;B=v+A**n
J=s.simplify(s.det(s.Matrix([A,B]).jacobian([u,v])))
check('all-parameter post-map is etale',J==1,True)
comp=[A.subs(u,w**e),B.subs(u,w**e)]
check('all-parameter finite-center label unchanged',s.simplify(s.det(s.Matrix(comp).jacobian([w,v]))-e*w**(e-1))==0,True)
d,z=s.symbols('d z',integer=True,positive=True)
for n0 in range(1,8):
 check('all-degree geometric Euler sum '+str(n0),s.expand((d-1)*sum(d**j for j in range(n0))-(d**n0-1))==0,True)
check('no external import or source cache',c['external_imports']==[] and c['new_primary_caches']==[])
out=dict(outcome='PASS',check_count=len(checks),algebraic_controls=alg,elapsed_seconds=time.perf_counter()-start,checks=checks,proof_sha256=sha(H/'ITERATED_DEFICIT_TRANSPORT_AND_FINITE_CENTER_LABELS.md'),checker_sha256=sha(Path(__file__)),context_sha256=sha(cp),registered_mathematical_checkers_replayed=[],independent_replay=False,formal_verification=False,external_review=False)
receipt.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['outcome','check_count','algebraic_controls','elapsed_seconds']}))
