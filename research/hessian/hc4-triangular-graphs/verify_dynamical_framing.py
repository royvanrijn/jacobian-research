"""First preregistration transcription controls for MG0; refuse registered replay."""
import hashlib,json,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
ctxpath=HERE/'sources/dynamical-framing-context.json'
receipt=HERE/'receipt-dynamical-framing.json'
known={e['id']:e for e in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']}
assert 'HC4MYMG0' not in known,'registered mathematical checker replay forbidden'
assert not receipt.exists(),'receipt already exists; do not replay'
ctx=json.loads(ctxpath.read_text());checks=[];algebraic=0;start=time.perf_counter()
def check(name,condition,algebra=False):
 global algebraic
 assert bool(condition),name
 checks.append(name);algebraic+=int(algebra)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
for p,d in ctx['retained_proof_and_source_sha256'].items():
 check('sealed '+p,sha(ROOT/p)==d)
pdf=ROOT/ctx['primary_source']['local_path']
check('published PDF pin',sha(pdf)==ctx['primary_source']['sha256'])
S,T,p,q,b=s.symbols('S T p q b')
M=s.symbols('M',integer=True,positive=True)
H=(p+q**M,b*(p+q**M)+q)
I=(S-(T-b*S)**M,T-b*S)
check('symbolic all-M target determinant one',s.simplify(s.det(s.Matrix(H).jacobian([p,q]))-1)==0,True)
back=[z.subs({S:H[0],T:H[1]},simultaneous=True) for z in I]
forward=[z.subs({p:I[0],q:I[1]},simultaneous=True) for z in H]
for j,z in enumerate([s.simplify(back[0]-p),s.simplify(back[1]-q),s.simplify(forward[0]-S),s.simplify(forward[1]-T)]):
 check('symbolic two-sided inverse component '+str(j),z==0,True)
check('literal old second generator recovered',s.expand(H[1]-b*H[0]-q)==0,True)
check('literal old first generator recovered',s.simplify(H[0]-(H[1]-b*H[0])**M-p)==0,True)
x,y,z=s.symbols('x y z')
forms=[(x+y,1),(x*x+x*y+y*y,1),(x**3-y**3,2),((x+2*y)**4,1)]
for j,(Q,b0) in enumerate(forms):
 degree=s.Poly(Q,x,y).total_degree();value=Q.subs({x:1,y:b0})
 check('nonzero frame direction '+str(j),value!=0,True)
 check('homogeneous ray '+str(j),s.expand(Q.subs({x:z,y:b0*z},simultaneous=True)-value*z**degree)==0,True)
 for m in [2,3]:
  D=m*degree;h=s.expand(Q**m);a=value**m
  check('framed leading degree '+str((j,m)),s.Poly(h,x,y).total_degree()==D,True)
  check('nonvanishing ray of full leading map '+str((j,m)),s.expand(h.subs({x:z,y:b0*z},simultaneous=True)-a*z**D)==0 and a!=0,True)
a=s.symbols('a',nonzero=True)
for D in [2,3,5]:
 c=s.Integer(1)
 for n in range(1,5):
  expected=a**sum(D**j for j in range(n-1))
  check('leading scalar recurrence '+str((D,n)),s.simplify(c-expected)==0,True)
  c=a*c**D
u,t=s.symbols('u t')
e,be=s.symbols('e be',integer=True,positive=True)
jac=s.det(s.Matrix([u**e,t]).jacobian([u,t]))
check('symbolic tame transverse differential',s.simplify(jac-e*u**(e-1))==0,True)
check('thinness coefficient simplification',s.simplify((1+e-1)/be-e/be)==0,True)
check('strict positive thinness',(e/be).is_positive,True)
expected_true={k for k,v in ctx['scope_flags'].items() if v}
check('scope flags present',len(expected_true)>=8)
check('no geometric theorem certified by samples',ctx['external_imports'][0]['verification']=='published theorem import;not certified by this checker')
out=dict(outcome='PASS',check_count=len(checks),algebraic_controls=algebraic,checks=checks,elapsed_seconds=time.perf_counter()-start,proof_sha256=sha(HERE/'INTEGER_DYNAMICAL_DEGREE_AND_POSITIVE_THINNESS.md'),checker_sha256=sha(Path(__file__)),context_sha256=sha(ctxpath),primary_pdf_sha256=sha(pdf),registered_mathematical_checkers_replayed=[],independent_replay=False,formal_verification=False,external_review=False)
receipt.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['outcome','check_count','algebraic_controls','elapsed_seconds']}))
