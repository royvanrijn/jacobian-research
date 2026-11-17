"""Exact transcription controls for PROOF_DRAFT.tex; not a universal proof."""
from pathlib import Path
import hashlib,json
from math import factorial,prod
import sympy as s
BASE=Path(__file__).resolve().parent
RESEARCH=BASE.parent.parent
REG=json.loads((RESEARCH/'MATH_STATUS.json').read_text())['entries']
REG={x['id']:x for x in REG}
AUDIT=json.loads((BASE/'audit.json').read_text())
checks=[]
def check(label,condition):
    if not condition:raise AssertionError(label)
    checks.append(label)
inputs={}
for deduction in AUDIT['deductions']:
    for claim_id in deduction['inputs']:
        claim=REG[claim_id]
        check('registered proved input '+claim_id,claim['state']=='proved')
        p=RESEARCH/claim['canonical_source']
        check('canonical source exists '+claim_id,p.is_file())
        inputs[claim_id]={'source':claim['canonical_source'],
            'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
# The degree controls do not claim exceptional invariance; invariant-degree
# lower bounds are supplied by the retained exact character theorem.
x,y,z,w=s.symbols('x y z w')
gb=s.groebner([x*w-y*z-1],x,y,z,w,order='grlex')
def normal(p):return s.Poly(gb.reduce(s.expand(p))[1],x,y,z,w)
for H in [2,8,12]:
    h=x**H+w**H+(x*w+y*z)**(H//2)
    for k in [1,2,3]:
        check(f'filtered degree product H{H} power{k}',normal(h**k).total_degree()==H*k)
    check(f'filtered composition H{H}',normal(2*h**3+3*h+1).total_degree()==3*H)
# Quadratic image collision identity and a cubic loop control. The latter
# prevents incorrectly excluding degree three by the arc argument.
t,u,A,B,C=s.symbols('t u A B C')
q=lambda v:A*v*v+B*v+C
check('quadratic collision factor',s.expand(q(t)-q(u)-(t-u)*(A*(t+u)+B))==0)
p=lambda v:v*v-s.Rational(1,2)+s.I*(v**3-v)
check('cubic loop joins endpoints',p(-1)==p(1))
check('cubic loop is nonconstant',p(0)!=p(1))
check('cubic loop regular at real-derivative zero',s.diff(p(t),t).subs(t,0)!=0)
# Positivity/nondegeneracy control on a real polynomial Gaussian seed.
mu=lambda k:s.factorial2(2*k-1) if k else s.Integer(1)
G=s.Matrix(3,3,lambda i,j:mu(i+j))
for k in [1,2,3]:check(f'real-seed Gaussian Gram minor{k}',G[:k,:k].det()>0)
# Check the radial constants against an independent exact SL2 Haar formula.
a,b,c,d=s.symbols('a b c d');r=x*w-y*z
phase=x**4+y*z+2*x*w+3
lift=x**4+r*y*z+2*r*x*w+3*r**2
marker=x*z+2;marker_lift=x*z+2*r
sub={x:a+s.I*b,y:c+s.I*d,z:-c+s.I*d,w:a-s.I*b}
L=s.expand(lift.subs(sub));M=s.expand(marker_lift.subs(sub))
def haar(poly):
    ans=0
    for (i,j,k,l),v in s.Poly(s.expand(poly),x,y,z,w).terms():
        if i==l and j==k:
            ans+=v*(-1)**j*s.Rational(factorial(i)*factorial(j),factorial(i+j+1))
    return s.simplify(ans)
def gauss(poly):
    ans=0
    for ex,v in s.Poly(s.expand(poly),a,b,c,d).terms():
        if all(e%2==0 for e in ex):
            ans+=v*prod(s.factorial2(e-1) if e else 1 for e in ex)
    return s.simplify(ans)
for marked in [False,True]:
    for m in range(4):
        hv=haar((marker if marked else 1)*phase**m)
        gv=gauss((M if marked else 1)*L**m)
        k=2*m+int(marked)
        check(f'radial identity marked{marked} m{m}',gv==2**k*factorial(k+1)*hv)
receipt={'status':'exact_controls_passed','check_count':len(checks),'checks':checks,
 'dependency_sources':inputs,
 'proof_sha256':hashlib.sha256((BASE/'PROOF_DRAFT.tex').read_bytes()).hexdigest(),
 'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'audit_sha256':hashlib.sha256((BASE/'audit.json').read_bytes()).hexdigest(),
 'universal_proofs_verified_by_samples':False,
 'topological_cover_classification_formally_verified':False,
 'exceptional_group_invariance_recomputed':False,
 'new_numeric_nilpotence_exponent_computed':False,
 'independent_review':False,'formal_verification':False,
 'registry_modified':False}
(BASE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('PASS',len(checks),'exact controls; universal claims remain written unregistered deductions.')
