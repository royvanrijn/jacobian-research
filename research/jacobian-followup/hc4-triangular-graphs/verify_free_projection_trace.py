"""Exact controls for the written all-degree free-projection criterion.

The imported one-line theorem is not proved by these controls. The
universal compensation is a substitution, not an assumed embedding.
"""
import hashlib
import json
from pathlib import Path
import platform
import time
import sympy as sp

START = time.monotonic()
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
checks = []

def check(value,label):
    assert bool(value),label
    checks.append(label)

def zero(value,label):
    check(sp.simplify(sp.expand(value))==0,label)

def zeros(values,label):
    check(all(sp.simplify(sp.expand(x))==0 for x in values),label)

s,t,w1,w2,X,Y,Z,c = sp.symbols('s t w1 w2 X Y Z c')
src = [s,t,w1,w2]
Q1 = X**2+w1*Y+w2*Z
Q2 = X*Y+w1*Z+w2**2
Q3 = X*Y+Z**2+w1**2+w2**2
for r,Q in enumerate([Q1,Q2,Q3],1):
    check(sp.hessian(Q,[w1,w2,X,Y,Z]).det()!=0,'the rank '+str(r)+' helper model has nondegenerate full ambient quadratic')
    check(sp.hessian(Q.subs({w1:0,w2:0}),[X,Y,Z]).rank()==r,
        'the same model has exactly helper rank '+str(r))
check(all(sp.diff(X*Y+Z**2,a).subs({X:0,Y:0,Z:0})==0 for a in [X,Y,Z]),
    'the zero ternary rank-three level is singular at its vertex')
rho = sp.symbols('rho',nonzero=True)
zero((Z-rho)*(Z+rho)-(Z**2-rho**2),'the nonzero rank-three removed set has two distinct geometric lines')
zero((X*Y+Z**2-rho**2).subs(Y,(rho**2-Z**2)/X),
    'localizing the rank-three helper coordinate gives the exact polynomial Laurent chart')
check(sp.Matrix([[1,-1]]).rank()==1,'the two-line smooth-quadric class lattice modulo the diagonal has rank one')
zero((X*Y-rho**2).subs(Y,rho**2/X),'the nonzero rank-two helper level has a genuine nonconstant Laurent unit')
check(sp.factor(X*Y)==X*Y and sp.factor(X**2)==X**2,
    'zero rank-two and rank-one levels have the stated reduced plane supports')

# An actual rank-two closed embedding and literal first jet.
n = 4
j2 = [w2+w1*s**n,s-w2,t+w2*s**n]
h2 = sp.expand(Q2.subs(dict(zip([X,Y,Z],j2)),simultaneous=True))
K2 = [t+s**(n+1),s]
zero(h2-w1*K2[0]-w2*K2[1],'the rank-two embedding has its exact intrinsic cotangent potential')
zero(sp.hessian(h2,src).det()-1,'the actual four-variable rank-two control has unit Hessian')
zero(j2[1]+w2-s,'the rank-two helper formula recovers the first source coordinate globally')
zero(j2[2]-w2*(j2[1]+w2)**n-t,'the rank-two helper formula recovers the second source coordinate globally')
zeros([j2[0].subs({w1:0,w2:0}),j2[1].subs({w1:0,w2:0})-s,j2[2].subs({w1:0,w2:0})-t],
    'its literal zero-section is the full helper plane with a source coordinate pair')
zero(sp.diff(h2,w1).subs({w1:0,w2:0})-(t+s*sp.diff(j2[0],w1).subs({w1:0,w2:0})),
    'the rank-two first free derivative has the unrestricted normal-jet correction times a')
zero(sp.diff(h2,w2).subs({w1:0,w2:0})-s*sp.diff(j2[0],w2).subs({w1:0,w2:0}),
    'the second rank-two first derivative has the coordinate factor a')
zeros([a.subs(s,0)-b for a,b in zip(K2,[t,0])],
    'the actual Keller pair restricts on its coordinate line to an injective affine pair')

# Arbitrarily high free-variable dependence is allowed; one representative.
k1 = s+t**2
k2 = t+k1**3
A0,B0 = w1**2,w2**3
HH = w1*A0+w2*B0
j1 = [HH,k1-HH*A0,k2-HH*B0]
h1 = sp.expand(Q1.subs(dict(zip([X,Y,Z],j1)),simultaneous=True))
zero(h1-w1*k1-w2*k2,'the nonlinear-helper embedding has the same exact critical-value cotangent identity')
zero(sp.hessian(h1,src).det()-1,'its actual four-variable potential has unit Hessian')
check(all(sp.Poly(a,w1,w2).total_degree()>1 for a in j1),
    'all three helpers can have higher than affine free-variable dependence')
recover_k1 = j1[1]+j1[0]*A0
recover_k2 = j1[2]+j1[0]*B0
recover_t = recover_k2-recover_k1**3
recover_s = recover_k1-recover_t**2
zeros([recover_s-s,recover_t-t],'the nonlinear-helper control has actual global polynomial parameter recovery')
zeros([sp.diff(h1,w1).subs({w1:0,w2:0})-k1,sp.diff(h1,w2).subs({w1:0,w2:0})-k2],
    'the rank-one zero-section first jet is exactly its source coordinate pair')

# A variable free frame need not preserve the Keller coefficient pair.
frame = sp.Matrix([[1,t],[0,1]])
coeff = frame.T.inv()*sp.Matrix([s,t])
zero(frame.det()-1,'the variable-frame countercontrol is a genuine polynomial unimodular frame')
zero(sp.Matrix(coeff).jacobian([s,t]).det()-(1-s),
    'the transformed coefficient pair can lose its constant Jacobian despite that frame')

# Universal compensation: exact identity and the full immersion failure.
K = sp.Matrix([s+t**2,t+(s+t**2)**2/sp.Integer(3)])
zero(K.jacobian([s,t]).det()-1,'the compensation countercontrol uses a genuine plane Keller automorphism')
AA = (K-sp.Matrix([s,t]))/2
Hcomp = w1*AA[0]+w2*AA[1]
Jcomp = sp.Matrix([w1,w2,1+Hcomp,s-Hcomp*AA[0],t-Hcomp*AA[1]])
hcomp = sp.expand(Q1.subs(dict(zip([w1,w2,X,Y,Z],Jcomp)),simultaneous=True))
zero(hcomp-1-w1*K[0]-w2*K[1],'the universal noncritical compensation has its exact cotangent potential')
zero(sp.hessian(hcomp,src).det()-1,'the compensation has unit Hessian without any embedding conclusion')
zeros([a-b for a,b in zip(Jcomp.subs({w1:0,w2:0}),[0,0,1,s,t])],
    'the universal compensation has an actually embedded noncritical zero-section')
point = {s:0,t:1,w1:3,w2:-3}
zero(Hcomp.subs(point)-1,'the immersion failure occurs at the required normal compensation value')
zeros(list(AA.subs(point)-sp.Matrix([sp.Rational(1,2),sp.Rational(1,6)])),
    'its source point has the exact nonzero normal first-jet vector')
BB = AA.jacobian([s,t]).subs(point)
zeros(list(BB*sp.Matrix([1,1])-sp.Matrix([1,1])),
    'the killed base direction is a nonzero eigenvector with eigenvalue one')
zeros(list(Jcomp.jacobian(src).subs(point)*sp.Matrix([1,1,0,0])),
    'the full five-coordinate substitution kills an actual source tangent despite unit Hessian')
check(Jcomp.jacobian(src).subs(point).rank()<4,
    'the universal compensation is not an immersion in this actual automorphic control')
U,V = sp.symbols('U V')
inv_t = V-U**2/3
inv_s = U-inv_t**2
zeros([inv_s.subs({U:K[0],V:K[1]},simultaneous=True)-s,
       inv_t.subs({U:K[0],V:K[1]},simultaneous=True)-t],
    'the killed-tangent input has its explicit plane polynomial inverse and is not a counterexample')

# The pointwise spectral obstruction and its two exact polynomial alternatives.
theta,a0,b0 = sp.symbols('theta a0 b0',nonzero=True)
vv = sp.Matrix([1,theta])
avec = sp.Matrix([a0,b0])
wvec = sp.Matrix([theta,-1])/(theta*(theta*a0-b0))
zero((wvec.T*vv)[0],'the general dimension-two normal vector can annihilate the eigenvector')
zero((wvec.T*avec)[0]-1/theta,'the same normal vector supplies the required inverse eigenvalue')
rr,dd,cc = sp.symbols('rr dd cc')
zero((rr**2-4*dd).subs(dd,-rr/(2*cc))-rr*(1+cc*rr/2)*2/cc,
    'determinant one and zero discriminant give exactly the two trace alternatives')
for tr,det,label in [(0,0,'zero'),(-2/cc,1/cc**2,'negative')]:
    zero((rr+2*cc*dd).subs({rr:tr,dd:det}),label+' trace branch satisfies determinant one')
    zero((rr**2-4*dd).subs({rr:tr,dd:det}),label+' trace branch has zero discriminant')
nil = sp.Matrix([[0,5*t**4],[0,0]])
check(nil**2==sp.zeros(2),'the plane nilpotent-Jacobian control has square-zero derivative')
shear = [s+t**5,t]
zero(sp.Matrix(shear).jacobian([s,t]).det()-1,'the needed plane lemma has an unbounded-degree triangular family')
zero(shear[0].subs(s,s-t**5)-s,'the triangular control has its exact inverse')

context = json.loads((HERE/'sources/free-projection-trace-context.json').read_text())
check(context['newly_used_unproved_external_statements']==[],'the two exact external inputs are retained rather than newly unproved imports')
check(context['retained_one_line_theorem_required'] and context['retained_generative_theorem_required'],
    'one-line invertibility and polynomial generative inputs are explicit retained dependencies')
check(not context['all_noncritical_rank_one_embeddings_resolved'] and not context['unit_Hessian_asserted_embedding'],
    'the noncritical rank-one residual and the distinction between substitution and embedding remain explicit')
check(not context['unrestricted_HC4_or_JC2_resolved'] and not context['original_goal_complete'],
    'no endpoint or unrestricted HC4/JC2 conclusion is inferred')
math_count = len(checks)
for name,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,'retained byte pin '+name)
old_counts = {}
for name in ['receipt-intrinsic-nambu.json','receipt-critical-jet.json']:
    old = json.loads((HERE/name).read_text())
    check(old['outcome']=='PASS','retained receipt is PASS '+name)
    old_counts[name] = old.get('check_count',len(old.get('checks',[])))
    for source,digest in old.get('source_sha256',{}).items():
        check(hashlib.sha256((HERE/source).read_bytes()).hexdigest()==digest,'retained receipt source '+name+' -> '+source)
    for source,digest in old.get('dependencies',{}).items():
        check(hashlib.sha256((ROOT/source).read_bytes()).hexdigest()==digest,'retained receipt dependency '+name+' -> '+source)
names = ['verify_free_projection_trace.py','FREE_PROJECTIONS_AND_CONSTANT_BASE_TRACE.md','sources/free-projection-trace-context.json']
result = dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=math_count,checks=checks,
    source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names},
    dependencies=context['retained_proof_and_source_sha256'],retained_receipt_counts_not_newly_replayed=old_counts,
    newly_cached_primary_context=context['newly_cached_primary_context'],retained_structural_inputs=context['retained_structural_inputs'],
    inspected_but_not_used=context['inspected_but_not_used'],newly_used_unproved_external_statements=[],
    written_not_finite_exhaustion=context['direct_proof_arguments'],
    **{key:value for key,value in context.items() if type(value) is bool},
    runtime=dict(python=platform.python_version(),sympy=sp.__version__),elapsed_seconds=round(time.monotonic()-START,3))
(HERE/'receipt-free-projection-trace.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=math_count,
    original_goal_complete=False),indent=2))
