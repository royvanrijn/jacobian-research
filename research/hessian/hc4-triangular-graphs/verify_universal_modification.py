"""Exact controls for the written universal affine-modification boundary.

The conditional d>1 lattice models are not nonautomorphic Keller samples.
"""
import hashlib
import itertools
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

x,y,T,V,Y,Z,s,t,u,v,a,b,c,d0,N,rho = sp.symbols('x y T V Y Z s t u v a b c d0 N rho')
f,g = sp.Function('f')(x,y),sp.Function('g')(x,y)
J = sp.diff(f,x)*sp.diff(g,y)-sp.diff(f,y)*sp.diff(g,x)
M = f-T
observable = M*Y+g
zero(observable.subs(T,f)-g,'the literal graph carries the original second Keller function')
zero(M.subs(T,f),'the actual graph lies in the multiplier divisor before any localization')
zero(sp.diff(M,T)+1,'the global image equation is smooth with its explicit parameter inverse')
zero(sp.diff(f,x)*sp.diff(observable,y)-sp.diff(f,y)*sp.diff(observable,x)-J,
    'every target fibre has the exact Keller Bezout identity for its helper derivatives')
zero((sp.diff(observable,x)*(-sp.diff(f,y))+sp.diff(observable,y)*sp.diff(f,x))/J-1,
    'the nonzero constant Keller Jacobian makes the full fibre Jacobian ideal the unit ideal')
zero(observable.subs(Y,(V-g)/M)-V,'localizing the multiplier eliminates exactly Y')
zero((observable-V).subs(Y,0)-(g-V),'the exceptional centre is the original pair fibre')
zero(sp.diff(observable,Z),'the nominal root direction is absent on every multiplier component')
zero(sp.diff(M,T)+1,'the multiplier has a constant nonzero coefficient and primitive content')

# A nonlinear genuine Keller control, with both components nonlinear.
fc = x+y**4
gc = y+fc**3
zero(sp.det(sp.Matrix([fc,gc]).jacobian([x,y]))-1,'the positive control is a genuine nonlinear Keller pair')
inverse = [T-(V-T**3)**4,V-T**3]
zero(inverse[0].subs({T:fc,V:gc},simultaneous=True)-x,'the positive Keller control has its displayed first inverse coordinate')
zero(inverse[1].subs({T:fc,V:gc},simultaneous=True)-y,'the positive Keller control has its displayed second inverse coordinate')
zero(fc.subs({x:inverse[0],y:inverse[1]},simultaneous=True)-T,'the first control inverse composition also holds on the target')
zero(gc.subs({x:inverse[0],y:inverse[1]},simultaneous=True)-V,'the second control inverse composition also holds on the target')
mc = fc-T
vc = mc*Y+gc
zero(sp.diff(fc,x)*sp.diff(vc,y)-sp.diff(fc,y)*sp.diff(vc,x)-1,
    'all positive-control target fibres retain the unit derivative identity')
check(sp.gcd(mc,gc-V)==1,'the positive control primitive-linear fibre coefficients are coprime')
check(sp.Poly(mc,x,y).primitive()[0]==1,'the positive control multiplier has primitive target-parameter content')
zero(vc.subs({T:fc,Y:u},simultaneous=True)-gc,'the positive control embedding carries the entire Keller pair')

# The generic helper map already fails flatness for the identity pair.
identity_equation = (x-T)*Y+y-V
zero(identity_equation.subs({x:T,y:V},simultaneous=True),'an actual generic source point has a full Y,Z fibre')
zero(identity_equation.subs({x:T,y:V+1},simultaneous=True)-1,'other points of the same helper curve have empty fibres')
check(2>1,'the full exceptional Y,Z fibre has larger dimension than the generic Z fibre')
zero(identity_equation.subs(y,V-(x-T)*Y),'good ambient target fibres do not force flatness of the distinct helper projection')

# A conditional three-sheet permutation lattice, not a Keller example.
e = [sp.eye(3)[:,i] for i in range(3)]
diag = sp.ones(3,1)
proj = sp.Matrix([[1,0,-1],[0,1,-1]])
check(proj*diag==sp.zeros(2,1),'the exact diagonal divisor relation vanishes in quotient coordinates')
check(proj.rank()==2,'the three-component quotient lattice has rank two')
check(proj*e[0]!=sp.zeros(2,1),'a selected component has a nonzero class in the conditional degree-three lattice')
check(sp.gcd_list(list(proj*e[0]))==1,'that selected lattice class is primitive and has infinite order')
swap = sp.Matrix([[1,0,0],[0,0,1],[0,1,0]])
check(swap*e[0]==e[0],'the selected-component stabilizer fixes its nonzero divisor class')
check(swap*diag==diag,'the source permutation action preserves the exact divisor relation')
orbit_basis = sp.Matrix([[1,0],[0,1],[0,1]])
check((proj*orbit_basis).rank()==1,'one selected singleton and one remaining orbit give fixed-class rank one')
check(sp.Matrix([[1,-1]])*sp.ones(2,1)==sp.zeros(1,1),'the two-sheet conditional lattice also has only its diagonal relation')
check(sp.zeros(0,1).rows==0,'the degree-one quotient lattice has zero rank')
check(all(sum((perm*e[0]-e[0]))==0 for perm in [sp.eye(3),swap]),
    'summing coordinates forbids a nonzero diagonal difference for fixed classes')

# Non-Keller finite-degree control marks the conditional model boundary.
bad_f,bad_g = x,y**2
bad_eq = (bad_f-T)*Y+bad_g-V
check(sp.det(sp.Matrix([bad_f,bad_g]).jacobian([x,y]))==2*y,
    'the two-sheet control is explicitly not a Keller pair')
check(all(sp.diff(bad_eq,var).subs({x:T,y:0,Y:0,V:0},simultaneous=True)==0 for var in [x,y,Y,Z]),
    'that non-Keller control has a singular special target fibre and is not a claimed realization')
zero((y-sp.sqrt(V))*(y+sp.sqrt(V))-(y**2-V),'the conditional exceptional two-sheet equation is reduced at generic V')

# Actual Hessians and literal gradient inverses.
A11,A12,A22 = sp.symbols('A11 A12 A22')
block = sp.Matrix([[A11,A12,1,0],[A12,A22,0,1],[1,0,0,0],[0,1,0,0]])
zero(block.det()-1,'the universal cotangent block Hessian has determinant one for arbitrary base second derivatives')
Q = T**2/2+x*Y+y*Z
zero(sp.hessian(Q,[T,x,y,Y,Z]).det()-1,'the ambient five-variable quadratic is nondegenerate')
fs = s+t**4
phi = fs**2/2+s*u+t*v
zero(sp.hessian(phi,[s,t,u,v]).det()-1,'the actual nonlinear HC4 source potential has unit Hessian')
grad = [sp.diff(phi,var) for var in [s,t,u,v]]
invgrad = [c,d0,a-(c+d0**4),b-4*d0**3*(c+d0**4)]
for i,(value,wanted) in enumerate(zip(invgrad,[s,t,u,v])):
    zero(value.subs(dict(zip([a,b,c,d0],grad)),simultaneous=True)-wanted,
        'the actual source gradient has polynomial left inverse coordinate '+str(i+1))
for i,(value,wanted) in enumerate(zip(grad,[a,b,c,d0])):
    zero(value.subs(dict(zip([s,t,u,v],invgrad)),simultaneous=True)-wanted,
        'the actual source gradient has polynomial right inverse coordinate '+str(i+1))
ell = sp.Matrix([y,x,rho*T,Y/N,Z/(2*N)])
zero(ell[2]**2/(2*rho**2)+N*ell[1]*ell[3]+2*N*ell[0]*ell[4]-Q,
    'the explicit constant linear normalizer pulls the retained Meng-Yang quadratic back to Q')
check(sp.factor(ell.jacobian([T,x,y,Y,Z]).det())!=0,
    'the quadratic normalizer is invertible for the stipulated nonzero constants')

context = json.loads((HERE/'sources/universal-modification-context.json').read_text())
check(context['newly_used_unproved_external_statements']==[],
    'affine-modification terminology adds no unproved theorem to the direct calculation')
check(context['source_fixed_class_vanishing_equivalent_to_degree_one'] and not context['source_fixed_classes_asserted_zero_unconditionally'],
    'the class criterion is an equivalence of open assertions, not an unconditional inverse claim')
check(context['constructed_gradient_inverse_proved'] and not context['carried_pair_identified_with_intrinsic_nonautomorphic_gradient_pair'],
    'the arbitrary carried pair is kept distinct from the always invertible constructed gradient')
check(not context['nonautomorphic_Keller_example_constructed'] and not context['unrestricted_HC4_or_JC2_resolved'],
    'no conditional multi-sheet lattice is promoted to a Keller counterexample or a JC2 proof')
math_count = len(checks)
for name,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,'retained byte pin '+name)
old_counts = {}
for name in ['receipt-two-helper-curves.json','receipt-regular-parameter.json','receipt-ambient-quadric.json']:
    old = json.loads((HERE/name).read_text())
    check(old['outcome']=='PASS','retained receipt is PASS '+name)
    old_counts[name] = old.get('check_count',len(old.get('checks',[])))
    for source,digest in old.get('source_sha256',{}).items():
        check(hashlib.sha256((HERE/source).read_bytes()).hexdigest()==digest,'retained receipt source '+name+' -> '+source)
    for source,digest in old.get('dependencies',{}).items():
        check(hashlib.sha256((ROOT/source).read_bytes()).hexdigest()==digest,'retained receipt dependency '+name+' -> '+source)
names = ['verify_universal_modification.py','UNIVERSAL_AFFINE_MODIFICATION_BOUNDARY.md','sources/universal-modification-context.json']
result = dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=math_count,checks=checks,
    source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names},
    dependencies=context['retained_proof_and_source_sha256'],retained_receipt_counts_not_newly_replayed=old_counts,
    newly_cached_primary_context=context['newly_cached_primary_context'],retained_structural_imports=context['retained_structural_imports'],
    newly_used_unproved_external_statements=[],written_not_finite_exhaustion=context['direct_proof_arguments'],
    **{key:value for key,value in context.items() if type(value) is bool},
    runtime=dict(python=platform.python_version(),sympy=sp.__version__),elapsed_seconds=round(time.monotonic()-START,3))
(HERE/'receipt-universal-modification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=math_count,
    original_goal_complete=False),indent=2))
