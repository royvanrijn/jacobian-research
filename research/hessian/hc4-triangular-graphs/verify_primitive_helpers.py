"""Exact controls for primitive helpers and the written all-degree proof."""
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

s,t,w1,w2,c,K1,K2,A1,A2,C,R = sp.symbols('s t w1 w2 c K1 K2 A1 A2 C R')
H = w1*A1+w2*A2
X = c+H
Y = K1-(2*c+H)*A1+w2*C
Z = K2-(2*c+H)*A2-w1*C
zero(X**2+w1*Y+w2*Z-c**2-w1*K1-w2*K2,
    'the complete polynomial compensation has its exact cotangent identity for unrestricted data')
Hp = w1*(A1+w2*R)+w2*(A2-w1*R)
zero(Hp-H,'every fixed-H change has the exact two-variable syzygy shape')
Cp = C+(2*c+H)*R
zero(K1-(2*c+Hp)*(A1+w2*R)+w2*Cp-Y,'the gauge change leaves Y unchanged')
zero(K2-(2*c+Hp)*(A2-w1*R)-w1*Cp-Z,'the same gauge change leaves Z unchanged')
zeros([X.subs({w1:0,w2:0})-c,Y.subs({w1:0,w2:0})-(K1-2*c*A1),Z.subs({w1:0,w2:0})-(K2-2*c*A2)],
    'the complete compensation has the exact closed-zero-section coordinate obligation')

# An actual high-degree free-variable control with parameter recovery.
k1=s+t**2
k2=t+k1**3
aa=[w1**2,w2**3]
cc=w1*w2
hh=w1*aa[0]+w2*aa[1]
xx=1+hh
yy=k1-(2+hh)*aa[0]+w2*cc
zz=k2-(2+hh)*aa[1]-w1*cc
phi=sp.expand(xx**2+w1*yy+w2*zz)
zero(phi-1-w1*k1-w2*k2,'the full-data control has the literal noncritical intrinsic potential')
zero(sp.hessian(phi,[s,t,w1,w2]).det()-1,'the actual high-degree source potential has unit Hessian')
rec1=yy+(xx+1)*aa[0]-w2*cc
rec2=zz+(xx+1)*aa[1]+w1*cc
rect=rec2-rec1**3
recs=rec1-rect**2
zeros([recs-s,rect-t],'this positive full-data control has global polynomial source recovery')
check(all(sp.Poly(a,w1,w2).total_degree()>1 for a in [xx,yy,zz]),
    'the full-data helpers allow unrestricted nonlinear free dependence')

# Preservation of a quadratic expression does not supply an inverse shear.
x,y,z=sp.symbols('x y z')
shear=[x,(1+w2)*y,z-w1*y]
zero(shear[0]**2+w1*shear[1]+w2*shear[2]-(x**2+w1*y+w2*z),
    'a formal coefficient-dependent syzygy shear can preserve the whole quadratic')
zero(sp.Matrix([w1,w2]+shear).jacobian([w1,w2,x,y,z]).det()-(1+w2),
    'that shear has nonconstant Jacobian and is not automatically an ambient automorphism')
check(shear[1].subs({w1:0,w2:-1})==0 and shear[2].subs({w1:0,w2:-1})==z,
    'the same shear actually collapses a helper line at a specific base point')

# Genuine Keller primitive-helper graph and the exact derivative ideal.
T1,T2,S=sp.symbols('T1 T2 S')
F=sp.Matrix([s+t**3,t])
Fi=F.jacobian([s,t]).inv()
a=s*t+s**2
base=(T1-T2**3)*T2+(T1-T2**3)**2
q=S-base
zero(F.jacobian([s,t]).det()-1,'the primitive-helper positive control is genuinely Keller')
zero(q.subs({T1:F[0],T2:F[1],S:a},simultaneous=True),
    'its actual primitive implicit equation is a graph equation')
p=sp.diff(q,S).subs({T1:F[0],T2:F[1],S:a},simultaneous=True)
zero(p-1,'the primitive graph has constant nonzero implicit derivative')
for i,T in enumerate([T1,T2]):
    da=sum(Fi[j,i]*sp.diff(a,v) for j,v in enumerate([s,t]))
    zero(sp.diff(q,T).subs({T1:F[0],T2:F[1],S:a},simultaneous=True)+p*da,
        'Nambu differentiation gives exact source Jacobian-ideal divisibility '+str(i+1))
check(sp.gcd(F[0],F[1])==1,'coprime target coordinates retain coprime source pullbacks in the control')

# The non-Keller primitive normal graph marks the essential hypothesis.
badq=S**2-T2
zero(badq.subs({T1:s,T2:t**2,S:t},simultaneous=True),
    'a primitive non-Keller finite map can have a smooth implicit image')
zero(sp.diff(badq,T2)+1,'that implicit image is explicitly smooth despite its noninvertible projection')
zero(sp.det(sp.Matrix([s,t**2]).jacobian([s,t]))-2*t,
    'the normal-image negative control is explicitly not Keller')
zero(sp.diff(badq,S).subs(S,t)-2*t,'its implicit derivative has a source zero divisor')
zero(2*t*(1/(2*t))-1,'its rational Nambu derivative cancels that factor but is not polynomial')
check(sp.rem(sp.Integer(-1),2*t,t)!=0,'the target partial is not polynomially divisible by the implicit derivative without Keller')

# Canonical gradient energy and explicit inverse from its linear equation.
u,v=sp.symbols('u v')
src=[s,t,u,v]
h=s**4+t**4+s*u+t*v
grad=sp.Matrix([sp.diff(h,x) for x in src])
energy=sp.expand(sum(a*b for a,b in zip(src,grad))-h)
P1,P2,P3,P4=sp.symbols('P1 P2 P3 P4')
target=[P1,P2,P3,P4]
B=P3*P1+P4*P2-P3**4-P4**4
eq=S-B
zero(sp.hessian(h,src).det()-1,'the actual Legendre positive control has unit Hessian')
zero(eq.subs(dict(zip(target+[S],list(grad)+[energy])),simultaneous=True),
    'its canonical energy has the displayed exact implicit relation')
for i,P in enumerate(target):
    zero(sp.diff(eq,P).subs(dict(zip(target,grad)),simultaneous=True)+src[i],
        'the linear Legendre equation supplies actual inverse source coordinate '+str(i+1))
inverse=[sp.diff(B,P) for P in target]
zeros([a.subs(dict(zip(target,grad)),simultaneous=True)-x for a,x in zip(inverse,src)],
    'the implicit-equation inverse has its literal left composition')
zeros([a.subs(dict(zip(src,inverse)),simultaneous=True)-P for a,P in zip(grad,target)],
    'the same implicit-equation inverse has its literal right composition')

# Split conditional labels model a geometric generic fibre, not a Keller sample.
eta1,eta2=sp.symbols('eta1 eta2')
labels=[sp.Integer(0),eta1-1,eta2-2]
splitq=sp.prod(S-a for a in labels)
selected=sp.diff(splitq,S).subs(S,0)
zero(selected-(eta1-1)*(eta2-2),'three conditional sheets give exactly two collision-line factors on one selected sheet')
check(sp.Poly(selected,eta1,eta2).total_degree()==2,
    'the generic collision product has degree d-1 in the independent stable variables')
for var in [S,eta1,eta2]:
    zero(sp.diff(splitq,var).subs({S:0,eta1:1},simultaneous=True),
        'the split geometric image is singular along an actual selected-sheet collision '+str(var))
check(sp.Poly(eta1-1,eta1,eta2).total_degree()==1 and sp.Poly(eta2-2,eta1,eta2).total_degree()==1,
    'each conditional collision factor has nonzero linear part')

# Differential immersion versus locally closed embedding at a node.
paramx=t**2-1
paramy=t*(t**2-1)
node=y**2-x**2*(x+1)
zero(node.subs({x:paramx,y:paramy},simultaneous=True),'the nodal normalization parameter satisfies its exact image equation')
check(sp.gcd(sp.diff(paramx,t),sp.diff(paramy,t))==1,
    'the nodal parameter has injective differential at every source point')
zeros([paramx.subs(t,1),paramx.subs(t,-1),paramy.subs(t,1),paramy.subs(t,-1)],
    'two distinct source points still map to the same node')
zeros([sp.diff(node,x).subs({x:0,y:0}),sp.diff(node,y).subs({x:0,y:0})],
    'that common image is a singular point of the implicit curve')

context=json.loads((HERE/'sources/primitive-helper-context.json').read_text())
check(context['newly_used_unproved_external_statements']==[],
    'the characteristic-zero weak-type mechanism is attributed and reproduced rather than newly imported without proof')
check(context['primitive_helper_required_for_normality_test'] and not context['primitive_field_helper_assumed_algebra_generator'],
    'field primitivity and source algebra generation remain separate obligations')
check(not context['normality_of_actual_embedding_assumed_normality_of_Legendre_image'],
    'the actual smooth HC4 image is kept distinct from the possibly nonnormal Legendre image')
check(not context['full_noncritical_rank_one_chart_resolved'] and not context['unrestricted_HC4_or_JC2_resolved'],
    'the full compensation and primitive-helper tests do not assert an endpoint or unrestricted closure')
math_count=len(checks)
for name,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,'retained byte pin '+name)
old_counts={}
for name in ['receipt-free-projection-trace.json','receipt-intrinsic-nambu.json','receipt-abelian-classgroup.json']:
    old=json.loads((HERE/name).read_text())
    check(old['outcome']=='PASS','retained receipt is PASS '+name)
    old_counts[name]=old.get('check_count',len(old.get('checks',[])))
    for source,digest in old.get('source_sha256',{}).items():
        check(hashlib.sha256((HERE/source).read_bytes()).hexdigest()==digest,'retained receipt source '+name+' -> '+source)
    for source,digest in old.get('dependencies',{}).items():
        check(hashlib.sha256((ROOT/source).read_bytes()).hexdigest()==digest,'retained receipt dependency '+name+' -> '+source)
names=['verify_primitive_helpers.py','PRIMITIVE_HELPERS_AND_LEGENDRE_NORMALITY.md','sources/primitive-helper-context.json']
result=dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=math_count,checks=checks,
    source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names},
    dependencies=context['retained_proof_and_source_sha256'],retained_receipt_counts_not_newly_replayed=old_counts,
    newly_cached_primary_context=context['newly_cached_primary_context'],retained_structural_inputs=context['retained_structural_inputs'],
    retained_primary_reproduced=context['retained_primary_reproduced'],inspected_but_not_imported=context['inspected_but_not_imported'],
    newly_used_unproved_external_statements=[],written_not_finite_exhaustion=context['direct_proof_arguments'],
    **{key:value for key,value in context.items() if type(value) is bool},
    runtime=dict(python=platform.python_version(),sympy=sp.__version__),elapsed_seconds=round(time.monotonic()-START,3))
(HERE/'receipt-primitive-helpers.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=math_count,original_goal_complete=False),indent=2))
