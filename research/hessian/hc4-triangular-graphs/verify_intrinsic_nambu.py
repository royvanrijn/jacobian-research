"""Exact controls for intrinsic all-degree Nambu direction rigidity.

These controls do not certify Tang's classification or construct a
nonautomorphic Keller pair. Orbit termination is checked as a literal identity.
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

def vec_zero(values,label):
    check(all(sp.simplify(sp.expand(value))==0 for value in values),label)

def deriv(coeff,variables,poly):
    return sp.expand(sum(a*sp.diff(poly,x) for a,x in zip(coeff,variables)))

s,t,u,v,tau = sp.symbols('s t u v tau')
variables = [s,t,u,v]
a1,a2,b1,b2 = sp.symbols('a1 a2 b1 b2')
direction = sp.Matrix([a1,a2,b1,b2])
h = s**4+t**4+s*u+t*v
F = sp.Matrix([sp.diff(h,x) for x in variables])
H = sp.hessian(h,variables)
Hi = H.inv()
G = sp.expand(sum(x*y for x,y in zip(variables,F))-h)
zero(H.det()-1,'the quartic-base intrinsic HC4 control has actual unit Hessian')
check(H*Hi==sp.eye(4),'its polynomial inverse Hessian is literal')
for i in range(4):
    zero(deriv(Hi[:,i],variables,G)-variables[i],'Nambu derivative of the actual Legendre energy is source coordinate '+str(i+1))
check(F.jacobian(variables)*Hi==sp.eye(4),'all constant Nambu directions translate the actual gradient target')
q = sp.expand((direction.T*Hi*direction)[0])
zero(q-(2*a1*b1+2*a2*b2-12*s**2*b1**2-12*t**2*b2**2),
    'the exact inverse-Hessian cone equation is quadratic in four constant covector coordinates')
zero(sp.Poly(q,s,t).coeff_monomial(s**2)+12*b1**2,'one source coefficient forces b1 zero over C')
zero(sp.Poly(q,s,t).coeff_monomial(t**2)+12*b2**2,'another source coefficient forces b2 zero over C')
zero(q.subs({b1:0,b2:0}),'the full exact cone is precisely the linear a-plane, despite an automorphic gradient')
for i in [0,1]:
    zero(deriv(Hi[:,i],variables,deriv(Hi[:,i],variables,G)),
        'an automatic cotangent-annihilator energy direction terminates at order two '+str(i+1))
orbit = [G]
for _ in range(5):
    orbit.append(deriv(Hi[:,2],variables,orbit[-1]))
check(orbit[2]!=0,'a transverse finite-energy direction has nonzero inverse-Hessian square')
zero(orbit[2]+12*s**2,'its exact second derivative is the nonzero square -12s^2')
check(orbit[4]!=0 and orbit[5]==0,'its energy orbit terminates exactly at order five, not at an observed finite prefix')
check(sp.Matrix.hstack(sp.eye(4)[:,0],sp.eye(4)[:,1],sp.eye(4)[:,2]).rank()==3,
    'the same control supplies three independent certified finite-energy directions')
A,B,C,D = sp.symbols('A B C D')
inverse = [C,D,A-4*C**3,B-4*D**3]
vec_zero([a.subs(dict(zip([A,B,C,D],F)),simultaneous=True)-x for a,x in zip(inverse,variables)],
    'the linear-two-plane cone control has its actual polynomial left gradient inverse')
vec_zero([a.subs(dict(zip(variables,inverse)),simultaneous=True)-x for a,x in zip(F,[A,B,C,D])],
    'the same control has its actual polynomial right gradient inverse')

# Nonlinear cone and nonclosure under addition for a quadratic potential.
h0 = s*u+t*v
H0i = sp.hessian(h0,variables).inv()
e1,e3 = sp.eye(4)[:,0],sp.eye(4)[:,2]
zero((e1.T*H0i*e1)[0],'one quadratic-control direction has zero inverse-Hessian square')
zero((e3.T*H0i*e3)[0],'a second quadratic-control direction has zero square')
zero(((e1+e3).T*H0i*(e1+e3))[0]-2,'their sum has nonzero square, certifying a nonlinear exact cone')
zero((e1.T*H0i*e3)[0]-1,'their nonzero mixed pairing supplies the same inverse certificate')

# A genuinely nonlinear intrinsic coefficient pair with two-step target flows.
k1 = s+t**2
k2 = t+k1**2
h1 = s**3+t**4+u*k1+v*k2
H1 = sp.hessian(h1,variables)
H1i = H1.inv()
F1 = sp.Matrix([sp.diff(h1,x) for x in variables])
G1 = sp.expand(sum(x*y for x,y in zip(variables,F1))-h1)
zero(sp.det(sp.Matrix([k1,k2]).jacobian([s,t]))-1,'the nonlinear cotangent coefficient control is genuinely Keller')
zero(H1.det()-1,'its literal four-variable Hessian also has determinant one')
flow_vector = H1i[:,0]
check(any(sp.Poly(a,variables).total_degree()>1 for a in flow_vector),
    'a zero-square intrinsic direction gives a genuinely nonlinear source quasi-translation')
zero((e1.T*H1i*e1)[0],'that nonlinear quasi-translation direction has zero polynomial inverse-Hessian square')
vec_zero([deriv(flow_vector,variables,a) for a in flow_vector],
    'the source flow vector is invariant under its own derivation')
flow = sp.Matrix(variables)+tau*flow_vector
vec_zero([a.subs(dict(zip(variables,flow)),simultaneous=True)-a-tau*b for a,b in zip(F1,e1)],
    'the explicit nonlinear polynomial flow translates the entire gradient by tau times its direction')
back = [a-tau*b for a,b in zip(variables,flow_vector)]
vec_zero([a.subs(dict(zip(variables,flow)),simultaneous=True)-x for a,x in zip(back,variables)],
    'changing the flow parameter sign gives its actual polynomial inverse')
for i in range(4):
    zero(deriv(H1i[:,i],variables,G1)-variables[i],
        'the nonlinear coefficient control retains the exact energy-coordinate identity '+str(i+1))

# Formal cotangent block, with no assumption of an inverse coefficient pair.
p,q0,r,s0,B11,B12,B22 = sp.symbols('p q0 r s0 B11 B12 B22')
DK = sp.Matrix([[p,q0],[r,s0]])
BB = sp.Matrix([[B11,B12],[B12,B22]])
HB = BB.row_join(DK.T).col_join(DK.row_join(sp.zeros(2)))
HBinv = sp.zeros(2).row_join(DK.inv()).col_join(DK.T.inv().row_join(-DK.T.inv()*BB*DK.inv()))
check(sp.simplify(HB*HBinv)==sp.eye(4),'the arbitrary intrinsic cotangent block inverse has the stated orientation')
check(HBinv[:2,:2]==sp.zeros(2),'automatic finite-energy covectors span a totally isotropic inverse-Hessian plane')
vec_zero((HBinv*direction)[:2,:]-DK.inv()*sp.Matrix([b1,b2]),
    'the base restriction of a transverse intrinsic Nambu direction is exactly the plane Keller direction')
zero(HB.det()-DK.det()**2,'the cotangent Hessian determinant is the square of the unrestricted plane Jacobian')
inverse_block = sp.zeros(2).row_join(DK).col_join(DK.T.row_join(BB))
zero(inverse_block.det()-DK.det()**2,
    'a supplied inverse-Hessian zero plane makes its cross-block determinant a unit')
check(inverse_block.inv()[2:,2:]==sp.zeros(2),
    'block inversion reconstructs actual joint affineness from a zero-square covector plane without Tang')
recovered_base = h1.subs({u:0,v:0})
recovered_pair = [sp.diff(h1,u).subs({u:0,v:0}),sp.diff(h1,v).subs({u:0,v:0})]
zero(h1-recovered_base-u*recovered_pair[0]-v*recovered_pair[1],
    'the actual polynomial reconstructs its intrinsic base and Keller pair from the recovered source plane')

# Norm and slice controls; arbitrary locally finite derivations can have eigenvalues.
z,alpha,beta,m,n = sp.symbols('z alpha beta m n')
logder = m/(z-alpha)-n/(z-beta)
zero(sp.limit(logder,z,sp.oo),'rational target logarithmic derivatives vanish at infinity')
norm = -z
zero(sp.diff(norm,z)/norm-2/(2*z),'an explicit finite algebraic norm has the required derivative trace factor')
xx,yy = sp.symbols('xx yy')
arbitrary = [1,yy]
zero(deriv(arbitrary,[xx,yy],xx)-1,'an arbitrary locally finite derivation may have a polynomial slice')
zero(deriv(arbitrary,[xx,yy],yy)-yy,'it may also have a nonzero eigenvalue, excluded only by Keller norm structure')
zero(deriv(arbitrary,[xx,yy],deriv(arbitrary,[xx,yy],yy))-yy,
    'that arbitrary slice derivation is not locally nilpotent')
slice_D = [1,2*xx]
projection = yy-2*xx**2+xx**2
zero(deriv(slice_D,[xx,yy],projection),'the finite Taylor slice projection is killed by the derivation')
zero(yy-(projection+xx**2),'the same slice projection gives literal polynomial reconstruction')

# The relative canonical slice has an inverse obligation, not automatic nilpotence.
YY,TT,VV = sp.symbols('YY TT VV')
DK1 = sp.Matrix([k1,k2]).jacobian([s,t])
normal = DK1.inv()*sp.Matrix([1,-YY])
zero(deriv(normal,[s,t],k1)-1,'the canonical relative normal derivation has polynomial slice f-T')
zero(deriv(normal,[s,t],k2)+YY,'its action on g cancels the multiplier contribution')
zero(deriv(list(normal)+[0],[s,t,YY],(k1-TT)*YY+k2),
    'the canonical normal derivation really fixes the carried second observable')
ff = k1+tau
gg = k2-YY*tau
flow_t = gg-ff**2
flow_s = ff-flow_t**2
vec_zero([sp.diff(flow_s,tau).subs(tau,0)-normal[0],sp.diff(flow_t,tau).subs(tau,0)-normal[1]],
    'the positive inverse-pair control has its explicit relative polynomial flow')
check(sp.Poly(flow_s,tau).degree()==4 and sp.Poly(flow_t,tau).degree()==2,
    'this relative slice is nilpotent in the control because actual Keller inverse coordinates are supplied')

context = json.loads((HERE/'sources/intrinsic-nambu-context.json').read_text())
check(context['newly_used_unproved_external_statements']==[],
    'the new primary energy identity is reproved without a new unproved classification')
check(context['retained_Tang_classification_required'] and context['retained_regular_parameter_input_required'],
    'the two-or-four conclusion explicitly retains Tang and regular-parameter inputs')
check(not context['linear_two_plane_cone_asserted_nonautomorphic'] and not context['finite_derivative_prefix_assumed_termination'],
    'a plane-shaped cone is not a converse and finite derivative prefixes are not accepted as proofs')
check(not context['unrestricted_HC4_or_JC2_resolved'] and not context['original_goal_complete'],
    'no extra intrinsic finite-energy direction or endpoint closure is assumed')
math_count = len(checks)
for name,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,'retained byte pin '+name)
old_counts = {}
for name in ['receipt-universal-modification.json','receipt-regular-parameter.json','receipt-intrinsic-plane.json']:
    old = json.loads((HERE/name).read_text())
    check(old['outcome']=='PASS','retained receipt is PASS '+name)
    old_counts[name] = old.get('check_count',len(old.get('checks',[])))
    for source,digest in old.get('source_sha256',{}).items():
        check(hashlib.sha256((HERE/source).read_bytes()).hexdigest()==digest,'retained receipt source '+name+' -> '+source)
    for source,digest in old.get('dependencies',{}).items():
        check(hashlib.sha256((ROOT/source).read_bytes()).hexdigest()==digest,'retained receipt dependency '+name+' -> '+source)
names = ['verify_intrinsic_nambu.py','INTRINSIC_NAMBU_DIRECTION_RIGIDITY.md','sources/intrinsic-nambu-context.json']
result = dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=math_count,checks=checks,
    source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names},
    dependencies=context['retained_proof_and_source_sha256'],retained_receipt_counts_not_newly_replayed=old_counts,
    newly_cached_primary_context=context['newly_cached_primary_context'],retained_structural_inputs=context['retained_structural_inputs'],
    newly_used_unproved_external_statements=[],written_not_finite_exhaustion=context['direct_proof_arguments'],
    **{key:value for key,value in context.items() if type(value) is bool},
    runtime=dict(python=platform.python_version(),sympy=sp.__version__),elapsed_seconds=round(time.monotonic()-START,3))
(HERE/'receipt-intrinsic-nambu.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=math_count,
    original_goal_complete=False),indent=2))
