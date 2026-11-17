"""Exact controls for the all-degree critical-determinant proof JC0/1."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as sp
START=time.monotonic()
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
checks=[]
def check(value,label):
    assert bool(value),label
    checks.append(label)
def zero(value,label):check(sp.cancel(sp.expand(value))==0,label)
def zeros(values,label):check(all(sp.cancel(sp.expand(v))==0 for v in values),label)
x,y,z,w,t,r,s,b,mu=sp.symbols('x y z w t r s b mu')
variables=[x,y,z,w]
qA=x*x+y*y+z*z+w*w
qB=x*x+2*y*y
outer=t**3/3+t*t/2+2*t
p=sp.diff(outer,t)
h=qA+outer.subs(t,qB)
hessian=sp.hessian(h,variables)
gradient=[sp.diff(h,v) for v in variables]
lam=[1,2,0,0]
for i in [0,1]:
    axis={variables[j]:0 for j in range(4) if j!=i}
    condition=p.subs(t,r)+sp.Rational(1,lam[i])
    axisdet=sp.expand(hessian.subs(axis).det())
    reduced=sp.rem(axisdet,variables[i]**2-r/lam[i],variables[i])
    Ci=sp.prod(1-sp.Rational(ll,lam[i]) for j,ll in enumerate(lam) if j!=i)
    expected=16*2*lam[i]*r*sp.diff(p,t).subs(t,r)*Ci
    zero(sp.rem(sp.expand(reduced-expected),condition,r),
         'the exact critical Hessian factor has its sign and dimension factor axis='+str(i))
    reduced_gradient=[]
    for component in gradient:
        reduced=sp.rem(sp.expand(component.subs(axis)),variables[i]**2-r/lam[i],variables[i])
        reduced_gradient.append(sp.rem(reduced,condition,r))
    zeros(reduced_gradient,'the normalized scalar roots give full critical points axis='+str(i))

# Polynomial coefficient rigidity at arbitrary root degree is proved in the note.
for degree in [1,2,4]:
    coeff=sp.symbols('a0:'+str(degree+1))
    pp=sum(coeff[j]*t**j for j in range(degree+1))
    expected=sum((j-degree)*coeff[j]*t**j for j in range(1,degree))-degree*(coeff[0]+1/mu)
    zero(t*sp.diff(pp,t)-degree*(pp+1/mu)-expected,
         'the all-root equal-value condition removes every intermediate coefficient degree='+str(degree))
    check(sp.degree(t*sp.diff(pp,t),t)==degree,
          'the root comparison has exactly the required outer derivative degree='+str(degree))

for count in [1,2,4]:
    lambdas=sp.symbols('l0:'+str(count),nonzero=True)
    P=sp.prod(1+ll*s for ll in lambdas)
    R=2*mu*(s-b)*sp.diff(P,s)-P.subs(s,b)
    check(sp.degree(R,s)==count,'the spectral comparison has the required degree count='+str(count))
    zero(sp.LC(sp.Poly(R,s))-2*mu*count*sp.LC(sp.Poly(P,s)),
         'the only possible scalar multiple has the exact leading coefficient count='+str(count))
    for i,ll in enumerate(lambdas):
        si=-1/ll
        Ci=sp.prod(1-lj/ll for j,lj in enumerate(lambdas) if j!=i)
        zero((si-b)*sp.diff(P,s).subs(s,si)+(1+ll*b)*Ci,
             'the axis Hessian condition matches the spectral root identity count='+str(count)+' axis='+str(i))
    zero((R-2*mu*count*P).subs(s,b)+(2*mu*count+1)*P.subs(s,b),
         'evaluation at the origin slope gives the nonzero characteristic-zero contradiction count='+str(count))

# Distinct defects which must not be silently omitted.
repeated=qA+(x*x+y*y)**2
point={x:sp.I/sp.sqrt(2),y:0,z:0,w:0}
zeros([sp.diff(repeated,v).subs(point) for v in variables],
      'a repeated nonzero eigenvalue gives an actual critical point on the critical quadric')
zero(sp.hessian(repeated,variables).det().subs(point),
     'the repeated-eigenvalue critical quadric has a killed tangent')
resonant=qA+qB*qB-qB
zero(sp.hessian(resonant,variables).det().subs(dict.fromkeys(variables,0)),
     'resonance at the origin is a degenerate critical point')
multiple=qA+x**6/3-x**4
point={x:1,y:0,z:0,w:0}
zeros([sp.diff(multiple,v).subs(point) for v in variables],
      'a multiple root of the outer derivative gives a genuine nonzero critical point')
zero(sp.hessian(multiple,variables).det().subs(point),
     'that repeated outer root makes the critical Hessian degenerate')

# The chain-rule obstruction survives a nonlinear polynomial source automorphism.
Phi=[x,y+x*x,z+x*y,w+z*z]
a,c,d,e=sp.symbols('a c d e')
inverse=[a,c-a*a,d-a*(c-a*a),e-(d-a*(c-a*a))**2]
zeros([Phi[j].subs(dict(zip(variables,inverse)),simultaneous=True)-[a,c,d,e][j]
       for j in range(4)],'the source reparameterization has an exact polynomial inverse')
check(sp.Matrix(Phi).jacobian(variables).det()==1,
      'its determinant factor is one everywhere')
hh=qA+qB*qB
composed=sp.expand(hh.subs(dict(zip(variables,Phi)),simultaneous=True))
native_hessian=sp.hessian(hh,variables)
composed_hessian=sp.hessian(composed,variables)
native_determinants=[]
for native in [[0,0,0,0],[sp.I/sp.sqrt(2),0,0,0],[0,sp.I/(2*sp.sqrt(2)),0,0]]:
    source=[v.subs(dict(zip([a,c,d,e],native)),simultaneous=True) if hasattr(v,'subs') else v for v in inverse]
    native_det=sp.simplify(native_hessian.subs(dict(zip(variables,native))).det())
    transformed_det=sp.simplify(composed_hessian.subs(dict(zip(variables,source))).det())
    zero(transformed_det-native_det,'nonlinear source mixing preserves the critical determinant at '+str(native))
    zeros([sp.diff(composed,v).subs(dict(zip(variables,source))) for v in variables],
          'the corresponding source point remains fully critical at '+str(native))
    native_determinants.append(native_det)
check(native_determinants==[16,32,-16],
      'even one low-degree outer polynomial has three incompatible native critical determinants')

# The excluded nonsemisimple family is an actual positive control in all outer degrees.
AA=sp.Matrix([[0,1,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]])
BB=sp.diag(1,0,0,0)
T=AA.inv()*BB
check(T!=sp.zeros(4) and T*T==sp.zeros(4),
      'the boundary pair is genuinely nonzero nilpotent and not diagonalizable')
H=(t**3+t+1)**2
positive=2*x*y+z*z+w*w+H.subs(t,x*x)
zero(sp.hessian(positive,variables).det()+16,
     'the nonsemisimple boundary has actual nonzero constant Hessian in degree twelve')
shift=sp.cancel((H.subs(t,x*x)-H.subs(t,0))/(2*x))
check(sp.denom(shift)==1,'the all-degree shear divides a polynomial by a certified factor')
zero(positive-(2*x*(y+shift)+z*z+w*w+H.subs(t,0)),
     'the positive boundary has the exact polynomial quadratic normal form')
grad_positive=[sp.diff(positive,v) for v in variables]
xx=c/2
grad_inverse=[xx,(a-2*xx*sp.diff(H,t).subs(t,xx*xx))/2,d/2,e/2]
zeros([grad_positive[i].subs(dict(zip(variables,grad_inverse)),simultaneous=True)-[a,c,d,e][i]
       for i in range(4)],'its high-degree gradient inverse is exact')
zero(sp.hessian(qA+7,variables).det()-16,
     'a constant graph remains the separate nondegenerate quadratic endpoint')

# Literal graph and quadratic-symmetry certificates.
X=sp.symbols('X')
P0=t**3+t+1
Q=qA+2*X*X-3*X+5
graph=P0.subs(t,qB)
outer_graph=2*P0*P0-3*P0+5
zero(Q.subs(X,graph)-(qA+outer_graph.subs(t,qB)),
     'the genuine high-degree graph has exactly the scalar-composite native potential')
check(sp.hessian(Q,variables+[X]).det()==64,
      'the graph ambient quadratic is nondegenerate')
check(sp.degree(outer_graph,t)==6,
      'the graph application keeps the arbitrary-degree outer square')
# Use the standard nondegenerate frame for an elementary quadric symmetry.
Qstandard=x*y+z*w+X*X
coef=z*z
E=[x,y-2*X*coef-x*coef**2,z,w,X+x*coef]
Einverse=[x,y+2*X*coef-x*coef**2,z,w,X-x*coef]
ambient=variables+[X]
zero(Qstandard.subs(dict(zip(ambient,E)),simultaneous=True)-Qstandard,
     'the saturation has a literal polynomial quadratic-preserving identity')
zeros([Einverse[i].subs(dict(zip(ambient,E)),simultaneous=True)-ambient[i]
       for i in range(5)],'the displayed ambient symmetry has its full polynomial inverse')

context=json.loads((HERE/'sources/spectral-pencil-context.json').read_text())
algebra_count=len(checks)
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,
          'unchanged prior bytes: '+path)
for entry in context['newly_cached_primary_context']:
    check(hashlib.sha256((HERE/'sources'/entry['retained_file']).read_bytes()).hexdigest()==entry['sha256'],
          'selected primary context is pinned: '+entry['retained_file'])
check(context['newly_used_unproved_external_statements']==[],
      'no quartic-classification or periodicity theorem is silently imported')
check(context['unrestricted_HC4_or_JC2_resolved'] is False,
      'the full HC4 and plane Keller problems remain open')
receipt=dict(outcome='PASS',claims=['HC4MYJC0','HC4MYJC1'],
 check_count=len(checks),algebraic_and_scope_controls=algebra_count+2,checks=checks,
 checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 proof_sha256=hashlib.sha256((HERE/'SEMISIMPLE_QUADRATIC_PENCILS_AND_ALL_SOURCE_ORBITS.md').read_bytes()).hexdigest(),
 source_sha256={'sources/spectral-pencil-context.json':hashlib.sha256((HERE/'sources/spectral-pencil-context.json').read_bytes()).hexdigest()},
 dependencies=context['retained_proof_and_source_sha256'],
 primary_context=context['newly_cached_primary_context'],
 assurance='Finite exact controls of the written all-degree proof; no independent replay, formal verification or external review.',
 all_source_orbit_semisimple_nonlinearity_excluded=True,
 polynomial_graph_and_parameterization_degrees_unbounded=True,
 all_nonsemisimple_pencils_resolved=False,
 full_w_dependent_escape_excluded=False,
 unrestricted_HC4_or_JC2_resolved=False,original_goal_complete=False,
 python=platform.python_version(),sympy=sp.__version__,elapsed_seconds=round(time.monotonic()-START,3))
(HERE/'receipt-spectral-pencils.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ['outcome','check_count','algebraic_and_scope_controls','elapsed_seconds','original_goal_complete']},indent=2))
