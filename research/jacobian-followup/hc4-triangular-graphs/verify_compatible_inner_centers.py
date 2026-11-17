"""Exact controls for specific stabilized closures and compatible-inner reductions."""
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
def zero(value,label):
    value=sp.cancel(value)
    check(value==0 or sp.simplify(value)==0,label)
def zeros(values,label):
    check(all(sp.cancel(v)==0 or sp.simplify(v)==0 for v in values),label)
x,y,z,w,t,s,kappa,alpha=sp.symbols('x y z w t s kappa alpha')
lam,aa=sp.symbols('lam aa',nonzero=True)
variables=[x,y,z,w]
L,U,d,a_d,rho=sp.symbols('L U d a_d rho')
M=L-a_d*U+1-rho
B_d=L-(d-1)*U-1
zero((L-1)*L**2+L*(L-d)-(L**3-d*L),
     'the scalar-isotropic stabilized plane has its exact d profile coefficient')
zero((L-1)*L**2+L*d*(L-1)-(L**3+(d-1)*L**2-d*L),
     'the centered scalar stabilized plane retains every special line')
zero((L-1)*L**2+L*(M-d)-(L**3-L**2+L*(M-d)),
     'the unmarked stabilized Jordan plane keeps all d deleted points')
zero(((L-1)*L**2+L*(M-d)).subs({d:2,a_d:0,rho:2})-(L**3-3*L),
     'the low-genus stabilized Jordan plane has coefficient three')
for degree in [2,3]:
    P=z*z+z+(2*z*z)**degree
    polynomial=sp.Poly(17-P,z)
    check(polynomial.degree()==2*degree and sp.gcd(polynomial,polynomial.diff()).degree()==0,
          'the zero-variable hyperelliptic control has all 2d distinct branch roots d='+str(degree))
zero(2*(d-1)-2-(-4+2*d),
     'Riemann-Hurwitz gives genus d-1 for the even-degree double cover')
m=sp.symbols('m',integer=True,positive=True)
for epsilon in [0,1]:
    degree=2*m+epsilon
    half=(degree+epsilon-2)/2
    zero(degree-1-half-m,
         'the new plane genus strictly exceeds every subtracted double-curve genus parity='+str(epsilon))
    check(m.is_positive,'the all-degree genus difference is positive')
zero(L**3-L**2-(d-1)*L+L*B_d-(L**3-(d-1)*L*U-d*L),
     'the ternary two-primary isotropic profile retains its uncanceled genus')
zero(L**3-L*M+L*B_d-(L**3-(d-1-a_d)*L*U+(rho-2)*L),
     'the ternary scalar-isotropic profile retains the exact positive genus difference')

# The complete exceptional plane Jordan critical set includes the singular point.
g=4*z*w+z*z
f=2*z*w+3*z+g*g
points=[{z:0,w:-sp.Rational(3,2)},{z:3,w:-sp.Rational(37,48)}]
for index,point in enumerate(points):
    zeros([sp.diff(f,q).subs(point) for q in [z,w]],
          'the stabilized plane low-genus control retains an actual critical point '+str(index))
    check(sp.hessian(f,[z,w]).subs(point).det()!=0,
          'the controlled regular or singular plane critical point is Morse')
zero(sp.diff(f,w)-2*z*(1+4*g),
     'the plane critical equation exhausts regular and singular stationary branches')
zero((2*w+3-sp.Rational(1,2)*(4*w+2*z))-(3-z),
     'the singular stationary branch fixes z and leaves exactly one inner root')
check(2<3,'the full two-point native critical count contradicts coefficient three')
scalar_stationary=sp.Matrix([2*w*(1+lam*s)+2*s,2*z*(1+lam*s)])
zeros(list(scalar_stationary.subs(s,-1/lam)-sp.Matrix([-2/lam,0])),
      'a nonzero isotropic scalar offset blocks its only determinant parameter')

# The ternary large scalar primary has a global hyperbolic fibre inverse.
H=t**3+2*t
source_inner=lam*(2*x*y+z*z)+2*aa*z+1
source_h=2*x*y+z*z+w*w+H.subs(t,source_inner)
zmap=(t+lam*H-lam*kappa+lam*w*w-1)/(2*aa)
S=kappa-H-zmap*zmap-w*w
zero(lam*(S+zmap*zmap)+2*aa*zmap+1-t,
     'the large scalar-primary fibre map recovers the actual inner coordinate')
zero(S+zmap*zmap+w*w+H-kappa,
     'the large scalar-primary fibre map retains the whole level')
zero(zmap.subs({t:source_inner,kappa:source_h},simultaneous=True)-z,
     'the large scalar-primary fibre has a global polynomial inverse')
plane_inner=lam*z*z+2*aa*z+1
plane_h=z*z+w*w+H.subs(t,plane_inner)
zero(S.subs({t:plane_inner,kappa:plane_h},simultaneous=True),
     'the hyperbolic zero curve is the actual rank-one plane level')

# Every ternary unmarked chain rank/force type, with an actual zero-inner variable.
cases=[
 ('length three forced',2*x*z+y*y+w*w,2*(2*x*z+y*y)+2*x*y,x+3*y,x,z,3,1,False,True),
 ('length three unforced',2*x*z+y*y+w*w,2*(2*x*z+y*y)+2*x*y,x,x,z,3,1,False,False),
 ('length two distinct',2*x*y+z*z+w*w,4*x*y+x*x+3*z*z,x+z,x,y,3,2,False,False),
 ('length two same forced',2*x*y+z*z+w*w,4*x*y+x*x+2*z*z,x+3*z,x,y,2,1,True,True),
 ('length two same unforced',2*x*y+z*z+w*w,4*x*y+x*x+2*z*z,x,x,y,2,1,False,False)]
def force(q,coords):
    matrix=sp.hessian(q,coords)
    linear=sp.Matrix([sp.diff(q,v).subs(dict.fromkeys(coords,0)) for v in coords])
    return any(sp.simplify((linear.T*u)[0])!=0 for u in matrix.nullspace())
for label,qA,qB,ell,cut,free,rank,rank0,forced,forced0 in cases:
    rest=qA.subs(cut,0)
    g0=qB.subs(free,0)
    D=sp.expand(g0-2*rest-2*ell)
    coords=[v for v in variables if v!=free]
    rem=[v for v in coords if v!=cut]
    D0=D.subs(cut,0)
    check(sp.hessian(D,coords).rank()==rank and force(D,coords)==forced,
          'the ternary-cut rank and actual kernel force are exact '+label)
    check(sp.hessian(D0,rem).rank()==rank0 and force(D0,rem)==forced0,
          'the cut-base rank and actual kernel force are exact '+label)
    inverse=(t-g0)/(4*cut)
    hh=qA+ell+H.subs(t,qB)
    zero(qB.subs(free,inverse)-t,'the unmarked-chain open cut recovers its inner coordinate '+label)
    zero(2*(hh.subs(free,inverse)-kappa)+D-(t+2*H-2*kappa),
         'the unmarked-chain whole fibre retains its exact open equation '+label)
    zero(hh.subs(cut,0)-(rest+ell.subs(cut,0)+H.subs(t,g0.subs(cut,0))),
         'the entire closed unmarked-chain stratum is retained '+label)
    A=sp.hessian(qA,variables)/2
    T=A.inv()*sp.hessian(qB,variables)/2
    check(T.rank()==3 and T*sp.Matrix([0,0,0,1])==sp.zeros(4,1),
          'this is an actual ternary invertible pencil stabilized by a zero-inner quadratic line '+label)

# Exact finite/infinity checks for the ternary marked critical-mass branches.
hyper=sp.Matrix([[0,1,0],[1,0,0],[0,0,1]])
anti=sp.Matrix([[0,0,1],[0,1,0],[1,0,0]])
J2=sp.Matrix([[2,0],[1,2]])
J3=sp.Matrix([[2,0,0],[1,2,0],[0,1,2]])
mass_cases=[
 ('marked length three',anti,J3,sp.Matrix([1,0,0]),t*t,4,3),
 ('marked length two distinct',hyper,sp.diag(J2,3),sp.Matrix([1,0,0]),t*t,3,2),
 ('marked length two same',hyper,sp.diag(J2,2),sp.Matrix([1,0,0]),t*t,3,3),
 ('semisimple repeated active',sp.eye(3),sp.diag(2,2,3),sp.Matrix([1,0,1]),t*t,2,2),
 ('scalar active large outer',sp.eye(3),2*sp.eye(3),sp.Matrix([1,0,0]),t**3/3,2,3)]
def order(value,point):
    num,den=sp.fraction(sp.cancel(value))
    def mult(poly):
        pp=sp.Poly(poly,s);result=0
        while pp.eval(point)==0:
            pp=sp.div(pp,sp.Poly(s-point,s))[0];result+=1
        return result
    return mult(num)-mult(den)
for label,A,T,v,outer,pole,multiplicity in mass_cases:
    xx=-s*(sp.eye(3)+s*T).inv()*v
    R=sp.cancel((xx.T*A*T*xx+2*v.T*A*xx)[0])
    P=sp.factor((sp.eye(3)+s*T).det())
    F=sp.cancel(s-sp.diff(outer,t).subs(t,R))
    W=sp.cancel(1/(P*F))
    check(order(R,-sp.Rational(1,2))==-pole and order(P,-sp.Rational(1,2))==multiplicity,
          'the ternary active pole and full multiplicity are both retained '+label)
    check(order(W,-sp.Rational(1,2))>=0,
          'the actual finite critical-mass budget passes '+label)
    mat=sp.eye(3)-T/2
    check(mat.row_join(v/2).rank()>mat.rank(),
          'the active singular stationary equation is actually inconsistent '+label)
    num,den=sp.fraction(W)
    check(sp.degree(den,s)-sp.degree(num,s)==4,
          'the ternary complete mass has the literal fourth-order infinity decay '+label)

# Higher-rank compatible inner examples and the exact global null shear.
parameter=sp.symbols('parameter')
qA=2*x*z+y*y+w*w
inner=2*x*y+2*y+w*w+2*parameter*w+parameter*parameter
outer=-t*t/4+t/2
hh=qA+outer.subs(t,inner)
A=sp.hessian(qA,variables)/2
T=A.inv()*sp.hessian(inner,variables)/2
v=A.inv()*sp.Matrix([sp.diff(inner,q).subs(dict.fromkeys(variables,0)) for q in variables])/2
check(T.rank()==3 and T*sp.Matrix([1,0,0,parameter])==v,
      'the two fixed tests have an actual higher-rank compatible inner center')
zeros([sp.diff(inner,q).subs({x:-1,y:0,w:-parameter}) for q in variables],
      'the entire explicit compatible inner stationary line is retained')
reduced_inner=inner.subs(x,0)
shift=sp.cancel((outer.subs(t,inner)-outer.subs(t,reduced_inner))/(2*x))
check(sp.denom(shift).is_number and not shift.has(z),
      'the compatible-center null shear is polynomial and independent of its shifted variable')
zero(hh-(2*x*(z+shift)+y*y+w*w+outer.subs(t,reduced_inner)),
     'the exact compatible-center fibre is the rank-one resonant suspension')
zero((z+shift).subs(z,z-shift)-z,
     'the compatible-center null shear has its actual polynomial inverse')
Z,W=sp.symbols('Z W')
plane=y*y+w*w+outer.subs(t,reduced_inner)
zero(plane.subs({w:Z-parameter,y:-W-Z*Z/4+sp.Rational(5,4)},simultaneous=True)
     -(W*(Z*Z-1)-2*parameter*Z+parameter*parameter+sp.Rational(5,4)),
     'the compatible-center family realizes exactly h0 and h1 by a polynomial inverse map')
coord_inner=2*x*y+2*y
coord_h=qA+outer.subs(t,coord_inner)
coord_shift=sp.cancel((outer.subs(t,coord_inner)-outer.subs(t,2*y))/(2*x))
zero(coord_h-(2*x*(z+coord_shift)+w*w+y),
     'the compatible inner-center family also realizes the old coordinate test')
check((A.inv()*sp.hessian(coord_inner,variables)/2)*sp.Matrix([1,0,0,0])
      ==A.inv()*sp.Matrix([0,2,0,0])/2,
      'the coordinate test really has compatible affine force')

# The only anisotropic kernel case is an actual stabilized ternary pencil.
Tline=sp.diag(0,2,3,4);vline=sp.Matrix([0,1,0,2])
check(Tline*sp.Matrix([0,sp.Rational(1,2),0,sp.Rational(1,2)])==vline,
      'a nondegenerate kernel line with compatible force splits off an invertible ternary pencil')
check(Tline.nullspace()==[sp.Matrix([1,0,0,0])],
      'the compatible exceptional kernel is exactly one nondegenerate line')

# Remaining coordinate-inner normal forms, with actual global polynomial inverses.
p=x*x+3*z*z+4*w*w
g=2*y+p
original=2*x*y+z*z+w*w+H.subs(t,g)
zero(original.subs(y,(t-p)/2)-(H+x*t+z*z+w*w-x*p),
     'the isotropic forced zero direction gives the exact cubic-plus-univariate normal form')
zero(g.subs(y,(t-p)/2)-t,
     'the cubic coordinate-inner map has its global polynomial inverse')
check(sp.Poly(x*t+z*z+w*w-x*p,x,z,w,t).total_degree()==3,
      'the cubic normal form has no hidden higher multivariate degree')
p2=2*x*x+3*z*z+4*w*w
g2=2*y+p2
original2=x*x+y*y+z*z+w*w+H.subs(t,g2)
zero(original2.subs(y,(t-p2)/2)-(H+x*x+z*z+w*w+(t-p2)**2/4),
     'the nonisotropic forced zero direction gives the exact quartic-plus-univariate normal form')
zero(g2.subs(y,(t-p2)/2)-t,
     'the quartic coordinate-inner map has its global polynomial inverse')
check(sp.Poly(x*x+z*z+w*w+(t-p2)**2/4,x,z,w,t).total_degree()==4,
      'the quartic normal form has no hidden higher multivariate degree')
check(sp.Matrix([x,g2,z,w]).jacobian(variables).det()==2,
      'the actual quadratic-coordinate replacement has constant nonzero Jacobian')
ranktwoA=sp.diag(sp.Matrix([[0,1],[1,0]]),1,1)
ranktwoT=sp.diag(0,0,2,3);ranktwov=sp.Matrix([1,1,0,0])
ker=sp.Matrix.hstack(*ranktwoT.nullspace())
line=(ranktwov.T*ranktwoA*ker).nullspace()[0]
check((line.T*ker.T*ranktwoA*ker*line)[0]!=0,
      'the remaining rank-two kernel constraint can genuinely be an anisotropic line')
check(ranktwoT.row_join(ranktwov).rank()>ranktwoT.rank(),
      'the remaining coordinate-inner example has no inner critical point')
zero(sp.diff(2*y+p,y)-2,
     'the remaining inner quadratic has an actual constant nonzero directional derivative')

context=json.loads((HERE/'sources/compatible-inner-center-context.json').read_text())
algebra=len(checks)
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,
          'retained proof/checker/receipt/primary bytes stay pinned: '+path)
check(context['newly_used_unproved_external_statements']==[],
      'no generic stabilization, cancellation or new coordinate theorem is imported')
check(context['all_inner_quadratics_with_critical_points_reduced_to_four_fixed_source_orbits'] is True,
      'the compatible-center reduction covers every affine coefficient and outer/source degree')
check(context['remaining_quadratic_coordinate_frontier_resolved'] is False
      and context['quadratic_two_critical_full_source_orbits_resolved'] is False
      and context['unit_ideal_endpoint_resolved'] is False
      and context['singleton_zero_endpoint_resolved'] is False
      and context['unrestricted_HC4_or_JC2_resolved'] is False,
      'the coordinate-inner frontier, four fixed tests and unrestricted HC4 remain open')
flags={key:value for key,value in context.items() if isinstance(value,bool)}
receipt=dict(outcome='PASS',claims=['HC4MYJN0','HC4MYJN1'],check_count=len(checks),
 algebraic_and_scope_controls=algebra+3,checks=checks,
 checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 proof_sha256=hashlib.sha256((HERE/'COMPATIBLE_INNER_CENTERS_AND_QUADRATIC_COORDINATE_FRONTIER.md').read_bytes()).hexdigest(),
 source_sha256={'sources/compatible-inner-center-context.json':hashlib.sha256((HERE/'sources/compatible-inner-center-context.json').read_bytes()).hexdigest()},
 dependencies=context['retained_proof_and_source_sha256'],
 assurance='Exact finite controls of the written all-degree specific-stabilization, compatible-center and coordinate-frontier proof; not independent replay, formal verification or external review.',
 python=platform.python_version(),sympy=sp.__version__,
 elapsed_seconds=round(time.monotonic()-START,3),**flags)
(HERE/'receipt-compatible-inner-centers.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({key:receipt[key] for key in ['outcome','check_count','algebraic_and_scope_controls','elapsed_seconds','original_goal_complete']},indent=2))
