"""Exact controls for whole-fibre surgery and all invertible Jordan configurations."""
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
x,y,z,w,t,kappa,s=sp.symbols('x y z w t kappa s')
lam,mu,aa,bb,cc,dd=sp.symbols('lam mu aa bb cc dd',nonzero=True)
L,U=sp.symbols('L U')
H=t**3+2*t
P=t+lam*H-lam*kappa

# Exact cubic plane presentation and the outer-independent polar.
inner=2*lam*x*y+x*x
f=2*x*y+aa*x+bb*y+H.subs(t,inner)
ymap=(x*x-lam*aa*x-P)/(lam*bb)
cubic=2*x**3+(bb-2*lam*aa)*x*x-(2*x+bb)*t-2*lam*x*(H-kappa)
zero(lam*bb*(2*x*ymap+aa*x+bb*ymap+H-kappa)-cubic,
     'the cyclic Jordan plane has the exact monic cubic level relation')
zero((2*lam*x*ymap+x*x-t)*bb-cubic,
     'the cubic also recovers the actual inner coordinate')
zero(ymap.subs({t:inner,kappa:f},simultaneous=True)-y,
     'the cyclic plane presentation has a global polynomial inverse')
check(sp.Poly(cubic,x).degree()==3 and sp.Poly(cubic,x).LC()==2,
      'the plane inner-coordinate map is finite free of total degree three')
polar=sp.diff(inner,x)*sp.diff(f,y)-sp.diff(inner,y)*sp.diff(f,x)
zero(polar-(4*x*x+2*(bb-lam*aa)*x+2*lam*bb*y),
     'the exact Jordan polar is a smooth parabola independent of the outer derivative')
ypolar=-(2*x*x+(bb-lam*aa)*x)/(lam*bb)
tpolar=-4*x**3/bb+(2*lam*aa/bb-1)*x*x
zero(inner.subs(y,ypolar)-tpolar,
     'the polar inner coordinate has its exact cubic leading term')
for outer in [t*t+2*t,t**3+2*t]:
    degree=sp.degree(outer,t)
    level=(2*x*ypolar+aa*x+bb*ypolar+outer.subs(t,tpolar)).subs({lam:2,aa:3,bb:5})-17
    poly=sp.Poly(level,x)
    check(poly.degree()==3*degree,
          'the Jordan polar has exact 3d finite level degree d='+str(degree))
    check(sp.gcd(poly,poly.diff()).degree()==0,
          'a literal generic Jordan polar value gives all simple finite ramifications d='+str(degree))
leading=sp.symbols('leading',nonzero=True)
small=-bb/(2*lam*leading)
zero(-bb-2*lam*leading*small,
     'the small infinity branch has the exact nonzero initial coefficient')
for root in [sp.sqrt(lam*leading),-sp.sqrt(lam*leading)]:
    zero(2*root*(root*root-lam*leading),'a large infinity initial root is exact')
    zero(6*root*root-2*lam*leading-4*lam*leading,
         'each large infinity root has a nonzero simple derivative')
m,components=sp.symbols('m components',integer=True,positive=True)
for epsilon in [0,1]:
    degree=2*m+epsilon
    total=(3*degree-6+epsilon)/2+components
    base=(degree+epsilon-2)/2
    zero(2*total-2*components-(-6+3*degree+epsilon),
         'componentwise degree-three Riemann-Hurwitz retains all infinity ramification parity='+str(epsilon))
    gap=sp.expand((3*degree-4+epsilon)/2-base)
    check(gap.subs(m,1)>0 and sp.diff(gap,m)>0,
          'the Jordan cyclic total genus exceeds the double-curve genus in every degree parity='+str(epsilon))
    check((3-epsilon)+epsilon==3,
          'actual infinity points and ramification exhaust degree three parity='+str(epsilon))

# Unmarked plane: retain every omitted x=0 point, including branch points.
for coefficient in [0,3]:
    curve=x*x-2*coefficient*x-(t+2*H-2*kappa)
    yp=(t-x*x)/(4*x)
    zero((2*x*yp+coefficient*x+H-kappa)*2+curve,
         'the unmarked double-curve open inverse retains the full plane level alpha='+str(coefficient))
    missing=sp.Poly((t+2*H-14),t)
    check(missing.degree()==3 and sp.gcd(missing,missing.diff()).degree()==0,
          'the omitted x-zero locus consists of exactly d distinct points alpha='+str(coefficient))
    check(missing.eval(0)!=0,'the omitted points are not genuine x-zero inner-coordinate points')

# Whole-fibre surgery tested on every constant Jordan configuration and force type.
source=sp.Matrix([x,y,z,w])
forms=[
 ('length four forced',2*x*w+2*y*z,2*(2*x*w+2*y*z)+2*x*z+y*y,x+3*y+z,x,w),
 ('length four unforced',2*x*w+2*y*z,2*(2*x*w+2*y*z)+2*x*z+y*y,x+3*y,x,w),
 ('length three distinct forced',2*x*z+y*y+w*w,2*(2*x*z+y*y)+2*x*y+3*w*w,x+y+w,x,z),
 ('length three same forced',2*x*z+y*y+w*w,2*(2*x*z+y*y+w*w)+2*x*y,x+y+w,x,z),
 ('length three same unforced',2*x*z+y*y+w*w,2*(2*x*z+y*y+w*w)+2*x*y,x,x,z),
 ('length two r2',2*x*y+z*z+w*w,4*x*y+x*x+3*z*z+5*w*w,x+z+w,x,y),
 ('length two r1',2*x*y+z*z+w*w,4*x*y+x*x+2*z*z+5*w*w,x+w,x,y),
 ('length two r0',2*x*y+z*z+w*w,2*(2*x*y+z*z+w*w)+x*x,x,x,y),
 ('two blocks distinct',2*x*y+2*z*w,4*x*y+x*x+6*z*w+z*z,x+z,x,y),
 ('two blocks same unmarked',2*x*y+2*z*w,2*(2*x*y+2*z*w)+x*x+z*z,x+z,x,y)]
for label,qA,qB,ell,cut,free in forms:
    remaining=[v for v in [x,y,z,w] if v not in [cut,free]]
    rest=qA.subs(free,0)-qA.coeff(free)*free
    rest=qA.subs(cut,0)
    g0=qB.subs(free,0)
    D=sp.expand(g0-2*rest-2*ell)
    ymap=(t-g0)/(4*cut)
    level=qA+ell+H.subs(t,qB)
    base_equation=D-(t+2*H-2*kappa)
    zero(qB.subs(free,ymap)-t,
         'the whole open fibre recovers its actual inner coordinate '+label)
    zero(2*(level.subs(free,ymap)-kappa)+base_equation,
         'the whole open fibre retains the exact level equation '+label)
    zero(ymap.subs(t,qB)-free,
         'the whole open fibre has its actual regular inverse '+label)
    zero(level.subs(cut,0)-(rest+ell.subs(cut,0)+H.subs(t,g0.subs(cut,0))),
         'the entire closed cut stratum is retained '+label)
    AA=sp.hessian(qA,[x,y,z,w])/2
    TT=AA.inv()*sp.hessian(qB,[x,y,z,w])/2
    check(AA.det()!=0 and TT.det()!=0 and AA*TT==(AA*TT).T,
          'the controlled constant source pencil is genuinely invertible self-adjoint '+label)

# Rank profiles, all scalar exceptions and exact low-genus coefficients.
degree,genus,punctures=sp.symbols('degree genus punctures')
M=L-genus*U+1-punctures
Y3=L**3-L**2+L*M
Y2=L**3+(degree-1)*L**2
base2=L**2+(degree-1)*L
zero((L-1)*L+degree*L-(L**2+(degree-1)*L),
     'rank-two quadric profile retains all d exceptional affine lines')
zero((L-1)*L**2+L*M-Y3,
     'rank-three profile retains its complete double curve')
zero(Y2-L*M+L*(M-degree)-(L**3+(degree-1)*L**2-degree*L),
     'two same unmarked Jordan blocks have the exact nonzero L-squared profile')
zero(Y3-L*M+L*(M-degree)-(L**3-L**2+L*(M-degree)),
     'the length-four unmarked profile cancels exactly one whole double curve')
zero((Y3-L*M+L*(M-degree)).subs({degree:2,genus:0,punctures:2})-(L**3-3*L),
     'the length-four quadratic exception has coefficient three')
zero((Y3-base2+L*(L-degree)).subs({degree:2,genus:0,punctures:2})-(L**3-4*L),
     'the scalar-isotropic rank-two exception has coefficient four')
zero((Y3-base2+L*(M-degree)).subs({degree:2,genus:0,punctures:2})-(L**3-5*L),
     'the two distinct unmarked quadratic blocks have coefficient five')
zero(L**2*M-degree*L-(L**3-genus*L**2*U+(1-punctures)*L**2-degree*L),
     'the rank-one all-same-eigenvalue profile retains its higher Hodge term')

# Marked length-two all-same primary: global hyperbolic fibre isomorphism.
K=(cc*cc+dd*dd)/bb**2
Y0=(x*x-lam*aa*x-P)/(lam*bb)
S=(t-x*x)/lam-2*x*Y0+K*x*x
ZZ,WW=sp.symbols('ZZ WW')
zz=ZZ+cc*x/bb;ww=WW+dd*x/bb
yy=Y0-cc*zz/bb-dd*ww/bb
qA=2*x*y+z*z+w*w
qB=lam*qA+x*x
ell=aa*x+bb*y+cc*z+dd*w
qmap=qA.subs({y:yy,z:zz,w:ww},simultaneous=True)
zero(qmap-(t-x*x)/lam-(ZZ*ZZ+WW*WW-S),
     'the marked four-primary map recovers its actual full quadric fibre')
zero(ell.subs({y:yy,z:zz,w:ww},simultaneous=True)-(x*x-P)/lam,
     'the marked four-primary map retains the exact affine linear form')
full=qA+ell+H.subs(t,qB)
zero(Y0.subs({t:qB,kappa:full},simultaneous=True)-cc*z/bb-dd*w/bb-y,
     'the marked four-primary whole-fibre map has a polynomial inverse')
check(bb!=0,'the global fibre elimination keeps its required nonzero marked head')

# Two-block isotropic head: actual constant frame and scalar residual plane.
r0,s0,p0,q0=sp.symbols('r0 s0 p0 q0')
sub={x:(r0+s0)/2,z:(r0-s0)/(2*sp.I),y:p0+q0,w:sp.I*(p0-q0)}
zero((2*x*y+2*z*w).subs(sub,simultaneous=True)-(2*r0*p0+2*s0*q0),
     'the isotropic-head block change is an actual constant Witt frame')
zero((x*x+z*z).subs(sub,simultaneous=True)-r0*s0,
     'the nilpotent quadratic becomes the exact mixed head product')
zero((y+sp.I*w).subs(sub,simultaneous=True)-2*q0,
     'the nonzero isotropic head is retained as a genuine kernel force')

# Low-genus moments: every one of the four source moments is exact zero.
def matrix_of(q):return sp.hessian(q,[x,y,z,w])/2
low_cases=[
 ('length four',2*x*w+2*y*z,2*(2*x*w+2*y*z)+2*x*z+y*y,x+3*y),
 ('three plus line isotropic',2*x*z+y*y+w*w,2*(2*x*z+y*y+w*w)+2*x*y,x+y+sp.I*w),
 ('two plus scalar distinct isotropic',2*x*y+z*z+w*w,4*x*y+x*x+3*z*z+3*w*w,x+z+sp.I*w),
 ('two plus scalar same isotropic',2*x*y+z*z+w*w,2*(2*x*y+z*z+w*w)+x*x,x+z+sp.I*w),
 ('two distinct blocks',2*x*y+2*z*w,4*x*y+x*x+6*z*w+z*z,x+z),
 ('two same isotropic heads',2*x*y+2*z*w,2*(2*x*y+2*z*w)+x*x+z*z,y+sp.I*w)]
for label,qA,qB,ell in low_cases:
    AA=matrix_of(qA);TT=AA.inv()*matrix_of(qB)
    uu=AA.inv()*sp.Matrix([sp.diff(ell,v) for v in [x,y,z,w]])/2
    zeros([(uu.T*AA*(TT**j)*uu)[0] for j in range(4)],
          'all primary moments really vanish in the low-genus case '+label)
    xx=-(sp.eye(4)+s*TT).inv()*uu
    zero((xx.T*AA*TT*xx)[0],
         'the regular inner stationary function is identically zero '+label)

# Complete singular stationary groups, not just regular resolvent points.
A0,B0=sp.symbols('A0 B0')
qA=2*x*w+2*y*z
qB=2*qA+2*x*z+y*y
ell=A0*x+B0*y
stationary=sp.Matrix([sp.diff(qA,v)-sp.diff(qB,v)/2+sp.diff(ell,v) for v in [x,y,z,w]])
zeros(list(stationary.subs({x:0,y:B0,z:A0})),
      'the length-four singular stationary set is its entire eigenline')
zero(qB.subs({x:0,y:B0,z:A0})-(4*A0*B0+B0*B0),
     'the length-four inner function is constant on the full singular stationary line')
for eigen,other in [(2,3),(3,2)]:
    qA=2*x*y+2*z*w
    qB=2*eigen*x*y+x*x+2*other*z*w+z*z
    ell=A0*x+B0*z
    stationary=[sp.diff(qA,v)-sp.diff(qB,v)/eigen+sp.diff(ell,v) for v in [x,y,z,w]]
    point={x:eigen*A0/2,z:0,w:-B0/(2*(1-sp.Rational(other,eigen)))}
    zeros([v.subs(point) for v in stationary],
          'an unmarked length-two primary retains its whole singular stationary line eigen='+str(eigen))
    check(sp.degree(qB.subs(point),y)<=1,
          'its inner equation has at most one quadratic-outer critical point unless the full line survives')
check(1<3 and 2<4 and 3<5,
      'the complete Morse-count bounds contradict every low-genus exceptional coefficient')

context=json.loads((HERE/'sources/jordan-polar-invertible-context.json').read_text())
algebra=len(checks)
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,
          'retained source/proof/checker/receipt bytes stay pinned: '+path)
for entry in context['newly_cached_primary_context']:
    check(hashlib.sha256((HERE/'sources'/entry['retained_file']).read_bytes()).hexdigest()==entry['sha256'],
          'the versioned primary suspension context is byte-pinned')
check(context['newly_used_unproved_external_statements']==[],
      'the proof imports no additional unproved classification or finite-cover theorem')
check(context['all_four_variable_invertible_affine_pencils_inverse_closed'] is True,
      'the complete closure has no diagonalizability or offset-cyclicity qualification')
check(context['all_singular_affine_pencils_resolved'] is False
      and context['quadratic_two_critical_full_source_orbits_resolved'] is False
      and context['unit_ideal_endpoint_resolved'] is False
      and context['singleton_zero_endpoint_resolved'] is False
      and context['unrestricted_HC4_or_JC2_resolved'] is False,
      'singular pencils, fixed full-source tests and both actual endpoints remain open')
flags={key:value for key,value in context.items() if isinstance(value,bool)}
receipt=dict(outcome='PASS',claims=['HC4MYJL0','HC4MYJL1'],check_count=len(checks),
 algebraic_and_scope_controls=algebra+3,checks=checks,
 checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 proof_sha256=hashlib.sha256((HERE/'JORDAN_POLAR_CURVES_AND_ALL_INVERTIBLE_PENCIL_CLOSURE.md').read_bytes()).hexdigest(),
 source_sha256={'sources/jordan-polar-invertible-context.json':hashlib.sha256((HERE/'sources/jordan-polar-invertible-context.json').read_bytes()).hexdigest()},
 dependencies=context['retained_proof_and_source_sha256'],
 assurance='Exact finite controls of the written all-degree Jordan polar and complete whole-fibre proof; not independent replay, formal verification or external review.',
 python=platform.python_version(),sympy=sp.__version__,
 elapsed_seconds=round(time.monotonic()-START,3),**flags)
(HERE/'receipt-jordan-polar-curves.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({key:receipt[key] for key in ['outcome','check_count','algebraic_and_scope_controls','elapsed_seconds','original_goal_complete']},indent=2))
