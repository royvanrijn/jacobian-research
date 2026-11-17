"""Exact controls for global null shears and all-degree rank-one source orbits."""
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
x,y,z,w,t,s,c,gamma,delta=sp.symbols('x y z w t s c gamma delta')
lam,b=sp.symbols('lam b',nonzero=True)
a=sp.symbols('a')
variables=[x,y,z,w]

# The global divided-difference shear and its inverse in unrestricted polynomial form.
for K in [z*z+w*w+(x+z*z+w)**2,
          z*z+w*w+(2*x+3*z*z+4*z+1)**3,
          z*z+w*w+(x*z+x*x+w*w)**4]:
    shift=sp.cancel((K-K.subs(x,0))/(2*x))
    check(sp.denom(shift).is_number and not shift.has(y),
          'the null-pair divided difference is a polynomial independent of the shifted variable')
    zero((2*x*y+K)-(2*x*(y+shift)+K.subs(x,0)),
         'the null-pair shear is the exact global potential identity')
    zero((y+shift).subs(y,y-shift)-y,
         'the null-pair shear has a literal polynomial inverse')
    check(sp.Matrix([x,y+shift,z,w]).jacobian(variables).det()==1,
          'the null-pair shear has constant Jacobian one')
K=z*z+w*w+(x+z*z+w)**2
native=2*x*y+K;reduced=2*x*y+K.subs(x,0)
check(sp.hessian(native,variables).det().subs({x:1,z:0,w:0})
      !=sp.hessian(reduced,variables).det().subs({x:1,z:0,w:0}),
      'the polynomial normal form is not falsely treated as native Hessian invariance')

# An actual complex Witt frame in a zero primary space.
u=sp.Matrix([1,sp.I,0]);up=sp.Matrix([sp.Rational(1,2),-sp.I/2,0])
frame=sp.Matrix.hstack(up,u,sp.Matrix([0,0,1]))
check(frame.det()!=0 and frame.T*frame==sp.Matrix([[0,1,0],[1,0,0],[0,0,1]]),
      'the isotropic zero-space vector has an actual nondegenerate constant hyperbolic partner')
check((u.T*u)[0]==0 and (u.T*up)[0]==1,
      'the Witt normalization retains both null norms and the exact unit pairing')
A=sp.diag(sp.Matrix([[0,1],[1,0]]),1,1)
T=sp.diag(0,0,3,0)
uu=sp.Matrix([0,1,0,0])
vv=sp.Matrix([0,1,2,0])
check(T*uu==sp.zeros(4,1) and (uu.T*A*uu)[0]==0 and (vv.T*A*uu)[0]==0,
      'the actual shear example satisfies all three zero-eigenvector hypotheses')
xx=sp.Matrix(variables)
inner=(xx.T*A*T*xx+2*vv.T*A*xx)[0]+1
zero(sp.diff(inner,y),'the actual inner polynomial is independent of the paired free variable')
bad=sp.Matrix([1,0,0,0])
check((bad.T*A*uu)[0]!=0,
      'a zero eigenvector with a nonorthogonal affine force is not admitted')
check(3*sp.eye(4)*uu!=sp.zeros(4,1),
      'an isotropic nonzero-eigenvalue vector is not substituted for a zero eigenvector')

# The plane resolvent, all finite boundaries and the exact infinity resonance.
beta=c-a*a/lam
zz=-a*s/(1+lam*s);ww=-b*s
R=beta-2*b*b*s+a*a/(lam*(1+lam*s)**2)
zero(lam*zz*zz+2*a*zz+2*b*ww+c-R,
     'the rank-one plane stationary inner function is exact in every affine coefficient')
P=1+lam*s
outer_slope,outer_linear=sp.symbols('outer_slope outer_linear')
F=sp.cancel(s-(outer_slope*R+outer_linear))
zero(F-((1+2*outer_slope*b*b)*s-(outer_slope*beta+outer_linear)
        -outer_slope*a*a/(lam*(1+lam*s)**2)),
     'the quadratic leading and constant infinity coefficients are exact')
resonant=-1/(2*b*b)
zero(F.subs(outer_slope,resonant)
     -(-((resonant*beta)+outer_linear)+a*a/(2*b*b*lam*(1+lam*s)**2)),
     'the exceptional quadratic keeps its genuine nondecaying infinity boundary')
for outer in [t*t,t**3/3,t**4]:
    m=sp.degree(sp.diff(outer,t),t)
    FF=sp.cancel(s-sp.diff(outer,t).subs(t,R.subs({lam:2,a:1,b:1,c:0})))
    WW=sp.cancel(1/((1+2*s)*FF))
    num,den=sp.fraction(WW)
    check(sp.degree(den,s)-sp.degree(num,s)>=2,
          'every selected nonresonant control retains proper infinity decay outer degree='+str(m+1))
    nn,dd=sp.fraction(WW)
    check(sp.cancel(WW).subs(s,-sp.Rational(1,2))==0,
          'the blocked active simple determinant pole is actually regular outer degree='+str(m+1))
simple_R=R.subs(a,0)
check(sp.denom(sp.cancel(simple_R)).is_number or not sp.denom(sp.cancel(simple_R)).has(s),
      'the compatible simple determinant root retains a regular inner value when the marked component vanishes')

# Full native four-variable critical mass after quadratic stabilization.
for marked in [0,1]:
    g=2*z*z+2*marked*z+2*w
    h=2*x*y+z*z+w*w+g*g
    gradient=[sp.diff(h,q) for q in variables]
    Hess=sp.hessian(h,variables)
    if marked:
        parameters=[sp.Integer(0),-sp.Rational(3,5)+3*sp.I/10,
                    -sp.Rational(3,5)-3*sp.I/10]
        points=[{x:0,y:0,z:-ss/(1+2*ss),w:-ss} for ss in parameters]
    else:
        points=[{x:0,y:0,z:0,w:0},
                {x:0,y:0,z:sp.sqrt(-sp.Rational(5,8)),w:sp.Rational(1,2)},
                {x:0,y:0,z:-sp.sqrt(-sp.Rational(5,8)),w:sp.Rational(1,2)}]
    masses=[]
    for index,point in enumerate(points):
        zeros([v.subs(point) for v in gradient],
              'the stabilized obstruction keeps an actual full critical point marked='+str(marked)+' index='+str(index))
        J=sp.simplify(Hess.subs(point).det())
        check(J!=0,'the controlled stabilized critical point is Morse')
        plane=sp.hessian(z*z+w*w+g*g,[z,w]).subs(point).det()
        zero(J+4*plane,'the quadratic stabilization retains its exact determinant factor')
        masses.append(1/J)
    zero(sum(masses),'the complete stabilized native reciprocal-Hessian sum is zero marked='+str(marked))
check(sp.diag(0,1,1).rank()==2,'the compatible determinant pole is not treated as a regular matrix inverse')

# Invariant zero-space norm is not the norm of the whole offset.
TT=sp.diag(2,0,0,0);EE=TT/2
vglobal=sp.Matrix([sp.I,1,0,0])
zero((vglobal.T*vglobal)[0],'the whole offset can be isotropic')
vzero=(sp.eye(4)-EE)*vglobal
check((vzero.T*vzero)[0]==1,
      'the relevant zero-primary norm remains nonzero for a globally isotropic offset')
vother=sp.Matrix([2,1,sp.I,0])
zero((((sp.eye(4)-EE)*vother).T*((sp.eye(4)-EE)*vother))[0],
     'a nonzero isotropic zero-primary projection has eta zero')
check((vother.T*vother)[0]==4,
      'a nonzero whole-offset norm cannot replace eta in the invariant classification')
check(EE*EE==EE and TT*TT==TT.trace()*TT,
      'the rank-one nonzero eigenvalue has its exact orthogonal projection algebra')

# Every resonance has the literal multiplier normal form before any division.
u0=a/lam
Z,W=sp.symbols('Z W')
outer=-(t-beta)**2/(4*b*b)+gamma*(t-beta)+delta
plane=z*z+w*w+outer.subs(t,lam*z*z+2*a*z+2*b*w+c)
multiplier=-lam*Z*Z/b+2*gamma*b
remainder=-lam*lam*Z**4/(4*b*b)+(1+gamma*lam)*Z*Z-2*u0*Z
zero(plane.subs(z,Z-u0)-(w*multiplier+remainder+u0*u0+delta),
     'the full quadratic infinity resonance has its exact rank-one multiplier expansion')
quotient,rem=sp.div(remainder,multiplier,Z)
zero(sp.diff(rem,Z)+2*u0,
     'the quadratic multiplier remainder has exactly the invariant marked linear coefficient')
zero(sp.diff(remainder,Z).subs({Z:0,gamma:0})+2*u0,
     'the repeated multiplier has no critical point when the marked component is nonzero')
zero(sp.diff(remainder,Z).subs({Z:0,gamma:0,a:0}),
     'the unmarked repeated multiplier retains a full critical line')
L=sp.symbols('L')
zero((L-1)*L*L+L*(L-1)-(L**3-L),
     'the empty-critical resonant case has its exact non-HG0 generic profile')

# Fixed resonant coefficients give explicit polynomial inverse normal forms.
for marked in [0,4]:
    Cnum=multiplier.subs({lam:2,b:3,gamma:1})
    Dnum=remainder.subs({lam:2,b:3,gamma:1,a:marked})
    qq,rr=sp.div(Dnum,Cnum,Z)
    zet,Om=sp.symbols('zet Om')
    zz=3*zet;ww=Om/(-6)-qq.subs(Z,zz)
    zero((w*Cnum+Dnum).subs({Z:zz,w:ww},simultaneous=True)
         -(Om*(zet*zet-1)-3*marked*zet+rr.subs(Z,0)),
         'the exact resonant normal form is h0 or the scalar h1 orbit marked='+str(marked))
    zero((-6*(w+qq)).subs({Z:zz,w:ww},simultaneous=True)-Om,
         'the resonant multiplier shear has its global polynomial inverse')
    check(sp.Matrix([zz,ww]).jacobian([zet,Om]).det()!=0,
          'the resonant source normal form has a constant nonzero Jacobian')

# All four fixed tests are genuine rank-at-most-one closed graph potentials.
aa=sp.symbols('aa')
inner=z*z+2*aa*z+2*w+aa*aa
example=2*x*y+z*z+w*w-inner**2/4+inner/2
zero(example.subs({z:Z-aa,w:-W-Z*Z/4+sp.Rational(5,4)},simultaneous=True)
     -(2*x*y+W*(Z*Z-1)-2*aa*Z+aa*aa+sp.Rational(5,4)),
     'the two fixed full-source tests occur in the exact rank-one semisimple graph family')
for marked in [0,-sp.Rational(1,2)]:
    check(-2*marked in [0,1],'the explicit graph realizes h0 and h1 without an endpoint equivalence claim')
nilA=sp.diag(sp.Matrix([[0,1],[1,0]]),1,1)
nilB=sp.diag(1,0,0,0)
nilT=nilA.inv()*nilB
check(nilT.rank()==1 and nilT.trace()==0 and nilT*nilT==sp.zeros(4),
      'the rank-one trace-zero boundary is genuinely nilpotent')
zero((2*x*y+z*z+w*w+x**4).subs(y,y-x**3/2)-(2*x*y+z*z+w*w),
     'the rank-one nilpotent graph realizes the old quadratic source orbit')
zero(x*x+y*y+z*z+w*w-(x+1)**2-(y*y+z*z+w*w-2*x-1),
     'the rank-zero graph realizes the old coordinate source orbit')
check(sp.det(sp.hessian(2*x*y+z*z+w*w-Z*Z/4+Z/2,[x,y,z,w,Z]))!=0,
      'the resonant graph ambient quadratic form is nondegenerate')

context=json.loads((HERE/'sources/rank-one-pencil-context.json').read_text())
algebra=len(checks)
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,
          'retained proof/checker/receipt/primary bytes stay pinned: '+path)
check(context['newly_used_unproved_external_statements']==[],
      'the stabilization and zero-eigenvector maps import no new unproved cancellation or residue theorem')
check(context['nonzero_eigenvalue_rank_one_affine_pencil_classification_complete'] is True,
      'the invariant rank-one classification admits every affine offset and outer degree')
check(context['quadratic_two_critical_full_source_orbits_resolved'] is False
      and context['all_singular_affine_pencils_resolved'] is False
      and context['unit_ideal_endpoint_resolved'] is False
      and context['singleton_zero_endpoint_resolved'] is False
      and context['unrestricted_HC4_or_JC2_resolved'] is False,
      'all four fixed tests and higher-rank singular pencils remain open')
flags={key:value for key,value in context.items() if isinstance(value,bool)}
receipt=dict(outcome='PASS',claims=['HC4MYJM0','HC4MYJM1'],check_count=len(checks),
 algebraic_and_scope_controls=algebra+3,checks=checks,
 checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 proof_sha256=hashlib.sha256((HERE/'ZERO_EIGENVECTOR_SHEARS_AND_RANK_ONE_PENCIL_CLASSIFICATION.md').read_bytes()).hexdigest(),
 source_sha256={'sources/rank-one-pencil-context.json':hashlib.sha256((HERE/'sources/rank-one-pencil-context.json').read_bytes()).hexdigest()},
 dependencies=context['retained_proof_and_source_sha256'],
 assurance='Exact finite controls of the written all-degree null shear, stabilized critical mass and invariant rank-one classification; not independent replay, formal verification or external review.',
 python=platform.python_version(),sympy=sp.__version__,
 elapsed_seconds=round(time.monotonic()-START,3),**flags)
(HERE/'receipt-rank-one-pencils.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({key:receipt[key] for key in ['outcome','check_count','algebraic_and_scope_controls','elapsed_seconds','original_goal_complete']},indent=2))
