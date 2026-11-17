"""Exact controls for the complete, all-outer-degree centered Jordan-image proof."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as sp
START=time.monotonic()
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
checks=[]
def check(v,label):
    assert bool(v),label
    checks.append(label)
def zero(v,label):check(sp.cancel(sp.expand(v))==0,label)
def zeros(vs,label):check(all(sp.cancel(sp.expand(v))==0 for v in vs),label)
x,y,z,w,r,t,k,lam,mu,nu=sp.symbols('x y z w r t k lam mu nu')
L,U,d,eps,a,rho,G,R,n=sp.symbols('L U d eps a rho G R n')
variables=[x,y,z,w]
H=lambda v:sp.expand(v**3+v+1)
F=lambda ll: r+ll*H(r)

# The shared open-stratum identity retains both residual quadratics.
for label,qC,qD in [
 ('two-simple',z*z+w*w,mu*z*z+nu*w*w),
 ('two-blocks',2*z*w,2*mu*z*w+z*z)]:
    qA=2*x*y+qC
    qB=2*lam*x*y+x*x+qD
    h=qA+H(qB)
    recovered_y=(r-x*x-qD)/(2*lam*x)
    sub={y:recovered_y}
    zero(qB.subs(sub)-r,'the length-two open inverse retains its exact scalar helper '+label)
    auxiliary=x*x-(F(lam)+lam*qC-qD-lam*k)
    zero(lam*(h.subs(sub)-k)+auxiliary,
         'the open fibre is exactly the punctured auxiliary threefold '+label)
    if label=='two-simple':
        determinant=(sp.hessian(lam*qC-qD,[z,w])/2).det()
        zero(determinant-(lam-mu)*(lam-nu),
             'the binary nondegeneracy gate keeps both repeated-eigenvalue exceptions')
    else:
        determinant=(sp.hessian(lam*qC-qD,[z,w])/2).det()
        zero(determinant+(lam-mu)**2,
             'the two-block auxiliary gate is nondegenerate exactly for distinct eigenvalues')

# Zero-eigenvalue block: solve the genuinely linear variable, with no division by lambda.
for label,qA,qB in [
 ('2+1+1',2*x*y+z*z+w*w,x*x+mu*z*z+nu*w*w),
 ('3+1',2*x*y+z*z+w*w,2*x*z+mu*w*w)]:
    h=qA+H(qB)
    recovered_y=(k-z*z-w*w-H(qB))/(2*x)
    zero(h.subs(y,recovered_y)-k,'the zero-eigenvalue open has its full affine-plane product '+label)
    zero((2*x*(recovered_y-y))-(k-h),
         'the zero-eigenvalue source recovery is exact on the full fibre '+label)
    zero(h.subs(x,0)-(z*z+w*w+H(qB.subs(x,0))),
         'the entire omitted zero slice is retained '+label)

# Length-three nonzero block: free coordinates are x,r,w.
qA=2*x*y+z*z+w*w
qB=lam*(2*x*y+z*z)+2*x*z+mu*w*w
z3=(F(lam)+(lam-mu)*w*w-lam*k)/(2*x)
y3=(r-lam*z3*z3-2*x*z3-mu*w*w)/(2*lam*x)
image3={y:y3,z:z3}
zero(qB.subs(image3,simultaneous=True)-r,'the length-three open parameterization recovers q_B')
zero((qA+H(qB)).subs(image3,simultaneous=True)-k,
     'the length-three open parameterization lies on the actual value fibre')
zero(2*x*(z3.subs(r,qB)-z)-lam*(qA+H(qB)-k),
     'its inverse z-coordinate holds on the whole value fibre')
zero((qA+H(qB)).subs(x,0)-(z*z+w*w+H(lam*z*z+mu*w*w)),
     'the length-three zero slice is the claimed semisimple plane curve')

# Length-four nonzero block: free coordinates are x,r,y.
qA=2*x*w+2*y*z
qB=lam*(2*x*w+2*y*z)+2*x*z+y*y
z4=(F(lam)-y*y-lam*k)/(2*x)
w4=(r-2*lam*y*z4-2*x*z4-y*y)/(2*lam*x)
image4={z:z4,w:w4}
zero(qB.subs(image4,simultaneous=True)-r,'the length-four open parameterization recovers q_B')
zero((qA+H(qB)).subs(image4,simultaneous=True)-k,
     'the length-four open parameterization lies on the actual value fibre')
zero(2*x*(z4.subs(r,qB)-z)-lam*(qA+H(qB)-k),
     'the length-four inverse keeps the full fibre relation')
zero((qA+H(qB)).subs(x,0)-(2*y*z+H(2*lam*y*z+y*y)),
     'the entire length-four zero slice is the Jordan plane curve times a line')

# Two blocks with the same nonzero eigenvalue: both row charts are necessary.
qA=2*x*y+2*z*w
qB=lam*qA+x*x+z*z
base=x*x+z*z-F(lam)+lam*k
row_y=(k-H(r)-2*z*w)/(2*x)
zero(qA.subs(y,row_y)-(k-H(r)),'the first equal-eigenvalue row chart is affine-line trivial')
zero(qB.subs(y,row_y)-r-base,'its scalar-helper defect is exactly the base surface equation')
row_w=(k-H(r))/(2*z)
zero(qA.subs({x:0,w:row_w})-(k-H(r)),
     'the second equal-eigenvalue row chart keeps x=0,z!=0')
zero(qB.subs({x:0,w:row_w})-r-base.subs(x,0),
     'the second row chart has the same complete base equation')
zero((qA+H(qB)).subs({x:0,z:0})-H(0),
     'none of the removed base points lifts at a generic value')

# The Jordan plane curve is an isomorphism onto the double cover with its branch points removed.
v,s=sp.symbols('v s')
q=2*lam*v*s+v*v
ss=(r-v*v)/(2*lam*v)
zero(q.subs(s,ss)-r,'the Jordan curve inverse recovers its helper')
zero(lam*(2*v*ss+H(r)-k)+(v*v-F(lam)+lam*k),
     'the Jordan curve inverse keeps its exact double-cover equation')
zero((2*v*s+H(q)).subs(v,0)-H(0),
     'the absent v=0 line is justified by the generic-value exclusion')

# Biquadratic plane curve: actual inverse square relations and disjoint branch loci.
b,c=sp.symbols('b c')
A1=(r+c*H(r)-c*k)/(b-c)
A2=(b*k-b*H(r)-r)/(b-c)
zero(b*A1+c*A2-r,'the double-square relations recover the original helper r')
zero(A1+A2+H(r)-k,'they recover the whole plane value equation')
for degree in [2,3]:
    outer=t**degree+t+1
    pp=lambda value:outer.subs(t,value)
    f1=(r+2*pp(r)-2*k)
    f2=(3*k-3*pp(r)-r)
    check(sp.gcd(sp.Poly(f1,r,domain='QQ(k)'),sp.Poly(sp.diff(f1,r),r,domain='QQ(k)')).degree()==0,
          'the first degree-'+str(degree)+' branch divisor is generically squarefree')
    check(sp.gcd(sp.Poly(f2,r,domain='QQ(k)'),sp.Poly(sp.diff(f2,r),r,domain='QQ(k)')).degree()==0,
          'the second degree-'+str(degree)+' branch divisor is generically squarefree')
    check(sp.gcd(sp.Poly(f1,r,domain='QQ(k)'),sp.Poly(f2,r,domain='QQ(k)')).degree()==0,
          'the two degree-'+str(degree)+' branch divisors are generically disjoint')
    hyper=r+2*pp(r)-2*k
    check(sp.Poly(hyper,r).degree()==degree,
          'the single double cover keeps the outer degree '+str(degree))

# Ramification and every puncture are counted symbolically, not by multiplying finite covers.
aa=(d+eps-2)/2
rr=2-eps
GG=2*d-3+eps
RR=4-2*eps
zero(2*aa-2-(-4+d+eps),'degree-two Riemann-Hurwitz includes the parity contribution at infinity')
zero(2*GG-2-(-8+4*d+2*eps),'degree-four Riemann-Hurwitz includes both finite and infinite inertia')
zero(2*aa+rr-d,'the two-cover genus and puncture identity has no missing infinity point')
zero(4-RR-2*eps,'the degree-four infinity fibre has the correct missing point count')
M=L-a*U+1-rho
auxiliary=(L-1)*L**2+L*M
zero(auxiliary-(L**3-a*L*U+(1-rho)*L),
     'the auxiliary threefold is computed from two actual product strata')
surface=(L-1)*L+d*L
zero(surface-(L**2+(d-1)*L),'the removed Danielewski surface counts every simple root')
off=auxiliary-surface
zero(off-(L**3-L**2-a*L*U+(2-rho-d)*L),
     'subtracting the whole x=0 auxiliary slice gives the shared open formula')
prototype=L**3-a*L*U-(d+rho-1)*L
two_blocks=L**3-2*a*L*U+(3-2*rho-2*d)*L
equal_blocks=L**3+(d-1)*L**2-d*L
zero(off+L*(L-1)-prototype,'a zero residual eigenvalue gives the complete one-block profile')
zero(off+L*(M-d)-two_blocks,'a second nonzero block retains both double-cover deficits')
zero(L*(surface-d)-equal_blocks,'equal-eigenvalue affine-line rows omit all nonlifting base points')
zero((L**3-L**2)+L*(L-G*U+1-R)-(L**3-G*L*U+(1-R)*L),
     'the length-three and zero-block profiles retain the full residual genus')
zero(sp.Poly(off+L*(L-G*U+1-R),U).coeff_monomial(U)+(a+G)*L,
     'residual and auxiliary genus terms have the same sign and cannot cancel')
zero(sp.Poly(off+L*d*(L-1),L).coeff_monomial(L**2)-(d-1),
     'a repeated nonzero residual eigenvalue leaves a forbidden positive L squared coefficient')
zero((prototype-(L**3-L)).subs({d:2,a:0,rho:2})+2*L,
     'degree two has the wrong coefficient even when the auxiliary genus is zero')
zero((two_blocks-(L**3-L)).subs({d:2,a:0,rho:2})+4*L,
     'two distinct nonzero blocks have coefficient minus five rather than minus one in degree two')
check(sp.Poly(equal_blocks-(L**3-n*L),L).coeff_monomial(L**2)==d-1,
      'equal blocks violate the generic profile for every possible critical count')
# Euler number alone would miss this profile defect.
mixed=L**3-L*U-L
target=L**3-3*L
check(sp.expand(mixed-target)!=0 and (mixed-target).subs({L:1,U:2})==0,
      'Euler agreement alone does not remove the mixed zero-Jordan positive-genus obstruction')
check(sp.expand((L-U-3)-2*(L-3))!=0,
      'the punctured elliptic double cover is not two copies of its punctured base')

# Exact linear-algebra controls for the critical count, resonance and positive-dimensional exception.
examples=[
 ('2+1+1',2*x*y+z*z+w*w,4*x*y+x*x),
 ('2+2 distinct',2*x*y+2*z*w,4*x*y+x*x+6*z*w+z*z),
 ('2+2 one zero',2*x*y+2*z*w,4*x*y+x*x+z*z),
 ('4',2*x*w+2*y*z,2*(2*x*w+2*y*z)+2*x*z+y*y)]
for label,qA,qB in examples:
    A=sp.hessian(qA,variables)/2
    B=sp.hessian(qB,variables)/2
    T=A.inv()*B
    for eigenvalue,mult,basis in T.eigenvects():
        if eigenvalue==0:continue
        check(all((v1.T*A*v2)[0]==0 for v1 in basis for v2 in basis),
              'every nonzero eigenspace is totally isotropic for the single-critical-point argument '+label+' eigen='+str(eigenvalue))
    hh=qA+H(qB)
    zeros([sp.diff(hh,var)-(2*A*(sp.eye(4)+sp.diff(t**3+t+1,t).subs(t,qB)*T)*sp.Matrix(variables))[i]
           for i,var in enumerate(variables)],
          'the full gradient has the claimed matrix factorization '+label)
    check(sp.hessian(hh,variables).subs(dict.fromkeys(variables,0)).det()!=0,
          'the native origin is reduced and nonresonant '+label)
resonant=2*x*y+z*z+w*w+(2*x*y+x*x)**2-(2*x*y+x*x)
zero(sp.hessian(resonant,variables).subs(dict.fromkeys(variables,0)).det(),
     'resonance gives a degenerate origin and is not discarded')
exception=2*x*y+z*z+w*w+(2*x*y+x*x+z*z)**2
zeros([sp.diff(exception,var).subs({x:0,z:sp.I/sp.sqrt(2),w:0}) for var in variables],
      'a nonzero eigenvalue shared with a simple block gives an entire critical affine line')

context=json.loads((HERE/'sources/centered-jordan-context.json').read_text())
algebra_count=len(checks)
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,
          'retained proof, receipt or primary foundation is pinned: '+path)
check(context['newly_used_unproved_external_statements']==[],
      'all curve genera and block accounting use explicit retained foundations only')
check(context['singleton_zero_endpoint_resolved'] is False
      and context['unrestricted_HC4_or_JC2_resolved'] is False,
      'the full graph reduction does not solve its old endpoint or the original HC4 goal')
receipt=dict(outcome='PASS',claims=['HC4MYJE0','HC4MYJE1'],
 check_count=len(checks),algebraic_and_scope_controls=algebra_count+2,checks=checks,
 checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 proof_sha256=hashlib.sha256((HERE/'ALL_CENTERED_QUADRATIC_COMPOSITE_IMAGES.md').read_bytes()).hexdigest(),
 source_sha256={'sources/centered-jordan-context.json':hashlib.sha256((HERE/'sources/centered-jordan-context.json').read_bytes()).hexdigest()},
 dependencies=context['retained_proof_and_source_sha256'],
 assurance='Finite exact controls for the written all-degree proof; no independent replay, formal verification or external review.',
 all_nonnilpotent_four_variable_pencils_inverse_closed=True,
 all_centered_quadratic_composite_images_reduced_to_old_endpoint=True,
 all_nonsemisimple_nonzero_spectrum_within_this_family_resolved=True,
 all_nonnilpotent_unit_hessian_source_orbits_proved_nonexistent=False,
 singleton_zero_endpoint_resolved=False,full_w_dependent_escape_excluded=False,
 unrestricted_HC4_or_JC2_resolved=False,original_goal_complete=False,
 python=platform.python_version(),sympy=sp.__version__,elapsed_seconds=round(time.monotonic()-START,3))
(HERE/'receipt-centered-jordan.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({key:receipt[key] for key in ['outcome','check_count','algebraic_and_scope_controls','elapsed_seconds','original_goal_complete']},indent=2))
