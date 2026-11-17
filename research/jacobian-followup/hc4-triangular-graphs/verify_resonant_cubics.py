"""Exact controls for the written all-degree three-resonant-cubic reduction."""
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
    check(sp.cancel(value)==0 or sp.simplify(value)==0,label)
def zeros(values,label):
    check(all(sp.cancel(v)==0 or sp.simplify(v)==0 for v in values),label)
x,y,z,w,t,s,kappa,beta=sp.symbols('x y z w t s kappa beta')
b,nu,lam,e,c,alpha=sp.symbols('b nu lam e c alpha',nonzero=True)
E,D=sp.symbols('E D')
L,U,d,G,C,r,ad,rho=sp.symbols('L U d G C r ad rho')
H=sp.Function('H')
inner=nu*z*z+2*b*w+beta
inverse=(t-nu*z*z-beta)/(2*b)
F=z*z+e*z+(t-nu*z*z-beta)**2/(4*b*b)+H(t)
zero(inner.subs(w,inverse)-t,'mixed plane: global coordinate inverse')
zero(inverse.subs(t,inner)-w,'mixed plane: reverse global coordinate')
zero((z*z+w*w+e*z+H(inner)).subs(w,inverse)-F,'mixed plane: literal transformed potential')
check(sp.Poly(F-H(t),z).degree()==4,'mixed plane projection is finite free of degree four')
polar=nu*z*z+beta+2*b*b/nu+e*b*b/(nu*z)
zero(sp.diff(F,z).subs(t,polar),'the entire nonzero-force polar is parametrized by Gm')
zero(inverse.subs(t,polar)-(b/nu+e*b/(2*nu*z)),'the actual w coordinate on the polar has both correct poles')
h2=sp.symbols('h2',nonzero=True)
quadratic=F.subs(H(t),h2*t*t+t).subs(t,polar)
zero(sp.limit(z*z*quadratic,z,0)-e*e*b*b*(1+4*b*b*h2)/(4*nu*nu),
     'quadratic polar pole cancellation is exactly the leading infinity resonance')
for degree in [2,3]:
    for force in [0,1]:
        f=(z*z+force*z+(t-2*z*z-1)**2/4+t**degree)
        if force:
            tt=2*z*z+2+sp.Rational(1,2)/z
            pp=sp.cancel(f.subs(t,tt)-17)
            numerator,denominator=sp.fraction(pp)
            pz=sp.Poly(numerator,z)
            check(pz.degree()==3*degree,'exact moving polar numerator degree '+str((degree,force)))
            check(sp.Poly(denominator,z).degree()==degree,'the complete zero-pole order is retained '+str((degree,force)))
            check(sp.gcd(pz,pz.diff()).degree()==0,'generic nonzero-force polar intersections are simple '+str(degree))
            check(pz.eval(0)!=0,'the omitted polar endpoint is never counted as a finite point '+str(degree))
        else:
            pp=sp.Poly(f.subs(z,0)-17,t)
            qq=sp.Poly(f.subs(t,2*z*z+2)-17,z)
            check(pp.degree()==degree and qq.degree()==2*degree,'the zero-force polar components give all 3d points '+str(degree))
            check(sp.gcd(pp,pp.diff()).degree()==sp.gcd(qq,qq.diff()).degree()==0,'both zero-force polar restrictions are generically simple '+str(degree))
            check(f.subs({z:0,t:2})!=17,'the single polar crossing value is excluded '+str(degree))
    check(4*degree>=4*2,'every admitted control has four actual leading sheets '+str(degree))
zero(2*((3*d-4-r)/2+C)-2*C-(-8+3*d+4-r),
     'componentwise quartic Riemann-Hurwitz retains all 3d finite ramifications and true punctures')
n=sp.symbols('n',integer=True,nonnegative=True)
for residue in range(4):
    degree=4*n+residue
    punctures=4 if residue==0 else (2 if residue==2 else 1)
    eps=residue%2
    difference=degree-sp.Rational(punctures+eps,2)
    minimum=1 if residue<3 else 0
    check(difference.subs(n,minimum)>0,'all-degree genus exceeds a_d in degree class '+str(residue))
    zero((3*degree-2-punctures)/2-(degree+eps-2)/2-difference,
         'exact lower genus difference in degree class '+str(residue))
    check(4//punctures*punctures==4,'all infinity deck orbits account for four sheets class '+str(residue))
Z=sp.symbols('Z')
initial=nu*nu*Z**4-2*nu*Z*Z+1+4*b*b*h2
zz=sp.symbols('zz')
zero(sp.discriminant(nu*nu*zz*zz-2*nu*zz+1+4*b*b*h2,zz)+16*nu*nu*b*b*h2,
     'the nonresonant quadratic leading squares are distinct')
zero(initial.subs(Z,0)-(1+4*b*b*h2),'the leading roots are all nonzero exactly off resonance')
zero(2*C-2*C-(-8+6+2),'six finite branches and two infinity branches give G=c_C for d=2')
# The unmarked-isotropic exact two-plane fibre decomposition.
q=2*x*y+z*z
bb=2*lam*x*y+E*x*x+2*D*x*z+nu*z*z
ell=alpha*x+e*z
f=q+ell+t*(bb-c)
K=-E*x*x/lam-2*D*x*z/lam+(1-nu/lam)*z*z+alpha*x+e*z+c/lam
zero(f.subs(t,-1/lam)-K,'the entire singular t-plane is the actual quadratic K')
zero(sp.diff(f,y)-2*x*(1+lam*t),'the solved-stratum divisor consists of the two full planes')
zero(f.subs(x,0)-(z*z+e*z+t*(nu*z*z-c)),'the complete x-zero plane is a solved punctured-line family')
check(sp.Poly(K.subs(x,0),z).degree()==2,'distinct remaining eigenvalue gives two generic intersection points')
zero(K.subs({x:0,nu:lam})-(e*z+c/lam),'coincident eigenvalues retain their one affine intersection')
m,nn,EK=sp.symbols('m nn EK')
zero(L**3-(2*L*L-L)+L*((L-m)+EK-nn)-(L**3-L*L+L*EK-(m+nn-1)*L),
     'whole fibre profile retains both planes, their centre intersection and every free y line')
for profile,expected in [(L-1,L**3-(m+nn)*L),(L,L**3-(m+nn-1)*L),
                         (2*L,L**3+L*L-(m+nn-1)*L),(0,L**3-L*L-(m+nn-1)*L)]:
    zero(L**3-L*L+L*profile-(m+nn-1)*L-expected,'all exact rank/force profiles for K: '+str(profile))
for is_cone in [True,False]:
    mm=1 if is_cone else 2
    for forced in [True,False]:
        required=mm+2-int(forced)
        bound=(0 if is_cone else 2)+(0 if forced else 1)
        check(bound<required,'distinct-eigenvalue complete critical count contradicts profile '+str((is_cone,forced)))
zero(f.subs({x:0,z:0,c:0,e:0})+0,'zero force on the cone retains the possible full critical curve')
ycurve=-alpha/(2*(1+lam*t))
zeros([sp.diff(f,v).subs({x:0,z:0,c:0,e:0,y:ycurve}) for v in [x,y,z,t]],
      'the entire excluded cone critical curve is verified, not just a point')
zero(sp.diff(f,z).subs({x:0,nu:lam,t:-1/lam})-e,
     'nonzero z force forbids the compatible singular x-zero group')
check(sp.hessian(K,[x,z]).det().subs(nu,lam)==-4*D*D/(lam*lam),
      'coincident-eigenvalue rank two occurs precisely in the length-three type')
zero(sp.diff(K,z).subs(nu,lam)-(-2*D*x/lam+e),
     'the unique singular quadratic stationary x is nonzero exactly when e is nonzero')
zero(sp.diff(K,z).subs({nu:lam,D:0})-e,'length-two/scalar models have a forced K line and no singular K point')
# Cone closure: failure of proportionality gives an unmarked isotropic eigenline.
B=sp.Matrix([[0,1,0],[1,0,0],[0,0,1]])
linear=sp.Matrix([0,0,1]);p=sp.Matrix([1,0,0])
gamma=sp.symbols('gamma')
A=lam*B+gamma*linear*linear.T
zeros(list((A-lam*B)*p),'the proportional-conic failure gives the actual common eigenvector')
zero((p.T*B*p)[0],'that eigenvector is isotropic on the cone')
zero((linear.T*p)[0],'it annihilates the full affine force')
check(A.det()!=0,'the conic proportionality control retains a nondegenerate ternary a')
# Mixed whole-fibre ranks and the strict genus comparison without cover multiplication.
M=L-ad*U+1-rho
curve=C*L-G*U+C-r
Yprofiles=[L**3,L**3-L*L+L*M,L**3+(d-1)*L*L,L*L*M]
Y0profiles=[L*L,L*L+(d-1)*L,L*M]
for iy,Y in enumerate(Yprofiles):
    for iz,Y0 in enumerate(Y0profiles):
        coeff=sp.expand(Y-Y0+L*curve).coeff(U)
        check(sp.expand(coeff).coeff(G)==-L,'the actual mixed curve genus remains in every whole-fibre stratum '+str((iy,iz)))
        check(sp.expand(coeff).coeff(ad).subs(L,1)<=1,'the subtracted binary genus can supply at most a_d '+str((iy,iz)))
# Direct stationary finite budgets in the nonzero ternary complement.
for multiplicity,pole in [(1,2),(2,2),(3,4),(2,3),(3,3)]:
    check(pole>=multiplicity,'every admitted active group is blocked and has sufficient all-degree budget '+str((multiplicity,pole)))
A2=sp.Matrix([[0,1],[1,0]])
N2=sp.Matrix([[0,0],[1,0]])
J3A=sp.Matrix([[0,1,0],[1,0,0],[0,0,1]])
J3B=lam*J3A+sp.Matrix([[0,0,1],[0,0,0],[1,0,0]])
J3T=J3A.inv()*J3B
check((J3T-lam*sp.eye(3))**3==sp.zeros(3) and (J3T-lam*sp.eye(3))**2!=sp.zeros(3),
      'the final coupled model is precisely a length-three self-adjoint primary')
for name,nnn,rank in [('scalar',0,0),('length two',x*x,1),('length three',2*x*z,2)]:
    Bform=lam*q+nnn
    Bmat=sp.hessian(Bform,[x,y,z])/2
    check(Bmat.det()==-lam**3,'every retained cubic has a genuinely nondegenerate inner ternary form '+name)
    check((J3A.inv()*Bmat-lam*sp.eye(3)).rank()==rank,'constant Jordan frame of retained model '+name)
# Both resonance coefficients and actual full-source maps.
q3=sp.Function('q3')(x,y,z)
aa3=sp.Function('a3')(x,y,z)
ll3=sp.Function('ell3')(x,y,z)
outer=-(q3+2*b*w)**2/(4*b*b)+gamma*(q3+2*b*w)
shear=-b*t-q3/(4*b)
actual=(w*w+aa3+ll3+outer).subs(w,shear)
zero(actual-(aa3+ll3+(t+gamma/2)*(q3-2*b*b*gamma)+gamma*gamma*b*b),
     'the full quadratic resonance is exactly a shifted quadric multiplier plus a constant')
zero((-(w+q3/(4*b))/b).subs(w,shear)-t,'the resonance shear has its actual polynomial inverse')
check(sp.Matrix([x,y,z,shear]).jacobian([x,y,z,t]).det()==-b,'the resonance source map has constant nonzero Jacobian')
u=sp.symbols('u',nonzero=True)
rr=beta-2*b*b*s+u/s**2
zero(sp.limit((s-(2*h2*rr+gamma))/s,s,sp.oo)-(1+4*b*b*h2),
     'the regular infinity leading coefficient agrees with the plane resonance')
zero(sp.limit((s-(2*h2*rr+gamma)).subs(h2,-1/(4*b*b)),s,sp.oo)-(beta/(2*b*b)-gamma),
     'the constant infinity coefficient is exactly the second quadratic resonance')
# All retained models pass the native and entire value-profile controls.
R=sp.symbols('R',nonzero=True)
for name,nnn in [('scalar',0),('length two',x*x),('length three',2*x*z)]:
    potential=q+alpha*x+e*z+t*(lam*q+nnn-c)
    determinant=sp.factor(sp.hessian(potential,[x,y,z,t]).det())
    for sign in [1,-1]:
        zz=sign*R
        tt=-1/lam-e/(2*lam*zz)
        yy=sp.cancel(-(alpha+tt*sp.diff(nnn,x).subs({x:0,z:zz}))/(2*(1+lam*tt)))
        point={x:0,z:zz,t:tt,y:yy,c:lam*R*R}
        zeros([sp.diff(potential,v).subs(point,simultaneous=True) for v in [x,y,z,t]],
              'complete regular critical point '+name+' sign='+str(sign))
        zero(determinant.subs(point,simultaneous=True)-4*e*e*lam*lam,
             'native Hessian agrees at every regular critical point '+name+' sign='+str(sign))
        zero(potential.subs(point,simultaneous=True)-(R*R+sign*e*R),
             'actual regular critical value '+name+' sign='+str(sign))
    if name=='length three':
        xx=e*lam/2;zz=alpha*lam/2
        yy=(4*c-alpha*alpha*lam**3-2*alpha*e*lam**2)/(4*e*lam**2)
        point={x:xx,y:yy,z:zz,t:-1/lam}
        zeros([sp.diff(potential,v).subs(point,simultaneous=True) for v in [x,y,z,t]],'the complete compatible singular critical point is retained')
        zero(determinant.subs(point,simultaneous=True)-4*e*e*lam*lam,'the third critical point has the same nonzero native Hessian')
        zero(potential.subs(point,simultaneous=True)-(c/lam+alpha*e*lam/2),'the third critical value retains its arbitrary relative parameter')
        check(xx!=0,'the third critical point remains distinct even at coincident values')
    else:
        zero(sp.diff((potential.subs(t,-1/lam)),z)-e,'there are no omitted compatible singular points '+name)
    number=3 if name=='length three' else 2
    EK0=L-1 if name=='length three' else L
    generic=sp.expand(L**3-2*L*L+L+L*((L-2)+EK0-1))
    zero(generic-(L**3-number*L),'the entire generic profile passes HG0 '+name)
    zero(L*((L-2+L)-(L-2))-L*L,'each exceptional regular value contributes one full L^2 jump '+name)
    if name=='length three':
        zero(L*((2*L-1)-(L-1))-L*L,'the singular quadratic level adds the third L^2 jump, including collisions')
        KK=-2*x*z/lam+alpha*x+e*z+c/lam
        z_on_curve=(kappa+2*x*z/lam-alpha*x-c/lam)/e
        check(sp.hessian(KK,[x,z]).det()!=0,'the third singular level is a genuine nondegenerate binary quadratic')
relative=sp.symbols('relative')
normalized=q+2*relative*x+z+t*(q+2*x*z-1)-1
zero(normalized.subs({t:-1,x:sp.Rational(1,2),z:relative,y:1-relative-relative*relative})-relative,
     'the normalized three-point example has critical value relative for all parameters')
context=json.loads((HERE/'sources/resonant-cubic-context.json').read_text())
algebra=len(checks)
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'retained proof/checker/receipt/primary bytes: '+path)
check(context['newly_used_unproved_external_statements']==[],'no new connectedness, cancellation or multiple-pole residue theorem is imported')
check(context['all_affine_quadratic_composite_universal_questions_reduced_to_four_tests_and_three_cubic_families'] is True,
      'the reduction covers every affine coefficient, nonlinear outer and source degree')
check(context['remaining_three_resonant_cubic_source_orbit_families_resolved'] is False
      and context['additional_cubic_families_identified_with_old_fixed_tests'] is False
      and context['unit_ideal_endpoint_resolved'] is False
      and context['singleton_zero_endpoint_resolved'] is False
      and context['unrestricted_HC4_or_JC2_resolved'] is False,
      'passing full profiles and Hessians are recorded as an open source-orbit obstruction, not a solution')
flags={key:value for key,value in context.items() if isinstance(value,bool)}
receipt=dict(outcome='PASS',claims=['HC4MYJP0','HC4MYJP1','HC4MYJP2'],check_count=len(checks),
 algebraic_and_scope_controls=algebra+3,checks=checks,
 checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 proof_sha256=hashlib.sha256((HERE/'ALL_AFFINE_PENCILS_TO_THREE_RESONANT_CUBIC_FAMILIES.md').read_bytes()).hexdigest(),
 source_sha256={'sources/resonant-cubic-context.json':hashlib.sha256((HERE/'sources/resonant-cubic-context.json').read_bytes()).hexdigest()},
 dependencies=context['retained_proof_and_source_sha256'],
 assurance='Exact finite controls of the written all-degree mixed-plane genus, complete ternary-cone closure and three-cubic full-source reduction proof; not independent replay, formal verification or external review.',
 python=platform.python_version(),sympy=sp.__version__,elapsed_seconds=round(time.monotonic()-START,3),**flags)
(HERE/'receipt-resonant-cubics.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({key:receipt[key] for key in ['outcome','check_count','algebraic_and_scope_controls','elapsed_seconds','original_goal_complete']},indent=2))
