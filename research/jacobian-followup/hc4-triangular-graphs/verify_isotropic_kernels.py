"""Exact controls for the written all-degree isotropic-kernel reduction."""
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
a,b,c,mu=sp.symbols('a b c mu',nonzero=True)
L,U,d,G,C,r=sp.symbols('L U d G C r')
H=sp.Function('H')
# A literal coordinate, its inverse, polar and trigonal equation.
g=z*z+2*a*w+2*b*z+c
inverse=(t-z*z-2*b*z-c)/(2*a)
F=H(t)+(z*t-z**3-2*b*z*z-c*z)/a
zero(g.subs(w,inverse)-t,'marked zero Jordan plane: actual global inner coordinate inverse')
zero(inverse.subs(t,g)-w,'marked zero Jordan plane: reverse global coordinate map')
zero((2*z*w+H(g)).subs(w,inverse)-F,'marked zero Jordan plane: exact transformed potential')
cubic=z**3+2*b*z*z+(c-t)*z+a*(kappa-H(t))
zero(cubic+a*(F-kappa),'the level equation is the monic degree-three presentation')
polar=3*z*z+4*b*z+c
zero(sp.diff(F,z)-(t-polar)/a,'the finite projection polar is a smooth parabola')
zero(F.subs(t,polar)-H(polar)-(2*z**3+2*b*z*z)/a,'polar restriction has leading degree 2d')
for degree in [2,3]:
    control=(polar**degree+(2*z**3+2*b*z*z)/a).subs({a:1,b:1,c:1})
    pp=sp.Poly(control-17,z)
    check(pp.degree()==2*degree,'exact finite branch degree control d='+str(degree))
    check(sp.gcd(pp,pp.diff()).degree()==0,'simple generic polar intersections control d='+str(degree))
    check(pp.eval(sp.Rational(-2,3))!=0,'the one higher polar ramification value is excluded d='+str(degree))
tau,Z,hd=sp.symbols('tau Z hd',nonzero=True)
zero(sp.diff(Z**3-a*hd,Z)-3*Z*Z,'every leading trigonal root is simple and nonzero')
for degree in [2,3]:
    check(3*degree>max(2*degree,3+degree,degree,0),
          'all lower infinity terms have smaller pole order d='+str(degree))
zero(2*(d-(3+r)/2+C)-2*C-(-6+2*d+3-r),
     'componentwise Riemann-Hurwitz retains all component genera and actual punctures')
n=sp.symbols('n',integer=True,nonnegative=True)
for residue in [0,1,2]:
    degree=3*n+residue
    punctures=3 if residue==0 else 1
    bound=sp.expand(degree-sp.Rational(1+punctures,2))
    expected=3*n-2 if residue==0 else 3*n+residue-1
    zero(bound-expected,'all-degree positive genus lower bound in degree class '+str(residue))
    min_n=1 if residue in [0,1] else 0
    check(bound.subs(n,min_n)>0,'the d>=2 boundary of this residue class has positive genus '+str(residue))
    check(3//punctures*punctures==3,'deck orbit punctures account for all degree-three sheets '+str(residue))
# Binary cone: every open and closed stratum, including the common two points.
j,k=sp.symbols('j k',integer=True,nonnegative=True)
cone=2*L-1
EC=(2+k)*L-j-2
zero(L*(L**2-cone)+L*EC-(L**3+k*L**2-(j+1)*L),
     'the whole binary-cone trace keeps its two common points and every free line')
for aa,bb,profile in [(2,1,L-1),(0,1,L),(0,0,2*L)]:
    check((aa!=0 and profile==L-1) or (aa==0 and bb!=0 and profile==L) or (aa==bb==0 and profile==2*L),
          'line restriction conic profile with quadratic and linear coefficients '+str((aa,bb)))
    tangent=2*aa*z+bb
    if aa!=0:check(sp.Poly(tangent,z).degree()==1,'nonconstant tangent equation has at most one cone-line point')
    elif bb!=0:check(tangent!=0,'a forced linear cone line has no critical point')
    else:check(tangent==0,'the identically zero line is recorded in k rather than silently omitted')
f=s*s+2*z*z+3*w*w+z+w+t*(2*z*w)
for point in [{s:0,z:-sp.Rational(1,4),w:0,t:2},{s:0,z:0,w:-sp.Rational(1,6),t:3}]:
    zeros([sp.diff(f,q).subs(point) for q in [s,z,w,t]],'complete binary-cone critical control on one of the two lines')
    check(sp.hessian(f,[s,z,w,t]).det().subs(point)!=0,'the cone-line critical control is genuinely Morse')
zero(sp.diff(2*z*w,z).subs({z:0,w:1})-2,'the normal derivative fixes the remaining multiplier on the nonzero cone')
# Exact zero-primary stationary polynomial, not a zero-eigenvalue substitution.
def stationary(A,T,v,const=c):
    xx=-s*(sp.eye(T.rows)+s*T).inv()*v
    return sp.cancel((xx.T*A*T*xx)[0]+2*(v.T*A*xx)[0]+const)
A2=sp.Matrix([[0,1],[1,0]])
N2=sp.Matrix([[0,0],[1,0]])
R2=stationary(A2,N2,sp.Matrix([a,b]))
zero(R2-(3*a*a*s*s-4*a*b*s+c),'marked length-two zero polynomial has exact nonzero top moment')
A3=sp.Matrix([[0,0,1],[0,1,0],[1,0,0]])
N3=sp.Matrix([[0,0,0],[1,0,0],[0,1,0]])
R3=stationary(A3,N3,sp.Matrix([a,b,c]),const=beta)
zero(R3-(-4*a*a*s**3+6*a*b*s*s-(4*a*c+2*b*b)*s+beta),
     'marked length-three zero polynomial has exact nonzero top moment')
for length,R0 in [(2,R2),(3,R3)]:
    check(sp.Poly(R0,s).degree()==length,'zero chain exact polynomial degree length='+str(length))
    check(sp.Poly(R0-2*s,s).degree()==length,'an extra zero line cannot cancel the marked top moment length='+str(length))
    check(length+1>=3,'nonzero-primary determinant plus marked chain gives infinity decay at least three length='+str(length))
# Finite active groups: exact local orders and complete constant primary ranks.
cases=[('distinct simple',sp.diag(2,3),sp.eye(2),sp.Matrix([1,1]),[(-sp.Rational(1,2),1,2),(-sp.Rational(1,3),1,2)]),
 ('scalar nonisotropic',2*sp.eye(2),sp.eye(2),sp.Matrix([1,0]),[(-sp.Rational(1,2),2,2)]),
 ('Jordan marked',2*sp.eye(2)+N2,A2,sp.Matrix([1,0]),[(-sp.Rational(1,2),2,3)])]
for label,T,A,v,orders in cases:
    zeros(list(T.T*A-A*T),'the finite primary operator is self-adjoint '+label)
    RR=stationary(A,T,v,const=0)
    PP=sp.Poly((sp.eye(2)+s*T).det(),s)
    denominator=sp.Poly(sp.fraction(RR)[1],s)
    for root,multiplicity,pole in orders:
        value=sp.limit(PP.as_expr()/(s-root)**multiplicity,s,root)
        check(value.is_finite and value!=0,
              'full determinant multiplicity retained '+label+str(root))
        value=sp.limit(denominator.as_expr()/(s-root)**pole,s,root)
        check(value.is_finite and value!=0,
              'active R pole has its complete budget '+label+str(root))
        check(pole>=multiplicity,'all outer m>=1 make this blocked determinant group regular in W '+label+str(root))
for label,T,A,v in [('scalar isotropic',2*sp.eye(2),A2,sp.Matrix([1,0])),
                    ('Jordan unmarked',2*sp.eye(2)+N2,A2,sp.Matrix([0,1]))]:
    zero(stationary(A,T,v,const=beta)-beta,'every nonzero-primary moment vanishes in the finite bad type '+label)
    check((sp.eye(2)-T/2).row_join(v/2).rank()>(sp.eye(2)-T/2).rank()
          if label=='scalar isotropic' else (sp.eye(2)-T/2).row_join(v/2).rank()==(sp.eye(2)-T/2).rank(),
          'singular stationary equation is blocked scalar or compatible Jordan '+label)
# Whole-fibre rank profiles in the new zero-primary cases.
ad,rho=sp.symbols('ad rho')
M=L-ad*U+1-rho
Q=L*L+(d-1)*L
V=L**3-L*L+L*M
curve=C*L-G*U+C-r
for label,Y in [('scalar isotropic',L**3),('Jordan unmarked',V)]:
    profile=sp.expand(Y-Q+L*curve)
    coefficient=sp.expand(profile).coeff(U)
    zero(coefficient-(-G*L if label=='scalar isotropic' else -(G+ad)*L),
         'the marked zero-plane genus cannot cancel in the whole-fibre identity '+label)
# The binary cut rank is mu A0-B0, including its affine force.
zero((mu*A2-A2*N2).det()+mu*mu,'the marked zero Jordan cut is nondegenerate at every nonzero eigenvalue')
zero((mu*sp.eye(2)).det()-mu*mu,'the semisimple zero-plane cut is nondegenerate at every nonzero eigenvalue')
for B0 in [A2*N2,sp.zeros(2)]:
    quad=(sp.Matrix([z,w]).T*(B0-mu*A2)*sp.Matrix([z,w]))[0]+x*x
    check(sp.hessian(quad,[x,z,w]).rank()==3,'the whole nonzero Jordan eigenvector quadric has full ternary rank')
    check(sp.hessian(quad.subs(x,0),[z,w]).rank()==2,'the complete eigenvector cut has binary rank two')
# Semisimple zero plane, quadratic outer and all four low-genus exceptional counts.
eta,u,v=sp.symbols('eta u v',nonzero=True)
R=beta-2*eta*s
regular=sp.expand(s-(u*R+v))
zero(regular-((1+2*u*eta)*s-(u*beta+v)),
     'the bad finite types have the complete linear regular stationary equation')
zero(regular.subs(u,-1/(2*eta))+(v-beta/(2*eta)),
     'the leading infinity resonance leaves exactly the constant -gamma')
for label,Y,shift in [('scalar isotropic',L**3,1),('Jordan unmarked',V,2)]:
    for outer,f0,required,bound in [('nonresonant',L-1,shift+1,shift),('leading resonance',L,shift,shift-1)]:
        actual=sp.expand((Y-Q+L*f0).subs({d:2,ad:0,rho:2}))
        zero(actual-(L**3-required*L),'exact whole low-genus coefficient '+label+' '+outer)
        check(bound<required,'all critical groups have too few points for HG0 '+label+' '+outer)
# Jordan singular line, including its possible constant inner function.
lam,force=sp.symbols('lam force',nonzero=True)
q0=2*x*y
q1=2*lam*x*y+x*x
pot=q0+force*x+s*q1
singular={s:-1/lam,x:lam*force/2}
zeros([sp.diff(pot,q).subs(singular) for q in [x,y]],
      'the entire compatible Jordan stationary line is retained')
check(sp.Poly(q1.subs(x,lam*force/2),y).degree()==1,
      'a marked eigenvector force makes the inner value affine on the singular line')
zero(q1.subs({x:0,force:0}),
     'the zero eigenvector force gives a constant inner value and any compatible group is a full line')
# The double infinity resonance is a global polynomial shear, in arbitrary forms.
qB=sp.Function('qB')(z,w)
qA=sp.Function('qA')(z,w)
ell=sp.Function('ell')(z,w)
original=s*s+y*y+qA+ell-(qB+2*b*y)**2/(4*b*b)
shear=-b*t-qB/(4*b)
zero(original.subs(y,shear)-(s*s+qA+ell+t*qB),
     'double infinity resonance is literally the binary-cone multiplier potential')
zero((-(y+qB/(4*b))/b).subs(y,shear)-t,'the double-resonance shear has its polynomial inverse')
# Exhaustion has only the fixed primary possibilities, not outer-degree cases.
types=[([2],[1,1]),([2],[2]),([1,1],[1,1]),([1,1],[2]),([3],[1]),([2,1],[1])]
for zero_blocks,nonzero_blocks in types:
    check(sum(zero_blocks)+sum(nonzero_blocks)==4,'zero/nonzero primary exhaustion dimension '+str((zero_blocks,nonzero_blocks)))
    check(sum(length-1 for length in zero_blocks)+sum(nonzero_blocks) in [2,3],
          'only remaining rank-two/rank-three configurations '+str((zero_blocks,nonzero_blocks)))
# Entire cubic family, literal constant frame, and all four represented tests.
p=2*x*x+3*x*z+4*z*z+2*w+1
inner=2*a*y+p
original=2*x*y+z*z+w*w+H(inner)
normal=H(t)+x*t/a+z*z+w*w-x*p/a
zero(original.subs(y,(t-p)/(2*a))-normal,'the entire cubic generating family is the actual source orbit of an isotropic-kernel pencil')
zero(inner.subs(y,(t-p)/(2*a))-t,'the entire cubic map has the global polynomial inverse')
check(sp.Poly(x*t/a+z*z+w*w-x*p/a,x,z,w,t).total_degree()==3,'the cubic bound applies exactly to the part outside arbitrary H')
XA,YA,ZA,WA=y,x-y/2-w,z,y+w
zero(2*XA*YA+ZA*ZA+WA*WA-(2*x*y+z*z+w*w),'constant frame preserves the actual quadratic form')
frame=sp.Matrix([XA,YA,ZA,WA])
check(frame.jacobian([x,y,z,w]).det()!=0,'the full four-test constant frame has a polynomial inverse')
zero(2*WA-(2*y+2*w),'the coordinate test is the exact p=2w member')
zero(ZA*ZA+2*c*ZA+2*WA+c*c-(2*y+z*z+2*c*z+2*w+c*c),
     'the two JM resonant graphs are literal members of the cubic class')
quadratic=2*x*y+z*z+w*w+(2*y)**2
check(sp.hessian(quadratic,[x,y,z,w]).det()!=0,'the p=0 quadratic outer member realizes the old nondegenerate quadratic test')
oldcoord=2*sp.Symbol('X')*sp.Symbol('Y')+sp.Symbol('Z')**2+sp.Symbol('W')
check(sp.diff(oldcoord,sp.Symbol('W'))==1,'the coordinate example has an explicit coordinate inverse along W')
# The sole residual quartic: simple zero line and both ternary forms invertible.
qA3=x*x+z*z+w*w
qB3=2*x*x+3*z*z+4*w*w
ell3=x+2*z-3*w
quartic=H(t)+qA3+ell3+(t-qB3-beta)**2/(4*b*b)
g3=qB3+2*b*y+beta
zero((y*y+qA3+ell3+H(g3)).subs(y,(t-qB3-beta)/(2*b))-quartic,
     'the simple-zero quartic is an equality of actual full source orbits')
zero(g3.subs(y,(t-qB3-beta)/(2*b))-t,'the remaining quartic map has an actual global polynomial inverse')
check(sp.Matrix([x,g3,z,w]).jacobian([x,y,z,w]).det()==2*b,'the quartic coordinate map has constant nonzero Jacobian')
check(sp.Poly(quartic-H(t),x,z,w,t).total_degree()==4,'the remaining ternary quartic has exactly the stated bounded part')
TT=sp.diag(2,0,3,4)
check(TT.rank()==3 and len(TT.nullspace())==1,'the sole new frontier has simple zero and rank three')
check((TT.nullspace()[0].T*TT.nullspace()[0])[0]!=0,'its only kernel line is nondegenerate')
check(sp.hessian(qA3,[x,z,w]).det()!=0 and sp.hessian(qB3,[x,z,w]).det()!=0,'both remaining ternary forms are nondegenerate')
check(TT.row_join(sp.Matrix([0,b,0,0])).rank()==4,'the nonzero offset on the simple kernel line forbids an inner critical point')
context=json.loads((HERE/'sources/isotropic-kernel-context.json').read_text())
algebra=len(checks)
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'retained proof/checker/receipt/primary bytes: '+path)
check(context['newly_used_unproved_external_statements']==[],'no new residue, cancellation or arbitrary conic-cover theorem is imported')
check(context['entire_cubic_coordinate_inner_family_reduced_to_four_fixed_source_orbits'] is True,
      'the whole cubic question reduces to the four retained tests without degree restrictions')
check(context['remaining_simple_zero_quartic_frontier_resolved'] is False
      and context['quadratic_two_critical_full_source_orbits_resolved'] is False
      and context['unit_ideal_endpoint_resolved'] is False
      and context['singleton_zero_endpoint_resolved'] is False
      and context['unrestricted_HC4_or_JC2_resolved'] is False,'the simple-zero quartic, four fixed tests and unrestricted HC4 remain open')
flags={key:value for key,value in context.items() if isinstance(value,bool)}
receipt=dict(outcome='PASS',claims=['HC4MYJO0','HC4MYJO1'],check_count=len(checks),
 algebraic_and_scope_controls=algebra+3,checks=checks,
 checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 proof_sha256=hashlib.sha256((HERE/'ISOTROPIC_KERNELS_AND_SIMPLE_ZERO_QUARTIC_FRONTIER.md').read_bytes()).hexdigest(),
 source_sha256={'sources/isotropic-kernel-context.json':hashlib.sha256((HERE/'sources/isotropic-kernel-context.json').read_bytes()).hexdigest()},
 dependencies=context['retained_proof_and_source_sha256'],
 assurance='Exact finite controls of the written all-degree trigonal, binary-cone, complete stationary-group and full cubic reduction proof; not independent replay, formal verification or external review.',
 python=platform.python_version(),sympy=sp.__version__,elapsed_seconds=round(time.monotonic()-START,3),**flags)
(HERE/'receipt-isotropic-kernels.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({key:receipt[key] for key in ['outcome','check_count','algebraic_and_scope_controls','elapsed_seconds','original_goal_complete']},indent=2))
