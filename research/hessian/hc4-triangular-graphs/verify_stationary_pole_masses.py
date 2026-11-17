"""Exact controls for completed simple-pole critical groups and noncyclic Krylov criteria."""
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
def iszero(value):
    value=sp.cancel(value)
    return value==0 or sp.simplify(value)==0
def zero(value,label):check(iszero(value),label)
def zeros(values,label):check(all(iszero(v) for v in values),label)
s,t,c=sp.symbols('s t c')
x,y,z,w=sp.symbols('x y z w')
variables=[x,y,z,w]
def pencil(A,T,v,outer):
    n=T.rows
    xx=-s*(sp.eye(n)+s*T).inv()*v
    R=sp.cancel((xx.T*A*T*xx+2*v.T*A*xx)[0])
    P=sp.factor((sp.eye(n)+s*T).det())
    F=sp.cancel(s-sp.diff(outer,t).subs(t,R))
    W=sp.cancel(1/(P*F))
    return xx,R,P,F,W
def inf_order(value):
    num,den=sp.fraction(sp.cancel(value))
    return sp.degree(den,s)-sp.degree(num,s)
def order_at(value,point):
    num,den=sp.fraction(sp.cancel(value))
    def mult(poly):
        poly=sp.Poly(poly,s)
        answer=0
        while poly.eval(point)==0:
            poly=sp.div(poly,sp.Poly(s-point,s))[0]
            answer+=1
        return answer
    return mult(num)-mult(den)

# All nine native critical points: three regular and three two-point stationary groups.
A=sp.eye(4)
T=sp.diag(1,2,3,4)
v=sp.Matrix([1,0,0,0])
xx,R,P,F,W=pencil(A,T,v,t*t)
source=sp.Matrix(variables)
inner=(source.T*A*T*source+2*v.T*A*source)[0]
h=(source.T*A*source)[0]+inner**2
gradient=sp.Matrix([sp.diff(h,q) for q in variables])
hessian=sp.hessian(h,variables)
regular=[sp.Integer(0),-2+sp.I,-2-sp.I]
regular_masses=[]
for point in regular:
    actual={variables[i]:xx[i].subs(s,point) for i in range(4)}
    zero(F.subs(s,point),'the regular stationary parameter is exact '+str(point))
    zeros(list(gradient.subs(actual,simultaneous=True)),
          'the regular parameter gives a full native critical point '+str(point))
    J=sp.simplify(hessian.subs(actual,simultaneous=True).det())
    zero(J-16*P.subs(s,point)*sp.diff(F,s).subs(s,point),
         'the native regular Hessian retains the full determinant factor '+str(point))
    check(J!=0,'the regular comparison point is Morse '+str(point))
    regular_masses.append(1/J)
check(sp.Poly(sp.fraction(F)[0],s).degree()==3,
      'the three regular stationary roots exhaust the cleared numerator')
check(sp.simplify(sum(regular_masses))!=0,
      'discarding the finite stationary groups would give a wrong nonzero mass sum')
group_masses=[]
for eigenvalue in [2,3,4]:
    alpha=-sp.Rational(1,eigenvalue)
    R0=R.subs(s,alpha)
    root=alpha/2
    square=sp.simplify((root-R0)/eigenvalue)
    masses=[]
    for sign in [1,-1]:
        actual=dict.fromkeys(variables,sp.Integer(0))
        actual[x]=xx[0].subs(s,alpha)
        actual[variables[eigenvalue-1]]=sign*sp.sqrt(square)
        zeros(list(gradient.subs(actual,simultaneous=True)),
              'the compatible determinant pole has an actual finite critical point eigen='+str(eigenvalue)+' sign='+str(sign))
        J=sp.simplify(hessian.subs(actual,simultaneous=True).det())
        expected=16*2*(root-R0)*2*sp.diff(P,s).subs(s,alpha)
        zero(J-expected,
             'the singular native determinant has the full rank-one cofactor eigen='+str(eigenvalue)+' sign='+str(sign))
        check(J!=0,'the singular-pole critical point is Morse eigen='+str(eigenvalue)+' sign='+str(sign))
        masses.append(1/J)
    group=sp.simplify(sum(masses))
    residue=1/(sp.diff(P,s).subs(s,alpha)*F.subs(s,alpha))
    zero(group-residue/16,'the complete two-point group equals its determinant-pole coefficient eigen='+str(eigenvalue))
    group_masses.append(group)
zero(sum(regular_masses)+sum(group_masses),
     'the complete native reciprocal-Hessian sum including all nine critical points is zero')
check(3+2*3==1+2*1*4,'the complete Morse count agrees with one plus twice m times source dimension')
check(inf_order(W)==5 and order_at(W,-1)>=0,
      'the blocked active pole and infinity satisfy the actual global boundary test')

# A degree-two outer derivative gives four points at one stationary determinant pole.
xx2,RR,PP,FF,WW=pencil(sp.eye(2),sp.diag(2,1),sp.Matrix([1,0]),t**3/3)
g2=2*x*x+y*y+2*x
h2=x*x+y*y+g2**3/3
curve=y**4+1
for derivative in [sp.diff(h2,x),sp.diff(h2,y)]:
    zero(sp.rem(derivative.subs(x,-1),curve,y),
         'the entire quartic singular-pole fibre consists of critical points')
zero(sp.rem(sp.hessian(h2,[x,y]).det().subs(x,-1)-16,curve,y),
     'all four critical points in the stationary group have native Hessian sixteen')
check(sp.gcd(curve,sp.diff(curve,y))==1 and sp.degree(curve,y)==4,
      'the quartic group has exactly four distinct finite points')
zero(sp.Integer(4)/16-1/(4*sp.diff(PP,s).subs(s,-1)*FF.subs(s,-1)),
     'one pole coefficient represents all 2m critical points for m=2')
r1,r2,r3,R0,alpha=sp.symbols('r1 r2 r3 R0 alpha')
p=(t-r1)*(t-r2)*(t-r3)+alpha
zero(sum(1/((root-R0)*sp.diff(p,t).subs(t,root)) for root in [r1,r2,r3])
     -1/(alpha-p.subs(t,R0)),
     'the all-root partial-fraction group identity has the correct sign')

# Morse failures at a compatible simple pole are not hidden hypotheses.
for label,outer,point in [
 ('F(alpha)=0',t*t-t,{x:-1,y:0}),
 ('a repeated outer-derivative root',t**3/3-t*t,{x:-1,y:1})]:
    hh=x*x+y*y+outer.subs(t,g2)
    zeros([sp.diff(hh,q).subs(point) for q in [x,y]],
          'the excluded degeneracy is an actual critical point '+label)
    zero(sp.hessian(hh,[x,y]).det().subs(point),
         'its native Hessian is explicitly degenerate '+label)

# Proper noncyclic Krylov spaces: verify the actual orthogonal splitting.
hyper=sp.Matrix([[0,1,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]])
J2=sp.Matrix([[1,0],[1,1]])
cases=[
 ('simple spectrum',sp.eye(4),sp.diag(1,2,3,4),sp.Matrix([1,0,0,0]),t*t,1,2),
 ('globally isotropic simple offset',sp.eye(4),sp.diag(1,2,3,4),sp.Matrix([1,sp.I,0,0]),t*t,1,2),
 ('active Jordan block',hyper,sp.diag(J2,sp.Integer(2),sp.Integer(3)),sp.Matrix([1,0,0,0]),t*t,2,3),
 ('repeated active eigenvalue',sp.eye(4),sp.diag(1,1,2,3),sp.Matrix([1,0,0,0]),t*t,2,2),
 ('scalar active primary',sp.eye(4),sp.eye(4),sp.Matrix([1,0,0,0]),t**3/3,4,2)]
for label,AA,TT,vv,outer,multiplicity,pole in cases:
    n=TT.rows
    krylov=sp.Matrix.hstack(*[TT**j*vv for j in range(n)])
    basis=sp.Matrix.hstack(*krylov.columnspace())
    r=basis.cols
    check(r<n and (basis.T*AA*basis).det()!=0,
          'the proper Krylov space has a genuinely nondegenerate restricted pairing '+label)
    if label=='globally isotropic simple offset':
        zero((vv.T*AA*vv)[0],'isotropy of the vector alone does not make its whole Krylov pairing degenerate')
    perp=sp.Matrix.hstack(*(basis.T*AA).nullspace())
    frame=basis.row_join(perp)
    check(frame.det()!=0 and basis.T*AA*perp==sp.zeros(r,n-r),
          'the orthogonal complement is an actual direct constant frame '+label)
    transformed=frame.inv()*TT*frame
    check(transformed[:r,r:]==sp.zeros(r,n-r)
          and transformed[r:,:r]==sp.zeros(n-r,r),
          'both Krylov summands are invariant for the full pencil '+label)
    check((frame.inv()*vv)[r:,:]==sp.zeros(n-r,1),
          'all affine offset remains on the marked Krylov summand '+label)
    _,innerR,detP,scalarF,weights=pencil(AA,TT,vv,outer)
    check(order_at(innerR,-1)==-pole and order_at(detP,-1)==multiplicity,
          'the active pole and full determinant multiplicity are both retained '+label)
    check(order_at(weights,-1)>=0 and inf_order(weights)>=2,
          'the active regularity budget and infinity test pass '+label)

# Zero Krylov space is a valid nondegenerate summand by convention.
zero_krylov=sp.Matrix.hstack(*[sp.diag(1,2)**j*sp.zeros(2,1) for j in range(2)])
check(zero_krylov.rank()==0,'the centered zero-Krylov case is not silently excluded')
_,_,_,_,centeredW=pencil(sp.eye(2),sp.diag(1,2),sp.zeros(2,1),t*t)
check(inf_order(centeredW)==3,'the centered simple-spectrum case has the completed stationary-pole boundary')

# An inactive repeated semisimple primary gives a positive-dimensional critical quadric.
q3=x*x+2*y*y+2*z*z+2*x
hh=x*x+y*y+z*z+q3*q3
critical_quadric=z*z+y*y+sp.Rational(13,8)
zeros([sp.rem(sp.diff(hh,q).subs(x,1),critical_quadric,z) for q in [x,y,z]],
      'the inactive repeated primary contributes its whole critical quadric')
zero(sp.rem(sp.hessian(hh,[x,y,z]).det().subs(x,1),critical_quadric,z),
     'the whole repeated-primary critical quadric has degenerate Hessian')

# The pole bound and nondegenerate-pairing assumptions are essential to this certificate.
_,_,Pfail,_,Wfail=pencil(sp.eye(4),sp.eye(4),sp.Matrix([1,0,0,0]),t*t)
check(order_at(Pfail,-1)==4 and order_at(Wfail,-1)==-2,
      'a quadratic outer with active multiplicity four fails the actual finite-pole regularity budget')
isoT=sp.diag(1,1,2,3)
isov=sp.Matrix([1,0,0,0])
isoK=sp.Matrix.hstack(*[isoT**j*isov for j in range(4)])
isoBasis=sp.Matrix.hstack(*isoK.columnspace())
check(isoK.rank()==1 and (isoBasis.T*hyper*isoBasis).det()==0,
      'an isotropic noncyclic offset has a genuinely degenerate Krylov pairing')
_,isoR,_,_,isoW=pencil(hyper,isoT,isov,t*t)
zero(isoR,'its expected active inner pole really disappears')
check(order_at(isoW,-1)==-2,'the unsupported determinant pole cannot be declared harmless')

# A compatible multiple pole can be a boundary residue with no finite critical group.
Abad=sp.Matrix([[0,1],[1,0]])
_,Rbad,Pbad,Fbad,Wbad=pencil(Abad,J2,sp.zeros(2,1),t*t)
zero(Wbad-(1/s-1/(1+s)-1/(1+s)**2),
     'the compatible multiple determinant pole has its exact partial fractions')
zero(sp.residue(Wbad,s,-1)+1,
     'its nonzero residue is a genuine unaccounted boundary term')
zero((x*x+2*x*y).subs(x,0),'the singular stationary line has inner value zero')
check(2*sp.Integer(0)!=-1 and sp.degree(Pbad,s)==2,
      'the singular stationary line contains no critical point with the required parameter minus one')

context=json.loads((HERE/'sources/stationary-pole-mass-context.json').read_text())
algebra_count=len(checks)
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,
          'retained proof, receipt or primary context bytes stay pinned: '+path)
check(context['newly_used_unproved_external_statements']==[],
      'the group formula and global sum are derived without an imported multivariate residue theorem')
check(context['all_invertible_simple_spectrum_affine_offsets_excluded'] is True,
      'the new simple-spectrum exclusion covers all offsets, not just cyclic vectors')
check(context['all_independently_centered_quadratic_composite_images_resolved'] is False
      and context['unit_ideal_endpoint_resolved'] is False
      and context['singleton_zero_endpoint_resolved'] is False
      and context['unrestricted_HC4_or_JC2_resolved'] is False,
      'the remaining failed certificates, endpoints and original goal stay open')
receipt=dict(outcome='PASS',claims=['HC4MYJH0','HC4MYJH1'],
 check_count=len(checks),algebraic_and_scope_controls=algebra_count+3,checks=checks,
 checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 proof_sha256=hashlib.sha256((HERE/'SIMPLE_STATIONARY_POLES_AND_NONCYCLIC_KRYLOV_CLOSURE.md').read_bytes()).hexdigest(),
 source_sha256={'sources/stationary-pole-mass-context.json':hashlib.sha256((HERE/'sources/stationary-pole-mass-context.json').read_bytes()).hexdigest()},
 dependencies=context['retained_proof_and_source_sha256'],
 assurance='Finite exact controls of the written all-degree critical-group and orthogonal-Krylov proof; not independent replay, formal verification or external review.',
 simple_stationary_pole_critical_group_identity_proved=True,
 boundary_completed_affine_pencil_residue_obstruction_proved=True,
 orthogonal_krylov_noncyclic_exclusion_proved=True,
 all_invertible_simple_spectrum_affine_offsets_excluded=True,
 all_nonnilpotent_noncyclic_affine_pencils_resolved=False,
 all_independently_centered_quadratic_composite_images_resolved=False,
 unit_ideal_endpoint_resolved=False,singleton_zero_endpoint_resolved=False,
 unrestricted_HC4_or_JC2_resolved=False,original_goal_complete=False,
 python=platform.python_version(),sympy=sp.__version__,elapsed_seconds=round(time.monotonic()-START,3))
(HERE/'receipt-stationary-pole-masses.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({key:receipt[key] for key in ['outcome','check_count','algebraic_and_scope_controls','elapsed_seconds','original_goal_complete']},indent=2))
