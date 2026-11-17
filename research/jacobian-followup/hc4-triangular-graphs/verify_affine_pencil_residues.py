"""Exact controls for affine-pencil partial fractions, cyclicity and nilpotent endpoint reduction."""
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
    reduced=sp.cancel(value)
    return reduced==0 or sp.simplify(reduced)==0
def zero(value,label):check(iszero(value),label)
def zeros(values,label):check(all(iszero(v) for v in values),label)
s,t,c=sp.symbols('s t c')
def inf_order(value):
    num,den=sp.fraction(sp.cancel(value))
    return sp.degree(den,s)-sp.degree(num,s)
def order_at(value,point):
    value=sp.cancel(value)
    num,den=sp.fraction(value)
    def multiplicity(poly):
        poly=sp.Poly(poly,s)
        out=0
        while poly.eval(point)==0:
            poly=sp.div(poly,sp.Poly(s-point,s))[0]
            out+=1
        return out
    return multiplicity(num)-multiplicity(den)
def pencil(A,T,v,constant,outer):
    dim=T.rows
    xx=sp.simplify(-s*(sp.eye(dim)+s*T).inv()*v)
    B=A*T
    R=sp.cancel((xx.T*B*xx+2*v.T*A*xx)[0]+constant)
    P=sp.factor((sp.eye(dim)+s*T).det())
    F=sp.cancel(s-sp.diff(outer,t).subs(t,R))
    W=None if F==0 else sp.cancel(1/(P*F))
    return xx,R,P,F,W

# A shifted nonzero Jordan block: the full determinant and every critical mass can be checked.
A=sp.Matrix([[0,1],[1,0]])
T=sp.Matrix([[1,0],[1,1]])
v=sp.Matrix([1,0])
xx,R,P,F,W=pencil(A,T,v,0,t*t)
zero(R-sp.diff(-s*s*(v.T*A*(sp.eye(2)+s*T).inv()*v)[0],s),
     'the rational inner value is the derivative of the stationary pencil value')
zeros(list((sp.eye(2)+s*T)*xx+s*v),
      'the rational source point solves the complete affine stationarity equations')
check(sp.Matrix.hstack(v,T*v).det()!=0 and T.det()!=0,
      'the shifted Jordan block is genuinely cyclic and invertible')
check(order_at(R,-1)==-3 and order_at(P,-1)==2,
      'the nonzero length-two block has inner pole order three and determinant order two')
check(order_at(W,-1)>=0 and inf_order(W)>=2,
      'all finite pencil poles and infinity satisfy the actual residue gate')
M=sp.eye(2)-T
check(M.rank()<M.row_join(v).rank(),
      'the singular affine stationarity equation has no missing source solution')
x,y=sp.symbols('x y')
source=sp.Matrix([x,y])
B=A*T
inner=(source.T*B*source+2*v.T*A*source)[0]
h=(source.T*A*source)[0]+inner*inner
inner_gradient=sp.Matrix([sp.diff(inner,q) for q in [x,y]])
zz=inner_gradient.subs({x:xx[0],y:xx[1]},simultaneous=True)
N=A+s*B
zero(sp.diff(R,s)+(zz.T*N.inv()*zz)[0]/2,
     'the derivative of the inner value has the correct negative half-factor')
critical_matrix=2*N+2*zz*zz.T
zero(critical_matrix.det()-4*A.det()*P*sp.diff(F,s),
     'the rank-one critical Hessian determinant has its exact factor and sign')
roots=[sp.Integer(0),sp.Integer(1),-1+sp.sqrt(2),-1-sp.sqrt(2)]
masses=[]
native_j=[]
gradient=sp.Matrix([sp.diff(h,q) for q in [x,y]])
native_hessian=sp.hessian(h,[x,y])
for point in roots:
    zero(F.subs(s,point),'the scalar critical parameter is exact '+str(point))
    actual={x:xx[0].subs(s,point),y:xx[1].subs(s,point)}
    zeros(list(gradient.subs(actual,simultaneous=True)),
          'the scalar root gives a full native critical point '+str(point))
    J=sp.simplify(native_hessian.subs(actual,simultaneous=True).det())
    zero(J-4*A.det()*P.subs(s,point)*sp.diff(F,s).subs(s,point),
         'the literal native Hessian agrees at the critical parameter '+str(point))
    check(J!=0,'the comparison critical point is Morse '+str(point))
    masses.append(1/J);native_j.append(J)
zero(sum(masses),'all reciprocal native critical Hessians sum to zero')
check(native_j[0]==native_j[1] and native_j[2]!=native_j[0],
      'two matching critical Hessians do not replace the whole-spectrum test')
partial=sum(1/(P.subs(s,p)*sp.diff(F,s).subs(s,p)*(s-p)) for p in roots)
zero(W-partial,'the rational function is exactly the sum of its simple-pole partial fractions')
zero(sum(1/(P.subs(s,p)*sp.diff(F,s).subs(s,p)) for p in roots),
     'the coefficient of one over s vanishes exactly')
check(sp.Poly(sp.fraction(F)[0],s).degree()==4,
      'no regular critical root is lost by the scalar reduction')

# A nonempty open four-variable independent-center example, without diagonalization assumptions in the theorem.
A4=sp.eye(4)
T4=sp.diag(1,2,3,4)
v4=sp.ones(4,1)
K=sp.Matrix.hstack(*[T4**j*v4 for j in range(4)])
check(K.det()==12 and T4.det()==24,'the four-variable cyclic and invertible gates are a nonempty open condition')
for degree in [2,3]:
    _,R4,P4,F4,W4=pencil(A4,T4,v4,1,t**degree+t)
    check(inf_order(W4)==5,'the invertible cyclic example has no infinity residue outer='+str(degree))
    check(all(order_at(R4,-sp.Rational(1,j))==-2 for j in [1,2,3,4]),
          'each simple nonzero primary has its exact double inner pole outer='+str(degree))
    check(all(order_at(W4,-sp.Rational(1,j))>=0 for j in [1,2,3,4]),
          'all finite boundary residues vanish outer='+str(degree))
    check(all((sp.eye(4)-T4/j).rank()<(sp.eye(4)-T4/j).row_join(v4).rank()
              for j in [1,2,3,4]),
          'none of the singular pencil parameters has an affine stationary solution outer='+str(degree))

# Nondiagonalizable cyclic length four and nilpotent infinity growth.
AA=sp.zeros(4)
NN=sp.zeros(4)
for j in range(4):
    AA[j,3-j]=1
    if j<3:NN[j+1,j]=1
vv=sp.Matrix([1,0,0,0])
for label,TT in [('nonzero',sp.eye(4)+NN),('nilpotent',NN)]:
    check(AA*TT==(AA*TT).T,'the cyclic chain is A-self-adjoint '+label)
    check(sp.Matrix.hstack(*[TT**j*vv for j in range(4)]).det()!=0,
          'the linear offset is cyclic through the whole chain '+label)
    _,RR,PP,FF,WW=pencil(AA,TT,vv,0,t*t)
    if label=='nonzero':
        check(order_at(RR,-1)==-5 and order_at(PP,-1)==4,
              'the length-four nonzero primary has pole order five against determinant order four')
        check(order_at(WW,-1)>=0 and inf_order(WW)==5,
              'the nonsemisimple invertible chain passes every residue boundary')
    else:
        zero(RR-5*s**4,'the nilpotent cyclic inner value has exact infinity degree four')
        check(PP==1 and sp.degree(FF,s)==4 and inf_order(WW)==4,
              'the nilpotent chain is excluded by polynomial degree rather than longer-jet checking')

# The exceptional simple-zero-eigenvalue quadratic outer case must not be silently included.
TS=sp.diag(0,1,2,3)
_,RS,PS,FS,WS=pencil(A4,TS,v4,sp.Rational(11,6),-t*t/4)
check(sp.Matrix.hstack(*[TS**j*v4 for j in range(4)]).det()!=0,
      'the double-infinity-resonance control is still cyclic')
zero(sp.limit(RS/s,s,sp.oo)+2,'the zero primary contributes the correct nonzero linear infinity slope')
zero(sp.limit(RS+2*s,s,sp.oo),'the constant infinity term is exactly resonant')
check(inf_order(FS)==2 and inf_order(WS)==1,
      'double infinity resonance leaves a genuine one-over-s boundary term')
check(all(order_at(WS,-sp.Rational(1,j))>=0 for j in [1,2,3]),
      'that unresolved control passes the finite boundaries but fails infinity')
_,_,_,FH,WH=pencil(A4,TS,v4,0,t**3)
check(inf_order(WH)>=2,'a higher nonlinear outer derivative restores the cyclic infinity obstruction')

# Nilpotent affine pencils: empty, singleton, multiple, and positive-dimensional critical strata.
vlin=sp.Matrix([sp.Rational(1,2),0,0,0])
Tzero=sp.zeros(4)
_,R0,_,F0,W0=pencil(A4,Tzero,vlin,0,-t*t)
zero(F0,'a completely resonant polynomial parameter equation is identically zero')
x,y,z,w,X=sp.symbols('x y z w X')
coord_h=y*y+z*z+w*w-2*x-1
Qcoord=x*x+y*y+z*z+w*w-X*X
zero(Qcoord.subs(X,x+1)-coord_h,'the fixed nilpotent affine graph is an actual coordinate potential')
recovered_x=(y*y+z*z+w*w-sp.Symbol('V')-1)/2
zero(coord_h.subs(x,recovered_x)-sp.Symbol('V'),
     'the coordinate graph has literal polynomial source recovery')
_,_,_,Fc,_=pencil(A4,Tzero,vlin,1,-t*t)
zero(Fc-2,'its stationary polynomial has no root and exactly the empty critical stratum')
_,_,_,Fl,_=pencil(A4,Tzero,vlin,1,t*t)
check(sp.degree(Fl,s)==1,'the polynomial degree-one branch has exactly one stationary parameter')
_,_,_,Fm,Wm=pencil(A4,Tzero,vlin,1,t**3)
check(sp.degree(Fm,s)==2 and inf_order(Wm)==2,
      'a higher polynomial stationary equation meets the complete residue obstruction')
zero((x*x+y*y+z*z+w*w-x*x)-(y*y+z*z+w*w),
     'the identically-zero parameter equation corresponds to a full critical line')
Ap=sp.Matrix([[0,1,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]])
Bp=sp.diag(1,0,0,0)
Tp=Ap.inv()*Bp
positive=2*x*y+z*z+w*w+(x*x)**3
check(Tp!=sp.zeros(4) and Tp*Tp==sp.zeros(4),
      'the nonlinear singleton control has a genuinely nonzero nilpotent inner pencil')
zero(sp.hessian(positive,[x,y,z,w]).det()+16,
     'the nonlinear singleton graph has actual unit Hessian outside the infinity gate')
_,_,_,Fp,Wp=pencil(Ap,Tp,sp.zeros(4,1),0,t**3)
zero(Fp-s,'the centered nilpotent singleton retains its exact stationary root')
check(inf_order(Wp)==1,'its nonzero infinity coefficient explains why the zero-sum test does not apply')
# Missing critical points at a singular pencil are a separate forbidden shortcut.
_,_,Pbad,Fbad,Wbad=pencil(sp.eye(1),sp.eye(1),sp.zeros(1,1),0,t*t)
check((sp.eye(1)-sp.eye(1)).row_join(sp.zeros(1,1)).rank()==0,
      'a centered singular pencil may have stationary solutions and cannot be discarded')
check(order_at(Wbad,-1)==-1,'that discarded singular pencil would contribute an unaccounted finite pole')

context=json.loads((HERE/'sources/affine-pencil-residue-context.json').read_text())
algebra_count=len(checks)
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,
          'retained endpoint, transfer or previous proof bytes remain pinned: '+path)
for entry in context['newly_cached_primary_context']:
    check(hashlib.sha256((HERE/'sources'/entry['retained_file']).read_bytes()).hexdigest()==entry['sha256'],
          'primary residue context bytes are pinned')
check(context['newly_used_unproved_external_statements']==[],
      'the partial-fraction theorem is proved directly rather than imported with missing hypotheses')
check(context['unit_ideal_endpoint_resolved'] is False and context['singleton_zero_endpoint_resolved'] is False,
      'both actual endpoint problems remain unresolved')
check(context['unrestricted_HC4_or_JC2_resolved'] is False,
      'the original all-embedding objective remains open')
receipt=dict(outcome='PASS',claims=['HC4MYJF0','HC4MYJF1','HC4MYJF2'],
 check_count=len(checks),algebraic_and_scope_controls=algebra_count+3,checks=checks,
 checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 proof_sha256=hashlib.sha256((HERE/'AFFINE_QUADRATIC_PENCIL_RESIDUES_AND_ENDPOINTS.md').read_bytes()).hexdigest(),
 source_sha256={'sources/affine-pencil-residue-context.json':hashlib.sha256((HERE/'sources/affine-pencil-residue-context.json').read_bytes()).hexdigest()},
 dependencies=context['retained_proof_and_source_sha256'],
 assurance='Finite exact controls for a written all-degree pencil/partial-fraction proof; not independent replay, formal verification or external review.',
 rational_affine_pencil_residue_obstruction_proved=True,
 cyclic_affine_pencil_source_orbits_excluded_under_stated_conditions=True,
 nilpotent_affine_pencil_family_reduced_to_two_old_endpoints=True,
 cyclic_double_infinity_resonance_resolved=False,
 all_independently_centered_quadratic_composite_images_resolved=False,
 unit_ideal_endpoint_resolved=False,singleton_zero_endpoint_resolved=False,
 unrestricted_HC4_or_JC2_resolved=False,original_goal_complete=False,
 python=platform.python_version(),sympy=sp.__version__,elapsed_seconds=round(time.monotonic()-START,3))
(HERE/'receipt-affine-pencil-residues.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({key:receipt[key] for key in ['outcome','check_count','algebraic_and_scope_controls','elapsed_seconds','original_goal_complete']},indent=2))
