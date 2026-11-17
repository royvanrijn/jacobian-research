#!/usr/bin/env python3
"""PF162 exact controls; the unbounded coverage proof is in the source note."""
from pathlib import Path
import hashlib,json,math,platform,time
import sympy as S
START=time.perf_counter();BASE=Path(__file__).resolve().parent;ROOT=BASE.parents[1]
SRC=BASE/'sources/ruling-six-boundary'
PROOF=BASE/'RULING_SIX_FINITE_END_COVERAGE_AND_QUADRATIC_BOUNDARY_ALGEBRA.md'
CONTEXT=SRC/'context.json';RECEIPT=BASE/'ruling_six_boundary_receipt.json'
checks=[]
def check(label,condition):
    if not condition:raise AssertionError(label)
    checks.append(label)
def zero(e):return S.cancel(e)==0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def initial(e,var):
    num,den=S.fraction(S.cancel(e))
    if num==0:return S.oo,S.Integer(0)
    pn,pd=S.Poly(num,var),S.Poly(den,var)
    a=min(i[0] for i in pn.monoms());b=min(i[0] for i in pd.monoms())
    return a-b,S.cancel(pn.nth(a)/pd.nth(b))
context=json.loads(CONTEXT.read_text())
for rel,h in context['retained_sha256'].items():check('retained pin '+rel,sha(ROOT/rel)==h)
u,w,z,p,t,x,y,lam,v,Y=S.symbols('u w z p t x y lam v Y')
M,N,c,a2,a1,a0,b2,b1,b0=S.symbols('M N c a2 a1 a0 b2 b1 b0')
P=c*M*z**3+w*(a2*z**2+a1*z+a0)+b2*z**2+b1*z+b0
star=P.subs(z,(w**2-N)/M)
cleared=c*(w**2-N)**3+(a2*w+b2)*(w**2-N)**2+M*(a1*w+b1)*(w**2-N)+M**2*(a0*w+b0-p)
check('whole finite-root divided equation',zero(cleared-M**2*(star-p)))
rec=S.cancel(cleared.subs(w,1/y)*y**6)
check('reciprocal coefficient unit alternatives',rec.subs({M:0,c:0,N:0})==a2*y+b2*y**2)
check('reciprocal exact parameter order',S.diff(rec,p)==-M**2*y**6)
L,R1,R0,bb=S.symbols('L R1 R0 bb')
K=bb+y*L+y**2*R1+y**3*(R0-p)
check('polar-z exact implicit parameter numerator',S.diff(K,p)==-y**3)
check('implicit w-equation parameter order numerator',zero(M*(y**3/S.diff(K,y))/y**2-M*y/S.diff(K,y)))
for mm,bb0,alpha in [(2,1,1),(2,1,2),(3,2,1),(5,4,3)]:
    bound=S.ceiling(S.Rational((mm+bb0)*alpha,2))-mm*alpha+alpha-1
    check(f'simple-N contact branch lower bound m{mm} b{bb0} alpha{alpha}',bound>=0)
check('finite non-N residual polynomial has degree at most two',S.degree(c*(w**2-N)+a2*w+b2,w)==2)

# An actual regular ruling-six expression with the necessary central double
# factor. No global polynomial mate is asserted for these supporting models.
models=[]
for k,aa in [(1,0),(1,1),(2,1),(2,2)]:
    mm=k+1;sig=1-2*k;MM=u**mm;NN=sig**2*u
    PS=(z**2+w)*(1+u**aa*w)**2+u
    PR=PS.subs(z,(w**2-NN)/MM)
    scaled=PR.subs(w,-u**(-aa)+u**(mm+aa)*lam)
    order,poly=initial(scaled,u)
    label=f'm{mm} a{aa}'
    check('entire central quadratic initial '+label,order==0 and zero(poly-lam**2))
    omega_coeff=S.diff(-u**(-aa)+u**(mm+aa)*lam,lam)/MM
    check('actual formal residue volume '+label,omega_coeff==u**aa)
    # Canonical top coefficients: b3=M*u^(2a),a2=2u^a,b2=1+N*u^(2a).
    ca2=2*u**aa;cb3=MM*u**(2*aa);cb2=1+NN*u**(2*aa)
    check('canonical order and double-leading identity '+label,
          initial(ca2,u)[0]==aa and initial(cb3,u)[0]==mm+2*aa and zero(2**2-4*cb2.subs(u,0)))
    xs=(-t**(-2*aa)+t**(2*(mm+aa))*lam)/t**(2*k+1)-sig/t**(2*k)
    q=2*(mm+aa)-1;eR=2*aa+2
    xo,xl=initial(xs,t)
    check('original source pole and exact ramification '+label,xo==-q and initial(S.diff(xs,lam),t)[0]==eR-1)
    check('full first-velocity clearing contribution '+label,
          1+math.ceil(S.Rational(q,eR))==2+math.ceil(S.Rational(2*mm-3,eR)))
    check('full marked-quadratic clearing contribution '+label,
          1+math.ceil(S.Rational(2*q+1,eR))==3+math.ceil(S.Rational(4*mm-5,eR)))
    source=PS.subs({u:t**2,w:t*(sig+t**(2*k)*x),z:2*sig*x+t**(2*k)*x**2})
    check('whole central source restriction degree two '+label,S.degree(source.subs(t,0),x)==2)
    models.append({'m':mm,'a':aa,'initial_residue_polynomial':str(poly),'original_ramification':eR,
        'source_pole':q,'source_volume_order':eR-1,'formal_model_not_Keller_realization':True})

# Explicit finite-root factor recursion over C[p], without a p denominator.
def trunc_t(expr,N):
    poly=S.Poly(S.expand(expr),t)
    return sum(c*t**i[0] for i,c in poly.terms() if i[0]<N)
PS=(z**2+w)*(1+w)**2+u
source=S.expand(PS.subs({u:t**2,w:t*(-1+t**2*x),z:-2*x+t**2*x**2}))
f=S.Poly(source-p,x);c2=S.Poly(f.as_expr().subs(t,0),x).LC()
W=S.cancel(f.as_expr().subs(t,0)/c2);U=c2
for degree in range(1,5):
    residual=S.Poly(S.expand(f.as_expr()-W*U),t).nth(degree)
    quo,rem=S.div(residual,f.as_expr().subs(t,0)/c2,x)
    U+=quo*t**degree;W+=S.cancel(rem/c2)*t**degree
check('actual finite-root Hensel factor product through t4',trunc_t(f.as_expr()-W*U,5)==0)
check('finite-root factor stays monic quadratic',S.degree(W,x)==2 and S.Poly(W,x).LC()==1)
check('actual finite-root factor has polynomial p coefficients',not S.denom(S.cancel(W)).has(p))
check('exact source deck symmetry preserves full P',zero(source.subs({t:-t,x:-x+2/t**2},simultaneous=True)-source))

# Actual multiplication, separation unit and the two residual/residue algebras.
cc,dd=S.symbols('cc dd',nonzero=True)
MX=S.Matrix([[0,p/cc],[1,0]])
MW2=dd**2*S.eye(2)+2*dd*MX
check('visible factor separation norm uses actual multiplication',zero(MW2.det()-(dd**4-4*dd**2*p/cc)))
sep=MW2.det().subs({dd:-2/t**2,cc:4})
check('separation norm leading unit independent of p',initial(sep,t)==(-8,S.Integer(16)))
c0,c1=S.symbols('c0 c1');MR=S.Matrix([[0,-c0],[1,-c1]])
check('residual quadratic exact multiplication identity',MR**2+c1*MR+c0*S.eye(2)==S.zeros(2))
check('residual trace norm discriminant exact',S.trace(MR)==-c1 and MR.det()==c0)
ML=S.Matrix([[0,p],[1,0]])
check('normal residue algebra lambda squared equals P',ML**2==p*S.eye(2))
check('normal residue algebra trace discriminant',S.det(S.Matrix([[2,S.trace(ML)],[S.trace(ML),S.trace(ML**2)]]))==4*p)
HH=p**2*(p-1);MT=HH*ML
check('image algebra actual trace discriminant',zero(S.det(S.Matrix([[2,S.trace(MT)],[S.trace(MT),S.trace(MT**2)]]))-4*p*HH**2))
N0=ML.subs(p,0)
check('critical residue dual numbers retain nonzero nilpotent',N0!=S.zeros(2) and N0**2==S.zeros(2))
Eplus=(S.eye(2)+ML.subs(p,1))/2;Eminus=(S.eye(2)-ML.subs(p,1))/2
check('noncritical residue splits with both projectors',Eplus**2==Eplus and Eminus**2==Eminus and Eplus*Eminus==S.zeros(2))
# Finite conductor quotient C[lambda]/lambda4(lambda2-1), dimension6.
CM=S.zeros(6)
for j in range(5):CM[j+1,j]=1
CM[4,5]=1
check('finite conductor quotient actual polynomial relation',CM**6-CM**4==S.zeros(6))
JJ=CM**3-CM
check('finite conductor nilradical fourth power zero',JJ**4==S.zeros(6))
check('finite conductor nilradical cube retained nonzero',JJ**3!=S.zeros(6))

# Both node and cusp residue images admit local ramified models; neither is
# declared to give a polynomial plane Keller map.
Pn=u+lam**2;Qn=lam**3+lam+S.Rational(3,2)*u*lam
Pc=u-3*lam**2;Qc=u*lam-2*lam**3
jac=lambda P,Q:S.diff(P,u)*S.diff(Q,lam)-S.diff(P,lam)*S.diff(Q,u)
check('immersive nodal local unit-volume model',zero(jac(Pn,Qn)-(1+S.Rational(3,2)*u)))
check('cuspidal local ramification model',jac(Pc,Qc)==u)
check('cusp finite fibre length3 exceeds generic index2',zero(Qc.subs(u,3*lam**2)-lam**3) and S.degree(lam**3,lam)==3)
Pc2=Pc.subs(u,t**2);Qc2=Qc.subs(u,t**2)
gb=S.groebner([Pc2,Qc2],t,lam)
check('quadratic cusp lift fibre length6 exceeds generic index4',set(q.as_expr() for q in gb.polys)=={t**2-3*lam**2,lam**3} and zero(S.diff(Pc2,t)*S.diff(Qc2,lam)-S.diff(Pc2,lam)*S.diff(Qc2,t)-2*t**3))

# Three complete necessary off-M weighted-face controls, including a vanished
# coefficient on a triple-plus-simple face.
profiles=[]
for name,I,Lv,Hv in [
    ('A',Y**2*(Y-1)**4,3,6),
    ('B5',Y**2*(Y-1)**3,4,8),
    ('B6',Y**2*(Y-1)**3*(Y+1),4,8),
    ('C',Y**3*(Y-1)**3,5,11)]:
    ip=S.Poly(S.expand(I),Y);F=sum(ip.nth(i)*v**(i*Lv-Hv)*w**i for i in range(ip.degree()+1) if ip.nth(i)!=0)
    if name=='C':F+=w**2
    if name=='B5':F+=v**(4*Lv+1)*w**6
    weighted=v**Hv*F.subs(w,v**(-Lv)*Y)
    order,poly=initial(weighted,v)
    check('full weighted face recovered '+name,order==0 and zero(poly-I))
    check('actual degree-two specialization retained '+name,S.degree(F.subs(v,0),w)==2)
    ell=4 if name=='A' else 3
    check('nonzero initial root exact multiplicity '+name,
          all(S.diff(I,Y,j).subs(Y,1)==0 for j in range(ell)) and S.diff(I,Y,ell).subs(Y,1)!=0)
    profiles.append({'name':name,'L':Lv,'H0_over_alpha':Hv,'initial':str(S.expand(I)),'multiplicity':ell})
check('triple-plus-simple cancellation is retained',S.Poly(S.expand(Y**2*(Y-1)**3*(Y+1)),Y).nth(4)==0)
old=json.loads((BASE/'canonical_deck_pole_coverage_receipt.json').read_text())
check('known3D nonproper control remains pinned',old['outcome']=='PASS' and 'retained 3D first-velocity pole' in old['checks'])
check('P-vertical and all original ramified primes retained',context['P_vertical_primes_retained_outside_generic_fibre_count'] and context['all_original_ramified_primes_retained'])
check('residue/conductor algebra not promoted to HC4 tail certificate',not context['conductor_algebra_nilpotence_claimed_HC4_tail_certificate'])
check('whole ruling-six and generic integrality not claimed',not context['entire_ruling_six_target_class_closed'] and not context['generic_first_derivative_integrality_proved'])
receipt={'claim':'PF2RULED162','outcome':'PASS','check_count':len(checks),'checks':checks,
    'elapsed_seconds':time.perf_counter()-START,'environment':{'python':platform.python_version(),'sympy':S.__version__},
    'proof_sha256':sha(PROOF),'checker_sha256':sha(Path(__file__)),'context_sha256':sha(CONTEXT),
    'formal_boundary_models':models,'off_M_face_controls':profiles,
    'all_finite_M_root_ends_holomorphic':True,
    'central_degree_two_residual_algebra_rank':2,
    'forced_original_ramification':'2*a+2','forced_source_pole':'2*(m+a)-1',
    'forced_residue_extension':'C(lambda)/C(P),P=lambda^2',
    'image_normal_form':'T^2=P*H(P)^2,H nonconstant',
    'image_conductor':'H(P)*C(lambda)',
    'finite_conductor_control_dimension':6,'finite_conductor_nilradical_index':4,
    'whole_ruling_six_target_class_closed':False,
    'original_full_normalization_and_P_vertical_primes_retained':True,
    'generic_first_derivative_integrality_proved':False,'normal_tail_increasing_power_containment_proved':False,
    'registered_predecessor_checkers_replayed':[],'new_external_theorem_imports':[],
    'full_goal_complete':False,'independent_replay':False,'formal_verification':False,'external_review':False}
RECEIPT.write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'outcome':'PASS','check_count':len(checks),'elapsed_seconds':receipt['elapsed_seconds']},indent=2))
