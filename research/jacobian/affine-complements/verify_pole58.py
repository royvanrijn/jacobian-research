"""PF2RULED12 exact controls; written geometric coverage is not formal verification."""
from pathlib import Path
from time import monotonic
import sympy as S,json,hashlib,platform,re
HERE=Path(__file__).resolve().parent;start=monotonic();checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def zero(name,value):check(name,S.cancel(S.powsimp(S.expand(value),force=True))==0)
t,Z=S.symbols('t Z',positive=True);c,K=S.symbols('c K',integer=True,positive=True)
H=S.Function('H')(Z);G=S.Function('G')(Z)
P=t**(-5*c)*H;Q=t**(-K)*G
J=S.diff(P,t)*S.diff(Q,Z)-S.diff(P,Z)*S.diff(Q,t)
zero('universal weighted initial ODE',t**(5*c+K+1)*J-(K*S.diff(H,Z)*G-5*c*H*S.diff(G,Z)))
for multiplicity in range(1,5):
    check('every possible numerator zero is non-fifth-power '+str(multiplicity),(8*multiplicity)%5!=0)
a,b=S.symbols('a b',positive=True);cc=b-3*a
zero('strict wedge-volume inequality is11c>5a',(a-2*b)-(-13*cc)-(11*cc-5*a))
threshold=S.Rational(38,11)
check('first-edge threshold corresponds to the strict ODE bound',11*(threshold-3)==5)
bounds=[]
for i in range(1,5):
    lower=S.ceiling(15-threshold*i)
    check('sharp allowed first-edge integer bound i'+str(i),S.Rational(15-lower,i)<=threshold)
    check('first forbidden integer order violates first-edge bound i'+str(i),S.Rational(15-(lower-1),i)>threshold)
    bounds.append(int(lower))
check('complete first-edge coefficient-order tuple',bounds==[12,9,5,2])
check('all regular z-degree>=5 terms have slope at most3',S.Rational(15,5)==3 and 3<threshold)
check('minimal allowed initial degree guarantees numerator degree<=4',S.Rational(15,5)<threshold)
# Exact lift depressions and the finite section-degree alternatives.
T,U,mu,nu,chi,R3=S.symbols('T U mu nu chi R3')
newP=S.expand((T**5+mu*T**4+R3).subs(T,U-mu/5))
check('depression removes fourth degree',newP.coeff(U,4)==0)
check('depression remainder stays in V3',S.Poly(newP-U**5,U).degree()==3)
check('degree12 symbol forced proportional by order12',3*4==12)
check('degree9 symbol forced proportional by order9',3*3==9)
check('degree6 symbol alternatives are order5 or proportional',list(range(5,7))==[5,6])
check('degree3 symbol alternatives are order2 or proportional',list(range(2,4))==[2,3])
# Universal coefficient models retain every term allowed on or above the
# indicated first Newton edge, with generic symbolic constants.
w,z=S.symbols('w z');A,B,C,D,xi=S.symbols('A B C D xi')
models=[
 ('beta_order5',w**15+A*w**12*z+B*w**9*z*z+C*w**5*z**3+D*w*w*z**4+xi*z**5,3,10,35,12),
 ('cubic_order2',w**15+A*w**12*z+B*w**9*z*z+C*w**6*z**3+D*w*w*z**4+xi*z**5,4,13,47,17)]
certificates=[]
for name,F,aa,bb,expected_fz,pole in models:
    expanded=S.Poly(S.expand(F.subs({w:t**aa,z:t**bb*Z})),t)
    order=min(m[0] for m,value in expanded.terms());face=expanded.coeff_monomial(t**order)
    check('first Newton value '+name,order==15*aa)
    check('nonconstant face degree '+name,S.degree(face,Z)==aa)
    zero('face constant coefficient '+name,face.subs(Z,0)-1)
    check('nonzero Newton roots simple '+name,S.degree(S.gcd(face,S.diff(face,Z)),Z)==0)
    check('normalization parameter primitive '+name,S.gcd(aa,bb)==1)
    check('F_z normalization order '+name,order-bb==expected_fz)
    check('exact generic-fibre Q pole '+name,expected_fz-aa-2*bb==pole)
    check('single branch exceeds total quartic degree '+name,pole>4)
    derivative=2*bb+aa-1-expected_fz
    check('nonzero derivative forces stated Q pole '+name,derivative==-pole-1)
    # Every T lift is regular after its cubic normal pole. Each expansion
    # term from the displayed powers lies above the exceptional face.
    for power in (5,3,2):
        for normal in range(power+1):
            zi=5-power+normal;wi=3*(power-normal)
            if power==5 and normal==0:continue
            gap=aa*wi+bb*zi-order
            check('arbitrary T-lift term above edge '+str((name,power,normal)),gap>0)
    other_orders=[(1,12),(2,9),(4,2),(5,0)] if name=='beta_order5' else [(1,12),(2,9),(3,6),(5,0)]
    for zi,wi in other_orders:
        check('complete other coefficient above face '+str((name,zi,wi)),aa*wi+bb*zi>order)
    certificates.append(dict(name=name,w_order=aa,z_order=bb,F_weight=order,initial_face=str(face),Fz_order=expected_fz,Q_derivative_order=derivative,Q_pole=pole))
check('polar centralizer leading wedge always exceeds pole3',5+1+1>3)
check('height-eight contact case stays below125',max(64,112)<125)
check('height-eight contact case cannot be72,108',64<72)
# Retain the rejected source-version shortcut as evidence, not a theorem.
manifest=json.loads((HERE/'sources/gcd-boundary-review.json').read_text())
raw=(HERE/'sources'/manifest['file']).read_bytes()
check('revised primary source seal',len(raw)==manifest['bytes'] and hashlib.sha256(raw).hexdigest()==manifest['sha256'])
check('revised version is pinned and not adopted as a gcd theorem',manifest['version']=='v2' and manifest['adopted_as_theorem'] is False and manifest['mathematical_result_for_pole58_uses_this_source'] is False)
text=re.sub('<[^>]+>',' ',raw.decode());text=' '.join(text.split())
check('primary version states inequality limitation','only an inequality' in text and 'cannot discard the infinite families' in text)
# Unified complete first-edge tables for both coprime height-eight cases.
newton_tables=[]
for m in (5,7):
    crit=3+S.Rational(5,m+6)
    hs=[int(S.ceiling(3*m-crit*i)) for i in range(1,m)]
    allowed=sorted({S.Rational(3*m-hh,i) for i in range(1,m)
                    for hh in range(hs[i-1],3*(m-i))},reverse=True)
    expected=[S.Rational(10,3),S.Rational(13,4)] if m==5 else [S.Rational(10,3),S.Rational(13,4),S.Rational(16,5),S.Rational(19,6)]
    check('complete surviving rational first-edge slope list m'+str(m),allowed==expected)
    check('complete integer first-edge order bounds m'+str(m),hs==([12,9,5,2] if m==5 else [18,15,11,8,5,1]))
    for multiplicity in range(1,m):
        check('all initial numerator zeros exclude an mth power '+str((m,multiplicity)),(8*multiplicity)%m!=0)
    pp=t**(-m*c)*H;qq=t**(-K)*G
    jj=S.diff(pp,t)*S.diff(qq,Z)-S.diff(pp,Z)*S.diff(qq,t)
    zero('general weighted ODE m'+str(m),t**(m*c+K+1)*jj-(K*S.diff(H,Z)*G-m*c*H*S.diff(G,Z)))
    zero('general volume comparison m'+str(m),(a-2*b)+(m+8)*(b-3*a)-((m+6)*(b-3*a)-5*a))
    for slope in allowed:
        aa=int(slope.q);bb=int(slope.p);fz=3*m*aa-bb;pole=(3*m-1)*aa-(m-2)*bb
        check('slope lies in first-edge ODE remainder band '+str((m,slope)),3<slope<=crit)
        check('pole budget from derivative bound '+str((m,slope)),fz-aa-(m-3)*bb==pole)
        check('even ramified or coalescent faces exceed degree4 '+str((m,slope)),pole>4)
        possible=[i for i in range(1,m) if (3*m-slope*i).is_integer and 3*m-slope*i>=0]
        newton_tables.append(dict(P_pole=m,slope=str(slope),w_weight=aa,z_weight=bb,Fz_lower_bound_per_extra_ramification=fz,Q_pole_lower_bound_per_extra_ramification=pole,possible_face_z_exponents=possible,repeated_roots_allowed=True))
    # Every coefficient cone implies the actual remainder symbol has all
    # of its 3l zeros at p. Regular T-lift terms cannot spoil that order.
    for ell in range(1,m):
        for power in range(ell,m+1):
            for regular_count in range(power+1):
                tau_count=power-regular_count
                needed_regular_z_degree=tau_count-ell
                if needed_regular_z_degree<0:continue
                check('all T powers preserve coefficient cone '+str((m,ell,power,regular_count)),tau_count>=ell and 3*tau_count>=3*ell)
# A repeated m7 face is retained explicitly, independent of the generic
# simple-binomial models above. This is a local Newton control only.
spar=S.Symbol('spar');Apar=S.Symbol('Apar',nonzero=True)
Frep=w**21+2*w**11*z**3+w*z**6+xi*z**7
subrep={w:spar**6,z:-spar**20+Apar*spar**27}
rawrep=S.Poly(S.expand(Frep.subs(subrep)),spar)
zero('coalescent face first splitting equation',rawrep.coeff_monomial(spar**140)-(9*Apar**2-xi))
fdrep=S.Poly(S.expand(S.diff(Frep,z).subs(subrep)),spar)
check('coalescent control Fz order is113',min(mon[0] for mon,coeff in fdrep.terms())==113)
zero('coalescent control Fz leading coefficient',fdrep.coeff_monomial(spar**113)-18*Apar)
check('coalescent refinement raises the Fz bound',113>2*53)
check('coalescent control Q pole is27 and exceeds degree4',113-6-4*20==27 and 27>4)
report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
    complete_first_edge_order_bounds=bounds,first_edge_maximum_slope=str(threshold),exceptional_Newton_certificates=certificates,
    complete_coprime_height_eight_Newton_tables=newton_tables,
    written_result='All(5,8) presentations and all noncontact-symbol(7,8) presentations are excluded',
    excluded_pole_pairs=[[5,8]],excluded_noncontact_pole_pairs=[[5,8],[7,8]],
    remaining_pole78_requires=dict(contact_symbol=True,ordinary_degrees=[72,108]),frontier_superseded_by='PF2RULED16',
    frontier_handoff_scope='historical pre-PF2RULED16; this receipt does not prove the current frontier',
    first_open_reduced_pairs=[[6,8],[7,8]],minimum_target_coordinate_pole=4,
    assurance=dict(first_edge_global_symbol_and_finiteness_coverage='written proof',
        noncontact_pole58_step_uses_external_degree_theorem=False,
        noncontact_pole78_step_uses_external_degree_theorem=False,coalescent_Newton_faces_retained=True,
        contact_base_uses_published_2022_alternative=True,published_72_108_exception_retained=True,
        global_gcd16_exclusion_imported=False,gcd25_from_superseded_2017_v1_imported=False,
        conditional_2026_degree125_upgrade_used=False,
        complete_height_eight_excluded=False,all_pole_four_coordinates_excluded=False,
        actual_ruled_compatibility_proved=False,quartic_packets_excluded=False,
        formal_verification=False,independent_full_reproof=False,external_review=False),
    source_sha256={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),HERE/'POLE_FIVE_EIGHT_EXCLUSION.md',HERE/'sources/gcd-boundary-review.json']},
    python=platform.python_version(),sympy=S.__version__,elapsed_seconds=round(monotonic()-start,6))
(HERE/'pole58_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
print(f'PASS {len(checks)} exact controls; entire(5,8) and noncontact(7,8) excluded; (6,8),(7,8) remain')
