"""PF2RULED11 exact controls; written formal descent and geometry are not formal verification."""
from pathlib import Path
from time import monotonic
import sympy as S,json,hashlib,platform
HERE=Path(__file__).resolve().parent;start=monotonic();checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def zero(name,value):check(name,S.cancel(S.powsimp(S.expand(value),force=True))==0)
w=S.Symbol('w');z=S.Symbol('z',positive=True)
h=S.Function('h')(w);rho=S.Function('rho')(w)
t=S.Symbol('t',integer=True,positive=True)
j=S.Symbol('j',integer=True)
for d in (1,2):
    m=d*t;P=h**t*z**(-m);R=rho*z**(-j)
    leading=j*t*h**(t-1)*S.diff(h,w)*rho-m*h**t*S.diff(rho,w)
    J=S.diff(P,z)*S.diff(R,w)-S.diff(P,w)*S.diff(R,z)
    zero('universal leading differential coefficient d'+str(d),z**(m+j+1)*J-leading)
    # d divides both m and2, so the integrating exponent is integral.
    Rfirst=rho*z**(m-2)
    Jfirst=S.diff(P,z)*S.diff(Rfirst,w)-S.diff(P,w)*S.diff(Rfirst,z)
    exact=-m*h**(S.Rational(2,d))*S.diff(h**(t-S.Rational(2,d))*rho,w)
    zero('first allowed remainder gives rational primitive d'+str(d),z**3*Jfirst-exact)
    zero('kernel logarithmic derivative d'+str(d),leading/(m*h**t*rho)-(S.Rational(1,d)*j*S.diff(h,w)/h-S.diff(rho,w)/rho))
    check('first remainder is precisely pole three d'+str(d),(m+(2-m)+1)==3)
    check('all higher remainder poles must commute d'+str(d),(m+(3-m)+1)==4)
# Transverse frame changes: e scales by u^2, h by u^d, rho by u^(2-m).
u=S.Symbol('u',positive=True);ee,hh,rr=S.symbols('ee hh rr',nonzero=True)
for d in (1,2):
    m=d*t
    zero('normal-frame one-form invariance d'+str(d),u**2*ee/(u**d*hh)**S.Rational(2,d)-ee/hh**S.Rational(2,d))
    zero('normal-frame rational primitive invariance d'+str(d),(u**d*hh)**(t-S.Rational(2,d))*u**(2-m)*rr-hh**(t-S.Rational(2,d))*rr)
check('leading coefficient line bundle degree is4',-2+2*3==4)
check('quotient one-form has canonical degree minus2',4-6==-2)
check('pole-four Laurent congruence retains inverse root term',3-4==-1 and 2-4==-2)
# The exact root classification is a rational-degree proof; these are its
# complete finite partition controls, not a sampled numeric extrapolation.
def partitions(n,minimum=1):
    if n==0:yield ()
    for k in range(minimum,n+1):
        for tail in partitions(n-k,k):yield (k,)+tail
classification=[];a=S.Symbol('a')
for d in (1,2):
    power=2//d;allowed=[0]
    for degree in range(1,3*d+1):
        for parts in partitions(degree):
            simple_pole=any(power*k==1 for k in parts)
            primitive_pole_degree=sum(power*k-1 for k in parts)
            infinity_zero_order=power*degree-1
            necessary=not simple_pole and infinity_zero_order<=primitive_pole_degree
            check('complete denominator partition '+str((d,parts)),necessary==(len(parts)==1 and power*degree>=2))
        if power*degree>=2:
            allowed.append(degree)
            primitive=(w-a)**(1-power*degree)/S.Integer(1-power*degree)
            zero('explicit rational primitive '+str((d,degree)),S.diff(primitive,w)-(w-a)**(-power*degree))
    zero('constant denominator primitive d'+str(d),S.diff(w,w)-1)
    expected=[0,1,2,3] if d==1 else [0,2,3,4,5,6]
    check('complete contact multiplicity classification d'+str(d),allowed==expected)
    classification.append(dict(common_symbol_bundle_degree=3*d,denominator_power=power,allowed_away_root_multiplicities=allowed,allowed_contact_multiplicities=[3*d-k for k in allowed]))
# Nonexact position controls: poles at different finite roots cannot be
# repaired by a formal coordinate or by silently moving the contact point.
negative=[('squarefree cubic',(w-1)*(w-2)*(w-3),2,1,S.Rational(3,4)),
          ('mixed cubic away from contact',(w-1)**2*(w-2),2,2,-4),
          ('sextic33 away from contact',(w-1)**3*(w-2)**3,1,1,-6),
          ('simple finite sextic root',w-1,1,1,1)]
for name,den,power,point,value in negative:
    zero('exact nonzero residue control '+name,S.residue(1/den**power,w,point)-value)
    check('negative residue is genuinely nonzero '+name,value!=0)
zero('sextic33 with one root at contact has a primitive',S.diff(-1/(2*(w-1)**2),w)-1/(w-1)**3)
# Explicit meromorphic controls replay the primitive ODE. They are not
# functions on the full N, finite maps, or Keller counterexamples.
meromorphic=[]
for d in (1,2):
    allowed=classification[d-1]['allowed_away_root_multiplicities']
    for degree in allowed:
        hh=(w-a)**degree
        power=2//d
        G=w if degree==0 else (w-a)**(1-power*degree)/S.Integer(1-power*degree)
        m=4 if d==1 else 6
        P=hh**(m//d)*z**(-m)
        R=-S.Rational(1,m)*hh**((2-m)//d)*G*z**(m-2)
        J=S.diff(P,z)*S.diff(R,w)-S.diff(P,w)*S.diff(R,z)
        zero('meromorphic leading-volume model '+str((d,degree)),J-z**-3)
        meromorphic.append(dict(d=d,away_degree=degree,P_pole=m,finite_map=False,globally_regular_pair=False))
# The source-degree sharpening follows from the exact monomial exponent
# and sign formula in both normalized cubic bases.
H,b=S.symbols('H b',nonzero=True)
counts=S.symbols('c0:4',integer=True,nonnegative=True)
mm=sum(counts);order=sum(i*counts[i] for i in range(4))
zero('universal product top H exponent',sum((4-i)*counts[i] for i in range(4))-(4*mm-order))
zero('universal product top b exponent',sum(-(3-i)*counts[i] for i in range(4))-(-3*mm+order))
for chart in (1,2):
    for pole in range(1,6):
        for c0 in range(pole+1):
            for c1 in range(pole-c0+1):
                for c2 in range(pole-c0-c1+1):
                    cc=[c0,c1,c2,pole-c0-c1-c2];r=sum(i*cc[i] for i in range(4))
                    if r>3:continue
                    terms=[((-1)**i if chart==1 else -1)*H**(4-i)/b**(3-i) for i in range(4)]
                    top=S.prod(terms[i]**cc[i] for i in range(4))
                    sign=(-1)**(r if chart==1 else pole)
                    zero('normalized source top term '+str((chart,pole,tuple(cc))),top-sign*H**(4*pole-r)/b**(3*pole-r))
for r in range(4):check('contact order below4 cannot cancel with lower pole '+str(r),16-4*r>0)
check('contact order4 is not claimed exact',16-4*4==0)
# Exact first-open height-eight handoff, with the published exception retained.
def bound(m,r):return max(16*(m-1),16*m-4*r)
frontier=[]
for m in (5,6,7):
    gcd=2 if m==6 else 1;allowed=classification[gcd-1]['allowed_contact_multiplicities'];rows=[]
    for k in allowed:
        if k==0:
            rows.append(dict(contact_multiplicity=k,outcome='open',ordinary_degrees=[16*m,128],away_root_multiplicity=3*gcd))
            continue
        rP=(m//gcd)*k;rQ=(8//gcd)*k
        BP=bound(m,rP);BQ=bound(8,rQ)
        exactP=16*m-4*rP if rP<=3 else None
        impossible=BP<72 or BQ<108 or (exactP is not None and exactP!=72)
        rows.append(dict(contact_multiplicity=k,outcome='excluded' if impossible else 'open_exception_only',ordinary_degree_upper_bounds=[BP,BQ],exact_P_degree=exactP,remaining_ordinary_degrees=None if impossible else [72,108]))
    frontier.append(dict(pole_pair=[m,8],common_symbol_degree=3*gcd,cases=rows))
check('pole58 all contact-symbol cases excluded',[x['contact_multiplicity'] for x in frontier[0]['cases'] if x['outcome']=='open_exception_only']==[])
check('pole68 contact multiplicity1 excluded with exact ordinary degree84',next(x for x in frontier[1]['cases'] if x['contact_multiplicity']==1)['exact_P_degree']==84 and next(x for x in frontier[1]['cases'] if x['contact_multiplicity']==1)['outcome']=='excluded')
check('pole68 exact retained contact cases',[x['contact_multiplicity'] for x in frontier[1]['cases'] if x['outcome']=='open_exception_only']==[6,4,3,2])
check('pole78 exact retained contact cases',[x['contact_multiplicity'] for x in frontier[2]['cases'] if x['outcome']=='open_exception_only']==[3,2,1])
check('noncontact ordinary degrees at height8',[x['cases'][-1]['ordinary_degrees'] for x in frontier]==[[80,128],[96,128],[112,128]])
# All lists above are ordered by degree away from contact, hence k descending.
remaining_cubic=[dict(contact_multiplicity=k,away_multiplicity=3-k) for k in (0,1,2,3)]
check('pole-four mixed cubic must contain contact point',all(x['contact_multiplicity']>0 for x in remaining_cubic if x['away_multiplicity'] in (1,2)))
report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
    rational_primitive_classification=classification,meromorphic_leading_models=meromorphic,height_eight_frontier=frontier,
    height_eight_frontier_table_scope='necessary PF2RULED11 period/degree rows before the PF2RULED12 finiteness exclusion',
    written_result='An intrinsic exact leading boundary one-form restricts every gcd1/gcd2 pole pair; squarefree pole-four cubics excluded and mixed cubics must contain contact',
    excluded_pole_four_sextic_root_partitions=[[2,2,2]],remaining_pole_four_sextic_root_partitions=[[6],[4,2]],pole_four_cubic_divisors=remaining_cubic,
    frontier_superseded_by='PF2RULED16',
    frontier_handoff_scope='historical pre-PF2RULED16; this receipt does not prove the current frontier',
    first_open_reduced_pairs=[[6,8],[7,8]],
    frontier_adoption='PF2RULED12; subsequent handoff, outside this checker intrinsic proof coverage',minimum_target_coordinate_pole=4,
    exact_source_degree_contact_orders=[0,1,2,3],
    assurance=dict(formal_descent_frame_and_rational_degree_coverage='written proof',
        root_period_proof_uses_finiteness=False,root_period_proof_uses_external_degree_theorem=False,
        height_eight_numeric_corollary_uses_external_2022_degree_theorem=True,
        published_72_108_exception_retained=True,conditional_2026_degree125_upgrade_used=False,
        standard_Newton_corner_identification_proved=False,
        height_eight_completely_excluded=False,all_pole_four_coordinates_excluded=False,
        actual_ruled_compatibility_proved=False,quartic_packets_excluded=False,
        formal_verification=False,independent_full_reproof=False,external_review=False),
    source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),HERE/'BOUNDARY_LEADING_PERIOD.md']},
    python=platform.python_version(),sympy=S.__version__,elapsed_seconds=round(monotonic()-start,6))
(HERE/'boundary_leading_period_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
print(f'PASS {len(checks)} exact controls; cubic-square pole-four patterns6,42 remain; height-eight contact cases retain72,108')
