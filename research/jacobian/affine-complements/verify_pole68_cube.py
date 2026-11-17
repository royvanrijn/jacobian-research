"""PF2RULED13 exact controls; necessary profiles do not construct a finite map."""
from pathlib import Path
from time import monotonic
import sympy as S,json,hashlib,platform
HERE=Path(__file__).resolve().parent;start=monotonic();checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def zero(name,value):check(name,S.cancel(S.powsimp(S.expand(value),force=True))==0)
t,Z=S.symbols('t Z',positive=True);c,K=S.symbols('c K',positive=True,integer=True)
H=S.Function('H')(Z);G=S.Function('G')(Z)
P=t**(-6*c)*H;Q=t**(-K)*G
J=S.diff(P,t)*S.diff(Q,Z)-S.diff(P,Z)*S.diff(Q,t)
zero('cube weighted initial ODE',t**(6*c+K+1)*J-(K*S.diff(H,Z)*G-6*c*H*S.diff(G,Z)))
# All nonconstant degree<=5 cube numerators have exactly one triple root.
def partitions(n,minimum=1):
    if n==0:yield ()
    for a in range(minimum,n+1):
        for tail in partitions(n-a,a):yield (a,)+tail
cube_types=[]
for degree in range(1,6):
    for parts in partitions(degree):
        if all(4*k%3==0 for k in parts):cube_types.append(parts)
check('complete cube numerator multiplicity classification',cube_types==[(3,)])
B=S.Symbol('B',nonzero=True);sface=(1+B*Z)/Z**2
zero('cube and quartic faces commute',6*sface**3*S.diff(sface**4,Z)-8*S.diff(sface**3,Z)*sface**4)
check('linear-cube first face has nonzero z term',S.expand((1+B*Z)**3).coeff(Z,1)==3*B)
crit=S.Rational(41,12)
check('strict initial threshold is41/12',12*(crit-3)==5)
check('only integer cube slopes above threshold and <=6',list(range(int(S.floor(crit))+1,7))==[4,5,6])
bounds=[int(S.ceiling(18-crit*i)) for i in range(1,6)]
check('complete intervening-band coefficient bounds',bounds==[15,12,8,5,1])
slopes=sorted({S.Rational(18-h,i) for i in range(1,6) for h in range(bounds[i-1],3*(6-i))},reverse=True)
check('complete intervening-band slopes',slopes==[S.Rational(17,5),S.Rational(10,3),S.Rational(13,4),S.Rational(16,5)])
band=[]
for slope in slopes:
    a,b=int(slope.q),int(slope.p);pole=17*a-4*b
    check('all intervening-band faces exceed degree4 '+str(slope),pole>4)
    band.append(dict(slope=str(slope),Q_pole_lower_bound=pole))
check('exact intervening-band pole list',[r['Q_pole_lower_bound'] for r in band]==[17,11,16,21])
# Global cube completion only uses symbol degrees and exact products.
T,A=S.symbols('T A')
zero('global cubic completion exact expansion',(T*T+A)**3-T**6-3*T**4*A-3*T*T*A*A-A**3)
check('cube completion remainder has maximum pole4',max(2+2,3)==4)
strict=[]
for j in range(3):
    slope=6-j;row=[max(0,19-i*slope) for i in range(2,6)]
    strict.append(dict(j=j,first_slope=slope,remainder_z2_to_z5_orders=row))
check('complete strict cube-remainder rows',[r['remainder_z2_to_z5_orders'] for r in strict]==[[7,1,0,0],[9,4,0,0],[11,7,3,0]])
for j in (0,1):
    eta_order=7+2*j;u_order=eta_order-6
    check('global X2 symbol order '+str(j),u_order==2*j+1)
    check('global Y3 symbol order '+str(j),min(j+u_order,6)>=3*j+1)
# Slope4: first cluster has total index3, next slope lies in[2,3].
check('slope4 cubic cluster contributes at least3',3*(17-4*4)==3)
next_slopes=sorted({max(S.Rational(6-h4,1),S.Rational(6-h5,2),S.Rational(6,3)) for h4 in range(3,8) for h5 in range(0,8)},reverse=True)
check('slope4 complete next-edge slope list',next_slopes==[3,S.Rational(5,2),2])
for slope in next_slopes:
    a,b=int(slope.q),int(slope.p)
    check('slope4 next-edge pole at least2 '+str(slope),5*a-b>=2)
check('slope4 high plus low exceeds quartic degree',3+2>4)
# Low-edge controls for slope5, retaining all relevant next coefficients.
w,z,C,D,xi=S.symbols('w z C D xi')
for h,a,b,expected_pole in [(0,1,3,-1),(1,1,2,0)]:
    F=(w**6+w*z)**3+C*w**h*z**4+D*z**5+xi*z**6
    f=S.Poly(S.expand(F.subs({w:t**a,z:t**b*Z})),t)
    value=min(mon[0] for mon,coef in f.terms());face=f.coeff_monomial(t**value)
    residual=S.cancel(face/Z**3)
    check('slope5 forbidden low face has a simple root '+str(h),S.degree(residual,Z)==1 and S.diff(residual,Z)!=0)
    check('slope5 low coefficient obstruction '+str(h),value-b-a-3*b==expected_pole)
# The two retained low alternatives, with ordinary generic fibre parameter.
F=(w**6+w*z)**3+D*z**5+xi*z**6
f=S.Poly(S.expand(F.subs({w:t**2,z:t**3*Z})),t);value=min(mon[0] for mon,coef in f.terms())
check('retained low cusp Newton order',value==15)
check('retained low cusp Fz and Q pole',value-3==12 and 12-2-3*3==1)
check('retained low cusp parameter primitive',S.gcd(2,3)==1)
low=1+C*Z+D*Z*Z+xi*Z**3
disc_low=S.Poly(S.discriminant(low,Z),xi)
zero('generic low cubic discriminant leading coefficient',disc_low.LC()+27)
check('generic low cubic parameter genuinely changes discriminant',disc_low.degree()==2)
check('retained low smooth branch Fz and Q pole',5-1-3==1)
# Norm/discriminant accounting fixes the high cubic jet orders.
check('slope6 entire cubic discriminant61',4+3+3*18==61)
check('slope5 high cluster low-factor contribution9',3*3==9)
check('slope5 high signed pole formula subtracts39',9-3-3*15==-39)
check('slope5 low cusp forces high discriminant42',(4-1)+39==42)
check('slope5 low three smooth points force high discriminant40',(4-3)+39==40)
aa,bb=S.symbols('aa bb')
zero('universal depressed cubic discriminant',S.discriminant(Z**3+aa*Z+bb,Z)+4*aa**3+27*bb**2)
check('discriminant61 not an uncancelled square or cube order',61%2!=0 and 61%3!=0)
allowed_t=[q for q in range(1,11) if q>6 and 2*q-19>0]
check('slope6 complete positive smooth-pole cancellation list',allowed_t==[10])
check('slope6 unique coefficient orders20,30',2*allowed_t[0]==20 and 3*allowed_t[0]==30)
check('slope6 smooth branch pole1',2*10-19==1)
check('slope6 double-factor discriminant21',61-2*20==21)
check('slope6 double factor ramifies',21%2==1)
for target in (40,42):
    candidates=[q for q in range(1,(target-1)//6+1) if q>5]
    check('all premature high-cubic cancellations violate finiteness '+str(target),all(2*q-13<0 for q in candidates))
check('uncancelled discriminant42 coefficient bounds',S.ceiling(S.Rational(42,3))==14 and S.ceiling(S.Rational(42,2))==21)
check('uncancelled discriminant40 coefficient bounds',S.ceiling(S.Rational(40,3))==14 and 40//2==20 and 40%3!=0)
check('high(3,20) cusp is primitive',S.gcd(3,20)==1)
# Three local curve controls check different/pole consistency only.
eps=S.Symbol('eps',nonzero=True)
F0=(z+w**6)**3-3*w**20*(z+w**6)+2*w**30+eps*w**31
H1=(z+w**5)**3+w**21;L1=z*z+w**3
H2=(z+w**5)**3+w**20;L2=z**3+w**3
models=[]
for name,F,degree,expected_disc,expected_poles in [
 ('slope6_two_infinity_points',F0,3,61,[3,1]),
 ('slope5_low_cusp',H1*L1,5,63,[1,1,1,1]),
 ('slope5_low_three_smooth',H2*L2,6,64,[1,1,1,1])]:
    disc=S.Poly(S.discriminant(F,z),w)
    order=min(mon[0] for mon,coef in disc.terms())
    const=S.Poly(F.subs(z,0),w);const_order=min(mon[0] for mon,coef in const.terms())
    local_degree=min(mon[0] for mon,coef in S.Poly(F.subs(w,0),z).terms())
    check('local curve discriminant order '+name,order==expected_disc)
    check('local curve projection degree '+name,local_degree==degree)
    check('local curve constant order18 '+name,const_order==18)
    check('local curve pole-budget consistency only '+name,order-degree-3*const_order==sum(expected_poles)==4)
    models.append(dict(name=name,local_projection_degree=degree,constant_order=const_order,full_discriminant_order=order,pole_partition=expected_poles,Q_constructed=False,is_global_finite_Wright_map=False))
for Hh,Ll,dh,dl in [(H1,L1,42,3),(H2,L2,40,6)]:
    hd=S.Poly(S.discriminant(Hh,z),w);ld=S.Poly(S.discriminant(Ll,z),w);res=S.Poly(S.resultant(Hh,Ll,z),w)
    check('local high cubic discriminant '+str(dh),min(mon[0] for mon,coef in hd.terms())==dh)
    check('local low factor discriminant '+str(dl),min(mon[0] for mon,coef in ld.terms())==dl)
    check('high-low resultant order9 '+str(dh),min(mon[0] for mon,coef in res.terms())==9)
# Explicit ramified profile of the coalescent slope6 control.
spar=S.Symbol('spar');Apar=S.Symbol('Apar',nonzero=True)
sub={w:spar**2,z:-spar**12+spar**20+Apar*spar**21}
curve=S.Poly(S.expand(F0.subs(sub)),spar)
zero('slope6 cusp first splitting relation',curve.coeff_monomial(spar**62)-(3*Apar**2+eps))
fd=S.Poly(S.expand(S.diff(F0,z).subs(sub)),spar)
check('slope6 cusp derivative order41',min(mon[0] for mon,coef in fd.terms())==41)
check('slope6 cusp Q pole3',41-2-3*12==3)
profiles=[
 dict(first_slope=6,beta_order=0,high_cubic_a_order=20,high_cubic_b_order=30,high_discriminant_order=61,Q_poles=[3,1]),
 dict(first_slope=5,beta_order=1,fifth_remainder_coefficient='unit',high_cubic_a_order_at_least=14,high_cubic_b_order_at_least=21,high_discriminant_order=42,Q_poles=[1,1,1,1]),
 dict(first_slope=5,beta_order=1,fifth_remainder_coefficient='vanishes',high_cubic_a_order_at_least=14,high_cubic_b_order=20,high_discriminant_order=40,Q_poles=[1,1,1,1])]
report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
    eliminated_intervening_band=band,strict_cube_remainder_orders=strict,
    excluded_perfect_cube_first_slopes=[4],remaining_first_slopes=[5,6],necessary_profiles=profiles,local_curve_controls=models,
    profile_gate_superseded_by='PF2RULED14',
    written_result='Retained PF2RULED13 calculation: noncontact(6,8) requires a global cubic completion and one of three precise necessary Weierstrass/pole profiles; slope4 is excluded',
    frontier_superseded_by='PF2RULED16',
    frontier_handoff_scope='historical pre-PF2RULED16; this receipt does not prove the current frontier',
    first_open_reduced_pairs=[[6,8],[7,8]],minimum_target_coordinate_pole=4,
    assurance=dict(full_first_edge_and_local_profile_coverage='written proof',
        global_cube_completion_and_discriminant_interpretation='written proof',
        new_noncontact_cube_analysis_uses_external_degree_theorem=False,
        profiles_are_necessary_not_realizability_certificates=True,
        assurance_scope='PF2RULED13 local calculation before the PF2RULED14 global exclusion',
        local_curve_controls_construct_Q=False,required_fibre_residues_proved_zero=False,
        any_profile_realizable_as_global_finite_map_proved=False,noncontact_pole68_excluded=False,
        complete_height_eight_excluded=False,all_pole_four_coordinates_excluded=False,
        actual_ruled_compatibility_proved=False,quartic_packets_excluded=False,
        formal_verification=False,independent_full_reproof=False,external_review=False),
    source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),HERE/'POLE_SIX_EIGHT_CUBE_REDUCTION.md']},
    python=platform.python_version(),sympy=S.__version__,elapsed_seconds=round(monotonic()-start,6))
(HERE/'pole68_cube_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
print(f'PASS {len(checks)} exact controls; retained cube-profile controls pass; global noncontact exclusion is PF2RULED14')
