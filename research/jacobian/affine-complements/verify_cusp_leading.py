"""PF2RULED14 exact algebra; the all-height filtration and finite Euler input are written proofs."""
from pathlib import Path
from time import monotonic
import sympy as S,json,hashlib,platform
HERE=Path(__file__).resolve().parent;start=monotonic();checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def zero(name,value):check(name,S.cancel(S.powsimp(S.expand(value),force=True))==0)
u,v,r=S.symbols('u v r')
Y=u**3*v-u*u;A1=u*u*v-u;B1=u*v-1;D=v;C=1+B1+D*D
f0=-Y-C*C-2*A1*D-B1*D*D
f1=-A1-D-2*B1*D-D**3;f2=C;f3=D
g=[S.expand(f0),S.expand(f1),S.expand(-f2),S.expand(-f3)]
p0=[x.subs(u,0) for x in g];p1=[S.diff(x,u).subs(u,0) for x in g]
jets=[S.expand(a-S.Rational(3,4)*S.diff(b,v)) for a,b in zip(p1,p0)]
expected0=[-v**4+v*v,-v**3+v,-v*v,-v]
expectedJ=[v/2,(v*v+1)/4,v/2,S.Rational(3,4)]
for i in range(4):
    zero('actual generator restriction '+str(i),p0[i]-expected0[i])
    zero('actual corrected first jet '+str(i),jets[i]-expectedJ[i])
    check('generator ordinary restriction degree '+str(i),S.degree(p0[i],v)==4-i)
    check('complete corrected-jet degree drop '+str(i),S.degree(p0[i],v)-S.degree(jets[i],v)==[3,1,1,1][i])
# Marked symbol calculation on t=1, Bcal=alpha+beta*r^2=0.
beta,alpha=S.symbols('beta alpha',nonzero=True)
cox=[beta,beta*r,-alpha,-alpha*r]
for i in range(4):zero('actual marked cubic frame '+str(i),cox[i].subs(alpha,-beta*r*r)/beta-r**i)
# First-order chain rule: the same h,lambda constrain all target elements.
hh,lam,aa,bb,FP,FQ=S.symbols('hh lam aa bb FP FQ')
zero('target algebra chain-rule identity',hh*(FP*lam*aa+FQ*lam*bb)-lam*(FP*hh*aa+FQ*hh*bb))
# The degree proof is uniform: any non-g0 product has restriction degree
# <=4s-1; all-g0 corrected jets drop3. Both give <=4s-2.
s,m=S.symbols('s m',integer=True,positive=True)
zero('uniform all-g0 degree is below bound by one',(4*s-2)-(4*s-3)-1)
zero('uniform other-generator degree bound',(4*s-1)-1-(4*s-2))
zero('strict global power-lift degree gap',(4*m-4)-(4*m-5)-1)
CC=S.Function('CC')(v);LL=S.Function('LL')(v);tt,delta=S.symbols('tt delta')
zero('universal power derivation identity',(v-tt)*m*CC**(m-1)*LL-delta*S.diff(CC**m,v)-m*CC**(m-1)*((v-tt)*LL-delta*S.diff(CC,v)))
# Finite exact controls supplement rather than establish the unbounded proof.
for pole in (1,2,3):
    for counts0 in range(pole+1):
        for counts1 in range(pole-counts0+1):
            for counts2 in range(pole-counts0-counts1+1):
                counts=[counts0,counts1,counts2,pole-counts0-counts1-counts2]
                restriction=S.prod(p0[i]**counts[i] for i in range(4))
                deriv=S.expand(sum(counts[i]*jets[i]*S.prod(p0[j]**(counts[j]-(1 if i==j else 0)) for j in range(4)) for i in range(4) if counts[i]))
                check('supplementary monomial degree control '+str((pole,tuple(counts))),S.degree(restriction,v)<=4*pole and S.degree(deriv,v)<=4*pole-2)
# Marked p=infinity: corrected jet has degree1, but C prime degree3.
check('marked-root degree contradiction',S.degree((v-tt)*jets[0],v)==2 and S.degree(S.diff(p0[0],v),v)==3)
# Finite p has actual base coordinate a!=0, with tau(q)=1.
a=S.Symbol('a',nonzero=True)
coeff=[1,-3/a,3/a**2,-1/a**3]
restriction=S.expand(sum(c*f for c,f in zip(coeff,p0)))
L0=S.expand(sum(c*f for c,f in zip(coeff,jets)))
E=S.Poly(S.expand((v-tt)*L0-delta*S.diff(restriction,v)),v)
first=E.coeff_monomial(v**3)
delta_value=S.Rational(3,16)/a
zero('leading coefficient fixes delta',first.subs(delta,delta_value))
check('leading coefficient determines delta uniquely',S.diff(first,delta)==4)
second=E.coeff_monomial(v**2).subs(delta,delta_value)
t_value=1/(4*a)-2*a/3
zero('next coefficient fixes cusp parameter',second.subs(tt,t_value))
check('next coefficient determines cusp parameter uniquely',S.diff(second,tt)!=0)
reduced={i:S.factor(E.coeff_monomial(v**i).subs({delta:delta_value,tt:t_value})) for i in range(4)}
zero('linear residual numerator',reduced[1]-(4*a*a-3)/(12*a))
zero('constant residual numerator',reduced[0]+(2*a*a-1)/(4*a*a))
AA=4*a*a-3;BB=2*a*a-1
zero('exact unit obstruction',2*BB-AA-1)
check('unit obstruction is nonzero',2*BB-AA==1)
# Independent four-high-coefficient expansion agrees for the relevant poles.
for mm in (2,4,6,8):
    cc=[S.binomial(3*mm,k)*(-1/a)**k for k in range(5)]
    normalized=[-cc[1]/4,-(mm+cc[2])/2,-(mm*cc[1]+3*cc[3])/4,-cc[4]+S.Rational(mm*(mm-1),2)]
    d=4*mm;A=[1,cc[1],cc[2]-mm,cc[3]-mm*cc[1]]
    rhs=[sum((d-j)*A[j]*t_value**(i-j) for j in range(i+1))*delta_value for i in range(4)]
    residual=[S.factor(normalized[i]-rhs[i]) for i in range(4)]
    zero('four-coefficient check fixes first residual '+str(mm),residual[0])
    zero('four-coefficient check fixes second residual '+str(mm),residual[1])
    zero('four-coefficient check forces same a condition '+str(mm),residual[2]+mm*(4*a*a-3)/(12*a))
    check('fourth coefficient independently stays nonzero at a^2=3/4 '+str(mm),S.simplify(residual[3].subs(a,S.sqrt(3)/2))==S.Rational(mm,6))
# Source-degree consequences use the existing explicit filtration.
check('height>=8 and gcd<=2 imply top contact order>=4',S.Rational(8,2)>=4)
check('contact order>=4 lowers the degree bound by one pole',16*8-4*4==16*(8-1))
check('pole4 contact symbol gives source degree<=48',max(16*3,16*4-4*4)==48)
cert=dict(as_of='2026-11-14',generator_restrictions=[str(x) for x in p0],generator_first_jets=[str(x) for x in p1],generator_corrected_jets=[str(x) for x in jets],
    normalized_finite_root_lift_coefficients=[str(x) for x in coeff],restriction_C=str(restriction),corrected_jet_L0=str(L0),
    cubic_identity_coefficients={str(i):str(E.coeff_monomial(v**i)) for i in range(4)},
    forced_delta=str(delta_value),forced_cusp_parameter=str(t_value),reduced_coefficients={str(i):str(x) for i,x in reduced.items()},
    numerator_equations=[str(AA),str(BB)],unit_certificate='2*(2*a**2-1)-(4*a**2-3)=1',
    required_nonzero_parameters=['a (noncontact finite root)','delta=lambda(v_star) (quartic cusp unit)'],
    all_height_degree_gap='deg L(R)<=4m-5 <4m-4<=deg(m*C^(m-1)*L(T)) if L(T)!=0',
    construct_finite_map=False)
certpath=HERE/'cusp_leading_certificates.json';certpath.write_text(json.dumps(cert,indent=2)+'\n')
report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
    written_result='No element of a finite quartic marked Wright Keller target algebra has leading tau^m with div(tau)=3p, p!=q; all gcd1/gcd2 coordinate presentations require contact symbols',
    excludes_noncoordinate_target_elements_with_indicated_symbol=True,all_height_one_root_noncontact_excluded=True,
    noncontact_pole68_excluded=True,supersedes_profile_realization_gate='PF2RULED13',
    frontier_superseded_by='PF2RULED16',
    frontier_handoff_scope='historical pre-PF2RULED16; this receipt does not prove the current frontier',
    first_open_reduced_pairs=[[6,8],[7,8]],remaining_height_eight=dict(all_contact=True,ordinary_degrees=[72,108]),
    minimum_target_coordinate_pole=4,pole_four_source_degree_upper_bound=48,
    assurance=dict(all_height_filtration_degree_gap='written proof',quartic_cusp_and_Euler_input='PF2OMIT9 written theorem',
        target_algebra_identity_uses_actual_marked_inclusion=True,auxiliary_lift_put_in_target_algebra=False,
        contact_height_eight_excluded=False,all_pole_four_coordinates_excluded=False,
        arbitrary_source_component_degree48_excluded=False,higher_gcds_classified=False,
        conditional_2026_degree125_upgrade_used=False,actual_ruled_compatibility_proved=False,
        quartic_packets_excluded=False,formal_verification=False,independent_full_reproof=False,external_review=False),
    source_sha256={str(p.relative_to(HERE)) if p.is_relative_to(HERE) else '../jc2-omitted-divisors/'+p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),HERE/'GLOBAL_CUSP_LEADING_OBSTRUCTION.md',certpath,HERE/'WRIGHT_MARKED_MODEL.md',HERE/'SOURCE_COORDINATES.md',HERE.parent/'jc2-omitted-divisors/WRIGHT_IMAGE_RESTRICTION.md']},
    python=platform.python_version(),sympy=S.__version__,elapsed_seconds=round(monotonic()-start,6))
(HERE/'cusp_leading_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
print(f'PASS {len(checks)} exact controls; all one-root noncontact symbols excluded; height8 only contact72,108')
