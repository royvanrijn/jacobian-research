"""PF2RULED19: generic-polar pole budgets and exact marked noncontact cusp sectors.
All-height geometry is written proof; exact finite identities, full section
jet blocks and two-root rational univariate representations are replayed here.
"""
from pathlib import Path
from time import monotonic
from math import gcd
import hashlib,json,platform
import sympy as S
import verify_contact_filtration as V
HERE=Path(__file__).resolve().parent
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def zero(name,value):check(name,S.cancel(S.expand(value))==0)
def partitions(n,start=1):
    if n==0:yield ()
    for k in range(start,n+1):
        for rest in partitions(n-k,k):yield (k,)+rest

def main():
    start=monotonic();m,l,J,I,d,ki=S.symbols('m l J I d ki',positive=True,integer=True)
    zero('local pole-sum formula',I-(m-2)*l-J-(I-m*l+2*l-J))
    zero('generic polar upper budget',m*(l-1)-(m-2)*l-J-(2*l-m-J))
    zero('common-symbol multiplicity threshold',2*m*ki/d-m-m*(2*ki/d-1))
    check('contact6 exceptional local pole budget is zero',2*6-6-(4*6-18)==0)
    check('contact1 local budget strictly negative',2*2-6-2<0)
    check('contact2 local budget strictly negative even at minimum E intersection',2*4-6-4<0)
    check('one simple common-nonic root cannot be finite',2*2-6<0 and 2*3-9<0)
    # Attained and strict generic-polar controls; these are local curves, not maps.
    z,w,c=S.symbols('z w c');polar_controls=[]
    for pole,order in [(6,2),(6,4),(6,6)]:
        ff=w**order-c*z**pole;disc=S.factor(S.discriminant(ff,w))
        zorder=S.Poly(disc,z).monoms()[0][0]
        check(f'balanced polar control order{order}',zorder==pole*(order-1))
        polar_controls.append(dict(numerator=str(ff),polar_intersection=zorder,
            Q_pole_sum_if_finite=zorder-(pole-2)*order,finite_map=False))
    check('strict polar bound control',5*5==25<6*5 and 25-(6-2)*6==1)
    polar_controls.append(dict(numerator='w^6+z^5-c*z^6',polar_intersection=25,upper_bound=30,Q_pole_sum_if_finite=1,finite_map=False))
    # Full partition census after finite-map root multiplicity restrictions.
    noncontact=sorted([list(reversed(p)) for p in partitions(9) if min(p)>=2 and not all(k%3==0 for k in p)],key=lambda p:(len(p),tuple(-k for k in p)))
    check('exact five noncontact noncube partitions',noncontact==[[7,2],[5,4],[5,2,2],[4,3,2],[3,2,2,2]])
    check('quartic pole degree supports at most four missing points',max(len(p) for p in partitions(4))==4)
    check('contact8 has a forbidden simple away root',9-8==1)
    # The retained elliptic contact3 sector has four distinct boundary roots.
    elliptic_profile=dict(contact=3,common_noncontact_multiplicities=[2,2,2],
        P_pole=6,P_restriction_degree=20,number_of_missing_points=4,Q_poles=[1,1,1,1],
        away_polar_intersections=[17,17,17],contact_polar_intersection=29,
        global_PQ_constructed=False,profile_excluded=False,
        profile_stage="PF2RULED19 necessary data before PF2RULED20")
    check('elliptic four-place pole budget',sum(elliptic_profile['Q_poles'])==4)
    check('elliptic away polar orders',4*(6-2)+1==17)
    check('elliptic contact polar order',6*(6-2)+(4*6-20)+1==29)
    # Every section of V6 and V9: certify the complete highest four jet coefficients.
    ap=V.powers({(0,0):1,(2,1):1},9);bp=V.powers({(0,0):-1,(2,1):-2,(4,2):-1,(0,1):-1},9);tp=V.powers({(0,0):1,(2,0):1,(4,1):1},18)
    jet_counts={}
    for pole in (6,9):
        count=0
        for jj in range(pole+1):
            for ii in range(pole+jj+1):
                source={(a+ii,b):V.F(v) for (a,b),v in V.mul(V.mul(ap[pole-jj],bp[jj]),tp[pole+jj-ii]).items()}
                cc=V.chart(source);check(f'whole actual section chart pole{pole} row{count}',all(a>=0 for a,b in cc))
                f0={b:v for (a,b),v in cc.items() if a==0};f1={b:v for (a,b),v in cc.items() if a==1}
                corrected=V.add(f1,{b-1:b*v for b,v in f0.items() if b},-V.F(3,4))
                power=2*(pole-jj)+ii;symbol=[V.F((-1)**(pole-jj)) if k==power else V.F(0) for k in range(5)]
                c0,c1,c2,c3,c4=symbol
                A=[c0,c1,c2-pole*c0,c3-pole*c1]
                B=[-c1/4,-(pole*c0+c2)/2,-(pole*c1+3*c3)/4,-c4+V.F(pole*(pole-1),2)*c0]
                for k in range(4):
                    check(f'full restriction high jet pole{pole} row{count} index{k}',f0.get(4*pole-k,0)==(-1)**pole*A[k])
                    check(f'full corrected high jet pole{pole} row{count} index{k}',corrected.get(4*pole-2-k,0)==(-1)**pole*B[k])
                count+=1
        check('complete section jet coverage pole'+str(pole),count==(pole+1)*(3*pole+2)//2)
        jet_counts[str(pole)]=count
    # Leading common nonic H=w^9+a w^8+b w^7+c w^6+d w^5+...
    a,b,cc,dd,t,delta,r=S.symbols('a b cc dd t delta r')
    h=1+a*r+b*r*r+cc*r**3+dd*r**4
    E=a*a*b+11*a*a+3*a*cc-4*b*b-24*b-36
    F=a**3*b+16*a**3+5*a*a*cc-4*a*b*b-30*a*b+4*a*dd-66*a-6*b*cc-18*cc
    deltav=-a/48;tv=(24*b+72-11*a*a)/(12*a)
    for pole,exponent in [(6,2),(9,3)]:
        pp=S.Poly(S.expand(h**exponent),r);cs=[pp.nth(k) for k in range(5)];D=4*pole
        AA=[1,cs[1],cs[2]-pole,cs[3]-pole*cs[1]]
        BB=[-cs[1]/4,-(pole+cs[2])/2,-(pole*cs[1]+3*cs[3])/4,-cs[4]+S.Rational(pole*(pole-1),2)]
        residual=[S.factor(BB[k]-(t*BB[k-1] if k else 0)-delta*(D-k)*AA[k]) for k in range(4)]
        reduced=[S.factor(x.subs({delta:deltav,t:tv})) for x in residual]
        zero('common nonic delta constraint pole'+str(pole),reduced[0])
        zero('common nonic cusp parameter pole'+str(pole),reduced[1])
        zero('common nonic third constraint pole'+str(pole),reduced[2]+(E/(2*a) if pole==6 else 3*E/(4*a)))
        zero('common nonic fourth constraint pole'+str(pole),reduced[3]+(F/(2*a) if pole==6 else 3*(F+a*E)/(4*a)))
    # Two-root sector: k is the smaller multiplicity, roots x,y are fixed marked coordinates.
    x,y=S.symbols('x y');sectors=[]
    parameters=[(2,21*y*y-4,1225*y**8-1568*y**6+5600*y**4-13440*y*y+8192),
                (4,15*y*y-8,25*y**8+2800*y**6-4480*y**4+1920*y*y-256)]
    for k,q2,p8 in parameters:
        HH=S.Poly((w-x)**k*(w-y)**(9-k),w);mp={a:HH.nth(8),b:HH.nth(7),cc:HH.nth(6),dd:HH.nth(5)}
        ee=S.expand(E.subs(mp));ff=S.expand(F.subs(mp));amean=mp[a]
        resultant=S.factor(S.resultant(ee,ff,x))
        check('complete resultant factorization sector'+str(k),S.cancel(resultant/(y*q2*p8)).is_Rational and resultant!=0)
        # y=0 is impossible: exact Bezout identity over Q[x].
        s0,t0,g0=S.gcdex(ee.subs(y,0),ff.subs(y,0),x)
        zero('y0 unit certificate sector'+str(k),s0*ee.subs(y,0)+t0*ff.subs(y,0)-1)
        check('y0 gcd is unit sector'+str(k),g0==1)
        # q2=0 forces the forbidden a=0. Exact Groebner basis is replayed,
        # with both input reduction and the explicit linear output retained.
        G=S.groebner([ee,ff,q2],x,y,order='lex')
        check('spurious quadratic forces zero cusp unit sector'+str(k),list(G)==[k*x+(9-k)*y,q2])
        check('spurious quadratic a vanishes sector'+str(k),G.reduce(amean)[1]==0)
        sub=S.subresultants(ee,ff,x);linear=next(q for q in reversed(sub) if S.degree(q,x)==1)
        lin=S.Poly(linear,x);lc,lc0=lin.nth(1),lin.nth(0)
        check('linear subresultant denominator is a unit sector'+str(k),S.gcd(lc,p8)==1)
        xp=S.rem(-lc0*S.invert(lc,p8,y),p8,y)
        zero('linear subresultant substitution sector'+str(k),S.rem(linear.subs(x,xp),p8,y))
        zero('first marked jet equation sector'+str(k),S.rem(ee.subs(x,xp),p8,y))
        zero('second marked jet equation sector'+str(k),S.rem(ff.subs(x,xp),p8,y))
        check('eight distinct candidate roots sector'+str(k),S.degree(p8,y)==8 and S.gcd(p8,S.diff(p8,y))==1)
        units={}
        for label,guard in [('a',S.rem(amean.subs(x,xp),p8,y)),('root_separation',xp-y),('y',y),('x',xp)]:
            uu,vv,gg=S.gcdex(guard,p8,y)
            check('marked guard is unit '+str(k)+' '+label,gg==1)
            zero('explicit marked guard Bezout '+str(k)+' '+label,uu*guard+vv*p8-1)
            units[label]={'guard':str(guard),'guard_coefficient':str(uu),'polynomial_coefficient':str(vv)}
        sectors.append(dict(multiplicities=[k,9-k],H=f'(w-x)^{k}*(w-y)^{9-k}',
            E=str(ee),F=str(ff),resultant=str(resultant),discarded_quadratic=str(q2),
            discarded_quadratic_Groebner_basis=list(map(str,G)),
            y0_Bezout_coefficients=[str(s0),str(t0)],
            P8=str(p8),x_as_polynomial_in_y=str(S.factor(xp)),linear_subresultant=str(S.factor(linear)),
            units=units,distinct_marked_coefficient_points=8,
            complete_leading_period_and_high_jet_sector=True,global_finite_map=False))
    frontier=[dict(contacts=[0],ordinary_degrees=[96,144]),
              dict(contacts=[3,4,5,6,7,9],ordinary_degrees=[80,128]),
              dict(contacts=[7,9],ordinary_degrees=[48,128])]
    cert=dict(as_of='2026-11-14',generic_polar_formula='sum nu=I(F,A_w)-(m-2)*l-I(F,E)<=2*l-m-I(F,E)',
        support_bound=4,common_symbol_root_rule='2*k>d at every root of common section in L^d',
        polar_controls=polar_controls,noncontact_nonic_partitions=noncontact,
        current_pole69_frontier=frontier,excluded_contacts=[1,2,8],excluded_pole69_ordinary_pairs=[[88,132],[72,108]],
        elliptic_contact3_necessary_profile=elliptic_profile,
        complete_section_jet_counts=jet_counts,
        noncontact_marked_cusp=dict(delta=str(deltav),cusp_parameter=str(tv),required_nonzero='a',E=str(E),F=str(F)),
        two_root_marked_sectors=sectors,total_distinct_marked_two_root_points=16,
        finite_map_constructed=False)
    path=HERE/'polar_cusp_sector_certificates.json';path.write_text(json.dumps(cert,indent=2)+'\n')
    inputs=[Path(__file__),HERE/'POLAR_CUSP_LEADING_SECTORS.md',path,HERE/'verify_contact_filtration.py',
        HERE/'CUBIC_COVER_BOUNDARY_PERIOD.md',HERE/'cubic_cover_period_receipt.json',HERE/'GLOBAL_CUSP_LEADING_OBSTRUCTION.md',
        HERE/'SOURCE_COORDINATES.md',HERE/'WRIGHT_MARKED_MODEL.md',HERE/'POLE_FIVE_EIGHT_EXCLUSION.md',
        HERE.parent/'jc2-omitted-divisors/WRIGHT_IMAGE_RESTRICTION.md']
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
        written_result='Generic polar budgets forbid simple common nonic roots, remove pole69 contact1/contact8 and72,108, and mark the two-root noncontact sector by sixteen exact leading candidates',
        all_height_generic_polar_budget_proved=True,all_height_common_root_multiplicity_rule_proved=True,
        root_support_bound_for_finite_quartic=4,current_pole69_frontier=frontier,
        current_noncontact_nonic_partitions=noncontact,pole69_72_108_excluded=True,
        marked_two_root_sector_exact_point_count=16,
        elliptic_profile_superseded_by='PF2RULED20',
        noncontact_partition_frontier_superseded_by='PF2RULED20',
        profile_scope='PF2RULED19 necessary data before the PF2RULED20 generic quartic gap',
        first_open_reduced_pairs=[[5,9],[6,9],[7,9],[8,9]],minimum_compatible_quartic_pole_height=9,
        assurance=dict(all_height_polar_and_finiteness_coverage='written proof',
            full_section_high_jet_coverage=True,marked_two_root_classification_is_only_leading_tests=True,
            period_solution_implies_finite_map=False,sixteen_points_construct_actual_maps=False,
            all_noncontact_pole69_excluded=False,elliptic_contact3_profile_excluded=False,
            all_pole69_excluded=False,full_degree9_period_position_classification=False,
            all_height_actual_72_108_exception_excluded=False,actual_ruled_compatibility_proved=False,
            quartic_packets_excluded=False,conditional_2026_degree125_upgrade_used=False,
            formal_verification=False,independent_full_reproof=False,external_review=False),
        source_sha256={str(f.relative_to(HERE)) if f.is_relative_to(HERE) else '../jc2-omitted-divisors/'+f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs},
        python=platform.python_version(),sympy=S.__version__,elapsed_seconds=round(monotonic()-start,6))
    (HERE/'polar_cusp_sector_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} exact controls; three pole69 source groups remain;16 two-root leading candidates; full goal active')

if __name__=='__main__':main()
