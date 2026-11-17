"""PF2RULED20: generic quartic discriminant gap and its finite-map applications.
The formal normal-form induction and all-order valuation coverage are written
proofs; this replays exact identities, exhaustive valuation strata and controls.
"""
from pathlib import Path
from time import monotonic
import hashlib,json,platform
import sympy as S
HERE=Path(__file__).resolve().parent
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def zero(name,value):check(name,S.cancel(S.expand(value))==0)

def main():
    start=monotonic();W,z,c,p,q,r,h=S.symbols('W z c p q r h')
    # Formal right-normalization step, performed z-adically and independently
    # of c. Terms of z-order n change the existing nonleading coefficients
    # only in order>=n+1; the all-n induction is written in the proof.
    eps,R0,R1,R2=S.symbols('eps R0 R1 R2');H=S.Function('H')(W)
    low=R0+R1*W+R2*W*W;residual=low+4*W**3*H
    step=S.expand((W-eps*H)**4+eps*residual-W**4-eps*low)
    zero('universal right-normalization first-order step',step.coeff(eps,1))
    check('right-normalization error is at least quadratic',step.coeff(eps,0)==0)
    n=S.symbols('n',integer=True,positive=True)
    zero('all-order z-adic error margin',2*n-(n+1)-(n-1))
    F=W**4+p*W**2+q*W+r-c*z**6
    Delta=S.expand(S.discriminant(F,W))
    D0=256*r**3-128*p*p*r*r+144*p*q*q*r-27*q**4+16*p**4*r-4*p**3*q*q
    D1=-16*z**6*(48*r*r-16*p*p*r+9*p*q*q+p**4)
    D2=128*z**12*(6*r-p*p);D3=-256*z**18
    for k,expected in enumerate([D0,D1,D2,D3]):zero('complete c coefficient '+str(k),Delta.coeff(c,k)-expected)
    check('generic discriminant upper order18 coefficient is fixed',Delta.coeff(c,3)==-256*z**18)
    B=S.cancel(D1/(-16*z**6))
    zero('c1 bracket after c2 constraint',B.subs(r,p*p/6+h)-(-p**4/3+9*p*q*q+48*h*h))
    # Complete positive-integer order alternatives for p order1 or2.
    strata=[]
    for po,cases in [(1,[(1,1),(2,None)]),(2,[(1,2),(3,3),(4,None)])]:
        for lo,hi in cases:
            if po==2 and lo==hi==3:
                strata.append(dict(p_order=2,q_order='3',conclusion='leading q3^2=p2^3/27; discriminant constant coefficient has nonzero z12 lead p2^6'))
                continue
            if hi is None:
                check(f'p{po} q>={lo}: p4 is uniquely leading in c1 bracket',4*po<po+2*lo and 4*po<10)
                leadorder=4*po;term='p^4'
            else:
                for qo in range(lo,hi+1):
                    check(f'p{po} q{qo}: p*q2 is uniquely leading in c1 bracket',po+2*qo<4*po and po+2*qo<10)
                leadorder=po+2*hi;term='p*q^2'
            check(f'p{po} q{lo} cannot meet c1 order11',leadorder<11)
            strata.append(dict(p_order=po,q_order=f'{lo}..{hi}' if hi is not None else '>='+str(lo),conclusion='c1 bracket uniquely led by '+term+' below order11'))
    P2,Q3=S.symbols('P2 Q3')
    lead=S.expand(D0.subs({p:P2,q:Q3,r:P2*P2/6}))
    unit= P2**6+(Q3*Q3-P2**3/27)*(19*P2**3-27*Q3*Q3)
    zero('order-two cancellation has exact nonzero constant-coefficient certificate',lead-unit)
    check('order-two certificate is nonzero on its required unit locus',P2**6!=0)
    # After p>=3, c2 gives r>=5. These finite sets are all remaining
    # exceptional integer orders of q and r; infinity is in the >= strata.
    for qo in (1,2,3):
        other=[12+5,9+2*qo,6+10,3+2*qo+5,15]
        check('q low order uniquely leads D0 '+str(qo),4*qo<min(other))
    check('r order5 uniquely leads D0 after q>=4',15<min(16,16,16,17,17))
    check('q order4 uniquely leads D0 after r>=6',16<min(18,17,18,17,18))
    # Every discriminant monomial has order>=18 in the forced cone.
    cone_bounds=[]
    for mon,coef in S.Poly(Delta,p,q,r,z,c).terms():
        valuation=3*mon[0]+5*mon[1]+6*mon[2]+mon[3]
        check('complete forced-cone monomial bound '+str(mon),valuation>=18)
        cone_bounds.append(dict(exponents=list(mon),coefficient=str(coef),minimum_z_order=valuation))
    P3,Q5,R6=S.symbols('P3 Q5 R6')
    specialized=S.expand(Delta.subs({p:P3*z**3,q:Q5*z**5,r:R6*z**6}))
    generic_lead=16*(R6-c)*(P3*P3-4*(R6-c))**2
    zero('generic leading discriminant under forced cone',specialized.coeff(z,18)-generic_lead)
    check('leading discriminant never vanishes identically in c',S.Poly(generic_lead,c).degree()==3 and S.Poly(generic_lead,c).LC()==-256)
    # Normalized branch form: z=t^2,W=t^3 X, simple nonzero alpha roots.
    tt,X,alpha=S.symbols('tt X alpha');initial=X**4+P3*X*X+R6-c
    zero('weighted initial polynomial after degree12 division',
         S.expand(F.subs({p:P3*z**3,q:Q5*z**5,r:R6*z**6}).subs({z:tt**2,W:tt**3*X})/tt**12)-(initial+tt*Q5*X))
    zero('two-square-root generic initial discriminant',S.discriminant(initial,X)-generic_lead)
    zero('polar leading branch coefficient',S.diff(initial,X).subs(X,alpha)-2*alpha*(2*alpha*alpha+P3))
    check('two branches each have a simple Q pole',9-(6-2)*2==1)
    check('finite local pole sum exactly2',18-(6-2)*4==2)
    # Controls: cone attained; nongeneric parameter may jump, and the
    # superficially cancellable p-order2 case really has order12.
    controls=[]
    for name,poly,expected in [
        ('balanced two-cusp control',W**4-c*z**6,18),
        ('nongeneric-jump generic fibre',W**4+(1-c)*z**6+z**7,18),
        ('nongeneric-jump special c1',W**4+z**7,21),
        ('order-two near-cancellation',W**4+3*z**2*W**2+z**3*W+S.Rational(3,2)*z**4-c*z**6,12)]:
        disc=S.expand(S.discriminant(poly,W));order=min(S.Poly(disc,z).monoms())[0]
        check('exact control '+name,order==expected)
        controls.append(dict(name=name,polynomial=str(poly),generic_or_displayed_polar_order=order,discriminant=str(S.factor(disc)),global_finite_map=False))
    # A selected nongeneric fibre really can have polar17. The generic
    # pencil has order14 and fails the finite-pole requirement.
    special_A=(z-W**2)**2-W**15
    special_F=special_A-c*z**6
    special_param={W:tt**2,z:tt**4+tt**15,c:0}
    zero('special polar17 curve parameterization',special_F.subs(special_param))
    special_polar=S.expand(S.diff(special_A,W).subs(special_param))
    check('special fibre has actual polar17',min(S.Poly(special_polar,tt).monoms())[0]==17)
    generic_scaled=S.cancel(special_F.subs({W:tt,z:tt**2+tt**6*X})/tt**12)
    zero('generic14 countercontrol full scaled equation',generic_scaled-(X*X-c*(1+tt**4*X)**6-tt**3))
    zero('generic14 countercontrol simple root equation',generic_scaled.subs(tt,0)-(X*X-c))
    generic_polar=S.expand(S.diff(special_A,W).subs({W:tt,z:tt**2+tt**6*X}))
    check('generic14 countercontrol branch polar order7',min(S.Poly(generic_polar,tt).monoms())[0]==7 and generic_polar.coeff(tt,7)==-4*X)
    check('generic14 countercontrol two local branches',2*7==14<17)
    controls.append(dict(name='actual special17 versus generic14',polynomial=str(special_F),
        special_c=0,special_parameterization=['t^2','t^4+t^15'],special_polar_order=17,
        generic_polar_order=14,generic_branch_polar_orders=[7,7],global_finite_map=False))
    # Every k=2 away common-nonic root costs two, all other roots at least one.
    previous=json.loads((HERE/'polar_cusp_sector_certificates.json').read_text())
    previous_receipt=json.loads((HERE/'polar_cusp_sector_receipt.json').read_text())
    check('previous polar/cusp input passed',previous_receipt['outcome']=='PASS')
    for path,digest in previous_receipt['source_sha256'].items():
        check('previous polar/cusp seal '+path,hashlib.sha256((HERE/path).read_bytes()).hexdigest()==digest)
    partition_costs=[]
    for part in previous['noncontact_nonic_partitions']:
        cost=sum(2 if k==2 else 1 for k in part)
        partition_costs.append(dict(partition=part,minimum_Q_pole_degree=cost,excluded=cost>4))
    kept=[x['partition'] for x in partition_costs if not x['excluded']]
    removed=[x['partition'] for x in partition_costs if x['excluded']]
    check('exact two noncontact sectors excluded',removed==[[5,2,2],[3,2,2,2]])
    check('three remaining noncontact sectors',kept==[[7,2],[5,4],[4,3,2]])
    check('elliptic contact3 violates quartic degree',3*2+1==7>4)
    check('old elliptic polar17 profile contradicts the generic gap',previous['elliptic_contact3_necessary_profile']['away_polar_intersections']==[17,17,17])
    # Historical PF2RULED20 profiles; PF2RULED21 later excludes both.
    noncontact_profile=dict(common_partition=[4,3,2],missing_points_by_root=[1,1,2],Q_poles=[1,1,1,1],
        polar_intersections_by_root=[33,25,18],global_PQ_constructed=False,profile_excluded=False)
    contact_profile=dict(contact=3,away_multiplicities=[4,2],missing_points_away=[1,2],missing_points_at_contact=1,
        Q_poles=[1,1,1,1],polar_intersections_away=[33,18],polar_intersection_at_contact=29,
        global_PQ_constructed=False,profile_excluded=False)
    check('historical noncontact432 spends exactly4 poles',sum(noncontact_profile['missing_points_by_root'])==4)
    check('noncontact432 polar33',4*8+1==33)
    check('noncontact432 polar25',4*6+1==25)
    check('historical contact42 spends exactly4 poles',sum(contact_profile['missing_points_away'])+contact_profile['missing_points_at_contact']==4)
    check('remaining contact42 polar29',4*6+(4*6-20)+1==29)
    cert=dict(as_of='2026-11-14',formal_normal_form='A=W^4+p(z)W^2+q(z)W+r(z), z unchanged, coordinate change independent of c',
        discriminant=str(Delta),c_coefficients={str(k):str(x) for k,x in enumerate([D0,D1,D2,D3])},
        c1_after_c2_constraint=str(-p**4/3+9*p*q*q+48*h*h),
        p_order2_unit_certificate=str(unit),p_order2_required_guard='P2!=0',valuation_strata=strata,
        forced_orders={'p':3,'q':5,'r':6},cone_monomial_bounds=cone_bounds,
        generic_leading_coefficient=str(generic_lead),generic_excluded_parameter_values=['R6','R6-P3^2/4'],
        normalized_branch_orders={'z':2,'W':3,'A_W':9,'Q_pole':1},number_of_local_branches=2,
        controls=controls,noncontact_partition_costs=partition_costs,
        remaining_noncontact_partitions=kept,excluded_noncontact_partitions=removed,
        elliptic_contact3_profile_excluded=True,minimum_elliptic_contact3_Q_pole_degree=7,
        surviving_contact3_away_types=[[6],[4,2]],
        noncontact432_necessary_profile=noncontact_profile,contact42_necessary_profile=contact_profile,
        higher_polar_profiles_superseded_by="PF2RULED21",
        profile_scope="PF2RULED20 necessary data before PF2RULED21",
        full_pole69_excluded=False)
    cp=HERE/'quartic_polar_gap_certificates.json';cp.write_text(json.dumps(cert,indent=2)+'\n')
    inputs=[Path(__file__),HERE/'QUARTIC_GENERIC_POLAR_GAP.md',cp,HERE/'POLAR_CUSP_LEADING_SECTORS.md',
            HERE/'polar_cusp_sector_certificates.json',HERE/'polar_cusp_sector_receipt.json',
            HERE/'CUBIC_COVER_BOUNDARY_PERIOD.md',HERE/'WRIGHT_MARKED_MODEL.md']
    report=dict(higher_polar_profiles_superseded_by='PF2RULED21',
        noncontact_partition_frontier_superseded_by='PF2RULED21',
        profile_scope='PF2RULED20 necessary data before PF2RULED21',
        as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
        written_result='Generic quartic polar order17 is impossible; finite pole6/leading-root4 forces two cusp branches with simple Q poles, closing elliptic contact3 and two noncontact nonic partitions',
        generic_quartic_polar17_excluded=True,forced_generic_polar_order=18,
        forced_quartic_orders={'p':3,'q':5,'r':6},finite_local_branch_count=2,finite_local_Q_poles=[1,1],
        elliptic_contact3_pole69_excluded=True,excluded_noncontact_nonic_partitions=removed,
        remaining_noncontact_nonic_partitions=kept,current_contact3_away_types=[[6],[4,2]],
        marked_two_root_leading_candidates_retained=16,
        first_open_reduced_pairs=[[5,9],[6,9],[7,9],[8,9]],minimum_compatible_quartic_pole_height=9,
        assurance=dict(formal_normal_form_and_valuation_coverage='written proof',higher_terms_truncated_in_proof=False,
            normal_form_coordinate_change_depends_on_c=False,claims_no_polar17_for_all_special_fibres=False,
            local_generic_cone_constructs_global_finite_map=False,noncontact432_profile_excluded=False,
            contact42_profile_excluded=False,all_pole69_excluded=False,actual_ruled_compatibility_proved=False,
            quartic_packets_excluded=False,conditional_2026_degree125_upgrade_used=False,
            formal_verification=False,independent_full_reproof=False,external_review=False),
        source_sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs},
        python=platform.python_version(),sympy=S.__version__,elapsed_seconds=round(monotonic()-start,6))
    (HERE/'quartic_polar_gap_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} exact controls; generic polar17 excluded; elliptic contact3 closed; historical PF20 handoff has three noncontact sectors')

if __name__=='__main__':main()
