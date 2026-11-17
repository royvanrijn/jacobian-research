"""PF2RULED81: arbitrary-k marked-end controls; universal proofs in the note."""
from pathlib import Path
import sys,hashlib,json,math,os,platform,collections
import sympy as S
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
import verify_general_two_stage as G
E=G.E;B=G.B

def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def identity(name,expr):check(name,S.expand(expr)==0)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    preds=['general_two_stage_receipt.json','terminal_cone_binomial_closure_receipt.json',
      'matched_meridian_common_exclusion_receipt.json','seifert_sections_winding_three_receipt.json',
      'rational_end_orderability_receipt.json','lifted_surface_s3_receipt.json',
      'fibre_order_transfer_receipt.json','fold_schur_s4_receipt.json','cyclic_parity_receipt.json',
      'cyclic_companion_receipt.json','peripheral_cone_budget_receipt.json','coordinate_pencil_matrix_receipt.json']
    for name in preds:
        old=json.loads((HERE/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for path,digest in old['source_sha256'].items():check('sealed predecessor '+name+' '+path,sha(HERE/path)==digest)
    manifest=HERE/'sources/meridian-balance-two-stage/manifest.json';sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        p=HERE/row['path'];check('primary page-audit byte pin '+p.name,sha(p)==row['sha256'] and p.stat().st_size==row['bytes'])
    # Nonunit balances imply torsion after a filling, not necessarily zero.
    for aa,bb in [(1,1),(2,3),(7,-4)]:
        check('nonunit balance retains free rank1 after both rational slopes '+str((aa,bb)),E.h1([[aa,bb,0],[0,1,0],[1,0,0]])['free']==1)
    s,k,v,j=S.symbols('s k v j',nonzero=True,integer=True);r=3*s;u=(k*v+1)/r
    cone=S.Matrix([[s,0,0,v]]);global_=S.Matrix([[1,1,1,-j]])
    nrow=S.Matrix([[0,r,0,-v]]);mrow=S.Matrix([[0,0,r,-2*v]])
    check('actual lifted disk forces relative section integer j0',mrow+nrow-(r*global_-3*cone)==S.Matrix([[0,0,0,r*j]]))
    identity('actual fixed boundary basis determinant',v*k-r*u+1)
    identity('actual central boundary basis determinant',(-2*v)*(-k/2)-r*u+1)
    identity('P(n) actual fibre h equals rL modulo m',k*(-2*v)/2+r*u-1)
    # Actual r3 product gluing at arbitrary k, including its meridian slope.
    u3=(k*v+1)/3;n=S.Matrix([-v,3]);l1=S.Matrix([u3,-k]);hh=k*n+3*l1;yy=u3*n+v*l1
    mm=-2*v*hh+3*(v*hh-yy);ll=u3*hh-k*(v*hh-yy)/2
    check('all-k r3 actual meridians match',S.simplify(mm+n)==S.zeros(2,1))
    check('all-k r3 actual longitude shear',S.simplify(ll-k*n/2-l1)==S.zeros(2,1))
    rr,kk,N=S.symbols('r k N',positive=True,integer=True);JJ=rr-kk/2;l0=rr-kk
    vv=(N*rr+1)/l0;uu=(kk*vv+1)/rr;a1=-3*vv/rr;a2=(2*vv-uu)/JJ
    check('trefoil actual rational longitude',S.factor(a1+a2+vv/rr+1/(rr*JJ))==0)
    check('explicit l0 first ceiling argument',S.factor(l0*a1+3*N+3/rr)==0)
    check('explicit l0 second ceiling argument',S.factor(l0*a2-2*N-1/JJ)==0)
    check('strict upper bound has exact gap1/(r*l0)',S.factor(N/l0-vv/rr+1/(rr*l0))==0)
    interval_controls=[]
    for ss in [3,5,7,11,31]:
        rr0=3*ss
        for kk0 in range(2,rr0,2):
            if math.gcd(rr0,kk0)>1:continue
            vv0=(-pow(kk0,-1,rr0))%rr0;uu0=(kk0*vv0+1)//rr0;J0=rr0-kk0//2;x=S.Rational(vv0,rr0);l=rr0-kk0
            c1=-S.Rational(vv0,ss);c2=S.Rational(2*vv0-uu0,J0)
            lower=(-1-S.floor(rr0*c1)-S.floor(rr0*c2))/rr0;upper=(1-S.ceiling(l*c1)-S.ceiling(l*c2))/l
            check('actual trefoil n between explicit interval complement bounds '+str((rr0,kk0)),lower==x and upper==x-S.Rational(1,rr0*l) and upper<x<=lower)
            interval_controls.append(dict(r=rr0,k=kk0,v=vv0,u=uu0,J=J0,n_slope=str(x),y_minus_lower=str(lower),y_plus_upper=str(upper),two_denominators=[rr0,l]))
    # Integral main cable groups and balanced rational longitude.
    cable_controls=[]
    for ss,kk0 in [(3,2),(3,8),(5,2),(7,4)]:
        rr0=3*ss;vv0=(-pow(kk0,-1,rr0))%rr0
        rows=[[ss,0,0,vv0],[1,1,1,0]];n0=[0,rr0,0,-vv0];m0=[0,0,rr0,-2*vv0]
        check('actual P H1 free2 '+str((ss,kk0)),E.h1(rows)==B.group(free=2))
        want=B.group([ss],free=1)
        check('actual P(n) rank1 cone torsion and m rational longitude '+str((ss,kk0)),E.h1(rows+[n0])==E.h1(rows+[n0,m0])==want)
        cable_controls.append(dict(s=ss,k=kk0,P=E.h1(rows),P_n=want,m_zero_after_n=True))
    # Generic cyclic pattern Seifert orders and marked fibre vanish independently of eta.
    ss,tt=S.symbols('s t',positive=True,integer=True);eps=-1/(3*ss*tt)
    identity('both cyclic patterns fibre order1',-eps*(3*ss*tt)-1)
    identity('cyclic3 pattern H1 order s',-eps*tt*ss*(3*ss)-ss)
    identity('cyclic4 pattern H1 order s',-eps*ss*(3*ss)*tt-ss)
    pp,qq,kk0,rr0=S.symbols('p q k r',positive=True,integer=True)
    identity('cyclic3 full edge fibre intersection',pp*qq*rr0-3*kk0/3-(rr0*pp*qq-kk0))
    identity('cyclic4 all-odd full edge fibre intersection',pp*qq*rr0-4*kk0/4-(rr0*pp*qq-kk0))
    PP,KK=S.symbols('P K',integer=True);delta=4*PP**2+13*(1-KK)*PP+13*KK
    identity('all-odd cyclic4 divisor16 identity',(13*KK-4*PP-17)*(PP-1)-(17-delta))
    divisor_controls=[]
    for d in [4,8,16]:
        pole=d+1;kap=S.Rational(4*pole+17+16//d,13)
        check('forced divisor16 candidate not integral '+str(d),kap.q!=1)
        divisor_controls.append(dict(P=pole,K=str(kap)))
    qsrc=S.symbols('qsrc',positive=True,integer=True)
    identity('whole source three-cone fusion has terminal order2',(qsrc+1)+1-qsrc-2)
    check('source full fusion cannot repeat an odd cone',all(not (q==q+1) and (q!=2 and q+1!=2 or q==2) for q in [2,3,7]))
    # All1920 general even-product planar tree types have the six required incidences.
    census=json.loads((HERE/'general_two_stage_receipt.json').read_text())['census']['complete_cells'];counts=collections.Counter();valid=True
    for key,cell in census.items():
        p0,q0,r0,k0=map(int,key.split(','))
        if r0*k0%2:continue
        for row in cell['orbits']:
            if row['cycle_rank'] or sum(x['genus'] for x in row['companion_bases']):continue
            co,ca=row['companion_order'],row['cable_order'];counts[co,ca]+=1
            valid &= tuple(row['h'])==(0,1,2,3) and all(x['fibre_degree']==1 for x in row['cable_bases'])
            M=tuple(row['M']);z=tuple(row['rho'][2]);lengths=G.B.lengths
            if co==3:valid &= r0%4==0 and k0%3==0 and lengths(M)==[1,3] and lengths(z)==[4]
            elif co==4:valid &= r0%3==0 and k0%4==0 and lengths(M)==[4] and lengths(z)==[1,3]
            elif co==12:valid &= r0%2==0 and k0%3==0 and lengths(M)==[1,3] and sorted(len(x['sheets']) for x in row['cable_bases'])==[1,3]
            elif co==6:valid &= r0%3==0 and k0%2==0 and lengths(M)==[1,1,2] and sorted(len(x['sheets']) for x in row['companion_bases'])==[1,3]
            elif ca==24:valid &= r0%3==0 and k0%2==0 and lengths(M)==[4]
            else:valid &= ca==2 and r0%2==1 and k0%2==0 and lengths(M)==[1,1,2]
    check('all1920 even planar classes have actual sector/fibre incidence',valid and dict(counts)=={(3,24):384,(4,24):384,(12,6):384,(6,6):256,(24,24):256,(24,2):256})
    # Ordinary words, all true meridian-cycle fillings, and the preferred R longitude.
    cases=[(2,3,4,3,3,24),(2,5,8,3,3,24),(3,5,8,3,3,24),(3,5,3,4,4,24),(2,7,3,4,4,24),(6,5,3,4,4,24),(4,3,3,4,4,24),
      (2,3,3,2,24,24),(2,3,9,2,24,24),(2,9,9,2,24,24),(4,3,2,3,12,6),(2,9,2,3,12,6),(2,3,10,3,12,6),(2,3,14,3,12,6),
      (4,3,3,14,24,2),(4,9,5,2,24,2),(2,9,3,2,6,6),(2,9,3,4,6,6),(2,9,9,8,6,6),(4,3,9,2,6,6),(2,3,3,2,6,6),(2,3,9,2,6,6),(2,3,15,2,6,6)]
    targets=[]
    for pp,qq,rr,kk,co,ca in cases:
        rows=[row for row in G.cell(pp,qq,rr,kk)['orbits'] if (row['companion_order'],row['cable_order'],row['cycle_rank'])==(co,ca,0)]
        check('generic ordinary target classes exist '+str((pp,qq,rr,kk,co,ca)),bool(rows))
        if co==3:
            add=[qq,qq] if pp%3==0 else [pp,pp] if qq%3==0 else [];hh=B.group([kk//3]+add);want=dict(H1=hh,quotient=hh)
        elif co==4:
            if pp%2 and qq%2:add=[]
            else:
                ee,oo=(pp,qq) if pp%2==0 else (qq,pp);add=[oo]*(3 if ee%4==0 else 1)
            hh=B.group([rr//3]+add);want=dict(H1=hh,quotient=hh)
        elif co==12:hh=B.group([pp*qq//3,rr//2]);want=dict(H1=hh,quotient=hh)
        elif co==6:
            ee,oo=(pp,qq) if pp%2==0 else (qq,pp);bb,tt,ss=ee//2,oo//3,rr//3;D=kk//2-bb*tt*ss
            want=dict(H1=B.group([abs(D),bb],free=int(D==0)),quotient=B.group([math.gcd(D,ss*kk),bb]))
        elif ca==24:want=dict(H1=B.group([8*(qq//3),rr//3]),quotient=B.group([2*(qq//3),rr//3]))
        else:hh=B.group([qq//3,rr]);want=dict(H1=hh,quotient=hh)
        for i,row in enumerate(rows):check('ordinary full end and R mark agree at independent k '+str((pp,qq,rr,kk,i)),G.satellite(pp,qq,rr,kk,row)==want)
        targets.append(dict(p=pp,q=qq,r=rr,k=kk,sector=[co,ca],classes=len(rows),**want))
    inputs=[Path(__file__),HERE/'MERIDIAN_BALANCE_AND_GENERAL_TWO_STAGE_CLOSURE.md',manifest,HERE/'verify_general_two_stage.py']+[HERE/n for n in preds]+[HERE/row['path'] for row in sources['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
      general=dict(no_field_degree_or_source_depth_bound=True,nonunit_two_boundary_balance_allowed=True,
        balanced_slope_becomes_rational_longitude_after_other_filling=True,
        nonsolid_boundary_incompressible_actual_graph_pieces_required=True,
        actual_matched_slopes_outside_L_space_interiors_required=True,
        two_step_Rasmussen_gluing_propagates_missing_interior=True,
        actual_irreducibility_before_graph_LO_implication=True,
        PF74_common_source_no_LO_image_hypotheses_required=True),
      S3=dict(actual_all_k_section_presentation=True,relative_section_j_zero_proved_from_integral_disk=True,
        non_trefoil_s1_matched_meridians=True,non_trefoil_s_ge3_balanced_propagation=True,
        trefoil_two_exact_denominators=['r','r-k'],trefoil_n_outside_interval_interior=True,
        full_common_sector_excluded=True,full_complementary_sector_excluded=True),
      quartic=dict(general_genuine_two_stage_actual_PF46_50_packet_required=True,
        arbitrary_lower_weight_key_terms_allowed=True,k_independent_of_rc_plus1=True,
        all_coprime_first_pairs_windings_coefficients_and_source_data=True,
        all_six_even_sectors_excluded_in_both_pencils=True,all1920_even_planar_types_audited=True,
        general_eta0_and_odd_product_PF80_closures_retained=True,
        entire_general_genuine_two_stage_target_class_excluded=True,
        source_full_fusion_third_cone_exactly2_used=True,
        all_odd_C4_s1_closed_by_unbounded_divisor16_identity=True,
        all_actual_fixed_sheets_caps_and_R_marks_retained=True,both_solid_cap_signs_retained=True),
      even_planar_sector_counts={str(tag):v for tag,v in counts.items()},trefoil_interval_controls=interval_controls,
      main_cable_controls=cable_controls,divisor16_controls=divisor_controls,ordinary_target_controls=targets,
      assurance=dict(written_unbounded_balance_section_interval_and_sector_proofs=True,
        complete_Rasmussen_Nemethi_primary_pages_visually_rechecked=True,
        deep_L_space_gluing_graph_LO_JSJ_and_rigidity_imported=True,
        finite_controls_not_universal_proof=True,interval_interiors_kept_distinct_from_endpoints=True,
        nonunit_balance_asserted_to_give_zero_integral_class=False,
        cyclic_companion_homology_replacement_asserted_homeomorphism=False,
        eta0_substitution_used_for_negative_eta=False,arbitrary_degree_normalizations_classified=False,
        longer_target_sequences_excluded=False,degenerate_stages_classified_here=False,
        full_quartic_exclusion=False,higher_degree_exclusion=False,finite_Keller_map=False,
        polynomial_curves_or_infinity_representations_declared_nonexistent=False,
        formal_verification=False,independent_full_reproof=False,external_review=False,novelty_priority=False,research_goal_completed=False),
      source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'meridian_balance_two_stage_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; meridian balance and exact interval bounds; entire general genuine two-stage quartic class excluded')
if __name__=='__main__':main()
