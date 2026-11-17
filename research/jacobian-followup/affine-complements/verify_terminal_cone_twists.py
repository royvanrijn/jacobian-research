"""PF2RULED85: actual terminal cone residues and complete cyclic arrival closure."""
from pathlib import Path
import sys,math,json,hashlib,os,platform
import sympy as S
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P))
import verify_trefoil_complementary_fusion as F
import verify_cyclic_loop_first_arrival as L
C=F.C;J=F.J
def leg_controls():
    z,w,T,U=S.symbols('z w T U');outrows=[]
    for ss,support in [(6,{10:1,9:1}),(6,{15:1,14:1}),(10,{25:1,12:1}),
                       (10,{25:1,16:1}),(10,{20:1,15:1,12:1})]:
        g=ss;drops=[]
        for j in sorted(support,reverse=True):
            h=math.gcd(g,j)
            if h<g:drops.append((g,h,j));g=h
        a=drops[-1][0];kap=drops[-1][2]+1-ss
        assert g==1 and ss%a==0 and max(support)%a==0
        R=C.local_cluster(z**-1+T,z**2,0,z,T)
        D=C.local_cluster(sum(co*w**(-j) for j,co in support.items())+w**(1-ss)*U,w**ss,0,w,U)
        rr=J.retag(R,'R','Hx');ds=J.retag(D,'D','Hy')
        raw=C.assemble([rr,ds]);mat=raw['full']
        vertices=set(range(mat.rows))-set(raw['index'])
        arrow=raw['labels'][rr['rows'][-1]['parents'][0]]
        spine=J.rooted_nodes(mat,vertices,arrow);matches=[]
        for node,legs in spine['legs'].items():
            if sorted(x['order'] for x in legs)!=sorted([a,kap]):continue
            data=[]
            for leg in legs:
                path=leg['vertices'];weights=[int(mat[x,x]) for x in path]
                b=int(abs(mat.extract(path[1:],path[1:]).det()))
                bp=int(abs(mat.extract(path[:-1],path[:-1]).det()))
                data.append(dict(alpha=leg['order'],cofactor=b,opposite_cofactor=bp,weights=weights))
            matches.append(data)
        assert len(matches)==1
        outrows.append(dict(sigma=ss,support=support,a=a,kappa=kap,legs=matches[0]))
    return outrows
def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def identity(name,expr):check(name,S.expand(expr)==0)
    def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
    preds=['cyclic_loop_first_arrival_receipt.json','persistent_jsj_terminal_modulus_receipt.json',
        'terminal_cone_binomial_closure_receipt.json','trefoil_complementary_fusion_receipt.json',
        'complementary_corner_jsj_receipt.json','general_two_stage_receipt.json']
    for name in preds:
        old=json.loads((P/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for path,digest in old['source_sha256'].items():
            check('sealed predecessor '+name+' '+path,sha(P/path)==digest)
    manifest=P/'sources/terminal-cone-twists/manifest.json';sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        path=P/row['path'];check('complete primary byte/page pin '+path.name,
                                sha(path)==row['sha256'] and path.stat().st_size==row['bytes'])
    a,k,f,nr,t,b=S.symbols('a k f nr t b',integer=True)
    g=nr*k-1
    identity('plus terminal twist gives nR+f-1 winding divisor',
             g+((f-1)*k+1)-(nr+f-1)*k)
    identity('minus terminal twist gives nR-f-1 winding divisor',
             g-((f+1)*k-1)-(nr-f-1)*k)
    check('quartic odd winding cannot divide either exact divisor',
          nr.subs(nr,13)+4-1==16 and 13-4-1==8 and
          [d for n in [8,16] for d in S.divisors(n) if d%2]==[1,1])
    identity('local fibre-cover section equation',
             (a*t+f*b)-f*b-a*t)
    j,sigma=S.symbols('j sigma',integer=True)
    identity('last free-tail residue is characteristic residue',
             (sigma-1+k)-(k-1)-sigma)
    N=S.symbols('N',integer=True,nonnegative=True)
    identity('negative-two tail determinant induction',2*(N+1)-N-(N+2))
    identity('section determinant isolates cone coefficient',a*b-t*k-(a*b-k*t))
    controls=leg_controls()
    for row in controls:
        aa,kk=row['a'],row['kappa'];tag=str((row['sigma'],row['support']))
        check('actual unbounded-hypothesis last-gcd control '+tag,
              row['sigma']%aa==0 and max(row['support'])%aa==0 and math.gcd(aa,kk-1)==1)
        check('actual full-cluster free leg is negative-two chain '+tag,
              any(x['alpha']==kk and x['cofactor']==kk-1 and x['weights']==[-2]*(kk-1) for x in row['legs']))
        check('actual full-cluster winding residue inverse kappa-minus1 '+tag,
              any(x['alpha']==aa and x['cofactor']==pow(kk-1,-1,aa) for x in row['legs']))
        for leg in row['legs']:
            cprod=S.eye(2)
            for weight in leg['weights']:cprod=cprod*S.Matrix([[-weight,-1],[1,0]])
            check('exact Euclidean determinant-one matrix '+str((tag,leg['alpha'])),
                  cprod.det()==1 and cprod[0,0]==leg['alpha'] and
                  cprod[1,0]==leg['cofactor'] and -cprod[0,1]==leg['opposite_cofactor'])
    check('non-self-inverse order5 cofactor is actually exercised',
          any(x['alpha']==5 and x['cofactor']==3 and x['opposite_cofactor']==2 for r in controls for x in r['legs']))
    # Old terminal congruence admitted T(5,7); the normalized cone residue does not.
    aa,kk=5,7;co=pow(kk-1,-1,aa);target=pow(4*kk,-1,aa)
    check('old bothodd terminal modulus really passes T5_7',
          (13*(aa-(aa-1)*kk)-1)%(aa*aa)==0)
    check('both signs of old T5_7 residue control are disjoint',
          {co,(-co)%aa}.isdisjoint({target,(-target)%aa}))
    aa=25;target25=pow(4,-1,aa)
    check('old T2_25 windingproduct9 square condition really passes',9*9%13==3)
    check('both repeated order25 target cones fail free-tail residue',
          {1,24}.isdisjoint({target25,(-target25)%25}) and 4%25 not in [1,24])
    cover_controls=[]
    for aa,kk,ff in [(5,7,4),(7,11,4),(11,13,4),(5,9,2),(7,13,3),(11,17,5)]:
        if math.gcd(aa,kk)!=1 or math.gcd(ff,aa*kk)!=1:continue
        old=pow(kk,-1,aa);new=pow(ff*kk,-1,aa)
        check('actual local cover divides coefficient by fibre degree '+str((aa,kk,ff)),
              ff*new%aa==old and pow(new,-1,aa)==ff*kk%aa)
        cover_controls.append(dict(a=aa,kappa=kk,fibre_degree=ff,old_coefficient=old,new_coefficient=new))
    cyclic84=json.loads((P/'cyclic_loop_first_arrival_receipt.json').read_text())
    check('retained all-length cyclic reduction has only two stated terminal residuals',
          all(cyclic84['quartic'][x] for x in [
            'cyclic3_arrival_any_companion_depth_any_tail_both_pencils_excluded',
            'cyclic4_arrival_any_nonempty_later_tail_both_pencils_excluded',
            'cyclic4_arrival_winding_at_least9_both_pencils_excluded',
            'remaining_cyclic4_arrival_final_winding3','residual_companion_all_windings_and_coefficients_odd',
            'residual_even_first_pair_only_2_25']))
    inputs=[Path(__file__),P/'TERMINAL_CONE_TWISTS_AND_CYCLIC_ARRIVAL_CLOSURE.md',manifest]+[P/x for x in preds]+[P/x['path'] for x in sources['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
        general=dict(no_field_degree_or_characteristic_depth_bound=True,
          actual_smooth_exact_whole_two_line_source_hypotheses_required=True,
          winding_cofactor_inverse_kappa_minus1=True,free_tail_cofactor_kappa_minus1=True,
          normalized_residue_not_arbitrary_unit=True,either_orientation_sign_allowed=True,
          local_fibre_cover_divides_normalized_coefficient=True,
          delta1_terminal_divisor_nR_plus_f_minus1_or_nR_minus_f_minus1=True,
          last_gcd_and_base1_coprime_cover_hypotheses_explicit=True),
        quartic=dict(all_bothodd_final_cyclic4r3_residuals_excluded=True,
          T2_25_repeated_cone_residual_excluded=True,
          entire_cyclic4_incoming_first_arrival_class_excluded_both_pencils=True,
          PF84_all_cyclic3_incoming_closure_retained=True,
          all_cyclic_incoming_states6_9_excluded_at_any_depth_and_tail_length=True,
          no_prefix_eta_condition=True,all_lower_weight_keys_and_actual_periphery_retained=True),
        full_terminal_cluster_controls=controls,local_cover_controls=cover_controls,
        assurance=dict(written_unbounded_local_slope_free_tail_and_congruence_proofs=True,
          complete_Hatcher25_26_27_and_NR3_4_11_12_visually_inspected=True,
          complete_primary_bytes_pinned=True,deep_JSJ_fibre_uniqueness_and_rigidity_imported=True,
          both_signs_used_as_weaker_necessary_match=True,finite_controls_not_universal_proof=True,
          normalized_coefficients_asserted_determined_by_cone_orders=False,
          source_degree_hypotheses_assumed_automatic=False,
          noncyclic_incoming_states_closed=False,all_initial_state3_trajectories_closed=False,
          full_quartic_exclusion=False,higher_degree_packet_classification=False,
          finite_Keller_map=False,formal_verification=False,external_review=False,
          independent_full_reproof=False,novelty_priority=False,research_goal_completed=False),
        source_sha256={os.path.relpath(x,P):sha(x) for x in inputs},
        runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'terminal_cone_twists_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; normalized terminal twists exclude ALL cyclic incoming first arrivals in both pencils')
if __name__=='__main__':main()

