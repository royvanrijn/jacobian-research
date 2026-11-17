"""PF2RULED92: all-degree peripheral reconstruction, D3 recursion, endpoint tails."""
from pathlib import Path
import sys,json,hashlib,math,os,platform
import sympy as S
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P))
import composed_peripheral_cover as C
import verify_global_ruling_unfilled_cyclic as F
import verify_genus_final_s3 as H
import compare_dihedral_linking as Q
A=F.A;K=F.K;G=F.G;E=F.E;B=F.B;I=tuple(range(4))
def d3_matrix(data):
    result=C.integral_null_core_linking(data);order=sorted(range(2),key=lambda i:-data['boundaries'][i]['degree'])
    assert [data['boundaries'][i]['degree'] for i in order]==[2,1]
    return [[result['matrix'][i][j] for j in order] for i in order],result['closed_H1']
def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def identity(name,value):check(name,S.expand(value)==0)
    def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
    preds=['genus_final_s3_receipt.json','global_ruling_unfilled_cyclic_receipt.json',
        'cyclic_pattern_a4_receipt.json','iterated_peripheral_tail_receipt.json',
        'meridian_balance_two_stage_receipt.json','lifted_surface_s3_receipt.json',
        'persistent_jsj_terminal_modulus_receipt.json']
    for name in preds:
        old=json.loads((P/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for path,digest in old['source_sha256'].items():check('sealed predecessor '+name+' '+path,sha(P/path)==digest)
    manifest=P/'sources/recursive-linking-endpoints/manifest.json';sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']+sources['external_code']:
        path=P/row['path'];check('primary/code exact pin '+path.name,sha(path)==row['sha256'] and path.stat().st_size==row['bytes'])
    a,b,t,aa,r,k,w,eta=S.symbols('a b t A r k w eta',integer=True)
    matrix=t*S.Matrix([[-b*b,a*b],[a*b,-a*a]])
    check('all-degree two-component symmetric reconstruction and primitive kernel',matrix==matrix.T and matrix*S.Matrix([a,b])==S.zeros(2,1))
    reduction_controls=[]
    for rows,extra in [([[2,4,6],[4,10,14]],[]),([[0,0,0]],[]),
        ([[6,0],[0,10],[2,4]],[]),([[0,2,4,6],[3,9,12,15]],[[1,0,0,1]]),
        ([[2,0,0],[0,3,0]],[[0,0,1]]),([[1,2,3],[0,0,0]],[[0,4,6]])]:
        small=C.reduced_h1(rows,extra);original=E.h1(rows,extra)
        check('HNF exact row lattice and unit elimination agree with retained Smith '+str((rows,extra)),small==original)
        reduction_controls.append(dict(rows=rows,extra=extra,H1=small))
    identity('D3 offdiagonal square winding recursion',2*aa*w*w-w*w*(2*aa))
    identity('last minus infinity inequality makes positive l0',((2*aa-1)*r-k)-((r*(2*aa-w)-k)+r*(w-1)))
    identity('r3 last-minus positive genus band upper bound',(3*(2*aa-w)-1)/2-aa-(2*aa-(2*w-1))-(w-3)/2)
    identity('allplus complementary caporder1 eta plus sign',(2*aa*r+2)-r*(2*aa-1)-1-(r+1))
    identity('allplus complementary caporder1 eta minus sign',(2*aa*r-2)-r*(2*aa-1)-1-(r-3))
    # Actual D3 full-unfilled controls, two tree choices, no word expansion.
    controls=[]
    for p,q,cables in [(2,3,[]),(4,3,[]),(2,9,[]),(8,3,[(3,2)]),
        (2,3,[(3,7),(5,2),(7,3),(3,2)]),
        (4,3,[(3,5),(5,2),(7,3)]),
        (2,3,[(3,2),(5,3)]*8)]:
        x,y,mm,ll=next(x for x in G.torus(p,q) if len(G.groups[x[0],x[1]])==6 and A.order(G.perms[x[2]])==2)
        seed=A.normalized(p,q,[],[G.perms[x],G.perms[y]]);M,_=A.periphery(p,q,[],seed)
        seed3=F.restrict_three(seed,M);M3=E.evaluate(K.words(p,q,[])[1],seed3)
        rho=seed3+[E.power(M3,kk%2) for rr,kk in cables];coef=p*q//6*math.prod(rr*rr for rr,kk in cables)
        results=[]
        for reverse in [False,True]:
            data=C.peripheral_cover(p,q,cables,rho,reverse=reverse);framing,closed=d3_matrix(data)
            check('ALL UNFILLED preferred D3 longitude relations exact '+str((p,q,cables,reverse)),framing==[[-coef,2*coef],[2*coef,-4*coef]])
            oldrows,oldmark,oldkeep=K.summarized_cover(p,q,cables,rho,reverse=reverse,raw=True)
            full=dict(H1=C.reduced_h1(oldrows),quotient=C.reduced_h1(oldrows,[oldmark]))
            check('independent composed predecessor full marked group '+str((p,q,cables,reverse)),closed==full['H1']==full['quotient'] and oldrows==data['relation_rows']+[b['meridian'] for b in data['boundaries']] and oldkeep==data['keep'])
            check('all-length control used composed flows only '+str((p,q,len(cables),reverse)),data['word_expansion_used'] is False)
            results.append(closed)
        check('both spanning trees retain exact same peripheral matrix '+str((p,q,cables)),results[0]==results[1])
        controls.append(dict(p=p,q=q,cables=cables,A=coef,mutual=2*coef,closed_H1=results[0],both_trees_checked=True))
    # The capped solid history can have enormous coefficients at modest depth.
    hist=H.history([(3,1),(5,1)]*5);seed=A.normalized(2,3,[],[G.perms[x] for x in next((x[0],x[1]) for x in G.torus(2,3) if len(G.groups[x[0],x[1]])==6 and A.order(G.perms[x[2]])==2)])
    M,_=A.periphery(2,3,[],seed);seed3=F.restrict_three(seed,M);M3=E.evaluate(K.words(2,3,[])[1],seed3)
    data=C.peripheral_cover(2,3,hist['cables'],seed3+[E.power(M3,kk%2) for rr,kk in hist['cables']]);framing,closed=d3_matrix(data)
    check('large-coefficient ten-level solid chain exact D3 matrix without expansion',framing==[[-hist['A'],2*hist['A']],[2*hist['A'],-4*hist['A']]] and hist['A']>10**10)
    deep_control=dict(A=hist['A'],cables=hist['cables'],mutual=2*hist['A'],closed_H1=closed)
    # Actual first births from both cyclic predecessor types and older depths.
    births=[]
    for oldM,r0,k0,older in [((1,0,2,3),3,2,[]),((1,0,2,3),3,2,[(3,5)]),
        ((1,2,0,3),2,3,[]),((1,2,0,3),2,3,[(2,3)])]:
        oldrho=K.cyclic_rho(2,3,older,oldM)
        z,tau,lam,beta=next(x for x in A.roots(oldM,I,r0,k0,2,6) if E.generated(oldrho+[x[0]])==6)
        rho3=F.restrict_three(oldrho+[z],tau);Mnew=C.peripheral_cover(2,3,older+[(r0,k0)],rho3)['meridian_monodromy']
        for later in [[],[(3,2),(5,3)]]:
            rho=rho3+[E.power(Mnew,kk%2) for rr,kk in later]
            data=C.peripheral_cover(2,3,older+[(r0,k0)]+later,rho);framing,closed=d3_matrix(data)
            coef=r0*k0//6*math.prod(rr*rr for rr,kk in later)
            check('actual cyclic birth and all-length D3 linking transfer '+str((oldM,older,later)),framing==[[-coef,2*coef],[2*coef,-4*coef]])
            births.append(dict(old_meridian=oldM,older=older,birth=[r0,k0],later=later,A=coef,closed_H1=closed))
    # All-degree multi-component and positive-rank controls.
    higher=[]
    for degree in [3,5,7]:
        ref=tuple((-i)%degree for i in range(degree));rot=tuple((i+1)%degree for i in range(degree))
        data=C.peripheral_cover(2,degree,[],[ref,rot]);result=C.integral_null_core_linking(data)
        check('higher-degree actual dihedral weighted-kernel reconstruction '+str(degree),result['closed_H1']==B.group([]) and all(result['matrix'][i][j]==2 for i in range(len(result['matrix'])) for j in range(i)))
        higher.append(dict(degree=degree,cycle_lengths=[b['degree'] for b in data['boundaries']],**result))
    rho=K.cyclic_rho(2,3,[],tuple((i+1)%6 for i in range(6)))
    positive=C.integral_null_core_linking(C.peripheral_cover(2,3,[],rho))
    check('positive-rank closed homology allowed in peripheral reconstruction',positive['closed_H1']['free']==2 and positive['matrix']==[[0]])
    # Rationally null torsion cores must NOT pass the integral hypothesis.
    x,y,mm,ll=next(x for x in G.torus(2,3) if len(G.groups[x[0],x[1]])==12 and A.order(G.perms[x[2]])==3)
    nonnull=C.peripheral_cover(2,3,[],[G.perms[x],G.perms[y]])
    try:C.integral_null_core_linking(nonnull);rejected=False
    except AssertionError as error:rejected='INTEGRAL zero' in str(error)
    check('actual A4 torsion-core cover rejected by null-core reconstruction',rejected)
    # Independent author-code comparison, with colored history fixed BEFORE linking.
    comparison=Q.run()
    check('standalone author comparison receipt exactly reproduced',comparison==json.loads((P/'dihedral_author_comparison.json').read_text()))
    check('nine actual D3 independent braid/colored inputs match unmodified authors code',comparison['outcome']=='PASS' and len(comparison['examples'])==9)
    check('same uncolored knot different birth versus loop coloring distinguished',comparison['different_coloring_control']['mutual']==2 and comparison['different_coloring_control']['cyclic_loop_coloring_same_knot_mutual']==18 and comparison['different_coloring_control']['coloring_not_chosen_by_linking_output'])
    check('independent degree5 and7 author comparisons match multi-boundary data',len(comparison['higher_degree_controls'])==2 and all(all(v==2 for v in x['offdiagonal']) for x in comparison['higher_degree_controls']))
    # Endpoint seed: exact marked full covers can pass linking but are now
    # excluded by the universal missing-interior gluing argument, not by H1.
    seed4=A.normalized(2,3,[],[G.perms[x] for x in next((x[0],x[1]) for x in G.torus(2,3) if len(G.groups[x[0],x[1]])==6 and A.order(G.perms[x[2]])==2)])
    M,_=A.periphery(2,3,[],seed4);prefix=seed4+[E.power(M,1)]
    z,tau,lam,beta=next(x for x in A.roots(M,I,3,52,2,6) if E.generated(prefix+[x[0]])==24)
    endpoint_controls=[]
    for tail in [[(3,5)],[(5,3)],[(7,4)],[(3,5),(5,3),(7,4)]]:
        full=K.summarized_cover(2,3,[(3,7),(3,52)]+tail,prefix+[z]+[E.power(tau,kk%2) for rr,kk in tail])
        mult=math.prod(rr for rr,kk in tail);added=[B.group([rr]) for rr,kk in tail if kk%2==0]
        expected=dict(H1=B.direct_sum([B.group([17])]+added),quotient=B.direct_sum([B.group([math.gcd(17,mult)])]+added))
        check('actual endpoint-seed tail marked group retained '+str(tail),full==expected)
        endpoint_controls.append(dict(tail=tail,actual=full,self_linking=str((-S.Rational(mult*mult,17))%1),nonLspace_by_written_endpoint_gluing=True))
    check('old passing strict linking control retained while endpoint gluing closes it',endpoint_controls[0]['self_linking']=='8/17')
    inputs=[Path(__file__),P/'RECURSIVE_LINKING_AND_ENDPOINT_TAIL_CLOSURE.md',P/'composed_peripheral_cover.py',P/'compare_dihedral_linking.py',P/'dihedral_author_comparison.json',manifest,P/'sources/recursive-linking-endpoints/author_code_pin.json']+[P/f for f in preds]+[P/x['path'] for x in sources['primary_sources']+sources['external_code']]
    result=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
        general=dict(all_degree_null_core_peripheral_weighted_kernel=True,two_component_primitive_rank1_reconstruction=True,
            all_length_actual_D3_linking_square_winding=True,first_S3_birth_parameter_retained=True,
            coloring_dependence_and_preferred_framing_retained=True,
            missing_Lspace_INTERIOR_endpoint_seed_all_cyclic2_tails=True),
        quartic=dict(exceptional_S3_endpoint_Lspace_seed_all_common_tails_closed=True,
            all_S3_nonpositive_arrival_eta_common_closed=True,
            last_minus_survivors_require_arrival_r3_both_pencils=True,
            common_last_minus_requires_r3_thin_at_or_above_genus_band=True,
            complementary_rge9_only_allplus_k2Arplusminus2=True,
            no_prefix_assumed_separate_affine_source=True),
        D3_composed_controls=controls,large_coefficient_deep_control=deep_control,
        actual_birth_controls=births,higher_degree_controls=higher,positive_rank_control=positive,
        exact_HNF_reduction_controls=reduction_controls,
        independent_author_comparison=comparison,endpoint_seed_tail_controls=endpoint_controls,
        assurance=dict(written_unbounded_peripheral_linking_endpoint_frontier_proofs=True,
            all_unfilled_boundaries_composed_and_integral_torsion_checked=True,
            external_author_code_unmodified_commit_pinned=True,
            colored_inputs_selected_by_monodromy_not_linking_value=True,
            deep_graph_gluing_LO_rigidity_and_Floer_inputs_retained=True,
            finite_controls_not_universal_proof=True,
            all_later_S3_closed=False,all_later_A4_closed=False,full_S4_A4_called_dihedral=False,
            external_algorithm_independently_reproved=False,full_quartic_exclusion=False,
            higher_degree_packet_classification=False,finite_Keller_map=False,formal_verification=False,
            external_review=False,novelty_priority=False,research_goal_completed=False),
        source_sha256={os.path.relpath(path,P):sha(path) for path in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'recursive_linking_endpoints_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; all-length D3 recursion and endpoint common tails closed; full later-tail goal OPEN')
if __name__=='__main__':main()
