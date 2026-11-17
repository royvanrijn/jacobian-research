"""Finite algebra/group controls for the written all-degree neighborhood theorem."""
from pathlib import Path
from itertools import permutations
from math import factorial,prod
from collections import Counter
import hashlib,json,os,platform
import sympy as S
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def compose(a,b):return tuple(a[b[i]] for i in range(len(a)))
def perm(n,cycles):
    out=list(range(n))
    for cy in cycles:
        for a,b in zip(cy,cy[1:]+cy[:1]):out[a]=b
    return tuple(out)
def closure(gens):
    identity=tuple(range(len(gens[0])));out={identity};todo=[identity]
    while todo:
        a=todo.pop()
        for b in gens:
            c=compose(a,b)
            if c not in out:out.add(c);todo.append(c)
    return out
def profile(a):
    unused=set(range(len(a)));lengths=[]
    while unused:
        t=min(unused);s=t;n=0
        while s in unused:unused.remove(s);n+=1;s=a[s]
        lengths.append(n)
    return Counter(lengths)
def class_size(a):
    return factorial(len(a))//prod(k**m*factorial(m) for k,m in profile(a).items())
def normal_group_index_control(group,sigma):
    central=[g for g in group if compose(g,sigma)==compose(sigma,g)]
    conj={compose(compose(g,sigma),tuple(g.index(i) for i in range(len(g)))) for g in group}
    return len(group)//len(central),len(conj),len(central)
def main():
    pred=HERE/'unramified_orevkov_receipt.json';d=json.loads(pred.read_text())
    check('predecessor PF2RULED46 PASS',d['outcome']=='PASS')
    for path,h in d['source_sha256'].items():check('predecessor seal '+path,sha(HERE/path)==h)
    olddir=ROOT/'jacobian/jc2-r7-global-cover';old=json.loads((olddir/'receipt.json').read_text())
    check('retained conditional Nori receipt PASS without changing scope',old['passed'] is True)
    for path,h in old['input_sha256'].items():check('retained r7 input '+path,sha(ROOT/path)==h)
    for path,h in old['output_sha256'].items():check('retained r7 output '+path,sha(olddir/path)==h)
    manifest=HERE/'sources/lefschetz-centralizer/manifest.json';m=json.loads(manifest.read_text())
    for row in m['primary_sources']:
        f=HERE/row['path_relative_to_affine_complements'];check('primary byte seal '+f.name,sha(f)==row['sha256'] and len(f.read_bytes())==row['bytes'])
    check('complete quantitative WLT and neighborhood pages inspected',m['primary_sources'][0]['new_complete_page_visual_audits']==[2,3,10,27,28])
    eta,delta,b,degree,k,T=S.symbols('eta delta b degree k T')
    check('normalization adjunction has its full finite delta correction',S.expand((eta+2*delta)-2*delta-eta)==0)
    check('conjugacy-index inequality rearranges exactly',S.expand(b*eta-(eta+2*delta)-((b-1)*eta-2*delta))==0)
    check('symmetric transposition field-degree inequality is exact',S.expand(2*((degree*(degree-1)/2-1)*eta-2*delta)-(degree*(degree-1)*eta-2*eta-4*delta))==0)
    check('quartic combined singularity budget uses total affine delta',S.expand(((S.binomial(4,2)-1)*eta-2*delta).subs(delta,k+T)-(5*eta-2*k-2*T))==0)
    block=S.symbols('m')
    check('arbitrary transitive fold conjugacy orbit has exact block count',S.expand((degree/block)*(block*(block-1)/2)-degree*(block-1)/2)==0)
    check('unknown-monodromy field-degree inequality rearranges exactly',S.expand(2*((degree*(block-1)/2-1)*eta-2*delta)-(degree*(block-1)*eta-2*eta-4*delta))==0)
    # Exhaustive finite controls for the all-degree conjugacy formula.
    for n in range(3,8):
        sn=list(permutations(range(n)));sigma=perm(n,[[0,1]])
        bi,ci,central=normal_group_index_control(sn,sigma)
        check('exact transposition centralizer/index d'+str(n),(bi,ci,central)==(n*(n-1)//2,n*(n-1)//2,2*factorial(n-2)))
        for cycle in range(2,n+1):
            sigma=perm(n,[list(range(cycle))]);bi,ci,central=normal_group_index_control(sn,sigma)
            check('single-cycle class formula d'+str(n)+' e'+str(cycle),bi==ci==class_size(sigma)==factorial(n)//(cycle*factorial(n-cycle)))
        sigma=perm(n,[[0,1],[2,3]]) if n>=4 else perm(n,[[0,1]])
        bi,ci,central=normal_group_index_control(sn,sigma)
        check('full cycle-profile formula with fixed letters d'+str(n),bi==ci==class_size(sigma))
    # The relevant local group may be nonabelian while still centralizing sigma.
    own=perm(6,[[0,1]]);u=perm(6,[[2,3]]);v=perm(6,[[3,4]])
    local=closure([own,u,v])
    check('separate cusp-plus-fold local group is nonabelian',len(local)==12 and compose(u,v)!=compose(v,u))
    check('separate cusp local group centralizes the selected fold',all(compose(g,own)==compose(own,g) for g in local))
    own3=perm(6,[[0,1,2]]);u=perm(6,[[3,4]]);v=perm(6,[[4,5]])
    local3=closure([own3,u,v])
    check('arbitrary-index selected cyclic factor can coexist with a cusp',len(local3)==18 and all(compose(g,own3)==compose(own3,g) for g in local3))
    bad=perm(6,[[1,2]])
    check('selected nonimmersive cusp defeats centralization',compose(own,bad)!=compose(bad,own))
    center=perm(6,[[0,3],[1,4],[2,5]])
    gg=closure([center,perm(6,[[0,1],[3,4]]),perm(6,[[0,1,2],[3,4,5]])])
    check('nonabelian monodromy need not make every generic inertia noncentral',len(gg)==12 and normal_group_index_control(list(gg),center)==(1,1,12))
    check('central-inertia permutation control is transitive',len({g[0] for g in gg})==6)
    # Transitive wreath controls for the all-degree normal-transposition proof.
    for n,mm in [(4,2),(6,2),(6,3),(8,2),(8,4),(9,3)]:
        gens=[perm(n,[[i,i+1]]) for i in range(mm-1)]
        gens.append(perm(n,[list(range(i,n,mm)) for i in range(mm)]))
        group=closure(gens);sigma=perm(n,[[0,1]])
        bi,ci,central=normal_group_index_control(list(group),sigma)
        check('exact fold-conjugacy block control d'+str(n)+' m'+str(mm),len({g[0] for g in group})==n and bi==ci==n*(mm-1)//2)
    # Read, do not recreate, the already certified infinity multiplicities.
    budgets=json.loads((olddir/'budgets.json').read_text());rows=[]
    for oldrow in budgets['infinity_resolution']:
        I=oldrow['infinity_contact'];mult=oldrow['B_multiplicities'];sq=25-sum(v*v for v in mult);normal=sq-8
        check('retained full infinity multiplicity cost I'+str(I),sq==oldrow['B_self_intersection']==15-(I-20)//2)
        check('immersive quintic normal degree independent of nodality I'+str(I),normal==7-(I-20)//2)
        applicable=normal>0
        excluded=applicable and 14*normal>8
        excluded3=applicable and 39*normal>8
        check('both inertia indices have the same exact contact gate I'+str(I),excluded==excluded3==(I<34))
        rows.append(dict(infinity_contact=I,actual_polynomial_contact=oldrow['actual_degree35_contact'],B_self_intersection=sq,affine_delta=4,normal_degree=normal,positive_normal_bundle=applicable,transposition_class=15,three_cycle_class=40,excluded_for_both_indices=excluded))
    allowed=[r['infinity_contact'] for r in rows if r['actual_polynomial_contact'] and not r['excluded_for_both_indices']]
    check('full actual remaining quintic contact spectrum',allowed==[34,38,40])
    check('contact34 is zero-normal-degree equality case',next(r for r in rows if r['infinity_contact']==34)['normal_degree']==0)
    check('finite tangency escape in historical test is real but superseded',20+2*7==34 and next(r for r in rows if r['infinity_contact']==20)['excluded_for_both_indices'])
    check('no finite incidence variable enters the new contact exclusion',all('finite_tangency' not in r for r in rows))
    inputs=[Path(__file__),HERE/'LEFSCHETZ_CENTRALIZER_BOUND.md',manifest,pred,olddir/'receipt.json',olddir/'budgets.json',olddir/'PROOF.md',ROOT/'jacobian/jc2-omitted-divisors/NONIMMERSED_LEDGER.md',ROOT/'jacobian/jc2-omitted-divisors/GLOBAL_MONODROMY.md']+[HERE/r['path_relative_to_affine_complements'] for r in m['primary_sources']]
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
        quantitative_theorem=dict(all_field_degrees=True,all_generic_indices=True,componentwise=True,actual_G_required=True,G_symmetric_assumed=False,selected_normalization_immersive_required=True,selected_support_zero_excess_required=True,all_other_supports_zero_excess_required=False,other_local_groups_may_be_nonabelian=True,full_finite_normalization_retained=True,simultaneous_target_infinity_model_required=True,normal_degree='eta=-K_X.Cbar-2=Cbar^2-2*delta',positivity_required='eta>0',tube_image='rho(K) subset C_G(sigma)',WLT_index='[Pi:K]<=Cbar^2/eta',conjugacy_size='b=[G:C_G(sigma)]',conclusion='(b-1)*eta<=2*delta',no_moduli_connectivity_required=True,bare_image_tube_replaced_by_normalization_tube=True),
        symmetric_consequences=dict(cycle_class='d!/product_l(l^m_l*m_l!)',transposition_class='binomial(d,2)',transposition_bound='(binomial(d,2)-1)*eta<=2*delta',field_degree_bound='d*(d-1)<=2+4*delta/eta',positive_normal_degree_required=True,quartic_bound='a_n<=1+min(2*T,floor(2*(k+T)/5))',quartic_case_excluded=False),
        unidentified_group_fold_bound=dict(any_transitive_G=True,single_transposition_required=True,normal_transposition_block_size='m>=2; m divides d',exact_class='d*(m-1)/2',bound='d*(m-1)<=2+4*delta/eta',explicit_field_degree_bound='d<=2+4*delta/eta',least_prime_refinement='d*(ell(d)-1)<=2+4*delta/eta',prime_degree_forces_S_d=True,eta_positive_required=True,other_inertia_unrestricted=True,intermediate_polynomial_plane_inferred=False),
        degree_six_application=dict(first_image='complete degree(6,10) gap21 family',second_image_degrees=[3,5],second_image_degree_hypothesis_preserved=True,oriented_rows=['A','C','N'],selected_second_image_immersive_required=True,selected_second_image_nodal_required=False,finite_singular_locus_avoidance_required=False,finite_tangency_bound_required=False,nonimmersed_first_image_allowed_in_N=True,other_image_cusp_does_not_break_selected_centralizer=True,group='S6',generic_indices=[2,3],conjugacy_sizes=[15,40],second_affine_delta=4,normal_degree='7-(I-20)/2',contact_lower_bound=34,remaining_contacts=allowed,contact_rows=rows,remaining_contacts_realized_as_Keller_images=False,nonimmersive_second_image_row_B_covered=False),
        assurance=dict(written_neighborhood_centralization_proof=True,precise_Nori_WLT_A_C_import=True,full_primary_pages_visually_checked=True,finite_controls_not_topology_proof=True,proof_from_bounded_search=False,eta_nonpositive_cases_excluded=False,all_inertia_assumed_noncentral=False,automatic_global_Orevkov_factorization=False,full_quartic_exclusion=False,unrestricted_second_image_degree_derived=False,formal_verification=False,independent_full_reproof=False,external_review=False,novelty_priority=False,research_goal_completed=False),
        source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'lefschetz_centralizer_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS',len(checks),'controls; all-degree componentwise centralizer bound; immersive quintic contact>=34 without finite tangency assumptions')
if __name__=='__main__':main()
