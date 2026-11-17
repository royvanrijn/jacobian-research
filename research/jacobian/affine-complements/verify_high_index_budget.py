"""PF112 exact controls supplement the written all-degree proof and classical imports."""
from pathlib import Path
import hashlib,json,math,platform
import sympy as S
import high_index_budget as H
import full_peripheral_linking as F
import whole_plane_peripheral as W
import normal_resolution_chern as N

P=Path(__file__).resolve().parent


def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def ident(name,value):check(name,S.factor(value)==0)
    def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
    for predecessor in ['whole_plane_peripheral_receipt.json','cyclic_quotient_morse_receipt.json','global_boundary_receipt.json']:
        old=json.loads((P/predecessor).read_text());check('predecessor PASS '+predecessor,old['outcome']=='PASS')
        for name,h in old['source_sha256'].items():check('predecessor seal '+predecessor+' '+name,sha(P/name)==h)
    manifest=json.loads((P/'sources/high-index-budget/manifest.json').read_text())
    for item in manifest['sources']:
        f=P/item['path'];check('exact primary/context '+f.name,sha(f)==item['sha256'] and f.stat().st_size==item['bytes'])
    check('entire published Chau paper freshly viewed',manifest['sources'][0]['fresh_complete_PDF_pages_rendered_and_visually_inspected']==[1,2,3,4,5,6] and 'PUBLISHED' in manifest['sources'][0]['version'])
    check('complete relevant Orevkov pages freshly viewed',manifest['sources'][1]['fresh_complete_PDF_pages_rendered_and_visually_inspected']==[7,8,9,10])
    check('classical inputs declared imports and maximal-index predecessor credited',manifest['primary_theorems_imported'] and not manifest['independent_whole_Chau_or_Orevkov_proofs'] and manifest['old_maximal_index_exclusion_credited'])
    check('no withdrawn claim or invented source symmetry',not manifest['withdrawn_Jelonek_2020_claim_imported'] and not manifest['source_Cstar_action_or_target_descent_assumed'])
    # Symbolic boundary arithmetic used in the written universal proof.
    d=S.symbols('d',integer=True,positive=True);extra,overflow=S.symbols('extra overflow',integer=True,nonnegative=True)
    ident('UNBOUNDED near-max contradiction consumes one excess and another prime',d-1-((d-2+overflow)+(2+extra))-(-1-overflow-extra))
    check('near-max contradiction is strictly negative at all nonnegative excess counts',(-1-overflow-extra).is_negative)
    ident('UNBOUNDED next strip uses exact two positive costs',d-1-((d-3)+1+1))
    ident('UNBOUNDED quotient-order improvement',(d-3)-1-(d-4))
    r=S.symbols('r',integer=True,positive=True)
    ident('UNBOUNDED high-index whole-prime budget',d-1-((d-3)+(r-1)+1)-(2-r))
    g=S.symbols('g',integer=True,positive=True)
    ident('UNBOUNDED hidden quotient plus excess reserve',d-1-((d-2)+g+1)-(-g))
    check('strict half condition retained at d6/e3',not H.high_index_reserve(6,3)['strict_high_index'])
    check('fielddegree4 not silently inserted into new cap',H.high_index_reserve(4,2)['necessary_budget_passes'])
    for degree in [5,6,7,17,49,1000]:
        check('near-max sector excluded '+str(degree),not H.high_index_reserve(degree,degree-2)['necessary_budget_passes'])
        check('general index cap is necessary only '+str(degree),H.universal_index_cap(degree)==degree-3)
    for degree in [7,8,23,1000]:
        p=H.next_strip_profile(degree)
        check('next strip complete numeric ledger '+str(degree),p['ramified_index']+p['other_index']+p['own_excess']+p['other_excess']==degree-1)
        check('next strip distinct immersed-collision image remains unproved as a realization '+str(degree),p['other_image_distinct'] and p['other_image_has_collision'] and p['ramified_nonimmersive_parameter_count']==1 and not p['realization_claimed'] and not p['source_smoothness_claimed'])
    check('new partial-link reserve is stronger than bare index cap',H.universal_index_cap(17)==14 and not H.high_index_reserve(17,14,2)['necessary_budget_passes'])
    check('index-cap equality is not falsely excluded without extra hidden generators',H.high_index_reserve(17,14,1)['necessary_budget_passes'])
    grouped=H.same_image_mass_filter(12,[(3,1),(3,1),(2,1)])
    check('image-mass obstruction excludes a class with NO individually high prime',all(2*e<=12 for e,f in [(3,1),(3,1),(2,1)]) and grouped['known_mass']==8 and grouped['own_excess_units_forced']==3 and not grouped['necessary_budget_passes'])
    grouped=H.same_image_mass_filter(12,[(3,1),(2,1),(2,1)])
    check('image-mass boundary equality is necessary not a realization',grouped['known_mass']==7 and grouped['necessary_budget_passes'] and not grouped['realization_claimed'])
    grouped=H.same_image_mass_filter(20,[(4,3)])
    check('high-image mass normalization argument permits arbitrary residue degrees',grouped['strict_high_mass'] and grouped['known_mass']==12 and grouped['known_residue_one_prime_count']==0 and grouped['own_excess_units_forced']==0 and grouped['outside_image_prime_forced'])
    grouped=H.same_image_mass_filter(20,[(4,3)],whole_nonproper_set_irreducible=True)
    check('irreducible whole nonproper set rejects high TOTAL mass even with no unit-residue prime',not grouped['necessary_budget_passes'])
    for degree in [5,6,17,40]:
        mass=degree//2
        grouped=H.same_image_mass_filter(degree,[(1,mass)],whole_nonproper_set_irreducible=True)
        check('irreducible total mass and normalized excess controls '+str(degree),grouped['irreducible_total_mass_cap']==mass and grouped['irreducible_residue_normalized_excess_lower_bound']==degree-1-mass)
    for bad in [1,4,True]:
        try:H.universal_index_cap(bad)
        except ValueError:check('universal cap domain guarded '+str(bad),True)
        else:raise AssertionError('unsupported cap domain')
    for bad in [5,6,False]:
        try:H.next_strip_profile(bad)
        except ValueError:check('strict next-strip domain guarded '+str(bad),True)
        else:raise AssertionError('unsupported next-strip domain')
    check('degree5 has no nontrivial quotient order after combined bound',H.universal_index_cap(5)-1==1)
    check('degree6 only possible quotient is A1 on index3',H.universal_index_cap(6)-1==2 and (-3)%2==1)
    # Genuine finite wrong-open example retains the failure of the WHOLE A2 premise.
    wrong=N.a1_wrong_open_model();v,t=wrong['v'],wrong['t']
    check('genuine degree12/e9 finite quotient control is not a Keller source',wrong['finite_degree']==12 and wrong['transverse_index']==9 and wrong['complement']=='C* times C' and not wrong['actual_Keller_source'])
    ident('genuine wrong-open injective cuspidal image P',wrong['P'].subs(v,0)-t**2)
    ident('genuine wrong-open injective cuspidal image Q',wrong['Q'].subs(v,0)-t**3)
    # Exact field and finite seed, not a scan over groups or cover degrees.
    x=S.symbols('x')
    check('GF16 defining polynomial irreducible over F2',S.Poly(x**4+x+1,x,modulus=2).is_irreducible)
    check('chosen GF16 element has exact order15',H.gf16_power(2,15)==1 and all(H.gf16_power(2,n)!=1 for n in [1,3,5]))
    check('all nonzero GF16 inverses exact',all(H.gf16_mul(a,H.gf16_inverse(a))==1 for a in range(1,16)))
    rho=H.seventeen_generators();identity=tuple(range(17));group={identity};queue=[identity]
    while queue:
        a=queue.pop()
        for b in rho:
            c=tuple(b[a[i]] for i in range(17))
            if c not in group:group.add(c);queue.append(c)
    check('projective seed has exact full PSL2F16 order',len(group)==16*(16*16-1)==4080)
    expectedF=S.Rational(6,17)*S.ones(3)+S.diag(0,-6,-6)
    seeds=[]
    for reverse in [False,True]:
        data=H.seventeen_cover(reverse=reverse);out=F.rational_core_linking(data);gen=W.core_generation(data)
        bs=data['boundaries'];rows=data['relation_rows'];ms=[b['meridian'] for b in bs];ls=[b['longitude'] for b in bs]
        check('exact17seed peripheral orbit types '+str(reverse),[(b['normal_index'],b['longitudinal_degree'],b['shear']) for b in bs]==[(15,1,6),(1,1,0),(1,1,0)])
        check('exact17seed full framed linking matrix '+str(reverse),S.Matrix(out['matrix'])==expectedF)
        check('exact17seed individual row defects all vanish '+str(reverse),out['core_orders']==[17]*3 and out['rational_row_denominators']==[17]*3 and out['cleared_relation_torsion_orders']==[1]*3)
        check('exact17seed still has hidden quotient '+str(reverse),gen['closed_H1']==dict(free=0,torsion=[51]) and gen['exterior_H1']==dict(free=3,torsion=[3]) and gen['unseen_quotient']==dict(free=0,torsion=[3]) and not gen['all_cores_generate'])
        for label,extra_rows,want in [('closed',ms,dict(free=0,torsion=[51])),('outside',ms+ls,dict(free=0,torsion=[3]))]:
            check('independent unreduced integer Smith17 '+str((reverse,label)),W.group_of_columns(S.Matrix(rows+extra_rows).T)==want)
        seeds.append(dict(rho=rho,peripheral=out,generation=gen,integer_relation_rows=rows,integer_meridian_rows=ms,integer_longitude_rows=ls,
                          meridian_monodromy=data['meridian_monodromy'],longitude_monodromy=data['longitude_monodromy']))
    check('both trees agree on complete17seed',seeds[0]['peripheral']==seeds[1]['peripheral'] and seeds[0]['generation']==seeds[1]['generation'])
    seed=seeds[0];mu=seed['meridian_monodromy'];lam=seed['longitude_monodromy']
    check('old b1 sector fails, new b2 sector succeeds',math.gcd(15,1-6)==5 and math.gcd(15,2-6)==1)
    try:W.linking_module(expectedF,generating_cores_certified=False)
    except ValueError:check('zero row defects cannot certify complete F reconstruction',True)
    else:raise AssertionError('false full generation accepted')
    cable_controls=[]
    for cables in [[(16,17)],[(31,17)],[(16,47),(31,77)],[(136,47)]]:
        data=H.seventeen_cover(cables);out=F.rational_core_linking(data);gen=W.core_generation(data);theory=H.seventeen_transfer(cables)
        check('full peripheral action preserved17 '+str(cables),data['meridian_monodromy']==mu and data['longitude_monodromy']==lam)
        check('direct cable linking equals universal formula17 '+str(cables),out['matrix']==theory['matrix'])
        check('direct core/group transfer17 '+str(cables),out['closed_H1']==theory['closed_H1'] and out['core_orders']==theory['core_orders'])
        check('direct unseen quotient and exterior torsion17 '+str(cables),gen['unseen_quotient']==theory['hidden_quotient'] and gen['exterior_H1']==theory['exterior_H1'])
        check('entire family excluded by all-degree index bound '+str(cables),theory['known_ramified_index']>theory['all_degree_index_cap'] and not H.high_index_reserve(17,15,gen['required_outside_prime_count'])['necessary_budget_passes'])
        cable_controls.append(dict(cables=cables,peripheral=out,generation=gen,unbounded_formula=theory))
    killed=next(z for z in cable_controls if z['cables']==[(136,47)])
    check('windings can kill every core while closed homology persists',killed['peripheral']['core_orders']==[1]*3 and killed['generation']['unseen_quotient']==dict(free=0,torsion=[51]))
    for bad in [[(16,16)],[(17,17)],[(16,32)],[(1,17)]]:
        try:H.seventeen_transfer(bad)
        except ValueError:check('out-of-sector cable rejected17 '+str(bad),True)
        else:raise AssertionError('unproved cable sector accepted')
    text=' '.join((P/'HIGH_INDEX_COLLISION_EXCESS_AND_NEAR_MAXIMAL_EXCLUSION.md').read_text().lower().split())
    for phrase in ['every generic degree d>=5','e <= d-3','floor(d/2)','another omitted prime has a different image','its own finite excess','the entire e=d-2','only to an algebraically smooth embedded a1','whole exceptional-value','2e>d','no historical priority','the next boundary strip','not an exclusion or a realization','order n<=d-4','degree-five actual normalizations are smooth','not merely orevkov','wrong open c* times c','remain open']:
        check('written scope '+phrase,phrase in text)
    files=['HIGH_INDEX_COLLISION_EXCESS_AND_NEAR_MAXIMAL_EXCLUSION.md','high_index_budget.py','verify_high_index_budget.py','sources/high-index-budget/manifest.json','whole_plane_peripheral_receipt.json','cyclic_quotient_morse_receipt.json','global_boundary_receipt.json']+[z['path'] for z in manifest['sources']]
    result=dict(claim='PF2RULED112',outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={f:sha(P/f) for f in dict.fromkeys(files)},
        theorem=dict(all_normal_actual_sources=True,all_field_degrees_main_high_index_theorem=True,
                     whole_marked_plane_and_all_omitted_primes_required=True,
                     high_index_residue_one_and_injective_normalization=True,
                     high_index_requires_positive_own_excess=True,
                     high_index_requires_distinct_nonproper_image=True,
                     high_index_e_plus_r_at_most_d_minus1=True,
                     arbitrary_residue_high_image_mass_normalization_injective=True,
                     high_image_mass_plus_unit_residue_count_at_most_d_minus2=True,
                     irreducible_whole_nonproper_TOTAL_mass_half_degree_bound=True,
                     irreducible_whole_nonproper_normalized_excess_lower_bound=True,
                     all_omitted_indices_at_most_d_minus3_for_d_at_least5=True,
                     complete_near_maximal_e_equals_d_minus2_class_excluded=True,
                     irreducible_whole_nonproper_set_half_degree_bound=True,
                     complete_next_strip_necessary_profile=True,
                     actual_quotient_order_at_most_d_minus4=True,
                     degree5_normalization_smooth=True,degree6_only_A1_quotients=True,
                     stronger_partial_link_excess_reserve=True,
                     complete_degree17_stated_colored_iterated_family_excluded=True,
                     arbitrary_cable_length_and_stated_parameters=True),
        exact17_seed=seed,cable_controls=cable_controls,
        assurance=dict(written_unbounded_proof=True,classical_primary_results_imported=True,
                       independently_reproved_Chau_and_Orevkov=False,finite_sweep_proves_theorem=False,
                       old_maximal_index_exclusion_credited=True,
                       singular_normal_sources_excluded_by_assumption=False,
                       contractible_component_alone_forbidden=False,
                       whole_nonproper_set_confused_with_branch_curve=False,
                       next_strip_all_models_excluded=False,numerical_profile_realized=False,
                       all_degree5_or_degree6_maps_excluded=False,all_degree17_maps_excluded=False,
                       uncolored_knots_without_representation_excluded=False,
                       arbitrary_finite_etale_open_suffices=False,
                       source_Cstar_action_or_target_descent_assumed=False,
                       withdrawn_Jelonek_2020_claim_imported=False,
                       global_finite_algebra_closed=False,formal_verification=False,external_review=False,
                       external_review_deferred=True,research_goal_completed=False),
        environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'high_index_budget_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',len(checks),'PF112 controls; all-degree e<=d-3 and complete near-maximal class excluded, full goal active')


if __name__=='__main__':main()
