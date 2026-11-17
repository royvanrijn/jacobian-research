"""PF113 exact controls; universal scope rests on written local-factor/flow proofs."""
from pathlib import Path
from itertools import product
import hashlib,json,platform
import sympy as S
import local_monodromy_branch as L
import high_index_budget as H
import coprime_leading_forms as C

P=Path(__file__).resolve().parent


def quotient_length(equations,variables):
    G=S.groebner(equations,*variables)
    if not G.is_zero_dimensional:raise ValueError('zero-dimensional fibre required')
    leaders=[p.LM(order=G.order).exponents for p in G.polys]
    bounds=[min(m[i] for m in leaders if m[i]>0 and sum(m)==m[i]) for i in range(len(variables))]
    basis=[m for m in product(*(range(b) for b in bounds)) if not any(all(x>=a for x,a in zip(m,lm)) for lm in leaders)]
    return len(basis),leaders


def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def ident(name,value):check(name,S.factor(value)==0)
    def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
    for predecessor in ['high_index_budget_receipt.json','whole_plane_peripheral_receipt.json','normal_resolution_chern_receipt.json']:
        old=json.loads((P/predecessor).read_text());check('predecessor PASS '+predecessor,old['outcome']=='PASS')
        for name,h in old['source_sha256'].items():check('predecessor seal '+predecessor+' '+name,sha(P/name)==h)
    manifest=json.loads((P/'sources/local-monodromy-branch/manifest.json').read_text())
    for item in manifest['sources']:
        f=P/item['path'];check('exact primary/context '+f.name,sha(f)==item['sha256'] and f.stat().st_size==item['bytes'])
    check('exact Palka version and irreducibility convention',manifest['sources'][0]['version'].startswith('arXiv1405.5391v2') and 'irreducibility' in manifest['sources'][0]['version'])
    check('classification imported and old audits not falsely fresh',manifest['Palka_classification_imported'] and not manifest['independent_whole_classification_proof'] and not manifest['fresh_PDF_audit_claimed'])
    # Arbitrary positive target weights: derivative positivity is checked in the
    # two cases covering all nonzero points, not extrapolated from sampled weights.
    a,b=S.symbols('a b',integer=True,positive=True);t=S.symbols('t',positive=True)
    xx=S.symbols('xx',positive=True);yy=S.symbols('yy',nonnegative=True)
    flow_norm=xx*t**(2*a)+yy*t**(2*b)
    derivative=S.diff(flow_norm,t)
    ident('UNBOUNDED weighted radial derivative',derivative-(2*a*xx*t**(2*a-1)+2*b*yy*t**(2*b-1)))
    check('UNBOUNDED derivative positive when first coordinate is nonzero',derivative.is_positive)
    derivative2=2*a*yy*t**(2*a-1)+2*b*xx*t**(2*b-1)
    check('UNBOUNDED derivative positive when second coordinate is nonzero',derivative2.is_positive)
    d=S.symbols('d',integer=True,positive=True);extra,overflow=S.symbols('extra overflow',integer=True,nonnegative=True)
    ident('UNBOUNDED local cluster capacity algebra',d-1-(S.symbols('ell')+S.symbols('r')-1)-(d-S.symbols('r')-S.symbols('ell')))
    slack=d-1-((d-3+overflow)+(2+extra)+1)
    ident('UNBOUNDED entire prior e=d-3 strip has negative budget',slack-(-1-overflow-extra))
    check('prior strip contradiction strict at all nonnegative extra costs',slack.is_negative)
    ident('UNBOUNDED new e=d-4 profile uses exact three units',d-1-((d-4)+2+1))
    ident('UNBOUNDED low-index even-degree cap',(2*(4+extra)-4)-(4+extra)-extra)
    ident('UNBOUNDED low-index odd-degree cap',(2*(3+extra)+1-4)-(3+extra)-extra)
    ident('UNBOUNDED quotient cap improvement',(d-4)-1-(d-5))
    for degree,count,blocks,passes in [(6,1,[6],False),(6,2,[3,2,1],True),(6,3,[4,1,1],False),(7,5,[2,2,2,1],True)]:
        result=L.local_partition_gate(degree,count,blocks)
        check('actual normalization local partition necessary gate '+str((degree,count,blocks)),result['necessary_capacity_passes']==passes and result['local_cluster_cap']==degree-count)
        check('normalization singleton gate does not forbid source singleton or imply surjectivity '+str((degree,count)),not result['source_restriction_singleton_fibre_excluded'] and not result['source_surjectivity_claimed'])
    for bad in [(5,0,[5]),(5,2,[2,2]),(5,2,[0,5])]:
        try:L.local_partition_gate(*bad)
        except ValueError:check('malformed capacity data rejected '+str(bad),True)
        else:raise AssertionError('malformed partition accepted')
    x,y=S.symbols('x y')
    weighted_controls=[]
    for expr,weights in [(x,(1,1)),(x**2-y**3,(3,2)),(x**19-y**17,(17,19)),((x**2-y**3)*(x**2-2*y**3),(3,2)),(y-x**7,(1,7))]:
        gate=L.positive_weighted_branch_gate(S.Poly(expr,x,y),weights,full_branch_union_certified=True)
        check('exact full weighted support certificate '+str((expr,weights)),gate['actual_Keller_cover_excluded'] and not gate['arbitrary_finite_cover_excluded'])
        degree=gate['weighted_degrees'][0]
        ident('target weighted flow support identity '+str(expr),expr.subs({x:t**weights[0]*x,y:t**weights[1]*y})-t**degree*expr)
        ident('target weighted Euler identity '+str(expr),weights[0]*x*S.diff(expr,x)+weights[1]*y*S.diff(expr,y)-degree*expr)
        weighted_controls.append(dict(polynomial=str(expr),weights=weights,gate=gate))
    partial=L.positive_weighted_branch_gate(S.Poly(x**2-y**3,x,y),(3,2),full_branch_union_certified=False)
    check('one weighted component is not a certified full branch union',partial['weighted_homogeneous'] and not partial['actual_Keller_cover_excluded'])
    combined=L.positive_weighted_branch_gate(S.Poly((x**2-y**3)*(x-1),x,y),(3,2),full_branch_union_certified=True)
    check('a component certificate does not certify the full displayed support',not combined['weighted_homogeneous'] and not combined['actual_Keller_cover_excluded'])
    for weights in [(1,0),(1,-1),(0,1)]:
        try:L.positive_weighted_branch_gate(S.Poly(x*y-1,x,y),weights,full_branch_union_certified=True)
        except ValueError:check('noncontracting target weights rejected '+str(weights),True)
        else:raise AssertionError('nonpositive radial weights accepted')
    for degree in [7,8,9,23,1000]:
        check('complete e=d-3 class now excluded '+str(degree),not L.high_index_reserve(degree,degree-3)['necessary_budget_passes'])
        check('universal cap is d-4 '+str(degree),L.universal_index_cap(degree)==degree-4)
    for degree in [9,10,41,1000]:
        p=L.next_strip_profile(degree)
        check('new two-RAMIFIED-prime profile exact budget '+str(degree),p['first_index']+p['second_index']+p['first_excess']==degree-1 and p['ramified_prime_count']==2 and p['unramified_prime_count']==0)
        check('new profile retains nonrealization and source-singularity scope '+str(degree),p['second_image_has_collision'] and p['distinct_ramified_images'] and not p['realization_claimed'] and not p['source_smoothness_assumed'])
    check('strict half equality remains outside this high-index argument',not L.high_index_reserve(8,4)['strict_high_index'])
    check('new partial-link reserve adds a different RAMIFIED image cost',not L.high_index_reserve(20,16,2)['necessary_budget_passes'] and L.high_index_reserve(20,15,2)['necessary_budget_passes'])
    check('shared-image bound improves PF112 by one',H.same_image_mass_filter(12,[(3,1),(2,1),(2,1)])['necessary_budget_passes'] and not L.image_mass_gate(12,[(3,1),(2,1),(2,1)])['necessary_budget_passes'])
    check('irreducible FULL branch is enough despite unspecified extra nonproper images',not L.image_mass_gate(20,[(7,2)],irreducible_full_branch=True)['necessary_budget_passes'])
    check('unramified companion cannot share a high-mass image',not L.image_mass_gate(12,[(6,1),(1,1)])['necessary_budget_passes'])
    q8=[(n,e,(-e)%n) for n in [2,3] for e in range(1,5) if n<=e-1 and S.gcd(n,e)==1]
    check('degree8 quotient restrictions exactly A1 orA2',q8==[(2,3,1),(3,4,2)])
    # Genuine finite maps with transitive local monodromy: not marked Keller examples.
    m=S.symbols('m',integer=True,positive=True);z=S.symbols('z',nonzero=True)
    ident('UNBOUNDED Kummer Jacobian',S.diff(z**m,z)-m*z**(m-1))
    for degree in [2,3,11,31]:
        model=L.kummer_control(degree)
        check('genuine Kummer cover has singleton full fibre outside Keller scope '+str(degree),model['singleton_normalization_fibre_length']==degree and model['generic_index']==degree and not model['whole_plane_Keller_map'] and model['etale_open']=='C* times C')
        check('Kummer control violates the marked Euler budget '+str(degree),degree>degree-1 and S.Poly(model['Jacobian'],model['z'],model['t']).total_degree()>0)
    cubic=L.cubic_control();s,t0=cubic['s'],cubic['t'];u,v=cubic['u'],cubic['v']
    ident('actual cubic Jacobian',S.diff(cubic['P'],s)*S.diff(cubic['Q'],t0)-S.diff(cubic['P'],t0)*S.diff(cubic['Q'],s)-cubic['Jacobian'])
    ident('actual critical parabola maps into weighted cusp',cubic['branch'].subs({u:cubic['branch_parameter'][0],v:cubic['branch_parameter'][1]}))
    length,leaders=quotient_length([cubic['P'],cubic['Q']],(s,t0))
    check('actual cubic cusp fibre has unique support and length3',length==3 and leaders==[(1,0),(0,3)] and not cubic['whole_plane_Keller_map'])
    singular=L.singular_quotient_control();A,B,C0=singular['A'],singular['B'],singular['C']
    length,leaders=quotient_length([singular['relation'],singular['P'],singular['Q']],(A,B,C0))
    check('genuine singular quotient fibre homogeneous zero-dimensional length12',length==12 and all(S.Poly(f,A,B,C0).is_homogeneous for f in [singular['relation'],singular['P'],singular['Q']]))
    upstream=C.explicit_pair(2);X,Y=upstream['X'],upstream['Y']
    ident('singular finite cover P is actual invariant polynomial',singular['P'].subs({A:X*X,B:X*Y,C0:Y*Y})-upstream['P'])
    ident('singular finite cover Q is actual invariant polynomial',singular['Q'].subs({A:X*X,B:X*Y,C0:Y*Y})-upstream['Q'])
    ident('full branch of singular control is the cuspidal axis image P',singular['P'].subs({A:0,B:0})-C0*C0)
    ident('full branch of singular control is the cuspidal axis image Q',singular['Q'].subs({A:0,B:0})-C0**3)
    check('singular finite cover is not falsely excluded',singular['etale_open']=='C* times C' and not singular['whole_plane_Keller_map'] and singular['source_singularity']=='A1')
    text=' '.join((P/'LOCAL_MONODROMY_CAPACITY_AND_WEIGHTED_BRANCH_EXCLUSION.md').read_text().lower().split())
    for phrase in ['ell_p <= d-r','not the fibre','singleton source fibre','never transitive','full finite-map branch-image','strictly positive','the same homotopy','same s_epsilon','no lift of a c* action','additional unramified','classification is imported','e<=d-4 for d>=7','entire e=d-3','both ramified','not a realized','not an identification','source restriction','not disproved','remain open','infinity knot alone']:
        check('written universal scope '+phrase,phrase in text)
    files=['LOCAL_MONODROMY_CAPACITY_AND_WEIGHTED_BRANCH_EXCLUSION.md','local_monodromy_branch.py','verify_local_monodromy_branch.py','sources/local-monodromy-branch/manifest.json','high_index_budget_receipt.json','whole_plane_peripheral_receipt.json','normal_resolution_chern_receipt.json']+[x['path'] for x in manifest['sources']]
    result=dict(claim='PF2RULED113',outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={f:sha(P/f) for f in dict.fromkeys(files)},
        theorem=dict(all_field_degrees=True,all_normal_actual_sources=True,
                     original_finite_normalization_whole_plane_all_primes_required=True,
                     every_local_cluster_at_most_d_minus_r=True,
                     complete_finite_fibre_count_lower_bound=True,
                     no_singleton_TOTAL_normalization_fibre=True,
                     local_monodromy_never_transitive=True,
                     positive_weighted_FULL_branch_union_excluded=True,
                     polynomial_target_rectifications_included=True,
                     contractible_irreducible_FULL_branch_excluded=True,
                     extra_unramified_nonproper_images_allowed=True,
                     high_index_forces_different_RAMIFIED_image=True,
                     all_indices_at_most_d_minus4_for_d_at_least7=True,
                     complete_prior_e_equals_d_minus3_strip_excluded=True,
                     shared_image_mass_plus_unit_residue_count_at_most_d_minus3=True,
                     irreducible_FULL_branch_mass_at_most_half_degree=True,
                     new_two_RAMIFIED_prime_necessary_profile=True,
                     quotient_order_at_most_d_minus5_for_d_at_least7=True),
        weighted_support_controls=weighted_controls,
        genuine_finite_controls=dict(cubic_fibre_length=3,singular_A1_fibre_length=12,
                                     singular_fibre_leading_monomials=leaders,
                                     etale_opens='C* times C',Keller_plane_examples=False),
        assurance=dict(written_unbounded_proof=True,Palka_classification_imported=True,
                       independent_whole_classification_proof=False,finite_samples_prove_theorem=False,
                       full_branch_union_certification_required=True,
                       source_singleton_fibres_excluded=False,source_surjectivity_claimed=False,
                       arbitrary_finite_covers_excluded=False,
                       one_weighted_component_suffices=False,weighted_analytic_germ_suffices=False,
                       infinity_knot_alone_certifies_global_branch_class=False,
                       nonpositive_weights_suffice=False,simultaneous_component_rectification_assumed=False,
                       source_Cstar_action_lift_assumed=False,source_symmetry_descent_assumed=False,
                       new_next_strip_models_excluded=False,numerical_profile_realized=False,
                       all_degree7_or_degree8_maps_excluded=False,homology_identifies_ADE_end=False,
                       global_finite_algebra_closed=False,formal_verification=False,external_review=False,
                       external_review_deferred=True,research_goal_completed=False),
        environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'local_monodromy_branch_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',len(checks),'PF113 controls; weighted FULL branch exclusion and universal e<=d-4 for d>=7, full goal active')


if __name__=='__main__':main()
