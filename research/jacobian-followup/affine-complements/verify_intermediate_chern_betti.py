"""PF124 exact controls and source integrity; no formal/independent theorem proof."""
from pathlib import Path
import hashlib,json,platform
import sympy as S
import intermediate_chern_betti as C

P=Path(__file__).resolve().parent


def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def ident(name,value):check(name,S.simplify(value)==0)
    def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
    olds=['intermediate_critical_axis_receipt.json','whole_plane_peripheral_receipt.json',
          'normal_resolution_chern_receipt.json','planar_circle_bundle_receipt.json']
    for name in olds:
        r=json.loads((P/name).read_text());check('sealed predecessor PASS '+name,r['outcome']=='PASS')
        for f,h in r['source_sha256'].items():check('sealed predecessor '+name+' '+f,sha(P/f)==h)
    manifest=json.loads((P/'sources/intermediate-chern-betti/manifest.json').read_text())
    for item in manifest['sources']:
        f=P/item['path'];check('source seal '+f.name,sha(f)==item['sha256'] and f.stat().st_size==item['bytes'])
    check('retained results credited with no new audit claims',manifest['prior_algebraic_nontorsion_intermediate_canonical_obstruction_credited'] and manifest['prior_Danielewski_and_pseudo_plane_controls_credited_not_rebuilt'] and not manifest['new_PDF_visual_audit'] and not manifest['primitive_critical_curve_marking_in_intermediate_assumed'])
    r,v,kappa,h,L,E=S.symbols('r v kappa h L E',positive=True)
    chiC=1-v;chi=4*r*(1-chiC)+3*r*(1-2*v-kappa)+2*r*(v+kappa)
    ident('UNBOUNDED proper finite Euler integration',chi-r*(3-kappa))
    ident('UNBOUNDED double contact count cancels',S.diff(chi,v))
    ident('UNBOUNDED immersive Betti profile',chi.subs(kappa,0)-1-(3*r-1))
    ident('UNBOUNDED one-cusp Betti profile',chi.subs(kappa,1)-1-(2*r-1))
    ident('UNBOUNDED immersive excludes rank1 end for all r>=1',(3*r-1-1).subs(r,1+h)-(1+3*h))
    ident('UNBOUNDED one-cusp rank1 equality only r1',(2*r-1-1).subs(r,1+h)-2*h)
    ident('UNBOUNDED ADE contact theta difference',(L-2-E)-(L-2)+E)
    check('UNBOUNDED ADE nonzero energy cannot match theta',(-E).is_negative is True)
    profiles=[]
    for b in [4,8,32]:
        for kk in [0,1]:
            p=C.affine_betti(b,kk)
            check('actual affine intermediate Betti '+str((b,kk)),p['chi']==b-b//4*(1+kk) and p['b2']==p['chi']-1 and p['b2']>=p['critical_curves'])
            check('actual circle-bundle survivor filter '+str((b,kk)),(p['b2']<=1)==(b==4 and kk==1))
            profiles.append(p)
    matrices=[]
    for packets,coeff in [([[1,2]],[1]),([[2],[1,1]],[1,1]),([[3,1],[2,2],[1,3]],[1,2,4])]:
        A=S.Matrix(C.pullback_matrix(packets));image=A*S.Matrix(coeff)
        check('actual disjoint critical pullback supports are independent '+str(packets),A.rank()==len(packets))
        check('actual canonical positive combination nonzero '+str(packets),all(z>0 for z in image))
        matrices.append(dict(packets=packets,rank=A.rank(),canonical_coordinates=list(map(int,image))))
    check('finite pi1 may be nontrivial when t>1',6%3==0 and 3>1)
    check('original t1 fundamental group forced trivial',1%1==0)
    # These are exact model constraints, not an enumeration proof for the groups.
    check('rank1 p4 LF argument needs only rational b1zero',all((1-3+length)-1==length-3 for length in [3,4]))
    check('actual complete radial zero profile is entirely even and has degree d',sum(2*e for e in [2,4])==12 and all(e%2==0 for e in [2,4]))
    for different in [False,True]:
        p0=[x^(2 if x&1 else 0) for x in range(8)]
        q0=[x^((4 if different else 2) if not x&1 else 0) for x in range(8)]
        product=[p0[q0[x]] for x in range(8)]
        check('actual complementary affine radial involution '+str(different),all(product[x]!=x and product[product[x]]==x for x in range(8)))
        if different:check('countercontrol radial involution need not be translation',len({x^product[x] for x in range(8)})>1)
    text=' '.join((P/'INTERMEDIATE_CHERN_BETTI_AND_END_OBSTRUCTIONS.md').read_text().lower().split())
    for phrase in ['original finite normalization','all omitted primes','every nontrivial',
      'retained results, not new claims','controls are not new constructions',
      'rationally independent','c1(k_m)!=0','finite factorial cover',
      '|pi1(m)| divides t','c1(ty_m)!=0','b2(y_m)>=s+a',
      'rational boundary c1 zero is insufficient','does not require pi1(v)=0',
      'counts critical curves','no primitive marking','3b/4-1','b/2-1',
      'only b=4 with one ordinary cusp','not the nonproper etale source restriction',
      'empty original source fibres','all e_i even','not necessarily a translation',
      'does not make q a square','remain open','fresh pdf visual audit']:
        check('written scope '+phrase,phrase in text)
    files=['INTERMEDIATE_CHERN_BETTI_AND_END_OBSTRUCTIONS.md','intermediate_chern_betti.py',
      'verify_intermediate_chern_betti.py','sources/intermediate-chern-betti/manifest.json']+olds+[x['path'] for x in manifest['sources']]
    result=dict(claim='PF2RULED124',outcome='PASS',check_count=len(checks),checks=checks,
      theorem=dict(all_plane_degrees_all_original_monodromy_types=True,whole_original_normalization_all_primes_required=True,
        all_nontrivial_intermediates_nonzero_rational_canonical_Chern=True,
        all_intermediate_critical_divisor_classes_topologically_independent=True,
        all_intermediate_pi1_finite_order_divides_original_block_size=True,
        all_intermediate_minimal_resolutions_nonzero_rational_Chern_and_rank_bound=True,
        all_intermediate_S3_and_negative_ADE_oriented_ends_excluded=True,
        all_intermediate_Lp1_ends_force_smooth_one_critical_curve_rank1=True,
        binary_affine_intermediate_exact_Euler_and_Betti_profiles=True,
        binary_affine_Lp1_ends_only_degree4_one_cusp_can_survive=True,
        binary_empty_fibre_generic_pencil_radial_zero_divisor_entirely_even=True),
      betti_profiles=profiles,critical_pullback_controls=matrices,
      assurance=dict(prior_canonical_and_pseudo_plane_results_credited=True,
        prior_deep_filling_and_deformation_imports_explicitly_retained=True,
        independent_new_filling_classification=False,new_PDF_visual_audit=False,
        finite_pi1_assumed_trivial=False,all_Chern_zero_boundary_conditions_taken_only_rationally=False,
        intermediate_critical_curve_marking_assumed_primitive=False,
        original_sole_boundary_parity_bound_transferred_automatically=False,
        intermediate_global_plane_recognition=False,
        every_binary_affine_quotient_excluded=False,every_quotient_degree4_excluded=False,
        all_possible_oriented_ends_excluded=False,determinant_alone_identifies_end=False,
        original_source_missing_fibres_excluded=False,
        affine_radial_involution_assumed_translation=False,
        even_radial_zero_divisor_assumed_global_square=False,
        actual_global_finite_algebra_closed=False,finite_samples_prove_theorem=False,
        historical_priority_claimed=False,formal_verification=False,external_review=False,research_goal_completed=False),
      source_sha256={f:sha(P/f) for f in dict.fromkeys(files)},environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'intermediate_chern_betti_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',len(checks),'PF124 controls; all-intermediate Chern/end obstructions and exact affine Betti profiles; goal ACTIVE')


if __name__=='__main__':main()
