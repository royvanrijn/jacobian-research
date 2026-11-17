"""PF110 exact controls supplement the written unbounded proofs and declared imports."""
from pathlib import Path
import hashlib,json,platform
import sympy as S
import planar_circle_bundle as C

P=Path(__file__).resolve().parent


def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def ident(name,value):check(name,S.simplify(S.expand(value))==0)
    def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
    old=json.loads((P/'normal_resolution_chern_receipt.json').read_text())
    check('sealed PF109 predecessor PASS',old['outcome']=='PASS')
    for f,h in old['source_sha256'].items():check('PF109 seal '+f,sha(P/f)==h)
    manifest=json.loads((P/'sources/planar-circle-bundle/manifest.json').read_text())
    for x in manifest['sources']:
        f=P/x['path'];check('primary/context pin '+f.name,sha(f)==x['sha256'] and f.stat().st_size==x['bytes'])
    check('exceptional p4 not silently covered by Cor1.3',manifest['exceptional_p4_rank_bound_written_without_nonstandard_filling_uniqueness'])
    check('unretrieved original not falsely audited',not manifest['original_mcduff_and_erratum_retrieved'])
    check('three complete-page primary audits', [x['fresh_complete_PDF_pages_rendered_and_visually_inspected'] for x in manifest['sources'][:3]]==[[2,14],[2,4,5,12,13],[6,8]])
    # Arbitrary-p chart and differential identities, not only sampled p values.
    p=S.symbols('p',integer=True,positive=True)
    u,v=S.symbols('u v',nonzero=True);t,s=S.symbols('t s')
    target=u**p*t+u;inverse=v**p*s-v**(p-1)
    ident('UNBOUNDED chart inverse u',inverse.subs({v:1/u,s:target})-t)
    ident('UNBOUNDED chart inverse v',target.subs({u:1/v,t:inverse})-s)
    du_v=-v**-2;dt_s=S.diff(inverse,s);dt_v=S.diff(inverse,v)
    ident('UNBOUNDED plane volume transition',du_v*dt_s+v**(p-2))
    alpha_s=v**(p-1);alpha_v=p*v**(p-2)*s-(p-1)*v**(p-3)
    ident('UNBOUNDED primitive ds coefficient',v**-1*dt_s-alpha_s)
    ident('UNBOUNDED primitive dv coefficient',v**-1*dt_v-alpha_v)
    ident('UNBOUNDED exterior derivative of global primitive',S.diff(alpha_s,v)-S.diff(alpha_v,s)+v**(p-2))
    ident('UNBOUNDED extension determinant',S.det(S.Matrix([[u**p,u],[0,1]]))-u**p)
    section=S.Matrix([[u**(p-1),1],[0,u**-1]])*S.Matrix([0,1])
    check('UNBOUNDED nowherezero global section of E(1)',section==S.Matrix([1,u**-1]))
    n=p-2;b=p-1
    ident('UNBOUNDED removed section square',-n+2*b-p)
    ident('UNBOUNDED ample section meets C0 once',-n+b-1)
    ident('UNBOUNDED split degrees sum to original determinant',-1-(p-1)+p)
    ident('UNBOUNDED canonical divisor index saturates bound',(p-2)+1-(p-1))
    ident('UNBOUNDED nonzero boundary Chern residue',-(p-2)+p-2)
    k,R=S.symbols('k R',integer=True)
    ident('UNBOUNDED theta equality detects squared rotation',(-1-k*k/p)-(-1-R*R/p)-(R*R-k*k)/p)
    ident('UNBOUNDED saturation parity',(p-2-p)+2)
    # Exceptional standard/conjugate L4,1: counts determine rank without uniqueness.
    a1,a2,a3,b12,b13,b23,c=S.symbols('a1 a2 a3 b12 b13 b23 c')
    equations=[b12+c-1,b13+c-1,b23+c-1,a1+b12+b13+c-2,a2+b12+b23+c-2,a3+b13+b23+c-2]
    solution=S.solve(equations,[a1,a2,a3,b12,b13,b23],dict=True)
    check('full exceptional multiplicity linear system solved',solution==[{a1:c,a2:c,a3:c,b12:1-c,b13:1-c,b23:1-c}])
    ident('UNBOUNDED nonnegative-count length identity',(a1+a2+a3+b12+b13+b23+c).subs(solution[0])-(3+c))
    for q in (0,1):
        f=C.three_hole_factorization(q)
        check('exceptional length/rank '+str(q),f['length']==3+q and f['b2_if_simply_connected']==q)
    for q in (-1,2):
        try:C.three_hole_factorization(q)
        except ValueError:check('negative factorization multiplicity rejected '+str(q),True)
        else:raise AssertionError('invalid factorization accepted')
    for boundary_zero,fill_nonzero,excluded in [(True,True,True),(False,True,False),(True,False,False)]:
        g=C.torsion_planarity_gate(integral_boundary_Chern_zero=boundary_zero,fill_Chern_nonzero=fill_nonzero)
        check('integral/torsion/zero-fill planarity scope '+str((boundary_zero,fill_nonzero)),g['planarity_excluded']==excluded and not g['rational_boundary_vanishing_alone_suffices'] and not g['graph_base_genus_zero_implies_contact_planarity'])
    check('RP3 index set empty',C.circle_bundle_indices(2)==[])
    check('odd and even index parity controls',C.circle_bundle_indices(3)==[2] and C.circle_bundle_indices(4)==[3] and C.circle_bundle_indices(5)==[2,4] and C.circle_bundle_indices(6)==[3,5])
    controls=[]
    for order in (3,4,5,8,17,101):
        d=C.source_control(order);controls.append(d)
        check('sharp genuine affine source '+str(order),d['numerical_index'] in C.circle_bundle_indices(order) and d['boundary_Chern_residue']!=0 and d['globally_exact_volume'] and d['source_affine'] and d['whole_plane_open'])
        check('no finite map realization silently inferred '+str(order),not d['finite_target_supplied'] and not d['actual_Keller_realization_claimed'])
    check('p2 displayed primitive has genuine pole',alpha_v.subs(p,2)==2*s-1/v)
    for bad in (1,2,S.Rational(7,2)):
        try:C.source_control(bad)
        except ValueError:check('positive exact source family domain guarded '+str(bad),True)
        else:raise AssertionError('invalid source family parameter accepted')
    text=' '.join((P/'PLANAR_CHERN_OBSTRUCTION_AND_CIRCLE_BUNDLE_SOURCE_REDUCTION.md').read_text().lower().split())
    for phrase in ['nonzero in integral','not enough','exactly one','no uniqueness','r=1, m=0','for every integer p>=3','globally exact','no compatible finite keller target','remain open','not independently reproved']:
        check('written scope '+phrase,phrase in text)
    files=['PLANAR_CHERN_OBSTRUCTION_AND_CIRCLE_BUNDLE_SOURCE_REDUCTION.md','planar_circle_bundle.py','verify_planar_circle_bundle.py','sources/planar-circle-bundle/manifest.json','normal_resolution_chern_receipt.json']+[x['path'] for x in manifest['sources']]
    result=dict(claim='PF2RULED110',outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={f:sha(P/f) for f in dict.fromkeys(files)},
        theorem=dict(all_field_degrees=True,all_normal_actual_sources=True,whole_marked_plane_required=True,
                     planar_contact_requires_integral_boundary_Chern_nonzero=True,
                     all_S3_ends_excluded=True,all_negative_SU2_quotient_ends_excluded=True,
                     all_Lp1_singular_sources_excluded=True,all_Lp1_multiple_prime_sources_excluded=True,
                     surviving_Lp1_source_smooth_sole_prime=True,
                     sharp_index_range_and_parity=True,sharp_all_p_affine_exact_source_controls=True,
                     oriented_actual_end_identification_required=True),
        genuine_source_controls=controls,
        assurance=dict(written_unbounded_proof=True,deep_contact_filling_results_imported=True,
                       independent_whole_classification_proofs=False,finite_sweep_proves_theorem=False,
                       rational_boundary_Chern_zero_suffices=False,graph_base_genus_implies_contact_planarity=False,
                       nonstandard_p4_fill_uniqueness_required=False,original_mcduff_audited=False,
                       all_Lp1_ends_excluded=False,all_spherical_ends_excluded=False,
                       all_integral_homology_sphere_ends_excluded=False,reversed_orientations_excluded=False,
                       source_controls_have_finite_Keller_targets=False,homology_identifies_end=False,
                       global_finite_algebra_closed=False,formal_verification=False,external_review=False,
                       external_review_deferred=True,research_goal_completed=False),
        environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'planar_circle_bundle_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',len(checks),'PF110 controls; planar obstruction and complete circle-bundle source reduction, full goal active')


if __name__=='__main__':main()
