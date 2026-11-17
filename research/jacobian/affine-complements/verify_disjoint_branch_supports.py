"""PF121 small controls; no group census or proof by finite samples."""
from pathlib import Path
import hashlib,json,platform
import sympy as S
import disjoint_branch_supports as L

P=Path(__file__).resolve().parent


def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def ident(name,value):check(name,S.expand(value)==0)
    def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
    olds=['primitive_product_action_receipt.json','simple_ramification_ledger_receipt.json',
          'half_support_inertia_receipt.json','local_monodromy_branch_receipt.json']
    for name in olds:
        r=json.loads((P/name).read_text());check('sealed predecessor PASS '+name,r['outcome']=='PASS')
        for f,h in r['source_sha256'].items():check('sealed predecessor '+name+' '+f,sha(P/f)==h)
    manifest=json.loads((P/'sources/disjoint-branch-supports/manifest.json').read_text())
    for item in manifest['sources']:
        f=P/item['path'];check('source seal '+f.name,sha(f)==item['sha256'] and f.stat().st_size==item['bytes'])
    check('primary web import source boundary',manifest['Guralnick_Magaard_CFSG_theorem_imported'] and
          manifest['narrow_all_element_inclusive_half_property_only'] and
          not manifest['original_Guralnick_Magaard_PDF_obtained_or_visually_audited'] and
          not manifest['old_OCR_subset_bound_used'] and not manifest['independent_CFSG_proof'])
    a,b=S.symbols('a b',positive=True)
    ident('UNBOUNDED shifted torus fibre Riemann-Hurwitz genus',-2*b+(a+1)*(b-1)-((a-1)*(b-1)-2))
    genera=[]
    for aa,bb in [(2,3),(3,4),(4,7),(5,9)]:
        genus=(aa-1)*(bb-1)//2;check('positive shifted torus-fibre genus '+str((aa,bb)),genus>=1)
        genera.append(dict(a=aa,b=bb,genus=genus))
    stars=[]
    for m,k,i,j in [(4,2,0,1),(5,2,2,4),(7,3,0,6),(8,2,3,6),(8,6,3,6),(10,4,1,9)]:
        family,sigs=L.subset_signatures(m,k,i,j)
        check('actual oriented-star family separates ALL letters '+str((m,k,i,j)),len(set(sigs))==m)
        check('actual distinguished signatures '+str((m,k,i,j)),all(sigs[i]) and not any(sigs[j]))
        stars.append(dict(m=m,k=k,family_size=len(family),distinct_signatures=len(set(sigs))))
    # Endpoint controls: the theorem deliberately excludes k=1,m-1.
    check('natural action actually permits disjoint nonidentity supports',{0,1}.isdisjoint({2,3}))
    A=S.symbols('A',positive=True)
    ident('UNBOUNDED common classical overlap polynomial',A*(A/2-1)-(A*A/2-A))
    h=S.symbols('h',nonnegative=True)
    poly=S.Poly(S.expand((A*A/2-A).subs(A,4+h)),h)
    check('UNBOUNDED strictly positive classical overlap for A>=4',all(c>=0 for c in poly.all_coeffs()) and poly.eval(0)>0)
    # Enumerate only small quadratic spaces, not groups or permutations of groups.
    form_controls=[]
    for aa in [3,4]:
        forms=L.minus_forms(aa);v,w=1,2
        p=L.form_permutation(aa,v);q=L.form_permutation(aa,w)
        sv,sw=L.support(p),L.support(q);overlap=2**(aa-1)*(2**(aa-2)-1)
        degree=2**(aa-1)*(2**aa-1);moved=2**(aa-1)*(2**(aa-1)-1)
        check('actual minus-form degree '+str(aa),len(forms)==degree)
        check('actual minus-form transvection support '+str(aa),len(sv)==len(sw)==moved)
        check('actual commuting minus-form supports overlap '+str(aa),len(sv&sw)==overlap>0)
        check('actual minus-form involutions '+str(aa),all(p[p[i]]==q[q[i]]==i for i in range(degree)))
        # Nonorthogonal transvections fail to commute, hence cannot be disjoint.
        r=L.form_permutation(aa,1<<aa)
        check('actual nonorthogonal form transvections do not commute '+str(aa),any(p[r[i]]!=r[p[i]] for i in range(degree)))
        form_controls.append(dict(a=aa,degree=degree,moved=moved,overlap=overlap))
    vector_controls=[]
    for aa,negative,epsilon in [(3,True,0),(4,True,0),(4,False,1)]:
        v=1 if negative else (1|(1<<aa));w=2|(2<<aa)
        points,p=L.vector_permutation(aa,negative,epsilon,v)
        same,q=L.vector_permutation(aa,negative,epsilon,w)
        sv,sw=L.support(p),L.support(q);overlap=2**(aa-1)*(2**(aa-2)-1)
        degree=(2**aa+1)*(2**(aa-1)-1) if negative else 2**(aa-1)*(2**aa-1)
        check('actual orthogonal-vector degree '+str((aa,negative)),len(points)==degree and points==same)
        check('actual orthogonal-vector transvection involutions '+str((aa,negative)),all(p[p[i]]==q[q[i]]==i for i in range(degree)))
        check('actual orthogonal-vector overlap '+str((aa,negative)),len(sv&sw)==overlap>0)
        check('actual orthogonal-vector commuting pair '+str((aa,negative)),all(p[q[i]]==q[p[i]] for i in range(degree)))
        vector_controls.append(dict(a=aa,negative=negative,degree=degree,moved=len(sv),overlap=overlap))
    # Affine controls prevent upgrading intersecting classical supports to ALL transvections.
    for aa in [3,5]:
        v=2;points=list(range(1<<aa))
        p=[x^(v if x&1 else 0) for x in points]
        q=[x^(v if not x&1 else 0) for x in points]
        sv,sw=L.support(p),L.support(q)
        check('actual affine complementary transvections survive '+str(aa),not(sv&sw) and sv|sw==set(points) and len(sv)==len(sw)==len(points)//2)
    text=' '.join((P/'DISJOINT_BRANCH_SUPPORTS_AND_PRIMITIVE_MONODROMY_REDUCTION.md').read_text().lower().split())
    for phrase in ['original finite normalization','all omitted primes','arbitrary indices and singular n',
      'two geometric branch meridians','disjoint moved-sheet supports','at most one omitted-curve germ',
      'one irreducible analytic branch','extends to p1','no simultaneous coordinate rectification',
      'membership signatures','endpoints k=1,m-1','at least half','all-element, inclusive-half',
      'not visually audited','ocr subset-parameter bound','triality','generated by linear transvections',
      'imprimitive cases','remain open','not an exclusion']:
        check('written scope '+phrase,phrase in text)
    files=['DISJOINT_BRANCH_SUPPORTS_AND_PRIMITIVE_MONODROMY_REDUCTION.md','disjoint_branch_supports.py',
      'verify_disjoint_branch_supports.py','sources/disjoint-branch-supports/manifest.json']+olds+[x['path'] for x in manifest['sources']]
    result=dict(claim='PF2RULED121',outcome='PASS',check_count=len(checks),checks=checks,
      theorem=dict(all_plane_degrees=True,whole_original_normalization_all_primes_required=True,
        universal_disjoint_nonidentity_BRANCH_meridian_pair=True,
        all_faithful_interior_subset_actions_excluded=True,subset_exclusion_CFSG_free=True,
        all_three_remaining_binary_classical_primitive_actions_excluded=True,
        primitive_only_natural_alternating_symmetric_or_binary_affine=True,
        binary_affine_linear_part_generated_by_transvections=True,
        binary_affine_exact_two_half_support_cusp_profiles_retained=True),
      genus_controls=genera,subset_controls=stars,minus_form_controls=form_controls,orthogonal_controls=vector_controls,
      assurance=dict(published_CFSG_classification_imported=True,
        GM_all_element_inclusive_half_primary_web_text_imported=True,
        original_GM_PDF_audited=False,original_GM_source_bytes_locally_pinned=False,
        review_record_misrepresented_as_original_source=False,OCR_subset_parameter_bound_used=False,
        independent_CFSG_classification_proved=False,finite_samples_prove_theorem=False,
        local_support_disjointness_requires_N_smooth=False,
        one_chosen_meridian_per_image_assumed_to_generate=False,
        simultaneous_rectification_of_components_assumed=False,
        binary_affine_family_excluded=False,natural_alternating_symmetric_actions_excluded=False,
        all_imprimitive_actions_excluded=False,intermediate_sources_assumed_planes=False,
        all_listed_containers_realizable=False,actual_global_finite_algebra_closed=False,
        historical_priority_claimed=False,formal_verification=False,external_review=False,research_goal_completed=False),
      source_sha256={f:sha(P/f) for f in dict.fromkeys(files)},
      environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'disjoint_branch_supports_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',len(checks),'PF121 controls; subset/classical families excluded, primitive reduction; goal ACTIVE')


if __name__=='__main__':main()
