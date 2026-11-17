"""Small PF118 controls and retained-source seals; no group/subgroup campaign."""
from pathlib import Path
from itertools import product
import hashlib,json,platform
import sympy as S
import half_support_inertia as H

P=Path(__file__).resolve().parent


def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def ident(name,z):check(name,S.simplify(S.expand(z))==0)
    def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
    olds=['high_index_budget_receipt.json','local_monodromy_branch_receipt.json',
          'residual_milnor_receipt.json','nilradical_trace_receipt.json']
    for name in olds:
        r=json.loads((P/name).read_text());check('sealed predecessor PASS '+name,r['outcome']=='PASS')
        for f,h in r['source_sha256'].items():check('sealed predecessor '+name+' '+f,sha(P/f)==h)
    manifest=json.loads((P/'sources/half-support-inertia/manifest.json').read_text())
    for item in manifest['sources']:
        f=P/item['path'];check('primary/source seal '+f.name,sha(f)==item['sha256'] and f.stat().st_size==item['bytes'])
    check('primary HTML audit, no fresh PDF claim',manifest['primary_HTML_sections_read'] and not manifest['fresh_PDF_audit_claimed'])
    # Reuse the completed KMFI95 census rather than enumerate its 620 subgroups again.
    M=P/'../monodromy-theory';old=json.loads((M/'local_image_receipt.json').read_text())
    check('retained exact KMFI95 census PASS',old['result']=='pass')
    for f,h in old['sha256'].items():check('sealed KMFI95 '+f,sha(M/f)==h)
    model=json.loads((M/'local_image_model.json').read_text())
    rows={x['kind']:x for x in model['complete_conjugacy_classes']}
    fixed={2:3,3:2,5:1,6:0,11:0}
    for order,mark in fixed.items():check('both exceptional11 cyclic marks '+str(order),rows['C'+str(order)]['marks']==[mark,mark])
    check('exceptional11 minimum support8 excludes both plane actions',11-max(fixed.values())==8 and H.support_obstruction(11,8)['excluded_for_plane_Keller'])
    # These dimension identities encode the written unbounded inequalities.
    q,h,a,b=S.symbols('q h a b',positive=True)
    ident('UNBOUNDED projective half-support difference',(q*h+1)-2*(h+1)-((q-2)*h-1))
    ident('UNBOUNDED binary projective strict difference',(2*h+1)-2*h-1)
    ident('UNBOUNDED affine half-support difference',q*h-2*h-(q-2)*h)
    ident('UNBOUNDED eigenspace merging',q**(a+b-1)+q-q**a-q**b-q*(q**(a-1)-1)*(q**(b-1)-1))
    profiles=[]
    for qq,nn in [(2,2),(2,3),(2,5),(3,3),(4,2),(5,3),(8,3),(16,2),(2,100)]:
        z=H.projective_bound(qq,nn);check('projective all-parameter bound transcription '+str((qq,nn)),z['half_support_excluded']);profiles.append(z)
    check('projective degree4 bound equality is retained',not H.projective_bound(3,2)['half_support_excluded'] and H.projective_bound(3,2)['support_lower']==2)
    for qq,nn in [(3,1),(3,2),(4,2),(5,1),(8,3),(17,10)]:check('nonbinary affine bound '+str((qq,nn)),H.affine_bound(qq,nn)['half_support_excluded'])
    # Exact sharp eigenline controls over prime fields; no sampled theorem proof.
    for qq,nn in [(2,4),(3,2),(3,3),(5,3)]:
        vectors=[v for v in product(range(qq),repeat=nn) if any(v)]
        def normal(v):
            lead=next(x for x in v if x);iv=pow(lead,-1,qq)
            return tuple(iv*x%qq for x in v)
        points=sorted({normal(v) for v in vectors})
        def image(v):
            w=list(v)
            if qq==2:w[1]=(w[1]+w[0])%qq
            else:w[0]=2*w[0]%qq
            return normal(w)
        moved=sum(image(v)!=v for v in points)
        bound=H.projective_bound(qq,nn)
        check('actual sharp projective fixed-space matrix '+str((qq,nn)),len(points)==bound['degree'] and len(points)-moved==bound['fixed_upper'])
    # AGL5(2) equality control: a genuine permutation, not a Keller realization.
    n=5;d=2**n;perm=[z ^ (((z>>0)&1)<<1) for z in range(d)]
    cyc=sorted(map(len,H.cycles(perm)))
    check('actual binary transvection cycle structure',cyc==[1]*16+[2]*8)
    profile=H.binary_affine_profile(n)
    check('binary affine32 is NOT excluded by strict bound',not H.support_obstruction(d,16)['excluded_for_plane_Keller'] and not profile['binary_affine_family_excluded'])
    check('binary affine32 necessary full-source ledger',2*profile['ramified_prime_count']+profile['unramified_primes_plus_normalized_excess']==d-1)
    check('binary affine double-value full fibre consumes all sheets',2*profile['ramified_prime_count']*profile['branch_normalization_max_fibre']==d)
    check('binary profile retains smooth N, no unramified prime over C and empty source double fibres',profile['normalization_smooth'] and not profile['unramified_prime_over_branch_allowed'] and profile['source_fibre_at_double_values']==0)
    for degree,mu in [(6,2),(6,3),(7,3),(11,2)]:check('small-support/equality action is not falsely excluded '+str((degree,mu)),not H.support_obstruction(degree,mu)['excluded_for_plane_Keller'])
    text=' '.join((P/'HALF_SUPPORT_INERTIA_AND_PLANE_MONODROMY_EXCLUSIONS.md').read_text().lower().split())
    for phrase in ['original finite normalization','every omitted prime','full finite',
        'not a new derivation','both exceptional','higher source dimension','no prime-degree',
        'binary means','not excluded','source fibre is empty','contact order is unrestricted',
        'published quartic','global finite-algebra','no formal verification']:
        check('written scope '+phrase,phrase in text)
    files=['HALF_SUPPORT_INERTIA_AND_PLANE_MONODROMY_EXCLUSIONS.md','half_support_inertia.py',
        'verify_half_support_inertia.py','sources/half-support-inertia/manifest.json']+olds+[x['path'] for x in manifest['sources']]
    result=dict(claim='PF2RULED118',outcome='PASS',check_count=len(checks),checks=checks,
      theorem=dict(all_plane_field_degrees=True,whole_plane_all_omitted_primes_required=True,
        exact_global_inertia_support_budget=True,nontrivial_inertia_moves_at_most_half=True,
        large_minimal_degree_actions_excluded=True,both_exceptional_PSL2_11_plane_actions_excluded=True,
        natural_projective_linear_point_and_hyperplane_actions_excluded_in_stated_range=True,
        nonbinary_affine_linear_actions_excluded=True,binary_affine_half_support_profile=True,
        binary_profile_normalization_smooth=True,binary_profile_double_values_force_empty_source_fibre=True),
      exceptional11=dict(degree=11,cyclic_fixed_counts=fixed,minimum_support=8,plane_excluded=True),
      projective_controls=profiles,binary_affine32=profile,
      assurance=dict(retained_global_ledger_and_PF113_used=True,primary_Lin_Zaidenberg_import_retained=True,
        elementary_unbounded_linear_algebra_proofs=True,retained_KMFI95_marks_reused=True,
        subgroup_census_rerun=False,finite_samples_prove_theorem=False,
        arbitrary_affine_dimension_exclusion=False,all_semilinear_actions_excluded=False,
        binary_affine_family_excluded=False,all_degree6_maps_excluded=False,
        actual_plane_binary_realization_constructed=False,source_surjectivity_assumed=False,
        global_trace_moment_regularity_proved=False,actual_global_finite_algebra_closed=False,
        historical_priority_claimed=False,formal_verification=False,external_review=False,
        research_goal_completed=False),source_sha256={f:sha(P/f) for f in dict.fromkeys(files)},
      environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'half_support_inertia_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',len(checks),'PF118 controls; all-degree monodromy classes excluded; full goal ACTIVE')


if __name__=='__main__':main()
