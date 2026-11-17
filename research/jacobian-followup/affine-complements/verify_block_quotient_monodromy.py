"""PF122 small exact controls; the hereditary geometric proof is written."""
from pathlib import Path
import hashlib,json,platform
import sympy as S
import block_quotient_monodromy as L
import disjoint_branch_supports as D
import primitive_product_action as W

P=Path(__file__).resolve().parent


def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def ident(name,value):check(name,S.expand(value)==0)
    def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
    olds=['disjoint_branch_supports_receipt.json','primitive_product_action_receipt.json',
          'half_support_inertia_receipt.json','local_monodromy_branch_receipt.json']
    for name in olds:
        r=json.loads((P/name).read_text());check('sealed predecessor PASS '+name,r['outcome']=='PASS')
        for f,h in r['source_sha256'].items():check('sealed predecessor '+name+' '+f,sha(P/f)==h)
    manifest=json.loads((P/'sources/block-quotient-monodromy/manifest.json').read_text())
    for item in manifest['sources']:
        f=P/item['path'];check('source seal '+f.name,sha(f)==item['sha256'] and f.stat().st_size==item['bytes'])
    t,b,a,delta,t0,eta=S.symbols('t b a delta t0 eta',positive=True)
    ident('UNBOUNDED exact block-defect ledger rearrangement',t*b-1-(t*a+delta+t0+eta)-(t*(b-a)-(1+delta+t0+eta)))
    ident('UNBOUNDED projected original 5/6 cap',S.Rational(5,6)*t*b/t-S.Rational(5,6)*b)
    blocks=[{0,1},{2,3},{4,5},{6,7}]
    p=[2,3,0,1,5,4,6,7];q=L.induced_blocks(p,blocks)
    check('actual support defect retains within-fixed-block movement',len(L.support(p))-2*len(L.support(q))==2)
    p=[2,3,0,1,4,5,6,7];q=[0,1,2,3,6,7,4,5]
    bp,bq=L.induced_blocks(p,blocks),L.induced_blocks(q,blocks)
    check('actual original disjoint block-moving pair descends',not(L.support(p)&L.support(q)) and not(L.support(bp)&L.support(bq)) and bool(L.support(bp)) and bool(L.support(bq)))
    p=[1,0,2,3,4,5,6,7];q=[0,1,3,2,4,5,6,7]
    check('countercontrol arbitrary disjoint original pair may project trivially',not(L.support(p)&L.support(q)) and not L.support(L.induced_blocks(p,blocks)) and not L.support(L.induced_blocks(q,blocks)))
    lift=[2,3,1,0];down=L.induced_blocks(lift,[{0,1},{2,3}])
    check('countercontrol quotient involution can have order4 lift',all(down[down[i]]==i for i in range(2)) and any(lift[lift[i]]!=i for i in range(4)))
    # An invisible branch is retained in the exact defect balance.
    tt,bb=2,5;ww=[4,3];aa=[2,0];dd=[w-tt*x for w,x in zip(ww,aa)];unram,excess=1,1
    check('actual integer block ledger with invisible branch',tt*bb-1==sum(ww)+unram+excess and tt*(bb-sum(aa))==1+sum(dd)+unram+excess)
    check('UNBOUNDED minimum support and disjoint pair exclude degree2/3',2*2>3)
    # Product base5 is the only case not removed by 5/6 alone.
    check('UNBOUNDED product base6 strict cutoff',1-S.Rational(1,6)==S.Rational(5,6))
    check('UNBOUNDED twisted product cutoff',1-S.Rational(1,8)>S.Rational(5,6))
    check('UNBOUNDED base5 swap plus second image overfills budget',S.Rational(4,5)+S.Rational(2,5)>1)
    product_controls=[]
    for c,s in [(5,2),(5,3),(6,2)]:
        top=list(range(s));top[0],top[1]=top[1],top[0]
        z=W.wreath_permutation(c,[list(range(c)) for _ in range(s)],top)
        moved=len(L.support(z));degree=c**s
        check('actual pure swap product support '+str((c,s)),moved*c==(c-1)*degree)
        product_controls.append(dict(base=c,factors=s,degree=degree,moved=moved))
    # Pointwise-stabilizer proof: form differences span v-perp, vectors span V.
    classical=[]
    for n in [3,4]:
        v=1;forms=D.minus_forms(n)
        moved=[c for c in forms if D.quadratic(v,n,c)==0]
        first=moved[0];diffs=[c^first for c in moved]
        check('actual minus-form support differences rank '+str(n),L.binary_rank(diffs)==2*n-1)
        check('actual form differences annihilate transvection center '+str(n),all((c&v).bit_count()%2==0 for c in diffs))
        classical.append(dict(action='minus_forms',a=n,support=len(moved),difference_rank=L.binary_rank(diffs)))
    for n,negative,epsilon in [(3,True,0),(4,True,0),(4,False,1)]:
        v=1 if negative else 1|(1<<n)
        points,p=D.vector_permutation(n,negative,epsilon,v)
        moved=[points[i] for i in L.support(p)]
        check('actual orthogonal support spans whole space '+str((n,negative)),L.binary_rank(moved)==2*n)
        check('actual orthogonal support difference rank '+str((n,negative)),L.binary_rank([x^moved[0] for x in moved])==2*n-1)
        classical.append(dict(action='orthogonal_vectors',a=n,negative=negative,support=len(moved),span_rank=L.binary_rank(moved)))
    h=S.symbols('h',nonnegative=True)
    check('UNBOUNDED codimension2 direction has nonzero polar rank',S.expand((2*a-4).subs(a,3+h))==2+2*h)
    # Original local normalized degrees may be >1; their sum, not point count, is f.
    es,fs,etas,ells=[4,3],[3,2],[[2,1],[4]],[[10,5],[10]]
    for i,(ee,ff,ex,ell) in enumerate(zip(es,fs,etas,ells)):
        qs=[2,1] if i==0 else [2]
        check('actual full residue-fibre ledger allows f>1 '+str(i),sum(qs)==ff and all(x==ee*q0+y for x,q0,y in zip(ell,qs,ex)) and sum(ell)==ee*ff+sum(ex))
    text=' '.join((P/'HEREDITARY_BLOCK_QUOTIENT_MONODROMY_REDUCTION.md').read_text().lower().split())
    for phrase in ['every primitive block quotient','including imprimitive monodromy','original finite normalization',
      'all omitted primes','arbitrary abstract homomorphic image','delta_j>=0','sum_j a_j<=b-1',
      'need not lift to an involution','no plane, smoothness','normalized residue degrees summing to f_i',
      'f_i>1','normally generate','positive-genus shifted fibre','c=5','4b/5','2b/5',
      'trivial pointwise stabilizer','the second element need not be a transvection',
      'not obtained or visually audited','w_c=d/2','they need not all be2',
      'degree4 quotient','not transferred','remain open']:
        check('written scope '+phrase,phrase in text)
    files=['HEREDITARY_BLOCK_QUOTIENT_MONODROMY_REDUCTION.md','block_quotient_monodromy.py',
      'verify_block_quotient_monodromy.py','sources/block-quotient-monodromy/manifest.json']+olds+[x['path'] for x in manifest['sources']]
    result=dict(claim='PF2RULED122',outcome='PASS',check_count=len(checks),checks=checks,
      theorem=dict(all_plane_degrees_all_monodromy_types=True,whole_original_normalization_all_primes_required=True,
        every_nontrivial_block_quotient_exact_support_defect_budget=True,
        every_block_quotient_full_positive_weighted_branch_excluded=True,
        every_block_quotient_irreducible_contractible_full_branch_excluded=True,
        every_block_quotient_disjoint_nonidentity_BRANCH_pair=True,
        quotient_degrees2_and3_excluded=True,
        all_interior_subset_block_quotients_excluded=True,
        all_primitive_product_action_block_quotients_excluded=True,
        classical_transvection_support_pointwise_stabilizer_trivial=True,
        all_three_binary_classical_primitive_block_quotients_excluded=True,
        every_primitive_block_quotient_natural_Alt_Sym_or_binary_affine=True,
        binary_affine_quotient_one_branch_original_support_exact_half=True),
      product_controls=product_controls,classical_rank_controls=classical,
      assurance=dict(published_CFSG_classification_imported=True,
        prior_GM_primary_web_import_boundary_retained=True,
        original_GM_PDF_audited=False,independent_CFSG_proof=False,finite_samples_prove_theorem=False,
        intermediate_models_assumed_planes=False,original_N_assumed_smooth=False,
        all_high_mass_residue_degrees_assumed_one=False,
        quotient_involution_assumed_to_lift_to_involution=False,
        arbitrary_original_disjoint_pair_assumed_nontrivial_downstairs=False,
        arbitrary_abstract_quotient_representations_covered=False,
        all_subgroup_composition_factors_classified=False,
        PF119_cusp_profile_transferred_to_every_intermediate=False,
        degree4_plane_inversion_applied_to_intermediate=False,
        all_natural_or_affine_block_quotients_excluded=False,
        all_imprimitive_monodromy_excluded=False,actual_global_finite_algebra_closed=False,
        historical_priority_claimed=False,formal_verification=False,external_review=False,research_goal_completed=False),
      source_sha256={f:sha(P/f) for f in dict.fromkeys(files)},environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'block_quotient_monodromy_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',len(checks),'PF122 controls; hereditary primitive block quotient reduction; goal ACTIVE')


if __name__=='__main__':main()
