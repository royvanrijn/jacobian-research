"""PF120 exact small controls and source pins; published CFSG theorem is imported."""
from pathlib import Path
from fractions import Fraction
import hashlib,json,platform
import sympy as S
import primitive_product_action as L

P=Path(__file__).resolve().parent


def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def ident(name,z):check(name,S.simplify(z)==0)
    def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
    olds=['torus_cusp_excess_receipt.json','local_monodromy_branch_receipt.json',
          'half_support_inertia_receipt.json','simple_ramification_ledger_receipt.json']
    for name in olds:
        r=json.loads((P/name).read_text());check('sealed predecessor PASS '+name,r['outcome']=='PASS')
        for f,h in r['source_sha256'].items():check('sealed predecessor '+name+' '+f,sha(P/f)==h)
    manifest=json.loads((P/'sources/primitive-product-action/manifest.json').read_text())
    for item in manifest['sources']:
        f=P/item['path'];check('source seal '+f.name,sha(f)==item['sha256'] and f.stat().st_size==item['bytes'])
    check('classification import and full-page audit declared',manifest['CFSG_classification_imported'] and manifest['author_PDF_pages_rendered_and_read']==[4,5] and not manifest['independent_CFSG_proof_claimed'])
    b,d=S.symbols('b d',positive=True)
    ident('UNBOUNDED twisted fixed-ratio maximum',b**2-8*(b-2)-(b-4)**2)
    check('UNBOUNDED 7/8 exceeds 5/6',S.Rational(7,8)>S.Rational(5,6))
    check('UNBOUNDED pure swap minimum equals involution cutoff',1-S.Rational(1,3)==S.Rational(2,3))
    examples=[]
    for bb,top,twist,fixed,invol in [(3,[1,0],None,3,True),(4,[1,0],0,2,False),
          (5,[1,0,2],None,25,True),(3,[1,2,0],None,3,False),
          (3,[1,0,3,2],None,9,True),(3,[1,0,2],2,3,True)]:
        base=[list(range(bb)) for _ in top]
        if twist is not None:base[twist][0],base[twist][1]=base[twist][1],base[twist][0]
        p=L.wreath_permutation(bb,base,top);deg=len(p);fix=sum(i==x for i,x in enumerate(p));order2=all(p[p[i]]==i for i in range(deg))
        check('actual Cartesian fixed count '+str((bb,top,twist)),fix==fixed and order2==invol)
        ratio=Fraction(fix,deg)
        check('actual uniform top-action dichotomy '+str((bb,top,twist)),(order2 and ratio==Fraction(1,bb)) or ratio<=Fraction(1,8))
        check('actual top-inertia excluded by correct strict bound '+str((bb,top,twist)),Fraction(deg-fix,deg)>=Fraction(5,6) or order2 and Fraction(deg-fix,deg)>=Fraction(2,3))
        examples.append(dict(base_degree=bb,factors=len(top),degree=deg,fixed=fix,involution=order2))
    A,B,C=S.symbols('A B C',positive=True);rows=L.classical_table(A,B,C,S)
    expected={'L_field2_points':1,'L_field3_points':(A-3)/2,'U4_field3':32,
      'Sp_field2_points':1,'Sp_plus_forms':A,'Sp_minus_forms':-A,
      'odd_O_field3_P1':(3*B-1)*(B-1)/2,'odd_O_field3_Nminus':(3*B**2-5*B-4)/2,
      'Oplus_field2_singular':A+1,'Ominus_field2_singular':-(A-1),
      'Ominus_field3_P1':(C-1)**2/2,'Oplus_field2_nonsingular':-A,
      'Ominus_field2_nonsingular':A,'Oplus_field3_N1':C*(C-3)/2,
      'Ominus_field3_N1':(C**2-C-4)/2}
    for key,(degree,mu) in rows.items():ident('UNBOUNDED Table2 2mu-degree '+key,2*mu-degree-expected[key])
    k=S.symbols('k',integer=True,nonnegative=True)
    for key,parameter,lower in [('L_field3_points',A,9),('odd_O_field3_P1',B,3),
        ('odd_O_field3_Nminus',B,3),('Ominus_field3_P1',C,3),
        ('Oplus_field3_N1',C,9),('Ominus_field3_N1',C,3)]:
        poly=S.Poly(S.expand(expected[key].subs(parameter,lower+k)),k)
        check('UNBOUNDED positive discarded row '+key,all(c>=0 for c in poly.all_coeffs()) and poly.eval(0)>0)
    check('small U4(2).2 is retained degree27 orthogonal row',rows['Ominus_field2_singular'][0].subs(A,4)==27 and rows['Ominus_field2_singular'][1].subs(A,4)==12)
    for m,kk in [(5,2),(7,3),(10,3),(12,4)]:
        z=L.binomial_composite_certificate(m,kk);check('exact interior binomial composite '+str((m,kk)),z['factor']*z['cofactor']==z['degree'])
    ident('UNBOUNDED nonalternating lower band even rows',7*(A-1)-3*(2*A-1)-(A-4))
    ident('UNBOUNDED nonalternating lower band odd row',7*A-3*(2*A+1)-(A-3))
    ident('UNBOUNDED high mass cost contradicts complete ledger',S.Rational(6,5)*d/2+S.Rational(3,7)*d-d-d/35)
    check('no third branch image can fit',3*S.Rational(3,7)>1)
    text=' '.join((P/'PRIMITIVE_PRODUCT_ACTION_EXCLUSION_AND_PRIME_DEGREE_REDUCTION.md').read_text().lower().split())
    for phrase in ['original finite normalization','all omitted primes','every affine divisorial inertia',
        'normally generate','not assumed one chosen meridian','no cfgs',
        'classification is not proved','triality','nonsimple socle','necessary containers',
        'prime generic degree','nonalternating socle','all branch meridians',
        'ordinary imprimitive wreath action','no intermediate field','compact exceptional',
        'remain open']:
        # Keep spelling controls independent of typography and accepted CFSG ordering.
        if phrase=='no cfgs':phrase='no cfsg'
        check('written scope '+phrase,phrase in text)
    files=['PRIMITIVE_PRODUCT_ACTION_EXCLUSION_AND_PRIME_DEGREE_REDUCTION.md','primitive_product_action.py',
      'verify_primitive_product_action.py','sources/primitive-product-action/manifest.json']+olds+[x['path'] for x in manifest['sources']]
    result=dict(claim='PF2RULED120',outcome='PASS',check_count=len(checks),checks=checks,
      theorem=dict(all_plane_degrees=True,whole_original_normalization_all_primes_required=True,
        every_inertia_support_strictly_less_five_sixths=True,involutory_inertia_support_strictly_less_two_thirds=True,
        all_primitive_Cartesian_product_actions_base_at_least3_excluded=True,
        product_action_exclusion_CFSG_free=True,complete_primitive_necessary_container_reduction=True,
        all_prime_degrees_force_natural_alternating_or_symmetric=True,
        nonalternating_primitive_at_most_two_branch_images=True,nonalternating_every_inertia_support_at_most_half=True),
      actual_wreath_controls=examples,
      assurance=dict(published_CFSG_classification_imported=True,full_PDF_table_visually_audited=True,
        finite_samples_prove_theorem=False,independent_CFSG_classification_proved=False,
        all_wreath_actions_including_imprimitive_excluded=False,all_intermediate_models_assumed_planes=False,
        one_meridian_per_image_assumed_to_generate=False,all_listed_containers_realizable=False,
        all_prime_degree_maps_invertible=False,binary_classical_families_excluded=False,
        alternating_subset_families_excluded=False,binary_affine_family_excluded=False,
        target_blowup_affine_boundary_hypotheses_preserved=False,
        actual_global_finite_algebra_closed=False,historical_priority_claimed=False,
        formal_verification=False,external_review=False,research_goal_completed=False),
      source_sha256={f:sha(P/f) for f in dict.fromkeys(files)},environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'primitive_product_action_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',len(checks),'PF120 controls; primitive product actions excluded and prime-degree reduction; goal ACTIVE')


if __name__=='__main__':main()
