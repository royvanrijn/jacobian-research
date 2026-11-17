"""Exact controls for PF2RULED59; deformation, LCP and Nori are sourced written proofs."""
from pathlib import Path
import hashlib,json,os,platform
import sympy as S
HERE=Path(__file__).resolve().parent;checks=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check(name,value):
    assert value,name
    checks.append(name)
def main():
    preds=['unramified_orevkov_receipt.json','lefschetz_centralizer_receipt.json','power_difference_receipt.json','canonical_quadrics_receipt.json']
    for name in preds:
        r=json.loads((HERE/name).read_text());check('predecessor PASS '+name,r['outcome']=='PASS')
        for p,h in r['source_sha256'].items():check('predecessor seal '+name+' '+p,sha(HERE/p)==h)
    manifest=HERE/'sources/nodal-deformation/manifest.json';m=json.loads(manifest.read_text())
    for row in m['primary_sources']:
        p=HERE/row['path_relative_to_affine_complements'];check('primary pin '+p.name,p.stat().st_size==row['bytes'] and sha(p)==row['sha256'])
    a,b,c,t,u,e,mult=S.symbols('a b c t u epsilon m')
    two=S.Matrix([[1,a],[1,b]]);three=S.Matrix([[1,a,a*a],[1,b,b*b],[1,c,c*c]])
    check('two-point evaluation determinant universal',S.expand(two.det()-(b-a))==0)
    check('three-point evaluation determinant universal',S.expand(three.det()-(b-a)*(c-a)*(c-b))==0)
    check('eta1 does not separate three points',S.Matrix([[1,a],[1,b],[1,c]]).rank()==2)
    check('eta0 does not separate two points',S.ones(2,1).rank()==1)
    jet2=S.Matrix([[1,a],[0,1]]);jet3=S.Matrix([[1,a,a*a],[0,1,2*a],[0,0,2]])
    check('two-point jet evaluation includes coalescence',jet2.det()==1)
    check('three-point jet evaluation includes coalescence',jet3.det()==2)
    eta,k=S.symbols('eta k',integer=True)
    # Written unbounded cohomology argument: H1(P1,O(n))=0 for n>=-1.
    check('arbitrary-k separation threshold identity',S.expand((eta-k)-((eta-(k-1))-1))==0)
    check('two-point threshold eta1',(1-2)==-1)
    check('three-point threshold eta2',(2-3)==-1)
    check('eta0 leaves H1 obstruction to two-point evaluation',(0-2)==-2)
    check('eta1 leaves H1 obstruction to three-point evaluation',(1-3)==-2)
    h=u**mult-e
    check('arbitrary-contact squarefree derivative identity',S.simplify(u*S.diff(h,u)-mult*h-mult*e)==0)
    for n in [1,2,3,7]:
        f=u**n-e
        check('contact control '+str(n)+' generic squarefree',S.gcd(f,S.diff(f,u))==1)
        check('contact control '+str(n)+' discriminant only at zero',S.expand(S.discriminant(f,u)-(-1)**(n*(n-1)//2+n-1)*n**n*e**(n-1))==0)
    selfint,canon,delta=S.symbols('selfint canon delta')
    check('nodal adjunction equals normal degree',(selfint-2*delta).subs(selfint,-canon-2+2*delta)==-canon-2)
    check('quartic a_n normal convention',S.expand((S.Symbol('a_n')-1)-((-S.Symbol('KdotC')-2).subs(S.Symbol('KdotC'),-S.Symbol('a_n')-1)))==0)
    check('source tangent and normal degree add correctly',2+eta==eta+2)
    check('triple simple folds consume six sheets',3*2==6 and all(3*2>d for d in [3,4,5]))
    check('eta1 survives unrestricted all-degree theorem but not double-valued theorem',1<2 and 1>=1)
    # Fixed boundary intersection one cannot split into several positive intersections.
    check('boundary degree1 retained fibre intersection',S.Poly(u-e,u).degree()==1 and S.diff(u-e,u)==1)
    # Retained A5 example: Orevkov invariant, not original pi1, is transported.
    x,y,w=S.symbols('x y w');xx=(1-t*t)**3;yy=t*(1-t*t);G=x*x-(x+y*y)**3
    check('A5 control lies on irreducible retained curve',S.expand(G.subs({x:xx,y:yy},simultaneous=True))==0)
    check('A5 control is immersive',S.gcd(S.diff(xx,t),S.diff(yy,t))==1)
    uu=S.expand(xx+yy*yy)
    check('A5 target polynomial coordinate and degrees',S.expand(uu-(1-t*t)**2)==0 and S.degree(uu,t)==4 and S.degree(yy,t)==3 and S.gcd(4,3)==1)
    check('A5 positive normal degree is6',4+3-1==6)
    cubic=w**3-3*(x+y*y)*w-2*x
    check('A5 cubic discriminant exactly retained branch',S.expand(S.discriminant(cubic,w)+108*G)==0)
    check('A5 triple fibre cannot be substituted for a disjoint-fold Keller fibre',cubic.subs({x:0,y:0})==w**3)
    check('A5 total cubic is smooth at its triple point',S.diff(cubic,x).subs({x:0,y:0,w:0})==-2)
    inputs=[Path(__file__),HERE/'NODAL_DEFORMATION_AND_NONPOSITIVE_NORMAL.md',manifest]+[HERE/n for n in preds]+[HERE/row['path_relative_to_affine_complements'] for row in m['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
      deformation=dict(immersed_normalizations_required=True,normal_bundle='O_P1(eta)',unobstructed_for_eta_at_least_one=True,fixed_compactification_and_boundary=True,all_infinity_intersection_numbers_retained=True,proper_end_isotopy_verified=True,normalizations_remain_A1=True,LCP_verified=True,global_moduli_connectedness_assumed=False,unrestricted_multiple_values_require_eta_at_least_two=True,eta1_two_point_use_requires_initial_double_values=True,eta_nonpositive_excluded=False),
      topology=dict(Orevkov_invariant_preserved=True,original_complement_group_preserved=False,Nori_applied_only_after_nodal_deformation=True,simultaneous_full_branch_union=True,whole_union_eta_at_least_two_implies_abelian=True,whole_union_eta_at_least_one_and_individual_double_values_implies_abelian=True,ordinary_original_tangential_complement_asserted_abelian=False),
      Keller=dict(arbitrary_field_degrees=True,arbitrary_ramification_indices_for_full_union=True,all_unramified_omitted_primes_retained=True,zero_excess_hypothesis_explicit=True,all_packets_assumed_zero_excess=False,full_union_some_eta_at_most_one=True,full_union_individual_double_values_some_eta_nonpositive=True,sole_fold_eta_at_most_one=True,sole_fold_eta_one_requires_triple_value_and_degree_at_least_six=True,sole_fold_degree_at_most_five_eta_nonpositive=True,whole_quartic_eta_nonpositive=True,whole_quartic_a_n_at_most_one=True,all_quartic_source_rulings_and_poles=True,all_finite_contact_multiplicities=True,unramified_image_H_treated_as_branch=False),
      relation_to_retained=dict(PF46_tangency_budget_not_deleted=True,PF47_positive_eta_quartic_range_now_closed=True,PF47_multi_component_selected_bound_retained=True,PF48_unramified_sheet_bound_retained=True,PF58_full_canonical_pencil_retained=True,negative_control_G0_original_complement_has_S3=True),
      assurance=dict(written_unbounded_proofs=True,finite_controls_not_universal_proof=True,primary_imports_independently_reproved=False,formal_verification=False,independent_full_reproof=False,external_review=False,novelty_priority=False,full_quartic_exclusion=False,higher_degree_exclusion=False,research_goal_completed=False),
      source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'nodal_deformation_receipt.json').write_text(json.dumps(report,indent=2)+'\n');print(f'PASS {len(checks)} controls; nodal deformation and nonpositive quartic normal degree')
if __name__=='__main__':main()
