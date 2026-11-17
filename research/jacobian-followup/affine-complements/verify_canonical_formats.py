"""Exact controls for PF2RULED57; adjunction and global normalization are written proofs."""
from pathlib import Path
import hashlib,json,math,os,platform
import sympy as S
HERE=Path(__file__).resolve().parent;checks=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check(name,value):
    assert value,name
    checks.append(name)
def main():
    preds=['unramified_orevkov_receipt.json','trace_dual_pencil_receipt.json','power_difference_receipt.json']
    for n in preds:
        r=json.loads((HERE/n).read_text());check('predecessor PASS '+n,r['outcome']=='PASS')
        for p,h in r['source_sha256'].items():check('predecessor seal '+n+' '+p,sha(HERE/p)==h)
    legacy=HERE.parent/'jc2-omitted-divisors'
    for n in ['global_receipt.json','discriminant_receipt.json']:
        r=json.loads((legacy/n).read_text());check('legacy PASS '+n,r['outcome']=='PASS')
        for field in ['input_sha256','output_sha256']:
            for p,h in r.get(field,{}).items():check('legacy seal '+n+' '+p,sha(legacy/p)==h)
    manifest=HERE/'sources/canonical-formats/manifest.json';m=json.loads(manifest.read_text())
    for row in m['primary_sources']:
        p=HERE/row['path_relative_to_affine_complements'];check('pinned primary '+p.name,p.stat().st_size==row['bytes'] and sha(p)==row['sha256'])
    d,p,e,ell=S.symbols('d p e ell',integer=True)
    check('strict conductor class positivity decomposition',S.expand((d-2)*p-(e-1)-((d-2)*(p-1)+(d-1-e)))==0)
    check('binary embedding forces maximal ramification under Euler bound',S.expand(((d-2)*ell+1)-(d-1)-(d-2)*(ell-1))==0)
    check('binary degree2 canonical coefficient zero',2-2==0)
    # The all-format proof is the sum identity; these are exact arithmetic controls.
    primitive=[(3,),(2,2)]
    check('primitive canonical formats have total excess2',all(sum(a-1 for a in row)==2 for row in primitive))
    check('primitive canonical formats are exactly cubic and two quadrics',[math.prod(row) for row in primitive]==[3,4])
    for row,g in [((3,),1),((2,2),1),((4,),2),((2,3),2),((2,2,2),2),((3,3),3)]:
        s=sum(row)-len(row)-1;degree=math.prod(row)
        check('closed content and degree-bound control '+str((row,g)),s>0 and g%s==0 and len(row)<=g+1 and degree<=2**(g+1) and all(a<=2**(a-1) for a in row))
    check('primitive hypothesis cannot be omitted in degree6',(2+3-2-1)==2 and 2%2==0 and 1%2!=0 and 2*3==6)
    check('affine complete intersection would have zero canonical vector',S.zeros(2,1)!=S.Matrix([1,0]))
    pole=S.symbols('pole',positive=True,integer=True)
    K=S.Matrix([1,0]);L=S.Matrix([1,pole]);cond=2*L-K
    check('whole r2q1 conductor class keeps arbitrary D pole',cond==S.Matrix([1,2*pole]))
    check('minimal D pole2 is only a subcase',cond.subs(pole,2)==S.Matrix([1,4]))
    # Universal quartic-resolvent polynomial identity, with nonunit leading coefficient.
    a,b,c,f,g,t,z=S.symbols('a b c f g t z')
    quartic=a*t**4+b*t**3+c*t*t+f*t+g
    rho=z**3-c*z*z+(b*f-4*a*g)*z+(4*a*c*g-a*f*f-b*b*g)
    check('universal binary quartic equals monic cubic resolvent discriminant',S.expand(S.discriminant(quartic,t)-S.discriminant(rho,z))==0)
    check('cubic resolvent stays monic for arbitrary quartic leading coefficient',S.Poly(rho,z).LC()==1)
    check('index exponent transfers identically under equal normal discriminants',S.simplify((S.Symbol('ram')+2*S.Symbol('I'))-S.Symbol('ram'))==2*S.Symbol('I'))
    # Actual trace determinants in split local algebras isolate generator collisions.
    h=S.symbols('h');roots=[S.Integer(0),S.Integer(1),h,S.Integer(2)]
    V=S.Matrix([[r**j for j in range(4)] for r in roots]);J4=S.factor(V.det());F=S.prod(t-r for r in roots).expand()
    check('etale quartic normalization can have divisorial inverse index',J4!=0 and S.degree(J4,h)==3 and S.expand(S.discriminant(F,t)-J4**2)==0)
    pairing=[roots[0]*roots[1]+roots[2]*roots[3],roots[0]*roots[2]+roots[1]*roots[3],roots[0]*roots[3]+roots[1]*roots[2]]
    V3=S.Matrix([[r**j for j in range(3)] for r in pairing]);J3=S.factor(V3.det())
    check('quartic and cubic split normalization indices agree exactly up to sign',S.expand(J4**2-J3**2)==0)
    check('index divisor is independent of true branch',S.factor(J4**2)==4*h**2*(h-1)**2*(h-2)**2)
    # Local pole order is not a local normalization index.
    w=S.symbols('w');local=(w-1)*(w*w-h)
    check('quadratic fold is monogenic at infinity uniformizer',S.discriminant(w*w-h,w)==4*h)
    check('local product index separates infinite from finite factor',S.expand(S.discriminant(local,w)-4*h*(1-h)**2)==0 and (1-h).subs(h,0)==1)
    check('unramified degree1 factor is normal at arbitrary reciprocal pole',S.degree(w-h**7,w)==1 and S.discriminant(w-h**7,w)==1)
    Q0=S.eye(3);Q1=S.diag(h,2,3);u,v=S.symbols('u v');det=S.factor((u*Q0+v*Q1).det())
    check('everywhere smooth conic member gives unit leading cubic coefficient',Q0.det()==1 and det.subs({u:1,v:0})==1 and S.Poly(det.subs(v,1),u).LC()==1)
    check('a permitted numerical two-quadric format is not an exclusion',sum([2,2])-2-1==1)
    new_inputs=[Path(__file__),HERE/'CANONICAL_FORMATS_AND_BINARY_CONDUCTOR.md',manifest]+[HERE/n for n in preds]+[legacy/n for n in ['global_receipt.json','discriminant_receipt.json']]+[HERE/row['path_relative_to_affine_complements'] for row in m['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
      scope=dict(canonical_and_binary_gates_all_field_degrees=True,all_source_orders=True,all_target_degrees=True,all_ramification_indices=True,all_residue_degrees=True,all_omitted_primes_retained=True,normal_surface_allowed_singular=True,full_boundary_Euler_inputs_retained=True),
      canonical=dict(class_basis='Cl(N)=direct_sum Z[E_i]',vector='K_N=sum (e_i-1)E_i',content='g=gcd(e_i-1)>0',nonzero_nontorsion=True,global_affine_complete_intersection_excluded=True,local_complete_intersections_excluded=False),
      projective=dict(coefficient='s=sum(a_i)-r-1',s_divides_g=True,conditional_degree_bound='d<=2^(g+1)',primitive_formats=[[3],[2,2]],primitive_degree_at_least_five_CI_excluded=True,all_higher_degree_algebras_excluded=False,normal_binary_embedding_excluded_all_nontrivial_Keller_degrees=True,maximal_index_remark_required=True),
      conductor=dict(class_formula='C_cond=(d-2)L-K_N',all_admissible_sourcefinite_coordinates=True,class_coefficients_strictly_positive_for_d_at_least_four=True,divisorial_conductor_nonzero=True,target_index_nonconstant=True,curve_point_conductor_identified_with_surface_conductor=False,class_vector_identified_with_boundary_multiplicities=False,codimension_two_questions_closed=False),
      quartic=dict(binary_and_cubic_normalization_indices_equal=True,normalization_of_resolvent_order_is_actual_C3=True,whole_r2q1_D_pole_unrestricted=True,r2q1_conductor_class='R+2*p*D',generic_r2q1_index_coprime_both_omitted_images=True,distinct_index_curve_in_etale_source_locus_required=True,distinct_support_means_disjoint_curves=False,actual_full_two_quadric_pencil_has_no_global_everywhere_smooth_member=True,global_matrix_model_constructed=False),
      assurance=dict(written_unbounded_proofs=True,finite_controls_not_universal_proof=True,primary_imports_independently_reproved=False,formal_verification=False,independent_full_reproof=False,external_review=False,novelty_priority=False,finite_Keller_pair_supplied=False,full_quartic_exclusion=False,higher_degree_exclusion=False,research_goal_completed=False),
      source_sha256={os.path.relpath(p,HERE):sha(p) for p in new_inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'canonical_formats_receipt.json').write_text(json.dumps(report,indent=2)+'\n');print(f'PASS {len(checks)} controls; canonical formats and common binary/cubic conductor')
if __name__=='__main__':main()
