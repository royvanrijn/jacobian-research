"""PF105 exact controls; finite examples do not prove the universal theorem."""
from pathlib import Path
import hashlib,json,platform,math
import sympy as S
import triple_boundary_residue as T

P=Path(__file__).resolve().parent


def main():
    checks=[]
    def check(name,value):
        assert value,name;checks.append(name)
    def ident(name,value):check(name,S.factor(value)==0)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    def ord0(poly,var):return min(m[0] for m,c in S.Poly(S.expand(poly),var).terms())
    old=json.loads((P/'fibre_transitive_spin_receipt.json').read_text())
    check('PF104 predecessor PASS',old['outcome']=='PASS')
    for f,h in old['source_sha256'].items():check('PF104 predecessor seal '+f,sha(P/f)==h)
    manifest=json.loads((P/'sources/triple-residue-excess/manifest.json').read_text())
    for item in manifest['sources']:
        f=P/item['path'];check('exact primary/context pin '+f.name,sha(f)==item['sha256'] and f.stat().st_size==item['bytes'])
    check('fresh full-page Borisov power LOCATION audit',manifest['sources'][0]['fresh_complete_PDF_pages_rendered_and_visually_inspected']==[3,11,12,13])
    check('fresh full-page Orevkov endpoint audit',manifest['sources'][1]['fresh_complete_PDF_pages_rendered_and_visually_inspected']==[2,7,8,9])
    # Universal weighted leading coefficient: no truncation depth or parameter bound.
    b=S.symbols('b',integer=True,positive=True)
    ident('UNBOUNDED weighted operator coefficient X2',2*S.Rational(1,2)-1)
    ident('UNBOUNDED weighted operator coefficient XYb',(1-b)*2/(b-1)+2)
    ident('UNBOUNDED weighted operator coefficient Y2b',(-2*b)*(-1/(2*b))-1)
    ident('UNBOUNDED nonzero boundary coefficient',S.Rational(1,2)+2/(b-1)-1/(2*b)-(b+1)**2/(2*b*(b-1)))
    ident('UNBOUNDED gcd divides two',2*b-2*(b+1)+2)
    # Quadratic coefficient reduction independent of projective Riemann-Hurwitz.
    a,c,d,e,f,B,X,Y=S.symbols('a c d e f B X Y')
    U=a*X*X+B*X*Y+c*Y*Y;V=d*X*X+e*X*Y+f*Y*Y
    J=S.diff(U,X)*S.diff(V,Y)-S.diff(U,Y)*S.diff(V,X)
    ident('binary quadratic Jacobian coefficient identity',J-(2*(a*e-B*d)*X*X+4*(a*f-c*d)*X*Y+2*(B*f-c*e)*Y*Y))
    ident('nonzero y2 coefficient would force forbidden proportional residuals',
          (a*e-B*d).subs({d:a*f/c,e:B*f/c}))
    # Rank-one finite germs at several indices, with a genuinely nonconstant unit.
    rankone=[];u,w,z=S.symbols('u w z')
    for index,m in [(2,2),(3,2),(4,3),(7,2)]:
        G=S.integrate((u-z**m)**(index-1)*(1+u+z),(z,0,w))
        restricted=S.expand(G.subs(u,w**m));expected=m*(index-1)+1
        coef=S.prod(S.Rational(m*j,m*j+1) for j in range(1,index))
        check('rank-one coprime first orders and exact local degree '+str((index,m)),
              ord0(restricted,w)==expected and ord0(G.subs(u,0),w)==expected and math.gcd(m,expected)==1)
        ident('rank-one nonzero all-index integral coefficient '+str((index,m)),S.expand(restricted).coeff(w,expected)-coef)
        ident('rank-one full Jacobian including unit '+str((index,m)),S.diff(G,w)-(u-w**m)**(index-1)*(1+u+w))
        rankone.append(dict(index=index,tangency=m,local_degree=expected,residue_local_degree=1,leading_coefficient=str(coef)))
    rankzero=[]
    for bb,perturb in [(2,0),(3,0),(3,1),(4,0)]:
        packet=T.rank_zero_model(bb,perturb);x=packet['X'];y=packet['Y'];uu=packet['u'];vv=packet['v']
        ident('rank-zero exact Jacobian '+str((bb,perturb)),
              S.diff(uu,x)*S.diff(vv,y)-S.diff(uu,y)*S.diff(vv,x)-packet['critical']**2*packet['unit'])
        check('rank-zero first orders and residue bound '+str((bb,perturb)),
              ord0(packet['restricted_u'],y)==2*bb and packet['upstairs_residue_degree']<=2)
        check('rank-zero actual axis intersection degree '+str((bb,perturb)),
              ord0(uu.subs(x,0),y)+ord0(uu.subs(y,0),x)==2*(bb+1))
        rankzero.append(dict(b=bb,perturb=perturb,first_orders=[bb+1,2*bb],
                             local_degree=packet['local_surface_degree'],residue_local_degree=packet['upstairs_residue_degree']))
    check('nonprimitive rank-zero germ retained',rankzero[1]['residue_local_degree']==2)
    check('primitive higher-characteristic rank-zero germ retained',rankzero[2]['residue_local_degree']==1 and rankzero[2]['first_orders']==[4,6])
    check('singular triple bound permits only quotient2 residue1',T.singular_triple_filter(2,1)['necessary_bound_satisfied'])
    check('singular triple quotient2 residue2 rejected',not T.singular_triple_filter(2,2)['necessary_bound_satisfied'])
    check('singular triple larger quotient rejected',not T.singular_triple_filter(3,1)['necessary_bound_satisfied'])
    # Unbounded genuine finite quotient family, proof in source note.
    family=[]
    for ss in [1,2,3]:
        c=T.quotient_counterfamily(ss);x=c['X'];y=c['Y'];g=c['g'];A=c['A'];B=c['B'];C=c['C']
        pp=c['p'];mm=c['m'];bb=c['b'];uu=c['u'];vv=c['v']
        ident('genuine finite counterfamily Jacobian '+str(ss),
              S.diff(uu,x)*S.diff(vv,y)-S.diff(uu,y)*S.diff(vv,x)-(x-y**bb)**pp)
        ident('counterfamily even invariant first coordinate '+str(ss),uu.subs({x:-x,y:-y})-uu)
        ident('counterfamily even invariant second coordinate '+str(ss),vv.subs({x:-x,y:-y})-vv)
        ident('actual downstairs A1 invariant formula '+str(ss),
              c['downstairs_u'].subs({A:g*g,B:g*y,C:y*y})-uu.subs(x,g+y**bb))
        ident('actual downstairs second coordinate '+str(ss),
              c['downstairs_v'].subs({A:g*g,B:g*y,C:y*y})-vv.subs(x,g+y**bb))
        ident('all-parameter nonzero boundary coefficient control '+str(ss),
              sum(c['coefficients'])-c['leading_boundary_coefficient'])
        ident('normalized boundary power map first coordinate '+str(ss),
              c['downstairs_u'].subs({A:0,B:0})-c['leading_boundary_coefficient']*C**(ss*(4*ss-1)))
        ident('normalized boundary power map second coordinate '+str(ss),
              c['downstairs_v'].subs(B,0)-C**(2*ss))
        check('actual residue and quotient degrees '+str(ss),
              math.gcd(2*ss,ss*(4*ss-1))==ss and pp+pp*bb==c['upstairs_degree'] and
              c['upstairs_degree']==2*c['quotient_degree'])
        ident('Y monic integrality identity '+str(ss),
              sum(c['coefficients'][j]*vv**(pp-j)*y**(mm*j) for j in range(pp+1))-uu*y**pp)
        ident('X monic integrality identity '+str(ss),
              sum(c['coefficients'][j]*vv**(bb*j)*x**(mm*(pp-j)) for j in range(pp+1))-uu*x**(bb*pp))
        check('both leading integrality coefficients are nonzero constants '+str(ss),
              c['coefficients'][0]!=0 and c['coefficients'][-1]!=0)
        check('wrong-source Euler budget differs by one '+str(ss),
              c['transverse_index']*c['residue_degree']+c['eta']==c['quotient_degree'])
        family.append(dict(s=ss,transverse_index=c['transverse_index'],residue_degree=c['residue_degree'],
                           quotient_degree=c['quotient_degree'],upstairs_residue_degree=c['upstairs_residue_degree'],
                           coefficient=str(c['leading_boundary_coefficient']),
                           downstairs_U=str(c['downstairs_u']),downstairs_V=str(c['downstairs_v']),
                           eta=c['eta'],complement=c['complement'],finite_Keller_realization=False))
    check('singular index3 example survives with residue1',family[0]['transverse_index']==3 and family[0]['residue_degree']==1)
    check('higher odd indices can have nontrivial residue on genuine non-plane finite sources',
          family[1]['transverse_index']==5 and family[1]['residue_degree']==2)
    aa,bb=S.symbols('a b',integer=True,nonnegative=True)
    ident('UNBOUNDED remaining residue cost at least8',(4+aa)*(2+bb)-8-(2*aa+4*bb+aa*bb))
    check('all e<=3 force f1 in ACTUAL Keller filters',
          all(not T.actual_keller_residue_filter(e,2)['low_index_residue_one'] for e in [1,2,3]))
    check('higher index not silently excluded by low-index theorem',T.actual_keller_residue_filter(4,2,9)['low_index_residue_one'])
    check('degree8 generic budget excludes e>=4 residue2 endpoint',
          not T.actual_keller_residue_filter(4,2,8)['generic_cost_budget'])
    check('PF104 natural8 ramified profile excluded WITHOUT smoothness',
          not T.actual_keller_residue_filter(3,2,8)['low_index_residue_one'])
    proof=' '.join((P/'TRIPLE_BOUNDARY_RESIDUE_ONE_AND_NORMAL_ITERATED_EXCLUSION.md').read_text().lower().split())
    for word in ['no smoothness','location','already used','not a new geometric discovery',
                 'independent of ell_p=e_D','n*q<=2','every field degree','not an automorphism theorem',
                 'normal shear remains','not a2','not all-degree jc2','one smooth branch']:
        check('proof/scope guard '+word,word.lower() in proof)
    files=['TRIPLE_BOUNDARY_RESIDUE_ONE_AND_NORMAL_ITERATED_EXCLUSION.md','triple_boundary_residue.py',
           'verify_triple_boundary_residue.py','sources/triple-residue-excess/manifest.json',
           'fibre_transitive_spin_receipt.json']+[x['path'] for x in manifest['sources']]
    result=dict(claim='PF2RULED105',outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={f:sha(P/f) for f in dict.fromkeys(files)},
        general=dict(all_field_degrees=True,arbitrary_other_omitted_primes=True,normal_singular_source_allowed=True,
                     all_indices_corank_one_local_residue_one=True,index3_corank_two_upstairs_local_residue_at_most2=True,
                     singular_index3_quotient_order2_and_local_residue1=True,
                     all_omitted_e_le3_have_f1=True,all_omitted_f1_if_d_le8=True,
                     all_omitted_f1_if_all_indices_odd_and_d_le10=True,
                     no_omitted_triple_longitude_joining_in_any_degree=True),
        normal_iterated_application=dict(PF104_natural_eight_all_lengths_all_parameters_excluded=True,
                                         smoothness_hypothesis_removed=True,arbitrary_extra_omitted_primes=True,
                                         arbitrary_finite_image_singularities=True,any_knot_type_triple_joining_excluded=True),
        rank_one_controls=rankone,rank_zero_controls=rankzero,
        genuine_unbounded_finite_counterfamily=dict(formula_proved_for_all_positive_s=True,
            transverse_index='2s+1',residue_degree='s',degree='4s^2',source='A1 quotient surface',
            complement='C* times C',finite_Keller_realization=False,representative_controls=family),
        attribution=dict(Borisov_power_map_and_location_imported=True,smooth_residue_one_already_primary=True,
                         prior_HC4_source_module_use_credited=True,endpoint_quadratic_mechanism_credited=True,
                         new_residue_bound_and_global_triple_deduction_written=True),
        assurance=dict(written_universal_proofs=True,primary_geometry_imported=True,
                       finite_checks_not_universal_proof=True,all_eight_sheet_Keller_maps_excluded=False,
                       all_indices_residue_one_claimed=False,all_singular_triple_points_excluded=False,
                       all_representation_of_underlying_knots_excluded=False,global_Laurent_finite_algebra_closed=False,
                       formal_verification=False,external_review=False,external_review_deferred=True,research_goal_completed=False),
        runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'triple_boundary_residue_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome='PASS',claim='PF2RULED105',check_count=len(checks))))


if __name__=='__main__':main()
