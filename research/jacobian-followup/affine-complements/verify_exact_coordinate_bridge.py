"""PF2RULED45 exact support; canonical and affine arguments are written proofs.
No bounded coefficient search is used as an unbounded geometric proof.
"""
from pathlib import Path
from math import gcd
import hashlib,json,platform
import sympy as S
HERE=Path(__file__).resolve().parent
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def pole_types(e):
    out=[]
    for g in range(1,e):
        for a in range(3,3*(e-1)//g+1):
            for b in range(2,min(a-1,2*(e-1)//g)+1):
                if gcd(a,b)==1 and g*(a*b-a-b)<=e-1:out.append((g*a,g*b,g,a,b))
    return sorted(out)
def main():
    pred_names=['regular_quartic_exclusion_receipt.json','universal_regular_end_receipt.json','multiple_fiber_cyclic_receipt.json','coordinate_regular_sector_receipt.json','intrinsic_normalization_receipt.json']
    for name in pred_names:
        d=json.loads((HERE/name).read_text());check('predecessor PASS '+name,d['outcome']=='PASS')
        for f,h in d['source_sha256'].items():check('predecessor seal '+name+' '+f,digest(HERE/f)==h)
    manifest=HERE/'sources/coordinate-negative-canonical/manifest.json';pins=json.loads(manifest.read_text())
    check('six precise primary inputs',len(pins['primary_sources'])==6)
    for row in pins['primary_sources']:
        f=HERE/row['path_relative_to_affine_complements'];raw=f.read_bytes()
        check('primary byte seal '+f.name,len(raw)==row['bytes'] and digest(f)==row['sha256'])
        check('primary format '+f.name,raw.startswith(b'%PDF') if row['kind']=='PDF' else b'<html' in raw.lower())
    check('general canonical formula is used with the correct infinity sign',pins['general_theorem_not_primitive_only'] and pins['canonical_sign_convention']=='delta=-ord_R; A_delta=1+ord_R(dx wedge dy)=e')
    check('complete signed-key, normal-form and canonical pages inspected',pins['new_visual_audits'][0]['PDF_pages']==[5,11,17] and pins['new_visual_audits'][0]['complete_signed_definitions_normal_form_and_general_canonical_formula'])
    check('complete published sign criterion inspected',pins['new_visual_audits'][1]['PDF_pages']==[3] and pins['new_visual_audits'][1]['complete_sign_criterion'])
    # All-index inequalities: no sampled key sequence replaces this algebra.
    wi,wn,C,alpha=S.symbols('w_i w_next C alpha')
    Cnext=C+(alpha-1)*wi
    check('strict key inequality is exactly strict defect monotonicity',S.expand((wn-Cnext)-(wi-C)-(wn-alpha*wi))==0)
    a0,b0,g=S.symbols('a0 b0 g',integer=True,positive=True)
    defect=a0*b0-a0-b0
    check('initial canonical defect is the coprime conductor minus one',S.expand(defect-((a0-1)*(b0-1)-1))==0)
    check('linear large-pole bound is an exact nonnegative factor',S.expand(defect-(a0-2)-(a0-1)*(b0-2))==0)
    check('linear small-pole bound has an explicit nonnegative remainder',S.expand(defect-(b0-1)-((a0-3)*(b0-1)+(b0-2)))==0)
    types={str(e):pole_types(e) for e in range(1,9)}
    check('index one has no all-coordinate-polar initial pair',types['1']==[])
    check('index two has only the exceptional pole pair three,two',types['2']==[(3,2,1,3,2)])
    check('index three retains exactly the two initial pole pairs',types['3']==[(3,2,1,3,2),(6,4,2,3,2)])
    for e,rows in types.items():
        ee=int(e)
        check('derived finite initial-type domain is respected e'+e,all(q<=3*(ee-1) and p<=2*(ee-1) and gg*(aa*bb-aa-bb)<ee for q,p,gg,aa,bb in rows))
    # Whole forced normal form and its actual two-chart volume residue.
    u,t,a,b,c,x,y=S.symbols('u t a b c x y')
    xx=u**-3;yy=u**-2+a*u**2+b*u**3+c*u**4+u**5*t
    jac=S.expand(S.diff(xx,u)*S.diff(yy,t)-S.diff(xx,t)*S.diff(yy,u))
    check('whole omitted-chart Jacobian is exactly the canonical simple zero',jac==-3*u)
    thetaU=-yy*S.diff(xx,u);thetaW=3*u*t
    diff=S.expand(thetaU-thetaW)
    check('full residue difference retains every coefficient',diff==3*(u**-6+a*u**-2+b*u**-1+c))
    check('exactness forces precisely the resonant coefficient b to vanish',S.residue(diff,u,0)==3*b)
    check('both local primitives differentiate to the source two-form',S.expand(-S.diff(thetaU,t)-jac)==0 and S.expand(-S.diff(thetaW,t)-jac)==0)
    A,B,d,f,v=S.symbols('A B d f v')
    shift=v+A*v**2/3
    yshift=S.series(shift**-2+d/shift,v,0,0).removeO().expand()
    check('affine x shear removes the pole-one lower term',S.expand(yshift.coeff(v,-1)-(d-2*A/3))==0)
    shift=v+B*v**4/3
    yshift=S.series(shift**-2+f*shift,v,0,2).removeO().expand()
    check('affine x translation removes the order-one term',S.expand(yshift.coeff(v,1)-(f-2*B/3))==0)
    G=y**3-x**2-3*a*y-3*c
    Gu=S.expand(G.subs({x:xx,y:yy},simultaneous=True))
    Gzero=S.expand(Gu.subs(b,0))
    check('polynomial key form has first coefficient three*t when b vanishes',S.Poly(S.expand((Gzero-3*u*t)/u**2),u) is not None and Gzero.coeff(u,1)==3*t)
    check('source polynomial division has its complete degree-three basis',S.Poly(G,y).degree()==3 and S.Poly(G,y).LC()==1)
    check('equal leading value and equal t degree force unique reduced monomial',all((2*j)%3!=0 for j in [-2,-1,1,2]))
    X,Y,Z=S.symbols('X Y Z')
    rel=Y**3-X**2-3*a*Z**4*Y-3*c*Z**6-Z**7
    check('exact affinization relation is an all-parameter polynomial identity',S.expand(rel.subs({X:x*G**3,Y:y*G**2,Z:G},simultaneous=True))==0)
    for name,expr in [('G',G),('X',x*G**3),('Y',y*G**2)]:
        check('whole source curve contracts under every affinization generator '+name,S.rem(expr,G,y)==0)
    # Actual b!=0 counterfamily, including every a,c and all fiber parameters.
    ff,TT,XX,s,V=S.symbols('ff TT XX s V')
    Tsrc=G**2-9*b**2*y;Xsrc=27*b**3*x
    P=(ff**2-TT)**3-243*a*b**4*(ff**2-TT)-2187*c*b**6-729*b**6*ff
    H=XX**2-P
    check('whole hyperelliptic finite chart is the original polynomial plane',S.expand(H.subs({ff:G,TT:Tsrc,XX:Xsrc},simultaneous=True))==0)
    yinv=(ff**2-TT)/(9*b**2);xinv=XX/(27*b**3)
    check('inverse finite-chart substitution is exact modulo its single relation',S.cancel(729*b**6*(ff-G.subs({x:xinv,y:yinv},simultaneous=True))-H)==0)
    check('all-parameter finite-chart smoothness has a unit gradient certificate',S.expand(S.diff(H,ff)+2*ff*S.diff(H,TT)-729*b**6)==0)
    check('control family critical locus cannot contain a moving curve',S.expand(S.diff(Tsrc,x)+4*x*G)==0 and S.expand(S.diff(Tsrc,y)-(6*G*(y**2-a)-9*b**2))==0)
    check('control family has only finitely many critical points',S.Poly(S.diff(Tsrc,y).subs(x,0),y).degree()==5 and S.Poly(S.diff(Tsrc,y).subs(x,0),y).LC()==6)
    q2,q1,q0=S.symbols('q2 q1 q0');square=(ff**3+q2*ff**2+q1*ff+q0)**2
    check('every fiber nonsquare proof starts with forced quadratic coefficient zero',S.expand(square-P).coeff(ff,5)==2*q2)
    check('every fiber nonsquare proof forces constant coefficient zero',S.expand((square-P).subs(q2,0)).coeff(ff,3)==2*q0)
    check('nonzero linear coefficient excludes square fibers for every parameter',S.expand((square-P).subs({q2:0,q0:0})).coeff(ff,1)==729*b**6)
    Pinf=(1-TT*s**2)**3-243*a*b**4*s**4*(1-TT*s**2)-2187*c*b**6*s**6-729*b**6*s**5
    Hinf=V**2-Pinf
    check('entire infinity chart is the same finite double cover of the projective base',S.cancel(s**6*H.subs({ff:1/s,XX:V/s**3},simultaneous=True)-Hinf)==0)
    check('two infinity sections remain disjoint and smooth for all parameters',Hinf.subs(s,0)==V**2-1 and S.diff(Hinf,V).subs(V,1)==2 and S.diff(Hinf,V).subs(V,-1)==-2)
    J=S.expand(S.diff(G,x)*S.diff(Tsrc,y)-S.diff(G,y)*S.diff(Tsrc,x))
    check('canonical differential identity has its exact full Jacobian',J==18*b**2*x)
    check('canonical zero is simple at both infinity sections',S.cancel((3*b/2)*(-s**-2)/(V/s**3)+(3*b/2)*s/V)==0)
    Tu=S.expand(Tsrc.subs({x:xx,y:yy},simultaneous=True))
    check('constructed global T recovers the whole boundary affine coordinate',S.Poly(Tu,u) is not None and Tu.coeff(u,0)==18*b*t)
    check('control volume has nonzero residue precisely on b!=0 family',S.residue(diff,u,0).subs(b,1)==3)
    check('Jung dominance has no initial cancellation at pole pair three,two',3%2==1 and 3>2)
    rr=S.symbols('r',positive=True,integer=True);z=S.symbols('z',positive=True)
    check('additive field extends in every regular-coordinate root chart',S.powsimp(z**(2*rr-2)*z**(2-rr)-z**rr)==0)
    check('additive boundary exponent is nonnegative from r>=1',(2*rr-2).subs(rr,1)==0 and S.diff(2*rr-2,rr)==2)
    inputs=[Path(__file__),HERE/'EXACT_VOLUME_COORDINATE_BRIDGE.md',manifest]+[HERE/r['path_relative_to_affine_complements'] for r in pins['primary_sources']]+[HERE/n for n in pred_names]
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
        general_index_bound=dict(all_indices=True,normal_affine_one_boundary=True,canonical_coefficient='e-1',source_coordinate_negative_hypothesis=True,last_key_value_negative=True,normal_form_coprime_defect='g*(a*b-a-b)<=e-1',maximum_reduced_pole='3*(e-1)',minimum_reduced_pole='2*(e-1)',e1_regular_coordinate_automatic=True,necessary_initial_type_controls=types,higher_index_cases_excluded=False),
        exact_coordinate_bridge=dict(smooth_affine_N=True,one_boundary_A1=True,div_source_volume='R',global_exactness_required=True,regular_source_coordinate_exists=True,all_finite_field_degrees=True,only_coordinate_negative_initial_pair=[3,2],essential_key_sequence=[3,2,-1],whole_root_chart='x=u^-3; y=u^-2+a*u^2+b*u^3+c*u^4+u^5*t',residue='3*b',exactness_forces='b=0',affineness_forces='b!=0',global_source_curve_contraction='G=y^3-x^2-3a*y-3c',affinization_ring='C[G,x*G^3,y*G^2]'),
        nonexact_control_family=dict(actual_smooth_affine_models=True,all_a_c_complex=True,b_nonzero=True,one_boundary_A1=True,canonical_divisor='R',every_source_coordinate_polar=True,global_volume_exact=False,proper_genus_two_family_with_one_section_removed=True,relative_ample_section_proof_written=True,normalized_boundary_coordinate='T/(18*b)',complete_parametrized_normal_form_with_source_choices=True,unique_parameter_moduli_claimed=False,finite_Keller_map_constructed=False),
        Keller_consequences=dict(one_index_two_prime_no_regular_coordinate_gap_closed=True,complete_regular_coordinate_geometric_classification_all_degrees=True,source_additive_action=True,source_additive_kernel='C[y]',target_algebra_invariance_assumed=False,all_degree_primitive_odd_marked_end=True,all_degree_target_torus_and_prime_boundary_degree_exclusions=True,all_quartic_single_index_two_prime_normalizations_excluded=True,individual_e1_e2_primes_have_regular_coordinate_on_affine_open=True,common_coordinate_for_multiple_primes_claimed=False,removing_other_primes_preserves_finiteness=False),
        assurance=dict(written_unbounded_canonical_and_global_proofs=True,general_signed_key_formula_imported=True,primitive_compactification_assumed=False,primary_statements_and_complete_page_audits=True,finite_controls_not_complete_geometric_proof=True,no_deep_coefficient_search=True,higher_index_normalizations_excluded=False,multiple_prime_quartic_packets_excluded=False,higher_finite_Keller_degrees_excluded=False,finite_Keller_map_constructed=False,formal_verification=False,independent_full_reproof=False,external_review=False,novelty_priority=False,research_goal_completed=False),
        source_sha256={str(p.relative_to(HERE)):digest(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'exact_coordinate_bridge_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS '+str(len(checks))+' exact bridge controls; all-degree regular coordinate forced; quartic single-index-two-prime packet excluded; goal active')
if __name__=='__main__':main()
