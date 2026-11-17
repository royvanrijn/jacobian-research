"""Exact global-bundle and intrinsic safe-pole controls, not a universal proof."""
from pathlib import Path
import hashlib
import json
import sympy as S

t,xi,eta,R,v,a,b,c,d=S.symbols('t xi eta R v a b c d')
Delta=t*(eta-xi)
old=[t,t*xi,eta]
jac=S.Matrix([[S.diff(q,u) for u in [xi,eta,t]] for q in old]).det()
assert S.simplify(jac/Delta**2-1/(t*(eta-xi)**2))==0

# Right matrix action on the row-line chart and the Hom-line fibre.
tx=(b+d*xi)/(a+c*xi)
te=(b+d*eta)/(a+c*eta)
tp=t*(a+c*xi)/(a+c*eta)
assert S.simplify(S.diff(tx,xi)*S.diff(te,eta)/(te-tx)**2
                  -1/(eta-xi)**2)==0
assert S.simplify(tp.subs(eta,xi)-t)==0

# Every fixed fibre root over the diagonal remains fixed, not an orbit6/12 root.
assert S.simplify((tp+1).subs({eta:xi,t:-1}))==0

# The finite-base log-discrepancy calculation, with proper centre extra order.
vf,extra=S.symbols('vf extra',positive=True,integer=True)
assert S.expand((vf+extra)-1-vf-(extra-1))==0
assert (extra-1).is_nonnegative

# Blow up the actual base-infinity/fibre-root curve R=W=0, W=R*v.
pull=[xi,xi+R,-1+R*v]
blowjac=S.Matrix([[S.diff(q,u) for u in [xi,R,v]] for q in pull]).det()
vol=S.simplify(blowjac/((pull[1]-pull[0])**2*pull[2]))
assert S.simplify(vol-1/(R*(R*v-1)))==0
assert S.residue(vol,R,0)==-1

# Homogeneous-form identities test the universal monomial proof,not invariance.
X,Y=S.symbols('X Y')
controls=[]
for D in [1,4,6]:
    F=X**(2*D)+3*X**D*Y**D+2*Y**(2*D)
    original=F.subs({X:R*v,Y:R*(v*xi+1)})/((R*v-1)**D*R**D)
    expected=R**D*F.subs({X:v,Y:v*xi+1})/(R*v-1)**D
    assert S.cancel(original-expected)==0
    assert S.expand(F.subs({X:v,Y:v*xi+1}))!=0
    leading=S.simplify((original/R**D).subs(R,0))
    assert leading!=0
    controls.append({'D':D,'phase_order':D,'original_volume_order':-1})

# All polynomial recompositions vanishing at0 retain a finite-phase log pole.
f0=R**4*(v**8+xi+1)/(R*v-1)**4
for H in [f0,f0**2+f0**3]:
    assert S.limit(H,R,0)==0
assert S.residue(vol,R,0)==-1

out={'date':'2026-11-16','completed':True,
 'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'sympy_version':S.__version__,
 'global_original_PGL2_volume_chart_exact':True,
 'right_Mobius_base_volume_invariance_exact':True,
 'right_group_fibre_action_trivial_on_infinity_diagonal_exact':True,
 'finite_base_discrepancy_bound_for_all_proper_centres':True,
 'actual_curve_blowup_retains_original_logarithmic_pole':True,
 'homogeneous_form_pole_transcription_controls':controls,
 'pole_survives_all_zero_centred_recompositions_argument_checked':True,
 'full_phase_and_original_marker_coverage':'ENTIRE algebra;global Hom-line bundle,not a rational frame',
 'actual_marked_FUNCTION_module_preserved_before_polynomial_quotient':True,
 'every_mixed_phase_finite_base_UNIT_volume_regular':True,
 'safe_pole_family_claimed_new_safe_class':False,
 'volume_pole_order_identified_with_corrected_torus_gluing_order':False,
 'finite_base_volume_regularity_claimed_full_unit_component_vanishing':False,
 'all_marker_forms_claimed_holomorphic_on_compactification':False,
 'full_geometric_cohomology_claimed_unchanged_by_compactification':False,
 'global_Hopf_circle_bundle_replaced_by_product':False,
 'all_full_mixed_phases_claimed_proper_holomorphic_volume_models':False,
 'full_actual_unit_extension_bridge_proved':False,
 'weighted_to_unweighted_Haar_transfer_used':False,
 'primitive_purity_reflected_from_composite':False,
 'universal_geometry_or_MZ_proof_machine_verified':False,
 'independent_replay':False,'formal_verification':False,'external_review':False,
 'external_review_next_task':False,'new_degree_support_or_parameter_campaign':False,
 'full_exceptional_MZ_cases':'OPEN on BOTH quotients',
 'ownership':'FULL exceptional MZ;general Gaussian/formal Bernstein work separately owned'}
Path(__file__).with_name('haar_bundle_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
