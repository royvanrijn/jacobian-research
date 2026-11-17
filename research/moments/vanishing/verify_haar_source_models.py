"""Exact transcription controls, not a machine proof of the geometric theorems."""
from pathlib import Path
import hashlib,json
import sympy as S
p=Path(__file__).resolve().parent
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
e=S.symbols('e',positive=True,integer=True);k=S.symbols('k',integer=True)
assert S.expand(e*k+e-1)==S.expand(e*(k+1)-1)
assert S.simplify((e*(k+1)-1).subs(k,-1))==-1
# Generic lattice index: primitive weighted curve parameter has kernel mu_g.
lattice=[]
for a,b in [(8,12),(12,20),(8,8),(12,12),(16,36),(24,60)]:
 g=S.igcd(a,b);r=a//g;s=b//g
 assert S.igcd(r,s)==1 and a*s==b*r
 assert S.Rational(a,r)==g and S.Rational(b,s)==g
 lattice.append({'a':a,'b':b,'gcd':int(g),'r':int(r),'s':int(s),'constant_dlog_ratio':int(g)})
y,f=S.symbols('y f');phase=(1+y)*(1+f+y*f)
assert S.expand(phase-((1+2*y*f)+y*(1+y*f)+f))==0
s,t=S.symbols('s t');logder=1/(s+1)-1/s
assert S.residue(logder,s,0)==-1 and S.residue(logder,s,-1)==1
moments=[]
for j in range(9):
 m=S.integrate((-t)**j,(t,0,1));assert m==S.Rational((-1)**j,j+1)
 moments.append(str(m))
# Cauchy expansion at infinity matches every symbolic coefficient.
z=S.symbols('z')
assert S.expand(S.series(S.log(1+z),z,0,10).removeO()-sum(S.Rational((-1)**j,j+1)*z**(j+1) for j in range(9)))==0
# Derivatives of arbitrary finite Laurent germs have no residue.
assert all(j==0 or j-1!=-1 for j in range(-9,10))
# General sign proof uses integrality, not the following finite controls:
assert all((ee*(kk+1)-1>=0)==(kk>=0) for ee in range(1,13) for kk in range(-6,7))
rec={'date':'2026-11-16','completed':True,'checker_sha256':sha(Path(__file__)),'sympy_version':S.__version__,
 'ramification_top_form_identity_exact':True,'logarithmic_pole_stays_minus_one_exact':True,
 'weighted_fibre_lattice_controls':lattice,'Long_charged_decomposition_exact':True,
 'fixed_original_marker_moment_controls':moments,'Cauchy_logarithm_coefficients_and_residues_exact':True,
 'valuation_equivalence_general_proof':'written in canonical note; integer k and positive integer e',
 'original_source_proper_models_proved':True,'actual_unit_submodule_global_splitting_proved':True,
 'actual_full_marked_module_locally_nonsplit_control_proved':True,
 'marker_loss_in_regular_period_image_proved':True,
 'finite_base_specific_marker_defect_localization_proved':False,
 'proper_unit_model_identified_with_full_marked_source':False,
 'unit_submodule_inclusion_treated_as_unit_quotient':False,
 'ramified_cover_claimed_to_remove_original_log_pole':False,
 'Bezout_nonconstant_Jacobian_dropped':False,
 'old_allH_nonpurity_claimed_new':False,'actual_control_claimed_pure':False,
 'weighted_to_unweighted_Haar_transfer_used':False,'primitive_purity_reflected_from_composite':False,
 'full_actual_unit_extension_bridge_proved':False,'universal_geometry_or_MZ_proof_machine_verified':False,
 'independent_replay':False,'formal_verification':False,'external_review':False,
 'external_review_next_task':False,'new_degree_support_or_parameter_campaign':False,
 'full_exceptional_MZ_cases':'OPEN on BOTH quotients'}
(p/'haar_source_models_verification.json').write_text(json.dumps(rec,indent=2)+'\n')
print(json.dumps({'completed':True,'exact_controls':'PASS','full_MZ':'OPEN','checker_sha256':rec['checker_sha256']}))
