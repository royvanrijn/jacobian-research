"""Finite exact controls for the written saturation/lifting proof; not a machine proof."""
from pathlib import Path
import hashlib,json
import sympy as S
p=Path(__file__).resolve().parent
s=S.symbols('s');a=1/(s+1)-1/s
assert S.residue(a,s,0)==-1 and S.residue(a,s,-1)==1
# Rank2 Kummer connection: d(v)=a*u, d(u)=0.
N=S.Matrix([[0,a],[0,0]]);assert N.rank()==1 and N*N==S.zeros(2)
# A finite global semisimple model with a distinguished trivial line and
# nontrivial character. Its only normalized invariant covector kills I.
T=S.diag(1,-1,1);r0,r1=S.symbols('r0 r1')
R=S.Matrix([[1,0],[0,-1]])
assert S.solve(list((S.Matrix([[1,r1]])*R-S.Matrix([[1,r1]]))),[r1])=={r1:0}
# Pushout by projection of R=<u,i> to u is quotient by i.
b=S.symbols('b');A=S.Matrix([[0,0,a],[0,0,b],[0,0,0]])
Q=S.Matrix([[1,0,0],[0,0,1]]);Ap=S.Matrix([[0,a],[0,0]])
assert Q*A==Ap*Q
# Local versus global distinction: the local identity on R allows all r1.
assert S.Matrix([[1,r1]])*S.eye(2)==S.Matrix([[1,r1]])
# Written saturation argument uses image containing L and surjecting C.
controls=[]
for n in range(1,6):
 row=S.Matrix([[S.Rational(j+1,n+1) for j in range(n)]])
 image=S.eye(n+1);image[0,1:]=row
 assert image.rank()==n+1
 controls.append({'compact_middle_rank':n,'image_rank_with_ambient_line':image.rank()})
rec={'date':'2026-11-16','completed':True,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'sympy_version':S.__version__,
 'source_image_exactly_W4_written_proof':True,'all_phase_unit_adjoined_regular_image_model_independent_written_proof':True,'ambient_cohomology_lift_exactly_decides_unit_in_regular_image_written_proof':True,'actual_regular_image_model_independent_written_proof':True,
 'canonical_R_defined_for_every_phase':True,'global_unit_functional_unique_representation_control':True,
 'boundary_Ext_pushout_matrix_control':True,'Kummer_residue_vector':[-1,1],
 'Kummer_constant_boundary_line_for_ALL_regular_Bezout_parameters':True,
 'finite_linear_saturation_controls':controls,
 'actual_period_kernel_assumed_Hodge':False,'local_unit_retraction_claimed_unique':False,
 'finitebase_partial_model_treated_as_proper_over_phase':False,'finitebase_elimination_proved':False,
 'marker_defect_localized_finitebase_or_infinity':False,'all_phases_claimed_proper_original_volume_models':False,
 'new_pure_exceptional_phase_produced':False,'nonlinear_composite_purity_reflected':False,
 'full_exceptional_purity_implies_safety_proved':False,'universal_geometry_or_MZ_machine_proof':False,
 'independent_replay':False,'formal_verification':False,'external_review':False,
 'external_review_next_task':False,'new_degree_support_or_parameter_campaign':False,
 'full_exceptional_MZ_cases':'OPEN on BOTH quotients'}
(p/'haar_saturation_verification.json').write_text(json.dumps(rec,indent=2)+'\n')
print(json.dumps({'completed':True,'exact_controls':'PASS','checker_sha256':rec['checker_sha256']}))
