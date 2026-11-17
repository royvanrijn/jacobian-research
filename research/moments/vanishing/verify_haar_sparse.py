"""Exact controls for the written logarithmic-jet proof; no theorem from samples."""
from pathlib import Path
import hashlib,json,itertools
import sympy as S
p=Path(__file__).resolve().parent;T=S.symbols('T')
supports=[[-9,4,7],[-5,-1,3,9],[-20,-3,1,7,9,25]];projectors=[]
for support in supports:
 polys=[]
 for i in support:
  L=S.prod((T-j)/S.Integer(i-j) for j in support if j!=i)
  assert S.degree(L,T)==len(support)-1
  assert [S.simplify(L.subs(T,j)) for j in support]==[int(j==i) for j in support]
  polys.append(str(S.expand(L)))
 projectors.append({'charges':support,'projectors':polys})
# Exact algebra of the simple-root strengthened inequality.
k,d,h,D,b,a,c=S.symbols('k d h D b a c')
assert S.expand((k-1)*(d-h)-(D*b+a*c-h))==S.expand((k-1)*d-(D*b+a*c+(k-2)*h))
# Sufficient ring cutoff follows using z+nminor<=nu-1.
z,n,nu,w,W,B=S.symbols('z n nu w W B')
assert S.expand((z+(1+W/w)*n+B/w).subs(z,nu-1-n))==S.expand(nu-1+(W/w)*n+B/w)
capacity_controls=[]
for divisor,minimum,bound in [(3,4,4),(6,6,6)]:
 for degree in [minimum,minimum+1,2*divisor,2*divisor+1,3*divisor,7*divisor+2]:
  mmax=degree//divisor
  cap=degree if mmax<=1 else degree-mmax
  assert cap>=bound
  capacity_controls.append({'orbit_half':divisor,'D':degree,'root_bound':mmax,'capacity_lower_bound':cap})
# Enumerate monomial COUNTS ONLY for a selected ring cutoff transcription.
# These are not phase or support searches and do not establish the universal proof.
count=0;nu0=3;w0=2;W0=7;B0=9
for z0,nm,n2,n7 in itertools.product(range(19),repeat=4):
 total=z0+nm+n2+n7
 if total>18:continue
 charge=-5*nm+2*n2+7*n7
 if abs(charge)>B0:continue
 paired=z0+min(nm,n2+n7)
 if w0*total>(W0+w0)*(nu0-1)+B0:assert paired>=nu0
 count+=1
# Componentwise orientation control CxC: both ideals nonzero, product zero.
eplus=(1,0);eminus=(0,1)
assert tuple(a*b for a,b in zip(eplus,eminus))==(0,0)
assert tuple(a*a for a in eplus)==eplus and tuple(a*a for a in eminus)==eminus
rec={'date':'2026-11-16','completed':True,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'sympy_version':S.__version__,
 'Vandermonde_projection_controls':projectors,'simple_root_improved_inequality_exact':True,
 'capacity_threshold_transcription_controls':capacity_controls,'ring_cutoff_balanced_count_controls':count,
 'nonreduced_or_componentwise_ring_orientation_retained':True,
 'logarithmic_derivation_bound_all_valuations_written_proof':True,
 'general_two_sign_capacity_proper_model_criterion_written_proof':True,
 'four_weight_2O_six_weight_2I_all_degree_classification_written_proof':True,
 'all_ring_J_nilpotence_iff_eventual_purity_written_proof':True,
 'all_original_marker_cutoff_written_proof':True,'all_H_mixed_nonpurity_written_proof':True,
 'minimum_failure_weights_after_every_auto_torus':{'2O':5,'2I':7},
 'new_pure_failure_produced':False,'full_MZ_cases_claimed_closed':False,
 'all_mixed_phases_claimed_proper_models':False,'primitive_purity_reflected_from_composite':False,
 'nilpotence_index_claimed_computed_without_certificate':False,'finitebase_volume_extension_claimed_full_unit_elimination':False,
 'actual_period_kernel_assumed_Hodge':False,'one_global_ring_sign_assumed':False,
 'weighted_to_unweighted_Haar_transfer_used':False,'universal_geometry_or_MZ_machine_proof':False,
 'independent_replay':False,'formal_verification':False,'external_review':False,
 'external_review_next_task':False,'new_degree_support_or_parameter_campaign':False,
 'full_exceptional_MZ_cases':'OPEN on BOTH quotients'}
(p/'haar_sparse_verification.json').write_text(json.dumps(rec,indent=2)+'\n')
print(json.dumps({'completed':True,'exact_controls':'PASS','balanced_count_controls':count,'checker_sha256':rec['checker_sha256']}))
