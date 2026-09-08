#!/usr/bin/env python3
"""Fixture preparation only. Never launches a reconstruction or evaluates a result."""
import hashlib,json,subprocess,datetime
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4]
CAS=Path(__file__).resolve().parent
PKG=ROOT/'research/artifacts/generated-results/elliptic-curves/det1092_blind_mw16_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,d):
    with p.open('x') as f:json.dump(d,f,indent=2,sort_keys=True);f.write('\n')
    p.chmod(0o444)
commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
parent_path=ROOT/'research/artifacts/generated-results/elliptic-curves/curve302_recovered_mw17_parent_v1.json'
parent=json.loads(parent_path.read_text())
assert not (PKG/'protocol.json').exists()
policy={
 'objective':'Recover a rational section over Q(t) independent of this retained rank16 core; known parent rank17 is a positive control, not a new rank claim.',
 'centre_selection':'LLL_gram columns verified unimodular; 512 SHA256 parity masks in reduced coordinates, first two digest bytes of det1092-blind-mw16-v1:j, j=0..511, big endian. Deterministic greedy moves +/-2ei and +/-2(ei+/-ej), strictly reducing norm, at most128 moves. Normalize original word sign. Deduplicate original parities; select up to32 with heights8..14 ordered by decreasing height, original L1 word size, original lexicographic word. No oracle, full-lattice coset data, or coefficient-fitted constants.',
 'RR':'For each selected centre T, h=height(T), use n=ceil(3h/4)+1 then n+1. Solve f0+f1*x(T)+f2*y(T)=0 in degrees(n,n-4,n-6); take each exact RREF kernel basis row in order, no fitted linear combinations. Remove T from the cubic-line intersection. Test discriminant square in Q(t) and return every distinct rational root of a split member.',
 'candidate_gate':'Exact equation membership and the K3 Shioda height Gram on the original parent; rank17 and positive Schur complement required. No finite-field or specialization independence proxy.',
 'stop':'First RR member giving exact rank17, after retaining both rational roots. Otherwise exhaust selected centres and both degree bounds. Preserve every failed stage; no refill/replacement/rerun.',
 'limits':{'parity_samples':512,'greedy_steps_per_sample':128,'selected_centres':32,'RR_bounds_per_centre':2,'CPU_seconds':180,'CPU_hard_seconds':185,'wall_seconds':240,'address_space_GiB':8,'per_file_MiB':128,'workers_concurrent':1},
 'allowed':['same original parent a-invariants','exact16 retained rational sections','principal16x16 height Gram independently recomputed from those sections','known full generic rank17','general RR, elliptic group law and K3 height formula','deterministic algorithm and CAS runtime'],
 'forbidden':['omitted section coordinates or other-basis words','omitted height pairings','full17 parent/basis/lattice artifacts during reconstruction','other arms or outputs','exceptional302 points, public rank31 set, exceptional14 complement, V1/V2/V3 data','historical winning centres/orbits/nets/covers or fitted selectors/coordinate transformations','specialization points or tests during reconstruction','network and repository access outside the explicit own-fixture and worker-file allowlist'],
 'classification_precedence':'All arms terminal and replay passing:17 successes BROAD_16_TO_17_RECOVERY;1..16 CORE_DEPENDENT_16_TO_17_RECOVERY;0 with complete execution BOUNDED_NO_RECOVERY;0 and any failure UNRESOLVED_DUE_TO_EXPERIMENT_FAILURE. Alternative quotient representatives are reported independently; ALTERNATIVE_QUOTIENT_RECOVERY reserved for a result whose only demonstrated successes are nonprimitive quotient multiples.',
 'scope':'Only the declared fresh RR split-component mechanism is tested. A miss cannot attribute failure to every recent geometric construction or imply nonexistence.'}
(PKG/'fixtures').mkdir(exist_ok=True);(PKG/'frozen_sources').mkdir()
sources=[parent_path,ROOT/'research/MATH_STATUS.json',ROOT/'research/elliptic-curves/README.md',ROOT/'AGENTS.md',ROOT/'research/elliptic-curves/AGENTS.md',ROOT/'research/elliptic-curves/cas/load_curve302_recovered_parent.sage']
for pattern in ['CURVE302_RECOVERED_MW17_PARENT_2026-09-07.md','CURVE302_PARENT_SEARCH_AND_GEOMETRY_2026-09-07.md','DET1092_RANK18_BASE_CHANGE_AND_INITIAL_UNLOCK_2026-09-08.md','DET1092_TWO_COVER_BRANCH_GATE_2026-09-08.md','DET1092_FIRST_UNLOCK_RR_NET_2026-09-08.md','DET1092_SINGLE_SEED_COVER_2026-09-08.md','DET1092_FIRST_UNLOCK_JACOBIAN_CLASS_2026-09-08.md','DET1092_RR_GENERIC_POINT_SPECIFICITY_2026-09-08.md','DET1092_SEED_KUMMER_COVER_2026-09-08.md','DET1092_RR_NET_REDUCIBLE_LOCUS_2026-09-08.md','RANK_JUMP_REASSESSMENT_2026-09-05.md']:
 sources.append(ROOT/'research/elliptic-curves/notes'/pattern)
for name in ['curve302_recovered_mw17_parent_proof_v1.json','curve302_parent_geometric_picard19_v1.json']:
 sources.append(parent_path.parent/name)
sources += [CAS/'worker.py',CAS/'isolate.py',Path(__file__)]
bindings={str(p.relative_to(ROOT)):sha(p) for p in sources}
for p in sources:
 dest=PKG/'frozen_sources'/p.relative_to(ROOT)
 dest.parent.mkdir(parents=True,exist_ok=True)
 dest.write_bytes(p.read_bytes());dest.chmod(0o444)
roster=[]
for omitted in range(17):
 keep=[j for j in range(17) if j!=omitted];aid='arm-%02d'%(omitted+1)
 fixture={'starting_commit':commit,'arm_id':aid,'a_invariants':parent['a_invariants'],
          'sections':[parent['basis_weierstrass_coordinates'][j] for j in keep],
          'gram':[[parent['generic_height_gram'][j][k] for k in keep] for j in keep],'policy':policy}
 fp=PKG/'fixtures'/(aid+'.json');write(fp,fixture)
 roster.append({'arm_id':aid,'omitted_basis_index_one_based':omitted+1,'retained_basis_indices_one_based':[j+1 for j in keep],'fixture':str(fp.relative_to(PKG)),'fixture_sha256':sha(fp)})
write(PKG/'roster.json',{'starting_commit':commit,'arms':roster,'immutable':True})
write(PKG/'protocol.json',{'starting_commit':commit,'frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'roster_sha256':sha(PKG/'roster.json'),'source_bindings':bindings,'policy':policy,'stages':['fixture_preparation','isolated_reconstruction_all17','freeze_primary_outputs','post_run_full_basis_evaluation','independent_replay','secondary_specialization_diagnostic'],'blinding':'Algorithm/process input blindness; designer has read canonical background. No claim of human ignorance of the known parent. Linux Landlock ABI7 allowlist with separate network/PID namespaces. Only own fixture and worker source from repository are readable.'})
status=subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True)
diff=subprocess.check_output(['git','diff','--','research/MATH_STATUS.json','research/elliptic-curves/README.md','research/elliptic-curves/notes'],cwd=ROOT,text=True)
(PKG/'preflight/starting_worktree.txt').write_text('Starting commit: '+commit+'\n'+status+'\n'+diff)
print('FROZEN',commit,'17 arms',sha(PKG/'protocol.json'))
