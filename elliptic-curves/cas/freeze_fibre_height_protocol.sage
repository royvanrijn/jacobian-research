#!/usr/bin/env sage-python
"""Extract only generic centres and exclusion j classes; no ranks/points enter selection.

Source files can contain outcomes. This one-time sanitizer emits only the
fixed twelve generic representatives, equation-derived j exclusions and pins.
"""
import json,gzip,hashlib,sys
from pathlib import Path
from sage.all import EllipticCurve,QQ,PolynomialRing
root=Path(__file__).resolve().parents[2];sys.path.insert(0,str(root/'elliptic-curves/cas'))
from importlib.machinery import SourceFileLoader
runner=SourceFileLoader('population',str(root/'elliptic-curves/cas/fibre_height_population.sage')).load_module()
bundle=root/'artifacts/generated-results/elliptic-curves/mw16_sensitivity_controls_v1.json.gz'
b=json.load(gzip.open(bundle,'rt'))
c=json.loads(b['files']['artifacts/local/elliptic-curves/mw16-sensitivity/generic-metrics.json']['text'])
r=next(r for r in c['results'] if r['parent_id']=='curve398-p16875')
s=next(s for s in r['settings'] if s['key']=='generic|metric:16|100000')
centres=[x['centre_construction']['representative'] for x in sorted(s['charts'],key=lambda c:c['mask'])]
paths=['elliptic-curves/data/icarm_mw16_parent_ladder_blind_inputs_v1.json','elliptic-curves/data/elkies_2026_r17_j_recognition_targets.json','artifacts/generated-results/elliptic-curves/icarm_mw16_nagao_finalist_specializations_h300_v1.json','artifacts/generated-results/elliptic-curves/a1_mw16_target_free_parameter_candidates_h300_v1.json','artifacts/generated-results/elkies-k3-r17-extreme-anchored-mw18-specializations-h1000-v1.json',str(bundle.relative_to(root)),'artifacts/generated-results/elkies-k3-r17-norm12-orbit07ca9-direct-fibration-v1.json']
excluded=set()
def add(model):excluded.add(str(EllipticCurve(QQ,model).j_invariant()))
for p in json.load(open(root/paths[0]))['parents']:add(p['target_short_model'])
for p in json.load(open(root/paths[1]))['targets']:add(p['ainvs'])
for path in paths[2:5]:
 d=json.load(open(root/path))
 for p in d.get('candidates',d.get('results',[])):
  if 'raw_short_model' in p:add(p['raw_short_model'])
cert=runner.read(runner.COVERS); R=PolynomialRing(QQ,'t')
for ch in cert['charts']:
 model=json.load(open(root/ch['direct_model']))['weierstrass_model'];A=R(model['A_coefficients_low_to_high']);B=R(model['B_coefficients_low_to_high'])
 for f in ch['fibres']:add([A(QQ(f['native_parameter'])),B(QQ(f['native_parameter']))])
h=cert['historical_rank28_anchor'];m=json.load(open(root/h['direct_model']));m=m.get('weierstrass_model',m)
add([R(m[key])(QQ(h['native_parameter'])) for key in ['A_coefficients_low_to_high','B_coefficients_low_to_high']])
paths.extend([ch['direct_model'] for ch in cert['charts']]);paths.append(h['direct_model'])
p=dict(nagao_keep_per_bucket=[16,8,4],nagao_bucket_width=1,schema='elliptic-curves.fibre-height-protocol.v1',mw16_presentation='a1-presentation-01',mw18_cover='07ca9-orbit-08c1e',parameter_height=6,arm_size=8,minimal_model_seconds=3,sample_height=1000,sample_seconds=2,sample_repeats=3,sample_section_indices=[0,1],search_height=100000,search_seconds=20,independence_prime_bound=1000,
 centres=dict(mw16=centres,mw18=[[int(i==j) for i in range(18)] for j in range(18)]),
 prime_blocks=dict(mw16=[list(x) for x in runner.source('elliptic-curves/cas/build_a1_mw16_target_free_parameter_candidates.sage').PRIME_BLOCKS],mw18=[list(x) for x in runner.source('elkies-k3/scripts/search_h92_q12o5867_rootless_nagao.py').DEFAULT_PRIME_BLOCKS]),
 excluded_j_invariants=sorted(excluded),source_hashes={p:runner.digest(root/p) for p in paths},
 exclusion_policy='Conservatively exclude entire j classes of known calibration equations and all retained prior prospective specializations. This also excludes rational twists of those curves.',
 population_policy='Fixed union of the identity and equation-derived improved H<=6 boxes; finite native coordinate only; exact Q-isomorphism deduplication before selection.',
 height_cost_policy='Equal-weight sum of empirical lower-is-better ranks for log2 j height, weighted integral global minimal Weierstrass size, and median repeated chart cost; no point outcome inputs.',
 nagao_policy='Existing three-stage capped Pareto policy and final ordering, applied separately to each coordinate population with width-one height buckets for H<=6, then merged and truncated to the common arm size.',
 centre_policy=dict(mw16='All twelve generic maximum-depth representatives for the fixed Gram, extracted without point outcomes from the certified calibration.',mw18='All eighteen singleton section centres; bounded pilot without a transferred sensitivity claim.'),
 interpretation='No gain in two finite arms is inconclusive about equivalence or usefulness. Report height separation, actual completed work and gain per worker time separately.')
if '--check' in sys.argv:
 assert runner.read(runner.PROTOCOL)==p, 'frozen protocol differs'
else:
 if runner.POPULATION.exists(): raise SystemExit('a population is already frozen; create a new version to change its protocol')
 runner.write(runner.PROTOCOL,p)
print('excluded',len(excluded),'centres',len(centres),flush=True)
