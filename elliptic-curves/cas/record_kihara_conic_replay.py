#!/usr/bin/env python3
"""Record the conic proof, primitive subgroup and pinned control comparison."""
import argparse,json,hashlib
from pathlib import Path
import certify_compact_r17_candidates as cert
import audit_inventory200_current_catalogue as catalogue
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';D=LOCAL/'kihara-conic-standalone-v1';OUT=ART/'kihara_conic_replay_v1.json'
def compute():
 stages={};paths=[Path(__file__).resolve()]
 for name,folder,log in [('construction','kihara-split-infinity-v1','worker.log'),('control','kihara-conic-control-v1','worker.log'),('independent','kihara-conic-standalone-v1','replay.log')]:
  d=LOCAL/folder;r=cert.read(d/'supervisor.json');assert r['outcome']=='completed' and r['returncode']==0 and r['failure_reason'] is None and cert.hashed(d/log)==r['log_sha256'];stages[name]=r;paths.extend([d/'protocol.json',d/'supervisor.json',d/log])
 transcript=(D/'replay.log').read_text();assert transcript.splitlines()==['PASS explicit conic, all infinity images, primitive rank13 subgroup and anti-height10','PASS balanced v774161 control,998-bit model and independent rational rank13']
 for p in [ROOT/'elliptic-curves/cas/verify_kihara_conic.sage',ART/'kihara_conic_replay_bundle_v1.json']:
  assert p.read_bytes()==(D/p.name).read_bytes();paths.append(p)
 control_path=ART/'kihara_conic_control_v1.json';control=cert.read(control_path);paths.append(control_path)
 database=catalogue.D/'database.json';metadata=catalogue.D/'metadata.json';index=ART/'new_high_rank_curve_index_v22.json';db=cert.read(database);inv=cert.read(index)
 assert cert.hashed(database)==cert.read(metadata)['sha256'] and len(db['curves'])==620 and len(inv['curves'])==201;paths.extend([database,metadata,index,Path(catalogue.__file__)])
 j=catalogue.j(control['curve']);public=[r['id'] for r in db['curves'] if catalogue.j(r['ainvs'])==j and cert.isomorphic(control['curve'],r['ainvs'])];local=[r['id'] for r in inv['curves'] if catalogue.j(r['curve'])==j and cert.isomorphic(control['curve'],r['curve'])]
 return {'schema':'elliptic-curves.kihara-conic-replay.v1','status':'PASS','generic_subgroup_rank':13,'generic_subgroup_saturated_in_its_rational_span':True,'ambient_generic_rank':'UNKNOWN','anti_invariant_section_height':'10','geometric_nontorsion_height_lower_bound':'6','base_chi':4,
 'control_rank_lower_bound':13,'control_model_bits':998,'control_catalogue_matches':public,'control_inventory_matches':local,'point_search_boxes':0,'stages':stages,'total_supervised_seconds':sum(r['wall_seconds'] for r in stages.values()),'transcript':transcript,
 'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},
 'scope':'One explicit rational conic base change over the first new Kihara K3 parent. The13-section subgroup is primitive by the deck action, anti-height10 and universal height lower6; full ambient rank remains unknown. One balanced control has13 independent points. Pinned catalogue/inventory comparisons are post-construction and do not prove literature novelty. No new parent count, near-record curve, searched rank gain, or superiority to the existing rank14 path is claimed. No point search or following population campaign.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args();r=compute()
 if args.check:assert r==cert.read(OUT)
 else:
  with OUT.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
 print('PASS primitive generic13 conic; one rank13 control;',r['total_supervised_seconds'],'seconds; matches',r['control_catalogue_matches'],r['control_inventory_matches'])
