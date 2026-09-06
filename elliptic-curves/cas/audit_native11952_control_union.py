#!/usr/bin/env python3
"""Exact point-cloud ranks for old blind generic43 and its union with metric49."""
from pathlib import Path
import certify_compact_r17_candidates as cert
from audit_recorded_point_mod2_rank_v2 import build,check
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves'
D=LOCAL/'native11952-control-union-v1'

def main():
    oldpath=ART/'half_lattice_search_ablation_rank29_holdout_blind_v1.json';metricpath=LOCAL/'native11952-metric49-tails-v1/combined-point-cloud-only.json';parentpath=LOCAL/'native11952-metric49-control-v1/candidate-00/result.json'
    old=next(r for r in cert.read(oldpath)['results'] if r['label']=='curve12-2024-rank29');metric=cert.read(metricpath);parent=cert.read(parentpath)
    arm=next(r for r in old['arms'] if r['id']=='generic-deepest43');points=[]
    for i in arm['candidate_point_indices']:
        p=old['candidate_points'][i]['point'];points.append({'x':str(cert.F(p['x'])/36),'y':str(cert.F(p['y'])/216)})
    for name,charts,state in [('old43',[{'search':{'finite_curve_points':points}}],parent['initial_state']),('union',metric['charts']+[{'search':{'finite_curve_points':points}}],metric['final_state'])]:
        path=D/(name+'-point-cloud-only.json');output=ART/f'native11952_{name}_control_mod2_v1.json'
        if path.exists():raise FileExistsError('preserve diagnostic input')
        checkpoint(path,{'status':'POINT_CLOUD_ONLY','curve':metric['curve'],'generic_points':metric['generic_points'],'family':'native11952-control','parameter':'generic-only-fixture','final_state':state,'charts':charts,'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),oldpath,metricpath,parentpath)},'claim_boundary':'Old blind generic43 discoveries, transported exactly, and optional union with completed metric49 cloud. Not a chronological search history or a new curve.'})
        build(path,output,997,cert.hashed(path));check(output)
if __name__=='__main__':main()
