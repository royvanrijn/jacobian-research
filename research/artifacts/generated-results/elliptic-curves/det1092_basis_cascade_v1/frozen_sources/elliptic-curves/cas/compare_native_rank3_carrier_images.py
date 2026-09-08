#!/usr/bin/env python3
"""Post-construction catalogue and prior-equation comparison for twelve images."""
import argparse,json
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
INPUT=ART/'native_rank3_carrier_images_v1.json';PROOF=ART/'native_rank3_carrier_images_replay_v1.json';OLD=ART/'full11952_late64_r17_results_v1.json';CAT=ROOT/'artifacts/local/elliptic-curves/retention24-current-catalogue-v1/database.json';OUT=ART/'native_rank3_carrier_images_comparison_v1.json'
def expected():
    proof=cert.read(PROOF)
    if proof['status']!='PASS' or proof['distinct_j_invariants']!=12 or any(cert.hashed(ROOT/n)!=h for n,h in proof['sources'].items()):raise ArithmeticError('all twelve exact image proofs required before comparison')
    data=cert.read(INPUT);old=cert.read(OLD);catalogue=cert.read(CAT)['curves'];previous=old['previous_equations']+[{'address':OLD.name+':'+r['id'],'curve':r['curve']} for r in old['curves']]
    if len(previous)!=917 or len(catalogue)!=593 or cert.hashed(CAT)!=old['catalogue']['raw_sha256']:raise ArithmeticError('pinned comparison population differs')
    rows=[]
    for r in data['rows']:
        model=tuple(map(cert.F,r['curve']));rows.append({'word':r['word'],'compact_parameter':r['compact_parameter'],'catalogue_matches':[q['id'] for q in catalogue if cert.isomorphic(model,q['ainvs'])],'previous_matches':[q['address'] for q in previous if cert.isomorphic(model,q['curve'])]})
    paths=[Path(__file__).resolve(),Path(cert.__file__),INPUT,PROOF,OLD,CAT]
    return {'schema':'elliptic-curves.native-rank3-carrier-image-comparison.v1','status':'PASS','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'rows':rows,'catalogue_equations':593,'previous_equations':917,'unmatched_images':sum(not r['catalogue_matches'] and not r['previous_matches'] for r in rows),'claim_boundary':'Exact rational-isomorphism comparisons after all twelve fixed constructions and independent arithmetic proofs completed. Absence from this pinned catalogue and917 prior equations is relative novelty only. No specialized independence, high-rank, record or universal novelty claim. These comparisons do not enter the separate prospective MW16 selection or point execution.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('image comparison replay differs')
    else:
        if OUT.exists():raise FileExistsError('preserve post-construction comparison')
        checkpoint(OUT,d)
    print('NATIVE CARRIER IMAGE COMPARISON',d['unmatched_images'],'UNMATCHED OF12',flush=True)
