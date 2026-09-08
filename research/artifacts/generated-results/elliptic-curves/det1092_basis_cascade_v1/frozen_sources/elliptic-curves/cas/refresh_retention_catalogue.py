#!/usr/bin/env python3
"""Exact equation-only comparison against the independently frozen593-row update."""
import argparse
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/retention24-current-catalogue-v1';DB=D/'database.json';INDEX=ART/'new_high_rank_curve_index_v7.json';BATCH=ART/'retention24_r17_results_v1.json';OUT=ART/'retention_refreshed_catalogue_comparison_v1.json'
def result():
    intake=cert.read(D/'intake.json');db=cert.read(DB);index=cert.read(INDEX);batch=cert.read(BATCH)
    if cert.hashed(DB)!=intake['sha256'] or len(db['curves'])!=593:raise ArithmeticError('pinned update differs')
    projection=[{k:r[k] for k in ('id','ainvs','rank_lower_bound','conductor')} for r in db['curves']]
    inventory=[{'id':r['id'],'matches':[q['id'] for q in projection if cert.isomorphic(r['curve'],q['ainvs'])]} for r in index['curves']]
    candidates=[{'id':r['id'],'matches':[q['id'] for q in projection if cert.isomorphic(r['curve'],q['ainvs'])]} for r in batch['curves']]
    return {'schema':'elliptic-curves.retention-refreshed-catalogue.v1','status':'PASS','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),Path(cert.__file__).resolve(),ROOT/'elliptic-curves/cas/elliptic_candidate_record.py',DB,D/'intake.json',INDEX,BATCH)},'catalogue_intake':intake,'equations':projection,'inventory_comparisons':inventory,'batch_comparisons':candidates,'maximum_reported_catalogue_rank':max(r['rank_lower_bound'] for r in projection),'claim_boundary':'Exact rational-isomorphism comparisons with593 pinned equations after terminal fixed selection/search/replays. Public rank and conductor metadata are reported, not independently recertified here. Catalogue absence is not universal novelty; existing frozen586-row proofs remain intact.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();r=result()
    if a.check:
        if cert.read(OUT)!=r:raise ArithmeticError('catalogue comparison differs')
    else:
        if OUT.exists():raise FileExistsError('preserve refreshed comparison')
        checkpoint(OUT,r)
    print('EXACT593 CATALOGUE COMPARISON',len(r['inventory_comparisons']),'inventory curves; matches',[q for q in r['inventory_comparisons'] if q['matches']],flush=True)
