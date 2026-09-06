#!/usr/bin/env python3
"""Bind the completed fixed own27-point follow-up without another search."""
import argparse
import json
from pathlib import Path
import certify_compact_r17_candidates as cert
import full11952_specialized_followup as batch
from research_runtime.store import checkpoint

ROOT=batch.ROOT
ART=batch.ART
CONTROL=batch.LOCAL/'full11952-specialized-followup-controller-v1'
OUT=ART/'full11952_specialized_followup_v1.json'


def expected():
    p=batch.protocol();paths={Path(__file__).resolve(),batch.BATCH/'protocol.json'}
    def read(path):paths.add(path);return cert.read(path)
    controller=read(CONTROL/'protocol.json');ledger=read(CONTROL/'ledger.json')
    if any(cert.hashed(ROOT/n)!=h for n,h in controller['sources'].items()):raise ArithmeticError('frozen controller sources differ')
    if ledger['status']!='PASS' or ledger['stronger_odd_prime_bounds'] or [r['name'] for r in ledger['rows']]!=['freeze','trial','clouds','geometry']:raise ArithmeticError('complete four-stage follow-up required')
    for row in ledger['rows']:
        s=row['supervision'];log=CONTROL/(row['name']+'.log');paths.add(log)
        if row['status']!='PASS' or s['outcome']!='completed' or s['returncode']!=0 or cert.hashed(log)!=s['log_sha256']:raise ArithmeticError('follow-up transcript differs')
    execution=read(batch.BATCH/'ledger.json');proof=read(batch.BATCH/'verification-ledger.json')
    roster=[r['id'] for r in p['rows']]
    if len(roster)!=1 or execution['status']!='PASS' or proof['status']!='PASS' or [r['id'] for r in execution['rows']]!=roster or [r['id'] for r in proof['rows']]!=roster:raise ArithmeticError('the one follow-up attempts and clouds required')
    if [r['id'] for r in execution['maps']]!=roster or any(r['status']!='PASS' for r in execution['maps']):raise ArithmeticError('the one prior geometry stages required')
    rows=[]
    for r,worker,check in zip(p['rows'],execution['rows'],proof['rows']):
        folder=batch.BATCH/r['id'];raw=read(folder/'result.json');maps=read(folder/'maps.json');cloud=read(ROOT/check['mod2_certificate']);odd=read(ROOT/check['modl_certificate'])
        if worker['status']!='PASS' or check['status']!='PASS' or cert.hashed(folder/'result.json')!=worker['result_sha256'] or worker['result_sha256']!=check['input_sha256'] or cert.hashed(ROOT/check['mod2_certificate'])!=check['mod2_sha256'] or cert.hashed(ROOT/check['modl_certificate'])!=check['modl_sha256']:raise ArithmeticError('completed point/cloud hashes differ')
        if cert.hashed(folder/'maps.json')!=raw['maps_sha256'] or len(maps['sample'])!=2048 or len(maps['rows'])!=49 or raw['initial_dimension']!=r['initial_rank']:raise ArithmeticError('fixed seed, sample or map binding differs')
        if raw['status']!='COMPLETE_DECLARED_ADAPTIVE_ATTEMPT' or len(raw['charts'])!=49 or any(c['search']['status']!='bounded_search_complete' or c['search']['height_bound']!=125000 or c['search']['timeout_seconds']!=10 for c in raw['charts']):raise ArithmeticError('complete fixed49 point exposure required')
        bounds={str(a['modulus']):a['finite_column_rank'] for a in odd['audits']}
        if cloud['rank_lower_bound']!=27 or bounds!={'3':27,'5':27} or check['rank_lower_bound']!=27 or bounds!=check['odd_modulus_lower_bounds']:raise ArithmeticError('observed matching27 bounds differ')
        rows.append({'id':r['id'],'family':r['family'],'parameter':r['parameter'],'initial_rank':r['initial_rank'],'rank_lower_bound':27,'odd_modulus_lower_bounds':bounds,'completed_boxes':49,'retained_points':len(cloud['points'])})
    return {'schema':'elliptic-curves.full11952-specialized-followup-result.v1','status':'PASS','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in sorted(paths)},'rows':rows,'declared_point_boxes':49,'completed_point_boxes':49,'total_retained_point_witnesses':sum(r['retained_points'] for r in rows),'claim_boundary':'One complete49-box specialized-parity follow-up of the one new full11952 curves with27 certified points. All2048-mask samples and49-map rosters precede points; exact histories, rational geometry and full-cloud proofs modulo2,3,5 pass. Every retained cloud still certifies27, giving no rank gain or new inventory entry. These are lower bounds, not exact ranks, upper bounds, saturation, covering optimality, point absence or a universal sensitivity claim.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');args=parser.parse_args();d=expected()
    if args.check:
        if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('completed follow-up summary differs')
    else:
        if OUT.exists():raise FileExistsError('preserve completed follow-up proof')
        checkpoint(OUT,d)
    print('ONE27 FOLLOWUP: ALL49 BOXES PASS;',d['total_retained_point_witnesses'],'WITNESSES; NO BOUND GAIN',flush=True)
