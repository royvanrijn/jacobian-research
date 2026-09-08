#!/usr/bin/env python3
"""Bind the fixed21 endpoint point campaign to its exact lower-bound evidence."""
import argparse,json
from collections import Counter
from pathlib import Path
import endpoint_specialized_parity_trial as batch
import complete_endpoint_specialized_trial as completion_driver
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ROOT=batch.ROOT;ART=batch.ART;OUT=ART/'endpoint_specialized_trial_v1.json'

def expected():
    p=batch.protocol();completion_driver.protocol();paths={Path(completion_driver.__file__).resolve(),batch.BATCH/'completion-v1/protocol.json',Path(__file__).resolve(),batch.BATCH/'protocol.json',batch.GATE,batch.INDEX}
    def read(path):paths.add(path);return cert.read(path)
    def terminal(path):
        s=read(path)
        if s['outcome']!='completed' or s['returncode']!=0:raise ArithmeticError('required completed supervision differs')
    completion=read(batch.BATCH/'completion-v1/ledger.json');terminal(batch.BATCH/'geometry.supervisor.json')
    if completion['status']!='PASS':raise ArithmeticError('complete endpoint roster and selected replays required')
    ledger={'status':'PASS','rows':[]}
    for i in range(21):
        row=read(batch.BATCH/'cloud-verification-v2'/f'{i:02}.json')
        if row['status']!='PASS' or len(row['rows'])!=1:raise ArithmeticError('complete21 individual cloud proofs required')
        ledger['rows'].extend(row['rows'])
    if ledger['status']!='PASS' or [r['id'] for r in ledger['rows']]!=[r['id'] for r in p['rows']]:raise ArithmeticError('all21 cloud rows required')
    rows=[]
    for index,(r,v) in enumerate(zip(p['rows'],ledger['rows'])):
        batch.configure(index);d=batch.D;raw=read(d/'result.json');maps=read(d/'maps.json');cloud=read(ROOT/v['mod2_certificate']);odd=read(ROOT/v['modl_certificate'])
        for label in ('maps','worker','mod2-build','mod2-check','modl-build','modl-check'):terminal(d/(label+'.supervisor.json'))
        replay=completion['rows'][index]
        if replay['id']!=r['id'] or replay['status']!='PASS':raise ArithmeticError('selected exact replay roster differs')
        for label in ('worker','replay'):
            supervision=ROOT/replay[label+'_supervision']
            if cert.hashed(supervision)!=replay[label+'_sha256']:raise ArithmeticError('selected supervision hash differs')
            terminal(supervision)
        if raw['initial_dimension']!=r['initial_rank'] or len(maps['sample'])!=256 or len(maps['rows'])!=12:raise ArithmeticError('endpoint seed/sample/map roster differs')
        if cert.hashed(d/'result.json')!=v['input_sha256'] or cert.hashed(ROOT/v['mod2_certificate'])!=v['mod2_sha256'] or cert.hashed(ROOT/v['modl_certificate'])!=v['modl_sha256']:raise ArithmeticError('point proof bindings differ')
        if cloud['rank_lower_bound']!=v['rank_lower_bound'] or {str(a['modulus']):a['finite_column_rank'] for a in odd['audits']}!=v['odd_modulus_lower_bounds']:raise ArithmeticError('finite point proofs disagree with ledger')
        attempted=len(raw['charts']);completed=sum(c['search']['status']=='bounded_search_complete' for c in raw['charts'])
        if attempted!=v['attempted_charts'] or completed!=v['completed_boxes'] or any(c['search']['height_bound']!=125000 or c['search']['timeout_seconds']!=10 for c in raw['charts']):raise ArithmeticError('recorded point budget differs')
        if raw['status']=='COMPLETE_DECLARED_ADAPTIVE_ATTEMPT':
            if attempted!=12:raise ArithmeticError('short complete endpoint trial')
        elif raw['status']=='TARGET_REACHED_PENDING_INDEPENDENT_REPLAY':
            if raw['rank_lower_bound']<22:raise ArithmeticError('unsupported target stop')
        else:raise ArithmeticError('nonterminal endpoint attempt')
        rows.append({**r,'rank_lower_bound':cloud['rank_lower_bound'],'odd_modulus_lower_bounds':v['odd_modulus_lower_bounds'],'retained_points':len(cloud['points']),'attempted_charts':attempted,'completed_boxes':completed,'status':raw['status'],'curve':raw['curve'],'point_certificates':{k:v[k] for k in ('mod2_certificate','mod2_sha256','modl_certificate','modl_sha256')}})
    return {'schema':'elliptic-curves.endpoint-specialized-trial.v1','status':'PASS','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in sorted(paths)},'rows':rows,'rank_lower_bound_counts':dict(sorted(Counter(r['rank_lower_bound'] for r in rows).items())),'declared_point_boxes':252,'attempted_point_boxes':sum(r['attempted_charts'] for r in rows),'completed_point_boxes':sum(r['completed_boxes'] for r in rows),'catalogue_snapshot_curves':593,'previous_address_equations':528,'claim_boundary':'All21 previously omitted nonsingular endpoints were chosen by exact equations before point searching. Each uses only its own independently certified11to17-point subset to define numerical specialized geometry. All256-mask samples and12-chart rosters were frozen before any point attempt. Trial clouds contain the certified seed plus every returned finite point; they are not an assertion that every generic section or possible rational point was returned. Endpoint equation nonisomorphism is certified separately by the prior exact audit. Finite subgroup lower bounds do not prove exact rank, rank upper bounds, saturation, absence of points, generic-section dependence, sampled-CVP optimality, universal sensitivity or universal novelty. High-rank inventory promotion requires a separate certificate and inventory replay.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('endpoint point aggregate differs')
    else:
        if OUT.exists():raise FileExistsError('preserve endpoint point aggregate')
        checkpoint(OUT,d)
    print('ENDPOINT TRIAL PASS',d['rank_lower_bound_counts'],d['completed_point_boxes'],'completed boxes',flush=True)
