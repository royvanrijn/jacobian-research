#!/usr/bin/env python3
"""Bind completed retained native19 exposures, proofs and measured cost."""
import argparse
from pathlib import Path
import retained_native19_trial_v3 as batch
from research_runtime.store import checkpoint
ROOT,ART=batch.ROOT,batch.ART
OUT=ART/'retained_native19_trial_v3_results_v1.json'

def expected():
    p=batch.protocol();paths={Path(__file__).resolve(),batch.BATCH/'protocol.json'}
    def read(q):paths.add(q);return batch.cert.read(q)
    run=read(batch.BATCH/'ledger.json');audit=read(batch.BATCH/'verification-ledger.json');geo=read(batch.BATCH/'geometry.supervisor.json')
    assert run['status']==audit['status']=='PASS' and geo['outcome']=='completed' and geo['returncode']==0
    assert len(run['maps'])==len(run['rows'])==len(audit['rows'])==len(p['rows'])==2
    assert batch.cert.hashed(batch.BATCH/'geometry.log')==geo['log_sha256'];paths.add(batch.BATCH/'geometry.log')
    rows=[];total=geo['wall_seconds']
    for metadata,execution,verified in zip(p['rows'],run['rows'],audit['rows']):
        ident=metadata['id'];assert execution['id']==verified['id']==ident and execution['status']==verified['status']=='PASS'
        D=batch.BATCH/ident;raw=read(D/'result.json');maps=read(D/'maps.json')
        cloud=read(ROOT/verified['mod2_certificate']);odd=read(ROOT/verified['modl_certificate'])
        assert batch.cert.hashed(D/'result.json')==execution['result_sha256']==verified['input_sha256']
        assert raw['maps_sha256']==batch.cert.hashed(D/'maps.json') and len(raw['charts'])==49 and len(maps['sample'])==2048
        assert raw['status']=='COMPLETE_DECLARED_ADAPTIVE_ATTEMPT' and raw['initial_dimension']==19
        assert cloud['rank_lower_bound']==raw['rank_lower_bound']==verified['rank_lower_bound']
        stage_records={}
        for name in ['maps','worker','replay','mod2-build','mod2-check','modl-build','modl-check']:
            q=read(D/(name+'.supervisor.json'));assert q['outcome']=='completed' and q['returncode']==0
            log=D/(name+'.log');paths.add(log);assert batch.cert.hashed(log)==q['log_sha256'];stage_records[name]=q
        for c in raw['charts']:
            assert c['search']['height_bound']==125000 and c['search']['timeout_seconds']==10
        completed=sum(c['search']['status']=='bounded_search_complete' for c in raw['charts'])
        cpu=[c['search']['search_cpu_ms'] for c in raw['charts']]
        seconds=sum(s['wall_seconds'] for s in stage_records.values());total+=seconds
        rows.append({'id':ident,'family':metadata['family'],'parameter':metadata['parameter'],
            'initial_rank_lower_bound':19,'rank_lower_bound':cloud['rank_lower_bound'],
            'odd_modulus_finite_ranks':{str(a['modulus']):a['finite_column_rank'] for a in odd['audits']},
            'attempted_boxes':49,'completed_boxes':completed,'retained_points':len(cloud['points']),
            'maps_seconds':stage_records['maps']['wall_seconds'],'point_seconds':stage_records['worker']['wall_seconds'],
            'local_verification_seconds':seconds-stage_records['maps']['wall_seconds']-stage_records['worker']['wall_seconds'],
            'gp_cpu_ms':sum(cpu) if all(type(c)==int for c in cpu) else None})
    earlier=[]
    for version in (1,2):
        folder=batch.LOCAL/('retained-native19-point-trial-v'+str(version))
        previous_protocol=read(folder/'protocol.json')
        assert all(batch.cert.hashed(ROOT/n)==h for n,h in {**previous_protocol['sources'],**previous_protocol['inputs'],**previous_protocol['seed_hashes']}.items())
        ledger=read(folder/'ledger.json');assert ledger['status']=='FAILED_OR_CENSORED' and not ledger['rows']
        exposure=0;cost=0
        for item in ledger['maps']:
            snapshot=read(folder/item['id']/'maps.json');assert not snapshot['rows']
            q=read(folder/item['id']/'maps.supervisor.json');cost+=q['wall_seconds']
            log=folder/item['id']/'maps.log';paths.add(log);assert batch.cert.hashed(log)==q['log_sha256']
            assert not (folder/item['id']/'result.json').exists()
        earlier.append({'version':version,'point_boxes':exposure,'preparation_supervised_seconds':cost,'status':'FAILED_OR_CENSORED'})
    control=read(batch.GATE);assert control['status']=='PASS' and control['rank_lower_bound']==28
    return {'schema':'elliptic-curves.retained-native19-followup-result.v3','status':'PASS_VERIFIED_DECLARED_ATTEMPTS',
        'sources':{str(q.relative_to(ROOT)):batch.cert.hashed(q) for q in sorted(paths)},'rows':rows,
        'attempted_boxes':98,'completed_boxes':sum(r['completed_boxes'] for r in rows),
        'total_supervised_seconds':total,'previous_failed_attempts':earlier,'trial_supervised_seconds_including_failed_preparations':total+sum(r['preparation_supervised_seconds'] for r in earlier),'accounting_boundary':'Supervised geometry, point workers, history, mod2/3/5 and exact map replays for this trial and prior V1/V2 preparations. Score/intake, method-control, report generation, later independent Sage and standalone proof costs are separate. This is not a complete all-session CPU measurement.','shared_geometry_replay_seconds':geo['wall_seconds'],
        'certified_added_directions':sum(r['rank_lower_bound']-19 for r in rows),
        'boundary':'Two fixed retained native19 curves, both complete maps before points,49 declared125000-height ten-second charts each, no rank stopping. Full histories, clouds and rational geometry replay. Lower bounds only; a null result is not exact rank, saturation, absence or a sensitivity theorem. No new parameter scan or automatic next wave.'}

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true');args=parser.parse_args();r=expected()
    if args.check:assert r==batch.cert.read(OUT)
    else:
        if OUT.exists():raise FileExistsError('preserve outcome')
        checkpoint(OUT,r)
    print(r['status'],r['completed_boxes'],'completed; gains',r['certified_added_directions'],r['rows'])
