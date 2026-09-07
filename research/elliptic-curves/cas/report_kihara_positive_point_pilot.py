#!/usr/bin/env python3
"""Completed six-fibre exposure, exact initial spans, gains and comparisons."""
import argparse,json
from pathlib import Path
import kihara_positive_point_pilot as batch
import certify_compact_r17_candidates as cert
import audit_inventory200_current_catalogue as catalogue
ROOT=batch.ROOT;ART=batch.ART;LOCAL=ROOT/'artifacts/local/elliptic-curves';OUT=ART/'kihara_positive_point_pilot_report_v1.json'
def compute():
    d=batch.BATCH;p=batch.protocol();assert cert.read(d/'ledger.json')['status']==cert.read(d/'verification-ledger.json')['status']=='PASS'
    paths={Path(__file__).resolve(),d/'protocol.json',d/'ledger.json',d/'verification-ledger.json',d/'verification-protocol.json'};stages=[]
    def supervision(path,expected='completed'):
        s=cert.read(path);assert s['outcome']==expected and (s['returncode']==0)==(expected=='completed');log=Path(s['log']);assert cert.hashed(log)==s['log_sha256'];paths.update([path,log]);stages.append({'stage':str(path.parent.relative_to(LOCAL))+'/'+path.name,'outcome':expected,'wall_seconds':s['wall_seconds']})
    for folder,expected in [('kihara-positive-fibre-intake-v1','completed'),('kihara-positive-seed-gaps-v1','completed'),('kihara-positive-seed-spans-v1','backend_failure'),('kihara-positive-seed-spans-v2','completed'),('kihara-positive-fibre-seeds-standalone-v1','completed')]:
        q=LOCAL/folder;supervision(q/'supervisor.json',expected);paths.update(q.glob('*.json'));pr=cert.read(q/'protocol.json')
        for name,h in pr['sources'].items():assert cert.hashed(ROOT/name)==h;paths.add(ROOT/name)
    seedproof=ART/'kihara_positive_fibre_seed_replay_v1.json';proof=ART/'kihara_positive_point_cloud_replay_v1.json';paths.update([seedproof,proof]);initial=cert.read(seedproof);independent=cert.read(proof);assert initial['status']==independent['status']=='PASS'
    dbpath=catalogue.D/'database.json';meta=catalogue.D/'metadata.json';index=ART/'new_high_rank_curve_index_v22.json';paths.update([dbpath,meta,index]);assert cert.hashed(dbpath)==cert.read(meta)['sha256'];db=cert.read(dbpath);inventory=cert.read(index);assert len(db['curves'])==620 and len(inventory['curves'])==201
    byj={};local={}
    for r in db['curves']:byj.setdefault(catalogue.j(r['ainvs']),[]).append(r)
    for r in inventory['curves']:local.setdefault(catalogue.j(r['curve']),[]).append(r)
    rows=[]
    for i,row in enumerate(p['rows']):
        folder=d/row['id'];seed=cert.read(folder/'seed.json');result=cert.read(folder/'result.json');cloudpath=ART/('kihara_positive_'+row['id'].replace('-','_')+'_cloud_v1.json');cloud=cert.read(cloudpath)
        assert result['status']=='COMPLETE_DECLARED_POINT_ATTEMPT' and len(result['charts'])==49 and cloud['status']=='PASS'
        assert cloud['source_sha256']==cert.hashed(folder/'result.json') and cloud['points']==result['retained_points'] and cloud['curve']==seed['curve']
        init=initial['rows'][i];full=independent['rows'][i];assert init['id']==full['id']==row['id'] and init['exact_displayed_span_rank']==row['initial_rank']==full['initial_rank'] and full['rank_lower_bound']==cloud['rank_lower_bound']
        assert [a['rank_lower_bound'] for a in cloud['audits']]==[a['rank_lower_bound'] for a in full['audits']]
        j=catalogue.j(seed['curve']);completed=sum(c['search']['status']=='bounded_search_complete' for c in result['charts'])
        rows.append({**row,'curve':seed['curve'],'points':len(cloud['points']),'rank_lower_bound':cloud['rank_lower_bound'],'discovered_rank_gain_lower_bound':cloud['rank_lower_bound']-row['initial_rank'],'completed_boxes':completed,'attempted_boxes':49,'incomplete_boxes':49-completed,'gp_cpu_milliseconds':sum(c['search']['search_cpu_ms'] or 0 for c in result['charts']),
            'rank_moduli':[a['modulus'] for a in full['audits']],
            'catalogue_matches':[r['id'] for r in byj.get(j,[]) if cert.isomorphic(seed['curve'],r['ainvs'])],
            'inventory_matches':[r['id'] for r in local.get(j,[]) if cert.isomorphic(seed['curve'],r['curve'])]})
        for stage in ('maps','worker','replay','cloud'):supervision(folder/(stage+'.supervisor.json'))
        paths.update([folder/'seed.json',folder/'maps.json',folder/'result.json',cloudpath])
    assert len(set(catalogue.j(r['curve']) for r in rows))==6
    assert [r['rank_lower_bound'] for r in rows]==[7,7,4,9,14,12]
    supervision(d/'geometry.supervisor.json');supervision(d/'standalone/supervisor.json');paths.add(d/'standalone/protocol.json')
    for name in ('kihara_positive_point_cloud_bundle_v1.json','kihara_positive_fibre_seed_bundle_v1.json'):
        path=ART/name;paths.add(path)
        for n,h in cert.read(path)['sources'].items():assert cert.hashed(ROOT/n)==h;paths.add(ROOT/n)
    for folder,script,bundle in [(d/'standalone','verify_kihara_positive_point_clouds.sage','kihara_positive_point_cloud_bundle_v1.json'),(LOCAL/'kihara-positive-fibre-seeds-standalone-v1','verify_kihara_positive_fibre_seeds.sage','kihara_positive_fibre_seed_bundle_v1.json')]:
        for path in [batch.CAS/script,ART/bundle]:assert path.read_bytes()==(folder/path.name).read_bytes();paths.add(path)
    for n,h in cert.read(d/'verification-protocol.json')['sources'].items():assert cert.hashed(ROOT/n)==h;paths.add(ROOT/n)
    return {'schema':'kihara-positive-point-pilot-report.v1','status':'PASS','rows':rows,'attempted_boxes':294,'completed_boxes':sum(r['completed_boxes'] for r in rows),'discovered_rank_gain_lower_bound':sum(r['discovered_rank_gain_lower_bound'] for r in rows),'new_Q_distinct_parents_in_this_stage':0,'inventory_additions':0,'stages':stages,'total_supervised_seconds_including_failed_span_attempt':sum(s['wall_seconds'] for s in stages),'sources':{str(path.relative_to(ROOT)):cert.hashed(path) for path in sorted(paths)},
        'scope':'Equal49-box fixed point exposure on six fibres of three previously certified new positive-ratio Kihara parents. Exact displayed specialization spans are6,7,4,9,12,12; two unit fibres lose generic directions, and both first-parent fibres have rational2-torsion. Empty-state explicit-centre backend searches avoid the adaptive no-two-torsion gate, without weakening it. Frozen numerical centre proposals have exact parity/norm/group/map checks. The rank4 seed uses additional representatives in existing cosets, not new parity classes. Complete returned clouds independently prove final lower bounds7,7,4,9,14,12 at3and5. Three discovered directions occur on two curves; no exact whole-curve rank, parent-superiority, new high-rank inventory or near-record claim. Pinned620 catalogue and201 inventory comparisons are post-search, not literature-wide novelty. No larger population or automatic next wave follows. Initial intake and the failed numerical-span source attempt are charged; previous parent-construction costs are not counted again.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();r=compute()
    if a.check:assert r==cert.read(OUT)
    else:
        with OUT.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print('PASS',r['completed_boxes'],'completed boxes;',r['discovered_rank_gain_lower_bound'],'certified directions;',r['total_supervised_seconds_including_failed_span_attempt'],'supervised seconds')
    print('Pinned matches',[(a['id'],a['catalogue_matches'],a['inventory_matches']) for a in r['rows']])
