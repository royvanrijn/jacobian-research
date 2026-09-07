#!/usr/bin/env python3
"""Bound the claims of the six-parent calibration by verified exposure."""
import argparse
from pathlib import Path
import mestre_parent_calibration as batch
import certify_compact_r17_candidates as cert
import audit_inventory200_current_catalogue as catalogue
from research_runtime.store import checkpoint
OUT=batch.ART/'mestre_parent_calibration_v1.json'

def compute():
    p=batch.protocol();d=batch.BATCH;ledger=cert.read(d/'ledger.json');verified=cert.read(d/'verification-ledger.json')
    independent=batch.ART/'mestre_parent_calibration_sage_replay_v1.json'
    span=batch.ART/'mestre_parent_calibration_input_span_v1.json'
    seed_clouds=batch.ART/'mestre_parent_calibration_seed_clouds_v1.json'
    checks=cert.read(d/'post-verification-ledger.json')
    if any(x['status']!='PASS' for x in (ledger,verified,checks,cert.read(independent),cert.read(span),cert.read(seed_clouds))):raise ArithmeticError('complete independent audits required')
    for stage in checks['stages']:
        if stage['status']!='PASS' or cert.hashed(batch.ROOT/stage['source'])!=stage['source_sha256']:raise ArithmeticError('independent audit source changed')
    if [r['id'] for r in ledger['rows']]!=[r['id'] for r in p['rows']] or len(verified['rows'])!=12:raise ArithmeticError('fixed12 roster differs')
    database=catalogue.D/'database.json';metadata=catalogue.D/'metadata.json';db=cert.read(database)
    if cert.hashed(database)!=cert.read(metadata)['sha256'] or db['count']!=620 or len(db['curves'])!=620:raise ArithmeticError('pinned catalogue differs')
    if len(db['curves'])!=620:raise ArithmeticError('620 pinned equations required')
    index=batch.ART/'new_high_rank_curve_index_v22.json';inventory=cert.read(index)
    byj={};local={}
    for r in db['curves']:byj.setdefault(catalogue.j(r['ainvs']),[]).append(r)
    for r in inventory['curves']:local.setdefault(catalogue.j(r['curve']),[]).append(r)
    paths=[Path(__file__).resolve(),Path(batch.__file__),Path(catalogue.__file__),database,metadata,index,independent,span,seed_clouds,
      d/'protocol.json',d/'intake-protocol.json',d/'intake-ledger.json',d/'ledger.json',d/'verification-ledger.json',d/'post-verification-ledger.json']
    rows=[];seen=set();audit_by_id={r['id']:r for r in verified['rows']}
    for row in p['rows']:
        folder=d/row['id'];seed=cert.read(folder/'seed.json');result=cert.read(folder/'result.json');v=audit_by_id[row['id']]
        if result['status']!='COMPLETE_DECLARED_ADAPTIVE_ATTEMPT' or len(result['charts'])!=49:raise ArithmeticError('all49 chart attempts required')
        maps=cert.read(folder/'maps.json')
        if len(maps['sample'])!=2047 or len(maps['rows'])!=49:raise ArithmeticError('fixed geometry differs')
        finite=cert.read(batch.ROOT/v['mod2_certificate']);odd=cert.read(batch.ROOT/v['modl_certificate'])
        if cert.hashed(batch.ROOT/v['mod2_certificate'])!=v['mod2_sha256'] or cert.hashed(batch.ROOT/v['modl_certificate'])!=v['modl_sha256']:raise ArithmeticError('full-cloud proof changed')
        lower=finite['rank_lower_bound'];bounds={str(a['modulus']):a['finite_column_rank'] for a in odd['audits']}
        if lower<result['rank_lower_bound']:raise ArithmeticError('full-cloud rank lost a worker direction')
        if any(b!=lower for b in bounds.values()):raise ArithmeticError('unequal finite ranks need separate interpretation')
        j=catalogue.j(seed['curve']);seen.add(j)
        timings={}
        for stage in ('intake','maps','worker','replay','mod2-build','mod2-check','modl-build','modl-check'):
            path=folder/(stage+'.supervisor.json');s=cert.read(path)
            if s['outcome']!='completed' or s['returncode']!=0:raise ArithmeticError('completed stage required')
            timings[stage]=s['wall_seconds'];paths.append(path)
        completed=sum(c['search']['status']=='bounded_search_complete' for c in result['charts'])
        rows.append({**row,'model':seed['curve'],'model_coefficient_bits':seed['model_coefficient_bits'],
          'supplied_image_subgroup_rank':11,'rank_lower_bound':lower,'discovered_rank_gain':lower-11,
          'retained_points':len(finite['points']),'odd_modulus_lower_bounds':bounds,
          'attempted_boxes':49,'completed_boxes':completed,'incomplete_boxes':49-completed,
          'gp_cpu_milliseconds':sum(c['search']['search_cpu_ms'] or 0 for c in result['charts']),
          'supervised_stage_seconds':timings,
          'catalogue_matches':[r['id'] for r in byj.get(j,[]) if cert.isomorphic(seed['curve'],r['ainvs'])],
          'inventory_matches':[r['id'] for r in local.get(j,[]) if cert.isomorphic(seed['curve'],r['curve'])],
          'rank_certificate_path':v['mod2_certificate']})
        paths.extend([folder/'seed.json',folder/'maps.json',folder/'result.json',batch.ROOT/v['mod2_certificate'],batch.ROOT/v['modl_certificate']])
    groups=[]
    for t in ('1','1009/101'):
        chosen=[r for r in rows if r['fibre_T']==t];gain=sum(r['discovered_rank_gain'] for r in chosen)
        worker=sum(r['supervised_stage_seconds']['worker'] for r in chosen)
        inclusive=sum(sum(r['supervised_stage_seconds'].values()) for r in chosen)
        groups.append({'fibre_T':t,'fibres':len(chosen),'discovered_rank_gains':gain,
          'completed_boxes':sum(r['completed_boxes'] for r in chosen),'point_worker_seconds':worker,
          'per_row_inclusive_seconds':inclusive,'gains_per_point_worker_second':gain/worker,
          'gains_per_per_row_inclusive_second':gain/inclusive})
    total=sum(sum(r['supervised_stage_seconds'].values()) for r in rows)
    total+=sum(r['supervision']['wall_seconds'] for r in checks['stages'] if r['name']!='cloud-audits')
    return {'schema':'elliptic-curves.mestre-parent-calibration.v1','status':'PASS','rows':rows,'paired_fibre_groups':groups,
      'parents':6,'fibres':12,'distinct_j_count':len(seen),'attempted_boxes':sum(r['attempted_boxes'] for r in rows),
      'completed_boxes':sum(r['completed_boxes'] for r in rows),'discovered_rank_gains':sum(r['discovered_rank_gain'] for r in rows),
      'supervised_seconds_including_independent_audits':total,'catalogue_count':620,'inventory_count':len(inventory['curves']),
      'sources':{str(q.relative_to(batch.ROOT)):cert.hashed(q) for q in paths},
      'scope':'Fixed twelve-fibre feasibility/visibility pilot across six verified Q-distinct parents, with49 charts each. Exact supplied-image relations prove the original14-image clouds span rank11; gains increase that subgroup rank. Exact generic ranks remain unknown, so these are not certified specialization jumps. Completed boxes trust pinned PARI transcripts. Comparison is within parent at two declared base parameters, not a height-matched score study or a parent superiority test. No scores or validation-prime scoring were used. Catalogue comparison is post-search and limited to the pinned620 equations and201 inventory equations, not literature-wide novelty. Lower bounds only; no record or near-record claim. No automatic follow-up campaign.'}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');a=parser.parse_args();result=compute()
    if a.check:assert result==cert.read(OUT)
    else:
        if OUT.exists():raise FileExistsError('preserve parent calibration report')
        checkpoint(OUT,result)
    print('PASS',result['completed_boxes'],'boxes;',result['discovered_rank_gains'],'certified gains;',result['distinct_j_count'],'distinct j',flush=True)
