#!/usr/bin/env python3
"""Post-outcome lookup in already completed score tables; no new traces."""
import argparse
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/corrected-mw16-higher-scores-v1'
OUT=ART/'strata60_retained_score_positions_v1.json'


def expected():
    selection=ART/'retained_mw16_score_strata_selection_v1.json'
    report_path=ART/'retained_mw16_score_strata_experiment_v1.json'
    audit_path=ART/'retained_mw16_score_strata_accounting_replay_v1.json'
    original_path=ART/'corrected_mw16_higher_selection_v1.json'
    selected=cert.read(selection);report=cert.read(report_path);audit=cert.read(audit_path);old=cert.read(original_path)
    if selected['status']!='PASS_FROZEN60_MATCHED_SELECTION' or report['status']!='COMPLETE_FIXED_COMPARISON' or audit['status']!='PASS':raise ArithmeticError('terminal independently audited comparison required')
    if old['protocol_sha256']!=cert.hashed(D/'protocol.json') or old['scores_sha256']!=cert.hashed(D/'result.json'):raise ArithmeticError('original retained scores changed')
    p=cert.read(D/'protocol.json');scores=cert.read(D/'result.json')
    if scores['status']!='COMPLETE_FROZEN10240' or len(scores['rows'])!=10240 or len(p['rows'])!=10240:raise ArithmeticError('complete old scalar population required')
    if cert.read(D/'controller/ledger.json')['status']!='PASS':raise ArithmeticError('original scalar replay required')
    pools={}
    for row,score in zip(p['rows'],scores['rows']):
        if row['id']!=score['id']:raise ArithmeticError('original score roster differs')
        pools.setdefault((row['family'],row['band']),[]).append({**row,**score})
    lookup={}
    for block,pool in pools.items():
        if len(pool)!=1024:raise ArithmeticError('equal original scalar pools required')
        pool.sort(key=lambda r:(-r['combined_late_units'],-r['combined_late_good'],r['denominator'],r['numerator']))
        for i,r in enumerate(pool):lookup[(*block,r['parameter'])]=(i+1,r,pool[5]['combined_late_units'])
    original_ids={r['id'] for r in old['selected']};rows=[]
    for row,outcome in zip(selected['selected'],report['rows']):
        if row['id']!=outcome['id']:raise ArithmeticError('completed comparison roster differs')
        record={k:row[k] for k in ('id','family','band','arm','parameter','retained_rank')}
        record.update(rank_lower_bound=outcome['rank_lower_bound'],certified_gain=outcome['certified_gain'],original_scalar_row=None)
        match=lookup.get((row['family'],row['band'],row['parameter']))
        if match:
            rank,prior,cutoff=match
            if not cert.isomorphic(row['model'],prior['model']):raise ArithmeticError('same parameter has different curve')
            record['original_scalar_row']={'id':prior['id'],'initial_combined_sign_pool_rank':prior['retained_index']+1,
                'late_rank_in_1024':rank,'late_selection_units':prior['combined_late_units'],
                'sixth_late_selection_units':cutoff,'late_score_gap_to_sixth_units':cutoff-prior['combined_late_units'],
                'original_finalist':prior['id'] in original_ids}
        rows.append(record)
    paths=[Path(__file__).resolve(),Path(cert.__file__),selection,report_path,audit_path,original_path,D/'protocol.json',D/'result.json',D/'controller/ledger.json']
    return {'schema':'elliptic-curves.strata60-retained-score-positions.v1','status':'PASS',
        'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'rows':rows,
        'overlap_with_original_scalar_population':sum(r['original_scalar_row'] is not None for r in rows),
        'claim_boundary':'Retrospective accounting for all60 completed matched curves against existing initial/late score tables. No additional prime trace, validation-prime result, parameter or point search is read or computed. The rank labels are joined only after the independently audited comparison ends. The initial rank refers to the signed retained block; the old initial and late positions use its two-sign family/band pool. A late shortlist miss locates a retained candidate; it does not prove another score rule, cutoff or point-budget allocation optimal.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=d:raise ArithmeticError('retrospective lookup replay differs')
    else:
        if OUT.exists():raise FileExistsError('preserve retrospective lookup')
        checkpoint(OUT,d)
    print('RETAINED SCORE POSITION AUDIT PASS',d['overlap_with_original_scalar_population'],'OF60 HAD EXISTING LATE SCORES',flush=True)
