#!/usr/bin/env python3
"""Freeze, execute and replay bounded blind chart-policy comparisons.

The policy input contains only exact curves, known points and a scoring Gram.
No missing/public target direction enters centre selection. Held-out identities
are fixed before any enumeration and are excluded from policy ranking.
"""
import argparse
from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import monotonic

from research_runtime.chart_policy import ChartPolicy, calibration_protocol, rank_calibration
from research_runtime.store import checkpoint, digest
from research_runtime.supervisor import capture, Limits
from run_mw_search import search_request
from pointed_quartic_search import ROOT

R17=ROOT/'elliptic-curves/data/r17_refresh_jump_ladder_blind_inputs_v1.json'
MW18=ROOT/'artifacts/generated-results/elkies-k3-r17-extreme-anchored-mw18-specializations-h1000-v1.json'


def freeze(path):
    if path.exists():raise FileExistsError('a frozen protocol is immutable')
    cases=[]; groups={17:[],18:[]}
    for row in json.loads(R17.read_text())['cases']:
        cases.append({'curve':row['short_model'],'points':row['generic_points'],'metric_gram':row['generic_height_gram']})
        groups[17].append(digest(cases[-1]))
    for row in json.loads(MW18.read_text())['candidates']:
        if not row['nagao']['is_certified_anchor']:continue
        points=row['specialized_points']['generic_R17']+[row['specialized_points']['cover_section']]
        # Euclidean scoring here is deliberately not called a canonical height.
        cases.append({'curve':row['raw_short_model'],'points':points,
                      'metric_gram':[[int(i==j) for j in range(18)] for i in range(18)]})
        groups[18].append(digest(cases[-1]))
    if list(map(len,groups.values()))!=[16,5]:raise ArithmeticError('blind control panel changed')
    keyed={digest(c):c for c in cases}
    held_out=[min(group) for group in groups.values()]
    policies=[ChartPolicy(chart_metric_kind=kind,chart_metric_weight=weight,diversity_window=2)
              for kind,weight in [('raw','1'),('metric','1/16'),('metric','1'),('metric','16')]]
    limits={'next_holes':12,'height':10000,'seconds_per_chart':'0.5','cvp_node_budget':100000,
            'cell_wall_seconds':240,'rss_bytes':1073741824,
            'extra_primes':[211,223,227,229,233,239,241,251]}
    sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in (R17,MW18)}
    protocol=calibration_protocol(panel=sorted(set(keyed)-set(held_out)),policies=policies,
                                 limits=limits,outcome_commitment=digest(sources),controls=held_out)
    record={'protocol':protocol,'cases':keyed,'inputs':sources,
            'panel_roles':{str(rank):ids for rank,ids in groups.items()},
            'blinding':'Outcome and target-point fields are excluded from every search request. These are historical controls, not unseen mathematical families.',
            'metric_boundary':'R17 uses its declared generic Gram; MW18 uses Euclidean scoring. Chart metric weights vary with identical centres within each case.'}
    checkpoint(path,record)
    return record


def run(document,output,verify=False):
    frozen=json.loads(document.read_text());protocol=frozen['protocol'];root=output.parent
    rows=[];measurements=[]
    order=protocol['panel']+protocol['held_out_controls']
    for case in order:
        for policy in protocol['policies']:
            key=digest(policy); stem=f'{case[:16]}-{key[:16]}'
            request={**frozen['cases'][case],'schema':'elliptic-curves.lazy-mw-search.v1','policy':policy,
                     **{k:protocol['limits'][k] for k in ('next_holes','height','seconds_per_chart','cvp_node_budget','extra_primes')}}
            path=root/'cells'/f'{stem}.json'; request_path=root/'requests'/f'{stem}.json'
            if verify:
                retained=json.loads(path.read_text())
                replay_path=root/'replay'/f'{stem}.json'
                result=capture([sys.executable,str(ROOT/'elliptic-curves/cas/run_mw_search.py'),
                    '--worker','--verify',str(path),'--output',str(replay_path)],
                    limits=Limits(protocol['limits']['cell_wall_seconds'],protocol['limits']['rss_bytes']),
                    log_path=root/'replay-logs'/f'{stem}.log')
                checked=json.loads(replay_path.read_text())
                supervision=result.supervision
            else:
                if not path.exists():
                    checkpoint(request_path,request)
                    result=capture([sys.executable,str(ROOT/'elliptic-curves/cas/run_mw_search.py'),
                        '--worker','--request',str(request_path),'--output',str(path)],
                        limits=Limits(protocol['limits']['cell_wall_seconds'],protocol['limits']['rss_bytes']),
                        log_path=root/'logs'/f'{stem}.log',check=False)
                    supervision=result.supervision
                    if result.returncode or not path.exists():
                        rows.append({'case':case,'policy':key,'status':'FAILED_CELL','supervision':supervision})
                        checkpoint(output,{'protocol_hash':protocol['protocol_hash'],'status':'INCOMPLETE','cells':rows})
                        print(f"CALIBRATION|cells={len(rows)}|FAILED",flush=True);continue
                checked=json.loads(path.read_text())
                if checked['request_hash']!=digest(request):raise ArithmeticError('cached cell differs from frozen request')
            if checked['cvp_status']!='COMPLETE_REQUEST' or len(checked['charts'])!=request['next_holes']:
                rows.append({'case':case,'policy':key,'status':'INCOMPLETE_CENTRES','result':str(path)})
                continue
            wall=sum(r['search']['wall_seconds'] for r in checked['charts'])
            row={'case':case,'policy':key,'protocol_hash':protocol['protocol_hash'],
                 'certified_independent_recoveries':checked['certified_rank_gain'],'wall_seconds':str(wall),
                 'complete_box_count':sum(r['search']['status']=='bounded_search_complete' for r in checked['charts']),
                 'status':checked['status'],'result':str(path.relative_to(root)),
                 'result_sha256':sha256(path.read_bytes()).hexdigest()}
            rows.append(row)
            if case in protocol['panel']:measurements.append(row)
            checkpoint(output,{'protocol_hash':protocol['protocol_hash'],'status':'RUNNING','cells':rows})
            print(f"CALIBRATION|cells={len(rows)}/{len(order)*len(protocol['policies'])}|gain={checked['certified_rank_gain']}",flush=True)
    complete=len(measurements)==len(protocol['panel'])*len(protocol['policies']) and all('certified_independent_recoveries' in r for r in rows)
    result={'schema':'elliptic-curves.chart-policy-sweep-result.v1','protocol_hash':protocol['protocol_hash'],
            'status':'COMPLETE' if complete else 'INCOMPLETE','cells':rows,
            'ranking':rank_calibration(protocol,measurements) if complete else None,
            'claim_boundary':'Bounded policy measurements only. Zero recovery and timeouts are not mathematical exclusions. A speed-only win does not establish the optimal recovery metric.'}
    checkpoint(output,result)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--freeze',type=Path);p.add_argument('--protocol',type=Path)
    p.add_argument('--output',type=Path);p.add_argument('--verify',action='store_true')
    a=p.parse_args()
    if a.freeze:
        result=freeze(a.freeze);print('FROZEN|'+result['protocol']['protocol_hash']);return
    if not a.protocol or not a.output:p.error('--protocol and --output are required for a sweep')
    run(a.protocol,a.output,a.verify)

if __name__=='__main__':main()
