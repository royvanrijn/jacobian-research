#!/usr/bin/env python3
"""Freeze the common initial setting with greatest certified five-curve gain.

Ties use height then the setting key. Runtime and prospective outcomes do not
enter selection. The input must have completed the same menu on all five
controls; each retained setting carries its full exact group certificate.
"""
from hashlib import sha256
import argparse
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]


def freeze(source):
    if source['status']!='COMPLETE' or source['declared_budget']['mode']!='initial':
        raise ArithmeticError('initial calibration menu is incomplete')
    results=source['results']
    if sorted(r['curve_id'] for r in results)!=[398,400,401,542,548]:
        raise ArithmeticError('five-curve control panel changed')
    common=set(s['key'] for s in results[0]['settings'])
    for r in results:
        common.intersection_update(s['key'] for s in r['settings'])
    ranking=[]
    for key in common:
        settings=[next(s for s in r['settings'] if s['key']==key) for r in results]
        if any(s['classification']['status']!='PASS_BASIS_EQUALS_DISCOVERED_GROUP' for s in settings):
            raise ArithmeticError('unclassified rank score')
        ranking.append((-sum(s['exact_quotient_rank_recovered'] for s in settings),settings[0]['height'],key))
    if not ranking: raise ArithmeticError('no common completed setting')
    ranking.sort(); selected=ranking[0][2]
    rows=[]
    for r in results:
        setting=next(s for s in r['settings'] if s['key']==selected)
        rows.append({**r,'settings':[setting],'current_basis':setting['current_basis'],
            'classification':setting['classification'],
            'exact_quotient_rank_recovered':setting['exact_quotient_rank_recovered']})
    setting=rows[0]['settings'][0]
    payload={**source,'results':rows,
        'declared_budget':{**source['declared_budget'],'specifications':[setting['specification']],
            'centres':[setting['centre']],'heights':[setting['height']]},
        'selection_rule':'maximize sum of exact quotient ranks on five controls; ties by height then key; runtime excluded',
        'selection_ranking':[{'key':key,'total_control_quotient_rank':-negative,'height':height} for negative,height,key in ranking],
        'total_control_quotient_rank':-ranking[0][0]}
    return payload


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True); p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    source=json.loads(args.input.read_text()); payload=freeze(source)
    payload['inputs']={str(path.resolve().relative_to(ROOT)):sha256(path.read_bytes()).hexdigest()
        for path in (args.input,Path(__file__),ROOT/'elliptic-curves/data/icarm_mw16_parent_ladder_blind_inputs_v1.json')}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    temporary=args.output.with_suffix('.tmp'); temporary.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n'); temporary.replace(args.output)
    print(f'FREEZE_SENSITIVITY|setting={payload["selection_ranking"][0]["key"]}|gain={payload["total_control_quotient_rank"]}',flush=True)


if __name__=='__main__': main()
