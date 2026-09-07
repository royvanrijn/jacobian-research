#!/usr/bin/env python3
"""Authorize the frozen prospective replay only after exact control replay."""
import argparse
from hashlib import sha256
import gzip
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]


def build(initial,adaptive,verified):
    if verified['status']!='PASS_EXACT_REPLAY':
        raise ArithmeticError('independent exact replay has not passed')
    if initial['status']!='COMPLETE' or adaptive['status']!='COMPLETE':
        raise ArithmeticError('control campaign incomplete')
    rows=adaptive['results']
    if sorted(r['curve_id'] for r in rows)!=[398,400,401,542,548]:
        raise ArithmeticError('five-curve calibration panel changed')
    total=sum(r['exact_quotient_rank_recovered'] for r in rows)
    if total<54:
        raise ArithmeticError('fewer than 54 certified control directions; prospectives remain gated')
    result=[]
    for row in rows:
        before=next(r for r in initial['results'] if r['curve_id']==row['curve_id'])
        settings=[{k:s[k] for k in ('centre','specification','height')} for r in (before,row) for s in r['settings']]
        result.append({'curve_id':row['curve_id'],'exact_quotient_rank_recovered':row['exact_quotient_rank_recovered'],
            'settings':settings,'initial_gain':before['exact_quotient_rank_recovered']})
    return {'schema':'elliptic-curves.mw16-sensitivity-prospective-gate.v1','status':'COMPLETE',
        'total_control_quotient_rank':total,'results':result,
        'pipeline':'Frozen initial centre/coordinate setting, then the declared bounded five-bit adaptive policy when the frozen gain trigger holds.',
        'adaptive_trigger':adaptive.get('declared_budget',{}).get('adaptive_trigger'),
        'claim_boundary':'This gate measures detector sensitivity; it supplies no prospective rank evidence.'}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('initial','adaptive','verified','bundle','output'): p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args()
    initial=json.loads(args.initial.read_text()); adaptive=json.loads(args.adaptive.read_text())
    verified=json.loads(args.verified.read_text()); bundle=json.loads(gzip.decompress(args.bundle.read_bytes()))
    if sha256(args.bundle.read_bytes()).hexdigest()!=verified['bundle_sha256']:
        raise ArithmeticError('independent replay bundle checksum changed')
    for path in (args.initial,args.adaptive):
        name=str(path.resolve().relative_to(ROOT))
        if sha256(path.read_bytes()).hexdigest()!=bundle['files'][name]['sha256']:
            raise ArithmeticError('calibration differs from exactly replayed evidence')
        if not any(c['campaign']==name for c in verified['campaigns']):
            raise ArithmeticError('calibration was not independently replayed')
    payload=build(initial,adaptive,verified)
    paths=(Path(__file__),args.initial,args.adaptive,args.verified,args.bundle)
    payload['inputs']={str(path.resolve().relative_to(ROOT)):sha256(path.read_bytes()).hexdigest() for path in paths}
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
    print(f'MW16_PROSPECTIVE_GATE|PASS|control_directions={payload["total_control_quotient_rank"]}',flush=True)


if __name__=='__main__': main()
