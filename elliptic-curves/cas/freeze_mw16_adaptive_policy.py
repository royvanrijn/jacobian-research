#!/usr/bin/env python3
"""Freeze the uniform rule: adapt iff 0 < initial certified gain < threshold.

The threshold is a calibration scheduling choice, not a bound on E(Q).
Every active control must have completed both declared centre constructions.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]


def freeze(initial,trials,threshold):
    if initial['status']!='COMPLETE' or trials['status']!='COMPLETE' or threshold<=1:
        raise ArithmeticError('incomplete calibration inputs or invalid trigger')
    eligible={r['parent_id'] for r in initial['results'] if 0<r['exact_quotient_rank_recovered']<threshold}
    if {r['parent_id'] for r in trials['results']}!=eligible:
        raise ArithmeticError('adaptive trial set differs from the uniform gain trigger')
    rows=[]
    for before in initial['results']:
        if before['parent_id'] in eligible:
            row=next(r for r in trials['results'] if r['parent_id']==before['parent_id'])
            if row['status']!='COMPLETE': raise ArithmeticError('active adaptive control incomplete')
        else:
            row={**before,'settings':[],'ranking':{},
                'skipped_adaptive_reason':f'initial certified gain outside 1..{threshold-1}'}
        rows.append(row)
    return {**trials,'results':rows,'declared_budget':{**trials['declared_budget'],
        'parent_ids':[r['parent_id'] for r in rows],
        'adaptive_trigger':{'initial_gain_minimum':1,'initial_gain_exclusive_maximum':threshold},
        'active_adaptive_parent_ids':sorted(eligible)},
        'total_control_quotient_rank':sum(r['exact_quotient_rank_recovered'] for r in rows),
        'claim_boundary':trials['claim_boundary']+['The gain trigger schedules adaptive work; it is not a rank upper bound.']}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--initial',type=Path,required=True);p.add_argument('--trials',type=Path,required=True)
    p.add_argument('--threshold',type=int,default=8);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    result=freeze(json.loads(args.initial.read_text()),json.loads(args.trials.read_text()),args.threshold)
    result['inputs']={str(path.resolve().relative_to(ROOT)):sha256(path.read_bytes()).hexdigest()
        for path in (Path(__file__),args.initial,args.trials)}
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(f'FREEZE_ADAPTIVE|total={result["total_control_quotient_rank"]}|active={len(result["declared_budget"]["active_adaptive_parent_ids"])}')


if __name__=='__main__': main()
