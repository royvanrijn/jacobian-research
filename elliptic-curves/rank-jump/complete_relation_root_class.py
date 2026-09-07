#!/usr/bin/env python3
"""Complete the unchanged relation-root computation with JSON-safe integers."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import retrospective as r
import relation_root_class as prior

OUTPUT=r.OUT/'rank_jump_relation_root_class_v2.json'
WORK=r.ROOT/'artifacts/local/rank-jump-relation-root-class-v2'


def compute():
    from sage.rings.integer import Integer
    def encode(value):
        if isinstance(value,Integer):return int(value)
        raise TypeError(type(value).__name__)
    result=json.loads(json.dumps(prior.compute(),default=encode))
    result['schema']='rank-jump.relation-root-class.v2'
    result['bindings'].update(prior.bindings([Path(__file__)]))
    result['completion_reason']='v1 arithmetic completed but JSON encoding rejected a Sage integer. Identical candidate, primes and arithmetic; convert Sage integers to Python integers only.'
    return result


def capture():
    WORK.mkdir(parents=True,exist_ok=True);path=WORK/'worker.json'
    if not path.exists():
        reason=None
        with (WORK/'worker.log').open('x') as log:
            try:
                p=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker'],stdout=log,stderr=log,timeout=30)
                if p.returncode:reason='Worker failure'
            except subprocess.TimeoutExpired:reason='Bounded worker timeout'
        if reason:
            r.write_new(OUTPUT,{'status':'UNKNOWN','reason':reason});return
    result=r.read(path);r.write_new(OUTPUT,result)
    print(result['candidate_status'],[(row['p'],row['f'],row['projection_valuation']) for row in result['rows'] if row['projection_valuation']%2],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker','check']);a=p.parse_args()
    if a.mode=='worker':r.write_new(WORK/'worker.json',compute())
    elif a.mode=='check':assert compute()==r.read(OUTPUT);print('PASS completed relation-root replay')
    else:capture()
