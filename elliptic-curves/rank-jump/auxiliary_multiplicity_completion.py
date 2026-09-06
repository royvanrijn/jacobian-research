#!/usr/bin/env python3
"""Read-only adapter for the frozen prime set's bad-reduction miss."""
import argparse
from pathlib import Path
import subprocess
import retrospective as r
import auxiliary_elliptic_multiplicity as producer

PROTOCOL=Path(__file__).with_name('AUXILIARY_MULTIPLICITY_COMPLETION_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_auxiliary_multiplicity_completion_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-auxiliary-multiplicity-completion-v1'


def compute():
    assert r.read(producer.OUTPUT)['rows'][1]['status']=='UNKNOWN'
    producer.PROTOCOL=PROTOCOL
    row=producer.compute(1)
    return {'schema':'rank-jump.auxiliary-multiplicity-completion.v1',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
            (Path(__file__),PROTOCOL,producer.OUTPUT)},'row':row}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker','check']);a=p.parse_args()
    if a.mode=='worker':r.write_new(WORK/'checkpoint.json',compute())
    elif a.mode=='check':assert r.read(OUTPUT)==compute();print('PASS swapped multiplicity completion')
    else:
        WORK.mkdir(parents=True,exist_ok=True);path=WORK/'checkpoint.json'
        if not path.exists():
            with (WORK/'worker.log').open('x') as log:
                try:
                    proc=subprocess.run(['sage','-python',str(Path(__file__).resolve()),'worker'],
                        cwd=r.ROOT,stdout=log,stderr=log,timeout=30)
                    reason=None if proc.returncode==0 else 'worker failure'
                except subprocess.TimeoutExpired:reason='30-second timeout'
                if reason and not path.exists():r.write_new(path,{'status':'UNKNOWN','reason':reason})
        data=r.read(path);r.write_new(OUTPUT,data);print(data.get('row',{}).get('status',data.get('status')))
