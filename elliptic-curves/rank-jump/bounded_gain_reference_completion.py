#!/usr/bin/env python3
"""Complete four inherited Artin entries; join conditional labels afterward."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import bounded_gain_reference as source
import fresh_governing_completion as repair

PROTOCOL=Path(__file__).with_name('BOUNDED_GAIN_REFERENCE_COMPLETION_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_bounded_gain_reference_completion_v1.json'
REPAIR=source.WORK/'reference-00-repair.json'
ADAPTER=source.WORK/'repair-input.json'


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
            (Path(__file__),PROTOCOL,source.OUTPUT,source.INPUT,Path(source.__file__),Path(repair.__file__),repair.base.LOCAL)}


def worker():
    from sage.all import pari
    source.configure();repair.PROTOCOL=PROTOCOL;repair.base.OUTPUT=ADAPTER
    pari.allocatemem(64000000,r.read(PROTOCOL)['limits']['pari_stack_bytes'],silent=True)
    old=r.read(source.OUTPUT)['stages'];adapter=r.read(ADAPTER)
    assert adapter=={'rows':[{'token':source.TOKEN,'factor':old['factor'],'local':old['local']}]}
    value=repair.repair(source.TOKEN,old['local'])
    r.write_new(REPAIR,{'bindings':bindings(),**value})


def capture():
    old=r.read(source.OUTPUT);stages=old['stages']
    if not ADAPTER.exists():r.write_new(ADAPTER,{'rows':[{'token':source.TOKEN,'factor':stages['factor'],'local':stages['local']}]})
    if not REPAIR.exists():
        error=None
        with (source.WORK/'repair.log').open('x') as log:
            try:
                p=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker'],stdout=log,stderr=log,
                                 timeout=r.read(PROTOCOL)['limits']['worker_seconds'])
                if p.returncode:error='Worker failure'
            except subprocess.TimeoutExpired:error='Bounded worker timeout'
        if error:r.write_new(REPAIR,{'bindings':bindings(),'status':'UNKNOWN','reason':error})
    result=r.read(REPAIR);assert result['bindings']==bindings()
    if result['status']=='PASS_REPAIRED':
        stages['local']={**stages['local'],**{k:v for k,v in result.items() if k!='bindings'}}
    r.write_new(OUTPUT,{'schema':'rank-jump.bounded-gain-reference-completion.v1',
        'token':source.TOKEN,'stages':stages,'repair':result,'bindings':bindings()})
    print(result['status'],result.get('minus_twist_CT_rank'),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['worker','capture','report']);args=p.parse_args()
    if args.mode=='report':source.OUTPUT=OUTPUT;source.report()
    else:globals()[args.mode]()
