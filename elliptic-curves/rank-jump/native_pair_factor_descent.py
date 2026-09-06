#!/usr/bin/env python3
"""Two fixed known-point pair models; reuse the existing conic/Jacobian descent."""
import argparse
from pathlib import Path
import subprocess
import retrospective as r
import minimal_native_block_carrier as pair

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'NATIVE_PAIR_FACTOR_DESCENT_PROTOCOL.json'
SOURCE=r.OUT/'rank_jump_native_genus_five_lift_inputs_v1.json'
PRIOR=r.OUT/'rank_jump_native_genus_five_lift_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-native-pair-factor-descent-v1'
OUTPUT=r.OUT/'rank_jump_native_pair_factor_descent_v1.json'


def input_path(i):return r.OUT/f'rank_jump_native_pair_factor_descent_{i}_inputs_v1.json'


def capture():
    src=r.read(SOURCE);lift=src['retained_lift']
    for i in range(2):
        r.write_new(input_path(i),{'schema':'rank-jump.native-pair-factor-descent-inputs.v1',
            'covers':[src['covers'][0],src['covers'][i+1]],
            'anchor':{'t':lift['t'],'u':lift['roots'][0],'v':lift['roots'][i+1]},
            'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (PROTOCOL,SOURCE,PRIOR,HERE/'minimal_native_block_carrier.py')}})


def worker(i):
    pair.INPUT=input_path(i);pair.PROTOCOL=PROTOCOL;pair.WORK=WORK/str(i)
    pair.worker()


def run():
    rows=[]
    for i in range(2):
        wd=WORK/str(i);wd.mkdir(parents=True,exist_ok=True)
        log=wd/'worker.log';execution=wd/'execution.json'
        if not log.exists():
            with log.open('x') as out:
                try:
                    p=subprocess.run(['/home/royvanrijn/.local/bin/sage','-python',str(Path(__file__).resolve()),'worker','--case',str(i)],stdout=out,stderr=out,timeout=60)
                    status={'status':'COMPLETE' if p.returncode==0 else 'FAILED','returncode':p.returncode}
                except subprocess.TimeoutExpired:status={'status':'TIMEOUT'}
            r.write_new(execution,status)
        row={'index':i,'execution':r.read(execution),'log':log.read_text()}
        for key in ('geometry','descent'):
            path=wd/f'{key}.json';row[key]=r.read(path) if path.exists() else {'status':'UNKNOWN'}
        rows.append(row);print(row['execution'],row['log'],flush=True)
    paths=[PROTOCOL,SOURCE,PRIOR,Path(__file__),HERE/'minimal_native_block_carrier.py',HERE/'retrospective.py']+[input_path(i) for i in range(2)]
    r.write_new(OUTPUT,{'schema':'rank-jump.native-pair-factor-descent.v1','pairs':rows,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths},
        'boundary':r.read(PROTOCOL)['boundary']})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','run','worker']);p.add_argument('--case',type=int);a=p.parse_args()
    if a.mode=='capture':capture()
    elif a.mode=='run':run()
    else:worker(a.case)
