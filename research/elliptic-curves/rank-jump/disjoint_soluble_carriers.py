#!/usr/bin/env python3
"""One complementary pair in the already certified modulo-23 obstruction."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import subprocess
import global_pair_solubility as g
import retrospective as r
from cover_experiment import sqrtq
from local_solubility_blocks import value

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'DISJOINT_SOLUBLE_CARRIERS_PROTOCOL.json'
INPUT=r.OUT/'rank_jump_disjoint_soluble_carriers_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_disjoint_soluble_carriers_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-disjoint-soluble-carriers-v1'


def capture():
    old=r.read(g.OLD_INPUT);results=r.read(g.OLD_RESULT);covers=old['covers']
    common,partner=8,1;mask=(1<<common)|(1<<partner)
    group=next(x for x in old['observed_groups'] if x['mask']>>common&1)
    t=F(group['published_parameter']);f=covers[common]['form']
    u=sqrtq(F(value(f,t.numerator,t.denominator),t.denominator**2));assert u is not None
    witnesses=[]
    for p in results['places']:
        w=next(w for w in p['maximal_proved_masks'] if w['mask']&mask==mask)
        witnesses.append({'prime':p['prime'],'base_point':w['rational_base_point']})
    c={'id':'complementary_AD','labels':[covers[common]['label'],covers[partner]['label']],
       'forms':[f,covers[partner]['form']],'single_conic_anchor':[str(t),str(u)],'old_local_witnesses':witnesses}
    bc=next(x for x in old['observed_groups'] if x['mask']&12==12)
    r.write_new(INPUT,{'schema':'rank-jump.disjoint-soluble-carriers-inputs.v1','cases':[c],
                       'known_BC_pair':{'labels':[covers[i]['label'] for i in [2,3]],'forms':[covers[i]['form'] for i in [2,3]],'parameter':bc['published_parameter']},
                       'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (PROTOCOL,g.OLD_INPUT,g.OLD_RESULT,Path(g.__file__))}})


def run():
    WORK.mkdir(parents=True,exist_ok=True);name='complementary_AD';log=WORK/(name+'.log')
    if not log.exists():
        with log.open('x') as f:
            try:
                p=subprocess.run(['sage','-python',str(Path(__file__).resolve()),'worker'],cwd=r.ROOT,stdout=f,stderr=f,timeout=60)
                status='COMPLETE' if p.returncode==0 else 'FAILED'
            except subprocess.TimeoutExpired:status='TIMEOUT'
        r.write_new(WORK/(name+'_execution.json'),{'status':status})
    row={'id':name,'execution':r.read(WORK/(name+'_execution.json')),'log':log.read_text()}
    for stage in ('geometry','local','descent'):
        path=WORK/(name+'_'+stage+'.json');row[stage]=r.read(path) if path.exists() else {'status':'UNKNOWN'}
    r.write_new(OUTPUT,{'schema':'rank-jump.disjoint-soluble-carriers.v1','rows':[row],
                        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (INPUT,PROTOCOL,Path(__file__),Path(g.__file__),HERE/'retrospective.py',HERE/'local_solubility_blocks.py')},
                        'boundary':r.read(PROTOCOL)['boundary']})
    print(row['execution']);print(row['log'])


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','run','worker']);a=p.parse_args()
    if a.mode=='capture':capture()
    elif a.mode=='run':run()
    else:
        g.INPUT=INPUT;g.PROTOCOL=PROTOCOL;g.WORK=WORK
        g.worker(0)
