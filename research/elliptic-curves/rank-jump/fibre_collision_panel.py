#!/usr/bin/env python3
"""Uniform bounded collision arithmetic on all nontrivial census blocks."""
import argparse
from itertools import combinations
from math import prod
from pathlib import Path
import subprocess
import retrospective as r

HERE=Path(__file__).resolve().parent
INPUT=r.OUT/'rank_jump_fibre_discrimination_v1.json'
PROTOCOL=HERE/'FIBRE_COLLISION_PANEL_PROTOCOL.json'
OUTPUT=r.OUT/'rank_jump_fibre_collision_panel_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-fibre-collision-panel-v1'


def targets():return [(key,b) for key,b in sorted(r.read(INPUT)['blocks'].items()) if b['compatible_cover_count']>=2]
def ev(q,x):return q[0]+q[1]*x+q[2]*x*x


def local(block):
    qs=[a['primitive_form'] for a in block['compatible_forms']];n=len(qs)
    potential=[];actual=[];collisions=[];witnesses=[]
    for p in [2]+r.primes(97):
        for chart,root in [('finite',x) for x in range(p)]+[('infinity',0)]:
            forms=qs if chart=='finite' else [q[::-1] for q in qs]
            inds=[i for i,q in enumerate(forms) if ev(q,root)%p==0]
            if len(inds)<2:continue
            masks=[(1<<i)|(1<<j) for i,j in combinations(inds,2)];potential.extend(masks)
            collisions.append({'prime':p,'chart':chart,'root':root,'indices':inds})
            if p<5 or len(inds)!=2:continue
            a=[ev(forms[i],root)//p%p for i in inds];d=[(forms[i][1]+2*forms[i][2]*root)%p for i in inds]
            if not all(d) or (a[0]*d[1]-a[1]*d[0])%p==0:continue
            for s in range(min(p,16)):
                t=root+p*s;values=[ev(q,t) for q in forms]
                units=[(v//p if i in inds else v)%p for i,v in enumerate(values)]
                unit=prod(units)%p
                if unit and pow(unit,(p-1)//2,p)==1:
                    mask=sum(1<<i for i in inds);actual.append(mask)
                    witnesses.append({'prime':p,'chart':chart,'root':root,'s':s,'mask':mask,'unit':unit});break
    return {'prime_bound':97,'collision_cells':collisions,'odd_pair_witnesses':witnesses,
            'potential_pair_span_interval':[r.rank(potential),n-1],
            'realizable_defect_span_interval':[r.rank(actual),n-1],
            'potential_pair_span_exact':r.rank(potential) if r.rank(potential)==n-1 else 'UNKNOWN',
            'realizable_defect_span_exact':r.rank(actual) if r.rank(actual)==n-1 else 'UNKNOWN'}


def factor_worker(index):
    from sage.all import ZZ
    key,block=targets()[index];wd=WORK/str(index)
    for j,pair in enumerate(block['pair_resultants']):
        path=wd/f'pair_{j}.json'
        if path.exists():continue
        n=abs(ZZ(pair['resultant']));fac=list(n.factor())
        assert prod(p**e for p,e in fac)==n and all(p.is_prime(proof=True) for p,e in fac)
        r.write_new(path,{'indices':pair['indices'],'resultant':pair['resultant'],'factors':[[str(p),int(e)] for p,e in fac]})


def run():
    rows=[]
    for i,(key,block) in enumerate(targets()):
        wd=WORK/str(i);wd.mkdir(parents=True,exist_ok=True)
        loc=wd/'local.json'
        if not loc.exists():r.write_new(loc,local(block))
        log=wd/'factor.log';state=wd/'execution.json'
        if not log.exists():
            with log.open('x') as out:
                try:
                    p=subprocess.run(['/home/royvanrijn/.local/bin/sage','-python',str(Path(__file__).resolve()),'factor','--index',str(i)],stdout=out,stderr=out,timeout=5)
                    execution={'status':'COMPLETE' if p.returncode==0 else 'FAILED','returncode':p.returncode}
                except subprocess.TimeoutExpired:execution={'status':'TIMEOUT'}
            r.write_new(state,execution)
        factors=[r.read(wd/f'pair_{j}.json') for j in range(len(block['pair_resultants'])) if (wd/f'pair_{j}.json').exists()]
        primes=sorted({int(p) for f in factors for p,e in f['factors']});complete=len(factors)==len(block['pair_resultants'])
        row={'block_key':key,'cover_count':block['compatible_cover_count'],'local':r.read(loc),'factor_execution':r.read(state),
             'pair_factorizations':factors,'pair_factor_count':len(factors),'pair_count':len(block['pair_resultants']),
             'full_support_status':'COMPLETE' if complete else 'UNKNOWN','certified_support_primes':list(map(str,primes)),
             'collision_support_size':len(primes) if complete else 'UNKNOWN'}
        rows.append(row)
        print(i,block['compatible_cover_count'],'covers',len(factors),'/',len(block['pair_resultants']),row['factor_execution']['status'],row['local']['potential_pair_span_interval'],row['local']['realizable_defect_span_interval'],flush=True)
    r.write_new(OUTPUT,{'schema':'rank-jump.fibre-collision-panel.v1','rows':rows,'layer':'solubility',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (INPUT,PROTOCOL,Path(__file__),HERE/'retrospective.py')},
        'boundary':r.read(PROTOCOL)['failure_semantics']})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['run','factor']);p.add_argument('--index',type=int);a=p.parse_args()
    if a.mode=='run':run()
    else:factor_worker(a.index)
