#!/usr/bin/env python3
"""Prove independent outside ramification without factoring any norm."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import three_radical_incidence as source

PROTOCOL=Path(__file__).with_name('THREE_RADICAL_PRIVATE_RAMIFICATION_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_three_radical_private_ramification_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-three-radical-private-ramification-v1'


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
            (Path(__file__),PROTOCOL,source.OUTPUT,Path(r.__file__))}


def worker(token):
    from sage.all import ZZ,prod
    row=next(x for x in r.read(source.OUTPUT)['rows'] if x['token']==token)
    assert row['status']=='PASS'
    norms=[ZZ(x['outside_discriminant_norm']) for x in row['rows']]
    total=prod(norms);records=[]
    for i,n in enumerate(norms):
        others=total//n;private=n
        while True:
            common=private.gcd(others)
            if common==1:break
            private//=common
        witness=not private.is_square()
        assert n%private==0 and private.gcd(others)==1
        records.append({'index':i,'generic_indices':row['rows'][i]['generic_indices'],
            'private_norm_factor':str(private),'private_factor_is_nonsquare':bool(witness)})
    rank=sum(x['private_factor_is_nonsquare'] for x in records)
    return {'token':token,'status':'PASS','bindings':bindings(),'rows':records,
        'candidate_count':len(records),'ramification_rank_lower_bound':rank,
        'candidate_span_unramified_dimension_upper_bound':len(records)-rank,
        'candidate_span_Selmer_dimension':0 if rank==len(records) else 'UNKNOWN',
        'new_Selmer_dimension_after_adding_generic_subgroup':0 if rank==len(records) else 'UNKNOWN'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker']);p.add_argument('--token');args=p.parse_args()
    if args.mode=='worker':r.write_new(WORK/f'{args.token}.json',worker(args.token))
    else:
        WORK.mkdir(parents=True,exist_ok=True);rows=[]
        for case in r.read(source.OUTPUT)['rows']:
            token=case['token'];path=WORK/f'{token}.json'
            if not path.exists():
                with (WORK/f'{token}.log').open('x') as log:
                    try:
                        proc=subprocess.run([sys.executable,__file__,'worker','--token',token],stdout=log,stderr=log,timeout=30)
                        reason=None if proc.returncode==0 else 'worker failure'
                    except subprocess.TimeoutExpired:reason='30-second timeout'
                if reason:r.write_new(path,{'token':token,'status':'UNKNOWN','reason':reason,'bindings':bindings()})
            row=r.read(path);assert row['bindings']==bindings();rows.append(row)
            print(token,row['status'],row.get('ramification_rank_lower_bound'),row.get('candidate_count'),flush=True)
        r.write_new(OUTPUT,{'schema':'rank-jump.three-radical-private-ramification.v1','bindings':bindings(),
            'status':'PASS' if all(x['status']=='PASS' for x in rows) else 'PARTIAL','rows':rows,
            'boundary':'Ranks concern outside ramification of the frozen additive generator spaces, not the complete Selmer or Mordell-Weil group.'})
