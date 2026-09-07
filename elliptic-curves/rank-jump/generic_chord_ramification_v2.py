#!/usr/bin/env python3
"""Replay the frozen chord experiment with overlapping-power refinement fixed."""
import argparse
from math import gcd
from pathlib import Path
import subprocess
import sys
import generic_chord_ramification as old
import retrospective as r

WORK=r.ROOT/'artifacts/local/rank-jump-generic-chord-ramification-v2'
OUT=r.OUT/'rank_jump_generic_chord_ramification_v2.json'


def insert(basis,n):
    n=abs(n)
    if n<=1:return
    i=0
    while i<len(basis):
        a=basis[i];d=gcd(a,n)
        if d==1:i+=1;continue
        if d==a:
            while n%a==0:n//=a
            if n==1:return
            # Revisit a: a=49,n=343 leaves n=7 and requires splitting a.
            continue
        basis.pop(i)
        insert(basis,d);insert(basis,a//d);insert(basis,n)
        return
    basis.append(n)
    assert len(basis)<=4096


def capture():
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for token in r.read(old.PROTOCOL)['cases']:
        path=WORK/f'{token}.json'
        if not path.exists():
            with (WORK/f'{token}.log').open('x') as log:
                try:
                    proc=subprocess.run([sys.executable,__file__,'worker','--token',token],stdout=log,stderr=log,timeout=60)
                    reason=None if proc.returncode==0 else 'worker failure'
                except subprocess.TimeoutExpired:reason='60-second timeout'
            if reason:r.write_new(path,{'status':'UNKNOWN','token':token,'reason':reason})
        row=r.read(path);rows.append(row);print(token,row['status'],flush=True)
    result={'schema':'rank-jump.generic-chord-ramification.v2','rows':rows,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
                    [Path(__file__),Path(old.__file__),old.PROTOCOL,old.PANEL,old.BOUNDARY,old.SUPPLEMENT]},
        'repair':'Revisit an atom after dividing its powers from the new norm; preserve v1 unchanged.'}
    prior=r.read(old.OUT)
    assert rows==prior['rows']
    result['arithmetic_rows_identical_to_v1']=True
    r.write_new(OUT,result)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker']);p.add_argument('--token');args=p.parse_args()
    old.insert=insert
    if args.mode=='capture':capture()
    else:
        row=old.worker(args.token)
        atoms=list(map(int,row['atoms']))
        assert all(gcd(a,b)==1 for i,a in enumerate(atoms) for b in atoms[i+1:])
        r.write_new(WORK/f'{args.token}.json',row)
