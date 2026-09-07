#!/usr/bin/env python3
"""Retain the failed-prime panel and rerun one fixed larger panel."""
import argparse
from pathlib import Path
import subprocess
import retrospective as r
import direct_relation_local_gate as d

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'DIRECT_RELATION_GOOD_PRIMES_PROTOCOL.json'
OUTPUT=r.OUT/'rank_jump_direct_relation_good_primes_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-direct-relation-good-primes-v1'


def compute():
    d.PROTOCOL=PROTOCOL
    data=d.compute()
    data['bindings'][str(Path(__file__).resolve().relative_to(r.ROOT))]=r.digest(Path(__file__).read_bytes())
    return data


def run():
    WORK.mkdir(parents=True,exist_ok=True);log=WORK/'worker.log';execution=WORK/'execution.json'
    if not log.exists():
        with log.open('x') as out:
            try:
                p=subprocess.run(['/home/royvanrijn/.local/bin/sage','-python',str(Path(__file__).resolve()),'worker'],stdout=out,stderr=out,timeout=60)
                status={'status':'COMPLETE' if p.returncode==0 else 'FAILED','returncode':p.returncode}
            except subprocess.TimeoutExpired:status={'status':'TIMEOUT'}
        r.write_new(execution,status)
    print(r.read(execution));print(log.read_text())


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['run','worker','check']);a=p.parse_args()
    if a.mode=='run':run()
    else:
        data=compute()
        if a.mode=='worker':r.write_new(OUTPUT,data)
        else:assert r.read(OUTPUT)==data
        for row in data['rows']:print(row['prime'],row['status'],[(x['index'],x['Fp_point_count']) for x in row.get('relations',[])],flush=True)
