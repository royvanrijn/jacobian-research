"""Independently replay one column of the retained rational lift matrix.

FLINT polynomial multiplication and addition, no Groebner membership call.
The large optional input is retained in the ignored local artifact cache.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
from flint import fmpq_mpoly_ctx,fmpq_mpoly

ROOT=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser();parser.add_argument('column',type=int,choices=range(1,7));args=parser.parse_args()
    meta=json.loads((ROOT/'lift_receipt.json').read_text())
    source=ROOT.parents[1]/meta['artifact']
    ctx=fmpq_mpoly_ctx.get(('D','B','c','a','b','z'))
    data=json.loads((ROOT/'critical_one_elimination.json').read_text())
    generators=[fmpq_mpoly(f.replace('**','^'),ctx) for f in data['equations']+['b*z-1']]
    hs=[r.split('=',1)[1] for r in (ROOT/'r7-critical-one.log').read_text().splitlines() if r.startswith('H[')]
    result=ctx.constant(0);rows=[];start=time.monotonic()
    digest=hashlib.sha256()
    with source.open('rb') as handle:
        for line in handle:
            digest.update(line)
            fields=line.decode().split(' ',3)
            if fields[0]!='CERT':continue
            i,j=int(fields[1]),int(fields[2])
            if j!=args.column:continue
            coefficient=fmpq_mpoly(fields[3].strip(),ctx)
            result+=coefficient*generators[i-1]
            rows.append(i)
            print('column',j,'row',i,'seconds',round(time.monotonic()-start,2),flush=True)
    assert digest.hexdigest()==meta['sha256']
    assert sorted(rows)==list(range(1,7))
    assert result==fmpq_mpoly(hs[args.column-1],ctx)
    receipt={'column':args.column,'literal_identity':True,'matrix_sha256':meta['sha256'],
             'generators_sha256':hashlib.sha256((ROOT/'critical_one_elimination.json').read_bytes()).hexdigest(),
             'targets_sha256':hashlib.sha256((ROOT/'r7-critical-one.log').read_bytes()).hexdigest(),
             'method':'exact FLINT rational multivariate multiplication and addition, no Groebner membership decision'}
    (ROOT/('lift-checked-%d.json'%args.column)).write_text(json.dumps(receipt,indent=2)+'\n')
    print('PASS LITERAL LIFT COLUMN',args.column,flush=True)


if __name__=='__main__':main()
