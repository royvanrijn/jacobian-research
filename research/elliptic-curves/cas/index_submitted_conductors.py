#!/usr/bin/env python3
"""Index complete conductor proofs separately from unresolved factorization work."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[2]
WORK=ROOT/'artifacts/local/elliptic-curves/submitted627-630-conductors-v1'
OUT=ROOT/'artifacts/generated-results/elliptic-curves/submitted627_630_conductors_v1'
SUMMARY=OUT/'summary.json'
CURRENT=ROOT/'elliptic-curves/data/submitted_conductors_current.json'


def read(path):
    return json.loads(path.read_text())


def hashed(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def expected():
    rows=[]
    for source in read(WORK/'protocol.json')['roster']:
        identifier=source['icarm_id'];certificate=OUT/f'curve{identifier}.json'
        row={'icarm_id':identifier,'inventory_id':source['inventory_id']}
        if certificate.exists():
            d=read(certificate)
            if d['status']!='PASS' or d['icarm_id']!=identifier or d['ainvs']!=source['ainvs']:
                raise ArithmeticError('conductor certificate identity mismatch')
            for path,digest in d['sources'].items():
                if hashed(ROOT/path)!=digest:
                    raise ArithmeticError('conductor proof source changed')
            row.update(status='EXACT',conductor=d['conductor'],bad_primes=d['bad_primes'],
                       certificate=str(certificate.relative_to(ROOT)),certificate_sha256=hashed(certificate))
        else:
            state=read(WORK/f'ecm_{identifier}/state.json')
            if math.prod(map(int,state['factors']))!=int(source['remaining_cofactor']):
                raise ArithmeticError('partial factorization product mismatch')
            unresolved=[]
            for text in state['factors']:
                n=int(text)
                witness=next((a for a in range(2,30) if math.gcd(a,n)==1 and pow(a,n-1,n)!=1),None)
                if witness is not None:
                    unresolved.append({'integer':text,'decimal_digits':len(text),
                        'fermat_compositeness_witness_base':witness,
                        'residue':str(pow(witness,n-1,n))})
            if not unresolved:
                raise ArithmeticError('unresolved row lacks a concrete composite obstruction')
            row.update(status='UNKNOWN',conductor=None,bad_primes=None,
                       remaining_composites=unresolved,partial_factors=state['factors'],
                       prior_conductor_upper_bound=source['prior_conductor_upper_bound'],
                       factorization_stage_status=state['status'])
        rows.append(row)
    return {'schema':'elliptic-curves.submitted-conductor-index.v1',
            'rows':rows,'exact_count':sum(r['status']=='EXACT' for r in rows),
            'unresolved_count':sum(r['status']=='UNKNOWN' for r in rows),
            'sources':{str(p.relative_to(ROOT)):hashed(p) for p in
                       (Path(__file__).resolve(),WORK/'protocol.json')},
            'claim_boundary':'Complete certified bad-prime lists and exact conductors only '
                'for EXACT rows. UNKNOWN rows supply no submission prime list. Broader '
                '201-curve factorization work is deferred while these four have priority. '
                'No website update or new rank claim.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');args=parser.parse_args()
    d=expected()
    for row in d['rows']:
        if row['status']=='EXACT':
            subprocess.run(['sage','-python',str(ROOT/'elliptic-curves/cas/certify_submitted627_630_conductors.sage'),
                            '--check','--id',str(row['icarm_id'])],check=True)
    current={'schema':'elliptic-curves.current-submitted-conductors.v1',
             'summary':str(SUMMARY.relative_to(ROOT)),
             'summary_sha256':hashlib.sha256((json.dumps(d,sort_keys=True,indent=2)+'\n').encode()).hexdigest()}
    if args.check:
        if read(SUMMARY)!=d or read(CURRENT)!=current:
            raise ArithmeticError('current conductor index differs')
    else:
        for path,value in [(SUMMARY,d),(CURRENT,current)]:
            if path.exists() and read(path)!=value:
                raise FileExistsError('preserve conductor index; make a new version')
            path.write_text(json.dumps(value,sort_keys=True,indent=2)+'\n')
    print('SUBMITTED CONDUCTORS:',d['exact_count'],'EXACT;',d['unresolved_count'],'UNKNOWN')
