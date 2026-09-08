#!/usr/bin/env python3
"""Exact homogeneous resultants bound possible prime-dependent model scalings."""
import argparse
import json
from pathlib import Path
import certify_compact_r17_candidates as cert
from mod2_reduction_independence import _is_prime
from research_runtime.store import checkpoint

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
INPUT=ART/'compact_six_r17_atlas_v1.json'
OUT=ART/'r17_scaling_prime_support_v1.json'
D=ROOT/'artifacts/local/elliptic-curves/r17-scaling-prime-support-v1'


def sources():
    paths=[Path(__file__).resolve(),INPUT,ROOT/'elliptic-curves/cas/mod2_reduction_independence.py']
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}


def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve fixed prime-support audit')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.r17-scaling-prime-support.v1',
        'sources':sources(),'families':6,'A_degree':8,'B_degree':12,
        'maximum_trial_prime':997,'wall_seconds':120,'rss_bytes':1073741824,
        'scope':'For the six compact binary coefficient pairs compute the exact20-by20 homogeneous Sylvester determinant. Trial-divide each nonzero determinant by primes through997, preserving the full residual cofactor as UNKNOWN if it is not1. No integer-factorization escalation, parameter scan, point search or score update. At a primitive parameter, a prime dividing both specialized coefficients must divide the homogeneous resultant. This bounds candidate model-scaling primes but does not assert that any such scaling occurs, nor local minimality at2or3, smaller coefficients, new curves or rank gains.'})


def sylvester(a,b):
    if len(a)!=9 or len(b)!=13:raise ArithmeticError('fixed binary degrees8 and12 required')
    rows=[]
    for k in range(12):rows.append([0]*k+a[::-1]+[0]*(11-k))
    for k in range(8):rows.append([0]*k+b[::-1]+[0]*(7-k))
    return rows


def determinant(rows):
    a=[r[:] for r in rows];sign=1;previous=1;n=len(a)
    for k in range(n-1):
        pivot=next((i for i in range(k,n) if a[i][k]),None)
        if pivot is None:return 0
        if pivot!=k:a[k],a[pivot]=a[pivot],a[k];sign=-sign
        current=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                value=a[i][j]*current-a[i][k]*a[k][j]
                if value%previous:raise ArithmeticError('nonexact Bareiss division')
                a[i][j]=value//previous
            a[i][k]=0
        previous=current
    return sign*a[-1][-1]


def expected():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('frozen resultant inputs changed')
    rows=[]
    for f in cert.read(INPUT)['families']:
        a=list(map(int,f['A_coefficients_low_to_high']));b=list(map(int,f['B_coefficients_low_to_high']))
        result=determinant(sylvester(a,b));remaining=abs(result);factors=[]
        if remaining:
            for prime in range(2,p['maximum_trial_prime']+1):
                if not _is_prime(prime):continue
                exponent=0
                while remaining%prime==0:remaining//=prime;exponent+=1
                if exponent:factors.append([prime,exponent])
        product=remaining
        for prime,exponent in factors:product*=prime**exponent
        if product!=abs(result):raise ArithmeticError('exact factor reconstruction differs')
        rows.append({'family':f['family'],'homogeneous_resultant':str(result),
                     'trial_prime_factors':factors,'remaining_cofactor':str(remaining),
                     'complete_prime_support':remaining==1,
                     'support_primes_at_least5':[q for q,e in factors if q>=5]})
    if len(rows)!=6:raise ArithmeticError('all six binary pairs required')
    return {'schema':'elliptic-curves.r17-scaling-prime-support-result.v1','status':'PASS_EXACT_BOUNDED_AUDIT',
            'sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'rows':rows,'claim_boundary':p['scope']}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','build','check']);args=parser.parse_args()
    if args.stage=='prepare':prepare()
    else:
        data=expected()
        if args.stage=='check':
            if cert.read(OUT)!=json.loads(json.dumps(data)):raise ArithmeticError('exact resultant support differs')
        else:
            if OUT.exists():raise FileExistsError('preserve prime-support certificate')
            checkpoint(OUT,data)
        print('EXACT RESULTANT PRIME SUPPORT',[(r['family'],r['trial_prime_factors'],r['remaining_cofactor']) for r in data['rows']],flush=True)
