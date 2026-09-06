#!/usr/bin/env python3
"""Complete bounded residue trees for the26 known non13 scaling-prime candidates."""
import argparse
import json
from math import comb
from pathlib import Path
import certify_compact_r17_candidates as cert
import audit_r17_scaling_prime_support as support
from research_runtime.store import checkpoint

ROOT=support.ROOT
ART=support.ART
D=ROOT/'artifacts/local/elliptic-curves/r17-other-small-prime-scalings-v1'
OUT=ART/'r17_other_small_prime_scalings_v1.json'


def sources():
    paths=[Path(__file__).resolve(),support.INPUT,support.OUT,support.D/'protocol.json',
           Path(support.__file__).resolve(),ROOT/'elliptic-curves/cas/verify_r17_scaling_prime_support.sage',
           support.D/'independent.supervisor.json']
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}


def roster():
    return [{'family':r['family'],'prime':p,'resultant_valuation':e}
            for r in cert.read(support.OUT)['rows'] for p,e in r['trial_prime_factors']
            if p>=5 and p!=13 and e>=4]


def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve fixed small-prime classification')
    s=cert.read(support.D/'independent.supervisor.json')
    if s['outcome']!='completed' or s['returncode']!=0 or len(roster())!=26:raise ArithmeticError('independent resultants and fixed26-pair gate required')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.r17-other-small-prime-scalings.v1',
        'sources':sources(),'pairs':roster(),'maximum_depth':6,'maximum_live_residues_per_chart':4096,
        'maximum_candidates_per_level':200000,'wall_seconds':300,'rss_bytes':2147483648,
        'gate':'The exact homogeneous resultants and independent polynomial checks show26 family/prime pairs other than13 with known prime at least5 and resultant valuation at least4. A primitive parameter with a removable short-model scale must satisfy this necessary resultant condition. Determine actual coefficient congruences at every one of these candidate pairs, without choosing parameter points or invoking factorization.',
        'scope':'Cover all affine residues and the infinity chart, normalizing the coordinate which is a p-adic unit. At depth k exclude a residue only when A or B fails its necessary congruence modulo p^min(k,4) or p^min(k,6). Certify an entire residue ball only by exact coefficients of A(r+p^k*z),B(r+p^k*z) being divisible by p^4,p^6. Lift every remaining residue by all p digits through depth6. Preserve UNKNOWN for any cap or surviving unresolved branch. No parameter scan, point search, larger-prime factorization or automatic deeper lifting; primes2,3 and unresolved large resultant factors remain outside this classification.'})


def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['pairs']!=roster():raise ArithmeticError('fixed small-prime inputs or roster changed')
    return p


def evaluate(coefficients,x,modulus):
    result=0
    for c in reversed(coefficients):result=(result*x+c)%modulus
    return result


def translated(coefficients,residue,step):
    return [sum(coefficients[i]*comb(i,j)*residue**(i-j)*step**j for i in range(j,len(coefficients))) for j in range(len(coefficients))]


def chart(a,b,prime,infinity,p):
    if infinity:a,b=a[::-1],b[::-1]
    current=[0] if infinity else list(range(prime));levels=[];balls=[]
    for depth in range(1,p['maximum_depth']+1):
        if len(current)>p['maximum_candidates_per_level']:
            return {'status':'UNKNOWN_CANDIDATE_CAP','levels':levels,'scale_balls':balls,'cap_depth':depth}
        modulus=prime**depth;ma=prime**min(depth,4);mb=prime**min(depth,6)
        excluded=[];admitted=[];live=[]
        for residue in current:
            if evaluate(a,residue,ma) or evaluate(b,residue,mb):
                excluded.append(residue);continue
            aa=translated(a,residue,modulus);bb=translated(b,residue,modulus)
            if all(c%prime**4==0 for c in aa) and all(c%prime**6==0 for c in bb):
                ball={'depth':depth,'residue':residue,'modulus':modulus,
                      'A_divided_coefficients':[str(c//prime**4) for c in aa],
                      'B_divided_coefficients':[str(c//prime**6) for c in bb]}
                admitted.append(residue);balls.append(ball)
            else:live.append(residue)
        levels.append({'depth':depth,'modulus':modulus,'excluded_residues':excluded,'admitted_residues':admitted,'unresolved_residues':live})
        if not live:return {'status':'COMPLETE_RESIDUE_CLASSIFICATION','levels':levels,'scale_balls':balls}
        if len(live)>p['maximum_live_residues_per_chart']:
            return {'status':'UNKNOWN_LIVE_RESIDUE_CAP','levels':levels,'scale_balls':balls}
        if depth==p['maximum_depth']:
            return {'status':'UNKNOWN_MAXIMUM_DEPTH','levels':levels,'scale_balls':balls}
        current=sorted(residue+modulus*digit for residue in live for digit in range(prime))
    raise AssertionError('unreachable residue loop')


def expected():
    p=protocol();families={r['family']:r for r in cert.read(support.INPUT)['families']};rows=[]
    for pair in p['pairs']:
        f=families[pair['family']];a=list(map(int,f['A_coefficients_low_to_high']));b=list(map(int,f['B_coefficients_low_to_high']))
        charts=[{'chart':'affine',**chart(a,b,pair['prime'],False,p)},
                {'chart':'infinity',**chart(a,b,pair['prime'],True,p)}]
        complete=all(c['status']=='COMPLETE_RESIDUE_CLASSIFICATION' for c in charts)
        count=sum(len(c['scale_balls']) for c in charts)
        status=('CLASSIFIED_SCALE_BALLS' if count else 'NO_REMOVABLE_SCALE') if complete else 'UNKNOWN_INCOMPLETE_CLASSIFICATION'
        rows.append({**pair,'status':status,'scale_balls':count,'charts':charts})
    return {'schema':'elliptic-curves.r17-other-small-prime-scalings-result.v1','status':'PASS_BOUNDED_EXACT_AUDIT',
            'sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'rows':rows,'claim_boundary':p['scope']}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','build','check']);args=parser.parse_args()
    if args.stage=='prepare':prepare()
    else:
        d=expected()
        if args.stage=='check':
            if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('complete residue trees differ')
        else:
            if OUT.exists():raise FileExistsError('preserve bounded small-prime classification')
            checkpoint(OUT,d)
        print('EXACT OTHER-PRIME CLASSIFICATIONS',[(r['family'],r['prime'],r['status'],r['scale_balls']) for r in d['rows']],flush=True)
