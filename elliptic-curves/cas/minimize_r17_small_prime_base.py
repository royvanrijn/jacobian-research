#!/usr/bin/env python3
"""Finite exact small-prime base-neighbour descent for integral R17 models.

Accept only strict decreases of max(|A_i|^3,|B_j|^2). This is a coefficient
height heuristic, not a proof of global minimality or a rank predictor.
"""
import argparse
from math import comb
from pathlib import Path
import certify_compact_r17_candidates as cert
from mod2_reduction_independence import _is_prime
from research_runtime.store import checkpoint, digest

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / 'artifacts/generated-results/elliptic-curves/r17_constant_scaling_audit_v1.json'
DIRECTORY = ROOT / 'artifacts/local/elliptic-curves/r17-small-prime-base-v1'


def transform(f, weight, m):
    a,b,c,d = m; result = [0]*(weight+1)
    for i,v in enumerate(f):
        if not v: continue
        for j in range(i+1):
            left = v*comb(i,j)*a**j*b**(i-j)
            if not left: continue
            for k in range(weight-i+1):
                result[j+k] += left*comb(weight-i,k)*c**k*d**(weight-i-k)
    return result


def multiply(m,n):
    a,b,c,d=m; e,f,g,h=n
    return [a*e+b*g,a*f+b*h,c*e+d*g,c*f+d*h]


def score(A,B):
    return max(max(abs(x) for x in A)**3,max(abs(x) for x in B)**2)


def valuation(n,p):
    if n==0: return 10**9
    e=0
    while n%p==0:n//=p;e+=1
    return e


def roots(A,B,p):
    aa,bb=[[x%p for x in f] for f in (A,B)]
    def vanishes(f,r,count):
        for j in range(count):
            if sum(f[i]*comb(i,j)*pow(r,i-j,p) for i in range(j,len(f)))%p:return False
        return True
    result=[r for r in range(p) if vanishes(aa,r,4) and vanishes(bb,r,6)]
    if all(x==0 for x in aa[-4:]) and all(x==0 for x in bb[-6:]):result.append('infinity')
    return result


def source_hashes():
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),INPUT)}


def prepare(directory):
    if (directory/'protocol.json').exists():raise FileExistsError('protocol already frozen')
    checkpoint(directory/'protocol.json',{'schema':'elliptic-curves.r17-small-prime-base.v1',
        'sources':source_hashes(),'families':['103b2','11952','074d9','07ca9','08234','08f72'],
        'prime_bound':997,'maximum_steps':128,'worker_wall_seconds':120,'worker_rss_bytes':1073741824,
        'maximum_concurrent_workers':2,
        'rule':'At every step test determinant-p neighbours at roots with multiplicities at least 4 in A and 6 in B, including infinity. Remove maximal common p scale of weights 4/6. Choose the smallest strictly decreasing exact max(abs(A)^3,abs(B)^2); ties by prime then root order.',
        'gate':'Begin with 103b2. Run the other five only if its final largest coefficient bit size decreases by at least 25 percent.',
        'scope':'Coordinate changes only; no parameter or point search, no global minimality assertion, and no generic rank increase.'})


def run(directory,family):
    protocol=cert.read(directory/'protocol.json')
    if protocol['sources']!=source_hashes() or family not in protocol['families']:raise ArithmeticError('protocol mismatch')
    output=directory/(family+'.json')
    if output.exists():raise FileExistsError('preserve previous attempt')
    row=next(r for r in cert.read(INPUT)['rows'] if r['family']==family)
    A,B=[[cert.F(x) for x in row[k]] for k in ('A_coefficients_low_to_high','B_coefficients_low_to_high')]
    if any(x.denominator!=1 for x in A+B):raise ArithmeticError('input not integral')
    A,B=list(map(int,A)),list(map(int,B));A += [0]*(9-len(A));B += [0]*(13-len(B))
    initial_A,initial_B=A[:],B[:]
    primes=[p for p in range(2,998) if _is_prime(p) and roots(A,B,p)]
    # Outside these primes, determinant-p changes and p rescaling are
    # invertible, preserving the absence of high-multiplicity common roots.
    M=[1,0,0,1];u=1;steps=[]
    result={'status':'RUNNING','family':family,'protocol_sha256':cert.hashed(directory/'protocol.json'),
        'source':row['source'],'source_sha256':row['source_sha256'],'input_sha256':cert.hashed(INPUT),
        'initial_primes':primes,'initial_constant_scale_u':row['scale_u'],'initial_bits':row['after_bits'],'steps':steps}
    for iteration in range(protocol['maximum_steps']):
        current=score(A,B);best=None
        for p in primes:
            for r in roots(A,B,p):
                m=[1,0,0,p] if r=='infinity' else [p,r,0,1]
                aa,bb=transform(A,8,m),transform(B,12,m)
                e=min(min(valuation(x,p) for x in aa)//4,min(valuation(x,p) for x in bb)//6)
                if e<1:continue
                aa=[x//p**(4*e) for x in aa];bb=[x//p**(6*e) for x in bb]
                candidate=score(aa,bb)
                if candidate<current and (best is None or candidate<best[0]):best=(candidate,p,r,e,m,aa,bb)
        if best is None:result['status']='COMPLETE_BOUNDED_LOCAL_DESCENT';break
        value,p,r,e,m,A,B=best;M=multiply(M,m);u*=p**e
        steps.append({'prime':p,'root':r,'scale_exponent':e,'before_score_bits':current.bit_length(),
            'after_score_bits':value.bit_length(),'model_digest':digest([A,B])})
        result.update(base_matrix_a_b_c_d=M,additional_scale_u=str(u),
            A_coefficients_low_to_high=list(map(str,A)),B_coefficients_low_to_high=list(map(str,B)))
        checkpoint(output,result)
        print('LOCAL BASE',family,len(steps),'p',p,'root',r,'bits',max(abs(x).bit_length() for x in A+B),flush=True)
    else:result['status']='STEP_LIMIT_REACHED'
    a,b,c,d=M;det=a*d-b*c
    if [x*u**4 for x in transform(A,8,[d,-b,-c,a])]!=[x*det**8 for x in initial_A]:raise ArithmeticError('inverse A identity failed')
    if [x*u**6 for x in transform(B,12,[d,-b,-c,a])]!=[x*det**12 for x in initial_B]:raise ArithmeticError('inverse B identity failed')
    result.update(base_matrix_a_b_c_d=M,additional_scale_u=str(u),
        total_scale_from_literal_source=str(cert.F(row['scale_u'])*u),
        A_coefficients_low_to_high=list(map(str,A)),B_coefficients_low_to_high=list(map(str,B)),
        after_bits=max(abs(x).bit_length() for x in A+B),inverse_coefficient_identities_verified=True)
    checkpoint(output,result);print('COMPLETE',family,row['after_bits'],'->',result['after_bits'],flush=True)


def check(path):
    row=cert.read(path);source=ROOT/row['source']
    if cert.hashed(source)!=row['source_sha256']:raise ArithmeticError('source changed')
    data=cert.read(source);model=data['representative' if row['family']=='074d9' else 'weierstrass_model']
    u=cert.F(row['total_scale_from_literal_source']);M=row['base_matrix_a_b_c_d']
    if M[0]*M[3]-M[1]*M[2]==0 or not u:raise ArithmeticError('singular coordinate change')
    for key,weight,power in [('A_coefficients_low_to_high',8,4),('B_coefficients_low_to_high',12,6)]:
        before=list(map(cert.F,model[key]));after=list(map(cert.F,row[key]))
        if transform(before,weight,M)!=[x*u**power for x in after]:raise ArithmeticError('literal coefficient identity failed')
    print('REPLAYED BASE IDENTITY',row['family'],row['after_bits'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','run','check'])
    p.add_argument('--directory',type=Path,default=DIRECTORY);p.add_argument('--family',default='103b2');p.add_argument('--input',type=Path)
    a=p.parse_args()
    if a.stage=='prepare':prepare(a.directory)
    elif a.stage=='run':run(a.directory,a.family)
    else:check(a.input)
