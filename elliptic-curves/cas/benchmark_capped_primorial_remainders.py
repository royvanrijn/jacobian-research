#!/usr/bin/env python3
"""Fixed actual-target calibration of capped, full-tree and scalar remainders."""
import argparse
from hashlib import sha256
from math import gcd,prod
from pathlib import Path
import time
import extend_small_conductor_norm_batch as batch
import capped_primorial_remainders as capped
from research_runtime.store import checkpoint

ROOT,ART,cert=batch.ROOT,batch.ART,batch.cert
SOURCE=Path(__file__).resolve()
PROTOCOL=ROOT/'artifacts/local/elliptic-curves/small-conductor-class-target-protected-v1/wave_003/protocol.json'
OUT=ART/'small_conductor_capped_remainders_benchmark_v1.json'


def calculate():
    for modulus in range(40):
        for n in [1,2,3,8,33]:
            values=[1+(i*i+7*i+modulus)%81 for i in range(n)]
            if capped.residues(values,modulus)!=[modulus%v for v in values]:raise ArithmeticError('edge case differs')
    p=cert.read(PROTOCOL);target=p['selection']['targets'][0];v1,v2=target['lattice_basis']
    c=list(map(int,cert.read(batch.forms.OUT)['reduced_binary_cubic_descending']))
    values=[];seen=set()
    for v in range(target['old_vmax']+1,p['vmax']+1):
        for slope in target['root_slopes_scaled']:
            center=slope*v//2**96
            for u in range(center-1,center+2):
                if (u,v) in seen or gcd(u,v)!=1:continue
                seen.add((u,v));m,n=u*v1[0]+v*v2[0],u*v1[1]+v*v2[1]
                if gcd(m,n)!=1:continue
                values.append(abs(((c[0]*m+c[1]*n)*m+c[2]*n*n)*m+c[3]*n*n*n))
    modulus=prod(batch.prior.primes_to(p['smooth_bound']))
    timings={};results=[]
    for name,method in [('full_tree',batch.residues),('capped_tree',capped.residues),('scalar',lambda vs,M:[M%v for v in vs])]:
        start=time.monotonic();results.append(method(values,modulus));timings[name]=time.monotonic()-start
    if not results[0]==results[1]==results[2]:raise ArithmeticError('actual target remainders differ')
    return {'schema':'elliptic-curves.capped-primorial-remainders-benchmark.v1','status':'PASS',
        'sources':{str(q.relative_to(ROOT)):cert.hashed(q) for q in [SOURCE,Path(batch.__file__).resolve(),Path(capped.__file__).resolve(),PROTOCOL,batch.forms.OUT]},
        'target_column':target['column'],'value_count':len(values),'modulus_bits':modulus.bit_length(),
        'edge_case_lists':200,'value_digest':sha256((','.join(map(str,values))).encode()).hexdigest(),
        'remainder_digest':sha256((','.join(map(str,results[0]))).encode()).hexdigest(),
        'exact_agreement':True,'wall_seconds':timings,'measured_full_over_capped_ratio':timings['full_tree']/timings['capped_tree'],
        'proof':'With cap=M+1, each node stores min(true product,cap). Capping commutes with positive multiplication. Every propagated remainder is<=M, so reducing it modulo a product>=cap leaves it unchanged. All leaf remainders are exactly M mod value.',
        'claim_boundary':'Exact arithmetic equivalence, with timing measured on one fixed actual target. No universal speedup or full-search runtime claim.'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--check',action='store_true');a=p.parse_args();r=calculate()
    if a.check:
        old=cert.read(OUT);r['wall_seconds']=old['wall_seconds'];r['measured_full_over_capped_ratio']=old['measured_full_over_capped_ratio']
        if r!=old:raise ArithmeticError('benchmark replay differs')
    else:
        if OUT.exists():raise FileExistsError('preserve benchmark')
        checkpoint(OUT,r)
    print('CAPPED/FULL/SCALAR EXACT AGREEMENT PASS',r['value_count'],r['wall_seconds'],flush=True)
