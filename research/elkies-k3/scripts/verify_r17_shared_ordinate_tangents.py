#!/usr/bin/env python3
"""Independent arithmetic replay of the frozen shared-ordinate tangent genus gate.

The characteristic-zero lower bound also uses the written odd-part reduction
argument. No Sage, producer factorization or proposed branch count is trusted.
"""
import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import importlib.util
from itertools import combinations
import json
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'artifacts/generated-results/elkies-k3-r17-correlated-genus-one-v1/input.json'
DEFAULT=ROOT/'artifacts/generated-results/elkies-k3-r17-shared-ordinate-tangents-v1'
SPEC=importlib.util.spec_from_file_location('tangent_finite_replay',ROOT/'elkies-k3/scripts/verify_r17_mestre_shared_twist.py')
V=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(V)
require=V.require


def ff_div(a,b,p):
    require(b,'zero divisor');a=list(a);q=[0]*max(0,len(a)-len(b)+1)
    inv=pow(b[-1],-1,p)
    while a and len(a)>=len(b):
        shift=len(a)-len(b);c=a[-1]*inv%p;q[shift]=c
        for j,x in enumerate(b):a[shift+j]=(a[shift+j]-c*x)%p
        a=V.trim(a)
    require(not a,'inexact finite polynomial division')
    return V.trim(q)


def odd_part(f,p):
    require(f and len(f)-1<p,'zero or too-high-degree branch polynomial')
    f=[x*pow(f[-1],-1,p)%p for x in f]
    c=V.ff_gcd(f,V.trim([i*x%p for i,x in enumerate(f)][1:]),p)
    w=ff_div(f,c,p);out=[1];i=1
    while len(w)>1:
        y=V.ff_gcd(w,c,p);z=ff_div(w,y,p)
        if i%2:out=V.ff_mul(out,z,p)
        w=y;c=ff_div(c,y,p);i+=1
    require(len(c)==1,'unhandled inseparable remainder')
    return [x*pow(out[-1],-1,p)%p for x in out]


def finite_witness(packet,pair,p):
    # Unlike an irreducibility witness, this bound permits degree drop.
    require(p>=2 and all(p%d for d in range(2,int(p**0.5)+1)),'nonprime reduction')
    def cv(xs):
        qs=[Fraction(x) for x in xs]
        require(all(x.denominator%p for x in qs),'nonintegral generic coefficient')
        return V.trim([x.numerator*pow(x.denominator,-1,p)%p for x in qs])
    mul=lambda a,b:V.ff_mul(a,b,p)
    sub=lambda a,b:V.ff_sub(a,b,p)
    scale=lambda a,k:V.trim([k*x%p for x in a])
    add=lambda a,b:sub(a,scale(b,-1))
    square=lambda a:mul(a,a)
    cube=lambda a:mul(square(a),a)
    A,B=cv(packet['A']),cv(packet['B']);i,j=pair
    r,yr=cv(packet['basis'][i]['x']),cv(packet['basis'][i]['y'])
    s,ys=cv(packet['basis'][j]['x']),cv(packet['basis'][j]['y'])
    require(yr and ys,'zero ordinate reduction')
    fr,fs=square(yr),square(ys)
    F,G=add(scale(square(r),3),A),add(scale(square(s),3),A)
    N=sub(mul(square(fr),cube(G)),mul(square(fs),cube(F)))
    M=sub(mul(mul(r,fr),square(G)),mul(mul(s,fs),square(F)))
    require(N,'zero tangent denominator reduction')
    U=sub(mul(r,N),scale(mul(mul(M,fr),G),3))
    W=sub(mul(s,N),scale(mul(mul(M,fs),F),3))
    evaluate=lambda x:add(add(cube(x),mul(mul(A,x),square(N))),mul(B,cube(N)))
    H,H2=evaluate(U),evaluate(W)
    require(mul(fs,H)==mul(fr,H2),'finite shared-ordinate identity')
    D=mul(H,N);odd=odd_part(D,p);finite=len(odd)-1;branch=finite+finite%2
    return {'raw_degree':len(D)-1,'odd_part':odd,'finite_branch_degree':finite,
            'total_branch_degree':branch,'genus_lower_bound':(branch-2)//2}


def verify(path):
    packet=json.loads((path/'input.json').read_text());result=json.loads((path/'result.json').read_text())
    require(packet['schema']=='r17-shared-ordinate-tangent-input-v1','input schema')
    source=json.loads(SOURCE.read_text())
    require(packet['source_sha256']==sha256(SOURCE.read_bytes()).hexdigest(),'generic source binding')
    require(result['input_sha256']==sha256((path/'input.json').read_bytes()).hexdigest(),'frozen input binding')
    require(all(packet[k]==source[k] for k in ('A','B','basis')),'generic-only source projection')
    V.V.source_projection(source)
    exact=V.V;A,B=exact.poly(packet['A']),exact.poly(packet['B'])
    for row in packet['basis']:
        x,y=exact.poly(row['x']),exact.poly(row['y'])
        require(exact.power(y,2)==exact.add(exact.add(exact.power(x,3),exact.mul(A,x)),B),'exact inherited section identity')
    pairs=[list(x) for x in combinations(range(17),2)]
    require(packet['pairs']==pairs and packet['primes']==[1009,1013,1019],'frozen basis-pair and prime panel')
    require([r['pair'] for r in result['rows']]==pairs,'complete unique basis-pair coverage')
    counts=Counter()
    for row in result['rows']:
        require(row['status']=='EXCLUDED_GENUS_GREATER_THAN_ONE','unresolved pair')
        witness=row['witness'];p=witness['prime'];require(p in packet['primes'],'unfrozen witness prime')
        calculated=finite_witness(packet,row['pair'],p)
        require(all(witness[k]==v for k,v in calculated.items()),'altered finite branch witness')
        require(calculated['total_branch_degree']>4,'no genus obstruction')
        counts[calculated['genus_lower_bound']]+=1
    return {'status':'PASS_136_SHARED_ORDINATE_TANGENT_GENUS_EXCLUSIONS','excluded_pairs':len(pairs),
            'genus_lower_bound_counts':dict(sorted(counts.items())),'minimum_genus_lower_bound':min(counts),
            'no_positive_candidate':True,'independence_not_assessed':True,'goal_complete':False,
            'scope':'One tangent from each unordered pair of the 17 published generic basis sections; both ordinate signs give the same branch class. Other words, iterates and identities remain outside this fixed gate.'}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input-dir',type=Path,default=DEFAULT);p.add_argument('--output',type=Path);args=p.parse_args()
    start=time.monotonic();result=verify(args.input_dir);result['elapsed_seconds']=round(time.monotonic()-start,3)
    if args.output:
        with args.output.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps(result,sort_keys=True))


if __name__=='__main__':main()
