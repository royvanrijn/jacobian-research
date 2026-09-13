#!/usr/bin/env sage-python
"""Genus gate for one shared-ordinate cubic tangent from every generic basis pair.

This is a fixed 136-pair identity gate, not a search for rational points.
Only the published generic model and basis are read. A finite polynomial's
odd part gives a lower bound for the characteristic-zero branch degree.
"""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import resource
import time

from sage.all import GF, PolynomialRing, QQ

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'artifacts/generated-results/elkies-k3-r17-correlated-genus-one-v1/input.json'
DEFAULT=ROOT/'artifacts/generated-results/elkies-k3-r17-shared-ordinate-tangents-v1'


def write_new(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x') as f:json.dump(data,f,indent=2,sort_keys=True);f.write('\n')


def odd_part(f):
    f=f.monic();c=f.gcd(f.derivative());w=f//c;out=f.parent()(1);i=1
    while w.degree()>0:
        y=w.gcd(c);z=w//y
        if i%2:out*=z
        w=y;c=c//y;i+=1
    if c.degree()>0:raise ValueError('inseparable remainder: defer this prime')
    return out.monic()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=DEFAULT);args=p.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU,(40,45));resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
    start=time.process_time();source=json.loads(SOURCE.read_text())
    packet={'schema':'r17-shared-ordinate-tangent-input-v1','source_sha256':sha256(SOURCE.read_bytes()).hexdigest(),
            'A':source['A'],'B':source['B'],'basis':source['basis'],'pairs':[list(x) for x in combinations(range(17),2)],
            'primes':[1009,1013,1019],'limits':{'cpu_seconds':40,'address_space_gib':4},
            'selection':'All unordered pairs of the 17 published generic basis sections, one cubic tangent each. Ordinate signs leave the branch squareclass unchanged. No exceptional point or specialization enters the packet.',
            'acceptance':'Retain a candidate for exact construction only if this finite-reduction gate cannot certify genus greater than one.'}
    write_new(args.output/'input.json',packet);rows=[]
    for i,j in packet['pairs']:
        attempts=[];row={'pair':[i,j],'status':'UNKNOWN'}
        for ell in packet['primes']:
            R=PolynomialRing(GF(ell),'t')
            try:
                cv=lambda xs:R([GF(ell)(QQ(x)) for x in xs])
                A,B=cv(packet['A']),cv(packet['B'])
                r,yr=cv(packet['basis'][i]['x']),cv(packet['basis'][i]['y'])
                s,ys=cv(packet['basis'][j]['x']),cv(packet['basis'][j]['y'])
                fr,fs=yr*yr,ys*ys;F,G=3*r*r+A,3*s*s+A
                assert fr==r**3+A*r+B and fs==s**3+A*s+B
                N=fr*fr*G**3-fs*fs*F**3
                M=r*fr*G**2-s*fs*F**2
                if not N or not yr or not ys:
                    attempts.append({'prime':ell,'status':'DEFER_ZERO_REDUCTION'});continue
                U=r*N-3*M*fr*G;V=s*N-3*M*fs*F
                H=U**3+A*U*N*N+B*N**3
                H2=V**3+A*V*N*N+B*N**3
                assert fs*H==fr*H2
                D=H*N
                if not D:
                    attempts.append({'prime':ell,'status':'DEFER_ZERO_BRANCH_POLYNOMIAL'});continue
                if D.degree()>=ell:
                    attempts.append({'prime':ell,'status':'DEFER_DEGREE_AT_LEAST_CHARACTERISTIC'});continue
                odd=odd_part(D);branch=int(odd.degree())+int(odd.degree()%2)
                witness={'prime':ell,'status':'FINITE_BRANCH_WITNESS','raw_degree':int(D.degree()),
                         'odd_part':[int(x) for x in odd.list()],
                         'finite_branch_degree':int(odd.degree()),'total_branch_degree':branch,
                         'genus_lower_bound':max(-1,(branch-2)//2)}
                attempts.append(witness)
                if branch>4:
                    row.update({'status':'EXCLUDED_GENUS_GREATER_THAN_ONE','witness':witness});break
            except (ValueError,ZeroDivisionError) as error:
                attempts.append({'prime':ell,'status':'DEFER_ARITHMETIC_REDUCTION','reason':str(error)})
        row['attempts']=attempts;rows.append(row)
        write_new(args.output/'checkpoints'/('%02d-%02d.json'%(i,j)),row)
    result={'schema':'r17-shared-ordinate-tangent-result-v1','input_sha256':sha256((args.output/'input.json').read_bytes()).hexdigest(),
            'rows':rows,'pair_count':len(rows),'excluded':sum(r['status']=='EXCLUDED_GENUS_GREATER_THAN_ONE' for r in rows),
            'unresolved_pairs':[r['pair'] for r in rows if r['status']=='UNKNOWN'],
            'genus_lower_bound_counts':dict(sorted(Counter(r['witness']['genus_lower_bound'] for r in rows if 'witness' in r).items())),
            'cpu_seconds':time.process_time()-start,
            'scope':'Only one shared-ordinate cubic tangent from each pair in the frozen published 17-section basis. No bound on other generic words, iterated tangents, other identities or arbitrary quadratic-cover ranks.'}
    write_new(args.output/'result.json',result)
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},sort_keys=True),flush=True)


if __name__=='__main__':main()
