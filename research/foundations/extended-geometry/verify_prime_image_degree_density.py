#!/usr/bin/env python3
"""Finite exact arithmetic controls for the separately proved prime-density corollary."""
if not __debug__:
    raise RuntimeError('Assertions must be enabled')
from pathlib import Path
from fractions import Fraction
from itertools import combinations
from math import gcd, isqrt, prod
import argparse, hashlib, json, sys
HERE=Path(__file__).resolve().parent

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def sieve(n):
    a=bytearray(b'\1')*(n+1);a[:2]=b'\0\0'
    for p in range(2,isqrt(n)+1):
        if a[p]:a[p*p:n+1:p]=b'\0'*((n-p*p)//p+1)
    return [p for p in range(2,n+1) if a[p]]
def factor_member(d):
    v=2*d-1
    return any(v%u==0 and v//u>=5 and v//u%4==1 for u in range(5,isqrt(v)+1,4))
def run():
    ledger_path=HERE/'prime_image_degree_density_source_ledger.json'
    ledger=json.loads(ledger_path.read_text())
    for r in ledger['retained_resources']:
        assert sha(HERE/r['retained_path'])==r['sha256'],r['retained_path']
    limit=20000;primes=sieve(2*limit);ps=[p for p in primes if p<=limit]
    family={};a=1
    while (2*a+3)*2+a+2<=limit:
        b=2*a+3;c=a+2
        for n in range(2,(limit-c)//b+1,2):
            d=b*n+c;assert 2*d-1==(2*n+1)*(2*a+3)
            family.setdefault(d,(a,n))
        a+=2
    for d in range(3,limit+1,2):assert (d in family)==factor_member(d),d
    qp={q for q in primes if q%4==1}
    prime_q_family={}
    for q in qp:
        a=(q-3)//2;c=(q+1)//2
        assert a>=1 and a%2==1 and c%2==1 and gcd(c,q)==1
        for p in ps:
            if p>c and p%q==c:
                n=(p-c)//q;assert n>=2 and n%2==0
                assert (2*a+3)*n+a+2==p
                assert (2*p-1)//q==2*n+1
                prime_q_family[p]=(a,n,q)
    for p in ps:
        div=[q for q in qp if q<2*p-1 and (2*p-1)%q==0]
        assert (p in prime_q_family)==bool(div)
        if div:assert p in family
        if p%10==3 and p>3:
            n=(p-3)//5;assert n>=2 and n%2==0 and 5*n+3==p
    assert 3 not in family and 13 in family and 23 in family
    assert 41 in family and 41 not in prime_q_family
    assert 37 not in family
    crt=[];qs=[5,13,17,29]
    for size in range(1,len(qs)+1):
        for subset in combinations(qs,size):
            M=prod(subset);tot=hit=0
            for r in range(1,M):
                if gcd(r,M)!=1:continue
                tot+=1
                hit+=any(r%q==(q+1)//2 for q in subset)
            phi=prod(q-1 for q in subset);miss=prod(q-2 for q in subset)
            assert tot==phi and hit==phi-miss
            density=Fraction(hit,tot)
            assert density==1-prod((Fraction(q-2,q-1) for q in subset),start=Fraction(1))
            crt.append({'q':list(subset),'modulus':M,'unit_residues':tot,'covered':hit,'density':str(density)})
    assert Fraction(1,4)==1-Fraction(5-2,5-1)
    assert Fraction(1,4)!=Fraction(1,5) # Incorrect integer-residue denominator.
    proof=HERE/'PRIME_IMAGE_DEGREE_DENSITY.md'
    return {'schema_version':1,'status':'passed','passed':True,'scope':'Finite exact parameter/divisor/CRT controls only. Classical fixed-PNT and written sequential-limit argument prove densityone; no finite sample proves an asymptotic.','odd_degree_limit':limit,'odd_degree_equivalence_cases':len(range(3,limit+1,2)),'prime_controls':len(ps),'full_family_prime_controls':sum(p in family for p in ps),'prime_q_subfamily_controls':len(prime_q_family),'crt_controls':crt,'negative_controls':['prime3 would require n0','prime41 belongs to full family but not prime-q subfamily','prime37 belongs to neither old parameter set','incorrect denominator q instead of q-1 rejected'],'sources':{p.name:sha(p) for p in [proof,Path(__file__),ledger_path]},'retained_inputs':ledger['retained_resources'],'python':sys.version.split()[0],'independent_replay':False,'external_review':False,'formal_verification':False}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write-receipt',type=Path,required=True);args=ap.parse_args()
    receipt=run();args.write_receipt.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'passed':True,'odd_degree_cases':receipt['odd_degree_equivalence_cases'],'prime_controls':receipt['prime_controls'],'CRT_cases':len(receipt['crt_controls'])}))
