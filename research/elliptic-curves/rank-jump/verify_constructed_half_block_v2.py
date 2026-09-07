#!/usr/bin/env python3
"""Add proved primality and complete prime-ideal replay to the half-block certificate."""
import argparse
from pathlib import Path
from math import prod
import retrospective as r
import verify_constructed_half_block as previous
OUTPUT=r.OUT/'rank_jump_constructed_half_block_verification_v2.json'

def compute():
    from sage.all import QQ,ZZ,PolynomialRing,pari
    result=previous.compute();assert result==r.read(previous.OUTPUT)
    data=r.read(previous.ARTIN);ref=r.read(previous.half.REFERENCE)
    R=PolynomialRing(QQ,'z');nf=pari.nfinit([pari(R(ref['cubic_ascending'])),ref['S_finite']])
    primes=sorted({row['p'] for row in data['prime_artin_certificates']})
    assert all(ZZ(p).is_prime(proof=True) for p in primes)
    lookup={};records={}
    for row in data['prime_artin_certificates']:
        P=next(P for P in pari.idealprimedec(nf,row['p']) if str(pari.idealhnf(nf,P))==row['hnf'])
        assert int(P[2])==row['e'] and int(P[3])==row['f']
        assert pari.idealval(nf,pari(row['uniformizer_GP']),P)==1
        lookup[row['hnf']]=P;records[row['hnf']]=row
    for col in data['columns']:
        I=pari(col['reduced_half_ideal_hnf']);norm=QQ(pari.idealnorm(nf,I));product=pari.idealhnf(nf,1)
        assert str(norm)==col['ideal_norm']
        for q in col['prime_factors']:
            P=lookup[q['prime_hnf']];e=q['valuation'];assert e>=0 and pari.idealval(nf,I,P)==e
            product=pari.idealmul(nf,product,pari.idealpow(nf,P,e))
        assert product==I
    result.update(schema='rank-jump.constructed-half-block-verification.v2',
        proved_rational_primes=primes,verified_prime_ideal_records=len(lookup),verified_complete_ideal_factorizations=len(data['columns']))
    result['bindings'].update({str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in [Path(__file__),previous.OUTPUT]})
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();result=compute()
    if a.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS elementary factor and proved prime-ideal factorizations')
