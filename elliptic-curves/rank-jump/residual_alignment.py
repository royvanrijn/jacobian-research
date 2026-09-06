#!/usr/bin/env python3
"""One point-free coordinate swap tests synchronization dependence at 2."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import retrospective as r
import residual_double_covers as original

PROTOCOL=Path(__file__).with_name('RESIDUAL_ALIGNMENT_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_residual_alignment_v1.json'


def dyadic(value):
    value=F(value)
    if not value:return {'zero':True,'square':True}
    n,d=value.numerator,value.denominator;v=0
    while n%2==0:n//=2;v+=1
    while d%2==0:d//=2;v-=1
    unit=n*pow(d,-1,8)%8
    return {'valuation':v,'odd_unit_mod8':unit,'square':v%2==0 and unit==1}


def compute():
    from sage.all import QQ,PolynomialRing
    data=r.read(original.OUTPUT);assert data['status']=='PASS'
    spec=r.read(PROTOCOL);h=[list(map(F,x['quartic_coefficients'])) for x in data['quartics']]
    h=[[c/F(d) for c in coeff] for coeff,d in zip(h,[100,25])]
    assert all(c.denominator==1 for coeff in h for c in coeff)
    U=PolynomialRing(QQ,'x');x=U.gen()
    params=[(s,1) for s in range(256)]+[(1,t) for t in range(0,256,2)]
    rows=[]
    for alignment in spec['alignments']:
        coeff=[h[0],h[1] if alignment=='identity' else list(reversed(h[1]))]
        polys=[sum(QQ(c)*x**(4-i) for i,c in enumerate(q)) for q in coeff]
        res=polys[0].resultant(polys[1]);assert res and all(f.degree()==4 and f.gcd(f.derivative())==1 for f in polys)
        signs=[]
        for e in [1,-1]:
            finite=original.local_test(coeff,2,3,e);witness=None
            for s,t in params:
                values=[e*original.evaluate(q,s,t) for q in coeff];certs=list(map(dyadic,values))
                if all(c['square'] and not c.get('zero') for c in certs):
                    witness={'parameter':[s,t],'values':list(map(str,values)),'square_certificates':certs};break
            if finite['survivor_count']==0:assert witness is None
            signs.append({'epsilon':e,'mod8_test':finite,'Q2_witness':witness,
                'Q2_solubility':'NO' if finite['survivor_count']==0 else 'YES' if witness else 'UNKNOWN'})
        rows.append({'alignment':alignment,'quartic_coefficients':[list(map(str,q)) for q in coeff],
            'resultant':str(res),'product_quotient_genus':3,'simultaneous_curve_genus':5,'signs':signs})
    # All primitive residue pairs prove the simpler parity obstruction too.
    for s in range(8):
        for t in range(8):
            if s%2:assert original.evaluate(h[0],s,t)%8==1
            if t%2:assert original.evaluate(h[1],s,t)%8==1
    return {'schema':'rank-jump.residual-alignment.v1',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),PROTOCOL,original.OUTPUT)},
        'rows':rows,'parity_obstruction_verified':True,
        'boundary':'Exact Q2 results for two chosen alignments. No rational point or full local solubility assertion for either genus-five curve.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();data=compute()
    if a.mode=='check':assert r.read(OUTPUT)==data;print('PASS residual alignment')
    else:r.write_new(OUTPUT,data)
    for row in data['rows']:print(row['alignment'],[(x['epsilon'],x['Q2_solubility'],x['Q2_witness']) for x in row['signs']])
