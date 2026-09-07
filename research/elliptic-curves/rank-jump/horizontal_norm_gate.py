#!/usr/bin/env python3
"""Coefficient-only exact local obstruction to rational horizontal pairs."""
import argparse
from pathlib import Path
import retrospective as r

PROTOCOL=Path(__file__).with_name('HORIZONTAL_NORM_GATE_PROTOCOL.json')
INPUT=r.OUT/'rank_jump_horizontal_norm_gate_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_horizontal_norm_gate_v1.json'


def valuation(x,p):
    x=r.F(x)
    if not x:raise ValueError('zero valuation is not finite')
    n,d=abs(x.numerator),x.denominator;e=0
    while n%p==0:n//=p;e+=1
    while d%p==0:d//=p;e-=1
    return e


def norm_gate(A,primes):
    a=-r.F(A)
    if a<0:return {'status':'EXCLUDED','obstruction':'real','minus_A':str(a)}
    if a==0:return {'status':'EXCLUDED','obstruction':'zero positive-definite norm allows only ell=z=0','minus_A':'0'}
    records=[]
    for p in primes:
        assert p>=5 and p%3==2 and pow(p-3,(p-1)//2,p)==p-1
        e=valuation(a,p)
        if e:records.append({'prime':p,'valuation_minus_A':e})
    odd=[row for row in records if row['valuation_minus_A']%2]
    return {'status':'EXCLUDED' if odd else 'UNKNOWN','minus_A':str(a),
            'nonzero_dictionary_valuations':records,'odd_inert_valuation_obstructions':odd}


def export():
    import bad_prime_support as bad
    rows=[{'case_index':i,'model':bad.cases()[i]['model']} for i in r.read(PROTOCOL)['cases']]
    r.write_new(INPUT,{'schema':'rank-jump.horizontal-norm-gate-inputs.v1','rows':rows})


def calculate():
    spec=r.read(PROTOCOL)['prime_dictionary']
    primes=[p for p in r.primes(spec['upper_bound']) if p>=spec['lower_bound'] and p%3==spec['congruence_mod_3']]
    rows=[]
    for row in r.read(INPUT)['rows']:
        assert set(row)=={'case_index','model'}
        short,_=r.short(row['model'],[]);A,B=map(r.F,short[3:])
        result=norm_gate(A,primes)
        rows.append({'case_index':row['case_index'],'short_A':str(A),'short_B':str(B),**result})
    controls=[]
    for c,n in [(1,5),(3,7)]:
        m=r.F(n*n-c*c-2,2);short,_=r.short([0,m,0,-m-3,c*c],[]);A=r.F(short[3])
        z=r.F(1,2)+m/3;ell=r.F(3,2)
        assert ell*ell+3*z*z==-A
        result=norm_gate(A,primes);assert result['status']=='UNKNOWN'
        controls.append({'c':c,'n':n,'short_A':str(A),'z':str(z),'ell':str(ell),'status':'EXACT_NORM_WITNESS',
            'dictionary_status':result['status']})
    return {'schema':'rank-jump.horizontal-norm-gate.v1',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),PROTOCOL,INPUT)},
        'primes':primes,'rows':rows,'controls':controls,
        'boundary':'An exclusion rules out every distinct rational horizontal pair in any Weierstrass model when completed ordinates are used. It does not exclude unequal-ordinate points sharing a quadratic-cover squareclass or any high-rank curve.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['export','build','check']);args=p.parse_args()
    if args.mode=='export':export()
    else:
        data=calculate()
        if args.mode=='check':assert r.read(OUTPUT)==data;print('PASS coefficient-only horizontal norm obstruction')
        else:r.write_new(OUTPUT,data)
        for row in data['rows']:print(row['case_index'],row['status'],row.get('odd_inert_valuation_obstructions'))
