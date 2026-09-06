#!/usr/bin/env python3
"""Independent Sage invariant and exact valuation replay of the coefficient gate."""
import argparse
from pathlib import Path
import retrospective as r
import horizontal_norm_gate as gate

OUTPUT=r.OUT/'rank_jump_horizontal_norm_gate_verification_v1.json'


def compute():
    from sage.all import QQ,ZZ,EllipticCurve
    data=r.read(gate.OUTPUT);models={x['case_index']:x['model'] for x in r.read(gate.INPUT)['rows']}
    for p in data['primes']:assert ZZ(p).is_prime(proof=True) and ZZ(-3).kronecker(p)==-1
    rows=[]
    for row in data['rows']:
        E=EllipticCurve(QQ,list(map(QQ,models[row['case_index']])))
        A=-E.c4()/48
        assert A==QQ(row['short_A'])
        observed=[{'prime':p,'valuation_minus_A':int((-A).valuation(p))} for p in data['primes'] if (-A).valuation(p)%2]
        assert observed==row['odd_inert_valuation_obstructions']
        assert (row['status']=='EXCLUDED')==bool(observed)
        rows.append({'case_index':row['case_index'],'c4':str(E.c4()),'obstructions':observed,'status':'PASS'})
    # Separate elementary local check at the actual obstruction prime.
    p=23
    zero_pairs=[(x,y) for x in range(p) for y in range(p) if (x*x+3*y*y)%p==0]
    assert zero_pairs==[(0,0)]
    return {'schema':'rank-jump.horizontal-norm-gate-verification.v1',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),gate.INPUT,gate.OUTPUT)},
        'rows':rows,'mod_23_zero_pairs':[list(x) for x in zero_pairs],
        'local_conclusion':'Any nonzero rational value ell^2+3*z^2 has even 23-adic valuation.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();data=compute()
    if args.mode=='check':assert r.read(OUTPUT)==data;print('PASS independent c4 and norm-gate verification')
    else:r.write_new(OUTPUT,data)
