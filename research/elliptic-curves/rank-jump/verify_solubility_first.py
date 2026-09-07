#!/usr/bin/env python3
"""Independent rational-square replay and polynomial lift-identity verification."""
from fractions import Fraction as F
import gzip
import json
from math import isqrt
from pathlib import Path
import retrospective as r
from cover_experiment import evaluate, mul, sub, trim
import solubility_first as s

OUTPUT = r.OUT/'rank_jump_solubility_first_verification_v1.json'


def add(a,b): return sub(a,[-x for x in b])
def scale(a,c): return trim([c*x for x in a])


def verify():
    result = r.read(s.OUTPUT); inp = r.read(s.INPUT)
    for path,sha in result['bindings'].items(): assert r.digest((r.ROOT/path).read_bytes())==sha
    geometry = json.loads(gzip.decompress(s.GEOMETRY.read_bytes()))
    assert len({c['label'] for c in geometry})==39119
    qs = {c['label']:[F(x,c['denominator_scale']**2) for x in c['integer_quadratic']] for c in geometry}
    # Two polynomial identities certify each selected map before specialization:
    # y0^2+q*y1^2 = x0^3+3q*x0*x1^2+A*x0+B,
    # 2*y0*y1 = 3*x0^2*x1+q*x1^3+A*x1.
    A,B = list(map(F,inp['A'])),list(map(F,inp['B']))
    for label,c in inp['split_lift_maps'].items():
        x0,x1,y0,y1 = (list(map(F,c[k+'_coefficients'])) for k in ('x0','x1','y0','y1'))
        q=qs[label]
        lhs=add(mul(y0,y0),mul(q,mul(y1,y1)))
        rhs=add(add(mul(mul(x0,x0),x0),scale(mul(q,mul(x0,mul(x1,x1))),3)),add(mul(A,x0),B))
        assert sub(lhs,rhs)==[0],label
        lhs=scale(mul(y0,y1),2)
        rhs=add(add(scale(mul(mul(x0,x0),x1),3),mul(q,mul(mul(x1,x1),x1))),mul(A,x1))
        assert sub(lhs,rhs)==[0],label
    tested=0
    for row in result['rows']:
        hits={};zeros=[];t=F(row['published_parameter'])
        assert t==-(F(row['compact_parameter'])+50)/26
        for label,q in qs.items():
            v=evaluate(q,t);tested+=1
            if not v:zeros.append(label)
            elif v>0:
                a,b=isqrt(v.numerator),isqrt(v.denominator)
                if a*a==v.numerator and b*b==v.denominator:hits[label]=str(F(a,b))
        assert hits=={h['label']:h['square_root'] for h in row['nonzero_square_hits']}
        assert zeros==row['branch_degeneracies']
    assert tested==1251808
    return {'schema':'rank-jump.solubility-first-verification.v1','status':'PASS',
            'result_sha256':r.digest(s.OUTPUT.read_bytes()), 'script_sha256':r.digest(Path(__file__).read_bytes()),
            'rational_square_evaluations':tested,'generic_lift_maps_verified':len(inp['split_lift_maps']),
            'boundary':'Independent rational-fraction square evaluation and generic lift identities. Finite Kummer rank is replayed separately by the producer check.'}


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args()
    data=verify()
    if a.mode=='build':r.write_new(OUTPUT,data)
    else:assert r.read(OUTPUT)==data
    print(data)
