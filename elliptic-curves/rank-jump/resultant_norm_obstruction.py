#!/usr/bin/env python3
"""Seven exact family resultants underlying the j-preimage sign obstruction."""
import argparse
from math import isqrt
from pathlib import Path
import retrospective as r
import second_generic_presentation as source

OUTPUT=r.OUT/'rank_jump_resultant_norm_obstruction_v1.json'


def compute():
    from sage.all import QQ,PolynomialRing
    from sage.env import SAGE_VERSION
    R=PolynomialRing(QQ,'s');data=r.read(source.INPUT);rows=[]
    for f in data['families']:
        A=R(f['A']);B=R(f['B']);assert A.degree()==8 and B.degree()==12
        result=A.resultant(B);assert result>0 and not result.is_square()
        n,d=int(result.numerator()),int(result.denominator())
        rows.append({'family':f['family'],'resultant':str(result),'numerator_floor_sqrt':str(isqrt(n)),
            'denominator_floor_sqrt':str(isqrt(d)),'rational_square':False})
    return {'schema':'rank-jump.resultant-norm-obstruction.v1','status':'PASS','families':rows,
        'norm_identity':'For H=b^2*A^3-a^3*B^2 of full degree24 and leading coefficient L: N_{Q[s]/H}(b*A/(a*B))=L^4/Res(A,B).',
        'software':{'sage':SAGE_VERSION},
        'bindings':source.bindings([Path(__file__),source.INPUT]),
        'boundary':'Seven fixed equation resultants, not a parameter search. The norm identity is proved in the note; a nonsquare resultant forbids a square scaling class on an irreducible residual preimage once only rational square-scaling factors have been removed.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS seven nonsquare family resultants')
