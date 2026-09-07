#!/usr/bin/env python3
"""Bounded metric relation proposal followed by exact group identity checking."""
from sage.all import EllipticCurve, QQ, RealField, matrix, vector
from importlib.machinery import SourceFileLoader
from fractions import Fraction
from pathlib import Path
import json
from hashlib import sha256
import argparse

ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
PUBLIC=ART/'inventory188_public28_reproduction_v1.json'
CLOUD=ART/'inventory188_fixed49_point_control_mod2_v1.json'
OUT=ART/'inventory188_recovered_direction_relation_v1.json'


def calculate():
    public=json.loads(PUBLIC.read_text());cloud=json.loads(CLOUD.read_text())
    assert cloud['curve']==public['curve'] and cloud['rank_lower_bound']==28
    assert cloud['independent_points'][:27]==public['points'][:27]
    E=EllipticCurve(QQ,public['curve'])
    basis=[E(list(map(QQ,p))) for p in public['points'][:27]]
    basis.append(E(list(map(QQ,public['transported_public_points'][26]))))
    P=E(list(map(QQ,cloud['independent_points'][27])))
    geometry=SourceFileLoader('relation_geometry',str(ROOT/'elliptic-curves/cas/prospective_half_lattice_v2.sage')).load_module()
    gram,_=geometry.canonical_height_gram(tuple(map(Fraction,public['curve'])),
        [tuple(map(Fraction,map(str,p.xy()))) for p in basis+[P]])
    real=RealField(384);G=matrix(real,[[str(x) for x in row] for row in gram])
    approx=G[:28,:28].solve_right(G.column(28)[:28])
    word=[QQ(Fraction(str(x)).limit_denominator(64)) for x in approx]
    denominator=QQ(1).denominator()
    from sage.all import lcm
    denominator=lcm(c.denominator() for c in word)
    integers=[int(c*denominator) for c in word]
    exact=sum((c*B for c,B in zip(integers,basis)),E(0))==denominator*P
    paths=[PUBLIC,CLOUD,Path(__file__).resolve(),ROOT/'elliptic-curves/cas/prospective_half_lattice_v2.sage']
    return {'schema':'elliptic-curves.inventory188-recovery-relation.v1',
            'status':'PASS_EXACT_GROUP_IDENTITY' if exact else 'UNKNOWN',
            'sources':{str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
            'curve':public['curve'],'basis':[list(map(str,p.xy())) for p in basis],
            'recovered_point':list(map(str,P.xy())), 'denominator':int(denominator),'integer_coefficients':integers,
            'maximum_coordinate_denominator_proposed':64,
            'scope':'One384-bit metric proposal per coefficient, denominator at most64, accepted only after exact rational group identity. Public points enter only after the completed search. No numerical relation is promoted without its exact group equation.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true');args=parser.parse_args()
    result=calculate()
    if args.check:assert result==json.loads(OUT.read_text())
    else:
        with OUT.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
    print(result['status'],'denominator',result['denominator'],'word',result['integer_coefficients'])
