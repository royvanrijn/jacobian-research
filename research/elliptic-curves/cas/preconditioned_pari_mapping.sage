"""Exact factor-free preconditioning, then PARI minimization at2 and3 only.

This is a coordinate policy, not a claim of global minimality. Run under an
external map-worker limit. All composed rational map identities are checked.
"""
import sys
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
from sage.all import pari

CAS = Path(__file__).resolve().parent
sys.path.insert(0, str(CAS))
from v3_warm_engine import load
from search_observability import multiply, transform


def mapping(model, points, centre):
    ff = load('preconditioned_factor_free', CAS/'factor_free_pari_mapping.sage')
    old = ff.mapping(model, points, centre)
    poly = lambda values: '+'.join(f'({F(x)})*x^{i}' for i, x in enumerate(values))
    C = '['+poly(old['reduced_P'])+','+poly(old['reduced_Q'])+']'
    result = pari('my(m1,m2,C1,C2);C1=hyperellminimalmodel('+C+',&m1,[2,3]);'
                  'C2=hyperellred(C1,&m2);[C2,m1,m2]')
    matrices = [tuple(F(str(m[i,j])) for i in range(2) for j in range(2))
                for m in (result[1][1], result[2][1])]
    first = tuple(map(F, old['matrix']))
    second = multiply(*matrices)
    M = multiply(first, second)
    P = [F(str(result[0][0].polcoef(i))) for i in range(5)]
    Q = [F(str(result[0][1].polcoef(i))) for i in range(3)]
    disc = [4*P[i]+sum(Q[j]*Q[i-j] for j in range(3) if 0 <= i-j < 3) for i in range(5)]
    raw = tuple(map(F, old['raw_coefficients']))
    transformed = transform(raw, M)
    k = next(i for i, x in enumerate(disc) if x)
    ratio = transformed[k]/disc[k]
    if (ratio <= 0 or any(a != ratio*b for a,b in zip(transformed,disc)) or
            isqrt(ratio.numerator)**2 != ratio.numerator or
            isqrt(ratio.denominator)**2 != ratio.denominator):
        raise ArithmeticError('preconditioned map identity failed')
    return dict(centre=centre, raw_coefficients=old['raw_coefficients'],
                denominator_clearing=old['denominator_clearing'],
                first_matrix=list(map(str,first)), second_matrix=list(map(str,second)),
                matrix=list(map(str,M)), reduced_P=list(map(str,P)), reduced_Q=list(map(str,Q)),
                discriminant_quartic=list(map(str,disc)), square_ratio=str(ratio),
                coordinate_policy={'kind':'raw','matrix':list(map(str,M))},
                minimization_primes=[2,3], claim_boundary='Exact coordinates; no global minimality claim.')


if __name__ == '__main__':
    import argparse, hashlib, json
    from research_runtime.store import checkpoint
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('preserve previous map result')
    raw = args.input.read_bytes()
    data = json.loads(raw)
    pari.allocatemem(256000000,silent=True)
    result = mapping(tuple(map(F,data['curve'])),
                     tuple(tuple(map(F,p)) for p in data['points']), data['centre'])
    checkpoint(args.output, {'status':'EXACT_PRECONDITIONED_MAP_CONSTRUCTED', 'mapping':result,
               'input_sha256':hashlib.sha256(raw).hexdigest(),
               'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    print('EXACT_PRECONDITIONED_MAP_CONSTRUCTED',flush=True)
