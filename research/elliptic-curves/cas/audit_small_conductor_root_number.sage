#!/usr/bin/env sage-python
"""Exact local root-number diagnostic on the already factored rank22 curve.

The global sign is analytic scheduling information, not an unconditional
rank-parity or new-point certificate. No L-series summation is attempted.
"""
import argparse
from pathlib import Path
import sys
from sage.all import QQ, ZZ, EllipticCurve, prod
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'elliptic-curves/cas'))
import certify_compact_r17_candidates as cert
INPUT = ROOT/'artifacts/generated-results/elliptic-curves/next12_rank22_exact_conductor_v1.json'

def build(output):
    if output.exists():
        raise FileExistsError('preserve root-number diagnostic')
    data = cert.read(INPUT)
    factors = [(ZZ(p), int(e)) for p, e in data['discriminant_factorization']]
    E = EllipticCurve(QQ, [QQ(a) for a in data['integral_model']])
    if prod(p**e for p, e in factors) != abs(E.discriminant()) or any(not p.is_prime(proof=True) for p, e in factors):
        raise ArithmeticError('complete proven factorization required')
    rows = [{'prime': str(p), 'root_number': int(E.root_number(p))} for p, e in factors]
    if any(r['root_number'] not in (-1, 1) for r in rows):
        raise ArithmeticError('invalid local sign')
    result = {'schema': 'elliptic-curves.small-conductor-root-number.v1',
              'sources': {str(p.relative_to(ROOT)): cert.hashed(p) for p in (INPUT, Path(__file__).resolve())},
              'infinite_root_number': -1, 'bad_prime_root_numbers': rows,
              'global_root_number': -int(prod(r['root_number'] for r in rows)),
              'claim_boundary': 'Exact functional-equation sign from local factors. It does not prove algebraic rank parity, an extra point, exact rank or an upper bound.'}
    cert.write(output, result)
    print('ROOT NUMBER', result['global_root_number'], rows, flush=True)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    build(p.parse_args().output.resolve())
