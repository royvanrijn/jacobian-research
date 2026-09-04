#!/usr/bin/env python3
"""Replay the pinned Nagao rank-17/20 point proofs without repeating discovery.

Only exact finite-group arithmetic on the recorded primes is used. This
does not rerun point search, saturation, numerical heights, or conductor
factorization; those have separate producers and certificates.
"""

from dataclasses import asdict
from fractions import Fraction as Q
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'elliptic-curves/cas'))
from mod2_reduction_independence import (
    combined_mod2_rank, mod2_reduction_signature,
    short_curve_has_no_rational_2_torsion_modular_certificate,
)

ARTIFACTS = ROOT / 'artifacts/generated-results/elliptic-curves'


def replay(model, basis, certificate, expected_rank, signatures_key):
    coefficients = tuple(map(Q, model))
    points = tuple((Q(row['jacobian_x']), Q(row['jacobian_y'])) for row in basis)
    if len(points) != expected_rank:
        raise ArithmeticError('wrong number of pinned points')
    signatures = []
    for pinned in certificate[signatures_key]:
        actual = mod2_reduction_signature(coefficients, points, pinned['prime'])
        # JSON normalization converts the immutable row tuples to lists.
        if json.loads(json.dumps(asdict(actual))) != pinned:
            raise ArithmeticError('finite-reduction signature changed')
        signatures.append(actual)
    if combined_mod2_rank(signatures, len(points)) != expected_rank:
        raise ArithmeticError('the points lack a full-rank mod-2 certificate')
    if not short_curve_has_no_rational_2_torsion_modular_certificate(
            coefficients, certificate['two_torsion_certificate_prime']):
        raise ArithmeticError('the rational 2-torsion exclusion failed')


def main():
    lower17 = json.loads((ARTIFACTS / 'elliptic_nagao_rank17_frontier_certificate.json').read_text())
    if len(lower17['certificates']) != 4:
        raise ArithmeticError('expected four rank-17 frontier curves')
    for row in lower17['certificates']:
        replay(row['short_weierstrass_coefficients'], row['saturated_basis'],
               row['finite_reduction_certificate'], 17, 'signatures')
        print(f"PASS Nagao u={row['parameter_u']}: rank >= 17")
    lower20 = json.loads((ARTIFACTS / 'elliptic_nagao_rank20_t5081_rank20_certificate.json').read_text())
    certificate = lower20['exact_rank_certificate']
    replay(lower20['candidate']['short_weierstrass_coefficients'],
           certificate['saturated_basis'], certificate, 20, 'finite_reduction_signatures')
    print('PASS Nagao T=5081/47: rank >= 20')


if __name__ == '__main__':
    main()
