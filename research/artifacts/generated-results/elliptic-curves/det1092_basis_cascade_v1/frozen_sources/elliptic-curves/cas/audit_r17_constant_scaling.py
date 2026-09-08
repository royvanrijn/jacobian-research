#!/usr/bin/env python3
"""Finite exact constant-scaling diagnostic for the six compiled R17 models.

Strip primes through 997 from the rational weighted coefficient content and
extract any remaining exact twelfth powers. No factorization, base-parameter
change, specialization sweep, or claim of a globally minimal family is made.
"""
from fractions import Fraction as F
from functools import reduce
from math import gcd, lcm
from pathlib import Path
import argparse
import certify_compact_r17_candidates as cert
from mod2_reduction_independence import _is_prime

ROOT = Path(__file__).resolve().parents[2]


def content(values):
    values = [abs(x) for x in values if x]
    return F(reduce(gcd, (x.numerator for x in values)), reduce(lcm, (x.denominator for x in values)))


def root12(value):
    lo, hi = 0, 1 << ((value.bit_length() + 11) // 12)
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid**12 <= value:
            lo = mid
        else:
            hi = mid
    return hi if hi**12 == value else lo if lo**12 == value else None


def scale_for(a, b):
    ca, cb = content(a), content(b)
    weighted = content([ca**3, cb**2]); num, den = weighted.numerator, weighted.denominator
    scale = F(1); exponents = []
    for p in range(2, 998):
        if not _is_prime(p):
            continue
        exponent = 0
        while num % p == 0:
            num //= p; exponent += 1
        while den % p == 0:
            den //= p; exponent -= 1
        power = exponent // 12
        if power:
            scale *= F(p)**power; exponents.append([p, exponent, power])
    numerator_root, denominator_root = root12(num), root12(den)
    if numerator_root is not None:
        scale *= numerator_root
    if denominator_root is not None:
        scale /= denominator_root
    return scale, {'weighted_content': str(weighted), 'small_prime_exponents': exponents,
                   'remaining_numerator_twelfth_root': numerator_root,
                   'remaining_denominator_twelfth_root': denominator_root}


def bits(values):
    return max(max(abs(q.numerator).bit_length(), q.denominator.bit_length()) for q in values)


def build(output):
    if output.exists():
        raise FileExistsError('use a new immutable diagnostic output')
    rows = []
    for family in ('074d9', '07ca9', '08234', '08f72', '11952', '103b2'):
        path = ROOT / 'artifacts/generated-results' / ('elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json' if family == '074d9'
                                                       else f'elkies-k3-r17-norm12-orbit{family}-direct-fibration-v1.json')
        data = cert.read(path); model = data['representative' if family == '074d9' else 'weierstrass_model']
        a, b = (list(map(F, model[key])) for key in ('A_coefficients_low_to_high', 'B_coefficients_low_to_high'))
        u, diagnostic = scale_for(a, b); new_a, new_b = [v/u**4 for v in a], [v/u**6 for v in b]
        if any(x*u**4 != y for x, y in zip(new_a, a)) or any(x*u**6 != y for x, y in zip(new_b, b)):
            raise ArithmeticError('constant Weierstrass identity failed')
        row = {'family': family, 'source': str(path.relative_to(ROOT)), 'source_sha256': cert.hashed(path),
               'scale_u': str(u), 'before_bits': bits(a+b), 'after_bits': bits(new_a+new_b),
               'A_coefficients_low_to_high': list(map(str, new_a)), 'B_coefficients_low_to_high': list(map(str, new_b)),
               'diagnostic': diagnostic}
        rows.append(row); print(family, row['before_bits'], '->', row['after_bits'], flush=True)
    cert.write(output, {'schema': 'elliptic-curves.r17-constant-scaling-audit.v1', 'rows': rows,
                       'checker_sha256': cert.hashed(Path(__file__).resolve()),
                       'scope': 'Exact constant change x=u^2 X, y=u^3 Y on six literal compiled models. No PGL2 optimization, minimality, generic-rank increase or new-curve claim.'})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', type=Path, required=True)
    build(parser.parse_args().output)
