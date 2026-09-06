#!/usr/bin/env python3
"""A coefficient-only Galois and rational-composition gate for original R17."""
import argparse
from fractions import Fraction as F
from math import comb, factorial, gcd
from pathlib import Path

import retrospective as r
import parameter_cover_capacity as capacity
import surface_discriminant_irreducibility as modular
import verify_surface_discriminant_modular as verifier

PROTOCOL = Path(__file__).with_name('SINGULAR_FIBRE_GALOIS_PROTOCOL.json')
OUTPUT = r.OUT / 'rank_jump_singular_fibre_galois_v1.json'


def trim(a):
    a = list(a)
    while len(a) > 1 and not a[-1]:
        a.pop()
    return a


def multiply(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return trim(out)


def remainder(a, b):
    a, b = trim(a), trim(b)
    assert b != [0]
    while a != [0] and len(a) >= len(b):
        shift, c = len(a)-len(b), a[-1]/b[-1]
        for i, v in enumerate(b):
            a[shift+i] -= c*v
        a = trim(a)
    return a


def polynomial_gcd(a, b):
    a, b = list(map(F, a)), list(map(F, b))
    while b != [0]:
        a, b = b, remainder(a, b)
    return [x/a[-1] for x in a]


def power_cycle_type(lengths, exponent):
    return sorted(d//gcd(d, exponent) for d in lengths
                  for _ in range(gcd(d, exponent)))


def primitive_gate(n, prime_cycle):
    assert prime_cycle > 1
    assert all(prime_cycle % d for d in range(2, prime_cycle))
    rows = [{'block_size': s, 'number_of_blocks': n//s,
             'excluded_by_prime_cycle': prime_cycle > s and prime_cycle > n//s}
            for s in range(2, n) if n % s == 0]
    return rows, all(x['excluded_by_prime_cycle'] for x in rows)


def compute():
    checked = verifier.compute()
    assert checked == r.read(verifier.OUTPUT)
    assert checked['status'] == 'PASS_IRREDUCIBLE_OVER_Q'
    row = next(x for x in checked['rows'] if x['prime'] == 181)
    cycles = row['factor_degrees']
    assert cycles == [1, 6, 17] and sum(cycles) == 24
    powered = power_cycle_type(cycles, 6)
    assert powered == [1]*7 + [17]
    blocks, primitive = primitive_gate(24, 17)
    assert primitive and 17 <= 24-3
    sign = (-1)**sum(d-1 for d in cycles)
    assert sign == -1

    inp = r.read(capacity.INPUT)
    A, B = [list(map(F, inp[k])) for k in ['A', 'B']]
    cubeA, squareB = multiply(multiply(A, A), A), multiply(B, B)
    delta = [-64*a-432*b for a, b in zip(cubeA, squareB, strict=True)]
    primitive_delta = list(map(F, r.read(modular.OUTPUT)['primitive_coefficients_ascending']))
    content = F(r.read(modular.OUTPUT)['integer_content'])
    assert delta == [content*x for x in primitive_delta]
    numerator = [(-48)**3*x for x in cubeA]
    common = polynomial_gcd(numerator, delta)
    assert common == [1] and len(numerator) == len(delta) == 25
    # No pole at infinity: both leading coefficients are nonzero.
    infinity_value = numerator[-1]/delta[-1]
    assert infinity_value

    paths = (Path(__file__), PROTOCOL, capacity.INPUT, modular.OUTPUT,
             verifier.OUTPUT, Path(verifier.__file__))
    return {
        'schema': 'rank-jump.singular-fibre-galois.v1',
        'bindings': {str(p.relative_to(r.ROOT)): r.digest(p.read_bytes()) for p in paths},
        'status': 'PASS_S24_AND_RATIONAL_J_INDECOMPOSABLE',
        'galois_certificate': {
            'degree': 24, 'transitive': True, 'frobenius_prime': 181,
            'frobenius_cycle_type': cycles, 'power': 6, 'powered_cycle_type': powered,
            'block_exclusions': blocks, 'primitive': primitive,
            'jordan_hypotheses': {'prime_cycle_length': 17, 'fixed_points': 7},
            'odd_frobenius_sign': sign, 'group': 'S24',
            'theorem': 'Jordan prime-cycle theorem, as stated in Jones, Theorem 1.1, arXiv:1209.5169'
        },
        'j_map': {
            'numerator_coefficients_ascending': list(map(str, numerator)),
            'denominator_coefficients_ascending': list(map(str, delta)),
            'exact_monic_gcd': list(map(str, common)), 'degree': 24,
            'finite_simple_poles': 24, 'value_at_infinity': str(infinity_value),
            'proper_rational_composition_over_Q': False,
            'proper_composition_over_Qbar': 'UNKNOWN',
            'proof': 'Every proper rational composition partitions the reduced pole divisor into equally sized fibres, contradicting primitive constant-field Galois action.'
        },
        'subset_field_degrees': [{'subset_size': k, 'degree_of_field_of_monic_factor': comb(24, k)}
                                 for k in [1, 2, 3, 4, 12]],
        'partition_orbit_sizes': [
            {'block_size': s, 'number_of_blocks': 24//s,
             'orbit_size': factorial(24)//(factorial(s)**(24//s)*factorial(24//s))}
            for s in range(2, 24) if 24 % s == 0],
        'rational_base_mobius_symmetry_preserving_singular_divisor': 'identity_only',
        'solvable_galois_constant_extension': 'singular polynomial remains irreducible; no nontrivial block system',
        'boundary': 'No Mordell-Weil or Selmer rank conclusion. No comparison of high and low specializations; this obstruction is shared by all fibres of the original R17 family.'
    }


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('mode', choices=['build', 'check'])
    args = p.parse_args()
    result = compute()
    if args.mode == 'check':
        assert r.read(OUTPUT) == result
        print('PASS S24, exact j pole divisor, and rational composition gate')
    else:
        r.write_new(OUTPUT, result)
        print(result['status'])
