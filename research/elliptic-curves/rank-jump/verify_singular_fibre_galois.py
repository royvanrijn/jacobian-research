#!/usr/bin/env python3
"""Independent Sage arithmetic replay of the S24/j-map proof inputs."""
import argparse
from pathlib import Path
import retrospective as r
import singular_fibre_galois as source
import parameter_cover_capacity as capacity
import verify_surface_discriminant_modular as modular

OUTPUT = r.OUT / 'rank_jump_singular_fibre_galois_verification_v1.json'


def compute():
    from sage.all import QQ, ZZ, PolynomialRing, Permutation, binomial, factorial
    data = r.read(source.OUTPUT)
    for name, digest in data['bindings'].items():
        assert r.digest((r.ROOT/name).read_bytes()) == digest
    assert modular.compute() == r.read(modular.OUTPUT)
    inp = r.read(capacity.INPUT)
    R = PolynomialRing(QQ, 'u')
    A, B = [R(list(map(QQ, inp[k]))) for k in ['A', 'B']]
    numerator, delta = (-48*A)**3, -16*(4*A**3+27*B**2)
    j = R.fraction_field()(numerator/delta)
    assert numerator.gcd(delta) == 1
    assert delta.gcd(delta.derivative()) == 1
    assert delta.degree() == numerator.degree() == 24
    assert j.numerator().degree() == j.denominator().degree() == 24
    assert list(map(str, numerator)) == data['j_map']['numerator_coefficients_ascending']
    assert list(map(str, delta)) == data['j_map']['denominator_coefficients_ascending']
    assert numerator.leading_coefficient()/delta.leading_coefficient() == QQ(data['j_map']['value_at_infinity'])

    # A representative of the certified Frobenius conjugacy class.
    g = Permutation([1] + list(range(3, 8)) + [2] + list(range(9, 25)) + [8])
    assert sorted(g.cycle_type()) == [1, 6, 17] and g.signature() == -1
    powered = Permutation(g.to_permutation_group_element()**6)
    assert sorted(powered.cycle_type()) == [1]*7 + [17]
    assert ZZ(17).is_prime() and 24-17 >= 3
    divisors = [s for s in ZZ(24).divisors() if 1 < s < 24]
    assert [int(s) for s in divisors] == [x['block_size'] for x in data['galois_certificate']['block_exclusions']]
    # A 17-cycle cannot act nontrivially on the blocks or within a block.
    assert all(factorial(s) % 17 and factorial(24//s) % 17 for s in divisors)
    for row in data['subset_field_degrees']:
        assert row['degree_of_field_of_monic_factor'] == binomial(24, row['subset_size'])
    for row in data['partition_orbit_sizes']:
        s, m = row['block_size'], row['number_of_blocks']
        assert row['orbit_size'] == factorial(24)/(factorial(s)**m*factorial(m))
    return {
        'schema': 'rank-jump.singular-fibre-galois-verification.v1',
        'bindings': {str(p.relative_to(r.ROOT)): r.digest(p.read_bytes()) for p in (Path(__file__), source.OUTPUT)},
        'status': 'PASS', 'exact_j_degree': 24, 'reduced_simple_pole_divisor': True,
        'frobenius_power_checked_by_permutation_arithmetic': True,
        'all_proper_block_sizes_excluded': True,
        'theorem_boundary': 'Jordan theorem and the written rational-composition proof are mathematical dependencies, not enumerated group searches.'
    }


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('mode', choices=['build', 'check'])
    args = p.parse_args()
    result = compute()
    if args.mode == 'check':
        assert r.read(OUTPUT) == result
        print('PASS independent exact j-map and permutation replay')
    else:
        r.write_new(OUTPUT, result)
        print(result['status'])
