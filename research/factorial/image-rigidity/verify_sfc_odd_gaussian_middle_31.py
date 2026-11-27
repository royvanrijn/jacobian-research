#!/usr/bin/env python3
"""Exact finite audit of the uniform middle-31 valuation-one proof.

The unbounded exponent assertion is the written digit/block argument.
This audits its fixed F_31 identity and its coefficient formulas at e=3.
"""

import json
import math
import platform
from fractions import Fraction as F

import verify_sfc_odd_gaussian_middle_layer as control
from verify_sfc_odd_gaussian_prime_divisors import factorial_tables, pure_hermite

P = 31
MOD = P * P


def multiply(a, b, cutoff):
    result = [0] * cutoff
    for i, x in enumerate(a):
        for j, y in enumerate(b[:cutoff - i]):
            if i + j < cutoff:
                result[i + j] = (result[i + j] + x * y) % P
    return result


def fixed_identity():
    denominator = [pow(math.factorial(j), -1, P) for j in range(16)]
    numerator = [0] * 32
    upper = [0] * 32
    for s in range(1, 16):
        numerator[31 - s] = (-1) ** (s - 1) * math.factorial(s - 1) % P
    for v in range(16, 31):
        upper[v] = (3 * pow(4, v, P)
                    * pow(math.factorial(31 - v) * math.factorial(2 * v - 31), -1, P)) % P
    quotient = [0] * 32
    for k in range(32):
        quotient[k] = (numerator[k] - sum(
            denominator[j] * quotient[k - j] for j in range(1, min(15, k) + 1)
        )) % P
    expected = [-c * pow(6, -1, P) % P for c in upper]
    expected[31] = 12
    assert quotient == expected
    assert multiply(denominator, quotient, 32) == numerator
    corrupted = quotient[:]
    corrupted[31] = 0
    assert multiply(denominator, corrupted, 32) != numerator
    return dict(denominator=denominator, numerator=numerator,
                upper_scalar=upper, quotient=quotient,
                boundary=quotient[31], omitted_boundary_rejected=True)


def leading_block_formula():
    wilson = math.factorial(30) % MOD
    harmonic = sum(pow(j, -1, P) for j in range(1, 16)) % P
    kappa = -pow(2, 30, MOD) * pow(wilson, -1, MOD) % MOD
    eta = kappa * (1 + P * harmonic) % MOD
    a = [1, 3 * pow(2, -1, MOD), pow(4, -1, MOD), pow(120, -1, MOD)]
    b = [1, -pow(3, -1, MOD), pow(60, -1, MOD)]
    a = [c * pow(kappa, j, MOD) % MOD for j, c in enumerate(a)]
    b = [c * pow(eta, j, MOD) % MOD for j, c in enumerate(b)]
    assert (harmonic, kappa, eta) == (19, 776, 404)
    assert [a, b] == [[1, 203, 628, 612], [1, 506, 542]]
    assert control.leading_rows_mod_square(3)["leading_rows"] == [a, b]
    av, bv = [control.evaluate(row, modulus=MOD) for row in (a, b)]
    assert av % P == bv % P == 0
    constants = [av // P, bv // P]
    assert constants == [29, 10]
    return dict(wilson=wilson, harmonic_half=harmonic, kappa=kappa, eta=eta,
                rows=[a, b], constants=constants)


def coefficient_controls():
    exponent = 3
    power = P ** exponent
    middle = 6 * power
    base = (power - 1) // 2
    child = power // P
    ell = (power - 1) // 30
    lam = F(ell, power)
    scale = (-1) ** ell * 2
    y = 7 * pow(scale, -1, P) % P
    factorial_values, factorial_units = factorial_tables(2 * middle + 2, P)
    plus = control.recurrence(middle + 1, P)
    minus = control.recurrence(middle - 1, P)
    assert plus == pure_hermite(middle + 1, P, factorial_values, factorial_units)
    assert minus == pure_hermite(middle - 1, P, factorial_values, factorial_units)
    pv, pu = plus
    mv, mu = minus
    inverse_base = pow(mu[base], -1, P)
    retained_upper, retained_lower, upper_groups = 0, 0, []
    lower_scalar, upper_values = {}, {}
    unit = {}
    for b in range(power):
        gap = pv[b] + b * lam
        if gap <= 1:
            upper_groups.append(b)
            shape = []
            for j, r in enumerate(range(b, len(pv), power)):
                assert pv[r] + r * lam == gap
                shape.append(pu[r] * pow(pu[b], -1, P) * pow(scale, -j, P) % P)
                retained_upper += 1
            expected = [1, 17, 8, 23] if b == 0 else [1, 21, 15]
            assert shape == expected
            if b:
                v = b // child
                child_ell = (child - 1) // 30
                assert pv[b] == -v * child_ell
                polynomial_coefficient = (3 * pow(4, v, P) * pow(
                    math.factorial(31 - v) * math.factorial(2 * v - 31), -1, P)) % P
                x_unit = (-1) ** child_ell * pow(2, -child, P) % P
                assert pu[b] == polynomial_coefficient * pow(x_unit, v, P) % P
            upper_values[int(gap * power)] = (
                pu[b] * control.evaluate(shape) * pow(y, b, P) % P)

        gap = mv[b] - mv[base] + (b - base) * lam
        if gap <= 1:
            shape = []
            for j, r in enumerate(range(b, len(mv), power)):
                assert mv[r] - mv[base] + (r - base) * lam == gap
                shape.append(mu[r] * pow(mu[b], -1, P) * pow(scale, -j, P) % P)
                retained_lower += 1
            expected = [1, 10, 15] if b >= base else [1, 30, 13]
            assert shape == expected
            scalar = mu[b] * inverse_base % P
            weight = int(gap * power)
            scalar = scalar * pow(y, b - base, P) % P
            if b >= base:
                unit[weight] = (unit.get(weight, 0) + scalar) % P
            else:
                lower_scalar[weight] = (lower_scalar.get(weight, 0) + scalar) % P

    assert upper_groups == [0] + [v * child for v in range(16, 31)]
    # Every scalar coefficient in (6), including those outside the cutoff.
    for c in range(base + 1):
        expected_value = -factorial_values[c]
        expected_unit = pow(2, -c, P) * pow(factorial_units[c], -1, P) % P
        assert mv[base + c] - mv[base] == expected_value
        assert mu[base + c] * inverse_base % P == expected_unit
    for d in range(1, base + 1):
        expected_value = exponent + factorial_values[d - 1]
        expected_unit = pow(-2, d - 1, P) * factorial_units[d - 1] % P
        assert mv[base - d] - mv[base] == expected_value
        assert mu[base - d] * inverse_base % P == expected_unit

    # Full division through and including weight one. The proof replaces this
    # exponent-sized calculation by the fixed identity in fixed_identity().
    normalized = control.divide_series(lower_scalar, unit, power + 1)
    assert all((3 * upper_values.get(w, 0) - 7 * 11 * normalized.get(w, 0)) % P == 0
               for w in range(power))
    assert normalized.get(power, 0) == 6
    assert 11 * normalized[power] % P == 4
    # Recover that boundary directly from the previously verified <1 inventory.
    earlier = control.check_layers(exponent)
    assert len(earlier["surviving_layers"]) == 15
    boundary = -sum(unit.get(power - v * child, 0)
                    * upper_values[v * child] * 3 * pow(7, -1, P)
                    for v in range(16, 31)) % P
    assert boundary == 4
    return dict(exponent=exponent, upper_groups=upper_groups,
                retained_upper_coefficients=retained_upper,
                retained_lower_coefficients=retained_lower,
                scalar_coefficients_compared=2 * base + 1,
                unit_terms=len(unit), normalized_scalar_boundary=normalized[power],
                normalized_lower_boundary=boundary,
                full_division_multiplied_back=True,
                independent_factorial_and_recurrence_arrays=True)


def main():
    identity = fixed_identity()
    leading = leading_block_formula()
    coefficients = coefficient_controls()
    a, b = leading["constants"]
    leading_projection = (3 * a - 7 * b) % P
    child_projection = -7 * coefficients["normalized_lower_boundary"] % P
    total = (leading_projection + child_projection) % P
    assert (leading_projection, child_projection, total) == (17, 3, 20)
    print(json.dumps(dict(
        status="PASS", python=platform.python_version(), prime=P,
        scope="fixed polynomial and prime-square identities; full e=3 coefficient controls through weight one",
        polynomial_identity=identity, leading_block_formula=leading,
        coefficient_controls=coefficients,
        compatibility=dict(leading=leading_projection, row_unit=child_projection, total=total),
        uniform_scope="e>=3 follows from the written digit cutoff, cancellation and block-recursion proof; e=1,2 use the existing prefix",
        assurance="internal exact replay; no independent review, external review or formal verification",
    ), indent=2))


if __name__ == "__main__":
    main()
