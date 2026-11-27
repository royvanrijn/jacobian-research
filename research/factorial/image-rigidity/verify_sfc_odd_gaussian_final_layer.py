#!/usr/bin/env python3
"""Exact controls for the next digit layer at final Hermite residues.

The uniform theorem is a written digit and valuation proof. This checker
audits complete coefficient inventories, normalized shape residues,
differential identities, the first-order common coordinate change, and
positive/negative compatibility controls at selected exact inputs.
"""

import json
import platform
from fractions import Fraction as F
from math import comb, factorial, isqrt

from verify_sfc_odd_gaussian_prime_divisors import factorial_tables, pure_hermite


def trim(poly):
    while poly and not poly[-1]:
        poly.pop()
    return poly


def add(a, b, p):
    return trim([((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % p
                 for i in range(max(len(a), len(b)))])


def multiply(a, b, p):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] = (out[i + j] + x * y) % p
    return trim(out)


def derivative(a, p):
    return trim([i * a[i] % p for i in range(1, len(a))])


def gcd(a, b, p):
    a, b = a[:], b[:]
    while b:
        while len(a) >= len(b):
            scalar = a[-1] * pow(b[-1], -1, p) % p
            offset = len(a) - len(b)
            for i, value in enumerate(b):
                a[i + offset] = (a[i + offset] - scalar * value) % p
            trim(a)
        a, b = b, a
    return [x * pow(a[-1], -1, p) % p for x in a] if a else []


def remainder(a, b, p):
    a = a[:]
    while len(a) >= len(b):
        scalar = a[-1] * pow(b[-1], -1, p) % p
        offset = len(a) - len(b)
        for i, value in enumerate(b):
            a[i + offset] = (a[i + offset] - scalar * value) % p
        trim(a)
    return a


def residue(value, p):
    return value.numerator % p * pow(value.denominator, -1, p) % p


def shapes(k):
    falling = [factorial(k) // factorial(k - j) for j in range(k + 1)]
    return ([F(falling[j], factorial(2 * j + 1)) for j in range(k + 1)],
            [F(falling[j], factorial(2 * j)) for j in range(k + 1)])


def check_case(k, p, e, claimed_character=None):
    assert k >= 1 and p > 4 * k + 2 and e >= 1
    assert all(p % divisor for divisor in range(2, isqrt(p) + 1))
    power = p ** e
    q = (2 * k + 1) * power
    values, units = factorial_tables(2 * q + 4, p)
    even = pure_hermite(q + 1, p, values, units)
    odd = pure_hermite(q + 2, p, values, units)
    left = (power + 1) // 2
    ell = (power - 1) // (p - 1)
    lam = F(ell, power)
    beta = F(e, 2) - lam / 2
    gamma = even[0][left] + left * lam
    assert gamma == -beta
    delta = beta if e == 1 else F(p + 1, 2 * p)
    off = 0 if e == 1 else (power // p + 1) // 2
    leading = {left + j * power for j in range(k + 1)}
    secondary = {off + j * power for j in range(k + 1)}
    children = {i + 1 for i in leading}
    threshold = gamma + delta

    def layer(poly, shift):
        return {i + shift for i, v in enumerate(poly[0])
                if v + (i + shift) * lam <= threshold}

    assert layer(even, 0) == leading | secondary
    assert layer(odd, 1) == leading | children | (secondary if e >= 2 else set())
    assert all(even[0][i] + i * lam == gamma for i in leading)
    assert all(even[0][i] + i * lam == threshold for i in secondary)

    s, even_shape = shapes(k)
    z_scale = (-1) ** ell * 2 % p

    def shape_coefficients(poly, index):
        v, unit = poly[0][index], poly[1][index]
        out = []
        for j in range(k + 1):
            i = index + j * power
            assert poly[0][i] == v - j * ell
            out.append(poly[1][i] * pow(unit, -1, p) * pow(z_scale, -j, p) % p)
        return out

    sr = [residue(x, p) for x in s]
    pr = [residue(x, p) for x in even_shape]
    assert shape_coefficients(even, left) == sr
    assert shape_coefficients(odd, left - 1) == [x * pow(2, j, p) % p for j, x in enumerate(sr)]
    assert shape_coefficients(even, off) == pr
    for i in leading:
        assert odd[0][i] == odd[0][i - 1]
        assert odd[1][i] * pow(odd[1][i - 1], -1, p) % p == 2
    chi = 1 if pow(2, (p - 1) // 2, p) == 1 else -1
    if claimed_character is not None:
        chi = claimed_character
    if e >= 2:
        assert shape_coefficients(odd, off - 1) == [x * pow(2, j, p) % p for j, x in enumerate(pr)]
        assert odd[0][off - 1] - odd[0][left - 1] == even[0][off] - even[0][left]
        ratio = (even[1][left] * odd[1][off - 1]
                 * pow(odd[1][left - 1] * even[1][off], -1, p)) % p
        assert ratio == chi % p

    doubled = [x * pow(2, j, p) % p for j, x in enumerate(sr)]
    p_doubled = [x * pow(2, j, p) % p for j, x in enumerate(pr)]
    ds, db = derivative(sr, p), derivative(doubled, p)
    assert gcd(sr, ds, p) == [1]
    common = gcd(sr, doubled, p)
    wronskian = add(multiply(db, sr, p), [-x % p for x in multiply(ds, doubled, p)], p)
    assert remainder(wronskian, multiply(common, common, p), p) == []
    if e == 1:
        test = multiply(db, pr, p)
        margin = 2 * (beta - lam / 2) - delta
    else:
        test = add(multiply(db, pr, p),
                   [(-chi * x) % p for x in multiply(ds, p_doubled, p)], p)
        margin = beta - lam / 2 - delta
    assert margin > 0
    obstruction = gcd(common, test, p)
    if e == 1 or chi == -1:
        assert obstruction == [1]
    else:
        # The residue-prime case must survive this compatibility test.
        assert obstruction == common
    return dict(k=k, prime=p, exponent=e, start=q, delta=str(delta),
                leading_residue_gcd=common, compatibility_gcd=obstruction,
                perturbation_margin=str(margin), character=chi)


def exact_differential_identities():
    for k in range(1, 10):
        s, p = shapes(k)
        for j in range(k + 1):
            assert p[j] == (2 * j + 1) * s[j]
            coefficient = (j - k) * s[j]
            if j < k:
                coefficient += (j + 1) * (4 * j + 6) * s[j + 1]
            assert coefficient == 0


def rational_valuation(value, p):
    numerator, denominator = abs(value.numerator), value.denominator
    result = 0
    while numerator % p == 0:
        numerator //= p
        result += 1
    while denominator % p == 0:
        denominator //= p
        result -= 1
    return result


def straightening_controls():
    count = 0
    for m in (1, 2, 3, 6, 8, 66, 67):
        def moment(power):
            result = F(1)
            for j in range(1, power, 2):
                result *= F(j, m)
            return result
        pure = {i: comb(m, i) * moment(m - i)
                for i in range(m + 1) if (m - i) % 2 == 0}
        actual = {}
        for i in range(m + 1):
            for j in range(m - i + 1):
                linear = m - i - j
                if linear % 2 == 0:
                    actual[i, j] = (F(factorial(m), factorial(i) * factorial(j) * factorial(linear))
                                    * moment(linear + 2 * j))
        difference = dict(actual)
        # Subtract (1-m*x*h) P_m(x+(1+x^2)h) exactly.
        for i, coefficient in pure.items():
            for j in range(i + 1):
                for t in range(j + 1):
                    value = coefficient * comb(i, j) * comb(j, t)
                    key = (i - j + 2 * t, j)
                    difference[key] = difference.get(key, F(0)) - value
                    key = (key[0] + 1, key[1] + 1)
                    difference[key] = difference.get(key, F(0)) + m * value
        difference = {key: value for key, value in difference.items() if value}
        assert all(j >= 2 for i, j in difference)
        alpha, beta = F(1, 26), F(6, 13)
        gamma = min(rational_valuation(c, 13) + i * alpha for i, c in pure.items())
        assert all(rational_valuation(c, 13) + i * alpha + j * beta
                   >= gamma + 2 * (beta - alpha)
                   for (i, j), c in difference.items())
        count += len(difference)
    return count


def main():
    exact_differential_identities()
    inputs = [(2, 13, 1), (2, 13, 2), (2, 13, 3), (2, 13, 4),
              (3, 37, 1), (3, 37, 2), (4, 89, 1), (4, 89, 2),
              (1, 7, 2), (2, 11, 2)]
    cases = [check_case(*item) for item in inputs]
    # Explicit same-residue incompatibility for the requested family.
    s, even_shape = shapes(2)
    evaluate = lambda poly, z: sum(c * z ** j for j, c in enumerate(poly))
    ds = [j * s[j] for j in range(1, len(s))]
    value = (2 * evaluate(ds, 4) * evaluate(even_shape, 2)
             + evaluate(ds, 2) * evaluate(even_shape, 4))
    assert value == F(254, 45) and residue(value, 13) == 12
    assert residue(-evaluate(even_shape, 2) / evaluate(ds, 2), 13) == 9
    assert residue(evaluate(even_shape, 4) / (2 * evaluate(ds, 4)), 13) == 4
    # Exact Wronskian identity: 90*(B'*A-A'*B)=(Z-2)^2+13*(Z+2).
    doubled = [c * 2 ** j for j, c in enumerate(s)]
    db = [j * doubled[j] for j in range(1, len(doubled))]
    wronskian = [F(0)] * 4
    for i, c in enumerate(db):
        for j, d in enumerate(s):
            wronskian[i + j] += c * d
    for i, c in enumerate(ds):
        for j, d in enumerate(doubled):
            wronskian[i + j] -= c * d
    assert [90 * x for x in wronskian] == [30, 9, 1, 0]
    # The positive-character example must remain unresolved by this layer.
    residue_control = check_case(4, 89, 2)
    assert len(residue_control["compatibility_gcd"]) > 1
    try:
        check_case(2, 13, 2, claimed_character=1)
    except AssertionError:
        pass
    else:
        raise AssertionError("A corrupted inter-row phase was accepted")
    print(json.dumps(dict(
        status="PASS", python=platform.python_version(),
        scope="finite exact controls of the written final-layer theorem; not an all-order experiment",
        cases=cases, straightening_remainder_coefficients=straightening_controls(),
        exact_13_certificate="254/45 = 12 modulo 13; incompatible lift corrections 9 and 4",
        exact_wronskian="90*(B'*A-A'*B)=(Z-2)^2+13*(Z+2)",
        retained_residue_control="q=9*89^e, e>=2: common residue survives this first compatibility layer",
        corruption_control="incorrect positive inter-row phase at p=13 rejected by exact coefficient residues",
        assurance="internal replay; no independent review or formal verification",
    ), indent=2))


if __name__ == "__main__":
    main()
