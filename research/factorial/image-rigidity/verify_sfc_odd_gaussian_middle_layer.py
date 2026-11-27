#!/usr/bin/env python3
"""Finite exact normalized-layer controls for n=6*31**e.

Inspects every pure odd coefficient below relative valuation one at e=2,3.
Divides by the entire leading lower-row unit, including its children.
This is not an all-exponent nonlifting certificate.
"""

import json
import platform
from fractions import Fraction as F

from verify_sfc_odd_gaussian_prime_divisors import factorial_tables, pure_hermite

P = 31
Z = 7
SHAPES = ([1, 17, 8, 23], [1, 10, 15])


def strip(number):
    valuation = 0
    while number % P == 0:
        number //= P
        valuation += 1
    return valuation, number


def recurrence(order, modulus):
    """Independent coefficient recurrence, without factorial tables."""
    half = (order - 1) // 2
    valuations, units = [0], [1]
    for r in range(half):
        vn, un = strip(2 * order * (half - r))
        vd, ud = strip((2 * r + 2) * (2 * r + 3))
        valuations.append(valuations[-1] + vn - vd)
        units.append(units[-1] * un * pow(ud, -1, modulus) % modulus)
    return valuations, units


def evaluate(poly, z=Z, modulus=P):
    return sum(c * pow(z, j, modulus) for j, c in enumerate(poly)) % modulus


def derivative(poly):
    return [j * poly[j] for j in range(1, len(poly))]


def divide_series(raw, unit, cutoff):
    """Formal division in the filtered associated graded ring at Z=7."""
    assert unit[0] == 1
    positive = [(w, c) for w, c in unit.items() if w and c]
    quotient = {}
    for w in range(cutoff):
        c = (raw.get(w, 0) - sum(
            d * quotient.get(w - u, 0) for u, d in positive if u <= w
        )) % P
        if c:
            quotient[w] = c
    # Multiplication back checks every retained weight, including zero layers.
    product = {}
    for w, c in quotient.items():
        for u, d in unit.items():
            if w + u < cutoff:
                product[w + u] = (product.get(w + u, 0) + c * d) % P
    assert all(product.get(w, 0) == raw.get(w, 0) for w in range(cutoff))
    return quotient


def check_layers(exponent, omit_child=False):
    power = P ** exponent
    middle = 6 * power
    ell = (power - 1) // (P - 1)
    lam = F(ell, power)
    base_minus = (power - 1) // 2
    scale = (-1) ** ell * 2
    # In the graded ring put p=t**power, Y=y*t**ell. Then Z=scale*y**power.
    # Z=7 is rational over F_31, so inverse Frobenius fixes it.
    y = Z * pow(scale, -1, P) % P
    values, units = factorial_tables(2 * middle + 2, P)
    rows, unit, inventories = [], {0: 1}, []
    for side, order, base in [(0, middle + 1, 0),
                              (1, middle - 1, base_minus)]:
        vals, coeffs = pure_hermite(order, P, values, units)
        assert (vals, coeffs) == recurrence(order, P)
        gamma = vals[base] + base * lam
        assert gamma == (0 if side == 0 else -F(exponent, 2) - lam / 2)
        row, retained, leading_groups = {}, 0, 0
        for b in range(power):
            gap = vals[b] + b * lam - gamma
            if gap >= 1:
                continue
            weight = gap * power
            assert weight.denominator == 1 and weight >= 0
            weight = int(weight)
            shape = []
            for j, r in enumerate(range(b, len(vals), power)):
                assert vals[r] + r * lam - gamma == gap
                shape.append(coeffs[r] * pow(coeffs[base], -1, P)
                             * pow(scale, -j, P) % P)
                retained += 1
            factor = pow(y, b - base, P)
            row[weight] = (row.get(weight, 0) + evaluate(shape) * factor) % P
            if side == 1 and b >= base:
                assert shape == [shape[0] * c % P for c in SHAPES[1]]
                unit[weight] = shape[0] * factor % P
                leading_groups += 1
            if omit_child and side == 1 and b == (P - 1) // 2 + 1:
                row[weight] = 0
        rows.append(row)
        inventories.append(dict(order=order, retained_coefficients=retained,
                                leading_unit_groups=leading_groups))

    normalized = divide_series(rows[1], unit, power)
    events = []
    transverse = []
    for w in sorted(set(rows[0]) | set(normalized)):
        vector = (rows[0].get(w, 0), normalized.get(w, 0))
        if not any(vector):
            continue
        determinant = (3 * vector[0] - 7 * vector[1]) % P
        event = dict(gap=str(F(w, power)), vector=list(vector),
                     determinant=determinant)
        events.append(event)
        if determinant:
            transverse.append(event)

    if omit_child:
        assert transverse and F(transverse[0]["gap"]) < F(945, 961)
        return transverse[0]

    assert [F(event["gap"]) for event in events[:15]] == [F(j, P) for j in range(16, 31)]
    assert all(event["determinant"] == 0 for event in events[:15])
    if exponent == 2:
        assert transverse == [dict(gap="945/961", vector=[3, 0], determinant=9)]
        assert len(events) == 16
        # eps=p^16*Y^-465, T=Y^31/p. Remove the fixed graded residue of T.
        formal = [[event["vector"][s] * pow(y, -v, P) % P for s in (0, 1)]
                  for v, event in enumerate(events[:15])]
        assert formal[:2] == [[5, 11], [6, 7]]
        assert [c * pow(10, -1, P) % P for c in formal[1]] == [13, 10]
        # Before dividing by 1+15*T+..., the second vector falsely separates.
        raw_second = rows[1][17 * power // P] * pow(y, -1, P) % P
        assert raw_second == 17 and (3 * 6 - 7 * raw_second) % P == 23
    else:
        assert exponent == 3 and not transverse and len(events) == 15
        formal = None
    return dict(exponent=exponent, middle=middle, cutoff="relative weight < 1",
                inventories=inventories, leading_unit_terms=len(unit),
                surviving_layers=events, formal_T_vectors=formal)


def factorial_unit_mod_square(number):
    """Wilson block recursion modulo p^2; runs without an exponent-sized table."""
    modulus = P * P
    prefixes = [1]
    for j in range(1, modulus + 1):
        prefixes.append(prefixes[-1] * (j if j % P else 1) % modulus)
    valuation, unit = 0, 1
    while number:
        blocks, tail = divmod(number, modulus)
        unit = unit * pow(prefixes[-1], blocks, modulus) * prefixes[tail] % modulus
        number //= P
        valuation += number
    return valuation, unit


def leading_rows_mod_square(exponent):
    modulus = P * P
    power = P ** exponent
    middle = 6 * power
    ell = (power - 1) // (P - 1)
    scale = (-1) ** ell * 2

    def coefficient(order, r):
        half = (order - 1) // 2
        vn, un = factorial_unit_mod_square(half)
        vd, ud = factorial_unit_mod_square(half - r)
        vf, uf = factorial_unit_mod_square(2 * r + 1)
        return vn - vd - vf, pow(2 * order, r, modulus) * un * pow(ud * uf, -1, modulus) % modulus

    rows = []
    for order, base, top in [(middle + 1, 0, 3),
                              (middle - 1, (power - 1) // 2, 2)]:
        bv, bu = coefficient(order, base)
        row = []
        for j in range(top + 1):
            v, u = coefficient(order, base + j * power)
            assert v - bv == -j * ell
            row.append(u * pow(bu, -1, modulus) * pow(scale, -j, modulus) % modulus)
        rows.append(row)
        if exponent <= 3:
            rv, ru = recurrence(order, modulus)
            assert row == [ru[base + j * power] * pow(ru[base], -1, modulus)
                           * pow(scale, -j, modulus) % modulus for j in range(top + 1)]
    assert rows == [[1, 203, 628, 612], [1, 506, 542]]
    a, b = rows
    wronskian = (evaluate(derivative(b), modulus=modulus) * evaluate(a, modulus=modulus)
                 - evaluate(derivative(a), modulus=modulus) * evaluate(b, modulus=modulus)) % modulus
    assert wronskian == 31 * 17
    return dict(exponent=exponent, leading_rows=rows, wronskian_at_7=wronskian)


def main():
    assert [evaluate(poly) for poly in SHAPES] == [0, 0]
    assert [evaluate(derivative(poly)) for poly in SHAPES] == [7, 3]
    cases = [check_layers(e) for e in (2, 3)]
    corrupted = check_layers(2, omit_child=True)
    print(json.dumps(dict(
        status="PASS", python=platform.python_version(), prime=P, common_residue=Z,
        scope="complete pure-array inventories below weight one at e=2,3; leading rows modulo 31^2 at e=2,...,6",
        cases=cases, leading_mod_square=[leading_rows_mod_square(e) for e in range(2, 7)],
        corruption_control=dict(omitted_first_secondary_child=corrupted,
                                omitted_row_unit="false second-layer determinant 23"),
        remaining_obligation="uniform all-exponent normalized-layer control and valuation-one compatibility",
        assurance="internal exact replay; no independent review or formal verification",
    ), indent=2))


if __name__ == "__main__":
    main()
