#!/usr/bin/env python3
"""Verify conductor-18 reconstruction certificates over Q with SageMath."""

from collections import defaultdict
import json
from math import comb
from pathlib import Path

from sage.all import PolynomialRing, QQ, matrix


def phase_terms(odd_terms, even_radial=0):
    # Keys are the angular exponent and the weight deficit.
    terms = {(2 * j - 1, 0): (-1)**j * comb(4, j) for j in range(5)}
    terms[-1, 2] = -1
    terms.update({(j, 1): value for j, value in enumerate(odd_terms)})
    terms[0, 2] = even_radial
    return {key: value for key, value in terms.items() if value}


def reconstruction_rows(phase, columns, final_shift, ring):
    coefficients = {(0, 0): ring.one()}
    for k in range(final_shift + 1):
        row = []
        for j, parity in columns:
            entry = sum(
                value / QQ(13 - parity + 18 * k - deficit)
                for (angle, deficit), value in coefficients.items()
                if angle == -j and deficit % 2 == parity
            )
            row.append(ring(entry))
        yield row
        if k == final_shift:
            break
        product = defaultdict(ring.zero)
        for (angle, deficit), value in coefficients.items():
            for (step, weight), coefficient in phase.items():
                # Remaining factors can each lower the angle by only one.
                if angle + step <= final_shift - k - 1:
                    product[angle + step, deficit + weight] += value * coefficient
        coefficients = {key: value for key, value in product.items() if value}


def verify_exceptional_radial_line():
    ring = PolynomialRing(QQ, "b")
    b = ring.gen()
    phase = phase_terms((b, 0, 0, 1))
    columns = [(j, 0) for j in range(7)] + [(j, 1) for j in range(6)]
    rows = list(reconstruction_rows(phase, columns, 13, ring))
    delta = matrix(ring, rows[:13]).det()
    replacement = matrix(ring, rows[:12] + rows[13:]).det()
    assert delta.degree() == 42 and delta(0) != 0
    assert delta.gcd(delta.derivative()) == 1
    assert replacement.degree() == 39
    gcd, left, right = delta.xgcd(replacement)
    assert gcd == 1 and left * delta + right * replacement == 1
    print("Exceptional radial line: 42 distinct roots, two minors with exact Bezout identity")


def normalize_radial_power(polynomial, odd, ring):
    terms = {}
    for (radial_degree, *others), value in polynomial.dict().items():
        assert radial_degree % 2 == odd
        terms[(radial_degree // 2, *others)] = value
    return ring(terms)


def solve_pivot(polynomial, variable):
    pivot = polynomial.coefficient({variable: 1})
    assert polynomial.degree(variable) == 1 and pivot and pivot.is_constant()
    solution = -(polynomial - pivot * variable) / QQ(pivot)
    assert variable not in solution.variables()
    return solution


def reduced_residual_rows(odd_degree):
    names = tuple(f"l{j}" for j in range(1, odd_degree + 1)) + ("v",)
    phase_ring = PolynomialRing(QQ, names=("u", *names))
    u, *others = phase_ring.gens()
    phase = phase_terms((u, *(u * value for value in others[:-1])), others[-1])

    ring = PolynomialRing(QQ, names=("r", *names), order="degrevlex")
    subring = PolynomialRing(QQ, names=("r", "l", "v"), order="degrevlex")
    r, l, v = subring.gens()
    reduction = ring.hom([r, l, *([0] * (odd_degree - 1)), v], subring)
    poles = []
    odd_coefficients = {}
    residuals = []
    columns = [(1, 0), (3, 0), (5, 0), (0, 1), (1, 1)]
    rows = reconstruction_rows(phase, columns, odd_degree + 5, phase_ring)
    for k, (a1, a3, a5, c0, c1) in enumerate(rows):
        if k == 0:
            assert a1 == a3 == a5 == c0 == c1 == 0
            continue
        coherent = normalize_radial_power(a1 - 2*a3 + a5, 0, ring)
        nuisance = [normalize_radial_power(column, 1, ring) for column in (c0, c1)]
        residual = coherent + sum(column * entry for column, entry in zip(nuisance, poles))
        if k <= 2:
            pivot = nuisance[k-1]
            assert pivot and pivot.is_constant() and not any(nuisance[k:])
            poles.append(-residual / QQ(pivot))
            continue

        residual = ring(residual.subs(odd_coefficients))
        if k <= odd_degree + 1:
            variable = ring.gen(k-1)
            odd_coefficients[variable] = solve_pivot(residual, variable)
            continue

        assert not set(residual.variables()) & odd_coefficients.keys()
        reduced = reduction(residual)
        assert reduced
        integral = reduced * reduced.denominator()
        residuals.append(integral / QQ(integral.content()))
    return residuals


def verify_odd_cubic_family():
    residuals = reduced_residual_rows(3)
    assert [(f.degree(), len(f.monomials())) for f in residuals] == [
        (4, 13), (5, 16), (6, 27), (7, 33),
    ]
    ring = residuals[0].parent()
    multipliers = ring.one().lift(ring.ideal(residuals))
    assert sum(g * f for g, f in zip(multipliers, residuals, strict=True)) == 1
    print("Odd cubic polynomial and even radial term: rows 5 through 8 certify the unit ideal")


def verify_odd_quartic_family():
    residuals = reduced_residual_rows(4)
    assert [(f.degree(), len(f.monomials())) for f in residuals] == [
        (5, 19), (6, 26), (7, 37), (8, 48),
    ]
    ring = residuals[0].parent()
    certificate = Path(__file__).with_name("cm_newton_conductor_18_quartic_bezout.json")
    multipliers = []
    for polynomial in json.loads(certificate.read_text(encoding="utf-8")):
        terms = {(i, j, k): QQ(numerator) / QQ(denominator)
                 for i, j, k, numerator, denominator in polynomial}
        assert len(terms) == len(polynomial)
        multipliers.append(ring(terms))
    assert sum(g * f for g, f in zip(multipliers, residuals, strict=True)) == 1
    print("Odd quartic polynomial and even radial term: rows 6 through 9 certify the unit ideal")


if __name__ == "__main__":
    verify_exceptional_radial_line()
    verify_odd_cubic_family()
    verify_odd_quartic_family()

