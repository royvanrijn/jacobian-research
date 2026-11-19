#!/usr/bin/env python3
"""Exact finite-field checks for the three-height quintic classification."""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
from functools import cache
from hashlib import sha256
from itertools import product
import json
from math import comb, factorial, isqrt
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory


PRIME = 2731
SMALL_ORDERS = tuple(range(2, 19, 2))
# Entries give the numerator's x-degree and the t-degree before removing t.
# The two entries in grade zero include the separate factor 1-t^2.
HEIGHTS = {
    "a": (0, 1), "b": (0, 3),
    "c": (1, 1), "d": (1, 3), "e": (1, 5),
    "f": (2, 1), "g": (2, 3),
    "h": (3, 1), "i": (3, 3), "j": (4, 1), "k": (5, 1),
}
HIGH_ORDERS = (546, 630, 650)
HERE = Path(__file__).resolve().parent


def peel_orders():
    """Check that each selected height isolates the asserted single face."""
    assert all(PRIME % d for d in range(2, isqrt(PRIME) + 1))
    support = {name: (degree, height + 2 * (degree == 0))
               for name, (degree, height) in HEIGHTS.items()}
    for name, order in zip(("e", "i", "k"), HIGH_ORDERS):
        assert order % 2 == 0 and 0 < order < PRIME
        assert order < PRIME - 1 <= 5 * order < 3 * PRIME - 1
        degree, height = support[name]
        if name == "e":
            assert degree == 1 and height == 5
            assert all(d == 0 or h < height or v == name
                       for v, (d, h) in support.items())
            face_height = height
        else:
            slope = Fraction(height - 5, degree)
            assert support["b"] == (0, 5)
            assert all(h < 5 + slope * d
                       for v, (d, h) in support.items() if v not in ("b", name))
            assert order % degree == 0
            assert comb(order, order // degree) % PRIME
            face_height = 5 + slope
        assert order * face_height + 1 == PRIME
        del support[name]


@cache
def integral(order, t_square_power, endpoint_power):
    """Integral of t^(order+2*t_square_power)*(1-t^2)^endpoint_power."""
    def inverse(i):
        return pow(i, -1, PRIME)

    value = inverse(order + 2 * t_square_power + 1)
    for step in range(1, endpoint_power + 1):
        value = value * (2 * step) * inverse(order + 2 * t_square_power + 2 * step + 1)
        value %= PRIME
    return value


def compositions(total, count):
    if count == 0:
        if total == 0:
            yield ()
    elif count == 1:
        yield (total,)
    else:
        for first in range(total + 1):
            for tail in compositions(total - first, count - 1):
                yield (first, *tail)


def positive_grades(remaining, grades, degree=2):
    if degree == len(grades):
        yield ()
    else:
        limit = remaining // degree if grades[degree] else 0
        for count in range(limit + 1):
            for rest in positive_grades(remaining - degree * count, grades, degree + 1):
                yield (count, *rest)


def chart_moments(orders, zero, one):
    variables = tuple(name for name in HEIGHTS if name not in zero | one)
    grades = [[] for _ in range(6)]
    for name, (degree, height) in HEIGHTS.items():
        if name not in zero:
            grades[degree].append((name, (height - 1) // 2))

    equations = []
    for order in orders:
        polynomial = defaultdict(int)
        for high_counts in positive_grades(order, grades):
            # The order factors have total x-degree order in the numerator.
            linear_count = order - sum(d * k for d, k in enumerate(high_counts, 2))
            negative_count = order - linear_count - sum(high_counts)
            counts = (negative_count, linear_count, *high_counts)
            options = [compositions(count, len(terms))
                       for count, terms in zip(counts, grades)]
            for powers_by_grade in product(*options):
                powers = {}
                t_square_power = 0
                denominator = 1
                for terms, powers_here in zip(grades, powers_by_grade):
                    for (name, height), power in zip(terms, powers_here):
                        powers[name] = power
                        t_square_power += height * power
                        denominator *= factorial(power)
                exponent = tuple(powers[name] for name in variables)
                coefficient = (factorial(order) // denominator) * integral(
                    order, t_square_power, negative_count
                )
                polynomial[exponent] = (polynomial[exponent] + coefficient) % PRIME
        equations.append({monomial: value for monomial, value in polynomial.items() if value})
    return variables, equations


def literal_moment(order, values):
    """Independent numeric check using two-variable polynomial multiplication."""
    terms = defaultdict(int)
    for name, (degree, height) in HEIGHTS.items():
        t_square_power = (height - 1) // 2
        terms[degree, t_square_power] += values[name]
        if degree == 0:
            terms[degree, t_square_power + 1] -= values[name]
    powers = {(0, 0): 1}
    for _ in range(order):
        next_powers = defaultdict(int)
        for (x_old, z_old), v_old in powers.items():
            for (x, z), value in terms.items():
                if x_old + x <= order:
                    key = (x_old + x, z_old + z)
                    next_powers[key] = (next_powers[key] + v_old * value) % PRIME
        powers = next_powers
    return sum(
        value * pow(order + 2 * z + 1, -1, PRIME)
        for (x, z), value in powers.items() if x == order
    ) % PRIME


def replay_chart(msolve, workspace, label, orders, zero, one):
    variables, equations = chart_moments(orders, zero, one)
    for multiplier in (2, 11):
        values = {
            name: 0 if name in zero else 1 if name in one else multiplier * (i + 2)
            for i, name in enumerate(HEIGHTS)
        }
        for order, polynomial in zip(orders, equations):
            value = 0
            for monomial, coefficient in polynomial.items():
                for name, power in zip(variables, monomial):
                    coefficient *= pow(values[name], power, PRIME)
                value += coefficient
            assert value % PRIME == literal_moment(order, values), (label, order)

    if any(set(eq) == {(0,) * len(variables)} for eq in equations):
        return {"chart": label, "orders": orders, "unit": "nonzero constant"}

    def render(eq):
        terms = []
        for monomial, value in eq.items():
            factors = [str(value)]
            factors.extend(f"{name}^{power}" for name, power in zip(variables, monomial) if power)
            terms.append("*".join(factors))
        return "+".join(terms)

    source = workspace / f"{label}.ms"
    output = workspace / f"{label}.out"
    source.write_text(",".join(variables) + f"\n{PRIME}\n"
                      + ",\n".join(render(eq) for eq in equations if eq) + "\n")
    run([msolve, "-g", "2", "-l", "2", "-t", "2", "-f", str(source), "-o", str(output)],
        capture_output=True, text=True, check=True)
    basis = [line.strip() for line in output.read_text().splitlines()
             if line.strip() and not line.startswith("#")]
    assert basis == ["[1]:"], (label, basis[:3])
    return {"chart": label, "orders": orders, "unit": "exact Groebner basis [1]",
            "input_sha256": sha256(source.read_bytes()).hexdigest()}


def fixed_charts(msolve):
    receipts = []
    with TemporaryDirectory() as temporary:
        workspace = Path(temporary)

        def check_charts(label, initial_zero, initial_one, leader_limits):
            zero = set(initial_zero)
            for leader, bound in leader_limits:
                orders = tuple(m for m in SMALL_ORDERS if m <= bound)
                receipts.append(replay_chart(msolve, workspace, label + leader, orders,
                                             zero, set(initial_one) | {leader}))
                zero.add(leader)

        check_charts("zero-", "abfghijk", "", (("c", 6), ("d", 4), ("e", 2)))
        check_charts("a-", "be", "a", tuple(zip("cdfghijk", (16, 14, 12, 2, 12, 6, 4, 10))))
        check_charts("b-", "eik", "b", tuple(zip("cdfghj", (18, 12, 10, 8, 12, 8))))
    return receipts


def dormant_example():
    """Check the two univariate integrals defining the first nonzero moment."""
    def integral_polynomial(r):
        return [comb(8 * r, j) * integral(10 * r, j, 8 * r) % PRIME
                for j in range(8 * r + 1)]

    def remainder(a, b):
        a = a[:]
        while len(a) >= len(b):
            shift = len(a) - len(b)
            factor = a[-1] * pow(b[-1], -1, PRIME)
            for j, value in enumerate(b):
                a[shift + j] = (a[shift + j] - factor * value) % PRIME
            while a and a[-1] == 0:
                a.pop()
        return a

    a, b = integral_polynomial(1), integral_polynomial(2)
    assert a[-1] and b[-1]
    while b:
        a, b = b, remainder(a, b)
    assert len(a) == 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--msolve", default="msolve",
                        help="path to the msolve executable (checked with 0.6.5)")
    parser.add_argument("--write-receipt", type=Path)
    args = parser.parse_args()
    peel_orders()
    dormant_example()
    charts = fixed_charts(args.msolve)
    proof = HERE / "THREE_HEIGHT_QUINTIC_IMAGE_SAFETY.md"
    receipt = {"schema_version": 1, "status": "passed", "prime": PRIME,
               "isolating_orders": HIGH_ORDERS, "charts": charts,
               "first_nineteen_example": "orders 10 and 20 have coprime defining factors",
               "proof_sha256": sha256(proof.read_bytes()).hexdigest(),
               "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
               "exact_linear_algebra": "msolve -g 2 -l 2; internal replay"}
    if args.write_receipt:
        args.write_receipt.write_text(json.dumps(receipt, indent=2) + "\n")
    print(f"PASS: {len(charts)} exact charts over F_{PRIME}, "
          "three isolated coefficients, and example gcd")


if __name__ == "__main__":
    main()
