#!/usr/bin/env python3
"""Exact finite-degree audit of the real-sheet spectrum theorem."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

w, t = sp.symbols("w t")


def extra_roots(degree: int) -> tuple[int, ...]:
    """The extra integral roots in the all-degree rational-fiber seed."""
    k = degree // 2
    if degree % 2:
        return (2,) + tuple(
            root for j in range(3, k + 2) for root in (j, 1 - j)
        )
    return (3, 4) + tuple(
        root for j in range(5, k + 3) for root in (j, 1 - j)
    )


def seed_data(degree: int) -> tuple[sp.Poly, sp.Rational]:
    """Return H_N and the rational slope making its pencil G_N/lambda_N."""
    roots = (0, -1) + extra_roots(degree)
    assert len(roots) == degree and len(set(roots)) == degree
    G = sp.Poly(sp.prod(w - root for root in roots), w, domain=sp.QQ)
    g0 = G.diff().eval(0)
    lam = g0 - G.diff().eval(1)
    assert lam != 0
    H = sp.Poly((G.as_expr() - g0 * w) / lam, w, domain=sp.QQ)
    slope = -sp.Rational(g0, lam)
    assert sp.Poly(H.as_expr() - slope * w, w) == sp.Poly(
        G.as_expr() / lam, w
    )
    return H, slope


def real_root_count(poly: sp.Poly) -> int:
    """Count real roots exactly; callers separately certify squarefreeness."""
    return int(poly.count_roots(-sp.oo, sp.oo))


def generic_rational_slope(
    H: sp.Poly, slope0: sp.Rational, degree: int
) -> tuple[sp.Rational, sp.Poly, sp.Poly]:
    """Find a nearby rational slope with distinct real critical values."""
    for exponent in range(3, 10):
        denominator = 10**exponent
        for numerator in (1, -1):
            slope = slope0 + sp.Rational(numerator, denominator)
            P = sp.Poly(H.as_expr() - slope * w, w, domain=sp.QQ)
            if P.gcd(P.diff()).degree() != 0:
                continue
            if real_root_count(P) != degree:
                continue
            derivative = P.diff()
            if derivative.gcd(derivative.diff()).degree() != 0:
                continue
            if real_root_count(derivative) != degree - 1:
                continue
            critical_discriminant = sp.Poly(
                sp.discriminant(P.as_expr() + t, w), t, domain=sp.QQ
            )
            if critical_discriminant.degree() != degree - 1:
                continue
            if critical_discriminant.gcd(critical_discriminant.diff()).degree() != 0:
                continue
            if real_root_count(critical_discriminant) != degree - 1:
                continue
            return slope, P, critical_discriminant
    raise AssertionError(f"no audited generic rational slope in degree {degree}")


def positive_ray_thresholds(
    critical_discriminant: sp.Poly, direction: int
) -> list[tuple[sp.Rational, sp.Rational]]:
    """Isolate q>0 for which Disc(P + direction*q)=0."""
    raw_intervals = critical_discriminant.intervals(
        eps=sp.Rational(1, 10**12)
    )
    thresholds: list[tuple[sp.Rational, sp.Rational]] = []
    for (left, right), multiplicity in raw_intervals:
        assert multiplicity == 1
        left = sp.Rational(left)
        right = sp.Rational(right)
        assert not (left <= 0 <= right)
        if direction == 1 and left > 0:
            thresholds.append((left, right))
        elif direction == -1 and right < 0:
            thresholds.append((-right, -left))
    thresholds.sort()
    for first, second in zip(thresholds, thresholds[1:]):
        assert first[1] < second[0]
    return thresholds


def simple_rational_between(
    lower: sp.Rational, upper: sp.Rational
) -> sp.Rational:
    """Choose a low-denominator rational strictly inside an exact interval."""
    assert lower < upper
    for denominator in range(1, 100_001):
        numerator = int(sp.floor(lower * denominator)) + 1
        candidate = sp.Rational(numerator, denominator)
        if candidate < upper:
            return candidate
    raise AssertionError("failed to find a simple rational interval sample")


def audit_degree(degree: int) -> tuple[sp.Rational, dict[int, sp.Rational]]:
    """Return and exactly certify rational t-witnesses for one degree."""
    H, slope0 = seed_data(degree)
    slope, P, critical_discriminant = generic_rational_slope(
        H, slope0, degree
    )
    direction = 1 if P.LC() > 0 else -1
    thresholds = positive_ray_thresholds(critical_discriminant, direction)
    assert len(thresholds) == degree // 2

    witnesses: dict[int, sp.Rational] = {degree: sp.Rational(0)}
    for crossed in range(1, len(thresholds) + 1):
        if crossed < len(thresholds):
            q = simple_rational_between(
                thresholds[crossed - 1][1], thresholds[crossed][0]
            )
        else:
            q = sp.floor(thresholds[-1][1]) + 1
        target_t = sp.Rational(direction) * q
        fiber = sp.Poly(P.as_expr() + target_t, w, domain=sp.QQ)
        assert fiber.gcd(fiber.diff()).degree() == 0
        count = real_root_count(fiber)
        expected = degree - 2 * crossed
        assert count == expected
        witnesses[count] = target_t

    assert set(witnesses) == set(range(degree % 2, degree + 1, 2))
    return slope, witnesses


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--max-degree",
        type=int,
        default=12,
        help="largest audited degree (default: 12)",
    )
    parser.add_argument(
        "--show-witnesses",
        action="store_true",
        help="print the exact rational (s,t) witness for every count",
    )
    args = parser.parse_args()
    if args.max_degree < 3:
        parser.error("--max-degree must be at least 3")

    for degree in range(3, args.max_degree + 1):
        slope, witnesses = audit_degree(degree)
        if args.show_witnesses:
            rendered = ", ".join(
                f"{count}:(s={slope},t={witnesses[count]})"
                for count in sorted(witnesses, reverse=True)
            )
            print(f"N={degree}: {rendered}")

    print(
        "PASS: exact rational witnesses realize every parity-compatible "
        f"real-sheet count in degrees 3 through {args.max_degree}"
    )
    print("PASS: each audited vertical chain crosses distinct simple critical values")


if __name__ == "__main__":
    main()
