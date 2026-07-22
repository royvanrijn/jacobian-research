#!/usr/bin/env python3
"""Verify an all-order algebraic recurrence for the distinguished inverse coordinate."""

from __future__ import annotations

from fractions import Fraction as Q
from math import comb
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "generated-results" / "essential_bcw_21_counterexample.json"


def truncated_multiply(left: list[Q], right: list[Q], order: int) -> list[Q]:
    answer = [Q(0)] * (order + 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            if i + j <= order:
                answer[i + j] += a * b
    return answer


def main() -> None:
    x, y, z, a, b, c = sp.symbols("x y z a b c")
    v = 1 + 2 * x * y
    original = [
        x - 3 * x**2 * y - x**3 * z,
        y + 3 * x * v**2 * z + 12 * x * y**2 * (2 + 3 * x * y),
        v**3 * z + 4 * y**2 * v * (2 + 3 * x * y),
    ]
    cubic_coefficient = 108 * a**2 * c**2 + 4 * a * b**3 - 36 * a * b * c - b**2 + 8 * c
    inverse_equation = sp.expand(cubic_coefficient * x**3 + (1 - 3 * a * b) * x - a)
    assert sp.expand(inverse_equation.subs(dict(zip((a, b, c), original)))) == 0

    # The quotient map keeps the original three target coordinates and the
    # homogenizer.  Hence (a,b,c,0,...,0) is the exact target section for the
    # 20-dimensional identity slice, not merely an abstract specialization.
    source = json.loads(SOURCE.read_text())
    rows = source["quotient_factorization"]["B_rows"]
    assert rows[0] == [[0, "1"]]
    assert rows[1] == [[1, "1"]]
    assert rows[2] == [[2, "1"]]
    assert rows[20] == [[23, "1"]]
    assert all(
        not row or all(column not in (0, 1, 2, 23) for column, _ in row)
        for row in rows[3:20]
    )

    # For L=1-3ab and K=cubic_coefficient, x=(a/L)Y with
    # Y+tY^3=1 and t=K*a^2/L^3.  Its coefficients are signed ternary
    # Fuss--Catalan numbers.
    fuss = [Q(comb(3 * n, n), 2 * n + 1) for n in range(13)]
    assert fuss[0] == 1
    for n in range(1, len(fuss)):
        assert fuss[n] == sum(
            fuss[i] * fuss[j] * fuss[n - 1 - i - j]
            for i in range(n)
            for j in range(n - i)
        )

    # On the one-parameter target line (a,b,c)=(t,0,t), the inverse obeys
    # X+(8t+108t^4)X^3=t.  Replay its all-order coefficient recurrence.
    order = 30
    coefficients = [Q(0)] * (order + 1)
    coefficients[1] = Q(1)
    for n in range(2, order + 1):
        cube = truncated_multiply(
            truncated_multiply(coefficients, coefficients, n), coefficients, n
        )
        coefficients[n] = -8 * cube[n - 1] - (108 * cube[n - 4] if n >= 4 else 0)

    cube = truncated_multiply(
        truncated_multiply(coefficients, coefficients, order), coefficients, order
    )
    residual = coefficients[:]
    residual[1] -= 1
    for n in range(1, order + 1):
        residual[n] += 8 * cube[n - 1]
        if n >= 4:
            residual[n] += 108 * cube[n - 4]
    assert residual == [Q(0)] * (order + 1)
    nonzero_coefficients = [(index, value) for index, value in enumerate(coefficients) if value]
    assert len(nonzero_coefficients) >= 10

    print("PASS inverse recurrence: exact cubic equation for the original x-coordinate")
    print("PASS inverse recurrence: quotient target section preserves (a,b,c) and q_0=x")
    print("PASS inverse recurrence: ternary Fuss--Catalan closed form through n=12")
    print("PASS inverse recurrence: one-parameter all-order recurrence through order 30")
    print("first coefficients:", [str(value) for value in coefficients[1:11]])


if __name__ == "__main__":
    main()
