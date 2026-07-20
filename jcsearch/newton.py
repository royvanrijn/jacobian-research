"""Newton-edge support utilities for the translated two-divisor chart.

For x=uv and y=u-(uv)^(-1), an ordinary monomial has the exact expansion

    x^a y^b = sum_k (-1)^k binomial(b,k)
                         u^(a+b-2k) v^(a-k).

The routines below use this identity to turn a proposed polynomial Newton
support into the smallest obvious Laurent ansatz containing it.  Cancellation
is still recomputed independently; the expansion is not assumed to describe
the whole polynomializable kernel.
"""
from __future__ import annotations

from math import gcd
import sympy as sp

from .translated import u, v


def monomial_laurent_terms(a: int, b: int):
    """Return {(u exponent, v exponent): coefficient} for x^a y^b."""
    if a < 0 or b < 0:
        raise ValueError("ordinary polynomial exponents must be nonnegative")
    return {
        (a + b - 2*k, a - k): (-1)**k * sp.binomial(b, k)
        for k in range(b + 1)
    }


def laurent_support(polynomial_support):
    """Union of exact Laurent supports of ordinary monomials."""
    result = set()
    for a, b in polynomial_support:
        result.update(monomial_laurent_terms(int(a), int(b)))
    return tuple(sorted(result))


def monomial_expansion(a: int, b: int):
    return sp.expand(sum(c*u**i*v**j for (i, j), c in
                         monomial_laurent_terms(a, b).items()))


def base_support(degree: int = 2):
    """All ordinary monomials of total degree at most ``degree``."""
    return {(a, b) for a in range(degree + 1)
            for b in range(degree + 1 - a)}


def weighted_edge(weight, top: int, width: int = 0):
    """Nonnegative lattice points in top-width <= r*a+s*b <= top."""
    r, s = map(int, weight)
    if r <= 0 or s <= 0 or gcd(r, s) != 1:
        raise ValueError("weight must be a positive primitive pair")
    if top < 0 or width < 0:
        raise ValueError("top and width must be nonnegative")
    points = set()
    for a in range(top // r + 1):
        for b in range(top // s + 1):
            value = r*a + s*b
            if top - width <= value <= top:
                points.add((a, b))
    return points


def total_degree_edge(top: int, width: int = 0):
    return weighted_edge((1, 1), top, width)


def leading_bracket_exponents(support_a, support_b):
    """Possible exponents in [A,B], ignoring coefficient cancellations."""
    possible = set()
    for i, j in support_a:
        for k, ell in support_b:
            if i*ell - j*k:
                possible.add((i + k - 1, j + ell - 1))
    return possible


def convex_lattice_points(vertices):
    """All integer points in/on a convex polygon with integral vertices."""
    vertices = tuple((int(a), int(b)) for a, b in vertices)
    if len(vertices) < 3:
        raise ValueError("at least three polygon vertices are required")
    twice_area = sum(a*d - b*c for (a, b), (c, d) in
                     zip(vertices, vertices[1:] + vertices[:1]))
    if twice_area == 0:
        raise ValueError("degenerate polygon")
    sign = 1 if twice_area > 0 else -1
    points = set()
    min_a, max_a = min(a for a, _ in vertices), max(a for a, _ in vertices)
    min_b, max_b = min(b for _, b in vertices), max(b for _, b in vertices)
    for a in range(min_a, max_a + 1):
        for b in range(min_b, max_b + 1):
            inside = True
            for (c, d), (e, f) in zip(vertices, vertices[1:] + vertices[:1]):
                cross = (e - c)*(b - d) - (f - d)*(a - c)
                if sign*cross < 0:
                    inside = False
                    break
            if inside:
                points.add((a, b))
    return tuple(sorted(points))


def reduced_72_108_polygons():
    """The two reduced `(8,28)` cases in Proposition 4.3 of arXiv:2204.14178."""
    return {
        "case_with_vertical_vertices": {
            "P_vertices": ((0, 0), (1, 0), (8, 14), (8, 16), (0, 8)),
            "Q_vertices": ((0, 0), (2, 1), (12, 21), (12, 24), (0, 12)),
        },
        "case_without_vertical_vertices": {
            "P_vertices": ((0, 0), (1, 0), (8, 14), (8, 16)),
            "Q_vertices": ((0, 0), (2, 1), (12, 21), (12, 24)),
        },
    }


def stage_d_families():
    """Small Newton/Puiseux-directed prototypes, not the full (72,108) chain."""
    low = base_support(2)
    common = {
        "total_ratio_2_3": (
            low | total_degree_edge(6),
            low | total_degree_edge(9),
            "total-degree edges 6 and 9"
        ),
        "weighted_3_2_ratio_2_3": (
            low | weighted_edge((3, 2), 12, 1),
            low | weighted_edge((3, 2), 18, 1),
            "thin (3,2)-weighted bands at 12 and 18"
        ),
        "weighted_2_3_ratio_2_3": (
            low | weighted_edge((2, 3), 12, 1),
            low | weighted_edge((2, 3), 18, 1),
            "thin (2,3)-weighted bands at 12 and 18"
        ),
        "corner_2_7_square_cube": (
            low | {(4, 14)},
            low | {(6, 21)},
            "square/cube of the primitive (2,7) corner"
        ),
        "corner_7_2_square_cube": (
            low | {(14, 4)},
            low | {(21, 6)},
            "coordinate-reversed square/cube of the primitive (2,7) corner"
        ),
    }
    return {name: {
        "P_polynomial_support": tuple(sorted(P)),
        "Q_polynomial_support": tuple(sorted(Q)),
        "description": description,
    } for name, (P, Q, description) in common.items()}
