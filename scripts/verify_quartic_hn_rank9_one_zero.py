#!/usr/bin/env python3
"""Exact audit of the rank-nine one-zero-Gale quartic HN gate.

This script starts from the 16-parameter triple-root normal form recorded in
QUARTIC_HN_WARING_RIGIDITY.md.  It verifies only the finite coefficient
identities used in Theorem 7.9:

* the K_(3,3) channel has already been checked by the companion script;
* the second-trace equations used to close the p0 and p3 charts;
* the six channel-separation monomials;
* the three higher trace pivots needed when p10 is the first active channel;
* the final quartic equation forcing p13=0.

No bounded coefficient search is involved.  All calculations are exact over
QQ with SymPy.
"""

from __future__ import annotations

import sympy as sp

x1, x2, y1, y2 = sp.symbols("x1 x2 y1 y2")
X = (x1, x2, y1, y2)
p = sp.symbols("p0:16")

# Triple-root linear normal form, already after the exact 19-rank linear gate.
E = (
    3 * p[0] * y2**4
    + 3 * p[1] * y1 * y2**3
    + 3 * p[2] * x2**3 * y1
    + 3 * p[3] * x2**4
    + 3 * p[4] * x1 * y2**3
    + 3 * p[5] * x1 * y1 * y2**2
    - p[5] * x2 * y2**3
    + 3 * p[6] * x1 * x2**2 * y1
    - p[6] * x2**3 * y2
    + 3 * p[7] * x1 * x2**3
    + 3 * p[8] * x1**2 * y2**2
    + 3 * p[9] * x1**2 * y1 * y2
    - 3 * p[9] * x1 * x2 * y2**2
    + 3 * p[10] * x1**2 * x2 * y1
    - 3 * p[10] * x1 * x2**2 * y2
    + 3 * p[11] * x1**2 * x2**2
    + 3 * p[12] * x1**3 * y2
    + 3 * p[13] * x1**3 * y1
    - 9 * p[13] * x1**2 * x2 * y2
    + 3 * p[14] * x1**3 * x2
    + 3 * p[15] * x1**4
) / 3

Q = sp.Matrix(
    [
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [1, 0, 0, 0],
        [0, 1, 0, 0],
    ]
)
C = x1**3
M = Q * sp.hessian(E, X)
N = Q * sp.hessian(C, X)

# Normalize to the channel p1=p5=p9=0.  The other channel is exchanged by
# x2 <-> y2.
CHANNEL = {p[index]: 0 for index in (1, 5, 9)}
M_CHANNEL = sp.simplify(M.subs(CHANNEL))


def coefficient(expression: sp.Expr, monomial: tuple[int, int, int, int]) -> sp.Expr:
    return sp.factor(sp.Poly(sp.expand(expression), *X).coeff_monomial(monomial))


def assert_equal(label: str, actual: sp.Expr, expected: sp.Expr) -> None:
    if sp.expand(actual - expected) != 0:
        raise AssertionError((label, sp.factor(actual), sp.factor(expected)))
    print("QUARTIC_HN_WARING_ONE_ZERO_COEFFICIENT_PASS", label, sp.factor(actual))


trace2 = sp.expand(sp.trace(M_CHANNEL**2))

TRACE2_TARGETS = (
    ("p0p2", (0, 1, 1, 2), 144 * p[0] * p[2]),
    ("p0p6", (1, 0, 1, 2), 48 * p[0] * p[6]),
    ("p0p10", (1, 0, 0, 3), -48 * p[0] * p[10]),
    ("p0p7", (1, 1, 0, 2), 144 * p[0] * p[7]),
    (
        "p0p11",
        (2, 0, 0, 2),
        12 * (4 * p[0] * p[11] - p[10] * p[4]),
    ),
    (
        "p0p3",
        (0, 2, 0, 2),
        36 * (8 * p[0] * p[3] + p[2] * p[4]),
    ),
    (
        "p13-first",
        (4, 0, 0, 0),
        4 * (3 * p[10] * p[12] + 2 * p[11] * p[8] + 9 * p[13] ** 2),
    ),
    ("p2p4", (1, 1, 1, 1), 72 * p[2] * p[4]),
    (
        "p3p4",
        (1, 2, 0, 1),
        48 * (p[2] * p[8] + 3 * p[3] * p[4]),
    ),
    (
        "p6-square",
        (0, 4, 0, 0),
        -4 * (3 * p[10] * p[2] - p[6] ** 2),
    ),
    (
        "p10p8",
        (3, 0, 0, 1),
        8 * (p[10] * p[8] + 3 * p[11] * p[4]),
    ),
    (
        "p3p8",
        (2, 2, 0, 0),
        12
        * (
            p[10] ** 2
            + 3 * p[12] * p[2]
            - 2 * p[13] * p[6]
            + 4 * p[3] * p[8]
        ),
    ),
    ("p6p4", (2, 0, 1, 1), 24 * p[4] * p[6]),
    ("p6p8", (3, 0, 1, 0), 8 * p[6] * p[8]),
)

for label, monomial, expected in TRACE2_TARGETS:
    assert_equal(label, coefficient(trace2, monomial), expected)

# The three p10 channel pivots occur in the coefficient of tau in trace(M+tau*N)^4,
# which is 4 trace(M^3 N) by cyclicity.
trace4_linear = sp.expand(4 * sp.trace(M_CHANNEL**3 * N))
TRACE4_LINEAR_TARGETS = (
    ("p2-square-p8", (3, 4, 0, 0), 432 * p[2] ** 2 * p[8]),
    ("p10-square-p4", (6, 0, 0, 1), 144 * p[10] ** 2 * p[4]),
    ("p10-square-p8", (7, 0, 0, 0), 48 * p[10] ** 2 * p[8]),
)
for label, monomial, expected in TRACE4_LINEAR_TARGETS:
    assert_equal(label, coefficient(trace4_linear, monomial), expected)

# Deepest branch p0=p2=p3=p6=p10=0 in addition to the channel normalization.
DEEP = CHANNEL | {p[index]: 0 for index in (0, 2, 3, 6, 10)}
M_DEEP = sp.simplify(M.subs(DEEP))
trace2_deep = sp.expand(sp.trace(M_DEEP**2))
trace4_deep = sp.expand(sp.trace(M_DEEP**4))

DEEP_TRACE2_TARGETS = (
    (
        "deep-quadratic",
        (4, 0, 0, 0),
        4 * (2 * p[11] * p[8] + 9 * p[13] ** 2),
    ),
    ("deep-p7p8", (3, 1, 0, 0), 24 * p[7] * p[8]),
    ("deep-p11p4", (3, 0, 0, 1), 24 * p[11] * p[4]),
    ("deep-p4p7", (2, 1, 0, 1), 72 * p[4] * p[7]),
)
for label, monomial, expected in DEEP_TRACE2_TARGETS:
    assert_equal(label, coefficient(trace2_deep, monomial), expected)

assert_equal(
    "deep-quartic",
    coefficient(trace4_deep, (8, 0, 0, 0)),
    4
    * (
        8 * p[11] ** 2 * p[8] ** 2
        + 108 * p[11] * p[13] ** 2 * p[8]
        + 81 * p[13] ** 4
    ),
)

print("QUARTIC_HN_WARING_RANK9_ONE_ZERO_GATE_PASS")
