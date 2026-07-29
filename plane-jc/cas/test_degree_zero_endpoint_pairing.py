#!/usr/bin/env python3
"""Regression for the degree-zero endpoint-pairing insufficiency theorem.

The proof in JC2_GLOBAL_COX_PACKET_ATTACK.md is uniform in n.  This script
checks representative members of the explicit family and the bounded
monomial normal forms of the graded source bridge; the bounded loop is a
regression, not the proof of infinitude.
"""
from __future__ import annotations

import hashlib
import itertools
from pathlib import Path

import sympy as sp

from endpoint_valuation_compiler import (
    endpoint_initials,
    monotone_reduce,
    monotone_polynomial_reduce,
    pole_matrix,
    pole_height,
    reducing_polynomial_shears,
    reducing_shears,
)


x, y, t = sp.symbols("x y t")
a, s = sp.symbols("a s")

COMPILER_PATH = Path(__file__).with_name("endpoint_valuation_compiler.py")
COMPILER_SHA256 = (
    "2b94709c916c287adad99597a8f4858f4a49bfb6566f73c1595687c2af52f597"
)
assert hashlib.sha256(COMPILER_PATH.read_bytes()).hexdigest() == COMPILER_SHA256


def laurent_poles(expression: sp.Expr) -> tuple[int, int]:
    """Return pole orders at t=0 and t=infinity for a Laurent polynomial."""
    initials = endpoint_initials(expression, t)
    return initials[0].pole, initials[1].pole


def reduce_bridge_monomial(s_power: int, a_power: int) -> tuple[int, int, int]:
    """Reduce s^i a^j using a*s^2=h.

    Return (remaining s exponent, remaining a exponent, h exponent).
    """
    h_power = min(s_power // 2, a_power)
    return s_power - 2 * h_power, a_power - h_power, h_power


# Every homogeneous bridge piece has one of the three normal-form generators
# s^d, a^m, or a^(m+1)*s.  Check all monomials in a generous finite box.
for s_power in range(25):
    for a_power in range(13):
        degree = s_power - 2 * a_power
        reduced_s, reduced_a, h_power = reduce_bridge_monomial(
            s_power, a_power
        )
        assert degree == reduced_s - 2 * reduced_a
        assert reduced_s < 2 or reduced_a == 0
        if degree >= 0:
            assert (reduced_s, reduced_a) == (degree, 0)
        elif degree % 2 == 0:
            assert (reduced_s, reduced_a) == (0, -degree // 2)
        else:
            assert (reduced_s, reduced_a) == (1, (1 - degree) // 2)
        assert h_power >= 0

# The odd-square generator identity in the bridge.
h_symbol = sp.symbols("h")
assert sp.expand((a * s) ** 2).subs(a * s**2, h_symbol) == a * h_symbol

cusp = y**2 - x**3

for n in range(1, 21):
    connector = x * y - y ** (n + 1) - 1
    packet_divisor = sp.expand(cusp * connector)

    # The connector component is smooth, the two components are coprime,
    # and their product is squarefree.
    connector_jacobian = sp.groebner(
        [connector, sp.diff(connector, x), sp.diff(connector, y)],
        x,
        y,
        domain=sp.QQ,
    )
    assert connector_jacobian.contains(sp.Integer(1))
    assert sp.gcd(cusp, connector) == 1
    assert sp.gcd(
        packet_divisor,
        sp.gcd(
            sp.diff(packet_divisor, x),
            sp.diff(packet_divisor, y),
        ),
    ) == 1

    # Normalizations:
    #   cusp:     (x,y)=(t^2,t^3);
    #   connector: (x,y)=(t^n+t^-1,t).
    assert sp.expand(cusp.subs({x: t**2, y: t**3})) == 0
    assert sp.cancel(
        connector.subs({x: t**n + t ** (-1), y: t})
    ) == 0

    # Pole pairs (pole(x), pole(y)) are (2,3) at the cusp endpoint and
    # (1,0), (n,1) at the two connector endpoints.  The connector matrix
    # is unimodular for every n, so a saturation/unimodularity prefilter
    # retains the entire infinite family.
    cusp_pole_pair = (2, 3)
    connector_pole_pairs = ((1, 0), (n, 1))
    determinant = (
        connector_pole_pairs[0][0] * connector_pole_pairs[1][1]
        - connector_pole_pairs[0][1] * connector_pole_pairs[1][0]
    )
    assert cusp_pole_pair == (2, 3)
    assert determinant == 1

# Raw pole pairs of polynomial coordinate generators are not invariant.
# The triangular coordinate change (X,Y)=(x,y+x^m) has Jacobian one and
# sends the connector parametrization to
# X=t^n+t^-1, Y=t+(t^n+t^-1)^m.
for n in range(1, 9):
    connector_x = t**n + t ** (-1)
    assert laurent_poles(connector_x) == (1, n)

    # The inverse shear u=x-y^n, v=y is an automorphism-minimal connector
    # presentation: u=t^-1 and v=t have the two axis pole pairs.
    reduced_connector_x = sp.expand(connector_x - t**n)
    assert reduced_connector_x == t ** (-1)
    assert laurent_poles(reduced_connector_x) == (1, 0)
    assert laurent_poles(t) == (0, 1)
    connector_reduction = monotone_reduce(connector_x, t, t)
    assert pole_matrix(
        connector_reduction.first,
        connector_reduction.second,
        t,
    ) == ((1, 0), (0, 1))
    assert connector_reduction.height == 2
    assert monotone_polynomial_reduce(connector_x, t, t).height == 2

    for m in range(2, 9):
        connector_y_sheared = t + connector_x**m
        assert laurent_poles(connector_y_sheared) == (m, m * n)
        assert sp.det(sp.Matrix([[1, 0], [m * x ** (m - 1), 1]])) == 1

        # The inverse triangular change recovers the original coordinate.
        recovered_y = sp.expand(
            connector_y_sheared - connector_x**m
        )
        assert recovered_y == t
        sheared_reduction = monotone_reduce(
            connector_x, connector_y_sheared, t
        )
        assert sheared_reduction.height == 2
        assert len(sheared_reduction.steps) == 2
        assert monotone_polynomial_reduce(
            connector_x, connector_y_sheared, t
        ).height == 2

        cusp_x = t**2
        cusp_y_sheared = t**3 + cusp_x**m
        assert laurent_poles(cusp_x) == (0, 2)
        assert laurent_poles(cusp_y_sheared) == (0, 2 * m)

# Equal valuation matrices do not determine whether a triangular reduction
# exists: the leading residues at all marked endpoints must be retained.
residue_x = t + t ** (-1)
for residue_parameter in (1, 2, 3):
    residue_y = (
        t**2
        + residue_parameter * t ** (-2)
        + t
    )
    assert laurent_poles(residue_x) == (1, 1)
    assert laurent_poles(residue_y) == (2, 2)

    for shear_coefficient in (1, 2, 3):
        remainder = sp.expand(
            residue_y - shear_coefficient * residue_x**2
        )
        zero_pole, infinity_pole = laurent_poles(remainder)
        assert (infinity_pole < 2) == (shear_coefficient == 1)
        assert (zero_pole < 2) == (
            shear_coefficient == residue_parameter
        )

    # Simultaneous quadratic cancellation occurs exactly when the two
    # endpoint residue ratios agree.
    simultaneous = [
        coefficient
        for coefficient in (1, 2, 3)
        if all(
            pole < 2
            for pole in laurent_poles(
                sp.expand(residue_y - coefficient * residue_x**2)
            )
        )
    ]
    assert simultaneous == ([1] if residue_parameter == 1 else [])

    initial_moves = reducing_shears(residue_x, residue_y, t)
    assert {
        move[2].coefficient for move in initial_moves
    } == {sp.Integer(1), sp.Integer(residue_parameter)}
    residue_reduction = monotone_reduce(residue_x, residue_y, t)
    assert residue_reduction.height == (
        2 if residue_parameter == 1 else 4
    )

# A complete polynomial shear can lower height even when its leading
# monomial cancellation is only height-neutral.  The compiler must retain
# that prefix and continue through strictly descending polynomial degrees.
polynomial_base = t ** (-1) + t**2
polynomial_target = (
    polynomial_base**3
    - t**6
    + t**5
    + polynomial_base**2
)
assert pole_matrix(polynomial_base, polynomial_target, t) == (
    (1, 3),
    (2, 5),
)
assert pole_height(polynomial_base, polynomial_target, t) == 11
assert reducing_shears(polynomial_base, polynomial_target, t) == ()

polynomial_moves = reducing_polynomial_shears(
    polynomial_base, polynomial_target, t
)
assert polynomial_moves
best_polynomial_move = polynomial_moves[0]
assert best_polynomial_move[2].target == 1
assert best_polynomial_move[2].terms == (
    (3, sp.Integer(1)),
    (2, sp.Integer(1)),
)
assert best_polynomial_move[2].height_before == 11
assert best_polynomial_move[2].height_after == 9
assert sp.expand(best_polynomial_move[1]) == -t**6 + t**5
assert monotone_polynomial_reduce(
    polynomial_base, polynomial_target, t
).height == 9

# If a shear makes one restricted coordinate zero, an endpoint-invisible
# constant translation restores a nonzero coordinate of pole height zero.
constant_base = t + t ** (-1)
constant_target = constant_base**2
constant_reduction = monotone_polynomial_reduce(
    constant_base, constant_target, t
)
assert constant_reduction.second == 1
assert constant_reduction.height == 2

# Bounded evidence for the remaining alternating peak theorem.  Enumerate
# all 49 positive-support pairs on exponents {-1,0,1}, then all alternating
# two-step monomial shears of degrees 1,2 and coefficients +/-1.  Every
# globally lowering path with a nondecreasing first step already has a
# lowering complete polynomial shear at its initial pair.
peak_polynomials = [
    sp.Add(
        *[
            t**exponent
            for index, exponent in enumerate((-1, 0, 1))
            if mask >> index & 1
        ]
    )
    for mask in range(1, 8)
]
peak_candidates = 0
terminal_peak_candidates = 0
for peak_first, peak_second in itertools.product(
    peak_polynomials, repeat=2
):
    initial_height = pole_height(peak_first, peak_second, t)
    for first_target, degree, coefficient in itertools.product(
        (0, 1), (1, 2), (-1, 1)
    ):
        if first_target == 0:
            middle_first = sp.expand(
                peak_first + coefficient * peak_second**degree
            )
            middle_second = peak_second
        else:
            middle_first = peak_first
            middle_second = sp.expand(
                peak_second + coefficient * peak_first**degree
            )
        if middle_first == 0 or middle_second == 0:
            continue
        if pole_height(middle_first, middle_second, t) < initial_height:
            continue

        for next_degree, next_coefficient in itertools.product(
            (1, 2), (-1, 1)
        ):
            if first_target == 0:
                final_first = middle_first
                final_second = sp.expand(
                    middle_second
                    + next_coefficient * middle_first**next_degree
                )
            else:
                final_first = sp.expand(
                    middle_first
                    + next_coefficient * middle_second**next_degree
                )
                final_second = middle_second
            if final_first == 0 or final_second == 0:
                continue
            if pole_height(final_first, final_second, t) >= initial_height:
                continue

            peak_candidates += 1
            if not reducing_polynomial_shears(
                peak_first, peak_second, t
            ):
                terminal_peak_candidates += 1

assert peak_candidates == 16
assert terminal_peak_candidates == 0

print(
    "PASS: the four degree-zero filters retain an infinite family, and "
    "coordinate reduction requires endpoint valuation residues"
)
