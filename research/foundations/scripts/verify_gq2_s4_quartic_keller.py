#!/usr/bin/env python3
"""Compile the named Roe--Turturean S_4 action into a Keller fiber.

The finite presentation has three marked S_4 quotient orbits over the fixed
tame frame, and Q_2 has three corresponding S_4 extensions.  The base
quadratic obstruction is one on x_0=(12)(34) and zero on the other two
choices.  The companion GP certificate computes the same normalized relative
Stiefel--Whitney bit on the polynomial side, naming the unique matching model.
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from jcsearch.keller_fiber import compile_polynomial_to_keller_fiber  # noqa: E402
from verify_finite_etale_keller_fibers import (  # noqa: E402
    check_scheme_reconstruction,
)
from verify_gq2_permutation_action import (  # noqa: E402
    CONVENTION,
    SCHEMA,
    conjugate,
    identity,
    multiply,
    verify_certificate,
)


sigma = [1, 0, 2, 3]
tau = [1, 2, 0, 3]
one = [0, 1, 2, 3]
all_wild_choices = [
    one,
    [1, 0, 3, 2],
    [2, 3, 0, 1],
    [3, 2, 1, 0],
]

reports = []
for x0 in all_wild_choices:
    certificate = {
        "schema": SCHEMA,
        "degree": 4,
        "conventions": CONVENTION,
        "generators": {
            "sigma": sigma,
            "tau": tau,
            "x0": x0,
            "x1": one,
        },
    }
    reports.append(verify_certificate(certificate, max_group_order=24))

assert all(report["admissible"] for report in reports)
assert [report["image_group_order"] for report in reports] == [6, 24, 24, 24]
assert [report["wild_normal_closure_order"] for report in reports] == [1, 4, 4, 4]
exact_image_reports = [
    report for report in reports if report["image_group_order"] == 24
]
assert len(exact_image_reports) == 3
assert all(report["orbit_sizes"] == [4] for report in exact_image_reports)
assert all(report["stabilizer_orders"] == [6] for report in exact_image_reports)

sigma_tuple = tuple(sigma)
tau_tuple = tuple(tau)
common_centralizer = [
    element
    for element in itertools.permutations(range(4))
    if conjugate(sigma_tuple, element) == sigma_tuple
    and conjugate(tau_tuple, element) == tau_tuple
]
assert common_centralizer == [tuple(one)]

# For the ramified simple V_4 module, Roe--Turturean (6.25) gives
#
#   Q_A^0(c) = q(c) + b_q(c, S^{-1}c).
#
# Here q is the unique S_3-invariant nonsingular quadratic form: it is one on
# every nonzero element of V_4.  The fixed double transposition is therefore
# the unique class with obstruction one.
one_tuple = identity(4)
wild_tuples = [tuple(choice) for choice in all_wild_choices[1:]]


def q_value(value: tuple[int, ...]) -> int:
    assert value in [tuple(choice) for choice in all_wild_choices]
    return int(value != one_tuple)


def polar(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    return q_value(multiply(left, right)) ^ q_value(left) ^ q_value(right)


candidate_obstruction_bits = [
    q_value(choice) ^ polar(choice, conjugate(choice, sigma_tuple))
    for choice in wild_tuples
]
assert candidate_obstruction_bits == [1, 0, 0]

# The companion GP certificate proves that LMFDB 2.1.4.8a1.2 is the unique
# quartic model with the same normalized relative Stiefel--Whitney bit one.
# This names the Appendix-D tuple x_0=(12)(34), without choosing an arbitrary
# bijection between the two remaining zero-obstruction classes.
T, S = sp.symbols("T S")
x, y, z = sp.symbols("x y z")
polynomial = T**4 + 4 * T**2 - 4 * T + 2
compilation = compile_polynomial_to_keller_fiber(
    polynomial,
    T,
    translation=1,
    inverse_variable=S,
    source_variables=(x, y, z),
)

assert compilation.seed == S**4 + 4 * S**3 + 10 * S**2 + 8 * S
assert compilation.target == (1, 0, sp.Rational(-3, 4))
assert compilation.inverse_polynomial == sp.expand(polynomial.subs(T, 1 + S))
assert compilation.coordinate_degrees == (7, 26, 24)
assert sp.factor(
    sp.Matrix(compilation.determinant_one_map).jacobian((x, y, z)).det()
) == 1
check_scheme_reconstruction(polynomial.subs(T, S), sp.Integer(1))

print("PASS: exactly three marked S_4 action orbits lie over the fixed tame frame")
print("PASS: candidate obstruction bits are [1,0,0], uniquely one at (12)(34)")
print("PASS: T^4+4T^2-4T+2 compiles at translation 1 to target (1,0,-3/4)")
print("PASS: the determinant-one map has coordinate degrees (7,26,24)")
print("PASS: its complete fiber is Spec(Q[T]/(T^4+4T^2-4T+2))")
print("MATCH: x_0=(12)(34) is the unique relative-SW2-one quartic action")
print("LIMIT: relative SW2 does not order the other two marked x_0 choices")
