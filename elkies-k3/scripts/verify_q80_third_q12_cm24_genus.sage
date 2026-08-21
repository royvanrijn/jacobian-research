#!/usr/bin/env sage
"""Verify generic fibers of the resolved CM24 third-q12 pencil have genus one.

This verifier specializes a generic new-base value, reduces at the split good
prime 73 (s=33, j=17), and asks Singular to normalize the residual cubic.  Two
independent new-base values give irreducible genus-one curves, certifying the
resolved D7 filtration and distinguishing it from the rejected unresolved
cusp jets.
"""

import argparse
from pathlib import Path

from sage.all import GF, PolynomialRing
from sage.interfaces.singular import singular


HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument("--new-base", type=int, default=7)
args = parser.parse_args()
load(str(HERE / "derive_q80_third_q12_cm24_pencil.sage"))

finite = GF(73)
finite_s = finite(33)
finite_j = finite(17)


def reduce_quadratic(value):
    value = K(value)
    return finite(value[0]) + finite_s * finite(value[1])


def reduce_relative(value):
    value = L(value)
    return reduce_quadratic(value[0]) + finite_j * reduce_quadratic(value[1])


new_base_value = L(args.new_base)
finite_plane = PolynomialRing(finite, names=("w", "x"))
finite_w, finite_x = finite_plane.gens()
finite_equation = finite_plane(0)
for x_degree, coefficient in enumerate(residual_cubic.list()):
    coefficient = new_old_base(coefficient)
    finite_coefficient = finite_plane(0)
    for w_degree, parameter_coefficient in enumerate(coefficient.list()):
        specialized = new_parameter_ring(parameter_coefficient)(new_base_value)
        finite_coefficient += reduce_relative(specialized) * finite_w**w_degree
    finite_equation += finite_coefficient * finite_x**x_degree

assert finite_equation.degree(finite_x) == 3
factorization = tuple(finite_equation.factor())
assert len(factorization) == 1 and factorization[0][1] == 1
singular.lib("normal.lib")
genus = int(finite_plane.ideal(finite_equation)._singular_().genus())

print(
    f"Q80THIRDCM24GENUS|prime=73|s=33|j=17|new_base={args.new_base}|"
    f"degrees_w,x={finite_equation.degree(finite_w)},{finite_equation.degree(finite_x)}|"
    f"irreducible=1|normalization_genus={genus}|"
    f"status={'PASS' if genus == 1 else 'FAIL_WRONG_GENUS'}",
    flush=True,
)
assert genus == 1
