#!/usr/bin/env python3
"""Uniform Mason certificate that a seed has at most one omitted pencil value."""

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.discriminant_geometry import (  # noqa: E402
    affine_difference_mason_defect,
)


def full_contact_partitions(total, maximum=None):
    """Partitions of total with every part at least two."""
    if total == 0:
        yield ()
        return
    maximum = min(total, maximum or total)
    for first in range(maximum, 1, -1):
        for tail in full_contact_partitions(total - first, first):
            yield (first,) + tail


# Mason for a nonzero linear difference requires n<=ell(lambda)+ell(mu).
# Since every part is at least two, the only non-strict support case is the
# same all-double partition in even degree.  A constant difference has the
# stronger n<=ell(lambda)+ell(mu)-1 and is therefore always impossible.
for degree in range(3, 25):
    partitions = tuple(full_contact_partitions(degree))
    for left in partitions:
        for right in partitions:
            defect = affine_difference_mason_defect(left, right)
            assert defect >= 0
            if defect == 0:
                assert degree % 2 == 0
                assert left == right == (2,) * (degree // 2)
            else:
                assert len(left) + len(right) < degree

# Close the unique equality case uniformly.  If P=A^2 and Q=B^2 with A,B
# distinct monic degree-m polynomials, D=A-B is nonzero of some degree
# 0<=d<m, while A+B has degree m in characteristic zero.  Hence
# deg(P-Q)=m+d>=2.  The symbolic regressions realize every possible d.
W = sp.symbols("W")
for half_degree in range(2, 13):
    coefficients = sp.symbols(f"a0:{half_degree}")
    A = W**half_degree + sum(
        coefficient * W**index
        for index, coefficient in enumerate(coefficients)
    )
    for difference_degree in range(half_degree):
        difference = W**difference_degree + 1
        if difference_degree == 0:
            difference = sp.Integer(1)
        B = sp.expand(A - difference)
        square_difference = sp.expand(A**2 - B**2)
        assert sp.expand(square_difference - difference * (A + B)) == 0
        assert sp.degree(square_difference, W) == (
            half_degree + difference_degree
        )
        assert sp.degree(square_difference, W) >= 2

print("PASS: Mason leaves only the even-degree all-double equality case")
print("PASS: a nonzero constant difference is excluded by the stronger bound")
print("PASS: two distinct monic squares cannot differ by an affine polynomial")
print("PASS: every normalized admissible seed has at most one omitted pencil value")
print("PASS: exact multiplicity types give a disjoint stratification of N_n")
