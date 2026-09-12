#!/usr/bin/env python3
"""Verify the algebraic identities of the weighted stable-moduli core.

The geometric assertions use D1/F2 and the marked Hurwitz quotient.  This
checker independently guards the explicit ingredients used in the finite-
type chart formulation:

* normalized rerooting and its inverse/composition laws;
* transport of the Hessian divisor and inverse pencil; and
* the minimal weighted coordinate profiles in degrees four through six.
"""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.weighted import WeightedSeedModel, w, x, y, z  # noqa: E402
from verify_all_degree_coefficient_tangents import (  # noqa: E402
    explicit_seed,
)


def total_degree(polynomial: sp.Expr) -> int:
    return sp.Poly(polynomial, x, y, z).total_degree()


def reroot(polynomial: sp.Expr, root: sp.Expr) -> sp.Expr:
    multiplier = -1 / (root * sp.diff(polynomial, w).subs(w, root))
    return sp.factor(multiplier * polynomial.subs(w, root * w))


def verify_groupoid_identities() -> None:
    a, c, q, s, t = sp.symbols("a c q s t", nonzero=True)

    # A generic sextic chart with two named nonzero roots.  The remaining
    # factor keeps the calculation off the symmetric special loci.
    raw = w**2 * (w - 1) * (w - a) * (w - c) * (1 + q * w)
    H = sp.factor(-raw / sp.diff(raw, w).subs(w, 1))
    assert sp.factor(H.subs(w, 0)) == 0
    assert sp.factor(sp.diff(H, w).subs(w, 0)) == 0
    assert sp.factor(H.subs(w, 1)) == 0
    assert sp.factor(sp.diff(H, w).subs(w, 1) + 1) == 0

    H_a = reroot(H, a)
    assert sp.factor(H_a.subs(w, 1)) == 0
    assert sp.factor(sp.diff(H_a, w).subs(w, 1) + 1) == 0

    # The root c of H becomes c/a in R_a H.  Rerooting there is the same as
    # rerooting H directly at c; rerooting at 1/a is the inverse arrow.
    assert sp.factor(H_a.subs(w, c / a)) == 0
    assert sp.factor(reroot(H_a, c / a) - reroot(H, c)) == 0
    assert sp.factor(reroot(H_a, 1 / a) - H) == 0

    multiplier = -1 / (a * sp.diff(H, w).subs(w, a))
    assert sp.factor(
        sp.diff(H_a, w, 2)
        - multiplier * a**2 * sp.diff(H, w, 2).subs(w, a * w)
    ) == 0

    inverse_a = H_a - s * w + t
    transported = multiplier * (
        H.subs(w, a * w)
        - (s / (multiplier * a)) * (a * w)
        + t / multiplier
    )
    assert sp.factor(inverse_a - transported) == 0


def verify_coordinate_profiles() -> None:
    expected = {
        4: (12, 11, 4),
        5: (17, 16, 4),
        6: (22, 21, 4),
    }
    for degree, profile in expected.items():
        primitive = explicit_seed(degree)
        model = WeightedSeedModel(sp.diff(primitive, w))
        mapping = model.mapping()
        assert tuple(total_degree(component) for component in mapping) == profile
        assert sp.factor(
            sp.Matrix(mapping).jacobian((x, y, z)).det()
        ) == 1


def main() -> None:
    verify_groupoid_identities()
    verify_coordinate_profiles()
    print(
        "PASS weighted stable-moduli chart: rerooting groupoid, "
        "pencil transport, and degree profiles"
    )


if __name__ == "__main__":
    main()
