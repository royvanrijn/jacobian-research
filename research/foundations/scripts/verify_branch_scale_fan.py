#!/usr/bin/env python3
"""Exact degree-six test of the branch-scale fan conjecture.

The radial calculation is toric: for three pair clusters the moving critical
values have weights 2*v(x), 2*v(y), 2*v(z), and the common refinement of
their pairwise ratio graphs is the A_2 braid fan.  The final calculation
checks a triple leading-value resonance and proves that the full stable-target
graph needs a further (non-radial) refinement.
"""

from __future__ import annotations

import itertools
import json
import sys
from fractions import Fraction
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "artifacts" / "generated-results" / "branch_scale_fan_degree6.json"

W, X, L, t = sp.symbols("W X L t")
x, y, z = sp.symbols("x y z")

CENTERS = (sp.Integer(0), sp.Integer(1), sp.Integer(3))
SCALES = (x, y, z)
LABELS = ("x", "y", "z")


def t_order(expression: sp.Expr) -> int:
    """Return the t-adic order of a nonzero polynomial/rational series."""

    numerator, denominator = sp.cancel(expression).as_numer_denom()

    def polynomial_order(polynomial: sp.Expr) -> int:
        poly = sp.Poly(polynomial, t)
        return min(monomial[0] for monomial, coefficient in poly.terms() if coefficient)

    return polynomial_order(numerator) - polynomial_order(denominator)


def lower_newton_root_valuations(polynomial: sp.Poly) -> list[Fraction]:
    """Extract root valuations from the lower Newton polygon."""

    points = []
    for exponent in range(polynomial.degree() + 1):
        coefficient = polynomial.nth(exponent)
        if coefficient != 0:
            points.append((exponent, t_order(coefficient)))

    hull: list[tuple[int, int]] = []
    for point in points:
        while len(hull) >= 2:
            x0, y0 = hull[-2]
            x1, y1 = hull[-1]
            x2, y2 = point
            left_slope = Fraction(y1 - y0, x1 - x0)
            right_slope = Fraction(y2 - y1, x2 - x1)
            if left_slope >= right_slope:
                hull.pop()
            else:
                break
        hull.append(point)

    valuations: list[Fraction] = []
    for (x0, y0), (x1, y1) in zip(hull, hull[1:]):
        valuation = -Fraction(y1 - y0, x1 - x0)
        valuations.extend([valuation] * (x1 - x0))
    return sorted(valuations)


def degree_six_family(
    x_value: sp.Expr = x,
    y_value: sp.Expr = y,
    z_value: sp.Expr = z,
) -> sp.Expr:
    """Three pair clusters centered at 0, 1, and 3."""

    return sp.expand(
        W
        * (W - x_value)
        * (W - 1)
        * (W - 1 - y_value)
        * (W - 3)
        * (W - 3 - z_value)
    )


family = degree_six_family()

# Each local pair polynomial is X(X-1).  The other four factors supply a
# nonzero center-dependent unit, so the moving critical values have scales
# x^2, y^2, z^2 and the displayed leading constants.
leading_constants: list[sp.Rational] = []
local_models: list[sp.Expr] = []
for center, scale in zip(CENTERS, SCALES):
    rescaled = sp.cancel(family.subs(W, center + scale * X) / scale**2)
    local_model = sp.factor(rescaled.subs({x: 0, y: 0, z: 0}))
    local_models.append(local_model)
    assert sp.solve(sp.diff(local_model, X), X) == [sp.Rational(1, 2)]
    leading_constants.append(sp.factor(local_model.subs(X, sp.Rational(1, 2))))

assert local_models == [
    9 * X * (X - 1),
    4 * X * (X - 1),
    36 * X * (X - 1),
]
assert leading_constants == [
    -sp.Rational(9, 4),
    -sp.Integer(1),
    -sp.Integer(9),
]

# The special derivative has exactly three cluster critical points with value
# zero and two global critical points with nonzero values.  Hence these are
# all of the moving critical values.
special_family = sp.factor(family.subs({x: 0, y: 0, z: 0}))
special_derivative = sp.factor(sp.diff(special_family, W))
assert special_derivative == (
    2 * W * (W - 3) * (W - 1) * (3 * W**2 - 8 * W + 3)
)
global_critical_points = sp.solve(3 * W**2 - 8 * W + 3, W)
assert len(global_critical_points) == 2
assert all(special_family.subs(W, point) != 0 for point in global_critical_points)

# An independent eliminant/Newton-polygon calculation recovers the two
# stationary valuations and the three moving valuations for every strict
# ordering of the cluster rates.
newton_traits = []
for exponents in itertools.permutations((1, 2, 3)):
    trait_family = degree_six_family(*(t**exponent for exponent in exponents))
    critical_value_eliminant = sp.Poly(
        sp.resultant(sp.diff(trait_family, W), L - trait_family, W),
        L,
    )
    root_valuations = lower_newton_root_valuations(critical_value_eliminant)
    assert root_valuations == [
        Fraction(0),
        Fraction(0),
        Fraction(2),
        Fraction(4),
        Fraction(6),
    ]
    newton_traits.append(
        {
            "scale_valuations": dict(zip(LABELS, exponents)),
            "critical_value_valuations": [str(value) for value in root_valuations],
        }
    )

# The common normalized graph of the three initial pairwise ratios has the
# braid-fan subdivision of the positive octant.  For an ordering
# v_{i0} <= v_{i1} <= v_{i2}, its rays are the full diagonal, the indicator
# of {i1,i2}, and the indicator of {i2}.  Every cone is unimodular.
fan_cones = []
for order in itertools.permutations(range(3)):
    full = (1, 1, 1)
    upper_pair = tuple(int(index in order[1:]) for index in range(3))
    upper_singleton = tuple(int(index == order[2]) for index in range(3))
    ray_matrix = sp.Matrix.hstack(
        sp.Matrix(full),
        sp.Matrix(upper_pair),
        sp.Matrix(upper_singleton),
    )
    assert abs(ray_matrix.det()) == 1
    fan_cones.append(
        {
            "weight_order": [LABELS[index] for index in order],
            "inequalities": [
                f"2v_{LABELS[order[0]]} <= 2v_{LABELS[order[1]]}",
                f"2v_{LABELS[order[1]]} <= 2v_{LABELS[order[2]]}",
            ],
            "rays": [full, upper_pair, upper_singleton],
        }
    )

assert {
    ray
    for cone in fan_cones
    for ray in map(tuple, cone["rays"])
} == {
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, 1, 0),
    (1, 0, 1),
    (0, 1, 1),
    (1, 1, 1),
}

# On a strict chamber the target points 0, lambda_x, lambda_y, lambda_z have
# pairwise-difference valuations min(2v_i,2v_j).  Thus the six fan chambers
# are exactly the six nonresonant caterpillar types of the stabilized target.
for exponents in itertools.permutations((1, 2, 3)):
    weights = tuple(2 * exponent for exponent in exponents)
    for left, right in itertools.combinations(range(3), 2):
        assert weights[left] != weights[right]
        assert min(weights[left], weights[right]) in (2, 4)


def critical_value_mod_t4(
    trait_family: sp.Expr,
    center: sp.Rational,
    first_displacement: sp.Rational,
) -> sp.Expr:
    """Critical value through t^3 for a pair cluster on a linear trait."""

    beta = sp.symbols("beta")
    critical_point = center + first_displacement * t + beta * t**2
    derivative_series = sp.series(
        sp.diff(trait_family, W).subs(W, critical_point),
        t,
        0,
        4,
    ).removeO().expand()
    beta_equation = next(
        derivative_series.coeff(t, exponent)
        for exponent in range(4)
        if beta in derivative_series.coeff(t, exponent).free_symbols
    )
    beta_value = sp.solve(beta_equation, beta)[0]
    critical_point = critical_point.subs(beta, beta_value)
    return sp.factor(
        sp.series(trait_family.subs(W, critical_point), t, 0, 4).removeO()
    )


# Pairwise Maxwell initial forms live on the equality walls.  Their possible
# cancellation is information not contained in the radial valuation vector.
maxwell_initial_forms = {
    "lambda_x-lambda_y": sp.factor(
        leading_constants[0] * x**2 - leading_constants[1] * y**2
    ),
    "lambda_x-lambda_z": sp.factor(
        leading_constants[0] * x**2 - leading_constants[2] * z**2
    ),
    "lambda_y-lambda_z": sp.factor(
        leading_constants[1] * y**2 - leading_constants[2] * z**2
    ),
}
expected_maxwell_initial_forms = {
    "lambda_x-lambda_y": -(3 * x - 2 * y) * (3 * x + 2 * y) / 4,
    "lambda_x-lambda_z": -9 * (x - 2 * z) * (x + 2 * z) / 4,
    "lambda_y-lambda_z": -(y - 3 * z) * (y + 3 * z),
}
assert all(
    sp.expand(maxwell_initial_forms[key] - value) == 0
    for key, value in expected_maxwell_initial_forms.items()
)

# At the triple leading resonance all three values equal -9/4*t^2.  The
# second-scale cross-ratio varies with the second jets A,B,C although the
# complete radial fan datum (including the leading scale ratios) is fixed.
A, B, C = sp.symbols("A B C")
triple_resonance_family = degree_six_family(
    t + A * t**2,
    sp.Rational(3, 2) * t + B * t**2,
    sp.Rational(1, 2) * t + C * t**2,
)
triple_values = [
    critical_value_mod_t4(
        triple_resonance_family,
        center,
        first_displacement,
    )
    for center, first_displacement in zip(
        CENTERS,
        (sp.Rational(1, 2), sp.Rational(3, 4), sp.Rational(1, 4)),
    )
]
assert all(
    sp.expand(value).coeff(t, 2) == -sp.Rational(9, 4)
    for value in triple_values
)

lambda_x_minus_y = sp.factor(triple_values[0] - triple_values[1])
lambda_z_minus_y = sp.factor(triple_values[2] - triple_values[1])
second_scale_cross_ratio = sp.factor(lambda_x_minus_y / lambda_z_minus_y)
assert lambda_x_minus_y == -3 * t**3 * (6 * A - 4 * B + 1) / 4
assert lambda_z_minus_y == 3 * t**3 * (2 * B - 6 * C + 1) / 2
assert sp.cancel(
    second_scale_cross_ratio
    + (6 * A - 4 * B + 1) / (2 * (2 * B - 6 * C + 1))
) == 0
assert second_scale_cross_ratio.subs({A: 0, B: 0, C: 0}) == -sp.Rational(1, 2)
assert second_scale_cross_ratio.subs({A: 1, B: 0, C: 0}) == -sp.Rational(7, 2)

artifact = {
    "experiment": "degree-six three-pair-cluster branch-scale fan",
    "family": "W(W-x)(W-1)(W-1-y)(W-3)(W-3-z)",
    "cluster_multiplicities": [2, 2, 2],
    "moving_critical_value_initial_terms": {
        "lambda_x": "-9/4*x^2",
        "lambda_y": "-y^2",
        "lambda_z": "-9*z^2",
    },
    "valuation_vectors": {
        "lambda_x": [2, 0, 0],
        "lambda_y": [0, 2, 0],
        "lambda_z": [0, 0, 2],
        "lambda_x/lambda_y": [2, -2, 0],
        "lambda_x/lambda_z": [2, 0, -2],
        "lambda_y/lambda_z": [0, 2, -2],
    },
    "radial_fan": {
        "hyperplanes": ["2v_x=2v_y", "2v_x=2v_z", "2v_y=2v_z"],
        "rays_in_positive_octant": [
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 1, 0],
            [1, 0, 1],
            [0, 1, 1],
            [1, 1, 1],
        ],
        "maximal_cones": fan_cones,
        "all_maximal_cones_unimodular": True,
    },
    "newton_traits": newton_traits,
    "maxwell_initial_forms": {
        key: str(value) for key, value in maxwell_initial_forms.items()
    },
    "triple_resonance": {
        "scale_jets": [
            "x=t+A*t^2",
            "y=3/2*t+B*t^2",
            "z=1/2*t+C*t^2",
        ],
        "common_leading_critical_value": "-9/4*t^2",
        "lambda_x-lambda_y": str(lambda_x_minus_y),
        "lambda_z-lambda_y": str(lambda_z_minus_y),
        "second_scale_cross_ratio": str(second_scale_cross_ratio),
        "specializations": {
            "A=B=C=0": "-1/2",
            "A=1,B=C=0": "-7/2",
        },
        "conclusion": (
            "the radial weighted braid fan is the first toric layer, but the "
            "full stable-target graph must refine its triple-resonance locus"
        ),
    },
}

expected_artifact = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
if "--emit-json" in sys.argv:
    print(expected_artifact, end="")
    raise SystemExit(0)
assert ARTIFACT.read_text() == expected_artifact, (
    f"{ARTIFACT.relative_to(ROOT)} is stale; regenerate it from this script"
)

print("PASS degree-six moving critical scales: x^2, y^2, z^2")
print("PASS radial normalized graph: six unimodular A2 braid-fan cones")
print("PASS eliminant Newton polygons: stationary 0,0 and moving 2a,2b,2c")
print("PASS triple resonance: fixed radial data, varying second-scale cross-ratio")
print("BRANCH_SCALE_FAN_DEGREE6_PASS")
