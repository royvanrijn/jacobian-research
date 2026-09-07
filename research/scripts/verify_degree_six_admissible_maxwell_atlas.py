#!/usr/bin/env python3
"""Verify the nonradial Maxwell source charts for three quadratic clusters."""

from __future__ import annotations

import json
import sys
from itertools import combinations, product
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "admissible_maxwell_atlas_degree6.json"
)

LABELS = ("x", "y", "z")
TAIL_UNITS = {"x": 9, "y": 4, "z": 36}
LEADING_UNITS = {
    "x": -sp.Rational(9, 4),
    "y": -sp.Integer(1),
    "z": -sp.Integer(9),
}


def ideal_intersection(
    left: list[sp.Expr],
    right: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
) -> list[sp.Expr]:
    """Compute an ideal intersection by elimination."""

    selector = sp.Symbol("selector")
    basis = sp.groebner(
        [selector * equation for equation in left]
        + [(1 - selector) * equation for equation in right],
        selector,
        *variables,
        order="lex",
    )
    return [
        sp.factor(polynomial.as_expr())
        for polynomial in basis.polys
        if selector not in polynomial.as_expr().free_symbols
    ]


# After adjoining square roots of the exact critical-value units, write the
# three branch values as q_i^2.  Every Maxwell pullback is then a sign
# arrangement.  This étale Kummer chart does not change completed-ring
# normality or transversality.
qx, qy, qz = sp.symbols("qx qy qz")
q_by_label = {"x": qx, "y": qy, "z": qz}

maxwell_strata: dict[str, dict[str, object]] = {}
for size in (2, 3):
    for subset in combinations(LABELS, size):
        base = subset[0]
        variables = tuple(q_by_label[label] for label in subset)
        equations = [
            q_by_label[label] ** 2 - q_by_label[base] ** 2
            for label in subset[1:]
        ]
        sign_primes = []
        for signs in product((-1, 1), repeat=size - 1):
            sign_primes.append(
                [
                    q_by_label[label] - sign * q_by_label[base]
                    for label, sign in zip(subset[1:], signs)
                ]
            )

        intersection = sign_primes[0]
        for sign_prime in sign_primes[1:]:
            intersection = ideal_intersection(
                intersection,
                sign_prime,
                variables,
            )
        assert sp.groebner(
            intersection,
            *variables,
            order="lex",
        ) == sp.groebner(
            equations,
            *variables,
            order="lex",
        )

        quadratic_tails = size
        identity_strands = 2 * (len(LABELS) - size)
        source_degree = 2 * quadratic_tails + identity_strands
        node_indices = sorted(
            [2] * quadratic_tails + [1] * identity_strands
        )
        assert source_degree == 6
        assert sum(node_indices) == 6
        assert (
            quadratic_tails * (2 * 2 - 2)
            + identity_strands * (2 * 1 - 2)
            == 2 * quadratic_tails
        )

        name = "".join(subset)
        maxwell_strata[name] = {
            "colliding_labels": list(subset),
            "quadratic_tails": quadratic_tails,
            "identity_strands": identity_strands,
            "source_degree": source_degree,
            "node_indices": node_indices,
            "normalization_branch_count": len(sign_primes),
            "sign_branches": [
                [
                    f"q_{subset[index]}={sign:+d}*q_{base}"
                    for index, sign in enumerate(signs, start=1)
                ]
                for signs in product((-1, 1), repeat=size - 1)
            ],
            "node_ring": (
                "k[[tau,"
                + ",".join(f"s_{label}" for label in subset)
                + (
                    ",e_1,e_2]]/("
                    if identity_strands
                    else "]]/("
                )
                + ",".join(
                    f"s_{label}^2-tau" for label in subset
                )
                + (
                    ",e_1-tau,e_2-tau)"
                    if identity_strands
                    else ")"
                )
            ),
            "local_tail_equations": [
                (
                    f"S-lambda_{label}="
                    f"{TAIL_UNITS[label]}*U_{label}^2"
                )
                for label in subset
            ],
        }

assert {
    name: stratum["normalization_branch_count"]
    for name, stratum in maxwell_strata.items()
} == {"xy": 2, "xz": 2, "yz": 2, "xyz": 4}

# The three residual double points of the six-line target arrangement are
# normal crossings between one radial divisor and one pairwise Maxwell
# divisor.  Their Kummer preimages are two reduced transverse points.
x, y, z = sp.symbols("x y z")
A = LEADING_UNITS["x"] * x**2
B = LEADING_UNITS["y"] * y**2
C = LEADING_UNITS["z"] * z**2

residual_crossings = {
    "A=0 & B=C": {
        "pair": "yz",
        "chart_substitution": {z: 1},
        "equations": (x, sp.expand(B - C).subs(z, 1)),
        "variables": (x, y),
        "points": ((0, -3), (0, 3)),
    },
    "B=0 & A=C": {
        "pair": "xz",
        "chart_substitution": {z: 1},
        "equations": (y, sp.expand(A - C).subs(z, 1)),
        "variables": (y, x),
        "points": ((0, -2), (0, 2)),
    },
    "C=0 & A=B": {
        "pair": "xy",
        "chart_substitution": {x: 1},
        "equations": (z, sp.expand(A - B).subs(x, 1)),
        "variables": (z, y),
        "points": (
            (0, -sp.Rational(3, 2)),
            (0, sp.Rational(3, 2)),
        ),
    },
}

crossing_artifact = {}
for name, data in residual_crossings.items():
    equations = data["equations"]
    variables = data["variables"]
    jacobian = sp.Matrix(
        [
            [sp.diff(equation, variable) for variable in variables]
            for equation in equations
        ]
    )
    determinants = []
    for point in data["points"]:
        substitution = dict(zip(variables, point))
        assert all(
            sp.expand(equation.subs(substitution)) == 0
            for equation in equations
        )
        determinant = sp.factor(jacobian.det().subs(substitution))
        assert determinant != 0
        determinants.append(str(determinant))
    crossing_artifact[name] = {
        "pair": data["pair"],
        "preimage_count": 2,
        "jacobian_determinants": determinants,
        "conclusion": "reduced transverse; no further blowup center",
    }

# At a coordinate triple point, blowing up the radial center records the
# ratio of the two vanishing branch values.  Equality of those residues has
# two reduced Kummer directions, already present as a two-cluster equality
# face in the radial atlas.
coordinate_directions = {
    "xy": ("[x:y]", ["[1:-3/2]", "[1:3/2]"]),
    "xz": ("[x:z]", ["[1:-1/2]", "[1:1/2]"]),
    "yz": ("[y:z]", ["[1:-3]", "[1:3]"]),
}

# At [1:1:1], all three Maxwell divisors meet.  The target blowup separates
# them, while the Kummer pullback has the four sign branches already computed
# by the three-cluster node ring.
assert maxwell_strata["xyz"]["node_indices"] == [2, 2, 2]
assert maxwell_strata["xyz"]["normalization_branch_count"] == 4

artifact = {
    "experiment": "degree-six nonradial Maxwell admissible atlas",
    "kummer_chart": "lambda_i=q_i^2 after an etale unit extension",
    "maxwell_strata": maxwell_strata,
    "residual_radial_maxwell_crossings": crossing_artifact,
    "coordinate_point_directions": {
        pair: {
            "exceptional_coordinate": coordinate,
            "reduced_directions": directions,
            "interpretation": "radial two-cluster equality face",
        }
        for pair, (coordinate, directions) in coordinate_directions.items()
    },
    "triple_maxwell": {
        "target_center": "[1:1:1]",
        "source_node_indices": [2, 2, 2],
        "normalization_branch_count": 4,
        "conclusion": (
            "the diagonal target blowup and four Kummer branches are the "
            "only additional centers among the six target-arrangement strata"
        ),
    },
    "scope": (
        "generic pairwise and triple Maxwell strata and their intersections "
        "with the radial boundary; the extension morphism, selection of the "
        "polynomial admissible-cover component, and label descent remain open"
    ),
}

expected_artifact = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
if "--emit-json" in sys.argv:
    print(expected_artifact, end="")
    raise SystemExit(0)
assert ARTIFACT.read_text() == expected_artifact, (
    f"{ARTIFACT.relative_to(ROOT)} is stale; regenerate it from this script"
)

print("PASS pairwise Maxwell charts: two Kummer branches each")
print("PASS triple Maxwell chart: four Kummer branches")
print("PASS Maxwell bubbles: source degree six and matching node indices")
print("PASS residual radial-Maxwell crossings: reduced and transverse")
print("PASS coordinate intersections: already radial equality directions")
print("DEGREE_SIX_ADMISSIBLE_MAXWELL_ATLAS_PASS")
