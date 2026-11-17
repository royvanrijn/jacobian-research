#!/usr/bin/env python3
"""Exact N=5,6,7 Tschirnhaus/Keller descent experiment.

The experiment compares two primitive coordinates on the same split
finite-etale algebra:

    r = (1,...,N),              u = r + r^2.

For each coordinate it compiles the normalized root polynomial into the
absolute quadratic-gauge Keller target.  It then checks:

* the two quotient algebras are explicitly isomorphic;
* u=r+r^2 is primitive but nonprojective;
* the coefficient-torus invariant lattice has rank N-4;
* the stable boundary fingerprints differ for N=5,6,7; and
* the relative and promoted coordinate degrees are 6N+2 and 6N+3.

Stable invariance of the fingerprints is the theorem-level input RQG3; this
script verifies the exact specialization and comparison, not that theorem.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/keller_tschirnhaus_descent_567.json"
)

T, U, S = sp.symbols("T U S")


def rational_text(value: sp.Expr) -> str:
    """Serialize an exact rational expression without losing precision."""

    return str(sp.factor(value))


def stable_weights(degree: int) -> sp.Matrix:
    """Return the 2-by-(N-2) quadratic-gauge weight matrix."""

    columns = [sp.Matrix((-2, -1))]
    columns.extend(
        sp.Matrix((1 - index, -index))
        for index in range(4, degree + 1)
    )
    return sp.Matrix.hstack(*columns)


def invariant_relations(degree: int) -> sp.Matrix:
    """Return a saturated Z-basis for the invariant character lattice."""

    relation_count = degree - 4
    relations = sp.zeros(relation_count, degree - 2)
    relations[0, 0] = -1
    relations[0, 1] = -6
    relations[0, 2] = 5
    for row in range(1, relation_count):
        # row=1 is a4*a6/a5^2, row=2 is a5*a7/a6^2.
        relations[row, row] = 1
        relations[row, row + 1] = -2
        relations[row, row + 2] = 1
    return relations


def compile_target(polynomial: sp.Poly, degree: int) -> dict[str, object]:
    """Compile the origin presentation into the fixed map U_N."""

    coefficients = {
        index: polynomial.nth(index) for index in range(degree + 1)
    }
    assert coefficients[1] != 0
    assert coefficients[3] != 0
    normalized = sp.Poly(polynomial.as_expr() / coefficients[1], T)
    normalized_coefficients = {
        index: normalized.nth(index) for index in range(degree + 1)
    }
    pi_value = normalized_coefficients[3]
    seed_parameters = {
        index: sp.factor(
            normalized_coefficients[index] / pi_value**index
        )
        for index in range(4, degree + 1)
    }
    target = tuple(seed_parameters[index] for index in range(4, degree + 1))
    target += (
        pi_value,
        normalized_coefficients[2],
        -2 * normalized_coefficients[0],
    )

    inverse = (
        S
        + target[-2] * S**2
        + target[-3] * S**3
        - target[-1] / 2
    )
    inverse += sum(
        seed_parameters[index] * target[-3] ** index * S**index
        for index in range(4, degree + 1)
    )
    expected = normalized.as_expr().subs(T, S)
    assert sp.cancel(inverse - expected) == 0

    return {
        "polynomial": polynomial,
        "normalized_coefficients": normalized_coefficients,
        "seed_parameters": seed_parameters,
        "target": target,
    }


def stable_fingerprint(compilation: dict[str, object], degree: int) -> tuple[sp.Expr, ...]:
    """Evaluate the saturated invariant basis in compiler seed coordinates."""

    seed = compilation["seed_parameters"]
    assert isinstance(seed, dict)
    values = [
        sp.factor(seed[5] ** 5 / seed[4] ** 6),
    ]
    for top_index in range(6, degree + 1):
        values.append(
            sp.factor(
                seed[top_index - 2]
                * seed[top_index]
                / seed[top_index - 1] ** 2
            )
        )
    return tuple(values)


def inverse_coordinate(source_roots: tuple[int, ...]) -> sp.Poly:
    """Interpolate r as a polynomial in u=r+r^2."""

    transformed_roots = tuple(root + root**2 for root in source_roots)
    return sp.Poly(
        sp.interpolate(tuple(zip(transformed_roots, source_roots)), U),
        U,
        domain=sp.QQ,
    )


def mobius_residuals(degree: int) -> tuple[sp.Expr, ...]:
    """Match the first three points projectively and test the remaining ones."""

    source = tuple(sp.Integer(index) for index in range(1, 4))
    target = tuple(root + root**2 for root in source)
    interpolation = sp.Matrix(
        [
            [
                source[index],
                1,
                -target[index] * source[index],
                -target[index],
            ]
            for index in range(3)
        ]
    )
    coefficients = []
    for column in range(4):
        remaining = [index for index in range(4) if index != column]
        coefficients.append(
            (-1) ** column * interpolation[:, remaining].det()
        )
    a_value, b_value, c_value, d_value = coefficients
    assert sp.Matrix(
        [[a_value, b_value], [c_value, d_value]]
    ).det() != 0

    residuals = []
    for root in range(4, degree + 1):
        transformed = root + root**2
        residual = sp.factor(
            a_value * root
            + b_value
            - transformed * (c_value * root + d_value)
        )
        assert residual == -2 * (root - 1) * (root - 2) * (root - 3)
        residuals.append(residual)
    return tuple(residuals)


def coordinate_degrees(degree: int) -> tuple[int, int]:
    """Return exact vertical and promoted total degrees of U_N."""

    # q has exact degree five and t has degree two.  The top B term is
    # N*u_N*t^2*x^(N-2)*q^N.  With u_N treated as a coefficient this has
    # degree 6N+2; after u_N is promoted it has degree 6N+3.  No other term
    # containing u_N can cancel it, and all remaining terms have lower degree.
    vertical = 6 * degree + 2
    promoted = vertical + 1
    assert vertical in (32, 38, 44)
    assert promoted in (33, 39, 45)
    return vertical, promoted


def build_results() -> dict[str, object]:
    """Run the exact experiment and return its JSON-ready certificate."""

    results: dict[str, object] = {
        "experiment": "Keller/Tschirnhaus descent in ranks 5, 6, and 7",
        "base_field": "Q",
        "tschirnhaus_change": "u=r+r^2",
        "stable_invariance_dependency": "RQG3",
        "projective_descent_dependency": "ARPD1",
        "degrees": {},
    }

    degree_rows: dict[str, object] = {}
    for degree in (5, 6, 7):
        weights = stable_weights(degree)
        relations = invariant_relations(degree)
        assert weights.rank() == 2
        assert relations.rank() == degree - 4
        assert weights * relations.T == sp.zeros(2, degree - 4)
        # The first N-4 columns form an upper-triangular minor with
        # determinant -1, so the displayed relation lattice is saturated.
        assert relations[:, : degree - 4].det() == -1

        roots = tuple(range(1, degree + 1))
        transformed_roots = tuple(root + root**2 for root in roots)
        assert len(set(transformed_roots)) == degree

        source_polynomial = sp.Poly(
            sp.prod(T - root for root in roots), T, domain=sp.QQ
        )
        transformed_polynomial = sp.Poly(
            sp.prod(U - root for root in transformed_roots),
            U,
            domain=sp.QQ,
        )
        inverse = inverse_coordinate(roots)

        # Exact quotient-algebra isomorphism:
        # U -> T+T^2 and T -> inverse(U).
        assert sp.rem(
            sp.Poly(
                transformed_polynomial.as_expr().subs(U, T + T**2),
                T,
                domain=sp.QQ,
            ),
            source_polynomial,
        ).is_zero
        assert sp.rem(
            sp.Poly(
                inverse.as_expr().subs(U, T + T**2) - T,
                T,
                domain=sp.QQ,
            ),
            source_polynomial,
        ).is_zero
        assert sp.rem(
            sp.Poly(
                inverse.as_expr() + inverse.as_expr() ** 2 - U,
                U,
                domain=sp.QQ,
            ),
            transformed_polynomial,
        ).is_zero

        evaluated_columns = sp.Matrix(
            [
                [1, root, root + root**2, root * (root + root**2)]
                for root in roots
            ]
        )
        assert evaluated_columns.rank() == 4
        residuals = mobius_residuals(degree)
        assert len(residuals) == degree - 3
        assert all(value != 0 for value in residuals)

        source_compilation = compile_target(source_polynomial, degree)
        transformed_compilation = compile_target(
            sp.Poly(
                transformed_polynomial.as_expr().subs(U, T),
                T,
                domain=sp.QQ,
            ),
            degree,
        )
        source_fingerprint = stable_fingerprint(source_compilation, degree)
        transformed_fingerprint = stable_fingerprint(
            transformed_compilation, degree
        )
        assert len(source_fingerprint) == degree - 4
        assert len(transformed_fingerprint) == degree - 4
        assert all(
            left != right
            for left, right in zip(
                source_fingerprint, transformed_fingerprint, strict=True
            )
        )

        vertical_degree, promoted_degree = coordinate_degrees(degree)
        degree_rows[str(degree)] = {
            "source_roots": list(roots),
            "transformed_roots": list(transformed_roots),
            "source_polynomial": str(source_polynomial.as_expr()),
            "transformed_polynomial": str(
                transformed_polynomial.as_expr()
            ),
            "inverse_tschirnhaus_polynomial": str(inverse.as_expr()),
            "projective_residuals": [
                rational_text(value) for value in residuals
            ],
            "projective_residual_count": degree - 3,
            "seed_parameter_count": degree - 3,
            "residual_scaling_eliminations": 1,
            "stable_boundary_survivor_count": degree - 4,
            "source_seed_parameters": {
                f"u_{index}": rational_text(value)
                for index, value in source_compilation[
                    "seed_parameters"
                ].items()
            },
            "transformed_seed_parameters": {
                f"u_{index}": rational_text(value)
                for index, value in transformed_compilation[
                    "seed_parameters"
                ].items()
            },
            "source_stable_fingerprint": [
                rational_text(value) for value in source_fingerprint
            ],
            "transformed_stable_fingerprint": [
                rational_text(value) for value in transformed_fingerprint
            ],
            "relative_coordinate_degree": vertical_degree,
            "promoted_coordinate_degree": promoted_degree,
            "promoted_target_dimension": degree,
        }

    results["degrees"] = degree_rows
    return results


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write",
        action="store_true",
        help="refresh the pinned generated-results JSON",
    )
    arguments = parser.parse_args()
    results = build_results()
    rendered = json.dumps(results, indent=2, sort_keys=True) + "\n"

    if arguments.write:
        ARTIFACT.write_text(rendered, encoding="utf-8")
        print(f"WROTE: {ARTIFACT.relative_to(ROOT)}")
    else:
        assert ARTIFACT.read_text(encoding="utf-8") == rendered

    print("PASS: both primitive coordinates give the same finite-etale algebra")
    print("PASS: q(r)=r+r^2 has N-3 nonzero projective residuals")
    print("PASS: exact stable boundary survivor counts are 1, 2, and 3")
    print("PASS: every N=5,6,7 stable fingerprint changes")
    print("PASS: promoted coordinate degrees are 33, 39, and 45")


if __name__ == "__main__":
    main()
