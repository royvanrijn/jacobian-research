#!/usr/bin/env python3
"""Exact line obstruction and bounded quadratic-curve search at sqrt(-31).

Work on the fixed-discriminant common-resolvent double cover

    Z^2 = 31 M(A, (A^2+31 R^2)/4, Pi^3)

through (A,R,Pi,Z)=(-8,2,5,1984).  The line calculation is exhaustive over
the algebraic closure: after recursively determining the degree-four square
root, the residual equations have no nonzero projective direction.  Exact
good-prime projective certificates also exclude degree-at-most-two curves
on all three coordinate-fixed slices.  The unrestricted quadratic
calculation is only a bounded search in integral coefficients.  It is
discovery evidence, not a proof that no rational curve exists and not an
infinitude theorem for Hasse failures.
"""

from __future__ import annotations

import argparse
import itertools
import json
import shutil
import subprocess
from fractions import Fraction
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "fixed_quintic_hasse_curve_search.json"
)

A, V, W, R, P, Z, t = sp.symbols("A V W R P Z t")
a, r, p = sp.symbols("a r p")

M = (
    3 * A**8
    - 24 * A**6 * V
    - 50 * A**6
    + 70 * A**4 * V**2
    + 270 * A**4 * V
    + 56 * A**4 * W
    + 275 * A**4
    - 76 * A**2 * V**3
    - 510 * A**2 * V**2
    - 288 * A**2 * V * W
    - 750 * A**2 * V
    - 360 * A**2 * W
    - 500 * A**2
    + 27 * V**4
    + 270 * V**3
    + 216 * V**2 * W
    + 675 * V**2
    + 1080 * V * W
    + 432 * W**2
)

FIXED_V = (A**2 + 31 * R**2) / 4
DOUBLE_COVER_RHS = sp.factor(31 * M.subs({V: FIXED_V, W: P**3}))
INTEGRAL_RHS = sp.Poly(
    sp.expand(256 * DOUBLE_COVER_RHS),
    A,
    R,
    P,
)
INTEGRAL_TERMS = [
    (monomial, int(coefficient))
    for monomial, coefficient in INTEGRAL_RHS.terms()
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quadratic-bound", type=int, default=2)
    parser.add_argument(
        "--write",
        action="store_true",
        help="write the default-bound generated artifact",
    )
    parser.add_argument("--emit-json", action="store_true")
    return parser.parse_args()


def primitive_residual(poly: sp.Expr) -> sp.Poly:
    cleared = sp.Poly(sp.together(poly), a, r, p).clear_denoms()[1]
    return cleared.primitive()[1]


def exact_line_obstruction() -> dict[str, object]:
    line_rhs = sp.Poly(
        sp.expand(
            DOUBLE_COVER_RHS.subs(
                {
                    A: -8 + a * t,
                    R: 2 + r * t,
                    P: 5 + p * t,
                }
            )
        ),
        t,
    )
    assert line_rhs.degree() == 8
    coefficients = [line_rhs.coeff_monomial(t**degree) for degree in range(9)]
    assert coefficients[0] == 1984**2

    square_root_coefficients = [sp.Integer(1984)]
    for degree in range(1, 5):
        earlier_products = sum(
            square_root_coefficients[index]
            * square_root_coefficients[degree - index]
            for index in range(1, degree)
        )
        square_root_coefficients.append(
            sp.factor(
                (coefficients[degree] - earlier_products)
                / (2 * square_root_coefficients[0])
            )
        )

    residuals: list[sp.Poly] = []
    for degree in range(5, 9):
        square_coefficient = sum(
            square_root_coefficients[index]
            * square_root_coefficients[degree - index]
            for index in range(max(0, degree - 4), min(4, degree) + 1)
        )
        residual = primitive_residual(coefficients[degree] - square_coefficient)
        assert residual.total_degree() == degree
        residuals.append(residual)

    chart_results: dict[str, object] = {}
    for chart_variable, remaining in (
        (a, (r, p)),
        (r, (a, p)),
        (p, (a, r)),
    ):
        chart_polynomials = [
            residual.as_expr().subs(chart_variable, 1)
            for residual in residuals
        ]
        basis = sp.groebner(
            chart_polynomials,
            *remaining,
            order="grevlex",
        )
        assert len(basis.polys) == 1
        assert basis.polys[0].as_expr() == 1
        chart_results[str(chart_variable)] = {
            "normalization": f"{chart_variable}=1",
            "groebner_basis": ["1"],
        }

    return {
        "base_ansatz": [
            "A(t)=-8+a*t",
            "R(t)=2+r*t",
            "Pi(t)=5+p*t",
        ],
        "square_root_degree": 4,
        "residual_degrees": [residual.total_degree() for residual in residuals],
        "residual_monomial_counts": [
            len(residual.terms()) for residual in residuals
        ],
        "projective_direction_charts": chart_results,
        "conclusion": (
            "no nonconstant base line through (-8,2,5) lifts to the "
            "fixed-sqrt(-31) common-resolvent double cover"
        ),
        "scope": (
            "exact only for base coordinates affine-linear in one parameter"
        ),
    }


def quadratic_coordinate_slice_obstructions() -> dict[str, object]:
    singular = shutil.which("Singular")
    if singular is None:
        raise SystemExit(
            "Singular is required for the quadratic coordinate-slice obstructions"
        )

    a1, a2, r1, r2, p1, p2 = sp.symbols("a1 a2 r1 r2 p1 p2")
    coordinate_polynomials = {
        A: -8 + a1 * t + a2 * t**2,
        R: 2 + r1 * t + r2 * t**2,
        P: 5 + p1 * t + p2 * t**2,
    }
    slice_specs = {
        "fixed_A": {
            "variables": (r1, r2, p1, p2),
            "substitution": {
                A: -8,
                R: coordinate_polynomials[R],
                P: coordinate_polynomials[P],
            },
            "base_ansatz": [
                "A(t)=-8",
                "R(t)=2+r1*t+r2*t^2",
                "Pi(t)=5+p1*t+p2*t^2",
            ],
            "conclusion": (
                "no nonconstant degree-at-most-two base curve through "
                "(-8,2,5) with A(t)=-8 lifts to the double cover"
            ),
        },
        "fixed_R": {
            "variables": (a1, a2, p1, p2),
            "substitution": {
                A: coordinate_polynomials[A],
                R: 2,
                P: coordinate_polynomials[P],
            },
            "base_ansatz": [
                "A(t)=-8+a1*t+a2*t^2",
                "R(t)=2",
                "Pi(t)=5+p1*t+p2*t^2",
            ],
            "conclusion": (
                "no nonconstant degree-at-most-two base curve through "
                "(-8,2,5) with R(t)=2 lifts to the double cover"
            ),
        },
        "fixed_Pi": {
            "variables": (a1, a2, r1, r2),
            "substitution": {
                A: coordinate_polynomials[A],
                R: coordinate_polynomials[R],
                P: 5,
            },
            "base_ansatz": [
                "A(t)=-8+a1*t+a2*t^2",
                "R(t)=2+r1*t+r2*t^2",
                "Pi(t)=5",
            ],
            "conclusion": (
                "no nonconstant degree-at-most-two base curve through "
                "(-8,2,5) with Pi(t)=5 lifts to the double cover"
            ),
        },
    }

    prime = 32003
    results: dict[str, object] = {}
    for slice_name, spec in slice_specs.items():
        variables = spec["variables"]
        weights = {
            variable: 2 if str(variable).endswith("2") else 1
            for variable in variables
        }
        quadratic_rhs = sp.Poly(
            sp.expand(DOUBLE_COVER_RHS.subs(spec["substitution"])),
            t,
        )
        assert quadratic_rhs.degree() == 16
        coefficients = [
            quadratic_rhs.coeff_monomial(t**degree)
            for degree in range(17)
        ]
        assert coefficients[0] == 1984**2

        square_root_coefficients = [sp.Integer(1984)]
        for degree in range(1, 9):
            earlier_products = sum(
                square_root_coefficients[index]
                * square_root_coefficients[degree - index]
                for index in range(1, degree)
            )
            square_root_coefficients.append(
                sp.cancel(
                    (coefficients[degree] - earlier_products)
                    / (2 * square_root_coefficients[0])
                )
            )

        residuals: list[sp.Poly] = []
        for degree in range(9, 17):
            square_coefficient = sum(
                square_root_coefficients[index]
                * square_root_coefficients[degree - index]
                for index in range(max(0, degree - 8), min(8, degree) + 1)
            )
            cleared = sp.Poly(
                sp.together(coefficients[degree] - square_coefficient),
                *variables,
            ).clear_denoms()[1]
            residual = cleared.primitive()[1]
            for monomial, _ in residual.terms():
                weighted_degree = sum(
                    exponent * weights[variable]
                    for variable, exponent in zip(variables, monomial)
                )
                assert weighted_degree == degree
            residuals.append(residual)

        chart_results: dict[str, object] = {}
        for chart_variable in variables:
            remaining = tuple(
                variable
                for variable in variables
                if variable != chart_variable
            )
            expressions = [
                str(
                    sp.expand(
                        residual.as_expr().subs(chart_variable, 1)
                    )
                ).replace("**", "^")
                for residual in residuals
            ]
            program = (
                f"ring q={prime},({','.join(map(str, remaining))}),dp;\n"
                "option(redSB);\n"
                f"ideal I={','.join(expressions)};\n"
                "ideal G=std(I);\n"
                "size(G);\n"
                "G[1];\n"
            )
            completed = subprocess.run(
                [singular, "-q"],
                input=program,
                text=True,
                capture_output=True,
                timeout=60,
                check=True,
            )
            assert "?" not in completed.stdout
            output = [
                line.strip()
                for line in completed.stdout.splitlines()
                if line.strip()
            ]
            assert output[-2:] == ["1", "1"]
            chart_results[str(chart_variable)] = {
                "normalization": f"{chart_variable}=1",
                "groebner_basis": ["1"],
            }

        results[slice_name] = {
            "base_ansatz": spec["base_ansatz"],
            "coefficient_weights": {
                str(variable): weight
                for variable, weight in weights.items()
            },
            "residual_degrees": list(range(9, 17)),
            "residual_monomial_counts": [
                len(residual.terms()) for residual in residuals
            ],
            "weighted_projective_charts": chart_results,
            "special_fiber": "empty",
            "conclusion": spec["conclusion"],
            "scope": f"exact on the {slice_name} coordinate slice",
        }

    return {
        "good_prime": prime,
        "slices": results,
        "characteristic_zero_argument": (
            "each weighted-homogeneous residual scheme is projective over "
            "Z; its empty good-prime fiber forces the generic fiber to be empty"
        ),
        "conclusion": (
            "a degree-at-most-two curve through the known point must "
            "genuinely vary A, R, and Pi"
        ),
        "scope": "exact on all three coordinate-fixed quadratic slices",
    }


def add_polynomials(left: list[int], right: list[int]) -> list[int]:
    result = [0] * max(len(left), len(right))
    for index, value in enumerate(left):
        result[index] += value
    for index, value in enumerate(right):
        result[index] += value
    return result


def multiply_polynomials(left: list[int], right: list[int]) -> list[int]:
    result = [0] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        if not left_value:
            continue
        for right_index, right_value in enumerate(right):
            if right_value:
                result[left_index + right_index] += left_value * right_value
    return result


def polynomial_powers(poly: list[int], exponent: int) -> list[list[int]]:
    result = [[1]]
    for _ in range(exponent):
        result.append(multiply_polynomials(result[-1], poly))
    return result


def evaluated_integral_rhs(
    coefficients: tuple[int, int, int, int, int, int],
) -> list[int]:
    a1, a2, r1, r2, p1, p2 = coefficients
    a_powers = polynomial_powers([-8, a1, a2], 8)
    r_powers = polynomial_powers([2, r1, r2], 8)
    p_powers = polynomial_powers([5, p1, p2], 6)
    result = [0] * 17
    for (a_exponent, r_exponent, p_exponent), coefficient in INTEGRAL_TERMS:
        term = multiply_polynomials(
            multiply_polynomials(
                a_powers[a_exponent],
                r_powers[r_exponent],
            ),
            p_powers[p_exponent],
        )
        for degree, value in enumerate(term):
            result[degree] += coefficient * value
    return result


def rational_square_root(poly: list[int]) -> list[Fraction] | None:
    # INTEGRAL_RHS=256*DOUBLE_COVER_RHS, so its known square root is
    # 16*1984=31744.
    root = [Fraction(31744)]
    assert poly[0] == root[0] ** 2
    for degree in range(1, 9):
        earlier_products = sum(
            root[index] * root[degree - index]
            for index in range(1, degree)
        )
        root.append(
            (Fraction(poly[degree]) - earlier_products) / (2 * root[0])
        )
    for degree in range(9, 17):
        square_coefficient = sum(
            root[index] * root[degree - index]
            for index in range(degree - 8, 9)
        )
        if square_coefficient != poly[degree]:
            return None
    return root


def has_line_image(coefficients: tuple[int, ...]) -> bool:
    a1, a2, r1, r2, p1, p2 = coefficients
    return (
        a1 * r2 - a2 * r1 == 0
        and a1 * p2 - a2 * p1 == 0
        and r1 * p2 - r2 * p1 == 0
    )


def bounded_quadratic_search(bound: int) -> dict[str, object]:
    assert bound > 0
    coefficient_range = range(-bound, bound + 1)
    tested = 0
    degenerate_line_images = 0
    candidates: list[dict[str, object]] = []
    for coefficients in itertools.product(coefficient_range, repeat=6):
        if coefficients[1] == coefficients[3] == coefficients[5] == 0:
            continue
        tested += 1
        if has_line_image(coefficients):
            degenerate_line_images += 1
            continue
        root = rational_square_root(evaluated_integral_rhs(coefficients))
        if root is None:
            continue
        candidates.append(
            {
                "coefficients": list(coefficients),
                "integral_square_root": [str(value) for value in root],
            }
        )

    return {
        "ansatz": [
            "A(t)=-8+a1*t+a2*t^2",
            "R(t)=2+r1*t+r2*t^2",
            "Pi(t)=5+p1*t+p2*t^2",
        ],
        "coefficient_order": ["a1", "a2", "r1", "r2", "p1", "p2"],
        "coefficient_box": [-bound, bound],
        "quadratic_tuples_tested": tested,
        "degenerate_line_images_removed": degenerate_line_images,
        "genuine_quadratic_tuples_tested": tested - degenerate_line_images,
        "square_pullback_candidates": candidates,
        "conclusion": (
            "no candidate in the stated coefficient box"
            if not candidates
            else "candidate pullbacks require exact arithmetic auditing"
        ),
        "scope": "bounded integral-coefficient experiment only",
    }


def local_chart_facts() -> dict[str, object]:
    quadratic_discriminant = sp.factor(A**2 - 4 * FIXED_V)
    assert quadratic_discriminant == -31 * R**2

    cubic_linear = A**2 - FIXED_V - 5
    cubic_constant = (4 * P**3 - FIXED_V * cubic_linear) / A
    x = sp.symbols("x")
    cubic = x**3 - A * x**2 + cubic_linear * x + cubic_constant
    known_cubic = sp.expand(cubic.subs({A: -8, R: 2, P: 5}))
    assert known_cubic == x**3 + 8 * x**2 + 12 * x + 8
    assert int(known_cubic.subs(x, 15)) % 31 == 0
    assert int(sp.diff(known_cubic, x).subs(x, 15)) % 31

    return {
        "quadratic_discriminant_identity": "A^2-4V=-31*R^2",
        "prime_2": (
            "for R nonzero, the quadratic splits over Q_2 because "
            "-31 is a square in Q_2"
        ),
        "prime_31": (
            "the known cubic has a simple root 15 modulo 31; this persists "
            "on a sufficiently small 31-adic parameter neighborhood"
        ),
        "remaining_gap": (
            "new primes ramifying in cubic specializations are not controlled "
            "by these two local observations"
        ),
    }


def build_artifact(bound: int) -> dict[str, object]:
    assert DOUBLE_COVER_RHS.subs({A: -8, R: 2, P: 5}) == 1984**2
    return {
        "schema_version": 1,
        "equation": (
            "Z^2=31*M(A,(A^2+31*R^2)/4,Pi^3)"
        ),
        "base_point": {
            "A": -8,
            "R": 2,
            "Pi": 5,
            "Z": 1984,
        },
        "exact_line_obstruction": exact_line_obstruction(),
        "exact_coordinate_slice_quadratic_obstructions": (
            quadratic_coordinate_slice_obstructions()
        ),
        "quadratic_search": bounded_quadratic_search(bound),
        "local_chart_facts": local_chart_facts(),
        "mathematical_status": {
            "line_ansatz": "exactly excluded",
            "coordinate_fixed_quadratic_ansatz": "exactly excluded",
            "quadratic_ansatz": "bounded experiment",
            "rational_curve": "open",
            "infinitely_many_hasse_failures": "open",
        },
        "reproducing_command": (
            ".venv/bin/python "
            "scripts/search_fixed_quintic_hasse_rational_curves.py"
        ),
    }


def main() -> None:
    args = parse_args()
    if args.quadratic_bound <= 0:
        raise SystemExit("--quadratic-bound must be positive")
    artifact = build_artifact(args.quadratic_bound)
    rendered = json.dumps(artifact, indent=2, sort_keys=True) + "\n"

    if args.emit_json:
        print(rendered, end="")
        return
    if args.write:
        if args.quadratic_bound != 2:
            raise SystemExit("--write is reserved for the default bound 2")
        ARTIFACT.write_text(rendered)
    elif args.quadratic_bound == 2:
        assert ARTIFACT.read_text() == rendered, (
            f"{ARTIFACT.relative_to(ROOT)} is stale; rerun with --write"
        )

    line = artifact["exact_line_obstruction"]
    coordinate_slices = artifact[
        "exact_coordinate_slice_quadratic_obstructions"
    ]
    quadratic = artifact["quadratic_search"]
    print("PASS: fixed-sqrt(-31) double-cover equation and base point are exact")
    print("PASS: all three projective line-direction charts have basis [1]")
    print(f"PASS: {line['conclusion']}")
    print(
        "PASS: all twelve coordinate-slice weighted-projective quadratic "
        f"charts have basis [1] modulo {coordinate_slices['good_prime']}"
    )
    print(f"PASS: {coordinate_slices['conclusion']}")
    print(
        "SEARCH: "
        f"{quadratic['genuine_quadratic_tuples_tested']} genuine quadratic "
        f"tuples in {quadratic['coefficient_box']}, "
        f"{len(quadratic['square_pullback_candidates'])} candidates"
    )
    print("OPEN: rational curves and Hasse-failure infinitude remain open")


if __name__ == "__main__":
    main()
