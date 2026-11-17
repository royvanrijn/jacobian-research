#!/usr/bin/env python3
"""Sparsify and profile the essential 21-variable cubic Keller collision.

The source artifact is already independently certified.  This script applies
five exact elementary linear conjugations, transports the stored collision,
and writes a sparser conjugate.  It also records exact generic-Jacobian and
vector-Waring diagnostics.  Singular is used only for the exact rank over
QQ(x); all conjugation and collision checks use rational arithmetic here.
"""

from __future__ import annotations

from fractions import Fraction as Q
from itertools import combinations_with_replacement
import json
from math import comb
from pathlib import Path
import re
import shutil
import subprocess

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "generated-results" / "essential_bcw_21_counterexample.json"
OUTPUT = ROOT / "artifacts" / "generated-results" / "low_complexity_bcw_21_counterexample.json"

# Each triple (i,j,a) denotes S=I+a E_ij and replaces F by S^-1 F S.
# The fourth move is term-neutral; it exposes the cancellation in the fifth.
SHEARS = (
    (4, 6, Q(-3, 2)),
    (5, 8, Q(-3)),
    (10, 8, Q(1)),
    (4, 3, Q(-3)),
    (13, 7, Q(6, 7)),
)

Exponent = tuple[int, ...]
Poly = dict[Exponent, Q]
Vector = list[Q]


def qtext(value: Q) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def decode_h(stored: dict[str, object]) -> list[Poly]:
    n = int(stored["dimension"])
    result: list[Poly] = []
    for component in stored["H"]:
        poly: Poly = {}
        for term in component:
            exponents = [0] * n
            for index, exponent in term["monomial"]:
                exponents[index] = exponent
            poly[tuple(exponents)] = Q(term["coefficient"])
        result.append(poly)
    return result


def encode_h(components: list[Poly]) -> list[list[dict[str, object]]]:
    return [
        [
            {
                "coefficient": qtext(coefficient),
                "monomial": [[i, exponent] for i, exponent in enumerate(exponents) if exponent],
            }
            for exponents, coefficient in sorted(poly.items(), reverse=True)
        ]
        for poly in components
    ]


def shear_conjugate(components: list[Poly], i: int, j: int, a: Q) -> list[Poly]:
    """Return S^-1 H(Sx) for S=I+a E_ij."""

    substituted: list[Poly] = []
    for poly in components:
        answer: Poly = {}
        for exponents, coefficient in poly.items():
            power = exponents[i]
            for moved in range(power + 1):
                new = list(exponents)
                new[i] -= moved
                new[j] += moved
                key = tuple(new)
                answer[key] = answer.get(key, Q(0)) + coefficient * comb(power, moved) * a**moved
        substituted.append({key: value for key, value in answer.items() if value})

    changed = substituted[i].copy()
    for exponents, coefficient in substituted[j].items():
        changed[exponents] = changed.get(exponents, Q(0)) - a * coefficient
    substituted[i] = {key: value for key, value in changed.items() if value}
    return substituted


def evaluate(poly: Poly, point: Vector) -> Q:
    return sum(
        coefficient
        * product(point[index] ** exponent for index, exponent in enumerate(exponents) if exponent)
        for exponents, coefficient in poly.items()
    )


def product(values) -> Q:
    answer = Q(1)
    for value in values:
        answer *= value
    return answer


def transport(vector: Vector, i: int, j: int, a: Q) -> Vector:
    """Apply S^-1 to a point, for S=I+a E_ij."""

    answer = vector.copy()
    answer[i] -= a * answer[j]
    return answer


def matrix_rank(matrix: list[list[Q]]) -> int:
    rows = [row.copy() for row in matrix]
    if not rows:
        return 0
    m, n = len(rows), len(rows[0])
    rank = 0
    for column in range(n):
        pivot = next((row for row in range(rank, m) if rows[row][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        scale = rows[rank][column]
        rows[rank] = [value / scale for value in rows[rank]]
        for row in range(rank + 1, m):
            if rows[row][column]:
                scale = rows[row][column]
                rows[row] = [left - scale * right for left, right in zip(rows[row], rows[rank])]
        rank += 1
    return rank


def matrix_product(left: list[list[Q]], right: list[list[Q]]) -> list[list[Q]]:
    n = len(left)
    answer = [[Q(0)] * n for _ in range(n)]
    nonzero_right = [[(j, value) for j, value in enumerate(row) if value] for row in right]
    for i, row in enumerate(left):
        for k, value in enumerate(row):
            if value:
                for j, right_value in nonzero_right[k]:
                    answer[i][j] += value * right_value
    return answer


def evaluated_jacobian(components: list[Poly], point: Vector) -> list[list[Q]]:
    n = len(components)
    answer = [[Q(0)] * n for _ in range(n)]
    for output, poly in enumerate(components):
        for exponents, coefficient in poly.items():
            for variable, exponent in enumerate(exponents):
                if not exponent:
                    continue
                term = coefficient * exponent
                for index, power in enumerate(exponents):
                    term *= point[index] ** (power - (index == variable))
                answer[output][variable] += term
    return answer


def exact_point_power_ranks(components: list[Poly]) -> list[int]:
    point = [Q(i * i + 3 * i + 5) for i in range(len(components))]
    jacobian = evaluated_jacobian(components, point)
    power = jacobian
    ranks: list[int] = []
    for _ in range(len(components)):
        ranks.append(matrix_rank(power))
        if ranks[-1] == 0:
            break
        power = matrix_product(power, jacobian)
    return ranks


def singular_generic_rank(components: list[Poly]) -> int:
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required for the exact QQ(x) rank audit")
    n = len(components)
    entries: list[str] = []
    for poly in components:
        for variable in range(n):
            terms: list[str] = []
            for exponents, coefficient in poly.items():
                if not exponents[variable]:
                    continue
                factors = [f"({qtext(coefficient)})", str(exponents[variable])]
                for index, exponent in enumerate(exponents):
                    exponent -= index == variable
                    if exponent:
                        factors.append(f"x{index}^{exponent}")
                terms.append("*".join(factors))
            entries.append("+".join(terms) or "0")
    program = (
        "ring r=0,(" + ",".join(f"x{i}" for i in range(n)) + "),dp;\n"
        + f"matrix J[{n}][{n}]=" + ",".join(entries) + ";\n"
        + "print(rank(J));\nquit;\n"
    )
    result = subprocess.run(
        [singular, "-q"], input=program, text=True, capture_output=True, check=True
    )
    integers = re.findall(r"(?m)^\s*(\d+)\s*$", result.stdout)
    if not integers:
        raise RuntimeError(f"could not parse Singular rank output: {result.stdout!r}")
    return int(integers[-1])


def flattening_ranks(components: list[Poly]) -> tuple[int, int]:
    n = len(components)
    cubic_monomials = sorted({exponents for poly in components for exponents in poly})
    component_matrix = sp.Matrix(
        [[sp.Rational(poly.get(exponents, 0).numerator, poly.get(exponents, 0).denominator)
          for exponents in cubic_monomials] for poly in components]
    )

    derivatives: list[dict[Exponent, Q]] = []
    for poly in components:
        for variable in range(n):
            derivative: Poly = {}
            for exponents, coefficient in poly.items():
                if exponents[variable]:
                    reduced = list(exponents)
                    reduced[variable] -= 1
                    key = tuple(reduced)
                    derivative[key] = derivative.get(key, Q(0)) + coefficient * exponents[variable]
            derivatives.append(derivative)
    quadratic_monomials = sorted({exponents for poly in derivatives for exponents in poly})
    jacobian_flattening = sp.Matrix(
        [[sp.Rational(poly.get(exponents, 0).numerator, poly.get(exponents, 0).denominator)
          for exponents in quadratic_monomials] for poly in derivatives]
    )
    return component_matrix.rank(), jacobian_flattening.rank()


def add_form(
    decomposition: dict[tuple[Q, ...], Vector], form: Vector, vector: Vector
) -> None:
    first = next(value for value in form if value)
    if first < 0:
        form = [-value for value in form]
        vector = [-value for value in vector]
    key = tuple(form)
    current = decomposition.setdefault(key, [Q(0)] * len(vector))
    decomposition[key] = [left + right for left, right in zip(current, vector)]


def waring_decomposition(components: list[Poly]) -> dict[tuple[Q, ...], Vector]:
    """Construct H=sum v_r ell_r^3 by polarizing its sparse monomials."""

    n = len(components)
    coefficient_vectors: dict[Exponent, Vector] = {}
    for output, poly in enumerate(components):
        for exponents, coefficient in poly.items():
            coefficient_vectors.setdefault(exponents, [Q(0)] * n)[output] = coefficient

    decomposition: dict[tuple[Q, ...], Vector] = {}
    for exponents, vector in coefficient_vectors.items():
        support = [(i, exponent) for i, exponent in enumerate(exponents) if exponent]
        if len(support) == 1:
            form = [Q(0)] * n
            form[support[0][0]] = 1
            add_form(decomposition, form, vector)
        elif len(support) == 2:
            squared = next(i for i, exponent in support if exponent == 2)
            linear = next(i for i, exponent in support if exponent == 1)
            plus = [Q(0)] * n
            minus = [Q(0)] * n
            pure = [Q(0)] * n
            plus[squared] = plus[linear] = 1
            minus[squared], minus[linear] = 1, -1
            pure[linear] = 1
            add_form(decomposition, plus, [value / 6 for value in vector])
            add_form(decomposition, minus, [-value / 6 for value in vector])
            add_form(decomposition, pure, [-value / 3 for value in vector])
        elif len(support) == 3:
            indices = [i for i, _ in support]
            for signs, scale in (
                ((1, 1, 1), Q(1, 24)),
                ((1, 1, -1), Q(-1, 24)),
                ((1, -1, 1), Q(-1, 24)),
                ((-1, 1, 1), Q(-1, 24)),
            ):
                form = [Q(0)] * n
                for index, sign in zip(indices, signs):
                    form[index] = sign
                add_form(decomposition, form, [scale * value for value in vector])
        else:
            raise AssertionError(f"unexpected cubic monomial support: {support}")

    decomposition = {
        form: vector for form, vector in decomposition.items() if any(vector)
    }
    reconstructed = [dict() for _ in range(n)]
    for form, vector in decomposition.items():
        support = [i for i, value in enumerate(form) if value]
        for i, j, k in combinations_with_replacement(support, 3):
            multiplicity = 1 if i == k else (3 if i == j or j == k else 6)
            exponents = [0] * n
            exponents[i] += 1
            exponents[j] += 1
            exponents[k] += 1
            scalar = multiplicity * form[i] * form[j] * form[k]
            for output, value in enumerate(vector):
                if value:
                    key = tuple(exponents)
                    reconstructed[output][key] = reconstructed[output].get(key, Q(0)) + value * scalar
    reconstructed = [
        {key: value for key, value in poly.items() if value} for poly in reconstructed
    ]
    assert reconstructed == components
    return decomposition


def main() -> None:
    stored = json.loads(SOURCE.read_text())
    assert stored["dimension"] == 21
    assert stored["jacobian_determinant"] == "1"
    original = decode_h(stored)
    components = original
    points = [[Q(value) for value in point] for point in stored["collision_points"]]
    common_image = [Q(value) for value in stored["common_image"]]

    trace = []
    for i, j, a in SHEARS:
        before = sum(len(poly) for poly in components)
        components = shear_conjugate(components, i, j, a)
        points = [transport(point, i, j, a) for point in points]
        common_image = transport(common_image, i, j, a)
        trace.append(
            {
                "i": i,
                "j": j,
                "a": qtext(a),
                "terms_before": before,
                "terms_after": sum(len(poly) for poly in components),
            }
        )

    assert sum(len(poly) for poly in original) == 60
    assert sum(len(poly) for poly in components) == 53
    assert all(sum(exponents) == 3 for poly in components for exponents in poly)
    assert len({tuple(point) for point in points}) == 3
    for point in points:
        image = [coordinate + evaluate(poly, point) for coordinate, poly in zip(point, components)]
        assert image == common_image

    distinct_monomials = len({exponents for poly in components for exponents in poly})
    generic_rank = singular_generic_rank(components)
    power_ranks = exact_point_power_ranks(components)
    component_rank, derivative_flattening_rank = flattening_ranks(components)
    decomposition = waring_decomposition(components)

    assert distinct_monomials == 45
    assert generic_rank == 18
    assert power_ranks == list(range(18, -1, -1))
    assert component_rank == 20
    assert derivative_flattening_rank == 42
    assert len(decomposition) == 164

    artifact = {
        "format": "low-complexity-essential-bcw-sparse-cubic-homogeneous-map-v1",
        "source": str(SOURCE.relative_to(ROOT)),
        "dimension": 21,
        "linear_part": "identity",
        "jacobian_determinant": "1",
        "jacobian_certificate": "linear conjugate of the certified essential BCW 21 map",
        "conjugation_shears": trace,
        "H": encode_h(components),
        "collision_points": [[qtext(value) for value in point] for point in points],
        "common_image": [qtext(value) for value in common_image],
        "statistics": {
            "component_monomial_terms": 53,
            "distinct_scalar_monomials": distinct_monomials,
            "nonzero_components": sum(bool(poly) for poly in components),
            "generic_rank_JH_over_QQ_x": generic_rank,
            "nilpotency_index_JH": 19,
            "exact_point_power_ranks": power_ranks,
            "component_flattening_rank": component_rank,
            "derivative_flattening_rank": derivative_flattening_rank,
            "vector_waring_rank_lower_bound": derivative_flattening_rank,
            "vector_waring_rank_constructive_upper_bound": len(decomposition),
        },
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")

    print("PASS low-complexity BCW 21: exact linear conjugation reduces 60 terms to 53")
    print("PASS low-complexity BCW 21: transported three-point collision remains exact")
    print("PASS low-complexity BCW 21: generic rank(JH)=18 and nilpotency index=19")
    print("PASS low-complexity BCW 21: vector-Waring rank is between 42 and 164")
    print(f"PASS low-complexity BCW 21: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
