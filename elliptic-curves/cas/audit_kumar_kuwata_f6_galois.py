#!/usr/bin/env python3
"""Compute the constant-field Galois eigenspaces in Kumar--Kuwata Example 9.1.

The paper's ancillary file gives a geometric rank-18 basis for

    Y^2 = X^3 - 33 X + t^6 + 8/t^6.

It is defined over Q(i, sqrt(2), 3^(1/4)).  This script reconstructs the
ten even and eight odd generators in the ancillary file, identifies their
Galois conjugates exactly by elliptic-curve arithmetic at two good rational
fibres, and computes all quadratic-character eigenspace dimensions.

The lattice enumeration uses PARI/GP's qfminim.  All subsequent arithmetic is
exact SymPy arithmetic in the degree-16 number field; no floating-point point
matching is used.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import platform
import re
import subprocess
from typing import Any, Iterable

import sympy as sp
from sympy import I, Matrix, QQ, Rational, eye, sqrt

from pari_bridge import pari_version


PRIMARY_SOURCE_ARXIV = "https://arxiv.org/abs/1409.2931"
PRIMARY_SOURCE_BUNDLE = "https://export.arxiv.org/e-print/1409.2931"
PRIMARY_SOURCE_BUNDLE_SHA256 = (
    "d3903ffa610826528ee87fa1872a46b7fccb1bcc2b0b2b6240bcdd14d71a73dd"
)
PRIMARY_SOURCE_FILE = "auxfiles/Example9.1.txt"
PRIMARY_SOURCE_FILE_SHA256 = (
    "53848ce6b2205353a03b134bf4086bf2ab4d6671d6134a6858800262b4444b2f"
)
REPRODUCE_COMMAND = (
    ".venv/bin/python elliptic-curves/cas/audit_kumar_kuwata_f6_galois.py "
    "--output artifacts/generated-results/elliptic_kumar_kuwata_f6_galois.json"
)


# Ordering from singular.tex, Example 9.1:
# P3,P4,P5,P6,sigma^2(P6),sigma^2(P5),sigma^2(P4),sigma^2(P3),P7,P8.
P_GRAM = [
    [4, 0, 0, -2, 1, 0, 0, -2, 0, 0],
    [0, 4, 0, -2, 1, 0, -2, 0, 0, 0],
    [0, 0, 4, -2, 1, -2, 0, 0, 0, 0],
    [-2, -2, -2, 4, -2, 1, 1, 1, 0, 0],
    [1, 1, 1, -2, 4, -2, -2, -2, -1, -1],
    [0, 0, -2, 1, -2, 4, 0, 0, 0, 0],
    [0, -2, 0, 1, -2, 0, 4, 0, 0, 2],
    [-2, 0, 0, 1, -2, 0, 0, 4, 2, 0],
    [0, 0, 0, 0, -1, 0, 0, 2, 4, 0],
    [0, 0, 0, 0, -1, 0, 2, 0, 0, 4],
]

# The printed Q1,...,Q8 matrix is E8 with diagonal 2 on the twisted F^(3).
# Pullback to F^(6) doubles the heights, as also seen in the ancillary
# rank-18 height matrix.
Q_GRAM = [
    [4, -2, 0, 0, 0, 0, 0, 0],
    [-2, 4, -2, 0, 0, 0, 0, 0],
    [0, -2, 4, -2, 0, 0, 0, 0],
    [0, 0, -2, 4, -2, 0, 0, 0],
    [0, 0, 0, -2, 4, -2, 0, -2],
    [0, 0, 0, 0, -2, 4, -2, 0],
    [0, 0, 0, 0, 0, -2, 4, 0],
    [0, 0, 0, 0, -2, 0, 0, 4],
]

P_NAMES = ["P3", "P4", "P5", "P6", "P6s", "P5s", "P4s", "P3s", "P7", "P8"]
Q_NAMES = [f"Q{j}" for j in range(1, 9)]


def _pari_matrix(matrix: list[list[int]]) -> str:
    return "[" + ";".join(",".join(str(x) for x in row) for row in matrix) + "]"


def minimal_vectors(gram: list[list[int]], norm: int = 4) -> list[tuple[int, ...]]:
    """Return all vectors of the requested norm, using exact PARI qfminim."""

    script = f"""
G={_pari_matrix(gram)};
r=qfminim(G,{norm});
V=r[3];
for(j=1,matsize(V)[2],print(Vec(V[,j])));
"""
    completed = subprocess.run(
        ["gp", "-fq"],
        input=script,
        text=True,
        capture_output=True,
        check=True,
    )
    representatives: list[tuple[int, ...]] = []
    for line in completed.stdout.splitlines():
        if not line.strip():
            continue
        values = tuple(int(token) for token in re.findall(r"-?\d+", line))
        if len(values) != len(gram):
            raise ValueError(f"unexpected qfminim vector: {line!r}")
        representatives.append(values)
    vectors = representatives + [tuple(-x for x in v) for v in representatives]
    expected_count = int(
        re.search(
            r"^\s*(\d+)",
            subprocess.run(
                ["gp", "-fq"],
                input=f"r=qfminim({_pari_matrix(gram)},{norm});print(r[1]);\n",
                text=True,
                capture_output=True,
                check=True,
            ).stdout,
        ).group(1)
    )
    if len(vectors) != expected_count or len(set(vectors)) != len(vectors):
        raise ValueError("qfminim root count or sign expansion was inconsistent")
    return sorted(vectors)


def _number_field() -> tuple[Any, Any, Any, Any]:
    field = QQ.algebraic_field(I, sqrt(2), 3 ** Rational(1, 4))
    ii = field.from_sympy(I)
    beta = field.from_sympy(sqrt(2))
    root3 = field.from_sympy(3 ** Rational(1, 4))
    return field, ii, beta, root3


Point = tuple[Any, Any] | None


def _sections(t_value: Any, ii: Any, beta: Any, root3: Any, field: Any) -> tuple[list[Point], list[Point]]:
    """Transcription of the P and Q sections in ancillary Example9.1.txt."""

    one = field.one
    two = field.convert(2)
    three = field.convert(3)
    four = field.convert(4)
    half = field.convert(Fraction(1, 2))
    t = field.convert(t_value)
    alpha = root3**2
    gamma = beta * root3
    w = (-one + ii * alpha) * half
    T = t**2
    sp = T + two / T
    sps = w**2 * T + two * w / T
    sm = T - two / T

    def p3(s: Any) -> Point:
        return (-s - 9, 3 * ii * alpha * (s + four))

    def p4(s: Any) -> Point:
        return (-s + 9, -3 * alpha * (s - four))

    def p5(s: Any) -> Point:
        return (
            -s + ii * alpha,
            gamma * alpha * (one + ii) * half * (-s + 2 * ii * alpha),
        )

    def p6(s: Any) -> Point:
        x = w * (
            -s + 3 * w * (one - ii) * half * gamma + 2 * alpha + 3 * ii
        )
        y = (
            ((w - one) * ii * half * gamma + 3 * half * (w - one) * (one + ii)) * s
            + (7 * alpha * ii - 3) * half * gamma
            + 3 * (ii + one) * alpha
        )
        return (x, y)

    p_points = [
        p3(sp),
        p4(sp),
        p5(sp),
        p6(sp),
        p6(sps),
        p5(sps),
        p4(sps),
        p3(sps),
        (3 - sp, -3 * sm),
        (-3 - sp, -3 * ii * sm),
    ]

    q1 = (
        -(w**2) * (11 * T**2 + 4 * w**2) / (2 * T),
        t * (-21 * beta * (2 * w + one) * T / 4),
    )
    q2 = (
        2 * w * (4 - 3 * beta) * T - 3 * (one - 2 * beta) + 4 * w**2 / T,
        t
        * (
            3 * (11 - 8 * beta) * T
            - 12 * w**2 * (3 - 2 * beta)
            + 6 * w * (4 - beta) / T
            + 6 * beta / T**2
        ),
    )
    q3 = (
        (2 * T**2 - 6 * T + one) / T,
        t * (-3 * (T**3 - 4 * T**2 + T - one) / T**2),
    )
    q4 = (
        2 * T - 3 * (one - 2 * beta) + 4 * (4 - 3 * beta) / T,
        t
        * (
            -3 * T
            + 6 * (one - 2 * beta)
            - 12 * (4 - 3 * beta) / T
            + (96 - 66 * beta) / T**2
        ),
    )
    q5 = (
        -w * T + 4 * w**2 * (4 - 3 * beta) / T,
        t * (-3 * w**2 * (one - 2 * beta) + (96 - 66 * beta) / T**2),
    )
    q6 = (
        2 * T + 3 * (one - 2 * beta) + 4 * (4 - 3 * beta) / T,
        t
        * (
            3 * T
            + 6 * (one - 2 * beta)
            + 12 * (4 - 3 * beta) / T
            + (96 - 66 * beta) / T**2
        ),
    )
    q7 = (
        (2 * T**2 + 6 * T + one) / T,
        t * 3 * (T**3 + 4 * T**2 + T + one) / T**2,
    )
    q8 = (-T - field.convert(11) / T, t * 21 * (one + 2 * w) / T**2)
    return p_points, [q1, q2, q3, q4, q5, q6, q7, q8]


def _on_curve(point: Point, t_value: Any, field: Any) -> bool:
    if point is None:
        return True
    x, y = point
    t = field.convert(t_value)
    return y**2 == x**3 - 33 * x + t**6 + field.convert(8) / t**6


def _negate(point: Point) -> Point:
    if point is None:
        return None
    return (point[0], -point[1])


def _add(left: Point, right: Point, field: Any) -> Point:
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2:
        if y1 + y2 == field.zero:
            return None
        slope = (3 * x1**2 - 33) / (2 * y1)
    else:
        slope = (y2 - y1) / (x2 - x1)
    x3 = slope**2 - x1 - x2
    return (x3, slope * (x1 - x3) - y1)


def _multiply(coefficient: int, point: Point, field: Any) -> Point:
    if coefficient < 0:
        return _multiply(-coefficient, _negate(point), field)
    answer: Point = None
    addend = point
    n = coefficient
    while n:
        if n & 1:
            answer = _add(answer, addend, field)
        addend = _add(addend, addend, field)
        n >>= 1
    return answer


def _linear_combination(points: list[Point], vector: Iterable[int], field: Any) -> Point:
    answer: Point = None
    for coefficient, point in zip(vector, points, strict=True):
        if coefficient:
            answer = _add(answer, _multiply(coefficient, point, field), field)
    return answer


def _action_matrix(
    basis_by_fibre: list[list[Point]],
    conjugates_by_fibre: list[list[Point]],
    roots: list[tuple[int, ...]],
    field: Any,
) -> Matrix:
    lookup: dict[tuple[Point, ...], tuple[int, ...]] = {}
    for vector in roots:
        key = tuple(
            _linear_combination(basis, vector, field) for basis in basis_by_fibre
        )
        if key in lookup and lookup[key] != vector:
            raise ValueError("the two-fibre root lookup was not injective")
        lookup[key] = vector
    columns: list[tuple[int, ...]] = []
    for column in range(len(basis_by_fibre[0])):
        key = tuple(points[column] for points in conjugates_by_fibre)
        if key not in lookup:
            raise ValueError(f"Galois conjugate of basis column {column} was not found")
        columns.append(lookup[key])
    return Matrix.hstack(*(Matrix(v) for v in columns))


def _block_diagonal(left: Matrix, right: Matrix) -> Matrix:
    rows = left.rows + right.rows
    columns = left.cols + right.cols
    answer = Matrix.zeros(rows, columns)
    answer[: left.rows, : left.cols] = left
    answer[left.rows :, left.cols :] = right
    return answer


def _fixed_dimension(actions: list[Matrix], signs: tuple[int, ...]) -> int:
    dimension = actions[0].rows
    equations = Matrix.vstack(
        *(action - sign * eye(dimension) for action, sign in zip(actions, signs, strict=True))
    )
    return dimension - equations.rank()


def _matrix_rows(matrix: Matrix) -> list[list[int]]:
    return [[int(matrix[row, column]) for column in range(matrix.cols)] for row in range(matrix.rows)]


def compute_report() -> dict[str, Any]:
    field, ii, beta, root3 = _number_field()
    fibres = [2, 3, 5, 7]
    original = [_sections(t, ii, beta, root3, field) for t in fibres]
    if not all(
        _on_curve(point, t, field)
        for t, pair in zip(fibres, original, strict=True)
        for part in pair
        for point in part
    ):
        raise ValueError("a transcribed source section failed its curve equation")

    transformations = {
        "complex_conjugation": (-ii, beta, root3),
        "quartic_rotation": (ii, beta, ii * root3),
        "sqrt2_flip": (ii, -beta, root3),
    }
    p_roots = minimal_vectors(P_GRAM)
    q_roots = minimal_vectors(Q_GRAM)
    p_basis_by_fibre = [pair[0] for pair in original]
    q_basis_by_fibre = [pair[1] for pair in original]

    p_actions: dict[str, Matrix] = {}
    q_actions: dict[str, Matrix] = {}
    for name, constants in transformations.items():
        conjugate = [_sections(t, *constants, field) for t in fibres]
        p_actions[name] = _action_matrix(
            p_basis_by_fibre, [pair[0] for pair in conjugate], p_roots, field
        )
        q_actions[name] = _action_matrix(
            q_basis_by_fibre, [pair[1] for pair in conjugate], q_roots, field
        )

    generator_order = list(transformations)
    full_actions = [
        _block_diagonal(p_actions[name], q_actions[name]) for name in generator_order
    ]
    j, rotation, beta_flip = full_actions
    identity = eye(18)
    relations = {
        "complex_conjugation_squared": j**2 == identity,
        "quartic_rotation_fourth_power": rotation**4 == identity,
        "sqrt2_flip_squared": beta_flip**2 == identity,
        "dihedral_relation": j * rotation * j == rotation**-1,
        "sqrt2_flip_commutes_with_complex_conjugation": beta_flip * j == j * beta_flip,
        "sqrt2_flip_commutes_with_quartic_rotation": beta_flip * rotation == rotation * beta_flip,
    }
    if not all(relations.values()):
        raise ValueError("reconstructed Galois matrices failed a field relation")

    gram = Matrix.diag(Matrix(P_GRAM), Matrix(Q_GRAM))
    isometries = {
        name: action.T * gram * action == gram
        for name, action in zip(generator_order, full_actions, strict=True)
    }
    if not all(isometries.values()):
        raise ValueError("reconstructed Galois matrix failed the height pairing")

    characters: list[dict[str, Any]] = []
    squareclasses = {
        (1, 1, 1): 1,
        (-1, 1, 1): -1,
        (1, -1, 1): 3,
        (-1, -1, 1): -3,
        (1, 1, -1): 2,
        (-1, 1, -1): -2,
        (1, -1, -1): 6,
        (-1, -1, -1): -6,
    }
    for sign_j in (-1, 1):
        for sign_rotation in (-1, 1):
            for sign_beta in (-1, 1):
                signs = (sign_j, sign_rotation, sign_beta)
                p_dimension = _fixed_dimension(
                    [p_actions[name] for name in generator_order], signs
                )
                q_dimension = _fixed_dimension(
                    [q_actions[name] for name in generator_order], signs
                )
                characters.append(
                    {
                        "signs": dict(zip(generator_order, signs, strict=True)),
                        "p_subspace_rank": p_dimension,
                        "q_subspace_rank": q_dimension,
                        "full_rank": p_dimension + q_dimension,
                        "base_curve_over_Q_t": signs == (1, 1, 1),
                        "quadratic_character_squareclass": squareclasses[signs],
                        "quadratic_subfield": (
                            "Q"
                            if signs == (1, 1, 1)
                            else f"Q(sqrt({squareclasses[signs]}))"
                        ),
                    }
                )
    characters.sort(key=lambda item: (-item["full_rank"], not item["base_curve_over_Q_t"], str(item["signs"])))
    rational_rank = next(item["full_rank"] for item in characters if item["base_curve_over_Q_t"])
    largest_twist = max(item["full_rank"] for item in characters if not item["base_curve_over_Q_t"])

    return {
        "audit_scope": "Kumar--Kuwata Example 9.1 geometric rank-18 basis",
        "source": {
            "title": (
                "Elliptic K3 surfaces associated with the product of two "
                "elliptic curves: Mordell-Weil lattices and their fields of definition"
            ),
            "authors": ["Abhinav Kumar", "Masato Kuwata"],
            "arxiv": PRIMARY_SOURCE_ARXIV,
            "bundle_url": PRIMARY_SOURCE_BUNDLE,
            "bundle_sha256": PRIMARY_SOURCE_BUNDLE_SHA256,
            "ancillary_file": PRIMARY_SOURCE_FILE,
            "ancillary_file_sha256": PRIMARY_SOURCE_FILE_SHA256,
        },
        "curve": "Y^2 = X^3 - 33 X + t^6 + 8/t^6",
        "polynomial_model_after_X_t2_Y_t3_scaling": "Y^2 = X^3 - 33*t^4*X + t^12 + 8",
        "geometric_rank": 18,
        "basis_field": "Q(i, sqrt(2), 3^(1/4))",
        "field_generators": {
            "complex_conjugation": "i -> -i; sqrt(2), 3^(1/4) fixed",
            "quartic_rotation": "3^(1/4) -> i*3^(1/4); i, sqrt(2) fixed",
            "sqrt2_flip": "sqrt(2) -> -sqrt(2); i, 3^(1/4) fixed",
        },
        "exact_method": {
            "matching_fibres": fibres,
            "p_minimal_vector_count": len(p_roots),
            "q_minimal_vector_count": len(q_roots),
            "point_matching": "exact degree-16 number-field arithmetic at all four fibres",
            "floating_point_matching": False,
            "finite_index_note": (
                "The published P(10)+Q(8) lattice has finite index in the "
                "saturated rank-18 Mordell-Weil lattice. Tensoring with Q "
                "therefore preserves every Galois-character dimension."
            ),
        },
        "generator_actions": {
            name: {
                "p_basis_order": P_NAMES,
                "p_matrix_rows": _matrix_rows(p_actions[name]),
                "q_basis_order": Q_NAMES,
                "q_matrix_rows": _matrix_rows(q_actions[name]),
            }
            for name in generator_order
        },
        "checks": {
            "source_point_specialization_count": 18 * len(fibres),
            "all_source_point_specializations_on_curve": True,
            "galois_group_relations": relations,
            "height_isometries": isometries,
        },
        "quadratic_character_eigenspaces": characters,
        "conclusion": {
            "rank_over_Q_t": rational_rank,
            "largest_quadratic_twist_rank_over_Q_t": largest_twist,
            "competitive_rank_threshold": 12,
            "competitive_for_specialization_search": max(rational_rank, largest_twist) >= 12,
            "specialization_search_performed": False,
            "mathematical_status": "exact finite-fibre computational audit",
            "limitation": (
                "The Galois action was identified in the published height lattices "
                "by exact matching at four rational fibres. Group relations and all "
                "height isometries were then checked exactly. This artifact does not "
                "claim a formal symbolic proof of all 54 section identities over K(t)."
            ),
        },
        "reproduce": REPRODUCE_COMMAND,
        "script_sha256": sha256(Path(__file__).resolve().read_bytes()).hexdigest(),
        "software": {
            "python": platform.python_version(),
            "sympy": sp.__version__,
            "pari_gp": pari_version(),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = compute_report()
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(payload, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
        print(f"wrote {args.output} ({sha256(payload.encode()).hexdigest()})")


if __name__ == "__main__":
    main()
