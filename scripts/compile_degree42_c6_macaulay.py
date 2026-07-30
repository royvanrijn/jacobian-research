#!/usr/bin/env python3
"""Compile the characteristic-zero degree-42 ``c6`` certificate.

This replaces a global support-saturation calculation by three fixed sparse
linear systems:

* ``w0*c6`` belongs to ``J6``;
* ``w2*c6`` belongs to ``J6``;
* a finite-support Macaulay dual functional kills ``J6`` but not ``c6``.

The modular solver uses many 31-bit primes and a two-sided block-Wiedemann
projection.  Rational reconstruction is never accepted without exact replay.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Sequence

import sympy as sp
import flint
from flint import nmod_mat


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

from jcsearch.sparse_macaulay import (  # noqa: E402
    IntegerPolynomial,
    RationalPolynomial,
    add_polynomial,
    annihilator_multipliers,
    block_wiedemann_solve,
    canonical_digest,
    crt_pair,
    deserialize_polynomial,
    monomials,
    multiply_polynomial,
    primitive_integer_polynomial,
    reconstruct_vector,
    serialize_polynomial,
)
from verify_degree42_order7_known_witness import C6  # noqa: E402
from verify_degree42_transported_27_normal_jets import (  # noqa: E402
    transformed_problem,
)


DEFAULT_OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree42_c6_macaulay_certificate.json"
)
VARIABLES = ("u", "v", "sync42p_w0", "sync42p_w1", "sync42p_w2")
NORMAL_CUTOFF = 6


def descending_31_bit_primes(count: int) -> tuple[int, ...]:
    """Return a deterministic schedule beginning at the Mersenne prime."""

    output = []
    cursor = 2**31
    for _index in range(count):
        cursor = int(sp.prevprime(cursor))
        output.append(cursor)
    return tuple(output)


PRIMES = descending_31_bit_primes(96)


def sympy_polynomial(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> RationalPolynomial:
    polynomial = sp.Poly(sp.cancel(expression), *variables, domain=sp.QQ)
    return {
        tuple(int(value) for value in monomial): Fraction(
            int(coefficient.p),
            int(coefficient.q),
        )
        for monomial, coefficient in polynomial.terms()
        if coefficient and sum(monomial[:2]) < NORMAL_CUTOFF
    }


def build_problem() -> dict[str, Any]:
    """Derive and normalize the fixed fiber directly from the Ritt model."""

    normals, bases, residuals, _defect = transformed_problem()
    pivot_three = residuals[4].subs(normals[2], 0)
    after_three = [
        item.subs(normals[2], pivot_three) for item in residuals
    ]
    pivot_four = after_three[10].subs(normals[3], 0)
    after_four = [
        item.subs(normals[3], pivot_four) for item in after_three
    ]
    pivot_five = after_four[16].subs(normals[4], 0)
    reduced = [
        item.subs(normals[4], pivot_five) for item in after_four
    ]
    u, v = sp.symbols("u v")
    variables = (u, v) + bases[3:]
    specialization = {
        normals[0]: u,
        normals[1]: v,
        bases[0]: 1,
        bases[1]: 2,
        bases[2]: 3,
    }
    generators: list[IntegerPolynomial] = []
    source_scales: list[list[int]] = []
    source_indices: list[int] = []
    for index, expression in enumerate(reduced):
        polynomial = sympy_polynomial(
            expression.subs(specialization, simultaneous=True),
            variables,
        )
        if not polynomial:
            continue
        primitive, scale = primitive_integer_polynomial(polynomial)
        generators.append(primitive)
        source_indices.append(index)
        source_scales.append([scale.numerator, scale.denominator])
    locals_ = {str(variable): variable for variable in variables}
    c6_expression = sp.sympify(C6.replace("^", "**"), locals=locals_)
    c6 = sympy_polynomial(c6_expression, variables)
    assert all(coefficient.denominator == 1 for coefficient in c6.values())
    problem = {
        "variables": list(VARIABLES),
        "normal_variables": ["u", "v"],
        "normal_cutoff": NORMAL_CUTOFF,
        "base_specialization": {
            "sync42p_e1": 1,
            "sync42p_e2": 2,
            "sync42p_t": 3,
        },
        "source_residual_indices_zero_based": source_indices,
        "source_scales": source_scales,
        "generators": [
            serialize_polynomial(generator) for generator in generators
        ],
        "c6": serialize_polynomial(c6),
    }
    problem["sha256"] = canonical_digest(problem)
    return problem


def integer_polynomial(data: Sequence[Sequence[object]]) -> IntegerPolynomial:
    decoded = deserialize_polynomial(data)
    assert all(value.denominator == 1 for value in decoded.values())
    return {
        monomial: value.numerator for monomial, value in decoded.items()
    }


def membership_columns(
    generators: list[IntegerPolynomial],
    normal_multiplier_degree: int,
    boundary_multiplier_degree: int,
) -> tuple[
    list[tuple[int, tuple[int, ...]]],
    list[dict[tuple[int, ...], int]],
]:
    normal = monomials(2, normal_multiplier_degree)
    boundary = monomials(3, boundary_multiplier_degree)
    labels: list[tuple[int, tuple[int, ...]]] = []
    columns: list[dict[tuple[int, ...], int]] = []
    generator_order = sorted(
        range(len(generators)),
        key=lambda index: (
            len(generators[index]),
            max(
                (sum(monomial) for monomial in generators[index]),
                default=0,
            ),
            index,
        ),
    )
    for generator_index in generator_order:
        generator = generators[generator_index]
        for normal_monomial in normal:
            for boundary_monomial in boundary:
                multiplier = normal_monomial + boundary_monomial
                product = multiply_polynomial(
                    {multiplier: Fraction(1)},
                    generator,
                    normal_variables=2,
                    normal_cutoff=NORMAL_CUTOFF,
                )
                if not product:
                    continue
                labels.append((generator_index, multiplier))
                columns.append(
                    {
                        monomial: coefficient.numerator
                        for monomial, coefficient in product.items()
                    }
                )
    return labels, columns


def columns_to_rows(
    columns: list[dict[tuple[int, ...], int]],
    targets: list[RationalPolynomial],
) -> tuple[
    list[tuple[int, ...]],
    list[dict[int, int]],
    list[list[int]],
]:
    row_labels = sorted(
        {
            monomial
            for column in columns
            for monomial in column
        }
        | {
            monomial
            for target in targets
            for monomial in target
        },
        key=lambda item: (sum(item), item),
    )
    row_index = {
        monomial: index for index, monomial in enumerate(row_labels)
    }
    rows = [{} for _label in row_labels]
    for column_index, column in enumerate(columns):
        for monomial, coefficient in column.items():
            rows[row_index[monomial]][column_index] = coefficient
    right_sides = [
        [
            int(target.get(monomial, Fraction(0)))
            for monomial in row_labels
        ]
        for target in targets
    ]
    return row_labels, rows, right_sides


def rank_profile(
    rows: list[dict[int, int]],
    right_sides: list[list[int]],
    prime: int,
) -> tuple[list[int], list[int]]:
    """Select basis columns and nonsingular pivot rows over one good prime."""

    row_count = len(rows)
    column_count = 1 + max(
        (column for row in rows for column in row),
        default=-1,
    )
    dense = [
        [row.get(column, 0) % prime for column in range(column_count)]
        + [right_side[row_index] % prime for right_side in right_sides]
        for row_index, row in enumerate(rows)
    ]
    reduced, rank = nmod_mat(dense, prime).rref()
    pivots: list[int] = []
    for row_index in range(rank):
        pivot = next(
            column
            for column in range(column_count + len(right_sides))
            if int(reduced[row_index, column])
        )
        pivots.append(pivot)
    if any(pivot >= column_count for pivot in pivots):
        raise RuntimeError("a target is outside the Macaulay column space")
    basis_columns = pivots
    # A full-column-rank matrix has a nonsingular row minor.  Find its row
    # profile by rref of the transpose.
    transposed = [
        [
            rows[row_index].get(column, 0) % prime
            for row_index in range(row_count)
        ]
        for column in basis_columns
    ]
    transposed_reduced, transposed_rank = nmod_mat(
        transposed,
        prime,
    ).rref()
    if transposed_rank != len(basis_columns):
        raise RuntimeError("rank-profile construction lost column rank")
    pivot_rows = []
    for row_index in range(transposed_rank):
        pivot_rows.append(
            next(
                column
                for column in range(row_count)
                if int(transposed_reduced[row_index, column])
            )
        )
    return basis_columns, pivot_rows


def restricted_square(
    rows: list[dict[int, int]],
    right_sides: list[list[int]],
    basis_columns: list[int],
    pivot_rows: list[int],
) -> tuple[list[dict[int, int]], list[list[int]]]:
    basis_lookup = {
        column: index for index, column in enumerate(basis_columns)
    }
    square = [
        {
            basis_lookup[column]: coefficient
            for column, coefficient in rows[row_index].items()
            if column in basis_lookup
        }
        for row_index in pivot_rows
    ]
    restricted_right = [
        [right_side[row_index] for row_index in pivot_rows]
        for right_side in right_sides
    ]
    return square, restricted_right


def exact_membership_replay(
    generators: list[IntegerPolynomial],
    labels: list[tuple[int, tuple[int, ...]]],
    basis_columns: list[int],
    coefficients: list[Fraction],
    target: RationalPolynomial,
) -> bool:
    residual = {monomial: -value for monomial, value in target.items()}
    for coordinate, column_index in zip(
        coefficients,
        basis_columns,
        strict=True,
    ):
        if not coordinate:
            continue
        generator_index, multiplier = labels[column_index]
        product = multiply_polynomial(
            {multiplier: Fraction(1)},
            generators[generator_index],
            normal_variables=2,
            normal_cutoff=NORMAL_CUTOFF,
        )
        add_polynomial(residual, product, coordinate)
    return not residual


def modular_reconstruct(
    square: list[dict[int, int]],
    right_sides: list[list[int]],
    primes: tuple[int, ...],
    block_size: int,
    seed: int,
    replay,
    minimum_primes: int = 8,
) -> tuple[list[list[Fraction]], list[dict[str, object]], int]:
    residues = [[0] * len(square) for _right in right_sides]
    modulus = 1
    transcripts: list[dict[str, object]] = []
    for prime_index, prime in enumerate(primes):
        attempt = 0
        while True:
            prime_seed = seed + 1009 * prime_index + attempt
            try:
                solutions, transcript = block_wiedemann_solve(
                    square,
                    right_sides,
                    prime,
                    block_size,
                    prime_seed,
                )
                break
            except ZeroDivisionError:
                attempt += 1
                if attempt == 8:
                    raise RuntimeError(
                        f"eight unlucky projections modulo {prime}"
                    )
        transcripts.append(
            transcript
            | {
                "prime": prime,
                "solution_sha256": canonical_digest(solutions),
            }
        )
        for right_index, solution in enumerate(solutions):
            residues[right_index] = [
                crt_pair(old, modulus, new, prime)
                for old, new in zip(
                    residues[right_index],
                    solution,
                    strict=True,
                )
            ]
        modulus *= prime
        if len(transcripts) < minimum_primes:
            continue
        reconstructed = [
            reconstruct_vector(vector, modulus) for vector in residues
        ]
        if all(vector is not None for vector in reconstructed):
            exact = [vector for vector in reconstructed if vector is not None]
            if replay(exact):
                return exact, transcripts, modulus
    raise RuntimeError("rational reconstruction did not stabilize")


def encode_multiplier_certificate(
    labels: list[tuple[int, tuple[int, ...]]],
    basis_columns: list[int],
    coefficients: list[Fraction],
) -> list[dict[str, object]]:
    output = []
    for column_index, coefficient in zip(
        basis_columns,
        coefficients,
        strict=True,
    ):
        if not coefficient:
            continue
        generator_index, monomial = labels[column_index]
        output.append(
            {
                "generator": generator_index,
                "monomial": list(monomial),
                "coefficient": (
                    coefficient.numerator
                    if coefficient.denominator == 1
                    else [coefficient.numerator, coefficient.denominator]
                ),
            }
        )
    return output


def dual_system(
    generators: list[IntegerPolynomial],
    c6: RationalPolynomial,
    boundary_degree: int,
) -> tuple[
    list[tuple[int, ...]],
    list[dict[int, int]],
    list[list[int]],
]:
    row_support = [
        normal + boundary
        for normal in monomials(2, NORMAL_CUTOFF - 1)
        for boundary in monomials(3, boundary_degree)
    ]
    support_index = {
        monomial: index for index, monomial in enumerate(row_support)
    }
    equations: list[dict[int, int]] = []
    for generator in generators:
        for multiplier in annihilator_multipliers(
            row_support,
            generator,
        ):
            product = multiply_polynomial(
                {multiplier: Fraction(1)},
                generator,
            )
            equation = {
                support_index[monomial]: coefficient.numerator
                for monomial, coefficient in product.items()
                if monomial in support_index and coefficient
            }
            if equation:
                equations.append(equation)
    normalization = {
        support_index[monomial]: coefficient.numerator
        for monomial, coefficient in c6.items()
    }
    equations.append(normalization)
    right_side = [0] * (len(equations) - 1) + [1]
    return row_support, equations, [right_side]


def compile_certificate(
    problem: dict[str, Any],
    block_size: int,
) -> dict[str, Any]:
    generators = [
        integer_polynomial(data) for data in problem["generators"]
    ]
    c6 = deserialize_polynomial(problem["c6"])
    w0_c6 = {
        monomial[:2]
        + (monomial[2] + 1,)
        + monomial[3:]: coefficient
        for monomial, coefficient in c6.items()
    }
    w2_c6 = {
        monomial[:-1] + (monomial[-1] + 1,): coefficient
        for monomial, coefficient in c6.items()
    }
    labels, columns = membership_columns(generators, 4, 2)
    _row_labels, rows, right_sides = columns_to_rows(
        columns,
        [w0_c6, w2_c6],
    )
    basis_columns, pivot_rows = rank_profile(
        rows,
        right_sides,
        PRIMES[0],
    )
    square, restricted_right = restricted_square(
        rows,
        right_sides,
        basis_columns,
        pivot_rows,
    )

    def membership_replay(vectors: list[list[Fraction]]) -> bool:
        return all(
            exact_membership_replay(
                generators,
                labels,
                basis_columns,
                vector,
                target,
            )
            for vector, target in zip(
                vectors,
                (w0_c6, w2_c6),
                strict=True,
            )
        )

    membership_vectors, membership_runs, membership_modulus = (
        modular_reconstruct(
            square,
            restricted_right,
            PRIMES,
            block_size,
            420006,
            membership_replay,
        )
    )

    dual_support, dual_rows, dual_right = dual_system(
        generators,
        c6,
        2,
    )
    dual_basis, dual_pivots = rank_profile(
        dual_rows,
        dual_right,
        PRIMES[0],
    )
    dual_square, dual_restricted_right = restricted_square(
        dual_rows,
        dual_right,
        dual_basis,
        dual_pivots,
    )

    def dual_replay(vectors: list[list[Fraction]]) -> bool:
        vector = vectors[0]
        full = [Fraction(0)] * len(dual_support)
        for coordinate, column in zip(
            vector,
            dual_basis,
            strict=True,
        ):
            full[column] = coordinate
        return all(
            sum(
                Fraction(coefficient) * full[column]
                for column, coefficient in row.items()
            )
            == Fraction(right)
            for row, right in zip(
                dual_rows,
                dual_right[0],
                strict=True,
            )
        )

    dual_vectors, dual_runs, dual_modulus = modular_reconstruct(
        dual_square,
        dual_restricted_right,
        PRIMES,
        block_size,
        420106,
        dual_replay,
    )
    functional_vector = [Fraction(0)] * len(dual_support)
    for coordinate, column in zip(
        dual_vectors[0],
        dual_basis,
        strict=True,
    ):
        functional_vector[column] = coordinate
    functional = {
        monomial: coefficient
        for monomial, coefficient in zip(
            dual_support,
            functional_vector,
            strict=True,
        )
        if coefficient
    }
    certificate = {
        "schema": "degree42-c6-macaulay-certificate.v1",
        "mathematical_scope": (
            "Exact characteristic-zero statements on the specialized "
            "(e1,e2,t)=(1,2,3) sixth normal jet: c6 is nonzero modulo J6 "
            "and is annihilated there by w0 and w2. This is not a global "
            "saturation computation, an order-seven lift, or an inverse "
            "calculation."
        ),
        "software_assumptions": {
            "python": platform.python_version(),
            "sympy": sp.__version__,
            "python_flint": flint.__version__,
        },
        "problem": problem,
        "macaulay": {
            "membership": {
                "rows": len(rows),
                "columns": len(columns),
                "rank": len(basis_columns),
                "normal_multiplier_degree": 4,
                "boundary_multiplier_degree": 2,
            },
            "dual": {
                "equations": len(dual_rows),
                "functional_coordinates": len(dual_support),
                "rank": len(dual_basis),
                "boundary_support_degree": 2,
            },
        },
        "modular": {
            "membership_runs": membership_runs,
            "membership_crt_modulus": str(membership_modulus),
            "dual_runs": dual_runs,
            "dual_crt_modulus": str(dual_modulus),
        },
        "annihilation_certificates": {
            "w0": encode_multiplier_certificate(
                labels,
                basis_columns,
                membership_vectors[0],
            ),
            "w2": encode_multiplier_certificate(
                labels,
                basis_columns,
                membership_vectors[1],
            ),
        },
        "nonmembership_functional": serialize_polynomial(functional),
    }
    certificate["certificate_sha256"] = canonical_digest(certificate)
    return certificate


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--block-size", type=int, default=8)
    parser.add_argument(
        "--problem-only",
        action="store_true",
        help="derive and print the normalized problem without solving",
    )
    args = parser.parse_args()
    problem = build_problem()
    if args.problem_only:
        print(json.dumps(problem, indent=2, sort_keys=True))
        return
    certificate = compile_certificate(problem, args.block_size)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    )
    print(
        json.dumps(
            {
                "certificate": str(args.output.relative_to(ROOT)),
                "certificate_sha256": certificate["certificate_sha256"],
                "problem_sha256": problem["sha256"],
            },
            sort_keys=True,
        )
    )
    print("PASS: reconstructed exact degree-42 c6 Macaulay certificate")


if __name__ == "__main__":
    main()
