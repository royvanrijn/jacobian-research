#!/usr/bin/env python3
"""Modular interpolation of one degree-eight order-five pivot chart."""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
from pathlib import Path

from sympy.polys.domains import GF, QQ
from sympy.polys.matrices.sdm import sdm_irref

from explore_degree_five_quantum_residue import add, scale
from interpolate_degree_seven_order_five_chart import (
    evaluate_bivariate,
    polynomial_interpolate,
    polynomial_evaluate,
    rational_line_fit,
    serialize_terms,
    singular_groebner,
)
from verify_degree_eight_relative_quantization_obstruction import (
    family_presentation,
)


PIVOT_COLUMNS = (
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17,
    18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
    33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
    48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62,
    63, 64, 65, 66, 67, 68, 71, 72, 73, 74, 75, 76, 77, 78, 79,
    80, 81, 82, 83, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
    96, 97, 98, 99, 100, 103, 107, 108, 111, 112, 114, 115,
)
DENOMINATOR_DEGREE = 50
NUMERATOR_DEGREES = (52, 54, 56)
TARGET_INDICES = (104, 105, 106)


def ambient_data():
    presentation = family_presentation(QQ, QQ.one, QQ.zero)
    monomials = sorted(
        set(presentation["constant"]).union(
            *(set(column) for column in presentation["strong_columns"])
        )
    )
    assert len(monomials) == 269
    return tuple(monomials[:104]), tuple(
        monomials[index] for index in TARGET_INDICES
    )


PIVOT_MONOMIALS, TARGET_MONOMIALS = ambient_data()


def residual_task(task: tuple[int, int, int]):
    prime, sigma_value, tau_value = task
    field = GF(prime)
    presentation = family_presentation(
        field, field(sigma_value), field(tau_value)
    )
    columns = presentation["strong_columns"]
    basis = [columns[index] for index in PIVOT_COLUMNS]
    rows = {}
    for row_index, monomial in enumerate(PIVOT_MONOMIALS):
        row = {
            column_index: coefficient
            for column_index, column in enumerate(basis)
            if (coefficient := column.get(monomial, field.zero))
        }
        coefficient = presentation["constant"].get(monomial, field.zero)
        if coefficient:
            row[104] = coefficient
        rows[row_index] = row
    reduced, pivots, _ = sdm_irref(rows)
    if pivots != list(range(104)):
        return sigma_value, tau_value, None
    coordinates = [
        reduced[index].get(104, field.zero) for index in range(104)
    ]
    residual = dict(presentation["constant"])
    for coefficient, column in zip(coordinates, basis, strict=True):
        residual = add(residual, scale(column, -coefficient))
    return sigma_value, tau_value, [
        int(residual.get(monomial, field.zero)) % prime
        for monomial in TARGET_MONOMIALS
    ]


def reconstruct_line(tau_value, records, prime):
    first_numerator, denominator = rational_line_fit(
        [(sigma, values[0]) for sigma, values in records],
        NUMERATOR_DEGREES[0],
        DENOMINATOR_DEGREE,
        prime,
    )
    numerators = [first_numerator]
    for target_index, degree in enumerate(NUMERATOR_DEGREES[1:], start=1):
        samples = [
            (
                sigma,
                values[target_index]
                * polynomial_evaluate(denominator, sigma, prime)
                % prime,
            )
            for sigma, values in records
        ]
        numerators.append(polynomial_interpolate(samples, degree, prime))
    return tau_value, denominator, numerators


def nested_reconstruction(lines, prime):
    usable = []
    for tau_value, records in sorted(lines.items()):
        try:
            usable.append(reconstruct_line(tau_value, records, prime))
        except ValueError:
            continue
    assert len(usable) >= max(NUMERATOR_DEGREES) + 4

    denominator = {}
    for sigma_degree in range(DENOMINATOR_DEGREE + 1):
        coefficients = polynomial_interpolate(
            [
                (tau_value, line_denominator[sigma_degree])
                for tau_value, line_denominator, _ in usable
            ],
            DENOMINATOR_DEGREE - sigma_degree,
            prime,
        )
        for tau_degree, coefficient in enumerate(coefficients):
            if coefficient:
                denominator[(sigma_degree, tau_degree)] = coefficient

    numerators = []
    for target_index, total_degree in enumerate(NUMERATOR_DEGREES):
        terms = {}
        for sigma_degree in range(total_degree + 1):
            coefficients = polynomial_interpolate(
                [
                    (tau_value, line_numerators[target_index][sigma_degree])
                    for tau_value, _, line_numerators in usable
                ],
                total_degree - sigma_degree,
                prime,
            )
            for tau_degree, coefficient in enumerate(coefficients):
                if coefficient:
                    terms[(sigma_degree, tau_degree)] = coefficient
        numerators.append(terms)
    return usable, denominator, numerators


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--sigma-samples", type=int, default=106)
    parser.add_argument("--tau-lines", type=int, default=60)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.prime <= max(args.sigma_samples, args.tau_lines):
        raise SystemExit("prime must exceed the interpolation grid dimensions")
    tasks = [
        (args.prime, sigma, tau)
        for tau in range(args.tau_lines)
        for sigma in range(args.sigma_samples)
    ]
    if args.jobs == 1:
        records = [residual_task(task) for task in tasks]
    else:
        with mp.get_context("spawn").Pool(args.jobs) as pool:
            records = pool.map(residual_task, tasks, chunksize=3)
    lines = {}
    chart_poles = []
    for sigma, tau, values in records:
        if values is None:
            chart_poles.append((sigma, tau))
        else:
            lines.setdefault(tau, []).append((sigma, values))

    usable, denominator, numerators = nested_reconstruction(lines, args.prime)
    validated = 0
    for sigma, tau, values in records:
        if values is None:
            continue
        denominator_value = evaluate_bivariate(
            denominator, sigma, tau, args.prime
        )
        if denominator_value == 0:
            continue
        for target, terms in zip(values, numerators, strict=True):
            assert evaluate_bivariate(terms, sigma, tau, args.prime) == (
                target * denominator_value % args.prime
            )
        validated += 1

    values, basis, basis_terms, digest = singular_groebner(
        numerators, denominator, args.prime
    )
    certificate = {
        "scope": "modular degree-eight pivot-chart interpolation",
        "prime": args.prime,
        "grid": {
            "sigma_samples": args.sigma_samples,
            "tau_lines": args.tau_lines,
            "points": len(tasks),
            "chart_poles": len(chart_poles),
            "usable_lines": len(usable),
            "validated_points": validated,
        },
        "degree_bounds": {
            "denominator": DENOMINATOR_DEGREE,
            "numerators": list(NUMERATOR_DEGREES),
        },
        "pivot_columns_zero_based": list(PIVOT_COLUMNS),
        "target_monomials": [list(monomial) for monomial in TARGET_MONOMIALS],
        "denominator_terms": serialize_terms(denominator),
        "numerator_terms": [serialize_terms(terms) for terms in numerators],
        "zero_scheme": {
            "raw_dimension": int(values["RAW_DIMENSION"]),
            "raw_vector_space_dimension": int(values["RAW_VDIM"]),
            "dimension": int(values["DIMENSION"]),
            "vector_space_dimension": int(values["VDIM"]),
            "groebner_basis": basis,
            "groebner_basis_terms": [
                serialize_terms(terms) for terms in basis_terms
            ],
        },
        "singular_program_sha256": digest,
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
    print(
        f"PASS: degree-eight chart reconstructed over GF({args.prime}); "
        f"validated={validated}"
    )
    print(
        "ZERO SCHEME: dimension={} vdim={} gb_size={}".format(
            values["DIMENSION"], values["VDIM"], values["GB_SIZE"]
        )
    )
    for index, polynomial in enumerate(basis, start=1):
        print(f"GB[{index}] = {polynomial}")


if __name__ == "__main__":
    main()
