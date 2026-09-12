#!/usr/bin/env python3
"""Probe one modular line for degree-nine pivot-chart degree bounds.

This is a bounded discovery helper.  It fixes the generic 149-column chart,
evaluates selected residual coordinates on one ``tau`` line, and determines
the smallest numerator degree compatible with denominator degree 72.  Every
fit is checked on unused samples.
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
from pathlib import Path

from sympy.polys.domains import GF, QQ
from sympy.polys.matrices.sdm import sdm_irref

from explore_degree_five_quantum_residue import (
    add,
    pi_power,
    poisson,
    scale,
    solve_affine,
    split_pair,
)
from interpolate_degree_seven_order_five_chart import rational_line_fit
from verify_degree_nine_relative_quantization_obstruction import (
    S2_SUPPORT,
    S4_SUPPORT,
    T2_SUPPORT,
    T4_SUPPORT,
    degree_nine_pair,
    family_presentation,
)
from verify_degree_seven_relative_quantization_obstruction import fifth_defect


PIVOT_COLUMNS = (
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18,
    19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,
    35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
    51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66,
    67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82,
    83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 96, 97, 98, 99, 100,
    101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 115,
    116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128,
    130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142,
    143, 144, 145, 146, 147, 150, 154, 158, 159, 161, 162, 163,
)
DENOMINATOR_DEGREE = 72
TARGET_INDICES = tuple(range(149, 161))


def ambient_data():
    presentation = family_presentation(QQ, QQ.one, QQ.zero)
    monomials = sorted(
        set(presentation["constant"]).union(
            *(set(column) for column in presentation["strong_columns"])
        )
    )
    assert len(monomials) == 371
    pivot_monomials = tuple(monomials[:149])
    basis = [presentation["strong_columns"][index] for index in PIVOT_COLUMNS]
    rows = {
        row_index: {
            column_index: coefficient
            for column_index, column in enumerate(basis)
            if (coefficient := column.get(monomial, QQ.zero))
        }
        for row_index, monomial in enumerate(pivot_monomials)
    }
    _, pivots, _ = sdm_irref(rows)
    assert pivots == list(range(149))
    return pivot_monomials, tuple(
        monomials[index] for index in TARGET_INDICES
    )


PIVOT_MONOMIALS, TARGET_MONOMIALS = ambient_data()


def pivot_data(field, sigma, tau):
    """Build only the 164-column prefix needed by the fixed chart."""

    S, T = degree_nine_pair(field, sigma, tau)
    correction_three = [
        poisson({monomial: field.one}, T) for monomial in S2_SUPPORT
    ]
    correction_three += [
        poisson(S, {monomial: field.one}) for monomial in T2_SUPPORT
    ]
    rhs_three = scale(pi_power(S, T, 3), -field.one / field(24))
    particular, kernel, rank_three = solve_affine(
        correction_three, rhs_three, field
    )
    assert rank_three == 227 and len(kernel) == 12
    base_pair = split_pair(particular, S2_SUPPORT, T2_SUPPORT)
    kernel_pairs = [
        split_pair(vector, S2_SUPPORT, T2_SUPPORT) for vector in kernel[:8]
    ]
    constant = fifth_defect(S, T, base_pair, field)
    correction_five = [
        poisson({monomial: field.one}, T) for monomial in S4_SUPPORT
    ]
    correction_five += [
        poisson(S, {monomial: field.one}) for monomial in T4_SUPPORT
    ]
    lower_prefix = []
    for basis_s, basis_t in kernel_pairs:
        diagonal = poisson(basis_s, basis_t)
        shifted = fifth_defect(
            S,
            T,
            (add(base_pair[0], basis_s), add(base_pair[1], basis_t)),
            field,
        )
        linear = add(
            add(shifted, constant, -field.one),
            diagonal,
            -field.one,
        )
        lower_prefix.extend((linear, diagonal))
    strong_prefix = correction_five + lower_prefix
    assert len(strong_prefix) == 164
    return strong_prefix, constant


def residual_task(task):
    prime, sigma_value, tau_value = task
    field = GF(prime)
    columns, constant = pivot_data(
        field, field(sigma_value), field(tau_value)
    )
    basis = [columns[index] for index in PIVOT_COLUMNS]
    rows = {}
    for row_index, monomial in enumerate(PIVOT_MONOMIALS):
        row = {
            column_index: coefficient
            for column_index, column in enumerate(basis)
            if (coefficient := column.get(monomial, field.zero))
        }
        coefficient = constant.get(monomial, field.zero)
        if coefficient:
            row[149] = coefficient
        rows[row_index] = row
    reduced, pivots, _ = sdm_irref(rows)
    if pivots != list(range(149)):
        return sigma_value, None
    coordinates = [reduced[index].get(149, field.zero) for index in range(149)]
    residual = dict(constant)
    for coefficient, column in zip(coordinates, basis, strict=True):
        residual = add(residual, scale(column, -coefficient))
    return sigma_value, [
        int(residual.get(monomial, field.zero)) % prime
        for monomial in TARGET_MONOMIALS
    ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=1009)
    parser.add_argument("--tau", type=int, default=0)
    parser.add_argument("--samples", type=int, default=166)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.prime <= args.samples:
        raise SystemExit("prime must exceed sample count")
    tasks = [(args.prime, sigma, args.tau) for sigma in range(args.samples)]
    if args.jobs == 1:
        records = [residual_task(task) for task in tasks]
    else:
        with mp.get_context("spawn").Pool(args.jobs) as pool:
            records = pool.map(residual_task, tasks, chunksize=2)
    records = [(sigma, values) for sigma, values in records if values is not None]
    results = []
    for target_index in range(len(TARGET_INDICES)):
        samples = [(sigma, values[target_index]) for sigma, values in records]
        fit = None
        for numerator_degree in range(68, 92):
            try:
                numerator, denominator = rational_line_fit(
                    samples, numerator_degree, DENOMINATOR_DEGREE, args.prime
                )
            except ValueError:
                continue
            fit = {
                "target_index": TARGET_INDICES[target_index],
                "numerator_degree_bound": numerator_degree,
                "actual_numerator_degree": max(
                    index for index, coefficient in enumerate(numerator) if coefficient
                ),
                "actual_denominator_degree": max(
                    index for index, coefficient in enumerate(denominator) if coefficient
                ),
            }
            break
        results.append(fit)
    certificate = {
        "scope": "one-line modular degree-bound probe; discovery only",
        "prime": args.prime,
        "tau": args.tau,
        "samples": len(records),
        "denominator_degree_bound": DENOMINATOR_DEGREE,
        "results": results,
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
    print(json.dumps(certificate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
