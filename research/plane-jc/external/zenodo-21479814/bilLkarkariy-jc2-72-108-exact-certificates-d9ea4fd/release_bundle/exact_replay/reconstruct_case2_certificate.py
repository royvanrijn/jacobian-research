#!/usr/bin/env python3
"""Reconstruct an explicit characteristic-zero Nullstellensatz certificate.

The four compact Case-2 residuals have degrees 3, 3, 4, 4.  We search for
multipliers whose products have degree at most 8, as described in the research
report, and solve the resulting Macaulay system directly over Q[u]/(H).
"""

from __future__ import annotations

import pickle
import time
from itertools import product
from pathlib import Path

import case2_exact_generate as c2
import exact_core as ec


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "case2_exact_certificate.pkl"
MACAULAY_DEGREE = 8


def monomials_upto(nvars: int, degree: int):
    out = []
    for exponents in product(range(degree + 1), repeat=nvars):
        if sum(exponents) <= degree:
            out.append(exponents)
    return sorted(out, key=lambda m: (sum(m), m))


def add_monomials(a, b):
    return tuple(x + y for x, y in zip(a, b))


def encode_k(value: ec.K):
    return {
        i: (int(value.p[i].p), int(value.p[i].q))
        for i in range(len(value.p))
        if value.p[i]
    }


def encode_pp(poly: ec.PP):
    return {m: encode_k(c) for m, c in poly.d.items()}


def build_system(generators):
    row_monomials = monomials_upto(3, MACAULAY_DEGREE)
    row_index = {m: i for i, m in enumerate(row_monomials)}
    columns = []
    for generator_index, generator in enumerate(generators):
        degree = max(map(sum, generator.d))
        for multiplier_monomial in monomials_upto(3, MACAULAY_DEGREE - degree):
            columns.append((generator_index, multiplier_monomial))

    rows = [dict() for _ in row_monomials]
    for column_index, (generator_index, multiplier_monomial) in enumerate(columns):
        for monomial, coefficient in generators[generator_index].d.items():
            row = row_index[add_monomials(monomial, multiplier_monomial)]
            rows[row][column_index] = coefficient

    rhs = [ec.K(0) for _ in rows]
    rhs[row_index[(0, 0, 0)]] = ec.K(1)
    return row_monomials, columns, rows, rhs


def solve_sparse(rows, rhs, ncols):
    """Sparse forward elimination and back substitution over ec.K."""
    rows = [dict(row) for row in rows]
    rhs = list(rhs)
    pivot_columns = []
    rank = 0
    started = time.monotonic()

    for column in range(ncols):
        pivot = next((i for i in range(rank, len(rows)) if rows[i].get(column)), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        rhs[rank], rhs[pivot] = rhs[pivot], rhs[rank]
        pivot_value = rows[rank][column]
        pivot_row = rows[rank]

        for i in range(rank + 1, len(rows)):
            value = rows[i].get(column)
            if not value:
                continue
            factor = value / pivot_value
            row = rows[i]
            for j, coefficient in pivot_row.items():
                updated = row.get(j, ec.K(0)) - factor * coefficient
                if updated:
                    row[j] = updated
                else:
                    row.pop(j, None)
            rhs[i] = rhs[i] - factor * rhs[rank]

        pivot_columns.append(column)
        rank += 1
        if rank % 10 == 0:
            elapsed = time.monotonic() - started
            nonzero = sum(len(row) for row in rows[rank:])
            print(f"rank={rank} remaining_nnz={nonzero} elapsed={elapsed:.1f}s", flush=True)
        if rank == len(rows):
            break

    for i in range(rank, len(rows)):
        if not rows[i] and rhs[i]:
            raise RuntimeError(f"inconsistent Macaulay row {i}")

    solution = [ec.K(0) for _ in range(ncols)]
    for i in range(rank - 1, -1, -1):
        column = pivot_columns[i]
        value = rhs[i]
        for j, coefficient in rows[i].items():
            if j != column and solution[j]:
                value -= coefficient * solution[j]
        solution[column] = value / rows[i][column]
    print(f"rank={rank}/{len(rows)} columns={ncols}")
    return solution, pivot_columns


def build_multipliers(generators, columns, solution):
    multipliers = [ec.PP() for _ in generators]
    for (generator_index, monomial), coefficient in zip(columns, solution):
        if coefficient:
            multipliers[generator_index] += ec.PP({monomial + (0, 0, 0): coefficient})
    return multipliers


def verify_identity(generators, multipliers):
    residual = ec.PP.const(-1)
    for generator, multiplier in zip(generators, multipliers):
        residual += generator * multiplier
    if residual:
        raise RuntimeError(f"nonzero certificate residual with {len(residual.d)} terms")


def main():
    generators = [c2.final[i] for i in c2.inds]
    row_monomials, columns, rows, rhs = build_system(generators)
    print(f"Macaulay degree {MACAULAY_DEGREE}: {len(rows)} rows x {len(columns)} columns")
    solution, pivots = solve_sparse(rows, rhs, len(columns))
    multipliers = build_multipliers(generators, columns, solution)
    verify_identity(generators, multipliers)
    payload = {
        "format": "jc2-case2-nullstellensatz-v1",
        "field_minpoly": encode_k(ec.K(ec.MOD + 0)),
        "field_minpoly_coefficients": [
            (int(ec.MOD[i].p), int(ec.MOD[i].q)) for i in range(len(ec.MOD))
        ],
        "generator_indices": tuple(c2.inds),
        "macaulay_degree": MACAULAY_DEGREE,
        "row_monomials": row_monomials,
        "columns": columns,
        "pivot_columns": pivots,
        "generators": [encode_pp(q) for q in generators],
        "multipliers": [encode_pp(q) for q in multipliers],
    }
    OUT.write_bytes(pickle.dumps(payload, protocol=4))
    print("CASE2_EXACT_CERTIFICATE_PASS")
    print(f"wrote {OUT} ({OUT.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
