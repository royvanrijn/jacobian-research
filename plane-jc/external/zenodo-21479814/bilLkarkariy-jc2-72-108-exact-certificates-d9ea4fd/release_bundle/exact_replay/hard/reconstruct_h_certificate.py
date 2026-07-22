#!/usr/bin/env python3
"""Reconstruct the fixed-support exact certificate h = sum(T_i E_i).

This uses the deterministic 1925-row minor selected modulo 71, rebuilds the
same matrix over Q, clears denominators row by row, and delegates the exact
integer solve to FLINT.  The final solution is checked against all 2010 scalar
coefficient equations before the text certificate is written.
"""

from __future__ import annotations

import math
import pickle
import re
import sys
import time
from pathlib import Path

import numpy as np
from flint import ctx, fmpq, fmpz_mat


ROOT = Path(__file__).resolve().parents[1]
HARD = ROOT / "hard"
SYSTEM = ROOT / "hne0_polred.pkl"
MODULAR_SUPPORT = ROOT / "hcert71.out"
MODULAR_MATRIX = ROOT / "fixed_matrix_p71.npz"
PIVOT_ROWS = ROOT / "pivot_scalar_rows_p71.npy"
SOLUTION_CACHE = HARD / "h_solution_exact.pkl"
CERTIFICATE = HARD / "h_certificate_exact.txt"
P = 71
MINPOLY = (26, 0, 3, 3, -1, 1)

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


def q(n=0, d=1):
    return fmpq(int(n), int(d))


def add_monomials(a, b):
    return tuple(x + y for x, y in zip(a, b))


def reduce_poly(coefficients):
    values = list(coefficients)
    if len(values) < 5:
        values.extend([q()] * (5 - len(values)))
    for degree in range(len(values) - 1, 4, -1):
        lead = values[degree]
        if not lead:
            continue
        values[degree] = q()
        shift = degree - 5
        for offset, coefficient in ((4, 1), (3, -3), (2, -3), (0, -26)):
            target = shift + offset
            if target >= len(values):
                values.extend([q()] * (target + 1 - len(values)))
            values[target] += lead * coefficient
    return values[:5]


def field_vector(serialized):
    out = [q() for _ in range(5)]
    for degree, (numerator, denominator) in serialized.items():
        out[int(degree)] = q(numerator, denominator)
    return out


def field_add(a, b):
    return [x + y for x, y in zip(a, b)]


def field_mul(a, b):
    product = [q() for _ in range(9)]
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    product[i + j] += x * y
    return reduce_poly(product)


def multiplication_matrix(a):
    matrix = [[q() for _ in range(5)] for _ in range(5)]
    for column in range(5):
        basis = [q() for _ in range(5)]
        basis[column] = q(1)
        product = field_mul(a, basis)
        for row in range(5):
            matrix[row][column] = product[row]
    return matrix


def load_structure():
    system = pickle.loads(SYSTEM.read_bytes())
    pattern = re.compile(r"^C\|(\d+)\|(\d+)\|(\d+)\|(\d+)\|")
    columns = []
    for line in MODULAR_SUPPORT.read_text().splitlines():
        match = pattern.match(line)
        if match:
            columns.append(tuple(map(int, match.groups())))
    if len(columns) != 385:
        raise RuntimeError(f"expected 385 multiplier monomials, got {len(columns)}")

    row_monomials = {(1, 0, 0)}
    for equation, a, b, c in columns:
        for monomial in system[equation - 1]:
            row_monomials.add(add_monomials((a, b, c), monomial))
    row_monomials = sorted(row_monomials, key=lambda m: (sum(m), m))
    if len(row_monomials) != 402:
        raise RuntimeError(f"expected 402 field rows, got {len(row_monomials)}")
    return system, columns, row_monomials


def build_selected_rational_rows(system, columns, row_monomials, pivot_scalar_rows):
    row_index = {m: i for i, m in enumerate(row_monomials)}
    selected = {}
    for selected_index, scalar_row in enumerate(map(int, pivot_scalar_rows)):
        field_row, component = divmod(scalar_row, 5)
        selected.setdefault(field_row, []).append((component, selected_index))

    rows = [dict() for _ in pivot_scalar_rows]
    for block_column, (equation, a, b, c) in enumerate(columns):
        for monomial, serialized_coefficient in system[equation - 1].items():
            field_row = row_index[add_monomials((a, b, c), monomial)]
            wanted = selected.get(field_row)
            if not wanted:
                continue
            block = multiplication_matrix(field_vector(serialized_coefficient))
            for component, selected_index in wanted:
                for basis_column in range(5):
                    value = block[component][basis_column]
                    if value:
                        rows[selected_index][block_column * 5 + basis_column] = value

    rhs = [q() for _ in rows]
    target_field_row = row_index[(1, 0, 0)]
    for field_row, wanted in selected.items():
        if field_row == target_field_row:
            for component, selected_index in wanted:
                if component == 0:
                    rhs[selected_index] = q(1)
    if sum(bool(value) for value in rhs) != 1:
        raise RuntimeError("the selected minor does not contain exactly one target scalar row")
    return rows, rhs


def mod_q(value, prime=P):
    return (int(value.p) % prime) * pow(int(value.q) % prime, -1, prime) % prime


def verify_modular_matrix(rows, rhs, pivot_scalar_rows):
    archived = np.load(MODULAR_MATRIX)
    expected_matrix = archived["A"][pivot_scalar_rows].astype(np.int64) % P
    expected_rhs = archived["b"][pivot_scalar_rows].astype(np.int64) % P
    for i, row in enumerate(rows):
        if mod_q(rhs[i]) != int(expected_rhs[i]):
            raise RuntimeError(f"target mismatch modulo {P} at row {i}")
        expected_nonzero = np.flatnonzero(expected_matrix[i])
        actual_nonzero = {j: mod_q(value) for j, value in row.items() if mod_q(value)}
        if len(expected_nonzero) != len(actual_nonzero):
            raise RuntimeError(f"nonzero-pattern mismatch modulo {P} at row {i}")
        for j, value in actual_nonzero.items():
            if value != int(expected_matrix[i, j]):
                raise RuntimeError(f"matrix mismatch modulo {P} at ({i},{j})")
    print("EXACT_MATRIX_REDUCES_TO_FIXED_P71_PASS", flush=True)


def clear_row_denominators(rows, rhs):
    size = len(rows)
    matrix = fmpz_mat(size, size)
    target = fmpz_mat(size, 1)
    for i, row in enumerate(rows):
        denominator = int(rhs[i].q)
        for value in row.values():
            denominator = math.lcm(denominator, int(value.q))
        for j, value in row.items():
            matrix[i, j] = int(value.p) * (denominator // int(value.q))
        target[i, 0] = int(rhs[i].p) * (denominator // int(rhs[i].q))
        if (i + 1) % 100 == 0:
            print(f"integer rows {i + 1}/{size}", flush=True)
    return matrix, target


def save_solution(solution):
    encoded = [
        (int(solution[i, 0].p), int(solution[i, 0].q))
        for i in range(solution.nrows())
    ]
    SOLUTION_CACHE.write_bytes(pickle.dumps(encoded, protocol=4))


def load_solution():
    encoded = pickle.loads(SOLUTION_CACHE.read_bytes())
    return [q(n, d) for n, d in encoded]


def verify_full_identity(system, columns, row_monomials, solution):
    # Put the whole solution over one denominator.  This turns the final replay
    # into integer dot products and avoids hundreds of thousands of enormous
    # rational gcd computations.
    common_denominator = 1
    for value in solution:
        common_denominator = math.lcm(common_denominator, int(value.q))
    numerators = [
        int(value.p) * (common_denominator // int(value.q)) for value in solution
    ]

    all_scalar_rows = np.arange(len(row_monomials) * 5, dtype=np.int64)
    rows, rhs = build_selected_rational_rows(
        system, columns, row_monomials, all_scalar_rows
    )
    scalar_products = 0
    for i, row in enumerate(rows):
        row_denominator = int(rhs[i].q)
        for value in row.values():
            row_denominator = math.lcm(row_denominator, int(value.q))
        left = 0
        for j, value in row.items():
            left += (
                int(value.p)
                * (row_denominator // int(value.q))
                * numerators[j]
            )
            scalar_products += 1
        right = (
            int(rhs[i].p)
            * (row_denominator // int(rhs[i].q))
            * common_denominator
        )
        if left != right:
            raise RuntimeError(f"full exact scalar residual is nonzero at row {i}")
    print(f"FULL_SCALAR_PRODUCTS={scalar_products}", flush=True)
    print("H_EXACT_IDENTITY_PASS", flush=True)


def format_q(value):
    numerator, denominator = int(value.p), int(value.q)
    return str(numerator) if denominator == 1 else f"{numerator}/{denominator}"


def format_field(value):
    terms = []
    for degree, coefficient in enumerate(value):
        if not coefficient:
            continue
        cq = format_q(coefficient)
        if degree == 0:
            terms.append(f"({cq})")
        elif degree == 1:
            terms.append(f"({cq})*w")
        else:
            terms.append(f"({cq})*w^{degree}")
    return "+".join(terms).replace("+-", "-") if terms else "0"


def write_certificate(columns, solution):
    lines = [
        "JC2_H_CERTIFICATE_V1",
        "minpoly=w^5-w^4+3*w^3+3*w^2+26",
        "identity=h-T1*E1-T2*E2-T3*E3-T4*E4=0",
        f"coefficients={len(columns)}",
    ]
    for block_column, (equation, a, b, c) in enumerate(columns):
        value = solution[block_column * 5 : (block_column + 1) * 5]
        lines.append(f"C|{equation}|{a}|{b}|{c}|{format_field(value)}")
    CERTIFICATE.write_text("\n".join(lines) + "\n")
    print(f"wrote {CERTIFICATE} ({CERTIFICATE.stat().st_size} bytes)", flush=True)


def main():
    ctx.threads = max(1, min(8, (getattr(__import__("os"), "cpu_count")() or 1)))
    system, columns, row_monomials = load_structure()
    pivot_scalar_rows = np.load(PIVOT_ROWS)

    if SOLUTION_CACHE.exists() and "--rebuild" not in sys.argv:
        print(f"loading cached exact solution {SOLUTION_CACHE}", flush=True)
        solution = load_solution()
    else:
        started = time.monotonic()
        rows, rhs = build_selected_rational_rows(
            system, columns, row_monomials, pivot_scalar_rows
        )
        print(f"built {len(rows)} exact selected rows in {time.monotonic()-started:.1f}s", flush=True)
        verify_modular_matrix(rows, rhs, pivot_scalar_rows)
        matrix, target = clear_row_denominators(rows, rhs)
        del rows, rhs
        print(f"starting FLINT exact solve {matrix.nrows()}x{matrix.ncols()}", flush=True)
        solution_matrix = matrix.solve(target)
        save_solution(solution_matrix)
        solution = [solution_matrix[i, 0] for i in range(solution_matrix.nrows())]
        print(f"exact solve finished in {time.monotonic()-started:.1f}s", flush=True)

    verify_full_identity(system, columns, row_monomials, solution)
    write_certificate(columns, solution)


if __name__ == "__main__":
    main()
