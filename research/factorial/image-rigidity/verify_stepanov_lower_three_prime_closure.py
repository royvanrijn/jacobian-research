#!/usr/bin/env python3
"""Exact guards for the written uniform-jet proof; finite guards are not its proof."""

import argparse
import json
from math import comb, factorial, prod
from pathlib import Path

import sympy as sp


def inverse_matrix(matrix, modulus, prime):
    size = len(matrix)
    rows = [[value % modulus for value in row]
            + [int(i == j) for j in range(size)]
            for i, row in enumerate(matrix)]
    for column in range(size):
        pivot = next(row for row in range(column, size)
                     if rows[row][column] % prime)
        rows[column], rows[pivot] = rows[pivot], rows[column]
        inverse = pow(rows[column][column], -1, modulus)
        rows[column] = [value * inverse % modulus for value in rows[column]]
        for row in range(size):
            if row != column:
                factor = rows[row][column]
                rows[row] = [(left - factor * right) % modulus
                             for left, right in zip(rows[row], rows[column])]
    assert all(rows[i][j] == int(i == j) for i in range(size) for j in range(size))
    return [row[size:] for row in rows]


def multiply(matrix, vector, modulus):
    return [sum(a * b for a, b in zip(row, vector)) % modulus for row in matrix]


def det_mod(matrix, prime):
    rows = [[value % prime for value in row] for row in matrix]
    result = 1
    for column in range(len(rows)):
        pivot = next((row for row in range(column, len(rows)) if rows[row][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            result = -result
            rows[column], rows[pivot] = rows[pivot], rows[column]
        value = rows[column][column]
        result = result * value % prime
        inverse = pow(value, -1, prime)
        for row in range(column + 1, len(rows)):
            factor = rows[row][column] * inverse % prime
            rows[row] = [(a - factor * b) % prime for a, b in zip(rows[row], rows[column])]
    return result % prime


def actual_kernel(prime, terminal, precision=6, u_degree=3):
    n = terminal - 1
    gamma = 3 * n + 1
    c = gamma + prime * n + prime
    d = c + 1
    assert n >= 2 and (gamma + 1) % prime == 0 and d % (prime * prime) == 0
    modulus = prime ** precision
    size = prime * (u_degree + 1)
    b = [1] + [0] * (size - 1)
    lower_b = None
    # These are coefficients of the actual polynomials in z, BEFORE cyclic
    # projection. The derivative retains the full z exponent k, not k mod p.
    for index in range(n + 1):
        if index == n:
            lower_b = b[:]
        b = [((b[k - 1] if k else 0)
              + (prime * index + gamma + prime - 1 - k) * b[k]) % modulus
             for k in range(size)]
    ls = [b]
    for j in range(prime - 2):
        previous = ls[-1]
        ls.append([((k - j) * previous[k]
                    - (previous[k - 1] if k else 0)) % modulus
                   for k in range(size)])
    matrices = []
    rhs = []
    for degree in range(u_degree + 1):
        matrices.append([[-ls[prime - 1 - m][prime * degree + row] % modulus
                          for m in range(1, prime)]
                         for row in range(1, prime)])
        rhs.append([-lower_b[prime * degree + row] % modulus
                    for row in range(1, prime)])
    minor = det_mod(matrices[0], prime)
    assert minor != 0
    inverse = inverse_matrix(matrices[0], modulus, prime)
    kernel = []
    for degree in range(u_degree + 1):
        vector = rhs[degree][:]
        for i in range(1, degree + 1):
            product = multiply(matrices[i], kernel[degree - i], modulus)
            vector = [(a - b) % modulus for a, b in zip(vector, product)]
        kernel.append(multiply(inverse, vector, modulus))
    for degree in range(u_degree + 1):
        lhs = [0] * (prime - 1)
        for i in range(degree + 1):
            product = multiply(matrices[i], kernel[degree - i], modulus)
            lhs = [(a + b) % modulus for a, b in zip(lhs, product)]
        assert lhs == rhs[degree]
    scalar = [((c + 2) * vector[0] + vector[1]) % modulus for vector in kernel]
    assert scalar[0] % (prime * prime) == prime
    assert scalar[1] % prime == prime - 1
    constant_modulus = prime ** min(n, precision)
    denominator = prod((2 - d + r) for r in range(prime) if r != prime - 2)
    expected_constant = (2 * d - prime) * pow(denominator, -1, constant_modulus) % constant_modulus
    assert scalar[0] % constant_modulus == expected_constant
    return {"p": prime, "N": terminal, "gamma": gamma, "D": d,
            "coefficient_modulus": modulus, "u_degree": u_degree,
            "minor_determinant_mod_p": minor,
            "normalized_kernel_coefficients": kernel,
            "S_coefficients": scalar,
            "constant_formula_modulus": constant_modulus,
            "constant_formula_value": expected_constant}


def modulo_p_and_u_squared(poly, x, u, prime):
    return {(ix, iu): int(coefficient) % prime
            for (ix, iu), coefficient in sp.Poly(sp.expand(poly), x, u).terms()
            if iu < 2 and int(coefficient) % prime}


def run():
    x, u, p, c = sp.symbols("x u p c")
    rising_next = (p - 2) * (p - 3) / 2
    h_next = p * (p - 1) / 2 - p + 2
    assert sp.expand(2 * c + 1 + rising_next - h_next - (2 * (c + 1) - p)) == 0
    structural = []
    for prime in [5, 7, 11, 13]:
        rho = 2 - u
        h = sum(x ** (prime - 1 - i) * rho ** i for i in range(prime)) - 1
        assert not modulo_p_and_u_squared((x - rho) * h - (x ** prime - x - u), x, u, prime)
        assert not modulo_p_and_u_squared(h.subs(x, 2) + 1, x, u, prime)
        top = -1
        next_coefficient = -rho
        s = (2 * (-1) + 1 + (prime - 2) * (prime - 3) // 2) * top - next_coefficient
        assert not modulo_p_and_u_squared(s + u, x, u, prime)
        structural.append({"p": prime, "rho_mod_u2": "2-u",
                           "H_u_at_two_mod_p_u2": -1, "S_mod_p_u2": "-u"})

    samples = [(5, 22), (5, 47), (5, 172), (5, 422), (5, 547),
               (7, 5), (7, 103), (7, 446), (7, 789), (7, 1818), (7, 2161),
               (11, pow(14, -1, 121)), (13, pow(16, -1, 169))]
    actual = [actual_kernel(prime, terminal) for prime, terminal in samples]

    residue_records = []
    for prime, n0 in [(5, 1), (7, 4)]:
        q = sum(comb(n0, k) * prod(range(n0 + 2 - k, n0 + 2)) * x ** (n0 - k)
                for k in range(n0 + 1))
        q_poly = sp.Poly(q, x, modulus=prime)
        assert sp.gcd(q_poly, q_poly.diff()).degree() == 0
        assert q_poly.eval(1) % prime != 0
        if prime == 5:
            assert q_poly == sp.Poly(x + 2, x, modulus=5)
        else:
            assert q_poly == sp.Poly((x - 3) * (x + 3) * (x * x - x + 3), x, modulus=7)
            assert all((r * r - r + 3) % 7 for r in range(7))
        compatibility = []
        for residue in range(1, prime):
            common = sp.gcd(q_poly, sp.Poly(x * x - x - 3 * residue, x, modulus=prime)).monic()
            compatibility.append({"D_over_p3_mod_p": residue,
                                  "gcd_coefficients_ascending_mod_p":
                                  [int(common.nth(k)) % prime for k in range(common.degree() + 1)]})
        residue_records.append({"p": prime, "n_mod_p": n0,
                                "lower_unit_polynomial": str(q_poly.as_expr()),
                                "compatibility": compatibility})
    assert [r["D_over_p3_mod_p"] for r in residue_records[0]["compatibility"]
            if len(r["gcd_coefficients_ascending_mod_p"]) > 1] == [2]
    assert [r["D_over_p3_mod_p"] for r in residue_records[1]["compatibility"]
            if len(r["gcd_coefficients_ascending_mod_p"]) > 1] == [2, 4, 6]
    residual_terminal = (1 + 6 * 343) * pow(10, -1, 2401) % 2401
    assert residual_terminal == 446

    # Exact Gram identity and its necessary coefficient-field distinction.
    gram_degree = 10
    for i in range(gram_degree + 1):
        for j in range(gram_degree + 1):
            assert factorial(i + j) == factorial(i) * factorial(j) * sum(
                comb(i, k) * comb(j, k) for k in range(min(i, j) + 1))
    complex_square = sp.Poly(sp.expand((x - 1 + sp.I) ** 2), x)
    complex_moment = sp.simplify(sum(complex_square.nth(j) * factorial(j)
                                   for j in range(complex_square.degree() + 1)))
    assert complex_moment == 0

    return {"status": "PASS", "sympy_version": sp.__version__,
            "scope": "exact guards for the uniform normalized-jet proof used in both closures; the separate Hensel verifier certifies the final gap-seven branch; finite examples do not prove all-order identities",
            "structural_identities": structural, "actual_series_certificates": actual,
            "residue_root_checks": residue_records,
            "gap7_residual_terminal_mod2401": residual_terminal,
            "gram_identity_degree_guard": gram_degree,
            "complex_positivity_countercontrol": {"g": "x-1+i", "Lambda_g_squared": 0}}


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run without -O or -OO: exact assertions must execute.")
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    serialized = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(serialized)
    print(serialized, end="")
