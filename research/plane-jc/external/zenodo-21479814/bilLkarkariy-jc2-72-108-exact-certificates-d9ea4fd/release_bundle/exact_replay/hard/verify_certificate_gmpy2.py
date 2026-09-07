#!/usr/bin/env python3
"""Independent gmpy2 replay of h = T1*E1 + ... + T4*E4."""

from __future__ import annotations

import pickle
import re
import sys
from pathlib import Path

from gmpy2 import lcm, mpq, mpz


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "hard" / "h_certificate_exact.txt"
SYSTEM = ROOT / "hne0_polred.pkl"
TERM = re.compile(r"^\((-?\d+)(?:/(\d+))?\)(?:\*w(?:\^(\d+))?)?$")

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


def zero():
    return mpq(0)


def reduce_field(coefficients):
    values = list(coefficients)
    if len(values) < 5:
        values.extend(zero() for _ in range(5 - len(values)))
    for degree in range(len(values) - 1, 4, -1):
        lead = values[degree]
        if not lead:
            continue
        values[degree] = zero()
        shift = degree - 5
        # w^5 = w^4 - 3*w^3 - 3*w^2 - 26.
        values[shift + 4] += lead
        values[shift + 3] -= 3 * lead
        values[shift + 2] -= 3 * lead
        values[shift] -= 26 * lead
    return tuple(values[:5])


def field_mul(left, right):
    product = [zero() for _ in range(9)]
    for i, a in enumerate(left):
        if a:
            for j, b in enumerate(right):
                if b:
                    product[i + j] += a * b
    return reduce_field(product)


def field_add(left, right):
    return tuple(a + b for a, b in zip(left, right))


def decode_system_coefficient(serialized):
    out = [zero() for _ in range(5)]
    for degree, (numerator, denominator) in serialized.items():
        out[int(degree)] = mpq(int(numerator), int(denominator))
    return tuple(out)


def parse_field(expression):
    if expression == "0":
        return tuple(zero() for _ in range(5))
    out = [zero() for _ in range(5)]
    for item in expression.split("+"):
        match = TERM.match(item)
        if not match:
            raise RuntimeError(f"cannot parse certificate coefficient: {item[:120]}")
        numerator = int(match.group(1))
        denominator = int(match.group(2) or 1)
        degree = int(match.group(3) or (1 if "*w" in item else 0))
        if not 0 <= degree < 5:
            raise RuntimeError("certificate coefficient is not reduced modulo the minpoly")
        out[degree] += mpq(numerator, denominator)
    return tuple(out)


def multiplication_matrix(value):
    matrix = [[zero() for _ in range(5)] for _ in range(5)]
    for column in range(5):
        basis = [zero() for _ in range(5)]
        basis[column] = mpq(1)
        product = field_mul(value, basis)
        for row in range(5):
            matrix[row][column] = product[row]
    return matrix


def main():
    system = pickle.loads(SYSTEM.read_bytes())
    certificate_terms = []
    supports = set()
    header = []
    with CERTIFICATE.open() as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\n")
            if line_number <= 4:
                header.append(line)
                continue
            parts = line.split("|", 5)
            if len(parts) != 6 or parts[0] != "C":
                raise RuntimeError(f"malformed certificate line {line_number}")
            equation, a, b, c = map(int, parts[1:5])
            support = (equation, a, b, c)
            if support in supports:
                raise RuntimeError(f"duplicate multiplier support {support}")
            supports.add(support)
            certificate_terms.append((support, parse_field(parts[5])))

    if header != [
        "JC2_H_CERTIFICATE_V1",
        "minpoly=w^5-w^4+3*w^3+3*w^2+26",
        "identity=h-T1*E1-T2*E2-T3*E3-T4*E4=0",
        "coefficients=385",
    ]:
        raise RuntimeError("unexpected exact-certificate header")
    if len(supports) != 385:
        raise RuntimeError(f"expected 385 coefficients, got {len(supports)}")

    # Normalize the 1925 huge rational coordinates to one denominator.
    common_denominator = mpz(1)
    for _, value in certificate_terms:
        for coordinate in value:
            common_denominator = lcm(common_denominator, coordinate.denominator)
    solution_numerators = []
    for _, value in certificate_terms:
        for coordinate in value:
            solution_numerators.append(
                coordinate.numerator * (common_denominator // coordinate.denominator)
            )

    row_monomials = {(1, 0, 0)}
    for (equation, a, b, c), _ in certificate_terms:
        for monomial in system[equation - 1]:
            row_monomials.add(
                (a + monomial[0], b + monomial[1], c + monomial[2])
            )
    row_monomials = sorted(row_monomials, key=lambda m: (sum(m), m))
    if len(row_monomials) != 402:
        raise RuntimeError(f"expected 402 coefficient rows, got {len(row_monomials)}")
    row_index = {m: i for i, m in enumerate(row_monomials)}
    scalar_rows = [dict() for _ in range(len(row_monomials) * 5)]
    field_products = 0
    for block_column, ((equation, a, b, c), _) in enumerate(certificate_terms):
        for monomial, serialized_coefficient in system[equation - 1].items():
            target = (a + monomial[0], b + monomial[1], c + monomial[2])
            block = multiplication_matrix(
                decode_system_coefficient(serialized_coefficient)
            )
            base_row = row_index[target] * 5
            base_column = block_column * 5
            for i in range(5):
                for j in range(5):
                    if block[i][j]:
                        scalar_rows[base_row + i][base_column + j] = block[i][j]
            field_products += 1

    target_scalar_row = row_index[(1, 0, 0)] * 5
    scalar_products = 0
    for row_number, row in enumerate(scalar_rows):
        row_denominator = mpz(1)
        for value in row.values():
            row_denominator = lcm(row_denominator, value.denominator)
        left = mpz(0)
        for column, value in row.items():
            left += (
                value.numerator
                * (row_denominator // value.denominator)
                * solution_numerators[column]
            )
            scalar_products += 1
        right = (
            row_denominator * common_denominator
            if row_number == target_scalar_row
            else mpz(0)
        )
        if left != right:
            raise RuntimeError(f"independent scalar residual is nonzero at row {row_number}")
    print(f"COEFFICIENT_FIELD_PRODUCTS={field_products}")
    print(f"SCALAR_PRODUCTS={scalar_products}")
    print("GMPY2_EXACT_PASS")


if __name__ == "__main__":
    main()
