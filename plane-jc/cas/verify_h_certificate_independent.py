#!/usr/bin/env python3
"""Independent verifier for h = sum(T_i E_i) over Q[w]/(w^5-w^4+3w^3+3w^2+26).

This checker deliberately imports none of the certificate-generation modules
and uses no CAS.  It consumes the pinned system pickle and the text multiplier
certificate from Zenodo 21479814 and checks the 2010 rational scalar rows with
gmpy2 integer/rational arithmetic.
"""

from __future__ import annotations

import argparse
import hashlib
import pickle
import re
import sys
from pathlib import Path

from gmpy2 import lcm, mpq, mpz


TERM = re.compile(r"^\((-?\d+)(?:/(\d+))?\)(?:\*w(?:\^(\d+))?)?$")
SYSTEM_SHA256 = "5a6e423d74ef09fc9c7a7282c500bda566018d7e56a93124665796bbe417cedf"
CERTIFICATE_SHA256 = "0e48ffab32469ef8405a6945b16cf1521ddeb3c592ae4e5051968110a4dc656a"


def sha256_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def zero():
    return mpq(0)


def reduce_field(coefficients):
    values = list(coefficients) + [zero()] * max(0, 5 - len(coefficients))
    for degree in range(len(values) - 1, 4, -1):
        lead = values[degree]
        if lead:
            values[degree] = zero()
            shift = degree - 5
            # w^5 = w^4 - 3 w^3 - 3 w^2 - 26
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
            raise RuntimeError(f"cannot parse coefficient: {item[:120]}")
        numerator = int(match.group(1))
        denominator = int(match.group(2) or 1)
        degree = int(match.group(3) or (1 if "*w" in item else 0))
        if not 0 <= degree < 5:
            raise RuntimeError("certificate coefficient is not field-reduced")
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("system", type=Path, help="hne0_polred.pkl")
    parser.add_argument("certificate", type=Path, help="h_certificate_exact.txt")
    args = parser.parse_args()
    if hasattr(sys, "set_int_max_str_digits"):
        sys.set_int_max_str_digits(0)

    system_hash = sha256_file(args.system)
    certificate_hash = sha256_file(args.certificate)
    if system_hash != SYSTEM_SHA256:
        raise RuntimeError(f"unexpected system SHA-256: {system_hash}")
    if certificate_hash != CERTIFICATE_SHA256:
        raise RuntimeError(f"unexpected certificate SHA-256: {certificate_hash}")

    system = pickle.loads(args.system.read_bytes())
    certificate_terms = []
    supports = set()
    header = []
    with args.certificate.open() as handle:
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

    expected_header = [
        "JC2_H_CERTIFICATE_V1",
        "minpoly=w^5-w^4+3*w^3+3*w^2+26",
        "identity=h-T1*E1-T2*E2-T3*E3-T4*E4=0",
        "coefficients=385",
    ]
    if header != expected_header or len(supports) != 385:
        raise RuntimeError("unexpected certificate header or support count")

    common_denominator = mpz(1)
    for _, value in certificate_terms:
        for coordinate in value:
            common_denominator = lcm(common_denominator, coordinate.denominator)
    solution_numerators = [
        coordinate.numerator * (common_denominator // coordinate.denominator)
        for _, value in certificate_terms
        for coordinate in value
    ]

    row_monomials = {(1, 0, 0)}
    for (equation, a, b, c), _ in certificate_terms:
        for monomial in system[equation - 1]:
            row_monomials.add((a + monomial[0], b + monomial[1], c + monomial[2]))
    row_monomials = sorted(row_monomials, key=lambda m: (sum(m), m))
    if len(row_monomials) != 402:
        raise RuntimeError(f"expected 402 coefficient rows, got {len(row_monomials)}")
    row_index = {monomial: i for i, monomial in enumerate(row_monomials)}
    scalar_rows = [dict() for _ in range(402 * 5)]
    field_products = 0
    for block_column, ((equation, a, b, c), _) in enumerate(certificate_terms):
        for monomial, serialized_coefficient in system[equation - 1].items():
            target = (a + monomial[0], b + monomial[1], c + monomial[2])
            block = multiplication_matrix(decode_system_coefficient(serialized_coefficient))
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
        right = row_denominator * common_denominator if row_number == target_scalar_row else 0
        if left != right:
            raise RuntimeError(f"nonzero scalar residual at row {row_number}")

    assert field_products == 13_410
    assert scalar_products == 335_250
    print(f"COEFFICIENT_FIELD_PRODUCTS={field_products}")
    print(f"SCALAR_PRODUCTS={scalar_products}")
    print("INDEPENDENT_H_CERTIFICATE_PASS")


if __name__ == "__main__":
    main()
