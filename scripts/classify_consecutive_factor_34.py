#!/usr/bin/env python3
"""Exact rank/contact and point-count checks for the (3,4) factor row.

The script has two independent jobs.

* It verifies the expanded motivic formula from
  ``extended-geometry/CONSECUTIVE_FACTOR_CONTACT_CLASSIFICATION.md`` after
  replacing every class by an F_5-point count, for representatives of all
  four middle-catalecticant ranks (and both rank-two boundary types).
* It exhausts every hyperplane over F_2 and F_3.  Those small fields are only
  an early filter: equality with q^7 occurs, so point counts alone do not
  prove the characteristic-zero theorem.

Coefficient vectors are in descending binary-monomial order.  Thus a
functional ``ell=(e_0,...,e_n)`` evaluates the coefficient vector of a binary
form by the ordinary dot product, and its middle catalecticant is Hankel.
"""

from collections import Counter, defaultdict
from itertools import product

import sympy as sp


z = sp.Symbol("z")


def projective_points(degree, q):
    """Canonical representatives for P^degree(F_q)."""

    for first_nonzero in range(degree + 1):
        for tail in product(range(q), repeat=degree - first_nonzero):
            yield (0,) * first_nonzero + (1,) + tail


def multiply(left, right, q):
    answer = [0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            answer[i + j] = (answer[i + j] + x * y) % q
    return tuple(answer)


def trim(polynomial, q):
    polynomial = list(polynomial)
    while polynomial and polynomial[0] % q == 0:
        polynomial.pop(0)
    return polynomial


def remainder(dividend, divisor, q):
    dividend, divisor = trim(dividend, q), trim(divisor, q)
    while dividend and len(dividend) >= len(divisor):
        scale = dividend[0] * pow(divisor[0], -1, q) % q
        for i, value in enumerate(divisor):
            dividend[i] = (dividend[i] - scale * value) % q
        dividend = trim(dividend, q)
    return dividend


def coprime(left, right, q):
    # Simultaneous leading-coefficient zero is the common root at infinity;
    # the Euclidean algorithm checks every finite root.
    if left[0] % q == right[0] % q == 0:
        return False
    left, right = trim(left, q), trim(right, q)
    while right:
        left, right = right, remainder(left, right, q)
    return len(left) == 1


def matrix_rank(matrix, q):
    matrix = [list(row) for row in matrix]
    rank = 0
    for column in range(len(matrix[0])):
        pivot = next(
            (row for row in range(rank, len(matrix)) if matrix[row][column] % q),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], -1, q)
        matrix[rank] = [(entry * inverse) % q for entry in matrix[rank]]
        for row in range(len(matrix)):
            if row != rank and matrix[row][column] % q:
                scale = matrix[row][column]
                matrix[row] = [
                    (x - scale * y) % q
                    for x, y in zip(matrix[row], matrix[rank])
                ]
        rank += 1
        if rank == len(matrix):
            break
    return rank


def catalecticant_rank(functional, p, q):
    return matrix_rank(
        [[functional[i + j] for j in range(p + 2)] for i in range(p + 1)],
        q,
    )


def characteristic_zero_contact_partition(functional):
    """Multiplicity partition of ell((X+zY)^7), including infinity."""

    polynomial = sp.Poly(
        sum(sp.binomial(7, i) * functional[i] * z**i for i in range(8)), z
    )
    partition = [7 - polynomial.degree()] if polynomial.degree() < 7 else []
    for factor, multiplicity in sp.factor_list(polynomial)[1]:
        partition.extend([multiplicity] * factor.degree())
    return tuple(sorted(partition, reverse=True))


def square_contraction(functional, divisor, q):
    """Return F |-> ell(divisor^2 F) as a coefficient vector."""

    square = multiply(divisor, divisor, q)
    output_length = len(functional) - len(square) + 1
    return tuple(
        sum(square[i] * functional[i + j] for i in range(len(square))) % q
        for j in range(output_length)
    )


def cubic_discriminant(functional, q):
    # ell((xX+yY)^3) has ordinary coefficients (e0,3e1,3e2,e3).
    a, b, c, d = (
        functional[0],
        3 * functional[1] % q,
        3 * functional[2] % q,
        functional[3],
    )
    return (
        b * b * c * c
        - 4 * a * c**3
        - 4 * b**3 * d
        - 27 * a * a * d * d
        + 18 * a * b * c * d
    ) % q


def cubic_contact_flags(functional, q):
    """Return membership in T (singular/zero), O_3, and O_0."""

    is_zero = all(entry % q == 0 for entry in functional)
    is_singular = is_zero or cubic_discriminant(functional, q) == 0
    is_nonzero_cube = (
        not is_zero and catalecticant_rank(functional, p=1, q=q) == 1
    )
    return is_singular, is_nonzero_cube, is_zero


def kernel_point_count(functional, p, q, projective):
    matrix = [
        [functional[i + j] for j in range(p + 2)] for i in range(p + 1)
    ]
    return sum(
        1
        for vector in projective[p]
        if all(
            sum(vector[i] * matrix[i][j] for i in range(p + 1)) % q == 0
            for j in range(p + 2)
        )
    )


def quintic_contact_counts(functional, q, projective):
    """The five p=2 loci K,Q,T,O_3,O_0 for one quintic functional."""

    kernel = kernel_point_count(functional, p=2, q=q, projective=projective)
    quadratic_zero = sum(
        all(entry == 0 for entry in square_contraction(functional, divisor, q))
        for divisor in projective[2]
    )
    singular = cube = zero = 0
    for divisor in projective[1]:
        flags = cubic_contact_flags(square_contraction(functional, divisor, q), q)
        singular += flags[0]
        cube += flags[1]
        zero += flags[2]
    return kernel, quadratic_zero, singular, cube, zero


def septic_relative_counts(functional, q, projective):
    """Point counts of every relative locus in formula (24) of the note."""

    kernel = kernel_point_count(functional, p=3, q=q, projective=projective)

    rel_kernel = rel_quadratic_zero = rel_singular = 0
    rel_cube = rel_zero = 0
    for linear_divisor in projective[1]:
        quintic = square_contraction(functional, linear_divisor, q)
        data = quintic_contact_counts(quintic, q, projective)
        rel_kernel += data[0]
        rel_quadratic_zero += data[1]
        rel_singular += data[2]
        rel_cube += data[3]
        rel_zero += data[4]

    quadratic_singular = quadratic_cube = quadratic_zero = 0
    for quadratic_divisor in projective[2]:
        flags = cubic_contact_flags(
            square_contraction(functional, quadratic_divisor, q), q
        )
        quadratic_singular += flags[0]
        quadratic_cube += flags[1]
        quadratic_zero += flags[2]

    cubic_zero = sum(
        all(entry == 0 for entry in square_contraction(functional, divisor, q))
        for divisor in projective[3]
    )
    return {
        "K": kernel,
        "K1": rel_kernel,
        "Q1": rel_quadratic_zero,
        "T1": rel_singular,
        "O31": rel_cube,
        "O01": rel_zero,
        "T2": quadratic_singular,
        "O32": quadratic_cube,
        "O02": quadratic_zero,
        "Q3": cubic_zero,
    }


def expanded_formula_count(functional, q, projective):
    loci = septic_relative_counts(functional, q, projective)
    answer = q**7 - q**5 - q**4 * loci["K"]
    answer += q * (loci["T1"] - loci["Q1"] - loci["T2"] + loci["Q3"])
    answer += q**2 * (loci["O32"] - loci["O31"])
    answer += q**3 * (loci["K1"] - loci["O01"] + loci["O02"])
    return answer, loci


def coprime_product_weights(q, projective):
    weights = Counter()
    for cubic in projective[3]:
        for quartic in projective[4]:
            if coprime(cubic, quartic, q):
                weights[multiply(cubic, quartic, q)] += 1
    assert sum(weights.values()) == q**6 * (q + 1)
    return weights


def complement_count(functional, q, product_weights):
    return sum(
        multiplicity
        for binary_product, multiplicity in product_weights.items()
        if sum(x * y for x, y in zip(functional, binary_product)) % q
    )


def verify_expanded_formula():
    q = 5  # The cubic contact classifier requires characteristic != 2,3.
    projective = {degree: tuple(projective_points(degree, q)) for degree in range(5)}
    product_weights = coprime_product_weights(q, projective)
    representatives = {
        "rank 1 / (7)": (1, 0, 0, 0, 0, 0, 0, 0),
        "rank 2 secant": (1, 0, 0, 0, 0, 0, 0, 1),
        "rank 2 tangent / (6,1)": (0, 1, 0, 0, 0, 0, 0, 0),
        "rank 3 / (5,2)": (0, 0, 1, 0, 0, 0, 0, 0),
        "rank 4": (1, 2, 3, 4, 1, 3, 2, 4),
    }
    rows = []
    for label, functional in representatives.items():
        direct = complement_count(functional, q, product_weights)
        formula, loci = expanded_formula_count(functional, q, projective)
        assert direct == formula
        rows.append(
            (
                label,
                catalecticant_rank(functional, 3, q),
                characteristic_zero_contact_partition(functional),
                direct,
                loci["K1"],
            )
        )
    return rows


def exhaustive_small_field_sweep(q):
    projective = {degree: tuple(projective_points(degree, q)) for degree in range(8)}
    product_weights = coprime_product_weights(q, projective)
    distribution = Counter()
    by_rank = defaultdict(Counter)
    affine_count = 0
    affine_by_rank = Counter()
    for functional in projective[7]:
        count = complement_count(functional, q, product_weights)
        rank = catalecticant_rank(functional, p=3, q=q)
        distribution[count] += 1
        by_rank[rank][count] += 1
        if count == q**7:
            affine_count += 1
            affine_by_rank[rank] += 1
    return {
        "hyperplanes": len(projective[7]),
        "distinct_counts": len(distribution),
        "minimum": min(distribution),
        "maximum": max(distribution),
        "affine_count": affine_count,
        "affine_by_rank": dict(sorted(affine_by_rank.items())),
        "rank_sizes": {
            rank: sum(counts.values()) for rank, counts in sorted(by_rank.items())
        },
    }


def main():
    rows = verify_expanded_formula()
    print("PASS (3,4): expanded motivic formula agrees with direct F_5 counts")
    for label, rank, contact, count, relative_kernel in rows:
        print(
            f"  {label:27s} rank={rank}, contact={contact}, "
            f"#U={count}, #K1={relative_kernel}"
        )

    expected = {
        2: {
            "hyperplanes": 255,
            "distinct_counts": 14,
            "minimum": 64,
            "maximum": 108,
            "affine_count": 0,
            "affine_by_rank": {},
            "rank_sizes": {1: 3, 2: 12, 3: 48, 4: 192},
        },
        3: {
            "hyperplanes": 3280,
            "distinct_counts": 29,
            "minimum": 1458,
            "maximum": 1980,
            "affine_count": 0,
            "affine_by_rank": {},
            "rank_sizes": {1: 4, 2: 36, 3: 324, 4: 2916},
        },
    }
    for q in (2, 3):
        summary = exhaustive_small_field_sweep(q)
        assert summary == expected[q], (q, summary)
        print(f"PASS exhaustive F_{q} sweep: {summary}")

    print("PASS: no F_2 or F_3 hyperplane has the affine-space count q^7")


if __name__ == "__main__":
    main()
