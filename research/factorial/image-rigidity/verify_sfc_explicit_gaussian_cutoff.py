#!/usr/bin/env python3
"""Explicit radius-bound arithmetic and integer normalization audit.

The cutoff uses the imported Basu--Roy Theorem 3 and the written proof of
compactified corner exclusion. This script verifies the integer arithmetic
and, with --audit, the full normalization for the scale-one four-column case.
It does not independently verify either imported theorem or the limit proof.
"""

import argparse
import json
from math import comb, factorial


def parameters(size, degree, maximum_order):
    assert size >= 2 and degree >= size - 1 and maximum_order >= size
    cofactor = factorial(size - 1) * degree ** ((size - 1) ** 2)
    coefficient_bound = ((2 * size ** 2 * cofactor) ** maximum_order
                         * (maximum_order * (1 + 2 * degree))
                         ** (maximum_order * degree))
    tau = coefficient_bound.bit_length() + 1
    variables = 2 * size
    polynomials = 3 * size + 2
    total_degree = maximum_order * (3 * degree + 1)
    count = (2 * total_degree + 1) * (2 * total_degree) ** (variables - 1)
    delta = variables * (2 * total_degree - 1) + 2
    # The extra bit(d+1) safely pads both the theorem statement and the
    # final estimate in the author's Section 5 proof.
    height = (2 * tau + count.bit_length()
              + (variables + 1) * (total_degree + 1).bit_length()
              + polynomials.bit_length() + 3)
    exponent = (4 * count * delta * height
                + 2 * (count + 1).bit_length() + (2 * variables).bit_length())
    return dict(size=size, degree=degree, maximum_order=maximum_order,
                variables=variables, polynomials=polynomials,
                polynomial_degree=total_degree, coefficient_bits=tau,
                N=count, Delta=delta, A=height, cutoff_exponent=exponent,
                cutoff=f"2^{exponent}")


def audit_four_column():
    import sympy as sp
    from flint import fmpz_mpoly_ctx
    from verify_sfc_gaussian_diagonal_nonlifting import count_vectors

    ring = fmpz_mpoly_ctx.get(("a0", "a1", "a2", "a3", "sigma", "r"))
    *coefficients, sigma, r = ring.gens()
    offsets = (0, 3, 4, 9)
    degree = 9
    matrix = sp.Matrix([[comb(x, i) for x in offsets] for i in range(4)])
    determinant = int(matrix.det())
    assert determinant
    adjugate = matrix.adjugate()
    basis = [sum(int(adjugate[j, i]) * coefficients[i] * sigma ** (3 - i)
                 for i in range(4)) for j in range(4)]
    rows = []
    for order in (1, 2, 3, 5):
        products = {}
        for length in {sum(x * a for x, a in zip(offsets, alpha))
                       for alpha in count_vectors(order)}:
            product = ring.constant(1)
            for j in range(1, length + 1):
                product *= order + j * sigma ** 2 * r
            for j in range(1, order * degree - length + 1):
                product *= order + j * sigma ** 2 * (1 - r)
            products[length] = product
        numerator = ring.constant(0)
        for alpha in count_vectors(order):
            weight = factorial(order)
            term = ring.constant(1)
            for p, a in zip(basis, alpha):
                weight //= factorial(a)
                term *= p ** a
            length = sum(x * a for x, a in zip(offsets, alpha))
            numerator += weight * term * products[length]
        terms = numerator.to_dict()
        assert min(e[4] for e in terms) >= 3 * order
        normalized = {e[:4] + (e[4] - 3 * order, e[5]): c
                      for e, c in terms.items()}
        # Replacing sigma by 1/w and multiplying by w^(2mD).
        actual_degree = max(sum(e[:4]) + 2 * order * degree - e[4] + e[5]
                            for e in normalized)
        assert actual_degree <= order * (3 * degree + 1)
        assert all(e[4] <= 2 * order * degree for e in normalized)
        bound = ((32 * 6 * degree ** 9) ** order
                 * (order * (1 + 2 * degree)) ** (order * degree))
        # The factor 2^m covers real/imaginary expansion of m coefficients.
        actual_norm_bound = 2 ** order * sum(abs(int(c)) for c in normalized.values())
        assert actual_norm_bound <= bound
        # At sigma=0, recover the Gaussian moment independently by count
        # vectors in the powers 1,z,z^2,z^3.
        gaussian = {}
        for alpha in count_vectors(order):
            power = sum(i * a for i, a in enumerate(alpha))
            if power % 2:
                continue
            weight = factorial(order)
            for a in alpha:
                weight //= factorial(a)
            double_factorial = 1
            for j in range(1, power, 2):
                double_factorial *= j
            value = weight * determinant ** order * double_factorial
            value *= order ** (order * degree - power // 2)
            gaussian[alpha] = value
        leading = {e[:4]: int(c) for e, c in normalized.items() if e[4] == 0}
        assert all(e[5] == 0 for e in normalized if e[4] == 0)
        assert leading == gaussian
        rows.append(dict(order=order, terms=len(normalized),
                         actual_w_degree=int(actual_degree),
                         coefficient_norm_bits=actual_norm_bound.bit_length()))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scale", type=int, default=1)
    parser.add_argument("--audit", action="store_true")
    args = parser.parse_args()
    if args.scale < 1:
        parser.error("The scale must be positive.")
    result = parameters(4, 9 * args.scale, 5)
    if args.scale == 1:
        assert result["cutoff_exponent"] == 373029362779174993920000137
    result.update(status="PASS", scale=args.scale,
                  scope="explicit sufficient corner height cutoff; extremely conservative; not a computed exceptional set",
                  assurance="internal arithmetic/normalization replay; imported radius theorem and compactified exclusion are written proof inputs")
    if args.audit:
        result["scale_one_full_normalization_audit"] = audit_four_column()
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
