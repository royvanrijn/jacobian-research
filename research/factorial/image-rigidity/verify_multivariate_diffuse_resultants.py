#!/usr/bin/env python3
"""Certify the Curie–Weiss consecutive pairs at finitely many orders.

The computations use exact integer arithmetic in python-flint. The
default range is 2 <= m <= 12 for the pair of orders m and m+1.
"""

from argparse import ArgumentParser
from math import comb

from flint import __version__, fmpz, fmpz_mpoly_ctx, fmpz_poly


def angle_polynomial(order, other_order, q, x):
    """Remove the odd factor z+1 and use x=z+1/z."""
    coefficients = [
        comb(order, k) * q ** (other_order * k * (order - k))
        for k in range(order + 1)
    ]
    if order % 2:
        quotient = [coefficients[0]]
        for coefficient in coefficients[1:-1]:
            quotient.append(coefficient - quotient[-1])
        rebuilt = [
            quotient[0],
            *(left + right for left, right in zip(quotient, quotient[1:])),
            quotient[-1],
        ]
        if rebuilt != coefficients:
            raise ArithmeticError("The odd member is not divisible by z+1")
        coefficients = quotient
    if coefficients != coefficients[::-1]:
        raise ArithmeticError("The coefficient list is not reciprocal")

    middle = len(coefficients) // 2
    result = coefficients[middle]
    previous, power_sum = 2, x
    for j in range(1, middle + 1):
        result += coefficients[middle - j] * power_sum
        previous, power_sum = power_sum, x * power_sum - previous
    return result


def univariate(polynomial):
    """Read a polynomial in q,x which is already independent of x."""
    terms = polynomial.to_dict()
    if any(x_power for _, x_power in terms):
        raise ArithmeticError("Unexpected dependence on x")
    degree = max((q_power for q_power, _ in terms), default=0)
    return fmpz_poly([terms.get((i, 0), 0) for i in range(degree + 1)])


def certify_no_interior_zeros(polynomial):
    """Return endpoint multiplicities, reduced degree and coefficient sign."""
    if not polynomial:
        raise ArithmeticError("The polynomial is identically zero")
    q = fmpz_poly([0, 1])
    at_zero = at_one = 0
    while polynomial[0] == 0:
        polynomial, remainder = divmod(polynomial, q)
        if remainder:
            raise ArithmeticError("Nonexact division by q")
        at_zero += 1
    while polynomial(1) == 0:
        polynomial, remainder = divmod(polynomial, q - 1)
        if remainder:
            raise ArithmeticError("Nonexact division by q-1")
        at_one += 1

    # The coefficient of v^k after substituting q=v/(1+v) and
    # multiplying by (1+v)^degree is a sum of a_i*binom(degree-i,k-i).
    # A convolution computes those coefficients multiplied by the
    # positive factors degree!*(degree-k)!.
    degree = polynomial.degree()
    factorials = [fmpz(1)]
    for k in range(1, degree + 1):
        factorials.append(k * factorials[-1])
    first = fmpz_poly([
        polynomial[i] * factorials[degree - i] for i in range(degree + 1)
    ])
    second = fmpz_poly([
        factorials[degree] // factorials[i] for i in range(degree + 1)
    ])
    product = first * second
    sign = 1 if polynomial[0] > 0 else -1
    if any(sign * product[k] <= 0 for k in range(degree + 1)):
        raise ArithmeticError("The transformed coefficients differ in sign")
    return at_zero, at_one, degree, sign


def check_controls():
    q = fmpz_poly([0, 1])
    positive = q**2 * (q - 1)**3 * (1 + q)**2
    if certify_no_interior_zeros(positive) != (2, 3, 2, 1):
        raise ArithmeticError("Endpoint or coefficient control failed")
    for bad in (fmpz_poly(), 2 * q - 1):
        try:
            certify_no_interior_zeros(bad)
        except ArithmeticError:
            pass
        else:
            raise ArithmeticError("A negative control passed")


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--max-order", type=int, default=12,
                        help="largest lower moment order (default 12)")
    args = parser.parse_args()
    if args.max_order < 2:
        parser.error("--max-order must be at least 2")
    check_controls()
    q, x = fmpz_mpoly_ctx.get(("q", "x"), "lex").gens()
    print(f"python-flint {__version__}", flush=True)
    for m in range(2, args.max_order + 1):
        lower = angle_polynomial(m, m + 1, q, x)
        upper = angle_polynomial(m + 1, m, q, x)
        resultant = univariate(lower.resultant(upper, "x"))
        at_zero, at_one, degree, sign = certify_no_interior_zeros(resultant)

        # z=-1 is the one root removed from the odd member. Its
        # occurrence in the even member is equivalent to x=-2 there.
        even = lower if m % 2 == 0 else upper
        even_at_minus_one = univariate(
            even.compose(q, q.context().constant(-2))
        )
        _, _, even_degree, _ = certify_no_interior_zeros(even_at_minus_one)
        print(
            f"m={m}: resultant q^{at_zero}(q-1)^{at_one}, "
            f"reduced degree {degree}, {degree + 1} strict "
            f"{'positive' if sign == 1 else 'negative'} signs; "
            f"even z=-1 reduced degree {even_degree}",
            flush=True,
        )
    print(f"Verified all {args.max_order - 1} pairs and three controls.")


if __name__ == "__main__":
    main()
