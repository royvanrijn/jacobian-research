#!/usr/bin/env python3
"""Exact arithmetic obstructions to copying the physical row-merger collision."""

from flint import (
    fmpq,
    fmpq_poly,
    fmpz,
    fmpz_mod_poly_ctx,
    fmpz_mpoly_ctx,
    fmpz_poly,
    fq_default_ctx,
    fq_default_poly_ctx,
    nmod_mat,
    nmod_poly,
)

from verify_physical_row_merge_interior_collision import (
    bernstein_coefficients,
    moments_with_field,
    multiply,
    raw_moment,
    require,
    sign,
    symmetric_moment,
    weighted_sum,
)


PRIME = 389
DEGREE = 102


def correlation_eliminant(cubic, quartic, variance):
    context = fmpz_mpoly_ctx.get(("q", "t"))

    def as_bivariate(coefficients):
        return context.from_dict(
            {
                (q_degree, t_degree): value
                for t_degree, coefficient in enumerate(coefficients)
                for q_degree, value in enumerate(coefficient)
                if value
            }
        )

    resultant = as_bivariate(cubic).resultant(as_bivariate(quartic), 1)
    require(not resultant.is_zero(), "Nonzero resultant")
    require(resultant.degrees() == (130, 0), "Resultant eliminates t")
    eliminant = fmpz_poly([resultant[j, 0] for j in range(131)])
    common = eliminant.gcd(variance)
    while common.degree() > 0:
        eliminant, remainder = divmod(eliminant, common)
        require(remainder.is_zero(), "Exact saturation by the difference variance")
        common = eliminant.gcd(variance)
    return eliminant // eliminant.content()


def certify_irreducible(polynomial):
    require(fmpz(PRIME).is_prime(), "Prime for the reduction")
    reduced = nmod_poly([int(c % PRIME) for c in polynomial], PRIME)
    require(polynomial.degree() == reduced.degree() == DEGREE, "Preserved eliminant degree")
    generator = nmod_poly([0, 1], PRIME)
    power = generator
    # The proper maximal divisors of 102 are 6, 34 and 51.
    for exponent in range(1, DEGREE + 1):
        power = power.pow_mod(PRIME, reduced)
        if exponent in (6, 34, 51):
            require((power - generator).gcd(reduced).degree() == 0, "Frobenius proper divisor check")
    require(power == generator, "Frobenius check at degree 102")
    return reduced


def certify_coefficient_direction(reduced, cubic, quartic, variance, second):
    modulus = fmpz_mod_poly_ctx(PRIME)([int(c) for c in reduced])
    field = fq_default_ctx(modulus=modulus, var="q")
    field_polynomials = fq_default_poly_ctx(field)

    def reduce_coefficient(coefficient):
        return field([int(c % PRIME) for c in coefficient])

    def reduce_polynomial(coefficients):
        return field_polynomials([reduce_coefficient(c) for c in coefficients])

    cubic_reduced = reduce_polynomial(cubic)
    quartic_reduced = reduce_polynomial(quartic)
    require(cubic_reduced.degree() == 3 and quartic_reduced.degree() == 4, "Preserved t-degrees")
    quadratic = quartic_reduced % cubic_reduced
    require(quadratic.degree() == 2, "Quadratic Euclidean remainder with unit leading coefficient")
    linear = cubic_reduced % quadratic
    require(linear.degree() == 1, "Linear Euclidean remainder with unit leading coefficient")
    field_parameter = -linear[0] / linear[1]
    require(
        cubic_reduced(field_parameter).is_zero() and quartic_reduced(field_parameter).is_zero(),
        "Common field parameter modulo 389",
    )
    variance_reduced = reduce_coefficient(variance)
    require(not variance_reduced.is_zero(), "Variance remains a unit")
    squared_phase = reduce_polynomial(second)(field_parameter) / variance_reduced

    # The phase coordinate is z = i*s, and z**2 = -squared_phase.
    require(
        (-squared_phase) ** ((PRIME**DEGREE - 1) // 2) == field(-1),
        "Nonsquare phase residue in the degree-102 field",
    )

    rows = []
    for j in range(10):
        for k in range(19 - 2 * j):
            entry = squared_phase**j * field_parameter**k
            coefficients = [int(c) for c in entry.to_list()]
            rows.append(coefficients + [0] * (DEGREE - len(coefficients)))
    require(len(rows) == 100, "All mixed powers of weighted degree at most 18")
    require(nmod_mat(rows, PRIME).rank() == 100, "Independent reductions of the 100 mixed powers")


def main():
    variance = symmetric_moment(0, 2, 0)
    difference_fourth = symmetric_moment(0, 4, 0)
    second = moments_with_field(2, 0)
    cubic = weighted_sum(
        (variance, moments_with_field(3, 0)),
        (-3, multiply(second, moments_with_field(1, 2))),
    )
    quartic = weighted_sum(
        (variance**2, moments_with_field(4, 0)),
        (-6 * variance, multiply(second, moments_with_field(2, 2))),
        (difference_fourth, multiply(second, second)),
    )
    eliminant = correlation_eliminant(cubic, quartic, variance)
    reduced = certify_irreducible(eliminant)
    print("Correlation: degree-102 eliminant irreducible modulo 389")

    left, right = fmpq(21, 25), fmpq(17, 20)
    signs = [
        sign(value)
        for value in bernstein_coefficients(fmpq_poly(eliminant), left, right)
        if value
    ]
    require(signs[0] * signs[-1] == -1, "Opposite endpoint signs for correlation eliminant")
    require(sum(first != second for first, second in zip(signs, signs[1:])) == 1, "One Bernstein sign change")
    print("Correlation: exactly one root in (21/25, 17/20)")

    numerator = raw_moment(2, 0, 0) * raw_moment(0, 2, 0)
    denominator = raw_moment(1, 1, 0)**2
    require(numerator.degree() == 8 and denominator.degree() == 14, "Nonconstant second-moment ratio")
    print("Moments: irrational second-moment ratio under independent column rescalings")

    certify_coefficient_direction(reduced, cubic, quartic, variance, second)
    print("Direction: 100 independent mixed powers t^k (s^2)^j with k + 2j <= 18")
    print("Direction: quadratic phase extension modulo 389, giving residue degree 204")


if __name__ == "__main__":
    main()

