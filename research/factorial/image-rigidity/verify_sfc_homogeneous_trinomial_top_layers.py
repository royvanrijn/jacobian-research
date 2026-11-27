#!/usr/bin/env python3
"""Verify the leading height layers and the integer-gap sign certificates."""

import argparse
from collections import defaultdict
from functools import lru_cache
from math import comb, lcm

from flint import fmpq, fmpq_mpoly_ctx, fmpz_mpoly_ctx
from sympy import bernoulli

from verify_sfc_homogeneous_trinomial_factored_sign import (
    RING as HEIGHT_RING,
    count_terms,
    expand_factors,
    factored_negative_remainder,
)

RING = fmpq_mpoly_ctx.get(("a0", "b0", "x", "y"))
A0, B0, X, Y = RING.gens()
LIMIT = 24
DIFFERENCE_ORDER = 18


@lru_cache(maxsize=None)
def bernoulli_rational(index):
    numerator, denominator = bernoulli(index).as_numer_denom()
    return fmpq(int(numerator), int(denominator))


def power_sum(length, power):
    # The Bernoulli convention here has B_1 = +1/2.
    return fmpq(1, power + 1) * sum(
        bernoulli_rational(j) * comb(power + 1, j) * length ** (power + 1 - j)
        for j in range(power + 1)
    )


def normalized_series(order, y_degree, limit, shift):
    middle, gap = A0 + shift, B0 + shift
    span = middle + gap
    result = [RING.constant(0) for _ in range(limit + 1)]
    for multiplier, offset, i, j in count_terms(order, y_degree, middle, span):
        log = [RING.constant(0)]
        for power in range(1, limit + 1):
            left = (
                fmpq(1, order**power) * power_sum(order * span - offset, power)
                + (i + j) * power_sum(span, power)
                + (y_degree - j) * power_sum(middle, power)
            )
            right = (
                fmpq(1, order**power) * power_sum(offset, power)
                + (order - y_degree - i) * power_sum(span, power)
                + (y_degree - j) * power_sum(gap, power)
            )
            log.append(fmpq((-1)**(power - 1), power) * (left * X**power + right * Y**power))

        series = [RING.constant(1)]
        for degree in range(1, limit + 1):
            # Differentiate exp(log) and compare coefficients.
            series.append(fmpq(1, degree) * sum(
                j * log[j] * series[degree - j] for j in range(1, degree + 1)
            ))
        for degree, coefficient in enumerate(series):
            result[degree] += multiplier * coefficient
    return result


def multiply_series(left, right, limit):
    return [
        sum((left[i] * right[degree-i] for i in range(degree+1)
             if left[i] and right[degree-i]), RING.constant(0))
        for degree in range(limit + 1)
    ]


def remainder_layers(limit=LIMIT, shift=1):
    a, b, c = (normalized_series(2, k, limit, shift) for k in range(3))
    d, f, g = (normalized_series(3, k, limit, shift) for k in (0, 2, 3))
    if any(poly[0] for poly in (a, b, c)) or any(any(poly[:2]) for poly in (d, f, g)):
        raise ArithmeticError("Unexpected terms below order four.")
    ccd = multiply_series(multiply_series(c, c, limit), d, limit)
    abg = multiply_series(multiply_series(a, b, limit), g, limit)
    acf = multiply_series(multiply_series(a, c, limit), f, limit)
    return {
        degree: (ccd[degree] + abg[degree] - acf[degree]).to_dict()
        for degree in range(4, limit + 1)
    }


def compare_expanded_remainders(layers):
    for span in range(2, 7):
        for middle in range(1, span):
            common, quotient = factored_negative_remainder(middle, span)
            full = (expand_factors(common) * quotient).to_dict()
            scale = 2**(4*span) * 3**(3*span)
            for degree, layer in layers.items():
                specialized = defaultdict(fmpq)
                for (i, j, p, q), coefficient in layer.items():
                    exponent = (7*span + 4*middle - p, 7*span - 4*middle - q)
                    specialized[exponent] += (
                        coefficient * scale * (middle - 1)**i * (span - middle - 1)**j
                    )
                expected = {
                    exponent: coefficient for exponent, coefficient in full.items()
                    if sum(exponent) == 14*span - degree
                }
                if {key: value for key, value in specialized.items() if value} != expected:
                    raise ArithmeticError(f"Layer {degree} disagrees at {(middle, span)}.")


@lru_cache(maxsize=None)
def power_in_binomial_basis(power):
    if power == 0:
        return (1,)
    previous = power_in_binomial_basis(power - 1)
    # z*binom(z,j) = j*binom(z,j) + (j+1)*binom(z,j+1).
    return tuple(j * (
        (previous[j] if j < len(previous) else 0)
        + (previous[j-1] if j else 0)
    ) for j in range(power + 1))


def verify_integer_gap_layers(layers):
    for order in range(5, LIMIT + 1):
        converted = defaultdict(fmpq)
        for (m, n, p, q), coefficient in layers[order].items():
            for i, left in enumerate(power_in_binomial_basis(m)):
                if left:
                    for j, right in enumerate(power_in_binomial_basis(n)):
                        if right:
                            converted[i, j, p, q] += coefficient * left * right
        negative = next((key for key, value in converted.items() if value < 0), None)
        if not converted or negative is not None:
            raise ArithmeticError(f"Integer-gap layer {order} failed at {negative}.")
    print(f"Integer-gap binomial coefficients are nonnegative through order {LIMIT}.")


def differences_at_zero(samples):
    while samples:
        yield samples[0]
        samples = [upper - lower for lower, upper in zip(samples, samples[1:])]


def verify_complete_gap_differences():
    def scaled_reversed_polynomial(a0, b0):
        middle = a0 + 1
        span = a0 + b0 + 2
        common, quotient = factored_negative_remainder(middle, span)
        full = expand_factors(common) * quotient
        scale = 432 ** (DIFFERENCE_ORDER - a0 - b0)
        return HEIGHT_RING.from_dict({
            (7*span + 4*middle - p, 7*span - 4*middle - q): scale * coefficient
            for (p, q), coefficient in full.to_dict().items()
        })

    rows = [
        list(differences_at_zero([
            scaled_reversed_polynomial(a0, b0) for b0 in range(DIFFERENCE_ORDER - a0 + 1)
        ]))
        for a0 in range(DIFFERENCE_ORDER + 1)
    ]
    checked = 0
    for j in range(DIFFERENCE_ORDER + 1):
        # Later columns use progressively fewer rows of the triangular grid.
        samples = [row[j] for row in rows[:DIFFERENCE_ORDER-j+1]]
        for i, difference in enumerate(differences_at_zero(samples)):
            negative = next((key for key, value in difference.to_dict().items() if value < 0), None)
            if negative is not None:
                raise ArithmeticError(f"Complete gap difference {(i, j)} failed at {negative}.")
            checked += 1
    print(f"Verified all {checked} complete gap differences of order at most {DIFFERENCE_ORDER}.")


def matched_cutoff(coefficients, degree, grow_a):
    # One group of negative monomials per dominating highest-degree monomial.
    groups = defaultdict(list)
    for (i, j), value in coefficients.items():
        if value < 0:
            target = (degree-j, j) if grow_a else (i, degree-i)
            ratio = -value / coefficients[target]
            groups[target].append((degree-i-j, ratio))

    def sufficient(gap):
        return all(
            sum(ratio / gap**drop for drop, ratio in group) < 1
            for group in groups.values()
        )

    low, high = 0, 1
    while not sufficient(high):
        high *= 2
    while high - low > 1:
        middle = (low + high) // 2
        if sufficient(middle):
            high = middle
        else:
            low = middle
    return high


def verify_fixed_deficit(deficit):
    # Unshifted gap variables are needed for the a^2 b^4 support bound.
    layers = remainder_layers(limit=deficit, shift=0)
    by_height_power = defaultdict(dict)
    for (i, j, p, q), value in layers[deficit].items():
        by_height_power[p, q][i, j] = value

    gap_ring = fmpz_mpoly_ctx.get(("a", "b"))
    a, b = gap_ring.gens()
    degree = 2*deficit
    max_a, max_b, checked = 0, 0, 0
    for p in range(deficit + 1):
        q = deficit - p
        coefficients = by_height_power[p, q]
        if any(i < 2 or j < 4 or i+j > degree for i, j in coefficients):
            raise ArithmeticError(f"Gap support bound failed at {(p, q)}.")
        if any(coefficients.get((i, degree-i), 0) <= 0 for i in range(2, degree-3)):
            raise ArithmeticError(f"Highest gap support or sign failed at {(p, q)}.")
        cutoff_a = matched_cutoff(coefficients, degree, True)
        cutoff_b = matched_cutoff(coefficients, degree, False)
        max_a, max_b = max(cutoff_a, max_a), max(cutoff_b, max_b)

        denominator = lcm(*(int(c.denominator) for c in coefficients.values()))
        polynomial = gap_ring.from_dict({
            key: int(value*denominator) for key, value in coefficients.items()
        })
        step_a = polynomial.compose(a+1, b) - polynomial
        step_b = polynomial.compose(a, b+1) - polynomial
        candidates = (
            ((a0, b0) for a0 in range(1, cutoff_a) for b0 in range(1, a0+1)),
            ((a0, b0) for b0 in range(1, cutoff_b) for a0 in range(1, b0)),
        )
        for triangle in candidates:
            for a0, b0 in triangle:
                if step_a(a0, b0) < 0 or step_b(a0, b0) < 0:
                    raise ArithmeticError(f"Negative gap difference at {(p, q, a0, b0)}.")
                checked += 1

    print(f"Deficit {deficit}: verified both normalized differences at every gap pair.")
    print(f"Matched cutoffs: {max_a} in a-dominant region, {max_b} in b-dominant region.")
    print(f"Checked {checked} finite coefficient-gap pairs.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixed-deficit", type=int, help="verify one degree in the reciprocal heights")
    args = parser.parse_args()
    if args.fixed_deficit is not None:
        if args.fixed_deficit < 5:
            parser.error("The fixed deficit must be at least 5.")
        verify_fixed_deficit(args.fixed_deficit)
        return
    layers = remainder_layers()
    if layers[4]:
        raise ArithmeticError("Order four did not cancel.")
    for degree in range(5, 12):
        coefficients = layers[degree].values()
        if not coefficients or min(coefficients) <= 0:
            raise ArithmeticError(f"Negative or empty shifted-gap layer at order {degree}.")
        print(f"Order {degree}: all {len(coefficients)} nonzero coefficients are positive.")
    witness = (4, 0, 0, 12)
    if layers[12].get(witness) != fmpq(-676639794713, 5224277606400):
        raise ArithmeticError("Shifted-gap obstruction at order twelve disagrees.")
    verify_integer_gap_layers(layers)
    compare_expanded_remainders(layers)
    print("Order twelve has a negative shifted-gap coefficient.")
    print("All layers agree with the expanded remainders through span 6.")
    verify_complete_gap_differences()


if __name__ == "__main__":
    main()
