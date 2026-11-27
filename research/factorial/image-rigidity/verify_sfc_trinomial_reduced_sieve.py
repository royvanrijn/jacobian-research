#!/usr/bin/env python3
"""Exact controls for the total-degree factor and reduced subresultant sieve.

The written polynomial identities prove the unbounded factor and divisor
criteria. The six finite profiles below test the formulas and candidate
elimination mechanism, without extending the existing span theorem.
"""

from fractions import Fraction
from math import comb, factorial, gcd, prod

import sympy as sp

from verify_sfc_homogeneous_trinomial_factored_sign import (
    LEFT, RIGHT, RING, coefficient, expand_factors, factored_negative_remainder,
)


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def rising(start, length):
    return prod(start + k for k in range(length))


def projected_coefficients(order, middle, span):
    """Expand the moments directly, independently of the factored summands."""
    first = [rising(LEFT + 1, e) * rising(RIGHT + 1, span - e)
             for e in (0, middle, span)]
    result = [RING.constant(0) for _ in range(order + 1)]
    for i in range(order + 1):
        for j in range(order - i + 1):
            k = order - i - j
            offset = middle * j + span * k
            weight = factorial(order) // (factorial(i) * factorial(j) * factorial(k))
            term = (weight * rising(order * LEFT + 1, offset)
                    * rising(order * RIGHT + 1, order * span - offset)
                    * first[2] ** (i + j) * (-1) ** k)
            for r in range(k + 1):
                result[j + r] += term * comb(k, r) * first[0] ** (k-r) * first[1] ** r
    return result


def verify_total_degree_factor():
    count = 0
    for middle, span in ((1, 2), (1, 3), (2, 3), (2, 4)):
        total = LEFT + RIGHT + span
        for order in range(2, 6):
            for k, direct in enumerate(projected_coefficients(order, middle, span)):
                quotient, remainder = divmod(direct, total)
                require(remainder == 0, "A projected moment lacks the total-degree factor.")
                if order <= 3:
                    common, factored = coefficient(order, k, middle, span)
                    require(direct == expand_factors(common) * factored,
                            "Direct and factored projected moments disagree.")
                require(quotient * total == direct, "Incorrect polynomial division.")
                count += 1

    # Check the underlying univariate kernel formula, including both endpoints.
    ell = sp.Symbol("ell")
    for order, span in ((1, 2), (2, 3), (3, 4), (5, 2)):
        for offset in range(order * span + 1):
            left = rising(order * ell + 1, offset) * rising(
                -order * ell - order * span + 1, order * span - offset)
            right = ((-1) ** (order * span - offset)
                     * (order * ell + offset) * rising(order * ell + 1, order * span - 1))
            require(sp.expand(left - right) == 0, "The total-degree-zero kernel identity failed.")
    return count


def verify_subresultant_identity():
    aa, bb, cc, dd, ee, ff, gg, t = sp.symbols("A B C D E F G t")
    ss = cc**2 * dd + aa * bb * gg - aa * cc * ff
    tt = cc**2 * ee - bb * cc * ff + (bb**2 - aa * cc) * gg
    quad = aa + bb * t + cc * t**2
    cubic = dd + ee * t + ff * t**2 + gg * t**3
    require(sp.rem(cc**2 * cubic - ss - tt * t, quad, t) == 0,
            "The cubic remainder identity failed.")
    resultant = sp.resultant(quad, cubic, t)
    require(sp.expand(cc**2 * resultant - (cc * ss**2 - bb * ss * tt + aa * tt**2)) == 0,
            "The resultant norm formula failed.")


def marked_constants(middle, span, height):
    """Compute both divisor constants without expanding a height polynomial."""
    a, b, d, h = middle, span - middle, span, height
    ell = Fraction(-1, 3)
    baseline = Fraction(a) + ell
    p, r = rising(h + 1, b), rising(baseline + 1, b)
    u, v, z = (rising(2*h + 1, 2*b),
               rising(2*baseline + 1, b) * rising(2*h + 1, b),
               rising(2*baseline + 1, 2*b))
    face = p*p*z - 2*p*r*v + r*r*u
    require(face > 0, "The initial marked face is not positive.")
    rh = 3 ** (2*b) * face
    require(rh.denominator == 1, "The weighted face clearing is not integral.")
    anchor = rising(3*h + 1, 3*d)
    t1, t2 = rising(ell + 1, a), rising(2*ell + 1, 2*a)
    endpoint = rising(ell + 1, d)
    cmark = t1**2 * t2 * face
    dmark = anchor * endpoint**3
    kmark = (rising(ell + 1, d)**3 * rising(ell + 1, a)**4
             * rising(h + 1, d) * rising(h + 1, b)**2
             * prod(2*ell + j for j in range(1, a + 1, 2)))
    total = ell + h + d
    degree_bound = 4*d - (a+1)//2 - 3
    new = 3**degree_bound * cmark**2 * dmark / (kmark * total**3)
    require(new > 0 and new.denominator == 1, "The reduced subresultant mark is not a positive integer.")
    b1 = prod(3*j - 1 for j in range(1, a + 1))
    b2 = prod(3*j - 2 for j in range(1, 2*a + 1))
    old = anchor**2 * int(rh)**3 * b1**6 * b2**3
    total_mark = 3*(h+d) - 1
    require(old % total_mark**5 == 0, "The fifth-power resultant factor did not cancel.")
    return old, old // total_mark**5, int(new)


def candidate_count(number):
    factors = sp.factorint(number)
    factors.pop(3, None)
    all_divisors = prod(int(e) + 1 for e in factors.values())
    character_sum = prod(int(e) + 1 if int(p) % 3 == 1 else int(int(e) % 2 == 0)
                         for p, e in factors.items())
    return (all_divisors + character_sum) // 2


def verify_sieve_controls():
    controls = (
        (1, 2, 0, (108, 48, 15)),
        (1, 3, 0, (1944, 1404, 216)),
        (2, 3, 0, (1836, 1296, 144)),
        (1, 3, 1, (5616, 2496, 346)),
        (2, 4, 0, (30780, 13680, 1404)),
        (1, 4, 0, (31104, 18144, 1800)),
    )
    evaluated = 0
    for a, d, height, expected in controls:
        common, quotient = factored_negative_remainder(a, d)
        reduced, remainder = divmod(quotient, (LEFT + RIGHT + d)**3)
        require(remainder == 0, "The reduced sign polynomial is not integral.")
        terms = {(int(i), int(j)): int(v) for (i, j), v in reduced.to_dict().items()}
        if (a, d) == (1, 2):
            require(terms.get((0, 0), 0) > 0 and all(v >= 0 for v in terms.values()),
                    "The reduced gap-induction base failed.")
        bound = 4*d - (a+1)//2 - 3
        require(max(i for i, _ in terms) <= bound, "Incorrect degree bound.")
        old, norm_reduced, new = marked_constants(a, d, height)
        direct = sum(v * Fraction(-1, 3)**i * height**j for (i, j), v in terms.items())
        require(direct * 3**bound == new, "Closed and expanded subresultant marks disagree.")
        intersection = gcd(norm_reduced, new)
        candidates = [int(z) for z in sp.divisors(intersection) if z % 3 == 1]
        actual = candidate_count(old), candidate_count(norm_reduced), len(candidates)
        require(actual == expected, "The exact candidate counts changed.")
        for z in candidates:
            ell = (z-1)//3
            value = sum(v * ell**i * height**j for (i, j), v in terms.items())
            require(value != 0, "A surviving translation requires the other cubic remainder.")
        evaluated += len(candidates)
        print(f"(a,d,H)={(a,d,height)}: {actual[0]} -> {actual[1]} -> {actual[2]}; all final candidates rejected.")

    # Positivity of B*G-C*F alone is already false for an admissible support.
    quadratic = projected_coefficients(2, 1, 2)
    cubic = projected_coefficients(3, 1, 2)
    aa, bb, cc = [int(v(0, 10)) for v in quadratic]
    dd, _, ff, gg = [int(v(0, 10)) for v in cubic]
    require(bb*gg - cc*ff < 0 and cc*cc*dd + aa*(bb*gg-cc*ff) > 0,
            "The cubic-cancellation negative control failed.")
    return evaluated


def main():
    moments = verify_total_degree_factor()
    verify_subresultant_identity()
    candidates = verify_sieve_controls()
    print(f"Passed {moments} exact projected-moment factor controls and the symbolic remainder/resultant identities.")
    print(f"Rejected all {candidates} final translations in six fixed-profile mechanism controls.")
    print("No unrestricted gap inequality or all-support candidate elimination is certified.")


if __name__ == "__main__":
    main()
