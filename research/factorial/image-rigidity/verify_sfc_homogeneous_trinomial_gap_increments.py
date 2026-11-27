#!/usr/bin/env python3
"""Replay factored SFC gap steps and the origin estimates."""

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from math import comb, factorial

from verify_sfc_homogeneous_trinomial_factored_sign import (
    LEFT,
    RIGHT,
    factored_negative_remainder,
)
from verify_sfc_homogeneous_trinomial_short_spans import initial_remainder


@lru_cache(maxsize=64)
def checked_quotient(a, d):
    b = d - a
    expected = Counter({(0, 1, j): 3 + 4 * (j <= a) for j in range(1, d + 1)})
    expected.update({(1, 1, j): 1 + 2 * (j <= b) for j in range(1, d + 1)})
    expected.update((0, 2, j) for j in range(1, a+1, 2))
    common, quotient = factored_negative_remainder(a, d)
    if common != expected:
        raise ArithmeticError(f"Incorrect common product at {(a, d)}.")
    return quotient, sum(n for (axis, _, _), n in common.items() if axis == 0)


def verify_increment(a, d, step_a):
    b = d - a
    previous, left_degree = checked_quotient(a, d)
    current, _ = checked_quotient(a + step_a, d + 1)
    if step_a:
        parity = int(a % 2 == 0)
        p, q = 4 - parity, 2
        numerator = 432 * (a + 1)**2 * (d + 1)**3 * (d + 3)
        denominator = 2**parity * a**2 * d**3 * (d + 2)
    else:
        p, q = 4, 4
        numerator = 432 * (b + 1)**4 * (d + 1)**3 * (d + 3)
        denominator = b**4 * d**3 * (d + 2)

    difference = denominator*current - numerator*LEFT**p*RIGHT**q*previous
    coefficients = difference.to_dict()
    negative = next((key for key, val in coefficients.items() if val < 0), None)
    if negative is not None or coefficients.get((0, 0), 0) <= 0:
        raise ArithmeticError(
            f"Factored increment failed at {(a, d, step_a)}. "
            f"Negative coefficient: {negative}. Nonpositive constant: "
            f"{coefficients.get((0, 0), 0) <= 0}."
        )

    # The least L exponent in the leading part is the pure H term
    # of (L+H)^4 ((d-2)L+(d+2)H), times the common monomial.
    old_l = 7*d + 4*a - 5 - left_degree
    old_h = max(map(sum, previous.to_dict())) - old_l
    if not previous.to_dict().get((old_l, old_h)) or (old_l+p, old_h+q) in coefficients:
        raise ArithmeticError(f"The leading coefficient ratio was not sharp at {(a, d, step_a)}.")


def origin_normalized(a, d):
    b = d - a
    k = comb(d, a)
    t = Fraction(1, comb(2*d, d))
    p, r = Fraction(k, comb(2*d, a)), Fraction(k, comb(2*d, b))
    s = Fraction(k*k, comb(2*d, 2*a))
    x = Fraction(1, comb(3*d, d))
    z = Fraction(k, comb(3*d, b))
    w0, w = Fraction(k*k, comb(3*d, 2*a)), Fraction(k*k, comb(3*d, 2*b))
    g = Fraction(k**3, comb(3*d, 3*a))
    h = Fraction(k, comb(3*d, d+a))

    aa, bb, cc = 2*(1-t), 2*(1-t+p-r), 1+s-2*r
    ff, gg = 3*(w0-w-2*h+2*z+x-1), g-3*w+3*z-1
    return aa * (bb*gg - cc*ff)


def verify_origin():
    for d in range(2, 16):
        for a in range(1, d):
            origin = origin_normalized(a, d)
            if origin <= 0:
                raise ArithmeticError(f"The origin inequality failed at {(a, d)}.")
            if d <= 6:
                scale = Fraction(
                    factorial(2*d)**2 * factorial(3*d) * factorial(d)**7,
                    comb(d, a)**4,
                )
                expanded = -initial_remainder(a, d).get((0, 0), 0)
                if origin * scale != int(expanded):
                    raise ArithmeticError(f"Origin normalization disagrees at {(a, d)}.")

    e = Fraction(1, 8192)
    if not (
        Fraction(2, 3)**26 < e
        and Fraction(2, 81) + e < Fraction(1, 40)
        and Fraction(2, 3)**8 < Fraction(1, 25)
        and Fraction(32, 141) < Fraction(23, 100)
        and Fraction(2, 3)**4 * Fraction(3, 2) * Fraction(15, 29)**2 < Fraction(2, 25)
    ):
        raise ArithmeticError("A large-span origin bound disagrees.")

    bounds = (
        (Fraction(7, 8), Fraction(187, 200), Fraction(17, 8), Fraction(28, 25)),
        (Fraction(3, 2)-2*e, Fraction(77, 100)-3*e, 3, 1+3*e),
        (Fraction(21, 16)-2*e, Fraction(65, 81)-3*e, Fraction(5, 2), 1+3*e),
        (Fraction(1, 2), Fraction(5, 9)-2*e, 1+2*e, Fraction(421, 900)),
        (Fraction(7, 8), Fraction(7, 9)-2*e, Fraction(8, 5)+2*e, Fraction(31, 25)),
        (Fraction(17, 16), Fraction(25, 27)-2*e, 2+2*e, Fraction(307, 243)),
    )
    if any(3*c*t <= b*g for c, t, b, g in bounds):
        raise ArithmeticError("The origin table does not give strict margins.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-span", type=int, default=34, help="maximum source span")
    args = parser.parse_args()
    if args.max_span < 2:
        parser.error("The maximum span must permit three offsets.")
    verify_origin()
    base, _ = checked_quotient(1, 2)
    base_coefficients = base.to_dict()
    if base_coefficients.get((0, 0), 0) <= 0 or min(base_coefficients.values()) < 0:
        raise ArithmeticError("The factored induction base failed.")
    print("Verified the origin table and every oriented pair through span 15.", flush=True)
    for d in range(2, args.max_span + 1):
        for a in range(1, d):
            verify_increment(a, d, False)
            verify_increment(a, d, True)
    print(f"Verified both factored increment bounds through source span {args.max_span}.")


if __name__ == "__main__":
    main()
