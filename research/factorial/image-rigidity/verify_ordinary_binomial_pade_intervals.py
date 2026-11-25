#!/usr/bin/env python3
"""Replay the interval bounds for the ordinary binomial Padé comparison.

The moderate-angle test checks intervals in d/b at analytic endpoints.
The small-angle checks use scalar inequalities. Neither enumerates supports.
"""

from decimal import Decimal, Context, ROUND_CEILING, ROUND_FLOOR, localcontext
from math import isqrt


PRECISION = 45
DOWN = Context(prec=PRECISION, rounding=ROUND_FLOOR)
UP = Context(prec=PRECISION, rounding=ROUND_CEILING)


class Interval:
    """Closed interval with outward rounding for basic operations."""

    def __init__(self, value, upper=None):
        if isinstance(value, Interval):
            self.lower, self.upper = value.lower, value.upper
        else:
            self.lower = Decimal(value)
            self.upper = self.lower if upper is None else Decimal(upper)

    def __add__(self, other):
        other = Interval(other)
        return Interval(DOWN.add(self.lower, other.lower),
                        UP.add(self.upper, other.upper))

    __radd__ = __add__

    def __neg__(self):
        return Interval(self.upper.copy_negate(), self.lower.copy_negate())

    def __sub__(self, other):
        return self + -Interval(other)

    def __rsub__(self, other):
        return Interval(other) + -self

    def __mul__(self, other):
        other = Interval(other)
        endpoints = [(a, b) for a in (self.lower, self.upper)
                     for b in (other.lower, other.upper)]
        return Interval(min(DOWN.multiply(a, b) for a, b in endpoints),
                        max(UP.multiply(a, b) for a, b in endpoints))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = Interval(other)
        assert other.lower > 0 or other.upper < 0
        return self * Interval(DOWN.divide(Decimal(1), other.upper),
                               UP.divide(Decimal(1), other.lower))

    def __rtruediv__(self, other):
        return Interval(other) / self

    def log(self):
        assert self.lower > 0
        with localcontext(Context(prec=PRECISION)):
            return Interval(self.lower.ln().next_minus(),
                            self.upper.ln().next_plus())

    def exp(self):
        with localcontext(Context(prec=PRECISION)):
            return Interval(self.lower.exp().next_minus(),
                            self.upper.exp().next_plus())

    def sqrt(self):
        return (self.log() / 2).exp()


def largest(left, right):
    left, right = Interval(left), Interval(right)
    return Interval(max(left.lower, right.lower), max(left.upper, right.upper))


def atan_reciprocal(denominator):
    x = Interval(1) / denominator
    xx = x * x
    term = x
    total = term
    for i in range(1, 35):
        term *= -xx
        total += term / (2 * i + 1)
    # The last term is positive, so the next term bounds the downward error.
    return Interval((total - term * xx / 71).lower, total.upper)


PI = 16 * atan_reciprocal(5) - 4 * atan_reciprocal(239)


def cosine_upper(angle):
    negative_square = -(angle * angle)
    term = total = Interval(1)
    for i in range(1, 21):
        term *= negative_square / ((2 * i) * (2 * i - 1))
        total += term
    # For angles in [0, pi], the tail starts negative and is alternating.
    return Interval(total.upper)


def verify_prime_estimates():
    """Check (T27) extensions at primes and at the left endpoints."""
    estimates = {21: "-.22", 34: "-.30", 43: "-.33", 100: "-.39",
                 117: "-.425", 137: "-.449", 150: "-.45", 343: "-.49"}
    total, count = Interval(0), 0
    for n in range(2, 3001):
        prime = all(n % p for p in range(2, isqrt(n) + 1))
        if not prime and n not in estimates and n != 300:
            continue
        log_n = Interval(n).log()
        if prime:
            total += log_n / (n - 1)
            count += 1
        for start, limit in estimates.items():
            if n >= start:
                assert (total - log_n).upper < Decimal(limit), (n, start)
        if 117 <= n <= 300:
            assert (count * log_n / n).upper < Decimal("1.224"), n
    assert (1 + Interval("1.2762") / Interval(300).log()).upper < Decimal("1.224")


def g_over_s(alpha):
    right_shift = alpha * (Interval(".5") + alpha / 14)
    left_shift = alpha * (Interval(".5") - alpha / 14)
    right = -(1 - right_shift).log() - right_shift
    left = left_shift - (1 - alpha) * (1 + left_shift / (1 - alpha)).log()
    return largest(right, left) / (alpha * alpha)


def verify_proper_coset_bounds():
    rows = [(".03", ".1", ".134", "-.30", 200, ".028", "2.34"),
            (".1", Interval(1) / 3, ".160", ".20", 200, ".39", "1.60"),
            (Interval(1) / 3, ".75", ".218", ".20", 275, "1.14", ".51")]
    for left, right, kappa, excess, n0, first_margin, last_margin in rows:
        left, right, kappa = map(Interval, (left, right, kappa))
        assert g_over_s(right).upper < kappa.lower
        for n, margin in ((n0, first_margin), (20000, last_margin)):
            scale = 4 * left * n
            bound = (scale.log() / 2 - Interval("1.016")
                     - kappa * (n * right / 64).sqrt() - Interval(excess)
                     - Interval("3.001") * Interval("1.256") / scale.sqrt())
            assert bound.lower > Decimal(margin)
    for n, gap, margin in ((200, 26, ".72"), (2000, 233, "1.04")):
        u = Interval(n).log()
        log_x = (2 * u - 2 * Interval(2).log()) / 3
        bound = (u - Interval("1.016") - Interval(".00384") * gap
                 - log_x + Interval(".22") - Interval("4.001")
                 * (log_x - u).exp() * (1 + Interval("1.2762") / log_x))
        assert bound.lower > Decimal(margin)
    u = Interval(200).log()
    log_x = (2 * u + Interval(2).log()) / 3
    bound = (u - Interval("1.016") - 2 * Interval(".00384") - log_x
             + Interval(".33") - Interval("3.001") * Interval("1.256")
             * (log_x - u).exp())
    assert bound.lower > Decimal(".029")


def entropy(w):
    return (1 + w).log() + w * (1 + 1 / w).log()


def angular_gain(cosine, alpha, baseline, gap, refine_beta):
    ratio = (1 - alpha) / alpha
    if cosine.upper < 0:
        w = Interval(ratio.lower)
        return largest(-Interval("1.01") + entropy(w),
                       (ratio * (1 - cosine)).log())
    w = Interval((ratio - cosine / alpha).lower)
    if w.lower > 0:
        return largest(-Interval("1.01") + entropy(w), w.log())
    if refine_beta:
        lo_ratio = cosine / (1 - alpha)
        if lo_ratio.lower > 1:
            loss = gap / alpha * (cosine - (1 - alpha) * (1 + lo_ratio.log()))
            return largest(baseline, -Interval("1.01") - loss)
    return baseline


def verify_moderate_angles(n_min, n_max, gap_max, bins, excess, min_margin,
                           refine_beta=False):
    worst = Decimal(9)
    for gap in range(2, gap_max + 1):
        cosines = [cosine_upper(PI * j / gap) for j in range(gap)]
        for left, right in bins:
            if (gap * (1 / left - 1)).upper < 5:
                continue
            n = Interval(largest(n_min, ((gap * left - 2) / 2).exp()).lower)
            if n.lower > n_max:
                continue
            x_lower = (largest(7, gap / right) * n).sqrt()
            x_upper = largest(gap * (1 + 1 / left), (gap * n / left).sqrt())
            norm = (x_upper.log() - Interval(excess) + Interval("2.001")
                    * (1 + Interval("1.2762") / x_lower.log()) * x_upper / n)
            baseline = -Interval("1.01") - g_over_s(right) * gap * right
            gain = (n - 1).log() + sum(
                angular_gain(c, right, baseline, gap, refine_beta)
                for c in cosines
            ) / gap
            margin = (gain - norm).lower
            assert margin > Decimal(min_margin), (gap, left.lower, margin)
            worst = min(worst, margin)
    display = worst.quantize(Decimal(".00000001"), rounding=ROUND_FLOOR)
    print(f"Moderate-angle margin for N in [{n_min}, {n_max}] exceeds {display}")


def main():
    verify_prime_estimates()
    verify_proper_coset_bounds()
    fine_bins = [(Interval(k) / 200, Interval(k + (1 if k < 10 else 2)) / 200)
                 for k in [*range(6, 10), *range(10, 150, 2)]]
    coarse_bins = [(Interval(k) / 100, Interval(k + 2) / 100)
                   for k in range(3, 75, 2)]
    verify_moderate_angles(2200, 3500, 76, fine_bins, ".39", ".07", True)
    verify_moderate_angles(3500, 20000, 110, coarse_bins, ".45", ".089")
    print("Prime sums and proper-coset endpoints verified.")


if __name__ == "__main__":
    main()

