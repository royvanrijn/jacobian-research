#!/usr/bin/env python3
"""Replay the interval bounds for the ordinary binomial Padé comparison.

The moderate-angle test checks intervals in d/b at analytic endpoints.
The small-angle checks use scalar inequalities. Neither enumerates supports.
"""

from decimal import Decimal, Context, ROUND_CEILING, ROUND_FLOOR, localcontext
from math import factorial, isqrt


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

    def __pow__(self, exponent):
        assert isinstance(exponent, int) and exponent >= 0
        result = Interval(1)
        for _ in range(exponent):
            result *= self
        return result

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
    """Check the small-argument prime estimates at primes and at the left endpoints."""
    estimates = {21: "-.22", 34: "-.30", 43: "-.33", 100: "-.39",
                 117: "-.425", 137: "-.449", 150: "-.45", 343: "-.49"}
    total, count = Interval(0), 0
    for n in range(2, 3001):
        prime = all(n % p for p in range(2, isqrt(n) + 1))
        if not prime and n not in estimates and n not in (40, 300):
            continue
        log_n = Interval(n).log()
        if prime:
            total += log_n / (n - 1)
            count += 1
        for start, limit in estimates.items():
            if n >= start:
                assert (total - log_n).upper < Decimal(limit), (n, start)
        if 40 <= n <= 300:
            assert (count * log_n / n).upper < Decimal("1.256"), n
        if 117 <= n <= 300:
            assert (count * log_n / n).upper < Decimal("1.224"), n
    # For x >= 300 the Dusart bound is smaller than both constants.
    assert (1 + Interval("1.2762") / Interval(300).log()).upper < Decimal("1.224")


def g_over_s(alpha):
    right_shift = alpha * (Interval(".5") + alpha / 14)
    left_shift = alpha * (Interval(".5") - alpha / 14)
    right = -(1 - right_shift).log() - right_shift
    left = left_shift - (1 - alpha) * (1 + left_shift / (1 - alpha)).log()
    return largest(right, left) / (alpha * alpha)


def verify_ratio_bounds():
    """Check both rows of the ratio lemma at their lowest terminal."""
    q = Interval
    rows = [
        (25, 3, ".224", "1.218", ".52", q(7) / 8, 10),
        (200, 2, ".368", "1.033", ".24", q(19) / 20, 13),
    ]
    for n, c, a, growth, tail, eta, gap in rows:
        a, growth, tail = map(q, (a, growth, tail))
        log_n = q(n).log()
        sigma = 2 * log_n + c
        assert (-q(c) / 2).exp().upper < a.lower
        assert ((log_n + q(c) / 2) / (n - 1)).exp().upper < growth.lower
        assert (q(1).exp() * a * growth**5 / 5).upper < tail.lower
        assert (2 * (1 - sigma / 4).exp()).upper < tail.lower
        short_terms = sum((a * growth**k)**k / factorial(k) for k in range(1, 5))
        coefficient_sum = 2 * short_terms + 2 * tail**5 / (1 - tail)
        assert coefficient_sum.upper < eta.lower
        assert sigma.lower > gap - 1
        annulus_margin = q(gap * (n - 1)) / (2 * (2*n - 1)) + (1 - eta).log()
        assert annulus_margin.lower > 0
    print("Both annulus rows verified at their first terminals.")


def verify_pade_loss_and_angle_guards():
    """Check the constants shared between the Padé lemma and Appendix D."""
    q = Interval
    caps = [
        (q(".03"), q(".128")),
        (q(".1"), q(".134")),
        (q(1) / 3, q(".160")),
        (q(".5"), q(".176")),
        (q(".75"), q(7) / 32),
        (q(".75"), q(".218")),
    ]
    for alpha, cap in caps:
        assert g_over_s(alpha).upper < cap.lower
    # arcsin(sqrt(.03)) < .111*pi/2 from the alternating sine series.
    t = q(".111") * PI / 2
    assert (t - t**3 / 6).lower > q(".03").sqrt().upper
    assert (1 / q(".97").sqrt()).upper < q("1.0156").lower
    assert (4 * q(2).log() / PI * q(".03").sqrt()).upper < q(".153").lower
    assert (2 * q("1.006") / PI).upper < q(".641").lower
    # The normalized beta maximum costs at most log(1+W)+1 < 5.
    assert ((q(100) / 3).log() + 1).upper < 5
    first_m = 2199
    stirling_slack = q(".01") * (first_m + 1) + (2 * PI / first_m).log() / 2 - 5
    assert stirling_slack.lower > 0
    assert (1 / q(first_m)).upper < q(".001").lower


def verify_intermediate_scalar_bounds():
    """Compare the bounds in the intermediate-order appendix."""
    q = Interval
    log2 = q(2).log()
    log3 = q(3).log()
    min_t = (q(100) / 3).log()

    def above(value, margin):
        assert value.lower > Decimal(margin), (value.lower, margin)

    def below(value, bound):
        assert value.upper < Decimal(bound), (value.upper, bound)

    u0 = q(2200).log()
    for left, right, kappa, excess, margin in [
        (q(".03"), q(".1"), q(".134"), q("-.30"), "1.35"),
        (q(".1"), q(1) / 3, q(".160"), q(".20"), ".75"),
        (q(1) / 3, q(".75"), q(".218"), q(".20"), ".14"),
    ]:
        assert g_over_s(right).upper < kappa.lower
        x = (q(2200) / (4 * left)).sqrt()
        gain = u0 - q(".006") - q(2) / 3 * (q("1.01") + kappa * (2 * u0 + 2))
        gain += ((1 - right) / right).log() / 3
        charge = x.log() + excess + q("3.001") * q("1.256") * x / 2200
        above(gain - charge, margin)
    for d in (2, 3):
        if d == 2:
            log_x = (2 * u0 + log2) / 3
            gain = u0 - q("1.016") - q(7) / 32 * q("1.5")
            charge = log_x - q(".33") + q("3.001") * q("1.256") * (log_x - u0).exp()
            above(gain - charge, ".95")
        else:
            log_x = (2 * u0 + u0.log()) / 4
            gain = u0 - q("1.016") - q(21) / 128 * u0
            charge = log_x - q(".22") + q("4.001") * q("1.256") * (log_x - u0).exp()
            above(gain - charge, "1.10")
    for n, dmax, margin in [(2200, 233, "1.08"), (20000, 565, ".68")]:
        u = q(n).log()
        log_x = (2 * u - 2 * log2) / 3
        gain = u - q("1.016") - q(".00384") * dmax
        charge = log_x - q(".22") + q("4.001") * (log_x - u).exp() * (1 + q("1.2762") / log_x)
        above(gain - charge, margin)

    below(-log2 - q(".33") + q("3.001") * q("1.256") / 2, ".862")
    above(-q(".006") + (((q(".97") / q(".03")).log() - q("1.02536")) / 2), "1.21")
    above(q(".9") * (min_t - q(".7237")) - q(".006") - q(".10868")
          - q(".111") * (q("3.18") + log2) - q(".153"), "1.80")

    def angle_constant(size, gap, small):
        g0 = q("1.01") + q(".00384") * gap
        if small:
            angular = 2 * log2 * q(".03") * (size * size - 1) / (3 * size)
            omitted = q(0)
        else:
            angular = q(".153")
            omitted = q(".111") * (g0 + log2)
        return (-q(".006") + (1 - q(1) / size) * (min_t - q(".7237"))
                + (2 * (2 * size / PI).log() - g0) / size - omitted - angular)

    for d in range(6, 19, 2):
        above(angle_constant(d // 2, d, True), "1.89")
    for d in range(4,20):
        above(angle_constant(d, d, d < 10), "2.24" if d < 10 else "2.41")
    for n, dmax, margin in [(2200,233,".17"),(20000,565,".08")]:
        gain = (q(".95") * (min_t - q(".7237")) - q(".002") - q(".05434")
                - q(".111") * (q("1.01") + q(".00384") * dmax + log2) - q(".153"))
        x = q("1.03") * (q(20) * n / q(".03")).sqrt()
        charge = (q("1.03").log() - q(".49") + q("2.001") * q("1.03")
                  * (1 + q("1.2762") / x.log()))
        above(gain - charge, margin)
    for x, correction, bound in [(343,".49","2.052"),(137,".449","2.177")]:
        charge = (q("1.03").log() - q(correction) + q("2.001") * q("1.03")
                  * (1 + q("1.2762") / q(x).log()))
        below(charge, bound)
    above((u0 - log3) / 2 - q(".2") - q("1.03").log() + q(".425")
          - q("2.001") * q("1.03") * q("1.224"), ".96")

    def min_upper(x, y):
        return q(min(x.lower, y.lower), min(x.upper, y.upper))

    table = [(".03",".036",".128",".086",".12"),
             (".036",".045",".129",".096",".20"),
             (".045",".07",".131",".120",".11"),
             (".07",".11",".134",".151",".21"),
             (".11",".17",".140",".189",".51"),
             (".17",".28",".150",".245",".67"),
             (".28",".45",".170",".315",".62"),
             (".45",".60",".190",".370",".53"),
             (".60",".75",".218",".420",".08")]
    high_worst = Decimal(9)
    for left,right,kappa,eta,margin in table:
        left,right,kappa,eta = map(q,(left,right,kappa,eta))
        assert g_over_s(right).upper < kappa.lower
        assert cosine_upper(PI*eta).upper < (1-right).lower
        for n,dmax in [(20000,257),(32768,300),(65536,565)]:
            u = q(n).log()
            sigma = 2*u+2
            s = min_upper(sigma,dmax*right)
            b = min_upper(sigma/(left*left),dmax/left)
            gain = (u - q("1.011") - kappa*(eta*s+right)
                    + q(".45")*largest(0,q("1.01")+((1-right)/right).log()))
            charge = (u+b.log())/2-q(".49")+q("2.46")*(b/n).sqrt()
            above(gain - charge, margin)
            high_worst = min(high_worst, (gain-charge).lower)
            above(1-2*kappa*eta-q(".5")-1/sigma, "0")
    u = q(20000).log()
    gain = u-q("1.011")-q(".218")*(q(".42")*q("6.75")+q(".75"))
    charge = (u+q(300).log())/2-q(".49")+q("2.46")*(q(300)/20000).sqrt()
    above(gain-charge, "0")
    display = high_worst.quantize(Decimal(".00000001"), rounding=ROUND_FLOOR)
    print(f"High-order moderate-angle table margin exceeds {display}.")
    print("Intermediate scalar and endpoint comparisons passed.")



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
    if not __debug__:
        raise RuntimeError("run the exact checks without Python's -O option")
    verify_ratio_bounds()
    verify_pade_loss_and_angle_guards()
    verify_prime_estimates()
    verify_intermediate_scalar_bounds()
    fine_bins = [(Interval(k) / 200, Interval(k + (1 if k < 10 else 2)) / 200)
                 for k in [*range(6, 10), *range(10, 150, 2)]]
    coarse_bins = [(Interval(k) / 100, Interval(k + 2) / 100)
                   for k in range(3, 75, 2)]
    verify_moderate_angles(2200, 3500, 76, fine_bins, ".39", ".07", True)
    verify_moderate_angles(3500, 20000, 110, coarse_bins, ".45", ".089")
    print("Prime sums and proper-coset endpoints verified.")


if __name__ == "__main__":
    main()

