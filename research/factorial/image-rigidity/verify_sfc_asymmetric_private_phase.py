"""Exact rational margins and outward interval replay for asymmetric private degrees."""
from fractions import Fraction
from functools import lru_cache
from math import factorial, prod

from flint import arb, ctx, fmpq

ctx.prec = 96
PI = arb.pi()


def ball(lower, upper):
    midpoint = (lower + upper) / 2
    radius = (upper - lower) / 2
    return arb(fmpq(midpoint.numerator, midpoint.denominator),
               fmpq(radius.numerator, radius.denominator))


def coefficient_bounds(degree):
    log_mean = arb(factorial(degree)).log() / degree
    return tuple((2 * (arb(factorial(k)).log() - k * log_mean)).exp().upper()
                 for k in range(1, degree))


COEFFICIENTS = {degree: coefficient_bounds(degree) for degree in range(3, 41)}


@lru_cache(maxsize=None)
def radial_data(degree, lower, upper):
    q = ball(lower, upper)
    coefficients = COEFFICIENTS[degree]
    amplitudes = []
    for k, coefficient in enumerate(coefficients, start=1):
        interval = arb(coefficient / 2, coefficient / 2)
        amplitudes.append(interval * (k * q / degree).exp())
    return q.exp(), tuple(amplitudes)


def projection_excludes(baseline, amplitudes, directions):
    value = baseline
    for amplitude, direction in zip(amplitudes, directions):
        value += amplitude * direction
    return value > 0 or value < 0


def exclude_all_radii(degree, lower, upper):
    theta = PI * ball(lower, upper)
    trig = [(theta * k).sin_cos() for k in range(degree + 1)]
    sines, cosines = zip(*trig)
    projections = []
    for k in (0, degree, *range(1, degree)):
        sine_terms = [sines[j - k] if j >= k else -sines[k - j]
                      for j in range(1, degree)]
        cosine_terms = [cosines[abs(j - k)] for j in range(1, degree)]
        projections.append((-sines[k], sines[degree - k], sine_terms))
        projections.append((cosines[k], cosines[degree - k], cosine_terms))

    stack = [(Fraction(-3, 2), Fraction(3, 2), 0)]
    boxes = 0
    while stack:
        lo, hi, depth = stack.pop()
        boxes += 1
        last, amplitudes = radial_data(degree, lo, hi)
        if any(projection_excludes(constant + last * endpoint, amplitudes, directions)
               for constant, endpoint, directions in projections):
            continue
        if depth == 10:
            return False, boxes
        middle = (lo + hi) / 2
        stack.extend(((lo, middle, depth + 1), (middle, hi, depth + 1)))
    return True, boxes


def anchor_intervals(degree):
    return [(max(Fraction(0), Fraction(4*k-1, 4*degree)),
             min(Fraction(1), Fraction(4*k+1, 4*degree)))
            for k in range(1, degree+1, 2)]


def check_pair(n):
    N = n + 1
    stack = [(max(a, c), min(b, d), 0) for a, b in anchor_intervals(n)
             for c, d in anchor_intervals(N) if max(a, c) <= min(b, d)]
    angle_boxes = radial_boxes = 0
    while stack:
        lo, hi, depth = stack.pop()
        angle_boxes += 1
        lower_excluded, visited = exclude_all_radii(n, lo, hi)
        radial_boxes += visited
        if lower_excluded:
            continue
        upper_excluded, visited = exclude_all_radii(N, lo, hi)
        radial_boxes += visited
        if upper_excluded:
            continue
        if depth >= 16:
            raise AssertionError((n, lo, hi))
        middle = (lo + hi) / 2
        stack.extend(((lo, middle, depth + 1), (middle, hi, depth + 1)))
    return angle_boxes, radial_boxes


def check_degree_uniform_margins():
    Q = Fraction
    e_lower = Q(2718, 1000)
    e_upper = Q(87, 32)
    e_partial = sum(Q(1, factorial(k)) for k in range(7))
    assert e_lower < e_partial
    assert e_partial + Q(8, 7 * factorial(7)) < e_upper
    assert Q(5, 12)**40 * e_lower**39 > 40
    assert e_upper**5 * 20**4 < Q(8, 5)**40
    assert Q(5, 12)**2 * Q(40, 39) * (e_upper / 2)**2 < Q(1, 3)

    left = [factorial(h)**2 * ((e_upper / 40)**2 * Q(40, 39))**h
            for h in range(1, 7)]
    right = [(Q(5, 12)**2 * Q(40, 39))**h /
             prod(Q(40 - j, 40)**2 for j in range(h))
             for h in range(1, 7)]
    limits = ((left, Q(2, 5), Q(1, 100), Q(3, 100)),
              (right, Q(1, 3), Q(9, 40), Q(7, 25)))
    for terms, q, modulus, slope in limits:
        assert sum(terms) + q**7 / (1 - q) < modulus
        assert sum(h*c for h, c in enumerate(terms, start=1)) + (
            q**7 * (7-6*q) / (1-q)**2
        ) < slope
    assert e_lower * (1 - Q(9, 40)) > 1 + Q(1, 100)
    assert e_lower * (1 - Q(1, 100)) > 1 + Q(9, 40)
    assert Q(3, 99) + Q(28, 100) / (1 - Q(9, 40)) < Q(1, 2)
    assert Q(7, 10)**2 < Q(1, 2)
    assert 1 + Q(3, 2) + Q(9, 8) > Q(10, 3)

    r = Q(3, 8)
    real_part_with_factor = r - 6*r**3 - 120*r**5
    imag_part_with_factor = r + 6*r**3 - 120*r**5
    assert real_part_with_factor < 0 and imag_part_with_factor < 0
    assert Q(7, 10)**2 < Q(1, 2) < Q(3, 4)**2
    assert 1 - 24*r**4 + Q(7, 10)*real_part_with_factor < 0
    assert 2*r*r + Q(3, 4)*imag_part_with_factor > 0


def main():
    check_degree_uniform_margins()
    for coefficients in COEFFICIENTS.values():
        assert sum(coefficients) < arb(fmpq(7, 10))
    angle_boxes = radial_boxes = 0
    for n in range(3, 40):
        angle, radial = check_pair(n)
        angle_boxes += angle
        radial_boxes += radial
    print(f"Certified degrees 3 to 40 in {angle_boxes} angular and {radial_boxes} radial boxes.")
    print("Checked all degree-uniform rational margins and the sharpness example.")


if __name__ == "__main__":
    main()
