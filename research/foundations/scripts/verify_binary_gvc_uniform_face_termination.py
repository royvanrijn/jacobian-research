#!/usr/bin/env python3
"""Exact regressions for uniform binary Hall and weighted-face termination.

The all-degree statements are proved in the accompanying note.  This script
checks the Hall formula on a broad finite window, exhausts the valuation
inequality on small unequal-weight lattice segments, and gives an exact
moment example in which a cancellation at the first moment is broken at a
prime dilation by the unique lower endpoint.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from math import comb, factorial


def hall_regression(limit: int = 12) -> None:
    for r in range(1, limit + 1):
        for d in range(r + 1, limit + 2):
            for mu in range(1, r + 1):
                for c in range(d + 1):
                    hall_fails = d - c < mu
                    formula = c >= d - mu + 1
                    assert hall_fails == formula


def line_points(u: int, v: int, weight: int) -> list[tuple[int, int]]:
    return [
        (a, b)
        for a in range(weight // u + 1)
        for b in range(weight // v + 1)
        if u * a + v * b == weight
    ]


def valuation_regression(limit_weight: int = 30) -> None:
    """Check the coefficient-independent inequalities in Theorem 3.1."""

    for u in range(1, 7):
        for v in range(1, 7):
            if u == v:
                continue
            for weight in range(1, limit_weight + 1):
                points = line_points(u, v, weight)
                if not points:
                    continue
                points.sort()
                for i in range(len(points)):
                    for j in range(i, len(points)):
                        for h in range(len(points)):
                            for k in range(h, len(points)):
                                lo = max(i, h)
                                hi = min(j, k)
                                if lo > hi:
                                    continue
                                intersection = points[lo : hi + 1]
                                alpha = min(intersection, key=lambda z: z[0] + z[1])
                                s = sum(alpha)
                                p = max(
                                    5,
                                    1
                                    + max(
                                        coordinate
                                        for point in points
                                        for coordinate in point
                                    ),
                                )
                                # It is enough to check all lattice points in
                                # the dilated intersection.  Non-p-multiples
                                # receive two Frobenius factors.
                                left_a, right_a = points[i], points[j]
                                left_b, right_b = points[h], points[k]
                                min_x = max(p * left_a[0], p * left_b[0])
                                max_x = min(p * right_a[0], p * right_b[0])
                                for x_exp in range(min_x, max_x + 1):
                                    remaining = p * weight - u * x_exp
                                    if remaining < 0 or remaining % v:
                                        continue
                                    y_exp = remaining // v
                                    if x_exp % p == 0 and y_exp % p == 0:
                                        beta = (x_exp // p, y_exp // p)
                                        lower_bound = sum(beta)
                                        if beta == alpha:
                                            assert lower_bound == s
                                        else:
                                            assert lower_bound >= s + 1
                                    else:
                                        factorial_value = x_exp // p + y_exp // p
                                        assert factorial_value + 2 >= s + 1


def multiply(
    left: dict[tuple[int, int], int],
    right: dict[tuple[int, int], int],
) -> dict[tuple[int, int], int]:
    answer: dict[tuple[int, int], int] = defaultdict(int)
    for (a, b), c in left.items():
        for (i, j), d in right.items():
            answer[a + i, b + j] += c * d
    return dict(answer)


def power(
    polynomial: dict[tuple[int, int], int], exponent: int
) -> dict[tuple[int, int], int]:
    answer = {(0, 0): 1}
    base = polynomial
    n = exponent
    while n:
        if n & 1:
            answer = multiply(answer, base)
        base = multiply(base, base)
        n //= 2
    return answer


def scalar_moment(
    operator: dict[tuple[int, int], int],
    polynomial: dict[tuple[int, int], int],
    exponent: int,
) -> int:
    op_power = power(operator, exponent)
    poly_power = power(polynomial, exponent)
    return sum(
        factorial(a) * factorial(b) * coefficient * poly_power.get((a, b), 0)
        for (a, b), coefficient in op_power.items()
    )


def translate(
    polynomial: dict[tuple[int, int], int], point: tuple[int, int]
) -> dict[tuple[int, int], int]:
    answer: dict[tuple[int, int], int] = defaultdict(int)
    zx, zy = point
    for (a, b), coefficient in polynomial.items():
        for i in range(a + 1):
            for j in range(b + 1):
                answer[i, j] += (
                    coefficient
                    * comb(a, i)
                    * comb(b, j)
                    * zx ** (a - i)
                    * zy ** (b - j)
                )
    return {exponent: coefficient for exponent, coefficient in answer.items() if coefficient}


def differential_value(
    operator: dict[tuple[int, int], int],
    polynomial: dict[tuple[int, int], int],
    exponent: int,
    point: tuple[int, int],
) -> int:
    op_power = power(operator, exponent)
    poly_power = power(polynomial, exponent)
    zx, zy = point
    answer = 0
    for (a, b), op_coefficient in op_power.items():
        for (i, j), poly_coefficient in poly_power.items():
            if i < a or j < b:
                continue
            falling_x = factorial(i) // factorial(i - a)
            falling_y = factorial(j) // factorial(j - b)
            answer += (
                op_coefficient
                * poly_coefficient
                * falling_x
                * falling_y
                * zx ** (i - a)
                * zy ** (j - b)
            )
    return answer


def ray_moment(
    operator: dict[tuple[int, int], int],
    shifted_polynomial: dict[tuple[int, int], int],
    delta: tuple[int, int],
    exponent: int,
) -> int:
    """Coefficient on the output ray exponent * delta.

    ``shifted_polynomial`` stores B after subtracting one copy of delta
    from every exponent.
    """

    op_power = power(operator, exponent)
    poly_power = power(shifted_polynomial, exponent)
    dx, dy = exponent * delta[0], exponent * delta[1]
    answer = 0
    for (a, b), coefficient in op_power.items():
        other = poly_power.get((a, b), 0)
        if not other:
            continue
        falling_x = factorial(a + dx) // factorial(dx)
        falling_y = factorial(b + dy) // factorial(dy)
        answer += falling_x * falling_y * coefficient * other
    return answer


def valuation(number: int, prime: int) -> int:
    assert number
    answer = 0
    while number % prime == 0:
        answer += 1
        number //= prime
    return answer


def prime_endpoint_example() -> None:
    # Weights (3,2), common weight 12.  The first moment cancels between
    # (4,0) and (2,3), but (4,0) is the unique least-ordinary-degree point
    # of the Newton-segment intersection.
    operator = {(0, 6): 1, (2, 3): 1, (4, 0): 1}
    polynomial = {(2, 3): -2, (4, 0): 1}
    assert scalar_moment(operator, polynomial, 1) == 0
    fifth = scalar_moment(operator, polynomial, 5)
    assert fifth != 0
    assert valuation(fifth, 5) == 4

    # The shifted-ray version has carrier delta=(1,2).  Its first ray
    # coefficient cancels, while the seventh again isolates (4,0).
    shifted_polynomial = {(2, 3): -1, (4, 0): 3}
    delta = (1, 2)
    assert ray_moment(operator, shifted_polynomial, delta, 1) == 0
    shifted_seventh = ray_moment(operator, shifted_polynomial, delta, 7)
    assert shifted_seventh != 0
    assert valuation(shifted_seventh, 7) == 4


def translated_multiradial_identity() -> None:
    operator = {(1, 0): 2, (0, 2): -1}
    polynomial = {(2, 1): 3, (0, 3): -2, (1, 0): 5}
    point = (2, -1)
    translated = translate(polynomial, point)
    for exponent in range(1, 5):
        assert scalar_moment(operator, translated, exponent) == differential_value(
            operator, polynomial, exponent, point
        )


def homogeneous_factorial_and_channel_regression() -> None:
    """Check the beta identity and the two compatibility warnings."""

    # C(U,V)=2 U^3-3 U^2 V+5 U V^2-7 V^3.
    polynomial = {(3, 0): 2, (2, 1): -3, (1, 2): 5, (0, 3): -7}
    for exponent in range(1, 6):
        expanded = power(polynomial, exponent)
        factorial_moment = sum(
            coefficient * factorial(a) * factorial(b)
            for (a, b), coefficient in expanded.items()
        )
        # Integral of U^a(1-U)^b is a!b!/(a+b+1)!.
        beta_moment = sum(
            Fraction(
                coefficient * factorial(a) * factorial(b),
                factorial(a + b + 1),
            )
            for (a, b), coefficient in expanded.items()
        )
        assert factorial_moment == factorial(3 * exponent + 1) * beta_moment

    # For G=S U+S^{-1}V, the S-constant term is zero in odd powers and
    # binom(2k,k) U^k V^k in power 2k.
    def channel_constant_term(exponent: int) -> dict[tuple[int, int], int]:
        if exponent % 2:
            return {}
        half = exponent // 2
        return {(half, half): comb(exponent, half)}

    first_channel = channel_constant_term(1)
    second_channel = channel_constant_term(2)
    assert not first_channel
    assert second_channel == {(1, 1): 2}
    assert second_channel != power(first_channel, 2)

    # The toric exponent map (a,b)->(a+b,b) changes factorial weights by
    # a vector-dependent factor even on one ordinary-degree line.
    def blowup_ratio(a: int, b: int) -> Fraction:
        old = factorial(a) * factorial(b)
        new = factorial(a + b) * factorial(b)
        return Fraction(new, old)

    assert blowup_ratio(2, 0) == 1
    assert blowup_ratio(1, 1) == 2


def minimal_bernstein_hall_circuit() -> None:
    """Verify Theorem 5.2 and Long's forbidden-character warning."""

    def pure_moment(a: int, b: int, c: int, d: int, exponent: int) -> Fraction:
        # Angular balance forces the two binomial selection counts to agree.
        return sum(
            Fraction(
                comb(exponent, k)
                * (a * c) ** (exponent - k)
                * (b * d) ** k,
                exponent + 1,
            )
            for k in range(exponent + 1)
        )

    samples = [
        (1, 1, 1, -1),
        (2, -3, 5, 7),
        (-4, 1, 3, 2),
    ]
    for a, b, c, d in samples:
        for exponent in range(1, 10):
            assert pure_moment(a, b, c, d, exponent) == Fraction(
                (a * c + b * d) ** exponent,
                exponent + 1,
            )

    # Long's normalized circuit has ac+bd=0.  Its multiplier Z^{-1}
    # lies outside the polynomial-multiplier cone and leaves a boundary
    # term (-1)^(m-1)/(m+1).
    for exponent in range(1, 15):
        assert pure_moment(1, 1, 1, -1, exponent) == 0
        long_mixed = Fraction((-1) ** (exponent - 1), exponent + 1)
        assert long_mixed


def primitive_cusp_parallelogram() -> None:
    """Verify the exact all-degree obstruction in Theorem 5.3."""

    def central(n: int) -> int:
        return comb(2 * n, n)

    def trinomial(n: int) -> int:
        return factorial(3 * n) // factorial(n) ** 3

    def obstruction(r: int, s: int) -> int:
        # Twice E_{r,s}, avoiding fractions.
        cr, cs = central(r), central(s)
        return (
            2 * (trinomial(r) - trinomial(s))
            + 9 * (cs - cr) * (cr + cs - 2)
        )

    assert obstruction(1, 2) == 48
    assert obstruction(1, 3) == -108
    assert obstruction(2, 3) == -156
    for r in range(1, 25):
        for s in range(1, 25):
            assert (obstruction(r, s) == 0) == (r == s)

    # Replay (5.22) directly after solving the first two moments.
    # The coefficient [U^k V^k](1+U+V+tUV)^m is computed sparsely.
    for r in range(1, 8):
        for s in range(1, 8):
            if r == s:
                continue
            cr, cs = central(r), central(s)
            t = Fraction(4, cr + cs - 4)
            q = Fraction(-factorial(r), factorial(s)) / t
            polynomial = {
                (0, 0): Fraction(1),
                (1, 0): Fraction(1),
                (0, 1): Fraction(1),
                (1, 1): t,
            }
            moments: list[Fraction] = []
            for exponent in range(1, 4):
                diagonal = power(polynomial, exponent)
                moment = sum(
                    Fraction(comb(exponent, k))
                    * q**k
                    * factorial(r * (exponent - k))
                    * factorial(s * k)
                    * diagonal.get((k, k), 0)
                    for k in range(exponent + 1)
                )
                moments.append(moment)
            assert moments[0] == 0
            assert moments[1] == 0
            assert 2 * moments[2] == factorial(r) ** 3 * obstruction(r, s)

    # Sparse affine dilations 1+U+V+t U^p V^q have no off-axis
    # contribution through moment two unless (p,q)=(1,1).
    for r in range(1, 12):
        for s in range(1, 12):
            for p in range(1, 6):
                for q_power in range(1, 6):
                    if (p, q_power) == (1, 1):
                        continue
                    t = Fraction(3, 7)
                    q_operator = Fraction(-factorial(r), factorial(s)) / t
                    polynomial = {
                        (0, 0): Fraction(1),
                        (1, 0): Fraction(1),
                        (0, 1): Fraction(1),
                        (p, q_power): t,
                    }
                    diagonal = power(polynomial, 2)
                    second = sum(
                        Fraction(comb(2, k))
                        * q_operator**k
                        * factorial(r * (2 - k))
                        * factorial(s * k)
                        * diagonal.get((p * k, q_power * k), 0)
                        for k in range(3)
                    )
                    expected = factorial(r) ** 2 * (
                        central(r) + central(s) - 4
                    )
                    assert second == expected
                    assert (second == 0) == (r == s == 1)


def five_channel_minor_inheritance_warning() -> None:
    """Verify the exact finite-prefix obstruction in Example 5.5.

    This is not a pure-moment-zero pair: its fourth scalar moment is
    nonzero.  It shows that vanishing through the first three moments of a
    five-channel convolution need not descend to any four-channel minor.
    """

    operator = {
        (2, 0): Fraction(1),
        (0, 3): Fraction(-1, 3),
    }
    polynomial = {
        (2, 0): Fraction(1),
        (0, 3): Fraction(1),
        (0, 0): Fraction(13, 30),
        (1, 0): Fraction(11, 2),
        (1, 3): Fraction(1),
    }
    moments = [
        scalar_moment(operator, polynomial, exponent)
        for exponent in range(1, 5)
    ]
    assert moments == [0, 0, 0, 1_205_760]

    support = tuple(polynomial)
    for omitted in support:
        minor = {
            exponent: coefficient
            for exponent, coefficient in polynomial.items()
            if exponent != omitted
        }
        minor_moments = [
            scalar_moment(operator, minor, exponent)
            for exponent in range(1, 4)
        ]
        assert any(minor_moments)


def unit_line_half_bridge_pivot() -> None:
    """Verify Theorem 5.6 and its fourth-moment obstruction."""

    def invariants(n: int) -> tuple[int, int, int]:
        central = comb(2 * n, n)
        trinomial = factorial(3 * n) // factorial(n) ** 3
        quadrinomial = factorial(4 * n) // factorial(n) ** 4
        return central, trinomial, quadrinomial

    def obstruction(n: int) -> int:
        central, trinomial, quadrinomial = invariants(n)
        return (
            2 * quadrinomial
            - 40 * central * trinomial
            + 48 * trinomial
            + 81 * central**3
            - 180 * central**2
            + 132 * central
            - 48
        )

    assert obstruction(2) == -480
    previous_ratio: Fraction | None = None
    for n in range(4, 31):
        central, trinomial, quadrinomial = invariants(n)
        ratio = Fraction(quadrinomial, central * trinomial)
        assert ratio > 20
        if previous_ratio is not None:
            assert ratio > previous_ratio
        previous_ratio = ratio
        assert obstruction(n) > 0

    # Replay the normalized half-bridge moments directly.
    for n in range(2, 15, 2):
        central, trinomial, _ = invariants(n)
        u = Fraction(central - 2, 4)
        v = u * (central - 2) - Fraction(
            trinomial - 9 * central + 12,
            18,
        )
        assert u and v
        operator = {
            (1, 0): Fraction(1),
            (0, n): Fraction(-1, factorial(n)),
        }
        polynomial = {
            (1, 0): Fraction(1),
            (0, n): Fraction(1),
            (0, 0): v,
            (0, n // 2): u,
            (1, n // 2): Fraction(1),
        }
        moments = [
            scalar_moment(operator, polynomial, exponent)
            for exponent in range(1, 5)
        ]
        assert moments[:3] == [0, 0, 0]
        assert moments[3] == Fraction(obstruction(n), 2)
        assert moments[3]


def eight_five_channel_obstructions(limit: int = 30) -> None:
    """Replay Corollary 5.10 on an exact endpoint-order window."""

    def weights(
        r: int, s: int
    ) -> tuple[dict[tuple[int, int], Fraction], dict[int, Fraction]]:
        operator_ratio = Fraction(-factorial(r), factorial(s))
        w = {
            (moment_order, k): Fraction(comb(moment_order, k))
            * operator_ratio**k
            * factorial(r * (moment_order - k))
            * factorial(s * k)
            for moment_order in range(2, 5)
            for k in range(moment_order + 1)
        }
        endpoint = {
            moment_order: sum(
                Fraction(comb(moment_order, k)) * w[moment_order, k]
                for k in range(moment_order + 1)
            )
            for moment_order in range(2, 5)
        }
        return w, endpoint

    def quadratic_cubic(
        w: dict[tuple[int, int], Fraction],
        endpoint: dict[int, Fraction],
        q_level: int,
        cubic_level: int,
    ) -> Fraction:
        u = -endpoint[2] / (2 * w[2, q_level])
        s_q = w[3, q_level] + w[3, q_level + 1]
        v = -(endpoint[3] + 6 * s_q * u) / (3 * w[3, cubic_level])
        a_q = (
            w[4, q_level]
            + 2 * w[4, q_level + 1]
            + w[4, q_level + 2]
        )
        return (
            endpoint[4]
            + 12 * a_q * u
            + 6 * w[4, 2 * q_level] * u**2
            + 12 * (w[4, cubic_level] + w[4, cubic_level + 1]) * v
        )

    def double_quadratic(
        w: dict[tuple[int, int], Fraction],
        endpoint: dict[int, Fraction],
        first_level: int,
        second_level: int,
    ) -> tuple[Fraction, Fraction]:
        a = 2 * w[2, first_level]
        b = 2 * w[2, second_level]
        c = 6 * (w[3, first_level] + w[3, first_level + 1])
        d = 6 * (w[3, second_level] + w[3, second_level + 1])
        determinant = a * d - b * c
        assert determinant
        u = (-endpoint[2] * d + b * endpoint[3]) / determinant
        v = (-a * endpoint[3] + c * endpoint[2]) / determinant
        a_first = (
            w[4, first_level]
            + 2 * w[4, first_level + 1]
            + w[4, first_level + 2]
        )
        a_second = (
            w[4, second_level]
            + 2 * w[4, second_level + 1]
            + w[4, second_level + 2]
        )
        obstruction = (
            endpoint[4]
            + 12 * a_first * u
            + 6 * w[4, 2 * first_level] * u**2
            + 12 * a_second * v
            + 6 * w[4, 2 * second_level] * v**2
            + 12 * w[4, first_level + second_level] * u * v
        )
        return determinant, obstruction

    quadratic_cubic_types = (
        (0, 0),
        (0, 1),
        (0, 2),
        (0, 3),
        (1, 0),
        (1, 1),
    )
    double_quadratic_types = ((0, 1), (0, 2))
    for r in range(1, limit + 1):
        for s in range(1, limit + 1):
            if r == s:
                continue
            w, endpoint = weights(r, s)
            for q_level, cubic_level in quadratic_cubic_types:
                assert quadratic_cubic(
                    w,
                    endpoint,
                    q_level,
                    cubic_level,
                )
            for first_level, second_level in double_quadratic_types:
                determinant, obstruction = double_quadratic(
                    w,
                    endpoint,
                    first_level,
                    second_level,
                )
                assert determinant
                assert obstruction

    w_12, endpoint_12 = weights(1, 2)
    assert quadratic_cubic(w_12, endpoint_12, 1, 1) == -240


def opposite_three_by_three_packet(limit: int = 20) -> None:
    """Check the exact central-binomial dichotomy in Theorem 7.5."""

    central = lambda n: comb(2 * n, n)
    for degree in range(3, limit + 1):
        for first in range(1, degree - 1):
            for second in range(1, degree - first):
                third = degree - first - second
                profile_product = (
                    central(first) * central(second) * central(third)
                )
                assert central(degree) > profile_product
                for endpoint_order in range(1, limit + 1):
                    ratio = Fraction(
                        central(degree * endpoint_order),
                        central(degree) * profile_product,
                    )
                    if endpoint_order == 1:
                        assert ratio == Fraction(1, profile_product)
                        assert ratio <= Fraction(1, 8)
                    else:
                        assert ratio > 1


def radial_digit_spectrum(limit: int = 100) -> None:
    """Check the unordered-pair rigidity in Lemma 7.4 bis."""

    for total in range(limit + 1):
        products: dict[int, tuple[int, int]] = {}
        for first in range(total // 2 + 1):
            second = total - first
            product = factorial(first) * factorial(second)
            assert product not in products
            products[product] = (first, second)


def main() -> None:
    hall_regression()
    valuation_regression()
    prime_endpoint_example()
    translated_multiradial_identity()
    homogeneous_factorial_and_channel_regression()
    minimal_bernstein_hall_circuit()
    primitive_cusp_parallelogram()
    five_channel_minor_inheritance_warning()
    unit_line_half_bridge_pivot()
    eight_five_channel_obstructions()
    radial_digit_spectrum()
    opposite_three_by_three_packet()
    print("PASS: uniform Hall, weighted-face, and beta-Hall circuit regressions")


if __name__ == "__main__":
    main()
