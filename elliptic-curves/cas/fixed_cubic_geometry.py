"""Exact generic geometry of the theta + u*theta^2 family (over Q).

This is not a rank predictor for its specializations. The rank computation
uses the written Shioda--Tate/height proof in RANK_JUMP_REASSESSMENT_2026-09-05.md.
All polynomial identities below are checked over Q[u], not at sampled u.
"""
from fractions import Fraction as Q
from math import isqrt


def rational(value):
    if isinstance(value, (float, bool)):
        raise ValueError("exact rational input required")
    return Q(value)


def poly(values):
    values = list(map(rational, values)) or [Q(0)]
    while len(values) > 1 and not values[-1]:
        values.pop()
    return tuple(values)


def add(*rows):
    out = [Q(0)] * max(map(len, rows))
    for row in rows:
        for i, value in enumerate(row):
            out[i] += value
    return poly(out)


def scale(row, value):
    return poly(value * x for x in row)


def mul(left, right):
    out = [Q(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i+j] += a*b
    return poly(out)


def power(row, exponent):
    out = (Q(1),)
    for _ in range(exponent):
        out = mul(out, row)
    return out


def remainder(left, right):
    if right == (0,):
        raise ZeroDivisionError("zero polynomial")
    left = poly(left)
    while left != (0,) and len(left) >= len(right):
        shift = len(left)-len(right)
        left = add(left, scale((Q(0),)*shift+right, -left[-1]/right[-1]))
    return left


def coprime(left, right):
    while right != (0,):
        left, right = right, remainder(left, right)
    return len(left) == 1 and left[0] != 0


def value(row, u):
    out = Q(0)
    for coefficient in reversed(row):
        out = out*u + coefficient
    return out


def family_polynomials(A, B):
    """Return a2(u), a4(u), a6(u) for y^2=x^3+a2*x^2+a4*x+a6."""
    A, B = map(rational, (A, B))
    return poly((0, 2*A)), poly((A, 3*B, A*A)), poly((B, 0, A*B, -B*B))


def evaluate_cubic(A, B, u, x):
    A, B, u, x = map(rational, (A, B, u, x))
    a2, a4, a6 = family_polynomials(A, B)
    return x**3 + value(a2, u)*x*x + value(a4, u)*x + value(a6, u)


def rational_square(value):
    value = rational(value)
    return value >= 0 and isqrt(value.numerator)**2 == value.numerator and isqrt(value.denominator)**2 == value.denominator


def generic_geometry(A, B):
    """Generic ranks for B != 0, 4*A^3+27*B^2 != 0, characteristic zero.

    Geometrically: rational elliptic surface, three I2 fibres and one I0*.
    The rank-one free Q-vector space has the quadratic character of B,
    witnessed by (B*u^2, sqrt(B)*(1+A*u^2+B*u^3)), of height 1/2.
    A nonsquare B therefore gives Q(u)-rank ZERO, not a transported high rank.
    No upper bound for any E_u(Q) is returned.
    """
    A, B = map(rational, (A, B))
    discriminant_f = -4*A**3-27*B**2
    if not B or not discriminant_f:
        raise ValueError("the 3I2+I0* proof requires B != 0 and a separable cubic")
    a2, a4, a6 = family_polynomials(A, B)
    D = poly((1, 0, A, B))
    discriminant = add(mul(power(a2, 2), power(a4, 2)), scale(power(a4, 3), -4),
                       scale(mul(power(a2, 3), a6), -4), scale(power(a6, 2), -27),
                       scale(mul(mul(a2, a4), a6), 18))
    c4 = scale(add(power(a2, 2), scale(a4, -3)), 16)
    x_section = poly((0, 0, B))
    section_rhs = add(power(x_section, 3), mul(a2, power(x_section, 2)), mul(a4, x_section), a6)
    derivative_at_section = add(scale(power(x_section, 2), 3), scale(mul(a2, x_section), 2), a4)
    checks = {
        "disc_F_equals_disc_f_times_D_squared": discriminant == scale(power(D, 2), discriminant_f),
        "D_separable": coprime(D, poly((0, 2*A, 3*B))),
        "c4_coprime_to_D": coprime(c4, D),
        "section_identity": section_rhs == scale(power(D, 2), B),
        "section_at_each_finite_node": derivative_at_section == mul(poly((A, 3*B)), D),
        "infinity_discriminant_order_six": 12-(len(discriminant)-1) == 6,
        "infinity_c4_order_at_least_two": 4-(len(c4)-1) >= 2,
    }
    if not all(checks.values()):
        raise ArithmeticError("generic-family polynomial proof failed")
    return {
        "schema": "elliptic-curves.fixed-cubic-generic-geometry.v1",
        "A": str(A), "B": str(B), "disc_f": str(discriminant_f),
        "checks": checks, "geometric_fibres": ["I2", "I2", "I2", "I0*"],
        "chi": 1, "geometric_picard_rank": 10, "fibre_root_rank": 7,
        "geometric_generic_rank": 1, "section_height": "1/2",
        "B_is_rational_square": rational_square(B),
        "arithmetic_generic_rank": int(rational_square(B)),
        "rank_scope": "Q(u); not a bound for any specialized E_u(Q)",
        "proof_note": "elliptic-curves/notes/RANK_JUMP_REASSESSMENT_2026-09-05.md#the-fixed-field-family-has-generic-rank-zero",
    }


def negative_square_section(A, B, v):
    """One rational point after u=-v^2; not a transported anchor class.

    The identity is F_u(-A*u-1/u) = -(1+A*u^2+B*u^3)^2/u^3.
    A rational witness does not imply independence from other supplied points.
    """
    A, B, v = map(rational, (A, B, v))
    generic_geometry(A, B)
    if not v:
        raise ValueError("this affine section requires v != 0")
    u = -v*v
    D = 1 + A*u*u + B*u**3
    if not D:
        raise ValueError("singular specialization")
    x, y = A*v*v + 1/(v*v), D/(v**3)
    if evaluate_cubic(A, B, u, x) != y*y:
        raise ArithmeticError("base-change section identity failed")
    return {"u": str(u), "x": str(x), "y": str(y)}


def alternating_rank_distribution(dimension):
    """Uniform alternating F2 matrices, NOT a proved elliptic-curve model.

    When a row/column is appended to a form of radical dimension k, rank
    increases by two with probability 1-2^-k, otherwise stays unchanged.
    No independence assumption across elliptic curves, p-value or rank bound.
    """
    if type(dimension) is not int or not 0 <= dimension <= 256:
        raise ValueError("dimension must be an integer between 0 and 256")
    distribution = {0: Q(1)}
    for n in range(dimension):
        nxt = {}
        for rank, probability in distribution.items():
            stay = Q(1, 2**(n-rank))
            nxt[rank] = nxt.get(rank, Q(0)) + probability*stay
            if stay != 1:
                nxt[rank+2] = nxt.get(rank+2, Q(0)) + probability*(1-stay)
        distribution = nxt
    return distribution
