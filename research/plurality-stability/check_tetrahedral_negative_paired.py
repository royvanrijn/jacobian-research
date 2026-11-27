#!/usr/bin/env python3
"""Certify the negative tetrahedral paired interface at correlation -1/3."""

from math import comb
import json
import platform

import flint
from flint import arb, ctx, fmpq


ORDER = 32
RADIUS = 12


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def multiply(left, right):
    result = [fmpq(0)] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            result[i + j] += x * y
    return result


def polynomial_envelope(a, b, variance_factor):
    """Bound (1+(at+b)^3)P(at+b)^2 exp(-t^2/(2*variance_factor))."""
    p = [fmpq(1), fmpq(19, 100), fmpq(11, 100), fmpq(3, 200)]
    f = multiply([fmpq(1), fmpq(0), fmpq(0), fmpq(1)], multiply(p, p))
    moments = [arb(1)] + [
        (arb(fmpq(3 * j * variance_factor, 8))).sqrt()**j
        for j in range(1, len(f))
    ]
    return sum((arb(c * comb(i, j) * a**j * b**(i-j)) * moments[j]
                for i, c in enumerate(f) for j in range(i + 1)), arb(0))


def approximation_error():
    """Check the Bernstein ellipse bound and the radial exterior bounds."""
    gamma = arb(1) / 5
    h = arb(9) / 16
    radius_envelope = polynomial_envelope(fmpq(1), fmpq(9, 2), 1)
    real_envelope = polynomial_envelope(fmpq(1), fmpq(0), 1)
    angle_envelope = polynomial_envelope(fmpq(8, 5), fmpq(0), 1)
    radial_exp = 8 * arb(27)**2 / 4 * (h * arb(81) / 4
                    + gamma**2 / (2 * (h - arb(3) / 8))).exp()
    radial_poly = 26 * radius_envelope * real_envelope * (
        h * arb(81) / 4 + 8 * gamma**2).exp()
    angle_gamma = gamma * arb(13) / 10
    angle_exp = 8 * RADIUS**2 * ((angle_gamma**2 + gamma**2) /
                  (4 * (h - arb(3) / 8 * arb(13) / 10))).exp()
    angle_poly = 26 * angle_envelope * real_envelope * (
        4 * (angle_gamma**2 + gamma**2)).exp()
    m = arb(4 * 10**11)
    require(radial_exp + radial_poly < m, "Radial ellipse bound failed")
    require(angle_exp + angle_poly < m, "Angular ellipse bound failed")
    require((arb(3) / 4).cosh() < arb(13) / 10, "Angular real-part bound failed")
    require((arb(3) / 2).cosh().sqrt() < arb(8) / 5, "Angular norm bound failed")
    quadrature = 32 * 4 * RADIUS**2 * m / arb(2)**(2 * ORDER)

    tail_envelope = polynomial_envelope(fmpq(1), fmpq(0), 2)
    require(tail_envelope < 100, "Outer polynomial bound failed")
    q = arb(RADIUS)**2
    exponential_tail = 1728 * (-q / 6).exp()
    polynomial_tail = 4 * 52 * tail_envelope**2 * (arb(8) / RADIUS) * (-q / 4).exp()
    return quadrature + exponential_tail + polynomial_tail


def integrals(theta, gamma):
    roots = [arb.legendre_p_root(ORDER, j, weight=True) for j in range(ORDER)]
    nodes = []
    h = arb(9) / 16
    for r0, wr0 in roots:
        r = RADIUS * (r0 + 1) / 2
        wr = RADIUS * wr0 / 2
        for theta0, wt0 in roots:
            c, s = (theta * theta0).cos(), (theta * theta0).sin()
            u, t = r * c, r * s
            weight = r * wr * theta * wt0 * (-h * r*r).exp()
            polynomial = 1 - 19*u/100 + 3*u*u/40 + 7*t*t/200 - u**3/100 - u*t*t/200
            nodes.append((u, t, weight * (gamma * u).exp(), polynomial, weight))

    half = arb(1) / 2
    root_half = half.sqrt()
    tau = arb(3) / 8
    standard = arb(87) / 80
    rotation = -arb(49) / 120
    result = arb(0)
    mean = arb(0)
    for i, (u, t, w, p, w0) in enumerate(nodes):
        for j in range(i + 1):
            v, z, q, p2, q0 = nodes[j]
            pp = p * p2
            uz, tv, tz = u*z, t*v, t*z
            self_dot = u*v + tz
            adjacent_dot = root_half * (uz + tv) + half * tz
            adjacent_alt = half * (adjacent_dot - (u - root_half*t)*(v - root_half*z))
            same = (-tau * self_dot).exp()
            adjacent = (-tau * adjacent_dot).exp()
            e0 = 2 * same - (2*standard + rotation*self_dot)*pp
            e1 = adjacent - (standard + rotation*adjacent_alt)*pp
            e2 = rotation * tz * pp
            summand = (6*e0*e0 + 24*e1*e1 + 6*e2*e2) * w * q
            result += summand if i == j else 2 * summand
            mean += (2 if i == j else 4) * (same + adjacent) * w0 * q0

    prefactor = 1 / (64 * (arb(8) / 9)**2 * (arb(1) / 3)**2 * theta**2)
    norm2_error = 36 * prefactor * approximation_error()

    mean_ellipse = 4 * arb(27)**2 / 4 * (h * arb(81) / 4).exp()
    require(mean_ellipse < 10**8 and 4 * RADIUS**2 < 10**8, "Mean ellipse bound failed")
    mean_error = (32 * 4 * RADIUS**2 * arb(10**8) / arb(2)**(2 * ORDER)
                  + arb(512) / 9 * (-arb(3) * RADIUS**2 / 8).exp())
    mean_factor = 1 / (16 * (arb(8) / 9) * theta**2)
    return prefactor * result, norm2_error, mean_factor * (mean + mean_error)


def certify():
    ctx.prec = 96
    theta = arb(2).sqrt().atan()
    gamma = theta.sin() * arb.pi().sqrt() / (8 * theta)
    require(0 < theta and theta < 1, "Invalid facet angle")
    require(0 < gamma and gamma < arb(1)/5, "Invalid Jensen linear coefficient")
    integral, error, mean_upper = integrals(theta, gamma)
    require(integral + error < (arb(33)/50)**2, "Rank-six residual bound failed")
    require(mean_upper < arb(784)/1000, "Normalized weight mean bound failed")
    b = arb(7) / 16
    z = gamma / (2*b.sqrt())
    reciprocal_upper = (1 + z*arb.pi().sqrt()*(z*z).exp()*(1 + z.erf())) / (2*b)
    require(reciprocal_upper < arb(1487)/1000, "Normalized reciprocal weight mean bound failed")
    angle_upper = 1 - 1/(mean_upper * reciprocal_upper)
    require(angle_upper < arb(143)/1000, "Squared weighted balance angle bound failed")
    require(arb(33)/50 + (3 - arb(33)/50)*(arb(143)/1000) < 1,
            "Balanced gap bound failed")
    return {
        "command": ".venv/bin/python research/plurality-stability/check_tetrahedral_negative_paired.py",
        "python": platform.python_version(),
        "python_flint": flint.__version__,
        "precision_bits": ctx.prec,
        "gauss_legendre_order": ORDER,
        "radial_cutoff": RADIUS,
        "rank_six_residual_hilbert_schmidt_squared_quadrature": integral.str(18),
        "absolute_quadrature_and_exterior_error": error.str(18),
        "rank_six_residual_hilbert_schmidt_norm_bound": "33/50",
        "normalized_weight_mean_upper": mean_upper.str(18),
        "normalized_reciprocal_weight_mean_upper": reciprocal_upper.str(18),
        "weighted_balance_angle_squared_upper": angle_upper.str(18),
    }


if __name__ == "__main__":
    print(json.dumps(certify(), indent=2))
