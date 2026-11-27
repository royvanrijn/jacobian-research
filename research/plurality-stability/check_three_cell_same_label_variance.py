"""Directed same-label variance absorption for three balanced cells.

The analytic reduction, including the positive cumulative-density
comparison, is in THREE_CELL_SAME_LABEL_VARIANCE_PENALTY.md.  All asserted
inequalities use Arb; floating point is used only in the saved summary.
"""

import argparse
from fractions import Fraction
import hashlib
import json
from math import comb, factorial
from pathlib import Path
import sys

import flint
from flint import arb, ctx

import check_three_cell_sharp_arc_row_extension as radial
from check_three_cell_global_scalar_through_half import radial_coefficients, upper_float
from check_three_cell_sharp_positive_arc_envelope import sharp_moments


MODE_ORDER = 160
ROW_ORDER = 32
RADIAL_LIMIT = 12
RADIAL_PANELS = 4096
PRECISION = 140


def ball(value):
    return arb(value.numerator) / value.denominator


def angular_certificate():
    order = 40
    mm, pp = sharp_moments(order)
    anchor = arb(3) / 2
    value, derivative = arb(0), arb(0)
    for n in range(1, order + 1):
        uniform = arb(comb(n, n // 2)) / 2**n if n % 2 == 0 else arb(0)
        coefficient = mm[n] + pp[n] - 2 * uniform
        value += coefficient * anchor**(n - 1) / factorial(n)
        if n >= 2:
            derivative += (n - 1) * coefficient * anchor**(n - 2) / factorial(n)
    error = 4 * anchor**order / factorial(order + 1) / (1 - anchor / (order + 2))
    derivative_error = 4 * order * anchor**(order - 1) / factorial(order + 1)
    derivative_error /= 1 - anchor * (order + 1) / (order * (order + 2))
    value += (-error).union(error)
    derivative += (-derivative_error).union(derivative_error)
    tangent = value - anchor * derivative
    assert derivative > 0 and tangent > arb(19) / 100
    # A larger proposed uniform constant is false at this actual anchor.
    assert value < arb(1) / 5
    coefficient = arb(81) * 19 / (64 * 100)
    minimum = arb(2) / 45 - 1 / (400 * coefficient)
    assert arb.pi()**2 < 10 and minimum > arb(17) / 500
    return {
        "exponential_anchor": "3/2",
        "series_order": order,
        "convex_tangent_lower": -upper_float(-tangent),
        "derivative_lower": -upper_float(-derivative),
        "averaged_coefficient_lower": -upper_float(-minimum),
        "claimed_averaged_coefficient": "17/500",
        "negative_control_one_fifth": "rejected at the anchor",
    }


def radial_moments(d, kappa, row):
    """GH8 interpolation, with the whole positive polynomial in each tail.

    For x >= R, Q(x) >= Q(R).  This improves the first-coefficient-only
    tail without changing the interpolation argument or discarding modes.
    """
    exponents = (1 + kappa, arb(1), 1 - kappa)
    orders = (2 * MODE_ORDER, MODE_ORDER, 0)
    output = [[arb(0) for _ in range(order + 1)] for order in orders]
    width = arb(RADIAL_LIMIT) / RADIAL_PANELS

    def endpoint(x):
        polynomial, first, second = row[-1], arb(0), arb(0)
        for coefficient in reversed(row[:-1]):
            second = second * x + 2 * first
            first = first * x + polynomial
            polynomial = polynomial * x + coefficient
        assert polynomial > 0
        factors = [8 * d * (-exponent * x**2).exp() / polynomial for exponent in exponents]
        return polynomial, first, second, factors

    left = arb(0)
    left_polynomial, _, _, left_factors = endpoint(left)
    for panel in range(RADIAL_PANELS):
        right = (panel + 1) * width
        right_polynomial, first, second, right_factors = endpoint(right)
        log_slope_bound = first / left_polynomial
        curvature_bound = second / left_polynomial
        errors = []
        for exponent, left_factor in zip(exponents, left_factors):
            assert exponent > 0
            derivative_bound = left_factor * (
                (2 * exponent * right + log_slope_bound)**2
                + 2 * exponent + curvature_bound + log_slope_bound**2
            )
            error = (derivative_bound * width**2 / 8).upper()
            errors.append((-error).union(error))
        left_power, right_power = left, right
        mass = right - left
        for n in range(2 * MODE_ORDER + 1):
            left_power *= left
            right_power *= right
            moment = (right_power - left_power) / (n + 2)
            left_weight = (right * mass - moment) / width
            right_weight = (moment - left * mass) / width
            for moments, at_left, at_right, error in zip(output, left_factors, right_factors, errors):
                if n < len(moments):
                    moments[n] += at_left * left_weight + at_right * right_weight + error * mass
            mass = moment
        left, left_polynomial, left_factors = right, right_polynomial, right_factors

    radius = arb(RADIAL_LIMIT)
    polynomial_at_radius, _, _, _ = endpoint(radius)
    for moments, exponent in zip(output, exponents):
        for n in range(len(moments)):
            decay = exponent - arb(max(0, n - 1)) / (2 * radius**2)
            assert decay > 0
            tail = 4 * d * radius**(n - 1) * (-exponent * radius**2).exp()
            tail /= polynomial_at_radius * decay
            moments[n] += arb(0).union(tail)
    tt, uu, (ww,) = output
    assert ww > 0
    return tt, uu, ww


def tail_bound(lo, hi, angular):
    rho = lo.union(hi)
    d = 1 - rho**2
    kappa = rho**2 / d
    rate = 2 * rho / d
    row, _ = radial_coefficients(d, 1 + kappa, rate, *angular, 1)
    first = MODE_ORDER + 1
    small_ratio = 4 * rate / (first + 1)
    eta = 1 + 19 * kappa / 10
    large_ratio = rate / eta
    assert 0 < small_ratio < 1 and 0 < large_ratio < 1
    small = 16 * (4 * rate)**first / (row[0] * factorial(first) * (1 - small_ratio))
    large = 3200 * large_ratio**first / (31 * d * eta * (1 - large_ratio))
    return (small + large).upper()


def moment_calibration():
    d, kappa = arb(4) / 5, arb(1) / 4
    tt, uu, ww = radial_moments(d, kappa, [arb(1)])
    tested = 0
    for moments, exponent, indices in (
        (tt, 1 + kappa, (0, 1, 2, 8, 32, 160, 320)),
        (uu, arb(1), (0, 1, 2, 8, 32, 160)),
        ([ww], 1 - kappa, (0,)),
    ):
        for n in indices:
            power = arb(n + 1) / 2
            exact = 4 * d * power.gamma() / exponent**power
            assert moments[n].contains(exact), (n, str(moments[n]), str(exact))
            tested += 1
    return {"constant_polynomial_closed_gaussian_moments": tested, "status": "pass"}


def slab_bound(lo, hi, angular):
    d = 1 - hi**2
    kappa = hi**2 / d
    row, modes = radial_coefficients(d, 1 + kappa, 2 * hi / d, *angular, ROW_ORDER)
    penalty = 17 * hi * arb.pi() * d / 1000
    normalization = arb.pi().sqrt() * d * d.sqrt() / 2
    row[0] += 2 * d * penalty / normalization
    modes[1] += penalty / normalization**2
    row = [value * (lo / hi)**n for n, value in enumerate(row, 1)]
    assert all(value > 0 for value in row) and all(value > 0 for value in modes)
    tt, uu, ww = radial_moments(d, kappa, row)
    finite = radial.gram_norm(modes, tt, uu, ww)
    contact = sum((mode * moment**2 for mode, moment in zip(modes, uu)), arb(0)) / ww
    flux_coefficients = [mode * moment / ww.sqrt() for mode, moment in zip(modes, uu)]
    flux_squared = arb(0)
    for i in range(len(modes)):
        for j in range(i + 1):
            factor = 1 if i == j else 2
            flux_squared += factor * flux_coefficients[i] * flux_coefficients[j] * (tt[i + j] - uu[i] * uu[j] / ww)
    # Q >= Q(0) gives a closed upper bound for int t^2/J_mix dnu.
    eta = 1 - kappa
    imbalance = 3 * d * arb.pi().sqrt() / (row[0] * eta**2 * eta.sqrt())
    leakage = abs((hi**2 - lo**2) / d * (imbalance / ww).sqrt()).upper()
    leak_loss = 2 * flux_squared.upper().sqrt() * leakage + contact * leakage**2
    tail = tail_bound(lo, hi, angular)
    return (finite + leak_loss + tail).upper(), finite, leak_loss, tail


def certificate(point_only=False):
    directory = Path(__file__).resolve().parent
    source_names = (
        Path(__file__).name,
        "THREE_CELL_SAME_LABEL_VARIANCE_PENALTY.md",
        "check_three_cell_sharp_arc_row_extension.py",
        "check_three_cell_sharp_positive_arc_envelope.py",
        "check_three_cell_global_scalar_through_half.py",
    )
    pins = {name: hashlib.sha256((directory / name).read_bytes()).hexdigest() for name in source_names}
    ctx.prec = PRECISION
    angular_result = angular_certificate()
    calibration = moment_calibration()
    angular = sharp_moments(MODE_ORDER)
    accepted, widths = [], set()

    def certify(left, right, depth):
        result = slab_bound(ball(left), ball(right), angular)
        if result[0] < arb(1999) / 1000:
            accepted.append(result)
            widths.add(right - left)
            return
        assert depth < 10, (str(left), str(right), str(result[0]))
        midpoint = (left + right) / 2
        certify(left, midpoint, depth + 1)
        certify(midpoint, right, depth + 1)

    if point_only:
        certify(Fraction(61, 100), Fraction(61, 100), 10)
    else:
        for numerator in range(580, 610):
            certify(Fraction(numerator, 1000), Fraction(numerator + 1, 1000), 0)
    assert all(hashlib.sha256((directory / name).read_bytes()).hexdigest() == digest for name, digest in pins.items()), "a proof or checker source changed during the run"
    return {
        "status": "pass",
        "precision_bits": ctx.prec,
        "angular": angular_result,
        "moment_calibration": calibration,
        "proof_source": source_names[1],
        "proof_source_sha256": pins[source_names[1]],
        "source_hashes": pins,
        "environment": {"python": sys.version.split()[0], "python_flint": flint.__version__},
        "assurance": {"internal_directed_verification": True, "independent_replay": False, "formal_verification": False, "external_review": False},
        "mode_order": MODE_ORDER,
        "row_order": ROW_ORDER,
        "radial_limit": RADIAL_LIMIT,
        "radial_panels_per_slab": RADIAL_PANELS,
        "rho_bounds": ["61/100", "61/100"] if point_only else ["29/50", "61/100"],
        "correlation_slabs": len(accepted),
        "slab_widths": [str(width) for width in sorted(widths)],
        "max_corrected_norm_upper": max(upper_float(row[0]) for row in accepted),
        "max_leakage_loss_upper": max(upper_float(row[2]) for row in accepted),
        "max_kernel_tail_upper": max(upper_float(row[3]) for row in accepted),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--point-only", action="store_true")
    args = parser.parse_args()
    print(json.dumps(certificate(point_only=args.point_only), indent=2))
