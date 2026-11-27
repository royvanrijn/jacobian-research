"""Directed one-direction and quartic comparison through rho = 3/5.

Run with python-flint. The radial moment enclosure is shared with the
sharp arc row checker. The slab transport and quartic test are proved in
THREE_CELL_QUARTIC_EXTENSION_THROUGH_THREE_FIFTHS.md.
"""

import argparse
from fractions import Fraction
import json
from math import factorial

from flint import arb, ctx

import check_three_cell_sharp_arc_row_extension as radial
from check_three_cell_global_scalar_through_half import radial_coefficients, upper_float
from check_three_cell_sharp_positive_arc_envelope import sharp_moments


MODE_ORDER = 136
ROW_ORDER = 24
RADIAL_LIMIT = 12
RADIAL_PANELS = 8192
COEFFICIENTS = (
    10000000000000000, -1022292365645760, -6375532343588090, -1640388695171383,
    -2685602519963690, -574400693853652, -593076167949204, -117437123001083,
    -91761013614007, -17179249713812, -11014364508137, -1967405144630,
    -1083150456344, -185639743597, -90341414090, -14918527767, -6547259999,
    -1045185363, -419766163, -64957442, -24142355, -3630008, -1259578,
    -184394, -60164,
)


def ball(value):
    return arb(value.numerator) / value.denominator


def polynomial(coefficients, x):
    value = arb(0)
    for coefficient in reversed(coefficients):
        value = value * x + coefficient
    return value


def mix_tail(lo, hi, angular):
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


def imbalance_moment(d, kappa, row):
    eta = 1 - kappa
    radius = arb(RADIAL_LIMIT)
    width = radius / 4096
    upper = arb(0)
    for index in range(4096):
        x = (index * width).union((index + 1) * width)
        integrand = 8 * d * x**4 * (-eta * x**2).exp() / polynomial(row, x)
        upper += width * integrand.upper()
    decay = eta - arb(3) / (2 * radius**2)
    assert decay > 0
    return (upper + 4 * d * radius**3 * (-eta * radius**2).exp() / (row[0] * decay)).upper()


def spectral_block(lo, hi, angular):
    d = 1 - hi**2
    kappa = hi**2 / d
    row, modes = radial_coefficients(d, 1 + kappa, 2 * hi / d, *angular, ROW_ORDER)
    row = [value * (lo / hi)**n for n, value in enumerate(row, 1)]
    tt, uu, ww = radial.radial_moments(d, kappa, 1 + kappa, row)
    size = len(modes)
    cov = [[tt[i + j] - uu[i] * uu[j] / ww for j in range(size)] for i in range(size)]

    alpha = [arb(value) / 10**16 for value in COEFFICIENTS]
    image = [sum((entry * value for entry, value in zip(line, alpha)), arb(0)) for line in cov]
    norm = sum((value * entry for value, entry in zip(alpha, image)), arb(0))
    assert norm > 0
    rayleigh = sum((mode * value**2 for mode, value in zip(modes, image)), arb(0)) / norm
    residual_coefficients = [mode * value for mode, value in zip(modes, image)]
    anchor = rayleigh.mid()
    for i, value in enumerate(alpha):
        residual_coefficients[i] -= anchor * value
    flux_coefficients = [mode * moment / ww.sqrt() for mode, moment in zip(modes, uu)]

    complement_squared = arb(0)
    residual_squared = arb(0)
    flux_squared = arb(0)
    for i in range(size):
        for j in range(i + 1):
            factor = 1 if i == j else 2
            entry = cov[i][j]
            complement_squared += factor * modes[i] * modes[j] * (entry - image[i] * image[j] / norm)**2
            residual_squared += factor * residual_coefficients[i] * residual_coefficients[j] * entry
            flux_squared += factor * flux_coefficients[i] * flux_coefficients[j] * entry

    leakage = abs((hi**2 - lo**2) / d * (imbalance_moment(d, kappa, row) / ww).sqrt()).upper()
    contact = sum((mode * moment**2 for mode, moment in zip(modes, uu)), arb(0)) / ww
    flux = flux_squared.upper().sqrt()
    loss = 2 * flux * leakage + contact * leakage**2
    tail = mix_tail(lo, hi, angular)
    e = (rayleigh - 2 + tail + loss).upper()
    residual = (residual_squared / norm).upper().sqrt().upper()
    delta = (2 - tail - complement_squared.upper().sqrt() - loss).lower()
    center = sum((value * moment for value, moment in zip(alpha, uu)), arb(0)) / ww
    return e, residual, delta, row, alpha, center, norm, kappa


def trial_geometry(lo, hi, row, alpha, center, norm, kappa):
    d_min = 1 - hi**2
    d_max = 1 - lo**2
    pole_radius = arb(1) / 4
    pole_slope = arb(13) / 10
    pole_factor = 4 * d_min / norm.sqrt()
    pole_width = pole_radius / 1000
    pole_lower = arb(2)
    for index in range(1000):
        x = (index * pole_width).union((index + 1) * pole_width)
        value = pole_factor * (polynomial(alpha, x) - center * (kappa * x**2).exp()) / polynomial(row, x)
        pole_lower = min(pole_lower, (value + pole_slope * x).lower())
    j_origin = (row[0] * (-pole_radius**2 / d_min).exp() / (4 * (d_min * d_max).sqrt())).lower()

    beta = 1 + kappa
    normalization = arb.pi().sqrt() / (2 * beta * beta.sqrt())
    width = arb(3) / 4096
    mass, second, mean = arb(0), arb(0), arb(0)
    for index in range(4096):
        x = (index * width).union((index + 1) * width)
        square = x**2
        denominator = polynomial(row, x)
        numerator = (-kappa * square).exp() * polynomial(alpha, x) - center
        mass += 2 / normalization * square * (-beta * square).exp() * width
        second += 32 * d_min**2 / (norm * normalization) * (-(1 - kappa) * square).exp() * (numerator / denominator)**2 * width
        mean += 8 * d_min / (norm.sqrt() * normalization) * x * (-square).exp() * numerator / denominator * width
    variance = (d_min / d_max * (mass * second - mean**2)).lower()
    return pole_lower, pole_slope, pole_radius, j_origin, variance


def covers_residuals(e, residual, delta, lo, hi, geometry):
    discriminant = residual**2 + delta * e
    if discriminant < 0:
        return True
    horizon = ((residual + discriminant.sqrt()) / delta).upper()
    if horizon <= 0:
        return False
    pole, slope, radius, origin_weight, variance = geometry
    if pole <= 0 or variance <= 0:
        return False
    d_min = 1 - hi**2
    density_bound = 16 * arb.pi() / (11 * lo * d_min)
    penalty_factor = hi * d_min / (2 * arb.pi())
    circumradius = (arb(2) / 3).sqrt()

    for index in range(64):
        left = index * horizon / 64
        right = (index + 1) * horizon / 64
        scalar = e + 2 * residual * right - delta * left**2
        if scalar < 0:
            continue
        denominator = 2 * origin_weight * pole**3 - 3 * slope * right**2
        cutoff = 3 * right**2 / (2 * origin_weight * pole**2)
        if denominator <= 0 or cutoff > radius:
            return False
        amplitude = 3 * circumradius * right**2 / denominator
        purity = arb(2) / 3 - density_bound * amplitude**2 * (1 + right**2)
        spread = variance.sqrt() - density_bound.sqrt() * right
        if purity <= 0 or spread <= 0 or not scalar < penalty_factor * purity * spread**2:
            return False
    return True


def certificate(point_only=False):
    ctx.prec = 140
    radial.MODE_ORDER = MODE_ORDER
    radial.RADIAL_LIMIT = RADIAL_LIMIT
    radial.RADIAL_PANELS = RADIAL_PANELS
    angular = sharp_moments(MODE_ORDER)
    widths = set()
    accepted = []

    def certify(left, right, depth):
        lo, hi = ball(left), ball(right)
        e, residual, delta, *trial = spectral_block(lo, hi, angular)
        if delta > 0:
            geometry = trial_geometry(lo, hi, *trial)
            if covers_residuals(e, residual, delta, lo, hi, geometry):
                accepted.append((e, residual, delta))
                widths.add(right - left)
                return
        assert depth < 6, (str(left), str(right), str(e), str(residual), str(delta))
        midpoint = (left + right) / 2
        certify(left, midpoint, depth + 1)
        certify(midpoint, right, depth + 1)

    if point_only:
        certify(Fraction(3, 5), Fraction(3, 5), 6)
    else:
        for numerator in range(580, 600):
            certify(Fraction(numerator, 1000), Fraction(numerator + 1, 1000), 0)

    return {
        "status": "pass",
        "precision_bits": ctx.prec,
        "radial_panels_per_slab": RADIAL_PANELS,
        "radial_limit": RADIAL_LIMIT,
        "mode_order": MODE_ORDER,
        "row_order": ROW_ORDER,
        "rho_bounds": ["3/5", "3/5"] if point_only else ["29/50", "3/5"],
        "correlation_slabs": len(accepted),
        "slab_widths": [str(width) for width in sorted(widths)],
        "max_scalar_e_upper": max(upper_float(e) for e, _, _ in accepted),
        "max_scalar_residual_upper": max(upper_float(residual) for _, residual, _ in accepted),
        "min_scalar_delta_lower": min(-upper_float(-delta) for _, _, delta in accepted),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--point-only", action="store_true")
    args = parser.parse_args()
    print(json.dumps(certificate(point_only=args.point_only), indent=2))
