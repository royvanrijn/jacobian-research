"""Directed radial bound for the unrestricted three-cell scalar comparison.

Run with python-flint. THREE_CELL_GLOBAL_SCALAR_THROUGH_HALF.md gives the
angular inequality, the exact moment formulas and both infinite tail bounds.
"""

import argparse
import json
from math import comb, factorial, inf, nextafter

from flint import arb, ctx


MODE_ORDER = 48
ROW_ORDER = 12
RADIAL_LIMIT = 10


def upper_float(value):
    return nextafter(float(value.upper()), inf)


def angular_moments(mode_order):
    pi = arb.pi()
    root_three = arb(3).sqrt()
    negative = [arb(1)]
    positive = [arb(1)]
    for mode in range(1, mode_order + 1):
        angle = mode * pi
        negative.append(4 / angle**2 if mode % 2 else arb(0))
        remainder = mode % 3
        sine = [arb(0), root_three / 2, -root_three / 2][remainder]
        positive.append(3 * sine / angle - arb(27) / (4 * angle**2) if remainder else arb(0))

    m_values = []
    p_values = []
    for n in range(mode_order + 1):
        weights = [arb(comb(n, j)) / 2**n for j in range(n + 1)]
        m_values.append(sum((w * negative[abs(n - 2 * j)] for j, w in enumerate(weights)), arb(0)))
        p_integral = sum((w * positive[abs(n - 2 * j)] for j, w in enumerate(weights)), arb(0))
        p_values.append((2 * p_integral + (-arb(1) / 2)**n) / 3)
    return m_values, p_values


def radial_coefficients(d, beta, rate, m_values, p_values, row_order):
    gamma = [arb(1), arb.pi().sqrt() / 2]
    for n in range(2, row_order + 1):
        gamma.append(n * gamma[n - 2] / 2)

    row = []
    for n in range(1, row_order + 1):
        power = beta ** (1 + n // 2)
        if n % 2:
            power *= beta.sqrt()
        row.append(rate**n * (m_values[n] - p_values[n]) * gamma[n] / (factorial(n) * power))
    modes = [rate**n * (m_values[n] + p_values[n]) / (2 * d * factorial(n))
             for n in range(len(m_values))]
    assert all(value > 0 for value in row)
    assert all(value > 0 for value in modes)
    return row, modes


def radial_moments(d, kappa, beta, row, panels, mode_order, radial_limit):
    tt = [arb(0) for _ in range(2 * mode_order + 1)]
    uu = [arb(0) for _ in range(mode_order + 1)]
    vv = arb(0)
    width = arb(radial_limit) / panels
    for index in range(panels):
        x = (index * width).union((index + 1) * width)
        square = x**2
        denominator = sum((value * x**n for n, value in enumerate(row)), arb(0))
        assert denominator > 0
        common = 8 * d * width / denominator
        t_factor = common * (-beta * square).exp()
        u_factor = common * (-square).exp()
        vv += common * (-(1 - kappa) * square).exp()
        power = arb(1)
        for n in range(2 * mode_order + 1):
            tt[n] += t_factor * power
            if n <= mode_order:
                uu[n] += u_factor * power
            power *= x

    def tail(n, eta):
        radius = arb(radial_limit)
        power = arb(max(0, n - 1)) / 2
        decay = eta - power / radius**2
        assert decay > 0
        upper = 4 * d * radius**(n - 1) * (-eta * radius**2).exp() / (row[0] * decay)
        return arb(0).union(upper)

    tt = [value + tail(n, beta) for n, value in enumerate(tt)]
    uu = [value + tail(n, arb(1)) for n, value in enumerate(uu)]
    vv += tail(0, 1 - kappa)
    assert vv > 0
    return tt, uu, vv


def finite_bound(rho, angular, panels, row_order, radial_limit):
    d = 1 - rho**2
    kappa = rho**2 / d
    beta = 1 + kappa
    rate = 2 * rho / d
    row, modes = radial_coefficients(d, beta, rate, *angular, row_order)
    mode_order = len(modes) - 1
    tt, uu, vv = radial_moments(d, kappa, beta, row, panels, mode_order, radial_limit)

    squared = arb(0)
    for i in range(mode_order + 1):
        for j in range(i + 1):
            covariance = tt[i + j] - uu[i] * uu[j] / vv
            multiplicity = 1 if i == j else 2
            squared += multiplicity * modes[i] * modes[j] * covariance**2

    return squared.sqrt(), row, d, kappa, rate


def operator_bound(rho, angular, panels):
    finite, row, _, kappa, _ = finite_bound(rho, angular, panels, ROW_ORDER, RADIAL_LIMIT)
    # The rational correlation intervals have rho <= 1/2, so rate/beta <= 1.
    beta = 1 + kappa
    half_order = ROW_ORDER // 2
    trace_tail = 4 * beta**(half_order - 1) / row[-1]
    trace_tail *= arb(factorial(MODE_ORDER - half_order + 1))
    trace_tail /= (half_order - 1) * factorial(MODE_ORDER)
    return finite, trace_tail


def certificate(point_only=False):
    ctx.prec = 100
    panels = 8192
    angular = angular_moments(MODE_ORDER)
    finite, tail = operator_bound(arb(1) / 2, angular, panels)
    assert finite + tail < arb(193) / 100
    result = {
        "status": "pass",
        "precision_bits": ctx.prec,
        "radial_panels_per_correlation_ball": panels,
        "mode_order": MODE_ORDER,
        "row_order": ROW_ORDER,
        "point": {
            "rho": "1/2",
            "finite_upper": upper_float(finite),
            "tail_upper": upper_float(tail),
            "total_upper": upper_float(finite + tail),
        },
    }
    if point_only:
        return result

    maximum_upper = 0.0
    for index in range(1760, 2000):
        rho = (arb(index) / 4000).union(arb(index + 1) / 4000)
        finite, tail = operator_bound(rho, angular, panels)
        assert finite + tail < arb(194) / 100, index
        maximum_upper = max(maximum_upper, upper_float(finite + tail))
    result["interval"] = {
        "rho_bounds": ["11/25", "1/2"],
        "correlation_slabs": 240,
        "slab_width": "1/4000",
        "total_upper": maximum_upper,
    }
    return result


def extension_certificate():
    ctx.prec = 100
    panels = 16384
    mode_order = 80
    row_order = 16
    radial_limit = 12
    angular = angular_moments(mode_order)
    maximum_upper = 0.0
    maximum_tail = 0.0
    for start, stop, denominator in [(2000, 2080, 4000), (4160, 4240, 8000)]:
        for index in range(start, stop):
            rho = (arb(index) / denominator).union(arb(index + 1) / denominator)
            finite, row, d, kappa, rate = finite_bound(rho, angular, panels, row_order, radial_limit)

            # Use J_1 for x <= 2 and the exponential lower bound for J for x >= 2.
            small_rate = 4 * rate
            small_ratio = small_rate / (mode_order + 2)
            assert small_ratio < 1
            small_tail = 16 * small_rate**(mode_order + 1)
            small_tail /= row[0] * factorial(mode_order + 1) * (1 - small_ratio)

            eta = 1 + 7 * kappa / 4
            large_ratio = rate / eta
            assert large_ratio < 1
            large_tail = 288 * large_ratio**(mode_order + 1) / (7 * d * eta * (1 - large_ratio))
            tail = small_tail + large_tail
            assert finite + tail < arb(199) / 100, (index, denominator)
            maximum_upper = max(maximum_upper, upper_float(finite + tail))
            maximum_tail = max(maximum_tail, upper_float(tail))

    return {
        "status": "pass",
        "precision_bits": ctx.prec,
        "radial_panels_per_correlation_ball": panels,
        "radial_limit": radial_limit,
        "mode_order": mode_order,
        "row_order": row_order,
        "interval": {
            "rho_bounds": ["1/2", "53/100"],
            "correlation_slabs": 160,
            "slab_widths": ["1/4000", "1/8000"],
            "tail_upper": maximum_tail,
            "total_upper": maximum_upper,
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--point-only", action="store_true")
    group.add_argument("--extension", action="store_true")
    args = parser.parse_args()
    output = extension_certificate() if args.extension else certificate(point_only=args.point_only)
    print(json.dumps(output, indent=2))
