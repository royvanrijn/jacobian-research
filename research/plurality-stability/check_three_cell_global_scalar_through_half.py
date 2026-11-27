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


def angular_moments():
    pi = arb.pi()
    root_three = arb(3).sqrt()
    negative = [arb(1)]
    positive = [arb(1)]
    for mode in range(1, MODE_ORDER + 1):
        angle = mode * pi
        negative.append(4 / angle**2 if mode % 2 else arb(0))
        remainder = mode % 3
        sine = [arb(0), root_three / 2, -root_three / 2][remainder]
        positive.append(3 * sine / angle - arb(27) / (4 * angle**2) if remainder else arb(0))

    m_values = []
    p_values = []
    for n in range(MODE_ORDER + 1):
        weights = [arb(comb(n, j)) / 2**n for j in range(n + 1)]
        m_values.append(sum((w * negative[abs(n - 2 * j)] for j, w in enumerate(weights)), arb(0)))
        p_integral = sum((w * positive[abs(n - 2 * j)] for j, w in enumerate(weights)), arb(0))
        p_values.append((2 * p_integral + (-arb(1) / 2)**n) / 3)
    return m_values, p_values


def radial_coefficients(d, beta, rate, m_values, p_values):
    gamma = [arb(1), arb.pi().sqrt() / 2]
    for n in range(2, ROW_ORDER + 1):
        gamma.append(n * gamma[n - 2] / 2)

    row = []
    for n in range(1, ROW_ORDER + 1):
        power = beta ** (1 + n // 2)
        if n % 2:
            power *= beta.sqrt()
        row.append(rate**n * (m_values[n] - p_values[n]) * gamma[n] / (factorial(n) * power))
    modes = [rate**n * (m_values[n] + p_values[n]) / (2 * d * factorial(n))
             for n in range(MODE_ORDER + 1)]
    assert all(value > 0 for value in row)
    assert all(value > 0 for value in modes)
    return row, modes


def radial_moments(d, kappa, beta, row, panels):
    tt = [arb(0) for _ in range(2 * MODE_ORDER + 1)]
    uu = [arb(0) for _ in range(MODE_ORDER + 1)]
    vv = arb(0)
    width = arb(RADIAL_LIMIT) / panels
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
        for n in range(2 * MODE_ORDER + 1):
            tt[n] += t_factor * power
            if n <= MODE_ORDER:
                uu[n] += u_factor * power
            power *= x

    def tail(n, eta):
        radius = arb(RADIAL_LIMIT)
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


def operator_bound(rho, angular, panels):
    d = 1 - rho**2
    kappa = rho**2 / d
    beta = 1 + kappa
    rate = 2 * rho / d
    row, modes = radial_coefficients(d, beta, rate, *angular)
    tt, uu, vv = radial_moments(d, kappa, beta, row, panels)

    squared = arb(0)
    for i in range(MODE_ORDER + 1):
        for j in range(i + 1):
            covariance = tt[i + j] - uu[i] * uu[j] / vv
            multiplicity = 1 if i == j else 2
            squared += multiplicity * modes[i] * modes[j] * covariance**2

    # The rational correlation intervals have rho <= 1/2, so rate/beta <= 1.
    half_order = ROW_ORDER // 2
    trace_tail = 4 * beta**(half_order - 1) / row[-1]
    trace_tail *= arb(factorial(MODE_ORDER - half_order + 1))
    trace_tail /= (half_order - 1) * factorial(MODE_ORDER)
    return squared.sqrt(), trace_tail


def certificate(point_only=False):
    ctx.prec = 100
    panels = 8192
    angular = angular_moments()
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--point-only", action="store_true")
    args = parser.parse_args()
    print(json.dumps(certificate(point_only=args.point_only), indent=2))
