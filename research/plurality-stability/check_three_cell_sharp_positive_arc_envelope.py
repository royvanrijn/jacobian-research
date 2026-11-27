"""Directed radial comparison with the sharp aligned-arc coefficient.

Run with python-flint. The angular proof and the exponential row bound
are in THREE_CELL_SHARP_POSITIVE_ARC_ENVELOPE.md. The finite Gram bound
reuses the through-half checker's radial moment enclosure.
"""

import json
from math import comb, factorial

from flint import arb, ctx

from check_three_cell_global_scalar_through_half import angular_moments, finite_bound, upper_float


def sharp_moments(mode_order):
    m_values, _ = angular_moments(mode_order)
    pi = arb.pi()
    integrals = [arb(1)]
    for mode in range(1, mode_order + 1):
        cos_value = arb(1) if mode % 3 == 0 else -arb(1) / 2
        numerator = 9 * cos_value - 3 - 6 * (-1)**mode
        integrals.append(numerator / (mode * pi)**2)

    p_values = []
    for n in range(mode_order + 1):
        terms = (arb(comb(n, j)) * integrals[abs(n - 2 * j)] for j in range(n + 1))
        p_values.append(sum(terms, arb(0)) / 2**n)
    return m_values, p_values


def certificate():
    ctx.prec = 100
    assert (arb.pi() / 10).cos()**2 > arb(9) / 10
    assert (2 * arb.pi() / 9).cos()**2 < arb(3) / 5
    panels = 16384
    mode_order = 96
    row_order = 24
    radial_limit = 12
    angular = sharp_moments(mode_order)
    maximum_upper = 0.0
    maximum_tail = 0.0

    for start, stop, denominator in [(2120, 2160, 4000), (4320, 4400, 8000)]:
        for index in range(start, stop):
            rho = (arb(index) / denominator).union(arb(index + 1) / denominator)
            finite, row, d, kappa, rate = finite_bound(rho, angular, panels, row_order, radial_limit)

            small_rate = 4 * rate
            small_ratio = small_rate / (mode_order + 2)
            assert small_ratio < 1
            small_tail = 16 * small_rate**(mode_order + 1)
            small_tail /= row[0] * factorial(mode_order + 1) * (1 - small_ratio)

            eta = 1 + 19 * kappa / 10
            large_ratio = rate / eta
            assert large_ratio < 1
            large_tail = 3200 * large_ratio**(mode_order + 1) / (31 * d * eta * (1 - large_ratio))
            tail = small_tail + large_tail
            assert finite + tail < arb(99) / 50, (index, denominator)
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
            "rho_bounds": ["53/100", "11/20"],
            "correlation_slabs": 120,
            "slab_widths": ["1/4000", "1/8000"],
            "tail_upper": maximum_tail,
            "total_upper": maximum_upper,
        },
    }


if __name__ == "__main__":
    print(json.dumps(certificate(), indent=2))
