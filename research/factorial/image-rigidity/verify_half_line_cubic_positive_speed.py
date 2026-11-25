"""Verify the compact field exclusion in the half-line cubic speed proof."""

from math import comb, factorial

from flint import acb, arb, ctx, fmpq


ctx.prec = 128
I = acb(0, 1)
TAIL_START = arb(16)
TOLERANCE = arb("1e-11")
GAUSSIAN_MOMENTS = [
    arb(2) ** arb(fmpq(k + 1, 2)) * arb(fmpq(k + 1, 2)).gamma()
    for k in range(7)
]


def center_moments(s, x, y):
    """Enclose derivatives through order three at the exact rational center."""
    s, x, y = map(arb, (s, x, y))
    z = acb(x, y)
    saddle = 2 * z / (1 + (1 + 4 * s * z).sqrt())
    height = saddle.imag.mid()
    end = acb(TAIL_START, height)
    descent = TAIL_START - x + s * (TAIL_START**2 - height**2)
    if not height > 0 or not descent > 0:
        raise ArithmeticError("invalid center contour")

    def exponent(u):
        return -u**2 / 2 - s * u**3 / 3 + z * u

    def integrate(k):
        def integrand(u, _analytic):
            # The entire polynomial-exponential integrand is analytic on all balls.
            return u**k * exponent(u).exp()

        vertical = acb.integral(
            integrand, 0, I * height, abs_tol=TOLERANCE, rel_tol=TOLERANCE,
            eval_limit=150, depth_limit=20,
        )
        horizontal = acb.integral(
            integrand, I * height, end, abs_tol=TOLERANCE, rel_tol=TOLERANCE,
            eval_limit=250, depth_limit=20,
        )
        start_size = (TAIL_START + height)
        remainder = exponent(end).real.exp() * sum(
            comb(k, m) * start_size**(k - m) * factorial(m) / descent**(m + 1)
            for m in range(k + 1)
        )
        tail = arb(0, remainder.abs_upper())
        value = vertical + horizontal
        return acb(value.real + tail, value.imag + tail)

    return height, [integrate(k) for k in range(4)]


def absolute_moment_bounds(height, s_min, s_max, x_max, y_min):
    """Bound each absolute contour integral uniformly over an entire cell."""
    s_min, s_max, x_max, y_min = map(arb, (s_min, s_max, x_max, y_min))
    linear = x_max + s_max * height**2
    maximum = 2 * linear / (1 + (1 + 4 * s_min * linear).sqrt())
    peak = (
        -maximum**2 / 2 + linear * maximum - s_min * maximum**3 / 3
        + height**2 / 2 - y_min * height
    )
    peak_radius = (maximum**2 + height**2).sqrt()
    vertical_decay = y_min - height / 2
    if vertical_decay > 0:
        vertical = [
            min(
                (height**(k + 1) / (k + 1)).abs_upper(),
                (factorial(k) / vertical_decay**(k + 1)).abs_upper(),
            )
            for k in range(7)
        ]
    else:
        vertical_peak = (height**2 / 2 - y_min * height).abs_upper()
        vertical = [
            (height**(k + 1) * vertical_peak.exp() / (k + 1)).abs_upper()
            for k in range(7)
        ]
    return [
        vertical[k] + peak.exp() * sum(
            comb(k, m) * peak_radius**(k - m) * GAUSSIAN_MOMENTS[m]
            for m in range(k + 1)
        )
        for k in range(7)
    ]


def exclude(bounds):
    """Return the function whose nonvanishing excludes the cell, if certified."""
    (s_min, s_max), (_, x_max), (y_min, _) = bounds
    s, x, y = [(lo + hi) / 2 for lo, hi in bounds]
    h_s, h_x, h_y = [arb((hi - lo) / 2) for lo, hi in bounds]
    h_z = (h_x**2 + h_y**2).sqrt()
    height, values = center_moments(s, x, y)
    if not all(value.is_finite() for value in values):
        return None
    value, first, second, third = values
    moment = absolute_moment_bounds(height, s_min, s_max, x_max, y_min)
    quadratic = (
        h_z**2 * moment[2] / 2 + h_z * h_s * moment[4] / 3
        + h_s**2 * moment[6] / 18
    )
    real_error = h_x * abs(first.real) + h_y * abs(first.imag) + h_s * abs(third.real) / 3
    imag_error = h_x * abs(first.imag) + h_y * abs(first.real) + h_s * abs(third.imag) / 3
    if (
        value.real.abs_lower() > (real_error + quadratic).abs_upper()
        or value.imag.abs_lower() > (imag_error + quadratic).abs_upper()
    ):
        return "G"

    z = acb(arb(x), arb(y))
    speed = 4 * z * first - second
    first_error = h_z * moment[2] + h_s * moment[4] / 3
    second_error = h_z * moment[3] + h_s * moment[5] / 3
    speed_error = 4 * h_z * abs(first) + 4 * (abs(z) + h_z) * first_error + second_error
    if speed.abs_lower() > speed_error.abs_upper():
        return "D"
    return None


def verify():
    initial = (
        (fmpq(0), fmpq(3, 8)),
        (fmpq(0), fmpq(4)),
        (fmpq(7, 5), fmpq(4)),
    )
    leaves = {"outside": 0, "G": 0, "D": 0}
    visits = 0

    def visit(bounds, depth):
        nonlocal visits
        (s_min, s_max), (x_min, x_max), (y_min, y_max) = bounds
        if x_min**2 + y_min**2 > 16:
            leaves["outside"] += 1
            return
        visits += 1
        certificate = exclude(bounds)
        if certificate is not None:
            leaves[certificate] += 1
            return
        if depth >= 33:
            raise ArithmeticError(f"uncertified cell at depth {depth}: {bounds}")
        widths = (8 * (s_max - s_min), x_max - x_min, y_max - y_min)
        axis = max(range(3), key=lambda j: widths[j])
        lo, hi = bounds[axis]
        mid = (lo + hi) / 2
        for interval in ((lo, mid), (mid, hi)):
            child = list(bounds)
            child[axis] = interval
            visit(tuple(child), depth + 1)

    visit(initial, 0)
    print(f"Verified {visits} non-skipped cells. Terminal cells {leaves}.")


if __name__ == "__main__":
    verify()
