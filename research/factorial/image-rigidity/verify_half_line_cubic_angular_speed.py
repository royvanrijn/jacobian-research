"""Interval replay for the imaginary order speed of half-line cubic zeros."""

from flint import arb, fmpq

from verify_half_line_cubic_positive_speed import (
    absolute_moment_bounds,
    center_moments,
)


def exclude(bounds):
    """Return the strict witness excluding an upper zero of nonpositive speed."""
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

    first_error = h_z * moment[2] + h_s * moment[4] / 3
    second_error = h_z * moment[3] + h_s * moment[5] / 3
    first_lower = max(arb(0), (first.abs_lower() - first_error.abs_upper()).lower())
    if 4 * arb(y_min) * first_lower**2 > second.imag.upper() + second_error.abs_upper():
        return "speed"
    return None


def verify():
    initial = (
        (fmpq(0), fmpq(3, 8)),
        (fmpq(0), fmpq(6)),
        (fmpq(7, 5), fmpq(6)),
    )
    leaves = {"outside": 0, "G": 0, "speed": 0}
    visits = 0

    def visit(bounds, depth):
        nonlocal visits
        (s_min, s_max), (x_min, x_max), (y_min, y_max) = bounds
        if x_min**2 + y_min**2 > 36:
            leaves["outside"] += 1
            return
        visits += 1
        certificate = exclude(bounds)
        if certificate is not None:
            leaves[certificate] += 1
            return
        if depth >= 45:
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
