"""Certify the low-field half-line cubic zeros for -1 <= T <= 1."""

from functools import lru_cache
from math import factorial

from flint import acb, arb, ctx, fmpq

from verify_half_line_cubic_zero_temperature import winding

ctx.prec = 256
TERMS_T = 80
TERMS_H = 420
FIELD_RADIUS = 20
TOP_SLOPE = fmpq(17, 20)
LEFT = fmpq(-1, 4)
RIGHT = fmpq(15)
STEPS = 128


def initial_series():
    one_third = arb(fmpq(1, 3))
    two_thirds = arb(fmpq(2, 3))
    moments = [
        one_third.gamma() / arb(3) ** two_thirds,
        two_thirds.gamma() / arb(3) ** one_third,
        arb(1),
    ]
    for k in range(3, 2 * (TERMS_T + 3)):
        moments.append((k - 2) * moments[k - 3])
    return tuple(
        tuple(
            moments[j + 2 * index] / (2**index * factorial(index))
            for index in range(TERMS_T + 3)
        )
        for j in (0, 1)
    )


INITIAL = initial_series()


@lru_cache(None)
def field_series(temperature):
    t = arb(temperature)
    coefficients = []
    for row in INITIAL:
        value = arb(0)
        for coefficient in reversed(row[:TERMS_T]):
            value = -t * value + coefficient
        tail = 2 * sum(row[TERMS_T:])
        coefficients.append(value + arb(0, tail.abs_upper()))
    coefficients.append((1 - t * coefficients[1]) / 2)
    for k in range(3, TERMS_H):
        coefficients.append(
            (coefficients[k - 3] - t * (k - 1) * coefficients[k - 1])
            / (k * (k - 1))
        )
    # Subsequent blocks of three shrink by at least 1/4 in this weighted norm.
    end = max(
        (abs(coefficients[k]) * FIELD_RADIUS**k).abs_upper()
        for k in range(TERMS_H - 3, TERMS_H)
    )
    tail = arb(0, end)
    derivative_tail = arb(0, (2 * (TERMS_H + 5) * end / FIELD_RADIUS).abs_upper())
    return coefficients, tail, derivative_tail


def evaluate(temperature, z):
    if not abs(z).abs_upper() < FIELD_RADIUS:
        raise ArithmeticError("Point outside the coefficient-tail disk")
    coefficients, tail, derivative_tail = field_series(temperature)
    value, derivative = acb(0), acb(0)
    for coefficient in reversed(coefficients):
        derivative = z * derivative + value
        value = z * value + coefficient
    return (
        acb(value.real + tail, value.imag + tail),
        acb(derivative.real + derivative_tail, derivative.imag + derivative_tail),
    )


def variation_bounds(temperature, z, value, derivative, temp_radius, field_radius):
    """Bound value, affine remainder, and derivative on a parameter cell."""
    t, dt, dh = arb(temperature), arb(temp_radius), arb(field_radius)
    p_abs, h_abs = abs(derivative), abs(z)
    norm = abs(value) + p_abs + 1
    temp_norm = h_abs + 1
    value_t = (t * derivative - z * value - 1) / 2
    derivative_t = (
        (t * z - 1) * value - (t * t + z) * derivative + t
    ) / 2
    temp_variation = ((temp_norm * dt).exp() - 1) * norm
    value_temp_error = dt * (
        abs(value_t) + (h_abs + 1) * temp_variation / 2 + dt * p_abs / 2
    )
    derivative_temp_error = dt * (
        abs(derivative_t) + temp_norm * temp_variation + dt * (h_abs + 2) * norm / 2
    )

    moved_norm = norm + temp_variation
    field_variation = (((h_abs + dh + 2) * dh).exp() - 1) * moved_norm
    second = (
        abs(z * value + 1 - t * derivative)
        + (h_abs + 1) * (temp_variation + field_variation)
        + dh * (moved_norm + field_variation)
        + dt * p_abs
    )
    bend = dh * derivative_temp_error + dh * dh * second / 2
    return (
        value_temp_error + dh * p_abs + bend,
        abs(value) + value_temp_error + bend,
        derivative_temp_error + dh * second,
    )


def certify_top():
    accepted = 0

    def certify_cell(tlo, thi, xlo, xhi, depth):
        nonlocal accepted
        if depth > 32:
            raise ArithmeticError("Top subdivision reached its limit")
        t, x = (tlo + thi) / 2, (xlo + xhi) / 2
        dt, dx = (thi - tlo) / 2, (xhi - xlo) / 2
        z = acb(arb(x), arb(9 + TOP_SLOPE * t))
        value, derivative = evaluate(t, z)
        error, _, _ = variation_bounds(t, z, value, derivative, dt, dx + TOP_SLOPE * dt)
        if error.abs_upper() < abs(value).abs_lower():
            accepted += 1
        elif dt > dx:
            certify_cell(tlo, t, xlo, xhi, depth + 1)
            certify_cell(t, thi, xlo, xhi, depth + 1)
        else:
            certify_cell(tlo, thi, xlo, x, depth + 1)
            certify_cell(tlo, thi, x, xhi, depth + 1)

    certify_cell(fmpq(-1), fmpq(1), LEFT, RIGHT, 0)
    print(f"Top boundary: {accepted} certified cells")


def certify_root_tubes():
    size = fmpq(1, 64)
    dt = fmpq(1, 2 * STEPS)
    accepted = 0
    for direction in (-1, 1):
        centers = [
            acb(arb(fmpq(x, 100)), arb(fmpq(y, 100)))
            for x, y in ((132, 393), (272, 626), (387, 818))
        ]
        for j in range(STEPS):
            t = direction * fmpq(2 * j + 1, 2 * STEPS)
            for k, z in enumerate(centers):
                # The midpoint Newton iterates only propose a Rouché center.
                for _ in range(5):
                    value, derivative = evaluate(t, z)
                    candidate = z - value / derivative
                    z = acb(candidate.real.mid(), candidate.imag.mid())
                value, derivative = evaluate(t, z)
                _, error, derivative_error = variation_bounds(
                    t, z, value, derivative, dt, size
                )
                if not error.abs_upper() < size * abs(derivative).abs_lower():
                    raise ArithmeticError(f"Rouché inequality failed at {t}, root {k}")
                low_height = z.imag - size
                if not (
                    z.real - size > LEFT
                    and z.real + size < RIGHT
                    and low_height > 0
                    and z.imag + size < 9 + TOP_SLOPE * (t - dt)
                ):
                    raise ArithmeticError(f"Root disk leaves contour at {t}, root {k}")
                speed_margin = 4 * low_height * (
                    abs(derivative).abs_lower() - derivative_error
                ) - 1
                if not speed_margin > 5:
                    raise ArithmeticError(f"Speed bound failed at {t}, root {k}")
                centers[k] = z
                accepted += 1
            if not all(
                centers[k + 1].imag - centers[k].imag > 2 * size for k in range(2)
            ):
                raise ArithmeticError(f"Imaginary projections overlap at {t}")
    print(f"Root tubes: {accepted} certified disks, speed margin > 5")


def main():
    ratio = fmpq(FIELD_RADIUS**3, TERMS_H * (TERMS_H - 1))
    ratio += fmpq(FIELD_RADIUS, TERMS_H)
    if ratio >= fmpq(1, 4):
        raise ArithmeticError("Invalid field tail ratio")
    count, segments = winding(LEFT, fmpq(0), RIGHT, fmpq(9))
    if count != 3:
        raise ArithmeticError(f"Expected initial winding 3, obtained {count}")
    print(f"Initial winding: {count}, {segments} certified segments")
    certify_top()
    certify_root_tubes()


if __name__ == "__main__":
    main()
