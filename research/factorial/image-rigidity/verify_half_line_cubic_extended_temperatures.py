"""Certify the low-field half-line cubic zeros for -2 <= T <= 2."""

from functools import lru_cache
from math import factorial

from flint import acb, acb_poly, arb, ctx, fmpq

ctx.prec = 512
TERMS_T = 220
TOP_SLOPE = fmpq(43, 40)
LEFT = fmpq(-1)
RIGHT = fmpq(30)
STEPS = 1024


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
            moments[j + 2 * index] / factorial(index)
            for index in range(TERMS_T + 3)
        )
        for j in (0, 1)
    )


INITIAL = initial_series()


@lru_cache(256)
def field_series(temperature, radius, terms):
    t = arb(temperature)
    coefficients = []
    for row in INITIAL:
        value = arb(0)
        for coefficient in reversed(row[:TERMS_T]):
            value = -t / 2 * value + coefficient
        tail = 2 * sum(row[TERMS_T:])
        coefficients.append(value + arb(0, tail.abs_upper()))
    coefficients.append((1 - t * coefficients[1]) / 2)
    for k in range(3, terms):
        coefficients.append(
            (coefficients[k - 3] - t * (k - 1) * coefficients[k - 1])
            / (k * (k - 1))
        )
    # Subsequent blocks of three shrink by at least 1/4 in this weighted norm.
    end = max(
        (abs(coefficients[k]) * radius**k).abs_upper()
        for k in range(terms - 3, terms)
    )
    tail = arb(0, end)
    derivative_tail = arb(0, (2 * (terms + 5) * end / radius).abs_upper())
    polynomial = acb_poly(coefficients)
    return polynomial, polynomial.derivative(), tail, derivative_tail


def evaluate(temperature, z):
    norm = abs(z).abs_upper()
    for radius, terms in ((20, 520), (24, 800), (36, 1050)):
        if norm < radius:
            break
    else:
        raise ArithmeticError("Point outside the coefficient-tail disk")
    polynomial, derivative_polynomial, tail, derivative_tail = field_series(
        temperature, radius, terms
    )
    value, derivative = polynomial(z), derivative_polynomial(z)
    return (
        acb(value.real + tail, value.imag + tail),
        acb(derivative.real + derivative_tail, derivative.imag + derivative_tail),
    )


def variation_bounds(temperature, z, value, derivative, temp_radius, field_radius):
    """Bound value, affine remainder, and derivative on a parameter cell."""
    t, dt, dh = arb(temperature), arb(temp_radius), arb(field_radius)
    p_abs, h_abs = abs(derivative), abs(z)
    norm = abs(value) + p_abs + 1
    temp_norm = 3 * h_abs / 2 + 3
    value_t = (t * derivative - z * value - 1) / 2
    derivative_t = (
        (t * z - 1) * value - (t * t + z) * derivative + t
    ) / 2
    temp_variation = ((temp_norm * dt).exp() - 1) * norm
    value_temp_error = dt * (
        abs(value_t) + (h_abs + 2) * temp_variation / 2 + dt * p_abs / 2
    )
    derivative_temp_error = dt * (
        abs(derivative_t) + temp_norm * temp_variation + dt * (h_abs + 4) * norm / 2
    )

    moved_norm = norm + temp_variation
    field_variation = (((h_abs + dh + 3) * dh).exp() - 1) * moved_norm
    second = (
        abs(z * value + 1 - t * derivative)
        + (h_abs + 2) * (temp_variation + field_variation)
        + dh * (moved_norm + field_variation)
        + dt * p_abs
    )
    for _ in range(2):
        f_move = dh*(p_abs+derivative_temp_error)+dh*dh*second/2
        p_move = dh*second
        second = (abs(z*value+1-t*derivative)
                  + (h_abs+dh)*(value_temp_error+f_move)
                  + 2*(derivative_temp_error+p_move)
                  + dh*abs(value)+dt*p_abs)
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
        if depth > 38:
            raise ArithmeticError("Top subdivision reached its limit")
        t, x = (tlo + thi) / 2, (xlo + xhi) / 2
        dt, dx = (thi - tlo) / 2, (xhi - xlo) / 2
        z = acb(arb(x), arb(15 + TOP_SLOPE * t))
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

    certify_cell(fmpq(-2), fmpq(2), fmpq(3), fmpq(15), 0)
    print(f"Top boundary: {accepted} certified cells")


def certify_root_tubes():
    size = fmpq(1, 128)
    dt = fmpq(1, STEPS)
    accepted = 0
    for direction in (-1, 1):
        centers = [
            acb(arb(fmpq(x, 100)), arb(fmpq(y, 100)))
            for x, y in ((132, 393), (272, 626), (387, 818), (488, 989), (581,1146), (668,1292), (749,1430))
        ]
        for j in range(STEPS):
            t = direction * fmpq(2 * j + 1, STEPS)
            for k, z in enumerate(centers):
                # The midpoint Newton iterates only propose a Rouché center.
                for _ in range(2):
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
                    and z.imag + size < 15 + TOP_SLOPE * (t - dt)
                ):
                    raise ArithmeticError(f"Root disk leaves contour at {t}, root {k}")
                speed_margin = 4 * low_height * (
                    abs(derivative).abs_lower() - derivative_error
                ) - 2
                if not speed_margin > 1:
                    raise ArithmeticError(f"Speed bound failed at {t}, root {k}")
                centers[k] = z
                accepted += 1
            if not all(
                centers[k + 1].imag - centers[k].imag > 2 * size for k in range(len(centers) - 1)
            ):
                raise ArithmeticError(f"Imaginary projections overlap at {t}")
    print(f"Root tubes: {accepted} certified disks, speed margin > 1")


def initial_winding():
    def sample(x, y):
        h = acb(arb(x), arb(y))
        value, derivative = evaluate(fmpq(0), h)
        midpoint = (value.real.mid().fmpq(), value.imag.mid().fmpq())
        mid_value = acb(arb(midpoint[0]), arb(midpoint[1]))
        error = abs(value - mid_value).abs_upper()
        return x, y, value, derivative, midpoint, error

    polygon = []

    def edge(start, end, depth):
        if depth >= 24:
            raise ArithmeticError("Initial winding subdivision reached its limit")
        x0, y0, value, derivative, midpoint, error = start
        x1, y1, _, _, _, end_error = end
        length = arb(abs(x1 - x0) + abs(y1 - y0))
        norm = abs(value) + abs(derivative) + 1
        variation = length * norm * (36 * length).exp()
        if (variation + error + end_error).abs_upper() < abs(value).abs_lower():
            polygon.append(midpoint)
            return
        center = sample((x0 + x1) / 2, (y0 + y1) / 2)
        edge(start, center, depth + 1)
        edge(center, end, depth + 1)

    corners = [
        sample(LEFT, fmpq(0)),
        sample(RIGHT, fmpq(0)),
        sample(RIGHT, fmpq(15)),
        sample(LEFT, fmpq(15)),
    ]
    corners.append(corners[0])
    for start, end in zip(corners, corners[1:]):
        edge(start, end, 0)
    polygon.append(corners[0][4])
    count = 0
    for (a, b), (c, d) in zip(polygon, polygon[1:]):
        x0, y0, x1, y1 = b, -a, d, -c
        determinant = x0 * y1 - x1 * y0
        if y0 <= 0 < y1 and determinant > 0:
            count += 1
        if y1 <= 0 < y0 and determinant < 0:
            count -= 1
    return count, len(polygon) - 1


def main():
    for radius, terms in ((20, 520), (24, 800), (36, 1050)):
        ratio = fmpq(radius**3, terms * (terms - 1)) + fmpq(2 * radius, terms)
        if ratio >= fmpq(1, 4):
            raise ArithmeticError("Invalid field tail ratio")
    count, segments = initial_winding()
    if count != 7:
        raise ArithmeticError(f"Expected initial winding 7, obtained {count}")
    print(f"Initial winding: {count}, {segments} certified segments", flush=True)
    certify_top()
    field_series.cache_clear()
    certify_root_tubes()


if __name__ == "__main__":
    main()
