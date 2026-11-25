"""Certify the bounded upper zero counts for the zero-temperature cubic."""

from typing import NamedTuple
from flint import acb, arb, ctx, fmpq

ctx.prec = 256
TERMS = 320
DISK_RADIUS = 18


class Sample(NamedTuple):
    x: fmpq
    y: fmpq
    value: acb
    derivative: acb
    midpoint: tuple
    error: arb


def initial_coefficients():
    one_third = arb(fmpq(1, 3))
    two_thirds = arb(fmpq(2, 3))
    coefficients = [
        one_third.gamma() / arb(3) ** two_thirds,
        two_thirds.gamma() / arb(3) ** one_third,
        arb(fmpq(1, 2)),
    ]
    for k in range(3, TERMS + 3):
        coefficients.append(coefficients[k - 3] / (k * (k - 1)))
    if fmpq(DISK_RADIUS**3, TERMS * (TERMS - 1)) >= fmpq(1, 2):
        raise ArithmeticError("Invalid geometric tail ratio")
    tail = 2 * sum(
        coefficients[k] * DISK_RADIUS**k for k in range(TERMS, TERMS + 3)
    )
    derivative_tail = 2 * sum(
        k * coefficients[k] * DISK_RADIUS ** (k - 1)
        for k in range(TERMS, TERMS + 3)
    )
    return coefficients[:TERMS], arb(0, tail.abs_upper()), arb(
        0, derivative_tail.abs_upper()
    )


COEFFICIENTS, TAIL, DERIVATIVE_TAIL = initial_coefficients()


def sample(x, y):
    if x * x + y * y > DISK_RADIUS**2:
        raise ArithmeticError("Point outside the certified disk")
    z = acb(arb(x), arb(y))
    value, derivative = acb(0), acb(0)
    for coefficient in reversed(COEFFICIENTS):
        derivative = derivative * z + value
        value = value * z + coefficient
    value = acb(value.real + TAIL, value.imag + TAIL)
    derivative = acb(
        derivative.real + DERIVATIVE_TAIL, derivative.imag + DERIVATIVE_TAIL
    )
    midpoint = (value.real.mid().fmpq(), value.imag.mid().fmpq())
    point = acb(arb(midpoint[0]), arb(midpoint[1]))
    return Sample(x, y, value, derivative, midpoint, abs(point - value).abs_upper())


def winding(left, bottom, right, top):
    corners = [
        sample(left, bottom),
        sample(right, bottom),
        sample(right, top),
        sample(left, top),
    ]
    corners.append(corners[0])
    polygon = []

    def certify_edge(start, end, depth=0):
        if depth >= 24:
            raise ArithmeticError("Unable to certify a boundary segment")
        length = arb(abs(start.x - end.x) + abs(start.y - end.y))
        # F'' = z F + 1 gives the 1-norm growth bound for (F, F', 1).
        norm = abs(start.value) + abs(start.derivative) + 1
        oscillation = length * norm * (DISK_RADIUS * length).exp()
        disk_radius = oscillation + start.error + end.error
        if disk_radius.abs_upper() < abs(start.value).abs_lower():
            polygon.append(start.midpoint)
            return
        center = sample((start.x + end.x) / 2, (start.y + end.y) / 2)
        certify_edge(start, center, depth + 1)
        certify_edge(center, end, depth + 1)

    for start, end in zip(corners, corners[1:]):
        certify_edge(start, end)
    polygon.append(corners[0].midpoint)

    # Rotate by -i and count positive-real-ray crossings with exact rationals.
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
    rectangles = [
        ("outer", ("0", "0", "16", "8"), 2),
        ("first", ("131/100", "392/100", "133/100", "394/100"), 1),
        ("second", ("271/100", "625/100", "273/100", "627/100"), 1),
    ]
    for name, coordinates, expected in rectangles:
        count, segments = winding(*(fmpq(value) for value in coordinates))
        if count != expected:
            raise ArithmeticError(f"{name}: expected winding {expected}, got {count}")
        print(f"{name}: winding {count}, certified segments {segments}")


if __name__ == "__main__":
    main()
