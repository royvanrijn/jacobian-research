"""Check the explicit coefficients and conductor examples in the family note."""

from collections import defaultdict
from fractions import Fraction
from math import comb, gcd


def multiply(first, second):
    coefficients = defaultdict(Fraction)
    for (i, j), left in first.items():
        for (k, ell), right in second.items():
            coefficients[i + k, j + ell] += left * right
    return coefficients


def two_root_coefficient(left_degree, right_degree, root, degree):
    first = max(0, degree - right_degree)
    final = min(left_degree, degree)
    return sum(
        comb(left_degree, i) * (-1) ** (left_degree - i)
        * comb(right_degree, degree - i) * (-root) ** (right_degree - degree + i)
        for i in range(first, final + 1)
    )


def order_nine_scalars(prime, root):
    k = (prime - 1) // 9
    assert prime == 9 * k + 1 and root % prime not in (0, 1)
    return tuple(
        two_root_coefficient(k * (2 * i % 9), k * (3 * i % 9),
                             root, prime - 1 - k * (5 * i % 9)) % prime
        for i in (1, 2, 4, 5, 7, 8)
    )


def coherent_dimension(order, character):
    return sum(character * a % order for a in (-4, 2, 3, -1)) // order - 1


def check_connected_family(iterate, coefficient):
    prime = 37
    order = (prime ** iterate - 1) // 4
    assert gcd(order, order - 4, 2, 3, order - 1) == 1
    character = order // 3 - 1
    target = prime ** (iterate - 1) * character % order
    assert gcd(character, order) == gcd(target, order) == 1
    assert coherent_dimension(order, character) == 2
    assert coherent_dimension(order, target) == 1
    assert all(coherent_dimension(order, prime ** j) == 1 for j in range(iterate))

    # Lift the weight equality terms, expanding the negative Laurent angles.
    phase = defaultdict(Fraction)
    for i in range(-4, 2):
        value = two_root_coefficient(2, 3, 3, i + 4)
        height = order - 2 * i
        if i >= 0:
            phase[i, height] += value
        else:
            pole_order = -i
            for k in range(pole_order + 1):
                phase[i, height - 2 * pole_order + 2 * k] += (
                    value * (-1) ** (pole_order + k) * comb(pole_order, k)
                )
    # Add 3y+7/2 at weight zero.
    phase[-1, 0] += Fraction(3, 2)
    phase[-1, 2] -= Fraction(3, 2)
    phase[0, 0] += Fraction(7, 2)
    square = multiply(phase, phase)
    fourth_power = multiply(square, square)
    integral = sum(value / (j + 1) for (i, j), value in fourth_power.items()
                   if i == 0 and j % 2 == 0)
    difference = prime ** iterate * integral - coefficient
    assert difference.denominator % prime != 0 and difference.numerator % prime == 0
    print(f"Connected order {order}, degree {order + 4}, normalized moment residue {coefficient % prime}")


def main():
    if not __debug__:
        raise SystemExit("Run with assertions enabled.")
    assert order_nine_scalars(37, 3) == (32, 21, 0, 0, 28, 24)
    assert order_nine_scalars(73, 5) == (11, 0, 22, 11, 0, 44)
    assert order_nine_scalars(73, 6) == (0, 22, 30, 30, 22, 0)
    k = (73 ** 2 - 1) // 9
    assert k % 16 == 0
    assert two_root_coefficient(2 * k, 3 * k, 5, 4 * k) % 73 == 11 ** 2 % 73
    coefficient = two_root_coefficient(8, 12, 3, 16)
    assert coefficient % 37 == 32
    for iterate in (2, 3, 4):
        check_connected_family(iterate, coefficient)
    print("Checked three order-nine Hasse–Witt rows.")


if __name__ == "__main__":
    main()
