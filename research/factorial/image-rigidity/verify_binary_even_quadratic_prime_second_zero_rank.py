"""Exact rational arithmetic in the uniform prime-second zero-rank lift."""

from fractions import Fraction
from math import comb, factorial


def verify():
    x = Fraction(359, 300)
    g = Fraction(4, 5)
    h = 0  # Value of g_h at the forced residue.
    h_derivative = (g - 1) / (2 * x**2)
    g_derivative = (g - h) / (2 * x)
    g_second = ((g_derivative - h_derivative) / x - (g - h) / x**2) / 2
    g1 = 1 + 6 * x * g - 3 * x * h
    g1_derivative = 3 * (g1 - g) / (2 * x)
    g1_second = 3 * ((g1_derivative - g_derivative) / x - (g1 - g) / x**2) / 2
    g2 = 1 + 5 * x * (4 * g1 - 3 * g)
    v = 1 + 2 * x - 4 * g**2
    assert 3 * g * v == 2

    for sign in (1, -1):
        z = -2 * sign * g
        dx = sign * (4 * g * x - 2 * x**2 / (g - 1))
        dz = -v - sign * g * dx / x
        w = x - 2 * g**2
        theta = -z * (1 + 4 * w + 24 * w**2) / 3 - Fraction(16, 45) * sign
        assert theta == Fraction(4, 45) * sign

        q0 = g + 2 * x * g_derivative
        q1 = g1 + 2 * x * g1_derivative
        f0_x = 2 * sign * z * (3 * g_derivative + 2 * x * g_second) + 2
        f0_z = 2 * z + 2 * sign * q0
        f1_x = -12 * z**2 - Fraction(4, 3) * sign * z * (
            3 * g1_derivative + 2 * x * g1_second
        ) - 2 - 24 * x
        f1_z = 4 * z**3 - 2 * z * (1 + 12 * x) - Fraction(4, 3) * sign * q1

        k0 = (
            z**3 + z * (1 + 6 * x)
            + 2 * sign * z**2 * (g + 5 * x * g_derivative + 2 * x**2 * g_second)
            + Fraction(2, 3) * sign * g1
        )
        k1 = (
            z**5 - z**3 * (1 + 20 * x) - z * (1 + 12 * x + 120 * x**2) / 2
            - Fraction(4, 3) * sign * z**2
            * (g1 + 5 * x * g1_derivative + 2 * x**2 * g1_second)
            - Fraction(4, 15) * sign * g2 + theta / 2
        )
        second0 = k0 + f0_x * dx + f0_z * dz + sign * g_second * dx**2
        second1 = (
            k1 + f1_x * dx + f1_z * dz
            - Fraction(2, 3) * sign * g1_second * dx**2 - 6 * dx * dz + 3 * z * dz**2
        )
        multiplier = 2 - 2 / g
        residual = second1 - multiplier * second0
        assert residual == -Fraction(553, 72) * sign
        print(f"C={sign}: exact second-digit residual {residual}")

    for prime, expected in ((7, 5), (79, 45)):
        x_residue = x.numerator * pow(x.denominator, -1, prime) % prime
        half_order = (prime - 1) // 2
        value = sum(
            comb(half_order, i) * factorial(2 * i) * pow(x_residue, i, prime)
            for i in range(half_order + 1)
        ) % prime
        assert value == expected and value
        print(f"p={prime}: g_h({x_residue}) = {value} mod {prime}")


if __name__ == "__main__":
    verify()
