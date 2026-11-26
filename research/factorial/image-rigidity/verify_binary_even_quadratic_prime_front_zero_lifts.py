"""Check the first p-adic lift of the three prime-front zero-rank residues."""

from functools import cache
from math import comb, factorial


@cache
def factorial_cached(n):
    return factorial(n)


def valuation_and_unit(value, prime):
    valuation = 0
    while value % prime == 0:
        value //= prime
        valuation += 1
    return valuation, value % prime


def moment_jet(order, prime, cutoff):
    # Keys are powers of alpha, x, z after kappa is set to -2.
    half_order = (prime - 1) // 2
    coefficients = {}
    minimum_exponent = 2 * order * prime
    for left in range(order + 1):
        choose_left = comb(order, left)
        for cross in range(order - left + 1):
            choose_cross = choose_left * comb(order - left, cross)
            for right in range(order - left - cross + 1):
                constant = order - left - cross - right
                integer_coefficient = (
                    choose_cross
                    * comb(order - left - cross, right)
                    * factorial_cached(2 * left + cross)
                    * factorial_cached(2 * right + cross)
                )
                valuation, unit = valuation_and_unit(integer_coefficient, prime)
                exponent = cross + 2 * (right + constant) + prime * valuation
                minimum_exponent = min(minimum_exponent, exponent)
                if exponent > cutoff:
                    continue
                unit *= pow(2, cross, prime) * pow(
                    -2, constant + half_order * valuation, prime
                )
                for x_power in range(right + 1):
                    z_power = cross + 2 * (right - x_power)
                    index = (exponent, x_power, z_power)
                    contribution = (
                        unit * comb(right, x_power) * pow(-2, x_power, prime)
                    )
                    coefficients[index] = (
                        coefficients.get(index, 0) + contribution
                    ) % prime
    return {key: value for key, value in coefficients.items() if value}, minimum_exponent


def value_and_gradient(coefficients, exponent, x, z, prime):
    value, derivative_x, derivative_z = 0, 0, 0
    for (power, x_power, z_power), coefficient in coefficients.items():
        if power != exponent:
            continue
        value += coefficient * pow(x, x_power, prime) * pow(z, z_power, prime)
        if x_power:
            derivative_x += (
                coefficient * x_power * pow(x, x_power - 1, prime)
                * pow(z, z_power, prime)
            )
        if z_power:
            derivative_z += (
                coefficient * z_power * pow(x, x_power, prime)
                * pow(z, z_power - 1, prime)
            )
    return value % prime, derivative_x % prime, derivative_z % prime


def factorial_polynomial(order, x, prime):
    return sum(
        comb(order, i) * factorial_cached(2 * i) * pow(x, i, prime)
        for i in range(order + 1)
    ) % prime


def verify():
    cases = (
        # p, position parameter, x residue, first order, alpha normalizers,
        # expected derivative and correction rows, and final residual.
        (11, 2, -6, 11, (23, 25, 27), ((1, 2, 2), (9, 1, 4), (3, 4, 1)), 7),
        (19, 1, -2, 18, (18, 39, 41), ((14, 0, 16), (18, 2, 5), (11, 16, 3)), 1),
        (37, 2, -9, 37, (75, 77, 79), ((18, 2, 20), (5, 15, 2), (34, 28, 20)), 36),
    )
    for prime, ell, x, start, exponents, expected_rows, expected_residual in cases:
        h = (prime - 1) // 2
        sign = (-1) ** h
        g = factorial_polynomial(h + 1, x, prime)
        z = -2 * sign * g % prime
        rows = []
        for index, exponent in enumerate(exponents):
            coefficients, raw_minimum = moment_jet(start + index, prime, exponent + 1)
            # Replacing kappa and the p-free integer units preserves both layers.
            assert exponent + 1 < raw_minimum + min(prime, 2 * ell + 3)
            assert min(power for power, _, _ in coefficients) == exponent
            leading, derivative_x, derivative_z = value_and_gradient(
                coefficients, exponent, x, z, prime
            )
            assert leading == 0
            correction, _, _ = value_and_gradient(
                coefficients, exponent + 1, x, z, prime
            )
            rows.append((derivative_x, derivative_z, correction))
        assert tuple(rows) == expected_rows
        (a, b, c), (d, e, f), (g, h, i) = rows
        determinant = (a * e - b * d) % prime
        assert determinant
        inverse = pow(determinant, -1, prime)
        dx = (-c * e + b * f) * inverse % prime
        dz = (-a * f + c * d) * inverse % prime
        residual = (g * dx + h * dz + i) % prime
        assert residual == expected_residual and residual
        print(f"p={prime}, ell={ell}, residues {(x % prime, z)}, "
              f"first corrections {(dx, dz)}, obstructing residual {residual}")


if __name__ == "__main__":
    verify()
