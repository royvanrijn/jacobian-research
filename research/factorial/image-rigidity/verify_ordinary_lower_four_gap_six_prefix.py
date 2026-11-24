"""Exact finite-field certificates for all lower-four gap-six terminals below 10000."""

from math import isqrt

from flint import nmod_poly


PRIME = 1_000_003
FIRST_TERMINAL = 2
LAST_TERMINAL = 9_999


def verify() -> None:
    factorial_limit = 10 * LAST_TERMINAL
    if factorial_limit >= PRIME:
        raise ArithmeticError("The factorial table would contain a nonunit.")
    if any(PRIME % divisor == 0 for divisor in range(2, isqrt(PRIME) + 1)):
        raise ArithmeticError("The specified modulus must be prime.")

    factorial = [1] * (factorial_limit + 1)
    inverse_factorial = [1] * (factorial_limit + 1)
    for index in range(1, factorial_limit + 1):
        factorial[index] = factorial[index - 1] * index % PRIME
    inverse_factorial[-1] = pow(factorial[-1], -1, PRIME)
    for index in range(factorial_limit, 0, -1):
        inverse_factorial[index - 1] = inverse_factorial[index] * index % PRIME

    def reciprocal(order: int) -> nmod_poly:
        common = factorial[order] * inverse_factorial[4 * order] % PRIME
        coefficients = [
            common * inverse_factorial[index] * inverse_factorial[order - index]
            * factorial[4 * order + 6 * index] % PRIME
            for index in range(order, -1, -1)
        ]
        result = nmod_poly(coefficients, PRIME)
        if result.degree() != order or result[order] != 1:
            raise ArithmeticError(f"Incorrect reciprocal polynomial at order {order}.")
        return result

    lower = reciprocal(FIRST_TERMINAL - 1)
    for terminal in range(FIRST_TERMINAL, LAST_TERMINAL + 1):
        upper = reciprocal(terminal)
        if lower.gcd(upper).degree() != 0:
            raise ArithmeticError(f"Nonconstant modular gcd at terminal {terminal}.")
        lower = upper

    checked = LAST_TERMINAL - FIRST_TERMINAL + 1
    print(f"prime {PRIME}")
    print(f"terminals {FIRST_TERMINAL} through {LAST_TERMINAL}")
    print(f"coprime pairs {checked}")


if __name__ == "__main__":
    verify()

