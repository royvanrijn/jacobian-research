"""Replay the formal coefficient identities for the tangent-conic (12,6) exclusion."""

from sympy import Poly, Rational, binomial, expand, symbols


Z, a, b, c, f, e, u, rho, t = symbols("Z a b c f e u rho t")
S = Z**2 + a
U = (b * Z + c) * S + f * Z + e
K = c + 3 * u / 2


def negative_coefficient(exponent, order):
    """Coefficient of Z**(-order) in (S**3 + U)**exponent."""
    result = 0
    leading_degree = int(6 * exponent)
    for index in range((leading_degree + order) // 3 + 1):
        # U has degree three, so indices above this bound cannot contribute.
        for (degree,), coefficient in Poly(U**index, Z).terms():
            difference = leading_degree - 6 * index + degree + order
            if difference < 0 or difference % 2:
                continue
            shift = difference // 2
            result += (
                binomial(exponent, index)
                * coefficient
                * binomial(3 * (exponent - index), shift)
                * a**shift
            )
    return expand(result)


def coefficient(order):
    return negative_coefficient(Rational(4, 3), order) + u * negative_coefficient(
        Rational(2, 3), order
    )


first = b * e + K * f - b**3 / 9
third = -a * b**3 + 3 * b**2 * f + 3 * b * c * K - 9 * f * e
fifth = b * (b**2 * (4 * c + 3 * u) - 18 * f**2)

assert expand(coefficient(1) - 4 * first / 9) == 0
assert expand(coefficient(3) + 2 * (18 * a * first + 2 * third) / 81) == 0
assert expand(
    243 * (coefficient(5) + b * coefficient(2) / 3)
    - (108 * a**2 - 36 * c) * first
    - 24 * a * third
    - fifth
) == 0

eliminate_e = expand(3 * b**2 * f + 3 * b * c * K - 9 * f * (b**2 / 9 - K * f / b))
f_squared = (b**2 * (4 * c + 3 * u) - rho / (b * t**3)) / 18
eliminate_f_squared = 2 * b**2 * f + b * K * (5 * c + 3 * u / 2) - rho * K / (
    2 * b**2 * t**3
)
assert expand(eliminate_e.subs(f**2, f_squared) - eliminate_f_squared) == 0

print("Sextic tangent conic: centered (8,6) formal coefficients and elimination verified.")
