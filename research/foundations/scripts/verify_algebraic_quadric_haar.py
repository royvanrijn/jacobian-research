#!/usr/bin/env python3
"""Exact algebraic Haar functional on UV + T^2 = 1.

Work over a characteristic-zero field.  For nonnegative exponents define

    L(U^a V^b T^c) = 0                                      if a != b,
    L((UV)^n T^(2j+1)) = 0,
    L((UV)^n T^(2j))
      = 2*4^n*n!*(2j)!*(n+j+1)! / (j!*(2n+2j+2)!).

No real integration is used to prove invariance.  The checker verifies that
the formula annihilates the quadric ideal and every image of the three
derivations

    D0 = U d/dU - V d/dV,
    D+ = U d/dT - 2T d/dV,
    D- = V d/dT - 2T d/dU.

The written proof explains why these identities and L(1)=1 uniquely force
the formula.  The final comparison with verify_long_xz_mathieu.py is only an
independent regression against the repository's older spherical checker.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb, factorial


Exponent = tuple[int, int, int]  # U, V, T
Polynomial = dict[Exponent, Fraction]
ReducedMonomial = tuple[str, int, int]  # U^r T^c, V^r T^c, or (UV)^r T^c


def add(*polynomials: Polynomial) -> Polynomial:
    """Add sparse polynomials over Q."""
    out: Polynomial = {}
    for polynomial in polynomials:
        for exponent, coefficient in polynomial.items():
            out[exponent] = out.get(exponent, Fraction(0)) + coefficient
            if not out[exponent]:
                del out[exponent]
    return out


def scale(coefficient: Fraction | int, polynomial: Polynomial) -> Polynomial:
    """Scale a sparse polynomial."""
    scalar = Fraction(coefficient)
    return {
        exponent: scalar * value
        for exponent, value in polynomial.items()
        if scalar * value
    }


def multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    """Multiply sparse polynomials."""
    out: Polynomial = {}
    for (au, av, at), left_coefficient in left.items():
        for (bu, bv, bt), right_coefficient in right.items():
            exponent = (au + bu, av + bv, at + bt)
            out[exponent] = (
                out.get(exponent, Fraction(0))
                + left_coefficient * right_coefficient
            )
    return {
        exponent: coefficient
        for exponent, coefficient in out.items()
        if coefficient
    }


def power(polynomial: Polynomial, exponent: int) -> Polynomial:
    """Raise a sparse polynomial to a nonnegative power."""
    assert exponent >= 0
    out = ONE
    base = polynomial
    order = exponent
    while order:
        if order & 1:
            out = multiply(out, base)
        base = multiply(base, base)
        order //= 2
    return out


def apolar_scalar(operator: Polynomial, polynomial: Polynomial) -> Fraction:
    """Evaluate a top constant-coefficient contraction at the origin."""
    return sum(
        operator_coefficient
        * polynomial.get(exponent, Fraction(0))
        * factorial(exponent[0])
        * factorial(exponent[1])
        * factorial(exponent[2])
        for exponent, operator_coefficient in operator.items()
    )


def double_factorial(value: int) -> int:
    result = 1
    for entry in range(value, 0, -2):
        result *= entry
    return result


ONE: Polynomial = {(0, 0, 0): Fraction(1)}
U: Polynomial = {(1, 0, 0): Fraction(1)}
V: Polynomial = {(0, 1, 0): Fraction(1)}
T: Polynomial = {(0, 0, 1): Fraction(1)}
T_SQUARED = multiply(T, T)
QUADRIC = add(multiply(U, V), T_SQUARED, scale(-1, ONE))


def reduce_monomial(a: int, b: int, c: int) -> dict[ReducedMonomial, Fraction]:
    """Reduce U^a V^b T^c to the three declared normal-form families.

    When a=b the balanced factor (UV)^a is retained.  When a>b (respectively
    b>a), expand the paired factor with UV=1-T^2 and obtain only U^r T^e
    (respectively V^r T^e) terms.
    """
    assert min(a, b, c) >= 0
    if a == b:
        return {("UV", a, c): Fraction(1)}

    paired = min(a, b)
    side = "U" if a > b else "V"
    weight = abs(a - b)
    return {
        (side, weight, c + 2 * index): Fraction((-1) ** index * comb(paired, index))
        for index in range(paired + 1)
    }


def balanced_moment(radial_degree: int, height_degree: int) -> Fraction:
    """Closed value of L((UV)^n T^c), derived from Lie recurrences."""
    assert radial_degree >= 0 and height_degree >= 0
    if height_degree % 2:
        return Fraction(0)
    j = height_degree // 2
    return Fraction(
        2
        * 4**radial_degree
        * factorial(radial_degree)
        * factorial(2 * j)
        * factorial(radial_degree + j + 1),
        factorial(j) * factorial(2 * radial_degree + 2 * j + 2),
    )


def reduced_moment(monomial: ReducedMonomial) -> Fraction:
    """Apply L to one reduced normal-form monomial."""
    side, degree, height_degree = monomial
    if side in {"U", "V"}:
        assert degree > 0
        return Fraction(0)
    assert side == "UV"
    return balanced_moment(degree, height_degree)


def haar_monomial(a: int, b: int, c: int) -> Fraction:
    """Apply the algebraic Haar functional to U^a V^b T^c."""
    return sum(
        coefficient * reduced_moment(monomial)
        for monomial, coefficient in reduce_monomial(a, b, c).items()
    )


def haar(polynomial: Polynomial) -> Fraction:
    """Apply the algebraic Haar functional to a sparse polynomial."""
    return sum(
        coefficient * haar_monomial(*exponent)
        for exponent, coefficient in polynomial.items()
    )


def d_zero_monomial(exponent: Exponent) -> Polynomial:
    """Apply D0=U*d_U-V*d_V to one monomial."""
    a, b, _ = exponent
    coefficient = a - b
    return {} if coefficient == 0 else {exponent: Fraction(coefficient)}


def d_plus_monomial(exponent: Exponent) -> Polynomial:
    """Apply D+=U*d_T-2T*d_V to one monomial."""
    a, b, c = exponent
    terms: list[Polynomial] = []
    if c:
        terms.append({(a + 1, b, c - 1): Fraction(c)})
    if b:
        terms.append({(a, b - 1, c + 1): Fraction(-2 * b)})
    return add(*terms)


def d_minus_monomial(exponent: Exponent) -> Polynomial:
    """Apply D-=V*d_T-2T*d_U to one monomial."""
    a, b, c = exponent
    terms: list[Polynomial] = []
    if c:
        terms.append({(a, b + 1, c - 1): Fraction(c)})
    if a:
        terms.append({(a - 1, b, c + 1): Fraction(-2 * a)})
    return add(*terms)


DERIVATIONS = {
    "D0": d_zero_monomial,
    "D+": d_plus_monomial,
    "D-": d_minus_monomial,
}


def derive(polynomial: Polynomial, derivation: str) -> Polynomial:
    """Apply one of the three infinitesimal so3 derivations."""
    operator = DERIVATIONS[derivation]
    return add(
        *(
            scale(coefficient, operator(exponent))
            for exponent, coefficient in polynomial.items()
        )
    )


def long_witness() -> tuple[Polynomial, Polynomial]:
    """Return Long's P=(1+U)(V-(2+U)T^2) and Q=U."""
    p = multiply(
        add(ONE, U),
        add(V, scale(-1, multiply(add(scale(2, ONE), U), T_SQUARED))),
    )
    return p, U


def to_sphere_laurent(polynomial: Polynomial) -> dict[tuple[int, int], Fraction]:
    """Eliminate V=(1-T^2)/U for comparison with the older checker."""
    out: dict[tuple[int, int], Fraction] = {}
    for (a, b, c), coefficient in polynomial.items():
        for index in range(b + 1):
            exponent = (a - b, c + 2 * index)
            out[exponent] = (
                out.get(exponent, Fraction(0))
                + coefficient * (-1) ** index * comb(b, index)
            )
    return {
        exponent: coefficient
        for exponent, coefficient in out.items()
        if coefficient
    }


def main() -> None:
    assert haar(ONE) == 1

    # The derivations preserve UV+T^2-1 as a polynomial identity, before
    # passing to the quotient.
    for derivation in DERIVATIONS:
        assert derive(QUADRIC, derivation) == {}

    # Exact exhaustive monomial audit.  These are regressions of identities
    # proved for arbitrary exponents in the canonical note.
    for a in range(7):
        for b in range(7):
            for c in range(9):
                monomial = {(a, b, c): Fraction(1)}
                assert haar(multiply(QUADRIC, monomial)) == 0
                for derivation in DERIVATIONS:
                    assert haar(derive(monomial, derivation)) == 0

    # The uniqueness recurrences forced by D+, the relation, and L(1)=1.
    # Odd moments start at L((UV)^n T)=0.  Even height moments and the radial
    # moments then satisfy the displayed beta recurrence without integration.
    for n in range(8):
        assert balanced_moment(n, 1) == 0
        for height_degree in range(9):
            assert balanced_moment(n + 1, height_degree) == (
                balanced_moment(n, height_degree)
                - balanced_moment(n, height_degree + 2)
            )
            assert (height_degree + 1) * balanced_moment(
                n + 1, height_degree
            ) == 2 * (n + 1) * balanced_moment(n, height_degree + 2)

    # Reynolds--apolar transfer on every ternary monomial through degree ten.
    # The written proof gives the identity in arbitrary dimension and degree.
    delta_symbol = add(scale(4, multiply(U, V)), T_SQUARED)
    for order in range(6):
        operator = power(delta_symbol, order)
        transfer_constant = (
            2**order
            * factorial(order)
            * double_factorial(2 * order + 1)
        )
        total_degree = 2 * order
        for a in range(total_degree + 1):
            for b in range(total_degree - a + 1):
                c = total_degree - a - b
                monomial = {(a, b, c): Fraction(1)}
                assert apolar_scalar(operator, monomial) == (
                    transfer_constant * haar(monomial)
                )

    # Long's all-order identities are checked in a useful exact range.
    p, q = long_witness()
    p_power = ONE
    for order in range(1, 21):
        p_power = multiply(p_power, p)
        expected = Fraction(
            4**order * factorial(order) ** 2,
            factorial(2 * order + 1),
        )
        assert haar(p_power) == 0
        assert haar(multiply(q, p_power)) == expected

    # Independent agreement with the existing phase/height checker.  This
    # comparison is deliberately last: none of the algebraic invariance
    # assertions above depends on real or compact integration.
    import verify_long_xz_mathieu as spherical_checker

    old_p, old_q = spherical_checker.long_so3_witness()
    assert to_sphere_laurent(p) == old_p
    assert to_sphere_laurent(q) == old_q
    for a in range(5):
        for b in range(5):
            for c in range(7):
                monomial = {(a, b, c): Fraction(1)}
                assert haar(monomial) == spherical_checker.sphere_integral(
                    to_sphere_laurent(monomial)
                )

    p_power = ONE
    old_power = {(0, 0): Fraction(1)}
    for order in range(1, 16):
        p_power = multiply(p_power, p)
        old_power = spherical_checker.sphere_multiply(old_power, old_p)
        assert haar(p_power) == spherical_checker.sphere_integral(old_power)
        assert haar(multiply(q, p_power)) == spherical_checker.sphere_integral(
            spherical_checker.sphere_multiply(old_q, old_power)
        )

    print("PASS quadric Haar: L(1)=1 and the quadric ideal is annihilated")
    print("PASS quadric Haar: D0,D+,D- invariance is exact and algebraic")
    print("PASS quadric Haar: Lie/quadric recurrences force the closed moments")
    print("PASS quadric Haar: Reynolds-apolar transfer through degree 10")
    print("PASS quadric Haar: Long moments through order 20")
    print("PASS quadric Haar: agreement with the spherical checker")


if __name__ == "__main__":
    main()
