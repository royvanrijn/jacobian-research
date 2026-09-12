#!/usr/bin/env python3
"""Exact algebraic Haar functional on k[SL2].

Coordinates are ordered as

        [a c]
    g = [b d],

with ad-bc=1.  The normalized functional is

    H(a^r b^s c^t d^u)
      = (-1)^s delta_(r,u) delta_(s,t) r!s!/(r+s+1)!.

The checker proves algebraically, in exact arithmetic, that the closed
formula annihilates (ad-bc-1)f and the images of the left and right sl2
derivations.  Its final comparisons with older compact-coordinate scripts
are independent regressions, not inputs to the invariance proof.
"""

from __future__ import annotations

from fractions import Fraction
from math import factorial


Exponent = tuple[int, int, int, int]  # a, b, c, d
Polynomial = dict[Exponent, Fraction]
VectorField = tuple[tuple[int, int, int], ...]  # scalar, differentiated, output


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
    for left_exponent, left_coefficient in left.items():
        for right_exponent, right_coefficient in right.items():
            exponent = tuple(
                left_exponent[index] + right_exponent[index]
                for index in range(4)
            )
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


ONE: Polynomial = {(0, 0, 0, 0): Fraction(1)}
A: Polynomial = {(1, 0, 0, 0): Fraction(1)}
B: Polynomial = {(0, 1, 0, 0): Fraction(1)}
C: Polynomial = {(0, 0, 1, 0): Fraction(1)}
D: Polynomial = {(0, 0, 0, 1): Fraction(1)}
DETERMINANT_RELATION = add(multiply(A, D), scale(-1, multiply(B, C)), scale(-1, ONE))


def haar_monomial(r: int, s: int, t: int, u: int) -> Fraction:
    """Apply the algebraic SL2 Haar functional to one monomial."""
    assert min(r, s, t, u) >= 0
    if r != u or s != t:
        return Fraction(0)
    return Fraction(
        (-1) ** s * factorial(r) * factorial(s),
        factorial(r + s + 1),
    )


def haar(polynomial: Polynomial) -> Fraction:
    """Apply the algebraic SL2 Haar functional to a sparse polynomial."""
    return sum(
        coefficient * haar_monomial(*exponent)
        for exponent, coefficient in polynomial.items()
    )


# Infinitesimal left multiplication g -> exp(tX)g.
LEFT_H: VectorField = (
    (1, 0, 0),
    (-1, 1, 1),
    (1, 2, 2),
    (-1, 3, 3),
)
LEFT_E: VectorField = ((1, 0, 1), (1, 2, 3))  # b*d_a+d*d_c
LEFT_F: VectorField = ((1, 1, 0), (1, 3, 2))  # a*d_b+c*d_d

# Infinitesimal right multiplication g -> g exp(tX).
RIGHT_H: VectorField = (
    (1, 0, 0),
    (1, 1, 1),
    (-1, 2, 2),
    (-1, 3, 3),
)
RIGHT_E: VectorField = ((1, 2, 0), (1, 3, 1))  # a*d_c+b*d_d
RIGHT_F: VectorField = ((1, 0, 2), (1, 1, 3))  # c*d_a+d*d_b

DERIVATIONS = {
    "L_H": LEFT_H,
    "L_E": LEFT_E,
    "L_F": LEFT_F,
    "R_H": RIGHT_H,
    "R_E": RIGHT_E,
    "R_F": RIGHT_F,
}


def derive_monomial(exponent: Exponent, vector_field: VectorField) -> Polynomial:
    """Apply a linear vector field to one monomial."""
    out: Polynomial = {}
    for scalar, differentiated, output_variable in vector_field:
        degree = exponent[differentiated]
        if not degree:
            continue
        result = list(exponent)
        result[differentiated] -= 1
        result[output_variable] += 1
        result_exponent = tuple(result)
        out[result_exponent] = (
            out.get(result_exponent, Fraction(0))
            + Fraction(scalar * degree)
        )
    return {
        result_exponent: coefficient
        for result_exponent, coefficient in out.items()
        if coefficient
    }


def derive(polynomial: Polynomial, derivation: str) -> Polynomial:
    """Apply one of the six left/right sl2 derivations."""
    vector_field = DERIVATIONS[derivation]
    return add(
        *(
            scale(coefficient, derive_monomial(exponent, vector_field))
            for exponent, coefficient in polynomial.items()
        )
    )


def long_witness() -> tuple[Polynomial, Polynomial]:
    """Return F=(1+c)(ad+b) and G=-c."""
    f = multiply(add(ONE, C), add(multiply(A, D), B))
    return f, scale(-1, C)


def substitute_quadric_monomial(
    u_degree: int, v_degree: int, t_degree: int
) -> Polynomial:
    """Pull a quadric monomial back along SL2/T -> {UV+T^2=1}.

    For g=[[a,c],[b,d]], conjugating diag(1,-1) gives

        U=-2ac, V=2bd, T=ad+bc,

    and U*V+T^2=(ad-bc)^2.
    """
    quadric_u = scale(-2, multiply(A, C))
    quadric_v = scale(2, multiply(B, D))
    quadric_t = add(multiply(A, D), multiply(B, C))
    return multiply(
        multiply(power(quadric_u, u_degree), power(quadric_v, v_degree)),
        power(quadric_t, t_degree),
    )


def main() -> None:
    assert haar(ONE) == 1

    # All six derivations preserve ad-bc-1 before quotienting.
    for derivation in DERIVATIONS:
        assert derive(DETERMINANT_RELATION, derivation) == {}

    # Exact exhaustive monomial audit of the Hopf ideal and left/right
    # infinitesimal invariance.  The canonical note proves the identities for
    # arbitrary exponents.
    for r in range(6):
        for s in range(6):
            for t in range(6):
                for u in range(6):
                    monomial = {(r, s, t, u): Fraction(1)}
                    assert haar(multiply(DETERMINANT_RELATION, monomial)) == 0
                    for derivation in DERIVATIONS:
                        assert haar(derive(monomial, derivation)) == 0

    # The Cartan equations force r=u and s=t.  On balanced monomials h[r,s],
    # a raising equation and ad-bc=1 give the two recurrences below, starting
    # from h[0,0]=1.  This derives rather than assumes the factorial formula.
    for r in range(8):
        for s in range(8):
            current = haar_monomial(r, s, s, r)
            next_r = haar_monomial(r + 1, s, s, r + 1)
            next_s = haar_monomial(r, s + 1, s + 1, r)
            assert (r + 1) * next_s + (s + 1) * next_r == 0
            assert next_r - next_s == current

    # Long's witness, now entirely inside k[SL2].
    f, g = long_witness()
    f_power = ONE
    for order in range(1, 21):
        f_power = multiply(f_power, f)
        assert haar(f_power) == 0
        assert haar(multiply(g, f_power)) == Fraction(
            (-1) ** (order - 1), order + 1
        )

    # The explicit SL2/T quadric quotient and Haar pullback identity.
    import verify_algebraic_quadric_haar as quadric

    pulled_relation = add(
        multiply(
            substitute_quadric_monomial(1, 0, 0),
            substitute_quadric_monomial(0, 1, 0),
        ),
        power(substitute_quadric_monomial(0, 0, 1), 2),
        scale(-1, ONE),
    )
    assert pulled_relation == multiply(
        DETERMINANT_RELATION,
        add(DETERMINANT_RELATION, scale(2, ONE)),
    )
    for u_degree in range(5):
        for v_degree in range(5):
            for t_degree in range(7):
                assert haar(
                    substitute_quadric_monomial(u_degree, v_degree, t_degree)
                ) == quadric.haar_monomial(u_degree, v_degree, t_degree)

    # Independent agreement with the repository's older beta checker.
    import verify_long_xz_mathieu as compact_checker

    old_power = compact_checker.ONE
    f_power = ONE
    for order in range(1, 16):
        old_power = compact_checker.multiply(old_power, compact_checker.f)
        f_power = multiply(f_power, f)
        assert haar(f_power) == compact_checker.integral_constant_term(old_power)
        assert haar(multiply(g, f_power)) == compact_checker.integral_constant_term(
            compact_checker.multiply(compact_checker.Z1_INV, old_power)
        )

    print("PASS SL2 Haar: H(1)=1 and the determinant ideal is annihilated")
    print("PASS SL2 Haar: all left/right sl2 derivations are annihilated")
    print("PASS SL2 Haar: Lie/determinant recurrences force the factorial formula")
    print("PASS SL2 Haar: Long moments through order 20")
    print("PASS SL2/T: explicit quadric relation and Haar pullback")
    print("PASS SL2 Haar: agreement with the compact beta checker")


if __name__ == "__main__":
    main()
