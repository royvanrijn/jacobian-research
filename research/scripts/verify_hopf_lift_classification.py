#!/usr/bin/env python3
"""Exact regressions for the Hopf-lift classification.

The all-order proof is in
extended-geometry/HOPF_LIFT_CLASSIFICATION.md.  This script does not search
V_d: it expands the defining one-variable integral, checks the jet
identities, uses exact SymPy arithmetic for the fixed (4,6) and (5,7)
residual ideal certificates, and constructs the exact (6,8) residual
system.  With ``--require-singular`` it also verifies the rational
modular Groebner/FGLM candidate support profile for that system.  Because
the residual ideal is nonhomogeneous, this optional modular reconstruction
is evidence rather than a deterministic ideal-equality certificate.  The
same option verifies one exact rational boundary exclusion and three exact
specialized ideal memberships, and classifies the exact H=0 slice.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import comb, factorial
import shutil
import subprocess


Polynomial = list[Fraction]
CUTOFF = 20


def trim(polynomial: Polynomial) -> Polynomial:
    result = polynomial[:]
    while result and result[-1] == 0:
        result.pop()
    return result


def multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    result = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return trim(result)


def add(left: Polynomial, right: Polynomial) -> Polynomial:
    result = [Fraction(0)] * max(len(left), len(right))
    for index, coefficient in enumerate(left):
        result[index] += coefficient
    for index, coefficient in enumerate(right):
        result[index] += coefficient
    return trim(result)


def scale(coefficient: Fraction, polynomial: Polynomial) -> Polynomial:
    return trim([coefficient * value for value in polynomial])


def subtract(left: Polynomial, right: Polynomial) -> Polynomial:
    return add(left, scale(Fraction(-1), right))


def divide_with_remainder(
    dividend: Polynomial,
    divisor: Polynomial,
) -> tuple[Polynomial, Polynomial]:
    remainder = trim(dividend)
    divisor = trim(divisor)
    assert divisor
    quotient = [Fraction(0)] * max(1, len(remainder) - len(divisor) + 1)
    while remainder and len(remainder) >= len(divisor):
        degree = len(remainder) - len(divisor)
        leading = remainder[-1] / divisor[-1]
        quotient[degree] += leading
        subtraction = [Fraction(0)] * degree + scale(leading, divisor)
        remainder = subtract(remainder, subtraction)
    return trim(quotient), remainder


def monic_gcd(left: Polynomial, right: Polynomial) -> Polynomial:
    a = trim(left)
    b = trim(right)
    while b:
        _, remainder = divide_with_remainder(a, b)
        a, b = b, remainder
    return scale(Fraction(1, 1) / a[-1], a)


def quadratic_resultant(
    a: Polynomial,
    b: Polynomial,
    c: Polynomial,
    d: Polynomial,
    e: Polynomial,
    f: Polynomial,
) -> Polynomial:
    """Resultant in y of a*y^2+b*y+c and d*y^2+e*y+f."""
    af_minus_cd = subtract(multiply(a, f), multiply(c, d))
    ae_minus_bd = subtract(multiply(a, e), multiply(b, d))
    bf_minus_ce = subtract(multiply(b, f), multiply(c, e))
    return subtract(
        power(af_minus_cd, 2),
        multiply(ae_minus_bd, bf_minus_ce),
    )


def power(polynomial: Polynomial, exponent: int) -> Polynomial:
    result = [Fraction(1)]
    base = polynomial
    value = exponent
    while value:
        if value % 2:
            result = multiply(result, base)
        base = multiply(base, base)
        value //= 2
    return result


def integrated_linear_pair(
    constant_part: Polynomial,
    height_part: Polynomial,
    order: int,
) -> Polynomial:
    """Return integral_0^1 (C(x)+v^2 D(x))^order dv."""
    result: Polynomial = []
    for height_count in range(order + 1):
        term = multiply(
            power(constant_part, order - height_count),
            power(height_part, height_count),
        )
        coefficient = Fraction(
            comb(order, height_count),
            2 * height_count + 1,
        )
        result = add(result, scale(coefficient, term))
    return result


def coefficient(polynomial: Polynomial, degree: int) -> Fraction:
    return polynomial[degree] if degree < len(polynomial) else Fraction(0)


def integrated_profile(profile: Polynomial, order: int) -> Polynomial:
    """Return H(X)=X^(r*m) int_0^1 R(v^2 X^2)^m dv without X^(r*m).

    The returned coefficient at index 2*k is [z^k]R(z)^m/(2*k+1).
    """
    profile_power = power(profile, order)
    result = [Fraction(0)] * (2 * len(profile_power) - 1)
    for degree, coefficient in enumerate(profile_power):
        result[2 * degree] = coefficient / (2 * degree + 1)
    return result


def shifted_coefficient(polynomial: Polynomial, shift: int, degree: int) -> Fraction:
    """Coefficient of u^degree in X^shift P(X), evaluated at X=1+u."""
    return sum(
        coefficient * comb(shift + exponent, degree)
        for exponent, coefficient in enumerate(polynomial)
        if shift + exponent >= degree
    )


def integral_at_one(profile: Polynomial, order: int) -> Fraction:
    return sum(integrated_profile(profile, order))


def odd_double_factorial(order: int) -> int:
    result = 1
    for value in range(1, order + 1, 2):
        result *= value
    return result


def check_profile(r: int, s: int, profile: Polynomial) -> None:
    assert s >= r
    for order in range(1, CUTOFF + 1):
        integrated = integrated_profile(profile, order)
        shift = r * order
        pure_degree = r * order
        adjacent_degree = pure_degree - 1
        expected = integral_at_one(profile, order)
        assert expected
        assert shifted_coefficient(integrated, shift, pure_degree) == 0
        assert (
            shifted_coefficient(integrated, shift, adjacent_degree)
            == expected
        )
        for multiplier_degree in range(1, pure_degree + 1):
            actual = shifted_coefficient(
                integrated,
                shift,
                pure_degree - multiplier_degree,
            )
            ladder = (
                comb(pure_degree - 1, multiplier_degree - 1) * expected
            )
            assert actual == ladder

        # The same Taylor argument classifies the profile exponent q.
        for q in range(max(0, r - 1), r + 2):
            pure = shifted_coefficient(integrated, q * order, pure_degree)
            adjacent = shifted_coefficient(
                integrated, q * order, adjacent_degree
            )
            assert (pure == 0 and adjacent != 0) == (q == r)


def check_high_rectangles(require_singular: bool = False) -> None:
    """Verify the exact residual ideals in the higher fixed rectangles."""
    import sympy as sym

    e, f, g = sym.symbols("E F G")
    constant_part = [sym.S.One, sym.S.One, e, f, g]
    height_part = [
        -sym.S.One,
        -3,
        -3 - 5 * e,
        -1 - 12 * e - 7 * f,
        -3 * (3 * e + 6 * f + 4 * e**2 + 3 * g),
        -2 * e - 15 * f - 24 * e**2 - 24 * g - 36 * e * f,
        (
            -4 * f
            - 15 * e**2
            - 21 * g
            - 78 * e * f
            - 20 * e**3
            - 48 * e * g
            - 27 * f**2
        ),
    ]

    def truncated_multiply(
        left: list[sym.Expr],
        right: list[sym.Expr],
        limit: int,
    ) -> list[sym.Expr]:
        result = [sym.S.Zero] * (min(limit, len(left) + len(right) - 2) + 1)
        for left_degree, left_coefficient in enumerate(left):
            for right_degree, right_coefficient in enumerate(right):
                if left_degree + right_degree <= limit:
                    result[left_degree + right_degree] += (
                        left_coefficient * right_coefficient
                    )
        return [sym.expand(value) for value in result]

    def truncated_power(
        polynomial: list[sym.Expr],
        exponent: int,
        limit: int,
    ) -> list[sym.Expr]:
        result = [sym.S.One]
        base = polynomial
        value = exponent
        while value:
            if value % 2:
                result = truncated_multiply(result, base, limit)
            value //= 2
            if value:
                base = truncated_multiply(base, base, limit)
        return result

    def symbolic_jet(
        order: int,
        current_constant_part: list[sym.Expr],
        current_height_part: list[sym.Expr],
    ) -> sym.Expr:
        result = sym.S.Zero
        for height_count in range(order + 1):
            term = truncated_multiply(
                truncated_power(
                    current_constant_part,
                    order - height_count,
                    order,
                ),
                truncated_power(current_height_part, height_count, order),
                order,
            )
            coefficient_value = sym.Rational(
                comb(order, height_count),
                2 * height_count + 1,
            )
            if order < len(term):
                result += coefficient_value * term[order]
        return sym.factor(result)

    jets = {
        order: symbolic_jet(
            order,
            constant_part,
            height_part,
        )
        for order in range(1, 11)
    }
    assert all(jets[order] == 0 for order in range(1, 7))

    residuals = [
        sym.cancel(jets[7] / sym.Rational(1024, 2145)),
        sym.cancel(jets[8] / sym.Rational(16384, 109395)),
        sym.cancel(jets[9] / sym.Rational(32768, 230945)),
        sym.cancel(jets[10] / sym.Rational(131072, 323323)),
    ]
    target_generators = [
        e**3 - 4 * e**2 - 8 * g,
        g**2,
        f * g,
        -5 * e**2 + 3 * f**2 - 10 * g,
        2 * e**2 + e * g + 4 * g,
        2 * e**2 + 3 * e * f + 4 * g,
    ]

    residual_basis = sym.groebner(residuals, g, f, e, order="grevlex")
    target_basis = sym.groebner(target_generators, g, f, e, order="grevlex")
    assert all(
        sym.expand(residual_basis.reduce(generator)[1]) == 0
        for generator in target_generators
    )
    assert all(
        sym.expand(target_basis.reduce(residual)[1]) == 0
        for residual in residuals
    )

    # The (5,7) rectangle: seven reconstruction jets and five residuals.
    h = sym.symbols("H")
    constant_part_57 = [sym.S.One, sym.S.One, e, f, g, h]
    height_part_57 = [
        -sym.S.One,
        -3,
        -3 - 5 * e,
        -1 - 12 * e - 7 * f,
        -3 * (3 * e + 6 * f + 4 * e**2 + 3 * g),
        -2 * e - 15 * f - 24 * e**2 - 24 * g - 36 * e * f - 11 * h,
        (
            -4 * f
            - 15 * e**2
            - 21 * g
            - 78 * e * f
            - 20 * e**3
            - 48 * e * g
            - 27 * f**2
            - 30 * h
        ),
        -3
        * (
            12 * e**3
            + 32 * e**2 * f
            + e**2
            + 18 * e * f
            + 36 * e * g
            + 20 * e * h
            + 21 * f**2
            + 24 * f * g
            + 2 * g
            + 9 * h
        ),
    ]
    jets_57 = {
        order: symbolic_jet(
            order,
            constant_part_57,
            height_part_57,
        )
        for order in range(1, 13)
    }
    assert all(jets_57[order] == 0 for order in range(1, 8))
    residuals_57 = [
        sym.Poly(jets_57[order], e, f, g, h).primitive()[1].as_expr()
        for order in range(8, 13)
    ]
    residual_basis_57 = sym.groebner(
        residuals_57,
        h,
        g,
        f,
        e,
        order="grevlex",
    )
    assert residual_basis_57.is_zero_dimensional
    lex_basis_57 = residual_basis_57.fglm(order="lex")
    target_generators_57 = [
        (
            323471 * e**5
            + 108000 * e**4
            + 178200 * e**3
            - 307800 * e**2 * f
            + 102600 * e * f**2
            + 583200 * e * f
            - 194400 * f**2
            + 388800 * h
        ),
        16751 * e**5 - 2700 * e**4 - 43200 * e * f**2 + 24300 * g**2,
        (
            -5584 * e**5
            - 7425 * e**4
            - 8100 * e**2 * f
            + 13500 * e * f**2
            + 24300 * f * g
        ),
        (
            -53551 * e**5
            - 15525 * e**4
            + 32400 * e**3
            + 81000 * e**2 * f
            - 27000 * e * f**2
            + 145800 * e * g
            + 97200 * f**2
        ),
        13 * e**5 + 120 * f**3,
        e**2 * (e**3 + 30 * f**2),
        e**3 * (21 * e**2 + 50 * f),
        e**6,
    ]
    target_basis_57 = sym.groebner(
        target_generators_57,
        h,
        g,
        f,
        e,
        order="lex",
    )
    assert all(
        sym.expand(lex_basis_57.reduce(generator)[1]) == 0
        for generator in target_generators_57
    )
    assert all(
        sym.expand(target_basis_57.reduce(residual)[1]) == 0
        for residual in residuals_57
    )

    # The (6,8) rectangle: derive the exact residual system through the
    # predicted cutoff 14.  The optional Singular check reconstructs the
    # candidate support profile but does not certify ideal equality.
    i = sym.symbols("I")
    reconstruction_variables = sym.symbols("B1:9")
    constant_part_68 = [sym.S.One, sym.S.One, e, f, g, h, i]
    unsolved_height_68 = [-sym.S.One, *reconstruction_variables]
    substitutions: dict[sym.Symbol, sym.Expr] = {}
    for order in range(1, 9):
        current_height = [
            sym.expand(sym.sympify(value).subs(substitutions))
            for value in unsolved_height_68
        ]
        jet = symbolic_jet(
            order,
            constant_part_68,
            current_height,
        )
        substitutions[reconstruction_variables[order - 1]] = sym.factor(
            sym.solve(jet, reconstruction_variables[order - 1])[0]
        )
    solved_height_68 = [
        sym.expand(sym.sympify(value).subs(substitutions))
        for value in unsolved_height_68
    ]
    jets_68 = {
        order: symbolic_jet(
            order,
            constant_part_68,
            solved_height_68,
        )
        for order in range(9, 15)
    }
    residuals_68 = [
        sym.Poly(jets_68[order], e, f, g, h, i).primitive()[1].as_expr()
        for order in range(9, 15)
    ]
    assert [
        len(sym.Poly(residual, e, f, g, h, i).terms())
        for residual in residuals_68
    ] == [19, 28, 37, 51, 64, 83]
    residual_jacobian = sym.Matrix(
        [
            [
                sym.diff(residual, variable).subs(
                    {e: 0, f: 0, g: 0, h: 0, i: 0}
                )
                for variable in (e, f, g, h, i)
            ]
            for residual in residuals_68
        ]
    )
    assert residual_jacobian.rank() == 1
    assert residual_jacobian[:, :4] == sym.zeros(6, 4)
    assert residual_jacobian[:, 4] != sym.zeros(6, 1)

    if require_singular:
        singular = shutil.which("Singular")
        if singular is None:
            raise RuntimeError(
                "--require-singular requested, but Singular is not on PATH"
            )
        residual_text = ",".join(
            str(residual).replace("**", "^") for residual in residuals_68
        )
        specialized_residual_text = ",".join(
            str(sym.expand(residual.subs({h: 0, i: 0}))).replace("**", "^")
            for residual in residuals_68
        )
        h_zero_residual_text = ",".join(
            str(sym.expand(residual.subs(h, 0))).replace("**", "^")
            for residual in residuals_68
        )
        singular_program = f"""
LIB "modstd.lib";
ring r=0,(I,H,G,F,E),dp;
option(redSB);
ideal J={residual_text};
ideal Boundary=J,84E+54F+5;
ideal BoundaryBasis=slimgb(Boundary);
int exact_ok=1;
if (reduce(1,BoundaryBasis)!=0) {{ exact_ok=0; }}
ideal K=modStd(J);
int candidate_ok=1;
if (vdim(K)!=32) {{ candidate_ok=0; }}
ring l=0,(I,H,G,F,E),lp;
ideal L=fglm(r,K);
if (size(L)!=18) {{ candidate_ok=0; }}
if (L[1]!=E8) {{ candidate_ok=0; }}
if (L[6]!=F5) {{ candidate_ok=0; }}
poly q=L[13];
q=subst(q,E,0);
q=subst(q,F,0);
if (q!=607500G3) {{ candidate_ok=0; }}
q=L[17];
q=subst(q,E,0);
q=subst(q,F,0);
q=subst(q,G,0);
if (q!=1944000000H2) {{ candidate_ok=0; }}
q=L[18];
q=subst(q,E,0);
q=subst(q,F,0);
q=subst(q,G,0);
q=subst(q,H,0);
if (q!=186624000000I) {{ candidate_ok=0; }}
ring s=0,(G,F,E),dp;
ideal J0={specialized_residual_text};
ideal T=E6,48F3+7E5,972G2+864F2E+29E5-108E4;
matrix M=lift(J0,T);
matrix Q=matrix(T)-matrix(J0)*M;
if (Q[1,1]!=0) {{ exact_ok=0; }}
if (Q[1,2]!=0) {{ exact_ok=0; }}
if (Q[1,3]!=0) {{ exact_ok=0; }}
ring h0=0,(I,G,F,E),dp;
option(redSB);
ideal HJ={h_zero_residual_text};
ideal HK=slimgb(HJ);
if (vdim(HK)!=17) {{ exact_ok=0; }}
ring hlex=0,(I,G,F,E),lp;
ideal HL=fglm(h0,HK);
if (size(HL)!=8) {{ exact_ok=0; }}
if (HL[1]!=E7) {{ exact_ok=0; }}
poly hq=HL[4];
hq=subst(hq,E,0);
if (hq!=165256200000F3) {{ exact_ok=0; }}
hq=HL[7];
hq=subst(hq,E,0);
hq=subst(hq,F,0);
if (hq!=38731921875G2) {{ exact_ok=0; }}
hq=HL[8];
hq=subst(hq,E,0);
hq=subst(hq,F,0);
hq=subst(hq,G,0);
if (hq!=165256200000I) {{ exact_ok=0; }}
"HOPF_68_EXACT_PARTIAL";
exact_ok;
"HOPF_68_CANDIDATE";
candidate_ok;
"""
        completed = subprocess.run(
            [singular, "-q"],
            input=singular_program,
            text=True,
            capture_output=True,
            check=True,
            timeout=120,
        )
        if (
            "HOPF_68_EXACT_PARTIAL\n1" not in completed.stdout
            or "HOPF_68_CANDIDATE\n1" not in completed.stdout
        ):
            raise AssertionError(
                "Singular did not reproduce the (6,8) candidate profile:\n"
                + completed.stdout
                + completed.stderr
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--require-singular",
        action="store_true",
        help="reproduce the rational (6,8) modular Groebner/FGLM profile",
    )
    arguments = parser.parse_args()

    # R(z)=(1-z)^s for several windings and endpoint multiplicities.
    for r in range(1, 5):
        for s in range(r, r + 3):
            profile = power([Fraction(1), Fraction(-1)], s)
            check_profile(r, s, profile)

    # A non-power member: R(z)=(1-z)^2(1+z).
    non_power = multiply(
        power([Fraction(1), Fraction(-1)], 2),
        [Fraction(1), Fraction(1)],
    )
    check_profile(1, 2, non_power)

    # The nonvanishing hypothesis is essential.  For
    # R(z)=(1-z)(z-1/5), the first detecting integral is zero, while every
    # even power has positive integral.
    odd_obstruction = multiply(
        [Fraction(1), Fraction(-1)],
        [Fraction(-1, 5), Fraction(1)],
    )
    assert integral_at_one(odd_obstruction, 1) == 0
    assert integral_at_one(odd_obstruction, 2) > 0

    # In the complete numerator class deg(C)<=1, deg(D)<=3, the first
    # three pure jets are triangular and force D=-C^3.
    samples = [Fraction(-2), Fraction(-1, 2), Fraction(1), Fraction(3, 2)]
    for a in samples:
        for b1 in samples:
            for b2 in samples:
                for b3 in samples:
                    constant_part = [Fraction(1), a]
                    height_part = [Fraction(-1), b1, b2, b3]
                    jet1 = coefficient(
                        integrated_linear_pair(constant_part, height_part, 1),
                        1,
                    )
                    jet2 = coefficient(
                        integrated_linear_pair(constant_part, height_part, 2),
                        2,
                    )
                    jet3 = coefficient(
                        integrated_linear_pair(constant_part, height_part, 3),
                        3,
                    )
                    assert jet1 == Fraction(3 * a + b1, 3)
                    assert jet2 == Fraction(
                        15 * a**2
                        + 10 * a * b1
                        + 3 * b1**2
                        + 4 * b2,
                        15,
                    )
                    assert jet3 == Fraction(
                        35 * a**3
                        + 35 * a**2 * b1
                        + 21 * a * b1**2
                        + 28 * a * b2
                        + 5 * b1**3
                        + 12 * b1 * b2
                        + 8 * b3,
                        35,
                    )

        forced_height = [
            Fraction(-1),
            -3 * a,
            -3 * a**2,
            -a**3,
        ]
        for order in range(1, CUTOFF + 1):
            moment = integrated_linear_pair(
                [Fraction(1), a],
                forced_height,
                order,
            )
            assert coefficient(moment, order) == 0
            assert coefficient(moment, order - 1) != 0

    # In the next complete rectangle deg(C)<=2, deg(D)<=4, four jets solve
    # D, the fifth leaves one false branch, and the sixth removes it.
    quadratic_samples = [Fraction(-3, 2), Fraction(-1, 3), Fraction(0), Fraction(2)]
    for a in samples:
        for e in quadratic_samples:
            constant_part = [Fraction(1), a, e]
            solved_height = [
                Fraction(-1),
                -3 * a,
                -3 * a**2 - 5 * e,
                -a * (a**2 + 12 * e),
                -3 * e * (3 * a**2 + 4 * e),
            ]
            moments = {
                order: integrated_linear_pair(
                    constant_part,
                    solved_height,
                    order,
                )
                for order in range(1, 7)
            }
            assert all(
                coefficient(moments[order], order) == 0
                for order in range(1, 5)
            )
            assert coefficient(moments[5], 5) == Fraction(
                256 * a * e * (a**2 + 12 * e),
                693,
            )
            assert coefficient(moments[6], 6) == Fraction(
                512
                * e
                * (4 * a**4 + 63 * a**2 * e + 20 * e**2),
                3003,
            )

        exceptional_e = -a**2 / 12
        exceptional_height = [
            Fraction(-1),
            -3 * a,
            -3 * a**2 - 5 * exceptional_e,
            -a * (a**2 + 12 * exceptional_e),
            -3
            * exceptional_e
            * (3 * a**2 + 4 * exceptional_e),
        ]
        exceptional_sixth = integrated_linear_pair(
            [Fraction(1), a, exceptional_e],
            exceptional_height,
            6,
        )
        assert coefficient(exceptional_sixth, 6) == Fraction(
            1280 * a**6,
            81081,
        )

    # In the (3,5) rectangle, five jets solve D and the next three reduce
    # to P6=P7=P8=0 in the two weighted parameters E,F.
    cubic_samples = [Fraction(-2, 3), Fraction(0), Fraction(1, 4), Fraction(2)]
    for e_value in cubic_samples:
        for f_value in cubic_samples:
            constant_part = [
                Fraction(1),
                Fraction(1),
                e_value,
                f_value,
            ]
            solved_height = [
                Fraction(-1),
                Fraction(-3),
                -3 - 5 * e_value,
                -1 - 12 * e_value - 7 * f_value,
                -3 * (3 * e_value + 6 * f_value + 4 * e_value**2),
                (
                    -2 * e_value
                    - 15 * f_value
                    - 24 * e_value**2
                    - 36 * e_value * f_value
                ),
            ]
            moments = {
                order: integrated_linear_pair(
                    constant_part,
                    solved_height,
                    order,
                )
                for order in range(1, 9)
            }
            assert all(
                coefficient(moments[order], order) == 0
                for order in range(1, 6)
            )
            p6 = (
                4 * f_value
                + 15 * e_value**2
                + 78 * e_value * f_value
                + 20 * e_value**3
                + 27 * f_value**2
            )
            p7 = (
                f_value
                + 4 * e_value**2
                + 24 * e_value * f_value
                + 8 * e_value**3
                + 12 * f_value**2
                + 8 * e_value**2 * f_value
            )
            p8 = (
                12 * f_value
                + 51 * e_value**2
                + 350 * e_value * f_value
                + 150 * e_value**3
                + 231 * f_value**2
                + 324 * e_value**2 * f_value
                + 24 * e_value**4
                + 90 * e_value * f_value**2
            )
            assert coefficient(moments[6], 6) == Fraction(512, 3003) * p6
            assert coefficient(moments[7], 7) == Fraction(4096, 2145) * p7
            assert coefficient(moments[8], 8) == Fraction(32768, 109395) * p8

    # Treat P6,P7,P8 as quadratics in F with coefficients in Q[E].
    p6_coefficients = (
        [Fraction(27)],
        [Fraction(4), Fraction(78)],
        [Fraction(0), Fraction(0), Fraction(15), Fraction(20)],
    )
    p7_coefficients = (
        [Fraction(12)],
        [Fraction(1), Fraction(24), Fraction(8)],
        [Fraction(0), Fraction(0), Fraction(4), Fraction(8)],
    )
    p8_coefficients = (
        [Fraction(231), Fraction(90)],
        [Fraction(12), Fraction(350), Fraction(324)],
        [
            Fraction(0),
            Fraction(0),
            Fraction(51),
            Fraction(150),
            Fraction(24),
        ],
    )
    resultant_67 = quadratic_resultant(
        *p6_coefficients,
        *p7_coefficients,
    )
    resultant_68 = quadratic_resultant(
        *p6_coefficients,
        *p8_coefficients,
    )
    q5 = [
        Fraction(7),
        Fraction(-156),
        Fraction(-1632),
        Fraction(4928),
        Fraction(-16896),
        Fraction(11520),
    ]
    q6 = [
        Fraction(75),
        Fraction(-1734),
        Fraction(-20117),
        Fraction(32592),
        Fraction(-169440),
        Fraction(64512),
        Fraction(6912),
    ]
    e_squared = [Fraction(0), Fraction(0), Fraction(1)]
    assert resultant_67 == scale(3, multiply(e_squared, q5))
    assert resultant_68 == scale(192, multiply(e_squared, q6))
    assert monic_gcd(q5, q6) == [Fraction(1)]

    # Universal triangular coefficient of b_m in [x^m]K_m.
    for order in range(1, CUTOFF + 1):
        height_part = [Fraction(-1)] + [Fraction(0)] * (order - 1) + [
            Fraction(1)
        ]
        response = integrated_linear_pair(
            [Fraction(1)],
            height_part,
            order,
        )
        expected = Fraction(
            2 ** (order - 1) * factorial(order),
            odd_double_factorial(2 * order + 1),
        )
        assert coefficient(response, order) == expected

    check_high_rectangles(require_singular=arguments.require_singular)

    print(
        "PASS Hopf lifts: full lower-jet ladder and q=r rigidity "
        f"through m={CUTOFF}"
    )
    print("PASS Hopf lifts: polynomial and non-power endpoint profiles")
    print("PASS Hopf lifts: quadratic odd-moment obstruction at m=1")
    print("PASS Hopf lifts: three-jet rigidity in the full (1,3) class")
    print("PASS Hopf lifts: six-jet rigidity in the full (2,4) class")
    print("PASS Hopf lifts: eight-jet rigidity in the full (3,5) class")
    print("PASS Hopf lifts: universal triangular coefficient through m=20")
    print("PASS Hopf lifts: ten-jet rigidity in the full (4,6) class")
    print("PASS Hopf lifts: twelve-jet rigidity in the full (5,7) class")
    if arguments.require_singular:
        print(
            "PASS Hopf lifts: exact partial certificate and modular "
            "candidate support profile in (6,8)"
        )
    else:
        print("PASS Hopf lifts: exact residual system in the (6,8) class")


if __name__ == "__main__":
    main()
