#!/usr/bin/env python3
"""Exclude every lower-face prolongation of the radial HC4 Schur pair.

For R=x^2+y^2+z^2, start from

    h6 = R^3/30,  h5 = t*R^2.

The degree-fourteen Schur face forces D_t^2 h4=16R.  This checker retains
the resulting completely generic h4, h3, and q2, together with an arbitrary
base quintic.  Two sparse coefficients of the four-variable Hessian
determinant are incompatible.  Hence the exceptional Schur pair does not
extend to a constant-Hessian collision, even with non-invariant lower data.
"""

from __future__ import annotations

import itertools

import sympy as sp


x, y, z, t, lam = sp.symbols("x y z t lambda")
variables = (x, y, z, t)
coefficient_variables = (lam, t, x)
radius = x**2 + y**2 + z**2


def generic_ternary_form(
    degree: int, prefix: str
) -> tuple[sp.Expr, tuple[sp.Symbol, ...]]:
    monomials = [
        x**i * y**j * z ** (degree - i - j)
        for i in range(degree + 1)
        for j in range(degree + 1 - i)
    ]
    coefficients = sp.symbols(f"{prefix}0:{len(monomials)}")
    return (
        sum(
            coefficient * monomial
            for coefficient, monomial in zip(coefficients, monomials)
        ),
        coefficients,
    )


base_quintic, base_quintic_coefficients = generic_ternary_form(5, "r")
mixed_cubic, mixed_cubic_coefficients = generic_ternary_form(3, "u")
base_quartic, base_quartic_coefficients = generic_ternary_form(4, "v")
mixed_quadratic, mixed_quadratic_coefficients = generic_ternary_form(
    2, "w"
)
base_cubic, base_cubic_coefficients = generic_ternary_form(3, "p")

delta = sp.symbols("delta")
ell_x, ell_y, ell_z = sp.symbols("ell_x ell_y ell_z")
kappa, cross_x, cross_y, cross_z = sp.symbols(
    "kappa cross_x cross_y cross_z"
)
base_quadratic_coefficients = sp.symbols("b0:6")
b0, b1, b2, b3, b4, b5 = base_quadratic_coefficients

h6 = radius**3 / 30
h5 = t * radius**2 + base_quintic
h4 = 8 * t**2 * radius + t * mixed_cubic + base_quartic
h3 = (
    delta * t**3
    + t**2 * (ell_x * x + ell_y * y + ell_z * z)
    + t * mixed_quadratic
    + base_cubic
)
q2 = (
    sp.Rational(1, 2) * kappa * t**2
    + t * (cross_x * x + cross_y * y + cross_z * z)
    + sp.Rational(1, 2)
    * (
        b0 * x**2
        + 2 * b1 * x * y
        + 2 * b2 * x * z
        + b3 * y**2
        + 2 * b4 * y * z
        + b5 * z**2
    )
)

pencil = sp.zeros(4)
for weight, homogeneous_part in (
    (4, h6),
    (3, h5),
    (2, h4),
    (1, h3),
    (0, q2),
):
    pencil += lam**weight * sp.hessian(homogeneous_part, variables)

# Only coefficients with y=z=0 are needed.  Differentiate first and then
# restrict to the x-axis, so transverse Hessian jets are retained.
axis_pencil = pencil.applyfunc(
    lambda entry: sp.expand(entry.subs({y: 0, z: 0}))
)
entry_terms = [
    [
        sp.Poly(
            axis_pencil[row, column], *coefficient_variables
        ).terms()
        for column in range(4)
    ]
    for row in range(4)
]


def permutation_sign(permutation: tuple[int, ...]) -> int:
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(4)
        for right in range(left + 1, 4)
    )
    return -1 if inversions % 2 else 1


def determinant_coefficient(
    target: tuple[int, int, int]
) -> sp.Expr:
    total = sp.Integer(0)
    for permutation in itertools.permutations(range(4)):
        partial: dict[tuple[int, ...], sp.Expr] = {
            (0, 0, 0): sp.Integer(permutation_sign(permutation))
        }
        for row, column in enumerate(permutation):
            next_partial: dict[tuple[int, ...], sp.Expr] = {}
            for left_exponents, left_coefficient in partial.items():
                for right_exponents, right_coefficient in entry_terms[
                    row
                ][column]:
                    exponents = tuple(
                        left + right
                        for left, right in zip(
                            left_exponents, right_exponents
                        )
                    )
                    if all(
                        exponent <= bound
                        for exponent, bound in zip(exponents, target)
                    ):
                        next_partial[exponents] = (
                            next_partial.get(exponents, 0)
                            + left_coefficient * right_coefficient
                        )
            partial = next_partial
        total += partial.get(target, 0)
    return sp.factor(total)


face_13 = determinant_coefficient((13, 1, 12))
face_11 = determinant_coefficient((11, 3, 8))

assert sp.expand(
    face_13 - sp.Rational(2, 25) * (3 * delta - 32)
) == 0
assert sp.expand(
    face_11 - sp.Rational(64, 25) * (99 * delta - 1040)
) == 0

all_lower_parameters = set(
    base_quintic_coefficients
    + mixed_cubic_coefficients
    + base_quartic_coefficients
    + mixed_quadratic_coefficients
    + base_cubic_coefficients
    + base_quadratic_coefficients
    + (
        ell_x,
        ell_y,
        ell_z,
        kappa,
        cross_x,
        cross_y,
        cross_z,
    )
)
assert not (face_13.free_symbols & all_lower_parameters)
assert not (face_11.free_symbols & all_lower_parameters)

forced_delta = sp.Rational(32, 3)
residual_11 = sp.factor(face_11.subs(delta, forced_delta))
assert residual_11 == sp.Rational(1024, 25)
assert residual_11 != 0

print("PASS: lambda^13*t*x^12 forces delta=32/3")
print("PASS: lambda^11*t^3*x^8 then equals 1024/25")
print("PASS: every arbitrary lower-form coefficient cancels")
print("SCOPE: the radial exceptional Schur pair cannot prolong to a collision")
