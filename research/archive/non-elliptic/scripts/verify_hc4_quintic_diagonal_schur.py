#!/usr/bin/env python3
"""Verify the diagonal rank-three sextic--quintic Schur obstruction.

The calculation keeps every base-only quintic, every remaining quartic
term, every cubic term, and an arbitrary quadratic form.  A truncated
sparse determinant extracts only the four coefficients used in the proof,
so all unused parameters remain genuinely generic without causing a dense
symbolic expansion.
"""

from __future__ import annotations

import itertools
import re
import shutil
import subprocess

import sympy as sp


x, y, z, t, lam = sp.symbols("x y z t lam")
a, b, c, delta = sp.symbols("a b c delta")
variables = (x, y, z, t)
determinant_variables = (lam, t, x, y, z)


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
base_cubic, base_cubic_coefficients = generic_ternary_form(3, "f")

ell_x, ell_y, ell_z = sp.symbols("ell_x ell_y ell_z")
kappa, cross_x, cross_y, cross_z = sp.symbols(
    "kappa cross_x cross_y cross_z"
)
base_quadratic_coefficients = sp.symbols("base0:6")
base0, base1, base2, base3, base4, base5 = (
    base_quadratic_coefficients
)

h6 = (x**6 + y**6 + z**6) / 30
h5 = t * (a * x**4 + b * y**4 + c * z**4) + base_quintic

# The lambda^14 Schur face forces precisely this t^2 coefficient.
h4 = (
    8 * t**2 * (a**2 * x**2 + b**2 * y**2 + c**2 * z**2)
    + t * mixed_cubic
    + base_quartic
)
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
        base0 * x**2
        + 2 * base1 * x * y
        + 2 * base2 * x * z
        + base3 * y**2
        + 2 * base4 * y * z
        + base5 * z**2
    )
)


# First verify the polynomial Schur norm itself.
sextic_block = sp.hessian(h6, (x, y, z))
quintic_cross = sp.Matrix(
    [
        sp.diff(a * x**4 + b * y**4 + c * z**4, variable)
        for variable in (x, y, z)
    ]
)
schur_numerator = sp.expand(
    (
        quintic_cross.T
        * sextic_block.adjugate()
        * quintic_cross
    )[0]
)
schur_quotient = 16 * (
    a**2 * x**2 + b**2 * y**2 + c**2 * z**2
)
assert sp.expand(
    schur_numerator - sextic_block.det() * schur_quotient
) == 0


# Before using the diagonal form, certify that it is forced.  For a generic
# quartic s, divisibility of
#
#   y^4*z^4*s_x^2 + x^4*z^4*s_y^2 + x^4*y^4*s_z^2
#
# by x^4*y^4*z^4 has radical equal to the ideal of all twelve mixed
# quartic coefficients.  The elementary proof is the valuation argument
# x^4 | s_x^2 => x^2 | s_x, and cyclically; the ideal calculation is an
# independent exact coefficient certificate.
generic_quartic, quartic_coefficients = generic_ternary_form(4, "q")
generic_gradient = sp.Matrix(
    [
        sp.diff(generic_quartic, variable)
        for variable in (x, y, z)
    ]
)
generic_schur_numerator = sp.expand(
    (
        generic_gradient.T
        * sextic_block.adjugate()
        * generic_gradient
    )[0]
)
bad_divisibility_equations = [
    coefficient
    for monomial, coefficient in sp.Poly(
        generic_schur_numerator, x, y, z
    ).terms()
    if not (
        monomial[0] >= 4
        and monomial[1] >= 4
        and monomial[2] >= 4
    )
]
pure_quartic_coefficients = {
    quartic_coefficients[0],
    quartic_coefficients[4],
    quartic_coefficients[14],
}
mixed_quartic_coefficients = tuple(
    coefficient
    for coefficient in quartic_coefficients
    if coefficient not in pure_quartic_coefficients
)
assert len(bad_divisibility_equations) == 66
assert len(mixed_quartic_coefficients) == 12
assert all(
    sp.expand(
        equation.subs(
            {
                coefficient: 0
                for coefficient in mixed_quartic_coefficients
            }
        )
    )
    == 0
    for equation in bad_divisibility_equations
)

singular = shutil.which("Singular")
if singular is None:
    raise RuntimeError("Singular is required for the exact radical check")


def singular_expression(expression: sp.Expr) -> str:
    return str(sp.expand(expression)).replace("**", "^")


quartic_program = f"""
ring rr=0,({",".join(map(str, quartic_coefficients))}),dp;
option(redSB);
ideal I={",".join(
    map(singular_expression, bad_divisibility_equations)
)};
ideal GI=slimgb(I);
ideal R={",".join(map(str, mixed_quartic_coefficients))};
ideal GR=std(R);
print("QUARTIC_BASE "+string(size(GI))+" "+string(size(reduce(I,GR))));
int i;
int exponent;
poly power;
for (i=1;i<=size(R);i++)
{{
  power=R[i];
  exponent=1;
  while (exponent<=4 && reduce(power,GI)!=0)
  {{
    power=power*R[i];
    exponent++;
  }}
  print(
    "QUARTIC_POWER "+string(i)+" "+string(exponent)+" "
    +string(reduce(power,GI)==0)
  );
}}
"""
quartic_output = subprocess.run(
    [singular, "-q"],
    input=quartic_program,
    text=True,
    capture_output=True,
    check=True,
    timeout=120,
).stdout
quartic_base = re.search(
    r"(?m)^QUARTIC_BASE (\d+) (\d+)$", quartic_output
)
assert quartic_base is not None
assert tuple(map(int, quartic_base.groups())) == (101, 0)
quartic_powers = re.findall(
    r"(?m)^QUARTIC_POWER (\d+) (\d+) ([01])$",
    quartic_output,
)
assert len(quartic_powers) == 12
assert all(success == "1" for _, _, success in quartic_powers)
assert max(int(exponent) for _, exponent, _ in quartic_powers) == 3


pencil = sp.zeros(4)
for weight, homogeneous_part in (
    (4, h6),
    (3, h5),
    (2, h4),
    (1, h3),
    (0, q2),
):
    pencil += lam**weight * sp.hessian(
        homogeneous_part, variables
    )


# Store every entry sparsely in (lambda,t,x,y,z).  The coefficient extractor
# multiplies term dictionaries while discarding exponents above the requested
# target.  This is an exact determinant expansion, not a numerical sample.
entry_terms = [
    [
        sp.Poly(
            sp.expand(pencil[row, column]),
            *determinant_variables,
        ).terms()
        for column in range(4)
    ]
    for row in range(4)
]


def permutation_sign(permutation: tuple[int, ...]) -> int:
    inversions = sum(
        permutation[i] > permutation[j]
        for i in range(4)
        for j in range(i + 1, 4)
    )
    return -1 if inversions % 2 else 1


def determinant_coefficient(
    target: tuple[int, int, int, int, int]
) -> sp.Expr:
    total = sp.Integer(0)
    for permutation in itertools.permutations(range(4)):
        partial: dict[tuple[int, ...], sp.Expr] = {
            (0, 0, 0, 0, 0): sp.Integer(
                permutation_sign(permutation)
            )
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


face_13 = determinant_coefficient((13, 1, 4, 4, 4))
face_11_x = determinant_coefficient((11, 3, 0, 4, 4))
face_11_y = determinant_coefficient((11, 3, 4, 0, 4))
face_11_z = determinant_coefficient((11, 3, 4, 4, 0))

assert sp.expand(
    face_13
    + 2 * (32 * a**3 + 32 * b**3 + 32 * c**3 - 3 * delta)
) == 0
assert sp.expand(
    face_11_x
    + 32 * a**2 * (32 * b**3 + 32 * c**3 - 3 * delta)
) == 0
assert sp.expand(
    face_11_y
    + 32 * b**2 * (32 * a**3 + 32 * c**3 - 3 * delta)
) == 0
assert sp.expand(
    face_11_z
    + 32 * c**2 * (32 * a**3 + 32 * b**3 - 3 * delta)
) == 0

# Eliminate delta by the lambda^13 face.  The remaining coefficients are
# fifth powers, so characteristic zero forces a=b=c=0.
forced_delta = sp.Rational(32, 3) * (a**3 + b**3 + c**3)
assert sp.expand(face_11_x.subs(delta, forced_delta) - 1024 * a**5) == 0
assert sp.expand(face_11_y.subs(delta, forced_delta) - 1024 * b**5) == 0
assert sp.expand(face_11_z.subs(delta, forced_delta) - 1024 * c**5) == 0

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
for coefficient in (face_13, face_11_x, face_11_y, face_11_z):
    assert not (coefficient.free_symbols & all_lower_parameters)

print("PASS: Schur divisibility forces every quartic to be diagonal")
print("PASS: the diagonal quartic has a polynomial sextic Schur norm")
print("PASS: lambda^13 forces delta=32*(a^3+b^3+c^3)/3")
print("PASS: lambda^11 leaves 1024*a^5, 1024*b^5, 1024*c^5")
print("PASS: all base-quintic, quartic, cubic, and quadratic data cancel")
print("SCOPE: the full Fermat-sextic quintic Schur-norm stratum is closed")
