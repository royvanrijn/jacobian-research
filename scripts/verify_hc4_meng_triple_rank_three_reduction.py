#!/usr/bin/env python3
"""Verify the exact algebra in the rank-three triple-layer HC(4) reduction.

The potential is

    psi = q_2 + h_3 + h_4 + h_6.

If Hess(h_6) has generic rank three, its constant kernel direction t and
the two highest determinant layers force

    D_t h_6 = D_t^2 h_4 = D_t^2 h_3 = 0.

The nonisotropic t-direction Schur-descends to HC(3).  In the isotropic
case psi=t*s(u)+phi(u), where s has degree at most three.  The coefficient
of t^2 is the bordered invariant

    -grad(s)^T adj(Hess(s)) grad(s).

This checker verifies that every cubic s with vanishing bordered invariant
has a nonzero constant direction along which s is independent.  After
that direction is selected, the shared cotangent and terminal HC(2)
identities complete the reduction.  The conceptual inputs not checked here
are Gordan--Noether, the elementary binary-cubic orbit classification,
the two-variable singular-Hessian classification, HC(3), HC(2), and Moh's
plane degree bound.
"""

from __future__ import annotations

import contextlib
import io
from pathlib import Path
import runpy

import sympy as sp


# Recheck the cotangent determinant and terminal HC(2) block used after the
# constant direction in s is found.
SHARED_CHECKER = Path(__file__).with_name(
    "verify_hc4_meng_dense_rank_three_sextic_reduction.py"
)
with contextlib.redirect_stdout(io.StringIO()):
    runpy.run_path(str(SHARED_CHECKER))


# 1. The two top spatial layers when Hess(h_6) has rank three.
lam = sp.symbols("lam")
c1, c2, c3 = sp.symbols("c1 c2 c3")
C6 = sp.diag(c1, c2, c3, 0)


def generic_symmetric(prefix: str) -> sp.Matrix:
    entries = sp.symbols(f"{prefix}0:10")
    return sp.Matrix(
        [
            [entries[0], entries[1], entries[2], entries[3]],
            [entries[1], entries[4], entries[5], entries[6]],
            [entries[2], entries[5], entries[7], entries[8]],
            [entries[3], entries[6], entries[8], entries[9]],
        ]
    )


H0 = generic_symmetric("h")
A3 = generic_symmetric("a")
B4 = generic_symmetric("b")
spatial = sp.Poly(
    (H0 + lam * A3 + lam**2 * B4 + lam**4 * C6).det(
        method="berkowitz"
    ),
    lam,
)
quotient_determinant = c1 * c2 * c3
assert sp.expand(
    spatial.coeff_monomial(lam**14)
    - quotient_determinant * B4[3, 3]
) == 0
assert sp.expand(
    spatial.coeff_monomial(lam**13)
    - quotient_determinant * A3[3, 3]
) == 0


# 2. The isotropic bordered coefficient.
x, y, m, t = sp.symbols("x y m t")
u = (x, y, m)
s_generic = sp.Function("s")(x, y, m)
gradient_generic = sp.Matrix([sp.diff(s_generic, variable) for variable in u])
hessian_generic = sp.hessian(s_generic, u)
bordered_generic = sp.zeros(4)
bordered_generic[0, 1:] = gradient_generic.T
bordered_generic[1:, 0] = gradient_generic
bordered_generic[1:, 1:] = t * hessian_generic
bordered_invariant = (
    gradient_generic.T * hessian_generic.adjugate() * gradient_generic
)[0]
assert sp.expand(
    bordered_generic.det(method="berkowitz")
    + t**2 * bordered_invariant
) == 0


# 3. Normalize a nonzero ternary cubic with singular Hessian to a binary
# cubic.  The rank-two binary orbits are x*y*(x+y) and x^2*y.  Add a
# general quadratic and linear part and reconstruct the coefficient ideal
# of the bordered invariant.
qxx, qxy, qxm, qyy, qym = sp.symbols("qxx qxy qxm qyy qym")
lx, ly, lm = sp.symbols("lx ly lm")
parameters = (qxx, qxy, qxm, qyy, qym, lx, ly, lm)
quadratic = (
    qxx * x**2
    + 2 * qxy * x * y
    + 2 * qxm * x * m
    + qyy * y**2
    + 2 * qym * y * m
) / 2
linear = lx * x + ly * y + lm * m


def bordered_coefficients(cubic: sp.Expr) -> list[sp.Expr]:
    polynomial = linear + quadratic + cubic
    gradient = sp.Matrix(
        [sp.diff(polynomial, variable) for variable in u]
    )
    hessian = sp.hessian(polynomial, u)
    invariant = sp.Poly(
        sp.expand((gradient.T * hessian.adjugate() * gradient)[0]),
        x,
        y,
        m,
    )
    coefficients: list[sp.Expr] = []
    for _, coefficient in invariant.terms():
        numerator = sp.together(coefficient).as_numer_denom()[0]
        numerator = sp.expand(numerator)
        if numerator != 0 and numerator not in coefficients:
            coefficients.append(numerator)
    return coefficients


rank_two_radical = [qxm, qym, lm]
rank_two_square = [
    qxm**2,
    qxm * qym,
    qxm * lm,
    qym**2,
    qym * lm,
    lm**2,
]
for cubic in (x * y * (x + y), x**2 * y):
    coefficient_ideal = sp.groebner(
        bordered_coefficients(cubic),
        *parameters,
        order="grevlex",
    )
    square_ideal = sp.groebner(
        rank_two_square,
        *parameters,
        order="grevlex",
    )
    for generator in coefficient_ideal.polys:
        assert square_ideal.reduce(generator.as_expr())[1] == 0
    for generator in square_ideal.polys:
        assert coefficient_ideal.reduce(generator.as_expr())[1] == 0


# 4. The rank-one binary orbit is x^3.  Its coefficient ideal has radical
#
#   P1 intersection P2,
#   P1=(qym,qyy,qxm*ly-qxy*lm),  P2=(lm,qym,qxm).
#
# Instead of trusting a primary-decomposition black box, verify both
# radical inclusions exactly.  The displayed four generators generate the
# intersection set-theoretically; the elementary branch split is described
# in the canonical note.
rank_one_coefficients = bordered_coefficients(x**3)
rank_one_ideal = sp.groebner(
    rank_one_coefficients,
    *parameters,
    order="grevlex",
)
intersection_generators = [
    qym,
    qyy * lm,
    qxm * ly - qxy * lm,
    qxm * qyy,
]
intersection_ideal = sp.groebner(
    intersection_generators,
    *parameters,
    order="grevlex",
)

# I is contained in the intersection ideal.
for equation in rank_one_coefficients:
    assert intersection_ideal.reduce(equation)[1] == 0

# The intersection ideal is contained in radical(I), with explicit powers.
radical_powers = (2, 3, 2, 3)
for generator, power in zip(intersection_generators, radical_powers):
    assert rank_one_ideal.reduce(sp.expand(generator**power))[1] == 0


# 5. If the cubic part is zero, write s=ell+a_2.  The degree-two and
# degree-zero pieces of the bordered invariant are the standard adjugate
# identities used to choose a constant direction.
u1, u2, u3, scale = sp.symbols("u1 u2 u3 scale")
ell1, ell2, ell3 = sp.symbols("ell1 ell2 ell3")
d11, d12, d13, d22, d23, d33 = sp.symbols(
    "d11 d12 d13 d22 d23 d33"
)
D = sp.Matrix(
    [
        [d11, d12, d13],
        [d12, d22, d23],
        [d13, d23, d33],
    ]
)
vector = sp.Matrix([u1, u2, u3])
ell_vector = sp.Matrix([ell1, ell2, ell3])
quadratic_form = (vector.T * D * vector)[0] / 2
quadratic_s = (ell_vector.T * vector)[0] + quadratic_form
quadratic_gradient = sp.Matrix(
    [sp.diff(quadratic_s, variable) for variable in (u1, u2, u3)]
)
quadratic_invariant = (
    quadratic_gradient.T * D.adjugate() * quadratic_gradient
)[0]
scaled_invariant = sp.Poly(
    sp.expand(
        quadratic_invariant.subs(
            {u1: scale * u1, u2: scale * u2, u3: scale * u3}
        )
    ),
    scale,
)
assert sp.expand(
    scaled_invariant.coeff_monomial(scale**2)
    - 2 * quadratic_form * D.det()
) == 0
assert sp.expand(
    scaled_invariant.coeff_monomial(1)
    - (ell_vector.T * D.adjugate() * ell_vector)[0]
) == 0


print("PASS: rank-three Hess(h_6) forces D_t^2 h_4=D_t^2 h_3=0")
print("PASS: the isotropic t^2 coefficient is the cubic bordered invariant")
print("PASS: both rank-two binary-cubic orbits force s to omit one variable")
print("PASS: the rank-one cubic radical has two constant-direction branches")
print("PASS: the quadratic-leading boundary has a constant direction")
print("PASS: the shared cotangent-lift and terminal HC(2) blocks hold")
print(
    "SCOPE: with the cited low-dimensional theorems, this excludes the "
    "full Hess(h_6)-rank-three cubic--quartic--sextic chart"
)
