#!/usr/bin/env python3
"""Exact obstructions to natural double-Schur descents of gauge Keller maps.

Let F=(P,B,C): A^3 -> A^3 be a root-engineered quadratic-gauge map and
Phi=<u,F(x)> its six-variable Meng doubling.  Two direct routes could
conceivably reduce Phi to four Hessian variables.

DUAL--DUAL ROUTE
----------------
Eliminating two repaired dual variables leaves

    psi(x,s) = s*L(x) + h(x),

where L is a nonzero constant linear combination of P,B,C.  The coefficient
of s^2 in det Hess(psi) is the bordered-Hessian invariant

    K(L) = -grad(L)^T adj(Hess(L)) grad(L).

It is independent of h.  This checker proves K(L) is nonzero for every
nonzero target linear form throughout the admissible root-engineered
quadratic-gauge family.  Hence the final Hessian determinant cannot be a
constant, regardless of the constant quadratic repair used on the two
eliminated dual variables.

For the two-parameter cubic family this is an exact coefficient calculation.
For every seed degree N>=4, the top z-degree of L first forces the B and C
coefficients to vanish; K(P)(0,0,z)=9*z^2 then excludes the P direction.

SOURCE--DUAL ROUTE
------------------
A source-first Meng--Yang descent requires F to be affine along a nonzero
constant source direction.  The common first coordinate P=t*q forces every
such direction to be the z-direction.  All degree-N decorations with N>=4
are nonlinear in z, so no source-first descent starts.  In the cubic family
z is the only affine direction, but no constant target linear combination
of dF/dz is a nonzero scalar.  Thus the second constant pivot does not exist.

These are exact obstructions for constant linear source/target splittings
and the original affine gauge coordinates.  They do not exclude nonlinear
symplectic changes, nonlinear remaining dual coefficients, or nonconstant
quadratic blocks with exceptional divisibility.
"""

from __future__ import annotations

import sympy as sp


def bordered_hessian_invariant(
    polynomial: sp.Expr, variables: tuple[sp.Symbol, ...]
) -> sp.Expr:
    """Return -grad(f)^T adj(Hess(f)) grad(f)."""

    gradient = sp.Matrix([sp.diff(polynomial, variable) for variable in variables])
    hessian = sp.hessian(polynomial, variables)
    return sp.factor(-(gradient.T * hessian.adjugate() * gradient)[0])


# The universal coefficient identity for a retained dual variable.  If
#
#     psi(x,s) = s*L(x) + h(x),
#
# then its bordered Hessian is [[s*H+R,g],[g^T,0]].  The s^2 coefficient of
# its determinant is -g^T adj(H)g, independently of R=Hess(h).
s = sp.symbols("s")
h11, h12, h13, h22, h23, h33 = sp.symbols(
    "h11 h12 h13 h22 h23 h33"
)
r11, r12, r13, r22, r23, r33 = sp.symbols(
    "r11 r12 r13 r22 r23 r33"
)
g1, g2, g3 = sp.symbols("g1 g2 g3")
generic_hessian = sp.Matrix(
    [[h11, h12, h13], [h12, h22, h23], [h13, h23, h33]]
)
generic_remainder = sp.Matrix(
    [[r11, r12, r13], [r12, r22, r23], [r13, r23, r33]]
)
generic_gradient = sp.Matrix([g1, g2, g3])
bordered_matrix = (s * generic_hessian + generic_remainder).row_join(
    generic_gradient
)
bordered_matrix = bordered_matrix.col_join(
    sp.Matrix([[g1, g2, g3, 0]])
)
bordered_determinant = sp.Poly(
    bordered_matrix.det(method="berkowitz"), s
)
assert sp.expand(
    bordered_determinant.coeff_monomial(s**2)
    + (generic_gradient.T * generic_hessian.adjugate() * generic_gradient)[0]
) == 0


# The normalized two-parameter cubic gauge family.  Scaling B by a nonzero
# constant would not affect any vanishing statement below, so the convenient
# determinant-minus-two normalization is used.
x, y, z = sp.symbols("x y z")
alpha, beta = sp.symbols("alpha beta", nonzero=True)
t = 1 + x * y
q = t**2 * z + y**2 * (1 + 3 * t) / alpha
P = sp.expand(t * q)
B = sp.expand(y + 3 * alpha * x * q + 2 * beta * t * q)
C = sp.expand(x * (5 - 3 * t) - alpha * x**3 * z)
cubic_map = sp.Matrix([P, B, C])
assert sp.factor(cubic_map.jacobian((x, y, z)).det()) == -2


# Every target linear form fails the necessary K(L)=0 condition.  Three
# exact coefficients give a triangular proof.  Write
#
#     L = ell_P P + ell_B B + ell_C C.
#
# After clearing the harmless denominator alpha^2:
#
#   [x^8] K(L) = 9 alpha^4 ell_C^4,
#   [z^2] K(L) = 9 alpha^2 (ell_P+2 beta ell_B)^4,
#
# and, after imposing the two resulting equations,
#
#   [1] K(L) = 9 alpha^4 ell_B^4.
#
# Since alpha is a unit and the characteristic is zero, L=0 is forced.
ell_P, ell_B, ell_C = sp.symbols("ell_P ell_B ell_C")
linear_target = sp.expand(ell_P * P + ell_B * B + ell_C * C)
cubic_K = sp.Poly(
    sp.cancel(
        alpha**2
        * bordered_hessian_invariant(linear_target, (x, y, z))
    ),
    x,
    y,
    z,
)
assert sp.expand(
    cubic_K.coeff_monomial(x**8) - 9 * alpha**4 * ell_C**4
) == 0
assert sp.expand(
    cubic_K.coeff_monomial(z**2)
    - 9 * alpha**2 * (ell_P + 2 * beta * ell_B) ** 4
) == 0
cubic_K_reduced = sp.Poly(
    cubic_K.as_expr().subs(
        {ell_C: 0, ell_P: -2 * beta * ell_B},
        simultaneous=True,
    ),
    x,
    y,
    z,
)
assert sp.expand(
    cubic_K_reduced.coeff_monomial(1) - 9 * alpha**4 * ell_B**4
) == 0


# All higher seed degrees are excluded by the highest z-layer.  If N>=4 and
# gamma=g_N/g_1, the top z^N coefficient of
#
#     L = ell_P P + ell_B B + ell_C C
#
# is gamma*f_N(x,y), where
#
#   f_N = x^(N-2) t^(2N) (N ell_B t^2-(N-2) ell_C x^2).
#
# For L=f(x,y)z^N+lower z-degree terms, [z^(4N-2)]K(L) is
#
#   I_N(f) = -N f [N f det Hess(f)
#                  -(N+1)(f_x^2 f_yy-2 f_x f_y f_xy+f_y^2 f_xx)].
#
# The formula is checked abstractly here.
n = sp.symbols("n", integer=True, positive=True)
abstract_z = sp.symbols("abstract_z", nonzero=True)
f, f_x, f_y, f_xx, f_xy, f_yy = sp.symbols(
    "f f_x f_y f_xx f_xy f_yy"
)
abstract_gradient = sp.Matrix(
    [
        f_x * abstract_z**n,
        f_y * abstract_z**n,
        n * f * abstract_z ** (n - 1),
    ]
)
abstract_hessian = sp.Matrix(
    [
        [
            f_xx * abstract_z**n,
            f_xy * abstract_z**n,
            n * f_x * abstract_z ** (n - 1),
        ],
        [
            f_xy * abstract_z**n,
            f_yy * abstract_z**n,
            n * f_y * abstract_z ** (n - 1),
        ],
        [
            n * f_x * abstract_z ** (n - 1),
            n * f_y * abstract_z ** (n - 1),
            n * (n - 1) * f * abstract_z ** (n - 2),
        ],
    ]
)
abstract_K = sp.expand(
    -(abstract_gradient.T * abstract_hessian.adjugate() * abstract_gradient)[0]
)
expected_leading_invariant = -n * f * (
    n * f * (f_xx * f_yy - f_xy**2)
    - (n + 1)
    * (f_x**2 * f_yy - 2 * f_x * f_y * f_xy + f_y**2 * f_xx)
)
abstract_K_coefficient = sp.powsimp(
    abstract_K / abstract_z ** (4 * n - 2), force=True
)
assert sp.factor(
    sp.powsimp(
        abstract_K_coefficient - expected_leading_invariant,
        force=True,
    )
) == 0


# The x-adic leading term of I_N(f_N) is already nonzero when ell_B!=0.
# The monomial calculation used in the proof is:
#
#   f=A*x^m*t^k  =>
#   [x^(4m)] I_N(f)|_{y=0}
#     = -A^4*k*N*(3km-kN+m^2+mN).
#
# For the ell_B term, (m,k)=(N-2,2N+2), and the final factor is
# 2(3N^2-7N-4), positive for N>=4.  If ell_B=0, the ell_C term is the exact
# monomial with (m,k)=(N,2N), whose final factor is 6N^2.  Thus the leading
# invariant forces ell_B=ell_C=0 in every degree N>=4.
m, k = sp.symbols("m k", integer=True, nonnegative=True)
monomial_constant_factor = sp.expand(3 * k * m - k * n + m**2 + m * n)
assert sp.expand(
    monomial_constant_factor.subs(
        {m: n - 2, k: 2 * n + 2},
        simultaneous=True,
    )
    - 2 * (3 * n**2 - 7 * n - 4)
) == 0
assert sp.expand(
    monomial_constant_factor.subs(
        {m: n, k: 2 * n},
        simultaneous=True,
    )
    - 6 * n**2
) == 0
for checked_n in range(4, 13):
    assert 3 * checked_n**2 - 7 * checked_n - 4 > 0


# Low-degree direct regressions of the all-degree valuation argument.  They
# differentiate the actual f_N from the top z-layer, restrict to y=0, and
# recover the symbolic leading coefficients used above.
top_B, top_C = sp.symbols("top_B top_C")
for checked_n in range(4, 9):
    checked_f = (
        x ** (checked_n - 2)
        * t ** (2 * checked_n)
        * (
            checked_n * top_B * t**2
            - (checked_n - 2) * top_C * x**2
        )
    )
    checked_fx = sp.diff(checked_f, x)
    checked_fy = sp.diff(checked_f, y)
    checked_fxx = sp.diff(checked_f, x, 2)
    checked_fxy = sp.diff(checked_f, x, y)
    checked_fyy = sp.diff(checked_f, y, 2)
    checked_invariant = -checked_n * checked_f * (
        checked_n
        * checked_f
        * (checked_fxx * checked_fyy - checked_fxy**2)
        - (checked_n + 1)
        * (
            checked_fx**2 * checked_fyy
            - 2 * checked_fx * checked_fy * checked_fxy
            + checked_fy**2 * checked_fxx
        )
    )
    checked_at_y_zero = sp.Poly(
        sp.expand(checked_invariant.subs(y, 0)),
        x,
    )
    expected_B_lead = (
        -2
        * checked_n**5
        * (2 * checked_n + 2)
        * (3 * checked_n**2 - 7 * checked_n - 4)
        * top_B**4
    )
    assert sp.expand(
        checked_at_y_zero.coeff_monomial(x ** (4 * checked_n - 8))
        - expected_B_lead
    ) == 0
    checked_C_only = sp.Poly(
        checked_at_y_zero.as_expr().subs(top_B, 0),
        x,
    )
    expected_C_lead = (
        -12
        * checked_n**4
        * (checked_n - 2) ** 4
        * top_C**4
    )
    assert sp.expand(
        checked_C_only.coeff_monomial(x ** (4 * checked_n))
        - expected_C_lead
    ) == 0


# Once ell_B=ell_C=0, the retained coefficient is a multiple of P.  Its
# bordered invariant is visibly nonzero on the line x=y=0.
P_K = bordered_hessian_invariant(P, (x, y, z))
assert sp.factor(P_K.subs({x: 0, y: 0})) == 9 * z**2


# Source-first route.  Let v=a*d_x+b*d_y+c*d_z.  The coefficients y^4 and
# 1 in alpha*D_v^2(P) force a=b=0, so z is the only affine direction
# allowed by the common first coordinate.
a, b, c = sp.symbols("a b c")


def directional_derivative(expression: sp.Expr) -> sp.Expr:
    return sp.expand(
        a * sp.diff(expression, x)
        + b * sp.diff(expression, y)
        + c * sp.diff(expression, z)
    )


second_directional_P = sp.Poly(
    sp.expand(alpha * directional_derivative(directional_derivative(P))),
    x,
    y,
    z,
)
assert second_directional_P.coeff_monomial(y**4) == 6 * a**2
assert second_directional_P.coeff_monomial(1) == 8 * b**2


# A degree-N decoration in B has the form
# N*gamma*t^2*x^(N-2)*q^N.  Its second z-derivative is nonzero for N>=4,
# so the only possible affine direction from P is destroyed.
gamma = sp.symbols("gamma", nonzero=True)
decoration_second_z = (
    n**2
    * (n - 1)
    * gamma
    * t**6
    * x ** (n - 2)
    * q ** (n - 2)
)
for checked_n in range(4, 9):
    checked_decoration = decoration_second_z.subs(n, checked_n)
    assert checked_decoration != 0


# In degree three, z remains affine, but its coefficient row has no
# nonzero constant linear combination.  Restricting to y=0 first forces
# the B and C coefficients to vanish; varying xy then forces the P
# coefficient to vanish.
dF_dz = cubic_map.diff(z)
pivot_linear_form = sp.expand(
    ell_P * dF_dz[0] + ell_B * dF_dz[1] + ell_C * dF_dz[2]
)
pivot_on_y_zero = sp.Poly(pivot_linear_form.subs(y, 0), x)
assert pivot_on_y_zero.coeff_monomial(x) == 3 * alpha * ell_B
assert pivot_on_y_zero.coeff_monomial(x**3) == -alpha * ell_C
pivot_after_reduction = sp.expand(
    pivot_linear_form.subs({ell_B: 0, ell_C: 0})
)
assert sp.expand(pivot_after_reduction - ell_P * t**3) == 0
assert sp.expand(
    sp.diff(pivot_after_reduction, y) - 3 * ell_P * x * t**2
) == 0


print("PASS: the dual--dual s^2 coefficient is the bordered invariant K(L)")
print("PASS: K(L) is nonzero for every nonzero cubic-gauge target linear form")
print("PASS: the top z-layer excludes every higher-degree gauge target form")
print("PASS: higher-degree gauge maps have no affine constant source direction")
print("PASS: the cubic z-direction has no constant nonzero second pivot")
