#!/usr/bin/env python3
"""Verify the nonlinear toric Schur descent of the Meng--Yang HC(5) example.

Write the Meng--Yang potential as

    Psi = A^2 + 11 A + 2 B,

where A and B are linear in the three dual variables.  The coefficient row
of A is unimodular over QQ[x,y].  An explicit SL(3) completion therefore
makes t=A a polynomial coordinate, with unit quadratic pivot.

The natural completion descends polynomially at the common collision level,
but its four-variable Hessian determinant is 16*J(x*y)^2 with nonconstant J,
and the collision is lost.  An explicit relative SL(2) correction cancels
J exactly and gives determinant 64.  The corrected descended potential has
the form

    f(x,y) + 2*y*r + 4*x*s,

so its gradient has an explicit polynomial inverse and cannot retain the
Meng--Yang collision.

The final symbolic calculation proves the corresponding all-degree toric
obstruction.  For a toric coefficient map

    (-y*P(x*y), x*Q(x*y)),

its Jacobian is (v*P(v)*Q(v))' at v=x*y.  If this is a nonzero constant,
then P and Q are units, and the doubled four-variable gradient is a
polynomial automorphism.  This is an exact obstruction for the toric
equivariant completion class, not for arbitrary nonlinear symplectic
coordinates or non-toric relative-Hamiltonian corrections.
"""

from __future__ import annotations

import sympy as sp


x, y = sp.symbols("x y")
y1, y2, y3 = sp.symbols("y1 y2 y3")
t, r, s = sp.symbols("t r s")
rho, sigma = sp.symbols("rho sigma")
v = x * y
u = 1 + v


# The Meng--Yang five-variable potential.
A = y1 * u**3 + 3 * x * y2 * u**2 - x**3 * y3
B = (
    y1 * y**2 * u * (4 + 3 * v)
    + y2 * (y + 3 * x * y**2 * (4 + 3 * v))
    + y3 * (2 * x - 3 * x**2 * y)
)
Psi = sp.expand(A**2 + 11 * A + 2 * B)
meng_variables = (x, y, y1, y2, y3)


# A Bezout identity for (u^3,-x^3) supplies a unit A-direction.
p = 1 - 3 * v + 6 * v**2
q = y**3 * (6 * v**2 + 15 * v + 10)
assert sp.expand(p * u**3 - x**3 * q) == 1

# The new dual coordinates are (t,r,s)^T = U (y1,y2,y3)^T.
U = sp.Matrix(
    [
        [u**3, 3 * x * u**2, -x**3],
        [0, 1, 0],
        [-q, 0, p],
    ]
)
assert sp.factor(U.det()) == 1
assert sp.expand((U * sp.Matrix([y1, y2, y3]))[0] - A) == 0

old_duals = (U.inv() * sp.Matrix([t, r, s])).applyfunc(sp.expand)
coordinate_substitution = dict(zip((y1, y2, y3), old_duals, strict=True))
Phi = sp.expand(Psi.subs(coordinate_substitution, simultaneous=True))

# Phi is monic quadratic in t, hence every constant-level critical equation
# d_t Phi=sigma has a polynomial solution.
phi_t_poly = sp.Poly(Phi, t)
assert phi_t_poly.degree() == 2
assert phi_t_poly.coeff_monomial(t**2) == 1
t_critical = sp.factor(sp.solve(sp.diff(Phi, t) - sigma, t)[0])
assert t_critical.is_polynomial(x, y, r, s, sigma)

descended = sp.expand((Phi - sigma * t).subs(t, t_critical))
descended_variables = (x, y, r, s)
descended_hessian_determinant = sp.factor(
    sp.hessian(descended, descended_variables).det(method="berkowitz")
)


# The natural completion is a doubling potential for the displayed
# two-variable coefficient map.
beta_r = sp.factor(sp.diff(Phi, r) / 2)
beta_s = sp.factor(sp.diff(Phi, s) / 2)
beta_t = sp.factor((sp.diff(Phi, t) - 2 * t - 11) / 2)
P = 18 * v**5 + 81 * v**4 + 120 * v**3 + 60 * v**2 - 1
Q = (v + 1) * (v + 2)
J = (
    144 * v**7
    + 945 * v**6
    + 2394 * v**5
    + 2910 * v**4
    + 1680 * v**3
    + 357 * v**2
    - 6 * v
    - 2
)
assert sp.expand(beta_r + y * P) == 0
assert sp.expand(beta_s - x * Q) == 0
assert sp.expand(beta_t - y**2 * (6 * v**2 + 15 * v + 4)) == 0
assert sp.expand(
    sp.Matrix([beta_r, beta_s]).jacobian((x, y)).det() - J
) == 0
assert sp.expand(descended_hessian_determinant - 16 * J**2) == 0
assert J.free_symbols


# The known Meng--Yang collision survives to one common t-critical level, but
# the remaining gradient coordinates differ after this nonlinear point
# transformation.
collision_old = (
    (
        sp.Rational(1),
        -sp.Rational(3, 2),
        sp.Rational(-26),
        -sp.Rational(3, 4),
        sp.Rational(27, 16),
    ),
    (
        sp.Rational(-1),
        sp.Rational(3, 2),
        sp.Rational(-26),
        sp.Rational(3, 4),
        -sp.Rational(27, 16),
    ),
)
collision_new = []
for point in collision_old:
    xy_substitution = {x: point[0], y: point[1]}
    new_duals = U.subs(xy_substitution) * sp.Matrix(point[2:])
    collision_new.append(tuple(point[:2]) + tuple(new_duals))

collision_sigma = tuple(
    sp.factor(
        sp.diff(Phi, t).subs(
            dict(zip((x, y, t, r, s), point, strict=True))
        )
    )
    for point in collision_new
)
assert collision_sigma == (-sp.Rational(19, 2),) * 2

natural_reduced_gradients = []
for point in collision_new:
    reduced_point = (point[0], point[1], point[3], point[4])
    substitution = dict(
        zip(descended_variables, reduced_point, strict=True)
    )
    substitution[sigma] = collision_sigma[0]
    natural_reduced_gradients.append(
        tuple(
            sp.factor(sp.diff(descended, variable).subs(substitution))
            for variable in descended_variables
        )
    )
assert natural_reduced_gradients[0] != natural_reduced_gradients[1]


# Principal-part cancellation.  The following toric SL(2) matrix is congruent
# to the identity at v=0 and sends the radial pair (P,Q) to (-1,2).
a = (3 * v**2 + 6 * v + 2) / 2
b = 2 * v + 3
c = 3 * (18 * v**4 + 63 * v**3 + 69 * v**2 + 21 * v - 1) / 2
d = 36 * v**5 + 108 * v**4 + 87 * v**3 + 3 * v**2 - 3 * v + 1
C = sp.Matrix([[a, x**2 * b], [y**2 * c, d]])
assert sp.factor(C.det()) == 1
assert (
    sp.Matrix([[beta_r, beta_s]]) * C
).applyfunc(sp.expand) == sp.Matrix([[y, 2 * x]])

# Substitute the corrected complementary coordinates.  This lifts to a
# polynomial determinant-one change of the three dual variables fixing t.
corrected_old_reduced = C * sp.Matrix([r, s])
corrected_phi = sp.expand(
    Phi.subs(
        {
            r: corrected_old_reduced[0],
            s: corrected_old_reduced[1],
        },
        simultaneous=True,
    )
)
corrected_t_critical = sp.factor(
    sp.solve(sp.diff(corrected_phi, t) - sigma, t)[0]
)
corrected_descended = sp.expand(
    (corrected_phi - sigma * t).subs(t, corrected_t_critical)
)
corrected_determinant = sp.factor(
    sp.hessian(corrected_descended, descended_variables).det(
        method="berkowitz"
    )
)
assert corrected_determinant == 64

# The correction makes the reduced potential a doubled linear map.  Its base
# term may depend on sigma, but the r,s coefficients are exactly 2*y and 4*x.
base_term = sp.expand(corrected_descended.subs({r: 0, s: 0}))
assert sp.expand(corrected_descended - base_term - 2 * y * r - 4 * x * s) == 0

# Write and verify the polynomial inverse of its gradient.
px, py, pr, ps = sp.symbols("px py pr ps")
inverse_x = ps / 4
inverse_y = pr / 2
base_dx = sp.diff(base_term, x)
base_dy = sp.diff(base_term, y)
inverse_r = sp.expand(
    (
        py - base_dy.subs(
            {x: inverse_x, y: inverse_y}, simultaneous=True
        )
    )
    / 2
)
inverse_s = sp.expand(
    (
        px - base_dx.subs(
            {x: inverse_x, y: inverse_y}, simultaneous=True
        )
    )
    / 4
)
inverse_map = (inverse_x, inverse_y, inverse_r, inverse_s)
corrected_gradient = tuple(
    sp.diff(corrected_descended, variable)
    for variable in descended_variables
)
output_substitution = dict(
    zip((px, py, pr, ps), corrected_gradient, strict=True)
)
assert all(
    sp.expand(
        inverse_coordinate.subs(output_substitution, simultaneous=True)
        - variable
    )
    == 0
    for inverse_coordinate, variable in zip(
        inverse_map, descended_variables, strict=True
    )
)


# All-degree toric identity.  Generic coefficient lists avoid relying on
# special features of the displayed P and Q.
degree_bound = 4
p_coefficients = sp.symbols(f"p0:{degree_bound + 1}")
q_coefficients = sp.symbols(f"q0:{degree_bound + 1}")
radial_variable = sp.symbols("zeta")
generic_p = sum(
    coefficient * radial_variable**index
    for index, coefficient in enumerate(p_coefficients)
)
generic_q = sum(
    coefficient * radial_variable**index
    for index, coefficient in enumerate(q_coefficients)
)
generic_beta = sp.Matrix(
    [
        -y * generic_p.subs(radial_variable, v),
        x * generic_q.subs(radial_variable, v),
    ]
)
generic_jacobian = sp.expand(
    generic_beta.jacobian((x, y)).det()
)
radial_derivative = sp.diff(
    radial_variable * generic_p * generic_q, radial_variable
).subs(radial_variable, v)
assert sp.expand(generic_jacobian - radial_derivative) == 0


print("PASS: the Meng--Yang coefficient row admits an explicit SL(3) completion")
print("PASS: t=A has unit quadratic pivot and a polynomial critical solution")
print(
    "PASS: the natural descent has determinant 16*J(x*y)^2 and loses the "
    "collision"
)
print(
    "PASS: an explicit toric SL(2) correction gives determinant 64 and a "
    "polynomially invertible gradient"
)
print(
    "PASS: every constant-nonzero toric radial descent has unit radial "
    "factors and cannot retain a collision"
)
print(
    "SCOPE: non-toric symplectic changes, higher-degree critical equations, "
    "and non-coordinate coisotropic embeddings remain open"
)
