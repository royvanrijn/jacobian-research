#!/usr/bin/env python3
"""Exact quartic-sheet weighted map and its C!=0 root reconstruction."""

import sympy as sp


x, y, z = sp.symbols("x y z")
A, B, C, W = sp.symbols("A B C W")

u = 1 + 3*x*y
gamma = 1 - 4*x*y - x**2*z
w_source = u*gamma

# Integer normalization of the seed p(w)=w-2w^3, H=(w^2-w^4)/2.
G = (
    sp.cancel((2*u + u**2 - 3*u**4*gamma**2)/x**2),
    sp.cancel((1 + u - 2*u**3*gamma**2)/x),
    sp.expand(x*gamma),
)
assert all(sp.denom(component) == 1 for component in G)

determinant = sp.factor(sp.Matrix(G).jacobian((x, y, z)).det())
degrees = tuple(sp.Poly(component, x, y, z).total_degree() for component in G)
assert determinant == -6
assert degrees == (12, 11, 4)

points = ((1, 0, 0), (-1, 0, 2))
target = (0, 0, 1)
for point in points:
    image = tuple(sp.expand(component.subs(dict(zip((x, y, z), point)))) for component in G)
    assert image == target
assert points[0] != points[1]


# The inverse polynomial has degree four:
#
#   E(W)=W^2-W^4-2BCW+AC^2.
#
# On the source, E(w)=0 and E'(w)=-2 gamma=-2C/x.
E = W**2 - W**4 - 2*B*C*W + A*C**2
dE = sp.diff(E, W)
assert sp.factor(E.subs({A: G[0], B: G[1], C: G[2], W: w_source})) == 0
assert sp.factor(dE.subs({A: G[0], B: G[1], C: G[2], W: w_source}) + 2*gamma) == 0


# Conversely, for C!=0 every simple root reconstructs exactly one source point.
gamma_root = B*C - W + 2*W**3
assert sp.factor(gamma_root + dE/2) == 0
x_root = C/gamma_root
u_root = W/gamma_root
v_root = (u_root - 1)/3
y_root = v_root/x_root
s_root = 1 - 4*v_root - gamma_root
z_root = s_root/x_root**2

reconstruction = {x: x_root, y: y_root, z: z_root}
for got, want in zip(G, (A, B, C)):
    numerator = sp.factor(sp.together(got.subs(reconstruction) - want).as_numer_denom()[0])
    remainder = sp.rem(sp.Poly(numerator, W), sp.Poly(E, W)).as_expr()
    assert sp.factor(remainder) == 0

# These identities recover the construction coordinates, proving uniqueness.
assert sp.factor((1 + 3*x_root*y_root) - u_root) == 0
assert sp.factor((1 - 4*x_root*y_root - x_root**2*z_root) - gamma_root) == 0
assert sp.factor(u_root*gamma_root - W) == 0

# No hidden reconstruction denominator occurs on C!=0: gamma_root=-E'/2, so
# the only excluded roots are exactly the repeated roots E'=0.
assert sp.factor(x_root + 2*C/dE) == 0


# A concrete target has four distinct roots, certifying generic degree four.
# At (A,B,C)=(1,0,1), E=-(W^4-W^2-1), which is squarefree.
quartic_control = sp.factor(E.subs({A: 1, B: 0, C: 1}))
assert quartic_control == -W**4 + W**2 + 1
assert sp.gcd(sp.Poly(quartic_control, W), sp.Poly(sp.diff(quartic_control, W), W)).degree() == 0

print("PASS: quartic weighted map is polynomial with degrees", degrees)
print("PASS: det(DG) =", determinant, "and the stored collision is exact")
print("PASS: E(W)=W^2-W^4-2BCW+AC^2 is the inverse quartic")
print("PASS: every simple root for C!=0 reconstructs uniquely")
print("PASS: repeated roots are exactly the C!=0 reconstruction poles")


# The quartic map is the tau=0 degree-drop specialization of the normalized
# degree-five seed surface at kappa=-5:
#
#   H=w^2(w-1)((kappa/2+2)w-kappa/2-3)=(w^2-w^4)/2.
#
# Specialize the uniform rank-two descent formulas and verify that the
# resulting four-dimensional map is polynomial symplectic and remains
# noninvertible.
X, Q, Z, E_adapted = sp.symbols("X Q Z E_adapted")
kappa = sp.Integer(-5)
tau = sp.Integer(0)
a = -(1 + kappa) / (2 + kappa)
H_drop = sp.factor(
    W**2
    * (W - 1)
    * (
        tau * W**2
        + (kappa / 2 - 2 * tau + 2) * W
        - kappa / 2
        + tau
        - 3
    )
)
assert a == -sp.Rational(4, 3)
assert sp.expand(H_drop - (W**2 - W**4) / 2) == 0

drop_shear = sp.factor(
    (
        12 * (kappa + 1) * tau**2
        - 18 * (kappa + 1) * (kappa + 6) * tau
        + 9 * (kappa**3 + 16 * kappa**2 + 52 * kappa + 72)
    )
    / (28 * (kappa + 2))
)
assert drop_shear == -sp.Rational(261, 28)

p_drop = sp.diff(H_drop, W)
q_drop = sp.expand(W * p_drop - H_drop)
W_drop = Z + drop_shear * Q**2
Y_drop = Q - X * W_drop / 3
source_v = -3 * X * Y_drop / (2 * a)
source_gamma = 1 - 3 * X * Q / 2
source_u = 1 + source_v
marked = sp.expand(source_u * source_gamma)
A_drop = (source_u + q_drop.subs(W, marked) / source_gamma**2) / X**2
B_drop = (1 + p_drop.subs(W, marked) / source_gamma) / X
S_drop = sp.cancel(-2 * a * A_drop / 3)
T_drop = sp.cancel(B_drop)
R_drop = sp.expand(2 * X * source_gamma)
assert all(
    not ({X, Q, Z} & sp.denom(component).free_symbols)
    for component in (S_drop, T_drop, R_drop)
)
assert R_drop == 2 * X - 3 * X**2 * Q


def quotient_bracket(left, right):
    return sp.cancel(
        -3
        * X**2
        * (
            sp.diff(left, X) * sp.diff(right, Z)
            - sp.diff(left, Z) * sp.diff(right, X)
        )
        + (6 * X * Q - 2)
        * (
            sp.diff(left, Q) * sp.diff(right, Z)
            - sp.diff(left, Z) * sp.diff(right, Q)
        )
    )


base_jacobian = sp.Matrix([S_drop, T_drop, R_drop]).jacobian((X, Q, Z))
assert sp.factor(base_jacobian.det()) == -1
assert quotient_bracket(S_drop, T_drop) == 1
assert quotient_bracket(R_drop, S_drop) == 0
assert quotient_bracket(R_drop, T_drop) == 0

# Recover the polynomial Hamiltonian completing D=E+f.  This is the same
# finite homotopy used by the uniform descent, now specialized before any
# large parameter-field simplification.
w_family = (base_jacobian.adjugate() / base_jacobian.det())[:, 2].applyfunc(
    sp.cancel
)
w_E = sp.Matrix([(1 + 3 * X * Q) / 2, -3 * Q**2, 9 * Q * Z / 2])
difference = (w_family - w_E).applyfunc(sp.cancel)


def quotient_hamiltonian(f):
    return sp.Matrix(
        [
            3 * X**2 * sp.diff(f, Z),
            (2 - 6 * X * Q) * sp.diff(f, Z),
            -3 * X**2 * sp.diff(f, X)
            + (6 * X * Q - 2) * sp.diff(f, Q),
        ]
    )


f_zero = sp.integrate(sp.cancel(difference[0] / (3 * X**2)), Z)
residual = (difference - quotient_hamiltonian(f_zero)).applyfunc(sp.cancel)
assert residual[0] == residual[1] == 0
inverse_x, rho = sp.symbols("inverse_x rho")
residual_vrho = sp.cancel(
    residual[2].subs(
        {
            X: 1 / inverse_x,
            Q: inverse_x * (2 - rho * inverse_x) / 3,
        }
    )
)
h_vrho = sp.integrate(residual_vrho / 3, inverse_x)
h_drop = sp.cancel(
    h_vrho.subs(
        {
            inverse_x: 1 / X,
            rho: 2 * X - 3 * X**2 * Q,
        }
    )
)
f_drop = sp.cancel(f_zero + h_drop)
assert sp.denom(f_drop) in (1, -1)
assert all(
    sp.cancel(difference[index] - quotient_hamiltonian(f_drop)[index]) == 0
    for index in range(3)
)
D_drop = E_adapted + f_drop


def adapted_poisson(left, right):
    base = quotient_bracket(left, right)
    e_vector = ((1 + 3 * X * Q) / 2, -3 * Q**2, 9 * Q * Z / 2)
    e_left = sp.diff(left, E_adapted)
    e_right = sp.diff(right, E_adapted)
    connection = sum(
        e_vector[index]
        * (
            e_left * sp.diff(right, variable)
            - e_right * sp.diff(left, variable)
        )
        for index, variable in enumerate((X, Q, Z))
    )
    return sp.cancel(base + connection)


assert adapted_poisson(D_drop, R_drop) == 1
assert adapted_poisson(S_drop, T_drop) == 1
assert adapted_poisson(R_drop, S_drop) == 0
assert adapted_poisson(R_drop, T_drop) == 0
assert adapted_poisson(D_drop, S_drop) == 0
assert adapted_poisson(D_drop, T_drop) == 0
adapted_output_jacobian = sp.factor(
    sp.Matrix([R_drop, T_drop, D_drop, S_drop])
    .jacobian((X, Q, Z, E_adapted))
    .det()
)
assert adapted_output_jacobian == -1

# The adapted Poisson coordinates are themselves a polynomial coordinate
# system on canonical A^4.  This makes the completed map an honest
# polynomial symplectic endomorphism of standard affine four-space.
x4, q4, p4, z4 = sp.symbols("x4 q4 p4 z4")
adapted_source = {
    X: x4,
    Q: q4,
    Z: 3 * x4**2 * p4 + (2 - 6 * x4 * q4) * z4,
    E_adapted: (1 + 3 * x4 * q4) * p4 / 2 - 3 * q4**2 * z4,
}
adapted_coordinate_jacobian = sp.factor(
    sp.Matrix(tuple(adapted_source.values()))
    .jacobian((x4, q4, p4, z4))
    .det()
)
assert adapted_coordinate_jacobian == -1
p4_inverse = -6 * E_adapted * X * Q + 2 * E_adapted + 3 * Z * Q**2
z4_inverse = -3 * E_adapted * X**2 + 3 * Z * X * Q / 2 + Z / 2
for got, expected in zip(
    tuple(adapted_source.values()),
    (X, Q, Z, E_adapted),
    strict=True,
):
    assert sp.expand(
        got.subs(
            {
                x4: X,
                q4: Q,
                p4: p4_inverse,
                z4: z4_inverse,
            },
            simultaneous=True,
        )
        - expected
    ) == 0
assert adapted_output_jacobian * adapted_coordinate_jacobian == 1

# Identify the base triple with the quartic weighted map after explicit
# diagonal and triangular source changes.  This transports the generic
# degree-four inverse equation and the stored collision.
weighted_substitution = {
    x: X,
    y: sp.Rational(3, 8) * Y_drop,
    z: W_drop / 2,
}
G_in_drop_coordinates = tuple(
    sp.expand(component.subs(weighted_substitution))
    for component in G
)
assert sp.expand(S_drop - sp.Rational(4, 9) * G_in_drop_coordinates[0]) == 0
assert sp.expand(T_drop - G_in_drop_coordinates[1]) == 0
assert sp.expand(R_drop - 2 * G_in_drop_coordinates[2]) == 0

# The two quartic collision points become two distinct points of the
# completed symplectic map after choosing E so that D=0.
drop_base_points = (
    (sp.Integer(1), sp.Integer(0), sp.Integer(0)),
    (-sp.Integer(1), -sp.Rational(4, 3), sp.Rational(144, 7)),
)
completed_points = []
for base_point in drop_base_points:
    base_substitution = dict(zip((X, Q, Z), base_point, strict=True))
    e_value = -sp.factor(f_drop.subs(base_substitution))
    completed_point = (*base_point, e_value)
    completed_points.append(completed_point)
    image = tuple(
        sp.factor(
            component.subs(
                dict(
                    zip(
                        (X, Q, Z, E_adapted),
                        completed_point,
                        strict=True,
                    )
                )
            )
        )
        for component in (R_drop, T_drop, D_drop, S_drop)
    )
    assert image == (2, 0, 0, 0)
assert completed_points[0] != completed_points[1]

drop_orders = (
    sp.Poly(S_drop, X, Q, Z).degree(Z),
    sp.Poly(T_drop, X, Q, Z).degree(Z),
)
drop_degrees = tuple(
    sp.Poly(component, X, Q, Z, E_adapted).total_degree()
    for component in (R_drop, T_drop, D_drop, S_drop)
)
canonical_drop_degrees = tuple(
    sp.Poly(
        sp.expand(component.subs(adapted_source)),
        x4,
        q4,
        p4,
        z4,
    ).total_degree()
    for component in (R_drop, T_drop, D_drop, S_drop)
)

print("PASS: (kappa,tau)=(-5,0) is the normalized quartic degree-drop seed")
print("PASS: the specialized rank-two Hamiltonian is polynomial")
print("PASS: all six completed Poisson brackets hold exactly over Q")
print("PASS: the adapted source coordinates are polynomial canonical coordinates")
print("PASS: the completed map has two distinct points over (2,0,0,0)")
print("PASS: left-right equivalence transports generic degree four")
print("DEGREE_DROP_Z_ORDERS=", drop_orders)
print("DEGREE_DROP_ADAPTED_TOTAL_DEGREES=", drop_degrees)
print("DEGREE_DROP_CANONICAL_TOTAL_DEGREES=", canonical_drop_degrees)
