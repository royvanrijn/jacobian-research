#!/usr/bin/env python3
"""Exact audit of the normalized PGL2 action and cubic dual-number tangent."""

import sympy as sp

a, b, c, d, e, y, z = sp.symbols("a b c d e y z")
factor_variables = (a, b, c, d, e)
middle = a * d + b * c
resultant = a**2 * e - a * b * d + b**2 * c

b_forward = 1 + a * y
c_forward = 1 - sp.Rational(3, 2) * a * y + a**2 * z
d_forward = sp.Rational(1, 2) * y - a * z + sp.Rational(3, 2) * a * y**2 - a**2 * y * z
e_forward = -2 * z + 4 * y**2 - 4 * a * y * z + 3 * a * y**3 - 2 * a**2 * y**2 * z
forward = {b: b_forward, c: c_forward, d: d_forward, e: e_forward}

y_inverse = 2 * b * d - a * e
z_inverse = 2 * d**2 + c * e + 6 * b * d**2 + 3 * b * c * e - sp.Rational(9, 2) * e

# E=T*d/dS, F=S*d/dT, H=T*d/dT-S*d/dS on binary forms.
raw_actions = {
    "E": sp.Matrix((b, 0, d, 2 * e, 0)),
    "F": sp.Matrix((0, a, 0, 2 * c, d)),
    "H": sp.Matrix((a, -b, 2 * c, 0, -2 * e)),
}


def derivative(polynomial, vector):
    return sp.expand(sum(
        sp.diff(polynomial, variable) * value
        for variable, value in zip(factor_variables, vector)
    ))


# If kappa is the raw variation of middle, the unique infinitesimal
# renormalization to middle=resultant=1 is L -> (1+t*kappa)L and
# Q -> (1-2*t*kappa)Q.
corrected_actions = {}
chart_actions = {}
for name, raw in raw_actions.items():
    assert derivative(resultant, raw) == 0
    kappa = derivative(middle, raw)
    correction = sp.Matrix((
        kappa * a, kappa * b, -2 * kappa * c, -2 * kappa * d, -2 * kappa * e,
    ))
    corrected = sp.expand(raw + correction)
    corrected_actions[name] = corrected
    assert sp.expand(derivative(middle, corrected).subs(forward)) == 0
    assert sp.expand(derivative(resultant, corrected).subs(forward)) == 0

    pulled = sp.Matrix((
        corrected[0], derivative(y_inverse, corrected), derivative(z_inverse, corrected),
    )).subs(forward).applyfunc(sp.factor)
    chart_actions[name] = pulled

    chart_jacobian = sp.Matrix((a, b_forward, c_forward, d_forward, e_forward)).jacobian((a, y, z))
    assert all(sp.expand(entry) == 0 for entry in chart_jacobian * pulled - corrected.subs(forward))


# P=p0*T^3+T^2*S+p2*T*S^2+p3*S^3 and dot(P)=rho(P)-kappa*P.
p0, p2, p3 = sp.symbols("p0 p2 p3")
target_variables = (p0, p2, p3)
target_actions = {
    "E": sp.Matrix((1 - 2 * p0 * p2, 3 * p3 - 2 * p2**2, -2 * p2 * p3)),
    "F": sp.Matrix((-3 * p0**2, 2 - 3 * p0 * p2, p2 - 3 * p0 * p3)),
    "H": sp.Matrix((2 * p0, -2 * p2, -4 * p3)),
}

G = sp.Matrix((
    sp.expand((a * c).subs(forward)),
    sp.expand((a * e + b * d).subs(forward)),
    sp.expand((b * e).subs(forward)),
))
for name in raw_actions:
    source_side = G.jacobian((a, y, z)) * chart_actions[name]
    target_side = target_actions[name].subs(dict(zip(target_variables, G)), simultaneous=True)
    assert all(sp.expand(entry) == 0 for entry in source_side - target_side)


# Conjugate to the announced coordinates:
# (a,y,z_chart)=(x,yy,-zz/2), (F1,F2,F3)=(p3,2*p2,2*p0).
x, yy, zz = sp.symbols("x yy zz")
source_substitution = {a: x, y: yy, z: -zz / 2}
announced_source_actions = {}
for name, vector in chart_actions.items():
    announced_source_actions[name] = sp.Matrix((
        vector[0].subs(source_substitution),
        vector[1].subs(source_substitution),
        -2 * vector[2].subs(source_substitution),
    )).applyfunc(sp.factor)

u, v, w = sp.symbols("u v w")
announced_target_actions = {
    "E": sp.Matrix((-u * v, 6 * u - v**2, 2 - v * w)),
    "F": sp.Matrix(((v - 3 * w * u) / 2, 4 - sp.Rational(3, 2) * w * v, -sp.Rational(3, 2) * w**2)),
    "H": sp.Matrix((-4 * u, -2 * v, 2 * w)),
}

one_plus_xy = 1 + x * yy
foundational = sp.Matrix((
    one_plus_xy**3 * zz + yy**2 * one_plus_xy * (4 + 3 * x * yy),
    yy + 3 * x * one_plus_xy**2 * zz + 3 * x * yy**2 * (4 + 3 * x * yy),
    2 * x - 3 * x**2 * yy - x**3 * zz,
))
source_variables = (x, yy, zz)
for name in raw_actions:
    source_side = foundational.jacobian(source_variables) * announced_source_actions[name]
    target_side = announced_target_actions[name].subs(
        {u: foundational[0], v: foundational[1], w: foundational[2]}, simultaneous=True,
    )
    assert all(sp.expand(entry) == 0 for entry in source_side - target_side)


deformation = sp.Matrix((
    sp.Rational(7, 12) * yy**2 + sp.Rational(1, 2) * x * yy * zz
    + sp.Rational(9, 4) * x * yy**3 + x**2 * yy**2 * zz
    + sp.Rational(3, 2) * x**2 * yy**4 + sp.Rational(1, 2) * x**3 * yy**3 * zz,
    x * zz + sp.Rational(11, 2) * x * yy**2 + 3 * x**2 * yy * zz
    + 6 * x**2 * yy**3 + 2 * x**3 * yy**2 * zz,
    0,
))

# H_def is not in the constant span of the normalized PGL2 directions.
alpha_e, alpha_f, alpha_h = sp.symbols("alpha_e alpha_f alpha_h")
orbit_combination = sp.zeros(3, 1)
for coefficient, name in zip((alpha_e, alpha_f, alpha_h), ("E", "F", "H")):
    orbit_combination += coefficient * foundational.jacobian(source_variables) * announced_source_actions[name]
equations = []
for entry in orbit_combination - deformation:
    equations.extend(sp.Poly(sp.expand(entry), *source_variables).coeffs())
matrix, vector = sp.linear_eq_to_matrix(equations, (alpha_e, alpha_f, alpha_h))
assert matrix.rank() < matrix.row_join(vector).rank()

# The full polynomial right-orbit tangent is all polynomial deformations for
# a Keller map: V=(DF)^(-1)H is polynomial.  Here it is divergence-free too.
jacobian = foundational.jacobian(source_variables)
assert jacobian.det() == -2
right_trivializer = (jacobian.adjugate() * deformation / jacobian.det()).applyfunc(sp.factor)
assert all(sp.expand(entry) == 0 for entry in jacobian * right_trivializer - deformation)
divergence = sum(
    sp.diff(right_trivializer[index], variable)
    for index, variable in enumerate(source_variables)
)
assert sp.expand(divergence) == 0

# The trivializer respects the original source weights (1,-1,-2), hence
# descends to the invariant plane k[xy,x^2*z].  Its highest ordinary-degree
# plane part has a compact factorization which rules out local nilpotence and
# local finiteness, so V_def is not itself an algebraic one-parameter-group
# generator.
weight_euler = sp.Matrix((x, -yy, -2 * zz))


def bracket(left, right, variables):
    return right.jacobian(variables) * left - left.jacobian(variables) * right


assert all(sp.expand(entry) == 0 for entry in bracket(weight_euler, right_trivializer, source_variables))

plane_u, plane_v = sp.symbols("plane_u plane_v")
induced_u = sp.expand(yy * right_trivializer[0] + x * right_trivializer[1])
induced_v = sp.expand(2 * x * zz * right_trivializer[0] + x**2 * right_trivializer[2])
to_plane = {yy: plane_u / x, zz: plane_v / x**2}
induced_u = sp.cancel(induced_u.subs(to_plane))
induced_v = sp.cancel(induced_v.subs(to_plane))
assert x not in induced_u.free_symbols and x not in induced_v.free_symbols


def highest_homogeneous(polynomial):
    poly = sp.Poly(polynomial, plane_u, plane_v)
    degree = max(sum(monomial) for monomial, _coefficient in poly.terms())
    part = sum(
        coefficient * plane_u**monomial[0] * plane_v**monomial[1]
        for monomial, coefficient in poly.terms()
        if sum(monomial) == degree
    )
    return degree, sp.factor(part)


degree_u, top_u = highest_homogeneous(induced_u)
degree_v, top_v = highest_homogeneous(induced_v)
common = sp.Rational(3, 2) * plane_u**4 * (3 * plane_u + plane_v) ** 2
assert degree_u == degree_v == 7
assert sp.expand(top_u + common * plane_u) == 0
assert sp.expand(top_v - common * (6 * plane_u + plane_v)) == 0

# If K=-u*d_u+(6u+v)*d_v, then K(common)=-2*common and
# (common*K)^n(u) is a nonzero scalar times common^n*u for every n.
k_on_common = sp.expand(
    -plane_u * sp.diff(common, plane_u)
    + (6 * plane_u + plane_v) * sp.diff(common, plane_v)
)
assert sp.expand(k_on_common + 2 * common) == 0

# Despite not being a one-parameter-group generator, V is tangent to an
# honest reduced curve of automorphisms.  Construct a vector potential
# V=(A3_y-A2_z,-A3_x,A2_x), polarize A3 in (x,y) over QQ[z] and A2 in
# (x,z) over QQ[y], and obtain a finite sum of locally nilpotent shears.
component_x, component_y, component_z = tuple(right_trivializer)
potential_2 = sp.integrate(component_z, x)
potential_3 = -sp.integrate(component_y, x) + sp.integrate(component_x.subs(x, 0), yy)
reconstructed = sp.Matrix((
    sp.diff(potential_3, yy) - sp.diff(potential_2, zz),
    -sp.diff(potential_3, x),
    sp.diff(potential_2, x),
))
assert all(sp.expand(entry) == 0 for entry in reconstructed - right_trivializer)


def polarize_binary(polynomial, first, second, invariant):
    """Write a binary polynomial as sum coefficient(invariant)*(first+lambda*second)^d."""
    binary_poly = sp.Poly(sp.expand(polynomial), first, second)
    maximum = max((sum(monomial) for monomial, _ in binary_poly.terms()), default=0)
    answer = []
    for degree in range(1, maximum + 1):
        coefficients = [
            binary_poly.coeff_monomial(first ** (degree - index) * second**index)
            for index in range(degree + 1)
        ]
        if not any(coefficients):
            continue
        lambdas = tuple(range(degree + 1))
        polarization_matrix = sp.Matrix([
            [sp.binomial(degree, index) * lam**index for lam in lambdas]
            for index in range(degree + 1)
        ])
        polarized = polarization_matrix.inv() * sp.Matrix(coefficients)
        for lam, coefficient in zip(lambdas, polarized):
            coefficient = sp.factor(coefficient)
            if coefficient == 0:
                continue
            assert coefficient.free_symbols <= {invariant}
            linear_form = first + lam * second
            answer.append((coefficient, linear_form, degree, lam))
    return answer


polarized_3 = polarize_binary(potential_3, x, yy, zz)
polarized_2 = polarize_binary(potential_2, x, zz, yy)
assert sp.expand(potential_3 - sum(coefficient * linear**degree for coefficient, linear, degree, _ in polarized_3)) == 0
assert sp.expand(potential_2 - sum(coefficient * linear**degree for coefficient, linear, degree, _ in polarized_2)) == 0
assert (len(polarized_3), len(polarized_2)) == (46, 87)

# Each summand gives g*(lambda*d_first-d_second) or
# g*(-lambda*d_first+d_second), where g is killed by the constant direction.
# Its exponential is therefore an elementary polynomial shear.  The product
# of these 133 shear flows is a reduced A1-curve alpha_t in Aut(A3) with
# alpha_0=id and alpha'_0=V_def.
for coefficient, linear, degree, lam in polarized_3:
    shear_coefficient = sp.expand(degree * coefficient * linear ** (degree - 1))
    assert sp.expand(lam * sp.diff(shear_coefficient, x) - sp.diff(shear_coefficient, yy)) == 0
for coefficient, linear, degree, lam in polarized_2:
    shear_coefficient = sp.expand(degree * coefficient * linear ** (degree - 1))
    assert sp.expand(-lam * sp.diff(shear_coefficient, x) + sp.diff(shear_coefficient, zz)) == 0

print("PASS: computed and normalized the three binary sl2 coefficient fields")
print("PASS: pulled the fields polynomially through the global A^3 chart")
print("PASS: computed the corresponding quadratic target cubic fields")
print("PASS: source and target fields are equivariant for G and the announced F")
print("PASS: H is not in the normalized PGL2 orbit distribution with constant coefficients")
print("PASS: H=DF*V for an explicit polynomial divergence-free vector field V")
print("PASS: V preserves the weight invariants but is neither locally nilpotent nor locally finite")
print("PASS: V is a sum of 133 LND shears and is tangent to a reduced automorphism curve")
