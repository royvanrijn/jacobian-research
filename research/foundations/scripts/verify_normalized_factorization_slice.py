#!/usr/bin/env python3
"""Exact certificate for the normalized linear-quadratic factorization bridge."""

import sympy as sp


a, b, c, d, e = sp.symbols("a b c d e")
y, z = sp.symbols("y z")

resultant = a**2 * e - a * b * d + b**2 * c
middle = a * d + b * c

# Ambient coefficient-resultant map.
theta = (
    a * c,
    middle,
    a * e + b * d,
    b * e,
    resultant,
)
ambient_jacobian = sp.factor(
    sp.Matrix(theta).jacobian((a, b, c, d, e)).det()
)
assert sp.factor(ambient_jacobian + resultant**2) == 0
print("PASS: det DTheta = -Res(L,Q)^2")

# The affine normalization is the unique representative of the projective open.
m, r, alpha, beta = sp.symbols("m r alpha beta")
lambda_normal = m / r
mu_normal = r / m**2
assert sp.cancel(lambda_normal * mu_normal * m - 1) == 0
assert sp.cancel(lambda_normal**2 * mu_normal * r - 1) == 0

m_scaled = alpha * beta * m
r_scaled = alpha**2 * beta * r
lambda_scaled = m_scaled / r_scaled
mu_scaled = r_scaled / m_scaled**2
assert sp.cancel(alpha * lambda_scaled - lambda_normal) == 0
assert sp.cancel(beta * mu_scaled - mu_normal) == 0
print("PASS: projective factor pairs have a unique representative with m=Res=1")

# Forward polynomial coordinates on the normalized slice.
b_forward = 1 + a * y
c_forward = 1 - sp.Rational(3, 2) * a * y + a**2 * z
d_forward = (
    sp.Rational(1, 2) * y
    - a * z
    + sp.Rational(3, 2) * a * y**2
    - a**2 * y * z
)
e_forward = (
    -2 * z
    + 4 * y**2
    - 4 * a * y * z
    + 3 * a * y**3
    - 2 * a**2 * y**2 * z
)
forward = {
    b: b_forward,
    c: c_forward,
    d: d_forward,
    e: e_forward,
}

assert sp.expand((middle - 1).subs(forward)) == 0
assert sp.expand((resultant - 1).subs(forward)) == 0
print("PASS: forward formulas land on Res=1 and [LQ]_(T^2S)=1")

# Residual torus equivariance, including the a=0 divisor.
scale = sp.symbols("scale", nonzero=True)
source_scale = {a: scale * a, y: y / scale, z: z / scale**2}
assert sp.cancel(b_forward.subs(source_scale) - b_forward) == 0
assert sp.cancel(c_forward.subs(source_scale) - c_forward) == 0
assert sp.cancel(d_forward.subs(source_scale) - d_forward / scale) == 0
assert sp.cancel(e_forward.subs(source_scale) - e_forward / scale**2) == 0

# Polynomial inverse coordinates.
y_inverse = 2 * b * d - a * e
z_inverse = 2 * d**2 + c * e + 6 * b * d**2 + 3 * b * c * e - sp.Rational(9, 2) * e

assert sp.expand(y_inverse.subs(forward) - y) == 0
assert sp.expand(z_inverse.subs(forward) - z) == 0
print("PASS: inverse after forward is the identity")

coefficient_scale = {
    a: scale * a,
    b: b,
    c: c,
    d: d / scale,
    e: e / scale**2,
}
assert sp.cancel(y_inverse.subs(coefficient_scale) - y_inverse / scale) == 0
assert sp.cancel(z_inverse.subs(coefficient_scale) - z_inverse / scale**2) == 0
print("PASS: polynomial coordinate certificate is equivariant for the residual torus")

# Verify the other composition in the coordinate ring of the complete intersection.
relations = (resultant - 1, middle - 1)
groebner = sp.groebner(relations, e, d, c, b, a, order="lex")
reverse = {
    y: y_inverse,
    z: z_inverse,
}
reconstructed = (
    1 + a * y_inverse,
    1 - sp.Rational(3, 2) * a * y_inverse + a**2 * z_inverse,
    (
        sp.Rational(1, 2) * y
        - a * z
        + sp.Rational(3, 2) * a * y**2
        - a**2 * y * z
    ).subs(reverse),
    (
        -2 * z
        + 4 * y**2
        - 4 * a * y * z
        + 3 * a * y**3
        - 2 * a**2 * y**2 * z
    ).subs(reverse),
)
for got, expected in zip(reconstructed, (b, c, d, e)):
    remainder = groebner.reduce(sp.expand(got - expected))[1]
    assert sp.expand(remainder) == 0
print("PASS: forward after inverse is the identity on the slice")

# The primitive vertical derivation along the projection to (a,b).  It is the
# cross product of the two (c,d,e)-gradients after removing their common
# factor a, and its coefficients are D(c),D(d),D(e).
vertical = {
    a: sp.Integer(0),
    b: sp.Integer(0),
    c: a**2,
    d: -a * b,
    e: -2 * b**2,
}


def D(polynomial):
    return sp.expand(
        sum(sp.diff(polynomial, variable) * value for variable, value in vertical.items())
    )


gradient_middle = sp.Matrix([sp.diff(middle, variable) for variable in (c, d, e)])
gradient_resultant = sp.Matrix(
    [sp.diff(resultant, variable) for variable in (c, d, e)]
)
cross_product = gradient_middle.cross(gradient_resultant)
assert cross_product == a * sp.Matrix([vertical[c], vertical[d], vertical[e]])
assert D(middle - 1) == 0
assert D(resultant - 1) == 0
assert all(D(D(variable)) == 0 for variable in (a, b, c, d, e))
assert D(y_inverse) == 0
assert groebner.reduce(sp.expand(D(z_inverse) - 1))[1] == 0

# The half-normalized generator comes directly from
# Q -> Q + t/2 * L*(a*T-2*b), and the normalization identity b-a*y=1
# integrates term by term to the alternative slice z_tilde=2*z.
def D0(polynomial):
    return sp.cancel(D(polynomial) / 2)


z_bezout = 4 * y_inverse * d + 2 * c * y_inverse**2 - e
assert groebner.reduce(sp.expand(b - a * y_inverse - 1))[1] == 0
assert D0(-e) == b**2
assert sp.expand(D0(4 * y_inverse * d) + 2 * a * b * y_inverse) == 0
assert sp.expand(D0(2 * c * y_inverse**2) - a**2 * y_inverse**2) == 0
assert groebner.reduce(sp.expand(D0(z_bezout) - 1))[1] == 0
assert groebner.reduce(sp.expand(z_bezout - 2 * z_inverse))[1] == 0
print("PASS: primitive vertical LND has invariant y and global polynomial slice z")

# Multiplication on the slice.
g_map = (
    sp.expand((a * c).subs(forward)),
    sp.expand((a * e + b * d).subs(forward)),
    sp.expand((b * e).subs(forward)),
)
g_jacobian = sp.factor(sp.Matrix(g_map).jacobian((a, y, z)).det())
assert g_jacobian == -1
assert sp.cancel(g_map[0].subs(source_scale) - scale * g_map[0]) == 0
assert sp.cancel(g_map[1].subs(source_scale) - g_map[1] / scale) == 0
assert sp.cancel(g_map[2].subs(source_scale) - g_map[2] / scale**2) == 0
print("PASS: normalized multiplication map has determinant -1")

# Linear equivalence with the foundational polynomial.
z1, z2, z3 = sp.symbols("z1 z2 z3")
u = 1 + z1 * z2
foundational = (
    u**3 * z3 + z2**2 * u * (4 + 3 * z1 * z2),
    z2 + 3 * z1 * u**2 * z3 + 3 * z1 * z2**2 * (4 + 3 * z1 * z2),
    2 * z1 - 3 * z1**2 * z2 - z1**3 * z3,
)
g_sub = tuple(
    sp.expand(component.subs({a: z1, y: z2, z: -z3 / 2}))
    for component in g_map
)
assert sp.expand(foundational[0] - g_sub[2]) == 0
assert sp.expand(foundational[1] - 2 * g_sub[1]) == 0
assert sp.expand(foundational[2] - 2 * g_sub[0]) == 0
print("PASS: foundational F = (G3,2G2,2G1) after z -> -z/2")

# Unequal-degree extension: a linear factor and a cubic factor.
T = sp.symbols("T")
a0, a1, b0, b1, b2, b3 = sp.symbols("a0 a1 b0 b1 b2 b3")
A13 = a0 * T + a1
B13 = b0 * T**3 + b1 * T**2 + b2 * T + b3
product13 = sp.Poly(sp.expand(A13 * B13), T).all_coeffs()
resultant13 = sp.factor(sp.resultant(A13, B13, T))
variables13 = (a0, a1, b0, b1, b2, b3)
relative13 = sp.Matrix((a0, a1, -b0, -b1, -b2, -b3))
tangent13 = sp.Matrix((*product13, resultant13)).jacobian(variables13) * relative13
tangent13 = tuple(sp.factor(entry) for entry in tangent13)
assert tangent13[:-1] == (0, 0, 0, 0, 0)
assert sp.factor(tangent13[-1] - 2 * resultant13) == 0
print("PASS: for degrees (1,3), the resultant detects relative scaling with weight 2")

# Equal degrees show the obstruction: the resultant has relative weight zero.
u0, u1, v0, v1 = sp.symbols("u0 u1 v0 v1")
A11 = u0 * T + u1
B11 = v0 * T + v1
product11 = sp.Poly(sp.expand(A11 * B11), T).all_coeffs()
resultant11 = sp.factor(sp.resultant(A11, B11, T))
variables11 = (u0, u1, v0, v1)
relative11 = sp.Matrix((u0, u1, -v0, -v1))
tangent11 = sp.Matrix((*product11, resultant11)).jacobian(variables11) * relative11
assert tuple(sp.factor(entry) for entry in tangent11) == (0, 0, 0, 0)
print("PASS: for equal degrees, multiplication and resultant retain the scaling kernel")

# Consecutive degrees have a unique normalization; a larger gap leaves roots.
lambda23 = m**2 / r
mu23 = r / m**3
assert sp.cancel(lambda23 * mu23 * m - 1) == 0
assert sp.cancel(lambda23**3 * mu23**2 * r - 1) == 0
print("PASS: consecutive degrees (2,3) again have a unique normalized representative")

# Direct finite-field audit of all three cubic hyperplane orbits.  We count
# the projective marked-factor open U_ell, canonically isomorphic to X_ell.
def projective_points(prime, dimension):
    """Canonical representatives of P^dimension(F_prime)."""
    points = []
    for pivot in range(dimension + 1):
        prefix = (0,) * pivot + (1,)
        for tail in __import__("itertools").product(range(prime), repeat=dimension - pivot):
            points.append(prefix + tail)
    return points


def mod_resultant_linear_quadratic(linear, quadratic, prime):
    aa, bb = linear
    cc, dd, ee = quadratic
    return (aa * aa * ee - aa * bb * dd + bb * bb * cc) % prime


def product_coefficients(linear, quadratic, prime):
    aa, bb = linear
    cc, dd, ee = quadratic
    return (
        aa * cc % prime,
        (aa * dd + bb * cc) % prime,
        (aa * ee + bb * dd) % prime,
        bb * ee % prime,
    )


# Representatives ell=m1-m2, ell=m1, ell=m0 have contact types
# (1,1,1), (2,1), and (3), respectively.
orbit_functionals = {
    "111": lambda coeffs, p: (coeffs[1] - coeffs[2]) % p,
    "21": lambda coeffs, p: coeffs[1] % p,
    "3": lambda coeffs, p: coeffs[0] % p,
}
# The chosen integral representative m1-m2 has bad reduction at 3 because its
# restriction to the twisted cubic acquires the scalar factor 3.  Use good
# primes for the three-orbit comparison.
for prime in (2, 5, 7, 11):
    linears = projective_points(prime, 1)
    quadratics = projective_points(prime, 2)
    counts = {}
    for name, functional in orbit_functionals.items():
        counts[name] = sum(
            mod_resultant_linear_quadratic(linear, quadratic, prime) != 0
            and functional(product_coefficients(linear, quadratic, prime), prime) != 0
            for linear in linears
            for quadratic in quadratics
        )
    assert counts == {
        "111": prime**3 - prime,
        "21": prime**3,
        "3": prime**3 - prime**2,
    }
print("PASS: all cubic hyperplane orbits have the stated point-count classes")
