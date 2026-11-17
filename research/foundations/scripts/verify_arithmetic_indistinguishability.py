#!/usr/bin/env python3
"""Exact checks for arithmetic indistinguishability of weighted seeds.

The first half verifies the degree-five involution that identifies the full
inverse pencils while moving the affine root mark.  The second half verifies
the classical degree-seven Davenport factorization and the point/line
Gassmann triple for the Fano plane.
"""

from itertools import product

import sympy as sp


# ---------------------------------------------------------------------------
# The degree-five marked/unmarked separation.

w, lam = sp.symbols("w lambda")
mu = sp.Rational(4, 5) - lam
beta = sp.Rational(4, 5)


def quintic(parameter):
    return sp.expand(
        w**2
        * (w - 1)
        * (3 * w**2 - (5 * parameter + 1) * w + 3 * parameter)
        / 60
    )


H = quintic(lam)
G = quintic(mu)
u = sp.factor(sp.diff(G.subs(w, beta - w) + H, w))
v = sp.factor((G.subs(w, beta - w) + H).subs(w, 0))

assert sp.factor(u - 2 * (5 * lam - 2) / 1875) == 0
assert sp.factor(v + 4 * (25 * lam + 8) / 46875) == 0
assert sp.factor(G.subs(w, beta - w) + H - u * w - v) == 0

s, t = sp.symbols("s t")
s_prime = s - u
t_prime = -t - v + beta * (s - u)
assert sp.factor(
    G.subs(w, beta - w) - s_prime * (beta - w) + t_prime
    + (H - s * w + t)
) == 0

# The direct C=0 fiber is one rational sheet plus D*x^2=1.  These scalings
# identify D for the two parameters, so they identify every boundary-fiber
# zeta function, not only the number of rational points.
A, B = sp.symbols("A B")
c_lam = (lam - 1) / 30
c_mu = (mu - 1) / 30
k_lam = sp.factor((-lam / 20) / c_lam)
k_mu = sp.factor((-mu / 20) / c_mu)
D_lam = B**2 / c_lam**2 - 4 * k_lam * A
D_mu = B**2 / c_mu**2 - 4 * k_mu * A
assert sp.factor(
    D_mu.subs({A: k_lam * A / k_mu, B: c_mu * B / c_lam}, simultaneous=True)
    - D_lam
) == 0


# ---------------------------------------------------------------------------
# Cassou-Nogues--Couveignes' degree-seven Davenport pair.

y, z, a, T = sp.symbols("y z a T")
minimal_a = a**2 + a + 2

g = (
    sp.Rational(1, 7) * y**7
    + (1 + a) * T * y**5
    + (1 + a) * T * y**4
    - (3 - 2 * a) * T**2 * y**3
    - 2 * (1 - 2 * a) * T**2 * y**2
    - sp.Rational(1, 28) * (5 + 3 * a) * (28 * T - 2 - 11 * a) * T**2 * y
    - (1 + a) * T**3
)


def conjugate_a(expression):
    """Apply a -> -1-a and reduce modulo a^2+a+2."""
    substituted = sp.expand(expression.subs(a, -1 - a))
    domain = sp.QQ.frac_field(T, y, z)
    return sp.rem(
        sp.Poly(substituted, a, domain=domain),
        sp.Poly(minimal_a, a, domain=domain),
    ).as_expr()


h = conjugate_a(g.subs(y, z))
factor = (
    y**3
    + a * y**2 * z
    + y**2 * z
    + a * y * z**2
    - z**3
    + 5 * T * y
    + 3 * a * T * y
    - 2 * T * z
    + 3 * a * T * z
    + 2 * a * T
    + T
)

groebner = sp.groebner(
    [minimal_a, factor], y, z, a, order="grlex", domain=sp.QQ.frac_field(T)
)
assert groebner.reduce(sp.expand(g - h))[1] == 0

# The Hessian divisors of g and h are not affinely equivalent when T != 0.
# Both Hessians have zero x^4 coefficient, so an affine equivalence has no
# translation.  The x^3 and x^2 coefficients would then force alpha^2=r and
# alpha^3=r, with r=-(1+a)/a.  Hence alpha=r=1, contrary to 2a+1 != 0.
g_second = sp.diff(g, y, 2)
h_second = sp.diff(h, z, 2)
assert sp.Poly(g_second, y).coeff_monomial(y**5) == 6
assert sp.Poly(h_second, z).coeff_monomial(z**5) == 6
assert sp.Poly(g_second, y).coeff_monomial(y**4) == 0
assert sp.Poly(h_second, z).coeff_monomial(z**4) == 0
assert sp.expand(
    sp.Poly(g_second, y).coeff_monomial(y**3) - 20 * (1 + a) * T
) == 0
assert sp.expand(
    sp.Poly(h_second, z).coeff_monomial(z**3) + 20 * a * T
) == 0
assert sp.expand(
    sp.Poly(g_second, y).coeff_monomial(y**2) - 12 * (1 + a) * T
) == 0
assert sp.expand(
    sp.Poly(h_second, z).coeff_monomial(z**2) + 12 * a * T
) == 0
assert sp.resultant(minimal_a, 2 * a + 1, a) != 0


# ---------------------------------------------------------------------------
# GL(3,2) on Fano points and lines.


def mat_vec(matrix, vector):
    bits = ((vector >> 0) & 1, (vector >> 1) & 1, (vector >> 2) & 1)
    output = []
    for row in range(3):
        output.append(
            sum(matrix[3 * row + col] * bits[col] for col in range(3)) % 2
        )
    return output[0] | (output[1] << 1) | (output[2] << 2)


def mat_mul(left, right):
    return tuple(
        sum(left[3 * row + k] * right[3 * k + col] for k in range(3)) % 2
        for row in range(3)
        for col in range(3)
    )


def transpose(matrix):
    return tuple(matrix[3 * col + row] for row in range(3) for col in range(3))


identity = (1, 0, 0, 0, 1, 0, 0, 0, 1)
matrices = []
for entries in product((0, 1), repeat=9):
    images = {mat_vec(entries, vector) for vector in range(1, 8)}
    if len(images) == 7 and 0 not in images:
        matrices.append(entries)
assert len(matrices) == 168

inverse = {}
for matrix in matrices:
    inverse[matrix] = next(
        candidate for candidate in matrices if mat_mul(matrix, candidate) == identity
    )


def point_permutation(matrix):
    return tuple(mat_vec(matrix, vector) - 1 for vector in range(1, 8))


def line_permutation(matrix):
    dual = transpose(inverse[matrix])
    return tuple(mat_vec(dual, covector) - 1 for covector in range(1, 8))


def fixed_points(permutation):
    return sum(index == image for index, image in enumerate(permutation))


for matrix in matrices:
    assert fixed_points(point_permutation(matrix)) == fixed_points(
        line_permutation(matrix)
    )

point_stabilizer = {m for m in matrices if point_permutation(m)[0] == 0}
line_stabilizer = {m for m in matrices if line_permutation(m)[0] == 0}
assert len(point_stabilizer) == len(line_stabilizer) == 24
for matrix in matrices:
    conjugate = {
        mat_mul(mat_mul(matrix, member), inverse[matrix])
        for member in point_stabilizer
    }
    assert conjugate != line_stabilizer


def permutation(*cycles):
    result = list(range(7))
    for cycle in cycles:
        for index, entry in enumerate(cycle):
            result[entry - 1] = cycle[(index + 1) % len(cycle)] - 1
    return tuple(result)


def compose(left, right):
    return tuple(left[right[index]] for index in range(7))


sigma_1 = permutation((1, 2), (3, 6))
sigma_half = permutation((2, 3), (4, 5))
sigma_0 = permutation((3, 4), (6, 7))

generated = {tuple(range(7))}
frontier = list(generated)
while frontier:
    current = frontier.pop()
    for generator in (sigma_1, sigma_half, sigma_0):
        candidate = compose(generator, current)
        if candidate not in generated:
            generated.add(candidate)
            frontier.append(candidate)
assert len(generated) == 168
sigma_infinity = compose(compose(sigma_1, sigma_half), sigma_0)
orbit = []
entry = 0
for _ in range(7):
    orbit.append(entry)
    entry = sigma_infinity[entry]
assert len(set(orbit)) == 7 and entry == 0

print("PASS: the quintic involution identifies the inverse pencils exactly")
print("PASS: the quintic C=0 fiber equations agree after boundary relabeling")
print("PASS: the explicit degree-seven Davenport factor divides g(y)-h(z)")
print("PASS: the Davenport Hessian divisors are not affinely equivalent for T != 0")
print("PASS: Fano point/line stabilizers are nonconjugate and Gassmann equivalent")
print("PASS: the displayed double-transposition branch cycles generate GL(3,2)")
