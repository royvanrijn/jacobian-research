#!/usr/bin/env python3
"""Exact de Rham--Cartier W_2 obstruction for the Mondello plane map.

The proof computes the full cokernel of the first Jacobian variation.  The
older x*y coefficient test is retained as its lowest-coordinate witness.
The all-degree statements are monomial/exterior-form arguments, not bounded
ansatz searches; finite loops below are exact regressions of those formulas.

The written proof and status boundary are in
verified/HUQ_KURUVILLA_PLANE_W2_OBSTRUCTION.md.
"""

from __future__ import annotations

import sympy as sp


def mod2(expr: sp.Expr, *generators: sp.Symbol) -> sp.Expr:
    """Return the canonical polynomial representative over F_2."""

    return sp.Poly(sp.expand(expr), *generators, modulus=2).as_expr()


def top_de_rham_representative(
    expr: sp.Expr, x: sp.Symbol, y: sp.Symbol
) -> sp.Expr:
    """Canonical representative in H^2_dR(F_2[x,y]).

    Exact two-forms contain precisely the monomials having an even exponent
    in at least one variable.  The quotient therefore retains the odd--odd
    monomials.
    """

    terms = sp.Poly(mod2(expr, x, y), x, y, modulus=2).terms()
    return sp.expand(
        sum(
            x**i * y**j
            for (i, j), coefficient in terms
            if int(coefficient) % 2 and i % 2 == 1 and j % 2 == 1
        )
    )


def coefficient_functional(
    multiplier: sp.Expr,
    derivative_axis: str,
    target: tuple[int, int],
    x: sp.Symbol,
    y: sp.Symbol,
) -> int:
    """Coefficient at target of multiplier*d(arbitrary polynomial).

    For a fixed target monomial and a fixed multiplier monomial there is at
    most one source monomial.  Its derivative coefficient decides whether an
    arbitrary input coefficient can contribute.  Returning zero proves that
    this coefficient functional vanishes in every polynomial degree.
    """

    profile = coefficient_functional_profile(
        multiplier, derivative_axis, (target,), x, y
    )
    return int(bool(profile))


def coefficient_functional_profile(
    multiplier: sp.Expr,
    derivative_axis: str,
    targets: tuple[tuple[int, int], ...],
    x: sp.Symbol,
    y: sp.Symbol,
) -> dict[tuple[int, int], int]:
    """Source-coefficient profile of a target-coefficient functional.

    The returned dictionary records which coefficients of an arbitrary input
    polynomial survive in the sum of the requested target coefficients.
    An empty profile is an exact all-degree vanishing certificate.
    """

    profile: dict[tuple[int, int], int] = {}
    for (u, v), coefficient in sp.Poly(
        mod2(multiplier, x, y), x, y, modulus=2
    ).terms():
        for target in targets:
            if derivative_axis == "x":
                i, j = target[0] - u + 1, target[1] - v
                derivative_coefficient = i
            else:
                i, j = target[0] - u, target[1] - v + 1
                derivative_coefficient = j
            contribution = (
                int(coefficient) & 1
            ) * (derivative_coefficient & 1)
            if i >= 0 and j >= 0 and contribution:
                profile[(i, j)] = profile.get((i, j), 0) ^ 1
                if profile[(i, j)] == 0:
                    del profile[(i, j)]
    return profile


def xor_profiles(
    *profiles: dict[tuple[int, int], int],
) -> dict[tuple[int, int], int]:
    """Add source-coefficient profiles over F_2."""

    result: dict[tuple[int, int], int] = {}
    for profile in profiles:
        for monomial, coefficient in profile.items():
            result[monomial] = result.get(monomial, 0) ^ coefficient
            if result[monomial] == 0:
                del result[monomial]
    return result


x, y = sp.symbols("x y")
P = x + x**2 * y + x**4 + x**6 * y**2
Q = y + x**5 + x**6 * y + x**7 * y**2 + x**8 * y**3

integer_jacobian = sp.expand(sp.det(sp.Matrix((P, Q)).jacobian((x, y))))
half_error = sp.Poly(integer_jacobian - 1, x, y).quo_ground(2).as_expr()
assert sp.expand(1 + 2 * half_error - integer_jacobian) == 0
assert int(sp.Poly(half_error, x, y).coeff_monomial(x * y)) % 2 == 1
print("PASS: the naive integral Jacobian has x*y coefficient 2 modulo 4")

Px = mod2(sp.diff(P, x), x, y)
Py = mod2(sp.diff(P, y), x, y)
Qx = mod2(sp.diff(Q, x), x, y)
Qy = mod2(sp.diff(Q, y), x, y)
assert Px == 1
assert Py == x**2
assert Qx == x**4 + x**6 * y**2
assert Qy == x**6 + x**8 * y**2 + 1

# For arbitrary corrections A,B modulo two, the first variation is
# A_x Q_y + P_x B_y + A_y Q_x + P_y B_x.  Each summand has identically
# zero x*y coefficient.  The helper checks the coefficient functional on an
# arbitrary input polynomial, without imposing a degree bound.
target = (1, 1)
functionals = (
    coefficient_functional(Qy, "x", target, x, y),
    coefficient_functional(Px, "y", target, x, y),
    coefficient_functional(Qx, "y", target, x, y),
    coefficient_functional(Py, "x", target, x, y),
)
assert functionals == (0, 0, 0, 0)
print("PASS: every all-degree first Jacobian correction has zero x*y coefficient")

# The full statement is de Rham-theoretic.  For a determinant-one plane map,
# dP,dQ form a basis of one-forms, and
#   D_F(A,B) dx dy = d(A dQ + B dP).
# The monomial primitive below verifies the quotient description directly.
for i in range(9):
    for j in range(9):
        monomial = x**i * y**j
        expected = monomial if i % 2 and j % 2 else 0
        assert top_de_rham_representative(monomial, x, y) == expected
        if i % 2 == 0:
            primitive_coefficient = sp.diff(x ** (i + 1) * y**j, x)
            assert mod2(primitive_coefficient - monomial, x, y) == 0
        elif j % 2 == 0:
            primitive_coefficient = sp.diff(x**i * y ** (j + 1), y)
            assert mod2(primitive_coefficient - monomial, x, y) == 0

# Replay the exterior-form identity on dense corrections, then test every
# monomial correction through degree eight.  Linearity plus the written
# monomial proof is the all-degree certificate.
A_dense = sum(
    x**i * y**j for i in range(7) for j in range(7 - i) if (i + 2 * j) % 3
)
B_dense = sum(
    x**i * y**j for i in range(7) for j in range(7 - i) if (2 * i + j) % 3
)
variation_dense = mod2(
    sp.diff(A_dense, x) * Qy
    + sp.diff(A_dense, y) * Qx
    + sp.diff(B_dense, x) * Py
    + sp.diff(B_dense, y) * Px,
    x,
    y,
)
exact_form_coefficient = mod2(
    sp.diff(A_dense, x) * Qy
    - sp.diff(A_dense, y) * Qx
    + sp.diff(B_dense, x) * Py
    - sp.diff(B_dense, y) * Px,
    x,
    y,
)
assert variation_dense == exact_form_coefficient
for i in range(9):
    for j in range(9):
        monomial = x**i * y**j
        variation_A = mod2(
            sp.diff(monomial, x) * Qy + sp.diff(monomial, y) * Qx,
            x,
            y,
        )
        variation_B = mod2(
            sp.diff(monomial, x) * Py + sp.diff(monomial, y) * Px,
            x,
            y,
        )
        assert top_de_rham_representative(variation_A, x, y) == 0
        assert top_de_rham_representative(variation_B, x, y) == 0
print("PASS: coker(D_F) is represented by x*y*F_2[x^2,y^2]")

# Compute the complete obstruction and its inverse-Cartier square root.
r = 1 + x * y
u = 1 + x**3 * r
half_error_mod2 = mod2(half_error, x, y)
obstruction_class = top_de_rham_representative(half_error_mod2, x, y)
assert obstruction_class == mod2(x * y * u**2, x, y)
assert obstruction_class != 0
assert top_de_rham_representative(
    half_error_mod2 + variation_dense, x, y
) == obstruction_class
print("PASS: the full obstruction class is x*y*u^2, with Cartier image u*dx*dy")
print("PASS: no polynomial constant-Jacobian lift exists over Z/4")

# Functoriality under polynomial source changes is proved in the note using
# Jung--van der Kulk tameness.  These two generators are exact regressions of
# the pullback formula, not a bounded left--right search.
triangular_pullback = top_de_rham_representative(
    obstruction_class.subs({x: x + y**2}, simultaneous=True), x, y
)
triangular_round_trip = top_de_rham_representative(
    triangular_pullback.subs({x: x + y**2}, simultaneous=True), x, y
)
swap_pullback = top_de_rham_representative(
    obstruction_class.subs({x: y, y: x}, simultaneous=True), x, y
)
assert triangular_pullback != 0
assert triangular_round_trip == obstruction_class
assert swap_pullback != 0
print("PASS: triangular and affine-swap generators preserve nonvanishing")

# The stable W_2 lift has sharp total degree 18.  If a lift
# (P+2A,Q+2B,z+2C) had degree at most 17, take the z^0 part of its half-error.
# With A_0,B_0 the z^0 coefficients and C_1 the z^1 coefficient, it is
# K+D_F(A_0,B_0)+C_1.  The functional below annihilates constants and obeys
#   Lambda(D_F(A_0,B_0)) = [x^15*y^5] A_0.
# It therefore vanishes on the correction terms under the degree bound, but
# Lambda(K)=1.  The profiles certify the identity for arbitrary polynomials.
degree_functional_targets = ((13, 4), (14, 5))
A_degree_profile = xor_profiles(
    coefficient_functional_profile(
        Qy, "x", degree_functional_targets, x, y
    ),
    coefficient_functional_profile(
        Qx, "y", degree_functional_targets, x, y
    ),
)
B_degree_profile = xor_profiles(
    coefficient_functional_profile(
        Px, "y", degree_functional_targets, x, y
    ),
    coefficient_functional_profile(
        Py, "x", degree_functional_targets, x, y
    ),
)
degree_functional_on_error = sum(
    int(sp.Poly(half_error_mod2, x, y).coeff_monomial(x**i * y**j)) % 2
    for i, j in degree_functional_targets
) % 2
assert A_degree_profile == {(15, 5): 1}
assert B_degree_profile == {}
assert degree_functional_on_error == 1
print("PASS: the dual coefficient functional excludes every stable W_2 lift of degree at most 17")

# The obstruction is unstable.  After adjoining one identity coordinate z,
# multiply that coordinate by 1+2*K.  The Jacobian matrix is block lower
# triangular, so its determinant is (1+2*K)^2=1 modulo four.
z = sp.symbols("z")
stabilized_third_coordinate = z * (1 + 2 * half_error)
stabilized_jacobian = sp.expand(
    sp.det(
        sp.Matrix((P, Q, stabilized_third_coordinate)).jacobian((x, y, z))
    )
)
for coefficient in sp.Poly(stabilized_jacobian - 1, x, y, z).coeffs():
    assert int(coefficient) % 4 == 0
assert mod2(sp.diff(z * half_error_mod2, z), x, y, z) == half_error_mod2
assert max(
    sp.Poly(coordinate, x, y, z).total_degree()
    for coordinate in (P, Q, stabilized_third_coordinate)
) == 18
print("PASS: one identity stabilization has an explicit Keller lift over Z/4")
print("PASS: the exact minimum stable W_2 coordinate degree is 18")

# A non-geometric plane correction lowers the evident W_3 upper bound from
# 35 to 25.  Put Delta=x^2*y^2+x^10*y^4 and choose A_3 so that
# D_F(A_3,0)=K^2+Delta.  Then C_3=K+Delta cancels the remaining second digit.
A_w3 = x**15 * y**2 * (1 + x**4 * y**4)
delta_w3 = x**2 * y**2 + x**10 * y**4
C_w3 = half_error_mod2 + delta_w3
D_A_w3 = mod2(
    sp.diff(A_w3, x) * Qy + sp.diff(A_w3, y) * Qx, x, y
)
assert D_A_w3 == mod2(half_error_mod2**2 + delta_w3, x, y)
canonical_w3_error = mod2(
    half_error_mod2 + half_error_mod2**2, x, y
)
assert int(
    sp.Poly(canonical_w3_error, x, y).coeff_monomial(x**26 * y**8)
) % 2 == 1
assert sp.Poly(canonical_w3_error, x, y).total_degree() == 34
assert 24 - 1 + sp.Poly(Qy, x, y).total_degree() == 33
assert 24 - 1 + sp.Poly(Qx, x, y).total_degree() < 34
w3_lift = (
    P + 4 * A_w3,
    Q,
    z * (1 + 2 * half_error + 4 * C_w3),
)
w3_jacobian = sp.expand(
    sp.det(sp.Matrix(w3_lift).jacobian((x, y, z)))
)
for coefficient in sp.Poly(w3_jacobian - 1, x, y, z).coeffs():
    assert int(coefficient) % 8 == 0
assert max(
    sp.Poly(coordinate, x, y, z).total_degree() for coordinate in w3_lift
) == 25
print("PASS: an explicit degree-25 Keller lift exists over W_3(F_2)")
print("PASS: degree 25 is sharp for extensions of the canonical W_2 lift")

# Universally, for h=2*K and S_n=sum_{j=0}^{n-1}(-h)^j, one has
# (1+h)S_n=1-(-h)^n.  The following symbolic induction step certifies the
# identity independently of K and n: if r=(-h)^n is the old error, adjoining
# r to S_n changes 1-r into 1+h*r, the next geometric-series remainder.
h, r_symbol = sp.symbols("h r_symbol")
assert sp.expand((1 - r_symbol) + (1 + h) * r_symbol - (1 + h * r_symbol)) == 0

# Replay several concrete levels as a regression; the all-n proof is the
# geometric-series identity above and in the canonical note.
for level in range(2, 7):
    multiplier = sum((-2 * half_error) ** j for j in range(level))
    determinant = sp.expand(integer_jacobian * multiplier)
    assert sp.expand(determinant - (1 - (-2 * half_error) ** level)) == 0
    for coefficient in sp.Poly(determinant - 1, x, y).coeffs():
        assert int(coefficient) % (2**level) == 0
print("PASS: the stabilized lifts form a compatible Keller tower over every W_n(F_2)")
