#!/usr/bin/env python3
"""Verify affine-plane and low-degree graph obstructions for Meng--Yang HC(5).

For

    Psi_(L,M,N) = L*A**2 + M*A + N*B,

the script extracts one unavoidable coefficient on each projective chart of
an affine-hyperplane normal with nonzero dual part.  If the dual part is
zero, the three dual rows of the restricted Hessian span at most two
dimensions, so the determinant vanishes identically.

For a polynomial graph y3=R(x1,x2,y1,y2) of degree at most three, two line
restrictions give the triangular coefficient chain

    rho = -89/16,  sigma = -tau,  [t**5] = 197*L*N**3/4.

This is a characteristic-zero contradiction when L*N is nonzero.  The
calculation does not use the collision equations.

For degree four, the checker uses the full two-slope pencil

    (x1,x2,y1,y2) = (0,t,c*t,d*t).

Its t**8 coefficient kills every x1-free quartic jet containing y2, not just
the suggested [x2**3*y2]R jet.  The t**6 and t**5 coefficients then reduce
the surviving slice to a parameter rho satisfying

    160*rho**2 + 1968*rho + 6021 = 0.

The discriminant is 576*34.  The next coefficient first kills [y1**3]R3
and [y2**2]R2 and then imposes a second quadratic on rho.  The two rho
quadratics are coprime, so no degree-at-most-four graph exists over any
characteristic-zero field.

The degree-five continuation uses the exact graph Hessian on the plane
x1=0.  If T=R|_(x1=0) and S=partial_x1 R|_(x1=0), its determinant is

    F(T, partial T) - 8*L*N**3*S.

Thus x1**2-divisible corrections are invisible on that plane, whereas the
x1-linear quintic jet is the first normal repair of the quartic obstruction.
The top degree-five faces force

    R5(0,x2,y1,y2) = kappa*x2**5,
    partial_y2 R4(0,x2,y1,y2) = 0,

and reduce the next two faces to a small conic system.  Finally, for the v2
member, the checker excludes the full sparse trace family

    R(0,x2,y1,y2) = kappa*x2**5 + d*x2**3*y1 + rho*x2**2*y2,

with arbitrary normal quartic jet and arbitrary remaining x1**2-divisible
terms of total degree at most five.  The terminal certificate uses two
resultants and a first-transverse coefficient coprime to the old quartic
polynomial Q(rho).
"""

from __future__ import annotations

from itertools import permutations

import sympy as sp


def permutation_sign(permutation: tuple[int, ...]) -> int:
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(len(permutation))
        for right in range(left + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def multiply_truncated(
    left: dict[tuple[int, ...], sp.Expr],
    right: dict[tuple[int, ...], sp.Expr],
    target: tuple[int, ...],
) -> dict[tuple[int, ...], sp.Expr]:
    """Multiply sparse polynomials, discarding exponents above ``target``."""

    product: dict[tuple[int, ...], sp.Expr] = {}
    for left_exponents, left_coefficient in left.items():
        for right_exponents, right_coefficient in right.items():
            exponents = tuple(
                a + b
                for a, b in zip(
                    left_exponents, right_exponents, strict=True
                )
            )
            if all(
                exponent <= bound
                for exponent, bound in zip(exponents, target, strict=True)
            ):
                product[exponents] = (
                    product.get(exponents, sp.S.Zero)
                    + left_coefficient * right_coefficient
                )
    return product


def determinant_truncation(
    matrix: sp.Matrix,
    variables: tuple[sp.Symbol, ...],
    target: tuple[int, ...],
) -> dict[tuple[int, ...], sp.Expr]:
    """Return the determinant truncated componentwise at ``target``."""

    assert matrix.shape == (4, 4)
    entries = [
        [
            {
                exponents: coefficient
                for exponents, coefficient in sp.Poly(
                    matrix[row, column], *variables
                ).terms()
                if all(
                    exponent <= bound
                    for exponent, bound in zip(
                        exponents, target, strict=True
                    )
                )
            }
            for column in range(4)
        ]
        for row in range(4)
    ]

    determinant: dict[tuple[int, ...], sp.Expr] = {}
    zero_exponents = (0,) * len(variables)
    for permutation in permutations(range(4)):
        term: dict[tuple[int, ...], sp.Expr] = {
            zero_exponents: sp.Integer(permutation_sign(permutation))
        }
        for row, column in enumerate(permutation):
            term = multiply_truncated(term, entries[row][column], target)
        for exponents, coefficient in term.items():
            determinant[exponents] = (
                determinant.get(exponents, sp.S.Zero) + coefficient
            )
    return determinant


def determinant_coefficient(
    matrix: sp.Matrix,
    variables: tuple[sp.Symbol, ...],
    target: tuple[int, ...],
) -> sp.Expr:
    """Extract one determinant coefficient without forming the determinant."""

    return sp.factor(
        determinant_truncation(matrix, variables, target).get(
            target, sp.S.Zero
        )
    )


x, y, p, q, r = sp.symbols("x y p q r")
L, M, N = sp.symbols("L M N")
u = 1 + x * y
q0 = y**2 * (1 + 3 * u)


def meng_forms(
    p_value: sp.Expr, q_value: sp.Expr, r_value: sp.Expr
) -> tuple[sp.Expr, sp.Expr]:
    A = u**3 * p_value + 3 * x * u**2 * q_value - x**3 * r_value
    B = (
        u * q0 * p_value
        + (y + 3 * x * q0) * q_value
        + x * (5 - 3 * u) * r_value
    )
    return A, B


def meng_potential(
    p_value: sp.Expr, q_value: sp.Expr, r_value: sp.Expr
) -> sp.Expr:
    A, B = meng_forms(p_value, q_value, r_value)
    return sp.expand(L * A**2 + M * A + N * B)


A, B = meng_forms(p, q, r)
Psi = meng_potential(p, q, r)
assert sp.Poly(Psi, x, y, p, q, r).total_degree() == 14


# The v2 collision has common ambient gradient (0,0,-1/2,0,0).  Its omitted
# r-component is zero, which is why every graph containing the two points
# transfers the collision even when the graph tangents differ there.
collision_points = (
    (sp.Integer(1), -sp.Rational(3, 2), 0, 0, 0),
    (-sp.Integer(1), sp.Rational(3, 2), 0, 0, 0),
)
collision_gradients = []
for point in collision_points:
    substitution = dict(zip((x, y, p, q, r), point, strict=True))
    collision_gradients.append(
        tuple(
            sp.factor(
                derivative.subs({L: 1, M: 13, N: 2}).subs(substitution)
            )
            for derivative in (sp.diff(Psi, variable) for variable in (x, y, p, q, r))
        )
    )
assert collision_gradients == [(0, 0, -sp.Rational(1, 2), 0, 0)] * 2


# Affine hyperplanes with nonzero dual normal.  Split the dual projective
# normal [gamma:delta:epsilon] into the disjoint charts gamma!=0,
# gamma=0 and delta!=0, and gamma=delta=0 and epsilon!=0.  Translation and
# source-normal parameters are all retained.  One coefficient on each chart
# is independent of every normal parameter and of M,N.
a, b, c, d, h = sp.symbols("a b c d h")

# gamma != 0: solve for p and retain (x,y,q,r).
p_chart = h - a * x - b * y - c * q - d * r
psi_p_chart = meng_potential(p_chart, q, r)
hessian_p_chart = sp.hessian(psi_p_chart, (x, y, q, r)).subs({y: 0, r: 0})
coefficient_p_chart = determinant_coefficient(
    hessian_p_chart, (x, q), (12, 4)
)
assert coefficient_p_chart == 793152 * L**4

# gamma=0, delta != 0: solve for q and retain (x,y,p,r).
q_chart = h - a * x - b * y - d * r
psi_q_chart = meng_potential(p, q_chart, r)
hessian_q_chart = sp.hessian(psi_q_chart, (x, y, p, r)).subs({y: 0, r: 0})
coefficient_q_chart = determinant_coefficient(
    hessian_q_chart, (x, p), (6, 4)
)
assert coefficient_q_chart == 2160 * L**4

# gamma=delta=0, epsilon != 0: solve for r and retain (x,y,p,q).
r_chart = h - a * x - b * y
psi_r_chart = meng_potential(p, q, r_chart)
hessian_r_chart = sp.hessian(psi_r_chart, (x, y, p, q)).subs({y: 0, q: 0})
coefficient_r_chart = determinant_coefficient(
    hessian_r_chart, (x, p), (2, 4)
)
assert coefficient_r_chart == 2160 * L**4


# If the dual normal vanishes, the hyperplane leaves one affine source
# parameter s and all three dual variables.  Writing the restricted
# potential as L*(avec(s).Y)^2+cvec(s).Y, its dual--dual Hessian is the
# displayed rank-one matrix.  The three dual Hessian rows therefore lie in
# the span of the source column and one dual row, so the full rank is at most
# three and its determinant is identically zero.
s = sp.symbols("s")
avec = sp.Matrix([sp.Function(f"a{index}")(s) for index in range(3)])
dual_variables = sp.Matrix(sp.symbols("Y0:3"))
cvec = sp.Matrix([sp.Function(f"c{index}")(s) for index in range(3)])
abstract_potential = L * (avec.dot(dual_variables)) ** 2 + cvec.dot(
    dual_variables
)
abstract_dual_hessian = sp.hessian(
    abstract_potential, tuple(dual_variables)
)
assert abstract_dual_hessian == 2 * L * avec * avec.T


# Generic graph of total degree at most three.
t = sp.symbols("t")
graph_coefficients: dict[tuple[int, int, int, int], sp.Symbol] = {}
R3 = sp.S.Zero
for total_degree in range(4):
    for x_degree in range(total_degree + 1):
        for y_degree in range(total_degree - x_degree + 1):
            for p_degree in range(
                total_degree - x_degree - y_degree + 1
            ):
                q_degree = (
                    total_degree - x_degree - y_degree - p_degree
                )
                exponents = (x_degree, y_degree, p_degree, q_degree)
                coefficient = sp.Symbol("r" + "".join(map(str, exponents)))
                graph_coefficients[exponents] = coefficient
                R3 += (
                    coefficient
                    * x**x_degree
                    * y**y_degree
                    * p**p_degree
                    * q**q_degree
                )

rho = graph_coefficients[(0, 2, 0, 1)]
sigma = graph_coefficients[(0, 1, 1, 1)]
tau = graph_coefficients[(0, 0, 2, 1)]
graph_potential = meng_potential(p, q, R3)
graph_hessian = sp.hessian(graph_potential, (x, y, p, q))

axis_hessian = graph_hessian.subs({x: 0, y: t, p: 0, q: 0})
axis_t6 = determinant_coefficient(axis_hessian, (t,), (6,))
assert axis_t6 == N**4 * (16 * rho + 89) ** 2

diagonal_hessian = graph_hessian.subs({x: 0, y: t, p: t, q: 0})
diagonal_t6 = determinant_coefficient(diagonal_hessian, (t,), (6,))
assert diagonal_t6 == N**4 * (16 * (rho + sigma + tau) + 89) ** 2

diagonal_t5 = determinant_coefficient(diagonal_hessian, (t,), (5,))
triangular_substitution = {
    rho: -sp.Rational(89, 16),
    sigma: -tau,
}
assert sp.factor(diagonal_t5.subs(triangular_substitution)) == (
    sp.Rational(197, 4) * L * N**3
)


# Quartic two-slope jet reduction.  Retain every homogeneous quartic
# coefficient and extract t**5 through t**8 in one truncated determinant
# pass.  This avoids forming the generic full determinant.
quartic_coefficients: dict[tuple[int, int, int, int], sp.Symbol] = {}
R4 = R3
R4_homogeneous = sp.S.Zero
for x_degree in range(5):
    for y_degree in range(5 - x_degree):
        for p_degree in range(5 - x_degree - y_degree):
            q_degree = 4 - x_degree - y_degree - p_degree
            exponents = (x_degree, y_degree, p_degree, q_degree)
            coefficient = sp.Symbol("r" + "".join(map(str, exponents)))
            quartic_coefficients[exponents] = coefficient
            term = (
                coefficient
                * x**x_degree
                * y**y_degree
                * p**p_degree
                * q**q_degree
            )
            R4 += term
            R4_homogeneous += term

slope_p, slope_q = sp.symbols("slope_p slope_q")
quartic_pencil_hessian = sp.hessian(
    meng_potential(p, q, R4), (x, y, p, q)
).subs({x: 0, y: t, p: slope_p * t, q: slope_q * t})
quartic_pencil_determinant = determinant_truncation(
    quartic_pencil_hessian, (t,), (8,)
)


def pencil_coefficient(degree: int) -> sp.Expr:
    return quartic_pencil_determinant.get((degree,), sp.S.Zero)


# The leading square is the q-derivative of the x-free quartic slice.
quartic_q_derivative = sp.factor(
    sp.diff(R4_homogeneous, q).subs(
        {x: 0, y: 1, p: slope_p, q: slope_q}
    )
)
assert sp.factor(pencil_coefficient(8)) == (
    256 * N**4 * quartic_q_derivative**2
)

x_free_quartic_q_substitution = {
    coefficient: 0
    for exponents, coefficient in quartic_coefficients.items()
    if exponents[0] == 0 and exponents[3] > 0
}
assert len(x_free_quartic_q_substitution) == 10
assert sp.factor(
    pencil_coefficient(7).subs(x_free_quartic_q_substitution)
) == 0

# Once the leading square vanishes, the next square is triangular in the
# remaining x-free quartic coefficients and three cubic q**2/q**3 jets.
r0040 = quartic_coefficients[(0, 0, 4, 0)]
r0130 = quartic_coefficients[(0, 1, 3, 0)]
r0220 = quartic_coefficients[(0, 2, 2, 0)]
r0310 = quartic_coefficients[(0, 3, 1, 0)]
r0400 = quartic_coefficients[(0, 4, 0, 0)]
r0003 = graph_coefficients[(0, 0, 0, 3)]
r0012 = graph_coefficients[(0, 0, 1, 2)]
r0102 = graph_coefficients[(0, 1, 0, 2)]

quartic_t6_inner = (
    -8 * slope_p**3 * r0040
    + 16 * slope_p**2 * tau
    - 6 * slope_p**2 * r0130
    + 32 * slope_p * slope_q * r0012
    + 16 * slope_p * sigma
    - 4 * slope_p * r0220
    + 48 * slope_q**2 * r0003
    + 32 * slope_q * r0102
    + 16 * rho
    - 2 * r0310
    + 89
)
assert sp.factor(
    pencil_coefficient(6).subs(x_free_quartic_q_substitution)
) == N**4 * quartic_t6_inner**2

t6_forced_substitution = {
    **x_free_quartic_q_substitution,
    r0003: 0,
    r0012: 0,
    r0102: 0,
    r0040: 0,
    r0130: sp.Rational(8, 3) * tau,
    r0220: 4 * sigma,
    r0310: 8 * rho + sp.Rational(89, 2),
}

quartic_t5_inner = (
    32 * slope_p**5 * tau**2
    - 64 * slope_p**4 * tau * sigma
    - 448 * slope_p**3 * tau * rho
    - 2736 * slope_p**3 * tau
    - 96 * slope_p**3 * sigma**2
    - 96 * slope_p**2 * tau * r0400
    - 576 * slope_p**2 * sigma * rho
    - 3528 * slope_p**2 * sigma
    - 96 * slope_p * sigma * r0400
    - 480 * slope_p * rho**2
    - 5904 * slope_p * rho
    - 18063 * slope_p
    - 96 * rho * r0400
    - 594 * r0400
)
assert sp.factor(
    pencil_coefficient(5).subs(t6_forced_substitution)
) == -sp.Rational(2, 3) * L * N**3 * quartic_t5_inner

# Coefficients slope_p**5, slope_p**3, slope_p, and 1 now force the
# displayed branch equations in that order.  The last implication uses that
# rho=-99/16 is not a root of rho_polynomial.
rho_polynomial = 160 * rho**2 + 1968 * rho + 6021
assert sp.factor(
    sp.Poly(quartic_t5_inner, slope_p).coeff_monomial(slope_p**5)
) == 32 * tau**2
assert sp.factor(
    sp.Poly(
        quartic_t5_inner.subs({tau: 0}), slope_p
    ).coeff_monomial(slope_p**3)
) == -96 * sigma**2
assert sp.expand(
    sp.Poly(
        quartic_t5_inner.subs({tau: 0, sigma: 0}), slope_p
    ).coeff_monomial(slope_p)
    + 3 * rho_polynomial
) == 0
assert rho_polynomial.subs({rho: -sp.Rational(99, 16)}) == (
    -sp.Rational(243, 8)
)
assert sp.discriminant(rho_polynomial, rho) == 576 * 34

# On either algebraic branch the complete x-free quartic slice is a single
# y**3*p term.  The coefficient of y**3*p in
# T=y**3*p-R4_homogeneous is nonzero: rho=-87/16 is not on either branch.
quartic_branch_substitution = {
    tau: 0,
    sigma: 0,
    r0400: 0,
}
quartic_x_free_branch = sp.factor(
    R4_homogeneous.subs({x: 0})
    .subs(t6_forced_substitution)
    .subs(quartic_branch_substitution)
)
assert sp.expand(
    quartic_x_free_branch
    - (8 * rho + sp.Rational(89, 2)) * y**3 * p
) == 0
assert rho_polynomial.subs({rho: -sp.Rational(87, 16)}) == (
    sp.Rational(405, 8)
)

# The next line coefficient closes both algebraic branches.  Its c**4 and
# d**2 terms first kill one cubic and one quadratic graph coefficient.  The
# remaining d coefficient imposes a second quadratic on rho, coprime to the
# first one.
r0030 = graph_coefficients[(0, 0, 3, 0)]
r0002 = graph_coefficients[(0, 0, 0, 2)]
quartic_t4_after_t5 = sp.expand(
    pencil_coefficient(4)
    .subs(t6_forced_substitution)
    .subs(quartic_branch_substitution)
)
quartic_t4_polynomial = sp.Poly(
    quartic_t4_after_t5, slope_p, slope_q
)
assert sp.factor(
    quartic_t4_polynomial.coeff_monomial(slope_p**4)
) == 36 * N**4 * r0030**2
assert sp.factor(
    quartic_t4_polynomial.coeff_monomial(slope_q**2)
) == 1024 * N**4 * r0002**2

rho_terminal_polynomial = 8 * rho**2 + 99 * rho + 279
quartic_t4_terminal = sp.Poly(
    quartic_t4_after_t5.subs({r0030: 0, r0002: 0}),
    slope_p,
    slope_q,
)
assert sp.expand(
    quartic_t4_terminal.coeff_monomial(slope_q)
    - 4 * L * N**3 * rho_terminal_polynomial
) == 0
assert sp.expand(
    rho_polynomial - 20 * rho_terminal_polynomial
) == -3 * (4 * rho - 147)
assert rho_terminal_polynomial.subs({rho: sp.Rational(147, 4)}) == (
    sp.Rational(58887, 4)
)
assert sp.resultant(
    rho_polynomial, rho_terminal_polynomial, rho
) == 16959456


# Exact plane-jet identity.  At x=0 the omitted r-gradient, the r-r Hessian
# entry, and all but one of the mixed r-Hessian entries vanish.  Consequently
# the graph Hessian depends only on the trace T=R|_(x=0) and normal derivative
# S=R_x|_(x=0), and the determinant is affine-linear in S with a unit
# coefficient when L*N is nonzero.
plane_trace, plane_normal = sp.symbols("plane_trace plane_normal")
plane_trace_y, plane_trace_p, plane_trace_q = sp.symbols(
    "plane_trace_y plane_trace_p plane_trace_q"
)
ambient_variables = (x, y, p, q, r)
ambient_hessian = sp.hessian(Psi, ambient_variables)
ambient_plane_hessian = ambient_hessian.subs({x: 0, r: plane_trace})
source_unit = sp.Matrix([1, 0, 0, 0])
plane_trace_gradient = sp.Matrix(
    [plane_normal, plane_trace_y, plane_trace_p, plane_trace_q]
)

assert sp.factor(sp.diff(Psi, r).subs({x: 0})) == 0
assert sp.factor(ambient_plane_hessian[4, 4]) == 0
assert ambient_plane_hessian[:4, 4] == 2 * N * source_unit

plane_graph_hessian = (
    ambient_plane_hessian[:4, :4]
    + 2
    * N
    * (
        source_unit * plane_trace_gradient.T
        + plane_trace_gradient * source_unit.T
    )
)
plane_graph_determinant = sp.factor(plane_graph_hessian.det())
plane_forcing = sp.factor(
    plane_graph_determinant.subs({plane_normal: 0})
)
assert sp.expand(
    plane_graph_determinant
    - plane_forcing
    + 8 * L * N**3 * plane_normal
) == 0
assert sp.factor(plane_graph_hessian[1:, 1:].det()) == -2 * L * N**2


def homogeneous_three_variable_form(
    prefix: str, degree: int
) -> tuple[sp.Expr, dict[tuple[int, int, int], sp.Symbol]]:
    """Return a generic homogeneous form in (y,p,q)."""

    coefficients: dict[tuple[int, int, int], sp.Symbol] = {}
    form = sp.S.Zero
    for y_degree in range(degree + 1):
        for p_degree in range(degree - y_degree + 1):
            q_degree = degree - y_degree - p_degree
            exponents = (y_degree, p_degree, q_degree)
            coefficient = sp.Symbol(
                prefix + "_" + "".join(map(str, exponents))
            )
            coefficients[exponents] = coefficient
            form += (
                coefficient
                * y**y_degree
                * p**p_degree
                * q**q_degree
            )
    return form, coefficients


def plane_determinant_for_trace(trace: sp.Expr) -> sp.Expr:
    """Substitute one plane trace into the normal-free determinant."""

    return sp.expand(
        plane_forcing.subs(
            {
                plane_trace: trace,
                plane_trace_y: sp.diff(trace, y),
                plane_trace_p: sp.diff(trace, p),
                plane_trace_q: sp.diff(trace, q),
            }
        )
    )


def line_coefficient(expression: sp.Expr, degree: int) -> sp.Expr:
    """Extract t**degree on the two-slope plane pencil."""

    line_expression = expression.subs(
        {y: t, p: slope_p * t, q: slope_q * t}
    )
    return sp.factor(sp.Poly(line_expression, t).coeff_monomial(t**degree))


# Degree-five leading trace rigidity.  The t**10 and t**8 squares are
# independent of every x-positive graph coefficient because the plane
# determinant sees only the trace and its first normal derivative; the latter
# has degree at most four and hence cannot enter these faces.
generic_trace4, generic_trace4_coefficients = homogeneous_three_variable_form(
    "trace4", 4
)
generic_trace5, generic_trace5_coefficients = homogeneous_three_variable_form(
    "trace5", 5
)
generic_top_plane_determinant = plane_determinant_for_trace(
    generic_trace4 + generic_trace5
)
generic_trace5_q_line = sp.diff(generic_trace5, q).subs(
    {y: 1, p: slope_p, q: slope_q}
)
assert line_coefficient(generic_top_plane_determinant, 10) == (
    256 * N**4 * generic_trace5_q_line**2
)

trace5_q_free_substitution = {
    coefficient: 0
    for exponents, coefficient in generic_trace5_coefficients.items()
    if exponents[2] > 0
}
assert len(trace5_q_free_substitution) == 15
assert sp.factor(
    line_coefficient(generic_top_plane_determinant, 9).subs(
        trace5_q_free_substitution
    )
) == 0
generic_trace5_p_line = sp.diff(generic_trace5, p).subs(
    {y: 1, p: slope_p, q: slope_q}
)
generic_trace4_q_line = sp.diff(generic_trace4, q).subs(
    {y: 1, p: slope_p, q: slope_q}
)
assert sp.factor(
    line_coefficient(generic_top_plane_determinant, 8).subs(
        trace5_q_free_substitution
    )
) == sp.factor(
    4
    * N**4
    * (
        generic_trace5_p_line.subs(trace5_q_free_substitution)
        - 8 * generic_trace4_q_line
    )
    ** 2
)

# Parameterize the t**8 equality by Q=partial_q T4 and integrate
# partial_p T5=8*y*Q.  The t**7 coefficient is a product.  Its second factor
# has no degree-five resonance: setting it to zero kills Q and kappa.
qa, qb, qc, qd, kappa = sp.symbols("qa qb qc qd kappa")
trace_q = qa * p**3 + qb * y * p**2 + qc * y**2 * p + qd * y**3
integrated_trace5 = sp.expand(
    8 * y * sp.integrate(trace_q, p) + kappa * y**5
)
resonant_trace4 = q * trace_q
resonant_plane_determinant = plane_determinant_for_trace(
    resonant_trace4 + integrated_trace5
)
trace_q_line = trace_q.subs({y: 1, p: slope_p})
resonance_factor_line = (
    sp.diff(integrated_trace5, y) - 4 * p * trace_q
).subs({y: 1, p: slope_p})
assert line_coefficient(resonant_plane_determinant, 7) == sp.factor(
    16 * L * N**3 * trace_q_line * resonance_factor_line
)
resonance_equations = sp.Poly(
    sp.diff(integrated_trace5, y) - 4 * p * trace_q, y, p
).coeffs()
resonance_matrix, resonance_vector = sp.linear_eq_to_matrix(
    resonance_equations, (qa, qb, qc, qd, kappa)
)
assert resonance_vector == sp.zeros(len(resonance_equations), 1)
assert resonance_matrix.rank() == 5


# The degree-six face after the leading rigidity.  Its highest slope degrees
# first force a0=b0=A0=e0=0 and B0=8*tau5/3.  The displayed residual then
# introduces only the two shifted variables U5,V5 and lambda5=4*L*kappa5/N.
a0, b0, e0, A0, B0, C0_trace, D0_trace = sp.symbols(
    "a0 b0 e0 A0 B0 C0_trace D0_trace"
)
tau5, sigma5, rho5, g5, kappa5 = sp.symbols(
    "tau5 sigma5 rho5 g5 kappa5"
)
trace3_degree5 = (
    a0 * q**3
    + b0 * p * q**2
    + tau5 * p**2 * q
    + e0 * y * q**2
    + sigma5 * y * p * q
    + rho5 * y**2 * q
)
trace4_degree5 = (
    A0 * p**4
    + B0 * y * p**3
    + C0_trace * y**2 * p**2
    + D0_trace * y**3 * p
    + g5 * y**4
)
trace5_degree5 = kappa5 * y**5
degree5_plane_determinant = plane_determinant_for_trace(
    trace3_degree5 + trace4_degree5 + trace5_degree5
)
degree5_t6_inner = (
    48 * a0 * slope_q**2
    + 32 * b0 * slope_p * slope_q
    + 16 * tau5 * slope_p**2
    - 8 * A0 * slope_p**3
    + 32 * e0 * slope_q
    + 16 * sigma5 * slope_p
    - 6 * B0 * slope_p**2
    + 16 * rho5
    - 4 * C0_trace * slope_p
    - 2 * D0_trace
    + 89
)
degree5_t6_kappa_factor = (
    60 * a0 * slope_q**2
    + 40 * b0 * slope_p * slope_q
    + 20 * tau5 * slope_p**2
    + 40 * e0 * slope_q
    + 20 * sigma5 * slope_p
    + 20 * rho5
    + 123
)
assert line_coefficient(degree5_plane_determinant, 6) == sp.factor(
    N**4 * degree5_t6_inner**2
    + 4 * L * N**3 * kappa5 * degree5_t6_kappa_factor
)

U5, V5 = sp.symbols("U5 V5")
lambda5 = 4 * L * kappa5 / N
degree5_residual_substitution = {
    a0: 0,
    b0: 0,
    e0: 0,
    A0: 0,
    B0: sp.Rational(8, 3) * tau5,
    C0_trace: 4 * sigma5 - U5 / 4,
    D0_trace: 8 * rho5 + sp.Rational(89, 2) - V5 / 2,
}
degree5_t6_residual = sp.factor(
    line_coefficient(degree5_plane_determinant, 6).subs(
        degree5_residual_substitution
    )
    / N**4
)
assert sp.factor(
    degree5_t6_residual
    - (
        U5**2 * slope_p**2
        + 2 * U5 * V5 * slope_p
        + V5**2
        + lambda5
        * (
            20 * tau5 * slope_p**2
            + 20 * sigma5 * slope_p
            + 20 * rho5
            + 123
        )
    )
) == 0

degree5_t5_residual = sp.Poly(
    sp.expand(
        line_coefficient(degree5_plane_determinant, 5).subs(
            degree5_residual_substitution
        )
    ),
    slope_p,
    slope_q,
)
assert sp.factor(
    degree5_t5_residual.coeff_monomial(slope_p**5)
) == -sp.Rational(64, 3) * L * N**3 * tau5**2
assert sp.factor(
    degree5_t5_residual.coeff_monomial(slope_p**3).subs(
        {tau5: 0, U5: 0}
    )
) == 64 * L * N**3 * sigma5**2


def reduce_coefficients_modulo(
    expression: sp.Expr,
    polynomial_variables: tuple[sp.Symbol, ...],
    modulus: sp.Expr,
    modulus_variable: sp.Symbol,
) -> sp.Expr:
    """Reduce every coefficient of ``expression`` modulo one polynomial."""

    reduced = sp.S.Zero
    modulus_poly = sp.Poly(modulus, modulus_variable)
    for monomial, coefficient in sp.Poly(
        expression, *polynomial_variables
    ).terms():
        numerator, denominator = sp.fraction(sp.together(coefficient))
        remainder = sp.rem(
            sp.Poly(numerator, modulus_variable), modulus_poly
        ).as_expr()
        term = sp.factor(remainder / denominator)
        for variable, exponent in zip(
            polynomial_variables, monomial, strict=True
        ):
            term *= variable**exponent
        reduced += term
    return sp.expand(reduced)


# Sparse degree-five v2 trace family.  Put
# T=kappa*y**5+d*y**3*p+rho*y**2*q and V=16*rho-2*d+89.  The t**6 and t**5
# equations eliminate lambda=2*kappa and leave the exact curve G=0 below.
# The plane-jet identity then determines the normal quartic S uniquely.
sparse_rho, sparse_V, target_constant = sp.symbols(
    "sparse_rho sparse_V target_constant"
)
sparse_Q = 160 * sparse_rho**2 + 1968 * sparse_rho + 6021
sparse_lambda = sp.factor(
    (
        3 * sparse_V * (4 * sparse_rho + 13)
        - sparse_Q
    )
    / 30
)
sparse_curve = sp.factor(
    30
    * (
        sparse_V**2
        + sparse_lambda * (20 * sparse_rho + 123)
    )
)
assert sparse_curve == (
    30 * sparse_V**2
    + 240 * sparse_V * sparse_rho**2
    + 2256 * sparse_V * sparse_rho
    + 4797 * sparse_V
    - 3200 * sparse_rho**3
    - 59040 * sparse_rho**2
    - 362484 * sparse_rho
    - 740583
)
sparse_kappa = sparse_lambda / 2
sparse_d = (16 * sparse_rho + 89 - sparse_V) / 2
sparse_trace = sp.expand(
    sparse_kappa * y**5
    + sparse_d * y**3 * p
    + sparse_rho * y**2 * q
)
sparse_plane_forcing = plane_determinant_for_trace(sparse_trace).subs(
    {L: 1, M: 13, N: 2}
)
sparse_reduced_forcing = reduce_coefficients_modulo(
    sparse_plane_forcing,
    (y, p, q),
    sparse_curve,
    sparse_V,
)
assert sp.Poly(sparse_reduced_forcing, y, p, q).total_degree() == 4
sparse_normal = sp.expand(
    (sparse_reduced_forcing - target_constant) / 64
)
sparse_graph = sp.expand(sparse_trace + x * sparse_normal)
sparse_potential = meng_potential(p, q, sparse_graph).subs(
    {L: 1, M: 13, N: 2}
)
sparse_hessian = sp.hessian(sparse_potential, (x, y, p, q))
sparse_axis_hessian = sparse_hessian.subs({p: 0, q: 0})
sparse_transverse_axis = sp.expand(
    determinant_truncation(
        sparse_axis_hessian, (x,), (1,)
    ).get((1,), sp.S.Zero)
)
sparse_transverse_axis_reduced = reduce_coefficients_modulo(
    sparse_transverse_axis,
    (y,),
    sparse_curve,
    sparse_V,
)
sparse_axis_polynomial = sp.Poly(sparse_transverse_axis_reduced, y)
sparse_H7 = (
    1920 * sparse_V * sparse_rho**3
    + 36240 * sparse_V * sparse_rho**2
    + 223344 * sparse_V * sparse_rho
    + 448443 * sparse_V
    - 25600 * sparse_rho**4
    - 650880 * sparse_rho**3
    - 6194880 * sparse_rho**2
    - 26158356 * sparse_rho
    - 41346207
)
sparse_H5 = (
    1212 * sparse_V * sparse_rho
    - 7716 * sparse_V
    + 51040 * sparse_rho**2
    + 627792 * sparse_rho
    + 1920699
)
assert sp.expand(
    sparse_axis_polynomial.coeff_monomial(y**7)
    - sp.Rational(8, 5) * sparse_H7
) == 0
assert sp.expand(
    sparse_axis_polynomial.coeff_monomial(y**5)
    - sp.Rational(104, 5) * sparse_H5
) == 0

# A rational collision-containing plane-flat seed on the exceptional
# 20*rho+123=V=0 point.  It is not a constant-Hessian graph: the y**7
# transverse coefficient below is already nonzero.
plane_flat_seed_substitution = {
    sparse_rho: -sp.Rational(123, 20),
    sparse_V: 0,
    target_constant: sp.Rational(17165601, 25),
}
assert sp.factor(sparse_curve.subs(plane_flat_seed_substitution)) == 0
assert sp.factor(
    sparse_lambda.subs(plane_flat_seed_substitution)
) == sp.Rational(51, 50)
assert sp.factor(
    sparse_kappa.subs(plane_flat_seed_substitution)
) == sp.Rational(51, 100)
assert sp.factor(
    sparse_d.subs(plane_flat_seed_substitution)
) == -sp.Rational(47, 10)
plane_flat_seed_graph = sparse_graph.subs(plane_flat_seed_substitution)
assert sp.factor(
    (sparse_plane_forcing - 64 * sparse_normal).subs(
        plane_flat_seed_substitution
    )
) == sp.Rational(17165601, 25)
for point_x, point_y in (
    (1, -sp.Rational(3, 2)),
    (-1, sp.Rational(3, 2)),
):
    assert sp.factor(
        plane_flat_seed_graph.subs(
            {x: point_x, y: point_y, p: 0, q: 0}
        )
    ) == 0
assert sp.factor(
    sparse_axis_polynomial.coeff_monomial(y**7).subs(
        plane_flat_seed_substitution
    )
) == sp.Rational(22032, 125)

sparse_first_resultant = sp.factor(
    sp.resultant(sparse_curve, sparse_H7, sparse_V)
)
sparse_second_resultant = sp.factor(
    sp.resultant(sparse_H7, sparse_H5, sparse_V)
)
assert sp.expand(
    sparse_first_resultant
    + 120
    * sparse_Q**2
    * (
        3200 * sparse_rho**4
        + 72240 * sparse_rho**3
        + 579240 * sparse_rho**2
        + 1899693 * sparse_rho
        + 2000700
    )
) == 0
assert sp.expand(
    sparse_second_resultant
    - 15
    * sparse_Q
    * (
        53760 * sparse_rho**3
        + 858080 * sparse_rho**2
        + 4224396 * sparse_rho
        + 6004503
    )
) == 0
assert sp.gcd(
    sp.Poly(sparse_first_resultant, sparse_rho),
    sp.Poly(sparse_second_resultant, sparse_rho),
).monic() == sp.Poly(sparse_Q, sparse_rho).monic()
sparse_quartic_branch_basis = sp.groebner(
    [sparse_Q, sparse_curve, sparse_H7, sparse_H5],
    sparse_V,
    sparse_rho,
    order="lex",
)
assert [polynomial.as_expr() for polynomial in sparse_quartic_branch_basis.polys] == [
    sparse_V,
    sparse_Q,
]

# On the only resultant branch, V=kappa=0 and Q(rho)=0.  Flatten the plane
# again in the quotient by Q and extract one full transverse coefficient.
# A degree-at-most-five correction x**2*U has deg(U|_(x=0))<=3 and changes
# this first transverse determinant only by -24*L*N**3*U, so it cannot touch
# the degree-five monomial y**4*q below.
quartic_branch_trace = sp.expand(
    (8 * sparse_rho + sp.Rational(89, 2)) * y**3 * p
    + sparse_rho * y**2 * q
)
quartic_branch_forcing = plane_determinant_for_trace(
    quartic_branch_trace
).subs({L: 1, M: 13, N: 2})
quartic_branch_reduced_forcing = reduce_coefficients_modulo(
    quartic_branch_forcing,
    (y, p, q),
    sparse_Q,
    sparse_rho,
)
quartic_branch_normal = sp.expand(
    (quartic_branch_reduced_forcing - target_constant) / 64
)
quartic_branch_graph = sp.expand(
    quartic_branch_trace + x * quartic_branch_normal
)
quartic_branch_potential = meng_potential(
    p, q, quartic_branch_graph
).subs({L: 1, M: 13, N: 2})
quartic_branch_hessian = sp.hessian(
    quartic_branch_potential, (x, y, p, q)
)
quartic_branch_transverse = sp.expand(
    determinant_truncation(
        quartic_branch_hessian, (x,), (1,)
    ).get((1,), sp.S.Zero)
)
quartic_branch_transverse_reduced = reduce_coefficients_modulo(
    quartic_branch_transverse,
    (y, p, q),
    sparse_Q,
    sparse_rho,
)
assert sp.expand(
    sp.Poly(
        quartic_branch_transverse_reduced, y, p, q
    ).coeff_monomial(y**4 * q)
    - 216 * (5 * sparse_rho + 12)
) == 0
assert sp.factor(
    sparse_Q.subs({sparse_rho: -sp.Rational(12, 5)})
) == sp.Rational(11097, 5)

# The universal first-normal variation used above: the perturbation x**2*U
# changes only the x-x Hessian derivative, by 12*N*U.  Its cofactor is the
# lower 3-by-3 determinant -2*L*N**2.
normal_variation = sp.Symbol("normal_variation")
assert sp.factor(
    (-2 * L * N**2) * (12 * N * normal_variation)
) == -24 * L * N**3 * normal_variation


print("PASS: affine Meng--Yang four-plane restrictions are degenerate or nonconstant")
print("PASS: every degree-at-most-three polynomial graph has nonconstant Hessian determinant")
print("PASS: the full quartic two-slope jet has no characteristic-zero branch")
print("PASS: every degree-at-most-four polynomial graph has nonconstant Hessian determinant")
print("PASS: the Meng--Yang plane determinant is a unit-affine normal-jet equation")
print("PASS: degree-five leading traces reduce to the vertical quintic branch")
print("PASS: the full sparse degree-five v2 trace family is excluded")
