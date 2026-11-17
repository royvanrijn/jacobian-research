#!/usr/bin/env python3
"""Exact audit of the canonical D5 two-mask blowdown.

This verifies the automorphic obstruction's algebraic identities and
exhausts the 72 coordinate assignments from the unchanged and two
ramification-incidence source charts.  The genus-two argument itself is
proved in the accompanying note.
"""

from __future__ import annotations

from itertools import permutations

import sympy as sp


U, V, A, B, S, T, Z, W = sp.symbols("U V A B S T Z W")
Delta = V**2 - 4 * U**5
mask_matrix = sp.Matrix([[V, 2 * U], [2 * U**4, V]])
assert sp.expand(mask_matrix.det() - Delta) == 0

mask_outputs = mask_matrix * sp.Matrix([A, B])
inverse_masks = sp.cancel(mask_matrix.adjugate() * sp.Matrix([S, T]) / Delta)
expected_inverse_masks = sp.Matrix(
    [(V * S - 2 * U * T) / Delta, (-2 * U**4 * S + V * T) / Delta]
)
assert all(
    sp.cancel(entry) == 0
    for entry in inverse_masks - expected_inverse_masks
)
print("PASS: the canonical two-mask blowdown has determinant Delta")

# Constant-linear base rigidity.  The displayed coefficient sequence proves
# p=cU and q=dV for every nonzero scalar multiple of Delta.
pU, pV, pZ, pW = sp.symbols("pU pV pZ pW")
qU, qV, qZ, qW, lam = sp.symbols("qU qV qZ qW lam")
p = pU * U + pV * V + pZ * Z + pW * W
q = qU * U + qV * V + qZ * Z + qW * W
error = sp.Poly(sp.expand(q**2 - 4 * p**5 - lam * Delta), U, V, Z, W)

assert error.coeff_monomial(Z**5) == -4 * pZ**5
assert error.coeff_monomial(W**5) == -4 * pW**5
after_p_masks = sp.Poly(
    error.as_expr().subs({pZ: 0, pW: 0}), U, V, Z, W
)
assert after_p_masks.coeff_monomial(Z**2) == qZ**2
assert after_p_masks.coeff_monomial(W**2) == qW**2
after_all_masks = sp.Poly(
    after_p_masks.as_expr().subs({qZ: 0, qW: 0}), U, V, Z, W
)
assert after_all_masks.coeff_monomial(V**5) == -4 * pV**5
after_pV = sp.Poly(after_all_masks.as_expr().subs(pV, 0), U, V)
assert after_pV.coeff_monomial(U**2) == qU**2
after_qU = sp.Poly(after_pV.as_expr().subs(qU, 0), U, V)
assert after_qU.coeff_monomial(V**2) == qV**2 - lam
assert after_qU.coeff_monomial(U**5) == 4 * (lam - pU**5)
print("PASS: constant-Jacobian linear base rows are only cusp scalings")

# With normalized base rows p=U,q=V, linear inverse divisibility makes the
# mask rows dependent on the base rows.
s_coefficients = sp.symbols("s0:4")
t_coefficients = sp.symbols("t0:4")
linear_variables = (U, V, Z, W)
linear_S = sum(
    coefficient * variable
    for coefficient, variable in zip(s_coefficients, linear_variables)
)
linear_T = sum(
    coefficient * variable
    for coefficient, variable in zip(t_coefficients, linear_variables)
)
first_numerator = sp.Poly(
    sp.expand(V * linear_S - 2 * U * linear_T),
    *linear_variables,
)
linear_solution = sp.linsolve(
    first_numerator.coeffs(),
    (*s_coefficients, *t_coefficients),
)
assert linear_solution == {
    (
        2 * t_coefficients[1],
        0,
        0,
        0,
        0,
        t_coefficients[1],
        0,
        0,
    )
}
print("PASS: linear adjugate divisibility forces a singular rechart")

# The generic fibre used in the all-degree rigidity proof is the odd-degree
# hyperelliptic curve y^2=4*x^5+delta.  Its quintic is squarefree for
# nonzero transcendental delta, and its smooth projective genus is two.
curve_x, curve_delta = sp.symbols("curve_x curve_delta")
curve_polynomial = 4 * curve_x**5 + curve_delta
assert sp.resultant(
    curve_polynomial,
    sp.diff(curve_polynomial, curve_x),
    curve_x,
) == 4**5 * 5**5 * curve_delta**4
assert (sp.degree(curve_polynomial, curve_x) - 1) // 2 == 2
print("PASS: the generic fixed-discriminant fibre is a smooth genus-two curve")

# The first nonautomorphic normalized-cusp chart.  Its desired log-crepant
# right side is polynomial but not divisible by Delta, whereas every
# Hamiltonian contraction determinant has a factor Delta.
old_u, old_v, mask_s, mask_t = sp.symbols(
    "old_u old_v mask_s mask_t"
)
old_delta = old_v**2 - 4 * old_u**5
cusp_p = old_v**2 + old_delta * mask_s
cusp_q = 2 * old_v**5 + old_delta * mask_t
cusp_quotient = sp.cancel(
    (cusp_q**2 - 4 * cusp_p**5) / old_delta
)
assert sp.denom(cusp_quotient) == 1

cusp_variables = (old_v, mask_s, mask_t)
contraction_vector = sp.Matrix(
    [sp.diff(cusp_p, variable) for variable in cusp_variables]
).cross(
    sp.Matrix(
        [sp.diff(cusp_q, variable) for variable in cusp_variables]
    )
)
assert all(
    sp.rem(entry, old_delta, old_v) == 0
    for entry in contraction_vector
)
cusp_basis = sp.groebner(
    [old_delta],
    old_u,
    old_v,
    mask_s,
    mask_t,
)
cusp_remainder = sp.factor(cusp_basis.reduce(cusp_quotient)[1])
assert cusp_remainder == -4 * old_v**5 * (
    5 * mask_s * old_v**3 - mask_t
)
print("PASS: the first nonautomorphic cusp contraction has a divisor mismatch")

# The minimal tangential chart h=V+S,
# p=h^2+a*Delta*T, q=2*h^5+b*Delta*T.  The T^5 and T^2 equations force
# a=b=0, leaving zero Jacobian.
normal_a, normal_b, log_scale = sp.symbols(
    "normal_a normal_b log_scale"
)
tangent_h = old_v + mask_s
tangent_p = tangent_h**2 + normal_a * old_delta * mask_t
tangent_q = 2 * tangent_h**5 + normal_b * old_delta * mask_t
tangent_map = sp.Matrix([tangent_p, tangent_q, old_u, mask_s])
tangent_jacobian = sp.factor(
    tangent_map.jacobian((old_u, old_v, mask_s, mask_t)).det()
)
assert sp.expand(
    tangent_jacobian
    - 2
    * tangent_h
    * old_delta
    * (5 * normal_a * tangent_h**3 - normal_b)
) == 0
tangent_error = sp.Poly(
    sp.expand(
        tangent_q**2
        - 4 * tangent_p**5
        - log_scale * old_delta * tangent_jacobian
    ),
    mask_t,
)
assert sp.expand(
    tangent_error.coeff_monomial(mask_t**5)
    + 4 * normal_a**5 * old_delta**5
) == 0
after_normal_a = sp.Poly(
    tangent_error.as_expr().subs(normal_a, 0),
    mask_t,
)
assert sp.expand(
    after_normal_a.coeff_monomial(mask_t**2)
    - normal_b**2 * old_delta**2
) == 0
assert tangent_jacobian.subs({normal_a: 0, normal_b: 0}) == 0
print("PASS: the minimal tangential log-crepant chart is degenerate")


sqrt5 = sp.sqrt(5)
alpha = (3 + sqrt5) / 2
beta = (3 - sqrt5) / 2


def dickson_data(a: sp.Symbol, u: sp.Expr) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    polynomial = sp.expand(a**5 - 5 * a**3 * u + 5 * a * u**2)
    derivative_factor = sp.expand(a**4 - 3 * a**2 * u + u**2)
    unramified_factor = sp.expand(a**2 - 4 * u)
    return polynomial, derivative_factor, unramified_factor


def divides(
    denominator: sp.Expr,
    numerator: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> bool:
    if numerator == 0:
        return True
    denominator_poly = sp.Poly(denominator, *variables, extension=sqrt5)
    numerator_poly = sp.Poly(numerator, *variables, extension=sqrt5)
    if numerator_poly.total_degree() < denominator_poly.total_degree():
        return False
    basis = sp.groebner([denominator], *variables, extension=sqrt5)
    return basis.reduce(numerator)[1] == 0


def screen_chart(
    outputs: tuple[sp.Expr, ...],
    names: tuple[str, ...],
    variables: tuple[sp.Symbol, ...],
    expected_branch: sp.Expr,
) -> tuple[int, int]:
    polynomial_assignments = 0
    determinant_assignments = 0
    for assignment in permutations(range(4)):
        base_u, base_v, mask_s, mask_t = (
            outputs[index] for index in assignment
        )
        denominator = sp.expand(base_v**2 - 4 * base_u**5)
        if sp.expand(denominator - expected_branch) == 0:
            determinant_assignments += 1
        first = sp.expand(base_v * mask_s - 2 * base_u * mask_t)
        second = sp.expand(
            -2 * base_u**4 * mask_s + base_v * mask_t
        )
        if divides(denominator, first, variables) and divides(
            denominator, second, variables
        ):
            polynomial_assignments += 1
    assert polynomial_assignments == 0, names
    assert determinant_assignments == 2, names
    return polynomial_assignments, determinant_assignments


# The unchanged primitive source chart.
a0, u0, x0, y0 = sp.symbols("a0 u0 x0 y0")
P0, Q0, C0 = dickson_data(a0, u0)
branch0 = sp.expand(C0 * Q0**2)
chart0 = (u0, P0, sp.expand(C0 * Q0 * x0), y0)

total_polynomial = 0
total_determinant = 0
counts = screen_chart(
    chart0,
    ("u", "P", "CQx", "y"),
    (a0, u0, x0, y0),
    branch0,
)
total_polynomial += counts[0]
total_determinant += counts[1]

# The two charts whose mask-zero section maps to R_plus or R_minus.
s0, t0, x1, y1 = sp.symbols("s0 t0 x1 y1")
for label, gamma in (("plus", alpha), ("minus", beta)):
    incidence_u = sp.expand(s0**2 / gamma + x1)
    incidence_P, incidence_Q, incidence_C = dickson_data(s0, incidence_u)
    incidence_branch = sp.expand(incidence_C * incidence_Q**2)
    incidence_outputs = (
        incidence_u,
        incidence_P,
        t0,
        sp.expand(incidence_C * incidence_Q * y1),
    )
    counts = screen_chart(
        incidence_outputs,
        (f"u_{label}", f"P_{label}", "t", "CQy"),
        (s0, t0, x1, y1),
        incidence_branch,
    )
    total_polynomial += counts[0]
    total_determinant += counts[1]

assert total_polynomial == 0
assert total_determinant == 6
print("PASS: none of the 72 primitive coordinate assignments is polynomial")
print("PASS: only six assignments pass the determinant-divisor prefilter")
print("PASS canonical D5 two-mask obstructions")
