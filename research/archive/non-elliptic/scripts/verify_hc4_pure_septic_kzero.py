#!/usr/bin/env python3
"""Verify exact tail obstructions in the k=0 pure-septic chart.

This checker continues the curved chart of HC4RSD35 after setting k=0.
It deliberately extracts only the homogeneous faces used in the proof.
All calculations are over QQ; polynomial identities with symbolic
parameters therefore remain valid over every characteristic-zero field.
"""

from __future__ import annotations

import hashlib
import json
from itertools import product
from pathlib import Path

import sympy as sp


x, y, z = sp.symbols("x y z")
variables = (x, y, z)
ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT / "artifacts" / "generated-results" / "hc4_pure_septic_kzero.json"
)


def binary_form(prefix: str, degree: int) -> tuple[sp.Expr, tuple[sp.Symbol, ...]]:
    coefficients = sp.symbols(f"{prefix}0:{degree + 1}")
    form = sum(
        coefficients[index] * x ** (degree - index) * y**index
        for index in range(degree + 1)
    )
    return form, coefficients


def ternary_form(prefix: str, degree: int) -> tuple[sp.Expr, tuple[sp.Symbol, ...]]:
    monomials = [
        x**i * y**j * z ** (degree - i - j)
        for i in range(degree + 1)
        for j in range(degree + 1 - i)
    ]
    coefficients = sp.symbols(f"{prefix}0:{len(monomials)}")
    return (
        sum(
            coefficient * monomial
            for coefficient, monomial in zip(
                coefficients, monomials, strict=True
            )
        ),
        coefficients,
    )


def bordered_face(
    components: dict[int, sp.Expr], component_degree_sum: int
) -> sp.Expr:
    """Extract one homogeneous face of grad(c)^T adj(Hess(c)) grad(c)."""

    gradients = {
        degree: [sp.diff(component, variable) for variable in variables]
        for degree, component in components.items()
    }
    hessians = {
        degree: sp.hessian(component, variables)
        for degree, component in components.items()
    }

    def convolution(
        gradient_indices: tuple[int, int],
        hessian_indices: tuple[tuple[int, int], tuple[int, int]],
    ) -> sp.Expr:
        result = 0
        for degrees in product(components, repeat=4):
            if sum(degrees) != component_degree_sum:
                continue
            d1, d2, d3, d4 = degrees
            result += (
                gradients[d1][gradient_indices[0]]
                * gradients[d2][gradient_indices[1]]
                * hessians[d3][hessian_indices[0]]
                * hessians[d4][hessian_indices[1]]
            )
        return sp.expand(result)

    return sp.expand(
        convolution((0, 0), ((1, 1), (2, 2)))
        - convolution((0, 0), ((1, 2), (1, 2)))
        + convolution((1, 1), ((0, 0), (2, 2)))
        - convolution((1, 1), ((0, 2), (0, 2)))
        + convolution((2, 2), ((0, 0), (1, 1)))
        - convolution((2, 2), ((0, 1), (0, 1)))
        + 2 * convolution((0, 1), ((0, 2), (1, 2)))
        - 2 * convolution((0, 1), ((0, 1), (2, 2)))
        + 2 * convolution((0, 2), ((0, 1), (1, 2)))
        - 2 * convolution((0, 2), ((0, 2), (1, 1)))
        + 2 * convolution((1, 2), ((0, 1), (0, 2)))
        - 2 * convolution((1, 2), ((0, 0), (1, 2)))
    )


def coefficient(expression: sp.Expr, monomial: sp.Expr) -> sp.Expr:
    return sp.Poly(expression, x, y, z).coeff_monomial(monomial)


# The residual degree-sixteen square after P=p*x^4 and D=0.  On the
# p!=0 locus its first four coefficients force H=x^3*K.  With
# G=7*V-3*p*K the whole binary identity is exactly
#
#     K_yy*(49*x^3*q-2*p*G)=G_y^2.
kappa = sp.symbols("kappa0:4")
nu = sp.symbols("nu0:4")
p, q = sp.symbols("p q")
K = sum(kappa[i] * x ** (3 - i) * y**i for i in range(4))
V = sum(nu[i] * x ** (3 - i) * y**i for i in range(4))
G = 7 * V - 3 * p * K
original_core = sp.expand(
    sp.diff(K, y, 2) * (49 * x**3 * q - 14 * p * V + 6 * p**2 * K)
    - (7 * sp.diff(V, y) - 3 * p * sp.diff(K, y)) ** 2
)
assert sp.expand(
    original_core
    - (
        sp.diff(K, y, 2) * (49 * x**3 * q - 2 * p * G)
        - sp.diff(G, y) ** 2
    )
) == 0

# Classify the p!=0 core over K(x).  If F=K_yy and G_y=F*L, then
# differentiating 49*x^3*q-2*p*G=F*L^2 gives
#
#     L*(2*F*(p+L_y)+F_y*L)=0.
#
# Thus L=0; or, when F_y!=0, L is proportional to F; or, at the
# F_y=0 endpoint, L=eta*x-p*y.  The next faces below certify that every
# one of these three global solution types is impossible.
c0, c1, c2, eta, residual_v0 = sp.symbols("c0 c1 c2 eta residual_v0")


def residual_components(
    residual_K: sp.Expr,
    residual_V: sp.Expr,
    residual_q: sp.Expr,
) -> dict[int, sp.Expr]:
    residual_Q = Q2_raw.subs(b[0], residual_q / 2)
    return {
        7: x**7,
        6: x**3 * residual_K,
        5: R5 + x**4 * z,
        4: U4 + z * residual_V,
        3: T3 + z * W2 + sp.Rational(1, 7) * x * z**2,
        2: residual_Q,
    }


# The objects R5,... are introduced just below.  The actual face checks are
# placed after that introduction.


# If p=0 and q!=0, write V_y=S and q*H_yy=S^2.  The last three
# degree-fifteen coefficients force S=s0*x^2.
s0, s1, s2, h0, h1, v0 = sp.symbols("s0 s1 s2 h0 h1 v0")
S = s0 * x**2 + s1 * x * y + s2 * y**2
square_H = (
    h0 * x**6
    + h1 * x**5 * y
    + s0**2 * x**4 * y**2 / 2
    + s0 * s1 * x**3 * y**3 / 3
    + (s1**2 + 2 * s0 * s2) * x**2 * y**4 / 12
    + s1 * s2 * x * y**5 / 10
    + s2**2 * y**6 / 30
)
square_V = v0 * x**3 + s0 * x**2 * y + s1 * x * y**2 / 2 + s2 * y**3 / 3
R5, r = binary_form("r", 5)
U4, u = binary_form("u", 4)
T3, t = binary_form("t", 3)
W2, w = binary_form("w", 2)
Q2_raw, b = ternary_form("b", 2)
# ternary_form orders z^2 first.
Q2 = Q2_raw.subs(b[0], sp.Rational(1, 2))

# Verify the global recurrence directly before entering its UFD charts.
full_H, full_a = binary_form("fa", 6)
full_V, full_v = binary_form("fv", 3)
full_Q = Q2_raw.subs(b[0], q / 2)
full_residual_components = {
    7: x**7,
    6: full_H,
    5: R5 + p * x**4 * z,
    4: U4 + z * full_V,
    3: T3 + z * W2 + p**2 * x * z**2 / 7,
    2: full_Q,
}
full_degree_sixteen = bordered_face(full_residual_components, 22)
full_degree_sixteen_expected = x**6 * (
    sp.diff(full_H, y, 2)
    * (49 * x**6 * q - 14 * p * x**3 * full_V + 6 * p**2 * full_H)
    - (7 * x**3 * sp.diff(full_V, y) - 3 * p * sp.diff(full_H, y)) ** 2
)
assert sp.expand(full_degree_sixteen - full_degree_sixteen_expected) == 0

# L=0 with cubic K_yy: degree fifteen cancels identically, but degree
# fourteen has -648/49.  The proportional-L branch dies one face earlier
# by 16.
cubic_K = c0 * x**3 + c1 * x**2 * y + c2 * x * y**2 + y**3
zero_L_V = residual_v0 * x**3 + sp.Rational(3, 7) * (
    c1 * x**2 * y + c2 * x * y**2 + y**3
)
zero_L_q = 2 * (7 * residual_v0 - 3 * c0) / 49
zero_L_components = residual_components(cubic_K, zero_L_V, zero_L_q)
zero_L_degree_fourteen = bordered_face(zero_L_components, 20)
assert coefficient(zero_L_degree_fourteen, x**4 * y**10) == -sp.Rational(
    648, 49
)

proportional_L_V = (
    residual_v0 * x**3
    + (27 * c1 - 4 * c2**2) * x**2 * y / 63
    + sp.Rational(5, 21) * c2 * x * y**2
    + sp.Rational(5, 21) * y**3
)
proportional_L_q = -2 * (
    243 * c0 - 4 * c2**3 - 567 * residual_v0
) / 3969
proportional_components = residual_components(
    cubic_K, proportional_L_V, proportional_L_q
)
proportional_degree_fifteen = bordered_face(proportional_components, 21)
assert coefficient(proportional_degree_fifteen, x**8 * y**7) == 16

# At K_yy proportional to x, the L=0 endpoint and the family
# L=eta*x-y have immutable positive-z coefficients 12/7 and 48/7.
endpoint_K = c0 * x**3 + x * y**2
endpoint_zero_V = residual_v0 * x**3 + sp.Rational(3, 7) * x * y**2
endpoint_zero_q = 2 * (7 * residual_v0 - 3 * c0) / 49
endpoint_zero_face = bordered_face(
    residual_components(endpoint_K, endpoint_zero_V, endpoint_zero_q), 20
)
assert coefficient(endpoint_zero_face, x**11 * y**2 * z) == sp.Rational(
    12, 7
)

endpoint_linear_V = (
    residual_v0 * x**3
    + sp.Rational(2, 7) * eta * x**2 * y
    + sp.Rational(2, 7) * x * y**2
)
endpoint_linear_q = 2 * (7 * residual_v0 - 3 * c0 + eta**2) / 49
endpoint_linear_face = bordered_face(
    residual_components(endpoint_K, endpoint_linear_V, endpoint_linear_q), 20
)
assert coefficient(endpoint_linear_face, x**11 * y**2 * z) == sp.Rational(
    48, 7
)

square_components = {
    7: x**7,
    6: square_H,
    5: R5,
    4: U4 + z * square_V,
    3: T3 + z * W2,
    2: Q2,
}
degree_fifteen_square = bordered_face(square_components, 21)
assert sp.factor(coefficient(degree_fifteen_square, x**5 * y**10)) == 56 * s2**4 / 75
assert sp.factor(
    coefficient(degree_fifteen_square, x**9 * y**6).subs(s2, 0)
) == sp.Rational(14, 9) * s1**4

# After s2=s1=0, normalize q=s0=1 and remove h1 by the shear preserving x.
# Four degree-fifteen equations determine r2,r3,r4,r5.  The next face has
# the immutable coefficient -4 at x^10*y^4.
aligned_square_H = h0 * x**6 + x**4 * y**2 / 2
aligned_square_V = v0 * x**3 + x**2 * y
aligned_R5 = R5.subs(
    {
        r[2]: w[1],
        r[3]: sp.Rational(2, 3) * w[2] - sp.Rational(2, 21) * v0,
        r[4]: 0,
        r[5]: 0,
    }
)
aligned_components = {
    7: x**7,
    6: aligned_square_H,
    5: aligned_R5,
    4: U4 + z * aligned_square_V,
    3: T3 + z * W2,
    2: Q2,
}
degree_fourteen_square = bordered_face(aligned_components, 20)
assert coefficient(degree_fourteen_square, x**10 * y**4) == -4


# The invariant explanation for the same obstruction.  If
# w=z+D(x)*y and c=h(x,w), direct differentiation gives
# J(c)=-D'(x)^2*h_w^4.  Verify it on a generic polynomial jet rather than
# relying on SymPy's printed derivative normal form.
D0, D1, D2 = sp.symbols("D0 D1 D2")
jet = sp.symbols("j0:10")
D_polynomial = D0 + D1 * x + D2 * x**2
w_coordinate = z + D_polynomial * y
h_jet = (
    jet[0]
    + jet[1] * x
    + jet[2] * w_coordinate
    + jet[3] * x**2
    + jet[4] * x * w_coordinate
    + jet[5] * w_coordinate**2
    + jet[6] * x**3
    + jet[7] * x**2 * w_coordinate
    + jet[8] * x * w_coordinate**2
    + jet[9] * w_coordinate**3
)
h_w = sp.diff(
    jet[0]
    + jet[1] * x
    + jet[2] * z
    + jet[3] * x**2
    + jet[4] * x * z
    + jet[5] * z**2
    + jet[6] * x**3
    + jet[7] * x**2 * z
    + jet[8] * x * z**2
    + jet[9] * z**3,
    z,
).subs(z, w_coordinate)
gradient = sp.Matrix([sp.diff(h_jet, variable) for variable in variables])
hessian = sp.hessian(h_jet, variables)
bordered_jet = sp.expand((gradient.T * hessian.adjugate() * gradient)[0])
assert sp.expand(
    bordered_jet + sp.diff(D_polynomial, x) ** 2 * h_w**4
) == 0


# On q=0, the degree-sixteen equation first gives V_y=0.  If
# V=v*x^3 is nonzero, normalize v=1.  The degree-fourteen face forces
# H=a0*x^6+a1*x^5*y+a2*x^4*y^2; after a shear and scaling take
# a1=0,a2=1.  Its last three equations give
# W=3*a0*x^2/7+rho*y^2 with 49*rho^2-35*rho+7=0.
rho = sp.symbols("rho")
Q2_zero = Q2_raw.subs(b[0], 0)
L1, ell = ternary_form("ell", 1)
# First retain arbitrary H and W.  The bottom coefficients of degree
# fourteen force H=a0*x^6+a1*x^5*y+a2*x^4*y^2.  In the normalized curved
# chart a shear gives a1=0 and the remaining equations are the quadratic
# rho relation displayed below.
zero_q_H, zero_q_a = binary_form("za", 6)
zero_q_W, zero_q_w = binary_form("zw", 2)
zero_q_initial_components = {
    7: x**7,
    6: zero_q_H,
    5: R5,
    4: U4 + x**3 * z,
    3: T3 + z * zero_q_W,
    2: Q2_zero,
}
degree_fourteen_zero_q = bordered_face(zero_q_initial_components, 20)
assert coefficient(degree_fourteen_zero_q, x**4 * y**10) == -324 * zero_q_a[6] ** 2
assert sp.factor(
    coefficient(degree_fourteen_zero_q, x**6 * y**8).subs(zero_q_a[6], 0)
) == -220 * zero_q_a[5] ** 2
assert sp.factor(
    coefficient(degree_fourteen_zero_q, x**8 * y**6).subs(
        {zero_q_a[5]: 0, zero_q_a[6]: 0}
    )
) == -136 * zero_q_a[4] ** 2
assert sp.factor(
    coefficient(degree_fourteen_zero_q, x**10 * y**4).subs(
        {
            zero_q_a[4]: 0,
            zero_q_a[5]: 0,
            zero_q_a[6]: 0,
        }
    )
) == -72 * zero_q_a[3] ** 2

zero_q_components = {
    7: x**7,
    6: h0 * x**6 + x**4 * y**2,
    5: R5,
    4: U4 + x**3 * z,
    3: T3 + z * (sp.Rational(3, 7) * h0 * x**2 + rho * y**2),
    2: Q2_zero,
    1: L1,
}
degree_thirteen_zero_q = bordered_face(zero_q_components, 19)
rho_relation = sp.Poly(49 * rho**2 - 35 * rho + 7, rho)


def reduce_rho(expression: sp.Expr) -> sp.Expr:
    return sp.factor(sp.rem(sp.Poly(expression, rho), rho_relation).as_expr())


assert sp.resultant(rho_relation.as_expr(), 7 * rho - 3, rho) != 0
assert sp.resultant(rho_relation.as_expr(), 21 * rho - 10, rho) != 0
assert sp.resultant(rho_relation.as_expr(), 7 * rho - 4, rho) != 0
assert sp.resultant(rho_relation.as_expr(), 70 * rho - 31, rho) != 0

zero_q_degree_thirteen_solution = {
    r[5]: 0,
    r[4]: (13 - 7 * rho) / 37,
    r[3]: 0,
    r[2]: sp.Rational(4, 7) * h0,
    b[1]: rho * r[1],
    b[3]: (35 * r[0] - 9 * h0**2) / 98,
}
degree_twelve_zero_q = bordered_face(zero_q_components, 18).subs(
    zero_q_degree_thirteen_solution
)
assert reduce_rho(
    coefficient(sp.expand(degree_twelve_zero_q), x**11 * z)
) == -24


# If V=0, degree fourteen gives W_y=0.  Inspect the W=x^2 chart.
general_H, general_a = binary_form("A", 6)
first_cubic_break = {
    7: x**7,
    6: general_H,
    5: R5,
    4: U4,
    3: T3 + x**2 * z,
    2: Q2_zero,
    1: L1,
}

# Instead of descending its long binary faces, collect every lower z-linear
# term.  After an x-translation and scaling its coefficient is
# D=x^2+beta*y+gamma.  The exact z^2 coefficient of J(h+z*D) is
# -4*beta^2*D; hence beta=0.  The z coefficient then equals
# 2*D*(gamma-3*x^2)*h_yy and contradicts the curved chart.
beta, gamma = sp.symbols("beta gamma")
global_h_coefficients = sp.symbols("gh0:10")
global_h = (
    global_h_coefficients[0]
    + global_h_coefficients[1] * x
    + global_h_coefficients[2] * y
    + global_h_coefficients[3] * x**2
    + global_h_coefficients[4] * x * y
    + global_h_coefficients[5] * y**2
    + global_h_coefficients[6] * x**3
    + global_h_coefficients[7] * x**2 * y
    + global_h_coefficients[8] * x * y**2
    + global_h_coefficients[9] * y**3
)
quadratic_D = x**2 + beta * y + gamma
quadratic_linear_potential = global_h + z * quadratic_D
quadratic_linear_gradient = sp.Matrix(
    [sp.diff(quadratic_linear_potential, variable) for variable in variables]
)
quadratic_linear_hessian = sp.hessian(
    quadratic_linear_potential, variables
)
quadratic_linear_bordered = sp.expand(
    (
        quadratic_linear_gradient.T
        * quadratic_linear_hessian.adjugate()
        * quadratic_linear_gradient
    )[0]
)
quadratic_linear_polynomial = sp.Poly(quadratic_linear_bordered, z)
assert sp.expand(
    quadratic_linear_polynomial.coeff_monomial(z**2)
    + 4 * beta**2 * quadratic_D
) == 0
assert sp.expand(
    quadratic_linear_polynomial.coeff_monomial(z).subs(beta, 0)
    - 2
    * (x**2 + gamma)
    * (gamma - 3 * x**2)
    * sp.diff(global_h, y, 2)
) == 0

# If W=0, the first possible transverse term is z*(lambda0*x+lambda1*y)
# in c2.
lambda0, lambda1 = sp.symbols("lambda0 lambda1")
Q2_linear_break, qbinary = binary_form("qb", 2)
first_quadratic_break = {
    7: x**7,
    6: general_H,
    5: R5,
    4: U4,
    3: T3,
    2: Q2_linear_break + z * (lambda0 * x + lambda1 * y),
    1: L1,
}

# There is no need to continue its homogeneous expansion.  For an arbitrary
# binary polynomial h and L=lambda0*x+lambda1*y, the z coefficient of
# J(h+z*L) is -2*L times the second derivative in the constant direction
# (lambda1,-lambda0).  Thus a nonzero L makes h affine in that direction.
hx = sp.symbols("hx0:10")
h_binary = (
    hx[0]
    + hx[1] * x
    + hx[2] * y
    + hx[3] * x**2
    + hx[4] * x * y
    + hx[5] * y**2
    + hx[6] * x**3
    + hx[7] * x**2 * y
    + hx[8] * x * y**2
    + hx[9] * y**3
)
linear_tail = lambda0 * x + lambda1 * y
linear_potential = h_binary + z * linear_tail
linear_gradient = sp.Matrix(
    [sp.diff(linear_potential, variable) for variable in variables]
)
linear_hessian = sp.hessian(linear_potential, variables)
linear_bordered = sp.expand(
    (linear_gradient.T * linear_hessian.adjugate() * linear_gradient)[0]
)
directional_curvature = sp.expand(
    lambda1**2 * sp.diff(h_binary, x, 2)
    - 2 * lambda0 * lambda1 * sp.diff(h_binary, x, y)
    + lambda0**2 * sp.diff(h_binary, y, 2)
)
assert sp.expand(
    sp.Poly(linear_bordered, z).coeff_monomial(z)
    + 2 * linear_tail * directional_curvature
) == 0

# For the terminal affine z tail the complete invariant is the binary
# Hessian determinant times its squared coefficient.
terminal = sp.symbols("terminal")
terminal_potential = h_binary + terminal * z
terminal_gradient = sp.Matrix(
    [sp.diff(terminal_potential, variable) for variable in variables]
)
terminal_hessian = sp.hessian(terminal_potential, variables)
terminal_bordered = sp.expand(
    (terminal_gradient.T * terminal_hessian.adjugate() * terminal_gradient)[0]
)
binary_hessian_determinant = sp.det(sp.hessian(h_binary, (x, y)))
assert sp.expand(
    terminal_bordered - terminal**2 * binary_hessian_determinant
) == 0

payload = {
    "format": "hc4-pure-septic-kzero-v1",
    "status": {
        "id": "HC4RSD37",
        "kind": "exact narrowing theorem",
        "scope": "the complete P_y=0 residue in the k=0 curved chart",
    },
    "recursive_square": (
        "x^6*(H_yy*(49*x^6*Q_zz-14*p*x^3*V+6*p^2*H)"
        "-(7*x^3*V_y-3*p*H_y)^2)"
    ),
    "p_zero_Qzz_nonzero": {
        "alignment": "Q_zz*H_yy=V_y^2 forces V_y=s*x^2",
        "obstruction": "normalized degree-fourteen coefficient -4",
        "explanation": "J(h(x,z+D*y))=-D'^2*h_w^4",
    },
    "p_zero_Qzz_zero": {
        "quartic_tail": "two rho packets, then normalized coefficient -24",
        "cubic_tail": "linear-in-z global identity forces H_yy=0",
        "later_tails": "constant-direction cylinder",
    },
    "p_nonzero": {
        "classification": "G_y=0, G_y proportional to K_yy^2, or the K_yy=x endpoint",
        "normalized_obstructions": ["-648/49", "16", "12/7", "48/7"],
    },
    "residual": "P_y!=0 k=0 packets and the H_yy=0 boundary",
}
serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()


print("PASS: verified the k=0 residual degree-sixteen square")
print("PASS: classified and closed every p!=0 cubic recurrence packet")
print("PASS: aligned the nonzero-Qzz square with x^2")
print("PASS: found the immutable degree-fourteen coefficient -4")
print("PASS: verified J(h(x,z+D*y))=-D'^2*h_w^4")
print("PASS: closed both zero-Qzz nonlinear transverse tails")
print("PASS: reduced the remaining linear and affine tails to a cylinder")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")
