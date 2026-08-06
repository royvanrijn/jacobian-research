#!/usr/bin/env python3
"""Verify the degree-eighteen square obstruction over the pure septic top.

HC4RSD34 leaves, on its curved chart,

    c_7 = x^7,
    c_6 = H_6(x,y) + k*x^5*z,
    c_5 = R_5(x,y) + z*P_4(x,y) + (2/7)*k^2*x^3*z^2.

This checker retains every coefficient of the quartic correction.  It
extracts the degree-eighteen part of the bordered invariant directly from
its four homogeneous derivative factors, verifies the global square
identity, and records the forced z^4 and z^3 coefficients.
"""

from __future__ import annotations

import hashlib
import json
from itertools import product
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT / "artifacts" / "generated-results" / "hc4_pure_septic_degree18.json"
)

x, y, z = sp.symbols("x y z")
variables = (x, y, z)


def bordered_face(
    components: dict[int, sp.Expr], component_degree_sum: int
) -> sp.Expr:
    """Extract one homogeneous face of grad(c)^T adj(Hess(c)) grad(c).

    A bordered-invariant monomial has two gradient factors and two Hessian
    factors.  If their source components have degrees d1,...,d4, its
    ordinary degree is d1+...+d4-6.  Convolving only the requested source
    degree sum avoids expanding irrelevant faces.
    """

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
        return result

    # The six diagonal and six mixed terms of the bordered invariant.
    result = (
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
    return sp.expand(result)


a = sp.symbols("a0:7")
p = sp.symbols("p0:5")
r = sp.symbols("r0:6")
u = sp.symbols("u0:15")
k = sp.symbols("k")

H = sum(a[index] * x ** (6 - index) * y**index for index in range(7))
P = sum(p[index] * x ** (4 - index) * y**index for index in range(5))
R = sum(r[index] * x ** (5 - index) * y**index for index in range(6))
quartic_monomials = [
    x**i * y**j * z ** (4 - i - j)
    for i in range(5)
    for j in range(5 - i)
]
C = sum(
    coefficient * monomial
    for coefficient, monomial in zip(u, quartic_monomials, strict=True)
)

components = {
    7: x**7,
    6: H + k * x**5 * z,
    5: R + z * P + sp.Rational(2, 7) * k**2 * x**3 * z**2,
    4: C,
}

degree_eighteen = bordered_face(components, 24)
H_y = sp.diff(H, y)
H_yy = sp.diff(H, y, 2)
P_y = sp.diff(P, y)
square = 7 * x**2 * P_y - 4 * k * H_y
schur_factor = (
    49 * x**4 * sp.diff(C, z, 2)
    - 42 * k * x**2 * P
    + 18 * k**2 * H
    - 6 * k**3 * x**5 * z
)
expected = x**8 * (H_yy * schur_factor - square**2)
assert sp.expand(degree_eighteen - expected) == 0

# Split a homogeneous quartic according to its passive z-degree.  The z^2
# coefficient of the identity kills z^4.  Its z coefficient then fixes the
# complete z^3 linear form to k^3*x/49.
U_coefficients = sp.symbols("A0:5")
V_coefficients = sp.symbols("B0:4")
D_coefficients = sp.symbols("D0:3")
ell_x, ell_y, q = sp.symbols("ell_x ell_y q")
U = sum(
    U_coefficients[index] * x ** (4 - index) * y**index
    for index in range(5)
)
V = sum(
    V_coefficients[index] * x ** (3 - index) * y**index
    for index in range(4)
)
D = sum(
    D_coefficients[index] * x ** (2 - index) * y**index
    for index in range(3)
)
split_C = (
    U
    + z * V
    + z**2 * D / 2
    + z**3 * (ell_x * x + ell_y * y)
    + q * z**4
)
split_expression = sp.expand(
    H_yy
    * (
        49 * x**4 * sp.diff(split_C, z, 2)
        - 42 * k * x**2 * P
        + 18 * k**2 * H
        - 6 * k**3 * x**5 * z
    )
    - square**2
)
split_polynomial = sp.Poly(split_expression, z)
assert sp.factor(split_polynomial.coeff_monomial(z**2)) == 588 * q * x**4 * H_yy
assert sp.expand(
    split_polynomial.coeff_monomial(z)
    - 6 * x**4 * H_yy * (49 * ell_x * x + 49 * ell_y * y - k**3 * x)
) == 0

forced_C = sp.expand(
    U + z * V + z**2 * D / 2 + sp.Rational(1, 49) * k**3 * x * z**3
)
forced_expression = sp.expand(
    H_yy
    * (
        49 * x**4 * sp.diff(forced_C, z, 2)
        - 42 * k * x**2 * P
        + 18 * k**2 * H
        - 6 * k**3 * x**5 * z
    )
    - square**2
)
binary_identity = sp.expand(
    H_yy * (49 * x**4 * D - 42 * k * x**2 * P + 18 * k**2 * H)
    - square**2
)
assert sp.expand(forced_expression - binary_identity) == 0

# On the genuinely moving chart k!=0, the coefficients below the x^4*D
# threshold are immutable.  They align the passive curvature with the pure
# root x before any Groebner splitting is needed.
binary_polynomial = sp.Poly(binary_identity, x, y)
assert sp.factor(binary_polynomial.coeff_monomial(y**10)) == -36 * a[6] ** 2 * k**2
assert sp.factor(
    binary_polynomial.coeff_monomial(x**2 * y**8).subs(a[6], 0)
) == -40 * a[5] ** 2 * k**2
aligned_substitution = {a[5]: 0, a[6]: 0}
assert sp.expand(
    binary_polynomial.coeff_monomial(x**4 * y**6).subs(aligned_substitution)
    + 8 * (-5 * a[4] * k + 14 * p[4]) * (-a[4] * k + 7 * p[4])
) == 0

# The second nonzero-curvature ratio fixes the next coefficient.  At the
# a4=0 endpoint, the following square fixes it with the third ratio.
second_ratio = {
    **aligned_substitution,
    p[4]: sp.Rational(5, 14) * a[4] * k,
}
assert sp.expand(
    binary_polynomial.coeff_monomial(x**5 * y**5).subs(second_ratio)
    + 18 * a[4] * k * (-5 * a[3] * k + 14 * p[3])
) == 0
rank_drop = {**aligned_substitution, a[4]: 0, p[4]: 0}
assert sp.expand(
    binary_polynomial.coeff_monomial(x**6 * y**4).subs(rank_drop)
    + 9 * (-2 * a[3] * k + 7 * p[3]) ** 2
) == 0


def solve_triangular_D(
    substitution: dict[sp.Symbol, sp.Expr], start_degree: int
) -> dict[sp.Symbol, sp.Expr]:
    """Solve the three triangular coefficients for D2,D1,D0."""

    solution: dict[sp.Symbol, sp.Expr] = {}
    for x_degree, variable in zip(
        range(start_degree, start_degree + 3),
        reversed(D_coefficients),
        strict=True,
    ):
        equation = binary_polynomial.coeff_monomial(
            x**x_degree * y ** (10 - x_degree)
        ).subs(substitution).subs(solution)
        roots = sp.solve(equation, variable)
        assert len(roots) == 1
        solution[variable] = sp.factor(roots[0])
    return solution


# Complete the degree-eighteen classification on k*a4!=0.  Each curvature
# ratio first determines D2,D1,D0.  The last two coefficients split into a
# generic resonant component and the double-root discriminant
# Delta=8*a2*a4-3*a3^2.
normalized = {a[0]: 0, a[1]: 0, a[5]: 0, a[6]: 0}
a2, a3, a4 = a[2], a[3], a[4]
p0, p1, p2, p3, p4 = p

ratio_A = {**normalized, p4: a4 * k / 7}
solution_A = solve_triangular_D(ratio_A, 6)
residual_A9 = sp.factor(
    binary_polynomial.coeff_monomial(x**9 * y).subs(ratio_A).subs(solution_A)
)
residual_A10 = sp.factor(
    binary_polynomial.coeff_monomial(x**10).subs(ratio_A).subs(solution_A)
)
A_left = -12 * a2 * a4 * k + 6 * a3**2 * k - 21 * a3 * p3 + 28 * a4 * p2
A_right = (
    20 * a2 * a3 * a4 * k
    - 28 * a2 * a4 * p3
    - 6 * a3**3 * k
    + 21 * a3**2 * p3
    - 28 * a3 * a4 * p2
    + 56 * a4**2 * p1
)
assert sp.cancel(residual_A9 + A_left * A_right / (8 * a4**3)) == 0

A_generic_p2 = 3 * (4 * a2 * a4 * k - 2 * a3**2 * k + 7 * a3 * p3) / (
    28 * a4
)
A_generic_linear = 2 * a2 * a3 * k - 7 * a2 * p3 + 14 * a4 * p1
assert sp.cancel(
    residual_A10.subs(p2, A_generic_p2)
    + A_generic_linear**2 / (4 * a4**2)
) == 0

A_exceptional_p1 = (
    -20 * a2 * a3 * a4 * k
    + 28 * a2 * a4 * p3
    + 6 * a3**3 * k
    - 21 * a3**2 * p3
    + 28 * a3 * a4 * p2
) / (56 * a4**2)
discriminant = 8 * a2 * a4 - 3 * a3**2
assert sp.cancel(
    residual_A10.subs(p1, A_exceptional_p1)
    - discriminant * A_left**2 / (192 * a4**4)
) == 0

ratio_B = {
    **normalized,
    p4: sp.Rational(5, 14) * a4 * k,
    p3: sp.Rational(5, 14) * a3 * k,
}
solution_B = solve_triangular_D(ratio_B, 6)
residual_B9 = sp.factor(
    binary_polynomial.coeff_monomial(x**9 * y).subs(ratio_B).subs(solution_B)
)
residual_B10 = sp.factor(
    binary_polynomial.coeff_monomial(x**10).subs(ratio_B).subs(solution_B)
)
B_left = -28 * a2 * a4 * k + 3 * a3**2 * k + 56 * a4 * p2
B_right = 32 * a2 * a3 * a4 * k - 3 * a3**3 * k - 56 * a3 * a4 * p2 + 112 * a4**2 * p1
assert sp.cancel(residual_B9 + B_left * B_right / (32 * a4**3)) == 0

B_generic_p2 = k * (28 * a2 * a4 - 3 * a3**2) / (56 * a4)
B_generic_linear = a2 * a3 * k + 28 * a4 * p1
assert sp.cancel(
    residual_B10.subs(p2, B_generic_p2)
    + B_generic_linear**2 / (16 * a4**2)
) == 0

B_exceptional_p1 = a3 * (
    -32 * a2 * a4 * k + 3 * a3**2 * k + 56 * a4 * p2
) / (112 * a4**2)
assert sp.cancel(
    residual_B10.subs(p1, B_exceptional_p1)
    - discriminant * B_left**2 / (768 * a4**4)
) == 0

# At a4=0, the a3!=0 chart is triangular and ends in one square.  At
# a3=0, the pure-fourth curvature endpoint has a unique D for arbitrary
# p0,p1,p2 and no residual equation.
rank_three = {
    **normalized,
    a4: 0,
    p4: 0,
    p3: sp.Rational(2, 7) * a3 * k,
}
solution_rank_three = solve_triangular_D(rank_three, 7)
rank_three_residual = sp.factor(
    binary_polynomial.coeff_monomial(x**10)
    .subs(rank_three)
    .subs(solution_rank_three)
)
assert sp.cancel(
    rank_three_residual
    + (6 * a2**2 * k - 14 * a2 * p2 + 21 * a3 * p1) ** 2
    / (9 * a3**2)
) == 0

rank_four = {**rank_three, a3: 0, p3: 0}
solution_rank_four = solve_triangular_D(rank_four, 8)
assert sp.cancel(binary_identity.subs(rank_four).subs(solution_rank_four)) == 0


payload = {
    "format": "hc4-pure-septic-degree18-v1",
    "status": {
        "id": "HC4RSD35",
        "kind": "exact narrowing theorem",
        "scope": "degree-eighteen descendants of the curved pure-septic chart",
        "result": (
            "the complete face is one Hessian-factor times a Schur factor "
            "minus a square; it fixes the quartic z^4 and z^3 tails and "
            "leaves one binary divisibility equation"
        ),
    },
    "degree_18": (
        "x^8*(H_yy*(49*x^4*(c4)_zz-42*k*x^2*P+18*k^2*H"
        "-6*k^3*x^5*z)-(7*x^2*P_y-4*k*H_y)^2)"
    ),
    "forced_quartic": (
        "c4=U4(x,y)+z*V3(x,y)+z^2*D2(x,y)/2+k^3*x*z^3/49"
    ),
    "binary_identity": (
        "H_yy*(49*x^4*D2-42*k*x^2*P+18*k^2*H)"
        "=(7*x^2*P_y-4*k*H_y)^2"
    ),
    "moving_alignment": {
        "forced_curvature": "a6=a5=0, so H_yy is divisible by x^2",
        "nonzero_a4_ratios": ["p4=a4*k/7", "p4=5*a4*k/14"],
        "second_ratio_descendant": "p3=5*a3*k/14",
        "a4_zero_endpoint": "p4=0 and p3=2*a3*k/7",
    },
    "moving_packet_classification": {
        "a4_nonzero_ratio_A": [
            "generic resonance A_left=A_generic_linear=0",
            "double-root discriminant with A_right=0",
        ],
        "a4_nonzero_ratio_B": [
            "generic resonance B_left=B_generic_linear=0",
            "double-root discriminant with B_right=0",
        ],
        "a4_zero_a3_nonzero": "triangular D and one terminal square",
        "a4_zero_a3_zero": "pure-fourth curvature with unique D",
        "packet_count": 6,
    },
    "residual": (
        "six exact moving-curvature packets followed by degree seventeen, "
        "the k=0 curved chart, and the passive-affine H_yy=0 boundary"
    ),
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: verified the complete pure-septic degree-eighteen face")
print("PASS: exposed the global Hessian-factor minus square obstruction")
print("PASS: killed quartic z^4 and fixed quartic z^3 exactly")
print("PASS: aligned every nonzero-k curvature packet with the pure root")
print("PASS: decomposed both ratios and both rank-drop endpoints")
print("THEOREM: the moving chart reduces to six exact degree-eighteen packets")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")
