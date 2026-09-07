#!/usr/bin/env python3
"""Close the k!=0 curved pure-septic packets left by HC4RSD35."""

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
    ROOT
    / "artifacts"
    / "generated-results"
    / "hc4_pure_septic_moving_closure.json"
)


def bordered_face(
    components: dict[int, sp.Expr], component_degree_sum: int
) -> sp.Expr:
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
        for d1, d2, d3, d4 in product(components, repeat=4):
            if d1 + d2 + d3 + d4 != component_degree_sum:
                continue
            result += (
                gradients[d1][gradient_indices[0]]
                * gradients[d2][gradient_indices[1]]
                * hessians[d3][hessian_indices[0]]
                * hessians[d4][hessian_indices[1]]
            )
        return result

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


a = sp.symbols("a0:7")
p = sp.symbols("p0:5")
r = sp.symbols("r0:6")
u = sp.symbols("u0:5")
v = sp.symbols("v0:4")
d = sp.symbols("d0:3")
t = sp.symbols("t0:10")
k = sp.symbols("k")


def binary(coefficients: tuple[sp.Symbol, ...], degree: int) -> sp.Expr:
    return sum(
        coefficients[index] * x ** (degree - index) * y**index
        for index in range(degree + 1)
    )


H = binary(a, 6)
P = binary(p, 4)
R = binary(r, 5)
U = binary(u, 4)
V = binary(v, 3)
D = binary(d, 2)
cubic_monomials = [
    x**i * y**j * z ** (3 - i - j)
    for i in range(4)
    for j in range(4 - i)
]
C3 = sum(
    value * monomial for value, monomial in zip(t, cubic_monomials, strict=True)
)

components = {
    7: x**7,
    6: H + k * x**5 * z,
    5: R + z * P + sp.Rational(2, 7) * k**2 * x**3 * z**2,
    4: U
    + z * V
    + z**2 * D / 2
    + sp.Rational(1, 49) * k**3 * x * z**3,
    3: C3,
}
degree_seventeen = bordered_face(components, 23)

common_alignment = {a[5]: 0, a[6]: 0}
ratio_A = {**common_alignment, p[4]: a[4] * k / 7}
ratio_B = {**common_alignment, p[4]: 5 * a[4] * k / 14}
assert sp.factor(
    coefficient(degree_seventeen, x**7 * y**10).subs(ratio_A)
) == sp.Rational(72, 7) * a[4] ** 3 * k**2
assert sp.factor(
    coefficient(degree_seventeen, x**7 * y**10).subs(ratio_B)
) == sp.Rational(18, 7) * a[4] ** 3 * k**2

rank_three = {
    **common_alignment,
    a[4]: 0,
    p[4]: 0,
    p[3]: sp.Rational(2, 7) * a[3] * k,
}
rank_three_expected = sp.Rational(24, 7) * a[3] ** 3 * k**2
rank_three_hits = [
    monomial
    for monomial, value in sp.Poly(degree_seventeen, x, y, z).terms()
    if sp.factor(value.subs(rank_three)) == rank_three_expected
]
assert rank_three_hits == [(10, 7, 0)]


# Pure-x^4 curvature endpoint.  Solve its unique degree-eighteen D first.
endpoint = {
    **rank_three,
    a[0]: 0,
    a[1]: 0,
    a[3]: 0,
    p[3]: 0,
}
H_endpoint = H.subs(endpoint)
P_endpoint = P.subs(endpoint)
endpoint_identity = sp.expand(
    sp.diff(H_endpoint, y, 2)
    * (49 * x**4 * D - 42 * k * x**2 * P_endpoint + 18 * k**2 * H_endpoint)
    - (7 * x**2 * sp.diff(P_endpoint, y) - 4 * k * sp.diff(H_endpoint, y)) ** 2
)
D_solution: dict[sp.Symbol, sp.Expr] = {}
for monomial, variable in zip((x**8 * y**2, x**9 * y, x**10), (d[2], d[1], d[0]), strict=True):
    equation = coefficient(endpoint_identity.subs(D_solution), monomial)
    roots = sp.solve(equation, variable)
    assert len(roots) == 1
    D_solution[variable] = sp.factor(roots[0])
assert sp.expand(endpoint_identity.subs(D_solution)) == 0

endpoint_degree_seventeen = sp.expand(
    degree_seventeen.subs(endpoint).subs(D_solution)
)
endpoint_polynomial = sp.Poly(endpoint_degree_seventeen, x, y, z)

# These are the two independent z-coefficients and the terminal binary
# coefficient.  They leave exactly the resonances p2=3*a2*k/7 and
# p2=4*a2*k/7, with p1=0 on the second.
expected_z_coefficients = {
    sp.factor(
        -12 * p[1] * (-3 * a[2] * k + 7 * p[2]) ** 2 / a[2]
    ),
    sp.factor(
        -12
        * (-4 * a[2] * k + 7 * p[2])
        * (-3 * a[2] * k + 7 * p[2]) ** 2
        / (7 * a[2])
    ),
}
actual_z_coefficients = {
    sp.factor(value)
    for monomial, value in endpoint_polynomial.terms()
    if monomial[2] == 1 and sp.factor(value) in expected_z_coefficients
}
assert actual_z_coefficients == expected_z_coefficients

expected_binary = sp.factor(
    140
    * r[5]
    * (-3 * a[2] * k + 7 * p[2])
    * (-a[2] * k + 2 * p[2])
    / a[2]
)
assert expected_binary in {
    sp.factor(value)
    for monomial, value in endpoint_polynomial.terms()
    if monomial[2] == 0
}

# Normalize a2=k=1.  Solve every degree-seventeen equation in the two
# endpoint resonances and expose the immutable degree-sixteen coefficient.
normalized_endpoint = {a[2]: 1, k: 1}
quadratic_form = sum(
    sp.symbols("b0:6")[index] * monomial
    for index, monomial in enumerate((z**2, y * z, y**2, x * z, x * y, x**2))
)
components_with_quadratic = {**components, 2: quadratic_form}
degree_sixteen = bordered_face(components_with_quadratic, 22)


E4_polynomial = sp.Poly(
    endpoint_degree_seventeen.subs(normalized_endpoint).subs(p[2], sp.Rational(4, 7)).subs(p[1], 0),
    x,
    y,
    z,
)

E4_solution = {
    r[5]: 0,
    r[4]: sp.Rational(4, 35),
    v[3]: sp.Rational(5, 14) * r[3],
    v[2]: (4 * p[0] + 5 * r[2]) / 14,
    t[0]: p[0] / 49,
    t[1]: (14 * v[1] - 5 * r[1]) / 49,
    t[4]: (7 * p[0] ** 2 - 5 * r[0] + 14 * v[0]) / 49,
}
assert sp.expand(E4_polynomial.as_expr().subs(E4_solution)) == 0
E4_degree_sixteen = sp.Poly(
    sp.expand(
        degree_sixteen.subs(endpoint)
        .subs(D_solution)
        .subs(normalized_endpoint)
        .subs({p[2]: sp.Rational(4, 7), p[1]: 0})
        .subs(E4_solution)
    ),
    x,
    y,
    z,
)
assert coefficient(E4_degree_sixteen.as_expr(), x**10 * y**6) == sp.Rational(
    256, 1225
)

E3_polynomial = sp.Poly(
    endpoint_degree_seventeen.subs(normalized_endpoint).subs(
        p[2], sp.Rational(3, 7)
    ),
    x,
    y,
    z,
)
E3_solution = {
    r[4]: (2 - 245 * p[1] * r[5]) / 14,
    v[3]: (2 * p[1] + 2 * r[3] - 245 * p[1] ** 2 * r[5]) / 14,
    t[0]: (4 * p[0] + 7 * p[1] ** 2) / 196,
    t[1]: -(
        28 * p[0] * p[1]
        + 147 * p[1] ** 2 * r[3]
        + 56 * p[1] * r[2]
        - 196 * p[1] * v[2]
        + 8 * r[1]
        - 28 * v[1]
    )
    / 196,
    t[4]: (
        28 * p[0] ** 2
        - 49 * p[1] ** 2 * r[2]
        - 42 * p[1] * r[1]
        + 98 * p[1] * v[1]
        - 20 * r[0]
        + 56 * v[0]
    )
    / 196,
}
assert sp.expand(E3_polynomial.as_expr().subs(E3_solution)) == 0
E3_degree_sixteen = sp.Poly(
    sp.expand(
        degree_sixteen.subs(endpoint)
        .subs(D_solution)
        .subs(normalized_endpoint)
        .subs(p[2], sp.Rational(3, 7))
        .subs(E3_solution)
    ),
    x,
    y,
    z,
)
assert coefficient(E3_degree_sixteen.as_expr(), x**11 * y**4 * z) == -sp.Rational(
    24, 49
)

payload = {
    "format": "hc4-pure-septic-moving-closure-v1",
    "status": {
        "id": "HC4RSD36",
        "kind": "exact closure theorem",
        "scope": "the six k!=0 curved packets left by HC4RSD35",
    },
    "degree_17_obstructions": {
        "a4_ratio_A": "72*a4^3*k^2/7",
        "a4_ratio_B": "18*a4^3*k^2/7",
        "a4_zero_a3_nonzero": "24*a3^3*k^2/7",
    },
    "pure_x4_resonances": {
        "ratios": ["p2=3*a2*k/7", "p2=4*a2*k/7 with p1=0"],
        "normalized_degree_16_obstructions": ["-24/49", "256/1225"],
    },
    "residual": "k=0 curved chart and passive-affine H_yy=0 boundary",
}
serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: killed all four nonzero-a4 packets at degree seventeen")
print("PASS: killed the a4=0,a3!=0 packet at degree seventeen")
print("PASS: killed both pure-x4 endpoint resonances at degree sixteen")
print("THEOREM: all six nonzero-k packets from HC4RSD35 are impossible")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")
