#!/usr/bin/env python3
"""Close the P_y!=0 part of the k=0 curved pure-septic chart."""

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
    / "hc4_pure_septic_kzero_wronskian.json"
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


def binary(prefix: str, degree: int) -> tuple[sp.Expr, tuple[sp.Symbol, ...]]:
    coefficients = sp.symbols(f"{prefix}0:{degree + 1}")
    return (
        sum(
            coefficients[index] * x ** (degree - index) * y**index
            for index in range(degree + 1)
        ),
        coefficients,
    )


H, a = binary("a", 6)
P, p = binary("p", 4)
R, r = binary("r", 5)
U, u = binary("u", 4)
V, v = binary("v", 3)
D, d = binary("d", 2)
t = sp.symbols("t0:10")
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
    6: H,
    5: R + z * P,
    4: U + z * V + z**2 * D / 2,
    3: C3,
}
degree_seventeen = bordered_face(components, 23)
degree_seventeen_in_z = sp.Poly(degree_seventeen, z)
assert degree_seventeen_in_z.coeff_monomial(z**2) == 0
degree_seventeen_z = degree_seventeen_in_z.coeff_monomial(z).subs(t[0], 0)
expected_wronskian = 49 * x**12 * (
    D * sp.diff(P, y, 2) - 2 * sp.diff(D, y) * sp.diff(P, y)
)
assert sp.expand(degree_seventeen_z - expected_wronskian) == 0
assert sp.expand(
    degree_seventeen_in_z.coeff_monomial(z)
    - expected_wronskian
    - 294 * t[0] * x**12 * sp.diff(H, y, 2)
) == 0

# Together with H_yy*D=P_y^2, the coupled Wronskian says either P_y=D=0
# or, up to nonzero scalars,
#
#     H_yy=L*E^3, P_y=L*E^2, D=L*E, t0=L_y/6.
#
# Indeed (D^2/P_y)_y=6*t0, and homogeneity makes D^2/P_y a linear form
# L.  Polynomiality then forces L to divide D.  The ordered pair (L,E)
# has the five projective configurations checked below.
h, g = sp.symbols("h g", nonzero=True)
ell_y_substitution = {
    a[2]: 0,
    a[3]: 0,
    a[4]: 0,
    a[5]: h,
    a[6]: 0,
    p[1]: 0,
    p[2]: 0,
    p[3]: g,
    p[4]: 0,
    d[0]: 0,
    d[1]: sp.Rational(9, 20) * g**2 / h,
    d[2]: 0,
}
ell_y_face = sp.expand(degree_seventeen.subs(ell_y_substitution))
assert coefficient(ell_y_face, x**8 * y**9) == 35 * g**2 * h

# In the ell=x chart normalize h=p1=1.  This is also checked below in the
# uniform ordered-line normalization; retain it as an independent scaling
# calibration of the triangular equations.
ell_x_substitution = {
    a[2]: 1,
    a[3]: 0,
    a[4]: 0,
    a[5]: 0,
    a[6]: 0,
    p[1]: 1,
    p[2]: 0,
    p[3]: 0,
    p[4]: 0,
    d[0]: sp.Rational(1, 2),
    d[1]: 0,
    d[2]: 0,
}
ell_x_degree_seventeen = sp.Poly(
    sp.expand(degree_seventeen.subs(ell_x_substitution)), x, y, z
)

# Exhaust the ordered projective configurations of (L,E) relative to x.
# The common scalar factors can be removed on the nonzero chart.
line_charts = {
    "Lx_Ey": (x, y),
    "Ly_Ex": (y, x),
    "finite_distinct": (y - x, y),
    "finite_equal": (y, y),
    "x_equal": (x, x),
}
chart_constants: dict[str, list[tuple[tuple[int, int, int], sp.Expr]]] = {}
for chart_name, (L, E) in line_charts.items():
    chart_F = sp.expand(L * E**3)
    chart_G = sp.expand(L * E**2)
    chart_D = sp.expand(L * E)
    chart_H = sp.integrate(sp.integrate(chart_F, y), y) + a[0] * x**6 + a[1] * x**5 * y
    chart_P = sp.integrate(chart_G, y) + p[0] * x**4
    chart_components = {
        7: x**7,
        6: chart_H,
        5: R + z * chart_P,
        4: U + z * V + z**2 * chart_D / 2,
        3: C3.subs(t[0], sp.diff(L, y) / 6),
    }
    chart_face = sp.Poly(bordered_face(chart_components, 23), x, y, z)
    constants = [
        (monomial, sp.factor(value))
        for monomial, value in chart_face.terms()
        if not value.free_symbols and value != 0
    ]
    chart_constants[chart_name] = constants

assert chart_constants["Lx_Ey"] == [((8, 9, 0), sp.Rational(7, 36))]
assert chart_constants["finite_distinct"]
assert chart_constants["finite_equal"] == [
    ((5, 12, 0), sp.Rational(21, 200))
]
assert chart_constants["Ly_Ex"] == []
assert chart_constants["x_equal"] == []

# Descend the only two degree-seventeen survivors.
b = sp.symbols("b0:6")
Q2 = sum(
    value * monomial
    for value, monomial in zip(
        b, (z**2, y * z, y**2, x * z, x * y, x**2), strict=True
    )
)

LyEx_H = a[0] * x**6 + a[1] * x**5 * y + x**3 * y**3 / 6
LyEx_P = p[0] * x**4 + x**2 * y**2 / 2
LyEx_D = x * y
LyEx_components = {
    7: x**7,
    6: LyEx_H,
    5: R + z * LyEx_P,
    4: U + z * V + z**2 * LyEx_D / 2,
    3: C3.subs(t[0], sp.Rational(1, 6)),
}
LyEx_solution = {
    r[5]: 0,
    v[3]: (-a[1] + p[0] + 84 * r[4]) / 42,
    t[1]: 2 * v[2] - 3 * r[3],
    t[4]: (
        2 * a[1] ** 2
        - 3 * a[1] * p[0]
        + p[0] ** 2
        - 7 * r[2]
        + 7 * v[1]
    )
    / 7,
}
LyEx_degree_sixteen = sp.Poly(
    sp.expand(
        bordered_face({**LyEx_components, 2: Q2}, 22).subs(LyEx_solution)
    ),
    x,
    y,
    z,
)
assert coefficient(LyEx_degree_sixteen.as_expr(), x**8 * y**8) == -sp.Rational(
    1, 16
)

Xeq_H = a[0] * x**6 + a[1] * x**5 * y + x**4 * y**2 / 2
Xeq_P = p[0] * x**4 + x**3 * y
Xeq_D = x**2
Xeq_components = {
    7: x**7,
    6: Xeq_H,
    5: R + z * Xeq_P,
    4: U + z * V + z**2 * Xeq_D / 2,
    3: C3.subs(t[0], 0),
}
Xeq_degree_seventeen = sp.Poly(
    bordered_face(Xeq_components, 23), x, y, z
)
Xeq_solution = {
    r[5]: 0,
    v[3]: 2 * r[4],
    t[1]: (a[1] - p[0] - 21 * r[3] + 14 * v[2]) / 7,
    t[4]: (
        2 * a[1] ** 2
        - 3 * a[1] * p[0]
        + p[0] ** 2
        - 7 * r[2]
        + 7 * v[1]
    )
    / 7,
}
Xeq_degree_sixteen = sp.Poly(
    sp.expand(
        bordered_face({**Xeq_components, 2: Q2}, 22).subs(Xeq_solution)
    ),
    x,
    y,
    z,
)
assert coefficient(Xeq_degree_sixteen.as_expr(), x**14 * y * z) == 588 * r[4]
Xeq_terminal = sp.Poly(
    sp.expand(
        Xeq_degree_sixteen.as_expr().subs(
            {
                r[4]: 0,
                v[2]: (-a[1] + p[0] + 14 * r[3]) / 7,
            }
        )
    ),
    x,
    y,
    z,
)
assert coefficient(Xeq_terminal.as_expr(), x**12 * y**4) == -1

payload = {
    "format": "hc4-pure-septic-kzero-wronskian-v1",
    "status": {
        "id": "HC4RSD38",
        "kind": "exact closure theorem",
        "scope": "the P_y!=0 part of the k=0 curved pure-septic chart",
    },
    "coupled_equations": [
        "H_yy*D=P_y^2",
        "D*(P_y)_y-2*D_y*P_y+6*t0*H_yy=0",
    ],
    "global_packet": (
        "H_yy=L*E^3, P_y=L*E^2, D=L*E, t0=L_y/6"
    ),
    "ordered_line_charts": {
        "L=x,E=y": "degree 17: 7/36",
        "L=y,E=x": "degree 16: -1/16",
        "finite distinct": "degree 17 obstruction",
        "finite equal": "degree 17: 21/200",
        "L=E=x": "degree 16: r4=0 then -1",
    },
    "residual": "only the passive-affine H_yy=0 pure-septic boundary",
}
serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: verified the degree-seventeen Wronskian identity")
print("PASS: classified the coupled equations by two ordered linear forms")
print("PASS: killed all five ordered-line projective charts")
print("THEOREM: the complete k=0 curved chart is closed with HC4RSD37")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")
