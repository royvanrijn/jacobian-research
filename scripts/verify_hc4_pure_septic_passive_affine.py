#!/usr/bin/env python3
"""Open the final passive-affine pure-septic boundary exactly."""

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
    / "hc4_pure_septic_passive_affine.json"
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


# Normalize c7=x^7 and the nonzero passive-affine sextic correction to
# c6=x^5*y.  For an arbitrary quintic C5 the degree-eighteen face is the
# passive Hessian determinant after one exact square completion.
c5_coefficients = sp.symbols("c0:21")
c5_monomials = [
    x**i * y**j * z ** (5 - i - j)
    for i in range(6)
    for j in range(6 - i)
]
C5 = sum(
    value * monomial
    for value, monomial in zip(c5_coefficients, c5_monomials, strict=True)
)
opening_components = {7: x**7, 6: x**5 * y, 5: C5}
degree_eighteen = bordered_face(opening_components, 24)
C5_hat = C5 - sp.Rational(2, 7) * x**3 * y**2
expected_degree_eighteen = 49 * x**12 * (
    sp.diff(C5_hat, y, 2) * sp.diff(C5_hat, z, 2)
    - sp.diff(C5_hat, y, z) ** 2
)
assert sp.expand(degree_eighteen - expected_degree_eighteen) == 0


# On the curved singular-passive-Hessian chart, use passive coordinates in
# which C5_hat=H5(x,y)+k*x^4*z.  Relative to the sextic linear form there
# are only the aligned chart L=y and the misaligned chart L=z.
h = sp.symbols("h0:6")
H5 = sum(h[index] * x ** (5 - index) * y**index for index in range(6))
k = sp.symbols("k")
c4_coefficients = sp.symbols("q0:15")
c4_monomials = [
    x**i * y**j * z ** (4 - i - j)
    for i in range(5)
    for j in range(5 - i)
]
C4 = sum(
    value * monomial
    for value, monomial in zip(c4_coefficients, c4_monomials, strict=True)
)


def degree_seventeen(linear_form: sp.Expr) -> sp.Expr:
    return bordered_face(
        {
            7: x**7,
            6: x**5 * linear_form,
            5: sp.Rational(2, 7) * x**3 * linear_form**2
            + H5
            + k * x**4 * z,
            4: C4,
        },
        23,
    )


H5_yy = sp.diff(H5, y, 2)
aligned_degree_seventeen = degree_seventeen(y)
assert sp.expand(
    aligned_degree_seventeen - 49 * x**12 * H5_yy * sp.diff(C4, z, 2)
) == 0

misaligned_forced_C4 = (
    sum(
        sp.symbols("u0:5")[index] * x ** (4 - index) * y**index
        for index in range(5)
    )
    + z
    * sum(
        sp.symbols("v0:4")[index] * x ** (3 - index) * y**index
        for index in range(4)
    )
    + sp.Rational(3, 7) * k * x**2 * z**2
    + sp.Rational(1, 49) * x * z**3
)
misaligned_degree_seventeen = degree_seventeen(z)
misaligned_schur_factor = (
    sp.diff(C4, z, 2)
    - sp.Rational(6, 7) * k * x**2
    - sp.Rational(6, 49) * x * z
)
assert sp.expand(
    misaligned_degree_seventeen
    - 49 * x**12 * H5_yy * misaligned_schur_factor
) == 0
assert sp.expand(
    sp.diff(misaligned_forced_C4, z, 2)
    - sp.Rational(6, 7) * k * x**2
    - sp.Rational(6, 49) * x * z
) == 0


# The aligned degree-sixteen descendant is another exact square recurrence.
u = sp.symbols("u0:5")
v = sp.symbols("v0:4")
U4 = sum(u[index] * x ** (4 - index) * y**index for index in range(5))
V3 = sum(v[index] * x ** (3 - index) * y**index for index in range(4))
t = sp.symbols("t0:10")
c3_monomials = [
    x**i * y**j * z ** (3 - i - j)
    for i in range(4)
    for j in range(4 - i)
]
C3 = sum(
    value * monomial for value, monomial in zip(t, c3_monomials, strict=True)
)
aligned_components = {
    7: x**7,
    6: x**5 * y,
    5: sp.Rational(2, 7) * x**3 * y**2 + H5 + k * x**4 * z,
    4: U4 + z * V3,
    3: C3,
}
aligned_degree_sixteen = bordered_face(aligned_components, 22)
E1 = t[4] * x + t[1] * y
aligned_expected = x**12 * (
    294 * t[0] * z * H5_yy
    + 14 * H5_yy * (7 * E1 - k**2 * x)
    - (7 * sp.diff(V3, y) - 3 * k * x**2) ** 2
)
assert sp.expand(aligned_degree_sixteen - aligned_expected) == 0


# In the misaligned chart the same face first removes the y^5 and y^4
# coefficients of H5, then leaves two exact ratios for v3 when h3!=0.
misaligned_components = {
    7: x**7,
    6: x**5 * z,
    5: sp.Rational(2, 7) * x**3 * z**2 + H5 + k * x**4 * z,
    4: misaligned_forced_C4,
    3: C3,
}
misaligned_degree_sixteen = bordered_face(misaligned_components, 22)
assert coefficient(misaligned_degree_sixteen, x**8 * y**8) == -25 * h[5] ** 2
assert sp.factor(
    coefficient(misaligned_degree_sixteen, x**10 * y**6).subs(h[5], 0)
) == -24 * h[4] ** 2
misaligned_ratio_equation = sp.factor(
    coefficient(misaligned_degree_sixteen, x**12 * y**4).subs(
        {h[4]: 0, h[5]: 0}
    )
)
assert sp.expand(
    misaligned_ratio_equation
    + 21 * (h[3] - 3 * v[3]) * (h[3] - 7 * v[3])
) == 0

# Normalize the curved misaligned endpoint further by h3=1, h2=k=0.
# The remaining degree-sixteen coefficients give v1=3*h1/7 and the
# following triangular solutions in the two v3 ratios.  Descend both to
# degree fifteen while retaining an arbitrary quadratic correction.
b = sp.symbols("b0:6")
Q2 = sum(
    value * monomial
    for value, monomial in zip(
        b, (z**2, y * z, y**2, x * z, x * y, x**2), strict=True
    )
)
misaligned_normalization = {
    h[2]: 0,
    h[3]: 1,
    h[4]: 0,
    h[5]: 0,
    k: 0,
    v[1]: sp.Rational(3, 7) * h[1],
    t[0]: 0,
    t[4]: (14 * v[0] - 5 * h[0]) / 49,
}
misaligned_ratio_A = {
    **misaligned_normalization,
    v[3]: sp.Rational(1, 3),
    v[2]: 0,
    t[1]: h[1] / 49,
}
misaligned_ratio_B = {
    **misaligned_normalization,
    v[3]: sp.Rational(1, 7),
    t[1]: h[1] / 49 + v[2] ** 2 / 3,
}
misaligned_constants: dict[str, list[tuple[tuple[int, int, int], sp.Expr]]] = {}
for ratio_name, ratio_substitution in (
    ("one_third", misaligned_ratio_A),
    ("one_seventh", misaligned_ratio_B),
):
    assert sp.expand(misaligned_degree_sixteen.subs(ratio_substitution)) == 0
    degree_fifteen_ratio = sp.Poly(
        sp.expand(
            bordered_face({**misaligned_components, 2: Q2}, 21).subs(
                ratio_substitution
            )
        ),
        x,
        y,
        z,
    )
    constants = [
        (monomial, sp.factor(value))
        for monomial, value in degree_fifteen_ratio.terms()
        if not value.free_symbols and value != 0
    ]
    misaligned_constants[ratio_name] = constants

assert misaligned_constants["one_third"] == [
    ((10, 4, 1), -sp.Rational(20, 21))
]
assert misaligned_constants["one_seventh"] == [
    ((10, 4, 1), sp.Rational(36, 7))
]


# Classify the aligned square globally.  Since F*Q=G^2 with degrees
# (3,1,2), unique factorization gives, up to nonzero scalars,
#
#     F=L*E^2, G=L*E, Q=L.
#
# Exhaust the same five ordered-line configurations relative to x and
# inspect their complete degree-fifteen descendants.
aligned_line_charts = {
    "Lx_Ey": (x, y),
    "Ly_Ex": (y, x),
    "finite_distinct": (y - x, y),
    "finite_equal": (y, y),
    "x_equal": (x, x),
}
aligned_chart_faces: dict[str, sp.Poly] = {}
aligned_chart_components: dict[str, dict[int, sp.Expr]] = {}
aligned_chart_special: dict[
    str, list[tuple[tuple[int, int, int], sp.Expr]]
] = {}
for chart_name, (L, E) in aligned_line_charts.items():
    chart_F = sp.expand(L * E**2)
    chart_G = sp.expand(L * E)
    chart_H = (
        sp.integrate(sp.integrate(chart_F, y), y)
        + h[0] * x**5
        + h[1] * x**4 * y
    )
    chart_V = (
        sp.integrate((chart_G + 3 * k * x**2) / 7, y)
        + v[0] * x**3
    )
    chart_E1 = L / 98 + k**2 * x / 7
    chart_C3 = (
        t[9] * x**3
        + t[8] * x**2 * y
        + t[6] * x * y**2
        + t[3] * y**3
        + z * (t[7] * x**2 + t[5] * x * y + t[2] * y**2)
        + z**2 * chart_E1
    )
    chart_components = {
        7: x**7,
        6: x**5 * y,
        5: sp.Rational(2, 7) * x**3 * y**2 + chart_H + k * x**4 * z,
        4: U4 + z * chart_V,
        3: chart_C3,
        2: Q2,
    }
    aligned_chart_components[chart_name] = chart_components
    aligned_chart_faces[chart_name] = sp.Poly(
        bordered_face(chart_components, 21), x, y, z
    )
    special = [
        (monomial, sp.factor(value))
        for monomial, value in aligned_chart_faces[chart_name].terms()
        if value != 0 and value.free_symbols <= {k}
    ]
    aligned_chart_special[chart_name] = special

assert ((14, 0, 1), sp.Rational(1, 7)) in aligned_chart_special["Lx_Ey"]
assert ((13, 1, 1), -sp.Rational(1, 7)) in aligned_chart_special["Ly_Ex"]
assert ((14, 0, 1), sp.Rational(1, 7)) in aligned_chart_special[
    "finite_distinct"
]
assert ((10, 5, 0), -sp.Rational(1, 42)) in aligned_chart_special[
    "finite_equal"
]
assert aligned_chart_special["x_equal"] == []

xeq_degree_fifteen_solution = {
    u[4]: 0,
    t[2]: (63 * k + 147 * u[3] - 10) / 686,
    b[0]: (
        -14 * h[1] * k
        + 3 * h[1]
        + 49 * k * v[0]
        + 49 * t[5]
        - 7 * u[2]
        - 14 * v[0]
    )
    / 343,
}
assert sp.expand(
    aligned_chart_faces["x_equal"].as_expr().subs(xeq_degree_fifteen_solution)
) == 0
ell0, ell1, ell2 = sp.symbols("ell0 ell1 ell2")
C1 = ell0 * z + ell1 * y + ell2 * x
xeq_degree_fourteen = sp.Poly(
    sp.expand(
        bordered_face(
            {**aligned_chart_components["x_equal"], 1: C1}, 20
        ).subs(xeq_degree_fifteen_solution)
    ),
    x,
    y,
    z,
)
xeq_u3 = sp.Rational(8, 49) - k
xeq_t3_roots = sp.solve(
    coefficient(xeq_degree_fourteen.as_expr(), x**13 * y).subs(u[3], xeq_u3),
    t[3],
)
assert len(xeq_t3_roots) == 1
xeq_t3 = sp.factor(xeq_t3_roots[0])
xeq_b1_roots = sp.solve(
    coefficient(xeq_degree_fourteen.as_expr(), x**14)
    .subs(u[3], xeq_u3)
    .subs(t[3], xeq_t3),
    b[1],
)
assert len(xeq_b1_roots) == 1
xeq_b1 = sp.factor(xeq_b1_roots[0])
xeq_degree_fourteen_solution = {u[3]: xeq_u3, t[3]: xeq_t3, b[1]: xeq_b1}
xeq_resonance = 490 * k**2 - 35 * k + 1
assert sp.expand(
    coefficient(xeq_degree_fourteen.as_expr(), x**12 * y**2).subs(
        xeq_degree_fourteen_solution
    )
    + xeq_resonance / 49
) == 0

xeq_degree_thirteen = sp.Poly(
    sp.expand(
        bordered_face(
            {**aligned_chart_components["x_equal"], 1: C1}, 19
        )
        .subs(xeq_degree_fifteen_solution)
        .subs(xeq_degree_fourteen_solution)
    ),
    x,
    y,
    z,
)
xeq_resonance_polynomial = sp.Poly(xeq_resonance, k)


def reduce_xeq_resonance(expression: sp.Expr) -> sp.Expr:
    numerator, denominator = sp.cancel(expression).as_numer_denom()
    inverse = sp.invert(sp.Poly(denominator, k), xeq_resonance_polynomial)
    return sp.factor(
        sp.rem(
            sp.Poly(numerator, k) * inverse,
            xeq_resonance_polynomial,
        ).as_expr()
    )


xeq_linear_obstruction = reduce_xeq_resonance(
    coefficient(xeq_degree_thirteen.as_expr(), x**11 * y * z)
)
assert sp.expand(
    xeq_linear_obstruction
    + 3 * (1449 * k - 125) / 24010
) == 0
assert sp.resultant(xeq_resonance, 1449 * k - 125, k) != 0


# The sole residue has two constant passive linear forms:
# c6=x^5*L and c5=2*x^3*L^2/7+x^4*M.  Split by their rank and inspect the
# next complete face with arbitrary c4.
tower_charts = {
    "independent": (y, z),
    "dependent": (y, 0),
    "L_zero": (0, y),
    "both_zero": (0, 0),
}
tower_faces: dict[str, sp.Expr] = {}
for tower_name, (tower_L, tower_M) in tower_charts.items():
    tower_components = {
        7: x**7,
        6: x**5 * tower_L,
        5: sp.Rational(2, 7) * x**3 * tower_L**2 + x**4 * tower_M,
        4: C4,
        3: C3,
    }
    tower_faces[tower_name] = bordered_face(tower_components, 22)

tower_shifts = {
    "independent": x * y**3 / 49 + sp.Rational(3, 7) * x**2 * y * z,
    "dependent": x * y**3 / 49,
    "L_zero": 0,
    "both_zero": 0,
}
for tower_name, shift in tower_shifts.items():
    shifted_quartic = C4 - shift
    passive_hessian_determinant = (
        sp.diff(shifted_quartic, y, 2) * sp.diff(shifted_quartic, z, 2)
        - sp.diff(shifted_quartic, y, z) ** 2
    )
    assert sp.expand(
        tower_faces[tower_name]
        - 49 * x**12 * passive_hessian_determinant
    ) == 0


# Passive singular-Hessian classification makes C4-shift a binary quartic
# Q4(x,N), plus x^3 times the complementary passive linear form.  Descend
# the finite relative directions of N against L and M.
A = sp.symbols("A0:5")


def binary_quartic(passive_form: sp.Expr) -> sp.Expr:
    return sum(
        A[index] * x ** (4 - index) * passive_form**index
        for index in range(5)
    )


tower_direction_charts = {
    "independent_Ny": (y, z, y),
    "independent_Nz": (y, z, z),
    "independent_Nsum": (y, z, y + z),
    "dependent_Ny": (y, 0, y),
    "dependent_Nz": (y, 0, z),
    "Lzero_Ny": (0, y, y),
    "Lzero_Nz": (0, y, z),
    "bothzero_Ny": (0, 0, y),
}
tower_degree_fifteen: dict[str, sp.Poly] = {}
tower_positive_z_counts: dict[str, int] = {}
for chart_name, (tower_L, tower_M, cylinder_N) in tower_direction_charts.items():
    base_name = (
        "independent"
        if tower_L == y and tower_M == z
        else "dependent"
        if tower_L == y
        else "L_zero"
        if tower_M == y
        else "both_zero"
    )
    transverse_coefficient = sp.symbols(f"s_{chart_name}")
    complementary_form = (
        z
        if cylinder_N == y
        else y
        if cylinder_N == z
        else z - y
    )
    chart_components = {
        7: x**7,
        6: x**5 * tower_L,
        5: sp.Rational(2, 7) * x**3 * tower_L**2 + x**4 * tower_M,
        4: tower_shifts[base_name]
        + binary_quartic(cylinder_N)
        + transverse_coefficient * x**3 * complementary_form,
        3: C3,
        2: Q2,
    }
    tower_degree_fifteen[chart_name] = sp.Poly(
        bordered_face(chart_components, 21), x, y, z
    )
    tower_positive_z_counts[chart_name] = sum(
        monomial[2] > 0
        for monomial, value in tower_degree_fifteen[chart_name].terms()
        if value != 0
    )

assert tower_positive_z_counts == {
    "independent_Ny": 3,
    "independent_Nz": 5,
    "independent_Nsum": 6,
    "dependent_Ny": 3,
    "dependent_Nz": 5,
    "Lzero_Ny": 3,
    "Lzero_Nz": 5,
    "bothzero_Ny": 3,
}


payload = {
    "format": "hc4-pure-septic-passive-affine-v1",
    "status": {
        "id": "HC4RSD39",
        "kind": "exact narrowing theorem",
        "scope": "the final H_yy=0 pure-septic boundary",
    },
    "degree_18": (
        "49*x^12*det Hess_(y,z)(c5-(2/7)*x^3*L^2)"
    ),
    "curved_shifted_quintic": {
        "aligned": (
            "c4_zz=0; the exact factor-minus-square recurrence has five "
            "ordered-line charts, all closed by degree 15 or 13"
        ),
        "misaligned": (
            "c4_zz=6*k*x^2/7+6*x*z/49; then h5=h4=0 and "
            "v3=h3/3 or h3/7, killed by -20/21 and 36/7"
        ),
    },
    "two_linear_form_tower": {
        "degree_16": (
            "49*x^12*det Hess_(y,z)(c4-shift), with shifts 0, "
            "x*y^3/49, or x*y^3/49+3*x^2*y*z/7"
        ),
        "quartic_packets": list(tower_direction_charts),
        "degree_15_positive_z_counts": tower_positive_z_counts,
    },
    "residual": (
        "eight explicit degree-fifteen quartic-cylinder direction packets"
    ),
}
serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: verified the shifted passive-Hessian opening")
print("PASS: split the curved shifted-quintic chart by alignment")
print("PASS: killed both misaligned ratios")
print("PASS: killed all five aligned ordered-line charts")
print("PASS: reduced the two-linear-form tower to eight quartic packets")
print("NARROWING: eight explicit degree-fifteen packets remain")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")
