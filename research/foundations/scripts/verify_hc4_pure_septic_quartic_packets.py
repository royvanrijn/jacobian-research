#!/usr/bin/env python3
"""Close the eight passive-affine pure-septic quartic packets exactly.

This continues HC4RSD39.  The eight nominal direction charts share one
quartic-polar rank-one equation.  Its zero stratum is either impossible or
a fixed cylinder, while its square-Hessian resonance has an immutable lower
face.  All calculations are over QQ and remain valid over every
characteristic-zero field.
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
    ROOT
    / "artifacts"
    / "generated-results"
    / "hc4_pure_septic_quartic_packets.json"
)


def bordered_face(
    components: dict[int, sp.Expr], component_degree_sum: int
) -> sp.Expr:
    """Extract one face of grad(c)^T*adj(Hess(c))*grad(c)."""

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


def face(
    components: dict[int, sp.Expr],
    component_degree_sum: int,
    substitutions: dict[sp.Symbol, sp.Expr],
) -> sp.Expr:
    return sp.expand(bordered_face(components, component_degree_sum).subs(substitutions))


def unique_solution(expression: sp.Expr, variable: sp.Symbol) -> sp.Expr:
    solutions = sp.solve(expression, variable)
    assert len(solutions) == 1
    return sp.factor(solutions[0])


# Generic lower corrections.  The monomial orders agree with HC4RSD39.
A = sp.symbols("A0:5")
t = sp.symbols("t0:10")
b = sp.symbols("b0:6")
ell = sp.symbols("ell0:3")


def binary_quartic(passive_form: sp.Expr) -> sp.Expr:
    return sum(
        A[index] * x ** (4 - index) * passive_form**index
        for index in range(5)
    )


C3 = sum(
    value * monomial
    for value, monomial in zip(
        t,
        (
            z**3,
            y * z**2,
            y**2 * z,
            y**3,
            x * z**2,
            x * y * z,
            x * y**2,
            x**2 * z,
            x**2 * y,
            x**3,
        ),
        strict=True,
    )
)
Q2 = sum(
    value * monomial
    for value, monomial in zip(
        b, (z**2, y * z, y**2, x * z, x * y, x**2), strict=True
    )
)
C1 = ell[0] * z + ell[1] * y + ell[2] * x


def packet_components(
    chart_name: str,
) -> tuple[dict[int, sp.Expr], sp.Symbol]:
    charts = {
        "independent_Ny": (y, z, y),
        "independent_Nz": (y, z, z),
        "independent_Nsum": (y, z, y + z),
        "dependent_Ny": (y, 0, y),
        "dependent_Nz": (y, 0, z),
        "Lzero_Ny": (0, y, y),
        "Lzero_Nz": (0, y, z),
        "bothzero_Ny": (0, 0, y),
    }
    L, M, N = charts[chart_name]
    shift = (
        x * y**3 / 49 + sp.Rational(3, 7) * x**2 * y * z
        if L == y and M == z
        else x * y**3 / 49
        if L == y
        else 0
    )
    transverse = sp.symbols(f"s_{chart_name}")
    complement = z if N == y else y if N == z else z - y
    components = {
        7: x**7,
        6: x**5 * L,
        5: sp.Rational(2, 7) * x**3 * L**2 + x**4 * M,
        4: shift + binary_quartic(N) + transverse * x**3 * complement,
        3: C3,
        2: Q2,
        1: C1,
    }
    return components, transverse


# Three transverse-direction packets first force A4=A3=0.  The remaining
# A2 is nonzero on the curved quartic locus and is normalized to one.
transverse_constants: dict[str, list[str]] = {}
for chart_name, t2_value in (
    ("independent_Nz", sp.Rational(1, 49)),
    ("dependent_Nz", sp.Integer(0)),
):
    components, transverse = packet_components(chart_name)
    degree_fifteen_solution = {
        transverse: sp.Rational(7, 2) * t[6],
        t[2]: t2_value,
        t[3]: 0,
    }
    assert face(components, 21, degree_fifteen_solution) == 0
    degree_fourteen = face(components, 20, degree_fifteen_solution)
    assert coefficient(degree_fourteen, x**8 * z**6) == -16 * A[4] ** 2
    assert coefficient(degree_fourteen, x**9 * z**5) == -24 * A[3] * A[4]
    assert sp.expand(
        coefficient(degree_fourteen, x**10 * z**4).subs(A[4], 0)
        + 12 * A[3] ** 2
    ) == 0

    normalized = {
        **degree_fifteen_solution,
        A[2]: 1,
        A[3]: 0,
        A[4]: 0,
    }
    normalized_fourteen = face(components, 20, normalized)
    assert sp.expand(
        coefficient(normalized_fourteen, x**12 * z**2)
        + 4 * (7 * t[1] - 2) * (7 * t[1] - 1)
    ) == 0

    # First root: an immutable degree-thirteen coefficient.
    root_one = {**normalized, t[1]: sp.Rational(1, 7)}
    b2_root_one = unique_solution(
        coefficient(face(components, 20, root_one), x**14), b[2]
    )
    degree_thirteen_one = face(
        components, 19, {**root_one, b[2]: b2_root_one}
    )
    assert coefficient(degree_thirteen_one, x**10 * y * z**2) == sp.Rational(
        12, 7
    )

    # Second root: the next equation gives t5=2*A1/7, and degree twelve
    # has a parameter-free obstruction.
    root_two = {
        **normalized,
        t[1]: sp.Rational(2, 7),
        t[5]: sp.Rational(2, 7) * A[1],
    }
    b2_root_two = unique_solution(
        coefficient(face(components, 20, root_two), x**14), b[2]
    )
    assert face(
        components,
        19,
        {**root_two, b[2]: b2_root_two},
    ) == 0
    degree_twelve_two = face(
        components,
        18,
        {**root_two, b[2]: b2_root_two},
    )
    assert coefficient(degree_twelve_two, x**7 * y**5) == -sp.Rational(
        12, 2401
    )
    transverse_constants[chart_name] = ["12/7", "-12/2401"]


# The diagonal passive direction has the same A4,A3 collapse.  Its two
# scalar roots are q=0 and q=-7, with the same lower obstructions.
components, transverse = packet_components("independent_Nsum")
sum_degree_fifteen = {
    transverse: -sp.Rational(7, 4) * t[4]
    + sp.Rational(7, 4) * t[5]
    - sp.Rational(7, 4) * t[6]
    + sp.Rational(1, 4),
    t[0]: t[2] - 2 * t[3] - sp.Rational(1, 49),
    t[1]: 2 * t[2] - 3 * t[3] - sp.Rational(2, 49),
}
assert face(components, 21, sum_degree_fifteen) == 0
sum_fourteen = face(components, 20, sum_degree_fifteen)
assert coefficient(sum_fourteen, x**8 * z**6) == -16 * A[4] ** 2
assert sp.expand(
    coefficient(sum_fourteen, x**10 * z**4).subs(A[4], 0)
    + 12 * A[3] ** 2
) == 0
sum_normalized = {
    **sum_degree_fifteen,
    A[2]: 1,
    A[3]: 0,
    A[4]: 0,
}
q_sum = 49 * t[2] - 147 * t[3] + 6
sum_fourteen_normalized = face(components, 20, sum_normalized)
assert sp.expand(
    coefficient(sum_fourteen_normalized, x**12 * z**2)
    + 4 * q_sum * (q_sum + 7) / 49
) == 0

sum_q_zero = {**sum_normalized, t[2]: 3 * t[3] - sp.Rational(6, 49)}
b0_sum_zero = unique_solution(
    coefficient(face(components, 20, sum_q_zero), x**14), b[0]
)
sum_thirteen_zero = face(
    components, 19, {**sum_q_zero, b[0]: b0_sum_zero}
)
assert coefficient(sum_thirteen_zero, x**10 * y**3) == sp.Rational(12, 7)

sum_q_minus_seven = {
    **sum_normalized,
    t[2]: 3 * t[3] - sp.Rational(13, 49),
    t[4]: (-4 * A[1] + 7 * t[5] + 7 * t[6] + 3) / 21,
}
b0_sum_minus_seven = unique_solution(
    coefficient(face(components, 20, sum_q_minus_seven), x**14), b[0]
)
assert face(
    components,
    19,
    {**sum_q_minus_seven, b[0]: b0_sum_minus_seven},
) == 0
sum_twelve_minus_seven = face(
    components,
    18,
    {**sum_q_minus_seven, b[0]: b0_sum_minus_seven},
)
assert coefficient(sum_twelve_minus_seven, x**7 * y**5) == -sp.Rational(
    12, 2401
)
transverse_constants["independent_Nsum"] = ["12/7", "-12/2401"]


# The remaining five charts share one quartic-polar split.  For
# independent_Ny put
#
#   r=s-7*b0, q=2*s-7*t5, p=49*t2-1.
#
# Its degree-fourteen equations are q^2=-28*A2*r,
# q*p=147*A3*r, p^2=-2058*A4*r.  Thus r=0 is the zero
# stratum, while r!=0 forces 3*A3^2=8*A2*A4.
components, transverse = packet_components("independent_Ny")
independent_y_fifteen = {t[0]: 0, t[1]: 0, t[4]: sp.Rational(1, 7)}
assert face(components, 21, independent_y_fifteen) == 0
independent_y_fourteen = face(components, 20, independent_y_fifteen)
r_expr = transverse - 7 * b[0]
q_expr = 2 * transverse - 7 * t[5]
p_expr = 49 * t[2] - 1
assert sp.expand(
    coefficient(independent_y_fourteen, x**14)
    + 28 * A[2] * r_expr
    + q_expr**2
) == 0
assert sp.expand(
    coefficient(independent_y_fourteen, x**13 * y)
    + sp.Rational(4, 7) * (147 * A[3] * r_expr - q_expr * p_expr)
) == 0
assert sp.expand(
    coefficient(independent_y_fourteen, x**12 * y**2)
    + sp.Rational(4, 49) * (2058 * A[4] * r_expr + p_expr**2)
) == 0

# Zero stratum: degree twelve successively kills A4,A3,A2, contradicting
# passive curvature.
independent_y_zero = {
    **independent_y_fifteen,
    transverse: 7 * b[0],
    t[5]: 2 * b[0],
    t[2]: sp.Rational(1, 49),
}
independent_y_zero_twelve = face(components, 18, independent_y_zero)
assert coefficient(independent_y_zero_twelve, x**6 * y**6) == -64 * A[4] ** 2
assert coefficient(
    independent_y_zero_twelve.subs(A[4], 0), x**9 * y**2 * z
) == sp.Rational(36, 7) * A[3]
assert coefficient(
    independent_y_zero_twelve.subs({A[4]: 0, A[3]: 0}),
    x**10 * y * z,
) == sp.Rational(12, 7) * A[2]

# Resonance: parameterize the square Hessian by nonzero r,q.  Degree
# thirteen kills p, and after its two remaining linear equations the next
# face is -3*q^2/(49*r).
r_parameter, q_parameter, p_parameter = sp.symbols(
    "r_parameter q_parameter p_parameter", nonzero=True
)
independent_y_resonance_general = {
    **independent_y_fifteen,
    transverse: r_parameter + 7 * b[0],
    t[5]: (2 * (r_parameter + 7 * b[0]) - q_parameter) / 7,
    t[2]: (p_parameter + 1) / 49,
    A[2]: -q_parameter**2 / (28 * r_parameter),
    A[3]: q_parameter * p_parameter / (147 * r_parameter),
    A[4]: -p_parameter**2 / (2058 * r_parameter),
}
independent_y_thirteen_general = face(
    components, 19, independent_y_resonance_general
)
assert coefficient(independent_y_thirteen_general, x**12 * z) == (
    -4 * r_parameter * p_parameter / 7
)
independent_y_resonance = {
    **independent_y_resonance_general,
    p_parameter: 0,
    t[2]: sp.Rational(1, 49),
    A[3]: 0,
    A[4]: 0,
}
independent_y_thirteen = face(components, 19, independent_y_resonance)
t3_resonance = unique_solution(
    coefficient(independent_y_thirteen, x**12 * y), t[3]
)
t6_resonance = unique_solution(
    coefficient(independent_y_thirteen, x**13), t[6]
)
independent_y_twelve = face(
    components,
    18,
    {
        **independent_y_resonance,
        t[3]: t3_resonance,
        t[6]: t6_resonance,
    },
)
assert coefficient(independent_y_twelve, x**10 * y * z) == (
    -3 * q_parameter**2 / (49 * r_parameter)
)


# The other four degree-fourteen systems are the same rank-one condition.
# On b0!=0 the parameter v below is killed at degree thirteen, leaving the
# pure-x square root u.  Each chart then has one immutable lower coefficient.
resonant_constants: dict[str, str] = {
    "independent_Ny": "-3*q^2/(49*r) at degree 12",
}
u, v = sp.symbols("u v", nonzero=True)


# dependent_Ny: b0=0 is affine-transverse; b0!=0 dies in degree twelve.
components, transverse = packet_components("dependent_Ny")
dependent_y_fifteen = {t[0]: 0, t[1]: 0, t[4]: 0}
assert face(components, 21, dependent_y_fifteen) == 0
dependent_y_fourteen = face(components, 20, dependent_y_fifteen)
q_dependent = 2 * transverse - 7 * t[5]
assert sp.expand(
    coefficient(dependent_y_fourteen, x**14)
    - (196 * A[2] * b[0] - q_dependent**2)
) == 0
assert sp.expand(
    coefficient(dependent_y_fourteen, x**13 * y)
    - 28 * (21 * A[3] * b[0] + q_dependent * t[2])
) == 0
assert sp.expand(
    coefficient(dependent_y_fourteen, x**12 * y**2)
    - 196 * (6 * A[4] * b[0] - t[2] ** 2)
) == 0
dependent_y_resonance = {
    **dependent_y_fifteen,
    t[2]: v,
    t[5]: (2 * transverse - u) / 7,
    A[2]: u**2 / (196 * b[0]),
    A[3]: -u * v / (21 * b[0]),
    A[4]: v**2 / (6 * b[0]),
}
dependent_y_thirteen = face(components, 19, dependent_y_resonance)
assert coefficient(dependent_y_thirteen, x**12 * z) == 196 * b[0] * v
dependent_y_pure = {
    **dependent_y_resonance,
    v: 0,
    t[2]: 0,
    A[3]: 0,
    A[4]: 0,
}
dependent_t3 = unique_solution(
    coefficient(face(components, 19, dependent_y_pure), x**12 * y), t[3]
)
dependent_t6 = unique_solution(
    coefficient(face(components, 19, dependent_y_pure), x**13), t[6]
)
dependent_y_twelve = face(
    components,
    18,
    {**dependent_y_pure, t[3]: dependent_t3, t[6]: dependent_t6},
)
assert coefficient(dependent_y_twelve, x**11 * z) == 4 * b[0] * u
resonant_constants["dependent_Ny"] = "4*b0*u at degree 12"


# Lzero_Ny: its resonance reaches degree ten and dies by -6*u^2.
components, transverse = packet_components("Lzero_Ny")
lzero_y_fifteen = {t[0]: 0, t[1]: 0, t[4]: 0}
assert face(components, 21, lzero_y_fifteen) == 0
lzero_y_fourteen = face(components, 20, lzero_y_fifteen)
assert coefficient(lzero_y_fourteen, x**14) == 49 * (
    4 * A[2] * b[0] - t[5] ** 2
)
assert coefficient(lzero_y_fourteen, x**13 * y) == 196 * (
    3 * A[3] * b[0] - t[2] * t[5]
)
assert coefficient(lzero_y_fourteen, x**12 * y**2) == 196 * (
    6 * A[4] * b[0] - t[2] ** 2
)
lzero_y_resonance = {
    **lzero_y_fifteen,
    t[5]: u,
    t[2]: v,
    A[2]: u**2 / (4 * b[0]),
    A[3]: u * v / (3 * b[0]),
    A[4]: v**2 / (6 * b[0]),
}
lzero_y_thirteen = face(components, 19, lzero_y_resonance)
assert coefficient(lzero_y_thirteen, x**12 * z) == 196 * b[0] * v
lzero_y_pure = {
    **lzero_y_resonance,
    v: 0,
    t[2]: 0,
    A[3]: 0,
    A[4]: 0,
    t[3]: 0,
    t[6]: (2 * b[0] + 7 * b[1] * u - transverse * u) / (14 * b[0]),
}
lzero_y_twelve = face(components, 18, lzero_y_pure)
transverse_lzero = 2 * b[0] / u
assert sp.factor(
    coefficient(lzero_y_twelve, x**11 * y)
    - 7 * u**2 * (-2 * b[0] + transverse * u) / b[0]
) == 0
b2_lzero = unique_solution(
    coefficient(lzero_y_twelve.subs(transverse, transverse_lzero), x**12),
    b[2],
)
lzero_y_after_twelve = {
    **lzero_y_pure,
    transverse: transverse_lzero,
    t[6]: sp.factor(lzero_y_pure[t[6]].subs(transverse, transverse_lzero)),
    b[2]: b2_lzero,
}
lzero_y_eleven = face(components, 17, lzero_y_after_twelve)
t7_lzero = unique_solution(
    coefficient(lzero_y_eleven, x**10 * y), t[7]
)
t7_lzero = sp.factor(t7_lzero.subs(transverse, transverse_lzero))
b3_lzero = unique_solution(
    coefficient(lzero_y_eleven.subs(t[7], t7_lzero), x**11), b[3]
)
b3_lzero = sp.factor(
    b3_lzero.subs({transverse: transverse_lzero, t[7]: t7_lzero})
)
lzero_y_ten = face(
    components,
    16,
    {**lzero_y_after_twelve, t[7]: t7_lzero, b[3]: b3_lzero},
)
assert coefficient(lzero_y_ten, x**8 * y**2) == -6 * u**2
resonant_constants["Lzero_Ny"] = "-6*u^2 at degree 10"


# bothzero_Ny: the last resonance is the nonlinear-coordinate square and
# survives until the immutable degree-eight coefficient -u^6/(4*b0^2).
components, transverse = packet_components("bothzero_Ny")
both_y_fifteen = {t[0]: 0, t[1]: 0, t[4]: 0}
assert face(components, 21, both_y_fifteen) == 0
both_y_fourteen = face(components, 20, both_y_fifteen)
assert coefficient(both_y_fourteen, x**14) == 49 * (
    4 * A[2] * b[0] - t[5] ** 2
)
assert coefficient(both_y_fourteen, x**13 * y) == 196 * (
    3 * A[3] * b[0] - t[2] * t[5]
)
assert coefficient(both_y_fourteen, x**12 * y**2) == 196 * (
    6 * A[4] * b[0] - t[2] ** 2
)
both_y_resonance = {
    **both_y_fifteen,
    t[5]: u,
    t[2]: v,
    A[2]: u**2 / (4 * b[0]),
    A[3]: u * v / (3 * b[0]),
    A[4]: v**2 / (6 * b[0]),
}
both_y_thirteen = face(components, 19, both_y_resonance)
assert coefficient(both_y_thirteen, x**12 * z) == 196 * b[0] * v
both_y_pure = {
    **both_y_resonance,
    v: 0,
    t[2]: 0,
    A[3]: 0,
    A[4]: 0,
    t[3]: 0,
    t[6]: b[1] * u / (2 * b[0]),
}
both_y_twelve = face(components, 18, both_y_pure)
assert sp.factor(
    coefficient(both_y_twelve, x**11 * y)
    - 7 * transverse * u**3 / b[0]
) == 0
both_y_solution = {
    **both_y_pure,
    transverse: 0,
    b[2]: b[1] ** 2 / (4 * b[0]),
    t[7]: 2 * A[1] * b[0] / u,
    b[3]: (-2 * A[1] * b[0] * b[1] + 2 * b[0] * t[8] * u)
    / u**2,
    b[4]: (
        -2 * A[1] * b[0] * b[1] ** 2
        + 2 * b[0] * b[1] * t[8] * u
        + ell[0] * u**3
    )
    / (2 * b[0] * u**2),
}
both_y_eight = face(components, 14, both_y_solution)
assert coefficient(both_y_eight, x**4 * y**4) == -u**6 / (4 * b[0] ** 2)
resonant_constants["bothzero_Ny"] = "-u^6/(4*b0^2) at degree 8"


# Lzero_Nz: the zero stratum dies by A4,A3,A2^2.  The resonance first
# kills its N^4 direction and then has -147*t5^4/(4*r^2).
components, transverse = packet_components("Lzero_Nz")
lzero_z_fifteen = {
    t[2]: 0,
    t[3]: 0,
    t[6]: sp.Rational(1, 7),
}
assert face(components, 21, lzero_z_fifteen) == 0
lzero_z_fourteen = face(components, 20, lzero_z_fifteen)
r_lzero_z_expr = 7 * b[2] - transverse
assert sp.expand(
    coefficient(lzero_z_fourteen, x**14)
    - 7 * (4 * A[2] * r_lzero_z_expr - 7 * t[5] ** 2)
) == 0
assert sp.expand(
    coefficient(lzero_z_fourteen, x**13 * z)
    - 28 * (3 * A[3] * r_lzero_z_expr - 7 * t[1] * t[5])
) == 0
assert sp.expand(
    coefficient(lzero_z_fourteen, x**12 * z**2)
    - 28 * (6 * A[4] * r_lzero_z_expr - 7 * t[1] ** 2)
) == 0
lzero_z_zero = {
    **lzero_z_fifteen,
    b[2]: transverse / 7,
    t[5]: 0,
    t[1]: 0,
}
lzero_z_zero_twelve = face(components, 18, lzero_z_zero)
assert coefficient(lzero_z_zero_twelve, x**6 * z**6) == -64 * A[4] ** 2
assert coefficient(
    lzero_z_zero_twelve.subs(A[4], 0), x**8 * z**4
) == -33 * A[3] ** 2
assert coefficient(
    lzero_z_zero_twelve.subs({A[4]: 0, A[3]: 0}), x**10 * z**2
) == -12 * A[2] ** 2

r_lzero_z = sp.symbols("r_lzero_z", nonzero=True)
lzero_z_resonance_general = {
    **lzero_z_fifteen,
    b[2]: (transverse + r_lzero_z) / 7,
    A[2]: 7 * t[5] ** 2 / (4 * r_lzero_z),
    A[3]: 7 * t[1] * t[5] / (3 * r_lzero_z),
    A[4]: 7 * t[1] ** 2 / (6 * r_lzero_z),
}
lzero_z_thirteen_general = face(components, 19, lzero_z_resonance_general)
assert coefficient(lzero_z_thirteen_general, x**12 * y) == (
    28 * r_lzero_z * t[1]
)
lzero_z_resonance = {
    **lzero_z_resonance_general,
    t[1]: 0,
    A[3]: 0,
    A[4]: 0,
    t[0]: -7 * t[5] ** 3 / (12 * r_lzero_z**2),
    t[4]: t[5] * (7 * b[1] - A[1]) / (2 * r_lzero_z),
}
lzero_z_twelve = face(components, 18, lzero_z_resonance)
assert coefficient(lzero_z_twelve, x**10 * z**2) == (
    -147 * t[5] ** 4 / (4 * r_lzero_z**2)
)
resonant_constants["Lzero_Nz"] = "-147*t5^4/(4*r^2) at degree 12"


# Exact affine-transverse cylinder lemma used on the zero strata of
# dependent_Ny, Lzero_Ny, and bothzero_Ny.  Put
#
#   c=h(x,y)+z*(a(x)*y+d(x)),  deg(a)<=1.
#
# The z^2 coefficient below first makes a constant.  If a!=0 it then makes
# d affine, and after an affine normalization D=y the z coefficient is
# -2*y*h_xx; the remaining square makes h depend on x only through x*y,
# which combines with y*z into one linear form.  If a=0,
# the z coefficient is d*h_yy*(d*d''-2*d'^2), so polynomial degree makes d
# constant.  The binary straight-level lemma then gives a fixed cylinder.
hx, hy, hxx, hxy, hyy = sp.symbols("hx hy hxx hxy hyy")
a0, a1, d0, d1, d2, d3 = sp.symbols("a0 a1 d0 d1 d2 d3")
a_linear = a0 + a1 * x
d_cubic = d0 + d1 * x + d2 * x**2 + d3 * x**3
D_affine_y = a_linear * y + d_cubic
Dx = sp.diff(D_affine_y, x)
Dy = sp.diff(D_affine_y, y)
Dxx = sp.diff(D_affine_y, x, 2)
Dxy = sp.diff(D_affine_y, x, y)
gradient_jet = sp.Matrix([hx + z * Dx, hy + z * Dy, D_affine_y])
hessian_jet = sp.Matrix(
    [
        [hxx + z * Dxx, hxy + z * Dxy, Dx],
        [hxy + z * Dxy, hyy, Dy],
        [Dx, Dy, 0],
    ]
)
affine_transverse_J = sp.Poly(
    sp.expand((gradient_jet.T * hessian_jet.adjugate() * gradient_jet)[0]),
    z,
)
expected_z2 = -D_affine_y * (
    -3 * y * a_linear * a1**2
    + 2 * a_linear**2 * sp.diff(d_cubic, x, 2)
    - 4 * a_linear * a1 * sp.diff(d_cubic, x)
    + d_cubic * a1**2
)
assert sp.expand(affine_transverse_J.coeff_monomial(z**2) - expected_z2) == 0
assert sp.expand(
    affine_transverse_J.coeff_monomial(z)
    .subs({a0: 0, a1: 0})
    - d_cubic
    * hyy
    * (d_cubic * sp.diff(d_cubic, x, 2) - 2 * sp.diff(d_cubic, x) ** 2)
) == 0

# The final affine normalization D=y is checked at the jet level.
normalized_gradient = sp.Matrix([hx, hy + z, y])
normalized_hessian = sp.Matrix(
    [[hxx, hxy, 0], [hxy, hyy, 1], [0, 1, 0]]
)
normalized_J = sp.Poly(
    sp.expand((normalized_gradient.T * normalized_hessian.adjugate() * normalized_gradient)[0]),
    z,
)
assert normalized_J.coeff_monomial(z) == -2 * y * hxx
q_jet, q_prime_jet = sp.symbols("q_jet q_prime_jet")
assert sp.expand(
    normalized_J.coeff_monomial(1).subs(
        {hxx: 0, hx: q_jet, hxy: q_prime_jet}
    )
    + (y * q_prime_jet - q_jet) ** 2
) == 0
assert sp.expand(
    affine_transverse_J.coeff_monomial(1).subs(
        {a0: 0, a1: 0, d1: 0, d2: 0, d3: 0}
    )
    - d0**2 * (hxx * hyy - hxy**2)
) == 0


payload = {
    "format": "hc4-pure-septic-quartic-packets-v1",
    "status": {
        "id": "HC4RSD40",
        "kind": "exact closure theorem",
        "scope": "the eight degree-fifteen packets left by HC4RSD39",
    },
    "common_quartic_polar_split": {
        "equation": "rank-one 2-by-2 system",
        "resonance": "3*A3^2=8*A2*A4, equivalently Q_NN is a square",
    },
    "transverse_direction_packets": transverse_constants,
    "resonant_obstructions": resonant_constants,
    "zero_strata": {
        "independent_Ny": "curvature contradiction at degree 12",
        "Lzero_Nz": "curvature contradiction at degree 12",
        "dependent_Ny": "fixed cylinder by the affine-transverse lemma",
        "Lzero_Ny": "fixed cylinder by the affine-transverse lemma",
        "bothzero_Ny": "fixed cylinder by the affine-transverse lemma",
    },
    "conclusion": (
        "all eight quartic packets are impossible or fixed cylinders; "
        "the scalar pure-septic branch is closed"
    ),
}
serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: exposed the common quartic-polar rank-one split")
print("PASS: closed all three transverse-direction packets")
print("PASS: closed all five aligned quartic packets")
print("PASS: verified the affine-transverse fixed-cylinder identity")
print("THEOREM: every scalar pure-septic direction is fixed")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")
