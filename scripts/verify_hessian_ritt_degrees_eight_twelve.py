#!/usr/bin/env python3
"""Exact Hessian-Ritt incidence and intersection atlas in degrees 8 and 12.

The coefficient calculations use the monic twice-integrated Hessian chart.
Singular is used only for exact minimal-prime counts and dimensions of the
degree-twelve pair and higher intersections.
"""

from __future__ import annotations

import re
import shutil
import subprocess

import sympy as sp


w, z = sp.symbols("w z")


def canonical_atlas(a: int, b: int, prefix: str = "") -> dict:
    """Return the canonical monic A o B chart and its residual equations."""

    degree = a * b
    coefficients = {j: sp.symbols(f"c{j}") for j in range(2, degree)}
    inner = {
        j: sp.symbols(f"{prefix}b{j}") for j in range(1, b)
    }
    outer = {
        m: sp.symbols(f"{prefix}a{m}") for m in range(1, a)
    }
    B = w**b + sum(inner[j] * w**j for j in inner)
    A = z**a + sum(outer[m] * z**m for m in outer)
    composition = sp.Poly(sp.expand(A.subs(z, B)), w)

    reconstruction: dict[sp.Symbol, sp.Expr] = {}
    used_degrees: list[int] = []
    for j in range(b - 1, 0, -1):
        coefficient_degree = degree - b + j
        equation = sp.expand(
            composition.nth(coefficient_degree).subs(reconstruction)
            - coefficients[coefficient_degree]
        )
        reconstruction[inner[j]] = sp.factor(sp.solve(equation, inner[j])[0])
        used_degrees.append(coefficient_degree)

    for m in range(a - 1, 0, -1):
        coefficient_degree = m * b
        equation = sp.expand(
            composition.nth(coefficient_degree).subs(reconstruction)
            - coefficients[coefficient_degree]
        )
        reconstruction[outer[m]] = sp.factor(sp.solve(equation, outer[m])[0])
        used_degrees.append(coefficient_degree)

    residuals: list[tuple[int, sp.Expr]] = []
    for coefficient_degree in range(2, degree):
        if coefficient_degree in used_degrees:
            continue
        numerator = sp.factor(
            sp.together(
                composition.nth(coefficient_degree).subs(reconstruction)
                - coefficients[coefficient_degree]
            ).as_numer_denom()[0]
        )
        residuals.append((coefficient_degree, numerator))

    parameters = list(inner.values()) + list(outer.values())
    coefficient_image = {
        coefficients[j]: composition.nth(j) for j in range(2, degree)
    }
    for _, residual in residuals:
        assert sp.expand(residual.subs(coefficient_image)) == 0
    assert len(reconstruction) == a + b - 2
    assert len(residuals) == degree - a - b

    return {
        "a": a,
        "b": b,
        "degree": degree,
        "parameters": parameters,
        "B": B,
        "A": A,
        "composition": composition,
        "coefficients": coefficients,
        "reconstruction": reconstruction,
        "residuals": residuals,
    }


def pulled_equations(source: dict, target: dict) -> list[sp.Expr]:
    substitution = {
        target["coefficients"][j]: source["composition"].nth(j)
        for j in range(2, source["degree"])
    }
    equations = []
    for _, residual in target["residuals"]:
        pulled = sp.factor(residual.subs(substitution))
        primitive = sp.primitive(sp.Poly(pulled, *source["parameters"]))[1]
        equations.append(primitive.as_expr())
    return equations


def singular_summary(
    variables: list[sp.Symbol], equations: list[sp.Expr]
) -> tuple[int, int]:
    singular = shutil.which("Singular")
    assert singular is not None, "degree-twelve minimal-prime audit requires Singular"
    serialized = [str(equation).replace("**", "^") for equation in equations]
    program = (
        'LIB "primdec.lib";\n'
        f'ring r=0,({",".join(map(str, variables))}),dp;\n'
        "option(redSB);\n"
        f'ideal I={",".join(serialized)};\n'
        'print("RESULT");\n'
        "dim(std(I));\n"
        "list L=minAssGTZ(I);\n"
        "size(L);\n"
    )
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=120,
    )
    match = re.search(r"RESULT\s+(\d+)\s+(\d+)", result.stdout)
    assert match is not None, result.stdout + result.stderr
    return int(match.group(1)), int(match.group(2))


def normalized_clean_seed(polynomial: sp.Expr) -> sp.Expr | None:
    polynomial = sp.expand(polynomial)
    tangent = sp.diff(polynomial, w).subs(w, 0)
    endpoint = sp.factor(polynomial.subs(w, 1) - polynomial.subs(w, 0) - tangent)
    if endpoint != 0:
        return None
    derivative_gap = sp.factor(sp.diff(polynomial, w).subs(w, 1) - tangent)
    if derivative_gap == 0:
        return None
    H = sp.factor(
        -(polynomial - polynomial.subs(w, 0) - tangent * w) / derivative_gap
    )
    K = sp.diff(H, w, 2)
    Q = sp.cancel(H / w**2)
    clean_quantities = (
        sp.discriminant(K, w),
        K.subs(w, 0),
        sp.discriminant(Q, w),
        K.subs(w, 1) + 2,
    )
    if any(quantity == 0 for quantity in clean_quantities):
        return None
    assert H.subs(w, 0) == 0
    assert sp.diff(H, w).subs(w, 0) == 0
    assert H.subs(w, 1) == 0
    assert sp.diff(H, w).subs(w, 1) == -1
    return H


# ---------------------------------------------------------------------------
# Degree eight.

D24 = canonical_atlas(2, 4, "u")
D42 = canonical_atlas(4, 2, "v")
assert [degree for degree, _ in D24["residuals"]] == [2, 3]
assert [degree for degree, _ in D42["residuals"]] == [3, 5]

c2, c3, c4, c5, c6, c7 = sp.symbols("c2 c3 c4 c5 c6 c7")
expected_24 = {
    2: (
        -512 * c2
        + 256 * c4 * c6
        - 64 * c4 * c7**2
        + 128 * c5**2
        - 256 * c5 * c6 * c7
        + 64 * c5 * c7**3
        - 64 * c6**3
        + 144 * c6**2 * c7**2
        - 60 * c6 * c7**4
        + 7 * c7**6
    ),
    3: (
        -128 * c3
        + 64 * c4 * c7
        + 64 * c5 * c6
        - 48 * c5 * c7**2
        - 48 * c6**2 * c7
        + 40 * c6 * c7**3
        - 7 * c7**5
    ),
}
expected_42 = {
    3: -2 * (256 * c3 - 128 * c4 * c7 + 20 * c6 * c7**3 - 7 * c7**5),
    5: -32 * c5 + 24 * c6 * c7 - 7 * c7**3,
}
for degree, residual in D24["residuals"]:
    assert sp.expand(residual - expected_24[degree]) == 0
for degree, residual in D42["residuals"]:
    assert sp.expand(residual - expected_42[degree]) == 0

u1, u2, u3, ua1 = D24["parameters"]
degree_eight_intersection = 8 * u1 - 4 * u2 * u3 + u3**3
pulled_42 = pulled_equations(D24, D42)
assert sp.factor(pulled_42[1] + degree_eight_intersection) == 0
assert sp.factor(
    pulled_42[0] - (u3**2 - u2) * degree_eight_intersection
) == 0

center = -u3 / 4
completed_quartic = sp.expand(D24["B"] + ua1 / 2)
centered_quartic = sp.Poly(sp.expand(completed_quartic.subs(w, w + center)), w)
assert centered_quartic.nth(3) == 0
assert sp.factor(8 * centered_quartic.nth(1) - degree_eight_intersection) == 0

# A rational clean point on the 2 o 2 o 2 component.
inner_8 = w**2 - 3 * w
middle_8 = inner_8**2 - 3 * inner_8
f8 = sp.expand(middle_8**2 - 100 * middle_8)
H8 = normalized_clean_seed(f8)
assert H8 == -w**2 * (w - 1) * (
    w**5 - 11 * w**4 + 37 * w**3 - 17 * w**2 - 189 * w + 519
) / 340


# ---------------------------------------------------------------------------
# Degree twelve incidence loci and pairwise intersections.

pairs = ((2, 6), (3, 4), (4, 3), (6, 2))
atlases = {pair: canonical_atlas(*pair, f"p{pair[0]}{pair[1]}_") for pair in pairs}
for a, b in pairs:
    atlas = atlases[(a, b)]
    assert len(atlas["parameters"]) == a + b - 2
    assert len(atlas["residuals"]) == 12 - a - b

pair_expectations = {
    ((2, 6), (6, 2)): (4, 1),
    ((2, 6), (4, 3)): (4, 1),
    ((2, 6), (3, 4)): (3, 1),
    ((3, 4), (6, 2)): (4, 1),
    ((4, 3), (6, 2)): (3, 1),
    ((3, 4), (4, 3)): (2, 2),
}
pair_equations: dict[tuple[tuple[int, int], tuple[int, int]], list[sp.Expr]] = {}
for ordered_pair, expected in pair_expectations.items():
    source_pair, target_pair = ordered_pair
    equations = pulled_equations(atlases[source_pair], atlases[target_pair])
    pair_equations[ordered_pair] = equations
    assert singular_summary(atlases[source_pair]["parameters"], equations) == expected

# The easy common-refinement charts collapse to compact equations.
source_34 = atlases[(3, 4)]
b1, b2, b3, a1, a2 = source_34["parameters"]
equations_34_62 = pair_equations[((3, 4), (6, 2))]
common_322 = b3**3 - 4 * b2 * b3 + 8 * b1
assert any(sp.factor(equation + common_322) == 0 for equation in equations_34_62)
for equation in equations_34_62:
    assert sp.rem(equation, common_322, b1) == 0

source_43 = atlases[(4, 3)]
q1, q2, d1, d2, d3 = source_43["parameters"]
equations_43_62 = pair_equations[((4, 3), (6, 2))]
common_left_ritt_1 = d3**3 - 4 * d2 * d3 + 8 * d1
common_left_ritt_2 = 8 * q2**3 - 36 * q1 * q2 + 27 * d3
G_left = sp.groebner(
    [common_left_ritt_1, common_left_ritt_2],
    *source_43["parameters"],
    order="grevlex",
)
for equation in equations_43_62:
    assert G_left.reduce(equation)[1] == 0

# The coprime 3 o 4 versus 4 o 3 intersection has exactly two components.
power_equations = [
    a2**2 - 3 * a1,
    3 * b3**2 - 8 * b2,
    b2**2 - 3 * b1 * b3 + 4 * a2,
]

shift, power_parameter = sp.symbols("shift power_parameter")
power_quartic = (w - shift) ** 4 + power_parameter * (w - shift)
power_constant = power_quartic.subs(w, 0)
power_B = sp.expand(power_quartic - power_constant)
power_A = sp.expand((z + power_constant) ** 3 - power_constant**3)
power_parameter_substitution = {
    b1: sp.Poly(power_B, w).nth(1),
    b2: sp.Poly(power_B, w).nth(2),
    b3: sp.Poly(power_B, w).nth(3),
    a1: sp.Poly(power_A, z).nth(1),
    a2: sp.Poly(power_A, z).nth(2),
}
for equation in power_equations:
    assert sp.expand(equation.subs(power_parameter_substitution)) == 0
for equation in pair_equations[((3, 4), (4, 3))]:
    assert sp.expand(equation.subs(power_parameter_substitution)) == 0

scale, translation = sp.symbols("scale translation", nonzero=True)
T3 = sp.chebyshevt(3, z)
T4 = sp.chebyshevt(4, w)
cheb_value = sp.chebyshevt(4, translation)
cheb_B = sp.expand(
    (sp.chebyshevt(4, scale * w + translation) - cheb_value)
    / (8 * scale**4)
)
cheb_A = sp.expand(
    (
        sp.chebyshevt(3, 8 * scale**4 * z + cheb_value)
        - sp.chebyshevt(3, cheb_value)
    )
    / (4 * (8 * scale**4) ** 3)
)
cheb_substitution = {
    b1: sp.Poly(cheb_B, w).nth(1),
    b2: sp.Poly(cheb_B, w).nth(2),
    b3: sp.Poly(cheb_B, w).nth(3),
    a1: sp.Poly(cheb_A, z).nth(1),
    a2: sp.Poly(cheb_A, z).nth(2),
}
for equation in pair_equations[((3, 4), (4, 3))]:
    assert sp.factor(equation.subs(cheb_substitution)) == 0

# Chebyshev lies on all four loci.  The power component meets either even
# factor-order locus only at its monomial specialization.
for target_pair in ((2, 6), (6, 2)):
    equations = pulled_equations(source_34, atlases[target_pair])
    for equation in equations:
        assert sp.factor(equation.subs(cheb_substitution)) == 0
    power_pullback = [
        sp.factor(equation.subs(power_parameter_substitution))
        for equation in equations
    ]
    assert all(expression.subs(power_parameter, 0) == 0 for expression in power_pullback)
    assert any(
        expression != 0
        and sp.factor(expression / power_parameter).is_polynomial(shift, power_parameter)
        for expression in power_pullback
    )

# Higher intersections, computed in the 3 o 4 chart.
I26 = pulled_equations(source_34, atlases[(2, 6)])
I43 = pulled_equations(source_34, atlases[(4, 3)])
I62 = pulled_equations(source_34, atlases[(6, 2)])
higher_expectations = [
    (I26 + I62, (3, 1)),       # 2o6 cap 3o4 is already inside 6o2
    (I26 + I43, (2, 1)),       # Chebyshev
    (I43 + I62, (2, 1)),       # Chebyshev
    (I26 + I43 + I62, (2, 1)), # all four
]
for equations, expected in higher_expectations:
    assert singular_summary(source_34["parameters"], equations) == expected


# ---------------------------------------------------------------------------
# Clean normalized witnesses for every degree-twelve component type.

def assert_clean(polynomial: sp.Expr) -> None:
    assert normalized_clean_seed(polynomial) is not None


# Common refinements 2o3o2, 2o2o3, and 3o2o2.
q_inner = w**2 - 3 * w
c_middle = q_inner**3 - 3 * q_inner**2 - 3 * q_inner
assert_clean(c_middle**2 + sp.Rational(196, 23) * c_middle)

c_inner = w**3 - 3 * w**2 - 3 * w
q_middle = c_inner**2 - 3 * c_inner
assert_clean(q_middle**2 - sp.Rational(1600, 31) * q_middle)

q_inner = w**2 - 3 * w
q_middle = q_inner**2 - 3 * q_inner
assert_clean(q_middle**3 - 3 * q_middle**2 - 700 * q_middle)

# Degree-six Ritt collision composed on the right and on the left.
right_quadratic = w**2 - 6 * w
right_ritt = (
    (right_quadratic + 5) ** 3
    - sp.Rational(775, 7) * (right_quadratic + 5)
) ** 2
assert_clean(right_ritt)

degree_six_ritt = ((w + 6) ** 3 - 6 * (w + 6)) ** 2
left_ritt = degree_six_ritt**2 - sp.Rational(4779325201, 21481) * degree_six_ritt
assert_clean(left_ritt)

# Power component.
power_witness = ((w - 1) ** 4 + sp.Rational(11, 2) * (w - 1)) ** 3
assert_clean(power_witness)

# Chebyshev component: translation=1 gives ten squarefree algebraic endpoint
# scales, all disjoint from every clean-boundary obstruction.
cheb_scale = sp.symbols("cheb_scale")
cheb_polynomial = sp.chebyshevt(12, cheb_scale * w + 1)
cheb_tangent = sp.diff(cheb_polynomial, w).subs(w, 0)
cheb_endpoint = sp.factor(
    cheb_polynomial.subs(w, 1)
    - cheb_polynomial.subs(w, 0)
    - cheb_tangent
)
cheb_endpoint_factor = sp.factor(cheb_endpoint / (8 * cheb_scale**2))
assert sp.degree(cheb_endpoint_factor, cheb_scale) == 10
assert sp.gcd(cheb_endpoint_factor, sp.diff(cheb_endpoint_factor, cheb_scale)) == 1
cheb_D = sp.factor(sp.diff(cheb_polynomial, w).subs(w, 1) - cheb_tangent)
cheb_K = sp.diff(cheb_polynomial, w, 2)
cheb_Q = sp.cancel(
    (cheb_polynomial - cheb_polynomial.subs(w, 0) - cheb_tangent * w) / w**2
)
cheb_bad_factors = (
    cheb_D,
    cheb_K.subs(w, 0),
    2 * cheb_D - cheb_K.subs(w, 1),
    sp.discriminant(cheb_K, w),
    sp.discriminant(cheb_Q, w),
)
for bad_factor in cheb_bad_factors:
    assert sp.gcd(cheb_endpoint_factor, bad_factor) == 1

print("PASS: the canonical Hessian atlas has exact codimension N-a-b in every factor pair")
print("PASS: degree eight has one excess 2o2o2 intersection of marked dimension two")
print("PASS: all six degree-twelve pair intersections and their dimensions are exact")
print("PASS: the 3o4/4o3 intersection has exactly power and Chebyshev components")
print("PASS: the all-four degree-twelve intersection is exactly the Chebyshev component")
print("PASS: every degree-eight and degree-twelve component meets the marked clean locus")
