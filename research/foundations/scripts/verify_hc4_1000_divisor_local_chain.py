#!/usr/bin/env python3
"""Verify the first divisor-local polynomiality chain in chart 1000.

This works only with the leading Laurent coefficient along L=0.  It avoids
forming the fully cancelled rational expressions for n2, n1, and n0.

Write

    A = f_aa*b + 2*(f_ab*h_a-f_aa*h_b),
    C = 4*f_aa*b^3 + 2*(f_ab*g_a-f_aa*g_b) - f_ab,
    B = f_aa^3*f_bbb - 5*f_aa^3 - 3*f_aa^2*f_ab*f_abb
        + 3*f_aa*f_aab*f_ab^2 - f_aaa*f_ab^3.

The leading Laurent coefficients of the four recursively forced pure-normal
third derivatives, in the order n3,n2,n1,n0, are

    -A^3*B/(4*f_aa^3),
    -A^2*C*B/(4*f_aa^3),
    -A*B*(3*C^2-2*f_aa*kappa)/(12*f_aa^3),
    -C*B*(C^2-2*f_aa*kappa)/(4*f_aa^3),

all multiplying L^-3.  Considered alone, these residues require either
B=0, or A=0 together with C=0 or C^2=2*f_aa*kappa.  Polynomiality of the
already-forced normal Hessian is stronger: its residues require

    A=0,  C^2=2*f_aa*kappa.

Only this unique branch is therefore relevant for a polynomial potential.

On A=0 and C^2=2*f_aa*kappa, the remaining Laurent numerators are four
explicit expressions K1,...,K4.  They admit derivative identities that
show, componentwise, that h_a != 0 forces a double divisibility p^2 | A
on every reduced irreducible caustic component p.
"""

from __future__ import annotations

import runpy

import sympy as sp


data = runpy.run_path("scripts/verify_hc4_x_caustic_formal_compatibility.py")

b = data["b"]
f_aa = data["f_aa"]
f_ab = data["f_ab"]
f_bb = data["f_bb"]
L = data["L"]
kappa = data["kappa"]
n0, n1, n2, n3 = (data[name] for name in ("n0", "n1", "n2", "n3"))
prolongation = data["prolongation"]
top_remainder = data["top_remainder"]

ell = sp.symbols("ell")
divisor_substitution = {f_bb: 5 * b + f_ab**2 / f_aa - ell / f_aa}


def coefficients(expression: sp.Expr, maximum_degree: int) -> list[sp.Expr]:
    """Return the ell coefficients after using ell=L."""

    converted = sp.expand(expression.subs(divisor_substitution))
    polynomial = sp.Poly(converted, ell, domain="EX")
    return [polynomial.nth(degree) for degree in range(maximum_degree + 1)]


def add(*vectors: list[sp.Expr]) -> list[sp.Expr]:
    size = max(map(len, vectors))
    return [
        sum(vector[index] if index < len(vector) else 0 for vector in vectors)
        for index in range(size)
    ]


def scale(vector: list[sp.Expr], scalar: sp.Expr) -> list[sp.Expr]:
    return [scalar * entry for entry in vector]


def multiply(
    left: list[sp.Expr], right: list[sp.Expr], maximum_degree: int
) -> list[sp.Expr]:
    product = [sp.Integer(0)] * (maximum_degree + 1)
    for left_degree, left_coefficient in enumerate(left):
        for right_degree, right_coefficient in enumerate(right):
            degree = left_degree + right_degree
            if degree <= maximum_degree:
                product[degree] += left_coefficient * right_coefficient
    return product


equations = [prolongation.nth(power) for power in range(4)]
pivots = [
    sp.diff(equations[0], n0),
    sp.diff(equations[1], n1),
    sp.diff(equations[2], n2),
    sp.diff(equations[3], n3),
]
remainders = [
    equations[0].subs({n0: 0, n1: 0, n2: 0, n3: 0}),
    equations[1].subs({n1: 0, n2: 0, n3: 0}),
    equations[2].subs({n2: 0, n3: 0}),
]
couplings = {
    (0, 1): sp.diff(equations[0], n1),
    (0, 2): sp.diff(equations[0], n2),
    (0, 3): sp.diff(equations[0], n3),
    (1, 2): sp.diff(equations[1], n2),
    (1, 3): sp.diff(equations[1], n3),
    (2, 3): sp.diff(equations[2], n3),
}

# n3 = N3/D3.
n3_numerator = -top_remainder
n3_denominator = 256 * L**3

maximum_degree = 15
series = {
    "N3": coefficients(n3_numerator, maximum_degree),
    "D3": coefficients(n3_denominator, maximum_degree),
}
for index, pivot in enumerate(pivots):
    series[f"p{index}"] = coefficients(pivot, maximum_degree)
for index, remainder in enumerate(remainders):
    series[f"r{index}"] = coefficients(remainder, maximum_degree)
for indices, coupling in couplings.items():
    series[f"a{indices[0]}{indices[1]}"] = coefficients(
        coupling, maximum_degree
    )

# n2 = N2/(p2*D3).
series["N2"] = scale(
    add(
        multiply(series["a23"], series["N3"], maximum_degree),
        multiply(series["r2"], series["D3"], maximum_degree),
    ),
    -1,
)
series["D2"] = multiply(series["p2"], series["D3"], maximum_degree)

# n1 = N1/(p1*p2*D3).
series["N1"] = scale(
    add(
        multiply(series["a12"], series["N2"], maximum_degree),
        multiply(
            multiply(series["a13"], series["N3"], maximum_degree),
            series["p2"],
            maximum_degree,
        ),
        multiply(
            multiply(series["r1"], series["p2"], maximum_degree),
            series["D3"],
            maximum_degree,
        ),
    ),
    -1,
)
series["D1"] = multiply(series["p1"], series["D2"], maximum_degree)

# n0 = N0/(p0*p1*p2*D3).
series["N0"] = scale(
    add(
        multiply(series["a01"], series["N1"], maximum_degree),
        multiply(
            multiply(series["a02"], series["N2"], maximum_degree),
            series["p1"],
            maximum_degree,
        ),
        multiply(
            multiply(
                series["a03"], series["N3"], maximum_degree
            ),
            multiply(series["p1"], series["p2"], maximum_degree),
            maximum_degree,
        ),
        multiply(series["r0"], series["D1"], maximum_degree),
    ),
    -1,
)
series["D0"] = multiply(series["p0"], series["D1"], maximum_degree)


A = b * f_aa - 2 * f_aa * data["h_b"] + 2 * f_ab * data["h_a"]
B = data["second_residue_factor"]
C = (
    4 * b**3 * f_aa
    - 2 * f_aa * data["g_b"]
    + 2 * f_ab * data["g_a"]
    - f_ab
)
kernel_derivative_L = (
    f_ab * data["total_derivative"](L, data["tangential_a"])
    - f_aa * data["total_derivative"](L, data["tangential_b"])
)
assert (
    sp.factor(
        f_aa * kernel_derivative_L
        - B
        + (f_aa * data["f_aab"] - data["f_aaa"] * f_ab) * L
    )
    == 0
)

expected = {
    "N3": (0, -64 * A**3 * B / f_aa**3),
    "D3": (3, sp.Integer(256)),
    "N2": (4, -117964800 * A**2 * C * B / f_aa**3),
    "D2": (7, sp.Integer(471859200)),
    "N1": (
        8,
        -54358179840000
        * A
        * B
        * (3 * C**2 - 2 * f_aa * kappa)
        / f_aa**3,
    ),
    "D1": (11, sp.Integer(652298158080000)),
    "N0": (
        12,
        -56358560858112000000
        * C
        * B
        * (C**2 - 2 * f_aa * kappa)
        / f_aa**3,
    ),
    "D0": (15, sp.Integer(225434243432448000000)),
}

for name, (order, leading_coefficient) in expected.items():
    assert all(sp.factor(entry) == 0 for entry in series[name][:order])
    assert sp.factor(series[name][order] - leading_coefficient) == 0

H1 = (
    2 * f_aa**3 * data["h_bb"]
    - f_aa**3
    - 4 * f_aa**2 * f_ab * data["h_ab"]
    - 4 * f_aa**2 * data["f_abb"] * data["h_a"]
    + 8 * f_aa * data["f_aab"] * f_ab * data["h_a"]
    + 2 * f_aa * f_ab**2 * data["h_aa"]
    - 4 * data["f_aaa"] * f_ab**2 * data["h_a"]
)
H2 = (
    f_aa**2 * data["h_ab"]
    - f_aa * data["f_aab"] * data["h_a"]
    - f_aa * f_ab * data["h_aa"]
    + data["f_aaa"] * f_ab * data["h_a"]
)
f_bbb_on_B = (
    5
    + 3 * f_ab * data["f_abb"] / f_aa
    - 3 * data["f_aab"] * f_ab**2 / f_aa**2
    + data["f_aaa"] * f_ab**3 / f_aa**3
)
assert (
    sp.factor(
        series["N3"][1].subs(data["f_bbb"], f_bbb_on_B)
        - 96 * A**2 * H1 / f_aa**3
    )
    == 0
)
assert (
    sp.factor(
        series["N3"][2].subs(data["f_bbb"], f_bbb_on_B)
        - 768 * data["h_a"] * A * H2 / f_aa**3
    )
    == 0
)

A_substitution = {data["h_b"]: b / 2 + f_ab * data["h_a"] / f_aa}
C_substitution = {
    data["g_b"]: 2 * b**3
    + f_ab * data["g_a"] / f_aa
    - f_ab / (2 * f_aa)
}
assert all(
    sp.factor(series["N3"][degree].subs(A_substitution)) == 0
    for degree in range(3)
)
assert all(
    sp.factor(series["N2"][degree].subs(A_substitution)) == 0
    for degree in (4, 5)
)
assert (
    sp.factor(
        series["N2"][6].subs(A_substitution)
        - 471859200 * data["h_a"] * C * H2 / f_aa**3
    )
    == 0
)

J1 = (
    2 * f_aa**3 * data["h_bb"]
    - f_aa**3
    - 4 * f_aa**2 * f_ab * data["h_ab"]
    + 4 * f_aa**2 * data["f_abb"] * data["h_a"]
    - 8 * f_aa * data["f_aab"] * f_ab * data["h_a"]
    + 2 * f_aa * f_ab**2 * data["h_aa"]
    + 4 * data["f_aaa"] * f_ab**2 * data["h_a"]
)
J2 = f_aa * data["h_aa"] + 2 * data["f_aaa"] * data["h_a"]
J3 = (
    6 * b**2 * f_aa**3
    - f_aa**3 * data["g_bb"]
    + 2 * f_aa**2 * f_ab * data["g_ab"]
    - 2 * f_aa**2 * data["f_abb"] * data["g_a"]
    + f_aa**2 * data["f_abb"]
    + 4 * f_aa * data["f_aab"] * f_ab * data["g_a"]
    - 2 * f_aa * data["f_aab"] * f_ab
    - f_aa * f_ab**2 * data["g_aa"]
    - 2 * data["f_aaa"] * f_ab**2 * data["g_a"]
    + data["f_aaa"] * f_ab**2
)
J4 = (
    21 * b * f_aa**2
    - 8 * f_aa * data["g_aa"]
    - 16 * data["f_aaa"] * data["g_a"]
    + 8 * data["f_aaa"]
)
AC_substitution = A_substitution | C_substitution
assert (
    sp.factor(
        series["N1"][9].subs(AC_substitution)
        - 54358179840000 * kappa * J1 / f_aa**2
    )
    == 0
)
assert (
    sp.factor(
        series["N1"][10].subs(AC_substitution)
        + 108716359680000 * kappa * J2 / f_aa**2
    )
    == 0
)
assert (
    sp.factor(
        series["N0"][13].subs(AC_substitution)
        + 112717121716224000000 * kappa * J3 / f_aa**2
    )
    == 0
)
assert (
    sp.factor(
        series["N0"][14].subs(AC_substitution)
        - 14089640214528000000 * kappa * J4 / f_aa**2
    )
    == 0
)

K1 = (
    2 * f_aa**3 * data["h_bb"]
    - f_aa**3
    - 4 * f_aa**2 * f_ab * data["h_ab"]
    - 2 * f_aa**2 * data["f_abb"] * data["h_a"]
    + 4 * f_aa * data["f_aab"] * f_ab * data["h_a"]
    + 2 * f_aa * f_ab**2 * data["h_aa"]
    - 2 * data["f_aaa"] * f_ab**2 * data["h_a"]
)
K2 = (
    4 * b**3 * f_aa**2 * data["h_aa"]
    + 8 * b**3 * f_aa * data["f_aaa"] * data["h_a"]
    - 12 * f_aa**2 * data["g_a"] * data["h_ab"]
    - 12 * f_aa**2 * data["g_ab"] * data["h_a"]
    - 2 * f_aa**2 * data["g_b"] * data["h_aa"]
    + 6 * f_aa**2 * data["h_ab"]
    - 4 * f_aa * data["f_aaa"] * data["g_b"] * data["h_a"]
    + 24 * f_aa * data["f_aab"] * data["g_a"] * data["h_a"]
    - 12 * f_aa * data["f_aab"] * data["h_a"]
    + 14 * f_aa * f_ab * data["g_a"] * data["h_aa"]
    + 12 * f_aa * f_ab * data["g_aa"] * data["h_a"]
    - 7 * f_aa * f_ab * data["h_aa"]
    - 20 * data["f_aaa"] * f_ab * data["g_a"] * data["h_a"]
    + 10 * data["f_aaa"] * f_ab * data["h_a"]
)
K3 = (
    4 * b**3 * f_aa**2 * data["f_aab"]
    - 4 * b**3 * f_aa * data["f_aaa"] * f_ab
    + 24 * b**2 * f_aa**3
    - 4 * f_aa**3 * data["g_bb"]
    - 2 * f_aa**2 * data["f_aab"] * data["g_b"]
    + 8 * f_aa**2 * f_ab * data["g_ab"]
    + 4 * f_aa**2 * data["f_abb"] * data["g_a"]
    - 2 * f_aa**2 * data["f_abb"]
    + 2 * f_aa * data["f_aaa"] * f_ab * data["g_b"]
    - 6 * f_aa * data["f_aab"] * f_ab * data["g_a"]
    + 3 * f_aa * data["f_aab"] * f_ab
    - 4 * f_aa * f_ab**2 * data["g_aa"]
    + 2 * data["f_aaa"] * f_ab**2 * data["g_a"]
    - data["f_aaa"] * f_ab**2
)
K4 = (
    84 * b**4 * f_aa**3
    - 32 * b**3 * f_aa**2 * data["g_aa"]
    - 64 * b**3 * f_aa * data["f_aaa"] * data["g_a"]
    + 32 * b**3 * f_aa * data["f_aaa"]
    - 42 * b * f_aa**3 * data["g_b"]
    + 42 * b * f_aa**2 * f_ab * data["g_a"]
    - 21 * b * f_aa**2 * f_ab
    + 96 * f_aa**2 * data["g_a"] * data["g_ab"]
    + 16 * f_aa**2 * data["g_aa"] * data["g_b"]
    - 48 * f_aa**2 * data["g_ab"]
    + 32 * f_aa * data["f_aaa"] * data["g_a"] * data["g_b"]
    - 16 * f_aa * data["f_aaa"] * data["g_b"]
    - 96 * f_aa * data["f_aab"] * data["g_a"] ** 2
    + 96 * f_aa * data["f_aab"] * data["g_a"]
    - 24 * f_aa * data["f_aab"]
    - 112 * f_aa * f_ab * data["g_a"] * data["g_aa"]
    + 56 * f_aa * f_ab * data["g_aa"]
    + 64 * data["f_aaa"] * f_ab * data["g_a"] ** 2
    - 64 * data["f_aaa"] * f_ab * data["g_a"]
    + 16 * data["f_aaa"] * f_ab
)
quadratic_C_substitution = A_substitution | {kappa: C**2 / (2 * f_aa)}
assert (
    sp.factor(
        series["N1"][9].subs(quadratic_C_substitution)
        - 108716359680000 * C**2 * K1 / f_aa**3
    )
    == 0
)
assert (
    sp.factor(
        series["N1"][10].subs(quadratic_C_substitution)
        + 54358179840000 * C * K2 / f_aa**3
    )
    == 0
)
assert (
    sp.factor(
        series["N0"][13].subs(quadratic_C_substitution)
        + 56358560858112000000 * C**2 * K3 / f_aa**3
    )
    == 0
)
assert (
    sp.factor(
        series["N0"][14].subs(quadratic_C_substitution)
        - 7044820107264000000 * C * K4 / f_aa**3
    )
    == 0
)

A_a = data["total_derivative"](A, data["tangential_a"]).subs(A_substitution)
A_b = data["total_derivative"](A, data["tangential_b"]).subs(A_substitution)
C_a = data["total_derivative"](C, data["tangential_a"])
C_b = data["total_derivative"](C, data["tangential_b"])
S_a = 2 * f_aa * C_a - C * data["f_aaa"]
S_b = 2 * f_aa * C_b - C * data["f_aab"]
R = (
    21 * b * f_aa**2
    - 8 * f_aa * data["g_aa"]
    + 8 * data["f_aaa"] * data["g_a"]
    - 4 * data["f_aaa"]
)
assert sp.factor(H2 + f_aa * A_a / 2) == 0
assert sp.factor(K1 + f_aa**2 * A_b + 2 * f_ab * H2) == 0
assert (
    sp.factor(
        K2
        - C * (f_aa * data["h_aa"] - data["f_aaa"] * data["h_a"])
        + 6 * (2 * data["g_a"] - 1) * H2
        - 3 * data["h_a"] * S_a
    )
    == 0
)
assert sp.factor(K3 - f_aa * S_b + f_ab * S_a) == 0
assert (
    sp.factor(K4 - C * R + 12 * (2 * data["g_a"] - 1) * S_a)
    == 0
)

line_slope = sp.symbols("line_slope")
line_h_substitution = {
    data["h_a"]: 0,
    data["h_ab"]: -line_slope * data["h_aa"],
    data["h_bb"]: sp.Rational(1, 2)
    + line_slope**2 * data["h_aa"],
}
line_characteristic = f_aa * line_slope + f_ab
assert (
    sp.factor(
        H2.subs(line_h_substitution)
        + f_aa * line_characteristic * data["h_aa"]
    )
    == 0
)
assert (
    sp.factor(
        K1.subs(line_h_substitution)
        - 2 * f_aa * line_characteristic**2 * data["h_aa"]
    )
    == 0
)
assert (
    sp.factor(
        K2.subs(line_h_substitution)
        - f_aa
        * data["h_aa"]
        * (C + 6 * (2 * data["g_a"] - 1) * line_characteristic)
    )
    == 0
)

E = (
    f_aa**2 * data["f_abb"]
    - 2 * f_aa * f_ab * data["f_aab"]
    + f_ab**2 * data["f_aaa"]
)
h_aa_compatible = -2 * data["f_aaa"] * data["h_a"] / f_aa
h_ab_compatible = (
    (f_aa * data["f_aab"] - 3 * f_ab * data["f_aaa"])
    * data["h_a"]
    / f_aa**2
)
h_bb_compatible = (
    sp.Rational(1, 2)
    + (data["f_abb"] * f_aa - f_ab * data["f_aab"])
    * data["h_a"]
    / f_aa**2
    + f_ab * h_ab_compatible / f_aa
)
assert (
    sp.factor(
        J1.subs(
            {
                data["h_aa"]: h_aa_compatible,
                data["h_ab"]: h_ab_compatible,
                data["h_bb"]: h_bb_compatible,
            }
        )
        - 6 * data["h_a"] * E
    )
    == 0
)

g_aa_compatible = (
    sp.Rational(21, 8) * b * f_aa
    - (2 * data["g_a"] - 1) * data["f_aaa"] / f_aa
)
g_ab_compatible = (
    (data["f_aab"] * f_aa - f_ab * data["f_aaa"])
    * (data["g_a"] - sp.Rational(1, 2))
    / f_aa**2
    + f_ab * g_aa_compatible / f_aa
)
g_bb_compatible = (
    6 * b**2
    + (data["f_abb"] * f_aa - f_ab * data["f_aab"])
    * (data["g_a"] - sp.Rational(1, 2))
    / f_aa**2
    + f_ab * g_ab_compatible / f_aa
)
assert (
    sp.factor(
        J3.subs(
            {
                data["g_aa"]: g_aa_compatible,
                data["g_ab"]: g_ab_compatible,
                data["g_bb"]: g_bb_compatible,
            }
        )
        + sp.Rational(3, 2) * (2 * data["g_a"] - 1) * E
    )
    == 0
)


def main() -> None:
    print("PASS: all four forced normal third derivatives have pole order at most 3")
    print("PASS: their leading numerators are A^3*B, A^2*C*B,")
    print("      A*B*(3*C^2-2*f_aa*kappa), C*B*(C^2-2*f_aa*kappa)")
    print("RESULT: normal-Hessian polynomiality leaves only")
    print("        A=0 and C^2=2*f_aa*kappa")
    print("PASS: B is f_aa times the kernel derivative of L modulo L")
    print("PASS: on B=0, the next n3 numerators are A^2*H1 and A*h_a*H2")
    print("PASS: on A=0, C^2=2*f_aa*kappa the four later poles are K1,...,K4")
    print("PASS: K1,...,K4 are derivative identities in A and C^2/f_aa")
    print("PASS: the h_a=0 line branch forces h_aa=h_ab=0 and h_bb=1/2")
    print("SCOPE: these are necessary divisor-local conditions, not inconsistency")


if __name__ == "__main__":
    main()
