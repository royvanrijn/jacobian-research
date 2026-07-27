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

all multiplying L^-3.  Therefore polynomial extendability requires either
B=0, or A=0 together with C=0 or C^2=2*f_aa*kappa.

On the characteristic branch B=0, the next two Laurent numerators of n3
factor as A^2*H1 and A*h_a*H2.  Thus its A != 0 part additionally requires
H1=0 and h_a*H2=0.

On A=C=0, the remaining two pole coefficients of n1 and n0 are nonzero
scalar multiples of kappa*J1, kappa*J2, kappa*J3, and kappa*J4.  Since the
Keller constant kappa is nonzero, this subbranch requires all four J's to
vanish.  Mixed-derivative compatibility then reduces J1 and J3 to multiples
of h_a*E and (2*g_a-1)*E for one common expression E.
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
    print("RESULT: polynomiality requires B=0, or")
    print("        A=0 and (C=0 or C^2=2*f_aa*kappa)")
    print("PASS: B is f_aa times the kernel derivative of L modulo L")
    print("PASS: on B=0, the next n3 numerators are A^2*H1 and A*h_a*H2")
    print("PASS: on A=C=0, the remaining pole numerators force J1=J2=J3=J4=0")
    print("PASS: mixed compatibility reduces this to h_a*E=(2*g_a-1)*E=0")
    print("SCOPE: these are necessary divisor-local conditions, not inconsistency")


if __name__ == "__main__":
    main()
