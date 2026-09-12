#!/usr/bin/env python3
"""Verify the double-conic normal layers and the first decic strata.

Normalize q=x*z-y^2 and use the SL2-equivariant harmonic splitting

    h5 = H5(f10) + q*H3(g6) + q^2*H1(k2).

The Hessian determinant has the unique harmonic expansion

    H9(Phi18) + q*H7(Phi14) + q^2*H5(Phi10)
        + q^3*H3(Phi6) + q^4*H1(Phi2).

This checker verifies the displayed transvectant formulae for the four
divisibility layers Phi18, Phi14, Phi10, Phi6.  Exact Clebsch--Gordan rank
checks make the finite interpolation injective on every cubic covariant
block.  It then uses exact rational modular standard bases in Singular to
exclude every binary-decic root partition with at most three support points,
every four-support partition at harmonic cross-ratio, eight complete
four-support partitions, and the generic point of the ninth partition.
"""

from __future__ import annotations

import argparse
import gc
import random
import shutil
import subprocess
from collections.abc import Iterable

import sympy as sp


x, y, z, s, t = sp.symbols("x y z s t")
lam, mu, nu = sp.symbols("lam mu nu")
q = x * z - y**2


def ternary_monomials(degree: int) -> list[sp.Expr]:
    return [
        x**a * y**b * z ** (degree - a - b)
        for a in range(degree + 1)
        for b in range(degree - a + 1)
    ]


def harmonic_lift_matrix(degree: int) -> tuple[list[sp.Expr], sp.Matrix]:
    """Return the lift S^(2d) -> ker(Box) in S^d(S^2)."""

    monomials = ternary_monomials(degree)
    rows: list[list[sp.Expr]] = []
    for index in range(2 * degree + 1):
        rows.append(
            [
                sp.Poly(
                    monomial.subs({x: s**2, y: s * t, z: t**2}), s, t
                ).coeff_monomial(s ** (2 * degree - index) * t**index)
                for monomial in monomials
            ]
        )
    if degree >= 2:
        for lower_monomial in ternary_monomials(degree - 2):
            rows.append(
                [
                    sp.Poly(
                        sp.diff(monomial, x, z)
                        - sp.diff(monomial, y, 2) / 4,
                        x,
                        y,
                        z,
                    ).coeff_monomial(lower_monomial)
                    for monomial in monomials
                ]
            )
    constraint_matrix = sp.Matrix(rows)
    assert constraint_matrix.rows == constraint_matrix.cols == len(monomials)
    right_hand_side = sp.zeros(len(monomials), 2 * degree + 1)
    right_hand_side[: 2 * degree + 1, :] = sp.eye(2 * degree + 1)
    return monomials, constraint_matrix.inv() * right_hand_side


LIFT_DATA = {
    degree: harmonic_lift_matrix(degree) for degree in (1, 3, 5, 7, 9)
}


def harmonic_lift(coefficients: list[sp.Expr], degree: int) -> sp.Expr:
    monomials, matrix = LIFT_DATA[degree]
    assert len(coefficients) == 2 * degree + 1
    return sp.expand(
        sum(
            monomials[row]
            * sum(
                matrix[row, column] * coefficients[column]
                for column in range(2 * degree + 1)
            )
            for row in range(len(monomials))
        )
    )


def binary_coefficients(polynomial: sp.Expr, degree: int) -> list[sp.Expr]:
    binary = sp.Poly(sp.expand(polynomial), s, t)
    return [
        binary.coeff_monomial(s ** (degree - index) * t**index)
        for index in range(degree + 1)
    ]


def harmonic_layers(
    f_coefficients: list[sp.Expr],
    g_coefficients: list[sp.Expr],
    k_coefficients: list[sp.Expr],
) -> dict[int, sp.Expr]:
    hessian = (
        lam * sp.hessian(harmonic_lift(f_coefficients, 5), (x, y, z))
        + mu * sp.hessian(q * harmonic_lift(g_coefficients, 3), (x, y, z))
        + nu * sp.hessian(q**2 * harmonic_lift(k_coefficients, 1), (x, y, z))
    )
    current = sp.expand(hessian.det())
    layers: dict[int, sp.Expr] = {}
    for ternary_degree in (9, 7, 5, 3):
        binary_degree = 2 * ternary_degree
        layer = sp.expand(current.subs({x: s**2, y: s * t, z: t**2}))
        layers[binary_degree] = layer
        if ternary_degree > 3:
            lifted = harmonic_lift(
                binary_coefficients(layer, binary_degree), ternary_degree
            )
            quotient, remainder = sp.div(
                sp.Poly(sp.expand(current - lifted), z), sp.Poly(q, z)
            )
            assert remainder.as_expr() == 0
            current = sp.expand(quotient.as_expr())
    return layers


def transvectant(first: sp.Expr, second: sp.Expr, order: int) -> sp.Expr:
    return sp.expand(
        sum(
            (-1) ** index
            * sp.binomial(order, index)
            * sp.diff(first, s, order - index, t, index)
            * sp.diff(second, s, index, t, order - index)
            for index in range(order + 1)
        )
    )


# A table row is (paired forms, first transvection, third form,
# second transvection, coefficient).  Thus ("ff",4,"g",2,c) means
# c*((f,f)_4,g)_2.  The unnormalized transvectant convention is used.
TABLE: dict[int, list[tuple[str, int, str, int, sp.Rational]]] = {
    18: [
        ("ff", 4, "f", 2, -sp.Rational(1, 56582064)),
        ("ff", 6, "f", 0, sp.Rational(1, 18860688)),
        ("ff", 2, "g", 2, sp.Rational(1, 12150)),
        ("ff", 4, "g", 0, -sp.Rational(19, 317520)),
        ("ff", 2, "k", 0, sp.Rational(16, 81)),
        ("gg", 0, "f", 2, -sp.Rational(8, 99)),
        ("gg", 2, "f", 0, sp.Rational(6, 55)),
        ("fg", 0, "k", 0, -sp.Rational(80)),
        ("gg", 0, "g", 0, sp.Rational(32)),
    ],
    14: [
        ("ff", 0, "f", 8, -sp.Rational(1, 64818903663360)),
        ("ff", 2, "f", 6, sp.Rational(191, 66292060564800)),
        ("ff", 2, "g", 4, sp.Rational(1, 10024560)),
        ("ff", 4, "g", 2, sp.Rational(1, 39584160)),
        ("ff", 6, "g", 0, -sp.Rational(19, 227026800)),
        ("ff", 2, "k", 2, sp.Rational(5, 4131)),
        ("ff", 4, "k", 0, sp.Rational(1, 95256)),
        ("gg", 0, "f", 4, -sp.Rational(1, 33660)),
        ("gg", 2, "f", 2, -sp.Rational(3, 4675)),
        ("gg", 4, "f", 0, sp.Rational(1, 360)),
        ("fg", 0, "k", 2, -sp.Rational(25, 51)),
        ("fg", 1, "k", 1, -sp.Rational(9, 14)),
        ("fg", 2, "k", 0, -sp.Rational(16, 945)),
        ("kk", 0, "f", 0, -sp.Rational(160)),
        ("gg", 0, "g", 2, -sp.Rational(27, 680)),
        ("gg", 0, "k", 0, sp.Rational(112)),
    ],
    10: [
        ("ff", 0, "f", 10, -sp.Rational(1, 608083813148797440)),
        ("ff", 2, "f", 8, sp.Rational(4399, 1303036742461708800)),
        ("ff", 2, "g", 6, sp.Rational(29, 113837724000)),
        ("ff", 4, "g", 4, sp.Rational(53, 367783416000)),
        ("ff", 6, "g", 2, -sp.Rational(31, 619783164000)),
        ("ff", 8, "g", 0, -sp.Rational(1, 3960744480)),
        ("ff", 4, "k", 2, sp.Rational(109, 136216080)),
        ("ff", 6, "k", 0, -sp.Rational(17, 227026800)),
        ("gg", 0, "f", 6, -sp.Rational(1, 68108040)),
        ("gg", 2, "f", 4, -sp.Rational(1, 126126000)),
        ("gg", 4, "f", 2, -sp.Rational(467, 26535600)),
        ("gg", 6, "f", 0, sp.Rational(29, 75600)),
        ("fg", 2, "k", 2, sp.Rational(62, 19305)),
        ("fg", 3, "k", 1, -sp.Rational(1, 975)),
        ("fg", 4, "k", 0, -sp.Rational(263, 300300)),
        ("kk", 0, "f", 2, -sp.Rational(32, 39)),
        ("kk", 2, "f", 0, -sp.Rational(40)),
        ("gg", 0, "g", 4, -sp.Rational(31, 409500)),
        ("gg", 0, "k", 2, sp.Rational(4, 143)),
        ("gg", 2, "k", 0, -sp.Rational(18, 275)),
        ("kk", 0, "g", 0, sp.Rational(144)),
    ],
    6: [
        ("ff", 2, "f", 10, sp.Rational(1, 77203623553056000)),
        ("ff", 4, "f", 8, sp.Rational(277, 19455313135370112000)),
        ("ff", 4, "g", 6, sp.Rational(101, 101949562915200)),
        ("ff", 6, "g", 4, sp.Rational(19, 35399153790000)),
        ("ff", 8, "g", 2, -sp.Rational(1, 668375631000)),
        ("ff", 10, "g", 0, -sp.Rational(1, 108020304000)),
        ("ff", 6, "k", 2, sp.Rational(19, 28605376800)),
        ("ff", 8, "k", 0, -sp.Rational(53, 118822334400)),
        ("gg", 0, "f", 8, -sp.Rational(41, 2022808788000)),
        ("gg", 2, "f", 6, sp.Rational(397, 374594220000)),
        ("gg", 4, "f", 4, sp.Rational(377, 2829103200)),
        ("fg", 4, "k", 2, sp.Rational(23, 7567560)),
        ("fg", 5, "k", 1, -sp.Rational(1, 249480)),
        ("fg", 6, "k", 0, -sp.Rational(829, 78586200)),
        ("kk", 0, "f", 4, -sp.Rational(13, 4158)),
        ("gg", 0, "g", 6, -sp.Rational(1747, 5096520000)),
        ("gg", 2, "g", 4, sp.Rational(211, 217800000)),
        ("gg", 2, "k", 2, sp.Rational(8, 1925)),
        ("gg", 4, "k", 0, -sp.Rational(307, 113400)),
        ("kk", 0, "g", 2, -sp.Rational(4, 15)),
        ("kk", 2, "g", 0, sp.Rational(8)),
        ("kk", 0, "k", 0, sp.Rational(96)),
    ],
}


def covariant_term(
    pair: str,
    first_order: int,
    third: str,
    second_order: int,
    forms: dict[str, sp.Expr],
) -> sp.Expr:
    paired = transvectant(forms[pair[0]], forms[pair[1]], first_order)
    return transvectant(paired, forms[third], second_order)


def term_multidegree(pair: str, third: str) -> tuple[int, int, int]:
    letters = pair + third
    return letters.count("f"), letters.count("g"), letters.count("k")


def formula_component(
    degree: int,
    multidegree: tuple[int, int, int],
    forms: dict[str, sp.Expr],
) -> sp.Expr:
    return sp.expand(
        sum(
            coefficient
            * covariant_term(pair, first_order, third, second_order, forms)
            for pair, first_order, third, second_order, coefficient in TABLE[degree]
            if term_multidegree(pair, third) == multidegree
        )
    )


MULTIDEGREES = [
    (3, 0, 0),
    (2, 1, 0),
    (2, 0, 1),
    (1, 2, 0),
    (1, 1, 1),
    (1, 0, 2),
    (0, 3, 0),
    (0, 2, 1),
    (0, 1, 2),
    (0, 0, 3),
]

EXPECTED_RANKS = {
    18: [2, 3, 2, 2, 1, 0, 1, 0, 0, 0],
    14: [2, 4, 2, 3, 3, 1, 1, 1, 0, 0],
    10: [2, 4, 2, 4, 3, 2, 1, 2, 1, 0],
    6: [2, 4, 2, 3, 3, 1, 2, 2, 2, 1],
}


def natural_candidates(
    degree: int,
    multidegree: tuple[int, int, int],
    forms: dict[str, sp.Expr],
) -> list[sp.Expr]:
    f_count, g_count, k_count = multidegree
    total_order = 10 * f_count + 6 * g_count + 2 * k_count
    if total_order < degree:
        return []
    total_contraction = (total_order - degree) // 2
    orders = {"f": 10, "g": 6, "k": 2}

    if f_count >= 2:
        paired = "f"
        third = "f" if f_count == 3 else ("g" if g_count else "k")
    elif g_count >= 2:
        paired = "g"
        third = "g" if g_count == 3 else ("f" if f_count else "k")
    elif k_count >= 2:
        paired = "k"
        third = "k" if k_count == 3 else ("f" if f_count else "g")
    else:
        candidates = []
        for first_order in range(7):
            second_order = total_contraction - first_order
            intermediate_order = 16 - 2 * first_order
            if 0 <= second_order <= min(intermediate_order, 2):
                candidates.append(
                    covariant_term(
                        "fg", first_order, "k", second_order, forms
                    )
                )
        return candidates

    candidates = []
    for first_order in range(0, orders[paired] + 1, 2):
        second_order = total_contraction - first_order
        intermediate_order = 2 * orders[paired] - 2 * first_order
        if 0 <= second_order <= min(intermediate_order, orders[third]):
            candidates.append(
                covariant_term(
                    paired + paired,
                    first_order,
                    third,
                    second_order,
                    forms,
                )
            )
    return candidates


def verify_covariant_layers(degrees: tuple[int, ...] = (18, 14, 10, 6)) -> None:
    random_generator = random.Random(73)
    samples: list[
        tuple[dict[str, sp.Expr], dict[int, sp.Expr]]
    ] = []
    for _ in range(2):
        f_coefficients = [
            sp.Integer(random_generator.randint(-2, 2)) for _ in range(11)
        ]
        g_coefficients = [
            sp.Integer(random_generator.randint(-2, 2)) for _ in range(7)
        ]
        k_coefficients = [
            sp.Integer(random_generator.randint(-2, 2)) for _ in range(3)
        ]
        forms = {
            "f": sum(
                f_coefficients[index] * s ** (10 - index) * t**index
                for index in range(11)
            ),
            "g": sum(
                g_coefficients[index] * s ** (6 - index) * t**index
                for index in range(7)
            ),
            "k": sum(
                k_coefficients[index] * s ** (2 - index) * t**index
                for index in range(3)
            ),
        }
        samples.append(
            (forms, harmonic_layers(f_coefficients, g_coefficients, k_coefficients))
        )

    for degree in degrees:
        for multidegree, expected_rank in zip(
            MULTIDEGREES, EXPECTED_RANKS[degree]
        ):
            evaluation_rows: list[list[sp.Expr]] = []
            for forms, layers in samples:
                actual = sp.Poly(layers[degree], lam, mu, nu).coeff_monomial(
                    lam ** multidegree[0]
                    * mu ** multidegree[1]
                    * nu ** multidegree[2]
                )
                predicted = formula_component(degree, multidegree, forms)
                assert sp.expand(actual - predicted) == 0

                candidates = natural_candidates(degree, multidegree, forms)
                coefficient_columns = [
                    binary_coefficients(candidate, degree)
                    for candidate in candidates
                ]
                for index in range(degree + 1):
                    evaluation_rows.append(
                        [column[index] for column in coefficient_columns]
                    )
            if expected_rank == 0:
                assert not evaluation_rows or not evaluation_rows[0]
            else:
                assert sp.Matrix(evaluation_rows).rank() == expected_rank
        print(f"PASS: Phi_{degree} transvectant identity and covariant ranks")


def canonical_lift(binary_form: sp.Expr) -> sp.Expr:
    """A convenient lift; changing it only translates the general q*G3."""

    polynomial = sp.Poly(sp.expand(binary_form), s, t)
    result = 0
    for index in range(11):
        coefficient = polynomial.coeff_monomial(s ** (10 - index) * t**index)
        half = index // 2
        if index % 2 == 0:
            monomial = x ** (5 - half) * z**half
        else:
            monomial = x ** (4 - half) * y * z**half
        result += coefficient * monomial
    return sp.expand(result)


g_coefficients = sp.symbols("g0:10")
line_x, line_y, line_z, inverse = sp.symbols("line_x line_y line_z inv")
cross_ratio = sp.symbols("cr")
cubic_monomials = [
    x**a * y**b * z ** (3 - a - b)
    for a in range(4)
    for b in range(4 - a)
]
general_cubic = sum(
    coefficient * monomial
    for coefficient, monomial in zip(g_coefficients, cubic_monomials)
)
singular_variables = g_coefficients + (line_x, line_y, line_z, inverse)


def singular_polynomial(expression: sp.Expr) -> str:
    _, polynomial = sp.Poly(
        expression, *singular_variables, domain=sp.QQ
    ).clear_denoms(convert=True)
    return str(sp.expand(polynomial.as_expr())).replace("**", "^")


def verify_decic(binary_form: sp.Expr, name: str) -> None:
    h5 = canonical_lift(binary_form) + q * general_cubic
    target = q**4 * (line_x * x + line_y * y + line_z * z)
    equations = [
        coefficient
        for _, coefficient in sp.Poly(
            sp.expand(sp.hessian(h5, (x, y, z)).det() - target), x, y, z
        ).terms()
    ]
    equation_source = ",".join(
        singular_polynomial(equation) for equation in equations
    )
    for chart in (line_x, line_y, line_z):
        source = "\n".join(
            [
                'LIB "modstd.lib";',
                f"ring r=0,({','.join(map(str, singular_variables))}),dp;",
                "option(redSB);",
                "ideal I="
                + equation_source
                + ","
                + singular_polynomial(inverse * chart - 1)
                + ";",
                "ideal J=modStd(I);",
                'print("NONZERO_LINE_CHART");',
                "J;",
                "exit;",
            ]
        )
        completed = subprocess.run(
            ["Singular", "-q"],
            input=source,
            text=True,
            capture_output=True,
            check=True,
            timeout=300,
        )
        assert "NONZERO_LINE_CHART\nJ[1]=1" in completed.stdout
    print(f"PASS: {name} has unit saturation in all residual-line charts")


def verify_generic_four_root(
    partition: tuple[int, int, int, int],
    charts: tuple[sp.Symbol, ...] = (line_x, line_y, line_z),
) -> None:
    """Test one complete four-root cross-ratio family."""

    first, second, third, fourth = partition
    binary_form = (
        s**first
        * t**second
        * (s - t) ** third
        * (s - cross_ratio * t) ** fourth
    )
    h5 = canonical_lift(binary_form) + q * general_cubic
    target = q**4 * (line_x * x + line_y * y + line_z * z)
    equations = [
        coefficient
        for _, coefficient in sp.Poly(
            sp.expand(sp.hessian(h5, (x, y, z)).det() - target), x, y, z
        ).terms()
    ]
    parameter_variables = (cross_ratio,) + singular_variables

    def parameter_polynomial(expression: sp.Expr) -> str:
        _, polynomial = sp.Poly(
            expression, *parameter_variables, domain=sp.QQ
        ).clear_denoms(convert=True)
        return str(sp.expand(polynomial.as_expr())).replace("**", "^")

    equation_source = ",".join(
        parameter_polynomial(equation) for equation in equations
    )
    for chart in charts:
        localization = inverse * chart * cross_ratio * (cross_ratio - 1) - 1
        source = "\n".join(
            [
                'LIB "modstd.lib";',
                f"ring r=0,({','.join(map(str, parameter_variables))}),dp;",
                "option(redSB);",
                "ideal I="
                + equation_source
                + ","
                + parameter_polynomial(localization)
                + ";",
                "ideal J=modStd(I);",
                'print("GENERIC_CROSS_RATIO_CHART");',
                "J;",
                "exit;",
            ]
        )
        completed = subprocess.run(
            ["Singular", "-q"],
            input=source,
            text=True,
            capture_output=True,
            check=True,
            timeout=1200,
        )
        assert "GENERIC_CROSS_RATIO_CHART\nJ[1]=1" in completed.stdout
        print(f"PASS: generic partition {partition}, chart {chart}, is a unit")


def verify_balanced_function_fields() -> None:
    """Prove generic emptiness of the two balanced four-root rows.

    This works over QQ(cross_ratio).  It does not identify the finite
    exceptional cross-ratio locus left by denominators of the certificate.
    """

    parameter_variables = singular_variables
    for partition in ((4, 2, 2, 2), (3, 3, 2, 2)):
        first, second, third, fourth = partition
        binary_form = (
            s**first
            * t**second
            * (s - t) ** third
            * (s - cross_ratio * t) ** fourth
        )
        h5 = canonical_lift(binary_form) + q * general_cubic
        target = q**4 * (line_x * x + line_y * y + line_z * z)
        equations = [
            coefficient
            for _, coefficient in sp.Poly(
                sp.expand(sp.hessian(h5, (x, y, z)).det() - target), x, y, z
            ).terms()
        ]

        def function_field_polynomial(expression: sp.Expr) -> str:
            return str(sp.expand(expression)).replace("**", "^")

        equation_source = ",".join(
            function_field_polynomial(equation) for equation in equations
        )
        for chart in (line_x, line_y, line_z):
            source = "\n".join(
                [
                    'LIB "ffmodstd.lib";',
                    f"ring r=(0,{cross_ratio}),"
                    + f"({','.join(map(str, parameter_variables))}),dp;",
                    "option(redSB);",
                    "ideal I="
                    + equation_source
                    + ","
                    + function_field_polynomial(inverse * chart - 1)
                    + ";",
                    "ideal J=ffmodStd(I);",
                    'print("FUNCTION_FIELD_CHART");',
                    "J;",
                    "exit;",
                ]
            )
            completed = subprocess.run(
                ["Singular", "-q"],
                input=source,
                text=True,
                capture_output=True,
                check=True,
                timeout=1200,
            )
            assert "FUNCTION_FIELD_CHART\nJ[1]=1" in completed.stdout, (
                completed.stdout + completed.stderr
            )
            print(
                f"PASS: QQ(cross_ratio) partition {partition}, "
                f"chart {chart}, is a unit"
            )


def partition_form(partition: Iterable[int], fourth_root: int | None = None) -> sp.Expr:
    parts = list(partition)
    if len(parts) == 1:
        return s**parts[0]
    if len(parts) == 2:
        return s ** parts[0] * t ** parts[1]
    if len(parts) == 3:
        return s ** parts[0] * t ** parts[1] * (s - t) ** parts[2]
    assert len(parts) == 4 and fourth_root is not None
    return (
        s ** parts[0]
        * t ** parts[1]
        * (s - t) ** parts[2]
        * (s - fourth_root * t) ** parts[3]
    )


SUPPORT_ONE_TWO = [
    (10,),
    (9, 1),
    (8, 2),
    (7, 3),
    (6, 4),
    (5, 5),
]

SUPPORT_THREE = [
    (8, 1, 1),
    (7, 2, 1),
    (6, 3, 1),
    (6, 2, 2),
    (5, 4, 1),
    (5, 3, 2),
    (4, 4, 2),
    (4, 3, 3),
]

HARMONIC_FOUR_SUPPORT = [
    (7, 1, 1, 1),
    (6, 2, 1, 1),
    (5, 3, 1, 1),
    (5, 2, 2, 1),
    (4, 4, 1, 1),
    (4, 3, 2, 1),
    (4, 2, 2, 2),
    (3, 3, 3, 1),
    (3, 3, 2, 2),
]

GENERIC_FOUR_GROUPS = {
    "generic-four-7111": (7, 1, 1, 1),
    "generic-four-6211": (6, 2, 1, 1),
    "generic-four-5311": (5, 3, 1, 1),
    "generic-four-5221": (5, 2, 2, 1),
    "generic-four-4411": (4, 4, 1, 1),
    "generic-four-4321": (4, 3, 2, 1),
    "generic-four-3331": (3, 3, 3, 1),
}


def verify_root_strata(group: str) -> None:
    assert shutil.which("Singular") is not None, "Singular is required"
    if group == "support-one-two":
        partitions = SUPPORT_ONE_TWO
    elif group == "support-three":
        partitions = SUPPORT_THREE
    elif group == "harmonic-four":
        partitions = HARMONIC_FOUR_SUPPORT
    else:
        raise AssertionError(f"unknown root-stratum group: {group}")

    for partition in partitions:
        if group == "harmonic-four":
            binary_form = partition_form(partition, fourth_root=2)
            name = f"harmonic four-root partition {partition}"
        else:
            binary_form = partition_form(partition)
            name = f"partition {partition}"
        verify_decic(binary_form, name)
        sp.core.cache.clear_cache()
        gc.collect()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--group",
        choices=(
            "all",
            "layers",
            "layer-18",
            "layer-14",
            "layer-10",
            "layer-6",
            "support-one-two",
            "support-three",
            "harmonic-four",
            "balanced-function-fields",
            "generic-four-4222-root-chart",
            *GENERIC_FOUR_GROUPS,
        ),
        default="all",
    )
    arguments = parser.parse_args()

    if arguments.group in ("all", "layers"):
        verify_covariant_layers()
        print("THEOREM: four normal layers are the displayed binary covariants")
    elif arguments.group.startswith("layer-"):
        degree = int(arguments.group.removeprefix("layer-"))
        verify_covariant_layers((degree,))
        print(f"THEOREM: Phi_{degree} is the displayed binary covariant")
    if arguments.group == "all":
        for group in ("support-one-two", "support-three", "harmonic-four"):
            verify_root_strata(group)
        for partition in GENERIC_FOUR_GROUPS.values():
            verify_generic_four_root(partition)
        verify_generic_four_root((4, 2, 2, 2), (line_x,))
        verify_balanced_function_fields()
    elif arguments.group in GENERIC_FOUR_GROUPS:
        verify_generic_four_root(GENERIC_FOUR_GROUPS[arguments.group])
    elif arguments.group == "balanced-function-fields":
        verify_balanced_function_fields()
    elif arguments.group == "generic-four-4222-root-chart":
        verify_generic_four_root((4, 2, 2, 2), (line_x,))
    elif arguments.group not in (
        "layers",
        "layer-18",
        "layer-14",
        "layer-10",
        "layer-6",
        "balanced-function-fields",
        "generic-four-4222-root-chart",
        *GENERIC_FOUR_GROUPS,
    ):
        verify_root_strata(arguments.group)
    if arguments.group == "all":
        print("THEOREM: support <= 3 and harmonic four-support decics are empty")
        print("THEOREM: eight complete four-support partitions are empty")
        print(
            "FRONTIER: only a finite exceptional cross-ratio locus in "
            "partition (3,3,2,2) remains among four-support decics"
        )
    elif arguments.group == "harmonic-four":
        print("THEOREM: all harmonic four-support decics are empty")
    elif arguments.group in GENERIC_FOUR_GROUPS:
        partition = GENERIC_FOUR_GROUPS[arguments.group]
        print(f"THEOREM: the generic four-support partition {partition} is empty")
    elif arguments.group == "balanced-function-fields":
        print("THEOREM: both balanced four-root rows are empty over QQ(lambda)")
        print("FRONTIER: only finite exceptional cross-ratio loci remain")
    elif arguments.group == "generic-four-4222-root-chart":
        print("THEOREM: the (4,2,2,2) double-root value chart is empty")


if __name__ == "__main__":
    main()
