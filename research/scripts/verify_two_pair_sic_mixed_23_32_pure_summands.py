#!/usr/bin/env python3
"""Exact pure-summand theorem for V_(2,3)+V_(3,2).

Under contraction-preserving SL_2,

    V_(2,3) = Sym^5 + Sym^3 + Sym^1.

The first contraction removes Sym^1.  On pure Sym^3, the second
contraction cuts out the rational normal cubic of pure cubes.  On pure
Sym^5, it cuts out the tangential variety of binary quintics L^4 M.

This checker verifies those statements and then inserts the resulting
normal forms into the full two-sided mixture with an arbitrary V_(3,2)
block.  Moments two through four remove every negative-block coefficient
that can participate in an unbounded central contraction.  Explicit
exponent cones prove eventual mixed vanishing.

The case in which both Sym^5 and Sym^3 components of the positive block
are nonzero remains open.
"""

from __future__ import annotations

import json
from math import comb, factorial
from pathlib import Path
import re
import shutil
import subprocess

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_mixed_23_32_pure_summands.json"
)

W, V, Z, Y = sp.symbols("W V Z Y")
PAIRING = W * Z + V * Y
VARIABLES = (W, V, Z, Y)


def contraction(expression: sp.Expr) -> sp.Expr:
    result = sp.Integer(0)
    for (w, v, z, y), coefficient in sp.Poly(
        sp.expand(expression), *VARIABLES
    ).terms():
        if w > z or v > y:
            continue
        result += (
            coefficient
            * sp.Rational(factorial(z), factorial(z - w))
            * sp.Rational(factorial(y), factorial(y - v))
            * Z ** (z - w)
            * Y ** (y - v)
        )
    return sp.expand(result)


def lowering(polynomial: sp.Expr) -> sp.Expr:
    return sp.expand(
        Y * sp.diff(polynomial, Z) - W * sp.diff(polynomial, V)
    )


def component_basis(k: int) -> list[sp.Expr]:
    """Return the Sym^(5-2k) divided-power basis inside V_(2,3)."""

    degree = 5 - 2 * k
    polynomial = sp.expand(
        V ** (2 - k) * Z ** (3 - k) * PAIRING**k
    )
    result: list[sp.Expr] = []
    for order in range(degree + 1):
        result.append(polynomial)
        polynomial = sp.expand(lowering(polynomial) / (order + 1))
    return result


def singular_prime_audit(
    singular: str,
    quadrics: tuple[sp.Expr, ...],
    variables: tuple[sp.Symbol, ...],
) -> dict[str, int | bool]:
    names = ",".join(str(variable) for variable in variables)
    expressions = ",".join(str(item).replace("**", "^") for item in quadrics)
    commands = f"""
LIB "primdec.lib";
ring r=0,({names}),dp;
ideal I={expressions};
ideal G=std(I);
list L=primdecGTZ(I);
ideal P=L[1][2];
ideal GP=std(P);
size(G);
dim(G);
size(L);
reduce(G,GP);
reduce(GP,G);
"""
    completed = subprocess.run(
        [singular, "-q"],
        input=commands,
        text=True,
        capture_output=True,
        check=True,
        timeout=120,
    )
    numeric = [
        int(line)
        for line in completed.stdout.splitlines()
        if re.fullmatch(r"-?[0-9]+", line.strip())
    ]
    assert numeric[:3] == [6, 3, 1]
    remainder_lines = [
        line for line in completed.stdout.splitlines() if line.startswith("_[")
    ]
    assert len(remainder_lines) == 12
    assert all(line.endswith("=0") for line in remainder_lines)
    return {
        "groebner_basis_size": numeric[0],
        "affine_dimension": numeric[1],
        "prime_components": numeric[2],
        "input_equals_prime_component": True,
    }


def main() -> None:
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required for the prime audit")

    bases = {
        "s": component_basis(0),
        "t": component_basis(1),
        "r": component_basis(2),
    }
    monomials = [(i, j) for i in range(3) for j in range(4)]
    columns = []
    for prefix in ("s", "t", "r"):
        for polynomial in bases[prefix]:
            expanded = sp.Poly(polynomial, *VARIABLES)
            columns.append(
                [
                    expanded.coeff_monomial(
                        W ** (2 - i) * V**i * Z ** (3 - j) * Y**j
                    )
                    for i, j in monomials
                ]
            )
    change_of_basis = sp.Matrix(
        12,
        12,
        lambda row, column: columns[column][row],
    )
    assert change_of_basis.det() == -22500

    r0, r1 = sp.symbols("r0 r1")
    R_component = r0 * bases["r"][0] + r1 * bases["r"][1]
    assert contraction(R_component) == 12 * (r0 * Z + r1 * Y)

    # Pure Sym^3: the second contraction is the rational-normal-cubic ideal.
    t = sp.symbols("t0:4")
    T_component = sum(
        coefficient * vector for coefficient, vector in zip(t, bases["t"])
    )
    T_second = contraction(T_component**2)
    t_minors = (
        t[0] * t[2] - t[1] ** 2,
        t[0] * t[3] - t[1] * t[2],
        t[1] * t[3] - t[2] ** 2,
    )
    assert sp.expand(
        T_second
        + 336
        * (
            t_minors[2] * Y**2
            + t_minors[1] * Y * Z
            + t_minors[0] * Z**2
        )
    ) == 0
    t_basis = sp.groebner(t_minors, *t, order="grevlex")
    assert len(t_basis) == 3
    assert all(t_basis.reduce(item)[1] == 0 for item in t_minors)
    lam, alpha, beta = sp.symbols("lam alpha beta")
    cube_parameterization = tuple(
        lam * alpha ** (3 - index) * beta**index
        for index in range(4)
    )
    assert all(
        sp.expand(minor.subs(dict(zip(t, cube_parameterization)))) == 0
        for minor in t_minors
    )
    assert sp.Matrix(cube_parameterization).jacobian(
        (lam, alpha, beta)
    ).subs({lam: 2, alpha: 3, beta: 5}).rank() == 2

    # Pure Sym^5: the second contraction is the tangent-quintic ideal.
    s = sp.symbols("s0:6")
    S_component = sum(
        coefficient * vector for coefficient, vector in zip(s, bases["s"])
    )
    S_second = contraction(S_component**2)
    s_quadrics = (
        s[0] * s[4] - 4 * s[1] * s[3] + 3 * s[2] ** 2,
        s[0] * s[5] - 3 * s[1] * s[4] + 2 * s[2] * s[3],
        s[1] * s[5] - 4 * s[2] * s[4] + 3 * s[3] ** 2,
    )
    assert sp.expand(
        S_second
        - 72
        * (
            s_quadrics[2] * Y**2
            + s_quadrics[1] * Y * Z
            + s_quadrics[0] * Z**2
        )
    ) == 0
    prime_audit = singular_prime_audit(singular, s_quadrics, s)

    a, b, c, d = sp.symbols("a b c d")
    tangent_parameterization = (
        a**4 * c,
        a**3 * (a * d + 4 * b * c) / 5,
        a**2 * b * (2 * a * d + 3 * b * c) / 5,
        a * b**2 * (3 * a * d + 2 * b * c) / 5,
        b**3 * (4 * a * d + b * c) / 5,
        b**4 * d,
    )
    tangent_substitution = dict(zip(s, tangent_parameterization))
    assert all(
        sp.expand(quadratic.subs(tangent_substitution)) == 0
        for quadratic in s_quadrics
    )
    assert sp.Matrix(tangent_parameterization).jacobian(
        (a, b, c, d)
    ).subs({a: 1, b: 2, c: 3, d: 5}).rank() == 3

    # The two nonzero orbit normal forms and their visible one-sided factors.
    sym3_normal = bases["t"][0]
    sym5_power_normal = bases["s"][0]
    sym5_tangent_normal = bases["s"][1]
    assert sp.expand(sym3_normal - V * Z**2 * PAIRING) == 0
    assert sym5_power_normal == V**2 * Z**3
    assert sp.expand(
        sym5_tangent_normal - V * Z**2 * (3 * V * Y - 2 * W * Z)
    ) == 0

    # Relevant V_(3,2) coefficients have delta in {-3,-2,-1}.
    relevant_pairs = [
        (i, j)
        for i in range(4)
        for j in range(3)
        if j - i in (-3, -2, -1)
    ]
    assert relevant_pairs == [
        (1, 0),
        (2, 0),
        (2, 1),
        (3, 0),
        (3, 1),
        (3, 2),
    ]
    u = sp.symbols("u0:6")
    relevant_C = sum(
        coefficient
        * W**i
        * V ** (3 - i)
        * Z**j
        * Y ** (2 - j)
        for coefficient, (i, j) in zip(u, relevant_pairs)
    )

    branch_data = {}
    for name, positive, expected in (
        (
            "Sym3_cube",
            sym3_normal,
            (
                2 * u[1] + 3 * u[4],
                u[3],
                (
                    4 * u[0] * u[3]
                    + 2 * u[1] ** 2
                    + 5 * u[1] * u[4]
                    + 5 * u[2] * u[3]
                    + 10 * u[3] * u[5]
                    + 5 * u[4] ** 2
                ),
            ),
        ),
        (
            "Sym5_tangent",
            sym5_tangent_normal,
            (
                u[1] - u[4],
                u[3],
                (
                    18 * u[0] * u[3]
                    + 9 * u[1] ** 2
                    + 15 * u[1] * u[4]
                    + 15 * u[2] * u[3]
                    + 40 * u[3] * u[5]
                    + 20 * u[4] ** 2
                ),
            ),
        ),
    ):
        mixture = sp.expand(positive + relevant_C)
        second = contraction(mixture**2)
        third = contraction(mixture**3)
        fourth = contraction(mixture**4)
        f2, f3, f4 = expected
        if name == "Sym3_cube":
            assert second == 24 * f2
            assert third == 10368 * f3 * Z
            assert fourth == 190080 * f4
        else:
            assert second == 24 * f2
            assert third == 6912 * f3 * Z
            assert fourth == 34560 * f4
        branch_basis = sp.groebner((f2, f3, f4), *u, order="grevlex")
        assert branch_basis.reduce(u[1] ** 2)[1] == 0
        assert branch_basis.reduce(u[3])[1] == 0
        assert branch_basis.reduce(u[4] ** 2)[1] == 0
        branch_data[name] = {
            "moment_2": str(f2),
            "moment_3_z": str(f3),
            "moment_4": str(f4),
            "groebner_basis": [
                str(item.as_expr()) for item in branch_basis.polys
            ],
            "forced_zero_coefficients": ["u1", "u3", "u4"],
            "residual_cutoff": (
                "m > delta_Q+3*max(0,-epsilon_Q)"
            ),
        }

    # Pure-fifth-power branch: moment two removes the unique delta=-3 corner.
    corner = W**3 * Y**2
    h = sp.symbols("h")
    assert contraction((sym5_power_normal + h * corner) ** 2) == 24 * h

    # Exact local obstruction on the remaining mixed-positive branch.
    x = sp.symbols("x0:12")
    monomial_basis = [
        W ** (2 - i) * V**i * Z ** (3 - j) * Y**j
        for i, j in monomials
    ]
    general_positive = sum(
        coefficient * monomial
        for coefficient, monomial in zip(x, monomial_basis)
    )
    mixed_base = V * Z**2 * (2 * V * Y + 3 * W * Z + 5 * V * Z)

    def coefficient_vector(polynomial: sp.Expr) -> sp.Matrix:
        expanded = sp.Poly(sp.expand(polynomial), *VARIABLES)
        return sp.Matrix(
            [
                expanded.coeff_monomial(monomial)
                for monomial in monomial_basis
            ]
        )

    base_vector = coefficient_vector(mixed_base)
    base_point = {
        variable: base_vector[index] for index, variable in enumerate(x)
    }
    moment_equations = []
    for order in range(1, 5):
        value = sp.Poly(contraction(general_positive**order), Z, Y)
        moment_equations.extend(
            value.coeff_monomial(Z ** (order - index) * Y**index)
            for index in range(order + 1)
        )
    jacobian = sp.Matrix(
        [
            [sp.diff(equation, variable) for variable in x]
            for equation in moment_equations
        ]
    ).subs(base_point)
    assert jacobian.rank() == 6

    incidence_directions = sp.Matrix.hstack(
        coefficient_vector(V**2 * Z**2 * Y),
        coefficient_vector(W * V * Z**3),
        coefficient_vector(V**2 * Z**3),
        coefficient_vector(lowering(mixed_base)),
    )
    assert incidence_directions.rank() == 4
    assert jacobian * incidence_directions == sp.zeros(jacobian.rows, 4)

    excess_1 = W * Z**2 * (3 * V * Y - W * Z) / 3
    excess_2 = (
        Y
        * (20 * V**2 * Y**2 - 97 * V * W * Y * Z + 37 * W**2 * Z**2)
        / 20
    )
    excess_vectors = sp.Matrix.hstack(
        coefficient_vector(excess_1),
        coefficient_vector(excess_2),
    )
    assert jacobian * excess_vectors == sp.zeros(jacobian.rows, 2)
    assert sp.Matrix.hstack(
        incidence_directions, excess_vectors
    ).rank() == 6

    alpha_excess, beta_excess = sp.symbols(
        "alpha_excess beta_excess"
    )
    excess = alpha_excess * excess_1 + beta_excess * excess_2
    second_rhs = []
    for order in range(1, 5):
        value = (
            -comb(order, 2)
            * contraction(excess**2 * mixed_base ** (order - 2))
            if order >= 2
            else sp.Integer(0)
        )
        expanded = sp.Poly(value, Z, Y)
        second_rhs.extend(
            expanded.coeff_monomial(Z ** (order - index) * Y**index)
            for index in range(order + 1)
        )
    second_rhs_vector = sp.Matrix(second_rhs)
    obstruction_polynomials = [
        sp.factor((functional.T * second_rhs_vector)[0])
        for functional in jacobian.T.nullspace()
    ]
    nonzero_obstructions = [
        polynomial
        for polynomial in obstruction_polynomials
        if polynomial != 0
    ]
    obstruction_basis = sp.groebner(
        nonzero_obstructions,
        alpha_excess,
        beta_excess,
        order="grevlex",
    )
    assert len(obstruction_basis) == 1
    assert obstruction_basis.reduce(beta_excess**2)[1] == 0
    assert obstruction_basis.reduce(alpha_excess**2)[1] != 0

    second_correction = W * Y * Z * (V * Y - W * Z) / 10
    third_correction = (
        W
        * Y
        * (
            4 * V * Y**2
            - 169 * V * Y * Z
            - 172 * W * Y * Z
            + 169 * W * Z**2
        )
        / 30280
    )
    for order in range(1, 5):
        second_coefficient = order * contraction(
            second_correction * mixed_base ** (order - 1)
        )
        if order >= 2:
            second_coefficient += comb(order, 2) * contraction(
                excess_1**2 * mixed_base ** (order - 2)
            )
        assert sp.expand(second_coefficient) == 0

        third_coefficient = order * contraction(
            third_correction * mixed_base ** (order - 1)
        )
        if order >= 2:
            third_coefficient += order * (order - 1) * contraction(
                excess_1
                * second_correction
                * mixed_base ** (order - 2)
            )
        if order >= 3:
            third_coefficient += comb(order, 3) * contraction(
                excess_1**3 * mixed_base ** (order - 3)
            )
        assert sp.expand(third_coefficient) == 0

    fourth_rhs = []
    for order in range(1, 5):
        value = sp.Integer(0)
        if order >= 2:
            value += order * (order - 1) * contraction(
                excess_1
                * third_correction
                * mixed_base ** (order - 2)
            )
            value += comb(order, 2) * contraction(
                second_correction**2 * mixed_base ** (order - 2)
            )
        if order >= 3:
            value += 3 * comb(order, 3) * contraction(
                excess_1**2
                * second_correction
                * mixed_base ** (order - 3)
            )
        if order >= 4:
            value += comb(order, 4) * contraction(
                excess_1**4 * mixed_base ** (order - 4)
            )
        expanded = sp.Poly(sp.expand(value), Z, Y)
        fourth_rhs.extend(
            expanded.coeff_monomial(Z ** (order - index) * Y**index)
            for index in range(order + 1)
        )
    fourth_rhs_vector = sp.Matrix(fourth_rhs)
    fourth_obstructions = [
        sp.factor((functional.T * fourth_rhs_vector)[0])
        for functional in jacobian.T.nullspace()
    ]
    assert fourth_obstructions == [
        sp.Rational(-533232, 64345),
        0,
        0,
        sp.Rational(21204864, 64345),
        0,
        0,
        0,
        0,
    ]

    tangent_kernel = sp.Matrix.hstack(*jacobian.nullspace())
    second_parameters = sp.symbols("a0:6")
    general_second_vector = (
        coefficient_vector(second_correction)
        + tangent_kernel * sp.Matrix(second_parameters)
    )
    general_second = sum(
        coefficient * monomial
        for coefficient, monomial in zip(
            general_second_vector,
            monomial_basis,
        )
    )
    general_third_rhs = []
    for order in range(1, 5):
        value = sp.Integer(0)
        if order >= 2:
            value += order * (order - 1) * contraction(
                excess_1
                * general_second
                * mixed_base ** (order - 2)
            )
        if order >= 3:
            value += comb(order, 3) * contraction(
                excess_1**3 * mixed_base ** (order - 3)
            )
        expanded = sp.Poly(sp.expand(value), Z, Y)
        general_third_rhs.extend(
            expanded.coeff_monomial(Z ** (order - index) * Y**index)
            for index in range(order + 1)
        )
    general_third_rhs_vector = sp.Matrix(general_third_rhs)
    general_third_obstructions = [
        sp.factor((functional.T * general_third_rhs_vector)[0])
        for functional in jacobian.T.nullspace()
    ]
    assert general_third_obstructions == [0] * 8
    general_third_vector, third_free_parameters = (
        jacobian.gauss_jordan_solve(-general_third_rhs_vector)
    )
    assert len(third_free_parameters) == 6
    general_third = sum(
        coefficient * monomial
        for coefficient, monomial in zip(
            general_third_vector,
            monomial_basis,
        )
    )
    general_fourth_rhs = []
    for order in range(1, 5):
        value = sp.Integer(0)
        if order >= 2:
            value += order * (order - 1) * contraction(
                excess_1
                * general_third
                * mixed_base ** (order - 2)
            )
            value += comb(order, 2) * contraction(
                general_second**2 * mixed_base ** (order - 2)
            )
        if order >= 3:
            value += 3 * comb(order, 3) * contraction(
                excess_1**2
                * general_second
                * mixed_base ** (order - 3)
            )
        if order >= 4:
            value += comb(order, 4) * contraction(
                excess_1**4 * mixed_base ** (order - 4)
            )
        expanded = sp.Poly(sp.expand(value), Z, Y)
        general_fourth_rhs.extend(
            expanded.coeff_monomial(Z ** (order - index) * Y**index)
            for index in range(order + 1)
        )
    general_fourth_rhs_vector = sp.Matrix(general_fourth_rhs)
    general_fourth_obstructions = [
        sp.factor((functional.T * general_fourth_rhs_vector)[0])
        for functional in jacobian.T.nullspace()
    ]
    fourth_polynomial_1 = (
        198510381 * second_parameters[5] ** 2
        + 6275108 * second_parameters[5]
        + 44436
    )
    fourth_polynomial_2 = (
        1311046029 * second_parameters[5] ** 2
        + 20142212 * second_parameters[5]
        + 220884
    )
    expected_fourth_obstructions = [
        -sp.Rational(12, 64345) * fourth_polynomial_1,
        0,
        0,
        sp.Rational(96, 64345) * fourth_polynomial_2,
        0,
        0,
        0,
        0,
    ]
    assert all(
        sp.expand(actual - expected) == 0
        for actual, expected in zip(
            general_fourth_obstructions,
            expected_fourth_obstructions,
        )
    )
    assert not any(
        parameter in obstruction.free_symbols
        for parameter in third_free_parameters
        for obstruction in general_fourth_obstructions
    )
    assert sp.gcd(fourth_polynomial_1, fourth_polynomial_2) == 1
    fourth_resultant = sp.resultant(
        fourth_polynomial_1,
        fourth_polynomial_2,
        second_parameters[5],
    )
    assert fourth_resultant == 2283980165392458318151680000

    payload = {
        "claim": (
            "SIC(2) holds on V_(2,3)+V_(3,2) whenever the nonzero "
            "positive block lies in one of the pure Sym^5 or Sym^3 "
            "summands and satisfies its pure premise"
        ),
        "claim_boundary": (
            "The positive block with both Sym^5 and Sym^3 components "
            "nonzero remains open"
        ),
        "clebsch_gordan": {
            "decomposition": "Sym^5 + Sym^3 + Sym^1",
            "change_of_basis_determinant": -22500,
            "first_contraction": "12*(r0*Z+r1*Y)",
        },
        "pure_Sym3": {
            "second_contraction_ideal": [str(item) for item in t_minors],
            "geometry": "rational normal cubic of binary cubes L^3",
            "normal_form": str(sym3_normal),
        },
        "pure_Sym5": {
            "second_contraction_ideal": [str(item) for item in s_quadrics],
            "prime_audit": prime_audit,
            "geometry": "tangential variety of binary quintics L^4*M",
            "power_normal_form": str(sym5_power_normal),
            "tangent_normal_form": str(sym5_tangent_normal),
        },
        "two_sided_branches": branch_data,
        "pure_power_branch": {
            "moment_2": "24*h",
            "forced_corner": "W^3*Y^2",
            "multiplier_bound": (
                "If w_Q is the central weight and epsilon_Q=deg_V-deg_Y, "
                "put R=max(0,w_Q-epsilon_Q).  A surviving term has "
                "m<=2*R+max(w_Q,0)."
            ),
        },
        "remaining_mixed_positive_local_gate": {
            "base_point": str(mixed_base),
            "moment_orders": [1, 2, 3, 4],
            "jacobian_rank": 6,
            "tangent_kernel_dimension": 6,
            "incidence_tangent_dimension": 4,
            "excess_directions": [str(excess_1), str(excess_2)],
            "second_order_obstruction_ideal": "(beta_excess^2)",
            "surviving_direction_lifts_through_order": 3,
            "second_correction": str(second_correction),
            "third_correction": str(third_correction),
            "fourth_order_obstruction": {
                "all_second_corrections_parameter_count": 6,
                "all_third_corrections_parameter_count": 6,
                "two_residual_quadratics": [
                    str(fourth_polynomial_1),
                    str(fourth_polynomial_2),
                ],
                "resultant": int(fourth_resultant),
                "conclusion": (
                    "The surviving excess direction has no fourth-order "
                    "lift for moments one through four, even after all "
                    "second- and third-order corrections"
                ),
            },
            "claim_boundary": (
                "This is a local calculation for moments one through "
                "four at one mixed-positive point, not a global "
                "classification of the mixed-positive branch"
            ),
        },
        "status": "exact characteristic-zero theorem",
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    print(
        "PASS V_(2,3)+V_(3,2) pure summands: CG classification, "
        "prime tangent quintics, central elimination, and multiplier cones"
    )


if __name__ == "__main__":
    main()
