#!/usr/bin/env python3
"""Exact generic local gate on the mixed-positive V_(2,3) branch.

The incidence family

    V*Z^2*(a*V*Y + b*W*Z + c*V*Z)

has one essential parameter on the chart b*(a-b) != 0.  A
contraction-preserving raising transformation removes c, and scaling
sets b=1.  This checker works over Q(u), where u=a/b.

For the first four full polynomial-valued contractions it:

* constructs the exact 14-by-12 moment Jacobian;
* separates the four incidence tangents from two excess directions;
* proves that the second excess is obstructed at order two;
* lets the first excess carry all six second- and all six third-order
  correction parameters;
* computes the two order-four compatibility quadratics and their
  resultant in Q[u].

The conclusion is a generic local jet exclusion.  It is not a global
classification of the mixed-positive locus and not an SIC theorem.
"""

from __future__ import annotations

import json
from math import comb
from pathlib import Path

import sympy as sp

from verify_two_pair_sic_mixed_23_32_pure_summands import (
    V,
    VARIABLES,
    W,
    Y,
    Z,
    contraction,
    lowering,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_mixed_23_32_generic_local_gate.json"
)


def main() -> None:
    u = sp.symbols("u")
    monomials = [
        W ** (2 - i) * V**i * Z ** (3 - j) * Y**j
        for i in range(3)
        for j in range(4)
    ]

    def coefficient_vector(polynomial: sp.Expr) -> sp.Matrix:
        expanded = sp.Poly(sp.expand(polynomial), *VARIABLES)
        return sp.Matrix(
            [expanded.coeff_monomial(monomial) for monomial in monomials]
        )

    def output_vector(values: list[sp.Expr]) -> sp.Matrix:
        result: list[sp.Expr] = []
        for order, value in enumerate(values, start=1):
            expanded = sp.Poly(sp.cancel(sp.expand(value)), Z, Y)
            result.extend(
                expanded.coeff_monomial(
                    Z ** (order - index) * Y**index
                )
                for index in range(order + 1)
            )
        return sp.Matrix(result)

    def raise_once(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            Z * sp.diff(polynomial, Y) - V * sp.diff(polynomial, W)
        )

    a, b, c, shift = sp.symbols("a b c shift")
    incidence_general = V * Z**2 * (
        a * V * Y + b * W * Z + c * V * Z
    )
    assert sp.expand(
        raise_once(incidence_general) - (a - b) * V**2 * Z**3
    ) == 0
    assert raise_once(raise_once(incidence_general)) == 0
    raised = sp.expand(
        incidence_general + shift * raise_once(incidence_general)
    )
    assert sp.cancel(
        raised.subs(shift, -c / (a - b))
        - V * Z**2 * (a * V * Y + b * W * Z)
    ) == 0

    base = V * Z**2 * (u * V * Y + W * Z)
    coefficients = sp.symbols("x0:12")
    general = sum(
        coefficient * monomial
        for coefficient, monomial in zip(coefficients, monomials)
    )
    base_vector = coefficient_vector(base)
    base_point = {
        coefficient: base_vector[index]
        for index, coefficient in enumerate(coefficients)
    }

    moment_equations: list[sp.Expr] = []
    for order in range(1, 5):
        value = sp.Poly(contraction(general**order), Z, Y)
        moment_equations.extend(
            value.coeff_monomial(Z ** (order - index) * Y**index)
            for index in range(order + 1)
        )
    jacobian = sp.Matrix(
        [
            [
                sp.diff(equation, coefficient)
                for coefficient in coefficients
            ]
            for equation in moment_equations
        ]
    ).subs(base_point)
    assert jacobian.shape == (14, 12)
    assert jacobian.rank() == 6

    incidence_directions = sp.Matrix.hstack(
        coefficient_vector(V**2 * Y * Z**2),
        coefficient_vector(W * V * Z**3),
        coefficient_vector(V**2 * Z**3),
        coefficient_vector(lowering(base)),
    )
    assert incidence_directions.rank() == 4
    assert all(
        sp.cancel(entry) == 0
        for entry in jacobian * incidence_directions
    )

    excess_1 = W * Z**2 * (3 * V * Y - W * Z) / 3
    excess_2 = (
        Y
        * (
            (u + 6) * V**2 * Y**2
            - (8 * u + 27) * V * W * Y * Z
            + (5 * u + 9) * W**2 * Z**2
        )
        / (u + 6)
    )
    excess_vectors = sp.Matrix.hstack(
        coefficient_vector(excess_1),
        coefficient_vector(excess_2),
    )
    assert all(
        sp.cancel(entry) == 0 for entry in jacobian * excess_vectors
    )
    assert sp.Matrix.hstack(
        incidence_directions,
        excess_vectors,
    ).rank() == 6

    left_kernel = jacobian.T.nullspace()
    assert len(left_kernel) == 8
    alpha, beta = sp.symbols("alpha beta")
    excess = alpha * excess_1 + beta * excess_2
    second_rhs = output_vector(
        [
            (
                -comb(order, 2)
                * contraction(excess**2 * base ** (order - 2))
                if order >= 2
                else sp.Integer(0)
            )
            for order in range(1, 5)
        ]
    )
    second_obstructions = [
        sp.factor(sp.cancel((functional.T * second_rhs)[0]))
        for functional in left_kernel
    ]
    second_polynomial_1 = (
        8 * u**4 + 44 * u**3 - 12 * u**2 - 189 * u - 243
    )
    second_polynomial_2 = (
        44 * u**5
        + 426 * u**4
        + 1612 * u**3
        + 2388 * u**2
        + 567 * u
        - 1215
    )
    second_denominator = (u + 6) ** 2 * (4 * u + 3)
    assert second_obstructions == [
        -8640 * beta**2 * second_polynomial_1 / second_denominator,
        0,
        0,
        -69120 * beta**2 * second_polynomial_2 / second_denominator,
        0,
        0,
        0,
        0,
    ]
    second_resultant = sp.resultant(
        second_polynomial_1,
        second_polynomial_2,
        u,
    )
    assert second_resultant == -10554055205970310272

    # Carry the surviving excess_1 direction with every correction.
    second_nonlinear = output_vector(
        [
            (
                comb(order, 2)
                * contraction(excess_1**2 * base ** (order - 2))
                if order >= 2
                else sp.Integer(0)
            )
            for order in range(1, 5)
        ]
    )
    second_vector, second_free = jacobian.gauss_jordan_solve(
        -second_nonlinear
    )
    assert len(second_free) == 6
    second_correction = sp.cancel(
        sum(
            coefficient * monomial
            for coefficient, monomial in zip(second_vector, monomials)
        )
    )

    third_nonlinear = output_vector(
        [
            (
                order
                * (order - 1)
                * contraction(
                    excess_1
                    * second_correction
                    * base ** (order - 2)
                )
                if order >= 2
                else sp.Integer(0)
            )
            + (
                comb(order, 3)
                * contraction(excess_1**3 * base ** (order - 3))
                if order >= 3
                else sp.Integer(0)
            )
            for order in range(1, 5)
        ]
    )
    third_obstructions = [
        sp.factor(sp.cancel((functional.T * third_nonlinear)[0]))
        for functional in left_kernel
    ]
    assert third_obstructions == [0] * 8
    third_vector, third_free = jacobian.gauss_jordan_solve(
        -third_nonlinear
    )
    assert len(third_free) == 6
    third_correction = sp.cancel(
        sum(
            coefficient * monomial
            for coefficient, monomial in zip(third_vector, monomials)
        )
    )

    fourth_nonlinear = output_vector(
        [
            (
                order
                * (order - 1)
                * contraction(
                    excess_1
                    * third_correction
                    * base ** (order - 2)
                )
                if order >= 2
                else sp.Integer(0)
            )
            + (
                comb(order, 2)
                * contraction(
                    second_correction**2 * base ** (order - 2)
                )
                if order >= 2
                else sp.Integer(0)
            )
            + (
                3
                * comb(order, 3)
                * contraction(
                    excess_1**2
                    * second_correction
                    * base ** (order - 3)
                )
                if order >= 3
                else sp.Integer(0)
            )
            + (
                comb(order, 4)
                * contraction(excess_1**4 * base ** (order - 4))
                if order >= 4
                else sp.Integer(0)
            )
            for order in range(1, 5)
        ]
    )
    fourth_obstructions = [
        sp.factor(sp.cancel((functional.T * fourth_nonlinear)[0]))
        for functional in left_kernel
    ]

    correction_parameter = second_free[5]
    fourth_polynomial_1 = (
        144 * correction_parameter**2 * u**7
        + 1512 * correction_parameter**2 * u**6
        + 5256 * correction_parameter**2 * u**5
        + 4482 * correction_parameter**2 * u**4
        - 20088 * correction_parameter**2 * u**3
        - 58563 * correction_parameter**2 * u**2
        - 61236 * correction_parameter**2 * u
        - 19683 * correction_parameter**2
        - 56 * correction_parameter * u**6
        - 282 * correction_parameter * u**5
        - 576 * correction_parameter * u**4
        - 4833 * correction_parameter * u**3
        - 9234 * correction_parameter * u**2
        - 4374 * correction_parameter * u
        + 12 * u**5
        + 18 * u**4
        - 270 * u**3
        - 243 * u**2
    )
    fourth_polynomial_2 = (
        792 * correction_parameter**2 * u**8
        + 11628 * correction_parameter**2 * u**7
        + 75672 * correction_parameter**2 * u**6
        + 272142 * correction_parameter**2 * u**5
        + 564300 * correction_parameter**2 * u**4
        + 611064 * correction_parameter**2 * u**3
        + 191241 * correction_parameter**2 * u**2
        - 183708 * correction_parameter**2 * u
        - 98415 * correction_parameter**2
        - 248 * correction_parameter * u**7
        + 1022 * correction_parameter * u**6
        + 22890 * correction_parameter * u**5
        + 70452 * correction_parameter * u**4
        + 52785 * correction_parameter * u**3
        - 19926 * correction_parameter * u**2
        - 21870 * correction_parameter * u
        + 222 * u**6
        + 1704 * u**5
        + 2952 * u**4
        + 54 * u**3
        - 1215 * u**2
    )
    chart_cubic = 2 * u**3 + 10 * u**2 + 21 * u + 9
    fourth_denominator = (
        (u + 6) ** 2 * (4 * u + 3) * chart_cubic
    )
    expected_fourth = [
        960 * fourth_polynomial_1 / fourth_denominator,
        0,
        0,
        7680 * fourth_polynomial_2 / fourth_denominator,
        0,
        0,
        0,
        0,
    ]
    fourth_differences = [
        sp.factor(sp.cancel(actual - expected))
        for actual, expected in zip(
            fourth_obstructions,
            expected_fourth,
        )
    ]
    assert fourth_differences == [0] * 8
    assert not any(
        parameter in obstruction.free_symbols
        for parameter in third_free
        for obstruction in fourth_obstructions
    )

    exceptional_polynomial = (
        14138 * u**6
        + 142955 * u**5
        + 483945 * u**4
        + 727020 * u**3
        + 540270 * u**2
        + 185004 * u
        + 18225
    )
    fourth_resultant = sp.factor(
        sp.resultant(
            fourth_polynomial_1,
            fourth_polynomial_2,
            correction_parameter,
        )
    )
    assert fourth_resultant == (
        648
        * u**8
        * (u + 6) ** 4
        * (4 * u + 3) ** 2
        * chart_cubic**2
        * exceptional_polynomial
    )
    assert sp.degree(exceptional_polynomial, u) == 6
    assert sp.Poly(chart_cubic, u).is_irreducible
    assert sp.Poly(exceptional_polynomial, u).is_irreducible
    exceptional_correction_coefficient = (
        20 * u**7
        + 1262 * u**6
        + 10842 * u**5
        + 31626 * u**4
        + 38295 * u**3
        + 21168 * u**2
        + 6075 * u
        + 2187
    )
    exceptional_correction_constant = (
        52 * u**6
        + 378 * u**5
        + 581 * u**4
        + 42 * u**3
        + 27 * u**2
        + 243 * u
    )
    fourth_subresultants = sp.subresultants(
        fourth_polynomial_1,
        fourth_polynomial_2,
        correction_parameter,
    )
    assert [sp.degree(item, correction_parameter) for item in fourth_subresultants] == [
        2,
        2,
        1,
        0,
    ]
    expected_linear_subresultant = (
        54
        * u**2
        * (u + 6)
        * (4 * u + 3)
        * chart_cubic
        * (
            exceptional_correction_coefficient * correction_parameter
            + exceptional_correction_constant
        )
    )
    assert sp.expand(
        fourth_subresultants[-2] - expected_linear_subresultant
    ) == 0
    assert sp.gcd(
        exceptional_correction_coefficient,
        exceptional_polynomial,
    ) == 1
    exceptional_correction_resultant = sp.resultant(
        exceptional_correction_coefficient,
        exceptional_polynomial,
        u,
    )
    assert (
        exceptional_correction_resultant
        == 98518095894778815317044086562668309937383692646528
    )

    def direct_fourth_gate(
        label: str,
        special_base: sp.Expr,
    ) -> dict[str, object]:
        """Recompute a rational exceptional point without Q(u) charts."""

        special_vector = coefficient_vector(special_base)
        special_point = {
            coefficient: special_vector[index]
            for index, coefficient in enumerate(coefficients)
        }
        # Rebuild rather than specialize the generic Jacobian when the
        # special point is not on the affine u-chart.
        special_jacobian = sp.Matrix(
            [
                [
                    sp.diff(equation, coefficient)
                    for coefficient in coefficients
                ]
                for equation in moment_equations
            ]
        ).subs(special_point)
        assert special_jacobian.rank() == 6

        special_incidence = sp.Matrix.hstack(
            coefficient_vector(V**2 * Y * Z**2),
            coefficient_vector(W * V * Z**3),
            coefficient_vector(V**2 * Z**3),
            coefficient_vector(lowering(special_base)),
        )
        assert special_incidence.rank() == 4
        special_kernel = special_jacobian.nullspace()
        combined = special_incidence
        special_excess_vectors: list[sp.Matrix] = []
        for vector in special_kernel:
            if combined.row_join(vector).rank() > combined.rank():
                special_excess_vectors.append(vector)
                combined = combined.row_join(vector)
        assert len(special_excess_vectors) == 2

        special_excess = [
            sp.factor(
                sum(
                    vector[index] * monomial
                    for index, monomial in enumerate(monomials)
                )
            )
            for vector in special_excess_vectors
        ]
        direct_alpha, direct_beta = sp.symbols(
            "direct_alpha direct_beta"
        )
        combined_excess = (
            direct_alpha * special_excess[0]
            + direct_beta * special_excess[1]
        )
        special_left_kernel = special_jacobian.T.nullspace()
        special_second_rhs = output_vector(
            [
                (
                    -comb(order, 2)
                    * contraction(
                        combined_excess**2
                        * special_base ** (order - 2)
                    )
                    if order >= 2
                    else sp.Integer(0)
                )
                for order in range(1, 5)
            ]
        )
        special_second_obstructions = [
            sp.factor(
                (functional.T * special_second_rhs)[0]
            )
            for functional in special_left_kernel
        ]
        nonzero_second = [
            obstruction
            for obstruction in special_second_obstructions
            if obstruction != 0
        ]
        special_second_basis = sp.groebner(
            nonzero_second,
            direct_alpha,
            direct_beta,
            order="grevlex",
        )
        assert [
            polynomial.as_expr()
            for polynomial in special_second_basis.polys
        ] == [direct_beta**2]

        surviving_excess = special_excess[0]
        special_second_nonlinear = output_vector(
            [
                (
                    comb(order, 2)
                    * contraction(
                        surviving_excess**2
                        * special_base ** (order - 2)
                    )
                    if order >= 2
                    else sp.Integer(0)
                )
                for order in range(1, 5)
            ]
        )
        special_second_vector, special_second_free = (
            special_jacobian.gauss_jordan_solve(
                -special_second_nonlinear
            )
        )
        assert len(special_second_free) == 6
        special_second_correction = sum(
            coefficient * monomial
            for coefficient, monomial in zip(
                special_second_vector,
                monomials,
            )
        )
        special_third_nonlinear = output_vector(
            [
                (
                    order
                    * (order - 1)
                    * contraction(
                        surviving_excess
                        * special_second_correction
                        * special_base ** (order - 2)
                    )
                    if order >= 2
                    else sp.Integer(0)
                )
                + (
                    comb(order, 3)
                    * contraction(
                        surviving_excess**3
                        * special_base ** (order - 3)
                    )
                    if order >= 3
                    else sp.Integer(0)
                )
                for order in range(1, 5)
            ]
        )
        assert all(
            (functional.T * special_third_nonlinear)[0] == 0
            for functional in special_left_kernel
        )
        special_third_vector, special_third_free = (
            special_jacobian.gauss_jordan_solve(
                -special_third_nonlinear
            )
        )
        assert len(special_third_free) == 6
        special_third_correction = sum(
            coefficient * monomial
            for coefficient, monomial in zip(
                special_third_vector,
                monomials,
            )
        )
        special_fourth_nonlinear = output_vector(
            [
                (
                    order
                    * (order - 1)
                    * contraction(
                        surviving_excess
                        * special_third_correction
                        * special_base ** (order - 2)
                    )
                    if order >= 2
                    else sp.Integer(0)
                )
                + (
                    comb(order, 2)
                    * contraction(
                        special_second_correction**2
                        * special_base ** (order - 2)
                    )
                    if order >= 2
                    else sp.Integer(0)
                )
                + (
                    3
                    * comb(order, 3)
                    * contraction(
                        surviving_excess**2
                        * special_second_correction
                        * special_base ** (order - 3)
                    )
                    if order >= 3
                    else sp.Integer(0)
                )
                + (
                    comb(order, 4)
                    * contraction(
                        surviving_excess**4
                        * special_base ** (order - 4)
                    )
                    if order >= 4
                    else sp.Integer(0)
                )
                for order in range(1, 5)
            ]
        )
        special_fourth_obstructions = [
            sp.factor(
                (functional.T * special_fourth_nonlinear)[0]
            )
            for functional in special_left_kernel
        ]
        nonzero_fourth = [
            obstruction
            for obstruction in special_fourth_obstructions
            if obstruction != 0
        ]
        special_fourth_basis = sp.groebner(
            nonzero_fourth,
            *list(special_second_free),
            order="grevlex",
        )
        assert [
            polynomial.as_expr()
            for polynomial in special_fourth_basis.polys
        ] == [1]
        return {
            "label": label,
            "base_point": str(special_base),
            "surviving_excess": str(surviving_excess),
            "second_order_basis": ["direct_beta**2"],
            "fourth_order_obstructions": [
                str(obstruction) for obstruction in nonzero_fourth
            ],
            "fourth_order_basis": ["1"],
        }

    rational_special_points = [
        direct_fourth_gate(
            "u=-6",
            V * Z**2 * (-6 * V * Y + W * Z),
        ),
        direct_fourth_gate(
            "u=-3/4",
            V * Z**2 * (-sp.Rational(3, 4) * V * Y + W * Z),
        ),
        direct_fourth_gate(
            "u=0",
            W * V * Z**3,
        ),
        direct_fourth_gate(
            "u=infinity",
            V**2 * Y * Z**2,
        ),
        direct_fourth_gate(
            "a=b=c=1 representative",
            V * Z**2 * (V * Y + W * Z + V * Z),
        ),
    ]

    payload = {
        "claim": (
            "On the chart b*(a-b)!=0, the incidence family reduces to "
            "B(u)=V*Z^2*(u*V*Y+W*Z).  Away from the displayed finite "
            "exceptional set, both transverse directions of the first "
            "four full moment equations are obstructed by deformation "
            "order four."
        ),
        "claim_boundary": (
            "This is an exact generic local jet theorem, not a global "
            "classification of the mixed-positive branch and not an "
            "all-order SIC theorem."
        ),
        "incidence_reduction": {
            "raising_operator": "Z*d/dY-V*d/dW",
            "c_shift": "c -> c+shift*(a-b)",
            "normalized_chart": "b*(a-b)!=0; u=a/b; c=0; b=1",
        },
        "moment_jacobian": {
            "shape": [14, 12],
            "generic_rank": 6,
            "incidence_tangent_dimension": 4,
            "transverse_excess_dimension": 2,
        },
        "second_order_gate": {
            "conclusion": "the beta direction is obstructed",
            "coefficient_polynomials": [
                str(second_polynomial_1),
                str(second_polynomial_2),
            ],
            "coefficient_resultant": int(second_resultant),
            "chart_denominator": str(second_denominator),
        },
        "fourth_order_gate": {
            "second_correction_parameter_count": 6,
            "third_correction_parameter_count": 6,
            "compatibility_polynomials": [
                str(fourth_polynomial_1),
                str(fourth_polynomial_2),
            ],
            "chart_denominator": str(fourth_denominator),
            "resultant_factorization": str(fourth_resultant),
            "degree_6_exceptional_polynomial": str(
                exceptional_polynomial
            ),
            "exceptional_sextic_common_correction": {
                "formula": (
                    "tau=-B(u)/A(u) in Q[u]/(S), with A and B below"
                ),
                "A": str(exceptional_correction_coefficient),
                "B": str(exceptional_correction_constant),
                "resultant_A_S": int(
                    exceptional_correction_resultant
                ),
            },
        },
        "remaining_special_ratios": {
            "affine_normalization_exceptions": ["b=0", "a=b"],
            "calculation_chart_exceptions": [
                "u=-6",
                "4*u+3=0",
                str(chart_cubic),
            ],
            "fourth_resultant_exceptions": [
                "u=0",
                str(exceptional_polynomial),
            ],
            "note": (
                "Some chart exceptions may be removable basis artifacts; "
                "each must be recomputed directly before being treated "
                "as a geometric component."
            ),
        },
        "direct_rational_special_points": rational_special_points,
        "status": "exact characteristic-zero generic local theorem",
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    print(
        "PASS generic mixed-positive V_(2,3) local gate: one quotient "
        "parameter, second-order beta obstruction, and factored "
        "fourth-order resultant"
    )


if __name__ == "__main__":
    main()
