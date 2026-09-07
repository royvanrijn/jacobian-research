#!/usr/bin/env python3
"""Specialize the F2 lower-Laurent linear maps at the exact carrier points.

The carrier Wronskian leaves one rational squarefree cofactor and two
conjugate double-root cofactors.  The existing endpoint/Schur presentation
belongs to the earliest descent-eight component, whose support condition is
that the cofactor have a movable double root.  This checker therefore does
two distinct things.

* It proves that the squarefree carrier is not a point of that component and
  must be routed to a later-defect compiler rather than substituted into a
  nonexistent ``w0`` parameter.
* Over ``K=QQ[rho]/(rho^2-3*rho+1)`` it specializes every exposed linear map
  on the double-root component: the zero-row triangular maps through descent
  35 and the first lower row at descent 37, the complete target cokernel, and
  the complete layer-zero Hermite cokernel.

The local descent-eight ratio satisfies ``27*y^2-9*y+1=0``.  It is linearly
disjoint from ``K``, so the remaining nonlinear forcing lives over a quartic
compositum with four embeddings.  This checker fixes the exact kernel and
cokernel bases into which that forcing must be evaluated.  It does not yet
form the nonlinear forcing circuits on layers 28 through 3 and obtains
neither a unit ideal nor a Keller map.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from math import factorial
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.matrices import DomainMatrix

from classify_f2_75_125_layers import band_factor_data, make_band
from reduce_f2_75_125_endpoint_system import normalized_layer_support_audit


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/jc2_f2_75_125_carrier_specializations.json"
)

RHO_SYMBOL = sp.Symbol("rho")
RHO_POLYNOMIAL = sp.Poly(RHO_SYMBOL**2 - 3 * RHO_SYMBOL + 1, RHO_SYMBOL)
RHO_FIELD = QQ.alg_field_from_poly(RHO_POLYNOMIAL, alias="rho")
RHO = RHO_FIELD.ext.as_expr()


def field_element(expression: sp.Expr | int) -> object:
    """Convert an expression in ``rho`` to the exact quadratic field."""

    expression = sp.cancel(sp.sympify(expression))
    numerator, denominator = sp.fraction(expression)

    def evaluate(polynomial_expression: sp.Expr) -> object:
        polynomial = sp.Poly(polynomial_expression, RHO_SYMBOL, domain=sp.QQ)
        value = RHO_FIELD.zero
        for coefficient in polynomial.all_coeffs():
            value = value * RHO_FIELD.unit + RHO_FIELD.convert(coefficient)
        return value

    return evaluate(numerator) / evaluate(denominator)


def field_matrix(rows: list[list[sp.Expr | int]]) -> DomainMatrix:
    """Build a dense exact matrix over ``QQ(rho)``."""

    if not rows:
        return DomainMatrix([], (0, 0), RHO_FIELD)
    column_count = len(rows[0])
    if any(len(row) != column_count for row in rows):
        raise ValueError("matrix rows have unequal lengths")
    converted = [
        [field_element(entry) for entry in row]
        for row in rows
    ]
    return DomainMatrix(converted, (len(rows), column_count), RHO_FIELD)


def rational_record(value: object) -> list[int]:
    return [int(value.numerator), int(value.denominator)]


def field_record(value: object) -> list[list[int]]:
    """Serialize ``a*rho+b`` as ``[b,a]`` rational pairs."""

    coefficients = list(value.to_list())
    coefficients = [RHO_FIELD.dom.zero] * (2 - len(coefficients)) + coefficients
    rho_coefficient, constant = coefficients
    return [rational_record(constant), rational_record(rho_coefficient)]


def matrix_record(matrix: DomainMatrix) -> list[list[list[list[int]]]]:
    return [
        [field_record(entry) for entry in row]
        for row in matrix.to_list()
    ]


def matrix_digest(matrix: DomainMatrix) -> str:
    encoded = json.dumps(
        {"shape": list(matrix.shape), "entries": matrix_record(matrix)},
        separators=(",", ":"),
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def reduced_kernel(matrix: DomainMatrix) -> DomainMatrix:
    kernel = matrix.nullspace()
    if kernel.shape[0] == 0:
        return kernel
    return kernel.rref()[0]


def reduced_left_cokernel(matrix: DomainMatrix) -> DomainMatrix:
    cokernel = matrix.transpose().nullspace()
    if cokernel.shape[0] == 0:
        return cokernel
    return cokernel.rref()[0]


def polynomial_coefficient_matrix(
    columns: list[sp.Expr], variable: sp.Symbol, maximum_degree: int
) -> DomainMatrix:
    rows: list[list[sp.Expr]] = [
        [sp.Integer(0) for _ in columns] for _ in range(maximum_degree + 1)
    ]
    for column, expression in enumerate(columns):
        polynomial = sp.Poly(sp.expand(expression), variable)
        if polynomial.degree() > maximum_degree:
            raise AssertionError("a polynomial column exceeded its ambient degree")
        for degree in range(maximum_degree + 1):
            rows[degree][column] = polynomial.coeff_monomial(variable**degree)
    return field_matrix(rows)


def evaluation_row(
    maximum_degree: int,
    point: sp.Expr,
    derivative_order: int,
) -> list[sp.Expr]:
    row: list[sp.Expr] = []
    for degree in range(maximum_degree + 1):
        if degree < derivative_order:
            row.append(sp.Integer(0))
        else:
            coefficient = factorial(degree) // factorial(degree - derivative_order)
            row.append(coefficient * point ** (degree - derivative_order))
    return row


def carrier_routing_audit() -> dict[str, object]:
    """Separate the squarefree carrier from the descent-eight component."""

    w = sp.Symbol("w")
    squarefree_R = (w**2 - 3 * w + 3) / 25
    squarefree_discriminant = sp.factor(sp.discriminant(squarefree_R, w))
    if squarefree_discriminant != -sp.Rational(3, 625):
        raise AssertionError("the squarefree carrier discriminant changed")

    double_R = (w - RHO_SYMBOL) ** 2 / (25 * (1 - RHO_SYMBOL) ** 2)
    reduced_double_discriminant = sp.rem(
        sp.Poly(sp.together(sp.discriminant(double_R, w)).as_numer_denom()[0], RHO_SYMBOL),
        RHO_POLYNOMIAL,
    ).as_expr()
    if reduced_double_discriminant != 0:
        raise AssertionError("the double carrier left the discriminant-zero stratum")

    rho_inverse = 3 - RHO_SYMBOL
    rho_minus_one_inverse = RHO_SYMBOL - 2
    for identity in (
        RHO_SYMBOL * rho_inverse - 1,
        (RHO_SYMBOL - 1) * rho_minus_one_inverse - 1,
    ):
        if sp.rem(sp.Poly(identity, RHO_SYMBOL), RHO_POLYNOMIAL).as_expr() != 0:
            raise AssertionError("a localized double-carrier unit changed")

    return {
        "squarefree_carrier": {
            "R": "(w^2-3*w+3)/25",
            "discriminant": "-3/625",
            "descent_eight_component_membership": False,
            "routing": (
                "branch before the movable-double-root specialization and "
                "compile the later first-defect strata"
            ),
        },
        "double_carrier": {
            "R": "(w-rho)^2/(25*(1-rho)^2)",
            "field": "K=QQ[rho]/(rho^2-3*rho+1)",
            "carrier_points": 2,
            "descent_eight_component_membership": True,
            "localized_units": {
                "rho_inverse": "3-rho",
                "rho_minus_one_inverse": "rho-2",
            },
        },
        "claim_boundary": (
            "the existing endpoint reduction is a component presentation, "
            "not a family into which the squarefree carrier can be substituted"
        ),
    }


def compositum_audit() -> dict[str, object]:
    """Compile the exact coefficient field needed by nonlinear forcing."""

    rho, y, theta = sp.symbols("rho y theta")
    rho_polynomial = rho**2 - 3 * rho + 1
    defect_polynomial = 27 * y**2 - 9 * y + 1
    quartic = sp.Poly(
        sp.resultant(
            rho_polynomial,
            defect_polynomial.subs(y, theta - rho),
            rho,
        ),
        theta,
        domain=sp.QQ,
    )
    expected = sp.Poly(
        729 * theta**4
        - 4860 * theta**3
        + 10341 * theta**2
        - 7470 * theta
        + 1756,
        theta,
        domain=sp.QQ,
    )
    if quartic != expected or not quartic.is_irreducible:
        raise AssertionError("the carrier/defect compositum changed")

    return {
        "carrier_polynomial": "rho^2-3*rho+1",
        "carrier_discriminant_square_class": 5,
        "descent_eight_polynomial": "27*y^2-9*y+1",
        "descent_eight_discriminant_square_class": -3,
        "linearly_disjoint": True,
        "compositum_degree": 4,
        "primitive_element": "theta=rho+y",
        "primitive_element_polynomial": str(quartic.as_expr()),
        "geometric_branch_count": 4,
        "linear_map_field": "QQ(rho); y enters only the nonlinear forcing/incidence rows",
    }


def normalized_q_operator(descent: int) -> tuple[DomainMatrix, dict[str, int]]:
    """Return the fixed-top q-only operator at one Laurent descent."""

    q_layer = 25 - descent
    q_band = make_band("Q", 125, 5, q_layer)
    q_u, q_vanishing, q_degree = band_factor_data(q_band)
    rows = [
        [sp.Integer(0) for _ in range(q_degree + 1)]
        for _ in range(q_degree + 3)
    ]
    for column in range(q_degree + 1):
        rows[column][column] = RHO_SYMBOL * (5 * column + q_u)
        rows[column + 1][column] = (
            (RHO_SYMBOL + 1) * (2 * q_layer - 5 * column - q_u)
            - 5 * q_vanishing * RHO_SYMBOL
        )
        rows[column + 2][column] = (
            5 * column + q_u + 5 * q_vanishing - 4 * q_layer
        )
    return field_matrix(rows), {
        "q_layer": q_layer,
        "q_u": q_u,
        "q_vanishing": q_vanishing,
        "q_degree": q_degree,
    }


def full_q_operator(
    descent: int,
    minimum_degree: int,
    maximum_degree: int,
    *,
    endpoint_domain: bool,
) -> DomainMatrix:
    """Build the actual Laurent-row map before common-factor division."""

    w = sp.Symbol("w")
    p_layer = 15
    q_layer = 25 - descent
    p_band = make_band("P", 75, 3, p_layer)
    q_band = make_band("Q", 125, 5, q_layer)
    p_u, _, _ = band_factor_data(p_band)
    q_u, q_vanishing, q_degree = band_factor_data(q_band)
    R = (w - RHO_SYMBOL) ** 2 / (25 * (1 - RHO_SYMBOL) ** 2)
    p_factor = (w - 1) ** 6 * R**3
    shift = (p_u + q_u - 1) // 5
    orders = list(range(q_degree + 1))
    if endpoint_domain and q_layer in (1, 13):
        orders.remove(0)

    columns: list[sp.Expr] = []
    for order in orders:
        q_factor = (w - 1) ** (q_vanishing + order)
        columns.append(
            sp.cancel(
                w**shift
                * (
                    p_layer
                    * p_factor
                    * (q_u * q_factor + 5 * w * sp.diff(q_factor, w))
                    - q_layer
                    * (p_u * p_factor + 5 * w * sp.diff(p_factor, w))
                    * q_factor
                )
            )
        )

    rows: list[list[sp.Expr]] = [
        [sp.Integer(0) for _ in columns]
        for _ in range(maximum_degree - minimum_degree + 1)
    ]
    for column, expression in enumerate(columns):
        polynomial = sp.Poly(sp.cancel(expression), w)
        for degree in range(minimum_degree, maximum_degree + 1):
            rows[degree - minimum_degree][column] = polynomial.coeff_monomial(
                w**degree
            )
    return field_matrix(rows)


def descent_linear_map_audit() -> dict[str, object]:
    """Specialize all zero-row Schur maps exposed by the current compiler."""

    descents = [*range(1, 36), 37]
    p_bands = {
        layer: make_band("P", 75, 3, layer) for layer in range(-75, 16)
    }
    q_bands = {
        layer: make_band("Q", 125, 5, layer) for layer in range(-125, 26)
    }
    support = normalized_layer_support_audit(p_bands, q_bands)
    support_by_layer = {
        int(row["layer"]): row for row in support["rows"]
    }

    records: list[dict[str, object]] = []
    for descent in descents:
        operator, band = normalized_q_operator(descent)
        kernel = reduced_kernel(operator)
        cokernel = reduced_left_cokernel(operator)
        expected_kernel = int(descent in (5, 10, 15, 20, 25))
        if kernel.shape[0] != expected_kernel:
            raise AssertionError("a specialized q-only kernel changed")
        if cokernel.shape[0] != 2 + expected_kernel:
            raise AssertionError("a specialized q-only cokernel changed")

        p_band = p_bands[15 - descent]
        layer = 40 - descent
        layer_row = support_by_layer[layer]
        minimum_degree, maximum_degree = layer_row["w_degree_interval"]
        actual_operator = full_q_operator(
            descent,
            int(minimum_degree),
            int(maximum_degree),
            endpoint_domain=True,
        )
        actual_cokernel = reduced_left_cokernel(actual_operator)
        endpoint_fixed_coordinate = int(25 - descent in (1, 13))
        expected_actual_rank = operator.rank() - endpoint_fixed_coordinate
        if actual_operator.rank() != expected_actual_rank:
            raise AssertionError("a full specialized Laurent-row rank changed")
        records.append(
            {
                "descent": descent,
                "Laurent_layer": layer,
                "P_band": 15 - descent,
                "P_follow_dimension": p_band.dimension,
                "Q_band": int(band["q_layer"]),
                "quotient_matrix_shape": list(operator.shape),
                "quotient_rank": operator.rank(),
                "Q_kernel_dimension": kernel.shape[0],
                "full_tangent_kernel_dimension": (
                    p_band.dimension + kernel.shape[0]
                ),
                "quotient_forcing_cokernel_dimension": cokernel.shape[0],
                "quotient_matrix_digest_sha256": matrix_digest(operator),
                "kernel_basis_digest_sha256": matrix_digest(kernel),
                "quotient_cokernel_basis_digest_sha256": matrix_digest(cokernel),
                "endpoint_fixed_Q_coordinate": bool(endpoint_fixed_coordinate),
                "full_Laurent_matrix_shape": list(actual_operator.shape),
                "full_Laurent_rank": actual_operator.rank(),
                "full_Laurent_forcing_cokernel_dimension": (
                    actual_cokernel.shape[0]
                ),
                "full_Laurent_matrix_digest_sha256": matrix_digest(
                    actual_operator
                ),
                "full_Laurent_cokernel_basis_digest_sha256": matrix_digest(
                    actual_cokernel
                ),
            }
        )

    upper = [row for row in records if int(row["descent"]) <= 11]
    coupled = [row for row in records if int(row["descent"]) >= 12]
    upper_rank = sum(int(row["full_Laurent_rank"]) for row in upper)
    upper_kernel = sum(int(row["Q_kernel_dimension"]) for row in upper)
    if (upper_rank, upper_kernel) != (134, 2):
        raise AssertionError("the specialized upper block changed")
    coupled_quotient_cokernel = sum(
        int(row["quotient_forcing_cokernel_dimension"]) for row in coupled
    )
    if coupled_quotient_cokernel != 53:
        raise AssertionError("the coupled raw Fitting-coordinate count changed")
    upper_full_cokernel = sum(
        int(row["full_Laurent_forcing_cokernel_dimension"]) for row in upper
    )
    coupled_full_cokernel = sum(
        int(row["full_Laurent_forcing_cokernel_dimension"]) for row in coupled
    )
    if (upper_full_cokernel, coupled_full_cokernel) != (219, 347):
        raise AssertionError("the full Laurent-row cokernel ledger changed")

    return {
        "field": "QQ(rho)",
        "operator": (
            "the normalized form of T_delta(q)="
            "5*C0*q'-(25-delta)*C0'*q"
        ),
        "arbitrary_P_follow": "q_follow=-3*C0^2*p",
        "rows": records,
        "endpoint_disjoint_descents_1_through_11": {
            "rank": upper_rank,
            "Q_kernel_dimension": upper_kernel,
            "full_Laurent_forcing_cokernel_dimension": upper_full_cokernel,
            "agreement_with_parametric_reducer": True,
        },
        "coupled_descents": [[12, 35], [37, 37]],
        "coupled_quotient_cokernel_after_common_factor": (
            coupled_quotient_cokernel
        ),
        "coupled_full_Laurent_forcing_cokernel_dimension": (
            coupled_full_cokernel
        ),
        "warning": (
            "the 53 quotient coordinates do not include divisibility/local-jet "
            "conditions.  Nonlinear forcing must be projected into the full "
            "347-dimensional coupled Laurent-row cokernel"
        ),
    }


def target_cokernel_audit() -> dict[str, object]:
    """Build the complete degree-33 target image and its cokernel."""

    w = sp.Symbol("w")
    E = (w - 1) * (w - RHO_SYMBOL)
    local_modulus = sp.expand(w**2 * (w - 1) ** 5 * (w - RHO_SYMBOL) ** 5)
    quotient_columns: list[sp.Expr] = []
    for degree in range(20):
        source = w**degree
        quotient_columns.append(
            sp.expand(
                5 * w * E * sp.diff(source, w)
                + (11 * E + 22 * w * sp.diff(E, w)) * source
            )
        )
    quotient_operator = polynomial_coefficient_matrix(quotient_columns, w, 21)
    if quotient_operator.shape != (22, 20) or quotient_operator.rank() != 20:
        raise AssertionError("the specialized target quotient operator changed")
    quotient_cokernel = reduced_left_cokernel(quotient_operator)
    if quotient_cokernel.shape != (2, 22):
        raise AssertionError("the target quotient cokernel changed")

    full_columns = [local_modulus * column for column in quotient_columns]
    full_image = polynomial_coefficient_matrix(full_columns, w, 33)
    full_cokernel = reduced_left_cokernel(full_image)
    if full_image.rank() != 20 or full_cokernel.shape != (14, 34):
        raise AssertionError("the complete target cokernel changed")

    local_rows = [
        evaluation_row(33, sp.Integer(0), order) for order in range(2)
    ] + [
        evaluation_row(33, sp.Integer(1), order) for order in range(5)
    ] + [
        evaluation_row(33, RHO_SYMBOL, order) for order in range(5)
    ]
    local_map = field_matrix(local_rows)
    if local_map.rank() != 12:
        raise AssertionError("the specialized target local-jet rank changed")
    if not local_map.matmul(full_image).is_zero_matrix:
        raise AssertionError("a target image column violated its local jets")
    if DomainMatrix.vstack(local_map, full_cokernel).rank() != 14:
        raise AssertionError("the local target jets left the complete cokernel")

    rho_local_map = field_matrix(
        [evaluation_row(33, RHO_SYMBOL, order) for order in range(5)]
    )
    if rho_local_map.rank() != 5:
        raise AssertionError("the movable target jet block changed")

    return {
        "ambient": "G_red in K[w]_(<=33)",
        "local_modulus": "D=w^2*(w-1)^5*(w-rho)^5",
        "quotient_operator": {
            "formula": (
                "N(S)=5*w*E*S'+(11*E+22*w*E')*S, "
                "E=(w-1)*(w-rho)"
            ),
            "matrix_shape": list(quotient_operator.shape),
            "rank": quotient_operator.rank(),
            "cokernel_dimension": quotient_cokernel.shape[0],
            "matrix_digest_sha256": matrix_digest(quotient_operator),
            "cokernel_basis_digest_sha256": matrix_digest(quotient_cokernel),
        },
        "complete_image": {
            "formula": "D*N(K[w]_(<=19)) in K[w]_(<=33)",
            "matrix_shape": list(full_image.shape),
            "rank": full_image.rank(),
            "cokernel_dimension": full_cokernel.shape[0],
            "matrix_digest_sha256": matrix_digest(full_image),
            "cokernel_basis_digest_sha256": matrix_digest(full_cokernel),
        },
        "local_jet_block": {
            "profile": "orders 0..1 at 0, 0..4 at 1, and 0..4 at rho",
            "rank": local_map.rank(),
            "contained_in_complete_cokernel": True,
        },
        "after_prior_control_and_fixed_endpoint_elimination": {
            "movable_rho_jets": rho_local_map.rank(),
            "quotient_fitting_coordinates": quotient_cokernel.shape[0],
            "remaining_target_coordinates": (
                rho_local_map.rank() + quotient_cokernel.shape[0]
            ),
        },
    }


def layer_zero_cokernel_audit() -> dict[str, object]:
    """Build the degree-33 layer-zero image and its Hermite cokernel."""

    w = sp.Symbol("w")
    modulus = sp.expand(w**3 * (w - 1) ** 6 * (w - RHO_SYMBOL) ** 6)
    image_columns = [sp.Integer(1)] + [
        modulus * w**degree for degree in range(19)
    ]
    image = polynomial_coefficient_matrix(image_columns, w, 33)
    cokernel = reduced_left_cokernel(image)
    if image.rank() != 20 or cokernel.shape != (14, 34):
        raise AssertionError("the specialized layer-zero image changed")

    h_at_zero = evaluation_row(33, sp.Integer(0), 0)
    hermite_rows = [
        evaluation_row(33, sp.Integer(0), order) for order in (1, 2)
    ]
    row = evaluation_row(33, sp.Integer(1), 0)
    hermite_rows.append([left - right for left, right in zip(row, h_at_zero)])
    hermite_rows.extend(
        evaluation_row(33, sp.Integer(1), order) for order in range(1, 6)
    )
    row = evaluation_row(33, RHO_SYMBOL, 0)
    hermite_rows.append([left - right for left, right in zip(row, h_at_zero)])
    hermite_rows.extend(
        evaluation_row(33, RHO_SYMBOL, order) for order in range(1, 6)
    )
    hermite_map = field_matrix(hermite_rows)
    if hermite_map.rank() != 14:
        raise AssertionError("the complete Hermite rank changed")
    if not hermite_map.matmul(image).is_zero_matrix:
        raise AssertionError("a layer-zero image column violated Hermite membership")
    if DomainMatrix.vstack(hermite_map, cokernel).rank() != 14:
        raise AssertionError("Hermite rows and the exact cokernel disagree")

    residual_rows = [
        [
            left - right
            for left, right in zip(
                evaluation_row(33, RHO_SYMBOL, 0), h_at_zero
            )
        ]
    ] + [
        evaluation_row(33, RHO_SYMBOL, order) for order in range(1, 6)
    ]
    residual_map = field_matrix(residual_rows)
    if residual_map.rank() != 6:
        raise AssertionError("the residual movable Hermite rank changed")

    return {
        "ambient": "H_red in K[w]_(<=33)",
        "image": "K*1 + w^3*(w-1)^6*(w-rho)^6*K[w]_(<=18)",
        "image_matrix_shape": list(image.shape),
        "image_rank": image.rank(),
        "cokernel_dimension": cokernel.shape[0],
        "image_matrix_digest_sha256": matrix_digest(image),
        "cokernel_basis_digest_sha256": matrix_digest(cokernel),
        "complete_Hermite_map": {
            "profile": (
                "H'(0),H''(0); H(1)-H(0),H^(1..5)(1); "
                "H(rho)-H(0),H^(1..5)(rho)"
            ),
            "matrix_shape": list(hermite_map.shape),
            "rank": hermite_map.rank(),
            "matrix_digest_sha256": matrix_digest(hermite_map),
            "same_row_space_as_cokernel": True,
        },
        "after_prior_control_and_fixed_endpoint_elimination": {
            "matrix_shape": list(residual_map.shape),
            "rank": residual_map.rank(),
            "matrix_digest_sha256": matrix_digest(residual_map),
            "remaining_layer_zero_coordinates": residual_map.rank(),
        },
    }


def build_payload() -> dict[str, object]:
    routing = carrier_routing_audit()
    descents = descent_linear_map_audit()
    target = target_cokernel_audit()
    layer_zero = layer_zero_cokernel_audit()
    remaining_coordinates = (
        int(
            target["after_prior_control_and_fixed_endpoint_elimination"][
                "remaining_target_coordinates"
            ]
        )
        + int(
            layer_zero["after_prior_control_and_fixed_endpoint_elimination"][
                "remaining_layer_zero_coordinates"
            ]
        )
    )
    if remaining_coordinates != 13:
        raise AssertionError("the specialized residual coordinate count changed")

    return {
        "schema": "plane-jc.f2-75-125-carrier-specializations.v1",
        "status": "exact-number-field-linear-maps;nonlinear-forcing-open",
        "carrier_routing": routing,
        "double_branch_coefficient_field": compositum_audit(),
        "specialized_zero_row_maps": descents,
        "specialized_target_cokernel": target,
        "specialized_layer_zero_cokernel": layer_zero,
        "combined_residual_coordinate_count": remaining_coordinates,
        "decision": {
            "unit_ideal_obtained": False,
            "counterexample_obtained": False,
            "double_branch_progress": (
                "all currently exposed kernel/cokernel maps are now fixed "
                "over QQ(rho), including the 53-dimensional factor-quotient "
                "target, the full 347-dimensional coupled Laurent cokernel, "
                "and the final 7+6 coordinates"
            ),
            "squarefree_branch_progress": (
                "proved that the rational carrier is outside the existing "
                "descent-eight component; it must be compiled from the later "
                "first-defect ledger rather than falsely specialized at w0"
            ),
            "next_exact_operation": (
                "use compile_f2_75_125_nonlinear_forcing.py to build the "
                "294+53 coupled equations and final 13 functionals, then "
                "test the localized circuit ideal by good-reduction methods"
            ),
        },
        "reproduction_command": (
            ".venv/bin/python plane-jc/cas/"
            "verify_f2_75_125_carrier_specializations.py"
        ),
        "software": {"sympy": sp.__version__},
    }


def artifact_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=ARTIFACT)
    parser.add_argument("--refresh", action="store_true")
    args = parser.parse_args()

    payload = build_payload()
    artifact = args.artifact.resolve()
    if args.refresh:
        artifact.parent.mkdir(parents=True, exist_ok=True)
        artifact.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        try:
            display = artifact.relative_to(ROOT)
        except ValueError:
            display = artifact
        print(f"WROTE {display}")
    else:
        expected = json.loads(artifact.read_text())
        current_claim = {key: value for key, value in payload.items() if key != "software"}
        pinned_claim = {key: value for key, value in expected.items() if key != "software"}
        if current_claim != pinned_claim:
            raise AssertionError(
                "the pinned carrier-specialization artifact is stale; inspect "
                "the change before using --refresh"
            )

    print("F2_CARRIER_SPECIALIZATION_ROUTING_PASS")
    print("F2_DOUBLE_BRANCH_QUADRATIC_FIELD_PASS")
    print("F2_SPECIALIZED_ZERO_ROW_KERNEL_COKERNEL_PASS")
    print("F2_SPECIALIZED_TARGET_COKERNEL_PASS")
    print("F2_SPECIALIZED_LAYER_ZERO_COKERNEL_PASS")
    print("F2_SQUAREFREE_LATER_DEFECT_ROUTE_OPEN")
    print("F2_DOUBLE_NONLINEAR_FORCING_HANDED_OFF_TO_PF2NF1")
    print(f"F2_CARRIER_SPECIALIZATION_ARTIFACT_SHA256={artifact_sha256(artifact)}")


if __name__ == "__main__":
    main()
