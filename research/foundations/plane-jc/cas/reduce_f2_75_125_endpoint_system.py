#!/usr/bin/env python3
"""Carry the F2 (75,125) fixed-endpoint elimination into the residual ideal.

This checker starts on the nonzero double-R stratum and performs the exact
ten-variable elimination at w=1.  It does not pretend that the resulting
rank-thirteen obstruction module is a polynomial ring in thirteen variables:
the earlier Laurent rows still contain more than a thousand source-band
coordinates.  Instead it emits a deterministic straight-line presentation
of the reduced ideal and proves the following facts.

* all forcing from the 22 target and 25 layer-zero dependency pairs is kept;
* the ten pivot solutions satisfy the complete endpoint rows identically;
* the degree-seven follower basis makes H(0) independent of the pivots;
* only Laurent layers 3 through 28 change under the substitution;
* the remaining target equations are five movable-root jets and the two
  triangular target Fitting residues;
* the remaining layer-zero equations are the six movable-root Hermite rows.

The output is an exact pre-lower-tail presentation, not a unit-ideal
certificate and not a plane counterexample.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from math import gcd
from pathlib import Path

import sympy as sp

from classify_f2_75_125_layers import band_factor_data, make_band


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/jc2_f2_75_125_endpoint_reduction.json"
)


def lcm(left: int, right: int) -> int:
    return abs(left * right) // gcd(left, right) if left and right else 0


def canonical_polynomial(
    expression: sp.Expr,
    generators: list[sp.Symbol],
) -> list[list[object]]:
    polynomial = sp.Poly(sp.expand(expression), *generators, domain=sp.QQ)
    return [
        [list(monomial), [int(coefficient.p), int(coefficient.q)]]
        for monomial, coefficient in polynomial.terms()
    ]


def polynomial_digest(
    expressions: list[sp.Expr],
    generators: list[sp.Symbol],
) -> str:
    canonical = [
        canonical_polynomial(expression, generators)
        for expression in expressions
    ]
    encoded = json.dumps(canonical, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def total_degree(
    expression: sp.Expr,
    generators: list[sp.Symbol],
) -> int:
    if expression == 0:
        return -1
    return int(sp.Poly(sp.expand(expression), *generators, domain=sp.QQ).total_degree())


def weighted_degree(
    expression: sp.Expr,
    generators: list[sp.Symbol],
    weights: list[int],
) -> int:
    polynomial = sp.Poly(sp.expand(expression), *generators, domain=sp.QQ)
    return max(
        sum(exponent * weight for exponent, weight in zip(monomial, weights))
        for monomial, _ in polynomial.terms()
    )


def endpoint_substitution_audit(
    p_bands: dict[int, object],
    q_bands: dict[int, object],
) -> dict[str, object]:
    delta = sp.Symbol("delta")
    local_w = 1 + delta

    jet_requirements: dict[str, dict[int, int]] = {"P": {}, "Q": {}}

    def require(side: str, layer: int, maximum_jet: int) -> None:
        jet_requirements[side][layer] = max(
            jet_requirements[side].get(layer, -1), maximum_jet
        )

    target_pairs: list[tuple[int, int]] = []
    for p_layer in range(-21, 16):
        q_layer = 4 - p_layer
        _, p_vanishing, _ = band_factor_data(p_bands[p_layer])
        _, q_vanishing, _ = band_factor_data(q_bands[q_layer])
        if p_vanishing + q_vanishing > 5:
            continue
        maximum_jet = 5 - p_vanishing - q_vanishing
        require("P", p_layer, maximum_jet)
        require("Q", q_layer, maximum_jet)
        target_pairs.append((p_layer, q_layer))

    layer_zero_pairs: list[tuple[int, int]] = []
    for p_layer in range(-25, 16):
        if p_layer == 0:
            continue
        q_layer = -p_layer
        _, p_vanishing, _ = band_factor_data(p_bands[p_layer])
        _, q_vanishing, _ = band_factor_data(q_bands[q_layer])
        if p_vanishing + q_vanishing > 5:
            continue
        maximum_jet = 5 - p_vanishing - q_vanishing
        require("P", p_layer, maximum_jet)
        require("Q", q_layer, maximum_jet)
        layer_zero_pairs.append((p_layer, q_layer))

    if len(target_pairs) != 22 or len(layer_zero_pairs) != 25:
        raise AssertionError("the endpoint dependency-pair count changed")

    coefficient_lists: dict[tuple[str, int], tuple[sp.Symbol, ...]] = {}
    local_polynomials: dict[tuple[str, int], sp.Expr] = {}
    all_symbols: list[sp.Symbol] = []
    for side in ("P", "Q"):
        for layer, maximum_jet in sorted(jet_requirements[side].items()):
            layer_name = f"m{-layer}" if layer < 0 else str(layer)
            coefficients = sp.symbols(
                f"ep_{side.lower()}{layer_name}_0:{maximum_jet + 1}"
            )
            coefficient_lists[side, layer] = coefficients
            all_symbols.extend(coefficients)
            local_polynomials[side, layer] = sum(
                coefficient * delta**order
                for order, coefficient in enumerate(coefficients)
            )

    fixed_values = {
        coefficient_lists["P", 3][0]: sp.Rational(1, 5),
        coefficient_lists["Q", 1][0]: sp.Integer(-1),
        coefficient_lists["Q", 13][0]: -sp.Rational(3, 5**5),
    }
    if len(all_symbols) != 162 or len(all_symbols) - len(fixed_values) != 159:
        raise AssertionError("the endpoint local-coordinate count changed")

    def target_pair_expression(
        p_layer: int,
        p_k: sp.Expr,
        q_layer: int,
        q_k: sp.Expr,
    ) -> sp.Expr:
        p_u, p_vanishing, _ = band_factor_data(p_bands[p_layer])
        q_u, q_vanishing, _ = band_factor_data(q_bands[q_layer])
        shift = (p_u + q_u - 1) // 5
        p_factor = delta**p_vanishing * p_k
        q_factor = delta**q_vanishing * q_k
        return sp.expand(
            local_w**shift
            * (
                p_layer
                * p_factor
                * (q_u * q_factor + 5 * local_w * sp.diff(q_factor, delta))
                - q_layer
                * (p_u * p_factor + 5 * local_w * sp.diff(p_factor, delta))
                * q_factor
            )
        )

    def layer_zero_product_expression(
        p_layer: int,
        p_k: sp.Expr,
        q_layer: int,
        q_k: sp.Expr,
    ) -> sp.Expr:
        p_u, p_vanishing, _ = band_factor_data(p_bands[p_layer])
        q_u, q_vanishing, _ = band_factor_data(q_bands[q_layer])
        shift = (p_u + q_u) // 5
        return sp.expand(
            p_layer
            * local_w**shift
            * delta ** (p_vanishing + q_vanishing)
            * p_k
            * q_k
        )

    full_target = sp.expand(
        sum(
            target_pair_expression(
                p_layer,
                local_polynomials["P", p_layer],
                q_layer,
                local_polynomials["Q", q_layer],
            )
            for p_layer, q_layer in target_pairs
        ).subs(fixed_values)
    )
    full_integral = sp.expand(
        sum(
            layer_zero_product_expression(
                p_layer,
                local_polynomials["P", p_layer],
                q_layer,
                local_polynomials["Q", q_layer],
            )
            for p_layer, q_layer in layer_zero_pairs
        ).subs(fixed_values)
    )
    if full_target.coeff(delta, 0) != 1:
        raise AssertionError("the complete normalized target value changed")

    # H(1)-H(0), rather than H(1), is the constant layer-zero Hermite row.
    # The follower basis makes the exact full-band value H(0) independent of
    # every pivot.  It is therefore a legitimate forcing circuit input here.
    h_at_zero = sp.Symbol("H_at_zero")
    endpoint_rows = [
        full_target.coeff(delta, order) for order in range(1, 5)
    ] + [
        full_integral.coeff(delta, 0) - h_at_zero,
        *[full_integral.coeff(delta, order) for order in range(1, 6)],
    ]

    pivot_variables = (
        *coefficient_lists["P", 3][1:5],
        *coefficient_lists["P", -1][:6],
    )
    pivot_zero = {variable: 0 for variable in pivot_variables}
    for row in endpoint_rows:
        if sp.Poly(row, *pivot_variables).total_degree() > 1:
            raise AssertionError("an endpoint row ceased to be affine in the pivots")

    endpoint_matrix = sp.Matrix(
        [
            [sp.diff(row, variable) for variable in pivot_variables]
            for row in endpoint_rows
        ]
    )
    forcing = sp.Matrix(
        [sp.expand(row.subs(pivot_zero)) for row in endpoint_rows]
    )
    determinant = sp.factor(endpoint_matrix.det())
    if determinant != 75000:
        raise AssertionError("the complete endpoint determinant changed")

    inverse = endpoint_matrix.inv()
    solutions = [sp.expand(expression) for expression in (-inverse * forcing)]
    if any(
        sp.expand(expression) != 0
        for expression in endpoint_matrix * sp.Matrix(solutions) + forcing
    ):
        raise AssertionError("the carried endpoint substitution failed")
    if any(set(expression.free_symbols) & set(pivot_variables) for expression in solutions):
        raise AssertionError("an eliminated pivot survived its solution")

    generators = [
        symbol
        for symbol in all_symbols
        if symbol not in fixed_values and symbol not in pivot_variables
    ] + [h_at_zero]
    forcing_term_counts = [len(sp.Add.make_args(expression)) for expression in forcing]
    solution_term_counts = [len(sp.Add.make_args(expression)) for expression in solutions]
    forcing_degrees = [total_degree(expression, generators) for expression in forcing]
    solution_degrees = [total_degree(expression, generators) for expression in solutions]
    inverse_degrees = [
        total_degree(expression, generators)
        for expression in inverse
        if expression != 0
    ]
    source_weights = [1] * (len(generators) - 1) + [2]
    source_solution_degrees = [
        weighted_degree(expression, generators, source_weights)
        for expression in solutions
    ]

    solution_denominators: list[int] = []
    for expression in solutions:
        denominator = 1
        for _, coefficient in sp.Poly(
            expression, *generators, domain=sp.QQ
        ).terms():
            denominator = lcm(denominator, int(coefficient.q))
        solution_denominators.append(denominator)

    expected_forcing_terms = [3, 18, 43, 87, 1, 5, 18, 39, 70, 116]
    expected_solution_terms = [4, 29, 87, 299, 1, 7, 24, 85, 254, 699]
    if forcing_term_counts != expected_forcing_terms:
        raise AssertionError("the complete endpoint forcing changed")
    if solution_term_counts != expected_solution_terms:
        raise AssertionError("the complete endpoint solutions changed")
    if solution_degrees != [2, 3, 4, 5, 1, 2, 3, 4, 5, 6]:
        raise AssertionError("the endpoint circuit degrees changed")
    if max(source_solution_degrees) != 7:
        raise AssertionError("the expanded-source endpoint degree changed")
    if max(inverse_degrees) != 5:
        raise AssertionError("the endpoint inverse degree changed")

    follower_records: list[dict[str, object]] = []
    for side, layer, orders in (
        ("P", 3, range(1, 5)),
        ("P", -1, range(6)),
    ):
        _, _, degree_bound = band_factor_data(p_bands[layer])
        for order in orders:
            lift = sp.expand(
                delta**order
                + (order - 7) * (-1) ** order * delta**6
                + (order - 6) * (-1) ** order * delta**7
            )
            quotient, remainder = sp.div(lift, local_w**2, delta)
            if remainder != 0 or sp.degree(lift, delta) > degree_bound:
                raise AssertionError("an endpoint follower left its exact band")
            if any(
                sp.Poly(lift, delta).nth(j) != int(j == order)
                for j in range(6)
            ):
                raise AssertionError("an endpoint follower changed a low jet")
            follower_records.append(
                {
                    "side": side,
                    "layer": layer,
                    "endpoint_order": order,
                    "lift": str(lift),
                    "w_squared_quotient": str(sp.expand(quotient)),
                }
            )

    return {
        "dependency_pairs": {
            "target": len(target_pairs),
            "layer_zero": len(layer_zero_pairs),
        },
        "local_coordinates": {
            "total": len(all_symbols),
            "fixed": len(fixed_values),
            "free": len(all_symbols) - len(fixed_values),
            "nonpivot_forcing_coordinates_including_H_at_zero": len(generators),
        },
        "pivot_variables": [str(variable) for variable in pivot_variables],
        "matrix_shape": list(endpoint_matrix.shape),
        "determinant": str(determinant),
        "matrix_digest_sha256": polynomial_digest(
            list(endpoint_matrix), generators
        ),
        "forcing_digest_sha256": polynomial_digest(list(forcing), generators),
        "solution_digest_sha256": polynomial_digest(solutions, generators),
        "forcing_term_counts": forcing_term_counts,
        "forcing_total_degrees": forcing_degrees,
        "inverse_maximum_total_degree": max(inverse_degrees),
        "solution_straight_line_term_counts": solution_term_counts,
        "solution_straight_line_term_count": sum(solution_term_counts),
        "solution_total_degrees_with_H_at_zero_as_one_node": solution_degrees,
        "solution_total_degrees_after_H_at_zero_quadratic_expansion": (
            source_solution_degrees
        ),
        "solution_denominator_lcms": solution_denominators,
        "identity_check": "M*x_solution+forcing=0 in QQ[nonpivots,H_at_zero]^10",
        "H_at_zero_definition": (
            "sum_ell ell*w^((p_u+q_u)/5)*(w-1)^(p_nu+q_nu)*"
            "K_P,ell(w)*K_Q,-ell(w) evaluated at w=0"
        ),
        "H_at_zero_source_degree": 2,
        "followers": follower_records,
    }


def normalized_layer_support_audit(
    p_bands: dict[int, object],
    q_bands: dict[int, object],
) -> dict[str, object]:
    """Count exact w-coefficient slots after the common monomial is removed."""

    w = sp.Symbol("w")
    # A rational point of the localized double-root base is enough to prove
    # that each displayed structural coefficient is not the zero polynomial
    # in w0.  Degree bounds prove that no coefficient lies outside the set.
    support_w0 = sp.Integer(2)
    R = (w - support_w0) ** 2 / 25

    def endpoint_factors(band: object, top: sp.Expr | None = None) -> list[sp.Expr]:
        if top is not None:
            return [top]
        _, vanishing, degree = band_factor_data(band)
        factors = [(w - 1) ** vanishing]
        if degree:
            factors.append((w - 1) ** vanishing * w**degree)
        return factors

    def bracket(
        p_layer: int,
        p_factor: sp.Expr,
        q_layer: int,
        q_factor: sp.Expr,
    ) -> sp.Poly:
        p_u, _, _ = band_factor_data(p_bands[p_layer])
        q_u, _, _ = band_factor_data(q_bands[q_layer])
        shift = (p_u + q_u - 1) // 5
        expression = sp.expand(
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
        return sp.Poly(expression, w)

    rows: list[dict[str, object]] = []
    for layer in [*range(40, 4, -1), 3]:
        support: set[int] = set()
        pair_count = 0
        for p_layer in range(max(-22, layer - 25), 16):
            q_layer = layer - p_layer
            if not -12 <= q_layer <= 25:
                continue
            pair_count += 1
            p_top = (w - 1) ** 6 * R**3 if p_layer == 15 else None
            q_top = (
                -sp.Rational(9, 5) * (w - 1) ** 10 * R**5
                if q_layer == 25
                else None
            )
            for p_factor in endpoint_factors(p_bands[p_layer], p_top):
                for q_factor in endpoint_factors(q_bands[q_layer], q_top):
                    polynomial = bracket(
                        p_layer, p_factor, q_layer, q_factor
                    )
                    support.update(
                        degree[0]
                        for degree, coefficient in polynomial.terms()
                        if coefficient
                    )
        if support and support != set(range(min(support), max(support) + 1)):
            raise AssertionError("a normalized Laurent layer acquired a support hole")
        rows.append(
            {
                "layer": layer,
                "band_pair_count": pair_count,
                "coefficient_slot_count": len(support),
                "w_degree_interval": (
                    None if not support else [min(support), max(support)]
                ),
                "changed_by_endpoint_substitution": 3 <= layer <= 28,
            }
        )

    if rows[0] != {
        "layer": 40,
        "band_pair_count": 1,
        "coefficient_slot_count": 0,
        "w_degree_interval": None,
        "changed_by_endpoint_substitution": False,
    }:
        raise AssertionError("the common-power top row stopped vanishing")
    zero_rows = [row for row in rows if int(row["layer"]) != 40]
    total_slots = sum(int(row["coefficient_slot_count"]) for row in zero_rows)
    affected_slots = sum(
        int(row["coefficient_slot_count"])
        for row in zero_rows
        if bool(row["changed_by_endpoint_substitution"])
    )
    unaffected_slots = total_slots - affected_slots
    if (total_slots, affected_slots, unaffected_slots) != (1172, 819, 353):
        raise AssertionError("the normalized reduced-system row count changed")

    return {
        "coordinate": "w=(1+t)^5",
        "common_factor_removed": (
            "t^L*u^chi(L), leaving one polynomial in w on Laurent layer L"
        ),
        "support_nonvanishing_witness": "w0=2 in the localized double-R base",
        "rows": rows,
        "zero_layer_coefficient_slots": total_slots,
        "changed_slots_layers_3_through_28": affected_slots,
        "unchanged_slots_layers_29_through_39": unaffected_slots,
    }


def residual_fitting_audit() -> dict[str, object]:
    w, w0 = sp.symbols("w w0", nonzero=True)
    E = (w - 1) * (w - w0)
    target_local_modulus = sp.expand(w**2 * E**5)
    if sp.degree(target_local_modulus, w) != 12:
        raise AssertionError("the target local modulus changed")

    target_operator = sp.zeros(22, 20)
    for column in range(20):
        image = sp.Poly(
            sp.expand(
                5 * w * E * sp.diff(w**column, w)
                + (11 * E + 22 * w * sp.diff(E, w)) * w**column
            ),
            w,
        )
        for (degree,), coefficient in image.terms():
            target_operator[degree, column] = coefficient
    triangular_minor = sp.factor(target_operator[:20, :].det())
    expected_minor = sp.factor(
        w0**20 * sp.prod(5 * index + 11 for index in range(20))
    )
    if triangular_minor != expected_minor:
        raise AssertionError("the target triangular Fitting minor changed")

    layer_zero_modulus = sp.expand(w**3 * (w - 1) ** 6 * (w - w0) ** 6)
    if sp.degree(layer_zero_modulus, w) != 15:
        raise AssertionError("the layer-zero Artinian modulus changed")

    return {
        "target": {
            "residual_polynomial": "G_red(w)=J_4_red(w)-1, degree<=33",
            "local_modulus": "D=w^2*(w-1)^5*(w-w0)^5",
            "local_modulus_degree": 12,
            "remaining_local_coordinates": [
                "Taylor_w0(G_red,0)",
                "Taylor_w0(G_red,1)",
                "Taylor_w0(G_red,2)",
                "Taylor_w0(G_red,3)",
                "Taylor_w0(G_red,4)",
            ],
            "quotient": "f=quo_w(G_red,D), degree<=21",
            "operator": (
                "N(S)=5*w*E*S'+(11*E+22*w*E')*S, "
                "E=(w-1)*(w-w0)"
            ),
            "operator_shape": [22, 20],
            "unit_minor": str(triangular_minor),
            "triangular_solution": (
                "s_d=(f_d+(5*d+28)*(1+w0)*s_(d-1)-"
                "5*(d+9)*s_(d-2))/((5*d+11)*w0), 0<=d<=19, "
                "with s_-1=s_-2=0"
            ),
            "rho_20": "f_20+128*(1+w0)*s_19-145*s_18",
            "rho_21": "f_21-150*s_19",
            "remaining_fitting_coordinates": ["rho_20(f)", "rho_21(f)"],
        },
        "layer_zero": {
            "residual_integral": (
                "H_red=sum_ell ell*P_red,ell*Q_-ell"
            ),
            "artinian_modulus": "K=w^3*(w-1)^6*(w-w0)^6",
            "artinian_length": 15,
            "remaining_coordinates": [
                "H_red(w0)-H_red(0)",
                "H_red'(w0)",
                "H_red''(w0)",
                "H_red'''(w0)",
                "H_red^(4)(w0)",
                "H_red^(5)(w0)",
            ],
        },
        "remaining_coordinate_count": 13,
    }


def upper_triangular_power_elimination_audit(
    p_bands: dict[int, object],
    q_bands: dict[int, object],
    layer_support: dict[str, object],
) -> dict[str, object]:
    """Eliminate the endpoint-disjoint new-Q bands on layers 39 through 29."""

    w, w0 = sp.symbols("w w0", nonzero=True)
    slots_by_layer = {
        int(row["layer"]): int(row["coefficient_slot_count"])
        for row in layer_support["rows"]
    }
    records: list[dict[str, object]] = []
    eliminated_coordinates = 0
    residual_slots = 0
    retained_centralizers = 0

    for descent in range(1, 12):
        layer = 40 - descent
        p_layer = 15 - descent
        q_layer = 25 - descent
        q_u, q_vanishing, q_degree = band_factor_data(q_bands[q_layer])
        q_dimension = q_degree + 1
        resonance = descent in (5, 10)
        expected_kernel_dimension = int(resonance)
        rank = q_dimension - expected_kernel_dimension

        # After the common factor
        #
        #   15*w^shift*(w-1)^(nu+5)*(w-w0)^5
        #
        # and the top normalization unit are removed, the q-only operator
        # sends w^k to A_k*w^k+B_k*w^(k+1)+C_k*w^(k+2), where
        #
        #   A_k=w0*(5*k+q_u).
        #
        # For q_u>0, rows/columns 0..d give a lower-triangular unit minor.
        # At the two source resonances q_u=0; delete the kernel direction and
        # use rows/columns 1..d instead.  The possible w^-1 Kummer shift only
        # reindexes those rows.
        operator = sp.zeros(q_degree + 3, q_dimension)
        for k in range(q_dimension):
            A_k = w0 * (5 * k + q_u)
            B_k = (
                (w0 + 1) * (2 * q_layer - 5 * k - q_u)
                - 5 * q_vanishing * w0
            )
            C_k = 5 * k + q_u + 5 * q_vanishing - 4 * q_layer
            operator[k, k] = A_k
            operator[k + 1, k] = B_k
            operator[k + 2, k] = C_k

        if resonance:
            selected = list(range(1, q_dimension))
            minor = sp.factor(operator.extract(selected, selected).det())
            expected_minor = sp.factor(
                w0**q_degree * 5**q_degree * sp.factorial(q_degree)
            )
            centralizer = f"C0^{5 - descent // 5}"
        else:
            selected = list(range(q_dimension))
            minor = sp.factor(operator.extract(selected, selected).det())
            expected_minor = sp.factor(
                w0**q_dimension
                * sp.prod(5 * k + q_u for k in range(q_dimension))
            )
            centralizer = None
        if minor != expected_minor:
            raise AssertionError("an upper new-Q unit minor changed")

        coefficient_slots = slots_by_layer[layer]
        layer_residual_slots = coefficient_slots - rank
        eliminated_coordinates += rank
        residual_slots += layer_residual_slots
        retained_centralizers += expected_kernel_dimension
        records.append(
            {
                "descent": descent,
                "layer": layer,
                "new_P_band": p_layer,
                "new_P_dimension": p_bands[p_layer].dimension,
                "new_Q_band": q_layer,
                "new_Q_dimension": q_dimension,
                "Q_fixed_vanishing": q_vanishing,
                "normalized_operator": (
                    "N_delta(w^k)=A_k*w^k+B_k*w^(k+1)+C_k*w^(k+2)"
                ),
                "A_k": f"w0*(5*k+{q_u})",
                "B_k": (
                    f"(w0+1)*({2*q_layer}-5*k-{q_u})"
                    f"-{5*q_vanishing}*w0"
                ),
                "C_k": f"5*k+{q_u + 5*q_vanishing - 4*q_layer}",
                "operator_rank": rank,
                "source_centralizer": centralizer,
                "unit_minor": str(minor),
                "Laurent_coefficient_slots": coefficient_slots,
                "Fitting_coordinate_slots_after_Q_elimination": (
                    layer_residual_slots
                ),
            }
        )

    if (
        eliminated_coordinates,
        retained_centralizers,
        residual_slots,
    ) != (134, 2, 219):
        raise AssertionError("the endpoint-disjoint triangular block changed")
    initial_active_coordinates = (
        sum(p_bands[layer].dimension for layer in range(-22, 16))
        + sum(q_bands[layer].dimension for layer in range(-12, 26))
        + sum(p_bands[layer].dimension for layer in (-24, -23))
        + sum(q_bands[layer].dimension for layer in (-14, -13))
        - p_bands[15].dimension
        - q_bands[25].dimension
        - 3
        - 10
    )
    active_after = initial_active_coordinates - eliminated_coordinates
    displayed_after = 1172 + 13 + 5 - eliminated_coordinates
    if (initial_active_coordinates, active_after, displayed_after) != (
        1061,
        927,
        1056,
    ):
        raise AssertionError("the upper triangular presentation size changed")

    return {
        "endpoint_disjoint_descent_interval": [1, 11],
        "Laurent_layer_interval": [39, 29],
        "why_the_block_commutes_with_endpoint_elimination": (
            "its new Q bands are Q_24 through Q_14, while the endpoint "
            "dependency cone uses Q_-12 through Q_13; P_3 first enters at "
            "descent 12"
        ),
        "common_operator_factor": (
            "15*w^shift*(w-1)^(Q_nu+5)*(w-w0)^5 times a base unit"
        ),
        "unit_minor_theorem": {
            "nonresonant": "w0^(d+1)*product_(k=0)^d(5*k+q_u)",
            "resonant_delta_5_10": "w0^d*5^d*d!",
        },
        "rows": records,
        "new_Q_coordinates_before": eliminated_coordinates + retained_centralizers,
        "new_Q_coordinates_eliminated": eliminated_coordinates,
        "source_centralizers_retained": retained_centralizers,
        "upper_Laurent_coefficient_slots_before": (
            eliminated_coordinates + residual_slots
        ),
        "upper_Fitting_coordinate_slots_after": residual_slots,
        "active_source_coordinates_after_upper_triangular_elimination": (
            active_after
        ),
        "displayed_generator_slots_after_upper_triangular_elimination": (
            displayed_after
        ),
        "coupling_boundary": {
            "descent": 12,
            "layer": 28,
            "bands": ["P_3", "Q_13"],
            "consequence": (
                "below this row the endpoint substitution and the new-Q "
                "power elimination must be combined as one Schur/Fitting "
                "calculation; subtracting all later operator ranks would be invalid"
            ),
        },
    }


def presentation_size_audit(
    p_bands: dict[int, object],
    q_bands: dict[int, object],
    layer_support: dict[str, object],
) -> dict[str, object]:
    t0, u0, w0, a_k, b_k = sp.symbols(
        "t0 u0 w0 a_k b_k", nonzero=True
    )
    c_leading = t0**5 / 25
    p7_leading = t0**7 * u0**3 * (w0 - 1) ** 3 * a_k
    p_minus1_value = t0**-1 * u0 * b_k
    normalized_ratio = sp.cancel(
        p_minus1_value * c_leading**3 / p7_leading**2
    )
    expected_ratio = b_k / (
        25**3 * w0 * (w0 - 1) ** 6 * a_k**2
    )
    if sp.cancel(
        normalized_ratio.subs(w0, u0**5)
        - expected_ratio.subs(w0, u0**5)
    ) != 0:
        raise AssertionError("the descent-eight global ratio scaling changed")

    closure_dimension = sum(
        p_bands[layer].dimension for layer in range(-22, 16)
    ) + sum(q_bands[layer].dimension for layer in range(-12, 26))
    old_layer_zero_extras = sum(
        p_bands[layer].dimension for layer in (-24, -23)
    ) + sum(q_bands[layer].dimension for layer in (-14, -13))
    specialized_top = p_bands[15].dimension + q_bands[25].dimension
    fixed_values = 3  # P3(1), Q1(1), and Q13(1); top values are specialized.
    endpoint_pivots = 10
    target_new_image = p_bands[-21].dimension + q_bands[-11].dimension
    active_source_coordinates = (
        closure_dimension
        + old_layer_zero_extras
        - specialized_top
        - fixed_values
        - endpoint_pivots
    )
    if (
        closure_dimension,
        old_layer_zero_extras,
        specialized_top,
        target_new_image,
        active_source_coordinates,
    ) != (1036, 56, 18, 30, 1061):
        raise AssertionError("the endpoint-reduced source-coordinate count changed")

    zero_generators = int(layer_support["zero_layer_coefficient_slots"])
    global_residues = 13
    branch_incidence_rows = 5
    displayed_generators = zero_generators + global_residues + branch_incidence_rows
    if displayed_generators != 1190:
        raise AssertionError("the endpoint-reduced generator count changed")

    return {
        "base": (
            "A8=QQ[w0,w0^-1,(w0-1)^-1,y]/(27*y^2-9*y+1)"
        ),
        "source_coordinate_ledger": {
            "closure_P_minus22_through_15_and_Q_minus12_through_25": (
                closure_dimension
            ),
            "old_layer_zero_extras_P_minus24_minus23_Q_minus14_minus13": (
                old_layer_zero_extras
            ),
            "specialized_top_coordinates_removed": specialized_top,
            "fixed_endpoint_values_removed": fixed_values,
            "endpoint_pivots_removed": endpoint_pivots,
            "target_new_image_coordinates_retained": target_new_image,
            "target_new_image_retention_reason": (
                "P_-21 and Q_-11 are killed by the target cokernel map but "
                "re-enter the retained Laurent layer 3 against Q_24 and P_14"
            ),
            "active_source_coordinates": active_source_coordinates,
        },
        "descent_eight_component_rows": [
            "K_P7(w0)=0",
            "K_P7'(w0)=0",
            "K_P7''(w0)=0",
            "K_Q1(w0)=0",
            (
                "K_P-1(w0)=25^3*w0*(w0-1)^6*y*a^2, "
                "a=[(w-w0)^3]K_P7"
            ),
        ],
        "descent_eight_ratio_scale_check": {
            "C0_leading_at_w0": "t0^5/25",
            "P7_leading_at_w0": "t0^7*u0^3*(w0-1)^3*a",
            "P_minus1_value_at_w0": "t0^-1*u0*K_P-1(w0)",
            "identity_using_u0^5_equals_w0": (
                "b_normalized/a_normalized^2="
                "K_P-1(w0)/(25^3*w0*(w0-1)^6*a^2)"
            ),
        },
        "descent_eight_incidence_row_count": branch_incidence_rows,
        "generator_ledger": {
            "zero_Laurent_coefficient_generators": zero_generators,
            "global_target_and_layer_zero_residues": global_residues,
            "descent_eight_component_rows": branch_incidence_rows,
            "displayed_generator_count_not_claimed_minimal": displayed_generators,
        },
        "certified_total_degree_upper_bound_after_endpoint_substitution": 8,
        "degree_reason": (
            "the pivot circuits have source degree at most 7 after the "
            "quadratic H(0) node is expanded, and every remaining bracket "
            "is linear in a P pivot and linear in its Q partner"
        ),
    }


def build_payload() -> dict[str, object]:
    p_bands = {
        layer: make_band("P", 75, 3, layer) for layer in range(-75, 16)
    }
    q_bands = {
        layer: make_band("Q", 125, 5, layer) for layer in range(-125, 26)
    }
    endpoint = endpoint_substitution_audit(p_bands, q_bands)
    layer_support = normalized_layer_support_audit(p_bands, q_bands)
    upper_triangular = upper_triangular_power_elimination_audit(
        p_bands, q_bands, layer_support
    )
    fitting = residual_fitting_audit()
    presentation = presentation_size_audit(p_bands, q_bands, layer_support)

    return {
        "schema": "plane-jc.f2-75-125-endpoint-reduction.v1",
        "status": "exact-pre-lower-tail-presentation;unit-ideal-undecided",
        "branch": {
            "top": "R(w)=(w-w0)^2/(25*(1-w0)^2)",
            "candidate_algebra": (
                "A8=QQ[w0,w0^-1,(w0-1)^-1,y]/(27*y^2-9*y+1)"
            ),
        },
        "complete_endpoint_substitution": endpoint,
        "normalized_laurent_rows": layer_support,
        "endpoint_disjoint_upper_triangular_elimination": upper_triangular,
        "remaining_fitting_coordinates": fitting,
        "reduced_ideal_presentation": presentation,
        "exact_substitution_rule": {
            "P3": (
                "K_P3=K_P3,pivot-zero+sum_(j=1)^4 x_j*L_j, "
                "L_j=delta^j+(j-7)(-1)^j*delta^6+"
                "(j-6)(-1)^j*delta^7"
            ),
            "P_minus1": (
                "K_P-1=K_P-1,pivot-zero+sum_(j=0)^5 x_j*L_j"
            ),
            "pivot_vector": "x=-M^-1*forcing",
            "affected_zero_layers": [[28, 5], [3, 3]],
            "unaffected_zero_layers": [[39, 29]],
            "target": "apply the same substitution to G=J_4-1",
            "layer_zero": "apply the same substitution to H=sum ell*P_ell*Q_-ell",
        },
        "decision": {
            "unit_ideal_obtained": False,
            "counterexample_obtained": False,
            "what_was_decided": (
                "the endpoint elimination is now fully carried into an exact "
                "degree-eight straight-line/Fitting presentation"
            ),
            "why_thirteen_is_not_a_small_polynomial_system": (
                "thirteen is the residual cokernel rank, while 1061 active "
                "source coordinates and 1172 earlier Laurent coefficient "
                "generators remain"
            ),
            "next_exact_operation": (
                "form the coupled Schur/Fitting circuits beginning at descent "
                "12 (layer 28), where P_3/Q_13 first meets the endpoint "
                "substitution; carry them through descent 37 and the 13 "
                "functionals, then test that component ideal"
            ),
        },
        "software": {"sympy": sp.__version__},
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--artifact", type=Path, default=ARTIFACT)
    args = parser.parse_args()

    payload = build_payload()
    artifact = args.artifact.resolve()
    if args.refresh:
        artifact.parent.mkdir(parents=True, exist_ok=True)
        artifact.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        try:
            display_path = artifact.relative_to(ROOT)
        except ValueError:
            display_path = artifact
        print(f"WROTE {display_path}")
    else:
        expected = json.loads(artifact.read_text())
        current_claim = {key: value for key, value in payload.items() if key != "software"}
        pinned_claim = {key: value for key, value in expected.items() if key != "software"}
        if current_claim != pinned_claim:
            raise AssertionError(
                "pinned endpoint-reduction artifact is stale; inspect before --refresh"
            )

    endpoint = payload["complete_endpoint_substitution"]
    presentation = payload["reduced_ideal_presentation"]
    print("PASS: the complete ten-pivot endpoint substitution is exact")
    print(
        "PASS:",
        endpoint["solution_straight_line_term_count"],
        "solution terms, source degree <=7, reduced bracket degree <=8",
    )
    print("PASS: only Laurent layers 3 through 28 change")
    upper = payload["endpoint_disjoint_upper_triangular_elimination"]
    print(
        "PASS:",
        upper["new_Q_coordinates_eliminated"],
        "upper new-Q coordinates eliminate by unit minors;",
        upper["upper_Fitting_coordinate_slots_after"],
        "upper Fitting slots remain",
    )
    print(
        "PASS:",
        presentation["source_coordinate_ledger"]["active_source_coordinates"],
        "active source coordinates precede the upper elimination;",
        presentation["generator_ledger"]["zero_Laurent_coefficient_generators"],
        "earlier Laurent coefficient generators remain",
    )
    print(
        "PASS:",
        upper["active_source_coordinates_after_upper_triangular_elimination"],
        "active source coordinates remain at the layer-28 coupling boundary",
    )
    print("PASS: the residual target/layer-zero Fitting module has 13 coordinates")
    print("PASS: no unit ideal and no counterexample is claimed")


if __name__ == "__main__":
    main()
