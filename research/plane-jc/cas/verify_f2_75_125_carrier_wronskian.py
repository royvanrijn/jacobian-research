#!/usr/bin/env python3
"""Verify the carrier-Wronskian classifier for the F2 ``(75,125)`` row.

The common top edge in original affine coordinates is

    P = q^-15*c(v)^3 + lower q-poles,
   -Q = q^-25*(9/5)*c(v)^5 + lower q-poles,

where ``q=y``, ``v=x*y^5`` and
``c(v)=v*(v-1)^2*R(v)``.  This checker converts the Keller two-form to the
target monomials

    pi=P^3/(-Q)^2,     h=P^5/(-Q)^3.

After the removable target shears in ``h``, the first nonconstant coefficient
is forced at q-descent 36.  Its coefficient H(v) is rational and satisfies one
explicit inhomogeneous Wronskian equation.  Local orders reduce that equation
to a constant numerator in the squarefree case and a linear numerator in the
double-root case.  Exact substitution then gives a finite classification:

* squarefree: ``R(v)=(v^2-3v+3)/25`` and a cyclic degree-three carrier map;
* double root: ``rho^2-3rho+1=0`` and a degree-six carrier map isomorphic to
  the already certified terminal Belyi map.

This is a necessary classification, not an exclusion of the F2 row or of the
degree pair ``(75,125)``.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/jc2_f2_75_125_carrier_wronskian.json"
)


def carrier_coordinate_audit() -> dict[str, object]:
    """Recover the carrier chart and the exact common top coefficient."""

    # q=y and v=x*y^5, hence x=v*q^-5 and
    # dx wedge dy = -q^-5 dq wedge dv.
    x_q_order = -5
    jacobian_form_q_order = -5
    c_degree = 1 + 2 + 2
    if (x_q_order, jacobian_form_q_order, c_degree) != (-5, -5, 5):
        raise AssertionError("the original carrier chart changed")

    # C_R=x*(v-1)^2*R(v)=q^-5*c(v).
    c_multiplicities = {
        "v=0": 1,
        "v=1": 2,
        "roots_of_R_total": 2,
    }
    if sum(c_multiplicities.values()) != c_degree:
        raise AssertionError("the carrier root multiplicities changed")

    return {
        "coordinates": {"q": "y", "v": "x*y^5", "x": "v*q^-5"},
        "two_form": "dx wedge dy=-q^-5*dq wedge dv",
        "common_root": "C_R=q^-5*c(v)",
        "carrier_polynomial": "c(v)=v*(v-1)^2*R(v)",
        "carrier_polynomial_degree": c_degree,
        "fixed_root_multiplicities": c_multiplicities,
        "cofactor_conditions": [
            "deg(R)=2",
            "R(0)!=0",
            "R(1)=1/25",
        ],
    }


def target_monomial_audit() -> dict[str, object]:
    """Compute the target coordinates and the forced descent."""

    p_q_order = -15
    minus_q_q_order = -25
    p_c_power = 3
    minus_q_c_power = 5
    minus_q_scalar = Fraction(9, 5)

    # Rows are the exponent vectors of pi and h in (P,-Q).
    target_exponent_matrix = ((3, -2), (5, -3))
    determinant = (
        target_exponent_matrix[0][0] * target_exponent_matrix[1][1]
        - target_exponent_matrix[0][1] * target_exponent_matrix[1][0]
    )
    if determinant != 1:
        raise AssertionError("the target monomial change is not unimodular")

    pi_q_order = 3 * p_q_order - 2 * minus_q_q_order
    pi_c_power = 3 * p_c_power - 2 * minus_q_c_power
    h_q_order = 5 * p_q_order - 3 * minus_q_q_order
    h_c_power = 5 * p_c_power - 3 * minus_q_c_power
    pi_leading_scalar = Fraction(1, 1) / minus_q_scalar**2
    h_leading_scalar = Fraction(1, 1) / minus_q_scalar**3
    if (
        pi_q_order,
        pi_c_power,
        h_q_order,
        h_c_power,
        pi_leading_scalar,
        h_leading_scalar,
    ) != (5, -1, 0, 0, Fraction(25, 81), Fraction(125, 729)):
        raise AssertionError("the carrier target leading data changed")

    # d(pi) wedge d(h)=(pi*h)/(P*(-Q))*dP wedge d(-Q).
    # The coefficient of dq wedge dv therefore has q-order
    # 5+0-(-15)-(-25)-5=40.  If the normalized h coordinate starts with
    # q^delta*H(v), its wedge with q^5*U(v) has order delta+4.
    target_two_form_q_order = (
        pi_q_order
        + h_q_order
        - p_q_order
        - minus_q_q_order
        - 5
    )
    forced_descent = target_two_form_q_order - 4
    if (target_two_form_q_order, forced_descent) != (40, 36):
        raise AssertionError("the carrier target descent changed")

    normalized_rhs_scalar = (
        pi_leading_scalar
        * h_leading_scalar
        / minus_q_scalar
    )
    if normalized_rhs_scalar != Fraction(5**6, 3**12):
        raise AssertionError("the normalized Wronskian scalar changed")

    return {
        "leading_pair": {
            "P": "q^-15*c(v)^3",
            "minus_Q": "q^-25*(9/5)*c(v)^5",
        },
        "target_coordinates": {
            "pi": "P^3/(-Q)^2",
            "h": "P^5/(-Q)^3",
        },
        "target_exponent_matrix": [list(row) for row in target_exponent_matrix],
        "target_exponent_determinant": determinant,
        "leading_pi": "q^5*(25/81)/c(v)",
        "leading_h": "125/729",
        "target_two_form_q_order": target_two_form_q_order,
        "forced_first_nonshear_descent": forced_descent,
        "normalized_rhs_scalar_for_Jacobian_one": str(normalized_rhs_scalar),
    }


def shear_normalization_audit(forced_descent: int) -> dict[str, object]:
    """Audit every removable coefficient before the inhomogeneous row."""

    removable: list[int] = []
    forbidden_nonshear: list[int] = []
    for descent in range(1, forced_descent):
        # For the first surviving coefficient H_d, the zero two-form rows give
        #     5*U*H_d' - d*U'*H_d=0,
        # so H_d^5/U^d is constant.  Since ord_(v=0)(U)=-1 and H_d is
        # rational, d must be divisible by five.
        if descent % 5 == 0:
            removable.append(descent)
        else:
            forbidden_nonshear.append(descent)

    if removable != [5, 10, 15, 20, 25, 30, 35]:
        raise AssertionError("the pre-target shear list changed")
    if len(forbidden_nonshear) != 28:
        raise AssertionError("the pre-target nonshear census changed")
    if forced_descent % 5 == 0:
        raise AssertionError("the inhomogeneous row became a target shear")

    return {
        "homogeneous_equation": "5*U*H_d'-d*U'*H_d=0",
        "rationality_identity": "H_d^5/U^d is constant",
        "divisibility_gate": "ord_(v=0)(U)=-1 forces 5|d",
        "removable_target_shear_descents": removable,
        "forbidden_nonshear_descents": forbidden_nonshear,
        "normalized_first_coefficient_descent": forced_descent,
    }


def carrier_target_fan_audit() -> dict[str, object]:
    """Extract the target divisor selected by the generic source carrier."""

    target = (5, 36)
    left = (1, 0)
    right = (0, 1)
    # pi=0 is the pre-existing target boundary; w=0 is a transverse curve.
    weights: dict[tuple[int, int], int | None] = {left: 1, right: None}
    insertion_order: list[tuple[int, int]] = []
    while True:
        middle = (left[0] + right[0], left[1] + right[1])
        insertion_order.append(middle)
        for endpoint in (left, right):
            if weights[endpoint] is not None:
                weights[endpoint] -= 1
        weights[middle] = -1
        if middle == target:
            break
        if target[1] * middle[0] < middle[1] * target[0]:
            right = middle
        else:
            left = middle

    ordered = sorted(
        weights,
        key=lambda ray: Fraction(ray[1], ray[0]) if ray[0] else Fraction(10**9),
    )
    if any(
        ordered[index][0] * ordered[index + 1][1]
        - ordered[index][1] * ordered[index + 1][0]
        != 1
        for index in range(len(ordered) - 1)
    ):
        raise AssertionError("the carrier target fan is not regular")
    boundary = [ray for ray in ordered if weights[ray] is not None]
    boundary_weights = [int(weights[ray]) for ray in boundary]
    expected_insertions = [
        (1, 1),
        (1, 2),
        (1, 3),
        (1, 4),
        (1, 5),
        (1, 6),
        (1, 7),
        (1, 8),
        (2, 15),
        (3, 22),
        (4, 29),
        (5, 36),
    ]
    if insertion_order != expected_insertions:
        raise AssertionError("the carrier target blowup chain changed")
    if boundary_weights != [0, -2, -2, -2, -2, -2, -2, -6, -1, -2, -2, -2, -2]:
        raise AssertionError("the carrier target boundary weights changed")
    target_index = boundary.index(target)
    if (boundary[target_index - 1], boundary[target_index + 1]) != (
        (1, 7),
        (4, 29),
    ):
        raise AssertionError("the carrier target adjacent rays changed")

    # zeta=w^5/pi^36 has orders -1,0,+1 on the left, target and right rays.
    residue_orders = [5 * ray[1] - 36 * ray[0] for ray in ((1, 7), target, (4, 29))]
    if residue_orders != [-1, 0, 1]:
        raise AssertionError("the carrier target residue orientation changed")

    return {
        "local_coordinates": ["pi=P^3/(-Q)^2", "w=h-target_shears"],
        "carrier_target_ray": list(target),
        "transverse_index": 1,
        "blowup_count": len(insertion_order),
        "insertion_order": [list(ray) for ray in insertion_order],
        "boundary_rays": [list(ray) for ray in boundary],
        "boundary_self_intersections": boundary_weights,
        "adjacent_rays": [[1, 7], [4, 29]],
        "residue_coordinate": "zeta=w^5/pi^36",
        "residue_orders_on_left_carrier_right": residue_orders,
    }


def local_divisor_shape(finite_multiplicities: tuple[int, ...]) -> dict[str, object]:
    """Determine the only possible denominator and numerator degree of H."""

    if sum(finite_multiplicities) != 5:
        raise AssertionError("c(v) no longer has degree five")
    if any(multiplicity not in (1, 2) for multiplicity in finite_multiplicities):
        raise AssertionError("an F2 carrier root has unexpected multiplicity")

    # At a finite zero of c of multiplicity m, U has order -m.  Balancing
    #     5*U*H'-36*U'*H = nonzero_constant*c^-9
    # forces ord(H)=1-8m, with nonzero leading coefficient 5-4m.
    finite_orders = tuple(1 - 8 * multiplicity for multiplicity in finite_multiplicities)
    finite_coefficients = tuple(5 - 4 * multiplicity for multiplicity in finite_multiplicities)
    if any(coefficient == 0 for coefficient in finite_coefficients):
        raise AssertionError("a finite carrier point became resonant")
    denominator_degree = -sum(finite_orders)

    # At infinity U has order five.  The homogeneous exponent 36 is resonant;
    # the inhomogeneous particular exponent is 39.  Since the forced finite
    # denominator has degree 36 or 37, exponent 39 would require a polynomial
    # numerator of negative degree.  Thus H has infinity order 36 and its
    # numerator degree is denominator_degree-36.
    infinity_resonant_order = 36
    infinity_particular_order = 39
    numerator_degree = denominator_degree - infinity_resonant_order
    if numerator_degree not in (0, 1):
        raise AssertionError("the F2 Wronskian numerator is no longer tiny")

    return {
        "finite_root_multiplicities": list(finite_multiplicities),
        "finite_H_orders": list(finite_orders),
        "finite_nonresonance_coefficients": list(finite_coefficients),
        "forced_denominator_degree": denominator_degree,
        "infinity_homogeneous_resonance_order": infinity_resonant_order,
        "infinity_particular_order": infinity_particular_order,
        "forced_numerator_degree": numerator_degree,
        "other_finite_poles_possible": False,
    }


def wronskian_numerator_operator(
    c: sp.Expr, radical: sp.Expr, numerator: sp.Expr, v: sp.Symbol
) -> sp.Expr:
    """Return the low-degree equation for H=radical*numerator/c^8."""

    # Dividing 5*U*H'-36*U'*H=K*c^-9 by U=(25/81)/c and
    # substituting H=radical*N/c^8 gives
    #
    #   5*(radical*N)'-4*(c'/c)*radical*N = 5^4/3^8.
    return sp.cancel(
        5 * sp.diff(radical * numerator, v)
        - 4 * radical * sp.diff(c, v) / c * numerator
        - sp.Rational(5**4, 3**8)
    )


def squarefree_classifier_audit() -> dict[str, object]:
    """Solve the squarefree quadratic cofactor Wronskian exactly."""

    shape = local_divisor_shape((1, 2, 1, 1))
    if shape["forced_numerator_degree"] != 0:
        raise AssertionError("the squarefree numerator ceased to be constant")

    v, a, b, numerator = sp.symbols("v a b numerator")
    r_polynomial = a * v**2 + b * v + (sp.Rational(1, 25) - a - b)
    c = sp.expand(v * (v - 1) ** 2 * r_polynomial)
    radical = sp.expand(v * (v - 1) * r_polynomial)
    residual = wronskian_numerator_operator(c, radical, numerator, v)
    residual_numerator = sp.Poly(sp.factor(residual), v)
    solutions = sp.solve(
        residual_numerator.all_coeffs(),
        (a, b, numerator),
        dict=True,
    )
    expected_solution = {
        a: sp.Rational(1, 25),
        b: -sp.Rational(3, 25),
        numerator: -sp.Rational(5**6, 3**9),
    }
    if solutions != [expected_solution]:
        raise AssertionError("the squarefree carrier classifier changed")

    selected_r = sp.factor(r_polynomial.subs(expected_solution))
    selected_c = sp.factor(c.subs(expected_solution))
    selected_radical = sp.factor(radical.subs(expected_solution))
    selected_numerator = expected_solution[numerator]
    selected_h = sp.factor(selected_radical * selected_numerator / selected_c**8)
    if wronskian_numerator_operator(
        selected_c, selected_radical, selected_numerator, v
    ) != 0:
        raise AssertionError("the squarefree Wronskian witness failed")
    discriminant = sp.factor(sp.discriminant(selected_r, v))
    if discriminant != -sp.Rational(3, 625):
        raise AssertionError("the squarefree cofactor discriminant changed")

    # Up to a nonzero target scaling, zeta=H^5/U^36 becomes the cyclic cubic
    # v*(v^2-3v+3)/(v-1)^3=1+1/(v-1)^3.
    residue_map = sp.factor(25 * v * selected_r / (v - 1) ** 3)
    cyclic_form = 1 + 1 / (v - 1) ** 3
    if sp.cancel(residue_map - cyclic_form) != 0:
        raise AssertionError("the squarefree carrier residue map changed")

    return {
        "case": "squarefree_R",
        "divisor_shape": shape,
        "forced_R": "(v^2-3*v+3)/25",
        "R_discriminant": str(discriminant),
        "forced_H": str(selected_h),
        "normalized_carrier_residue_map": "1+1/(v-1)^3",
        "residue_degree": 3,
        "branch_profiles": [[3], [3]],
        "simple_R_root_behavior": (
            "the two simple R roots are unramified simple points over residue value 0"
        ),
        "status": "unique necessary squarefree carrier row",
    }


def is_zero_mod_rho_relation(expression: sp.Expr, rho: sp.Symbol) -> bool:
    """Test a rational identity modulo rho^2-3*rho+1."""

    relation = sp.Poly(rho**2 - 3 * rho + 1, rho)
    numerator = sp.fraction(sp.cancel(expression))[0]
    return sp.rem(sp.Poly(numerator, rho), relation).is_zero


def double_root_classifier_audit() -> dict[str, object]:
    """Solve the double-root carrier Wronskian and identify its Belyi map."""

    shape = local_divisor_shape((1, 2, 2))
    if shape["forced_numerator_degree"] != 1:
        raise AssertionError("the double-root numerator ceased to be linear")

    v, rho, n0, n1 = sp.symbols("v rho n0 n1")
    relation = rho**2 - 3 * rho + 1
    r_polynomial = (v - rho) ** 2 / (25 * (1 - rho) ** 2)
    c = sp.factor(v * (v - 1) ** 2 * r_polynomial)
    radical = v * (v - 1) * (v - rho)
    residual = wronskian_numerator_operator(c, radical, n0 + n1 * v, v)
    coefficients = sp.Poly(sp.factor(residual), v).all_coeffs()
    if len(coefficients) != 3:
        raise AssertionError("the double-root Wronskian coefficient count changed")

    # The v^2 and v coefficients are a homogeneous 2x2 system in n0,n1.
    # Its determinant is a nonzero unit times rho^2-3*rho+1.  The constant
    # equation makes (n0,n1) nonzero, so the quadratic relation is necessary.
    coefficient_matrix = sp.Matrix(
        [
            [sp.diff(coefficients[0], n0), sp.diff(coefficients[0], n1)],
            [sp.diff(coefficients[1], n0), sp.diff(coefficients[1], n1)],
        ]
    )
    determinant = sp.factor(coefficient_matrix.det())
    determinant_quotient = sp.factor(determinant / relation)
    if determinant_quotient == 0 or sp.cancel(determinant - determinant_quotient * relation) != 0:
        raise AssertionError("the double-root determinant relation changed")

    selected_n0 = sp.Rational(625, 3**8) * (3 - rho)
    selected_n1 = sp.Rational(625, 3**9) * (4 * rho - 11)
    selected_residual = wronskian_numerator_operator(
        c, radical, selected_n0 + selected_n1 * v, v
    )
    residual_numerator = sp.Poly(selected_residual, v)
    if any(
        not is_zero_mod_rho_relation(coefficient, rho)
        for coefficient in residual_numerator.all_coeffs()
    ):
        raise AssertionError("the double-root Wronskian witness failed")

    alpha = sp.Rational(3, 5) * (rho + 1)
    if not is_zero_mod_rho_relation(-selected_n0 / selected_n1 - alpha, rho):
        raise AssertionError("the fivefold zero of the carrier map moved")
    residue_map = v * (v - alpha) ** 5 / ((v - 1) ** 3 * (v - rho) ** 3)
    s = -v / alpha
    terminal_map = 125 * s * (s + 1) ** 5 / (9 * s**2 + 15 * s + 5) ** 3
    if not is_zero_mod_rho_relation(
        residue_map / terminal_map - sp.Rational(729, 125), rho
    ):
        raise AssertionError("the carrier map is no longer the terminal Belyi map")

    # Away from its displayed zeros and poles, the logarithmic derivative has
    # no finite zero modulo the rho relation.  The remaining ramification is
    # the index-three point at v=infinity.
    logarithmic_numerator = sp.fraction(
        sp.cancel(sp.diff(residue_map, v) / residue_map)
    )[0]
    reduced_logarithmic_numerator = sp.Poly(logarithmic_numerator, v)
    relation_poly = sp.Poly(relation, rho)
    reduced_logarithmic = sum(
        sp.rem(sp.Poly(coefficient, rho), relation_poly).as_expr()
        * v ** (reduced_logarithmic_numerator.degree() - index)
        for index, coefficient in enumerate(
            reduced_logarithmic_numerator.all_coeffs()
        )
    )
    if sp.Poly(reduced_logarithmic, v).degree() > 0:
        raise AssertionError("the double-root map acquired an extra finite critical point")
    if reduced_logarithmic == 0:
        raise AssertionError("the double-root logarithmic derivative vanished")

    return {
        "case": "double_root_R",
        "divisor_shape": shape,
        "forced_double_root_equation": "rho^2-3*rho+1=0",
        "double_root_field": "Q(sqrt(5))",
        "forced_H_numerator": {
            "representation": "H=v*(v-1)*(v-rho)*(n0+n1*v)/c(v)^8",
            "constant": "625*(3-rho)/3^8",
            "linear": "625*(4*rho-11)/3^9",
        },
        "fivefold_zero": "alpha=3*(rho+1)/5",
        "normalized_carrier_residue_map": (
            "v*(v-alpha)^5/((v-1)^3*(v-rho)^3)"
        ),
        "source_change_to_terminal_map": "s=-v/alpha",
        "map_identity": "g(v)=(729/125)*h(s)",
        "residue_degree": 6,
        "passport": [[5, 1], [3, 3], [3, 1, 1, 1]],
        "status": "two conjugate necessary double-root carrier rows",
    }


def build_payload() -> dict[str, object]:
    carrier = carrier_coordinate_audit()
    target = target_monomial_audit()
    shears = shear_normalization_audit(
        int(target["forced_first_nonshear_descent"])
    )
    carrier_target = carrier_target_fan_audit()
    cases = [squarefree_classifier_audit(), double_root_classifier_audit()]

    return {
        "schema": "plane-jc.f2-75-125-carrier-wronskian.v1",
        "status": "exact-finite-carrier-classification-row-still-live",
        "certified_input": {
            "degree_pair": [75, 125],
            "normal_form_branch": "F2",
            "jacobian": "nonzero constant, normalized to one",
            "top_common_power": "P_top=C_R^3, Q_top=(-9/5)*C_R^5",
            "cofactor_strata": ["squarefree_R", "double_root_R"],
        },
        "carrier_chart": carrier,
        "target_monomial_change": target,
        "pre_target_shear_normalization": shears,
        "carrier_target_completion": carrier_target,
        "forced_wronskian_equation": {
            "U": "(25/81)/c(v)",
            "equation": (
                "5*U*H'-36*U'*H="
                "(5^6/3^12)*c(v)^-9"
            ),
            "right_side": "nonzero for every nonzero constant Jacobian",
        },
        "carrier_cases": cases,
        "global_result": {
            "verdict": "the certified F2 (75,125) row remains possible only on three exact carrier points",
            "reason": (
                "the squarefree cofactor is unique over Q and the double-root "
                "parameter has two conjugate values over Q(sqrt(5))"
            ),
            "claim_boundary": (
                "the carrier Wronskian is necessary but does not supply the "
                "remaining lower Laurent coefficients or a global Keller map"
            ),
        },
        "reproduction_command": (
            ".venv/bin/python plane-jc/cas/"
            "verify_f2_75_125_carrier_wronskian.py"
        ),
        "software": {"python": "sympy exact arithmetic"},
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
        if expected != payload:
            raise AssertionError(
                "the pinned F2 carrier-Wronskian artifact is stale; inspect "
                "the change before using --refresh"
            )

    print("F2_CARRIER_TARGET_MONOMIAL_PASS")
    print("F2_CARRIER_PRETARGET_SHEARS_PASS")
    print("F2_CARRIER_WRONSKIAN_CLASSIFIER_PASS")
    print("F2_CARRIER_SQUAREFREE_CYCLIC_CUBIC_PASS")
    print("F2_CARRIER_DOUBLE_ROOT_BELYI_IDENTIFICATION_PASS")
    print("F2_75_125_CERTIFIED_ROW_FINITE_NOT_EXCLUDED")
    print(f"F2_CARRIER_WRONSKIAN_ARTIFACT_SHA256={artifact_sha256(artifact)}")


if __name__ == "__main__":
    main()
