#!/usr/bin/env sage-python
"""Fit one native alternate-Q80 norm-eight pencil through curve 12's Z^12 basis.

The exact curve-12 quotient certificate supplies twelve displayed complement
directions.  For the equation-cheapest norm-eight trace in the native 11952
chart, this replay solves the one-parameter residual-chord pencil incidence
equation for each direction and proves the resulting quartic cover is smooth,
irreducible, genus one, and carries the displayed rational witness.

This is a fitted positive-control signature, not evidence that any one fixed
pencil member splits at a different fibre.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import runpy
import sys

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ
from sage.env import SAGE_VERSION


sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[2]
DIRECT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
PRIORITY = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.tsv"
PRIORITY_CERT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.json"
SWEEP = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-database-sweep-v1.json"
PUBLIC = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-public-fibres-v1.json"
QUOTIENT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-curve12-alternate-q80-quotient-v1.json"
HELPER = ROOT / "elkies-k3/scripts/construct_elkies_2026_bisections.sage"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-curve12-norm8-incidence-v1.json"
CURVE_ID = 12
CHART = "norm12-orbit-11952"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def rational_text(value) -> str:
    value = QQ(value)
    if value.denominator() == 1:
        return str(value.numerator())
    return f"{value.numerator()}/{value.denominator()}"


def polynomial_coefficients(polynomial):
    if not polynomial:
        return ["0"]
    return [rational_text(polynomial[index]) for index in range(polynomial.degree() + 1)]


def polynomial(ring, values):
    return ring([QQ(value) for value in values])


def rational_function(ring, record):
    return ring.fraction_field()(
        polynomial(ring, record["numerator_coefficients_low_to_high"])
        / polynomial(ring, record["denominator_coefficients_low_to_high"])
    )


def short_invariants(ainvs):
    a1, a2, a3, a4, a6 = map(QQ, ainvs)
    b2 = a1**2 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3**2 + 4 * a6
    c4 = b2**2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    return a1, a3, b2, -c4 / 48, -c6 / 864


def build():
    direct = json.loads(DIRECT.read_text())
    priority_certificate = json.loads(PRIORITY_CERT.read_text())
    sweep = json.loads(SWEEP.read_text())
    public = json.loads(PUBLIC.read_text())
    quotient = json.loads(QUOTIENT.read_text())
    helper = runpy.run_path(str(HELPER))
    if direct["sections"]["status"] != "PASS_EXACT_SATURATED_RANK17_BASIS":
        raise ArithmeticError("the native 11952 generic basis is not saturated")
    if priority_certificate["status"] != "PASS_EXACT_COMPLETE_ALTERNATE_NORM8_PENCIL_PRIORITY":
        raise ArithmeticError("the native norm-eight priority table is not complete")
    if int(priority_certificate["class_count"]) != 63917:
        raise ArithmeticError("the native norm-eight class count changed")
    if quotient["displayed_exceptional_quotient"]["quotient"] != "Z^12":
        raise ArithmeticError("curve 12 lost its displayed Z^12 quotient")

    with PRIORITY.open(newline="") as source:
        first = next(csv.DictReader(source, delimiter="\t"))
    if int(first["priority_rank"]) != 1:
        raise ArithmeticError("the equation-cheapest trace changed")
    trace_word = [ZZ(value) for value in first["section_basis_w"].split()]
    if len(trace_word) != 17:
        raise ArithmeticError("the norm-eight trace word has wrong length")

    hit = next(
        record
        for record in sweep["rational_j_hits_and_twists"]
        if int(record["curve_id"]) == CURVE_ID and record["representative"] == CHART
    )
    native = next(
        record for record in hit["native_chart_twists"] if record["chart"] == CHART
    )
    parameter = QQ(native["native_parameter"]["numerator"]) / QQ(
        native["native_parameter"]["denominator"]
    )
    projective_scale_q = QQ(native["twist"]["quadratic_twist_parameter_q"])
    projective_scale_s = QQ(native["twist"]["qq_isomorphism_scale_s_with_s_squared_q"])
    if projective_scale_s**2 != projective_scale_q:
        raise ArithmeticError("the curve-12 projective isomorphism scale is inconsistent")

    public_record = next(record for record in public["records"] if int(record["id"]) == CURVE_ID)
    a1, a3, b2, target_A, target_B = short_invariants(public_record["ainvs"])
    ring = PolynomialRing(QQ, "u")
    field = ring.fraction_field()
    A = polynomial(ring, direct["weierstrass_model"]["A_coefficients_low_to_high"])
    B = polynomial(ring, direct["weierstrass_model"]["B_coefficients_low_to_high"])
    discriminant = ring(-16 * (4 * A**3 + 27 * B**2))
    fibre_A = A(parameter)
    fibre_B = B(parameter)
    scale_q = target_B * fibre_A / (fibre_B * target_A)
    if not scale_q.is_square():
        raise ArithmeticError("curve 12 has a nontrivial affine-chart twist")
    scale_s = scale_q.sqrt()
    parameter_denominator = QQ(native["native_parameter"]["denominator"])
    if projective_scale_q * parameter_denominator**4 != scale_q:
        raise ArithmeticError("curve-12 projective and affine q scales disagree")
    if projective_scale_s * parameter_denominator**2 not in (scale_s, -scale_s):
        raise ArithmeticError("curve-12 projective and affine s scales disagree")
    if target_A != scale_q**2 * fibre_A or target_B != scale_q**3 * fibre_B:
        raise ArithmeticError("curve 12 is not the asserted native 11952 fibre")

    generic_curve = EllipticCurve(field, [A, B])
    basis = [
        generic_curve(
            rational_function(ring, record["X"]),
            rational_function(ring, record["Y"]),
        )
        for record in direct["sections"]["records"]
    ]
    trace = sum(
        (int(coefficient) * point for coefficient, point in zip(trace_word, basis)),
        generic_curve(0),
    )
    frame = helper["trace_chord_frame"](trace[0], trace[1], ring)
    h, Nx, Ny, M0 = (frame[key] for key in ("h", "Nx", "Ny", "M0"))
    if h.degree() != 2:
        raise ArithmeticError("the selected norm-eight trace is not finite-pole genus-one")

    target_curve = EllipticCurve(QQ, [target_A, target_B])
    fibre_curve = EllipticCurve(QQ, [A(parameter), B(parameter)])
    public_target_points = [
        target_curve(
            QQ(x) + b2 / 12,
            QQ(y) + (a1 * QQ(x) + a3) / 2,
        )
        for x, y in public_record["points"]
    ]
    public_fibre_points = [
        fibre_curve(point[0] / scale_q, point[1] / scale_s**3)
        for point in public_target_points
    ]
    complement_labels = quotient["displayed_exceptional_quotient"][
        "free_basis_modulo_specialized_generic"
    ]
    complement_indices = [int(label[1:]) for label in complement_labels]
    if complement_indices != [2, 11, 4, 3, 6, 8, 17, 10, 28, 24, 19, 15]:
        raise ArithmeticError("the exact curve-12 complement basis changed")

    trace_x = QQ(trace[0](parameter))
    trace_y = QQ(trace[1](parameter))
    h_value = QQ(h(parameter))
    if not h_value:
        raise ArithmeticError("the selected trace has a pole at curve 12")
    records = []
    normalized_quartics = []
    for direction_number, (label, public_index) in enumerate(
        zip(complement_labels, complement_indices), start=1
    ):
        target = public_fibre_points[public_index - 1]
        target_x, target_y = target[0], target[1]
        if target_x == trace_x:
            raise ArithmeticError(f"vertical incidence degeneracy for {label}")
        slope = (target_y + trace_y) / (target_x - trace_x)
        pencil_parameter = (slope * h_value - QQ(M0(parameter))) / h_value**2
        M = M0 + pencil_parameter * h**2
        data = helper["chord_data_from_slope_numerator"](
            h,
            Nx,
            Ny,
            M,
            A,
            B,
            discriminant,
            ring,
            field,
            expected_q_degree=4,
        )
        q = data["q"]
        if not q.is_irreducible() or q.gcd(discriminant).degree() or q.gcd(h).degree():
            raise ArithmeticError(f"{label} failed the smooth irreducible quartic gate")
        normalized_quartics.append(q.monic())

        sum_x = data["sum_x"]
        product_x = data["product_x"]
        cover_coordinate = (2 * target_x - QQ(sum_x(parameter))) / h_value
        line_value = QQ(M(parameter)) / h_value * target_x - (
            trace_y + QQ(M(parameter)) / h_value * trace_x
        )
        if cover_coordinate**2 != QQ(q(parameter)):
            raise ArithmeticError(f"{label} has no exact cover witness")
        if line_value != target_y:
            raise ArithmeticError(f"{label} misses the residual chord")
        if target_x**2 - QQ(sum_x(parameter)) * target_x + QQ(product_x(parameter)):
            raise ArithmeticError(f"{label} misses the residual quadratic")
        if not (
            data["x0"].degree() <= 4
            and data["x1"].degree() <= 2
            and data["y0"].degree() <= 6
            and data["y1"].degree() <= 4
        ):
            raise ArithmeticError(f"{label} lift exceeds the integral degree bounds")

        records.append(
            {
                "direction_number": direction_number,
                "quotient_basis_direction": label,
                "source_public_point_index_one_based": public_index,
                "pencil_parameter_lambda": rational_text(pencil_parameter),
                "slope_numerator_M_coefficients_low_to_high": polynomial_coefficients(M),
                "branch_quartic_coefficients_low_to_high": polynomial_coefficients(q),
                "branch_quartic_degree": 4,
                "branch_quartic_irreducible_and_squarefree_over_Q": True,
                "branch_quartic_coprime_to_surface_discriminant_and_trace_denominator": True,
                "curve12_rational_cover_witness": {
                    "u": rational_text(parameter),
                    "s": rational_text(cover_coordinate),
                },
                "fibre_point": [rational_text(target_x), rational_text(target_y)],
                "residual_quadratic_sum_x_coefficients_low_to_high": polynomial_coefficients(sum_x),
                "residual_quadratic_product_x_coefficients_low_to_high": polynomial_coefficients(product_x),
                "lift_degree_bounds": {"x0": 4, "x1": 2, "y0": 6, "y1": 4},
                "exact_incidence_and_lift_identities_verified": True,
            }
        )

    if len(set(normalized_quartics)) != 12:
        raise ArithmeticError("two fitted directions produced the same branch quartic")
    return {
        "schema": "elkies-k3.r17-norm12-curve12-norm8-incidence.v1",
        "status": "PASS_EXACT_CURVE12_TWELVE_DIRECTION_NORM8_GENUS_ONE_INCIDENCE",
        "curve": {
            "icarm_id": 12,
            "native_chart": CHART,
            "native_parameter": rational_text(parameter),
            "displayed_exceptional_quotient": "Z^12",
        },
        "complete_norm_eight_frontier": {
            "class_count": 63917,
            "selected_priority_rank": 1,
            "selected_orbit_mask": int(first["orbit_mask"]),
            "selected_trace_word_in_saturated_generic_basis": [int(value) for value in trace_word],
            "selected_trace_height": 8,
            "selected_trace_denominator_degree": 2,
            "selected_trace_h_coefficients_low_to_high": polynomial_coefficients(h),
            "selected_trace_Nx_coefficients_low_to_high": polynomial_coefficients(Nx),
            "selected_trace_Ny_coefficients_low_to_high": polynomial_coefficients(Ny),
            "least_slope_M0_coefficients_low_to_high": polynomial_coefficients(M0),
        },
        "incidence_signature": {
            "directions": records,
            "successful_directions": len(records),
            "all_twelve_normalized_branch_quartics_distinct": True,
            "pairwise_distinct_quartics_give_pairwise_distinct_squareclasses": (
                "Each branch polynomial is irreducible of degree four; distinct monic "
                "representatives cannot have square ratio in QQ(u)."
            ),
        },
        "generic_cover_consequence": {
            "base_curve_genus": 1,
            "base_change_degree": 2,
            "four_simple_branch_points_on_smooth_fibres": True,
            "lift_disjoint_from_zero": True,
            "anti_invariant_height": 16,
            "generic_mw_rank_lower_bound_on_each_fitted_cover": 18,
        },
        "claim_boundary": {
            "proved": [
                "all twelve displayed quotient-basis directions lie on an exact fitted member of one native norm-eight genus-one pencil",
                "each fitted member has an irreducible squarefree quartic branch polynomial and an exact curve-12 rational witness",
                "the twelve branch quartics are pairwise distinct squareclasses",
            ],
            "not_proved": [
                "that one fixed pencil member exposes all twelve directions",
                "that these fitted covers split at any other specialization",
                "that the displayed subgroup is the full Mordell--Weil group",
            ],
        },
        "inputs": {
            relative(path): digest(path)
            for path in (DIRECT, PRIORITY, PRIORITY_CERT, SWEEP, PUBLIC, QUOTIENT, HELPER)
        },
        "software_assumptions": {"sage_version": SAGE_VERSION},
        "reproducing_command": (
            "sage -python elkies-k3/scripts/"
            "certify_r17_norm12_curve12_norm8_incidence.sage"
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    output = args.output.resolve()
    serialized = json.dumps(build(), indent=2, sort_keys=True) + "\n"
    if args.check:
        if not output.exists() or output.read_text() != serialized:
            raise ArithmeticError("stored norm-eight incidence certificate differs from replay")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
    print(
        "R17CURVE12NORM8|directions=12|genus=1|anti_invariant_height=16|"
        f"status=PASS|output={relative(output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
