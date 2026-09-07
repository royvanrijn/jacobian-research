#!/usr/bin/env sage-python
"""Fit a native norm-eight genus-one pencil through four ICARM quotients.

This extends the curve-12 positive control to curve 395 in the 11952 chart
and curves 363, 364, and 378 in the 08f72 chart.  For every preferred basis
direction of the exact displayed exceptional quotient, it fits the cheapest
norm-eight residual-chord pencil and checks the incidence, cover witness,
quartic smoothness/irreducibility, and integral lift degree bounds exactly.
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
DIRECT = {
    "norm12-orbit-11952": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json",
    "norm12-orbit-08f72": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit08f72-direct-fibration-v1.json",
}
PRIORITY = {
    "norm12-orbit-11952": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.tsv",
    "norm12-orbit-08f72": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-08f72-alternate-norm8-pencil-priority-v1.tsv",
}
PRIORITY_CERT = {
    "norm12-orbit-11952": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.json",
    "norm12-orbit-08f72": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-08f72-alternate-norm8-pencil-priority-v1.json",
}
SWEEP = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-database-sweep-v1.json"
PUBLIC = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-public-fibres-v1.json"
QUOTIENT_AUDIT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-native-icarm-quotient-audit-v1.json"
HELPER = ROOT / "elkies-k3/scripts/construct_elkies_2026_bisections.sage"
CURVE12_REPLAY = ROOT / "elkies-k3/scripts/certify_r17_norm12_curve12_norm8_incidence.sage"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-norm8-incidence-v1.json"
CURVE_CHARTS = {
    395: "norm12-orbit-11952",
    363: "norm12-orbit-08f72",
    364: "norm12-orbit-08f72",
    378: "norm12-orbit-08f72",
}


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


def load_chart(chart, helper):
    direct = json.loads(DIRECT[chart].read_text())
    priority_certificate = json.loads(PRIORITY_CERT[chart].read_text())
    if direct["sections"]["status"] != "PASS_EXACT_SATURATED_RANK17_BASIS":
        raise ArithmeticError(f"{chart} generic basis is not saturated")
    if priority_certificate["status"] != "PASS_EXACT_COMPLETE_ALTERNATE_NORM8_PENCIL_PRIORITY":
        raise ArithmeticError(f"{chart} norm-eight priority table is not complete")
    if int(priority_certificate["class_count"]) != 63917:
        raise ArithmeticError(f"{chart} norm-eight class count changed")
    with PRIORITY[chart].open(newline="") as source:
        first = next(csv.DictReader(source, delimiter="\t"))
    if int(first["priority_rank"]) != 1:
        raise ArithmeticError(f"{chart} equation-cheapest trace changed")
    trace_word = [ZZ(value) for value in first["section_basis_w"].split()]
    if len(trace_word) != 17:
        raise ArithmeticError(f"{chart} norm-eight trace word has wrong length")

    ring = PolynomialRing(QQ, "u")
    field = ring.fraction_field()
    A = polynomial(ring, direct["weierstrass_model"]["A_coefficients_low_to_high"])
    B = polynomial(ring, direct["weierstrass_model"]["B_coefficients_low_to_high"])
    discriminant = ring(-16 * (4 * A**3 + 27 * B**2))
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
        raise ArithmeticError(f"{chart} selected trace is not finite-pole genus-one")
    return {
        "direct": direct,
        "first": first,
        "trace_word": trace_word,
        "ring": ring,
        "field": field,
        "A": A,
        "B": B,
        "discriminant": discriminant,
        "trace": trace,
        "h": h,
        "Nx": Nx,
        "Ny": Ny,
        "M0": M0,
    }


def fit_fibre(curve_id, chart, chart_data, sweep, public, quotient_audit, helper):
    direct = chart_data["direct"]
    first = chart_data["first"]
    trace_word = chart_data["trace_word"]
    ring, field = chart_data["ring"], chart_data["field"]
    A, B = chart_data["A"], chart_data["B"]
    discriminant = chart_data["discriminant"]
    trace = chart_data["trace"]
    h, Nx, Ny, M0 = (chart_data[key] for key in ("h", "Nx", "Ny", "M0"))

    hit = next(
        record
        for record in sweep["rational_j_hits_and_twists"]
        if int(record["curve_id"]) == curve_id and record["representative"] == chart
    )
    native = next(record for record in hit["native_chart_twists"] if record["chart"] == chart)
    parameter = QQ(native["native_parameter"]["numerator"]) / QQ(
        native["native_parameter"]["denominator"]
    )
    projective_scale_q = QQ(native["twist"]["quadratic_twist_parameter_q"])
    projective_scale_s = QQ(native["twist"]["qq_isomorphism_scale_s_with_s_squared_q"])
    if projective_scale_s**2 != projective_scale_q:
        raise ArithmeticError(f"curve {curve_id} projective isomorphism scale is inconsistent")

    public_record = next(record for record in public["records"] if int(record["id"]) == curve_id)
    quotient_record = next(record for record in quotient_audit["fibres"] if int(record["curve_id"]) == curve_id)
    if quotient_record["native_chart"] != chart:
        raise ArithmeticError(f"curve {curve_id} quotient/native chart mismatch")
    labels = quotient_record["displayed_exceptional_quotient"]["preferred_public_quotient_basis"]
    expected_rank = int(quotient_record["displayed_exceptional_quotient"]["free_rank"])
    if len(labels) != expected_rank:
        raise ArithmeticError(f"curve {curve_id} preferred quotient basis has wrong length")
    indices = [int(label[1:]) for label in labels]

    a1, a3, b2, target_A, target_B = short_invariants(public_record["ainvs"])
    fibre_A, fibre_B = A(parameter), B(parameter)
    scale_q = target_B * fibre_A / (fibre_B * target_A)
    if not scale_q.is_square():
        raise ArithmeticError(f"curve {curve_id} has a nontrivial affine-chart twist")
    scale_s = scale_q.sqrt()
    parameter_denominator = QQ(native["native_parameter"]["denominator"])
    if projective_scale_q * parameter_denominator**4 != scale_q:
        raise ArithmeticError(f"curve {curve_id} projective and affine q scales disagree")
    if projective_scale_s * parameter_denominator**2 not in (scale_s, -scale_s):
        raise ArithmeticError(f"curve {curve_id} projective and affine s scales disagree")
    if target_A != scale_q**2 * fibre_A or target_B != scale_q**3 * fibre_B:
        raise ArithmeticError(f"curve {curve_id} is not the asserted native fibre")

    target_curve = EllipticCurve(QQ, [target_A, target_B])
    fibre_curve = EllipticCurve(QQ, [fibre_A, fibre_B])
    public_target_points = [
        target_curve(QQ(x) + b2 / 12, QQ(y) + (a1 * QQ(x) + a3) / 2)
        for x, y in public_record["points"]
    ]
    public_fibre_points = [
        fibre_curve(point[0] / scale_q, point[1] / scale_s**3)
        for point in public_target_points
    ]

    trace_x, trace_y = QQ(trace[0](parameter)), QQ(trace[1](parameter))
    h_value = QQ(h(parameter))
    if not h_value:
        raise ArithmeticError(f"selected trace has a pole at curve {curve_id}")
    records = []
    normalized_quartics = []
    for direction_number, (label, public_index) in enumerate(zip(labels, indices), start=1):
        target = public_fibre_points[public_index - 1]
        target_x, target_y = target[0], target[1]
        if target_x == trace_x:
            raise ArithmeticError(f"curve {curve_id} vertical incidence degeneracy for {label}")
        slope = (target_y + trace_y) / (target_x - trace_x)
        pencil_parameter = (slope * h_value - QQ(M0(parameter))) / h_value**2
        M = M0 + pencil_parameter * h**2
        data = helper["chord_data_from_slope_numerator"](
            h, Nx, Ny, M, A, B, discriminant, ring, field, expected_q_degree=4
        )
        q = data["q"]
        if not q.is_irreducible() or q.gcd(discriminant).degree() or q.gcd(h).degree():
            raise ArithmeticError(f"curve {curve_id} {label} failed the smooth irreducible quartic gate")
        normalized_quartics.append(q.monic())
        sum_x, product_x = data["sum_x"], data["product_x"]
        cover_coordinate = (2 * target_x - QQ(sum_x(parameter))) / h_value
        line_value = QQ(M(parameter)) / h_value * target_x - (
            trace_y + QQ(M(parameter)) / h_value * trace_x
        )
        if cover_coordinate**2 != QQ(q(parameter)):
            raise ArithmeticError(f"curve {curve_id} {label} has no exact cover witness")
        if line_value != target_y:
            raise ArithmeticError(f"curve {curve_id} {label} misses the residual chord")
        if target_x**2 - QQ(sum_x(parameter)) * target_x + QQ(product_x(parameter)):
            raise ArithmeticError(f"curve {curve_id} {label} misses the residual quadratic")
        if not (
            data["x0"].degree() <= 4
            and data["x1"].degree() <= 2
            and data["y0"].degree() <= 6
            and data["y1"].degree() <= 4
        ):
            raise ArithmeticError(f"curve {curve_id} {label} lift exceeds integral degree bounds")
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
                "rational_cover_witness": {
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
    if len(set(normalized_quartics)) != len(records):
        raise ArithmeticError(f"curve {curve_id} produced duplicate normalized quartics")
    return {
        "curve_id": curve_id,
        "native_chart": chart,
        "native_parameter": rational_text(parameter),
        "displayed_exceptional_quotient": f"Z^{expected_rank}",
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
            "all_normalized_branch_quartics_distinct_within_fibre": True,
        },
        "generic_cover_consequence": {
            "base_curve_genus": 1,
            "base_change_degree": 2,
            "four_simple_branch_points_on_smooth_fibres": True,
            "lift_disjoint_from_zero": True,
            "anti_invariant_height": 16,
            "generic_mw_rank_lower_bound_on_each_fitted_cover": 18,
        },
    }


def build():
    sweep = json.loads(SWEEP.read_text())
    public = json.loads(PUBLIC.read_text())
    quotient_audit = json.loads(QUOTIENT_AUDIT.read_text())
    helper = runpy.run_path(str(HELPER))
    chart_data = {chart: load_chart(chart, helper) for chart in sorted(set(CURVE_CHARTS.values()))}
    fibres = [
        fit_fibre(curve_id, chart, chart_data[chart], sweep, public, quotient_audit, helper)
        for curve_id, chart in CURVE_CHARTS.items()
    ]
    total_directions = sum(record["incidence_signature"]["successful_directions"] for record in fibres)
    return {
        "schema": "elkies-k3.r17-norm12-icarm-norm8-incidence.v1",
        "status": "PASS_EXACT_FOUR_FIBRE_NORM8_GENUS_ONE_INCIDENCE",
        "fibres": fibres,
        "summary": {
            "curve_count": len(fibres),
            "successful_quotient_basis_directions": total_directions,
            "all_fitted_covers_genus_one": True,
            "all_anti_invariant_lifts_height": 16,
        },
        "claim_boundary": {
            "proved": [
                "every preferred displayed quotient-basis direction on curves 395,363,364,378 lies on an exact fitted member of the native equation-cheapest norm-eight genus-one pencil",
                "every fitted branch quartic is irreducible, squarefree, avoids the surface discriminant and trace poles, and has an exact rational witness at its source fibre",
            ],
            "not_proved": [
                "that one fixed pencil member exposes all directions",
                "that any fitted cover splits at another specialization",
                "that any displayed subgroup is the full Mordell--Weil group",
            ],
        },
        "inputs": {
            relative(path): digest(path)
            for path in (
                *DIRECT.values(),
                *PRIORITY.values(),
                *PRIORITY_CERT.values(),
                SWEEP,
                PUBLIC,
                QUOTIENT_AUDIT,
                HELPER,
                CURVE12_REPLAY,
            )
        },
        "software_assumptions": {"sage_version": SAGE_VERSION},
        "reproducing_command": (
            "sage -python elkies-k3/scripts/"
            "certify_r17_norm12_icarm_norm8_incidence.sage"
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    output = arguments.output.resolve()
    serialized = json.dumps(build(), indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not output.exists() or output.read_text() != serialized:
            raise ArithmeticError("stored four-fibre norm-eight incidence certificate differs from replay")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
    payload = json.loads(serialized)
    print(
        "R17ICARMNORM8|fibres={}|directions={}|genus=1|anti_invariant_height=16|status=PASS|output={}".format(
            payload["summary"]["curve_count"],
            payload["summary"]["successful_quotient_basis_directions"],
            relative(output),
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
