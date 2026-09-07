#!/usr/bin/env sage-python
"""Audit repeated late public labels on the five 074d9 fibres.

Labels P18 through P22 occur on all five fibres.  Curves 351, 356, 376, and
377 are the fixed reconstruction set and curve 385 is the untouched holdout.
For each label this script tests:

* the unique cubic polynomial x(u) through the four reconstruction fibres;
* every total-degree-three rational interpolation x=N(u)/D(u), with
  deg(N)+deg(D)=3 and a unique four-point projective solve;
* a quadratic inverse relation u=c0+c1*x+c2*x^2, overdetermined on the four
  reconstruction points before the holdout is inspected;
* the unique cubic inverse relation u=c0+c1*x+c2*x^2+c3*x^3;
* the unique quartic x(u) through all five fibres, followed by the exact
  Weierstrass-square test needed for an integral section.

The cubic inverse template has exactly four free coefficients, so the fifth
fibre is a genuine holdout.  A fully general cubic in x with affine-linear
coefficient functions has six free coefficients and is not identifiable from
four fibres; it is explicitly excluded rather than post-hoc normalized.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import runpy
import sys

from sage.all import PolynomialRing, QQ, matrix, vector
from sage.env import SAGE_VERSION


sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[2]
LINEAGE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json"
)
PUBLIC = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-public-fibres-v1.json"
)
CROSS_HELPER = (
    ROOT / "elkies-k3/scripts/certify_r17_074d9_cross_fibre_bisection_transfer.sage"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-074d9-late-point-interpolation-v1.json"
)

CHART = "norm12-orbit-074d9"
TRAIN_IDS = (351, 356, 376, 377)
HOLDOUT_ID = 385
TARGET_IDS = (*TRAIN_IDS, HOLDOUT_ID)
REPEATED_LABELS = tuple(range(18, 23))
STATUS = "PASS_EXACT_BOUNDED_LATE_POINT_INTERPOLATION_HOLDOUT"
PROTOCOL = "R17074D9LATEPOINTS"


def digest(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def rational_text(value) -> str:
    value = QQ(value)
    if value.denominator() == 1:
        return str(value.numerator())
    return f"{value.numerator()}/{value.denominator()}"


def polynomial_text(polynomial) -> list[str]:
    if not polynomial:
        return ["0"]
    return [rational_text(polynomial[index]) for index in range(polynomial.degree() + 1)]


def lagrange(points, ring):
    u = ring.gen()
    answer = ring.zero()
    for index, (abscissa, ordinate) in enumerate(points):
        basis = ring.one()
        denominator = QQ.one()
        for other_index, (other_abscissa, _other_ordinate) in enumerate(points):
            if other_index == index:
                continue
            basis *= u - other_abscissa
            denominator *= abscissa - other_abscissa
        answer += ordinate / denominator * basis
    return ring(answer)


def solve_inverse_relation(samples, degree):
    coefficient_matrix = matrix(
        QQ, [[x_value**power for power in range(degree + 1)] for _u, x_value in samples]
    )
    target = vector(QQ, [u_value for u_value, _x in samples])
    if coefficient_matrix.rank() < degree + 1:
        return None
    if len(samples) == degree + 1:
        coefficients = coefficient_matrix.solve_right(target)
        consistent = True
    else:
        square = coefficient_matrix[: degree + 1, :]
        coefficients = square.solve_right(target[: degree + 1])
        consistent = coefficient_matrix * coefficients == target
    return {
        "coefficients_low_to_high": [rational_text(value) for value in coefficients],
        "reconstruction_set_consistent": bool(consistent),
        "coefficients": coefficients,
    }


def rational_interpolants(samples, holdout, ring):
    u = ring.gen()
    results = []
    for numerator_degree in range(4):
        denominator_degree = 3 - numerator_degree
        rows = []
        for parameter, x_value in samples:
            rows.append(
                [parameter**power for power in range(numerator_degree + 1)]
                + [
                    -x_value * parameter**power
                    for power in range(denominator_degree + 1)
                ]
            )
        kernel = matrix(QQ, rows).right_kernel()
        record = {
            "numerator_degree_bound": numerator_degree,
            "denominator_degree_bound": denominator_degree,
            "kernel_dimension": int(kernel.dimension()),
            "unique_projective_reconstruction": kernel.dimension() == 1,
        }
        if kernel.dimension() != 1:
            record["holdout_hit"] = False
            record["status"] = "NONUNIQUE_OR_DEGENERATE_RECONSTRUCTION"
            results.append(record)
            continue
        row = vector(QQ, kernel.basis()[0])
        row_entries = list(row)
        numerator = ring(row_entries[: numerator_degree + 1])
        denominator = ring(row_entries[numerator_degree + 1 :])
        common = numerator.gcd(denominator)
        if common.degree() > 0:
            numerator //= common
            denominator //= common
        first = next(
            (value for value in reversed(list(numerator) + list(denominator)) if value),
            None,
        )
        if first is None:
            raise ArithmeticError("a rational interpolation kernel is zero")
        numerator /= first
        denominator /= first
        if any(not denominator(parameter) for parameter, _x in samples):
            record["holdout_hit"] = False
            record["status"] = "DENOMINATOR_ZERO_ON_RECONSTRUCTION_SET"
            results.append(record)
            continue
        holdout_parameter, holdout_x = holdout
        denominator_holdout = denominator(holdout_parameter)
        holdout_hit = bool(
            denominator_holdout
            and numerator(holdout_parameter) == holdout_x * denominator_holdout
        )
        monic_denominator = denominator / denominator.leading_coefficient()
        denominator_is_square = bool(
            monic_denominator.degree() % 2 == 0 and monic_denominator.is_square()
        )
        record.update(
            {
                "status": "EXACT_RECONSTRUCTION_AND_HOLDOUT_TEST",
                "numerator_coefficients_low_to_high": polynomial_text(numerator),
                "denominator_coefficients_low_to_high": polynomial_text(denominator),
                "reduced_denominator_is_scalar_times_polynomial_square": denominator_is_square,
                "holdout_denominator_nonzero": bool(denominator_holdout),
                "holdout_hit": holdout_hit,
                "holdout_residual": rational_text(
                    numerator(holdout_parameter) - holdout_x * denominator_holdout
                ),
            }
        )
        results.append(record)
    return results


def build_payload():
    lineage = json.loads(LINEAGE.read_text())
    public = json.loads(PUBLIC.read_text())
    cross_helper = runpy.run_path(str(CROSS_HELPER))
    if lineage["representative"]["chart"] != CHART:
        raise ArithmeticError("the lineage representative changed")
    parameters, isomorphisms = cross_helper["target_parameters"](lineage)
    public_by_id = {int(record["id"]): record for record in public["records"]}
    ring = PolynomialRing(QQ, "u")
    representative = lineage["representative"]
    A = ring([QQ(value) for value in representative["A_coefficients_low_to_high"]])
    B = ring([QQ(value) for value in representative["B_coefficients_low_to_high"]])

    chart_points = {}
    occurrence_inventory = {}
    for curve_id in TARGET_IDS:
        _invariants, target_curve, target_points = cross_helper["public_short_points"](
            public_by_id[curve_id]
        )
        target_iso = isomorphisms[curve_id]
        scale_q = QQ(target_iso["q"])
        scale_s = QQ(target_iso["s_with_s_squared_q"])
        parameter = parameters[curve_id]
        points = []
        for point in target_points:
            x_value = QQ(point[0] / scale_q)
            y_value = QQ(point[1] / scale_s**3)
            if y_value**2 != x_value**3 + A(parameter) * x_value + B(parameter):
                raise ArithmeticError("a late public point missed the 074d9 chart")
            points.append((x_value, y_value))
        chart_points[curve_id] = points
        for point_index in range(18, len(points) + 1):
            occurrence_inventory.setdefault(point_index, []).append(curve_id)

    label_results = []
    any_holdout_hit = False
    any_integral_section = False
    for point_index in REPEATED_LABELS:
        samples = [
            (parameters[curve_id], chart_points[curve_id][point_index - 1][0])
            for curve_id in TRAIN_IDS
        ]
        holdout = (
            parameters[HOLDOUT_ID],
            chart_points[HOLDOUT_ID][point_index - 1][0],
        )
        cubic_x = lagrange(samples, ring)
        cubic_x_residual = cubic_x(holdout[0]) - holdout[1]
        rational_results = rational_interpolants(samples, holdout, ring)

        quadratic = solve_inverse_relation(samples, 2)
        if quadratic is None:
            quadratic_record = {"status": "DEGENERATE_RECONSTRUCTION"}
        else:
            coefficients = quadratic.pop("coefficients")
            holdout_residual = (
                sum(coefficients[power] * holdout[1] ** power for power in range(3))
                - holdout[0]
            )
            quadratic_record = {
                **quadratic,
                "template": "u=c0+c1*x+c2*x^2",
                "holdout_residual": rational_text(holdout_residual),
                "holdout_hit": bool(
                    quadratic["reconstruction_set_consistent"] and not holdout_residual
                ),
            }

        cubic = solve_inverse_relation(samples, 3)
        if cubic is None:
            cubic_record = {"status": "DEGENERATE_RECONSTRUCTION"}
        else:
            coefficients = cubic.pop("coefficients")
            holdout_residual = (
                sum(coefficients[power] * holdout[1] ** power for power in range(4))
                - holdout[0]
            )
            cubic_record = {
                **cubic,
                "template": "u=c0+c1*x+c2*x^2+c3*x^3",
                "holdout_residual": rational_text(holdout_residual),
                "holdout_hit": not holdout_residual,
            }

        all_samples = samples + [holdout]
        quartic_x = lagrange(all_samples, ring)
        rhs = quartic_x**3 + A * quartic_x + B
        rhs_is_square = bool(rhs.is_square())
        section_record = {
            "template": "integral section with deg(x)<=4 on the polynomial 074d9 chart",
            "x_coefficients_low_to_high": polynomial_text(quartic_x),
            "weierstrass_rhs_is_polynomial_square": rhs_is_square,
            "section_exists": False,
        }
        if rhs_is_square:
            y_polynomial = rhs.sqrt()
            observed_y = [
                chart_points[curve_id][point_index - 1][1] for curve_id in TARGET_IDS
            ]
            plus = all(
                y_polynomial(parameters[curve_id]) == observed
                for curve_id, observed in zip(TARGET_IDS, observed_y)
            )
            minus = all(
                -y_polynomial(parameters[curve_id]) == observed
                for curve_id, observed in zip(TARGET_IDS, observed_y)
            )
            section_record.update(
                {
                    "y_coefficients_low_to_high": polynomial_text(y_polynomial),
                    "all_five_ordinate_signs_consistent": bool(plus or minus),
                    "section_exists": bool(plus or minus),
                }
            )
        label_hit = bool(
            not cubic_x_residual
            or any(record.get("holdout_hit") for record in rational_results)
            or quadratic_record.get("holdout_hit")
            or cubic_record.get("holdout_hit")
        )
        any_holdout_hit |= label_hit
        any_integral_section |= section_record["section_exists"]
        label_results.append(
            {
                "label": f"P{point_index}",
                "reconstruction_curve_ids": list(TRAIN_IDS),
                "holdout_curve_id": HOLDOUT_ID,
                "cubic_polynomial_x_of_u": {
                    "x_coefficients_low_to_high": polynomial_text(cubic_x),
                    "holdout_residual": rational_text(cubic_x_residual),
                    "holdout_hit": not cubic_x_residual,
                },
                "total_degree_three_rational_x_of_u": rational_results,
                "quadratic_inverse_relation": quadratic_record,
                "cubic_inverse_relation": cubic_record,
                "all_five_integral_section_test": section_record,
                "any_predeclared_holdout_hit": label_hit,
            }
        )
        print(
            f"{PROTOCOL}|label=P{point_index}|holdout_hit={label_hit}"
            f"|integral_section={section_record['section_exists']}",
            flush=True,
        )

    inputs = {
        relative(path): digest(path) for path in (LINEAGE, PUBLIC, CROSS_HELPER)
    }
    return {
        "schema": "elkies-k3.r17-074d9-late-point-interpolation.v1",
        "status": STATUS,
        "claim": (
            "Exact fixed-holdout rejection of the declared low-degree section and "
            "degree-two/degree-three x-relation templates for labels P18 through P22."
        ),
        "claim_boundary": (
            "A miss excludes only the displayed predeclared templates in the 074d9 "
            "coordinate. It does not exclude higher-degree rational functions, general "
            "multisections, a different base coordinate, isogenies, or historical label "
            "meaning. The fully general affine-coefficient cubic is not tested because "
            "four reconstruction fibres do not identify its six coefficients."
        ),
        "protocol": {
            "chart": CHART,
            "reconstruction_curve_ids": list(TRAIN_IDS),
            "holdout_curve_id": HOLDOUT_ID,
            "labels": [f"P{index}" for index in REPEATED_LABELS],
            "holdout_was_not_used_for_reconstruction": True,
        },
        "late_label_occurrence_inventory": {
            f"P{point_index}": curve_ids
            for point_index, curve_ids in sorted(occurrence_inventory.items())
        },
        "results": label_results,
        "summary": {
            "tested_label_count": len(label_results),
            "predeclared_holdout_hit_count": sum(
                bool(record["any_predeclared_holdout_hit"]) for record in label_results
            ),
            "integral_section_hit_count": sum(
                bool(record["all_five_integral_section_test"]["section_exists"])
                for record in label_results
            ),
            "all_declared_templates_miss": not any_holdout_hit and not any_integral_section,
            "lattice_trace_identification": (
                "NOT_TRIGGERED_NO_HIT"
                if not any_holdout_hit and not any_integral_section
                else "REQUIRED_FOR_REPORTED_HIT"
            ),
        },
        "nonidentifiable_template": {
            "template": (
                "x^3+(a0+a1*u)*x^2+(b0+b1*u)*x+(c0+c1*u)=0"
            ),
            "free_coefficient_count": 6,
            "reconstruction_equation_count": 4,
            "status": "NOT_IDENTIFIABLE_FROM_FOUR_RECONSTRUCTION_FIBRES",
        },
        "generation": {
            "command": (
                "sage -python elkies-k3/scripts/"
                "audit_r17_074d9_late_point_interpolation.sage"
            ),
            "checker_sha256": digest(Path(__file__)),
            "inputs": inputs,
        },
        "software_assumptions": {
            "sage_version": SAGE_VERSION,
            "required_features": ["exact QQ linear algebra", "exact polynomial squares"],
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    payload = build_payload()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not arguments.output.is_file() or arguments.output.read_text() != rendered:
            raise SystemExit("stale 074d9 late-point interpolation certificate")
        terminal = "PASS"
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered)
        terminal = "WROTE"
    print(
        f"{PROTOCOL}|labels={payload['summary']['tested_label_count']}"
        f"|holdout_hits={payload['summary']['predeclared_holdout_hit_count']}"
        f"|sections={payload['summary']['integral_section_hit_count']}"
        f"|status={terminal}|output={arguments.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
