#!/usr/bin/env sage-python
"""Specialize exact MW16 bases at the Nagao-prefilter finalists.

This is the bridge between the cheap local sieve and bounded half-lattice
recovery.  It reconstructs each parent quartic once, evaluates the selected
saturated generic MW16 basis at every finalist and groups exact Q-isomorphic
raw fibres.  Global minimization is intentionally deferred until a candidate
returns an independent direction, immediately before its residual-Selmer gate.

No point search, public complement, target jump label, or Selmer calculation
is performed here.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import sys

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "elliptic-curves/data/icarm_mw16_parent_ladder_blind_inputs_v1.json"
NAGAO = ROOT / "artifacts/generated-results/elliptic-curves/icarm_mw16_parent_nagao_prefilter_h300_v1.json"
MODEL = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
TABLE = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.tsv"
SCREEN = ROOT / "elkies-k3/scripts/screen_icarm_curve398_norm8_a1_fibrations.sage"
CHORD = ROOT / "elkies-k3/scripts/construct_elkies_2026_bisections.sage"
PREPARE = ROOT / "elliptic-curves/cas/prepare_icarm_mw16_parent_ladder_inputs.sage"
OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/icarm_mw16_nagao_finalist_specializations_h300_v1.json"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def qtext(value) -> str:
    value = QQ(value)
    return (
        str(value.numerator())
        if value.denominator() == 1
        else f"{value.numerator()}/{value.denominator()}"
    )


def point_record(point):
    if point.is_zero():
        return {"infinity": True}
    return {"x": qtext(point[0]), "y": qtext(point[1])}


def candidate_id(parent_id: str, numerator: int, denominator: int) -> str:
    sign = "m" if numerator < 0 else "p"
    return f"{parent_id}-{sign}{abs(numerator)}d{denominator}"


def invert_mobius(function, value, ring):
    numerator = ring(function.numerator())
    denominator = ring(function.denominator())
    if numerator.degree() > 1 or denominator.degree() > 1:
        raise ArithmeticError("stored old-section base map is not Mobius")
    n0, n1 = numerator[0], numerator[1]
    d0, d1 = denominator[0], denominator[1]
    bottom = value * d1 - n1
    if not bottom:
        raise ZeroDivisionError("selected section meets quartic infinity")
    answer = (n0 - value * d0) / bottom
    if function(answer) != value:
        raise ArithmeticError("stored Mobius inversion failed")
    return QQ(answer)


def specialize(
    *,
    parameter,
    source_points,
    base_maps,
    old_ring,
    h,
    nx,
    m0,
    quartic,
    child_a,
    child_b,
    prepare,
):
    specialized_a = QQ(child_a(parameter))
    specialized_b = QQ(child_b(parameter))
    if 4 * specialized_a**3 + 27 * specialized_b**2 == 0:
        return None, "SINGULAR_FIBRE"
    fixed_m = m0 + parameter * h**2
    fixed_quartic = old_ring(
        [QQ(quartic[index](parameter)) for index in range(5)]
    )
    sum_x = old_ring((fixed_m**2 - nx) // h**2)
    quartic_points = []
    try:
        for source_point, base_map in zip(source_points, base_maps):
            old_parameter = invert_mobius(base_map, parameter, old_ring)
            x_value = QQ(source_point[0](old_parameter))
            y_value = QQ(source_point[1](old_parameter))
            h_value = QQ(h(old_parameter))
            if not h_value:
                return None, "SELECTED_SECTION_MEETS_QUARTIC_INFINITY"
            w_value = (2 * x_value - QQ(sum_x(old_parameter))) / h_value
            if w_value**2 != fixed_quartic(old_parameter):
                raise ArithmeticError("selected old section missed the finalist quartic")
            quartic_points.append((old_parameter, w_value))
    except ZeroDivisionError:
        return None, "SELECTED_SECTION_SPECIALIZATION_UNDEFINED"

    t0, w0 = quartic_points[0]
    if not w0:
        return None, "POINTED_QUARTIC_ORIGIN_IS_BRANCH_POINT"
    shift_ring = PolynomialRing(QQ, "z")
    z = shift_ring.gen()
    shifted = shift_ring(fixed_quartic(t0 + z))
    ee, dd, cc, bb, aa = [QQ(shifted[index]) for index in range(5)]
    if ee != w0**2:
        raise ArithmeticError("pointed quartic constant term changed")
    a1g = dd / w0
    a2g = cc - dd**2 / (4 * w0**2)
    a3g = 2 * w0 * bb
    a4g = -4 * w0**2 * aa
    a6g = a2g * a4g
    b2g = a1g**2 + 4 * a2g
    b4g = a1g * a3g + 2 * a4g
    b6g = a3g**2 + 4 * a6g
    c4g = b2g**2 - 24 * b4g
    c6g = -b2g**3 + 36 * b2g * b4g - 216 * b6g
    if 81 * (-c4g / 48) != specialized_a or 729 * (-c6g / 864) != specialized_b:
        raise ArithmeticError("pointed quartic normalization missed the raw short model")

    raw_curve = EllipticCurve(QQ, [specialized_a, specialized_b])
    raw_points = []
    for old_parameter, w_value in quartic_points[1:]:
        zz = old_parameter - t0
        if not zz:
            return None, "GENERIC_SECTION_COLLIDES_WITH_SELECTED_ZERO"
        x_general = (2 * w0 * (w_value + w0) + dd * zz) / zz**2
        y_general = (
            4 * w0**2 * (w_value + w0)
            + 2 * w0 * dd * zz
            + (2 * w0 * cc - dd**2 / (2 * w0)) * zz**2
        ) / zz**3
        raw_points.append(
            raw_curve(
                9 * (x_general + b2g / 12),
                27 * (y_general + (a1g * x_general + a3g) / 2),
            )
        )
    if len(raw_points) != 16:
        raise ArithmeticError("finalist MW16 point list is incomplete")
    return {
        "raw_short_model": [
            "0",
            "0",
            "0",
            qtext(specialized_a),
            qtext(specialized_b),
        ],
        "raw_generic_points": [point_record(point) for point in raw_points],
        "j_invariant": qtext(raw_curve.j_invariant()),
    }, None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--nagao", type=Path, default=NAGAO)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    blind = json.loads(args.input.read_text())
    nagao = json.loads(args.nagao.read_text())
    if blind.get("status") != "PASS_EXACT_COMPLEMENT_BLIND_NINE_PARENT_INPUTS":
        raise ArithmeticError("MW16 parent inputs are not passing")
    if nagao.get("status") != "PASS_BOUNDED_NINE_PRESENTATION_NAGAO_PREFILTER":
        raise ArithmeticError("Nagao prefilter is not passing")
    if nagao["next_gate"]["stage"] != "bounded_half_lattice_jump_recovery":
        raise ArithmeticError("Nagao artifact does not authorize this stage")

    screen = SourceFileLoader("mw16_finalist_screen", str(SCREEN)).load_module()
    chord = SourceFileLoader("mw16_finalist_chord", str(CHORD)).load_module()
    prepare = SourceFileLoader("mw16_finalist_prepare", str(PREPARE)).load_module()
    source = json.loads(MODEL.read_text())
    table = screen.load_rows(TABLE)
    old_ring = PolynomialRing(QQ, "t")
    old_field = old_ring.fraction_field()
    old_a = old_ring(source["weierstrass_model"]["A_coefficients_low_to_high"])
    old_b = old_ring(source["weierstrass_model"]["B_coefficients_low_to_high"])
    old_curve = EllipticCurve(old_field, [old_a, old_b])
    old_basis = tuple(
        old_curve(
            screen.polynomial_from_record(record["X"], old_ring, QQ),
            screen.polynomial_from_record(record["Y"], old_ring, QQ),
        )
        for record in source["sections"]["records"]
    )
    parent_by_id = {row["parent_id"]: row for row in blind["parents"]}
    records = []
    failures = []
    for nagao_parent in nagao["presentations"]:
        parent = parent_by_id[nagao_parent["parent_id"]]
        trace_vector = vector(ZZ, parent["source_marking"]["trace_section_basis_w"])
        if list(trace_vector) != list(
            screen.parse_vector(
                table[int(parent["priority_rank"]) - 1]["section_basis_w"]
            )
        ):
            raise ArithmeticError("parent trace no longer matches the priority table")
        trace = sum(
            (
                coefficient * point
                for coefficient, point in zip(trace_vector, old_basis)
                if coefficient
            ),
            old_curve(0),
        )
        h, nx, m0, quartic, child_a, child_b = prepare.child_geometry(
            trace, old_a, old_b, old_ring, chord, screen
        )
        section_vectors = (
            vector(
                ZZ,
                parent["source_marking"][
                    "new_zero_source_section_basis_coordinates"
                ],
            ),
        ) + tuple(
            vector(ZZ, row)
            for row in parent["source_marking"][
                "generic_source_section_basis_coordinates"
            ]
        )
        source_points = tuple(
            sum(
                (
                    coefficient * point
                    for coefficient, point in zip(section_vector, old_basis)
                    if coefficient
                ),
                old_curve(0),
            )
            for section_vector in section_vectors
        )
        base_maps = tuple(
            old_field(
                old_ring(
                    [
                        QQ(value)
                        for value in record[
                            "numerator_coefficients_low_to_high"
                        ]
                    ]
                )
            )
            / old_field(
                old_ring(
                    [
                        QQ(value)
                        for value in record[
                            "denominator_coefficients_low_to_high"
                        ]
                    ]
                )
            )
            for record in parent["source_marking"]["base_maps_lambda_of_old_t"]
        )
        if len(source_points) != 17 or len(base_maps) != 17:
            raise ArithmeticError("parent section context is incomplete")
        for finalist in nagao_parent["finalists"]:
            numerator, denominator = map(int, finalist["projective_pair"])
            if denominator == 0:
                failures.append(
                    {
                        "parent_id": parent["parent_id"],
                        "parameter": "infinity",
                        "reason": "INFINITY_CHART_NOT_IMPLEMENTED",
                    }
                )
                continue
            parameter = QQ(numerator) / QQ(denominator)
            specialization, failure = specialize(
                parameter=parameter,
                source_points=source_points,
                base_maps=base_maps,
                old_ring=old_ring,
                h=h,
                nx=nx,
                m0=m0,
                quartic=quartic,
                child_a=child_a,
                child_b=child_b,
                prepare=prepare,
            )
            if failure is not None:
                failures.append(
                    {
                        "parent_id": parent["parent_id"],
                        "parameter": finalist["parameter"],
                        "reason": failure,
                    }
                )
                continue
            identifier = candidate_id(parent["parent_id"], numerator, denominator)
            records.append(
                {
                    "candidate_id": identifier,
                    "parent_id": parent["parent_id"],
                    "parent_curve_id": int(parent["curve_id"]),
                    "parameter": finalist["parameter"],
                    "projective_pair": [numerator, denominator],
                    "nagao": finalist,
                    "generic_height_gram": parent["generic_height_gram"],
                    "generic_rank": 16,
                    **specialization,
                }
            )
        print(
            f"MW16SPECIALIZE|parent={parent['parent_id']}|"
            f"requested={len(nagao_parent['finalists'])}|status=COMPLETE",
            flush=True,
        )

    groups = []
    for record in records:
        curve = EllipticCurve(QQ, [QQ(value) for value in record["raw_short_model"]])
        for group in groups:
            if (
                record["j_invariant"] == group["j_invariant"]
                and curve.is_isomorphic(group["representative_curve"])
            ):
                group["candidate_ids"].append(record["candidate_id"])
                break
        else:
            groups.append(
                {
                    "j_invariant": record["j_invariant"],
                    "representative_curve": curve,
                    "candidate_ids": [record["candidate_id"]],
                }
            )
    exact_classes = []
    for index, group in enumerate(groups, 1):
        class_id = f"QISO-{index:04d}"
        candidate_ids = sorted(group["candidate_ids"])
        exact_classes.append(
            {
                "class_id": class_id,
                "j_invariant": group["j_invariant"],
                "candidate_ids": candidate_ids,
            }
        )
        for record in records:
            if record["candidate_id"] in candidate_ids:
                record["q_isomorphism_class_id"] = class_id
    duplicate_groups = [row for row in exact_classes if len(row["candidate_ids"]) > 1]
    payload = {
        "schema": "elliptic-curves.icarm-mw16-nagao-finalist-specializations.v1",
        "status": "PASS_EXACT_MW16_NAGAO_FINALIST_SPECIALIZATIONS",
        "requested_finalist_count": sum(
            len(row["finalists"]) for row in nagao["presentations"]
        ),
        "successful_specialization_count": len(records),
        "structural_failure_count": len(failures),
        "exact_q_isomorphism_class_count": len(groups),
        "q_isomorphism_classes": exact_classes,
        "duplicate_q_isomorphism_groups": duplicate_groups,
        "structural_failures": failures,
        "candidates": records,
        "next_gate": {
            "stage": "bounded_half_lattice_jump_recovery",
            "authorized_search": (
                "complete exact maximum-depth MW16 stratum at reduced-coordinate "
                "height at most 100000 and at most 15 seconds per chart"
            ),
            "after_positive_recovery": (
                "compute the global minimal model, transport MW16 plus the recovered "
                "points, then run the complete residual 2-Selmer gate before any "
                "adaptive or unrestricted expensive continuation"
            ),
        },
        "inputs": {
            relative(path): digest(path)
            for path in (
                args.input,
                args.nagao,
                MODEL,
                TABLE,
                SCREEN,
                CHORD,
                PREPARE,
                Path(__file__),
            )
        },
        "claim_boundary": [
            "Every successful row has sixteen exact specialized generic points on its raw short model.",
            "The generic height form is the saturated function-field MW16 form; specialization independence is checked at the next gate.",
            "Structural failures are excluded from this chart implementation, not asserted to have low rank.",
            "Exact Q-isomorphism deduplication uses the raw fibres; global minimization is deferred until after positive half-lattice recovery.",
            "No point search, public-complement comparison, rank jump, minimization, or Selmer computation occurs here.",
        ],
        "reproducing_command": (
            "sage -python elliptic-curves/cas/specialize_icarm_mw16_nagao_finalists.sage --check"
        ),
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.is_file() or args.output.read_text() != serialized:
            raise ArithmeticError("stored finalist specialization artifact differs from replay")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    print(
        f"MW16SPECIALIZE|requested={payload['requested_finalist_count']}|"
        f"successful={len(records)}|q_isomorphism_classes={len(groups)}|"
        f"failures={len(failures)}|output={relative(args.output)}|status=PASS",
        flush=True,
    )


if __name__ == "__main__":
    main()
