#!/usr/bin/env sage-python
"""Fit and freeze one target-blind 074d9 norm-eight pencil across record fibres.

The selected trace is the first representative under a complete norm-eight
parity-coset enumeration ranked only by equation-basis coefficient cost.  No
target fibre enters that selection.  Twelve pencil members are fitted through
the displayed quotient basis at curve 356, frozen, and evaluated at curves
351, 376, 377, and 385.  The experiment is then repeated from curve 385.

Curve 12 belongs to the inequivalent 11952 fibration.  A specialization of a
074d9 base cover at "curve 12" is therefore undefined; the output records this
as not applicable instead of converting it into a false zero-transfer datum.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import runpy
import sys

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, matrix, vector
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
ENUMERATOR = ROOT / "elkies-k3/scripts/enumerate_rootless_bisection_orbits.sage"
CHORD_HELPER = ROOT / "elkies-k3/scripts/construct_elkies_2026_bisections.sage"
CROSS_HELPER = (
    ROOT / "elkies-k3/scripts/certify_r17_074d9_cross_fibre_bisection_transfer.sage"
)
QUOTIENT_HELPER = (
    ROOT / "elliptic-curves/scripts/evaluate_elkies_2026_bisections_at_controls.py"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-074d9-norm8-cross-fibre-transfer-v1.json"
)

CHART = "norm12-orbit-074d9"
TARGET_IDS = (351, 356, 376, 377, 385)
SOURCE_IDS = (356, 385)
STATUS = "PASS_EXACT_CANONICAL_NORM8_CROSS_FIBRE_TRANSFER"
PROTOCOL = "R17074D9NORM8TRANSFER"


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


def point_text(point) -> list[str]:
    return [rational_text(point[0]), rational_text(point[1])]


def parity_mask(value) -> int:
    return sum((int(entry) & 1) << index for index, entry in enumerate(value))


def word_score(value):
    oriented = min(tuple(map(int, value)), tuple(-int(entry) for entry in value))
    return (
        sum(abs(entry) for entry in oriented),
        sum(bool(entry) for entry in oriented),
        max(abs(entry) for entry in oriented),
        oriented,
    )


def canonical_norm8_trace(gram, enumerator):
    short_change = gram.LLL_gram().transpose()
    if abs(short_change.det()) != 1:
        raise ArithmeticError("the LLL coordinate change is not unimodular")
    short_gram = short_change * gram * short_change.transpose()

    def representative_key(short_tuple):
        return word_score(vector(ZZ, short_tuple) * short_change)

    streaming = enumerator["streaming_short_vectors"](
        short_gram, bound=8, representative_key=representative_key
    )
    excluded = (
        streaming["masks_by_norm"][2]
        | streaming["masks_by_norm"][4]
        | streaming["masks_by_norm"][6]
    )
    candidates = set(streaming["representatives"]) - excluded
    if len(candidates) != 63925:
        raise ArithmeticError("the published-R17 norm-eight class count changed")
    rows = []
    seen_masks = set()
    for short_mask in candidates:
        short_vector = vector(ZZ, streaming["representatives"][short_mask])
        section_vector = short_vector * short_change
        if section_vector * gram * section_vector != 8:
            raise ArithmeticError("a selected representative lost norm eight")
        mask = parity_mask(section_vector)
        if mask in seen_masks:
            raise ArithmeticError("the norm-eight enumeration repeated a parity class")
        seen_masks.add(mask)
        score = word_score(section_vector)
        rows.append(
            {
                "orbit_mask": mask,
                "short_orbit_mask": int(short_mask),
                "word": score[3],
                "score": score,
                "minimum_unoriented_count": int(
                    streaming["unoriented_multiplicities"][short_mask]
                ),
            }
        )
    rows.sort(key=lambda row: (row["score"], row["orbit_mask"]))
    return rows[0], rows[:64], streaming


def target_data(lineage, public, cross_helper):
    parameters, isomorphisms = cross_helper["target_parameters"](lineage)
    public_by_id = {int(record["id"]): record for record in public["records"]}
    result = {}
    for curve_id in TARGET_IDS:
        invariants, target_curve, target_points = cross_helper["public_short_points"](
            public_by_id[curve_id]
        )
        target_iso = isomorphisms[curve_id]
        scale_q = QQ(target_iso["q"])
        scale_s = QQ(target_iso["s_with_s_squared_q"])
        if scale_s**2 != scale_q:
            raise ArithmeticError("a 074d9 target scale is inconsistent")
        result[curve_id] = {
            "parameter": parameters[curve_id],
            "isomorphism": target_iso,
            "scale_q": scale_q,
            "scale_s": scale_s,
            "target_curve": target_curve,
            "target_points": target_points,
            "invariants": invariants,
        }
    return result


def fit_source_covers(
    *, source_id, trace_rank, trace, frame, A, B, discriminant, ring, field, targets, chord_helper
):
    source = targets[source_id]
    parameter = source["parameter"]
    scale_q = source["scale_q"]
    scale_s = source["scale_s"]
    target_curve = source["target_curve"]
    target_points = source["target_points"]
    fibre_curve = EllipticCurve(QQ, [A(parameter), B(parameter)])
    trace_point = fibre_curve(trace[0](parameter), trace[1](parameter))
    h, Nx, Ny, M0 = (frame[key] for key in ("h", "Nx", "Ny", "M0"))
    if h.degree() != 2 or not h(parameter):
        raise ArithmeticError("the canonical norm-eight trace is not regular at a source")
    covers = []
    for direction, public_index in enumerate(range(18, 30), start=1):
        public_point = target_points[public_index - 1]
        point = fibre_curve(public_point[0] / scale_q, public_point[1] / scale_s**3)
        if point[0] == trace_point[0]:
            raise ArithmeticError("a norm-eight source incidence is vertical")
        slope = (point[1] + trace_point[1]) / (point[0] - trace_point[0])
        h_value = QQ(h(parameter))
        pencil_parameter = (slope * h_value - QQ(M0(parameter))) / h_value**2
        M = M0 + pencil_parameter * h**2
        data = chord_helper["chord_data_from_slope_numerator"](
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
        q = ring(data["q"])
        if q.gcd(q.derivative()).degree() or q.gcd(discriminant).degree() or q.gcd(h).degree():
            raise ArithmeticError("a fitted norm-eight cover is not smooth and transverse")
        witness = (2 * point[0] - data["sum_x"](parameter)) / h_value
        if witness**2 != q(parameter):
            raise ArithmeticError("the fitted source point missed its quartic cover")
        if (
            point[0] ** 2
            - data["sum_x"](parameter) * point[0]
            + data["product_x"](parameter)
        ):
            raise ArithmeticError("the fitted source point missed the residual quadratic")
        covers.append(
            {
                "cover_id": f"trace-{trace_rank}-from-{source_id}-P{public_index}",
                "source_curve_id": source_id,
                "source_quotient_basis_direction": f"P{public_index}",
                "source_direction_number": direction,
                "pencil_parameter_lambda": rational_text(pencil_parameter),
                "slope_numerator_M_coefficients_low_to_high": polynomial_text(M),
                "branch_quartic_coefficients_low_to_high": polynomial_text(q),
                "branch_quartic_irreducible_over_Q": bool(q.is_irreducible()),
                "branch_quartic_squarefree_and_coprime_to_discriminant_and_trace_denominator": True,
                "source_cover_witness": {
                    "u": rational_text(parameter),
                    "s": rational_text(witness),
                    "public_point": f"P{public_index}",
                },
                "_data": data,
            }
        )
    monic_quartics = {
        tuple(q / q.leading_coefficient())
        for q in (
            ring(cover["_data"]["q"]) for cover in covers
        )
    }
    if len(monic_quartics) != 12:
        raise ArithmeticError("two source directions fitted the same norm-eight cover")
    return covers


def evaluate_frozen_covers(covers, trace, A, B, ring, targets, cross_helper):
    split_records = {curve_id: [] for curve_id in TARGET_IDS}
    evaluations = []
    for cover in covers:
        data = cover["_data"]
        q = ring(data["q"])
        cover_evaluations = []
        for curve_id in TARGET_IDS:
            parameter = targets[curve_id]["parameter"]
            value = QQ(q(parameter))
            square_root = cross_helper["exact_square_root"](value)
            if square_root is None:
                cover_evaluations.append(
                    {
                        "curve_id": curve_id,
                        "status": "EXACT_NONSQUARE",
                        "q_at_fibre": rational_text(value),
                    }
                )
                continue
            positive = (
                QQ(data["x0"](parameter) + data["x1"](parameter) * square_root),
                QQ(data["y0"](parameter) + data["y1"](parameter) * square_root),
            )
            negative = (
                QQ(data["x0"](parameter) - data["x1"](parameter) * square_root),
                QQ(data["y0"](parameter) - data["y1"](parameter) * square_root),
            )
            fibre_curve = EllipticCurve(QQ, [A(parameter), B(parameter)])
            positive_point = fibre_curve(positive)
            negative_point = fibre_curve(negative)
            if positive_point + negative_point != fibre_curve(
                trace[0](parameter), trace[1](parameter)
            ):
                raise ArithmeticError("a frozen norm-eight split lost its trace")
            record = {
                "cover_id": cover["cover_id"],
                "label": cover["cover_id"],
                "source_curve_id": cover["source_curve_id"],
                "source_quotient_basis_direction": cover[
                    "source_quotient_basis_direction"
                ],
                "q_at_fibre": rational_text(value),
                "canonical_square_root": rational_text(square_root),
                "positive_chart_point": point_text(positive),
                "negative_chart_point": point_text(negative),
                "exact_cover_and_branch_identities_verified": True,
            }
            split_records[curve_id].append(record)
            cover_evaluations.append(
                {
                    "curve_id": curve_id,
                    "status": "SPLIT_OVER_Q",
                    "q_at_fibre": rational_text(value),
                    "canonical_square_root": rational_text(square_root),
                }
            )
        public_cover = {key: value for key, value in cover.items() if key != "_data"}
        public_cover["frozen_evaluations"] = cover_evaluations
        evaluations.append(public_cover)
    return split_records, evaluations


def transfer_row(source_id, quotient_results):
    return {
        str(target_id): {
            "status": "EXACT",
            "split_cover_count": quotient_results[target_id]["split_count"],
            "T_rank": quotient_results[target_id]["class_span_rank"],
        }
        for target_id in TARGET_IDS
    } | {
        "12": {
            "status": "NOT_APPLICABLE_DIFFERENT_FIBRATION",
            "split_cover_count": None,
            "T_rank": None,
            "reason": (
                "curve 12 is a 11952 fibre and has no 074d9 base parameter; "
                "cross-fibration intersection requires an explicit K3 birational transport"
            ),
        }
    }


def screen_cheapest_finite_traces(
    *, cheapest_64, curve, multiples, A, B, discriminant, ring, field,
    targets, lineage, public, chord_helper, cross_helper, quotient_helper,
):
    """Run the target-blind cheapest-sixteen finite-pole trace screen."""

    selected = []
    skipped = []
    for priority_rank, trace_record in enumerate(cheapest_64, start=1):
        trace = cross_helper["trace_from_word"](
            trace_record["word"], curve, multiples
        )
        frame = chord_helper["trace_chord_frame"](trace[0], trace[1], ring)
        if frame["h"].degree() != 2:
            skipped.append(
                {
                    "priority_rank": priority_rank,
                    "orbit_mask": trace_record["orbit_mask"],
                    "trace_denominator_degree": int(frame["h"].degree()),
                }
            )
            continue
        selected.append((priority_rank, trace_record, trace, frame))
        if len(selected) == 16:
            break
    if len(selected) != 16:
        raise ArithmeticError("the cheapest-64 list did not contain sixteen finite-pole traces")

    aggregated = {
        source_id: {target_id: [] for target_id in TARGET_IDS}
        for source_id in SOURCE_IDS
    }
    trace_summaries = []
    for screen_rank, (priority_rank, trace_record, trace, frame) in enumerate(
        selected, start=1
    ):
        source_summaries = []
        for source_id in SOURCE_IDS:
            covers = fit_source_covers(
                source_id=source_id,
                trace_rank=priority_rank,
                trace=trace,
                frame=frame,
                A=A,
                B=B,
                discriminant=discriminant,
                ring=ring,
                field=field,
                targets=targets,
                chord_helper=chord_helper,
            )
            split_records, evaluations = evaluate_frozen_covers(
                covers, trace, A, B, ring, targets, cross_helper
            )
            split_counts = {
                target_id: len(split_records[target_id]) for target_id in TARGET_IDS
            }
            if split_counts[source_id] != 12:
                raise ArithmeticError("a cheapest-sixteen fitted diagonal lost a direction")
            for target_id in TARGET_IDS:
                if target_id != source_id:
                    aggregated[source_id][target_id].extend(split_records[target_id])
            source_summaries.append(
                {
                    "source_curve_id": source_id,
                    "split_counts": {
                        str(target_id): split_counts[target_id]
                        for target_id in TARGET_IDS
                    },
                    "cover_fingerprints": [
                        {
                            "cover_id": cover["cover_id"],
                            "branch_quartic_sha256": hashlib.sha256(
                                json.dumps(
                                    cover["branch_quartic_coefficients_low_to_high"],
                                    separators=(",", ":"),
                                ).encode()
                            ).hexdigest(),
                            "off_diagonal_split_targets": [
                                item["curve_id"]
                                for item in cover["frozen_evaluations"]
                                if item["curve_id"] != source_id
                                and item["status"] == "SPLIT_OVER_Q"
                            ],
                        }
                        for cover in evaluations
                    ],
                }
            )
        trace_summaries.append(
            {
                "screen_rank": screen_rank,
                "global_priority_rank": priority_rank,
                "orbit_mask": trace_record["orbit_mask"],
                "equation_basis_word": list(trace_record["word"]),
                "trace_denominator_h_coefficients_low_to_high": polynomial_text(
                    frame["h"]
                ),
                "source_results": source_summaries,
            }
        )
        print(
            f"{PROTOCOL}|stage=screen16|completed={screen_rank}/16"
            f"|priority_rank={priority_rank}",
            flush=True,
        )

    isomorphisms = {
        curve_id: targets[curve_id]["isomorphism"] for curve_id in TARGET_IDS
    }
    quotient_results = {}
    matrix_rows = {}
    for source_id in SOURCE_IDS:
        quotient_results[source_id] = cross_helper["resolve_quotients"](
            lineage,
            public,
            aggregated[source_id],
            isomorphisms,
            quotient_helper,
        )
        matrix_rows[source_id] = {}
        for target_id in TARGET_IDS:
            if target_id == source_id:
                matrix_rows[source_id][target_id] = {
                    "status": "EXACT_FITTED_DIAGONAL",
                    "split_cover_count": 16 * 12,
                    "T_rank": 12,
                }
            else:
                result = quotient_results[source_id][target_id]
                matrix_rows[source_id][target_id] = {
                    "status": "EXACT",
                    "split_cover_count": result["split_count"],
                    "T_rank": result["class_span_rank"],
                }
        matrix_rows[source_id][12] = {
            "status": "NOT_APPLICABLE_DIFFERENT_FIBRATION",
            "split_cover_count": None,
            "T_rank": None,
        }

    return {
        "selection": {
            "requested_trace_count": 16,
            "selection_rule": (
                "first sixteen finite-pole traces in the complete target-blind "
                "equation-cost ordering"
            ),
            "selected_global_priority_ranks": [row[0] for row in selected],
            "skipped_earlier_nonfinite_traces": skipped,
        },
        "traces": trace_summaries,
        "transfer_matrix": {
            str(source_id): {
                str(target_id): matrix_rows[source_id][target_id]
                for target_id in (*TARGET_IDS, 12)
            }
            for source_id in SOURCE_IDS
        },
        "off_diagonal_quotient_results": {
            str(source_id): [
                {"curve_id": target_id, **quotient_results[source_id][target_id]}
                for target_id in TARGET_IDS
                if target_id != source_id
            ]
            for source_id in SOURCE_IDS
        },
        "record_corridor": {
            "forward_356_to_385": matrix_rows[356][385],
            "reverse_385_to_356": matrix_rows[385][356],
            "at_least_one_record_anchored_cover": bool(
                matrix_rows[356][385]["split_cover_count"]
                or matrix_rows[385][356]["split_cover_count"]
            ),
        },
    }


def build_payload():
    lineage = json.loads(LINEAGE.read_text())
    public = json.loads(PUBLIC.read_text())
    if lineage["representative"]["chart"] != CHART:
        raise ArithmeticError("the lineage representative changed")
    chord_helper = runpy.run_path(str(CHORD_HELPER))
    cross_helper = runpy.run_path(str(CROSS_HELPER))
    quotient_helper = runpy.run_path(str(QUOTIENT_HELPER))
    enumerator = runpy.run_path(str(ENUMERATOR))
    gram = matrix(ZZ, lineage["generic_basis"]["height_gram"])
    selected, cheapest_64, streaming = canonical_norm8_trace(gram, enumerator)
    ring, field, A, B, discriminant, curve, _basis, multiples = cross_helper[
        "build_exact_context"
    ](lineage)
    trace = cross_helper["trace_from_word"](selected["word"], curve, multiples)
    frame = chord_helper["trace_chord_frame"](trace[0], trace[1], ring)
    if frame["h"].degree() != 2:
        raise ArithmeticError("the target-blind first trace does not define a finite norm-eight pencil")
    targets = target_data(lineage, public, cross_helper)

    source_results = {}
    transfer_matrix = {}
    for source_id in SOURCE_IDS:
        covers = fit_source_covers(
            source_id=source_id,
            trace_rank=1,
            trace=trace,
            frame=frame,
            A=A,
            B=B,
            discriminant=discriminant,
            ring=ring,
            field=field,
            targets=targets,
            chord_helper=chord_helper,
        )
        split_records, evaluations = evaluate_frozen_covers(
            covers, trace, A, B, ring, targets, cross_helper
        )
        quotient_results = cross_helper["resolve_quotients"](
            lineage,
            public,
            split_records,
            {curve_id: targets[curve_id]["isomorphism"] for curve_id in TARGET_IDS},
            quotient_helper,
        )
        if (
            quotient_results[source_id]["split_count"] != 12
            or quotient_results[source_id]["class_span_rank"] != 12
        ):
            raise ArithmeticError("the fitted diagonal did not recover all twelve source directions")
        source_results[source_id] = {
            "source_curve_id": source_id,
            "fitted_cover_count": len(covers),
            "covers": evaluations,
            "quotient_transfer_results": [
                {"curve_id": curve_id, **quotient_results[curve_id]}
                for curve_id in TARGET_IDS
            ],
        }
        transfer_matrix[source_id] = transfer_row(source_id, quotient_results)
        print(
            f"{PROTOCOL}|source={source_id}|splits="
            + ",".join(
                f"{target}:{quotient_results[target]['split_count']}"
                for target in TARGET_IDS
            )
            + "|ranks="
            + ",".join(
                f"{target}:{quotient_results[target]['class_span_rank']}"
                for target in TARGET_IDS
            ),
            flush=True,
        )

    forward_record = transfer_matrix[356]["385"]
    reverse_record = transfer_matrix[385]["356"]
    shared_forward = [
        cover["cover_id"]
        for cover in source_results[356]["covers"]
        if next(
            record for record in cover["frozen_evaluations"] if record["curve_id"] == 385
        )["status"]
        == "SPLIT_OVER_Q"
    ]
    shared_reverse = [
        cover["cover_id"]
        for cover in source_results[385]["covers"]
        if next(
            record for record in cover["frozen_evaluations"] if record["curve_id"] == 356
        )["status"]
        == "SPLIT_OVER_Q"
    ]
    cheapest_16_screen = screen_cheapest_finite_traces(
        cheapest_64=cheapest_64,
        curve=curve,
        multiples=multiples,
        A=A,
        B=B,
        discriminant=discriminant,
        ring=ring,
        field=field,
        targets=targets,
        lineage=lineage,
        public=public,
        chord_helper=chord_helper,
        cross_helper=cross_helper,
        quotient_helper=quotient_helper,
    )
    inputs = {
        relative(path): digest(path)
        for path in (
            LINEAGE,
            PUBLIC,
            ENUMERATOR,
            CHORD_HELPER,
            CROSS_HELPER,
            QUOTIENT_HELPER,
        )
    }
    return {
        "schema": "elkies-k3.r17-074d9-norm8-cross-fibre-transfer.v1",
        "status": STATUS,
        "claim": (
            "Exact target-blind selection, source incidence fits, frozen off-diagonal "
            "specializations, and displayed-quotient transfer ranks for one canonical "
            "074d9 norm-eight trace."
        ),
        "claim_boundary": (
            "The diagonal members are deliberately fitted. Off-diagonal rational splits "
            "and quotient coordinates are exact. A nonsquare is decided by exact integer "
            "square tests. Transfers are defined only among fibres of the same 074d9 "
            "fibration. Quotients refer to displayed public subgroups and are not full "
            "Mordell--Weil rank upper bounds."
        ),
        "canonical_trace": {
            "selection_is_target_independent": True,
            "complete_minimum_norm_eight_class_count": 63925,
            "excluded_lower_norm_coset_count": len(
                streaming["masks_by_norm"][2]
                | streaming["masks_by_norm"][4]
                | streaming["masks_by_norm"][6]
            ),
            "selection_score": [
                "coefficient_l1",
                "support_count",
                "maximum_absolute_coefficient",
                "lexicographic_oriented_equation_basis_word",
            ],
            "selected_priority_rank": 1,
            "selected_orbit_mask": selected["orbit_mask"],
            "selected_equation_basis_word": list(selected["word"]),
            "selected_minimum_unoriented_count": selected[
                "minimum_unoriented_count"
            ],
            "trace_denominator_h_coefficients_low_to_high": polynomial_text(frame["h"]),
            "trace_Nx_coefficients_low_to_high": polynomial_text(frame["Nx"]),
            "trace_Ny_coefficients_low_to_high": polynomial_text(frame["Ny"]),
            "least_slope_numerator_M0_coefficients_low_to_high": polynomial_text(
                frame["M0"]
            ),
            "cheapest_64": [
                {
                    "priority_rank": rank,
                    "orbit_mask": record["orbit_mask"],
                    "equation_basis_word": list(record["word"]),
                    "score": list(record["score"][:3]),
                    "minimum_unoriented_count": record[
                        "minimum_unoriented_count"
                    ],
                }
                for rank, record in enumerate(cheapest_64, start=1)
            ],
        },
        "source_experiments": [source_results[source_id] for source_id in SOURCE_IDS],
        "transfer_matrix": {
            str(source_id): transfer_matrix[source_id] for source_id in SOURCE_IDS
        },
        "record_corridor": {
            "forward_356_to_385": forward_record,
            "reverse_385_to_356": reverse_record,
            "forward_shared_cover_ids": shared_forward,
            "reverse_shared_cover_ids": shared_reverse,
            "at_least_one_record_anchored_norm8_cover": bool(
                shared_forward or shared_reverse
            ),
        },
        "different_family_control": {
            "curve_id": 12,
            "status": "NOT_APPLICABLE_DIFFERENT_FIBRATION",
            "required_next_input": (
                "an explicit birational common-K3 transport that turns a frozen 074d9 "
                "cover into intersection data on the 11952 curve-12 fibre"
            ),
        },
        "cheapest_16_finite_trace_screen": cheapest_16_screen,
        "generation": {
            "command": (
                "sage -python elkies-k3/scripts/"
                "certify_r17_074d9_norm8_cross_fibre_transfer.sage"
            ),
            "checker_sha256": digest(Path(__file__)),
            "inputs": inputs,
        },
        "software_assumptions": {
            "sage_version": SAGE_VERSION,
            "required_features": [
                "exact lattice enumeration",
                "exact elliptic function fields and group law",
                "PARI height relations checked by exact addition",
            ],
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
            raise SystemExit("stale 074d9 norm-eight transfer certificate")
        terminal = "PASS"
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered)
        terminal = "WROTE"
    print(
        f"{PROTOCOL}|forward={payload['record_corridor']['forward_356_to_385']['T_rank']}"
        f"|reverse={payload['record_corridor']['reverse_385_to_356']['T_rank']}"
        f"|corridor={payload['record_corridor']['at_least_one_record_anchored_norm8_cover']}"
        f"|status={terminal}|output={arguments.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
