#!/usr/bin/env sage-python
"""Scan the pinned 474-curve ICARM snapshot against the six norm-12 j-classes.

The 43 compiled norm-twelve fibrations form six exact rational-PGL2 classes.
It is therefore enough to solve one projective degree-24 j-preimage equation
per target curve and class.  A rational hit is transported to every native
chart in its class, where the specialized short Weierstrass model is compared
with the target over QQ and its exact quadratic-twist class is recorded.

The original 2026-09-01 database response had SHA-256
``18699517...`` and contained curves 1 through 474.  The live append-only
endpoint has since grown.  For a stable replay we select exactly ids 1..474
and require the hash-pinned projection (id, curve key, a-invariants, creation
time) below.  Rank improvements made after the cutoff are rolled back from
the public history when rank metadata is reported; they do not affect the
equation scan.  The default replay uses the exact 474-row projection already
stored in the certificate; ``--live-source`` is an explicit drift audit.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from math import gcd
from pathlib import Path
import sys
from urllib.request import urlopen

from sage.all import PolynomialRing, QQ, ZZ, lcm
from sage.env import SAGE_VERSION


sys.set_int_max_str_digits(0)

ROOT = Path(__file__).resolve().parents[2]
ATLAS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-record-lineage-atlas-v1.json"
)
LINEAGE = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/icarm_wgxli_rank17_lineage_v1.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-database-sweep-v1.json"
)
DATABASE_URL = "https://elliptic-rank.icarm.cloud/database.json"
PINNED_DATABASE_DATE = "2026-09-01"
PINNED_DATABASE_CUTOFF = "2026-09-01 00:00:00"
PINNED_DATABASE_COUNT = 474
PINNED_DATABASE_SHA256 = (
    "18699517c2969c8c3a250ae612d5caae9fb23c379fe054ba3c7fdf2ec2a83e50"
)
PINNED_EQUATION_PROJECTION_SHA256 = (
    "6c88a323d4e893072fc88613b002ad01f4bb7f6f2f5a059c77383fb92462adfb"
)
EXPECTED_REPRESENTATIVES = (
    "norm12-orbit-074d9",
    "norm12-orbit-08234",
    "norm12-orbit-0e80b",
    "norm12-orbit-11952",
    "norm12-orbit-07ca9",
    "norm12-orbit-08f72",
)
WGXLI_COMPONENTS = (
    (351, 356, 376, 377, 385),
    (363, 364, 378),
    (389, 390, 391),
    (393,),
    (395,),
)
EXPECTED_CLASS_HIT_COUNTS = {
    "norm12-orbit-074d9": 5,
    "norm12-orbit-08234": 54,
    "norm12-orbit-0e80b": 2,
    "norm12-orbit-11952": 2,
    "norm12-orbit-07ca9": 3,
    "norm12-orbit-08f72": 3,
}
RATIONAL_ROOT_OBSTRUCTION_PRIMES = (
    5,
    7,
    11,
    13,
    17,
    19,
    23,
    29,
    31,
    37,
    41,
    43,
    47,
    53,
    59,
    61,
    67,
    71,
    73,
    79,
    83,
    89,
    97,
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


def rational_text(value) -> str:
    value = QQ(value)
    if value.denominator() == 1:
        return str(value.numerator())
    return f"{value.numerator()}/{value.denominator()}"


def polynomial_from_record(record, ring, key):
    return ring([QQ(value) for value in record["weierstrass_model"][key]])


def weierstrass_invariants(ainvs):
    if len(ainvs) == 2:
        ainvs = [0, 0, 0, ainvs[0], ainvs[1]]
    if len(ainvs) != 5:
        raise ArithmeticError("an ICARM record has neither two nor five a-invariants")
    a1, a2, a3, a4, a6 = map(QQ, ainvs)
    b2 = a1**2 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3**2 + 4 * a6
    b8 = a1**2 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3**2 - a4**2
    c4 = b2**2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    delta = -b2**2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
    if not delta or c4**3 - c6**2 != 1728 * delta:
        raise ArithmeticError("an ICARM target is singular or violates the invariant identity")
    return {
        "c4": c4,
        "c6": c6,
        "delta": delta,
        "j": c4**3 / delta,
        "A": -c4 / 48,
        "B": -c6 / 864,
    }


def pinned_rank_lower_bound(record) -> int:
    rank = int(record["rank_lower_bound"])
    later_improvements = sorted(
        (
            item
            for item in record.get("history", [])
            if item.get("kind") == "rank_improved"
            and item.get("at", "") >= PINNED_DATABASE_CUTOFF
        ),
        key=lambda item: item["at"],
    )
    if later_improvements:
        rank = int(later_improvements[0]["old_rank"])
    return rank


def equation_projection(records) -> list[dict]:
    return [
        {
            "id": int(record["id"]),
            "curve_key": record["curve_key"],
            "ainvs": record["ainvs"],
            "created_at": record["created_at"],
        }
        for record in records
    ]


def projection_digest(projection) -> str:
    raw = json.dumps(projection, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def load_pinned_records(
    database_path: Path | None,
    *,
    live_source: bool,
    stored_certificate: Path,
):
    if database_path is not None:
        raw = database_path.read_bytes()
        payload = json.loads(raw)
        source_records = payload["curves"]
    elif live_source:
        with urlopen(DATABASE_URL, timeout=120) as response:
            raw = response.read()
        payload = json.loads(raw)
        source_records = payload["curves"]
    else:
        if not stored_certificate.is_file():
            raise ArithmeticError(
                "offline replay requires the stored certificate; use --database "
                "or --live-source only to reconstruct it"
            )
        stored = json.loads(stored_certificate.read_text())
        if stored.get("schema") != "elkies-k3.r17-norm12-icarm-database-sweep.v1":
            raise ArithmeticError("stored ICARM sweep has an unknown schema")
        stored_snapshot = stored.get("snapshot")
        if not isinstance(stored_snapshot, dict):
            raise ArithmeticError("stored ICARM sweep has no source projection")
        stored_rows = stored_snapshot.get("curves")
        if not isinstance(stored_rows, list):
            raise ArithmeticError("stored ICARM sweep source projection is malformed")
        source_records = [
            {
                "id": record["id"],
                "curve_key": record["curve_key"],
                "ainvs": record["ainvs"],
                "rank_lower_bound": record["snapshot_rank_lower_bound"],
                "submitter": record.get("submitter"),
                "created_at": record["created_at"],
                "history": [],
            }
            for record in stored_rows
        ]
    records_by_id = {int(record["id"]): record for record in source_records}
    expected_ids = set(range(1, PINNED_DATABASE_COUNT + 1))
    if not expected_ids.issubset(records_by_id):
        missing = sorted(expected_ids - set(records_by_id))
        raise ArithmeticError(f"the public source no longer contains pinned ids: {missing}")
    records = [records_by_id[curve_id] for curve_id in sorted(expected_ids)]
    if any(record["created_at"] >= PINNED_DATABASE_CUTOFF for record in records):
        raise ArithmeticError("a pinned curve id has a creation time after the cutoff")
    projection = equation_projection(records)
    observed_projection_digest = projection_digest(projection)
    if observed_projection_digest != PINNED_EQUATION_PROJECTION_SHA256:
        raise ArithmeticError(
            "the pinned ICARM curve-equation projection changed: "
            f"{observed_projection_digest} != {PINNED_EQUATION_PROJECTION_SHA256}"
        )
    return records


def primitive_projective_polynomial(poly, degree: int) -> list[int]:
    coefficients = [QQ(poly[index]) for index in range(degree + 1)]
    denominator_lcm = lcm(value.denominator() for value in coefficients)
    integers = [int(value * denominator_lcm) for value in coefficients]
    content = 0
    for value in integers:
        content = gcd(content, abs(value))
    if not content:
        raise ArithmeticError("cannot normalize the zero projective equation")
    integers = [value // content for value in integers]
    first_nonzero = next(value for value in reversed(integers) if value)
    if first_nonzero < 0:
        integers = [-value for value in integers]
    return integers


def coefficient_digest(coefficients: list[int]) -> str:
    raw = json.dumps(coefficients, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def modular_no_projective_root_prime(primitive: list[int]) -> int | None:
    for prime in RATIONAL_ROOT_OBSTRUCTION_PRIMES:
        residues = [value % prime for value in primitive]
        if not residues[-1]:
            continue
        if all(
            sum(
                coefficient * pow(value, index, prime)
                for index, coefficient in enumerate(residues)
            )
            % prime
            for value in range(prime)
        ):
            return prime
    return None


def normalize_projective(numerator, denominator) -> tuple[int, int]:
    numerator = int(ZZ(numerator))
    denominator = int(ZZ(denominator))
    common = gcd(abs(numerator), abs(denominator))
    if not common:
        raise ArithmeticError("the zero projective point is invalid")
    numerator //= common
    denominator //= common
    if denominator < 0 or (denominator == 0 and numerator < 0):
        numerator = -numerator
        denominator = -denominator
    return numerator, denominator


def rational_projective_roots(poly, degree: int):
    primitive = primitive_projective_polynomial(poly, degree)
    obstruction_prime = modular_no_projective_root_prime(primitive)
    roots = []
    decision_method = "projective modular no-root obstruction"
    if obstruction_prime is None:
        roots = [
            normalize_projective(QQ(root).numerator(), QQ(root).denominator())
            for root, _multiplicity in poly.roots(QQ)
        ]
        if not primitive[-1]:
            roots.append((1, 0))
        roots = sorted(set(roots), key=lambda item: (item[1] == 0, QQ(item[0]) / item[1] if item[1] else 0))
        decision_method = "exact QQ rational-root fallback"
    return {
        "primitive": primitive,
        "primitive_equation_sha256": coefficient_digest(primitive),
        "finite_degree": int(poly.degree()),
        "projective_degree": degree,
        "decision_method": decision_method,
        "modular_no_projective_root_prime": obstruction_prime,
        "roots": roots,
    }


def homogeneous_evaluate(poly, numerator, denominator, degree: int):
    return QQ(
        sum(
            poly[index] * numerator**index * denominator ** (degree - index)
            for index in range(degree + 1)
        )
    )


def exact_rational_nth_root(value, exponent: int):
    value = QQ(value)
    if not value:
        return QQ(0)
    if value < 0 and exponent % 2 == 0:
        return None
    numerator_root, numerator_exact = ZZ(value.numerator()).nth_root(
        exponent, truncate_mode=True
    )
    denominator_root, denominator_exact = ZZ(value.denominator()).nth_root(
        exponent, truncate_mode=True
    )
    if not numerator_exact or not denominator_exact:
        return None
    return QQ(numerator_root) / QQ(denominator_root)


def twist_record(source_A, source_B, target):
    source_A = QQ(source_A)
    source_B = QQ(source_B)
    target_A = target["A"]
    target_B = target["B"]
    source_delta = -16 * (4 * source_A**3 + 27 * source_B**2)
    if not source_delta:
        raise ArithmeticError("a rational j-preimage specialized to a singular fibre")

    if target["j"] == 0:
        if source_A or target_A or not source_B or not target_B:
            raise ArithmeticError("the j=0 specialization has inconsistent short invariants")
        ratio = target_B / source_B
        q_value = exact_rational_nth_root(ratio, 3)
        if q_value is None:
            return {
                "status": "SAME_J_NON_QUADRATIC_SEXTIC_TWIST",
                "B_target_over_B_fibre": rational_text(ratio),
            }
    elif target["j"] == 1728:
        if source_B or target_B or not source_A or not target_A:
            raise ArithmeticError("the j=1728 specialization has inconsistent short invariants")
        ratio = target_A / source_A
        q_value = exact_rational_nth_root(ratio, 2)
        if q_value is None:
            return {
                "status": "SAME_J_NON_QUADRATIC_QUARTIC_TWIST",
                "A_target_over_A_fibre": rational_text(ratio),
            }
    else:
        if not source_A or not source_B or not target_A or not target_B:
            raise ArithmeticError("a generic-j specialization has a zero short coefficient")
        q_value = target_B * source_A / (source_B * target_A)

    if target_A != q_value**2 * source_A or target_B != q_value**3 * source_B:
        raise ArithmeticError("the computed quadratic-twist parameter fails the model identities")
    result = {
        "status": (
            "QQ_ISOMORPHIC_UNTWISTED"
            if q_value.is_square()
            else "NONTRIVIAL_QUADRATIC_TWIST"
        ),
        "quadratic_twist_parameter_q": rational_text(q_value),
        "model_identities": "A_target=q^2*A_fibre; B_target=q^3*B_fibre",
    }
    if q_value.is_square():
        result["qq_isomorphism_scale_s_with_s_squared_q"] = rational_text(
            q_value.sqrt()
        )
        result["qq_isomorphism"] = "x_target=q*x_fibre; y_target=s^3*y_fibre"
    return result


def transport_projective(point, matrix_entries):
    numerator, denominator = point
    a, b, c, d = map(ZZ, matrix_entries)
    return normalize_projective(
        a * numerator + b * denominator,
        c * numerator + d * denominator,
    )


def projective_parameter_record(point):
    numerator, denominator = point
    return {
        "numerator": numerator,
        "denominator": denominator,
        "affine_parameter": (
            rational_text(QQ(numerator) / denominator) if denominator else "infinity"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, default=ATLAS)
    parser.add_argument("--database", type=Path)
    parser.add_argument(
        "--live-source",
        action="store_true",
        help="recover ids 1..474 from the current append-only endpoint",
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.database is not None and args.live_source:
        parser.error("--database and --live-source are mutually exclusive")
    args.atlas = args.atlas.resolve()
    args.output = args.output.resolve()

    atlas = json.loads(args.atlas.read_text())
    if atlas["atlas"]["chart_count"] != 43:
        raise ArithmeticError("the source atlas no longer contains exactly 43 charts")
    classes = atlas["atlas"]["pgl2_equivalence_classes"]
    representatives = tuple(record["representative"] for record in classes)
    if representatives != EXPECTED_REPRESENTATIVES:
        raise ArithmeticError("the six rational-PGL2 representatives changed")
    if [len(record["members"]) for record in classes] != [8, 5, 18, 6, 2, 4]:
        raise ArithmeticError("the six rational-PGL2 class sizes changed")

    lineage = json.loads(LINEAGE.read_text())
    lineage_snapshot = lineage["snapshot"]
    if (
        lineage_snapshot["database_curve_count"] != PINNED_DATABASE_COUNT
        or lineage_snapshot["database_sha256"] != PINNED_DATABASE_SHA256
        or lineage_snapshot["date"] != PINNED_DATABASE_DATE
    ):
        raise ArithmeticError("the original pinned ICARM snapshot provenance changed")

    records = load_pinned_records(
        args.database.resolve() if args.database else None,
        live_source=args.live_source,
        stored_certificate=args.output,
    )
    records_by_id = {int(record["id"]): record for record in records}
    target_invariants = {
        curve_id: weierstrass_invariants(record["ainvs"])
        for curve_id, record in records_by_id.items()
    }
    distinct_j_count = len({record["j"] for record in target_invariants.values()})

    ring = PolynomialRing(QQ, "u")
    charts_by_label = {
        record["label"]: record for record in atlas["atlas"]["charts"]
    }
    chart_models = {
        label: (
            polynomial_from_record(
                record, ring, "A_coefficients_low_to_high"
            ),
            polynomial_from_record(
                record, ring, "B_coefficients_low_to_high"
            ),
        )
        for label, record in charts_by_label.items()
    }
    representative_maps = {}
    for label in representatives:
        A, B = chart_models[label]
        representative_maps[label] = {
            "A": A,
            "B": B,
            "c4_cubed": (-48 * A) ** 3,
            "delta": -16 * (4 * A**3 + 27 * B**2),
        }

    decision_cache = {}
    decision_ledger = []
    rational_hits = []
    obstruction_counts = Counter()
    fallback_count = 0
    miss_count = 0
    class_hit_count = 0
    chart_twist_counts = Counter()

    for record_index, record in enumerate(records, start=1):
        curve_id = int(record["id"])
        target = target_invariants[curve_id]
        target_j = target["j"]
        curve_decisions = []
        for class_record in classes:
            representative = class_record["representative"]
            cache_key = (representative, target_j)
            decision = decision_cache.get(cache_key)
            if decision is None:
                source_map = representative_maps[representative]
                equation = ring(
                    target_j.denominator() * source_map["c4_cubed"]
                    - target_j.numerator() * source_map["delta"]
                )
                if not equation:
                    raise ArithmeticError("a nonconstant j-map became a constant target map")
                decision = rational_projective_roots(equation, 24)
                decision_cache[cache_key] = decision

            roots = decision["roots"]
            compact_decision = {
                "representative": representative,
                "primitive_equation_sha256": decision[
                    "primitive_equation_sha256"
                ],
                "finite_degree": decision["finite_degree"],
                "decision_method": decision["decision_method"],
                "modular_no_projective_root_prime": decision[
                    "modular_no_projective_root_prime"
                ],
                "rational_projective_root_count": len(roots),
            }
            if roots:
                compact_decision["rational_projective_roots"] = [
                    projective_parameter_record(point) for point in roots
                ]
            curve_decisions.append(compact_decision)

            if decision["modular_no_projective_root_prime"] is not None:
                obstruction_counts[decision["modular_no_projective_root_prime"]] += 1
            else:
                fallback_count += 1
            if not roots:
                miss_count += 1
                continue

            class_hit_count += 1
            for root in roots:
                chart_twists = []
                for member in class_record["members"]:
                    label = member["label"]
                    member_point = transport_projective(
                        root,
                        member[
                            "representative_to_member_pgl2_matrix_a_b_c_d"
                        ],
                    )
                    A_member, B_member = chart_models[label]
                    source_A = homogeneous_evaluate(
                        A_member, member_point[0], member_point[1], 8
                    )
                    source_B = homogeneous_evaluate(
                        B_member, member_point[0], member_point[1], 12
                    )
                    twist = twist_record(source_A, source_B, target)
                    chart_twist_counts[twist["status"]] += 1
                    chart_twists.append(
                        {
                            "chart": label,
                            "frame_class": charts_by_label[label]["frame_class"],
                            "native_parameter": projective_parameter_record(
                                member_point
                            ),
                            "twist": twist,
                        }
                    )
                rational_hits.append(
                    {
                        "curve_id": curve_id,
                        "snapshot_rank_lower_bound": pinned_rank_lower_bound(record),
                        "submitter": record.get("submitter"),
                        "representative": representative,
                        "representative_frame_class": charts_by_label[representative][
                            "frame_class"
                        ],
                        "representative_parameter": projective_parameter_record(root),
                        "primitive_preimage_equation_coefficients_u_degree_0_through_24": [
                            str(value) for value in decision["primitive"]
                        ],
                        "native_chart_twists": chart_twists,
                    }
                )

        decision_ledger.append(
            {
                "curve_id": curve_id,
                "snapshot_rank_lower_bound": pinned_rank_lower_bound(record),
                "classes": curve_decisions,
            }
        )
        if record_index % 50 == 0 or record_index == len(records):
            print(
                f"progress={record_index}/{len(records)}|class_hits={class_hit_count}",
                flush=True,
            )

    if miss_count + class_hit_count != PINNED_DATABASE_COUNT * len(classes):
        raise ArithmeticError("the class decision count is inconsistent")

    hit_ids = sorted({record["curve_id"] for record in rational_hits})
    class_hit_counts = Counter(
        record["representative"] for record in rational_hits
    )
    if dict(class_hit_counts) != EXPECTED_CLASS_HIT_COUNTS:
        raise ArithmeticError("the pinned class-hit distribution changed")
    if len(rational_hits) != class_hit_count or len(hit_ids) != class_hit_count:
        raise ArithmeticError(
            "a pinned curve acquired multiple rational roots or multiple hit classes"
        )
    if {273, 302, 398}.intersection(hit_ids):
        raise ArithmeticError("a pinned rank-30 or rank-31 curve acquired a hit")
    if chart_twist_counts != Counter({"QQ_ISOMORPHIC_UNTWISTED": 376}):
        raise ArithmeticError("the pinned native-chart twist distribution changed")
    wgxli_component_results = []
    for component in WGXLI_COMPONENTS:
        component_hits = [curve_id for curve_id in component if curve_id in hit_ids]
        wgxli_component_results.append(
            {
                "curve_ids": list(component),
                "rational_j_hit_curve_ids": component_hits,
                "all_miss": not component_hits,
            }
        )
    if any(
        record["rational_j_hit_curve_ids"] != record["curve_ids"]
        for record in wgxli_component_results
    ):
        raise ArithmeticError("a pinned wgxli component is not completely recognized")

    snapshot_curves = [
        {
            "id": int(record["id"]),
            "curve_key": record["curve_key"],
            "ainvs": record["ainvs"],
            "snapshot_rank_lower_bound": pinned_rank_lower_bound(record),
            "submitter": record.get("submitter"),
            "created_at": record["created_at"],
        }
        for record in records
    ]
    payload = {
        "schema": "elkies-k3.r17-norm12-icarm-database-sweep.v1",
        "status": "PASS_EXACT_COMPLETE_PINNED_ICARM_J_PREIMAGE_AND_TWIST_SWEEP",
        "outcome": {
            "pinned_curve_count": PINNED_DATABASE_COUNT,
            "distinct_target_j_invariant_count": distinct_j_count,
            "rational_pgl2_j_map_class_count": len(classes),
            "class_preimage_decision_count": PINNED_DATABASE_COUNT * len(classes),
            "class_rational_hit_count": class_hit_count,
            "class_hit_counts_by_representative": dict(class_hit_counts),
            "class_miss_count": miss_count,
            "hit_curve_ids": hit_ids,
            "rational_hit_root_count": len(rational_hits),
            "native_chart_twist_counts": dict(sorted(chart_twist_counts.items())),
            "wgxli_components": wgxli_component_results,
        },
        "snapshot": {
            "original_database_date": PINNED_DATABASE_DATE,
            "original_raw_database_sha256": PINNED_DATABASE_SHA256,
            "original_curve_count": PINNED_DATABASE_COUNT,
            "original_curve_ids": "1..474",
            "equation_projection_fields": [
                "id",
                "curve_key",
                "ainvs",
                "created_at",
            ],
            "equation_projection_sha256": PINNED_EQUATION_PROJECTION_SHA256,
            "recovery_source": DATABASE_URL,
            "recovery_rule": (
                "select ids 1 through 474, require creation before 2026-09-01, "
                "then require the pinned equation-projection SHA-256"
            ),
            "curves": snapshot_curves,
        },
        "j_map_classes": [
            {
                "representative": record["representative"],
                "size": len(record["members"]),
                "frame_class": charts_by_label[record["representative"]][
                    "frame_class"
                ],
                "members": [member["label"] for member in record["members"]],
            }
            for record in classes
        ],
        "decision_summary": {
            "modular_obstruction_prime_counts": {
                str(prime): count for prime, count in sorted(obstruction_counts.items())
            },
            "exact_QQ_fallback_count": fallback_count,
            "exact_QQ_fallback_miss_count": sum(
                1
                for record in decision_ledger
                for decision in record["classes"]
                if decision["decision_method"] == "exact QQ rational-root fallback"
                and not decision["rational_projective_root_count"]
            ),
        },
        "rational_j_hits_and_twists": rational_hits,
        "decision_ledger": decision_ledger,
        "method": {
            "preimage_test": (
                "For each target and class representative form the primitive projective "
                "degree-24 equation denominator(j_target)*c4(u)^3-"
                "numerator(j_target)*Delta(u)=0. A prime with no P1(F_p) root "
                "proves a QQ miss; otherwise factor exactly over QQ."
            ),
            "class_transport": (
                "Transport every representative projective root through each stored exact "
                "representative-to-member PGL2(Q) matrix."
            ),
            "twist_test": (
                "Compare short models by A_target=q^2*A_fibre and "
                "B_target=q^3*B_fibre. A rational square q gives a QQ-isomorphism; "
                "a nonsquare q is a nontrivial quadratic twist. At j=0 or 1728, "
                "first test whether the sextic or quartic twist is quadratic."
            ),
        },
        "claim_boundary": {
            "proved": [
                "all 474 curve equations in the pinned 2026-09-01 ICARM snapshot are decided against all six rational-PGL2 j-map classes",
                "every rational j-preimage found by the sweep has an exact native-chart twist classification",
                "every declared miss has either a projective finite-field no-root witness or an exact QQ rational-root factorization",
            ],
            "not_proved": [
                "absence from rootless fibrations outside the certified 43-chart norm-twelve atlas",
                "any Mordell-Weil rank upper bound",
                "that later ICARM curves with ids above 474 miss these six classes",
            ],
        },
        "inputs": {
            relative(args.atlas): digest(args.atlas),
            relative(LINEAGE): digest(LINEAGE),
            "ICARM_database_original_raw_sha256": PINNED_DATABASE_SHA256,
            "ICARM_ids_1_through_474_equation_projection_sha256": PINNED_EQUATION_PROJECTION_SHA256,
        },
        "software_assumptions": {
            "sage_version": SAGE_VERSION,
            "required_features": [
                "exact QQ polynomial factorization",
                "exact modular projective-root evaluation",
                "exact rational square and nth-power tests",
            ],
        },
        "reproducing_command": (
            "sage -python elkies-k3/scripts/certify_r17_norm12_icarm_database_sweep.sage"
        ),
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != serialized:
            raise ArithmeticError("stored complete ICARM sweep differs from exact replay")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    print(
        "R17ICARMSWEEP|curves={}|j_classes={}|decisions={}|hit_curves={}|"
        "class_hits={}|native_twists={}|status=PROVED|output={}".format(
            PINNED_DATABASE_COUNT,
            len(classes),
            PINNED_DATABASE_COUNT * len(classes),
            ",".join(map(str, hit_ids)) if hit_ids else "none",
            class_hit_count,
            ",".join(f"{key}:{value}" for key, value in sorted(chart_twist_counts.items())),
            relative(args.output),
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
