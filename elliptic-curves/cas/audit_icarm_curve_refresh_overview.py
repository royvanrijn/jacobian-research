#!/usr/bin/env python3
"""Freeze and audit the ICARM ids-475-through-573 refresh tranche.

The exact 43-chart j-preimage decisions come from the separately certified
573-curve norm-twelve sweep. This companion artifact retains the public rows
needed for offline follow-up, checks every displayed point and discriminant,
and independently proves the displayed rank lower bounds for the
rank-at-least-24 priority tranche by finite good-reduction quotients.

Construction labels are deterministic commentary tags, not mathematical
provenance theorems. Atlas misses are only misses from the six certified
norm-twelve j-map classes.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
import sys


sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[2]
SWEEP = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-database-sweep-v2.json"
OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/icarm_curve_refresh_475_573_overview_v1.json"
MOD3_SOURCE = ROOT / "elliptic-curves/cas/search_mestre_root_tuple_scale_max200.py"
DATABASE_URL = "https://elliptic-rank.icarm.cloud/database.json"
DATABASE_RAW_SHA256 = "e57d991894722f0e5ab2f548b77f09064a46ec926c93ef3730f47685e016aab0"
SWEEP_SHA256 = "77a3c051111e7ead5ee2a6f88df4a975c2f5bdb87be1bfe4d88b195f293da50c"
MOD3_SOURCE_SHA256 = "405a2b9f7653c89af0e3e6caf2e77765cb4bfc88fccf88edffa67d3435aebf24"
SNAPSHOT_DATE = "2026-09-04"
SNAPSHOT_COUNT = 573
FIRST_NEW_ID = 475
LAST_NEW_ID = 573
PRIORITY_RANK_LOWER_BOUND = 24

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mod2_reduction_independence import (  # noqa: E402
    combined_mod2_rank,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
    finite_curve_points,
)
from search_mestre_root_tuple_scale_max200 import (  # noqa: E402
    mod3_independence_certificate,
)


Q = Fraction


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_digest(value) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return sha256(raw).hexdigest()


def normalize_ainvs(ainvs) -> tuple[Fraction, ...]:
    if len(ainvs) == 2:
        ainvs = [0, 0, 0, *ainvs]
    if len(ainvs) != 5:
        raise ArithmeticError("an ICARM curve has neither two nor five a-invariants")
    return tuple(Q(value) for value in ainvs)


def invariants(ainvs):
    a1, a2, a3, a4, a6 = normalize_ainvs(ainvs)
    b2 = a1 * a1 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3 * a3 + 4 * a6
    b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
    c4 = b2 * b2 - 24 * b4
    c6 = -(b2**3) + 36 * b2 * b4 - 216 * b6
    delta = -(b2 * b2 * b8) - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
    if not delta or c4**3 - c6**2 != 1728 * delta:
        raise ArithmeticError("a refreshed ICARM curve is singular or inconsistent")
    return {
        "a1": a1,
        "a2": a2,
        "a3": a3,
        "a4": a4,
        "a6": a6,
        "b2": b2,
        "delta": delta,
        "short_A": -c4 / 48,
        "short_B": -c6 / 864,
    }


def point_on_general_curve(inv, point) -> bool:
    x_value, y_value = map(Q, point)
    return (
        y_value**2 + inv["a1"] * x_value * y_value + inv["a3"] * y_value
        == x_value**3 + inv["a2"] * x_value**2 + inv["a4"] * x_value + inv["a6"]
    )


def short_point(inv, point):
    x_value, y_value = map(Q, point)
    return (
        x_value + inv["b2"] / 12,
        y_value + (inv["a1"] * x_value + inv["a3"]) / 2,
    )


def primes_up_to(bound: int):
    primes = []
    for candidate in range(2, bound + 1):
        if all(candidate % prime for prime in primes if prime * prime <= candidate):
            primes.append(candidate)
    return primes


def reduce_fraction(value: Fraction, prime: int) -> int:
    value = Q(value)
    if value.denominator % prime == 0:
        raise ValueError("nonintegral coefficient at the selected prime")
    return value.numerator * pow(value.denominator, -1, prime) % prime


def torsion_gcd_certificate(coefficients):
    running_gcd = 0
    records = []
    for prime in primes_up_to(500):
        if prime <= 3:
            continue
        try:
            coefficient_a = reduce_fraction(coefficients[3], prime)
            coefficient_b = reduce_fraction(coefficients[4], prime)
        except ValueError:
            continue
        if (-16 * (4 * coefficient_a**3 + 27 * coefficient_b**2)) % prime == 0:
            continue
        order = len(finite_curve_points(coefficient_a, coefficient_b, prime))
        records.append({"prime": prime, "group_order": order})
        running_gcd = gcd(running_gcd, order)
        if running_gcd == 1:
            return {"gcd": 1, "good_reduction_group_orders": records}
    raise ArithmeticError("the priority curve torsion gcd did not close through 500")


def signature_record(signature):
    return {
        "prime": signature.prime,
        "group_order": signature.group_order,
        "doubled_subgroup_order": signature.doubled_subgroup_order,
        "quotient_dimension": signature.quotient_dimension,
        "rows": [list(row) for row in signature.rows],
    }


def independence_certificate(record, inv):
    points = [short_point(inv, point) for point in record["points"]]
    coefficients = (Q(0), Q(0), Q(0), inv["short_A"], inv["short_B"])
    signatures = find_mod2_reduction_certificate(coefficients, points, prime_bound=500)
    rank = combined_mod2_rank(signatures, len(points))
    no_two_torsion_prime = find_two_torsion_certificate_prime(coefficients, prime_bound=500)
    mod3 = None
    certified_rank = rank
    if rank != len(points):
        mod3 = mod3_independence_certificate(
            coefficients, tuple(points), prime_bound=499
        )
        certified_rank = int(mod3["combined_exact_rank_over_F3"])
    if certified_rank != len(points):
        raise ArithmeticError(
            f"curve {record['id']} priority independence did not close mod 2 or 3"
        )
    return {
        "status": "PROVED_DISPLAYED_POINT_INDEPENDENCE",
        "displayed_point_count": len(points),
        "proved_displayed_subgroup_rank": certified_rank,
        "no_rational_2_torsion_prime": no_two_torsion_prime,
        "mod2_reduction_certificate": {
            "prime_bound": 500,
            "combined_exact_rank_over_F2": rank,
            "signatures": [signature_record(signature) for signature in signatures],
        },
        "mod3_fallback_certificate": mod3,
        "torsion_gcd_certificate": torsion_gcd_certificate(coefficients),
        "proof": (
            "Full F2 rank in a product of E(F_p)/2E(F_p), together with an "
            "irreducible mod-p 2-division cubic, proves Z-independence by infinite "
            "descent. Coprime good-reduction group orders prove trivial rational torsion."
        ),
    }


def construction_tag(commentary) -> str:
    text = (commentary or "").lower()
    if "2-neighbor fibration" in text:
        return "DECLARED_ELKIES_R17_TWO_NEIGHBOR"
    if (
        "rank-17 elliptic k3" in text
        or "rank-17 k3 fibration" in text
        or "rank-17 elliptic fibration" in text
    ):
        return "DECLARED_ELKIES_PUBLISHED_R17"
    if "z/2z rank-9 family" in text:
        return "DECLARED_ELKIES_KLAGSBRUN_Z2_RANK9"
    if "mestre" in text or "fermigier" in text or "quartic family" in text:
        return "DECLARED_MESTRE_FERMIGIER"
    if text.strip():
        return "OTHER_PUBLIC_COMMENTARY"
    return "NO_PUBLIC_CONSTRUCTION_COMMENTARY"


SOURCE_FIELDS = (
    "id",
    "curve_key",
    "ainvs",
    "rank_lower_bound",
    "torsion",
    "naive_height",
    "faltings_height",
    "conductor",
    "bad_primes",
    "discriminant",
    "regulator",
    "points",
    "submitter",
    "commentary",
    "created_at",
    "updated_at",
)


def source_projection(record):
    return {key: record.get(key) for key in SOURCE_FIELDS}


def load_records(database: Path | None, output: Path):
    if database is not None:
        raw = database.read_bytes()
        if sha256(raw).hexdigest() != DATABASE_RAW_SHA256:
            raise ArithmeticError("the supplied 573-curve database response changed")
        payload = json.loads(raw)
        if int(payload["count"]) != SNAPSHOT_COUNT:
            raise ArithmeticError("the supplied database count is not 573")
        all_records = {int(record["id"]): record for record in payload["curves"]}
        if set(all_records) != set(range(1, SNAPSHOT_COUNT + 1)):
            raise ArithmeticError("the supplied database ids are not exactly 1 through 573")
        return [
            source_projection(all_records[curve_id])
            for curve_id in range(FIRST_NEW_ID, LAST_NEW_ID + 1)
        ]
    if not output.is_file():
        raise ArithmeticError("offline replay requires the stored overview or --database")
    stored = json.loads(output.read_text())
    if stored.get("schema") != "elliptic-curves.icarm-refresh-overview.v1":
        raise ArithmeticError("the stored overview schema changed")
    return stored["snapshot"]["records"]


def build(records):
    if [int(record["id"]) for record in records] != list(range(FIRST_NEW_ID, LAST_NEW_ID + 1)):
        raise ArithmeticError("the refresh tranche is not the exact 475..573 id interval")
    if digest(SWEEP) != SWEEP_SHA256:
        raise ArithmeticError("the exact 573-curve atlas sweep changed")
    if digest(MOD3_SOURCE) != MOD3_SOURCE_SHA256:
        raise ArithmeticError("the exact mod-3 fallback implementation changed")
    sweep = json.loads(SWEEP.read_text())
    if sweep["outcome"]["pinned_curve_count"] != SNAPSHOT_COUNT:
        raise ArithmeticError("the atlas sweep is not the 573-curve refresh")
    hits = {
        int(record["curve_id"]): record
        for record in sweep["rational_j_hits_and_twists"]
        if int(record["curve_id"]) >= FIRST_NEW_ID
    }
    ledger = {int(record["curve_id"]): record for record in sweep["decision_ledger"]}

    rank_distribution = Counter()
    construction_distribution = Counter()
    hit_distribution = Counter()
    overview = []
    priority = []
    conductor_missing = []
    for record in records:
        curve_id = int(record["id"])
        rank = int(record["rank_lower_bound"])
        if len(record["points"]) != rank:
            raise ArithmeticError(f"curve {curve_id} point count differs from its lower bound")
        inv = invariants(record["ainvs"])
        if any(not point_on_general_curve(inv, point) for point in record["points"]):
            raise ArithmeticError(f"curve {curve_id} has an off-curve displayed point")
        if inv["delta"] != Q(record["discriminant"]):
            raise ArithmeticError(f"curve {curve_id} stored discriminant changed")
        tag = construction_tag(record.get("commentary"))
        rank_distribution[rank] += 1
        construction_distribution[tag] += 1
        if record.get("conductor") is None:
            conductor_missing.append(curve_id)

        hit = hits.get(curve_id)
        if hit is None:
            decisions = ledger[curve_id]["classes"]
            if any(decision["rational_projective_root_count"] for decision in decisions):
                raise ArithmeticError(f"curve {curve_id} lost a recorded atlas hit")
            atlas = {
                "status": "EXACT_MISS_FROM_ALL_SIX_CERTIFIED_NORM12_J_CLASSES",
                "miss_certificates": [
                    {
                        "representative": decision["representative"],
                        "method": decision["decision_method"],
                        "modular_no_projective_root_prime": decision[
                            "modular_no_projective_root_prime"
                        ],
                        "primitive_equation_sha256": decision["primitive_equation_sha256"],
                    }
                    for decision in decisions
                ],
            }
        else:
            twist_statuses = Counter(
                item["twist"]["status"] for item in hit["native_chart_twists"]
            )
            if set(twist_statuses) != {"QQ_ISOMORPHIC_UNTWISTED"}:
                raise ArithmeticError(f"curve {curve_id} has a nontrivial chart twist")
            hit_distribution[hit["representative"]] += 1
            atlas = {
                "status": "EXACT_UNTWISTED_NORM12_ATLAS_FIBRE",
                "representative": hit["representative"],
                "frame_class": hit["representative_frame_class"],
                "representative_parameter": hit["representative_parameter"],
                "native_chart_count": len(hit["native_chart_twists"]),
                "native_twist_status_counts": dict(twist_statuses),
            }

        item = {
            "curve_id": curve_id,
            "rank_lower_bound": rank,
            "displayed_point_count": len(record["points"]),
            "all_displayed_points_on_curve_exactly": True,
            "stored_discriminant_matches_exact_invariant": True,
            "public_torsion_field": record.get("torsion"),
            "conductor_recorded": record.get("conductor") is not None,
            "construction_commentary_tag": tag,
            "construction_tag_is_only_a_text_classification": True,
            "atlas": atlas,
        }
        if rank >= PRIORITY_RANK_LOWER_BOUND:
            item["independence_and_torsion_certificate"] = independence_certificate(record, inv)
            priority.append(item)
        overview.append(item)

    priority.sort(key=lambda item: (-item["rank_lower_bound"], item["curve_id"]))
    atlas_priority = [
        item["curve_id"]
        for item in priority
        if item["atlas"]["status"] == "EXACT_UNTWISTED_NORM12_ATLAS_FIBRE"
    ]
    non_atlas_priority = [
        item["curve_id"]
        for item in priority
        if item["atlas"]["status"] != "EXACT_UNTWISTED_NORM12_ATLAS_FIBRE"
    ]
    independently_reproved = [
        item["curve_id"]
        for item in priority
        if item["independence_and_torsion_certificate"]["status"]
        == "PROVED_DISPLAYED_POINT_INDEPENDENCE"
    ]
    independence_open = [
        item["curve_id"]
        for item in priority
        if item["curve_id"] not in independently_reproved
    ]
    return {
        "schema": "elliptic-curves.icarm-refresh-overview.v1",
        "status": "PASS_EXACT_OVERVIEW_OF_ICARM_CURVES_475_THROUGH_573",
        "summary": {
            "previous_snapshot_curve_count": 474,
            "refreshed_snapshot_curve_count": SNAPSHOT_COUNT,
            "new_curve_count": len(records),
            "new_curve_id_interval": [FIRST_NEW_ID, LAST_NEW_ID],
            "rank_lower_bound_distribution": {
                str(rank): count for rank, count in sorted(rank_distribution.items(), reverse=True)
            },
            "all_displayed_points_checked_on_curve": True,
            "all_stored_discriminants_checked": True,
            "new_norm12_atlas_hit_count": len(hits),
            "new_norm12_atlas_miss_count": len(records) - len(hits),
            "new_atlas_hits_by_representative": dict(sorted(hit_distribution.items())),
            "priority_rule": f"rank lower bound at least {PRIORITY_RANK_LOWER_BOUND}",
            "priority_curve_ids": [item["curve_id"] for item in priority],
            "priority_atlas_fibre_ids": atlas_priority,
            "priority_non_atlas_ids": non_atlas_priority,
            "priority_rank_lower_bounds_independently_reproved_curve_ids": independently_reproved,
            "priority_curves_requiring_mod3_fallback": [542],
            "priority_independence_not_closed_curve_ids": independence_open,
            "all_priority_curves_have_trivial_torsion_certificates": True,
            "curves_missing_public_conductor": conductor_missing,
            "construction_commentary_tag_counts": dict(sorted(construction_distribution.items())),
        },
        "priority_curves": priority,
        "curve_overview": overview,
        "snapshot": {
            "date": SNAPSHOT_DATE,
            "source": DATABASE_URL,
            "raw_database_sha256": DATABASE_RAW_SHA256,
            "record_projection_fields": list(SOURCE_FIELDS),
            "new_records_projection_sha256": canonical_digest(records),
            "records": records,
        },
        "claim_boundary": {
            "proved": [
                "all 99 newly appended equations and displayed points are internally exact",
                "all 99 stored discriminants match the displayed Weierstrass equations",
                "each of the 99 curves has an exact hit or miss against all six certified norm-twelve j-map classes",
                "the displayed public points on every rank-at-least-24 priority curve are independent by exact finite reductions modulo 2 or 3",
                "every rank-at-least-24 priority curve has trivial rational torsion",
            ],
            "not_proved": [
                "that an atlas miss is absent from every rank-17 or rootless K3 fibration",
                "that any displayed subgroup is the full Mordell-Weil group",
                "an exact rank upper bound for any refreshed curve",
                "that a construction commentary tag independently verifies the stated provenance",
                "independence of the non-priority displayed point lists by this local replay",
            ],
        },
        "inputs": {
            relative(SWEEP): SWEEP_SHA256,
            relative(MOD3_SOURCE): MOD3_SOURCE_SHA256,
            "ICARM_database_raw_sha256": DATABASE_RAW_SHA256,
        },
        "reproducing_command": ".venv/bin/python elliptic-curves/cas/audit_icarm_curve_refresh_overview.py",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    output = args.output.resolve()
    records = load_records(args.database.resolve() if args.database else None, output)
    serialized = json.dumps(build(records), indent=2, sort_keys=True) + "\n"
    if args.check:
        if not output.is_file() or output.read_text() != serialized:
            raise ArithmeticError("the stored ICARM refresh overview changed")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
    summary = json.loads(serialized)["summary"]
    print(
        "ICARMREFRESH|new_curves={}|atlas_hits={}|priority={}|"
        "priority_atlas={}|priority_nonatlas={}|status=PROVED|output={}".format(
            summary["new_curve_count"],
            summary["new_norm12_atlas_hit_count"],
            ",".join(map(str, summary["priority_curve_ids"])),
            ",".join(map(str, summary["priority_atlas_fibre_ids"])),
            ",".join(map(str, summary["priority_non_atlas_ids"])),
            relative(output),
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
