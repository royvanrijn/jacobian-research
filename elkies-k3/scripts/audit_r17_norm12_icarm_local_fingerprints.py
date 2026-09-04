#!/usr/bin/env python3
"""Compute exact local/Nagao fingerprints for all 69 recognized ICARM fibres.

The six exact j-map classes supply the native family models and parameters.
For each pinned fibre this replay records the native local symbol at every
prime in three fixed blocks, checks it against the target curve whenever the
certified QQ-isomorphism has usable reduction, and adds an exact empirical
within-family incidence rank.  Nagao scores remain heuristic features.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "elkies-k3/scripts"))

from search_h92_q12o5867_rootless_nagao import (  # noqa: E402
    Candidate,
    DEFAULT_PRIME_BLOCKS,
    FamilyModel,
    SCORE_SCALE,
    build_residue_tables,
    local_symbol_record,
    projective_index,
    score_block,
)


Q = Fraction
ATLAS = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-record-lineage-atlas-v1.json"
SWEEP = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-database-sweep-v1.json"
PUBLIC = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-public-fibres-v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-local-fingerprints-v1.json"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def short_coefficients(ainvs: list[str]) -> tuple[Q, Q]:
    a1, a2, a3, a4, a6 = map(Q, ainvs)
    b2 = a1 * a1 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3 * a3 + 4 * a6
    c4 = b2 * b2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    return -c4 / 48, -c6 / 864


def reduce_fraction(value: Q, prime: int) -> int:
    if value.denominator % prime == 0:
        raise ZeroDivisionError
    return value.numerator % prime * pow(value.denominator % prime, -1, prime) % prime


def target_trace(ainvs: list[str], prime: int) -> tuple[int | None, int | None]:
    coefficient_a_q, coefficient_b_q = short_coefficients(ainvs)
    try:
        coefficient_a = reduce_fraction(coefficient_a_q, prime)
        coefficient_b = reduce_fraction(coefficient_b_q, prime)
    except ZeroDivisionError:
        return None, None
    if (4 * coefficient_a**3 + 27 * coefficient_b**2) % prime == 0:
        return None, None
    characters = [-1] * prime
    characters[0] = 0
    for value in range(1, prime):
        characters[value * value % prime] = 1
    trace = -sum(
        characters[(x**3 + coefficient_a * x + coefficient_b) % prime]
        for x in range(prime)
    )
    return trace, prime + 1 - trace


def qtext_mean(values: list[int]) -> str:
    value = Q(sum(values), len(values))
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def qtext_median(values: list[int]) -> str:
    ordered = sorted(values)
    middle = len(ordered) // 2
    value = (
        Q(ordered[middle])
        if len(ordered) % 2
        else Q(ordered[middle - 1] + ordered[middle], 2)
    )
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def build() -> dict[str, object]:
    atlas = json.loads(ATLAS.read_text())
    sweep = json.loads(SWEEP.read_text())
    public = json.loads(PUBLIC.read_text())
    if public["status"] != "PASS_PINNED_PUBLIC_POINT_PROJECTION_FOR_69_RECOGNIZED_FIBRES":
        raise ArithmeticError("the pinned public-fibre input is not certified")

    chart_records = {
        record["label"]: record for record in atlas["atlas"]["charts"]
    }
    hit_records = {
        int(record["curve_id"]): record
        for record in sweep["rational_j_hits_and_twists"]
    }
    public_records = {int(record["id"]): record for record in public["records"]}
    representatives = sorted({record["representative"] for record in hit_records.values()})
    if len(representatives) != 6:
        raise ArithmeticError("the six-class partition changed")

    models: dict[str, FamilyModel] = {}
    tables = {}
    rejected = {}
    for representative in representatives:
        chart = chart_records[representative]
        weierstrass = chart["weierstrass_model"]
        model = FamilyModel(
            source=ATLAS.resolve(),
            source_sha256=digest(ATLAS),
            a_coefficients=tuple(Q(value) for value in weierstrass["A_coefficients_low_to_high"]),
            b_coefficients=tuple(Q(value) for value in weierstrass["B_coefficients_low_to_high"]),
            a_degree=8,
            b_degree=12,
            coordinate=f"native parameter of {representative}",
            coefficient_source_keys=(
                "A_coefficients_low_to_high",
                "B_coefficients_low_to_high",
            ),
        )
        if len(model.a_coefficients) != 9 or len(model.b_coefficients) != 13:
            raise ArithmeticError(f"{representative} is not an (8,12) family model")
        models[representative] = model
        tables[representative], rejected[representative] = build_residue_tables(
            model, DEFAULT_PRIME_BLOCKS
        )

    rows: list[dict[str, object]] = []
    for curve_id in sorted(hit_records):
        hit = hit_records[curve_id]
        record = public_records[curve_id]
        representative = hit["representative"]
        parameter_record = hit["representative_parameter"]
        numerator = int(parameter_record["numerator"])
        denominator = int(parameter_record["denominator"])
        native_twist = next(
            twist
            for twist in hit["native_chart_twists"]
            if twist["chart"] == representative
        )["twist"]
        if native_twist["status"] != "QQ_ISOMORPHIC_UNTWISTED":
            raise ArithmeticError(f"curve {curve_id} is not untwisted in its representative chart")
        scale_s = Q(native_twist["qq_isomorphism_scale_s_with_s_squared_q"])

        candidate = Candidate(
            numerator=numerator,
            denominator=denominator,
            height=max(abs(numerator), denominator),
        )
        inverse_cache: dict[tuple[int, int], int | None] = {}
        for block in tables[representative]:
            candidate = score_block(candidate, block, inverse_cache)

        local_records: list[dict[str, object]] = []
        checked = 0
        skipped_scale = 0
        for block_number, block in enumerate(tables[representative], start=1):
            for prime, table in block.items():
                index = projective_index(numerator, denominator, prime)
                symbol = table[index]
                target_ap, target_count = target_trace(record["ainvs"], prime)
                scale_usable = (
                    scale_s.numerator % prime != 0
                    and scale_s.denominator % prime != 0
                )
                if symbol.good_reduction and target_ap is not None and scale_usable:
                    if symbol.trace != target_ap or symbol.point_count != target_count:
                        raise ArithmeticError(
                            f"local trace mismatch for curve {curve_id} at p={prime}"
                        )
                    checked += 1
                elif not scale_usable:
                    skipped_scale += 1

                good_population = [entry for entry in table if entry.good_reduction]
                scores = [entry.contribution_units for entry in good_population]
                row = local_symbol_record(symbol)
                row.update(
                    {
                        "prime": prime,
                        "block": block_number,
                        "target_a_p": target_ap,
                        "target_point_count": target_count,
                        "target_comparison_checked": bool(
                            symbol.good_reduction and target_ap is not None and scale_usable
                        ),
                        "family_good_parameter_count": len(good_population),
                        "family_score_mean_units_1e12": qtext_mean(scores),
                        "family_upper_tail_count_including_ties": (
                            sum(value >= symbol.contribution_units for value in scores)
                            if symbol.good_reduction
                            else None
                        ),
                    }
                )
                local_records.append(row)

        rows.append(
            {
                "curve_id": curve_id,
                "representative": representative,
                "family_frame_class": hit["representative_frame_class"],
                "parameter": parameter_record["affine_parameter"],
                "projective_pair": [numerator, denominator],
                "projective_height": candidate.height,
                "snapshot_rank_lower_bound": int(hit["snapshot_rank_lower_bound"]),
                "public_displayed_point_count": len(record["points"]),
                "bad_prime_count": len(record["bad_primes"]),
                "conductor_decimal_digits": len(str(record["conductor"]).lstrip("-")),
                "block_score_units_1e12": list(candidate.block_score_units),
                "total_score_units_1e12": candidate.total_score_units,
                "good_prime_count": candidate.good_primes,
                "bad_reduction_prime_count": candidate.bad_primes,
                "target_trace_comparisons_checked": checked,
                "comparisons_skipped_for_nonunit_isomorphism_scale": skipped_scale,
                "local_symbols": local_records,
            }
        )

    class_summaries = []
    for representative in representatives:
        group = [row for row in rows if row["representative"] == representative]
        ranks = [int(row["snapshot_rank_lower_bound"]) for row in group]
        scores = [int(row["total_score_units_1e12"]) for row in group]
        class_summaries.append(
            {
                "representative": representative,
                "family_frame_class": group[0]["family_frame_class"],
                "recognized_public_hit_count": len(group),
                "rank_lower_bound_minimum": min(ranks),
                "rank_lower_bound_maximum": max(ranks),
                "rank_lower_bound_mean": qtext_mean(ranks),
                "rank_lower_bound_multiset": sorted(ranks),
                "nagao_score_units_1e12_minimum": min(scores),
                "nagao_score_units_1e12_median": qtext_median(scores),
                "nagao_score_units_1e12_maximum": max(scores),
                "nagao_score_units_1e12_mean": qtext_mean(scores),
                "rejected_model_primes": list(rejected[representative]),
            }
        )

    return {
        "schema": "elkies-k3.r17-norm12-icarm-local-fingerprints.v1",
        "status": "PASS_EXACT_LOCAL_FINGERPRINTS_FOR_ALL_69_RECOGNIZED_FIBRES",
        "summary": {
            "recognized_fibres": len(rows),
            "representative_families": len(representatives),
            "prime_blocks": [list(block) for block in DEFAULT_PRIME_BLOCKS],
            "score_scale": SCORE_SCALE,
            "all_usable_target_trace_comparisons_passed": True,
        },
        "class_summaries": class_summaries,
        "fibres": rows,
        "feature_table_columns": [
            "parameter",
            "representative/family_frame_class",
            "snapshot_rank_lower_bound",
            "recognized_public_hit_count (class-level selection/effort proxy only)",
            "local a_p vector and within-family upper-tail incidences",
            "three Nagao block scores",
        ],
        "claim_boundary": {
            "proved": [
                "the recorded local point counts and Frobenius traces are exact",
                "all locally usable native-family traces equal the corresponding target-curve traces",
                "all 69 exact-sweep hits receive the same predeclared local feature map",
            ],
            "heuristic_only": [
                "Nagao score as a rank-jump predictor",
                "recognized public hit count as a proxy for historical search effort",
                "any family-quality versus search-effort interpretation",
            ],
        },
        "inputs": {
            relative(ATLAS): digest(ATLAS),
            relative(SWEEP): digest(SWEEP),
            relative(PUBLIC): digest(PUBLIC),
        },
        "reproducing_command": (
            ".venv/bin/python elkies-k3/scripts/"
            "audit_r17_norm12_icarm_local_fingerprints.py"
        ),
        "software_assumptions": {"python": sys.version.split()[0]},
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    output = args.output.resolve()
    serialized = json.dumps(build(), indent=2, sort_keys=True) + "\n"
    if args.check:
        if not output.exists() or output.read_text() != serialized:
            raise ArithmeticError("stored local fingerprints differ from exact replay")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
    print(
        "R17ICARMLOCAL|fibres=69|families=6|prime_blocks=3|"
        f"status=PASS|output={relative(output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
