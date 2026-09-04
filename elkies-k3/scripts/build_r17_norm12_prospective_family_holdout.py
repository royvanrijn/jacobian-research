#!/usr/bin/env python3
"""Freeze an ordinary-parameter R17 cohort with whole-PGL2-family holdouts.

This command is intentionally label-free.  It samples rational parameters by
counter-mode SHA-256 draws from a new projective-height shell, commits the six
exact PGL2 family identities, and leaves every feature and search outcome null.
The downstream protocol opens development families first and forbids opening
either holdout family until the predictor artifact itself is hash-frozen.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ATLAS = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-record-lineage-atlas-v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-prospective-ordinary-family-holdout-v1.json"
DEFAULT_SALT = "r17-norm12-ordinary-family-holdout-v1"
PRIME_BLOCKS = (
    (19, 41, 43, 61, 71, 73, 79, 83),
    (89, 107, 113, 127, 131, 137, 139, 151),
    (157, 163, 167, 173, 179, 181, 191, 193, 197),
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_digest(value) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def parameter_text(numerator: int, denominator: int) -> str:
    return str(numerator) if denominator == 1 else f"{numerator}/{denominator}"


def homogeneous_value(coefficients: list[str], numerator: int, denominator: int) -> int:
    degree = 24
    return sum(
        int(coefficient) * numerator**index * denominator ** (degree - index)
        for index, coefficient in enumerate(coefficients)
    )


def j_fingerprint(chart: dict[str, object], numerator: int, denominator: int) -> str | None:
    j_map = chart["normalized_j_map"]
    assert isinstance(j_map, dict)
    j_numerator = homogeneous_value(
        j_map["numerator_coefficients_low_to_high"], numerator, denominator
    )
    j_denominator = homogeneous_value(
        j_map["denominator_coefficients_low_to_high"], numerator, denominator
    )
    if j_denominator == 0:
        return None
    value = Fraction(j_numerator, j_denominator)
    return sha256(f"{value.numerator}/{value.denominator}".encode()).hexdigest()


def draw_pair(
    family: str, counter: int, salt: str, minimum_height: int, maximum_height: int
) -> tuple[int, int]:
    block = sha256(f"{salt}|parameter|{family}|{counter}".encode()).digest()
    numerator = int.from_bytes(block[:16], "big") % (2 * maximum_height + 1)
    numerator -= maximum_height
    denominator = int.from_bytes(block[16:], "big") % maximum_height + 1
    return numerator, denominator


def family_inventory(atlas: dict[str, object]):
    atlas_body = atlas["atlas"]
    if atlas_body["pgl2_equivalence_class_count"] != 6:
        raise ArithmeticError("the exact atlas no longer has six PGL2 classes")
    charts = {row["label"]: row for row in atlas_body["charts"]}
    families = []
    for equivalence_class in atlas_body["pgl2_equivalence_classes"]:
        representative = equivalence_class["representative"]
        chart = charts[representative]
        members = [row["label"] for row in equivalence_class["members"]]
        if representative not in members:
            raise ArithmeticError("a PGL2 representative is absent from its class")
        if not all(row["solve_certificate"]["identity_verified"] for row in equivalence_class["members"]):
            raise ArithmeticError("a PGL2 family member lacks an exact identity")
        families.append(
            {
                "representative": representative,
                "frame_class": chart["frame_class"],
                "member_charts": members,
                "member_count": len(members),
                "representative_j_map_sha256": chart["normalized_j_map"]["sha256"],
                "chart": chart,
            }
        )
    if sum(row["member_count"] for row in families) != 43:
        raise ArithmeticError("the six PGL2 families no longer partition 43 charts")
    return sorted(families, key=lambda row: row["representative"])


def build(
    *,
    minimum_height: int,
    maximum_height: int,
    rows_per_family: int,
    holdout_family_count: int,
    salt: str,
):
    if not (30_000 < minimum_height <= maximum_height):
        raise ValueError("the prospective shell must begin beyond height 30000")
    if rows_per_family < 1:
        raise ValueError("rows_per_family must be positive")
    atlas = json.loads(ATLAS.read_text())
    families = family_inventory(atlas)
    if not (1 <= holdout_family_count < len(families)):
        raise ValueError("holdout_family_count must leave development families")

    split_order = sorted(
        (row["representative"] for row in families),
        key=lambda family: sha256(f"{salt}|family-split|{family}".encode()).digest(),
    )
    holdout_families = set(split_order[:holdout_family_count])
    seen_j: set[str] = set()
    rows = []
    structural_rejections = Counter()
    for family_record in families:
        family = family_record["representative"]
        selected = 0
        counter = 0
        seen_pairs: set[tuple[int, int]] = set()
        while selected < rows_per_family:
            numerator, denominator = draw_pair(
                family, counter, salt, minimum_height, maximum_height
            )
            selection_counter = counter
            counter += 1
            height = max(abs(numerator), denominator)
            if height < minimum_height:
                structural_rejections["below_shell"] += 1
                continue
            if gcd(abs(numerator), denominator) != 1:
                structural_rejections["nonprimitive_pair"] += 1
                continue
            pair = (numerator, denominator)
            if pair in seen_pairs:
                structural_rejections["duplicate_parameter"] += 1
                continue
            fingerprint = j_fingerprint(
                family_record["chart"], numerator, denominator
            )
            if fingerprint is None:
                structural_rejections["singular_fibre"] += 1
                continue
            if fingerprint in seen_j:
                structural_rejections["duplicate_fibre_across_cohort"] += 1
                continue
            seen_pairs.add(pair)
            seen_j.add(fingerprint)
            selected += 1
            row_key = {
                "family": family,
                "projective_pair": [numerator, denominator],
            }
            rows.append(
                {
                    "sample_id": canonical_digest(row_key)[:24],
                    "family": family,
                    "frame_class": family_record["frame_class"],
                    "pgl2_family_members": family_record["member_charts"],
                    "outer_split": (
                        "locked_family_holdout"
                        if family in holdout_families
                        else "prospective_development"
                    ),
                    "parameter": parameter_text(numerator, denominator),
                    "projective_pair": [numerator, denominator],
                    "projective_height": height,
                    "selection_counter": selection_counter,
                    "selection_lane": "ordinary_counter_hash_draw_no_score",
                    "j_invariant_sha256": fingerprint,
                    "features": None,
                    "search_outcome": None,
                    "search_status": "NOT_OPENED",
                }
            )

    family_counts = Counter(row["family"] for row in rows)
    split_counts = Counter(row["outer_split"] for row in rows)
    if set(family_counts.values()) != {rows_per_family}:
        raise ArithmeticError("ordinary sampling is not balanced by whole family")
    if len({row["sample_id"] for row in rows}) != len(rows):
        raise ArithmeticError("sample identifiers are not unique")
    if len({row["j_invariant_sha256"] for row in rows}) != len(rows):
        raise ArithmeticError("the prospective cohort contains duplicate fibres")
    if any(row["features"] is not None or row["search_outcome"] is not None for row in rows):
        raise ArithmeticError("a prospective label was opened during commitment")

    development_families = [
        family for family in split_order if family not in holdout_families
    ]
    locked_families = [family for family in split_order if family in holdout_families]
    payload = {
        "schema": "elkies-k3.r17-norm12-prospective-ordinary-family-holdout.v1",
        "status": "FROZEN_UNOPENED_PROSPECTIVE_WHOLE_FAMILY_HOLDOUT",
        "commitment": {
            "frozen_before_any_search_outcome": True,
            "feature_values_materialized": False,
            "search_labels_materialized": False,
            "selection_used_public_rank_or_hit_labels": False,
            "historical_search_denominator_known": False,
            "salt": salt,
        },
        "population": {
            "parameter_domain": (
                "primitive (a,b), b>0, with minimum_height<=max(abs(a),b)<=maximum_height"
            ),
            "minimum_height": minimum_height,
            "maximum_height": maximum_height,
            "rows_per_family": rows_per_family,
            "family_count": len(families),
            "scheduled_parameter_count": len(rows),
            "selection": (
                "first accepted counter-mode SHA-256 draws; reject only nonprimitive/out-of-shell "
                "pairs, singular fibres, and duplicate exact j-invariants"
            ),
            "ordinary_means": (
                "no Nagao, conductor, public-rank, cover-visibility, or point-search score "
                "enters parameter selection"
            ),
            "structural_rejections_before_commitment": dict(
                sorted(structural_rejections.items())
            ),
            "family_counts": dict(sorted(family_counts.items())),
            "split_counts": dict(sorted(split_counts.items())),
        },
        "pgl2_family_split": {
            "split_unit": "entire exact rational-PGL2 equivalence class of j-maps",
            "selection": "smallest salted family-level SHA-256 keys are locked holdouts",
            "family_order_by_frozen_hash": split_order,
            "prospective_development_families": development_families,
            "locked_holdout_families": locked_families,
            "holdout_family_count": holdout_family_count,
            "holdout_parameter_count": holdout_family_count * rows_per_family,
            "development_parameter_count": (
                len(families) - holdout_family_count
            ) * rows_per_family,
            "member_charts_by_representative": {
                row["representative"]: row["member_charts"] for row in families
            },
            "leakage_rule": (
                "No fibre, transformed parameter, feature, label, threshold, or fitted value "
                "from any member chart of a locked family may enter model fitting or tuning."
            ),
        },
        "frozen_search_protocol": {
            "protocol_id": "r17-norm12-uniform-extra-direction-search-v1",
            "execution_order": [
                "run all prospective-development families under this unchanged protocol",
                "freeze the complete predictor artifact and its SHA-256",
                "open both locked PGL2 families without retraining, threshold changes, or depth changes",
            ],
            "per_parameter_stages": [
                {
                    "stage": "exact_specialization_gate",
                    "rule": (
                        "specialize the representative degree-(8,12) model and a saturated "
                        "generic MW17 basis; verify every section and an exact finite-reduction "
                        "rank-17 certificate on the same global minimal curve"
                    ),
                    "failure": "CENSORED_STRUCTURAL; never a rank label",
                },
                {
                    "stage": "fixed_local_features",
                    "prime_blocks": [list(block) for block in PRIME_BLOCKS],
                    "rule": (
                        "compute every usable a_p exactly; bad or model-degenerate primes have "
                        "the predeclared zero standardized contribution"
                    ),
                    "mathematical_status": "heuristic features only",
                },
                {
                    "stage": "uniform_bounded_point_search",
                    "engine": "Sage eclib mwrank_MordellWeil.search",
                    "global_minimal_model_required": True,
                    "search_height": 12,
                    "maxr": 32,
                    "pp": 0,
                    "wall_clock_limit_seconds": 300,
                    "adaptive_depth": False,
                    "presieve": False,
                    "timeout_or_backend_failure": "CENSORED_OPERATIONAL",
                },
                {
                    "stage": "exact_positive_confirmation",
                    "rule": (
                        "transport each returned point back exactly, verify the curve equation, "
                        "and retain a positive only after exact finite-quotient escape from MW17"
                    ),
                    "positive_label": "CERTIFIED_DISPLAYED_QUOTIENT_GAIN_AT_LEAST_ONE",
                    "completed_miss_label": "BOUNDED_PROTOCOL_NO_GAIN_FOUND",
                },
            ],
            "allocation_rule": "identical stages, bounds, and time limit for every scheduled row",
            "no_promotion_rules": [
                "A completed bounded miss is not a rank-17 or quotient-zero theorem.",
                "A timeout is censored and is not counted as a completed miss.",
                "Nagao scores may rank no row and may not change search depth in this cohort.",
                "Public recognized fibres are diagnostics only and may not fit or tune the predictor.",
            ],
        },
        "families": [
            {key: value for key, value in row.items() if key != "chart"}
            for row in families
        ],
        "rows": rows,
        "inputs": {relative(ATLAS): digest(ATLAS)},
        "generation": {
            "script": relative(Path(__file__)),
            "script_sha256": digest(Path(__file__)),
            "command": (
                ".venv/bin/python elkies-k3/scripts/"
                "build_r17_norm12_prospective_family_holdout.py"
            ),
        },
        "proof_boundary": [
            "The exact PGL2 family partition, parameter draws, and absence of duplicate sampled fibres are certified.",
            "The cohort supplies a known prospective search denominator only after protocol execution; it does not recover the unknown historical denominator.",
            "No rank, quotient, point-search, or predictor-performance outcome has been computed for any sampled row.",
            "Whole-family locking blocks direct fibre-level family leakage, but cannot erase human feature-design choices made before this commitment.",
        ],
    }
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-height", type=int, default=30_001)
    parser.add_argument("--maximum-height", type=int, default=60_000)
    parser.add_argument("--rows-per-family", type=int, default=256)
    parser.add_argument("--holdout-family-count", type=int, default=2)
    parser.add_argument("--salt", default=DEFAULT_SALT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = build(
        minimum_height=args.minimum_height,
        maximum_height=args.maximum_height,
        rows_per_family=args.rows_per_family,
        holdout_family_count=args.holdout_family_count,
        salt=args.salt,
    )
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != serialized:
            raise SystemExit("stale prospective whole-family holdout commitment")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    split = payload["pgl2_family_split"]
    print(
        "R17PROSPECTIVEFAMILY"
        f"|rows={payload['population']['scheduled_parameter_count']}"
        f"|development_families={len(split['prospective_development_families'])}"
        f"|holdout_families={len(split['locked_holdout_families'])}"
        f"|holdouts={','.join(split['locked_holdout_families'])}"
        "|status=FROZEN_UNOPENED",
        flush=True,
    )


if __name__ == "__main__":
    main()
