#!/usr/bin/env python3
"""Build the target-screen atlas stratified by actual generic root rank."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
TABLE = GENERATED / "elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.tsv"
STRATA = GENERATED / "elkies-k3-icarm-11952-norm8-low-root-strata-v1.json"
OLD_ATLAS = GENERATED / "elkies-k3-icarm-11952-norm8-a1-mw16-atlas-v1.json"
FIXED_MW15 = GENERATED / "elkies-k3-icarm-curve302-273-fixed-2a1-mw15-screen-v1.json"
DEFAULT_OUTPUT = GENERATED / "elkies-k3-icarm-11952-norm8-low-root-atlas-v2.json"
TARGETS = (302, 273, 542, 548, 398, 399, 400, 403, 401, 402, 10)
COMPILED_TARGETS = {398, 400, 401, 542, 548}


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def paths(curve_id: int) -> tuple[Path, Path, Path]:
    if curve_id == 398:
        modular = GENERATED / "elkies-k3-curve398-11952-norm8-a1-modular-screen-v1.json"
        exact = GENERATED / "elkies-k3-curve398-11952-norm8-a1-exact-survivors-v1.json"
    else:
        prefix = GENERATED / f"elkies-k3-icarm-curve{curve_id}-11952-norm8-a1"
        modular = Path(f"{prefix}-modular-screen-v1.json")
        exact = Path(f"{prefix}-exact-survivors-v1.json")
    compiled = GENERATED / f"elkies-k3-icarm-curve{curve_id}-11952-norm8-a1-compiled-survivors-v1.json"
    return modular, exact, compiled


def qq_hit(record: dict) -> bool:
    return any(
        bool(
            specialization.get("isomorphic_to_target_over_Q")
            or specialization.get("isomorphic_to_curve398_over_Q")
        )
        for specialization in record.get("specializations", [])
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    stratum_certificate = json.loads(STRATA.read_text())
    if stratum_certificate.get("status") != "PASS_EXACT_COMPLETE_NORM8_LOW_ROOT_STRATIFICATION":
        raise ArithmeticError("low-root stratum certificate is not passing")
    strata = {
        int(row["minimum_unoriented_split_member_count"]): row
        for row in stratum_certificate["strata"]
    }
    if set(strata) != set(range(1, 9)):
        raise ArithmeticError("expected the complete m=1,...,8 stratum range")

    multiplicity_by_rank = {}
    with TABLE.open(newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            rank = int(row["priority_rank"])
            if rank in multiplicity_by_rank:
                raise ArithmeticError("duplicate priority rank")
            multiplicity_by_rank[rank] = int(row["minimal_unoriented_count"])
    if set(multiplicity_by_rank) != set(range(1, 63918)):
        raise ArithmeticError("priority table is not complete")
    if Counter(multiplicity_by_rank.values()) != Counter(
        {multiplicity: row["class_count"] for multiplicity, row in strata.items()}
    ):
        raise ArithmeticError("priority-table and stratum counts differ")

    old_atlas = json.loads(OLD_ATLAS.read_text())
    if old_atlas.get("status") != "PASS_EXACT_COMPLETE_PRIORITY_ICARM_A1_MW16_ATLAS":
        raise ArithmeticError("historical target atlas is not passing")
    if old_atlas["scope"]["target_curve_order"] != list(TARGETS):
        raise ArithmeticError("historical target order changed")
    old_target_by_id = {
        int(row["curve_id"]): row for row in old_atlas["targets"]
    }

    input_paths = [Path(__file__).resolve(), TABLE, STRATA, OLD_ATLAS, FIXED_MW15]
    target_rows = []
    positive_control_ranks = []
    for curve_id in TARGETS:
        modular_path, exact_path, compiled_path = paths(curve_id)
        modular = json.loads(modular_path.read_text())
        exact = json.loads(exact_path.read_text())
        input_paths.extend((modular_path, exact_path))
        search = modular["search"]
        if int(search["priority_table_class_count"]) != 63917:
            raise ArithmeticError(f"curve {curve_id}: modular screen is incomplete")
        survivors = list(map(int, search["survivor_priority_ranks"]))
        witnesses = search.get("first_exclusion_prime_by_priority_rank")
        if witnesses is not None:
            if len(witnesses) != 63917:
                raise ArithmeticError(f"curve {curve_id}: exclusion vector is incomplete")
            if {index + 1 for index, value in enumerate(witnesses) if value is None} != set(survivors):
                raise ArithmeticError(f"curve {curve_id}: survivors and witnesses do not partition")
        elif int(search["excluded_count"]) + len(survivors) != 63917:
            raise ArithmeticError(f"curve {curve_id}: legacy screen does not partition the atlas")
        exact_records = {int(record["priority_rank"]): record for record in exact.get("records", [])}
        if set(exact_records) != set(survivors):
            raise ArithmeticError(f"curve {curve_id}: exact stage does not resolve every survivor")
        hit_ranks = sorted(rank for rank, record in exact_records.items() if qq_hit(record))
        compiled_ranks = []
        if curve_id in COMPILED_TARGETS:
            compiled = json.loads(compiled_path.read_text())
            input_paths.append(compiled_path)
            if compiled.get("status") != "PASS_EXACT_COMPILED_A1_MW16_FIBRATIONS":
                raise ArithmeticError(f"curve {curve_id}: compiled certificate is not passing")
            compiled_ranks = sorted(int(row["priority_rank"]) for row in compiled["fibrations"])
            if compiled_ranks != hit_ranks:
                raise ArithmeticError(f"curve {curve_id}: exact/compiled hit ranks differ")
            for row in compiled["fibrations"]:
                if (
                    row["equation"]["fibre_configuration"] != "I2 at infinity + 22 I1"
                    or int(row["generic_mordell_weil"]["rank"]) != 16
                    or row["generic_mordell_weil"]["height_gram_determinant"] != "474"
                ):
                    raise ArithmeticError(f"curve {curve_id}: compiled MW16 invariant changed")
        elif hit_ranks:
            raise ArithmeticError(f"curve {curve_id}: uncompiled QQ-isomorphic hit")
        if any(multiplicity_by_rank[rank] != 1 for rank in hit_ranks):
            raise ArithmeticError(f"curve {curve_id}: a target hit lies outside the true A1 stratum")
        if curve_id == 398:
            positive_control_ranks = hit_ranks

        per_stratum = []
        for multiplicity, stratum in sorted(strata.items()):
            class_ranks = {
                rank for rank, value in multiplicity_by_rank.items()
                if value == multiplicity
            }
            stratum_survivors = sorted(class_ranks & set(survivors))
            stratum_hits = sorted(class_ranks & set(hit_ranks))
            excluded = len(class_ranks) - len(stratum_survivors)
            if any(rank not in exact_records for rank in stratum_survivors):
                raise ArithmeticError("a stratum survivor is not exactly resolved")
            outcome = "HIT" if stratum_hits else "MISS_EXACT_COMPLETE_STRATUM"
            per_stratum.append(
                {
                    "minimum_unoriented_split_member_count": multiplicity,
                    "root_lattice": stratum["root_lattice"],
                    "root_rank": int(stratum["root_rank"]),
                    "geometric_mw_rank_at_rho_19": int(
                        stratum["geometric_mw_rank_at_rho_19"]
                    ),
                    "class_count": len(class_ranks),
                    "modular_excluded_count": excluded,
                    "exact_survivor_count": len(stratum_survivors),
                    "exact_survivor_priority_ranks": stratum_survivors,
                    "qq_isomorphic_hit_count": len(stratum_hits),
                    "qq_isomorphic_hit_priority_ranks": stratum_hits,
                    "outcome": outcome,
                }
            )
        target_rows.append(
            {
                "curve_id": curve_id,
                "snapshot_rank_lower_bound": int(
                    old_target_by_id[curve_id]["snapshot_rank_lower_bound"]
                ),
                "total_class_count": 63917,
                "modular_excluded_count": int(search["excluded_count"]),
                "exact_survivor_count": len(survivors),
                "qq_isomorphic_hit_count": len(hit_ranks),
                "qq_isomorphic_hit_priority_ranks": hit_ranks,
                "strata": per_stratum,
                "outcome": "HIT" if hit_ranks else "MISS_EXACT_COMPLETE_NORM8_LOW_ROOT_ATLAS",
            }
        )

    if positive_control_ranks != [16875, 63669]:
        raise ArithmeticError("curve-398 positive control changed")
    curve302 = next(row for row in target_rows if row["curve_id"] == 302)
    if curve302["outcome"] != "MISS_EXACT_COMPLETE_NORM8_LOW_ROOT_ATLAS":
        raise ArithmeticError("curve 302 no longer misses the complete norm-eight atlas")
    if any(row["outcome"] != "MISS_EXACT_COMPLETE_STRATUM" for row in curve302["strata"]):
        raise ArithmeticError("curve 302 has a low-root stratum hit")

    fixed_mw15 = json.loads(FIXED_MW15.read_text())
    if fixed_mw15.get("status") != "PASS_EXACT_BOUNDED_FIXED_2A1_MW15_TARGET_SCREEN":
        raise ArithmeticError("fixed-corridor MW15 control is not passing")
    fixed302 = next(row for row in fixed_mw15["records"] if int(row["curve_id"]) == 302)
    if fixed302["status"] != "PASS_MODULAR_EXCLUDED_FIXED_2A1_FIBRATION":
        raise ArithmeticError("fixed-corridor curve-302 outcome changed")

    payload = {
        "schema": "elkies-k3.icarm-11952-norm8-low-root-atlas.v2",
        "status": "PASS_EXACT_COMPLETE_PRIORITY_ICARM_NORM8_LOW_ROOT_ATLAS",
        "source_chart": "norm12-orbit-11952 alternate-Q80 rootless/MW17",
        "scope": {
            "target_curve_order": list(TARGETS),
            "target_curve_count": len(TARGETS),
            "classes_per_target": 63917,
            "target_fibration_pairs_screened": 63917 * len(TARGETS),
            "root_rank_range": [1, 8],
            "geometric_mw_rank_range": [9, 16],
            "curve302_complete_norm8_miss": True,
            "curve398_positive_control_priority_ranks": positive_control_ranks,
        },
        "stratum_totals": [strata[index] for index in sorted(strata)],
        "targets": target_rows,
        "additional_fixed_corridor_control": {
            "scope": "one independent source-identified 2A1/MW15 fibration; not an atlas",
            "curve302_first_exclusion_prime": int(fixed302["first_exclusion_prime"]),
            "curve302_outcome": fixed302["status"],
        },
        "correction": (
            "The prior v1 atlas uniformly called all 63,917 classes A1/MW16. "
            "The exact split-member and discriminant certificates show that only 1,266 "
            "classes have A1/MW16; 8,410 have 2A1/MW15, and the remaining classes "
            "continue through 8A1/MW9. The target-j computations themselves covered "
            "all classes and remain valid."
        ),
        "curve302_conclusion": (
            "Curve 302 has no rational target parameter in any of the eight complete "
            "root-rank strata of the committed 63,917-class norm-eight old-degree-two "
            "atlas. In particular it misses all 1,266 A1/MW16 and all 8,410 "
            "2A1/MW15 classes in that atlas."
        ),
        "proof_boundary": (
            "Complete means the minimum-norm-eight, old-degree-two residual-chord "
            "translation layer on source chart 11952. This is not the global fibration "
            "atlas of X948: A2/MW15, other old degrees, other trace norms, and other "
            "rootless source charts remain outside the theorem. A bounded miss or this "
            "layer-complete miss is not a nonexistence theorem for a parent of curve 302."
        ),
        "inputs": {relative(path): digest(path) for path in input_paths},
        "reproducing_command": (
            "python3 elkies-k3/scripts/build_icarm_norm8_low_root_atlas.py"
        ),
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output_path = args.output.resolve()
    if args.check:
        if not output_path.is_file() or output_path.read_text() != serialized:
            raise ArithmeticError("stored low-root ICARM atlas differs from replay")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(
        "ICARMLOWROOTATLAS|targets={}|pairs={}|curve302=MISS|mw16={}|mw15={}|control398={}|status={}|output={}".format(
            len(TARGETS),
            63917 * len(TARGETS),
            strata[1]["class_count"],
            strata[2]["class_count"],
            ",".join(map(str, positive_control_ranks)),
            payload["status"],
            relative(output_path),
        )
    )


if __name__ == "__main__":
    main()
