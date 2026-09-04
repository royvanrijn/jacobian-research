#!/usr/bin/env python3
"""Build the compact exact ICARM A1/MW16 atlas from per-target certificates."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
DEFAULT_OUTPUT = GENERATED / "elkies-k3-icarm-11952-norm8-a1-mw16-atlas-v1.json"

# User-directed order: priority pair, named follow-ups, then the remaining
# unexplained rank >= 24 misses.  Curve 398 is retained as the positive control.
TARGETS = (302, 273, 542, 548, 398, 399, 400, 403, 401, 402, 10)
COMPILED_TARGETS = {398, 400, 401, 542, 548}


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def paths(curve_id: int):
    if curve_id == 398:
        modular = GENERATED / "elkies-k3-curve398-11952-norm8-a1-modular-screen-v1.json"
        exact = GENERATED / "elkies-k3-curve398-11952-norm8-a1-exact-survivors-v1.json"
    else:
        prefix = GENERATED / f"elkies-k3-icarm-curve{curve_id}-11952-norm8-a1"
        modular = Path(f"{prefix}-modular-screen-v1.json")
        exact = Path(f"{prefix}-exact-survivors-v1.json")
    compiled = GENERATED / f"elkies-k3-icarm-curve{curve_id}-11952-norm8-a1-compiled-survivors-v1.json"
    return modular, exact, compiled


def target_rank(curve_id: int, modular: dict, exact: dict) -> int:
    for document in (exact, modular):
        target = document.get("target")
        if isinstance(target, dict) and "snapshot_rank_lower_bound" in target:
            return int(target["snapshot_rank_lower_bound"])
    if curve_id == 398:
        return 30
    raise ArithmeticError(f"curve {curve_id}: rank lower bound missing")


def exact_qq_hit_count(curve_id: int, exact: dict) -> int:
    if "qq_isomorphic_candidate_count" in exact:
        return int(exact["qq_isomorphic_candidate_count"])
    key = "isomorphic_to_curve398_over_Q" if curve_id == 398 else "isomorphic_to_target_over_Q"
    return sum(
        bool(specialization.get(key))
        for record in exact.get("records", [])
        for specialization in record.get("specializations", [])
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    input_paths = [Path(__file__).resolve()]
    rows = []
    hit_fibrations = []
    for curve_id in TARGETS:
        modular_path, exact_path, compiled_path = paths(curve_id)
        modular = json.loads(modular_path.read_text())
        exact = json.loads(exact_path.read_text())
        input_paths.extend((modular_path, exact_path))
        search = modular["search"]
        class_count = int(search["priority_table_class_count"])
        survivors = list(map(int, search["survivor_priority_ranks"]))
        if class_count != 63917:
            raise ArithmeticError(f"curve {curve_id}: incomplete A1 class count")
        qq_hits = exact_qq_hit_count(curve_id, exact)
        compiled_count = 0
        compiled_ranks = []
        if curve_id in COMPILED_TARGETS:
            compiled = json.loads(compiled_path.read_text())
            input_paths.append(compiled_path)
            if compiled.get("status") != "PASS_EXACT_COMPILED_A1_MW16_FIBRATIONS":
                raise ArithmeticError(f"curve {curve_id}: compiled status changed")
            compiled_count = int(compiled["compiled_fibration_count"])
            compiled_ranks = [int(row["priority_rank"]) for row in compiled["fibrations"]]
            if compiled_count != qq_hits or len(compiled_ranks) != compiled_count:
                raise ArithmeticError(f"curve {curve_id}: exact/compiled hit count mismatch")
            for fibration in compiled["fibrations"]:
                generic = fibration["generic_mordell_weil"]
                equation = fibration["equation"]
                if (
                    int(generic["rank"]) != 16
                    or generic["height_gram_determinant"] != "474"
                    or equation["fibre_configuration"] != "I2 at infinity + 22 I1"
                ):
                    raise ArithmeticError(f"curve {curve_id}: compiled fibration invariant changed")
                hit_fibrations.append({
                    "curve_id": curve_id,
                    "priority_rank": int(fibration["priority_rank"]),
                    "orbit_hex": fibration["orbit_hex"],
                    "complete_old_degree_one_section_count": int(
                        fibration["complete_old_degree_one_section_count"]
                    ),
                    "generic_mw_rank": 16,
                    "generic_mw_height_determinant": "474",
                    "fibre_configuration": "I2 at infinity + 22 I1",
                    "target_rank_lower_bound": target_rank(curve_id, modular, exact),
                    "specialization_rank_jump_lower_bound": (
                        target_rank(curve_id, modular, exact) - 16
                    ),
                })
        elif qq_hits:
            raise ArithmeticError(f"curve {curve_id}: uncompiled QQ-isomorphic hit")

        rank = target_rank(curve_id, modular, exact)
        rows.append({
            "curve_id": curve_id,
            "snapshot_rank_lower_bound": rank,
            "complete_a1_class_count": class_count,
            "modular_survivor_count": len(survivors),
            "modular_survivor_priority_ranks": survivors,
            "exact_qq_isomorphic_count": qq_hits,
            "compiled_fibration_count": compiled_count,
            "compiled_priority_ranks": compiled_ranks,
            "outcome": "HIT" if compiled_count else "MISS_COMPLETE_A1_LAYER",
        })

    if len(hit_fibrations) != 9:
        raise ArithmeticError(f"expected nine compiled fibrations, got {len(hit_fibrations)}")
    hit_curve_ids = [row["curve_id"] for row in rows if row["outcome"] == "HIT"]
    miss_curve_ids = [row["curve_id"] for row in rows if row["outcome"] != "HIT"]
    payload = {
        "schema": "elkies-k3.icarm-11952-norm8-a1-mw16-atlas.v1",
        "status": "PASS_EXACT_COMPLETE_PRIORITY_ICARM_A1_MW16_ATLAS",
        "source_chart": "norm12-orbit-11952 alternate Q80 rootless/MW17",
        "pipeline": [
            "X948 rootless frame",
            "complete minimum-norm-eight translation classes",
            "old-degree-two residual-chord A1/MW16 pencil",
            "projective modular target-j screen",
            "exact QQ factorization of survivors",
            "twist-sensitive QQ isomorphism",
            "exact equation and saturated generic-MW compile",
        ],
        "scope": {
            "target_curve_count": len(TARGETS),
            "target_curve_order": list(TARGETS),
            "classes_per_target": 63917,
            "target_fibration_pairs_screened": 63917 * len(TARGETS),
            "hit_curve_ids": hit_curve_ids,
            "miss_curve_ids": miss_curve_ids,
            "compiled_fibration_count": len(hit_fibrations),
        },
        "targets": rows,
        "compiled_fibrations": hit_fibrations,
        "proof_boundary": (
            "MISS_COMPLETE_A1_LAYER excludes rational target parameters only in the complete "
            "63,917-class minimum-norm-eight A1/MW16 layer on source chart 11952. HIT "
            "certifies the displayed polynomial pencil, I2+22I1 fibres, saturated generic "
            "MW16 height lattice, rational parameter, and QQ target isomorphism. Target ranks "
            "and displayed jumps are lower bounds from the pinned ICARM snapshot, not exact-rank proofs."
        ),
        "inputs": {relative(path): digest(path) for path in input_paths},
        "reproducing_command": (
            "python3 elkies-k3/scripts/build_icarm_a1_mw16_atlas.py"
        ),
    }
    output_path = args.output.resolve()
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not output_path.is_file() or output_path.read_text() != serialized:
            raise ArithmeticError("stored A1/MW16 atlas differs from replay")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(
        f"ICARMA1ATLAS|targets={len(TARGETS)}|pairs={63917 * len(TARGETS)}"
        f"|hit_curves={len(hit_curve_ids)}|fibrations={len(hit_fibrations)}"
        f"|status={payload['status']}|output={relative(output_path)}"
    )


if __name__ == "__main__":
    main()
