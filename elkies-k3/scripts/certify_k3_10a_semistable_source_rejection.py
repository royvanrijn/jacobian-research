#!/usr/bin/env python3
"""Aggregate the two-auxiliary ideal-source rejection for the det-750 MW17 K3."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "artifacts/generated-results"
ADAPTERS = (
    GEN / "elkies-k3-k3-10a14a46c14b3150-source-search-target-v1.json",
    GEN / "elkies-k3-k3-10a14a46c14b3150-source-search-target-partner2-v1.json",
)
SEARCHES = (
    GEN / "elkies-k3-k3-10a14a46c14b3150-semistable-mw0-2-sources-large-a-v1.json",
    GEN / "elkies-k3-k3-10a14a46c14b3150-semistable-mw0-2-sources-large-a-partner2-v1.json",
)
PARETO = GEN / "elkies-k3-rank7-surface-pareto-v1.json"
DEFAULT_OUTPUT = GEN / "elkies-k3-k3-10a14a46c14b3150-semistable-source-rejection-v1.json"
SURFACE_ID = "K3-10a14a46c14b3150"
TARGET_ID = "K3-10a14a46c14b3150-F001"
AMBIENTS = ["D24", "D16_E8", "3E8", "A24", "A17_E7", "A15_D9"]


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    adapters = [json.loads(path.read_text()) for path in ADAPTERS]
    searches = [json.loads(path.read_text()) for path in SEARCHES]
    pareto = json.loads(PARETO.read_text())
    surface = next(row for row in pareto["surfaces"] if row["surface_id"] == SURFACE_ID)

    if [row["auxiliary"]["partner_index_one_based"] for row in adapters] != [1, 2]:
        raise ArithmeticError("the two auxiliary classes are not both represented")
    if any(row["surface_id"] != SURFACE_ID for row in adapters):
        raise ArithmeticError("adapter surface mismatch")
    frame_hashes = {row["frame"]["gram_sha256"] for row in adapters}
    if len(frame_hashes) != 1:
        raise ArithmeticError("adapters do not share the same target frame")
    target = adapters[0]["frame"]
    if not (
        target["frame_id"] == TARGET_ID
        and target["root_rank"] == 0
        and target["mw_rank_for_rho_19"] == 17
        and target["determinant"] == 750
    ):
        raise ArithmeticError("target is no longer the rootless MW17 frame")

    summaries = []
    for adapter, search in zip(adapters, searches):
        scope = search["search_scope"]
        if search["status"] != "PASS_EXACT_SEARCH_NO_SUCCESS_HIT":
            raise ArithmeticError("a determinant-750 ideal source has appeared")
        if not scope["complete"] or not scope["all_a_only"]:
            raise ArithmeticError("source search is not complete in the declared cut")
        if scope["source_root_rank"] != [15, 17] or scope["source_support_count"] != [1, 3]:
            raise ArithmeticError("source search cut changed")
        if scope["rooted_niemeier_ambients"] != AMBIENTS:
            raise ArithmeticError("large-ambient list changed")
        accounting = search["accounting"]
        if accounting["distinct_reduced_gram_sources"] or accounting["success_condition_hits"]:
            raise ArithmeticError("negative source accounting is inconsistent")
        summaries.append(
            {
                "partner_index_one_based": adapter["auxiliary"]["partner_index_one_based"],
                "auxiliary_gram_sha256": adapter["auxiliary"]["gram_sha256"],
                "ambient_embedding_counts": {
                    row["ambient_label"]: row["totals"].get("complete_auxiliary_embeddings", 0)
                    for row in accounting["ambient_searches"]
                },
                "admissible_sources": 0,
            }
        )

    if not (
        surface["maximum_generic_mw_rank"] == 17
        and surface["moduli"]["rationality"] is True
        and surface["moduli"]["genus"] == 0
        and surface["moduli"]["t_arithmetic"]["curve_label"] == "Gamma_0(3)"
    ):
        raise ArithmeticError("surface Pareto arithmetic changed")

    payload = {
        "schema": "elkies-k3.k3-10a-semistable-source-rejection.v1",
        "status": "PASS_COMPLETE_TWO_AUXILIARY_LARGE_AMBIENT_IDEAL_SOURCE_REJECTION",
        "surface": {
            "surface_id": SURFACE_ID,
            "determinant": 750,
            "moduli_genus": 0,
            "moduli_rational": True,
            "modular_curve": "Gamma_0(3)",
        },
        "target": {
            "frame_id": TARGET_ID,
            "rootless": True,
            "mw_rank_for_rho_19": 17,
            "frame_gram_sha256": target["gram_sha256"],
        },
        "ideal_source_cut": {
            "mw_rank": [0, 2],
            "semistable": True,
            "reducible_support_count": [1, 3],
            "rooted_niemeier_ambients": AMBIENTS,
            "auxiliary_classes": 2,
            "searches": summaries,
        },
        "decision": {
            "action": "DEMOTE_BEHIND_DETERMINANT_500_CANDIDATE",
            "reason": (
                "The rational genus-zero moduli and rootless MW17 target remain attractive, "
                "but neither auxiliary class has a semistable MW0-2 source in the complete "
                "declared six-large-ambient cut."
            ),
        },
        "proof_boundary": {
            "proved": (
                "Both catalogue auxiliary classes have been searched completely in the "
                "declared six-large-ambient, all-A, root-rank-15--17, support-one--three cut, "
                "and no admissible source Gram occurs."
            ),
            "not_proved": (
                "This is not a complete fibration classification. Smaller-component "
                "Niemeier ambients, nonsemistable sources, MW3 sources, and sources with "
                "four or more reducible supports remain open."
            ),
        },
        "inputs": {
            **{relative(path): digest(path) for path in ADAPTERS + SEARCHES},
            relative(PARETO): digest(PARETO),
        },
        "reproduce": "python3 elkies-k3/scripts/certify_k3_10a_semistable_source_rejection.py",
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output_path = arguments.output.resolve()
    if arguments.check:
        if not output_path.exists() or output_path.read_text() != serialized:
            raise SystemExit(f"stale artifact: {output_path}")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print("K310AIDEALSOURCE|auxiliaries=2|sources=0|status=PASS")


if __name__ == "__main__":
    main()
