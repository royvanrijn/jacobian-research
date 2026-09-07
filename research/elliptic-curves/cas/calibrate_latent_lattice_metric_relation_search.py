#!/usr/bin/env python3
"""Freeze the bounded metric/relation/global-replay Phase-0 search.

No wgxli target is loaded.  Training R17 height forms define a secondary
metric proposal score; the held-out rank-25 public cloud is searched without
its embedding.  Exact integral lifting, primitive/global replay thresholds,
and withheld overlap evaluation are separate stages.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import importlib
import json
from pathlib import Path
import platform
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC = ROOT / "elliptic-curves"
sys.path[:0] = [str(ELLIPTIC), str(ELLIPTIC / "cas")]

from latent_lattice import (  # noqa: E402
    EllipticCurve,
    bounded_metric_relation_search,
    build_relation_complex,
    height_gram,
)


ARTIFACTS = ROOT / "artifacts/generated-results/elliptic-curves"
CONSENSUS = ARTIFACTS / "latent_lattice_relation_consensus_v1.json"
TRUTH = ARTIFACTS / "latent_lattice_calibration_truth_v1.json"
OUTPUT = ARTIFACTS / "latent_lattice_metric_relation_search_v1.json"
HELD_LABEL = "rank_at_least_25"
SEARCH_BOUNDS = {
    "source_center_limit": 20,
    "center_pair_limit": 256,
    "seed_edges_per_center": 4,
    "initial_states_per_center_pair": 64,
    "minimum_states_per_center_pair": 2,
    "beam_width": 500,
    "maximum_steps": 80,
    "maximum_exact_lift_attempts": 500,
    "maximum_embeddings": 1,
    "minimum_global_replay_rays": 100,
    "reseed_after_skips": 20,
    "reseed_source_limit": 16,
    "reseed_target_limit": 4,
    "reseed_state_limit": 32,
    "norm_log_tolerance": 0.45,
    "angle_tolerance": 0.16,
    "angle_hard_tolerance": 0.36,
}


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    consensus = json.loads(CONSENSUS.read_text())
    truth = json.loads(TRUTH.read_text())
    held = next(
        record
        for record in consensus["leave_one_fibre_out"]
        if record["held_out_label"] == HELD_LABEL
    )
    control = next(
        record for record in consensus["controls"] if record["label"] == HELD_LABEL
    )
    source = build_relation_complex(held["training_two_of_three_vectors"])
    target = build_relation_complex(
        tuple(
            tuple(map(int, row.split()))
            for row in control["ambient_short_vector_coordinate_rows"]
        )
    )

    training_forms = []
    for record in truth["positive_controls"]:
        if record["label"] == HELD_LABEL:
            continue
        form = np.asarray(record["canonical_height_gram"], dtype=float)
        training_forms.append(form / np.trace(form))
    inferred_form = sum(training_forms) / len(training_forms)
    module = importlib.import_module("elkies_rank25")
    target_form = height_gram(
        EllipticCurve(tuple(module.GENERAL_WEIERSTRASS_COEFFICIENTS)),
        tuple(module.POINTS),
        digits=60,
    )
    ledger = bounded_metric_relation_search(
        source,
        target,
        inferred_form,
        target_form,
        target_center_limit=len(target.vertices),
        **SEARCH_BOUNDS,
    )
    ledger_record = {
        "source_center_count": ledger.source_center_count,
        "target_center_count": ledger.target_center_count,
        "initial_state_count": ledger.initial_state_count,
        "maximum_beam_population": ledger.maximum_beam_population,
        "expanded_state_count": ledger.expanded_state_count,
        "exact_lift_attempt_count": ledger.exact_lift_attempt_count,
        "maximum_source_rank_reached": ledger.maximum_source_rank_reached,
        "maximum_mapped_vertex_count": ledger.maximum_mapped_vertex_count,
        "maximum_matched_relation_count": ledger.maximum_matched_relation_count,
        "maximum_global_replay_ray_count": ledger.maximum_global_replay_ray_count,
        "exact_lifts_below_global_support": ledger.exact_lifts_below_global_support,
        "maximum_metric_reseed_count": ledger.maximum_metric_reseed_count,
        "accepted_embedding_count": len(ledger.embeddings),
    }
    passed = bool(ledger.embeddings)
    library_sources = tuple(sorted((ELLIPTIC / "latent_lattice").glob("*.py")))
    payload = {
        "schema": "elliptic-curves.latent-lattice-metric-relation-search.v1",
        "status": (
            "PASS_BLIND_R17_RECOVERY" if passed
            else "FAIL_BLIND_R17_RECOVERY_GATE_CLOSED"
        ),
        "scope": "Phase-0 rank-25 R17 held-out control only; no wgxli target is loaded",
        "held_out_label": HELD_LABEL,
        "source_training_core": {
            "ray_count": len(source.vertices),
            "rank": 17,
            "ternary_relation_count": len(source.ternary_relations),
        },
        "target_full_cloud": {
            "ray_count": len(target.vertices),
            "ambient_rank": 25,
            "ternary_relation_count": len(target.ternary_relations),
        },
        "search_bounds": SEARCH_BOUNDS,
        "search_ledger": ledger_record,
        "supervised_validator_reference": {
            "artifact": str(
                (ARTIFACTS / "latent_lattice_hypergraph_matcher_v1.json").relative_to(
                    ROOT
                )
            ),
            "truth_global_replay_ray_count": 238,
            "truth_global_replay_rank": 17,
        },
        "gate_decision": (
            "CLOSED. The exact validator sees a 238-ray primitive truth embedding, "
            "but this blind finite search reaches at most "
            f"{ledger.maximum_global_replay_ray_count} replayed rays and accepts none."
        ),
        "proof_boundary": (
            "Relation construction, integral lifts, rational ranks, and global ray "
            "replay are exact. Height forms and metric pruning are numerical. The "
            "238-ray reference is supervised post-search calibration evidence. This "
            "bounded failure proves inadequacy of this proposal generator, not "
            "nonexistence of a common R17 subgroup."
        ),
        "inputs": {
            str(path.relative_to(ROOT)): digest(path)
            for path in (
                CONSENSUS,
                TRUTH,
                ELLIPTIC / "cas/elkies_rank25.py",
                *library_sources,
                Path(__file__).resolve(),
            )
        },
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
        },
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != rendered:
            raise SystemExit("latent-lattice metric-relation artifact is stale")
        print(
            "LATENTMETRIC|check=PASS|"
            f"sha256={sha256(rendered.encode()).hexdigest()}"
        )
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(
        f"LATENTMETRIC|status={payload['status']}|output={args.output}|"
        f"sha256={sha256(rendered.encode()).hexdigest()}"
    )


if __name__ == "__main__":
    main()
