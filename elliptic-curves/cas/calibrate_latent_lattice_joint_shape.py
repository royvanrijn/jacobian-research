#!/usr/bin/env python3
"""Audit primitive Hermite extremality on the blind R17 candidate ledger.

This is a Phase-0 control calculation only.  Candidate embeddings and ambient
height Grams come from the truth-free joint ledger.  Published R17 indices are
read only after all primitive-closure signatures and rankings are fixed.
"""

from __future__ import annotations

import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path
import platform
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC = ROOT / "elliptic-curves"
sys.path.insert(0, str(ELLIPTIC))

from latent_lattice import primitive_hermite_signatures  # noqa: E402


ARTIFACTS = ROOT / "artifacts/generated-results/elliptic-curves"
INPUT_LEDGER = ARTIFACTS / "latent_lattice_joint_fingerprint_ledger_v1.json.gz"
TRUTH_DIAGNOSTIC = ARTIFACTS / "latent_lattice_joint_fingerprints_v1.json"
OUTPUT = ARTIFACTS / "latent_lattice_joint_shape_v1.json"
SHAPE_LEDGER = ARTIFACTS / "latent_lattice_joint_shape_ledger_v1.json.gz"
BOUNDS = {
    "digits": 80,
    "maximum_vectors": 100_000,
    "batch_size": 64,
    "top_diagnostics": 16,
}


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--shape-ledger", type=Path, default=SHAPE_LEDGER)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    blind = json.loads(gzip.decompress(INPUT_LEDGER.read_bytes()))
    truth_document = json.loads(TRUTH_DIAGNOSTIC.read_text())
    truth_indices = {
        record["label"]: int(record["truth_source_index"])
        for record in truth_document["controls"]
    }
    shape_fibres = []
    controls = []
    exact_maximum_count = 0
    for fibre in blind["fibres"]:
        label = fibre["label"]
        signatures = primitive_hermite_signatures(
            fibre["ambient_height_gram"],
            [candidate["basis_rows"] for candidate in fibre["candidates"]],
            digits=BOUNDS["digits"],
            maximum_vectors=BOUNDS["maximum_vectors"],
            batch_size=BOUNDS["batch_size"],
            timeout=600,
        )
        ranked = sorted(
            range(len(signatures)),
            key=lambda index: (
                -float(signatures[index].hermite.log_hermite_invariant),
                index,
            ),
        )
        truth_index = truth_indices[label]
        truth_rank = ranked.index(truth_index)
        exact_maximum_count += int(truth_rank == 0)
        controls.append(
            {
                "label": label,
                "candidate_count": len(signatures),
                "primitive_candidate_count": sum(
                    signature.saturation_index == 1 for signature in signatures
                ),
                "truth_source_index": truth_index,
                "truth_descending_hermite_rank": truth_rank,
                "truth_signature": signatures[truth_index].to_record(),
                "selected_source_index": ranked[0],
                "selected_signature": signatures[ranked[0]].to_record(),
                "top_candidates": [
                    {
                        "source_index": index,
                        **signatures[index].to_record(),
                    }
                    for index in ranked[: BOUNDS["top_diagnostics"]]
                ],
            }
        )
        shape_fibres.append(
            {
                "label": label,
                "candidates": [
                    {"source_index": index, **signature.to_record()}
                    for index, signature in enumerate(signatures)
                ],
            }
        )
        print(
            f"LATENTSHAPEPROGRESS|label={label}|candidates={len(signatures)}|"
            f"truth_rank={truth_rank}",
            flush=True,
        )

    shape_payload = {
        "schema": "elliptic-curves.latent-lattice-joint-shape-ledger.v1",
        "scope": "Blind R17 control candidates only; no wgxli target is loaded",
        "bounds": BOUNDS,
        "fibres": shape_fibres,
    }
    shape_rendered = (
        json.dumps(shape_payload, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode()
    shape_bytes = gzip.compress(shape_rendered, compresslevel=9, mtime=0)
    shape_hash = sha256(shape_bytes).hexdigest()
    status = (
        "PASS_PRIMITIVE_HERMITE_EXTREMAL_SELECTOR"
        if exact_maximum_count == len(controls)
        else "FAIL_PRIMITIVE_HERMITE_EXTREMAL_SELECTOR_GATE_CLOSED"
    )
    library_sources = tuple(sorted((ELLIPTIC / "latent_lattice").glob("*.py")))
    payload = {
        "schema": "elliptic-curves.latent-lattice-joint-shape.v1",
        "status": status,
        "scope": "Phase-0 R17 rank-25--28 controls only; no wgxli target is loaded",
        "algorithm": {
            "exact_step": (
                "Smith index plus primitive row closure by double integer kernel "
                "for every candidate"
            ),
            "numerical_step": (
                "80-digit restricted Gram; PARI LLL and complete minimum search "
                "through the shortest reduced-basis diagonal"
            ),
            "blind_selector": "maximum scale-free Hermite invariant per fibre",
            "bounds": BOUNDS,
        },
        "controls": controls,
        "exact_maximum_selection_count": exact_maximum_count,
        "shape_ledger": {
            "path": str(args.shape_ledger.relative_to(ROOT)),
            "sha256": shape_hash,
            "compressed_bytes": len(shape_bytes),
            "uncompressed_bytes": len(shape_rendered),
        },
        "gate_decision": (
            "OPEN. The primitive-Hermite selector chooses all withheld controls."
            if exact_maximum_count == len(controls)
            else "CLOSED. Hermite extremality narrows the controls but is not an "
            f"exact selector ({exact_maximum_count}/{len(controls)})."
        ),
        "proof_boundary": (
            "Candidate matrices, Smith indices, primitive closures, and candidate "
            "identities are exact. Height Grams, minima, determinants, Hermite "
            "invariants, and their ordering are numerical at the declared precision. "
            "Published truth indices are postselection diagnostics only."
        ),
        "inputs": {
            str(path.relative_to(ROOT)): digest(path)
            for path in (
                INPUT_LEDGER,
                TRUTH_DIAGNOSTIC,
                *library_sources,
                Path(__file__).resolve(),
            )
        },
        "software": {"python": platform.python_version(), "numpy": np.__version__},
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if (
            not args.output.exists()
            or args.output.read_text() != rendered
            or not args.shape_ledger.exists()
            or args.shape_ledger.read_bytes() != shape_bytes
        ):
            raise SystemExit("latent-lattice joint-shape artifact is stale")
        print(f"LATENTSHAPE|check=PASS|sha256={sha256(rendered.encode()).hexdigest()}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.shape_ledger.parent.mkdir(parents=True, exist_ok=True)
    args.shape_ledger.write_bytes(shape_bytes)
    args.output.write_text(rendered)
    print(
        f"LATENTSHAPE|status={status}|output={args.output}|"
        f"sha256={sha256(rendered.encode()).hexdigest()}"
    )


if __name__ == "__main__":
    main()
