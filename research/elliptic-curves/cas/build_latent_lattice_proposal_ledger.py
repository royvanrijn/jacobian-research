#!/usr/bin/env python3
"""Build a bounded, replayable latent-subspace proposal ledger for one fibre.

This is a high-recall candidate generator, not a selector and not a proof of a
generic Mordell--Weil subgroup.  It consumes only a public curve model and its
displayed independent points.  Withheld calibration embeddings must be scored
in a separate process after this artifact has been frozen.
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
    build_relation_complex,
    enumerate_short_vectors,
    modular_row_space_key,
    recombined_core_extension_search,
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--curve-module", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--height-bound", type=float, required=True)
    parser.add_argument("--dimension", type=int, required=True)
    parser.add_argument("--digits", type=int, default=80)
    parser.add_argument("--maximum-lines", type=int, default=100_000)
    parser.add_argument("--seed-edges", type=int, default=3_000)
    parser.add_argument("--anchor-count", type=int, default=20)
    parser.add_argument("--enclosure-count", type=int, default=20)
    parser.add_argument("--inner-count", type=int, default=2)
    parser.add_argument("--include-complex", action="store_true")
    parser.add_argument("--include-stage-proposals", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    module = importlib.import_module(args.curve_module)
    source = Path(module.__file__).resolve()
    model = tuple(module.GENERAL_WEIERSTRASS_COEFFICIENTS)
    points = tuple(module.POINTS)
    curve = EllipticCurve(model)
    records = enumerate_short_vectors(
        curve,
        points,
        height_bound=args.height_bound,
        digits=args.digits,
        maximum_lines=args.maximum_lines,
        materialize_points=True,
    )
    metadata = [
        {
            "integral": bool(record.arithmetic.get("integral", False)),
            "total_bits": int(record.arithmetic.get("total_bits", 0)),
            "x_denominator_bits": int(
                record.arithmetic.get("x_denominator_bits", 0)
            ),
            "y_denominator_bits": int(
                record.arithmetic.get("y_denominator_bits", 0)
            ),
        }
        for record in records
    ]
    complex_ = build_relation_complex(
        [record.coordinates for record in records], metadata
    )
    ledger = recombined_core_extension_search(
        records,
        complex_,
        dimension=args.dimension,
        seed_edges=args.seed_edges,
        anchor_count=args.anchor_count,
        enclosure_count=args.enclosure_count,
        inner_count=args.inner_count,
    )
    def proposal_records(proposals):
        answer = []
        for rank, proposal in enumerate(proposals):
            record = proposal.to_record()
            record["blind_rank"] = rank
            record["modular_row_space_keys"] = {
                str(prime): [
                    list(row)
                    for row in modular_row_space_key(proposal.basis_rows, prime)
                ]
                for prime in (1_000_003, 1_000_033)
            }
            answer.append(record)
        return answer

    refined = proposal_records(ledger.refined_proposals)
    stage_proposals = None
    if args.include_stage_proposals:
        stage_proposals = {
            "direct": proposal_records(ledger.direct_proposals),
            "enclosures": proposal_records(ledger.enclosure_proposals),
        }
    relation_record = complex_.to_record(include_relations=args.include_complex)
    if not args.include_complex:
        relation_record.pop("vertices", None)
        relation_record.pop("metadata", None)
    payload = {
        "schema": "elliptic-curves.latent-lattice-proposal-ledger.v1",
        "status": "PASS_BOUNDED_PROPOSAL_LEDGER",
        "role": (
            "High-recall blind proposal ledger. It contains no generic points, "
            "withheld embeddings, family labels, or truth-derived scores."
        ),
        "label": args.label,
        "public_input": {
            "curve_module": args.curve_module,
            "source": str(source.relative_to(ROOT)),
            "source_sha256": digest(source),
            "displayed_rank": len(points),
        },
        "enumeration": {
            "height_bound": args.height_bound,
            "digits": args.digits,
            "maximum_lines": args.maximum_lines,
            "primitive_unoriented_lines": len(records),
            "complete_within_bound": len(records) < args.maximum_lines,
        },
        "relation_complex": relation_record,
        "search": ledger.summary_record(),
        "stage_proposals": stage_proposals,
        "refined_proposals": refined,
        "proof_boundary": (
            "Point arithmetic and relation identities are exact. Canonical "
            "heights, arithmetic scores, numerical span masks, and proposal "
            "ordering are numerical. Two-prime row-space keys are only fast "
            "matching filters; selected survivors require exact rational-rank "
            "and primitive-closure replay."
        ),
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
        },
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(
        f"LATENTLEDGER|label={args.label}|dimension={args.dimension}|"
        f"lines={len(records)}|candidates={len(refined)}|"
        f"output={args.output}|sha256={sha256(rendered.encode()).hexdigest()}"
    )


if __name__ == "__main__":
    main()
