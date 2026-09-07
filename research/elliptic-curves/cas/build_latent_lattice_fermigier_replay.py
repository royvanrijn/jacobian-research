#!/usr/bin/env python3
"""Build the truth-free ICARM-245 rank-12 cross-bound replay ledger."""

from __future__ import annotations

import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path
import platform
import sys


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC = ROOT / "elliptic-curves"
sys.path[:0] = [str(ELLIPTIC), str(ELLIPTIC / "cas")]

from icarm_curve245 import (  # noqa: E402
    GENERAL_WEIERSTRASS_COEFFICIENTS,
    POINTS,
)
from latent_lattice import (  # noqa: E402
    EllipticCurve,
    FiniteQuotientBlock,
    build_relation_complex,
    candidate_finite_signature,
    enumerate_short_vectors,
    height_gram,
    primitive_hermite_signatures,
    recombined_core_extension_search,
    repeated_cross_bound_intersection_ledger,
)


ARTIFACTS = ROOT / "artifacts/generated-results/elliptic-curves"
FINITE = ARTIFACTS / "latent_lattice_finite_calibration_v1.json"
OUTPUT = ARTIFACTS / "latent_lattice_fermigier_replay_v1.json.gz"
BOUNDS = {
    "dimension": 12,
    "height_bounds": [28.0, 29.0],
    "height_digits": 80,
    "seed_edges": 3_000,
    "anchor_count": 500,
    "enclosure_codimension": 3,
    "left_count": 200,
    "right_count": 200,
    "maximum_candidates": 128,
    "hermite_maximum_vectors": 100_000,
}


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def finite_block(record: dict[str, object]) -> FiniteQuotientBlock:
    return FiniteQuotientBlock(
        reduction_prime=int(record["reduction_prime"]),
        relation_prime=int(record["relation_prime"]),
        group_order=int(record["group_order"]),
        multiple_subgroup_order=int(record["multiple_subgroup_order"]),
        quotient_dimension=int(record["quotient_dimension"]),
        rows=tuple(tuple(map(int, row)) for row in record["rows"]),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    finite_document = json.loads(FINITE.read_text())
    finite_control = next(
        record
        for record in finite_document["controls"]
        if record["label"] == "ICARM_245_Fermigier_negative_control"
    )
    development = tuple(
        finite_block(record) for record in finite_control["development_blocks"]
    )
    held_out = tuple(
        finite_block(record) for record in finite_control["held_out_blocks"]
    )
    curve = EllipticCurve(tuple(GENERAL_WEIERSTRASS_COEFFICIENTS))
    ledgers = {}
    states = {}
    for bound in BOUNDS["height_bounds"]:
        records = enumerate_short_vectors(
            curve,
            POINTS,
            height_bound=bound,
            digits=BOUNDS["height_digits"],
            maximum_lines=100_000,
            materialize_points=True,
        )
        complex_ = build_relation_complex([record.coordinates for record in records])
        ledgers[bound] = recombined_core_extension_search(
            records,
            complex_,
            dimension=BOUNDS["dimension"],
            seed_edges=BOUNDS["seed_edges"],
            anchor_count=BOUNDS["anchor_count"],
            enclosure_codimension=BOUNDS["enclosure_codimension"],
            enclosure_count=0,
            inner_count=2,
        )
        states[bound] = (records, complex_)
        print(
            f"FERMIGIERREPLAYPROGRESS|bound={bound}|rays={len(records)}|"
            f"enclosures={len(ledgers[bound].enclosure_proposals)}",
            flush=True,
        )
    records, complex_ = states[28.0]
    repeated = repeated_cross_bound_intersection_ledger(
        records,
        complex_,
        ledgers[28.0].enclosure_proposals,
        ledgers[29.0].enclosure_proposals,
        target_dimension=BOUNDS["dimension"],
        left_count=BOUNDS["left_count"],
        right_count=BOUNDS["right_count"],
        maximum_candidates=BOUNDS["maximum_candidates"],
    )
    ambient_gram = height_gram(curve, POINTS, digits=BOUNDS["height_digits"])
    shape = primitive_hermite_signatures(
        ambient_gram,
        [proposal.primitive_basis_rows for proposal in repeated.proposals],
        digits=BOUNDS["height_digits"],
        maximum_vectors=BOUNDS["hermite_maximum_vectors"],
        batch_size=64,
        timeout=600,
    )
    candidates = []
    for index, (proposal, signature) in enumerate(zip(repeated.proposals, shape)):
        development_signature = candidate_finite_signature(
            proposal.primitive_basis_rows,
            complex_,
            finite_blocks=development,
        )
        held_out_signature = candidate_finite_signature(
            proposal.primitive_basis_rows,
            complex_,
            finite_blocks=held_out,
        )
        candidates.append(
            {
                "source_index": index,
                "proposal": proposal.to_record(),
                "primitive_hermite_signature": signature.to_record(),
                "development_finite_signature": development_signature.to_record(),
                "held_out_finite_signature": held_out_signature.to_record(),
            }
        )
    library_sources = tuple(sorted((ELLIPTIC / "latent_lattice").glob("*.py")))
    payload = {
        "schema": "elliptic-curves.latent-lattice-fermigier-replay.v1",
        "scope": "Truth-free Phase-0 ICARM-245 control; no wgxli target is loaded",
        "bounds": BOUNDS,
        "short_vector_counts": {
            str(bound): len(states[bound][0]) for bound in BOUNDS["height_bounds"]
        },
        "enclosure_counts": {
            str(bound): len(ledgers[bound].enclosure_proposals)
            for bound in BOUNDS["height_bounds"]
        },
        "repeated_intersection_summary": repeated.summary_record(),
        "ambient_height_gram": ambient_gram,
        "candidates": candidates,
        "proof_boundary": (
            "Candidate coordinates, cross-bound rational intersections, primitive "
            "closures, finite-code restrictions, and relation counts are exact. "
            "Canonical heights, proposal ordering, and Hermite signatures are "
            "numerical at the declared precision. This ledger contains no truth basis."
        ),
        "inputs": {
            str(path.relative_to(ROOT)): digest(path)
            for path in (
                FINITE,
                ELLIPTIC / "cas/icarm_curve245.py",
                *library_sources,
                Path(__file__).resolve(),
            )
        },
        "software": {"python": platform.python_version()},
    }
    rendered = (
        json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode()
    compressed = gzip.compress(rendered, compresslevel=9, mtime=0)
    if args.check:
        if not args.output.exists() or args.output.read_bytes() != compressed:
            raise SystemExit("latent-lattice Fermigier replay artifact is stale")
        print(f"FERMIGIERREPLAY|check=PASS|sha256={sha256(compressed).hexdigest()}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(compressed)
    print(
        f"FERMIGIERREPLAY|candidates={len(candidates)}|output={args.output}|"
        f"sha256={sha256(compressed).hexdigest()}"
    )


if __name__ == "__main__":
    main()
