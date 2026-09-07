#!/usr/bin/env python3
"""Export root-adapted child frames embedded in a neighbor-search artifact.

status: ACTIVE_SEARCH
claim: exact mechanical extraction with source hashes; no new lattice claim
inputs: a root-adapted neighbor-search artifact containing child frames
outputs: child frame files and a generated manifest
"""

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, required=True)
parser.add_argument("--child-mw-rank", type=int, required=True)
parser.add_argument("--output-prefix", type=Path, required=True)
parser.add_argument("--manifest", type=Path, required=True)
args = parser.parse_args()

source = args.input.resolve()
prefix = args.output_prefix.resolve()
manifest = args.manifest.resolve()
payload = json.loads(source.read_text())
assert payload["status"] in {
    "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS",
    "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS_TARGET_FILTERED",
}

outputs = []
for record in payload["neighbors"]:
    if int(record["child_mw_rank"]) != args.child_mw_rank:
        continue
    frame = record.get("child_root_adapted_frame")
    if frame is None:
        continue
    candidate = {
        "q": int(record["q"]),
        "old_fibre_degree": int(record["old_fiber_degree"]),
        "orbit_index": int(record["orbit_index"]),
    }
    path = Path(
        f"{prefix}-q{candidate['q']}o{candidate['orbit_index']}-frame.txt"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "# exact root-adapted child frame exported from neighbor search\n"
        + "\n".join(" ".join(map(str, row)) for row in frame)
        + "\n"
    )
    outputs.append({
        "candidate_id": candidate,
        "child_root_data": record["child_root_data"],
        "child_mw_rank": int(record["child_mw_rank"]),
        "frame": str(path.relative_to(ROOT)),
        "frame_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "source_neighbor_record": record,
    })

result = {
    "schema": "elkies-k3.h3-neighbor-child-frame-export.v1",
    "status": "PASS_EXACT_NEIGHBOR_CHILD_FRAME_EXPORT",
    "child_mw_rank": args.child_mw_rank,
    "count": len(outputs),
    "outputs": outputs,
    "input": str(source.relative_to(ROOT)),
    "input_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
}
manifest.parent.mkdir(parents=True, exist_ok=True)
manifest.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
print(f"CHILDFRAMES|MW={args.child_mw_rank}|count={len(outputs)}|output={manifest}")
