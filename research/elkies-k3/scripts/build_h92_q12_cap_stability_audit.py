#!/usr/bin/env python3
"""Compare two bounded q12 neighbor shells by exact rootless fibre identity.

status: ACTIVE_SEARCH
claim: exact stability comparison of the retained rootless sets at two MW caps
inputs: two root-adapted Weyl-neighbor JSON artifacts
outputs: caller-selected compact generated audit
"""

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--smaller", type=Path, required=True)
parser.add_argument("--larger", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path):
    path = path.resolve()
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def rootless_fibres(data):
    return sorted(
        tuple(item["fiber"])
        for item in data["neighbors"]
        if item["child_root_data"] == [0, 0, 1] and item["child_mw_rank"] == 17
    )


paths = (args.smaller.resolve(), args.larger.resolve())
smaller, larger = (json.loads(path.read_text()) for path in paths)
small_summary, large_summary = smaller["summaries"][0], larger["summaries"][0]
assert small_summary["q"] == large_summary["q"] == 12
assert small_summary["factor_order"] == large_summary["factor_order"] == [6, 2]
assert small_summary["mw_vector_cap"] < large_summary["mw_vector_cap"]
small_set, large_set = rootless_fibres(smaller), rootless_fibres(larger)
assert len(small_set) == len(set(small_set))
assert len(large_set) == len(set(large_set))

payload = {
    "schema": "elkies-k3.h3-q12-cap-stability-audit.v1",
    "status": "PASS_EXACT_BOUNDED_Q12_ROOTLESS_SET_STABLE_ACROSS_DOUBLED_CAP",
    "q": 12,
    "old_fibre_degree": 2,
    "smaller_cap": small_summary["mw_vector_cap"],
    "larger_cap": large_summary["mw_vector_cap"],
    "full_mw_vector_count": large_summary["mw_pari_vector_count"],
    "smaller_primitive_candidate_count": small_summary["primitive_neighbors"],
    "larger_primitive_candidate_count": large_summary["primitive_neighbors"],
    "smaller_rootless_count": len(small_set),
    "larger_rootless_count": len(large_set),
    "rootless_fibre_sets_equal": small_set == large_set,
    "new_rootless_fibres_in_larger_sample": [list(item) for item in sorted(set(large_set) - set(small_set))],
    "stable_rootless_fibres": [list(item) for item in small_set],
    "proof_boundary": (
        "The equality is exact for the two stored bounded samples. The larger cap covers "
        f"{large_summary['mw_vector_cap']} of {large_summary['mw_pari_vector_count']} MW vectors; "
        "this is not a complete q12-shell or global optimality theorem."
    ),
    "inputs": {
        "paths": [relative(path) for path in paths],
        "sha256": {relative(path): sha256(path) for path in paths},
    },
}
assert payload["rootless_fibre_sets_equal"]
output = args.output.resolve()
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(output)
