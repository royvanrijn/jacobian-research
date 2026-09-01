#!/usr/bin/env sage
"""Build the contiguous, hash-pinned N(24A1) octad-completion manifest.

By default this discovers every generated completion shard with schema v1,
sorts by its declared half-open prefix range, proves that the ranges form one
contiguous interval starting at zero, and pins every artifact hash and exact
local accounting value.  The manifest is the shared input contract for the
multi-shard Weyl-M24 canonicalizer and the surface-first catalogue.

status: EXACT_CONTIGUOUS_SHARD_MANIFEST
"""

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
GENERATED = ROOT / "artifacts/generated-results"
PATTERN = "elkies-k3-24a1-octad-rank7-completion-*-v1.json"
DEFAULT_OUTPUT = GENERATED / "elkies-k3-24a1-octad-completion-manifest-v1.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build(paths):
    loaded = []
    for path in paths:
        payload = json.loads(path.read_text())
        assert payload["schema"] == (
            "elkies-k3.24a1-octad-rank7-completion-shard.v1"
        )
        assert payload["status"] == (
            "PASS_EXACT_DECLARED_24A1_OCTAD_COMPLETION_SHARD"
        )
        parameters = payload["parameters"]
        loaded.append(
            (
                parameters["prefix_start_zero_based_inclusive"],
                parameters["prefix_stop_zero_based_exclusive"],
                path,
                payload,
            )
        )
    loaded.sort(key=lambda item: (item[0], item[1]))
    assert loaded and loaded[0][0] == 0
    assert all(
        left[1] == right[0] for left, right in zip(loaded, loaded[1:])
    )
    determinant_bounds = {
        payload["parameters"]["determinant_bound"]
        for unused_start, unused_stop, unused_path, payload in loaded
    }
    prefix_hashes = {
        payload["input"]["prefix_artifact_sha256"]
        for unused_start, unused_stop, unused_path, payload in loaded
    }
    assert len(determinant_bounds) == len(prefix_hashes) == 1

    determinant_distribution = Counter()
    mw_distribution = Counter()
    closure_distribution = Counter()
    local_records = 0
    compatible_records = 0
    records = []
    for index, (start, stop, path, payload) in enumerate(loaded):
        local_records += len(payload["orbits"])
        compatible = sum(
            record["k3_discriminant_gate"]["matching_even_ternary_genera"] > 0
            for record in payload["orbits"]
        )
        compatible_records += compatible
        determinant_distribution.update(
            record["determinant"] for record in payload["orbits"]
        )
        mw_distribution.update(
            record["mordell_weil_rank"] for record in payload["orbits"]
        )
        closure_distribution.update(
            record["primitive_closure_index"] for record in payload["orbits"]
        )
        records.append(
            {
                "shard_index_zero_based": index,
                "artifact": str(path.relative_to(ROOT)),
                "sha256": digest(path),
                "prefix_start_zero_based_inclusive": start,
                "prefix_stop_zero_based_exclusive": stop,
                "prefixes_processed": stop - start,
                "shard_local_residual_m24_records": len(payload["orbits"]),
                "shard_local_k3_compatible_genus_records": compatible,
            }
        )
    return {
        "schema": "elkies-k3.24a1-octad-completion-manifest.v1",
        "status": "PASS_EXACT_CONTIGUOUS_24A1_OCTAD_COMPLETION_SHARD_MANIFEST",
        "proof_scope": {
            "proved": (
                "Every listed completion artifact has the required exact shard "
                "schema/status and pinned hash; its declared prefix ranges form "
                "one gap-free, overlap-free half-open interval starting at zero."
            ),
            "not_proved": (
                "The manifest covers only the listed positive-octad completion "
                "frontier and does not promote it to a complete 24A1 auxiliary census."
            ),
        },
        "parameters": {
            "prefix_start_zero_based_inclusive": loaded[0][0],
            "prefix_stop_zero_based_exclusive": loaded[-1][1],
            "determinant_bound": next(iter(determinant_bounds)),
            "prefix_artifact_sha256": next(iter(prefix_hashes)),
        },
        "accounting": {
            "completion_shards": len(records),
            "prefixes_processed": loaded[-1][1] - loaded[0][0],
            "shard_local_residual_m24_records": local_records,
            "shard_local_k3_compatible_genus_records": compatible_records,
            "primitive_closure_index_distribution": {
                str(key): value for key, value in sorted(closure_distribution.items())
            },
            "determinant_distribution": {
                str(key): value for key, value in sorted(determinant_distribution.items())
            },
            "mordell_weil_rank_distribution": {
                str(key): value for key, value in sorted(mw_distribution.items())
            },
        },
        "shards": records,
    }


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--shard", type=Path, action="append")
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()
paths = arguments.shard or sorted(GENERATED.glob(PATTERN))
result = build(paths)
encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if not arguments.output.exists() or arguments.output.read_text() != encoded:
        raise SystemExit("24A1 octad completion manifest is stale")
else:
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(encoded)
print(
    "OCTADMANIFEST|shards={}|prefixes={}:{}|records={}|status=PASS_EXACT".format(
        result["accounting"]["completion_shards"],
        result["parameters"]["prefix_start_zero_based_inclusive"],
        result["parameters"]["prefix_stop_zero_based_exclusive"],
        result["accounting"]["shard_local_residual_m24_records"],
    )
)
