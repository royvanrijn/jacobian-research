#!/usr/bin/env python3
"""Certify adjacent coprime-stride NS0030 fibre scans as one full census."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PREFIX = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0030-source-ansatz-mod5-pilot100k-v1.json"
DEFAULT_SUFFIX = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0030-source-ansatz-mod5-suffix290625-v1.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0030-source-ansatz-mod5-complete-v1.json"
ACCOUNTING_KEYS = (
    "normalized_A_samples",
    "branch_eligible_with_signs",
    "hermite_compatible_with_signs",
    "exact_prescribed_orders",
    "squarefree_examples_with_signs",
)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prefix", type=Path, default=DEFAULT_PREFIX)
parser.add_argument("--suffix", type=Path, default=DEFAULT_SUFFIX)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--examples", type=int, default=100)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()

input_paths = [args.prefix.resolve(), args.suffix.resolve()]
payloads = [json.loads(path.read_text()) for path in input_paths]
for payload in payloads:
    if payload["schema"] != "elkies-k3.lattice-foundry-ns0030-source-ansatz-modp.v1":
        raise ValueError("unexpected NS0030 source-ansatz schema")
for key in ("prime", "ansatz", "source"):
    if payloads[0][key] != payloads[1][key]:
        raise ValueError(f"NS0030 scan inputs disagree on {key}")

per_pair = int(payloads[0]["scan"]["normalized_A_polynomials_per_support_pair"])
prime = int(payloads[0]["prime"])
if per_pair != prime**8:
    raise ValueError("unexpected normalized A population")

records_by_payload = []
for payload in payloads:
    records = {
        (int(row["lambda"]), int(row["mu"])): row
        for row in payload["scan"]["support_pairs"]
    }
    records_by_payload.append(records)
if set(records_by_payload[0]) != set(records_by_payload[1]):
    raise ValueError("support-pair populations differ")

combined_records = []
for support_pair in sorted(records_by_payload[0]):
    prefix = records_by_payload[0][support_pair]
    suffix = records_by_payload[1][support_pair]
    stride = int(prefix["sample_stride"])
    if int(suffix["sample_stride"]) != stride or math.gcd(stride, per_pair) != 1:
        raise ValueError("stride is not a shared permutation of the coefficient population")
    prefix_count = int(prefix["samples_consumed"])
    suffix_count = int(suffix["samples_consumed"])
    prefix_offset = int(prefix["sample_offset"])
    suffix_offset = int(suffix["sample_offset"])
    if (prefix_offset + prefix_count * stride) % per_pair != suffix_offset:
        raise ValueError("suffix does not begin after the prefix along the stride permutation")
    if prefix_count + suffix_count != per_pair:
        raise ValueError("stride segments do not cover the complete coefficient population")
    accounting = {
        key: int(prefix["accounting"][key]) + int(suffix["accounting"][key])
        for key in ACCOUNTING_KEYS
    }
    combined_records.append(
        {
            "lambda": support_pair[0],
            "mu": support_pair[1],
            "samples_consumed": per_pair,
            "exhausted": True,
            "sample_stride": stride,
            "sample_offset": prefix_offset,
            "certified_adjacent_segments": [
                {
                    "artifact": relative(input_paths[0]),
                    "offset": prefix_offset,
                    "samples": prefix_count,
                },
                {
                    "artifact": relative(input_paths[1]),
                    "offset": suffix_offset,
                    "samples": suffix_count,
                },
            ],
            "accounting": accounting,
        }
    )

accounting = {
    key: sum(row["accounting"][key] for row in combined_records)
    for key in ACCOUNTING_KEYS
}
examples = (payloads[0]["examples"] + payloads[1]["examples"])[: args.examples]
accounting["stored_examples"] = len(examples)
has_examples = bool(accounting["squarefree_examples_with_signs"])
if has_examples != bool(examples):
    raise ValueError("stored examples do not witness the aggregate squarefree count")

output = {
    "schema": "elkies-k3.lattice-foundry-ns0030-source-ansatz-modp-combined.v1",
    "status": (
        "PASS_EXACT_EXHAUSTIVE_MODULAR_SOURCE_FIBRE_ANSATZ_WITH_EXAMPLES"
        if has_examples
        else "PASS_EXACT_EXHAUSTIVE_MODULAR_SOURCE_FIBRE_ANSATZ_EMPTY"
    ),
    "prime": prime,
    "inputs": [
        {
            "artifact": relative(path),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }
        for path in input_paths
    ],
    "scan": {
        "support_pairs": combined_records,
        "exhausted": True,
        "normalized_A_polynomials_per_support_pair": per_pair,
        "certified_disjoint_union": (
            "For each support pair the two segments are adjacent intervals in "
            "one stride permutation; gcd(stride,p^8)=1 and their lengths sum to p^8."
        ),
    },
    "accounting": accounting,
    "ansatz": payloads[0]["ansatz"],
    "examples": examples,
    "source": payloads[0]["source"],
    "proof_boundary": {
        "proved": (
            "The two certified stride segments exhaust every normalized A "
            "polynomial exactly once for every displayed ordered support pair. "
            + (
                "Every stored example has the exact prescribed finite-field fibre profile."
                if has_examples
                else "No branch in this complete normalized finite-field fibre chart has the exact prescribed orders and squarefree residual cubic."
            )
        ),
        "not_proved": (
            "A finite-field chart result is not characteristic-zero nonexistence, "
            "a rational source marking, or a physical neighbour corridor."
        ),
    },
    "reproduce": (
        "python3 elkies-k3/scripts/combine_lattice_foundry_ns0030_source_ansatz_modp.py"
    ),
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
output_path = args.output.resolve()
if args.check:
    if output_path.read_text() != serialized:
        raise SystemExit("combined NS0030 source-ansatz artifact is stale")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)

print(
    "FOUNDRYNS0030COMBINE|"
    f"pairs={len(combined_records)}|samples={accounting['normalized_A_samples']}|"
    f"compatible={accounting['hermite_compatible_with_signs']}|"
    f"squarefree={accounting['squarefree_examples_with_signs']}|status=PASS",
    flush=True,
)
print(f"OUTPUT|{output_path}", flush=True)
