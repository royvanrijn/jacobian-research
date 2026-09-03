#!/usr/bin/env python3
"""Apply targeted extra-prime screens to surviving norm-four cover pairs."""

from __future__ import annotations

import argparse
from collections import defaultdict
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
DEFAULT_CANDIDATES = (
    GENERATED
    / "elkies-k3-r17-norm12-11952-genus3-cover-collision-intersection-v1.json"
)
DEFAULT_PRIMES = (43, 47, 53)
DEFAULT_OUTPUT = (
    GENERATED
    / "elkies-k3-r17-norm12-11952-genus3-cover-collision-targeted-v1.json"
)


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def cover_key(coefficients: list[int], prime: int):
    values = [int(value) % prime for value in coefficients]
    if len(values) != 3 or values[2] == 0:
        return None
    if (values[1] * values[1] - 4 * values[0] * values[2]) % prime == 0:
        return None
    inverse = pow(values[2], -1, prime)
    return (
        tuple(value * inverse % prime for value in values),
        pow(values[2], (prime - 1) // 2, prime),
    )


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--candidates", type=Path, default=DEFAULT_CANDIDATES)
parser.add_argument("--screen", type=Path, action="append")
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

candidate_path = args.candidates.resolve()
candidate_payload = json.loads(candidate_path.read_text())
if candidate_payload.get("schema") != (
    "elkies-k3.r17-norm12-direct-genus3-cover-collision-intersection.v1"
):
    raise ValueError("unexpected candidate schema")
pairs = {
    tuple(map(int, record["trace_masks"]))
    for record in candidate_payload["surviving_distinct_norm_four_trace_pairs"]
}

screen_paths = (
    [path.resolve() for path in args.screen]
    if args.screen
    else [
        GENERATED
        / f"elkies-k3-r17-norm12-11952-genus3-normalization-four-pair-traces-p{prime}-v1.json"
        for prime in DEFAULT_PRIMES
    ]
)
common = set(pairs)
prime_records = []
for path in screen_paths:
    payload = json.loads(path.read_text())
    if payload.get("schema") != (
        "elkies-k3.r17-norm12-direct-genus3-normalization-modp-search.v1"
    ):
        raise ValueError(f"unexpected screen schema: {path}")
    prime = int(payload["prime"])
    traces_by_key = defaultdict(set)
    for survivor in payload["survivors"]:
        key = cover_key(
            survivor["reduced_quadratic_coefficients_low_to_high"], prime
        )
        if key is not None:
            traces_by_key[key].add(int(survivor["translation_orbit_mask"]))
    present = {
        pair
        for pair in pairs
        if any(pair[0] in masks and pair[1] in masks for masks in traces_by_key.values())
    }
    common &= present
    prime_records.append(
        {
            "prime": prime,
            "processed_trace_count": int(payload["search"]["processed_trace_count"]),
            "input_candidate_pair_count": len(pairs),
            "pair_count_at_prime": len(present),
            "intersection_count_through_prime": len(common),
            "surviving_pairs_at_prime": [list(pair) for pair in sorted(present)],
        }
    )

output_path = args.output if args.output.is_absolute() else ROOT / args.output
output_path.parent.mkdir(parents=True, exist_ok=True)
inputs = [candidate_path] + screen_paths
result = {
    "schema": "elkies-k3.r17-norm12-direct-genus3-targeted-cover-collision.v1",
    "status": (
        "PASS_TARGETED_GOOD_REDUCTION_NO_COVER_COLLISION"
        if not common
        else "SURVIVING_TARGETED_MODULAR_COVER_COLLISIONS"
    ),
    "source_label": "norm12-orbit-11952",
    "input_candidate_pairs": [list(pair) for pair in sorted(pairs)],
    "prime_records": prime_records,
    "surviving_pairs": [list(pair) for pair in sorted(common)],
    "proof_boundary": (
        "The input pairs are exactly the survivors of the five-prime full "
        "affine census. An absent pair cannot be a characteristic-zero collision "
        "with simultaneous good integral reduction in the displayed charts at "
        "that prime. Parameter denominators and boundary charts prevent this "
        "targeted no-hit from being a characteristic-zero nonexistence theorem."
    ),
    "inputs": {relative(path): digest(path) for path in inputs},
    "reproducing_command": (
        "python3 "
        "elkies-k3/scripts/analyze_r17_norm12_direct_genus3_targeted_collisions.py"
    ),
}
output_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
print(
    "R17GENUS3TARGETED"
    f"|primes={','.join(str(record['prime']) for record in prime_records)}"
    f"|input_pairs={len(pairs)}|surviving_pairs={len(common)}"
    f"|output={relative(output_path)}|status={result['status']}",
    flush=True,
)
