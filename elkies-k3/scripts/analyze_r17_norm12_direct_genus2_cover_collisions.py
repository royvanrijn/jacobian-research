#!/usr/bin/env python3
"""Intersect exact mod-p norm-six cover-collision sieves.

The input artifacts enumerate all affine two-node residual-chord slopes for
each norm-six trace that has good reduction in the displayed chart.  This
checker retains the scalar squareclass and intersects both

* norm-six versus complete smooth-atlas cover pairs; and
* pairs of distinct norm-six trace classes sharing a cover.

An empty intersection is an exact simultaneous-good-reduction obstruction,
not a characteristic-zero nonexistence theorem: a rational solution whose
parameters or trace chart have bad reduction at one of the selected primes is
outside the conclusion.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PRIMES = (23, 29, 31, 17, 37)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-11952-genus2-cover-collision-intersection-v1.json"
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
    inverse = pow(values[2], -1, prime)
    atom = tuple(value * inverse % prime for value in values)
    scalar_character = pow(values[2], (prime - 1) // 2, prime)
    return atom, scalar_character


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--input",
    type=Path,
    action="append",
    help="norm-six mod-p artifact; defaults to the five committed full screens",
)
parser.add_argument("--source-label", default="norm12-orbit-11952")
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

if args.input:
    input_paths = [path.resolve() for path in args.input]
else:
    input_paths = [
        ROOT
        / "artifacts/generated-results"
        / f"elkies-k3-r17-norm12-11952-genus2-normalization-full-p{prime}-direct-v1.json"
        for prime in DEFAULT_PRIMES
    ]

smooth_common = None
singular_common = None
prime_records = []
for path in input_paths:
    payload = json.loads(path.read_text())
    if payload.get("schema") != "elkies-k3.r17-norm12-direct-genus2-normalization-modp-search.v1":
        raise ValueError(f"unexpected schema: {path}")
    prime = int(payload["prime"])
    if prime < 3:
        raise ValueError("expected an odd prime")

    smooth_pairs = set()
    traces_by_cover = {}
    for survivor in payload["survivors"]:
        trace_mask = int(survivor["translation_orbit_mask"])
        for smooth_mask in survivor["smooth_cover_match_masks"]:
            smooth_mask = int(smooth_mask)
            if smooth_mask != trace_mask:
                smooth_pairs.add((trace_mask, smooth_mask))
        key = cover_key(
            survivor["reduced_quadratic_coefficients_low_to_high"], prime
        )
        if key is not None:
            traces_by_cover.setdefault(key, set()).add(trace_mask)

    singular_pairs = set()
    for trace_masks in traces_by_cover.values():
        singular_pairs.update(combinations(sorted(trace_masks), 2))

    smooth_common = (
        smooth_pairs if smooth_common is None else smooth_common & smooth_pairs
    )
    singular_common = (
        singular_pairs
        if singular_common is None
        else singular_common & singular_pairs
    )
    prime_records.append(
        {
            "prime": prime,
            "processed_trace_count": int(payload["search"]["processed_trace_count"]),
            "survivor_count": int(payload["survivor_count"]),
            "smooth_pair_count_at_prime": len(smooth_pairs),
            "distinct_norm_six_pair_count_at_prime": len(singular_pairs),
            "smooth_pair_intersection_count_through_prime": len(smooth_common),
            "distinct_norm_six_pair_intersection_count_through_prime": len(
                singular_common
            ),
        }
    )

smooth_common = smooth_common or set()
singular_common = singular_common or set()
output_path = args.output if args.output.is_absolute() else ROOT / args.output
output_path.parent.mkdir(parents=True, exist_ok=True)
result = {
    "schema": "elkies-k3.r17-norm12-direct-genus2-cover-collision-intersection.v1",
    "status": (
        "PASS_SIMULTANEOUS_GOOD_REDUCTION_NO_COVER_COLLISION"
        if not smooth_common and not singular_common
        else "SURVIVING_MODULAR_COVER_COLLISIONS"
    ),
    "source_label": args.source_label,
    "primes_in_intersection_order": [record["prime"] for record in prime_records],
    "prime_records": prime_records,
    "surviving_smooth_vs_norm_six_pairs": [
        {
            "norm_six_trace_mask": trace_mask,
            "smooth_bisection_mask": smooth_mask,
        }
        for trace_mask, smooth_mask in sorted(smooth_common)
    ],
    "surviving_distinct_norm_six_trace_pairs": [
        {"trace_masks": [left, right]}
        for left, right in sorted(singular_common)
    ],
    "proof_boundary": (
        "Every pair in the intersection has the same full finite-field cover "
        "squareclass, including its scalar character, at every displayed prime. "
        "An empty intersection excludes a characteristic-zero collision only when "
        "both traces, their affine slope parameters, and the smooth cover when "
        "applicable have simultaneous good reduction in all displayed charts. "
        "Bad-reduction and parameter-at-infinity cases remain open."
    ),
    "inputs": {relative(path): digest(path) for path in input_paths},
    "reproducing_command": (
        "python3 elkies-k3/scripts/analyze_r17_norm12_direct_genus2_cover_collisions.py "
        f"--source-label {args.source_label} "
        + " ".join(f"--input {relative(path)}" for path in input_paths)
        + f" --output {relative(output_path)}"
    ),
}
output_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
print(
    "R17GENUS2COVERINTERSECTION"
    f"|primes={','.join(map(str, result['primes_in_intersection_order']))}"
    f"|smooth_pairs={len(smooth_common)}"
    f"|norm6_pairs={len(singular_common)}"
    f"|output={relative(output_path)}|status={result['status']}",
    flush=True,
)
