#!/usr/bin/env python3
"""Intersect exact mod-p cover-collision sieves involving norm-four traces.

For each prime this checker compares every rational-normalization survivor in
the finite three-node (norm-four) chord chart with

* the complete smooth rational-bisection atlas;
* every survivor in the finite two-node (norm-six) chord chart; and
* survivors attached to a distinct norm-four trace class.

The scalar squareclass is part of every cover key.  An empty intersection is
only a simultaneous-good-reduction obstruction in the displayed affine
charts; it is not a characteristic-zero nonexistence theorem.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PRIMES = (17, 23, 29, 31, 37)
GENERATED = ROOT / "artifacts/generated-results"
DEFAULT_SMOOTH = (
    GENERATED
    / "elkies-k3-r17-norm12-11952-alternate-bisections-full-v1.json"
)
DEFAULT_OUTPUT = (
    GENERATED
    / "elkies-k3-r17-norm12-11952-genus3-cover-collision-intersection-v1.json"
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


def smooth_branch_records(path: Path):
    """Stream the large pretty-printed atlas, retaining only mask and branch."""

    pending = None
    with path.open() as stream:
        iterator = iter(stream)
        for line in iterator:
            if '"numerator_coefficients": [' in line:
                values = []
                for coefficient_line in iterator:
                    value = coefficient_line.strip().rstrip(",")
                    if value == "]":
                        break
                    values.append(int(json.loads(value)))
                if len(values) == 3:
                    pending = values
            elif pending is not None and '"label": ' in line:
                label = json.loads(line.split(":", 1)[1].strip().rstrip(","))
                yield int(label.rsplit("-", 1)[1], 16), pending
                pending = None


def cover_key(coefficients: list[int], prime: int):
    values = [int(value) % prime for value in coefficients]
    if len(values) != 3 or values[2] == 0:
        return None
    if (values[1] * values[1] - 4 * values[0] * values[2]) % prime == 0:
        return None
    inverse = pow(values[2], -1, prime)
    atom = tuple(value * inverse % prime for value in values)
    scalar_character = pow(values[2], (prime - 1) // 2, prime)
    return atom, scalar_character


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--prime",
    type=int,
    action="append",
    help="prime to intersect; defaults to 17,23,29,31,37",
)
parser.add_argument("--smooth-covers", type=Path, default=DEFAULT_SMOOTH)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

primes = tuple(args.prime or DEFAULT_PRIMES)
if len(set(primes)) != len(primes):
    parser.error("primes must be distinct")

smooth_path = args.smooth_covers.resolve()
smooth_records = list(smooth_branch_records(smooth_path))
if len(smooth_records) != 39147:
    raise ValueError(f"expected 39147 smooth covers, found {len(smooth_records)}")

common_smooth = None
common_norm_six = None
common_norm_four = None
prime_records = []
input_paths = [smooth_path]
for prime in primes:
    norm_four_path = (
        GENERATED
        / f"elkies-k3-r17-norm12-11952-genus3-normalization-full-p{prime}-v1.json"
    )
    norm_six_path = (
        GENERATED
        / f"elkies-k3-r17-norm12-11952-genus2-normalization-full-p{prime}-direct-v1.json"
    )
    input_paths.extend((norm_four_path, norm_six_path))
    norm_four = json.loads(norm_four_path.read_text())
    norm_six = json.loads(norm_six_path.read_text())
    if norm_four.get("schema") != (
        "elkies-k3.r17-norm12-direct-genus3-normalization-modp-search.v1"
    ):
        raise ValueError(f"unexpected norm-four schema: {norm_four_path}")
    if norm_six.get("schema") != (
        "elkies-k3.r17-norm12-direct-genus2-normalization-modp-search.v1"
    ):
        raise ValueError(f"unexpected norm-six schema: {norm_six_path}")
    if int(norm_four["prime"]) != prime or int(norm_six["prime"]) != prime:
        raise ValueError(f"prime mismatch at {prime}")

    smooth_by_key = defaultdict(set)
    for mask, coefficients in smooth_records:
        key = cover_key(coefficients, prime)
        if key is not None:
            smooth_by_key[key].add(mask)

    norm_six_by_key = defaultdict(set)
    for survivor in norm_six["survivors"]:
        key = cover_key(
            survivor["reduced_quadratic_coefficients_low_to_high"], prime
        )
        if key is not None:
            norm_six_by_key[key].add(int(survivor["translation_orbit_mask"]))

    norm_four_by_key = defaultdict(set)
    for survivor in norm_four["survivors"]:
        key = cover_key(
            survivor["reduced_quadratic_coefficients_low_to_high"], prime
        )
        if key is not None:
            norm_four_by_key[key].add(int(survivor["translation_orbit_mask"]))

    smooth_pairs = set()
    norm_six_pairs = set()
    norm_four_pairs = set()
    for key, norm_four_masks in norm_four_by_key.items():
        for norm_four_mask in norm_four_masks:
            smooth_pairs.update(
                (norm_four_mask, smooth_mask)
                for smooth_mask in smooth_by_key.get(key, ())
                if smooth_mask != norm_four_mask
            )
            norm_six_pairs.update(
                (norm_four_mask, norm_six_mask)
                for norm_six_mask in norm_six_by_key.get(key, ())
                if norm_six_mask != norm_four_mask
            )
        norm_four_pairs.update(combinations(sorted(norm_four_masks), 2))

    common_smooth = smooth_pairs if common_smooth is None else common_smooth & smooth_pairs
    common_norm_six = (
        norm_six_pairs
        if common_norm_six is None
        else common_norm_six & norm_six_pairs
    )
    common_norm_four = (
        norm_four_pairs
        if common_norm_four is None
        else common_norm_four & norm_four_pairs
    )
    prime_records.append(
        {
            "prime": prime,
            "norm_four_processed_trace_count": int(
                norm_four["search"]["processed_trace_count"]
            ),
            "norm_four_survivor_count": int(norm_four["survivor_count"]),
            "norm_six_processed_trace_count": int(
                norm_six["search"]["processed_trace_count"]
            ),
            "norm_four_vs_smooth_pair_count_at_prime": len(smooth_pairs),
            "norm_four_vs_norm_six_pair_count_at_prime": len(norm_six_pairs),
            "distinct_norm_four_pair_count_at_prime": len(norm_four_pairs),
            "norm_four_vs_smooth_intersection_count_through_prime": len(
                common_smooth
            ),
            "norm_four_vs_norm_six_intersection_count_through_prime": len(
                common_norm_six
            ),
            "distinct_norm_four_intersection_count_through_prime": len(
                common_norm_four
            ),
        }
    )

common_smooth = common_smooth or set()
common_norm_six = common_norm_six or set()
common_norm_four = common_norm_four or set()
output_path = args.output if args.output.is_absolute() else ROOT / args.output
output_path.parent.mkdir(parents=True, exist_ok=True)
result = {
    "schema": "elkies-k3.r17-norm12-direct-genus3-cover-collision-intersection.v1",
    "status": (
        "PASS_SIMULTANEOUS_GOOD_REDUCTION_NO_COVER_COLLISION"
        if not common_smooth and not common_norm_six and not common_norm_four
        else "SURVIVING_MODULAR_COVER_COLLISIONS"
    ),
    "source_label": "norm12-orbit-11952",
    "primes_in_intersection_order": list(primes),
    "smooth_atlas_record_count": len(smooth_records),
    "prime_records": prime_records,
    "surviving_norm_four_vs_smooth_pairs": [
        {"norm_four_trace_mask": left, "smooth_bisection_mask": right}
        for left, right in sorted(common_smooth)
    ],
    "surviving_norm_four_vs_norm_six_pairs": [
        {"norm_four_trace_mask": left, "norm_six_trace_mask": right}
        for left, right in sorted(common_norm_six)
    ],
    "surviving_distinct_norm_four_trace_pairs": [
        {"trace_masks": [left, right]}
        for left, right in sorted(common_norm_four)
    ],
    "proof_boundary": (
        "Every retained pair has the same full finite-field quadratic cover "
        "squareclass, including its scalar character, at every displayed prime. "
        "An empty intersection is only a simultaneous-good-reduction obstruction "
        "for the processed traces in the finite affine slope charts. Trace or "
        "parameter bad reduction and parameter-at-infinity charts remain open."
    ),
    "inputs": {relative(path): digest(path) for path in input_paths},
    "reproducing_command": (
        "python3 "
        "elkies-k3/scripts/analyze_r17_norm12_direct_genus3_cover_collisions.py"
    ),
}
output_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
print(
    "R17GENUS3COVERINTERSECTION"
    f"|primes={','.join(map(str, primes))}"
    f"|smooth_pairs={len(common_smooth)}"
    f"|norm6_pairs={len(common_norm_six)}"
    f"|norm4_pairs={len(common_norm_four)}"
    f"|output={relative(output_path)}|status={result['status']}",
    flush=True,
)
