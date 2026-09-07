#!/usr/bin/env python3
"""Exclude three-cover collisions in the processed R17 low-genus charts.

The existing collision ledgers compare distinct trace masks.  This checker
adds the missing same-mask audit: at each prime it groups the complete smooth
atlas and every affine genus-two/genus-three normalization survivor by the
full scalar-sensitive quadratic cover key.  A characteristic-zero triple in
one trace mask would force that mask to occur with multiplicity at least three
at every good prime.

Combining an empty same-mask triple intersection with the certified distinct-
mask pair ledgers excludes a three-cover collision in the displayed layers.
Boundary charts, bad-reduction denominators, and arithmetic genus at least
four remain open.
"""

# status: UNPROMOTED_RESULT
# claim: no three-cover collision in the processed affine low-genus 11952 layers
# inputs: pinned smooth, genus-two, and genus-three normalization ledgers
# outputs: artifacts/generated-results/elkies-k3-r17-norm12-11952-low-genus-rank3-cover-collision-v1.json

from __future__ import annotations

import argparse
from collections import defaultdict
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
DEFAULT_PRIMES = (17, 23, 29, 31, 37)
SMOOTH = GENERATED / "elkies-k3-r17-norm12-11952-alternate-bisections-full-v1.json"
GENUS2_COLLISIONS = (
    GENERATED
    / "elkies-k3-r17-norm12-11952-genus2-cover-collision-intersection-v1.json"
)
GENUS3_COLLISIONS = (
    GENERATED
    / "elkies-k3-r17-norm12-11952-genus3-cover-collision-intersection-v1.json"
)
GENUS3_TARGETED = (
    GENERATED
    / "elkies-k3-r17-norm12-11952-genus3-cover-collision-targeted-v1.json"
)
DEFAULT_OUTPUT = (
    GENERATED
    / "elkies-k3-r17-norm12-11952-low-genus-rank3-cover-collision-v1.json"
)


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def smooth_records(path: Path):
    """Stream the large pretty-printed atlas, retaining only mask and branch."""

    pending = None
    with path.open() as source:
        iterator = iter(source)
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


def cover_key(coefficients, prime):
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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    smooth = list(smooth_records(SMOOTH))
    if len(smooth) != 39147:
        raise ValueError(f"expected 39147 smooth covers, found {len(smooth)}")
    g2_collision = json.loads(GENUS2_COLLISIONS.read_text())
    g3_collision = json.loads(GENUS3_COLLISIONS.read_text())
    g3_targeted = json.loads(GENUS3_TARGETED.read_text())
    if g2_collision.get("status") != "PASS_SIMULTANEOUS_GOOD_REDUCTION_NO_COVER_COLLISION":
        raise ValueError("the distinct-mask genus-two collision gate is not closed")
    if g3_targeted.get("status") != "PASS_TARGETED_GOOD_REDUCTION_NO_COVER_COLLISION":
        raise ValueError("the distinct-mask genus-three targeted gate is not closed")
    if g3_collision.get("status") != "SURVIVING_MODULAR_COVER_COLLISIONS":
        raise ValueError("the genus-three input survivor ledger changed")
    if g3_targeted.get("input_candidate_pairs") != [
        record["trace_masks"]
        for record in g3_collision["surviving_distinct_norm_four_trace_pairs"]
    ]:
        raise ValueError("the targeted genus-three gate does not cover every survivor")

    common_pair_masks = None
    common_triple_masks = None
    prime_records = []
    input_paths = [
        SMOOTH,
        GENUS2_COLLISIONS,
        GENUS3_COLLISIONS,
        GENUS3_TARGETED,
    ]
    for prime in DEFAULT_PRIMES:
        g2_path = (
            GENERATED
            / f"elkies-k3-r17-norm12-11952-genus2-normalization-full-p{prime}-direct-v1.json"
        )
        g3_path = (
            GENERATED
            / f"elkies-k3-r17-norm12-11952-genus3-normalization-full-p{prime}-v1.json"
        )
        input_paths.extend((g2_path, g3_path))
        g2 = json.loads(g2_path.read_text())
        g3 = json.loads(g3_path.read_text())
        if int(g2["prime"]) != prime or int(g3["prime"]) != prime:
            raise ValueError(f"prime mismatch at p={prime}")

        groups = defaultdict(set)
        for mask, coefficients in smooth:
            key = cover_key(coefficients, prime)
            if key is not None:
                groups[(mask, key)].add(("smooth", mask))
        for survivor in g2["survivors"]:
            key = cover_key(
                survivor["reduced_quadratic_coefficients_low_to_high"], prime
            )
            if key is not None:
                groups[(int(survivor["translation_orbit_mask"]), key)].add(
                    (
                        "genus2",
                        int(survivor["trace_index"]),
                        tuple(map(int, survivor["l0_l1"])),
                    )
                )
        for survivor in g3["survivors"]:
            key = cover_key(
                survivor["reduced_quadratic_coefficients_low_to_high"], prime
            )
            if key is not None:
                groups[(int(survivor["translation_orbit_mask"]), key)].add(
                    (
                        "genus3",
                        int(survivor["trace_index"]),
                        tuple(map(int, survivor["m0_m1_m2"])),
                    )
                )

        pair_masks = {
            mask for (mask, _), objects in groups.items() if len(objects) >= 2
        }
        triple_masks = {
            mask for (mask, _), objects in groups.items() if len(objects) >= 3
        }
        pair_group_count = sum(len(objects) >= 2 for objects in groups.values())
        triple_group_count = sum(len(objects) >= 3 for objects in groups.values())
        common_pair_masks = (
            pair_masks
            if common_pair_masks is None
            else common_pair_masks & pair_masks
        )
        common_triple_masks = (
            triple_masks
            if common_triple_masks is None
            else common_triple_masks & triple_masks
        )
        prime_records.append(
            {
                "prime": prime,
                "genus2_survivor_count": int(g2["survivor_count"]),
                "genus3_survivor_count": int(g3["survivor_count"]),
                "same_mask_cover_groups_of_size_at_least_2": pair_group_count,
                "same_mask_cover_groups_of_size_at_least_3": triple_group_count,
                "same_mask_pair_trace_mask_count": len(pair_masks),
                "same_mask_triple_trace_mask_count": len(triple_masks),
                "pair_trace_mask_intersection_count_through_prime": len(
                    common_pair_masks
                ),
                "triple_trace_mask_intersection_count_through_prime": len(
                    common_triple_masks
                ),
            }
        )

    common_pair_masks = common_pair_masks or set()
    common_triple_masks = common_triple_masks or set()
    if common_triple_masks:
        status = "SURVIVING_SAME_MASK_LOW_GENUS_TRIPLE_CANDIDATES"
    else:
        status = "PASS_PROCESSED_AFFINE_LOW_GENUS_NO_THREE_COVER_COLLISION"
    payload = {
        "schema": "elkies-k3.r17-norm12-low-genus-rank3-cover-collision.v1",
        "status": status,
        "source_label": "norm12-orbit-11952",
        "layers": [
            "complete smooth rational bisection atlas",
            "finite affine arithmetic-genus-two normalization chart",
            "finite affine arithmetic-genus-three normalization chart",
        ],
        "primes_in_intersection_order": list(DEFAULT_PRIMES),
        "prime_records": prime_records,
        "surviving_same_mask_pair_trace_masks": sorted(common_pair_masks),
        "surviving_same_mask_triple_trace_masks": sorted(common_triple_masks),
        "distinct_mask_pair_gate": {
            "genus2_status": g2_collision["status"],
            "genus3_initial_status": g3_collision["status"],
            "genus3_targeted_status": g3_targeted["status"],
            "genus3_targeted_surviving_pairs": g3_targeted["surviving_pairs"],
        },
        "mathematical_consequence": (
            "No scalar-sensitive quadratic cover occurs three times among the "
            "processed smooth, genus-two, and genus-three affine layers with "
            "simultaneous good integral reduction in the displayed charts. Thus "
            "these layers cannot supply three independent directions on one twist."
        ),
        "proof_boundary": (
            "An empty same-mask trace intersection is a necessary-condition "
            "obstruction, combined with the separately certified distinct-mask pair "
            "obstructions. It is not a global twist-rank upper bound. Rational "
            "solutions nonintegral at a selected prime, parameter-at-infinity or "
            "other boundary charts, and bisections of arithmetic genus at least four "
            "remain open."
        ),
        "inputs": {relative(path): digest(path) for path in input_paths},
        "reproducing_command": (
            "python3 elkies-k3/scripts/analyze_r17_norm12_low_genus_rank3_collisions.py"
        ),
    }
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = args.output if args.output.is_absolute() else ROOT / args.output
    if args.check:
        if output.read_text() != encoded:
            raise SystemExit("stored artifact differs from replay")
        print(f"PASS check {relative(output)}")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(encoded)
    print(
        "R17LOWGENUSRANK3"
        f"|same_mask_pairs={len(common_pair_masks)}"
        f"|same_mask_triples={len(common_triple_masks)}"
        f"|output={relative(output)}|status={status}"
    )


if __name__ == "__main__":
    main()
