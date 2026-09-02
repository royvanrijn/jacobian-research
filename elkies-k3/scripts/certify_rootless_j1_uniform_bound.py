#!/usr/bin/env python3
"""Certify the finite quadratic-form part of the rootless J1 bound.

status: ACTIVE_PROOF
claim: For the pinned Picard-rank-19 H3 surface, the discriminant-form
  isometry group has order eight and the image of the Hodge isometry group is
  {+1,-1}.  Braun--Kimura--Watari Proposition C' therefore bounds the J1
  multiplicity of either rootless J2 class by four.  Since the rootless J2
  classification has exactly two classes, the rootless J1 count is in [2,8].
inputs: artifacts/generated-results/elkies-k3-rootless-j2-niemeier-first.json,
  artifacts/generated-results/elkies-k3-rootless-j2-niemeier-controls.json
outputs: artifacts/generated-results/elkies-k3-rootless-j1-uniform-bound-v1.json
supersedes/superseded-by: none

The rank-three Hodge-rigidity argument and the application of Proposition C'
are mathematical parts of Theorem H2a in RANK_MUTATION_AND_LIFT_THEOREMS.md.
This dependency-free replay checks the remaining exact finite arithmetic and
locks the two complete J2 inputs.  It deliberately does not turn the smaller
frame-automorphism double quotients into J1 bounds.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CLASSIFICATION = (
    ROOT
    / "artifacts/generated-results/elkies-k3-rootless-j2-niemeier-first.json"
)
CONTROLS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-rootless-j2-niemeier-controls.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-rootless-j1-uniform-bound-v1.json"
)
MODULUS = 948


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def form_isometry_units(q: Fraction) -> list[int]:
    """Return units u with q*u^2=q in Q/2Z for a cyclic generator."""

    result = []
    for unit in range(MODULUS):
        if math.gcd(unit, MODULUS) != 1:
            continue
        difference = q * (unit * unit - 1)
        if difference.denominator == 1 and difference.numerator % 2 == 0:
            result.append(unit)
    return result


def left_cosets(group: list[int], subgroup: list[int]) -> list[list[int]]:
    unseen = set(group)
    cosets = []
    while unseen:
        representative = min(unseen)
        coset = sorted(
            {(element * representative) % MODULUS for element in subgroup}
        )
        assert set(coset) <= set(group)
        unseen.difference_update(coset)
        cosets.append(coset)
    return cosets


def build_certificate() -> dict:
    classification = load_json(CLASSIFICATION)
    controls = load_json(CONTROLS)

    assert classification["schema"] == "elkies-k3.rootless-j2-niemeier-first.v1"
    assert classification["status"] == "PASS_COMPLETE_ROOTLESS_J2_CLASSIFICATION"
    assert len(classification["rootless_classes"]) == 2
    assert all(row["determinant"] == MODULUS for row in classification["rootless_classes"])

    assert controls["schema"] == "elkies-k3.rootless-j2-niemeier-controls.v1"
    assert len(controls["controls"]) == 2
    assert [row["label"] for row in controls["controls"]] == [
        "published_R17",
        "alternate_Q80",
    ]

    hodge_image = [1, MODULUS - 1]
    rows = []
    common_group = None
    for control in controls["controls"]:
        q = Fraction(control["frame_discriminant_generator_q"])
        group = form_isometry_units(q)
        assert hodge_image == [unit for unit in group if unit in hodge_image]
        cosets = left_cosets(group, hodge_image)
        assert len(group) == 8
        assert len(cosets) == 4
        if common_group is None:
            common_group = group
        assert group == common_group
        rows.append(
            {
                "j2_class": control["label"],
                "cyclic_discriminant_form_generator_q": str(q),
                "orthogonal_group_units_mod_948": group,
                "hodge_image_units_mod_948": hodge_image,
                "hodge_left_cosets": cosets,
                "uniform_j1_multiplicity_upper_bound": len(cosets),
                "frame_automorphism_discriminant_image_recorded_but_not_used": control[
                    "frame_automorphism_discriminant_image"
                ],
            }
        )

    rootless_j2_count = len(classification["rootless_classes"])
    per_class_bound = max(row["uniform_j1_multiplicity_upper_bound"] for row in rows)
    return {
        "schema": "elkies-k3.rootless-j1-uniform-bound.v1",
        "status": "PASS_EXACT_ROOTLESS_J1_UNIFORM_BOUND_NOT_CLASSIFICATION",
        "inputs": {
            str(CLASSIFICATION.relative_to(ROOT)): sha256(CLASSIFICATION),
            str(CONTROLS.relative_to(ROOT)): sha256(CONTROLS),
        },
        "theorem_inputs": {
            "picard_rank": 19,
            "transcendental_rank": 3,
            "hodge_isometry_group": ["+identity", "-identity"],
            "hodge_discriminant_image_units_mod_948": hodge_image,
            "literature_bound": (
                "Braun--Kimura--Watari, arXiv:1312.4421, Proposition C'"
            ),
        },
        "finite_form_rows": rows,
        "conclusion": {
            "complete_rootless_j2_class_count": rootless_j2_count,
            "uniform_j1_multiplicity_upper_bound_per_j2_class": per_class_bound,
            "rootless_j1_class_count_lower_bound": rootless_j2_count,
            "rootless_j1_class_count_upper_bound": rootless_j2_count * per_class_bound,
        },
        "scope": {
            "proved": (
                "The exact finite-form quotient in the uniform J1 multiplicity "
                "bound has four cosets for each of the two complete rootless J2 "
                "classes. Together with the rank-three Hodge-rigidity lemma, the "
                "rootless J1 count is at least two and at most eight."
            ),
            "not_proved": (
                "The exact J1 count, representatives for its surface-automorphism "
                "orbits, or a frame-dependent improvement from the recorded frame "
                "automorphism images."
            ),
        },
    }


def canonical_bytes(value: dict) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    payload = canonical_bytes(build_certificate())
    if args.check:
        assert args.output.read_bytes() == payload
        print(f"PASS byte-identical {args.output.relative_to(ROOT)}")
        return

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(payload)
    print(f"wrote {args.output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
