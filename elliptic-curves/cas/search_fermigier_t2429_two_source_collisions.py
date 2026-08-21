#!/usr/bin/env python3
"""Search auxiliary-slice collisions from the new T=2429/6 rank-14 lead.

The unused-slope genus-two pass through Fermigier E22 produced the low-
conductor fibre ``T=2429/6``.  Exact height triage found two independent
exceptional directions there: the forced point ``x=13973/6`` and the new
height-50000 point ``x=17731/1654``.  Through each source point this script
constructs the two genus-one slices ``x=+/-T+b``, maps their generic-section
intersections to the pointed auxiliary elliptic curve, and enumerates every
signed support-at-most-two seed combination.

The target of this bounded pass is a parameter appearing from both distinct
source directions with distinct quartic abscissas.  Exact membership is
checked throughout.  No completeness or rank claim follows from no collision.
"""

from __future__ import annotations

import argparse
from dataclasses import replace
from fractions import Fraction
import hashlib
from itertools import combinations, product
import json
from pathlib import Path

from ek_k3 import rational_to_string
from fermigier_mestre import FermigierMestreFamily
from search_fermigier_rank22_accidental_slices import Slice, slice_polynomial
import search_fermigier_rank22_auxiliary_orbits as auxiliary_tools
from search_nagao_section7_auxiliary_jacobians import (
    weierstrass_add,
    weierstrass_multiply,
)


Q = Fraction
SOURCE_PARAMETER = Q(2429, 6)
SOURCES = (
    ("forced-P13", (Q(13973, 6), Q(35843568404, 9))),
    (
        "H50000-independent",
        (Q(17731, 1654), Q(2059430795581229, 18466083)),
    ),
)


def projective_height(value: Fraction) -> int:
    value = Q(value)
    return max(abs(value.numerator), value.denominator)


def low_weight_points(auxiliary, seeds, *, maximum_support: int):
    coefficients = auxiliary.weierstrass_coefficients
    seen = set()

    def emit(point, label):
        if point is None:
            return None
        key = Q(point[0]), Q(point[1])
        if key in seen:
            return None
        seen.add(key)
        return key, label

    for index, point in enumerate(seeds, 1):
        for scalar in (1, -1, 2, -2):
            emitted = emit(
                weierstrass_multiply(coefficients, point, scalar),
                f"{scalar:+d}P{index}",
            )
            if emitted is not None:
                yield emitted
    for left in range(len(seeds)):
        for right in range(left+1, len(seeds)):
            for sign in (1, -1):
                candidate = weierstrass_add(
                    coefficients,
                    seeds[left],
                    weierstrass_multiply(coefficients, seeds[right], sign),
                )
                emitted = emit(
                    candidate,
                    f"P{left+1}{'+' if sign == 1 else '-'}P{right+1}",
                )
                if emitted is not None:
                    yield emitted
    for support in range(3, maximum_support+1):
        for indices in combinations(range(len(seeds)), support):
            for signs in product((1, -1), repeat=support):
                candidate = None
                for index, sign in zip(indices, signs, strict=True):
                    candidate = weierstrass_add(
                        coefficients,
                        candidate,
                        weierstrass_multiply(coefficients, seeds[index], sign),
                    )
                label = "+".join(
                    f"{'-' if sign < 0 else ''}P{index+1}"
                    for index, sign in zip(indices, signs, strict=True)
                )
                emitted = emit(candidate, label)
                if emitted is not None:
                    yield emitted


def digest(parameters) -> str:
    value = hashlib.sha256()
    for parameter in sorted(parameters, key=lambda t: (projective_height(t), t)):
        value.update((rational_to_string(parameter)+"\n").encode())
    return value.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "artifacts/local/elliptic-curves/"
            "fermigier-t2429-two-source-collisions-v1.json"
        ),
    )
    parser.add_argument("--maximum-support", type=int, default=3)
    arguments = parser.parse_args()
    repository = Path(__file__).resolve().parents[2]
    output = arguments.output if arguments.output.is_absolute() else repository/arguments.output

    # The reusable pointed-quartic helpers take their source parameter from
    # this module constant.  This process owns no other slice calculation, so
    # pin it to the new anchor before constructing any auxiliary object.
    auxiliary_tools.T0 = SOURCE_PARAMETER

    candidates = {}
    slice_records = []
    for source_label, source_point in SOURCES:
        if source_point[1]**2 != FermigierMestreFamily.quartic_value(
            SOURCE_PARAMETER, source_point[0]
        ):
            raise AssertionError("a pinned T=2429/6 source missed the quartic")
        for slope in (-1, 1):
            intercept = source_point[0]-Q(slope)*SOURCE_PARAMETER
            data = Slice(
                accidental_label=source_label,
                source_point=source_point,
                slope=slope,
                intercept=intercept,
                coefficients=slice_polynomial(slope, intercept),
            )
            auxiliary = auxiliary_tools.make_auxiliary_slice(data)
            slice_id = f"{source_label}_{'p1' if slope == 1 else 'm1'}"
            auxiliary = replace(auxiliary, slice_id=slice_id)
            seeds, intersections = auxiliary_tools.exact_seed_points(auxiliary)
            generic_parameters = {
                abs(parameter)
                for parameter, _ in auxiliary_tools.generic_intersections(auxiliary)
            }
            accepted = 0
            for point, combination in low_weight_points(
                auxiliary, seeds, maximum_support=arguments.maximum_support
            ):
                inverse = auxiliary.inverse(point)
                if inverse is None:
                    continue
                signed_parameter, ordinate = inverse
                parameter = abs(signed_parameter)
                if parameter in (Q(0), SOURCE_PARAMETER) or parameter in generic_parameters:
                    continue
                if FermigierMestreFamily.discriminant_factor(parameter) == 0:
                    continue
                quartic_x = auxiliary.original_x(signed_parameter)
                if ordinate**2 != FermigierMestreFamily.quartic_value(
                    parameter, quartic_x
                ):
                    raise AssertionError("an auxiliary inverse missed the Fermigier fibre")
                candidates.setdefault(parameter, []).append(
                    {
                        "source_label": source_label,
                        "slice_id": slice_id,
                        "combination": combination,
                        "signed_parameter": signed_parameter,
                        "quartic_x": quartic_x,
                        "quartic_y": ordinate,
                    }
                )
                accepted += 1
            slice_records.append(
                {
                    "slice_id": slice_id,
                    "source_label": source_label,
                    "slope": slope,
                    "seed_count": len(seeds),
                    "generic_intersection_count": len(intersections),
                    "accepted_incidence_count": accepted,
                }
            )
            print(
                f"T2429AUX|slice={slice_id}|seeds={len(seeds)}|"
                f"accepted={accepted}|unique_global={len(candidates)}",
                flush=True,
            )

    collisions = []
    for parameter, rows in candidates.items():
        sources = {row["source_label"] for row in rows}
        abscissas = {row["quartic_x"] for row in rows}
        if len(sources) < 2 or len(abscissas) < 2:
            continue
        representatives = []
        seen_sources = set()
        for row in rows:
            if row["source_label"] in seen_sources:
                continue
            seen_sources.add(row["source_label"])
            representatives.append(
                {
                    key: rational_to_string(value) if isinstance(value, Q) else value
                    for key, value in row.items()
                }
            )
        collisions.append(
            {
                "parameter": rational_to_string(parameter),
                "projective_height": projective_height(parameter),
                "distinct_source_count": len(sources),
                "distinct_abscissa_count": len(abscissas),
                "representatives": representatives,
            }
        )
    collisions.sort(key=lambda row: (row["projective_height"], Q(row["parameter"])))
    payload = {
        "schema_version": 1,
        "artifact_kind": "bounded_fermigier_two_source_auxiliary_collision_search",
        "status": "collision_found" if collisions else "no_collision_in_declared_ball",
        "source_parameter": rational_to_string(SOURCE_PARAMETER),
        "source_points": [
            {
                "label": label,
                "x": rational_to_string(point[0]),
                "y": rational_to_string(point[1]),
            }
            for label, point in SOURCES
        ],
        "combination_ball": (
            "singletons, doubles, and signed seed sums through support "
            f"{arguments.maximum_support}"
        ),
        "slice_records": slice_records,
        "unique_parameter_count": len(candidates),
        "unique_parameter_sha256": digest(candidates),
        "collision_count": len(collisions),
        "collisions": collisions,
        "claim_boundary": "exact bounded search; no Mordell-Weil completeness or rank claim",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
    print(
        f"T2429AUX|slices=4|unique={len(candidates)}|"
        f"collisions={len(collisions)}|output={output}|status=PASS_BOUNDED",
        flush=True,
    )


if __name__ == "__main__":
    main()
