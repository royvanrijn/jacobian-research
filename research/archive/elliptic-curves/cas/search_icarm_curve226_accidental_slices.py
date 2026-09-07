#!/usr/bin/env python3
"""Bounded accidental-slice deformation audit around ICARM curve 226.

The submitted family label decodes to the exact six-root tuple below.  This
script searches the anchor quartic and all ``x=+/-T+n`` slices through its
bounded non-generic abscissas.  It carries the slice ordinates to every
cross-source collision before computing a numerical height rank; a
candidate-centered quartic search can otherwise miss these transported
points.

This is a bounded experiment.  Numerical height ranks are not Mordell--Weil
certificates, and no absence statement extends beyond the declared boxes.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
import json
from pathlib import Path

from mestre_root_tuples import SixRootMestreConstruction
from search_fermigier_rank22_accidental_slices import (
    canonical_signless_points,
    search_polynomial,
)
from search_icarm_curve243_accidental_slices import conductor_record, qtext
from search_mestre_root_tuple_scale import (
    primitive_visible_points,
    quartic_point_to_jacobian,
)
from search_six_root_low_conductor_centers import (
    bivariate_quartic_coefficients,
    slice_polynomial,
)
from triage_nagao_rank13_finalists import height_matrix_replay, stable_height_rank


Q = Fraction
ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = (
    ROOT / "artifacts/local/elliptic-curves"
    / "icarm226-accidental-slices-h200000-v1.json"
)
ANCHOR = Q(10167, 350)
ROOTS = tuple(map(Q, (-138, -90, -60, -12, 138, 162)))
EXPECTED_SINGULAR = (Q(-69), Q(69))
EXPECTED_REGULAR = (
    Q(10167, 2800),
    Q(10167, 1225),
    Q(10167, 550),
    Q(111837, 2450),
    Q(10167, 100),
    Q(40668, 175),
)


def build_payload(*, height: int, timeout: float, stack_bytes: int) -> dict:
    construction = SixRootMestreConstruction(ROOTS)
    if not construction.is_quartic_family:
        raise AssertionError("the curve-226 roots left the Mestre quartic locus")
    anchor_points, anchor_search = search_polynomial(
        construction.quartic_coefficients(ANCHOR),
        height_bound=height,
        timeout=timeout,
        stack_bytes=stack_bytes,
    )
    if anchor_search["status"] != "completed":
        raise RuntimeError("anchor search did not complete")
    anchor_x = sorted(
        {point[0] for point in canonical_signless_points(anchor_points)}
    )
    visible_x = {
        root + sign * ANCHOR
        for root in construction.roots
        for sign in (-1, 1)
    }
    if len(set(anchor_x).intersection(visible_x)) != 12:
        raise AssertionError("the visible anchor abscissas changed")
    accidentals = tuple(value for value in anchor_x if value not in visible_x)
    if len(accidentals) != 19:
        raise AssertionError("the bounded accidental set changed")

    bivariate = bivariate_quartic_coefficients(construction)
    sources: dict[Fraction, set[tuple[int, int]]] = defaultdict(set)
    transported: dict[Fraction, list[tuple[Fraction, Fraction]]] = defaultdict(list)
    slices = []
    for source_index, source_x in enumerate(accidentals):
        for slope in (-1, 1):
            intercept = source_x-slope*ANCHOR
            polynomial = slice_polynomial(bivariate, Q(slope), intercept)
            if len(polynomial) != 5:
                raise AssertionError("a +/-1 slice did not become quartic")
            points, search = search_polynomial(
                polynomial,
                height_bound=height,
                timeout=timeout,
                stack_bytes=stack_bytes,
            )
            if search["status"] != "completed":
                raise RuntimeError(f"slice {source_index}/{slope} did not complete")
            signless = canonical_signless_points(points)
            for parameter, ordinate in signless:
                sources[parameter].add((source_index, slope))
                transported[parameter].append(
                    (Q(slope)*parameter+intercept, ordinate)
                )
            slices.append(
                {
                    "source_index": source_index,
                    "slope": slope,
                    "intercept": qtext(intercept),
                    "parameter_count": len(signless),
                }
            )

    collision_parameters = tuple(
        parameter
        for parameter, signed_sources in sorted(sources.items())
        if parameter != ANCHOR
        and len({index for index, _ in signed_sources}) >= 2
    )
    if collision_parameters != tuple(sorted(EXPECTED_SINGULAR+EXPECTED_REGULAR)):
        raise AssertionError("the bounded collision set changed")

    collision_rows = []
    for parameter in collision_parameters:
        singular = construction.primitive_discriminant_value(parameter) == 0
        row = {
            "parameter": qtext(parameter),
            "distinct_accidental_sources": len(
                {index for index, _ in sources[parameter]}
            ),
            "signed_sources": [list(value) for value in sorted(sources[parameter])],
            "singular_quartic": singular,
        }
        if not singular:
            quartic_points = list(primitive_visible_points(construction, parameter))
            quartic_points.extend(transported[parameter])
            images = []
            seen_x = set()
            for point in quartic_points:
                image = quartic_point_to_jacobian(construction, parameter, point)
                if image[0] in seen_x:
                    continue
                seen_x.add(image[0])
                images.append(image)
            height_rows = height_matrix_replay(
                construction.primitive_jacobian_coefficients(parameter),
                images,
                precisions=(96, 192),
                timeout=timeout,
                stack_bytes=stack_bytes,
            )
            row.update(
                {
                    "transported_signed_points": len(transported[parameter]),
                    "distinct_jacobian_x": len(images),
                    "stable_numerical_rank": stable_height_rank(height_rows),
                    "height_replay": list(height_rows),
                    "conductor": conductor_record(
                        construction,
                        parameter,
                        timeout=timeout,
                        stack_bytes=stack_bytes,
                    ),
                }
            )
            if row["stable_numerical_rank"] != 11:
                raise AssertionError("a collision numerical rank changed")
        collision_rows.append(row)

    return {
        "schema_version": 1,
        "status": "bounded_exact_accidental_slice_experiment",
        "icarm_curve": 226,
        "anchor_parameter": qtext(ANCHOR),
        "roots": [qtext(value) for value in construction.roots],
        "height_bound": height,
        "anchor_signed_points": anchor_search["signed_point_count"],
        "anchor_distinct_abscissas": len(anchor_x),
        "visible_abscissas": 12,
        "bounded_accidental_abscissas": [qtext(value) for value in accidentals],
        "slice_count": len(slices),
        "slices": slices,
        "cross_source_collisions": collision_rows,
        "conclusion": (
            "the two all-source collisions are singular; every regular "
            "collision has stable numerical rank 11 after carrying the exact "
            "slice points, so this bounded lane is stopped"
        ),
        "rank_claim": None,
        "reproduction": (
            "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
            "elliptic-curves/cas/search_icarm_curve226_accidental_slices.py"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--height", type=int, default=200_000)
    parser.add_argument("--timeout", type=float, default=30.0)
    parser.add_argument("--stack-bytes", type=int, default=1_000_000_000)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    payload = build_payload(
        height=arguments.height,
        timeout=arguments.timeout,
        stack_bytes=arguments.stack_bytes,
    )
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
    print(
        "ICARM226SLICE|"
        f"accidentals={len(payload['bounded_accidental_abscissas'])}|"
        f"slices={payload['slice_count']}|"
        f"collisions={len(payload['cross_source_collisions'])}|"
        f"output={arguments.output}|status=PASS",
        flush=True,
    )


if __name__ == "__main__":
    main()
