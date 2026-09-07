#!/usr/bin/env python3
"""Bounded accidental-slice deformation audit around ICARM curve 243.

Curve 243 was submitted as the ``T=3895/6`` specialization of the exact
six-root Mestre family recorded below.  This experiment searches the anchor
quartic and then every ``x=+/-T+n`` quartic slice through the non-generic
anchor abscissas.  Cross-source parameter collisions force two distinct
accidental quartic directions on a new fibre.

All searches are bounded experiments.  A collision is not a Mordell--Weil
independence certificate, and the conductor computation supplies no rank
claim.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
import json
from pathlib import Path
import subprocess

from mestre_root_tuples import SixRootMestreConstruction
from search_fermigier_rank22_accidental_slices import (
    canonical_signless_points,
    search_polynomial,
)
from search_six_root_low_conductor_centers import (
    bivariate_quartic_coefficients,
    slice_polynomial,
)


Q = Fraction
ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = (
    ROOT / "artifacts/local/elliptic-curves"
    / "icarm243-accidental-slices-h200000-v1.json"
)
ANCHOR = Q(3895, 6)
ROOTS = (
    Q(-1455, 4),
    Q(2955, 4),
    Q(1437, 2),
    Q(-1149, 4),
    Q(-1851, 4),
    Q(-687, 2),
)
EXPECTED_COLLISIONS = (Q(27265, 144), Q(15580, 7))


def qtext(value: Fraction) -> str:
    value = Q(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def conductor_record(
    construction: SixRootMestreConstruction,
    parameter: Fraction,
    *,
    timeout: float,
    stack_bytes: int,
) -> dict:
    coefficients = construction.primitive_jacobian_coefficients(parameter)
    vector = ",".join(
        f"({value.numerator}/{value.denominator})" for value in coefficients
    )
    program = (
        "default(realprecision,60);"
        f"E=ellminimalmodel(ellinit([{vector}]));"
        'N=ellglobalred(E)[1];print("N|",N);'
        'print("LOGN|",log(N));print("ROOT|",ellrootno(E));quit\n'
    )
    process = subprocess.run(
        ["gp", "-q", "-s", str(stack_bytes)],
        input=program,
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    if process.returncode != 0 or "***" in process.stderr:
        raise RuntimeError(f"PARI/GP failed: {process.stderr.strip()}")
    rows = {}
    for line in process.stdout.splitlines():
        if "|" in line:
            key, value = line.strip().split("|", 1)
            rows[key] = value
    return {
        "conductor": rows["N"],
        "log_conductor": rows["LOGN"],
        "root_number": int(rows["ROOT"]),
        "below_182_72": float(rows["LOGN"]) < 182.72,
    }


def build_payload(
    *,
    height: int,
    timeout: float,
    stack_bytes: int,
) -> dict:
    construction = SixRootMestreConstruction(ROOTS)
    if not construction.is_quartic_family:
        raise AssertionError("the curve-243 roots left the Mestre quartic locus")

    anchor_points, anchor_search = search_polynomial(
        construction.quartic_coefficients(ANCHOR),
        height_bound=height,
        timeout=timeout,
        stack_bytes=stack_bytes,
    )
    if anchor_search["status"] != "completed":
        raise RuntimeError("anchor quartic search did not complete")
    anchor_x = sorted(
        {point[0] for point in canonical_signless_points(anchor_points)}
    )
    generic_x = {
        root + sign * ANCHOR
        for root in construction.roots
        for sign in (-1, 1)
    }
    if len(set(anchor_x).intersection(generic_x)) != 12:
        raise AssertionError("the twelve visible anchor abscissas changed")
    accidentals = tuple(value for value in anchor_x if value not in generic_x)
    if len(accidentals) != 16:
        raise AssertionError("the bounded curve-243 accidental set changed")

    bivariate = bivariate_quartic_coefficients(construction)
    parameter_sources: dict[Fraction, set[tuple[int, int]]] = defaultdict(set)
    slices = []
    for source_index, source_x in enumerate(accidentals):
        for slope in (-1, 1):
            intercept = source_x - slope * ANCHOR
            polynomial = slice_polynomial(bivariate, Q(slope), intercept)
            if len(polynomial) != 5:
                raise AssertionError("a +/-1 slice did not reduce to a quartic")
            points, search = search_polynomial(
                polynomial,
                height_bound=height,
                timeout=timeout,
                stack_bytes=stack_bytes,
            )
            if search["status"] != "completed":
                raise RuntimeError(
                    f"slice {source_index}/{slope} did not complete"
                )
            parameters = sorted(
                {point[0] for point in canonical_signless_points(points)}
            )
            for parameter in parameters:
                parameter_sources[parameter].add((source_index, slope))
            slices.append(
                {
                    "source_index": source_index,
                    "slope": slope,
                    "intercept": qtext(intercept),
                    "parameter_count": len(parameters),
                }
            )

    collisions = []
    for parameter, signed_sources in sorted(parameter_sources.items()):
        source_indices = sorted({source for source, _ in signed_sources})
        if parameter == ANCHOR or len(source_indices) < 2:
            continue
        collisions.append(
            {
                "parameter": qtext(parameter),
                "signed_sources": [list(row) for row in sorted(signed_sources)],
                "distinct_accidental_sources": len(source_indices),
                "conductor": conductor_record(
                    construction,
                    parameter,
                    timeout=timeout,
                    stack_bytes=stack_bytes,
                ),
            }
        )
    if tuple(Q(row["parameter"]) for row in collisions) != EXPECTED_COLLISIONS:
        raise AssertionError("the bounded cross-source collision set changed")

    return {
        "schema_version": 1,
        "status": "bounded_exact_accidental_slice_experiment",
        "icarm_curve": 243,
        "anchor_parameter": qtext(ANCHOR),
        "roots": [qtext(value) for value in construction.roots],
        "height_bound": height,
        "anchor_signed_points": anchor_search["signed_point_count"],
        "anchor_distinct_abscissas": len(anchor_x),
        "visible_abscissas": 12,
        "bounded_accidental_abscissas": [qtext(value) for value in accidentals],
        "slice_count": len(slices),
        "slices": slices,
        "cross_source_collisions": collisions,
        "conclusion": (
            "two non-anchor parameters force two distinct accidental quartic "
            "directions in the declared boxes; both conductors exceed the target"
        ),
        "rank_claim": None,
        "reproduction": (
            "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
            "elliptic-curves/cas/search_icarm_curve243_accidental_slices.py"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--height", type=int, default=200_000)
    parser.add_argument("--timeout", type=float, default=20.0)
    parser.add_argument("--stack-bytes", type=int, default=1_000_000_000)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    payload = build_payload(
        height=arguments.height,
        timeout=arguments.timeout,
        stack_bytes=arguments.stack_bytes,
    )
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        "ICARM243SLICE|"
        f"accidentals={len(payload['bounded_accidental_abscissas'])}|"
        f"slices={payload['slice_count']}|"
        f"collisions={len(payload['cross_source_collisions'])}|"
        f"output={arguments.output}|status=PASS",
        flush=True,
    )


if __name__ == "__main__":
    main()
