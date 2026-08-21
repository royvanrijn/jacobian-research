#!/usr/bin/env python3
"""Search unused affine slopes through Fermigier E22 accidental points.

Previous auxiliary-Jacobian searches used only lines ``x=+/-T+b``.  For any
other rational slope, a line through an accidental point on the E22 fibre
still meets each of the thirteen generic section lines at a rational
parameter.  Those intersections give an exact seed subgroup on the pointed
genus-one slice.

The cancellation to genus one occurs only for slopes ``+/-1``.  The smallest
unused integral slopes instead give degree-six, genus-two curves.  This
bounded experiment searches those curves directly.  Its primary gate is a
collision of the same Fermigier parameter across distinct accidental source
points: such a collision preserves more than one exceptional section and is
more relevant to rank 21 than merely widening a one-section height box.

All curve membership and inverse maps are checked exactly.  The calculation
is not a complete Mordell--Weil enumeration and a collision is not by itself
an independence certificate.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any

from ek_k3 import rational_to_string
from fermigier_mestre import FermigierMestreFamily
from search_fermigier_published_pair_fiber_products import (
    PRIMARY_ARTIFACT,
    published_accidentals,
    triage_specialization,
)
from search_fermigier_rank22_accidental_slices import (
    FERMIGIER_BIVARIATE_COEFFICIENTS,
    T0,
    conductor_probe,
    poly_add,
    poly_multiply,
    search_polynomial,
    trim_polynomial,
    quartic_group_pullback,
)


Q = Fraction
DEFAULT_SLOPES = (-3, -2, 0, 2, 3)


def projective_height(value: Fraction) -> int:
    value = Q(value)
    return max(abs(value.numerator), value.denominator)


def parameter_digest(parameters: list[Fraction]) -> str:
    digest = hashlib.sha256()
    for parameter in parameters:
        digest.update((rational_to_string(parameter) + "\n").encode())
    return digest.hexdigest()


@dataclass
class Incidence:
    source_label: str
    slice_id: str
    combination: str
    signed_parameter: Fraction
    quartic_x: Fraction
    quartic_y: Fraction

    def record(self) -> dict[str, Any]:
        return {
            "source_label": self.source_label,
            "slice_id": self.slice_id,
            "combination": self.combination,
            "signed_parameter": rational_to_string(self.signed_parameter),
            "quartic_x": rational_to_string(self.quartic_x),
            "quartic_y": rational_to_string(self.quartic_y),
        }


@dataclass(frozen=True)
class GenusTwoSlice:
    slice_id: str
    source_label: str
    source_point: tuple[Fraction, Fraction]
    slope: int
    intercept: Fraction
    coefficients: tuple[Fraction, ...]

    def original_x(self, parameter: Fraction) -> Fraction:
        return Q(self.slope)*Q(parameter)+self.intercept


def general_slice_polynomial(
    slope: int, intercept: Fraction
) -> tuple[Fraction, ...]:
    linear = (Q(intercept), Q(slope))
    x_power = (Q(1),)
    answer = (Q(0),)
    for coefficient_polynomial in FERMIGIER_BIVARIATE_COEFFICIENTS:
        answer = poly_add(
            answer, poly_multiply(coefficient_polynomial, x_power)
        )
        x_power = poly_multiply(x_power, linear)
    answer = trim_polynomial(answer)
    if len(answer)-1 != 6:
        raise AssertionError("an unused affine slope did not give degree six")
    return answer


def unused_slices(
    accidentals: tuple[tuple[str, tuple[Fraction, Fraction]], ...],
    slopes: tuple[int, ...],
):
    for label, source in accidentals:
        source_x, source_y = map(Q, source)
        for slope in slopes:
            intercept = source_x-Q(slope)*T0
            slope_label = f"m{abs(slope)}" if slope < 0 else f"p{slope}"
            yield GenusTwoSlice(
                slice_id=f"{label.lower()}_{slope_label}",
                source_label=label,
                source_point=(source_x, source_y),
                slope=slope,
                intercept=intercept,
                coefficients=general_slice_polynomial(slope, intercept),
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--slopes",
        default=",".join(map(str, DEFAULT_SLOPES)),
        help="comma-separated unused integral slopes",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "artifacts/local/elliptic-curves/"
            "fermigier-e22-multislope-collisions-v1.json"
        ),
    )
    parser.add_argument("--height", type=int, default=200_000)
    parser.add_argument("--timeout", type=float, default=10.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument("--triage", action="store_true")
    parser.add_argument("--conductor-timeout", type=float, default=20.0)
    parser.add_argument("--rank-timeout", type=float, default=30.0)
    arguments = parser.parse_args()
    slopes = tuple(dict.fromkeys(map(int, arguments.slopes.split(","))))
    if not slopes or any(slope in (-1, 1) for slope in slopes):
        raise SystemExit("supply at least one slope, excluding the prior +/-1 lanes")

    repository = Path(__file__).resolve().parents[2]
    output = arguments.output
    if not output.is_absolute():
        output = repository/output
    primary = json.loads(
        (repository/"artifacts/generated-results"/PRIMARY_ARTIFACT).read_text()
    )
    accidentals = published_accidentals(primary)

    parameters: dict[Fraction, list[Incidence]] = {}
    slice_records = []
    totals = Counter()
    for slice_data in unused_slices(accidentals, slopes):
        points, search = search_polynomial(
            slice_data.coefficients,
            height_bound=arguments.height,
            timeout=arguments.timeout,
            stack_bytes=arguments.stack_bytes,
        )
        accepted = 0
        seen_signed_x: set[tuple[Fraction, Fraction]] = set()
        for signed_parameter, ordinate in points:
            totals["searched_points"] += 1
            signed_parameter, ordinate = Q(signed_parameter), Q(ordinate)
            parameter = abs(signed_parameter)
            if parameter in (Q(0), abs(T0)):
                totals["source_or_zero"] += 1
                continue
            if FermigierMestreFamily.discriminant_factor(parameter) == 0:
                totals["singular"] += 1
                continue
            quartic_x = slice_data.original_x(signed_parameter)
            signed_key = signed_parameter, quartic_x
            if signed_key in seen_signed_x:
                totals["ordinate_sign_duplicate"] += 1
                continue
            seen_signed_x.add(signed_key)
            if ordinate**2 != FermigierMestreFamily.quartic_value(
                parameter, quartic_x
            ):
                raise AssertionError("a genus-two point missed the Fermigier quartic")
            generic_x = {
                point[0]
                for point in FermigierMestreFamily.known_quartic_points(parameter)
            }
            if quartic_x in generic_x:
                totals["generic_intersection"] += 1
                continue
            incidence = Incidence(
                source_label=slice_data.source_label,
                slice_id=slice_data.slice_id,
                combination="hyperellratpoints",
                signed_parameter=signed_parameter,
                quartic_x=quartic_x,
                quartic_y=ordinate,
            )
            parameters.setdefault(parameter, []).append(incidence)
            accepted += 1
            totals["accepted_incidences"] += 1
        slice_records.append(
            {
                "slice_id": slice_data.slice_id,
                "source_label": slice_data.source_label,
                "slope": slice_data.slope,
                "polynomial_degree": len(slice_data.coefficients)-1,
                "search": search,
                "accepted_incidence_count": accepted,
            }
        )
        print(
            f"MULTISLOPE|slice={slice_data.slice_id}|"
            f"points={len(points)}|search={search['status']}|"
            f"accepted={accepted}|unique_global={len(parameters)}",
            flush=True,
        )

    collision_records = []
    for parameter, incidences in parameters.items():
        by_source: dict[str, Incidence] = {}
        for incidence in incidences:
            by_source.setdefault(incidence.source_label, incidence)
        if len(by_source) < 2:
            continue
        collision_records.append(
            {
                "parameter": rational_to_string(parameter),
                "projective_height": projective_height(parameter),
                "distinct_source_count": len(by_source),
                "distinct_slice_count": len({row.slice_id for row in incidences}),
                "incidence_count": len(incidences),
                "source_representatives": [
                    incidence.record() for incidence in by_source.values()
                ],
            }
        )
    collision_records.sort(
        key=lambda row: (
            -row["distinct_source_count"],
            row["projective_height"],
            Q(row["parameter"]),
        )
    )
    ordered_parameters = sorted(
        parameters,
        key=lambda value: (projective_height(value), value),
    )
    triage_records = []
    if arguments.triage:
        for parameter in ordered_parameters:
            conductor = conductor_probe(
                parameter,
                timeout=arguments.conductor_timeout,
                stack_bytes=arguments.stack_bytes,
            )
            row: dict[str, Any] = {
                "parameter": rational_to_string(parameter),
                "conductor_probe": conductor,
            }
            if conductor.get("below_strict_log_conductor_target"):
                forced = []
                seen_x: set[Fraction] = set()
                for incidence in parameters[parameter]:
                    if incidence.quartic_x in seen_x:
                        continue
                    seen_x.add(incidence.quartic_x)
                    pullback = quartic_group_pullback(
                        parameter,
                        (incidence.quartic_x, incidence.quartic_y),
                    )
                    if pullback is None:
                        raise AssertionError("a multislope point pulled back to zero")
                    forced.append(
                        {
                            "quartic_x": rational_to_string(incidence.quartic_x),
                            "quartic_z": rational_to_string(incidence.quartic_y),
                            "basepoint_group_pullback": {
                                "jacobian_x": rational_to_string(pullback[0]),
                                "jacobian_y": rational_to_string(pullback[1]),
                            },
                        }
                    )
                row["rank_triage"] = triage_specialization(
                    {
                        "parameter_t": rational_to_string(parameter),
                        "forced_points": forced,
                    },
                    search_timeout=arguments.rank_timeout,
                    height_timeout=arguments.rank_timeout,
                    precisions=(72, 120),
                    stack_bytes=arguments.stack_bytes,
                    saturation_timeout=arguments.rank_timeout,
                    certificate_prime_bound=2_000,
                )
            triage_records.append(row)
            print(
                f"MULTISLOPE|triage_T={parameter}|"
                f"logN={conductor.get('log_conductor')}|"
                f"rank={(row.get('rank_triage') or {}).get('full_pool_stable_numerical_rank')}",
                flush=True,
            )
    payload = {
        "schema_version": 1,
        "artifact_kind": "bounded_fermigier_e22_multislope_collision_search",
        "status": (
            "bounded_search_with_distinct_source_collisions"
            if collision_records
            else "bounded_search_no_distinct_source_collision"
        ),
        "claim_scope": {
            "exact": "slice identities, seed group law, inverse maps, and Fermigier quartic membership",
            "bounded": "the declared slopes and signed support-at-most-two seed combinations only",
            "not_claimed": "Mordell-Weil completeness, section independence, or rank 21",
        },
        "parameters": {
            "slopes": list(slopes),
            "source_count": len(accidentals),
            "height_bound": arguments.height,
            "timeout_seconds_per_slice": arguments.timeout,
            "search": "PARI hyperellratpoints on degree-six affine slices",
        },
        "slice_records": slice_records,
        "totals": {
            **dict(totals),
            "slice_count": len(slice_records),
            "unique_parameter_count": len(parameters),
            "unique_parameter_sha256": parameter_digest(ordered_parameters),
            "distinct_source_collision_count": len(collision_records),
            "maximum_distinct_source_count": max(
                (row["distinct_source_count"] for row in collision_records),
                default=1,
            ),
        },
        "candidate_parameters": [
            {
                "parameter": rational_to_string(parameter),
                "projective_height": projective_height(parameter),
                "incidences": [row.record() for row in parameters[parameter]],
            }
            for parameter in ordered_parameters
        ],
        "collisions": collision_records,
        "triage": triage_records,
        "reproduction": {
            "command": (
                "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
                "elliptic-curves/cas/search_fermigier_rank22_multislope_collisions.py"
            )
        },
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
    print(
        "MULTISLOPE|"
        f"slices={len(slice_records)}|searched_points={totals['searched_points']}|"
        f"unique_parameters={len(parameters)}|collisions={len(collision_records)}|"
        f"max_sources={payload['totals']['maximum_distinct_source_count']}|"
        f"output={output}|status=PASS_BOUNDED",
        flush=True,
    )


if __name__ == "__main__":
    main()
