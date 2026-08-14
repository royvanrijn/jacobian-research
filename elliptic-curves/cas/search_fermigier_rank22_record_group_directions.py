#!/usr/bin/env python3
"""Manufacture Fermigier slices from the record quartic's exact group law.

The rank-22 record fibre at ``T0=39508/39`` is itself a pointed genus-one
quartic.  Taking the first visible quartic point as origin gives an explicit
birational map to a generalized Weierstrass curve.  This script transports
the canonical published rank-22 basis through that map, then enumerates the
global-sign quotient of coefficient vectors in ``{-1,0,1}^22`` of l1 norm at
most two.

The 22 singleton vectors are exact transport calibrations.  Every weight-two
vector produces a new record-fibre quartic abscissa, after excluding both the
27 abscissas in the H=10^6 replay and all 22 published preimages.  For every
remaining abscissa ``x0`` the two slices

    x = T + (x0-T0),       x = -T + (x0+T0)

are searched once at parameter height 50,000.  Every prior Fermigier
parameter in the pinned manifest, the earlier auxiliary-orbit tranche, and
the H=50,000 pair-product stage is removed.  Exact conductors precede any
rank work; only completed conductors below 182.72 receive H=50,000 point and
height-rank triage.  Stable numerical rank at least 21 immediately triggers
the existing saturation and finite-reduction certificate path.

This is a finite group-orbit and slice-height computation, not a complete
rational-point or Mordell--Weil enumeration.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import shlex
import subprocess
import sys
import time
from typing import Any, Iterable, Sequence

import sympy as sp

from ek_k3 import rational_to_string
from fermigier_mestre import FermigierMestreFamily
from pari_bridge import pari_version
from search_fermigier_published_pair_fiber_products import (
    EXPECTED_PUBLISHED_PREIMAGE_SHA256,
    PRIMARY_ARTIFACT,
    extract_parameter_values,
    generic_abscissas,
    polynomial_digest,
    published_preimage_digest,
    rational_digest,
    sha256_file,
    triage_specialization,
)
from search_fermigier_published_pair_fiber_products_h50000 import (
    EXPECTED_H50000_RESULT_SHA256,
    h50000_result_digest,
)
from search_fermigier_rank22_accidental_slices import (
    T0,
    canonical_signless_points,
    conductor_probe,
    point_record,
    poly_evaluate,
    quartic_group_pullback,
    search_polynomial,
    short_add,
    short_negate,
    slice_polynomial,
)
from search_fermigier_rank22_auxiliary_orbits import prior_parameter_manifest
from search_nagao_section7_auxiliary_jacobians import (
    translate_polynomial,
    weierstrass_add,
    weierstrass_multiply,
)
from triage_nagao_rank13_finalists import point_digest, point_on_short_curve


Q = Fraction
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

SLICE_HEIGHT = 50_000
EXPECTED_BASIS_COUNT = 22
EXPECTED_ORBIT_VECTOR_COUNT = 484
EXPECTED_NEW_DIRECTION_COUNT = 462
EXPECTED_SLICE_CALL_COUNT = 924
EXPECTED_TRANSPORT_SOURCE_SHA256 = (
    "8650d548b7c101514d3403a46aa8ea49cef071b231d86df2e3950be84a868011"
)
EXPECTED_PUBLISHED_BASIS_POINT_SHA256 = (
    "4b2b89f3f432be8599b6ab5109c1221af52d1aa96c05eba1229983d929c9e727"
)
EXPECTED_RECORD_X_SHA256 = (
    "57dea60915bd86a6594dddc7691ef17d0812d74132ff197d2162f87bd7535ff9"
)
EXPECTED_PUBLISHED_PREIMAGE_X_SHA256 = (
    "19a7f4b94ba385bb2f8d71d49d43c50558f5b9bfb1d400edd713b9e15f754c6d"
)
EXPECTED_KNOWN_RECORD_X_UNION_SHA256 = (
    "a692902016852128dac0997705d5463a4d94f48a0c4089d324f0c110f8ac8167"
)
EXPECTED_AUXILIARY_ORBIT_ARTIFACT_SHA256 = (
    "1008336232ac65bb2bace6ff7008ffb18c1b491c6f0295e16fda14949d5b94d6"
)
EXPECTED_PRIOR_PARAMETER_COUNT = 1_239
EXPECTED_PRIOR_PARAMETER_SHA256 = (
    "9482e61650aa8bb1fd45c3765e5db92c1474090faee8d831e0d73cee4fc864c4"
)
AUXILIARY_ORBIT_ARTIFACT = "elliptic_fermigier_rank22_auxiliary_orbits.json"
H50000_PAIR_ARTIFACT = (
    "elliptic_fermigier_published_pair_fiber_products_h50000.json"
)


def projective_height(value: Fraction) -> int:
    value = Q(value)
    return max(abs(value.numerator), value.denominator)


def fraction_pair(value: dict[str, Any], prefix: str) -> tuple[Fraction, Fraction]:
    return Q(value[f"{prefix}_x"]), Q(value[f"{prefix}_y"])


def transport_source_digest(rows: Sequence[dict[str, Any]]) -> str:
    digest = hashlib.sha256()
    for row in rows:
        digest.update(
            (
                f"{row['label']}|{Q(row['quartic_preimage']['x'])}|"
                f"{Q(row['quartic_preimage']['z'])}|"
                f"{Q(row['short_point']['jacobian_x'])}|"
                f"{Q(row['short_point']['jacobian_y'])}\n"
            ).encode()
        )
    return digest.hexdigest()


def point_pair_record(point: tuple[Fraction, Fraction]) -> dict[str, str | bool]:
    return {
        "x": rational_to_string(point[0]),
        "y": rational_to_string(point[1]),
        "exact_membership_checked": True,
    }


@dataclass(frozen=True)
class RecordQuarticAuxiliary:
    base_point: tuple[Fraction, Fraction]
    quartic_coefficients: tuple[Fraction, ...]
    shifted_coefficients: tuple[Fraction, ...]
    weierstrass_coefficients: tuple[Fraction, ...]

    @classmethod
    def construct(cls) -> "RecordQuarticAuxiliary":
        base = tuple(Q(value) for value in FermigierMestreFamily.known_quartic_points(T0)[0])
        coefficients = tuple(
            reversed(FermigierMestreFamily.quartic_coefficients(T0))
        )
        shifted = translate_polynomial(coefficients, base[0])
        if len(shifted) != 5 or shifted[0] != base[1] ** 2:
            raise AssertionError("the selected record origin missed the quartic")
        _, d_value, c_value, b_value, a_value = shifted
        q_value = base[1]
        weierstrass = (
            d_value / q_value,
            c_value - d_value**2 / (4 * q_value**2),
            2 * q_value * b_value,
            -4 * q_value**2 * a_value,
            a_value * (d_value**2 - 4 * q_value**2 * c_value),
        )
        return cls(base, coefficients, shifted, weierstrass)

    def quartic_value(self, x_value: Fraction) -> Fraction:
        answer = Q(0)
        for coefficient in reversed(self.quartic_coefficients):
            answer = answer * Q(x_value) + coefficient
        return answer

    def forward(
        self, point: tuple[Fraction, Fraction]
    ) -> tuple[Fraction, Fraction]:
        x_value, ordinate = (Q(value) for value in point)
        if ordinate**2 != self.quartic_value(x_value):
            raise ValueError("the point missed the record quartic")
        u_value = x_value - self.base_point[0]
        if u_value == 0:
            raise ValueError("the pointed record origin maps to infinity")
        _, d_value, c_value, _, _ = self.shifted_coefficients
        q_value = self.base_point[1]
        image_x = (
            2 * q_value * (ordinate + q_value) + d_value * u_value
        ) / u_value**2
        image_y = (
            4 * q_value**2 * (ordinate + q_value)
            + 2 * q_value * (d_value * u_value + c_value * u_value**2)
            - d_value**2 * u_value**2 / (2 * q_value)
        ) / u_value**3
        image = image_x, image_y
        if not point_on_short_curve(self.weierstrass_coefficients, image):
            raise AssertionError("the pointed-quartic forward map failed exactly")
        return image

    def inverse(
        self, point: tuple[Fraction, Fraction] | None
    ) -> tuple[Fraction, Fraction] | None:
        if point is None:
            return self.base_point
        if not point_on_short_curve(self.weierstrass_coefficients, point):
            raise ValueError("the point missed the auxiliary Weierstrass curve")
        image_x, image_y = point
        if image_y == 0:
            return None
        _, d_value, c_value, _, _ = self.shifted_coefficients
        q_value = self.base_point[1]
        u_value = (
            4 * q_value**2 * (image_x + c_value) - d_value**2
        ) / (2 * q_value * image_y)
        if u_value == 0:
            return None
        ordinate = (
            (image_x * u_value**2 - d_value * u_value) / (2 * q_value)
            - q_value
        )
        answer = self.base_point[0] + u_value, ordinate
        if ordinate**2 != self.quartic_value(answer[0]):
            raise AssertionError("the pointed-quartic inverse map failed exactly")
        return answer


def orbit_vectors() -> tuple[tuple[int, ...], ...]:
    vectors: list[tuple[int, ...]] = []
    for first in range(EXPECTED_BASIS_COUNT):
        vector = [0] * EXPECTED_BASIS_COUNT
        vector[first] = 1
        vectors.append(tuple(vector))
    for first in range(EXPECTED_BASIS_COUNT):
        for second in range(first + 1, EXPECTED_BASIS_COUNT):
            for second_sign in (-1, 1):
                vector = [0] * EXPECTED_BASIS_COUNT
                vector[first] = 1
                vector[second] = second_sign
                vectors.append(tuple(vector))
    if len(vectors) != EXPECTED_ORBIT_VECTOR_COUNT:
        raise AssertionError("the sign-quotiented weight-two orbit changed")
    return tuple(vectors)


def vector_id(vector: Sequence[int]) -> str:
    terms = []
    for index, coefficient in enumerate(vector, start=1):
        if coefficient:
            terms.append(f"{'p' if coefficient > 0 else 'm'}{index:02d}")
    return "_".join(terms)


def auxiliary_combination(
    coefficients: Sequence[Fraction],
    basis: Sequence[tuple[Fraction, Fraction]],
    vector: Sequence[int],
) -> tuple[Fraction, Fraction] | None:
    answer = None
    for point, scalar in zip(basis, vector, strict=True):
        if scalar:
            answer = weierstrass_add(
                coefficients,
                answer,
                weierstrass_multiply(coefficients, point, scalar),
            )
    return answer


def short_combination(
    coefficients: Sequence[Fraction],
    basis: Sequence[tuple[Fraction, Fraction]],
    vector: Sequence[int],
) -> tuple[Fraction, Fraction] | None:
    answer = None
    for point, scalar in zip(basis, vector, strict=True):
        if scalar not in (-1, 0, 1):
            raise ValueError("this orbit only supports coefficients -1,0,1")
        if scalar:
            answer = short_add(
                coefficients,
                answer,
                point if scalar == 1 else short_negate(point),
            )
    return answer


def direction_digest(records: Iterable[dict[str, Any]]) -> str:
    digest = hashlib.sha256()
    for record in sorted(records, key=lambda item: item["direction_id"]):
        digest.update(
            (
                f"{record['direction_id']}|{record['quartic_x']}|"
                f"{record['quartic_z']}|"
                f"{','.join(map(str, record['coefficient_vector']))}\n"
            ).encode()
        )
    return digest.hexdigest()


def load_transport_source(
    primary: dict[str, Any], auxiliary: RecordQuarticAuxiliary
) -> tuple[
    tuple[tuple[Fraction, Fraction], ...],
    tuple[tuple[Fraction, Fraction], ...],
    list[dict[str, Any]],
]:
    rows = primary["published_point_preimages"]
    if len(rows) != EXPECTED_BASIS_COUNT or [row["label"] for row in rows] != [
        f"P{index}" for index in range(1, 23)
    ]:
        raise AssertionError("the canonical published preimage ordering changed")
    if transport_source_digest(rows) != EXPECTED_TRANSPORT_SOURCE_SHA256:
        raise AssertionError("the exact published transport source changed")
    certificate = primary["record_fiber_height_1000000_replay"][
        "exact_rank_certificate"
    ]["canonical_published_point_finite_reduction_certificate"]
    if (
        certificate["status"] != "certified"
        or certificate["combined_exact_rank_over_F2"] != 22
        or certificate["point_sha256"]
        != EXPECTED_PUBLISHED_BASIS_POINT_SHA256
    ):
        raise AssertionError("the canonical published rank-22 certificate changed")

    auxiliary_basis = []
    short_basis = []
    records = []
    for row in rows:
        quartic_point = (
            Q(row["quartic_preimage"]["x"]),
            Q(row["quartic_preimage"]["z"]),
        )
        short_point = (
            Q(row["short_point"]["jacobian_x"]),
            Q(row["short_point"]["jacobian_y"]),
        )
        if quartic_group_pullback(T0, quartic_point) != short_point:
            raise AssertionError("a published preimage lost its short group image")
        image = auxiliary.forward(quartic_point)
        if auxiliary.inverse(image) != quartic_point:
            raise AssertionError("a published basis point failed the birational round trip")
        auxiliary_basis.append(image)
        short_basis.append(short_point)
        records.append(
            {
                "label": row["label"],
                "quartic_preimage": {
                    "x": rational_to_string(quartic_point[0]),
                    "z": rational_to_string(quartic_point[1]),
                },
                "short_group_point": point_record(short_point),
                "auxiliary_group_point": point_pair_record(image),
                "exact_forward_inverse_round_trip_checked": True,
                "exact_short_group_coordinate_checked": True,
            }
        )
    if point_digest(short_basis) != EXPECTED_PUBLISHED_BASIS_POINT_SHA256:
        raise AssertionError("the transported published point digest changed")
    return tuple(auxiliary_basis), tuple(short_basis), records


def known_record_abscissas(primary: dict[str, Any]) -> tuple[set[Fraction], dict[str, Any]]:
    exact_certificate = primary["record_fiber_height_1000000_replay"][
        "exact_rank_certificate"
    ]
    replay_x = {
        Q(record["quartic_x"])
        for record in exact_certificate["record_signless_points"]
    }
    published_x = {
        Q(record["quartic_preimage"]["x"])
        for record in primary["published_point_preimages"]
    }
    if (
        len(replay_x) != 27
        or rational_digest(sorted(replay_x)) != EXPECTED_RECORD_X_SHA256
        or len(published_x) != 22
        or rational_digest(sorted(published_x))
        != EXPECTED_PUBLISHED_PREIMAGE_X_SHA256
        or len(replay_x | published_x) != 32
        or rational_digest(sorted(replay_x | published_x))
        != EXPECTED_KNOWN_RECORD_X_UNION_SHA256
    ):
        raise AssertionError("the exact known record-fibre abscissas changed")
    return replay_x | published_x, {
        "H1000000_abscissa_count": len(replay_x),
        "H1000000_abscissa_sha256": rational_digest(sorted(replay_x)),
        "published_preimage_abscissa_count": len(published_x),
        "published_preimage_abscissa_sha256": rational_digest(sorted(published_x)),
        "union_count": len(replay_x | published_x),
        "union_sha256": rational_digest(sorted(replay_x | published_x)),
    }


def generate_directions(
    auxiliary: RecordQuarticAuxiliary,
    auxiliary_basis: Sequence[tuple[Fraction, Fraction]],
    short_basis: Sequence[tuple[Fraction, Fraction]],
    known_x: set[Fraction],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    vectors = orbit_vectors()
    short_coefficients = FermigierMestreFamily.coefficients(T0)
    all_x: set[Fraction] = set()
    directions = []
    exceptional = 0
    known_exclusions = 0
    for vector in vectors:
        auxiliary_point = auxiliary_combination(
            auxiliary.weierstrass_coefficients, auxiliary_basis, vector
        )
        inverse = auxiliary.inverse(auxiliary_point)
        if inverse is None:
            exceptional += 1
            continue
        x_value, ordinate = inverse
        expected_short = short_combination(short_coefficients, short_basis, vector)
        if quartic_group_pullback(T0, inverse) != expected_short:
            raise AssertionError("an orbit inverse lost its exact short group coordinate")
        if x_value in all_x:
            raise AssertionError("two sign-quotiented orbit vectors shared an abscissa")
        all_x.add(x_value)
        if x_value in known_x:
            known_exclusions += 1
            continue
        if sum(abs(value) for value in vector) != 2:
            raise AssertionError("a new direction did not come from a pair combination")
        directions.append(
            {
                "direction_id": vector_id(vector),
                "coefficient_vector": list(vector),
                "quartic_x": rational_to_string(x_value),
                "quartic_z": rational_to_string(ordinate),
                "projective_height": projective_height(x_value),
                "exact_auxiliary_inverse_checked": True,
                "exact_short_group_combination_checked": True,
            }
        )
    directions.sort(
        key=lambda record: (
            record["projective_height"],
            Q(record["quartic_x"]),
            record["direction_id"],
        )
    )
    if (
        exceptional != 0
        or known_exclusions != 22
        or len(all_x) != EXPECTED_ORBIT_VECTOR_COUNT
        or len(directions) != EXPECTED_NEW_DIRECTION_COUNT
    ):
        raise AssertionError("the exact low-weight record orbit changed")
    return directions, {
        "global_sign_quotient": True,
        "coefficient_alphabet": [-1, 0, 1],
        "maximum_l1_norm": 2,
        "singleton_transport_calibrations": 22,
        "sign_quotiented_vector_count": len(vectors),
        "exceptional_inverse_count": exceptional,
        "known_abscissa_exclusions": known_exclusions,
        "unique_orbit_abscissa_count": len(all_x),
        "genuinely_new_pair_combination_abscissa_count": len(directions),
        "new_direction_sha256": direction_digest(directions),
        "minimum_new_abscissa_projective_height": min(
            record["projective_height"] for record in directions
        ),
        "maximum_new_abscissa_projective_height": max(
            record["projective_height"] for record in directions
        ),
    }


def load_prior_parameters(
    artifact_directory: Path,
) -> tuple[set[Fraction], dict[str, Any]]:
    base, base_record = prior_parameter_manifest(artifact_directory)
    auxiliary_path = artifact_directory / AUXILIARY_ORBIT_ARTIFACT
    if sha256_file(auxiliary_path) != EXPECTED_AUXILIARY_ORBIT_ARTIFACT_SHA256:
        raise AssertionError("the earlier auxiliary-orbit artifact changed")
    auxiliary_values = extract_parameter_values(json.loads(auxiliary_path.read_text()))

    pair_path = artifact_directory / H50000_PAIR_ARTIFACT
    pair_artifact = json.loads(pair_path.read_text())
    pair_result_digest = h50000_result_digest(pair_artifact["pair_searches"])
    if pair_result_digest != EXPECTED_H50000_RESULT_SHA256:
        raise AssertionError("the exact H=50000 pair-product result changed")
    pair_seen = {
        abs(Q(incidence["canonical_parameter_t"]))
        for row in pair_artifact["pair_searches"]
        for incidence in row["search"]["incidences"]
    }
    parameters = base | auxiliary_values | pair_seen
    digest = rational_digest(sorted(parameters))
    if (
        len(parameters) != EXPECTED_PRIOR_PARAMETER_COUNT
        or digest != EXPECTED_PRIOR_PARAMETER_SHA256
    ):
        raise AssertionError("the terminal prior Fermigier parameter set changed")
    return parameters, {
        "base_manifest_parameter_count": len(base),
        "base_manifest_parameter_sha256": base_record["parameter_sha256"],
        "base_manifest": base_record,
        "auxiliary_orbit_artifact": AUXILIARY_ORBIT_ARTIFACT,
        "auxiliary_orbit_artifact_sha256": EXPECTED_AUXILIARY_ORBIT_ARTIFACT_SHA256,
        "auxiliary_orbit_extracted_parameter_count": len(auxiliary_values),
        "H50000_pair_artifact": H50000_PAIR_ARTIFACT,
        "H50000_exact_pair_result_sha256": pair_result_digest,
        "H50000_seen_parameters": [
            rational_to_string(value) for value in sorted(pair_seen)
        ],
        "terminal_prior_parameter_count": len(parameters),
        "terminal_prior_parameter_sha256": digest,
    }


def classify_slice_point(
    *,
    direction: dict[str, Any],
    slope: int,
    intercept: Fraction,
    signed_parameter: Fraction,
    ordinate: Fraction,
    prior_parameters: set[Fraction],
) -> dict[str, Any]:
    signed_parameter = Q(signed_parameter)
    canonical_parameter = abs(signed_parameter)
    x_value = Q(slope) * signed_parameter + Q(intercept)
    record: dict[str, Any] = {
        "signed_parameter_t": rational_to_string(signed_parameter),
        "canonical_parameter_t": rational_to_string(canonical_parameter),
        "quartic_x": rational_to_string(x_value),
        "quartic_z": rational_to_string(abs(ordinate)),
    }
    if ordinate**2 != FermigierMestreFamily.quartic_value(
        canonical_parameter, x_value
    ):
        raise AssertionError("a returned slice point missed the canonical fibre")
    if canonical_parameter == abs(T0):
        record["classification"] = "record-source-calibration-excluded"
        return record
    if canonical_parameter == 0:
        record["classification"] = "zero-parameter-excluded"
        return record
    if FermigierMestreFamily.discriminant_factor(canonical_parameter) == 0:
        record["classification"] = "singular-fibre-excluded"
        return record
    if canonical_parameter in prior_parameters:
        record["classification"] = "prior-parameter-excluded"
        return record
    if x_value in generic_abscissas(canonical_parameter):
        record["classification"] = "generic-section-collision-excluded"
        return record
    quartic_point = x_value, abs(ordinate)
    pullback = quartic_group_pullback(canonical_parameter, quartic_point)
    if pullback is None:
        raise AssertionError("a nongeneric forced point mapped to the origin")
    record.update(
        {
            "classification": "genuinely-new-forced-fibre",
            "source_direction_id": direction["direction_id"],
            "source_record_quartic_x": direction["quartic_x"],
            "slice_slope": slope,
            "slice_intercept": rational_to_string(intercept),
            "basepoint_group_pullback": point_record(pullback),
            "exact_quartic_membership_checked": True,
        }
    )
    return record


def search_direction(
    direction: dict[str, Any],
    *,
    height_bound: int,
    timeout: float,
    stack_bytes: int,
    prior_parameters: set[Fraction],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    source_x = Q(direction["quartic_x"])
    source_z = Q(direction["quartic_z"])
    rows = []
    qualifying = []
    for slope in (-1, 1):
        intercept = source_x - slope * T0
        coefficients = slice_polynomial(slope, intercept)
        if poly_evaluate(coefficients, T0) != source_z**2:
            raise AssertionError("a group-manufactured slice missed its source point")
        raw_points, search = search_polynomial(
            coefficients,
            height_bound=height_bound,
            timeout=timeout,
            stack_bytes=stack_bytes,
        )
        incidences = [
            classify_slice_point(
                direction=direction,
                slope=slope,
                intercept=intercept,
                signed_parameter=parameter,
                ordinate=ordinate,
                prior_parameters=prior_parameters,
            )
            for parameter, ordinate in canonical_signless_points(raw_points)
        ]
        calibrations = [
            incidence
            for incidence in incidences
            if incidence["classification"]
            == "record-source-calibration-excluded"
        ]
        if search["status"] == "completed" and len(calibrations) != 1:
            raise AssertionError("a completed slice lost its exact T0 calibration")
        qualifying.extend(
            incidence
            for incidence in incidences
            if incidence["classification"] == "genuinely-new-forced-fibre"
        )
        rows.append(
            {
                "slice_id": f"{direction['direction_id']}_{'m1' if slope == -1 else 'p1'}",
                "slope": slope,
                "intercept": rational_to_string(intercept),
                "quartic_polynomial_sha256": polynomial_digest(coefficients),
                "search": search,
                "record_T0_calibration_count": len(calibrations),
                "incidences": incidences,
                "qualifying_new_parameter_count": len(
                    {row["canonical_parameter_t"] for row in qualifying}
                ),
            }
        )
    return {
        **direction,
        "slice_searches": rows,
    }, qualifying


def aggregate_candidates(
    qualifying: Sequence[dict[str, Any]],
) -> list[dict[str, Any]]:
    candidates: dict[Fraction, dict[str, Any]] = {}
    for incidence in qualifying:
        parameter = Q(incidence["canonical_parameter_t"])
        candidate = candidates.setdefault(
            parameter,
            {
                "parameter_t": rational_to_string(parameter),
                "signed_parameters": set(),
                "source_direction_ids": set(),
                "forced_points_by_x": {},
            },
        )
        candidate["signed_parameters"].add(incidence["signed_parameter_t"])
        candidate["source_direction_ids"].add(incidence["source_direction_id"])
        point = {
            "quartic_x": incidence["quartic_x"],
            "quartic_z": incidence["quartic_z"],
            "source_direction_id": incidence["source_direction_id"],
            "source_record_quartic_x": incidence["source_record_quartic_x"],
            "slice_slope": incidence["slice_slope"],
            "slice_intercept": incidence["slice_intercept"],
            "basepoint_group_pullback": incidence["basepoint_group_pullback"],
            "exact_membership_checked": True,
        }
        existing = candidate["forced_points_by_x"].get(point["quartic_x"])
        if existing is None:
            candidate["forced_points_by_x"][point["quartic_x"]] = point
        elif existing["quartic_z"] != point["quartic_z"]:
            raise AssertionError("one forced abscissa acquired inconsistent ordinates")
    answer = []
    for parameter, candidate in sorted(
        candidates.items(), key=lambda item: (projective_height(item[0]), item[0])
    ):
        points = sorted(
            candidate["forced_points_by_x"].values(),
            key=lambda point: Q(point["quartic_x"]),
        )
        answer.append(
            {
                "parameter_t": candidate["parameter_t"],
                "projective_height": projective_height(parameter),
                "signed_parameters": sorted(candidate["signed_parameters"]),
                "source_direction_ids": sorted(candidate["source_direction_ids"]),
                "distinct_source_direction_count": len(
                    candidate["source_direction_ids"]
                ),
                "distinct_forced_quartic_abscissa_count": len(points),
                "forced_points": points,
            }
        )
    return answer


def slice_result_digest(direction_rows: Sequence[dict[str, Any]]) -> str:
    digest = hashlib.sha256()
    for direction in direction_rows:
        for row in direction["slice_searches"]:
            search = row["search"]
            digest.update(
                (
                    f"{row['slice_id']}|{row['quartic_polynomial_sha256']}|"
                    f"{search['status']}|{search.get('signed_point_count')}|"
                    f"{row['record_T0_calibration_count']}|"
                    f"{len(row['incidences'])}\n"
                ).encode()
            )
    return digest.hexdigest()


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--slice-height", type=int, default=SLICE_HEIGHT)
    parser.add_argument("--slice-timeout", type=float, default=15.0)
    parser.add_argument("--conductor-timeout", type=float, default=15.0)
    parser.add_argument("--specialization-timeout", type=float, default=30.0)
    parser.add_argument("--height-timeout", type=float, default=20.0)
    parser.add_argument("--saturation-timeout", type=float, default=30.0)
    parser.add_argument("--certificate-prime-bound", type=int, default=2_000)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=root
        / "artifacts"
        / "generated-results"
        / "elliptic_fermigier_rank22_record_group_directions.json",
    )
    return parser


def write_artifact(path: Path, artifact: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")


def main() -> None:
    args = build_parser().parse_args()
    if args.slice_height != SLICE_HEIGHT:
        raise SystemExit("this bounded tranche is pinned at slice H=50000")
    if min(
        args.slice_timeout,
        args.conductor_timeout,
        args.specialization_timeout,
        args.height_timeout,
        args.saturation_timeout,
    ) <= 0 or max(
        args.slice_timeout,
        args.conductor_timeout,
        args.specialization_timeout,
        args.height_timeout,
        args.saturation_timeout,
    ) > 60:
        raise SystemExit("all subprocess caps must lie in (0,60]")

    root = Path(__file__).resolve().parents[2]
    artifact_directory = root / "artifacts" / "generated-results"
    primary_path = artifact_directory / PRIMARY_ARTIFACT
    primary = json.loads(primary_path.read_text())
    if published_preimage_digest(primary) != EXPECTED_PUBLISHED_PREIMAGE_SHA256:
        raise AssertionError("the exact published accidental preimages changed")

    auxiliary = RecordQuarticAuxiliary.construct()
    auxiliary_basis, short_basis, basis_records = load_transport_source(
        primary, auxiliary
    )
    known_x, known_record = known_record_abscissas(primary)
    directions, orbit_record = generate_directions(
        auxiliary, auxiliary_basis, short_basis, known_x
    )
    prior_parameters, prior_record = load_prior_parameters(artifact_directory)

    script_path = Path(__file__).resolve()
    artifact: dict[str, Any] = {
        "schema_version": 1,
        "status": "in-progress exact record-group manufactured-direction search",
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": "182.72",
            "hit": False,
        },
        "source": {
            "primary_artifact": PRIMARY_ARTIFACT,
            "primary_artifact_sha256_observed": sha256_file(primary_path),
            "published_accidental_preimage_sha256": published_preimage_digest(primary),
            "transport_source_sha256": transport_source_digest(
                primary["published_point_preimages"]
            ),
            "certified_published_rank22_basis_point_sha256": point_digest(
                short_basis
            ),
            "known_record_abscissas": known_record,
        },
        "pointed_quartic_map": {
            "record_parameter_t": rational_to_string(T0),
            "origin_quartic_point": {
                "x": rational_to_string(auxiliary.base_point[0]),
                "z": rational_to_string(auxiliary.base_point[1]),
            },
            "quartic_coefficients_ascending_in_x": [
                rational_to_string(value)
                for value in auxiliary.quartic_coefficients
            ],
            "shifted_coefficients_ascending_in_u_x_minus_origin": [
                rational_to_string(value)
                for value in auxiliary.shifted_coefficients
            ],
            "generalized_weierstrass_coefficients_a1_a2_a3_a4_a6": [
                rational_to_string(value)
                for value in auxiliary.weierstrass_coefficients
            ],
            "formula": (
                "for f(x0+u)=q^2+d*u+c*u^2+b*u^3+a*u^4, "
                "X=(2q(z+q)+du)/u^2 and "
                "Y=(4q^2(z+q)+2q(du+cu^2)-d^2u^2/(2q))/u^3"
            ),
            "exact_inverse_formula_checked_for_every_basis_and_orbit_point": True,
        },
        "transported_rank22_basis": basis_records,
        "orbit_generation": orbit_record,
        "prior_decontamination": prior_record,
        "parameters": {
            "coefficient_alphabet": [-1, 0, 1],
            "maximum_l1_norm": 2,
            "global_sign_quotient": True,
            "new_direction_count": len(directions),
            "slopes_per_direction": [-1, 1],
            "declared_slice_call_count": len(directions) * 2,
            "slice_height": args.slice_height,
            "slice_timeout_seconds": args.slice_timeout,
            "conductor_timeout_seconds": args.conductor_timeout,
            "specialization_height": SLICE_HEIGHT,
            "specialization_timeout_seconds": args.specialization_timeout,
            "height_timeout_seconds": args.height_timeout,
            "height_precisions": [72, 120],
            "saturation_timeout_seconds": args.saturation_timeout,
            "certificate_prime_bound": args.certificate_prime_bound,
            "stack_bytes": args.stack_bytes,
            "no_retries": True,
            "checkpoint_after_each_direction": True,
        },
        "direction_searches": [],
        "candidates": [],
        "execution": {
            "phase": "slice-search-in-progress",
            "directions_completed": 0,
            "slice_calls_completed_or_attempted": 0,
        },
        "software": {
            "python": platform.python_version(),
            "sympy": sp.__version__,
            "pari_gp": pari_version(),
        },
        "reproducing_command": " ".join(
            shlex.quote(part) for part in [sys.executable, *sys.argv]
        ),
        "script_sha256": sha256_file(script_path),
    }
    write_artifact(args.output, artifact)

    qualifying_all: list[dict[str, Any]] = []
    started = time.monotonic()
    for index, direction in enumerate(directions, start=1):
        direction_row, qualifying = search_direction(
            direction,
            height_bound=args.slice_height,
            timeout=args.slice_timeout,
            stack_bytes=args.stack_bytes,
            prior_parameters=prior_parameters,
        )
        artifact["direction_searches"].append(direction_row)
        qualifying_all.extend(qualifying)
        artifact["execution"].update(
            {
                "directions_completed": index,
                "slice_calls_completed_or_attempted": index * 2,
                "last_direction_id": direction["direction_id"],
                "wall_seconds_so_far": time.monotonic() - started,
            }
        )
        write_artifact(args.output, artifact)
        if index % 25 == 0 or index == len(directions):
            print(
                f"directions {index}/{len(directions)} "
                f"qualifying_incidences={len(qualifying_all)}",
                flush=True,
            )

    candidates = aggregate_candidates(qualifying_all)
    artifact["candidates"] = candidates
    artifact["execution"]["phase"] = "conductor-first"
    write_artifact(args.output, artifact)
    for candidate in candidates:
        parameter = Q(candidate["parameter_t"])
        candidate["conductor_probe"] = conductor_probe(
            parameter,
            timeout=args.conductor_timeout,
            stack_bytes=args.stack_bytes,
        )
        write_artifact(args.output, artifact)
        if candidate["conductor_probe"].get(
            "below_strict_log_conductor_target"
        ):
            try:
                candidate["rank_triage"] = triage_specialization(
                    candidate,
                    search_timeout=args.specialization_timeout,
                    height_timeout=args.height_timeout,
                    precisions=(72, 120),
                    stack_bytes=args.stack_bytes,
                    saturation_timeout=args.saturation_timeout,
                    certificate_prime_bound=args.certificate_prime_bound,
                )
            except subprocess.TimeoutExpired as error:
                candidate["rank_triage"] = {
                    "status": "timeout-no-retry",
                    "error": str(error)[:1000],
                }
            except (RuntimeError, AssertionError, ValueError) as error:
                candidate["rank_triage"] = {
                    "status": "error-no-retry",
                    "error": str(error)[:1000],
                }
            write_artifact(args.output, artifact)

    slice_rows = [
        row
        for direction in artifact["direction_searches"]
        for row in direction["slice_searches"]
    ]
    classifications = Counter(
        incidence["classification"]
        for row in slice_rows
        for incidence in row["incidences"]
    )
    rank_records = [
        candidate["rank_triage"]
        for candidate in candidates
        if "rank_triage" in candidate
        and "full_pool_stable_numerical_rank" in candidate["rank_triage"]
    ]
    artifact["outcome"] = {
        "declared_direction_count": EXPECTED_NEW_DIRECTION_COUNT,
        "declared_slice_call_count": EXPECTED_SLICE_CALL_COUNT,
        "directions_completed": len(artifact["direction_searches"]),
        "slice_calls_attempted": len(slice_rows),
        "slice_calls_completed": sum(
            row["search"]["status"] == "completed" for row in slice_rows
        ),
        "slice_calls_timed_out_or_errored": sum(
            row["search"]["status"] != "completed" for row in slice_rows
        ),
        "record_T0_calibrated_slices": sum(
            row["record_T0_calibration_count"] == 1 for row in slice_rows
        ),
        "incidence_classification_counts": dict(sorted(classifications.items())),
        "genuinely_new_forced_fibres": len(candidates),
        "completed_conductors": sum(
            candidate.get("conductor_probe", {}).get("status") == "completed"
            for candidate in candidates
        ),
        "subtarget_conductors": sum(
            candidate.get("conductor_probe", {}).get(
                "below_strict_log_conductor_target"
            )
            is True
            for candidate in candidates
        ),
        "rank_triage_count": len(rank_records),
        "maximum_stable_numerical_rank": max(
            (
                record["full_pool_stable_numerical_rank"]
                for record in rank_records
            ),
            default=None,
        ),
        "exact_slice_result_sha256": slice_result_digest(
            artifact["direction_searches"]
        ),
        "wall_seconds": time.monotonic() - started,
    }
    artifact["target"]["hit"] = any(
        candidate.get("rank_triage", {})
        .get("finite_reduction_attempt", {})
        .get("certified_algebraic_rank_lower_bound", 0)
        >= 21
        and candidate.get("conductor_probe", {}).get(
            "below_strict_log_conductor_target"
        )
        for candidate in candidates
    )
    if not artifact["target"]["hit"]:
        artifact["target"]["reason"] = (
            "no group-manufactured subtarget fibre received an exact rank-21 certificate"
        )
    artifact["status"] = (
        "completed exact rank22 record-group orbit and bounded H=50000 manufactured-slice search"
    )
    artifact["execution"].update(
        {
            "phase": "complete",
            "wall_seconds": time.monotonic() - started,
            "owned_processes_remaining": 0,
        }
    )
    artifact["completed_at_utc"] = datetime.now(timezone.utc).isoformat()
    write_artifact(args.output, artifact)


if __name__ == "__main__":
    main()
