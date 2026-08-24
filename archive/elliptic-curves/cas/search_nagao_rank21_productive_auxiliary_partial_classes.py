#!/usr/bin/env python3
"""Search one- and two-prime finite classes on the productive rank-21 slice.

The companion finite-class calculation forces the three discriminant-power
conditions at 13, 37, and 83 simultaneously.  Those shortest lifts are much
too tall.  This calculation weakens the local requirement in a controlled
way: for every attainable target point in one finite group it retains the 32
shortest full-product representatives, and for every attainable target pair
it retains the eight shortest representatives.  Each representative comes
from the complete 20*48*90 finite product reached with coefficients in
[-8,8]^8.

Every retained vector is lifted with exact rational group law.  The resulting
specialization and forced quartic abscissa are checked exactly, the published
fiber and generic visible sections are removed, and the remaining parameters
receive the same rigorous small-prime radical upper proxy as the parent
orbit.  The result is exhaustive only for this declared representative
policy, not for arbitrary kernel lifts of a finite class.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import sys
import time
from typing import Any, Iterable, Sequence

from ek_k3 import primes_up_to, rational_to_string
from search_nagao_rank21_productive_auxiliary_local_classes import (
    COEFFICIENT_BOUND,
    EXPECTED_INPUT_SCRIPT_SHA256,
    EXPECTED_INPUT_SHA256,
    FiniteCurve,
    LOCAL_PRIMES,
    TARGET_RESIDUES,
    load_input,
    parameter_residues,
    shortest_state_vectors,
)
from search_nagao_rank21_productive_auxiliary_orbit import (
    PRODUCTIVE_INTERCEPT,
    exact_radical_proxy,
    file_sha256,
    visible_abscissas,
)
from search_nagao_section7_auxiliary_jacobians import (
    weierstrass_add,
    weierstrass_multiply,
)


Q = Fraction
PARENT_ENGINE_SHA256 = (
    "25e6798a801d5f56700f3be68ff82c41a2968cc9829783c9fdaf42ef0a7d0136"
)
DEFAULT_OUTPUT = Path(
    "artifacts/generated-results/"
    "elliptic_nagao_rank21_productive_auxiliary_partial_classes.json"
)
SINGLE_REPRESENTATIVES_PER_TARGET_POINT = 32
PAIR_REPRESENTATIVES_PER_TARGET_PAIR = 8
PROXY_PRIME_BOUND = 1000
PROXY_GATE = 190.0


def vector_priority(vector: Sequence[int]) -> tuple[Any, ...]:
    return (
        sum(abs(value) for value in vector),
        max(abs(value) for value in vector),
        tuple(vector),
    )


def selected_vectors(
    states: dict[tuple[int, ...], tuple[int, ...]],
    residue_maps: Sequence[Sequence[int | None]],
) -> tuple[dict[tuple[int, ...], set[str]], dict[str, Any]]:
    selected: dict[tuple[int, ...], set[str]] = defaultdict(set)
    single_group_counts: dict[str, int] = {}
    pair_group_counts: dict[str, int] = {}

    for index in range(len(LOCAL_PRIMES)):
        groups: dict[int, list[tuple[Any, ...]]] = defaultdict(list)
        for state, vector in states.items():
            residue = residue_maps[index][state[index]]
            if residue is None or residue not in TARGET_RESIDUES[index]:
                continue
            groups[state[index]].append((*vector_priority(vector), vector))
        label = f"p{LOCAL_PRIMES[index]}"
        single_group_counts[label] = len(groups)
        for records in groups.values():
            records.sort()
            for *_, vector in records[:SINGLE_REPRESENTATIVES_PER_TARGET_POINT]:
                selected[tuple(vector)].add(label)

    for left, right in ((0, 1), (0, 2), (1, 2)):
        groups: dict[tuple[int, int], list[tuple[Any, ...]]] = defaultdict(list)
        for state, vector in states.items():
            residues = (
                residue_maps[left][state[left]],
                residue_maps[right][state[right]],
            )
            if any(residue is None for residue in residues):
                continue
            if (
                residues[0] not in TARGET_RESIDUES[left]
                or residues[1] not in TARGET_RESIDUES[right]
            ):
                continue
            groups[(state[left], state[right])].append(
                (*vector_priority(vector), vector)
            )
        label = f"p{LOCAL_PRIMES[left]}_p{LOCAL_PRIMES[right]}"
        pair_group_counts[label] = len(groups)
        for records in groups.values():
            records.sort()
            for *_, vector in records[:PAIR_REPRESENTATIVES_PER_TARGET_PAIR]:
                selected[tuple(vector)].add(label)

    return dict(selected), {
        "single_target_finite_point_counts": single_group_counts,
        "pair_target_finite_point_counts": pair_group_counts,
        "unique_selected_vector_count": len(selected),
    }


def lift_vector(
    auxiliary: Any,
    basis: Sequence[tuple[Fraction, Fraction]],
    vector: Sequence[int],
) -> tuple[Fraction, Fraction] | None:
    point = None
    for basis_point, scalar in zip(basis, vector, strict=True):
        point = weierstrass_add(
            auxiliary.weierstrass_coefficients,
            point,
            weierstrass_multiply(
                auxiliary.weierstrass_coefficients, basis_point, scalar
            ),
        )
    inverse = auxiliary.inverse(point)
    if inverse is None or inverse[0] == 0:
        return None
    signed_parameter = inverse[0]
    parameter = abs(signed_parameter)
    forced_x = -signed_parameter + PRODUCTIVE_INTERCEPT
    return parameter, forced_x


def digest_strings(values: Iterable[str]) -> str:
    return hashlib.sha256("\n".join(values).encode()).hexdigest()


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=root / DEFAULT_OUTPUT)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    root = Path(__file__).resolve().parents[2]
    started = time.monotonic()
    parent_engine = (
        root
        / "elliptic-curves"
        / "cas"
        / "search_nagao_rank21_productive_auxiliary_local_classes.py"
    )
    if file_sha256(parent_engine) != PARENT_ENGINE_SHA256:
        raise AssertionError("the finite-class engine changed")
    _, auxiliary, basis = load_input(root)
    curves = tuple(
        FiniteCurve(auxiliary.weierstrass_coefficients, prime)
        for prime in LOCAL_PRIMES
    )
    residue_maps = tuple(parameter_residues(curve, auxiliary) for curve in curves)
    states, depth_counts = shortest_state_vectors(curves, basis)
    if len(states) != 20 * 48 * 90:
        raise AssertionError("the complete finite product changed")
    vectors, selection = selected_vectors(states, residue_maps)

    incidences: dict[tuple[Fraction, Fraction], dict[str, Any]] = {}
    exceptional_vector_count = 0
    for vector, labels in sorted(vectors.items(), key=lambda item: vector_priority(item[0])):
        lifted = lift_vector(auxiliary, basis, vector)
        if lifted is None:
            exceptional_vector_count += 1
            continue
        parameter, forced_x = lifted
        key = parameter, forced_x
        record = incidences.setdefault(
            key,
            {
                "parameter": parameter,
                "forced_quartic_x": forced_x,
                "vectors": [],
                "selection_labels": set(),
            },
        )
        record["vectors"].append(tuple(vector))
        record["selection_labels"].update(labels)

    primes = primes_up_to(PROXY_PRIME_BOUND)
    parameter_records: dict[Fraction, dict[str, Any]] = {}
    published_parameter = auxiliary.base_parameter
    generic_incidence_count = 0
    published_incidence_count = 0
    for record in incidences.values():
        parameter = record["parameter"]
        forced_x = record["forced_quartic_x"]
        generic = forced_x in visible_abscissas(parameter)
        if generic:
            generic_incidence_count += 1
        if parameter == published_parameter:
            published_incidence_count += 1
        if generic or parameter == published_parameter:
            continue
        aggregate = parameter_records.setdefault(
            parameter,
            {
                "parameter": parameter,
                "incidences": [],
                "selection_labels": set(),
            },
        )
        aggregate["incidences"].append(record)
        aggregate["selection_labels"].update(record["selection_labels"])

    for record in parameter_records.values():
        record["proxy"] = exact_radical_proxy(record["parameter"], primes)
    ordered = sorted(
        parameter_records.values(),
        key=lambda record: (
            record["proxy"]["log_radical_upper_proxy"],
            max(abs(record["parameter"].numerator), record["parameter"].denominator),
            record["parameter"],
        ),
    )
    below_gate = [
        record
        for record in ordered
        if record["proxy"]["log_radical_upper_proxy"] < PROXY_GATE
    ]
    if below_gate:
        raise AssertionError("a partial-class lift crossed the conductor gate")

    def stored_record(record: dict[str, Any]) -> dict[str, Any]:
        return {
            "parameter": rational_to_string(record["parameter"]),
            "selection_labels": sorted(record["selection_labels"]),
            "forced_incidence_count": len(record["incidences"]),
            "forced_quartic_x": [
                rational_to_string(item["forced_quartic_x"])
                for item in record["incidences"][:4]
            ],
            **record["proxy"],
        }

    script_path = Path(__file__).resolve()
    artifact = {
        "schema_version": 1,
        "status": "complete_partial_finite_class_representative_screen",
        "target_hit": False,
        "source": {
            "productive_orbit_artifact_sha256": EXPECTED_INPUT_SHA256,
            "productive_orbit_script_sha256": EXPECTED_INPUT_SCRIPT_SHA256,
            "finite_class_engine_sha256": PARENT_ENGINE_SHA256,
            "auxiliary_basis_count": len(basis),
        },
        "finite_product": {
            "primes": list(LOCAL_PRIMES),
            "group_orders": [len(curve.points) for curve in curves],
            "target_residue_unions": [sorted(values) for values in TARGET_RESIDUES],
            "coefficient_range_per_coordinate": [-COEFFICIENT_BOUND, COEFFICIENT_BOUND],
            "state_counts_after_each_basis_coordinate": list(depth_counts),
            "complete_state_count": len(states),
        },
        "selection": {
            "single_representatives_per_target_finite_point": (
                SINGLE_REPRESENTATIVES_PER_TARGET_POINT
            ),
            "pair_representatives_per_target_finite_pair": (
                PAIR_REPRESENTATIVES_PER_TARGET_PAIR
            ),
            **selection,
            "selected_vector_sha256": digest_strings(
                ",".join(map(str, vector)) for vector in sorted(vectors)
            ),
        },
        "exact_lifts": {
            "exceptional_vector_count": exceptional_vector_count,
            "distinct_parameter_x_incidence_count": len(incidences),
            "generic_visible_incidence_count": generic_incidence_count,
            "published_fiber_incidence_count": published_incidence_count,
            "decontaminated_distinct_parameter_count": len(parameter_records),
            "parameter_sha256": digest_strings(
                rational_to_string(parameter) for parameter in sorted(parameter_records)
            ),
            "proxy_prime_bound": PROXY_PRIME_BOUND,
            "strict_proxy_gate": PROXY_GATE,
            "proxy_below_gate_count": len(below_gate),
            "minimum_log_radical_upper_proxy": (
                ordered[0]["proxy"]["log_radical_upper_proxy"] if ordered else None
            ),
            "top_32": [stored_record(record) for record in ordered[:32]],
        },
        "conclusion": {
            "target_hit": False,
            "exact_conductor_call_count": 0,
            "reason": (
                "After exact published-fiber and generic-section decontamination, "
                "no declared one- or two-prime representative has radical upper "
                "proxy below 190."
            ),
            "scope_warning": (
                "This closes only the stated top-32/top-8 representative policy; "
                "other kernel lifts of the same finite classes remain open."
            ),
        },
        "reproduction": {
            "command": (
                "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
                "elliptic-curves/cas/"
                "search_nagao_rank21_productive_auxiliary_partial_classes.py"
            ),
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "script_sha256": file_sha256(script_path),
            "wall_seconds": time.monotonic() - started,
            "no_subprocesses": True,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(
        f"states={len(states)} vectors={len(vectors)} "
        f"parameters={len(parameter_records)} "
        f"min_proxy={artifact['exact_lifts']['minimum_log_radical_upper_proxy']:.12f}",
        flush=True,
    )
    print(f"wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
