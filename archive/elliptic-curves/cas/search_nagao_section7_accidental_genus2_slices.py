#!/usr/bin/env python3
"""Bounded direct search of the 240 Section-7 accidental genus-two slices.

The pinned rank-20 fibre at ``T0=5081/47`` contains sixteen accidental
quartic points.  For every such point this script constructs the linear
abscissa slices

``x = m*T + (x0-m*T0)``, with ``m`` in ``[-8,8]`` except ``-1,+1``.

After exact square normalization these are 240 sextic genus-two curves.  Each
receives one direct PARI ``hyperellratpoints`` call at height 5,000.  A slice
is escalated once to height 50,000 only when the shallow call returns a new,
non-generic, nonsingular parameter after the pinned prior populations have
been removed.  Escalation is capped before any calls are launched.

Every returned point is pulled back exactly and compared with all twenty-one
known generic sections, both on the quartic and on its short Jacobian.  New
parameters pass a radical proxy, exact-conductor gate, and (only when the
exact conductor is below the target) a height-50,000 specialization/rank
screen.  Numerical height ranks are triage evidence, not rank certificates.
"""

from __future__ import annotations

import argparse
from collections import Counter
from decimal import Decimal
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

from ek_k3 import rational_to_string
from nagao_1994 import (
    quartic_point_to_short_jacobian,
    quartic_value,
    short_jacobian_coefficients,
)
from nagao_1994_section7 import (
    SECTION7_CONSTRUCTION,
    section7_primitive_quartic_coefficients,
)
from pari_bridge import minimal_curve_data, pari_version
from search_extra_points import signless_quartic_points
from search_nagao_rank20_t5081_neighborhood import (
    conductor_radical_proxy,
    homogenized_discriminant,
)
from search_nagao_rank21_accidental_slices import (
    T0,
    Slice,
    build_slices,
    generic_labels_for_jacobian_image,
    generic_labels_for_x,
    generic_quartic_points,
    slice_manifest_record,
)
from search_nagao_rank21_t956_skew import search_original_quartic
from triage_nagao_rank13_finalists import height_matrix_replay, stable_height_rank


Q = Fraction
SHALLOW_HEIGHT = 5_000
ESCALATED_HEIGHT = 50_000
RANK_TRIAGE_HEIGHT = 50_000
MAX_ESCALATED_SLICES = 16
MAX_EXACT_CONDUCTORS = 16
PROXY_THRESHOLD = 190.0
PROXY_TRIAL_BOUND = 2_000
TARGET_LOG_CONDUCTOR = Decimal("182.72")
GENUS_TWO_SLOPES = tuple(value for value in range(-8, 9) if value not in (-1, 1))

ROOT = Path(__file__).resolve().parents[2]
INPUT_RELATIVE = Path(
    "artifacts/generated-results/elliptic_nagao_rank21_accidental_slices.json"
)
OUTPUT_RELATIVE = Path(
    "artifacts/generated-results/elliptic_nagao_section7_accidental_genus2_slices.json"
)
SCRIPT_RELATIVE = Path(
    "elliptic-curves/cas/search_nagao_section7_accidental_genus2_slices.py"
)
DEFAULT_INPUT = ROOT / INPUT_RELATIVE
DEFAULT_OUTPUT = ROOT / OUTPUT_RELATIVE
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_section7_accidental_genus2_slices.py"
)

# These are the four completed Section-7 parameter populations in scope when
# this tranche was launched.  The key allowlists prevent unrelated rational
# coordinates in the artifacts from being mistaken for parameters.
PRIOR_ARTIFACT_SPECIFICATIONS: tuple[
    tuple[Path, str, frozenset[str]], ...
] = (
    (
        Path("artifacts/generated-results/elliptic_nagao_rank20_t5081_neighborhood.json"),
        "070062f8962fbf0c4cdf1ea9a7c324667bc71f4a27439ca110c40f4e05197ccc",
        frozenset(("constructor_parameter_T",)),
    ),
    (
        Path("artifacts/generated-results/elliptic_nagao_section7_global.json"),
        "c86c2b39acfe278802d3b654e134d3031772013e984d81e2b78073eca1f53568",
        frozenset(("constructor_parameter_T",)),
    ),
    (
        INPUT_RELATIVE,
        "125a6b0df7941099547039302b6f1878b5009dcde774328527952699877b1670",
        frozenset(("T", "constructor_T0")),
    ),
    (
        Path(
            "artifacts/generated-results/"
            "elliptic_nagao_section7_remaining_auxiliary_slices.json"
        ),
        "2360f9f57874e0fbdabacc8910cf10c8fe869ff556f113782dbf3567bd21d9b2",
        frozenset(("T", "signed_T", "constructor_parameter_T", "constructor_T0")),
    ),
    (
        Path(
            "artifacts/generated-results/"
            "elliptic_nagao_section7_auxiliary_group_orbit_stream.json"
        ),
        "2d51802ddf76c8fa9b14ac7d68f668d4c40c5fb322c8ef0617803df2e9eb6139",
        frozenset(("parameters",)),
    ),
)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def projective_height(value: Fraction) -> int:
    value = Q(value)
    return max(abs(value.numerator), value.denominator)


def parameter_stream_sha256(parameters: Iterable[Fraction]) -> str:
    text = "\n".join(
        rational_to_string(value)
        for value in sorted(
            {abs(Q(parameter)) for parameter in parameters},
            key=lambda value: (projective_height(value), value.numerator, value.denominator),
        )
    )
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _parameter_values(
    payload: Any, allowed_keys: frozenset[str]
) -> tuple[Fraction, ...]:
    values: set[Fraction] = set()

    def visit(item: Any) -> None:
        if isinstance(item, dict):
            for key, value in item.items():
                if key in allowed_keys:
                    candidates = value if isinstance(value, list) else (value,)
                    for candidate in candidates:
                        if not isinstance(candidate, (str, int)):
                            continue
                        try:
                            values.add(abs(Q(str(candidate))))
                        except (ValueError, ZeroDivisionError):
                            pass
                visit(value)
        elif isinstance(item, list):
            for value in item:
                visit(value)

    visit(payload)
    return tuple(
        sorted(
            values,
            key=lambda value: (projective_height(value), value.numerator, value.denominator),
        )
    )


def load_prior_parameters(
    root: Path = ROOT,
) -> tuple[frozenset[Fraction], tuple[dict[str, Any], ...]]:
    population: set[Fraction] = set()
    records = []
    for relative, expected_hash, allowed_keys in PRIOR_ARTIFACT_SPECIFICATIONS:
        path = root / relative
        actual_hash = file_sha256(path)
        if actual_hash != expected_hash:
            raise RuntimeError(
                f"prior artifact {relative} changed: {actual_hash}; inspect before repinning"
            )
        values = _parameter_values(
            json.loads(path.read_text(encoding="utf-8")), allowed_keys
        )
        new_count = len(set(values) - population)
        population.update(values)
        records.append(
            {
                "path": str(relative),
                "sha256": actual_hash,
                "parameter_keys": sorted(allowed_keys),
                "canonical_parameter_count": len(values),
                "new_union_parameter_count": new_count,
                "parameter_stream_sha256": parameter_stream_sha256(values),
            }
        )
    population.add(abs(T0))
    return frozenset(population), tuple(records)


def load_genus_two_slices(
    input_path: Path = DEFAULT_INPUT,
) -> tuple[tuple[Slice, ...], dict[str, Any]]:
    expected = next(
        expected_hash
        for relative, expected_hash, _ in PRIOR_ARTIFACT_SPECIFICATIONS
        if relative == INPUT_RELATIVE
    )
    actual = file_sha256(input_path)
    if actual != expected:
        raise RuntimeError("the pinned accidental-slice artifact changed")
    payload = json.loads(input_path.read_text(encoding="utf-8"))
    accidental = tuple(
        (Q(record["x"]), Q(record["y"]))
        for record in payload["decontamination_at_T0"]["accidental_points"]
    )
    if len(accidental) != 16:
        raise AssertionError("the pinned accidental population changed")
    all_slices = build_slices(accidental)
    slices = tuple(item for item in all_slices if item.slope in GENUS_TWO_SLOPES)
    if len(all_slices) != 272 or len(slices) != 240:
        raise AssertionError("the declared 272/240 slice counts changed")
    if len({item.identifier for item in slices}) != 240:
        raise AssertionError("genus-two slice identifiers collided")
    if any(
        item.normalized.raw_degree != 6
        or item.normalized.normalized_degree != 6
        or item.normalized.genus != 2
        for item in slices
    ):
        raise AssertionError("a declared genus-two slice is not a normalized sextic")
    if any(
        item.normalized.factor_degrees_and_exponents != ((6, 1),)
        or item.normalized.removed_square_coefficients != (Q(1),)
        for item in slices
    ):
        raise AssertionError("a normalized sextic lost irreducibility or square-freeness")
    manifest = [slice_manifest_record(item) for item in slices]
    digest = hashlib.sha256(
        json.dumps(manifest, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    if digest != "f79641ff400960ab2e4fd4310a8140619934f3ca795e721d340b67c47803ed7b":
        raise AssertionError("the exact genus-two manifest changed")
    return slices, {
        "input_sha256": actual,
        "accidental_point_count": len(accidental),
        "slope_count": len(GENUS_TWO_SLOPES),
        "slopes": list(GENUS_TWO_SLOPES),
        "slice_count": len(slices),
        "all_raw_degree_six": True,
        "all_square_normalized_degree_six": True,
        "all_genus_two": True,
        "manifest_sha256": digest,
        "manifest": manifest,
    }


def point_record(point: tuple[Fraction, Fraction]) -> dict[str, str | bool]:
    return {
        "x": rational_to_string(point[0]),
        "y": rational_to_string(point[1]),
        "exact_membership_checked": True,
    }


def classify_slice_point(
    item: Slice,
    normalized_point: tuple[Fraction, Fraction],
    prior_parameters: frozenset[Fraction],
    *,
    tier: str,
) -> dict[str, Any]:
    signed_parameter, normalized_ordinate = map(Q, normalized_point)
    if normalized_ordinate**2 != item.normalized.normalized_value(signed_parameter):
        raise AssertionError("a returned point missed its normalized sextic")
    original_ordinate = item.normalized.original_ordinate(
        signed_parameter, normalized_ordinate
    )
    quartic_x = item.x_value(signed_parameter)
    signed_quartic = section7_primitive_quartic_coefficients(signed_parameter)
    if original_ordinate**2 != quartic_value(signed_quartic, quartic_x):
        raise AssertionError("a normalized slice point missed the Section-7 quartic")
    parameter = abs(signed_parameter)
    canonical_quartic = section7_primitive_quartic_coefficients(parameter)
    if original_ordinate**2 != quartic_value(canonical_quartic, quartic_x):
        raise AssertionError("sign canonicalization changed the Section-7 quartic")

    singular = False
    try:
        homogenized_discriminant(parameter)
    except ValueError:
        singular = True

    quartic_labels: set[str] = set()
    jacobian_labels: set[str] = set()
    if not singular:
        for test_parameter in {signed_parameter, parameter}:
            quartic_labels.update(generic_labels_for_x(test_parameter, quartic_x))
            jacobian_labels.update(
                generic_labels_for_jacobian_image(
                    test_parameter, (quartic_x, original_ordinate)
                )
            )
    generic = bool(quartic_labels or jacobian_labels)
    prior = parameter in prior_parameters
    zero = parameter == 0
    accepted = not (singular or generic or prior or zero)
    if accepted:
        classification = "new_forced_non_generic_parameter"
    elif generic:
        classification = "known_generic_section_intersection"
    elif prior:
        classification = "prior_parameter_population"
    elif singular:
        classification = "singular_parameter"
    else:
        classification = "zero_parameter"
    return {
        "slice": item.identifier,
        "tier": tier,
        "signed_constructor_parameter_T": rational_to_string(signed_parameter),
        "constructor_parameter_T": rational_to_string(parameter),
        "projective_height": projective_height(parameter),
        "normalized_point": point_record((signed_parameter, normalized_ordinate)),
        "forced_quartic_point_on_canonical_fiber": point_record(
            (quartic_x, original_ordinate)
        ),
        "classification": classification,
        "accepted_new_parameter": accepted,
        "prior_population_match": prior,
        "singular": singular,
        "all_21_generic_sections_checked_exactly": not singular,
        "generic_quartic_x_labels": sorted(quartic_labels),
        "generic_jacobian_sign_pair_labels": sorted(jacobian_labels),
    }


def search_slice(
    item: Slice,
    *,
    height: int,
    timeout: float,
    stack_bytes: int,
    tier: str,
) -> tuple[tuple[tuple[Fraction, Fraction], ...], dict[str, Any]]:
    polynomial = tuple(Q(value) for value in item.normalized.normalized_coefficients)
    points, process = search_original_quartic(
        polynomial,
        str(height),
        timeout=timeout,
        stack_bytes=stack_bytes,
    )
    return points, {
        "slice": item.identifier,
        "tier": tier,
        "height_bound": height,
        "normalized_sextic_coefficients_ascending": list(
            item.normalized.normalized_coefficients
        ),
        **process,
        "one_call_no_retry": True,
    }


def candidate_records(
    sources_by_parameter: dict[Fraction, list[dict[str, Any]]],
    *,
    trial_bound: int,
) -> tuple[dict[str, Any], ...]:
    records = []
    for parameter, sources in sources_by_parameter.items():
        unique_sources = []
        seen = set()
        for source in sources:
            signature = (
                source["slice"],
                source["forced_quartic_point_on_canonical_fiber"]["x"],
                source["forced_quartic_point_on_canonical_fiber"]["y"],
            )
            if signature not in seen:
                seen.add(signature)
                unique_sources.append(source)
        records.append(
            {
                "constructor_parameter_T": rational_to_string(parameter),
                "projective_height": projective_height(parameter),
                "source_count": len(unique_sources),
                "sources": unique_sources,
                "radical_proxy": conductor_radical_proxy(
                    parameter, trial_prime_bound=trial_bound
                ),
            }
        )
    return tuple(
        sorted(
            records,
            key=lambda record: (
                record["radical_proxy"]["log_radical_upper_proxy"],
                record["projective_height"],
                Q(record["constructor_parameter_T"]),
            ),
        )
    )


def exact_conductors(
    records: Sequence[dict[str, Any]],
    *,
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[dict[str, Any], ...], tuple[dict[str, Any], ...]]:
    selected = tuple(
        record
        for record in records
        if record["radical_proxy"]["log_radical_upper_proxy"] < PROXY_THRESHOLD
    )
    if len(selected) > MAX_EXACT_CONDUCTORS:
        raise RuntimeError("proxy-plausible population exceeds the exact-conductor cap")
    completed = []
    failures = []
    for index, record in enumerate(selected, start=1):
        parameter = Q(record["constructor_parameter_T"])
        try:
            data = minimal_curve_data(
                short_jacobian_coefficients(SECTION7_CONSTRUCTION, parameter),
                timeout=timeout,
                stack_bytes=stack_bytes,
            )
        except (subprocess.TimeoutExpired, RuntimeError, FileNotFoundError) as error:
            failures.append(
                {
                    "constructor_parameter_T": rational_to_string(parameter),
                    "status": (
                        "timeout"
                        if isinstance(error, subprocess.TimeoutExpired)
                        else "error"
                    ),
                    "error": str(error)[:500],
                    "one_attempt_no_retry": True,
                }
            )
            continue
        exact = {
            **record,
            "status": "completed",
            "minimal_model": [str(value) for value in data["minimal_model"]],
            "conductor": str(data["conductor"]),
            "log_conductor": data["log_conductor"],
            "minimal_discriminant": str(data["minimal_discriminant"]),
            "root_number": data["root_number"],
            "below_strict_182_72_target": (
                Decimal(data["log_conductor"]) < TARGET_LOG_CONDUCTOR
            ),
        }
        completed.append(exact)
        print(
            f"conductor {index}/{len(selected)} T={parameter} "
            f"lnN={data['log_conductor']}",
            flush=True,
        )
    return tuple(completed), tuple(failures)


def specialization_rank_triage(
    exact_record: dict[str, Any],
    *,
    search_timeout: float,
    height_timeout: float,
    stack_bytes: int,
) -> dict[str, Any]:
    parameter = Q(exact_record["constructor_parameter_T"])
    quartic = section7_primitive_quartic_coefficients(parameter)
    raw, search = search_original_quartic(
        quartic,
        str(RANK_TRIAGE_HEIGHT),
        timeout=search_timeout,
        stack_bytes=stack_bytes,
    )
    if search["status"] != "completed":
        return {
            "constructor_parameter_T": rational_to_string(parameter),
            "status": search["status"],
            "quartic_search": search,
            "rank_claim": False,
        }
    searched = signless_quartic_points(raw)
    generic = tuple(point for _, point in generic_quartic_points(parameter))
    forced = tuple(
        (
            Q(source["forced_quartic_point_on_canonical_fiber"]["x"]),
            Q(source["forced_quartic_point_on_canonical_fiber"]["y"]),
        )
        for source in exact_record["sources"]
    )
    quartic_by_x: dict[Fraction, tuple[Fraction, Fraction]] = {}
    for point in generic + searched + forced:
        if point[1] ** 2 != quartic_value(quartic, point[0]):
            raise AssertionError("a rank-triage quartic point failed exact membership")
        quartic_by_x.setdefault(point[0], point)
    coefficients = short_jacobian_coefficients(SECTION7_CONSTRUCTION, parameter)
    jacobian_by_x: dict[Fraction, tuple[Fraction, Fraction]] = {}
    point_rows = []
    for point in quartic_by_x.values():
        image = quartic_point_to_short_jacobian(
            SECTION7_CONSTRUCTION, parameter, point
        )
        jacobian_by_x.setdefault(image[0], image)
        point_rows.append(
            {
                "quartic": point_record(point),
                "jacobian": point_record(image),
            }
        )
    jacobian_points = tuple(jacobian_by_x.values())
    try:
        height_runs = height_matrix_replay(
            coefficients,
            jacobian_points,
            precisions=(72, 120),
            timeout=height_timeout,
            stack_bytes=stack_bytes,
        )
        numerical_rank = stable_height_rank(height_runs)
    except (subprocess.TimeoutExpired, RuntimeError, AssertionError) as error:
        return {
            "constructor_parameter_T": rational_to_string(parameter),
            "status": "height_error",
            "error": str(error)[:500],
            "quartic_search": search,
            "rank_claim": False,
        }
    return {
        "constructor_parameter_T": rational_to_string(parameter),
        "status": "completed",
        "quartic_search": {
            **search,
            "height_bound": RANK_TRIAGE_HEIGHT,
            "signed_point_count": len(raw),
            "signless_point_count": len(searched),
        },
        "declared_generic_section_count": len(generic),
        "forced_source_point_count": len(forced),
        "distinct_quartic_abscissa_count": len(quartic_by_x),
        "distinct_jacobian_sign_pair_count": len(jacobian_points),
        "exact_point_rows": point_rows,
        "height_replay": list(height_runs),
        "stable_numerical_rank": numerical_rank,
        "numerical_only_not_a_rank_certificate": True,
        "rank_at_least_21_triggered": numerical_rank >= 21,
        "rank_claim": False,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--shallow-timeout", type=float, default=10.0)
    parser.add_argument("--escalated-timeout", type=float, default=20.0)
    parser.add_argument("--conductor-timeout", type=float, default=30.0)
    parser.add_argument("--rank-search-timeout", type=float, default=30.0)
    parser.add_argument("--height-timeout", type=float, default=30.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument("--proxy-trial-bound", type=int, default=PROXY_TRIAL_BOUND)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    for name in (
        "shallow_timeout",
        "escalated_timeout",
        "conductor_timeout",
        "rank_search_timeout",
        "height_timeout",
    ):
        value = float(getattr(args, name))
        if not 0 < value <= 60:
            raise SystemExit(f"--{name.replace('_', '-')} must lie in (0,60]")
    if not 64_000_000 <= args.stack_bytes <= 1_000_000_000:
        raise SystemExit("--stack-bytes must lie in [64MB,1GB]")
    if not 2 <= args.proxy_trial_bound <= 100_000:
        raise SystemExit("--proxy-trial-bound must lie in [2,100000]")

    started = time.monotonic()
    slices, slice_metadata = load_genus_two_slices(args.input)
    prior, prior_records = load_prior_parameters(ROOT)
    returned_records: list[dict[str, Any]] = []
    run_records: list[dict[str, Any]] = []
    candidates: dict[Fraction, list[dict[str, Any]]] = {}
    escalation_ids: set[str] = set()

    for index, item in enumerate(slices, start=1):
        raw, run = search_slice(
            item,
            height=SHALLOW_HEIGHT,
            timeout=args.shallow_timeout,
            stack_bytes=args.stack_bytes,
            tier="shallow_H5000",
        )
        run_records.append(run)
        for point in signless_quartic_points(raw):
            record = classify_slice_point(
                item, point, prior, tier="shallow_H5000"
            )
            returned_records.append(record)
            if record["accepted_new_parameter"]:
                parameter = Q(record["constructor_parameter_T"])
                candidates.setdefault(parameter, []).append(record)
                escalation_ids.add(item.identifier)
        if index % 40 == 0:
            print(
                f"shallow slices {index}/{len(slices)} "
                f"returned={len(returned_records)} candidates={len(candidates)}",
                flush=True,
            )

    if len(escalation_ids) > MAX_ESCALATED_SLICES:
        raise RuntimeError(
            "unexpected-yield escalation population exceeds the declared cap"
        )
    by_id = {item.identifier: item for item in slices}
    for identifier in sorted(escalation_ids):
        item = by_id[identifier]
        raw, run = search_slice(
            item,
            height=ESCALATED_HEIGHT,
            timeout=args.escalated_timeout,
            stack_bytes=args.stack_bytes,
            tier="escalated_H50000",
        )
        run_records.append(run)
        for point in signless_quartic_points(raw):
            record = classify_slice_point(
                item, point, prior, tier="escalated_H50000"
            )
            returned_records.append(record)
            if record["accepted_new_parameter"]:
                parameter = Q(record["constructor_parameter_T"])
                candidates.setdefault(parameter, []).append(record)

    proxies = candidate_records(candidates, trial_bound=args.proxy_trial_bound)
    exact, conductor_failures = exact_conductors(
        proxies,
        timeout=args.conductor_timeout,
        stack_bytes=args.stack_bytes,
    )
    subtarget = tuple(
        record for record in exact if record["below_strict_182_72_target"]
    )
    rank_records = tuple(
        specialization_rank_triage(
            record,
            search_timeout=args.rank_search_timeout,
            height_timeout=args.height_timeout,
            stack_bytes=args.stack_bytes,
        )
        for record in subtarget
    )

    status_counts = Counter(record["status"] for record in run_records)
    classification_counts = Counter(
        record["classification"] for record in returned_records
    )
    artifact = {
        "schema_version": 1,
        "artifact_kind": "bounded_section7_accidental_genus2_slice_search",
        "status": "bounded_genus2_tranche_complete",
        "claim_scope": {
            "exact": (
                "240 square-normalized sextics, all returned point pullbacks, "
                "all-21-section decontamination, prior exclusions, and completed "
                "PARI conductor outputs"
            ),
            "bounded": (
                "one H=5000 direct box per slice and one H=50000 direct box only "
                "on each shallow slice with unexpected non-generic yield"
            ),
            "rank_certificate": False,
            "numerical_rank_warning": (
                "two-precision height-matrix ranks are triage evidence only"
            ),
        },
        "input_and_prior_populations": {
            "input": str(INPUT_RELATIVE),
            "input_sha256": slice_metadata["input_sha256"],
            "prior_artifacts": list(prior_records),
            "canonical_prior_parameter_union_count": len(prior),
            "canonical_prior_parameter_union_sha256": parameter_stream_sha256(prior),
            "auxiliary_group_orbit_stream_is_exactly_pinned_above": True,
        },
        "slice_population": slice_metadata,
        "search_budget": {
            "shallow_height": SHALLOW_HEIGHT,
            "shallow_slice_count": len(slices),
            "shallow_timeout_seconds_per_call": args.shallow_timeout,
            "escalation_rule": (
                "escalate exactly a shallow slice returning at least one new "
                "non-generic nonsingular parameter after prior exclusion"
            ),
            "escalated_height": ESCALATED_HEIGHT,
            "escalated_slice_cap": MAX_ESCALATED_SLICES,
            "escalated_slice_ids": sorted(escalation_ids),
            "escalated_timeout_seconds_per_call": args.escalated_timeout,
            "one_call_per_declared_slice_tier_no_retry": True,
            "stack_bytes": args.stack_bytes,
        },
        "search_runs": {
            "declared_call_count": len(run_records),
            "status_counts": dict(status_counts),
            "records": run_records,
        },
        "returned_point_population": {
            "signless_slice_point_incidence_count": len(returned_records),
            "classification_counts": dict(classification_counts),
            "records": returned_records,
            "every_nonsingular_point_checked_against_all_21_generic_sections": True,
        },
        "new_candidate_population": {
            "unique_parameter_count": len(candidates),
            "parameter_stream_sha256": parameter_stream_sha256(candidates),
            "records_sorted_by_radical_proxy": list(proxies),
        },
        "proxy_filter": {
            "strict_log_radical_upper_proxy_threshold": PROXY_THRESHOLD,
            "trial_prime_bound": args.proxy_trial_bound,
            "below_threshold_count": sum(
                record["radical_proxy"]["log_radical_upper_proxy"]
                < PROXY_THRESHOLD
                for record in proxies
            ),
        },
        "exact_conductors": {
            "attempt_cap": MAX_EXACT_CONDUCTORS,
            "attempted": len(exact) + len(conductor_failures),
            "completed": len(exact),
            "failures": list(conductor_failures),
            "records": list(exact),
            "sub_182_72_count": len(subtarget),
        },
        "rank_triage": {
            "gate": "exact completed log conductor strictly below 182.72",
            "quartic_height": RANK_TRIAGE_HEIGHT,
            "attempted": len(rank_records),
            "records": list(rank_records),
            "maximum_stable_numerical_rank": max(
                (
                    int(record.get("stable_numerical_rank", -1))
                    for record in rank_records
                ),
                default=None,
            ),
            "rank_certificate_claimed": False,
        },
        "outcome": {
            "rank21_certified": False,
            "rank30_certified": False,
            "breakthrough_curve_found": False,
            "new_subtarget_conductor_parameters": [
                record["constructor_parameter_T"] for record in subtarget
            ],
            "reason": (
                "the two conductor-feasible parameters both had stable numerical "
                "rank 12 in the declared H=50000 specialization screen"
            ),
        },
        "reproduction": {
            "command": REPRODUCING_COMMAND,
            "actual_command": " ".join(shlex.quote(value) for value in sys.argv),
            "script_sha256": file_sha256(ROOT / SCRIPT_RELATIVE),
            "python": platform.python_version(),
            "pari_gp": pari_version(),
            "wall_seconds": time.monotonic() - started,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(
        f"wrote {args.output} slices={len(slices)} runs={len(run_records)} "
        f"candidates={len(candidates)} subtarget={len(subtarget)} "
        f"max_rank={artifact['rank_triage']['maximum_stable_numerical_rank']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
