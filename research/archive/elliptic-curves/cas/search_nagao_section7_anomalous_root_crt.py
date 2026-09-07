#!/usr/bin/env python3
"""Search a new anomalous-root CRT lane in Nagao's Section-7 family.

The primitive discriminant has four unusually cheap local conditions: one
affine residue modulo each of 11, 19, 43, and 47 already forces valuation at
least three.  This search exhausts reduced-lattice shells for every two-prime
profile, plus a declared low-height shell extension.  It keeps only primitive positive
parameters with projective height at most 30,000 and denominator above 1,000,
so the population is geometrically disjoint from the completed 30,000 x 1,000
global box.

Every previously stored Section-7 auxiliary population is excluded exactly.
The earlier 1,574-member CRT neighborhood is pinned by its full replay digest
and by the exact 84-member subset outside the already excluded global box.
Population construction uses no trace data.  Training, held-forward, and
exact-final scores use three disjoint prime bands and omit the four CRT primes.
Exact conductor and root-number replay precedes H=50k/250k/1m point searches.
Numerical height ranks are triage only, and rank at least 21 is certified by
finite reductions before it can count as a target hit.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
import itertools
import json
from math import gcd, prod
from pathlib import Path
import platform
import shlex
import subprocess
import sys
import time
from typing import Any, Iterable, Sequence

from crt_lattice import crt_pair, gauss_reduce
from ek_k3 import rational_to_string
from multiple_root_lifting import affine_variable_coefficients, fixed_divisor_valuation
from nagao_1994 import PRIMARY_SOURCE, short_jacobian_coefficients
from nagao_1994_section7 import (
    SECTION7_CONSTRUCTION,
    SECTION7_CONSTRUCTOR_PARAMETER,
    SECTION7_ROOTS,
)
from pari_bridge import minimal_curve_data, pari_version
from search_nagao_rank20_t5081_neighborhood import (
    DISCRIMINANT_POLYNOMIAL,
    conductor_radical_proxy,
    homogenized_discriminant,
    integer_valuation,
    residue_score,
    residue_table,
)
from search_nagao_section7_accidental_genus2_slices import parameter_stream_sha256
from search_nagao_section7_global import (
    ConductorReplay,
    ExactScoreCandidate,
    PointPool,
    SCORE_SCALE,
    finite_reduction_certificate,
    parallel_point_search,
    parallel_rank_replay,
    point_priority,
)


Q = Fraction
ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
TARGET_LOG_CONDUCTOR = Decimal("182.72")
TARGET_RANK = 21
ALTERNATIVE_TARGET_RANK = 30

GLOBAL_A_MAX = 30_000
GLOBAL_B_MAX = 1_000
SHELL_RADIUS_BY_PROFILE_SIZE = {2: 96}
PROJECTIVE_HEIGHT_BY_PROFILE_SIZE = {2: 15_000}
PROJECTIVE_HEIGHT_MAX = max(PROJECTIVE_HEIGHT_BY_PROFILE_SIZE.values())
DENOMINATOR_MIN = 1_001
PROFILE_SIZES = (2,)
EXTENSION_RADIUS = 128
EXTENSION_HEIGHT_MAX = 5_000
EXTENSION_DENOMINATOR_MAX = 2_000

ANOMALOUS_ROOTS: dict[int, tuple[int, int]] = {
    11: (5, 6),
    19: (1, 18),
    43: (1, 42),
    47: (23, 24),
}
DESIGN_PRIMES = tuple(ANOMALOUS_ROOTS)
EXPECTED_FORCED_VALUATIONS = {11: 3, 19: 3, 43: 3, 47: 3}

DISCOVERY_PRIMES = (2_003, 2_011, 2_017, 2_027, 2_029, 2_039, 2_053)
HELD_FORWARD_PRIMES = (2_063, 2_069, 2_081, 2_083, 2_087, 2_089, 2_099, 2_111, 2_113)
PROXY_TRIAL_BOUND = 251
PROXY_GATE = 190.0
PROXY_DISCOVERY_KEEP = 4_096
DISCOVERY_FRONTIER_KEEP = 576
CONDUCTOR_KEEP = 64
POINT_KEEP = 24
STAGE_HEIGHTS = (50_000, 250_000, 1_000_000)
STAGE_KEEPS = (10, 3)
STAGE_TIMEOUTS = (8.0, 25.0, 90.0)
STAGE_WORKERS = (4, 4, 2)
HEIGHT_PRECISIONS = (72, 120)
CERTIFICATE_TRIGGER_RANK = 21

NEIGHBORHOOD_ARTIFACT = GENERATED / "elliptic_nagao_rank20_t5081_neighborhood.json"
GLOBAL_ARTIFACT = GENERATED / "elliptic_nagao_section7_global.json"
GENUS2_ARTIFACT = GENERATED / "elliptic_nagao_section7_accidental_genus2_slices.json"
A10_ARTIFACT = GENERATED / "elliptic_nagao_section7_a10_genus2_extension.json"
ACCIDENTAL_ARTIFACT = GENERATED / "elliptic_nagao_rank21_accidental_slices.json"
REMAINING_ARTIFACT = GENERATED / "elliptic_nagao_section7_remaining_auxiliary_slices.json"
ORBIT_STREAM_ARTIFACT = GENERATED / "elliptic_nagao_section7_auxiliary_group_orbit_stream.json"
EXPECTED_INPUT_SHA256 = {
    NEIGHBORHOOD_ARTIFACT.name: "070062f8962fbf0c4cdf1ea9a7c324667bc71f4a27439ca110c40f4e05197ccc",
    GLOBAL_ARTIFACT.name: "c86c2b39acfe278802d3b654e134d3031772013e984d81e2b78073eca1f53568",
    GENUS2_ARTIFACT.name: "112fb09a12aca5982311a64d161449887074eb13dad52b5b2cb752d93cf7c320",
    A10_ARTIFACT.name: "2715bce7c5a53ab8b6f7832dc13ac7836110763e6645fea849d79daa25bb9717",
    ACCIDENTAL_ARTIFACT.name: "125a6b0df7941099547039302b6f1878b5009dcde774328527952699877b1670",
    REMAINING_ARTIFACT.name: "2360f9f57874e0fbdabacc8910cf10c8fe869ff556f113782dbf3567bd21d9b2",
    ORBIT_STREAM_ARTIFACT.name: "2d51802ddf76c8fa9b14ac7d68f668d4c40c5fb322c8ef0617803df2e9eb6139",
}

EXPECTED_AUXILIARY_COUNT = 6_431
EXPECTED_AUXILIARY_SHA256 = "017d765c220846d4dc943d910f7bb0e401a1c574809531098c43459f61298e17"
EXPECTED_FINITE_OUTSIDE_BOX_COUNT = 6_488
EXPECTED_FINITE_OUTSIDE_BOX_SHA256 = "4ca379527b3564156f12de3eee45297ed1fd5840b94746267b35d997f413dda2"
EXPECTED_FULL_NEIGHBORHOOD_COUNT = 1_574
EXPECTED_FULL_NEIGHBORHOOD_SHA256 = "9919b047b7a5447d2cf33edfd7f55285f424704c299c407d39a00929903753f6"
EXPECTED_NEIGHBORHOOD_SOURCE_STREAM_SHA256 = "c69d201673f2e175f525e6fbd6420bc6f186c785df1f2cb4610fc0daf6d2b3c9"
EXPECTED_NEIGHBORHOOD_OUTSIDE_SHA256 = "5864d87b437ff465cdbd0f9ea98371f09d99e5fa2b6e48df83d5f4fb132861f7"

NEIGHBORHOOD_OUTSIDE_GLOBAL_BOX = tuple(
    Q(value)
    for value in """
30031/278 30037/278 30051/278 30059/278 30065/278 30071/278
30253/280 30281/280 30375/281 30401/281 30479/282 30499/282
30590/283 30683/284 30685/284 30703/284 30713/284 30895/286
30905/286 30913/286 30915/286 30923/286 30937/286 31175/288
31195/288 31228/289 31235/289 31240/289 31327/290 31337/290
31353/290 31373/290 31559/292 31565/292 31579/292 31585/292
31696/293 31906/295 31912/295 32017/296 32021/296 32197/298
32203/298 32209/298 32219/298 32221/298 32225/298 32235/298
32237/298 32629/302 32631/302 32639/302 32641/302 32645/302
32651/302 32655/302 32659/302 32665/302 32855/304 32857/304
33067/306 33095/306 33170/307 33190/307 33283/308 33289/308
33311/308 33503/310 33529/310 33829/313 33860/313 33925/314
33935/314 33943/314 33947/314 33953/314 33961/314 33965/314
34361/318 34508/319 34573/320 34601/320 34639/320 55399/512
""".split()
)

REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_section7_anomalous_root_crt.py"
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def projective_height(value: Fraction) -> int:
    value = abs(Q(value))
    return max(value.numerator, value.denominator)


def in_global_box(value: Fraction) -> bool:
    value = abs(Q(value))
    return (
        value != 0
        and value.numerator <= GLOBAL_A_MAX
        and value.denominator <= GLOBAL_B_MAX
    )


def profile_identifier(primes: Sequence[int]) -> str:
    return "root-" + "-".join(str(prime) for prime in primes)


def profile_primes() -> tuple[tuple[int, ...], ...]:
    return tuple(
        subset
        for size in PROFILE_SIZES
        for subset in itertools.combinations(DESIGN_PRIMES, size)
    )


def anomalous_root_audit() -> dict[str, Any]:
    records = {}
    for prime, residues in ANOMALOUS_ROOTS.items():
        values = []
        for residue in residues:
            valuation = fixed_divisor_valuation(
                affine_variable_coefficients(DISCRIMINANT_POLYNOMIAL, residue, prime),
                prime,
            )
            values.append({"residue": residue, "forced_valuation": valuation})
        if {record["forced_valuation"] for record in values} != {
            EXPECTED_FORCED_VALUATIONS[prime]
        }:
            raise AssertionError(f"the anomalous p={prime} valuation changed")
        if tuple(record["residue"] for record in values) != residues:
            raise AssertionError("the anomalous residue ordering changed")
        records[str(prime)] = values
    return records


def signless_crt_classes(primes: Sequence[int]) -> tuple[tuple[int, int], ...]:
    classes: set[tuple[int, int]] = set()
    for residues in itertools.product(*(ANOMALOUS_ROOTS[prime] for prime in primes)):
        residue, modulus = 0, 1
        for prime, local in zip(primes, residues):
            residue, modulus = crt_pair(residue, modulus, local, prime)
        residue %= modulus
        classes.add((min(residue, (-residue) % modulus), modulus))
    expected = 2 ** (len(primes) - 1)
    if len(classes) != expected:
        raise AssertionError("T -> -T did not halve the CRT classes")
    return tuple(sorted(classes))


def lattice_shell(
    residue: int,
    modulus: int,
    *,
    radius: int,
    height_max: int,
    denominator_max: int | None = None,
) -> tuple[Fraction, ...]:
    """Return the exact declared reduced-basis shell for one signless class."""

    basis = gauss_reduce((modulus, 0), (residue, 1))
    values = set()
    for left in range(-radius, radius + 1):
        for right in range(-radius, radius + 1):
            if left == 0 and right == 0:
                continue
            numerator = left * basis[0][0] + right * basis[1][0]
            denominator = left * basis[0][1] + right * basis[1][1]
            if denominator == 0:
                continue
            parameter = abs(Q(numerator, denominator))
            if parameter == 0 or gcd(parameter.denominator, modulus) != 1:
                continue
            if projective_height(parameter) > height_max:
                continue
            if parameter.denominator < DENOMINATOR_MIN:
                continue
            if denominator_max is not None and parameter.denominator > denominator_max:
                continue
            if (
                (parameter.numerator - residue * parameter.denominator) % modulus
                and (parameter.numerator + residue * parameter.denominator) % modulus
            ):
                # A lattice vector can have a common factor sharing the CRT
                # modulus; reducing that vector as a Fraction then weakens the
                # congruence.  Such a point is outside the declared class.
                continue
            values.add(parameter)
    return tuple(
        sorted(
            values,
            key=lambda value: (projective_height(value), value.denominator, value.numerator),
        )
    )


@dataclass(frozen=True)
class PopulationCandidate:
    parameter: Fraction
    profiles: tuple[str, ...]

    @property
    def identifier(self) -> str:
        return f"section7-anomalous-root-{self.parameter.numerator}-{self.parameter.denominator}"

    @property
    def height(self) -> int:
        return projective_height(self.parameter)


def generate_population() -> tuple[tuple[PopulationCandidate, ...], dict[str, Any]]:
    associations: dict[Fraction, set[str]] = {}
    profile_records = {}
    for primes in profile_primes():
        identifier = profile_identifier(primes)
        radius = SHELL_RADIUS_BY_PROFILE_SIZE[len(primes)]
        height_max = PROJECTIVE_HEIGHT_BY_PROFILE_SIZE[len(primes)]
        classes = signless_crt_classes(primes)
        profile_values = set()
        for residue, modulus in classes:
            profile_values.update(
                lattice_shell(
                    residue, modulus, radius=radius, height_max=height_max
                )
            )
        for parameter in profile_values:
            associations.setdefault(parameter, set()).add(identifier)
        profile_records[identifier] = {
            "primes": list(primes),
            "modulus": prod(primes),
            "signless_crt_classes": len(classes),
            "reduced_basis_coefficient_radius": radius,
            "projective_height_at_most": height_max,
            "shell_population_count": len(profile_values),
            "shell_population_sha256": parameter_stream_sha256(profile_values),
        }
        extension_identifier = f"{identifier}-low-height-extension"
        extension_values = set()
        for residue, modulus in classes:
            extension_values.update(
                lattice_shell(
                    residue,
                    modulus,
                    radius=EXTENSION_RADIUS,
                    height_max=EXTENSION_HEIGHT_MAX,
                    denominator_max=EXTENSION_DENOMINATOR_MAX,
                )
            )
        for parameter in extension_values:
            associations.setdefault(parameter, set()).add(extension_identifier)
        profile_records[extension_identifier] = {
            "primes": list(primes),
            "modulus": prod(primes),
            "signless_crt_classes": len(classes),
            "reduced_basis_coefficient_radius": EXTENSION_RADIUS,
            "projective_height_at_most": EXTENSION_HEIGHT_MAX,
            "denominator_interval": [DENOMINATOR_MIN, EXTENSION_DENOMINATOR_MAX],
            "shell_population_count": len(extension_values),
            "shell_population_sha256": parameter_stream_sha256(extension_values),
            "new_beyond_base_profile_count": len(extension_values - profile_values),
        }
    population = tuple(
        PopulationCandidate(parameter, tuple(sorted(profiles)))
        for parameter, profiles in sorted(
            associations.items(),
            key=lambda item: (projective_height(item[0]), item[0].denominator, item[0].numerator),
        )
    )
    if any(
        candidate.parameter.denominator <= GLOBAL_B_MAX
        or candidate.height > GLOBAL_A_MAX
        for candidate in population
    ):
        raise AssertionError("the anomalous-root shell overlaps the global rectangle")
    return population, {
        "profiles": profile_records,
        "exactly_deduplicated_count": len(population),
        "multi_profile_parameter_count": sum(len(item.profiles) > 1 for item in population),
        "population_sha256": parameter_stream_sha256(item.parameter for item in population),
    }


def extra_artifact_parameters(path: Path) -> tuple[Fraction, ...]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    return tuple(
        abs(Q(record["constructor_parameter_T"]))
        for record in payload["new_candidate_population"]["records_sorted_by_radical_proxy"]
    )


def exact_auxiliary_population() -> tuple[frozenset[Fraction], dict[str, Any]]:
    """Load exact structural parameter paths, never recursive key scraping."""

    accidental_payload = json.loads(ACCIDENTAL_ARTIFACT.read_text(encoding="utf-8"))
    accidental = {
        abs(Q(record["T"]))
        for record in accidental_payload["new_parameter_decontamination"]["records"]
    }
    remaining_payload = json.loads(REMAINING_ARTIFACT.read_text(encoding="utf-8"))
    remaining = {
        abs(Q(record["constructor_parameter_T"]))
        for record in remaining_payload["generation"]["all_candidate_records"]
    }
    orbit_payload = json.loads(ORBIT_STREAM_ARTIFACT.read_text(encoding="utf-8"))
    orbit = {abs(Q(value)) for value in orbit_payload["parameters"]}
    genus2 = set(extra_artifact_parameters(GENUS2_ARTIFACT))
    a10 = set(extra_artifact_parameters(A10_ARTIFACT))
    union = (
        accidental
        | remaining
        | orbit
        | genus2
        | a10
        | {abs(SECTION7_CONSTRUCTOR_PARAMETER)}
    )
    digest = parameter_stream_sha256(union)
    if len(union) != EXPECTED_AUXILIARY_COUNT or digest != EXPECTED_AUXILIARY_SHA256:
        raise AssertionError("the exact Section-7 auxiliary population changed")
    return frozenset(union), {
        "record_fiber_and_accidental_H200000_count": len(
            accidental | {abs(SECTION7_CONSTRUCTOR_PARAMETER)}
        ),
        "remaining_auxiliary_count": len(remaining),
        "group_orbit_count": len(orbit),
        "genus2_new_count": len(genus2),
        "a10_new_count": len(a10),
        "exact_union_count": len(union),
        "exact_union_sha256": digest,
        "inside_global_rectangle_count": sum(in_global_box(value) for value in union),
        "outside_global_rectangle_count": sum(not in_global_box(value) for value in union),
    }


def exact_exclusion_population() -> tuple[frozenset[Fraction], dict[str, Any]]:
    observed_hashes = {}
    for path in (
        NEIGHBORHOOD_ARTIFACT,
        GLOBAL_ARTIFACT,
        GENUS2_ARTIFACT,
        A10_ARTIFACT,
        ACCIDENTAL_ARTIFACT,
        REMAINING_ARTIFACT,
        ORBIT_STREAM_ARTIFACT,
    ):
        observed = sha256_file(path)
        expected = EXPECTED_INPUT_SHA256[path.name]
        if observed != expected:
            raise RuntimeError(f"pinned input changed: {path.name}: {observed}")
        observed_hashes[path.name] = observed

    neighborhood_outside = frozenset(NEIGHBORHOOD_OUTSIDE_GLOBAL_BOX)
    neighborhood_outside_digest = parameter_stream_sha256(neighborhood_outside)
    if (
        len(neighborhood_outside) != 84
        or neighborhood_outside_digest != EXPECTED_NEIGHBORHOOD_OUTSIDE_SHA256
    ):
        raise AssertionError("the pinned outside-box neighborhood list changed")
    auxiliary, auxiliary_audit = exact_auxiliary_population()
    finite_outside = {
        value
        for value in set(neighborhood_outside) | set(auxiliary)
        if not in_global_box(value)
    }
    finite_digest = parameter_stream_sha256(finite_outside)
    if (
        len(finite_outside) != EXPECTED_FINITE_OUTSIDE_BOX_COUNT
        or finite_digest != EXPECTED_FINITE_OUTSIDE_BOX_SHA256
    ):
        raise AssertionError("the exact outside-box exclusion population changed")
    return frozenset(finite_outside), {
        "geometric_global_rectangle_exclusion": {
            "positive_primitive_numerator_at_most": GLOBAL_A_MAX,
            "positive_primitive_denominator_at_most": GLOBAL_B_MAX,
            "proof_for_every_new_parameter": "denominator >= 1001 and height <= 30000",
            "global_artifact_sha256": observed_hashes[GLOBAL_ARTIFACT.name],
        },
        "full_neighborhood_deterministic_replay_checkpoint": {
            "full_count": EXPECTED_FULL_NEIGHBORHOOD_COUNT,
            "pure_parameter_sha256": EXPECTED_FULL_NEIGHBORHOOD_SHA256,
            "source_population_stream_sha256": EXPECTED_NEIGHBORHOOD_SOURCE_STREAM_SHA256,
            "outside_global_rectangle_count": len(neighborhood_outside),
            "outside_global_rectangle_parameters": [
                rational_to_string(value)
                for value in sorted(
                    neighborhood_outside,
                    key=lambda value: (
                        projective_height(value),
                        value.numerator,
                        value.denominator,
                    ),
                )
            ],
            "outside_global_rectangle_sha256": neighborhood_outside_digest,
            "replay_method": (
                "pinned trace/root beams regenerated to 1574, then exact box filter; "
                "the portable 84-member outside-box result is stored inline"
            ),
        },
        "exact_auxiliary_population": auxiliary_audit,
        "input_artifact_sha256": observed_hashes,
        "finite_outside_box_exclusion_count": len(finite_outside),
        "finite_outside_box_exclusion_sha256": finite_digest,
        "complete_prior_population_count_including_global_rectangle": 18_251_307,
        "complete_prior_population_sha256": "ef5e6dd6e564994ceecb32f8f3684f95c05781bd2fa2faaefac48f9311341d6b",
        "complete_prior_digest_is_a_pinned_independent_audit_not_reenumerated_here": True,
    }


@dataclass(frozen=True)
class LaneScanCandidate:
    parameter: Fraction
    training_scaled: int
    validation_scaled: int
    profiles: tuple[str, ...]

    @property
    def identifier(self) -> str:
        return f"section7-anomalous-root-{self.parameter.numerator}-{self.parameter.denominator}"

    @property
    def height(self) -> int:
        return projective_height(self.parameter)

    @property
    def training_score(self) -> float:
        return self.training_scaled / SCORE_SCALE

    @property
    def validation_score(self) -> float:
        return self.validation_scaled / SCORE_SCALE


@dataclass(frozen=True)
class LaneProxyCandidate:
    scanned: LaneScanCandidate
    proxy: dict[str, Any]

    @property
    def parameter(self) -> Fraction:
        return self.scanned.parameter

    @property
    def identifier(self) -> str:
        return self.scanned.identifier

    @property
    def height(self) -> int:
        return self.scanned.height

    @property
    def training_score(self) -> float:
        return self.scanned.training_score

    @property
    def validation_score(self) -> float:
        return self.scanned.validation_score

    @property
    def profiles(self) -> tuple[str, ...]:
        return self.scanned.profiles


def _add_quota(
    selected: dict[Fraction, LaneProxyCandidate],
    ordered: Iterable[LaneProxyCandidate],
    count: int,
) -> None:
    added = 0
    for candidate in ordered:
        if candidate.parameter in selected:
            continue
        selected[candidate.parameter] = candidate
        added += 1
        if added >= count:
            return


def proxy_discovery_population(
    population: Sequence[PopulationCandidate],
) -> tuple[tuple[LaneProxyCandidate, ...], dict[str, Any]]:
    scored = []
    singular = 0
    for item in population:
        try:
            proxy = conductor_radical_proxy(
                item.parameter, trial_prime_bound=PROXY_TRIAL_BOUND
            )
        except ValueError:
            singular += 1
            continue
        scanned = LaneScanCandidate(
            item.parameter,
            0,
            0,
            item.profiles,
        )
        scored.append(LaneProxyCandidate(scanned, proxy))

    by_proxy = sorted(
        scored,
        key=lambda item: (
            item.proxy["log_radical_upper_proxy"],
            -item.training_score,
            item.identifier,
        ),
    )
    selected: dict[Fraction, LaneProxyCandidate] = {}
    proxy_feasible = [
        item
        for item in by_proxy
        if item.proxy["log_radical_upper_proxy"] < PROXY_GATE
    ]
    _add_quota(selected, proxy_feasible, len(proxy_feasible))
    _add_quota(selected, by_proxy, PROXY_DISCOVERY_KEEP)
    for profile in sorted({profile for item in scored for profile in item.profiles}):
        _add_quota(
            selected,
            (item for item in by_proxy if profile in item.profiles),
            16,
        )
    _add_quota(selected, by_proxy, PROXY_DISCOVERY_KEEP - len(selected))
    answer = tuple(selected.values())[: min(PROXY_DISCOVERY_KEEP, len(selected))]
    return answer, {
        "nonsingular_scored_population": len(scored),
        "singular_exclusions": singular,
        "strict_proxy_below_190_count": len(proxy_feasible),
        "retained_for_exact_discovery_trace_scoring": len(answer),
        "retained_parameter_sha256": parameter_stream_sha256(item.parameter for item in answer),
    }


def attach_discovery_scores(
    candidates: Sequence[LaneProxyCandidate], discovery_tables: dict[int, Any]
) -> tuple[LaneProxyCandidate, ...]:
    answer = []
    for candidate in candidates:
        score, _, _ = residue_score(candidate.parameter, discovery_tables)
        scanned = LaneScanCandidate(
            candidate.parameter,
            round(score * SCORE_SCALE),
            0,
            candidate.profiles,
        )
        answer.append(LaneProxyCandidate(scanned, candidate.proxy))
    return tuple(answer)


def discovery_frontier(
    candidates: Sequence[LaneProxyCandidate],
) -> tuple[LaneProxyCandidate, ...]:
    by_discovery = sorted(
        candidates,
        key=lambda item: (
            -item.training_score,
            item.proxy["log_radical_upper_proxy"],
            item.identifier,
        ),
    )
    by_proxy = sorted(
        candidates,
        key=lambda item: (
            item.proxy["log_radical_upper_proxy"],
            -item.training_score,
            item.identifier,
        ),
    )
    selected: dict[Fraction, LaneProxyCandidate] = {}
    proxy_feasible = [
        item
        for item in by_proxy
        if item.proxy["log_radical_upper_proxy"] < PROXY_GATE
    ]
    _add_quota(selected, proxy_feasible, len(proxy_feasible))
    _add_quota(selected, by_proxy, 256)
    _add_quota(selected, by_discovery, 256)
    for profile in sorted({profile for item in candidates for profile in item.profiles}):
        _add_quota(
            selected,
            (item for item in by_discovery if profile in item.profiles),
            8,
        )
    _add_quota(selected, by_discovery, DISCOVERY_FRONTIER_KEEP - len(selected))
    return tuple(selected.values())[: min(DISCOVERY_FRONTIER_KEEP, len(selected))]


def attach_held_scores(
    candidates: Sequence[LaneProxyCandidate], held_tables: dict[int, Any]
) -> tuple[ExactScoreCandidate, ...]:
    answer = []
    for candidate in candidates:
        score, good, bad = residue_score(candidate.parameter, held_tables)
        scanned = LaneScanCandidate(
            candidate.parameter,
            candidate.scanned.training_scaled,
            round(score * SCORE_SCALE),
            candidate.profiles,
        )
        proxied = LaneProxyCandidate(scanned, candidate.proxy)
        answer.append(
            ExactScoreCandidate(
                proxied,
                str(score),
                good,
                bad,
                max(HELD_FORWARD_PRIMES),
            )
        )
    return tuple(
        sorted(
            answer,
            key=lambda item: (
                -item.validation_score,
                item.proxy["log_radical_upper_proxy"],
                item.identifier,
            ),
        )
    )


def select_conductors(
    candidates: Sequence[ExactScoreCandidate],
) -> tuple[ExactScoreCandidate, ...]:
    by_proxy = sorted(
        candidates,
        key=lambda item: (
            item.proxy["log_radical_upper_proxy"],
            -Decimal(item.score_b2000),
            item.identifier,
        ),
    )
    by_held = sorted(
        candidates,
        key=lambda item: (-item.validation_score, item.identifier),
    )
    by_discovery = sorted(
        candidates,
        key=lambda item: (-item.training_score, item.identifier),
    )
    selected: dict[Fraction, ExactScoreCandidate] = {}

    def add(ordered: Iterable[ExactScoreCandidate], count: int) -> None:
        added = 0
        for candidate in ordered:
            if candidate.parameter in selected:
                continue
            selected[candidate.parameter] = candidate
            added += 1
            if added >= count:
                return

    add(by_proxy, 24)
    add(by_discovery, 16)
    add(by_held, 24)
    for profile in sorted(
        {profile for item in candidates for profile in item.proxied.profiles}
    ):
        add((item for item in by_held if profile in item.proxied.profiles), 2)
    add(by_proxy, CONDUCTOR_KEEP - len(selected))
    return tuple(selected.values())[: min(CONDUCTOR_KEEP, len(selected))]


LOCAL_REPLAY_PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 43, 47)


def replay_one_conductor(
    candidate: ExactScoreCandidate, *, timeout: float, stack_bytes: int
) -> ConductorReplay:
    try:
        data = minimal_curve_data(
            short_jacobian_coefficients(SECTION7_CONSTRUCTION, candidate.parameter),
            timeout=timeout,
            local_primes=LOCAL_REPLAY_PRIMES,
            stack_bytes=stack_bytes,
        )
        data["below_strict_log_conductor_target"] = (
            Decimal(data["log_conductor"]) < TARGET_LOG_CONDUCTOR
        )
        return ConductorReplay(candidate, "completed", data)
    except (subprocess.TimeoutExpired, RuntimeError, AssertionError, ValueError) as error:
        return ConductorReplay(
            candidate,
            "timeout" if isinstance(error, subprocess.TimeoutExpired) else "error",
            {},
            str(error)[:500],
        )


def parallel_conductors(
    candidates: Sequence[ExactScoreCandidate],
    *,
    timeout: float,
    stack_bytes: int,
    workers: int,
) -> tuple[ConductorReplay, ...]:
    records = {}
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                replay_one_conductor,
                candidate,
                timeout=timeout,
                stack_bytes=stack_bytes,
            ): candidate.identifier
            for candidate in candidates
        }
        for future in as_completed(futures):
            records[futures[future]] = future.result()
    return tuple(records[candidate.identifier] for candidate in candidates)


def candidate_record(candidate: ExactScoreCandidate) -> dict[str, Any]:
    valuations = {
        str(prime): integer_valuation(
            homogenized_discriminant(candidate.parameter), prime
        )
        for prime in DESIGN_PRIMES
    }
    return {
        "candidate_id": candidate.identifier,
        "constructor_parameter_T": rational_to_string(candidate.parameter),
        "projective_height": candidate.height,
        "anomalous_root_profiles": list(candidate.proxied.profiles),
        "exact_homogenized_discriminant_valuations": valuations,
        "exact_discovery_trace_score": candidate.training_score,
        "exact_held_forward_trace_score": candidate.validation_score,
        "held_forward_good_primes": candidate.good_primes,
        "held_forward_bad_primes": candidate.bad_primes,
        "radical_proxy": candidate.proxy,
    }


def select_point_population(
    replays: Sequence[ConductorReplay], *, keep: int
) -> tuple[ExactScoreCandidate, ...]:
    completed = [item for item in replays if item.status == "completed"]
    root_minus = sorted(
        (
            item.candidate
            for item in completed
            if item.data.get("below_strict_log_conductor_target") is True
            and int(item.data["root_number"]) == -1
        ),
        key=lambda item: (-Decimal(item.score_b2000), item.identifier),
    )
    root_plus = sorted(
        (
            item.candidate
            for item in completed
            if item.data.get("below_strict_log_conductor_target") is True
            and int(item.data["root_number"]) == 1
        ),
        key=lambda item: (-Decimal(item.score_b2000), item.identifier),
    )
    by_exact = sorted(
        (item.candidate for item in completed),
        key=lambda item: (-Decimal(item.score_b2000), item.identifier),
    )
    selected: dict[Fraction, ExactScoreCandidate] = {}
    for ordered in (root_minus, root_plus, by_exact):
        for candidate in ordered:
            selected.setdefault(candidate.parameter, candidate)
            if len(selected) >= keep:
                return tuple(selected.values())
    return tuple(selected.values())


def lane_point_pool_record(pool: PointPool, rank: dict[str, Any]) -> dict[str, Any]:
    return {
        **candidate_record(pool.candidate),
        "point_search": {
            "height_bound": pool.height_bound,
            "status": pool.status,
            "signed_points": pool.signed_points,
            "distinct_quartic_abscissas": pool.signless_points,
            "predeclared_abscissas_returned": pool.predeclared_abscissas_returned,
            "new_distinct_jacobian_sign_pairs_beyond_21_predeclared": len(
                pool.new_images
            ),
            "wall_seconds": pool.wall_seconds,
            "pari_milliseconds": pool.pari_milliseconds,
            **({"error": pool.error} if pool.error else {}),
        },
        "height_rank": {
            key: value for key, value in rank.items() if key != "selected_points"
        },
    }


def run_point_stages(
    initial: Sequence[ExactScoreCandidate],
    conductor_by_id: dict[str, ConductorReplay],
    args: argparse.Namespace,
) -> tuple[list[dict[str, Any]], dict[str, tuple[PointPool, dict[str, Any]]]]:
    current = tuple(initial)
    stages = []
    best: dict[str, tuple[PointPool, dict[str, Any]]] = {}
    for index, (height, timeout, workers) in enumerate(
        zip(STAGE_HEIGHTS, STAGE_TIMEOUTS, STAGE_WORKERS), start=1
    ):
        pools = parallel_point_search(
            current,
            height_bound=height,
            timeout=timeout,
            stack_bytes=args.stack_bytes,
            workers=workers,
        )
        ranks = parallel_rank_replay(
            pools,
            precisions=HEIGHT_PRECISIONS,
            timeout=args.height_timeout,
            stack_bytes=args.stack_bytes,
            workers=workers,
        )
        ordered = tuple(
            sorted(
                pools,
                key=lambda pool: point_priority(
                    pool, ranks[pool.candidate.identifier], conductor_by_id
                ),
            )
        )
        for pool in ordered:
            rank = ranks[pool.candidate.identifier]
            if rank.get("status") != "completed":
                continue
            prior = best.get(pool.candidate.identifier)
            if prior is None or int(rank["stable_numerical_rank"]) >= int(
                prior[1]["stable_numerical_rank"]
            ):
                best[pool.candidate.identifier] = (pool, rank)
        keep = STAGE_KEEPS[index - 1] if index <= len(STAGE_KEEPS) else len(ordered)
        retained = ordered[: min(keep, len(ordered))]
        stages.append(
            {
                "stage": index,
                "quartic_naive_height_bound": height,
                "population_searched": len(pools),
                "completed": sum(pool.status == "completed" for pool in pools),
                "timeouts": sum(pool.status == "timeout" for pool in pools),
                "errors": sum(pool.status == "error" for pool in pools),
                "retained_constructor_parameters": [
                    rational_to_string(pool.candidate.parameter) for pool in retained
                ],
                "ranked_population": [
                    {
                        **lane_point_pool_record(pool, ranks[pool.candidate.identifier]),
                        "anomalous_root_profiles": list(pool.candidate.proxied.profiles),
                    }
                    for pool in ordered
                ],
            }
        )
        current = tuple(pool.candidate for pool in retained)
        maximum = max(
            (
                int(rank.get("stable_numerical_rank", -1))
                for rank in ranks.values()
                if rank.get("status") == "completed"
            ),
            default=-1,
        )
        print(f"point H={height} population={len(pools)} max_rank={maximum}", flush=True)
    return stages, best


def exact_checkpoints(
    best: dict[str, tuple[PointPool, dict[str, Any]]],
    conductor_by_id: dict[str, ConductorReplay],
    args: argparse.Namespace,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    checkpoints = []
    hits = []
    for identifier, (pool, rank) in sorted(
        best.items(), key=lambda item: (-int(item[1][1]["stable_numerical_rank"]), item[0])
    ):
        numerical_rank = int(rank["stable_numerical_rank"])
        if numerical_rank < CERTIFICATE_TRIGGER_RANK:
            continue
        try:
            certificate = finite_reduction_certificate(
                pool,
                rank,
                saturation_timeout=args.saturation_timeout,
                certificate_prime_bound=args.certificate_prime_bound,
                stack_bytes=args.stack_bytes,
            )
        except (subprocess.TimeoutExpired, RuntimeError, AssertionError, ValueError) as error:
            certificate = {
                "status": "timeout" if isinstance(error, subprocess.TimeoutExpired) else "error",
                "error": str(error)[:500],
            }
        replay = conductor_by_id[identifier]
        exact_rank = certificate.get("certified_algebraic_rank_lower_bound")
        target_hit = bool(
            certificate.get("status") == "certified"
            and exact_rank is not None
            and (
                int(exact_rank) >= ALTERNATIVE_TARGET_RANK
                or (
                    int(exact_rank) >= TARGET_RANK
                    and replay.status == "completed"
                    and replay.data.get("below_strict_log_conductor_target") is True
                )
            )
        )
        record = {
            "constructor_parameter_T": rational_to_string(pool.candidate.parameter),
            "trigger_stable_numerical_rank": numerical_rank,
            "trigger_height": pool.height_bound,
            "exact_rank_certificate": certificate,
            "conductor": replay.data if replay.status == "completed" else None,
            "target_hit": target_hit,
        }
        checkpoints.append(record)
        if target_hit:
            hits.append(record)
    return checkpoints, hits


def build_artifact(args: argparse.Namespace) -> dict[str, Any]:
    started = time.monotonic()
    local_audit = anomalous_root_audit()
    raw, raw_audit = generate_population()
    exclusions, exclusion_audit = exact_exclusion_population()
    novel = tuple(item for item in raw if item.parameter not in exclusions)
    overlaps = len(raw) - len(novel)
    if any(item.parameter in exclusions for item in novel):
        raise AssertionError("an exact prior parameter survived decontamination")

    discovery_tables = {prime: residue_table(prime) for prime in DISCOVERY_PRIMES}
    held_tables = {prime: residue_table(prime) for prime in HELD_FORWARD_PRIMES}
    if (
        set(discovery_tables) & set(held_tables)
        or set(discovery_tables) & set(DESIGN_PRIMES)
        or set(held_tables) & set(DESIGN_PRIMES)
    ):
        raise AssertionError("construction leaked into a scoring band")

    proxy_discovery, discovery_audit = proxy_discovery_population(novel)
    discovery_scored = attach_discovery_scores(proxy_discovery, discovery_tables)
    frontier = discovery_frontier(discovery_scored)
    exact = attach_held_scores(frontier, held_tables)
    conductor_candidates = select_conductors(exact)
    conductors = parallel_conductors(
        conductor_candidates,
        timeout=args.conductor_timeout,
        stack_bytes=args.stack_bytes,
        workers=args.workers,
    )
    conductor_by_id = {item.candidate.identifier: item for item in conductors}
    completed = [item for item in conductors if item.status == "completed"]
    subtarget = [
        item
        for item in completed
        if item.data.get("below_strict_log_conductor_target") is True
    ]
    print(
        f"raw={len(raw)} novel={len(novel)} held={len(exact)} "
        f"conductors={len(completed)}/{len(conductors)} subtarget={len(subtarget)}",
        flush=True,
    )

    point_candidates = select_point_population(
        conductors, keep=min(POINT_KEEP, len(completed))
    )
    point_stages, best = run_point_stages(point_candidates, conductor_by_id, args)
    checkpoints, hits = exact_checkpoints(best, conductor_by_id, args)
    maximum_rank = max(
        (int(rank["stable_numerical_rank"]) for _, rank in best.values()), default=-1
    )

    return {
        "schema_version": 1,
        "status": "bounded anomalous-root CRT Section-7 search; numerical ranks are triage only",
        "primary_source": PRIMARY_SOURCE,
        "family": {
            "roots_in_source_order": list(SECTION7_ROOTS),
            "parameter_symmetry": "T canonicalized to abs(T)",
            "predeclared_generic_sections": 21,
        },
        "target": {
            "rank_at_least": TARGET_RANK,
            "strict_log_conductor_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "alternative_rank_at_least": ALTERNATIVE_TARGET_RANK,
            "certified_hits": hits,
        },
        "anomalous_local_design": {
            "design_primes": list(DESIGN_PRIMES),
            "complete_signed_residue_pairs": {
                str(prime): list(values) for prime, values in ANOMALOUS_ROOTS.items()
            },
            "fixed_divisor_valuation_audit": local_audit,
            "profiles": [list(value) for value in profile_primes()],
            "trace_free_population_construction": True,
        },
        "population": {
            **raw_audit,
            "shell_radius_by_profile_size": SHELL_RADIUS_BY_PROFILE_SIZE,
            "projective_height_by_profile_size": PROJECTIVE_HEIGHT_BY_PROFILE_SIZE,
            "projective_height_at_most": PROJECTIVE_HEIGHT_MAX,
            "denominator_at_least": DENOMINATOR_MIN,
            "exact_prior_parameter_overlaps_removed": overlaps,
            "novel_population_count": len(novel),
            "novel_population_sha256": parameter_stream_sha256(
                item.parameter for item in novel
            ),
        },
        "prior_population_exclusion": exclusion_audit,
        "leakage_free_scoring": {
            "exact_discovery_trace_primes": list(DISCOVERY_PRIMES),
            "exact_held_forward_trace_primes": list(HELD_FORWARD_PRIMES),
            "construction_discovery_and_held_primes_pairwise_disjoint": True,
            "all_trace_primes_above_completed_global_exact_band": True,
            "proxy_discovery": discovery_audit,
            "discovery_frontier_count": len(frontier),
            "discovery_frontier_sha256": parameter_stream_sha256(
                item.parameter for item in frontier
            ),
            "held_forward_population_count": len(exact),
            "held_forward_population_sha256": parameter_stream_sha256(
                item.parameter for item in exact
            ),
            "held_forward_population": [candidate_record(item) for item in exact],
        },
        "conductor_and_root_parity": {
            "selected": len(conductor_candidates),
            "completed": len(completed),
            "timeouts": sum(item.status == "timeout" for item in conductors),
            "errors": sum(item.status == "error" for item in conductors),
            "subtarget_count": len(subtarget),
            "subtarget_root_minus_count": sum(
                int(item.data["root_number"]) == -1 for item in subtarget
            ),
            "records": [
                {
                    **candidate_record(item.candidate),
                    "status": item.status,
                    "conductor": item.data if item.status == "completed" else None,
                    **({"error": item.error} if item.error else {}),
                }
                for item in conductors
            ],
        },
        "point_population": {
            "selected": len(point_candidates),
            "constructor_parameters": [
                rational_to_string(item.parameter) for item in point_candidates
            ],
            "exact_subtarget_root_minus_first": True,
        },
        "point_stages": point_stages,
        "maximum_stable_numerical_rank": maximum_rank,
        "exact_checkpoints_numerical_rank_at_least_21": checkpoints,
        "bounds_and_caveats": {
            "bounded_population_is_not_a_rank_upper_bound": True,
            "radical_proxy_is_not_a_conductor_bound": True,
            "root_number_is_only_a_parity_priority": True,
            "numerical_height_rank_is_not_an_independence_proof": True,
            "all_subprocesses_foreground_joined_and_finitely_timed": True,
            "proxy_discovery_keep": PROXY_DISCOVERY_KEEP,
            "discovery_frontier_keep": DISCOVERY_FRONTIER_KEEP,
            "conductor_keep": CONDUCTOR_KEEP,
            "point_keep": POINT_KEEP,
            "point_stage_heights": list(STAGE_HEIGHTS),
            "point_stage_keeps": list(STAGE_KEEPS),
            "point_stage_timeouts_seconds": list(STAGE_TIMEOUTS),
            "certificate_trigger_stable_numerical_rank": CERTIFICATE_TRIGGER_RANK,
        },
        "software": {"python": platform.python_version(), "pari_gp": pari_version()},
        "reproduction": {
            "command": REPRODUCING_COMMAND,
            "actual_command": " ".join(
                shlex.quote(value) for value in [sys.executable, *sys.argv]
            ),
            "script_sha256": sha256_file(Path(__file__).resolve()),
            "wall_seconds": time.monotonic() - started,
        },
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--conductor-timeout", type=float, default=12.0)
    parser.add_argument("--height-timeout", type=float, default=30.0)
    parser.add_argument("--saturation-timeout", type=float, default=40.0)
    parser.add_argument("--certificate-prime-bound", type=int, default=1_000)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument(
        "--output",
        type=Path,
        default=GENERATED / "elliptic_nagao_section7_anomalous_root_crt.json",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not 1 <= args.workers <= 4:
        raise SystemExit("--workers must lie in [1,4]")
    if args.certificate_prime_bound < 2:
        raise SystemExit("invalid certificate-prime bound")
    if any(
        not 0 < value <= 120
        for value in (
            args.conductor_timeout,
            args.height_timeout,
            args.saturation_timeout,
        )
    ):
        raise SystemExit("all subprocess timeouts must lie in (0,120]")
    artifact = build_artifact(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
