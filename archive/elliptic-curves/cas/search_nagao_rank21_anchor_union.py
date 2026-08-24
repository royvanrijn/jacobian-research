#!/usr/bin/env python3
"""Broad anchor/root/trace-union search in Nagao's rank-21 family.

The population is a declared union of four independent trace-tail CRT beams,
root-power beams using *all* p-adic balls at the requested valuation, and
unconditioned smooth/even-denominator shells around the published and three
historical finalist parameters.  The known fibers are calibration only.

Every prime used to construct a trace or root condition is omitted from the
B=200 and exact PARI B=2000 scores.  Exact conductor and root-number replay
precedes point search, with conductor-feasible root-number -1 curves first.
Numerical height ranks are triage only.  Stable rank at least 18 is replayed
through saturation and exact finite reductions; the stated target still
requires a certified rank at least 21 and strict log(N)<182.72.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from math import gcd, log
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import sys
from typing import Any, Iterable, Sequence

from crt_lattice import crt_pair, gauss_reduce
from ek_k3 import primes_up_to, rational_to_string
from multiple_root_lifting import affine_variable_coefficients, fixed_divisor_valuation
from nagao_1994 import (
    PRIMARY_SOURCE,
    RANK21_CONSTRUCTION,
    RANK21_CONSTRUCTOR_PARAMETER,
    short_jacobian_coefficients,
)
from pari_bridge import minimal_curve_data, pari_version
from search_nagao_rank21_fingerprint_crt import (
    DISCRIMINANT_POLYNOMIAL,
    KNOWN_LEAD_PARAMETERS,
    conductor_radical_proxy,
    homogenized_discriminant,
    integer_valuation,
)
from search_nagao_rank21_unbiased import (
    PointPool,
    PrefilterCandidate,
    build_residue_tables,
    finite_reduction_certificate,
    parallel_point_search,
    parallel_rank_replay,
    point_pool_record,
    pool_priority,
    projective_index,
    residue_score,
)


Q = Fraction
REPOSITORY = Path(__file__).resolve().parents[2]
TARGET_LOG_CONDUCTOR = Decimal("182.72")
OLD_BOX_A_MAX = 10_000
OLD_BOX_B_MAX = 100
TRACE_COUNT = 5
TRACE_BEAM_WIDTH = 3_500
ROOT_BEAM_WIDTH = 500
PROXY_LIMIT = Decimal("205")
PROXY_SURVIVOR_CAP = 4_000
B200_KEEP = 500
CONDUCTOR_KEEP = 48
POINT_KEEP = 28
STAGE_HEIGHTS = (50_000, 250_000, 1_000_000)
STAGE_KEEPS = (12, 4)
STAGE_TIMEOUTS = (10.0, 30.0, 90.0)
CHECKPOINT_RANK = 18
TARGET_RANK = 21
ROOT_TARGETS = ((5, 8), (7, 6), (11, 5), (13, 4), (17, 3), (23, 3), (37, 2), (83, 2))
ROOT_PRIMES = tuple(prime for prime, _ in ROOT_TARGETS)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_rank21_anchor_union.py"
)


@dataclass(frozen=True)
class Anchor:
    label: str
    parameter: Fraction
    status: str


ANCHORS = (
    Anchor("published-rank21", RANK21_CONSTRUCTOR_PARAMETER, "Nagao rank at least 21"),
    Anchor("historical-1393", Q(1393, 108), "historical finalist; bounded stable rank 17"),
    Anchor("historical-1649", Q(1649, 6), "historical finalist; bounded stable rank 17"),
    Anchor("historical-6629", Q(6629, 174), "finite-reduction rank at least 18"),
)

# (prime, anchor residue, anchor trace, affine size of trace<=anchor-trace union).
# These are selected automatically by ``learn_trace_fingerprint`` and pinned
# so a change in the family or scoring implementation cannot silently alter
# the designed population.
EXPECTED_TRACE_FINGERPRINTS = {
    "published-rank21": (
        (67, 53, -16, 4),
        (83, 24, -18, 4),
        (97, 20, -18, 8),
        (103, 92, -19, 4),
        (127, 106, -22, 6),
    ),
    "historical-1393": (
        (19, 18, -8, 2),
        (73, 69, -17, 2),
        (103, 52, -17, 10),
        (139, 112, -23, 4),
        (149, 6, -24, 8),
    ),
    "historical-1649": (
        (53, 1, -14, 10),
        (67, 18, -14, 14),
        (79, 51, -15, 12),
        (127, 42, -19, 24),
        (167, 80, -24, 12),
    ),
    "historical-6629": (
        (107, 67, -20, 6),
        (113, 92, -21, 10),
        (157, 39, -25, 4),
        (181, 139, -25, 10),
        (191, 37, -27, 8),
    ),
}
TRACE_DESIGN_PRIMES = tuple(
    sorted({item[0] for values in EXPECTED_TRACE_FINGERPRINTS.values() for item in values})
)
DESIGN_PRIMES = tuple(sorted(set(TRACE_DESIGN_PRIMES) | set(ROOT_PRIMES)))


PRIOR_ARTIFACTS = {
    "elliptic_nagao_rank21_neighborhood.json": "7d59fe9a91c0f3e46604794e8931ae27e26eeea1ebf252176dffd6be8d6010fe",
    "elliptic_nagao_rank21_mutations.json": "0c0347246a1f154831c77c8f1334eb202493a0f9501caada255877a5429dd330",
    "elliptic_nagao_rank21_unbiased.json": "5bf7406855af5ec39b269fa4105c9225adb4a10d13fab5480b15264cc3e8fe1d",
    "elliptic_nagao_rank21_fingerprint_crt.json": "31ff52558551f5799633f7d3767af017a036180e34cb730caa50dc4527f547b1",
    "elliptic_nagao_rank21_neighbor_triage.json": "9b3151444a649fe1bc5a52c58a170a362f60f97bc17612ee5c8503c99c0d00bc",
    "elliptic_nagao_rank21_historical_finalists.json": "90fc658cdb7c39c96317ee888be1364b8c9f368859230e25161dc45cd6a3cec7",
}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_rational(text: str) -> Fraction | None:
    try:
        return abs(Q(text))
    except (ValueError, ZeroDivisionError):
        return None


def prior_parameter_exclusions() -> tuple[frozenset[Fraction], dict[str, Any]]:
    """Load every constructor parameter stored by the pinned prior lanes."""

    excluded = {abs(Q(value)) for value in KNOWN_LEAD_PARAMETERS}
    excluded.update(abs(anchor.parameter) for anchor in ANCHORS)
    per_artifact = {}

    def visit(value: Any) -> None:
        if isinstance(value, dict):
            for key, item in value.items():
                if (
                    isinstance(item, str)
                    and (
                        key.startswith("constructor_parameter")
                        or key in {"calibration_parameter", "lead_parameter"}
                    )
                ):
                    parsed = parse_rational(item)
                    if parsed is not None:
                        excluded.add(parsed)
                visit(item)
        elif isinstance(value, list):
            for item in value:
                visit(item)

    result_dir = REPOSITORY / "artifacts/generated-results"
    for name, expected_sha in PRIOR_ARTIFACTS.items():
        path = result_dir / name
        observed_sha = sha256_file(path)
        if observed_sha != expected_sha:
            raise AssertionError(f"prior exclusion artifact changed: {name}")
        before = len(excluded)
        visit(json.loads(path.read_text(encoding="utf-8")))
        per_artifact[name] = {
            "sha256": observed_sha,
            "new_exact_parameter_exclusions": len(excluded) - before,
        }
    return frozenset(excluded), {
        "pinned_artifacts": per_artifact,
        "total_exact_parameter_exclusions": len(excluded),
        "old_box_excluded_separately": [OLD_BOX_A_MAX, OLD_BOX_B_MAX],
    }


def learn_trace_fingerprint(anchor: Anchor, tables: dict[int, Any]) -> tuple[tuple[int, int, int, int], ...]:
    ranked = []
    for prime, table in tables.items():
        if prime < 11:
            continue
        index = projective_index(anchor.parameter.numerator, anchor.parameter.denominator, prime)
        symbol = table[index]
        if not symbol.good_reduction or symbol.ellap is None or symbol.ellap >= 0:
            continue
        union_size = sum(
            candidate.projective_index < prime
            and candidate.good_reduction
            and candidate.ellap is not None
            and candidate.ellap <= symbol.ellap
            for candidate in table
        )
        information_weight = symbol.contribution * log(prime / union_size)
        ranked.append((-information_weight, prime, index, int(symbol.ellap), union_size))
    answer = tuple(sorted(item[1:] for item in sorted(ranked)[:TRACE_COUNT]))
    if answer != EXPECTED_TRACE_FINGERPRINTS[anchor.label]:
        raise AssertionError(f"trace fingerprint changed for {anchor.label}")
    return answer


@dataclass(frozen=True)
class TraceChoice:
    prime: int
    residue: int
    ellap: int
    threshold: int
    contribution: float


@dataclass(frozen=True)
class RootBall:
    prime: int
    residue: int
    modulus: int
    forced_valuation: int


@dataclass(frozen=True)
class BeamState:
    residue: int
    modulus: int
    trace_choices: tuple[TraceChoice, ...]
    root_balls: tuple[RootBall, ...]
    representative: Fraction
    basis: tuple[tuple[int, int], tuple[int, int]]
    objective: float


def polynomial_value_mod(coefficients: Sequence[int], value: int, modulus: int) -> int:
    answer = 0
    for coefficient in reversed(coefficients):
        answer = (answer * value + coefficient) % modulus
    return answer


def root_ball_union(prime: int, target_valuation: int) -> tuple[RootBall, ...]:
    """Classify all disjoint affine balls forcing the requested valuation."""

    active = [
        (residue, prime)
        for residue in range(prime)
        if polynomial_value_mod(DISCRIMINANT_POLYNOMIAL, residue, prime) == 0
    ]
    completed = []
    while active:
        residue, modulus = active.pop()
        forced = fixed_divisor_valuation(
            affine_variable_coefficients(DISCRIMINANT_POLYNOMIAL, residue, modulus),
            prime,
        )
        if forced >= target_valuation:
            completed.append(RootBall(prime, residue, modulus, forced))
            continue
        if modulus > prime ** (target_valuation + 2):
            raise AssertionError("a p-adic root branch did not stabilize")
        lifted_modulus = modulus * prime
        for lift in range(prime):
            lifted = residue + lift * modulus
            if polynomial_value_mod(DISCRIMINANT_POLYNOMIAL, lifted, lifted_modulus) == 0:
                active.append((lifted, lifted_modulus))
    return tuple(sorted(completed, key=lambda ball: (ball.modulus, ball.residue)))


def gauss_shell(
    residue: int,
    modulus: int,
    *,
    radius: int = 4,
    limit: int = 6,
) -> tuple[tuple[Fraction, tuple[int, int], tuple[tuple[int, int], tuple[int, int]]], ...]:
    basis = gauss_reduce((modulus, 0), (residue, 1))
    found: dict[Fraction, tuple[int, int]] = {}
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
            if (
                (parameter.numerator - residue * parameter.denominator) % modulus
                and (parameter.numerator + residue * parameter.denominator) % modulus
            ):
                continue
            pair = (left, right)
            if parameter not in found or pair < found[parameter]:
                found[parameter] = pair
    ordered = sorted(
        found.items(),
        key=lambda item: (
            max(item[0].numerator, item[0].denominator),
            item[0].numerator,
            item[0].denominator,
            item[1],
        ),
    )
    return tuple((parameter, pair, basis) for parameter, pair in ordered[:limit])


def make_state(
    residue: int,
    modulus: int,
    traces: tuple[TraceChoice, ...],
    roots: tuple[RootBall, ...],
) -> BeamState:
    shell = gauss_shell(residue, modulus, radius=3, limit=1)
    if not shell:
        shell = gauss_shell(residue, modulus, radius=7, limit=1)
    if not shell:
        raise AssertionError("CRT lattice had no affine-unit representative")
    representative, _, basis = shell[0]
    height = max(representative.numerator, representative.denominator)
    trace_reward = sum(choice.contribution for choice in traces)
    root_savings = sum((ball.forced_valuation - 1) * log(ball.prime) for ball in roots)
    return BeamState(
        residue,
        modulus,
        traces,
        roots,
        representative,
        basis,
        20 * log(height) - 2 * trace_reward - root_savings,
    )


def retain_beam(states: Iterable[BeamState], width: int) -> tuple[BeamState, ...]:
    best: dict[tuple[int, int], BeamState] = {}
    for state in states:
        key = state.residue, state.modulus
        prior = best.get(key)
        if prior is None or (
            state.objective,
            max(state.representative.numerator, state.representative.denominator),
            state.residue,
        ) < (
            prior.objective,
            max(prior.representative.numerator, prior.representative.denominator),
            prior.residue,
        ):
            best[key] = state
    return tuple(
        sorted(
            best.values(),
            key=lambda state: (
                state.objective,
                max(state.representative.numerator, state.representative.denominator),
                state.residue,
            ),
        )[:width]
    )


def build_trace_beam(
    fingerprint: Sequence[tuple[int, int, int, int]],
    tables: dict[int, Any],
    *,
    width: int,
) -> tuple[tuple[BeamState, ...], tuple[BeamState, ...], tuple[dict[str, Any], ...]]:
    states = (make_state(0, 1, (), ()),)
    prefix = states
    audit = []
    for stage, (prime, _, threshold, expected_union) in enumerate(fingerprint, start=1):
        choices = tuple(
            TraceChoice(prime, symbol.projective_index, int(symbol.ellap), threshold, symbol.contribution)
            for symbol in tables[prime]
            if symbol.projective_index < prime
            and symbol.good_reduction
            and symbol.ellap is not None
            and symbol.ellap <= threshold
        )
        if len(choices) != expected_union:
            raise AssertionError("a pinned trace union changed")
        expanded = []
        for state in states:
            for choice in choices:
                residue, modulus = crt_pair(state.residue, state.modulus, choice.residue, prime)
                expanded.append(make_state(residue, modulus, state.trace_choices + (choice,), state.root_balls))
        states = retain_beam(expanded, width)
        if stage == 3:
            prefix = states
        audit.append(
            {
                "prime": prime,
                "anchor_trace_threshold": threshold,
                "affine_symbol_union_size": len(choices),
                "expanded": len(expanded),
                "retained": len(states),
            }
        )
    return states, prefix, tuple(audit)


def extend_root_beam(
    states: Sequence[BeamState], balls: Sequence[RootBall], *, width: int
) -> tuple[BeamState, ...]:
    expanded = []
    for state in states:
        for ball in balls:
            if gcd(state.modulus, ball.modulus) != 1:
                continue
            residue, modulus = crt_pair(state.residue, state.modulus, ball.residue, ball.modulus)
            expanded.append(make_state(residue, modulus, state.trace_choices, state.root_balls + (ball,)))
    return retain_beam(expanded, width)


def build_strata(
    trace_data: dict[str, tuple[tuple[BeamState, ...], tuple[BeamState, ...], tuple[dict[str, Any], ...]]],
    root_unions: dict[int, tuple[RootBall, ...]],
) -> dict[str, tuple[BeamState, ...]]:
    strata = {}
    profiles = {
        "published-rank21": ((13, 17, 23), (5, 13)),
        "historical-1393": ((5, 11, 23), (7, 11)),
        "historical-1649": ((5, 7, 11), (23, 83)),
        "historical-6629": ((5, 7, 23), (11, 23)),
    }
    for label, (full, prefix, _) in trace_data.items():
        strata[f"{label}-trace-full"] = full
        # Each full beam also receives its strongest family-wide alternative
        # root union, while prefix beams receive multi-prime saving profiles.
        primary_prime = profiles[label][0][0]
        strata[f"{label}-trace-full-root-{primary_prime}"] = extend_root_beam(
            full, root_unions[primary_prime], width=ROOT_BEAM_WIDTH
        )
        for profile_index, primes in enumerate(profiles[label], start=1):
            states = prefix
            for prime in primes:
                states = extend_root_beam(states, root_unions[prime], width=ROOT_BEAM_WIDTH)
            strata[f"{label}-trace-prefix-root-profile-{profile_index}"] = states
    return strata


@dataclass(frozen=True)
class GeneratedCandidate:
    parameter: Fraction
    stratum: str
    trace_conditions: tuple[tuple[int, int, int], ...]
    root_conditions: tuple[tuple[int, int, int, int], ...]
    proxy: dict[str, Any]

    @property
    def identifier(self) -> str:
        return f"anchor-union-{self.parameter.numerator}-{self.parameter.denominator}"

    @property
    def height(self) -> int:
        return max(self.parameter.numerator, self.parameter.denominator)


def outside_old_box(parameter: Fraction) -> bool:
    parameter = abs(Q(parameter))
    return not (
        parameter.numerator <= OLD_BOX_A_MAX
        and parameter.denominator <= OLD_BOX_B_MAX
    )


def nearest_integer(value: Fraction) -> int:
    quotient, remainder = divmod(value.numerator, value.denominator)
    return quotient + (2 * remainder >= value.denominator)


def smooth_denominators() -> tuple[int, ...]:
    values = set()
    for power_two in range(14):
        for power_three in range(9):
            for extra in (1, 5, 7, 11, 13):
                denominator = 2**power_two * 3**power_three * extra
                if 101 <= denominator <= 5_000:
                    values.add(denominator)
    return tuple(sorted(values))


def raw_parameters(strata: dict[str, tuple[BeamState, ...]]) -> tuple[tuple[Fraction, str, tuple[tuple[int, int, int], ...], tuple[tuple[int, int, int, int], ...]], ...]:
    raw = []
    for label, states in strata.items():
        for state in states:
            shell = gauss_shell(state.residue, state.modulus, radius=4, limit=6)
            indices = (0, 2, 5) if "trace-full" in label and "root" not in label else (0, 3)
            for index in indices:
                if index >= len(shell):
                    continue
                parameter = shell[index][0]
                traces = tuple(
                    (
                        choice.prime,
                        projective_index(parameter.numerator, parameter.denominator, choice.prime),
                        choice.threshold,
                    )
                    for choice in state.trace_choices
                )
                roots = tuple(
                    (
                        ball.prime,
                        parameter.numerator
                        * pow(parameter.denominator, -1, ball.modulus)
                        % ball.modulus,
                        ball.modulus,
                        ball.forced_valuation,
                    )
                    for ball in state.root_balls
                )
                raw.append((parameter, label, traces, roots))

    # Unconditioned controls: smooth denominators and all 16-divisible
    # denominators near each anchor.  These are intentionally not trace-scored
    # during construction and therefore test the beam's selection assumptions.
    for anchor in ANCHORS:
        for denominator in smooth_denominators():
            center = nearest_integer(anchor.parameter * denominator)
            for offset in range(-96, 97, 8):
                numerator = center + offset
                if numerator <= 0:
                    continue
                raw.append((Q(numerator, denominator), f"{anchor.label}-smooth-denominator", (), ()))
        for denominator in range(112, 2_001, 16):
            center = nearest_integer(anchor.parameter * denominator)
            for offset in range(-80, 81, 16):
                numerator = center + offset
                if numerator <= 0:
                    continue
                if numerator % 2 == 0:
                    numerator += 1
                raw.append((Q(numerator, denominator), f"{anchor.label}-even-denominator", (), ()))
    return tuple(raw)


def generate_candidates(
    strata: dict[str, tuple[BeamState, ...]],
    exclusions: frozenset[Fraction],
    *,
    proxy_limit: Decimal,
    survivor_cap: int,
) -> tuple[tuple[GeneratedCandidate, ...], dict[str, Any]]:
    raw = raw_parameters(strata)
    unique_raw: dict[Fraction, tuple[str, tuple[tuple[int, int, int], ...], tuple[tuple[int, int, int, int], ...]]] = {}
    old_box = prior = zero = 0
    for parameter, stratum, traces, roots in raw:
        parameter = abs(Q(parameter))
        if parameter == 0:
            zero += 1
            continue
        if not outside_old_box(parameter):
            old_box += 1
            continue
        if parameter in exclusions:
            prior += 1
            continue
        record = (stratum, traces, roots)
        existing = unique_raw.get(parameter)
        if existing is None or (len(roots), len(traces), stratum) > (
            len(existing[2]), len(existing[1]), existing[0]
        ):
            unique_raw[parameter] = record

    retained = []
    proxy_rejected = singular = 0
    for parameter, (stratum, traces, roots) in unique_raw.items():
        try:
            proxy = conductor_radical_proxy(parameter)
        except ValueError:
            singular += 1
            continue
        if Decimal(str(proxy["log_radical_upper_proxy"])) >= proxy_limit:
            proxy_rejected += 1
            continue
        retained.append(GeneratedCandidate(parameter, stratum, traces, roots, proxy))

    retained.sort(key=lambda candidate: (candidate.proxy["log_radical_upper_proxy"], candidate.height, candidate.identifier))
    # Exact per-stratum quotas prevent dense smooth shells from erasing CRT arms.
    selected: dict[Fraction, GeneratedCandidate] = {}
    labels = sorted({candidate.stratum for candidate in retained})
    quota = max(1, min(30, survivor_cap // max(1, len(labels))))
    for label in labels:
        for candidate in (item for item in retained if item.stratum == label):
            selected[candidate.parameter] = candidate
            if sum(item.stratum == label for item in selected.values()) >= quota:
                break
    for candidate in retained:
        if len(selected) >= survivor_cap:
            break
        selected.setdefault(candidate.parameter, candidate)
    answer = tuple(sorted(selected.values(), key=lambda candidate: (candidate.proxy["log_radical_upper_proxy"], candidate.height, candidate.identifier)))
    if any(candidate.parameter in exclusions or not outside_old_box(candidate.parameter) for candidate in answer):
        raise AssertionError("a prior parameter leaked into the anchor-union population")
    digest = hashlib.sha256()
    for candidate in answer:
        digest.update(f"{candidate.parameter}|{candidate.stratum}|{candidate.proxy['log_radical_upper_proxy']!r}\n".encode())
    return answer, {
        "raw_records": len(raw),
        "exactly_deduplicated_outside_prior_population": len(unique_raw),
        "excluded_old_box": old_box,
        "excluded_exact_prior_fibers": prior,
        "excluded_zero": zero,
        "singular": singular,
        "proxy_rejected": proxy_rejected,
        "below_proxy_before_cap": len(retained),
        "retained_after_stratified_cap": len(answer),
        "stratum_counts": {label: sum(candidate.stratum == label for candidate in answer) for label in labels},
        "survivor_stream_sha256": digest.hexdigest(),
    }


@dataclass(frozen=True)
class ExactCandidate:
    prefilter: PrefilterCandidate
    generated: GeneratedCandidate
    exact_score_b2000: str
    exact_good_primes: int
    exact_bad_primes: int
    exact_last_prime: int
    omitted_design_primes: int

    @property
    def parameter(self) -> Fraction:
        return self.prefilter.parameter

    @property
    def identifier(self) -> str:
        return self.generated.identifier


def prefilter_candidates(candidates: Sequence[GeneratedCandidate], tables: dict[int, Any], *, keep: int) -> tuple[tuple[PrefilterCandidate, GeneratedCandidate], ...]:
    scored = []
    for candidate in candidates:
        score, good, bad = residue_score(candidate.parameter.numerator, candidate.parameter.denominator, tables)
        scored.append((PrefilterCandidate(candidate.parameter.numerator, candidate.parameter.denominator, score, good, bad), candidate))
    scored.sort(key=lambda item: (-item[0].residue_score_b200, item[1].proxy["log_radical_upper_proxy"], item[1].height, item[1].identifier))
    return tuple(scored[: min(keep, len(scored))])


def gp_rational(value: Fraction) -> str:
    return f"({value.numerator}/{value.denominator})"


def exact_decontaminated_scores(
    candidates: Sequence[tuple[PrefilterCandidate, GeneratedCandidate]],
    *,
    cutoff: int,
    batch_size: int,
    timeout: float,
    stack_bytes: int,
) -> tuple[ExactCandidate, ...]:
    executable = shutil.which("gp")
    if executable is None:
        raise FileNotFoundError("PARI/GP executable 'gp' was not found")
    omitted_test = "||".join(f"p=={prime}" for prime in DESIGN_PRIMES)
    records = {}
    last_prime = primes_up_to(cutoff)[-1]
    for start in range(0, len(candidates), batch_size):
        batch = candidates[start : start + batch_size]
        commands = ["default(realprecision,80);"]
        for index, (prefilter, _) in enumerate(batch):
            coefficients = short_jacobian_coefficients(RANK21_CONSTRUCTION, prefilter.parameter)
            vector = ",".join(gp_rational(value) for value in coefficients)
            commands.extend(
                (
                    f"E=ellminimalmodel(ellinit([{vector}]));",
                    "S=0;USED=0;BAD=0;OMITTED=0;",
                    f"forprime(p=5,{cutoff},if({omitted_test},OMITTED++,if(valuation(E.disc,p)>0,BAD++,A=ellap(E,p);S+=(2-A)/(p+1-A)*log(p);USED++)));",
                    f'print("ROW|{index}|",S,"|",USED,"|",BAD,"|",OMITTED);',
                )
            )
        commands.append("quit")
        result = subprocess.run(
            [executable, "-q", "-s", str(stack_bytes)],
            input="\n".join(commands) + "\n",
            text=True,
            capture_output=True,
            timeout=timeout,
        )
        if result.returncode != 0 or "***" in result.stderr:
            raise RuntimeError(f"PARI exact score failed: {result.stderr.strip()}")
        observed = 0
        for line in result.stdout.splitlines():
            if not line.startswith("ROW|"):
                continue
            _, index_text, score, used, bad, omitted = line.split("|")
            prefilter, generated = batch[int(index_text)]
            records[generated.identifier] = (score, int(used), int(bad), int(omitted), prefilter, generated)
            observed += 1
        if observed != len(batch):
            raise RuntimeError("PARI omitted exact score rows")
    answer = tuple(
        ExactCandidate(prefilter, generated, score, used, bad, last_prime, omitted)
        for score, used, bad, omitted, prefilter, generated in records.values()
    )
    return tuple(sorted(answer, key=lambda candidate: (-Decimal(candidate.exact_score_b2000), candidate.generated.proxy["log_radical_upper_proxy"], candidate.generated.height, candidate.identifier)))


@dataclass(frozen=True)
class ConductorReplay:
    candidate: ExactCandidate
    status: str
    data: dict[str, Any]
    error: str | None = None


def replay_one_conductor(candidate: ExactCandidate, *, timeout: float, stack_bytes: int) -> ConductorReplay:
    try:
        data = minimal_curve_data(
            short_jacobian_coefficients(RANK21_CONSTRUCTION, candidate.parameter),
            timeout=timeout,
            stack_bytes=stack_bytes,
            local_primes=(2, 3, *ROOT_PRIMES),
        )
        data["below_strict_log_conductor_target"] = Decimal(data["log_conductor"]) < TARGET_LOG_CONDUCTOR
        return ConductorReplay(candidate, "completed", data)
    except (subprocess.TimeoutExpired, RuntimeError, AssertionError, ValueError) as error:
        return ConductorReplay(candidate, "timeout" if isinstance(error, subprocess.TimeoutExpired) else "error", {}, str(error)[:500])


def parallel_conductors(candidates: Sequence[ExactCandidate], *, timeout: float, stack_bytes: int, workers: int) -> tuple[ConductorReplay, ...]:
    by_id = {}
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(replay_one_conductor, candidate, timeout=timeout, stack_bytes=stack_bytes): candidate.identifier
            for candidate in candidates
        }
        for future in as_completed(futures):
            by_id[futures[future]] = future.result()
    return tuple(by_id[candidate.identifier] for candidate in candidates)


def select_conductor_population(candidates: Sequence[ExactCandidate], *, keep: int) -> tuple[ExactCandidate, ...]:
    selected = {candidate.identifier: candidate for candidate in candidates[:24]}
    by_proxy = sorted(candidates, key=lambda candidate: (candidate.generated.proxy["log_radical_upper_proxy"], -Decimal(candidate.exact_score_b2000), candidate.identifier))
    for candidate in by_proxy[:16]:
        selected[candidate.identifier] = candidate
    for stratum in sorted({candidate.generated.stratum for candidate in candidates}):
        leader = next(candidate for candidate in candidates if candidate.generated.stratum == stratum)
        selected[leader.identifier] = leader
    for candidate in candidates:
        if len(selected) >= keep:
            break
        selected.setdefault(candidate.identifier, candidate)
    return tuple(sorted(selected.values(), key=lambda candidate: (-Decimal(candidate.exact_score_b2000), candidate.identifier))[:keep])


def select_point_population(replays: Sequence[ConductorReplay], *, keep: int) -> tuple[ExactCandidate, ...]:
    completed = [replay for replay in replays if replay.status == "completed"]
    root_minus = sorted(
        (replay for replay in completed if replay.data.get("below_strict_log_conductor_target") and replay.data.get("root_number") == -1),
        key=lambda replay: (-Decimal(replay.candidate.exact_score_b2000), replay.candidate.identifier),
    )
    feasible_plus = sorted(
        (replay for replay in completed if replay.data.get("below_strict_log_conductor_target") and replay.data.get("root_number") == 1),
        key=lambda replay: (-Decimal(replay.candidate.exact_score_b2000), replay.candidate.identifier),
    )
    selected = {replay.candidate.identifier: replay.candidate for replay in root_minus[:keep]}
    for replay in feasible_plus[:8]:
        if len(selected) >= keep:
            break
        selected.setdefault(replay.candidate.identifier, replay.candidate)
    for replay in completed:
        if len(selected) >= keep:
            break
        selected.setdefault(replay.candidate.identifier, replay.candidate)
    root_minus_ids = {replay.candidate.identifier for replay in root_minus}
    return tuple(sorted(selected.values(), key=lambda candidate: (candidate.identifier not in root_minus_ids, -Decimal(candidate.exact_score_b2000), candidate.identifier))[:keep])


def generated_record(candidate: ExactCandidate) -> dict[str, Any]:
    return {
        "candidate_id": candidate.identifier,
        "constructor_parameter": rational_to_string(candidate.parameter),
        "height": candidate.generated.height,
        "stratum": candidate.generated.stratum,
        "trace_conditions": [list(value) for value in candidate.generated.trace_conditions],
        "root_conditions": [list(value) for value in candidate.generated.root_conditions],
        "radical_proxy": candidate.generated.proxy,
        "exact_small_prime_discriminant_valuations": {
            str(prime): integer_valuation(homogenized_discriminant(candidate.parameter), prime)
            for prime in (2, 3, *ROOT_PRIMES)
        },
        "leakage_free_b200": {
            "score": candidate.prefilter.residue_score_b200,
            "good_primes_used": candidate.prefilter.residue_good_primes,
            "bad_primes_skipped": candidate.prefilter.residue_bad_primes,
            "design_primes_omitted": list(DESIGN_PRIMES),
        },
        "exact_pari_leakage_free_b2000": {
            "score": candidate.exact_score_b2000,
            "good_primes_used": candidate.exact_good_primes,
            "bad_primes_skipped": candidate.exact_bad_primes,
            "design_primes_omitted_count": candidate.omitted_design_primes,
            "last_numerical_prime": candidate.exact_last_prime,
        },
    }


def staged_points(candidates: Sequence[ExactCandidate], args: argparse.Namespace) -> tuple[list[dict[str, Any]], dict[str, tuple[PointPool, dict[str, Any]]]]:
    retained = tuple(candidates)
    best = {}
    stages = []
    for index, (height, timeout) in enumerate(zip(STAGE_HEIGHTS, STAGE_TIMEOUTS), start=1):
        pools = parallel_point_search(retained, height_bound=height, timeout=timeout, stack_bytes=args.stack_bytes, workers=args.workers)
        ranks = parallel_rank_replay(pools, precisions=(72, 120), timeout=args.height_timeout, stack_bytes=args.stack_bytes, workers=args.workers)
        ranked = sorted(((pool, ranks[pool.candidate.identifier]) for pool in pools), key=lambda item: pool_priority(item[0], item[1]))
        for pool, rank in ranked:
            if rank.get("status") == "completed":
                prior = best.get(pool.candidate.identifier)
                if prior is None or int(rank["stable_numerical_rank"]) >= int(prior[1]["stable_numerical_rank"]):
                    best[pool.candidate.identifier] = (pool, rank)
        keep = len(ranked) if index == len(STAGE_HEIGHTS) else STAGE_KEEPS[index - 1]
        retained = tuple(pool.candidate for pool, _ in ranked[:keep])
        stages.append(
            {
                "stage": index,
                "quartic_naive_height_bound": height,
                "population_searched": len(ranked),
                "completed": sum(pool.status == "completed" for pool, _ in ranked),
                "timeouts": sum(pool.status == "timeout" for pool, _ in ranked),
                "errors": sum(pool.status == "error" for pool, _ in ranked),
                "retained_constructor_parameters": [rational_to_string(candidate.parameter) for candidate in retained],
                "ranked_population": [
                    {**point_pool_record(pool, rank, include_points=index == len(STAGE_HEIGHTS)), "anchor_union_generation": generated_record(pool.candidate)}
                    for pool, rank in ranked
                ],
            }
        )
        top_rank = max((int(rank["stable_numerical_rank"]) for _, rank in ranked if rank.get("status") == "completed"), default=-1)
        print(f"stage H={height} population={len(ranked)} max_rank={top_rank}", flush=True)
    return stages, best


def exact_checkpoints(best: dict[str, tuple[PointPool, dict[str, Any]]], conductor_by_id: dict[str, ConductorReplay], args: argparse.Namespace) -> list[dict[str, Any]]:
    records = []
    for identifier, (pool, rank) in sorted(best.items(), key=lambda item: (-int(item[1][1]["stable_numerical_rank"]), item[0])):
        numerical_rank = int(rank["stable_numerical_rank"])
        if numerical_rank < CHECKPOINT_RANK:
            continue
        replay = conductor_by_id.get(identifier)
        if replay is None or replay.status != "completed":
            continue
        certificate = finite_reduction_certificate(pool, rank, saturation_timeout=args.saturation_timeout, certificate_prime_bound=args.certificate_prime_bound, stack_bytes=args.stack_bytes)
        exact_rank = certificate["certified_algebraic_rank_lower_bound"]
        target_hit = bool(exact_rank is not None and exact_rank >= TARGET_RANK and replay.data["below_strict_log_conductor_target"])
        record = {
            "constructor_parameter": rational_to_string(pool.candidate.parameter),
            "stable_numerical_rank": numerical_rank,
            "deepest_completed_height": pool.height_bound,
            "conductor": replay.data,
            "exact_rank_certificate": certificate,
            "target_rank21_under_log_conductor_hit": target_hit,
        }
        records.append(record)
        print(f"checkpoint T={record['constructor_parameter']} numerical={numerical_rank} exact={exact_rank} target={target_hit}", flush=True)
    return records


def build_artifact(args: argparse.Namespace) -> dict[str, Any]:
    exclusions, exclusion_audit = prior_parameter_exclusions()
    tables_all = build_residue_tables(200)
    trace_data = {}
    trace_audit = {}
    for anchor in ANCHORS:
        fingerprint = learn_trace_fingerprint(anchor, tables_all)
        trace_data[anchor.label] = build_trace_beam(fingerprint, tables_all, width=args.trace_beam_width)
        trace_audit[anchor.label] = {
            "constructor_parameter": rational_to_string(anchor.parameter),
            "source_status": anchor.status,
            "fingerprint": [list(value) for value in fingerprint],
            "beam_stages": list(trace_data[anchor.label][2]),
        }
    root_unions = {prime: root_ball_union(prime, target) for prime, target in ROOT_TARGETS}
    expected_counts = {5: 16, 7: 16, 11: 6, 13: 6, 17: 6, 23: 12, 37: 3, 83: 6}
    if {prime: len(balls) for prime, balls in root_unions.items()} != expected_counts:
        raise AssertionError("a full p-adic root-ball union changed")
    strata = build_strata(trace_data, root_unions)
    generated, generation_audit = generate_candidates(strata, exclusions, proxy_limit=args.proxy_limit, survivor_cap=args.proxy_survivor_cap)
    scoring_tables = {prime: table for prime, table in tables_all.items() if prime not in DESIGN_PRIMES}
    prefiltered = prefilter_candidates(generated, scoring_tables, keep=args.b200_keep)
    exact = exact_decontaminated_scores(prefiltered, cutoff=2_000, batch_size=args.score_batch_size, timeout=args.score_timeout, stack_bytes=args.stack_bytes)
    conductor_candidates = select_conductor_population(exact, keep=args.conductor_keep)
    conductors = parallel_conductors(conductor_candidates, timeout=args.conductor_timeout, stack_bytes=args.stack_bytes, workers=args.workers)
    conductor_by_id = {replay.candidate.identifier: replay for replay in conductors}
    point_candidates = select_point_population(conductors, keep=args.point_keep)
    point_stages, best = staged_points(point_candidates, args)
    checkpoints = exact_checkpoints(best, conductor_by_id, args)
    hits = [record for record in checkpoints if record["target_rank21_under_log_conductor_hit"]]
    return {
        "schema_version": 1,
        "status": "bounded multi-anchor root/trace-union search; numerical ranks are triage only",
        "primary_source": PRIMARY_SOURCE,
        "target": {
            "rank_at_least": TARGET_RANK,
            "strict_log_conductor_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "certified_hits": hits,
        },
        "prior_population_exclusion": exclusion_audit,
        "trace_design": trace_audit,
        "design_primes_omitted_from_scores": list(DESIGN_PRIMES),
        "root_ball_unions": {
            str(prime): {
                "target_valuation": dict(ROOT_TARGETS)[prime],
                "complete_disjoint_affine_union": [ball.__dict__ for ball in balls],
            }
            for prime, balls in root_unions.items()
        },
        "generation": generation_audit,
        "exact_b2000_population": [generated_record(candidate) for candidate in exact],
        "conductor_population": [
            {
                **generated_record(replay.candidate),
                "status": replay.status,
                "conductor": replay.data if replay.status == "completed" else None,
                "error": replay.error,
            }
            for replay in conductors
        ],
        "point_population_selection": {
            "count": len(point_candidates),
            "root_minus_one_and_subtarget_count": sum(
                conductor_by_id[candidate.identifier].data.get("root_number") == -1
                and conductor_by_id[candidate.identifier].data.get("below_strict_log_conductor_target") is True
                for candidate in point_candidates
                if conductor_by_id[candidate.identifier].status == "completed"
            ),
            "constructor_parameters": [rational_to_string(candidate.parameter) for candidate in point_candidates],
        },
        "point_stages": point_stages,
        "exact_checkpoints_stable_numerical_rank_at_least_18": checkpoints,
        "bounds_and_caveats": {
            "old_positive_box_excluded": {"numerator_at_most": OLD_BOX_A_MAX, "denominator_at_most": OLD_BOX_B_MAX},
            "proxy_limit": str(args.proxy_limit),
            "proxy_survivor_cap": args.proxy_survivor_cap,
            "b200_keep": args.b200_keep,
            "conductor_keep": args.conductor_keep,
            "point_keep": args.point_keep,
            "stage_heights": list(STAGE_HEIGHTS),
            "stage_keeps": list(STAGE_KEEPS),
            "stage_timeouts_seconds": list(STAGE_TIMEOUTS),
            "bounded_search_is_not_a_rank_upper_bound": True,
            "root_number_priority_is_a_parity_heuristic": True,
            "all_subprocesses_synchronous_with_finite_timeouts": True,
        },
        "software": {"python": platform.python_version(), "pari_gp": pari_version()},
        "reproducing_command": REPRODUCING_COMMAND,
        "actual_command": " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv]),
        "script_sha256": sha256_file(Path(__file__).resolve()),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trace-beam-width", type=int, default=TRACE_BEAM_WIDTH)
    parser.add_argument("--proxy-limit", type=Decimal, default=PROXY_LIMIT)
    parser.add_argument("--proxy-survivor-cap", type=int, default=PROXY_SURVIVOR_CAP)
    parser.add_argument("--b200-keep", type=int, default=B200_KEEP)
    parser.add_argument("--score-batch-size", type=int, default=50)
    parser.add_argument("--score-timeout", type=float, default=40.0)
    parser.add_argument("--conductor-keep", type=int, default=CONDUCTOR_KEEP)
    parser.add_argument("--conductor-timeout", type=float, default=60.0)
    parser.add_argument("--point-keep", type=int, default=POINT_KEEP)
    parser.add_argument("--height-timeout", type=float, default=30.0)
    parser.add_argument("--saturation-timeout", type=float, default=30.0)
    parser.add_argument("--certificate-prime-bound", type=int, default=1_000)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPOSITORY / "artifacts/generated-results/elliptic_nagao_rank21_anchor_union.json",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not 1 <= args.workers <= 4:
        raise SystemExit("--workers must lie in [1,4]")
    if any(value < 1 for value in (args.trace_beam_width, args.proxy_survivor_cap, args.b200_keep, args.conductor_keep, args.point_keep)):
        raise SystemExit("all population bounds must be positive")
    if any(not 0 < value <= 120 for value in (args.score_timeout, args.conductor_timeout, args.height_timeout, args.saturation_timeout)):
        raise SystemExit("all subprocess timeouts must lie in (0,120]")
    artifact = build_artifact(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
