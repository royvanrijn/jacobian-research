#!/usr/bin/env python3
"""Bounded projective CRT/Gauss constructor for q12/orbit5867 parameters.

High-Nagao symbols are selected from exact discovery-prime P^1 tables.  A
bounded beam combines one projective line modulo each prime, exact Gauss
reduction finds short vectors in the resulting congruence lattice, and bounded
reduced-basis combinations generate rational parameters far outside a height
box.  Candidates are reranked on disjoint held-out primes and deduplicated
against the prior H=10000 scan.

This is a heuristic parameter constructor, not rank evidence.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
import json
from math import gcd, lcm, log
from pathlib import Path
import shlex
import sys
from time import perf_counter
from typing import Iterable, Sequence


SCRIPTS = Path(__file__).resolve().parent
ROOT = SCRIPTS.parents[1]
sys.path.insert(0, str(ROOT / "elliptic-curves/cas"))

from crt_lattice import crt_pair, gauss_reduce  # noqa: E402
from search_h92_q12o5867_rootless_nagao import (  # noqa: E402
    Candidate,
    DEFAULT_PRIME_BLOCKS,
    SCORE_SCALE,
    build_residue_tables,
    is_prime,
    load_family_model,
    score_block,
)


DEFAULT_PRIOR = (
    ROOT
    / "artifacts/local/elkies-k3/q12o5867-rootless-nagao-cpp-h10000.json"
)
DEFAULT_OUTPUT = (
    ROOT / "artifacts/local/elkies-k3/q12o5867-rootless-nagao-crt-gauss.json"
)
DEFAULT_VALIDATION_PRIMES = tuple(prime for prime in range(199, 500) if is_prime(prime))


Vector = tuple[int, int]


@dataclass(frozen=True)
class ProjectiveChoice:
    prime: int
    projective_index: int
    score_units: int
    trace: int

    @property
    def label(self) -> str:
        return "infinity" if self.projective_index == self.prime else str(self.projective_index)


@dataclass(frozen=True)
class BeamState:
    finite_residue: int
    finite_modulus: int
    infinity_modulus: int
    choices: tuple[ProjectiveChoice, ...]
    local_score_units: int
    reduced_basis: tuple[Vector, Vector]
    representative: tuple[int, int, int]
    objective: float

    @property
    def modulus(self) -> int:
        return self.finite_modulus * self.infinity_modulus


def projective_choice_groups(
    construction_blocks: Sequence[dict[int, tuple[object, ...]]], choices_per_prime: int
) -> tuple[tuple[ProjectiveChoice, ...], ...]:
    if choices_per_prime < 1:
        raise ValueError("choices_per_prime must be positive")
    groups = []
    for block in construction_blocks:
        for prime, table in block.items():
            ranked = sorted(
                (symbol for symbol in table if symbol.good_reduction),
                key=lambda symbol: (-symbol.contribution_units, symbol.projective_index),
            )[:choices_per_prime]
            if not ranked:
                raise ValueError(f"discovery prime p={prime} has no good symbols")
            groups.append(
                tuple(
                    ProjectiveChoice(
                        prime=prime,
                        projective_index=symbol.projective_index,
                        score_units=symbol.contribution_units,
                        trace=int(symbol.trace),
                    )
                    for symbol in ranked
                )
            )
    return tuple(groups)


def congruence_basis(
    finite_residue: int, finite_modulus: int, infinity_modulus: int
) -> tuple[Vector, Vector]:
    """Basis for a=R*b mod F and b=0 mod I, with gcd(F,I)=1."""

    if finite_modulus < 1 or infinity_modulus < 1:
        raise ValueError("projective CRT moduli must be positive")
    if gcd(finite_modulus, infinity_modulus) != 1:
        raise ValueError("finite and infinity moduli must be coprime")
    finite_residue %= finite_modulus
    return (finite_modulus, 0), (
        finite_residue * infinity_modulus,
        infinity_modulus,
    )


def normalize_and_validate(
    vector: Vector,
    finite_residue: int,
    finite_modulus: int,
    infinity_modulus: int,
) -> tuple[int, int, int] | None:
    numerator, denominator = vector
    if denominator == 0:
        return None
    common = gcd(abs(numerator), abs(denominator))
    if common:
        numerator //= common
        denominator //= common
    if denominator < 0:
        numerator, denominator = -numerator, -denominator
    if gcd(denominator, finite_modulus) != 1:
        return None
    if (numerator - finite_residue * denominator) % finite_modulus:
        return None
    if denominator % infinity_modulus:
        return None
    if gcd(numerator, infinity_modulus) != 1:
        return None
    return numerator, denominator, max(abs(numerator), denominator)


def short_representatives(
    finite_residue: int,
    finite_modulus: int,
    infinity_modulus: int,
    *,
    coefficient_radius: int,
) -> tuple[tuple[int, int, int], ...]:
    if coefficient_radius < 1:
        raise ValueError("coefficient_radius must be positive")
    basis = gauss_reduce(*congruence_basis(finite_residue, finite_modulus, infinity_modulus))
    answers: dict[tuple[int, int], tuple[int, int, int]] = {}
    for left in range(-coefficient_radius, coefficient_radius + 1):
        for right in range(-coefficient_radius, coefficient_radius + 1):
            if left == 0 and right == 0:
                continue
            normalized = normalize_and_validate(
                (
                    left * basis[0][0] + right * basis[1][0],
                    left * basis[0][1] + right * basis[1][1],
                ),
                finite_residue,
                finite_modulus,
                infinity_modulus,
            )
            if normalized is not None:
                answers[normalized[:2]] = normalized
    return tuple(
        sorted(
            answers.values(),
            key=lambda item: (item[2], abs(item[0]), item[1], item[0]),
        )
    )


def extend_state(
    state: BeamState,
    choice: ProjectiveChoice,
    *,
    height_weight: float,
    beam_coefficient_radius: int,
) -> BeamState | None:
    finite_residue = state.finite_residue
    finite_modulus = state.finite_modulus
    infinity_modulus = state.infinity_modulus
    if choice.projective_index == choice.prime:
        infinity_modulus *= choice.prime
    else:
        finite_residue, finite_modulus = crt_pair(
            finite_residue,
            finite_modulus,
            choice.projective_index,
            choice.prime,
        )
    representatives = short_representatives(
        finite_residue,
        finite_modulus,
        infinity_modulus,
        coefficient_radius=beam_coefficient_radius,
    )
    if not representatives:
        return None
    representative = representatives[0]
    local_score_units = state.local_score_units + choice.score_units
    objective = local_score_units / SCORE_SCALE - height_weight * log(
        max(2, representative[2])
    )
    return BeamState(
        finite_residue=finite_residue,
        finite_modulus=finite_modulus,
        infinity_modulus=infinity_modulus,
        choices=state.choices + (choice,),
        local_score_units=local_score_units,
        reduced_basis=gauss_reduce(
            *congruence_basis(finite_residue, finite_modulus, infinity_modulus)
        ),
        representative=representative,
        objective=objective,
    )


def beam_combine_projective(
    groups: Sequence[Sequence[ProjectiveChoice]],
    *,
    beam_width: int,
    height_weight: float,
    beam_coefficient_radius: int,
) -> tuple[BeamState, ...]:
    if beam_width < 1:
        raise ValueError("beam_width must be positive")
    initial = BeamState(0, 1, 1, (), 0, ((1, 0), (0, 1)), (0, 1, 1), 0.0)
    states = (initial,)
    used_primes: set[int] = set()
    for group in groups:
        primes = {choice.prime for choice in group}
        if len(primes) != 1:
            raise ValueError("each beam group must contain one prime")
        prime = next(iter(primes))
        if prime in used_primes:
            raise ValueError("a discovery prime occurs twice")
        used_primes.add(prime)
        expanded: dict[tuple[int, int, int], BeamState] = {}
        for state in states:
            for choice in group:
                candidate = extend_state(
                    state,
                    choice,
                    height_weight=height_weight,
                    beam_coefficient_radius=beam_coefficient_radius,
                )
                if candidate is None:
                    continue
                key = (
                    candidate.finite_residue,
                    candidate.finite_modulus,
                    candidate.infinity_modulus,
                )
                previous = expanded.get(key)
                if previous is None or (
                    -candidate.objective,
                    candidate.representative[2],
                    key,
                ) < (
                    -previous.objective,
                    previous.representative[2],
                    key,
                ):
                    expanded[key] = candidate
        states = tuple(
            sorted(
                expanded.values(),
                key=lambda candidate: (
                    -candidate.objective,
                    candidate.representative[2],
                    candidate.finite_residue,
                    candidate.infinity_modulus,
                ),
            )[:beam_width]
        )
        if not states:
            raise RuntimeError(f"beam became empty at p={prime}")
    return states


def enumerate_beam_parameters(
    states: Sequence[BeamState], *, coefficient_radius: int, minimum_height: int
) -> tuple[tuple[int, int, int, int], ...]:
    """Return unique (a,b,height,state-index) outside the declared old box."""

    answers: dict[tuple[int, int], tuple[int, int, int, int]] = {}
    for state_index, state in enumerate(states):
        for numerator, denominator, height in short_representatives(
            state.finite_residue,
            state.finite_modulus,
            state.infinity_modulus,
            coefficient_radius=coefficient_radius,
        ):
            if height < minimum_height:
                continue
            key = (numerator, denominator)
            candidate = (numerator, denominator, height, state_index)
            previous = answers.get(key)
            if previous is None or state_index < previous[3]:
                answers[key] = candidate
    return tuple(
        sorted(answers.values(), key=lambda item: (item[2], abs(item[0]), item[1], item[0]))
    )


def prior_pairs(path: Path) -> set[tuple[int, int]]:
    document = json.loads(path.read_text())
    return {tuple(map(int, record["projective_pair"])) for record in document["finalists"]}


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def score_population(
    parameters: Iterable[tuple[int, int, int, int]],
    construction_blocks: Sequence[dict[int, tuple[object, ...]]],
    validation_tables: dict[int, tuple[object, ...]],
) -> list[dict[str, object]]:
    construction_tables = {
        prime: table for block in construction_blocks for prime, table in block.items()
    }
    construction_cache: dict[tuple[int, int], int | None] = {}
    validation_cache: dict[tuple[int, int], int | None] = {}
    records = []
    for numerator, denominator, height, state_index in parameters:
        seed = Candidate(numerator, denominator, height)
        construction = score_block(seed, construction_tables, construction_cache)
        validation = score_block(seed, validation_tables, validation_cache)
        records.append(
            {
                "parameter": f"{numerator}/{denominator}",
                "projective_pair": [numerator, denominator],
                "projective_height": height,
                "beam_state_index": state_index,
                "construction_score_units_1e12": construction.total_score_units,
                "construction_score": construction.total_score_units / SCORE_SCALE,
                "construction_good_prime_count": construction.good_primes,
                "construction_bad_reduction_prime_count": construction.bad_primes,
                "validation_score_units_1e12": validation.total_score_units,
                "validation_score": validation.total_score_units / SCORE_SCALE,
                "validation_good_prime_count": validation.good_primes,
                "validation_bad_reduction_prime_count": validation.bad_primes,
            }
        )
    construction_order = sorted(
        records,
        key=lambda row: (
            -row["construction_score_units_1e12"],
            -row["construction_good_prime_count"],
            row["construction_bad_reduction_prime_count"],
            row["projective_height"],
            row["projective_pair"][1],
            row["projective_pair"][0],
        ),
    )
    validation_order = sorted(
        records,
        key=lambda row: (
            -row["validation_score_units_1e12"],
            -row["validation_good_prime_count"],
            row["validation_bad_reduction_prime_count"],
            row["projective_height"],
            row["projective_pair"][1],
            row["projective_pair"][0],
        ),
    )
    previous_signature = None
    competition_rank = 0
    for position, record in enumerate(construction_order, 1):
        signature = (
            record["construction_score_units_1e12"],
            record["construction_good_prime_count"],
            record["construction_bad_reduction_prime_count"],
        )
        if signature != previous_signature:
            competition_rank = position
            previous_signature = signature
        record["construction_rank"] = competition_rank
    previous_signature = None
    competition_rank = 0
    for position, record in enumerate(validation_order, 1):
        signature = (
            record["validation_score_units_1e12"],
            record["validation_good_prime_count"],
            record["validation_bad_reduction_prime_count"],
        )
        if signature != previous_signature:
            competition_rank = position
            previous_signature = signature
        record["validation_rank"] = competition_rank
        record["worst_rank"] = max(record["construction_rank"], competition_rank)
        record["rank_sum"] = record["construction_rank"] + competition_rank
    return sorted(
        records,
        key=lambda row: (
            row["worst_rank"],
            row["rank_sum"],
            row["projective_height"],
            row["projective_pair"][1],
            row["projective_pair"][0],
        ),
    )


def exact_specialization(model: object, record: dict[str, object]) -> dict[str, object]:
    """Specialize the exact rational short model and clear denominators."""

    numerator, denominator = map(int, record["projective_pair"])
    if denominator == 0:
        raise ValueError("the exact affine adapter cannot specialize infinity")
    parameter = Fraction(numerator, denominator)

    def evaluate(coefficients: Sequence[Fraction]) -> Fraction:
        value = Fraction(0)
        for coefficient in reversed(coefficients):
            value = value * parameter + coefficient
        return value

    coefficient_a = evaluate(model.a_coefficients)
    coefficient_b = evaluate(model.b_coefficients)
    discriminant_core = 4 * coefficient_a**3 + 27 * coefficient_b**2
    if discriminant_core == 0:
        raise ValueError("constructed parameter specializes to a singular curve")
    scale = lcm(coefficient_a.denominator, coefficient_b.denominator)
    integral_a = coefficient_a * scale**4
    integral_b = coefficient_b * scale**6
    if integral_a.denominator != 1 or integral_b.denominator != 1:
        raise AssertionError("denominator clearing did not produce an integral model")
    integral_a_value = integral_a.numerator
    integral_b_value = integral_b.numerator
    integral_discriminant = -16 * (
        4 * integral_a_value**3 + 27 * integral_b_value**2
    )
    if integral_discriminant == 0:
        raise AssertionError("integral specialization became singular")
    if (
        (-48 * integral_a_value) ** 3 - (-864 * integral_b_value) ** 2
        != 1728 * integral_discriminant
    ):
        raise AssertionError("integral c4/c6/discriminant identity failed")

    def rational_text(value: Fraction) -> str:
        return f"{value.numerator}/{value.denominator}"

    return {
        "parameter": record["parameter"],
        "projective_pair": record["projective_pair"],
        "source_short_model": {
            "A": rational_text(coefficient_a),
            "B": rational_text(coefficient_b),
            "discriminant": rational_text(-16 * discriminant_core),
            "nonsingular": True,
        },
        "denominator_cleared_integral_short_model": {
            "equation": "y^2=x^3+A*x+B",
            "A": str(integral_a_value),
            "B": str(integral_b_value),
            "c4": str(-48 * integral_a_value),
            "c6": str(-864 * integral_b_value),
            "discriminant": str(integral_discriminant),
            "clearing_scale": str(scale),
            "isomorphism": "x_integral=scale^2*x_source; y_integral=scale^3*y_source",
            "nonsingular": True,
            "minimality_claimed": False,
        },
    }


def primes_in_interval(lower: int, upper: int) -> tuple[int, ...]:
    if lower < 5 or upper < lower:
        raise ValueError("invalid prime interval")
    return tuple(prime for prime in range(lower, upper + 1) if is_prime(prime))


def state_record(state: BeamState) -> dict[str, object]:
    return {
        "finite_residue": str(state.finite_residue),
        "finite_modulus": str(state.finite_modulus),
        "infinity_modulus": str(state.infinity_modulus),
        "total_modulus": str(state.modulus),
        "local_score_units_1e12": state.local_score_units,
        "objective": state.objective,
        "gauss_reduced_basis": [list(vector) for vector in state.reduced_basis],
        "shortest_valid_representative": list(state.representative),
        "choices": [
            {
                "prime": choice.prime,
                "projective_index": choice.projective_index,
                "label": choice.label,
                "a_p": choice.trace,
                "score_units_1e12": choice.score_units,
            }
            for choice in state.choices
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Bounded projective CRT/Gauss Nagao constructor for q12/orbit5867."
    )
    parser.add_argument("--prior", type=Path, default=DEFAULT_PRIOR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--choices-per-prime", type=int, default=6)
    parser.add_argument("--beam-width", type=int, default=512)
    parser.add_argument("--height-weight", type=float, default=0.025)
    parser.add_argument("--beam-coefficient-radius", type=int, default=2)
    parser.add_argument("--coefficient-radius", type=int, default=10)
    parser.add_argument("--minimum-height", type=int, default=10001)
    parser.add_argument("--finalists", type=int, default=1000)
    parser.add_argument("--extra-construction-prime-min", type=int)
    parser.add_argument("--extra-construction-prime-max", type=int)
    parser.add_argument("--validation-prime-min", type=int, default=199)
    parser.add_argument("--validation-prime-max", type=int, default=499)
    parser.add_argument("--specialize-top", type=int, default=0)
    args = parser.parse_args()
    if args.finalists < 1 or args.minimum_height < 1 or args.specialize_top < 0:
        raise SystemExit("finalist count and minimum height must be positive")
    if (args.extra_construction_prime_min is None) != (
        args.extra_construction_prime_max is None
    ):
        raise SystemExit("both extra construction prime bounds are required together")

    started = perf_counter()
    model = load_family_model()
    base_blocks, base_rejected = build_residue_tables(
        model, DEFAULT_PRIME_BLOCKS
    )
    extra_requested: tuple[int, ...] = ()
    extra_blocks: tuple[dict[int, tuple[object, ...]], ...] = ()
    extra_rejected: tuple[dict[str, object], ...] = ()
    if args.extra_construction_prime_min is not None:
        extra_requested = primes_in_interval(
            args.extra_construction_prime_min,
            args.extra_construction_prime_max,
        )
        extra_blocks, extra_rejected = build_residue_tables(model, (extra_requested,))
    construction_blocks = base_blocks + extra_blocks
    validation_requested = primes_in_interval(
        args.validation_prime_min, args.validation_prime_max
    )
    construction_primes = {
        prime for block in construction_blocks for prime in block
    }
    prime_overlap = construction_primes.intersection(validation_requested)
    if prime_overlap:
        raise ValueError(
            f"validation interval overlaps usable construction primes: {sorted(prime_overlap)}"
        )
    validation_blocks, validation_rejected = build_residue_tables(
        model, (validation_requested,)
    )
    groups = projective_choice_groups(construction_blocks, args.choices_per_prime)
    beam_started = perf_counter()
    states = beam_combine_projective(
        groups,
        beam_width=args.beam_width,
        height_weight=args.height_weight,
        beam_coefficient_radius=args.beam_coefficient_radius,
    )
    beam_seconds = perf_counter() - beam_started
    parameters = enumerate_beam_parameters(
        states,
        coefficient_radius=args.coefficient_radius,
        minimum_height=args.minimum_height,
    )
    old_pairs = prior_pairs(args.prior)
    novel_parameters = tuple(
        parameter for parameter in parameters if parameter[:2] not in old_pairs
    )
    population_text = "".join(
        f"{numerator}/{denominator}\n"
        for numerator, denominator, _, _ in novel_parameters
    ).encode()
    scored = score_population(
        novel_parameters, construction_blocks, validation_blocks[0]
    )
    robust_finalists = scored[: args.finalists]
    construction_extreme = sorted(
        scored,
        key=lambda row: (
            row["construction_rank"],
            row["validation_rank"],
            row["projective_height"],
        ),
    )[: args.finalists]
    validation_extreme = sorted(
        scored,
        key=lambda row: (
            row["validation_rank"],
            row["construction_rank"],
            row["projective_height"],
        ),
    )[: args.finalists]
    exact_specializations = [
        exact_specialization(model, record)
        for record in robust_finalists[: args.specialize_top]
    ]
    runtime = perf_counter() - started

    output = {
        "schema": "h92-q12o5867-rootless-projective-crt-gauss-nagao-v2",
        "status": "PASS_BOUNDED_HEURISTIC_PROJECTIVE_CRT_GAUSS_CONSTRUCTOR",
        "proof_boundary": (
            "This is a bounded Nagao-score constructor and disjoint-prime rerank. "
            "Scores, beam survival, and CRT profiles are not rank evidence. Exact "
            "specialization, when requested, does not evaluate sections or search points."
        ),
        "model_sha256": model.source_sha256,
        "prior_h10000": {
            "path": str(args.prior.resolve()),
            "sha256": file_sha256(args.prior),
        },
        "bounds": {
            "choices_per_prime": args.choices_per_prime,
            "beam_width": args.beam_width,
            "height_weight": args.height_weight,
            "beam_coefficient_radius": args.beam_coefficient_radius,
            "final_coefficient_radius": args.coefficient_radius,
            "minimum_height": args.minimum_height,
            "base_construction_primes": [
                prime for block in base_blocks for prime in block
            ],
            "extra_construction_requested_primes": list(extra_requested),
            "extra_construction_usable_primes": [
                prime for block in extra_blocks for prime in block
            ],
            "all_construction_primes": [choice_group[0].prime for choice_group in groups],
            "validation_requested_primes": list(validation_requested),
            "validation_usable_primes": list(validation_blocks[0]),
            "base_construction_rejected_primes": list(base_rejected),
            "extra_construction_rejected_primes": list(extra_rejected),
            "validation_rejected_primes": list(validation_rejected),
            "construction_validation_disjoint": True,
        },
        "construction": {
            "congruences": (
                "finite choices impose a=rb mod p; infinity choices impose b=0 mod p"
            ),
            "beam_state_count": len(states),
            "beam_states": [state_record(state) for state in states],
            "enumerated_parameter_count_before_prior_dedup": len(parameters),
            "prior_h10000_pair_count": len(old_pairs),
            "novel_parameter_count": len(novel_parameters),
            "novel_parameter_population_sha256": sha256(population_text).hexdigest(),
            "minimum_observed_height": min(parameter[2] for parameter in novel_parameters),
            "maximum_observed_height": max(parameter[2] for parameter in novel_parameters),
        },
        "ranking": {
            "method": (
                "competition ranks tie equal score/good/bad signatures; minimize "
                "max(construction_rank,validation_rank), then rank sum"
            ),
            "scored_population_count": len(scored),
            "stored_per_order": args.finalists,
            "robust_finalists": robust_finalists,
            "construction_extreme": construction_extreme,
            "validation_extreme": validation_extreme,
        },
        "exact_specializations": exact_specializations,
        "point_search_launched": False,
        "runtime": {
            "beam_seconds": beam_seconds,
            "total_seconds": runtime,
        },
        "reproducing_command": shlex.join(sys.argv),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        f"PASS states={len(states)} novel={len(novel_parameters)} "
        f"top_robust={robust_finalists[0]['parameter']} "
        f"top_validation={validation_extreme[0]['parameter']} "
        f"specialized={len(exact_specializations)} "
        f"seconds={runtime:.3f} output={args.output}"
    )


if __name__ == "__main__":
    main()
