#!/usr/bin/env python3
"""Search higher-weight mod-2 classes on the rank-29 record curve.

The earlier alternate-cover pass exhausts positive subset sums of Hamming
weight two and three.  This pass is disjoint from it: every class mask has
weight at least four.  It also uses a freedom absent from that pass.  Replacing
``P_i`` by ``-P_i`` changes a representative by ``2*P_i`` and hence preserves
its class in ``E(Q)/2E(Q)``.  A deterministic beam/local search optimizes these
signed representatives and their class masks simultaneously.

The objective is exact and chart-aware.  For a base point ``Q`` it forms the
58 alternate parameters of ``+/- P_i``, then scores the best affine
normalizations in a small declared pool.  Five Hamming-weight bands are
searched separately so that coordinate growth does not make the retained
tranche collapse to weight four.  Global negation is canonicalized; bitwise
complements are *not* identified because they are generally different mod-2
classes.

Each retained cover receives bounded offset, affine cross-ratio, three-point
Mobius, and skew-denominator searches.  PARI calls are one-shot foreground
process groups with strict timeouts and no retries.  Every returned affine
point is mapped and checked over ``QQ``.  Known seed/companion relations and
PARI-proposed relations are replayed with exact Fraction group arithmetic;
anything left unresolved immediately receives an exact finite-reduction
rank-30 test.

This is a bounded point search, not a rank upper bound.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import hashlib
from itertools import combinations
import json
from pathlib import Path
import platform
import random
import re
import shlex
import sys
import time
from typing import Any, Iterable, Sequence

from alternate_quartic_covers import (
    AlternateQuarticCover,
    alternate_cover,
    mobius_preimage,
    short_add,
    three_point_mobius_matrix,
)
from elkies_klagsbrun_rank29 import (
    GENERAL_WEIERSTRASS_COEFFICIENTS,
    PUBLISHED_POINTS,
    from_short_point,
    point_on_general_curve,
    point_on_short_curve,
    published_short_points,
    short_weierstrass_coefficients,
    to_short_point,
)
from mod2_reduction_independence import (
    Mod2ReductionSignature,
    combined_mod2_rank,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)
from pari_bridge import pari_version
from search_elkies_klagsbrun_rank30 import (
    RationalPoint,
    affine_substitute,
    exact_linear_combination,
    point_add,
    point_negate,
    poly_evaluate,
)
from search_elkies_klagsbrun_rank30_alternate_covers import (
    alternate_parameter,
    parameter_height,
)
from search_extra_points import gp_rational, gp_vector
from search_nagao_rank21_t956_skew import run_gp_once, search_original_quartic
from search_nagao_u42_skew_height import map_chart_point, transform_binary_quartic


Q = Fraction
POINT_COUNT = len(PUBLISHED_POINTS)
ALL_MASK = (1 << POINT_COUNT) - 1
SHORT_PUBLIC_POINTS = published_short_points()
SIGNED_SHORT_PUBLIC_POINTS = tuple(
    oriented
    for point in SHORT_PUBLIC_POINTS
    for oriented in (point, (point[0], -point[1]))
)
DEFAULT_SEED = 20260814
WEIGHT_BANDS = ((4, 6), (7, 10), (11, 17), (18, 25), (26, 29))
DEFAULT_OUTPUT = Path(
    "artifacts/generated-results/"
    "elliptic_elkies_klagsbrun_rank30_higher_weight_covers.json"
)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/"
    "search_elkies_klagsbrun_rank30_higher_weight_covers.py"
)


@dataclass(frozen=True)
class SignedClassState:
    """A signed representative of one 29-bit mod-2 class.

    A set bit has coefficient ``+1`` or ``-1`` according to
    ``negative_mask``.  The least-index nonzero coefficient is always positive,
    which identifies a representative with its global negative.
    """

    mask: int
    negative_mask: int

    def __post_init__(self) -> None:
        if not 0 < self.mask <= ALL_MASK:
            raise ValueError("a class mask must be a nonzero 29-bit integer")
        if self.negative_mask & ~self.mask:
            raise ValueError("negative signs may occur only on selected indices")
        least_bit = self.mask & -self.mask
        if self.negative_mask & least_bit:
            raise ValueError("global sign is not canonical")

    @property
    def weight(self) -> int:
        return self.mask.bit_count()

    @property
    def identifier(self) -> str:
        return (
            f"w{self.weight:02d}_m{self.mask:08x}_"
            f"n{self.negative_mask:08x}"
        )

    @property
    def coefficients(self) -> tuple[int, ...]:
        return tuple(
            0
            if not (self.mask >> index) & 1
            else (-1 if (self.negative_mask >> index) & 1 else 1)
            for index in range(POINT_COUNT)
        )


@dataclass(frozen=True)
class Evaluation:
    state: SignedClassState
    score: tuple[int, ...]
    best_anchor_indices: tuple[int, int]
    base_coordinate_bits: int


@dataclass(frozen=True)
class CoverCandidate:
    evaluation: Evaluation
    base_short_point: RationalPoint
    parameters: tuple[Fraction, ...]

    @property
    def state(self) -> SignedClassState:
        return self.evaluation.state

    @property
    def cover(self) -> AlternateQuarticCover:
        return alternate_cover(short_weierstrass_coefficients(), self.base_short_point)


@dataclass(frozen=True)
class CoverChartPlan:
    identifier: str
    kind: str
    candidate_index: int
    polynomial: tuple[Fraction, ...]
    height_specification: str
    center: Fraction | None = None
    scale: Fraction | None = None
    matrix: tuple[int, int, int, int] | None = None
    compression_score: tuple[int, ...] = ()


def canonical_state(mask: int, negative_mask: int) -> SignedClassState:
    """Return the unique global-sign normalization of signed coefficients."""

    mask = int(mask)
    negative_mask = int(negative_mask) & mask
    if not 0 < mask <= ALL_MASK:
        raise ValueError("a class mask must be a nonzero 29-bit integer")
    least_bit = mask & -mask
    if negative_mask & least_bit:
        negative_mask ^= mask
    return SignedClassState(mask, negative_mask)


def state_indices(state: SignedClassState) -> tuple[int, ...]:
    return tuple(index for index in range(POINT_COUNT) if (state.mask >> index) & 1)


def signed_short_sum(state: SignedClassState) -> RationalPoint:
    """Construct the representative exactly on the integral short model."""

    coefficients = short_weierstrass_coefficients()
    answer = None
    for index, coefficient in enumerate(state.coefficients):
        if coefficient == 0:
            continue
        point = SHORT_PUBLIC_POINTS[index]
        if coefficient < 0:
            point = point[0], -point[1]
        answer = short_add(coefficients, answer, point)
    if answer is None:
        raise AssertionError("a nonzero class representative vanished")
    if not point_on_short_curve(answer):
        raise AssertionError("a signed subset sum left the short curve")
    return answer


def signed_public_short_points() -> tuple[RationalPoint, ...]:
    return SIGNED_SHORT_PUBLIC_POINTS


def signed_parameter_label(index: int) -> str:
    point_index, sign_index = divmod(index, 2)
    return f"p{point_index + 1:02d}{'+' if sign_index == 0 else '-'}"


def cover_parameters(base_point: RationalPoint) -> tuple[Fraction, ...]:
    answer = []
    for point in signed_public_short_points():
        try:
            answer.append(alternate_parameter(base_point, point))
        except ValueError as error:
            raise AssertionError(
                "a higher-weight representative collided with a public point"
            ) from error
    if len(answer) != 2 * POINT_COUNT:
        raise AssertionError("the signed public parameter count changed")
    return tuple(answer)


def rational_bit_height(value: Fraction) -> int:
    return parameter_height(Q(value)).bit_length()


def height_signature(values: Iterable[Fraction]) -> tuple[int, ...]:
    heights = sorted(rational_bit_height(value) for value in values)
    if len(heights) < 20:
        raise ValueError("a compression score needs at least twenty values")
    positions = (4, 9, len(heights) // 3, len(heights) // 2, 5 * len(heights) // 6, len(heights) - 1)
    return tuple(heights[index] for index in positions) + (
        sum(heights[2 : min(16, len(heights))]),
    )


def best_affine_normalization(
    parameters: Sequence[Fraction], *, pool_size: int
) -> tuple[tuple[int, ...], tuple[int, int]]:
    """Score the best ``t=center+scale*s`` from a declared small pool."""

    if not 2 <= pool_size <= len(parameters):
        raise ValueError("the affine normalization pool size is invalid")
    pool = tuple(
        sorted(
            range(len(parameters)),
            key=lambda index: (parameter_height(parameters[index]), index),
        )[:pool_size]
    )
    best: tuple[tuple[int, ...], tuple[int, int]] | None = None
    for first, second in combinations(pool, 2):
        center = Q(parameters[first])
        scale = Q(parameters[second]) - center
        if scale == 0:
            continue
        normalized = tuple((Q(value) - center) / scale for value in parameters)
        candidate = height_signature(normalized), (first, second)
        if best is None or candidate < best:
            best = candidate
    if best is None:
        raise AssertionError("all affine normalization anchors coincided")
    return best


def evaluate_state(state: SignedClassState, *, pool_size: int = 5) -> Evaluation:
    if state.weight <= 3:
        raise ValueError("weight-zero through weight-three classes are excluded")
    base = signed_short_sum(state)
    parameters = cover_parameters(base)
    score, anchors = best_affine_normalization(parameters, pool_size=pool_size)
    coordinate_bits = max(
        abs(base[0].numerator).bit_length(),
        base[0].denominator.bit_length(),
        abs(base[1].numerator).bit_length(),
        base[1].denominator.bit_length(),
    )
    return Evaluation(state, score + (coordinate_bits,), anchors, coordinate_bits)


def build_candidate(evaluation: Evaluation) -> CoverCandidate:
    base = signed_short_sum(evaluation.state)
    parameters = cover_parameters(base)
    replay = evaluate_state(evaluation.state)
    if replay != evaluation:
        raise AssertionError("a selected beam evaluation was not deterministic")
    return CoverCandidate(evaluation, base, parameters)


def random_state(
    rng: random.Random, lower_weight: int, upper_weight: int
) -> SignedClassState:
    weight = rng.randint(lower_weight, upper_weight)
    indices = rng.sample(range(POINT_COUNT), weight)
    mask = sum(1 << index for index in indices)
    negative = sum(1 << index for index in indices if rng.getrandbits(1))
    return canonical_state(mask, negative)


def mutate_state(
    state: SignedClassState,
    rng: random.Random,
    lower_weight: int,
    upper_weight: int,
) -> SignedClassState:
    """Apply one deterministic-RNG local move within a weight band."""

    selected = list(state_indices(state))
    unselected = [index for index in range(POINT_COUNT) if not (state.mask >> index) & 1]
    move = rng.randrange(5)
    mask = state.mask
    negative = state.negative_mask
    if move == 0 and len(selected) >= 2:
        # Sign motion stays in the same mod-2 class.
        index = rng.choice(selected)
        negative ^= 1 << index
    elif move == 1 and len(selected) >= 3:
        for index in rng.sample(selected, 2):
            negative ^= 1 << index
    elif move == 2 and selected and unselected:
        removed = rng.choice(selected)
        added = rng.choice(unselected)
        mask ^= (1 << removed) | (1 << added)
        negative &= ~(1 << removed)
        if rng.getrandbits(1):
            negative |= 1 << added
    elif move == 3 and state.weight < upper_weight and unselected:
        added = rng.choice(unselected)
        mask |= 1 << added
        if rng.getrandbits(1):
            negative |= 1 << added
    elif state.weight > lower_weight and selected:
        removed = rng.choice(selected)
        mask &= ~(1 << removed)
        negative &= ~(1 << removed)
    elif selected:
        index = rng.choice(selected)
        negative ^= 1 << index
    answer = canonical_state(mask, negative)
    if not lower_weight <= answer.weight <= upper_weight:
        raise AssertionError("a mutation escaped its declared weight band")
    return answer


def _beam(
    evaluations: Sequence[Evaluation], *, lower: int, upper: int, width: int
) -> tuple[Evaluation, ...]:
    best_by_mask: dict[int, Evaluation] = {}
    for evaluation in evaluations:
        if not lower <= evaluation.state.weight <= upper:
            continue
        prior = best_by_mask.get(evaluation.state.mask)
        if prior is None or (evaluation.score, evaluation.state.identifier) < (
            prior.score,
            prior.state.identifier,
        ):
            best_by_mask[evaluation.state.mask] = evaluation
    return tuple(
        sorted(
            best_by_mask.values(),
            key=lambda item: (item.score, item.state.identifier),
        )[:width]
    )


def search_signed_classes(
    *,
    seed: int,
    evaluation_budget: int,
    beam_width: int,
    rounds: int,
    mutations_per_state: int,
) -> tuple[tuple[Evaluation, ...], tuple[dict[str, Any], ...]]:
    """Run the pinned stratified beam/local search exactly once."""

    if evaluation_budget < len(WEIGHT_BANDS) * 20:
        raise ValueError("the evaluation budget is too small for five bands")
    if min(beam_width, rounds, mutations_per_state) <= 0:
        raise ValueError("beam parameters must be positive")
    rng = random.Random(seed)
    base_budget, remainder = divmod(evaluation_budget, len(WEIGHT_BANDS))
    all_evaluations: list[Evaluation] = []
    band_records: list[dict[str, Any]] = []
    for band_index, (lower, upper) in enumerate(WEIGHT_BANDS):
        budget = base_budget + (band_index < remainder)
        evaluations: dict[SignedClassState, Evaluation] = {}

        def add(state: SignedClassState) -> bool:
            if state in evaluations or len(evaluations) >= budget:
                return False
            evaluations[state] = evaluate_state(state)
            return True

        # Cyclic windows provide deterministic structure; random signs stop
        # them from being only the positive subset representatives.
        for weight in sorted({lower, (lower + upper) // 2, upper}):
            for start in range(POINT_COUNT):
                indices = tuple((start + offset) % POINT_COUNT for offset in range(weight))
                mask = sum(1 << index for index in indices)
                negative = sum(1 << index for index in indices if rng.getrandbits(1))
                add(canonical_state(mask, negative))
                if len(evaluations) >= min(budget, 2 * beam_width):
                    break
            if len(evaluations) >= min(budget, 2 * beam_width):
                break
        while len(evaluations) < min(budget, 2 * beam_width):
            add(random_state(rng, lower, upper))

        for _ in range(rounds):
            current = _beam(
                tuple(evaluations.values()), lower=lower, upper=upper, width=beam_width
            )
            for parent in current:
                for _ in range(mutations_per_state):
                    add(mutate_state(parent.state, rng, lower, upper))
                    if len(evaluations) >= budget:
                        break
                if len(evaluations) >= budget:
                    break
            if len(evaluations) >= budget:
                break
        # Use random injections to make the exact declared budget independent
        # of duplicate local proposals.
        while len(evaluations) < budget:
            add(random_state(rng, lower, upper))

        ordered = tuple(
            sorted(
                evaluations.values(),
                key=lambda item: (item.score, item.state.identifier),
            )
        )
        all_evaluations.extend(ordered)
        band_records.append(
            {
                "weight_lower": lower,
                "weight_upper": upper,
                "evaluated_signed_representative_count": len(ordered),
                "distinct_class_mask_count": len({item.state.mask for item in ordered}),
                "best_state": ordered[0].state.identifier,
                "best_score": list(ordered[0].score),
            }
        )
        print(
            f"beam band {lower:02d}-{upper:02d}: "
            f"evaluated={len(ordered)} best={ordered[0].score}",
            flush=True,
        )
    return tuple(all_evaluations), tuple(band_records)


def hamming_distance(left: int, right: int) -> int:
    return (int(left) ^ int(right)).bit_count()


def select_diverse_evaluations(
    evaluations: Sequence[Evaluation], *, count: int, minimum_distance: int
) -> tuple[Evaluation, ...]:
    """Retain a weight-stratified, class-mask-diverse top tranche."""

    if count < len(WEIGHT_BANDS):
        raise ValueError("at least one retained cover per weight band is required")
    if minimum_distance < 0:
        raise ValueError("the Hamming distance must be nonnegative")
    best_by_mask: dict[int, Evaluation] = {}
    for evaluation in evaluations:
        prior = best_by_mask.get(evaluation.state.mask)
        if prior is None or (evaluation.score, evaluation.state.identifier) < (
            prior.score,
            prior.state.identifier,
        ):
            best_by_mask[evaluation.state.mask] = evaluation
    ordered = tuple(
        sorted(
            best_by_mask.values(),
            key=lambda item: (item.score, item.state.identifier),
        )
    )
    selected: list[Evaluation] = []

    def eligible(evaluation: Evaluation, distance: int) -> bool:
        return evaluation not in selected and all(
            hamming_distance(evaluation.state.mask, prior.state.mask) >= distance
            for prior in selected
        )

    quota, extra = divmod(count, len(WEIGHT_BANDS))
    for band_index, (lower, upper) in enumerate(WEIGHT_BANDS):
        band_quota = quota + (band_index < extra)
        band = tuple(
            item for item in ordered if lower <= item.state.weight <= upper
        )
        for distance in (minimum_distance, max(0, minimum_distance // 2), 0):
            for evaluation in band:
                if sum(lower <= item.state.weight <= upper for item in selected) >= band_quota:
                    break
                if eligible(evaluation, distance):
                    selected.append(evaluation)
            if sum(lower <= item.state.weight <= upper for item in selected) >= band_quota:
                break
    for distance in (minimum_distance, max(0, minimum_distance // 2), 0):
        for evaluation in ordered:
            if len(selected) >= count:
                break
            if eligible(evaluation, distance):
                selected.append(evaluation)
        if len(selected) >= count:
            break
    if len(selected) != count:
        raise AssertionError("the diverse selector could not fill its tranche")
    return tuple(selected)


def coefficient_bit_size(polynomial: Sequence[Fraction]) -> int:
    return max(
        max(abs(Q(value).numerator).bit_length(), Q(value).denominator.bit_length())
        for value in polynomial
    )


def optimized_affine_pairs(
    candidate: CoverCandidate, *, pool_size: int, count: int
) -> tuple[tuple[tuple[int, ...], int, int], ...]:
    parameters = candidate.parameters
    pool = tuple(
        sorted(
            range(len(parameters)),
            key=lambda index: (parameter_height(parameters[index]), index),
        )[:pool_size]
    )
    ranked = []
    for first, second in combinations(pool, 2):
        center = parameters[first]
        scale = parameters[second] - center
        if scale == 0:
            continue
        normalized = tuple((value - center) / scale for value in parameters)
        polynomial = affine_substitute(candidate.cover.coefficients, center, scale)
        score = height_signature(normalized) + (coefficient_bit_size(polynomial),)
        ranked.append((score, first, second))
    ranked.sort(key=lambda item: (item[0], item[1], item[2]))
    return tuple(ranked[:count])


def optimized_mobius_charts(
    candidate: CoverCandidate,
    *,
    affine_pairs: Sequence[tuple[tuple[int, ...], int, int]],
    pool_size: int,
    count: int,
) -> tuple[tuple[tuple[int, ...], tuple[int, int, int, int], tuple[int, int, int]], ...]:
    parameters = candidate.parameters
    raw_pool = sorted(
        range(len(parameters)),
        key=lambda index: (parameter_height(parameters[index]), index),
    )
    preferred = []
    for _, first, second in affine_pairs:
        preferred.extend((first, second))
    preferred.extend(raw_pool)
    pool = []
    for index in preferred:
        if index not in pool:
            pool.append(index)
        if len(pool) >= pool_size:
            break
    ranked = []
    seen: set[tuple[int, int, int, int]] = set()
    for indices in combinations(pool, 3):
        values = tuple(parameters[index] for index in indices)
        for rotation in range(3):
            rotated_values = values[rotation:] + values[:rotation]
            rotated_indices = indices[rotation:] + indices[:rotation]
            if len(set(rotated_values)) != 3:
                continue
            matrix = three_point_mobius_matrix(*rotated_values)
            if matrix in seen:
                continue
            seen.add(matrix)
            preimages = tuple(
                preimage
                for value in parameters
                if (preimage := mobius_preimage(matrix, value)) is not None
            )
            transformed = transform_binary_quartic(candidate.cover.coefficients, matrix)
            score = height_signature(preimages) + (coefficient_bit_size(transformed),)
            ranked.append((score, matrix, rotated_indices))
    ranked.sort(key=lambda item: (item[0], item[1], item[2]))
    return tuple(ranked[:count])


def build_cover_chart_plans(
    candidates: Sequence[CoverCandidate],
    *,
    offset_count: int,
    affine_count: int,
    mobius_count: int,
    affine_pool_size: int,
    mobius_pool_size: int,
    offset_height: int,
    affine_height: int,
    mobius_height: int,
    skew_numerator_bound: int,
    skew_denominator_bound: int,
) -> tuple[CoverChartPlan, ...]:
    plans: list[CoverChartPlan] = []
    for candidate_index, candidate in enumerate(candidates):
        cover = candidate.cover
        offset_ranked = []
        for parameter_index in sorted(
            range(len(candidate.parameters)),
            key=lambda index: (parameter_height(candidate.parameters[index]), index),
        ):
            center = candidate.parameters[parameter_index]
            polynomial = affine_substitute(cover.coefficients, center, Q(1))
            offset_ranked.append(
                (coefficient_bit_size(polynomial), parameter_index, polynomial)
            )
        offset_ranked.sort(key=lambda item: (item[0], item[1]))
        for rank, (bits, parameter_index, polynomial) in enumerate(
            offset_ranked[:offset_count], 1
        ):
            plans.append(
                CoverChartPlan(
                    identifier=(
                        f"c{candidate_index + 1:02d}_offset{rank}_"
                        f"{signed_parameter_label(parameter_index)}"
                    ),
                    kind="integer_offset",
                    candidate_index=candidate_index,
                    polynomial=polynomial,
                    height_specification=str(offset_height),
                    center=candidate.parameters[parameter_index],
                    scale=Q(1),
                    compression_score=(bits,),
                )
            )

        affine_pairs = optimized_affine_pairs(
            candidate, pool_size=affine_pool_size, count=affine_count
        )
        for rank, (score, first, second) in enumerate(affine_pairs, 1):
            center = candidate.parameters[first]
            scale = candidate.parameters[second] - center
            plans.append(
                CoverChartPlan(
                    identifier=(
                        f"c{candidate_index + 1:02d}_affine{rank}_"
                        f"{signed_parameter_label(first)}_"
                        f"{signed_parameter_label(second)}"
                    ),
                    kind="two_point_affine_cross_ratio",
                    candidate_index=candidate_index,
                    polynomial=affine_substitute(cover.coefficients, center, scale),
                    height_specification=str(affine_height),
                    center=center,
                    scale=scale,
                    compression_score=score,
                )
            )

        mobius = optimized_mobius_charts(
            candidate,
            affine_pairs=affine_pairs,
            pool_size=mobius_pool_size,
            count=mobius_count,
        )
        for rank, (score, matrix, indices) in enumerate(mobius, 1):
            transformed = transform_binary_quartic(cover.coefficients, matrix)
            label = "_".join(signed_parameter_label(index) for index in indices)
            plans.append(
                CoverChartPlan(
                    identifier=f"c{candidate_index + 1:02d}_mobius{rank}_{label}",
                    kind="three_point_mobius_cross_ratio",
                    candidate_index=candidate_index,
                    polynomial=transformed,
                    height_specification=str(mobius_height),
                    matrix=matrix,
                    compression_score=score,
                )
            )
            plans.append(
                CoverChartPlan(
                    identifier=f"c{candidate_index + 1:02d}_skew{rank}_{label}",
                    kind="mobius_skew_denominator_box",
                    candidate_index=candidate_index,
                    polynomial=transformed,
                    height_specification=(
                        f"[{skew_numerator_bound},{skew_denominator_bound}]"
                    ),
                    matrix=matrix,
                    compression_score=score,
                )
            )
    return tuple(plans)


def map_plan_points(
    plan: CoverChartPlan,
    candidate: CoverCandidate,
    raw_points: Iterable[tuple[Fraction, Fraction]],
) -> tuple[tuple[RationalPoint, ...], int]:
    """Map every finite returned point exactly; report Mobius pole count."""

    images = []
    pole_count = 0
    cover = candidate.cover
    for raw_parameter, raw_ordinate in raw_points:
        raw = Q(raw_parameter), Q(raw_ordinate)
        if poly_evaluate(plan.polynomial, raw[0]) != raw[1] ** 2:
            raise AssertionError("PARI returned a point off a transformed quartic")
        if plan.matrix is not None:
            original = map_chart_point(raw, plan.matrix)
            if original is None:
                pole_count += 1
                continue
        else:
            if plan.center is None or plan.scale is None:
                raise AssertionError("an affine plan omitted its center or scale")
            original = plan.center + plan.scale * raw[0], raw[1]
        if cover.value(original[0]) != original[1] ** 2:
            raise AssertionError("a chart map missed the original alternate quartic")
        short_point = cover.cover_point_to_curve(original)
        if not point_on_short_curve(short_point):
            raise AssertionError("an alternate-cover image missed the short curve")
        point = from_short_point(short_point)
        if not point_on_general_curve(point):
            raise AssertionError("short-model inverse transport missed the curve")
        images.append(point)
    return tuple(images), pole_count


def state_relation(state: SignedClassState) -> tuple[int, ...]:
    return state.coefficients


def known_relation_map(candidates: Sequence[CoverCandidate]) -> dict[RationalPoint, tuple[int, ...]]:
    """Relations for +/- public seeds and their ``Q-R`` companions."""

    answer: dict[RationalPoint, tuple[int, ...]] = {}
    for candidate in candidates:
        base_general = from_short_point(candidate.base_short_point)
        base_relation = state_relation(candidate.state)
        if exact_linear_combination(base_relation) != base_general:
            raise AssertionError("a signed base relation failed exact replay")
        for index, point in enumerate(PUBLISHED_POINTS):
            for sign in (1, -1):
                oriented = point if sign > 0 else point_negate(point)
                relation = [0] * POINT_COUNT
                relation[index] = sign
                relation_tuple = tuple(relation)
                answer.setdefault(oriented, relation_tuple)
                companion = point_add(base_general, point_negate(oriented))
                if companion is None:
                    continue
                companion_relation = tuple(
                    left - right
                    for left, right in zip(base_relation, relation_tuple)
                )
                answer.setdefault(companion, companion_relation)
    return answer


def discover_relations_batch(
    points: Sequence[RationalPoint], *, timeout: float, stack_bytes: int
) -> tuple[tuple[tuple[int, ...] | None, ...], dict[str, Any]]:
    """Propose several relations in one height-matrix process, then replay."""

    if not points:
        return (), {"status": "not_needed", "retried": False, "wall_seconds": 0.0}
    curve = ",".join(gp_rational(value) for value in GENERAL_WEIERSTRASS_COEFFICIENTS)
    basis = ",".join(gp_vector(value) for value in PUBLISHED_POINTS)
    commands = [
        "default(realprecision,160);",
        f"E=ellinit([{curve}]);",
        f"B=[{basis}];",
        "H=ellheightmatrix(E,B);",
    ]
    for index, point in enumerate(points):
        commands.extend(
            (
                f"Q={gp_vector(point)};",
                "V=vector(#B,j,ellheight(E,B[j],Q))~;",
                "C=round(matsolve(H,V));",
                "S=[0];for(j=1,#B,S=elladd(E,S,ellmul(E,B[j],C[j])));",
                f'print("RELATION_{index} ",Vec(C)," EXACT ",S==Q);',
            )
        )
    commands.append("quit")
    output, process = run_gp_once(
        "\n".join(commands) + "\n", timeout=timeout, stack_bytes=stack_bytes
    )
    process = {**process, "timeout_seconds": timeout, "retried": False}
    if output is None:
        return tuple(None for _ in points), process
    relations: list[tuple[int, ...] | None] = []
    for index, point in enumerate(points):
        match = re.search(
            rf"^RELATION_{index} \[(.*?)\] EXACT ([01])$", output, re.MULTILINE
        )
        if match is None or match.group(2) != "1":
            relations.append(None)
            continue
        relation = tuple(int(value.strip()) for value in match.group(1).split(","))
        if len(relation) != POINT_COUNT:
            raise AssertionError("PARI returned a relation of the wrong length")
        if exact_linear_combination(relation) != point:
            raise AssertionError("a proposed relation failed exact Fraction replay")
        relations.append(relation)
    return tuple(relations), process


def signature_record(signature: Mod2ReductionSignature) -> dict[str, Any]:
    return {
        "prime": signature.prime,
        "group_order": signature.group_order,
        "doubled_subgroup_order": signature.doubled_subgroup_order,
        "quotient_dimension": signature.quotient_dimension,
        "rows": [list(row) for row in signature.rows],
    }


def point_record(point: RationalPoint) -> dict[str, str]:
    return {"x": str(point[0]), "y": str(point[1])}


def evaluation_payload(evaluation: Evaluation) -> dict[str, Any]:
    state = evaluation.state
    return {
        "class_mask_hex": f"0x{state.mask:08x}",
        "negative_mask_hex": f"0x{state.negative_mask:08x}",
        "weight": state.weight,
        "score": list(evaluation.score),
        "best_anchor_indices": list(evaluation.best_anchor_indices),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--evaluation-budget", type=int, default=2000)
    parser.add_argument("--beam-width", type=int, default=20)
    parser.add_argument("--beam-rounds", type=int, default=5)
    parser.add_argument("--mutations-per-state", type=int, default=4)
    parser.add_argument("--cover-count", type=int, default=10)
    parser.add_argument("--minimum-mask-distance", type=int, default=4)
    parser.add_argument("--offset-count", type=int, default=1)
    parser.add_argument("--affine-count", type=int, default=2)
    parser.add_argument("--mobius-count", type=int, default=1)
    parser.add_argument("--affine-pool-size", type=int, default=12)
    parser.add_argument("--mobius-pool-size", type=int, default=9)
    parser.add_argument("--offset-height", type=int, default=100_000)
    parser.add_argument("--affine-height", type=int, default=50_000)
    parser.add_argument("--mobius-height", type=int, default=20_000)
    parser.add_argument("--skew-numerator-bound", type=int, default=100_000_000)
    parser.add_argument("--skew-denominator-bound", type=int, default=50)
    parser.add_argument("--chart-timeout", type=float, default=4.0)
    parser.add_argument("--relation-timeout", type=float, default=60.0)
    parser.add_argument("--relation-batch-size", type=int, default=16)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument("--certificate-prime-bound", type=int, default=1000)
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not 100 <= args.evaluation_budget <= 10_000:
        raise SystemExit("--evaluation-budget must lie in [100,10000]")
    if min(
        args.beam_width,
        args.beam_rounds,
        args.mutations_per_state,
        args.cover_count,
        args.offset_count,
        args.affine_count,
        args.mobius_count,
        args.affine_pool_size,
        args.mobius_pool_size,
        args.offset_height,
        args.affine_height,
        args.mobius_height,
        args.skew_numerator_bound,
        args.skew_denominator_bound,
        args.relation_batch_size,
    ) <= 0:
        raise SystemExit("counts, heights, and beam parameters must be positive")
    if not 5 <= args.cover_count <= 25:
        raise SystemExit("--cover-count must lie in [5,25]")
    if not 2 <= args.affine_pool_size <= 20:
        raise SystemExit("--affine-pool-size must lie in [2,20]")
    if not 3 <= args.mobius_pool_size <= 14:
        raise SystemExit("--mobius-pool-size must lie in [3,14]")
    if not 0 < args.chart_timeout <= 60 or not 0 < args.relation_timeout <= 60:
        raise SystemExit("subprocess timeouts must lie in (0,60]")
    if args.stack_bytes < 64_000_000:
        raise SystemExit("--stack-bytes must be at least 64MB")
    if args.certificate_prime_bound < 3:
        raise SystemExit("--certificate-prime-bound must be at least 3")


def main() -> None:
    args = build_parser().parse_args()
    validate_args(args)
    started = time.monotonic()
    evaluations, band_records = search_signed_classes(
        seed=args.seed,
        evaluation_budget=args.evaluation_budget,
        beam_width=args.beam_width,
        rounds=args.beam_rounds,
        mutations_per_state=args.mutations_per_state,
    )
    if any(evaluation.state.weight <= 3 for evaluation in evaluations):
        raise AssertionError("the search overlapped a prior weight<=3 class")
    selected_evaluations = select_diverse_evaluations(
        evaluations,
        count=args.cover_count,
        minimum_distance=args.minimum_mask_distance,
    )
    candidates = tuple(build_candidate(evaluation) for evaluation in selected_evaluations)
    plans = build_cover_chart_plans(
        candidates,
        offset_count=args.offset_count,
        affine_count=args.affine_count,
        mobius_count=args.mobius_count,
        affine_pool_size=args.affine_pool_size,
        mobius_pool_size=args.mobius_pool_size,
        offset_height=args.offset_height,
        affine_height=args.affine_height,
        mobius_height=args.mobius_height,
        skew_numerator_bound=args.skew_numerator_bound,
        skew_denominator_bound=args.skew_denominator_bound,
    )
    evaluation_digest_payload = [
        evaluation_payload(evaluation)
        for evaluation in sorted(evaluations, key=lambda item: item.state.identifier)
    ]
    evaluation_sha256 = hashlib.sha256(
        json.dumps(evaluation_digest_payload, separators=(",", ":")).encode()
    ).hexdigest()
    plan_payload = [
        {
            "id": plan.identifier,
            "kind": plan.kind,
            "candidate_index": plan.candidate_index,
            "height": plan.height_specification,
            "coefficients": [str(value) for value in plan.polynomial],
            "center": None if plan.center is None else str(plan.center),
            "scale": None if plan.scale is None else str(plan.scale),
            "matrix": None if plan.matrix is None else list(plan.matrix),
        }
        for plan in plans
    ]
    plan_sha256 = hashlib.sha256(
        json.dumps(plan_payload, separators=(",", ":")).encode()
    ).hexdigest()

    returned: dict[RationalPoint, set[str]] = {}
    run_records = []
    for run_index, plan in enumerate(plans, 1):
        raw_points, process = search_original_quartic(
            plan.polynomial,
            plan.height_specification,
            timeout=args.chart_timeout,
            stack_bytes=args.stack_bytes,
        )
        mapped, pole_count = map_plan_points(
            plan, candidates[plan.candidate_index], raw_points
        )
        for point in mapped:
            returned.setdefault(point, set()).add(plan.identifier)
        run_records.append(
            {
                "id": plan.identifier,
                "kind": plan.kind,
                "candidate_index_one_based": plan.candidate_index + 1,
                "height_specification": plan.height_specification,
                "compression_score": list(plan.compression_score),
                **process,
                "mapped_affine_image_count_including_duplicates": len(mapped),
                "mobius_pole_point_count": pole_count,
                "all_mapped_images_checked_exactly": True,
            }
        )
        if run_index % 10 == 0 or process["status"] != "completed":
            print(
                f"higher-weight charts {run_index}/{len(plans)}; "
                f"status={process['status']} unique_images={len(returned)}",
                flush=True,
            )

    relation_map = known_relation_map(candidates)
    classifications: dict[RationalPoint, dict[str, Any]] = {}
    unresolved = []
    for point in sorted(returned, key=lambda item: (item[0], item[1])):
        relation = relation_map.get(point)
        if relation is None:
            unresolved.append(point)
            continue
        if exact_linear_combination(relation) != point:
            raise AssertionError("a known seed/companion relation failed replay")
        classifications[point] = {
            "classification": "exact_seed_or_companion_in_rank29_subgroup",
            "published_basis_relation": list(relation),
            "exact_fraction_group_law_replay": True,
        }

    relation_processes = []
    still_unresolved = []
    for batch_start in range(0, len(unresolved), args.relation_batch_size):
        batch = unresolved[batch_start : batch_start + args.relation_batch_size]
        relations, process = discover_relations_batch(
            batch,
            timeout=args.relation_timeout,
            stack_bytes=args.stack_bytes,
        )
        relation_processes.append(
            {
                "batch_start_zero_based": batch_start,
                "point_count": len(batch),
                **process,
            }
        )
        for point, relation in zip(batch, relations):
            if relation is None:
                still_unresolved.append(point)
            else:
                classifications[point] = {
                    "classification": "exactly_in_published_rank29_subgroup",
                    "published_basis_relation": list(relation),
                    "exact_fraction_group_law_replay": True,
                }

    target_hit = False
    two_torsion_prime = find_two_torsion_certificate_prime(
        short_weierstrass_coefficients(), prime_bound=args.certificate_prime_bound
    )
    for point in still_unresolved:
        augmented = SHORT_PUBLIC_POINTS + (to_short_point(point),)
        signatures = find_mod2_reduction_certificate(
            short_weierstrass_coefficients(),
            augmented,
            prime_bound=args.certificate_prime_bound,
        )
        binary_rank = combined_mod2_rank(signatures, len(augmented))
        independent = binary_rank == 30
        classifications[point] = {
            "classification": (
                "exact_independent_30th_point"
                if independent
                else "unresolved_after_exact_relation_and_mod2_search"
            ),
            "augmented_mod2_rank": binary_rank,
            "two_torsion_certificate_prime": two_torsion_prime,
            "certificate_prime_bound": args.certificate_prime_bound,
            "signatures": [signature_record(signature) for signature in signatures],
        }
        if independent:
            target_hit = True
            print(
                "RANK30_DIRECTION exact finite-reduction rank 30 at "
                f"x={point[0]}",
                flush=True,
            )

    point_records = []
    for point in sorted(returned, key=lambda item: (item[0], item[1])):
        record = {
            **point_record(point),
            "source_charts": sorted(returned[point]),
            "exact_curve_membership_checked": True,
            **classifications[point],
        }
        point_records.append(record)
    if len(point_records) != len(returned):
        raise AssertionError("not every returned image received a classification")

    selected_records = []
    for index, candidate in enumerate(candidates, 1):
        state = candidate.state
        selected_records.append(
            {
                "index_one_based": index,
                "identifier": state.identifier,
                "class_mask_hex": f"0x{state.mask:08x}",
                "class_mask_binary_29_bits": f"{state.mask:029b}",
                "negative_mask_hex": f"0x{state.negative_mask:08x}",
                "weight": state.weight,
                "subset_indices_one_based": [index + 1 for index in state_indices(state)],
                "signed_coefficients": list(state.coefficients),
                "base_Q_short": point_record(candidate.base_short_point),
                "base_Q_general": point_record(from_short_point(candidate.base_short_point)),
                "beam_score": list(candidate.evaluation.score),
                "beam_best_anchor_labels": [
                    signed_parameter_label(value)
                    for value in candidate.evaluation.best_anchor_indices
                ],
                "quartic_coefficients_ascending": [
                    str(value) for value in candidate.cover.coefficients
                ],
                "all_58_signed_public_parameters_checked": True,
            }
        )

    script_path = Path(__file__).resolve()
    actual_command = " ".join(
        shlex.quote(part) for part in [sys.executable, *sys.argv]
    )
    completed_count = sum(record["status"] == "completed" for record in run_records)
    timeout_count = sum(record["status"] == "timeout" for record in run_records)
    artifact = {
        "schema_version": 1,
        "artifact_kind": "bounded_higher_weight_mod2_cover_rank30_search",
        "status": (
            "exact_rank30_target_hit"
            if target_hit
            else "bounded_search_no_certified_30th_point"
        ),
        "claim_scope": {
            "exact": (
                "class-mask exclusion, signed representatives, group sums, "
                "quartics, chart maps, curve memberships, subgroup relation "
                "replays, and any finite-reduction independence certificate"
            ),
            "bounded": (
                "deterministic beam/local evaluations and exactly the completed "
                "PARI chart boxes; no rank upper bound"
            ),
        },
        "quotient_semantics": {
            "dimension": 29,
            "class_mask": "support of coefficients modulo 2",
            "prior_exhaustive_classes_excluded": "every mask of Hamming weight <=3",
            "minimum_evaluated_weight": min(item.state.weight for item in evaluations),
            "signed_representative_rule": (
                "coefficient signs may change because +/-P_i differ by 2P_i"
            ),
            "global_negation_canonicalization": (
                "the least-index nonzero coefficient is positive"
            ),
            "complement_rule": (
                "bitwise complements are retained as distinct: v and "
                "v XOR (2^29-1) are not generally equal in E(Q)/2E(Q)"
            ),
        },
        "reproduction": {
            "command": REPRODUCING_COMMAND,
            "actual_command": actual_command,
            "python": platform.python_version(),
            "pari_gp": pari_version(),
            "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
        },
        "selection": {
            "seed": args.seed,
            "evaluation_budget": args.evaluation_budget,
            "beam_width": args.beam_width,
            "beam_rounds": args.beam_rounds,
            "mutations_per_state": args.mutations_per_state,
            "weight_bands": [list(band) for band in WEIGHT_BANDS],
            "band_records": list(band_records),
            "evaluated_manifest_sha256": evaluation_sha256,
            "evaluated_manifest_stored_inline": False,
            "retained_cover_count": len(candidates),
            "minimum_mask_hamming_distance_requested": args.minimum_mask_distance,
            "retained_covers": selected_records,
        },
        "search_budget": {
            "offset_charts_per_cover": args.offset_count,
            "affine_cross_ratio_charts_per_cover": args.affine_count,
            "mobius_cross_ratio_charts_per_cover": args.mobius_count,
            "mobius_skew_boxes_per_cover": args.mobius_count,
            "affine_normalization_pool_size": args.affine_pool_size,
            "mobius_normalization_pool_size": args.mobius_pool_size,
            "offset_height": args.offset_height,
            "affine_height": args.affine_height,
            "mobius_height": args.mobius_height,
            "skew_numerator_bound": args.skew_numerator_bound,
            "skew_denominator_interval": [1, args.skew_denominator_bound],
            "chart_timeout_seconds_each": args.chart_timeout,
            "relation_timeout_seconds_each_batch": args.relation_timeout,
            "relation_batch_size": args.relation_batch_size,
            "stack_bytes_each": args.stack_bytes,
            "certificate_prime_bound": args.certificate_prime_bound,
            "one_pass_no_retry": True,
        },
        "chart_manifest": {
            "declared_count": len(plans),
            "sha256": plan_sha256,
            "stored_inline": False,
            "kind_counts": {
                kind: sum(plan.kind == kind for plan in plans)
                for kind in sorted({plan.kind for plan in plans})
            },
        },
        "search_result": {
            "completed_chart_count": completed_count,
            "timed_out_chart_count": timeout_count,
            "other_noncompleted_chart_count": len(plans) - completed_count - timeout_count,
            "all_declared_runs_completed": completed_count == len(plans),
            "run_records": run_records,
            "unique_exact_mapped_image_count": len(returned),
            "relation_processes": relation_processes,
            "point_records": point_records,
            "exact_seed_or_companion_relation_count": sum(
                record["classification"]
                == "exact_seed_or_companion_in_rank29_subgroup"
                for record in point_records
            ),
            "exact_other_rank29_relation_count": sum(
                record["classification"] == "exactly_in_published_rank29_subgroup"
                for record in point_records
            ),
            "unresolved_after_relation_and_mod2_count": sum(
                record["classification"]
                == "unresolved_after_exact_relation_and_mod2_search"
                for record in point_records
            ),
            "certified_independent_30th_point_count": sum(
                record["classification"] == "exact_independent_30th_point"
                for record in point_records
            ),
            "rank30_target_hit": target_hit,
            "wall_seconds": time.monotonic() - started,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {args.output}")
    print(
        f"completed={completed_count}/{len(plans)} "
        f"unique_images={len(returned)} target_hit={str(target_hit).lower()}",
        flush=True,
    )


if __name__ == "__main__":
    main()
