#!/usr/bin/env python3
"""Direct quotient-aware search beyond the Fermigier rank-20 anchor.

General 2-descent backends are pathological on this specialization: PARI
spends hours in the cubic-field class group and eclib overflows native search
bounds.  This runner instead searches many exact degree-two coordinate systems
of the same elliptic curve and asks a cheaper question first: did any returned
point escape the known rank-20 subgroup in small finite quotients?

The full workflow is deliberately staged:

``smoke``
    Search the original quartic in two skew boxes and a few low-weight
    alternate covers.  This is the safe default.
``frontier``
    Scan the ``2^20-1`` represented nonzero mod-2 classes in Gray-code order,
    retaining the covers whose known points have the smallest exact coordinate
    heights.  The scan is resumable from an atomic checkpoint.
``search``
    Search the retained covers from a completed frontier checkpoint.
``all``
    Build/resume the full frontier and then search it.

Every quartic point and every map is checked over ``QQ``.  Height pairings are
used only to propose subgroup relations, which are replayed with exact
``Fraction`` arithmetic.  Finite-quotient escape is exact evidence, but a rank
claim is emitted only after ``build_finite_quotient_certificate`` certifies the
augmented point set.
"""

from __future__ import annotations

from research_runtime.supervisor import Limits, capture, capture_record, captured_run, run as supervised_run

import argparse
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import heapq
import itertools
import json
from math import isqrt
import os
from pathlib import Path
import re
import shlex
import shutil
import signal
import subprocess
import sys
import time
from typing import Any, Iterable, Sequence

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

from alternate_quartic_covers import (
    AlternateQuarticCover,
    alternate_cover,
    short_subset_sum,
)
from elliptic_candidate_record import (
    WeierstrassChange,
    build_finite_quotient_certificate,
    fraction_text,
    is_on_weierstrass_curve,
    model_from_record,
    point_from_record,
    point_sequence_sha256,
    source_point_to_target,
    verify_finite_quotient_certificate,
)
from fermigier_mestre import FermigierMestreFamily
from overcomplete_quotient_bank import (
    OvercompleteQuotientBank,
    build_overcomplete_quotient_bank,
    evaluate_candidate_batches_against_bank,
    load_overcomplete_quotient_bank,
    save_overcomplete_quotient_bank,
)
from search_nagao_u135_alternate_covers import best_cross_ratio_charts
from search_nagao_u42_skew_height import (
    centered_unimodular_matrix,
    exact_linear_combination,
    map_chart_point,
    short_add as fast_short_add,
    transform_binary_quartic,
)


Q = Fraction
RationalPoint = tuple[Fraction, Fraction]
ROOT = Path(__file__).resolve().parents[2]
ANCHOR_ADAPTER_U = Q(28_917, 20)
ANCHOR_T = 2 * ANCHOR_ADAPTER_U
EXPECTED_STRONG_BASIS_SHA256 = (
    "1c637a6a592604fd8aee99a17c152c82739854155324e98c7009810ded5f5a73"
)
PLAN_SCHEMA = "elliptic-curves.fermigier-rank20-direction-plan.v1"
RESULT_SCHEMA = "elliptic-curves.fermigier-rank20-direction-search.v1"
DEFAULT_CANDIDATE_INPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/"
    "elliptic_curve_candidate_fermigier_mestre_v1_u28917_20.json"
)
DEFAULT_PLAN = (
    ROOT
    / "artifacts/local/elliptic-curves/"
    "fermigier_rank20_direction_plan.json"
)
DEFAULT_BANK = (
    ROOT
    / "artifacts/local/elliptic-curves/"
    "fermigier_rank20_overcomplete_quotient_bank.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/local/elliptic-curves/"
    "fermigier_rank20_direction_search.json"
)


@dataclass(frozen=True)
class SearchBox:
    identifier: str
    numerator_bound: int
    denominator_lower: int
    denominator_upper: int

    def __post_init__(self) -> None:
        if not 1 <= self.denominator_lower <= self.denominator_upper:
            raise ValueError("invalid denominator interval")
        if self.numerator_bound < self.denominator_upper:
            raise ValueError("PARI requires numerator bound >= denominator upper bound")

    @property
    def gp_height(self) -> str:
        if self.denominator_lower == 1:
            return f"[{self.numerator_bound},{self.denominator_upper}]"
        return (
            f"[{self.numerator_bound},"
            f"[{self.denominator_lower},{self.denominator_upper}]]"
        )


# These boxes are disjoint in denominator and complementary to the previous
# balanced H<=2,000,000, denominator<=13,000 ratpoints run.
ORIGINAL_SKEW_BOXES = (
    SearchBox("d000001_000010", 10_000_000_000, 1, 10),
    SearchBox("d000011_000100", 1_000_000_000, 11, 100),
    SearchBox("d000101_001000", 100_000_000, 101, 1_000),
    SearchBox("d001001_013000", 10_000_000, 1_001, 13_000),
    SearchBox("d013001_050000", 2_000_000, 13_001, 50_000),
    SearchBox("d050001_200000", 500_000, 50_001, 200_000),
)


@dataclass(frozen=True)
class AnchorData:
    model: tuple[Fraction, Fraction, Fraction, Fraction, Fraction]
    basis: tuple[RationalPoint, ...]
    basis_sha256: str
    known_signless_points: frozenset[RationalPoint]
    known_quartic_abscissas: tuple[Fraction, ...]
    candidate_record: dict[str, Any]


@dataclass(frozen=True)
class FrontierEntry:
    maximum_parameter_bit_length: int
    sum_parameter_bit_lengths: int
    mask: int

    @property
    def score(self) -> tuple[int, int, int]:
        return (
            self.maximum_parameter_bit_length,
            self.sum_parameter_bit_lengths,
            self.mask,
        )

    @property
    def subset_indices(self) -> tuple[int, ...]:
        return tuple(
            index for index in range(self.mask.bit_length()) if self.mask >> index & 1
        )

    def to_record(self) -> dict[str, Any]:
        return {
            "maximum_parameter_bit_length": self.maximum_parameter_bit_length,
            "sum_parameter_bit_lengths": self.sum_parameter_bit_lengths,
            "mask": hex(self.mask),
            "subset_weight": self.mask.bit_count(),
            "subset_indices_one_based": [index + 1 for index in self.subset_indices],
        }

    @classmethod
    def from_record(cls, record: dict[str, Any]) -> "FrontierEntry":
        return cls(
            maximum_parameter_bit_length=int(
                record["maximum_parameter_bit_length"]
            ),
            sum_parameter_bit_lengths=int(record["sum_parameter_bit_lengths"]),
            mask=int(str(record["mask"]), 16),
        )


@dataclass(frozen=True)
class CoverSearchPlan:
    frontier: FrontierEntry
    cover: AlternateQuarticCover
    charts: tuple[Any, ...]

    @property
    def identifier(self) -> str:
        return "q_" + "_".join(
            str(index + 1) for index in self.frontier.subset_indices
        )


@dataclass(frozen=True)
class StageSettings:
    cover_limit: int
    charts_per_cover: int
    original_box_limit: int
    original_center_count: int
    original_center_shifts: tuple[int, ...]
    pilot_height: int
    escalation_height: int
    escalation_charts: int
    deep_height: int
    deep_charts: int
    skew_cover_count: int
    skew_box_limit: int

    def __post_init__(self) -> None:
        if self.cover_limit <= 0 or self.charts_per_cover <= 0:
            raise ValueError("cover and chart counts must be positive")
        if self.original_box_limit < 0 or self.original_center_count < 0:
            raise ValueError("original-search counts must be nonnegative")
        if not self.original_center_shifts:
            raise ValueError("at least one original centered-chart shift is required")
        if self.pilot_height <= 0 or self.escalation_height <= 0 or self.deep_height <= 0:
            raise ValueError("point-search heights must be positive")
        if self.escalation_charts < 0 or self.deep_charts < 0:
            raise ValueError("escalation/deep chart counts must be nonnegative")
        if self.skew_cover_count < 0 or self.skew_box_limit < 0:
            raise ValueError("skew-search counts must be nonnegative")


def projective_height(value: Fraction) -> int:
    value = Q(value)
    return max(abs(value.numerator), value.denominator)


def canonical_sign(point: RationalPoint) -> RationalPoint:
    """Choose one exact representative from ``{P,-P}`` on a short model."""

    x_value, y_value = (Q(value) for value in point)
    return (x_value, y_value) if y_value >= 0 else (x_value, -y_value)


def point_record(point: RationalPoint | None) -> dict[str, str] | None:
    if point is None:
        return None
    return {"x": fraction_text(point[0]), "y": fraction_text(point[1])}


def point_from_xy_record(record: dict[str, Any]) -> RationalPoint:
    return Q(record["x"]), Q(record["y"])


def rational_square_root(value: Fraction) -> Fraction:
    value = Q(value)
    if value < 0:
        raise ValueError("a negative rational is not a rational square")
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    if numerator**2 != value.numerator or denominator**2 != value.denominator:
        raise ValueError("the rational is not a square")
    return Q(numerator, denominator)


def _atomic_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def _direct_quartic_images(abscissas: Iterable[Fraction]) -> set[RationalPoint]:
    images: set[RationalPoint] = set()
    for abscissa in abscissas:
        square = FermigierMestreFamily.quartic_value(ANCHOR_T, Q(abscissa))
        ordinate = rational_square_root(square)
        if ordinate == 0:
            continue
        image = FermigierMestreFamily.quartic_point_to_jacobian(
            ANCHOR_T, (Q(abscissa), ordinate)
        )
        images.add(canonical_sign(image))
    return images


def load_anchor(path: Path = DEFAULT_CANDIDATE_INPUT) -> AnchorData:
    record = json.loads(path.read_text())
    identity = record["identity"]["canonical_parameter"]
    if identity["name"] != "u" or Q(identity["value"]) != ANCHOR_ADAPTER_U:
        raise AssertionError("the candidate artifact has the wrong Fermigier identity")

    model = model_from_record(
        record["models"]["legacy_normalized_short_jacobian"]["coefficients"]
    )
    expected_model = FermigierMestreFamily.coefficients(ANCHOR_T)
    if model != expected_model or any(model[:3]):
        raise AssertionError("the pinned legacy short model changed")

    basis = tuple(
        point_from_record(item)
        for item in record["bounded_saturation_status"]["returned_legacy_basis"]
    )
    basis_sha256 = point_sequence_sha256(basis)
    recorded_hash = record["bounded_saturation_status"][
        "returned_legacy_point_sequence_sha256"
    ]
    if (
        len(basis) != 20
        or basis_sha256 != EXPECTED_STRONG_BASIS_SHA256
        or recorded_hash != basis_sha256
        or any(not is_on_weierstrass_curve(model, point) for point in basis)
    ):
        raise AssertionError("the bounded-saturation rank-20 basis changed")

    change = WeierstrassChange.from_values(
        record["exact_transformations"][
            "canonical_to_legacy_normalized_short"
        ]["change_u_r_s_t"]
    )
    known: set[RationalPoint] = {canonical_sign(point) for point in basis}
    for item in record["complete_point_pool"]["difference_pool"]:
        canonical = point_from_record(item["canonical_point"])
        legacy = source_point_to_target(canonical, change)
        if not is_on_weierstrass_curve(model, legacy):
            raise AssertionError("a transported historical point left the short curve")
        known.add(canonical_sign(legacy))

    abscissas = tuple(
        Q(value) for value in record["complete_point_pool"]["abscissas"]
    )
    known.update(_direct_quartic_images(abscissas))
    visible = FermigierMestreFamily.known_quartic_points(ANCHOR_T)
    known.update(
        canonical_sign(
            FermigierMestreFamily.quartic_point_to_jacobian(ANCHOR_T, point)
        )
        for point in visible
        if point[1]
    )
    if any(not is_on_weierstrass_curve(model, point) for point in known):
        raise AssertionError("the decontamination set contains an off-curve point")

    return AnchorData(
        model=model,
        basis=basis,
        basis_sha256=basis_sha256,
        known_signless_points=frozenset(known),
        known_quartic_abscissas=abscissas,
        candidate_record=record,
    )


def score_base_point(
    coefficient_a: Fraction,
    basis: Sequence[RationalPoint],
    base_point: RationalPoint,
) -> tuple[int, int]:
    """Score one exact degree-two chart by known parameter bit lengths."""

    x_base, y_base = base_point
    bit_lengths: list[int] = []
    for x_value, y_value in basis:
        if (x_value, y_value) == base_point:
            continue
        if x_value == x_base:
            if y_value != -y_base or y_base == 0:
                raise AssertionError("an unexpected exceptional cover parameter")
            parameter = -(3 * x_base**2 + coefficient_a) / (2 * y_base)
        else:
            parameter = (y_value + y_base) / (x_value - x_base)
        bit_lengths.append(projective_height(parameter).bit_length())
    if not bit_lengths:
        raise AssertionError("a cover has no finite known parameters")
    return max(bit_lengths), sum(bit_lengths)


def _heap_item(entry: FrontierEntry) -> tuple[int, int, int, FrontierEntry]:
    return (
        -entry.maximum_parameter_bit_length,
        -entry.sum_parameter_bit_lengths,
        -entry.mask,
        entry,
    )


def _push_frontier(
    heap: list[tuple[int, int, int, FrontierEntry]],
    entry: FrontierEntry,
    retain_count: int,
) -> None:
    item = _heap_item(entry)
    if len(heap) < retain_count:
        heapq.heappush(heap, item)
        return
    if entry.score < heap[0][3].score:
        heapq.heapreplace(heap, item)


def select_frontier_entries(
    frontier: Sequence[FrontierEntry],
    weight_frontier: Sequence[FrontierEntry],
    *,
    global_count: int,
    maximum_count: int,
) -> tuple[FrontierEntry, ...]:
    selected = {entry.mask: entry for entry in sorted(frontier, key=lambda e: e.score)[:global_count]}
    for entry in weight_frontier:
        selected.setdefault(entry.mask, entry)
    return tuple(sorted(selected.values(), key=lambda entry: entry.score)[:maximum_count])


def _frontier_record(
    anchor: AnchorData,
    *,
    status: str,
    target: int,
    total: int,
    processed: int,
    previous_gray: int,
    current: RationalPoint | None,
    heap: Sequence[tuple[int, int, int, FrontierEntry]],
    weight_best: dict[int, FrontierEntry],
    retain_count: int,
    global_count: int,
    maximum_count: int,
    progress_interval: int,
) -> dict[str, Any]:
    frontier = tuple(sorted((item[3] for item in heap), key=lambda entry: entry.score))
    weights = tuple(weight_best[weight] for weight in sorted(weight_best))
    selected = select_frontier_entries(
        frontier,
        weights,
        global_count=global_count,
        maximum_count=maximum_count,
    )
    return {
        "schema": PLAN_SCHEMA,
        "status": status,
        "candidate_key": "fermigier-mestre-v1:u=28917/20",
        "tuple_parameter_T": fraction_text(ANCHOR_T),
        "basis_count": len(anchor.basis),
        "basis_sha256": anchor.basis_sha256,
        "configuration": {
            "total_nonzero_classes": total,
            "scan_target": target,
            "retain_count": retain_count,
            "global_count": global_count,
            "maximum_count": maximum_count,
            "progress_interval": progress_interval,
        },
        "state": {
            "processed_integer": processed,
            "previous_gray": previous_gray,
            "current_point": point_record(current),
        },
        "frontier": [entry.to_record() for entry in frontier],
        "weight_frontier": [entry.to_record() for entry in weights],
        "selected": [entry.to_record() for entry in selected],
        "scope": (
            "all represented nonzero mod-2 classes"
            if target == total
            else f"Gray-code prefix through integer {target}; exploratory only"
        ),
    }


def scan_frontier(
    anchor: AnchorData,
    checkpoint: Path,
    *,
    retain_count: int = 64,
    global_count: int = 40,
    maximum_count: int = 60,
    progress_interval: int = 65_536,
    scan_limit: int = 0,
    resume: bool = True,
) -> tuple[FrontierEntry, ...]:
    """Scan exact represented classes, atomically checkpointing resumable state."""

    if not 0 < global_count <= retain_count or maximum_count < global_count:
        raise ValueError("invalid frontier retention counts")
    if progress_interval <= 0 or scan_limit < 0:
        raise ValueError("progress and scan bounds must be nonnegative")
    total = (1 << len(anchor.basis)) - 1
    target = total if scan_limit == 0 else min(total, scan_limit)
    coefficient_a = Q(anchor.model[3])

    processed = 0
    previous_gray = 0
    current: RationalPoint | None = None
    heap: list[tuple[int, int, int, FrontierEntry]] = []
    weight_best: dict[int, FrontierEntry] = {}

    if resume and checkpoint.exists():
        prior = json.loads(checkpoint.read_text())
        configuration = prior.get("configuration", {})
        expected = {
            "total_nonzero_classes": total,
            "scan_target": target,
            "retain_count": retain_count,
            "global_count": global_count,
            "maximum_count": maximum_count,
            "progress_interval": progress_interval,
        }
        if (
            prior.get("schema") != PLAN_SCHEMA
            or prior.get("basis_sha256") != anchor.basis_sha256
            or configuration != expected
        ):
            raise ValueError("the existing frontier checkpoint has different inputs")
        if prior["status"] in {"complete", "bounded-complete"}:
            return tuple(FrontierEntry.from_record(item) for item in prior["selected"])
        state = prior["state"]
        processed = int(state["processed_integer"])
        previous_gray = int(state["previous_gray"])
        current_record = state["current_point"]
        current = None if current_record is None else point_from_xy_record(current_record)
        for item in prior["frontier"]:
            _push_frontier(heap, FrontierEntry.from_record(item), retain_count)
        weight_best = {
            int(item["subset_weight"]): FrontierEntry.from_record(item)
            for item in prior["weight_frontier"]
        }
        if previous_gray != (processed ^ (processed >> 1)):
            raise AssertionError("the checkpoint Gray-code state is inconsistent")
        expected_current = (
            None
            if previous_gray == 0
            else short_subset_sum(
                anchor.model,
                anchor.basis,
                tuple(
                    index
                    for index in range(len(anchor.basis))
                    if previous_gray >> index & 1
                ),
            )
        )
        if current != expected_current:
            raise AssertionError("the checkpoint current point does not match its Gray word")
        if current is not None and not is_on_weierstrass_curve(anchor.model, current):
            raise AssertionError("the checkpoint current point left the curve")
        print(
            f"F20SEARCH|stage=frontier|status=resume|processed={processed}|target={target}",
            flush=True,
        )

    started = time.monotonic()

    def save(status: str) -> None:
        _atomic_json(
            checkpoint,
            _frontier_record(
                anchor,
                status=status,
                target=target,
                total=total,
                processed=processed,
                previous_gray=previous_gray,
                current=current,
                heap=heap,
                weight_best=weight_best,
                retain_count=retain_count,
                global_count=global_count,
                maximum_count=maximum_count,
                progress_interval=progress_interval,
            ),
        )

    # Keep the last durable checkpoint transactionally consistent.  A SIGINT
    # may arrive between arbitrary Python bytecodes, so the interrupt handler
    # deliberately leaves this last completed checkpoint untouched rather than
    # serializing partially committed in-memory frontier state.
    if not checkpoint.exists():
        save("partial")
    last_checkpoint_processed = processed

    try:
        for integer in range(processed + 1, target + 1):
            gray = integer ^ (integer >> 1)
            changed = gray ^ previous_gray
            if changed == 0 or changed & (changed - 1):
                raise AssertionError("consecutive Gray words did not differ in one bit")
            index = changed.bit_length() - 1
            summand = anchor.basis[index]
            if not (gray >> index & 1):
                summand = summand[0], -summand[1]
            next_current = fast_short_add(coefficient_a, current, summand)
            if next_current is None:
                raise AssertionError("a represented nonzero class summed to infinity")
            maximum_bits, sum_bits = score_base_point(
                coefficient_a, anchor.basis, next_current
            )
            entry = FrontierEntry(maximum_bits, sum_bits, gray)
            _push_frontier(heap, entry, retain_count)
            weight = gray.bit_count()
            prior_weight = weight_best.get(weight)
            if prior_weight is None or entry.score < prior_weight.score:
                weight_best[weight] = entry
            current = next_current
            processed = integer
            previous_gray = gray
            if processed % progress_interval == 0:
                if not is_on_weierstrass_curve(anchor.model, current):
                    raise AssertionError("the Gray-code group walk left the curve")
                save("partial")
                last_checkpoint_processed = processed
                print(
                    "F20SEARCH|stage=frontier|status=progress|"
                    f"processed={processed}|target={target}|"
                    f"elapsed={time.monotonic()-started:.3f}",
                    flush=True,
                )
    except KeyboardInterrupt:
        print(
            "F20SEARCH|stage=frontier|status=interrupted|"
            f"durable_processed={last_checkpoint_processed}",
            flush=True,
        )
        raise

    if current is not None and not is_on_weierstrass_curve(anchor.model, current):
        raise AssertionError("the completed Gray-code group walk left the curve")
    save("complete" if target == total else "bounded-complete")
    final = json.loads(checkpoint.read_text())
    print(
        "F20SEARCH|stage=frontier|status=complete|"
        f"processed={processed}|selected={len(final['selected'])}|"
        f"elapsed={time.monotonic()-started:.3f}",
        flush=True,
    )
    return tuple(FrontierEntry.from_record(item) for item in final["selected"])


def low_weight_frontier(
    anchor: AnchorData,
    *,
    maximum_weight: int = 2,
    retain_count: int = 12,
) -> tuple[FrontierEntry, ...]:
    """Cheap deterministic plan used by the smoke stage."""

    if maximum_weight <= 0 or retain_count <= 0:
        raise ValueError("low-weight bounds must be positive")
    entries: list[FrontierEntry] = []
    coefficient_a = Q(anchor.model[3])
    for weight in range(1, min(maximum_weight, len(anchor.basis)) + 1):
        for indices in itertools.combinations(range(len(anchor.basis)), weight):
            base_point = short_subset_sum(anchor.model, anchor.basis, indices)
            if base_point is None:
                raise AssertionError("a nonempty rank-20 subset summed to infinity")
            maximum_bits, sum_bits = score_base_point(
                coefficient_a, anchor.basis, base_point
            )
            mask = sum(1 << index for index in indices)
            entries.append(FrontierEntry(maximum_bits, sum_bits, mask))
    entries.sort(key=lambda entry: entry.score)
    selected = {entry.mask: entry for entry in entries[:retain_count]}
    for weight in range(1, maximum_weight + 1):
        best = next((entry for entry in entries if entry.mask.bit_count() == weight), None)
        if best is not None:
            selected.setdefault(best.mask, best)
    return tuple(sorted(selected.values(), key=lambda entry: entry.score)[:retain_count])


def build_cover_plans(
    anchor: AnchorData,
    entries: Sequence[FrontierEntry],
    *,
    charts_per_cover: int,
    cover_limit: int,
) -> tuple[CoverSearchPlan, ...]:
    if charts_per_cover <= 0 or cover_limit <= 0:
        raise ValueError("cover and chart counts must be positive")
    plans: list[CoverSearchPlan] = []
    for position, entry in enumerate(entries[:cover_limit], start=1):
        base_point = short_subset_sum(
            anchor.model, anchor.basis, entry.subset_indices
        )
        if base_point is None:
            raise AssertionError("a selected cover base point vanished")
        cover = alternate_cover(anchor.model, base_point)
        charts = best_cross_ratio_charts(
            cover, anchor.basis, count=charts_per_cover
        )
        plans.append(CoverSearchPlan(entry, cover, charts))
        print(
            "F20SEARCH|stage=plan|status=cover|"
            f"position={position}|count={min(len(entries), cover_limit)}|"
            f"id=q_{'_'.join(str(i+1) for i in entry.subset_indices)}|"
            f"max_bits={entry.maximum_parameter_bit_length}",
            flush=True,
        )
    return tuple(plans)


def gp_rational(value: Fraction | int) -> str:
    value = Q(value)
    return f"({value.numerator}/{value.denominator})"


def gp_point(point: RationalPoint) -> str:
    return f"[{gp_rational(point[0])},{gp_rational(point[1])}]"


def gp_curve(model: Sequence[Fraction]) -> str:
    return ",".join(gp_rational(value) for value in model)


def gp_quartic(coefficients_ascending: Sequence[Fraction]) -> str:
    if len(coefficients_ascending) != 5:
        raise ValueError("five ascending quartic coefficients are required")
    return "+".join(
        f"{gp_rational(coefficient)}*x^{power}"
        for power, coefficient in enumerate(coefficients_ascending)
    )




def run_gp_once(program, *, timeout, stack_bytes):
    if timeout<=0 or stack_bytes<8_000_000:raise ValueError("timeout and stack bounds must be positive")
    executable=shutil.which('gp')
    if executable is None:return None,{'status':'unavailable','wall_seconds':0.0}
    record=capture_record([executable,'-q','-s',str(stack_bytes)],input_text=program,
        limits=Limits(timeout,max(512_000_000,2*stack_bytes),pari_stack_bytes=stack_bytes))
    status=record['outcome']
    fatal=[line for line in record['stderr'].splitlines() if '***' in line and "Warning:" not in line]
    if status=='strict_wall_timeout':return None,{**record,'status':'timeout'}
    if status!='completed':return None,{**record,'status':status}
    if fatal or record['returncode']:
        return None,{**record,'status':'pari_error','error':' '.join(fatal or record['stderr'].splitlines())[:2000]}
    return record['stdout'],{**record,'status':'completed'}


def parse_gp_points(text: str) -> tuple[RationalPoint, ...]:
    pairs = re.findall(
        r"\[(-?\d+(?:/\d+)?),\s*(-?\d+(?:/\d+)?)\]", text
    )
    return tuple((Q(x_value), Q(y_value)) for x_value, y_value in pairs)


def search_quartic(
    coefficients_ascending: Sequence[Fraction],
    height_specification: str,
    *,
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[RationalPoint, ...], dict[str, Any]]:
    program = "\n".join(
        (
            f"Q={gp_quartic(coefficients_ascending)};gettime();",
            f"R=hyperellratpoints(Q,{height_specification});",
            'print("PARI_MILLISECONDS ",gettime());',
            'print("POINTS_BEGIN");print(R);print("POINTS_END");',
            "quit",
        )
    ) + "\n"
    output, process = run_gp_once(
        program, timeout=timeout, stack_bytes=stack_bytes
    )
    record = {
        **process,
        "height_specification": height_specification,
        "timeout_seconds": timeout,
        "pari_stack_bytes": stack_bytes,
        "retried": False,
    }
    if output is None:
        return (), record
    marker = re.search(r"POINTS_BEGIN\n(.*?)\nPOINTS_END", output, re.DOTALL)
    milliseconds = re.search(r"^PARI_MILLISECONDS (\d+)$", output, re.MULTILINE)
    if marker is None or milliseconds is None:
        raise AssertionError("PARI omitted point-search output markers")
    points = parse_gp_points(marker.group(1))
    record.update(
        {
            "pari_milliseconds": int(milliseconds.group(1)),
            "signed_point_count": len(points),
            "distinct_abscissa_count": len({point[0] for point in points}),
        }
    )
    return points, record


def _assigned(index: int, *, shard_index: int, shard_count: int) -> bool:
    return index % shard_count == shard_index


def absorb_curve_points(
    anchor: AnchorData,
    discoveries: dict[RationalPoint, set[str]],
    points: Iterable[RationalPoint],
    source: str,
) -> tuple[int, int]:
    before = len(discoveries)
    known_hits = 0
    for point in points:
        point = (Q(point[0]), Q(point[1]))
        if not is_on_weierstrass_curve(anchor.model, point):
            raise AssertionError("a search map returned an off-curve point")
        normalized = canonical_sign(point)
        if normalized in anchor.known_signless_points:
            known_hits += 1
            continue
        discoveries.setdefault(normalized, set()).add(source)
    return len(discoveries) - before, known_hits


def map_original_quartic_points(
    points: Iterable[RationalPoint],
) -> tuple[RationalPoint, ...]:
    mapped: list[RationalPoint] = []
    for point in points:
        if point[1] == 0:
            continue
        image = FermigierMestreFamily.quartic_point_to_jacobian(ANCHOR_T, point)
        mapped.append(image)
    return tuple(mapped)


def run_original_searches(
    anchor: AnchorData,
    settings: StageSettings,
    discoveries: dict[RationalPoint, set[str]],
    *,
    timeout: float,
    stack_bytes: int,
    shard_index: int,
    shard_count: int,
) -> list[dict[str, Any]]:
    quartic = tuple(reversed(FermigierMestreFamily.quartic_coefficients(ANCHOR_T)))
    runs: list[dict[str, Any]] = []
    job_index = 0
    for box in ORIGINAL_SKEW_BOXES[: settings.original_box_limit]:
        assigned = _assigned(job_index, shard_index=shard_index, shard_count=shard_count)
        job_index += 1
        if not assigned:
            continue
        raw, process = search_quartic(
            quartic,
            box.gp_height,
            timeout=timeout,
            stack_bytes=stack_bytes,
        )
        mapped = map_original_quartic_points(raw)
        source = f"original:skew:{box.identifier}"
        gained, known_hits = absorb_curve_points(
            anchor, discoveries, mapped, source
        )
        record = {
            "kind": "original-skew",
            "job_index": job_index - 1,
            "source": source,
            "numerator_absolute_bound": box.numerator_bound,
            "denominator_lower_bound": box.denominator_lower,
            "denominator_upper_bound": box.denominator_upper,
            "mapped_exact_curve_points_with_sign": len(mapped),
            "known_or_prior_hits_with_sign": known_hits,
            "new_global_signless_curve_points": gained,
            **process,
        }
        runs.append(record)
        print(
            "F20SEARCH|stage=original|"
            f"source={source}|status={process['status']}|"
            f"discoveries={len(discoveries)}",
            flush=True,
        )

    centers = sorted(
        set(anchor.known_quartic_abscissas),
        key=lambda value: (-projective_height(value), value),
    )[: settings.original_center_count]
    for center_position, center in enumerate(centers):
        for shift in settings.original_center_shifts:
            assigned = _assigned(
                job_index, shard_index=shard_index, shard_count=shard_count
            )
            job_index += 1
            if not assigned:
                continue
            matrix = centered_unimodular_matrix(center, shift)
            transformed = transform_binary_quartic(quartic, matrix)
            raw, process = search_quartic(
                transformed,
                str(settings.pilot_height),
                timeout=timeout,
                stack_bytes=stack_bytes,
            )
            poles = 0
            original_points: list[RationalPoint] = []
            for transformed_point in raw:
                point = map_chart_point(transformed_point, matrix)
                if point is None:
                    poles += 1
                    continue
                original_points.append(point)
            mapped = map_original_quartic_points(original_points)
            source = f"original:center:{center_position:02d}:shift:{shift}"
            gained, known_hits = absorb_curve_points(
                anchor, discoveries, mapped, source
            )
            record = {
                "kind": "original-centered-chart",
                "job_index": job_index - 1,
                "source": source,
                "center": fraction_text(center),
                "shift": shift,
                "matrix_a_b_c_d": list(matrix),
                "points_at_chart_pole": poles,
                "mapped_exact_curve_points_with_sign": len(mapped),
                "known_or_prior_hits_with_sign": known_hits,
                "new_global_signless_curve_points": gained,
                **process,
            }
            runs.append(record)
            print(
                "F20SEARCH|stage=original-chart|"
                f"source={source}|status={process['status']}|"
                f"discoveries={len(discoveries)}",
                flush=True,
            )
    return runs


def run_cover_chart(
    plan: CoverSearchPlan,
    chart: Any,
    height_specification: str,
    *,
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[RationalPoint, ...], dict[str, Any]]:
    transformed = transform_binary_quartic(plan.cover.coefficients, chart.matrix)
    raw, process = search_quartic(
        transformed,
        height_specification,
        timeout=timeout,
        stack_bytes=stack_bytes,
    )
    mapped: list[RationalPoint] = []
    poles = 0
    for transformed_point in raw:
        cover_point = map_chart_point(transformed_point, chart.matrix)
        if cover_point is None:
            poles += 1
            continue
        mapped.append(plan.cover.cover_point_to_curve(cover_point))
    return tuple(mapped), {
        **process,
        "points_at_chart_pole": poles,
        "mapped_exact_curve_points_with_sign": len(mapped),
    }


def run_alternate_cover_searches(
    anchor: AnchorData,
    plans: Sequence[CoverSearchPlan],
    settings: StageSettings,
    discoveries: dict[RationalPoint, set[str]],
    *,
    timeout: float,
    stack_bytes: int,
    shard_index: int,
    shard_count: int,
) -> list[dict[str, Any]]:
    runs: list[dict[str, Any]] = []
    pilot_yields: list[tuple[int, CoverSearchPlan, Any]] = []
    assigned_plans = [
        plan
        for index, plan in enumerate(plans)
        if _assigned(index, shard_index=shard_index, shard_count=shard_count)
    ]

    def execute(
        plan: CoverSearchPlan,
        chart: Any,
        *,
        stage: str,
        height_specification: str,
    ) -> int:
        mapped, process = run_cover_chart(
            plan,
            chart,
            height_specification,
            timeout=timeout,
            stack_bytes=stack_bytes,
        )
        source = (
            f"alternate:{plan.identifier}:{stage}:"
            f"{'_'.join(str(index+1) for index in chart.basis_indices)}"
        )
        gained, known_hits = absorb_curve_points(
            anchor, discoveries, mapped, source
        )
        runs.append(
            {
                "kind": "alternate-cover",
                "source": source,
                "stage": stage,
                "cover_id": plan.identifier,
                "cover_mask": hex(plan.frontier.mask),
                "cover_subset_indices_one_based": [
                    index + 1 for index in plan.frontier.subset_indices
                ],
                "cover_score": {
                    "maximum_parameter_bit_length": plan.frontier.maximum_parameter_bit_length,
                    "sum_parameter_bit_lengths": plan.frontier.sum_parameter_bit_lengths,
                },
                "normalizing_basis_indices_one_based": [
                    index + 1 for index in chart.basis_indices
                ],
                "matrix_a_b_c_d": list(chart.matrix),
                "known_or_prior_hits_with_sign": known_hits,
                "new_global_signless_curve_points": gained,
                **process,
            }
        )
        return gained

    for position, plan in enumerate(assigned_plans, start=1):
        for chart in plan.charts:
            gained = execute(
                plan,
                chart,
                stage="pilot",
                height_specification=str(settings.pilot_height),
            )
            pilot_yields.append((gained, plan, chart))
        print(
            "F20SEARCH|stage=alternate-pilot|"
            f"position={position}|count={len(assigned_plans)}|"
            f"id={plan.identifier}|discoveries={len(discoveries)}",
            flush=True,
        )

    pilot_yields.sort(
        key=lambda item: (-item[0], item[1].frontier.score, item[2].score)
    )
    escalation_yields: list[tuple[int, CoverSearchPlan, Any]] = []
    for _, plan, chart in pilot_yields[: settings.escalation_charts]:
        gained = execute(
            plan,
            chart,
            stage="escalation",
            height_specification=str(settings.escalation_height),
        )
        escalation_yields.append((gained, plan, chart))
    escalation_yields.sort(
        key=lambda item: (-item[0], item[1].frontier.score, item[2].score)
    )
    for _, plan, chart in escalation_yields[: settings.deep_charts]:
        execute(
            plan,
            chart,
            stage="deep",
            height_specification=str(settings.deep_height),
        )

    skew_plan_ids = {
        plan.identifier for plan in plans[: settings.skew_cover_count]
    }
    for plan in assigned_plans:
        if plan.identifier not in skew_plan_ids:
            continue
        chart = plan.charts[0]
        for box in ORIGINAL_SKEW_BOXES[: settings.skew_box_limit]:
            execute(
                plan,
                chart,
                stage=f"skew-{box.identifier}",
                height_specification=box.gp_height,
            )
    return runs


def relation_proposals(
    model: Sequence[Fraction],
    basis: Sequence[RationalPoint],
    points: Sequence[RationalPoint],
    *,
    timeout: float,
    stack_bytes: int,
    batch_size: int = 32,
) -> tuple[tuple[tuple[int, ...] | None, bool, str], ...]:
    """Propose height-pairing relations and replay every success exactly."""

    if batch_size <= 0:
        raise ValueError("the relation batch size must be positive")
    if not points:
        return ()
    answers: list[tuple[tuple[int, ...] | None, bool, str]] = []
    coefficient_a = Q(model[3])
    for start in range(0, len(points), batch_size):
        batch = points[start : start + batch_size]
        commands = [
            "default(realprecision,120);",
            f"E=ellinit([{gp_curve(model)}]);",
            f"B=[{','.join(gp_point(point) for point in basis)}];",
            "H=ellheightmatrix(E,B);",
        ]
        for index, point in enumerate(batch):
            commands.extend(
                (
                    f"Q={gp_point(point)};",
                    "V=vector(#B,j,ellheight(E,B[j],Q))~;",
                    "C=round(matsolve(H,V));",
                    "S=[0];for(j=1,#B,S=elladd(E,S,ellmul(E,B[j],C[j])));",
                    f'print("RELATION_{index} ",Vec(C)," EXACT ",S==Q);',
                )
            )
        commands.append("quit")
        output, process = run_gp_once(
            "\n".join(commands) + "\n",
            timeout=timeout,
            stack_bytes=stack_bytes,
        )
        if output is None:
            answers.extend((None, False, process["status"]) for _ in batch)
            continue
        for index, point in enumerate(batch):
            match = re.search(
                rf"^RELATION_{index} \[(.*?)\] EXACT ([01])$",
                output,
                re.MULTILINE,
            )
            if match is None:
                answers.append((None, False, "missing_output"))
                continue
            body = match.group(1).strip()
            relation = (
                ()
                if not body
                else tuple(int(value.strip()) for value in body.split(","))
            )
            exact = match.group(2) == "1"
            if len(relation) != len(basis):
                raise AssertionError("PARI returned a relation with the wrong length")
            if exact and exact_linear_combination(
                coefficient_a, basis, relation
            ) != point:
                raise AssertionError("a height-proposed relation failed exact replay")
            answers.append((relation, exact, "completed"))
        print(
            "F20SEARCH|stage=relations|status=batch|"
            f"processed={min(start+len(batch),len(points))}|count={len(points)}|"
            f"process={process['status']}",
            flush=True,
        )
    if len(answers) != len(points):
        raise AssertionError("relation batching lost a candidate")
    return tuple(answers)


def load_or_build_bank(
    anchor: AnchorData,
    path: Path,
    *,
    prime_bound: int,
    row_target: int,
) -> tuple[OvercompleteQuotientBank, str]:
    if path.exists():
        bank = load_overcomplete_quotient_bank(
            path, model=anchor.model, known_points=anchor.basis
        )
        if (
            bank.prime_bound == prime_bound
            and bank.row_target_per_modulus == row_target
        ):
            return bank, "loaded"
    print(
        "F20SEARCH|stage=quotient-bank|status=start|"
        f"prime_bound={prime_bound}|row_target={row_target}",
        flush=True,
    )
    bank = build_overcomplete_quotient_bank(
        anchor.model,
        anchor.basis,
        moduli=(2, 3, 5),
        prime_bound=prime_bound,
        row_target_per_modulus=row_target,
        require_full_baseline_rank=True,
        progress=True,
    )
    save_overcomplete_quotient_bank(bank, path)
    return bank, "built"


def escaping_labels_from_profile(record: dict[str, Any] | None) -> tuple[str, ...]:
    if record is None:
        return ()
    if "escaping_labels" in record:
        return tuple(sorted(str(value) for value in record["escaping_labels"]))
    if record.get("profile") is None:
        return ()
    labels: set[str] = set()
    for profile in record["profile"]["profiles"]:
        labels.update(str(value) for value in profile["individually_escaping_labels"])
    return tuple(sorted(labels))


def attempt_augmented_certificates(
    anchor: AnchorData,
    labelled_points: dict[str, RationalPoint],
    labels: Sequence[str],
    *,
    prime_bound: int,
    individual_limit: int,
    pair_limit: int,
) -> dict[str, Any]:
    """Promote quotient escape only after replayed full-column certificates."""

    attempts: list[dict[str, Any]] = []
    successful_individual: list[dict[str, Any]] = []
    successful_pairs: list[dict[str, Any]] = []
    selected_labels = tuple(labels[:individual_limit])

    def try_points(
        current_labels: tuple[str, ...], current_points: tuple[RationalPoint, ...]
    ) -> dict[str, Any] | None:
        augmented = anchor.basis + current_points
        expected_rank = len(augmented)
        for modulus in (2, 3, 5):
            started = time.monotonic()
            try:
                certificate = build_finite_quotient_certificate(
                    anchor.model,
                    augmented,
                    relation_prime=modulus,
                    prime_bound=prime_bound,
                )
            except Exception as error:  # preserve a failed backend as evidence
                attempts.append(
                    {
                        "labels": list(current_labels),
                        "modulus": modulus,
                        "status": "error",
                        "error_type": type(error).__name__,
                        "error": str(error),
                        "seconds": time.monotonic() - started,
                    }
                )
                continue
            rank = int(certificate["combined_rank_over_relation_field"])
            certified = bool(certificate["certified_independent"])
            attempts.append(
                {
                    "labels": list(current_labels),
                    "modulus": modulus,
                    "status": "completed",
                    "combined_rank": rank,
                    "expected_rank": expected_rank,
                    "certified_independent": certified,
                    "certificate_primes": certificate["certificate_primes"],
                    "seconds": time.monotonic() - started,
                }
            )
            print(
                "F20SEARCH|stage=certificate|"
                f"labels={','.join(current_labels)}|modulus={modulus}|"
                f"rank={rank}|expected={expected_rank}|certified={str(certified).lower()}",
                flush=True,
            )
            if certified and rank == expected_rank:
                verify_finite_quotient_certificate(
                    anchor.model, augmented, certificate
                )
                return {
                    "labels": list(current_labels),
                    "certified_rank_lower_bound": expected_rank,
                    "relation_modulus": modulus,
                    "certificate": certificate,
                }
        return None

    for label in selected_labels:
        success = try_points((label,), (labelled_points[label],))
        if success is not None:
            successful_individual.append(success)

    pair_attempts = 0
    for left_index, left in enumerate(selected_labels):
        for right in selected_labels[left_index + 1 :]:
            if pair_attempts >= pair_limit:
                break
            pair_attempts += 1
            success = try_points(
                (left, right),
                (labelled_points[left], labelled_points[right]),
            )
            if success is not None:
                successful_pairs.append(success)
                break
        if pair_attempts >= pair_limit or successful_pairs:
            break

    certified_rank = 20
    if successful_individual:
        certified_rank = max(
            certified_rank,
            max(item["certified_rank_lower_bound"] for item in successful_individual),
        )
    if successful_pairs:
        certified_rank = max(
            certified_rank,
            max(item["certified_rank_lower_bound"] for item in successful_pairs),
        )
    return {
        "attempted_labels": list(selected_labels),
        "individual_limit": individual_limit,
        "pair_attempt_limit": pair_limit,
        "attempts": attempts,
        "successful_individual_certificates": successful_individual,
        "successful_pair_certificates": successful_pairs,
        "certified_rank_lower_bound_after_search": certified_rank,
        "claim_boundary": (
            "rank is promoted only by a verified full-column finite-quotient "
            "infinite-descent certificate"
        ),
    }


def resolved_settings(args: argparse.Namespace) -> StageSettings:
    smoke = args.stage == "smoke"

    def choose(value: int | None, smoke_default: int, full_default: int) -> int:
        return value if value is not None else (smoke_default if smoke else full_default)

    shifts = (
        tuple(int(value) for value in args.original_center_shifts.split(","))
        if args.original_center_shifts is not None
        else ((0,) if smoke else (-1, 0, 1))
    )
    return StageSettings(
        cover_limit=choose(args.cover_limit, 8, 60),
        charts_per_cover=choose(args.charts_per_cover, 2, 3),
        original_box_limit=choose(args.original_box_limit, 2, len(ORIGINAL_SKEW_BOXES)),
        original_center_count=choose(args.original_center_count, 3, 12),
        original_center_shifts=shifts,
        pilot_height=choose(args.pilot_height, 50_000, 50_000),
        escalation_height=choose(args.escalation_height, 250_000, 250_000),
        escalation_charts=choose(args.escalation_charts, 0, 16),
        deep_height=choose(args.deep_height, 1_000_000, 1_000_000),
        deep_charts=choose(args.deep_charts, 0, 4),
        skew_cover_count=choose(args.skew_cover_count, 2, 8),
        skew_box_limit=choose(args.skew_box_limit, 1, 2),
    )


def load_completed_plan(
    path: Path,
    anchor: AnchorData,
    *,
    allow_bounded: bool,
) -> tuple[FrontierEntry, ...]:
    record = json.loads(path.read_text())
    allowed = {"complete"} | ({"bounded-complete"} if allow_bounded else set())
    if (
        record.get("schema") != PLAN_SCHEMA
        or record.get("basis_sha256") != anchor.basis_sha256
        or record.get("status") not in allowed
    ):
        raise ValueError("a completed matching frontier plan is required")
    return tuple(FrontierEntry.from_record(item) for item in record["selected"])


def candidate_sort_key(point: RationalPoint) -> tuple[int, int, RationalPoint]:
    return projective_height(point[0]), projective_height(point[1]), point


def reproducing_command() -> str:
    return " ".join(shlex.quote(value) for value in (sys.executable, *sys.argv))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        choices=("smoke", "frontier", "search", "all"),
        default="smoke",
    )
    parser.add_argument("--candidate-input", type=Path, default=DEFAULT_CANDIDATE_INPUT)
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--bank-cache", type=Path, default=DEFAULT_BANK)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--frontier-retain", type=int, default=64)
    parser.add_argument("--frontier-global", type=int, default=40)
    parser.add_argument("--frontier-maximum", type=int, default=60)
    parser.add_argument("--frontier-progress-interval", type=int, default=65_536)
    parser.add_argument(
        "--frontier-scan-limit",
        type=int,
        default=0,
        help="test-only Gray prefix; zero scans all represented classes",
    )
    parser.add_argument("--no-resume", action="store_true")
    parser.add_argument("--allow-bounded-plan", action="store_true")
    parser.add_argument("--smoke-maximum-weight", type=int, default=2)
    parser.add_argument("--cover-limit", type=int)
    parser.add_argument("--charts-per-cover", type=int)
    parser.add_argument("--original-box-limit", type=int)
    parser.add_argument("--original-center-count", type=int)
    parser.add_argument("--original-center-shifts")
    parser.add_argument("--pilot-height", type=int)
    parser.add_argument("--escalation-height", type=int)
    parser.add_argument("--escalation-charts", type=int)
    parser.add_argument("--deep-height", type=int)
    parser.add_argument("--deep-charts", type=int)
    parser.add_argument("--skew-cover-count", type=int)
    parser.add_argument("--skew-box-limit", type=int)
    parser.add_argument("--search-timeout", type=float, default=60.0)
    parser.add_argument("--relation-timeout", type=float, default=180.0)
    parser.add_argument("--relation-batch-size", type=int, default=32)
    parser.add_argument("--stack-bytes", type=int, default=1_000_000_000)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--skip-original", action="store_true")
    parser.add_argument("--skip-alternate", action="store_true")
    parser.add_argument("--skip-bank", action="store_true")
    parser.add_argument("--skip-certification", action="store_true")
    parser.add_argument("--bank-prime-bound", type=int, default=2_000)
    parser.add_argument("--bank-row-target", type=int, default=48)
    parser.add_argument("--bank-candidate-limit", type=int, default=64)
    parser.add_argument("--bank-batch-size", type=int, default=8)
    parser.add_argument("--certificate-prime-bound", type=int, default=2_000)
    parser.add_argument("--certificate-individual-limit", type=int, default=8)
    parser.add_argument("--certificate-pair-limit", type=int, default=12)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if (
        args.shard_count <= 0
        or not 0 <= args.shard_index < args.shard_count
        or args.search_timeout <= 0
        or args.relation_timeout <= 0
        or args.stack_bytes < 8_000_000
    ):
        raise SystemExit("invalid shard, timeout, or PARI stack configuration")

    anchor = load_anchor(args.candidate_input)
    print(
        "F20SEARCH|version=1|stage=input|"
        f"known_rank=20|basis_sha256={anchor.basis_sha256}|"
        f"known_signless_points={len(anchor.known_signless_points)}|"
        f"shard={args.shard_index}/{args.shard_count}",
        flush=True,
    )

    if args.stage in {"frontier", "all"}:
        entries = scan_frontier(
            anchor,
            args.plan,
            retain_count=args.frontier_retain,
            global_count=args.frontier_global,
            maximum_count=args.frontier_maximum,
            progress_interval=args.frontier_progress_interval,
            scan_limit=args.frontier_scan_limit,
            resume=not args.no_resume,
        )
        if args.stage == "frontier":
            return
    elif args.stage == "search":
        entries = load_completed_plan(
            args.plan,
            anchor,
            allow_bounded=args.allow_bounded_plan,
        )
    else:
        entries = low_weight_frontier(
            anchor,
            maximum_weight=args.smoke_maximum_weight,
            retain_count=args.cover_limit or 8,
        )

    settings = resolved_settings(args)
    plans = (
        ()
        if args.skip_alternate
        else build_cover_plans(
            anchor,
            entries,
            charts_per_cover=settings.charts_per_cover,
            cover_limit=settings.cover_limit,
        )
    )
    discoveries: dict[RationalPoint, set[str]] = {}
    original_runs = (
        []
        if args.skip_original
        else run_original_searches(
            anchor,
            settings,
            discoveries,
            timeout=args.search_timeout,
            stack_bytes=args.stack_bytes,
            shard_index=args.shard_index,
            shard_count=args.shard_count,
        )
    )
    alternate_runs = (
        []
        if args.skip_alternate
        else run_alternate_cover_searches(
            anchor,
            plans,
            settings,
            discoveries,
            timeout=args.search_timeout,
            stack_bytes=args.stack_bytes,
            shard_index=args.shard_index,
            shard_count=args.shard_count,
        )
    )

    candidates = tuple(sorted(discoveries, key=candidate_sort_key))
    proposals = relation_proposals(
        anchor.model,
        anchor.basis,
        candidates,
        timeout=args.relation_timeout,
        stack_bytes=args.stack_bytes,
        batch_size=args.relation_batch_size,
    )
    unresolved: list[RationalPoint] = []
    unresolved_labels: list[str] = []
    candidate_records: list[dict[str, Any]] = []
    for index, (point, proposal) in enumerate(zip(candidates, proposals, strict=True)):
        relation, exact, status = proposal
        label = f"candidate_{index:04d}"
        candidate_records.append(
            {
                "label": label,
                **point_record(point),
                "sources": sorted(discoveries[point]),
                "relation_process_status": status,
                "exact_relation_in_known_rank20_subgroup": exact,
                "basis_relation": list(relation) if exact and relation is not None else None,
                "fraction_group_law_replay": exact,
            }
        )
        if not exact:
            unresolved.append(point)
            unresolved_labels.append(label)

    bank_record: dict[str, Any] | None = None
    bank_status = "skipped"
    labelled_bank_points: dict[str, RationalPoint] = {}
    if unresolved and not args.skip_bank:
        bank_points = tuple(unresolved[: args.bank_candidate_limit])
        bank_labels = tuple(unresolved_labels[: args.bank_candidate_limit])
        labelled_bank_points = dict(zip(bank_labels, bank_points, strict=True))
        bank, bank_status = load_or_build_bank(
            anchor,
            args.bank_cache,
            prime_bound=args.bank_prime_bound,
            row_target=args.bank_row_target,
        )
        bank_record = evaluate_candidate_batches_against_bank(
            bank,
            anchor.basis,
            tuple(labelled_bank_points.values()),
            candidate_labels=tuple(labelled_bank_points),
            batch_size=args.bank_batch_size,
        )
        bank_record["candidate_inventory"] = [
            {"label": label, **point_record(point)}
            for label, point in labelled_bank_points.items()
        ]
        bank_record["cache"] = {
            "path": str(args.bank_cache),
            "status": bank_status,
            "baseline_ranks": bank.baseline_ranks,
        }

    escaping_labels = escaping_labels_from_profile(bank_record)
    certification = {
        "status": "skipped",
        "certified_rank_lower_bound_after_search": 20,
        "reason": "no exact quotient-escaping candidate was selected",
    }
    if (
        escaping_labels
        and labelled_bank_points
        and not args.skip_certification
    ):
        certification = {
            "status": "completed",
            **attempt_augmented_certificates(
                anchor,
                labelled_bank_points,
                escaping_labels,
                prime_bound=args.certificate_prime_bound,
                individual_limit=args.certificate_individual_limit,
                pair_limit=args.certificate_pair_limit,
            ),
        }

    certified_rank = int(
        certification.get("certified_rank_lower_bound_after_search", 20)
    )
    script_path = Path(__file__).resolve()
    output = args.output
    if args.shard_count > 1 and output == DEFAULT_OUTPUT:
        output = output.with_name(
            f"{output.stem}.shard_{args.shard_index:03d}_of_"
            f"{args.shard_count:03d}{output.suffix}"
        )
    artifact = {
        "schema": RESULT_SCHEMA,
        "status": "bounded direct Fermigier rank-20 direction search complete",
        "candidate": {
            "candidate_key": "fermigier-mestre-v1:u=28917/20",
            "tuple_parameter_T": fraction_text(ANCHOR_T),
            "short_model": [fraction_text(value) for value in anchor.model],
            "root_number": anchor.candidate_record["models"]["global_minimal"][
                "root_number"
            ],
            "certified_rank_lower_bound_before_search": 20,
            "basis_sha256": anchor.basis_sha256,
        },
        "scope": {
            "stage": args.stage,
            "shard_index": args.shard_index,
            "shard_count": args.shard_count,
            "frontier_entry_count": len(entries),
            "cover_plan_count": len(plans),
            "settings": {
                "cover_limit": settings.cover_limit,
                "charts_per_cover": settings.charts_per_cover,
                "original_box_limit": settings.original_box_limit,
                "original_center_count": settings.original_center_count,
                "original_center_shifts": list(settings.original_center_shifts),
                "pilot_height": settings.pilot_height,
                "escalation_height": settings.escalation_height,
                "escalation_charts": settings.escalation_charts,
                "deep_height": settings.deep_height,
                "deep_charts": settings.deep_charts,
                "skew_cover_count": settings.skew_cover_count,
                "skew_box_limit": settings.skew_box_limit,
                "timeout_seconds_per_search": args.search_timeout,
            },
            "bounded_search_not_complete": True,
        },
        "input": {
            "candidate_artifact": str(args.candidate_input),
            "candidate_artifact_sha256": hashlib.sha256(
                args.candidate_input.read_bytes()
            ).hexdigest(),
            "known_signless_decontamination_count": len(
                anchor.known_signless_points
            ),
            "frontier_plan": str(args.plan) if args.stage != "smoke" else None,
        },
        "cover_plans": [
            {
                **plan.frontier.to_record(),
                "id": plan.identifier,
                "base_point": point_record(plan.cover.base_point),
                "quartic_coefficients_ascending": [
                    fraction_text(value) for value in plan.cover.coefficients
                ],
                "charts": [
                    {
                        "basis_indices_one_based": [
                            index + 1 for index in chart.basis_indices
                        ],
                        "matrix_a_b_c_d": list(chart.matrix),
                        "mean_log10_known_height": chart.mean_log_height,
                        "median_log10_known_height": chart.median_log_height,
                        "maximum_log10_known_height": chart.maximum_log_height,
                    }
                    for chart in plan.charts
                ],
            }
            for plan in plans
        ],
        "search_runs": {
            "original": original_runs,
            "alternate": alternate_runs,
        },
        "point_triage": {
            "discovered_signless_points_outside_predeclared_set": len(candidates),
            "exactly_replayed_in_known_rank20_subgroup": sum(
                bool(record["exact_relation_in_known_rank20_subgroup"])
                for record in candidate_records
            ),
            "unresolved_after_relation_replay": len(unresolved),
            "candidates": candidate_records,
        },
        "overcomplete_quotient_bank": bank_record,
        "augmented_exact_certification": certification,
        "outcome": {
            "certified_rank_lower_bound_after_search": certified_rank,
            "rank_gain_certified": certified_rank > 20,
            "rank22_gain_certified": certified_rank > 21,
            "root_number_plus_one_parity_heuristic": (
                "even rank is expected conjecturally; no parity statement is used "
                "in the exact certificate"
            ),
        },
        "claim_boundary": (
            "bounded point-search absence is never promoted to a rank upper bound; "
            "only verified augmented finite-quotient certificates promote rank"
        ),
        "reproducing_command": reproducing_command(),
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
    }
    _atomic_json(output, artifact)
    print(
        "F20SEARCH|stage=done|status=complete|"
        f"candidates={len(candidates)}|unresolved={len(unresolved)}|"
        f"escaping={len(escaping_labels)}|certified_rank={certified_rank}|"
        f"output={output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
