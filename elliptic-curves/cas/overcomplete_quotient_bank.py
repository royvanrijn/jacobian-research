#!/usr/bin/env python3
"""Reusable overcomplete finite-quotient banks for elliptic-curve searches.

A compact independence certificate normally stops as soon as the known
columns have full rank.  That is ideal for certifying a basis, but a matrix
with exactly ``r`` independent rows cannot expose a new ``(r+1)``-st column.
This module deliberately retains additional good-reduction quotient rows.

The resulting bank is a search filter, not by itself a Mordell--Weil rank
certificate.  A positive marginal dimension is exact finite-quotient
evidence that a candidate escapes the image of the known subgroup in the
retained product.  Final rank claims should still be replayed with
``build_finite_quotient_certificate`` on the augmented point set.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Iterable, Sequence

from elliptic_candidate_record import (
    finite_quotient_signature,
    fraction_text,
    matrix_rank_and_pivots_mod_prime,
    model_from_record,
    point_sequence_sha256,
    primes_up_to,
)
from finite_quotient_escape import (
    MultiModulusEscapeProfile,
    QuotientBlock,
    analyze_multi_modulus_escape,
)


Q = Fraction
RationalPoint = tuple[Fraction, Fraction]
SCHEMA = "elliptic-curves.overcomplete-quotient-bank.v1"


@dataclass(frozen=True)
class QuotientBankEntry:
    modulus: int
    prime: int
    group_order: int
    multiple_subgroup_order: int
    quotient_dimension: int
    known_rows: tuple[tuple[int, ...], ...]

    @classmethod
    def from_signature(
        cls, modulus: int, signature: dict[str, Any]
    ) -> "QuotientBankEntry":
        rows = tuple(tuple(int(value) for value in row) for row in signature["rows"])
        quotient_dimension = int(signature["quotient_dimension"])
        if quotient_dimension <= 0 or len(rows) != quotient_dimension:
            raise ValueError("only nontrivial finite quotient signatures belong in a bank")
        return cls(
            modulus=int(modulus),
            prime=int(signature["prime"]),
            group_order=int(signature["group_order"]),
            multiple_subgroup_order=int(signature["multiple_subgroup_order"]),
            quotient_dimension=quotient_dimension,
            known_rows=rows,
        )

    def to_record(self) -> dict[str, Any]:
        return {
            "modulus": self.modulus,
            "prime": self.prime,
            "group_order": self.group_order,
            "multiple_subgroup_order": self.multiple_subgroup_order,
            "quotient_dimension": self.quotient_dimension,
            "known_rows": [list(row) for row in self.known_rows],
        }

    @classmethod
    def from_record(cls, record: dict[str, Any]) -> "QuotientBankEntry":
        answer = cls(
            modulus=int(record["modulus"]),
            prime=int(record["prime"]),
            group_order=int(record["group_order"]),
            multiple_subgroup_order=int(record["multiple_subgroup_order"]),
            quotient_dimension=int(record["quotient_dimension"]),
            known_rows=tuple(
                tuple(int(value) for value in row) for row in record["known_rows"]
            ),
        )
        if answer.quotient_dimension <= 0 or len(answer.known_rows) != answer.quotient_dimension:
            raise ValueError("a serialized quotient-bank entry is malformed")
        return answer


@dataclass(frozen=True)
class ModulusBank:
    modulus: int
    entries: tuple[QuotientBankEntry, ...]
    baseline_rank: int
    row_count: int
    scanned_prime_bound: int
    stopped_after_row_target: bool
    no_torsion_witness: tuple[int, int] | None

    def __post_init__(self) -> None:
        if any(entry.modulus != self.modulus for entry in self.entries):
            raise ValueError("a modulus bank mixed relation moduli")
        if self.row_count != sum(entry.quotient_dimension for entry in self.entries):
            raise ValueError("the stored quotient row count is inconsistent")

    def to_record(self) -> dict[str, Any]:
        return {
            "modulus": self.modulus,
            "baseline_rank": self.baseline_rank,
            "row_count": self.row_count,
            "entry_count": len(self.entries),
            "scanned_prime_bound": self.scanned_prime_bound,
            "stopped_after_row_target": self.stopped_after_row_target,
            "no_rational_modulus_torsion_witness": (
                None
                if self.no_torsion_witness is None
                else {
                    "prime": self.no_torsion_witness[0],
                    "group_order": self.no_torsion_witness[1],
                }
            ),
            "entries": [entry.to_record() for entry in self.entries],
        }

    @classmethod
    def from_record(cls, record: dict[str, Any]) -> "ModulusBank":
        witness_record = record.get("no_rational_modulus_torsion_witness")
        witness = (
            None
            if witness_record is None
            else (int(witness_record["prime"]), int(witness_record["group_order"]))
        )
        return cls(
            modulus=int(record["modulus"]),
            entries=tuple(
                QuotientBankEntry.from_record(entry) for entry in record["entries"]
            ),
            baseline_rank=int(record["baseline_rank"]),
            row_count=int(record["row_count"]),
            scanned_prime_bound=int(record["scanned_prime_bound"]),
            stopped_after_row_target=bool(record["stopped_after_row_target"]),
            no_torsion_witness=witness,
        )


@dataclass(frozen=True)
class OvercompleteQuotientBank:
    model: tuple[Fraction, Fraction, Fraction, Fraction, Fraction]
    known_point_count: int
    known_point_sequence_sha256: str
    moduli: tuple[ModulusBank, ...]
    prime_bound: int
    row_target_per_modulus: int

    def __post_init__(self) -> None:
        if len(self.model) != 5 or any(self.model[:3]):
            raise ValueError("the quotient bank currently requires a short model")
        if self.known_point_count <= 0:
            raise ValueError("a quotient bank needs a nonempty known point set")
        modulus_values = tuple(bank.modulus for bank in self.moduli)
        if len(set(modulus_values)) != len(modulus_values):
            raise ValueError("the quotient bank repeated a relation modulus")
        for bank in self.moduli:
            if any(
                len(row) != self.known_point_count
                for entry in bank.entries
                for row in entry.known_rows
            ):
                raise ValueError("a quotient-bank row has the wrong known width")
            rows = [row for entry in bank.entries for row in entry.known_rows]
            replay_rank, _ = matrix_rank_and_pivots_mod_prime(
                rows, self.known_point_count, bank.modulus
            )
            if replay_rank != bank.baseline_rank:
                raise ValueError("a quotient-bank baseline rank is inconsistent")
            if (
                bank.no_torsion_witness is not None
                and bank.no_torsion_witness[1] % bank.modulus == 0
            ):
                raise ValueError("the stored torsion witness is invalid")

    def to_record(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "claim_boundary": (
                "positive marginal quotient rank is exact escape evidence; "
                "final Mordell--Weil rank claims require an augmented full certificate"
            ),
            "model_coefficients_a1_a2_a3_a4_a6": [
                fraction_text(value) for value in self.model
            ],
            "known_point_count": self.known_point_count,
            "known_point_sequence_sha256": self.known_point_sequence_sha256,
            "prime_bound": self.prime_bound,
            "row_target_per_modulus": self.row_target_per_modulus,
            "moduli": [bank.to_record() for bank in self.moduli],
        }

    @classmethod
    def from_record(
        cls,
        record: dict[str, Any],
        *,
        model: Sequence[Fraction | int],
        known_points: Sequence[RationalPoint],
    ) -> "OvercompleteQuotientBank":
        if record.get("schema") != SCHEMA:
            raise ValueError("unknown overcomplete quotient-bank schema")
        expected_model = model_from_record(model)
        recorded_model = model_from_record(
            record["model_coefficients_a1_a2_a3_a4_a6"]
        )
        if recorded_model != expected_model:
            raise ValueError("the quotient-bank curve model does not match")
        expected_hash = point_sequence_sha256(tuple(known_points))
        if (
            int(record["known_point_count"]) != len(known_points)
            or record["known_point_sequence_sha256"] != expected_hash
        ):
            raise ValueError("the quotient-bank known basis does not match")
        return cls(
            model=expected_model,
            known_point_count=len(known_points),
            known_point_sequence_sha256=expected_hash,
            moduli=tuple(ModulusBank.from_record(item) for item in record["moduli"]),
            prime_bound=int(record["prime_bound"]),
            row_target_per_modulus=int(record["row_target_per_modulus"]),
        )

    @property
    def baseline_ranks(self) -> dict[int, int]:
        return {bank.modulus: bank.baseline_rank for bank in self.moduli}


def _combined_rank(
    entries: Iterable[QuotientBankEntry], column_count: int, modulus: int
) -> int:
    rows = [row for entry in entries for row in entry.known_rows]
    rank, _ = matrix_rank_and_pivots_mod_prime(rows, column_count, modulus)
    return rank


def build_overcomplete_quotient_bank(
    model: Sequence[Fraction | int],
    known_points: Sequence[RationalPoint],
    *,
    moduli: Sequence[int] = (2, 3, 5),
    prime_bound: int = 2_000,
    row_target_per_modulus: int = 48,
    require_full_baseline_rank: bool = True,
    progress: bool = False,
) -> OvercompleteQuotientBank:
    """Build deterministic extra-row finite-quotient banks.

    For each relation modulus, every nontrivial quotient encountered is kept
    until both the known columns have full rank and the requested row target
    has been reached.  The retained rows are intentionally not greedily
    minimized.
    """

    coefficients = model_from_record(model)
    points = tuple((Q(x), Q(y)) for x, y in known_points)
    if len(coefficients) != 5 or any(coefficients[:3]):
        raise ValueError("the quotient bank currently requires a short model")
    if not points:
        raise ValueError("at least one known point is required")
    if prime_bound < 3 or row_target_per_modulus <= 0:
        raise ValueError("prime and row bounds must be positive")
    relation_moduli = tuple(int(value) for value in moduli)
    if not relation_moduli or len(set(relation_moduli)) != len(relation_moduli):
        raise ValueError("provide distinct relation moduli")

    banks: list[ModulusBank] = []
    all_primes = primes_up_to(prime_bound)
    for modulus in relation_moduli:
        entries: list[QuotientBankEntry] = []
        row_count = 0
        rank = 0
        witness: tuple[int, int] | None = None
        last_scanned = 2
        stopped = False
        for reduction_prime in all_primes:
            last_scanned = reduction_prime
            if reduction_prime in (2, modulus):
                continue
            try:
                signature = finite_quotient_signature(
                    coefficients, points, reduction_prime, modulus
                )
            except ValueError:
                continue
            group_order = int(signature["group_order"])
            if witness is None and group_order % modulus:
                witness = (reduction_prime, group_order)
            if int(signature["quotient_dimension"]) <= 0:
                continue
            entry = QuotientBankEntry.from_signature(modulus, signature)
            entries.append(entry)
            row_count += entry.quotient_dimension
            rank = _combined_rank(entries, len(points), modulus)
            if progress:
                print(
                    "QTBANK|"
                    f"modulus={modulus}|prime={reduction_prime}|rows={row_count}|"
                    f"baseline_rank={rank}",
                    flush=True,
                )
            if rank == len(points) and row_count >= row_target_per_modulus:
                stopped = True
                break
        if require_full_baseline_rank and rank != len(points):
            raise RuntimeError(
                f"mod-{modulus} quotient bank reached rank {rank}, expected "
                f"{len(points)}, through reduction prime {last_scanned}"
            )
        banks.append(
            ModulusBank(
                modulus=modulus,
                entries=tuple(entries),
                baseline_rank=rank,
                row_count=row_count,
                scanned_prime_bound=last_scanned,
                stopped_after_row_target=stopped,
                no_torsion_witness=witness,
            )
        )

    return OvercompleteQuotientBank(
        model=coefficients,
        known_point_count=len(points),
        known_point_sequence_sha256=point_sequence_sha256(points),
        moduli=tuple(banks),
        prime_bound=prime_bound,
        row_target_per_modulus=row_target_per_modulus,
    )


def save_overcomplete_quotient_bank(
    bank: OvercompleteQuotientBank, path: Path
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(bank.to_record(), indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def load_overcomplete_quotient_bank(
    path: Path,
    *,
    model: Sequence[Fraction | int],
    known_points: Sequence[RationalPoint],
) -> OvercompleteQuotientBank:
    return OvercompleteQuotientBank.from_record(
        json.loads(path.read_text()), model=model, known_points=known_points
    )


def evaluate_candidates_against_bank(
    bank: OvercompleteQuotientBank,
    known_points: Sequence[RationalPoint],
    candidates: Sequence[RationalPoint],
    *,
    candidate_labels: Sequence[str] | None = None,
) -> dict[str, Any]:
    """Replay candidate columns in every usable retained quotient.

    A reduction prime is skipped when a candidate denominator is not
    invertible there.  The exact known-column prefix is checked against the
    cached bank before the resulting block is admitted.
    """

    known = tuple((Q(x), Q(y)) for x, y in known_points)
    extra = tuple((Q(x), Q(y)) for x, y in candidates)
    if len(known) != bank.known_point_count:
        raise ValueError("the supplied known point count differs from the bank")
    if point_sequence_sha256(known) != bank.known_point_sequence_sha256:
        raise ValueError("the supplied known point ordering differs from the bank")
    labels = (
        tuple(str(label) for label in candidate_labels)
        if candidate_labels is not None
        else tuple(f"candidate_{index:04d}" for index in range(len(extra)))
    )
    if len(labels) != len(extra) or len(set(labels)) != len(labels):
        raise ValueError("candidate labels must be distinct and match the candidates")
    if not extra:
        return {
            "candidate_count": 0,
            "candidate_labels": [],
            "profile": None,
            "availability": {},
            "claim_boundary": "no candidates were supplied",
        }

    blocks: list[QuotientBlock] = []
    availability: dict[str, dict[str, int]] = {}
    all_points = known + extra
    for modulus_bank in bank.moduli:
        used = 0
        skipped = 0
        rows = 0
        for entry in modulus_bank.entries:
            try:
                signature = finite_quotient_signature(
                    bank.model, all_points, entry.prime, entry.modulus
                )
            except ValueError:
                skipped += 1
                continue
            replay_rows = tuple(
                tuple(int(value) for value in row) for row in signature["rows"]
            )
            known_prefix = tuple(
                tuple(row[: len(known)]) for row in replay_rows
            )
            if known_prefix != entry.known_rows:
                raise AssertionError(
                    f"cached quotient prefix changed at p={entry.prime}, "
                    f"ell={entry.modulus}"
                )
            blocks.append(
                QuotientBlock.build(
                    modulus=entry.modulus,
                    rows=replay_rows,
                    column_count=len(all_points),
                    source=f"p={entry.prime}",
                )
            )
            used += 1
            rows += len(replay_rows)
        availability[str(modulus_bank.modulus)] = {
            "retained_entries": len(modulus_bank.entries),
            "usable_entries": used,
            "skipped_for_candidate_denominators": skipped,
            "usable_rows": rows,
        }
    if not blocks:
        raise RuntimeError("no retained quotient prime was usable for these candidates")

    profile: MultiModulusEscapeProfile = analyze_multi_modulus_escape(
        blocks,
        known_column_count=len(known),
        candidate_labels=labels,
    )
    return {
        "candidate_count": len(extra),
        "candidate_labels": list(labels),
        "profile": profile.to_record(),
        "availability": availability,
        "claim_boundary": (
            "escape is exact finite-quotient evidence; non-escape proves nothing, "
            "and final rank promotion requires an augmented full certificate"
        ),
    }


def evaluate_candidate_batches_against_bank(
    bank: OvercompleteQuotientBank,
    known_points: Sequence[RationalPoint],
    candidates: Sequence[RationalPoint],
    *,
    candidate_labels: Sequence[str] | None = None,
    batch_size: int = 8,
) -> dict[str, Any]:
    """Evaluate bounded candidate batches to avoid denominator-union collapse.

    ``finite_quotient_signature`` needs every supplied affine denominator to be
    invertible at the reduction prime.  Passing a large heterogeneous point
    set at once can therefore discard most of an otherwise useful bank.  This
    wrapper preserves exact within-batch marginal-rank analysis while keeping
    that union small.  Cross-batch independence is intentionally left to the
    final augmented certificate stage.
    """

    if batch_size <= 0:
        raise ValueError("the quotient-bank candidate batch size must be positive")
    extra = tuple((Q(x), Q(y)) for x, y in candidates)
    labels = (
        tuple(str(label) for label in candidate_labels)
        if candidate_labels is not None
        else tuple(f"candidate_{index:04d}" for index in range(len(extra)))
    )
    if len(labels) != len(extra) or len(set(labels)) != len(labels):
        raise ValueError("candidate labels must be distinct and match the candidates")
    batches: list[dict[str, Any]] = []
    escaping: set[str] = set()
    for start in range(0, len(extra), batch_size):
        current_points = extra[start : start + batch_size]
        current_labels = labels[start : start + batch_size]
        try:
            result = evaluate_candidates_against_bank(
                bank,
                known_points,
                current_points,
                candidate_labels=current_labels,
            )
        except RuntimeError as error:
            result = {
                "candidate_count": len(current_points),
                "candidate_labels": list(current_labels),
                "profile": None,
                "availability": {},
                "status": "no_usable_retained_prime",
                "error": str(error),
            }
        result["batch_index"] = len(batches)
        result["candidate_offset"] = start
        profile = result.get("profile")
        if profile is not None:
            for modulus_profile in profile["profiles"]:
                escaping.update(
                    str(label)
                    for label in modulus_profile["individually_escaping_labels"]
                )
        batches.append(result)
    return {
        "evaluation_mode": "bounded-candidate-batches",
        "candidate_count": len(extra),
        "candidate_labels": list(labels),
        "batch_size": batch_size,
        "batch_count": len(batches),
        "batches": batches,
        "escaping_labels": sorted(escaping),
        "claim_boundary": (
            "within-batch quotient escape is exact; cross-batch independence is "
            "not inferred and requires an augmented full certificate"
        ),
    }
