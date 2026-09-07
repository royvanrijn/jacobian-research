#!/usr/bin/env python3
"""Exhaust the cheap multiple-root CRT classes in Fermigier's family.

The default local conditions deliberately force large powers in the
degree-20 discriminant factor ``H(T)`` at much lower congruence cost than a
naive modulus ``p^k`` would suggest.  Every one of the 144 local-symbol
combinations is combined by exact CRT.  In each resulting two-dimensional
lattice, a declared coefficient box in a Gauss-reduced basis is enumerated
and its shortest rational representatives are retained.

The lowest-height bounded pool is ranked by a good-reduction heuristic.
PARI/GP computes exact conductor and local-reduction data for a deterministic
mixture of the lowest-height and highest-score candidates.  Neither the score
nor an optional numerical height-pairing determinant implies a rank.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction
import hashlib
from itertools import product
import json
from math import gcd, prod
from pathlib import Path
import platform
import shlex
import sys
from typing import Any, Iterable

from crt_lattice import crt_pair, short_rational_representatives
from ek_k3 import fraction_mod, rational_to_string, valuation
from fermigier_mestre import (
    DISCRIMINANT_FACTOR_COEFFICIENTS,
    FermigierMestreFamily,
    ROOTS,
)
from multiple_root_lifting import (
    affine_variable_coefficients,
    fixed_divisor_valuation,
)
from pari_bridge import minimal_curve_data, pari_version
from search_record_residue_class import build_score_tables, score_rational


TARGET_LOG_CONDUCTOR = Decimal("182.72")


@dataclass(frozen=True)
class LocalConstraintGroup:
    """A union of congruence classes with one uniform local guarantee."""

    prime: int
    modulus: int
    residues: tuple[int, ...]
    forced_h_valuation: int
    reduction: str
    presented_model_scaling: int

    @property
    def expected_minimal_discriminant_valuation(self) -> int:
        return self.forced_h_valuation - 12 * self.presented_model_scaling


DEFAULT_GROUP_ORDER = (7, 11, 17, 19, 37)
CONSTRAINT_GROUPS = {
    7: LocalConstraintGroup(
        prime=7,
        modulus=49,
        residues=(7, 14, 21, 28, 35, 42),
        forced_h_valuation=18,
        reduction="split multiplicative after one minimalizing 7-scaling",
        presented_model_scaling=1,
    ),
    11: LocalConstraintGroup(
        prime=11,
        modulus=11,
        residues=(0, 5, 6),
        forced_h_valuation=4,
        reduction="split multiplicative",
        presented_model_scaling=0,
    ),
    17: LocalConstraintGroup(
        prime=17,
        modulus=17,
        residues=(1, 16),
        forced_h_valuation=4,
        reduction="split multiplicative",
        presented_model_scaling=0,
    ),
    19: LocalConstraintGroup(
        prime=19,
        modulus=19,
        residues=(5, 14),
        forced_h_valuation=4,
        reduction="split multiplicative",
        presented_model_scaling=0,
    ),
    37: LocalConstraintGroup(
        prime=37,
        modulus=37,
        residues=(4, 33),
        forced_h_valuation=3,
        reduction="split multiplicative",
        presented_model_scaling=0,
    ),
}


def parse_groups(value: str) -> tuple[int, ...]:
    """Parse a comma-separated subset of the implemented local groups."""

    if value.strip().lower() == "all":
        return DEFAULT_GROUP_ORDER
    if not value.strip():
        return ()
    try:
        return tuple(int(item) for item in value.split(","))
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "groups must be comma-separated primes"
        ) from error


def choose_groups(
    included: Iterable[int], omitted: Iterable[int]
) -> tuple[LocalConstraintGroup, ...]:
    """Return groups in canonical order after validating selection options."""

    included_tuple = tuple(included)
    omitted_tuple = tuple(omitted)
    if len(set(included_tuple)) != len(included_tuple):
        raise ValueError("--groups contains a duplicate prime")
    if len(set(omitted_tuple)) != len(omitted_tuple):
        raise ValueError("--omit-groups contains a duplicate prime")
    unknown = (set(included_tuple) | set(omitted_tuple)) - set(CONSTRAINT_GROUPS)
    if unknown:
        raise ValueError(f"unknown constraint groups: {sorted(unknown)}")
    selected = set(included_tuple) - set(omitted_tuple)
    if not selected:
        raise ValueError("at least one constraint group must remain")
    return tuple(
        CONSTRAINT_GROUPS[prime]
        for prime in DEFAULT_GROUP_ORDER
        if prime in selected
    )


def exact_group_certificate(group: LocalConstraintGroup) -> dict[str, Any]:
    """Recheck fixed divisors and clean split reduction for one group."""

    residue_certificates: list[dict[str, Any]] = []
    for residue in group.residues:
        fixed_valuation = fixed_divisor_valuation(
            affine_variable_coefficients(
                DISCRIMINANT_FACTOR_COEFFICIENTS,
                residue,
                group.modulus,
            ),
            group.prime,
        )
        if fixed_valuation < group.forced_h_valuation:
            raise AssertionError(
                f"T={residue} mod {group.modulus} lost its v_{group.prime}(H) "
                "guarantee"
            )

        parameter = Fraction(residue)
        invariants = FermigierMestreFamily.invariants(parameter)
        raw_c4_valuation = valuation(invariants["c4"], group.prime)
        raw_c6_valuation = valuation(invariants["c6"], group.prime)
        if group.presented_model_scaling:
            coefficients = FermigierMestreFamily.coefficients(parameter)
            scaling = group.prime**group.presented_model_scaling
            coefficient_a = coefficients[3] / scaling**4
            coefficient_b = coefficients[4] / scaling**6
            character_sum = 0
            for x_value in range(group.prime):
                rhs = (
                    x_value**3
                    + fraction_mod(coefficient_a, group.prime) * x_value
                    + fraction_mod(coefficient_b, group.prime)
                )
                symbol = pow(rhs % group.prime, (group.prime - 1) // 2, group.prime)
                character_sum += (
                    -1 if symbol == group.prime - 1 else symbol
                )
            trace = -character_sum
            if trace != 1:
                raise AssertionError("the scaled p=7 fiber is no longer split")
            minimal_c4_valuation = raw_c4_valuation - 4
        else:
            local = FermigierMestreFamily.local_data(
                residue % group.prime, group.prime
            )
            if local.split_multiplicative is not True:
                raise AssertionError(
                    f"the p={group.prime} residue is no longer clean split"
                )
            trace = local.trace
            minimal_c4_valuation = raw_c4_valuation
        if minimal_c4_valuation != 0:
            raise AssertionError("a declared clean class has nonunit minimal c4")

        residue_certificates.append(
            {
                "residue": residue,
                "modulus": group.modulus,
                "fixed_divisor_valuation": fixed_valuation,
                "presented_c4_valuation": raw_c4_valuation,
                "presented_c6_valuation": raw_c6_valuation,
                "minimal_c4_valuation": minimal_c4_valuation,
                "minimal_trace": trace,
            }
        )
    return {
        **asdict(group),
        "expected_minimal_discriminant_valuation_at_least": (
            group.expected_minimal_discriminant_valuation
        ),
        "union_density": rational_to_string(
            Fraction(len(group.residues), group.modulus)
        ),
        "reciprocal_union_density": rational_to_string(
            Fraction(group.modulus, len(group.residues))
        ),
        "residue_certificates": residue_certificates,
    }


def crt_classes(
    groups: tuple[LocalConstraintGroup, ...]
) -> tuple[dict[str, Any], ...]:
    """Return every independent local-symbol combination and its CRT class."""

    classes: list[dict[str, Any]] = []
    for class_index, residues in enumerate(
        product(*(group.residues for group in groups)), start=1
    ):
        crt_residue = 0
        crt_modulus = 1
        for group, residue in zip(groups, residues, strict=True):
            crt_residue, crt_modulus = crt_pair(
                crt_residue,
                crt_modulus,
                residue,
                group.modulus,
            )
        classes.append(
            {
                "class_index": class_index,
                "residues": {
                    str(group.prime): residue
                    for group, residue in zip(groups, residues, strict=True)
                },
                "crt_residue": crt_residue,
                "crt_modulus": crt_modulus,
            }
        )
    expected = prod(len(group.residues) for group in groups)
    if len(classes) != expected or len(
        {item["crt_residue"] for item in classes}
    ) != expected:
        raise AssertionError("CRT class enumeration is not exhaustive and unique")
    return tuple(classes)


def _negated_residues(
    residues: dict[str, int], groups: tuple[LocalConstraintGroup, ...]
) -> dict[str, int]:
    return {
        str(group.prime): (-residues[str(group.prime)]) % group.modulus
        for group in groups
    }


def _exact_candidate(
    numerator: int,
    denominator: int,
    groups: tuple[LocalConstraintGroup, ...],
    path: dict[str, Any],
) -> dict[str, Any] | None:
    parameter = Fraction(numerator, denominator)
    discriminant_factor = FermigierMestreFamily.discriminant_factor(parameter)
    if discriminant_factor == 0:
        return None
    actual_valuations: dict[str, int] = {}
    checked_residues: dict[str, int] = {}
    for group in groups:
        if gcd(denominator, group.modulus) != 1:
            raise AssertionError("a constrained denominator is not a local unit")
        residue = fraction_mod(parameter, group.modulus)
        expected_residue = path["canonical_residues"][str(group.prime)]
        if residue != expected_residue:
            raise AssertionError("a rational representative lost its CRT class")
        if residue not in group.residues:
            raise AssertionError("sign normalization left the allowed residue union")
        checked_residues[str(group.prime)] = residue
        actual = valuation(discriminant_factor, group.prime)
        if actual < group.forced_h_valuation:
            raise AssertionError(
                f"candidate lost the v_{group.prime}(H) guarantee"
            )
        actual_valuations[str(group.prime)] = actual
    return {
        "t": rational_to_string(parameter),
        "numerator": numerator,
        "denominator": denominator,
        "height": max(abs(numerator), denominator),
        "residues": checked_residues,
        "h_valuations": actual_valuations,
        "forcing_paths": [path],
    }


def enumerate_candidates(
    classes: tuple[dict[str, Any], ...],
    groups: tuple[LocalConstraintGroup, ...],
    *,
    coefficient_radius: int,
    representatives_per_class: int,
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    """Enumerate the bounded lattice neighborhoods and deduplicate ``T <-> -T``."""

    records: dict[tuple[int, int], dict[str, Any]] = {}
    bounded_representatives = 0
    retained_representatives = 0
    singular_rejections = 0
    symmetry_or_duplicate_merges = 0
    maximum_in_box = (2 * coefficient_radius + 1) ** 2 - 1
    for crt_class in classes:
        all_in_box = short_rational_representatives(
            crt_class["crt_residue"],
            crt_class["crt_modulus"],
            coefficient_radius=coefficient_radius,
            limit=maximum_in_box,
        )
        bounded_representatives += len(all_in_box)
        for representative in all_in_box[:representatives_per_class]:
            retained_representatives += 1
            numerator = representative.numerator
            canonical_residue = crt_class["crt_residue"]
            canonical_residues = dict(crt_class["residues"])
            sign_normalized = numerator < 0
            if sign_normalized:
                numerator = -numerator
                canonical_residue = (-canonical_residue) % crt_class["crt_modulus"]
                canonical_residues = _negated_residues(canonical_residues, groups)
            path = {
                "source_class_index": crt_class["class_index"],
                "source_crt_residue": crt_class["crt_residue"],
                "canonical_crt_residue": canonical_residue,
                "crt_modulus": crt_class["crt_modulus"],
                "source_residues": crt_class["residues"],
                "canonical_residues": canonical_residues,
                "sign_normalized": sign_normalized,
            }
            if (
                numerator - canonical_residue * representative.denominator
            ) % crt_class["crt_modulus"]:
                raise AssertionError("sign-normalized CRT identity failed")
            key = (numerator, representative.denominator)
            if key in records:
                symmetry_or_duplicate_merges += 1
                if path not in records[key]["forcing_paths"]:
                    records[key]["forcing_paths"].append(path)
                continue
            record = _exact_candidate(
                numerator,
                representative.denominator,
                groups,
                path,
            )
            if record is None:
                singular_rejections += 1
                continue
            records[key] = record

    ordered = sorted(
        records.values(),
        key=lambda item: (
            item["height"],
            item["numerator"],
            item["denominator"],
        ),
    )
    for rank, record in enumerate(ordered, start=1):
        record["height_rank"] = rank
    counts = {
        "crt_classes_visited": len(classes),
        "bounded_representatives_in_reduced_basis_boxes": bounded_representatives,
        "short_representatives_retained_before_deduplication": (
            retained_representatives
        ),
        "singular_representatives_rejected": singular_rejections,
        "sign_symmetry_or_duplicate_merges": symmetry_or_duplicate_merges,
        "unique_nonsingular_representatives": len(ordered),
    }
    return ordered, counts


def score_height_pool(
    height_ranked: list[dict[str, Any]],
    *,
    height_pool: int,
    score_bound: int,
    score: str,
) -> list[dict[str, Any]]:
    """Score exactly the declared lowest-height prefix."""

    tables = build_score_tables(score_bound, score)
    pool = height_ranked[:height_pool]
    for record in pool:
        record["score"] = score_rational(
            record["numerator"], record["denominator"], tables
        )
    ranked = sorted(
        pool,
        key=lambda item: (
            -item["score"]["value"],
            item["height"],
            item["numerator"],
            item["denominator"],
        ),
    )
    for rank, record in enumerate(ranked, start=1):
        record["score_rank_within_height_pool"] = rank
    return ranked


def select_for_pari(
    score_ranked: list[dict[str, Any]], count: int
) -> list[tuple[dict[str, Any], str]]:
    """Alternate deterministically between lowest height and highest score."""

    if count <= 0:
        return []
    height_ranked = sorted(score_ranked, key=lambda item: item["height_rank"])
    sources = ((height_ranked, "lowest height"), (score_ranked, "highest score"))
    selected: list[tuple[dict[str, Any], str]] = []
    seen: set[tuple[int, int]] = set()
    index = 0
    while len(selected) < min(count, len(score_ranked)):
        progressed = False
        for records, reason in sources:
            if index >= len(records):
                continue
            record = records[index]
            key = (record["numerator"], record["denominator"])
            if key in seen:
                continue
            selected.append((record, reason))
            seen.add(key)
            progressed = True
            if len(selected) == count:
                break
        index += 1
        if not progressed and index >= len(score_ranked):
            break
    return selected


def retained_records(
    score_ranked: list[dict[str, Any]], keep: int
) -> list[dict[str, Any]]:
    """Retain the union of the best height and score prefixes."""

    height_ranked = sorted(score_ranked, key=lambda item: item["height_rank"])
    retained_keys = {
        (record["numerator"], record["denominator"])
        for record in score_ranked[:keep] + height_ranked[:keep]
    }
    return [
        record
        for record in score_ranked
        if (record["numerator"], record["denominator"]) in retained_keys
    ]


def verify_pari_local_data(
    record: dict[str, Any], groups: tuple[LocalConstraintGroup, ...]
) -> None:
    """Check PARI's minimal local data against every engineered condition."""

    reductions = record["pari"]["local_reduction"]
    checks: dict[str, bool] = {}
    for group in groups:
        local = reductions[str(group.prime)]
        expected_discriminant = (
            record["h_valuations"][str(group.prime)]
            - 12 * group.presented_model_scaling
        )
        valid = (
            local["conductor_exponent"] == 1
            and local["minimal_c4_valuation"] == 0
            and local["minimal_discriminant_valuation"] == expected_discriminant
            and local["ellap"] == 1
        )
        checks[str(group.prime)] = valid
        if not valid:
            raise AssertionError(
                f"PARI local data contradicted the p={group.prime} certificate"
            )
    record["pari_local_constraints_verified"] = checks


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--groups", type=parse_groups, default=DEFAULT_GROUP_ORDER)
    parser.add_argument("--omit-groups", type=parse_groups, default=())
    parser.add_argument("--coefficient-radius", type=int, default=12)
    parser.add_argument("--representatives-per-class", type=int, default=12)
    parser.add_argument(
        "--height-pool",
        type=int,
        default=512,
        help="number of lowest-height representatives admitted to score ranking",
    )
    parser.add_argument("--score-bound", type=int, default=200)
    parser.add_argument(
        "--score",
        choices=("fermigier-good", "nagao-log"),
        default="fermigier-good",
    )
    parser.add_argument("--keep", type=int, default=24)
    parser.add_argument("--pari-count", type=int, default=4)
    parser.add_argument("--pari-timeout", type=float, default=45.0)
    parser.add_argument("--pari-stack-bytes", type=int, default=256_000_000)
    parser.add_argument(
        "--point-determinant-count",
        type=int,
        default=1,
        help="supply twelve mapped sections for this many PARI calls",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            root
            / "artifacts"
            / "generated-results"
            / "elliptic_fermigier_multiple_root_crt.json"
        ),
    )
    return parser


def validate_args(args: argparse.Namespace) -> tuple[LocalConstraintGroup, ...]:
    try:
        groups = choose_groups(args.groups, args.omit_groups)
    except ValueError as error:
        raise SystemExit(str(error)) from error
    positive = {
        "--coefficient-radius": args.coefficient_radius,
        "--representatives-per-class": args.representatives_per_class,
        "--height-pool": args.height_pool,
        "--score-bound": args.score_bound,
        "--keep": args.keep,
        "--pari-timeout": args.pari_timeout,
        "--pari-stack-bytes": args.pari_stack_bytes,
    }
    for name, value in positive.items():
        if value <= 0:
            raise SystemExit(f"{name} must be positive")
    if args.score_bound < 5:
        raise SystemExit("--score-bound must be at least 5")
    if args.pari_count < 0 or args.point_determinant_count < 0:
        raise SystemExit("PARI candidate counts must be nonnegative")
    if args.point_determinant_count > args.pari_count:
        raise SystemExit("--point-determinant-count cannot exceed --pari-count")
    return groups


def main() -> None:
    args = build_parser().parse_args()
    groups = validate_args(args)
    group_certificates = [exact_group_certificate(group) for group in groups]
    classes = crt_classes(groups)
    expected_class_count = prod(len(group.residues) for group in groups)
    height_ranked, enumeration = enumerate_candidates(
        classes,
        groups,
        coefficient_radius=args.coefficient_radius,
        representatives_per_class=args.representatives_per_class,
    )
    if not height_ranked:
        raise SystemExit("the bounded lattice enumeration found no candidate")
    score_ranked = score_height_pool(
        height_ranked,
        height_pool=args.height_pool,
        score_bound=args.score_bound,
        score=args.score,
    )

    pari_errors: list[dict[str, str]] = []
    pari_selected = select_for_pari(score_ranked, args.pari_count)
    for selected_index, (record, reason) in enumerate(pari_selected):
        parameter = Fraction(record["numerator"], record["denominator"])
        known_points = None
        if selected_index < args.point_determinant_count:
            known_points = FermigierMestreFamily.known_jacobian_points(parameter)[1:]
        record["pari_selection_reason"] = reason
        try:
            record["pari"] = minimal_curve_data(
                FermigierMestreFamily.coefficients(parameter),
                timeout=args.pari_timeout,
                known_points=known_points,
                local_primes=tuple(group.prime for group in groups),
                stack_bytes=args.pari_stack_bytes,
            )
            record["below_log_conductor_target"] = (
                Decimal(record["pari"]["log_conductor"])
                < TARGET_LOG_CONDUCTOR
            )
            record["target_status"] = (
                "conductor inequality only; rank was not computed"
                if record["below_log_conductor_target"]
                else "not a target hit; rank was not computed"
            )
            verify_pari_local_data(record, groups)
        except Exception as error:
            pari_errors.append({"t": record["t"], "error": str(error)})

    retained = retained_records(score_ranked, args.keep)
    retained_keys = {
        (record["numerator"], record["denominator"]) for record in retained
    }
    for record, _ in pari_selected:
        key = (record["numerator"], record["denominator"])
        if key not in retained_keys:
            retained.append(record)
            retained_keys.add(key)
    retained.sort(
        key=lambda item: (
            item["score_rank_within_height_pool"],
            item["height_rank"],
        )
    )

    score_tables = build_score_tables(args.score_bound, args.score)
    for record in retained:
        record["score"] = score_rational(
            record["numerator"],
            record["denominator"],
            score_tables,
            include_traces=True,
        )

    completed = [record for record, _ in pari_selected if "pari" in record]
    completed.sort(key=lambda item: Decimal(item["pari"]["log_conductor"]))
    conductor_summary = {
        "calls_requested": args.pari_count,
        "calls_completed": len(completed),
        "below_log_conductor_bound": sum(
            record["below_log_conductor_target"] for record in completed
        ),
        "best": (
            {
                "t": completed[0]["t"],
                "height": completed[0]["height"],
                "log_conductor": completed[0]["pari"]["log_conductor"],
            }
            if completed
            else None
        ),
    }

    command = " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv])
    script_path = Path(__file__).resolve()
    artifact = {
        "schema_version": 1,
        "status": (
            "bounded exhaustive-CRT experiment; the lattice neighborhoods and "
            "score pool are bounded, scores are heuristic, and no rank is "
            "computed or inferred"
        ),
        "target": {
            "rank_at_least": 21,
            "log_conductor_strict_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "alternative_rank_at_least": 30,
            "hits": [],
            "explanation": (
                "conductor-only successes are not target hits without certified "
                "independent rational points"
            ),
        },
        "family": {
            "name": "normalized Fermigier--Mestre rank-at-least-12 family",
            "root_tuple": list(ROOTS),
            "source": "https://matwbn.icm.edu.pl/ksiazki/aa/aa82/aa8243.pdf",
        },
        "constraint_groups": group_certificates,
        "method": {
            "crt": (
                "all local-symbol combinations; exact CRT with pairwise-coprime "
                "group moduli"
            ),
            "lattice_neighborhood": (
                "all coefficient pairs in [-coefficient_radius,coefficient_radius]^2 "
                "in an exact Gauss-reduced lattice basis, followed by the configured "
                "shortest-per-class cutoff"
            ),
            "parameter_symmetry": (
                "T and -T give the same Jacobian coefficients and are deduplicated; "
                "every source CRT path is retained"
            ),
            "selection_order": (
                "exact local verification for every short representative; score only "
                "the lowest-height pool; retain the union of height and score prefixes"
            ),
            "score": args.score,
            "score_scope": (
                "good reduction only at numerical primes 5<=p<=score_bound; "
                "denominator and bad-reduction primes are excluded"
            ),
            "pari_subset": (
                "deterministic alternation between lowest height and highest score; "
                "conductor and local reductions only, never ellrank"
            ),
        },
        "parameters": {
            key: (list(value) if isinstance(value, tuple) else str(value))
            if isinstance(value, (tuple, Path))
            else value
            for key, value in vars(args).items()
        },
        "enumeration": {
            "expected_crt_classes": expected_class_count,
            "crt_modulus": classes[0]["crt_modulus"],
            **enumeration,
            "height_pool_scored": len(score_ranked),
            "records_output": len(retained),
        },
        "conductor_summary": conductor_summary,
        "pari_errors": pari_errors,
        "software": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "pari_gp": pari_version(),
        },
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "reproducing_command": command,
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
        "candidates": retained,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")

    print(f"wrote {args.output}")
    print(
        f"classes={len(classes)}/{expected_class_count} "
        f"unique={enumeration['unique_nonsingular_representatives']} "
        f"height_pool={len(score_ranked)} retained={len(retained)}"
    )
    for record, reason in pari_selected:
        summary = (
            f"t={record['t']} height={record['height']} "
            f"height_rank={record['height_rank']} "
            f"score_rank={record['score_rank_within_height_pool']} "
            f"selected={reason!r}"
        )
        if "pari" in record:
            summary += f" logN={record['pari']['log_conductor']}"
        else:
            summary += " logN=unavailable"
        print(summary)


if __name__ == "__main__":
    main()
