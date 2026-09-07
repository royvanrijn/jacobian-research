#!/usr/bin/env sage-python
"""Sweep exact half-lattice depth spectra on the custom rootless MW17 frames.

Corollary S4 of ``SPECIALIZATION_QUOTIENT_AND_RANK_JUMP_THEOREMS.md``
identifies the minimum norm in each class of ``M/2M`` with four times the
squared distance from the corresponding midpoint to ``M``.  This script
applies that finite statement to every rootless rank-17 frame attached to the
legacy custom lattices ``NS0002`` through ``NS0048``.  The published
``NS0001`` lattice is deliberately kept out of the target batch.

The primary closest-vector enumeration uses fplll double-double arithmetic.
Every returned norm is recomputed over the integers.  Every deepest class and
a deterministic stride through all other classes are independently repeated
with 256-bit MPFR Gram--Schmidt arithmetic.  Work is checkpointed per frame.

This is a lattice-geometric scheduling certificate.  It proves no rational
K3 realization, specialized point, Selmer bound, or rank jump.
"""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
from typing import Sequence

from fpylll import Enumeration, FPLLL, GSO, IntegerMatrix
from sage.all import ZZ, matrix
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
CATALOGUE = (
    ROOT / "artifacts/generated-results/elkies-k3-rank7-auxiliary-catalogue-v1.json"
)
CLASSIFIER = (
    ROOT
    / "artifacts/generated-results/elkies-k3-rank19-arithmetic-marking-classifier-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-custom-ns-rootless-half-lattice-sweep-v1.json"
)
DEFAULT_CHECKPOINT = (
    ROOT
    / "artifacts/generated-results/checkpoints"
    / "elkies-k3-custom-ns-rootless-half-lattice-sweep-v1.checkpoint.json"
)

DEGREE = 2
RANK = 17
AUDIT_PRECISION = 256
AUDIT_STRIDE = 4093
EXPECTED_CUSTOM_NS_COUNT = 47
EXPECTED_ELIGIBLE_NS_COUNT = 32
EXPECTED_INELIGIBLE_NS_COUNT = 15
EXPECTED_FRAME_COUNT = 136
CHECKSUM_MODULI = (1_000_000_007, 1_000_000_009)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def gram_digest(gram: Sequence[Sequence[int]]) -> str:
    text = "\n".join(" ".join(map(str, row)) for row in gram) + "\n"
    return hashlib.sha256(text.encode()).hexdigest()


def canonical_hash(value) -> str:
    text = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(text.encode()).hexdigest()


def rational_string(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def exact_norm(vector: Sequence[int], gram: Sequence[Sequence[int]]) -> int:
    dimension = len(vector)
    return sum(
        int(vector[left]) * int(gram[left][right]) * int(vector[right])
        for left in range(dimension)
        for right in range(dimension)
    )


class CosetOracle:
    """Floating CVP enumeration with exact congruence and norm checks."""

    def __init__(
        self,
        gram: Sequence[Sequence[int]],
        *,
        float_type: str,
        precision: int | None = None,
    ) -> None:
        self.gram = tuple(tuple(map(int, row)) for row in gram)
        self.dimension = len(self.gram)
        if precision is not None:
            FPLLL.set_precision(int(precision))
        self.gso = GSO.Mat(
            IntegerMatrix.from_matrix(self.gram),
            gram=True,
            float_type=float_type,
            update=True,
        )
        self.mu = tuple(
            tuple(
                self.gso.get_mu(row, column) if row > column else 0.0
                for column in range(self.dimension)
            )
            for row in range(self.dimension)
        )
        # The 0/1 residue divided by two is a uniform initial target lift.
        self.distance_bound = (
            sum(abs(value) for row in self.gram for value in row) / 4.0 + 1.0
        )

    def solve(self, mask: int) -> tuple[int, tuple[int, ...], float]:
        residue = tuple(
            (int(mask) >> index) & 1 for index in range(self.dimension)
        )
        target = [
            -(
                residue[index]
                + sum(
                    residue[row] * self.mu[row][index]
                    for row in range(index + 1, self.dimension)
                )
            )
            / DEGREE
            for index in range(self.dimension)
        ]
        solutions = Enumeration(self.gso).enumerate(
            0,
            self.dimension,
            self.distance_bound,
            0,
            target=target,
        )
        if not solutions:
            raise ArithmeticError("CVP enumeration returned no solution")
        reported_distance, coordinates = solutions[0]
        closest = tuple(int(round(value)) for value in coordinates)
        if any(
            abs(value - integer) > 1.0e-7
            for value, integer in zip(coordinates, closest)
        ):
            raise ArithmeticError("CVP enumeration returned nonintegral coordinates")
        representative = tuple(
            residue[index] + DEGREE * closest[index]
            for index in range(self.dimension)
        )
        if any(
            (representative[index] - residue[index]) % DEGREE
            for index in range(self.dimension)
        ):
            raise ArithmeticError("CVP representative is in the wrong parity class")
        norm = exact_norm(representative, self.gram)
        error = abs(DEGREE * DEGREE * float(reported_distance) - norm)
        if error > 1.0e-6 or norm < 0 or norm % 2:
            raise ArithmeticError(f"invalid CVP norm={norm}, error={error}")
        return norm, representative, error


def audit_mask(frame_id: str, mask: int) -> bool:
    residue = int(hashlib.sha256(frame_id.encode()).hexdigest()[:8], 16) % AUDIT_STRIDE
    return int(mask) % AUDIT_STRIDE == residue


def sweep_frame(task: dict) -> dict:
    frame_id = task["frame_id"]
    gram = task["gram"]
    if gram_digest(gram) != task["gram_sha256"]:
        raise ArithmeticError(f"{frame_id}: Gram hash changed")

    primary = CosetOracle(gram, float_type="dd")
    histogram: Counter[int] = Counter()
    records = []
    maximum_error = 0.0
    checksums = {modulus: 0 for modulus in CHECKSUM_MODULI}
    stride_masks = []
    for mask in range(1 << RANK):
        norm, representative, error = primary.solve(mask)
        histogram[norm] += 1
        records.append((norm, mask, representative))
        maximum_error = max(maximum_error, error)
        for modulus in CHECKSUM_MODULI:
            checksums[modulus] = (
                checksums[modulus] + (mask + 1) * (norm + 1)
            ) % modulus
        if audit_mask(frame_id, mask):
            stride_masks.append(mask)

    maximum = max(histogram)
    deepest = [row for row in records if row[0] == maximum]
    audit_targets = sorted(set(stride_masks) | {row[1] for row in deepest})
    audit = CosetOracle(
        gram,
        float_type="mpfr",
        precision=AUDIT_PRECISION,
    )
    by_mask = {mask: (norm, representative) for norm, mask, representative in records}
    audit_error = 0.0
    for mask in audit_targets:
        norm, representative, error = audit.solve(mask)
        expected_norm, unused_primary_representative = by_mask[mask]
        if norm != expected_norm or exact_norm(representative, gram) != norm:
            raise ArithmeticError(f"{frame_id}: MPFR audit minimum mismatch at {mask}")
        audit_error = max(audit_error, error)

    deepest.sort(key=lambda row: row[1])
    deepest_records = [
        {
            "mask": mask,
            "representative": list(map(int, representative)),
        }
        for unused_norm, mask, representative in deepest
    ]
    rho_2_squared = Fraction(maximum, 4)
    return {
        "frame_id": frame_id,
        "ns_id": task["ns_id"],
        "surface_id": task["surface_id"],
        "determinant": int(task["determinant"]),
        "gram_sha256": task["gram_sha256"],
        "rank": RANK,
        "translation_cosets": 1 << RANK,
        "minimum_norm_histogram": {
            str(norm): count for norm, count in sorted(histogram.items())
        },
        "maximum_mu_2": maximum,
        "maximum_midpoint_depth_squared": rational_string(rho_2_squared),
        "rho_2_squared": rational_string(rho_2_squared),
        "deepest_class_count": len(deepest),
        "deepest_masks": [mask for unused_norm, mask, unused_rep in deepest],
        "deepest_masks_sha256": canonical_hash(
            [mask for unused_norm, mask, unused_rep in deepest]
        ),
        "deepest_minimum_representatives_sha256": canonical_hash(deepest_records),
        "deepest_minimum_representative_examples": deepest_records[:4],
        "numerical_certificate": {
            "all_returned_norms_recomputed_over_ZZ": True,
            "primary_backend": "fplll GSO dd",
            "maximum_primary_distance_to_exact_norm_error": maximum_error,
            "audit_backend": f"fplll GSO mpfr {AUDIT_PRECISION}-bit",
            "all_deepest_classes_cross_precision_audited": True,
            "deterministic_nonmaximum_audit_stride": AUDIT_STRIDE,
            "cross_precision_masks_audited": len(audit_targets),
            "maximum_audit_distance_to_exact_norm_error": audit_error,
            **{
                f"checksum_mod_{modulus}": checksums[modulus]
                for modulus in CHECKSUM_MODULI
            },
        },
    }


def custom_ns_selection(catalogue: dict) -> tuple[list[dict], list[dict]]:
    eligible = []
    ineligible = []
    custom_surfaces = [
        surface
        for surface in catalogue["surfaces"]
        if surface.get("legacy_ns_ids")
        and surface["legacy_ns_ids"] != ["NS0001"]
    ]
    if len(custom_surfaces) != EXPECTED_CUSTOM_NS_COUNT:
        raise ArithmeticError("custom NS catalogue cardinality changed")
    for surface in custom_surfaces:
        if len(surface["legacy_ns_ids"]) != 1:
            raise ArithmeticError("a custom surface has ambiguous legacy NS labels")
        ns_id = surface["legacy_ns_ids"][0]
        frames = [
            frame
            for frame in surface["frames"]
            if int(frame["mw_rank_for_rho_19"]) == RANK
            and int(frame["root_rank"]) == 0
            and frame["root_type"] == "0"
            and int(frame["signed_root_count"]) == 0
        ]
        if not frames:
            ineligible.append(
                {
                    "ns_id": ns_id,
                    "surface_id": surface["surface_id"],
                    "determinant": int(surface["determinant"]),
                    "maximum_catalogued_mw_rank": max(
                        int(frame["mw_rank_for_rho_19"])
                        for frame in surface["frames"]
                    ),
                    "reason": (
                        "NO_ATTACHED_ROOTLESS_MW17_FRAME; a rooted rank-17 frame "
                        "Gram is not the Mordell--Weil height lattice"
                    ),
                }
            )
            continue
        for frame in frames:
            gram = frame["gram"]
            gram_matrix = matrix(ZZ, gram)
            intrinsics = frame.get("rootless_intrinsics") or {}
            if (
                gram_matrix.dimensions() != (RANK, RANK)
                or not gram_matrix.is_symmetric()
                or not gram_matrix.is_positive_definite()
                or any(int(gram_matrix[index, index]) % 2 for index in range(RANK))
                or int(gram_matrix.det()) != int(surface["determinant"])
                or int(intrinsics.get("minimum_squared_norm", 0)) < 4
                or gram_digest(gram) != frame["gram_sha256"]
            ):
                raise ArithmeticError(f"{frame['frame_id']}: invalid rootless MW17 input")
            eligible.append(
                {
                    "ns_id": ns_id,
                    "surface_id": surface["surface_id"],
                    "determinant": int(surface["determinant"]),
                    "frame_id": frame["frame_id"],
                    "gram": gram,
                    "gram_sha256": frame["gram_sha256"],
                }
            )
    eligible.sort(key=lambda row: (row["ns_id"], row["frame_id"]))
    ineligible.sort(key=lambda row: row["ns_id"])
    if (
        len({row["ns_id"] for row in eligible}) != EXPECTED_ELIGIBLE_NS_COUNT
        or len(ineligible) != EXPECTED_INELIGIBLE_NS_COUNT
        or len(eligible) != EXPECTED_FRAME_COUNT
    ):
        raise ArithmeticError("eligible custom rootless sweep dimensions changed")
    return eligible, ineligible


def arithmetic_rows(classifier: dict) -> dict[str, dict]:
    return {row["surface_id"]: row for row in classifier["candidates"]}


def ns_summaries(frame_rows: list[dict], classifier_by_surface: dict[str, dict]) -> list[dict]:
    grouped: dict[str, list[dict]] = {}
    for row in frame_rows:
        grouped.setdefault(row["ns_id"], []).append(row)
    result = []
    for ns_id, rows in sorted(grouped.items()):
        maxima = Counter(int(row["maximum_mu_2"]) for row in rows)
        best = max(maxima)
        rho_2_squared = Fraction(best, 4)
        surface_id = rows[0]["surface_id"]
        arithmetic = classifier_by_surface[surface_id]
        result.append(
            {
                "ns_id": ns_id,
                "surface_id": surface_id,
                "determinant": rows[0]["determinant"],
                "rootless_mw17_frame_count": len(rows),
                "maximum_mu_2_histogram_across_frames": {
                    str(norm): count for norm, count in sorted(maxima.items())
                },
                "largest_maximum_mu_2": best,
                "largest_rho_2_squared": rational_string(rho_2_squared),
                "frames_attaining_largest_hole": [
                    row["frame_id"] for row in rows if row["maximum_mu_2"] == best
                ],
                "arithmetic_marking": {
                    "classification": arithmetic["classification"],
                    "phase_2_certificate_status": arithmetic.get(
                        "phase_2_certificate_status"
                    ),
                    "different_ns_foundry_equation_eligible": arithmetic[
                        "different_ns_foundry_equation_eligible"
                    ],
                    "next_arithmetic_gate": arithmetic.get("next_arithmetic_gate"),
                },
            }
        )
    return result


def checkpoint_configuration(tasks: list[dict]) -> dict:
    return {
        "catalogue": relative(CATALOGUE),
        "catalogue_sha256": digest(CATALOGUE),
        "classifier": relative(CLASSIFIER),
        "classifier_sha256": digest(CLASSIFIER),
        "frame_ids": [task["frame_id"] for task in tasks],
        "frame_gram_sha256": {
            task["frame_id"]: task["gram_sha256"] for task in tasks
        },
        "degree": DEGREE,
        "rank": RANK,
        "cosets_per_frame": 1 << RANK,
        "audit_precision_bits": AUDIT_PRECISION,
        "audit_stride": AUDIT_STRIDE,
    }


def load_checkpoint(path: Path, configuration: dict) -> dict[str, dict]:
    if not path.is_file():
        return {}
    payload = json.loads(path.read_text())
    if (
        payload.get("schema")
        != "elkies-k3.custom-ns-rootless-half-lattice-sweep-checkpoint.v1"
        or payload.get("configuration") != configuration
    ):
        raise ArithmeticError("checkpoint configuration does not match this sweep")
    return {row["frame_id"]: row for row in payload["completed_frames"]}


def write_checkpoint(path: Path, configuration: dict, completed: dict[str, dict]) -> None:
    payload = {
        "schema": "elkies-k3.custom-ns-rootless-half-lattice-sweep-checkpoint.v1",
        "configuration": configuration,
        "completed_frames": [completed[key] for key in sorted(completed)],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    parser.add_argument("--workers", type=int, default=min(4, os.cpu_count() or 1))
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.workers <= 0:
        parser.error("--workers must be positive")

    catalogue = json.loads(CATALOGUE.read_text())
    classifier = json.loads(CLASSIFIER.read_text())
    if catalogue.get("schema") != "elkies-k3.rank7-auxiliary-catalogue.v1":
        raise ArithmeticError("unexpected auxiliary catalogue schema")
    if classifier.get("schema") != "elkies-k3.rank19-arithmetic-marking-classifier.v1":
        raise ArithmeticError("unexpected arithmetic classifier schema")
    tasks, ineligible = custom_ns_selection(catalogue)
    configuration = checkpoint_configuration(tasks)
    completed = load_checkpoint(args.checkpoint.resolve(), configuration)
    unknown_completed = set(completed) - {task["frame_id"] for task in tasks}
    if unknown_completed:
        raise ArithmeticError(f"checkpoint contains unknown frames: {sorted(unknown_completed)}")
    pending = [task for task in tasks if task["frame_id"] not in completed]

    if pending:
        with ProcessPoolExecutor(
            max_workers=args.workers,
            mp_context=multiprocessing.get_context("fork"),
        ) as executor:
            futures = {executor.submit(sweep_frame, task): task for task in pending}
            for future in as_completed(futures):
                task = futures[future]
                row = future.result()
                completed[row["frame_id"]] = row
                write_checkpoint(args.checkpoint.resolve(), configuration, completed)
                print(
                    "CUSTOMNSHALF|completed={}/{}|ns={}|frame={}|max_mu2={}|deepest={}".format(
                        len(completed),
                        len(tasks),
                        row["ns_id"],
                        row["frame_id"],
                        row["maximum_mu_2"],
                        row["deepest_class_count"],
                    ),
                    flush=True,
                )
    if len(completed) != len(tasks):
        raise ArithmeticError("custom NS sweep ended with incomplete frame coverage")

    frame_rows = [completed[task["frame_id"]] for task in tasks]
    classifier_by_surface = arithmetic_rows(classifier)
    summaries = ns_summaries(frame_rows, classifier_by_surface)
    global_maxima = Counter(row["maximum_mu_2"] for row in frame_rows)
    global_best = max(global_maxima)
    global_rho_2_squared = Fraction(global_best, 4)
    payload = {
        "schema": "elkies-k3.custom-ns-rootless-half-lattice-sweep.v1",
        "status": "PASS_COMPLETE_CUSTOM_NS_ROOTLESS_MW17_HALF_LATTICE_SWEEP",
        "scope": {
            "custom_ns_definition": (
                "legacy custom labels NS0002 through NS0048 in the rank-seven "
                "auxiliary catalogue; published NS0001 is excluded"
            ),
            "custom_ns_count": EXPECTED_CUSTOM_NS_COUNT,
            "eligible_custom_ns_count": len(summaries),
            "ineligible_custom_ns_count": len(ineligible),
            "rootless_mw17_frame_count": len(frame_rows),
            "parity_cosets_per_frame": 1 << RANK,
            "parity_cosets_total": len(frame_rows) * (1 << RANK),
        },
        "theorem_application": {
            "source": "elkies-k3/SPECIALIZATION_QUOTIENT_AND_RANK_JUMP_THEOREMS.md#5-half-lattice-depth-and-covering-radius",
            "corollary": "S4",
            "identity": "D_2(c)=mu_2(c)/4 and rho_2(M)=sqrt(max_c mu_2(c))/2",
            "interpretation": (
                "This is an exact discrete old-lattice midpoint-hole invariant. "
                "It can schedule a later equation-level half-lattice search but is "
                "not a success probability or rank-jump predictor."
            ),
        },
        "accounting": {
            "maximum_mu_2_histogram_across_all_frames": {
                str(norm): count for norm, count in sorted(global_maxima.items())
            },
            "largest_maximum_mu_2": global_best,
            "largest_rho_2_squared": rational_string(global_rho_2_squared),
            "frames_attaining_global_largest_hole": [
                row["frame_id"]
                for row in frame_rows
                if row["maximum_mu_2"] == global_best
            ],
            "ns_attaining_global_largest_hole": sorted(
                {
                    row["ns_id"]
                    for row in frame_rows
                    if row["maximum_mu_2"] == global_best
                }
            ),
            "published_r17_reference_not_swept": {
                "ns_id": "NS0001",
                "maximum_mu_2": 12,
                "deepest_class_count": 43,
                "rho_2_squared": "3",
                "source": "elkies-k3/SPECIALIZATION_QUOTIENT_AND_RANK_JUMP_THEOREMS.md",
            },
        },
        "custom_ns_summaries": summaries,
        "ineligible_custom_ns": ineligible,
        "frame_spectra": frame_rows,
        "certificate_method": {
            "primary": (
                "complete fplll double-double CVP enumeration of all 2^17 parity "
                "classes with exact integral norm and congruence recomputation"
            ),
            "audit": (
                "independent 256-bit MPFR GSO repetition of every deepest class "
                "and a deterministic stride through every frame"
            ),
            "checkpoint_granularity": "one complete frame spectrum",
            "input_gram_validation": (
                "rank 17, symmetric, positive definite, even, determinant equal to "
                "the NS determinant, rootless minimum at least four, and hash match"
            ),
        },
        "proof_boundary": (
            "The sweep computes S4 only for attached rootless MW17 frame lattices. "
            "It does not construct a K3 equation or QQ marking, specialize a family, "
            "find a point, calculate L/M or a residual Selmer group, prove a rank "
            "jump, or turn midpoint depth into prospective enrichment. Arithmetic "
            "marking exclusions continue to override geometric scheduling."
        ),
        "inputs": {
            relative(path): digest(path)
            for path in (Path(__file__).resolve(), CATALOGUE, CLASSIFIER)
        },
        "software": {
            "sage_version": SAGE_VERSION,
            "fpylll_precision_for_audit_bits": AUDIT_PRECISION,
        },
        "reproducing_command": (
            "sage -python elkies-k3/scripts/sweep_custom_ns_half_lattice_depths.sage"
        ),
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output_path = args.output.resolve()
    if args.check:
        if not output_path.is_file() or output_path.read_text() != serialized:
            raise ArithmeticError("stored custom NS half-lattice sweep differs from replay")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    if args.checkpoint.resolve().is_file():
        args.checkpoint.resolve().unlink()
    print(
        "CUSTOMNSHALF|ns={}|frames={}|cosets={}|max_mu2={}|status={}|output={}".format(
            len(summaries),
            len(frame_rows),
            len(frame_rows) * (1 << RANK),
            global_best,
            payload["status"],
            relative(output_path),
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
