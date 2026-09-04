#!/usr/bin/env sage -python
"""Blind half-lattice replay on frozen subgroup-only 273/302/control inputs."""

from __future__ import annotations

import argparse
from decimal import Decimal, getcontext
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import subprocess
import sys
from typing import Sequence

from fpylll import Enumeration, GSO, IntegerMatrix
from sage.all import EllipticCurve, QQ, pari


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC = ROOT / "elliptic-curves"
CAS = ELLIPTIC / "cas"
INPUT = ELLIPTIC / "data/half_lattice_heldout_subgroup_inputs_v1.json"
ENGINE = CAS / "half_lattice_fake_descent_replay.sage"
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "half_lattice_heldout_273_302_blind_v1.json"
)
sys.path[:0] = [str(ELLIPTIC), str(CAS)]

from mod2_reduction_independence import (  # noqa: E402
    combined_mod2_rank,
    find_mod2_reduction_certificate,
)
from search_nagao_u135_alternate_covers import relation_proposals  # noqa: E402


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def binary_rank(values) -> int:
    pivots: dict[int, int] = {}
    for value in values:
        value = int(value)
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    return len(pivots)


def restricted_signature_rank(signatures, columns: Sequence[int]) -> int:
    values = []
    for signature in signatures:
        for row in signature.rows:
            values.append(
                sum((int(row[column]) & 1) << offset for offset, column in enumerate(columns))
            )
    return binary_rank(values)


class CosetOracle:
    def __init__(self, gram: Sequence[Sequence[int]]) -> None:
        self.gram = tuple(tuple(int(value) for value in row) for row in gram)
        self.dimension = len(self.gram)
        self.gso = GSO.Mat(
            IntegerMatrix.from_matrix(self.gram),
            gram=True,
            float_type="dd",
            update=True,
        )
        self.mu = tuple(
            tuple(self.gso.get_mu(i, j) if i > j else 0.0 for j in range(self.dimension))
            for i in range(self.dimension)
        )
        self.distance_bound = sum(abs(value) for row in self.gram for value in row) / 4 + 1.0

    def solve(self, mask: int):
        residue = tuple((mask >> index) & 1 for index in range(self.dimension))
        target = [
            -(residue[i] + sum(residue[j] * self.mu[j][i] for j in range(i + 1, self.dimension))) / 2
            for i in range(self.dimension)
        ]
        solutions = Enumeration(self.gso).enumerate(
            0, self.dimension, self.distance_bound, 0, target=target
        )
        if not solutions:
            raise ArithmeticError("rounded-form CVP returned no solution")
        reported, coordinates = solutions[0]
        closest = tuple(int(round(value)) for value in coordinates)
        representative = tuple(residue[index] + 2 * closest[index] for index in range(self.dimension))
        norm = sum(
            representative[i] * self.gram[i][j] * representative[j]
            for i in range(self.dimension)
            for j in range(self.dimension)
        )
        if abs(4 * float(reported) - norm) > 1.0e-6:
            raise ArithmeticError("rounded-form CVP exact norm check failed")
        return norm, representative


def height_gram(model, points):
    getcontext().prec = 110
    pari.default("realprecision", 110)
    raw = pari(EllipticCurve(QQ, list(model))).ellheightmatrix([list(point) for point in points])
    n = len(points)
    gram = tuple(
        tuple(Decimal(str(raw[i, j])) for j in range(n)) for i in range(n)
    )
    if max(abs(gram[i][j] - gram[j][i]) for i in range(n) for j in range(n)) > Decimal("1e-90"):
        raise ArithmeticError("canonical-height Gram is unexpectedly asymmetric")
    return gram


def quadratic(gram, vector) -> Decimal:
    return sum(
        Decimal(vector[i]) * gram[i][j] * Decimal(vector[j])
        for i in range(len(vector))
        for j in range(len(vector))
    )


def point_record(point):
    return {"x": str(point[0]), "y": str(point[1])}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--height-bound", type=int, default=100_000)
    parser.add_argument("--timeout-seconds", type=float, default=15.0)
    parser.add_argument("--stack-bytes", type=int, default=1_000_000_000)
    parser.add_argument("--only", help="optional curve/config substring for diagnostics")
    args = parser.parse_args()

    frozen = json.loads(args.input.read_text())
    if frozen.get("status") != "FROZEN_SUBGROUP_ONLY_INPUTS_NO_HELDOUT_COORDINATES":
        raise ValueError("subgroup-only input is not frozen")
    if frozen["boundary"]["output_contains_heldout_coordinates"] is not False:
        raise ValueError("subgroup input does not enforce the held-out boundary")
    engine = SourceFileLoader("half_lattice_engine_for_heldout", str(ENGINE)).load_module()

    results = []
    for curve in frozen["curves"]:
        model = tuple(Fraction(value) for value in curve["short_model"])
        for configuration in curve["configurations"]:
            key = f"{curve['label']}/{configuration['id']}"
            if args.only and args.only not in key:
                continue
            basis = tuple(
                (Fraction(point[0]), Fraction(point[1]))
                for point in configuration["starting_subgroup_points"]
            )
            dimension = len(basis)
            baseline_signatures = find_mod2_reduction_certificate(
                model, basis, prime_bound=1000
            )
            if combined_mod2_rank(baseline_signatures, dimension) != dimension:
                raise ArithmeticError(f"{key}: starting subgroup independence did not replay")
            gram = height_gram(model, basis)
            rounded = tuple(
                tuple(int((value * Decimal(1_000_000)).to_integral_value()) for value in row)
                for row in gram
            )
            oracle = CosetOracle(rounded)
            ranked = []
            for mask in range(1 << dimension):
                unused_norm, representative = oracle.solve(mask)
                ranked.append((quadratic(gram, representative) / 4, mask, representative))
            ranked.sort(key=lambda item: (-item[0], item[1]))
            selected = ranked[:43]

            rounded_check = tuple(
                tuple(int((value * Decimal(100_000)).to_integral_value()) for value in row)
                for row in gram
            )
            check_oracle = CosetOracle(rounded_check)
            representative_disagreements = sum(
                check_oracle.solve(mask)[1] != representative
                for unused_depth, mask, representative in selected
            )

            discoveries: dict[tuple[Fraction, Fraction], set[int]] = {}
            cover_rows = []
            for position, (depth, mask, representative) in enumerate(selected, 1):
                outcome = engine.run_quartic_search(
                    mask=mask,
                    representative=representative,
                    short_model=model,
                    generic_points=basis,
                    height_bound=args.height_bound,
                    timeout_seconds=args.timeout_seconds,
                    stack_bytes=args.stack_bytes,
                )
                for point in outcome.curve_points:
                    discoveries.setdefault(point, set()).add(mask)
                record = outcome.record
                cover_rows.append(
                    {
                        "position": position,
                        "mask": mask,
                        "hex": f"0x{mask:0{(dimension + 3) // 4}x}",
                        "depth": str(depth),
                        "representative": list(representative),
                        "status": record["status"],
                        "finite_point_count": len(outcome.curve_points),
                        "raw_coefficient_bits": record["raw_rational_coefficient_maximum_bits"],
                        "integral_coefficient_bits": record["integral_model_maximum_coefficient_bits"],
                        "reduced_coefficient_bits": (
                            record["reduced_model"]["maximum_coefficient_bits"]
                            if record["status"] == "bounded_search_complete"
                            else None
                        ),
                        "modular_density_product": (
                            record["local_stage"]["joint_independent_density_product"]
                            if record["status"] == "bounded_search_complete"
                            else None
                        ),
                        "search_milliseconds": record.get("search_milliseconds"),
                        "wall_seconds": record["wall_seconds"],
                    }
                )
                print(
                    f"HELDOUTHALF|case={key}|cover={position}/43|mask={mask:#x}|"
                    f"status={record['status']}|points={len(outcome.curve_points)}",
                    flush=True,
                )

            basis_with_signs = {
                signed for point in basis for signed in (point, (point[0], -point[1]))
            }
            candidates = tuple(
                sorted(
                    (point for point in discoveries if point not in basis_with_signs),
                    key=lambda point: (
                        max(abs(point[0].numerator).bit_length(), point[0].denominator.bit_length()),
                        max(abs(point[1].numerator).bit_length(), point[1].denominator.bit_length()),
                        point,
                    ),
                )
            )
            proposals = (
                relation_proposals(
                    model, basis, candidates, timeout=60.0, stack_bytes=args.stack_bytes
                )
                if candidates
                else ()
            )
            unexplained = tuple(
                point for point, (unused_relation, exact) in zip(candidates, proposals) if not exact
            )
            combined_signatures = (
                find_mod2_reduction_certificate(
                    model, basis + unexplained, prime_bound=1500
                )
                if unexplained
                else baseline_signatures
            )
            combined_rank = combined_mod2_rank(
                combined_signatures, dimension + len(unexplained)
            )
            columns = list(range(dimension))
            current_rank = restricted_signature_rank(combined_signatures, columns)
            selected_offsets = []
            for offset in range(len(unexplained)):
                trial = columns + [dimension + offset]
                trial_rank = restricted_signature_rank(combined_signatures, trial)
                if trial_rank > current_rank:
                    selected_offsets.append(offset)
                    columns = trial
                    current_rank = trial_rank
            if current_rank != combined_rank:
                raise ArithmeticError(f"{key}: greedy quotient basis lost rank")
            unexplained_index = {point: index for index, point in enumerate(unexplained)}
            candidate_rows = []
            for point, (relation, exact) in zip(candidates, proposals):
                if exact:
                    continue
                offset = unexplained_index[point]
                candidate_rows.append(
                    {
                        "point": point_record(point),
                        "source_masks": sorted(discoveries[point]),
                        "selected_for_independent_quotient_basis": offset in selected_offsets,
                    }
                )
            results.append(
                {
                    "curve": curve["label"],
                    "configuration": configuration["id"],
                    "dimension": dimension,
                    "included_public_indices_one_based": configuration[
                        "included_public_indices_one_based"
                    ],
                    "selected_top_class_count": len(selected),
                    "top43_representative_disagreements_scale_1e5_vs_1e6": representative_disagreements,
                    "cover_records": cover_rows,
                    "blind_result": {
                        "distinct_nonbasis_candidates": len(candidates),
                        "unexplained_candidate_count": len(unexplained),
                        "finite_mod2_rank_lower_bound": combined_rank,
                        "finite_mod2_quotient_gain": combined_rank - dimension,
                        "selected_independent_candidate_count": len(selected_offsets),
                        "candidate_points": candidate_rows,
                    },
                }
            )

    payload = {
        "schema": "elliptic-curves.half-lattice-heldout-273-302-blind.v1",
        "status": "PASS_BOUNDED_HELDOUT_SUBGROUP_SEARCH",
        "blindness_boundary": {
            "search_input": str(args.input.relative_to(ROOT)),
            "search_imported_public_curve_modules": False,
            "search_loaded_heldout_point_coordinates": False,
            "verification_fixture": "separate verifier not yet run",
        },
        "input_hashes": {
            str(args.input.relative_to(ROOT)): digest(args.input),
            str(ENGINE.relative_to(ROOT)): digest(ENGINE),
            str(Path(__file__).resolve().relative_to(ROOT)): digest(Path(__file__).resolve()),
        },
        "declared_budget": {
            "top_classes_per_configuration": 43,
            "height_bound_each": args.height_bound,
            "timeout_seconds_each": args.timeout_seconds,
            "stack_bytes_each": args.stack_bytes,
            "single_pass_no_retry": True,
        },
        "configuration_count": len(results),
        "results": results,
        "claim_boundary": [
            "Starting subgroup and discovered-point independence are exact finite-reduction results.",
            "The rounded canonical-height CVPs are reproducibly checked at the displayed scales but are numerical evidence.",
            "All point-search misses are bounded, not nonexistence claims.",
            "No family or K3 provenance is inferred from recovery performance.",
        ],
        "software": {
            "sage_python": sys.executable,
            "pari_gp": subprocess.run(
                ["gp", "-fq"], input="print(version());quit\n", text=True, capture_output=True
            ).stdout.strip(),
        },
        "reproducing_command": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elliptic-curves/cas/replay_half_lattice_heldout_subgroups.sage"
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        "HELDOUTHALF|status=PASS|gains="
        + ",".join(str(row["blind_result"]["finite_mod2_quotient_gain"]) for row in results)
        + f"|output={display_path(args.output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
