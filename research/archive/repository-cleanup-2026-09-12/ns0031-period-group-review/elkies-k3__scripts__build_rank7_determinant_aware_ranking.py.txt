#!/usr/bin/env python3
"""Build the fail-closed determinant-aware rootless-MW17 factory ranking.

This is a cheap pre-solver pass over the exact rank-seven auxiliary catalogue.
It rejects surfaces without an actual rootless rank-17 frame, applies the
Blichfeldt necessary determinant bound, audits discriminant-form length, and
checks an explicit even rank-17 Gram witness before using any discovery score.

The surviving surfaces are placed in determinant regimes rather than ordered
by raw determinant.  Arithmetic field evidence, source-equation precursors,
and compiler corridors are typed readiness gates.  Rank-jump mechanism
evidence is reported separately and raw multisection counts are never used in
the priority key.  A proved obstruction to the full rational rank-19 marking
removes a surface from the arithmetic candidate queue while preserving it in a
separate exact-rejection ledger.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CATALOGUE = (
    ROOT / "artifacts/generated-results/elkies-k3-rank7-auxiliary-catalogue-v1.json"
)
SURFACE_LEDGER = (
    ROOT / "artifacts/generated-results/elkies-k3-rank7-surface-pareto-v1.json"
)
EVIDENCE = (
    ROOT
    / "elkies-k3/data/lattice-foundry/determinant-aware-ranking-evidence-v1.json"
)
ARITHMETIC_CLASSIFIER = (
    ROOT
    / "artifacts/generated-results/elkies-k3-rank19-arithmetic-marking-classifier-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-rank7-determinant-aware-ranking-v1.json"
)

RANK = 17
ROOTLESS_MINIMUM = 4
K3_DISCRIMINANT_LENGTH_MAXIMUM = 3
UNKNOWN_TIER = 3


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def determinant_bareiss(values: list[list[int]]) -> int:
    """Return the exact determinant using fraction-free elimination."""

    matrix = [list(map(int, row)) for row in values]
    size = len(matrix)
    if any(len(row) != size for row in matrix):
        raise ValueError("Gram matrix is not square")
    sign = 1
    previous = 1
    for pivot_index in range(size - 1):
        if matrix[pivot_index][pivot_index] == 0:
            swap = next(
                (
                    index
                    for index in range(pivot_index + 1, size)
                    if matrix[index][pivot_index] != 0
                ),
                None,
            )
            if swap is None:
                return 0
            matrix[pivot_index], matrix[swap] = matrix[swap], matrix[pivot_index]
            sign *= -1
        pivot = matrix[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = (
                    matrix[row][column] * pivot
                    - matrix[row][pivot_index] * matrix[pivot_index][column]
                )
                if numerator % previous:
                    raise ArithmeticError("Bareiss division was not exact")
                matrix[row][column] = numerator // previous
        previous = pivot
        for row in range(pivot_index + 1, size):
            matrix[row][pivot_index] = 0
    return sign * matrix[-1][-1]


def blichfeldt_data() -> dict:
    gamma_upper = (2 / math.pi) * math.gamma(2 + RANK / 2) ** (2 / RANK)
    determinant_lower = (ROOTLESS_MINIMUM / gamma_upper) ** RANK
    return {
        "rank": RANK,
        "required_minimum_squared_norm": ROOTLESS_MINIMUM,
        "hermite_constant_upper_bound": gamma_upper,
        "gram_determinant_lower_bound_real": determinant_lower,
        "gram_determinant_lower_bound_integer": math.ceil(determinant_lower),
        "formula": (
            "gamma_n <= (2/pi)*Gamma(2+n/2)^(2/n), hence "
            "det(Gram) >= (minimum/gamma_n)^n"
        ),
        "status_boundary": (
            "Necessary only. Passing the bound does not construct a lattice or "
            "prove that a required discriminant form is realized."
        ),
    }


def discriminant_length(surface: dict) -> int:
    invariants = surface["surface_key"]["ns_discriminant_form_key"]["invariants"]
    if math.prod(map(int, invariants)) != int(surface["determinant"]):
        raise ValueError(f"discriminant invariants changed for {surface['surface_id']}")
    return sum(int(value) > 1 for value in invariants)


def even_rootless_witnesses(surface: dict) -> list[dict]:
    witnesses = []
    determinant = int(surface["determinant"])
    for frame in surface["frames"]:
        if int(frame["mw_rank_for_rho_19"]) != RANK:
            continue
        gram = frame["gram"]
        exact_determinant = determinant_bareiss(gram)
        symmetric = all(
            gram[row][column] == gram[column][row]
            for row in range(RANK)
            for column in range(RANK)
        )
        even = all(int(gram[index][index]) % 2 == 0 for index in range(RANK))
        intrinsics = frame.get("rootless_intrinsics") or {}
        rootless = (
            int(frame["root_rank"]) == 0
            and frame["root_type"] == "0"
            and int(frame["signed_root_count"]) == 0
            and int(intrinsics.get("minimum_squared_norm", 0)) >= ROOTLESS_MINIMUM
        )
        if not (
            len(gram) == RANK
            and symmetric
            and even
            and rootless
            and exact_determinant == determinant
        ):
            continue
        witnesses.append(
            {
                "frame_id": frame["frame_id"],
                "rank": RANK,
                "even": True,
                "rootless": True,
                "minimum_squared_norm": int(intrinsics["minimum_squared_norm"]),
                "determinant": exact_determinant,
                "norm_four_unoriented_pairs": intrinsics.get(
                    "norm_four_unoriented_pairs"
                ),
                "gram_sha256": frame["gram_sha256"],
            }
        )
    return witnesses


def default_coordinate(name: str) -> dict:
    return {
        "tier": UNKNOWN_TIER,
        "status": f"UNKNOWN_NO_{name.upper()}_EVIDENCE_ATTACHED",
        "evidence": [],
    }


def arithmetic_coordinate(
    surface_row: dict, overlay: dict | None, classifier_row: dict
) -> dict:
    if overlay is not None:
        return dict(overlay)
    if (
        classifier_row.get("phase_2_certificate_status")
        == "UNRESOLVED_FOR_EXPLICIT_REASON"
    ):
        return {
            "tier": 1,
            "status": "PASS_EXACT_MARKED_CURVE_NONCM_RATIONAL_LOCUS_UNRESOLVED",
            "evidence": [
                replay["path"]
                for replay in classifier_row.get("certificate_replay", [])
            ],
            "boundary": classifier_row["next_arithmetic_gate"],
            "phase_2_certificate_status": classifier_row[
                "phase_2_certificate_status"
            ],
            "candidate_eligible": True,
        }
    gate = surface_row["moduli"]["t_arithmetic"]["pre_solver_gate"]
    if gate["equation_solver_may_launch"]:
        return {
            "tier": 1,
            "status": "PASS_EXACT_T_ARITHMETIC_CURVE_IDENTIFICATION_ONLY",
            "evidence": [],
            "boundary": (
                "The moduli curve is identified, but no rational elliptic source "
                "marking is inferred."
            ),
        }
    return {
        "tier": 2,
        "status": gate["status"],
        "evidence": [],
        "boundary": "The T-arithmetic attempt is exact but field identification is open.",
    }


def determinant_regime(
    determinant: int,
    theoretical_minimum: int,
    first_imported_rootless: int,
    reference_determinant: int,
    short_vector_benchmark_pass: bool,
) -> dict:
    if determinant < theoretical_minimum:
        return {"tier": 4, "label": "INFEASIBLE_BELOW_BLICHFELDT_BOUND"}
    if determinant < first_imported_rootless:
        return {
            "tier": 0,
            "label": "LOW_DETERMINANT_DESIGN_FRONTIER_IF_EXPLICITLY_REALIZED",
        }
    if determinant <= reference_determinant:
        return {"tier": 1, "label": "OBSERVED_DENSE_ROOTLESS_SWEET_SPOT"}
    if short_vector_benchmark_pass:
        return {"tier": 2, "label": "ABOVE_BAND_SHORT_VECTOR_QUALITY_EXCEPTION"}
    return {"tier": 3, "label": "ABOVE_BAND_SPARSE_PRESSURE"}


def build(
    catalogue: dict,
    surface_ledger: dict,
    evidence: dict,
    arithmetic_classifier: dict,
) -> dict:
    assert catalogue["schema"] == "elkies-k3.rank7-auxiliary-catalogue.v1"
    assert surface_ledger["schema"] == "elkies-k3.rank7-surface-pareto.v1"
    assert evidence["schema"] == "elkies-k3.determinant-aware-ranking-evidence.v1"
    assert arithmetic_classifier["schema"] == (
        "elkies-k3.rank19-arithmetic-marking-classifier.v1"
    )
    surface_rows = {row["surface_id"]: row for row in surface_ledger["surfaces"]}
    assert set(surface_rows) == {row["surface_id"] for row in catalogue["surfaces"]}
    overlays = {row["surface_id"]: row for row in evidence["surfaces"]}
    if not set(overlays) <= set(surface_rows):
        raise ValueError("evidence names a surface absent from the catalogue")

    for row in overlays.values():
        for coordinate in (
            "arithmetic_field_of_definition",
            "source_equation_precursor",
            "compiler_corridor",
            "rank_jump_mechanism",
        ):
            if not 0 <= int(row[coordinate]["tier"]) <= UNKNOWN_TIER:
                raise ValueError(f"bad evidence tier for {row['surface_id']}:{coordinate}")
            for evidence_path in row[coordinate].get("evidence", []):
                if not (ROOT / evidence_path).is_file():
                    raise FileNotFoundError(evidence_path)
        candidate_eligible = row["arithmetic_field_of_definition"].get(
            "candidate_eligible", True
        )
        if not isinstance(candidate_eligible, bool):
            raise ValueError(
                f"bad arithmetic candidate flag for {row['surface_id']}"
            )

    blichfeldt = blichfeldt_data()
    theoretical_minimum = blichfeldt["gram_determinant_lower_bound_integer"]
    candidate_witnesses = {
        surface["surface_id"]: even_rootless_witnesses(surface)
        for surface in catalogue["surfaces"]
    }
    arithmetic_rows = {
        row["surface_id"]: row for row in arithmetic_classifier["candidates"]
    }
    expected_arithmetic_ids = {
        surface_id for surface_id, witnesses in candidate_witnesses.items() if witnesses
    }
    if set(arithmetic_rows) != expected_arithmetic_ids:
        raise ValueError("arithmetic classifier does not cover the exact candidate set")
    realized_determinants = [
        int(surface["determinant"])
        for surface in catalogue["surfaces"]
        if candidate_witnesses[surface["surface_id"]]
    ]
    first_imported_rootless = min(realized_determinants)
    reference_determinant = int(
        evidence["policy"]["published_r17_reference_determinant"]
    )
    reference_pairs = int(
        evidence["policy"]["published_r17_reference_norm_four_unoriented_pairs"]
    )

    rejected = []
    arithmetic_rejected = []
    candidates = []
    for surface in catalogue["surfaces"]:
        surface_id = surface["surface_id"]
        determinant = int(surface["determinant"])
        length = discriminant_length(surface)
        witnesses = candidate_witnesses[surface_id]
        has_rank17 = any(
            int(frame["mw_rank_for_rho_19"]) == RANK for frame in surface["frames"]
        )
        filters = {
            "rootless_mw17_frame": {
                "pass": bool(witnesses),
                "status": (
                    "PASS_EXACT_ROOTLESS_MW17_FRAME_WITNESS"
                    if witnesses
                    else "REJECT_NO_EXPLICIT_ROOTLESS_MW17_FRAME"
                ),
            },
            "hermite_blichfeldt_feasibility": {
                "pass": determinant >= theoretical_minimum,
                "status": (
                    "PASS_NECESSARY_BLICHFELDT_BOUND"
                    if determinant >= theoretical_minimum
                    else "REJECT_BELOW_NECESSARY_BLICHFELDT_BOUND"
                ),
                "determinant": determinant,
                "minimum_integer_determinant": theoretical_minimum,
            },
            "discriminant_form_length": {
                "pass": length <= K3_DISCRIMINANT_LENGTH_MAXIMUM,
                "status": (
                    "PASS_K3_RANK_THREE_COMPLEMENT_LENGTH_BOUND"
                    if length <= K3_DISCRIMINANT_LENGTH_MAXIMUM
                    else "REJECT_DISCRIMINANT_LENGTH_EXCEEDS_THREE"
                ),
                "length": length,
                "maximum": K3_DISCRIMINANT_LENGTH_MAXIMUM,
            },
            "even_rank17_required_discriminant_form": {
                "pass": bool(witnesses),
                "status": (
                    "PASS_EXACT_EVEN_ROOTLESS_RANK17_CATALOGUE_WITNESS"
                    if witnesses
                    else "REJECT_NO_EVEN_ROOTLESS_RANK17_CATALOGUE_WITNESS"
                ),
                "witnesses": witnesses,
            },
        }
        if not all(row["pass"] for row in filters.values()):
            rejected.append(
                {
                    "surface_id": surface_id,
                    "determinant": determinant,
                    "maximum_catalogued_mw_rank": max(
                        int(frame["mw_rank_for_rho_19"])
                        for frame in surface["frames"]
                    ),
                    "has_rank17_frame_before_full_witness_audit": has_rank17,
                    "failed_filters": [
                        name for name, row in filters.items() if not row["pass"]
                    ],
                }
            )
            continue

        overlay = overlays.get(surface_id, {})
        classifier_row = arithmetic_rows[surface_id]
        arithmetic = arithmetic_coordinate(
            surface_rows[surface_id],
            overlay.get("arithmetic_field_of_definition"),
            classifier_row,
        )
        marking_classification = classifier_row["classification"]
        overlay_eligible = arithmetic.get("candidate_eligible", True)
        if (marking_classification == "ARITHMETICALLY_EXCLUDED") != (
            overlay_eligible is False
        ):
            raise ValueError(
                f"arithmetic evidence/classifier disagreement for {surface_id}"
            )
        if marking_classification == "ARITHMETICALLY_EXCLUDED":
            arithmetic_rejected.append(
                {
                    "surface_id": surface_id,
                    "legacy_ns_ids": surface["legacy_ns_ids"],
                    "determinant": determinant,
                    "status": marking_classification,
                    "evidence": [
                        "artifacts/generated-results/elkies-k3-rank19-arithmetic-marking-classifier-v1.json",
                        *arithmetic.get("evidence", []),
                    ],
                    "reason": (
                        "A proved obstruction to a full QQ-rational rank-19 "
                        "Neron--Severi marking excludes the required saturated "
                        "arithmetic MW17 endpoint."
                    ),
                }
            )
            continue
        precursor = dict(
            overlay.get(
                "source_equation_precursor",
                default_coordinate("source_equation_precursor"),
            )
        )
        corridor = dict(
            overlay.get("compiler_corridor", default_coordinate("compiler_corridor"))
        )
        rank_jump = dict(
            overlay.get(
                "rank_jump_mechanism", default_coordinate("rank_jump_mechanism")
            )
        )
        rank_jump.setdefault("raw_multisection_counts_used", False)
        if rank_jump["raw_multisection_counts_used"]:
            raise ValueError("raw multisection counts may not enter rank-jump scoring")

        best_pairs = max(
            int(row["norm_four_unoriented_pairs"])
            for row in witnesses
            if row["norm_four_unoriented_pairs"] is not None
        )
        benchmark_pass = best_pairs >= reference_pairs
        regime = determinant_regime(
            determinant,
            theoretical_minimum,
            first_imported_rootless,
            reference_determinant,
            benchmark_pass,
        )
        readiness = {
            "arithmetic_field_of_definition": arithmetic,
            "source_equation_precursor": precursor,
            "compiler_corridor": corridor,
        }
        strict_passes = {
            "arithmetic_field_of_definition": (
                marking_classification == "ARITHMETICALLY_POSSIBLE"
                and int(arithmetic["tier"]) == 0
            ),
            "source_equation_precursor": int(precursor["tier"]) <= 1,
            "compiler_corridor": int(corridor["tier"]) <= 1,
        }
        missing_count = sum(not value for value in strict_passes.values())
        source = surface_rows[surface_id]["easiest_known_source"]
        priority_key = (
            int(regime["tier"]),
            missing_count,
            sum(int(row["tier"]) for row in readiness.values()),
            int(arithmetic["tier"]),
            int(precursor["tier"]),
            int(corridor["tier"]),
            0 if benchmark_pass else 1,
            -best_pairs,
            int(source["mw_rank"]),
            int(source["reducible_fibre_support_count"]),
            surface_id,
        )
        candidate_record = {
                "surface_id": surface_id,
                "legacy_ns_ids": surface["legacy_ns_ids"],
                "determinant": determinant,
                "hard_theoretical_filters": filters,
                "determinant_regime": regime,
                "determinant_coordinates": {
                    "blichfeldt_lower_boundary_integer": theoretical_minimum,
                    "ratio_to_blichfeldt_boundary": determinant / theoretical_minimum,
                    "first_imported_rootless_determinant": first_imported_rootless,
                    "published_r17_reference_determinant": reference_determinant,
                    "raw_determinant_is_not_a_minimization_objective": True,
                },
                "rootless_short_vector_quality": {
                    "best_norm_four_unoriented_pairs": best_pairs,
                    "published_r17_benchmark": reference_pairs,
                    "passes_published_r17_short_vector_benchmark": benchmark_pass,
                    "boundary": (
                        "This is an exact lattice-density coordinate, not a rank-jump "
                        "mechanism score."
                    ),
                },
                "readiness_filters": readiness,
                "arithmetic_marking_classification": marking_classification,
                "readiness_strict_passes": strict_passes,
                "readiness_missing_gate_count": missing_count,
                "expensive_equation_scoring_eligible": missing_count == 0,
                "rank_jump_mechanism_coordinate": rank_jump,
                "easiest_known_source": {
                    "status": source["status"],
                    "root_type": source["root_type"],
                    "mw_rank": int(source["mw_rank"]),
                    "reducible_fibre_support_count": int(
                        source["reducible_fibre_support_count"]
                    ),
                },
                "priority_key": list(priority_key[:-1]),
                "_priority_key": priority_key,
            }
        if classifier_row.get("phase_2_certificate_status") is not None:
            candidate_record.update(
                {
                    "phase_2_certificate_status": classifier_row[
                        "phase_2_certificate_status"
                    ],
                    "next_arithmetic_gate": classifier_row["next_arithmetic_gate"],
                }
            )
        candidates.append(candidate_record)

    candidates.sort(key=lambda row: row["_priority_key"])
    for rank, row in enumerate(candidates, 1):
        row["determinant_aware_rank"] = rank
        del row["_priority_key"]

    mechanism_order = sorted(
        candidates,
        key=lambda row: (
            int(row["rank_jump_mechanism_coordinate"]["tier"]),
            row["determinant_aware_rank"],
            row["surface_id"],
        ),
    )
    expensive_queue = [
        row["surface_id"]
        for row in candidates
        if row["expensive_equation_scoring_eligible"]
    ]
    next_gate_queue = [
        {
            "surface_id": row["surface_id"],
            "determinant": row["determinant"],
            "determinant_regime": row["determinant_regime"]["label"],
            "missing_gates": [
                name
                for name, passed in row["readiness_strict_passes"].items()
                if not passed
            ],
            **(
                {"next_arithmetic_gate": row["next_arithmetic_gate"]}
                if row.get("phase_2_certificate_status") is not None
                else {}
            ),
        }
        for row in candidates
        if not row["expensive_equation_scoring_eligible"]
    ]
    regime_counts = Counter(row["determinant_regime"]["label"] for row in candidates)

    return {
        "schema": "elkies-k3.rank7-determinant-aware-ranking.v1",
        "status": "PASS_FAIL_CLOSED_THEORY_FILTERS_AND_TYPED_READINESS_RANKING",
        "proof_boundary": {
            "proved": (
                "Every retained row has an explicit even rootless rank-17 Gram of the "
                "required catalogue discriminant form, determinant at least the "
                "Blichfeldt necessary boundary, and discriminant length at most three."
            ),
            "not_proved": (
                "The Blichfeldt bound is not sufficient for lattice existence. The "
                "imported catalogue is not complete in any determinant band. A formal "
                "or finite p-adic precursor is not a QQ model. Short-vector density and "
                "priority rank do not predict arithmetic rank jumps. Only surfaces with "
                "a separately proved rational-marking obstruction are removed by the "
                "arithmetic rejection gate, and UNKNOWN rows are never equation-scoring "
                "eligible."
            ),
        },
        "theoretical_boundary": blichfeldt,
        "determinant_policy": {
            "first_imported_rootless_determinant": first_imported_rootless,
            "published_r17_reference_determinant": reference_determinant,
            "observed_sweet_spot_closed_interval": [
                first_imported_rootless,
                reference_determinant,
            ],
            "unpopulated_design_frontier_closed_interval": [
                theoretical_minimum,
                first_imported_rootless - 1,
            ],
            "boundary": (
                "The interval below the first imported rootless determinant is a design "
                "frontier, not a nonexistence range. Determinant is not minimized inside "
                "a regime; readiness and exact short-vector quality decide priority."
            ),
        },
        "objective_policy": {
            "primary_priority_key": [
                "determinant-regime tier",
                "number of missing strict readiness gates",
                "sum of readiness evidence tiers",
                "arithmetic field-of-definition tier",
                "source-equation precursor tier",
                "compiler-corridor tier",
                "published-R17 short-vector benchmark pass",
                "negative exact norm-four pair count",
                "easiest known source MW rank",
                "source support count",
            ],
            "excluded": [
                "raw determinant minimization",
                "raw bisection counts",
                "raw trisection counts",
                "raw quadrisection counts",
                "rank-jump mechanism coordinate (reported and ranked separately)",
            ],
        },
        "accounting": {
            "catalogue_surfaces": len(catalogue["surfaces"]),
            "theory_rejected_before_scoring": len(rejected),
            "arithmetic_marking_rejected_before_scoring": len(arithmetic_rejected),
            "theory_feasible_with_explicit_rootless_mw17_witness": (
                len(candidates) + len(arithmetic_rejected)
            ),
            "arithmetic_candidates_after_marking_rejections": len(candidates),
            "expensive_equation_scoring_eligible": len(expensive_queue),
            "determinant_regimes": dict(sorted(regime_counts.items())),
        },
        "expensive_equation_scoring_queue": expensive_queue,
        "next_readiness_gate_queue": next_gate_queue,
        "rank_jump_mechanism_priority_order": [
            row["surface_id"] for row in mechanism_order
        ],
        "theory_rejections": rejected,
        "arithmetic_marking_rejections": arithmetic_rejected,
        "candidates": candidates,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalogue", type=Path, default=CATALOGUE)
    parser.add_argument("--surface-ledger", type=Path, default=SURFACE_LEDGER)
    parser.add_argument("--evidence", type=Path, default=EVIDENCE)
    parser.add_argument(
        "--arithmetic-classifier", type=Path, default=ARITHMETIC_CLASSIFIER
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    paths = {
        "catalogue": arguments.catalogue.resolve(),
        "surface_ledger": arguments.surface_ledger.resolve(),
        "evidence": arguments.evidence.resolve(),
        "arithmetic_classifier": arguments.arithmetic_classifier.resolve(),
    }
    payloads = {name: json.loads(path.read_text()) for name, path in paths.items()}
    catalogue_hash = digest(paths["catalogue"])
    if payloads["surface_ledger"]["inputs"].get(relative(paths["catalogue"])) != catalogue_hash:
        raise SystemExit("surface ledger does not match the catalogue hash")
    result = build(
        payloads["catalogue"],
        payloads["surface_ledger"],
        payloads["evidence"],
        payloads["arithmetic_classifier"],
    )
    result["inputs"] = {relative(path): digest(path) for path in paths.values()}
    result["reproduce"] = (
        "python3 elkies-k3/scripts/build_rank7_determinant_aware_ranking.py"
    )
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit("determinant-aware ranking artifact is stale")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    accounting = result["accounting"]
    leader = result["candidates"][0]
    print(
        "DETERRANK|catalogue={}|theory_pass={}|expensive={}|leader={}|det={}|status=PASS".format(
            accounting["catalogue_surfaces"],
            accounting["theory_feasible_with_explicit_rootless_mw17_witness"],
            accounting["expensive_equation_scoring_eligible"],
            leader["surface_id"],
            leader["determinant"],
        )
    )


if __name__ == "__main__":
    main()
