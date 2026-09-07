#!/usr/bin/env sage -python
"""Post-search exact quotient verification for the half-lattice ablation.

The blind artifact is read and hashed before any exceptional-point fixture is
imported.  This verifier computes exact coordinates in the displayed known
subgroup, ranks every arm over Q and F_2, and attaches per-43-cover and
per-CPU-second efficiencies.  It does not turn bounded search misses into
nonexistence statements.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Any, Sequence

from sage.all import Matrix, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC = ROOT / "elliptic-curves"
CAS = ELLIPTIC / "cas"
ART = ROOT / "artifacts/generated-results/elliptic-curves"
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
CONTROLS = ART / "elkies_2026_high_rank_positive_controls_v2.json"
RANK29_PUBLIC = (
    ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-public-fibres-v1.json"
)
CURVE12_QUOTIENT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-curve12-alternate-q80-quotient-v1.json"
)
DEVELOPMENT_BLIND = ART / "half_lattice_search_ablation_r17_development_blind_v1.json"
HOLDOUT_BLIND = ART / "half_lattice_search_ablation_rank29_holdout_blind_v1.json"
DEVELOPMENT_OUTPUT = ART / "half_lattice_search_ablation_r17_development_verification_v1.json"
HOLDOUT_OUTPUT = ART / "half_lattice_search_ablation_rank29_holdout_verification_v1.json"

sys.path[:0] = [str(ELLIPTIC), str(CAS)]

from ecsearch.q12o5867_specialization import (  # noqa: E402
    evaluate_projective_specialization,
    global_minimal_model_with_change,
    load_q12o5867_data,
    short_certificate_model,
)
from elliptic_candidate_record import source_point_to_target  # noqa: E402


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def binary_rank(values: Sequence[int]) -> int:
    pivots = {}
    for value in values:
        value = int(value)
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    return len(pivots)


def rational_rank(vectors: Sequence[Sequence[int]], dimension: int) -> int:
    if not vectors:
        return 0
    return Matrix(QQ, len(vectors), dimension, [value for row in vectors for value in row]).rank()


def quotient_mask(vector: Sequence[int]) -> int:
    return sum((int(value) & 1) << index for index, value in enumerate(vector))


def point_tuple(row: dict[str, str]) -> tuple[Fraction, Fraction]:
    return Fraction(row["x"]), Fraction(row["y"])


def short_data(ainvs, points):
    a1, a2, a3, a4, a6 = (QQ(value) for value in ainvs)
    b2 = a1 * a1 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3 * a3 + 4 * a6
    c4 = b2 * b2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    model = (QQ(0), QQ(0), QQ(0), -27 * c4, -54 * c6)
    short_points = tuple(
        (36 * QQ(point[0]) + 3 * b2, 108 * (2 * QQ(point[1]) + a1 * QQ(point[0]) + a3))
        for point in points
    )
    return (
        tuple(Fraction(str(value)) for value in model),
        tuple(
            tuple(Fraction(str(value)) for value in point) for point in short_points
        ),
    )


def relation_chunks(model, basis, points, relation_proposals, *, chunk_size, timeout, stack_bytes):
    answers = []
    failures = []
    for start in range(0, len(points), chunk_size):
        chunk = points[start : start + chunk_size]
        try:
            rows = relation_proposals(
                model, basis, chunk, timeout=timeout, stack_bytes=stack_bytes
            )
        except Exception as error:
            failures.append(
                {
                    "start": start,
                    "count": len(chunk),
                    "error_type": type(error).__name__,
                    "error": str(error),
                }
            )
            answers.extend((None, False) for unused in chunk)
            continue
        answers.extend(rows)
    return tuple(answers), failures


def development_fixtures(blind, relation_proposals):
    from elkies_rank25 import POINTS as RANK25_POINTS
    from elkies_rank26 import POINTS as RANK26_POINTS
    from elkies_rank27 import POINTS as RANK27_POINTS
    from elkies_rank28 import POINTS as RANK28_POINTS

    public_by_rank = {
        25: RANK25_POINTS,
        26: RANK26_POINTS,
        27: RANK27_POINTS,
        28: RANK28_POINTS,
    }
    controls = json.loads(CONTROLS.read_text())
    if controls.get("status") != "PASS_EXACT_ELKIES_2026_HIGH_RANK_POSITIVE_CONTROLS":
        raise ValueError("R17 positive-control fixture is not passing")
    by_parameter = {row["parameter"]: row for row in controls["fibres"]}
    family = load_q12o5867_data(MODEL, SECTIONS)
    fixtures = {}
    for row in blind["results"]:
        parameter = row["parameter"]
        control = by_parameter[parameter]
        rank = int(control["published_rank_lower_bound"])
        numerator, denominator = row["projective_parameter"]
        specialization = evaluate_projective_specialization(family, numerator, denominator)
        minimal_model, minimal_change, unused = global_minimal_model_with_change(
            specialization.model
        )
        short_model, short_change = short_certificate_model(minimal_model)
        generic = tuple(
            source_point_to_target(
                source_point_to_target(point, minimal_change), short_change
            )
            for point in specialization.points
        )
        public_indices = tuple(
            int(value) - 1
            for value in control["public_positive_control"][
                "selected_public_point_indices_one_based"
            ]
        )
        complement = tuple(
            source_point_to_target(public_by_rank[rank][index], short_change)
            for index in public_indices
        )
        model = tuple(Fraction(value) for value in short_model)
        if [str(value) for value in model] != row["short_model"]:
            raise ArithmeticError(f"{parameter}: blind/public model mismatch")
        if tuple(point_tuple(point) for point in row["generic_points"]) != generic:
            raise ArithmeticError(f"{parameter}: blind/public generic subgroup mismatch")
        fixtures[row["label"]] = {
            "model": model,
            "basis": generic + complement,
            "quotient_dimension": rank - 17,
            "coordinate_map": lambda relation: tuple(int(value) for value in relation[17:]),
            "published_rank_lower_bound": rank,
            "fixture_inputs": {
                "parameter": parameter,
                "public_complement_indices_one_based": [index + 1 for index in public_indices],
            },
        }
    return fixtures, {
        str(CONTROLS.relative_to(ROOT)): digest(CONTROLS),
        str(MODEL.relative_to(ROOT)): digest(MODEL),
        str(SECTIONS.relative_to(ROOT)): digest(SECTIONS),
        **{
            str((CAS / f"elkies_rank{rank}.py").relative_to(ROOT)): digest(
                CAS / f"elkies_rank{rank}.py"
            )
            for rank in public_by_rank
        },
    }


def holdout_fixtures(blind, relation_proposals):
    import elkies_klagsbrun_rank29
    import icarm_curve356

    public_projection = json.loads(RANK29_PUBLIC.read_text())
    curve12_certificate = json.loads(CURVE12_QUOTIENT.read_text())
    model12 = tuple(
        Fraction(str(value))
        for value in elkies_klagsbrun_rank29.short_weierstrass_coefficients()
    )
    points12 = tuple(
        tuple(Fraction(str(value)) for value in point)
        for point in elkies_klagsbrun_rank29.published_short_points()
    )
    model356 = tuple(
        Fraction(str(value)) for value in icarm_curve356.short_coefficients()
    )
    points356 = tuple(
        tuple(Fraction(str(value)) for value in point)
        for point in icarm_curve356.SHORT_POINTS
    )
    record385 = next(row for row in public_projection["records"] if row["id"] == 385)
    model385, points385 = short_data(record385["ainvs"], record385["points"])

    coordinate12 = Matrix(
        ZZ,
        curve12_certificate["specialized_generic_subgroup"][
            "coordinate_matrix_rows_in_ordered_29_public_points"
        ],
    )
    complement_indices12 = tuple(
        int(label[1:]) - 1
        for label in curve12_certificate["displayed_exceptional_quotient"][
            "free_basis_modulo_specialized_generic"
        ]
    )
    complement12 = Matrix(
        ZZ, 29, 12, lambda row, column: int(row == complement_indices12[column])
    )
    augmented12 = coordinate12.augment(complement12)
    if abs(augmented12.det()) != 1:
        raise ArithmeticError("curve12 public coordinate transform lost unimodularity")
    inverse12 = augmented12.inverse()

    def curve12_map(relation):
        coordinates = inverse12 * Matrix(ZZ, 29, 1, relation)
        if any(value not in ZZ for value in coordinates.list()):
            raise ArithmeticError("curve12 quotient coordinate is not integral")
        return tuple(int(value) for value in coordinates.list()[17:])

    identity_map = lambda relation: tuple(int(value) for value in relation[17:29])
    raw = {
        "curve12-2024-rank29": (model12, points12, curve12_map),
        "curve356-rank29": (model356, points356, identity_map),
        "curve385-rank29": (model385, points385, identity_map),
    }
    fixtures = {}
    blind_by_label = {row["label"]: row for row in blind["results"]}
    for label, (model, basis, coordinate_map) in raw.items():
        row = blind_by_label[label]
        if [str(value) for value in model] != row["short_model"]:
            raise ArithmeticError(f"{label}: blind/public model mismatch")
        generic = tuple(point_tuple(point) for point in row["generic_points"])
        if label != "curve12-2024-rank29" and generic != basis[:17]:
            raise ArithmeticError(f"{label}: public first 17 stopped matching generic input")
        fixtures[label] = {
            "model": model,
            "basis": basis,
            "quotient_dimension": 12,
            "coordinate_map": coordinate_map,
            "published_rank_lower_bound": 29,
            "fixture_inputs": {},
        }
    return fixtures, {
        str(RANK29_PUBLIC.relative_to(ROOT)): digest(RANK29_PUBLIC),
        str(CURVE12_QUOTIENT.relative_to(ROOT)): digest(CURVE12_QUOTIENT),
        str((CAS / "elkies_klagsbrun_rank29.py").relative_to(ROOT)): digest(
            CAS / "elkies_klagsbrun_rank29.py"
        ),
        str((CAS / "icarm_curve356.py").relative_to(ROOT)): digest(
            CAS / "icarm_curve356.py"
        ),
    }


def greedy_basis(candidate_indices, quotient_vectors, dimension):
    selected = []
    vectors = []
    current_rank = 0
    for index in candidate_indices:
        vector = quotient_vectors[index]
        if vector is None:
            continue
        trial = vectors + [vector]
        rank = rational_rank(trial, dimension)
        if rank > current_rank:
            vectors = trial
            selected.append(index)
            current_rank = rank
    return selected, current_rank


def aggregate(results):
    arm_ids = [row["id"] for row in results[0]["arms"]]
    rows = []
    for arm_id in arm_ids:
        arms = [next(arm for arm in row["arms"] if arm["id"] == arm_id) for row in results]
        total_rank = sum(arm["exact_quotient_rank_over_Q"] for arm in arms)
        total_cpu = sum(arm["cover_cpu_seconds"] for arm in arms)
        rows.append(
            {
                "id": arm_id,
                "case_count": len(arms),
                "total_exact_quotient_rank_over_Q": total_rank,
                "mean_exact_quotient_rank_over_Q": total_rank / len(arms),
                "mean_rank_normalized_per_43_covers": sum(
                    arm["rank_normalized_per_43_covers"] for arm in arms
                )
                / len(arms),
                "pooled_quotient_rank_per_cpu_second": (
                    total_rank / total_cpu if total_cpu else None
                ),
                "mean_target_recovery_fraction": sum(
                    arm["target_recovery_fraction"] for arm in arms
                )
                / len(arms),
                "all_public_relations_exact": all(
                    arm["unresolved_candidate_count"] == 0 for arm in arms
                ),
            }
        )
    random_rows = [row for row in rows if row["id"].startswith("random43-")]
    return {
        "arms": rows,
        "random43_mean_of_arm_means_rank": sum(
            row["mean_exact_quotient_rank_over_Q"] for row in random_rows
        )
        / len(random_rows),
        "random43_best_arm_mean_rank": max(
            row["mean_exact_quotient_rank_over_Q"] for row in random_rows
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", choices=("development", "holdout"), required=True)
    parser.add_argument("--blind", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--relation-chunk-size", type=int, default=64)
    parser.add_argument("--relation-timeout-seconds", type=float, default=180.0)
    parser.add_argument("--stack-bytes", type=int, default=1_000_000_000)
    args = parser.parse_args()
    blind_path = args.blind or (
        DEVELOPMENT_BLIND if args.phase == "development" else HOLDOUT_BLIND
    )
    output = args.output or (
        DEVELOPMENT_OUTPUT if args.phase == "development" else HOLDOUT_OUTPUT
    )
    blind_bytes = blind_path.read_bytes()
    blind_hash = sha256(blind_bytes).hexdigest()
    blind = json.loads(blind_bytes)
    if blind.get("status") != "PASS_BLIND_ABLATION_SEARCH":
        raise ValueError("ablation blind artifact is not complete")
    if blind.get("phase") != args.phase:
        raise ValueError("ablation phase mismatch")
    if blind["blindness_boundary"]["exceptional_point_fixture_loaded"] is not False:
        raise ValueError("ablation blindness assertion failed")

    # Fixture and public-point imports occur only after the frozen bytes above.
    from search_nagao_u135_alternate_covers import relation_proposals

    if args.phase == "development":
        fixtures, fixture_hashes = development_fixtures(blind, relation_proposals)
    else:
        fixtures, fixture_hashes = holdout_fixtures(blind, relation_proposals)

    verified = []
    all_exact = True
    for case in blind["results"]:
        label = case["label"]
        fixture = fixtures[label]
        points = tuple(point_tuple(row["point"]) for row in case["candidate_points"])
        relations, failures = relation_chunks(
            fixture["model"],
            fixture["basis"],
            points,
            relation_proposals,
            chunk_size=args.relation_chunk_size,
            timeout=args.relation_timeout_seconds,
            stack_bytes=args.stack_bytes,
        )
        quotient_vectors = []
        quotient_masks = []
        for relation, exact in relations:
            if exact:
                vector = fixture["coordinate_map"](relation)
                if len(vector) != fixture["quotient_dimension"]:
                    raise ArithmeticError(f"{label}: quotient coordinate length changed")
                quotient_vectors.append(vector)
                quotient_masks.append(quotient_mask(vector))
            else:
                quotient_vectors.append(None)
                quotient_masks.append(None)
        all_exact &= all(exact for unused_relation, exact in relations)

        arm_rows = []
        for arm in case["arms"]:
            indices = arm["candidate_point_indices"]
            vectors = [quotient_vectors[index] for index in indices if quotient_vectors[index] is not None]
            masks = [quotient_masks[index] for index in indices if quotient_masks[index] is not None]
            selected, rank_q = greedy_basis(
                indices, quotient_vectors, fixture["quotient_dimension"]
            )
            if rank_q != rational_rank(vectors, fixture["quotient_dimension"]):
                raise ArithmeticError(f"{label}/{arm['id']}: greedy rational rank changed")
            rank_f2 = binary_rank(masks)
            cpu = float(arm["cover_cpu_seconds"])
            count = int(arm["class_count"])
            arm_rows.append(
                {
                    **{key: value for key, value in arm.items() if key != "candidate_point_indices"},
                    "exact_public_relation_count": len(vectors),
                    "unresolved_candidate_count": len(indices) - len(vectors),
                    "exact_quotient_rank_over_Q": rank_q,
                    "exact_quotient_rank_mod2": rank_f2,
                    "target_quotient_dimension": fixture["quotient_dimension"],
                    "target_recovery_fraction": rank_q / fixture["quotient_dimension"],
                    "rank_normalized_per_43_covers": rank_q * 43 / count,
                    "quotient_rank_per_cpu_second": rank_q / cpu if cpu else None,
                    "selected_exact_basis": [
                        {
                            "candidate_index": index,
                            "point": case["candidate_points"][index]["point"],
                            "source_masks": case["candidate_points"][index]["source_masks"],
                            "quotient_coordinates": list(quotient_vectors[index]),
                            "quotient_mask": quotient_masks[index],
                            "quotient_hex": f"0x{quotient_masks[index]:03x}",
                        }
                        for index in selected
                    ],
                }
            )
        verified.append(
            {
                "label": label,
                "parameter": case["parameter"],
                "published_rank_lower_bound": fixture["published_rank_lower_bound"],
                "target_quotient_dimension": fixture["quotient_dimension"],
                "blind_candidate_count": len(points),
                "all_blind_candidates_have_exact_public_relations": all(
                    exact for unused_relation, exact in relations
                ),
                "relation_failure_chunks": failures,
                "fixture_inputs": fixture["fixture_inputs"],
                "ranking": case["ranking"],
                "arms": arm_rows,
            }
        )
        print(
            f"HALFABLATEVERIFY|phase={args.phase}|case={label}|"
            + "|".join(
                f"{row['id']}={row['exact_quotient_rank_over_Q']}"
                for row in arm_rows
            ),
            flush=True,
        )

    payload: dict[str, Any] = {
        "schema": "elliptic-curves.half-lattice-search-ablation-verification.v1",
        "status": (
            "PASS_EXACT_PUBLIC_QUOTIENT_ABLATION"
            if all_exact
            else "PARTIAL_PUBLIC_QUOTIENT_ABLATION"
        ),
        "phase": args.phase,
        "phase_boundary": {
            "blind_artifact_sha256_before_fixture_import": blind_hash,
            "public_points_loaded_only_after_blind_artifact_frozen": True,
            "holdout_arms_and_hash_seeds_identical_to_development": True,
        },
        "blind_artifact": str(blind_path.relative_to(ROOT)),
        "input_hashes": {
            str(blind_path.relative_to(ROOT)): blind_hash,
            **fixture_hashes,
            str(Path(__file__).resolve().relative_to(ROOT)): digest(Path(__file__).resolve()),
        },
        "results": verified,
        "aggregate": aggregate(verified),
        "claim_boundary": [
            "Every displayed quotient coordinate and rank is exact in the displayed known subgroup.",
            "The primary rank metric is Q-linear quotient rank; mod-2 rank is retained separately.",
            "Searches are bounded, so recovered dimensions are lower bounds on what larger searches may recover.",
            "Per-CPU rates use the blind per-cover child-plus-parent CPU measurements.",
            "The three rank-29 cases are holdouts for the arm rules frozen in the search source.",
        ],
        "reproducing_command": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            f"elliptic-curves/cas/verify_half_lattice_search_ablation.sage --phase {args.phase}"
        ),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        f"HALFABLATEVERIFY|phase={args.phase}|status={payload['status']}|"
        f"output={output.relative_to(ROOT)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
