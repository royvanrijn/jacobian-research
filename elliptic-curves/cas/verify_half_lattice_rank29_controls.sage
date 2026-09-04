#!/usr/bin/env sage -python
"""Verification-only quotient audit for the blind rank-29 control replay."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys

from sage.all import Matrix, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
BLIND = ROOT / "artifacts/generated-results/elliptic-curves/half_lattice_rank29_controls_blind_v1.json"
PUBLIC = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-public-fibres-v1.json"
CURVE12 = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-curve12-alternate-q80-quotient-v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/half_lattice_rank29_controls_verification_v1.json"
sys.path.insert(0, str(CAS))


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def binary_rank(values) -> int:
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
    return model, short_points


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blind", type=Path, default=BLIND)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--timeout-seconds", type=float, default=180.0)
    parser.add_argument("--stack-bytes", type=int, default=1_000_000_000)
    args = parser.parse_args()

    blind_bytes = args.blind.read_bytes()
    blind_hash = sha256(blind_bytes).hexdigest()
    blind = json.loads(blind_bytes)
    if blind.get("status") != "PASS_BOUNDED_BLIND_RANK29_CONTROL_SEARCH":
        raise ValueError("rank29 blind artifact did not finish")
    if blind["blindness_boundary"]["search_loaded_exceptional_coordinates"] is not False:
        raise ValueError("rank29 blind boundary was not preserved")

    # Fixture imports occur only after the blind bytes and digest are frozen.
    import elkies_klagsbrun_rank29
    import icarm_curve356
    from search_nagao_u135_alternate_covers import relation_proposals

    public_projection = json.loads(PUBLIC.read_text())
    curve12_certificate = json.loads(CURVE12.read_text())
    model12 = tuple(Fraction(str(value)) for value in elkies_klagsbrun_rank29.short_weierstrass_coefficients())
    points12 = tuple(tuple(Fraction(str(value)) for value in point) for point in elkies_klagsbrun_rank29.published_short_points())
    model356 = tuple(Fraction(str(value)) for value in icarm_curve356.short_coefficients())
    points356 = tuple(tuple(Fraction(str(value)) for value in point) for point in icarm_curve356.SHORT_POINTS)
    record385 = next(record for record in public_projection["records"] if record["id"] == 385)
    model385_q, points385_q = short_data(record385["ainvs"], record385["points"])
    model385 = tuple(Fraction(str(value)) for value in model385_q)
    points385 = tuple(tuple(Fraction(str(value)) for value in point) for point in points385_q)

    coordinate12 = Matrix(
        ZZ,
        curve12_certificate["specialized_generic_subgroup"]["coordinate_matrix_rows_in_ordered_29_public_points"],
    )
    complement_indices12 = tuple(
        int(label[1:]) - 1
        for label in curve12_certificate["displayed_exceptional_quotient"]["free_basis_modulo_specialized_generic"]
    )
    complement12 = Matrix(
        ZZ, 29, 12, lambda row, column: int(row == complement_indices12[column])
    )
    augmented12 = coordinate12.augment(complement12)
    if abs(augmented12.det()) != 1:
        raise ArithmeticError("curve12 public quotient coordinates lost unimodularity")
    inverse12 = augmented12.inverse()

    fixtures = {
        "curve12-2024-rank29": (model12, points12),
        "curve356-rank29": (model356, points356),
        "curve385-rank29": (model385, points385),
    }
    results = []
    all_exact = True
    for result in blind["results"]:
        model, public_points = fixtures[result["label"]]
        candidate_records = result["blind_result"]["candidate_points"]
        candidates = tuple(
            (Fraction(row["point"]["x"]), Fraction(row["point"]["y"]))
            for row in candidate_records
        )
        proposals = relation_proposals(
            model, public_points, candidates,
            timeout=args.timeout_seconds, stack_bytes=args.stack_bytes,
        ) if candidates else ()
        relation_rows = []
        quotient_masks = []
        selected_masks = []
        for index, ((relation, exact), blind_row) in enumerate(zip(proposals, candidate_records)):
            all_exact &= exact
            if exact and result["label"] == "curve12-2024-rank29":
                coordinates = inverse12 * Matrix(ZZ, 29, 1, relation)
                if any(value not in ZZ for value in coordinates.list()):
                    raise ArithmeticError("curve12 quotient coordinate stopped being integral")
                quotient = tuple(int(value) for value in coordinates.list()[17:])
            elif exact:
                quotient = tuple(int(value) for value in relation[17:29])
            else:
                quotient = None
            mask = (
                sum((value & 1) << offset for offset, value in enumerate(quotient))
                if quotient is not None else None
            )
            if mask is not None:
                quotient_masks.append(mask)
                if blind_row["selected_for_independent_quotient_basis"]:
                    selected_masks.append(mask)
            relation_rows.append(
                {
                    "candidate_index": index,
                    "exact_relation_in_displayed_public_basis": exact,
                    "relation": list(relation) if exact else None,
                    "quotient_coordinates": list(quotient) if quotient is not None else None,
                    "quotient_mask": mask,
                    "quotient_hex": f"0x{mask:03x}" if mask is not None else None,
                    "source_half_class_masks": blind_row["source_masks"],
                    "selected_by_blind_finite_reduction": blind_row[
                        "selected_for_independent_quotient_basis"
                    ],
                }
            )
        exact_rank = binary_rank(quotient_masks)
        selected_rank = binary_rank(selected_masks)
        if result["blind_result"]["finite_reduction_certificate_valid"]:
            if exact_rank != result["blind_result"]["finite_mod2_quotient_gain"]:
                raise ArithmeticError(f"{result['label']}: blind and fixture quotient ranks differ")
            if selected_rank != exact_rank:
                raise ArithmeticError(f"{result['label']}: blind selected basis lost quotient rank")
        results.append(
            {
                "label": result["label"],
                "displayed_public_rank": 29,
                "displayed_exceptional_quotient_rank": 12,
                "all_blind_candidates_have_exact_public_basis_relations": all(
                    row["exact_relation_in_displayed_public_basis"] for row in relation_rows
                ),
                "exact_exceptional_quotient_rank_recovered": exact_rank,
                "recovery_fraction": f"{exact_rank}/12",
                "blind_selected_exact_quotient_rank": selected_rank,
                "relations": relation_rows,
            }
        )
        print(f"RANK29HALFVERIFY|case={result['label']}|rank={exact_rank}/12", flush=True)

    payload = {
        "schema": "elliptic-curves.half-lattice-rank29-controls-verification.v1",
        "status": "PASS_EXACT_RANK29_CONTROL_RELATIONS" if all_exact else "PARTIAL",
        "blind_artifact": str(args.blind.relative_to(ROOT)),
        "blind_artifact_sha256_before_fixture_import": blind_hash,
        "verification_boundary": "Full displayed point fixtures were loaded only after the blind artifact was frozen.",
        "results": results,
        "claim_boundary": [
            "All recorded relations are replayed by exact Fraction group arithmetic.",
            "Quotient ranks are exact inside the independently certified displayed rank-29 subgroups.",
            "Search incompleteness is bounded; no exact Mordell-Weil rank is asserted.",
        ],
        "input_hashes": {
            str(args.blind.relative_to(ROOT)): blind_hash,
            str(PUBLIC.relative_to(ROOT)): digest(PUBLIC),
            str(CURVE12.relative_to(ROOT)): digest(CURVE12),
            str((CAS / "elkies_klagsbrun_rank29.py").relative_to(ROOT)): digest(CAS / "elkies_klagsbrun_rank29.py"),
            str((CAS / "icarm_curve356.py").relative_to(ROOT)): digest(CAS / "icarm_curve356.py"),
            str(Path(__file__).resolve().relative_to(ROOT)): digest(Path(__file__).resolve()),
        },
        "reproducing_command": "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python elliptic-curves/cas/verify_half_lattice_rank29_controls.sage",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"RANK29HALFVERIFY|status={payload['status']}|output={args.output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
