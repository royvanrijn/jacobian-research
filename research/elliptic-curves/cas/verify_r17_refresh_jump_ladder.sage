#!/usr/bin/env sage-python
"""Open the public complements only after freezing the blind ladder result."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys

from sage.all import Matrix, QQ


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
BLIND = ROOT / "artifacts/generated-results/elliptic-curves/r17_refresh_jump_ladder_blind_v1.json"
PROTOCOL = ROOT / "artifacts/generated-results/elliptic-curves/r17_refresh_jump_ladder_protocol_v1.json"
OVERVIEW = ROOT / "artifacts/generated-results/elliptic-curves/icarm_curve_refresh_475_573_overview_v1.json"
QUOTIENTS = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-refresh-priority-quotients-v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/r17_refresh_jump_ladder_verification_v1.json"

sys.path[:0] = [str(ROOT / "elliptic-curves"), str(CAS)]


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def short_model_and_points(public_record):
    a1, a2, a3, a4, a6 = (Fraction(value) for value in public_record["ainvs"])
    b2 = a1 * a1 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3 * a3 + 4 * a6
    c4 = b2 * b2 - 24 * b4
    c6 = -(b2**3) + 36 * b2 * b4 - 216 * b6
    model = (Fraction(0), Fraction(0), Fraction(0), -c4 / 48, -c6 / 864)
    points = tuple(
        (
            Fraction(x_value) + b2 / 12,
            Fraction(y_value) + (a1 * Fraction(x_value) + a3) / 2,
        )
        for x_value, y_value in public_record["points"]
    )
    return model, points


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blind", type=Path, default=BLIND)
    parser.add_argument("--protocol", type=Path, default=PROTOCOL)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--timeout-seconds", type=float, default=180.0)
    parser.add_argument("--stack-bytes", type=int, default=1_000_000_000)
    args = parser.parse_args()

    # This digest and all boundary checks happen before either truth artifact is
    # opened or its helper modules are imported.
    blind_bytes = args.blind.read_bytes()
    blind_hash_before_truth = sha256(blind_bytes).hexdigest()
    blind = json.loads(blind_bytes)
    protocol = json.loads(args.protocol.read_text())
    if blind.get("status") != "PASS_COMPLETE_BLIND_RECOVERY_BEFORE_PUBLIC_COMPLEMENT":
        raise ArithmeticError("the blind ladder result is not complete")
    if blind["blindness_boundary"]["public_complement_loaded"]:
        raise ArithmeticError("the blind runner crossed the complement boundary")
    if blind["protocol_sha256"] != digest(args.protocol):
        raise ArithmeticError("the blind runner used another protocol")
    if protocol.get("status") != "FROZEN_BEFORE_BLIND_RECOVERY":
        raise ArithmeticError("the protocol was not frozen")

    # Public complement opening begins here.
    overview = json.loads(OVERVIEW.read_text())
    quotient = json.loads(QUOTIENTS.read_text())
    from search_nagao_u135_alternate_covers import relation_proposals

    public = {int(row["id"]): row for row in overview["snapshot"]["records"]}
    truth = {int(row["curve_id"]): row for row in quotient["fibres"]}
    if set(truth) != {int(row["curve_id"]) for row in blind["results"]}:
        raise ArithmeticError("blind and truth case inventories differ")

    results = []
    all_exact = True
    for blind_row in blind["results"]:
        curve_id = int(blind_row["curve_id"])
        model, public_points = short_model_and_points(public[curve_id])
        blind_model = tuple(Fraction(value) for value in next(
            row for row in json.loads(
                (ROOT / "elliptic-curves/data/r17_refresh_jump_ladder_blind_inputs_v1.json").read_text()
            )["cases"] if int(row["curve_id"]) == curve_id
        )["short_model"])
        if model != blind_model:
            raise ArithmeticError(f"curve {curve_id}: public and blind models differ")
        final_basis = tuple(
            (Fraction(row["x"]), Fraction(row["y"]))
            for row in blind_row["final_basis"]
        )
        proposals = relation_proposals(
            model,
            public_points,
            final_basis,
            timeout=args.timeout_seconds,
            stack_bytes=args.stack_bytes,
        )
        relation_rows = []
        quotient_rows = []
        for basis_index, (relation, exact) in enumerate(proposals):
            all_exact &= bool(exact)
            if exact:
                relation = tuple(map(int, relation))
                quotient_coordinates = relation[17:]
                quotient_rows.append(quotient_coordinates)
            else:
                relation = None
                quotient_coordinates = None
            relation_rows.append(
                {
                    "blind_basis_index": basis_index,
                    "exact_relation_in_opened_public_basis": bool(exact),
                    "public_basis_coordinates": list(relation) if relation is not None else None,
                    "P18_and_beyond_coordinates": (
                        list(quotient_coordinates)
                        if quotient_coordinates is not None
                        else None
                    ),
                }
            )
        exact_public_quotient_rank = (
            int(Matrix(QQ, quotient_rows).rank()) if quotient_rows else 0
        )
        blind_gain = int(
            blind_row["exact_quotient_rank_recovered_before_public_complement"]
        )
        row_exact = all(row["exact_relation_in_opened_public_basis"] for row in relation_rows)
        if row_exact and exact_public_quotient_rank != blind_gain:
            raise ArithmeticError(
                f"curve {curve_id}: blind gain {blind_gain} != public quotient rank "
                f"{exact_public_quotient_rank}"
            )
        true_jump = int(
            truth[curve_id]["displayed_exceptional_quotient"]["free_rank"]
        )
        if blind_gain > true_jump and row_exact:
            raise ArithmeticError(f"curve {curve_id}: recovered rank exceeds displayed jump")
        results.append(
            {
                "curve_id": curve_id,
                "representative_class": blind_row["representative_class"],
                "exact_quotient_rank_recovered_before_public_complement": blind_gain,
                "all_final_blind_basis_points_in_opened_public_subgroup": row_exact,
                "exact_rank_of_blind_basis_images_in_opened_P18_and_beyond_complement": (
                    exact_public_quotient_rank if row_exact else None
                ),
                "true_displayed_jump_opened_after_blind_freeze": true_jump,
                "recovery_fraction": f"{blind_gain}/{true_jump}",
                "initial_exact_recovered_rank": int(
                    blind_row["initial"]["exact_quotient_rank_recovered"]
                ),
                "adaptive_incremental_exact_recovered_rank": int(
                    blind_row.get("adaptive", {}).get(
                        "exact_incremental_quotient_rank_recovered", 0
                    )
                ),
                "attempted_chart_count": int(blind_row["attempted_chart_count"]),
                "timeout_chart_count": int(blind_row.get("timeout_chart_count", 0)),
                "pari_failure_chart_count": int(
                    blind_row.get("pari_failure_chart_count", 0)
                ),
                "opened_public_relations": relation_rows,
            }
        )
        print(
            f"R17JUMPLADDERVERIFY|curve={curve_id}|recovered={blind_gain}|"
            f"jump={true_jump}|relations={'PASS' if row_exact else 'UNKNOWN'}",
            flush=True,
        )

    status = (
        "PASS_ALL_BLIND_RECOVERY_RANKS_EXACT_IN_OPENED_DISPLAYED_QUOTIENTS"
        if all_exact
        else "UNKNOWN_BLIND_POINT_OUTSIDE_OR_UNRESOLVED_IN_DISPLAYED_SUBGROUP"
    )
    payload = {
        "schema": "elliptic-curves.r17-refresh-jump-ladder-verification.v1",
        "status": status,
        "phase_boundary": {
            "blind_artifact_sha256_before_public_complement_import": blind_hash_before_truth,
            "public_complement_opened_only_after_blind_status_and_hash_were_fixed": True,
        },
        "results": results,
        "pre_search_exclusion": protocol["pre_search_exclusion"],
        "inputs": {
            relative(args.blind): blind_hash_before_truth,
            relative(args.protocol): digest(args.protocol),
            relative(OVERVIEW): digest(OVERVIEW),
            relative(QUOTIENTS): digest(QUOTIENTS),
            relative(Path(__file__).resolve()): digest(Path(__file__).resolve()),
        },
        "generation": {
            "command": (
                "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
                "elliptic-curves/cas/verify_r17_refresh_jump_ladder.sage"
            )
        },
        "claim_boundary": [
            "Every displayed-subgroup relation is replayed by exact rational group law.",
            "Recovered ranks are exact ranks of the discovered subgroup modulo specialized MW17.",
            "True jumps are ranks of certified displayed-subgroup quotients, not full Mordell-Weil rank assertions.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"R17JUMPLADDERVERIFY|status={status}|output={relative(args.output)}")


if __name__ == "__main__":
    main()
