#!/usr/bin/env sage -python
"""Verification-only replay for the frozen 273/302/245 held-out experiment.

The blind artifact is required to exist before this program imports any full
public point fixture.  Exact height-pairing proposals are replayed in Fraction
group arithmetic by ``relation_proposals``.  Since each starting subgroup is
spanned by a coordinate subset of the displayed public basis, restriction of
an exact relation to the omitted coordinates gives its held-out quotient
class modulo 2.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
BLIND = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "half_lattice_heldout_273_302_blind_v1.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "half_lattice_heldout_273_302_verification_v1.json"
)
sys.path.insert(0, str(CAS))


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


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


def point_from_record(record):
    return Fraction(record["x"]), Fraction(record["y"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blind", type=Path, default=BLIND)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--chunk-size", type=int, default=96)
    parser.add_argument("--timeout-seconds", type=float, default=180.0)
    parser.add_argument("--stack-bytes", type=int, default=1_000_000_000)
    args = parser.parse_args()

    # Freeze the verification boundary before loading the fixture modules.
    blind_bytes = args.blind.read_bytes()
    blind_hash = sha256(blind_bytes).hexdigest()
    blind = json.loads(blind_bytes)
    if blind.get("status") != "PASS_BOUNDED_HELDOUT_SUBGROUP_SEARCH":
        raise ValueError("the blind artifact is absent or did not finish")
    if blind["blindness_boundary"]["search_loaded_heldout_point_coordinates"] is not False:
        raise ValueError("the blind artifact did not preserve the held-out boundary")

    import icarm_curve245
    import icarm_curve273
    import icarm_curve302
    from search_nagao_u135_alternate_covers import relation_proposals

    fixtures = {
        "curve273": icarm_curve273,
        "curve302": icarm_curve302,
        "curve245-adverse-control": icarm_curve245,
    }
    verification_rows = []
    all_exact = True
    for result in blind["results"]:
        module = fixtures[result["curve"]]
        full_basis = tuple(module.SHORT_POINTS)
        model = tuple(module.short_coefficients())
        included = {index - 1 for index in result["included_public_indices_one_based"]}
        omitted = tuple(index for index in range(len(full_basis)) if index not in included)
        target_rank = len(omitted)
        candidate_records = result["blind_result"]["candidate_points"]
        candidates = tuple(point_from_record(row["point"]) for row in candidate_records)
        relation_rows = []
        quotient_masks = []
        processed = 0
        # The small 273/302 candidate sets are verified in full.  Curve 245
        # returns thousands of easy combinations, so stop only once its entire
        # possible displayed held-out quotient has been certified.
        while processed < len(candidates):
            stop = min(len(candidates), processed + args.chunk_size)
            proposals = relation_proposals(
                model,
                full_basis,
                candidates[processed:stop],
                timeout=args.timeout_seconds,
                stack_bytes=args.stack_bytes,
            )
            for local_index, (relation, exact) in enumerate(proposals):
                candidate_index = processed + local_index
                all_exact &= exact
                if exact:
                    quotient_mask = sum(
                        (int(relation[index]) & 1) << offset
                        for offset, index in enumerate(omitted)
                    )
                    quotient_masks.append(quotient_mask)
                else:
                    quotient_mask = None
                relation_rows.append(
                    {
                        "candidate_index": candidate_index,
                        "exact_relation_in_displayed_public_basis": exact,
                        "heldout_quotient_mask": quotient_mask,
                        "heldout_quotient_hex": (
                            f"0x{quotient_mask:0{max(1, (target_rank + 3) // 4)}x}"
                            if quotient_mask is not None
                            else None
                        ),
                        "relation": list(relation) if exact else None,
                        "selected_by_blind_finite_reduction": candidate_records[candidate_index][
                            "selected_for_independent_quotient_basis"
                        ],
                        "source_half_class_masks": candidate_records[candidate_index]["source_masks"],
                    }
                )
            processed = stop
            if result["curve"] == "curve245-adverse-control" and binary_rank(quotient_masks) == target_rank:
                break

        exact_masks = [
            row["heldout_quotient_mask"]
            for row in relation_rows
            if row["exact_relation_in_displayed_public_basis"]
        ]
        heldout_rank = binary_rank(exact_masks)
        blind_selected_masks = [
            row["heldout_quotient_mask"]
            for row in relation_rows
            if row["exact_relation_in_displayed_public_basis"]
            and row["selected_by_blind_finite_reduction"]
        ]
        verification_rows.append(
            {
                "curve": result["curve"],
                "configuration": result["configuration"],
                "dimension": result["dimension"],
                "displayed_public_rank": len(full_basis),
                "heldout_dimension": target_rank,
                "candidate_count_in_blind_artifact": len(candidates),
                "candidate_count_verified": processed,
                "all_processed_candidates_have_exact_public_basis_relations": all(
                    row["exact_relation_in_displayed_public_basis"] for row in relation_rows
                ),
                "exact_heldout_quotient_rank_recovered": heldout_rank,
                "exact_heldout_recovery_fraction": f"{heldout_rank}/{target_rank}",
                "blind_selected_exact_heldout_rank": binary_rank(blind_selected_masks),
                "blind_finite_reduction_field_valid": result["curve"] != "curve245-adverse-control",
                "blind_finite_reduction_note": (
                    None
                    if result["curve"] != "curve245-adverse-control"
                    else "INVALID: simultaneous denominators of thousands of candidates excluded every tested small reduction prime; the negative blind gain is not a rank statement"
                ),
                "relations": relation_rows,
            }
        )
        print(
            f"HELDOUTVERIFY|case={result['curve']}/{result['configuration']}|"
            f"exact_rank={heldout_rank}/{target_rank}|processed={processed}/{len(candidates)}",
            flush=True,
        )

    payload = {
        "schema": "elliptic-curves.half-lattice-heldout-273-302-verification.v1",
        "status": "PASS_EXACT_HELDOUT_RELATION_VERIFICATION" if all_exact else "PARTIAL",
        "blind_artifact": str(args.blind.relative_to(ROOT)),
        "blind_artifact_sha256_before_fixture_import": blind_hash,
        "verification_boundary": {
            "full_public_point_fixtures_loaded_only_after_blind_artifact_frozen": True,
            "curve245_verification_may_stop_after_full_possible_heldout_span": True,
        },
        "results": verification_rows,
        "claim_boundary": [
            "Every exact relation was replayed in exact Fraction elliptic-curve arithmetic.",
            "Recovered ranks are ranks in the displayed public free lattice modulo the selected coordinate subgroup and modulo 2.",
            "They do not prove saturation in the full Mordell-Weil group or any family provenance.",
            "The blind finite-reduction gain fields for curve245 are invalid diagnostics, explicitly superseded here.",
        ],
        "input_hashes": {
            str(args.blind.relative_to(ROOT)): blind_hash,
            str((CAS / "icarm_curve245.py").relative_to(ROOT)): digest(CAS / "icarm_curve245.py"),
            str((CAS / "icarm_curve273.py").relative_to(ROOT)): digest(CAS / "icarm_curve273.py"),
            str((CAS / "icarm_curve302.py").relative_to(ROOT)): digest(CAS / "icarm_curve302.py"),
            str(Path(__file__).resolve().relative_to(ROOT)): digest(Path(__file__).resolve()),
        },
        "reproducing_command": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elliptic-curves/cas/verify_half_lattice_heldout_subgroups.sage"
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        "HELDOUTVERIFY|status=" + payload["status"] + f"|output={args.output.relative_to(ROOT)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
