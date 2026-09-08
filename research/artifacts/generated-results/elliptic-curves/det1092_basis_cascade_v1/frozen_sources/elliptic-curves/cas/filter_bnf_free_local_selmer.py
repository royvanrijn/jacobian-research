#!/usr/bin/env python3
"""Apply certified local Kummer conditions to a BNF-free image manifest.

``evaluate_bnf_free_signature_map.py`` evaluates global squareclass
candidates in a fixed product of local squareclass coordinates.  This script
is the next, deliberately separate step: it tests those local vectors against
a supplied basis for the *full local Kummer image*.  Known rational points are
checked to lie in that space, but are never used as a substitute for it.

The output remains a bookkeeping manifest.  A global S-class quotient bound
and the norm condition of the descent are still necessary before its surviving
directions are a certified residual 2-Selmer basis.

As a separate strict route, an audit from
``audit_bnf_free_two_cover_reduction.py`` can filter candidate covers at its
selected finite places.  Only candidates with a proved ``Q_p`` point at every
listed place survive that route; an obstruction rejects and an inconclusive
place is kept separate rather than admitted as a local-Selmer class.

``audit_bnf_free_local_kummer_coverage.py`` supplies another strict partial
route: at each odd/real place where known points are certified to span the
full local Kummer image, the candidate must reduce into their local span.
Other places remain explicitly unresolved.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


CAS_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(CAS_ROOT))

from residual_selmer_quotient import F2Error, _integer, audit_manifest


PROTOCOL = "BNFFREELOCAL"


def f2_reduce(mask: int, basis: list[int]) -> int:
    """Reduce a packed vector against an echelon GF(2) basis."""

    pivots: dict[int, int] = {}
    for vector in basis:
        reduced = vector
        while reduced:
            pivot = reduced.bit_length() - 1
            if pivot in pivots:
                reduced ^= pivots[pivot]
            else:
                pivots[pivot] = reduced
                break
    reduced = mask
    while reduced:
        pivot = reduced.bit_length() - 1
        if pivot not in pivots:
            break
        reduced ^= pivots[pivot]
    return reduced


def validate_condition_map(record: dict, local_dimension: int) -> list[int]:
    if record.get("schema") != "elliptic-curves.bnf-free-local-kummer-map.v1":
        raise F2Error("unexpected local Kummer-map schema")
    if _integer(record.get("local_dimension"), "local_dimension") != local_dimension:
        raise F2Error("local Kummer-map dimension does not match the image manifest")
    raw_basis = record.get("allowed_local_images")
    if not isinstance(raw_basis, list):
        raise F2Error("allowed_local_images must be a list")
    limit = 1 << local_dimension
    basis = [_integer(value, "allowed_local_images entry") for value in raw_basis]
    if any(value < 0 or value >= limit for value in basis):
        raise F2Error("allowed local image is outside the local coordinate space")
    return basis


def cover_audit_index(record: dict) -> dict[str, dict]:
    if record.get("protocol") != "BNFFREECOVERLOCAL-v1":
        raise F2Error("unexpected BNF-free cover-local audit")
    covers = record.get("covers")
    if not isinstance(covers, list):
        raise F2Error("cover-local audit lacks a cover list")
    index = {}
    for cover in covers:
        if not isinstance(cover, dict) or "label" not in cover:
            raise F2Error("cover-local audit contains an invalid cover record")
        label = str(cover["label"])
        if label in index:
            raise F2Error("cover-local audit has duplicate labels")
        places = cover.get("finite_places")
        if not isinstance(places, list) or not places:
            raise F2Error("cover-local audit must record at least one finite place per cover")
        index[label] = cover
    return index


def cover_audit_filter(images: dict, cover_audit: dict):
    """Keep only candidates proved locally soluble at every audited finite place."""
    index = cover_audit_index(cover_audit)
    survivors = []
    obstructed = []
    inconclusive = []
    for candidate in images.get("candidate_images", []):
        if not isinstance(candidate, dict):
            raise F2Error("candidate_images must contain objects")
        label = str(candidate["label"])
        cover = index.get(label)
        if cover is None:
            inconclusive.append({"label": label, "reason": "MISSING_COVER_LOCAL_AUDIT"})
            continue
        classifications = [
            str(place.get("classification", ""))
            for place in cover["finite_places"]
            if isinstance(place, dict)
        ]
        if any(item.startswith("PROVED_NO_QP_POINT") for item in classifications):
            obstructed.append({"label": label, "finite_place_classifications": classifications})
        elif classifications and all(item.startswith("PROVED_QP_POINT") for item in classifications):
            survivors.append(candidate)
        else:
            inconclusive.append({"label": label, "finite_place_classifications": classifications})
    return survivors, obstructed, inconclusive


def coverage_filter(images: dict, coverage: dict):
    if coverage.get("protocol") != "BNFFREELOCALCOVERAGE-v1":
        raise F2Error("unexpected BNF-free local-coverage audit")
    local_dimension = _integer(images["local_dimension"], "local_dimension")
    coverage_dimension = _integer(
        coverage.get("signature_local_dimension"), "coverage signature_local_dimension"
    )
    if local_dimension < coverage_dimension:
        raise F2Error("image manifest has fewer local coordinates than the coverage audit")
    covered_known = coverage.get("known_mw_local_images")
    if not isinstance(covered_known, list):
        raise F2Error("local-coverage audit lacks known MW local images")
    image_known = images.get("known_mw_images", [])
    if len(covered_known) != len(image_known):
        raise F2Error("local-coverage audit has a different known-MW image count")
    for expected, actual in zip(covered_known, image_known):
        if not isinstance(expected, dict) or not isinstance(actual, dict):
            raise F2Error("known-MW local-image records must be objects")
        expected_mask = _integer(expected["local"], "coverage known MW local image")
        actual_mask = _integer(actual["local"], "image known MW local image")
        mask_limit = (1 << coverage_dimension) - 1
        if str(expected["label"]) != str(actual["label"]) or expected_mask != actual_mask & mask_limit:
            raise F2Error("local-coverage audit does not match the known MW image manifest")

    spaces = []
    for record in coverage.get("odd_places", []):
        if record.get("classification") == "CERTIFIED_FULL_LOCAL_KUMMER_IMAGE_COVERAGE":
            spaces.append((f"p={record['rational_prime']}", [int(value) for value in record["coordinate_indices"]]))
    real = coverage.get("real_place", {})
    if real.get("classification") == "CERTIFIED_FULL_REAL_KUMMER_IMAGE_COVERAGE":
        spaces.append(("real", [int(value) for value in real["coordinate_indices"]]))
    for _, indices in spaces:
        if not indices or any(index < 0 or index >= coverage_dimension for index in indices):
            raise F2Error("local-coverage audit contains invalid coordinate indices")

    known_masks = [_integer(item["local"], "known MW local image") for item in image_known]
    survivors = []
    rejected = []
    for candidate in images.get("candidate_images", []):
        if not isinstance(candidate, dict):
            raise F2Error("candidate_images must contain objects")
        mask = _integer(candidate["local"], "candidate local image")
        obstructions = []
        for place, indices in spaces:
            basis = []
            for known_mask in known_masks:
                projected = sum(
                    ((known_mask >> index) & 1) << position
                    for position, index in enumerate(indices)
                )
                if projected:
                    basis.append(projected)
            projected_candidate = sum(
                ((mask >> index) & 1) << position
                for position, index in enumerate(indices)
            )
            if f2_reduce(projected_candidate, basis):
                obstructions.append(place)
        if obstructions:
            rejected.append({"label": str(candidate["label"]), "certified_local_obstructions": obstructions})
        else:
            survivors.append(candidate)
    return survivors, rejected, spaces


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--images", type=Path, required=True)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--local-kummer-map", type=Path)
    source.add_argument("--cover-local-audit", type=Path)
    source.add_argument("--local-coverage-audit", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    images = json.loads(args.images.read_text())
    if not isinstance(images, dict):
        raise F2Error("images input must be a JSON object")

    if args.cover_local_audit:
        cover_audit = json.loads(args.cover_local_audit.read_text())
        if not isinstance(cover_audit, dict):
            raise F2Error("cover-local audit must be a JSON object")
        survivors, obstructed, inconclusive = cover_audit_filter(images, cover_audit)
        quotient_input = dict(images)
        quotient_input["candidate_images"] = survivors
        quotient = audit_manifest(quotient_input)
        output = {
            "protocol": "BNFFREELOCAL-v1",
            "local_source": "selected_finite_two_cover_reduction_audit",
            "candidate_count": len(images.get("candidate_images", [])),
            "finite_local_survivor_count": len(survivors),
            "locally_obstructed_candidates": obstructed,
            "locally_inconclusive_candidates": inconclusive,
            "post_local_quotient": quotient,
            "status": (
                "SELECTED_FINITE_LOCAL_FILTER_ONLY: survivors have proved Q_p "
                "points at the audited finite places, but still require every "
                "other local condition and the global descent certification."
            ),
        }
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
        print(
            f"{PROTOCOL}|stage=cover_local_complete|candidates={output['candidate_count']}"
            f"|finite_local_survivors={output['finite_local_survivor_count']}"
            f"|obstructed={len(obstructed)}|inconclusive={len(inconclusive)}"
            f"|residual_signature_rank={quotient['candidate_residual_rank']}"
            "|status=SELECTED_FINITE_LOCAL_FILTER_ONLY",
            flush=True,
        )
        return

    if args.local_coverage_audit:
        coverage = json.loads(args.local_coverage_audit.read_text())
        if not isinstance(coverage, dict):
            raise F2Error("local-coverage audit must be a JSON object")
        survivors, rejected, spaces = coverage_filter(images, coverage)
        quotient_input = dict(images)
        quotient_input["candidate_images"] = survivors
        quotient = audit_manifest(quotient_input)
        output = {
            "protocol": "BNFFREELOCAL-v1",
            "local_source": "certified_odd_and_real_kummer_coverage",
            "candidate_count": len(images.get("candidate_images", [])),
            "coverage_local_survivor_count": len(survivors),
            "locally_rejected_candidates": rejected,
            "certified_covered_places": [place for place, _ in spaces],
            "post_local_quotient": quotient,
            "status": (
                "PARTIAL_CERTIFIED_LOCAL_FILTER_ONLY: rejection is certified at "
                "the listed covered places, while all uncovered places remain "
                "required local-descent work."
            ),
        }
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
        print(
            f"{PROTOCOL}|stage=coverage_complete|candidates={output['candidate_count']}"
            f"|coverage_local_survivors={output['coverage_local_survivor_count']}"
            f"|rejected={len(rejected)}|covered_places={len(spaces)}"
            f"|residual_signature_rank={quotient['candidate_residual_rank']}"
            "|status=PARTIAL_CERTIFIED_LOCAL_FILTER_ONLY",
            flush=True,
        )
        return

    local_map = json.loads(args.local_kummer_map.read_text())
    if not isinstance(local_map, dict):
        raise F2Error("local Kummer map must be a JSON object")
    local_dimension = _integer(images["local_dimension"], "local_dimension")
    basis = validate_condition_map(local_map, local_dimension)
    local_limit = 1 << local_dimension

    def checked_local_mask(record: dict, description: str) -> int:
        mask = _integer(record["local"], description)
        if mask < 0 or mask >= local_limit:
            raise F2Error(f"{description} is outside the local coordinate space")
        return mask

    for known in images.get("known_mw_images", []):
        if not isinstance(known, dict):
            raise F2Error("known_mw_images must contain objects")
        local = checked_local_mask(known, "known MW local image")
        if f2_reduce(local, basis):
            raise F2Error("a known Mordell--Weil image fails the supplied local condition")

    survivors = []
    rejected = []
    for candidate in images.get("candidate_images", []):
        if not isinstance(candidate, dict):
            raise F2Error("candidate_images must contain objects")
        local = checked_local_mask(candidate, "candidate local image")
        remainder = f2_reduce(local, basis)
        if remainder:
            rejected.append(
                {
                    "label": str(candidate["label"]),
                    "local_obstruction": remainder,
                }
            )
        else:
            survivors.append(candidate)

    quotient_input = dict(images)
    quotient_input["candidate_images"] = survivors
    quotient = audit_manifest(quotient_input)
    output = {
        "protocol": "BNFFREELOCAL-v1",
        "local_kummer_map": {
            "schema": local_map["schema"],
            "method": local_map.get("method"),
            "basis_size": len(basis),
            "basis_rank": len(basis) - sum(
                1
                for index, vector in enumerate(basis)
                if not f2_reduce(vector, basis[:index])
            ),
        },
        "candidate_count": len(images.get("candidate_images", [])),
        "local_survivor_count": len(survivors),
        "locally_rejected_candidates": rejected,
        "post_local_quotient": quotient,
        "status": (
            "LOCAL_FILTER_ONLY: surviving candidates still require the global "
            "S-class quotient bound and the descent norm condition."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        f"{PROTOCOL}|stage=complete|candidates={output['candidate_count']}"
        f"|local_survivors={output['local_survivor_count']}"
        f"|residual_signature_rank="
        f"{quotient['candidate_residual_rank']}|status=LOCAL_FILTER_ONLY",
        flush=True,
    )


if __name__ == "__main__":
    main()
