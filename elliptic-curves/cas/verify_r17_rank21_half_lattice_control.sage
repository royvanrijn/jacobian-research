#!/usr/bin/env sage -python
"""Verification-only quotient audit for the blind R17 rank-21 replay."""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
BLIND = ROOT / "artifacts/generated-results/elliptic-curves/half_lattice_r17_rank21_blind_v1.json"
FIXTURE = ROOT / "artifacts/generated-results/elliptic-curves/icarm_curve394_rank21_v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/half_lattice_r17_rank21_verification_v1.json"
sys.path[:0] = [str(ROOT / "elliptic-curves"), str(CAS)]

from ecsearch.q12o5867_specialization import evaluate_projective_specialization, global_minimal_model_with_change, load_q12o5867_data, short_certificate_model
from elliptic_candidate_record import source_point_to_target
from mod2_reduction_independence import find_mod2_reduction_certificate
from search_nagao_u135_alternate_covers import relation_proposals


MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"


def binary_rank(values):
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


def restricted_rank(signatures, columns):
    return binary_rank(
        sum((int(row[column]) & 1) << offset for offset, column in enumerate(columns))
        for signature in signatures for row in signature.rows
    )


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def main():
    blind_bytes = BLIND.read_bytes()
    blind_hash = sha256(blind_bytes).hexdigest()
    blind = json.loads(blind_bytes)
    if blind.get("status") != "PASS_BOUNDED_BLIND_R17_RANK21_CONTROL":
        raise ValueError("blind rank21 replay did not finish")
    if blind["blindness_boundary"]["loaded_public_rank21_points"] is not False:
        raise ValueError("blind rank21 boundary failed")
    fixture = json.loads(FIXTURE.read_text())

    family = load_q12o5867_data(MODEL, SECTIONS)
    specialization = evaluate_projective_specialization(family, 3, 8)
    minimal_model, change1, unused = global_minimal_model_with_change(specialization.model)
    short_model, change2 = short_certificate_model(minimal_model)
    if [str(value) for value in short_model] != blind["short_model"]:
        raise ArithmeticError("blind and public rank21 models differ")
    generic = tuple(
        source_point_to_target(source_point_to_target(point, change1), change2)
        for point in specialization.points
    )
    public = tuple(
        source_point_to_target((Fraction(point[0]), Fraction(point[1])), change2)
        for point in fixture["public_point_replay"]["points"]
    )
    signatures = find_mod2_reduction_certificate(short_model, generic + public, prime_bound=1800)
    columns = list(range(17))
    rank = restricted_rank(signatures, columns)
    if rank != 17:
        raise ArithmeticError("rank21 verification lost generic independence")
    complement_indices = []
    for index in range(len(public)):
        trial = columns + [17 + index]
        trial_rank = restricted_rank(signatures, trial)
        if trial_rank > rank:
            complement_indices.append(index)
            columns, rank = trial, trial_rank
    if rank != 21 or len(complement_indices) != 4:
        raise ArithmeticError("rank21 public complement did not have dimension four")
    basis = generic + tuple(public[index] for index in complement_indices)
    candidate_rows = blind["blind_result"]["candidate_points"]
    points = tuple(
        (Fraction(row["point"]["x"]), Fraction(row["point"]["y"])) for row in candidate_rows
    )
    relations = relation_proposals(short_model, basis, points, timeout=180.0, stack_bytes=1_000_000_000)
    if not all(exact for unused_relation, exact in relations):
        raise ArithmeticError("a blind rank21 candidate missed the exact public subgroup")
    masks = []
    by_center = defaultdict(set)
    exact_rows = []
    for source, (relation, exact) in zip(candidate_rows, relations):
        quotient = relation[17:]
        mask = sum((int(value) & 1) << offset for offset, value in enumerate(quotient))
        masks.append(mask)
        for center in source["source_masks"]:
            by_center[int(center)].add(mask)
        exact_rows.append(
            {
                "source_half_classes": [f"0x{int(center):05x}" for center in source["source_masks"]],
                "relation": list(relation),
                "quotient_coordinates": list(quotient),
                "quotient_hex": f"0x{mask:x}",
            }
        )
    recovered = binary_rank(masks)
    payload = {
        "schema": "elliptic-curves.half-lattice-r17-rank21-verification.v1",
        "status": "PASS_EXACT_R17_RANK21_HELDOUT_RELATIONS",
        "blind_artifact_sha256_before_fixture_load": blind_hash,
        "public_exceptional_quotient_dimension": 4,
        "exact_blind_recovered_quotient_dimension": recovered,
        "full_public_quotient_recovered": recovered == 4,
        "public_complement_indices_one_based": [index + 1 for index in complement_indices],
        "productive_centers": [
            {"half_class": f"0x{center:05x}", "quotient_span": binary_rank(values)}
            for center, values in sorted(by_center.items())
        ],
        "verified_candidate_sample_count": len(points),
        "relations": exact_rows,
        "claim_boundary": ["The quotient dimension, relations, and recovered span are exact inside the displayed public rank-21 subgroup.", "The underlying point search remains bounded."],
        "input_hashes": {
            str(BLIND.relative_to(ROOT)): blind_hash,
            str(FIXTURE.relative_to(ROOT)): digest(FIXTURE),
            str(Path(__file__).resolve().relative_to(ROOT)): digest(Path(__file__).resolve()),
        },
        "reproducing_command": "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python elliptic-curves/cas/verify_r17_rank21_half_lattice_control.sage",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"R17RANK21VERIFY|status=PASS|recovered={recovered}/4|output={OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
