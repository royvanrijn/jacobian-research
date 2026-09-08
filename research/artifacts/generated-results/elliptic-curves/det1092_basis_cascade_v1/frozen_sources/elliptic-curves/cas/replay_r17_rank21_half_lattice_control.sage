#!/usr/bin/env sage -python
"""Blind fixed-generic-deepest43 replay at the R17 rank-21 control t=3/8."""

from __future__ import annotations

from decimal import Decimal
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
ENGINE = CAS / "half_lattice_fake_descent_replay.sage"
OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/half_lattice_r17_rank21_blind_v1.json"
sys.path[:0] = [str(ROOT / "elliptic-curves"), str(CAS)]

from ecsearch.q12o5867_specialization import evaluate_projective_specialization, global_minimal_model_with_change, load_q12o5867_data, short_certificate_model
from elliptic_candidate_record import source_point_to_target
from mod2_reduction_independence import combined_mod2_rank, find_mod2_reduction_certificate
from search_nagao_u135_alternate_covers import relation_proposals


MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


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


def main():
    engine = SourceFileLoader("half_lattice_rank21_engine", str(ENGINE)).load_module()
    family = load_q12o5867_data(MODEL, SECTIONS)
    specialization = evaluate_projective_specialization(family, 3, 8)
    minimal_model, change1, unused = global_minimal_model_with_change(specialization.model)
    short_model, change2 = short_certificate_model(minimal_model)
    generic = tuple(
        source_point_to_target(source_point_to_target(point, change1), change2)
        for point in specialization.points
    )
    signatures0 = find_mod2_reduction_certificate(short_model, generic, prime_bound=800)
    if combined_mod2_rank(signatures0, 17) != 17:
        raise ArithmeticError("rank21 control generic specialization lost independence")

    generic_oracle = engine.CosetOracle(engine.GENERIC_GRAM)
    deepest = []
    for mask in range(1 << 17):
        norm, representative, unused_error = generic_oracle.solve(mask)
        if norm == 12:
            deepest.append(mask)
    if len(deepest) != 43:
        raise ArithmeticError("generic deepest count changed")
    height_gram = engine.canonical_height_gram(short_model, generic)
    rounded = tuple(
        tuple(int((value * Decimal(1_000_000)).to_integral_value()) for value in row)
        for row in height_gram
    )
    oracle = engine.CosetOracle(rounded)
    discoveries = {}
    covers = []
    for position, mask in enumerate(deepest, 1):
        unused_norm, representative, unused_error = oracle.solve(mask)
        depth = engine.quadratic_decimal(height_gram, representative) / 4
        outcome = engine.run_quartic_search(
            mask=mask, representative=representative, short_model=short_model,
            generic_points=generic, height_bound=100_000, timeout_seconds=15.0,
            stack_bytes=1_000_000_000,
        )
        for point in outcome.curve_points:
            discoveries.setdefault(point, set()).add(mask)
        record = outcome.record
        covers.append(
            {
                "mask": mask,
                "hex": f"0x{mask:05x}",
                "specialized_depth": str(depth),
                "status": record["status"],
                "finite_curve_point_count": len(outcome.curve_points),
                "reduced_coefficient_bits": record.get("reduced_model", {}).get("maximum_coefficient_bits"),
                "modular_density_product": record.get("local_stage", {}).get("joint_independent_density_product"),
                "search_milliseconds": record.get("search_milliseconds"),
            }
        )
        print(f"R17RANK21HALF|cover={position}/43|mask={mask:#x}|points={len(outcome.curve_points)}", flush=True)

    basis_signs = {signed for point in generic for signed in (point, (point[0], -point[1]))}
    candidates = tuple(sorted((point for point in discoveries if point not in basis_signs), key=str))
    proposals = relation_proposals(short_model, generic, candidates, timeout=120.0, stack_bytes=1_000_000_000)
    unexplained = tuple(point for point, (unused_relation, exact) in zip(candidates, proposals) if not exact)
    # Thousands of affine hits collectively exclude every small reduction
    # prime through their denominators.  Keep a deterministic verification
    # sample instead of emitting a bogus negative finite-reduction gain.
    verification_sample = unexplained[:256]
    candidate_rows = [
        {
            "point": {"x": str(point[0]), "y": str(point[1])},
            "source_masks": sorted(discoveries[point]),
        }
        for point in verification_sample
    ]
    payload = {
        "schema": "elliptic-curves.half-lattice-r17-rank21-blind.v1",
        "status": "PASS_BOUNDED_BLIND_R17_RANK21_CONTROL",
        "parameter": "3/8",
        "blindness_boundary": {"loaded_public_rank21_points": False},
        "declared_budget": {"generic_deepest_classes": 43, "height_bound_each": 100000, "timeout_seconds_each": 15.0},
        "short_model": [str(value) for value in short_model],
        "generic_points": [{"x": str(point[0]), "y": str(point[1])} for point in generic],
        "cover_records": covers,
        "blind_result": {
            "distinct_nonbasis_candidates": len(candidates),
            "unexplained_candidate_count": len(unexplained),
            "verification_sample_count": len(verification_sample),
            "finite_mod2_quotient_gain": None,
            "finite_reduction_certificate_valid": False,
            "candidate_points": candidate_rows,
        },
        "claim_boundary": ["Every returned point satisfies the curve equation exactly.", "The numerical generic-relation failure is not an independence proof; the separate fixture verifier is authoritative.", "Search misses are bounded; the rank-21 public fixture was not loaded."],
        "input_hashes": {
            str(MODEL.relative_to(ROOT)): digest(MODEL),
            str(SECTIONS.relative_to(ROOT)): digest(SECTIONS),
            str(ENGINE.relative_to(ROOT)): digest(ENGINE),
            str(Path(__file__).resolve().relative_to(ROOT)): digest(Path(__file__).resolve()),
        },
        "reproducing_command": "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python elliptic-curves/cas/replay_r17_rank21_half_lattice_control.sage",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"R17RANK21HALF|status=PASS|gain=verification_required|output={OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
