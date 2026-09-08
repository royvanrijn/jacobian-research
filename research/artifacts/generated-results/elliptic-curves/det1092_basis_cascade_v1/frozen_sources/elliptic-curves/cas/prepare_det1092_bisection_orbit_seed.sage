#!/usr/bin/env sage-python
"""Specialize only the redacted determinant-1092 generic MW17 parent."""

from __future__ import annotations

import argparse
import sys
from fractions import Fraction as F
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ


CAS = Path(__file__).resolve().parent
sys.path.insert(0, str(CAS))
import certify_compact_r17_candidates as cert
import det1092_bisection_orbit_holdout as campaign
from research_runtime.memory_store import MemoryFactStore
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as ReductionCache
from research_runtime.search_state import raw_state
from research_runtime.store import checkpoint


def main(index: int) -> None:
    campaign.install_execution_guard()
    row, directory = campaign.configure(index)
    output = directory / "seed.json"
    if output.exists():
        raise FileExistsError("preserve generic-only specialized seed")
    parent = cert.read(campaign.LOCAL / "parent-sections.json")
    t = QQ(row["parameter"])
    ring = PolynomialRing(QQ, "t")
    evaluate = lambda item: ring(item["numerator"])(t) / ring(item["denominator"])(t)
    a = [evaluate(item) for item in parent["a_invariants"]]
    if a[:3] != [0, 0, 0]:
        raise ArithmeticError("redacted generic determinant-reduced parent is not short")
    q = t.denominator()
    model = [0, 0, 0, a[3] * q**8, a[4] * q**12]
    if list(map(str, model)) != row["model"]:
        raise ArithmeticError("frozen held-out model is not the exact homogeneous specialization")
    curve = EllipticCurve(QQ, model)
    points = [curve([evaluate(point[0]) * q**4, evaluate(point[1]) * q**6])
              for point in parent["basis_weierstrass_coordinates"]]
    if len(points) != campaign.INITIAL_RANK or any(point.is_zero() for point in points):
        raise ArithmeticError("generic MW17 specialization failed")
    fractions = tuple(tuple(F(str(value)) for value in point.xy()) for point in points)
    rational_model = tuple(F(str(value)) for value in model)
    state = raw_state(rational_model, fractions, cache=ReductionCache(MemoryFactStore()), prime_bound=1000)
    if (tuple(tuple(F(value) for value in point) for point in state.basis) != fractions
            or state.rank != campaign.INITIAL_RANK):
        raise ArithmeticError("specialized generic points do not reproduce the exact seed")
    rank_proof = campaign.checked_rank(rational_model, fractions, state.reductions.primes, state.no_two_torsion_prime)
    checkpoint(output, {
        "schema": "elliptic-curves.det1092-bisection-orbit-seed.v1", "family": row["family"], "parameter": row["parameter"],
        "curve": list(map(str, model)), "points": [list(map(str, point)) for point in fractions],
        "generic_points": [list(map(str, point)) for point in fractions], "rank_certificate": rank_proof,
        "parent_redaction_sha256": campaign.sha(campaign.LOCAL / "parent-sections.json"),
        "candidate_row": row, "homogeneous_scale": str(q**2),
        "coordinate_transport": "certified generic determinant-reduced short chart, then (x,y) -> (q^4 x,q^6 y)",
    })
    print("BISECTION GENERIC MW17 SEED", row["id"], flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", required=True, type=int)
    main(parser.parse_args().index)
