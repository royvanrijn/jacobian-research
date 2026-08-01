#!/usr/bin/env python3
"""Verify the exact sign involution between the two Case-1 branches.

The archived systems use coefficient dictionaries over the pinned degree-35
first-block field.  The branch involution fixes h and u3, negates u1 and u2,
and rescales the displayed equations by recorded field units ±1.  Because the
field coefficients themselves are unchanged, this transport can be checked
exactly on the serialized rational data without a CAS dependency.
"""
from __future__ import annotations

import pickle
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / (
    "plane-jc/external/zenodo-21479814/"
    "bilLkarkariy-jc2-72-108-exact-certificates-d9ea4fd/"
    "release_bundle/exact_replay"
)


def negate_coefficient(coefficient: dict[int, tuple[int, int]]):
    return {exponent: (-numerator, denominator) for exponent, (numerator, denominator) in coefficient.items()}


def transport_polynomial(polynomial, row_scale: int):
    transported = {}
    for monomial, coefficient in polynomial.items():
        sign = row_scale * (-1) ** (monomial[1] + monomial[2])
        transported[monomial] = coefficient if sign == 1 else negate_coefficient(coefficient)
    return transported


def check_pair(left_name: str, right_name: str, row_scales: tuple[int, ...]) -> None:
    left = pickle.loads((ROOT / left_name).read_bytes())
    right = pickle.loads((ROOT / right_name).read_bytes())
    assert len(left) == len(right) == len(row_scales)
    for index, (source, target, row_scale) in enumerate(zip(left, right, row_scales), start=1):
        assert transport_polynomial(source, row_scale) == target, (
            f"branch transport mismatch in {left_name}, row {index}"
        )


check_pair(
    "case1_branch1_after_w.pkl",
    "case1_branch2_after_w.pkl",
    (1, 1, -1, -1, -1, -1, -1),
)
check_pair(
    "hne0_polred.pkl",
    "hne0_branch2_polred.pkl",
    (1, 1, -1, -1, -1, -1),
)

print("CASE1_FULL_SYSTEM_SYMMETRY_PASS")
print("CASE1_HARD_RESIDUAL_SYMMETRY_PASS")
print("CASE1_BRANCH_SYMMETRY_PASS")
