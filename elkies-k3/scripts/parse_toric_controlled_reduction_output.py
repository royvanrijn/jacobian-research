#!/usr/bin/env python3
"""Strict parser for one ToricControlledReduction ``readfile`` output line.

This module deliberately knows nothing about the R17 surface or the expected
Frobenius polynomial.  The mathematical verifier imports it and separately
checks every parsed field against the exact model.
"""

from __future__ import annotations

import json


def parse_readfile_output(raw_text: str) -> dict:
    """Parse the eight colon-delimited fields emitted by ``readfile.exe``."""

    raw_line = raw_text.strip()
    if not raw_line or "\n" in raw_line or "\r" in raw_line:
        raise ValueError("expected exactly one nonempty output line")
    parts = raw_line.split(":")
    if len(parts) != 8:
        raise ValueError("unexpected ToricControlledReduction output format")
    label = parts[0]
    if not label:
        raise ValueError("empty ToricControlledReduction label")
    monomials = json.loads(parts[1])
    coefficients = json.loads(parts[2])
    halfspace_a = json.loads(parts[3])
    halfspace_b = json.loads(parts[4])
    prime = int(parts[5])
    hodge_numbers = json.loads(parts[6])
    frobenius_coefficients = json.loads(parts[7])
    if not all(isinstance(row, list) and len(row) == 3 for row in monomials):
        raise ValueError("monomial list is not a three-column integer matrix")
    if not all(isinstance(row, list) and len(row) == 3 for row in halfspace_a):
        raise ValueError("half-space matrix is not three-column")
    if len(monomials) != len(coefficients):
        raise ValueError("monomial and coefficient lengths differ")
    if len(halfspace_a) != len(halfspace_b):
        raise ValueError("half-space matrix and vector lengths differ")
    for values in (
        monomials,
        coefficients,
        halfspace_a,
        halfspace_b,
        hodge_numbers,
        frobenius_coefficients,
    ):
        flattened = (
            [entry for row in values for entry in row]
            if values and isinstance(values[0], list)
            else values
        )
        if not all(isinstance(value, int) and not isinstance(value, bool) for value in flattened):
            raise ValueError("noninteger datum in ToricControlledReduction output")
    return {
        "label": label,
        "monomials": [tuple(row) for row in monomials],
        "coefficients": coefficients,
        "halfspace_A": [tuple(row) for row in halfspace_a],
        "halfspace_b": halfspace_b,
        "prime": prime,
        "hodge_numbers": hodge_numbers,
        "frobenius_coefficients": frobenius_coefficients,
    }

