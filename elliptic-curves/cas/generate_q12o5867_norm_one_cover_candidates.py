#!/usr/bin/env sage
"""Generate bounded exact norm-one cubic classes for q12o5867 covers.

For every primitive projective power-basis direction
``gamma=a+b*theta+c*theta^2`` in the declared box, form

    alpha = gamma^3 / Norm(gamma).

Then ``Norm(alpha)=1`` exactly, so ``alpha`` passes the norm condition for the
standard odd-degree 2-descent.  The corresponding intersections of quadrics
can be built immediately by ``build_bnf_free_two_covers.py`` and subjected to
strict local-solubility tests.  Most bounded classes are expected to fail a
local condition; generation alone is neither Selmer membership nor evidence
for a rational point.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import product
import json
from math import gcd
from pathlib import Path
import sys


PROTOCOL = "Q12NORMONE"
SIGNATURE_SCHEMA = "elliptic-curves.bnf-free-signature-map.v1"
OUTPUT_SCHEMA = "elliptic-curves.bnf-free-norm-filtered-squareclass-candidates.v1"


def normalize_direction(values: tuple[int, int, int]) -> tuple[int, int, int] | None:
    if values == (0, 0, 0):
        return None
    common = gcd(gcd(abs(values[0]), abs(values[1])), abs(values[2]))
    values = tuple(value // common for value in values)
    first = next(value for value in values if value)
    if first < 0:
        values = tuple(-value for value in values)
    return values


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--signature-map", type=Path, required=True)
    parser.add_argument("--coefficient-bound", type=int, default=3)
    parser.add_argument("--max-candidates", type=int, default=500)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.coefficient_bound < 1 or args.max_candidates < 1:
        parser.error("coefficient and candidate bounds must be positive")

    from sage.all import NumberField, PolynomialRing, QQ

    signature = json.loads(args.signature_map.read_text())
    if signature.get("schema") != SIGNATURE_SCHEMA:
        raise ValueError("expected a BNF-free signature map")
    coefficients = [QQ(Fraction(value)) for value in signature["defining_polynomial_ascending"]]
    if len(coefficients) != 4 or coefficients[-1] != 1:
        raise ValueError("signature map must define a monic cubic")
    ring = PolynomialRing(QQ, "x")
    x = ring.gen()
    field = NumberField(
        sum(coefficient * x**index for index, coefficient in enumerate(coefficients)),
        "theta",
    )
    theta = field.gen()

    directions = {
        normalized
        for values in product(
            range(-args.coefficient_bound, args.coefficient_bound + 1), repeat=3
        )
        if (normalized := normalize_direction(values)) is not None
    }
    candidates = []
    seen_alpha = set()
    for a, b, c in sorted(directions):
        gamma = field(a + b * theta + c * theta**2)
        norm = QQ(gamma.norm())
        if not norm:
            continue
        alpha = gamma**3 / norm
        if QQ(alpha.norm()) != 1:
            raise ArithmeticError("norm-one normalization failed")
        alpha_coefficients = list(alpha.polynomial().list())
        alpha_coefficients += [QQ(0)] * (3 - len(alpha_coefficients))
        key = tuple(alpha_coefficients[:3])
        if key in seen_alpha:
            continue
        seen_alpha.add(key)
        candidates.append(
            {
                "label": f"norm-one-gamma-{a}_{b}_{c}",
                "generator_coefficients": [str(value) for value in key],
                "norm": "1",
                "source_gamma_coefficients": [a, b, c],
                "normalization": "alpha=gamma^3/Norm(gamma)",
            }
        )
        if len(candidates) >= args.max_candidates:
            break

    output = {
        "schema": OUTPUT_SCHEMA,
        "status": "exact_norm_one_candidates_not_local_selmer_certificate",
        "field_polynomial_ascending": [str(value) for value in coefficients],
        "source_signature_map": str(args.signature_map.resolve()),
        "coefficient_bound": args.coefficient_bound,
        "primitive_projective_direction_count": len(directions),
        "candidate_truncated": len(candidates) < len(directions),
        "candidates": candidates,
        "claim_boundary": [
            "Every displayed alpha has exact field norm one.",
            "Norm one is necessary but not sufficient for a 2-Selmer class.",
            "Local solubility and a rational point on the cover remain separate gates.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        f"{PROTOCOL}|stage=complete|directions={len(directions)}"
        f"|candidates={len(candidates)}|status=EXACT_NORM_ONE_NOT_SELMER",
        flush=True,
    )


if __name__ == "__main__":
    main()
