#!/usr/bin/env sage
"""Turn exact closed relation-generator products into K(S,2) candidates.

The input is a per-run principal-relation ledger from one of the BNF-free
collectors.  A *combination* of relations whose factor-base parity is
supported entirely at S has even valuation outside S, hence its product of
stored principal generators is an explicit element of K(S,2).  We find a
GF(2) basis of such combinations by eliminating the non-S columns.  This is
strictly stronger than retaining only individual S-supported rows.  No class
group, BNF, or square-root ideal calculation is used here.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

from sage.all import NumberField, PolynomialRing, QQ


PROTOCOL = "BNFFREESQ"


def rational(value: str) -> QQ:
    fraction = Fraction(value)
    return QQ(fraction.numerator) / QQ(fraction.denominator)


def s_supported_dependencies(rows: list[int], s_mask: int) -> list[int]:
    """Return an independent basis of combinations supported outside no S-column.

    Each returned packed mask indexes ``rows``.  Sparse elimination is on the
    projection to the complement of S, while the provenance mask retains the
    exact principal-generator products needed for a reproducible candidate.
    """

    pivots: dict[int, tuple[int, int]] = {}
    dependencies = []
    for index, full_row in enumerate(rows):
        row = full_row & ~s_mask
        provenance = 1 << index
        while row:
            pivot = row.bit_length() - 1
            previous = pivots.get(pivot)
            if previous is None:
                pivots[pivot] = (row, provenance)
                break
            row ^= previous[0]
            provenance ^= previous[1]
        if row == 0:
            dependencies.append(provenance)
    return dependencies


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--relation-ledger", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-candidates", type=int, default=100)
    parser.add_argument(
        "--individual-s-supported-only",
        action="store_true",
        help="legacy diagnostic: do not eliminate non-S columns across rows",
    )
    args = parser.parse_args()
    if args.max_candidates < 1:
        raise ValueError("--max-candidates must be positive")

    ledger = json.loads(args.relation_ledger.read_text())
    if ledger.get("schema") != "elliptic-curves.bnf-free-principal-relation-ledger.v1":
        raise ValueError("unexpected relation-ledger schema")
    coefficients = [rational(value) for value in ledger["defining_polynomial_ascending"]]
    if len(coefficients) != 4 or coefficients[-1] != 1:
        raise ValueError("relation ledger needs a monic cubic polynomial")
    ring = PolynomialRing(QQ, "x")
    x = ring.gen()
    polynomial = sum(coefficients[index] * x**index for index in range(4))
    field = NumberField(polynomial, "theta")
    theta = field.gen()

    s_columns = {int(index) for index in ledger.get("S_columns", ())}
    s_mask = sum(1 << index for index in s_columns)
    relations = list(ledger.get("closed_relations", ()))
    relation_masks = [int(relation["fb_parity_mask_hex"], 16) for relation in relations]
    if args.individual_s_supported_only:
        dependencies = [
            1 << relation_index
            for relation_index, mask in enumerate(relation_masks)
            if mask & ~s_mask == 0
        ]
        generation_mode = "individual_S_supported_rows_only"
    else:
        dependencies = s_supported_dependencies(relation_masks, s_mask)
        generation_mode = "non_S_projection_kernel"
    candidates = []
    skipped_non_s = 0
    for mask in relation_masks:
        if mask & ~s_mask:
            skipped_non_s += 1
    for dependency_index, dependency in enumerate(dependencies):
        relation_indices = [
            relation_index
            for relation_index in range(len(relations))
            if dependency >> relation_index & 1
        ]
        combined_mask = 0
        product = field(1)
        generator_indices = []
        for relation_index in relation_indices:
            relation = relations[relation_index]
            combined_mask ^= relation_masks[relation_index]
            for generator_index in relation["generator_indices"]:
                generator = ledger["generators"][generator_index]
                coordinates = [rational(value) for value in generator["power_basis"]]
                product *= sum(
                    coordinates[index] * theta**index for index in range(3)
                )
                generator_indices.append(generator_index)
        if combined_mask & ~s_mask:
            raise ArithmeticError("non-S elimination produced an invalid dependency")
        support = [
            index for index in range(combined_mask.bit_length())
            if combined_mask >> index & 1
        ]
        coefficients_out = list(product.polynomial().list())
        coefficients_out += [QQ(0)] * (3 - len(coefficients_out))
        candidate = {
            "label": f"s-supported-kernel-{dependency_index}",
            "generator_coefficients": [str(value) for value in coefficients_out[:3]],
            "source_relation_indices": relation_indices,
            "source_kind": (
                relations[relation_indices[0]]["kind"]
                if len(relation_indices) == 1
                else "s-supported-relation-kernel"
            ),
            "factor_base_support": support,
            "generator_indices": generator_indices,
        }
        if len(relation_indices) == 1:
            # Retain the original, convenient provenance field for downstream
            # users that consume a single closed relation at a time.
            candidate["source_relation"] = relation_indices[0]
        candidates.append(candidate)
        if len(candidates) >= args.max_candidates:
            break

    output = {
        "schema": "elliptic-curves.bnf-free-squareclass-candidates.v1",
        "status": "explicit_KS2_candidates_not_local_selmer_certificate",
        "field_polynomial_ascending": [str(value) for value in coefficients],
        "candidates": candidates,
        "skipped_non_S_supported_relations": skipped_non_s,
        "s_supported_kernel_dimension": len(dependencies),
        "candidate_generation": generation_mode,
        "candidate_truncated": len(candidates) < len(dependencies),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        f"{PROTOCOL}|stage=complete|candidates={len(candidates)}"
        f"|s_supported_kernel_dimension={len(dependencies)}"
        f"|skipped_non_S_supported_relations={skipped_non_s}"
        f"|status=EXPLICIT_KS2_CANDIDATES_NOT_SELMER",
        flush=True,
    )


if __name__ == "__main__":
    main()
