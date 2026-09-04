#!/usr/bin/env sage
"""Add the canonical rational-prime principal relations to a BNF-free ledger.

For every rational prime ``p`` represented in a complete factor base, the
principal ideal ``(p)`` supplies an exact relation among *all* prime ideals
above ``p``.  These relations are free: their generators are the rational
integers themselves, so they require neither BNF nor class-group arithmetic.
They should be present before interpreting a mod-two relation quotient.

The input ledger is never modified.  The output keeps the original relation
rows and appends one stored generator and one verified principal row per
rational prime.  Re-running on an already augmented ledger is idempotent.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

from sage.all import NumberField, PolynomialRing, QQ, ZZ, pari


PROTOCOL = "BNFFREECANON"
SCHEMA = "elliptic-curves.bnf-free-principal-relation-ledger.v1"
CANONICAL_SOURCE = "canonical_rational_prime_principal"


def rational(value: str) -> QQ:
    value = Fraction(value)
    return QQ(value.numerator) / QQ(value.denominator)


def coefficient_field(ledger: dict):
    coefficients = [rational(value) for value in ledger["defining_polynomial_ascending"]]
    if len(coefficients) != 4 or coefficients[-1] != 1:
        raise ValueError("ledger must define a monic cubic polynomial")
    if ledger.get("curve_preset") == "elkies-2026-rank28":
        from run_elkies_2026_rank28_s_class_pari import validate_inputs

        _, _, certified_coefficients, factor_hint_primes = validate_inputs()
        if coefficients != [QQ(value) for value in certified_coefficients]:
            raise ValueError("Elkies rank-28 ledger has the wrong defining cubic")
        pari.addprimes(factor_hint_primes)
    else:
        declared_primes = [
            ZZ(value) for value in ledger.get("selmer_rational_primes", ())
        ]
        if any(value < 2 or not value.is_prime() for value in declared_primes):
            raise ValueError("ledger selmer_rational_primes are not proved prime")
        pari.addprimes(declared_primes)
    ring = PolynomialRing(QQ, "x")
    x = ring.gen()
    return NumberField(sum(value * x**index for index, value in enumerate(coefficients)), "theta")


def prime_signature(prime) -> tuple[str, int, int, int]:
    return (
        str(prime.pari_hnf()),
        int(prime.norm()),
        int(prime.residue_class_degree()),
        int(prime.smallest_integer()),
    )


def canonical_row(field, factor_base_column: dict, rational_prime: int) -> int:
    """Return the exact valuation-parity row of the principal ideal ``(p)``."""

    row = 0
    for prime, exponent in field.ideal(ZZ(rational_prime)).factor():
        if int(exponent) & 1:
            column = factor_base_column.get(prime_signature(prime))
            if column is None:
                raise ValueError(
                    f"factor base omits a prime ideal above canonical prime {rational_prime}"
                )
            row |= 1 << column
    return row


def existing_canonical_rows(ledger: dict) -> dict[int, dict]:
    records = {}
    for relation in ledger.get("closed_relations", ()):
        if not isinstance(relation, dict):
            raise ValueError("closed relation is not an object")
        if relation.get("source") != CANONICAL_SOURCE:
            continue
        rational_prime = int(relation["rational_prime"])
        if rational_prime in records:
            raise ValueError(f"duplicate canonical relation for rational prime {rational_prime}")
        records[rational_prime] = relation
    return records


def augment(ledger: dict) -> tuple[dict, int, int]:
    if ledger.get("schema") != SCHEMA:
        raise ValueError("unexpected relation-ledger schema")
    if not isinstance(ledger.get("factor_base"), list):
        raise ValueError("ledger factor_base must be a list")
    if not isinstance(ledger.get("generators"), list) or not isinstance(
        ledger.get("closed_relations"), list
    ):
        raise ValueError("ledger generators and closed_relations must be lists")

    field = coefficient_field(ledger)
    factor_base_column = {}
    rational_primes = set()
    for index, record in enumerate(ledger["factor_base"]):
        if not isinstance(record, dict):
            raise ValueError("factor-base entry is not an object")
        signature = (
            str(record["hnf"]),
            int(record["norm"]),
            int(record["residue_degree"]),
            int(record["rational_prime"]),
        )
        if signature in factor_base_column:
            raise ValueError("factor base has duplicate prime ideals")
        factor_base_column[signature] = index
        rational_primes.add(signature[3])

    existing = existing_canonical_rows(ledger)
    added = 0
    for rational_prime in sorted(rational_primes):
        row = canonical_row(field, factor_base_column, rational_prime)
        present = existing.get(rational_prime)
        if present is not None:
            if int(present.get("fb_parity_mask_hex", "0"), 16) != row:
                raise ValueError(
                    f"stored canonical relation for {rational_prime} has the wrong parity row"
                )
            continue
        generator_index = len(ledger["generators"])
        ledger["generators"].append(
            {
                "power_basis": [str(rational_prime), "0", "0"],
                "source": CANONICAL_SOURCE,
                "rational_prime": rational_prime,
                "norm": str(ZZ(rational_prime) ** field.degree()),
            }
        )
        ledger["closed_relations"].append(
            {
                "fb_parity_mask_hex": hex(row),
                "generator_indices": [generator_index],
                "kind": CANONICAL_SOURCE,
                "source": CANONICAL_SOURCE,
                "rational_prime": rational_prime,
            }
        )
        added += 1

    ledger["canonical_principal_relations"] = {
        "method": "exact factorization of the rational principal ideals (p)",
        "source": CANONICAL_SOURCE,
        "rational_primes": sorted(rational_primes),
        "completed_relation_count": len(rational_primes),
        "added_in_this_augmentation": added,
        "principal_generators_stored": True,
    }
    return ledger, added, len(rational_primes)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--relation-ledger", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.relation_ledger.resolve() == args.output.resolve():
        raise ValueError("write the augmented ledger to a distinct path")
    ledger = json.loads(args.relation_ledger.read_text())
    if not isinstance(ledger, dict):
        raise ValueError("relation ledger must be a JSON object")
    augmented, added, total = augment(ledger)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(augmented, indent=2, sort_keys=True) + "\n")
    print(
        f"{PROTOCOL}|stage=complete|rational_prime_relations={total}"
        f"|added={added}|closed_relations={len(augmented['closed_relations'])}"
        "|status=EXACT_PRINCIPAL_ROWS",
        flush=True,
    )


if __name__ == "__main__":
    main()
