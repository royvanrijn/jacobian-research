#!/usr/bin/env sage
"""Select repeated residual ideals for an adaptive Minkowski special-q pass.

The useful vertices in a large-prime forest are not fresh special ideals but
outside ideals already hit by several exact partial relations.  This tool
reads a mergeable relation ledger, ranks those vertices by exact incidence,
reconstructs their degree-one residues, and emits a replayable seed list for
``run_fermigier_rank20_minkowski_specialq.py``.

Selection is graph scheduling only.  It proves no new principal relation,
class-group bound, local condition, or Selmer statement.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from sage.all import GF, NumberField, PolynomialRing, QQ, ZZ, pari


PROTOCOL = "BNFFREEFEEDBACK"
SCHEMA = "elliptic-curves.bnf-free-principal-relation-ledger.v1"


def rational(value: str):
    parsed = Fraction(value)
    return QQ(parsed.numerator) / QQ(parsed.denominator)


def field_from_ledger(ledger: dict):
    coefficients = [
        rational(value) for value in ledger["defining_polynomial_ascending"]
    ]
    if len(coefficients) != 4 or coefficients[-1] != 1:
        raise ValueError("ledger must define a monic cubic polynomial")
    declared_primes = [
        ZZ(value) for value in ledger.get("selmer_rational_primes", ())
    ]
    if any(value < 2 or not value.is_prime() for value in declared_primes):
        raise ValueError("ledger selmer_rational_primes are not proved prime")
    pari.addprimes(declared_primes)
    ring = PolynomialRing(QQ, "x")
    x = ring.gen()
    polynomial = sum(
        coefficient * x**index for index, coefficient in enumerate(coefficients)
    )
    field = NumberField(polynomial, "theta")
    reported = ledger.get("field_discriminant")
    if reported is not None and ZZ(reported) != field.discriminant():
        raise ValueError("ledger field discriminant does not match its polynomial")
    return field, polynomial


def select_feedback_vertices(
    ledger: dict,
    *,
    minimum_occurrences: int,
    maximum_rational_prime: int,
    maximum_seeds: int,
) -> list[dict]:
    if ledger.get("schema") != SCHEMA:
        raise ValueError("unexpected relation-ledger schema")
    if not isinstance(ledger.get("partial_relations"), list):
        raise ValueError("ledger does not retain exact partial relations")
    if not isinstance(ledger.get("generators"), list):
        raise ValueError("ledger does not retain exact generators")

    source_rational_primes = {
        int(rational_prime)
        for generator in ledger["generators"]
        for rational_prime, _ in generator.get("source_special_ideals", ())
    }
    incidences = Counter()
    for relation_index, relation in enumerate(ledger["partial_relations"]):
        if not isinstance(relation, dict):
            raise ValueError(f"partial relation {relation_index} is not an object")
        for vertex in relation.get("large_prime_vertices", ()):
            if not isinstance(vertex, list) or len(vertex) != 2:
                raise ValueError(
                    f"partial relation {relation_index} has an invalid vertex"
                )
            rational_prime, ideal_hnf = int(vertex[0]), str(vertex[1])
            if rational_prime not in source_rational_primes:
                incidences[(rational_prime, ideal_hnf)] += 1

    field, polynomial = field_from_ledger(ledger)
    theta = field.gen()
    candidates = []
    for (rational_prime, ideal_hnf), occurrence_count in sorted(
        incidences.items(), key=lambda item: (-item[1], item[0][0], item[0][1])
    ):
        if occurrence_count < minimum_occurrences:
            continue
        if rational_prime > maximum_rational_prime:
            continue
        prime = next(
            (
                candidate
                for candidate in field.primes_above(ZZ(rational_prime))
                if str(candidate.pari_hnf()) == ideal_hnf
            ),
            None,
        )
        if prime is None:
            raise ValueError("a repeated vertex cannot be reconstructed")
        if int(prime.residue_class_degree()) != 1:
            continue
        residues = [
            ZZ(value)
            for value in polynomial.change_ring(GF(rational_prime)).roots(
                multiplicities=False
            )
            if int((theta - ZZ(value)).valuation(prime)) > 0
        ]
        if len(residues) != 1:
            raise ArithmeticError("degree-one prime did not have one matching residue")
        candidates.append(
            {
                "rational_prime": rational_prime,
                "residue": int(residues[0]),
                "prime_ideal_hnf": ideal_hnf,
                "partial_relation_incidence_count": occurrence_count,
            }
        )
        if len(candidates) >= maximum_seeds:
            break
    return candidates


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--relation-ledger", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--minimum-occurrences", type=int, default=2)
    parser.add_argument("--maximum-rational-prime", type=int, default=1 << 40)
    parser.add_argument("--maximum-seeds", type=int, default=20)
    args = parser.parse_args()
    if args.minimum_occurrences < 2:
        raise ValueError("--minimum-occurrences must be at least two")
    if args.maximum_rational_prime < 2 or args.maximum_seeds < 1:
        raise ValueError("feedback bounds must be positive")

    payload = args.relation_ledger.read_bytes()
    ledger = json.loads(payload)
    if not isinstance(ledger, dict):
        raise ValueError("relation ledger must be a JSON object")
    selected = select_feedback_vertices(
        ledger,
        minimum_occurrences=args.minimum_occurrences,
        maximum_rational_prime=args.maximum_rational_prime,
        maximum_seeds=args.maximum_seeds,
    )
    seed_specials = ",".join(
        f"{record['rational_prime']}:{record['residue']}" for record in selected
    )
    output = {
        "schema": "elliptic-curves.bnf-free-minkowski-feedback-specials.v1",
        "status": "EXACT_GRAPH_FEEDBACK_SELECTION_NOT_AN_ARITHMETIC_BOUND",
        "relation_ledger": {
            "path": str(args.relation_ledger.resolve()),
            "sha256": sha256(payload).hexdigest(),
        },
        "minimum_occurrences": args.minimum_occurrences,
        "maximum_rational_prime": args.maximum_rational_prime,
        "maximum_seeds": args.maximum_seeds,
        "selected_special_ideals": selected,
        "seed_specials_argument": seed_specials,
        "recommended_special_ideal_mode": (
            "cycle-pairs" if len(selected) >= 3 else "single"
        ),
        "recommended_pair_cycle_length": len(selected) if len(selected) >= 3 else None,
        "claim_boundary": (
            "The selected ideals are exact repeated graph vertices. Selection "
            "does not close a relation and supplies no class-group or Selmer bound."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        f"{PROTOCOL}|stage=complete|selected={len(selected)}"
        f"|seed_specials={seed_specials}|output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
