#!/usr/bin/env sage
"""Certify a mod-two S-class quotient from a BNF-free relation ledger.

The ledger contains exact principal-ideal relations, but a relation-rank
plateau alone says nothing about completeness.  This audit supplies the
missing finite-generation gate.  It checks every stored principal-generator
product and then uses either:

* the unconditional Minkowski prime-ideal generation bound; or
* Bach's ``12 (log |Delta_K|)^2`` bound, conditional on ERH/GRH.

When the audited factor base contains every prime ideal through the selected
bound, the quotient of its GF(2) vector space by the verified principal rows
and the declared S-primes surjects onto ``Cl(O_K[S^-1]) / 2``.  Its dimension
is consequently an upper bound.  No full BNF, regulator, or odd class-group
part is computed.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from sage.all import (
    NumberField,
    PolynomialRing,
    QQ,
    RealIntervalField,
    ZZ,
    factorial,
    pari,
    prime_range,
)


PROTOCOL = "BNFFREECLASS"
SCHEMA = "elliptic-curves.bnf-free-principal-relation-ledger.v1"
CANONICAL_SOURCE = "canonical_rational_prime_principal"


def rational(value: str) -> QQ:
    value = Fraction(value)
    return QQ(value.numerator) / QQ(value.denominator)


def prime_signature(prime) -> tuple[str, int, int, int]:
    """Stable enough identity for a prime ideal in the ledger's field model."""

    return (
        str(prime.pari_hnf()),
        int(prime.norm()),
        int(prime.residue_class_degree()),
        int(prime.smallest_integer()),
    )


def packed_rank(rows: list[int]) -> int:
    pivots: dict[int, int] = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            previous = pivots.get(pivot)
            if previous is None:
                pivots[pivot] = row
                break
            row ^= previous
    return len(pivots)


def coefficient_field(ledger: dict):
    coefficients = [
        rational(value) for value in ledger["defining_polynomial_ascending"]
    ]
    if len(coefficients) != 4 or coefficients[-1] != 1:
        raise ValueError("ledger must define a monic cubic polynomial")
    if ledger.get("curve_preset") == "elkies-2026-rank28":
        from run_elkies_2026_rank28_s_class_pari import validate_inputs

        _, _, certified_coefficients, factor_hint_primes = validate_inputs()
        if coefficients != [QQ(value) for value in certified_coefficients]:
            raise ValueError("Elkies rank-28 ledger has the wrong defining cubic")
        pari.addprimes(factor_hint_primes)
    ring = PolynomialRing(QQ, "x")
    x = ring.gen()
    return NumberField(sum(value * x**index for index, value in enumerate(coefficients)), "theta")


def required_prime_ideals(field, rational_bound: int, s_rational_primes: set[int]):
    rational_primes = set(int(q) for q in prime_range(2, rational_bound + 1))
    rational_primes.update(s_rational_primes)
    primes = []
    for rational_prime in sorted(rational_primes):
        primes.extend(field.primes_above(ZZ(rational_prime)))
    return tuple(primes)


def generator(field, record: dict):
    coordinates = [rational(value) for value in record["power_basis"]]
    if len(coordinates) != 3:
        raise ValueError("principal generator must have three power-basis coordinates")
    theta = field.gen()
    return sum(value * theta**index for index, value in enumerate(coordinates))


def bounds(field) -> tuple[int, int]:
    """Return safe integer upper bounds for Minkowski and Bach generation."""

    degree = field.degree()
    if degree != 3:
        raise ValueError("this audit currently supports cubic ledgers only")
    _, complex_places = field.signature()
    discriminant = abs(ZZ(field.discriminant()))
    real_interval = RealIntervalField(256)
    discriminant_interval = real_interval(discriminant)
    minkowski = (
        real_interval(factorial(degree))
        / real_interval(degree) ** degree
        * (real_interval(4) / real_interval.pi()) ** complex_places
        * discriminant_interval.sqrt()
    )
    bach = real_interval(12) * discriminant_interval.log() ** 2
    # Use the upper endpoint, so finite-precision interval arithmetic cannot
    # accidentally round a generating bound downward.
    return ZZ(minkowski.upper().ceil()), ZZ(bach.upper().ceil())


def canonical_principal_relation_audit(field, ledger: dict, factor_base_column: dict) -> dict:
    """Verify optional stored ``(p)`` rows, without making them a validity gate.

    These rows are cheap but often omitted by generic relation collectors. A
    class-quotient bound remains valid without them, just weaker; this record
    makes the distinction visible in the certificate.
    """

    expected_primes = sorted({signature[3] for signature in factor_base_column})
    found = {}
    for relation_index, relation in enumerate(ledger["closed_relations"]):
        if relation.get("source") != CANONICAL_SOURCE:
            continue
        rational_prime = int(relation.get("rational_prime", -1))
        if rational_prime in found:
            raise ValueError(f"duplicate canonical rational-prime row for {rational_prime}")
        indices = [int(index) for index in relation["generator_indices"]]
        if len(indices) != 1:
            raise ValueError(f"canonical row {rational_prime} must retain one generator")
        stored = generator(field, ledger["generators"][indices[0]])
        if stored != field(rational_prime):
            raise ValueError(f"canonical row {rational_prime} does not retain generator p")
        expected_row = 0
        for prime, exponent in field.ideal(ZZ(rational_prime)).factor():
            if int(exponent) & 1:
                column = factor_base_column.get(prime_signature(prime))
                if column is None:
                    raise ValueError(
                        f"canonical row {rational_prime} uses a prime outside the factor base"
                    )
                expected_row |= 1 << column
        actual_row = int(relation["fb_parity_mask_hex"], 16)
        if actual_row != expected_row:
            raise ValueError(f"canonical row {rational_prime} has the wrong valuation parity")
        found[rational_prime] = relation_index

    metadata = ledger.get("canonical_principal_relations")
    metadata_present = metadata is not None
    if metadata_present:
        if not isinstance(metadata, dict):
            raise ValueError("canonical_principal_relations metadata is not an object")
        recorded_primes = sorted(int(value) for value in metadata.get("rational_primes", ()))
        if recorded_primes != expected_primes:
            raise ValueError("canonical-principal metadata does not cover the factor-base primes")
        if int(metadata.get("completed_relation_count", -1)) != len(expected_primes):
            raise ValueError("canonical-principal metadata has the wrong row count")
        if not bool(metadata.get("principal_generators_stored")):
            raise ValueError("canonical-principal metadata does not attest stored generators")

    found_primes = sorted(found)
    missing = sorted(set(expected_primes).difference(found))
    extras = sorted(set(found).difference(expected_primes))
    if extras:
        raise ValueError("canonical-principal rows include a rational prime outside the factor base")
    if metadata_present and missing:
        raise ValueError("canonical-principal metadata is present but rows are missing")
    return {
        "status": (
            "COMPLETE_AND_VERIFIED"
            if not missing
            else "NOT_INCLUDED"
        ),
        "metadata_present": metadata_present,
        "factor_base_rational_prime_count": len(expected_primes),
        "verified_row_count": len(found_primes),
        "missing_rational_primes": missing,
        "interpretation": (
            "The canonical (p) rows are exact free principal relations. Their "
            "absence does not invalidate the class-quotient upper bound, but "
            "it can make that bound much weaker."
        ),
    }


def audit(ledger: dict, *, assume_erh: bool, verify_relations: bool) -> dict:
    if ledger.get("schema") != SCHEMA:
        raise ValueError("unexpected relation-ledger schema")
    required = (
        "factor_base_bound",
        "factor_base_completion",
        "selmer_rational_primes",
        "factor_base",
        "S_columns",
        "generators",
        "closed_relations",
    )
    missing = [key for key in required if key not in ledger]
    if missing:
        raise ValueError(
            "ledger lacks certification metadata: " + ", ".join(missing)
        )
    field = coefficient_field(ledger)
    reported_discriminant = ledger.get("field_discriminant")
    if reported_discriminant is not None and ZZ(reported_discriminant) != field.discriminant():
        raise ValueError("ledger field_discriminant does not match the defining polynomial")

    factor_base_bound = int(ledger["factor_base_bound"])
    completion = ledger["factor_base_completion"]
    if not isinstance(completion, dict) or int(
        completion.get("all_prime_ideals_above_rational_primes_through", -1)
    ) != factor_base_bound:
        raise ValueError("ledger does not attest a complete rational-prime factor base")
    factor_base_materialized = bool(
        completion.get("materialized_complete_factor_base")
    )
    s_rational_primes = {int(value) for value in ledger["selmer_rational_primes"]}
    if set(int(value) for value in completion.get("extra_declared_S_rational_primes", ())) != s_rational_primes:
        raise ValueError("ledger S-prime metadata is inconsistent")

    required_by_signature = {}
    if factor_base_materialized:
        required_by_signature = {
            prime_signature(prime): prime
            for prime in required_prime_ideals(
                field, factor_base_bound, s_rational_primes
            )
        }
    factor_base_records = ledger["factor_base"]
    if not isinstance(factor_base_records, list):
        raise ValueError("factor_base must be a list")
    factor_base_signatures = []
    factor_base_primes = []
    for record in factor_base_records:
        if not isinstance(record, dict):
            raise ValueError("factor-base entry is not an object")
        signature = (
            str(record["hnf"]),
            int(record["norm"]),
            int(record["residue_degree"]),
            int(record["rational_prime"]),
        )
        if factor_base_materialized:
            prime = required_by_signature.get(signature)
        else:
            prime = next(
                (
                    candidate
                    for candidate in field.primes_above(ZZ(signature[3]))
                    if prime_signature(candidate) == signature
                ),
                None,
            )
        if prime is None:
            raise ValueError("factor-base prime cannot be reconstructed from its ledger signature")
        factor_base_signatures.append(signature)
        factor_base_primes.append(prime)
    if len(set(factor_base_signatures)) != len(factor_base_signatures):
        raise ValueError("factor base has duplicate prime ideals")
    factor_base_column = {
        signature: index for index, signature in enumerate(factor_base_signatures)
    }
    missing_primes = set(required_by_signature).difference(factor_base_signatures)
    if factor_base_materialized and missing_primes:
        raise ValueError(
            f"factor base omits {len(missing_primes)} required prime ideals"
        )

    s_columns = [int(index) for index in ledger["S_columns"]]
    if len(s_columns) != len(set(s_columns)) or any(
        index < 0 or index >= len(factor_base_primes) for index in s_columns
    ):
        raise ValueError("invalid S_columns")
    actual_s = {factor_base_signatures[index] for index in s_columns}
    if factor_base_materialized:
        expected_s = {
            prime_signature(prime)
            for rational_prime in s_rational_primes
            for prime in field.primes_above(ZZ(rational_prime))
        }
        s_columns_complete = actual_s == expected_s
    else:
        s_columns_complete = False
    if factor_base_materialized and not s_columns_complete:
        raise ValueError("S_columns are not exactly the declared Selmer prime ideals")

    relation_rows = []
    canonical_relation_rows = []
    relation_count = len(ledger["closed_relations"])
    for relation_index, relation in enumerate(ledger["closed_relations"]):
        if not isinstance(relation, dict):
            raise ValueError("closed relation is not an object")
        row = int(relation["fb_parity_mask_hex"], 16)
        if row >> len(factor_base_primes):
            raise ValueError(f"relation {relation_index} exceeds factor-base dimension")
        indices = [int(index) for index in relation["generator_indices"]]
        if len(indices) != len(set(indices)) or any(
            index < 0 or index >= len(ledger["generators"]) for index in indices
        ):
            raise ValueError(f"relation {relation_index} has invalid generator indices")
        if verify_relations:
            product = field(1)
            for index in indices:
                product *= generator(field, ledger["generators"][index])
            # Factoring this one principal ideal is sparse in the relation
            # support.  Iterating its factors is far cheaper than querying a
            # valuation at every prime in an ERH-sized factor base.
            computed = 0
            for prime, exponent in field.ideal(product).factor():
                if int(exponent) & 1:
                    column = factor_base_column.get(prime_signature(prime))
                    if column is None:
                        raise ValueError(
                            f"relation {relation_index} has odd support outside the factor base"
                        )
                    computed |= 1 << column
            if computed != row:
                raise ValueError(
                    f"relation {relation_index} does not match its principal-generator product"
                )
        relation_rows.append(row)
        if relation.get("source") == CANONICAL_SOURCE:
            canonical_relation_rows.append(row)

    relation_rank = packed_rank(relation_rows)
    quotient_rank = packed_rank(
        relation_rows + [(1 << column) for column in s_columns]
    )
    quotient_dimension = len(factor_base_primes) - quotient_rank
    canonical_relation_rank = packed_rank(canonical_relation_rows)
    canonical_quotient_rank = packed_rank(
        canonical_relation_rows + [(1 << column) for column in s_columns]
    )
    canonical_quotient_dimension = len(factor_base_primes) - canonical_quotient_rank
    canonical_audit = canonical_principal_relation_audit(
        field, ledger, factor_base_column
    )
    minkowski_bound, bach_bound = bounds(field)

    if not factor_base_materialized:
        classification = "UNCERTIFIED_INCOMPLETE_FACTOR_BASE"
        certification = {
            "status": "UNCERTIFIED_RELATION_STABILIZATION",
            "method": "none",
            "remaining_dimension_upper_bound": None,
            "hypothesis": None,
            "reason": "FACTOR_BASE_NOT_MATERIALIZED_THROUGH_DECLARED_BOUND",
        }
    elif factor_base_bound >= minkowski_bound:
        classification = "CERTIFIED_UNCONDITIONAL_MINKOWSKI"
        certification = {
            "status": "CERTIFIED",
            "method": "Minkowski prime-ideal generation bound plus verified principal relation rows",
            "remaining_dimension_upper_bound": quotient_dimension,
            "hypothesis": None,
        }
    elif assume_erh and factor_base_bound >= bach_bound:
        classification = "CERTIFIED_UNDER_ERH_BACH"
        certification = {
            "status": "CERTIFIED_UNDER_HYPOTHESIS",
            "method": "Bach 12(log |Delta_K|)^2 prime-ideal generation bound plus verified principal relation rows",
            "remaining_dimension_upper_bound": quotient_dimension,
            "hypothesis": "ERH/GRH for the Dedekind zeta function of the cubic field",
        }
    else:
        classification = "UNCERTIFIED_FACTOR_BASE"
        reason = (
            "ERH_BACH_BOUND_NOT_ACCEPTED"
            if factor_base_bound >= bach_bound
            else "FACTOR_BASE_BELOW_BACH_BOUND"
        )
        certification = {
            "status": "UNCERTIFIED_RELATION_STABILIZATION",
            "method": "none",
            "remaining_dimension_upper_bound": None,
            "hypothesis": None,
            "reason": reason,
        }

    real_places, complex_places = field.signature()
    # For a cubic field, -1 is a non-square torsion unit, so the finite-S
    # unit contribution modulo squares is r1+r2+#S.
    s_unit_squareclass_dimension = int(
        real_places + complex_places + len(s_columns)
    )
    output = {
        "protocol": "BNFFREECLASS-v1",
        "classification": classification,
        "collector": {
            "status": ledger.get("status"),
            "curve_preset": ledger.get("curve_preset"),
            "special_ideal_mode": ledger.get("special_ideal_mode"),
            "special_residue_degree": ledger.get("special_residue_degree"),
            "special_primes_in_factor_base": ledger.get(
                "special_primes_in_factor_base"
            ),
            "sampled_generator_count": len(ledger.get("generators", ())),
            "noncanonical_closed_relation_count": sum(
                relation.get("source") != CANONICAL_SOURCE
                for relation in ledger.get("closed_relations", ())
            ),
            "factor_hint_certificate": ledger.get("factor_hint_certificate"),
            "early_quotient": ledger.get("collection_early_quotient"),
        },
        "defining_polynomial_ascending": ledger["defining_polynomial_ascending"],
        "field_discriminant": str(field.discriminant()),
        "field_degree": int(field.degree()),
        "signature": [int(real_places), int(complex_places)],
        "factor_base_bound": factor_base_bound,
        "factor_base_materialized_complete": factor_base_materialized,
        "factor_base_dimension": len(factor_base_primes),
        "declared_S_prime_ideal_count": len(s_columns),
        "declared_S_prime_ideal_columns_complete": s_columns_complete,
        "declared_S_rational_primes": sorted(s_rational_primes),
        "relation_count": relation_count,
        "relation_rank": relation_rank,
        "factor_base_quotient_dimension": quotient_dimension,
        "relation_source_rank_analysis": {
            "canonical_rational_principal_relation_count": len(canonical_relation_rows),
            "canonical_rational_principal_relation_rank": canonical_relation_rank,
            "factor_base_quotient_dimension_after_canonical_rows_and_S": canonical_quotient_dimension,
            "quotient_dimension_improvement_from_noncanonical_rows": (
                canonical_quotient_dimension - quotient_dimension
            ),
            "interpretation": (
                "This isolates the reduction supplied by collected noncanonical "
                "relations after the free exact (p) rows and S-columns have "
                "already been inserted."
            ),
        },
        "minkowski_generation_bound": int(minkowski_bound),
        "bach_erh_generation_bound": int(bach_bound),
        "principal_relations_verified": verify_relations,
        "canonical_rational_principal_relation_audit": canonical_audit,
        "class_quotient_certification": certification,
        "s_unit_squareclass_dimension": s_unit_squareclass_dimension,
        "k_s_2_dimension_upper_bound": (
            s_unit_squareclass_dimension + quotient_dimension
            if certification["remaining_dimension_upper_bound"] is not None
            else None
        ),
        "interpretation": (
            "The stated bound concerns Cl(O_K[S^-1])/2 and, when certified, "
            "the resulting K(S,2) envelope. It does not impose the norm or "
            "local conditions of elliptic 2-descent."
        ),
        "selmer_claim": {
            "completed": False,
            "residual_two_selmer_quotient_dimension": None,
            "all_local_solubility_conditions_completed": False,
            "expensive_search_authorized": False,
        },
    }
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--relation-ledger", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--assume-erh",
        action="store_true",
        help="permit the explicitly conditional Bach factor-base generation gate",
    )
    parser.add_argument(
        "--skip-principal-relation-verification",
        action="store_true",
        help="diagnostic only; never emits a certificate",
    )
    args = parser.parse_args()
    if args.skip_principal_relation_verification:
        raise ValueError(
            "a quotient certification requires principal-generator verification"
        )
    ledger = json.loads(args.relation_ledger.read_text())
    if not isinstance(ledger, dict):
        raise ValueError("relation ledger must be a JSON object")
    output = audit(ledger, assume_erh=args.assume_erh, verify_relations=True)
    output["relation_ledger"] = {
        "path": str(args.relation_ledger.resolve()),
        "sha256": sha256(args.relation_ledger.read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        f"{PROTOCOL}|stage=complete|classification={output['classification']}"
        f"|factor_base_dimension={output['factor_base_dimension']}"
        f"|relation_rank={output['relation_rank']}"
        f"|factor_base_quotient_dimension={output['factor_base_quotient_dimension']}"
        f"|bach_erh_generation_bound={output['bach_erh_generation_bound']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
