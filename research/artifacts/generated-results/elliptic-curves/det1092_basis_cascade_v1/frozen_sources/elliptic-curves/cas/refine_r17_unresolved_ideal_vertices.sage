#!/usr/bin/env sage
"""Refine cached R17 norm cofactors into exact residual-ideal incidence edges.

The input is a bounded Minkowski relation ledger.  A product-tree GCD pass
finds proved rational prime factors shared by its unresolved norm cofactors.
For each affected generator this program removes every exact factor-base
valuation and every proved shared prime-ideal valuation from the principal
ideal.  The unfactored remainder is retained as a reduced ideal HNF vertex;
it is never declared prime and its norm is never required to factor.

An incidence dependency is therefore an exact ideal-class relation even when
every individual edge has an unfactored tail.  Vertices exactly equivalent to
one of the certified point half-ideals are killed from the start.  Without a
factor-base generation certificate, any resulting dimension is still only a
bounded quotient presentation and not an S-class or Selmer upper bound.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from hashlib import sha256
import json
from pathlib import Path
import time

from sage.all import NumberField, PolynomialRing, QQ, ZZ, pari

from run_fermigier_rank20_fixedfb_quadratic_specialq import (
    SparseLargePrimeEliminator,
)
from run_fermigier_rank20_minkowski_specialq import (
    select_batch_gcd_records,
    split_shared_cofactors,
)


ROOT = Path(__file__).resolve().parents[2]
PRESSURE = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-kummer-classgroup-pressure-v1.json"
)
PRESSURE_STATUS = "PROVED_KUMMER_FORCED_CUBIC_CLASS_GROUP_2RANK_LOWER_BOUNDS"
SCHEMA = "elliptic-curves.r17-unresolved-ideal-vertices.v1"
PROTOCOL = "R17IDEALVERT"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def insert_row(pivots: dict[int, int], row: int) -> bool:
    while row:
        pivot = row.bit_length() - 1
        previous = pivots.get(pivot)
        if previous is None:
            pivots[pivot] = row
            return True
        row ^= previous
    return False


def parity_vertices(vertices):
    counts = Counter(vertices)
    return tuple(sorted(vertex for vertex, count in counts.items() if count & 1))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-inputs", type=int, default=10000)
    parser.add_argument(
        "--max-shared-records",
        type=int,
        default=0,
        help="zero processes every record with a proper shared factor",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.max_inputs <= 1:
        parser.error("--max-inputs must exceed one")
    if args.max_shared_records < 0:
        parser.error("--max-shared-records must be nonnegative")
    if args.output.exists() and not args.overwrite:
        raise FileExistsError(args.output)

    started = time.monotonic()
    ledger = json.loads(args.input.read_text())
    if ledger.get("schema") != "elliptic-curves.bnf-free-principal-relation-ledger.v1":
        raise ValueError("unexpected relation-ledger schema")
    if ledger.get("status") != "exact_minkowski_ideal_relations_not_class_group_completion":
        raise ValueError("the relation ledger is not an exact Minkowski checkpoint")
    if not ledger.get("special_primes_in_factor_base"):
        raise ValueError("this refinement requires every source special ideal in the factor base")
    if ledger.get("closed_relations") or ledger.get("partial_relations"):
        raise ValueError("merge the pre-existing exact relations before this zero-relation refinement")

    coefficients = [ZZ(value) for value in ledger["defining_polynomial_ascending"]]
    if len(coefficients) != 4 or coefficients[-1] != 1:
        raise ValueError("the input does not define a monic cubic")
    ring = PolynomialRing(QQ, "x")
    x = ring.gen()
    polynomial = sum(value * x**index for index, value in enumerate(coefficients))
    declared_s = [ZZ(value) for value in ledger["selmer_rational_primes"]]
    pari.addprimes(declared_s)
    field = NumberField(polynomial, "theta")
    theta = field.gen()
    nf = field.pari_nf()

    factor_base_records = ledger["factor_base"]
    factor_base = [None] * len(factor_base_records)
    by_rational_prime = defaultdict(list)
    # The Minkowski ledger fixes list order as column order.
    for index, record in enumerate(factor_base_records):
        rational_prime = ZZ(record["rational_prime"])
        expected_hnf = record["hnf"]
        match = next(
            (
                prime
                for prime in field.primes_above(rational_prime)
                if str(prime.pari_hnf()) == expected_hnf
            ),
            None,
        )
        if match is None:
            raise ArithmeticError("a factor-base prime ideal no longer reconstructs")
        factor_base[index] = match
        by_rational_prime[rational_prime].append((index, match))
    if any(prime is None for prime in factor_base):
        raise ArithmeticError("the factor-base reconstruction is incomplete")

    pressure = json.loads(PRESSURE.read_text())
    if pressure.get("status") != PRESSURE_STATUS:
        raise ArithmeticError("the point half-ideal certificate is not passing")
    pressure_record = None
    for record in pressure["curves"]:
        certified_polynomial = ring(record["two_division_cubic"].replace("z", "x"))
        if certified_polynomial == polynomial:
            pressure_record = record
            break
    if pressure_record is None:
        raise ValueError("the cubic does not match a certified R17 record fibre")
    curve_id = int(pressure_record["curve_id"])

    unit_hnf = str(field.ideal(1).pari_hnf())

    def reduced_ideal_record(ideal):
        source_hnf = ideal.pari_hnf()
        reduced = pari.idealred(nf, [source_hnf, 1])
        reduced_hnf = pari.idealhnf(nf, reduced[0])
        multiplier = reduced[1]
        expected = pari.idealmul(nf, reduced_hnf, multiplier)
        if str(pari.idealhnf(nf, expected)) != str(source_hnf):
            raise ArithmeticError("idealred residual-ideal identity failed")
        return str(reduced_hnf), str(multiplier)

    killed_reduced_hnfs = {unit_hnf}
    killed_half_ideal_records = []
    for record in pressure_record["point_half_ideals"]:
        half_ideal = field.ideal(pari(record["half_ideal_hnf"]))
        reduced_hnf, multiplier = reduced_ideal_record(half_ideal)
        killed_reduced_hnfs.add(reduced_hnf)
        killed_half_ideal_records.append(
            {
                "label": record["label"],
                "half_ideal_hnf": record["half_ideal_hnf"],
                "reduced_ideal_hnf": reduced_hnf,
                "reduction_multiplier": multiplier,
            }
        )

    canonical_pivots = {}
    for rational_prime, entries in sorted(by_rational_prime.items()):
        row = 0
        for index, prime in entries:
            exponent = int(field.ideal(rational_prime).valuation(prime))
            if exponent & 1:
                row ^= 1 << index
        insert_row(canonical_pivots, row)
    quotient_pivots = dict(canonical_pivots)
    for column in ledger["S_columns"]:
        insert_row(quotient_pivots, 1 << int(column))
    baseline_rank = len(quotient_pivots)
    baseline_dimension = len(factor_base) - baseline_rank

    unresolved = ledger["unresolved_cofactors"]
    selected = select_batch_gcd_records(unresolved, args.max_inputs)
    cofactors = [int(record["cofactor"]) for _index, record in selected]
    parts, gcd_statistics = split_shared_cofactors(cofactors, "product-tree")
    shared = [
        (source_index, record, factor_parts)
        for (source_index, record), factor_parts in zip(selected, parts)
        if len(factor_parts) > 1
    ]
    if args.max_shared_records:
        shared = shared[: args.max_shared_records]

    sparse = SparseLargePrimeEliminator()
    partial_edges = []
    closed_relations = []
    skipped_no_proved_shared_prime = 0
    killed_tail_count = 0
    factor_base_tail_count = 0

    for shared_index, (source_index, record, factor_parts) in enumerate(shared, 1):
        prime_parts = []
        for part in sorted(set(map(ZZ, factor_parts))):
            if part.is_pseudoprime() and part.is_prime():
                prime_parts.append(part)
        known_residual = [ZZ(value) for value in record["known_odd_residual_primes"]]
        split_primes = sorted(set(prime_parts + known_residual))
        if not prime_parts:
            skipped_no_proved_shared_prime += 1
            continue

        generator_index = int(record["generator_index"])
        generator = ledger["generators"][generator_index]
        coordinates = [QQ(value) for value in generator["power_basis"]]
        alpha = sum(
            coordinates[index] * theta**index for index in range(len(coordinates))
        )
        principal = field.ideal(alpha)
        residual = principal
        row = 0
        removed_factor_ideals = []

        for rational_prime, entries in by_rational_prime.items():
            for column, prime in entries:
                valuation = int(alpha.valuation(prime))
                if valuation < 0:
                    raise ArithmeticError("a cached Minkowski generator is not integral")
                if not valuation:
                    continue
                residual /= prime**valuation
                if valuation & 1:
                    row ^= 1 << column
                removed_factor_ideals.append(
                    {
                        "kind": "factor_base",
                        "column": column,
                        "rational_prime": str(rational_prime),
                        "prime_ideal_hnf": str(prime.pari_hnf()),
                        "valuation": valuation,
                    }
                )

        large_vertices = []
        for rational_prime in split_primes:
            for prime in field.primes_above(rational_prime):
                valuation = int(alpha.valuation(prime))
                if valuation < 0:
                    raise ArithmeticError("a split residual valuation is negative")
                if not valuation:
                    continue
                residual /= prime**valuation
                if valuation & 1:
                    large_vertices.append(("prime_ideal", str(prime.pari_hnf())))
                removed_factor_ideals.append(
                    {
                        "kind": "proved_shared_prime",
                        "rational_prime": str(rational_prime),
                        "prime_ideal_hnf": str(prime.pari_hnf()),
                        "valuation": valuation,
                    }
                )

        if not residual.is_integral():
            raise ArithmeticError("the exact residual ideal is not integral")
        reconstructed = residual
        for removed in removed_factor_ideals:
            prime = field.ideal(pari(removed["prime_ideal_hnf"]))
            reconstructed *= prime ** int(removed["valuation"])
        if reconstructed != principal:
            raise ArithmeticError("the residual-ideal factorization does not replay")

        tail_hnf, tail_multiplier = reduced_ideal_record(residual)
        tail_kind = "opaque_residual_ideal"
        factor_base_tail = next(
            (
                index
                for index, prime in enumerate(factor_base)
                if str(prime.pari_hnf()) == tail_hnf
            ),
            None,
        )
        if tail_hnf in killed_reduced_hnfs:
            killed_tail_count += 1
            tail_kind = "known_mw29_or_unit"
        elif factor_base_tail is not None:
            factor_base_tail_count += 1
            row ^= 1 << factor_base_tail
            tail_kind = "factor_base_prime"
        else:
            large_vertices.append(("residual_ideal", tail_hnf))

        vertices = parity_vertices(large_vertices)
        edge = {
            "source_unresolved_index": source_index,
            "generator_index": generator_index,
            "fb_parity_mask_hex": hex(row),
            "large_prime_vertices": [list(vertex) for vertex in vertices],
            "proved_shared_rational_primes": [str(value) for value in prime_parts],
            "factor_parts": [str(value) for value in factor_parts],
            "removed_factor_ideals": removed_factor_ideals,
            "residual_ideal_hnf": str(residual.pari_hnf()),
            "reduced_tail_hnf": tail_hnf,
            "tail_reduction_multiplier": tail_multiplier,
            "tail_classification": tail_kind,
        }
        partial_edges.append(edge)
        cycle, provenance = sparse.add(vertices, row, generator_index)
        if cycle is not None:
            gained = insert_row(quotient_pivots, cycle)
            closed_relations.append(
                {
                    "fb_parity_mask_hex": hex(cycle),
                    "generator_indices": sorted(provenance),
                    "increased_quotient_rank": gained,
                }
            )
        if shared_index % 500 == 0:
            print(
                f"{PROTOCOL}|curve={curve_id}|processed={shared_index}/{len(shared)}"
                f"|edges={len(partial_edges)}|cycles={len(closed_relations)}"
                f"|rank_gain={len(quotient_pivots)-baseline_rank}",
                flush=True,
            )

    output = {
        "schema": SCHEMA,
        "status": "BOUNDED_EXACT_RESIDUAL_IDEAL_INCIDENCE_NOT_CLASS_GROUP_COMPLETION",
        "curve_id": curve_id,
        "input": {"path": str(args.input.resolve()), "sha256": digest(args.input)},
        "pressure_certificate": {
            "path": str(PRESSURE.relative_to(ROOT)),
            "sha256": digest(PRESSURE),
            "known_half_ideal_count": len(killed_half_ideal_records),
        },
        "settings": {
            "max_inputs": args.max_inputs,
            "max_shared_records": args.max_shared_records,
            "selection": "lowest cofactor bit length, then generator index",
            "gcd_engine": "product-tree without generic cofactor factorization",
        },
        "factor_base_width": len(factor_base),
        "baseline_rank_after_canonical_and_S_rows": baseline_rank,
        "baseline_bounded_quotient_dimension": baseline_dimension,
        "selected_unresolved_count": len(selected),
        "proper_shared_factor_record_count": len(shared),
        "gcd_statistics": gcd_statistics,
        "skipped_no_proved_shared_prime": skipped_no_proved_shared_prime,
        "killed_tail_count": killed_tail_count,
        "factor_base_tail_count": factor_base_tail_count,
        "partial_edges": partial_edges,
        "closed_relations": closed_relations,
        "large_prime_elimination": {
            "vertex_count": len(sparse.vertex_columns),
            "edge_count": sparse.edge_count,
            "rank": len(sparse.pivots),
            "dependency_count": sparse.dependency_count,
            "nullity": sparse.edge_count - len(sparse.pivots),
        },
        "bounded_quotient_relation_rank_gain": len(quotient_pivots) - baseline_rank,
        "bounded_quotient_dimension_after_cycles": (
            len(factor_base) - len(quotient_pivots)
        ),
        "killed_half_ideal_reductions": killed_half_ideal_records,
        "elapsed_seconds": time.monotonic() - started,
        "claim_boundary": [
            "Every retained edge is verified as an exact ideal factorization with an unfactored residual ideal HNF.",
            "Only exact sparse-incidence dependencies create factor-base relations.",
            "The factor base has no generation certificate, so its displayed quotient dimension is not an S-class or Selmer upper bound.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    temporary.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    temporary.replace(args.output)
    print(
        f"{PROTOCOL}|curve={curve_id}|selected={len(selected)}|shared={len(shared)}"
        f"|edges={len(partial_edges)}|cycles={len(closed_relations)}"
        f"|rank_gain={len(quotient_pivots)-baseline_rank}"
        f"|dimension={len(factor_base)-len(quotient_pivots)}"
        f"|seconds={output['elapsed_seconds']:.3f}|output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
