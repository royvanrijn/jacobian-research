#!/usr/bin/env sage
"""Add missing principal (p) rows to a cached targeted R17 ideal closure.

This is a cache-only postprocessor: it performs no norm factorization and no
new lattice search.  It replays every cached target ideal factorization, then
adds

    (p) = product_{P|p} P^e(P/p)

for each outside rational prime actually used by the sparse incidence graph.
The output remains a bounded relation certificate, not a class-group or
Selmer upper bound.
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


INPUT_SCHEMA = "elliptic-curves.r17-targeted-residual-ideal-closure.v1"
INPUT_STATUS = "BOUNDED_EXACT_MW29_RELATIVE_IDEAL_RELATIONS_NOT_CLASS_GROUP_COMPLETION"
SOURCE_SCHEMA = "elliptic-curves.r17-unresolved-ideal-vertices.v1"
SCHEMA = "elliptic-curves.r17-targeted-residual-ideal-closure-canonical.v1"
PROTOCOL = "R17IDEALCANON"


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
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.output.exists() and not args.overwrite:
        raise FileExistsError(args.output)

    started = time.monotonic()
    targeted = json.loads(args.input.read_text())
    if targeted.get("schema") != INPUT_SCHEMA or targeted.get("status") != INPUT_STATUS:
        raise ValueError("unexpected targeted-closure input")
    source_path = Path(targeted["input"]["path"])
    if digest(source_path) != targeted["input"]["sha256"]:
        raise ArithmeticError("the residual-ideal source hash changed")
    source = json.loads(source_path.read_text())
    if source.get("schema") != SOURCE_SCHEMA:
        raise ValueError("unexpected residual-ideal source schema")
    ledger_path = Path(source["input"]["path"])
    if digest(ledger_path) != source["input"]["sha256"]:
        raise ArithmeticError("the Minkowski ledger hash changed")
    ledger = json.loads(ledger_path.read_text())

    ring = PolynomialRing(QQ, "x")
    x = ring.gen()
    polynomial = sum(
        ZZ(value) * x**index
        for index, value in enumerate(ledger["defining_polynomial_ascending"])
    )
    pari.addprimes([ZZ(value) for value in ledger["selmer_rational_primes"]])
    field = NumberField(polynomial, "theta")

    factor_base = []
    factor_base_by_hnf = {}
    by_rational_prime = defaultdict(list)
    for index, record in enumerate(ledger["factor_base"]):
        rational_prime = ZZ(record["rational_prime"])
        expected_hnf = record["hnf"]
        prime = next(
            (
                candidate
                for candidate in field.primes_above(rational_prime)
                if str(candidate.pari_hnf()) == expected_hnf
            ),
            None,
        )
        if prime is None:
            raise ArithmeticError("a factor-base ideal no longer reconstructs")
        factor_base.append(prime)
        factor_base_by_hnf[expected_hnf] = index
        by_rational_prime[rational_prime].append((index, prime))

    quotient_pivots = {}
    for rational_prime, entries in sorted(by_rational_prime.items()):
        row = 0
        rational_ideal = field.ideal(rational_prime)
        for index, prime in entries:
            if int(rational_ideal.valuation(prime)) & 1:
                row ^= 1 << index
        insert_row(quotient_pivots, row)
    for column in ledger["S_columns"]:
        insert_row(quotient_pivots, 1 << int(column))
    baseline_rank = len(quotient_pivots)

    sparse = SparseLargePrimeEliminator()
    references = {}
    cycles = []
    next_edge_id = 0

    def add_edge(vertices, row, reference):
        nonlocal next_edge_id
        edge_id = next_edge_id
        next_edge_id += 1
        references[edge_id] = reference
        cycle, provenance = sparse.add(parity_vertices(vertices), row, edge_id)
        if cycle is not None:
            cycles.append(
                {
                    "fb_parity_mask_hex": hex(cycle),
                    "edge_references": [references[index] for index in sorted(provenance)],
                    "increased_quotient_rank": insert_row(quotient_pivots, cycle),
                }
            )

    outside_hnf_to_rational_prime = {}
    for index, edge in enumerate(source["partial_edges"]):
        for removed in edge["removed_factor_ideals"]:
            outside_hnf_to_rational_prime[removed["prime_ideal_hnf"]] = int(
                removed["rational_prime"]
            )
        add_edge(
            [tuple(vertex) for vertex in edge["large_prime_vertices"]],
            int(edge["fb_parity_mask_hex"], 16),
            {"kind": "source_edge", "index": index},
        )

    known_by_label = {
        record["label"]: record
        for record in source["killed_half_ideal_reductions"]
    }

    def replay_target_edge(edge, source_ideal):
        beta = sum(
            (
                QQ(value) * field.gen() ** index
                for index, value in enumerate(edge["beta_power_basis"])
            ),
            field(0),
        )
        quotient = field.ideal(pari(edge["quotient_ideal_hnf"]))
        if source_ideal * quotient != field.ideal(beta):
            raise ArithmeticError("a cached target identity no longer replays")
        rebuilt = field.ideal(1)
        row = 0
        vertices = []
        for record in edge["prime_ideal_factorization"]:
            rational_prime = ZZ(record["rational_prime"])
            if not rational_prime.is_prime():
                raise ArithmeticError("a cached rational factor is not prime")
            prime = field.ideal(pari(record["prime_ideal_hnf"]))
            valuation = int(record["valuation"])
            rebuilt *= prime**valuation
            outside_hnf_to_rational_prime[record["prime_ideal_hnf"]] = int(
                rational_prime
            )
            if valuation & 1:
                column = factor_base_by_hnf.get(record["prime_ideal_hnf"])
                if column is None:
                    vertices.append(("prime_ideal", record["prime_ideal_hnf"]))
                else:
                    row ^= 1 << column
        if rebuilt != quotient:
            raise ArithmeticError("cached prime-ideal factors do not rebuild the quotient")
        if hex(row) != edge["fb_parity_mask_hex"]:
            raise ArithmeticError("cached target factor-base row changed")
        if [list(vertex) for vertex in parity_vertices(vertices)] != edge["outside_vertices"]:
            raise ArithmeticError("cached target outside vertices changed")
        return vertices, row

    for index, edge in enumerate(targeted["known_MW29_target_edges"]):
        source_ideal = field.ideal(
            pari(known_by_label[edge["label"]]["reduced_ideal_hnf"])
        )
        vertices, row = replay_target_edge(edge, source_ideal)
        add_edge(vertices, row, {"kind": "known_MW29_target_edge", "index": index})

    for index, edge in enumerate(targeted["tail_target_edges"]):
        source_ideal = field.ideal(pari(edge["source_reduced_ideal_hnf"]))
        vertices, row = replay_target_edge(edge, source_ideal)
        add_edge(
            [("residual_ideal", edge["source_reduced_ideal_hnf"])] + vertices,
            row,
            {"kind": "residual_ideal_target_edge", "index": index},
        )

    if cycles:
        raise ArithmeticError("the cached target stage did not have zero dependencies")

    used_outside_hnfs = {
        vertex[1]
        for edge in source["partial_edges"]
        for vertex in edge["large_prime_vertices"]
        if vertex[0] == "prime_ideal"
    }
    used_outside_hnfs.update(
        vertex[1]
        for edge in targeted["known_MW29_target_edges"]
        + targeted["tail_target_edges"]
        for vertex in edge["outside_vertices"]
    )
    if any(hnf not in outside_hnf_to_rational_prime for hnf in used_outside_hnfs):
        raise ArithmeticError("an outside ideal has no rational-prime label")

    canonical_edges = []
    for rational_prime in sorted(
        {outside_hnf_to_rational_prime[hnf] for hnf in used_outside_hnfs}
    ):
        rational_ideal = field.ideal(rational_prime)
        row = 0
        vertices = []
        factors = []
        for prime in field.primes_above(rational_prime):
            exponent = int(rational_ideal.valuation(prime))
            hnf = str(prime.pari_hnf())
            column = factor_base_by_hnf.get(hnf)
            if exponent & 1:
                if column is None:
                    vertices.append(("prime_ideal", hnf))
                else:
                    row ^= 1 << column
            factors.append(
                {
                    "prime_ideal_hnf": hnf,
                    "ramification_index": exponent,
                    "factor_base_column": column,
                }
            )
        edge_index = len(canonical_edges)
        canonical_edges.append(
            {
                "rational_prime": str(rational_prime),
                "fb_parity_mask_hex": hex(row),
                "outside_vertices": [list(vertex) for vertex in parity_vertices(vertices)],
                "prime_ideal_factorization": factors,
                "verified_identity": "(p) = product_{P|p} P^e(P/p)",
            }
        )
        add_edge(
            vertices,
            row,
            {"kind": "canonical_outside_rational_prime", "index": edge_index},
        )

    output = {
        "schema": SCHEMA,
        "status": "BOUNDED_EXACT_CANONICAL_AUGMENTATION_NOT_CLASS_GROUP_COMPLETION",
        "curve_id": targeted["curve_id"],
        "input": {"path": str(args.input.resolve()), "sha256": digest(args.input)},
        "residual_ideal_source": {
            "path": str(source_path.resolve()),
            "sha256": digest(source_path),
        },
        "source_minkowski_ledger": {
            "path": str(ledger_path.resolve()),
            "sha256": digest(ledger_path),
        },
        "factor_base_width": len(factor_base),
        "baseline_rank_after_canonical_and_S_rows": baseline_rank,
        "baseline_bounded_quotient_dimension": len(factor_base) - baseline_rank,
        "replayed_source_edge_count": len(source["partial_edges"]),
        "replayed_known_MW29_target_edge_count": len(targeted["known_MW29_target_edges"]),
        "replayed_tail_target_edge_count": len(targeted["tail_target_edges"]),
        "canonical_outside_rational_prime_edges": canonical_edges,
        "closed_relations": cycles,
        "large_prime_elimination": {
            "vertex_count": len(sparse.vertex_columns),
            "edge_count": sparse.edge_count,
            "rank": len(sparse.pivots),
            "dependency_count": sparse.dependency_count,
            "nullity": sparse.edge_count - len(sparse.pivots),
        },
        "bounded_quotient_relation_rank_gain": len(quotient_pivots) - baseline_rank,
        "bounded_quotient_dimension_after_cycles": len(factor_base) - len(quotient_pivots),
        "elapsed_seconds": time.monotonic() - started,
        "claim_boundary": [
            "All cached target factorizations and all principal outside-(p) rows replay exactly.",
            "Only sparse-incidence dependencies create factor-base relations.",
            "The factor base has no generation certificate, so the displayed quotient dimension is not an S-class or Selmer upper bound.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    temporary.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    temporary.replace(args.output)
    print(
        f"{PROTOCOL}|curve={targeted['curve_id']}|canonical={len(canonical_edges)}"
        f"|cycles={len(cycles)}|rank_gain={len(quotient_pivots)-baseline_rank}"
        f"|dimension={len(factor_base)-len(quotient_pivots)}"
        f"|seconds={output['elapsed_seconds']:.3f}|output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
