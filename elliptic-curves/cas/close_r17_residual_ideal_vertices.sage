#!/usr/bin/env sage
"""Target cached R17 residual ideals with exact directional ideal minima.

The preceding residual-ideal refiner turns an unfactored norm cofactor into an
exact ideal-class vertex.  Blind short-vector sampling leaves almost all such
vertices as leaves.  This second stage instead samples *inside each leaf
ideal* I.  For beta in I it factors the much smaller integral quotient

    J = (beta) / I,

and records the exact class relation [I] + [J] = 0 over F_2.  For a certified
MW29 half-ideal, [I] is killed and the same calculation directly supplies the
relative relation [J] = 0.  Sparse dependencies are accepted only after the
ideal factorization (beta) = I*J and every prime-ideal valuation replay.

This remains a bounded factor-base calculation.  It proves exact relations,
but it is not an S-class or Selmer upper bound without a separate generation
certificate.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from math import gcd
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import time

from sage.all import (
    ComplexField,
    NumberField,
    PolynomialRing,
    QQ,
    RealField,
    ZZ,
    factor,
    matrix,
    pari,
    vector,
)

from run_fermigier_rank20_fixedfb_quadratic_specialq import (
    SparseLargePrimeEliminator,
)


SCHEMA = "elliptic-curves.r17-targeted-residual-ideal-closure.v1"
INPUT_SCHEMA = "elliptic-curves.r17-unresolved-ideal-vertices.v1"
INPUT_STATUS = "BOUNDED_EXACT_RESIDUAL_IDEAL_INCIDENCE_NOT_CLASS_GROUP_COMPLETION"
LEDGER_SCHEMA = "elliptic-curves.bnf-free-principal-relation-ledger.v1"
PROTOCOL = "R17IDEALCLOSE"


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


def projective_coordinate_key(element):
    """Identify nonzero field elements modulo rational scaling."""

    coordinates = [QQ(value) for value in element.list()]
    denominator = ZZ(1)
    for value in coordinates:
        denominator = denominator.lcm(value.denominator())
    integral = [ZZ(value * denominator) for value in coordinates]
    content = 0
    for value in integral:
        content = gcd(content, abs(int(value)))
    if not content:
        raise ArithmeticError("the zero element has no projective key")
    primitive = [value // content for value in integral]
    first_nonzero = next(value for value in primitive if value)
    if first_nonzero < 0:
        primitive = [-value for value in primitive]
    return tuple(map(str, primitive))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--max-tail-ideals",
        type=int,
        default=0,
        help="zero targets every opaque residual ideal",
    )
    parser.add_argument("--max-samples-per-tail", type=int, default=2)
    parser.add_argument("--max-tail-factor-bits", type=int, default=120)
    parser.add_argument("--max-samples-per-known-ideal", type=int, default=1)
    parser.add_argument("--max-known-factor-bits", type=int, default=230)
    parser.add_argument(
        "--candidate-engine", choices=("idealmin", "lll"), default="lll"
    )
    parser.add_argument("--lattice-combination-bound", type=int, default=2)
    parser.add_argument("--lll-scale", type=int, default=10**30)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.max_tail_ideals < 0:
        parser.error("--max-tail-ideals must be nonnegative")
    if args.max_samples_per_tail < 1 or args.max_samples_per_known_ideal < 1:
        parser.error("sample counts must be positive")
    if args.max_tail_factor_bits < 2 or args.max_known_factor_bits < 2:
        parser.error("factor bit bounds must exceed one")
    if args.lattice_combination_bound < 1 or args.lll_scale < 1:
        parser.error("lattice bounds must be positive")
    if args.output.exists() and not args.overwrite:
        raise FileExistsError(args.output)

    started = time.monotonic()
    source = json.loads(args.input.read_text())
    if source.get("schema") != INPUT_SCHEMA or source.get("status") != INPUT_STATUS:
        raise ValueError("unexpected residual-ideal input")
    ledger_path = Path(source["input"]["path"])
    if digest(ledger_path) != source["input"]["sha256"]:
        raise ArithmeticError("the source Minkowski ledger hash changed")
    ledger = json.loads(ledger_path.read_text())
    if ledger.get("schema") != LEDGER_SCHEMA:
        raise ValueError("unexpected source relation-ledger schema")

    coefficients = [ZZ(value) for value in ledger["defining_polynomial_ascending"]]
    ring = PolynomialRing(QQ, "x")
    x = ring.gen()
    polynomial = sum(value * x**index for index, value in enumerate(coefficients))
    declared_s = [ZZ(value) for value in ledger["selmer_rational_primes"]]
    pari.addprimes(declared_s)
    field = NumberField(polynomial, "theta")
    nf = field.pari_nf()

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
            raise ArithmeticError("a factor-base prime ideal no longer reconstructs")
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

    signature = tuple(map(int, field.signature()))
    if signature == (1, 1):
        directions = ([1, 1], [10, 1], [100, 1], [1000, 1])
    elif signature == (3, 0):
        directions = (
            [1, 1, 1],
            [10, 1, 1],
            [100, 1, 1],
            [1, 10, 10],
            [1, 100, 100],
            [10, 10, 1],
        )
    else:
        raise ArithmeticError("unexpected cubic signature")

    real_field = RealField(256)
    complex_field = ComplexField(256)
    roots = list(
        polynomial.change_ring(complex_field).roots(multiplicities=False)
    )
    real_roots = sorted(root.real() for root in roots if root.imag() == 0)
    complex_root = next((root for root in roots if root.imag() > 0), None)
    sqrt_two = real_field(2).sqrt()
    if signature == (1, 1) and (len(real_roots) != 1 or complex_root is None):
        raise ArithmeticError("the complex cubic embeddings did not reconstruct")
    if signature == (3, 0) and len(real_roots) != 3:
        raise ArithmeticError("the real cubic embeddings did not reconstruct")

    def evaluate_at_root(element, root):
        return sum(
            complex_field(coefficient) * root**power
            for power, coefficient in enumerate(element.list())
        )

    def minkowski_coordinates(element):
        if signature == (3, 0):
            return [
                real_field(evaluate_at_root(element, root).real())
                for root in real_roots
            ]
        value = evaluate_at_root(element, complex_root)
        return [
            real_field(evaluate_at_root(element, real_roots[0]).real()),
            sqrt_two * real_field(value.real()),
            sqrt_two * real_field(value.imag()),
        ]

    primes_above_cache = {}

    def primes_above(rational_prime):
        rational_prime = ZZ(rational_prime)
        if rational_prime not in primes_above_cache:
            primes_above_cache[rational_prime] = tuple(
                field.primes_above(rational_prime)
            )
        return primes_above_cache[rational_prime]

    def verify_equivalent_reduction(record, original_key, reduced_key, multiplier_key):
        original = field.ideal(pari(record[original_key]))
        reduced = field.ideal(pari(record[reduced_key]))
        multiplier = pari(record[multiplier_key])
        rebuilt = pari.idealmul(nf, reduced.pari_hnf(), multiplier)
        if str(pari.idealhnf(nf, rebuilt)) != str(original.pari_hnf()):
            raise ArithmeticError("a cached ideal-reduction identity failed")
        return reduced

    def directional_candidates(ideal, maximum_samples, excluded_keys=()):
        candidates = {}
        failures = []
        for direction in directions:
            try:
                beta = field(pari.idealmin(nf, ideal.pari_hnf(), direction))
                quotient = field.ideal(beta) / ideal
                if not quotient.is_integral() or quotient * ideal != field.ideal(beta):
                    raise ArithmeticError("directional ideal minimum does not replay")
                key = projective_coordinate_key(beta)
                if key in excluded_keys:
                    continue
                current = candidates.get(key)
                candidate = (ZZ(quotient.norm()), beta, quotient, list(direction))
                if current is None or candidate[0] < current[0]:
                    candidates[key] = candidate
            except Exception as error:
                failures.append({"direction": list(direction), "error": str(error)})
        ordered = sorted(
            candidates.values(),
            key=lambda item: (item[0].nbits(), item[0], str(item[1])),
        )[:maximum_samples]
        return ordered, failures

    def lattice_candidates(ideal, maximum_samples, excluded_keys=()):
        basis = tuple(ideal.basis())
        embedding = matrix(
            real_field, [minkowski_coordinates(element) for element in basis]
        )
        if signature == (1, 1):
            twists = ((-5, 0), (0, 0), (5, 0))
        else:
            twists = ((-3, 0), (0, -3), (0, 0), (0, 3), (3, 0))
        bound = args.lattice_combination_bound
        combinations = [
            coefficients
            for coefficients in product(range(-bound, bound + 1), repeat=3)
            if any(coefficients)
        ]
        candidates = {}
        failures = []
        for u, v in twists:
            if signature == (3, 0):
                weights = [
                    real_field(10) ** u,
                    real_field(10) ** v,
                    real_field(10) ** (-u - v),
                ]
            else:
                exponent = u + v
                weights = [
                    real_field(10) ** (-2 * exponent),
                    real_field(10) ** exponent,
                    real_field(10) ** exponent,
                ]
            lattice = matrix(
                ZZ,
                [
                    [
                        ZZ(
                            (
                                args.lll_scale
                                * weights[column]
                                * embedding[row, column]
                            ).round()
                        )
                        for column in range(3)
                    ]
                    for row in range(3)
                ],
            )
            try:
                _reduced, transform = lattice.LLL(transformation=True)
                for combination in combinations:
                    original = vector(ZZ, combination) * transform
                    beta = sum(
                        (
                            original[index] * basis[index]
                            for index in range(3)
                        ),
                        field(0),
                    )
                    if not beta:
                        continue
                    quotient = field.ideal(beta) / ideal
                    if not quotient.is_integral() or quotient * ideal != field.ideal(beta):
                        raise ArithmeticError("LLL ideal sample does not replay")
                    key = projective_coordinate_key(beta)
                    if key in excluded_keys:
                        continue
                    norm = ZZ(quotient.norm())
                    candidate = (norm, beta, quotient, [u, v])
                    current = candidates.get(key)
                    if current is None or candidate[0] < current[0]:
                        candidates[key] = candidate
            except Exception as error:
                failures.append({"twist": [u, v], "error": str(error)})
        ordered = sorted(
            candidates.values(),
            key=lambda item: (item[0].nbits(), item[0], str(item[1])),
        )[:maximum_samples]
        return ordered, failures

    candidate_function = (
        directional_candidates if args.candidate_engine == "idealmin" else lattice_candidates
    )

    def factor_quotient(ideal, beta, quotient, maximum_bits):
        norm = ZZ(quotient.norm())
        if norm <= 0:
            raise ArithmeticError("an integral quotient ideal has nonpositive norm")
        if norm.nbits() > maximum_bits:
            return None
        rational_factorization = [(ZZ(p), int(e)) for p, e in factor(norm)]
        if any(not p.is_prime() for p, _ in rational_factorization):
            raise ArithmeticError("an alleged quotient-norm factor is not prime")
        row = 0
        vertices = []
        ideal_factors = []
        reconstructed = field.ideal(1)
        reconstructed_norm_exponents = defaultdict(int)
        for rational_prime, norm_exponent in rational_factorization:
            for prime in primes_above(rational_prime):
                valuation = int(quotient.valuation(prime))
                if valuation < 0:
                    raise ArithmeticError("an integral quotient has negative valuation")
                if not valuation:
                    continue
                reconstructed *= prime**valuation
                reconstructed_norm_exponents[rational_prime] += (
                    valuation * int(prime.residue_class_degree())
                )
                hnf = str(prime.pari_hnf())
                column = factor_base_by_hnf.get(hnf)
                if valuation & 1:
                    if column is None:
                        vertices.append(("prime_ideal", hnf))
                    else:
                        row ^= 1 << column
                ideal_factors.append(
                    {
                        "rational_prime": str(rational_prime),
                        "prime_ideal_hnf": hnf,
                        "residue_degree": int(prime.residue_class_degree()),
                        "valuation": valuation,
                        "factor_base_column": column,
                    }
                )
            if reconstructed_norm_exponents[rational_prime] != norm_exponent:
                raise ArithmeticError("prime-ideal valuations do not recover the norm")
        if reconstructed != quotient or quotient * ideal != field.ideal(beta):
            raise ArithmeticError("the targeted quotient factorization does not replay")
        return {
            "beta_power_basis": [str(QQ(value)) for value in beta.list()],
            "beta_norm": str(beta.norm()),
            "quotient_ideal_hnf": str(quotient.pari_hnf()),
            "quotient_norm": str(norm),
            "quotient_norm_bit_length": int(norm.nbits()),
            "rational_norm_factorization": [
                [str(p), exponent] for p, exponent in rational_factorization
            ],
            "prime_ideal_factorization": ideal_factors,
            "fb_parity_mask_hex": hex(row),
            "outside_vertices": [list(vertex) for vertex in parity_vertices(vertices)],
            "verified_identity": "(beta) = source_ideal * quotient_ideal",
        }

    sparse = SparseLargePrimeEliminator()
    edge_references = {}
    cycles = []
    next_edge_id = 0
    outside_hnf_to_rational_prime = {}

    def add_edge(vertices, row, reference):
        nonlocal next_edge_id
        edge_id = next_edge_id
        next_edge_id += 1
        edge_references[edge_id] = reference
        vertices = parity_vertices(vertices)
        cycle, provenance = sparse.add(vertices, row, edge_id)
        if cycle is not None:
            gained = insert_row(quotient_pivots, cycle)
            cycles.append(
                {
                    "fb_parity_mask_hex": hex(cycle),
                    "edge_references": [edge_references[index] for index in sorted(provenance)],
                    "increased_quotient_rank": gained,
                }
            )

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
    if cycles:
        raise ArithmeticError("the source claimed no cycles but its edge system has one")

    known_target_edges = []
    known_failures = []
    for point_index, record in enumerate(source["killed_half_ideal_reductions"]):
        ideal = verify_equivalent_reduction(
            record,
            "half_ideal_hnf",
            "reduced_ideal_hnf",
            "reduction_multiplier",
        )
        candidates, failures = candidate_function(
            ideal, args.max_samples_per_known_ideal, ()
        )
        if failures:
            known_failures.append({"label": record["label"], "failures": failures})
        for sample_index, (_norm, beta, quotient, direction) in enumerate(candidates):
            factored = factor_quotient(
                ideal, beta, quotient, args.max_known_factor_bits
            )
            if factored is None:
                continue
            edge_record = {
                "label": record["label"],
                "source_ideal_role": "certified_MW29_half_ideal_killed_in_relative_quotient",
                "source_reduced_ideal_hnf": str(ideal.pari_hnf()),
                "direction": direction,
                **factored,
            }
            edge_index = len(known_target_edges)
            known_target_edges.append(edge_record)
            for ideal_factor in factored["prime_ideal_factorization"]:
                outside_hnf_to_rational_prime[
                    ideal_factor["prime_ideal_hnf"]
                ] = int(ideal_factor["rational_prime"])
            add_edge(
                [tuple(vertex) for vertex in factored["outside_vertices"]],
                int(factored["fb_parity_mask_hex"], 16),
                {"kind": "known_MW29_target_edge", "index": edge_index},
            )

    opaque_records = [
        (index, edge)
        for index, edge in enumerate(source["partial_edges"])
        if edge["tail_classification"] == "opaque_residual_ideal"
    ]
    if args.max_tail_ideals:
        opaque_records = opaque_records[: args.max_tail_ideals]
    tail_target_edges = []
    tail_failures = []
    skipped_tail_factor_bound = 0
    for target_index, (source_edge_index, record) in enumerate(opaque_records, 1):
        ideal = verify_equivalent_reduction(
            record,
            "residual_ideal_hnf",
            "reduced_tail_hnf",
            "tail_reduction_multiplier",
        )
        generator = ledger["generators"][int(record["generator_index"])]
        alpha = sum(
            (
                QQ(value) * field.gen() ** index
                for index, value in enumerate(generator["power_basis"])
            ),
            field(0),
        )
        multiplier = field(pari(record["tail_reduction_multiplier"]))
        replay_beta = alpha / multiplier
        replay_quotient = field.ideal(replay_beta) / ideal
        if not replay_quotient.is_integral():
            raise ArithmeticError("the source relation does not transport to the reduced tail")
        excluded_keys = {projective_coordinate_key(replay_beta)}
        candidates, failures = candidate_function(
            ideal, args.max_samples_per_tail, excluded_keys
        )
        if failures:
            tail_failures.append(
                {"source_edge_index": source_edge_index, "failures": failures}
            )
        for sample_index, (_norm, beta, quotient, direction) in enumerate(candidates):
            factored = factor_quotient(
                ideal, beta, quotient, args.max_tail_factor_bits
            )
            if factored is None:
                skipped_tail_factor_bound += 1
                continue
            edge_record = {
                "source_edge_index": source_edge_index,
                "source_ideal_role": "opaque_residual_ideal_vertex",
                "source_reduced_ideal_hnf": str(ideal.pari_hnf()),
                "direction": direction,
                **factored,
            }
            edge_index = len(tail_target_edges)
            tail_target_edges.append(edge_record)
            for ideal_factor in factored["prime_ideal_factorization"]:
                outside_hnf_to_rational_prime[
                    ideal_factor["prime_ideal_hnf"]
                ] = int(ideal_factor["rational_prime"])
            add_edge(
                [("residual_ideal", str(ideal.pari_hnf()))]
                + [tuple(vertex) for vertex in factored["outside_vertices"]],
                int(factored["fb_parity_mask_hex"], 16),
                {"kind": "residual_ideal_target_edge", "index": edge_index},
            )
        if target_index % 500 == 0:
            print(
                f"{PROTOCOL}|curve={source['curve_id']}|targeted={target_index}/{len(opaque_records)}"
                f"|target_edges={len(tail_target_edges)}|cycles={len(cycles)}"
                f"|rank_gain={len(quotient_pivots)-baseline_rank}",
                flush=True,
            )

    used_outside_hnfs = {
        vertex[1]
        for edge in source["partial_edges"]
        for vertex in edge["large_prime_vertices"]
        if vertex[0] == "prime_ideal"
    }
    used_outside_hnfs.update(
        vertex[1]
        for edge in known_target_edges + tail_target_edges
        for vertex in edge["outside_vertices"]
    )
    missing_prime_labels = sorted(
        hnf for hnf in used_outside_hnfs if hnf not in outside_hnf_to_rational_prime
    )
    if missing_prime_labels:
        raise ArithmeticError("an outside prime-ideal vertex lacks its rational prime")
    outside_rational_primes = sorted(
        {outside_hnf_to_rational_prime[hnf] for hnf in used_outside_hnfs}
    )
    canonical_outside_edges = []
    for rational_prime in outside_rational_primes:
        row = 0
        vertices = []
        rational_ideal = field.ideal(rational_prime)
        factors = []
        for prime in primes_above(rational_prime):
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
        edge_index = len(canonical_outside_edges)
        canonical_outside_edges.append(
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
        "status": "BOUNDED_EXACT_MW29_RELATIVE_IDEAL_RELATIONS_NOT_CLASS_GROUP_COMPLETION",
        "curve_id": source["curve_id"],
        "input": {"path": str(args.input.resolve()), "sha256": digest(args.input)},
        "source_minkowski_ledger": {
            "path": str(ledger_path.resolve()),
            "sha256": digest(ledger_path),
        },
        "settings": {
            "max_tail_ideals": args.max_tail_ideals,
            "max_samples_per_tail": args.max_samples_per_tail,
            "max_tail_factor_bits": args.max_tail_factor_bits,
            "max_samples_per_known_ideal": args.max_samples_per_known_ideal,
            "max_known_factor_bits": args.max_known_factor_bits,
            "candidate_engine": args.candidate_engine,
            "lattice_combination_bound": args.lattice_combination_bound,
            "lll_scale": args.lll_scale,
            "directions": [list(direction) for direction in directions],
        },
        "factor_base_width": len(factor_base),
        "baseline_rank_after_canonical_and_S_rows": baseline_rank,
        "baseline_bounded_quotient_dimension": len(factor_base) - baseline_rank,
        "source_edge_count": len(source["partial_edges"]),
        "known_MW29_ideal_count": len(source["killed_half_ideal_reductions"]),
        "known_MW29_target_edge_count": len(known_target_edges),
        "opaque_tail_ideal_count_available": sum(
            edge["tail_classification"] == "opaque_residual_ideal"
            for edge in source["partial_edges"]
        ),
        "opaque_tail_ideal_count_targeted": len(opaque_records),
        "tail_target_edge_count": len(tail_target_edges),
        "skipped_tail_sample_above_factor_bound": skipped_tail_factor_bound,
        "directional_failures": {
            "known_MW29": known_failures,
            "opaque_tails": tail_failures,
        },
        "known_MW29_target_edges": known_target_edges,
        "tail_target_edges": tail_target_edges,
        "canonical_outside_rational_prime_edges": canonical_outside_edges,
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
            "Every target edge verifies (beta) = source_ideal * quotient_ideal and all prime-ideal valuations exactly.",
            "Certified MW29 half-ideal sources are killed before sparse elimination; opaque residual sources remain explicit vertices.",
            "Only exact sparse-incidence dependencies create factor-base relations.",
            "The factor base has no generation certificate, so the displayed quotient dimension is not an S-class or Selmer upper bound.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    temporary.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    temporary.replace(args.output)
    print(
        f"{PROTOCOL}|curve={source['curve_id']}|known_edges={len(known_target_edges)}"
        f"|tail_edges={len(tail_target_edges)}|cycles={len(cycles)}"
        f"|rank_gain={len(quotient_pivots)-baseline_rank}"
        f"|dimension={len(factor_base)-len(quotient_pivots)}"
        f"|seconds={output['elapsed_seconds']:.3f}|output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
