#!/usr/bin/env sage-python
"""Exhaust canonical-height shells on the promoted rank-nine paired base.

status: ACTIVE_SEARCH
claim: bounded exhaustive split-bisection sieve with exact specialization certificates
inputs: promoted-base lattice certificate, published R17 sections, complete bisection batch
outputs: artifacts/generated-results/elkies-2026-rank9-paired-base-sieve.json
supersedes: none; exhausts shells only partially sampled by search_elkies_2026_rank9_paired_base.sage

For each distinct finite parameter in the declared height shell, a modular
bitset sieve reduces the 39,118 non-defining square tests to exact integer
square tests.  The sieve has no false negatives: a rational square remains a
quadratic residue (or zero) at every usable prime.  At the validation height,
all nineteen known points are specialized and their independence is certified
by an exact finite-quotient witness.  Larger discovery shells may use
``--certify split``; every split hit is then fully materialized and certified,
while no rank claim is made for no-hit fibres.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import importlib.util
import json
from pathlib import Path
import shlex
import sys
from tempfile import NamedTemporaryFile
from time import perf_counter

from sage.all import EllipticCurve, QQ, RealField, ZZ, matrix


ROOT = Path(__file__).resolve().parents[2]
SEARCH_SCRIPT = ROOT / "elkies-k3/scripts/search_elkies_2026_rank9_paired_base.sage"
LATTICE = ROOT / "artifacts/generated-results/elkies-2026-rank9-paired-base-search.json"
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
BISECTIONS = ROOT / "artifacts/generated-results/elkies-2026-equation-bisections-full.json"
FINITE_HELPER = ROOT / "elliptic-curves/cas/elliptic_candidate_record.py"
OUTPUT = ROOT / "artifacts/generated-results/elkies-2026-rank9-paired-base-sieve.json"
TARGET_MASKS = (42110, 43109)
# These primes were selected from the exact denominator profile of the full
# batch.  None divides a coefficient denominator, so every table contributes
# a genuine quadratic-residue filter.
DEFAULT_PRIMES = (131, 137, 151, 157, 163, 167, 181, 191, 197, 199, 211, 227)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_loader(name, SourceFileLoader(name, str(path)))
    if specification is None or specification.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


def atomic_write(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile("w", dir=path.parent, prefix=path.name + ".", delete=False) as stream:
        temporary = Path(stream.name)
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
    temporary.replace(path)


def normalized_homogeneous_quadratic(coefficients):
    values = tuple(QQ(value) for value in coefficients)
    denominator = ZZ(1)
    for value in values:
        denominator = denominator.lcm(value.denominator())
    integers = tuple(ZZ(value * denominator) for value in values)
    return denominator, integers


def exact_square_value(denominator, coefficients, numerator, parameter_denominator):
    c0, c1, c2 = coefficients
    homogeneous = (
        c0 * parameter_denominator**2
        + c1 * numerator * parameter_denominator
        + c2 * numerator**2
    )
    square_integer = denominator * homogeneous
    if square_integer < 0 or not square_integer.is_square():
        return None
    square_root = QQ(square_integer.sqrt()) / (denominator * parameter_denominator)
    value = QQ(homogeneous) / (denominator * parameter_denominator**2)
    if square_root**2 != value:
        raise ArithmeticError("integer square reconstruction failed")
    return value, square_root


def build_residue_bitsets(quadratics, primes):
    all_bits = (1 << len(quadratics)) - 1
    tables = {}
    for prime in primes:
        residues = [0] * (prime + 1)
        for index, (denominator, coefficients) in enumerate(quadratics):
            bit = 1 << index
            if denominator % prime == 0:
                for projective_value in range(prime + 1):
                    residues[projective_value] |= bit
                continue
            c0, c1, c2 = (int(value % prime) for value in coefficients)
            factor = int(denominator % prime)
            for value in range(prime):
                evaluation = factor * (c0 + c1 * value + c2 * value * value) % prime
                if evaluation == 0 or pow(evaluation, (prime - 1) // 2, prime) == 1:
                    residues[value] |= bit
            evaluation = factor * c2 % prime
            if evaluation == 0 or pow(evaluation, (prime - 1) // 2, prime) == 1:
                residues[prime] |= bit
        if any(value & ~all_bits for value in residues):
            raise ArithmeticError("residue bitset exceeded the bisection universe")
        tables[prime] = residues
    return tables, all_bits


def sieve_indices(numerator, denominator, tables, all_bits):
    survivors = all_bits
    for prime, residues in tables.items():
        if denominator % prime:
            value = (numerator % prime) * pow(denominator % prime, -1, prime) % prime
        else:
            value = prime
        survivors &= residues[value]
        if survivors == 0:
            break
    while survivors:
        bit = survivors & -survivors
        yield bit.bit_length() - 1
        survivors ^= bit


def compact_certificate(search, helper, model, points, prime_bound):
    certificate = search.best_finite_certificate(helper, model, points, prime_bound)
    successful = certificate["independent_subset_certificate"]
    if successful is None:
        return {
            "certified_rank_lower_bound": 0,
            "best_relation_prime": certificate["best_relation_prime"],
            "independent_point_indices_zero_based": [],
            "successful_certificate": None,
        }
    return {
        "certified_rank_lower_bound": certificate["certified_rank_lower_bound"],
        "best_relation_prime": certificate["best_relation_prime"],
        "independent_point_indices_zero_based": certificate[
            "independent_point_indices_zero_based"
        ],
        "successful_certificate": successful,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lattice", type=Path, default=LATTICE)
    parser.add_argument("--model", type=Path, default=MODEL)
    parser.add_argument("--sections", type=Path, default=SECTIONS)
    parser.add_argument("--bisections", type=Path, default=BISECTIONS)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--height-bound", type=int, default=60)
    parser.add_argument("--height-floor", type=int, default=0)
    parser.add_argument(
        "--gram-scale",
        type=int,
        default=100000000,
        help="integer scale for qfminim; lower values can be needed for large shells",
    )
    parser.add_argument("--primes", default=",".join(map(str, DEFAULT_PRIMES)))
    parser.add_argument("--certificate-prime-bound", type=int, default=300)
    parser.add_argument("--certify", choices=("all", "split"), default="all")
    parser.add_argument(
        "--store",
        choices=("all", "hits"),
        default=None,
        help="store every fibre, or only split/collision hits plus a digest of the full shell",
    )
    parser.add_argument("--checkpoint-every", type=int, default=100)
    args = parser.parse_args()
    if not 0 <= args.height_floor < args.height_bound:
        parser.error("require 0 <= height-floor < height-bound")
    primes = tuple(int(value) for value in args.primes.split(",") if value)
    if any(value < 3 or not ZZ(value).is_prime() for value in primes):
        parser.error("all sieve moduli must be odd primes")
    store_policy = args.store or ("all" if args.certify == "all" else "hits")
    if args.certify == "all" and store_policy != "all":
        parser.error("--certify all requires --store all so every certificate is retained")

    started = perf_counter()
    search = load_module("elkies_rank9_search_helpers", SEARCH_SCRIPT)
    finite_helper = load_module("elkies_rank9_finite_helper", FINITE_HELPER)
    lattice_document = json.loads(args.lattice.read_text())
    model_document = json.loads(args.model.read_text())
    section_document = json.loads(args.sections.read_text())
    bisection_document = json.loads(args.bisections.read_text())
    lattice = lattice_document["base_lattice"]
    field = RealField(420)
    lll_gram = matrix(field, [[field(value) for value in row] for row in lattice["lll_reduced_canonical_height_gram"]])
    if args.gram_scale < 1000000:
        parser.error("gram-scale below 1e6 is not supported by the rounding guard")
    scale = ZZ(args.gram_scale)
    integral_gram = matrix(ZZ, [[ZZ((scale * value).round()) for value in row] for row in lll_gram.rows()])
    third_curve = EllipticCurve(QQ, lattice["third_quotient_model_a1_a2_a3_a4_a6"])
    paired_curve = EllipticCurve(QQ, list(search.PAIRED_MODEL))
    lll_basis = [third_curve(QQ(point[0]), QQ(point[1])) for point in lattice["lll_basis"]]
    to_paired, to_cover = search.paired_cover_maps(third_curve, paired_curve)
    enumerated, pareto, enumeration_summary = search.enumerate_base_points(
        third_curve,
        lll_basis,
        lll_gram,
        integral_gram,
        args.height_bound,
        scale,
        to_paired,
        to_cover,
    )
    shell = [record for record in enumerated if field(record["canonical_height"]) > args.height_floor]
    shell.sort(key=lambda record: (field(record["canonical_height"]), record["t"]))

    all_records = bisection_document["bisections"]
    by_mask = {int(record["lattice_orbit_mask"]): record for record in all_records}
    other_records = [
        record for record in all_records if int(record["lattice_orbit_mask"]) not in TARGET_MASKS
    ]
    quadratics = [
        normalized_homogeneous_quadratic(record["residual_chord"]["q_coefficients"])
        for record in other_records
    ]
    tables, all_bits = build_residue_bitsets(quadratics, primes)

    a_coefficients = tuple(QQ(value) for value in model_document["A_coefficients_low_to_high"])
    b_coefficients = tuple(QQ(value) for value in model_document["B_coefficients_low_to_high"])
    results = []
    shell_hasher = sha256()
    modular_survivor_tests = 0
    exact_square_tests = 0
    zero_count = 0
    split_fibre_count = 0
    total_extra_splits = 0
    for ordinal, candidate in enumerate(shell, start=1):
        t_value = QQ(candidate["t"])
        numerator = ZZ(t_value.numerator())
        denominator = ZZ(t_value.denominator())
        survivor_indices = tuple(sieve_indices(numerator, denominator, tables, all_bits))
        modular_survivor_tests += len(survivor_indices)
        split = []
        zeros = []
        for index in survivor_indices:
            exact_square_tests += 1
            exact = exact_square_value(*quadratics[index], numerator, denominator)
            if exact is None:
                continue
            q_value, square_root = exact
            record = other_records[index]
            mask = int(record["lattice_orbit_mask"])
            entry = {
                "lattice_orbit_mask": mask,
                "orbit_hex": f"0x{mask:05x}",
                "q_value": search.rational_text(q_value),
                "square_root": search.rational_text(square_root),
            }
            if q_value == 0:
                zeros.append(entry)
                zero_count += 1
            else:
                split.append((entry, record, square_root))

        should_materialize = args.certify == "all" or bool(split)
        row = {
            **candidate,
            "shell_ordinal": ordinal,
            "modular_survivor_count": len(survivor_indices),
            "exact_square_test_count": len(survivor_indices),
            "other_split_bisection_count": len(split),
            "branch_collision_count": len(zeros),
            "branch_collisions": zeros,
            "certification_policy": args.certify,
        }
        split_fibre_count += bool(split)
        total_extra_splits += len(split)
        shell_hasher.update(
            (json.dumps(
                {
                    "t": candidate["t"],
                    "lll_coefficients": candidate["lll_coefficients"],
                    "canonical_height": candidate["canonical_height"],
                    "split_masks": [entry[0]["lattice_orbit_mask"] for entry in split],
                    "collision_masks": [entry["lattice_orbit_mask"] for entry in zeros],
                },
                sort_keys=True,
                separators=(",", ":"),
            ) + "\n").encode()
        )
        if should_materialize:
            coefficient_a = search.evaluate(a_coefficients, t_value)
            coefficient_b = search.evaluate(b_coefficients, t_value)
            model = tuple(search.to_fraction(value) for value in (0, 0, 0, coefficient_a, coefficient_b))
            row["specialized_model_a1_a2_a3_a4_a6"] = [
                search.rational_text(value) for value in (0, 0, 0, coefficient_a, coefficient_b)
            ]
            curve = EllipticCurve(QQ, [coefficient_a, coefficient_b])
            generic_points = search.reconstruct_sections(section_document, t_value, curve)
            selected_points = []
            for mask in TARGET_MASKS:
                q_value = search.evaluate(by_mask[mask]["residual_chord"]["q_coefficients"], t_value)
                if q_value == 0 or not q_value.is_square():
                    raise ArithmeticError("a defining paired cover stopped splitting")
                selected_points.append(search.lifted_point(by_mask[mask], t_value, q_value.sqrt(), curve))
            known_points = generic_points + selected_points
            known_certificate = compact_certificate(
                search, finite_helper, model, known_points, args.certificate_prime_bound
            )
            extra_points = []
            seen = {(point[0], point[1]) for point in known_points}
            materialized_split = []
            for entry, record, square_root in split:
                point = search.lifted_point(record, t_value, square_root, curve)
                if (point[0], point[1]) not in seen:
                    extra_points.append(point)
                    seen.add((point[0], point[1]))
                materialized_split.append({**entry, "point": search.point_text(point)})
            all_points = known_points + extra_points
            total_certificate = (
                known_certificate
                if not extra_points
                else compact_certificate(
                    search, finite_helper, model, all_points, args.certificate_prime_bound
                )
            )
            row.update(
                {
                    "known_point_count": 19,
                    "known_rank_certificate": known_certificate,
                    "split_bisections": materialized_split,
                    "materialized_distinct_extra_point_count": len(extra_points),
                    "extra_points": [search.point_text(point) for point in extra_points],
                    "total_displayed_point_count": len(all_points),
                    "total_rank_certificate": total_certificate,
                }
            )
        else:
            row.update(
                {
                    "specialized_model_a1_a2_a3_a4_a6": None,
                    "known_point_count": 19,
                    "known_rank_certificate": None,
                    "split_bisections": [entry for entry, _record, _root in split],
                    "materialized_distinct_extra_point_count": 0,
                    "extra_points": [],
                    "total_displayed_point_count": 19,
                    "total_rank_certificate": None,
                }
            )
        if store_policy == "all" or split or zeros:
            results.append(row)
        if ordinal % args.checkpoint_every == 0:
            print(
                "ELKIES2026R9PAIRSIEVE|"
                f"progress={ordinal}/{len(shell)}|split_hits={split_fibre_count}|"
                f"exact_tests={exact_square_tests}",
                flush=True,
            )

    rank_counts = {}
    for row in results:
        certificate = row["total_rank_certificate"]
        if certificate is not None:
            rank = str(certificate["certified_rank_lower_bound"])
            rank_counts[rank] = rank_counts.get(rank, 0) + 1
    result = {
        "schema": "elkies-k3.elkies-2026-rank9-paired-base-sieve.v1",
        "status": "PASS_BOUNDED_EXHAUSTIVE_SPLIT_BISECTION_SIEVE",
        "inputs": {
            display_path(Path(__file__).resolve()): digest(Path(__file__).resolve()),
            display_path(SEARCH_SCRIPT): digest(SEARCH_SCRIPT),
            display_path(args.lattice): digest(args.lattice),
            display_path(args.model): digest(args.model),
            display_path(args.sections): digest(args.sections),
            display_path(args.bisections): digest(args.bisections),
            display_path(FINITE_HELPER): digest(FINITE_HELPER),
        },
        "bounds": {
            "canonical_height_floor_exclusive": args.height_floor,
            "canonical_height_bound_inclusive": args.height_bound,
            "integral_gram_scale": int(scale),
            "certification_policy": args.certify,
            "storage_policy": store_policy,
            "certificate_prime_bound": args.certificate_prime_bound,
            "sieve_primes": list(primes),
        },
        "enumeration": {
            **enumeration_summary,
            "shell_distinct_finite_t_values": len(shell),
            "pareto_frontier": pareto,
        },
        "sieve": {
            "nondefining_quadratic_count": len(other_records),
            "naive_exact_square_test_count": len(shell) * len(other_records),
            "modular_survivor_exact_test_count": modular_survivor_tests,
            "exact_square_test_count": exact_square_tests,
            "branch_collision_count": zero_count,
            "fibres_with_extra_splits": split_fibre_count,
            "total_extra_split_bisections": total_extra_splits,
            "full_shell_ledger_sha256": shell_hasher.hexdigest(),
            "stored_fibre_count": len(results),
            "no_false_negative_argument": (
                "For q(t)=M/(D*b^2), rational squareness is equivalent to D*M being an integer "
                "square. Its homogeneous reduction is therefore a square or zero at every sieve prime."
            ),
        },
        "certified_rank_lower_bound_counts": dict(sorted(rank_counts.items(), key=lambda item: int(item[0]))),
        "fibres": results,
        "runtime_seconds": perf_counter() - started,
        "reproducing_command": shlex.join(sys.argv),
        "proof_boundary": (
            "Enumeration and square tests are exhaustive only inside the declared canonical-height shell. "
            "Every nonzero split is exact. A rank lower bound is asserted only where a successful exact "
            "finite-quotient certificate is stored. No rank upper bound or completeness beyond the shell is inferred."
        ),
    }
    atomic_write(args.output, result)
    print(
        "ELKIES2026R9PAIRSIEVE|"
        f"shell_t={len(shell)}|extra_split_fibres={result['sieve']['fibres_with_extra_splits']}|"
        f"exact_tests={exact_square_tests}|status={result['status']}|output={display_path(args.output)}"
    )


if __name__ == "__main__":
    main()
