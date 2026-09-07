#!/usr/bin/env sage-python
"""Invert rational intersection-one V4 products in a native priority prefix.

The exact 1,024-prefix contains 10,362 pairs of rational bisections with
intersection number one; the complete 39,147-class atlas contains 4,358,409.
Their fibre products are genus-one V4 bases with a QQ-point, and their product
characters are explicit quartics.  This script tests all selected product
characters against every one of the 63,917 minimum-norm-eight genus-one
bisection pencils.

The first prime scans every projective pencil parameter and uses a hash table
of all target quartics.  Later primes revisit only surviving (trace,target)
pairs.  A pair is rejected only by projective coefficient inequality at a
good reduction; all survivors are resolved over QQ, including the rational
constant squareclass.  A hit supplies the third nontrivial V4 character and
therefore the exact rank decomposition 17+1+1+1.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
import runpy

import numpy as np
from sage.all import GF, PolynomialRing, QQ, ZZ, lcm, matrix, vector
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = Path(__file__).resolve().parent
BISECTIONS = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-11952-alternate-bisections-cheapest-1024-v1.json"
)
DIRECT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
)
NORM8_TABLE = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.tsv"
)
INVERSION_SCRIPT = SCRIPTS / "search_r17_norm12_11952_product_bisection_inversion.sage"
SHORTLIST_SCRIPT = SCRIPTS / "select_r17_norm12_11952_v4_pair_shortlist.sage"
CHORD_SCRIPT = SCRIPTS / "construct_elkies_2026_bisections.sage"
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-11952-v4-intersection-one-product-inversion-v1.json"
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def rational_text(value) -> str:
    value = QQ(value)
    return str(value.numerator()) if value.denominator() == 1 else f"{value.numerator()}/{value.denominator()}"


def coefficient_text(polynomial) -> list[str]:
    return [rational_text(polynomial[index]) for index in range(polynomial.degree() + 1)]


def primitive_coefficients(polynomial, length=5) -> tuple[int, ...]:
    values = [QQ(polynomial[index]) for index in range(length)]
    denominator = lcm([value.denominator() for value in values])
    integers = [ZZ(denominator * value) for value in values]
    content = ZZ(0)
    for value in integers:
        content = gcd(int(content), abs(int(value)))
    if content:
        integers = [value // content for value in integers]
    pivot = next(value for value in integers if value)
    if pivot < 0:
        integers = [-value for value in integers]
    return tuple(map(int, integers))


def projective_key(values):
    values = tuple(values)
    pivot = next((value for value in values if value), None)
    if pivot is None:
        return None
    inverse = pivot**-1
    return tuple(value * inverse for value in values)


def target_key(coefficients, field, inverted=False):
    values = [field(value) for value in coefficients]
    if inverted:
        values.reverse()
    return projective_key(values)


def pencil_keys(context, word, inverse):
    trace = inverse["trace_from_word"](context, word, context["multiples"])
    if trace.is_zero():
        raise ArithmeticError("norm-eight trace reduced to zero")
    pencil = inverse["trace_pencil_family"](context, trace)
    family = pencil["q_lambda_family"]
    field = context["coefficient_field"]
    keys = set()
    zero_vector = False
    for parameter in field:
        q = family[4]
        for coefficient in reversed(family[:4]):
            q = q * parameter + coefficient
        key = projective_key(tuple(q[index] for index in range(5)))
        if key is None:
            zero_vector = True
        else:
            keys.add(key)
    infinity_key = projective_key(tuple(family[4][index] for index in range(5)))
    if infinity_key is None:
        zero_vector = True
    else:
        keys.add(infinity_key)
    return pencil["chart"], keys, zero_vector


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bisections", type=Path, default=BISECTIONS)
    parser.add_argument("--direct", type=Path, default=DIRECT)
    parser.add_argument("--norm8-table", type=Path, default=NORM8_TABLE)
    parser.add_argument("--pool-size", type=int, default=1024)
    parser.add_argument("--primes", default="131,137,151")
    parser.add_argument("--progress-every", type=int, default=1000)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    primes = [int(value) for value in args.primes.split(",") if value]
    if args.pool_size < 2 or not primes or any(not ZZ(prime).is_prime() for prime in primes):
        parser.error("require pool-size >=2 and a nonempty prime list")

    bisection_path = args.bisections.resolve()
    direct_path = args.direct.resolve()
    norm8_path = args.norm8_table.resolve()
    batch = json.loads(bisection_path.read_text())
    direct = json.loads(direct_path.read_text())
    if batch.get("schema") != "elkies-k3.bisection-extension-input.v1":
        raise ValueError("unexpected bisection schema")
    if direct.get("status") != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
        raise ValueError("unexpected direct-model status")
    source_pool_size = len(batch["bisections"])
    records = sorted(batch["bisections"], key=lambda row: int(row["priority_rank"]))
    records = records[: args.pool_size]
    if len(records) != args.pool_size or [int(row["priority_rank"]) for row in records] != list(range(1, args.pool_size + 1)):
        raise ArithmeticError("bisection pool is not a complete priority prefix")

    selector = runpy.run_path(str(SHORTLIST_SCRIPT))
    inverse = runpy.run_path(str(INVERSION_SCRIPT))
    helper = runpy.run_path(str(CHORD_SCRIPT))
    ring = PolynomialRing(QQ, "u")
    gram = matrix(ZZ, direct["frame_certificate"]["frame_gram"])
    vector_array = np.array(
        [record["direct_alternate_w"] for record in records], dtype=np.int64
    )
    gram_array = np.array(gram, dtype=np.int64)
    crude_pairing_bound = (
        17
        * int(np.max(np.abs(vector_array))) ** 2
        * int(np.max(np.abs(gram_array)))
    )
    if crude_pairing_bound >= 2**63:
        raise ArithmeticError("int64 pairing bound is unsafe")
    covector_array = vector_array @ gram_array
    branches = [selector["branch_polynomial"](record, ring) for record in records]
    norms = np.einsum("ij,ij->i", covector_array, vector_array)
    if not np.all(norms == 10):
        raise ArithmeticError("bisection pool contains a non-norm-ten class")

    target_pairs = []
    target_primitives = []
    target_manifest_hasher = sha256()
    pairing_block_size = 512
    for right_start in range(0, len(records), pairing_block_size):
        right_stop = min(len(records), right_start + pairing_block_size)
        pairings = covector_array[right_start:right_stop] @ vector_array.T
        for offset, pairing_row in enumerate(pairings):
            right_index = right_start + offset
            for left_index in np.flatnonzero(pairing_row[:right_index] == 7):
                left_index = int(left_index)
                product = ring(branches[left_index] * branches[right_index])
                if product.degree() != 4 or not product.is_squarefree():
                    raise ArithmeticError(
                        "intersection-one product is not squarefree quartic"
                    )
                pair_key = (
                    f"{records[left_index]['label']}:"
                    f"{records[right_index]['label']}"
                )
                primitive = primitive_coefficients(product)
                target_pairs.append((left_index, right_index))
                target_primitives.append(primitive)
                target_manifest_hasher.update(
                    json.dumps(
                        [
                            pair_key,
                            [
                                int(records[left_index]["priority_rank"]),
                                int(records[right_index]["priority_rank"]),
                            ],
                            primitive,
                        ],
                        separators=(",", ":"),
                    ).encode()
                )
                target_manifest_hasher.update(b"\n")
    if args.pool_size == 1024 and len(target_pairs) != 10362:
        raise ArithmeticError(f"expected 10,362 intersection-one targets, found {len(target_pairs)}")
    complete_native_atlas = args.pool_size == source_pool_size == 39147
    if complete_native_atlas and len(target_pairs) != 4358409:
        raise ArithmeticError(
            f"expected 4,358,409 complete-atlas targets, found {len(target_pairs)}"
        )
    if not target_pairs:
        raise ArithmeticError("the selected prefix has no intersection-one targets")
    target_manifest_sha256 = target_manifest_hasher.hexdigest()

    def exact_target_record(target_index):
        left_index, right_index = target_pairs[target_index]
        product = ring(branches[left_index] * branches[right_index])
        return {
            "index": int(target_index),
            "pair_key": (
                f"{records[left_index]['label']}:{records[right_index]['label']}"
            ),
            "record_indices": [left_index, right_index],
            "labels": [records[left_index]["label"], records[right_index]["label"]],
            "priority_ranks": [
                int(records[left_index]["priority_rank"]),
                int(records[right_index]["priority_rank"]),
            ],
            "product_quartic_coefficients_low_to_high": coefficient_text(product),
            "primitive_product_coefficients_low_to_high": list(
                target_primitives[target_index]
            ),
        }

    with norm8_path.open(newline="") as stream:
        norm8_rows = list(csv.DictReader(stream, delimiter="\t"))
    if len(norm8_rows) != 63917:
        raise ArithmeticError("norm-eight table is incomplete")
    words = [inverse["parse_vector"](row["section_basis_w"]) for row in norm8_rows]
    maximum_coefficient = max(abs(int(value)) for word in words for value in word)
    contexts = {
        prime: inverse["modular_context"](prime, direct, helper, maximum_coefficient)
        for prime in primes
    }

    first_prime = primes[0]
    first_field = contexts[first_prime]["coefficient_field"]
    first_lookups = {"finite": defaultdict(list), "inverted_at_infinity": defaultdict(list)}
    for target_index, coefficients in enumerate(target_primitives):
        first_lookups["finite"][target_key(coefficients, first_field)].append(target_index)
        first_lookups["inverted_at_infinity"][target_key(coefficients, first_field, True)].append(target_index)

    candidates = set()
    trace_failures = Counter()
    zero_vector_traces = []
    for trace_index, word in enumerate(words):
        try:
            chart, keys, zero_vector = pencil_keys(contexts[first_prime], word, inverse)
        except (ArithmeticError, ValueError, ZeroDivisionError) as error:
            trace_failures[(first_prime, type(error).__name__)] += 1
            candidates.update((trace_index, target_index) for target_index in range(len(target_pairs)))
            continue
        if zero_vector:
            zero_vector_traces.append({"prime": first_prime, "trace_index": trace_index})
            candidates.update((trace_index, target_index) for target_index in range(len(target_pairs)))
            continue
        lookup = first_lookups[chart]
        for key in keys:
            candidates.update((trace_index, target_index) for target_index in lookup.get(key, ()))
        if args.progress_every and (trace_index + 1) % args.progress_every == 0:
            print(
                f"ALLV4INVERSIONPROGRESS|prime={first_prime}|done={trace_index + 1}"
                f"|traces={len(words)}|candidate_pairs={len(candidates)}",
                flush=True,
            )
    prime_candidate_counts = [{"prime": first_prime, "candidate_pair_count": len(candidates)}]

    for prime in primes[1:]:
        if not candidates:
            prime_candidate_counts.append({"prime": prime, "candidate_pair_count": 0, "skipped_after_empty": True})
            continue
        field = contexts[prime]["coefficient_field"]
        target_keys = {
            chart: [
                target_key(coefficients, field, chart == "inverted_at_infinity")
                for coefficients in target_primitives
            ]
            for chart in ("finite", "inverted_at_infinity")
        }
        by_trace = defaultdict(set)
        for trace_index, target_index in candidates:
            by_trace[trace_index].add(target_index)
        retained = set()
        for trace_index, target_indices in by_trace.items():
            try:
                chart, keys, zero_vector = pencil_keys(contexts[prime], words[trace_index], inverse)
            except (ArithmeticError, ValueError, ZeroDivisionError) as error:
                trace_failures[(prime, type(error).__name__)] += 1
                retained.update((trace_index, target_index) for target_index in target_indices)
                continue
            if zero_vector:
                zero_vector_traces.append({"prime": prime, "trace_index": trace_index})
                retained.update((trace_index, target_index) for target_index in target_indices)
                continue
            keys_for_chart = target_keys[chart]
            retained.update(
                (trace_index, target_index)
                for target_index in target_indices
                if keys_for_chart[target_index] in keys
            )
        candidates = retained
        prime_candidate_counts.append({"prime": prime, "candidate_pair_count": len(candidates)})
        print(
            f"ALLV4INVERSIONPRIME|prime={prime}|candidate_pairs={len(candidates)}",
            flush=True,
        )

    exact_context = inverse["exact_model_and_basis"](direct, helper)
    candidates_by_trace = defaultdict(list)
    for trace_index, target_index in sorted(candidates):
        candidates_by_trace[trace_index].append(target_index)
    exact_resolutions = []
    hits = []
    for trace_index, target_indices in candidates_by_trace.items():
        trace = inverse["trace_from_word"](exact_context, words[trace_index])
        pencil = inverse["trace_pencil_family"](exact_context, trace)
        resolutions = []
        trace_hits = []
        for target_index in target_indices:
            target_record = exact_target_record(target_index)
            target = ring(
                [
                    QQ(value)
                    for value in target_record[
                        "product_quartic_coefficients_low_to_high"
                    ]
                ]
            )
            if pencil["chart"] == "inverted_at_infinity":
                target = inverse["reciprocal_polynomial"](target, 4, ring)
            roots, infinity_matches, common = inverse["exact_parameters"](
                pencil["q_lambda_family"], target, ring
            )
            target_result = {
                "target_index": int(target_index),
                "pair_key": target_record["pair_key"],
                "coefficient_gcd_degree": int(common.degree()),
                "rational_finite_parameters": [rational_text(root) for root in roots],
                "infinity_matches": bool(infinity_matches),
                "squareclass_hits": [],
            }
            for parameter in roots:
                q = sum(
                    (
                        parameter**power * pencil["q_lambda_family"][power]
                        for power in range(5)
                    ),
                    ring.zero(),
                )
                multiplier = inverse["squareclass_multiplier"](q, target)
                if multiplier.is_square():
                    hit = {
                        "target_index": int(target_index),
                        "pair_key": target_record["pair_key"],
                        "parameter": rational_text(parameter),
                        "construction_chart": pencil["chart"],
                        "square_multiplier": rational_text(multiplier.sqrt()),
                        "branch_quartic_coefficients_low_to_high": coefficient_text(q),
                    }
                    target_result["squareclass_hits"].append(hit)
                    trace_hits.append(hit)
            if infinity_matches:
                multiplier = inverse["squareclass_multiplier"](
                    pencil["q_lambda_family"][4], target
                )
                if multiplier.is_square():
                    hit = {
                        "target_index": int(target_index),
                        "pair_key": target_record["pair_key"],
                        "parameter": "infinity",
                        "construction_chart": pencil["chart"],
                        "square_multiplier": rational_text(multiplier.sqrt()),
                        "branch_quartic_coefficients_low_to_high": coefficient_text(
                            pencil["q_lambda_family"][4]
                        ),
                    }
                    target_result["squareclass_hits"].append(hit)
                    trace_hits.append(hit)
            resolutions.append(target_result)
        exact_resolutions.append(
            {
                "trace_index": trace_index,
                "priority_rank": int(norm8_rows[trace_index]["priority_rank"]),
                "orbit_mask": int(norm8_rows[trace_index]["orbit_mask"]),
                "target_indices": target_indices,
                "targets": resolutions,
            }
        )
        for hit in trace_hits:
            target = exact_target_record(hit["target_index"])
            left_index, right_index = target["record_indices"]
            v4_point = selector["recover_unique_affine_intersection"](
                records[left_index], records[right_index], ring
            )
            hits.append(
                {
                    "trace_index": trace_index,
                    "norm8_priority_rank": int(norm8_rows[trace_index]["priority_rank"]),
                    "norm8_orbit_mask": int(norm8_rows[trace_index]["orbit_mask"]),
                    "norm8_section_basis_w": list(map(int, words[trace_index])),
                    "v4_point": v4_point,
                    "new_character_height_on_product_twist": 8,
                    "pulled_new_height_gram_on_v4_base": [[24, 0, 0], [0, 24, 0], [0, 0, 16]],
                    "generic_rank_lower_bound": 20,
                    **hit,
                    "target": target,
                }
            )

    status = (
        "PASS_EXACT_RANK20_V4_PRODUCT_CHARACTER_HITS"
        if hits
        else (
            "PASS_EXACT_NO_NORM8_PRODUCT_CHARACTER_HIT_IN_COMPLETE_RATIONAL_V4_SET"
            if complete_native_atlas
            else "PASS_EXACT_NO_NORM8_PRODUCT_CHARACTER_HIT_IN_RATIONAL_V4_PREFIX"
        )
    )
    output = args.output if args.output.is_absolute() else ROOT / args.output
    payload = {
        "schema": "elkies-k3.r17-norm12-11952-all-rational-v4-product-inversion.v1",
        "status": status,
        "scope": {
            "native_bisection_priority_prefix": args.pool_size,
            "source_bisection_record_count": source_pool_size,
            "complete_native_rational_bisection_atlas": complete_native_atlas,
            "intersection_one_rational_v4_target_count": len(target_pairs),
            "norm8_pencil_class_count": len(norm8_rows),
            "total_target_trace_pairs": len(target_pairs) * len(norm8_rows),
        },
        "target_manifest_sha256": target_manifest_sha256,
        "primes": primes,
        "prime_candidate_counts": prime_candidate_counts,
        "trace_failure_histogram": {
            f"p{prime}:{name}": count
            for (prime, name), count in sorted(trace_failures.items())
        },
        "zero_coefficient_vector_traces": zero_vector_traces,
        "final_modular_candidate_pair_count": len(candidates),
        "exact_resolutions": exact_resolutions,
        "hit_count": len(hits),
        "hits": hits,
        "inputs": {
            relative(path): digest(path)
            for path in (
                bisection_path,
                direct_path,
                norm8_path,
                INVERSION_SCRIPT,
                SHORTLIST_SCRIPT,
                CHORD_SCRIPT,
            )
        },
        "software_assumptions": {"sage": SAGE_VERSION},
        "reproducing_command": (
            "sage -python "
            "elkies-k3/scripts/search_r17_norm12_11952_all_rational_v4_product_inversion.sage "
            f"--bisections {relative(bisection_path)} "
            f"--direct {relative(direct_path)} "
            f"--norm8-table {relative(norm8_path)} "
            f"--pool-size {args.pool_size} --primes {args.primes} "
            f"--output {relative(output)}"
        ),
        "proof_boundary": (
            f"All {len(target_pairs):,} intersection-one targets in the "
            + (
                "complete native 39,147-class rational-bisection atlas"
                if complete_native_atlas
                else f"exact native {args.pool_size:,}-class priority prefix"
            )
            + " and all 63,917 minimum-norm-eight bisection pencils are covered. A "
            "no-hit result excludes product-character sections in this integral "
            "coboundary layer only; non-coboundary twist sections and rational V4 "
            + (
                "pairs outside the rational-bisection atlas remain open. "
                if complete_native_atlas
                else "pairs outside the prefix remain open. "
            )
            + "A hit gives a rational genus-one "
            "V4 base and three orthogonal nontrivial-character sections."
        ),
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not output.exists() or output.read_text() != serialized:
            raise ArithmeticError("stored all-target inversion differs from replay")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
    print(
        "ALLV4INVERSION"
        f"|targets={len(target_pairs)}|traces={len(norm8_rows)}"
        f"|modular_candidates={len(candidates)}|hits={len(hits)}"
        f"|status={status}|output={relative(output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
