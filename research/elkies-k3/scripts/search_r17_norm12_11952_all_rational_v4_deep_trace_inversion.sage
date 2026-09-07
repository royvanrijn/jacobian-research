#!/usr/bin/env sage-python
"""Compare every rational intersection-one V4 product with all 49 deep traces.

The complete native rational-bisection atlas contains 4,358,409 pairs with
intersection number one.  Their product characters are quartics and their V4
fibre products are genus-one curves with a rational point.  The separate
norm-eight inversion covers the minimum-trace coboundary layer.  This script
hashes every one of these product quartics against every squareclass in the
one-parameter families belonging to the 49 minimum-norm-twelve trace parities.

Finite-field equality retains the scalar squareclass.  A pair is rejected
only at a prime where its target and the deep trace have good reduction; all
modular survivors are resolved by exact polynomial division over QQ.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from hashlib import sha256
import json
from pathlib import Path
import runpy

import numpy as np
from sage.all import GF, PolynomialRing, QQ, ZZ, matrix, vector
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = Path(__file__).resolve().parent
DIRECT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
FULL_BISECTIONS = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisections-full-v1.json"
COLLISIONS = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisection-collisions-full-v1.json"
PRIORITY = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisection-priority-v1.tsv"
PARITY = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-product-tate-parity-v1.json"
DEEP_SCRIPT = SCRIPTS / "search_r17_norm12_11952_product_deep_trace_inversion.sage"
INVERSION_SCRIPT = SCRIPTS / "search_r17_norm12_11952_product_bisection_inversion.sage"
HELPER = SCRIPTS / "construct_elkies_2026_bisections.sage"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-deep-trace-inversion-full-v1.json"


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def rational_text(value) -> str:
    value = QQ(value)
    return str(value.numerator()) if value.denominator() == 1 else f"{value.numerator()}/{value.denominator()}"


def coefficient_text(polynomial, length: int) -> list[str]:
    return [rational_text(polynomial[index]) for index in range(length)]


def stream_branches(path: Path):
    """Read only each branch triple and its following label from the large JSON."""

    pending = None
    with path.open() as stream:
        iterator = iter(stream)
        for line in iterator:
            if '"numerator_coefficients": [' in line:
                values = []
                for coefficient_line in iterator:
                    text = coefficient_line.strip().rstrip(",")
                    if text == "]":
                        break
                    values.append(QQ(json.loads(text)))
                if len(values) == 3:
                    # The content and its sign are part of the rational
                    # squareclass.  Never primitive-normalize this triple.
                    pending = tuple(values)
            elif pending is not None and '"label": ' in line:
                label = json.loads(line.split(":", 1)[1].strip().rstrip(","))
                yield label, pending
                pending = None


def read_priority(path: Path):
    rows = []
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            rows.append(
                {
                    "priority_rank": int(row["priority_rank"]),
                    "label": f"alternate-orbit-{int(row['orbit_mask']):05x}",
                    "orbit_mask": int(row["orbit_mask"]),
                    "direct_w": tuple(map(int, row["direct_alternate_w"].split())),
                }
            )
    rows.sort(key=lambda row: row["priority_rank"])
    if [row["priority_rank"] for row in rows] != list(range(1, len(rows) + 1)):
        raise ArithmeticError("priority table is not a complete initial interval")
    return rows


def scalar_hash(values, prime: int, inverses, powers):
    """Hash projective coefficients together with the pivot square character."""

    values = [int(value) % prime for value in values]
    pivot = next((value for value in values if value), None)
    if pivot is None:
        return None
    normalized = [(value * inverses[pivot]) % prime for value in values]
    nonsquare = 0 if pow(pivot, (prime - 1) // 2, prime) == 1 else 1
    return sum(value * powers[index] for index, value in enumerate(normalized)) + nonsquare * powers[5]


def vector_hash(values, prime: int, inverses, powers):
    """Vectorized version of scalar_hash for an n-by-5 int64 array."""

    values = np.remainder(values, prime).astype(np.int64, copy=False)
    pivots = np.zeros(values.shape[0], dtype=np.int64)
    unset = np.ones(values.shape[0], dtype=bool)
    for column in range(5):
        take = unset & (values[:, column] != 0)
        pivots[take] = values[take, column]
        unset[take] = False
    if np.any(unset):
        raise ArithmeticError("a primitive target quartic vanished modulo p")
    normalized = np.remainder(values * inverses[pivots, None], prime)
    hashes = normalized @ powers[:5]
    character = np.array(
        [0 if pow(int(value), (prime - 1) // 2, prime) == 1 else 1 for value in pivots],
        dtype=np.int64,
    )
    return hashes + character * powers[5]


def deep_keys(prime, direct, deep_records, inverse, helper, deep):
    maximum = max(abs(int(value)) for row in deep_records for value in row["section_basis_w"])
    context = inverse["modular_context"](prime, direct, helper, maximum)
    inverses = np.array([0] + [pow(value, -1, prime) for value in range(1, prime)], dtype=np.int64)
    powers = np.array([prime**index for index in range(6)], dtype=np.int64)
    by_chart_and_hash = defaultdict(lambda: defaultdict(list))
    bad = {}
    zero_parameters = {}
    for deep_index, record in enumerate(deep_records):
        try:
            pencil = deep["deep_trace_family"](
                context,
                vector(ZZ, record["section_basis_w"]),
                inverse,
                helper,
            )
        except (ArithmeticError, ValueError, ZeroDivisionError) as error:
            bad[deep_index] = type(error).__name__
            continue
        ring = context["ring"]
        parameters = list(context["coefficient_field"]) + ["infinity"]
        for parameter in parameters:
            q = deep["evaluated_family"](pencil["q_lambda_family"], parameter, ring)
            parameter_text = "infinity" if parameter == "infinity" else str(int(parameter))
            if not q:
                zero_parameters.setdefault(deep_index, []).append(parameter_text)
                continue
            factorization = q.factor()
            reduced = ring(factorization.unit())
            for factor, exponent in factorization:
                if int(exponent) % 2:
                    reduced *= factor
            if reduced.degree() != 4:
                continue
            key = scalar_hash([reduced[index] for index in range(5)], prime, inverses, powers)
            by_chart_and_hash[pencil["chart"]][key].append((deep_index, parameter_text))
    return context, by_chart_and_hash, bad, zero_parameters, inverses, powers


def branch_array(rows, branches, prime):
    result = np.empty((len(rows), 3), dtype=np.int64)
    for index, row in enumerate(rows):
        coefficients = branches[row["label"]]
        reduced = []
        for value in coefficients:
            value = QQ(value)
            if int(value.denominator()) % prime == 0:
                raise ZeroDivisionError("branch coefficient has bad reduction")
            reduced.append(
                int(value.numerator())
                * pow(int(value.denominator()) % prime, -1, prime)
                % prime
            )
        if not any(reduced):
            raise ZeroDivisionError("an entire branch polynomial vanishes modulo p")
        result[index] = reduced
    return result


def reduce_values(values, prime):
    reduced = []
    for value in values:
        value = QQ(value)
        if int(value.denominator()) % prime == 0:
            raise ZeroDivisionError("rational polynomial has bad reduction")
        reduced.append(
            int(value.numerator())
            * pow(int(value.denominator()) % prime, -1, prime)
            % prime
        )
    return reduced


def product_coefficients(left, right):
    return np.column_stack(
        (
            left[:, 0] * right[:, 0],
            left[:, 0] * right[:, 1] + left[:, 1] * right[:, 0],
            left[:, 0] * right[:, 2] + left[:, 1] * right[:, 1] + left[:, 2] * right[:, 0],
            left[:, 1] * right[:, 2] + left[:, 2] * right[:, 1],
            left[:, 2] * right[:, 2],
        )
    )


def scan_first_prime(rows, vectors, gram_array, branch_modp, key_maps, prime, inverses, powers):
    candidates = set()
    target_count = 0
    direct_keys = key_maps.get("finite", {})
    inverted_keys = key_maps.get("inverted_at_infinity", {})
    direct_hashes = np.array(sorted(direct_keys), dtype=np.int64)
    inverted_hashes = np.array(sorted(inverted_keys), dtype=np.int64)
    covectors = vectors @ gram_array
    block = 512
    for right_start in range(0, len(rows), block):
        right_stop = min(len(rows), right_start + block)
        pairings = covectors[right_start:right_stop] @ vectors.T
        for offset, pairing_row in enumerate(pairings):
            right_index = right_start + offset
            left_indices = np.flatnonzero(pairing_row[:right_index] == 7)
            if not len(left_indices):
                continue
            target_count += len(left_indices)
            right = np.repeat(branch_modp[right_index][None, :], len(left_indices), axis=0)
            products = product_coefficients(branch_modp[left_indices], right)
            hashes = vector_hash(products, prime, inverses, powers)
            if len(direct_hashes):
                for position in np.flatnonzero(np.isin(hashes, direct_hashes)):
                    for deep_index, unused_parameter in direct_keys[int(hashes[position])]:
                        candidates.add((int(left_indices[position]), right_index, deep_index))
            if len(inverted_hashes):
                reverse_hashes = vector_hash(products[:, ::-1], prime, inverses, powers)
                for position in np.flatnonzero(np.isin(reverse_hashes, inverted_hashes)):
                    for deep_index, unused_parameter in inverted_keys[int(reverse_hashes[position])]:
                        candidates.add((int(left_indices[position]), right_index, deep_index))
    return target_count, candidates


def filter_candidates(candidates, branch_modp, key_maps, prime, inverses, powers, bad, zeros):
    kept = set()
    for left_index, right_index, deep_index in candidates:
        if deep_index in bad or deep_index in zeros:
            kept.add((left_index, right_index, deep_index))
            continue
        products = product_coefficients(
            branch_modp[left_index][None, :], branch_modp[right_index][None, :]
        )
        finite_hash = int(vector_hash(products, prime, inverses, powers)[0])
        inverted_hash = int(vector_hash(products[:, ::-1], prime, inverses, powers)[0])
        if (
            deep_index in {item[0] for item in key_maps.get("finite", {}).get(finite_hash, ())}
            or deep_index in {item[0] for item in key_maps.get("inverted_at_infinity", {}).get(inverted_hash, ())}
        ):
            kept.add((left_index, right_index, deep_index))
    return kept


def exact_resolve(candidates, rows, branches, deep_records, exact_context, inverse, helper, deep):
    records = []
    hits = []
    unresolved_symbolic = []
    ring = exact_context["ring"]
    grouped = defaultdict(list)
    for left_index, right_index, deep_index in sorted(candidates):
        grouped[deep_index].append((left_index, right_index))
    for deep_index, pairs in grouped.items():
        record = deep_records[deep_index]
        pencil = deep["deep_trace_family"](
            exact_context,
            vector(ZZ, record["section_basis_w"]),
            inverse,
            helper,
        )
        for left_index, right_index in pairs:
            target = ring(branches[rows[left_index]["label"]]) * ring(
                branches[rows[right_index]["label"]]
            )
            if pencil["chart"] == "inverted_at_infinity":
                target = inverse["reciprocal_polynomial"](target, 4, ring)
            parameters, unused = deep["exact_parameters"](pencil, target, exact_context)
            if parameters is None:
                unresolved_symbolic.append(
                    {
                        "deep_trace_index": deep_index,
                        "pair_key": f"{rows[left_index]['label']}:{rows[right_index]['label']}",
                    }
                )
                parameters = []
            row_hits = []
            for parameter in parameters:
                q = deep["evaluated_family"](pencil["q_lambda_family"], parameter, ring)
                square_root = deep["squarefactor_match"](q, target)
                if square_root is not None and square_root is not False:
                    hit = {
                        "deep_trace_index": deep_index,
                        "deep_orbit_mask": int(record["orbit_mask"]),
                        "pair_key": f"{rows[left_index]['label']}:{rows[right_index]['label']}",
                        "priority_ranks": [left_index + 1, right_index + 1],
                        "lambda": rational_text(parameter),
                        "construction_chart": pencil["chart"],
                        "target_quartic_coefficients_low_to_high": coefficient_text(target, 5),
                        "square_factor_coefficients_low_to_high": coefficient_text(square_root, 3),
                    }
                    row_hits.append(hit)
                    hits.append(hit)
            infinity_root = deep["squarefactor_match"](
                pencil["q_lambda_family"][4], target
            )
            if infinity_root is not None and infinity_root is not False:
                hit = {
                    "deep_trace_index": deep_index,
                    "deep_orbit_mask": int(record["orbit_mask"]),
                    "pair_key": f"{rows[left_index]['label']}:{rows[right_index]['label']}",
                    "priority_ranks": [left_index + 1, right_index + 1],
                    "lambda": "infinity",
                    "construction_chart": pencil["chart"],
                    "target_quartic_coefficients_low_to_high": coefficient_text(target, 5),
                    "square_factor_coefficients_low_to_high": coefficient_text(infinity_root, 3),
                }
                row_hits.append(hit)
                hits.append(hit)
            records.append(
                {
                    "deep_trace_index": deep_index,
                    "pair_key": f"{rows[left_index]['label']}:{rows[right_index]['label']}",
                    "rational_divisibility_parameter_count": len(parameters),
                    "squareclass_hit_count": len(row_hits),
                }
            )
    return records, hits, unresolved_symbolic


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primes", default="131,137")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    primes = [int(value) for value in args.primes.split(",") if value]
    if not primes or any(not ZZ(prime).is_prime() for prime in primes):
        parser.error("--primes must be a nonempty prime list")

    direct = json.loads(DIRECT.read_text())
    parity = json.loads(PARITY.read_text())
    collision = json.loads(COLLISIONS.read_text())
    if direct.get("status") != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
        raise ArithmeticError("direct model is not certified")
    if parity.get("status") != "PASS_EXACT_PRODUCT_TATE_PARITY_REDUCTION":
        raise ArithmeticError("deep parity input is not certified")
    if collision.get("collision_count") != 0 or collision.get("input_bisection_count") != 39147:
        raise ArithmeticError("complete rational-bisection atlas is not certified")
    if digest(FULL_BISECTIONS) != collision["input"]["sha256"]:
        raise ArithmeticError("large bisection atlas differs from its compact certificate")

    rows = read_priority(PRIORITY)
    branches = dict(stream_branches(FULL_BISECTIONS))
    if len(rows) != 39147 or len(branches) != 39147 or set(branches) != {row["label"] for row in rows}:
        raise ArithmeticError("priority vectors and streamed branches do not cover the same atlas")
    vectors = np.array([row["direct_w"] for row in rows], dtype=np.int64)
    gram_array = np.array(direct["frame_certificate"]["frame_gram"], dtype=np.int64)
    if not np.all(np.einsum("ij,ij->i", vectors @ gram_array, vectors) == 10):
        raise ArithmeticError("a rational bisection vector has the wrong norm")

    deep_records = parity["invariant_trace_parity"]["deep_norm12_classes"]
    if len(deep_records) != 49:
        raise ArithmeticError("expected 49 deep trace parities")
    deep = runpy.run_path(str(DEEP_SCRIPT))
    inverse = runpy.run_path(str(INVERSION_SCRIPT))
    helper = runpy.run_path(str(HELPER))
    exact_context = inverse["exact_model_and_basis"](direct, helper)
    control = deep["build_positive_control"](
        exact_context, deep_records[0], inverse, helper
    )
    control_deep_index = next(
        index
        for index, record in enumerate(deep_records)
        if int(record["orbit_mask"]) == int(control["orbit_mask"])
    )

    prime_records = []
    control_records = []
    candidates = None
    target_count = None
    for prime in primes:
        unused_context, key_maps, bad, zeros, inverses, powers = deep_keys(
            prime, direct, deep_records, inverse, helper, deep
        )
        if candidates is None and (bad or zeros):
            raise ArithmeticError(
                "the first prime must give good, nonzero reductions for all deep traces"
            )
        control_values = reduce_values(
            control["branch_quartic_coefficients_low_to_high"], prime
        )
        expected_parameter = str(
            int(GF(prime)(QQ(control["pencil_parameter_lambda"])))
        )
        control_matches = []
        for chart, values in (
            ("finite", control_values),
            ("inverted_at_infinity", list(reversed(control_values))),
        ):
            key = scalar_hash(values, prime, inverses, powers)
            for deep_index, parameter in key_maps.get(chart, {}).get(key, ()):
                if deep_index == control_deep_index:
                    control_matches.append((chart, parameter))
        if not any(parameter == expected_parameter for unused, parameter in control_matches):
            raise ArithmeticError("streaming target hash failed the deep positive control")
        control_records.append(
            {
                "prime": prime,
                "deep_trace_index": control_deep_index,
                "expected_parameter_residue": expected_parameter,
                "matching_chart_parameter_pairs": [list(item) for item in control_matches],
            }
        )
        branch_modp = branch_array(rows, branches, prime)
        before = None if candidates is None else len(candidates)
        if candidates is None:
            target_count, candidates = scan_first_prime(
                rows, vectors, gram_array, branch_modp, key_maps, prime, inverses, powers
            )
        else:
            candidates = filter_candidates(
                candidates, branch_modp, key_maps, prime, inverses, powers, bad, zeros
            )
        prime_records.append(
            {
                "prime": prime,
                "candidate_count_before": before,
                "deep_trace_bad_reduction_count": len(bad),
                "deep_zero_family_parameter_count": sum(map(len, zeros.values())),
                "deep_quartic_squareclass_key_count": sum(len(keys) for keys in key_maps.values()),
                "candidate_count_after": len(candidates),
            }
        )
        if not candidates:
            break

    if target_count != 4358409:
        raise ArithmeticError(f"expected 4,358,409 intersection-one targets, found {target_count}")
    exact_records, hits, unresolved_symbolic = exact_resolve(
        candidates, rows, branches, deep_records, exact_context, inverse, helper, deep
    )
    status = (
        "PASS_EXACT_DEEP_TRACE_PRODUCT_CHARACTER_HITS"
        if hits
        else (
            "INCOMPLETE_SYMBOLIC_DEEP_TRACE_FAMILIES"
            if unresolved_symbolic
            else "PASS_EXACT_NO_DEEP_TRACE_HIT_IN_COMPLETE_RATIONAL_V4_SET"
        )
    )
    payload = {
        "schema": "elkies-k3.r17-norm12-11952-all-rational-v4-deep-trace-inversion.v1",
        "status": status,
        "scope": {
            "complete_native_rational_bisection_count": len(rows),
            "intersection_one_rational_v4_target_count": target_count,
            "deep_norm12_trace_class_count": len(deep_records),
            "target_trace_comparison_count": target_count * len(deep_records),
        },
        "primes": prime_records,
        "positive_control": {
            "label": control["label"],
            "deep_orbit_mask": int(control["orbit_mask"]),
            "pencil_parameter_lambda": control["pencil_parameter_lambda"],
            "branch_quartic_coefficients_low_to_high": control[
                "branch_quartic_coefficients_low_to_high"
            ],
            "modular_recoveries": control_records,
            "status": "PASS_STREAMING_HASH_POSITIVE_CONTROL_RECOVERED",
        },
        "final_modular_candidate_count": len(candidates),
        "exact_resolutions": exact_records,
        "unresolved_symbolic_families": unresolved_symbolic,
        "hit_count": len(hits),
        "hits": hits,
        "inputs": {
            relative(path): digest(path)
            for path in (DIRECT, FULL_BISECTIONS, COLLISIONS, PRIORITY, PARITY, DEEP_SCRIPT, INVERSION_SCRIPT, HELPER)
        },
        "software_assumptions": {"sage": SAGE_VERSION, "numpy": np.__version__},
        "reproducing_command": (
            "sage -python elkies-k3/scripts/"
            "search_r17_norm12_11952_all_rational_v4_deep_trace_inversion.sage "
            f"--primes {','.join(map(str, primes))} --output {relative(args.output)}"
        ),
        "proof_boundary": (
            "Every one of the 4,358,409 intersection-one rational V4 product "
            "characters is compared with every parameter in all 49 deep norm-twelve "
            "trace families. Modular keys retain the rational scalar squareclass, "
            "and every modular survivor is resolved over QQ. A no-hit result excludes "
            "only this deep integral coboundary layer; non-coboundary product-twist "
            "sections and rational V4 bases outside the native rational-bisection "
            "atlas remain open."
        ),
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = args.output.resolve()
    if args.check:
        if not output.exists() or output.read_text() != serialized:
            raise ArithmeticError("stored all-target deep inversion differs from replay")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
    print(
        "ALLV4DEEPINVERSION|"
        f"targets={target_count}|deep_traces={len(deep_records)}|"
        f"modular_candidates={len(candidates)}|hits={len(hits)}|"
        f"status={status}|output={relative(output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
