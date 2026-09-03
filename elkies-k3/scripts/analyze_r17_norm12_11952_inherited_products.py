#!/usr/bin/env python3
"""Analyze exact product characters among inherited norm-12 covers.

The squareclass calculation is performed in a formal rational-function field
with the displayed base variable renamed to one common symbol.  A separate
Neron--Severi intersection gate records whether an older character curve is
actually a degree-two cover of the alternate-Q80 base; without that gate a
formal coefficient match cannot be used as a rank certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[2]
INHERITED = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-inherited-bisection-covers-v1.json"
DIRECT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
HIDDEN_INHERITED = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-103b2-inherited-bisection-covers-v1.json"
HIDDEN_DIRECT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit103b2-direct-fibration-v1.json"
RANK28 = ROOT / "artifacts/generated-results/elkies-k3-r17-rank28-genus-one-bisection-pilot-v1.json"
SPLITTING = ROOT / "artifacts/generated-results/elkies-k3-r17-genus-one-bisection-splitting-search-v1.json"
Q103B2 = ROOT / "artifacts/generated-results/elkies-k3-norm12-orbit-103b2-twist-section-v1.json"
GRAM = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
HASHER = ROOT / "elkies-k3/scripts/hash_bisection_extensions.py"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-inherited-product-characters-v1.json"
HIDDEN_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-103b2-inherited-product-characters-v1.json"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def load_hasher():
    specification = importlib.util.spec_from_file_location("bisection_extension_hasher", HASHER)
    module = importlib.util.module_from_spec(specification)
    assert specification.loader is not None
    specification.loader.exec_module(module)
    return module


def branch(coefficients) -> dict[str, list[str]]:
    return {
        "numerator_coefficients": [str(value) for value in coefficients],
        "denominator_coefficients": ["1"],
    }


def key_tuple(key: dict) -> tuple:
    constant = key["constant_squareclass"]
    factors = key["monic_odd_factors_ascending_coefficients"]
    return (
        int(constant["sign"]),
        tuple(int(prime) for prime in constant["odd_primes"]),
        tuple(tuple(str(value) for value in factor) for factor in factors),
    )


def key_payload(key: tuple) -> dict:
    sign, primes, factors = key
    return {
        "constant_squareclass": {"sign": sign, "odd_primes": list(primes)},
        "monic_odd_factors_ascending_coefficients": [list(factor) for factor in factors],
    }


def multiply_keys(left: tuple, right: tuple) -> tuple:
    left_sign, left_primes, left_factors = left
    right_sign, right_primes, right_factors = right
    primes = set(left_primes)
    primes.symmetric_difference_update(right_primes)
    factors = set(left_factors)
    factors.symmetric_difference_update(right_factors)
    return left_sign * right_sign, tuple(sorted(primes)), tuple(sorted(factors))


def key_digest(key: tuple) -> str:
    return hashlib.sha256(
        json.dumps(key_payload(key), sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def load_gram() -> list[list[int]]:
    return [
        [int(value) for value in line.split()]
        for line in GRAM.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]


def pairing(left: list[int], right: list[int], gram: list[list[int]]) -> int:
    return sum(
        left[i] * gram[i][j] * right[j]
        for i in range(17) for j in range(17)
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-label",
        choices=("norm12-orbit-11952", "norm12-orbit-103b2"),
        default="norm12-orbit-11952",
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()

    if arguments.source_label == "norm12-orbit-11952":
        inherited_path = INHERITED
        direct_path = DIRECT
        expected_status = "PASS_EXACT_121_INHERITED_ALTERNATE_Q80_BISECTION_COVERS"
        expected_count = 121
        expected_older_degrees = (8, 9)
        output = arguments.output or DEFAULT_OUTPUT
        schema = "elkies-k3.r17-norm12-11952-inherited-product-characters.v1"
        target_description = "alternate-Q80"
    else:
        inherited_path = HIDDEN_INHERITED
        direct_path = HIDDEN_DIRECT
        expected_status = "PASS_EXACT_82_INHERITED_103B2_MARKING_BISECTION_COVERS"
        expected_count = 82
        expected_older_degrees = (12, 0)
        output = arguments.output or HIDDEN_OUTPUT
        schema = "elkies-k3.r17-norm12-103b2-inherited-product-characters.v1"
        target_description = "hidden-103b2"

    inherited = json.loads(inherited_path.read_text())
    direct = json.loads(direct_path.read_text())
    rank28 = json.loads(RANK28.read_text())
    splitting = json.loads(SPLITTING.read_text())
    q103b2 = json.loads(Q103B2.read_text())
    if inherited.get("status") != expected_status:
        raise ValueError("the inherited-cover input has the wrong exact certificate status")
    if len(inherited["bisections"]) != expected_count:
        raise ValueError(f"expected exactly {expected_count} inherited covers")

    hasher = load_hasher()
    variable = sp.Symbol("z")
    individual: dict[str, tuple] = {}
    labels_by_key: dict[tuple, list[str]] = {}
    for record in inherited["bisections"]:
        label = str(record["label"])
        coefficients = record["canonical_squareclass"]["q_coefficients_low_to_high"]
        key = key_tuple(hasher.extension_key(branch(coefficients), variable))
        individual[label] = key
        labels_by_key.setdefault(key, []).append(label)
    if len(labels_by_key) != expected_count:
        raise ArithmeticError("the inherited covers are not pairwise-distinct squareclasses")

    product_buckets: dict[tuple, list[tuple[str, str]]] = {}
    inherited_matches = []
    ordered_labels = sorted(individual)
    for left_index, left in enumerate(ordered_labels):
        for right in ordered_labels[left_index + 1:]:
            product = multiply_keys(individual[left], individual[right])
            product_buckets.setdefault(product, []).append((left, right))
            if product in labels_by_key:
                inherited_matches.append({
                    "pair": [left, right],
                    "matching_inherited_characters": labels_by_key[product],
                    "product_extension_sha256": key_digest(product),
                })
    product_collisions = [
        {
            "product_extension_sha256": key_digest(key),
            "pairs": [list(pair) for pair in pairs],
            "pair_count": len(pairs),
        }
        for key, pairs in product_buckets.items() if len(pairs) > 1
    ]
    product_collisions.sort(key=lambda item: item["product_extension_sha256"])

    older_catalog = []
    trace = rank28["traces"][0]
    for target in trace["targets"]:
        coefficients = target["branch_polynomial_q_coefficients_low_to_high"]
        key = key_tuple(hasher.extension_key(branch(coefficients), variable))
        older_catalog.append({
            "label": f"rank28-{target['target_label']}",
            "source": relative(RANK28),
            "formal_extension_key": key,
            "formal_extension_sha256": key_digest(key),
            "source_base": "published R17 t-line",
        })
    q103_key = key_tuple(hasher.extension_key(branch(q103b2["q_coefficients_low_to_high"]), variable))
    older_catalog.append({
        "label": "q_103b2",
        "source": relative(Q103B2),
        "formal_extension_key": q103_key,
        "formal_extension_sha256": key_digest(q103_key),
        "source_base": "published R17 t-line",
    })

    formal_matches = []
    for older in older_catalog:
        pairs = product_buckets.get(older["formal_extension_key"], [])
        if pairs:
            formal_matches.append({
                "catalog_label": older["label"],
                "pairs": [list(pair) for pair in pairs],
                "formal_extension_sha256": older["formal_extension_sha256"],
            })
    for older in older_catalog:
        older.pop("formal_extension_key")

    gram = load_gram()
    target_w = [int(value) for value in direct["divisor"]["pinned_trace_vector_w"]]
    norm8_w = [int(value) for value in trace["pinned_rank17_w"]]
    q103_record = next(
        record for record in splitting["construction"]["records"]
        if record.get("label") == "norm12-orbit-103b2"
    )
    q103_w = [int(value) for value in q103_record["pinned_rank17_w"]]
    if pairing(target_w, target_w, gram) != 12:
        raise ArithmeticError("target divisor trace does not have norm 12")
    if pairing(norm8_w, norm8_w, gram) != 8:
        raise ArithmeticError("rank-28 pencil trace does not have norm 8")
    if pairing(q103_w, q103_w, gram) != 12:
        raise ArithmeticError("q_103b2 trace does not have norm 12")
    rank28_dot = pairing(target_w, norm8_w, gram)
    q103_dot = pairing(target_w, q103_w, gram)
    rank28_target_degree = 2 * 2 + 2 * 3 - rank28_dot
    q103_target_degree = 3 * 2 + 2 * 3 - q103_dot
    if (rank28_target_degree, q103_target_degree) != expected_older_degrees:
        raise ArithmeticError("older-character target-base degree calculation changed")

    result = {
        "schema": schema,
        "status": "PASS_EXACT_NO_INHERITED_PRODUCT_CHARACTER_CLOSURE",
        "inputs": {
            relative(path): digest(path)
            for path in (inherited_path, direct_path, RANK28, SPLITTING, Q103B2, GRAM, HASHER)
        },
        "inherited_character_group": {
            "character_count": len(individual),
            "distinct_character_count": len(labels_by_key),
            "pair_product_count": len(ordered_labels) * (len(ordered_labels) - 1) // 2,
            "distinct_pair_product_count": len(product_buckets),
            "pair_product_collision_count": len(product_collisions),
            "pair_product_collisions": product_collisions,
            "matches_another_inherited_character_count": len(inherited_matches),
            "matches_another_inherited_character": inherited_matches,
        },
        "older_published_base_catalog": {
            "formal_variable_rename_comparison_count": len(older_catalog),
            "characters": older_catalog,
            "formal_product_match_count": len(formal_matches),
            "formal_product_matches": formal_matches,
            "compatibility_gate": {
                "pairing_convention": "(a,b,w).(a',b',w')=a*b'+b*a'-<w,w'>_R17",
                "target_fibre": target_description,
                "target_fibre_divisor": {"a": 3, "b": 2, "trace_norm": 12},
                "rank28_pencil": {
                    "source_divisor": {"a": 2, "b": 2, "trace_norm": 8},
                    "trace_pairing_with_target": rank28_dot,
                    "degree_over_target_u_line": rank28_target_degree,
                    "is_quadratic_character_of_target_base": rank28_target_degree == 2,
                },
                "q_103b2": {
                    "source_divisor": {"a": 3, "b": 2, "trace_norm": 12},
                    "trace_pairing_with_target": q103_dot,
                    "degree_over_target_u_line": q103_target_degree,
                    "is_quadratic_character_of_target_base": q103_target_degree == 2,
                },
                "consequence": (
                    "The eleven rank-28 quartics and q_103b2 are characters over the old "
                    f"published t-line. Their source curves have degrees {rank28_target_degree} "
                    f"and {q103_target_degree} over the {target_description} u-line, not degree "
                    "2, so renaming t to u is only a formal coefficient comparison and cannot "
                    "supply the third V4 character."
                ),
            },
        },
        "proof_boundary": (
            f"All {len(product_buckets)} products of the {expected_count} exact inherited "
            f"{target_description}-base quadratic "
            "characters are compared in QQ(u)^*/QQ(u)^{*2}, including the rational "
            "constant squareclass. Older published-base characters are also compared "
            "after a formal variable rename, but the exact NS degree gate prevents "
            "treating those old-base equations as target-base twist characters."
        ),
        "reproducing_command": (
            ".venv/bin/python elkies-k3/scripts/analyze_r17_norm12_11952_inherited_products.py "
            f"--source-label {arguments.source_label} --output {relative(output)}"
        ),
    }
    if inherited_matches:
        result["status"] = "PASS_EXACT_INHERITED_PRODUCT_CHARACTER_CLOSURE_FOUND"
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not output.exists() or output.read_text() != serialized:
            raise ValueError("stored inherited-product artifact differs from replay")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
    print(
        "ALTPRODUCTS|characters={}|products={}|distinct_products={}|"
        "product_collisions={}|inherited_closures={}|formal_old_matches={}|status={}".format(
            len(individual), len(ordered_labels) * (len(ordered_labels) - 1) // 2,
            len(product_buckets), len(product_collisions), len(inherited_matches),
            len(formal_matches), result["status"],
        )
    )


if __name__ == "__main__":
    main()
