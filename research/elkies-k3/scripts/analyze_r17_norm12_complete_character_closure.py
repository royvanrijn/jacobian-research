#!/usr/bin/env python3
"""Search a complete norm-12 bisection frame for three-character closure.

Every accepted cover has a squareclass whose geometric branch divisor has
degree two.  Its finite irreducible-factor support consequently has size at
most two.  A relation q_i*q_j=q_k can therefore occur only as either an edge
and its two singleton endpoints, or as a triangle in the graph of two-factor
supports.  This gives an exhaustive exact search without forming all O(n^2)
pair products.  Rational constant squareclasses remain part of every key.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import importlib.util
import json
from pathlib import Path
from typing import Any

import sympy as sp


ROOT = Path(__file__).resolve().parents[2]
HASHER = ROOT / "elkies-k3/scripts/hash_bisection_extensions.py"
GRAM = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
RANK28 = ROOT / "artifacts/generated-results/elkies-k3-r17-rank28-genus-one-bisection-pilot-v1.json"
SPLITTING = ROOT / "artifacts/generated-results/elkies-k3-r17-genus-one-bisection-splitting-search-v1.json"
Q103B2 = ROOT / "artifacts/generated-results/elkies-k3-norm12-orbit-103b2-twist-section-v1.json"

CONFIG = {
    "norm12-orbit-11952": {
        "input": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisections-full-v1.json",
        "direct": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json",
        "output": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-complete-character-closure-v1.json",
        "status": "PASS_EXACT_COMPLETE_ALTERNATE_BISECTION_EQUATIONS",
        "count": 39147,
        "older_degrees": (8, 9),
        "description": "alternate-Q80",
    },
    "norm12-orbit-103b2": {
        "input": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-103b2-bisections-full-v1.json",
        "direct": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit103b2-direct-fibration-v1.json",
        "output": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-103b2-complete-character-closure-v1.json",
        "status": "PASS_EXACT_COMPLETE_103B2_HIDDEN_BISECTION_EQUATIONS",
        "count": 39120,
        "older_degrees": (12, 0),
        "description": "hidden-103b2",
    },
    "norm12-orbit-08f72": {
        "input": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-08f72-alternate-bisections-full-v1.json",
        "direct": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit08f72-direct-fibration-v1.json",
        "output": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-08f72-complete-character-closure-v1.json",
        "status": "PASS_EXACT_COMPLETE_ALTERNATE_BISECTION_EQUATIONS",
        "count": 39147,
        "older_degrees": (8, 12),
        "description": "alternate-Q80 orbit-08f72",
    },
    "norm12-orbit-08ab4": {
        "input": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-08ab4-alternate-bisections-full-v1.json",
        "direct": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit08ab4-direct-fibration-v1.json",
        "output": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-08ab4-complete-character-closure-v1.json",
        "status": "PASS_EXACT_COMPLETE_ALTERNATE_BISECTION_EQUATIONS",
        "count": 39147,
        "older_degrees": (10, 10),
        "description": "alternate-Q80 orbit-08ab4",
    },
    "norm12-orbit-091e4": {
        "input": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-091e4-alternate-bisections-full-v1.json",
        "direct": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit091e4-direct-fibration-v1.json",
        "output": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-091e4-complete-character-closure-v1.json",
        "status": "PASS_EXACT_COMPLETE_ALTERNATE_BISECTION_EQUATIONS",
        "count": 39147,
        "older_degrees": (9, 12),
        "description": "alternate-Q80 orbit-091e4",
    },
}

Atom = tuple[str, ...]
Key = frozenset[Atom]
Signature = frozenset[Atom]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def load_module(path: Path, name: str):
    specification = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(specification)
    assert specification.loader is not None
    specification.loader.exec_module(module)
    return module


def atoms(extension_key: dict[str, Any]) -> Key:
    result: set[Atom] = set()
    constant = extension_key["constant_squareclass"]
    if int(constant["sign"]) == -1:
        result.add(("sign", "-1"))
    result.update(("prime", str(prime)) for prime in constant["odd_primes"])
    result.update(
        ("polynomial", *(str(value) for value in factor))
        for factor in extension_key["monic_odd_factors_ascending_coefficients"]
    )
    return frozenset(result)


def polynomial_signature(key: Key) -> Signature:
    return frozenset(atom for atom in key if atom[0] == "polynomial")


def branch(coefficients: list[Any]) -> dict[str, list[str]]:
    return {
        "numerator_coefficients": [str(value) for value in coefficients],
        "denominator_coefficients": ["1"],
    }


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


def exact_relations(
    labels_by_key: dict[Key, str], keys_by_signature: dict[Signature, list[Key]]
) -> tuple[list[list[str]], dict[str, int]]:
    """Return every three-key relation, using the degree-two support graph."""

    relations: set[tuple[str, str, str]] = set()
    singleton = {
        next(iter(signature)): signature
        for signature in keys_by_signature if len(signature) == 1
    }
    edges = {signature for signature in keys_by_signature if len(signature) == 2}

    endpoint_candidate_count = 0
    for edge in edges:
        left_atom, right_atom = tuple(edge)
        left_signature = singleton.get(left_atom)
        right_signature = singleton.get(right_atom)
        if left_signature is None or right_signature is None:
            continue
        for left in keys_by_signature[left_signature]:
            for right in keys_by_signature[right_signature]:
                endpoint_candidate_count += 1
                third = left ^ right
                if third in labels_by_key:
                    relations.add(tuple(sorted((
                        labels_by_key[left], labels_by_key[right], labels_by_key[third]
                    ))))

    adjacency: dict[Atom, set[Atom]] = defaultdict(set)
    for edge in edges:
        left_atom, right_atom = tuple(edge)
        adjacency[left_atom].add(right_atom)
        adjacency[right_atom].add(left_atom)
    triangle_support_count = 0
    triangle_key_candidate_count = 0
    for left_atom in sorted(adjacency):
        for right_atom in sorted(atom for atom in adjacency[left_atom] if left_atom < atom):
            common = adjacency[left_atom].intersection(adjacency[right_atom])
            for third_atom in sorted(atom for atom in common if right_atom < atom):
                triangle_support_count += 1
                left_signature = frozenset((left_atom, right_atom))
                right_signature = frozenset((left_atom, third_atom))
                third_signature = frozenset((right_atom, third_atom))
                for left in keys_by_signature[left_signature]:
                    for right in keys_by_signature[right_signature]:
                        triangle_key_candidate_count += 1
                        third = left ^ right
                        if third in labels_by_key and polynomial_signature(third) == third_signature:
                            relations.add(tuple(sorted((
                                labels_by_key[left], labels_by_key[right], labels_by_key[third]
                            ))))
    return [list(relation) for relation in sorted(relations)], {
        "singleton_support_count": len(singleton),
        "two_factor_support_count": len(edges),
        "endpoint_key_candidate_count": endpoint_candidate_count,
        "triangle_support_count": triangle_support_count,
        "triangle_key_candidate_count": triangle_key_candidate_count,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-label", choices=tuple(CONFIG), default="norm12-orbit-11952")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    configuration = CONFIG[arguments.source_label]
    input_path = configuration["input"]
    direct_path = configuration["direct"]
    output = arguments.output or configuration["output"]

    payload = json.loads(input_path.read_text())
    if payload.get("status") != configuration["status"]:
        raise ValueError("complete cover artifact has the wrong exact status")
    records = payload.get("bisections")
    if not isinstance(records, list) or len(records) != configuration["count"]:
        raise ValueError("complete cover artifact has the wrong record count")

    hasher = load_module(HASHER, "bisection_extension_hasher")
    variable = sp.Symbol(str(payload.get("base_parameter", "u")))
    labels_by_key: dict[Key, str] = {}
    keys_by_signature: dict[Signature, list[Key]] = defaultdict(list)
    factor_count_histogram: dict[int, int] = defaultdict(int)
    for record in records:
        declared_branch, _ = hasher.branch_from_record(record, variable)
        extension = hasher.extension_key(declared_branch, variable)
        if hasher.is_trivial_squareclass(extension):
            raise ArithmeticError(f"{record['label']}: split squareclass")
        divisor = hasher.geometric_branch_divisor(declared_branch, variable)
        if divisor["geometric_degree"] != 2:
            raise ArithmeticError(f"{record['label']}: geometric branch degree is not two")
        key = atoms(extension)
        if key in labels_by_key:
            raise ArithmeticError("complete character input contains an equal-cover collision")
        signature = polynomial_signature(key)
        if not 1 <= len(signature) <= 2:
            raise ArithmeticError("degree-two cover has an unexpected finite-factor support")
        labels_by_key[key] = str(record["label"])
        keys_by_signature[signature].append(key)
        factor_count_histogram[len(signature)] += 1

    relations, search_counts = exact_relations(labels_by_key, keys_by_signature)

    rank28 = json.loads(RANK28.read_text())
    q103b2 = json.loads(Q103B2.read_text())
    older_catalog: list[tuple[str, Key]] = []
    trace = rank28["traces"][0]
    for target in trace["targets"]:
        extension = hasher.extension_key(
            branch(target["branch_polynomial_q_coefficients_low_to_high"]), variable
        )
        older_catalog.append((f"rank28-{target['target_label']}", atoms(extension)))
    q103_extension = hasher.extension_key(branch(q103b2["q_coefficients_low_to_high"]), variable)
    older_catalog.append(("q_103b2", atoms(q103_extension)))

    formal_matches = []
    for catalog_label, catalog_key in older_catalog:
        pairs: set[tuple[str, str]] = set()
        for key, label in labels_by_key.items():
            other = key ^ catalog_key
            if other in labels_by_key and label < labels_by_key[other]:
                pairs.add((label, labels_by_key[other]))
        if pairs:
            formal_matches.append({
                "catalog_label": catalog_label,
                "pairs": [list(pair) for pair in sorted(pairs)],
            })

    direct = json.loads(direct_path.read_text())
    splitting = json.loads(SPLITTING.read_text())
    gram = load_gram()
    target_w = [int(value) for value in direct["divisor"]["pinned_trace_vector_w"]]
    norm8_w = [int(value) for value in trace["pinned_rank17_w"]]
    q103_w = [int(value) for value in next(
        record for record in splitting["construction"]["records"]
        if record.get("label") == "norm12-orbit-103b2"
    )["pinned_rank17_w"]]
    rank28_degree = 2 * 3 + 2 * 2 - pairing(target_w, norm8_w, gram)
    q103_degree = 3 * 2 + 2 * 3 - pairing(target_w, q103_w, gram)
    if (rank28_degree, q103_degree) != configuration["older_degrees"]:
        raise ArithmeticError("older-character target-base degree calculation changed")

    result = {
        "schema": "elkies-k3.r17-norm12-complete-character-closure.v1",
        "status": (
            "PASS_EXACT_COMPLETE_THREE_CHARACTER_CLOSURE_FOUND"
            if relations else "PASS_EXACT_NO_COMPLETE_THREE_CHARACTER_CLOSURE"
        ),
        "source_label": arguments.source_label,
        "inputs": {
            relative(path): digest(path)
            for path in (input_path, direct_path, RANK28, SPLITTING, Q103B2, GRAM, HASHER)
        },
        "character_count": len(labels_by_key),
        "distinct_character_count": len(labels_by_key),
        "finite_factor_support_size_histogram": {
            str(size): count for size, count in sorted(factor_count_histogram.items())
        },
        "support_graph_search": search_counts,
        "three_character_relation_count": len(relations),
        "three_character_relations": relations,
        "older_published_base_catalog": {
            "formal_variable_rename_comparison_count": len(older_catalog),
            "formal_product_match_count": sum(len(item["pairs"]) for item in formal_matches),
            "formal_product_matches": formal_matches,
            "compatibility_gate": {
                "target_fibre": configuration["description"],
                "rank28_source_curve_degree_over_target": rank28_degree,
                "q_103b2_source_curve_degree_over_target": q103_degree,
                "both_are_quadratic_target_characters": rank28_degree == q103_degree == 2,
            },
        },
        "proof_boundary": (
            "The support-graph search is exhaustive for q_i*q_j=q_k because every exact "
            "character has one or two finite irreducible polynomial factors and the rational "
            "constant atoms are checked in the full squareclass key. A found relation would "
            "still require the separate V4 base, section-height, and saturation certification. "
            "Formal old-base matches are not target-base rank certificates unless the exact "
            "Neron--Severi degree gate is two."
        ),
        "reproducing_command": (
            ".venv/bin/python elkies-k3/scripts/analyze_r17_norm12_complete_character_closure.py "
            f"--source-label {arguments.source_label} --output {relative(output)}"
        ),
    }
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not output.exists() or output.read_text() != serialized:
            raise ValueError("stored complete-character artifact differs from replay")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
    print(
        "NORM12CLOSURE|source={}|characters={}|relations={}|formal_old_matches={}|status={}".format(
            arguments.source_label, len(labels_by_key), len(relations),
            sum(len(item["pairs"]) for item in formal_matches), result["status"],
        )
    )


if __name__ == "__main__":
    main()
