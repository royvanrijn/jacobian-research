#!/usr/bin/env python3
"""Stream a complete quadratic-cover atlas and search exact character closure.

The equation-level bisection artifacts are intentionally large because they
retain the lifted sections.  Loading a complete atlas merely to inspect its
quadratic branch polynomials can therefore require several gigabytes.  This
checker reads only the exact branch triple, label, orbit mask, and fixed-frame
vector from each record.

Every accepted branch divisor has geometric degree two.  Its squareclass has
one or two irreducible geometric atoms (the point at infinity is an atom when
the polynomial degree is odd).  Consequently a relation q_i*q_j=q_k can only
be an edge with its two singleton endpoints or a triangle of two-atom
supports.  The support search below is exhaustive and the rational constant
factor is tested by an exact rational-square predicate without factoring large
integer contents.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
import json
from math import gcd, isqrt, lcm
from pathlib import Path
from typing import Iterable

import sympy as sp


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PRIORITY = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-08f72-alternate-bisection-priority-v1.tsv"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-08f72-streaming-character-closure-v1.json"
)
RANK28 = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-rank28-genus-one-bisection-pilot-v1.json"
)
Q103B2 = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-norm12-orbit-103b2-twist-section-v1.json"
)

Atom = tuple[object, ...]
Support = frozenset[Atom]


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def relative(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


def parse_array(iterator: Iterable[str]) -> list[object]:
    values = []
    for line in iterator:
        text = line.strip().rstrip(",")
        if text == "]":
            return values
        values.append(json.loads(text))
    raise ValueError("unterminated JSON array")


def stream_records(path: Path):
    """Yield the four small fields needed from every pretty-printed record."""

    fixed_vector = None
    branch = None
    label = None
    with path.open() as stream:
        iterator = iter(stream)
        for line in iterator:
            if '"alternate_rank17_w": [' in line:
                fixed_vector = tuple(map(int, parse_array(iterator)))
            elif (
                '"numerator_coefficients": [' in line
                or '"q_coefficients_low_to_high": [' in line
            ):
                values = parse_array(iterator)
                if len(values) <= 3:
                    branch = tuple(Fraction(str(value)) for value in values)
            elif branch is not None and '"label": ' in line:
                label = str(json.loads(line.split(":", 1)[1].strip().rstrip(",")))
            elif label is not None and '"lattice_orbit_mask": ' in line:
                raw_mask = json.loads(line.split(":", 1)[1].strip().rstrip(","))
                mask = int(raw_mask, 0) if isinstance(raw_mask, str) else int(raw_mask)
                if fixed_vector is None:
                    raise ValueError(f"{path}: {label} has no fixed-frame vector")
                yield {
                    "label": label,
                    "orbit_mask": mask,
                    "fixed_vector": fixed_vector,
                    "branch": branch,
                }
                fixed_vector = None
                branch = None
                label = None
    if any(value is not None for value in (fixed_vector, branch, label)):
        raise ValueError(f"{path}: incomplete final bisection record")


def primitive(
    coefficients: tuple[object, ...], allowed_degrees: tuple[int, ...] | None = (1, 2)
) -> tuple[tuple[int, ...], Fraction]:
    values = [Fraction(str(value)) for value in coefficients]
    while values and values[-1] == 0:
        values.pop()
    if allowed_degrees is not None and len(values) - 1 not in allowed_degrees:
        raise ValueError(f"branch polynomial has degree {len(values) - 1}")
    denominator = lcm(*(value.denominator for value in values))
    integral = [
        value.numerator * (denominator // value.denominator) for value in values
    ]
    content = 0
    for value in integral:
        content = gcd(content, abs(value))
    if not content:
        raise ValueError("zero branch polynomial")
    values = [value // content for value in integral]
    first = next(value for value in values if value)
    scale = Fraction(content, denominator)
    if first < 0:
        values = [-value for value in values]
        scale = -scale
    return tuple(values), scale


def multiply(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    result = [0] * (len(left) + len(right) - 1)
    for i, x_value in enumerate(left):
        for j, y_value in enumerate(right):
            result[i + j] += x_value * y_value
    return tuple(result)


def linear_atom(root: Fraction) -> Atom:
    # denominator*u-numerator, normalized by its first nonzero coefficient.
    polynomial, unused_scale = primitive((-root.numerator, root.denominator))
    return ("L", *polynomial)


def branch_character(coefficients: tuple[int, ...]) -> tuple[Support, Fraction, int, bool]:
    polynomial, scale = primitive(coefficients)
    degree = len(polynomial) - 1
    atoms: list[Atom] = []
    irreducible = False
    if degree == 1:
        atoms.append(("L", *polynomial))
        atoms.append(("I",))
    else:
        constant, linear, quadratic = polynomial
        discriminant = linear * linear - 4 * quadratic * constant
        if discriminant == 0:
            raise ValueError("quadratic branch polynomial is not squarefree")
        root_discriminant = isqrt(discriminant) if discriminant > 0 else -1
        if root_discriminant >= 0 and root_discriminant**2 == discriminant:
            roots = {
                Fraction(-linear + root_discriminant, 2 * quadratic),
                Fraction(-linear - root_discriminant, 2 * quadratic),
            }
            if len(roots) != 2:
                raise ValueError("quadratic branch polynomial has a repeated root")
            atoms.extend(linear_atom(root) for root in sorted(roots))
        else:
            irreducible = True
            atoms.append(("Q", *polynomial))

    finite_polynomial = (1,)
    for atom in atoms:
        if atom[0] != "I":
            finite_polynomial = multiply(finite_polynomial, tuple(map(int, atom[1:])))
    ratio = None
    for actual, reconstructed in zip(polynomial, finite_polynomial):
        if reconstructed:
            ratio = Fraction(actual, reconstructed)
            break
        if actual:
            raise ArithmeticError("factor support does not reconstruct the branch")
    if ratio is None or any(
        Fraction(actual) != ratio * reconstructed
        for actual, reconstructed in zip(polynomial, finite_polynomial)
    ):
        raise ArithmeticError("factor support is not proportional to the branch")
    scalar = Fraction(scale) * ratio
    return frozenset(atoms), scalar, degree, irreducible


def catalog_character(coefficients: list[object]) -> tuple[Support, Fraction]:
    """Factor one committed squarefree character without factoring its content."""

    values = [Fraction(str(value)) for value in coefficients]
    while values and not values[-1]:
        values.pop()
    if not values:
        raise ValueError("zero catalog character")
    variable = sp.Symbol("u")
    expression = sum(
        sp.Rational(value.numerator, value.denominator) * variable**index
        for index, value in enumerate(values)
    )
    coefficient, factors = sp.factor_list(expression, variable)
    scalar = Fraction(int(coefficient.p), int(coefficient.q))
    atoms = set()
    factored_degree = 0
    for factor, exponent in factors:
        if exponent != 1:
            raise ValueError("committed catalog character is not squarefree")
        polynomial = sp.Poly(factor, variable, domain=sp.QQ)
        rational_coefficients = [
            Fraction(int(value.p), int(value.q))
            for value in reversed(polynomial.all_coeffs())
        ]
        denominator = lcm(*(value.denominator for value in rational_coefficients))
        integral = tuple(
            value.numerator * (denominator // value.denominator)
            for value in rational_coefficients
        )
        normalized, factor_scale = primitive(integral, allowed_degrees=None)
        scalar *= factor_scale / denominator
        degree = len(normalized) - 1
        atom_type = "L" if degree == 1 else "Q" if degree == 2 else "P"
        atoms.add((atom_type, *normalized))
        factored_degree += degree
    if factored_degree != len(values) - 1:
        raise ArithmeticError("catalog factor degrees do not reconstruct the character")
    if factored_degree % 2:
        atoms.add(("I",))
    return frozenset(atoms), scalar


def rational_square(value: Fraction) -> bool:
    if value <= 0:
        return False
    numerator_root = isqrt(value.numerator)
    denominator_root = isqrt(value.denominator)
    return (
        numerator_root**2 == value.numerator
        and denominator_root**2 == value.denominator
    )


def scalar_text(value: Fraction) -> str:
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"{value.numerator}/{value.denominator}"
    )


def atom_payload(atom: Atom) -> list[object]:
    return list(atom)


def support_payload(support: Support) -> list[list[object]]:
    return [atom_payload(atom) for atom in sorted(support)]


def priority_rows(path: Path) -> dict[int, tuple[int, ...]]:
    result = {}
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            mask = int(row["orbit_mask"])
            if mask in result:
                raise ValueError(f"{path}: duplicate orbit mask {mask}")
            result[mask] = tuple(map(int, row["historical_alternate_w"].split()))
    return result


def collision_groups(records: list[dict]) -> list[list[dict]]:
    groups: list[list[dict]] = []
    for record in records:
        for group in groups:
            if rational_square(record["scalar"] / group[0]["scalar"]):
                group.append(record)
                break
        else:
            groups.append([record])
    return [group for group in groups if len(group) > 1]


def relation_search(by_support: dict[Support, list[dict]]):
    singleton = {
        next(iter(support)): support
        for support in by_support
        if len(support) == 1
    }
    edges = {support for support in by_support if len(support) == 2}
    relations: set[tuple[str, str, str]] = set()
    endpoint_candidates = 0
    for edge in edges:
        left_atom, right_atom = tuple(edge)
        left_support = singleton.get(left_atom)
        right_support = singleton.get(right_atom)
        if left_support is None or right_support is None:
            continue
        for edge_record in by_support[edge]:
            for left_record in by_support[left_support]:
                for right_record in by_support[right_support]:
                    endpoint_candidates += 1
                    labels = {
                        edge_record["label"], left_record["label"], right_record["label"]
                    }
                    if len(labels) == 3 and rational_square(
                        edge_record["scalar"]
                        * left_record["scalar"]
                        * right_record["scalar"]
                    ):
                        relations.add(tuple(sorted(labels)))

    adjacency: dict[Atom, set[Atom]] = defaultdict(set)
    for edge in edges:
        left_atom, right_atom = tuple(edge)
        adjacency[left_atom].add(right_atom)
        adjacency[right_atom].add(left_atom)
    triangle_supports = 0
    triangle_candidates = 0
    for left_atom in sorted(adjacency):
        for right_atom in sorted(
            atom for atom in adjacency[left_atom] if left_atom < atom
        ):
            for third_atom in sorted(
                atom
                for atom in adjacency[left_atom].intersection(adjacency[right_atom])
                if right_atom < atom
            ):
                triangle_supports += 1
                supports = (
                    frozenset((left_atom, right_atom)),
                    frozenset((left_atom, third_atom)),
                    frozenset((right_atom, third_atom)),
                )
                for first in by_support[supports[0]]:
                    for second in by_support[supports[1]]:
                        for third in by_support[supports[2]]:
                            triangle_candidates += 1
                            labels = {first["label"], second["label"], third["label"]}
                            if len(labels) == 3 and rational_square(
                                first["scalar"] * second["scalar"] * third["scalar"]
                            ):
                                relations.add(tuple(sorted(labels)))
    return sorted(relations), {
        "singleton_support_count": len(singleton),
        "two_atom_support_count": len(edges),
        "endpoint_scalar_candidate_count": endpoint_candidates,
        "triangle_support_count": triangle_supports,
        "triangle_scalar_candidate_count": triangle_candidates,
    }


def committed_catalog():
    rank28 = json.loads(RANK28.read_text())
    q103b2 = json.loads(Q103B2.read_text())
    result = []
    for target in rank28["traces"][0]["targets"]:
        support, scalar = catalog_character(
            target["branch_polynomial_q_coefficients_low_to_high"]
        )
        result.append((f"rank28-{target['target_label']}", support, scalar))
    support, scalar = catalog_character(q103b2["q_coefficients_low_to_high"])
    result.append(("q_103b2", support, scalar))
    return result


def catalog_product_matches(by_support: dict[Support, list[dict]], catalog):
    matches = []
    for catalog_label, target_support, target_scalar in catalog:
        pairs = set()
        for left_support, left_records in by_support.items():
            right_support = left_support ^ target_support
            right_records = by_support.get(right_support)
            if right_records is None:
                continue
            for left in left_records:
                for right in right_records:
                    if left["label"] >= right["label"]:
                        continue
                    if rational_square(
                        left["scalar"] * right["scalar"] / target_scalar
                    ):
                        pairs.add((left["label"], right["label"]))
        if pairs:
            matches.append(
                {
                    "catalog_label": catalog_label,
                    "pairs": [list(pair) for pair in sorted(pairs)],
                }
            )
    return matches


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, action="append", required=True)
    parser.add_argument("--priority-table", type=Path)
    parser.add_argument("--source-label", default="norm12-orbit-08f72")
    parser.add_argument("--expected-count", type=int, default=39147)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    if arguments.expected_count <= 0:
        parser.error("--expected-count must be positive")

    input_paths = [path.resolve() for path in arguments.input]
    priority_path = (
        None if arguments.priority_table is None else arguments.priority_table.resolve()
    )
    expected = {} if priority_path is None else priority_rows(priority_path)
    if priority_path is not None and len(expected) != arguments.expected_count:
        raise ValueError(
            f"priority table has {len(expected)} rows, expected {arguments.expected_count}"
        )

    seen_labels = set()
    seen_masks = set()
    by_support: dict[Support, list[dict]] = defaultdict(list)
    degree_histogram = Counter()
    irreducible_count = 0
    for path in input_paths:
        for raw in stream_records(path):
            label = raw["label"]
            mask = raw["orbit_mask"]
            if label in seen_labels or mask in seen_masks:
                raise ValueError(f"duplicate label or orbit mask at {label}")
            if priority_path is not None and (
                mask not in expected or raw["fixed_vector"] != expected[mask]
            ):
                raise ValueError(f"{label}: priority-table orbit attachment mismatch")
            support, scalar, degree, irreducible = branch_character(raw["branch"])
            if len(support) not in (1, 2):
                raise ArithmeticError(f"{label}: unexpected branch support size")
            record = {
                "label": label,
                "orbit_mask": mask,
                "fixed_vector": raw["fixed_vector"],
                "scalar": scalar,
            }
            by_support[support].append(record)
            seen_labels.add(label)
            seen_masks.add(mask)
            degree_histogram[degree] += 1
            irreducible_count += int(irreducible)

    if len(seen_masks) != arguments.expected_count:
        raise ValueError(
            f"atlas has {len(seen_masks)} distinct masks, expected {arguments.expected_count}"
        )
    if priority_path is not None and seen_masks != set(expected):
        missing = sorted(set(expected) - seen_masks)
        extra = sorted(seen_masks - set(expected))
        raise ValueError(
            f"atlas is not complete: missing={missing[:10]}, extra={extra[:10]}"
        )

    collisions = []
    for support, records in by_support.items():
        for group in collision_groups(records):
            collisions.append(
                {
                    "support": support_payload(support),
                    "records": [
                        {
                            "label": record["label"],
                            "orbit_mask": record["orbit_mask"],
                            "fixed_frame_vector": list(record["fixed_vector"]),
                            "scalar_representative": scalar_text(record["scalar"]),
                        }
                        for record in group
                    ],
                }
            )
    collisions.sort(key=lambda item: [record["label"] for record in item["records"]])

    relations, relation_counts = relation_search(by_support)
    catalog = committed_catalog()
    catalog_matches = catalog_product_matches(by_support, catalog)
    catalog_match_count = sum(len(item["pairs"]) for item in catalog_matches)
    status_parts = []
    if collisions:
        status_parts.append("EQUAL_COVER_COLLISIONS_FOUND")
    if relations:
        status_parts.append("THREE_CHARACTER_CLOSURE_FOUND")
    if catalog_match_count:
        status_parts.append("COMMITTED_CATALOG_PRODUCT_MATCH_FOUND")
    status = (
        "PASS_EXACT_" + "_AND_".join(status_parts)
        if status_parts
        else "PASS_EXACT_NO_EQUAL_COVER_OR_THREE_CHARACTER_CLOSURE"
    )

    output_path = arguments.output.resolve()
    inputs = input_paths + ([] if priority_path is None else [priority_path]) + [
        RANK28,
        Q103B2,
    ]
    payload = {
        "schema": "elkies-k3.r17-norm12-streaming-quadratic-character-closure.v1",
        "status": status,
        "source_label": arguments.source_label,
        "character_count": len(seen_labels),
        "priority_orbit_count": len(expected),
        "distinct_polynomial_support_count": len(by_support),
        "branch_polynomial_degree_histogram": {
            str(key): value for key, value in sorted(degree_histogram.items())
        },
        "irreducible_quadratic_count": irreducible_count,
        "equal_cover_collision_count": len(collisions),
        "equal_cover_collisions": collisions,
        "support_graph_search": relation_counts,
        "three_character_relation_count": len(relations),
        "three_character_relations": relations,
        "committed_character_catalog": {
            "comparison_count": len(catalog),
            "formal_variable_rename_product_match_count": catalog_match_count,
            "formal_variable_rename_product_matches": catalog_matches,
        },
        "inputs": {relative(path): digest(path) for path in inputs},
        "proof_boundary": (
            "Every supplied exact degree-two branch record is streamed, with complete "
            "orbit coverage checked against the declared count"
            + (
                " and priority table. "
                if priority_path is not None
                else ". "
            )
            + "Equal covers are grouped by identical geometric support "
            "and exact rational-square scalar ratio. Because every support has size one or "
            "two, the endpoint and triangle search exhausts q_i*q_j=q_k while retaining the "
            "rational constant squareclass. Every pair product is also compared exactly "
            "against the eleven committed rank-28 characters and q_103b2 by support xor "
            "and a factorless rational-square test. Those catalog comparisons are formal "
            "variable renames until a target-base degree-two compatibility gate is proved. "
            "A collision or relation still requires the separate anti-invariant height and "
            "base certification."
        ),
        "reproducing_command": (
            ".venv/bin/python "
            "elkies-k3/scripts/search_r17_norm12_quadratic_character_closure_streaming.py "
            + f"--source-label {arguments.source_label} "
            + (
                ""
                if priority_path is None
                else f"--priority-table {relative(priority_path)} "
            )
            + " ".join(f"--input {relative(path)}" for path in input_paths)
            + f" --output {relative(output_path)}"
        ),
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not output_path.exists() or output_path.read_text() != serialized:
            raise ValueError("stored streaming closure artifact differs from replay")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(
        "R17STREAMCHAR"
        f"|source={arguments.source_label}|characters={len(seen_labels)}"
        f"|supports={len(by_support)}|collisions={len(collisions)}"
        f"|relations={len(relations)}|catalog_matches={catalog_match_count}|status={status}"
        f"|output={relative(output_path)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
