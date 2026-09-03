#!/usr/bin/env sage-python
"""Certify ADE mass strata as asymptotic Kneser-neighbour scores.

This is a bounded calibration, not the missing local-density inversion.  It
uses the complete determinant-78 census to aggregate all realized ADE strata,
computes the representation rows of rank at most four, and solves exact
nonnegative rational LPs with matching primal and dual certificates.  It also
records the exact determinant-948 rootless mass, a determinant-950 rootless
mass lower bound, the genus masses, and the first feasible finite-p line
counts for the three rank-17 controls.

The script deliberately fails closed: the rank-at-most-four determinant-78
LP has a positive optimum, so it does not claim to recover mu_0=0 from local
moments.  Nor does it call the census-realized ADE list "locally admissible".
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

from sage.all import (
    CartanMatrix,
    MixedIntegerLinearProgram,
    QQ,
    RootSystem,
    ZZ,
    matrix,
    vector,
)
from sage.env import SAGE_VERSION
from sage.quadratic_forms.genera.genus import Genus


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = Path(__file__).resolve()
DET78_GRAM = ROOT / "elkies-k3/data/lattice/e6_rank4_det78_frame.txt"
DET948_GRAM = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
DET78_CENSUS = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-e6-rank4-det78-niemeier-frames-v1.json"
)
DET948_ROOTLESS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-rootless-j2-niemeier-first.json"
)
DET950_FOUNDRY = (
    ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-v1.json"
)
DET950_NEW = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-ns0024-new-rootless-source-route-v1.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-ade-neighbor-mass-score-v1.json"
)

RANK_THREE_MOMENT_TYPES = ("A1", "2A1", "A2", "3A1", "A1+A2", "A3")
RANK_FOUR_MOMENT_TYPES = ("4A1", "2A1+A2", "2A2", "A1+A3", "A4", "D4")
MOMENT_TYPES = RANK_THREE_MOMENT_TYPES + RANK_FOUR_MOMENT_TYPES
TYPE_PATTERN = re.compile(r"(?:(\d+))?([ADE])(\d+)")
RELATION_CACHE = {}


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rational(value):
    value = QQ(value)
    return {
        "denominator": int(value.denominator()),
        "numerator": int(value.numerator()),
        "text": str(value),
    }


def load_matrix(path):
    return matrix(
        ZZ,
        [
            [ZZ(entry) for entry in line.split()]
            for line in Path(path).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def one_proper_spinor_genus(genus):
    ambient, kernel = genus._proper_spinor_kernel()
    count = ZZ(ambient.quotient(kernel).order())
    assert count == 1
    return int(count)


def component_list(root_type):
    if root_type == "0":
        return []
    answer = []
    for term in root_type.split("+"):
        match = TYPE_PATTERN.fullmatch(term)
        assert match is not None, term
        multiplicity = int(match.group(1) or 1)
        answer.extend(
            [(match.group(2), int(match.group(3)))] * multiplicity
        )
    return answer


def canonical_root_type(components):
    counts = Counter(components)
    terms = []
    for component in sorted(counts, key=lambda item: (item[0], item[1])):
        multiplicity = counts[component]
        prefix = str(multiplicity) if multiplicity > 1 else ""
        terms.append(f"{prefix}{component[0]}{component[1]}")
    return "+".join(terms) if terms else "0"


def abstract_ade_types(maximum_rank=17):
    """Return every abstract ADE root type of rank at most maximum_rank."""

    irreducibles = [("A", rank) for rank in range(1, maximum_rank + 1)]
    irreducibles += [("D", rank) for rank in range(4, maximum_rank + 1)]
    irreducibles += [("E", rank) for rank in (6, 7, 8)]
    states = {(0, ())}
    for component in irreducibles:
        next_states = set(states)
        component_rank = component[1]
        for old_rank, old_components in states:
            new_components = list(old_components)
            new_rank = old_rank
            while new_rank + component_rank <= maximum_rank:
                new_components.append(component)
                new_rank += component_rank
                next_states.add((new_rank, tuple(new_components)))
        states = next_states
    rows = [
        {
            "rank": rank,
            "root_type": canonical_root_type(components),
        }
        for rank, components in states
    ]
    rows.sort(key=lambda row: (row["rank"], row["root_type"]))
    assert len(rows) == 3768
    return rows


def irreducible_relations(component):
    """Return bit sets for orthogonal and inner-product-minus-one roots."""

    if component in RELATION_CACHE:
        return RELATION_CACHE[component]
    root_lattice = RootSystem(list(component)).root_lattice()
    cartan = CartanMatrix(list(component))
    roots = [vector(root.to_vector()) for root in root_lattice.roots()]
    orthogonal = []
    minus_one = []
    for left in roots:
        orthogonal_mask = 0
        minus_one_mask = 0
        for index, right in enumerate(roots):
            pairing = left * cartan * right
            if pairing == 0:
                orthogonal_mask |= 1 << index
            elif pairing == -1:
                minus_one_mask |= 1 << index
        orthogonal.append(orthogonal_mask)
        minus_one.append(minus_one_mask)
    RELATION_CACHE[component] = (orthogonal, minus_one)
    return RELATION_CACHE[component]


def target_relations(root_type):
    """Return the root masks and irreducible blocks of an ADE root system."""

    components = component_list(root_type)
    pieces = [irreducible_relations(component) for component in components]
    sizes = [len(piece[0]) for piece in pieces]
    root_count = sum(sizes)
    all_roots = (1 << root_count) - 1
    orthogonal = []
    minus_one = []
    blocks = []
    offset = 0
    for (local_orthogonal, local_minus_one), size in zip(pieces, sizes):
        local_block = ((1 << size) - 1) << offset
        other_components = all_roots ^ local_block
        blocks.append((offset, size, local_block))
        orthogonal.extend(
            other_components | (mask << offset) for mask in local_orthogonal
        )
        minus_one.extend(mask << offset for mask in local_minus_one)
        offset += size

    return all_roots, orthogonal, minus_one, blocks


def rank_three_counts_on_mask(allowed, orthogonal, minus_one):
    """Count the rank-at-most-three embeddings inside ``allowed`` roots."""

    root_count = allowed.bit_count()
    two_a1 = 0
    a2 = 0
    three_a1 = 0
    a1_a2 = 0
    a3 = 0
    first_candidates = allowed
    while first_candidates:
        first_bit = first_candidates & -first_candidates
        first_candidates -= first_bit
        first = first_bit.bit_length() - 1
        orthogonal_to_first = orthogonal[first] & allowed
        minus_one_to_first = minus_one[first] & allowed
        two_a1 += orthogonal_to_first.bit_count()
        a2 += minus_one_to_first.bit_count()
        candidates = orthogonal_to_first
        while candidates:
            bit = candidates & -candidates
            second = bit.bit_length() - 1
            candidates -= bit
            three_a1 += (
                orthogonal_to_first & orthogonal[second] & allowed
            ).bit_count()
        candidates = minus_one_to_first
        while candidates:
            bit = candidates & -candidates
            second = bit.bit_length() - 1
            candidates -= bit
            a1_a2 += (
                orthogonal[first] & orthogonal[second] & allowed
            ).bit_count()
            a3 += (
                minus_one[first] & orthogonal[second] & allowed
            ).bit_count()

    answer = (root_count, two_a1, a2, three_a1, a1_a2, a3)
    return dict(zip(RANK_THREE_MOMENT_TYPES, answer))


def low_rank_embedding_counts(root_type):
    """Count embeddings of every ADE root lattice of rank at most four."""

    all_roots, orthogonal, minus_one, blocks = target_relations(root_type)
    answer = rank_three_counts_on_mask(all_roots, orthogonal, minus_one)

    rank_four = {name: 0 for name in RANK_FOUR_MOMENT_TYPES}

    # An irreducible simply-laced Weyl group is transitive on its roots and
    # on ordered root pairs of inner product -1.  Fix one representative in
    # each irreducible target block, count in its exact orthogonal mask, and
    # multiply by the orbit size.  This avoids materializing huge tuple sets.
    for offset, size, block in blocks:
        representative = offset
        complement = rank_three_counts_on_mask(
            orthogonal[representative], orthogonal, minus_one
        )
        rank_four["4A1"] += size * complement["3A1"]
        rank_four["2A1+A2"] += size * complement["A1+A2"]
        rank_four["A1+A3"] += size * complement["A3"]

        adjacent = minus_one[representative] & block
        if adjacent:
            second_bit = adjacent & -adjacent
            second = second_bit.bit_length() - 1
            ordered_a2_in_block = sum(
                (minus_one[index] & block).bit_count()
                for index in range(offset, offset + size)
            )
            double_complement = rank_three_counts_on_mask(
                orthogonal[representative] & orthogonal[second],
                orthogonal,
                minus_one,
            )
            rank_four["2A2"] += (
                ordered_a2_in_block * double_complement["A2"]
            )

    # Connected source diagrams must land in one irreducible target block.
    # Count A4 by fixing its second node; count D4 by fixing its centre.
    for center in range(all_roots.bit_count()):
        first_candidates = minus_one[center]
        while first_candidates:
            first_bit = first_candidates & -first_candidates
            first_candidates -= first_bit
            first = first_bit.bit_length() - 1
            third_candidates = minus_one[center] & orthogonal[first]
            while third_candidates:
                third_bit = third_candidates & -third_candidates
                third_candidates -= third_bit
                third = third_bit.bit_length() - 1
                rank_four["A4"] += (
                    minus_one[third]
                    & orthogonal[center]
                    & orthogonal[first]
                ).bit_count()
        rank_four["D4"] += rank_three_counts_on_mask(
            minus_one[center], orthogonal, minus_one
        )["3A1"]

    answer.update(rank_four)
    return answer


def exact_lp(moment_names, types, normalized_moments):
    """Maximize the rootless fraction and emit exact primal/dual witnesses."""

    labels = ["0"] + sorted(types)
    columns = {
        label: [QQ(1)]
        + [QQ(0) if label == "0" else QQ(types[label][name]) for name in moment_names]
        for label in labels
    }
    right_hand_side = [QQ(1)] + [QQ(normalized_moments[name]) for name in moment_names]

    primal = MixedIntegerLinearProgram(maximization=True, solver="PPL")
    variables = primal.new_variable(nonnegative=True)
    for row in range(len(right_hand_side)):
        primal.add_constraint(
            sum(columns[label][row] * variables[index] for index, label in enumerate(labels))
            == right_hand_side[row]
        )
    primal.set_objective(variables[0])
    primal_optimum = QQ(primal.solve())
    primal_values = {
        label: QQ(primal.get_values(variables[index]))
        for index, label in enumerate(labels)
        if QQ(primal.get_values(variables[index])) != 0
    }
    assert primal_values["0"] == primal_optimum
    for row, expected in enumerate(right_hand_side):
        assert sum(
            columns[label][row] * value for label, value in primal_values.items()
        ) == expected

    dual = MixedIntegerLinearProgram(maximization=False, solver="PPL")
    positive = dual.new_variable(nonnegative=True)
    negative = dual.new_variable(nonnegative=True)
    dual_variables = [positive[row] - negative[row] for row in range(len(right_hand_side))]
    for label in labels:
        target = QQ(1) if label == "0" else QQ(0)
        dual.add_constraint(
            sum(
                columns[label][row] * dual_variables[row]
                for row in range(len(right_hand_side))
            )
            >= target
        )
    dual.set_objective(
        sum(right_hand_side[row] * dual_variables[row] for row in range(len(right_hand_side)))
    )
    dual_optimum = QQ(dual.solve())
    dual_values = [
        QQ(dual.get_values(positive[row])) - QQ(dual.get_values(negative[row]))
        for row in range(len(right_hand_side))
    ]
    assert dual_optimum == primal_optimum
    assert sum(a * b for a, b in zip(right_hand_side, dual_values)) == primal_optimum
    for label in labels:
        target = QQ(1) if label == "0" else QQ(0)
        assert sum(a * b for a, b in zip(columns[label], dual_values)) >= target

    row_labels = ["total_mass"] + list(moment_names)
    return {
        "candidate_root_system_count_including_zero": len(labels),
        "conclusion": (
            "inconclusive_positive_rootless_fraction_allowed"
            if primal_optimum > 0
            else "rootless_mass_forced_zero"
        ),
        "dual_certificate": {
            label: rational(value) for label, value in zip(row_labels, dual_values)
        },
        "maximum_rootless_fraction": rational(primal_optimum),
        "moment_rows": list(moment_names),
        "primal_certificate": {
            label: rational(value) for label, value in sorted(primal_values.items())
        },
        "strong_duality_checked_exactly": True,
    }


def genus_record(gram):
    genus = Genus(gram)
    return genus, {
        "determinant": int(gram.det()),
        "mass": rational(genus.mass()),
        "proper_spinor_genus_count": one_proper_spinor_genus(genus),
        "rank": gram.nrows(),
    }


def odd_rank_neighbor_count(rank, prime):
    assert rank % 2 == 1
    return sum(ZZ(prime) ** exponent for exponent in range(rank - 1))


def build():
    det78_census = json.loads(DET78_CENSUS.read_text())
    det948_rootless = json.loads(DET948_ROOTLESS.read_text())
    det950_foundry = json.loads(DET950_FOUNDRY.read_text())
    det950_new = json.loads(DET950_NEW.read_text())

    det78 = load_matrix(DET78_GRAM)
    det948 = load_matrix(DET948_GRAM)
    det950 = matrix(ZZ, det950_new["new_rootless_frame"]["gram"])
    genus78, genus_record78 = genus_record(det78)
    genus948, genus_record948 = genus_record(det948)
    genus950, genus_record950 = genus_record(det950)

    assert genus_record78["determinant"] == 78
    assert genus_record948["determinant"] == 948
    assert genus_record950["determinant"] == 950
    assert genus78.mass() == QQ(det78_census["accounting"]["target_genus_mass"])

    type_masses = defaultdict(lambda: QQ(0))
    type_class_counts = Counter()
    embedding_counts = {}
    moments = {name: QQ(0) for name in MOMENT_TYPES}
    census_mass = QQ(0)
    for frame in det78_census["frames"]:
        root_type = frame["root_type"]
        if root_type not in embedding_counts:
            embedding_counts[root_type] = low_rank_embedding_counts(root_type)
        assert embedding_counts[root_type]["A1"] == frame["signed_root_count"]
        weight = QQ(1) / ZZ(frame["automorphism_group_order"])
        census_mass += weight
        type_masses[root_type] += weight
        type_class_counts[root_type] += 1
        for name in MOMENT_TYPES:
            moments[name] += weight * embedding_counts[root_type][name]
    assert census_mass == genus78.mass()
    assert len(det78_census["frames"]) == 1549
    assert len(type_masses) == 621
    normalized_moments = {name: value / census_mass for name, value in moments.items()}

    diagonal_orders = {
        "A1": 2,
        "2A1": 8,
        "A2": 12,
        "3A1": 48,
        "A1+A2": 24,
        "A3": 48,
        "4A1": 384,
        "2A1+A2": 96,
        "2A2": 288,
        "A1+A3": 96,
        "A4": 240,
        "D4": 1152,
    }
    for name in MOMENT_TYPES:
        assert low_rank_embedding_counts(name)[name] == diagonal_orders[name]

    ade_rows = abstract_ade_types()
    realized = set(type_masses)
    assert realized.issubset({row["root_type"] for row in ade_rows})

    det78_strata = [
        {
            "class_count": type_class_counts[root_type],
            "low_rank_embedding_counts": embedding_counts[root_type],
            "mass": rational(type_masses[root_type]),
            "root_type": root_type,
        }
        for root_type in sorted(type_masses)
    ]

    lp_rows = [
        exact_lp(("A1",), embedding_counts, normalized_moments),
        exact_lp(("A1", "2A1", "A2"), embedding_counts, normalized_moments),
        exact_lp(RANK_THREE_MOMENT_TYPES, embedding_counts, normalized_moments),
        exact_lp(MOMENT_TYPES, embedding_counts, normalized_moments),
    ]
    assert all(QQ(row["maximum_rootless_fraction"]["text"]) > 0 for row in lp_rows)

    rootless948 = sum(
        (
            QQ(1) / ZZ(row["automorphism_group_order"])
            for row in det948_rootless["rootless_classes"]
        ),
        QQ(0),
    )
    assert len(det948_rootless["rootless_classes"]) == 2
    assert rootless948 == QQ(3) / 4

    ns0024 = next(row for row in det950_foundry["ns_classes"] if row["ns_id"] == "NS0024")
    catalog950 = [
        row for row in ns0024["frames"] if row["root_type"] == "0"
    ]
    assert {row["frame_id"] for row in catalog950} == {
        "NS0024-F001",
        "NS0024-F002",
        "NS0024-F005",
    }
    known950_orders = [
        ZZ(row["rootless_intrinsics"]["automorphism_group_order"])
        for row in catalog950
    ]
    known950_orders.append(ZZ(det950_new["new_rootless_frame"]["automorphism_group_order"]))
    known950_norm_four = [
        ZZ(row["rootless_intrinsics"]["norm_four_vectors"]) for row in catalog950
    ]
    known950_norm_four.append(ZZ(det950_new["new_rootless_frame"]["norm_four_vectors"]))
    assert len(set(known950_norm_four)) == 4
    rootless950_lower = sum((QQ(1) / order for order in known950_orders), QQ(0))
    assert rootless950_lower == QQ(7) / 4

    fraction948 = rootless948 / genus948.mass()
    fraction950_lower = rootless950_lower / genus950.mass()
    return {
        "abstract_ADE_rank_at_most_17": {
            "complete_abstract_list": ade_rows,
            "count_including_zero": len(ade_rows),
            "local_admissibility_status": (
                "NOT_COMPUTED: this complete abstract list is a certified superset; "
                "integral local embedding filters at 2 and determinant primes remain required"
            ),
        },
        "controls": {
            "determinant_78": {
                **genus_record78,
                "exact_rootless_mass": rational(0),
                "exact_rootless_neighbor_fraction": rational(0),
                "first_good_prime_test": {
                    "exact_rootless_lines": 0,
                    "prime": 5,
                    "total_isotropic_lines": int(odd_rank_neighbor_count(17, 5)),
                    "reason": "the complete genus census is rootful",
                },
                "realized_ADE_strata": det78_strata,
                "realized_ADE_type_count": len(det78_strata),
                "rootless_mass_source": "complete 1549-class mass-closed census",
            },
            "determinant_948": {
                **genus_record948,
                "asymptotic_rootless_neighbor_fraction": rational(fraction948),
                "exact_rootless_mass": rational(rootless948),
                "first_good_prime_line_count": {
                    "exact_frequency_status": "NOT_COMPUTED",
                    "prime": 5,
                    "total_isotropic_lines": int(odd_rank_neighbor_count(17, 5)),
                },
                "rootless_mass_source": "complete two-class rootless Niemeier classification",
            },
            "determinant_950": {
                **genus_record950,
                "asymptotic_rootless_neighbor_fraction_lower_bound": rational(fraction950_lower),
                "exact_rootless_mass_status": "NOT_COMPUTED",
                "first_good_prime_line_count": {
                    "exact_frequency_status": "NOT_COMPUTED",
                    "prime": 3,
                    "total_isotropic_lines": int(odd_rank_neighbor_count(17, 3)),
                },
                "known_distinct_rootless_classes": len(known950_orders),
                "known_rootless_automorphism_orders": list(map(int, known950_orders)),
                "known_rootless_mass_lower_bound": rational(rootless950_lower),
                "known_rootless_norm_four_counts": list(map(int, known950_norm_four)),
                "rootless_mass_source": "three foundry classes plus one distinct exact route class",
            },
        },
        "determinant_78_low_rank_LP": {
            "candidate_set": (
                "the 621 census-realized nonzero ADE types plus zero; this is an "
                "optimistically pruned subset of any valid locally admissible list"
            ),
            "exact_LPs": lp_rows,
            "moment_rows": {
                name: {
                    "normalized_average": rational(normalized_moments[name]),
                    "weighted_total": rational(moments[name]),
                }
                for name in MOMENT_TYPES
            },
            "result": (
                "Rank-at-most-four ADE averages do not force mu_0=0, even after "
                "using the complete census to remove every unrealized nonzero type. "
                "Substantially higher ADE rows or full triangular inversion are required."
            ),
            "true_census_rootless_mass": rational(0),
        },
        "inputs": {
            relative(path): "sha256:" + digest(path)
            for path in (
                DET78_GRAM,
                DET948_GRAM,
                DET78_CENSUS,
                DET948_ROOTLESS,
                DET950_FOUNDRY,
                DET950_NEW,
            )
        },
        "proof_boundary": (
            "The asymptotic conversion from stratum mass to neighbor frequency is a "
            "theorem of Chenevier after fixing a compatible spinor displacement. "
            "This artifact supplies exact mass inputs and a low-rank LP calibration. "
            "It does not compute the complete locally admissible ADE list, higher local "
            "representation densities, determinant-950 exact rootless mass, or the "
            "determinant-948/950 finite-p rootless-line frequencies."
        ),
        "reproduce": {
            "check": "sage -python elkies-k3/scripts/certify_ade_neighbor_mass_score.sage --check",
            "generate": "sage -python elkies-k3/scripts/certify_ade_neighbor_mass_score.sage",
        },
        "sage_version": SAGE_VERSION,
        "schema": "elkies-k3.ade-neighbor-mass-score.v1",
        "script": relative(SCRIPT),
        "status": "PASS_EXACT_ADE_MASS_SCORE_CALIBRATION_WITH_FAIL_CLOSED_BOUNDARY",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    rendered = json.dumps(build(), indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != rendered:
            raise SystemExit(f"stale artifact: run {relative(SCRIPT)}")
        print(f"OK: {relative(output)}")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(rendered)
    print(f"wrote {relative(output)}")


if __name__ == "__main__":
    main()
