#!/usr/bin/env python3
"""Conditional permutation gluing for the F2 terminal A6 packet.

The certified input is the degree-six branch-cycle passport

    (5,1) | (3,3) | (3,1,1,1)

with product one and monodromy A6.  The Kummer-orbit theorem also leaves two
simple cofactor orbits on the squarefree-R stratum, but it does *not* prove a
branch-cycle model for those orbits.

This checker therefore audits the following explicit conditional model:

* the A6 triple fixes every sheet outside its six-sheet packet;
* each of the two simple cofactor orbits supplies one additional simple
  branch value, hence one transposition;
* the five branch cycles can be ordered as the A6 triple followed by the two
  transpositions, with no other branch cycles on that compact target curve.

It exhaustively enumerates that model.  If at least one spectator sheet is
required, product one, transitivity, and Riemann--Hurwitz force degree seven.
There are 30 normalized tuples after fixing the core letters and five-cycle,
and six simultaneous-conjugacy classes; every tuple has genus zero and
generates S7.  The certified endpoint/interior markings give a further
conditional filter: requiring the connector anchor to avoid both source
endpoints leaves three classes of signature (5,3,1).  Thus even the marked
conditional model survives.

The checker also constructs paired-star witnesses in every remaining degree.
They show that an A6 packet plus 2*k simple cycles is compatible with a
connected genus-zero degree-(6+k) cover.  This witness family is exact
permutation theory, but its identification with F2 boundary components is
not proved and is not an exclusion of (75,125).
"""

from __future__ import annotations

import argparse
from itertools import combinations, permutations
import json
from pathlib import Path
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/jc2_f2_a6_simple_spectator_gluing.json"
)

Permutation = tuple[int, ...]
BranchTuple = tuple[Permutation, ...]


def identity(degree: int) -> Permutation:
    return tuple(range(degree))


def compose(left: Permutation, right: Permutation) -> Permutation:
    """Return ``left`` after ``right``."""

    return tuple(left[right[index]] for index in range(len(left)))


def inverse(permutation: Permutation) -> Permutation:
    result = [0] * len(permutation)
    for index, image in enumerate(permutation):
        result[image] = index
    return tuple(result)


def conjugate(element: Permutation, permutation: Permutation) -> Permutation:
    return compose(compose(element, permutation), inverse(element))


def permutation_product(branch_cycles: Sequence[Permutation]) -> Permutation:
    if not branch_cycles:
        raise ValueError("at least one branch cycle is required")
    result = identity(len(branch_cycles[0]))
    for branch_cycle in branch_cycles:
        result = compose(result, branch_cycle)
    return result


def cycle_type(permutation: Permutation) -> tuple[int, ...]:
    seen: set[int] = set()
    lengths: list[int] = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        point = start
        length = 0
        while point not in seen:
            seen.add(point)
            point = permutation[point]
            length += 1
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def cycle_length_at(permutation: Permutation, point: int) -> int:
    image = permutation[point]
    length = 1
    while image != point:
        image = permutation[image]
        length += 1
    return length


def ramification_contribution(permutation: Permutation) -> int:
    return len(permutation) - len(cycle_type(permutation))


def transposition(degree: int, first: int, second: int) -> Permutation:
    result = list(range(degree))
    result[first], result[second] = result[second], result[first]
    return tuple(result)


def all_transpositions(degree: int) -> tuple[Permutation, ...]:
    return tuple(
        transposition(degree, first, second)
        for first, second in combinations(range(degree), 2)
    )


def extend(permutation: Permutation, degree: int, offset: int = 0) -> Permutation:
    """Embed a permutation in consecutive letters starting at ``offset``."""

    if offset + len(permutation) > degree:
        raise ValueError("the requested permutation embedding is too large")
    result = list(range(degree))
    for index, image in enumerate(permutation):
        result[offset + index] = offset + image
    return tuple(result)


def generated_group(generators: Sequence[Permutation]) -> set[Permutation]:
    degree = len(generators[0])
    result = {identity(degree)}
    frontier = [identity(degree)]
    moves = (*generators, *(inverse(generator) for generator in generators))
    while frontier:
        current = frontier.pop()
        for move in moves:
            candidate = compose(move, current)
            if candidate not in result:
                result.add(candidate)
                frontier.append(candidate)
    return result


def generated_orbit(generators: Sequence[Permutation], start: int) -> set[int]:
    orbit = {start}
    frontier = [start]
    while frontier:
        point = frontier.pop()
        for generator in generators:
            image = generator[point]
            if image not in orbit:
                orbit.add(image)
                frontier.append(image)
    return orbit


def canonical_under_simultaneous_conjugacy(
    branch_cycles: BranchTuple,
) -> BranchTuple:
    degree = len(branch_cycles[0])
    return min(
        tuple(conjugate(relabeling, cycle) for cycle in branch_cycles)
        for relabeling in permutations(range(degree))
    )


def nontrivial_cycles(permutation: Permutation) -> list[list[int]]:
    """Return nontrivial cycles in stable one-based notation."""

    seen: set[int] = set()
    cycles: list[list[int]] = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        point = start
        cycle: list[int] = []
        while point not in seen:
            seen.add(point)
            cycle.append(point + 1)
            point = permutation[point]
        if len(cycle) > 1:
            smallest = min(cycle)
            pivot = cycle.index(smallest)
            cycle = cycle[pivot:] + cycle[:pivot]
            cycles.append(cycle)
    return sorted(cycles)


def terminal_a6_triples() -> tuple[BranchTuple, ...]:
    """Reproduce the five terminal triples with the five-cycle fixed."""

    sigma_zero = (1, 2, 3, 4, 0, 5)
    compatible: list[BranchTuple] = []
    for sigma_infinity in permutations(range(6)):
        if cycle_type(sigma_infinity) != (3, 3):
            continue
        sigma_third = inverse(compose(sigma_zero, sigma_infinity))
        if cycle_type(sigma_third) != (3, 1, 1, 1):
            continue
        triple = (sigma_zero, sigma_infinity, sigma_third)
        group = generated_group(triple[:2])
        if len(group) != 360:
            raise AssertionError("a compatible terminal triple lost A6")
        if permutation_product(triple) != identity(6):
            raise AssertionError("the terminal meridian product changed")
        compatible.append(triple)
    if len(compatible) != 5:
        raise AssertionError("the terminal passport no longer has five triples")
    return tuple(compatible)


def riemann_hurwitz_source_genus(
    degree: int,
    branch_cycles: Sequence[Permutation],
) -> int | None:
    """Return the source genus for a cover of P1, or ``None`` if impossible."""

    ramification = sum(
        ramification_contribution(branch_cycle)
        for branch_cycle in branch_cycles
    )
    twice_genus = 2 - 2 * degree + ramification
    if twice_genus < 0 or twice_genus % 2:
        return None
    return twice_genus // 2


def two_simple_cycle_degree_table(
    triples: Sequence[BranchTuple],
) -> list[dict[str, int | None]]:
    """Audit all product-one transposition pairs in degrees six through eight."""

    rows: list[dict[str, int | None]] = []
    for degree in range(6, 9):
        normalized = 0
        transitive = 0
        for triple in triples:
            extended = tuple(extend(cycle, degree) for cycle in triple)
            for first in all_transpositions(degree):
                for second in all_transpositions(degree):
                    branch_cycles = (*extended, first, second)
                    if permutation_product(branch_cycles) != identity(degree):
                        continue
                    normalized += 1
                    if len(generated_orbit(branch_cycles, 0)) == degree:
                        transitive += 1
        genus = riemann_hurwitz_source_genus(
            degree,
            (*tuple(extend(cycle, degree) for cycle in triples[0]),
             *all_transpositions(degree)[:1],
             *all_transpositions(degree)[:1]),
        )
        rows.append(
            {
                "degree": degree,
                "product_one_normalized_tuples": normalized,
                "transitive_normalized_tuples": transitive,
                "riemann_hurwitz_genus": genus,
            }
        )
    expected = [
        {
            "degree": 6,
            "product_one_normalized_tuples": 75,
            "transitive_normalized_tuples": 75,
            "riemann_hurwitz_genus": 1,
        },
        {
            "degree": 7,
            "product_one_normalized_tuples": 105,
            "transitive_normalized_tuples": 30,
            "riemann_hurwitz_genus": 0,
        },
        {
            "degree": 8,
            "product_one_normalized_tuples": 140,
            "transitive_normalized_tuples": 0,
            "riemann_hurwitz_genus": None,
        },
    ]
    if rows != expected:
        raise AssertionError("the two-simple-cycle degree table changed")
    return rows


def degree_seven_candidates(
    triples: Sequence[BranchTuple],
) -> tuple[BranchTuple, ...]:
    """Enumerate the connected degree-seven candidates exhaustively."""

    candidates: list[BranchTuple] = []
    for triple in triples:
        extended = tuple(extend(cycle, 7) for cycle in triple)
        for first in all_transpositions(7):
            for second in all_transpositions(7):
                branch_cycles = (*extended, first, second)
                if permutation_product(branch_cycles) != identity(7):
                    continue
                if len(generated_orbit(branch_cycles, 0)) != 7:
                    continue
                if first != second:
                    raise AssertionError("product one no longer pairs the cycles")
                if first[6] == 6 or first[6] >= 6:
                    raise AssertionError("the new sheet is not used")
                if riemann_hurwitz_source_genus(7, branch_cycles) != 0:
                    raise AssertionError("a degree-seven candidate lost genus zero")
                group = generated_group(branch_cycles)
                if len(group) != 5040:
                    raise AssertionError("a connected candidate no longer generates S7")
                candidates.append(branch_cycles)
    if len(candidates) != 30:
        raise AssertionError("the labelled degree-seven count changed")
    return tuple(candidates)


def unique_external_sheet(branch_cycles: BranchTuple) -> int:
    fixed = [
        point
        for point in range(len(branch_cycles[0]))
        if all(cycle[point] == point for cycle in branch_cycles[:3])
    ]
    if len(fixed) != 1:
        raise AssertionError("the embedded A6 packet has no unique exterior sheet")
    return fixed[0]


def conjugacy_class_records(
    candidates: Sequence[BranchTuple],
) -> list[dict[str, object]]:
    classes = sorted(
        {
            canonical_under_simultaneous_conjugacy(candidate)
            for candidate in candidates
        }
    )
    if len(classes) != 6:
        raise AssertionError("the degree-seven Nielsen class count changed")

    records: list[dict[str, object]] = []
    for branch_cycles in classes:
        external = unique_external_sheet(branch_cycles)
        connector = branch_cycles[3]
        anchor = connector[external]
        if branch_cycles[4] != connector:
            raise AssertionError("the simple connector is no longer paired")
        terminal_cycle_lengths = [
            cycle_length_at(cycle, anchor) for cycle in branch_cycles[:3]
        ]
        records.append(
            {
                "branch_cycles_one_based": [
                    nontrivial_cycles(cycle) for cycle in branch_cycles
                ],
                "external_sheet_one_based": external + 1,
                "core_anchor_one_based": anchor + 1,
                "anchor_cycle_lengths_at_terminal_values": terminal_cycle_lengths,
                "anchor_specializations": {
                    "over_zero": (
                        "interior_attachment_s=-1"
                        if terminal_cycle_lengths[0] == 5
                        else "source_endpoint_s=0"
                    ),
                    "over_infinity": (
                        "one_of_two_interior_denominator_attachments"
                    ),
                    "over_125_over_729": (
                        "source_endpoint_s=infinity"
                        if terminal_cycle_lengths[2] == 3
                        else "one_of_three_simple_interior_points"
                    ),
                },
            }
        )
    profile_counts: dict[str, int] = {}
    for record in records:
        profile = ",".join(
            str(value)
            for value in record["anchor_cycle_lengths_at_terminal_values"]
        )
        profile_counts[profile] = profile_counts.get(profile, 0) + 1
    if profile_counts != {"1,3,3": 1, "5,3,1": 3, "5,3,3": 2}:
        raise AssertionError("the attachment profile census changed")
    return records


def marked_incidence_filter(
    records: Sequence[dict[str, object]],
) -> dict[str, object]:
    """Refine the six classes by the certified endpoint/interior markings."""

    signatures = [
        tuple(record["anchor_cycle_lengths_at_terminal_values"])
        for record in records
    ]
    interior_over_zero = [signature for signature in signatures if signature[0] == 5]
    interior_over_third = [signature for signature in signatures if signature[2] == 1]
    avoids_both_source_endpoints = [
        signature
        for signature in signatures
        if signature[0] == 5 and signature[2] == 1
    ]
    assert len(interior_over_zero) == 5
    assert len(interior_over_third) == 3
    assert len(avoids_both_source_endpoints) == 3
    assert set(avoids_both_source_endpoints) == {(5, 3, 1)}
    return {
        "certified_markings": {
            "over_zero": (
                "the fixed point is the source endpoint s=0; the five-cycle "
                "is the interior attachment s=-1"
            ),
            "over_infinity": (
                "both three-cycles are the two interior denominator-root "
                "attachments"
            ),
            "over_125_over_729": (
                "the three-cycle is the source endpoint s=infinity; the fixed "
                "points are the three simple interior roots"
            ),
        },
        "class_counts": {
            "anchor_in_interior_attachment_over_zero": len(interior_over_zero),
            "anchor_in_an_interior_attachment_over_infinity": len(records),
            "anchor_avoids_source_endpoint_over_third_value": len(
                interior_over_third
            ),
            "anchor_avoids_both_source_endpoints": len(
                avoids_both_source_endpoints
            ),
        },
        "strongest_naive_marked_filter_survivors": {
            "signature": [5, 3, 1],
            "simultaneous_conjugacy_classes": 3,
            "verdict": "survives",
        },
        "scope_warning": (
            "a connector sheet specializing to a marked residue point is not "
            "a proof that the spectator boundary component glues there; that "
            "requires a toroidal source-node compatibility theorem"
        ),
    }


def paired_star_witness(
    triple: BranchTuple,
    remaining_degree: int,
) -> BranchTuple:
    """One A6 block plus paired connectors to every new sheet."""

    if remaining_degree < 1:
        raise ValueError("the paired-star witness needs a spectator sheet")
    degree = 6 + remaining_degree
    branch_cycles: list[Permutation] = [
        extend(cycle, degree) for cycle in triple
    ]
    for new_sheet in range(6, degree):
        connector = transposition(degree, 0, new_sheet)
        branch_cycles.extend((connector, connector))
    result = tuple(branch_cycles)
    if permutation_product(result) != identity(degree):
        raise AssertionError("the paired-star meridian product is not one")
    if len(generated_orbit(result, 0)) != degree:
        raise AssertionError("the paired-star witness is not transitive")
    if riemann_hurwitz_source_genus(degree, result) != 0:
        raise AssertionError("the paired-star witness is not genus zero")
    return result


def full_symmetric_witness_criterion(
    a6_group: Iterable[Permutation],
    remaining_degree: int,
) -> None:
    """Check the elementary conjugation proof that the star generates S_n."""

    degree = 6 + remaining_degree
    core_group = tuple(a6_group)
    if {element[0] for element in core_group} != set(range(6)):
        raise AssertionError("the core group is not transitive")
    for new_sheet in range(6, degree):
        connector = transposition(degree, 0, new_sheet)
        conjugate_edges = set()
        for element in core_group:
            extended = extend(element, degree)
            conjugated = conjugate(extended, connector)
            moved = tuple(
                point for point in range(degree) if conjugated[point] != point
            )
            conjugate_edges.add(moved)
        expected = {(core_sheet, new_sheet) for core_sheet in range(6)}
        normalized = {tuple(sorted(edge)) for edge in conjugate_edges}
        if normalized != expected:
            raise AssertionError("the connector conjugates missed a core sheet")
    # The derived transpositions join every vertex to the core.  A connected
    # graph's edge transpositions generate its full symmetric group.


def double_packet_witness(triple: BranchTuple) -> BranchTuple:
    """Two disjoint A6 packets joined by a repeated cross transposition."""

    degree = 12
    first = tuple(extend(cycle, degree, 0) for cycle in triple)
    second = tuple(extend(cycle, degree, 6) for cycle in triple)
    connector = transposition(degree, 0, 6)
    branch_cycles = (*first, *second, connector, connector)
    if permutation_product(branch_cycles) != identity(degree):
        raise AssertionError("the double-packet meridian product is not one")
    if len(generated_orbit(branch_cycles, 0)) != degree:
        raise AssertionError("the double-packet witness is not transitive")
    if riemann_hurwitz_source_genus(degree, branch_cycles) != 0:
        raise AssertionError("the double-packet witness is not genus zero")
    return branch_cycles


def double_packet_full_symmetric_criterion(
    a6_group: Iterable[Permutation],
) -> None:
    """Check that conjugated cross edges generate the full S12."""

    core_group = tuple(a6_group)
    representatives = {
        image: next(element for element in core_group if element[0] == image)
        for image in range(6)
    }
    connector = transposition(12, 0, 6)
    derived_edges: set[tuple[int, int]] = set()
    for left_image in range(6):
        left = extend(representatives[left_image], 12, 0)
        for right_image in range(6):
            right = extend(representatives[right_image], 12, 6)
            relabeling = compose(left, right)
            conjugated = conjugate(relabeling, connector)
            edge = tuple(
                point for point in range(12) if conjugated[point] != point
            )
            if len(edge) != 2:
                raise AssertionError("a conjugated connector is not an edge")
            derived_edges.add(edge)
    expected = {
        (left, right) for left in range(6) for right in range(6, 12)
    }
    if derived_edges != expected:
        raise AssertionError("the cross-packet connector missed a derived edge")
    # The complete bipartite edge graph is connected, so its edge
    # transpositions generate S12.


def build_payload() -> dict[str, object]:
    triples = terminal_a6_triples()
    degree_table = two_simple_cycle_degree_table(triples)
    candidates = degree_seven_candidates(triples)
    class_records = conjugacy_class_records(candidates)
    marked_filter = marked_incidence_filter(class_records)

    core_group = generated_group(triples[0][:2])
    remaining_degree_rows: list[dict[str, int | bool]] = []
    for remaining_degree in range(1, 20):
        witness = paired_star_witness(triples[0], remaining_degree)
        full_symmetric_witness_criterion(core_group, remaining_degree)
        degree = 6 + remaining_degree
        remaining_degree_rows.append(
            {
                "remaining_degree": remaining_degree,
                "global_degree": degree,
                "simple_cycle_count": 2 * remaining_degree,
                "total_ramification": sum(
                    ramification_contribution(cycle) for cycle in witness
                ),
                "product_one": True,
                "transitive": True,
                "source_genus": 0,
            }
        )

    double_witness = double_packet_witness(triples[0])
    double_packet_full_symmetric_criterion(core_group)
    return {
        "schema": "plane-jc.f2-a6-simple-spectator-gluing.v2",
        "status": "conditional-finite-enumeration-not-f2-spectator-classification",
        "certified_terminal_row_input_and_permutation_replay": {
            "geometry_source": "plane-jc/F2_TERMINAL_RESIDUE_COVER.md",
            "source_ray": [12, -17],
            "target_ray": [5, 2],
            "transverse_index": 1,
            "residue_degree": 6,
            "terminal_degree": 6,
            "branch_values": ["0", "infinity", "125/729"],
            "passport": [[5, 1], [3, 3], [3, 1, 1, 1]],
            "ramification_contributions": [4, 4, 2],
            "total_ramification": 10,
            "residue_different_packet": [4, 2, 2, 2],
            "transverse_different": 0,
            "target_center": "boundary_at_infinity",
            "compatible_triples_with_fixed_five_cycle": len(triples),
            "generated_group": "A6",
            "meridian_product": "one",
        },
        "conditional_bridge_assumptions": [
            "the terminal A6 triple embeds on six global sheets and fixes all remaining sheets",
            "each of the two squarefree-R simple Kummer orbits supplies exactly one separate simple branch value with transposition monodromy",
            "the A6 triple is followed by those two branch cycles in one compact target-P1 meridian system",
            "there are no other branch cycles on that target curve",
            "a spectator means at least one sheet outside the terminal six-sheet packet",
        ],
        "two_simple_cycle_degree_audit": degree_table,
        "squarefree_conditional_result": {
            "forced_global_degree": 7,
            "remaining_degree": 1,
            "source_genus": 0,
            "normalized_factorizations_fixed_core_and_five_cycle": len(
                candidates
            ),
            "simultaneous_conjugacy_classes_with_branch_labels_fixed": len(
                class_records
            ),
            "fully_sheet_labeled_factorizations": 5040 * len(class_records),
            "generated_group": "S7",
            "class_records": class_records,
            "verdict": "survives_as_six_exact_permutation_gluings",
        },
        "marked_terminal_incidence_filter": marked_filter,
        "strictly_disjoint_spectator_verdict": (
            "a transposition supported away from the A6 block cannot make the "
            "global action transitive"
        ),
        "all_remaining_degree_witness": {
            "formula": (
                "append (1,j),(1,j) for each new sheet j; then "
                "N=6+k, branch_count=3+2k, ramification=10+2k=2N-2"
            ),
            "generated_group": "S_(6+k)",
            "proof_criterion": (
                "A6 is transitive on its six letters; conjugating each star "
                "connector gives all core-to-new transpositions, whose "
                "connected edge graph generates the full symmetric group"
            ),
            "checked_remaining_degrees": remaining_degree_rows,
            "degree_25_diagnostic": {
                "warning": (
                    "25 is the common-edge polynomial degree, not a certified "
                    "global geometric degree"
                ),
                "if_global_degree_were_25": 25,
                "required_simple_cycles": 38,
                "deficit_from_one_cycle_per_two_R_orbits": 36,
                "verdict": "additional cycles would be required, not a contradiction",
            },
        },
        "double_packet_diagnostic": {
            "degree": 12,
            "core_ramification": 20,
            "simple_cross_connectors": 2,
            "total_ramification": sum(
                ramification_contribution(cycle) for cycle in double_witness
            ),
            "product_one": True,
            "transitive": True,
            "source_genus": 0,
            "generated_group": "S12 by the same conjugated-edge argument",
            "scope": (
                "an abstract witness if another global row supplies the two "
                "connectors; the double-root F2 stratum itself supplies no "
                "certified simple spectator cycles"
            ),
        },
        "claim_boundary": {
            "proved": (
                "the finite conditional enumeration and the abstract "
                "paired-connector permutation witnesses"
            ),
            "not_proved": [
                "that a simple R Kummer orbit is a branch value or a transposition",
                "that the spectator rows lie over the terminal target component",
                "that the relevant global cover curve is this compact P1 cover",
                "that no purity-forced affine or other boundary branch cycles occur",
                "that common-edge degree 25 equals geometric field degree",
                "exclusion of (75,125) or of the F2 family",
            ],
        },
        "reproduction_command": (
            ".venv/bin/python "
            "plane-jc/cas/verify_f2_a6_simple_spectator_gluing.py"
        ),
        "software": {"python": "standard library"},
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--artifact", type=Path, default=ARTIFACT)
    args = parser.parse_args()

    payload = build_payload()
    artifact = args.artifact.resolve()
    if args.refresh:
        artifact.parent.mkdir(parents=True, exist_ok=True)
        artifact.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        try:
            display_path = artifact.relative_to(ROOT)
        except ValueError:
            display_path = artifact
        print(f"WROTE {display_path}")
    else:
        expected = json.loads(artifact.read_text())
        if expected != payload:
            raise AssertionError(
                "the pinned F2 A6 simple-spectator artifact is stale; "
                "inspect before --refresh"
            )

    print("F2_A6_TERMINAL_TRIPLES=5")
    print("F2_A6_TWO_SIMPLE_DEGREE7_NORMALIZED=30")
    print("F2_A6_TWO_SIMPLE_DEGREE7_CONJUGACY_CLASSES=6")
    print("F2_A6_MARKED_INTERIOR_FILTER_CONJUGACY_CLASSES=3")
    print("F2_A6_TWO_SIMPLE_DEGREE7_MONODROMY=S7")
    print("F2_A6_TWO_SIMPLE_DEGREE7_RIEMANN_HURWITZ_GENUS=0")
    print("F2_A6_PAIRED_STAR_WITNESSES_REMAINING_DEGREE=1..19")
    print("F2_A6_SIMPLE_SPECTATOR_GLUING_CONDITIONAL_PASS")


if __name__ == "__main__":
    main()
