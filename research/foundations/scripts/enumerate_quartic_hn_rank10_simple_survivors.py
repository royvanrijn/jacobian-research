#!/usr/bin/env python3
"""Enumerate simple rank-four OP-QHNW10 abstract survivors on ten points.

The input is the complete rank-four, nine-element catalogue from
``matroid-database==0.3``.  Every simple rank-four matroid on ten elements has
a deletion to a simple rank-four matroid on nine elements, so the search can
be performed by enumerating the simple single-element extensions of those
catalogue entries.

Rank-preserving single-element extensions are encoded exactly by modular cuts
of the flat lattice.  In rank four, simplicity leaves Boolean choices only on
rank-two and rank-three flats.  Z3 enumerates the residual modular cuts after
the cyclic-complement conditions have forced the necessary hyperplanes.  No
finite-field representability assumption is used.  Every resulting extension
is independently checked and the labelled outputs are quotiented by abstract
matroid isomorphism.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import lzma
import sys
from collections import Counter
from itertools import combinations
from pathlib import Path

import z3

from verify_quartic_hn_rank10_matroid_survivors import (
    assert_basis_exchange,
    bases_from_rankline,
    hyperplanes,
    rank_table,
    revlex_four_sets,
)


CATALOGUE_LINE_COUNT = 190_214
CATALOGUE_SIMPLE_COUNT = 185_981
SOURCE_WHEEL_SHA256 = "85d0304575784ceb4797014fb5a761443d2a0bbf01aafc38ef293d2c64a1b5ce"


def is_simple(n: int, bases: list[int]) -> bool:
    support = 0
    pairs: set[int] = set()
    for basis in bases:
        support |= basis
        elements = [element for element in range(n) if basis & (1 << element)]
        pairs.update((1 << left) | (1 << right) for left, right in combinations(elements, 2))
    return support == (1 << n) - 1 and len(pairs) == n * (n - 1) // 2


def closures(n: int, ranks: list[int]) -> list[int]:
    result: list[int] = []
    full = (1 << n) - 1
    for mask in range(1 << n):
        closure = mask
        for element in range(n):
            bit = 1 << element
            if not mask & bit and ranks[mask | bit] == ranks[mask]:
                closure |= bit
        result.append(closure & full)
    return result


def is_cyclic(mask: int, ranks: list[int]) -> bool:
    rank = ranks[mask]
    return all(
        ranks[mask ^ (1 << element)] == rank
        for element in range(len(ranks).bit_length() - 1)
        if mask & (1 << element)
    )


def flat_data(n: int, ranks: list[int], closure: list[int]) -> tuple[list[int], list[int]]:
    flats = [mask for mask in range(1 << n) if closure[mask] == mask]
    return (
        [flat for flat in flats if ranks[flat] == 2],
        [flat for flat in flats if ranks[flat] == 3],
    )


def forced_planes_from_new_bases(
    n: int, ranks: list[int], closure: list[int]
) -> set[int]:
    """Planes that must enter the cut to block a bad basis I+e."""
    full = (1 << n) - 1
    forced: set[int] = set()
    for triple in combinations(range(n), 3):
        mask = sum(1 << element for element in triple)
        if ranks[mask] == 3 and is_cyclic(full ^ mask, ranks):
            plane = closure[mask]
            if ranks[plane] != 3:
                raise AssertionError((triple, plane, ranks[plane]))
            forced.add(plane)
    return forced


def close_forced_cut(
    ranks: list[int],
    lines: list[int],
    planes: list[int],
    forced_planes: set[int],
) -> tuple[set[int], set[int], str | None]:
    """Close forced planes under modular pairs, or return an exact conflict."""
    selected_lines: set[int] = set()
    selected_planes = set(forced_planes)
    changed = True
    while changed:
        changed = False
        for left, right in combinations(sorted(selected_planes), 2):
            intersection = left & right
            if ranks[intersection] == 2:
                line = intersection
                if line not in selected_lines:
                    selected_lines.add(line)
                    changed = True
        for line in tuple(selected_lines):
            for plane in planes:
                if line & ~plane == 0 and plane not in selected_planes:
                    selected_planes.add(plane)
                    changed = True

    for line in selected_lines:
        for plane in selected_planes:
            if line & ~plane and ranks[line & plane] == 1:
                return selected_lines, selected_planes, "line_plane_conflict"
    for left, right in combinations(sorted(selected_lines), 2):
        if ranks[left & right] + ranks[left | right] == 4:
            return selected_lines, selected_planes, "line_line_conflict"
    if any(plane.bit_count() >= 7 for plane in selected_planes):
        return selected_lines, selected_planes, "selected_large_plane"
    return selected_lines, selected_planes, None


def cut_contains(
    flat: int,
    ranks: list[int],
    selected_lines: set[int],
    selected_planes: set[int],
) -> bool:
    rank = ranks[flat]
    if rank <= 1:
        return False
    if rank == 2:
        return flat in selected_lines
    if rank == 3:
        return flat in selected_planes
    if rank == 4:
        return True
    raise AssertionError((flat, rank))


def extension_rank_without_new(
    mask: int,
    ranks: list[int],
    closure: list[int],
    selected_lines: set[int],
    selected_planes: set[int],
) -> int:
    rank = ranks[mask]
    return rank if cut_contains(closure[mask], ranks, selected_lines, selected_planes) else rank + 1


def old_basis_complements_pass(
    n: int,
    bases: list[int],
    ranks: list[int],
    closure: list[int],
    selected_lines: set[int],
    selected_planes: set[int],
) -> bool:
    full = (1 << n) - 1
    for basis in bases:
        old_complement = full ^ basis
        total_rank = extension_rank_without_new(
            old_complement, ranks, closure, selected_lines, selected_planes
        )
        new_is_coloop = ranks[old_complement] < total_rank
        old_coloop = any(
            extension_rank_without_new(
                old_complement ^ (1 << element),
                ranks,
                closure,
                selected_lines,
                selected_planes,
            )
            < total_rank
            for element in range(n)
            if old_complement & (1 << element)
        )
        if not (new_is_coloop or old_coloop):
            return False
    return True


def z3_cut_models(
    n: int,
    bases: list[int],
    ranks: list[int],
    closure: list[int],
    lines: list[int],
    planes: list[int],
    forced_planes: set[int],
) -> list[tuple[set[int], set[int]]]:
    line_vars = {line: z3.Bool(f"l_{line}") for line in lines}
    plane_vars = {plane: z3.Bool(f"p_{plane}") for plane in planes}
    solver = z3.Solver()

    def flat_var(flat: int) -> z3.BoolRef:
        rank = ranks[flat]
        if rank <= 1:
            return z3.BoolVal(False)
        if rank == 2:
            return line_vars[flat]
        if rank == 3:
            return plane_vars[flat]
        if rank == 4:
            return z3.BoolVal(True)
        raise AssertionError((flat, rank))

    for plane in forced_planes:
        solver.add(plane_vars[plane])
    for plane in planes:
        if plane.bit_count() >= 7:
            solver.add(z3.Not(plane_vars[plane]))
    for line in lines:
        for plane in planes:
            if line & ~plane == 0:
                solver.add(z3.Implies(line_vars[line], plane_vars[plane]))
            elif ranks[line & plane] == 1:
                solver.add(z3.Not(z3.And(line_vars[line], plane_vars[plane])))
    for left, right in combinations(lines, 2):
        if ranks[left & right] + ranks[left | right] == 4:
            solver.add(z3.Not(z3.And(line_vars[left], line_vars[right])))
    for left, right in combinations(planes, 2):
        intersection = left & right
        if ranks[intersection] == 2:
            solver.add(z3.Implies(z3.And(plane_vars[left], plane_vars[right]), line_vars[intersection]))

    def rank_with_new(old_mask: int) -> z3.ArithRef:
        return z3.IntVal(ranks[old_mask]) + z3.If(flat_var(closure[old_mask]), 0, 1)

    full = (1 << n) - 1
    for basis in bases:
        old_complement = full ^ basis
        total_rank = rank_with_new(old_complement)
        coloop_conditions: list[z3.BoolRef] = [
            z3.IntVal(ranks[old_complement]) < total_rank
        ]
        coloop_conditions.extend(
            rank_with_new(old_complement ^ (1 << element)) < total_rank
            for element in range(n)
            if old_complement & (1 << element)
        )
        solver.add(z3.Or(*coloop_conditions))

    variables = [line_vars[line] for line in lines] + [plane_vars[plane] for plane in planes]
    result: list[tuple[set[int], set[int]]] = []
    while solver.check() == z3.sat:
        model = solver.model()
        selected_lines = {
            line for line, variable in line_vars.items() if z3.is_true(model.eval(variable, model_completion=True))
        }
        selected_planes = {
            plane for plane, variable in plane_vars.items() if z3.is_true(model.eval(variable, model_completion=True))
        }
        result.append((selected_lines, selected_planes))
        solver.add(
            z3.Or(
                *(variable != model.eval(variable, model_completion=True) for variable in variables)
            )
        )
    return result


def extension_bases(
    n: int,
    old_bases: list[int],
    ranks: list[int],
    closure: list[int],
    selected_lines: set[int],
    selected_planes: set[int],
) -> list[int]:
    result = list(old_bases)
    new_bit = 1 << n
    for triple in combinations(range(n), 3):
        mask = sum(1 << element for element in triple)
        if ranks[mask] == 3 and not cut_contains(
            closure[mask], ranks, selected_lines, selected_planes
        ):
            result.append(mask | new_bit)
    return sorted(result)


def rankline_from_bases(n: int, bases: list[int]) -> str:
    basis_set = set(bases)
    return "".join(
        "*" if sum(1 << element for element in subset) in basis_set else "0"
        for subset in revlex_four_sets(n)
    )


def assert_simple_survivor(n: int, bases: list[int]) -> None:
    if not is_simple(n, bases):
        raise AssertionError("extension is not simple")
    assert_basis_exchange(n, bases)
    ranks = rank_table(n, bases)
    for plane in hyperplanes(n, ranks):
        if n - plane.bit_count() < 3:
            raise AssertionError(("small cocircuit", plane))
    full = (1 << n) - 1
    for basis in bases:
        complement = full ^ basis
        rank = ranks[complement]
        if all(
            ranks[complement ^ (1 << element)] == rank
            for element in range(n)
            if complement & (1 << element)
        ):
            raise AssertionError(("cyclic complement", basis, complement))


def element_invariants(n: int, bases: list[int], ranks: list[int]) -> list[tuple[object, ...]]:
    planes = hyperplanes(n, ranks)
    basis_degrees = [sum(bool(basis & (1 << element)) for basis in bases) for element in range(n)]
    pair_degrees = {
        (left, right): sum(bool(basis & (1 << left)) and bool(basis & (1 << right)) for basis in bases)
        for left, right in combinations(range(n), 2)
    }
    return [
        (
            basis_degrees[element],
            tuple(sorted(pair_degrees[tuple(sorted((element, other)))] for other in range(n) if other != element)),
            tuple(sorted(plane.bit_count() for plane in planes if plane & (1 << element))),
        )
        for element in range(n)
    ]


def partial_incidence_profile(family: list[int], assigned: list[int]) -> Counter[int]:
    profile: Counter[int] = Counter()
    for member in family:
        pattern = sum(1 << index for index, element in enumerate(assigned) if member & (1 << element))
        profile[pattern] += 1
    return profile


def partial_plane_profile(planes: list[int], assigned: list[int]) -> Counter[tuple[int, int]]:
    profile: Counter[tuple[int, int]] = Counter()
    for plane in planes:
        pattern = sum(1 << index for index, element in enumerate(assigned) if plane & (1 << element))
        profile[(plane.bit_count(), pattern)] += 1
    return profile


def are_isomorphic(n: int, left_bases: list[int], right_bases: list[int]) -> bool:
    if len(left_bases) != len(right_bases):
        return False
    left_ranks, right_ranks = rank_table(n, left_bases), rank_table(n, right_bases)
    left_planes, right_planes = hyperplanes(n, left_ranks), hyperplanes(n, right_ranks)
    if Counter(plane.bit_count() for plane in left_planes) != Counter(plane.bit_count() for plane in right_planes):
        return False
    left_invariants = element_invariants(n, left_bases, left_ranks)
    right_invariants = element_invariants(n, right_bases, right_ranks)
    if Counter(left_invariants) != Counter(right_invariants):
        return False

    source_order = sorted(
        range(n),
        key=lambda element: (Counter(left_invariants)[left_invariants[element]], left_invariants[element], element),
    )
    assigned_left: list[int] = []
    assigned_right: list[int] = []
    used_right: set[int] = set()

    def search(depth: int) -> bool:
        if depth == n:
            return True
        source = source_order[depth]
        assigned_left.append(source)
        left_basis_profile = partial_incidence_profile(left_bases, assigned_left)
        left_plane_profile = partial_plane_profile(left_planes, assigned_left)
        for target in range(n):
            if target in used_right or right_invariants[target] != left_invariants[source]:
                continue
            assigned_right.append(target)
            if (
                partial_incidence_profile(right_bases, assigned_right) == left_basis_profile
                and partial_plane_profile(right_planes, assigned_right) == left_plane_profile
            ):
                used_right.add(target)
                if search(depth + 1):
                    return True
                used_right.remove(target)
            assigned_right.pop()
        assigned_left.pop()
        return False

    return search(0)


def enumerate_survivors(database_root: Path) -> dict[str, object]:
    path = database_root / "allr4n09.txt.xz"
    if not path.is_file():
        raise FileNotFoundError(path)

    scan_counts: Counter[str] = Counter()
    residual: list[dict[str, object]] = []
    line_count = 0
    simple_count = 0
    with lzma.open(path, "rt") as handle:
        for local_index, line in enumerate(handle):
            line_count += 1
            rankline = line.strip()
            bases = bases_from_rankline(9, rankline)
            if not is_simple(9, bases):
                continue
            simple_count += 1
            ranks = rank_table(9, bases)
            closure = closures(9, ranks)
            lines, planes = flat_data(9, ranks, closure)
            if any(plane.bit_count() >= 8 for plane in planes):
                scan_counts["coloop_deletion"] += 1
                continue
            forced_planes = forced_planes_from_new_bases(9, ranks, closure)
            selected_lines, selected_planes, conflict = close_forced_cut(
                ranks, lines, planes, forced_planes
            )
            if conflict is not None:
                scan_counts[conflict] += 1
                continue
            minimal_passes = old_basis_complements_pass(
                9, bases, ranks, closure, selected_lines, selected_planes
            )
            scan_counts["minimal_survivor" if minimal_passes else "optional_cut_needed"] += 1
            residual.append(
                {
                    "catalogue_index": local_index,
                    "rankline": rankline,
                    "minimal_cut_passes": minimal_passes,
                }
            )
            if simple_count % 20_000 == 0:
                print(
                    f"QHNW10_SIMPLE_SCAN lines={line_count} simple={simple_count} residual={len(residual)}",
                    file=sys.stderr,
                    flush=True,
                )

    if line_count != CATALOGUE_LINE_COUNT or simple_count != CATALOGUE_SIMPLE_COUNT:
        raise AssertionError((line_count, simple_count))

    representatives: list[dict[str, object]] = []
    labelled_model_count = 0
    model_counts_by_deletion: dict[int, int] = {}
    for residual_record in residual:
        local_index = int(residual_record["catalogue_index"])
        rankline = str(residual_record["rankline"])
        bases = bases_from_rankline(9, rankline)
        ranks = rank_table(9, bases)
        closure = closures(9, ranks)
        lines, planes = flat_data(9, ranks, closure)
        forced_planes = forced_planes_from_new_bases(9, ranks, closure)
        models = z3_cut_models(9, bases, ranks, closure, lines, planes, forced_planes)
        model_counts_by_deletion[local_index] = len(models)
        print(
            f"QHNW10_SIMPLE_RESIDUAL index={local_index} models={len(models)}",
            file=sys.stderr,
            flush=True,
        )
        for selected_lines, selected_planes in models:
            candidate = extension_bases(
                9, bases, ranks, closure, selected_lines, selected_planes
            )
            assert_simple_survivor(10, candidate)
            labelled_model_count += 1
            provenance = {
                "deletion_catalogue_index": local_index,
                "selected_rank2_flats": sorted(selected_lines),
                "selected_rank3_flats": sorted(selected_planes),
            }
            match = next(
                (
                    representative
                    for representative in representatives
                    if are_isomorphic(10, candidate, list(representative["bases"]))
                ),
                None,
            )
            if match is None:
                representatives.append(
                    {
                        "bases": candidate,
                        "rankline": rankline_from_bases(10, candidate),
                        "labelled_extension_models": 1,
                        "deletion_catalogue_indices": [local_index],
                        "first_provenance": provenance,
                    }
                )
            else:
                match["labelled_extension_models"] = int(match["labelled_extension_models"]) + 1
                source_indices = list(match["deletion_catalogue_indices"])
                if local_index not in source_indices:
                    source_indices.append(local_index)
                    match["deletion_catalogue_indices"] = sorted(source_indices)

    records: list[dict[str, object]] = []
    for survivor_index, representative in enumerate(representatives):
        records.append(
            {
                "survivor_index": survivor_index,
                "rankline": representative["rankline"],
                "bases": representative["bases"],
                "basis_count": len(representative["bases"]),
                "labelled_extension_models": representative["labelled_extension_models"],
                "deletion_catalogue_indices": representative["deletion_catalogue_indices"],
                "first_provenance": representative["first_provenance"],
            }
        )

    return {
        "schema": "quartic-hn-rank10-simple-survivors-v1",
        "source": {
            "package": "matroid-database==0.3",
            "source_wheel_sha256": SOURCE_WHEEL_SHA256,
            "catalogue_file": path.name,
            "catalogue_file_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "catalogue_line_count": line_count,
            "simple_deletion_count": simple_count,
        },
        "scope": {
            "rank": 4,
            "total_elements": 10,
            "simple": True,
            "characteristic_zero_representability_imposed": False,
            "finite_field_representability_imposed": False,
        },
        "counts": {
            "scan_outcomes": dict(sorted(scan_counts.items())),
            "residual_deletion_types": len(residual),
            "labelled_extension_models": labelled_model_count,
            "abstract_isomorphism_types": len(records),
            "models_by_deletion_catalogue_index": {
                str(index): count for index, count in sorted(model_counts_by_deletion.items())
            },
        },
        "residual_deletions": residual,
        "records": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--database-root",
        type=Path,
        required=True,
        help="path containing allr4n09.txt.xz from matroid-database==0.3",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("artifacts/generated-results/quartic_hn_rank10_simple_survivors.json"),
    )
    arguments = parser.parse_args()
    result = enumerate_survivors(arguments.database_root.resolve())
    output = arguments.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(f"QHNW10_SIMPLE_RESIDUAL_DELETIONS={result['counts']['residual_deletion_types']}")
    print(f"QHNW10_SIMPLE_LABELLED_EXTENSIONS={result['counts']['labelled_extension_models']}")
    print(f"QHNW10_SIMPLE_ABSTRACT_SURVIVORS={result['counts']['abstract_isomorphism_types']}")
    print(f"QHNW10_SIMPLE_CENSUS_ARTIFACT={output}")


if __name__ == "__main__":
    main()
