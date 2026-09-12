#!/usr/bin/env python3
"""Enumerate the complete OP-QHNW10 Gale-loop survivor slice.

This is a catalogue extraction script, not the checked-in theorem replay.
It reads the rank-four revlex files from ``matroid-database==0.3``, filters
simple nonzero-column simplifications, assigns parallel multiplicities and
one or more Gale loops to reach ten elements, and quotients the surviving
colourings by the automorphism group of the simplification.

The output is an exact frozen JSON census.  Characteristic-zero realization
ideals and Hessian-trace obstructions are handled by a separate verifier.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import lzma
import sys
from collections import Counter, defaultdict
from itertools import combinations, permutations, product
from pathlib import Path

from verify_quartic_hn_rank10_matroid_survivors import (
    bases_from_rankline,
    hyperplanes,
    rank_table,
)


CATALOGUE_COUNTS = {4: 1, 5: 5, 6: 23, 7: 108, 8: 940, 9: 190_214}
SOURCE_WHEEL_SHA256 = "85d0304575784ceb4797014fb5a761443d2a0bbf01aafc38ef293d2c64a1b5ce"


def catalogue_sha256(paths: list[Path]) -> str:
    digest = hashlib.sha256()
    for path in paths:
        digest.update(path.name.encode("ascii"))
        digest.update(path.read_bytes())
    return digest.hexdigest()


def is_simple(n: int, bases: list[int]) -> bool:
    element_support = 0
    independent_pairs: set[tuple[int, int]] = set()
    for basis in bases:
        element_support |= basis
        elements = [element for element in range(n) if basis & (1 << element)]
        independent_pairs.update(combinations(elements, 2))
    return (
        element_support == (1 << n) - 1
        and len(independent_pairs) == n * (n - 1) // 2
    )


def weight_options(n: int, loops: int) -> tuple[tuple[int, ...], ...]:
    nonloop_total = 10 - loops
    return tuple(
        weights
        for weights in product((1, 2, 3), repeat=n)
        if sum(weights) == nonloop_total
        and sum(weight == 3 for weight in weights) <= 1
    )


def satisfies_weighted_constraints(
    n: int,
    bases: list[int],
    ranks: list[int],
    planes: list[int],
    weights: tuple[int, ...],
) -> bool:
    for plane in planes:
        if sum(
            weights[element]
            for element in range(n)
            if not plane & (1 << element)
        ) < 3:
            return False

    full = (1 << n) - 1
    for basis in bases:
        complement_support = full
        for element in range(n):
            if basis & (1 << element) and weights[element] == 1:
                complement_support ^= 1 << element
        complement_rank = ranks[complement_support]
        has_coloop = any(
            weights[element] - bool(basis & (1 << element)) == 1
            and ranks[complement_support ^ (1 << element)] < complement_rank
            for element in range(n)
        )
        if not has_coloop:
            return False
    return True


def automorphisms(n: int, bases: list[int]) -> tuple[tuple[int, ...], ...]:
    basis_set = set(bases)
    basis_degrees = [
        sum(bool(basis & (1 << element)) for basis in bases)
        for element in range(n)
    ]
    pair_degrees = {
        (left, right): sum(
            bool(basis & (1 << left)) and bool(basis & (1 << right))
            for basis in bases
        )
        for left, right in combinations(range(n), 2)
    }
    invariant_cells: dict[tuple[int, tuple[int, ...]], list[int]] = defaultdict(list)
    for element in range(n):
        pair_profile = tuple(
            sorted(
                pair_degrees[tuple(sorted((element, other)))]
                for other in range(n)
                if other != element
            )
        )
        invariant_cells[(basis_degrees[element], pair_profile)].append(element)
    cells = tuple(tuple(cell) for cell in invariant_cells.values())

    result: list[tuple[int, ...]] = []
    for cell_images in product(*(permutations(cell) for cell in cells)):
        image = list(range(n))
        for cell, targets in zip(cells, cell_images, strict=True):
            for source, target in zip(cell, targets, strict=True):
                image[source] = target
        mapped_bases = {
            sum(1 << image[element] for element in range(n) if basis & (1 << element))
            for basis in bases
        }
        if mapped_bases == basis_set:
            result.append(tuple(image))
    if not result:
        raise AssertionError("identity automorphism missing")
    return tuple(result)


def canonical_weight(
    weights: tuple[int, ...], group: tuple[tuple[int, ...], ...]
) -> tuple[int, ...]:
    images = []
    for automorphism in group:
        image = [0] * len(weights)
        for source, target in enumerate(automorphism):
            image[target] = weights[source]
        images.append(tuple(image))
    return min(images)


def enumerate_survivors(database_root: Path) -> dict[str, object]:
    catalogue_paths = [
        database_root / f"allr4n{n:02d}.txt.xz" for n in range(4, 10)
    ]
    for path in catalogue_paths:
        if not path.is_file():
            raise FileNotFoundError(path)

    coloured: dict[
        tuple[int, int, int, str], set[tuple[int, ...]]
    ] = defaultdict(set)
    labelled_counts: Counter[int] = Counter()
    simple_counts: Counter[int] = Counter()
    line_counts: Counter[int] = Counter()

    options = {
        (n, loops): weight_options(n, loops)
        for n in range(4, 10)
        for loops in range(1, 7)
    }

    for n, path in zip(range(4, 10), catalogue_paths, strict=True):
        with lzma.open(path, "rt") as handle:
            for local_index, line in enumerate(handle):
                line_counts[n] += 1
                rankline = line.strip()
                bases = bases_from_rankline(n, rankline)
                if not is_simple(n, bases):
                    continue
                simple_counts[n] += 1
                ranks = rank_table(n, bases)
                planes = hyperplanes(n, ranks)
                for loops in range(1, 7):
                    for weights in options[n, loops]:
                        if satisfies_weighted_constraints(
                            n, bases, ranks, planes, weights
                        ):
                            coloured[(loops, n, local_index, rankline)].add(weights)
                            labelled_counts[loops] += 1
        print(
            f"QHNW10_LOOP_CATALOGUE_SCAN n={n} "
            f"lines={line_counts[n]} simple={simple_counts[n]}",
            file=sys.stderr,
            flush=True,
        )

    if dict(line_counts) != CATALOGUE_COUNTS:
        raise AssertionError((line_counts, CATALOGUE_COUNTS))

    automorphism_cache: dict[tuple[int, int], tuple[tuple[int, ...], ...]] = {}
    records: list[dict[str, object]] = []
    coloured_counts: Counter[int] = Counter()
    for (loops, n, local_index, rankline), weights_set in sorted(coloured.items()):
        cache_key = (n, local_index)
        if cache_key not in automorphism_cache:
            bases = bases_from_rankline(n, rankline)
            automorphism_cache[cache_key] = automorphisms(n, bases)
        group = automorphism_cache[cache_key]
        representatives = sorted(
            {canonical_weight(weights, group) for weights in weights_set}
        )
        coloured_counts[loops] += len(representatives)
        records.append(
            {
                "loops": loops,
                "simple_points": n,
                "catalogue_index": local_index,
                "rankline": rankline,
                "automorphism_group_order": len(group),
                "weight_representatives": [
                    "".join(str(weight) for weight in weights)
                    for weights in representatives
                ],
            }
        )

    return {
        "schema": "quartic-hn-rank10-loop-survivors-v1",
        "source": {
            "package": "matroid-database==0.3",
            "source_wheel_sha256": SOURCE_WHEEL_SHA256,
            "catalogue_files_sha256": catalogue_sha256(catalogue_paths),
            "catalogue_line_counts": {
                str(n): line_counts[n] for n in range(4, 10)
            },
        },
        "scope": {
            "rank": 4,
            "total_elements": 10,
            "minimum_loops": 1,
            "characteristic_zero_representability_imposed": False,
            "finite_field_representability_imposed": False,
        },
        "counts": {
            "labelled_survivors_by_loops": {
                str(loops): labelled_counts[loops]
                for loops in sorted(labelled_counts)
            },
            "coloured_isomorphism_types_by_loops": {
                str(loops): coloured_counts[loops]
                for loops in sorted(coloured_counts)
            },
            "catalogue_keys": len(records),
            "underlying_simplifications": len(automorphism_cache),
        },
        "records": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--database-root",
        type=Path,
        required=True,
        help="path containing allr4n04.txt.xz through allr4n09.txt.xz",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "artifacts/generated-results/quartic_hn_rank10_loop_survivors.json"
        ),
    )
    arguments = parser.parse_args()
    result = enumerate_survivors(arguments.database_root.resolve())
    output = arguments.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    labelled = sum(result["counts"]["labelled_survivors_by_loops"].values())
    coloured = sum(result["counts"]["coloured_isomorphism_types_by_loops"].values())
    print(f"QHNW10_LOOP_LABELLED_SURVIVORS={labelled}")
    print(f"QHNW10_LOOP_COLOURED_SURVIVORS={coloured}")
    print(f"QHNW10_LOOP_CATALOGUE_KEYS={result['counts']['catalogue_keys']}")
    print(
        "QHNW10_LOOP_UNDERLYING_SIMPLIFICATIONS="
        f"{result['counts']['underlying_simplifications']}"
    )
    print(f"QHNW10_LOOP_CENSUS_ARTIFACT={output}")


if __name__ == "__main__":
    main()
