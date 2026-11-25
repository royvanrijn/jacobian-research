#!/usr/bin/env python3
"""Count split Newton polygons satisfying the simultaneous CM slope inequalities.

The relative rank is h, the CM degree is twice e, and the dimension is e*h.
A signature records the smaller Hodge dimension in each of the e pairs.
Numerical polygons are quotiented by pair permutation and conjugation within
each pair. No assertion about realization by abelian varieties is used.
"""

import argparse
from bisect import bisect_right
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations_with_replacement, groupby
import json
from math import comb, gcd, lcm, prod
from pathlib import Path


MAX_DIMENSION = 30
CERTIFICATE = Path(__file__).with_name("cm_newton_polygon_census.json")


def interior_slopes(height, unit):
    """Enumerate run multisets for slopes strictly between zero and one."""
    types = [
        (denominator, numerator, unit * numerator // denominator)
        for denominator in range(2, height + 1)
        for numerator in range(1, denominator)
        if gcd(numerator, denominator) == 1
    ]
    denominators = [denominator for denominator, _, _ in types]

    def visit(start, remaining, dimension, slopes):
        yield dimension, tuple(sorted(slopes))
        stop = bisect_right(denominators, remaining)
        for index in range(start, stop):
            denominator, numerator, value = types[index]
            if dimension + numerator <= height // 2:
                yield from visit(
                    index, remaining - denominator, dimension + numerator,
                    slopes + (value,) * denominator,
                )

    yield from visit(0, height, 0, ())


def components(height, unit):
    """Complete the interior runs by arbitrary zero and one slopes."""
    for interior_dimension, interior in interior_slopes(height, unit):
        vacancies = height - len(interior)
        max_ones = min(vacancies, height // 2 - interior_dimension)
        for ones in range(max_ones + 1):
            dimension = interior_dimension + ones
            slopes = (0,) * (vacancies - ones) + interior + (unit,) * ones
            if 2 * dimension == height:
                dual = tuple(unit - value for value in reversed(slopes))
                weight = 1 if slopes == dual else 2
                if slopes > dual:
                    continue
            else:
                weight = 1
            yield dimension, slopes, len(interior), weight


def deficits(slopes, unit):
    """Compute floor(top j slopes) minus ceil(bottom j slopes)."""
    smallest = 0
    largest = 0
    result = []
    for low, high in zip(slopes[:len(slopes) // 2], reversed(slopes)):
        smallest += low
        largest += high
        result.append(largest // unit - (smallest + unit - 1) // unit)
    return tuple(result)


@dataclass
class Profile:
    dimension: int
    deficits: tuple
    count: int = 0
    max_interior: int = 0


def component_profiles(height):
    unit = lcm(*range(1, height + 1))
    profiles = {}
    oriented = Counter()
    for dimension, slopes, interior, weight in components(height, unit):
        deficit = deficits(slopes, unit)
        key = dimension, deficit
        if key not in profiles:
            profiles[key] = Profile(dimension, deficit)
        profile = profiles[key]
        profile.count += 1
        profile.max_interior = max(interior, profile.max_interior)
        oriented[dimension] += weight
    return [profiles[key] for key in sorted(profiles)], oriented


def oriented_component_counts(bound):
    """Independent unbounded-knapsack count before quotienting conjugation."""
    counts = [[0] * (bound + 1) for _ in range(bound + 1)]
    counts[0][0] = 1
    for denominator in range(1, bound + 1):
        for numerator in range(denominator + 1):
            if gcd(numerator, denominator) != 1:
                continue
            for height in range(denominator, bound + 1):
                for dimension in range(numerator, height + 1):
                    counts[height][dimension] += counts[height - denominator][dimension - numerator]
    return counts


def passes(signature, profiles):
    return not any(signature) or all(sum(column) >= 1 for column in zip(*profiles))


def forced_ordinary(height, signature):
    mixed = tuple(dimension for dimension in signature if dimension)
    return (
        not mixed
        or mixed == (1,)
        or height == 4 and mixed == (2,)
        or height == 2 and mixed == (1, 1)
    )


@dataclass
class PolygonCounts:
    before: int = 0
    allowed: int = 0
    max_interior: int = 0


def signature_counts(height, pairs, profiles):
    rows = defaultdict(PolygonCounts)
    components_per_dimension = Counter()
    for profile in profiles:
        components_per_dimension[profile.dimension] += profile.count

    for group in combinations_with_replacement(profiles, pairs):
        signature = tuple(item.dimension for item in group)
        row = rows[signature]
        weight = 1
        for item, copies in groupby(group):
            repeats = sum(1 for _ in copies)
            weight *= comb(item.count + repeats - 1, repeats)
        row.before += weight
        if passes(signature, [item.deficits for item in group]):
            row.allowed += weight
            row.max_interior = max(row.max_interior, sum(item.max_interior for item in group))

    assert len(rows) == comb(pairs + height // 2, height // 2)
    for signature, row in sorted(rows.items()):
        assert row.before == prod(
            comb(components_per_dimension[dimension] + copies - 1, copies)
            for dimension, copies in Counter(signature).items()
        )
        assert 0 < row.allowed <= row.before
        assert (row.allowed == 1) == (row.max_interior == 0)
        assert (row.allowed == 1) == forced_ordinary(height, signature)
        yield [height, pairs, list(signature), row.before, row.allowed, height * pairs - row.max_interior]


def check_examples():
    defective = (0, 1, 1, 2)
    assert deficits(defective, 2) == (1, 0)
    assert passes((2,), [(1,)]) and not passes((2,), [(1, 0)])
    ordinary = (1,)
    supersingular = (-1,)
    assert passes((0,), [(0,)])
    assert not passes((1, 1), [ordinary, supersingular])
    assert passes((1, 1, 1), [ordinary, ordinary, supersingular])
    assert not passes((1, 1, 1), [ordinary, supersingular, supersingular])


def census(bound):
    check_examples()
    oriented = oriented_component_counts(bound)
    rows = []
    for height in range(1, bound + 1):
        profiles, counts = component_profiles(height)
        for dimension in range(height // 2 + 1):
            assert counts[dimension] == oriented[height][dimension], (height, dimension)
        for pairs in range(1, bound // height + 1):
            rows.extend(signature_counts(height, pairs, profiles))
    return {
        "max_dimension": bound,
        "conventions": (
            "Numerical split Newton polygons with integral slope-run endpoints, modulo "
            "permutation and conjugation of pairs. Survival means satisfying the simultaneous "
            "deficit inequalities. The p-rank entry is a minimum over surviving polygons."
        ),
        "columns": [
            "relative_rank", "conjugate_pairs", "smaller_hodge_dimensions",
            "polygons_before_inequalities", "surviving_polygons", "minimum_surviving_p_rank",
        ],
        "rows": rows,
    }


def surviving_polygons(height, signature):
    """Generate exact representatives for a requested signature."""
    unit = lcm(*range(1, height + 1))
    catalogues = {dimension: [] for dimension in signature}
    for dimension, slopes, _, _ in components(height, unit):
        if dimension in catalogues:
            catalogues[dimension].append(slopes)
    signature_groups = sorted(Counter(signature).items())

    def visit(index, blocks):
        if index == len(signature_groups):
            if passes(signature, [deficits(slopes, unit) for slopes in blocks]):
                yield [[str(Fraction(slope, unit)) for slope in slopes] for slopes in blocks]
            return
        dimension, copies = signature_groups[index]
        for group in combinations_with_replacement(catalogues[dimension], copies):
            yield from visit(index + 1, blocks + group)

    yield from visit(0, ())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--json", action="store_true", help="emit the full census")
    modes.add_argument(
        "--polygons", metavar="H:D,D,...",
        help="emit surviving polygons for a relative rank and signature, one per line",
    )
    args = parser.parse_args()
    if not __debug__:
        parser.error("run without -O so the audit assertions execute")
    if args.polygons:
        try:
            relative_rank, dimensions = args.polygons.split(":", 1)
            height = int(relative_rank)
            signature = tuple(sorted(map(int, dimensions.split(","))))
            assert 1 <= height * len(signature) <= MAX_DIMENSION
            assert 0 <= min(signature) and max(signature) <= height // 2
        except (AssertionError, ValueError):
            parser.error("expected H:D,D,... with 0 <= D <= H/2 and dimension at most 30")
        for polygon in surviving_polygons(height, signature):
            print(json.dumps(polygon))
        return

    result = census(MAX_DIMENSION)
    if args.json:
        rows = result.pop("rows")
        print(json.dumps(result, indent=2).removesuffix("\n}") + ',\n  "rows": [')
        print(",\n".join("    " + json.dumps(row) for row in rows))
        print("  ]\n}")
        return
    assert result == json.loads(CERTIFICATE.read_text(encoding="utf-8")), "census certificate changed"
    rows = result["rows"]
    print(f"Verified {len(rows)} signatures through dimension {MAX_DIMENSION}.")
    print(f"Numerical polygons before inequalities: {sum(row[3] for row in rows)}.")
    print(f"Numerical polygons satisfying all inequalities: {sum(row[4] for row in rows)}.")
    print(f"Signatures with only the ordinary polygon surviving: {sum(row[4] == 1 for row in rows)}.")
    print("Independent component counts, signature counts and negative controls passed.")


if __name__ == "__main__":
    main()
