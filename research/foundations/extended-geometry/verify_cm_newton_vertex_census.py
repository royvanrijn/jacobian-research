#!/usr/bin/env python3
"""Verify the historical full numerical census by integral convex paths.

Endpoint counts suffice for feasibility and minimum p-rank, but full counts
also require the rounded profiles of the other convex paths. This transfer
uses integer vertices and exact multiset weights, independently of
the reduced-slope-run enumerator that produced the retained census.
"""

import argparse
from collections import Counter, defaultdict
from functools import cache
from itertools import combinations_with_replacement, groupby
import json
from math import comb
from pathlib import Path


CERTIFICATE = Path(__file__).with_name("cm_newton_polygon_census.json")


def convex_paths(height, dimension):
    """Generate all strictly convex turn sequences with integral vertices."""
    def extend(x, y, previous, path):
        if x == height:
            yield path
            return
        for new_x in range(x + 1, height + 1):
            width = new_x - x
            low = dimension if new_x == height else max(y, dimension - (height - new_x))
            high = min(dimension, y + width)
            for new_y in range(low, high + 1):
                rise = new_y - y
                if previous and rise * previous[0] <= previous[1] * width:
                    continue
                if new_x < height and rise * (height - new_x) >= (dimension - new_y) * width:
                    continue
                yield from extend(new_x, new_y, (width, rise), (*path, (new_x, new_y)))

    yield from extend(0, 0, None, ((0, 0),))


def dual_vertices(path, height, dimension):
    return tuple((height - x, height - x - dimension + y) for x, y in reversed(path))


def rounded_prefixes(path):
    result = [0]
    for (x, y), (stop_x, stop_y) in zip(path, path[1:]):
        width = stop_x - x
        rise = stop_y - y
        for i in range(1, width + 1):
            result.append(y + (rise * i + width - 1) // width)
    return result


def endpoints(path):
    zero = path[1][0] if path[1][1] == 0 else 0
    before, end = path[-2:]
    one = end[0] - before[0] if end[0] - before[0] == end[1] - before[1] else 0
    return zero, one


@cache
def profiles(height, dimension):
    counts = Counter()
    minima = {}
    for path in convex_paths(height, dimension):
        dual = dual_vertices(path, height, dimension)
        if 2 * dimension == height and path > dual:
            continue
        indices = range(1, height // 2 + 1)
        prefixes = rounded_prefixes(path)
        paired = rounded_prefixes(dual)
        delta = tuple(prefixes[j] + paired[j] - j + min(j, dimension) for j in indices)
        z, o = endpoints(path)
        expected = tuple(int(z < j < height - o) + int(o < j < height - z) for j in indices)
        assert all(a >= b for a, b in zip(delta, expected))
        counts[delta] += 1
        minima[delta] = min(minima.get(delta, height), z + o)
    return tuple((delta, counts[delta], minima[delta]) for delta in sorted(counts))


def transfer(height, pairs):
    indices = range(1, height // 2 + 1)
    options = [(m, size, rank, tuple(min(j, m) - value for j, value in zip(indices, delta)))
               for m in range(height // 2 + 1) for delta, size, rank in profiles(height, m)]
    rows = defaultdict(lambda: [0, 0, height * pairs])
    for group in combinations_with_replacement(range(len(options)), pairs):
        choices = [options[index] for index in group]
        signature = tuple(m for m, _, _, _ in choices)
        weight = 1
        for index, repeated in groupby(group):
            multiplicity = sum(1 for _ in repeated)
            weight *= comb(options[index][1] + multiplicity - 1, multiplicity)
        row = rows[signature]
        row[0] += weight
        if not any(signature) or all(sum(column) >= 1 for column in zip(*(item[3] for item in choices))):
            row[1] += weight
            row[2] = min(row[2], sum(item[2] for item in choices))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path, default=CERTIFICATE)
    args = parser.parse_args()
    if not __debug__:
        parser.error("run without -O so the verification assertions execute")
    data = json.loads(args.certificate.read_text())
    rows = data["rows"]
    expected_keys = [(h, e, m) for h in range(1, data["max_dimension"] + 1)
                     for e in range(1, data["max_dimension"] // h + 1)
                     for m in combinations_with_replacement(range(h // 2 + 1), e)]
    assert len(expected_keys) == len(rows), "wrong number of signature rows"
    for (h, e), block in groupby(zip(expected_keys, rows), key=lambda entry: entry[0][:2]):
        computed = transfer(h, e)
        for key, (_, _, signature, *expected) in block:
            assert key == (h, e, tuple(signature)), "wrong signature key"
            actual = computed.pop(tuple(signature))
            assert actual == expected, (h, e, signature, actual, expected)
        assert not computed, (h, e, "extra signature rows")
    assert transfer(2, 3)[1, 1, 1] == [4, 2, 4]
    assert transfer(4, 1)[2,] == [4, 1, 4]
    print(f"Verified all {len(rows)} full rows through dimension {data['max_dimension']}.")
    print("Counts before and after both roundings, and minimum p-ranks, all agree.")
    print(f"Totals: {sum(row[3] for row in rows):,} candidates, {sum(row[4] for row in rows):,} survivors.")


if __name__ == "__main__":
    main()
