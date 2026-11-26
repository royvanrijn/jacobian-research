#!/usr/bin/env python3
"""Exact strand restrictions and Wirtinger input for the interleaved 8-plat."""

import argparse
from collections import defaultdict
from math import gcd


def twist_rows():
    changes = {
        1: {2: 3, 6: 3},
        2: {3: 3},
        4: {3: 3},
        6: {3: 3},
        7: {4: -3},
        8: {3: 3},
        10: {5: 3},
        11: {4: 3},
        12: {5: 3},
        33: {2: 3, 6: 3},
    }
    return [
        [(i, changes.get(row, {}).get(i, 4)) for i in range(2 if row % 2 else 1, 8, 2)]
        for row in range(1, 34)
    ]


def trace_components(rows):
    colors = 'AABBAABB'
    labels = list(range(8))
    subwords = defaultdict(list)
    column_exponents = defaultdict(list)
    for row in rows:
        for i, exponent in row:
            column_exponents[i].append(exponent)
            left, right = labels[i - 1:i + 1]
            color = colors[left]
            if color == colors[right]:
                j = 1 + sum(colors[k] == color for k in labels[:i - 1])
                subwords[color].append((j, exponent))
            if exponent % 2:
                labels[i - 1:i + 1] = right, left

    adjacent = [{i ^ 1} for i in range(8)]
    for i in range(0, 8, 2):
        a, b = labels[i:i + 2]
        adjacent[a].add(b)
        adjacent[b].add(a)
    unseen = set(range(8))
    components = []
    while unseen:
        stack = [min(unseen)]
        component = set()
        while stack:
            vertex = stack.pop()
            if vertex not in unseen:
                continue
            unseen.remove(vertex)
            component.add(vertex + 1)
            stack.extend(adjacent[vertex])
        components.append(sorted(component))
    return components, dict(subwords), [gcd(*column_exponents[i]) for i in (2, 4, 6)]


def wirtinger_relators(rows):
    """Signed letters for the disk-meridian presentation from the oriented height."""
    labels = [1, -1, 2, -2, 3, -3, 4, -4]
    n = 4
    relators = []
    for row in rows:
        for i, exponent in row:
            for _ in range(abs(exponent)):
                a, b = labels[i - 1:i + 1]
                n += 1
                if exponent > 0:
                    relators.append([-n, a, b, -a])
                    labels[i - 1:i + 1] = n, a
                else:
                    relators.append([-n, -b, a, b])
                    labels[i - 1:i + 1] = b, n
    relators.extend(labels[i:i + 2] for i in range(0, 8, 2))
    return n, relators


def write_gap_trial(n, relators):
    print(f'F:=FreeGroup({n});; U:=GeneratorsOfGroup(F);;')
    print(f'w:={relators};;')
    print('R:=List(w, w0 -> Product(List(w0, x -> U[AbsInt(x)]^SignInt(x))));;')
    print(r'''
Print("GAP ", GAPInfo.Version, "\n");
for s in [[1,2,4],[2,3,4],[1,2,3,4]] do
    C:=CosetTableFromGensAndRels(U,R,U{s}:silent:=true,max:=20000);;
    if C=fail then
        Print(s," cap reached\n");
    else
        Print(s," index ",Length(C[1]),"\n");
    fi;
od;
QUIT;
''')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gap-input', action='store_true', help='print the bounded GAP trial')
    args = parser.parse_args()
    rows = twist_rows()
    assert len(rows) == 33
    assert all(abs(e) >= 3 for row in rows for _, e in row)
    components, subwords, even_gcds = trace_components(rows)
    assert components == [[1, 2, 5, 6], [3, 4, 7, 8]]
    assert subwords == {'A': [(2, 4), (2, -3), (3, 4)], 'B': [(1, 4), (1, 4), (2, 3)]}
    assert even_gcds == [1, 1, 1]
    n, relators = wirtinger_relators(rows)
    if args.gap_input:
        write_gap_trial(n, relators)
    else:
        print(f'components={components}; subbraids={subwords}; even-column gcds={even_gcds}')
        print(f'{len(rows)} twist rows; Wirtinger input {n} generators and {len(relators)} relators')
