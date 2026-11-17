#!/usr/bin/env python3
"""Signed reflection lifts in any characteristic, with replayable Smith form."""


def relator_matrix(degree, edges, relators):
    """Integer translations of words in involutive affine reflections.

    Permutations are lists of images. A row variable a_i multiplies the
    ordered root e_p-e_q, p<q. Negative word letters equal positive ones.
    """
    count = len(edges)
    rows = []
    for word in relators:
        perm = list(range(degree))
        block = [[0] * count for _ in range(degree)]
        for letter in word:
            index = abs(letter) - 1
            p, q = edges[index]
            assert 0 <= p < q < degree
            block[perm[p]][index] += 1
            block[perm[q]][index] -= 1
            perm[p], perm[q] = perm[q], perm[p]
        if perm != list(range(degree)):
            raise ValueError("word does not project to the identity")
        rows.extend(row for row in block if any(row))
    return rows


def smith_certificate(matrix, columns=None):
    """Return a Smith diagonal and elementary unimodular row/column steps.

    Each step is swap, add an integer multiple of a different row/column,
    or negate a row. A verifier can replay these without a Smith algorithm.
    """
    a = [list(row) for row in matrix]
    m, n = len(a), len(a[0]) if a else columns
    assert n is not None and all(len(row) == n for row in a)
    steps = []

    def swap(axis, i, j):
        if i == j:
            return
        if axis == "r":
            a[i], a[j] = a[j], a[i]
        else:
            for row in a:
                row[i], row[j] = row[j], row[i]
        steps.append([axis, "swap", i, j])

    def add(axis, target, source, multiple):
        if not multiple:
            return
        assert target != source
        if axis == "r":
            a[target] = [x + multiple * y for x, y in zip(a[target], a[source])]
        else:
            for row in a:
                row[target] += multiple * row[source]
        steps.append([axis, "add", target, source, multiple])

    diagonal = []
    for k in range(min(m, n)):
        nonzero = [(abs(a[i][j]), i, j) for i in range(k, m)
                   for j in range(k, n) if a[i][j]]
        if not nonzero:
            break
        _, i, j = min(nonzero)
        swap("r", k, i)
        swap("c", k, j)
        while True:
            # Euclid first down the pivot column, then across the pivot row.
            changed = False
            for i in range(k + 1, m):
                if a[i][k]:
                    add("r", i, k, -(a[i][k] // a[k][k]))
                    if a[i][k]:
                        swap("r", i, k)
                    changed = True
                    break
            if changed:
                continue
            for j in range(k + 1, n):
                if a[k][j]:
                    add("c", j, k, -(a[k][j] // a[k][k]))
                    if a[k][j]:
                        swap("c", j, k)
                    changed = True
                    break
            if changed:
                continue
            offending = next(((i, j) for i in range(k + 1, m)
                              for j in range(k + 1, n)
                              if a[i][j] % a[k][k]), None)
            if offending is None:
                break
            add("r", k, offending[0], 1)
        if a[k][k] < 0:
            a[k] = [-x for x in a[k]]
            steps.append(["r", "neg", k])
        diagonal.append(a[k][k])
    assert all(diagonal[i + 1] % diagonal[i] == 0 for i in range(len(diagonal) - 1))
    assert all(a[i][j] == (diagonal[i] if i == j and i < len(diagonal) else 0)
               for i in range(m) for j in range(n))
    return {"rows": m, "columns": n, "diagonal": diagonal, "steps": steps,
            "characteristic_zero_dimension": n - len(diagonal),
            "maximum_dimension_all_characteristics": n - diagonal.count(1)}


def dimension_in_characteristic(diagonal, columns, prime):
    return columns - sum(x % prime != 0 for x in diagonal)
