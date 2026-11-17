#!/usr/bin/env python3
"""Oriented SL2 transvection colourings and restricted scalar lift kernels."""

from itertools import product
from verify_integer_lifts import nullspace


def multiply(a, b, p):
    return tuple(sum(a[2*i+k] * b[2*k+j] for k in range(2)) % p for i in range(2) for j in range(2))


def inverse(a, p):
    assert (a[0]*a[3]-a[1]*a[2]) % p == 1
    return (a[3], -a[1] % p, -a[2] % p, a[0])


def vector_image(a, u, p):
    return ((a[0]*u[0]+a[1]*u[1]) % p, (a[2]*u[0]+a[3]*u[1]) % p)


def roots_and_tables(p):
    def canonical(u):
        return min(u, tuple(-x % p for x in u))
    roots = sorted({canonical(u) for u in product(range(p), repeat=2) if any(u)})
    # Normalize first seed to the positive transvection with root (1,0).
    roots.remove(canonical((1, 0)))
    roots.insert(0, canonical((1, 0)))
    index = {u: i for i, u in enumerate(roots)}
    matrices = [((1-a*b) % p, a*a % p, -b*b % p, (1+a*b) % p) for a, b in roots]
    forwards = [[index[canonical(vector_image(inverse(o, p), u, p))] for u in roots] for o in matrices]
    backwards = [[index[canonical(vector_image(o, u, p))] for u in roots] for o in matrices]
    return roots, matrices, forwards, backwards


def oriented_diagram(pd):
    n = len(pd)
    parent = list(range(2*n+1))
    def root(e):
        while parent[e] != e:
            e = parent[e]
        return e
    for a, b, c, d in pd:
        assert c == a % (2*n) + 1
        assert d == b % (2*n) + 1 or b == d % (2*n) + 1
        parent[root(b)] = root(d)
    arc_ids = {r: i for i, r in enumerate(sorted({root(e) for e in range(1,2*n+1)}))}
    arc_of = {e: arc_ids[root(e)] for e in range(1, 2*n+1)}
    assert len(arc_ids) == n
    return [(arc_of[b], arc_of[a], arc_of[c], int(d == b % (2*n) + 1)) for a,b,c,d in pd]


def extend(n, crossings, seeds, assignment, forward, backward):
    labels = [None] * n
    for s, value in zip(seeds, assignment):
        labels[s] = value
    changed = True
    while changed:
        changed = False
        for o,a,b,sign in crossings:
            if labels[o] is None:
                continue
            direct, reverse = (forward, backward) if sign else (backward, forward)
            if labels[a] is not None:
                value = direct[labels[o]][labels[a]]
                if labels[b] is None:
                    labels[b] = value
                    changed = True
                elif labels[b] != value:
                    return None
            elif labels[b] is not None:
                labels[a] = reverse[labels[o]][labels[b]]
                changed = True
    assert all(v is not None for v in labels)
    return labels


def lift_kernel(crossings, labels, roots, matrices, p):
    n = len(labels)
    rows = []
    for o,a,b,sign in crossings:
        word = [-(o+1),a+1,o+1,-(b+1)] if sign else [o+1,a+1,-(o+1),-(b+1)]
        prefix = (1,0,0,1)
        block = [[0]*n for _ in range(2)]
        for letter in word:
            i = abs(letter)-1
            value = vector_image(prefix, roots[labels[i]], p)
            for j in range(2):
                block[j][i] += value[j] if letter > 0 else -value[j]
            base = matrices[labels[i]]
            prefix = multiply(prefix, base if letter > 0 else inverse(base,p), p)
        assert prefix == (1,0,0,1)
        rows.extend(block)
    return nullspace(rows, n, p)
