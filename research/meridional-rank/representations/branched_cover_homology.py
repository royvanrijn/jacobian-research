#!/usr/bin/env python3
"""Integral first homology from a transposition-coloured meridional diagram.

Use the complete lifted edge presentation, with each meridian cycle filled.
This implementation does not build the scalar root-lift matrix.
"""

from integer_lifts import smith_certificate


def cover_homology(degree, transpositions, relators):
    generators = len(transpositions)
    edge_count = generators * degree
    parent = list(range(degree))

    def root(v):
        while parent[v] != v:
            v = parent[v]
        return v

    def image(i, v):
        a, b = transpositions[i]
        return b if v == a else a if v == b else v

    tree = []
    for i, (a, b) in enumerate(transpositions):
        if root(a) != root(b):
            parent[root(a)] = root(b)
            tree.append(i * degree + a)
    assert len(tree) == degree - 1, "monodromy is not transitive"
    rows = []
    for word in relators:
        for start in range(degree):
            row = [0] * edge_count
            v = start
            for letter in word:
                i = abs(letter) - 1
                if letter > 0:
                    row[i * degree + v] += 1
                    v = image(i, v)
                else:
                    v = image(i, v)
                    row[i * degree + v] -= 1
            assert v == start
            rows.append(row)
    # The filling slope over a fixed sheet has length one; over a moved
    # pair it has length two. Both sorts of fillings are indispensable.
    for i, (a, b) in enumerate(transpositions):
        row = [0] * edge_count
        row[i * degree + a] = row[i * degree + b] = 1
        rows.append(row)
        for v in range(degree):
            if v not in (a, b):
                row = [0] * edge_count
                row[i * degree + v] = 1
                rows.append(row)
    retained = [j for j in range(edge_count) if j not in tree]
    presentation = [[row[j] for j in retained] for row in rows]
    smith = smith_certificate(presentation, len(retained))
    return {"cover_degree": degree, "tree_edges": tree, "edge_generators_after_tree": len(retained),
            "presentation": presentation, "smith": smith,
            "H1_free_rank": len(retained) - len(smith["diagonal"]),
            "H1_torsion_factors": [d for d in smith["diagonal"] if d > 1]}


def betti_in_prime(data, prime):
    return data["H1_free_rank"] + sum(d % prime == 0 for d in data["H1_torsion_factors"])
