#!/usr/bin/env python3
"""Separate replay of the pinned release's 243 pseudo-reflection certificates.

No author checker code is imported. Check conjugation by multiplication,
so inversion is unnecessary. Knot-name/PD identity remains source imported.
"""

from collections import Counter
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def column_rank(matrix, prime):
    pivots = {}
    for column in zip(*matrix):
        v = [x % prime for x in column]
        for pivot in sorted(pivots):
            factor = v[pivot]
            v = [(x - factor * y) % prime for x, y in zip(v, pivots[pivot])]
        pivot = next((i for i, x in enumerate(v) if x), None)
        if pivot is not None:
            inverse = pow(v[pivot], -1, prime)
            pivots[pivot] = [x * inverse % prime for x in v]
    return len(pivots)


def multiply(a, b, prime):
    columns = list(zip(*b))
    return [[sum(x * y for x, y in zip(row, col)) % prime for col in columns] for row in a]


def check(c):
    pd, prime, matrices = c["pd"], c["field"], c["matrices"]
    assert prime > 1 and all(prime % d for d in range(2, prime))
    n = 2 * len(pd)
    assert len(matrices) == n and sorted(e for crossing in pd for e in crossing) == sorted(list(range(1, n + 1)) * 2)
    assert all(len(a) == 4 and all(len(row) == 4 for row in a) for a in matrices)
    matrices = [[[v % prime for v in row] for row in a] for a in matrices]
    parent = list(range(n))

    def root(e):
        while parent[e] != e:
            e = parent[e]
        return e

    for a, b, c1, d in pd:
        assert c1 == a % n + 1
        assert d == b % n + 1 or b == d % n + 1
        assert matrices[b - 1] == matrices[d - 1]
        parent[root(b - 1)] = root(d - 1)
        A, O, C = matrices[a - 1], matrices[b - 1], matrices[c1 - 1]
        if d == b % n + 1:
            assert multiply(A, O, prime) == multiply(O, C, prime)
        else:
            assert multiply(O, A, prime) == multiply(C, O, prime)
    assert len({root(i) for i in range(n)}) == len(pd)
    differences = []
    for a in matrices:
        assert column_rank(a, prime) == 4
        delta = [[v - int(i == j) for j, v in enumerate(row)] for i, row in enumerate(a)]
        assert column_rank(delta, prime) == 1
        differences.extend(delta)
    assert column_rank(differences, prime) == 4
    seen = {root(e - 1) for e in c["wirtinger_seeds"]}
    assert len(seen) == 4
    changed = True
    while changed:
        changed = False
        for a, b, c1, d in pd:
            o, x, y = root(b - 1), root(a - 1), root(c1 - 1)
            if o in seen and ((x in seen) != (y in seen)):
                seen.update([x, y])
                changed = True
    assert seen == {root(i) for i in range(n)}


def main():
    source = HERE / "sources/current13/data/certificates.json"
    certs = json.loads(source.read_text())
    assert len(certs) == 243 and len({c["name"] for c in certs}) == 243
    for c in certs:
        check(c)
    # Corruptions must invalidate the hypotheses needed for the proof.
    controls = []
    for kind in ["matrix", "seeds", "identity"]:
        c = json.loads(json.dumps(certs[0]))
        if kind == "matrix":
            c["matrices"][0][0][0] += 1
        elif kind == "seeds":
            c["wirtinger_seeds"] = c["wirtinger_seeds"][:3]
        else:
            c["matrices"] = [[[int(i == j) for j in range(4)] for i in range(4)] for _ in c["matrices"]]
        try:
            check(c)
        except AssertionError:
            controls.append(kind)
        else:
            raise AssertionError(f"corruption accepted: {kind}")
    report = {"status": "PASS", "certificates": len(certs),
              "field_counts": dict(Counter(c["field"] for c in certs)),
              "rejected_corruptions": controls,
              "input_sha256": {str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in [source, Path(__file__), HERE / "sources/current13/manifest.json"]},
              "claim_boundary": "Separate standard-library replay of algebra and four-seed upper certificates passes. PD-to-knot-name identity and the October KnotInfo uncertainty entries are imported from the pinned authors' release, not checked against an independent current KnotInfo copy. Not a full paper audit or priority certification."}
    (HERE / "current13_certificate_verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
