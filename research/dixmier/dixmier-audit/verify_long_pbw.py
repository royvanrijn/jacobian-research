#!/usr/bin/env python3
"""Exact Long map and uncorrected normal PBW audit; standard library only."""
from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path

import tempfile as _check_tempfile

_CHECK_OUTPUT_DIRECTORY = None

def _check_output(original):
    """Write new runs to a temporary directory; preserve deposited receipts."""
    global _CHECK_OUTPUT_DIRECTORY
    if _CHECK_OUTPUT_DIRECTORY is None:
        _CHECK_OUTPUT_DIRECTORY = Path(_check_tempfile.mkdtemp(prefix="math-check-"))
        print("Verification output:", _CHECK_OUTPUT_DIRECTORY)
    return _CHECK_OUTPUT_DIRECTORY / Path(original).name


import hashlib
import itertools
import json
import platform
import time

N = 4
ZERO = (0,) * N

def const(c):
    return {ZERO: F(c)} if c else {}

def var(i):
    e = [0] * N
    e[i] = 1
    return {tuple(e): F(1)}

def add(*ps):
    out = {}
    for p in ps:
        for e, c in p.items():
            out[e] = out.get(e, F(0)) + c
    return {e: c for e, c in out.items() if c}

def scale(p, c):
    return {e: v * c for e, v in p.items() if v * c}

def mul(p, q):
    out = {}
    for e, c in p.items():
        for f, d in q.items():
            g = tuple(a + b for a, b in zip(e, f))
            out[g] = out.get(g, F(0)) + c * d
    return {e: c for e, c in out.items() if c}

def power(p, n):
    out = const(1)
    for _ in range(n):
        out = mul(out, p)
    return out

def diff(p, i):
    out = {}
    for e, c in p.items():
        if e[i]:
            f = list(e)
            f[i] -= 1
            out[tuple(f)] = c * e[i]
    return out

def poisson(p, q):
    return add(*[r for a, b in ((2, 0), (3, 1)) for r in
                 (mul(diff(p, a), diff(q, b)),
                  scale(mul(diff(p, b), diff(q, a)), -1))])

def prod(xs):
    out = 1
    for x in xs:
        out *= x
    return out

def evaluate(p, point):
    return sum((c * prod(F(v) ** k for v, k in zip(point, e))
                for e, c in p.items()), F(0))

def pbw(p, q):
    """Product with [p,x]=[z,q]=h; keys are (x,q,p,z,h) exponents."""
    out = {}
    for e, c in p.items():
        for f, d in q.items():
            for i in range(min(e[2], f[0]) + 1):
                ci = comb(e[2], i) * factorial(f[0]) // factorial(f[0] - i)
                for j in range(min(e[3], f[1]) + 1):
                    cj = comb(e[3], j) * factorial(f[1]) // factorial(f[1] - j)
                    g = (e[0]+f[0]-i, e[1]+f[1]-j,
                         e[2]+f[2]-i, e[3]+f[3]-j, i+j)
                    out[g] = out.get(g, F(0)) + c*d*ci*cj
    return {e: c for e, c in out.items() if c}

def digest(p):
    rows = [[list(e), str(c)] for e, c in sorted(p.items())]
    return hashlib.sha256(json.dumps(rows, separators=(",", ":")).encode()).hexdigest()

def require(condition, label):
    if not condition:
        raise AssertionError(label)

def long_images():
    """Return expanded classical images in source order x,q,p,z."""
    x, q, p, z = (var(i) for i in range(N))
    one = const(1)
    a = add(one, scale(mul(x, q), -3))
    beta = add(scale(mul(power(x, 2), p), 3), scale(mul(a, z), 2), scale(power(q, 2), -9))
    y = add(q, scale(mul(x, beta), -F(1, 3)))
    u = mul(x, y)
    w = add(one, u)
    R = add(scale(x, 2), scale(mul(power(x, 2), q), -3))
    S = add(y, scale(mul(mul(x, power(w, 2)), beta), 3),
            scale(mul(mul(x, power(y, 2)), add(const(4), scale(u, 3))), 3))
    T = scale(add(mul(power(w, 3), beta),
                  mul(mul(power(y, 2), w), add(const(4), scale(u, 3)))), -F(1, 2))
    H = add(scale(mul(power(y, 4), add(scale(power(u, 2), 18), scale(u, 78), const(125))), F(1,20)),
            scale(mul(mul(beta, power(y,2)), add(power(u,3), scale(power(u,2),5), scale(u,10), const(-5))), F(3,10)),
            scale(mul(power(beta,2), add(scale(u,9),const(2))), -F(1,6)),
            scale(mul(power(x,2),power(beta,3)), -F(1,6)))
    D = add(scale(mul(add(one,scale(mul(x,q),3)),p),F(1,2)), scale(mul(power(q,2),z),-3),H)
    return [R, T, D, S]

def main():
    start = time.monotonic()
    x, q, p, z = (var(i) for i in range(N))
    images = long_images()
    labels = ["R", "T", "D", "S"]
    source = [x, q, p, z]
    relations = []
    for i, j in itertools.combinations(range(N), 2):
        expected = poisson(source[i], source[j])
        require(poisson(images[i], images[j]) == expected, "Poisson relation " + labels[i] + labels[j])
        defect = add(pbw(images[i], images[j]), scale(pbw(images[j], images[i]),-1),
                     {e+(1,): -c for e,c in expected.items()})
        require(all(e[4] >= 2 for e in defect), "commutator first order")
        if (labels[i], labels[j]) == ("R", "S"):
            compact = scale(mul(power(x,6), mul(add(scale(mul(q,x),3),const(-2)),
                                              add(scale(mul(q,x),9),const(-4)))), -6)
            require(defect == {e+(2,):c for e,c in compact.items()}, "compact RS defect")
        specialized = {}
        for e,c in defect.items():
            g=e[:4]
            specialized[g]=specialized.get(g,F(0))+c
        specialized={e:c for e,c in specialized.items() if c}
        first=next(iter(sorted(defect.items())),None)
        relations.append({"pair": [labels[i], labels[j]], "poisson": "PASS",
                          "pbw_defect_terms": len(defect),
                          "pbw_h_orders": sorted({e[4] for e in defect}),
                          "pbw_defect_sha256": digest(defect),
                          "pbw_h1_defect_terms": len(specialized),
                          "pbw_h1_defect_sha256": digest(specialized),
                          "witness": [list(first[0]),str(first[1])] if first else None})
    points = [(0,0,F(1,24),-F(1,8)), (1,F(2,3),F(247,96),-F(89,64)),
              (-1,-F(2,3),F(247,96),-F(89,64))]
    target = (F(0),F(1,8),F(0),F(0))
    require(len(set(points))==3, "distinct points")
    for pt in points:
        require(tuple(evaluate(im,pt) for im in images)==target, "collision")
    sizes = {name:{"terms":len(im),"degree":max(map(sum,im)),"sha256":digest(im)}
             for name,im in zip(labels,images)}
    require([(sizes[n]["terms"],sizes[n]["degree"]) for n in labels]
            == [(2,3),(47,15),(139,23),(22,11)], "published size controls")
    require(any(r["pbw_h1_defect_terms"] for r in relations), "uncorrected PBW failure")
    here=Path(__file__).resolve().parent
    root=here.parents[2]
    claims=[]  # Mathematical controls do not depend on registry bookkeeping.
    receipt={"date":"2026-11-14","status":"PASS","python":platform.python_version(),
             "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             "arithmetic":"fractions.Fraction; complete sparse coefficient support",
             "source":"https://arxiv.org/html/2608.23777v1#S3",
             "classical_sizes":sizes,"relations":relations,
             "collision":[[str(t) for t in pt] for pt in points],
             "target":[str(t) for t in target],"claims_reviewed":claims,
             "compact_RS_defect":"-6*h^2*x^6*(3*q*x-2)*(9*q*x-4)",
             "boundary":"Exact identities and uncorrected normal PBW failure only; no general nonliftability, independent replay of archived proofs, or DC2 conclusion.",
             "elapsed_seconds":round(time.monotonic()-start,3)}
    (_check_output(here/"receipt.json")).write_text(json.dumps(receipt,indent=2)+"\n")
    print("PASS: all six Long Poisson relations, three-point collision and source-size controls")
    for r in relations:
        print("PBW",r["pair"],"h orders",r["pbw_h_orders"],"nonzero terms at h=1",r["pbw_h1_defect_terms"])
    print("Receipt:",here/"receipt.json")

if __name__=="__main__":
    main()
