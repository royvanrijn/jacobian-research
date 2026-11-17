#!/usr/bin/env python3
"""Exact finite witnesses for the unrestricted polynomial fiber homotopy.

Standard-library Fraction arithmetic only. The all-orders and all-rank
statements are proved in PROOF.md; this checks the retained seed and explicit
third/fifth-order corrections without solving or enlarging any degree box.
"""
from collections import defaultdict
from fractions import Fraction as Q
from math import comb
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
import json

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "dixmier-residue-transport"

def clean(p):
    return {m: c for m, c in p.items() if c}

def add(*ps):
    out = defaultdict(Q)
    for p in ps:
        for m, c in p.items():
            out[m] += c
    return clean(out)

def scale(p, c):
    return clean({m: c * v for m, v in p.items()})

def mul(p, q):
    out = defaultdict(Q)
    for m, c in p.items():
        for n, d in q.items():
            out[tuple(i + j for i, j in zip(m, n))] += c * d
    return clean(out)

def dz(p, r=1):
    for _ in range(r):
        p = {(i, j, k-1): k*c for (i, j, k), c in p.items() if k}
    return p

def integral(p):
    return {(i, j, k+1): c/Q(k+1) for (i, j, k), c in p.items()}

def delta(p, r=1):
    for _ in range(r):
        out = defaultdict(Q)
        for (i, j, k), c in p.items():
            out[i+1, j, k] += c*(3*i-6*j)
            if j:
                out[i, j-1, k] += 2*j*c
        p = clean(out)
    return p

def pi(p, q, r=1):
    return add(*(scale(mul(delta(dz(p, r-j), j),
                           dz(delta(q, r-j), j)), (-1)**j*comb(r,j))
                 for j in range(r+1)))

def linear(S, T, u, v):
    return add(pi(u,T), pi(S,v))

def homotopy(S, T, F):
    I = integral(F)
    return mul(dz(S),I), mul(dz(T),I)

def read(rows):
    out = {}
    for m, c in rows:
        assert len(m) in (3,4)
        if len(m) == 4:
            assert m[3] == 0
        assert m[1] >= 0 and m[2] >= 0
        out[tuple(m[:3])] = Q(c)
    return out

def serial(p):
    return [[list(m),str(c)] for m,c in sorted(p.items())]

def period(p, lam):
    return sum(c*lam.get(m,Q(0)) for m,c in p.items())

def stats(p):
    return {"terms": len(p), "Z_order": max((m[2] for m in p),default=-1),
            "min_X": min((m[0] for m in p),default=0)}

def main():
    data = json.loads((SOURCE/"inputs.json").read_text())
    S,T,s,t = (read(data[k]) for k in ("S","T","s2","t2"))
    lam = read(json.loads((SOURCE/"combined_cocycle.json").read_text())["functional"])
    assert pi(S,T) == {(0,0,0): Q(1)}
    lower_defect = scale(pi(S,T,3),Q(1,24))
    assert not add(linear(S,T,s,t),lower_defect)
    # A fresh third-order solution from the universal formula.
    fresh = tuple(scale(p,-1) for p in homotopy(S,T,lower_defect))
    assert not add(linear(S,T,*fresh),lower_defect)
    fifth = add(pi(s,t),scale(pi(s,T,3),Q(1,24)),
                scale(pi(S,t,3),Q(1,24)),scale(pi(S,T,5),Q(1,1920)))
    assert period(fifth,lam) == 1
    u,v = (scale(p,-1) for p in homotopy(S,T,fifth))
    image = linear(S,T,u,v)
    assert not add(image,fifth)
    assert period(image,lam) == -1
    assert stats(u)["Z_order"] > 1 or stats(v)["Z_order"] > 0
    # Arbitrarily deep Laurent poles are covered by the proof; these are
    # sign/implementation controls, not an all-pole finite enumeration.
    controls = [{(0,0,0):Q(1)}, {(0,0,7):Q(3,5)},
                {(-100000,2,3):Q(7,11),(5,1,0):Q(-2)}]
    for F in controls:
        assert linear(S,T,*homotopy(S,T,F)) == F
    witness = {"S":serial(S), "T":serial(T), "s":serial(s), "t":serial(t),
               "O5":serial(fifth), "u":serial(u), "v":serial(v)}
    (_check_output(HERE/"escape_witness.json")).write_text(json.dumps(witness,indent=2)+"\n")
    files = [SOURCE/"inputs.json",SOURCE/"combined_cocycle.json",
             HERE/"verify_fiber_homotopy.py",HERE/"escape_witness.json"]
    receipt = {
        "claim":"QFL1",
        "status":"PASS",
        "scope":"Exact retained-seed pair relation through hbar^5; universal written homotopy gives every order. No four-generator finite Weyl endomorphism is asserted.",
        "canonical_pair":True, "retained_lower_equation":True,
        "fresh_third_order_correction":[stats(p) for p in fresh],
        "fifth_order_defect":stats(fifth),
        "escape_correction":[stats(u),stats(v)],
        "fifth_order_equation_zero":True,
        "residue_on_defect":"1", "residue_on_correction_image":"-1",
        "Laurent_controls":len(controls),
        "hashes":{str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
        "coverage_boundary":"The all-orders fiber theorem and the all-rank full-frame de Rham theorem use written proofs. This finite replay is neither a formal verification nor an arithmetic-support or algebraization certificate."
    }
    (_check_output(HERE/"receipt.json")).write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt,indent=2))

if __name__ == "__main__":
    main()
