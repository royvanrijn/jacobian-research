#!/usr/bin/env python3
"""Replay ONLY the published Xns+(121) finite localization obstructions.

This does not certify the global Selmer bound, arithmetic interpretations,
p-adic precision, lifting lemma, or complete rational-point theorem.
"""
import argparse
import ast
from hashlib import sha1, sha256
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT/"elliptic-curves/data/external/xns121"
BLOBS = {
    "even_selmer_0.m": "170a6c803a136f7bb0b1a20e0a9cd94f1897bbea",
    "even_selmer_2.m": "9a5a6bd26724066b1df2613a3d5044b829933707",
    "even_selmer_4.m": "6e6a3af97b0b483c09c719c146b6e7f68bcbe32f",
    "even_selmer_m2.m": "5815d2135690e896dc1bc6592397bf9bfc34696c",
    "odd_formal_-1_p160_data.m": "f6d404a1e51c88591901e8d392a5981bea330f04",
    "odd_formal_1_p160_data.m": "39219d12268400bdd13b9d15e14a3535b8d48e34",
    "odd_formal_3_p160_data.m": "88fa7a6b242a567b9e4b9f7731df0e244525bd0a",
}


def require(ok, reason):
    if not ok:
        raise ValueError(reason)


def rank(rows, prime=11):
    require(bool(rows) and bool(rows[0]), "empty matrix")
    width = len(rows[0])
    require(all(len(r) == width and all(type(x) is int for x in r) for r in rows), "invalid matrix")
    a = [[x % prime for x in r] for r in rows]
    count = 0
    for j in range(width):
        pivot = next((k for k in range(count, len(a)) if a[k][j]), None)
        if pivot is None:
            continue
        a[count], a[pivot] = a[pivot], a[count]
        inv = pow(a[count][j], -1, prime)
        a[count] = [x*inv % prime for x in a[count]]
        for k in range(count+1, len(a)):
            s = a[k][j]
            a[k] = [(x-s*y) % prime for x, y in zip(a[k], a[count])]
        count += 1
    return count


def parse_assignments(text):
    # A data-only parser: never eval Magma or source code.
    text = re.sub(r"//[^\n]*", "", text)
    matches = list(re.finditer(r"([A-Za-z_]\w*)\s*:=\s*(.*?);", text, flags=re.S))
    values = {}
    for m in matches:
        require(m[1] not in values, "duplicate assignment")
        value = ast.literal_eval(m[2].strip())
        require(isinstance(value, (int, list)), "unsupported data")
        values[m[1]] = value
    rest = re.sub(r"([A-Za-z_]\w*)\s*:=\s*(.*?);", "", text, flags=re.S)
    require(not rest.strip(), "unparsed content")
    return values


def check_even(d, expected):
    S, U, T = d["SelmerLocalRows"], d["FormalRows"], d["LocalTorsionRows"]
    require(len(S) == expected and rank(S) == expected, "global row independence")
    require(all(len(r) == 36 for r in S+U+T) and len(U) == len(T) == 3, "even matrix dimensions")
    rows = []
    for j in (0, 1):
        s = [r[12*j:12*j+12] for r in S]
        u, t = U[j][12*j:12*j+12], T[j][12*j:12*j+12]
        require(rank(s) == expected, "localization not injective")
        require(rank([u]) == rank([t]) == 1 and rank([u, t]) == 2, "formal/torsion lines")
        require(rank(s+[t]) == expected, "torsion not in Selmer localization")
        require(rank(s+[u]) == expected+1, "formal direction not excluded")
        # dim(S intersect (U+T)) = dim S + dim(U+T) - dim(S+U+T).
        require(expected+2-rank(s+[u, t]) == 1, "intersection exceeds torsion")
        rows.append({"place_index": j+1, "selmer_dimension": expected,
                     "localization_injective": True, "formal_direction_excluded": True,
                     "intersection_with_formal_plus_torsion_dimension": 1})
    return rows


def check_odd(d):
    require(d["SelmerDimension"] == 1, "odd global dimension")
    S, U = d["SelmerRamifiedRows"], d["DirectionRows"]
    require(len(S) == 1 and len(U) == 2 and all(len(r) == 44 for r in S+U), "odd dimensions")
    rows = []
    for j in (0, 1):
        s, u = S[0][22*j:22*j+22], U[j][22*j:22*j+22]
        require(rank([s]) == rank([u]) == 1 and rank([s, u]) == 2, "odd direction not excluded")
        rows.append({"place_index": j+1, "selmer_dimension": 1,
                     "localization_injective": True, "formal_direction_excluded": True})
    return rows


def verify():
    results = []
    for name, blob in BLOBS.items():
        b = (DATA/name).read_bytes()
        require(sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest() == blob, "source blob changed")
        d = parse_assignments(b.decode())
        if name.startswith("even"):
            expected = 5 if name in ("even_selmer_0.m", "even_selmer_m2.m") else 4
            checks = check_even(d, expected)
        else:
            checks = check_odd(d)
        results.append({"file": name, "git_blob": blob, "sha256": sha256(b).hexdigest(), "checks": checks})
    require(sum(len(r["checks"]) for r in results) == 14, "coverage")
    return {"schema": "xns121-localization-replay-v1", "status": "PASS",
            "rational_residue_classes": 7, "class_place_checks": 14, "results": results,
            "boundary": "Only finite F11 matrices replayed independently. Global Selmer completeness, arithmetic identification and formal lifting remain external, unverified dependencies; full rational-point theorem is not certified."}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    result = verify()
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True)
            stream.write("\n")
    print(json.dumps({"status": result["status"], "class_place_checks": result["class_place_checks"]}))
