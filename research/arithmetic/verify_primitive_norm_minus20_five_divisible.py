#!/usr/bin/env python3
"""Exact finite support for CN2FiveDivisibleOperator and CN2EffectiveAll."""
if not __debug__:
    raise RuntimeError("This certificate requires Python assertions enabled.")
from collections import Counter
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


def matmul(a,b):
    return ((a[0][0]*b[0][0]+a[0][1]*b[1][0],
             a[0][0]*b[0][1]+a[0][1]*b[1][1]),
            (a[1][0]*b[0][0]+a[1][1]*b[1][0],
             a[1][0]*b[0][1]+a[1][1]*b[1][1]))


def val(n,q):
    if n == 0:
        return 1000
    result = 0
    while n%q == 0:
        n//=q
        result+=1
    return result


def endpoint_distances(g,b,q):
    def distance(rowshift,colshift,detshift):
        smallest = min(val(g[i][j],q)+rowshift[i]+colshift[j]
                       for i in range(2) for j in range(2))
        return b+detshift-2*smallest
    return (distance((0,0),(0,0),0),
            distance((0,0),(0,1),1),
            distance((0,-1),(0,0),-1),
            distance((0,-1),(0,1),0))


def lattice_shells(q=5,maximum=5):
    records = []
    for b in range(maximum+1):
        distances = Counter()
        selected = radical = 0
        for a in range(b+1):
            c = b-a
            for t in range(q**a):
                M = ((q**a,t),(0,q**c))
                for s in list(range(q))+[None]:
                    h = ((1,0),(s,1)) if s is not None else ((0,1),(1,0))
                    g = matmul(M,h)
                    if g[1][0]%q:
                        continue
                    four = endpoint_distances(g,b,q)
                    if (four[0] == four[3] == 0 or
                            four[1] == four[2] == 0):
                        edge_distance = 0
                    else:
                        edge_distance = 1+min(four)
                    assert edge_distance <= b
                    distances[edge_distance] += 1
                    is_radical = g[0][0]%q == 0 and g[1][1]%q == 0
                    if b:
                        assert is_radical == (edge_distance < b)
                    if is_radical:
                        radical += 1
                    else:
                        selected += 1
        expected = {0:1, **{k:2*q**k for k in range(1,b+1)}}
        assert distances == expected
        d = 1+2*sum(q**k for k in range(1,b+1))
        assert sum(distances.values()) == d
        if b:
            assert selected == 2*q**b
            assert radical == 1+2*sum(q**k for k in range(1,b))
        records.append({"b": b, "integral_cosets":d,
                        "edge_distance_histogram":dict(distances),
                        "nonradical_cosets":selected,
                        "radical_cosets":radical})
    return records


def tree_intersections(q=5,radius=5):
    # Build both halves of a rooted edge in a (q+1)-regular tree.
    endpoints = [(0,1)]
    incident = [[0],[0]]
    levels = [0]
    frontier = [(0,0),(1,0)]
    nextvertex = 2
    for depth in range(1,radius+1):
        newfrontier = []
        for v,parentedge in frontier:
            for _ in range(q):
                w = nextvertex
                nextvertex += 1
                k = len(endpoints)
                endpoints.append((v,w))
                levels.append(depth)
                incident[v].append(k)
                incident.append([k])
                newfrontier.append((w,k))
        frontier = newfrontier
    neighbors = []
    for k,(v,w) in enumerate(endpoints):
        neighbors.append(set(incident[v]+incident[w])-{k})
    for k in range(1,radius):
        for i,level in enumerate(levels):
            if level >= radius:
                continue
            actual = sum(levels[j] == k for j in neighbors[i])
            if k == 1:
                expected = (2*q if level == 0 else q-1 if level == 1
                            else 1 if level == 2 else 0)
            else:
                expected = (q if level == k-1 else q-1 if level == k
                            else 1 if level == k+1 else 0)
            assert actual == expected
    return {"tree_edges":len(endpoints), "radius":radius,
            "intersection_recurrence_checked_through":radius-1}


def constants():
    assert all(F(b**6,5**b) < F(7,5)**6 for b in range(1,5))
    assert F(5,4)**6 < 5
    assert 40**12 < 2**64
    upper = F(1,16)+F(117,200)+F(384,2**21)
    assert upper < F(3,4)
    assert F(3,4*486) > F(1,1000)
    return {"b_times_5_to_minus_b_over_6_upper":"7/5",
            "normalized_spectral_error_upper":"21*n^(-1/3)",
            "sieve_upper_at_n_2_64":str(upper),
            "margin_below_three_quarters":str(F(3,4)-upper),
            "claimed_count_denominator":1000}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).with_name('primitive_norm_minus20_five_divisible_checks.json'))
    args = parser.parse_args()
    result = {"checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "lattice_shells":lattice_shells(),
              "tree_intersections":tree_intersections(),
              "constants":constants(),
              "scope":"Finite local lattice cosets through b=5, independent tree intersection counts and exact infinite-inequality constants. Does not itself prove the global or full-space spectral argument."}
    out=args.output
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
