#!/usr/bin/env python3
"""Exact retrospective geometry of 37 retained bisections. No searches."""
import argparse
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
import json
from math import gcd, lcm
from pathlib import Path
import retrospective as r
from cover_experiment import evaluate, mul, sub, sqrtq

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE / "BRANCH_BLOCK_PROTOCOL.json"
INPUT = r.OUT / "rank_jump_branch_block_inputs_v1.json"
OUTPUT = r.OUT / "rank_jump_branch_blocks_v1.json"
CHECKPOINT = r.ROOT / "artifacts/local/rank-jump-branch-blocks-v1"


def bindings():
    return {p.name: r.digest(p.read_bytes()) for p in
            (PROTOCOL, Path(__file__), HERE / "retrospective.py",
             HERE / "cover_experiment.py")}


def capture():
    protocol = r.read(PROTOCOL)
    raw = (r.ROOT / protocol["census"]).read_bytes()
    census = json.loads(raw)
    atlas_raw = (r.ROOT / protocol["atlas"]).read_bytes()
    assert r.digest(atlas_raw) == protocol["atlas_sha256"]
    assert census["generation"]["inputs"][protocol["atlas"]] == protocol["atlas_sha256"]
    fibres = census["fibres"]
    assert len(fibres) == 5 and sum(len(f["hits"]) for f in fibres) == 37
    labels = {h["label"] for f in fibres for h in f["hits"]}
    atlas = json.loads(atlas_raw)
    selected = [b for b in atlas["bisections"] if b["label"] in labels]
    assert len(selected) == len(labels)
    assert all(b["residual_chord"]["construction_chart"] == "finite" for b in selected)
    # Re-reading closes a concurrent-mutation window without touching the source.
    assert (r.ROOT / protocol["census"]).read_bytes() == raw
    r.write_new(INPUT, {
        "schema": "rank-jump.branch-block-inputs.v1",
        "bindings": bindings(), "census_sha256": r.digest(raw),
        "atlas_sha256": r.digest(atlas_raw),
        "fibres": [{"parameter": f["parameter"], "hits": f["hits"],
                    "public_complement": f["public_complement"],
                    "rank_result": f["rank_result"],
                    "split_class_span": f["split_class_span"]} for f in fibres],
        "covers": sorted(selected, key=lambda b: b["label"]),
        "scope": "Oracle-selected retrospective diagnostics; exact relation certificates inherited from pinned census."
    })


def factors(q):
    """Monic irreducible factors of a squarefree polynomial of degree 1 or 2."""
    q = list(map(F, q))
    assert len(q) in (2, 3) and q[-1]
    if len(q) == 2:
        return [(q[0] / q[1], F(1))]
    c, b, a = q
    disc = b*b - 4*a*c
    assert disc
    root = sqrtq(disc)
    if root is None:
        return [(c/a, b/a, F(1))]
    return sorted([((b-root)/(2*a), F(1)), ((b+root)/(2*a), F(1))])


def character_data(polynomials):
    fs = [factors(q) for q in polynomials]
    universe = sorted(set(f for row in fs for f in row))
    rows = [sum(1 << universe.index(f) for f in row) for row in fs]
    infinity = [int((len(q)-1) % 2) for q in polynomials]
    if any(infinity):
        rows = [v | (bit << len(universe)) for v, bit in zip(rows, infinity)]
    n = r.rank(rows)
    b = sum(len(f)-1 for f in universe) + int(any(infinity))
    genus = F(1) + F(2**n, 4) * (b-4)
    assert n > 0 and genus.denominator == 1 and genus >= 0
    return {
        "geometric_character_rank": n,
        "branch_points": b,
        "cover_degree": 2**n,
        "genus": int(genus),
        "monic_branch_factors_ascending": [list(map(str, f)) for f in universe],
        "factor_degrees": [len(f)-1 for f in universe],
        "branch_character_masks": rows,
        "infinity_branched": bool(any(infinity)),
    }


def qrank(rows):
    rows = [list(map(F, row)) for row in rows]
    if not rows:
        return 0
    pivot = 0
    for j in range(len(rows[0])):
        i = next((i for i in range(pivot, len(rows)) if rows[i][j]), None)
        if i is None:
            continue
        rows[pivot], rows[i] = rows[i], rows[pivot]
        d = rows[pivot][j]
        rows[pivot] = [x/d for x in rows[pivot]]
        for i in range(pivot+1, len(rows)):
            d = rows[i][j]
            rows[i] = [x-d*y for x,y in zip(rows[i], rows[pivot])]
        pivot += 1
    return pivot


def solve_coordinates(fibre):
    """Solve exact stored relations for split points modulo generic points."""
    q = fibre["public_complement"]["dimension"]
    k = len(fibre["hits"])
    relation = fibre["rank_result"]["relation_basis"]
    assert relation["all_relations_verified_by_exact_group_addition"] is True
    ordered = relation["ordered_points"]
    assert ordered[:17] == [f"generic-P{i}" for i in range(1, 18)]
    assert ordered[17:17+q] == fibre["public_complement"]["ordered_basis_labels"]
    assert ordered[17+q:] == ["split-"+h["label"] for h in fibre["hits"]]
    rr = [list(map(F, row)) for row in relation["relations"]]
    assert len(rr) == k and all(len(row) == 17+q+k for row in rr)
    # C * (split points) = -D * (public complement), modulo generic points.
    augmented = [row[17+q:] + [-x for x in row[17:17+q]] for row in rr]
    for j in range(k):
        i = next(i for i in range(j, k) if augmented[i][j])
        augmented[j], augmented[i] = augmented[i], augmented[j]
        d = augmented[j][j]
        augmented[j] = [x/d for x in augmented[j]]
        for i in range(k):
            if i != j:
                d = augmented[i][j]
                augmented[i] = [x-d*y for x,y in zip(augmented[i], augmented[j])]
    coords = [row[k:] for row in augmented]
    assert all(sum(rr[i][17+q+j]*coords[j][s] for j in range(k)) ==
               -rr[i][17+s] for i in range(k) for s in range(q))
    return coords


def specialization_blocks(coords):
    """Projective direction groups and the exact kernel modulo generic M_Q."""
    groups = {}
    for i, row in enumerate(coords):
        scale = next((v for v in row if v), F(1))
        key = tuple(str(v/scale) for v in row)
        groups.setdefault(key, []).append(i)
    a = [list(row) for row in zip(*coords)]
    k = len(coords)
    pivots = []
    for j in range(k):
        p = len(pivots)
        i = next((i for i in range(p, len(a)) if a[i][j]), None)
        if i is None:
            continue
        a[p], a[i] = a[i], a[p]
        d = a[p][j]
        a[p] = [x/d for x in a[p]]
        for i in range(len(a)):
            if i != p:
                d = a[i][j]
                a[i] = [x-d*y for x,y in zip(a[i], a[p])]
        pivots.append(j)
    kernel = []
    for j in range(k):
        if j in pivots:
            continue
        v = [F(0)]*k
        v[j] = F(1)
        for i,p in enumerate(pivots):
            v[p] = -a[i][j]
        den = lcm(*(x.denominator for x in v))
        ints = [int(x*den) for x in v]
        divisor = gcd(*ints)
        ints = [x//divisor for x in ints]
        assert all(sum(F(ints[i])*coords[i][s] for i in range(k)) == 0
                   for s in range(len(coords[0])))
        kernel.append(ints)
    assert qrank(kernel) == k-qrank(coords)
    return {
        "rational_line_groups": [{"normalized_coordinates": list(key), "indices": value,
                                  "zero_direction": all(x == "0" for x in key)}
                                 for key,value in groups.items()],
        "kernel_dimension_modulo_generic": len(kernel),
        "kernel_basis_integer_coefficients": kernel,
        "kernel_scope": "Linear combinations lie in M tensor Q; not asserted to vanish in E or belong integrally to M.",
    }


def analyze(fibre, covers):
    qs = []
    hits = fibre["hits"]
    t = F(fibre["parameter"])
    for h in hits:
        c = covers[h["label"]]
        q = list(map(F, c["residual_chord"]["q_coefficients"]))
        assert len(q) == 3
        quadratic = c["quadratic_cover"]
        a, b, d = [list(map(F, quadratic[key])) for key in
                   ("leading_coefficients", "linear_coefficients", "constant_coefficients")]
        trace_h = list(map(F, c["trace_section"]["h_coefficients"]))
        assert sub(mul(b,b), [4*x for x in mul(a,d)]) == mul(mul(trace_h,trace_h),q)
        # The pinned census evaluates the affine branch polynomial at t.
        value = evaluate(q, t)
        assert value == F(h["q_value"]) != 0
        assert F(h["canonical_positive_square_root"])**2 == value
        qs.append(q)
    geometry = character_data(qs)
    coords = solve_coordinates(fibre)
    exact_rank = qrank(coords)
    finite_rank = r.rank([r.pack(h["finite_quotient_class_modulo_generic_17"]["coordinates_over_f2"]) for h in hits])
    assert exact_rank >= finite_rank
    assert finite_rank == fibre["split_class_span"]["dimension_modulo_generic_17"]
    pairs = []
    for i,j in combinations(range(len(qs)), 2):
        info = character_data([qs[i], qs[j]])
        pairs.append({"indices": [i,j], "character_rank": info["geometric_character_rank"],
                      "branch_points": info["branch_points"], "genus": info["genus"],
                      "specialized_quotient_rank": qrank([coords[i], coords[j]])})
    triple_counts = Counter()
    candidates = []
    for indices in combinations(range(len(qs)), 3):
        info = character_data([qs[i] for i in indices])
        triple_counts[(info["geometric_character_rank"], info["branch_points"], info["genus"])] += 1
        if info["genus"] <= 1:
            candidates.append({"indices": list(indices), "geometry": info,
                               "specialized_quotient_rank": qrank([coords[i] for i in indices])})
    return {
        "parameter": str(t), "labels": [h["label"] for h in hits],
        "geometry": geometry,
        "q_coefficients_ascending": [list(map(str,q)) for q in qs],
        "specialized_quotient_coordinates_over_q": [list(map(str,row)) for row in coords],
        "exact_specialized_quotient_rank": exact_rank,
        "specialization_blocks": specialization_blocks(coords),
        "finite_quotient_rank_lower_bound": finite_rank,
        "known_total_quotient_rank": fibre["public_complement"]["dimension"],
        "pairs": pairs,
        "triple_geometry_counts": [{"character_rank": n, "branch_points": b,
                                    "genus": g, "count": count}
                                   for (n,b,g),count in sorted(triple_counts.items())],
        "low_genus_triples": candidates,
        "rational_point_on_normalization": "Certified by nonzero square values of every quadratic at the retained parameter.",
        "arithmetic_character_rank_equals_geometric": True,
        "relation_scope": "Exact Q-coordinates derived from the previously group-verified relation certificate; whole-curve upper rank remains UNKNOWN.",
    }


def build(check=False):
    inp = r.read(INPUT)
    assert inp["bindings"] == bindings()
    covers = {c["label"]: c for c in inp["covers"]}
    rows = []
    for i, f in enumerate(inp["fibres"]):
        row = analyze(f, covers)
        rows.append(row)
        if not check:
            CHECKPOINT.mkdir(parents=True, exist_ok=True)
            dest = CHECKPOINT / f"case-{i}.json"
            if dest.exists():
                assert r.read(dest) == row
            else:
                r.write_new(dest, row)
    out = {"schema": "rank-jump.branch-blocks.v1",
           "input_sha256": r.digest(INPUT.read_bytes()),
           "bindings": bindings(), "fibres": rows,
           "low_genus_triple_count": sum(len(f["low_genus_triples"]) for f in rows),
           "scope": "Only the 37 retained retrospective hits; no parameter or point searches."}
    if check:
        assert r.read(OUTPUT) == out
        print("PASS exact branch-block replay")
    else:
        r.write_new(OUTPUT, out)
    for row in rows:
        print(row["parameter"], row["geometry"]["geometric_character_rank"],
              row["geometry"]["branch_points"], row["geometry"]["genus"],
              row["exact_specialized_quotient_rank"], row["triple_geometry_counts"])


def verify():
    """Independent Sage polynomial gcd/factor and rational matrix checks."""
    from sage.all import QQ, PolynomialRing, matrix
    inp, out = r.read(INPUT), r.read(OUTPUT)
    assert inp["bindings"] == bindings()
    assert out["input_sha256"] == r.digest(INPUT.read_bytes())
    ring = PolynomialRing(QQ, "t")
    for f, row in zip(inp["fibres"], out["fibres"]):
        polys = [ring([QQ(x) for x in q]) for q in row["q_coefficients_ascending"]]
        assert all(p.degree() == 2 and p.is_squarefree() for p in polys)
        factorizations = [list(p.factor()) for p in polys]
        union = {str(z.monic()): z.degree() for factors_ in factorizations for z,e in factors_ if e % 2}
        assert sum(union.values()) == row["geometry"]["branch_points"]
        assert all(polys[i].gcd(polys[j]).degree() ==
                   4-next(p["branch_points"] for p in row["pairs"] if p["indices"] == [i,j])
                   for i,j in combinations(range(len(polys)), 2))
        rel = matrix(QQ, f["rank_result"]["relation_basis"]["relations"])
        q = f["public_complement"]["dimension"]
        computed = -rel[:,17+q:].inverse() * rel[:,17:17+q]
        assert computed == matrix(QQ, row["specialized_quotient_coordinates_over_q"])
        assert computed.rank() == row["exact_specialized_quotient_rank"]
        print("PASS Sage branch factors, pair gcds and relation matrix", row["parameter"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("capture", "build", "check", "verify"))
    mode = parser.parse_args().mode
    if mode == "capture":
        capture()
    elif mode == "verify":
        verify()
    else:
        build(mode == "check")
