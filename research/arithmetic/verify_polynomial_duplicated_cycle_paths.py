#!/usr/bin/env python3
"""Exact path-normalization controls for two duplicated cycles and one T."""
from collections import Counter, deque
from itertools import product
from math import comb, factorial
from pathlib import Path
import hashlib
import json
import platform

import verify_polynomial_first_deficit as u

if not __debug__:
    raise RuntimeError("Duplicated-cycle path checker requires assertions")
HERE = Path(__file__).resolve().parent
c, m, l, h = u.c, u.m, u.l, u.h


def supports(g):
    return [set(cy) for cy in c.cycles(g, False)]


def forest(state):
    A, B, T = state
    N = len(A)
    ac, bc, tc = (supports(g) for g in state)
    assert len(ac) == len(bc) == 2 and len(tc) == 1 and len(tc[0]) == 2
    a, = {len(x) for x in ac}
    b, = {len(x) for x in bc}
    assert a >= 2 and b >= 2 and N == 2*(a+b-1)
    assert c.index(A)+c.index(B)+c.index(T) == N-1
    assert len(c.cycles(m.product(state))) == 1
    edges = []
    for g in (A, B):
        for cy in c.cycles(g, False):
            edges.extend((cy[0], x) for x in cy[1:])
    adjacency = [set() for _ in range(N)]
    for x, y in edges:
        adjacency[x].add(y); adjacency[y].add(x)
    remaining = set(range(N)); components = []
    while remaining:
        start = min(remaining); remaining.remove(start)
        seen = {start}; queue = [start]
        for x in queue:
            for y in adjacency[x] & remaining:
                remaining.remove(y); seen.add(y); queue.append(y)
        components.append(seen)
    assert len(components) == 2 and len(edges) == N-2
    assert all(len(x & tc[0]) == 1 for x in components)
    R = m.compose(A, B)
    assert {frozenset(cy) for cy in c.cycles(R)} == set(map(frozenset, components))
    intersections = [(i, j, next(iter(x & y))) for i, x in enumerate(ac)
                     for j, y in enumerate(bc) if x & y]
    assert all(len(x & y) <= 1 for x in ac for y in bc)
    da = [sum(i == k for i, j, x in intersections) for k in range(2)]
    db = [sum(j == k for i, j, x in intersections) for k in range(2)]
    outside = set(range(N))-set.union(*ac, *bc)
    assert len(intersections) == 2+len(outside) and len(outside) <= 1
    if outside:
        kind = "path"
        assert sorted(da) == sorted(db) == [1, 2]
    elif da == db == [1, 1]:
        kind = "matching"
    elif sorted(da) == [0, 2]:
        kind = "A-star"
        assert db == [1, 1]
    else:
        kind = "B-star"
        assert da == [1, 1] and sorted(db) == [0, 2]
    return dict(a=a, b=b, N=N, ac=ac, bc=bc, da=da, db=db, kind=kind,
                outside=outside, components=components, R=R, edge_count=len(edges))


def distance(g, x, y):
    z = x
    for k in range(1, len(g)+1):
        z = g[z]
        if z == y:
            return k
    raise AssertionError("points not in one cycle")


def path_data(state):
    d = forest(state)
    assert d["kind"] == "path"
    al, ai = (d["da"].index(k) for k in (1, 2))
    bl, bi = (d["db"].index(k) for k in (1, 2))
    p, = d["ac"][al] & d["bc"][bi]
    q, = d["ac"][ai] & d["bc"][bi]
    r, = d["ac"][ai] & d["bc"][bl]
    alpha, beta = distance(state[0], q, r), distance(state[1], p, q)
    assert 1 <= alpha < d["a"] and 1 <= beta < d["b"]
    return dict(d, alpha=alpha, beta=beta, p=p, q=q, r=r,
                al=al, ai=ai, bl=bl, bi=bi)


def conjugated(state, g):
    return tuple(m.conjugate(g, factor) for factor in state)


def power(g, k):
    out = tuple(range(len(g)))
    for _ in range(k):
        out = m.compose(g, out)
    return out


def rotate(state, k):
    """Q_1^(2k), followed by simultaneous relabeling by (AB)^k."""
    assert k >= 0
    Rk = power(m.compose(state[0], state[1]), k)
    actual = conjugated(m.braid(state, (1, 1)*k), Rk)
    assert actual[:2] == state[:2]
    assert actual[2] == m.conjugate(Rk, state[2])
    return actual


def rotate_to(state, start, target):
    R = m.compose(state[0], state[1]); x = start
    for k in range(len(R)):
        if x == target:
            return rotate(state, k)
        x = R[x]
    raise AssertionError("target in another AB component")


def capture(state):
    """Two positive moves on B,T; independently check the literal identity."""
    A, B, T = state
    points, = supports(T)
    fixed = [x for x in points if B[x] == x]
    assert len(fixed) == 1
    v, = fixed; u0, = points-{v}
    assert B[u0] != u0
    new = m.braid(state, (2, 2))
    expected = (A, m.conjugate(T, B),
                m.perm(len(A), [(u0, m.inverse(B)[u0])]))
    assert new == expected and m.product(new) == m.product(state)
    return new


def path_seed(a, b, alpha=1, beta=1):
    assert a >= 2 and b >= 2 and 1 <= alpha < a and 1 <= beta < b
    ai = tuple(range(a)); nxt = a
    bi = (0,)+tuple(range(nxt, nxt+b-1)); nxt += b-1
    p, r = bi[b-beta], ai[alpha]
    al = (p,)+tuple(range(nxt, nxt+a-1)); nxt += a-1
    bl = (r,)+tuple(range(nxt, nxt+b-1)); nxt += b-1
    N = nxt+1
    state = (m.perm(N, [ai, al]), m.perm(N, [bi, bl]), m.perm(N, [(nxt, 0)]))
    d = path_data(state)
    assert (d["alpha"], d["beta"]) == (alpha, beta)
    return state


def matching_seed(a, b, phase):
    assert a >= 2 and b >= 2
    degree = a+b-1
    assert 0 <= phase < degree
    N = 2*degree
    ac, bc = [], []
    for root in (0, degree):
        ac.append(tuple(range(root, root+a)))
        bc.append((root,)+tuple(range(root+a, root+degree)))
    A, B = m.perm(N, ac), m.perm(N, bc)
    y = power(m.compose(A, B), phase)[degree]
    state = A, B, m.perm(N, [(0, y)])
    assert forest(state)["kind"] == "matching"
    return state


def matching_to_path(state):
    d = forest(state)
    assert d["kind"] == "matching"
    Aonly = set.union(*d["ac"])-set.union(*d["bc"])
    Bonly = set.union(*d["bc"])-set.union(*d["ac"])
    for k in range(d["a"]+d["b"]-1):
        candidate = rotate(state, k)
        ends, = supports(candidate[2])
        if len(ends & Aonly) == len(ends & Bonly) == 1:
            result = capture(candidate)
            assert forest(result)["kind"] == "path"
            return result
    raise AssertionError("aligned matching has no opposite-color capture")


def star_to_path(state):
    d = forest(state)
    if d["kind"] == "B-star":
        changed = star_to_path(m.braid(state, (1,)))
        result = m.braid(changed, (-1,))
        assert forest(result)["kind"] == "path"
        return result
    assert d["kind"] == "A-star"
    isolated = d["ac"][d["da"].index(0)]
    ends, = supports(state[2])
    start, = ends-isolated
    Bonly = set.union(*d["bc"])-set.union(*d["ac"])
    result = capture(rotate_to(state, start, min(Bonly)))
    assert forest(result)["kind"] == "path"
    return result


def set_beta(state, target):
    old = path_data(state)
    assert 1 <= target < old["b"]
    outside, = old["outside"]
    ends, = supports(state[2]); start, = ends-{outside}
    cut = capture(rotate_to(state, start, old["p"]))
    d = forest(cut)
    assert d["kind"] == "A-star"
    assert cut[0] == state[0]
    isolated = old["ac"][old["al"]]
    ends, = supports(cut[2]); start, = ends-isolated
    # Retain q on the internal B-cycle and choose a new attachment target.
    y = power(m.inverse(cut[1]), target)[old["q"]]
    assert y not in set.union(*d["ac"])
    result = capture(rotate_to(cut, start, y))
    new = path_data(result)
    assert (new["alpha"], new["beta"]) == (old["alpha"], target)
    assert result[0] == state[0]
    # The B leaf cycle, and every edge outside this chosen B cycle, survives.
    leaf = old["bc"][old["bl"]]
    assert all(result[1][x] == state[1][x] for x in leaf)
    return result


def set_alpha(state, target):
    old = path_data(state)
    assert 1 <= target < old["a"]
    swapped = m.braid(state, (1,))
    data = path_data(swapped)
    assert (data["alpha"], data["beta"]) == (
        old["b"]-old["beta"], old["a"]-old["alpha"])
    changed = set_beta(swapped, old["a"]-target)
    result = m.braid(changed, (-1,))
    new = path_data(result)
    assert (new["alpha"], new["beta"]) == (target, old["beta"])
    return result


def normalization_controls():
    rows = []
    domains = list(product(range(2, 8), repeat=2))+[(2, 17), (4, 13), (7, 19), (16, 17)]
    for a, b in domains:
        matching = 0
        for phase in range(1, a+b-1):
            state = matching_to_path(matching_seed(a, b, phase))
            assert path_data(state)["a"] == a
            matching += 1
        # The phase-zero matching really is the excluded quadratic lift.
        aligned = matching_seed(a, b, 0)
        for k in range(a+b-1):
            trial = rotate(aligned, k)
            ea, eb = (set.union(*supports(trial[i])) for i in (0, 1))
            ends, = supports(trial[2])
            assert not (len(ends & (ea-eb)) == len(ends & (eb-ea)) == 1)
        if max(a, b) < 8:
            distances = list(product(range(1, a), range(1, b)))
        else:
            distances = sorted(set(product((1, a-1), (1, b-1))))
        beta_changes = alpha_changes = star_anchors = 0
        for alpha, beta in distances:
            state = path_seed(a, b, alpha, beta)
            targets = range(1, b) if max(a, b) < 8 else (1, b-1)
            for target in targets:
                set_beta(state, target); beta_changes += 1
            targets = range(1, a) if max(a, b) < 8 else (1, a-1)
            for target in targets:
                set_alpha(state, target); alpha_changes += 1
            normal = set_alpha(set_beta(state, 1), 1)
            nd = path_data(normal)
            assert (nd["alpha"], nd["beta"]) == (1, 1)
            ends, = supports(normal[2]); start, = ends-nd["outside"]
            normal = rotate_to(normal, start, nd["q"])
            label = lambda z: c.cyclic_label(c.v.one.fixed_product(z))
            assert label(normal) == label(path_seed(a, b))
        state = path_seed(a, b)
        d = path_data(state)
        star = capture(rotate_to(state, d["q"], d["p"]))
        fs = forest(star)
        isolated = fs["ac"][fs["da"].index(0)]
        main = next(comp for comp in fs["components"] if comp != isolated)
        ends, = supports(star[2]); start, = ends-isolated
        for target in sorted(main):
            trial = rotate_to(star, start, target)
            star_to_path(trial)
            star_to_path(m.braid(trial, (1,)))
            star_anchors += 2
        rows.append(dict(a=a, b=b, degree=2*(a+b-1), matching_escapes=matching,
                         path_distance_pairs=len(distances), beta_changes=beta_changes,
                         alpha_changes=alpha_changes, star_anchor_controls=star_anchors,
                         aligned_matching_excluded=True))
    return rows


def group_control(seed, primitive, a, b):
    N = len(seed[0]); base, blocks, edges = l.group_from_transposition(seed)
    assert edges == blocks*base*(base-1)//2
    if primitive:
        assert (base, blocks) == (N, 1)
        order, deck, all_blocks = factorial(N), 1, [1, N]
        triple_count = None
    else:
        assert (base, blocks) == (2, a+b-1)
        A, B = tuple(x[0] for x in l.to_wreath(seed, blocks))[:2]
        comm = m.compose(m.compose(A, B), m.compose(m.inverse(A), m.inverse(B)))
        triples = supports(comm)
        assert len(triples) == 1 and len(triples[0]) == 3
        start = tuple(sorted(triples[0])); seen = {start}; queue = deque([start])
        while queue:
            tri = queue.popleft()
            for g in (A, B, m.inverse(A), m.inverse(B)):
                new = tuple(sorted(g[x] for x in tri))
                if new not in seen:
                    seen.add(new); queue.append(new)
        triple_count = len(seen)
        assert triple_count == comb(blocks, 3)
        quotient_order = factorial(blocks)//(2 if a % 2 and b % 2 else 1)
        order, deck, all_blocks = 2**blocks*quotient_order, 2, [1, blocks, N]
    assert c.block_divisors(seed) == all_blocks and c.stabilizer_size(seed) == deck
    return dict(actual_group_order=order, normal_transposition_base=base,
                canonical_block_count=blocks, normal_transposition_edges=edges,
                deck=deck, all_block_count_divisors=all_blocks,
                quotient_all_three_subsets=triple_count)


def complete_components():
    domains = [(2, 2), (3, 3), (4, 4), (4, 5), (4, 7), (5, 5), (5, 6), (6, 6),
               (4, 12), (7, 9)]
    rows = []
    for a, b in domains:
        N, _, _, passport = u.first_deficit(((a, a), (b, b)), 1)
        degree, D = a+b-1, (2*a-1)*(2*b-1)
        placements = 3*(2-(a == b))
        assert placements == c.type_placements(passport)
        total = c.character_total(passport)
        assert total == c.comparison_total(passport) == degree**2*D
        components = []; union = set(); literal_total = 0
        for primitive, raw in ((True, path_seed(a, b)), (False, matching_seed(a, b, 0))):
            seed = c.v.one.fixed_product(raw)
            control = group_control(seed, primitive, a, b)
            orbit, generators = u.absolute_colored_closure(seed)
            assert not union & orbit; union.update(orbit)
            expected = degree*(D-1)//2 if primitive else degree
            assert len(orbit) == expected
            kinds = {forest(state)["kind"] for state in orbit}
            assert kinds == ({"path", "A-star", "B-star", "matching"} if primitive else {"matching"})
            kind_counts = Counter(forest(state)["kind"] for state in orbit)
            expected_counts = (
                {"matching": degree*(degree-1)//2,
                 "A-star": (a-1)*(a+2*b-2)//2,
                 "B-star": (b-1)*(b+2*a-2)//2,
                 "path": (a-1)*(b-1)*(N-1)}
                if primitive else {"matching": degree})
            assert kind_counts == expected_counts
            literal = len(orbit)*N//control["deck"]; literal_total += literal
            components.append(dict(control, component="primitive" if primitive else "imprimitive",
                                   literal_seed=seed, fixed_pattern_absolute_orbit=len(orbit),
                                   full_absolute_orbit=placements*len(orbit),
                                   full_literal_orbit=placements*literal,
                                   colored_generators=generators, core_kinds=sorted(kinds),
                                   core_counts=dict(kind_counts),
                                   orbit_sha256=c.digest(orbit)))
        assert literal_total == total
        rows.append(dict(a=a, b=b, degree=N, full_type_placements=placements,
                         ordered_character_total=total, full_literal_domain=placements*total,
                         components=components))
        print(json.dumps(dict(a=a, b=b, degree=N,
                              fixed_pattern_absolute_sizes=[x["fixed_pattern_absolute_orbit"]
                                                            for x in components])), flush=True)
    return rows


def negative_controls():
    rejected = 0
    for function, args in ((path_seed, (1, 3)), (path_seed, (3, 3, 3, 1)),
                           (matching_seed, (3, 3, 5)),
                           (matching_to_path, (matching_seed(4, 5, 0),))):
        try:
            function(*args)
        except AssertionError:
            rejected += 1
        else:
            raise AssertionError("invalid input or aligned escape accepted")
    return dict(rejected=rejected, aligned_imprimitive_escape_rejected=True)


def main():
    sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    inputs = ["POLYNOMIAL_DUPLICATED_CYCLE_PATH_NORMALIZATION.md",
              "POLYNOMIAL_FIRST_DEFICIT_COMPONENT_REDUCTION.md",
              "POLYNOMIAL_PASSPORT_COMPONENTS.md", "POLYNOMIAL_PASSPORT_LIFT_TRANSPORT.md",
              "verify_polynomial_first_deficit.py", "verify_polynomial_stable_hermite.py",
              "verify_polynomial_phase_interpolation.py", "verify_polynomial_involution_components.py",
              "verify_polynomial_lift_transport.py", "verify_polynomial_passport_components.py",
              "sources/jones-zvonkin-braid-cacti-2002.pdf",
              "sources/goulden-jackson-cacti-1992.pdf",
              "sources/wajnryb-two-special-fibers-1995.pdf"]
    inputs.extend(str(Path(module.__file__).resolve().relative_to(HERE))
                  for module in (c.m, c.h, c.v, c.v.one))
    controls = normalization_controls()
    components = complete_components()
    result = dict(
        schema_version=1, kind="duplicated_cycle_path_normalization", result="PASS",
        python=platform.python_version(), checker_sha256=sha(Path(__file__).resolve()),
        input_sha256={path: sha(HERE/path) for path in dict.fromkeys(inputs)},
        normalization_controls=controls, complete_component_controls=components,
        complete_passport_count=len(components),
        complete_literal_domains=sum(x["full_literal_domain"] for x in components),
        negative_controls=negative_controls(),
        scope="Written X1-X4 classify every a,b>=2 into exactly two full components. "
              "Finite controls check literal Q^2 identities, twists with relabeling, all four core shapes, "
              "opposite-color capture, both independent path-distance changes and unchanged complementary cycles. "
              "Ten complete component closures exhaust independent hook totals. "
              "No general residue-polynomial irreducibility or arbitrary-passport classification is inferred.",
        independent_replay=False, formal_verification=False, external_review=False)
    (HERE/"polynomial_duplicated_cycle_paths_verification.json").write_text(
        json.dumps(result, indent=2)+"\n")
    print("PASS: all four cores and two distance controls; "+
          str(len(components))+" complete passports; "+
          str(result["complete_literal_domains"])+" literal tuples represented")


if __name__ == "__main__":
    main()
