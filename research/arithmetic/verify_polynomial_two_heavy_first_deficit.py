#!/usr/bin/env python3
"""Exact transposition hubs and forest moves for the two-heavy first deficit."""
from collections import Counter, deque
from itertools import combinations, product
from math import factorial
from pathlib import Path
import hashlib
import json
import platform

import verify_polynomial_two_cycle_pair as y

if not __debug__:
    raise RuntimeError("Two-heavy checker requires active assertions")
HERE = Path(__file__).resolve().parent
x, c, m, l, h, u = y.x, y.c, y.m, y.l, y.h, y.u


def components(N, edges):
    adj = [set() for _ in range(N)]
    for a, b in edges:
        adj[a].add(b); adj[b].add(a)
    remaining = set(range(N)); answer = []
    while remaining:
        start = min(remaining); remaining.remove(start)
        seen = {start}; queue = [start]
        for a in queue:
            for b in adj[a] & remaining:
                remaining.remove(b); seen.add(b); queue.append(b)
        answer.append(seen)
    return answer


def forest(state):
    N = len(state[0]); ac, bc = x.supports(state[0]), x.supports(state[1])
    s, p = len(state)-2, len(ac)+len(bc)
    assert s >= 1 and s == p-3
    assert c.index(state[0]) >= 2 and c.index(state[1]) >= 2
    assert all(c.index(T) == 1 for T in state[2:])
    assert sum(c.index(g) for g in state) == N-1
    assert len(c.cycles(m.product(state))) == 1
    edges = []
    for g in state[:2]:
        for cy in c.cycles(g, False):
            edges.extend((cy[0], a) for a in cy[1:])
    cs = components(N, edges)
    R = m.compose(*state[:2])
    assert {frozenset(z) for z in cs} == {frozenset(z) for z in c.cycles(R)}
    assert len(cs) == s+1
    outside = set(range(N))-set.union(*(ac+bc))
    matrix = [[len(a & b) for b in bc] for a in ac]
    assert all(v <= 1 for row in matrix for v in row)
    e = sum(map(sum, matrix))
    assert e == len(outside)+2
    all_edges = edges+[tuple(next(iter(x.supports(T)))) for T in state[2:]]
    assert len(all_edges) == N-1 and len(components(N, all_edges)) == 1
    return dict(N=N, s=s, p=p, ac=ac, bc=bc, matrix=matrix, e=e,
                outside=outside, components=cs, R=R)


def gather_last(state, at):
    assert 2 <= at < len(state)
    word = tuple(range(at+1, len(state)))
    old_others = state[:at]+state[at+1:]
    changed = m.braid(state, word)
    assert changed[:-1] == old_others
    return changed, word


def gather_first(state, at):
    assert 2 <= at < len(state)
    word = tuple(range(at, 2, -1))
    changed = m.braid(state, word)
    assert changed[:2] == state[:2] and changed[2] == state[at]
    return changed, word


def rotate_last(state, target):
    C = m.product(state[:-1]); ends, = x.supports(state[-1])
    chosen = next(z for z in c.cycles(C) if target in z)
    start, = ends & set(chosen)
    point = start
    for k in range(len(C)):
        if point == target:
            break
        point = C[point]
    else:
        raise AssertionError("target not reached")
    count = len(state)-1
    twist = tuple(range(1, count))*count
    word = twist*k
    changed = x.conjugated(m.braid(state, word), m.power(C, k))
    assert changed[:-1] == state[:-1]
    assert changed[-1] == m.conjugate(m.power(C, k), state[-1])
    assert target in next(iter(x.supports(changed[-1])))
    return changed


def hub(state, root):
    """Keep both heavy factors literal; make all T edges meet root."""
    original_heavy = state[:2]
    d = forest(state)
    H = next(z for z in d["components"] if root in z)
    degree = lambda z: sum(bool(next(iter(x.supports(T))) & H) for T in z[2:])
    while degree(state) < d["s"]:
        at = next(i for i in range(2, len(state))
                  if not (next(iter(x.supports(state[i]))) & H))
        old = degree(state)
        state, _ = gather_last(state, at)
        state = rotate_last(state, root)
        assert degree(state) == old+1
    # Once the quotient tree is a star, every other cut side is one heavy component.
    for _ in range(d["s"]):
        at = next((i for i in range(2, len(state))
                   if root not in next(iter(x.supports(state[i])))), None)
        if at is None:
            break
        state, _ = gather_last(state, at)
        state = rotate_last(state, root)
    assert all(root in next(iter(x.supports(T))) for T in state[2:])
    assert state[:2] == original_heavy
    ends = [next(iter(x.supports(T)))-{root} for T in state[2:]]
    assert all(len(z) == 1 for z in ends)
    outer_components = [next(i for i, C0 in enumerate(d["components"]) if z <= C0) for z in ends]
    assert len(set(outer_components)) == d["s"]
    return state


def capture_first(state):
    changed = x.capture(state[:3])+state[3:]
    assert m.product(changed) == m.product(state)
    return changed


def make_outside(state):
    d = forest(state)
    if d["outside"]:
        return state
    assert d["p"] >= 5 and d["e"] == 2
    isolated = next((i for i, row in enumerate(d["matrix"]) if sum(row) == 0), None)
    if isolated is None:
        changed = make_outside(m.braid(state, (1,)))
        result = m.braid(changed, (-1,))
        assert forest(result)["outside"]
        return result
    H = d["ac"][isolated]
    root = min(set(range(d["N"]))-H)
    state = hub(state, root)
    at = next(i for i in range(2, len(state))
              if next(iter(x.supports(state[i]))) & H)
    state, _ = gather_last(state, at)
    Bonly = set.union(*d["bc"])-set.union(*d["ac"])
    assert Bonly
    state = rotate_last(state, min(Bonly))
    state, _ = gather_first(state, len(state)-1)
    result = capture_first(state)
    assert len(forest(result)["outside"]) == 1
    return result


def delete_edge(state, point):
    d = forest(state)
    assert d["outside"] and point in set.union(*d["ac"]) & set.union(*d["bc"])
    outside = min(d["outside"])
    H = next(z for z in d["components"] if point in z)
    root = min(H-{point})
    state = hub(state, root)
    at = next(i for i in range(2, len(state))
              if outside in next(iter(x.supports(state[i]))))
    state, _ = gather_last(state, at)
    state = rotate_last(state, point)
    state, _ = gather_first(state, len(state)-1)
    result = capture_first(state)
    nd = forest(result)
    assert (nd["e"], len(nd["outside"])) == (d["e"]-1, len(d["outside"])-1)
    return result


def leaf_transfer(state, leaf, target_b, slot=0):
    d = forest(state)
    assert d["outside"] and leaf in d["ac"] and target_b in d["bc"]
    ai, bj = d["ac"].index(leaf), d["bc"].index(target_b)
    assert sum(d["matrix"][ai]) == 1
    old_b = d["matrix"][ai].index(1)
    H = next(z for z in d["components"] if leaf <= z)
    assert H-leaf
    root, outside = min(H-leaf), min(d["outside"])
    point, = leaf & d["bc"][old_b]
    state = hub(state, root)
    at = next(i for i in range(2, len(state)) if outside in next(iter(x.supports(state[i]))))
    state, _ = gather_last(state, at)
    state = rotate_last(state, point)
    state, _ = gather_first(state, len(state)-1)
    swap0 = state[2]
    cut = capture_first(state)
    mapped_b = [{swap0[z] for z in B0} for B0 in d["bc"]]
    free = sorted(mapped_b[bj]-set.union(*d["ac"]))
    assert free and 0 <= slot < len(free)
    cut, _ = gather_last(cut, 2)
    cut = rotate_last(cut, free[slot])
    cut, _ = gather_first(cut, len(cut)-1)
    swap1 = cut[2]
    result = capture_first(cut)
    new_b = [{swap1[z] for z in B0} for B0 in mapped_b]
    expected = [row[:] for row in d["matrix"]]
    expected[ai][old_b] = 0; expected[ai][bj] = 1
    assert [[len(A0 & B0) for B0 in new_b] for A0 in d["ac"]] == expected
    nd = forest(result)
    assert result[0] == state[0] and (nd["e"], len(nd["outside"])) == (d["e"], len(d["outside"]))
    return result


def independent_rotation(state, at):
    """Pure loop rotates A,B by just one AB-component cycle, fixing every T."""
    d = forest(state); outside = d["outside"]
    assert len(outside) == 1
    v, = outside
    assert all(v in next(iter(x.supports(T))) for T in state[2:])
    other, = next(iter(x.supports(state[at])))-{v}
    component = next(z for z in d["components"] if other in z)
    L = len(component)
    gathered, g = gather_first(state, at)
    prefix_twist = (1, 2)*3*(L+1)
    inverse_g = tuple(-z for z in g[::-1])
    word = g+prefix_twist+inverse_g+(-1, -1)*(L+1)
    result = m.braid(state, word)
    rotation = tuple(d["R"][z] if z in component else z for z in range(d["N"]))
    assert result[:2] == tuple(m.conjugate(rotation, z) for z in state[:2])
    assert result[2:] == state[2:]
    assert m.labels(len(state), word) == tuple(range(len(state)))
    return result, word


def forest_seed(aa, bb, edges):
    assert len(edges) >= 2 and len(set(edges)) == len(edges)
    ac = [[k for k, (i, j) in enumerate(edges) if i == at] for at in range(len(aa))]
    bc = [[k for k, (i, j) in enumerate(edges) if j == at] for at in range(len(bb))]
    nxt = len(edges)
    for cycles, lengths in ((ac, aa), (bc, bb)):
        for cy, length in zip(cycles, lengths):
            assert len(cy) <= length
            cy.extend(range(nxt, nxt+length-len(cy)))
            nxt = max(nxt, max(cy)+1)
    N = sum(aa+bb)-2
    assert nxt+len(edges)-2 == N
    A, B = m.perm(N, ac), m.perm(N, bc)
    cs = c.cycles(m.compose(A, B))
    root = cs[0][0]
    state = (A, B)+tuple(m.perm(N, [(root, cy[0])]) for cy in cs[1:])
    forest(state)
    return state


def probe():
    domains = [((2, 2, 2), (2, 2)), ((2, 3, 4), (2, 3)),
               ((2, 2, 2, 2), (2, 2)), ((2, 2, 3), (2, 2, 2)),
               ((3, 4, 5, 6), (2, 3, 4)), ((3,), (2, 2, 2, 2)),
               ((2, 2, 2, 2), (3,)), ((2,)*6, (3, 3))]
    rows = []
    for aa, bb in domains:
        s = len(aa)+len(bb)-3
        seed = h.tree_seed((aa, bb), s)
        # Deliberately vary the ordering/topology before each constructive operation.
        two_edges = [(0, 0), (min(1, len(aa)-1), min(1, len(bb)-1))]
        if len(set(two_edges)) != 2:
            raise AssertionError("invalid two-edge fixture")
        states = [seed, forest_seed(aa, bb, two_edges)]
        words = u.colored_generators(tuple(tuple(m.cycle_lengths(z)) for z in seed))
        for k in range(5):
            states.append(m.braid(states[-1], words[(3*k+1) % len(words)]))
        checks = created_outside = deleted_edges = 0
        for raw in states:
            created_outside += not forest(raw)["outside"]
            state = make_outside(raw)
            while forest(state)["e"] > 3:
                fd = forest(state)
                point = min(set.union(*fd["ac"]) & set.union(*fd["bc"]))
                state = delete_edge(state, point)
                deleted_edges += 1
            d = forest(state); v, = d["outside"]
            star = hub(state, v)
            for at in range(2, len(star)):
                independent_rotation(star, at); checks += 1
            for at in range(3, len(star)):
                assert m.braid(star, (at,)*3) == star
                assert m.labels(len(star), (at,)*3) != tuple(range(len(star)))
                checks += 1
            A, B, T = star[:3]
            crossed = m.braid(star, (2, 1))
            assert crossed[:3] == (T, m.conjugate(T, A), m.conjugate(T, B))
            assert m.braid(crossed, (1, 2))[:2] == (crossed[1], crossed[2])
            checks += 1
            for ai, row in enumerate(d["matrix"]):
                if sum(row) != 1:
                    continue
                for bj, B0 in enumerate(d["bc"]):
                    available = len(B0)-sum(R0[bj] for R0 in d["matrix"])+(row[bj] == 1)
                    if available:
                        leaf_transfer(state, d["ac"][ai], B0); checks += 1
            swapped = m.braid(state, (1,))
            sd = forest(swapped)
            for ai, row in enumerate(sd["matrix"]):
                if sum(row) != 1:
                    continue
                for bj, B0 in enumerate(sd["bc"]):
                    available = len(B0)-sum(R0[bj] for R0 in sd["matrix"])+(row[bj] == 1)
                    if available:
                        changed = leaf_transfer(swapped, sd["ac"][ai], B0)
                        forest(m.braid(changed, (-1,)))
                        checks += 1
        rows.append(dict(profiles=[aa, bb], states=len(states), checks=checks,
                         outside_creations=created_outside, deleted_edges=deleted_edges))
        print(json.dumps(rows[-1]), flush=True)
    return rows


def three_edge_graph_controls():
    rows = []
    for a, b in product(range(2, 6), repeat=2):
        possible = list(product(range(a), range(b)))
        for pattern in ("two", "three", "alternating"):
            caps = ([2]*(a+b) if pattern == "two" else [3]*(a+b) if pattern == "three"
                    else [2+(i % 2) for i in range(a+b)])
            def degrees(es):
                out = [0]*(a+b)
                for i, j in es:
                    out[i] += 1; out[a+j] += 1
                return out
            states = {frozenset(es) for es in combinations(possible, 3)
                      if all(d <= cap for d, cap in zip(degrees(es), caps))}
            canonical = frozenset(((0, 0), (1, 0), (1, 1)))
            assert canonical in states
            seen = {canonical}; queue = deque([canonical]); moves = 0
            while queue:
                es = queue.popleft(); deg = degrees(es)
                for i, j in es:
                    if deg[i] == 1:
                        for k in range(b):
                            if k != j and deg[a+k] < caps[a+k]:
                                new = (es-{(i, j)}) | {(i, k)}
                                assert new in states
                                moves += 1
                                if new not in seen:
                                    seen.add(new); queue.append(new)
                    if deg[a+j] == 1:
                        for k in range(a):
                            if k != i and deg[k] < caps[k]:
                                new = (es-{(i, j)}) | {(k, j)}
                                assert new in states
                                moves += 1
                                if new not in seen:
                                    seen.add(new); queue.append(new)
            assert seen == states
            rows.append(dict(A_cycles=a, B_cycles=b, capacities=caps,
                             three_edge_forests=len(states), directed_leaf_moves=moves))
    return rows


def pure_closure(seed):
    r = len(seed)
    words = u.colored_generators(tuple((str(i),) for i in range(r)))
    assert all(m.labels(r, w) == tuple(range(r)) for w in words)
    seen = {c.cyclic_label(seed)}; queue = deque(seen)
    while queue:
        state = queue.popleft()
        for word in words:
            new = c.cyclic_label(m.braid(state, word))
            if new not in seen:
                seen.add(new); queue.append(new)
    return seen, len(words)


def complete_components():
    domains = [((2, 2, 2), (2, 2)), ((2, 2, 2), (2, 3)),
               ((2, 2, 2, 2), (3,)), ((2, 2, 2, 2), (2, 2)),
               ((2, 2, 2), (2, 2, 2))]
    rows = []
    for aa, bb in domains:
        s = len(aa)+len(bb)-3
        N, _, _, passport = u.first_deficit((aa, bb), s)
        seed = h.tree_seed((aa, bb), s)
        orbit, generators = pure_closure(seed)
        total = c.character_total(passport)
        assert total == c.comparison_total(passport) == N*len(orbit)
        placements = c.type_placements(passport)
        assert placements == (s+2)*(s+1)//(1+(aa == bb))
        group = y.primitive_group(seed)
        rows.append(dict(profiles=[aa, bb], simple_factors=s, degree=N,
                         pure_generators=generators, fixed_order_pure_absolute_size=len(orbit),
                         fixed_order_pure_literal_size=total, full_type_placements=placements,
                         full_literal_domain=placements*total, literal_seed=seed, group=group,
                         orbit_sha256=c.digest(orbit)))
        print(json.dumps(dict(profiles=[aa, bb], simple_factors=s, degree=N,
                              pure_absolute_size=len(orbit))), flush=True)
    return rows


def negative_controls():
    state = forest_seed((2, 2, 2), (2, 2), [(0, 0), (1, 0), (1, 1)])
    d = forest(state); v, = d["outside"]
    star = hub(state, v)
    at = 2
    ends, = x.supports(star[at]); other, = ends-{v}
    L = len(next(z for z in d["components"] if other in z))
    gathered, g = gather_first(star, at)
    wrong = g+(1, 2)*3*L+tuple(-z for z in g[::-1])+(-1, -1)*L
    assert m.braid(star, wrong)[2:] != star[2:]
    rejected = 0
    for fn, args in ((forest_seed, ((2, 2, 2), (2, 2), [(0, 0)])),
                     (forest_seed, ((2, 2), (2, 2, 2), [(0, 0), (0, 1), (0, 2)])),
                     (delete_edge, (forest_seed((2, 2, 2), (2, 2), [(0, 0), (1, 1)]), 0))):
        try:
            fn(*args)
        except AssertionError:
            rejected += 1
        else:
            raise AssertionError("invalid forest/edge move accepted")
    return dict(rejected=rejected, incorrect_component_rotation_exponent_rejected=True)


def main():
    controls = probe()
    graph = three_edge_graph_controls()
    orbits = complete_components()
    inputs = [
        "POLYNOMIAL_TWO_HEAVY_FIRST_DEFICIT_COMPONENTS.md",
        "POLYNOMIAL_TWO_CYCLE_PAIR_COMPONENTS.md",
        "POLYNOMIAL_FIRST_DEFICIT_COMPONENT_REDUCTION.md",
        "POLYNOMIAL_DUPLICATED_CYCLE_PATH_NORMALIZATION.md",
        "POLYNOMIAL_SINGLE_CYCLE_TAIL_COMPONENTS.md",
        "POLYNOMIAL_PASSPORT_COMPONENTS.md", "POLYNOMIAL_PASSPORT_LIFT_TRANSPORT.md",
        "verify_polynomial_two_cycle_pair.py", "verify_polynomial_duplicated_cycle_paths.py",
        "verify_polynomial_first_deficit.py", "verify_polynomial_stable_hermite.py",
        "verify_polynomial_phase_interpolation.py", "verify_polynomial_involution_components.py",
        "verify_polynomial_lift_transport.py", "verify_polynomial_passport_components.py",
        "sources/jones-zvonkin-braid-cacti-2002.pdf",
        "sources/goulden-jackson-cacti-1992.pdf",
        "sources/wajnryb-two-special-fibers-1995.pdf",
        "sources/chen-gendron-prado-tahar-isoresidual-curves-2026.pdf"]
    inputs.extend(str(Path(module.__file__).resolve().relative_to(HERE))
                  for module in (c.m, c.h, c.v, c.v.one))
    sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    result = dict(
        schema_version=1, kind="two_heavy_complete_first_deficit", result="PASS",
        python=platform.python_version(), checker_sha256=sha(Path(__file__).resolve()),
        input_sha256={p: sha(HERE/p) for p in dict.fromkeys(inputs)},
        constructive_controls=controls, three_edge_graph_controls=graph,
        complete_component_controls=orbits, complete_passport_count=len(orbits),
        complete_literal_domains=sum(z["full_literal_domain"] for z in orbits),
        negative_controls=negative_controls(),
        scope="Written Z1-Z5 classify the complete positive first deficit with two non-transposition factors "
              "and s=p-3>=1 transpositions. Everyp>=5 has one primitiveS_N component, also under fixed-order pure braids; "
              "p4 is Y's two-cycle classification with its sole quadratic exception. Universal hub and forest moves, "
              "three-edge graph normalization, component-phase loop and labeled one-heavy cabling prove the uniform claim. "
              "Finite controls are not extrapolated. First-deficit k>=3, deeper and zero-simple passports remain open.",
        independent_replay=False, formal_verification=False, external_review=False)
    (HERE/"polynomial_two_heavy_first_deficit_verification.json").write_text(json.dumps(result, indent=2)+"\n")
    print("PASS:hub/edge/phase controls;48 complete three-edge graphs;"+
          str(len(orbits))+" complete pure-braid passports;"+
          str(result["complete_literal_domains"])+" literal tuples represented")


if __name__ == "__main__":
    main()
