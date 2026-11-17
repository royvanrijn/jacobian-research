"""Exact controls for coupled restrictions; universal arguments are in the TeX note."""
from collections import deque
from pathlib import Path
import hashlib
import json
import math
import sympy as s

root = Path(__file__).resolve().parent

def swap_edge(q, edge, radical=0):
    a, b = edge
    out = []
    for point in range(4 << q):
        x, z = point >> q, point & ((1 << q)-1)
        moving = x == a or x == b
        xx = b if x == a else a if x == b else x
        zz = z ^ (radical if moving else 0)
        out.append((xx << q) | zz)
    return tuple(out)

def compose(p, q):
    return tuple(p[q[x]] for x in range(len(p)))

def orbit(gens, start=0):
    seen, queue = {start}, deque([start])
    while queue:
        x = queue.popleft()
        for g in gens:
            y = g[x]
            if y not in seen:
                seen.add(y)
                queue.append(y)
    return seen

def commuting_maps(gens):
    n = len(gens[0])
    answers = []
    # Any commuting permutation of a transitive action is determined by c(0).
    for target in range(n):
        values = {0: target}
        queue = deque([0])
        valid = True
        while queue and valid:
            x = queue.popleft()
            for g in gens:
                y, image = g[x], g[values[x]]
                if y in values:
                    if values[y] != image:
                        valid = False
                        break
                else:
                    values[y] = image
                    queue.append(y)
        if valid and len(values) == n and len(set(values.values())) == n:
            answers.append(tuple(values[x] for x in range(n)))
    return set(answers)

group_controls = []
for q in range(4):
    n = 4 << q
    gens = [swap_edge(q,e) for e in [(0,1),(1,2),(2,3)]]
    gens += [swap_edge(q,(0,1),1 << j) for j in range(q)]
    assert len(orbit(gens)) == n
    centralizer = commuting_maps(gens)
    translations = {tuple((x << q) | (z ^ c)
                          for x in range(4) for z in range(1 << q))
                    for c in range(1 << q)}
    assert centralizer == translations
    assert all(compose(c,c) == tuple(range(n)) for c in centralizer)
    for g in gens:
        assert sum(g[x] != x for x in range(n)) == n//2
        assert compose(g,g) == tuple(range(n))
        for point in range(n):
            xx, zz = point >> q, point & ((1 << q)-1)
            if (g[point] >> q) == xx:
                assert g[point] == point  # no motion inside fixed base blocks
    group_controls.append({"radical_dimension":q,"permutation_degree":n,
                           "group_order_formula":24*(1 << (3*q)),
                           "centralizer_size":len(centralizer),
                           "centralizer_exponent":1 if q == 0 else 2,
                           "noncyclic_deck_obstruction":q >= 2,
                           "all_tested_inertias_move_exactly_half":True})

z = s.symbols("z")
def char_dimension(n, trace_counts):
    t = s.symbols("T")
    a, b = s.Integer(1), t
    if n == 0:
        ch = a
    elif n == 1:
        ch = b
    else:
        for _ in range(2,n+1):
            a,b = b,s.expand(t*b-a)
        ch = b
    value = s.simplify(sum(count*ch.subs(t,tr) for tr,count in trace_counts)
                       /sum(count for _,count in trace_counts))
    assert value.is_Integer
    return int(value)

phi = (1+s.sqrt(5))/2
ico_traces = [(2,1),(-2,1),(0,30),(1,20),(-1,20),
              (phi,12),(-phi,12),(phi-1,12),(1-phi,12)]
oct_traces = [(2,1),(-2,1),(0,18),(1,8),(-1,8),
              (s.sqrt(2),6),(-s.sqrt(2),6)]
ico = {n:char_dimension(n,ico_traces) for n in range(20)}
octa = {n:char_dimension(n,oct_traces) for n in range(16)}
assert {n:v for n,v in ico.items() if v} == {0:1,12:1}
assert {n:v for n,v in octa.items() if v} == {0:1,8:1,12:1}
def capacity_lower(n, orbit_size):
    m = n//orbit_size
    return n//2 if m <= 1 else n//2-m
# These are boundary controls. The all-degree inequality is proved in the note.
assert capacity_lower(20,12) == 10
assert capacity_lower(22,12) == 11
assert capacity_lower(24,12) == 10
assert capacity_lower(16,6) == 6
assert capacity_lower(18,6) == 6

orbifold_controls = []
for signature, degree in [((2,2,2),4),((2,2,5),10),
                          ((2,3,3),12),((2,3,4),24),((2,3,5),60)]:
    assert 2*degree-sum(degree-degree//m for m in signature) == 2
    # Pullback fibres split into m reduced components at each of degree/m points.
    euler_cover = 2+sum((degree//m)*(m-1) for m in signature)
    assert euler_cover == 2*degree
    orbifold_controls.append({"signature":signature,"cover_degree":degree,
                              "pullback_euler":euler_cover,
                              "all_pullback_components_reduced":True})
# Bad orbifold boundary: these groups are cyclic and cannot use the good-cover argument.
bad_spindles = [{"orders":[a,b],"group_order":math.gcd(a,b)}
                for a,b in [(2,3),(2,4),(3,6)]]

note = root/"INWARD_COUPLED_RESTRICTIONS_2026_11_16.tex"
dependencies = [
root.parent/"affine-complements/HEREDITARY_EULER_EXCESS_AND_QUOTIENT_DEGREE_BOUNDS.md",
root.parent/"affine-complements/MINIMAL_DEGREE_PSEUDOPLANE_COMPRESSION_AND_DECK_DICHOTOMY.md",
root.parent/"hc4-triangular-graphs/EXACT_VOLUME_AND_COMPLETE_RULING_COMPONENTS.md",
root.parent/"moment-vanishing-audit/LOGARITHMIC_JETS_AND_SPARSE_ALL_RING_SAFETY.md",
root.parent/"meridional-rank/LIFT_DIMENSION_AND_RADICAL_TOWERS.tex"
]
receipt = {
"date":"2026-11-16","status":"finite_exact_controls_pass",
"scope":"Group centralizers and character gaps are exact finite calculations; all-degree conclusions rest on the written arguments and adopted geometric inputs.",
"radical_tower_controls":group_controls,
"icosahedral_nonzero_fixed_dimensions_below_20":{n:v for n,v in ico.items() if v},
"octahedral_nonzero_fixed_dimensions_below_16":{n:v for n,v in octa.items() if v},
"orbifold_cover_controls":orbifold_controls,"bad_spindle_controls":bad_spindles,
"proof_sha256":hashlib.sha256(note.read_bytes()).hexdigest(),
"checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
"dependency_sha256":{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in dependencies},
"external_review":False,"formal_verification":False,"registry_promotion":False,
"novelty_claim":False,"agent_messages_sent":0,
"critical_scope":"H_q exclusion requires every generic branch inertia to be a lifted transposition T_(e,r); other inertia classes are not excluded."
}
out=root/"inward-coupled-restrictions-controls-2026-11-16.json"
out.write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps(receipt,indent=2))
