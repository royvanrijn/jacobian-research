"""Small exact controls for the written all-degree arguments; not theorem inference."""
from collections import deque
from pathlib import Path
import hashlib
import itertools
import json

ROOT = Path(__file__).resolve().parent
def parity(x):
    return x.bit_count() & 1

def pair(x, y, k):
    return sum(((x >> i) & 1) * ((y >> (i ^ 1)) & 1) for i in range(k)) & 1

def q0(x, k):
    return sum(((x >> i) & 1) * ((x >> (i + 1)) & 1) for i in range(0, k, 2)) & 1

def affine_tv(u, a, k):
    return tuple(x ^ (u if pair(u, x, k) ^ a else 0) for x in range(1 << k))

def orbit(perms, start=0):
    reached, queue = {start}, deque([start])
    while queue:
        x = queue.popleft()
        for p in perms:
            y = p[x]
            if y not in reached:
                reached.add(y)
                queue.append(y)
    return reached

def rank(vectors):
    pivots = {}
    for x in vectors:
        while x:
            b = x.bit_length() - 1
            if b in pivots:
                x ^= pivots[b]
            else:
                pivots[b] = x
                break
    return len(pivots)

def refinements(directions, offsets, k):
    return [L for L in range(1 << k)
            if all((q0(u, k) ^ parity(L & u)) == (1 ^ a)
                   for u, a in zip(directions, offsets))]

def bases(k):
    first = [1 << i for i in range(k)]
    second = first[:]
    for i in range(k - 1):
        second[i + 1] ^= second[i]
    third = second[:]
    for i in range(k - 1, 0, -1):
        third[i - 1] ^= third[i]
    return [first, second, third]

controls = []
for k in [2, 4, 6]:
    cases = 0
    max_orbit = 0
    for basis in bases(k):
        assert rank(basis) == k
        for offsets in itertools.product([0, 1], repeat=k):
            perms = [affine_tv(u, a, k) for u, a in zip(basis, offsets)]
            sols = refinements(basis, offsets, k)
            assert len(sols) == 1
            L = sols[0]
            values = [q0(x, k) ^ parity(L & x) for x in range(1 << k)]
            assert set(values) == {0, 1}
            for p in perms:
                assert all(p[p[x]] == x for x in range(1 << k))
                assert sum(p[x] == x for x in range(1 << k)) == (1 << (k - 1))
                assert all(values[x] == values[p[x]] for x in range(1 << k))
            size = len(orbit(perms))
            assert size < (1 << k)
            max_orbit = max(max_orbit, size)
            cases += 1
    controls.append({"dimension": k, "basis_offset_cases": cases,
                     "largest_zero_orbit": max_orbit,
                     "every_case_has_unique_preserved_refinement": True})

# The extra-meridian lower bound can be attained in dimension two.
k = 2
directions, offsets = [1, 2, 1], [0, 0, 1]
perms = [affine_tv(u, a, k) for u, a in zip(directions, offsets)]
assert len(orbit(perms)) == 4
assert refinements(directions, offsets, k) == []

# Inconsistency of the refinement equations is necessary, not sufficient.
k = 4
directions, offsets = [1, 2, 4, 8, 1], [0, 0, 0, 0, 1]
perms = [affine_tv(u, a, k) for u, a in zip(directions, offsets)]
assert refinements(directions, offsets, k) == []
assert len(orbit(perms)) == 4

# Every inconsistent example has an origin-independent violating dependency.
dependency_checks = 0
for k, directions, offsets in [(2, [1, 2, 1], [0, 0, 1]),
                              (4, [1, 2, 4, 8, 1], [0, 0, 0, 0, 1])]:
    violating = []
    for mask in range(1, 1 << len(directions)):
        inds = [i for i in range(len(directions)) if (mask >> i) & 1]
        direction_sum = 0
        for i in inds:
            direction_sum ^= directions[i]
        if direction_sum:
            continue
        value = sum(1 ^ offsets[i] for i in inds) & 1
        for a, i in enumerate(inds):
            for j in inds[a + 1:]:
                value ^= pair(directions[i], directions[j], k)
        if value:
            violating.append(mask)
            for z in range(1 << k):
                shifted = [a ^ pair(u, z, k) for u, a in zip(directions, offsets)]
                shifted_value = sum(1 ^ shifted[i] for i in inds) & 1
                for a, i in enumerate(inds):
                    for j in inds[a + 1:]:
                        shifted_value ^= pair(directions[i], directions[j], k)
                assert shifted_value == 1
                dependency_checks += 1
    assert violating

def cycle(n, entries):
    p = list(range(n))
    for x, y in zip(entries, entries[1:] + entries[:1]):
        p[x] = y
    return tuple(p)

def cycles(p):
    seen, total = set(), 0
    for x in range(len(p)):
        if x not in seen:
            total += 1
            y = x
            while y not in seen:
                seen.add(y)
                y = p[y]
    return total

def product(perms):
    result = tuple(range(len(perms[0])))
    for p in perms:
        result = tuple(p[result[x]] for x in range(len(result)))
    return result

tuple_controls = []
for name, n, perms, transitive in [
    ("four_transpositions_six_sheets", 6, [cycle(6,[i,i+1]) for i in range(4)], False),
    ("five_transpositions_six_sheets", 6, [cycle(6,[i,i+1]) for i in range(5)], True),
    ("three_threecycles_six_sheets", 6,
     [cycle(6,[0,1,2]),cycle(6,[2,3,4]),cycle(6,[4,5,0])], True),
]:
    size = len(orbit(perms))
    assert (size == n) == transitive
    index = sum(n - cycles(p) for p in perms)
    infinity_cycles = cycles(product(perms))
    if transitive:
        twice_genus = index - n - infinity_cycles + 2
        assert twice_genus >= 0 and twice_genus % 2 == 0
    else:
        twice_genus = None
    tuple_controls.append({"name": name, "zero_orbit": size,
                           "index_sum": index, "infinity_cycles": infinity_cycles,
                           "twice_genus_if_transitive": twice_genus})

proof = ROOT / "BREAKTHROUGH_SEARCH_TOPOLOGICAL_COUPLING_2026_11_16.tex"
receipt = {
    "status": "all_exact_controls_pass",
    "scope": "Finite independent controls only; all-degree written proofs remain the evidence.",
    "basis_offset_controls": controls,
    "total_basis_offset_cases": sum(c["basis_offset_cases"] for c in controls),
    "transitive_extra_meridian_control": {"dimension": 2, "meridians": 3, "orbit": 4},
    "refinement_inconsistency_not_sufficient_control": {"dimension": 4, "meridians": 5, "orbit": 4},
    "origin_shift_dependency_checks": dependency_checks,
    "permutation_tuple_controls": tuple_controls,
    "proof_sha256": hashlib.sha256(proof.read_bytes()).hexdigest(),
    "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "no_agent_messages": True,
    "formal_verification": False,
    "external_review": False
}
out = ROOT / "topological-coupling-controls-2026-11-16.json"
out.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
print(json.dumps(receipt, indent=2))
