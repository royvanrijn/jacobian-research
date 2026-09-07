#!/usr/bin/env python3
"""Exact finite model of J[2] and the difference of two E[4] extensions."""
import argparse
from itertools import permutations, product
from pathlib import Path
import retrospective as r
import local_collision as lc
from cubic_bridge import Cubic

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE/"TORSION_DIFFERENCE_PROTOCOL.json"
FOUR = r.OUT/"rank_jump_four_division_separation_v1.json"
OUTPUT = r.OUT/"rank_jump_torsion_difference_v1.json"


def act(M, v):
    return lc.lift(v, M)


def compose(A, B):
    return tuple(act(A, b) for b in B)


def add(A, B):
    return tuple(a ^ b for a, b in zip(A, B))


def canonical_subset(s):
    assert s.bit_count() % 2 == 0
    return min(s, s ^ 63)


def inject(v):
    return canonical_subset((3 if v & 1 else 0) ^ (12 if v & 2 else 0))


def section(v):
    t = next(t for t in range(8) if t.bit_count() % 2 == 0
             and lc.lift(t, [1, 2, 3]) == v)
    return canonical_subset(sum(((t >> i) & 1) << (2*i) for i in range(3)))


def encode_state(k, v):
    return canonical_subset(inject(k) ^ section(v))


def permute_subset(s, perm, signs):
    out = 0
    for i in range(3):
        for bit in range(2):
            if s >> (2*i + bit) & 1:
                j = perm[i]
                out ^= 1 << (2*j + (bit ^ ((signs >> j) & 1)))
    return canonical_subset(out)


def weil(a, b):
    return ((a & 1)*((b >> 1) & 1) + ((a >> 1) & 1)*(b & 1)) % 2


N = tuple(tuple(t if weil(t, v) else 0 for v in (1, 2)) for t in (1, 2, 3))


def sign_matrix(signs):
    return tuple(lc.lift(signs, [M[j] for M in N]) for j in range(2))


def lifted_act(M, x):
    # Any 0/1 lift of a GL2(F2) matrix is invertible modulo four.
    return tuple(sum(((M[j] >> i) & 1)*x[j] for j in range(2)) % 4 for i in range(2))


def baer_coordinates(x, y):
    assert all((a-b) % 2 == 0 for a, b in zip(x, y))
    k = sum((((y[i]-x[i])//2) % 2) << i for i in range(2))
    v = sum((x[i] % 2) << i for i in range(2))
    return k | (v << 2)


def finite_model():
    states = {encode_state(k, v): k | (v << 2) for k in range(4) for v in range(4)}
    assert len(states) == 16
    matrices, records = [], []
    checked_baer = 0
    for perm in permutations(range(3)):
        g = (perm[0]+1, perm[1]+1)
        assert act(g, 3) == perm[2]+1
        for i in range(3):
            assert compose(g, N[i]) == compose(N[perm[i]], g)
        for signs in range(8):
            c = sign_matrix(signs)
            cg = compose(c, g)
            M = (g[0], g[1], cg[0] | (g[0] << 2), cg[1] | (g[1] << 2))
            for subset, state in states.items():
                assert states[permute_subset(subset, perm, signs)] == act(M, state)
            for x in product(range(4), repeat=2):
                for k in range(4):
                    y = tuple((x[i]+2*((k >> i) & 1)) % 4 for i in range(2))
                    Ax, Ay = lifted_act(g, x), lifted_act(g, y)
                    cay = act(c, sum((Ay[i] % 2) << i for i in range(2)))
                    By = tuple((Ay[i]+2*((cay >> i) & 1)) % 4 for i in range(2))
                    assert baer_coordinates(Ax, By) == act(M, baer_coordinates(x, y))
                    checked_baer += 1
            parity = signs.bit_count() % 2
            eta = signs ^ (7 if parity else 0)
            assert eta.bit_count() % 2 == 0
            assert c == add((1, 2) if parity else (0, 0), sign_matrix(eta))
            matrices.append(M)
            records.append({"permutation": list(perm), "target_sign_mask": signs,
                            "extension_cocycle_matrix_columns": c, "J2_matrix_columns": M,
                            "norm_sign": parity, "eta_sign_mask": eta})
    assert len(set(matrices)) == 48
    assert all(compose(A, B) in matrices for A in matrices for B in matrices)
    # Commutant as the kernel of exact linear equations in sixteen entries.
    equations = []
    elementary = [tuple((1 << i) if j == k else 0 for j in range(4))
                  for k in range(4) for i in range(4)]
    for M in matrices:
        images = [add(compose(X, M), compose(M, X)) for X in elementary]
        for j in range(4):
            for i in range(4):
                equations.append(r.pack([(X[j] >> i) & 1 for X in images]))
    centralizer = lc.orthogonal(equations, 16)
    ident = (1, 2, 4, 8)
    nilpotent = (0, 0, 1, 2)
    pack_matrix = lambda M: sum(v << (4*j) for j, v in enumerate(M))
    assert centralizer == lc.canonical([pack_matrix(ident), pack_matrix(nilpotent)])
    assert compose(nilpotent, nilpotent) == (0, 0, 0, 0)
    idempotents = []
    for mask in range(4):
        packed = lc.lift(mask, centralizer)
        M = tuple((packed >> (4*j)) & 15 for j in range(4))
        if compose(M, M) == M:
            idempotents.append(packed)
    assert sorted(idempotents) == [0, pack_matrix(ident)]
    return {"transvection_N_columns": N, "actions": records,
            "root_action_checks": 48*16, "Baer_action_checks": checked_baer,
            "commuting_algebra_basis_packed": centralizer,
            "commuting_algebra": "F2[epsilon]/(epsilon^2)",
            "commuting_idempotents_packed": sorted(idempotents),
            "module_status": "INDECOMPOSABLE_NONSPLIT_SELF_EXTENSION_OF_V"}


def build(check=False):
    finite = finite_model()
    inp = r.read(lc.INPUT)
    four = r.read(FOUR)
    assert four["input_sha256"] == r.digest(lc.INPUT.read_bytes())
    A, B = map(r.F, inp["anchor"]["short_model_ainvariants"][3:])
    K = Cubic(A, B)
    rows = []
    for u in r.read(PROTOCOL)["parameters"]:
        witness = next(x for x in four["rows"] if int(x["u"]) == u)
        D = 1+A*u*u+B*u**3
        gamma = K.sub(K.one, K.scale(K.theta, u))
        eta = K.scale(gamma, D)
        assert K.norm(gamma) == D and K.norm(eta) == D**4
        assert witness["relative_degree_over_anchor_four_division_field"] == 8
        rows.append({"u": u, "gamma_coordinates": list(map(str, gamma)),
                     "eta_coordinates": list(map(str, eta)), "D": str(D),
                     "simple_collision_prime": witness["simple_good_collision_prime"],
                     "J2_splitting_field_degree_over_Q": 48,
                     "J2_galois_image": "(C2)^3 semidirect S3",
                     "J2_module_indecomposable": True,
                     "new_CT_value_computed": False})
    paths = (PROTOCOL, Path(__file__), FOUR, lc.INPUT, HERE/"retrospective.py",
             HERE/"local_collision.py", HERE/"cubic_bridge.py")
    out = {"schema": "rank-jump.torsion-difference.v1",
           "bindings": {str(p.relative_to(r.ROOT)): r.digest(p.read_bytes()) for p in paths},
           "finite_model": finite, "cases": rows,
           "scope": "Exact underlying torsion-extension model. CT comparison requires the decorated local Baer-sum structure; the module alone is not a solubility or rank predictor."}
    # Normalize tuples so construction and replay compare the same JSON object.
    import json
    out = json.loads(json.dumps(out))
    if check:
        assert r.read(OUTPUT) == out
        print("PASS torsion difference, sign decomposition and commuting algebra")
    else:
        r.write_new(OUTPUT, out)
    print("root checks", finite["root_action_checks"], "Baer checks", finite["Baer_action_checks"])
    print("commutant", finite["commuting_algebra"], "fibre fields", [48]*len(rows))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("build", "check"))
    build(parser.parse_args().mode == "check")
