#!/usr/bin/env sage-python
"""Prospective, bounded same-cover test on generic R17 genus-one pencils.

No exceptional points or fibre parameters enter the frozen arithmetic packet.
The result concerns the selected complete pencils, not all K3 bisections.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from importlib.machinery import SourceFileLoader
from importlib.util import module_from_spec, spec_from_loader
import itertools
import json
from pathlib import Path
import time

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, matrix, vector

ROOT = Path(__file__).resolve().parents[2]
CHORD = ROOT / "elkies-k3/scripts/construct_elkies_2026_bisections.sage"
DATA = ROOT / "elkies-k3/data/fibrations"
DEFAULT = ROOT / "artifacts/generated-results/elkies-k3-r17-correlated-genus-one-v1"


def load_chord():
    loader = SourceFileLoader("correlated_chord", str(CHORD))
    spec = spec_from_loader(loader.name, loader)
    result = module_from_spec(spec)
    loader.exec_module(result)
    return result


def write_new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def coeffs(poly, size=None):
    return [str(poly[i]) for i in range(size if size is not None else poly.degree() + 1)]


def rows(mat):
    return [[str(x) for x in row] for row in mat.rows()]


def freeze(output):
    chord = load_chord()
    model_path = DATA / "elkies_2026_published_r17_model.json"
    sections_path = DATA / "elkies_2026_published_r17_sections.json"
    id_path = ROOT / "artifacts/generated-results/elkies-2026-published-r17-target.json"
    pin_path = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
    model = json.loads(model_path.read_text())
    data = json.loads(sections_path.read_text())
    R = PolynomialRing(QQ, "t")
    A = R(model["A_coefficients_low_to_high"])
    B = R(model["B_coefficients_low_to_high"])
    basis = chord.reconstruct_basis(R, A, B, data)
    identification = json.loads(id_path.read_text())["pinned_identification"]
    assert identification["gram_identity_orientation"] == "M^T*Gpub*M=Gpinned"
    M = matrix(QQ, identification["basis_change_matrix"])
    pinned = chord.load_matrix(pin_path)
    G = M.inverse().transpose() * pinned * M.inverse()
    assert M.transpose() * G * M == pinned and G.det() == 948
    choices = {}
    for i, j in itertools.combinations(range(17), 2):
        for sign in (1, -1):
            word = [0] * 17
            word[i], word[j] = 1, sign
            w = vector(ZZ, word)
            if w * G * w != 8:
                continue
            parity = tuple(x % 2 for x in word)
            choices.setdefault(parity, word)
    selected = sorted(choices.values())[:64]
    packet = {
        "schema": "r17-correlated-genus-one-input-v1",
        "scope": "At most 64 parity-distinct norm-eight words supported on two published generic basis vectors; lexicographic order; positive second coefficient preferred.",
        "selection_pool_count": len(choices),
        "selected_words": selected,
        "A": coeffs(A), "B": coeffs(B),
        "basis": [{"x": coeffs(x), "y": coeffs(y)} for x, y in basis],
        "gram": rows(G),
        "source_hashes": {str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest()
                          for p in (model_path, sections_path, id_path, pin_path, CHORD)},
        "quarantine": "Only model coefficients, generic section coordinates and the generic height identification are exported. Public specializations, exceptional points and quotient labels are not worker inputs. Field allowlisting is not OS isolation.",
        "limits": {"pencils": 64, "screen_primes": [101, 103, 107],
                   "cpu_seconds": 120, "address_space_gib": 4},
        "acceptance": ["same quadratic extension over the fixed original t-line",
                       "two exact section maps", "independence modulo inherited subgroup",
                       "rational base or genus-one base with certified positive rank"],
    }
    write_new(output / "input.json", packet)
    print(json.dumps({"frozen_pencils": len(selected), "selection_pool_count": len(choices)}), flush=True)


def normalize(v):
    pivot = next(x for x in v if x)
    return tuple(int(x / pivot) for x in v)


def veronese(v):
    return all(v[i] * v[j + 1] == v[i + 1] * v[j]
               for i, j in itertools.combinations(range(4), 2))


def run(output):
    import resource
    packet = json.loads((output / "input.json").read_text())
    resource.setrlimit(resource.RLIMIT_CPU, (120, 125))
    resource.setrlimit(resource.RLIMIT_AS, (4 * 1024**3, 4 * 1024**3))
    started = time.process_time()
    chord = load_chord()
    R = PolynomialRing(QQ, "t")
    t = R.gen()
    K = R.fraction_field()
    A, B = R(packet["A"]), R(packet["B"])
    E = EllipticCurve(K, [A, B])
    basis = [E(K(R(row["x"])), K(R(row["y"]))) for row in packet["basis"]]
    G = matrix(QQ, packet["gram"])
    pencils = []
    for index, word in enumerate(packet["selected_words"]):
        w = vector(ZZ, word)
        assert w * G * w == 8
        trace = sum((int(n) * p for n, p in zip(word, basis)), E(0))
        X, Y = trace[0], trace[1]
        frame = chord.trace_chord_frame(X, Y, R)
        h, Nx, Ny, M = (frame[key] for key in ("h", "Nx", "Ny", "M0"))
        if h.degree() != 2:
            raise ArithmeticError("Selected trace needs an unimplemented infinity chart; do not silently omit it")
        numerators = [M**4 - 6*M**2*Nx - 8*M*Ny - 3*Nx**2 - 4*A*h**4,
                      (4*M**3 - 12*M*Nx - 8*Ny)*h**2,
                      6*(M**2-Nx)*h**4, 4*M*h**6, h**8]
        qs = []
        for numerator in numerators:
            q, rem = numerator.quo_rem(h**6)
            assert not rem and q.degree() <= 4
            qs.append(q)
        branch = matrix(QQ, 5, 5, lambda i, j: qs[j][i])
        assert branch.det()
        record = {"index": index, "word": word,
                  "h": coeffs(h), "Nx": coeffs(Nx), "Ny": coeffs(Ny),
                  "M0": coeffs(M), "branch_matrix": rows(branch)}
        pencils.append(record)
        write_new(output / "pencils" / f"{index:03d}.json", record)
    write_new(output / "pencils.json", pencils)
    print(json.dumps({"constructed": len(pencils), "cpu_seconds": time.process_time()-started}), flush=True)
    matrices = [matrix(QQ, row["branch_matrix"]) for row in pencils]
    remaining = set(itertools.combinations(range(len(pencils)), 2))
    witnesses = []
    for prime in packet["limits"]["screen_primes"]:
        F = GF(prime)
        images = {}
        usable = []
        for i, mat in enumerate(matrices):
            try:
                reduced = matrix(F, mat)
            except (ValueError, ZeroDivisionError, TypeError):
                continue
            if not reduced.det():
                continue
            usable.append(i)
            for value in list(F) + [None]:
                vv = vector(F, [1, value, value**2, value**3, value**4]) if value is not None else vector(F, [0,0,0,0,1])
                images.setdefault(normalize(reduced * vv), set()).add(i)
        meet = set()
        for bucket in images.values():
            meet.update(itertools.combinations(sorted(bucket), 2))
        usable = set(usable)
        rejected = sorted(pair for pair in remaining if set(pair) <= usable and pair not in meet)
        remaining.difference_update(rejected)
        record = {"prime": prime, "usable": sorted(usable), "rejected_pairs": [list(p) for p in rejected],
                  "remaining_pairs": len(remaining)}
        witnesses.append(record)
        write_new(output / f"prime-{prime}.json", record)
        print(json.dumps({"prime": prime, "rejected": len(rejected), "remaining": len(remaining)}), flush=True)
    # Rational-normal quartic membership reduces to a univariate polynomial gcd.
    # No exceptional specialization or bounded parameter search is involved.
    S = PolynomialRing(QQ, "lambda")
    z = S.gen()
    exact = []
    for i, j in sorted(remaining):
        transport = matrices[j].inverse() * matrices[i]
        values = transport * vector(S, [1, z, z*z, z**3, z**4])
        minors = [values[a]*values[b+1]-values[a+1]*values[b]
                  for a, b in itertools.combinations(range(4), 2)]
        divisor = S(0)
        for minor in minors:
            divisor = divisor.gcd(minor)
        if not divisor:
            raise ArithmeticError("Entire branch images coincide: requires a separate pencil-equivalence audit")
        divisor = divisor.monic()
        roots = divisor.roots(QQ, multiplicities=False)
        candidates = []
        for root in roots:
            values_at_root = vector(QQ, [f(root) for f in values])
            assert veronese(values_at_root)
            other = None if not values_at_root[0] else values_at_root[1]/values_at_root[0]
            q = R(list(matrices[i] * vector(QQ, [1, root, root**2, root**3, root**4])))
            squarefree = q.degree() in (3,4) and q.gcd(q.derivative()).degree() == 0
            candidates.append({"lambda_i": str(root), "lambda_j": None if other is None else str(other),
                               "quartic": coeffs(q), "genus_one": squarefree,
                               "projective_scale": str(values_at_root[0]) if other is not None else str(values_at_root[4])})
        infinity = transport.column(4)
        record = {"pair": [i,j], "minor_gcd": coeffs(divisor),
                  "infinity_meets": bool(veronese(infinity)), "rational_candidates": candidates}
        exact.append(record)
        write_new(output / "pairs" / f"{i:03d}-{j:03d}.json", record)
    result = {"schema": "r17-correlated-genus-one-result-v1",
              "input_sha256": sha256((output/"input.json").read_bytes()).hexdigest(),
              "pencils": len(pencils), "pair_count": len(pencils)*(len(pencils)-1)//2,
              "modular_witnesses": witnesses, "exact_pairs": exact,
              "cpu_seconds": time.process_time()-started,
              "positive_result": "UNPROVEN",
              "boundary": "A projective branch match still requires its rational constant squareclass, exact maps, independence, and infinitely many rational base points. This finite bank is not a global nonexistence theorem."}
    write_new(output / "result.json", result)
    print(json.dumps({"completed_pairs": result["pair_count"], "exact_pairs": len(exact),
                      "rational_matches": sum(len(r["rational_candidates"]) for r in exact),
                      "cpu_seconds": result["cpu_seconds"]}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("freeze", "run"))
    parser.add_argument("--output", type=Path, default=DEFAULT)
    args = parser.parse_args()
    (freeze if args.mode == "freeze" else run)(args.output)
