#!/usr/bin/env sage -python
"""Bounded transfer of retained chord/Mestre/local-sieve methods to Hill's X1152.

freeze selects only generic equation/basis words, before any experiment.
run has a 120 CPU-second ceiling and checkpoints each constructed pencil.
No known exceptional fibre, target points, or rank labels select any word.
"""
import argparse
from hashlib import sha256
import importlib.util
import itertools
import json
from pathlib import Path
import resource
import time

from sage.all import Conic, EllipticCurve, GF, PolynomialRing, QQ, ZZ, gcd, lcm, matrix, prime_range, vector

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elkies-k3/data/external/hill2026/data/parent-rank17.json"
CHORD = ROOT / "elkies-k3/scripts/construct_elkies_2026_bisections.sage"


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    # This producer never overwrites a frozen input or existing certificate.
    with path.open("x") as stream:
        json.dump(data, stream, indent=2, sort_keys=True)
        stream.write("\n")


def coeffs(poly):
    return [str(x) for x in poly.list()]


def load_chord():
    from importlib.machinery import SourceFileLoader
    loader = SourceFileLoader("hill1152_retained_chord", str(CHORD))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    return mod


def build(data):
    R = PolynomialRing(QQ, "U")
    K = R.fraction_field()
    A = R(data["jacobian_model"]["a4_low_to_high"])
    B = R(data["jacobian_model"]["a6_low_to_high"])
    E = EllipticCurve(K, [A, B])
    P = [E(K(R(r["x_numerator"]))/R(r["x_denominator"]),
           K(R(r["y_numerator"]))/R(r["y_denominator"]))
         for r in data["neighbor_sections"]]
    G = matrix(QQ, data["exact_generic_section_basis"]["height_matrix"])
    assert len(P) == 17 and G.det() == 1152 and G.is_positive_definite()
    return R, K, A, B, E, P, G


def select_words(G, norm):
    # One deliberately selected representative per parity; no assertion that
    # the resulting finite bank exhausts the lattice or all bisections.
    choices = {}
    for i in range(17):
        w = [0]*17
        w[i] = 1
        if G[i, i] == norm:
            choices[tuple(w)] = w
    for i, j in itertools.combinations(range(17), 2):
        for sign in (1, -1):
            if G[i, i]+G[j, j]+2*sign*G[i, j] == norm:
                w = [0]*17
                w[i], w[j] = 1, sign
                choices.setdefault(tuple(x % 2 for x in w), w)
    return sorted(choices.values())


def freeze(out):
    data = json.loads(SOURCE.read_text())
    _, _, _, _, _, _, G = build(data)
    write(out/"input.json", {
        "schema": "hill1152-method-transfer-input-v1",
        "source": str(SOURCE.relative_to(ROOT)),
        "source_sha256": digest(SOURCE),
        "chord_sha256": digest(CHORD),
        "producer_sha256": digest(Path(__file__)),
        "selection": "All lexicographically sorted, parity-distinct representatives of norm 8 or 10 supported on at most two supplied generic basis vectors; sign + preferred. No exceptional fibres.",
        "words8": select_words(G, 8),
        "words10": select_words(G, 10),
        "screen_primes": [101, 103, 107, 109, 113, 127, 131],
        "local_prime_depths": [[2, 3], [3, 2], [5, 2], [7, 2], [11, 2]],
        "mestre_prime_bound_exclusive": 2000,
        "cpu_seconds": 120, "address_space_bytes": 4294967296,
    })


def compose(f, a):
    return f.numerator()(a)/f.denominator()(a)


def mestre_gate(A, B, packet):
    R = A.parent()
    D = -5*A*B*(9261*B**2+400*A**3)
    delta = 4*A**3+27*B**2
    x1, x2 = -21*B/(5*A), -21*B/(20*A)
    c1, c2 = 1/(25*A**2), 1/(200*A**2)
    assert D*c1*c1 == x1**3+A*x1+B
    assert D*c2*c2 == x2**3+A*x2+B
    xsum = D*((c2-c1)/(x2-x1))**2-x1-x2
    assert xsum == -QQ(7)/15*B/A-QQ(20)/81*A**2/B
    fields = {"A": {}, "B": {}}
    squarefree_prime = None
    for p in prime_range(5, packet["mestre_prime_bound_exclusive"]):
        F = PolynomialRing(GF(p), "t")
        a, b, d, e = [F(f) for f in (A, B, D, delta)]
        for name, f, src in [("A", a, A), ("B", b, B)]:
            if f.degree() != src.degree():
                continue
            if "irreducible_prime" not in fields[name] and f.is_irreducible():
                fields[name]["irreducible_prime"] = int(p)
            congruence = p % 12 == 11 if name == "A" else p % 4 == 3
            if congruence and "local_place" not in fields[name]:
                roots = sorted(int(z) for z in f.roots(multiplicities=False) if f.derivative()(z))
                if roots:
                    fields[name]["local_place"] = {"prime": int(p), "root": roots[0]}
        if (squarefree_prime is None and d.degree() == 44 and e.degree() == 24
                and d.is_squarefree() and e.is_squarefree() and d.gcd(e).degree() == 0):
            squarefree_prime = int(p)
        if squarefree_prime and all(len(v) == 2 for v in fields.values()):
            break
    passed = squarefree_prime is not None and all(len(v) == 2 for v in fields.values())
    return {"status": "PASS" if passed else "UNKNOWN", "fields": fields,
            "smooth_disjoint_prime": squarefree_prime, "cover_degree": int(D.degree()),
            "fixed_u2_genus": 21 if squarefree_prime else None,
            "all_rational_auxiliary_genus_lower_bound": 9 if passed else None,
            "scope": "Only Mestre's retained identity on this exact equation. The valuation proof in the canonical note is required; not a universal cover obstruction."}


def norm_projective(v):
    p = next(x for x in v if x)
    return tuple(int(x/p) for x in v)


def local_binary_values(q, p, depth):
    # Clear rational denominators by a square, then remove even common
    # p-adic valuation. This retains the actual constant squareclass.
    den = ZZ(lcm([c.denominator() for c in q.list()]))
    a = [ZZ(c*den**2) for c in q.list()]
    v = min(c.valuation(p) for c in a if c)
    a = [c//ZZ(p)**(2*(v//2)) for c in a]
    a += [ZZ(0)]*(3-len(a))
    modulus = p**depth
    squares = {n*n % modulus for n in range(modulus)}
    # P^1(Z/p^depth): z=1, or t=1 and z divisible by p.
    addresses = [(t, 1) for t in range(modulus)]
    addresses += [(1, p*z) for z in range(modulus//p)]
    return {i for i, (t, z) in enumerate(addresses)
            if int((a[2]*t*t+a[1]*t*z+a[0]*z*z) % modulus) in squares}


def run(out):
    packet = json.loads((out/"input.json").read_text())
    assert packet["source_sha256"] == digest(SOURCE)
    assert packet["chord_sha256"] == digest(CHORD)
    assert packet["producer_sha256"] == digest(Path(__file__))
    resource.setrlimit(resource.RLIMIT_CPU, (packet["cpu_seconds"], packet["cpu_seconds"]+5))
    resource.setrlimit(resource.RLIMIT_AS, (packet["address_space_bytes"], packet["address_space_bytes"]))
    started = time.process_time()
    data = json.loads(SOURCE.read_text())
    R, K, A, B, E, basis, G = build(data)
    assert packet["words8"] == select_words(G, 8)
    assert packet["words10"] == select_words(G, 10)
    chord = load_chord()
    traces = [sum((int(n)*P for n, P in zip(w, basis)), E(0))
              for w in packet["words8"]+packet["words10"]]
    # A common Mobius chart moves every trace pole away from infinity.
    c = next(QQ(c) for c in [0, 1, -1, 2, -2, 3]
             if all(P[0].denominator()(c) != 0 for P in traces))
    u = R.gen()
    old_t = c+1/u
    a, b = R(u**8*A(old_t)), R(u**12*B(old_t))
    delta = 4*a**3+27*b**2
    assert delta.degree() == 24 and delta.is_squarefree()
    write(out/"mestre.json", mestre_gate(A, B, packet))
    pencils, covers = [], []
    for index, (word, trace) in enumerate(zip(packet["words8"]+packet["words10"], traces)):
        norm = int(vector(QQ, word)*G*vector(QQ, word))
        x, y = u**4*compose(trace[0], old_t), u**6*compose(trace[1], old_t)
        frame = chord.trace_chord_frame(x, y, R)
        h, Nx, Ny, M = [frame[k] for k in ("h", "Nx", "Ny", "M0")]
        assert h.degree() == (norm-4)//2
        common = {"word": word, "h": coeffs(h), "Nx": coeffs(Nx),
                  "Ny": coeffs(Ny), "M0": coeffs(M)}
        if norm == 8:
            numerators = [M**4-6*M**2*Nx-8*M*Ny-3*Nx**2-4*a*h**4,
                (4*M**3-12*M*Nx-8*Ny)*h**2,
                6*(M**2-Nx)*h**4, 4*M*h**6, h**8]
            qs = []
            for numerator in numerators:
                q, rem = numerator.quo_rem(h**6)
                assert not rem and q.degree() <= 4
                qs.append(q)
            branch = matrix(QQ, 5, 5, lambda i, j: qs[j][i])
            assert branch.det()
            row = {**common, "index": len(pencils),
                   "branch_matrix": [[str(v) for v in line] for line in branch.rows()]}
            pencils.append(row)
            write(out/"pencils"/f"{len(pencils)-1:03d}.json", row)
        else:
            raw_q, remainder = (M**4-6*M**2*Nx-8*M*Ny-3*Nx**2-4*a*h**4).quo_rem(h**6)
            assert not remainder
            if raw_q.degree() not in (1, 2) or not raw_q.is_squarefree():
                row = {**common, "label": f"generic-norm10-{len(covers):02d}",
                       "status": "DEGENERATE_RESIDUAL", "q": coeffs(raw_q),
                       "rational_conic": False, "conic_point": None}
                covers.append(row)
                write(out/"covers"/f"{len(covers)-1:03d}.json", row)
                continue
            d = chord.chord_data_from_slope_numerator(h, Nx, Ny, M, a, b, delta, R, K, int(raw_q.degree()))
            q = d["q"]
            x0, x1, y0, y1 = [d[k] for k in ("x0", "x1", "y0", "y1")]
            xr = y0**2/(x1**2*q)-2*x0
            yr = -y1+y0*(x0-xr)/(x1*q)
            assert yr and q*yr**2 == xr**3+a*xr+b
            S = PolynomialRing(QQ, ["T", "Z", "S"])
            T, Z, ss = S.gens()
            conic = Conic(QQ, ss**2-q[2]*T**2-q[1]*T*Z-q[0]*Z**2)
            soluble, point = conic.has_rational_point(point=True)
            witness = [str(z) for z in point] if soluble else None
            if soluble:
                t0, z0, s0 = map(QQ, point)
                assert s0*s0 == q[2]*t0*t0+q[1]*t0*z0+q[0]*z0*z0
            row = {**common, "label": f"generic-norm10-{len(covers):02d}",
                   "q": coeffs(q), "branch_at_infinity": bool(q.degree() == 1), "smooth_branch": bool(d["branch_fibres_smooth"]),
                   "rational_conic": bool(soluble), "conic_point": witness,
                   "nonzero_anti_invariant": True,
                   "x0": coeffs(x0), "x1": coeffs(x1),
                   "y0": coeffs(y0), "y1": coeffs(y1)}
            covers.append(row)
            write(out/"covers"/f"{len(covers)-1:03d}.json", row)
        print(json.dumps({"trace": index+1, "norm": norm, "cpu_seconds": round(time.process_time()-started, 3)}), flush=True)
    write(out/"pencils.json", pencils)
    matrices = [matrix(QQ, p["branch_matrix"]) for p in pencils]
    remaining = set(itertools.combinations(range(len(pencils)), 2))
    screens = []
    for p in packet["screen_primes"]:
        F = GF(p)
        images = {}
        for i, mat in enumerate(matrices):
            try:
                mat = matrix(F, mat)
            except (ValueError, ZeroDivisionError, TypeError):
                continue
            if not mat.det():
                continue
            images[i] = {norm_projective(mat*vector(F, [1, z, z*z, z**3, z**4]))
                         for z in F}
            images[i].add(norm_projective(mat.column(4)))
            assert len(images[i]) == p+1
        rejected = sorted((i, j) for i, j in remaining
                          if i in images and j in images and images[i].isdisjoint(images[j]))
        remaining.difference_update(rejected)
        screens.append({"prime": p, "usable": sorted(images), "rejected_pairs": [list(v) for v in rejected],
                        "remaining_pairs": len(remaining)})
        write(out/f"prime-{p}.json", screens[-1])
    # Known covers enter only as separately labelled comparison/local-sieve
    # controls, after generic selection and construction are sealed.
    for i in [4339, 36773]:
        d = json.loads((SOURCE.parent/f"cover-{i}.json").read_text())
        original = R(d["quadratic_basechange"]["branch_coefficients"])
        q = R(u**2*original(old_t))
        covers.append({"label": f"published-cover-{i}", "q": coeffs(q)})
    candidates = set(itertools.combinations(range(len(covers)), 2))
    distinct = {(i, j) for i, j in candidates
                if covers[i].get("status") != "DEGENERATE_RESIDUAL" and covers[j].get("status") != "DEGENERATE_RESIDUAL"
                and R(covers[i]["q"]).monic() != R(covers[j]["q"]).monic()
                and R(covers[i]["q"]).gcd(R(covers[j]["q"])).degree() == 0}
    local_sieves = []
    for p, depth in packet["local_prime_depths"]:
        images = [local_binary_values(R(row["q"]), p, depth) for row in covers]
        rejected = sorted((i, j) for i, j in distinct if images[i].isdisjoint(images[j]))
        distinct.difference_update(rejected)
        local_sieves.append({"prime": p, "depth": depth, "rejected_pairs": [list(v) for v in rejected],
                             "remaining_pairs": len(distinct)})
    result = {"schema": "hill1152-method-transfer-result-v1",
        "input_sha256": digest(out/"input.json"),
        "common_chart": {"old_t": str(old_t), "c": str(c), "A": coeffs(a), "B": coeffs(b)},
        "pencil_count": len(pencils), "pencil_pairs": len(pencils)*(len(pencils)-1)//2,
        "projective_screens": screens, "unresolved_projective_pairs": [list(v) for v in sorted(remaining)],
        "quadratic_covers": covers, "local_sieves": local_sieves,
        "locally_surviving_disjoint_pairs": [list(v) for v in sorted(distinct)],
        "cpu_seconds": round(time.process_time()-started, 6),
        "boundary": "Finite generic word bank. Projective exclusion is necessary-condition arithmetic, not a universal no-go. Conic solubility is certified; local pair survival is not a rational point or a rank statement. Same-cover two-gain objective stays open."}
    write(out/"result.json", result)
    print(json.dumps({"pencils": len(pencils), "unresolved_pencil_pairs": len(remaining),
                      "constructed_quadratic_covers": len(covers)-2,
                      "rational_conics": sum(r.get("rational_conic", False) for r in covers),
                      "locally_surviving_pairs": len(distinct), "cpu_seconds": result["cpu_seconds"]}), flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=["freeze", "run"])
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    (freeze if args.mode == "freeze" else run)(args.output)
