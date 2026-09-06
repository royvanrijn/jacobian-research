#!/usr/bin/env python3
"""Certified local connecting classes for the u=2 gluing isogenies."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import retrospective as r
import affine_selmer as af
import local_collision as lc
import selmer_comparison as sc

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE/"JACOBIAN_LOCAL_PROTOCOL.json"
OUTPUT = r.OUT/"rank_jump_u2_isogeny_local_conditions_v1.json"
WORK = r.ROOT/"artifacts/local/rank-jump-u2-isogeny-local-v1"


def bindings():
    paths = (PROTOCOL, Path(__file__), af.INPUT, sc.OUTPUT, af.LOCAL)
    return {str(p.relative_to(r.ROOT)): r.digest(p.read_bytes()) for p in paths}


def inputs():
    u = r.read(PROTOCOL)["parameter_u"]
    case = next(c for c in r.read(af.INPUT)["cases"] if c["u"] == u)
    comparison = next(c for c in r.read(sc.OUTPUT)["cases"] if c["u"] == u)
    return case, comparison


def root_class(nf, th, A, B, u, p, approximation, chars):
    from sage.all import QQ, pari
    a = QQ(approximation)
    value = a**3+A*a+B
    derivative = 3*a*a+A
    h = int(value.valuation(p))
    d = int(derivative.valuation(p))
    assert h > 2*d
    accuracy = h-d
    # The numerator has integral coefficients in the root and theta.
    # Changing the root by O(p^accuracy) changes the quotient by
    # O(p^(accuracy-2*d)) in every component.
    e = (th*th+pari(a)*th+pari(A+a*a))/pari(derivative)
    beta = 1+pari(u)*(pari(a)+th)+e*pari(-u*a+u*u*(A+a*a))
    valuations = [int(pari.idealval(nf, beta, P)) for P in chars.primes]
    margins = [int(P[2])*(accuracy-2*d)-v-2*int(P[2])*(1 if p == 2 else 0)
               for P, v in zip(chars.primes, valuations)]
    assert all(m > 0 for m in margins)
    coefficients = [str(QQ(pari.polcoef(pari.lift(beta), i))) for i in range(3)]
    return {"root_approximation": str(a), "v_f": h, "v_derivative": d,
            "root_accuracy_lower_bound": accuracy,
            "class_coordinates": coefficients, "class_ideal_valuations": valuations,
            "square_stability_margins": margins,
            "signature": list(chars.signature(beta))}


def summarize(place, rows, local, h, real=False):
    signatures = list(map(r.pack, [x["signature"] for x in rows]))
    D = r.basis(local["sum_basis"])
    assert all(r.reduce(v, D) == 0 for v in signatures)
    if len(rows) == 3:
        assert signatures[0] ^ signatures[1] ^ signatures[2] == 0
    t = r.rank(signatures)
    ell = len(local["L0_basis"])
    Jtorsion = 2*h-t
    correction = -2 if real else 2 if place == 2 else 0
    middle = 2*ell-t
    assert middle == Jtorsion+correction
    assert middle == len(local["sum_basis"])-t+len(local["intersection_basis"])
    return {"place": place, "root_classes": rows,
            "rational_elliptic_two_torsion_dimension": h,
            "connecting_rank": t, "rational_Jacobian_two_torsion_dimension": Jtorsion,
            "left_local_dimension": len(local["sum_basis"]),
            "right_local_dimension": len(local["intersection_basis"]),
            "Jacobian_local_Kummer_dimension": middle}


def arithmetic(verify=False):
    from sage.all import QQ, AA, PolynomialRing, pari
    from sage.version import version
    sys.path.insert(0, str(r.ROOT/"elliptic-curves/cas"))
    from research_runtime.local_kummer import LocalSquareclasses
    case, comparison = inputs()
    A, B = map(QQ, case["anchor_model"][3:])
    u = case["u"]
    R = PolynomialRing(QQ, "z")
    z = R.gen()
    f = z**3+A*z+B
    primes = [x["place"] for x in case["local"] if x["place"] != "infinity"]
    nf = pari.nfinit([pari(f), primes])
    th = pari.Mod(pari(z), pari(f))
    retained = r.read(OUTPUT) if verify else None
    result = []
    WORK.mkdir(parents=True, exist_ok=True)
    for entry in case["local"]:
        p = entry["place"]
        local = next(x for x in comparison["local"] if x["place"] == p)
        previous = next(x for x in retained["places"] if x["place"] == p) if verify else None
        if p == "infinity":
            roots = sorted(f.roots(AA, multiplicities=False))
            assert len(roots) == 3
            signs = [int(1-u*a < 0) for a in roots]
            rows = []
            for i in range(3):
                signature = [sum(signs[k] for k in range(3) if k != i) % 2
                             if j == i else signs[3-i-j] for j in range(3)]
                rows.append({"root_index": i, "gamma_signs": signs, "signature": signature})
            row = summarize(p, rows, local, 2, real=True)
        else:
            chars = LocalSquareclasses(nf, p)
            h = len(chars.primes)-1
            expected_roots = (0, 1, 3)[h]
            precision = min(128, max(16, 4*int(f.discriminant().valuation(p))+16))
            approximations = ([x["root_approximation"] for x in previous["root_classes"]]
                              if verify else [str(pari.lift(x)) for x in pari.polrootspadic(pari(f), p, precision)])
            assert len(approximations) == expected_roots
            rows = [root_class(nf, th, A, B, u, p, a, chars) for a in approximations]
            for i, a in enumerate(rows):
                for b in rows[i+1:]:
                    assert (QQ(a["root_approximation"])-QQ(b["root_approximation"])).valuation(p) < min(
                        a["root_accuracy_lower_bound"], b["root_accuracy_lower_bound"])
            row = summarize(p, rows, local, h)
            row["requested_root_precision"] = precision
        row["bindings"] = bindings()
        if verify:
            assert row == previous
        else:
            checkpoint = WORK/f"place-{p}.json"
            if checkpoint.exists():
                assert r.read(checkpoint) == row
            else:
                r.write_new(checkpoint, row)
        result.append(row)
        print("PASS place" if verify else "checkpoint", p, "connecting", row["connecting_rank"],
              "J local dimension", row["Jacobian_local_Kummer_dimension"], flush=True)
    out = {"schema": "rank-jump.u2-isogeny-local.v1", "bindings": bindings(),
           "software": {"sage": version, "pari": str(pari("version()"))},
           "u": u, "places": result,
           "scope": "Exact local connecting classes and dimensions; no complete middle local basis or independent CT entry has been computed."}
    if not verify:
        r.write_new(OUTPUT, out)


def check():
    out = r.read(OUTPUT)
    assert out["bindings"] == bindings()
    _, comparison = inputs()
    for row in out["places"]:
        local = next(x for x in comparison["local"] if x["place"] == row["place"])
        rebuilt = summarize(row["place"], row["root_classes"], local,
                            row["rational_elliptic_two_torsion_dimension"],
                            real=row["place"] == "infinity")
        assert all(row[k] == v for k, v in rebuilt.items())
        for root in row["root_classes"]:
            if row["place"] != "infinity":
                assert root["v_f"] > 2*root["v_derivative"]
                assert all(m > 0 for m in root["square_stability_margins"])
    print("PASS local connecting ranks and exact-sequence dimensions")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("capture", "worker", "check", "verify"))
    mode = parser.parse_args().mode
    if mode == "capture":
        WORK.mkdir(parents=True, exist_ok=True)
        p = subprocess.run(["sage", "-python", str(Path(__file__).resolve()), "worker"],
                           cwd=r.ROOT, capture_output=True, text=True,
                           timeout=r.read(PROTOCOL)["limits"]["worker_seconds"])
        with (WORK/"capture.log").open("x") as log:
            log.write(p.stdout+p.stderr)
        print(p.stdout+p.stderr)
        assert p.returncode == 0
    elif mode == "check":
        check()
    else:
        arithmetic(verify=mode == "verify")
