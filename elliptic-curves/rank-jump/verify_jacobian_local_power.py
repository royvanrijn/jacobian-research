#!/usr/bin/env sage -python
"""Independent PARI local-power replay of every connecting-kernel relation."""
from sage.all import QQ, PolynomialRing, pari
import retrospective as r
import jacobian_local_conditions as jlc


def verify():
    data = r.read(jlc.OUTPUT)
    assert data["bindings"] == jlc.bindings()
    case, _ = jlc.inputs()
    A, B = map(QQ, case["anchor_model"][3:])
    R = PolynomialRing(QQ, "z")
    z = R.gen()
    f = z**3+A*z+B
    primes = [p["place"] for p in data["places"] if p["place"] != "infinity"]
    nf = pari.nfinit([pari(f), primes])
    th = pari.Mod(pari(z), pari(f))
    checks = 0
    for local in data["places"]:
        p = local["place"]
        if p == "infinity":
            continue
        roots = local["root_classes"]
        elements = [sum(pari(QQ(c))*th**i for i, c in enumerate(a["class_coordinates"]))
                    for a in roots]
        ideals = list(pari.idealprimedec(nf, p))
        for mask in range(1 << len(roots)):
            value, signature = pari(1), 0
            for i, element in enumerate(elements):
                if mask >> i & 1:
                    value *= element
                    signature ^= r.pack(roots[i]["signature"])
            value *= pari.denominator(pari.nfalgtobasis(nf, value))**2
            square = all(pari.nfislocalpower(nf, P, value, 2) for P in ideals)
            assert square == (signature == 0)
            checks += 1
        print("PASS independent connecting-kernel power test", p, flush=True)
    print("PASS", checks, "local-power comparisons")


if __name__ == "__main__":
    verify()
