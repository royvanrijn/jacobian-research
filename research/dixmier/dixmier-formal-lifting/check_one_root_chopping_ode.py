#!/usr/bin/env python3
"""Exact local one-root chopping ODE/order/face controls; no Weyl pair claim."""
from fractions import Fraction as Q
from pathlib import Path
from math import gcd
import argparse
import hashlib
import json
import platform
import tempfile

HERE = Path(__file__).resolve().parent
# Commuting variables x,y,C. C is kept indeterminate throughout.
ZERO = (0, 0, 0)
ONE = {ZERO: Q(1)}
X = {(1, 0, 0): Q(1)}
Y = {(0, 1, 0): Q(1)}
C = {(0, 0, 1): Q(1)}


def clean(p):
    return {e: c for e, c in p.items() if c}


def add(*polys):
    out = {}
    for p in polys:
        for e, c in p.items():
            out[e] = out.get(e, Q(0)) + c
    return clean(out)


def scale(p, c):
    return clean({e: c * a for e, a in p.items()})


def mul(a, b):
    out = {}
    for e, c in a.items():
        for f, d in b.items():
            ef = tuple(x + y for x, y in zip(e, f))
            out[ef] = out.get(ef, Q(0)) + c * d
    return clean(out)


def power(p, n):
    out = ONE
    for _ in range(n):
        out = mul(out, p)
    return out


def derivative(p, index):
    out = {}
    for e, c in p.items():
        if e[index]:
            ef = list(e)
            ef[index] -= 1
            out[tuple(ef)] = c * e[index]
    return out


def operator(g, H, l, v):
    return add(scale(mul(derivative(H, 1), g), -l),
               scale(mul(H, derivative(g, 1)), v))


def face(p, index):
    degree = max(e[index] for e in p)
    return degree, {tuple(0 if j == index else k for j, k in enumerate(e)): c
                    for e, c in p.items() if e[index] == degree}


def poisson(a, b):
    return add(mul(derivative(a, 0), derivative(b, 1)),
               scale(mul(derivative(a, 1), derivative(b, 0)), -1))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    u = add(ONE, Y)
    g = power(u, 4)
    H = add(scale(power(u, 5), Q(1, 2)), mul(C, power(u, 6)))
    H1 = add(scale(power(u, 20), Q(1, 16)),
             scale(mul(C, power(u, 21)), Q(1, 2)),
             scale(mul(power(C, 2), power(u, 22)), Q(3, 2)),
             scale(mul(power(C, 3), power(u, 23)), 2))
    assert operator(g, H, 2, 3) == power(g, 2)
    assert operator(g, power(u, 6), 2, 3) == {}
    assert H1 == add(power(H, 4), scale(mul(power(C, 4), power(g, 6)), -1))
    assert operator(g, H1, 2, 12) == scale(mul(power(H, 3), power(g, 2)), 4)
    assert face(H, 1) == (6, C)
    assert face(H1, 1) == (23, scale(power(C, 3), 2))
    # Complete polynomial kernel for this control follows diagonally:
    # L_v(u^j)=(d*v-l*j)u^(j+d-1). Check the only zero through degree80.
    assert [j for j in range(81) if 4 * 3 - 2 * j == 0] == [6]
    for j in [0, 1, 5, 6, 7, 12]:
        expected = scale(power(u, j + 3), 12 - 2 * j)
        assert operator(g, power(u, j), 2, 3) == expected

    # Center u=1+y, then use X temporarily as an indeterminate common
    # translation parameter. Reconstruct both translated polynomials.
    def substitute_y(poly, target):
        return add(*(scale(mul(mul(power(X, e[0]), power(target, e[1])),
                               power(C, e[2])), coefficient)
                     for e, coefficient in poly.items()))
    def coefficient_y(poly, degree):
        return {(e[0], 0, e[2]): coefficient for e, coefficient in poly.items()
                if e[1] == degree}
    centered_g = substitute_y(g, add(Y, scale(ONE, -1)))
    centered_H1 = substitute_y(H1, add(Y, scale(ONE, -1)))
    assert centered_g == power(Y, 4)
    assert coefficient_y(centered_H1, 22) == scale(power(C, 2), Q(3, 2))
    translated_g = substitute_y(centered_g, add(Y, X))
    translated_H1 = substitute_y(centered_H1, add(Y, X))
    topg, subg = coefficient_y(translated_g, 4), coefficient_y(translated_g, 3)
    toph, subh = coefficient_y(translated_H1, 23), coefficient_y(translated_H1, 22)
    # Difference subh/(23*toph)-subg/(4*topg)=3/(92*C).
    numerator = add(scale(mul(subh, topg), 4), scale(mul(subg, toph), -23))
    denominator = scale(mul(toph, topg), 92)
    assert numerator == scale(power(C, 2), 6)
    assert mul(scale(C, 92), numerator) == scale(denominator, 3)

    p, q, q1, q2, n0, n1 = 8, 12, 6, 23, 2, 4
    assert 1 < 2 < 4 and gcd(2, 3) == 1
    assert n1 == p // gcd(p, q1) and n1 * q1 // p == 3
    assert q1 > q * (n0 - 1) - p + 1
    M0 = p + q
    M1 = q1 - (q * (n0 - 1) - p)
    M2 = q2 - (n1 * q1 - M1)
    assert (M0, M1, M2) == (20, 2, 1)
    assert n0 * q - q1 == 18 and n1 * q1 - q2 == 1
    assert n0 * q - q1 + n1 * q1 - q2 == p + q - 1 == 19
    assert [gcd(p, q), gcd(gcd(p, q), q1), gcd(gcd(gcd(p, q), q1), q2)] == [4, 2, 1]

    a = add(ONE, X)
    Pstar = mul(power(a, 4), power(u, 8))
    Qstar = mul(power(a, 6), power(u, 12))
    Ustar = mul(power(a, 3), H)
    Vstar = mul(power(a, 12), H1)
    assert Vstar == add(power(Ustar, 4), scale(mul(power(C, 4), power(Pstar, 3)), -1))
    # These polynomial controls also satisfy the local face-ODE Poisson
    # identities; their first chop is deliberately NOT the required one.
    assert poisson(Ustar, Pstar) == scale(Qstar, 2)
    assert poisson(Vstar, Pstar) == scale(mul(Qstar, power(Ustar, 3)), 8)
    assert poisson(Qstar, Pstar) == {}
    assert add(power(Qstar, 2), scale(power(Pstar, 3), -1)) == {} and Ustar
    expected_x = [(4, power(g, 2)), (6, power(g, 3)), (3, H), (12, H1)]
    expected_y = [(8, power(a, 4)), (12, power(a, 6)),
                  (6, mul(C, power(a, 3))),
                  (23, scale(mul(power(C, 3), power(a, 12)), 2))]
    for poly, ex, ey in zip([Pstar, Qstar, Ustar, Vstar], expected_x, expected_y):
        assert face(poly, 0) == ex and face(poly, 1) == ey
        assert max(e[1] - e[0] for e in poly) == ey[0]

    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    receipt = {
        "status": "PASS",
        "claim": "Local one-root ODE homogeneous branch and next chopping/order/axis-face countercontrol",
        "environment": {"python": platform.python_version(), "arithmetic": "Independent Python integers/Fraction polynomials in x,y,C; C indeterminate"},
        "primary_source": {"url": "https://arxiv.org/html/2410.06959v5", "version_date": "2026-01-19", "printed_locations": ["Theorem4.2 prelude p44", "equation4.13 p45", "Case2b p46", "equation4.14 p47", "induction caseb p49", "Lemma6.1 p64"]},
        "parameters": {"l": 2, "d": 4, "n": 2, "m": 3, "beta": 1, "u": "1+y", "C": "indeterminate; order statements require C nonzero"},
        "identities": ["-2H'g+3Hg'=g^2", "H1=H^4-C^4g^6", "-2H1'g+12H1g'=4H^3g^2"],
        "H": "u^5/2+C*u^6", "H1": "u^20/16+C*u^21/2+3*C^2*u^22/2+2*C^3*u^23",
        "orders": [p, q, q1, q2], "chopping_exponents": [n0, n1],
        "M_values": [M0, M1, M2], "distances": [18, 1], "total_distance": 19,
        "sequential_gcds": [4, 2, 1],
        "axis_faces_and_bold_order": True,
        "simultaneous_centering_difference": "3/(92*C), nonzero for C!=0; invariant under common translations",
        "second_local_chop_polynomial_identity": True,
        "initial_chop_not_realized": "Qstar^2-Pstar^3=0 differs from nonzeroUstar",
        "hashes": {
            "ONE_ROOT_CHOPPING_ODE_HOMOGENEOUS_BRANCH.md": sha(HERE / "ONE_ROOT_CHOPPING_ODE_HOMOGENEOUS_BRANCH.md"),
            "check_one_root_chopping_ode.py": sha(Path(__file__)),
            "sources/zheglov-2410.06959v5.pdf": sha(HERE / "sources/zheglov-2410.06959v5.pdf"),
        },
        "scope": "Countercontrol to deriving a pure-power polynomial from the displayed local one-root ODE and next ODE/order/face constraints. Not a simultaneous Weyl pair, not a counterexample to Theorem4.2 under all global hypotheses, and not a main-proof or DC1 refutation.",
        "remaining_obligation": "Show an additional actual-pair/normal-form condition forces C=0, or repair the all-stage induction/terminal-edge claim while retaining C.",
        "formal_verification": False, "external_review": False,
    }
    out = args.output_dir or Path(tempfile.mkdtemp(prefix="one-root-chopping-ode-"))
    out.mkdir(parents=True, exist_ok=True)
    path = out / "one_root_chopping_ode_controls.json"
    path.write_text(json.dumps(receipt, indent=2) + "\n")
    print("PASS all-C one-root ODE, exact next chop, next ODE, orders/gcds/distances and individual axis faces")
    print("No simultaneous Weyl-pair claim; receipt:", path)


if __name__ == "__main__":
    main()
