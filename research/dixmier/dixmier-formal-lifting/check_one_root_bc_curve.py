#!/usr/bin/env python3
"""Exact one-root BC-curve/Hensel/trace/head controls; no ODO admission claim."""
from fractions import Fraction as Q
from pathlib import Path
from math import gcd
import argparse
import hashlib
import json
import platform
import tempfile

HERE = Path(__file__).resolve().parent
# Generic Laurent monomials in three variables. At first these are z,R,c;
# the independent curve-polynomial test reuses them as X,Y,c.
ONE = {(0, 0, 0): Q(1)}


def clean(poly):
    return {e: c for e, c in poly.items() if c}


def add(*polys):
    out = {}
    for poly in polys:
        for e, c in poly.items():
            out[e] = out.get(e, Q(0)) + c
    return clean(out)


def scale(poly, c):
    return clean({e: a * c for e, a in poly.items()})


def mul(a, b):
    out = {}
    for e, c in a.items():
        for f, d in b.items():
            ef = tuple(i + j for i, j in zip(e, f))
            out[ef] = out.get(ef, Q(0)) + c * d
    return clean(out)


def power(poly, n):
    out = ONE
    for _ in range(n):
        out = mul(out, poly)
    return out


def shift(poly, exponent):
    return {tuple(i + j for i, j in zip(e, exponent)): c for e, c in poly.items()}


def coefficient(poly, index, degree):
    return {tuple(0 if j == index else n for j, n in enumerate(e)): c
            for e, c in poly.items() if e[index] == degree}


def derivative(poly, index):
    out = {}
    for e, c in poly.items():
        if e[index]:
            f = list(e)
            f[index] -= 1
            out[tuple(f)] = c * e[index]
    return out


def substitute(poly, index, value):
    out = {}
    for e, c in poly.items():
        assert e[index] >= 0
        residual = list(e)
        residual[index] = 0
        out = add(out, scale(shift(power(value, e[index]), residual), c))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    z = {(1, 0, 0): Q(1)}
    R = {(0, 1, 0): Q(1)}
    c = {(0, 0, 1): Q(1)}
    XX = {(-8, 0, 0): Q(1)}
    YY = add({(-12, 0, 0): Q(1)}, {(6, 0, 1): Q(1, 2)}, {(7, 1, 0): Q(1)})
    ZZ = add(power(YY, 2), scale(power(XX, 3), -1))
    AA = add(power(ZZ, 2), scale(mul(power(c, 2), YY), -1))
    FF = add(power(AA, 2), scale(mul(mul(c, power(XX, 2)), ZZ), -1))
    G = shift(FF, (22, 0, 0))
    assert all(e[0] >= 0 for e in G)
    expected0 = mul(power(c, 2), add(scale(power(R, 2), 16), scale(ONE, -1)))
    expected1 = add(scale(mul(c, power(R, 3)), 32), scale(mul(c, R), -2))
    expected2 = scale(power(R, 4), 16)
    assert coefficient(G, 0, 0) == expected0
    assert coefficient(G, 0, 1) == expected1
    assert coefficient(G, 0, 2) == expected2
    r0 = scale(ONE, Q(1, 4))
    assert substitute(expected0, 1, r0) == {}
    assert substitute(coefficient(derivative(G, 1), 0, 0), 1, r0) == scale(power(c, 2), 8)
    assert substitute(expected1, 1, r0) == {}
    # The simple derivative forces R1=0. Verify the next coefficient
    # using R=1/4-z^2/(128c^2); c remains indeterminate and nonzero.
    r2 = add(r0, {(2, 0, -2): Q(-1, 128)})
    G2 = substitute(G, 1, r2)
    assert all(not coefficient(G2, 0, j) for j in [0, 1, 2])

    # The finite truncation supplies all coefficients affecting the
    # displayed two highest chopped terms, but is not asserted to solveF.
    yshort = substitute(YY, 1, r0)
    q1 = add(power(yshort, 2), scale(power(XX, 3), -1))
    q2 = add(power(q1, 4), scale(mul(power(c, 4), power(XX, 3)), -1))
    assert min(e[0] for e in q1) == -6
    assert coefficient(q1, 0, -6) == c
    assert coefficient(q1, 0, -5) == scale(ONE, Q(1, 2))
    assert coefficient(q1, 0, -4) == {}
    assert min(e[0] for e in q2) == -23
    assert coefficient(q2, 0, -23) == scale(power(c, 3), 2)
    assert coefficient(q2, 0, -22) == scale(power(c, 2), Q(3, 2))
    assert coefficient(yshort, 0, 8) == {}

    # Independent polynomial variables X,Y,c for monicity and BC top.
    x = {(1, 0, 0): Q(1)}
    y = {(0, 1, 0): Q(1)}
    zz = add(power(y, 2), scale(power(x, 3), -1))
    curve = add(power(add(power(zz, 2), scale(mul(power(c, 2), y), -1)), 2),
                scale(mul(mul(c, power(x, 2)), zz), -1))
    assert max(e[1] for e in curve) == 8
    assert coefficient(curve, 1, 8) == ONE
    assert coefficient(curve, 1, 7) == {}
    assert max(8 * e[0] + 12 * e[1] for e in curve) == 96
    top = {e: a for e, a in curve.items() if 8 * e[0] + 12 * e[1] == 96}
    assert top == power(zz, 4)
    assert gcd(8, 7) == 1
    assert [gcd(8, 12), gcd(gcd(8, 12), 6), gcd(gcd(gcd(8, 12), 6), 23)] == [4, 2, 1]
    assert 2 * 12 - 6 == 18 and 4 * 6 - 23 == 1 and 18 + 1 == 8 + 12 - 1

    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    receipt = {
        "status": "PASS",
        "claim": "Explicit monic one-place BC-type curve retains the nonzero homogeneous parameter and trace constraint",
        "environment": {"python": platform.python_version(), "arithmetic": "Independent integers/Fraction Laurent polynomials; c indeterminate"},
        "curve": "Z=Y^2-X^3; F=(Z^2-c^2Y)^2-cX^2Z",
        "hypotheses": "Characteristic-zero field, c nonzero; finite constant extension for eighth roots of unity",
        "BC_data": {"degree_Y": 8, "weighted_top": "(Y^2-X^3)^4", "weights_X_Y": [8, 12], "Y7_coefficient": "0"},
        "normalized_equation": {"G0": "c^2*(16R^2-1)", "G1": "32cR^3-2cR", "G2": "16R^4", "simple_root": "R0=1/4", "derivative": "8c^2", "R1": "0", "R2": "-1/(128c^2)"},
        "branch": "X=t^8; Y=t^12+(c/2)t^-6+(1/4)t^-7+0*t^-8-t^-9/(128c^2)+O(t^-10)",
        "all_orders_branch": "Written formal implicit-function recursion; not finite extrapolation",
        "one_place_and_trace": "Written eight-distinct-conjugates argument from exponent-7, monic degree8 and Vieta",
        "chopping": {"orders": [8, 12, 6, 23], "q1_leading": "c", "q2_leading": "2c^3", "q2_next": "3c^2/2", "distances": [18, 1], "sequential_gcds": [4, 2, 1]},
        "hashes": {
            "ONE_ROOT_BC_CURVE_AND_ADMISSION_GAP.md": sha(HERE / "ONE_ROOT_BC_CURVE_AND_ADMISSION_GAP.md"),
            "check_one_root_bc_curve.py": sha(Path(__file__)),
            "sources/zheglov-2410.06959v5.pdf": sha(HERE / "sources/zheglov-2410.06959v5.pdf"),
        },
        "scope": "Algebraic curve and commutative formal Laurent data only; no chosen Sato datum, commuting ordinary-differential reconstruction, string-pair reconstruction or simultaneous polynomial Weyl admission. No DC1 theorem or refutation.",
        "formal_verification": False, "external_review": False,
    }
    out = args.output_dir or Path(tempfile.mkdtemp(prefix="one-root-bc-curve-"))
    out.mkdir(parents=True, exist_ok=True)
    path = out / "one_root_bc_curve_controls.json"
    path.write_text(json.dumps(receipt, indent=2) + "\n")
    print("PASS exact BC polynomial, simple Hensel branch coefficients, trace/head controls")
    print("No differential-operator admission claim; receipt:", path)


if __name__ == "__main__":
    main()
