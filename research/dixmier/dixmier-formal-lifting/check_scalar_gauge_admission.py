#!/usr/bin/env python3
"""Finite Laurent-Ore controls for scalar gauges and the subprincipal corner."""
from fractions import Fraction as Q
from pathlib import Path
from math import comb
import argparse
import hashlib
import itertools
import json
import platform
import tempfile

HERE = Path(__file__).resolve().parent
ONE = {(0, 0): Q(1)}
U = {(1, 0): Q(1)}
D = {(0, 1): Q(1)}


def clean(p):
    return {e: c for e, c in p.items() if c}


def add(*polys):
    out = {}
    for p in polys:
        for e, c in p.items():
            out[e] = out.get(e, Q(0)) + c
    return clean(out)


def scale(p, c):
    return clean({e: a * c for e, a in p.items()})


def falling(a, n):
    out = 1
    for j in range(n):
        out *= a - j
    return out


def ore(a, b):
    out = {}
    for (i, j), c in a.items():
        for (k, l), d in b.items():
            for n in range(j + 1):
                e = (i + k - n, j + l - n)
                out[e] = out.get(e, Q(0)) + c * d * comb(j, n) * falling(k, n)
    return clean(out)


def power(p, n):
    out = ONE
    for _ in range(n):
        out = ore(out, p)
    return out


def comm(a, b):
    return add(ore(a, b), scale(ore(b, a), -1))


def gauge(poly, gamma):
    shift = add(D, {(-1, 0): gamma})
    out = {}
    for (i, j), c in poly.items():
        out = add(out, scale(ore({(i, 0): Q(1)}, power(shift, j)), c))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    assert comm(D, U) == ONE
    # Separate Weyl-symbol convention control: average all six words
    # with two positions and two derivatives.
    weyl = {}
    for positions in itertools.combinations(range(4), 2):
        word = ONE
        for i in range(4):
            word = ore(word, U if i in positions else D)
        weyl = add(weyl, scale(word, Q(1, 6)))
    assert weyl == {(2, 2): Q(1), (1, 1): Q(2), (0, 0): Q(1, 2)}

    corner_controls = []
    for a, b, c, d in [(8, 4, 12, 6), (6, 4, 12, 8), (4, 6, 6, 9)]:
        rho = Q(c, a)
        assert a * d == b * c and a != b
        for A, B in [(Q(49, 3), Q(73, 2)), (Q(2), Q(-3))]:
            P = {(a, b): Q(1), (a - 1, b - 1): A}
            QQ = {(c, d): Q(1), (c - 1, d - 1): B}
            AW = A - Q(a * b, 2)
            BW = B - Q(c * d, 2)
            expected = (a - b) * (rho * AW - BW)
            actual = comm(QQ, P).get((a + c - 2, b + d - 2), Q(0))
            assert actual == expected
            shifts = []
            for gamma in [Q(-2), Q(1, 2), Q(3)]:
                Pg, Qg = gauge(P, gamma), gauge(QQ, gamma)
                assert Pg.get((a - 1, b - 1), Q(0)) == A + b * gamma
                assert Qg.get((c - 1, d - 1), Q(0)) == B + d * gamma
                assert rho * (AW + b * gamma) - (BW + d * gamma) == rho * AW - BW
                assert gauge(comm(QQ, P), gamma) == comm(Qg, Pg)
                shifts.append(str(gamma))
            corner_controls.append({"corners": [a, b, c, d], "A": str(A), "B": str(B),
                                    "commutator_coefficient": str(actual), "gamma": shifts})
    genuine_pair = []
    for gamma in [Q(-2), Q(1, 2), Q(1), Q(3)]:
        gauged_D = gauge(D, gamma)
        assert gauged_D == {(0, 1): Q(1), (-1, 0): gamma}
        assert comm(gauged_D, U) == ONE
        assert any(i < 0 for i, j in gauged_D)
        genuine_pair.append({"gamma": str(gamma), "relation_preserved": True,
                             "D_loses_polynomial_admission": True})
    for m in [-2, 1, 3]:
        poly = {(3, 2): Q(1), (1, 0): Q(2)}
        assert gauge(poly, Q(m)) == ore(ore({(-m, 0): Q(1)}, poly), {(m, 0): Q(1)})

    def mm(a, b):
        return [[add(*(ore(a[i][k], b[k][j]) for k in range(2)))
                 for j in range(2)] for i in range(2)]
    def madd(a, b, factor=Q(1)):
        return [[add(a[i][j], scale(b[i][j], factor)) for j in range(2)]
                for i in range(2)]
    eye = [[ONE, {}], [{}, ONE]]
    xpos = [[U, {}], [{}, U]]
    deriv = [[D, {}], [{}, D]]
    G = [[ONE, U], [{}, ONE]]
    Ginv = [[ONE, scale(U, -1)], [{}, ONE]]
    assert mm(G, Ginv) == eye and mm(Ginv, G) == eye
    qgood = mm(mm(Ginv, deriv), G)
    assert qgood == [[D, ONE], [{}, D]]
    assert mm(mm(Ginv, xpos), G) == xpos
    assert madd(mm(qgood, xpos), mm(xpos, qgood), -1) == eye
    pole = {(-1, 0): Q(1)}
    Gbad = [[U, {}], [{}, pole]]
    Gbadinv = [[pole, {}], [{}, U]]
    qbad = mm(mm(Gbadinv, deriv), Gbad)
    assert qbad == [[add(D, pole), {}], [{}, add(D, scale(pole, -1))]]
    assert madd(mm(qbad, xpos), mm(xpos, qbad), -1) == eye
    assert any(i < 0 for row in qbad for poly in row for i, j in poly)

    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    receipt = {
        "status": "PASS",
        "claim": "Finite Ore controls for formal scalar-unit corner shifts and a written rational scalar-gauge polynomial-admission obstruction",
        "environment": {"python": platform.python_version(), "arithmetic": "Independent Python integers/Fraction finite Ore multiplication over Laurent coefficients"},
        "Weyl_ordering_control": "Sym(u^2D^2)=u^2D^2+2uD+1/2",
        "corner_controls": corner_controls,
        "genuine_pair_controls": genuine_pair,
        "integer_conjugation_identity": [-2, 1, 3],
        "written_global_theorem": "For P,Q in A1(k), [Q,P]=kappa nonzero, rational f with bothf^-1Pf,f^-1Qf polynomial must be constant; proof uses invariant O,fO,intersection and finite quotient traces.",
        "quotient_dimensions": ["deg numerator", "deg denominator"],
        "written_matrix_position_extension": "For a polynomial Mat_r(A1) pair with scalar nonzero commutator, any rational position matrix gauge preserving polynomial admission lies in GL_r(k[x]); same invariant intersection and two finite quotient traces.",
        "matrix_controls": {"positive_unimodular": "[[1,u],[0,1]] sends DI to DI+E12", "negative_nonunimodular": "diag(u,u^-1) sends DI to diag(D+u^-1,D-u^-1) and loses polynomial admission"},
        "hashes": {
            "SCALAR_GAUGE_SUBPRINCIPAL_AND_POLYNOMIAL_ADMISSION.md": sha(HERE / "SCALAR_GAUGE_SUBPRINCIPAL_AND_POLYNOMIAL_ADMISSION.md"),
            "check_scalar_gauge_admission.py": sha(Path(__file__)),
        },
        "scope": "Corner controls are not Weyl pairs. Rational scalar and one-position matrix admission theorems are proved in writing, not inferred from tests. No absolute zero-subprincipal theorem, nonrational formal/algebraic gauge classification, multivariable extension or full Weyl-division gauge exclusion.",
        "formal_verification": False, "external_review": False,
    }
    out = args.output_dir or Path(tempfile.mkdtemp(prefix="scalar-gauge-admission-"))
    out.mkdir(parents=True, exist_ok=True)
    destination = out / "scalar_gauge_admission_controls.json"
    destination.write_text(json.dumps(receipt, indent=2) + "\n")
    print("PASS exact scalar-gauge corner, Weyl-ordering and genuine-pair polynomial-loss controls")
    print("Receipt:", destination)


if __name__ == "__main__":
    main()
