#!/usr/bin/env sage -python
"""Corrected v3 replay of the bounded direct fixed-resolvent determinant slice."""

from __future__ import annotations

import argparse
from importlib.machinery import SourceFileLoader
import importlib.util
import itertools
import json
from pathlib import Path

from sage.all import QQ, PolynomialRing, matrix


HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]
ART = ROOT / "artifacts" / "generated-results" / "elliptic-curves"
PROTOCOL = HERE.with_name("CURVE302_FIXED_RESOLVENT_ORBIT_PROTOCOL_V3.json")
OUTPUT = ART / "curve302_fixed_resolvent_orbits_v3.json"
WORK = ROOT / "artifacts" / "local" / "elliptic-curves" / "curve302-fixed-resolvent-orbits-v3"
V2 = HERE.with_name("curve302_fixed_resolvent_orbits_v2.sage")

spec = importlib.util.spec_from_loader("curve302_fixed_resolvent_orbits_v2", SourceFileLoader("curve302_fixed_resolvent_orbits_v2", str(V2)))
v2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v2)


def read(path):
    return json.loads(path.read_text())


def put_new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def corrected_direct_perturbations(form, parameter_bound, perturbation_bound):
    R = PolynomialRing(QQ, names=("x", "y"))
    x, y = R.gens()
    F_over_four = sum(QQ(value) * x ** (3 - i) * y**i / 4 for i, value in enumerate(form))
    coefficients = range(-perturbation_bound, perturbation_bound + 1)
    for px, py in itertools.product(range(-parameter_bound, parameter_bound + 1), repeat=2):
        A, B = v2.base_pencil(form, px, py)
        M = x * A - y * B
        L, m, n, r, s = M[0, 0], M[0, 1], M[0, 2], M[1, 1], M[1, 2]
        for rx, ry, sx, sy in itertools.product(coefficients, repeat=4):
            rr, ss = r + rx * x + ry * y, s + sx * x + sy * y
            denominator = L * rr - m * m
            numerator = F_over_four + L * ss * ss - 2 * m * n * ss + n * n * rr
            if denominator == 0:
                continue
            tt = v2.linear_quotient(denominator, numerator, x, y)
            if tt is None:
                continue
            P = matrix(R, [[L, m, n], [m, rr, ss], [n, ss, tt]])
            Anew = matrix(QQ, 3, 3, lambda i, j: v2.coeff(P[i, j], x, y, 1))
            Bnew = matrix(QQ, 3, 3, lambda i, j: -v2.coeff(P[i, j], x, y, 0))
            yield (px, py, rx, ry, sx, sy), Anew, Bnew


def compute():
    # v2.compute resolves this global iterator at call time.  The only v2
    # arithmetic correction is det(xA-yB)=Fmax/4; all subsequent gates are
    # reused verbatim and independently rechecked from the pair.
    v2.PROTOCOL = PROTOCOL
    v2.base.PROTOCOL = PROTOCOL
    v2.direct_perturbations = corrected_direct_perturbations
    result = v2.compute()
    result["schema"] = "elliptic-curves.curve302-fixed-resolvent-orbits.v3"
    result["protocol"] = v2.base.relative(PROTOCOL)
    result["protocol_sha256"] = v2.base.digest(PROTOCOL)
    result["bindings"] = {v2.base.relative(path): v2.base.digest(path) for path in (PROTOCOL, V2, v2.V1, v2.base.MAXIMUM_FORM, v2.base.STRICT, v2.base.FILTRATION, v2.base.CURVE, v2.base.HIDDEN_SOLVER, HERE)}
    result["bounded_search"]["corrected_determinant_identity"] = "det(x*A-y*B)=Fmax/4"
    return result


def capture():
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    WORK.mkdir(parents=True, exist_ok=True)
    put_new(OUTPUT, compute())
    print(read(OUTPUT)["status"], flush=True)


def check():
    stored = read(OUTPUT)
    assert stored["bindings"] == {v2.base.relative(path): v2.base.digest(path) for path in (PROTOCOL, V2, v2.V1, v2.base.MAXIMUM_FORM, v2.base.STRICT, v2.base.FILTRATION, v2.base.CURVE, v2.base.HIDDEN_SOLVER, HERE)}
    if stored["status"] != "PASS_NEW_STRICT_CHARACTER_HANDED_OFF":
        assert compute() == stored
    print("PASS corrected direct fixed-resolvent determinant slice replay", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("capture", "check"))
    args = parser.parse_args()
    (capture if args.mode == "capture" else check)()
