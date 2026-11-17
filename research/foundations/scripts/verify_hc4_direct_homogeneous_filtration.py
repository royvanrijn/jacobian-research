#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "generated-results" / "hc4_direct_homogeneous_filtration.json"
EXPECTED_OUT_SHA256 = "9731d3636a0aee74b9573e051be3ef55a0682d227915aff03ad60b039806651a"


def audit_existing() -> None:
    actual = hashlib.sha256(OUT.read_bytes()).hexdigest()
    assert actual == EXPECTED_OUT_SHA256, (
        f"committed filtration artifact drifted: expected {EXPECTED_OUT_SHA256}, "
        f"got {actual}"
    )
    payload = json.loads(OUT.read_text(encoding="utf-8"))
    assert payload == {
        "scope": "direct HC4 top homogeneous filtration",
        "status": "verified",
        "identities": {
            "first": "[t] det M = c1 det(A0)",
            "second_when_c1_zero": (
                "[t^2] det M = c2 det(A0) - b1^T adj(A0) b1"
            ),
        },
    }
    print(
        "PASS: committed HC4 direct-filtration artifact is intact; "
        "no symbolic replay"
    )


parser = argparse.ArgumentParser()
parser.add_argument(
    "--audit-existing-only",
    action="store_true",
    help="validate the committed artifact without symbolic replay or rewriting it",
)
arguments = parser.parse_args()
if arguments.audit_existing_only:
    audit_existing()
    raise SystemExit(0)

import sympy as sp

OUT.parent.mkdir(parents=True, exist_ok=True)

t = sp.symbols("t")
a11, a12, a13, a22, a23, a33 = sp.symbols("a11 a12 a13 a22 a23 a33")
u11, u12, u13, u22, u23, u33 = sp.symbols("u11 u12 u13 u22 u23 u33")
v11, v12, v13, v22, v23, v33 = sp.symbols("v11 v12 v13 v22 v23 v33")
b11, b12, b13, b21, b22, b23, c1, c2 = sp.symbols(
    "b11 b12 b13 b21 b22 b23 c1 c2"
)

A0 = sp.Matrix([[a11, a12, a13], [a12, a22, a23], [a13, a23, a33]])
A1 = sp.Matrix([[u11, u12, u13], [u12, u22, u23], [u13, u23, u33]])
A2 = sp.Matrix([[v11, v12, v13], [v12, v22, v23], [v13, v23, v33]])
b1 = sp.Matrix([b11, b12, b13])
b2 = sp.Matrix([b21, b22, b23])

At = A0 + t * A1 + t**2 * A2
bt = t * b1 + t**2 * b2
ct = t * c1 + t**2 * c2
M = sp.Matrix.vstack(
    sp.Matrix.hstack(At, bt),
    sp.Matrix.hstack(bt.T, sp.Matrix([[ct]])),
)

poly = sp.Poly(sp.expand(M.det()), t)
detA = sp.factor(A0.det())
coef1 = sp.factor(poly.coeff_monomial(t))
coef2 = sp.factor(poly.coeff_monomial(t**2))

assert sp.expand(coef1 - c1 * detA) == 0
expected2 = c2 * detA - (b1.T * A0.adjugate() * b1)[0]
assert sp.expand(coef2.subs(c1, 0) - expected2) == 0

result = {
    "scope": "direct HC4 top homogeneous filtration",
    "status": "verified",
    "identities": {
        "first": "[t] det M = c1 det(A0)",
        "second_when_c1_zero": "[t^2] det M = c2 det(A0) - b1^T adj(A0) b1",
    },
}
OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
