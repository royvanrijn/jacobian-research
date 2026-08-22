#!/usr/bin/env sage -python
"""Exact H92 identity check for a reconstructed height-46 P2 section.

The initial implementation accepted only the multi-prime CRT representation.
The p-adic Newton lift records the same section more compactly as
``x=X/Z^2, y=Y/Z^3``; accepting both formats keeps the equation-level
transport independent of the reconstruction method.
"""

from pathlib import Path
import json
from importlib.machinery import SourceFileLoader

from sage.all import PolynomialRing, QQ
import argparse


ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--input",
    type=Path,
    default=ROOT / "artifacts/generated-results/elkies-k3-h92-p2-half-crt.json",
)
args = parser.parse_args()
payload = json.loads(args.input.read_text())
if not payload["complete"]:
    raise RuntimeError("the P2 reconstruction is incomplete")
anchor = SourceFileLoader(
    "h92_p2_anchor", str(ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage")
).load_module()
r, s = anchor.EXPECTED_H92
_, formulas = anchor.parse_h92(ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt")
A1, A, B1, B, B2 = (QQ(value(r, s)) for value in formulas)

ring = PolynomialRing(QQ, "t")
t = ring.gen()
field = ring.fraction_field()


def polynomial(coefficients):
    return sum(field(QQ(value)) * t**index for index, value in enumerate(coefficients))


if payload["schema"] == "elkies-k3.h92-p2-half-crt.v1":
    def coordinate(name):
        return polynomial(payload[name]["numerator"]) / polynomial(payload[name]["denominator"])
elif payload["schema"] == "elkies-k3.h92-p2-hensel-lift.v1":
    z = polynomial(payload["Z"])
    def coordinate(name):
        return polynomial(payload[name.upper()]) / z**(2 if name == "x" else 3)
else:
    raise ValueError("unsupported P2 reconstruction schema {}".format(payload["schema"]))


x, y = coordinate("x"), coordinate("y")
A_t = A1*t**3 + A*t**4
B_t = B1*t**5 + B*t**6 + B2*t**7
assert y**2 == x**3 + A_t*x + B_t
assert (x.numerator().degree(), x.denominator().degree()) == (46, 42)
assert (y.numerator().degree(), y.denominator().degree()) == (69, 63)
print("H92P2COORD|identity=PASS|x_degrees=46,42|y_degrees=69,63|schema={}".format(payload["schema"]))
