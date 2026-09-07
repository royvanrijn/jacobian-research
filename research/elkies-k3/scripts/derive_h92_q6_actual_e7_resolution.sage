#!/usr/bin/env sage -python
"""Reject an untransported standard E7 chart as an H92 resolution map.

The formal E7 model ``Y^2=U^3+U*Z^3`` is useful for local class-group
bookkeeping, but it is not automatically an actual coordinate chart of the
H92 K3. This diagnostic begins with the exact H92 short Weierstrass germ at
the finite additive fibre and replays the first two ordinary blow-ups that
the standard-form calculation would use. The second Z-chart is smooth because
the exact coefficient B1 is nonzero. Hence the standard chart's ordinary node
cannot be treated as a transported H92 blow-up chart without an explicit
analytic coordinate map.

This is a positive safety certificate: it prevents a high-degree vertical
condition from being certified using a merely Kodaira-compatible chart.
"""

import argparse
import hashlib
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-resolution-diagnostic.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def strict(ring, value, substitutions, exceptional):
    transformed = ring(value(*substitutions))
    quotient, remainder = transformed.quo_rem(exceptional**2)
    assert not remainder
    return quotient


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

anchor = SourceFileLoader("h92_actual_e7_diagnostic_anchor", str(ANCHOR)).load_module()
r, s = anchor.EXPECTED_H92
_, formulas = anchor.parse_h92(H92)
A1, A, B1, B, B2 = (QQ(value(r, s)) for value in formulas)
assert B1

ring = PolynomialRing(QQ, names=("Z", "U", "Y"))
Z, U, Y = ring.gens()
f0 = Y**2 - U**3 - (A1 * Z**3 + A * Z**4) * U - (
    B1 * Z**5 + B * Z**6 + B2 * Z**7
)
f1_z = strict(ring, f0, (Z, Z * U, Z * Y), Z)
f2_z = strict(ring, f1_z, (Z, Z * U, Z * Y), Z)
assert f1_z == Y**2 - Z * U**3 - A1 * Z**2 * U - A * Z**3 * U - B1 * Z**3 - B * Z**4 - B2 * Z**5
assert f2_z == Y**2 - Z**2 * U**3 - A1 * Z * U - A * Z**2 * U - B1 * Z - B * Z**2 - B2 * Z**3
gradient = {name: str(f2_z.derivative(variable)(0, 0, 0))
            for name, variable in (("Z", Z), ("U", U), ("Y", Y))}
assert QQ(gradient["Z"]) == -B1
assert QQ(gradient["Z"])

payload = {
    "schema": "elkies-k3.h92-q6-actual-e7-resolution-diagnostic.v1",
    "status": "PASS_REJECTS_UNTRANSPORTED_STANDARD_E7_CHART",
    "inputs": {"h92_source": {"path": str(H92.relative_to(ROOT)), "sha256": digest(H92)}},
    "finite_h92_germ": str(f0),
    "standard_sequence_first_charts": {"blow1_Z": str(f1_z), "blow2_Z": str(f2_z)},
    "B1": str(B1),
    "blow2_Z_gradient_at_origin": gradient,
    "conclusion": (
        "The blow2_Z origin is smooth in the exact H92 germ. Therefore the "
        "ordinary node in the standard E7 chart is not an actual H92 chart "
        "under this untransported coordinate identification."
    ),
    "compiler_instruction": (
        "Do not use the formal E7 valuation or Cartier atlas to certify a "
        "high-degree H92 condition until an explicit coordinate transport from "
        "the H92 germ to that formal model has been supplied and verified."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92ACTUALE7|blow2_Z_dZ={}|status=PASS_REJECTS_UNTRANSPORTED_STANDARD_E7_CHART".format(
        -B1
    ),
    flush=True,
)
