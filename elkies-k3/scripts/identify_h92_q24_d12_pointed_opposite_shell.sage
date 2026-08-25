#!/usr/bin/env sage -python
"""Match the exact q24 pointed-opposite D12 section against known shell points."""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--output",
    type=Path,
    default=GENERATED / "elkies-k3-h3-q24-d12-pointed-opposite-shell-match.json",
)
args = parser.parse_args()

POINT = LOCAL / "q24-d12-pointed-opposite-section-qq.json"
A11 = LOCAL / "q24-d12-to-a11-orbit42-resolved-rr-qq.json"
ZERO = LOCAL / "q24-orbit42-rational-zero-pole-sections-qq.json"
SPINOR = LOCAL / "q24-orbit42-spinor-zero-pole-sections-qq.json"
INPUTS = (POINT, A11, ZERO, SPINOR)
for path in INPUTS:
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

point = json.loads(POINT.read_text())
a11 = json.loads(A11.read_text())
zero = json.loads(ZERO.read_text())
spinor = json.loads(SPINOR.read_text())
assert point["status"] == "PASS_EXACT_Q24_D12_POINTED_OPPOSITE_SECTION_QQ"
assert zero["status"] == "PASS_EXACT_Q42_RATIONAL_ZERO_POLE_SECTIONS_QQ"
assert spinor["status"] == "PASS_EXACT_Q42_SPINOR_ZERO_POLE_SECTIONS_QQ"

UQ = PolynomialRing(QQ, "u")
KU = UQ.fraction_field()
VQ = PolynomialRing(QQ, "V")
KV = VQ.fraction_field()
u_of_V = KV(a11["coordinate_change"]["u_of_V"])
x_scale = KV(a11["coordinate_change"]["x_scale"])
y_scale = KV(a11["coordinate_change"]["y_scale"])


def evaluate_u(value):
    value = KU(value)
    return KV(VQ(value.numerator())(u_of_V)) / KV(VQ(value.denominator())(u_of_V))


section = point["section"]
X = VQ([QQ(value) for value in section["X_coefficients_low_to_high"]])
Y = VQ([QQ(value) for value in section["Y_coefficients_low_to_high"]])
Z = VQ([QQ(value) for value in section["Z_coefficients_low_to_high"]])
x = KV(X) / KV(Z**2)
y = KV(Y) / KV(Z**3)

matches = []
for kind, rows in (("identity", zero["sections"]), ("spinor", spinor["sections"])):
    for index, row in enumerate(rows):
        xu = UQ([QQ(value) for value in row["x_coefficients_low_to_high"]])
        yu = UQ([QQ(value) for value in row["y_coefficients_low_to_high"]])
        candidate_x = x_scale * evaluate_u(KU(xu))
        candidate_y = y_scale * evaluate_u(KU(yu))
        if x == candidate_x and y in (candidate_y, -candidate_y):
            matches.append(
                {
                    "shell_kind": kind,
                    "section_index": index,
                    "same_y_sign": y == candidate_y,
                    "pair_index": row.get("pair_index"),
                    "stored_sign": row.get("sign"),
                }
            )

payload = {
    "schema": "elkies-k3.h3-q24-d12-pointed-opposite-shell-match.v1",
    "status": (
        "PASS_EXACT_Q24_D12_POINTED_OPPOSITE_SHELL_MATCH"
        if matches else
        "Q24_D12_POINTED_OPPOSITE_NOT_IN_KNOWN_ZERO_POLE_SHELL"
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
    "match_count": len(matches),
    "matches": matches,
    "proof_boundary": (
        "Exact QQ(V) coordinate comparison with all eighteen identity and two "
        "spinor zero-pole points after the certified shell-to-minimal change."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q24D12POINTEDSHELL|matches={}|labels={}|status={}".format(
        len(matches),
        ",".join(f"{row['shell_kind']}:{row['section_index']}" for row in matches) or "none",
        payload["status"],
    ),
    flush=True,
)
print(f"OUTPUT|{args.output.resolve()}", flush=True)
