#!/usr/bin/env sage -python
"""Reject the tempting old-(t,x,y)-monomial E7 quotient for the q=6 third divisor.

The actual vertical target for the degree-44 third transported section is the
anti-nef E7 cycle c with C_E7*c=22*e5, hence its complete ideal must have
colength c.C_E7.c/2=363.  It is tempting to impose that cycle using only
monomials t^i x^a y^b in the singular old Weierstrass coordinates.  The
actual H92 valuation atlas lets us test that proposal exactly.  It gives a
monomial quotient of length 529, not 363.

Therefore old-coordinate valuations alone do not supply the desired complete
ideal: the missing generators/relations must come from actual resolved-chart
transition functions.  This is a regression guard against turning a
component-order calculation into a false condition matrix.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
ATLAS = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-valuation-atlas.json"
TARGET = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-third-e7-local-target.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-third-rejected-old-monomial-e7-ideal.json"
ATLAS_SHA256 = "ae7eb1e79a2fb41ab05d0092cb8c04663307bd8aa1e0cb562a8d3a014e94f451"
TARGET_SHA256 = "a0699e4ec75930cc93a9706ddf96f4ffc744954809e53a24029bd8c6668843f7"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--atlas", type=Path, default=ATLAS)
parser.add_argument("--target", type=Path, default=TARGET)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

if args.atlas == ATLAS:
    assert digest(args.atlas) == ATLAS_SHA256
if args.target == TARGET:
    assert digest(args.target) == TARGET_SHA256
atlas = json.loads(args.atlas.read_text())
target = json.loads(args.target.read_text())
assert atlas["status"] == "PASS_EXACT_H92_E7_VALUATION_ATLAS"
assert target["status"] == "PASS_EXACT_Q6_THIRD_E7_LATTICE_TARGET"

orders = [entry["old_coordinate_orders"] for entry in atlas["entries"]]
v_t = vector(ZZ, [entry["t"] for entry in orders])
v_x = vector(ZZ, [entry["x"] for entry in orders])
v_y = vector(ZZ, [entry["y"] for entry in orders])
c = -vector(ZZ, target["resolved_exceptional_coefficients"])
assert c == vector(ZZ, (22, 44, 66, 44, 33, 33, 55))
cartan = matrix(ZZ, [
    [2, 0, 0, -1, 0, 0, 0],
    [0, 2, 0, 0, -1, 0, -1],
    [0, 0, 2, -1, 0, -1, -1],
    [-1, 0, -1, 2, 0, 0, 0],
    [0, -1, 0, 0, 2, 0, 0],
    [0, 0, -1, 0, 0, 2, 0],
    [0, -1, -1, 0, 0, 0, 2],
])
assert cartan*c == vector(ZZ, (0, 0, 0, 0, 22, 0, 0))
expected_colength = ZZ(c*cartan*c)//2
assert expected_colength == 363


def in_naive_ideal(exponents):
    i, a, b = exponents
    return min(i*v_t+a*v_x+b*v_y-c) >= 0


# The minimum old-coordinate t valuation is one, and x,y have strictly
# positive valuations.  Hence any nonmember has i<33, a<17, b in {0,1};
# this finite box counts the full quotient in the Weierstrass normal form.
candidate_generators = [
    (i, a, b)
    for b in range(2)
    for a in range(18)
    for i in range(34)
    if in_naive_ideal((i, a, b))
]
minimal_generators = [
    exponent for exponent in candidate_generators
    if not any(
        predecessor != exponent
        and all(left <= right for left, right in zip(predecessor, exponent))
        for predecessor in candidate_generators
    )
]
quotient_basis = [
    (i, a, b)
    for b in range(2)
    for a in range(17)
    for i in range(33)
    if not in_naive_ideal((i, a, b))
]
naive_colength = len(quotient_basis)
assert naive_colength == 529
assert naive_colength != expected_colength

payload = {
    "schema": "elkies-k3.h92-q6-third-rejected-old-monomial-e7-ideal.v1",
    "status": "PASS_REJECTS_OLD_MONOMIAL_E7_COMPLETE_IDEAL",
    "inputs": {
        "actual_valuation_atlas": {"path": str(args.atlas.relative_to(ROOT)), "sha256": digest(args.atlas)},
        "third_lattice_target": {"path": str(args.target.relative_to(ROOT)), "sha256": digest(args.target)},
    },
    "cycle": [int(value) for value in c],
    "cartan_boundary": [int(value) for value in cartan*c],
    "complete_ideal_expected_colength": int(expected_colength),
    "naive_old_monomial_ideal": {
        "membership": "t^i*x^a*y^b with b in {0,1} and i*v_t+a*v_x+b*v_y >= cycle componentwise",
        "minimal_generators": [[int(value) for value in item] for item in minimal_generators],
        "quotient_colength": int(naive_colength),
    },
    "conclusion": "529 != 363, so a condition matrix built from this old-coordinate monomial ideal would be wrong. Derive actual chart transition generators before forming the finite quotient.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6THIRDREJECTMONOMIAL|naive_colength=529|expected=363|"
    "status=PASS_REJECTS_OLD_MONOMIAL_E7_COMPLETE_IDEAL",
    flush=True,
)
