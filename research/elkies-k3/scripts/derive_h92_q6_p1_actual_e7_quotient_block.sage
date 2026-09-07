#!/usr/bin/env sage -python
"""Construct the first q=6 E7 finite quotient in actual H92 coordinates.

Put ``T=t`` and translate the old coordinate by the exact branch jet

    U=x-c2*T^2-c3*T^3,  Y=y,

where c2=-B1/A1 and c3=-(c2^3+A*c2+B)/A1.  The actual resolved valuation
atlas gives

    v(T)=(2,2,4,3,1,2,3),
    v(U)=(2,4,6,4,3,3,5).

Their sum c has E7-Cartan boundary e1+2e5 and colength six.  Directly in
the translated *H92* local ring, its complete ideal is

    (U^2, T*U, T^4, T*Y, U*Y),

with quotient basis 1,T,T^2,T^3,U,Y.  The exact H92 equation reduces to
zero in this quotient. Combined with the corrected actual +/-P1 marked-chart
calculation, this is a finite q=6 E7 quotient block. It is not an all-edge
transition module: the global q=6 cover still requires its resolved-chart
gluing evaluator.
"""

import argparse
import hashlib
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import PolynomialRing, QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
ATLAS = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-valuation-atlas.json"
MARKED = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-p1-actual-e7-marked-module-corrected.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-p1-actual-e7-quotient-block.json"
ATLAS_SHA256 = "ae7eb1e79a2fb41ab05d0092cb8c04663307bd8aa1e0cb562a8d3a014e94f451"
MARKED_SHA256 = "4a94a5aca8686fbb666b5bb26b6f784eca33079d8922da0766ec5bd0ae2a4ba8"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--atlas", type=Path, default=ATLAS)
parser.add_argument("--marked-module", type=Path, default=MARKED)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

if args.atlas == ATLAS:
    assert digest(args.atlas) == ATLAS_SHA256
if args.marked_module == MARKED:
    assert digest(args.marked_module) == MARKED_SHA256
atlas = json.loads(args.atlas.read_text())
marked = json.loads(args.marked_module.read_text())
assert atlas["status"] == "PASS_EXACT_H92_E7_VALUATION_ATLAS"
assert marked["status"] == "PASS_EXACT_P1_ACTUAL_E7_MARKED_MODULE_CORRECTED"
assert marked["marked_component"] == "E7_5"

anchor = SourceFileLoader("h92_q6_actual_e7_quotient_anchor", str(ANCHOR)).load_module()
r, s = anchor.EXPECTED_H92
_, formulas = anchor.parse_h92(H92)
A1, A, B1, B, B2 = (QQ(value(r, s)) for value in formulas)
c2 = -B1/A1
c3 = -(c2**3+A*c2+B)/A1

ring = PolynomialRing(QQ, names=("T", "U", "Y"), order="degrevlex")
T, U, Y = ring.gens()
x_branch = c2*T**2+c3*T**3
relation = Y**2-(x_branch+U)**3-(A1*T**3+A*T**4)*(x_branch+U)-(
    B1*T**5+B*T**6+B2*T**7
)
ambient_complete_ideal = ring.ideal((U**2, T*U, T**4, T*Y, U*Y, Y**2))
assert relation.reduce(ambient_complete_ideal.groebner_basis()) == 0
surface_ideal = ring.ideal((relation, U**2, T*U, T**4, T*Y, U*Y))
assert surface_ideal.vector_space_dimension() == 6
quotient_basis = (ring(1), T, T**2, T**3, U, Y)
for value in quotient_basis:
    assert value.reduce(surface_ideal.groebner_basis()) == value

# This U valuation is read from the actual edge maps.  The only non-obvious
# cancellation occurs on E7_5: there t=Z*U0*(Z-A1/B1), x=Z^2*U0^2*(Z-A1/B1),
# and c2*(-A1/B1)=1, so x-c2*t^2 has exactly one additional Z factor.
v_t = vector(ZZ, atlas["old_base_fibre_multiplicities"])
v_u = vector(ZZ, (2, 4, 6, 4, 3, 3, 5))
assert v_t == vector(ZZ, (2, 2, 4, 3, 1, 2, 3))
cycle = v_t+v_u
cartan = matrix(ZZ, [
    [2, 0, 0, -1, 0, 0, 0],
    [0, 2, 0, 0, -1, 0, -1],
    [0, 0, 2, -1, 0, -1, -1],
    [-1, 0, -1, 2, 0, 0, 0],
    [0, -1, 0, 0, 2, 0, 0],
    [0, 0, -1, 0, 0, 2, 0],
    [0, -1, -1, 0, 0, 0, 2],
])
assert cycle == vector(ZZ, (4, 6, 10, 7, 4, 5, 8))
assert cartan*cycle == vector(ZZ, (1, 0, 0, 0, 2, 0, 0))
assert cycle*cartan*cycle//2 == 6

# Every ideal generator reaches the required divisorial cycle.  In addition
# to proving membership, the matching colength proves the complete ideal has
# no unrecorded generator in the actual rational-double-point local ring.
v_y = vector(ZZ, [entry["old_coordinate_orders"]["y"] for entry in atlas["entries"]])
for valuation in (2*v_u, v_t+v_u, 4*v_t, v_t+v_y, v_u+v_y, 2*v_y):
    assert min(valuation-cycle) >= 0

payload = {
    "schema": "elkies-k3.h92-q6-p1-actual-e7-quotient-block.v1",
    "status": "PASS_EXACT_Q6_P1_ACTUAL_E7_QUOTIENT_BLOCK",
    "inputs": {
        "h92_source": {"path": str(H92.relative_to(ROOT)), "sha256": digest(H92)},
        "actual_valuation_atlas": {"path": str(args.atlas.relative_to(ROOT)), "sha256": digest(args.atlas)},
        "actual_marked_module": {"path": str(args.marked_module.relative_to(ROOT)), "sha256": digest(args.marked_module)},
    },
    "translated_coordinates": {"T": "t", "U": "x-c2*t^2-c3*t^3", "Y": "y", "c2": str(c2), "c3": str(c3)},
    "actual_exceptional_orders": {"T": [int(value) for value in v_t], "U": [int(value) for value in v_u], "Y": [int(value) for value in v_y]},
    "cycle": [int(value) for value in cycle],
    "cartan_boundary": [int(value) for value in cartan*cycle],
    "complete_ideal": ["U^2", "T*U", "T^4", "T*Y", "U*Y"],
    "quotient_basis": ["1", "T", "T^2", "T^3", "U", "Y"],
    "quotient_dimension": 6,
    "relation_reduces_to_zero": True,
    "compiler_instruction": "Reduce cleared q=6 E7 numerators in this actual H92 quotient before stacking the resulting block with E8 and smooth collision blocks.",
    "boundary": "This supplies one actual E7 finite quotient and the corrected marked-chart branch. Its all-edge transition evaluator on the bounded global coefficient ambient remains the next assembly step.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6P1ACTUALE7QUOTIENT|length=6|boundary=e1+2e5|"
    "status=PASS_EXACT_Q6_P1_ACTUAL_E7_QUOTIENT_BLOCK",
    flush=True,
)
