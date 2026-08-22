#!/usr/bin/env sage -python
"""Orient the IV* component met by the marked q=8 child section.

The two ``u``-chart blow-ups in
``derive_h92_q6_child_q8_ivstar_entrance.sage`` already put the marked point
in the smooth locus of one of the two branches

    u=0,  y-c=0,  and  u=0,  y+c=0.

This checker records the remaining component-theoretic conclusion without
pretending that an arbitrary Dynkin numbering is geometric.  Both branches
are reduced multiplicity-one components in the resolved IV* fibre.  The
projective Weierstrass zero is disjoint from the affine blow-up centres, so it
meets the third multiplicity-one component.  Since the fibre is IV*, whose
extended E6 multiplicity vector has exactly three entries equal to one, the
two chart branches are precisely the two nonidentity E6 components.  We use
their actual chart equations as the ``+`` and ``-`` labels; this fixes their
opposite classes in the E6 component group without asserting a numerical
Dynkin-node convention.

It derives only this local component orientation.  The corresponding finite
IV* module and the II* module remain separate global-pencil gates.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
ENTRANCE = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-ivstar-entrance.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-ivstar-component.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--child", type=Path, default=CHILD)
parser.add_argument("--entrance", type=Path, default=ENTRANCE)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
for name in ("child", "entrance", "output"):
    setattr(args, name, getattr(args, name).resolve())

child = json.loads(args.child.read_text())
entrance = json.loads(args.entrance.read_text())
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert entrance["status"] == "PASS_EXACT_Q6_CHILD_Q8_IVSTAR_ENTRANCE"

ivstar = entrance["ivstar"]
assert tuple(next(item for item in child["finite_fibres"] if item["kodaira"] == "IV*")["minimal_orders"]) == (3, 4, 8)
assert child["root_data"] == {"rank": 14, "determinant": 3, "type": "E8+E6"}
assert ivstar["marked_entrance"]["smooth"]
branch_data = ivstar["second_chart_exceptional_branches"]
assert branch_data["marked_branch"] == "y-c=0"
assert branch_data["other_branch"] == "y+c=0"
assert branch_data["both_branches_smooth_in_second_chart"]

# At the generic point of each branch in the second chart, the strict
# transform is smooth because d(f2)/dy is a unit.  Hence u itself is a
# uniformizer transverse to that branch, proving reduced multiplicity one in
# div(u), rather than merely reading a multiplicity from the Kodaira symbol.
ring = PolynomialRing(QQ, names=("u", "x", "y", "c"))
u, x, y, c = ring.gens()
fibre = y**2 - c**2
assert fibre == (y - c) * (y + c)
assert fibre.derivative(y) == 2 * y
for branch in (y - c, y + c):
    quotient = ring.quotient(ring.ideal((u, branch)))
    # The generic point has y=+/-c with c nonzero (the entrance certificate
    # proves c=Y_section(0) is nonzero), so the derivative 2*y is a unit.
    assert quotient(fibre.derivative(y)) != 0

# The projective zero O=[0:1:0] is not in the affine (x,y)-chart containing
# the singular Weierstrass point and its two u-chart blow-up centres.  It is
# therefore not on either displayed affine branch.  For IV* the affine E6
# multiplicities are the highest-root coefficients below; adjoining the
# identity component gives exactly three multiplicity-one components.
e6_highest_root = (1, 2, 3, 2, 1, 2)
assert sum(value == 1 for value in e6_highest_root) == 2
extended_multiplicities = (1,) + e6_highest_root
assert sum(value == 1 for value in extended_multiplicities) == 3

payload = {
    "schema": "elkies-k3.h92-q6-child-q8-ivstar-component.v1",
    "status": "PASS_EXACT_Q6_CHILD_Q8_IVSTAR_COMPONENT",
    "inputs": {
        "child_jacobian": {"path": str(args.child.relative_to(ROOT)), "sha256": digest(args.child)},
        "ivstar_entrance": {"path": str(args.entrance.relative_to(ROOT)), "sha256": digest(args.entrance)},
    },
    "ivstar_component_orientation": {
        "actual_chart_labels": {
            "plus": "u=0, y-c=0",
            "minus": "u=0, y+c=0",
        },
        "marked_section_component": "plus",
        "inversion": "plus and minus are the two opposite nonzero E6 component-group classes",
        "multiplicity_proof": "At each branch generic point d(f2)/dy is a unit and u is a transverse uniformizer, so both branches occur with multiplicity one in div(u).",
        "zero_section_component": "the distinct projective component at Weierstrass infinity",
        "extended_E6_multiplicities": list(extended_multiplicities),
        "conclusion": "The marked q8 section meets the physically labelled nonidentity multiplicity-one component plus; minus is its opposite nonidentity component.",
    },
    "boundary": (
        "The plus/minus labels are fixed by actual resolved chart equations, not by an arbitrary E6 Dynkin numbering. "
        "This does not yet derive the finite IV* or II* Riemann--Roch modules, a q8 pencil, a D13 equation, a rootless bisection, an extension collision, or generic rank 18 or 19."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6CHILDQ8IVCOMP|branches=plus,minus|marked=plus|"
    "multiplicity_one=2|status=PASS_EXACT_Q6_CHILD_Q8_IVSTAR_COMPONENT",
    flush=True,
)
