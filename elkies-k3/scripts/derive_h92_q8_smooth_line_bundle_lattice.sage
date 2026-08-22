#!/usr/bin/env sage -python
"""Identify the q=8 smooth-collision line-bundle lattice from its divisor.

The source-nef q=8 class restricts to 9*O+9*(-P1) on the old generic
fibre.  Its remaining divisor is supported on the two additive fibres,
together with -11 times a chosen fibre at infinity.  Thus at each of the four
smooth transverse O.(-P1) collisions it has no vertical modification.

The local q=8 line-bundle lattice is consequently the regular algebra frame
in q=(m-y(P1)/x(P1))/h and X=h^2*x, namely

    <1,q,...,q^9,X,X*q,...,X*q^7>.

This turns the existing h-principal-part map into the actual smooth condition:
all negative h-principal parts must vanish.  It does not handle the E7/E8
resolved conditions or prove a global two-dimensional kernel.
"""

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GENERIC = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-generic-rr-ambient.json"
FRAME = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-smooth-collision-frame.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-smooth-line-bundle-lattice.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--generic", type=Path, default=GENERIC)
parser.add_argument("--frame", type=Path, default=FRAME)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

generic = json.loads(args.generic.read_text())
frame = json.loads(args.frame.read_text())
assert generic["status"] == "PASS_EXACT_Q8_GENERIC_RR_AMBIENT"
assert frame["status"] == "PASS_EXACT_Q8_SMOOTH_COLLISION_FRAME"
assert generic["generic_fibre_divisor"] == "9*O + 9*(-P1)"
vertical = generic["source_q8_lattice_selection"]["vertical_difference_D_minus_9O_minus_9minusP1"]
assert vertical == [-11, 0, 2, 3, 4, 6, 5, 5, 6, -4, -5, -7, -10, -8, -6, -4, -2, 0, 0]
assert len(frame["regular_degree_18_frame"]) == 18
assert frame["coordinates"]["q"] == "(m-y(P1)/x(P1))/h"
assert frame["coordinates"]["X"] == "h^2*x"

payload = {
    "schema": "elkies-k3.h92-q8-smooth-line-bundle-lattice.v1",
    "status": "PASS_EXACT_Q8_SMOOTH_LINE_BUNDLE_LATTICE",
    "inputs": {
        "generic_q8_divisor": {"path": str(args.generic.relative_to(ROOT)), "sha256": digest(args.generic)},
        "smooth_collision_frame": {"path": str(args.frame.relative_to(ROOT)), "sha256": digest(args.frame)},
    },
    "local_divisor": {
        "horizontal_part": "9*O+9*(-P1)",
        "vertical_difference": vertical,
        "support_argument": (
            "The non-fibre entries are E7/E8 component terms and are supported "
            "only at the additive fibres. The remaining -11*F may be represented "
            "by a fibre away from any chosen smooth collision. Hence the q8 divisor "
            "has no vertical modification in each h-adic collision neighbourhood."
        ),
    },
    "line_bundle_lattice": {
        "coordinates": frame["coordinates"],
        "regular_basis": frame["regular_degree_18_frame"],
        "condition": "all negative h-principal parts in the q,X regular frame vanish",
        "principal_part_map_use": (
            "The complete h^-15 principal-part map is the q8 smooth condition "
            "after setting its target coordinates to zero."
        ),
    },
    "boundary": (
        "This identifies the smooth q8 line-bundle lattice only. The E7/E8 "
        "resolved quotient conditions, their compatible common kernel, the D13 "
        "equation, rootless equation, bisection covers, collisions, and rank "
        "claims remain unproved."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("H92Q8SMOOTHLATTICE|basis=18|principal_parts=must_vanish|status=PASS_EXACT_Q8_SMOOTH_LINE_BUNDLE_LATTICE", flush=True)
