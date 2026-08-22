#!/usr/bin/env sage -python
"""Derive the resolved E7 local target for the third H3 q=6 marked divisor.

The Weyl transport fixes the third marked class as

    43*O + (22*(-P1)-P2) + V,

where the E7 part of ``V`` is integral but large.  This script carries that
specific class into the formal E7 blow-up numbering. It proves that its
exceptional cycle is exactly 22 times the q=6 marked P1 cycle, while checking
separately that the reconstructed P2 section specializes to the smooth affine
component of the old III* fibre.

That separation matters. It identifies the 22-fold *exceptional cycle*, but
does not identify the whole third local module with a power of the P1 branch
module: the integral vertical divisor V needs its own chart trivialization and
P2 still imposes a smooth-chart condition in a global degree-44 computation.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
WEYL = ROOT / "artifacts/generated-results/elkies-k3-h3-q6-weyl-section-transport.json"
P2 = ROOT / "artifacts/generated-results/elkies-k3-h92-p2-hensel-100003-p1024.json"
ATLAS = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-e7-valuation-atlas.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-third-e7-local-target.json"
WEYL_SHA256 = "c4b7e38f0ea9fc3f748200ca9923ea3ffe5c0028c979e5f81be6507954d7c822"
P2_SHA256 = "e02e2803387d3a7f53907f548b275bb592d366f653f630f6ba8c9ef2611f3e37"
ATLAS_SHA256 = "0a408e706820a62a4c2e3290e615f5883808803565912b5b3e3e0867761f5f58"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--weyl", type=Path, default=WEYL)
parser.add_argument("--p2", type=Path, default=P2)
parser.add_argument("--atlas", type=Path, default=ATLAS)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

if args.weyl == WEYL:
    assert digest(args.weyl) == WEYL_SHA256
if args.p2 == P2:
    assert digest(args.p2) == P2_SHA256
if args.atlas == ATLAS:
    assert digest(args.atlas) == ATLAS_SHA256
weyl = json.loads(args.weyl.read_text())
p2 = json.loads(args.p2.read_text())
atlas = json.loads(args.atlas.read_text())
assert weyl["status"] == "PASS_EXACT_Q6_WEYL_SECTION_TRANSPORT"
assert p2["schema"] == "elkies-k3.h92-p2-hensel-lift.v1" and p2["complete"]
assert atlas["status"] == "PASS_FORMAL_E7_VALUATION_ATLAS"

# The marked source-to-resolution map is fixed geometrically by div(Z), as
# in derive_h92_q8_e7_local_target.sage.  It is repeated here so this target
# remains independently reproducible from its three pinned inputs.
source_to_resolved = (1, 6, 4, 3, 7, 2, 5)
cartan = matrix(ZZ, [
    [2, 0, 0, -1, 0, 0, 0],
    [0, 2, 0, 0, -1, 0, -1],
    [0, 0, 2, -1, 0, -1, -1],
    [-1, 0, -1, 2, 0, 0, 0],
    [0, -1, 0, 0, 2, 0, 0],
    [0, 0, -1, 0, 0, 2, 0],
    [0, -1, -1, 0, 0, 0, 2],
])
assert cartan * vector(ZZ, atlas["old_base_fibre_multiplicities"]) == vector(ZZ, [1] + [0] * 6)

third = weyl["third_vertical_correction"]
assert third["old_generic_restriction"] == "43*O + (22*(-P1)-P2)"
source_coefficients = vector(ZZ, third["coordinates"][:7])
assert source_coefficients == vector(ZZ, (-22, -33, -44, -66, -55, -44, -33))
resolved_coefficients = vector(ZZ, [0] * 7)
for source_index, resolved_label in enumerate(source_to_resolved):
    resolved_coefficients[resolved_label - 1] = source_coefficients[source_index]
assert resolved_coefficients == vector(ZZ, (-22, -44, -66, -44, -33, -33, -55))

# On a resolution, the exceptional degree vector of O(sum c_i E_i) is
# -Cartan*c.  It is exactly the degree-22 marked component target; this is
# derived from V, not guessed from the number 22 or the E7 diagram.
resolved_degrees = -cartan * resolved_coefficients
assert resolved_degrees == vector(ZZ, (0, 0, 0, 0, 22, 0, 0))
assert vector(ZZ, third["old_E7_component_intersections"]["vertical_correction"]) == vector(ZZ, (0, 0, 0, 0, 0, 0, 22))
assert vector(ZZ, third["old_E7_component_intersections"]["horizontal_part"]) == vector(ZZ, [0] * 7)

q6_degrees = vector(ZZ, (0, 0, 0, 0, 1, 0, 0))
q6_cycle = -cartan.inverse() * q6_degrees
assert q6_cycle == vector(QQ, (-1, -2, -3, -2, QQ(-3) / 2, QQ(-3) / 2, QQ(-5) / 2))
assert resolved_coefficients == 22 * q6_cycle
assert -cartan * (22 * q6_cycle) == resolved_degrees

# The exact P2 model is in the old local base coordinate t=Z.  Its t=0 value
# lies on the smooth locus y^2=x^3 of the singular Weierstrass fibre, because
# x(0) and y(0) are both nonzero.  Replacing it by a singular-branch module
# would therefore be a false local identification.
z0 = QQ(p2["Z"][0])
x0 = QQ(p2["X"][0]) / z0**2
y0 = QQ(p2["Y"][0]) / z0**3
assert z0 and x0 and y0 and y0**2 == x0**3

payload = {
    "schema": "elkies-k3.h92-q6-third-e7-local-target.v1",
    "status": "PASS_EXACT_Q6_THIRD_E7_LATTICE_TARGET",
    "inputs": {
        "weyl_transport": {"path": str(args.weyl.relative_to(ROOT)), "sha256": digest(args.weyl)},
        "p2": {"path": str(args.p2.relative_to(ROOT)), "sha256": digest(args.p2)},
        "valuation_atlas": {"path": str(args.atlas.relative_to(ROOT)), "sha256": digest(args.atlas)},
    },
    "generic_restriction": "43*O + (22*(-P1)-P2)",
    "source_to_resolved_component_map": list(source_to_resolved),
    "resolved_exceptional_coefficients": [int(value) for value in resolved_coefficients],
    "resolved_component_degrees": [int(value) for value in resolved_degrees],
    "exceptional_cycle": {
        "q6_marked_cycle": [str(value) for value in q6_cycle],
        "identity": "third_E7_cycle = 22*q6_marked_cycle",
        "third_cycle": [str(value) for value in resolved_coefficients],
    },
    "p2_e7_specialization": {
        "base_denominator_at_Z0_nonzero": bool(z0),
        "x_at_Z0_nonzero": bool(x0),
        "y_at_Z0_nonzero": bool(y0),
        "smooth_affine_cubic_check": bool(y0**2 == x0**3),
        "conclusion": "P2 is a smooth-affine E7 contribution, not a singular-branch contribution.",
    },
    "compiler_instruction": (
        "The displayed exceptional cycle is in formal E7 numbering only. The equality "
        "with 22 times the q6 "
        "marked cycle is an exceptional-cycle identity only, not a license to replace "
        "the third module by a P1 branch-module power. Separately impose the smooth "
        "P2 chart condition. Do not infer either condition from component labels."
    ),
    "boundary": (
        "This derives only the old E7 lattice target for the degree-44 marked divisor. "
        "It does not exhibit a coordinate transport to an H92 resolution, construct "
        "the complete high-degree ambient, a smooth P2 quotient, E8/finite conditions, "
        "a two-dimensional kernel, or child coordinates."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6THIRDE7|degrees=0,0,0,0,22,0,0|cycle=22*q6|"
    "p2=smooth_affine|status=PASS_EXACT_Q6_THIRD_E7_LATTICE_TARGET",
    flush=True,
)
