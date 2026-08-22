#!/usr/bin/env sage -python
"""Transport the exact finite q8 coefficient module to q_regular=q-R/Nx.

The finite q-frame certificate gives the diagonal module

    B in (f_IV),   C_old in (f_II^2*f_IV^3)

for C_old+B*q.  Since q_regular=q-R/Nx, a pair C+B*q_regular has
C_old=C-B*R/Nx.  At the additive CRT modulus Nx is a unit.  This script
therefore exports the exact q_regular module

    < (f_IV, lift(f_IV*R/Nx)), (0, f_II^2*f_IV^3) >,

in (B,C) order.  It is the local finite condition in a frame that is already
globally base-regular.  No infinity condition or global pencil is asserted.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ, matrix


ROOT = Path(__file__).resolve().parents[2]
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
FINITE = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-finite-q-module-qq.json"
NORMALIZER = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-q-pole-normalization-crt.json"
REGULAR = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-q-regular-frame.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-q-regular-finite-module-qq.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--finite", type=Path, default=FINITE)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

child = json.loads(CHILD.read_text())
finite_path = args.finite.resolve()
finite = json.loads(finite_path.read_text())
normalizer_record = json.loads(NORMALIZER.read_text())
regular = json.loads(REGULAR.read_text())
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert finite["status"] == "PASS_EXACT_FINITE_Q_CONDITIONS"
assert normalizer_record["status"] == "PASS_EXACT_CRT_PRINCIPAL_PART_NORMALIZATION"
assert normalizer_record["complete"] and normalizer_record["exact_check"]
assert regular["status"] == "PASS_EXACT_Q_REGULAR_FRAME"

ring = PolynomialRing(QQ, "T")
ii = ring(next(item for item in child["finite_fibres"] if item["kodaira"] == "II*")["factor"])
iv = ring(next(item for item in child["finite_fibres"] if item["kodaira"] == "IV*")["factor"])
nef_ivstar = finite["module"]["exact_basis"][0][0] == "1"
assert finite["module"]["smith_degrees"] == ([0, 4] if nef_ivstar else [1, 5])
modulus = ii**2 * iv**(2 if nef_ivstar else 3)
assert modulus.degree() == (4 if nef_ivstar else 5) and ii.gcd(iv) in QQ

# The section data are deliberately read through the normalizer certificate's
# pinned input, while Nx itself is rebuilt from the marking it references.
marking_path = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-marking.json"
marking = json.loads(marking_path.read_text())
section = marking["selected_q8"]["relative_child_section_standard_jacobian_coordinates"]
nx = polynomial(ring, section["x_numerator_coefficients_low_to_high"])
normalizer = polynomial(ring, normalizer_record["normalizer"]["R_coefficients_low_to_high"])
assert nx.degree() == 96 and normalizer.degree() == 95 and nx.gcd(modulus) in QQ

r_modulus = (normalizer * nx.inverse_mod(modulus)).mod(modulus)
first_b = ring.one() if nef_ivstar else iv
first_c_lift = (first_b * r_modulus).mod(modulus)
assert (first_c_lift - first_b*r_modulus) % modulus == 0

# In (B,C) order, C_old=C-B*R/Nx.  The two displayed columns generate
# exactly B in (iv) and C_old in (modulus).  Their determinant has the old
# finite index degree 1+5=6.
module_matrix = matrix(ring, [[first_b, 0], [first_c_lift, modulus]])
assert module_matrix.det() == first_b*modulus
assert (first_c_lift - first_b*r_modulus) % modulus == 0

payload = {
    "schema": "elkies-k3.h92-q6-child-q8-q-regular-finite-module-qq.v1",
    "status": "PASS_EXACT_Q_REGULAR_FINITE_MODULE",
    "inputs": {
        "child": digest(CHILD),
        "finite_q_module": digest(finite_path),
        "normalizer": digest(NORMALIZER),
        "regular_frame": digest(REGULAR),
    },
    "frame_change": {
        "old": "C_old+B*q",
        "new": "C+B*q_regular",
        "q_regular": "q-R/Nx",
        "relation": "C_old=C-B*R/Nx",
    },
    "additive_CRT": {
        "modulus": "f_II*^2*f_IV*^2" if nef_ivstar else "f_II*^2*f_IV*^3",
        "degree": int(modulus.degree()),
        "Nx_is_a_unit": True,
        "R_over_Nx_modulus_coefficients_low_to_high": [str(value) for value in r_modulus.list()],
    },
    "module": {
        "coordinate_order": ["B", "C"],
        "basis": [
            ["1", "lift(R/Nx)"] if nef_ivstar else ["f_IV*", "lift(f_IV*R/Nx)"],
            ["0", "f_II*^2*f_IV*^2"] if nef_ivstar else ["0", "f_II*^2*f_IV*^3"],
        ],
        "first_C_lift_coefficients_low_to_high": [str(value) for value in first_c_lift.list()],
        "determinant": "f_II*^2*f_IV*^2" if nef_ivstar else "f_IV* * f_II*^2*f_IV*^3",
        "finite_codimension": int(modulus.degree() + first_b.degree()),
        "smith_degrees": [0, int(modulus.degree())] if nef_ivstar else [1, 5],
    },
    "boundary": (
        "This transports the complete finite additive module into the globally "
        "base-regular q frame. It does not impose the remaining smooth or "
        "infinity conditions, construct a q8 pencil/D13/rootless equation, "
        "or supply bisections, extension collisions, or rank."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q6CHILDQREGFINITE|CRT_degree={}|codimension={}|"
    "status=PASS_EXACT_Q_REGULAR_FINITE_MODULE".format(
        modulus.degree(), payload["module"]["finite_codimension"]
    ),
    flush=True,
)
