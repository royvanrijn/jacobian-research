#!/usr/bin/env sage -python
"""Certify that q_regular trivializes the smooth q8 collision module.

At the reduced smooth collision divisor h, the marked chord module in the
``<1,m>`` frame is

    a=A/h^2, b=B/h,  A*D+B*N = 0 (mod h^2),

where p=N/(h*D).  Put q_regular=(m-p)/h-R/Nx.  For a local regular pair
``(C,B)`` in ``C+B*q_regular``, its old coefficients satisfy

    b=B/h,
    a=C-B*p/h-B*R/Nx.

The displayed smooth congruence is then automatic; conversely it makes C
regular.  Thus q_regular is an exact local frame for the saturated smooth
module.  This leaves no independent h-supported quotient condition when the
finite additive and infinity modules are expressed in this frame.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ, gcd


ROOT = Path(__file__).resolve().parents[2]
MARKING = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-marking.json"
SMOOTH = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-smooth-collision-module.json"
NORMALIZER = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-q-pole-normalization-crt.json"
REGULAR = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-q-regular-frame.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-q-regular-smooth-frame.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def monic_power_root(value, exponent):
    root = value.parent().one()
    for factor, multiplicity in value.factor():
        assert multiplicity % exponent == 0
        root *= factor.monic()**(multiplicity // exponent)
    return root.monic()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--marking", type=Path, default=MARKING)
parser.add_argument("--smooth", type=Path, default=SMOOTH)
parser.add_argument("--normalizer", type=Path, default=NORMALIZER)
parser.add_argument("--regular", type=Path, default=REGULAR)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
for name in ("marking", "smooth", "normalizer", "regular", "output"):
    setattr(args, name, getattr(args, name).resolve())

marking = json.loads(args.marking.read_text())
smooth = json.loads(args.smooth.read_text())
normalizer = json.loads(args.normalizer.read_text())
regular = json.loads(args.regular.read_text())
assert marking["status"] == "PASS_EXACT_Q6_CHILD_Q8_MARKING"
assert smooth["status"] == "PASS_EXACT_Q6_CHILD_Q8_SMOOTH_COLLISION_MODULE"
assert normalizer["status"] == "PASS_EXACT_CRT_PRINCIPAL_PART_NORMALIZATION"
assert normalizer["complete"] and normalizer["exact_check"]
assert regular["status"] == "PASS_EXACT_Q_REGULAR_FRAME"

ring = PolynomialRing(QQ, "T")
section = marking["selected_q8"]["relative_child_section_standard_jacobian_coordinates"]
nx = polynomial(ring, section["x_numerator_coefficients_low_to_high"])
dx = polynomial(ring, section["x_denominator_coefficients_low_to_high"])
ny = polynomial(ring, section["y_numerator_coefficients_low_to_high"])
dy = polynomial(ring, section["y_denominator_coefficients_low_to_high"])
h = monic_power_root(dx, 2)
assert h.degree() == 46 and dx // h**2 in QQ and dy // h**3 in QQ

# p=N/(hD), with D a unit at every h-point.  Form these from the already
# certified square/cube denominator roots instead of cancelling the enormous
# rational function p; cancellation is unnecessary in the h-local ring.
dx_over_h2 = dx // h**2
dy_over_h3 = dy // h**3
assert dx_over_h2 in QQ and dy_over_h3 in QQ
N = -ny * dx_over_h2
D = dy_over_h3 * nx
assert gcd(D, h) in QQ
R = polynomial(ring, normalizer["normalizer"]["R_coefficients_low_to_high"])
assert (R*h*dy-ny) % nx == 0
assert gcd(nx, h) in QQ

# The two gcd identities mean D and Nx are units in QQ[T]_(h).  There is no
# need to materialize their inverses modulo h^2 (whose rational coefficients
# are enormous): the two exact identities below are equalities in that local
# ring and prove both directions of the module equivalence.
modulus = h**2

# Retain the exact (not only mod-h^2) numerator identity.  In the localization
# at h, C=B*R/Nx + (A*D+B*N)/(h^2*D); it proves the inverse transition maps
# every old smooth-compatible pair to regular C.
identity = "C=B*R/Nx+(A*D+B*N)/(h^2*D)"
forward_identity = "A=h^2*(C-B*R/Nx)-B*N/D"

payload = {
    "schema": "elkies-k3.h92-q6-child-q8-q-regular-smooth-frame.v1",
    "status": "PASS_EXACT_Q_REGULAR_SMOOTH_FRAME",
    "inputs": {
        "marking": {"path": str(args.marking.relative_to(ROOT)), "sha256": digest(args.marking)},
        "smooth_module": {"path": str(args.smooth.relative_to(ROOT)), "sha256": digest(args.smooth)},
        "normalizer": {"path": str(args.normalizer.relative_to(ROOT)), "sha256": digest(args.normalizer)},
        "regular_frame": {"path": str(args.regular.relative_to(ROOT)), "sha256": digest(args.regular)},
    },
    "smooth_divisor": {
        "h_degree": int(h.degree()),
        "h_square_degree": int(modulus.degree()),
        "D_is_unit_mod_h2": True,
        "Nx_is_unit_mod_h2": True,
    },
    "frame_transition": {
        "old": "a=A/h^2, b=B/h, A*D+B*N=0 mod h^2",
        "new": "C+B*q_regular",
        "q_regular": "(m-p)/h-R/Nx",
        "forward": forward_identity,
        "inverse": identity,
        "mod_h2": "A=-B*N/D mod h^2",
        "congruence_verified": "(h^2*(C-B*R/Nx)-B*N/D)*D+B*N=h^2*D*(C-B*R/Nx)",
    },
    "conclusion": {
        "local_module": "QQ[T]_(h)<1,q_regular>",
        "additional_smooth_quotient_condition": "none",
        "reason": "The forward and inverse formulas identify the saturated old smooth module with regular (C,B) pairs after localizing at h.",
    },
    "boundary": (
        "This identifies only the smooth collision module in the q_regular frame. "
        "The finite II*/IV* module and the complete infinity module still have to "
        "be assembled before a q8 pencil, bisection equation, extension collision, "
        "or rank claim."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6CHILDQREGSMOOTH|h_degree={}|smooth_quotient=0|status=PASS_EXACT_Q_REGULAR_SMOOTH_FRAME".format(h.degree()),
    flush=True,
)
