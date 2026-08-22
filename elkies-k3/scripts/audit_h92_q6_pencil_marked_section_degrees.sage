#!/usr/bin/env sage -python
"""Audit the first H92 q=6 pencil on its three marked old sections.

The first exact q=6 pencil has parameter ``T=(a1+b1*m)/(a0+b0*m)`` where
``m=(y-y(P1))/(x-x(P1))``.  This script prevents a tempting but false shortcut
in the bisection search: the old zero and ``-P1`` are sections of this pencil,
while ``P1`` maps to the q=6 base with degree 22.  Thus none of these three
marked old sections is a q=6 bisection.

This is only a marked-section restriction audit.  It neither enumerates
rootless bisections nor constructs a quadratic extension.
"""

import argparse
import hashlib
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
RR = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-global-rr.json"
SECTION = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-marked-section-degrees.json"
SECTION_SHA256 = "c323bf6346bb239934a5a2d8b1a3f4067e70e993d2e4eb32aaa30f469fca6397"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--rr", type=Path, default=RR)
parser.add_argument("--section", type=Path, default=SECTION)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
for name in ("rr", "section", "output"):
    setattr(args, name, getattr(args, name).resolve())

rr = json.loads(args.rr.read_text())
assert rr["status"] == "PASS_EXACT_GLOBAL_RR_KERNEL"
assert rr["kernel"]["dimension"] == 2
assert digest(args.section) == SECTION_SHA256
section = json.loads(args.section.read_text())

u_ring = PolynomialRing(QQ, "u")
u = u_ring.gen()
u_field = u_ring.fraction_field()
h = polynomial(u_ring, section["structured_denominator"]["Z4_coefficients"])
x_p = u_field(polynomial(
    u_ring, section["x_entrance_base"]["numerator_coefficients"]
)) / u_field(polynomial(
    u_ring, section["x_entrance_base"]["denominator_coefficients"]
))
y_p = u_field(polynomial(
    u_ring, section["y_entrance_base"]["numerator_coefficients"]
)) / u_field(polynomial(
    u_ring, section["y_entrance_base"]["denominator_coefficients"]
))


def section_function(entry):
    A = polynomial(u_ring, entry["A_coefficients_low_to_high"])
    B = polynomial(u_ring, entry["B_coefficients_low_to_high"])
    return u_field(A) / u_field(h**2), u_field(B) / u_field(h)


(a0, b0), (a1, b1) = tuple(section_function(entry) for entry in rr["kernel"]["sections"])
assert a0 * b1 - a1 * b0

# Both the old zero and -P1 are poles of m.  Taking the leading m coefficient
# in the pencil ratio consequently restricts T to b1/b0 on either section.
pole_restriction = b1 / b0
assert pole_restriction == u

# At P1 the chord has its tangent value.  This gives the complete restriction
# of T to P1 as a rational function in the old entrance parameter u.
anchor = SourceFileLoader("h92_q6_marked_degree_anchor", str(ANCHOR)).load_module()
r92, s92 = anchor.EXPECTED_H92
_, h92_formulas = anchor.parse_h92(H92)
A1, A, unused_B1, unused_B, unused_B2 = (
    QQ(value(r92, s92)) for value in h92_formulas
)
old_a = A1 / u**3 + A / u**4
m_at_p1 = (3 * x_p**2 + old_a) / (2 * y_p)
t_at_p1 = (a1 + b1 * m_at_p1) / (a0 + b0 * m_at_p1)
assert t_at_p1.numerator().degree() == 22
assert t_at_p1.denominator().degree() == 22
assert t_at_p1.derivative().numerator().degree() == 42

payload = {
    "schema": "elkies-k3.h92-q6-marked-section-degrees.v1",
    "status": "PASS_EXACT_Q6_MARKED_SECTION_RESTRICTIONS",
    "inputs": {
        "global_rr": {"path": str(args.rr.relative_to(ROOT)), "sha256": digest(args.rr)},
        "p1_section": {"path": str(args.section.relative_to(ROOT)), "sha256": digest(args.section)},
    },
    "pencil_parameter": "T=(a1+b1*m)/(a0+b0*m), m=(y-y(P1))/(x-x(P1))",
    "restrictions": {
        "old_zero": {
            "reason": "m has a pole at the old zero, so T restricts to b1/b0",
            "parameter": str(pole_restriction),
            "degree": 1,
        },
        "minus_P1": {
            "reason": "m has a pole at -P1, so T restricts to b1/b0",
            "parameter": str(pole_restriction),
            "degree": 1,
        },
        "P1": {
            "reason": "m(P1)=(3*x(P1)^2+a)/(2*y(P1)) is the tangent chord slope",
            "numerator_degree": int(t_at_p1.numerator().degree()),
            "denominator_degree": int(t_at_p1.denominator().degree()),
            "map_degree": 22,
            "derivative_numerator_degree": int(t_at_p1.derivative().numerator().degree()),
        },
    },
    "conclusion": "The old zero, -P1, and P1 have q=6-base degrees 1, 1, and 22 respectively; none is a bisection.",
    "boundary": "This rules out only these three marked old sections as immediate q=6 bisections. It does not enumerate all multisections, construct a rootless bisection equation, or compute a quadratic extension or collision.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6MARKEDDEGREES|O=1|-P1=1|P1=22|status=PASS_EXACT_Q6_MARKED_SECTION_RESTRICTIONS",
    flush=True,
)
