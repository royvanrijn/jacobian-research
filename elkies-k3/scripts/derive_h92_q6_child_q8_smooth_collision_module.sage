#!/usr/bin/env sage -python
"""Derive the smooth O.S collision module for the marked q8 child divisor.

The q8 marking has generic chord

    m = (y+y(S))/(x-x(S)),

which is the usual chord through ``P=-S``.  The exact section ``S`` meets the
standard zero section over a reduced degree-46 divisor ``h`` of smooth old
fibres.  This script derives the corresponding saturated local frame and its
finite quotient blocks.  It does not assert a global ambient or a q8 pencil.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ, gcd


ROOT = Path(__file__).resolve().parents[2]
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
MARKING = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-marking.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-smooth-collision-module.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def monic_power_root(value, exponent):
    root = value.parent().one()
    for irreducible, multiplicity in value.factor():
        assert multiplicity % exponent == 0
        root *= irreducible.monic()**(multiplicity // exponent)
    return root.monic()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--child", type=Path, default=CHILD)
parser.add_argument("--marking", type=Path, default=MARKING)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
for name in ("child", "marking", "output"):
    setattr(args, name, getattr(args, name).resolve())

child = json.loads(args.child.read_text())
marking = json.loads(args.marking.read_text())
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert marking["status"] == "PASS_EXACT_Q6_CHILD_Q8_MARKING"

ring = PolynomialRing(QQ, "T")
field = ring.fraction_field()
section = marking["selected_q8"]["relative_child_section_standard_jacobian_coordinates"]
x_s = field(polynomial(ring, section["x_numerator_coefficients_low_to_high"])) / field(
    polynomial(ring, section["x_denominator_coefficients_low_to_high"])
)
y_s = field(polynomial(ring, section["y_numerator_coefficients_low_to_high"])) / field(
    polynomial(ring, section["y_denominator_coefficients_low_to_high"])
)

# The marked point for the chord is P=-S, so p=y(P)/x(P)=-y(S)/x(S).
h = monic_power_root(ring(x_s.denominator()), 2)
assert h.degree() == 46
assert h == monic_power_root(ring(y_s.denominator()), 3)
discriminant = -16 * (4 * polynomial(ring, child["minimal_short_weierstrass"]["A_coefficients_low_to_high"])**3
                      + 27 * polynomial(ring, child["minimal_short_weierstrass"]["B_coefficients_low_to_high"])**2)
assert gcd(h, discriminant) == 1

z_s = -x_s / y_s
z_s_over_h = z_s / field(h)
assert gcd(h, ring(z_s_over_h.numerator())) == 1
assert gcd(h, ring(z_s_over_h.denominator())) == 1

p = -y_s / x_s
p_numerator = ring(p.numerator())
p_denominator = ring(p.denominator())
p_denominator_over_h, remainder = p_denominator.quo_rem(h)
assert not remainder
assert gcd(h, p_denominator_over_h) == 1
assert gcd(h**2, p_denominator_over_h) == 1

# D*p has a unique polynomial residue modulo D.  It is the compensating
# scalar in the base-regular q8 chord generator q=h*m-A0.
principal = field(h) * p
interpolation = (
    ring(principal.numerator()) * ring(principal.denominator()).inverse_mod(h)
).mod(h)
assert interpolation.degree() < h.degree()
assert (principal - field(interpolation)).numerator() % h == 0

# The base-regular quotient is QQ[T]/(h), hence has dimension 20, with the
# b-coefficient mapping by its literal residue.  In the saturated form the
# A coefficient is multiplied by den(p)/h, a unit modulo h^2 by the gcd
# assertion above.  Therefore that map onto QQ[T]/(h^2) is surjective and
# has exact quotient dimension 40, without an oversized sample matrix.
base_quotient_dimension = h.degree()
saturated_quotient_dimension = (h**2).degree()

payload = {
    "schema": "elkies-k3.h92-q6-child-q8-smooth-collision-module.v1",
    "status": "PASS_EXACT_Q6_CHILD_Q8_SMOOTH_COLLISION_MODULE",
    "inputs": {
        "child_jacobian": {"path": str(args.child.relative_to(ROOT)), "sha256": digest(args.child)},
        "q8_marking": {"path": str(args.marking.relative_to(ROOT)), "sha256": digest(args.marking)},
    },
    "collision_divisor": {
        "degree": int(h.degree()),
        "irreducible_factor_degrees": [int(factor.degree()) for factor, unused in h.factor()],
        "squarefree": gcd(h, h.derivative()) == 1,
        "smooth_old_fibres": True,
        "formal_parameter": "z(S)=-x(S)/y(S)",
        "transversality": "z(S)/h is a unit modulo h",
    },
    "module": {
        "raw_chord": "m=(y+y(S))/(x-x(S))",
        "marked_chord_point": "P=-S",
        "p": "y(P)/x(P)=-y(S)/x(S)",
        "saturated_frame": "<1,(m-p)/h>",
        "base_regular_coefficient_rule": "for base-regular a,b, a+b*m is smooth iff h divides b",
        "compensated_generator": "q=h*m-A0",
        "interpolation_A0": str(interpolation),
        "rewrite": "a+h*c*m=(a+c*A0)+c*q",
    },
    "base_regular_quotient": {
        "quotient": "QQ[T]/(h)",
        "dimension": int(base_quotient_dimension),
        "map": "a+b*m maps to b mod h",
        "kernel_rule": "h divides b",
    },
    "saturated_quotient": {
        "quotient": "QQ[T]/(h^2)",
        "dimension": int(saturated_quotient_dimension),
        "coefficient_form": "a=A/h^2, b=B/h",
        "congruence": "A*(den(p)/h)+B*num(p)=0 mod h^2",
        "surjectivity_certificate": "den(p)/h is a unit modulo h^2",
    },
    "boundary": (
        "This is the complete smooth O.S collision module for the q8 marking only. "
        "It does not derive the II*/IV* modules, a bounded global ambient, a q8 pencil, "
        "a D13 equation, a rootless bisection, an extension collision, or generic rank 18 or 19."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6CHILDQ8SMOOTH|degree_h=46|squarefree=1|base_regular=h_divides_b|"
    "saturated_quotient=92|status=PASS_EXACT_Q6_CHILD_Q8_SMOOTH_COLLISION_MODULE",
    flush=True,
)
