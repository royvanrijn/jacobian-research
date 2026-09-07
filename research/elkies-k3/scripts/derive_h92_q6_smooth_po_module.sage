#!/usr/bin/env sage -python
"""Derive the smooth ``P1.O`` collision module for the first H3 neighbour.

The marked section has affine denominator

    x(P1)=N_x/(c_x*u^4*h(u)^2),   y(P1)=N_y/(c_y*u^6*h(u)^3),

where ``h`` is squarefree of degree four.  At a root of ``h`` the formal
parameter ``z_P=-x(P1)/y(P1)`` is a unit times ``h``.  Thus ``P1`` (and
``-P1``) meets ``O`` transversely on each of these four *smooth* fibres.
The raw chord

    m=(y-y(P1))/(x-x(P1))

has one unpaired vertical pole there.  Put ``p=y(P1)/x(P1)=-1/z_P``.  The
complete local lattice is the saturated frame

    <1, (m-p)/h>.

For a *base-regular* coefficient pair ``a+b*m``, this restricts to ``h | b``:
no scalar ``a`` can cancel the pole of ``b*p`` otherwise.  Writing ``b=h*c``
then gives

    a+b*m = (a+c*A) + c*(h*m-A),

where ``A`` is the unique degree-<4 interpolation residue of ``h*y(P1)/x(P1)``
modulo ``h``.  More generally, for the global coefficients
``a=A/h^2, b=B/h``, saturated membership is the exact congruence

    A*(den(p)/h) + B*num(p) = 0 mod h^2.

This script certifies both presentations and exports their exact quotient
condition blocks.  It deliberately does not supply the global degree bounds
or any condition at the III* and II* fibres.
"""

import argparse
import hashlib
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import PolynomialRing, QQ, gcd, vector


ROOT = Path(__file__).resolve().parents[2]
SECTION = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
SECTION_SHA256 = "c323bf6346bb239934a5a2d8b1a3f4067e70e993d2e4eb32aaa30f469fca6397"
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-smooth-po-module.json"

exec(compile(CORE.read_text(), str(CORE), "exec"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def smooth_po_chord_condition(ambient, h):
    """Return the base-regular sublattice block imposing ``h | b``.

    ``ambient`` consists of ``('a', exponent)`` and ``('b', exponent)``
    labels.  Negative exponents are legitimate only after a caller has proved
    its Laurent ambient; here they reduce uniquely because ``h(0) != 0``.
    """
    ring = h.parent()
    u = ring.gen()
    quotient_basis = tuple(u**degree for degree in range(h.degree()))

    def residue(entry):
        kind, exponent = entry
        exponent = int(exponent)
        if kind == "a":
            return vector(QQ, len(quotient_basis))
        if kind != "b":
            raise ValueError("unknown smooth-collision coefficient kind {}".format(kind))
        if exponent >= 0:
            remainder = (u**exponent).mod(h)
        else:
            remainder = (u ** (-exponent)).inverse_mod(h)
        return vector(
            QQ,
            [remainder.monomial_coefficient(monomial) for monomial in quotient_basis],
        )

    return quotient_condition(
        "H92 smooth P1.O base-regular chord collisions",
        ambient,
        residue,
        quotient_basis,
        "P1 meets O transversely at the four squarefree roots of h; h | b",
    )


def smooth_po_saturated_condition(ambient, h, p_numerator, p_denominator_over_h):
    """Return the saturated ``a=A/h^2, b=B/h`` collision block.

    An ``('A', i)`` label denotes ``A=u^i`` and an ``('B', i)`` label
    denotes ``B=u^i``.  The local coefficient of ``1`` in the frame
    ``<1,(m-p)/h>`` is ``a+b*p``; clearing its unit denominator gives the
    displayed congruence modulo ``h^2``.
    """
    ring = h.parent()
    u = ring.gen()
    modulus = h**2
    quotient_basis = tuple(u**degree for degree in range(modulus.degree()))

    def residue(entry):
        kind, exponent = entry
        exponent = int(exponent)
        if kind == "A":
            value = u**exponent * p_denominator_over_h
        elif kind == "B":
            value = u**exponent * p_numerator
        else:
            raise ValueError("unknown saturated smooth coefficient kind {}".format(kind))
        remainder = value.mod(modulus)
        return vector(
            QQ,
            [remainder.monomial_coefficient(monomial) for monomial in quotient_basis],
        )

    return quotient_condition(
        "H92 saturated smooth P1.O chord collisions",
        ambient,
        residue,
        quotient_basis,
        "a=A/h^2, b=B/h in <1,(m-y_P/x_P)/h>",
    )


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

assert digest(SECTION) == SECTION_SHA256
section = json.loads(SECTION.read_text())
assert section["status"] == "PASS_EXACT_H92_P1"

ring = PolynomialRing(QQ, "u")
u = ring.gen()
field = ring.fraction_field()
h = polynomial(ring, section["structured_denominator"]["Z4_coefficients"])
x_num = polynomial(ring, section["x_entrance_base"]["numerator_coefficients"])
x_den = polynomial(ring, section["x_entrance_base"]["denominator_coefficients"])
y_num = polynomial(ring, section["y_entrance_base"]["numerator_coefficients"])
y_den = polynomial(ring, section["y_entrance_base"]["denominator_coefficients"])

assert h.degree() == 4
assert h(0) == 1
assert gcd(h, h.derivative()) == 1
z = x_den.sqrt()
assert z**2 == x_den
scale = z / (u**2 * h)
assert scale in QQ and scale
assert y_den == (y_den.leading_coefficient() / z.leading_coefficient()**3) * z**3

# z_P=-x_P/y_P is h times a unit at every root of h.  These gcd assertions
# are the exact transverse-intersection certificate in the formal O-chart.
x_p = field(x_num) / field(x_den)
y_p = field(y_num) / field(y_den)
z_p_over_h = -x_p / y_p / field(h)
assert z_p_over_h.numerator().mod(h)
assert z_p_over_h.denominator().mod(h)
assert gcd(h, ring(z_p_over_h.numerator())) == 1
assert gcd(h, ring(z_p_over_h.denominator())) == 1

principal = field(h) * y_p / x_p
assert gcd(ring(principal.denominator()), h) == 1
A = (ring(principal.numerator()) * ring(principal.denominator()).inverse_mod(h)).mod(h)
assert A.degree() < h.degree()
assert (principal - field(A)).numerator() % h == 0

p = y_p / x_p
p_numerator = ring(p.numerator())
p_denominator = ring(p.denominator())
p_denominator_over_h, p_remainder = p_denominator.quo_rem(h)
assert not p_remainder
assert gcd(h, p_denominator_over_h) == 1

# A regression block intentionally includes a Laurent coefficient.  It checks
# that reduction modulo h, not a root-by-root numerical approximation, is used.
regression_ambient = (("a", -2), ("a", 0), ("b", -1), ("b", 0), ("b", 1), ("b", 2), ("b", 3))
block = smooth_po_chord_condition(regression_ambient, h)
assert block["matrix"].rank() == 4
assert block["matrix"].right_kernel().dimension() == 3

saturated_ambient = tuple(
    [("A", exponent) for exponent in range(1, 9)]
    + [("B", exponent) for exponent in range(3, 5)]
)
saturated_block = smooth_po_saturated_condition(
    saturated_ambient, h, p_numerator, p_denominator_over_h
)
assert saturated_block["matrix"].rank() == 8
assert saturated_block["matrix"].right_kernel().dimension() == 2

payload = {
    "schema": "elkies-k3.h92-q6-smooth-po-module.v1",
    "status": "PASS_EXACT_SMOOTH_PO_CHORD_MODULE",
    "inputs": {
        "marked_section": {"path": str(SECTION.relative_to(ROOT)), "sha256": SECTION_SHA256},
    },
    "collision_divisor": {
        "h": str(h),
        "degree": int(h.degree()),
        "squarefree": True,
        "formal_parameter": "z_P=-x(P1)/y(P1)",
        "transversality": "z_P/h is a unit modulo h",
    },
    "module": {
        "raw_chord": "m=(y-y(P1))/(x-x(P1))",
        "saturated_frame": "<1,(m-y(P1)/x(P1))/h>",
        "base_regular_coefficient_rule": "for base-regular a,b, a+b*m is smooth iff h divides b",
        "compensated_generator": "q=h*m-A",
        "interpolation_A": str(A),
        "rewrite": "a+h*c*m=(a+c*A)+c*q",
    },
    "regression_block": {
        "ambient": [[kind, exponent] for kind, exponent in regression_ambient],
        "quotient_basis": list(block["quotient_basis"]),
        "matrix": [[str(value) for value in row] for row in block["matrix"].rows()],
        "rank": int(block["matrix"].rank()),
        "kernel_dimension": int(block["matrix"].right_kernel().dimension()),
    },
    "saturated_regression_block": {
        "ambient": [[kind, exponent] for kind, exponent in saturated_ambient],
        "coefficient_form": "a=A/h^2, b=B/h",
        "congruence": "A*(den(y(P1)/x(P1))/h)+B*num(y(P1)/x(P1))=0 mod h^2",
        "quotient_basis": list(saturated_block["quotient_basis"]),
        "matrix": [[str(value) for value in row] for row in saturated_block["matrix"].rows()],
        "rank": int(saturated_block["matrix"].rank()),
        "kernel_dimension": int(saturated_block["matrix"].right_kernel().dimension()),
    },
    "boundary": (
        "This is the complete smooth P1.O contribution only. It does not prove "
        "a bounded global Laurent ambient, the E7 representative conversion, "
        "or an h0(D)=2 calculation."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6SMOOTHPO|degree_h=4|squarefree=1|base_regular=h_divides_b|"
    "saturated_rows=8|saturated_kernel=2|status=PASS_EXACT_SMOOTH_PO_CHORD_MODULE",
    flush=True,
)
