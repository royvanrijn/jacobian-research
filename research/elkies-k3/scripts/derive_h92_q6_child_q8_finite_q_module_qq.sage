#!/usr/bin/env sage -python
"""Compute the exact finite q8-child coefficient module in the q=(m-p)/h frame.

The saturated smooth frame is ``<1,q>`` with
``q=(m-p)/h``, where ``p=-y(S)/x(S)`` and ``h`` is the reduced O.S collision
divisor.  At the two additive fibres, h is a unit.  This removes the costly
global interpolation of A0 modulo h^2: a finite coefficient pair is simply
``C+B*q``.  The complete II* and oriented IV* quotient rows have exact rank
six on the ten-dimensional CRT residue space.  In fact their kernel is the
diagonal module ``< (f_IV,0), (0,f_II^2*f_IV^3) >``.  No condition at
infinity is included, so this is not yet a q8 pencil.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ, matrix


ROOT = Path(__file__).resolve().parents[2]
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
MARKING = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-marking.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-finite-q-module-qq.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def monic_power_root(value, exponent):
    root = value.parent().one()
    for irreducible, multiplicity in value.factor():
        assert multiplicity % exponent == 0
        root *= irreducible.monic() ** (multiplicity // exponent)
    return root.monic()


def reduce_rational(value, modulus):
    ring = modulus.parent()
    numerator, denominator = ring(value.numerator()), ring(value.denominator())
    assert denominator.gcd(modulus) in QQ
    return (numerator * denominator.inverse_mod(modulus)).mod(modulus)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--nef-ivstar", action="store_true",
    help="use the Weyl-nef IV* ideal (u^2,X,Y) instead of the dominant D13 arm ideal",
)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

child = json.loads(CHILD.read_text())
marking = json.loads(MARKING.read_text())
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert marking["status"] == "PASS_EXACT_Q6_CHILD_Q8_MARKING"

ring = PolynomialRing(QQ, "T")
T = ring.gen()
field = ring.fraction_field()
section = marking["selected_q8"]["relative_child_section_standard_jacobian_coordinates"]
sx = field(polynomial(ring, section["x_numerator_coefficients_low_to_high"])) / field(
    polynomial(ring, section["x_denominator_coefficients_low_to_high"])
)
sy = field(polynomial(ring, section["y_numerator_coefficients_low_to_high"])) / field(
    polynomial(ring, section["y_denominator_coefficients_low_to_high"])
)
h = monic_power_root(ring(sx.denominator()), 2)
assert h.degree() == 46
assert h == monic_power_root(ring(sy.denominator()), 3)


def local_data(kodaira):
    fibre = next(item for item in child["finite_fibres"] if item["kodaira"] == kodaira)
    factor = ring(fibre["factor"])
    assert factor.degree() == 1
    point = -factor[0] / factor[1]
    u_ring = PolynomialRing(QQ, "u")
    u = u_ring.gen()
    u_field = u_ring.fraction_field()

    def translate(value):
        return u_field(u_ring(value.numerator()(point + u))) / u_field(
            u_ring(value.denominator()(point + u))
        )

    return factor, point, u_ring, u, translate


ii_factor, ii_point, ii_ring, ii_u, ii_translate = local_data("II*")
iv_factor, iv_point, iv_ring, iv_u, iv_translate = local_data("IV*")
modulus = ii_factor**2 * iv_factor**(2 if args.nef_ivstar else 3)
assert modulus.degree() == (4 if args.nef_ivstar else 5) and ii_factor.gcd(iv_factor) in QQ

# At II*, q has zero image in QQ[u]/(u^2,X,Y): m and p have the same
# cusp restriction.  Thus the two II* rows are exactly C=0 mod u^2.
ii_modulus = ii_u**2

if args.nef_ivstar:
    # The nef IV* target has the same colength-two cusp quotient as II*.
    # Since p is the cusp residue of m and h is a unit, q=(m-p)/h vanishes
    # in QQ[u]/(u^2,X,Y) at both additive fibres.
    iv_modulus = iv_u**2
else:
    # At the dominant D13 IV* target, Y=-c*u^2, uX=X^2=u^3=0.  Hence
    # q = c/(x(S)(0)h)*u^2 + (-y(S)(0)/(x(S)(0)^2 h))*X.
    B_curve = polynomial(ring, child["minimal_short_weierstrass"]["B_coefficients_low_to_high"])
    iv_b = iv_ring(B_curve(iv_point + iv_u))
    unit_b = iv_b // iv_u**4
    c = QQ(unit_b(0)).sqrt()
    assert c and c**2 == unit_b(0)
    iv_x, iv_y = iv_translate(sx), iv_translate(sy)
    iv_h = iv_ring(h(iv_point + iv_u))
    assert iv_h(0) and iv_x(0) and iv_y(0)
    iv_modulus = iv_u**3
    q_u = reduce_rational(c / (iv_x(0) * iv_h(0)) * iv_u**2, iv_modulus)
    q_x = -QQ(iv_y(0)) / (QQ(iv_x(0))**2 * QQ(iv_h(0)))
    assert q_x


def jets(value, modulus, count):
    return tuple(reduce_rational(value, modulus)[index] for index in range(count))


pairs = [(T**degree, ring.zero()) for degree in range(modulus.degree())]
pairs += [(ring.zero(), T**degree) for degree in range(modulus.degree())]
rows = []
for jet in range(2):
    rows.append([
        jets(ii_translate(field(C)), ii_modulus, 2)[jet]
        for B, C in pairs
    ])
if args.nef_ivstar:
    for jet in range(2):
        rows.append([jets(iv_translate(field(C)), iv_modulus, 2)[jet] for B, C in pairs])
else:
    for jet in range(3):
        rows.append([
            (jets(iv_translate(field(C)), iv_modulus, 3)[jet]
             + jets(iv_translate(field(B)), iv_modulus, 3)[jet] * q_u[0]
             + (jets(iv_translate(field(B)), iv_modulus, 3)[0] * q_u[jet]
                if jet else 0))
            for B, C in pairs
        ])
    rows[-1] = [
        jets(iv_translate(field(C)), iv_modulus, 3)[2]
        + jets(iv_translate(field(B)), iv_modulus, 3)[0] * q_u[2]
        for B, C in pairs
    ]
    rows.append([jets(iv_translate(field(B)), iv_u, 1)[0] * q_x for B, C in pairs])

condition_matrix = matrix(QQ, rows)
expected_rows, expected_columns = (4, 8) if args.nef_ivstar else (6, 10)
assert condition_matrix.dimensions() == (expected_rows, expected_columns)
# Columns B_0,C_0,...,C_4 exhibit an exact nonsingular minor: the five C
# columns are Hermite interpolation at the II*/IV* points, and B_0 supplies
# the independent IV* X quotient direction.  This avoids an impractically
# large rational nullspace computation while proving the full row rank.
pivot_columns = tuple(range(modulus.degree(), 2*modulus.degree())) if args.nef_ivstar else (0, 5, 6, 7, 8, 9)
pivot_minor = condition_matrix.matrix_from_columns(pivot_columns)
assert pivot_minor.det() != 0

# The IV* X row forces B(iv)=0.  For B divisible by iv, B*q has zero IV*
# u-jet because q_u is a multiple of u^2.  The remaining II*/IV* rows then
# force C to vanish modulo ii^2*iv^3.  The displayed diagonal submodule lies
# in the finite kernel and has the same degree-six index, hence is the entire
# module (not merely a Smith-profile calculation).
if args.nef_ivstar:
    b_vector = [QQ.one()] + [QQ.zero()] * (2*modulus.degree()-1)
    assert condition_matrix * matrix(QQ, 2*modulus.degree(), 1, b_vector) == matrix(QQ, expected_rows, 1, 0)
else:
    iv_vector = [iv_factor[index] if index <= iv_factor.degree() else QQ.zero()
                 for index in range(5)] + [QQ.zero()] * 5
    assert condition_matrix * matrix(QQ, 10, 1, iv_vector) == matrix(QQ, 6, 1, 0)
assert ii_factor.degree() == iv_factor.degree() == 1

payload = {
    "schema": "elkies-k3.h92-q6-child-q8-finite-q-module-qq.v2",
    "status": "PASS_EXACT_FINITE_Q_CONDITIONS",
    "inputs": {"child": digest(CHILD), "marking": digest(MARKING)},
    "smooth_frame": "q=(m-p)/h; finite coefficient pair C+B*q",
    "finite_conditions": {
        "CRT_modulus_degree": int(modulus.degree()),
        "rows": int(condition_matrix.nrows()),
        "rank": int(expected_rows),
        "residue_ambient_dimension": int(condition_matrix.ncols()),
        "residue_kernel_dimension": int(expected_columns-expected_rows),
        "nonzero_minor_columns": ["C{}".format(index) for index in range(modulus.degree())] if args.nef_ivstar else ["B0", "C0", "C1", "C2", "C3", "C4"],
        "ii_star": "q=0 in QQ[u]/(u^2,X,Y)",
        "iv_star": "q=0 in QQ[u]/(u^2,X,Y)" if args.nef_ivstar else "q=c*u^2/(x(S)(0)h)+(-y(S)(0)/(x(S)(0)^2h))*X",
    },
    "module": {
        "rank": int(2),
        "finite_codimension": int(expected_rows),
        "determinant_degree": int(expected_rows),
        "smith_degrees": [int(0), int(modulus.degree())] if args.nef_ivstar else [int(1), int(5)],
        "exact_basis": [["1", "0"], ["0", "f_II*^2*f_IV*^2"]] if args.nef_ivstar else [["f_IV*", "0"], ["0", "f_II*^2*f_IV*^3"]],
    },
    "boundary": (
        "This is the complete exact finite module in the q frame. It has no "
        "infinity condition and therefore does not construct a q8 pencil, D13 "
        "equation, rootless bisection, collision, or rank statement."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6CHILDQ8FINITEQQ|nef_ivstar={}|rows={}|rank={}|module_rank=2|determinant_degree={}"
    "|status=PASS_EXACT_FINITE_Q_CONDITIONS".format(
        int(args.nef_ivstar), expected_rows, expected_rows, expected_rows
    ),
    flush=True,
)
