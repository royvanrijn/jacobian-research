#!/usr/bin/env sage -python
"""Evaluate the corrected child q=8 chord in the II*/IV* vertical quotients.

The q=8 marking on the explicit q6 child has generic basis ``<1,m>``, with

    m=(y+y(S))/(x-x(S)).

The vertical II* and IV* ideals have now been determined in the actual
Weierstrass charts.  This checker reduces ``m`` in their finite quotient
rings and exports the resulting linear conditions on a local coefficient
pair ``a(u)+b(u)m``.  It does not choose a global base-function ambient or
assemble the two fibres with the degree-46 smooth collision condition.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
MARKING = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-marking.json"
IISTAR = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-iistar-vertical-ideal.json"
IVSTAR = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-ivstar-orientation.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-additive-chord-blocks.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def rational_function(field, ring, data, numerator_key, denominator_key):
    return field(polynomial(ring, data[numerator_key])) / field(
        polynomial(ring, data[denominator_key])
    )


def reduce_rational(value, modulus):
    """Reduce a QQ(u) element modulo a polynomial at which its denominator is a unit."""

    ring = modulus.parent()
    numerator = ring(value.numerator())
    denominator = ring(value.denominator())
    assert denominator.gcd(modulus) in QQ
    return (numerator * denominator.inverse_mod(modulus)).mod(modulus)


def coefficient_vector(value, degree):
    return vector(QQ, [value[index] for index in range(degree)])


def local_block(residue, dimension, name):
    """Matrix for a_0+...+a_{d-1}u^(d-1)+b(u)*residue in QQ[u]/u^d."""

    ring = residue.parent()
    u = ring.gen()
    modulus = u**dimension
    labels = tuple(
        [("a", index) for index in range(dimension)]
        + [("b", index) for index in range(dimension)]
    )
    columns = []
    for kind, exponent in labels:
        value = u**exponent if kind == "a" else u**exponent * residue
        columns.append(coefficient_vector(value.mod(modulus), dimension))
    return {
        "name": name,
        "basis": ["u^{}".format(index) for index in range(dimension)],
        "ambient": [[kind, exponent] for kind, exponent in labels],
        "matrix": matrix(QQ, dimension, len(columns), lambda row, column: columns[column][row]),
    }


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--child", type=Path, default=CHILD)
parser.add_argument("--marking", type=Path, default=MARKING)
parser.add_argument("--iistar", type=Path, default=IISTAR)
parser.add_argument("--ivstar", type=Path, default=IVSTAR)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
for name in ("child", "marking", "iistar", "ivstar", "output"):
    setattr(args, name, getattr(args, name).resolve())

child = json.loads(args.child.read_text())
marking = json.loads(args.marking.read_text())
iistar = json.loads(args.iistar.read_text())
ivstar = json.loads(args.ivstar.read_text())
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert marking["status"] == "PASS_EXACT_Q6_CHILD_Q8_MARKING"
assert iistar["status"] == "PASS_EXACT_Q6_CHILD_Q8_IISTAR_VERTICAL_IDEAL"
assert ivstar["status"] in {
    "PASS_EXACT_Q6_CHILD_Q8_IVSTAR_ORIENTATION",
    "PASS_EXACT_Q6_CHILD_Q8_IVSTAR_VERTICAL_IDEAL_PAIR",
}
assert iistar["complete_ideal"]["generators"] == ["Y", "X", "u^2"]
if "selected_q8_ivstar_module" in ivstar:
    iv_generators = ivstar["selected_q8_ivstar_module"]["complete_ideal_generators"]
else:
    candidates = ivstar["orientation_candidates"]
    assert candidates and all(item["generators"] == candidates[0]["generators"] for item in candidates)
    iv_generators = candidates[0]["generators"]
nef_ivstar = iv_generators == ["Y", "X", "u^2"]
if not nef_ivstar:
    assert iv_generators == ["Y+c*u^2", "u*X", "X^2", "u^3"]

T_ring = PolynomialRing(QQ, "T")
T = T_ring.gen()
T_field = T_ring.fraction_field()
A = polynomial(T_ring, child["minimal_short_weierstrass"]["A_coefficients_low_to_high"])
B = polynomial(T_ring, child["minimal_short_weierstrass"]["B_coefficients_low_to_high"])
section_data = marking["selected_q8"]["relative_child_section_standard_jacobian_coordinates"]
sx = rational_function(
    T_field, T_ring, section_data,
    "x_numerator_coefficients_low_to_high", "x_denominator_coefficients_low_to_high",
)
sy = rational_function(
    T_field, T_ring, section_data,
    "y_numerator_coefficients_low_to_high", "y_denominator_coefficients_low_to_high",
)
assert sy**2 == sx**3 + A*sx + B


def local_section(fibre_name, expected_orders):
    fibre = next(item for item in child["finite_fibres"] if item["kodaira"] == fibre_name)
    factor = T_ring(fibre["factor"])
    assert factor.degree() == 1
    base_point = -factor[0] / factor[1]
    u_ring = PolynomialRing(QQ, "u")
    u = u_ring.gen()
    u_field = u_ring.fraction_field()

    def translate(value):
        return u_field(u_ring(value.numerator()(base_point + u))) / u_field(
            u_ring(value.denominator()(base_point + u))
        )

    x_u, y_u = translate(sx), translate(sy)
    assert x_u.valuation() == y_u.valuation() == 0
    assert (x_u(0), y_u(0)) != (0, 0)
    A_u = u_ring(A(base_point + u))
    B_u = u_ring(B(base_point + u))
    a, remainder = A_u.quo_rem(u**expected_orders[0])
    assert not remainder and a(0)
    b, remainder = B_u.quo_rem(u**expected_orders[1])
    assert not remainder and b(0)
    return u_ring, u, x_u, y_u, a, b


# II*: I=(u^2,X,Y), so only the two-jet at the cusp X=Y=0 remains.
ii_ring, ii_u, ii_x, ii_y, ii_a, ii_b = local_section("II*", (4, 5))
ii_modulus = ii_u**2
ii_residue = reduce_rational(-ii_y/ii_x, ii_modulus)
ii_block = local_block(ii_residue, 2, "II* chord quotient R/(u^2,X,Y)")
assert ii_block["matrix"].rank() == 2
assert ii_residue[0]

# IV*: the dominant D13 target has the arm ideal
# (Y+c*u^2,uX,X^2,u^3); the nef target has the invariant ideal (u^2,X,Y).
iv_ring, iv_u, iv_x, iv_y, iv_a, iv_b = local_section("IV*", (3, 4))
if nef_ivstar:
    iv_modulus = iv_u**2
    iv_u_residue = reduce_rational(-iv_y/iv_x, iv_modulus)
    iv_u_block = local_block(iv_u_residue, 2, "IV* quotient R/(u^2,X,Y)")
    iv_labels = iv_u_block["ambient"]
    iv_matrix = iv_u_block["matrix"]
    x_coefficient = None
    assert iv_matrix.rank() == 2
else:
    c = QQ(iv_b(0)).sqrt()
    assert c and c**2 == iv_b(0)
    iv_modulus = iv_u**3
    ratio = reduce_rational(-iv_y/iv_x, iv_modulus)
    inverse_x_constant = QQ(iv_x(0))**-1
    x_coefficient = -QQ(iv_y(0)) * inverse_x_constant**2
    u2_correction = c * inverse_x_constant
    iv_u_residue = iv_ring(ratio + u2_correction*iv_u**2).mod(iv_modulus)
    assert iv_x(0) and iv_y(0)
    assert x_coefficient
    iv_u_block = local_block(iv_u_residue, 3, "IV* u-part quotient")
    iv_labels = iv_u_block["ambient"]
    iv_matrix = matrix(QQ, 4, len(iv_labels), 0)
    iv_matrix[:3, :] = iv_u_block["matrix"]
    for column, (kind, exponent) in enumerate(iv_labels):
        if kind == "b" and exponent == 0:
            iv_matrix[3, column] = x_coefficient
    assert iv_matrix.rank() == 4

payload = {
    "schema": "elkies-k3.h92-q6-child-q8-additive-chord-blocks.v1",
    "status": "PASS_EXACT_Q6_CHILD_Q8_ADDITIVE_CHORD_BLOCKS",
    "inputs": {
        "child_jacobian": {"path": str(args.child.relative_to(ROOT)), "sha256": digest(args.child)},
        "q8_marking": {"path": str(args.marking.relative_to(ROOT)), "sha256": digest(args.marking)},
        "iistar_ideal": {"path": str(args.iistar.relative_to(ROOT)), "sha256": digest(args.iistar)},
        "ivstar_orientation": {"path": str(args.ivstar.relative_to(ROOT)), "sha256": digest(args.ivstar)},
    },
    "chord": "m=(y+y(S))/(x-x(S)), S=corrected marked q8 section",
    "ii_star": {
        "ideal": "(u^2,X,Y)",
        "quotient_basis": ["1", "u"],
        "m_residue_at_X=Y=0_mod_u2": str(ii_residue),
        "coefficient_block": {
            **{key: value for key, value in ii_block.items() if key != "matrix"},
            "matrix": [[str(value) for value in row] for row in ii_block["matrix"].rows()],
            "rank": int(ii_block["matrix"].rank()),
        },
    },
    "iv_star": {
        "ideal": "(u^2,X,Y)" if nef_ivstar else "(Y+c*u^2,u*X,X^2,u^3)",
        "quotient_basis": ["1", "u"] if nef_ivstar else ["1", "u", "u^2", "X"],
        "m_residue": str(iv_u_residue) if nef_ivstar else "{} + ({})*X modulo (u^3,uX,X^2,Y+c*u^2)".format(iv_u_residue, x_coefficient),
        "u_residue": str(iv_u_residue),
        "X_coefficient": None if x_coefficient is None else str(x_coefficient),
        "coefficient_block": {
            "basis": ["1", "u", "u^2", "X"],
            "ambient": iv_labels,
            "matrix": [[str(value) for value in row] for row in iv_matrix.rows()],
            "rank": int(iv_matrix.rank()),
        },
    },
    "boundary": (
        "These are the exact finite additive chord blocks for local coefficient "
        "jets. A global base-function envelope and its compatibility with the "
        "degree-46 smooth collision module have not been derived, so no q8 pencil, "
        "rootless equation, bisection cover, collision, or rank statement follows."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    ("H92Q6CHILDQ8ADDITIVE|II_rows=2|IV_rows={}|II_rank=2|IV_rank={}|"
     "status=PASS_EXACT_Q6_CHILD_Q8_ADDITIVE_CHORD_BLOCKS").format(
         iv_matrix.nrows(), iv_matrix.rank()
     ),
    flush=True,
)
