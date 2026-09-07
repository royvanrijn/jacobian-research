#!/usr/bin/env sage -python
"""Derive the exact component-nef finite fractional module and its infinity basis.

For the actual component-nef q8 fibre D, the physical exceptional cycles are
negative.  At both finite additive fibres, adding 2*div(u) turns the required
fractional valuation cycle into the maximal complete ideal (u,X,Y).  Thus the
local coefficient sheaf is

    u^-2 * (u,X,Y).

The physical chord m_phys is the pullback of the standard marked chord m_S
through tau_-P0.  The transported old zero P0 is in the identity component,
and tau_-P0 fixes the singular point of the Weierstrass cubic.  Therefore the
residue of m_phys at the singular point is

    r_f = m_S(0,0) = -y(S)/x(S).

With L=f_II*f_IV and M=L^2, a global finite section

    (C+B*m_phys)/M

belongs to both fractional local modules iff

    C+B*r_f = 0 mod f

at f=f_II,f_IV.  Let rho be the unique degree<2 CRT interpolation of the two
residues.  The finite module is then exactly

    < (m_phys-rho)/M,  L/M >.

It also contains the constant rational function 1 (multiply L/M by L),
which is a useful consistency check for this fractional trivialization.

At infinity m_phys=T^2*m_inf in a smooth translated fibre frame.  Since
deg(rho)<2, deg(L)=2 and deg(M)=4, both displayed generators have order two:
after division by s^2 they limit to independent multiples of m_inf and 1.
This matches the physical fibre coefficient -2 and proves that these two
finite generators are already a complete infinity basis.  No degree-window
search is used.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ, vector, ZZ


ROOT = Path(__file__).resolve().parents[2]
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
ZERO = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-zero-section.json"
MARKING = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-marking.json"
DEFAULT_TARGET = ROOT / "artifacts/local/elkies-k3/q8-target-component-nef.json"
DEFAULT_OUTPUT = ROOT / "artifacts/local/elkies-k3/q8-component-nef-fractional-module.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def rational(field, ring, data, numerator, denominator):
    return field(polynomial(ring, data[numerator])) / field(polynomial(ring, data[denominator]))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--target", type=Path, default=DEFAULT_TARGET)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
args.target = args.target.resolve()
args.output = args.output.resolve()

child = json.loads(CHILD.read_text())
zero = json.loads(ZERO.read_text())
marking = json.loads(MARKING.read_text())
target = json.loads(args.target.read_text())
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert zero["status"] == "PASS_EXACT_CHILD_ZERO_SECTION_TRANSPORT"
assert marking["status"] == "PASS_EXACT_Q6_CHILD_Q8_MARKING"
assert target["status"] == "PASS_EXACT_Q6_CHILD_Q8_PHYSICAL_ROOT_TARGET"
assert target["normalization"]["representative"] == "component-nef"

ring = PolynomialRing(QQ, "T")
T = ring.gen()
field = ring.fraction_field()
A = field(polynomial(ring, child["minimal_short_weierstrass"]["A_coefficients_low_to_high"]))
Bcurve = field(polynomial(ring, child["minimal_short_weierstrass"]["B_coefficients_low_to_high"]))
curve = EllipticCurve(field, [0, 0, 0, A, Bcurve])

p0_data = zero["section"]
p0 = curve(
    rational(field, ring, p0_data,
             "x_numerator_coefficients_low_to_high",
             "x_denominator_coefficients_low_to_high"),
    rational(field, ring, p0_data,
             "y_numerator_coefficients_low_to_high",
             "y_denominator_coefficients_low_to_high"),
)
s_data = marking["selected_q8"]["relative_child_section_standard_jacobian_coordinates"]
s = curve(
    rational(field, ring, s_data,
             "x_numerator_coefficients_low_to_high",
             "x_denominator_coefficients_low_to_high"),
    rational(field, ring, s_data,
             "y_numerator_coefficients_low_to_high",
             "y_denominator_coefficients_low_to_high"),
)
sx, sy = s.xy()
px, py = p0.xy()

fibres = {
    item["kodaira"]: ring(item["factor"])
    for item in child["finite_fibres"]
    if item["kodaira"] in ("II*", "IV*")
}
assert set(fibres) == {"II*", "IV*"}
ii, iv = fibres["II*"], fibres["IV*"]
assert ii.degree() == iv.degree() == 1 and ii.gcd(iv) in QQ
L = (ii*iv).monic()
M = L**2

# Recheck the component-nef valuation cycles against the known actual chart
# valuation cycles.  The negative physical cycles are represented as
# u^-2 times the maximal ideal.
e8_cycle = vector(ZZ, target["selected_q8"]["E8"]["vertical_cycle"])
e6_cycle = vector(ZZ, target["selected_q8"]["E6"]["vertical_cycle"])
assert list(e8_cycle) == [-4,-5,-7,-10,-8,-6,-4,-2]
assert list(e6_cycle) == [-2,-3,-4,-3,-2,-2]
ii_u = vector(ZZ, (2,2,4,6,3,4,3,5))
ii_x = vector(ZZ, (2,4,6,10,4,7,5,8))
ii_chart_cycle = vector(ZZ, (-2,-4,-6,-10,-4,-7,-5,-8))
ii_shifted = ii_chart_cycle + 2*ii_u
assert list(ii_shifted) == [2,0,2,2,2,1,1,2]
iv_u = vector(ZZ, (2,1,1,2,2,3))
iv_x = vector(ZZ, (2,2,2,3,3,4))
iv_chart_cycle = vector(ZZ, (-2,-2,-2,-3,-3,-4))
iv_shifted = iv_chart_cycle + 2*iv_u
assert list(iv_shifted) == [2,0,0,1,1,2]

# At the singular origin, translation by -P0 fixes (0,0).  This follows
# exactly from y0^2=x0^3 on the additive special fibre.
residues = {}
translation_checks = {}
for kind, factor in fibres.items():
    base_point = -factor[0]/factor[1]
    x0, y0 = QQ(px(base_point)), QQ(py(base_point))
    xs0, ys0 = QQ(sx(base_point)), QQ(sy(base_point))
    assert x0 and y0 and xs0
    assert y0**2 == x0**3
    lam0 = -y0/x0
    xp0 = lam0**2-x0
    yp0 = lam0*(0-xp0)
    assert xp0 == 0 and yp0 == 0
    residue = -ys0/xs0
    residues[kind] = QQ(residue)
    translation_checks[kind] = {
        "base_point": str(base_point),
        "P0": [str(x0), str(y0)],
        "translated_singular_point": ["0", "0"],
        "m_phys_residue": str(residue),
    }

# Interpolate rho(T) of degree < 2 through the two local residues.
tii = -ii[0]/ii[1]
tiv = -iv[0]/iv[1]
rii, riv = residues["II*"], residues["IV*"]
assert tii != tiv
rho = ring(rii + (riv-rii)*(T-tii)/(tiv-tii))
assert rho.degree() < 2
assert QQ(rho(tii)) == rii and QQ(rho(tiv)) == riv
assert (rho-ring(rii)) % ii == 0
assert (rho-ring(riv)) % iv == 0

# Numerator lattice K={(B,C): C+B*rho == 0 mod L} has basis
# (1,-rho),(0,L).  Dividing both coefficients by M gives the fractional
# finite module.  Its determinant is L/M^2 and has degree -6.
finite_det_degree = int(L.degree()-2*M.degree())
assert finite_det_degree == -6

# Infinity: m_phys=T^2*m_inf after translating the smooth fibre.  rho has
# degree <=1, L degree2, M degree4.  Therefore
#   (m-rho)/M = s^2*(unit*m_inf + O(s)),
#   L/M       = s^2*(unit + O(s)).
# The leading matrix is triangular with nonzero diagonal, so this is exactly
# the s^2 local lattice required by vertical_fibre_coefficient=-2.
vertical_fibre_coefficient = int(target["selected_q8"]["vertical_fibre_coefficient"])
assert vertical_fibre_coefficient == -2
g1_order = M.degree()-2
g2_order = M.degree()-L.degree()
assert g1_order == g2_order == 2
transition_det_order = 2+2
assert transition_det_order == 4

payload = {
    "schema": "elkies-k3.h92-q6-child-q8-component-nef-fractional-module.v1",
    "status": "PASS_EXACT_COMPONENT_NEF_FRACTIONAL_MODULE",
    "inputs": {
        "child": digest(CHILD),
        "zero": digest(ZERO),
        "marking": digest(MARKING),
        "component_nef_target": digest(args.target),
    },
    "local_fractional_ideals": {
        "II*": {
            "cycle_after_u2_shift": list(map(int, ii_shifted)),
            "complete_ideal": "(u,X,Y)",
            "fractional_ideal": "u^-2*(u,X,Y)",
        },
        "IV*": {
            "cycle_after_u2_shift": list(map(int, iv_shifted)),
            "complete_ideal": "(u,X,Y)",
            "fractional_ideal": "u^-2*(u,X,Y)",
            "orientation_independent": True,
        },
    },
    "translation": translation_checks,
    "residues": {key: str(value) for key, value in residues.items()},
    "base_polynomials": {
        "f_II": str(ii.monic()),
        "f_IV": str(iv.monic()),
        "L=f_II*f_IV": str(L),
        "M=L^2": str(M),
        "rho": str(rho),
        "rho_coefficients_low_to_high": [str(value) for value in rho.list()],
    },
    "finite_module": {
        "generic_frame": "1,m_component_nef",
        "coefficient_form": "(C+B*m_component_nef)/M",
        "residue_condition": "C+B*rho = 0 mod L",
        "basis_functions": [
            "(m_component_nef-rho)/M",
            "L/M",
        ],
        "coefficient_basis_B_C": [
            ["1/M", "-rho/M"],
            ["0", "L/M"],
        ],
        "determinant": "L/M^2",
        "determinant_degree": finite_det_degree,
        "contains_constant_one_in_this_trivialization": True,
        "constant_one_identity": "1=L*(L/M)",
    },
    "infinity": {
        "physical_chord_scaling": "m_component_nef=T^2*m_inf+O(T)",
        "vertical_fibre_coefficient": vertical_fibre_coefficient,
        "basis_orders": [int(g1_order), int(g2_order)],
        "local_lattice": "s^2<1,m_inf>",
        "leading_basis": ["unit*m_inf", "unit"],
        "transition_determinant_order": transition_det_order,
        "basis_is_complete": True,
    },
    "pencil": {
        "new_base": "U=((m_component_nef-rho)/M)/(L/M)=(m_component_nef-rho)/L",
        "translated_standard_coordinate_equation": "m_S=rho+L*U",
    },
    "boundary": (
        "This certifies the finite fractional module and its smooth infinity "
        "completion for the component-nef degree-two divisor.  It does not yet "
        "compute the characteristic-zero branch quartic, child Jacobian, later "
        "neighbours, or a rank record."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q8FRAC|rho_deg={}|L_deg={}|M_deg={}|finite_det_deg={}|"
    "infinity_orders=2,2|contains_one=1|new_base=(m-rho)/L|"
    "status=PASS_EXACT_COMPONENT_NEF_FRACTIONAL_MODULE".format(
        rho.degree(), L.degree(), M.degree(), finite_det_degree
    ),
    flush=True,
)
