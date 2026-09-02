#!/usr/bin/env sage-python
"""Certify a bounded rational-surface quadratic-base-change search.

status: ACTIVE_PROOF
claim: exact E6+A1 quadratic-base-change family and bounded degree-(2,2) search
inputs: repeated-fibre rational-base-change audit certificate
outputs: elkies-k3-rational-surface-quadratic-rank-search-v1.json

The search starts with the one-modulus E6+A1 rational elliptic surfaces

    y^2 = x^3 + (u-3)x + c^2*u^2 + u - 2,

whose rational section (-1,c*u) generates the geometric rank-one
Mordell--Weil lattice [1/6].  It branches a quadratic cover at the I2 fibre
u=0 and at u=lambda, and solves a complete degree-(2,2) twist-section ansatz.

Two rational solution curves result.  The preferred curve has especially
small formulas.  Its quadratic pullback is a non-isotrivial one-modulus K3
family with fibres 2IV*+I4+4I1, two independent rational sections, height
matrix diag(1/3,3), and generic Picard rank 19.  The script also imports the
already-certified Golay/NS0031 rational-surface descents as controls.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import (
    CartanMatrix,
    EllipticCurve,
    Infinity,
    PolynomialRing,
    QQ,
    ZZ,
    ceil,
    matrix,
    vector,
    zero_matrix,
)


ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "artifacts/generated-results"
DEFAULT_OUTPUT = GEN / "elkies-k3-rational-surface-quadratic-rank-search-v1.json"
CONTROL = GEN / "elkies-k3-repeated-fibre-rational-base-change-audit-v1.json"


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compact_height(curve, point, multiplier, chi):
    """Recover a Shioda height after clearing all component groups."""

    coordinate = (multiplier * point)[0]
    numerator = coordinate.numerator()
    denominator = coordinate.denominator()
    if not denominator.is_square():
        raise ArithmeticError("x-coordinate denominator is not a square")
    finite_poles = denominator.degree() // 2
    infinity_poles = ceil((numerator.degree() - 2 * chi) / 2)
    intersection = max(finite_poles, infinity_poles)
    height = (2 * chi + 2 * intersection) / QQ(multiplier**2)
    return {
        "multiplier": multiplier,
        "x_numerator_degree": int(numerator.degree()),
        "x_denominator_degree": int(denominator.degree()),
        "intersection_with_zero": int(intersection),
        "height": str(height),
    }


# Work over QQ(k)(u).
K_RING = PolynomialRing(QQ, "k")
k = K_RING.gen()
K = K_RING.fraction_field()
U_RING = PolynomialRing(K, "u")
u = U_RING.gen()


def preferred_family():
    c = K(2 * k / (3 * k**2 - 4))
    lam = K(-(k - 2) * (k + 2) * (3 * k**2 - 4) / 4)
    z = K(2 / (3 * k**2 - 4))
    p = z**2
    q = K(k**2 / (3 * k**2 - 4))
    r = K(-1)
    s = z**3
    v = K(2 * (k**2 + 2) / (3 * k**2 - 4) ** 2)
    w = K(0)
    return c, lam, (p, q, r), (s, v, w)


def secondary_family():
    h = K((3 * k**2 + 16) / (6 * k))
    m = K((16 - 3 * k**2) / (6 * k))
    if 3 * (h**2 - m**2) != 16:
        raise ArithmeticError("secondary conic parameterization failed")
    c = K(m / (h**2 - 8))
    lam = K((8 - h**2) / 3)
    p = K(64 / (h**2 * (h**2 - 8) ** 2))
    q = K(8 * (3 * h**2 + 8) / (3 * h**2 * (h**2 - 8)))
    r = K(2)
    s = K(-512 / (h**3 * (h**2 - 8) ** 3))
    v = K(-32 * (9 * h**2 + 16) / (3 * h**3 * (h**2 - 8) ** 2))
    w = K(-3 * (h**2 + 8) / (h * (h**2 - 8)))
    return c, lam, (p, q, r), (s, v, w), h, m


def verify_twist_family(label, c, lam, x_coefficients, y_coefficients):
    p, q, r = x_coefficients
    s, v, w = y_coefficients
    d = u * (u - lam)
    x_twist = p * u**2 + q * u + r
    y_twist = s * u**2 + v * u + w
    cubic = x_twist**3 + (u - 3) * x_twist + c**2 * u**2 + u - 2
    if d * y_twist**2 != cubic:
        raise ArithmeticError(f"{label}: twist-section identity failed")
    if (c * u) ** 2 != (-1) ** 3 + (u - 3) * (-1) + c**2 * u**2 + u - 2:
        raise ArithmeticError(f"{label}: invariant rational-surface section failed")

    base_discriminant = 4 * (u - 3) ** 3 + 27 * (c**2 * u**2 + u - 2) ** 2
    expected_quadratic = (
        27 * c**4 * u**2 + (54 * c**2 + 4) * u - (108 * c**2 + 9)
    )
    if base_discriminant != u**2 * expected_quadratic:
        raise ArithmeticError(f"{label}: rational-surface discriminant factorization failed")

    return {
        "label": label,
        "parameters": {"c": str(c), "lambda": str(lam)},
        "rational_surface": {
            "equation": "y^2=x^3+(u-3)*x+c^2*u^2+u-2",
            "fibre_profile": "IV*+I2+2I1",
            "root_lattice": "E6+A1",
            "root_rank": 7,
            "mordell_weil_group": "Z",
            "mordell_weil_lattice_gram": [["1/6"]],
            "generator": {"x": "-1", "y": "c*u"},
            "torsion": "trivial",
        },
        "quadratic_cover": {
            "squareclass": "d(u)=u*(u-lambda)",
            "map": "u=lambda/(1-t^2)",
            "branch_values": ["u=0 (I2)", "u=lambda (smooth generically)"],
        },
        "twist_section": {
            "model": "d(u)*y^2=x^3+(u-3)*x+c^2*u^2+u-2",
            "x_coefficients_low_to_high": [str(r), str(q), str(p)],
            "y_coefficients_low_to_high": [str(w), str(v), str(s)],
        },
        "base_changed_k3": {
            "equation": (
                "Y^2=X^3+H^3*(lambda-3H)*X+"
                "H^4*(c^2*lambda^2+lambda*H-2H^2), H=1-t^2"
            ),
            "fibre_profile": "2IV*+I4+4I1",
            "root_lattice": "2E6+A3",
            "root_rank": 15,
            "invariant_section": {"X": "-H^2", "Y": "c*lambda*H^2"},
            "anti_invariant_section": {
                "X": "p*lambda^2+q*lambda*H+r*H^2",
                "Y": "lambda*t*(s*lambda^2+v*lambda*H+w*H^2)",
            },
            "rank_sum_lower_bound": 2,
            "picard_rank_lower_bound": 19,
        },
    }


preferred = preferred_family()
preferred_record = verify_twist_family("E6A1-clean-rminus1", *preferred)
secondary = secondary_family()
secondary_record = verify_twist_family("E6A1-secondary-r2", *secondary[:4])
secondary_record["parameters"]["auxiliary_h"] = str(secondary[4])
secondary_record["parameters"]["auxiliary_m"] = str(secondary[5])
secondary_record["parameters"]["conic"] = "3*(h^2-m^2)=16"


# Verify the polynomial K3 sections for the preferred family over QQ(k)(t).
T_RING = PolynomialRing(K, "t")
t = T_RING.gen()
H = 1 - t**2
c, lam, (p, q, r), (s, v, w) = preferred
A_k3 = H**3 * (lam - 3 * H)
B_k3 = H**4 * (c**2 * lam**2 + lam * H - 2 * H**2)
P_x = -H**2
P_y = c * lam * H**2
Q_x = p * lam**2 + q * lam * H + r * H**2
Q_y = lam * t * (s * lam**2 + v * lam * H + w * H**2)
for name, x_coordinate, y_coordinate in (
    ("invariant", P_x, P_y),
    ("anti-invariant", Q_x, Q_y),
):
    if y_coordinate**2 != x_coordinate**3 + A_k3 * x_coordinate + B_k3:
        raise ArithmeticError(f"preferred K3 {name} section failed")
if P_x(-t) != P_x or P_y(-t) != P_y or Q_x(-t) != Q_x or Q_y(-t) != -Q_y:
    raise ArithmeticError("deck-character identities failed")


# Complete degree-(2,2) ansatz elimination.  The constant coefficient forces
# x(0) in {-1,2}.  These are the two rational branches retained above.
ELIM = PolynomialRing(QQ, names=("z", "L", "q"), order="lex")
z0, L0, q0 = ELIM.gens()

e3_minus = 5 * L0**3 * z0**6 + 9 * L0**2 * q0 * z0**4 + 3 * L0 * q0**2 * z0**2 - 12 * L0 * z0**4 - q0**3 - 12 * q0 * z0**2 + 8 * z0**2
e1_minus = L0 * (L0 * z0**2 + q0 - 2 * z0) ** 2 * (L0 * z0**2 + q0 + 2 * z0) ** 2
resultant_minus = e3_minus.resultant(e1_minus, q0).factor()
expected_minus = (
    4096
    * L0**3
    * z0**8
    * (3 * L0 * z0**2 - 4 * z0 + 1) ** 2
    * (3 * L0 * z0**2 + 4 * z0 + 1) ** 2
)
if resultant_minus.prod() != expected_minus:
    raise ArithmeticError("r=-1 ansatz resultant changed")

e3_plus = 5 * L0**3 * z0**6 + 9 * L0**2 * q0 * z0**4 + 3 * L0 * q0**2 * z0**2 + 24 * L0 * z0**4 - q0**3 + 24 * q0 * z0**2 + 8 * z0**2
e1_plus = (
    3 * L0**5 * z0**8
    + 12 * L0**4 * q0 * z0**6
    + 18 * L0**3 * q0**2 * z0**4
    + 48 * L0**3 * z0**6
    + 12 * L0**2 * q0**3 * z0**2
    + 96 * L0**2 * q0 * z0**4
    + 3 * L0 * q0**4
    + 48 * L0 * q0**2 * z0**2
    + 192 * L0 * z0**4
    + 192 * q0 * z0**2
    + 64 * z0**2
)
resultant_plus = e3_plus.resultant(e1_plus, q0).factor()
expected_plus = 4096 * z0**6 * (27 * L0**3 * z0**2 - 72 * L0**2 * z0**2 + 64)
if resultant_plus.prod() != expected_plus:
    raise ArithmeticError("r=2 ansatz resultant changed")


# Exact height and non-torsion calibration at k=1.  The component-group
# exponent is lcm(3,4)=12 for 2IV*+I4.
RT = PolynomialRing(QQ, "tt")
tt = RT.gen()
FT = RT.fraction_field()
k1 = QQ(1)
c1 = QQ(c(k1))
lam1 = QQ(lam(k1))
p1, q1, r1 = (QQ(value(k1)) for value in (p, q, r))
s1, v1, w1 = (QQ(value(k1)) for value in (s, v, w))
H1 = 1 - tt**2
curve_k3 = EllipticCurve(
    FT,
    [
        H1**3 * (lam1 - 3 * H1),
        H1**4 * (c1**2 * lam1**2 + lam1 * H1 - 2 * H1**2),
    ],
)
point_invariant = curve_k3(-H1**2, c1 * lam1 * H1**2)
point_anti = curve_k3(
    p1 * lam1**2 + q1 * lam1 * H1 + r1 * H1**2,
    lam1 * tt * (s1 * lam1**2 + v1 * lam1 * H1 + w1 * H1**2),
)
height_invariant = compact_height(curve_k3, point_invariant, 12, 2)
height_anti = compact_height(curve_k3, point_anti, 12, 2)
if height_invariant["height"] != "1/3" or height_anti["height"] != "3":
    raise ArithmeticError("preferred-family height calibration changed")

# Pin the integral Neron--Severi lattice, not merely its rational height
# decomposition.  The basis is
#
#   O,F,E6_1[1..6],E6_2[1..6],A3[1..3],P0,P1.
#
# P0 meets the 4/3 component in each IV* and the middle I4 component;
# P1 meets only the middle I4 component.  Character orthogonality and the
# I4 local cross correction give P0.P1=1.
ns = zero_matrix(ZZ, 19)
ns[0, 0] = -2
ns[0, 1] = ns[1, 0] = 1
e6 = CartanMatrix(["E", 6])
a3 = CartanMatrix(["A", 3])
ns[2:8, 2:8] = -e6
ns[8:14, 8:14] = -e6
ns[14:17, 14:17] = -a3
P0_index, P1_index = 17, 18
for section_index in (P0_index, P1_index):
    ns[section_index, section_index] = -2
    ns[1, section_index] = ns[section_index, 1] = 1
ns[P0_index, P1_index] = ns[P1_index, P0_index] = 1
ns[P0_index, 2] = ns[2, P0_index] = 1
ns[P0_index, 8] = ns[8, P0_index] = 1
ns[P0_index, 15] = ns[15, P0_index] = 1
ns[P1_index, 15] = ns[15, P1_index] = 1
if ns.det() != 36 or ns.elementary_divisors()[-2:] != [3, 12]:
    raise ArithmeticError("generic Neron--Severi Gram matrix changed")

# Even overlattices of this displayed lattice correspond to isotropic
# subgroups of its discriminant form.  In Smith coordinates Z/3 + Z/12 the
# form is represented by [[4/3,1/3],[1/3,-1/4]].  Its only nonzero isotropic
# elements generate two order-three subgroups, so index three is the only
# possible saturation defect.
discriminant_isotropic = []
for first in range(3):
    for second in range(12):
        if first == 0 and second == 0:
            continue
        norm = (
            QQ(4) * first**2 / 3
            + QQ(2) * first * second / 3
            - QQ(second**2) / 4
        )
        if norm in ZZ and ZZ(norm) % 2 == 0:
            discriminant_isotropic.append((first, second))
if discriminant_isotropic != [(0, 4), (0, 8), (1, 4), (2, 8)]:
    raise ArithmeticError("generic Neron--Severi isotropic classes changed")

# Eliminate both possible index-three Mordell--Weil saturations on one good
# specialization.  If a generic class (a*P0+b*P1)/3 existed, it would extend
# over the good (k,t)=(1,3) fibre.  The four projective combinations modulo
# three are all non-divisible in this exact elliptic curve over QQ.
t3 = QQ(3)
H3 = 1 - t3**2
curve_t3 = EllipticCurve(
    QQ,
    [
        H3**3 * (lam1 - 3 * H3),
        H3**4 * (c1**2 * lam1**2 + lam1 * H3 - 2 * H3**2),
    ],
)
point_invariant_t3 = curve_t3(-H3**2, c1 * lam1 * H3**2)
point_anti_t3 = curve_t3(
    p1 * lam1**2 + q1 * lam1 * H3 + r1 * H3**2,
    lam1 * t3 * (s1 * lam1**2 + v1 * lam1 * H3 + w1 * H3**2),
)
saturation_combinations = {
    "P0": point_invariant_t3,
    "P1": point_anti_t3,
    "P0+P1": point_invariant_t3 + point_anti_t3,
    "P0-P1": point_invariant_t3 - point_anti_t3,
}
if any(point.is_divisible_by(3) for point in saturation_combinations.values()):
    raise ArithmeticError("generic Mordell--Weil lattice is not 3-saturated")

# A good fibre of the twist at (k,u)=(1,1) proves that the anti-invariant
# generic section is non-torsion.
u1 = QQ(1)
d1 = u1 * (u1 - lam1)
x1 = p1 * u1**2 + q1 * u1 + r1
y1 = s1 * u1**2 + v1 * u1 + w1
twist_fibre = EllipticCurve(
    QQ,
    [
        d1**2 * (u1 - 3),
        d1**3 * (c1**2 * u1**2 + u1 - 2),
    ],
)
twist_point = twist_fibre(d1 * x1, d1**2 * y1)
if twist_point.order() != Infinity:
    raise ArithmeticError("specialized anti-invariant point is torsion")


# A marked-fibre cross-ratio surrogate: if alpha,beta are the residual I1
# values then the product of their two squared preimage coordinates is
# Q(lambda)/Q(0).  Its nonconstant formula proves non-isotriviality of the
# marked K3 family.
residual = lambda value: (
    27 * c**4 * value**2 + (54 * c**2 + 4) * value - (108 * c**2 + 9)
)
configuration_invariant = K(residual(lam) / residual(K(0)))
expected_configuration_invariant = K(
    (2 * k**2 + 1) * (9 * k**2 - 20) ** 2 / (9 * (3 * k**2 + 4) ** 2)
)
if configuration_invariant != expected_configuration_invariant:
    raise ArithmeticError("marked fibre-configuration invariant changed")
if configuration_invariant.derivative() == 0:
    raise ArithmeticError("preferred K3 family became isotrivial")

preferred_record["base_changed_k3"].update(
    {
        "mordell_weil_rank_exact_generic": 2,
        "mordell_weil_height_gram": [["1/3", "0"], ["0", "3"]],
        "mordell_weil_determinant": "1",
        "neron_severi_determinant_absolute": "36",
        "picard_rank_exact_generic": 19,
        "nonisotrivial_marked_configuration_invariant": str(configuration_invariant),
        "rank_decomposition": "rank E(QQ(k)(u))=1 plus rank E^(d)(QQ(k)(u))=1",
        "neron_severi_smith_nontrivial": [3, 12],
        "mordell_weil_saturation": "exact; no index-three overlattice",
    }
)
preferred_record["arithmetic_open"] = {
    "rational_parameter": "k in QQ",
    "excluded_rational_values": ["-2", "0", "2"],
    "residual_I1_discriminant": "16*(9*k^4+12*k^2+16)^3/(3*k^2-4)^6",
    "branch_smoothness_value": "-(2*k^2+1)*(9*k^2-20)^2/(3*k^2-4)^2",
    "comment": (
        "For rational k outside {-2,0,2}, the displayed denominators, branch degree, "
        "and generic rational-surface/K3 fibre profiles are valid; further special-fibre "
        "rank jumps are not excluded."
    ),
}
preferred_record["height_calibration_k_equals_1"] = {
    "invariant": height_invariant,
    "anti_invariant": height_anti,
    "twist_good_fibre_u_equals_1": {
        "curve_a_invariants": [str(value) for value in twist_fibre.a_invariants()],
        "point": [str(twist_point[0]), str(twist_point[1])],
        "order": "infinite",
    },
    "saturation_fibre_t_equals_3": {
        "curve_a_invariants": [str(value) for value in curve_t3.a_invariants()],
        "curve_discriminant": str(curve_t3.discriminant()),
        "points": {
            label: [str(point[0]), str(point[1])]
            for label, point in saturation_combinations.items()
        },
        "all_nondivisible_by_3": True,
    },
}


if not CONTROL.exists():
    raise FileNotFoundError(CONTROL)
control_payload = json.loads(CONTROL.read_text())
control_rows = control_payload["exact_promoted_model_certificates"]
golay_control = next(row for row in control_rows if row["base_field"] == "QQ")
ns0031_control = next(row for row in control_rows if row["base_field"] == "GF(7)")


catalogue = [
    {
        "id": "Kimura-E7",
        "moduli_dimension": 1,
        "equation": "y^2=x^3+(u-3)x+(a*u-2)",
        "rational_section_subfamily": "a=-c^2-2, P=(-a,c*(a-1))",
        "fibre_profile": "III*+3I1",
        "root_lattice": "E7",
        "mw_lattice": [["1/2"]],
        "torsion": "trivial",
        "priority": "high roots, but the E6+A1 seed has a simpler rational generator",
    },
    {
        "id": "Kimura-D7",
        "moduli_dimension": 1,
        "equation": "y^2=x^3+u^2*(u^2-3)x+u^3*(a*u^3+u^2-2)",
        "fibre_profile": "I3*+3I1",
        "root_lattice": "D7",
        "mw_lattice": [["1/4"]],
        "torsion": "trivial generically",
        "priority": "rank-budget seed; a compact QQ generator was not compiled here",
    },
    {
        "id": "Kimura-E6A1-rational-generator",
        "moduli_dimension": 1,
        "equation": "y^2=x^3+(u-3)x+(c^2*u^2+u-2)",
        "section": "(-1,c*u)",
        "fibre_profile": "IV*+I2+2I1",
        "root_lattice": "E6+A1",
        "mw_lattice": [["1/6"]],
        "torsion": "trivial",
        "priority": "selected seed",
    },
    {
        "id": "Golay-rational-control",
        "moduli_dimension": 0,
        "equation_reference": relative(CONTROL),
        "fibre_profile": golay_control["quotient"]["fibre_profile"],
        "root_lattice": "A5+A2",
        "mw_lattice": [["1/2"]],
        "torsion": "Z/3",
        "rank_decomposition": "1 invariant + 1 anti-invariant",
        "status": "exact QQ control; source NS determinant is 20, not Golay 720",
    },
    {
        "id": "NS0031-mod7-control",
        "moduli_dimension": 0,
        "equation_reference": relative(CONTROL),
        "base_field": "GF(7)",
        "fibre_profile": ns0031_control["quotient"]["fibre_profile"],
        "root_lattice": "A7",
        "mw_rank_from_shioda_tate": 1,
        "rank_decomposition": "one nonzero invariant trace + one nonzero anti-invariant direction",
        "status": "finite-field precursor, not an E/QQ(u) entry",
    },
]


payload = {
    "schema": "elkies-k3.rational-surface-quadratic-rank-search.v1",
    "status": "PASS_EXACT_NEW_ONE_MODULUS_PICARD19_QUADRATIC_BASE_CHANGE",
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/certify_rational_surface_quadratic_rank_search.sage"
    ),
    "inputs": {relative(CONTROL): digest(CONTROL)},
    "low_complexity_rational_surface_catalogue": catalogue,
    "bounded_degree_2_2_twist_section_elimination": {
        "ansatz": "x=p*u^2+q*u+r, y=s*u^2+v*u+w, p=z^2, s=z^3",
        "constant_equation": "(r-2)*(r+1)^2=0",
        "r_minus_1_resultant": str(resultant_minus),
        "r_2_resultant": str(resultant_plus),
        "scope": (
            "Complete for nonzero leading coefficient and nonzero lambda in this polynomial "
            "degree-(2,2) ansatz; it is not a complete search over rational-function sections."
        ),
    },
    "new_families": [preferred_record, secondary_record],
    "pareto_selection": {
        "winner": "E6A1-clean-rminus1",
        "objectives": {
            "generic_picard_rank": 19,
            "rank_sum": 2,
            "rational_moduli_dimension": 1,
            "rational_sections": 2,
            "quadratic_map_degree": 2,
            "base_equation_monomial_count": 6,
            "twist_section_degree_pair": [2, 2],
            "arithmetic_parameter": "k in QQ outside an explicit finite bad locus",
        },
        "reason": (
            "It has the same exact rank and Picard outcome as the secondary component, "
            "but substantially smaller c, lambda, and twist-section formulas."
        ),
    },
    "proof_boundary": {
        "proved": (
            "The rational-surface and twist section identities, quadratic map, deck characters, "
            "degree-(2,2) ansatz resultants, fibre-root budget, exact specialized heights, "
            "non-torsion specialization, nonconstant marked fibre invariant, generic Picard rank "
            "19, and generic rank decomposition 1+1 are exact."
        ),
        "not_proved": (
            "The catalogue is a low-complexity rank-seven seed list, not the full Oguiso--Shioda "
            "classification.  The ansatz does not search higher-degree or rational-function "
            "twist sections.  No claim is made that rank sum two is globally maximal among all "
            "quadratic base changes, or that every rational specialization retains exact rank two."
        ),
    },
}


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()
output_path = arguments.output.resolve()
encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if not output_path.exists() or output_path.read_text() != encoded:
        raise SystemExit(f"stale artifact: {output_path}")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(encoded)

print(
    "RESSQBC|catalogue=5|controls=Golay,NS0031|"
    "new=E6A1-clean-rminus1,E6A1-secondary-r2|"
    "winner=rho19_MW2_heightdet1|status=PASS",
    flush=True,
)
print(f"OUTPUT|{output_path}", flush=True)
