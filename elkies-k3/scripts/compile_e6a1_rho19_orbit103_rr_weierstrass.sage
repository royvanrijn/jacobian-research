#!/usr/bin/env sage-python
"""Compile the orbit-103 E6+A1 neighbor to an exact Weierstrass model.

status: ACTIVE_PROOF
claim: resolved RR pencil for P0+P1+A3_2 and its 2D5+2A2/MW3 Jacobian
inputs: elkies-k3-e6a1-rho19-genuine-q2-neighbors-v1.json
outputs: elkies-k3-e6a1-rho19-orbit103-rr-weierstrass-v1.json

The old generic-fibre divisor is P0+P1.  A normalized secant function has
those two simple horizontal poles.  Its only vertical pole is the middle
component of the old I4 fibre, as shown by exact central/side resolution
charts at infinity.  The certified nef-isotropic class has h0=2, so this
function together with 1 is the complete resolved Riemann--Roch pencil.
Elimination gives a binary quartic with rational points at t=+/-1 and hence
the explicit Jacobian elliptic K3 recorded below.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import (
    EllipticCurve,
    PolynomialRing,
    QQ,
    ZZ,
    factor,
    matrix,
    vector,
)


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
GEN = ROOT / "artifacts/generated-results"
SOURCE = GEN / "elkies-k3-e6a1-rho19-genuine-q2-neighbors-v1.json"
DEFAULT_OUTPUT = GEN / "elkies-k3-e6a1-rho19-orbit103-rr-weierstrass-v1.json"

_compiler_path = HERE / "elliptic_neighbor_compiler.sage"
exec(compile(_compiler_path.read_text(), str(_compiler_path), "exec"), globals())


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def poly_coefficients(polynomial):
    polynomial = polynomial.parent()(polynomial)
    return [str(polynomial[i]) for i in range(polynomial.degree() + 1)]


def rational_leading_term_at_zero(value, uniformizer):
    """Return (order, leading coefficient) in a univariate local parameter."""
    field = uniformizer.parent().fraction_field()
    value = field(value)
    numerator = uniformizer.parent()(value.numerator())
    denominator = uniformizer.parent()(value.denominator())
    numerator_order = numerator.valuation()
    denominator_order = denominator.valuation()
    numerator_lead = (numerator / uniformizer**numerator_order)(0)
    denominator_lead = (denominator / uniformizer**denominator_order)(0)
    return numerator_order - denominator_order, numerator_lead / denominator_lead


if not SOURCE.exists():
    raise FileNotFoundError(SOURCE)
source = json.loads(SOURCE.read_text())
if source.get("status") != (
    "PASS_EXACT_COMPLETE_GENUINE_Q2_CENSUS_AND_18_NEF_MW3_FRAMES"
):
    raise ArithmeticError("genuine q=2 source certificate is not exact")
target = source["preferred_equation_compiler_target"]
if target["orbit"] != 103 or target["root_type"] != "2A2+2D5":
    raise ArithmeticError("preferred orbit-103 target changed")
expected_divisor = [0] * 19
expected_divisor[15] = expected_divisor[17] = expected_divisor[18] = 1
if target["divisor_in_ns_basis"] != expected_divisor:
    raise ArithmeticError("orbit 103 is no longer P0+P1+A3_2")
orbit_record = next(
    item for item in source["nef_rank_three_neighbors"] if item["orbit"] == 103
)
if not orbit_record["nef"] or orbit_record["physical_reduction"]["sequence"]:
    raise ArithmeticError("orbit 103 lost its nef empty-reduction certificate")


# Old K3 and marked sections over QQ(k)(t).
K_RING = PolynomialRing(QQ, "k")
k_polynomial = K_RING.gen()
K = K_RING.fraction_field()
T_RING = PolynomialRing(K, "t")
t = T_RING.gen()
FT = T_RING.fraction_field()
k = K(k_polynomial)

D_parameter = 3 * k**2 - 4
c = 2 * k / D_parameter
lam = -(k**2 - 4) * D_parameter / 4
p = 4 / D_parameter**2
q = k**2 / D_parameter
s_coefficient = 8 / D_parameter**3
v = 2 * (k**2 + 2) / D_parameter**2
H = 1 - t**2

old_a = H**3 * (lam - 3 * H)
old_b = H**4 * (c**2 * lam**2 + lam * H - 2 * H**2)
old_curve = EllipticCurve(FT, [old_a, old_b])
x0 = -H**2
y0 = c * lam * H**2
x1 = p * lam**2 + q * lam * H - H**2
y1 = lam * t * (s_coefficient * lam**2 + v * lam * H)
P0 = old_curve(x0, y0)
P1 = old_curve(x1, y1)
for label, point in (("P0", P0), ("P1", P1)):
    if point[1] ** 2 != point[0] ** 3 + old_a * point[0] + old_b:
        raise ArithmeticError(f"old section {label} changed")


# The secant line through P0,P1 has a coefficient pole at t=2/k.  Multiplying
# the degree-two function by L=t-2/k cancels it.  With affine old coordinates
# x,y, the second RR generator is
#
# z = L*(y+y0+m*(x-x0))/((x-x0)*(x-x1)),  m=(y1-y0)/(x1-x0).
L = t - 2 / k
delta_x = T_RING(x1 - x0)
delta_y = T_RING(y1 - y0)
mL = T_RING(L * delta_y / delta_x)
expected_mL = T_RING(
    2 / k
    * (
        t**3
        + (-k**3 / 4 + k / 2) * t**2
        + (k**2 / 2 - 3) * t
        + k / 2
    )
)
if mL != expected_mL:
    raise ArithmeticError("normalized secant coefficient changed")
if delta_x != T_RING(k**2 * (k**2 - 4) / 4 * (t - 2 / k) * (t + 2 / k)):
    raise ArithmeticError("P0/P1 x-collision divisor changed")

# At both IV* fibres, P0 has orders (2,2), P1 is a smooth affine point, and
# L,m are units.  For every exceptional divisorial valuation of
# y^2=x^3+h^3*a*x+h^4*b, the tropical equality gives
# min(v(y),2v(h)) >= min(v(x),2v(h)).  Hence the displayed secant quotient is
# regular on all six E6 exceptional components.  The exact hypotheses are
# checked here at both h=t-1 and h=t+1.
ivstar_inputs = []
for label, place in (("t=1", t - 1), ("t=-1", t + 1)):
    record = {
        "place": label,
        "ord_a": int(T_RING(old_a).valuation(place)),
        "ord_b": int(T_RING(old_b).valuation(place)),
        "ord_x_P0": int(T_RING(x0).valuation(place)),
        "ord_y_P0": int(T_RING(y0).valuation(place)),
        "ord_x_P1": int(T_RING(x1).valuation(place)),
        "ord_y_P1": int(T_RING(y1).valuation(place)),
        "ord_L": int(T_RING(L).valuation(place)),
        "ord_mL": int(T_RING(mL).valuation(place)),
    }
    if tuple(record[key] for key in ("ord_a", "ord_b")) != (3, 4):
        raise ArithmeticError("old IV* coefficient orders changed")
    if tuple(record[key] for key in ("ord_x_P0", "ord_y_P0")) != (2, 2):
        raise ArithmeticError("P0 IV* entrance orders changed")
    if tuple(record[key] for key in ("ord_x_P1", "ord_y_P1")) != (0, 0):
        raise ArithmeticError("P1 is no longer smooth at IV*")
    if record["ord_L"] != 0 or record["ord_mL"] != 0:
        raise ArithmeticError("normalized secant is not an IV* unit")
    ivstar_inputs.append(record)


# Resolve the old I4 fibre at infinity.  Put s=1/t and use minimal coordinates
# x_s=s^4*x, y_s=s^6*y.  The node is (-1,0).  The side charts have
# x_s=-1+s*xi, y_s=s*eta; the middle chart has s^2 instead.  Exact leading
# terms prove that z is regular on the two side components and has a simple
# pole on A3_2.
XI_ETA_RING = PolynomialRing(K, names=("xi", "eta"))
xi_polynomial, eta_polynomial = XI_ETA_RING.gens()
XI_ETA_FIELD = XI_ETA_RING.fraction_field()
S_RING = PolynomialRing(XI_ETA_FIELD, "s")
s = S_RING.gen()
FS = S_RING.fraction_field()
xi = XI_ETA_FIELD(xi_polynomial)
eta = XI_ETA_FIELD(eta_polynomial)
k_local = XI_ETA_FIELD(k)
D_local = 3 * k_local**2 - 4
c_local = 2 * k_local / D_local
lam_local = -(k_local**2 - 4) * D_local / 4
p_local = 4 / D_local**2
q_local = k_local**2 / D_local
sc_local = 8 / D_local**3
v_local = 2 * (k_local**2 + 2) / D_local**2
Hm = s**2 - 1
infinity_a = Hm**3 * (lam_local * s**2 - 3 * Hm)
infinity_b = Hm**4 * (
    c_local**2 * lam_local**2 * s**4
    + lam_local * Hm * s**2
    - 2 * Hm**2
)
x0s = -Hm**2
y0s = c_local * lam_local * s**2 * Hm**2
x1s = p_local * lam_local**2 * s**4 + q_local * lam_local * Hm * s**2 - Hm**2
y1s = lam_local * (sc_local * lam_local**2 * s**5 + v_local * lam_local * Hm * s**3)
t_local = 1 / s
mL_local = 2 / k_local * (
    t_local**3
    + (-k_local**3 / 4 + k_local / 2) * t_local**2
    + (k_local**2 / 2 - 3) * t_local
    + k_local / 2
)
ms = FS(s**3 * mL_local / (1 - 2 * s / k_local))

infinity_charts = {}
for label, exponent in (("side_A3_1_or_A3_3", 1), ("middle_A3_2", 2)):
    xs = -1 + s**exponent * xi
    ys = s**exponent * eta
    equation = ys**2 - xs**3 - infinity_a * xs - infinity_b
    numerator = ys + y0s + ms * (xs - x0s)
    denominator = (xs - x0s) * (xs - x1s)
    z_local = FS(s * (1 - 2 * s / k_local) * numerator / denominator)
    equation_order, equation_lead = rational_leading_term_at_zero(equation, s)
    z_order, z_lead = rational_leading_term_at_zero(z_local, s)
    infinity_charts[label] = {
        "equation_order": int(equation_order),
        "equation_initial_form": str(factor(equation_lead)),
        "z_order": int(z_order),
        "z_initial_form": str(factor(z_lead)),
    }
if infinity_charts["side_A3_1_or_A3_3"]["z_order"] != 0:
    raise ArithmeticError("secant acquired a pole on a side I4 component")
if infinity_charts["middle_A3_2"]["z_order"] != -1:
    raise ArithmeticError("secant lost its simple A3_2 pole")

resolved_rr = {
    "divisor": "P0+P1+A3_2",
    "h0": 2,
    "basis": ["1", "z"],
    "generic_horizontal_poles": ["P0 with multiplicity one", "P1 with multiplicity one"],
    "finite_base_regularization": "multiply the secant by L=t-2/k",
    "ivstar_tropical_inputs": ivstar_inputs,
    "infinity_resolved_charts": infinity_charts,
    "completeness_reason": (
        "The source certificate proves D primitive nef isotropic, hence h0(D)=2; "
        "the two displayed independent global sections therefore form all of H0(D)."
    ),
}


# Eliminate x after imposing the pencil parameter z.  Clearing L gives
# y=(z*(x-x0)*(x-x1)-L*line)/L.  The resulting quartic polynomial in x is
# divisible by the two known pole factors, leaving a quadratic.
Z_RING = PolynomialRing(K, "z")
z_polynomial = Z_RING.gen()
KZ = Z_RING.fraction_field()
TZ_RING = PolynomialRing(KZ, "t")
tz = TZ_RING.gen()
XZ_RING = PolynomialRing(TZ_RING, "x")
x = XZ_RING.gen()
kz = KZ(k)
Dz = 3 * kz**2 - 4
cz = 2 * kz / Dz
lamz = -(kz**2 - 4) * Dz / 4
pz = 4 / Dz**2
qz = kz**2 / Dz
scz = 8 / Dz**3
vz = 2 * (kz**2 + 2) / Dz**2
Hz = 1 - tz**2
az = Hz**3 * (lamz - 3 * Hz)
bz = Hz**4 * (cz**2 * lamz**2 + lamz * Hz - 2 * Hz**2)
x0z = -Hz**2
y0z = cz * lamz * Hz**2
x1z = pz * lamz**2 + qz * lamz * Hz - Hz**2
y1z = lamz * tz * (scz * lamz**2 + vz * lamz * Hz)
Lz = tz - 2 / kz
mLz = TZ_RING(
    2 / kz
    * (
        tz**3
        + (-kz**3 / 4 + kz / 2) * tz**2
        + (kz**2 / 2 - 3) * tz
        + kz / 2
    )
)
pole_denominator = (x - x0z) * (x - x1z)
scaled_line = Lz * y0z + mLz * (x - x0z)
cleared = (
    (z_polynomial * pole_denominator - scaled_line) ** 2
    - Lz**2 * (x**3 + az * x + bz)
)
quadratic, remainder = cleared.quo_rem(pole_denominator)
if remainder or quadratic.degree() != 2:
    raise ArithmeticError("resolved pencil did not eliminate to a quadratic in x")
radicand = TZ_RING(quadratic.discriminant())
quartic_z, square_factor = squarefree_binary_quartic(radicand, TZ_RING)
if square_factor**2 != Lz**2 or quartic_z.degree() != 4:
    raise ArithmeticError("pencil discriminant did not reduce to the expected quartic")


# Put the two D5 fibres symmetrically at r=+/-2.
# r = k(k^2-4)z/2-k, so z=2(r+k)/(k(k^2-4)).
R_RING = PolynomialRing(K, "r")
r = R_RING.gen()
FR = R_RING.fraction_field()
TR2_RING = PolynomialRing(FR, "t")
t2 = TR2_RING.gen()
z_substitution = 2 * (r + k) / (k * (k**2 - 4))
quartic_r = TR2_RING(
    sum(FR(quartic_z[index](z=z_substitution)) * t2**index for index in range(5))
)

w_plus = FR((r + 2) ** 2 / (k * (k - 2)))
w_minus = FR((r - 2) ** 2 / (k * (k + 2)))
if quartic_r(1) != w_plus**2 or quartic_r(-1) != w_minus**2:
    raise ArithmeticError("rational IV* points on the quartic changed")

standard_a_z, standard_b_z, unused_quartic_discriminant = (
    binary_quartic_jacobian_coefficients(quartic_z)
)
minimalizing_unit = k * (k**2 - 4)

A4 = (
    k**4 * r**4
    + 16 * k * (k**2 + 12) * r**3
    + 8 * k**2 * (k**2 + 152) * r**2
    + 64 * k * (31 * k**2 - 12) * r
    + 16 * k**2 * (61 * k**2 - 48)
)
B6 = (
    k**6 * r**6
    + 24 * k**3 * (k**2 + 12) * r**5
    + 12 * (k**6 + 160 * k**4 + 192 * k**2 + 1152) * r**4
    + k * (3072 * k**4 + 35072 * k**2 + 27648) * r**3
    + (1488 * k**6 + 91776 * k**4 + 4608 * k**2 - 55296) * r**2
    + k * (85632 * k**4 - 13824 * k**2 - 110592) * r
    + k**2 * (26560 * k**4 - 4608 * k**2 - 55296)
)
child_a = R_RING(-27 * (r**2 - 4) ** 2 * A4)
child_b = R_RING(54 * (r**2 - 4) ** 3 * B6)
if FR(minimalizing_unit**4 * standard_a_z(z=z_substitution)) != FR(child_a):
    raise ArithmeticError("clean child A coefficient changed")
if FR(minimalizing_unit**6 * standard_b_z(z=z_substitution)) != FR(child_b):
    raise ArithmeticError("clean child B coefficient changed")

child = EllipticCurve(FR, [child_a, child_b])

# Two explicit rational points obtained from the binary-quartic covariants.
X_plus = R_RING(
    3
    * (r + 2)
    * (
        k**2 * r**3
        + (-2 * k**2 - 16 * k + 96) * r**2
        + (-20 * k**2 + 320 * k + 192) * r
        + 232 * k**2
        + 192 * k
    )
)
Y_plus = R_RING(
    216
    * (r + k)
    * (r + 2) ** 2
    * (
        (k - 2) ** 2 * r**2
        + (-12 * k**2 + 32 * k + 48) * r
        + 52 * k**2
        + 80 * k
        + 16
    )
)
X_minus = R_RING(
    3
    * (r - 2)
    * (
        k**2 * r**3
        + (2 * k**2 - 16 * k - 96) * r**2
        + (-20 * k**2 - 320 * k + 192) * r
        - 232 * k**2
        + 192 * k
    )
)
Y_minus = R_RING(
    -216
    * (r + k)
    * (r - 2) ** 2
    * (
        (k + 2) ** 2 * r**2
        + (12 * k**2 + 32 * k - 48) * r
        + 52 * k**2
        - 80 * k
        + 16
    )
)
Q_plus = child(X_plus, Y_plus)
Q_minus = child(X_minus, Y_minus)
for label, point in (("Q_plus", Q_plus), ("Q_minus", Q_minus)):
    if point[1] ** 2 != point[0] ** 3 + child_a * point[0] + child_b:
        raise ArithmeticError(f"transported child point {label} changed")


# Exact fibre classification.  The reduced cubic discriminant suffices for
# factor orders; the full elliptic discriminant differs by -16.
reduced_discriminant = R_RING(4 * child_a**3 + 27 * child_b**2)
constant = 2**11 * 3**13 * k**3 * (k**2 - 4) ** 2 * (k**2 + QQ(4) / 3)
residual_quartic = R_RING(
    reduced_discriminant
    / (constant * (r + k) ** 3 * (r**2 - 4) ** 7)
)
if residual_quartic.degree() != 4:
    raise ArithmeticError("child residual I1 polynomial changed")
if any(
    residual_quartic.gcd(factor_value).degree() > 0
    for factor_value in (r + k, r - 2, r + 2)
):
    raise ArithmeticError("generic child fibre factors collided")

classification = classify_finite_short_weierstrass_fibres(
    R_RING, child_a, child_b
)
finite_symbols = sorted(
    (item["kodaira"], item["degree"]) for item in classification["finite_fibres"]
)
if finite_symbols != [("I1", 4), ("I1*", 1), ("I1*", 1), ("I3", 1)]:
    raise ArithmeticError(f"child finite fibre profile changed: {finite_symbols}")
if classification["infinity_boundary"]["normalized_orders"] != (0, 0, 3):
    raise ArithmeticError("child infinity fibre is no longer I3")
infinity_kodaira = kodaira_data_from_short_orders(0, 0, 3)
if infinity_kodaira != (2, 3, 3, "I3"):
    raise ArithmeticError("infinity Kodaira classification changed")

finite_root_rank = classification["finite_root_rank"]
finite_euler = classification["finite_euler_number"]
finite_root_determinant = classification["finite_root_determinant"]
if (finite_root_rank, finite_euler, finite_root_determinant) != (12, 21, 48):
    raise ArithmeticError("finite child root data changed")
if (finite_root_rank + 2, finite_euler + 3, finite_root_determinant * 3) != (
    14,
    24,
    144,
):
    raise ArithmeticError("global child root/Euler data changed")


payload = {
    "schema": "elkies-k3.e6a1-rho19-orbit103-rr-weierstrass.v1",
    "status": "PASS_EXACT_RESOLVED_RR_QUARTIC_AND_WEIERSTRASS",
    "inputs": {relative(SOURCE): digest(SOURCE)},
    "parameter_exclusions": [
        "k=0",
        "k=2",
        "k=-2",
        "3*k^2=4",
        "3*k^2+4=0",
        "zeros of the displayed residual-quartic discriminant",
    ],
    "source_k3": {
        "equation": "Y^2=X^3+H^3*(lambda-3H)*X+H^4*(c^2*lambda^2+lambda*H-2H^2)",
        "H": "1-t^2",
        "D": "3*k^2-4",
        "c": "2*k/D",
        "lambda": "-(k^2-4)*D/4",
        "P0": [str(x0), str(y0)],
        "P1": [str(x1), str(y1)],
    },
    "resolved_riemann_roch": resolved_rr,
    "pencil_parameter": {
        "z": "(t-2/k)*(y+y0+m*(x-x0))/((x-x0)*(x-x1))",
        "m": "(y1-y0)/(x1-x0)",
        "normalized_secant_L_times_m": str(mL),
        "new_base_normalization": "r=k*(k^2-4)*z/2-k",
    },
    "elimination": {
        "quadratic_in_old_x": [str(quadratic[index]) for index in range(3)],
        "radicand_square_factor": str(square_factor),
        "binary_quartic_in_t_coefficients_after_r_normalization": poly_coefficients(quartic_r),
        "rational_points": [
            {"t": "1", "w": str(w_plus)},
            {"t": "-1", "w": str(w_minus)},
        ],
        "jacobian_identification": (
            "The binary quartic has a QQ(k,r)-point, so the genus-one pencil "
            "is isomorphic to its displayed Jacobian after choosing that point as origin."
        ),
    },
    "child_weierstrass": {
        "equation": "y^2=x^3-27*(r^2-4)^2*A4*x+54*(r^2-4)^3*B6",
        "A4_coefficients_low_to_high": poly_coefficients(A4),
        "B6_coefficients_low_to_high": poly_coefficients(B6),
        "short_a_coefficients_low_to_high": poly_coefficients(child_a),
        "short_b_coefficients_low_to_high": poly_coefficients(child_b),
        "explicit_rational_points": {
            "Q_plus": [str(X_plus), str(Y_plus)],
            "Q_minus": [str(X_minus), str(Y_minus)],
        },
    },
    "fibre_certificate": {
        "finite_fibres": [
            {
                **item,
                "factor": str(item["factor"]),
                "raw_orders": list(item["raw_orders"]),
                "minimal_orders": list(item["minimal_orders"]),
            }
            for item in classification["finite_fibres"]
        ],
        "infinity": {
            "orders": [0, 0, 3],
            "kodaira": "I3",
        },
        "residual_I1_quartic_coefficients_low_to_high": poly_coefficients(residual_quartic),
        "profile": "2I1*+2I3+4I1",
        "root_type": "2D5+2A2",
        "root_data": [14, 92, 144],
        "euler_number": 24,
    },
    "mordell_weil": {
        "geometric_rank": 3,
        "torsion": "trivial",
        "height_gram": target["mw_height_gram"],
        "regulator": "1/4",
        "proof": (
            "The equation realizes the certified orbit-103 pencil on the same "
            "rho=19 K3; Shioda--Tate with root rank 14 gives rank 3, and the "
            "saturated frame certificate supplies torsion and the height lattice."
        ),
        "arithmetic_field_warning": (
            "The geometric rank-three lattice is exact. The two displayed "
            "QQ(k,r)-points are explicit, but this certificate does not claim "
            "that they generate the full arithmetic Mordell--Weil group."
        ),
    },
    "proof_boundary": {
        "proved": (
            "Complete resolved H0 basis, exact elimination, rational quartic "
            "origin, clean Weierstrass equation, fibre profile, geometric MW "
            "rank three, torsion, and height-lattice isometry."
        ),
        "open": (
            "An explicit three-section arithmetic basis over QQ(k)(r), its "
            "field-of-definition decomposition, and coefficient-optimality "
            "against the A7+D7 secondary pencil are not asserted."
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
    "E6A1O103|RR=2|quartic=4|fibres=2I1*+2I3+4I1|MW=3|status=PASS_EXACT",
    flush=True,
)
print(f"OUTPUT|{output_path}", flush=True)
