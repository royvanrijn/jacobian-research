#!/usr/bin/env sage-python
"""Certify orbit-103 arithmetic descent and audit the orbit-96 tangent trace.

status: ACTIVE_PROOF
claim: orbit-103 arithmetic MW rank two and fail-closed A7+D7 comparison gate
inputs: elkies-k3-e6a1-rho19-orbit103-rr-weierstrass-v1.json
outputs: elkies-k3-e6a1-rho19-orbit103-arithmetic-orbit96-audit-v1.json

The two polynomial orbit-103 sections are rational and independent.  A third
geometric direction comes from the two points at infinity of the old-base
binary quartic; its ordinate is proportional to sqrt(-3), and conjugation is
elliptic negation.  The geometric rank-three certificate then forces generic
arithmetic rank exactly two over QQ(k)(r).

The abstract orbit-96 record has horizontal trace 2*P0.  The obvious resolved
tangent function on the displayed equation is compiled exactly here.  Its
minimal Jacobian has the old 2IV*+I4+4I1 fingerprint, not A7+D7.  Thus the
current abstract E6 root basis is not enough to make an equation-level
coefficient comparison; an explicit physical E6 marking is a required gate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, factor


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
GEN = ROOT / "artifacts/generated-results"
SOURCE = GEN / "elkies-k3-e6a1-rho19-orbit103-rr-weierstrass-v1.json"
Q2_SOURCE = GEN / "elkies-k3-e6a1-rho19-genuine-q2-neighbors-v1.json"
DEFAULT_OUTPUT = (
    GEN / "elkies-k3-e6a1-rho19-orbit103-arithmetic-orbit96-audit-v1.json"
)

_compiler_path = HERE / "elliptic_neighbor_compiler.sage"
exec(compile(_compiler_path.read_text(), str(_compiler_path), "exec"), globals())


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial_coefficients(polynomial):
    polynomial = polynomial.parent()(polynomial)
    return [str(polynomial[index]) for index in range(polynomial.degree() + 1)]


if not SOURCE.exists() or not Q2_SOURCE.exists():
    raise FileNotFoundError("orbit-103 equation and q=2 census artifacts are required")
source = json.loads(SOURCE.read_text())
q2_source = json.loads(Q2_SOURCE.read_text())
if source.get("status") != "PASS_EXACT_RESOLVED_RR_QUARTIC_AND_WEIERSTRASS":
    raise ArithmeticError("orbit-103 equation source is not exact")
if q2_source.get("status") != (
    "PASS_EXACT_COMPLETE_GENUINE_Q2_CENSUS_AND_18_NEF_MW3_FRAMES"
):
    raise ArithmeticError("q=2 neighbor source is not exact")
if source["mordell_weil"]["geometric_rank"] != 3:
    raise ArithmeticError("orbit-103 geometric rank changed")
if source["mordell_weil"]["torsion"] != "trivial":
    raise ArithmeticError("orbit-103 torsion changed")


# Reconstruct the clean orbit-103 equation and its two rational sections.
K_RING = PolynomialRing(QQ, "k")
k_polynomial = K_RING.gen()
K = K_RING.fraction_field()
R_RING = PolynomialRing(K, "r")
r = R_RING.gen()
FR = R_RING.fraction_field()
k = FR(k_polynomial)

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
child = EllipticCurve(FR, [child_a, child_b])

X_plus = R_RING(
    3 * (r + 2) * (
        k**2 * r**3
        + (-2 * k**2 - 16 * k + 96) * r**2
        + (-20 * k**2 + 320 * k + 192) * r
        + 232 * k**2 + 192 * k
    )
)
Y_plus = R_RING(
    216 * (r + k) * (r + 2) ** 2 * (
        (k - 2) ** 2 * r**2
        + (-12 * k**2 + 32 * k + 48) * r
        + 52 * k**2 + 80 * k + 16
    )
)
X_minus = R_RING(
    3 * (r - 2) * (
        k**2 * r**3
        + (2 * k**2 - 16 * k - 96) * r**2
        + (-20 * k**2 - 320 * k + 192) * r
        - 232 * k**2 + 192 * k
    )
)
Y_minus = R_RING(
    -216 * (r + k) * (r - 2) ** 2 * (
        (k + 2) ** 2 * r**2
        + (12 * k**2 + 32 * k - 48) * r
        + 52 * k**2 - 80 * k + 16
    )
)
Q_plus = child(X_plus, Y_plus)
Q_minus = child(X_minus, Y_minus)
if Q_plus.is_zero() or Q_minus.is_zero() or Q_plus == Q_minus or Q_plus == -Q_minus:
    raise ArithmeticError("the two rational orbit-103 directions lost independence data")

# Both rational sections have P.O=0.  Each meets one I1* spinor component and
# the nonidentity component of each I3.  The standard local corrections are
# 5/4, 2/3, 2/3, hence height 4-5/4-2/3-2/3=17/12.
if (
    X_minus.valuation(r - 2), Y_minus.valuation(r - 2),
    X_plus.valuation(r + 2), Y_plus.valuation(r + 2),
) != (1, 2, 1, 2):
    raise ArithmeticError("I1* spinor entrances of Q+/- changed")
if X_plus(r=2) == 0 or X_minus(r=-2) == 0:
    raise ArithmeticError("Q+/- left the identity component at the other I1*")
if X_plus.degree() != 4 or X_minus.degree() != 4:
    raise ArithmeticError("Q+/- acquired an intersection with the old zero")
if Y_plus.valuation(r + k) != 1 or Y_minus.valuation(r + k) != 1:
    raise ArithmeticError("Q+/- lost the finite I3 node entrance")
finite_i3_node_x = 3 * k**2 * (k**2 - 4) ** 2
if X_plus(r=-k) != finite_i3_node_x or X_minus(r=-k) != finite_i3_node_x:
    raise ArithmeticError("Q+/- no longer pass through the finite I3 node")
if X_plus[4] != 3 * k**2 or X_minus[4] != 3 * k**2:
    raise ArithmeticError("Q+/- lost the infinity-I3 node entrance")
rational_section_height = QQ(4) - QQ(5) / 4 - 2 * QQ(2) / 3
if rational_section_height != QQ(17) / 12:
    raise ArithmeticError("Q+/- height calculation changed")


# The old-base quartic has leading coefficient
# -48*(r+k)^2/(k^2*(k^2-4)^2).  Its two infinity points are defined over
# j^2=-3.  The covariant image on the clean child is the following point.
X_delta = R_RING(
    (-6 * k**2 - 12) * r**4
    - 96 * k * r**3
    + (-48 * k**2 - 672) * r**2
    - 1152 * k * r
    - 480 * k**2 - 192
)
Y_delta_coefficient = R_RING(
    (-18 * k**2 - 24) * r**6
    + (-36 * k**3 + 144 * k) * r**5
    + (648 * k**2 + 3168) * r**4
    + (672 * k**3 + 12672 * k) * r**3
    + (14112 * k**2 + 12672) * r**2
    + (4032 * k**3 + 20736 * k) * r
    + 8064 * k**2 - 1536
)
if R_RING(-3 * Y_delta_coefficient**2) != R_RING(
    X_delta**3 + child_a * X_delta + child_b
):
    raise ArithmeticError("the sqrt(-3) orbit-103 direction left the child")
if not X_delta or not Y_delta_coefficient:
    raise ArithmeticError("the anti-invariant direction became zero")

# Conjugation j -> -j fixes X_delta and negates its ordinate, hence sends the
# section to its elliptic inverse.  It supplies a nonzero anti-invariant line.
# The geometric rank is three.  The two rational points are independent: they
# have the same positive height, and dependence in the torsion-free MW group
# would force Q_plus=+/-Q_minus, already excluded above.  Therefore the fixed
# subspace has dimension exactly two and the generic arithmetic rank is two.
geometric_rank = 3
rational_independent_rank = 2
anti_invariant_rank_lower_bound = 1
arithmetic_rank = geometric_rank - anti_invariant_rank_lower_bound
if arithmetic_rank != rational_independent_rank:
    raise ArithmeticError("orbit-103 Galois rank decomposition changed")


# Audit the tempting orbit-96 tangent trace.  Orbit 96 has abstract horizontal
# trace 2*P0, but its E6 root coordinates have not been matched to resolved
# equation components.  Compile the unique generic-fibre tangent function
# (both signs give the same Jacobian) with the apparent t=-1 regularizer.
T_RING = PolynomialRing(K, "t")
t = T_RING.gen()
FT = T_RING.fraction_field()
k_t = FT(k_polynomial)
D_parameter = 3 * k_t**2 - 4
c = 2 * k_t / D_parameter
lam = -(k_t**2 - 4) * D_parameter / 4
H = 1 - t**2
old_a = H**3 * (lam - 3 * H)
old_b = H**4 * (c**2 * lam**2 + lam * H - 2 * H**2)
x0 = -H**2
y0 = c * lam * H**2
tangent_slope = FT((3 * x0**2 + old_a) / (2 * y0))
if tangent_slope != ((-3 * k_t**2 / 4 + 1) / k_t) * (t - 1) * (t + 1):
    raise ArithmeticError("P0 tangent slope changed")

Z_RING = PolynomialRing(K, "z")
z = Z_RING.gen()
KZ = Z_RING.fraction_field()
TZ_RING = PolynomialRing(KZ, "t")
tz = TZ_RING.gen()
XZ_RING = PolynomialRing(TZ_RING, "x")
x = XZ_RING.gen()
kz = KZ(k_polynomial)
Dz = 3 * kz**2 - 4
cz = 2 * kz / Dz
lamz = -(kz**2 - 4) * Dz / 4
Hz = 1 - tz**2
az = Hz**3 * (lamz - 3 * Hz)
bz = Hz**4 * (cz**2 * lamz**2 + lamz * Hz - 2 * Hz**2)
x0z = -Hz**2
y0z = cz * lamz * Hz**2
mz = (3 * x0z**2 + az) / (2 * y0z)
regularizer = tz + 1
pole_factor = x - x0z

standard_models = []
for sign in (1, -1):
    # sign=+1 has poles at +P0; sign=-1 has poles at -P0.
    scaled_line = sign * regularizer * (y0z + mz * pole_factor)
    cleared = (
        (z * pole_factor**2 - scaled_line) ** 2
        - regularizer**2 * (x**3 + az * x + bz)
    )
    quadratic, remainder = cleared.quo_rem(pole_factor**2)
    if remainder or quadratic.degree() != 2:
        raise ArithmeticError("orbit-96 tangent trace did not leave a quadratic")
    quartic, square_factor = squarefree_binary_quartic(
        TZ_RING(quadratic.discriminant()), TZ_RING
    )
    if quartic.degree() != 3 or square_factor**2 != tz**2:
        raise ArithmeticError("orbit-96 tangent trace quartic changed")
    standard_models.append(binary_quartic_jacobian_coefficients(quartic)[:2])
if standard_models[0] != standard_models[1]:
    raise ArithmeticError("the two tangent orientations lost their common Jacobian")

standard_a, standard_b = standard_models[0]
minimalizing_unit = KZ(z / (6 * (z + 1)))
tangent_a = Z_RING(minimalizing_unit**4 * standard_a)
tangent_b = Z_RING(minimalizing_unit**6 * standard_b)
C2 = Z_RING(z**2 - kz**4 / 4 + 4 * kz**2 / 3 - QQ(7) / 3)
C4 = Z_RING(
    z**4
    + (-3 * kz**4 / 8 + 2 * kz**2 - 4) * z**2
    - kz**6 / 8 + 11 * kz**4 / 8 - 4 * kz**2 + 3
)
expected_tangent_a = Z_RING(-3 * (z**2 - 1) ** 3 * C2)
expected_tangent_b = Z_RING(-2 * (z**2 - 1) ** 4 * C4)
if tangent_a != expected_tangent_a or tangent_b != expected_tangent_b:
    raise ArithmeticError("clean tangent-trace equation changed")

tangent_classification = classify_finite_short_weierstrass_fibres(
    Z_RING, tangent_a, tangent_b
)
finite_profile = sorted(
    (item["kodaira"], item["degree"])
    for item in tangent_classification["finite_fibres"]
)
if finite_profile != [("I1", 4), ("IV*", 1), ("IV*", 1)]:
    raise ArithmeticError(f"unexpected tangent-trace fibres: {finite_profile}")
if tangent_classification["infinity_boundary"]["normalized_orders"] != (0, 0, 4):
    raise ArithmeticError("tangent trace lost its I4 fibre at infinity")
if (
    tangent_classification["finite_root_rank"],
    tangent_classification["finite_euler_number"],
    tangent_classification["finite_root_determinant"],
) != (12, 20, 9):
    raise ArithmeticError("tangent trace is no longer the 2E6+A3 frame")

orbit96 = q2_source["secondary_fibre_simple_compiler_target"]
if orbit96["orbit"] != 96 or orbit96["root_type"] != "A7+D7":
    raise ArithmeticError("secondary abstract target changed")
if orbit96["divisor_complexity"] != {"max_abs": 2, "l1": 12}:
    raise ArithmeticError("orbit-96 source complexity changed")
orbit103 = q2_source["preferred_equation_compiler_target"]
if orbit103["divisor_complexity"] != {"max_abs": 1, "l1": 3}:
    raise ArithmeticError("orbit-103 source complexity changed")


payload = {
    "schema": "elkies-k3.e6a1-rho19-orbit103-arithmetic-orbit96-audit.v1",
    "status": "PASS_EXACT_ARITHMETIC_RANK2_AND_FAIL_CLOSED_A7D7_GATE",
    "inputs": {
        relative(SOURCE): digest(SOURCE),
        relative(Q2_SOURCE): digest(Q2_SOURCE),
    },
    "orbit103_arithmetic": {
        "base_field": "QQ(k)(r)",
        "geometric_rank": geometric_rank,
        "arithmetic_rank": arithmetic_rank,
        "rational_independent_points": ["Q_plus", "Q_minus"],
        "rational_point_height": str(rational_section_height),
        "independence_proof": (
            "Q_plus and Q_minus have equal positive height 17/12. In the "
            "torsion-free geometric MW group, dependence would therefore force "
            "Q_plus=+/-Q_minus, contradicted by their exact coordinates."
        ),
        "third_geometric_direction": {
            "field": "QQ(k,sqrt(-3))(r)",
            "X_coefficients_low_to_high": polynomial_coefficients(X_delta),
            "Y_over_sqrt_minus3_coefficients_low_to_high": polynomial_coefficients(
                Y_delta_coefficient
            ),
            "galois_action": "sqrt(-3)->-sqrt(-3) sends Q_delta to -Q_delta",
            "source": "the two projective points at old-base infinity on the binary quartic",
        },
        "rank_decomposition": "3 geometric = 2 invariant + 1 anti-invariant",
        "conclusion": (
            "There is no third independent generator over QQ(k)(r); the third "
            "geometric direction is anti-invariant over QQ(sqrt(-3))."
        ),
    },
    "orbit96_audit": {
        "abstract_target": {
            "orbit": 96,
            "root_type": "A7+D7",
            "mw_height_gram": orbit96["mw_height_gram"],
            "source_divisor_complexity": orbit96["divisor_complexity"],
        },
        "naive_equation_tangent_trace": {
            "function": (
                "(t+1)*(y-y0-m_tan*(x-x0))/(x-x0)^2, with the opposite "
                "orientation giving the same Jacobian"
            ),
            "equation": "y^2=x^3-3*(z^2-1)^3*C2(z)*x-2*(z^2-1)^4*C4(z)",
            "C2_coefficients_low_to_high": polynomial_coefficients(C2),
            "C4_coefficients_low_to_high": polynomial_coefficients(C4),
            "fibre_profile": "2IV*+I4+4I1",
            "root_type": "2E6+A3",
            "root_data": [15, 156, 36],
            "rejection": "This is not the abstract A7+D7 orbit-96 fibration.",
        },
        "missing_gate": (
            "The abstract E6 simple-root coordinates must be identified with a "
            "complete resolved physical E6 component marking on the equation."
        ),
    },
    "coefficient_comparison": {
        "orbit103_source_divisor_complexity": orbit103["divisor_complexity"],
        "orbit96_abstract_source_divisor_complexity": orbit96["divisor_complexity"],
        "orbit103_reducible_fibre_count": 4,
        "orbit96_expected_reducible_fibre_count": 2,
        "weierstrass_degree_floor": {"A": 8, "B": 12},
        "status": "NOT_COMPARABLE_WITHOUT_A_PHYSICALLY_MARKED_A7D7_EQUATION",
        "warning": (
            "The displayed low-coefficient tangent-trace equation is a negative "
            "control, not an A7+D7 competitor. It must not be used to claim "
            "coefficient optimality."
        ),
    },
    "proof_boundary": {
        "proved": (
            "Exact orbit-103 arithmetic rank two, explicit anti-invariant third "
            "geometric direction, and exact rejection of the naive orbit-96 "
            "tangent trace by its fibre fingerprint."
        ),
        "open": (
            "A physically resolved E6 marking and a genuine orbit-96 A7+D7 "
            "Weierstrass equation are still required for coefficient optimality."
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
    "E6A1O103ARITH|rank_Qkr=2|rank_geom=3|anti=sqrt(-3)|"
    "orbit96_tangent=REJECT_2E6+A3|status=PASS_EXACT",
    flush=True,
)
print(f"OUTPUT|{output_path}", flush=True)
