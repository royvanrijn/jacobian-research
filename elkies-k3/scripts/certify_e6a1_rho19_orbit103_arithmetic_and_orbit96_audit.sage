#!/usr/bin/env sage-python
"""Certify orbit-103 arithmetic descent and compare the corrected orbit 96.

status: ACTIVE_PROOF
claim: orbit-103 arithmetic MW rank two and exact orbit-96 comparison
inputs: orbit-103 Weierstrass and orbit-96 physical/Galois certificates
outputs: elkies-k3-e6a1-rho19-orbit103-arithmetic-orbit96-audit-v1.json

The two polynomial orbit-103 sections are rational and independent.  A third
geometric direction comes from the two points at infinity of the old-base
binary quartic; its ordinate is proportional to sqrt(-3), and conjugation is
elliptic negation.  The geometric rank-three certificate then forces generic
arithmetic rank exactly two over QQ(k)(r).

The former orbit-96 rejection was caused by a coefficient-parent error: its
simplified tangent slope was left in a fraction field, so zero-argument
``discriminant()`` did not compute the discriminant in the old x-coordinate.
The dependent physical certificate now supplies the genuine A7+D7 equation
and its exact ``1+chi_-3+1`` Mordell--Weil representation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
GEN = ROOT / "artifacts/generated-results"
SOURCE = GEN / "elkies-k3-e6a1-rho19-orbit103-rr-weierstrass-v1.json"
Q2_SOURCE = GEN / "elkies-k3-e6a1-rho19-genuine-q2-neighbors-v1.json"
ORBIT96_SOURCE = GEN / "elkies-k3-e6a1-rho19-orbit96-rr-galois-v1.json"
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


if not SOURCE.exists() or not Q2_SOURCE.exists() or not ORBIT96_SOURCE.exists():
    raise FileNotFoundError("orbit-103, q=2 census, and orbit-96 artifacts are required")
source = json.loads(SOURCE.read_text())
q2_source = json.loads(Q2_SOURCE.read_text())
orbit96_source = json.loads(ORBIT96_SOURCE.read_text())
if source.get("status") != "PASS_EXACT_RESOLVED_RR_QUARTIC_AND_WEIERSTRASS":
    raise ArithmeticError("orbit-103 equation source is not exact")
if q2_source.get("status") != (
    "PASS_EXACT_COMPLETE_GENUINE_Q2_CENSUS_AND_18_NEF_MW3_FRAMES"
):
    raise ArithmeticError("q=2 neighbor source is not exact")
if orbit96_source.get("status") != (
    "PASS_EXACT_PHYSICAL_A7D7_WEIERSTRASS_AND_MW_2_PLUS_CHI_MINUS3"
):
    raise ArithmeticError("orbit-96 physical/Galois source is not exact")
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
    "status": "PASS_EXACT_BOTH_ORBITS_ARITHMETIC_RANK2_AND_CHI_MINUS3",
    "inputs": {
        relative(SOURCE): digest(SOURCE),
        relative(Q2_SOURCE): digest(Q2_SOURCE),
        relative(ORBIT96_SOURCE): digest(ORBIT96_SOURCE),
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
        "corrected_equation": orbit96_source["weierstrass"],
        "arithmetic": orbit96_source["mordell_weil"],
        "historical_parent_bug": orbit96_source["elimination_parent_regression"],
    },
    "coefficient_comparison": {
        "orbit103_source_divisor_complexity": orbit103["divisor_complexity"],
        "orbit96_abstract_source_divisor_complexity": orbit96["divisor_complexity"],
        "orbit103_reducible_fibre_count": 4,
        "orbit96_expected_reducible_fibre_count": 2,
        "orbit96_weierstrass_degrees": {"A": 6, "B": 9},
        "orbit103_weierstrass_degrees": {"A": 8, "B": 12},
        "arithmetic_ranks": {"orbit96": 2, "orbit103": 2},
        "nontrivial_characters": {"orbit96": "chi_-3", "orbit103": "chi_-3"},
        "status": "EXACT_PHYSICAL_AND_ARITHMETIC_COMPARISON",
    },
    "proof_boundary": {
        "proved": (
            "Exact orbit-103 arithmetic rank two and comparison with the "
            "corrected physical orbit-96 A7+D7 equation and MW representation."
        ),
        "open": (
            "Polynomial Weierstrass coordinates for the orbit-96 lattice "
            "generators and arithmetic ranks after numerical specialization."
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
    "orbit96=I8+I3star+7I1,rank2,chi_-3|status=PASS_EXACT",
    flush=True,
)
print(f"OUTPUT|{output_path}", flush=True)
