#!/usr/bin/env sage
"""Compile the exact q12/o5867 smooth chord pencil over QQ.

For the certified horizontal section P with P.O=10 and fibre twist -4, use
only the global smooth chord module for D=O+P-4F.  The ambient degrees are
deg(AA)<=16 and deg(BB)<=4.  The congruence AA*X=BB*Y mod Z^2 is a direct
22-to-2 linear calculation.  The resulting pencil is compiled to its binary
quartic and minimal Jacobian without resolved vertical rows, elimination, or
a Groebner basis.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import PolynomialRing, QQ, ZZ, inverse_mod, matrix


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
TARGET = LOCAL / "q12o5867-target-horizontal-qq.json"
SURFACE = LOCAL / "q4o164-q8o376-smooth-rr-qq.json"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=LOCAL / "q12o5867-smooth-rr-qq.json",
    )
    return parser.parse_args()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficients(poly):
    return [str(value) for value in poly.list()]


def rational_bits(values):
    answer = 0
    for value in values:
        value = QQ(value)
        answer = max(
            answer,
            abs(ZZ(value.numerator())).nbits(),
            ZZ(value.denominator()).nbits(),
        )
    return int(answer)


args = parse_args()
started = time.monotonic()
target_artifact = json.loads(TARGET.read_text(encoding="utf-8"))
surface_artifact = json.loads(SURFACE.read_text(encoding="utf-8"))
assert target_artifact["status"] == "PASS_EXACT_QQ_Q12O5867_TARGET_HORIZONTAL_SECTION"
assert target_artifact["target"]["P_dot_O"] == 10
assert target_artifact["target"]["height"] == "22"
assert target_artifact["target"]["equation_component_profile"] == [1, 1, 1, 1]

R = PolynomialRing(QQ, "v")
v = R.gen()
K = R.fraction_field()
A = R(surface_artifact["child"]["minimal_A_coefficients_low_to_high"])
B = R(surface_artifact["child"]["minimal_B_coefficients_low_to_high"])
target = target_artifact["target"]
X = R(target["x_numerator_coefficients_low_to_high"])
denominator_x = R(target["x_denominator_coefficients_low_to_high"])
Y = R(target["y_numerator_coefficients_low_to_high"])
denominator_y = R(target["y_denominator_coefficients_low_to_high"])
assert (X.degree(), denominator_x.degree()) == (24, 20)
assert (Y.degree(), denominator_y.degree()) == (35, 30)

# Recover the common monic denominator Z exactly.
Z = R.one()
for factor, exponent in denominator_x.factor():
    assert int(exponent) % 2 == 0
    Z *= factor.monic() ** (int(exponent) // 2)
assert Z.degree() == 10
assert Z ** 2 == denominator_x and Z ** 3 == denominator_y
assert Y ** 2 == X ** 3 + A * X * Z ** 4 + B * Z ** 6

# D=O+P-4F: deg(AA)<=2(P.O)-4=16 and deg(BB)<=P.O-2-4=4.
aa_degree = 16
bb_degree = 4
ambient = [
    (v ** degree, R.zero()) for degree in range(aa_degree + 1)
] + [
    (R.zero(), v ** degree) for degree in range(bb_degree + 1)
]
assert len(ambient) == 22
collision_modulus = Z ** 2
assert collision_modulus.degree() == 20
remainders = [R((AA * X - BB * Y) % collision_modulus) for AA, BB in ambient]
condition_matrix = matrix(QQ, [
    [remainder[degree] for remainder in remainders]
    for degree in range(collision_modulus.degree())
])
assert condition_matrix.nrows() == 20 and condition_matrix.ncols() == 22
assert condition_matrix.rank() == 20
kernel = condition_matrix.right_kernel().basis_matrix()
assert kernel.nrows() == kernel.rank() == 2

pairs = []
for row in kernel.rows():
    AA = sum((row[index] * ambient[index][0] for index in range(22)), R.zero())
    BB = sum((row[index] * ambient[index][1] for index in range(22)), R.zero())
    assert AA.degree() <= aa_degree and BB.degree() <= bb_degree
    assert (AA * X - BB * Y) % collision_modulus == 0
    pairs.append((AA, BB))
AA0, BB0 = pairs[0]
AA1, BB1 = pairs[1]

# Compile the degree-two chord directly over QQ(u).
UQ = PolynomialRing(QQ, "u")
u = UQ.gen()
KU = UQ.fraction_field()
VU = PolynomialRing(KU, "v")
vv = VU.gen()


def lift(poly):
    return VU([KU(value) for value in R(poly).list()])


aa = lift(AA0) + KU(u) * lift(AA1)
bb = lift(BB0) + KU(u) * lift(BB1)
X_u, Y_u, Z_u, A_u = map(lift, (X, Y, Z, A))
raw = (
    aa ** 4 - 6 * X_u * aa ** 2 * bb ** 2 + 8 * Y_u * aa * bb ** 3
    - 3 * X_u ** 2 * bb ** 4 - 4 * A_u * bb ** 4 * Z_u ** 4
)
after_collision, remainder = raw.quo_rem(Z_u ** 4)
assert not remainder

# This is a bivariate factorization of a degree-bounded chord radicand, not a
# surface elimination.  Retain precisely the odd squareclass.
UV = PolynomialRing(UQ, "v")
after_polynomial = UV(after_collision)
factorization = after_polynomial.factor()
factor_list = list(factorization)
odd_factors = [factor for factor, exponent in factor_list if int(exponent) % 2]
quartic = UV(factorization.unit())
for factor in odd_factors:
    quartic *= factor
square_factor = UV.one()
for factor, exponent in factor_list:
    square_factor *= factor ** (int(exponent) // 2)
assert after_polynomial == quartic * square_factor ** 2
assert quartic.degree() == 4
quartic_coefficients = list(quartic)
quartic_coefficients += [UQ.zero()] * (5 - len(quartic_coefficients))
e, d, c, b, a = quartic_coefficients

# Classical binary-quartic Jacobian.
I = UQ(12 * a * e - 3 * b * d + c ** 2)
J = UQ(72 * a * c * e + 9 * b * c * d - 27 * a * d ** 2 - 27 * b ** 2 * e - 2 * c ** 3)
A_child = UQ(-27 * I)
B_child = UQ(-27 * J)

# Exact global minimality at every finite base factor.
common = A_child.gcd(B_child)
removed_scalings = []
for factor, unused_exponent in common.factor():
    order_A = A_child.valuation(factor)
    order_B = B_child.valuation(factor)
    scaling_order = min(order_A // 4, order_B // 6)
    if scaling_order:
        A_child //= factor ** (4 * scaling_order)
        B_child //= factor ** (6 * scaling_order)
        removed_scalings.append((factor, int(scaling_order)))
assert A_child.degree() <= 8 and B_child.degree() <= 12
assert all(
    min(A_child.valuation(factor) // 4, B_child.valuation(factor) // 6) == 0
    for factor, unused in A_child.gcd(B_child).factor()
)
assert not (A_child.degree() <= 4 and B_child.degree() <= 6)

Delta_child = UQ(-16 * (4 * A_child ** 3 + 27 * B_child ** 2))
delta_factorization = list(Delta_child.factor())
finite_fibres = []
root_rank = 0
for factor, exponent in delta_factorization:
    order_delta = int(exponent)
    order_A = int(A_child.valuation(factor))
    order_B = int(B_child.valuation(factor))
    if order_delta == 1:
        kodaira = "I1"
        fibre_root_rank = 0
    elif order_A == 0 and order_B == 0:
        kodaira = f"I{order_delta}"
        fibre_root_rank = order_delta - 1
    elif order_A >= 1 and order_B == 1 and order_delta == 2:
        kodaira = "II"
        fibre_root_rank = 0
    else:
        raise ArithmeticError(("unclassified finite fibre", factor, order_A, order_B, order_delta))
    root_rank += fibre_root_rank * factor.degree()
    finite_fibres.append({
        "factor_coefficients_low_to_high": coefficients(factor),
        "factor_degree": int(factor.degree()),
        "orders_A_B_Delta": [order_A, order_B, order_delta],
        "kodaira": kodaira,
        "root_rank_contribution": int(fibre_root_rank * factor.degree()),
    })

infinity_orders = [
    int(8 - A_child.degree()),
    int(12 - B_child.degree()),
    int(24 - Delta_child.degree()),
]
if infinity_orders[2] == 0:
    infinity_kodaira = "smooth"
elif infinity_orders[2] == 1:
    infinity_kodaira = "I1"
elif infinity_orders[0] == 0 and infinity_orders[1] == 0:
    infinity_kodaira = f"I{infinity_orders[2]}"
    root_rank += infinity_orders[2] - 1
elif infinity_orders[0] >= 1 and infinity_orders[1] == 1 and infinity_orders[2] == 2:
    infinity_kodaira = "II"
else:
    raise ArithmeticError(("unclassified infinity fibre", infinity_orders))

euler_number = sum(item["factor_degree"] * item["orders_A_B_Delta"][2] for item in finite_fibres) + infinity_orders[2]
assert euler_number == 24
assert root_rank == 0
assert all(item["kodaira"] in ("I1", "II") for item in finite_fibres)
assert infinity_kodaira in ("smooth", "I1", "II")

quartic_values = []
for coefficient in quartic_coefficients:
    quartic_values.extend(UQ(coefficient).list())
jacobian_values = A_child.list() + B_child.list() + Delta_child.list()
payload = {
    "schema": "h92-q12o5867-smooth-rr-qq-v1",
    "reproducing_command": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/compile_h92_q12o5867_smooth_rr_qq.sage "
        "--output artifacts/local/elkies-k3/q12o5867-smooth-rr-qq.json"
    ),
    "inputs": {
        "target": {"path": str(TARGET.relative_to(ROOT)), "sha256": sha256(TARGET)},
        "surface": {"path": str(SURFACE.relative_to(ROOT)), "sha256": sha256(SURFACE)},
    },
    "divisor": {
        "class": "O+P-4F",
        "P_dot_O": 10,
        "fibre_twist": -4,
        "vertical_support": 0,
        "vertical_layers": 0,
    },
    "smooth_RR": {
        "ambient_dimension": 22,
        "collision_modulus_degree": 20,
        "condition_rank": int(condition_matrix.rank()),
        "h0": int(kernel.nrows()),
        "AA_degree_bound": aa_degree,
        "BB_degree_bound": bb_degree,
        "basis_pairs": [
            {
                "AA_coefficients_low_to_high": coefficients(AA),
                "BB_coefficients_low_to_high": coefficients(BB),
            }
            for AA, BB in pairs
        ],
        "no_resolved_vertical_rows": True,
    },
    "binary_quartic": {
        "coefficients_in_old_v_low_to_high": [coefficients(value) for value in quartic_coefficients],
        "factor_degrees_and_exponents_after_collision": [
            [int(factor.degree()), int(exponent)] for factor, exponent in factor_list
        ],
        "raw_after_collision_degree": int(after_polynomial.degree()),
        "square_factor_degree": int(square_factor.degree()),
        "maximum_rational_bits": rational_bits(quartic_values),
    },
    "child": {
        "minimal_A_coefficients_low_to_high": coefficients(A_child),
        "minimal_B_coefficients_low_to_high": coefficients(B_child),
        "discriminant_coefficients_low_to_high": coefficients(Delta_child),
        "degrees_A_B_Delta": [int(A_child.degree()), int(B_child.degree()), int(Delta_child.degree())],
        "removed_nonminimal_finite_scalings": [
            {"factor_coefficients_low_to_high": coefficients(factor), "order": order}
            for factor, order in removed_scalings
        ],
        "finite_fibres": finite_fibres,
        "infinity": {"orders_A_B_Delta": infinity_orders, "kodaira": infinity_kodaira},
        "euler_number": int(euler_number),
        "root_rank": int(root_rank),
        "ADE": "rootless",
        "MW_rank_if_rho19": 17,
        "maximum_A_B_Delta_rational_bits": rational_bits(jacobian_values),
    },
    "method": {
        "direct_smooth_chord_module_only": True,
        "resolved_vertical_rows": False,
        "groebner_or_surface_elimination": False,
        "stopped_before_17_section_endpoint_basis": True,
    },
    "runtime_seconds": time.monotonic() - started,
    "proof_boundary": "Exact QQ h0=2 pencil, binary quartic, minimal Jacobian, and fibre/rootless classification; the 17-section endpoint basis is not attempted.",
    "status": "PASS_EXACT_QQ_Q12O5867_SMOOTH_RR_ROOTLESS_JACOBIAN",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print(
    "Q12O5867SMOOTHQQ|ambient=22|rank={}|h0={}|quartic_degree={}|"
    "degrees_A_B_Delta={}|fibres={}|euler={}|root_rank={}|bits={}|runtime={:.3f}|"
    "status={}|output={}".format(
        condition_matrix.rank(), kernel.nrows(), quartic.degree(),
        payload["child"]["degrees_A_B_Delta"],
        [(item["factor_degree"], item["kodaira"]) for item in finite_fibres]
        + [([1], infinity_kodaira)] if infinity_kodaira != "smooth" else
        [(item["factor_degree"], item["kodaira"]) for item in finite_fibres],
        euler_number, root_rank, payload["child"]["maximum_A_B_Delta_rational_bits"],
        payload["runtime_seconds"], payload["status"], args.output.resolve(),
    )
)
