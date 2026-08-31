#!/usr/bin/env sage
"""Compile the fixed-corridor terminal A1 fibration from the q12 endpoint.

For the exact rootless section P constructed by
``construct_h92_fixed_final_a1_horizontal_from_q12_endpoint_qq.sage``, use

    D = O + P - 2F,   P.O=6.

The smooth chord module has deg(AA)<=10, deg(BB)<=2.  This is a 14-to-2
linear calculation with no vertical component rows and no Groebner basis.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import PolynomialRing, QQ, ZZ, matrix


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--target", type=Path, default=LOCAL / "fixed-final-a1-horizontal-from-q12-endpoint-qq.json")
parser.add_argument("--surface", type=Path, default=LOCAL / "q12o5867-smooth-rr-qq.json")
parser.add_argument("--output", type=Path, default=LOCAL / "fixed-final-a1-reverse-rr-qq.json")
parser.add_argument("--target-status", default="PASS_EXACT_QQ_FIXED_FINAL_A1_REVERSE_HORIZONTAL")
parser.add_argument("--surface-status", default="PASS_EXACT_QQ_Q12O5867_SMOOTH_RR_ROOTLESS_JACOBIAN")
parser.add_argument("--edge", default="rootless/MW17 reverse to A1/MW16")
parser.add_argument("--expected-i2-count", type=int, default=1)
parser.add_argument("--expected-ade", default="A1")
parser.add_argument("--expected-mw-rank", type=int, default=16)
parser.add_argument("--result-status", default="PASS_EXACT_QQ_FIXED_FINAL_A1_REVERSE_RR_JACOBIAN")
args = parser.parse_args()
TARGET = args.target.resolve()
SURFACE = args.surface.resolve()
OUTPUT = args.output.resolve()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficients(poly):
    return [str(value) for value in poly.list()]


def rational_bits(values):
    answer = 0
    for value in values:
        value = QQ(value)
        answer = max(answer, abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())
    return int(answer)


started = time.monotonic()
target_artifact = json.loads(TARGET.read_text())
surface_artifact = json.loads(SURFACE.read_text())
assert target_artifact["status"] == args.target_status
assert surface_artifact["status"] == args.surface_status
assert target_artifact["fixed_edge"]["P_dot_O"] == 6

R = PolynomialRing(QQ, "u")
u = R.gen()
K = R.fraction_field()
A = R(surface_artifact["child"]["minimal_A_coefficients_low_to_high"])
B = R(surface_artifact["child"]["minimal_B_coefficients_low_to_high"])
section = target_artifact["section"]
X = R(section["x_numerator_coefficients_low_to_high"])
Y = R(section["y_numerator_coefficients_low_to_high"])
Z = R(section["Z_coefficients_low_to_high"])
assert (X.degree(), Y.degree(), Z.degree()) == (16, 24, 6)
assert Y ** 2 == X ** 3 + A * X * Z ** 4 + B * Z ** 6

# For D=O+P+kF, the compact chord bounds are
# deg(AA)<=2(P.O)+k and deg(BB)<=P.O-2+k.  Here k=-2.
aa_degree = 10
bb_degree = 2
ambient = [(u ** degree, R.zero()) for degree in range(aa_degree + 1)] + [
    (R.zero(), u ** degree) for degree in range(bb_degree + 1)
]
collision_modulus = Z ** 2
assert len(ambient) == 14 and collision_modulus.degree() == 12
remainders = [R((AA * X - BB * Y) % collision_modulus) for AA, BB in ambient]
condition_matrix = matrix(QQ, [
    [remainder[degree] for remainder in remainders]
    for degree in range(collision_modulus.degree())
])
assert condition_matrix.dimensions() == (12, 14)
assert condition_matrix.rank() == 12
# Sage 10.9 routes a dense-QQ kernel through IML even when ``flint`` is
# requested.  Clear the one global denominator first so the documented FLINT
# integer-kernel implementation is actually selected.  ``computed`` avoids an
# unnecessary echelon normalization.
integer_conditions, unused_denominator = condition_matrix._clear_denom()
kernel = integer_conditions.right_kernel_matrix(algorithm="flint", basis="computed").change_ring(QQ)
assert kernel.nrows() == kernel.rank() == 2

pairs = []
for row in kernel.rows():
    AA = sum((row[index] * ambient[index][0] for index in range(14)), R.zero())
    BB = sum((row[index] * ambient[index][1] for index in range(14)), R.zero())
    assert (AA * X - BB * Y) % collision_modulus == 0
    pairs.append((AA, BB))
AA0, BB0 = pairs[0]
AA1, BB1 = pairs[1]

# Compile the chord radicand with s as the new base parameter.
S = PolynomialRing(QQ, "s")
s = S.gen()
KS = S.fraction_field()
US = PolynomialRing(KS, "u")


def lift(poly):
    return US([KS(value) for value in R(poly).list()])


aa = lift(AA0) + KS(s) * lift(AA1)
bb = lift(BB0) + KS(s) * lift(BB1)
Xu, Yu, Zu, Au = map(lift, (X, Y, Z, A))
raw = aa ** 4 - 6 * Xu * aa ** 2 * bb ** 2 + 8 * Yu * aa * bb ** 3 - 3 * Xu ** 2 * bb ** 4 - 4 * Au * bb ** 4 * Zu ** 4
after_collision, remainder = raw.quo_rem(Zu ** 4)
assert not remainder

# Strip the forced square by multivariate gcd rather than factoring the large
# bivariate radicand.  This is exact in QQ[s,u].
SU = PolynomialRing(QQ, names=("s", "u"))
ss, uu = SU.gens()
after = SU(after_collision)
square_factor = after.gcd(after.derivative(ss)).gcd(after.derivative(uu))
quartic, remainder = after.quo_rem(square_factor ** 2)
assert not remainder
assert quartic.degree(uu) == 4
assert quartic.gcd(quartic.derivative(ss)).gcd(quartic.derivative(uu)).is_constant()

U_over_S = PolynomialRing(S, "u")
quartic_univariate = U_over_S(quartic)
square_factor_univariate = U_over_S(square_factor)
quartic_coefficients = [S(quartic_univariate[degree]) for degree in range(5)]
e, d, c, b, a = quartic_coefficients
I = S(12 * a * e - 3 * b * d + c ** 2)
J = S(72 * a * c * e + 9 * b * c * d - 27 * a * d ** 2 - 27 * b ** 2 * e - 2 * c ** 3)
A_child = S(-27 * I)
B_child = S(-27 * J)

# Remove any finite fourth/sixth common scaling.  Factoring only gcd(A,B) is
# tiny compared with factoring the chord radicand or the surface ideal.
removed_scalings = []
for factor, unused in A_child.gcd(B_child).factor():
    order = min(A_child.valuation(factor) // 4, B_child.valuation(factor) // 6)
    if order:
        A_child //= factor ** (4 * order)
        B_child //= factor ** (6 * order)
        removed_scalings.append((factor, int(order)))
assert A_child.degree() <= 8 and B_child.degree() <= 12
print(
    "FIXEDFINALA1RR_STAGE|"
    f"square_factor_total_degree={square_factor.total_degree()}|"
    f"quartic_degrees_s_u=({quartic.degree(ss)},{quartic.degree(uu)})|"
    f"preminimal_degrees=({A_child.degree()},{B_child.degree()})|"
    f"removed_scalings={[(factor.degree(), order) for factor, order in removed_scalings]}",
    flush=True,
)
assert not (A_child.degree() <= 4 and B_child.degree() <= 6)

Delta = S(-16 * (4 * A_child ** 3 + 27 * B_child ** 2))
assert Delta.degree() <= 24

# Semistable classification without factoring the degree-24 discriminant.
# The repeated part gives the expected I2 supports; the residual factor is
# squarefree and disjoint from A and B, hence gives geometric I1 fibres.
repeated = Delta.gcd(Delta.derivative()).monic()
assert repeated.degree() == args.expected_i2_count
residual, remainder = Delta.quo_rem(repeated ** 2)
assert not remainder and residual.degree() == 24 - 2 * args.expected_i2_count
assert repeated.is_squarefree()
assert residual.gcd(residual.derivative()).is_constant()
assert repeated.gcd(residual).is_constant()
assert repeated.gcd(A_child).is_constant() and repeated.gcd(B_child).is_constant()
assert residual.gcd(A_child).is_constant() and residual.gcd(B_child).is_constant()
infinity_orders = [int(8 - A_child.degree()), int(12 - B_child.degree()), int(24 - Delta.degree())]
assert infinity_orders == [0, 0, 0]

quartic_values = []
for coefficient in quartic_coefficients:
    quartic_values.extend(coefficient.list())
jacobian_values = A_child.list() + B_child.list() + Delta.list()
payload = {
    "schema": "elkies-k3.fixed-reverse-smooth-rr-qq.v1",
    "status": args.result_status,
    "reproducing_command": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/compile_h92_fixed_final_a1_reverse_rr_qq.sage"
    ),
    "inputs": {
        "target": {"path": str(TARGET.relative_to(ROOT)), "sha256": sha256(TARGET)},
        "surface": {"path": str(SURFACE.relative_to(ROOT)), "sha256": sha256(SURFACE)},
    },
    "divisor": {"class": "O+P-2F", "P_dot_O": 6, "fibre_twist": -2, "vertical_support": 0},
    "smooth_RR": {
        "ambient_dimension": 14,
        "collision_modulus_degree": 12,
        "condition_rank": int(condition_matrix.rank()),
        "h0": int(kernel.nrows()),
        "AA_degree_bound": aa_degree,
        "BB_degree_bound": bb_degree,
        "basis_pairs": [
            {"AA_coefficients_low_to_high": coefficients(AA), "BB_coefficients_low_to_high": coefficients(BB)}
            for AA, BB in pairs
        ],
        "no_resolved_vertical_rows": True,
    },
    "binary_quartic": {
        "coefficients_in_old_u_low_to_high": [coefficients(value) for value in quartic_coefficients],
        "square_factor_coefficients_in_old_u_low_to_high": [
            coefficients(S(square_factor_univariate[degree]))
            for degree in range(square_factor_univariate.degree() + 1)
        ],
        "square_factor_total_degree": int(square_factor.total_degree()),
        "maximum_rational_bits": rational_bits(quartic_values),
    },
    "child": {
        "minimal_A_coefficients_low_to_high": coefficients(A_child),
        "minimal_B_coefficients_low_to_high": coefficients(B_child),
        "discriminant_coefficients_low_to_high": coefficients(Delta),
        "degrees_A_B_Delta": [int(A_child.degree()), int(B_child.degree()), int(Delta.degree())],
        "finite_fibres": [
            {"factor_coefficients_low_to_high": coefficients(repeated), "factor_degree": int(repeated.degree()), "orders_A_B_Delta": [0, 0, 2], "kodaira": "I2", "root_rank_contribution": int(repeated.degree())},
            {"factor_coefficients_low_to_high": coefficients(residual), "factor_degree": int(residual.degree()), "orders_A_B_Delta": [0, 0, 1], "kodaira": "I1", "root_rank_contribution": 0},
        ],
        "infinity": {"orders_A_B_Delta": infinity_orders, "kodaira": "smooth"},
        "euler_number": 24,
        "root_rank": args.expected_i2_count,
        "ADE": args.expected_ade,
        "MW_rank_if_rho19": args.expected_mw_rank,
        "removed_nonminimal_finite_scalings": [
            {"factor_coefficients_low_to_high": coefficients(factor), "order": order}
            for factor, order in removed_scalings
        ],
        "maximum_A_B_Delta_rational_bits": rational_bits(jacobian_values),
    },
    "method": {
        "direct_smooth_chord_module_only": True,
        "multivariate_gcd_square_stripping": True,
        "full_discriminant_factorization": False,
        "groebner_or_surface_elimination": False,
        "runtime_seconds": time.monotonic() - started,
    },
    "proof_boundary": "Exact QQ h0=2, binary quartic, minimal semistable Jacobian and fibre classification. Pointing and component marking are separate gates.",
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "FIXEDREVERSERR|edge={}|ambient=14|rank={}|h0={}|quartic=4|degrees={}|I2={}|I1={}|"
    "root_rank={}|bits={}|seconds={:.3f}|status={}|output={}".format(
        args.edge,
        condition_matrix.rank(), kernel.nrows(), payload["child"]["degrees_A_B_Delta"],
        args.expected_i2_count, residual.degree(), args.expected_i2_count,
        payload["child"]["maximum_A_B_Delta_rational_bits"], payload["method"]["runtime_seconds"],
        payload["status"], OUTPUT,
    ), flush=True,
)
