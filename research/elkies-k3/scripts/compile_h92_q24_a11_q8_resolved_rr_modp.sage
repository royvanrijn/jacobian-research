#!/usr/bin/env sage -python
"""Compile the marked A11 q8 resolved RR plane modulo the pinned prime.

Translate the elliptic group origin to the exact pinned A11 zero.  The second
horizontal point is then the already marked modular section
``D=P12-O_pinned``, with pole order six and denominator Z.

For the divisor ``O+D-2F`` the complete chord-frame ambient is

    a = AA/Z^2,  deg(AA) <= 10,
    b = BB/Z,    deg(BB) <= 2,

in ``a+b*m``, where ``m=(y-y(D))/(x-x(D))``.  Thus the ambient has dimension
14.  Regularity at the six smooth collision points is the single polynomial
congruence ``AA*X=BB*Y mod Z^2`` (twelve scalar rows).  The expected kernel is
the full two-dimensional H0 plane.

The script compiles that same plane to its squarefree binary quartic and
classifies the modular Jacobian.  It is a pinned-good-reduction construction
regression, not a characteristic-zero equation certificate.  No Groebner
basis or nonlinear section solve is used.
"""

import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q24-d12-to-a11-orbit42-resolved-rr-qq.json"
MARKED = LOCAL / "q24-a11-q8-horizontal-points-mod100003.json"
OUTPUT = LOCAL / "q24-a11-q8-resolved-rr-mod100003.json"

for path in (MODEL, MARKED):
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

model = json.loads(MODEL.read_text())
marked = json.loads(MARKED.read_text())
assert model["status"] == "PASS_EXACT_Q24_D12_Q6_A11_COMPONENT_VALUATION_RR"
assert marked["status"] == "PASS_Q24_A11_Q8_HORIZONTAL_POINTS_RECONSTRUCTION_MODP"

p = ZZ(marked["prime"])
assert p == 100003
F = GF(p)
FV = PolynomialRing(F, "V")
V = FV.gen()


def red_q(value):
    value = QQ(value)
    if value.denominator() % p == 0:
        raise ArithmeticError("bad exact denominator at pinned prime")
    return F(value.numerator()) / F(value.denominator())


A = FV([red_q(value) for value in model["child"]["minimal_A_coefficients_low_to_high"]])
B = FV([red_q(value) for value in model["child"]["minimal_B_coefficients_low_to_high"]])
difference = marked["marked_difference_P12_minus_Opinned"]["section"]
X = FV(difference["X_coefficients_low_to_high"])
Y = FV(difference["Y_coefficients_low_to_high"])
Z = FV(difference["Z_coefficients_low_to_high"])
if (X.degree(), Y.degree(), Z.degree()) != (16, 24, 6) or Z.leading_coefficient() != 1:
    raise ArithmeticError("marked difference has the wrong resolved degree profile")
if Y**2 != X**3 + A * X * Z**4 + B * Z**6:
    raise ArithmeticError("marked difference misses the reduced A11 equation")

# -------------------------------------------------------------------------
# Complete 14 -> 2 collision RR kernel.
# -------------------------------------------------------------------------
aa_degree = 10
bb_degree = 2
ambient = [
    (V**degree, FV.zero()) for degree in range(aa_degree + 1)
] + [
    (FV.zero(), V**degree) for degree in range(bb_degree + 1)
]
collision_modulus = Z**2
remainders = [FV((AA * X - BB * Y) % collision_modulus) for AA, BB in ambient]
rows = [
    [remainder[degree] for remainder in remainders]
    for degree in range(collision_modulus.degree())
]
condition_matrix = matrix(F, rows)
rank = condition_matrix.rank()
kernel = condition_matrix.right_kernel().basis_matrix()
if condition_matrix.nrows() != 12 or condition_matrix.ncols() != 14:
    raise ArithmeticError("resolved RR matrix dimensions changed")
if rank != 12 or kernel.nrows() != 2:
    raise ArithmeticError(f"resolved RR rank/kernel is {rank}/{kernel.nrows()}, expected 12/2")


def pair_from_row(row):
    AA = sum(row[index] * V**index for index in range(aa_degree + 1))
    offset = aa_degree + 1
    BB = sum(row[offset + index] * V**index for index in range(bb_degree + 1))
    if (AA * X - BB * Y) % collision_modulus:
        raise ArithmeticError("displayed RR basis misses collision congruence")
    return FV(AA), FV(BB)


rr_pairs = [pair_from_row(row) for row in kernel.rows()]

# -------------------------------------------------------------------------
# Compile the exact modular plane to a binary quartic over F(U).
# -------------------------------------------------------------------------
FU = PolynomialRing(F, "U")
U0 = FU.gen()
KU = FU.fraction_field()
RV = PolynomialRing(KU, "V")
W = RV.gen()


def lift_v(poly):
    return RV([KU(value) for value in FV(poly).list()])


AA0, BB0 = rr_pairs[0]
AA1, BB1 = rr_pairs[1]
z_lift, x_lift, y_lift, a_lift = map(lift_v, (Z, X, Y, A))
a0 = lift_v(AA0) / z_lift**2
b0 = lift_v(BB0) / z_lift
a1 = lift_v(AA1) / z_lift**2
b1 = lift_v(BB1) / z_lift
determinant = a0 * b1 - a1 * b0
if not determinant:
    raise ArithmeticError("RR basis is dependent in the chord direction")

chord = (a1 - KU(U0) * a0) / (KU(U0) * b0 - b1)
x_point = x_lift / z_lift**2
y_point = y_lift / z_lift**3
radicand = (
    chord**4
    - 6 * x_point * chord**2
    - 8 * y_point * chord
    - 3 * x_point**2
    - 4 * a_lift
)
rad_num = RV(radicand.numerator())
rad_den = RV(radicand.denominator())
raw_branch = rad_num * rad_den
factorization = raw_branch.factor()
odd_branch = RV(factorization.unit())
branch_factor_profile = []
for factor, exponent in factorization:
    exponent = int(exponent)
    branch_factor_profile.append((str(factor), int(factor.degree()), exponent))
    if exponent % 2:
        odd_branch *= factor
if odd_branch.degree() != 4:
    raise ArithmeticError(
        f"resolved RR square class has degree {odd_branch.degree()}, expected 4; "
        f"profile={branch_factor_profile}"
    )

quartic = odd_branch
coefficients = list(quartic.list()) + [KU.zero()] * 5
e, d, c, b, a = coefficients[:5]
Iinv = 12 * a * e - 3 * b * d + c**2
Jinv = 72 * a * c * e + 9 * b * c * d - 27 * a * d**2 - 27 * b**2 * e - 2 * c**3
jacA = -27 * Iinv
jacB = -27 * Jinv
jacDelta = -16 * (4 * jacA**3 + 27 * jacB**2)


def kodaira_data(ord_a, ord_b, ord_delta):
    if ord_a == 0 or ord_b == 0:
        n = int(ord_delta)
        return n - 1, n, n, f"I{n}"
    if ord_delta == 2:
        return 0, 1, 2, "II"
    if ord_delta == 3:
        return 1, 2, 3, "III"
    if ord_delta == 4:
        return 2, 3, 4, "IV"
    if ord_delta == 6 and ord_a >= 2 and ord_b >= 3:
        return 4, 4, 6, "I0*"
    if ord_delta >= 7 and ord_a == 2 and ord_b == 3:
        n = int(ord_delta - 6)
        return n + 4, 4, n + 6, f"I{n}*"
    if ord_delta == 8:
        return 6, 3, 8, "IV*"
    if ord_delta == 9:
        return 7, 2, 9, "III*"
    if ord_delta == 10:
        return 8, 1, 10, "II*"
    raise ArithmeticError((ord_a, ord_b, ord_delta))


def valuation(value, factor):
    return int(value.numerator().valuation(factor) - value.denominator().valuation(factor))


factors = set()
for value in (jacA, jacB, jacDelta):
    for polynomial in (value.numerator(), value.denominator()):
        factors.update(factor for factor, unused in polynomial.factor())

finite_fibres = []
root_rank = 0
root_determinant = 1
euler = 0
for factor in sorted(factors, key=str):
    raw = (valuation(jacA, factor), valuation(jacB, factor), valuation(jacDelta, factor))
    scale = min(raw[0] // 4, raw[1] // 6, raw[2] // 12)
    orders = (raw[0] - 4 * scale, raw[1] - 6 * scale, raw[2] - 12 * scale)
    if orders[2] <= 0:
        continue
    fibre_rank, determinant_local, euler_local, kind = kodaira_data(*orders)
    degree = int(factor.degree())
    root_rank += degree * fibre_rank
    root_determinant *= determinant_local**degree
    euler += degree * euler_local
    finite_fibres.append(
        {"factor": str(factor), "degree": degree, "orders_A_B_Delta": list(orders), "kind": kind}
    )

raw_infinity = tuple(
    int(value.denominator().degree() - value.numerator().degree())
    for value in (jacA, jacB, jacDelta)
)
scale_infinity = min(raw_infinity[0] // 4, raw_infinity[1] // 6, raw_infinity[2] // 12)
infinity_orders = tuple(
    raw_infinity[index] - (4, 6, 12)[index] * scale_infinity for index in range(3)
)
infinity_kind = "smooth"
if infinity_orders[2] > 0:
    fibre_rank, determinant_local, euler_local, infinity_kind = kodaira_data(*infinity_orders)
    root_rank += fibre_rank
    root_determinant *= determinant_local
    euler += euler_local

status = (
    "PASS_MODP_A11_Q8_RESOLVED_RR_2A5"
    if (root_rank, root_determinant, euler) == (10, 36, 24)
    else "FAIL_MODP_A11_Q8_RESOLVED_RR_CHILD_SIGNATURE"
)
payload = {
    "schema": "elkies-k3.h3-q24-a11-q8-resolved-rr-modp.v1",
    "status": status,
    "prime": int(p),
    "inputs": {
        str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in (MODEL, MARKED)
    },
    "zero_translation": "equation group translated by -O_pinned",
    "horizontal_section": "P12-O_pinned",
    "resolved_RR": {
        "ambient_dimension": 14,
        "AA_degree_bound": aa_degree,
        "BB_degree_bound": bb_degree,
        "collision_modulus_degree": int(collision_modulus.degree()),
        "condition_rank": int(rank),
        "kernel_dimension": int(kernel.nrows()),
        "basis_rows": [list(map(int, row)) for row in kernel.rows()],
        "pairs": [
            {
                "AA_coefficients_low_to_high": list(map(int, AA.list())),
                "BB_coefficients_low_to_high": list(map(int, BB.list())),
            }
            for AA, BB in rr_pairs
        ],
    },
    "quartic": {
        "degree": int(quartic.degree()),
        "coefficients_low_to_high": [str(value) for value in quartic.list()],
        "branch_factor_profile": branch_factor_profile,
    },
    "jacobian": {
        "A": str(jacA),
        "B": str(jacB),
        "finite_fibres": finite_fibres,
        "infinity_orders_A_B_Delta": list(infinity_orders),
        "infinity_kind": infinity_kind,
        "root_rank": int(root_rank),
        "root_determinant": int(root_determinant),
        "euler_number": int(euler),
        "MW_rank_if_rho19": int(17 - root_rank),
    },
    "large_Groebner_required": False,
    "proof_boundary": (
        "Exact pinned-prime arithmetic constructs the complete 14-to-2 collision RR kernel, "
        "the binary quartic, and the modular Jacobian fibre signature. The horizontal section, "
        "plane, quartic, minimized Jacobian, and marking still require characteristic-zero replay."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    f"A11Q8RRMOD|prime={p}|ambient=14|rank={rank}|kernel={kernel.nrows()}|"
    f"quartic_degree={quartic.degree()}|root_rank={root_rank}|root_det={root_determinant}|"
    f"euler={euler}|MW={17-root_rank}|status={status}",
    flush=True,
)
print(f"OUTPUT|{OUTPUT}", flush=True)
