#!/usr/bin/env sage -python
"""Compile the physical q4/orbit208 2A5-to-3A3 pencil over QQ.

status: ACTIVE_PROOF
claim: exact QQ H0 plane, quartic, Jacobian, and 3A3 fibre profile
outputs: artifacts/local/elkies-k3/q24-2a5-physical-q4o208-rr-qq.json

The fibre has the literal divisor presentation

    D = O + P1229 + C10 + C8.

Write a generic numerator as ``a(T)+b*m`` with deg(a)<=2, b constant, and
``m=(y+Y_P)/(x-X_P)``, the chord through -P1229.  The common denominator is
the product of the two old I6 parameters.  On every nonidentity component at
an I6 root r, m has the constant residue

    c_r = Y_P(r)/(x_node(r)-X_P(r)).

Thus two exact rows cut the four-dimensional numerator space to H0=2.  The
old I6 factors occur squared in the chord branch numerator and are divided
exactly; no Groebner basis or polynomial factor search is used.
"""

import hashlib
import json
import time
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
SURFACE = LOCAL / "q24-a11-to-2a5-q8-resolved-rr-qq.json"
P1229_PATH = LOCAL / "q24-2a5-p1229-scaled-x-qq.json"
ROUTE = GENERATED / "elkies-k3-h3-a5a5-physical-q4o208-to-pinned-r17-certificate.json"
OUTPUT = LOCAL / "q24-2a5-physical-q4o208-rr-qq.json"

started = time.monotonic()


def log(stage, **fields):
    suffix = "|".join(f"{key}={value}" for key, value in fields.items())
    print(
        f"A5Q4O208RRQQ|stage={stage}|elapsed={time.monotonic()-started:.3f}"
        + (f"|{suffix}" if suffix else ""),
        flush=True,
    )


surface = json.loads(SURFACE.read_text())
p1229 = json.loads(P1229_PATH.read_text())
route = json.loads(ROUTE.read_text())
assert surface["status"] == "PASS_EXACT_Q24_A11_Q8_2A5_RESOLVED_RR"
assert p1229["status"] == "PASS_EXACT_QQ_P1229_POLYNOMIAL_SECTION"
assert route["status"] == "PASS_EXACT_PHYSICAL_Q4O208_3A3_TO_PINNED_R17"
assert route["candidate_id"] == {"q": 4, "old_fibre_degree": 2, "orbit_index": 208}
assert route["compiler_profile"]["literal_special_I4"]["identity"] == (
    "F_q4 = old_zero + P1229 + old_A11_component_10 + old_A11_component_8"
)

R = PolynomialRing(QQ, "T")
T = R.gen()
A = R([QQ(value) for value in surface["child"]["minimal_A_coefficients_low_to_high"]])
B = R([QQ(value) for value in surface["child"]["minimal_B_coefficients_low_to_high"]])
P = p1229["P1229"]
PX = R([QQ(value) for value in P["X_coefficients_low_to_high"]])
PY = R([QQ(value) for value in P["Y_coefficients_low_to_high"]])
assert (PX.degree(), PY.degree()) == (4, 6)
assert PY**2 == PX**3 + A * PX + B
log("LOAD", P1229_bits=P["height_profile"]["maximum_rational_bits"])


# The two exact old-I6 roots are already certified by the parent equation.
i6_factors = []
for record in surface["child"]["discriminant_factorization"]:
    if int(record["multiplicity"]) != 6:
        continue
    factor = R(sage_eval(record["factor"], locals={"T": T}))
    assert factor.degree() == 1
    i6_factors.append(factor.monic())
assert len(i6_factors) == 2
i6_roots = sorted((-factor[0] for factor in i6_factors))


def node_x(root):
    ar = QQ(A(root))
    br = QQ(B(root))
    assert ar
    value = -3 * br / (2 * ar)
    assert value**3 + ar * value + br == 0
    assert 3 * value**2 + ar == 0
    return value


residues = []
condition_rows = []
for root in i6_roots:
    node = node_x(root)
    assert PX(root) != node
    residue = QQ(PY(root) / (node - PX(root)))
    residues.append(residue)
    condition_rows.append((QQ.one(), root, root**2, residue))

conditions = matrix(QQ, condition_rows)
assert conditions.rank() == 2
tail = conditions.matrix_from_columns([2, 3])
assert tail.det()
kernel_rows = []
for free_column in (0, 1):
    solution = tail.solve_right(-conditions.column(free_column))
    row = vector(QQ, [0, 0, solution[0], solution[1]])
    row[free_column] = 1
    kernel_rows.append(row)
kernel = matrix(QQ, kernel_rows)
assert conditions * kernel.transpose() == matrix(QQ, 2, 2)
assert kernel.rank() == 2

# Pinned-prime comparison with the independently resolved modular plane.
p = ZZ(103)
Fp = GF(p)


def reduce_q(value):
    value = QQ(value)
    assert value.denominator() % p
    return Fp(value.numerator()) / Fp(value.denominator())


kernel_mod = matrix(Fp, [[reduce_q(value) for value in row] for row in kernel.rows()])
pinned_mod = matrix(Fp, [[1, 0, 92, 56], [0, 1, 89, 58]])
assert kernel_mod.row_space() == pinned_mod.row_space()
assert sorted(reduce_q(root) for root in i6_roots) == [Fp(68), Fp(89)]
assert sorted(reduce_q(value) for value in residues) == [Fp(44), Fp(53)]
log("RR", ambient=4, local_rank=2, kernel=2, status="PASS_EXACT")


# Compile the exact binary quartic over QQ(U).
UR = PolynomialRing(QQ, "U")
U0 = UR.gen()
KU = UR.fraction_field()
S = PolynomialRing(KU, "T")


def lift(poly):
    return S([KU(value) for value in R(poly).list()])


pairs = []
for row in kernel.rows():
    apoly = R(row[0] + row[1] * T + row[2] * T**2)
    pairs.append((lift(apoly), KU(row[3])))
(a0, b0), (a1, b1) = pairs
denominator = KU(U0) * b0 - b1
assert denominator
chord = (a1 - KU(U0) * a0) / denominator
radicand = (
    chord**4 - 6 * lift(PX) * chord**2 - 8 * lift(PY) * chord
    - 3 * lift(PX)**2 - 4 * lift(A)
)
numerator = S(radicand.numerator())
radicand_denominator = S(radicand.denominator())
old_square = prod((lift(factor) for factor in i6_factors), S.one())
quartic_numerator, remainder = numerator.quo_rem(old_square**2)
if remainder or quartic_numerator.degree() != 4:
    raise ArithmeticError(
        f"branch numerator is not the two old I6 squares times a quartic: "
        f"remainder={bool(remainder)}, quotient_degree={quartic_numerator.degree()}"
    )
quartic = S([
    KU(value) / KU(radicand_denominator[0])
    for value in quartic_numerator.list()
])
assert radicand == S(old_square)**2 * quartic
assert quartic.degree() == 4
log("QUARTIC", degree=4, removed="two_exact_old_I6_squares")

coefficients = list(quartic.list()) + [KU.zero()] * 5
e, d, c, b, a = coefficients[:5]
I = 12 * a * e - 3 * b * d + c**2
J = 72 * a * c * e + 9 * b * c * d - 27 * a * d**2 - 27 * b**2 * e - 2 * c**3
A_raw = -27 * I
B_raw = -27 * J

# The quartic has the scalar denominator ``denominator^4``.  Hence I and J
# have denominator powers 8 and 12, respectively; equivalently the standard
# Weierstrass gauge uses u=denominator^2.
A_child_value = KU(A_raw * denominator**8)
B_child_value = KU(B_raw * denominator**12)
if A_child_value.denominator() != 1 or B_child_value.denominator() != 1:
    raise ArithmeticError("q4 Jacobian retains a non-gauge base denominator")
A_child = UR(A_child_value)
B_child = UR(B_child_value)
Delta_child = UR(-16 * (4 * A_child**3 + 27 * B_child**2))
log(
    "CHILD_DEGREES",
    A=A_child.degree(), B=B_child.degree(), Delta=Delta_child.degree(),
)
assert (A_child.degree(), B_child.degree(), Delta_child.degree()) == (8, 12, 24)

# Three finite I4 fibres.  Use gcd/squarefree structure,
# avoiding a large full factorization of the degree-24 discriminant.
repeated_gcd = Delta_child.gcd(Delta_child.derivative()).monic()
repeated_support = repeated_gcd.squarefree_part().monic()
assert repeated_support.degree() == 3
assert repeated_gcd == repeated_support**3
nodal, remainder = Delta_child.quo_rem(repeated_support**4)
assert not remainder and nodal.degree() == 12 and nodal.is_squarefree()
assert A_child.gcd(repeated_support) == 1
assert B_child.gcd(repeated_support) == 1
infinity_orders = (
    8 - A_child.degree(),
    12 - B_child.degree(),
    24 - Delta_child.degree(),
)
assert infinity_orders == (0, 0, 0)
log("CHILD", finite="3I4+12I1", infinity="smooth", ADE="3A3", euler=24)


def rational_bits(values):
    values = [QQ(value) for value in values]
    return int(max(
        max(abs(value.numerator()).nbits(), value.denominator().nbits())
        for value in values
    ))


payload = {
    "schema": "elkies-k3.q24-2a5-physical-q4o208-rr-qq.v1",
    "status": "PASS_EXACT_QQ_PHYSICAL_Q4O208_3A3_RR_AND_JACOBIAN",
    "divisor": {
        "identity": "F_q4=O+P1229+C10+C8",
        "horizontal": "P1229",
        "P_dot_O": 0,
        "old_fibre_degree": 2,
        "q": 4,
    },
    "resolved_RR": {
        "ambient_basis": ["1", "T", "T^2", "m"],
        "chord": "m=(y+Y_P1229)/(x-X_P1229)",
        "old_I6_roots": [str(value) for value in i6_roots],
        "old_I6_chord_residues": [str(value) for value in residues],
        "condition_rows": [[str(value) for value in row] for row in conditions.rows()],
        "condition_rank": int(conditions.rank()),
        "kernel_basis": [[str(value) for value in row] for row in kernel.rows()],
        "kernel_dimension": int(kernel.nrows()),
        "mod103_kernel_row_space_match": True,
        "maximum_kernel_rational_bits": rational_bits(kernel.list()),
    },
    "quartic": {
        "coefficients_in_old_T_low_to_high": [str(value) for value in quartic.list()],
        "degree": 4,
        "exact_square_removal": [str(factor) + "^2" for factor in i6_factors],
    },
    "child": {
        "minimal_A_coefficients_low_to_high": [str(value) for value in A_child.list()],
        "minimal_B_coefficients_low_to_high": [str(value) for value in B_child.list()],
        "degrees_A_B_Delta": [8, 12, 24],
        "finite_repeated_support": str(repeated_support),
        "finite_nodal_factor": str(nodal),
        "finite_fibres": [
            {"kodaira": "I4", "count": 3},
            {"kodaira": "I1", "count": 12},
        ],
        "infinity": {"kodaira": "smooth", "orders_A_B_Delta": [0, 0, 0]},
        "ADE": "3A3",
        "root_rank": 9,
        "MW_rank_if_rho19": 8,
        "euler_number": 24,
    },
    "method": {
        "large_Groebner_required": False,
        "full_discriminant_factorization_required": False,
        "runtime_seconds": float(time.monotonic() - started),
    },
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (SURFACE, P1229_PATH, ROUTE)],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (SURFACE, P1229_PATH, ROUTE)
        },
    },
    "proof_boundary": (
        "Exact QQ H0 plane, binary quartic, globally minimal Jacobian, three finite I4 "
        "fibres, twelve finite I1 fibres, Euler number 24, and "
        "3A3/MW8 conditional on rho=19. The equation-effective C5 pointing and full "
        "old-curve-to-child marking are supplied lattice-wise by the route certificate "
        "and remain to be attached to the quartic/Jacobian equation."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A5Q4O208RRQQ|ambient=4|rank=2|kernel=2|quartic=4|fibres=3I4+12I1|"
    "kernel_bits={}|status={}|output={}".format(
        payload["resolved_RR"]["maximum_kernel_rational_bits"],
        payload["status"], OUTPUT,
    ),
    flush=True,
)
