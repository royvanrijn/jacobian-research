#!/usr/bin/env sage -python
"""Exact resolved-RR lift for A11/MW6 --q8 orbit12--> 2A5/MW7.

The exact translated section D=P12-O_pinned has degrees (16,24,6).  For the
divisor O+D-2F, write a generic-fibre function as ``a+b*m`` with
``m=(y-y(D))/(x-x(D))`` and

    a=AA/Z^2, deg AA<=10;  b=BB/Z, deg BB<=2.

The 14-dimensional ambient is cut out by the smooth-collision congruence
``AA*X=BB*Y mod Z^2``.  Rather than form a large rational nullspace, invert X
modulo Z^2 and impose only that the degree-11 coefficient of
``BB*Y/X mod Z^2`` vanish.  This is the bidirectional coefficient recurrence:
one linear condition on the three BB coefficients leaves the exact H0 plane.

The same plane is compiled to a binary quartic and globally minimal short
Jacobian.  No Groebner basis or nonlinear characteristic-zero solve is used.
"""

import hashlib
import json
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, inverse_mod, matrix


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q24-d12-to-a11-orbit42-resolved-rr-qq.json"
DIFFERENCE = LOCAL / "q24-a11-q8-difference-section-qq.json"
MODULAR_RR = LOCAL / "q24-a11-q8-resolved-rr-mod100003.json"
OUTPUT = LOCAL / "q24-a11-to-2a5-q8-resolved-rr-qq.json"

for path in (MODEL, DIFFERENCE, MODULAR_RR):
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

started = time.monotonic()


def log(stage, **fields):
    suffix = "|".join(f"{key}={value}" for key, value in fields.items())
    print(
        f"A11Q8RRQQ|stage={stage}|elapsed={time.monotonic()-started:.3f}"
        + (f"|{suffix}" if suffix else ""),
        flush=True,
    )


model = json.loads(MODEL.read_text())
difference = json.loads(DIFFERENCE.read_text())
modular_rr = json.loads(MODULAR_RR.read_text())
assert model["status"] == "PASS_EXACT_Q24_D12_Q6_A11_COMPONENT_VALUATION_RR"
assert difference["status"] == "PASS_EXACT_Q24_A11_Q8_DIFFERENCE_SECTION_QQ"
assert modular_rr["status"] == "PASS_MODP_A11_Q8_RESOLVED_RR_2A5"

VQ = PolynomialRing(QQ, "V")
V = VQ.gen()
KV = VQ.fraction_field()
A = VQ([QQ(value) for value in model["child"]["minimal_A_coefficients_low_to_high"]])
B = VQ([QQ(value) for value in model["child"]["minimal_B_coefficients_low_to_high"]])
section = difference["section"]
X = VQ([QQ(value) for value in section["X_coefficients_low_to_high"]])
Y = VQ([QQ(value) for value in section["Y_coefficients_low_to_high"]])
Z = VQ([QQ(value) for value in section["Z_coefficients_low_to_high"]])
if (X.degree(), Y.degree(), Z.degree()) != (16, 24, 6) or Z.leading_coefficient() != 1:
    raise ArithmeticError("exact translated section has wrong degree profile")
if Y**2 != X**3 + A * X * Z**4 + B * Z**6:
    raise ArithmeticError("exact translated section misses A11")
log("LOAD", section_degrees="16,24,6")

# -------------------------------------------------------------------------
# Exact 14 -> 2 resolved RR plane by one coefficient recurrence.
# -------------------------------------------------------------------------
collision_modulus = Z**2
log("INVERSE_START", modulus_degree=collision_modulus.degree())
X_inverse = VQ(inverse_mod(X, collision_modulus))
if (X * X_inverse) % collision_modulus != 1:
    raise ArithmeticError("X inverse modulo Z^2 failed")
log("INVERSE_DONE")

bb_basis = (VQ.one(), V, V**2)
base_residue = VQ((Y * X_inverse) % collision_modulus)
aa_residues = [VQ((BB * base_residue) % collision_modulus) for BB in bb_basis]
top_coefficients = [QQ(AA[11]) for AA in aa_residues]
nonzero_indices = [index for index, value in enumerate(top_coefficients) if value]
if not nonzero_indices:
    raise ArithmeticError("degree-11 recurrence did not leave a two-plane")
pivot = min(
    nonzero_indices,
    key=lambda index: max(
        abs(ZZ(top_coefficients[index].numerator())).nbits(),
        abs(ZZ(top_coefficients[index].denominator())).nbits(),
    ),
)
bb_kernel_rows = []
for free_index in range(3):
    if free_index == pivot:
        continue
    row = [QQ.zero()] * 3
    row[free_index] = QQ.one()
    row[pivot] = -top_coefficients[free_index] / top_coefficients[pivot]
    if sum(top_coefficients[index] * row[index] for index in range(3)):
        raise ArithmeticError("explicit one-row recurrence kernel failed")
    bb_kernel_rows.append(row)

rr_pairs = []
ambient_rows = []
for row in bb_kernel_rows:
    BB = VQ(sum(row[index] * bb_basis[index] for index in range(3)))
    AA = VQ((BB * Y * X_inverse) % collision_modulus)
    if AA.degree() > 10 or BB.degree() > 2:
        raise ArithmeticError("resolved RR degree bound failed")
    if (AA * X - BB * Y) % collision_modulus:
        raise ArithmeticError("resolved RR collision congruence failed")
    rr_pairs.append((AA, BB))
    ambient_rows.append(
        [AA[index] for index in range(11)] + [BB[index] for index in range(3)]
    )
log("RR", ambient=14, collision_rank=12, kernel=2, status="PASS_EXACT")

# Good-reduction plane comparison, allowing the unavoidable GL2 basis change.
p = ZZ(modular_rr["prime"])
Fp = GF(p)


def reduce_q(value):
    value = QQ(value)
    if value.denominator() % p == 0:
        raise ArithmeticError("bad denominator in RR plane reduction")
    return Fp(value.numerator()) / Fp(value.denominator())


exact_plane_mod = matrix(Fp, [[reduce_q(value) for value in row] for row in ambient_rows])
modular_plane = matrix(Fp, modular_rr["resolved_RR"]["basis_rows"])
if exact_plane_mod.row_space() != modular_plane.row_space():
    raise ArithmeticError("exact RR plane misses pinned modular plane")
plane_transport_mod = exact_plane_mod.solve_left(modular_plane)
if plane_transport_mod.det() == 0:
    raise ArithmeticError("modular RR plane transport is singular")

# -------------------------------------------------------------------------
# Binary quartic over QQ(U), primitive normalization, and exact child.
# -------------------------------------------------------------------------
TQ = PolynomialRing(QQ, "T")
T = TQ.gen()
KT = TQ.fraction_field()
WV = PolynomialRing(KT, "V")
W = WV.gen()


AA0, BB0 = rr_pairs[0]
AA1, BB1 = rr_pairs[1]
if AA0 * BB1 - AA1 * BB0 == 0:
    raise ArithmeticError("exact RR basis is chord-dependent")

# Form the chord square class fraction-free over QQ[T][V].  With
# N=AA1-T*AA0 and Db=T*BB0-BB1, one has m=N/(Z*Db).  Its fourth-power
# denominator is already a square, so only the displayed numerator matters.
WVT = PolynomialRing(TQ, "Vb")
Vb = WVT.gen()


def lift_nested(poly):
    return WVT([TQ(value) for value in VQ(poly).list()])


AA0b, BB0b, AA1b, BB1b, Zb, Xb, Yb, Ab = map(
    lift_nested, (AA0, BB0, AA1, BB1, Z, X, Y, A)
)
N_bi = AA1b - WVT(T) * AA0b
Db_bi = WVT(T) * BB0b - BB1b
raw_nested = (
    N_bi**4
    - 6 * Xb * N_bi**2 * Db_bi**2
    - 8 * Yb * N_bi * Db_bi**3
    - 3 * Xb**2 * Db_bi**4
    - 4 * Ab * Zb**4 * Db_bi**4
)
raw_branch = WV([KT(coefficient) for coefficient in raw_nested.list()])
rad_num_degree = int(raw_branch.degree())
rad_den_degree = int(4 * (Z.degree() + max(BB0.degree(), BB1.degree())))
log("SQUARECLASS_START", numerator_degree=rad_num_degree, denominator_degree=rad_den_degree)
quartic_nested, remainder = raw_nested.quo_rem(Zb**6)
if remainder or quartic_nested.degree() != 4:
    raise ArithmeticError(
        f"fraction-free numerator is not Z^6 times a quartic: "
        f"remainder={bool(remainder)}, quotient_degree={quartic_nested.degree()}"
    )
quartic = WV([KT(coefficient) for coefficient in quartic_nested.list()])
z_lift = WV([KT(value) for value in Z.list()])
decomposition = ((z_lift, 6), (quartic, 1))
if raw_branch != z_lift**6 * quartic:
    raise ArithmeticError("exact Z^6 quartic factorization failed")
log("SQUARECLASS_DONE", factor="Z^6*quartic", quartic_degree=4)

e, d, c, b, a = quartic.list()
I = 12 * a * e - 3 * b * d + c**2
J = 72 * a * c * e + 9 * b * c * d - 27 * a * d**2 - 27 * b**2 * e - 2 * c**3
if I.denominator() != 1 or J.denominator() != 1:
    raise ArithmeticError("q8 quartic invariants retain base denominators")
A_child = TQ(-27 * I)
B_child = TQ(-27 * J)
Delta_child = TQ(-16 * (4 * A_child**3 + 27 * B_child**2))
if (A_child.degree(), B_child.degree(), Delta_child.degree()) != (8, 12, 24):
    raise ArithmeticError("2A5 child is not a globally minimal elliptic K3 model")

log("FACTOR_CHILD_START")
delta_factors = sorted(
    ((factor.monic(), int(multiplicity)) for factor, multiplicity in Delta_child.factor()),
    key=lambda row: (row[1], row[0].degree(), str(row[0])),
)
factor_degree_profile = sorted((factor.degree(), multiplicity) for factor, multiplicity in delta_factors)
if factor_degree_profile != [(1, 6), (1, 6), (12, 1)]:
    raise ArithmeticError(f"child discriminant profile is {factor_degree_profile}, expected 2I6+12I1")
nodal_factor = next(factor for factor, multiplicity in delta_factors if multiplicity == 1)
if not nodal_factor.is_squarefree():
    raise ArithmeticError("degree-12 nodal support is not squarefree")
if A_child.gcd(Delta_child) != 1 or B_child.gcd(Delta_child) != 1:
    raise ArithmeticError("2A5 multiplicative fibres collide with c4/c6")
log("CHILD", fibres="2I6+12I1", ADE="2A5", euler=24, MW_if_rho19=7, status="PASS_EXACT")

# Pinned-prime child regression.  A GL2 change of the RR plane may change the
# displayed new-base coordinate, so compare invariant fibre data, not raw rows.
TF = PolynomialRing(Fp, "T")
A_mod = TF([reduce_q(value) for value in A_child.list()])
B_mod = TF([reduce_q(value) for value in B_child.list()])
D_mod = TF(-16 * (4 * A_mod**3 + 27 * B_mod**2))
if sorted((factor.degree(), int(multiplicity)) for factor, multiplicity in D_mod.factor()) != [
    (1, 1), (1, 6), (1, 6), (3, 1), (4, 1), (4, 1)
]:
    raise ArithmeticError("exact 2A5 child misses pinned modular fibre pattern")


def text_hash(rows):
    blob = json.dumps(rows, separators=(",", ":"), ensure_ascii=True).encode()
    return hashlib.sha256(blob).hexdigest()


input_hashes = {
    str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
    for path in (MODEL, DIFFERENCE, MODULAR_RR)
}
payload = {
    "schema": "elkies-k3.h3-q24-a11-to-2a5-q8-resolved-rr-qq.v1",
    "status": "PASS_EXACT_Q24_A11_Q8_2A5_RESOLVED_RR",
    "inputs": {"paths": list(input_hashes), "sha256": input_hashes},
    "selected_section": {
        "relation": "P12-O_pinned after translating O_pinned to the Weierstrass origin",
        "P_dot_O": 6,
        "I12_component_depth": 6,
        "minimal_equation_degrees": [16, 24, 6],
        "exact_weierstrass_identity": True,
    },
    "resolved_RR": {
        "divisor_formula": "O+(P12-O_pinned)-2F in the translated group frame",
        "ambient_dimension": 14,
        "AA_degree_bound": 10,
        "BB_degree_bound": 2,
        "collision_modulus": "Z^2",
        "collision_rank": 12,
        "recurrence": "AA=BB*Y/X mod Z^2; coefficient degree 11 vanishes",
        "recurrence_rank": 1,
        "kernel_dimension": 2,
        "h0": 2,
        "basis_pairs": [
            {
                "AA_coefficients_low_to_high": [str(value) for value in AA.list()],
                "BB_coefficients_low_to_high": [str(value) for value in BB.list()],
            }
            for AA, BB in rr_pairs
        ],
        "kernel_basis_sha256": text_hash([
            [[str(value) for value in AA.list()], [str(value) for value in BB.list()]]
            for AA, BB in rr_pairs
        ]),
        "mod_100003_plane_matches": True,
        "mod_100003_GL2_transport": [list(map(int, row)) for row in plane_transport_mod.rows()],
    },
    "quartic": {
        "degree": 4,
        "raw_radicand_degrees": [rad_num_degree, rad_den_degree],
        "squarefree_decomposition": [
            [int(factor.degree()), int(multiplicity)] for factor, multiplicity in decomposition
        ],
        "coefficients_in_T_low_to_high": [str(value) for value in quartic.list()],
        "I": str(I),
        "J": str(J),
    },
    "child": {
        "minimal_A_coefficients_low_to_high": [str(value) for value in A_child.list()],
        "minimal_B_coefficients_low_to_high": [str(value) for value in B_child.list()],
        "discriminant_factorization": [
            {"factor": str(factor), "multiplicity": multiplicity}
            for factor, multiplicity in delta_factors
        ],
        "finite_fibres_geometric": [
            {"kodaira": "I6", "count": 2},
            {"kodaira": "I1", "count": 12},
        ],
        "infinity_kind": "smooth",
        "root_lattice": "2A5",
        "root_rank": 10,
        "root_determinant": 36,
        "euler_number": 24,
        "MW_rank_if_rho19": 7,
    },
    "verification": {
        "exact_section_identity": True,
        "exact_collision_congruence": True,
        "exact_H0_dimension": True,
        "exact_quartic_squareclass": True,
        "exact_minimal_jacobian": True,
        "exact_2A5_fibre_classification": True,
        "mod_100003_plane_regression": True,
        "mod_100003_child_regression": True,
    },
    "large_Groebner_required": False,
    "next_required": "transport an explicit old-A11 component as the chosen 2A5 zero and close the full marked NS handoff",
    "proof_boundary": (
        "The translated section, complete two-dimensional H0 plane, quartic, globally minimal Jacobian, "
        "2A5 fibre classification, Euler number, and MW rank conditional on rho=19 are exact over QQ. "
        "The orbit12 orientation is bound by the pinned mod-100003 plane regression. The explicit child zero, "
        "component marking, and bidirectional NS-to-equation transport remain to be attached before promotion."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
log("RESULT", status=payload["status"])
print(f"OUTPUT|{OUTPUT}", flush=True)
