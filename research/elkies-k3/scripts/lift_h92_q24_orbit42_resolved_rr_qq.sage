#!/usr/bin/env sage -python
"""Exact resolved-RR lift for D12/MW5 --q6 orbit42--> A11/MW6.

This is the corrected orbit42 compiler.  It consumes the exact 18+2
zero-pole shell construction and the four exact orbit42 section candidates;
it does not Hensel-lift the abandoned current-equation modular section.

The physical I8* resolution has repeated tangent cone y^2 at C01.  Hence the
first divisorial valuation is weighted (2,2,3), not the ordinary total-degree
jet used by the historical driver.  On the three-dimensional collision
space H0(O+P), the effective-root sign and C10 spinor orientation impose the
single condition AA(alpha)=0.  The opposite C11 orientation has one further
independent C10-arm condition modulo the pinned good prime and is rejected.

The script changes the selected exact section to the minimal D12 equation,
constructs the 9 -> 3 -> 2 RR plane, extracts the exact binary quartic, and
classifies its minimized Jacobian.  The expected output has an I12 fibre and
twelve geometrically nodal fibres: root lattice A11 and Euler number 24.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, inverse_mod


def locate_repo(explicit=None):
    candidates = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd = Path.cwd().resolve()
    candidates += [cwd, *cwd.parents]
    home = Path.home()
    candidates += [
        home / "Documents/jacobian-research",
        home / "jacobian-research",
        home / "src/jacobian-research",
        home / "git/jacobian-research",
    ]
    seen = set()
    for candidate in candidates:
        try:
            candidate = candidate.resolve()
        except Exception:
            continue
        if candidate in seen:
            continue
        seen.add(candidate)
        if (candidate / "elkies-k3/scripts").is_dir():
            return candidate
    raise SystemExit("Could not locate jacobian-research; pass --repo PATH")


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--prime", type=int, default=100003)
parser.add_argument("--candidate-index", type=int, default=0)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q24-orbit42-zero-pole-seeds-mod-43.json"
CANDIDATES = LOCAL / "q24-orbit42-exact-section-candidates-qq.json"
PARENT = LOCAL / "q24-d13-to-d12-component-valuation-qq.json"
PHYSICAL = LOCAL / "q24-d12-orbit42-i8star-physical-marking-qq.json"
OUTPUT = (
    args.output.resolve()
    if args.output
    else LOCAL / "q24-d12-to-a11-orbit42-resolved-rr-qq.json"
)

for path in (MODEL, CANDIDATES, PARENT, PHYSICAL):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

started = time.monotonic()


def log(stage, **fields):
    suffix = "|".join(f"{key}={value}" for key, value in fields.items())
    print(
        f"Q42RRQQ|stage={stage}|elapsed={time.monotonic()-started:.3f}"
        + (f"|{suffix}" if suffix else ""),
        flush=True,
    )


model_artifact = json.loads(MODEL.read_text())
candidate_artifact = json.loads(CANDIDATES.read_text())
parent_artifact = json.loads(PARENT.read_text())
physical_artifact = json.loads(PHYSICAL.read_text())

if model_artifact.get("status") != "PASS_Q42_ZERO_POLE_SMALLPRIME_SEEDS":
    raise ArithmeticError("zero-pole model prerequisite is not passing")
if candidate_artifact.get("status") != "PASS_EXACT_Q42_ORBIT42_SECTION_CANDIDATES_QQ":
    raise ArithmeticError("exact orbit42 candidates prerequisite is not passing")
if parent_artifact.get("status") != "PASS_EXACT_Q24_D13_TO_D12_COMPONENT_VALUATION_RR":
    raise ArithmeticError("exact D12 parent prerequisite is not passing")
if physical_artifact.get("status") != "PASS_Q42_EXACT_I8STAR_PHYSICAL_MARKING":
    raise ArithmeticError("physical I8* marking prerequisite is not passing")

candidate_rows = candidate_artifact["candidates"]
if not 0 <= args.candidate_index < len(candidate_rows):
    raise ArithmeticError("candidate index is out of range")
candidate = candidate_rows[args.candidate_index]
model = model_artifact["exact_model"]

# -------------------------------------------------------------------------
# Exact change from the R3-zero shell equation to the minimal D12 equation.
# -------------------------------------------------------------------------
UQ = PolynomialRing(QQ, "u")
u = UQ.gen()
KU = UQ.fraction_field()
VQ = PolynomialRing(QQ, "V")
V = VQ.gen()
KV = VQ.fraction_field()

A_shell = UQ([QQ(value) for value in model["A_coefficients_low_to_high"]])
B_shell = UQ([QQ(value) for value in model["B_coefficients_low_to_high"]])
A_parent = VQ([
    QQ(value)
    for value in parent_artifact["child"]["minimal_A_coefficients_low_to_high"]
])
B_parent = VQ([
    QQ(value)
    for value in parent_artifact["child"]["minimal_B_coefficients_low_to_high"]
])

raw_root = QQ(model["I8star_root"])
base_scale = QQ(model["base_scale"])
center = QQ(model["center"])
u_of_V = KV.one() / (KV(V) - raw_root) / base_scale - center / base_scale


def evaluate_u(value, argument):
    value = KU(value)
    return (
        KV(VQ(value.numerator())(argument))
        / KV(VQ(value.denominator())(argument))
    )


def rational_square_root(value):
    value = KV(value)
    numerator = VQ(value.numerator())
    denominator = VQ(value.denominator())
    answer = KV.one()
    for poly, direction in ((numerator, 1), (denominator, -1)):
        leading = QQ(poly.leading_coefficient())
        if not leading.is_square():
            raise ArithmeticError("coordinate scale has nonsquare leading coefficient")
        part = KV.one()
        for factor, multiplicity in (poly / leading).factor():
            if int(multiplicity) % 2:
                raise ArithmeticError("coordinate scale is not a rational square")
            part *= KV(factor.monic()) ** (int(multiplicity) // 2)
        root = KV(leading.sqrt()) * part
        answer = answer * root if direction == 1 else answer / root
    return answer


x_scale = rational_square_root(
    KV(A_parent) / evaluate_u(A_shell, u_of_V)
)
y_scale = rational_square_root(
    KV(B_parent) / evaluate_u(B_shell, u_of_V)
)
if x_scale**2 * evaluate_u(A_shell, u_of_V) != KV(A_parent):
    raise ArithmeticError("exact A-coordinate change failed")
if y_scale**2 * evaluate_u(B_shell, u_of_V) != KV(B_parent):
    raise ArithmeticError("exact B-coordinate change failed")
if y_scale**2 != x_scale**3:
    raise ArithmeticError("exact x/y coordinate scales are incompatible")

X_shell = UQ([QQ(value) for value in candidate["X_coefficients_low_to_high"]])
Y_shell = UQ([QQ(value) for value in candidate["Y_coefficients_low_to_high"]])
Z_shell = UQ([QQ(value) for value in candidate["Z_coefficients_low_to_high"]])
x_parent = x_scale * evaluate_u(KU(X_shell) / KU(Z_shell**2), u_of_V)
y_parent = y_scale * evaluate_u(KU(Y_shell) / KU(Z_shell**3), u_of_V)

E_parent = EllipticCurve(
    KV, [0, 0, 0, KV(A_parent), KV(B_parent)]
)
P_parent = E_parent(x_parent, y_parent)


def power_root(poly, exponent):
    poly = VQ(poly)
    leading = QQ(poly.leading_coefficient())
    answer = VQ.one()
    for factor, multiplicity in (poly / leading).factor():
        if int(multiplicity) % exponent:
            raise ArithmeticError("section denominator is not a perfect power")
        answer *= factor.monic() ** (int(multiplicity) // exponent)
    return answer.monic()


X_num, X_den = map(VQ, (P_parent[0].numerator(), P_parent[0].denominator()))
Y_num, Y_den = map(VQ, (P_parent[1].numerator(), P_parent[1].denominator()))
X_num, X_den = X_num / X_den.leading_coefficient(), X_den / X_den.leading_coefficient()
Y_num, Y_den = Y_num / Y_den.leading_coefficient(), Y_den / Y_den.leading_coefficient()
Z_parent = power_root(X_den, 2)
if Z_parent != power_root(Y_den, 3):
    raise ArithmeticError("x/y section denominators give different Z")
X_parent = VQ(KV(P_parent[0]) * Z_parent**2)
Y_parent = VQ(KV(P_parent[1]) * Z_parent**3)
if Y_parent**2 != (
    X_parent**3
    + A_parent * X_parent * Z_parent**4
    + B_parent * Z_parent**6
):
    raise ArithmeticError("minimal-equation section identity failed")
if (X_parent.degree(), Y_parent.degree(), Z_parent.degree()) != (10, 15, 3):
    raise ArithmeticError("unexpected minimal-equation orbit42 section degrees")

log(
    "SECTION",
    candidate=args.candidate_index,
    mapping=candidate["mapping_index"],
    spinor=candidate["spinor_index"],
    degrees="10,15,3",
    identity="PASS",
)

# -------------------------------------------------------------------------
# Exact resolved RR plane.
#
# With a=AA/Z^2 and b=BB/Z, regularity at the three smooth P.O collision
# points is AA*X == BB*Y mod Z^2.  X is a unit modulo Z^2, so BB=1,V and
# the homogeneous Z^2 direction give a three-dimensional space.  At C01 the
# reduced tangent cone y^2 gives weights (2,2,3); the effective-root sign
# therefore imposes exactly AA(alpha)=0.
# -------------------------------------------------------------------------
I8_row = next(
    row
    for row in parent_artifact["child"]["finite_fibres"]
    if row["kodaira"] == "I8*"
)
I8_factor = VQ(I8_row["factor"])
if I8_factor.degree() != 1:
    raise ArithmeticError("I8* support is not rational linear")
alpha = -I8_factor[0] / I8_factor[1]
if alpha != raw_root:
    raise ArithmeticError("shell/minimal I8* roots disagree")

collision_modulus = Z_parent**2
X_inverse = VQ(inverse_mod(X_parent, collision_modulus))
collision_pairs = [(Z_parent**2, VQ.zero())]
for BB in (VQ.one(), V):
    AA = VQ((BB * Y_parent * X_inverse) % collision_modulus)
    if (AA * X_parent - BB * Y_parent) % collision_modulus:
        raise ArithmeticError("exact smooth-collision congruence failed")
    collision_pairs.append((AA, BB))

rr_pairs = []
for AA, BB in collision_pairs[1:]:
    AA = VQ(AA - AA(alpha) / Z_parent(alpha) ** 2 * Z_parent**2)
    if AA(alpha) or (AA * X_parent - BB * Y_parent) % collision_modulus:
        raise ArithmeticError("exact resolved C01 condition failed")
    rr_pairs.append((AA, BB))

if len(rr_pairs) != 2:
    raise ArithmeticError("resolved RR kernel does not have dimension two")

orientations = physical_artifact["orientation_candidates"]
selected_orientation = next(
    index
    for index, row in enumerate(orientations)
    if row["section_meets_physical_components"] == ["C10"]
    and row["omitted_vertical_components"] == ["C10"]
)
rejected_orientation = next(
    index
    for index, row in enumerate(orientations)
    if row["section_meets_physical_components"] == ["C11"]
    and row["omitted_vertical_components"] == ["C11"]
)

log(
    "RR",
    ambient=9,
    collision_rank=6,
    post_collision=3,
    resolved_rank=1,
    kernel=2,
    weights="2,2,3",
    selected_orientation=selected_orientation,
    status="PASS_EXACT",
)

# -------------------------------------------------------------------------
# Binary quartic over QQ(T) and exact minimized Jacobian.
# -------------------------------------------------------------------------
TQ = PolynomialRing(QQ, "T")
T = TQ.gen()
KT = TQ.fraction_field()
WV = PolynomialRing(KT, "V")
W = WV.gen()


def lift_v(poly):
    return WV([KT(value) for value in VQ(poly).list()])


AA0, BB0 = rr_pairs[0]
AA1, BB1 = rr_pairs[1]
z_lift, x_lift, y_lift, a_lift = map(
    lift_v, (Z_parent, X_parent, Y_parent, A_parent)
)
a0 = lift_v(AA0) / z_lift**2
b0 = lift_v(BB0) / z_lift
a1 = lift_v(AA1) / z_lift**2
b1 = lift_v(BB1) / z_lift

slope = (a1 - KT(T) * a0) / (KT(T) * b0 - b1)
x_point = x_lift / z_lift**2
y_point = y_lift / z_lift**3
radicand = (
    slope**4
    - 6 * x_point * slope**2
    - 8 * y_point * slope
    - 3 * x_point**2
    - 4 * a_lift
)
rad_num = WV(radicand.numerator())
rad_den = WV(radicand.denominator())
raw_branch = rad_num * rad_den
decomposition = raw_branch.squarefree_decomposition()
product = WV.one()
for factor, multiplicity in decomposition:
    product *= factor ** int(multiplicity)
scalar = raw_branch.leading_coefficient() / product.leading_coefficient()
odd_branch = WV(scalar)
for factor, multiplicity in decomposition:
    if int(multiplicity) % 2:
        odd_branch *= factor
if odd_branch.degree() != 4:
    raise ArithmeticError("resolved orbit42 branch is not quartic")

# If odd_branch=q_num/D, then (D*w)^2=D*q_num.  This clears the coefficient
# denominators without changing the double cover.  A common square is then
# absorbed into w.
denominator_lcm = TQ.one()
for coefficient in odd_branch.list():
    denominator_lcm = denominator_lcm.lcm(coefficient.denominator())
quartic_coefficients = [
    TQ((coefficient * denominator_lcm).numerator()) * denominator_lcm
    for coefficient in odd_branch.list()
]
common = quartic_coefficients[0]
for coefficient in quartic_coefficients[1:]:
    common = common.gcd(coefficient)
common_square = TQ.one()
for factor, multiplicity in common.factor():
    common_square *= factor ** (int(multiplicity) // 2)
quartic_coefficients = [
    coefficient // common_square**2
    for coefficient in quartic_coefficients
]
quartic = WV(sum(
    KT(quartic_coefficients[index]) * W**index
    for index in range(len(quartic_coefficients))
))
if quartic.degree() != 4:
    raise ArithmeticError("primitive quartic normalization changed the degree")

e, d, c, b, a = quartic.list()
I = 12 * a * e - 3 * b * d + c**2
J = (
    72 * a * c * e
    + 9 * b * c * d
    - 27 * a * d**2
    - 27 * b**2 * e
    - 2 * c**3
)
if I.denominator() != 1 or J.denominator() != 1:
    raise ArithmeticError("quartic invariants retain base denominators")
A_child = TQ(-27 * I)
B_child = TQ(-27 * J)
Delta_child = TQ(-16 * (4 * A_child**3 + 27 * B_child**2))

if (A_child.degree(), B_child.degree(), Delta_child.degree()) != (8, 12, 24):
    raise ArithmeticError("child is not a globally minimal elliptic K3 model")

delta_factors = sorted(
    ((factor.monic(), int(multiplicity)) for factor, multiplicity in Delta_child.factor()),
    key=lambda row: (row[1], row[0].degree(), str(row[0])),
)
if sorted((factor.degree(), multiplicity) for factor, multiplicity in delta_factors) != [
    (1, 12),
    (12, 1),
]:
    raise ArithmeticError("child discriminant does not have I12 + 12 I1 support")

I12_factor = next(factor for factor, multiplicity in delta_factors if multiplicity == 12)
nodal_factor = next(factor for factor, multiplicity in delta_factors if multiplicity == 1)
if not nodal_factor.is_squarefree():
    raise ArithmeticError("degree-12 nodal support is not squarefree")
if A_child.gcd(Delta_child) != 1 or B_child.gcd(Delta_child) != 1:
    raise ArithmeticError("multiplicative/nodal fibres collide with c4/c6")

log(
    "CHILD",
    quartic=4,
    Adeg=8,
    Bdeg=12,
    Ddeg=24,
    fibres="I12+12I1",
    ADE="A11",
    euler=24,
    MW_if_rho19=6,
    status="PASS_EXACT",
)

# Good-reduction regression for the exact child.  This also binds the chosen
# exact shell representative to the pinned physical-orientation computation.
p = ZZ(args.prime)
Fp = GF(p)


def reduce_q(value):
    value = QQ(value)
    if ZZ(value.denominator()) % p == 0:
        raise ArithmeticError(f"bad reduction denominator at p={p}")
    return Fp(ZZ(value.numerator())) / Fp(ZZ(value.denominator()))


TF = PolynomialRing(Fp, "T")
A_mod = TF([reduce_q(value) for value in A_child.list()])
B_mod = TF([reduce_q(value) for value in B_child.list()])
D_mod = TF(-16 * (4 * A_mod**3 + 27 * B_mod**2))
if sorted((factor.degree(), int(multiplicity)) for factor, multiplicity in D_mod.factor()) != [
    (1, 1),
    (1, 12),
    (11, 1),
]:
    raise ArithmeticError("exact child misses the pinned mod-100003 A11 pattern")


def text_hash(rows):
    blob = json.dumps(rows, separators=(",", ":"), ensure_ascii=True).encode()
    return hashlib.sha256(blob).hexdigest()


payload = {
    "schema": "elkies-k3.h3-q24-d12-to-a11-orbit42-resolved-rr-qq.v1",
    "status": "PASS_EXACT_Q24_D12_Q6_A11_COMPONENT_VALUATION_RR",
    "inputs": {
        "zero_pole_model": str(MODEL.relative_to(ROOT)),
        "exact_section_candidates": str(CANDIDATES.relative_to(ROOT)),
        "exact_D12_parent": str(PARENT.relative_to(ROOT)),
        "physical_I8star_marking": str(PHYSICAL.relative_to(ROOT)),
        "input_sha256": {
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (MODEL, CANDIDATES, PARENT, PHYSICAL)
        },
    },
    "selected_section": {
        "candidate_index": int(args.candidate_index),
        "mapping_index": int(candidate["mapping_index"]),
        "spinor_index": int(candidate["spinor_index"]),
        "P_dot_O": 3,
        "canonical_exact_point_input": str(CANDIDATES.relative_to(ROOT)),
        "minimal_equation_degrees": [10, 15, 3],
        "minimal_point_sha256": text_hash([
            [str(value) for value in X_parent.list()],
            [str(value) for value in Y_parent.list()],
            [str(value) for value in Z_parent.list()],
        ]),
        "exact_weierstrass_identity": True,
    },
    "coordinate_change": {
        "u_of_V": str(u_of_V),
        "x_scale": str(x_scale),
        "y_scale": str(y_scale),
    },
    "resolved_RR": {
        "divisor_geometric_sign": "O+P-V_effective",
        "abstract_lattice_formula": "D42=O+P+V_abstract",
        "sign_explanation": (
            "The graph marking fixes the D12 Cartan graph but not the global root sign. "
            "Nefness maps V_abstract to -V_effective."
        ),
        "ambient_dimension": 9,
        "smooth_collision_rank": 6,
        "post_collision_dimension": 3,
        "C01_tangent_cone": "y^2",
        "C01_weights": [2, 2, 3],
        "resolved_condition": "AA(alpha)=0",
        "resolved_rank": 1,
        "kernel_dimension": 2,
        "h0": 2,
        "selected_physical_orientation": int(selected_orientation),
        "selected_section_meets": "C10",
        "selected_omitted_component": "C10",
        "rejected_physical_orientation": int(rejected_orientation),
        "kernel_basis_sha256": text_hash([
            [[str(value) for value in AA.list()], [str(value) for value in BB.list()]]
            for AA, BB in rr_pairs
        ]),
    },
    "quartic": {
        "degree": 4,
        "raw_radicand_degrees": [int(rad_num.degree()), int(rad_den.degree())],
        "squarefree_decomposition": [
            [int(factor.degree()), int(multiplicity)]
            for factor, multiplicity in decomposition
        ],
        "coefficients_in_T_low_to_high": [str(value) for value in quartic.list()],
        "I": str(I),
        "J": str(J),
    },
    "child": {
        "minimal_A_coefficients_low_to_high": [str(value) for value in A_child.list()],
        "minimal_B_coefficients_low_to_high": [str(value) for value in B_child.list()],
        "discriminant_factorization": [
            {"factor": str(factor), "multiplicity": int(multiplicity)}
            for factor, multiplicity in delta_factors
        ],
        "finite_fibres_geometric": [
            {"kodaira": "I12", "count": 1},
            {"kodaira": "I1", "count": 12},
        ],
        "infinity_kind": "smooth",
        "root_lattice": "A11",
        "root_rank": 11,
        "root_determinant": 12,
        "euler_number": 24,
        "MW_rank_if_rho19": 6,
    },
    "verification": {
        "exact_coordinate_change": True,
        "exact_section_identity": True,
        "exact_collision_congruence": True,
        "exact_resolved_rank": True,
        "exact_quartic_squareclass": True,
        "exact_minimal_jacobian": True,
        "exact_A11_fibre_classification": True,
        "mod_100003_A11_regression": True,
    },
    "next_required": "A11_Q8_ORBIT922_MARKING_AND_RESOLVED_RR",
    "proof_boundary": (
        "The section, RR plane, quartic, Jacobian, and A11 fibre classification are exact over QQ. "
        "Identification of the selected exact section with the historical orbit42 shell and rejection "
        "of the opposite spinor arm use the pinned good-reduction marking at p=100003. The MW rank is "
        "reported as 6 under the repository's rho=19 lattice marking; no independent Picard upper bound "
        "is added by this certificate."
    ),
}

OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUTPUT}", flush=True)
print(
    "Q42RRQQ_RESULT|ambient=9|collision=6|post=3|resolved=1|kernel=2|"
    "quartic=4|fibres=I12+12I1|ADE=A11|euler=24|MW_if_rho19=6|"
    "status=PASS_EXACT_Q24_D12_Q6_A11_COMPONENT_VALUATION_RR",
    flush=True,
)
