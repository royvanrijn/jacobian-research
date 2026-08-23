#!/usr/bin/env sage -python
'''Test the H3 q24 candidate local module using a reduced 16-variable RR system.

The smooth-collision condition is

    A*X - B*Y == 0 mod Z^2,

with deg(A)<=41 and deg(B)<=15. Since X is invertible modulo Z^2, A is uniquely
determined by B. Only the six coefficients of degrees 42..47 need to vanish, so
the global system reduces from 58 variables to 16 variables:

    16 - 6 = 10.

At the I9* place put

    m  = (y+yP)/(x-xP),
    m0 = m(u,0,0) = -yP/xP,
    q  = m-m0 = (xP*y+yP*x)/(xP*(x-xP)).

The proposed local module is <u^8,q>. Within the marked-chord ambient this is
equivalent to requiring a+b*m0 to vanish modulo u^8. The modular mode also
proves row-space equality with the historical 58-variable eight-cusp-jet
system. The --exact mode solves the same 16 -> 10 -> 2 problem over QQ and
attempts the exact quartic/Jacobian compilation.
'''

import argparse
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix


def locate_repo(explicit=None):
    candidates = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd = Path.cwd().resolve()
    candidates += [cwd, *cwd.parents]
    home = Path.home()
    candidates += [
        home / "Documents" / "jacobian-research",
        home / "jacobian-research",
        home / "src" / "jacobian-research",
        home / "git" / "jacobian-research",
        home / "projects" / "jacobian-research",
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
        if (
            (candidate / "elkies-k3/scripts").is_dir()
            and (candidate / "artifacts/generated-results").is_dir()
        ):
            return candidate
    raise SystemExit("Could not locate jacobian-research; pass --repo PATH")


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--prime", type=int, default=100003)
parser.add_argument("--exact", action="store_true", help="solve over QQ")
parser.add_argument(
    "--rr-only",
    action="store_true",
    help="stop after constructing the exact/modular two-row RR pencil",
)
parser.add_argument("--output", type=Path)
parser.add_argument(
    "--no-assert-d12",
    action="store_true",
    help="write diagnostics instead of failing if the candidate child is not D12",
)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
LOCAL = ROOT / "artifacts/local/elkies-k3"
GEN = ROOT / "artifacts/generated-results"
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"
if not CORE.exists():
    raise SystemExit(f"Missing compiler core: {CORE}")
exec(compile(CORE.read_text(), str(CORE), "exec"))

q8_candidates = [
    LOCAL / "q8-corrected2cover-qq-child.json",
    GEN / "elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
]
Q8 = next(
    (
        path
        for path in q8_candidates
        if path.exists()
        and json.loads(path.read_text()).get("status")
        == "PASS_EXACT_CORRECTED_Q8_D13_CHILD"
    ),
    None,
)
if Q8 is None:
    raise SystemExit("No passing exact corrected D13 q8 child artifact")

if args.exact:
    p = None
    F = QQ
    field_label = "QQ"
else:
    p = ZZ(args.prime)
    F = GF(p)
    field_label = f"GF({p})"

R = PolynomialRing(F, "U")
U = R.gen()
K = R.fraction_field()


def base_coefficient(value):
    value = QQ(value)
    if args.exact:
        return value
    denominator = ZZ(value.denominator())
    if denominator % p == 0:
        raise ZeroDivisionError(
            f"coefficient denominator {denominator} is divisible by {p}"
        )
    return F(ZZ(value.numerator())) / F(denominator)


def base_poly(values):
    return R([base_coefficient(value) for value in values])


q8 = json.loads(Q8.read_text())
child = q8["child"]
Acurve = base_poly(child["minimal_A_coefficients_low_to_high"])
Bcurve = base_poly(child["minimal_B_coefficients_low_to_high"])
Delta = -16 * (4 * Acurve**3 + 27 * Bcurve**2)

EXACT = LOCAL / "q8-q24-horizontal-section-qq.json"
section_source = None
if EXACT.exists():
    exact = json.loads(EXACT.read_text())
    if exact.get("status") == "PASS_EXACT_Q24_HORIZONTAL_SECTION":
        sec = exact["section"]
        Z = base_poly(sec["Z_coefficients_low_to_high"])
        X = base_poly(sec["X_coefficients_low_to_high"])
        Y = base_poly(sec["Y_coefficients_low_to_high"])
        section_source = str(EXACT.relative_to(ROOT))

if section_source is None:
    if args.exact:
        raise SystemExit(
            f"Exact mode requires {EXACT} with status PASS_EXACT_Q24_HORIZONTAL_SECTION"
        )
    MOD = LOCAL / f"q24-degree46-direct-global-mod-{p}.json"
    if not MOD.exists():
        raise SystemExit(f"Missing q24 section: {EXACT} and {MOD}")
    modular = json.loads(MOD.read_text())
    assert modular["status"] == "PASS_MODULAR_Q24_FROM_DIRECT_DEGREE46_BRIDGE"
    sec = modular["section_mod_p"]
    Z = R([F(int(v)) for v in sec["Z_coefficients_low_to_high"]])
    X = R([F(int(v)) for v in sec["X_coefficients_low_to_high"]])
    Y = R([F(int(v)) for v in sec["Y_coefficients_low_to_high"]])
    section_source = str(MOD.relative_to(ROOT))

assert (Z.degree(), X.degree(), Y.degree()) == (24, 52, 78)
assert Z.is_monic()
assert X.gcd(Z).degree() == 0
assert Y.gcd(Z).degree() == 0
assert Y**2 == X**3 + Acurve * X * Z**4 + Bcurve * Z**6

xP = K(X) / K(Z**2)
yP = K(Y) / K(Z**3)
assert yP**2 == xP**3 + K(Acurve) * xP + K(Bcurve)

print(
    "Q24RED16_INPUT|"
    f"field={field_label}|section={section_source}|"
    "Z=24|X=52|Y=78|status=PASS",
    flush=True,
)

# ---------------------------------------------------------------------------
# 1. Eliminate A from the smooth-collision block.
# ---------------------------------------------------------------------------
modulus = Z**2
assert modulus.degree() == 48
Xinv = X.inverse_mod(modulus)

A_columns = [((U**i) * Y * Xinv) % modulus for i in range(16)]
high_rows = tuple(range(42, 48))
H = matrix(
    F,
    6,
    16,
    lambda row, col: A_columns[col][high_rows[row]],
)
high_rank = H.rank()
K10 = H.right_kernel().basis_matrix()
print(
    "Q24RED16_COLLISION|"
    f"field={field_label}|Bcols=16|rows=6|rank={high_rank}|"
    f"kernel={K10.nrows()}|"
    f"status={'PASS' if (high_rank, K10.nrows()) == (6, 10) else 'UNEXPECTED'}",
    flush=True,
)
if (high_rank, K10.nrows()) != (6, 10):
    raise ArithmeticError("reduced collision system is not 16 -> 10")


def AB_from_Brow(row):
    Bcoef = R(list(row))
    Acoef = (Bcoef * Y * Xinv) % modulus
    assert Acoef.degree() <= 41
    assert (Acoef * X - Bcoef * Y) % modulus == 0
    return Acoef, Bcoef, K(Acoef) / K(Z**2), K(Bcoef) / K(Z)


pairs10 = [AB_from_Brow(row) for row in K10.rows()]

# ---------------------------------------------------------------------------
# 2. Candidate local module <u^8,m-m0>.
# ---------------------------------------------------------------------------
i9 = next(item for item in child["finite_fibres"] if item["kodaira"] == "I9*")
RQ = PolynomialRing(QQ, "U")
fQ = RQ(str(i9["factor"]))
f = R([base_coefficient(c) for c in fQ.list()])
assert f.degree() == 1
alpha = -f[0] / f[1]

S = PolynomialRing(F, "u")
u = S.gen()
KS = S.fraction_field()


def shift_poly(poly):
    return S(R(poly)(alpha + u))


def shift_rf(value):
    value = K(value)
    numerator = shift_poly(R(value.numerator()))
    denominator = shift_poly(R(value.denominator()))
    if not denominator[0]:
        raise ZeroDivisionError("local denominator is not a unit at I9*")
    return KS(numerator) / KS(denominator)


def jet(value, length):
    value = KS(value)
    numerator = S(value.numerator())
    denominator = S(value.denominator())
    if not denominator[0]:
        raise ZeroDivisionError("jet denominator is not a unit")
    remainder = (numerator * denominator.inverse_mod(u**length)) % (u**length)
    return [remainder[i] for i in range(length)]


assert (
    shift_poly(Acurve).valuation(),
    shift_poly(Bcurve).valuation(),
    shift_poly(Delta).valuation(),
) == (2, 3, 15)

xloc = shift_rf(xP)
yloc = shift_rf(yP)
assert xloc.valuation() >= 0 and yloc.valuation() >= 0
x0 = F(xloc(0))
y0 = F(yloc(0))
assert x0 or y0
assert y0**2 == x0**3
assert x0 != 0
m0 = -yloc / xloc
assert m0.valuation() >= 0

# Exact local algebra identity.
XY = PolynomialRing(KS, names=("x", "y"))
xpoly, ypoly = XY.gens()
XYF = XY.fraction_field()
xx, yy = XYF(xpoly), XYF(ypoly)
xlocF, ylocF, m0F = XYF(xloc), XYF(yloc), XYF(m0)
mexpr = (yy + ylocF) / (xx - xlocF)
Lexpr = xlocF * yy + ylocF * xx
assert mexpr - m0F == Lexpr / (xlocF * (xx - xlocF))


def local_scalar(pair):
    unused_A, unused_B, a, b = pair
    bloc = shift_rf(b)
    assert bloc.valuation() >= 0
    return shift_rf(a) + bloc * m0


local10 = [local_scalar(pair) for pair in pairs10]
J8_values = [jet(value, 8) for value in local10]
J8 = matrix(F, 8, 10, lambda row, col: J8_values[col][row])
J8_rank = J8.rank()
K2_in_K10 = J8.right_kernel().basis_matrix()
B2 = K2_in_K10 * K10
if (J8_rank, K2_in_K10.nrows(), B2.rank()) != (8, 2, 2):
    raise ArithmeticError("candidate local module does not cut 10 -> 2")

print(
    "Q24RED16_LOCAL|"
    f"field={field_label}|module=<u^8,m-m0>|"
    "identity=(xP*y+yP*x)/(xP*(x-xP))|"
    f"rank={J8_rank}|kernel=2|status=PASS",
    flush=True,
)

# Ninth jet selects the distinguished Theta line.
J9_values = [jet(value, 9) for value in local10]
J9 = matrix(F, 9, 10, lambda row, col: J9_values[col][row])
Ktheta_in_K10 = J9.right_kernel().basis_matrix()
if (J9.rank(), Ktheta_in_K10.nrows()) != (9, 1):
    raise ArithmeticError("ninth scalar jet did not select one line")
Btheta = Ktheta_in_K10 * K10

# ---------------------------------------------------------------------------
# 3. In modular mode, prove exact equivalence with the original 58 columns.
# ---------------------------------------------------------------------------
equivalence = {
    "checked": False,
    "same_B_rowspace": None,
    "same_theta_line": None,
}
if not args.exact:
    ambient58 = [("A", i) for i in range(42)] + [("B", i) for i in range(16)]
    collision_columns = [
        (U**i * X) % modulus if kind == "A" else (-U**i * Y) % modulus
        for kind, i in ambient58
    ]
    Cfull = matrix(F, 48, 58, lambda row, col: collision_columns[col][row])
    assert Cfull.rank() == 48

    Jcolumns8 = []
    Jcolumns9 = []
    for kind, i in ambient58:
        if kind == "A":
            scalar = shift_rf(K(U**i) / K(Z**2))
        else:
            scalar = shift_rf(K(U**i) / K(Z)) * m0
        Jcolumns8.append(jet(scalar, 8))
        Jcolumns9.append(jet(scalar, 9))

    Mfull8 = Cfull.stack(
        matrix(F, 8, 58, lambda row, col: Jcolumns8[col][row])
    )
    Kfull2 = Mfull8.right_kernel().basis_matrix()
    assert Mfull8.rank() == 56 and Kfull2.dimensions() == (2, 58)
    Bprojection2 = matrix(F, [list(row[42:58]) for row in Kfull2.rows()])
    same_two_space = Bprojection2.stack(B2).rank() == 2
    assert same_two_space

    Mfull9 = Cfull.stack(
        matrix(F, 9, 58, lambda row, col: Jcolumns9[col][row])
    )
    Kfulltheta = Mfull9.right_kernel().basis_matrix()
    assert Mfull9.rank() == 57 and Kfulltheta.dimensions() == (1, 58)
    Bprojection_theta = matrix(
        F, [list(row[42:58]) for row in Kfulltheta.rows()]
    )
    same_theta_line = Bprojection_theta.stack(Btheta).rank() == 1
    assert same_theta_line

    equivalence = {
        "checked": True,
        "full_rank": int(Mfull8.rank()),
        "full_kernel_dimension": int(Kfull2.nrows()),
        "same_B_rowspace": True,
        "same_theta_line": True,
    }
    print(
        "Q24RED16_EQUIV|full58_rank=56|kernel=2|"
        "same_B_rowspace=1|same_theta_line=1|status=PASS",
        flush=True,
    )

pairs2 = [AB_from_Brow(row) for row in B2.rows()]


def coefficient_list(poly):
    poly = R(poly)
    if args.exact:
        return [str(c) for c in poly.list()]
    return [int(c) for c in poly.list()]


base_payload = {
    "schema": "elkies-k3.h3-q24-d12-reduced16.v1",
    "field": field_label,
    "mode": "exact" if args.exact else "modular",
    "prime": None if args.exact else int(p),
    "inputs": {
        "q8_child": str(Q8.relative_to(ROOT)),
        "q24_section": section_source,
    },
    "reduced_collision": {
        "B_dimension": 16,
        "high_degree_rows": list(high_rows),
        "rank": int(high_rank),
        "kernel_dimension": int(K10.nrows()),
    },
    "local_module": {
        "module": "<u^8,m-m0>",
        "identity": "m-m0=(xP*y+yP*x)/(xP*(x-xP))",
        "xP_is_unit": True,
        "rank": int(J8_rank),
        "kernel_dimension": int(K2_in_K10.nrows()),
    },
    "equivalence": equivalence,
    "pencil": [
        {
            "A_coefficients_low_to_high": coefficient_list(Acoef),
            "B_coefficients_low_to_high": coefficient_list(Bcoef),
        }
        for Acoef, Bcoef, unused_a, unused_b in pairs2
    ],
    "theta": {
        "A_coefficients_low_to_high": coefficient_list(
            AB_from_Brow(Btheta[0])[0]
        ),
        "B_coefficients_low_to_high": coefficient_list(R(list(Btheta[0]))),
    },
    "proof_boundary": (
        "This constructs the 2D space cut out by the candidate local module. "
        "Promotion still requires proving the component valuations of "
        "(u^8,xP*y+yP*x) on the actual resolved I9* atlas."
    ),
}

if args.rr_only:
    base_payload["status"] = "PASS_REDUCED16_RR_PENCIL"
    OUT = (
        args.output.resolve()
        if args.output
        else LOCAL
        / ("q24-d12-reduced16-qq-rr.json" if args.exact else f"q24-d12-reduced16-mod-{p}-rr.json")
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(base_payload, indent=2, sort_keys=True) + "\n")
    print(f"OUTPUT|{OUT}", flush=True)
    print(
        f"Q24RED16_RESULT|field={field_label}|16-6=10|10-8=2|"
        "stage=RR|status=PASS_REDUCED16_RR_PENCIL",
        flush=True,
    )
    raise SystemExit(0)

# ---------------------------------------------------------------------------
# 4. Compile the two-row pencil.
# ---------------------------------------------------------------------------
VR = PolynomialRing(F, "V")
V = VR.gen()
VF = VR.fraction_field()
UR = PolynomialRing(VF, "U")
UK = UR.fraction_field()


def lift_poly(poly):
    poly = R(poly)
    return UR([VF(c) for c in poly.list()])


def lift_rf(value):
    value = K(value)
    return UK(lift_poly(R(value.numerator()))) / UK(
        lift_poly(R(value.denominator()))
    )


a0, b0 = lift_rf(pairs2[0][2]), lift_rf(pairs2[0][3])
a1, b1 = lift_rf(pairs2[1][2]), lift_rf(pairs2[1][3])
xPV, yPV = lift_rf(xP), lift_rf(yP)
AV, BV = lift_poly(Acurve), lift_poly(Bcurve)

denominator = b1 - VF(V) * b0
if not denominator:
    raise ArithmeticError("both pencil generators have zero chord coefficient")
mvalue = -(a1 - VF(V) * a0) / denominator

XR = PolynomialRing(UK, "x")
xvar = XR.gen()
yline = XR(mvalue) * (xvar - XR(xPV)) - XR(yPV)
relation = yline**2 - xvar**3 - XR(AV) * xvar - XR(BV)
quadratic, remainder = relation.quo_rem(xvar - XR(xPV))
assert not remainder and quadratic.degree() == 2
discriminant = UK(quadratic[1] ** 2 - 4 * quadratic[2] * quadratic[0])

quartic, square_factor = squarefree_binary_quartic(discriminant, UR)
quartic_degree = int(quartic.degree())
if quartic_degree not in (3, 4):
    base_payload.update(
        {
            "status": "REJECT_REDUCED16_LOCAL_MODULE_NOT_GENUS_ONE",
            "quartic_degree": quartic_degree,
            "diagnosis": (
                "The reduced 16-variable calculation is correct and exactly "
                "matches the provisional eight ordinary cusp jets, but that "
                "two-space has squarefree branch degree {} and therefore is "
                "not the q24 genus-one pencil."
            ).format(quartic_degree),
        }
    )
    OUT = (
        args.output.resolve()
        if args.output
        else LOCAL
        / (
            "q24-d12-reduced16-qq-rejected.json"
            if args.exact
            else f"q24-d12-reduced16-mod-{p}-rejected.json"
        )
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(base_payload, indent=2, sort_keys=True) + "\n")
    print(
        "Q24RED16_QUARTIC|"
        f"degree={quartic_degree}|expected=3or4|"
        "status=REJECT_LOCAL_MODULE",
        flush=True,
    )
    print(f"OUTPUT|{OUT}", flush=True)
    print(
        "Q24RED16_RESULT|"
        f"field={field_label}|16-6=10|10-8=2|"
        f"quartic={quartic_degree}|"
        "status=REJECT_REDUCED16_LOCAL_MODULE_NOT_GENUS_ONE",
        flush=True,
    )
    raise SystemExit(3)
I, J = binary_quartic_invariants(quartic)
jacA = VF(-27) * VF(I)
jacB = VF(-27) * VF(J)
classification = classify_finite_short_weierstrass_fibres(VR, jacA, jacB)

finite = [
    {
        "factor": str(item["factor"]),
        "degree": int(item["degree"]),
        "minimal_orders": list(map(int, item["minimal_orders"])),
        "kodaira": item["kodaira"],
    }
    for item in classification["finite_fibres"]
]
root_rank = int(classification["finite_root_rank"])
euler = int(classification["finite_euler_number"])
root_det = int(classification["finite_root_determinant"])
inf_orders = tuple(
    map(int, classification["infinity_boundary"]["normalized_orders"])
)
inf_kind = "smooth"
if inf_orders[2] > 0:
    inf_rank, inf_euler, inf_det, inf_kind = kodaira_data_from_short_orders(
        *inf_orders
    )
    root_rank += int(inf_rank)
    euler += int(inf_euler)
    root_det *= int(inf_det)

is_d12 = (root_rank, root_det) == (12, 4)
base_payload.update(
    {
        "status": (
            "PASS_REDUCED16_EXACT_D12_CANDIDATE"
            if args.exact and is_d12
            else "PASS_REDUCED16_EQUIVALENT_PROVISIONAL_D12"
            if is_d12
            else "REDUCED16_CHILD_NOT_D12"
        ),
        "quartic_degree": int(quartic.degree()),
        "quartic_coefficients_low_to_high": [str(c) for c in quartic.list()],
        "jacobian_raw": {"A": str(jacA), "B": str(jacB)},
        "child": {
            "finite_fibres": finite,
            "infinity_orders": list(inf_orders),
            "infinity_kind": inf_kind,
            "root_rank": int(root_rank),
            "root_determinant": int(root_det),
            "euler": int(euler),
        },
    }
)

OUT = (
    args.output.resolve()
    if args.output
    else LOCAL
    / ("q24-d12-reduced16-qq.json" if args.exact else f"q24-d12-reduced16-mod-{p}.json")
)
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(base_payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUT}", flush=True)
print(
    "Q24RED16_RESULT|"
    f"field={field_label}|16-6=10|10-8=2|quartic={quartic.degree()}|"
    f"root_rank={root_rank}|root_det={root_det}|euler={euler}|"
    f"D12={int(is_d12)}|status={base_payload['status']}",
    flush=True,
)
if not is_d12 and not args.no_assert_d12:
    raise SystemExit(2)
