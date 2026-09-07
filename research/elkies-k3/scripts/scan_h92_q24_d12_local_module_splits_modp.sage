#!/usr/bin/env sage -python
"""Search simple two-generator local modules for the H3 q24 -> D12 pencil.

The globally correct reduced collision system is

    B dimension 16 --6 conditions--> dimension 10.

At the I9* point write

    m0(u) = m(u,0,0) = -yP/xP,
    q = m-m0 = (xP*y+yP*x)/(xP*(x-xP)).

For every pair (r,s), this script tests the local module

    M(r,s) = <u^r, u^s*q>.

For f=a+b*m=(a+b*m0)+b*q, membership means

    a+b*m0 == 0 mod u^r,
    b       == 0 mod u^s.

Only pairs whose local condition has rank 8 on the ten-dimensional collision
space are compiled. Their squarefree branch degree is then measured. Degree
3/4 is genus one; degree 18 is the already-rejected eight-ordinary-jet model
M(8,0). A D12 result remains diagnostic until the module is verified on the
actual resolved I9* component atlas.
"""

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
parser.add_argument("--max-order", type=int, default=12)
parser.add_argument(
    "--factor",
    action="store_true",
    help="factor each squarefree branch polynomial over GF(p)(V); slower",
)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
LOCAL = ROOT / "artifacts/local/elkies-k3"
GEN = ROOT / "artifacts/generated-results"
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"
exec(compile(CORE.read_text(), str(CORE), "exec"))

p = ZZ(args.prime)
F = GF(p)
R = PolynomialRing(F, "U")
U = R.gen()
K = R.fraction_field()


def red_q(value):
    value = QQ(value)
    denominator = ZZ(value.denominator())
    if denominator % p == 0:
        raise ZeroDivisionError(
            f"coefficient denominator {denominator} is divisible by {p}"
        )
    return F(ZZ(value.numerator())) / F(denominator)


def red_poly(values):
    return R([red_q(value) for value in values])


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

q8 = json.loads(Q8.read_text())
child = q8["child"]
Acurve = red_poly(child["minimal_A_coefficients_low_to_high"])
Bcurve = red_poly(child["minimal_B_coefficients_low_to_high"])
Delta = -16 * (4 * Acurve**3 + 27 * Bcurve**2)

EXACT = LOCAL / "q8-q24-horizontal-section-qq.json"
MOD = LOCAL / f"q24-degree46-direct-global-mod-{p}.json"
section_source = None
if EXACT.exists():
    exact = json.loads(EXACT.read_text())
    if exact.get("status") == "PASS_EXACT_Q24_HORIZONTAL_SECTION":
        sec = exact["section"]
        Z = red_poly(sec["Z_coefficients_low_to_high"])
        X = red_poly(sec["X_coefficients_low_to_high"])
        Y = red_poly(sec["Y_coefficients_low_to_high"])
        section_source = str(EXACT.relative_to(ROOT))

if section_source is None:
    if not MOD.exists():
        raise SystemExit(f"Missing q24 section: {EXACT} and {MOD}")
    modular = json.loads(MOD.read_text())
    if modular.get("status") != "PASS_MODULAR_Q24_FROM_DIRECT_DEGREE46_BRIDGE":
        raise SystemExit(f"{MOD} is not a passing q24 section artifact")
    sec = modular["section_mod_p"]
    Z = R([F(int(value)) for value in sec["Z_coefficients_low_to_high"]])
    X = R([F(int(value)) for value in sec["X_coefficients_low_to_high"]])
    Y = R([F(int(value)) for value in sec["Y_coefficients_low_to_high"]])
    section_source = str(MOD.relative_to(ROOT))

assert (Z.degree(), X.degree(), Y.degree()) == (24, 52, 78)
assert Z.is_monic()
assert X.gcd(Z).degree() == 0
assert Y.gcd(Z).degree() == 0
assert Y**2 == X**3 + Acurve * X * Z**4 + Bcurve * Z**6

xP = K(X) / K(Z**2)
yP = K(Y) / K(Z**3)
assert yP**2 == xP**3 + K(Acurve) * xP + K(Bcurve)

# ---------------------------------------------------------------------------
# Reduced smooth-collision system: 16 B coefficients -> dimension 10.
# ---------------------------------------------------------------------------
modulus = Z**2
Xinv = X.inverse_mod(modulus)
A_columns = [((U**i) * Y * Xinv) % modulus for i in range(16)]
H = matrix(F, 6, 16, lambda row, col: A_columns[col][42 + row])
assert H.rank() == 6
K10 = H.right_kernel().basis_matrix()
assert K10.dimensions() == (10, 16)


def AB_from_Brow(row):
    Bcoef = R(list(row))
    Acoef = (Bcoef * Y * Xinv) % modulus
    assert Acoef.degree() <= 41
    assert (Acoef * X - Bcoef * Y) % modulus == 0
    return Acoef, Bcoef, K(Acoef) / K(Z**2), K(Bcoef) / K(Z)


pairs10 = [AB_from_Brow(row) for row in K10.rows()]

# ---------------------------------------------------------------------------
# Local I9* expansions.
# ---------------------------------------------------------------------------
i9 = next(item for item in child["finite_fibres"] if item["kodaira"] == "I9*")
RQ = PolynomialRing(QQ, "U")
fQ = RQ(str(i9["factor"]))
f = R([red_q(coefficient) for coefficient in fQ.list()])
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
    return [remainder[index] for index in range(length)]


assert (
    shift_poly(Acurve).valuation(),
    shift_poly(Bcurve).valuation(),
    shift_poly(Delta).valuation(),
) == (2, 3, 15)

xloc = shift_rf(xP)
yloc = shift_rf(yP)
assert xloc.valuation() == 0 and F(xloc(0)) != 0
assert yloc.valuation() >= 0
m0 = -yloc / xloc
assert m0.valuation() >= 0

max_order = ZZ(args.max_order)
if max_order < 0:
    raise ValueError("--max-order must be nonnegative")

scalar_local = []
q_coefficient_local = []
for unused_A, unused_B, a, b in pairs10:
    bloc = shift_rf(b)
    assert bloc.valuation() >= 0
    scalar_local.append(shift_rf(a) + bloc * m0)
    q_coefficient_local.append(bloc)

scalar_jets = [jet(value, max_order) for value in scalar_local]
coefficient_jets = [jet(value, max_order) for value in q_coefficient_local]

# New-base rings used for candidate branch calculations.
VR = PolynomialRing(F, "V")
V = VR.gen()
VF = VR.fraction_field()
UR = PolynomialRing(VF, "U")
UK = UR.fraction_field()


def lift_poly(poly):
    poly = R(poly)
    return UR([VF(coefficient) for coefficient in poly.list()])


def lift_rf(value):
    value = K(value)
    return UK(lift_poly(R(value.numerator()))) / UK(
        lift_poly(R(value.denominator()))
    )


xPV = lift_rf(xP)
yPV = lift_rf(yP)
AV = lift_poly(Acurve)
BV = lift_poly(Bcurve)


def canonical_two_space(B2):
    echelon = matrix(F, B2).echelon_form()
    if echelon.dimensions() != (2, 16) or echelon.rank() != 2:
        raise ArithmeticError("candidate B-space is not two-dimensional")
    return echelon, tuple(int(value) for value in echelon.list())


def compile_candidate(B2):
    pairs2 = [AB_from_Brow(row) for row in B2.rows()]
    a0, b0 = lift_rf(pairs2[0][2]), lift_rf(pairs2[0][3])
    a1, b1 = lift_rf(pairs2[1][2]), lift_rf(pairs2[1][3])

    # Our chord is (y+yP)/(x-xP), i.e. the standard chord through -P.
    mvalue = pencil_chord_solution(a0, b0, a1, b1, VF(V))
    radicand = chord_discriminant(xPV, -yPV, AV, mvalue)
    quartic, square_factor = squarefree_binary_quartic(radicand, UR)
    degree = int(quartic.degree())

    record = {
        "squarefree_degree": degree,
        "radicand_numerator_degree": int(UR(radicand.numerator()).degree()),
        "radicand_denominator_degree": int(UR(radicand.denominator()).degree()),
        "quartic_coefficients_low_to_high": [str(value) for value in quartic.list()],
        "factor_degrees": [],
        "root_rank": None,
        "root_determinant": None,
        "euler": None,
        "D12": False,
    }

    if args.factor:
        record["factor_degrees"] = [
            [int(factor.degree()), int(exponent)]
            for factor, exponent in quartic.factor()
        ]

    if degree not in (3, 4):
        return record

    invariant_i, invariant_j = binary_quartic_invariants(quartic)
    jacA = VF(-27) * VF(invariant_i)
    jacB = VF(-27) * VF(invariant_j)
    classification = classify_finite_short_weierstrass_fibres(VR, jacA, jacB)

    root_rank = int(classification["finite_root_rank"])
    root_det = int(classification["finite_root_determinant"])
    euler = int(classification["finite_euler_number"])
    infinity_orders = tuple(
        map(int, classification["infinity_boundary"]["normalized_orders"])
    )
    infinity_kind = "smooth"
    if infinity_orders[2] > 0:
        inf_rank, inf_euler, inf_det, infinity_kind = (
            kodaira_data_from_short_orders(*infinity_orders)
        )
        root_rank += int(inf_rank)
        root_det *= int(inf_det)
        euler += int(inf_euler)

    record.update(
        {
            "jacobian_raw": {"A": str(jacA), "B": str(jacB)},
            "finite_fibres": [
                {
                    "factor": str(item["factor"]),
                    "degree": int(item["degree"]),
                    "minimal_orders": list(map(int, item["minimal_orders"])),
                    "kodaira": item["kodaira"],
                }
                for item in classification["finite_fibres"]
            ],
            "infinity_orders": list(infinity_orders),
            "infinity_kind": infinity_kind,
            "root_rank": root_rank,
            "root_determinant": root_det,
            "euler": euler,
            "D12": (root_rank, root_det) == (12, 4),
        }
    )
    return record


spaces = {}
all_pairs = []
for r in range(int(max_order) + 1):
    for s in range(int(max_order) + 1):
        rows = []
        row_labels = []
        for order in range(r):
            rows.append([scalar_jets[column][order] for column in range(10)])
            row_labels.append(["scalar", order])
        for order in range(s):
            rows.append(
                [coefficient_jets[column][order] for column in range(10)]
            )
            row_labels.append(["q_coefficient", order])

        local = matrix(F, rows) if rows else matrix(F, 0, 10)
        rank = int(local.rank())
        kernel_dimension = 10 - rank
        pair_record = {
            "r": r,
            "s": s,
            "condition_rows": len(rows),
            "rank": rank,
            "kernel_dimension": kernel_dimension,
            "row_labels": row_labels,
            "space_id": None,
        }
        all_pairs.append(pair_record)

        if kernel_dimension != 2:
            continue

        K2_in_K10 = local.right_kernel().basis_matrix()
        if K2_in_K10.dimensions() != (2, 10):
            raise ArithmeticError("rank-eight local block has wrong kernel")
        B2 = K2_in_K10 * K10
        canonical, key = canonical_two_space(B2)

        if key not in spaces:
            candidate = compile_candidate(canonical)
            candidate.update(
                {
                    "space_id": len(spaces) + 1,
                    "B_echelon_rows": [
                        [int(value) for value in row]
                        for row in canonical.rows()
                    ],
                    "modules": [],
                }
            )
            spaces[key] = candidate

        candidate = spaces[key]
        pair_record["space_id"] = candidate["space_id"]
        candidate["modules"].append([r, s])

        print(
            "Q24MODSCAN_CANDIDATE|"
            f"r={r}|s={s}|space={candidate['space_id']}|"
            f"branch={candidate['squarefree_degree']}|"
            f"root_rank={candidate['root_rank'] if candidate['root_rank'] is not None else 'NA'}|"
            f"root_det={candidate['root_determinant'] if candidate['root_determinant'] is not None else 'NA'}|"
            f"D12={int(candidate['D12'])}|"
            f"status={'D12' if candidate['D12'] else 'GENUS_ONE_OTHER' if candidate['squarefree_degree'] in (3,4) else 'NOT_GENUS_ONE'}",
            flush=True,
        )

candidates = sorted(
    spaces.values(),
    key=lambda item: (
        0 if item["D12"] else 1,
        abs(item["squarefree_degree"] - 4),
        item["squarefree_degree"],
        item["space_id"],
    ),
)
d12_candidates = [item for item in candidates if item["D12"]]
genus_one_candidates = [
    item for item in candidates if item["squarefree_degree"] in (3, 4)
]

rejected_80 = next(
    (item for item in candidates if [8, 0] in item["modules"]),
    None,
)
if rejected_80 is None:
    raise ArithmeticError("M(8,0) did not produce a two-dimensional candidate")
if rejected_80["squarefree_degree"] != 18:
    raise ArithmeticError(
        "M(8,0) no longer reproduces the observed degree-18 rejection"
    )

payload = {
    "schema": "elkies-k3.h3-q24-local-module-split-scan-modp.v1",
    "status": (
        "FOUND_D12_LOCAL_MODULE_CANDIDATE"
        if d12_candidates
        else "FOUND_GENUS_ONE_LOCAL_MODULE_CANDIDATE"
        if genus_one_candidates
        else "NO_GENUS_ONE_SPLIT_MODULE_UP_TO_BOUND"
    ),
    "proof_boundary": (
        "This is a modular discovery scan in the verified reduced collision "
        "space. A successful module still requires valuation verification on "
        "the actual resolved I9* component atlas before it can certify h0(D)."
    ),
    "prime": int(p),
    "max_order": int(max_order),
    "inputs": {
        "q8_child": str(Q8.relative_to(ROOT)),
        "q24_section": section_source,
    },
    "collision": {
        "B_dimension": 16,
        "rank": 6,
        "kernel_dimension": 10,
    },
    "rejected_eight_jet_model": {
        "module": [8, 0],
        "space_id": rejected_80["space_id"],
        "squarefree_degree": rejected_80["squarefree_degree"],
    },
    "pair_tests": all_pairs,
    "candidate_spaces": candidates,
    "D12_candidate_count": len(d12_candidates),
    "genus_one_candidate_count": len(genus_one_candidates),
}
OUT = (
    args.output.resolve()
    if args.output
    else LOCAL / f"q24-local-module-split-scan-mod-{p}.json"
)
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(f"OUTPUT|{OUT}", flush=True)
print(
    "Q24MODSCAN_RESULT|"
    f"prime={p}|max_order={max_order}|"
    f"spaces={len(candidates)}|genus_one={len(genus_one_candidates)}|"
    f"D12={len(d12_candidates)}|"
    f"best_branch={candidates[0]['squarefree_degree'] if candidates else 'NA'}|"
    f"status={payload['status']}",
    flush=True,
)
