#!/usr/bin/env sage -python
'''
Search non-split local HNF modules for the H3 q24 -> D12 neighbour.

We already know:
  * smooth collision: 16 B coefficients -> rank 6 -> dimension 10;
  * any correct remaining local block must cut this to dimension 2;
  * split modules <u^r, u^s q>, r+s=8, give branch degrees
      20,16,12,8,4,6,10,14,18,
    and the sole genus-one split (r,s)=(4,4) is a D13 self-neighbour.

Here
    m  = (y+yP)/(x-xP)
    m0 = m(u,0,0) = -yP/xP
    q  = m-m0.

A general full-rank local submodule can be put in HNF form

    M(r,s,c) = < u^r, c(u) + u^s q >,

with r+s=8 and c considered modulo u^r.

For f = a + b*m = scalar + b*q, scalar=a+b*m0, membership is

    b = u^s h
    scalar - c*h == 0 mod u^r.

For each tested c we:
  1. build the local condition on the verified 10D post-collision space;
  2. require rank 8 / kernel dimension 2;
  3. compile that exact 2-space to its chord double cover;
  4. measure the squarefree branch degree;
  5. if degree 3/4, classify its Jacobian and test D12 root data.

This is a discovery scan, not a resolved-I9* proof.
'''

import argparse
import json
import random
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
parser.add_argument(
    "--random-per-split",
    type=int,
    default=24,
    help="number of deterministic random c(u) polynomials per r+s=8 split",
)
parser.add_argument(
    "--small-lambda",
    type=int,
    default=8,
    help="test monomial shears lambda*u^k for +/-lambda up to this bound",
)
parser.add_argument("--factor", action="store_true")
parser.add_argument("--stop-on-d12", action="store_true")
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
# Reduced smooth collision: 16 B coefficients -> dimension 10.
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
# Actual I9* local coordinate.
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
    if length == 0:
        return []
    value = KS(value)
    numerator = S(value.numerator())
    denominator = S(value.denominator())
    if not denominator[0]:
        raise ZeroDivisionError("jet denominator is not a unit")
    remainder = (numerator * denominator.inverse_mod(u**length)) % (u**length)
    return [remainder[index] for index in range(length)]


Al = shift_poly(Acurve)
Bl = shift_poly(Bcurve)
Dl = shift_poly(Delta)
assert (Al.valuation(), Bl.valuation(), Dl.valuation()) == (2, 3, 15)

xloc = shift_rf(xP)
yloc = shift_rf(yP)
assert xloc.valuation() == 0 and F(xloc(0)) != 0
assert yloc.valuation() >= 0
m0 = -yloc / xloc
assert m0.valuation() >= 0

MAX = 8
scalar_local = []
b_local = []
for unused_A, unused_B, a, b in pairs10:
    bloc = shift_rf(b)
    assert bloc.valuation() >= 0
    scalar_local.append(shift_rf(a) + bloc * m0)
    b_local.append(bloc)

scalar_jets = [jet(value, MAX) for value in scalar_local]
b_jets = [jet(value, MAX) for value in b_local]

# ---------------------------------------------------------------------------
# Branch/Jacobian compiler.
# ---------------------------------------------------------------------------
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


def canonical_B_space(B2):
    E = matrix(F, B2).echelon_form()
    assert E.dimensions() == (2, 16) and E.rank() == 2
    return E, tuple(int(value) for value in E.list())


def compile_space(B2):
    pairs2 = [AB_from_Brow(row) for row in B2.rows()]
    a0, b0 = lift_rf(pairs2[0][2]), lift_rf(pairs2[0][3])
    a1, b1 = lift_rf(pairs2[1][2]), lift_rf(pairs2[1][3])

    # m=(y+yP)/(x-xP) is the standard chord based at -P.
    mvalue = pencil_chord_solution(a0, b0, a1, b1, VF(V))
    radicand = chord_discriminant(xPV, -yPV, AV, mvalue)
    quartic, square_factor = squarefree_binary_quartic(radicand, UR)
    degree = int(quartic.degree())

    record = {
        "squarefree_degree": degree,
        "root_rank": None,
        "root_determinant": None,
        "euler": None,
        "D12": False,
        "factor_degrees": [],
    }
    if args.factor:
        record["factor_degrees"] = [
            [int(g.degree()), int(e)] for g, e in quartic.factor()
        ]
    if degree not in (3, 4):
        return record

    I, J = binary_quartic_invariants(quartic)
    jacA = VF(-27) * VF(I)
    jacB = VF(-27) * VF(J)
    classification = classify_finite_short_weierstrass_fibres(VR, jacA, jacB)

    root_rank = int(classification["finite_root_rank"])
    root_det = int(classification["finite_root_determinant"])
    euler = int(classification["finite_euler_number"])
    infinity_orders = tuple(
        map(int, classification["infinity_boundary"]["normalized_orders"])
    )
    infinity_kind = "smooth"
    if infinity_orders[2] > 0:
        ir, ie, idet, infinity_kind = kodaira_data_from_short_orders(*infinity_orders)
        root_rank += int(ir)
        root_det *= int(idet)
        euler += int(ie)

    record.update(
        {
            "root_rank": root_rank,
            "root_determinant": root_det,
            "euler": euler,
            "D12": (root_rank, root_det) == (12, 4),
            "infinity_orders": list(infinity_orders),
            "infinity_kind": infinity_kind,
            "finite_fibres": [
                {
                    "factor": str(item["factor"]),
                    "degree": int(item["degree"]),
                    "minimal_orders": list(map(int, item["minimal_orders"])),
                    "kodaira": item["kodaira"],
                }
                for item in classification["finite_fibres"]
            ],
        }
    )
    return record


# ---------------------------------------------------------------------------
# Candidate c(u) generation.
# ---------------------------------------------------------------------------
def trunc_rf(value, n):
    vals = jet(value, n)
    return S(vals)


def add_candidate(store, poly, label, r):
    poly = S(poly) % (u**r) if r else S.zero()
    key = tuple(int(poly[i]) for i in range(r))
    if key not in store:
        store[key] = (poly, label)


rng = random.Random(202608231409)
spaces = {}
attempts = []
found_d12 = False

for r in range(9):
    s = 8 - r
    cands = {}
    add_candidate(cands, S.zero(), "split", r)

    # Small monomial shears.
    for k in range(r):
        for lam in range(1, int(args.small_lambda) + 1):
            add_candidate(cands, F(lam) * u**k, f"+{lam}u^{k}", r)
            add_candidate(cands, -F(lam) * u**k, f"-{lam}u^{k}", r)

    # Geometry-derived truncated local series and shifts.
    if r:
        derived = {
            "m0": trunc_rf(m0, r),
            "xP": trunc_rf(xloc, r),
            "yP": trunc_rf(yloc, r),
            "A/u2": (Al // (u**2)) % (u**r),
            "B/u3": (Bl // (u**3)) % (u**r),
        }
        # Unit inverses where available.
        if xloc(0):
            derived["1/xP"] = trunc_rf(1 / xloc, r)
        if yloc(0):
            derived["1/yP"] = trunc_rf(1 / yloc, r)
        for name, poly in list(derived.items()):
            add_candidate(cands, poly, name, r)
            add_candidate(cands, -poly, "-" + name, r)
            for shift in range(1, r):
                add_candidate(cands, u**shift * poly, f"u^{shift}*{name}", r)

    # Deterministic random HNF shears.
    for index in range(int(args.random_per_split)):
        coeffs = [F(rng.randrange(int(p))) for unused in range(r)]
        poly = S(coeffs)
        add_candidate(cands, poly, f"random{index}", r)

    print(
        f"Q24SHEAR_SPLIT|r={r}|s={s}|candidates={len(cands)}|status=START",
        flush=True,
    )

    for key, (c, label) in cands.items():
        # Local HNF membership:
        #   b = u^s h
        #   scalar - c*h == 0 mod u^r.
        rows = []
        row_labels = []

        # b divisible by u^s.
        for order in range(s):
            rows.append([b_jets[column][order] for column in range(10)])
            row_labels.append(("b", order))

        # Build h=b/u^s to enough order for c*h mod u^r.
        # For each of the ten basis functions, first s coefficients must vanish
        # only after imposing the whole matrix; h coefficients use the shifted
        # formal series b_{s+j}.
        for order in range(r):
            row = []
            for column in range(10):
                value = scalar_jets[column][order]
                # coefficient of u^order in c * (b/u^s)
                correction = F(0)
                for ck in range(min(order, r - 1) + 1):
                    bj = order - ck + s
                    if bj < MAX:
                        correction += c[ck] * b_jets[column][bj]
                row.append(value - correction)
            rows.append(row)
            row_labels.append(("scalar-c*h", order))

        local = matrix(F, rows)
        rank = int(local.rank())
        if rank != 8:
            attempts.append(
                {
                    "r": r, "s": s, "label": label, "c": str(c),
                    "rank": rank, "kernel": 10-rank,
                    "space_id": None,
                }
            )
            continue

        K2 = local.right_kernel().basis_matrix()
        B2 = K2 * K10
        canonical, space_key = canonical_B_space(B2)

        if space_key not in spaces:
            compiled = compile_space(canonical)
            compiled.update(
                {
                    "space_id": len(spaces) + 1,
                    "B_echelon_rows": [
                        [int(value) for value in row] for row in canonical.rows()
                    ],
                    "modules": [],
                }
            )
            spaces[space_key] = compiled

        result = spaces[space_key]
        module_record = {
            "r": r,
            "s": s,
            "label": label,
            "c": str(c),
        }
        result["modules"].append(module_record)
        attempts.append(
            {
                **module_record,
                "rank": 8,
                "kernel": 2,
                "space_id": result["space_id"],
                "branch": result["squarefree_degree"],
                "D12": result["D12"],
            }
        )

        # Only print novel spaces; repeated c's often yield the same 2-space.
        if len(result["modules"]) == 1:
            print(
                "Q24SHEAR_SPACE|"
                f"space={result['space_id']}|r={r}|s={s}|"
                f"label={label}|c={c}|"
                f"branch={result['squarefree_degree']}|"
                f"root_rank={result['root_rank'] if result['root_rank'] is not None else 'NA'}|"
                f"root_det={result['root_determinant'] if result['root_determinant'] is not None else 'NA'}|"
                f"D12={int(result['D12'])}|"
                f"status={'D12' if result['D12'] else 'GENUS_ONE_OTHER' if result['squarefree_degree'] in (3,4) else 'NOT_GENUS_ONE'}",
                flush=True,
            )

        if result["D12"]:
            found_d12 = True
            print(
                "Q24SHEAR_HIT|"
                f"space={result['space_id']}|r={r}|s={s}|"
                f"label={label}|c={c}|branch={result['squarefree_degree']}|"
                "root_rank=12|root_det=4|status=FOUND_D12",
                flush=True,
            )
            if args.stop_on_d12:
                break

    if found_d12 and args.stop_on_d12:
        break

candidate_spaces = sorted(
    spaces.values(),
    key=lambda item: (
        0 if item["D12"] else 1,
        0 if item["squarefree_degree"] in (3,4) else 1,
        abs(item["squarefree_degree"] - 4),
        item["space_id"],
    ),
)
d12_spaces = [item for item in candidate_spaces if item["D12"]]
genus_one_spaces = [
    item for item in candidate_spaces if item["squarefree_degree"] in (3,4)
]

# Regression: split (4,4,c=0) should be the D13 self-neighbour seen already.
split44 = next(
    (
        item
        for item in candidate_spaces
        if any(
            module["r"] == 4 and module["s"] == 4 and module["label"] == "split"
            for module in item["modules"]
        )
    ),
    None,
)
if split44 is None:
    raise ArithmeticError("missing split (4,4) regression")
if split44["squarefree_degree"] != 4 or split44["root_rank"] != 13:
    raise ArithmeticError(
        "split (4,4) no longer reproduces the genus-one D13 self-neighbour"
    )

payload = {
    "schema": "elkies-k3.h3-q24-sheared-local-module-scan-modp.v1",
    "status": (
        "FOUND_D12_SHEARED_LOCAL_MODULE"
        if d12_spaces
        else "FOUND_ONLY_NON_D12_GENUS_ONE_SHEARS"
        if genus_one_spaces
        else "NO_GENUS_ONE_SHEARED_MODULE_IN_SEARCH_SET"
    ),
    "proof_boundary": (
        "Discovery scan only. A D12 hit identifies a candidate HNF local "
        "module but must still be checked on the actual resolved I9* charts "
        "against the exact q24 vertical component cycle."
    ),
    "prime": int(p),
    "inputs": {
        "q8_child": str(Q8.relative_to(ROOT)),
        "q24_section": section_source,
    },
    "collision": {"B_dimension": 16, "rank": 6, "kernel_dimension": 10},
    "search": {
        "random_per_split": int(args.random_per_split),
        "small_lambda": int(args.small_lambda),
        "attempt_count": len(attempts),
        "unique_two_spaces": len(candidate_spaces),
    },
    "split44_regression": {
        "space_id": split44["space_id"],
        "branch": split44["squarefree_degree"],
        "root_rank": split44["root_rank"],
        "root_determinant": split44["root_determinant"],
    },
    "D12_space_count": len(d12_spaces),
    "genus_one_space_count": len(genus_one_spaces),
    "candidate_spaces": candidate_spaces,
    "attempts": attempts,
}

OUT = (
    args.output.resolve()
    if args.output
    else LOCAL / f"q24-sheared-local-module-scan-mod-{p}.json"
)
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(f"OUTPUT|{OUT}", flush=True)
print(
    "Q24SHEAR_RESULT|"
    f"prime={p}|attempts={len(attempts)}|spaces={len(candidate_spaces)}|"
    f"genus_one={len(genus_one_spaces)}|D12={len(d12_spaces)}|"
    f"status={payload['status']}",
    flush=True,
)
