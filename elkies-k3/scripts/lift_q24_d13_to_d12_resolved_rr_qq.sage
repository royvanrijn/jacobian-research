#!/usr/bin/env sage -python
"""
Exact H92 q24 D13 -> D12 resolved-Riemann--Roch compiler.

This script has one job:

    exact q24 section on the exact D13 parent
        -> exact resolved RR pencil
        -> exact genus-one quartic
        -> exact minimal D12/MW5 child equation.

It does NOT discover q24, Hensel-lift q24, use modular q24 coefficients, or
search over local modules.

Inputs
------
  artifacts/local/elkies-k3/q8-q24-horizontal-section-qq.json
      status = PASS_EXACT_Q24_HORIZONTAL_SECTION

  exact corrected D13 parent:
      artifacts/local/elkies-k3/q8-corrected2cover-qq-child.json
      or generated-results fallback.

Selected resolved divisor datum
-------------------------------
The already-certified orbit-85 q24 vertical profile gives the nonzero
infinitely-near point cluster on the long I9* branch

    C01:2, C02:2, C04:2, C06:3.

This script independently reconstructs the actual I9* blow-up tree over QQ and
requires those centres to occur on one nested branch.  It then imposes those
orders in the actual strict-transform surface local rings.

Output
------
  artifacts/local/elkies-k3/q24-d13-to-d12-resolved-rr-qq.json

Terminal success
----------------
  Q24D12QQ_RESULT|...|status=PASS_EXACT_Q24_D13_TO_D12_RESOLVED_RR
"""

import argparse
import json
import time
from pathlib import Path

from sage.all import (
    QQ, ZZ, PolynomialRing, gcd, identity_matrix, lcm, matrix
)


# ---------------------------------------------------------------------------
# Paths / arguments.
# ---------------------------------------------------------------------------
def locate_repo(explicit=None):
    candidates = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd = Path.cwd().resolve()
    candidates += [cwd, *cwd.parents]
    h = Path.home()
    candidates += [
        h / "Documents" / "jacobian-research",
        h / "jacobian-research",
        h / "src" / "jacobian-research",
        h / "git" / "jacobian-research",
        h / "projects" / "jacobian-research",
    ]
    seen = set()
    for c in candidates:
        try:
            c = c.resolve()
        except Exception:
            continue
        if c in seen:
            continue
        seen.add(c)
        if (
            (c / "elkies-k3/scripts").is_dir()
            and (c / "artifacts/generated-results").is_dir()
        ):
            return c
    raise SystemExit("Could not locate jacobian-research; pass --repo PATH")


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--section", type=Path)
parser.add_argument("--parent", type=Path)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
LOCAL = ROOT / "artifacts/local/elkies-k3"
GEN = ROOT / "artifacts/generated-results"
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"

SECTION = (
    args.section.resolve()
    if args.section
    else LOCAL / "q8-q24-horizontal-section-qq.json"
)

parent_candidates = (
    [args.parent.resolve()]
    if args.parent
    else [
        LOCAL / "q8-corrected2cover-qq-child.json",
        GEN / "elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
    ]
)
PARENT = next(
    (
        p for p in parent_candidates
        if p.exists()
        and json.loads(p.read_text()).get("status")
        == "PASS_EXACT_CORRECTED_Q8_D13_CHILD"
    ),
    None,
)

OUTPUT = (
    args.output.resolve()
    if args.output
    else LOCAL / "q24-d13-to-d12-resolved-rr-qq.json"
)

if not CORE.exists():
    raise SystemExit(f"Missing compiler core: {CORE}")
if not SECTION.exists():
    raise SystemExit(f"Missing exact q24 section: {SECTION}")
if PARENT is None:
    raise SystemExit("No passing exact corrected D13 parent artifact found")

core = {}
exec(compile(CORE.read_text(), str(CORE), "exec"), core)
squarefree_binary_quartic = core["squarefree_binary_quartic"]
binary_quartic_invariants = core["binary_quartic_invariants"]
classify_finite_short_weierstrass_fibres = core[
    "classify_finite_short_weierstrass_fibres"
]
kodaira_data_from_short_orders = core["kodaira_data_from_short_orders"]

started = time.monotonic()


def stamp(stage, **kwargs):
    fields = "|".join(f"{k}={v}" for k, v in kwargs.items())
    print(
        f"Q24D12QQ|stage={stage}|elapsed={time.monotonic()-started:.3f}"
        + (f"|{fields}" if fields else ""),
        flush=True,
    )


# ---------------------------------------------------------------------------
# 1. Exact parent and exact q24 section.
# ---------------------------------------------------------------------------
parent = json.loads(PARENT.read_text())
exact = json.loads(SECTION.read_text())

if exact.get("status") != "PASS_EXACT_Q24_HORIZONTAL_SECTION":
    raise SystemExit(
        f"{SECTION} does not have status PASS_EXACT_Q24_HORIZONTAL_SECTION"
    )
if exact.get("verification", {}).get("exact_weierstrass_identity") is not True:
    raise SystemExit(f"{SECTION} does not certify the exact Weierstrass identity")

child = parent["child"]
R = PolynomialRing(QQ, "U")
U = R.gen()
K = R.fraction_field()

A = R([QQ(v) for v in child["minimal_A_coefficients_low_to_high"]])
B = R([QQ(v) for v in child["minimal_B_coefficients_low_to_high"]])

sec = exact["section"]
Z = R([QQ(v) for v in sec["Z_coefficients_low_to_high"]])
X = R([QQ(v) for v in sec["X_coefficients_low_to_high"]])
Y = R([QQ(v) for v in sec["Y_coefficients_low_to_high"]])

assert (Z.degree(), X.degree(), Y.degree()) == (24, 52, 78)
assert Z.leading_coefficient() == 1
assert gcd(X, Z) == 1
assert gcd(Y, Z) == 1

section_identity = Y**2 - X**3 - A * X * Z**4 - B * Z**6
if section_identity:
    raise ArithmeticError(
        f"exact q24 artifact fails parent identity, degree={section_identity.degree()}"
    )

xP = K(X) / K(Z**2)
yP = K(Y) / K(Z**3)
assert yP**2 == xP**3 + K(A) * xP + K(B)

stamp(
    "INPUT",
    section=SECTION.relative_to(ROOT),
    parent=PARENT.relative_to(ROOT),
    Z=24,
    X=52,
    Y=78,
    identity="PASS",
)


# ---------------------------------------------------------------------------
# Helpers for coefficient growth.
# ---------------------------------------------------------------------------
def primitive_poly(poly):
    """Scale a QQ multivariate polynomial to a primitive ZZ polynomial."""
    ring = poly.parent()
    poly = ring(poly)
    if not poly:
        return poly
    den = ZZ.one()
    for c in poly.coefficients():
        den = lcm(den, ZZ(QQ(c).denominator()))
    q = ring(poly * den)
    nums = [abs(ZZ(QQ(c).numerator())) for c in q.coefficients() if c]
    content = gcd(nums) if nums else ZZ.one()
    if content > 1:
        q = ring(q / content)
    # Deterministic sign.
    coeffs = q.coefficients()
    if coeffs and coeffs[0] < 0:
        q = -q
    return q


# ---------------------------------------------------------------------------
# 2. Exact I9* germ and exact blow-up resolution over QQ.
# ---------------------------------------------------------------------------
Delta = -16 * (4 * A**3 + 27 * B**2)
i9 = next(item for item in child["finite_fibres"] if item["kodaira"] == "I9*")
f = R(str(i9["factor"]))
if f.degree() != 1:
    raise ArithmeticError("D13 I9* factor is not rational linear")
alpha = -QQ(f[0]) / QQ(f[1])

T = PolynomialRing(QQ, "t")
t = T.gen()
Al1 = T(A(alpha + t))
Bl1 = T(B(alpha + t))
Dl1 = T(Delta(alpha + t))
local_orders = (Al1.valuation(), Bl1.valuation(), Dl1.valuation())
if local_orders != (2, 3, 15):
    raise ArithmeticError(f"unexpected D13 I9* orders {local_orders}")

S = PolynomialRing(QQ, names=("u", "x", "y"), order="degrevlex")
u, x, y = S.gens()
Al = S(A(alpha + u))
Bl = S(B(alpha + u))
surface0 = primitive_poly(y**2 - x**3 - Al * x - Bl)
assert surface0(0, 0, 0) == 0


def hessian_matrix(poly, point):
    return matrix(
        QQ,
        [[poly.derivative(a, b)(*point) for b in (u, x, y)] for a in (u, x, y)],
    )


def shifted_polynomial(poly, point):
    a, b, c = map(QQ, point)
    return S(poly(a + u, b + x, c + y))


def order_at_point(poly, point):
    shifted = shifted_polynomial(poly, point)
    if not shifted:
        return 10**9
    return min(sum(exp) for exp, coef in shifted.dict().items() if coef)


def tangent_component_count(poly, point, multiplicity):
    if multiplicity != 2:
        return 1, None
    rank = int(hessian_matrix(poly, point).rank())
    if rank not in (1, 2, 3):
        raise ArithmeticError(f"unexpected tangent quadratic rank {rank}")
    return (2 if rank == 2 else 1), rank


def divide_power(poly, exceptional, power):
    q = S(poly)
    for _ in range(int(power)):
        q, rem = q.quo_rem(exceptional)
        if rem:
            raise ArithmeticError("required exceptional power does not divide")
    return primitive_poly(q)


def canonical_projective(direction):
    vals = [QQ(v) for v in direction]
    idx = next((i for i, v in enumerate(vals) if v), None)
    if idx is None:
        raise ArithmeticError("zero projective tangent direction")
    inv = vals[idx] ** -1
    return tuple(QQ(v * inv) for v in vals)


def singular_points_on_exceptional(poly, exceptional):
    ideal = S.ideal(
        [
            exceptional,
            poly,
            poly.derivative(u),
            poly.derivative(x),
            poly.derivative(y),
        ]
    )
    if ideal.is_one():
        return []
    if ideal.dimension() != 0:
        raise ArithmeticError(
            "exceptional singular locus is not zero-dimensional: "
            + repr(ideal.groebner_basis())
        )
    try:
        variety = ideal.variety(ring=QQ)
    except TypeError:
        variety = ideal.variety()
    points = []
    for sol in variety:
        point = tuple(QQ(sol.get(v, QQ.zero())) for v in (u, x, y))
        if point not in points:
            points.append(point)
    points.sort()
    for pt in points:
        assert exceptional(*pt) == 0
        assert poly(*pt) == 0
        assert all(poly.derivative(v)(*pt) == 0 for v in (u, x, y))
    return points


def chart_substitutions(point, kind):
    a, b, c = map(QQ, point)
    if kind == "u":
        return (a + u, b + u * x, c + u * y), u
    if kind == "x":
        return (a + x * u, b + x, c + x * y), x
    if kind == "y":
        return (a + y * u, b + y * x, c + y), y
    raise ValueError(kind)


def direction_from_chart(kind, point):
    pu, px, py = map(QQ, point)
    if kind == "u":
        assert pu == 0
        return canonical_projective((1, px, py))
    if kind == "x":
        assert px == 0
        return canonical_projective((pu, 1, py))
    assert kind == "y" and py == 0
    return canonical_projective((pu, px, 1))


center_records = []
leaf_count = 0
next_center = 0
MAX_CENTERS = 20


def blow_center(poly, point, path, depth):
    global next_center, leaf_count
    if len(center_records) >= MAX_CENTERS:
        raise RuntimeError("exact I9* resolution exceeded max centres")

    multiplicity = order_at_point(poly, point)
    if multiplicity < 2:
        raise ArithmeticError(("attempted blow-up at smooth point", path, point))

    components, tangent_rank = tangent_component_count(poly, point, multiplicity)
    label = f"C{next_center+1:02d}"
    next_center += 1

    record = {
        "label": label,
        "path": path,
        "depth": int(depth),
        "point": tuple(map(QQ, point)),
        "multiplicity": int(multiplicity),
        "geometric_exceptional_components": int(components),
        "tangent_rank": tangent_rank,
        "children": [],
    }
    center_records.append(record)

    candidates = {}
    for kind in ("u", "x", "y"):
        subs, e = chart_substitutions(point, kind)
        transformed = S(poly(*subs))
        strict = divide_power(transformed, e, multiplicity)
        for q in singular_points_on_exceptional(strict, e):
            direction = direction_from_chart(kind, q)
            candidates.setdefault(direction, []).append((kind, q, strict))

    preference = {"u": 0, "x": 1, "y": 2}
    for direction in sorted(candidates):
        reps = candidates[direction]
        reps.sort(key=lambda item: preference[item[0]])
        kind, q, strict = reps[0]
        direction_text = ",".join(str(v) for v in direction)
        child_path = f"{path}/{label}:{kind}:{direction_text}"
        record["children"].append(
            {
                "direction": direction,
                "selected_chart": kind,
                "path": child_path,
            }
        )
        blow_center(strict, q, child_path, depth + 1)

    if not candidates:
        leaf_count += 1


blow_center(surface0, (QQ.zero(), QQ.zero(), QQ.zero()), "root", 0)

center_count = len(center_records)
geometric_components = sum(
    int(r["geometric_exceptional_components"]) for r in center_records
)
split_centers = [
    r["label"]
    for r in center_records
    if r["geometric_exceptional_components"] == 2
]

if center_count != 12 or geometric_components != 13:
    raise ArithmeticError(
        f"exact I9* resolution gave centres/components "
        f"{center_count}/{geometric_components}, expected 12/13"
    )

records = {r["label"]: r for r in center_records}

# C01,C02,C04,C06 are labels in the modular DFS chronology, not intrinsic
# geometric names. Over QQ sibling branches can be enumerated differently.
# Intrinsically these are the first four centres on the long D13 arm.
parent_of = {}
for parent in center_records:
    for child_info in parent["children"]:
        target_path = str(child_info["path"])
        matches = [
            rec["label"] for rec in center_records
            if str(rec["path"]) == target_path
        ]
        if len(matches) != 1:
            raise ArithmeticError(
                f"could not identify unique child for {parent['label']} path={target_path}: {matches}"
            )
        child_label = matches[0]
        if child_label in parent_of and parent_of[child_label] != parent["label"]:
            raise ArithmeticError(f"centre {child_label} has multiple parents")
        parent_of[child_label] = parent["label"]

max_depth = max(int(r["depth"]) for r in center_records)
deepest = [r for r in center_records if int(r["depth"]) == max_depth]
if not deepest:
    raise ArithmeticError("exact I9* resolution has no deepest centre")

def ancestor_chain(label):
    chain = []
    cur = str(label)
    while True:
        chain.append(cur)
        if cur not in parent_of:
            break
        cur = parent_of[cur]
    chain.reverse()
    return chain

deep_chains = [ancestor_chain(r["label"]) for r in deepest]
if any(len(chain) < 4 for chain in deep_chains):
    raise ArithmeticError(f"D13 long arm is unexpectedly short: {deep_chains}")

prefix4 = tuple(deep_chains[0][:4])
if any(tuple(chain[:4]) != prefix4 for chain in deep_chains[1:]):
    raise ArithmeticError(
        f"deep D13 branches disagree before selected q24 cluster: {deep_chains}"
    )

SELECTED_CLUSTER = tuple(zip(prefix4, (2, 2, 2, 3)))
MODULAR_CLUSTER_LABELS = ("C01", "C02", "C04", "C06")

selected_depths = tuple(int(records[name]["depth"]) for name, _ in SELECTED_CLUSTER)
if selected_depths != (0, 1, 2, 3):
    raise ArithmeticError(
        f"selected long-arm depths are {selected_depths}, expected (0,1,2,3)"
    )

print(
    "Q24D12QQ_CLUSTER_MAP|"
    + "|".join(
        f"{mod}->{exact}:{order}"
        for mod, (exact, order)
        in zip(MODULAR_CLUSTER_LABELS, SELECTED_CLUSTER)
    )
    + "|status=PASS_INTRINSIC_LONG_ARM",
    flush=True,
)

stamp(
    "I9_RESOLUTION",
    base=alpha,
    orders="2,3,15",
    centers=center_count,
    components=geometric_components,
    split=",".join(split_centers),
    cluster=",".join(f"{name}:{order}" for name, order in SELECTED_CLUSTER),
    status="PASS",
)


# ---------------------------------------------------------------------------
# 3. Exact smooth-collision RR reduction: 16 -> 10.
#
# B has degree <=15.  A is determined mod Z^2 by
#     A X = B Y (mod Z^2).
# Requiring deg(A)<=41 gives six high-degree conditions.
# ---------------------------------------------------------------------------
modulus = Z**2
stamp("COLLISION_INVERSE", status="BEGIN")
Xinv = X.inverse_mod(modulus)
stamp("COLLISION_INVERSE", status="PASS")

A_columns = [((U**i) * Y * Xinv) % modulus for i in range(16)]
H = matrix(
    QQ,
    6,
    16,
    lambda row, col: A_columns[col][42 + row],
)
if H.rank() != 6:
    raise ArithmeticError(f"exact reduced collision rank {H.rank()}, expected 6")
K10 = H.right_kernel().basis_matrix()
if K10.dimensions() != (10, 16):
    raise ArithmeticError(
        f"exact reduced collision kernel {K10.dimensions()}, expected (10,16)"
    )


def AB_from_Brow(row):
    Bcoef = R(list(row))
    Acoef = (Bcoef * Y * Xinv) % modulus
    assert Acoef.degree() <= 41
    assert (Acoef * X - Bcoef * Y) % modulus == 0
    return Acoef, Bcoef


post_pairs = [AB_from_Brow(row) for row in K10.rows()]
stamp(
    "SMOOTH_COLLISION",
    ambient=16,
    rank=6,
    post_collision=10,
    equivalent_full_ambient=58,
    status="PASS",
)


# ---------------------------------------------------------------------------
# 4. Build exact common-unit local numerators at I9*.
# ---------------------------------------------------------------------------
Zl = S(Z(alpha + u))
Xl = S(X(alpha + u))
Yl = S(Y(alpha + u))
if Zl(0, 0, 0) == 0 or Xl(0, 0, 0) == 0:
    raise ArithmeticError("q24 section is not a unit in the selected I9* chart")

# m=(y+yP)/(x-xP)
#   m_num = y Z^3 + Y
#   m_den = Z (x Z^2 - X)
m_num = S(y * Zl**3 + Yl)
m_den = S(Zl * (x * Zl**2 - Xl))
common_den = S(Zl**2 * m_den)
if common_den(0, 0, 0) == 0:
    raise ArithmeticError("marked-chord common denominator is not a local unit")

local_numerators = []
for AA, BB in post_pairs:
    AAl = S(AA(alpha + u))
    BBl = S(BB(alpha + u))
    numerator = S(AAl * m_den + BBl * Zl * m_num)
    local_numerators.append(primitive_poly(numerator))

if len(local_numerators) != 10:
    raise AssertionError("post-collision local basis is not 10-dimensional")

stamp("LOCAL_NUMERATORS", basis=10, common_den_unit=1, status="PASS")


# ---------------------------------------------------------------------------
# 5. Actual local-surface quotient conditions along the resolved branch.
# ---------------------------------------------------------------------------
def order_at_origin(poly):
    poly = S(poly)
    if not poly:
        return 10**9
    return min(sum(exp) for exp, coef in poly.dict().items() if coef)


def shifted(poly, point):
    a, b, c = map(QQ, point)
    return S(poly(a + u, b + x, c + y))


def monomials_exact(total):
    result = []
    for i in range(total + 1):
        for j in range(total - i + 1):
            k = total - i - j
            result.append(u**i * x**j * y**k)
    return result


def monomials_below(total):
    return [
        mon
        for degree in range(total)
        for mon in monomials_exact(degree)
    ]


def local_order_matrix(basis, surface_eq, point, required_order):
    required_order = int(required_order)
    ss = shifted(surface_eq, point)
    ideal = S.ideal([ss] + monomials_exact(required_order))
    gb = ideal.groebner_basis()
    mons = monomials_below(required_order)
    remainders = []
    for poly in basis:
        rem = shifted(poly, point).reduce(gb)
        if rem and max(sum(exp) for exp in rem.dict()) >= required_order:
            raise ArithmeticError(
                "truncated local remainder escaped maximal-ideal quotient"
            )
        remainders.append(rem)
    M = matrix(
        QQ,
        len(mons),
        len(basis),
        lambda row, col: remainders[col].monomial_coefficient(mons[row]),
    )
    return M, mons


def canonical_after_condition(poly, surface_eq, point, required_order):
    a, b, c = map(QQ, point)
    sp = shifted(poly, point)
    ss = shifted(surface_eq, point)
    rem = sp.reduce(S.ideal([ss]).groebner_basis())
    if order_at_origin(rem) < required_order:
        raise ArithmeticError(
            f"surface-normal representative has order {order_at_origin(rem)}, "
            f"expected >= {required_order}"
        )
    unshifted = S(rem(u - a, x - b, y - c))
    return primitive_poly(unshifted)


def child_chart(parent, child):
    target = str(child["path"])
    matches = [
        str(item["selected_chart"])
        for item in parent["children"]
        if str(item["path"]) == target
    ]
    if len(matches) != 1:
        raise ArithmeticError(
            f"could not identify unique exact chart "
            f"{parent['label']} -> {child['label']}: {matches}"
        )
    return matches[0]


current_basis = list(local_numerators)
current_surface = surface0
transform = identity_matrix(QQ, 10)
condition_ledger = []

for step_index, (center_name, required_order) in enumerate(SELECTED_CLUSTER):
    record = records[center_name]
    point = record["point"]

    M, quotient_monomials = local_order_matrix(
        current_basis, current_surface, point, required_order
    )
    local_rank = int(M.rank())
    kernel = M.right_kernel().basis_matrix()
    before = len(current_basis)
    after = int(kernel.nrows())

    if before - local_rank != after:
        raise ArithmeticError("local rank/nullity mismatch")

    new_basis = []
    for row in kernel.rows():
        combination = S(
            sum(row[j] * current_basis[j] for j in range(before))
        )
        new_basis.append(
            canonical_after_condition(
                combination, current_surface, point, required_order
            )
        )

    transform = kernel * transform
    current_basis = new_basis
    cumulative_codim = 10 - after

    condition_ledger.append(
        {
            "center": center_name,
            "additional_order": int(required_order),
            "quotient_rows": int(M.nrows()),
            "local_rank": local_rank,
            "dimension_before": before,
            "dimension_after": after,
            "cumulative_codimension": cumulative_codim,
            "quotient_monomials": [str(mon) for mon in quotient_monomials],
        }
    )

    print(
        "Q24D12QQ_CENTER|"
        f"center={center_name}|order={required_order}|before={before}|"
        f"rank={local_rank}|after={after}|codim={cumulative_codim}|status=PASS",
        flush=True,
    )

    if step_index + 1 < len(SELECTED_CLUSTER):
        next_name = SELECTED_CLUSTER[step_index + 1][0]
        next_record = records[next_name]
        kind = child_chart(record, next_record)
        subs, exceptional = chart_substitutions(point, kind)

        current_surface = divide_power(
            S(current_surface(*subs)),
            exceptional,
            int(record["multiplicity"]),
        )

        current_basis = [
            divide_power(S(poly(*subs)), exceptional, required_order)
            for poly in current_basis
        ]

        print(
            "Q24D12QQ_BLOWUP|"
            f"from={center_name}|to={next_name}|chart={kind}|"
            f"surface_mult={record['multiplicity']}|"
            f"section_divide={required_order}|status=PASS",
            flush=True,
        )

final_dimension = len(current_basis)
resolved_codim = 10 - final_dimension
if final_dimension != 2 or resolved_codim != 8:
    raise ArithmeticError(
        f"resolved RR dimension mismatch: codim={resolved_codim}, "
        f"h0={final_dimension}; expected 8/2"
    )

B2 = transform * K10
if B2.dimensions() != (2, 16) or B2.rank() != 2:
    raise ArithmeticError(
        f"final exact B-space has dimensions/rank "
        f"{B2.dimensions()}/{B2.rank()}, expected (2,16)/2"
    )

final_pairs = [AB_from_Brow(row) for row in B2.rows()]

stamp(
    "RESOLVED_RR",
    post_collision=10,
    resolved_codim=8,
    h0=2,
    status="PASS",
)


# ---------------------------------------------------------------------------
# 6. Compile exact degree-two pencil -> genus-one quartic.
# ---------------------------------------------------------------------------
VR = PolynomialRing(QQ, "V")
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


def pair_rf(pair):
    AA, BB = pair
    return K(AA) / K(Z**2), K(BB) / K(Z)


a0_raw, b0_raw = pair_rf(final_pairs[0])
a1_raw, b1_raw = pair_rf(final_pairs[1])

a0, b0 = lift_rf(a0_raw), lift_rf(b0_raw)
a1, b1 = lift_rf(a1_raw), lift_rf(b1_raw)
xPV, yPV = lift_rf(xP), lift_rf(yP)
AV, BV = lift_poly(A), lift_poly(B)

den = b1 - VF(V) * b0
if not den:
    raise ArithmeticError("exact resolved q24 pencil has degenerate chord direction")

# New base V = f1/f0, so f1 - V*f0 = 0.
mval = -(a1 - VF(V) * a0) / den

XR = PolynomialRing(UK, "x")
xx = XR.gen()
yline = XR(mval) * (xx - XR(xPV)) - XR(yPV)
relation = yline**2 - xx**3 - XR(AV) * xx - XR(BV)
quadratic, remainder = relation.quo_rem(xx - XR(xPV))
if remainder or quadratic.degree() != 2:
    raise ArithmeticError("exact q24 chord elimination did not leave a quadratic")

disc = UK(quadratic[1] ** 2 - 4 * quadratic[2] * quadratic[0])
quartic, square_factor = squarefree_binary_quartic(disc, UR)
quartic_degree = int(quartic.degree())
if quartic_degree not in (3, 4):
    raise ArithmeticError(
        f"resolved exact q24 pencil is not genus one: "
        f"squarefree degree={quartic_degree}"
    )

stamp("QUARTIC", degree=quartic_degree, genus=1, status="PASS")


# ---------------------------------------------------------------------------
# 7. Binary-quartic invariants -> exact minimal D12 child.
# ---------------------------------------------------------------------------
I, J = binary_quartic_invariants(quartic)
jacA = VF(-27) * VF(I)
jacB = VF(-27) * VF(J)

classification = classify_finite_short_weierstrass_fibres(VR, jacA, jacB)

finite_data = [
    {
        "factor": str(item["factor"]),
        "degree": int(item["degree"]),
        "raw_orders": list(map(int, item["raw_orders"])),
        "scaling": int(item["scaling"]),
        "minimal_orders": list(map(int, item["minimal_orders"])),
        "kodaira": item["kodaira"],
    }
    for item in classification["finite_fibres"]
]

root_rank = int(classification["finite_root_rank"])
root_euler = int(classification["finite_euler_number"])
root_det = int(classification["finite_root_determinant"])

infinity = classification["infinity_boundary"]
infinity_orders = tuple(map(int, infinity["normalized_orders"]))
infinity_kind = "smooth"

if infinity_orders[2] > 0:
    irank, ieuler, idet, infinity_kind = kodaira_data_from_short_orders(
        *infinity_orders
    )
    root_rank += int(irank)
    root_euler += int(ieuler)
    root_det *= int(idet)

minimal_A = classification["finite_minimization"]["minimal_a"]
minimal_B = classification["finite_minimization"]["minimal_b"]
minimal_Delta = classification["finite_minimization"]["minimal_discriminant"]

if minimal_A.degree() > 8 or minimal_B.degree() > 12 or minimal_Delta.degree() > 24:
    raise ArithmeticError(
        "minimal q24 child violates elliptic K3 degree bounds: "
        f"{minimal_A.degree()}/{minimal_B.degree()}/{minimal_Delta.degree()}"
    )

if minimal_Delta != -16 * (4 * minimal_A**3 + 27 * minimal_B**2):
    raise ArithmeticError("minimal child discriminant does not replay exactly")

if root_euler != 24:
    raise ArithmeticError(f"child Euler number {root_euler}, expected 24")
if (root_rank, root_det) != (12, 4):
    raise ArithmeticError(
        f"child root data {(root_rank, root_det)}, expected D12=(12,4)"
    )

mw_rank = 19 - 2 - root_rank
if mw_rank != 5:
    raise ArithmeticError(f"child MW rank {mw_rank}, expected 5")

stamp(
    "CHILD",
    root="D12",
    root_rank=root_rank,
    root_det=root_det,
    euler=root_euler,
    MW=mw_rank,
    infinity=f"{infinity_orders},{infinity_kind}",
    status="PASS",
)


# ---------------------------------------------------------------------------
# 8. Write one production artifact.
# ---------------------------------------------------------------------------
def qlist(poly):
    return [str(v) for v in poly.list()]


payload = {
    "schema": "elkies-k3.h92-q24-d13-to-d12-resolved-rr-qq.v1",
    "status": "PASS_EXACT_Q24_D13_TO_D12_RESOLVED_RR",
    "inputs": {
        "q24_section": str(SECTION.relative_to(ROOT)),
        "d13_parent": str(PARENT.relative_to(ROOT)),
    },
    "q24_section": {
        "Z_degree": 24,
        "X_degree": 52,
        "Y_degree": 78,
        "exact_parent_identity": True,
    },
    "old_I9star_resolution": {
        "base": str(alpha),
        "local_orders_A_B_Delta": list(map(int, local_orders)),
        "blowup_centers": center_count,
        "geometric_exceptional_components": geometric_components,
        "split_centers": split_centers,
    },
    "resolved_rr": {
        "generic_fibre": "L(O+P)",
        "marked_chord": "m=(y+yP)/(x-xP)",
        "smooth_collision": {
            "reduced_B_dimension": 16,
            "condition_rank": 6,
            "post_collision_dimension": 10,
            "equivalent_full_ambient_dimension": 58,
        },
        "cluster": [
            {
                "canonical_modular_center": modular_name,
                "exact_QQ_center": name,
                "additional_order": int(order),
                "depth": int(records[name]["depth"]),
            }
            for modular_name, (name, order)
            in zip(MODULAR_CLUSTER_LABELS, SELECTED_CLUSTER)
        ],
        "condition_ledger": condition_ledger,
        "resolved_codimension": 8,
        "kernel_dimension": 2,
        "h0": 2,
        "kernel_B_coefficients_low_to_high": [
            [str(v) for v in row] for row in B2.rows()
        ],
        "pencil": [
            {
                "A_coefficients_low_to_high": qlist(AA),
                "B_coefficients_low_to_high": qlist(BB),
            }
            for AA, BB in final_pairs
        ],
        "new_base": "V=f1/f0",
    },
    "quartic": {
        "squarefree_degree": quartic_degree,
        "coefficients_in_U_low_to_high": [str(v) for v in quartic.list()],
        "binary_quartic_I": str(I),
        "binary_quartic_J": str(J),
    },
    "jacobian_raw": {
        "A": str(jacA),
        "B": str(jacB),
        "equation": "y^2=x^3+A(V)*x+B(V)",
    },
    "child": {
        "minimal_A_coefficients_low_to_high": qlist(minimal_A),
        "minimal_B_coefficients_low_to_high": qlist(minimal_B),
        "minimal_Delta_coefficients_low_to_high": qlist(minimal_Delta),
        "finite_fibres": finite_data,
        "infinity": {
            "minimal_orders": list(infinity_orders),
            "kodaira": infinity_kind,
        },
        "root_lattice": "D12",
        "root_rank": root_rank,
        "root_determinant": root_det,
        "root_euler": root_euler,
        "MW_rank_if_rho_19": mw_rank,
        "equation": "y^2=x^3+A_D12(V)*x+B_D12(V)",
    },
    "verification": {
        "exact_q24_parent_identity": True,
        "exact_I9star_resolution": True,
        "resolved_h0_two": True,
        "genus_one_quartic": True,
        "exact_minimal_discriminant_identity": True,
        "euler_24": True,
        "D12_root_data": True,
        "MW_rank_5_if_rho_19": True,
    },
    "boundary": (
        "This certifies the characteristic-zero equation arrow for the selected "
        "q24 orbit-85 resolved cluster: D13/MW4 -> D12/MW5.  The downstream "
        "choice of D12 zero/MW basis and the q6 D12->A11 resolved pencil are "
        "separate certificates."
    ),
}

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(f"OUTPUT|{OUTPUT}", flush=True)
print(
    "Q24D12QQ_RESULT|"
    f"h0=2|quartic={quartic_degree}|root_rank={root_rank}|"
    f"root_det={root_det}|euler={root_euler}|MW={mw_rank}|"
    "status=PASS_EXACT_Q24_D13_TO_D12_RESOLVED_RR",
    flush=True,
)
