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
from math import factorial
from pathlib import Path

from sage.all import (
    GF, QQ, ZZ, PolynomialRing, gcd, identity_matrix, lcm, matrix
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
pencil_chord_solution = core["pencil_chord_solution"]
chord_discriminant = core["chord_discriminant"]

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

# Exact QQ DFS labels are not intrinsic.  Match the exact resolution tree
# to the certified GF(100003) tree by reducing every projective tangent
# direction.  GF(100003) is used ONLY for geometric identification.
GEOM_PRIME = ZZ(100003)
MOD_RESOLUTION_PATH = LOCAL / "q24-i9star-resolution-mod-100003.json"
MOD_CLUSTER_PATH = LOCAL / "q24-i9star-effective-cluster-mod-100003.json"
MOD_RR_PATH = LOCAL / "q24-d12-resolved-cluster-rr-mod-100003.json"

for geom_path in (MOD_RESOLUTION_PATH, MOD_CLUSTER_PATH):
    if not geom_path.exists():
        raise SystemExit(
            f"Missing geometric prerequisite {geom_path}; "
            "cannot identify exact infinitely-near centres safely"
        )

mod_resolution = json.loads(MOD_RESOLUTION_PATH.read_text())
mod_cluster = json.loads(MOD_CLUSTER_PATH.read_text())

assert mod_resolution["status"] == "PASS_EXPLICIT_MODP_I9STAR_D13_COMPONENT_RESOLUTION"
assert mod_cluster["status"] in (
    "PASS_H3_Q24_EFFECTIVE_I9STAR_CLUSTER",
    "CANDIDATE_H3_Q24_EFFECTIVE_I9STAR_CLUSTER",
)

Fp_geom = GF(GEOM_PRIME)

def _red_geom(q):
    q = QQ(q)
    d = ZZ(q.denominator())
    if d % GEOM_PRIME == 0:
        raise ZeroDivisionError(
            f"geometry denominator divisible by {GEOM_PRIME}: {q}"
        )
    return Fp_geom(ZZ(q.numerator())) / Fp_geom(d)

def _canon_mod_direction(direction):
    vals = [_red_geom(v) for v in direction]
    pivot = next((i for i,v in enumerate(vals) if v), None)
    if pivot is None:
        raise ArithmeticError("zero projective direction")
    scale = vals[pivot]**-1
    return tuple(int(v*scale) for v in vals)

exact_by_path = {str(r["path"]): r for r in center_records}
mod_centers = list(mod_resolution["centers"])
mod_by_path = {str(r["path"]): r for r in mod_centers}

exact_root = exact_by_path["root"]
mod_root = mod_by_path["root"]

exact_to_mod = {str(exact_root["label"]): str(mod_root["label"])}
queue = [(exact_root, mod_root)]

while queue:
    ep, mp = queue.pop(0)

    mod_edges = {}
    for edge in mp.get("children", []):
        kind = str(edge["selected_chart"])
        direction = tuple(map(int, edge["direction"]))
        child_path = (
            f"{mp['path']}/{mp['label']}:{kind}:"
            + ",".join(map(str, direction))
        )
        if child_path not in mod_by_path:
            raise ArithmeticError(f"missing modular child path {child_path}")
        mod_edges[(kind, direction)] = mod_by_path[child_path]

    for edge in ep.get("children", []):
        kind = str(edge["selected_chart"])
        direction = _canon_mod_direction(edge["direction"])
        key = (kind, direction)
        if key not in mod_edges:
            raise ArithmeticError(
                f"no modular counterpart for exact edge "
                f"{ep['label']}:{kind}:{direction}; available={sorted(mod_edges)}"
            )
        exact_child_path = str(edge["path"])
        if exact_child_path not in exact_by_path:
            raise ArithmeticError(f"missing exact child path {exact_child_path}")
        ec = exact_by_path[exact_child_path]
        mc = mod_edges[key]
        exact_to_mod[str(ec["label"])] = str(mc["label"])
        queue.append((ec, mc))

if len(exact_to_mod) != len(center_records):
    raise ArithmeticError(
        f"incomplete exact/modular centre map "
        f"{len(exact_to_mod)}/{len(center_records)}"
    )

mod_to_exact = {}
for exact_label, mod_label in exact_to_mod.items():
    if mod_label in mod_to_exact and mod_to_exact[mod_label] != exact_label:
        raise ArithmeticError(f"noninjective modular centre map at {mod_label}")
    mod_to_exact[mod_label] = exact_label

mod_plan = [
    (str(row["center"]), int(row["additional_point_order"]))
    for row in mod_cluster["common_nonzero_centre_plan"]
]
assert mod_plan == [("C01",2),("C02",2),("C04",2),("C06",3)]

SELECTED_CLUSTER = tuple(
    (mod_to_exact[mod_label], order)
    for mod_label, order in mod_plan
)

print(
    "Q24D12QQ_CLUSTER_MAP|"
    + "|".join(
        f"{mod}->{exact}:{order}"
        for (mod,order),(exact,unused)
        in zip(mod_plan, SELECTED_CLUSTER)
    )
    + "|method=REDUCE_PROJECTIVE_DIRECTIONS_MOD_100003"
    + "|status=PASS_CERTIFIED_GEOMETRY_MAP",
    flush=True,
)

# Optional modular resolved-RR ledger.  If present, every exact local rank
# must replay the modular one.
expected_mod_ledger = {}
if MOD_RR_PATH.exists():
    try:
        mod_rr = json.loads(MOD_RR_PATH.read_text())
        expected_mod_ledger = {
            str(row["center"]): row
            for row in mod_rr.get("resolved_cluster", {}).get(
                "condition_ledger", []
            )
        }
        print(
            "Q24D12QQ_MODULAR_SHADOW|"
            f"centers={len(expected_mod_ledger)}|"
            f"artifact={MOD_RR_PATH.relative_to(ROOT)}|status=LOADED",
            flush=True,
        )
    except Exception as exc:
        print(
            "Q24D12QQ_MODULAR_SHADOW|"
            f"status=IGNORED|reason={type(exc).__name__}:{exc}",
            flush=True,
        )
        expected_mod_ledger = {}

stamp(
    "I9_RESOLUTION",
    base=alpha,
    orders="2,3,15",
    centers=center_count,
    components=geometric_components,
    split=",".join(split_centers),
    cluster=",".join(
        f"{exact_to_mod[name]}->{name}:{order}"
        for name,order in SELECTED_CLUSTER
    ),
    status="PASS",
)

def primitive_integer_basis(M):
    rows = []
    for row in M.rows():
        den = ZZ.one()
        for value in row:
            den = lcm(den, ZZ(QQ(value).denominator()))
        ints = [ZZ(QQ(value) * den) for value in row]
        content = ZZ.zero()
        for value in ints:
            content = gcd(content, abs(value))
        if content > 1:
            ints = [value // content for value in ints]
        pivot = next((value for value in ints if value), ZZ.zero())
        if pivot < 0:
            ints = [-value for value in ints]
        rows.append([QQ(value) for value in ints])
    return matrix(QQ, rows)


# ---------------------------------------------------------------------------
# 3. Exact smooth-collision RR reduction: 16 -> 10.
#
# B has degree <=15.  A is determined mod Z^2 by
#     A X = B Y (mod Z^2).
# Requiring deg(A)<=41 gives six high-degree conditions.
# ---------------------------------------------------------------------------
modulus = Z**2
stamp("COLLISION_INVERSE", status="BEGIN")
INV_CACHE = LOCAL / "q24-xinv-mod-z2-qq.json"
Xinv = None

if INV_CACHE.exists():
    try:
        cached = json.loads(INV_CACHE.read_text())
        candidate = R([QQ(v) for v in cached["coefficients_low_to_high"]])
        if (X * candidate) % modulus == R.one():
            Xinv = candidate
            stamp(
                "COLLISION_INVERSE",
                method="CACHE",
                degree=Xinv.degree(),
                status="PASS",
            )
    except Exception:
        Xinv = None

if Xinv is None:
    # modulus = Z^2.  First invert modulo degree-24 Z, then Newton-lift:
    #
    #   q2 = q * (2-X*q)
    #
    # Since 1-X*q = 0 mod Z, its new error is its square, hence 0 mod Z^2.
    stamp("COLLISION_INVERSE_MOD_Z", degree=Z.degree(), status="BEGIN")
    invZ = X.inverse_mod(Z)
    stamp("COLLISION_INVERSE_MOD_Z", degree=Z.degree(), status="PASS")

    stamp("COLLISION_INVERSE_LIFT_Z2", status="BEGIN")
    Xinv = (invZ * (2 - X * invZ)) % modulus
    if (X * Xinv) % modulus != R.one():
        raise ArithmeticError("Newton lift did not produce X^-1 mod Z^2")
    stamp(
        "COLLISION_INVERSE_LIFT_Z2",
        degree=Xinv.degree(),
        status="PASS",
    )

    INV_CACHE.write_text(json.dumps({
        "schema": "elkies-k3.q24-xinv-mod-z2-qq.v1",
        "status": "PASS_EXACT_XINV_MOD_Z2",
        "modulus": "Z^2",
        "Z_degree": int(Z.degree()),
        "X_degree": int(X.degree()),
        "inverse_degree": int(Xinv.degree()),
        "coefficients_low_to_high": [str(v) for v in Xinv.list()],
    }, indent=2, sort_keys=True) + "\n")

    stamp(
        "COLLISION_INVERSE",
        method="MOD_Z_NEWTON_LIFT",
        cache=INV_CACHE.relative_to(ROOT),
        status="PASS",
    )

# Compute Y/X modulo Z^2 ONCE.  The previous implementation repeated the
# huge exact Y*Xinv product sixteen times.
stamp("COLLISION_COLUMNS", method="ONE_PRODUCT_PLUS_SHIFTS", status="BEGIN")
YXinv = (Y * Xinv) % modulus

def multiply_by_U_mod_monic(poly, monic_modulus):
    """Fast U*poly mod a monic modulus; input already has degree < deg(modulus)."""
    q = U * poly
    d = monic_modulus.degree()
    if q.degree() >= d:
        # U*poly has degree at most d, hence exactly one reduction is enough.
        q -= q[d] * monic_modulus
    assert q.degree() < d
    return R(q)

A_columns = [YXinv]
for _ in range(1, 16):
    A_columns.append(multiply_by_U_mod_monic(A_columns[-1], modulus))

stamp("COLLISION_COLUMNS", columns=16, status="PASS")

H = matrix(
    QQ,
    6,
    16,
    lambda row, col: A_columns[col][42 + row],
)
# Do NOT ask Sage for a generic QQ right-kernel here: the 6x16 entries
# contain enormous rationals and generic elimination creates massive
# intermediate coefficients.
#
# Use a good finite prime only to locate an invertible 6-column minor.
# The actual solve and kernel are then constructed and verified over QQ.
stamp("COLLISION_KERNEL", method="MODULAR_PIVOTS_EXACT_6X6", status="BEGIN")

pivot_prime = ZZ(100003)
Fp = GF(pivot_prime)

def reduce_q_mod_p(value):
    value = QQ(value)
    den = ZZ(value.denominator())
    if den % pivot_prime == 0:
        raise ZeroDivisionError(
            f"collision pivot denominator divisible by {pivot_prime}"
        )
    return Fp(ZZ(value.numerator())) / Fp(den)

Hmod = matrix(
    Fp, H.nrows(), H.ncols(),
    [reduce_q_mod_p(v) for v in H.list()]
)

if Hmod.rank() != 6:
    raise ArithmeticError(
        f"collision matrix rank mod {pivot_prime} is {Hmod.rank()}, expected 6"
    )

pivots = tuple(int(i) for i in Hmod.pivots())
if len(pivots) != 6:
    raise ArithmeticError(f"expected 6 collision pivots, got {pivots}")

free = tuple(i for i in range(16) if i not in pivots)
if len(free) != 10:
    raise ArithmeticError(f"expected 10 collision free columns, got {free}")

stamp(
    "COLLISION_KERNEL_PIVOTS",
    prime=pivot_prime,
    pivots=",".join(map(str, pivots)),
    free=",".join(map(str, free)),
    status="PASS",
)

P = H.matrix_from_columns(pivots)
Ffree = H.matrix_from_columns(free)

stamp("COLLISION_KERNEL_SOLVE", size="6x6_by_10", status="BEGIN")
pivot_values = P.solve_right(-Ffree)   # 6 x 10, exact QQ
stamp("COLLISION_KERNEL_SOLVE", size="6x6_by_10", status="PASS")

rows = []
for j in range(10):
    row = [QQ.zero()] * 16
    row[free[j]] = QQ.one()
    for i, col in enumerate(pivots):
        row[col] = pivot_values[i, j]
    rows.append(row)

K10 = matrix(QQ, rows)
stamp("COLLISION_KERNEL_NORMALIZE", status="BEGIN")
K10 = primitive_integer_basis(K10)
stamp("COLLISION_KERNEL_NORMALIZE", status="PASS")

if K10.dimensions() != (10, 16) or K10.rank() != 10:
    raise ArithmeticError(
        f"constructed collision kernel has {K10.dimensions()}/rank={K10.rank()}"
    )

# Full exact certificate: modular arithmetic only selected the pivot minor.
if H * K10.transpose() != matrix(QQ, 6, 10):
    raise ArithmeticError("constructed collision kernel fails exact QQ replay")

stamp(
    "COLLISION_KERNEL",
    rank=6,
    kernel=10,
    exact_replay=1,
    status="PASS",
)


def AB_from_Brow(row):
    Bcoef = R(list(row))

    # By linearity:
    #   B*Y/X = sum_i b_i * (U^i*Y/X)
    # and A_columns already contains those sixteen reduced polynomials.
    Acoef = R.zero()
    for i, coefficient in enumerate(row):
        if coefficient:
            Acoef += coefficient * A_columns[i]

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
# IMPORTANT: do not truncate at the original I9* point.
#
# Later infinitely-near centres can have nonzero affine coordinates.  After
# blow-up and translation, terms of high degree in the original chart can
# contribute to low local order.  Therefore root-level total-degree
# truncation is not compatible with this resolved RR calculation.
#
# Keep the functions exact, while retaining the safe optimization that
# shifts only the elementary 16-column basis and applies K10 afterwards.

stamp(
    "LOCAL_EXACT",
    method="EXACT_SHIFTED_FUNCTIONS_WITH_RECURRENCE",
    status="BEGIN",
)

def shift_u_exact(poly):
    # Exact Horner evaluation f(alpha+u), avoiding generic substitution.
    result = S.zero()
    au = S(alpha + u)
    for coefficient in reversed(R(poly).list()):
        result = result * au + QQ(coefficient)
    return S(result)


Zl = shift_u_exact(Z)
Xl = shift_u_exact(X)
Yl = shift_u_exact(Y)

if Zl(0, 0, 0) == 0 or Xl(0, 0, 0) == 0:
    raise ArithmeticError(
        "q24 section is not a unit in the selected I9* chart"
    )

Z2l = S(Zl * Zl)
Z3l = S(Z2l * Zl)

# m=(y+yP)/(x-xP)
m_num = S(y * Z3l + Yl)
m_den = S(Zl * (x * Z2l - Xl))

if m_den(0, 0, 0) == 0:
    raise ArithmeticError(
        "marked-chord denominator is not a local unit"
    )

stamp(
    "LOCAL_EXACT_BASIS",
    method="ONE_A_SHIFT_PLUS_EXACT_RECURRENCE",
    columns=16,
    status="BEGIN",
)

# A_columns[i+1] = U*A_columns[i] - A_columns[i][47]*Z^2
# exactly, because Z^2 is monic of degree 48.
A_local_columns = [shift_u_exact(A_columns[0])]
B_local_columns = [S.one()]
alpha_plus_u = S(alpha + u)

for i in range(15):
    A_local_columns.append(
        S(
            alpha_plus_u * A_local_columns[-1]
            - QQ(A_columns[i][47]) * Z2l
        )
    )
    B_local_columns.append(
        S(alpha_plus_u * B_local_columns[-1])
    )

assert len(A_local_columns) == 16
assert len(B_local_columns) == 16

stamp(
    "LOCAL_EXACT_BASIS",
    columns=16,
    status="PASS",
)

# Apply the 10x16 exact collision kernel only after constructing the
# elementary local columns.  This avoids ten independent huge Taylor shifts.
stamp(
    "LOCAL_NUMERATOR_BASIS",
    method="UNIVARIATE_TRIPLE_RECURRENCE",
    columns=16,
    status="BEGIN",
)

# Every local numerator is linear in x,y:
#
#     N_i = C_i(u) + D_i(u)*x + E_i(u)*y.
#
# Do all heavy arithmetic in a dedicated univariate QQ[u] ring and only
# convert the final ten functions to S=QQ[u,x,y].

LU = PolynomialRing(QQ, "w")
w = LU.gen()

def S_u_to_LU(poly):
    poly = S(poly)
    terms = {}
    for exponent, coefficient in poly.dict().items():
        if exponent[1] != 0 or exponent[2] != 0:
            raise ArithmeticError(
                "expected polynomial involving only local base coordinate u"
            )
        terms[int(exponent[0])] = QQ(coefficient)
    return LU(terms)

def LU_to_S(poly):
    poly = LU(poly)
    return S(sum(
        QQ(coefficient) * u**i
        for i, coefficient in enumerate(poly.list())
        if coefficient
    ))

z  = S_u_to_LU(Zl)
xs = S_u_to_LU(Xl)
ys = S_u_to_LU(Yl)
a0 = S_u_to_LU(A_local_columns[0])

stamp(
    "LOCAL_NUMERATOR_COMMON",
    method="UNIVARIATE_PRODUCTS",
    status="BEGIN",
)

# Common univariate products.
z2 = z*z
z3 = z2*z
z4 = z3*z
z5 = z4*z

xz = xs*z
yz = ys*z

stamp(
    "LOCAL_NUMERATOR_COMMON",
    status="PASS",
)

# N0 = A0*m_den + Z*m_num
#
# m_den = x Z^3 - X Z
# Z*m_num = y Z^4 + Y Z
#
# hence
#   constant = YZ - A0*XZ
#   x coeff  = A0*Z^3
#   y coeff  = Z^4
N_columns = [
    (
        yz - a0*xz,
        a0*z3,
        z4,
    )
]

# Z^2*m_den = x Z^5 - X Z^3.
correction = (
    -xs*z3,   # constant
    z5,       # x
    LU.zero() # y
)

aw = LU(alpha) + w

for i in range(15):
    c = QQ(A_columns[i][47])
    C, D, E = N_columns[-1]

    N_columns.append(
        (
            aw*C - c*correction[0],
            aw*D - c*correction[1],
            aw*E,  # correction has no y term
        )
    )

if len(N_columns) != 16:
    raise AssertionError("elementary numerator recurrence has wrong size")

stamp(
    "LOCAL_NUMERATOR_RECURRENCE",
    columns=16,
    status="PASS",
)

# Apply K10 while STILL univariate.
local_numerators = []

for row in K10.rows():
    C = LU.zero()
    D = LU.zero()
    E = LU.zero()

    for i, coefficient in enumerate(row):
        if coefficient:
            q = QQ(coefficient)
            Ci, Di, Ei = N_columns[i]
            C += q*Ci
            D += q*Di
            E += q*Ei

    # Only now enter the trivariate ring.
    numerator = (
        LU_to_S(C)
        + x*LU_to_S(D)
        + y*LU_to_S(E)
    )

    local_numerators.append(S(numerator))

if len(local_numerators) != 10:
    raise AssertionError(
        "post-collision local basis is not 10-dimensional"
    )

stamp(
    "LOCAL_NUMERATOR_BASIS",
    method="UNIVARIATE_TRIPLE_RECURRENCE",
    columns=16,
    basis=10,
    status="PASS",
)

stamp(
    "LOCAL_EXACT",
    basis=10,
    common_den_unit=1,
    status="PASS",
)


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



def _monomial_exponent(mon):
    data = S(mon).dict()
    if len(data) != 1:
        raise ArithmeticError(f"expected monomial, got {mon}")
    return next(iter(data.keys()))

def _taylor_coefficient(poly, point, exponent):
    derivative = S(poly)
    denominator = 1
    for variable, degree in zip((u,x,y), exponent):
        for unused in range(int(degree)):
            derivative = derivative.derivative(variable)
        denominator *= factorial(int(degree))
    return QQ(derivative(*point)) / QQ(denominator)

def _low_vector_at(poly, point, mons):
    return [
        _taylor_coefficient(
            poly, point, _monomial_exponent(mon)
        )
        for mon in mons
    ]

def local_order_matrix(basis, surface_eq, point, required_order):
    required_order = int(required_order)
    if required_order not in (2,3):
        raise ArithmeticError(
            f"local quotient supports only orders 2/3, got {required_order}"
        )

    mons = monomials_below(required_order)

    # Certify surface multiplicity two using only its 2-jet.
    surf_mons = monomials_below(3)
    surf = _low_vector_at(surface_eq, point, surf_mons)
    surf_degrees = [sum(_monomial_exponent(mon)) for mon in surf_mons]
    if any(surf[i] for i,d in enumerate(surf_degrees) if d < 2):
        raise ArithmeticError("surface order dropped below two")
    if not any(surf[i] for i,d in enumerate(surf_degrees) if d == 2):
        raise ArithmeticError("surface quadratic tangent cone vanished")

    raw = matrix(
        QQ,
        len(mons),
        len(basis),
        lambda row,col:
            _taylor_coefficient(
                basis[col], point, _monomial_exponent(mons[row])
            ),
    )

    if required_order == 2:
        return raw, mons

    # Mod m^3 the only relation is the quadratic tangent cone.
    q = _low_vector_at(surface_eq, point, mons)
    pivot = next((i for i,value in enumerate(q) if value), None)
    if pivot is None:
        raise ArithmeticError("surface relation vanished mod m^3")

    L = matrix(QQ, len(mons)-1, len(mons))
    outrow = 0
    for i in range(len(mons)):
        if i == pivot:
            continue
        L[outrow,i] = q[pivot]
        L[outrow,pivot] = -q[i]
        outrow += 1
    return L*raw, mons

def canonical_after_condition(poly, surface_eq, point, required_order):
    required_order = int(required_order)

    if required_order == 2:
        mons = monomials_below(2)
        if any(_low_vector_at(poly, point, mons)):
            raise ArithmeticError(
                "order-2 kernel vector has nonzero low Taylor coefficient"
            )
        return S(poly)

    if required_order == 3:
        mons = monomials_below(3)
        v = _low_vector_at(poly, point, mons)
        q = _low_vector_at(surface_eq, point, mons)
        pivot = next((i for i,value in enumerate(q) if value), None)
        if pivot is None:
            raise ArithmeticError("missing quadratic surface relation")
        scalar = QQ(v[pivot]) / QQ(q[pivot])
        if any(v[i] != scalar*q[i] for i in range(len(mons))):
            raise ArithmeticError(
                "order-3 kernel vector not proportional to surface 2-jet"
            )
        result = S(poly - scalar*surface_eq)
        if any(_low_vector_at(result, point, mons)):
            raise ArithmeticError(
                "surface subtraction failed to produce order >=3"
            )
        return result

    raise ArithmeticError(
        f"canonical representative supports only 2/3, got {required_order}"
    )


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


def fast_right_kernel_QQ(M, label, pivot_prime=100003):
    # Right kernel of a small huge-coefficient QQ matrix.
    # GF(p) selects independent rows/columns only; the actual solve and
    # verification are exact over QQ.
    Fp_local = GF(ZZ(pivot_prime))

    def red(v):
        v = QQ(v)
        d = ZZ(v.denominator())
        if d % pivot_prime == 0:
            raise ZeroDivisionError(
                f"{label}: denominator divisible by pivot prime {pivot_prime}"
            )
        return Fp_local(ZZ(v.numerator())) / Fp_local(d)

    Mmod = matrix(
        Fp_local, M.nrows(), M.ncols(),
        [red(v) for v in M.list()]
    )
    rank = int(Mmod.rank())

    if rank == 0:
        return 0, identity_matrix(QQ, M.ncols())

    pivot_cols = tuple(int(i) for i in Mmod.pivots())
    if len(pivot_cols) != rank:
        raise ArithmeticError(
            f"{label}: modular pivot-column count {len(pivot_cols)} != rank {rank}"
        )

    Pmod = Mmod.matrix_from_columns(pivot_cols)
    pivot_rows = tuple(int(i) for i in Pmod.transpose().pivots())
    if len(pivot_rows) != rank:
        raise ArithmeticError(
            f"{label}: modular pivot-row count {len(pivot_rows)} != rank {rank}"
        )

    free_cols = tuple(i for i in range(M.ncols()) if i not in pivot_cols)

    print(
        "Q24D12QQ_LOCAL_PIVOTS|"
        f"center={label}|rank={rank}|"
        f"rows={','.join(map(str,pivot_rows))}|"
        f"cols={','.join(map(str,pivot_cols))}|"
        f"free={','.join(map(str,free_cols))}|status=PASS",
        flush=True,
    )

    P = M.matrix_from_rows_and_columns(pivot_rows, pivot_cols)

    if not free_cols:
        K = matrix(QQ, 0, M.ncols())
        return rank, K

    Ffree = M.matrix_from_rows_and_columns(pivot_rows, free_cols)

    print(
        "Q24D12QQ_LOCAL_SOLVE|"
        f"center={label}|size={rank}x{rank}_by_{len(free_cols)}|status=BEGIN",
        flush=True,
    )
    pivot_values = P.solve_right(-Ffree)
    print(
        "Q24D12QQ_LOCAL_SOLVE|"
        f"center={label}|size={rank}x{rank}_by_{len(free_cols)}|status=PASS",
        flush=True,
    )

    rows = []
    for j, free_col in enumerate(free_cols):
        row = [QQ.zero()] * M.ncols()
        row[free_col] = QQ.one()
        for i, pivot_col in enumerate(pivot_cols):
            row[pivot_col] = pivot_values[i, j]
        rows.append(row)

    K = primitive_integer_basis(matrix(QQ, rows))
    if K.dimensions() != (len(free_cols), M.ncols()):
        raise ArithmeticError(f"{label}: constructed kernel has wrong dimensions")
    if M * K.transpose() != matrix(QQ, M.nrows(), K.nrows()):
        raise ArithmeticError(f"{label}: exact local kernel replay failed")

    return rank, K


def divide_power_rr(poly, exceptional, power):
    q = S(poly)
    for unused in range(int(power)):
        q, rem = q.quo_rem(exceptional)
        if rem:
            raise ArithmeticError("RR exceptional power does not divide")
    return S(q)


current_basis = list(local_numerators)

# The full exact surface already certified the blow-up tree.  For the RR
# quotient along the selected four centres, work in the exact truncated
# local ring S / m^9.
current_surface = surface0
transform = identity_matrix(QQ, 10)
condition_ledger = []

for step_index, (center_name, required_order) in enumerate(SELECTED_CLUSTER):
    record = records[center_name]
    point = record["point"]

    print(
        "Q24D12QQ_CENTER_BEGIN|"
        f"center={center_name}|order={required_order}|"
        f"basis={len(current_basis)}|status=BEGIN",
        flush=True,
    )
    M, quotient_monomials = local_order_matrix(
        current_basis, current_surface, point, required_order
    )
    print(
        "Q24D12QQ_LOCAL_MATRIX|"
        f"center={center_name}|rows={M.nrows()}|cols={M.ncols()}|status=PASS",
        flush=True,
    )
    local_rank, kernel = fast_right_kernel_QQ(M, center_name)
    before = len(current_basis)
    after = int(kernel.nrows())

    mod_label = exact_to_mod[center_name]
    expected = expected_mod_ledger.get(mod_label)
    if expected is not None:
        exp = (
            int(expected["dimension_before"]),
            int(expected["local_rank"]),
            int(expected["dimension_after"]),
        )
        got = (before, local_rank, after)
        print(
            "Q24D12QQ_MODULAR_REPLAY|"
            f"mod_center={mod_label}|qq_center={center_name}|"
            f"qq={got[0]},{got[1]},{got[2]}|"
            f"mod={exp[0]},{exp[1]},{exp[2]}|"
            f"status={'PASS' if got==exp else 'MISMATCH'}",
            flush=True,
        )
        if got != exp:
            raise ArithmeticError(
                f"exact/modular RR mismatch at {mod_label}->{center_name}: "
                f"QQ={got}, MOD={exp}"
            )

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

        print(
            "Q24D12QQ_BLOWUP_BEGIN|"
            f"from={center_name}|to={next_name}|chart={kind}|"
            f"basis={len(current_basis)}|status=BEGIN",
            flush=True,
        )
        current_surface = divide_power_rr(
            S(current_surface(*subs)),
            exceptional,
            int(record["multiplicity"]),
        )

        current_basis = [
            divide_power_rr(S(poly(*subs)), exceptional, required_order)
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

B2 = primitive_integer_basis(transform * K10)
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

stamp(
    "PENCIL_ELIMINATION",
    method="DIRECT_CHORD_DISCRIMINANT",
    status="BEGIN",
)

mval = pencil_chord_solution(a0, b0, a1, b1, VF(V))
disc = chord_discriminant(xPV, -yPV, AV, mval)

stamp(
    "PENCIL_ELIMINATION",
    method="DIRECT_CHORD_DISCRIMINANT",
    status="PASS",
)

stamp("QUARTIC_SQUAREFREE", status="BEGIN")
quartic, square_factor = squarefree_binary_quartic(disc, UR)
stamp("QUARTIC_SQUAREFREE", status="PASS")
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

stamp("CHILD_CLASSIFY", status="BEGIN")
classification = classify_finite_short_weierstrass_fibres(VR, jacA, jacB)
stamp("CHILD_CLASSIFY", status="PASS")

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
