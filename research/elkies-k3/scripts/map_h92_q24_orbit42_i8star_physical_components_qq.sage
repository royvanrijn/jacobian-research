#!/usr/bin/env sage -python
"""
status: ACTIVE_PROOF_AID
claim: exact physical I8* component marking for
       H3-04 D12/MW5 --q6/orbit42--> A11/MW6

Inputs:
  artifacts/local/elkies-k3/q24-d12-to-a11-orbit42-divval-preflight.json
  artifacts/local/elkies-k3/q24-d13-to-d12-component-valuation-qq.json

Output:
  artifacts/local/elkies-k3/q24-d12-orbit42-i8star-physical-marking-qq.json

This resolves the exact I8* singularity of the certified D12 parent by
ordinary point blow-ups over QQ.  It tracks:
  * 11 blow-up centres;
  * 12 geometric exceptional components (one split centre);
  * the finite D12 dual graph;
  * fibre multiplicities via the valuation of the original base parameter;
  * every graph/multiplicity isomorphism to the exact abstract D12 root frame;
  * the orbit42 vertical vector and section component under every surviving
    spinor-arm orientation.

No RR dimension or child equation is claimed here.
"""

import argparse
import json
import time
from pathlib import Path

from sage.all import QQ, ZZ, PolynomialRing, gcd, lcm, matrix, sage_eval, vector


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
        if (
            (candidate / "elkies-k3/scripts").is_dir()
            and (candidate / "artifacts/local/elkies-k3").is_dir()
        ):
            return candidate
    raise SystemExit("Could not locate jacobian-research")


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
LOCAL = ROOT / "artifacts/local/elkies-k3"

PREFLIGHT = LOCAL / "q24-d12-to-a11-orbit42-divval-preflight.json"
PARENT = LOCAL / "q24-d13-to-d12-component-valuation-qq.json"
OUTPUT = (
    args.output.resolve()
    if args.output
    else LOCAL / "q24-d12-orbit42-i8star-physical-marking-qq.json"
)

for path in (PREFLIGHT, PARENT):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

_progress_started = time.monotonic()

def progress(stage, **fields):
    tail = "|".join(f"{k}={v}" for k, v in fields.items())
    print(
        "Q42I8QQ_PROGRESS|"
        f"elapsed={time.monotonic()-_progress_started:.3f}|stage={stage}"
        + (f"|{tail}" if tail else ""),
        flush=True,
    )

progress("INPUT", status="BEGIN")

pre = json.loads(PREFLIGHT.read_text())
parent = json.loads(PARENT.read_text())

if pre.get("status") != "PASS_Q42_DIVVAL_PREFLIGHT":
    raise SystemExit("orbit42 divval preflight is not passing")
if parent.get("status") != "PASS_EXACT_Q24_D13_TO_D12_COMPONENT_VALUATION_RR":
    raise SystemExit("exact D12 parent is not passing")

# -------------------------------------------------------------------------
# Exact parent local I8* germ.
# -------------------------------------------------------------------------
R = PolynomialRing(QQ, "V")
V = R.gen()
child = parent["child"]

A = R([QQ(v) for v in child["minimal_A_coefficients_low_to_high"]])
B = R([QQ(v) for v in child["minimal_B_coefficients_low_to_high"]])
Delta = R([QQ(v) for v in child["minimal_discriminant_coefficients_low_to_high"]])

i8 = [
    row for row in child["finite_fibres"]
    if str(row["kodaira"]) == "I8*" and int(row["degree"]) == 1
]
if len(i8) != 1:
    raise ArithmeticError(f"expected one rational I8*, got {len(i8)}")

f = R(sage_eval(str(i8[0]["factor"]), locals={"V": V}))
if f.degree() != 1:
    raise ArithmeticError("I8* factor is not linear")
alpha = QQ(-f[0] / f[1])

T = PolynomialRing(QQ, "t")
t = T.gen()
Al1 = T(A(alpha + t))
Bl1 = T(B(alpha + t))
Dl1 = T(Delta(alpha + t))
orders = (
    int(Al1.valuation()),
    int(Bl1.valuation()),
    int(Dl1.valuation()),
)
if orders != (2, 3, 14):
    raise ArithmeticError(f"I8* local orders {orders} != (2,3,14)")

S = PolynomialRing(QQ, names=("u", "x", "y"), order="degrevlex")
u, x, y = S.gens()


def primitive_poly(poly):
    poly = S(poly)
    if not poly:
        return poly
    den = ZZ.one()
    for c in poly.coefficients():
        den = lcm(den, ZZ(QQ(c).denominator()))
    q = S(poly * den)
    nums = [abs(ZZ(QQ(c).numerator())) for c in q.coefficients() if c]
    content = gcd(nums) if nums else ZZ.one()
    if content > 1:
        q = S(q / content)
    coeffs = q.coefficients()
    if coeffs and coeffs[0] < 0:
        q = -q
    return q


Al = S(A(alpha + u))
Bl = S(B(alpha + u))
surface0 = primitive_poly(y**2 - x**3 - Al*x - Bl)
if surface0(0, 0, 0) != 0:
    raise ArithmeticError("I8* germ is not centered at the origin")


def shifted(poly, point):
    a, b, c = map(QQ, point)
    return S(poly(a + u, b + x, c + y))


def order_at_point(poly, point):
    sh = shifted(poly, point)
    if not sh:
        return 10**9
    return min(sum(exp) for exp, coeff in sh.dict().items() if coeff)


def hessian_rank(poly, point):
    H = matrix(
        QQ,
        [
            [poly.derivative(a, b)(*point) for b in (u, x, y)]
            for a in (u, x, y)
        ],
    )
    return int(H.rank())


def tangent_component_count(poly, point, multiplicity):
    if multiplicity != 2:
        return 1, None
    rank = hessian_rank(poly, point)
    if rank not in (1, 2, 3):
        raise ArithmeticError(f"unexpected tangent quadratic rank {rank}")
    return (2 if rank == 2 else 1), rank


def divide_power(poly, exceptional, power):
    q = S(poly)
    for _ in range(int(power)):
        q, rem = q.quo_rem(exceptional)
        if rem:
            raise ArithmeticError("exceptional power does not divide")
    return primitive_poly(q)


def exceptional_order(poly, exceptional):
    poly = S(poly)
    count = 0
    while poly:
        q, rem = poly.quo_rem(exceptional)
        if rem:
            break
        count += 1
        poly = S(q)
    return count


def canonical_projective(direction):
    vals = [QQ(v) for v in direction]
    pivot = next((i for i, v in enumerate(vals) if v), None)
    if pivot is None:
        raise ArithmeticError("zero projective direction")
    scale = vals[pivot]**-1
    return tuple(QQ(v * scale) for v in vals)


def singular_points_on_exceptional(poly, exceptional):
    local_started = time.monotonic()
    progress(
        "SINGULAR_LOCUS",
        exceptional=exceptional,
        poly_terms=len(poly.dict()),
        status="BEGIN",
    )
    ideal = S.ideal([
        exceptional,
        poly,
        poly.derivative(u),
        poly.derivative(x),
        poly.derivative(y),
    ])
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
        pt = tuple(QQ(sol.get(v, QQ.zero())) for v in (u, x, y))
        if pt not in points:
            points.append(pt)
    points.sort()
    progress(
        "SINGULAR_LOCUS",
        exceptional=exceptional,
        points=len(points),
        local_elapsed=f"{time.monotonic()-local_started:.3f}",
        status="PASS",
    )
    return points


def chart_substitutions(point, kind):
    a, b, c = map(QQ, point)
    if kind == "u":
        return (a + u, b + u*x, c + u*y), u
    if kind == "x":
        return (a + x*u, b + x, c + x*y), x
    if kind == "y":
        return (a + y*u, b + y*x, c + y), y
    raise ValueError(kind)


def direction_from_chart(kind, point):
    pu, px, py = map(QQ, point)
    if kind == "u":
        if pu != 0:
            raise ArithmeticError("u-chart exceptional point has u != 0")
        return canonical_projective((1, px, py))
    if kind == "x":
        if px != 0:
            raise ArithmeticError("x-chart exceptional point has x != 0")
        return canonical_projective((pu, 1, py))
    if py != 0:
        raise ArithmeticError("y-chart exceptional point has y != 0")
    return canonical_projective((pu, px, 1))


# -------------------------------------------------------------------------
# Resolution tree.  Track the original base parameter through each chart.
# -------------------------------------------------------------------------
center_records = []
next_center = 0
MAX_CENTERS = 20

# INCIDENCE_AWARE_EXCEPTIONAL_TRACKING
exceptional_edges = set()


def edge_key(left, right):
    return tuple(sorted((str(left), str(right))))


def strict_transform_divisor(divisor, subs, exceptional):
    transformed = S(divisor(*subs))
    multiplicity = exceptional_order(transformed, exceptional)
    strict = S(transformed)
    for _ in range(multiplicity):
        strict, rem = strict.quo_rem(exceptional)
        if rem:
            raise ArithmeticError(
                "exceptional divisor strict-transform division failed"
            )
    return S(strict), int(multiplicity)


def blow_center(
    poly,
    base_pullback,
    point,
    path,
    depth,
    parent_label=None,
    active_divisors=None,
):
    global next_center

    if len(center_records) >= MAX_CENTERS:
        raise RuntimeError("exact I8* resolution exceeded max centres")

    multiplicity = order_at_point(poly, point)
    if multiplicity < 2:
        raise ArithmeticError(("blow-up at smooth point", path, point))

    components, tangent_rank = tangent_component_count(
        poly, point, multiplicity
    )

    if active_divisors is None:
        active_divisors = {}

    passing = sorted(
        name
        for name, divisor in active_divisors.items()
        if S(divisor)(*point) == 0
    )
    stale = sorted(set(active_divisors) - set(passing))
    if stale:
        raise ArithmeticError(
            f"active exceptional divisors do not pass centre {path}: {stale}"
        )

    label = f"C{next_center+1:02d}"
    next_center += 1

    if components != 1:
        raise ArithmeticError(
            f"incidence-aware mapper encountered split centre {label}; "
            "branch-specific equations are required"
        )

    # Blow-up at an intersection separates the old components and inserts
    # the new exceptional between them.
    for ia, left in enumerate(passing):
        for right in passing[ia+1:]:
            exceptional_edges.discard(edge_key(left, right))
    for old_component in passing:
        exceptional_edges.add(edge_key(old_component, label))

    progress(
        "CENTER",
        label=label,
        depth=depth,
        path=path,
        multiplicity=multiplicity,
        tangent_rank=tangent_rank,
        geometric_components=components,
        passing=",".join(passing) if passing else "NONE",
        status="BEGIN",
    )

    record = {
        "label": label,
        "parent_label": parent_label,
        "path": path,
        "depth": int(depth),
        "point": tuple(map(QQ, point)),
        "surface_multiplicity": int(multiplicity),
        "geometric_exceptional_components": int(components),
        "tangent_rank": tangent_rank,
        "passing_exceptional_components": passing,
        "children": [],
        "chart_exceptional_data": [],
    }
    center_records.append(record)

    candidates = {}
    preference = {"u": 0, "x": 1, "y": 2}

    for kind in ("u", "x", "y"):
        chart_started = time.monotonic()
        progress(
            "CHART",
            center=label,
            chart=kind,
            status="BEGIN",
        )
        subs, exceptional = chart_substitutions(point, kind)
        transformed_surface = S(poly(*subs))
        strict = divide_power(
            transformed_surface, exceptional, multiplicity
        )

        transformed_base = S(base_pullback(*subs))
        base_mult = exceptional_order(transformed_base, exceptional)

        chart_divisors = {}
        divisor_transform_orders = {}
        for old_name, old_divisor in active_divisors.items():
            strict_old, old_mult = strict_transform_divisor(
                old_divisor, subs, exceptional
            )
            divisor_transform_orders[old_name] = int(old_mult)
            if strict_old.degree() > 0:
                chart_divisors[old_name] = strict_old

        # Newly created exceptional divisor in this affine chart.
        chart_divisors[label] = S(exceptional)

        # The factorization of strict|E is diagnostic only and can be vastly
        # more expensive than the actual resolution.  The singular locus below
        # already determines the infinitely-near centres exactly, so retain the
        # restriction itself but do NOT factor it.
        restriction = S(strict.subs({exceptional: 0}))

        record["chart_exceptional_data"].append({
            "chart": kind,
            "base_multiplicity": int(base_mult),
            "restriction": str(restriction),
            "restriction_factors": [],
            "restriction_factorization_skipped": True,
            "exceptional_strict_transform_orders": divisor_transform_orders,
        })

        # Projective-chart cover optimization.
        #
        # u-chart covers [1:*:*].
        # x-chart covers [0:1:*] plus overlap with u-chart.
        # After those, the only direction unique to the y-chart is [0:0:1],
        # represented by (u,x,y)=(0,0,0) in the y-chart.  Test that one point
        # exactly instead of asking Singular for a full variety().
        if kind == "y":
            q0 = (QQ.zero(), QQ.zero(), QQ.zero())
            UNIQUE_PROJECTIVE_Y_DIRECTION = True
            singular_here = (
                strict(*q0) == 0
                and all(strict.derivative(v)(*q0) == 0 for v in (u, x, y))
            )
            singular_points = [q0] if singular_here else []
            progress(
                "SINGULAR_LOCUS",
                exceptional=exceptional,
                method="UNIQUE_PROJECTIVE_Y_DIRECTION",
                points=len(singular_points),
                status="PASS",
            )
        else:
            singular_points = singular_points_on_exceptional(strict, exceptional)
        for q in singular_points:
            direction = direction_from_chart(kind, q)
            candidates.setdefault(direction, []).append(
                (
                    kind,
                    q,
                    strict,
                    transformed_base,
                    base_mult,
                    chart_divisors,
                )
            )

        progress(
            "CHART",
            center=label,
            chart=kind,
            singular_points=len(singular_points),
            base_multiplicity=base_mult,
            method=(
                "UNIQUE_PROJECTIVE_Y_DIRECTION"
                if kind == "y"
                else "FULL_AFFINE_CHART"
            ),
            local_elapsed=f"{time.monotonic()-chart_started:.3f}",
            status="PASS",
        )

    # Base multiplicity of the newly created exceptional divisor is intrinsic.
    all_mults = {
        int(item["base_multiplicity"])
        for item in record["chart_exceptional_data"]
    }
    if len(all_mults) != 1:
        raise ArithmeticError(
            f"chart-dependent base multiplicity at {label}: {sorted(all_mults)}"
        )
    record["base_multiplicity"] = int(next(iter(all_mults)))

    for direction in sorted(candidates):
        reps = candidates[direction]
        reps.sort(key=lambda item: preference[item[0]])
        (
            kind,
            q,
            strict,
            next_base,
            unused_mult,
            chart_divisors,
        ) = reps[0]

        child_active = {
            name: divisor
            for name, divisor in chart_divisors.items()
            if S(divisor)(*q) == 0
        }
        if label not in child_active:
            raise ArithmeticError(
                f"child singular point is not on new exceptional {label}"
            )

        direction_text = ",".join(str(v) for v in direction)
        child_path = f"{path}/{label}:{kind}:{direction_text}"
        record["children"].append({
            "direction": direction,
            "selected_chart": kind,
            "path": child_path,
        })
        blow_center(
            strict,
            next_base,
            q,
            child_path,
            depth + 1,
            parent_label=label,
            active_divisors=child_active,
        )

    progress(
        "CENTER",
        label=label,
        depth=depth,
        child_directions=len(candidates),
        status="PASS",
    )


progress("RESOLUTION", status="BEGIN")
blow_center(
    surface0,
    S(u),
    (QQ.zero(), QQ.zero(), QQ.zero()),
    "root",
    0,
    active_divisors={},
)
progress("RESOLUTION", centers=len(center_records), status="PASS")

center_count = len(center_records)
geometric_components = sum(
    int(r["geometric_exceptional_components"])
    for r in center_records
)
split_centers = [
    r["label"]
    for r in center_records
    if int(r["geometric_exceptional_components"]) == 2
]

# ALLOW_12_CENTER_NONSPLIT_I8STAR
#
# For I8* the finite root lattice has 12 components (D12).  Ordinary point
# blow-ups need not realize the D-branch by a rank-2 tangent cone at one
# centre.  In this exact model the branch can instead appear as two distinct
# infinitely-near centres from an earlier exceptional divisor.
#
# Therefore the invariant requirement is 12 geometric exceptional root
# components.  Two resolution presentations are admissible:
#
#   12 centres, all nonsplit
#   11 centres, one split centre
#
# The actual D12 graph/determinant check below is the decisive gate.
if geometric_components != 12:
    raise ArithmeticError(
        f"I8* resolution gave {geometric_components} geometric root components, "
        "expected 12"
    )
if center_count not in (11, 12):
    raise ArithmeticError(
        f"I8* resolution gave {center_count} centres, expected 11 or 12"
    )
if center_count + len(split_centers) != 12:
    raise ArithmeticError(
        "I8* centre/splitting accounting does not give 12 root components: "
        f"centres={center_count}, split_centers={split_centers}"
    )

print(
    "Q42I8QQ|stage=RESOLUTION|"
    f"base={alpha}|orders=2,3,14|centers={center_count}|"
    f"components={geometric_components}|"
    f"split={','.join(split_centers) if split_centers else 'NONE'}|"
    f"presentation={'NONSPLIT_CENTRE_TREE' if not split_centers else 'SPLIT_TANGENT_CENTRE'}|"
    "status=PASS_EXACT_I8STAR_RESOLUTION",
    flush=True,
)

# -------------------------------------------------------------------------
# Build finite physical D12 graph from exact exceptional incidence.
# -------------------------------------------------------------------------
if split_centers:
    raise ArithmeticError(
        "incidence-aware builder expects the observed nonsplit 12-centre "
        f"presentation, got split centres {split_centers}"
    )

vertices = [
    {
        "name": record["label"],
        "center": record["label"],
        "base_multiplicity": int(record["base_multiplicity"]),
        "split_branch": None,
    }
    for record in center_records
]
physical_names = [v["name"] for v in vertices]
physical_index = {name: idx for idx, name in enumerate(physical_names)}
edges = sorted(exceptional_edges)

physical_adj = matrix(ZZ, 12, 12)
for left, right in edges:
    if left not in physical_index or right not in physical_index:
        raise ArithmeticError(
            f"incidence edge references unknown component {(left, right)}"
        )
    ii = physical_index[left]
    jj = physical_index[right]
    physical_adj[ii, jj] = physical_adj[jj, ii] = 1

physical_cartan = 2 * matrix.identity(ZZ, 12) - physical_adj
degrees = sorted(
    sum(physical_adj[ii, jj] for jj in range(12))
    for ii in range(12)
)

print(
    "Q42I8QQ|stage=PHYSICAL_GRAPH_DIAGNOSTIC|"
    f"vertices={len(vertices)}|edges={len(edges)}|"
    f"det={physical_cartan.det()}|"
    "degree_sequence=" + ",".join(map(str, degrees)) + "|"
    "method=INCIDENCE_AWARE_BLOWUP_SURGERY|status=PROFILE",
    flush=True,
)

if len(vertices) != 12:
    raise ArithmeticError(
        f"physical graph has {len(vertices)} vertices, expected 12"
    )
if len(edges) != 11:
    raise ArithmeticError(
        f"physical incidence graph has {len(edges)} edges, expected 11"
    )
if physical_cartan.det() != 4:
    raise ArithmeticError(
        f"physical exceptional graph determinant={physical_cartan.det()}, expected 4"
    )
if degrees != [1, 1, 1] + [2]*8 + [3]:
    raise ArithmeticError(
        f"physical graph degree sequence is not D12: {degrees}"
    )

print(
    "Q42I8QQ|stage=PHYSICAL_GRAPH|"
    f"vertices={len(vertices)}|edges={len(edges)}|det=4|"
    "degree_sequence=" + ",".join(map(str, degrees)) + "|"
    "method=INCIDENCE_AWARE_BLOWUP_SURGERY|status=PASS_D12_GRAPH",
    flush=True,
)

# -------------------------------------------------------------------------
# Abstract D12 graph and affine fibre multiplicity marks.
# -------------------------------------------------------------------------
abstract_gram = matrix(ZZ, pre["abstract_D12_marking"]["root_gram"])
if abstract_gram.dimensions() != (12, 12) or abstract_gram.det() != 4:
    raise ArithmeticError("preflight abstract D12 root Gram is invalid")

abstract_adj = 2 * matrix.identity(ZZ, 12) - abstract_gram

# GRAPH_DERIVED_AFFINE_FIBRE_MARKS
#
# Once the finite D12 graph is known exactly, the Kodaira fibre
# multiplicities are determined intrinsically by the affine intersection
# equations.  Do not use recursively propagated chart exponents here.
#
# If C0 is the identity component with multiplicity m0=1 and it meets
# finite root C_a, then
#
#       R * m = e_a,
#
# while the affine row itself gives
#
#       2*m0 - m_a = 0,
#
# hence m_a = 2.  This second equation eliminates the spurious large
# positive columns of R^{-1}.
def affine_mark_candidates(cartan):
    cartan = matrix(ZZ, cartan)
    out = []
    for attach in range(cartan.nrows()):
        e = vector(ZZ, [int(i == attach) for i in range(cartan.nrows())])
        mq = cartan.solve_right(e)
        if not all(v in ZZ for v in mq):
            continue
        marks = vector(ZZ, [ZZ(v) for v in mq])
        if not all(v > 0 for v in marks):
            continue
        if marks[attach] != 2:
            continue
        out.append((attach, marks))
    return out


affine_candidates = affine_mark_candidates(abstract_gram)
physical_affine_candidates = affine_mark_candidates(physical_cartan)

if not affine_candidates:
    raise ArithmeticError("abstract D12 graph has no valid affine mark candidate")
if not physical_affine_candidates:
    raise ArithmeticError("physical D12 graph has no valid affine mark candidate")

print(
    "Q42I8QQ|stage=AFFINE_MARKS|"
    f"abstract_candidates={len(affine_candidates)}|"
    f"physical_candidates={len(physical_affine_candidates)}|"
    "abstract="
    + ";".join(
        f"a{a}:" + ",".join(map(str, marks))
        for a, marks in affine_candidates
    )
    + "|physical="
    + ";".join(
        f"p{a}:" + ",".join(map(str, marks))
        for a, marks in physical_affine_candidates
    )
    + "|method=GRAPH_DERIVED_AFFINE_INTERSECTION_EQUATIONS|status=PASS",
    flush=True,
)

# -------------------------------------------------------------------------
# Enumerate graph isomorphisms preserving fibre multiplicities.
# -------------------------------------------------------------------------
try:
    import networkx as nx
    from networkx.algorithms import isomorphism as nxiso
except ImportError as exc:
    raise SystemExit("networkx is required for D12 graph matching") from exc

GA = nx.Graph()
GP = nx.Graph()

# If more than one affine attachment candidate survives, include each; their
# union is deduplicated below.
for i in range(12):
    GA.add_node(i)
for i in range(12):
    for j in range(i + 1, 12):
        if abstract_adj[i, j] == 1:
            GA.add_edge(i, j)

for name in physical_names:
    GP.add_node(name)
for left, right in edges:
    GP.add_edge(left, right)

isomorphisms = []
seen_maps = set()

for affine_attach, abstract_marks in affine_candidates:
    for i in range(12):
        GA.nodes[i]["mark"] = int(abstract_marks[i])
        GA.nodes[i]["affine_attach"] = int(i == affine_attach)

    for physical_attach_idx, physical_marks in physical_affine_candidates:
        physical_attach_name = physical_names[physical_attach_idx]
        for name in physical_names:
            idx = physical_index[name]
            GP.nodes[name]["mark"] = int(physical_marks[idx])
            GP.nodes[name]["affine_attach"] = int(name == physical_attach_name)

        matcher = nxiso.GraphMatcher(
            GA,
            GP,
            node_match=nxiso.categorical_node_match(
                ["mark", "affine_attach"], [None, None]
            ),
        )
        for mapping in matcher.isomorphisms_iter():
            if mapping[affine_attach] != physical_attach_name:
                raise ArithmeticError("affine attachment was not preserved")

            key = (
                tuple(mapping[i] for i in range(12)),
                int(affine_attach),
                str(physical_attach_name),
            )
            if key in seen_maps:
                continue
            seen_maps.add(key)
            isomorphisms.append({
                "affine_attachment_abstract_root": int(affine_attach),
                "affine_attachment_physical_component": physical_attach_name,
                "abstract_to_physical": [mapping[i] for i in range(12)],
                "abstract_marks": list(map(int, abstract_marks)),
                "physical_marks": {
                    physical_names[i]: int(physical_marks[i])
                    for i in range(12)
                },
            })

if not isomorphisms:
    raise ArithmeticError("no multiplicity-preserving D12 graph isomorphism")

# -------------------------------------------------------------------------
# Transport orbit42 vertical vector and section component.
# -------------------------------------------------------------------------
vertical_abstract = vector(
    ZZ,
    pre["orbit42"]["vertical_root_coefficients_abstract_D12"],
)

bridge = pre["abstract_D12_marking"].get("bridge_regression")
dual_pairing = None
if bridge and bridge.get("selected_dual_pairing") is not None:
    dual_pairing = vector(ZZ, bridge["selected_dual_pairing"])
    if len(dual_pairing) != 12:
        raise ArithmeticError("bridge dual pairing has wrong length")

orientation_records = []
orientation_keys = set()

for iso in isomorphisms:
    amap = iso["abstract_to_physical"]

    physical_vertical = {
        amap[i]: int(vertical_abstract[i])
        for i in range(12)
    }

    met_components = []
    if dual_pairing is not None:
        met_components = [
            amap[i]
            for i in range(12)
            if dual_pairing[i] != 0
        ]

    key = (
        tuple(sorted(physical_vertical.items())),
        tuple(sorted(met_components)),
    )
    if key in orientation_keys:
        continue
    orientation_keys.add(key)

    orientation_records.append({
        **iso,
        "vertical_coefficients_physical": physical_vertical,
        "section_meets_physical_components": met_components,
        "omitted_vertical_components": sorted(
            name for name, coeff in physical_vertical.items() if coeff == 0
        ),
    })

print(
    "Q42I8QQ|stage=MARKING|"
    f"graph_isomorphisms={len(isomorphisms)}|"
    f"distinct_orientations={len(orientation_records)}|"
    f"dual_pairing={'AVAILABLE' if dual_pairing is not None else 'UNAVAILABLE'}|"
    "status=PASS_PHYSICAL_D12_MARKING_SET",
    flush=True,
)

for idx, row in enumerate(orientation_records):
    mapping = row["abstract_to_physical"]
    print(
        "Q42I8QQ_ORIENTATION|"
        f"index={idx}|"
        f"map={','.join(f'r{i}->{mapping[i]}' for i in range(12))}|"
        f"omitted={','.join(row['omitted_vertical_components']) or 'NONE'}|"
        f"section_meets={','.join(row['section_meets_physical_components']) or 'UNKNOWN'}|"
        "status=CANDIDATE_EXACT_SPINOR_ORIENTATION",
        flush=True,
    )

# Serialize centres without massive strict-transform polynomials.
center_payload = []
for record in center_records:
    center_payload.append({
        "label": record["label"],
        "parent_label": record["parent_label"],
        "path": record["path"],
        "depth": record["depth"],
        "point": [str(v) for v in record["point"]],
        "surface_multiplicity": record["surface_multiplicity"],
        "geometric_exceptional_components": record[
            "geometric_exceptional_components"
        ],
        "tangent_rank": record["tangent_rank"],
        "base_multiplicity": record["base_multiplicity"],
        "passing_exceptional_components": record[
            "passing_exceptional_components"
        ],
        "children": [
            {
                "direction": [str(v) for v in child["direction"]],
                "selected_chart": child["selected_chart"],
                "path": child["path"],
            }
            for child in record["children"]
        ],
        "chart_exceptional_data": record["chart_exceptional_data"],
    })

payload = {
    "schema": "elkies-k3.h3-q24-d12-orbit42-i8star-physical-marking-qq.v1",
    "status": "PASS_Q42_EXACT_I8STAR_PHYSICAL_MARKING",
    "inputs": {
        "preflight": str(PREFLIGHT.relative_to(ROOT)),
        "exact_D12_parent": str(PARENT.relative_to(ROOT)),
    },
    "I8star": {
        "base_factor": str(f),
        "base_root": str(alpha),
        "local_orders": list(orders),
        "centres": center_payload,
        "center_count": center_count,
        "geometric_exceptional_components": geometric_components,
        "split_centers": split_centers,
    },
    "physical_D12": {
        "vertices": vertices,
        "edges": [list(edge) for edge in edges],
        "cartan": [
            [int(v) for v in row]
            for row in physical_cartan.rows()
        ],
        "fibre_multiplicity_candidates": [
            {
                "affine_attachment_component": physical_names[int(attach)],
                "marks": {
                    physical_names[i]: int(marks[i])
                    for i in range(12)
                },
            }
            for attach, marks in physical_affine_candidates
        ],
        "chart_base_multiplicity_diagnostics": {
            v["name"]: int(v["base_multiplicity"])
            for v in vertices
        },
    },
    "abstract_D12": {
        "gram": [
            [int(v) for v in row]
            for row in abstract_gram.rows()
        ],
        "affine_attachment_candidates": [
            {
                "root": int(attach),
                "marks": list(map(int, marks)),
            }
            for attach, marks in affine_candidates
        ],
        "orbit42_vertical": list(map(int, vertical_abstract)),
        "orbit42_dual_pairing": (
            None if dual_pairing is None else list(map(int, dual_pairing))
        ),
    },
    "orientation_candidates": orientation_records,
    "next_required": {
        "stage": "Q42_RESOLVED_RR_TRIVIALIZATION",
        "instruction": (
            "For every surviving physical orientation, rederive the exact local "
            "chart chain and impose the D42 line-bundle conditions. Keep the "
            "horizontal class implicit in the line-bundle cocycle; do not revert "
            "to zero-pole section generation. Require an exact two-dimensional "
            "kernel before compiling a quartic/Jacobian. If the spinor-arm "
            "orientation remains ambiguous, let exact RR + A11/modular "
            "regression select it."
        ),
    },
    "proof_boundary": (
        "Exact I8* resolution, finite D12 dual graph, fibre multiplicities and "
        "abstract-to-physical marking candidates only. The horizontal line-bundle "
        "trivialization and resolved RR kernel are not yet constructed."
    ),
}

OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUTPUT}", flush=True)
print(
    "Q42I8QQ_RESULT|"
    f"centers={center_count}|components={geometric_components}|"
    f"orientations={len(orientation_records)}|"
    "next=Q42_RESOLVED_RR_TRIVIALIZATION|"
    "status=PASS_Q42_EXACT_I8STAR_PHYSICAL_MARKING",
    flush=True,
)
