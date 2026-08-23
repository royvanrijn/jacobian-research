#!/usr/bin/env sage -python
"""
Recover the effective D13 chamber of the H3 q24 divisor and compare its exact
vertical correction with the actual resolved I9* component atlas.

This is deliberately a diagnostic bridge:
  deterministic D13 basis -> q24-dominant/effective D13 chamber
                           -> resolved I9* exceptional components.

Key observation:
`vertical_root_coefficients` are NOT component multiplicities in the effective
I9* chamber.  They live in the deterministic D13 simple-root basis constructed
by close_h92_q8_q24_by_q6_translation.sage.  We first Weyl-reduce the FULL q24
child-frame vector cD[2:] (including its MW lift) until its pairings with all
13 D13 simple roots are nonnegative.  The inverse Weyl image of those roots is
then the q24 effective simple system.  Effective fibre components are their
negatives.

We then:
  * express the exact root correction in that effective component basis;
  * compute the vanishing cycle forced by negative component coefficients;
  * test anti-nefness and the Du Val colength Z*C*Z/2;
  * rebuild the ordinary-blowup incidence graph from the existing modular
    I9* resolution artifact;
  * extract v_E(u), v_E(x), v_E(y) from the tracked origin maps;
  * enumerate graph identifications and a monomial approximation to the
    corresponding complete ideal.

A colength of 8 is the expected missing local RR codimension after the
16 -> 10 smooth-collision reduction.
"""

import argparse
import contextlib
import io
import json
import sys
from itertools import combinations, product
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix, vector


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
parser.add_argument("--monomial-u", type=int, default=18)
parser.add_argument("--monomial-x", type=int, default=8)
parser.add_argument("--monomial-y", type=int, default=4)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
SCRIPTS = ROOT / "elkies-k3/scripts"
LOCAL = ROOT / "artifacts/local/elkies-k3"
CLOSE = SCRIPTS / "close_h92_q8_q24_by_q6_translation.sage"
RESOLUTION_SCRIPT = SCRIPTS / "derive_h92_q24_i9star_resolution_modp.sage"
RESOLUTION = LOCAL / f"q24-i9star-resolution-mod-{args.prime}.json"

if not CLOSE.exists():
    raise SystemExit(f"Missing {CLOSE}")
if not RESOLUTION.exists():
    raise SystemExit(
        f"Missing {RESOLUTION}\n"
        f"Run first:\n  sage -python {RESOLUTION_SCRIPT.relative_to(ROOT)} "
        f"--prime {args.prime}"
    )

# ---------------------------------------------------------------------------
# Execute the exact q24 translation in a private scope so we use the precise
# deterministic D13 basis in which vr was solved.
# ---------------------------------------------------------------------------
saved_argv = list(sys.argv)
scope = {"__name__": "__embedded_q24_chamber_probe__"}
captured = io.StringIO()
try:
    sys.argv = [str(CLOSE)]
    with contextlib.redirect_stdout(captured):
        exec(compile(CLOSE.read_text(), str(CLOSE), "exec"), scope)
finally:
    sys.argv = saved_argv

required = ("adapted", "cD", "vr", "vf")
missing = [name for name in required if name not in scope]
if missing:
    raise SystemExit("q24 translation scope lacks: " + ",".join(missing))

adapted = matrix(ZZ, scope["adapted"])
cD = vector(ZZ, scope["cD"])
vr = vector(ZZ, scope["vr"])
vf = ZZ(scope["vf"])
assert adapted.dimensions() == (17, 17)
assert len(cD) == 19
assert len(vr) == 13
root = adapted[:13, :13]
assert root.det() == 4

l_original = vector(ZZ, cD[2:])
assert len(l_original) == 17
root_norm = ZZ(l_original * adapted * l_original)

# ---------------------------------------------------------------------------
# Weyl-reduce the FULL q24 frame vector to the dominant D13 chamber.
# Reflection in simple root e_i: l -> l - <l,e_i> e_i.
# ---------------------------------------------------------------------------
l = vector(ZZ, l_original)
weyl = []
for iteration in range(10000):
    pairings = vector(ZZ, l * adapted[:, :13])
    negative = [i for i, value in enumerate(pairings) if value < 0]
    if not negative:
        break
    i = negative[0]
    pairing = ZZ(pairings[i])
    before_norm = ZZ(l * adapted * l)
    l[i] -= pairing
    assert ZZ(l * adapted * l) == before_norm
    weyl.append((i, int(pairing)))
else:
    raise ArithmeticError("D13 Weyl reduction did not terminate")

dominant_pairings = vector(ZZ, l * adapted[:, :13])
assert all(value >= 0 for value in dominant_pairings)
assert ZZ(l * adapted * l) == root_norm

# w(l_original)=l.  Effective simple roots for the ORIGINAL q24 class are
# beta_j = w^{-1}(e_j).  Apply inverse reflections in reverse order.
B_eff_rows = []
for j in range(13):
    beta = vector(ZZ, [ZZ(i == j) for i in range(17)])
    for i, unused_pairing in reversed(weyl):
        pairing = ZZ(beta * adapted.column(i))
        beta[i] -= pairing
    assert beta * adapted * beta == 2
    B_eff_rows.append(beta)

B_eff = matrix(ZZ, [list(row) for row in B_eff_rows])
C_eff = B_eff * adapted * B_eff.transpose()
assert C_eff == root

replayed = vector(ZZ, l_original * adapted * B_eff.transpose())
assert replayed == dominant_pairings

print(
    "Q24CHAMBER|"
    f"weyl_steps={len(weyl)}|"
    f"dominant_pairings={','.join(map(str, dominant_pairings))}|"
    f"positive_hits={sum(value > 0 for value in dominant_pairings)}|"
    f"sum_hits={sum(dominant_pairings)}|status=PASS",
    flush=True,
)

# ---------------------------------------------------------------------------
# Express exact deterministic-basis vertical root correction in the effective
# simple-root system.  Effective geometric components are C_j=-beta_j.
#
#   rootpart = sum a_j beta_j = sum (-a_j) C_j
#
# so divisor coefficient along C_j is d_j=-a_j.
# Negative d_j forces vanishing v_E(f) >= a_j.
# ---------------------------------------------------------------------------
vr17 = vector(ZZ, list(vr) + [0] * 4)
a_q = matrix(QQ, B_eff.transpose()).solve_right(vector(QQ, vr17))
assert all(value in ZZ for value in a_q)
a = vector(ZZ, a_q)
assert vector(ZZ, a) * B_eff == vr17

component_coefficients = vector(ZZ, [-value for value in a])
required_vanishing = vector(ZZ, [max(ZZ(0), value) for value in a])

cycle_pairings = vector(ZZ, required_vanishing * C_eff)
anti_nef = all(value >= 0 for value in cycle_pairings)
cycle_square_positive = ZZ(required_vanishing * C_eff * required_vanishing)
colength_q = QQ(cycle_square_positive) / 2
colength_integral = colength_q in ZZ
colength = ZZ(colength_q) if colength_integral else None

print(
    "Q24EFFECTIVE_VERTICAL|"
    f"fibre_twist_elsewhere={vf}|"
    f"root_coeff_effective_beta={','.join(map(str, a))}|"
    f"component_coeff={','.join(map(str, component_coefficients))}|"
    f"required_vanishing={','.join(map(str, required_vanishing))}|"
    f"cycle_pairings={','.join(map(str, cycle_pairings))}|"
    f"anti_nef={int(anti_nef)}|"
    f"cycle_square={cycle_square_positive}|"
    f"duval_colength={colength if colength is not None else colength_q}|"
    f"expected=8|status={'HIT_COLENGTH8' if colength == 8 else 'DIAGNOSTIC'}",
    flush=True,
)

# ---------------------------------------------------------------------------
# Rebuild the actual blow-up graph from the existing resolution artifact.
# ---------------------------------------------------------------------------
resolution = json.loads(RESOLUTION.read_text())
if not str(resolution.get("status", "")).startswith("PASS_"):
    raise SystemExit(f"Resolution artifact is not passing: {resolution.get('status')}")
centers = resolution["centers"]
if len(centers) != 13:
    raise ArithmeticError(f"expected 13 blowup centers, got {len(centers)}")

vertices = []
edges = set()
self_intersection = {}
active_history = {}

def edge_key(a, b):
    return tuple(sorted((a, b)))

for record in centers:
    new = record["label"]
    active = list(record.get("active_components", []))
    active_history[new] = active
    vertices.append(new)

    # Blowing up an intersection separates every pair of components through
    # that center, and the new exceptional meets all their strict transforms.
    for left, right in combinations(active, 2):
        edges.discard(edge_key(left, right))
    for old in active:
        if old.startswith("E"):
            if old not in self_intersection:
                raise ArithmeticError(f"unknown active exceptional {old}")
            self_intersection[old] -= 1
        edges.add(edge_key(new, old))
    self_intersection[new] = -1

exceptional_edges = {
    edge for edge in edges if edge[0].startswith("E") and edge[1].startswith("E")
}
all_minus_two = all(self_intersection[name] == -2 for name in vertices)
tree_like = len(exceptional_edges) == 12

# Chart Cartan from incidence, independently of the abstract D13 Cartan.
index = {name: i for i, name in enumerate(vertices)}
C_chart = matrix(ZZ, 13, 13)
for i in range(13):
    C_chart[i, i] = 2
for left, right in exceptional_edges:
    i, j = index[left], index[right]
    C_chart[i, j] = C_chart[j, i] = -1

print(
    "Q24ATLAS_GRAPH|"
    f"components={len(vertices)}|edges={len(exceptional_edges)}|"
    f"all_minus_two={int(all_minus_two)}|tree_like={int(tree_like)}|"
    f"det={C_chart.det()}|"
    f"status={'PASS_D13_GRAPH' if all_minus_two and tree_like and C_chart.det()==4 else 'CHECK_GRAPH'}",
    flush=True,
)

# ---------------------------------------------------------------------------
# Extract valuations of original local u,x,y along each newly-created
# exceptional from the u-chart's tracked origin map.
# ---------------------------------------------------------------------------
p = ZZ(args.prime)
F = GF(p)
S3 = PolynomialRing(F, names=("u", "x", "y"), order="degrevlex")
uu, xx, yy = S3.gens()

def variable_order(poly, variable_index):
    poly = S3(poly)
    if not poly:
        return 10**9
    return min(
        exponent[variable_index]
        for exponent, coefficient in poly.dict().items()
        if coefficient
    )

coordinate_values = {}
for record in centers:
    u_chart = next(
        (chart for chart in record.get("charts", []) if chart.get("chart") == "u"),
        None,
    )
    if u_chart is None:
        raise ArithmeticError(f"{record['label']} has no u blow-up chart")
    origin_map = [S3(text) for text in u_chart["origin_map"]]
    values = tuple(variable_order(expr, 0) for expr in origin_map)
    coordinate_values[record["label"]] = values

print(
    "Q24ATLAS_VALUATIONS|"
    + "|".join(
        f"{name}={coordinate_values[name][0]},{coordinate_values[name][1]},{coordinate_values[name][2]}"
        for name in vertices
    )
    + "|status=PASS",
    flush=True,
)

# ---------------------------------------------------------------------------
# Enumerate graph isomorphisms abstract effective D13 -> chart D13.
# Generic backtracking; D13 should leave only its small diagram ambiguity.
# ---------------------------------------------------------------------------
abstract_adj = {
    i: {j for j in range(13) if i != j and C_eff[i, j] == -1}
    for i in range(13)
}
chart_adj = {
    name: {
        other
        for other in vertices
        if other != name and edge_key(name, other) in exceptional_edges
    }
    for name in vertices
}

abstract_signature = {
    i: (len(abstract_adj[i]), sorted(len(abstract_adj[j]) for j in abstract_adj[i]))
    for i in range(13)
}
chart_signature = {
    name: (len(chart_adj[name]), sorted(len(chart_adj[j]) for j in chart_adj[name]))
    for name in vertices
}

candidates = {
    i: [name for name in vertices if chart_signature[name] == abstract_signature[i]]
    for i in range(13)
}
order = sorted(range(13), key=lambda i: (len(candidates[i]), -len(abstract_adj[i])))

isomorphisms = []
def backtrack(position, mapping, used):
    if position == len(order):
        isomorphisms.append(dict(mapping))
        return
    i = order[position]
    for name in candidates[i]:
        if name in used:
            continue
        ok = True
        for j, other in mapping.items():
            if ((j in abstract_adj[i]) != (other in chart_adj[name])):
                ok = False
                break
        if not ok:
            continue
        mapping[i] = name
        used.add(name)
        backtrack(position + 1, mapping, used)
        used.remove(name)
        del mapping[i]

backtrack(0, {}, set())
if not isomorphisms:
    print("Q24ATLAS_ISO|count=0|status=NO_GRAPH_ISOMORPHISM", flush=True)
else:
    print(f"Q24ATLAS_ISO|count={len(isomorphisms)}|status=PASS", flush=True)

# ---------------------------------------------------------------------------
# For every orientation, map the candidate vanishing cycle to chart order and
# find minimal monomials u^a x^b y^c whose exceptional valuations dominate it.
# This is only the monomial subideal; branch-sensitive binomials may be needed.
# ---------------------------------------------------------------------------
orientation_records = []
for iso_index, iso in enumerate(isomorphisms, start=1):
    chart_cycle = {iso[i]: int(required_vanishing[i]) for i in range(13)}

    monomials = []
    for au in range(args.monomial_u + 1):
        for ax in range(args.monomial_x + 1):
            for ay in range(args.monomial_y + 1):
                if au == ax == ay == 0:
                    continue
                good = True
                vals = {}
                for name in vertices:
                    vu, vx, vy = coordinate_values[name]
                    value = au * vu + ax * vx + ay * vy
                    vals[name] = value
                    if value < chart_cycle[name]:
                        good = False
                        break
                if good:
                    monomials.append((au, ax, ay, vals))

    # Minimal under monomial divisibility.
    minimal = []
    for item in monomials:
        au, ax, ay, unused_vals = item
        dominated = False
        for other in monomials:
            bu, bx, by, unused_other_vals = other
            if (bu, bx, by) == (au, ax, ay):
                continue
            if bu <= au and bx <= ax and by <= ay:
                dominated = True
                break
        if not dominated:
            minimal.append(item)

    minimal.sort(key=lambda item: (sum(item[:3]), item[:3]))
    minimal_labels = [
        "u^{}*x^{}*y^{}".format(au, ax, ay)
        for au, ax, ay, unused in minimal
    ]

    print(
        "Q24ATLAS_ORIENTATION|"
        f"i={iso_index}|"
        f"map={';'.join(f'R{k+1}->{iso[k]}' for k in range(13))}|"
        f"cycle={','.join(f'{name}:{chart_cycle[name]}' for name in vertices)}|"
        f"minimal_monomials={','.join(minimal_labels[:16]) or 'NONE'}|"
        f"count={len(minimal)}|status=PASS",
        flush=True,
    )

    orientation_records.append(
        {
            "index": iso_index,
            "effective_root_to_chart": {f"R{i+1}": iso[i] for i in range(13)},
            "chart_cycle": chart_cycle,
            "minimal_monomials": minimal_labels,
        }
    )

payload = {
    "schema": "elkies-k3.h3-q24-effective-d13-chamber-cycle-modp.v1",
    "status": (
        "HIT_EXPECTED_LOCAL_COLENGTH_8"
        if colength == 8
        else "PASS_EFFECTIVE_D13_CHAMBER_DIAGNOSTIC"
    ),
    "prime": int(p),
    "q24": {
        "weyl_step_count": len(weyl),
        "weyl_steps_1_based": [[i + 1, pairing] for i, pairing in weyl],
        "dominant_simple_pairings": list(map(int, dominant_pairings)),
        "vertical_fibre_coefficient_moved_elsewhere": int(vf),
        "deterministic_root_coefficients": list(map(int, vr)),
        "effective_beta_coefficients": list(map(int, a)),
        "effective_component_coefficients": list(map(int, component_coefficients)),
        "required_vanishing_cycle": list(map(int, required_vanishing)),
        "cycle_pairings": list(map(int, cycle_pairings)),
        "anti_nef": bool(anti_nef),
        "cycle_square_positive": int(cycle_square_positive),
        "duval_colength": int(colength) if colength is not None else str(colength_q),
    },
    "resolution": {
        "artifact": str(RESOLUTION.relative_to(ROOT)),
        "all_exceptionals_minus_two": bool(all_minus_two),
        "exceptional_edge_count": len(exceptional_edges),
        "chart_cartan_determinant": int(C_chart.det()),
        "coordinate_valuations_u_x_y": {
            name: list(map(int, coordinate_values[name])) for name in vertices
        },
        "graph_isomorphism_count": len(isomorphisms),
        "orientations": orientation_records,
    },
    "boundary": (
        "This identifies the q24 effective D13 chamber and the divisorial "
        "vanishing cycle. The monomial list is only a subideal diagnostic; "
        "a complete resolved marked-chord module may require branch-sensitive "
        "binomial generators and chart-overlap conditions."
    ),
}
OUT = (
    args.output.resolve()
    if args.output
    else LOCAL / f"q24-effective-d13-chamber-cycle-mod-{p}.json"
)
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUT}", flush=True)
print(
    "Q24CHAMBER_RESULT|"
    f"weyl_steps={len(weyl)}|positive_hits={sum(v>0 for v in dominant_pairings)}|"
    f"anti_nef={int(anti_nef)}|colength={colength if colength is not None else colength_q}|"
    f"atlas_iso={len(isomorphisms)}|"
    f"status={payload['status']}",
    flush=True,
)
