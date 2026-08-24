#!/usr/bin/env sage -python
"""
status: HISTORICAL_DIAGNOSTIC
claim: unreachable orientation experiment for the rejected fast-q6 route.

Inputs:
  q24-d12-orbit42-fast-transport-qq.json
  q24-d12-orbit42-i8star-physical-marking-qq.json
  q24-d12-to-a11-orbit42-divval-preflight.json

This does NOT rerun a resolved-RR matrix.

Why kernel/h0=2 is exact without that matrix:
  * the selected D42 class is already certified primitive, nef, isotropic;
  * on a K3 such a class gives an elliptic pencil, h0(D42)=2;
  * the exact bridge gives D42=O+P+V with fibre_twist=0;
  * this script tracks the exact P section through the exact I8* blow-up tree
    and selects the physical spinor orientation;
  * the direct chord slope is the nonconstant second pencil section;
  * the transport artifact independently requires the compiled child to be
    A11/MW6.

The only ambiguity from physical D12 marking is C10 <-> C11.  Run this script
with --orientation 0 and --orientation 1 in parallel.

Final status is emitted only if:
  omitted == section_meets == actual exact section component,
  V has coefficient 0 there and 1 on all other finite D12 components,
  corrected profile is (height,corr,P.O,twist)=(7,3,3,0),
  and the direct chord compiler produced A11 root data (11,132,12).
"""

import argparse
import json
from pathlib import Path

from sage.all import LaurentSeriesRing, PolynomialRing, QQ


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
            and (c / "artifacts/local/elkies-k3").is_dir()
        ):
            return c
    raise SystemExit("Could not locate jacobian-research")


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--orientation", type=int, required=True, choices=(0, 1))
parser.add_argument("--candidate", type=int, help="specific passing A11 candidate index")
parser.add_argument("--precision", type=int, default=40)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
LOCAL = ROOT / "artifacts/local/elkies-k3"

TRANSPORT = LOCAL / "q24-d12-orbit42-fast-transport-qq.json"
MARKING = LOCAL / "q24-d12-orbit42-i8star-physical-marking-qq.json"
PREFLIGHT = LOCAL / "q24-d12-to-a11-orbit42-divval-preflight.json"

for path in (TRANSPORT, MARKING, PREFLIGHT):
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

transport = json.loads(TRANSPORT.read_text())
marking = json.loads(MARKING.read_text())
pre = json.loads(PREFLIGHT.read_text())

assert transport["status"] == "PASS_Q42_FAST_A11_CANDIDATES"
assert marking["status"] == "PASS_Q42_EXACT_I8STAR_PHYSICAL_MARKING"
assert pre["status"] == "PASS_Q42_DIVVAL_PREFLIGHT"

# Corrected orbit42 profile: hard anti-regression gates.
o42 = pre["orbit42"]
assert o42["mw_projection"] == [-1, 0, -1, -1, 0]
assert o42["height"] == "7"
assert o42["local_correction"] == "3"
assert int(o42["P_dot_O"]) == 3
assert int(o42["fibre_twist"]) == 0
assert o42["child_root_data"] == [11, 132, 12]
assert int(o42["child_mw_rank"]) == 6

orientations = marking["orientation_candidates"]
if len(orientations) != 2:
    raise ArithmeticError(
        f"expected exactly two physical spinor orientations, got {len(orientations)}"
    )
orientation = orientations[args.orientation]

omitted = list(orientation["omitted_vertical_components"])
meets = list(orientation["section_meets_physical_components"])

if len(omitted) != 1 or len(meets) != 1:
    raise ArithmeticError(
        f"orientation {args.orientation} does not have one omitted/met component"
    )
if omitted != meets:
    raise ArithmeticError(
        f"orientation {args.orientation}: omitted={omitted} != section_meets={meets}"
    )

expected_component = omitted[0]
vertical = {
    str(k): int(v)
    for k, v in orientation["vertical_coefficients_physical"].items()
}
if vertical.get(expected_component) != 0:
    raise ArithmeticError("omitted spinor component has nonzero V coefficient")
if sorted(vertical.values()) != [0] + [1]*11:
    raise ArithmeticError(
        "physical V is not exactly eleven coefficient-1 roots plus one omitted arm"
    )


# ---------------------------------------------------------------------------
# Exact section through the exact blow-up tree.
# ---------------------------------------------------------------------------
candidates = list(transport["passing_A11_candidates"])
if args.candidate is not None:
    candidates = [
        c for c in candidates
        if int(c["candidate_index"]) == int(args.candidate)
    ]
if not candidates:
    raise SystemExit("no passing A11 candidate selected")

R = PolynomialRing(QQ, "V")
V = R.gen()

alpha = QQ(marking["I8star"]["base_root"])
centres = list(marking["I8star"]["centres"])
by_label = {str(r["label"]): r for r in centres}
by_path = {str(r["path"]): r for r in centres}

root_rows = [r for r in centres if r["parent_label"] is None]
if len(root_rows) != 1:
    raise ArithmeticError("physical resolution has no unique root centre")
root = root_rows[0]


def child_record(parent, edge):
    path = str(edge["path"])
    if path not in by_path:
        raise ArithmeticError(f"missing child path {path}")
    return by_path[path]


LS = LaurentSeriesRing(QQ, "s", default_prec=int(args.precision))
s = LS.gen()


def eval_poly_series(poly, arg):
    out = LS(0)
    for c in reversed(list(poly)):
        out = out*arg + LS(QQ(c))
    return out


def series_constant(value):
    value = LS(value)
    if value.valuation() < 0:
        return None
    return QQ(value[0])


def equal_constant(value, target):
    c = series_constant(value)
    return c is not None and c == QQ(target)


def inverse_chart_series(coords, point, kind):
    """
    Lift a section through a point blow-up.

    Input coords are the old local coordinates as Laurent series in the
    original base parameter s.  Output is the selected affine blow-up chart.
    """
    us, xs, ys = map(LS, coords)
    a, b, c = map(QQ, point)

    du = us - a
    dx = xs - b
    dy = ys - c

    if kind == "u":
        if not du:
            raise ArithmeticError("section is identically in u=0 at blow-up")
        return (du, dx/du, dy/du)
    if kind == "x":
        if not dx:
            raise ArithmeticError("section is identically in x=0 at blow-up")
        return (du/dx, dx, dy/dx)
    if kind == "y":
        if not dy:
            raise ArithmeticError("section is identically in y=0 at blow-up")
        return (du/dy, dx/dy, dy)
    raise ValueError(kind)


def matches_point(coords, point):
    return all(
        equal_constant(value, target)
        for value, target in zip(coords, point)
    )


def section_component(XP, YP, ZP):
    base = LS(alpha) + s
    Zv = eval_poly_series(ZP, base)
    Xv = eval_poly_series(XP, base)
    Yv = eval_poly_series(YP, base)

    if not Zv:
        raise ArithmeticError(
            "marked P section has Z identically zero in I8* local series"
        )

    xsec = Xv / Zv**2
    ysec = Yv / Zv**3

    coords = (s, xsec, ysec)
    current = root
    visited = []

    while True:
        point = tuple(QQ(v) for v in current["point"])
        if not matches_point(coords, point):
            raise ArithmeticError(
                f"section does not pass expected blow-up centre {current['label']}"
            )

        visited.append(str(current["label"]))
        hits = []

        for edge in current.get("children", []):
            kind = str(edge["selected_chart"])
            try:
                lifted = inverse_chart_series(coords, point, kind)
            except (ArithmeticError, ZeroDivisionError):
                continue

            child = child_record(current, edge)
            child_point = tuple(QQ(v) for v in child["point"])
            if matches_point(lifted, child_point):
                hits.append((child, lifted, kind))

        if not hits:
            # The strict transform exits the unresolved singular chain here,
            # hence meets the exceptional component created at this centre.
            return str(current["label"]), visited, coords

        # A section is irreducible and must choose one infinitely-near centre.
        unique = {}
        for child, lifted, kind in hits:
            unique[str(child["label"])] = (child, lifted, kind)
        if len(unique) != 1:
            raise ArithmeticError(
                f"section lift is ambiguous at {current['label']}: "
                f"{sorted(unique)}"
            )

        current, coords, unused_kind = next(iter(unique.values()))


results = []

for cand in candidates:
    idx = int(cand["candidate_index"])
    XP = R([QQ(v) for v in cand["X"]])
    YP = R([QQ(v) for v in cand["Y"]])
    ZP = R([QQ(v) for v in cand["Z"]])

    assert ZP.degree() == 3
    assert XP.degree() <= 10
    assert YP.degree() <= 15

    actual_component, visited, final_coords = section_component(XP, YP, ZP)

    orientation_match = actual_component == expected_component

    print(
        "Q42RRQQ|"
        f"orientation={args.orientation}|candidate={idx}|"
        "ambient=SKIPPED_PRIMITIVE_NEF_ISOTROPIC|"
        "resolved=PHYSICAL_COMPONENT_MARKING|kernel=2|h0=2|"
        f"actual_component={actual_component}|"
        f"expected_component={expected_component}|"
        f"path={','.join(visited)}|"
        f"status={'PASS_ORIENTATION' if orientation_match else 'REJECT_ORIENTATION'}",
        flush=True,
    )

    if not orientation_match:
        continue

    child = cand["child"]
    if (
        int(child["root_rank"]) != 11
        or int(child["root_det"]) != 12
        or int(child["euler"]) != 24
        or int(child["MW_rank_if_rho19"]) != 6
    ):
        raise ArithmeticError(
            f"candidate {idx} passed orientation but is not exact A11/MW6"
        )

    if int(cand["quartic"]["degree"]) not in (3, 4):
        raise ArithmeticError("passing child quartic degree is not 3/4")

    results.append({
        "candidate_index": idx,
        "orientation": int(args.orientation),
        "actual_section_component": actual_component,
        "expected_section_component": expected_component,
        "visited_centres": visited,
        "omitted_equals_section_meets": True,
        "vertical_coefficients_physical": vertical,
        "rr": {
            "method": "PRIMITIVE_NEF_ISOTROPIC_K3_PLUS_DIRECT_CHORD",
            "ambient_matrix": "NOT_NEEDED",
            "kernel_dimension": 2,
            "h0": 2,
            "fibre_twist": 0,
        },
        "quartic_degree": int(cand["quartic"]["degree"]),
        "child_root_rank": int(child["root_rank"]),
        "child_root_det": int(child["root_det"]),
        "child_euler": int(child["euler"]),
        "child_MW_rank": int(child["MW_rank_if_rho19"]),
    })


OUT = (
    args.output.resolve()
    if args.output
    else LOCAL / (
        f"q24-d12-orbit42-fast-orientation-{args.orientation}-qq.json"
    )
)

status = (
    "PASS_EXACT_Q24_D12_Q6_A11_DIRECT_CHORD_COMPONENT_MARKING"
    if results
    else "Q42_FAST_ORIENTATION_REJECTED"
)

payload = {
    "schema": "elkies-k3.h3-q24-d12-orbit42-fast-orientation-qq.v1",
    "status": status,
    "orientation": int(args.orientation),
    "orientation_record": orientation,
    "results": results,
    "proof": {
        "corrected_profile": True,
        "primitive_nef_isotropic_D42": True,
        "K3_h0_D42": 2,
        "exact_D42_formula": "D42 = O + P + V",
        "fibre_twist": 0,
        "exact_section_lift_through_I8star": True,
        "omitted_equals_section_meets": bool(results),
        "direct_chord_child_A11": bool(results),
    },
}

OUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUT}", flush=True)
print(
    "Q42RRQQ_RESULT|"
    f"orientation={args.orientation}|passes={len(results)}|"
    f"status={status}",
    flush=True,
)

if not results:
    raise SystemExit(3)
