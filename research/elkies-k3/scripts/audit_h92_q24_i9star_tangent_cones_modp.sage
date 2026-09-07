#!/usr/bin/env sage -python
"""
Audit the exceptional divisors created by the current q24 I9* blow-up tree.

The existing resolver labels one blow-up center as one E## object.  That is
not generally the same thing as one irreducible exceptional curve: the
projectivized tangent cone of a nonordinary center may be reducible.

For every center and every affine blow-up chart this script restricts the
stored strict transform to the new exceptional coordinate and factors that
exceptional slice.  `max_reduced_factors` is a practical lower bound / usually
the exact number of reduced irreducible components created by that blow-up.

If 12 blow-up centers give a total of 13 reduced tangent-cone components, the
current "expected 13 centers" criterion is simply the wrong invariant.
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, ZZ


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
    ]
    for candidate in candidates:
        try:
            candidate = candidate.resolve()
        except Exception:
            continue
        if (candidate / "elkies-k3/scripts").is_dir():
            return candidate
    raise SystemExit("Could not locate jacobian-research")


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--prime", type=int, default=100003)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
path = ROOT / "artifacts/local/elkies-k3" / (
    f"q24-i9star-resolution-mod-{args.prime}.json"
)
if not path.exists():
    raise SystemExit(f"Missing {path}")

data = json.loads(path.read_text())
p = ZZ(args.prime)
F = GF(p)
S = PolynomialRing(F, names=("u", "x", "y"), order="degrevlex")
u, x, y = S.gens()
gens = {"u": u, "x": x, "y": y}

center_records = []
component_lower_bound = 0
split_centers = []

for center in data["centers"]:
    per_chart = []
    max_reduced = 0
    for chart in center.get("charts", []):
        kind = chart["chart"]
        exceptional = gens[kind]
        strict = S(chart["strict_transform"])
        slice_poly = S(strict.subs({exceptional: F(0)}))

        if not slice_poly:
            factor_record = "ZERO"
            distinct = 0
            multiplicities = []
        elif slice_poly.total_degree() == 0:
            factor_record = str(slice_poly)
            distinct = 0
            multiplicities = []
        else:
            fac = slice_poly.factor()
            factors = []
            multiplicities = []
            for factor, exponent in fac:
                factors.append(str(factor))
                multiplicities.append(int(exponent))
            distinct = len(factors)
            factor_record = "*".join(
                f"({factor})^{exponent}"
                for factor, exponent in fac
            ) or str(fac.unit())

        max_reduced = max(max_reduced, distinct)
        per_chart.append(
            {
                "chart": kind,
                "slice": str(slice_poly),
                "factorization": factor_record,
                "distinct_nonconstant_factors": distinct,
                "multiplicities": multiplicities,
            }
        )

        print(
            "Q24TANGENT_CHART|"
            f"center={center['label']}|chart={kind}|"
            f"slice={slice_poly}|factors={factor_record}|"
            f"distinct={distinct}|status=PASS",
            flush=True,
        )

    # A nonzero irreducible/repeated tangent-cone component may disappear in
    # one affine chart, so take the maximum seen among u/x/y charts.
    if max_reduced == 0:
        max_reduced = 1
    component_lower_bound += max_reduced
    if max_reduced > 1:
        split_centers.append(center["label"])

    center_records.append(
        {
            "label": center["label"],
            "ordinary_double_point": bool(center["ordinary_double_point"]),
            "multiplicity": int(center["multiplicity"]),
            "active_components": center.get("active_components", []),
            "max_reduced_factors": max_reduced,
            "charts": per_chart,
        }
    )

    print(
        "Q24TANGENT_CENTER|"
        f"center={center['label']}|"
        f"ordinary={int(center['ordinary_double_point'])}|"
        f"mult={center['multiplicity']}|"
        f"max_reduced_factors={max_reduced}|"
        f"status={'SPLIT' if max_reduced > 1 else 'ONE_REDUCED_COMPONENT'}",
        flush=True,
    )

print(
    "Q24TANGENT_RESULT|"
    f"centers={len(data['centers'])}|"
    f"reduced_component_lower_bound={component_lower_bound}|"
    f"split_centers={','.join(split_centers) or 'NONE'}|"
    f"expected_D13=13|"
    f"status={'EXPLAINS_D13' if component_lower_bound == 13 else 'NEEDS_MORE_GEOMETRY'}",
    flush=True,
)

out = ROOT / "artifacts/local/elkies-k3" / (
    f"q24-i9star-tangent-cone-audit-mod-{p}.json"
)
out.write_text(
    json.dumps(
        {
            "schema": "elkies-k3.h3-q24-i9star-tangent-cone-audit.v1",
            "source_status": data.get("status"),
            "prime": int(p),
            "center_count": len(data["centers"]),
            "reduced_component_lower_bound": component_lower_bound,
            "split_centers": split_centers,
            "centers": center_records,
        },
        indent=2,
        sort_keys=True,
    )
    + "\n"
)
print(f"OUTPUT|{out}", flush=True)
