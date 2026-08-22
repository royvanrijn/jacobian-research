#!/usr/bin/env sage -python
"""Derive valuations in the formal resolved E7 normal-form model.

The formal III* model is ``Y^2=U^3+U*Z^3``. The component labels in its
Kodaira diagram do *not* determine the local orders of ``Z``, ``U``, and
``Y``: the first three exceptional curves are nonreduced in the intermediate
charts. This script records that normal-form chart data. It is not a
transported H92 atlas; see ``derive_h92_q6_actual_e7_resolution.sage``.

The atlas is an input to a later high-degree Riemann--Roch module computation.
It does not turn the third marked NS divisor into old Weierstrass coordinates,
and therefore makes no claim about a global kernel or a child equation.
"""

import argparse
import json
from pathlib import Path

from sage.all import ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-e7-valuation-atlas.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

# Each entry lists the pullback monomial exponents of (old Z, old U, old Y)
# in its displayed resolved chart.  ``chart_orders`` gives the orders of that
# chart's variables at the generic point of the named reduced exceptional
# curve.  A zero order denotes a unit there.  In particular, the values 2 for
# the first three Z/U exceptional parameters come from the chart equations
# Y^2=Z*(unit), not from a Kodaira label.
atlas = (
    ("E7_1", "E1, blow1_Z generic", (2, 0, 1), ((1, 0, 0), (1, 1, 0), (1, 0, 1))),
    ("E7_2", "E2, blow2_Z generic", (2, 0, 1), ((1, 0, 0), (2, 1, 0), (2, 0, 1))),
    ("E7_3", "E3, blow3_U generic", (0, 2, 1), ((1, 2, 0), (1, 3, 0), (1, 4, 1))),
    ("E7_4", "N3u, blow-up of blow3_U origin", (1, 1, 1), ((1, 2, 0), (1, 3, 0), (1, 4, 1))),
    ("E7_5", "N2, blow-up of blow2_Z origin", (1, 1, 1), ((1, 0, 0), (2, 1, 0), (2, 0, 1))),
    ("E7_6", "N3g, blow-up of blow3_U (Z,U,Y)=(-1,0,0)", (0, 1, 1), ((1, 2, 0), (1, 3, 0), (1, 4, 1))),
    ("E7_7", "N3z, blow-up of blow3_Z origin", (1, 1, 1), ((2, 1, 0), (3, 2, 0), (4, 2, 1))),
)


def order(chart_orders, monomial):
    return sum(ZZ(left) * ZZ(right) for left, right in zip(chart_orders, monomial))


entries = []
for label, chart, chart_orders, pullbacks in atlas:
    values = tuple(order(chart_orders, monomial) for monomial in pullbacks)
    z_order, u_order, y_order = values
    # The leading terms of the pulled-back defining equation must have a
    # repeated smallest order.  This is a local sanity check on the chart map,
    # rather than an inference from the E7 diagram.
    relation_orders = (2 * y_order, 3 * u_order, u_order + 3 * z_order)
    assert relation_orders.count(min(relation_orders)) >= 2
    entries.append({
        "component": label,
        "resolved_chart": chart,
        "chart_variable_orders": list(chart_orders),
        "old_coordinate_pullback_monomials": {
            name: list(monomial)
            for name, monomial in zip(("Z", "U", "Y"), pullbacks)
        },
        "old_coordinate_orders": {name: int(value) for name, value in zip(("Z", "U", "Y"), values)},
        "defining_equation_term_orders": {
            "Y^2": int(relation_orders[0]),
            "U^3": int(relation_orders[1]),
            "U*Z^3": int(relation_orders[2]),
        },
    })

# The orders of the old base coordinate Z are the III* fibre multiplicities.
base_multiplicities = vector(ZZ, [entry["old_coordinate_orders"]["Z"] for entry in entries])
cartan = matrix(ZZ, [
    [2, 0, 0, -1, 0, 0, 0],
    [0, 2, 0, 0, -1, 0, -1],
    [0, 0, 2, -1, 0, -1, -1],
    [-1, 0, -1, 2, 0, 0, 0],
    [0, -1, 0, 0, 2, 0, 0],
    [0, 0, -1, 0, 0, 2, 0],
    [0, -1, -1, 0, 0, 0, 2],
])
# The affine old-fibre component meets E7_1 once, hence the displayed vector
# satisfies C*m=(1,0,...,0).  This verifies the chart-derived multiplicities
# against the independently recorded resolution graph.
assert cartan * base_multiplicities == vector(ZZ, [1] + [0] * 6)
assert tuple(base_multiplicities) == (2, 2, 4, 3, 1, 2, 3)

payload = {
    "schema": "elkies-k3.h92-q6-e7-valuation-atlas.v1",
    "status": "PASS_FORMAL_E7_VALUATION_ATLAS",
    "equation": "Y^2=U^3+U*Z^3",
    "coordinates": "The old H92 local coordinates (Z,U,Y), with Z the old base parameter at the III* fibre.",
    "entries": entries,
    "old_base_fibre_multiplicities": [int(value) for value in base_multiplicities],
    "cartan_check": "E7 Cartan times multiplicities equals the affine attachment vector (1,0,0,0,0,0,0).",
    "boundary": "This is formal normal-form chart data, not a transported H92 chart. It is not a component-label substitution and cannot itself impose an H92 branch module or compute h0.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6E7ATLAS|components=7|Z_orders={}|status=PASS_FORMAL_E7_VALUATION_ATLAS".format(
        tuple(base_multiplicities)
    ),
    flush=True,
)
