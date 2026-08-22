#!/usr/bin/env sage -python
"""Assemble the endpoint-compatible q=8 coefficient ambient on H92.

For each generic q=8 basis element ``G=x^a*m^b``, two actual resolved
calculations provide complementary bounds for a coefficient
``u^i/h(u)^k``:

* the II* E8 module requires ``i >= e(G)``;
* the actual marked E7 frame allows a pole ``G/t^d`` and hence requires
  ``i <= 4*k+d`` at ``u=infinity``, since ``deg(h)=4``.

Choosing the least nonnegative ``k`` with ``e <= 4*k+d`` gives a finite,
canonical endpoint envelope.  Its basis is the set of all

    u^i/h(u)^k * G,       e <= i <= 4*k+d.

It is a controlled ambient for the remaining q=8 condition compiler, not
the claimed global Riemann--Roch space: smooth P1.O collision conditions and
the other resolved E7 chart transitions may shrink or require enlargement of
this seed before a complete two-dimensional kernel can be asserted.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
P1 = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
GENERIC = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-generic-rr-ambient.json"
E8 = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e8-ambient-weights.json"
E7 = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-actual-e7-marked-frame.json"
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-endpoint-rr-ambient.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--p1", type=Path, default=P1)
parser.add_argument("--generic", type=Path, default=GENERIC)
parser.add_argument("--e8", type=Path, default=E8)
parser.add_argument("--e7", type=Path, default=E7)
parser.add_argument(
    "--extra-h-power", type=int, default=0,
    help="raise each minimal endpoint denominator power by this nonnegative amount",
)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
if args.extra_h_power < 0:
    raise ValueError("extra-h-power must be nonnegative")

p1 = json.loads(args.p1.read_text())
generic = json.loads(args.generic.read_text())
e8 = json.loads(args.e8.read_text())
e7 = json.loads(args.e7.read_text())
exec(compile(CORE.read_text(), str(CORE), "exec"))
assert p1["status"] == "PASS_EXACT_H92_P1"
assert generic["status"] == "PASS_EXACT_Q8_GENERIC_RR_AMBIENT"
assert e8["status"] == "PASS_EXACT_Q8_E8_AMBIENT_WEIGHTS"
assert e7["status"] == "PASS_EXACT_Q8_ACTUAL_E7_MARKED_FRAME"

u_ring = PolynomialRing(QQ, "u")
h = u_ring([QQ(value) for value in p1["structured_denominator"]["Z4_coefficients"]])
assert h.degree() == 4 and h(0) == 1 and h.leading_coefficient()

def key(entry):
    return (entry["kind"], int(entry["x_power"]), int(entry["m_power"]))

generic_basis = {key(entry): entry for entry in generic["basis"]}
e8_floors = {
    key(entry["basis"]): int(entry["minimal_u_power"])
    for entry in e8["basis_weight_floors"]
}
e7_frames = {
    key(entry["basis"]): int(entry["t_denominator_power"])
    for entry in e7["normalized_q8_generators"]
}
assert set(generic_basis) == set(e8_floors) == set(e7_frames)

families = []
ambient_basis = []
for entry in generic["basis"]:
    label = key(entry)
    e8_floor = e8_floors[label]
    e7_pole = e7_frames[label]
    interval = endpoint_coefficient_interval(e8_floor, e7_pole, h.degree())
    h_power = interval["denominator_power"]+args.extra_h_power
    upper_u_power = int(h.degree()*h_power+e7_pole)
    assert e8_floor == interval["u_power_lower"] <= upper_u_power
    family = {
        "generic_basis": {"kind": label[0], "x_power": label[1], "m_power": label[2]},
        "e8_minimal_u_power": e8_floor,
        "e7_allowed_t_denominator_power": e7_pole,
        "h_power": h_power,
        "u_power_range": [int(e8_floor), upper_u_power],
    }
    families.append(family)
    for u_power in range(e8_floor, upper_u_power+1):
        ambient_basis.append({
            **family["generic_basis"],
            "u_power": int(u_power),
            "h_power": h_power,
            "coefficient": "u^{}/h(u)^{}".format(u_power, h_power),
        })

assert len(families) == 18
assert len(ambient_basis) == 54+72*args.extra_h_power
assert [family["h_power"] for family in families[:10]] == [2+args.extra_h_power, 2+args.extra_h_power, 3+args.extra_h_power, 3+args.extra_h_power, 4+args.extra_h_power, 4+args.extra_h_power, 5+args.extra_h_power, 5+args.extra_h_power, 6+args.extra_h_power, 6+args.extra_h_power]
assert [family["h_power"] for family in families[10:]] == [2+args.extra_h_power, 2+args.extra_h_power, 3+args.extra_h_power, 3+args.extra_h_power, 4+args.extra_h_power, 4+args.extra_h_power, 5+args.extra_h_power, 5+args.extra_h_power]

status = (
    "PASS_EXACT_Q8_ENDPOINT_RR_AMBIENT"
    if args.extra_h_power == 0 else "PASS_EXACT_Q8_ENLARGED_ENDPOINT_RR_AMBIENT"
)
boundary = (
    "This is the least-h-denominator ambient compatible with the actual "
    "E8 module and the actual E7 marked frame. It is not yet a complete "
    "q8 cover: impose smooth collision conditions and all remaining E7 "
    "chart/gluing conditions before computing a kernel or child equation."
    if args.extra_h_power == 0 else
    "This is an endpoint-compatible enlarged ambient, not yet a complete q8 "
    "cover: impose smooth collision conditions and all remaining E7 "
    "chart/gluing conditions before computing a kernel or child equation."
)

payload = {
    "schema": "elkies-k3.h92-q8-endpoint-rr-ambient.v1",
    "status": status,
    "inputs": {
        "p1": {"path": str(args.p1.relative_to(ROOT)), "sha256": digest(args.p1)},
        "compiler_core": {"path": str(CORE.relative_to(ROOT)), "sha256": digest(CORE)},
        "generic_ambient": {"path": str(args.generic.relative_to(ROOT)), "sha256": digest(args.generic)},
        "actual_e8_floors": {"path": str(args.e8.relative_to(ROOT)), "sha256": digest(args.e8)},
        "actual_e7_marked_frame": {"path": str(args.e7.relative_to(ROOT)), "sha256": digest(args.e7)},
    },
    "base_coordinate": "u=1/t",
    "collision_polynomial": {"h": str(h), "degree": int(h.degree()), "h_at_E8": str(h(0))},
    "coefficient_rule": "u^i/h(u)^k * x^a*m^b with e8_floor<=i<=4*k+e7_pole",
    "families": families,
    "ambient_basis": ambient_basis,
    "ambient_dimension": len(ambient_basis),
    "endpoint_proof": {
        "E8": "h(0) is a unit, so i>=e8_floor is exactly the actual E8 order condition.",
        "E7_marked": "At u=infinity, u^i/h^k has t-order 4*k-i; i<=4*k+e7_pole makes it at least t^(-e7_pole).",
    },
    "boundary": boundary,
}
if args.extra_h_power:
    payload["enlargement"] = {
        "extra_h_power": args.extra_h_power,
        "rule": "replace each least endpoint denominator power k by k+r while retaining e8_floor<=i<=4*(k+r)+e7_pole",
        "base_ambient_dimension": 54,
        "enlarged_ambient_dimension": len(ambient_basis),
        "boundary": "This enlargement relaxes no E8 floor or marked-E7 inequality; its remaining E7 edges still require resolved quotient conditions.",
    }
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
if args.extra_h_power == 0:
    print(
        "H92Q8ENDPOINTAMBIENT|families=18|basis=54|h_degree=4|"
        "status=PASS_EXACT_Q8_ENDPOINT_RR_AMBIENT",
        flush=True,
    )
else:
    print(
        "H92Q8ENDPOINTAMBIENT|families=18|basis={}|h_degree=4|extra_h={}|"
        "status={}".format(len(ambient_basis), args.extra_h_power, status),
        flush=True,
    )
