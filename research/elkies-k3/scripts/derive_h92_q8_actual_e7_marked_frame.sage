#!/usr/bin/env sage -python
"""Normalize the q=8 generic basis at the actual marked E7 chart.

On the actual ``E7_2--E7_5`` chart, the corrected q=6 marked frame proves

    n=Z*m/t = unit/W,

near ``-P1``, where ``W`` is the transverse marked-point parameter.  The
q=8 E7 twist has Cartier factor ``g=Z^6*Y^5`` there.  Since ``t/Z``,
``Y``, and ``x/t^2`` are units at this point, the following generators obey

    g*(m^b/t^6) = unit*n^b,
    g*(x*m^b/t^8) = unit*n^b.

They therefore lie in the ninth q=6 marked module after multiplication by
the actual q=8 gluing factor, and have marked pole orders b (at most nine).
The result is local to this actual chart: it gives the E7 marked-point side
of a q=8 coefficient ambient, not compatible global functions.
"""

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
AMBIENT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-generic-rr-ambient.json"
MARKED = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-p1-actual-e7-marked-module-corrected.json"
GLUING = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-actual-e7-gluing.json"
TRACE = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-p1-actual-e7-trace.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-actual-e7-marked-frame.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--ambient", type=Path, default=AMBIENT)
parser.add_argument("--marked", type=Path, default=MARKED)
parser.add_argument("--gluing", type=Path, default=GLUING)
parser.add_argument("--trace", type=Path, default=TRACE)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

ambient = json.loads(args.ambient.read_text())
marked = json.loads(args.marked.read_text())
gluing = json.loads(args.gluing.read_text())
trace = json.loads(args.trace.read_text())
assert ambient["status"] == "PASS_EXACT_Q8_GENERIC_RR_AMBIENT"
assert marked["status"] == "PASS_EXACT_P1_ACTUAL_E7_MARKED_MODULE_CORRECTED"
assert gluing["status"] == "PASS_EXACT_Q8_ACTUAL_E7_GLUING"
assert trace["status"] == "PASS_EXACT_P1_ACTUAL_E7_TRACE"
assert marked["chart"] == "E7_2--E7_5"
assert marked["marked_component"] == "E7_5"
assert marked["unit_coefficients_at_minus_P1"]["(Z*m/t)*W"] != "0"
assert trace["node_chart"]["-P1"]["U"] != "0"
assert trace["node_chart"]["-P1"]["Y"] != "0"

edge = next(
    entry for entry in gluing["actual_edge_chart_gluing"]
    if entry["name"] == "E7_2--E7_5"
)
assert edge["w_cartier_equation"] == "Z^6*Y^5"
assert edge["comparison_condition"] == "(Z^6*Y^5)*f belongs to (q6_marked_module)^9"

# The exact actual-chart pullbacks from the corrected q6 marked-frame
# certificate show t/Z, Y, and x/t^2 are units at -P1.  Thus g/t^6 is a
# unit, and m=(t/Z)*n is a unit multiple of the corrected generator n.
pullbacks = marked["local_coordinates"]["old_coordinate_pullbacks"]
assert pullbacks["t"] == "Z*U*(Z-A1/B1)"
assert pullbacks["x"] == "Z^2*U^2*(Z-A1/B1)"
unit_data = marked["unit_coefficients_at_minus_P1"]
assert unit_data["t_over_Z"] != "0"
assert unit_data["(Z*m/t)*W"] != "0"

generators = []
for entry in ambient["basis"]:
    x_power = int(entry["x_power"])
    m_power = int(entry["m_power"])
    denominator_power = 6+2*x_power
    assert x_power in (0, 1)
    assert m_power <= 9
    if x_power == 0:
        numerator = "m^{}".format(m_power)
    else:
        numerator = "x*m^{}".format(m_power)
    generators.append({
        "basis": {"kind": entry["kind"], "x_power": x_power, "m_power": m_power},
        "normalized_generator": "{}/t^{}".format(numerator, denominator_power),
        "t_denominator_power": denominator_power,
        "q6_ninth_power_comparison": (
            "(Z^6*Y^5)*({}/t^{}) = unit*(Z*m/t)^{}".format(
                numerator, denominator_power, m_power
            )
        ),
        "marked_pole_order_at_minus_P1": m_power,
    })

assert [row["normalized_generator"] for row in generators[:10]] == [
    "m^{}/t^6".format(power) for power in range(10)
]
assert [row["normalized_generator"] for row in generators[10:]] == [
    "x*m^{}/t^8".format(power) for power in range(8)
]
assert max(row["marked_pole_order_at_minus_P1"] for row in generators) == 9

payload = {
    "schema": "elkies-k3.h92-q8-actual-e7-marked-frame.v1",
    "status": "PASS_EXACT_Q8_ACTUAL_E7_MARKED_FRAME",
    "inputs": {
        "generic_ambient": {"path": str(args.ambient.relative_to(ROOT)), "sha256": digest(args.ambient)},
        "q6_actual_marked_frame": {"path": str(args.marked.relative_to(ROOT)), "sha256": digest(args.marked)},
        "q8_actual_e7_gluing": {"path": str(args.gluing.relative_to(ROOT)), "sha256": digest(args.gluing)},
        "actual_p1_trace": {"path": str(args.trace.relative_to(ROOT)), "sha256": digest(args.trace)},
    },
    "actual_chart": {
        "name": "E7_2--E7_5",
        "marked_component": "E7_5",
        "transverse_parameter": marked["local_coordinates"]["transverse_at_minus_P1"],
        "q6_frame": "Z*m/t=unit/W; m=(t/Z)*(Z*m/t)",
        "q8_gluing_factor": "g=Z^6*Y^5",
        "unit_facts": ["t/Z", "Y", "x/t^2", "(Z*m/t)*W"],
    },
    "normalized_q8_generators": generators,
    "conclusion": (
        "Every listed q8 generic-basis element has an exact marked-chart "
        "normalization whose multiplication by g lies in the ninth q6 marked "
        "module. Its pole order at -P1 is the stated m-power."
    ),
    "boundary": (
        "This is only the marked smooth-point and one edge-chart portion of "
        "the q8 E7 condition. It does not prove compatibility at the other "
        "E7 edges, a finite quotient, global base bounds, h0=2, or a pencil."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q8ACTUALE7MARKED|basis=18|m_denominator=6|xm_denominator=8|"
    "status=PASS_EXACT_Q8_ACTUAL_E7_MARKED_FRAME",
    flush=True,
)
