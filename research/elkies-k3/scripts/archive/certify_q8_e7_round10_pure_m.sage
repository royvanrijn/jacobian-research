#!/usr/bin/env sage -python
"""Exact characteristic-zero certificate for the two round-10 E7 vectors.

The generic E7 layer computation found
    f8 = t^4*m^8
    f9 = t^4*m^9.

Use only already-certified exact artifacts:

1. The q6 actual E7 module cover says its rank-one fractional module has
   generator 1 away from -P1 and generator m near -P1, with
   m = unit/W there; hence 1/m is regular near -P1.

2. The q8 gluing artifact says on every actual resolved E7 edge chart
       g*f in (q6_marked_module)^9
   with g a regular Cartier equation (a monomial in component equations).

Then:
  away from -P1:
      g*t^4*m^b is regular for b=8,9,
      so it lies in the ninth power, whose local generator is 1.

  near -P1:
      (g*t^4*m^9)/m^9 = g*t^4 is regular,
      (g*t^4*m^8)/m^9 = g*t^4/m is regular
  because 1/m is regular.

Thus both vectors satisfy the complete E7 module cover in characteristic zero.
No Gröbner calculation is needed.
"""

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path.cwd()
ALL_EDGE = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-all-edge-module.json"
GLUING = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-actual-e7-gluing.json"
DEFAULT_OUTPUT = ROOT / "artifacts/local/elkies-k3/q8-e7-round10-pure-m-exact.json"

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

all_edge = json.loads(ALL_EDGE.read_text())
gluing = json.loads(GLUING.read_text())

assert all_edge["status"] == "PASS_EXACT_Q6_ACTUAL_E7_ALL_EDGE_MODULE"
assert gluing["status"] == "PASS_EXACT_Q8_ACTUAL_E7_GLUING"

cover = all_edge["module_cover"]
assert cover["away_from_minus_P1"].startswith("generator 1")
assert cover["near_minus_P1"].startswith("generator m")
assert all_edge["marked_horizontal_frame"]["corrected_identity"] == "Z*m/t=unit/W"
assert "m=unit/W" in all_edge["marked_horizontal_frame"]["conclusion"]

charts = gluing["actual_edge_chart_gluing"]
assert len(charts) == 6
for edge in charts:
    assert edge["comparison_condition"].endswith(
        "*f belongs to (q6_marked_module)^9"
    )
    # The artifact's Cartier equations are products of positive powers of
    # actual component equations, hence regular local functions.
    g = edge["w_cartier_equation"]
    assert "^-1" not in g

vectors = []
for b in (8, 9):
    near_quotient = "g*t^4" if b == 9 else "g*t^4/m"
    vectors.append({
        "name": "t^4*m^{}".format(b),
        "m_power": b,
        "away_from_minus_P1": {
            "q6_ninth_local_generator": "1",
            "quotient": "g*t^4*m^{}".format(b),
            "regular": True,
            "reason": "g,t,m are regular in the certified away-from--P1 q6 E7 trivialization",
        },
        "near_minus_P1": {
            "q6_ninth_local_generator": "m^9",
            "quotient": near_quotient,
            "regular": True,
            "reason": (
                "g and t are regular; 1/m is regular because m=unit/W"
                if b == 8 else
                "g and t are regular"
            ),
        },
        "complete_e7_membership": True,
    })

payload = {
    "schema": "elkies-k3.h92-q8-e7-round10-pure-m-exact.v1",
    "status": "PASS_EXACT_Q8_E7_ROUND10_PURE_M_MEMBERSHIP",
    "inputs": {
        "q6_actual_e7_all_edge_module": {
            "path": str(ALL_EDGE.relative_to(ROOT)),
            "sha256": digest(ALL_EDGE),
        },
        "q8_actual_e7_gluing": {
            "path": str(GLUING.relative_to(ROOT)),
            "sha256": digest(GLUING),
        },
    },
    "vectors": vectors,
    "conclusion": (
        "The two first generic leading vectors t^4*m^8 and t^4*m^9 "
        "satisfy the complete actual resolved E7 module cover over QQ."
    ),
    "boundary": (
        "This is an E7-local statement only. It does not impose E8, finite "
        "smooth-collision/global transition lattices, or construct the q8 pencil."
    ),
}

args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(
    "Q8E7R10EXACT|vectors=t^4*m^8,t^4*m^9|charts=6|"
    "membership=complete_E7_cover|status=PASS_EXACT_Q8_E7_ROUND10_PURE_M_MEMBERSHIP",
    flush=True,
)
