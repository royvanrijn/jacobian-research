#!/usr/bin/env sage -python
"""Record actual H92 (t,x,y) pullbacks on the six resolved E7 edge charts."""

import argparse
import hashlib
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import PolynomialRing, QQ

ROOT = Path(__file__).resolve().parents[2]
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
RESOLUTION = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-resolution-full.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-chart-pullbacks.json"
RESOLUTION_SHA256 = "14378f4718d3fbe781d5b351ba4943a962a430d9827c0ce285cd1125a9e8c500"

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def strict(ring, value, substitutions, exceptional):
    transformed = ring(value(*substitutions))
    quotient, remainder = transformed.quo_rem(exceptional**2)
    assert not remainder
    return quotient

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--resolution", type=Path, default=RESOLUTION)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
if args.resolution == RESOLUTION:
    assert digest(args.resolution) == RESOLUTION_SHA256
resolution = json.loads(args.resolution.read_text())
assert resolution["status"] == "PASS_EXACT_H92_E7_BLOWUP_TREE"

anchor = SourceFileLoader("h92_actual_e7_pullback_anchor", str(ANCHOR)).load_module()
r, s = anchor.EXPECTED_H92
_, formulas = anchor.parse_h92(H92)
A1, A, B1, B, B2 = (QQ(value(r, s)) for value in formulas)

ring = PolynomialRing(QQ, names=("Z", "U", "Y"))
Z, U, Y = ring.gens()
f0 = Y**2 - U**3 - (A1 * Z**3 + A * Z**4) * U - (B1 * Z**5 + B * Z**6 + B2 * Z**7)
f1_z = strict(ring, f0, (Z, Z * U, Z * Y), Z)
f2_u = strict(ring, f1_z, (U * Z, U, U * Y), U)
f3_u = strict(ring, f2_u, (U * Z, U, U * Y), U)
f3_z = strict(ring, f2_u, (Z, Z * U, Z * Y), Z)

second = -A1 / B1
third = -QQ(1) / A1
edge_equations = {
    "E7_2--E7_5": strict(ring, ring(f2_u(Z + second, U, Y)), (Z, Z * U, Z * Y), Z),
    "E7_1--E7_4": strict(ring, f3_u, (U * Z, U, U * Y), U),
    "E7_4--E7_3": strict(ring, f3_u, (Z, Z * U, Z * Y), Z),
    "E7_3--E7_7": strict(ring, f3_z, (U * Z, U, U * Y), U),
    "E7_7--E7_2": strict(ring, f3_z, (Z, Z * U, Z * Y), Z),
    "E7_3--E7_6": strict(ring, ring(f3_u(Z + third, U, Y)), (Z, Z * U, Z * Y), Z),
}
assert {name: str(value) for name, value in edge_equations.items()} == resolution["edge_charts"]

maps = {
    "E7_2--E7_5": ((Z * U * (Z + second), Z**2 * U**2 * (Z + second), Z**3 * U**2 * (Z + second) * Y), Z**6 * U**4 * (Z + second)**2),
    "E7_1--E7_4": ((U**3 * Z, U**4 * Z, U**6 * Z * Y), U**12 * Z**2),
    "E7_4--E7_3": ((Z**3 * U**2, Z**4 * U**3, Z**6 * U**4 * Y), Z**12 * U**8),
    "E7_3--E7_7": ((U**3 * Z**2, U**5 * Z**3, U**7 * Z**4 * Y), U**14 * Z**8),
    "E7_7--E7_2": ((Z**3 * U, Z**5 * U**2, Z**7 * U**2 * Y), Z**14 * U**4),
    "E7_3--E7_6": ((Z**2 * U**2 * (Z + third), Z**3 * U**3 * (Z + third), Z**5 * U**4 * (Z + third) * Y), Z**10 * U**8 * (Z + third)**2),
}

records = []
for name, (pullback, total_factor) in maps.items():
    if f0(*pullback) != total_factor * edge_equations[name]:
        raise ValueError("pullback factor mismatch in {}".format(name))
    records.append({
        "name": name,
        "surface_equation": str(edge_equations[name]),
        "old_coordinate_pullback": {key: str(value) for key, value in zip(("t", "x", "y"), pullback)},
        "old_coordinate_values_at_chart_origin": {key: str(value(0, 0, 0)) for key, value in zip(("t", "x", "y"), pullback)},
        "total_transform_factor": str(total_factor),
    })

payload = {
    "schema": "elkies-k3.h92-q6-actual-e7-chart-pullbacks.v1",
    "status": "PASS_EXACT_H92_E7_CHART_PULLBACKS",
    "inputs": {"actual_resolution": {"path": str(args.resolution.relative_to(ROOT)), "sha256": digest(args.resolution)}},
    "node_translations": {"second_U_nonzero": str(second), "third_U_generic": str(third)},
    "charts": records,
    "compiler_instruction": "Substitute these maps before evaluating an old-model basis element or a marked-point DAG in the local chart ring.",
    "boundary": "This supplies pullbacks, not the marked-chord evaluation or finite quotient conditions.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("H92ACTUALE7PULLBACK|charts=6|status=PASS_EXACT_H92_E7_CHART_PULLBACKS", flush=True)
