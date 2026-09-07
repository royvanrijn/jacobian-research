#!/usr/bin/env sage
"""Move the exact Q80 coefficient parameter so the CM24 point is t=0.

This small adapter lets the modular marked-jet fitter use the reconstructed
global coefficient line.  Its coordinates are centered exactly as in the
formal CM24 system, and ``u=u_CM24+t`` fixes the cross-prime marking of the
local parameter.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--parameter",
    type=Path,
    default=ROOT / "artifacts/generated-results/"
    "q80-cm24-slope-8-87-qq-PDQE-parameter.json",
)
parser.add_argument(
    "--family",
    type=Path,
    default=ROOT / "artifacts/generated-results/"
    "q80-cm24-slope-8-87-qq-surface-family.json",
)
parser.add_argument(
    "--output",
    type=Path,
    default=ROOT / "artifacts/generated-results/"
    "q80-cm24-slope-8-87-qq-local-parameter.json",
)
args = parser.parse_args()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


parameter_payload = json.loads(args.parameter.read_text())
family_payload = json.loads(args.family.read_text())
if parameter_payload.get("schema") != "q80-cm24-qq-PDQE-parameter-v1":
    raise ValueError("unexpected Q80 parameter schema")
if family_payload.get("schema") != "q80-qq-unmarked-surface-family-v1":
    raise ValueError("unexpected Q80 family schema")

source_ring = PolynomialRing(QQ, "u")
u = source_ring.gen()
local_ring = PolynomialRing(QQ, "t")
t = local_ring.gen()
u_cm24 = QQ(family_payload["cm24_anchor"]["parameter"])

centers = {
    "D": QQ(-1)/2,
    "P": QQ(9)/4,
    "Q": QQ(-9)/4,
    "E": QQ(-27)/32,
}
source_names = {"D": "d", "P": "p", "Q": "q", "E": "e"}
functions = {}
for centered_name, source_name in source_names.items():
    record = parameter_payload["original_functions"][source_name]
    numerator = source_ring(record["numerator"])
    denominator = source_ring(record["denominator"])
    local_numerator = local_ring(numerator(u_cm24+t))
    local_denominator = local_ring(denominator(u_cm24+t))
    centered_numerator = local_numerator-centers[centered_name]*local_denominator
    common = centered_numerator.gcd(local_denominator)
    centered_numerator //= common
    local_denominator //= common
    if local_denominator == 0 or centered_numerator.gcd(local_denominator) != 1:
        raise ArithmeticError(f"bad local function {centered_name}")
    if centered_numerator(0) != 0 or local_denominator(0) == 0:
        raise ArithmeticError(f"{centered_name} is not regular and centered at CM24")
    functions[centered_name] = {
        "value": f"({centered_numerator})/({local_denominator})",
        "numerator": str(centered_numerator),
        "denominator": str(local_denominator),
        "degrees": [int(centered_numerator.degree()), int(local_denominator.degree())],
    }

if functions["P"]["degrees"][0] < 1:
    raise ArithmeticError("P is not a local parameter at CM24")

output = {
    "schema": "q80-cm24-formal-branch-parameter-v1",
    "status": "PASS_EXACT_CM24_CENTERED_PARAMETER",
    "slope": "8/87",
    "global_to_local": f"u={u_cm24}+t",
    "cm24_global_parameter": str(u_cm24),
    "functions": functions,
    "inputs": {
        "parameter": {"path": str(args.parameter.relative_to(ROOT)), "sha256": sha256(args.parameter)},
        "family": {"path": str(args.family.relative_to(ROOT)), "sha256": sha256(args.family)},
    },
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
print(
    "Q80LOCALPARAMETER|coordinate=t|cm24=t0|P_valuation=1|"
    f"output={args.output}|status=PASS_EXACT_CM24_CENTERED_PARAMETER",
    flush=True,
)
