#!/usr/bin/env python3
"""Adapt the exact generic-Jacobian interpolator to a common-producer prime."""

import argparse
import contextlib
import hashlib
import io
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CORE = ROOT / "elkies-k3/scripts/interpolate_q80_third_q12_jacobian_mod19_quadratic.sage"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
args.input = args.input.resolve()
args.output = args.output.resolve()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


batch = json.loads(args.input.read_text())
if batch.get("status") != "PASS_EXACT_THIRD_Q12_WEIERSTRASS_SAMPLE_BATCH_COMMON_PRODUCER":
    raise ValueError("input is not a certified common-producer mapped batch")
specialization = batch["specialization"]
prime = int(specialization["prime"])
match = re.fullmatch(r"r\^2 \+ (\d+)\*r \+ (\d+)", specialization["extension_modulus"])
if match is None:
    raise ValueError("cannot parse quadratic extension modulus")
linear, constant = map(int, match.groups())

source = CORE.read_text()
old_paths = (
    'INPUT = ROOT / "artifacts/generated-results/q80-third-q12-p19-weierstrass-sample-batch.json"\n'
    'OUTPUT = ROOT / "artifacts/generated-results/q80-third-q12-p19-jacobian-interpolated.json"'
)
new_paths = f"INPUT = Path({str(args.input)!r})\nOUTPUT = Path({str(args.output)!r})"
if source.count(old_paths) != 1:
    raise ArithmeticError("immutable Jacobian core path contract changed")
source = source.replace(old_paths, new_paths, 1)
source = source.replace(
    'EXPECTED = "PASS_EXACT_THIRD_Q12_WEIERSTRASS_SAMPLE_BATCH_MOD19_QUADRATIC"',
    'EXPECTED = "PASS_EXACT_THIRD_Q12_WEIERSTRASS_SAMPLE_BATCH_COMMON_PRODUCER"',
    1,
)
old_field = (
    'base_finite = GF(19)\n'
    'modulus_ring = PolynomialRing(base_finite, "m")\n'
    'm = modulus_ring.gen()\n'
    'finite = GF(19**2, "r", modulus=m**2 + 12 * m + 3)'
)
new_field = (
    f'base_finite = GF({prime})\n'
    'modulus_ring = PolynomialRing(base_finite, "m")\n'
    'm = modulus_ring.gen()\n'
    f'finite = GF({prime}**2, "r", modulus=m**2 + {linear} * m + {constant})'
)
if source.count(old_field) != 1:
    raise ArithmeticError("immutable Jacobian core field contract changed")
source = source.replace(old_field, new_field, 1)

saved_name = "PASS_EXACT_INTERPOLATED_THIRD_Q12_JACOBIAN_MOD19_QUADRATIC"
temporary_name = "PASS_EXACT_INTERPOLATED_THIRD_Q12_JACOBIAN_COMMON_PRODUCER"
source = source.replace(saved_name, temporary_name)
source = source.replace(
    '"specialization": {"u": "-2", "prime": 19, "extension_modulus": "r^2+12*r+3"}',
    '"specialization": batch_specialization',
    1,
)
source = source.replace(
    '"reproduce": "sage -python elkies-k3/scripts/interpolate_q80_third_q12_jacobian_mod19_quadratic.sage"',
    '"reproduce": "adapter fills this field"',
    1,
)

namespace = {
    "__file__": str(CORE),
    "__name__": "__main__",
    "batch_specialization": specialization,
}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(source, str(CORE), "exec"), namespace)

output = json.loads(args.output.read_text())
if output.get("status") != temporary_name:
    raise ArithmeticError("adapted generic Jacobian did not pass")
output["schema"] = "elkies-k3.q80-third-q12-jacobian-interpolated-modp2.v2"
output["status"] = temporary_name
output["specialization"] = specialization
output["input"] = {"path": str(args.input.relative_to(ROOT)), "sha256": sha256(args.input)}
output["worker"] = {
    "core": {"path": str(CORE.relative_to(ROOT)), "sha256": sha256(CORE)},
    "adapter": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
}
output["claim_boundary"] = {
    "proved": [
        f"generic long Weierstrass Jacobian over GF({prime}^2)(V) in the pinned Laurent gauge",
        "exact discriminant and j identities",
        f"{batch['training_count']} training and {batch['held_out_count']} held-out mapped-fibre replays",
    ],
    "not_proved": [
        "generic interpolation of the forward/inverse maps",
        "global minimality and A5+A3+3A1 fibre marking",
        "a characteristic-zero coefficient reconstruction",
    ],
}
output["reproduce"] = (
    "sage -python elkies-k3/scripts/interpolate_q80_third_q12_jacobian_modp2.py "
    f"--input {args.input} --output {args.output}"
)
args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
degrees = tuple(
    output["weierstrass"][name]["degrees_numerator_denominator"]
    for name in ("a1", "a2", "a3", "a4", "a6")
)
print(
    f"Q80THIRDQ12COMMONJACOBIAN|prime={prime}|training={batch['training_count']}|"
    f"heldout={batch['held_out_count']}|degrees={degrees}|"
    "status=PASS_EXACT_INTERPOLATED_THIRD_Q12_JACOBIAN_COMMON_PRODUCER"
)
