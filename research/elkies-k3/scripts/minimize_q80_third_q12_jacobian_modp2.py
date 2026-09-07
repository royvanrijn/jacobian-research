#!/usr/bin/env python3
"""Adapt the exact third-q12 global minimizer to a common-producer prime."""

import argparse
import contextlib
import hashlib
import io
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CORE = ROOT / "elkies-k3/scripts/minimize_q80_third_q12_jacobian_mod19_quadratic.sage"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
args.input = args.input.resolve()
args.output = args.output.resolve()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


jacobian = json.loads(args.input.read_text())
expected_input = "PASS_EXACT_INTERPOLATED_THIRD_Q12_JACOBIAN_COMMON_PRODUCER"
if jacobian.get("status") != expected_input:
    raise ValueError("input is not a certified common-producer Jacobian")
specialization = jacobian["specialization"]
prime = int(specialization["prime"])
match = re.fullmatch(r"r\^2 \+ (\d+)\*r \+ (\d+)", specialization["extension_modulus"])
if match is None:
    raise ValueError("cannot parse quadratic extension modulus")
linear, constant = map(int, match.groups())

source = CORE.read_text()
old_paths = (
    'INPUT = ROOT / "artifacts/generated-results/q80-third-q12-p19-jacobian-interpolated.json"\n'
    'OUTPUT = ROOT / "artifacts/generated-results/q80-third-q12-p19-jacobian-minimal.json"'
)
new_paths = f"INPUT = Path({str(args.input)!r})\nOUTPUT = Path({str(args.output)!r})"
if source.count(old_paths) != 1:
    raise ArithmeticError("immutable minimizer path contract changed")
source = source.replace(old_paths, new_paths, 1)
source = source.replace(
    'EXPECTED = "PASS_EXACT_INTERPOLATED_THIRD_Q12_JACOBIAN_MOD19_QUADRATIC"',
    f'EXPECTED = "{expected_input}"', 1,
)
old_field = (
    'base_finite = GF(19)\nmodulus_ring = PolynomialRing(base_finite, "m")\n'
    'm = modulus_ring.gen()\nfinite = GF(19**2, "r", modulus=m**2 + 12 * m + 3)'
)
new_field = (
    f'base_finite = GF({prime})\nmodulus_ring = PolynomialRing(base_finite, "m")\n'
    f'm = modulus_ring.gen()\nfinite = GF({prime}**2, "r", modulus=m**2 + {linear} * m + {constant})'
)
if source.count(old_field) != 1:
    raise ArithmeticError("immutable minimizer field contract changed")
source = source.replace(old_field, new_field, 1)
status = "PASS_EXACT_MINIMAL_THIRD_Q12_JACOBIAN_AND_FIBRES_COMMON_PRODUCER"
source = source.replace(
    '"PASS_EXACT_MINIMAL_THIRD_Q12_JACOBIAN_AND_FIBRES_MOD19_QUADRATIC"',
    f'"{status}"',
)
source = source.replace(
    '"specialization": {"u": "-2", "prime": 19, "extension_modulus": "r^2+12*r+3"}',
    '"specialization": common_specialization', 1,
)

namespace = {
    "__file__": str(CORE),
    "__name__": "__main__",
    "common_specialization": specialization,
}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(source, str(CORE), "exec"), namespace)

output = json.loads(args.output.read_text())
if output.get("status") != status:
    raise ArithmeticError("adapted global minimizer did not pass")
output["schema"] = "elkies-k3.q80-third-q12-jacobian-minimal-modp2.v2"
output["status"] = status
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
        f"explicit integral short Weierstrass K3 model over GF({prime}^2)(V)",
        "exact long/short scaling and inverse",
        "minimal multiplicative fibre configuration I6+I4+3I2+8I1",
        "root lattice A5+A3+3A1",
    ],
    "not_proved": [
        "transported old-component and zero-section marking at this prime",
        "composition of the generic maps with this minimal gauge",
        "a characteristic-zero coefficient reconstruction",
    ],
}
output["reproduce"] = (
    "sage -python elkies-k3/scripts/minimize_q80_third_q12_jacobian_modp2.py "
    f"--input {args.input} --output {args.output}"
)
args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
print(
    f"Q80THIRDQ12COMMONMINIMAL|prime={prime}|degrees=8,12|Delta_degree=24|"
    "fibres=I6+I4+3I2+8I1|roots=A5+A3+3A1|"
    "status=PASS_EXACT_MINIMAL_THIRD_Q12_JACOBIAN_AND_FIBRES_COMMON_PRODUCER"
)
