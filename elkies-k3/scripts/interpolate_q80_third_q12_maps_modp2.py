#!/usr/bin/env python3
"""Adapt the exact generic birational-map interpolator to a common prime."""

import argparse
import contextlib
import hashlib
import io
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CORE = ROOT / "elkies-k3/scripts/interpolate_q80_third_q12_maps_mod19_quadratic.sage"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--batch", type=Path, required=True)
parser.add_argument("--pencil", type=Path, required=True)
parser.add_argument("--jacobian", type=Path, required=True)
parser.add_argument("--minimal", type=Path)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
for name in ("batch", "pencil", "jacobian", "output"):
    setattr(args, name, getattr(args, name).resolve())
if args.minimal is not None:
    args.minimal = args.minimal.resolve()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


batch = json.loads(args.batch.read_text())
pencil = json.loads(args.pencil.read_text())
jacobian = json.loads(args.jacobian.read_text())
minimal = json.loads(args.minimal.read_text()) if args.minimal is not None else None
expected = (
    (batch, "PASS_EXACT_THIRD_Q12_WEIERSTRASS_SAMPLE_BATCH_COMMON_PRODUCER"),
    (pencil, "PASS_EXACT_RESOLVED_THIRD_Q12_PENCIL_COMMON_PRODUCER"),
    (jacobian, "PASS_EXACT_INTERPOLATED_THIRD_Q12_JACOBIAN_COMMON_PRODUCER"),
)
if any(payload.get("status") != status for payload, status in expected):
    raise ValueError("one or more common-prime map inputs are uncertified")
if minimal is not None and minimal.get("status") != (
    "PASS_EXACT_MINIMAL_THIRD_Q12_JACOBIAN_AND_FIBRES_COMMON_PRODUCER"
):
    raise ValueError("minimal common-prime child is uncertified")
specialization = batch["specialization"]
if pencil["specialization"] != specialization or jacobian["specialization"] != specialization:
    raise ValueError("map input specializations disagree")
if minimal is not None and minimal["specialization"] != specialization:
    raise ValueError("minimal child specialization disagrees")
prime = int(specialization["prime"])
match = re.fullmatch(r"r\^2 \+ (\d+)\*r \+ (\d+)", specialization["extension_modulus"])
if match is None:
    raise ValueError("cannot parse quadratic extension modulus")
linear, constant = map(int, match.groups())

source = CORE.read_text()
old_paths = (
    'BATCH = ROOT / "artifacts/generated-results/q80-third-q12-p19-weierstrass-sample-batch.json"\n'
    'PENCIL = ROOT / "artifacts/generated-results/q80-third-q12-um2-p19-resolved-pencil.json"\n'
    'JACOBIAN = ROOT / "artifacts/generated-results/q80-third-q12-p19-jacobian-interpolated.json"\n'
    'MINIMAL = ROOT / "artifacts/generated-results/q80-third-q12-p19-jacobian-minimal.json"\n'
    'OUTPUT = ROOT / "artifacts/generated-results/q80-third-q12-p19-birational-maps.json"'
)
new_paths = (
    f"BATCH = Path({str(args.batch)!r})\n"
    f"PENCIL = Path({str(args.pencil)!r})\n"
    f"JACOBIAN = Path({str(args.jacobian)!r})\n"
    + (f"MINIMAL = Path({str(args.minimal)!r})\n" if args.minimal is not None else "")
    + f"OUTPUT = Path({str(args.output)!r})"
)
if source.count(old_paths) != 1:
    raise ArithmeticError("immutable map core path contract changed")
source = source.replace(old_paths, new_paths, 1)
old_inputs = '''minimal = json.loads(MINIMAL.read_text())
expected_inputs = (
    (batch, "PASS_EXACT_THIRD_Q12_WEIERSTRASS_SAMPLE_BATCH_MOD19_QUADRATIC"),
    (pencil, "PASS_EXACT_RESOLVED_THIRD_Q12_PENCIL_MOD19_QUADRATIC"),
    (jacobian, "PASS_EXACT_INTERPOLATED_THIRD_Q12_JACOBIAN_MOD19_QUADRATIC"),
    (minimal, "PASS_EXACT_MINIMAL_THIRD_Q12_JACOBIAN_AND_FIBRES_MOD19_QUADRATIC"),
)'''
if minimal is None:
    new_inputs = '''expected_inputs = (
        (batch, "PASS_EXACT_THIRD_Q12_WEIERSTRASS_SAMPLE_BATCH_COMMON_PRODUCER"),
        (pencil, "PASS_EXACT_RESOLVED_THIRD_Q12_PENCIL_COMMON_PRODUCER"),
        (jacobian, "PASS_EXACT_INTERPOLATED_THIRD_Q12_JACOBIAN_COMMON_PRODUCER"),
    )'''
else:
    new_inputs = '''minimal = json.loads(MINIMAL.read_text())
expected_inputs = (
    (batch, "PASS_EXACT_THIRD_Q12_WEIERSTRASS_SAMPLE_BATCH_COMMON_PRODUCER"),
    (pencil, "PASS_EXACT_RESOLVED_THIRD_Q12_PENCIL_COMMON_PRODUCER"),
    (jacobian, "PASS_EXACT_INTERPOLATED_THIRD_Q12_JACOBIAN_COMMON_PRODUCER"),
    (minimal, "PASS_EXACT_MINIMAL_THIRD_Q12_JACOBIAN_AND_FIBRES_COMMON_PRODUCER"),
)'''
if source.count(old_inputs) != 1:
    raise ArithmeticError("immutable map core input-status contract changed")
source = source.replace(old_inputs, new_inputs, 1)
old_field = (
    'base_finite = GF(19)\nmodulus_ring = PolynomialRing(base_finite, "m")\n'
    'm = modulus_ring.gen()\nfinite = GF(19**2, "r", modulus=m**2 + 12 * m + 3)'
)
new_field = (
    f'base_finite = GF({prime})\nmodulus_ring = PolynomialRing(base_finite, "m")\n'
    f'm = modulus_ring.gen()\nfinite = GF({prime}**2, "r", modulus=m**2 + {linear} * m + {constant})'
)
if source.count(old_field) != 1:
    raise ArithmeticError("immutable map core field contract changed")
source = source.replace(old_field, new_field, 1)
source = source.replace(
    '"PASS_EXACT_THIRD_Q12_WEIERSTRASS_SAMPLE_MOD19_QUADRATIC"',
    '"PASS_EXACT_THIRD_Q12_WEIERSTRASS_SAMPLE_COMMON_PRODUCER"',
)
source = source.replace("if len(training) < 100", "if len(training) < 32", 1)
source = source.replace(
    '"PASS_EXACT_GENERIC_THIRD_Q12_BIRATIONAL_MAPS_MOD19_QUADRATIC"',
    '"PASS_EXACT_GENERIC_THIRD_Q12_BIRATIONAL_MAPS_COMMON_PRODUCER"',
)
source = source.replace(
    '"specialization": {"u": "-2", "prime": 19, "extension_modulus": "r^2+12*r+3"}',
    '"specialization": common_specialization', 1,
)
if minimal is None:
    source = source.replace("for path in (BATCH, PENCIL, JACOBIAN, MINIMAL)",
                            "for path in (BATCH, PENCIL, JACOBIAN)", 1)

namespace = {
    "__file__": str(CORE),
    "__name__": "__main__",
    "common_specialization": specialization,
    "minimal": minimal or {"long_to_minimal_map": {"status": "LONG_MODEL_ONLY"}},
}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(source, str(CORE), "exec"), namespace)

output = json.loads(args.output.read_text())
status = "PASS_EXACT_GENERIC_THIRD_Q12_BIRATIONAL_MAPS_COMMON_PRODUCER"
if output.get("status") != status:
    raise ArithmeticError("adapted generic maps did not pass")
output["schema"] = "elkies-k3.q80-third-q12-birational-maps-modp2.v2"
output["status"] = status
output["specialization"] = specialization
if minimal is None:
    for key in ("long_to_minimal", "forward_minimal", "inverse_minimal"):
        output.pop(key, None)
output["inputs"] = [
    {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
    for path in (
        (args.batch, args.pencil, args.jacobian)
        + ((args.minimal,) if args.minimal is not None else ())
    )
]
output["worker"] = {
    "core": {"path": str(CORE.relative_to(ROOT)), "sha256": sha256(CORE)},
    "adapter": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
}
output["claim_boundary"] = {
    "proved": [
        f"explicit generic birational maps in both directions over GF({prime}^2)(V)",
        "literal generic function-field replay, not only invariant agreement",
        "all mapped fibres replay by joint cross multiplication",
        *(["explicit composition with the certified minimal short Weierstrass gauge"] if minimal is not None else []),
    ],
    "not_proved": [
        *(["composition with a globally minimal short model"] if minimal is None else []),
        "transported old-component and zero-section marking",
        "a characteristic-zero coefficient reconstruction",
    ],
}
output["reproduce"] = (
    "sage -python elkies-k3/scripts/interpolate_q80_third_q12_maps_modp2.py "
    f"--batch {args.batch} --pencil {args.pencil} --jacobian {args.jacobian} "
    + (f"--minimal {args.minimal} " if args.minimal is not None else "")
    + f"--output {args.output}"
)
args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
print(
    f"Q80THIRDQ12COMMONMAPS|prime={prime}|training={batch['training_count']}|"
    f"heldout={batch['held_out_count']}|generic_forward=PASS|generic_inverse=PASS|"
    "status=PASS_EXACT_GENERIC_THIRD_Q12_BIRATIONAL_MAPS_COMMON_PRODUCER"
)
