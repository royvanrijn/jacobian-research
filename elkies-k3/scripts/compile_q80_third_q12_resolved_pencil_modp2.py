#!/usr/bin/env python3
"""Adapt the immutable p=19 connected-quotient compiler to any quadratic prime."""

import argparse
import contextlib
import hashlib
import io
import json
import re
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CORE = ROOT / "elkies-k3/scripts/compile_q80_third_q12_um2_p19_resolved_pencil.sage"
LATTICE = ROOT / "artifacts/generated-results/q80-d7d5-mw5-height-lattice.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--surface", type=Path, required=True)
parser.add_argument("--horizontal", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
args.surface = args.surface.resolve()
args.horizontal = args.horizontal.resolve()
args.output = args.output.resolve()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


surface = json.loads(args.surface.read_text())
horizontal = json.loads(args.horizontal.read_text())
if horizontal.get("status") != "PASS_EXACT_THIRD_Q12_HORIZONTAL_FROM_COMMON_CLOSURE_PRODUCER":
    raise ValueError("common-producer horizontal is not certified")
if len(horizontal["pairwise_producer"]["candidates"]) != 1:
    raise ValueError("resolved compiler requires one Frobenius orbit representative")
specialization = horizontal["specialization"]
prime = int(specialization["prime"])
if specialization["common_extension_degree"] != 2:
    raise ValueError("resolved compiler currently requires a quadratic representative")
if surface["parameters"][0]["u"] != specialization["u"]:
    raise ValueError("surface/horizontal u mismatch")

candidate = horizontal["pairwise_producer"]["candidates"][0]


def coefficient_string(coordinates):
    a, b = map(int, coordinates)
    pieces = []
    if b:
        pieces.append("r" if b == 1 else f"{b}*r")
    if a:
        pieces.append(str(a))
    return " + ".join(pieces) if pieces else "0"


normalized_horizontal = {
    "schema": "elkies-k3.q80-po0-rur-third-q12-modp.v1",
    "status": "PASS_EXACT_MODP2_THIRD_Q12_HORIZONTAL_FROBENIUS_ORBIT",
    "third_q12": {
        "candidates_up_to_sign": [
            {
                "x": {
                    "numerator_coefficients_low_to_high": [
                        coefficient_string(value)
                        for value in candidate["x"]["numerator_coefficients_low_to_high"]
                    ],
                    "denominator_coefficients_low_to_high": [
                        coefficient_string(value)
                        for value in candidate["x"]["denominator_coefficients_low_to_high"]
                    ],
                },
                "y": {
                    "numerator_coefficients_low_to_high": [
                        coefficient_string(value)
                        for value in candidate["y"]["numerator_coefficients_low_to_high"]
                    ],
                    "denominator_coefficients_low_to_high": [
                        coefficient_string(value)
                        for value in candidate["y"]["denominator_coefficients_low_to_high"]
                    ],
                },
            }
        ]
    },
}

modulus_text = specialization["common_extension_modulus"]
modulus_match = re.fullmatch(r"x\^2 \+ (?:(\d+)\*)?x \+ (\d+)", modulus_text)
if modulus_match is None:
    raise ValueError(f"cannot parse monic quadratic modulus {modulus_text!r}")
linear = int(modulus_match.group(1) or 1)
constant = int(modulus_match.group(2))

local_root = ROOT / "artifacts/local/elkies-k3"
local_root.mkdir(parents=True, exist_ok=True)
with tempfile.TemporaryDirectory(dir=local_root) as temporary_directory:
    temporary = Path(temporary_directory)
    normalized_path = temporary / "horizontal.json"
    core_output = temporary / "resolved.json"
    normalized_path.write_text(json.dumps(normalized_horizontal, indent=2, sort_keys=True) + "\n")

    source = CORE.read_text()
    path_pattern = re.compile(
        r"SURFACE = ROOT / \(.*?DEFAULT_OUTPUT = ROOT / \(.*?\n\)",
        re.DOTALL,
    )
    replacement = (
        f"SURFACE = Path({str(args.surface)!r})\n"
        f"HORIZONTAL = Path({str(normalized_path)!r})\n"
        f"LATTICE = Path({str(LATTICE)!r})\n"
        f"DEFAULT_OUTPUT = Path({str(core_output)!r})"
    )
    source, path_replacements = path_pattern.subn(replacement, source, count=1)
    source, prime_replacements = source.replace("prime = 19", f"prime = {prime}", 1), source.count("prime = 19")
    modulus_old = "r_modulus**2 + 12 * r_modulus + 3"
    modulus_new = f"r_modulus**2 + {linear} * r_modulus + {constant}"
    if path_replacements != 1 or prime_replacements != 1 or source.count(modulus_old) != 1:
        raise ArithmeticError("immutable core no longer matches the adapter contract")
    source = source.replace(modulus_old, modulus_new, 1)
    pole_old = '''h_factors = tuple(x_denominator.factor())
if len(h_factors) != 1 or int(h_factors[0][1]) != 2:
    raise ArithmeticError("horizontal x denominator is not one square")
h = h_factors[0][0].monic()'''
    pole_new = '''h_factors = tuple(x_denominator.factor())
if any(int(exponent) != 2 for factor, exponent in h_factors):
    raise ArithmeticError("horizontal x denominator is not a square")
h = base_ring.one()
for factor, exponent in h_factors:
    h *= factor.monic()
if h.degree() != 2:
    raise ArithmeticError("horizontal pole divisor does not have degree two")'''
    if source.count(pole_old) != 1:
        raise ArithmeticError("immutable core pole check no longer matches the adapter contract")
    source = source.replace(pole_old, pole_new, 1)

    saved_argv = sys.argv
    namespace = {"__file__": str(CORE), "__name__": "__main__"}
    sys.argv = [str(CORE), "--output", str(core_output)]
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            exec(compile(source, str(CORE), "exec"), namespace)
    finally:
        sys.argv = saved_argv
    compiled = json.loads(core_output.read_text())

if compiled.get("status") != "PASS_EXACT_RESOLVED_THIRD_Q12_PENCIL_MOD19_QUADRATIC":
    raise ArithmeticError("adapted core did not pass its exact connected-quotient gates")
compiled["schema"] = "elkies-k3.q80-third-q12-resolved-pencil-modp2.v2"
compiled["status"] = "PASS_EXACT_RESOLVED_THIRD_Q12_PENCIL_COMMON_PRODUCER"
compiled["specialization"] = {
    "u": specialization["u"],
    "prime": prime,
    "extension_modulus": f"r^2 + {linear}*r + {constant}",
}
compiled["inputs"] = {
    "surface": {"path": str(args.surface.relative_to(ROOT)), "sha256": sha256(args.surface)},
    "horizontal": {
        "path": str(args.horizontal.relative_to(ROOT)),
        "sha256": sha256(args.horizontal),
    },
    "lattice": {"path": str(LATTICE.relative_to(ROOT)), "sha256": sha256(LATTICE)},
    "core": {"path": str(CORE.relative_to(ROOT)), "sha256": sha256(CORE)},
    "adapter": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
}
compiled["claim_boundary"] = {
    "proved": [
        "common-producer horizontal literal replay on the parent",
        "complete finite-prime Smith saturation and shifted-Popov ambient",
        "resolved D7 complete ideal and finite D5 connected quotient",
        "rank-five gate with a two-dimensional pencil",
        "exact moving equation of degrees (2,9,3)",
    ],
    "not_proved": [
        "generic genus one or a minimal child Jacobian at this prime",
        "canonical PGL2/Weierstrass alignment with p=19",
        "a characteristic-zero lift",
    ],
}
compiled["reproduce"] = (
    "sage -python elkies-k3/scripts/compile_q80_third_q12_resolved_pencil_modp2.py "
    f"--surface {args.surface} --horizontal {args.horizontal} --output {args.output}"
)
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(compiled, indent=2, sort_keys=True) + "\n")
print(
    f"Q80THIRDQ12COMMONPENCIL|u={specialization['u']}|prime={prime}|"
    "smith=0,0,6|ambient=7|gate_rank=5|kernel=2|moving_degrees=2,9,3|"
    "status=PASS_EXACT_RESOLVED_THIRD_Q12_PENCIL_COMMON_PRODUCER"
)
