#!/usr/bin/env sage-python
"""Specialize coefficients in an NS0007 reduced modular system.

The four ``a2`` entries remain the primary slice.  Repeated ``--fix`` options
can additionally pin node-jet variables, for example
``--fix si_0=1 --fix sl_0=3``.  This supports exact complexity probes before
committing to a large finite-field slice census.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "artifacts/local/elkies-k3/ns0007-pole0-reduced-modp/p7-lambda2.json"
DEFAULT_OUTPUT_DIR = ROOT / "artifacts/local/elkies-k3/ns0007-pole0-reduced-modp/a2-slices"


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
parser.add_argument(
    "--a2",
    required=True,
    help="four comma-separated values or * for a2_4,a2_3,a2_2,a2_1",
)
parser.add_argument(
    "--fix",
    action="append",
    default=[],
    help="additional variable specialization NAME=VALUE (repeatable)",
)
parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()

metadata_path = args.input.resolve()
metadata = json.loads(metadata_path.read_text())
if metadata.get("schema") != "elkies-k3.lattice-foundry-ns0007-pole0-reduced-modp-system.v1":
    raise ValueError("unexpected reduced-system schema")
prime = int(metadata["prime"])
field = GF(prime)
tokens = args.a2.split(",")
if len(tokens) != 4:
    raise SystemExit("--a2 requires four values")
values = [None if token.strip() == "*" else field(int(token)) for token in tokens]

system_path = ROOT / metadata["system"]["msolve_input"]
if hashlib.sha256(system_path.read_bytes()).hexdigest() != metadata["system"]["msolve_sha256"]:
    raise ArithmeticError("reduced-system digest mismatch")
lines = system_path.read_text().splitlines()
old_names = lines[0].split(",")
old_ring = PolynomialRing(field, names=old_names, order="degrevlex")
equation_text = "\n".join(lines[2:]).replace("^", "**")
old_equations = [old_ring(piece.strip()) for piece in equation_text.split(",") if piece.strip()]

fixed_names = ["a2_4", "a2_3", "a2_2", "a2_1"]
fixed = {
    name: value
    for name, value in zip(fixed_names, values)
    if value is not None
}
extra_fixed = {}
for assignment in args.fix:
    if assignment.count("=") != 1:
        raise SystemExit("--fix requires NAME=VALUE")
    name, value_text = assignment.split("=", 1)
    name = name.strip()
    if name not in old_names:
        raise SystemExit(f"unknown --fix variable {name!r}")
    if name in fixed_names:
        raise SystemExit(f"specialize {name} through --a2, not --fix")
    if name in extra_fixed:
        raise SystemExit(f"duplicate --fix variable {name!r}")
    extra_fixed[name] = field(int(value_text))
fixed.update(extra_fixed)
new_names = [name for name in old_names if name not in fixed]
new_ring = PolynomialRing(field, names=new_names, order="degrevlex")
new_generators = new_ring.gens_dict()
images = [fixed[name] if name in fixed else new_generators[name] for name in old_names]
specialize = old_ring.hom(images, new_ring)
equations = []
inconsistent_constant = False
for equation in old_equations:
    specialized = new_ring(specialize(equation))
    if specialized.is_constant():
        if specialized:
            inconsistent_constant = True
            break
        continue
    equations.append(specialized)
if inconsistent_constant:
    equations = [new_ring.one()]

tag = "-".join("x" if value is None else str(int(value)) for value in values)
extra_tag = "".join(
    f"-{name}-{int(value)}" for name, value in sorted(extra_fixed.items())
)
output_dir = args.output_dir.resolve()
output_dir.mkdir(parents=True, exist_ok=True)
stem = f"p{prime}-lambda{metadata['lambda']}-a2-{tag}{extra_tag}"
msolve_path = output_dir / f"{stem}.ms"
output_metadata_path = output_dir / f"{stem}.json"
msolve_text = ",".join(new_names) + "\n" + str(prime) + "\n"
msolve_text += ",\n".join(str(equation).replace("**", "^") for equation in equations) + "\n"
output = {
    "schema": "elkies-k3.lattice-foundry-ns0007-pole0-reduced-a2-slice-modp.v1",
    "status": "PASS_EXACT_A2_SPECIALIZATION",
    "input": {
        "metadata": relative(metadata_path),
        "metadata_sha256": hashlib.sha256(metadata_path.read_bytes()).hexdigest(),
        "system": relative(system_path),
        "system_sha256": hashlib.sha256(system_path.read_bytes()).hexdigest(),
    },
    "prime": prime,
    "lambda": metadata["lambda"],
    "a2_coefficients_high_to_low": [
        None if value is None else int(value) for value in values
    ],
    "additional_fixed_variables": {
        name: int(value) for name, value in sorted(extra_fixed.items())
    },
    "system": {
        "variables": new_names,
        "variable_count": len(new_names),
        "equation_count": len(equations),
        "inconsistent_constant": inconsistent_constant,
        "equation_total_degrees": [int(equation.degree()) for equation in equations],
        "equation_term_counts": [len(equation.monomials()) for equation in equations],
        "msolve_input": relative(msolve_path),
        "msolve_sha256": hashlib.sha256(msolve_text.encode()).hexdigest(),
    },
    "proof_boundary": (
        "Exact substitution fixes the displayed a2 polynomial in the reduced "
        "finite-field NS0007 chart and any explicitly named node jets. Solver "
        "output still requires independent "
        "decoding and geometric verification."
    ),
}
metadata_text = json.dumps(output, indent=2, sort_keys=True) + "\n"
if args.check:
    if msolve_path.read_text() != msolve_text or output_metadata_path.read_text() != metadata_text:
        raise SystemExit("NS0007 a2 slice is stale")
else:
    msolve_path.write_text(msolve_text)
    output_metadata_path.write_text(metadata_text)

print(
        "FOUNDRYNS0007A2SLICE|"
        f"p={prime}|lambda={metadata['lambda']}|a2={tag}|"
        f"extra_fixed={len(extra_fixed)}|"
    f"vars={len(new_names)}|eqs={len(equations)}|"
    f"max_degree={max(equation.degree() for equation in equations)}|status=PASS",
    flush=True,
)
print(f"MSOLVE_INPUT|{msolve_path}", flush=True)
print(f"OUTPUT|{output_metadata_path}", flush=True)
