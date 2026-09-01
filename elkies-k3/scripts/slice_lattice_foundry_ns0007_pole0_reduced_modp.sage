#!/usr/bin/env sage-python
"""Specialize the four a2 coefficients in an NS0007 reduced modular system."""

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
output_dir = args.output_dir.resolve()
output_dir.mkdir(parents=True, exist_ok=True)
stem = f"p{prime}-lambda{metadata['lambda']}-a2-{tag}"
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
        "finite-field NS0007 chart. Solver output still requires independent "
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
    f"vars={len(new_names)}|eqs={len(equations)}|"
    f"max_degree={max(equation.degree() for equation in equations)}|status=PASS",
    flush=True,
)
print(f"MSOLVE_INPUT|{msolve_path}", flush=True)
print(f"OUTPUT|{output_metadata_path}", flush=True)
