#!/usr/bin/env sage-python
"""Export a ToricControlledReduction input for one R17 singleton twist.

For a squarefree quadratic ``d`` the regular model is

    d(t)*y^2 = x^3 + A(t)*x + B(t).

Its Newton polytope has vertices ``(0,0,0)``, ``(12,0,0)``, ``(0,3,0)``,
``(0,0,2)``, and ``(2,0,2)``.  The independent certifier re-derives the
input from the rational source artifacts before accepting the raw output.
"""

import argparse
from hashlib import sha256
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BISECTIONS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisections-cheapest-1024-v1.json"
)
DEFAULT_MODEL = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
)


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def ntl_matrix(rows) -> str:
    return "[" + "".join("[" + " ".join(map(str, row)) + "]" for row in rows) + "]"


def ntl_vector(values) -> str:
    return "[" + " ".join(map(str, values)) + "]"


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--label", required=True)
parser.add_argument("--prime", type=int, required=True)
parser.add_argument("--bisections", type=Path, default=DEFAULT_BISECTIONS)
parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
parser.add_argument(
    "--output-dir",
    type=Path,
    default=ROOT / "artifacts/local/elkies-k3/r17-singleton-toric-frobenius",
)
args = parser.parse_args()

prime = int(args.prime)
field = GF(prime)
if prime < 5 or not field.is_prime_field():
    raise ValueError("--prime must be an odd prime at least five")

bisections = json.loads(args.bisections.read_text())
if bisections.get("schema") != "elkies-k3.bisection-extension-input.v1":
    raise ValueError("unexpected bisection schema")
by_label = {record["label"]: record for record in bisections["bisections"]}
if args.label not in by_label:
    raise ValueError("label absent from bisection input")

Rt = PolynomialRing(field, "t")
d = Rt(
    [
        field(QQ(value))
        for value in by_label[args.label]["branch"]["numerator_coefficients"]
    ]
)
model = json.loads(args.model.read_text())
if model.get("status") != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
    raise ValueError("unexpected direct-model status")
weierstrass = model["weierstrass_model"]
A = Rt([field(QQ(value)) for value in weierstrass["A_coefficients_low_to_high"]])
B = Rt([field(QQ(value)) for value in weierstrass["B_coefficients_low_to_high"]])
base_discriminant = -field(16) * (field(4) * A**3 + field(27) * B**2)
if (
    d.degree() != 2
    or not d.is_squarefree()
    or A.degree() != 8
    or B.degree() != 12
    or base_discriminant.degree() != 24
    or not base_discriminant.is_squarefree()
    or d.gcd(base_discriminant).degree()
):
    raise ArithmeticError("declared prime is not good for the singleton twist")

terms = {(0, 3, 0): -field.one()}
for index, coefficient in enumerate(d):
    if coefficient:
        terms[(index, 0, 2)] = coefficient
for index, coefficient in enumerate(-A):
    if coefficient:
        terms[(index, 1, 0)] = coefficient
for index, coefficient in enumerate(-B):
    if coefficient:
        terms[(index, 0, 0)] = coefficient
monomials = sorted(terms)
coefficients = [int(terms[monomial]) for monomial in monomials]

required_vertices = {
    (0, 0, 0),
    (12, 0, 0),
    (0, 3, 0),
    (0, 0, 2),
    (2, 0, 2),
}
if not required_vertices.issubset(terms):
    raise ArithmeticError("support does not span the expected Newton polytope")
halfspace_A = (
    (0, 0, 1),
    (0, 1, 0),
    (-1, -4, -5),
    (0, -2, -3),
    (1, 0, 0),
)
halfspace_b = (0, 0, 12, 6, 0)

output_dir = args.output_dir.resolve() / args.label / f"p{prime}"
output_dir.mkdir(parents=True, exist_ok=True)
input_path = output_dir / "toric-controlled-reduction.input"
label = f"r17-singleton-{args.label}-p{prime}"
line = ":".join(
    (
        label,
        ntl_matrix(monomials),
        ntl_vector(coefficients),
        ntl_matrix(halfspace_A),
        ntl_vector(halfspace_b),
        str(prime),
    )
)
input_path.write_text(line + "\n")

record = {
    "schema": "elkies-k3.r17-singleton-toric-frobenius-input.v1",
    "status": "PASS_EXACT_TORIC_FROBENIUS_INPUT_EXPORT",
    "label": args.label,
    "prime": prime,
    "surface_equation": "d(t)*y^2-x^3-A(t)*x-B(t)",
    "quadratic_coefficients_low_to_high_mod_p": [int(value) for value in d],
    "support_monomials_t_x_y": [list(value) for value in monomials],
    "support_coefficients_mod_p": coefficients,
    "newton_halfspace_A": [list(row) for row in halfspace_A],
    "newton_halfspace_b": list(halfspace_b),
    "newton_vertices": [list(value) for value in sorted(required_vertices)],
    "toric_input": {
        "path": str(input_path.relative_to(ROOT)),
        "sha256": digest(input_path),
    },
    "inputs": {
        str(path.resolve().relative_to(ROOT)): digest(path)
        for path in (args.bisections, args.model)
    },
}
record_path = output_dir / "input-certificate.json"
record_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
print(
    f"R17SINGLETONTORICINPUT|label={args.label}|p={prime}|terms={len(monomials)}"
    f"|input={input_path}|certificate={record_path}|status=PASS_EXPORTED",
    flush=True,
)

