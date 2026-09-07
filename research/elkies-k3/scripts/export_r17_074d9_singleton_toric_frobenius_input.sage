#!/usr/bin/env sage-python
"""Export toric controlled-reduction inputs for the four 074d9 record twists.

status: ACTIVE_COMPILER
claim: exact sparse toric-hypersurface input at one declared good prime

For a quadratic character ``q(u)``, use the regular twist presentation

    q(u)*y^2 = x^3 + A(u)*x + B(u).

It is birational over the generic fibre to

    Y^2 = X^3 + q(u)^2*A(u)*X + q(u)^3*B(u)

through ``X=q*x, Y=q^2*y``.  The Newton polytope has vertices

    (0,0,0), (12,0,0), (0,3,0), (0,0,2), (2,0,2)

and half-space description

    u,x,y >= 0,  u+4*x+5*y <= 12,  2*x+3*y <= 6.

The output line is accepted by the open-source ToricControlledReduction
``examples/readfile.cpp`` driver.  Nondegeneracy, Frobenius computation, and
the conversion of its output into a Picard/Mordell--Weil bound are separate
certificates.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_COVERS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-074d9-cross-fibre-bisection-transfer-v1.json"
)
DEFAULT_MODEL = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json"
)
RECORD_LABELS = (
    "074d9-orbit-04b07",
    "074d9-orbit-11a44",
    "074d9-orbit-11279",
    "074d9-orbit-080fa",
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
parser.add_argument("--label", choices=RECORD_LABELS, required=True)
parser.add_argument("--prime", type=int, required=True)
parser.add_argument("--covers", type=Path, default=DEFAULT_COVERS)
parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
parser.add_argument(
    "--output-dir",
    type=Path,
    default=ROOT / "artifacts/local/elkies-k3/r17-074d9-singleton-toric-frobenius",
)
args = parser.parse_args()

prime = int(args.prime)
field = GF(prime)
if prime < 5 or not field.is_prime_field():
    raise ValueError("--prime must be an odd prime at least five")

covers = json.loads(args.covers.read_text())
if covers.get("status") != "PASS_EXACT_COMPLETE_074D9_CROSS_FIBRE_BISECTION_TRANSFER":
    raise ValueError("unexpected 074d9 cover-certificate status")
by_label = {
    row["label"]: row
    for fibre in covers["fibres"]
    for row in fibre["records"]
    if row["label"] in RECORD_LABELS
}
if tuple(label for label in RECORD_LABELS if label not in by_label):
    raise ArithmeticError("the cover certificate lost a record-specific twist")

model = json.loads(args.model.read_text())
if model.get("status") != "PROVED_EXACT_LINEAGE_REALIZATION_AND_DISPLAYED_QUOTIENTS":
    raise ValueError("unexpected 074d9 lineage-model status")
representative = model["representative"]

base_ring = PolynomialRing(field, "u")
q = base_ring(
    [field(QQ(value)) for value in by_label[args.label]["branch_quadratic_coefficients_low_to_high"]]
)
A = base_ring(
    [field(QQ(value)) for value in representative["A_coefficients_low_to_high"]]
)
B = base_ring(
    [field(QQ(value)) for value in representative["B_coefficients_low_to_high"]]
)
base_discriminant = -field(16) * (field(4) * A**3 + field(27) * B**2)
if (
    q.degree() != 2
    or not q.is_squarefree()
    or A.degree() != 8
    or B.degree() != 12
    or base_discriminant.degree() != 24
    or not base_discriminant.is_squarefree()
    or q.gcd(base_discriminant).degree()
):
    raise ArithmeticError("declared prime is not good for the singleton-twist model")

# Variables are ordered (u,x,y).  Coefficients are represented in 0,...,p-1.
terms = {(0, 3, 0): -field.one()}
for index, coefficient in enumerate(q):
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

tag = args.label.removeprefix("074d9-orbit-")
output_dir = args.output_dir.resolve() / tag / f"p{prime}"
output_dir.mkdir(parents=True, exist_ok=True)
input_path = output_dir / "toric-controlled-reduction.input"
driver_label = f"r17-074d9-singleton-{tag}-p{prime}"
line = ":".join(
    (
        driver_label,
        ntl_matrix(monomials),
        ntl_vector(coefficients),
        ntl_matrix(halfspace_A),
        ntl_vector(halfspace_b),
        str(prime),
    )
)
input_path.write_text(line + "\n")

record = {
    "schema": "elkies-k3.r17-074d9-singleton-toric-frobenius-input.v1",
    "status": "PASS_EXACT_TORIC_FROBENIUS_INPUT_EXPORT",
    "proof_boundary": (
        "This certifies the finite-field coefficients, good-reduction gates, full "
        "Newton polytope, and ToricControlledReduction input syntax. Toric "
        "nondegeneracy and Frobenius output are separate certificates."
    ),
    "label": args.label,
    "prime": prime,
    "surface_equation": "q(u)*y^2-x^3-A(u)*x-B(u)",
    "generic_fibre_isomorphism_to_short_twist": "X=q*x, Y=q^2*y",
    "quadratic_coefficients_low_to_high_mod_p": [int(value) for value in q],
    "support_monomials_u_x_y": [list(value) for value in monomials],
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
        for path in (args.covers, args.model)
    },
}
record_path = output_dir / "input-certificate.json"
record_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
print(
    f"R17074D9TORICINPUT|label={args.label}|p={prime}|terms={len(monomials)}"
    f"|input={input_path}|certificate={record_path}|status=PASS_EXPORTED",
    flush=True,
)
