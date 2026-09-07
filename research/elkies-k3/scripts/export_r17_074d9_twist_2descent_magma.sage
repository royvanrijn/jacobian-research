#!/usr/bin/env sage-python
"""Export a fail-closed Magma 2-descent job for one 074d9 singleton twist.

Magma's function-field ``TwoSelmerGroup`` is currently implemented for
``F_p(u)`` in odd characteristic, not for ``QQ(u)``.  At a good prime the
result supplies an unconditional arithmetic rank upper bound for the
characteristic-zero Mordell--Weil group.  This exporter fixes every
coefficient and emits the exact job; execution and interpretation of the raw
Magma output are separate certificate steps.
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


def magma_polynomial(polynomial, variable: str = "u") -> str:
    terms = []
    for degree, coefficient in enumerate(polynomial):
        if not coefficient:
            continue
        if degree == 0:
            monomial = str(int(coefficient))
        elif degree == 1:
            monomial = f"{int(coefficient)}*{variable}"
        else:
            monomial = f"{int(coefficient)}*{variable}^{degree}"
        terms.append(monomial)
    return " + ".join(terms) if terms else "0"


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--label", choices=RECORD_LABELS, required=True)
parser.add_argument("--prime", type=int, required=True)
parser.add_argument(
    "--allow-nongood-reduction",
    action="store_true",
    help="export a discovery-only job when singular surface fibres collide mod p",
)
parser.add_argument("--covers", type=Path, default=DEFAULT_COVERS)
parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
parser.add_argument(
    "--output-dir",
    type=Path,
    default=ROOT / "artifacts/local/elkies-k3/r17-074d9-twist-2descent",
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

model = json.loads(args.model.read_text())
if model.get("status") != "PROVED_EXACT_LINEAGE_REALIZATION_AND_DISPLAYED_QUOTIENTS":
    raise ValueError("unexpected 074d9 lineage-model status")
representative = model["representative"]

ring = PolynomialRing(field, "u")
q = ring(
    [
        field(QQ(value))
        for value in by_label[args.label]["branch_quadratic_coefficients_low_to_high"]
    ]
)
A = ring(
    [field(QQ(value)) for value in representative["A_coefficients_low_to_high"]]
)
B = ring(
    [field(QQ(value)) for value in representative["B_coefficients_low_to_high"]]
)
base_discriminant = -field(16) * (field(4) * A**3 + field(27) * B**2)
coefficient_model_valid = (
    q.degree() == 2
    and q.is_squarefree()
    and A.degree() == 8
    and B.degree() == 12
    and base_discriminant != 0
)
good_reduction = (
    coefficient_model_valid
    and base_discriminant.degree() == 24
    and base_discriminant.is_squarefree()
    and q.gcd(base_discriminant).degree() == 0
)
if not coefficient_model_valid or (
    not good_reduction and not args.allow_nongood_reduction
):
    raise ArithmeticError("declared prime is not permitted for the singleton twist")

tag = args.label.removeprefix("074d9-orbit-")
output_dir = args.output_dir.resolve() / tag / f"p{prime}"
output_dir.mkdir(parents=True, exist_ok=True)
job_path = output_dir / "two-descent.m"
lines = [
    f'print "R17074D9_2DESCENT|label={args.label}|prime={prime}";',
    f"F := GF({prime});",
    "K<u> := FunctionField(F);",
    f"q := {magma_polynomial(q)};",
    f"A := {magma_polynomial(A)};",
    f"B := {magma_polynomial(B)};",
    "E := EllipticCurve([K | 0, 0, 0, q^2*A, q^3*B]);",
    "R<x> := PolynomialRing(K);",
    "two_factors := Factorization(x^3 + q^2*A*x + q^3*B);",
    'print "TWO_DIVISION_FACTOR_DEGREES", [ Degree(row[1]) : row in two_factors ];',
    "S2, selmer_map := TwoSelmerGroup(E);",
    'print "TWO_SELMER_INVARIANTS", Invariants(S2);',
    'print "TWO_SELMER_ORDER", #S2;',
    'print "TWO_SELMER_NGENS", Ngens(S2);',
    "two_torsion_dimension := #two_factors - 1;",
    'print "TWO_TORSION_DIMENSION", two_torsion_dimension;',
    'print "MW_RANK_UPPER_FROM_TWO_SELMER", Ngens(S2)-two_torsion_dimension;',
]
job_path.write_text("\n".join(lines) + "\n")
record = {
    "schema": "elkies-k3.r17-074d9-twist-2descent-input.v1",
    "status": (
        "PASS_EXACT_MAGMA_2DESCENT_INPUT_EXPORT"
        if good_reduction
        else "PASS_EXACT_MAGMA_2DESCENT_DISCOVERY_INPUT_EXPORT"
    ),
    "proof_boundary": (
        "This certifies the good reduction and exact Magma input only. The raw "
        "Magma output and any rank conclusion require separate verification."
    ),
    "label": args.label,
    "prime": prime,
    "good_reduction": bool(good_reduction),
    "descent_role": "rank_upper_bound" if good_reduction else "discovery_only",
    "short_twist_equation": "y^2=x^3+q(u)^2*A(u)*x+q(u)^3*B(u)",
    "quadratic_coefficients_low_to_high_mod_p": [int(value) for value in q],
    "A_coefficients_low_to_high_mod_p": [int(value) for value in A],
    "B_coefficients_low_to_high_mod_p": [int(value) for value in B],
    "magma_job": {
        "path": str(job_path.relative_to(ROOT)),
        "sha256": digest(job_path),
    },
    "inputs": {
        str(path.resolve().relative_to(ROOT)): digest(path)
        for path in (args.covers, args.model)
    },
}
record_path = output_dir / "input-certificate.json"
record_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
print(
    f"R17074D92DESCENTINPUT|label={args.label}|p={prime}|job={job_path}"
    f"|certificate={record_path}|status=PASS_EXPORTED",
    flush=True,
)
