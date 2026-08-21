#!/usr/bin/env sage
"""CRT-reconstruct the generic q=80 CM24 slope-1/12 sextic space.

Inputs are canonical RREF relation spaces produced by
``extend_q80_cm24_branch_modp.sage``.  RREF makes each matrix entry a rational
function of the characteristic-zero jet coefficients, so coefficientwise CRT
and rational reconstruction are meaningful whenever the pivot pattern is
stable.  One declared split prime is excluded from reconstruction and used as
a blind exact validation modulus.
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, QQ, ZZ, crt


parser = argparse.ArgumentParser()
parser.add_argument("inputs", nargs="+")
parser.add_argument("--validation", required=True)
parser.add_argument("--output")
arguments = parser.parse_args()


def load_artifact(path):
    payload = json.loads(Path(path).read_text())
    if (
        payload.get("schema") != "q80-cm24-split-prime-formal-branch-v1"
        or payload.get("kind") != "canonical_centered_relation_space"
        or payload.get("characteristic_zero_slope") != "1/12"
        or payload.get("relation_degree") != 6
    ):
        raise ValueError(f"unexpected relation-space artifact: {path}")
    return payload


training = [load_artifact(path) for path in arguments.inputs]
validation = load_artifact(arguments.validation)
all_payloads = training+[validation]
reference = training[0]
for payload in all_payloads[1:]:
    for key in ("centered_variables", "monomial_exponents", "pivot_columns"):
        if payload[key] != reference[key]:
            raise ValueError(
                f"prime {payload['prime']} does not have the generic {key}"
            )
    if len(payload["rref_basis"]) != len(reference["rref_basis"]):
        raise ValueError(f"prime {payload['prime']} has exceptional kernel dimension")

primes = [ZZ(payload["prime"]) for payload in training]
if len(set(primes+[ZZ(validation["prime"])])) != len(primes)+1:
    raise ValueError("training and validation primes must be distinct")
modulus = ZZ.prod(primes)
rows = len(reference["rref_basis"])
columns = len(reference["rref_basis"][0])
reconstructed = []
failures = []
for row in range(rows):
    reconstructed_row = []
    for column in range(columns):
        residues = [
            ZZ(payload["rref_basis"][row][column])
            for payload in training
        ]
        combined = ZZ(crt(residues, primes))
        try:
            value = QQ(combined.rational_reconstruction(modulus))
        except ArithmeticError:
            failures.append((row, column, "no_rational_reconstruction"))
            value = None
        if value is not None:
            expected = GF(validation["prime"])(
                validation["rref_basis"][row][column]
            )
            try:
                actual = GF(validation["prime"])(value)
            except ZeroDivisionError:
                failures.append((row, column, "validation_denominator_zero"))
            else:
                if actual != expected:
                    failures.append((row, column, "blind_validation_mismatch"))
        reconstructed_row.append(value)
    reconstructed.append(reconstructed_row)

nonzero = [value for row in reconstructed for value in row if value]
maximum_numerator = max((abs(value.numerator()) for value in nonzero), default=0)
maximum_denominator = max((value.denominator() for value in nonzero), default=1)
status = "PASS" if not failures else "NEEDS_MORE_PRIMES"
print(
    f"Q80CM24CRT|training_primes={','.join(map(str, primes))}|"
    f"modulus={modulus}|validation_prime={validation['prime']}|"
    f"rows={rows}|columns={columns}|failures={len(failures)}|"
    f"max_numerator={maximum_numerator}|max_denominator={maximum_denominator}|"
    f"status={status}",
    flush=True,
)
if failures:
    print(f"Q80CM24CRT|first_failures={tuple(failures[:20])}", flush=True)

if arguments.output and not failures:
    payload = {
        "schema": "q80-cm24-slope-1-12-sextic-space-qq-v1",
        "status": "exact_rational_reconstruction_with_blind_modular_validation",
        "characteristic_zero_slope": "1/12",
        "normalization": "P=P_CM24+h",
        "centered_variables": reference["centered_variables"],
        "relation_degree": 6,
        "monomial_exponents": reference["monomial_exponents"],
        "pivot_columns": reference["pivot_columns"],
        "rref_basis": [
            [str(value) for value in row] for row in reconstructed
        ],
        "training_primes": list(map(int, primes)),
        "crt_modulus": str(modulus),
        "blind_validation_prime": int(validation["prime"]),
        "formal_order_at_each_prime": reference["order"],
        "maximum_absolute_numerator": str(maximum_numerator),
        "maximum_denominator": str(maximum_denominator),
        "caveat": (
            "The rational coefficient matrix is reconstructed exactly from "
            "finite fields and passes one excluded prime. Direct substitution "
            "into a characteristic-zero algebraic family is still pending."
        ),
    }
    Path(arguments.output).write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
