#!/usr/bin/env sage
"""Recover short integral sextics from modular q=80 relation spaces.

The canonical RREF of the 15-dimensional relation space has very large
rational entries.  Instead of reconstructing that ill-conditioned basis,
this script forms the full integer congruence lattice

    {v in Z^210 : v mod M lies in the modular relation space}

for a CRT modulus M, LLL-reduces it, and retains the 15 shortest independent
vectors which also pass a completely excluded validation prime.
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, Matrix, QQ, ZZ, crt, vector


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
        or payload.get("relation_degree") != 6
    ):
        raise ValueError(f"unexpected relation artifact: {path}")
    return payload


training = [load_artifact(path) for path in arguments.inputs]
validation = load_artifact(arguments.validation)
reference = training[0]
for payload in training[1:]+[validation]:
    for key in ("centered_variables", "monomial_exponents", "pivot_columns"):
        if payload[key] != reference[key]:
            raise ValueError(f"prime {payload['prime']} changed {key}")
    if len(payload["rref_basis"]) != 15:
        raise ValueError(f"prime {payload['prime']} is an exceptional reduction")

primes = [ZZ(payload["prime"]) for payload in training]
modulus = ZZ.prod(primes)
pivots = tuple(reference["pivot_columns"])
column_count = len(reference["monomial_exponents"])
free_columns = tuple(column for column in range(column_count) if column not in pivots)

basis_rows = []
for row in range(15):
    values = [ZZ.zero() for _ in range(column_count)]
    for column in range(column_count):
        residues = [
            ZZ(payload["rref_basis"][row][column]) for payload in training
        ]
        values[column] = ZZ(crt(residues, primes))
        if values[column] > modulus//2:
            values[column] -= modulus
    basis_rows.append(values)
for column in free_columns:
    values = [ZZ.zero() for _ in range(column_count)]
    values[column] = modulus
    basis_rows.append(values)
lattice_basis = Matrix(ZZ, basis_rows)
assert lattice_basis.nrows() == lattice_basis.ncols() == column_count

print(
    f"Q80CM24SHORT|stage=lll_start|training_primes={','.join(map(str, primes))}|"
    f"modulus={modulus}|dimension={column_count}",
    flush=True,
)
reduced = lattice_basis.LLL(delta=0.75)
print("Q80CM24SHORT|stage=lll_complete", flush=True)


def lies_in_space(row, payload):
    field = GF(payload["prime"])
    basis = Matrix(field, payload["rref_basis"])
    candidate = vector(field, list(row))
    return basis.stack(Matrix(field, [candidate])).rank() == 15


selected = []
selected_matrix = Matrix(QQ, 0, column_count)
diagnostics = []
for row in sorted(reduced.rows(), key=lambda value: value.norm()):
    if not all(lies_in_space(row, payload) for payload in training):
        raise AssertionError("LLL returned a vector outside the CRT lattice")
    validation_ok = lies_in_space(row, validation)
    trial = selected_matrix.stack(Matrix(QQ, [row]))
    independent = trial.rank() > selected_matrix.rank()
    diagnostics.append((row.norm()**2, max(map(abs, row)), validation_ok, independent))
    if validation_ok and independent:
        selected.append(vector(ZZ, row))
        selected_matrix = trial
        if len(selected) == 15:
            break

status = "PASS" if len(selected) == 15 else "NEEDS_LARGER_MODULUS"
maximum_coefficient = max(
    (max(map(abs, row)) for row in selected), default=ZZ.zero()
)
print(
    f"Q80CM24SHORT|training_primes={','.join(map(str, primes))}|"
    f"modulus={modulus}|validation_prime={validation['prime']}|"
    f"selected={len(selected)}|max_coefficient={maximum_coefficient}|"
    f"first_diagnostics={tuple(diagnostics[:20])}|status={status}",
    flush=True,
)

if arguments.output and len(selected) == 15:
    payload = {
        "schema": "q80-cm24-slope-1-12-short-sextic-space-qq-v1",
        "status": "integral_CRT_lattice_recovery_with_blind_modular_validation",
        "characteristic_zero_slope": "1/12",
        "normalization": "P=P_CM24+h",
        "centered_variables": reference["centered_variables"],
        "relation_degree": 6,
        "monomial_exponents": reference["monomial_exponents"],
        "integral_basis": [list(map(int, row)) for row in selected],
        "training_primes": list(map(int, primes)),
        "crt_modulus": str(modulus),
        "blind_validation_prime": int(validation["prime"]),
        "maximum_absolute_coefficient": str(maximum_coefficient),
        "formal_order_at_each_prime": reference["order"],
        "caveat": (
            "The short integral space is recovered from congruence lattices "
            "and passes an excluded split prime. A direct characteristic-zero "
            "formal or global substitution certificate is still required."
        ),
    }
    Path(arguments.output).write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
