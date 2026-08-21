#!/usr/bin/env sage
"""CRT-reconstruct the slope-8/87 q80 quartic relation ideal over QQ.

The split-prime artifacts contain the canonical RREF of the 22-dimensional
space of centered quartic relations in D,P,Q,E.  Canonical pivots make every
matrix entry a well-defined rational number.  We reconstruct from eight large
primes, verify two withheld primes, and evaluate all 22 relations on the exact
order-28 characteristic-zero jet.
"""

import argparse
import glob
import hashlib
import json
from pathlib import Path

from sage.all import *


parser = argparse.ArgumentParser()
parser.add_argument("--input-glob", action="append", required=True)
parser.add_argument("--validation-glob", action="append", default=[])
parser.add_argument(
    "--qq-series",
    default="artifacts/generated-results/q80-cm24-slope-8-87-qq-order28-series.json",
)
parser.add_argument(
    "--output",
    default="artifacts/generated-results/q80-cm24-slope-8-87-qq-quartic-ideal.json",
)
arguments = parser.parse_args()


def expand(patterns):
    return tuple(
        sorted({Path(name) for pattern in patterns for name in glob.glob(pattern)})
    )


input_paths = expand(arguments.input_glob)
validation_paths = expand(arguments.validation_glob)
if not input_paths:
    raise ValueError("no CRT relation artifacts matched")
if set(input_paths) & set(validation_paths):
    raise ValueError("CRT and withheld validation inputs must be disjoint")


def relation_record(path):
    payload = json.loads(path.read_text())
    if payload.get("schema") != "q80-cm24-split-prime-formal-branch-v1":
        raise ValueError(f"unexpected schema in {path}")
    if payload.get("kind") != "canonical_centered_relation_space":
        raise ValueError(f"expected a relation space in {path}")
    if payload.get("characteristic_zero_slope") != "8/87":
        raise ValueError(f"wrong branch in {path}")
    if tuple(payload.get("centered_variables", ())) != ("D", "P", "Q", "E"):
        raise ValueError(f"wrong centered coordinates in {path}")
    if payload.get("relation_degree") != 4:
        raise ValueError(f"wrong relation degree in {path}")
    matrix = tuple(tuple(ZZ(value) for value in row) for row in payload["rref_basis"])
    if len(matrix) != 22 or any(len(row) != 70 for row in matrix):
        raise ValueError(f"expected a 22 by 70 RREF matrix in {path}")
    return {
        "path": path,
        "payload": payload,
        "prime": ZZ(payload["prime"]),
        "matrix": matrix,
        "pivots": tuple(map(ZZ, payload["pivot_columns"])),
        "monomials": tuple(
            tuple(map(ZZ, exponents))
            for exponents in payload["monomial_exponents"]
        ),
    }


crt_records = tuple(relation_record(path) for path in input_paths)
validation_records = tuple(relation_record(path) for path in validation_paths)
all_records = crt_records+validation_records
reference = all_records[0]
for record in all_records[1:]:
    if record["pivots"] != reference["pivots"]:
        raise ValueError(f"pivot pattern differs in {record['path']}")
    if record["monomials"] != reference["monomials"]:
        raise ValueError(f"monomial order differs in {record['path']}")
if len(set(record["prime"] for record in all_records)) != len(all_records):
    raise ValueError("duplicate prime across CRT/validation artifacts")

primes = tuple(record["prime"] for record in crt_records)
modulus = prod(primes)


def reconstruct(row, column):
    residue = ZZ(
        CRT_list(
            [record["matrix"][row][column] for record in crt_records],
            list(primes),
        )
    )
    return QQ(residue.rational_reconstruction(modulus))


qq_matrix = tuple(
    tuple(reconstruct(row, column) for column in range(70))
    for row in range(22)
)
qq_ring = PolynomialRing(QQ, names=("D", "P", "Q", "E"), order="degrevlex")
D, P, Q, E = qq_ring.gens()
monomials = tuple(
    qq_ring.monomial(*exponents) for exponents in reference["monomials"]
)
relations = tuple(
    sum(
        (coefficient*monomial for coefficient, monomial in zip(row, monomials)),
        qq_ring.zero(),
    )
    for row in qq_matrix
)

# Verify canonical RREF exactly over QQ.
matrix = Matrix(QQ, qq_matrix)
if matrix.rref() != matrix:
    raise ArithmeticError("reconstructed relation matrix is not RREF")
if tuple(matrix.pivots()) != reference["pivots"]:
    raise ArithmeticError("reconstructed RREF has the wrong pivots")

# Verify every CRT and withheld finite-field artifact coefficientwise.
for record in all_records:
    finite = GF(record["prime"])
    for row in range(22):
        for column in range(70):
            coefficient = qq_matrix[row][column]
            reduced = finite(coefficient.numerator())/finite(coefficient.denominator())
            if reduced != finite(record["matrix"][row][column]):
                raise AssertionError(
                    f"QQ RREF misses {record['path']} at ({row},{column})"
                )

# Check all relations on the exact centered QQ jet.
series_path = Path(arguments.qq_series)
series_payload = json.loads(series_path.read_text())
if series_payload.get("schema") != "q80-cm24-qq-surface-series-v1":
    raise ValueError("unexpected exact-series schema")
order = ZZ(series_payload["order"])
centers = {
    "D": -QQ(1)/2,
    "P": QQ(9)/4,
    "Q": -QQ(9)/4,
    "E": -QQ(27)/32,
}
series = {}
for name in ("D", "P", "Q", "E"):
    values = [QQ(value) for value in series_payload["series"][name]][:order]
    values[0] -= centers[name]
    series[name] = tuple(values)


def multiply(left, right):
    result = [QQ.zero()]*order
    for first_index, first_value in enumerate(left):
        for second_index, second_value in enumerate(right[:order-first_index]):
            result[first_index+second_index] += first_value*second_value
    return tuple(result)


powers = {}
for name in ("D", "P", "Q", "E"):
    powers[(name, 0)] = (QQ.one(),)+(QQ.zero(),)*(order-1)
    for exponent in range(1, 5):
        powers[(name, exponent)] = multiply(
            powers[(name, exponent-1)], series[name]
        )
for relation_index, relation in enumerate(relations):
    residual = [QQ.zero()]*order
    for exponents, coefficient in relation.dict().items():
        term = powers[("D", exponents[0])]
        for name, exponent in zip(("P", "Q", "E"), exponents[1:]):
            term = multiply(term, powers[(name, exponent)])
        for index, value in enumerate(term):
            residual[index] += coefficient*value
    if any(residual):
        first_nonzero = next(index for index, value in enumerate(residual) if value)
        raise AssertionError(
            f"relation {relation_index} misses exact QQ jet at order {first_nonzero}"
        )


def metadata(record):
    path = record["path"]
    return {
        "path": str(path),
        "prime": int(record["prime"]),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


output_payload = {
    "schema": "q80-cm24-qq-quartic-ideal-v1",
    "status": "cross_prime_crt_reconstruction_with_withheld_and_exact_jet_checks",
    "scope": "candidate_global_quartic_ideal_pending_exact_parameter_substitution",
    "slope": "8/87",
    "centered_coordinates": {
        "D": "d+1/2",
        "P": "p-9/4",
        "Q": "q+9/4",
        "E": "e+27/32",
    },
    "relation_degree": 4,
    "matrix_dimensions": [22, 70],
    "pivot_columns": [int(value) for value in reference["pivots"]],
    "monomial_exponents": [list(map(int, value)) for value in reference["monomials"]],
    "rref_basis": [[str(value) for value in row] for row in qq_matrix],
    "relations": [str(value) for value in relations],
    "crt_modulus_bits": int(modulus.nbits()),
    "crt_inputs": [metadata(record) for record in crt_records],
    "withheld_inputs": [metadata(record) for record in validation_records],
    "exact_series": {"path": str(series_path), "sha256": hashlib.sha256(series_path.read_bytes()).hexdigest()},
    "exact_series_orders_checked": int(order),
}
output_path = Path(arguments.output)
output_path.parent.mkdir(parents=True, exist_ok=True)
output_path.write_text(json.dumps(output_payload, indent=2, sort_keys=True)+"\n")

largest_bits = max(
    max(abs(value.numerator()).nbits(), value.denominator().nbits())
    for row in qq_matrix
    for value in row
)
print(
    "Q80SURFACEIDEALCRT|slope=8/87|relations=22|monomials=70|"
    f"crt_primes={len(crt_records)}|withheld_primes={len(validation_records)}|"
    f"modulus_bits={modulus.nbits()}|largest_numden_bits={largest_bits}|"
    f"exact_orders={order}|output={output_path}|"
    "status=PASS_CROSS_PRIME_EXACT_QUARTIC_RECONSTRUCTION",
    flush=True,
)
