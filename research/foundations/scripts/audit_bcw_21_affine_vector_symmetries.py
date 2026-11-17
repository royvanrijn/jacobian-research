#!/usr/bin/env python3
"""Exclude affine-vector-field translation symmetries of the 21D BCW map.

For the cubic homogeneous part H, an infinitesimal affine translation
V(x)=A*x+b preserves H only if JH(x)V(x)=0. Homogeneity separates this into
JH(x)b=0 and JH(x)A*x=0. The first space is already the constant Jacobian
kernel; this audit proves the coefficient system for A has full rank over Q.
"""

from fractions import Fraction
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "generated-results" / "essential_bcw_21_counterexample.json"
PRIME = 1_000_003


stored = json.loads(SOURCE.read_text())
dimension = stored["dimension"]
assert dimension == 21

# Unknown j*n+k is A[j,k]. Expand (d H_i / d x_j) A[j,k] x_k.
equations: dict[tuple[int, tuple[int, ...]], dict[int, Fraction]] = {}
for output_index, component in enumerate(stored["H"]):
    for term in component:
        coefficient = Fraction(term["coefficient"])
        exponents = [0] * dimension
        for variable_index, exponent in term["monomial"]:
            exponents[variable_index] = exponent
        for derivative_index, exponent in enumerate(exponents):
            if exponent == 0:
                continue
            derivative_coefficient = coefficient * exponent
            derivative_exponents = exponents[:]
            derivative_exponents[derivative_index] -= 1
            for multiplier_index in range(dimension):
                result_exponents = derivative_exponents[:]
                result_exponents[multiplier_index] += 1
                key = output_index, tuple(result_exponents)
                unknown = derivative_index * dimension + multiplier_index
                row = equations.setdefault(key, {})
                row[unknown] = row.get(unknown, Fraction(0)) + derivative_coefficient


def residue(value: Fraction) -> int:
    assert value.denominator % PRIME
    return value.numerator * pow(value.denominator, -1, PRIME) % PRIME


rows = [
    {column: residue(value) for column, value in row.items() if value}
    for row in equations.values()
]
rows = [row for row in rows if row]
pivots: dict[int, dict[int, int]] = {}
for source_row in rows:
    row = dict(source_row)
    while row:
        pivot = min(row)
        if pivot not in pivots:
            inverse = pow(row[pivot], -1, PRIME)
            pivots[pivot] = {
                column: value * inverse % PRIME
                for column, value in row.items()
                if value % PRIME
            }
            break
        factor = row[pivot]
        for column, value in pivots[pivot].items():
            reduced = (row.get(column, 0) - factor * value) % PRIME
            if reduced:
                row[column] = reduced
            else:
                row.pop(column, None)

unknowns = dimension * dimension
assert len(rows) == 2484
assert len(pivots) == unknowns

# Full rank modulo a good prime exhibits a nonzero maximal minor over Q.
print("PASS 21D affine symmetry audit: 2484 exact coefficient equations")
print("PASS 21D affine symmetry audit: rank 441 modulo 1000003")
print("PASS 21D affine symmetry audit: JH(x) A x = 0 forces A = 0 over Q")
print("PASS 21D affine symmetry audit: no nonzero affine translation symmetry")
