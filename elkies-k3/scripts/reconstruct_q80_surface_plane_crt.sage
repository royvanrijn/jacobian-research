#!/usr/bin/env sage
"""CRT-reconstruct the slope-8/87 q80 ``(D,Q)`` plane equation.

The inputs are canonical centered relation spaces produced at split primes by
``extend_q80_cm24_branch_modp.sage``.  For each prime we form the four-variable
candidate ideal, eliminate ``P,E``, and normalize the unique ``(D,Q)``
eliminant to have coefficient one at ``D^10``.  Coefficients reconstructed
from the declared CRT inputs are checked at every reconstruction prime, at
separate withheld validation primes, and against the available exact QQ jet.

This certifies a cross-prime finite-jet equation.  It is not a global identity
until a rational parametrization or direct substitution into the marked
surface family proves it.
"""

import argparse
import glob
import hashlib
import json
from pathlib import Path

from sage.all import CRT_list, GF, PolynomialRing, QQ, ZZ, prod


parser = argparse.ArgumentParser()
parser.add_argument("--input-glob", action="append", required=True)
parser.add_argument("--validation-glob", action="append", default=[])
parser.add_argument(
    "--qq-series",
    default="artifacts/generated-results/q80-cm24-slope-8-87-qq-order28-series.json",
)
parser.add_argument(
    "--output",
    default="artifacts/generated-results/q80-cm24-slope-8-87-qq-DQ-plane.json",
)
arguments = parser.parse_args()


def expand(patterns):
    return tuple(sorted({Path(name) for pattern in patterns for name in glob.glob(pattern)}))


input_paths = expand(arguments.input_glob)
validation_paths = expand(arguments.validation_glob)
if not input_paths:
    raise ValueError("no CRT relation artifacts matched")
if set(input_paths) & set(validation_paths):
    raise ValueError("CRT and withheld validation inputs must be disjoint")


def plane_record(path):
    payload = json.loads(path.read_text())
    if payload.get("schema") != "q80-cm24-split-prime-formal-branch-v1":
        raise ValueError(f"unexpected schema in {path}")
    if payload.get("kind") != "canonical_centered_relation_space":
        raise ValueError(f"expected a relation space in {path}")
    if payload.get("characteristic_zero_slope") != "8/87":
        raise ValueError(f"wrong branch in {path}")
    if tuple(payload.get("centered_variables", ())) != ("D", "P", "Q", "E"):
        raise ValueError(f"wrong centered coordinates in {path}")
    prime = ZZ(payload["prime"])
    finite = GF(prime)
    ring = PolynomialRing(finite, names=("D", "P", "Q", "E"), order="degrevlex")
    D, P, Q, E = ring.gens()
    monomials = tuple(
        ring.monomial(*map(ZZ, exponents))
        for exponents in payload["monomial_exponents"]
    )
    relations = tuple(
        sum(
            (finite(coefficient)*monomial
             for coefficient, monomial in zip(row, monomials)),
            ring.zero(),
        )
        for row in payload["rref_basis"]
    )
    ideal = ring.ideal(relations)
    if ideal.dimension() != 1:
        raise ValueError(f"relation ideal in {path} does not define a curve")
    elimination = ideal.elimination_ideal((P, E))
    plane_ring = PolynomialRing(finite, names=("D", "Q"), order="degrevlex")
    generators = tuple(plane_ring(str(value)) for value in elimination.gens())
    if len(generators) != 1:
        raise ValueError(f"expected one plane eliminant in {path}")
    polynomial = generators[0]
    leading = polynomial[(10, 0)]
    if not leading:
        raise ValueError(f"plane eliminant in {path} lacks D^10")
    polynomial /= leading
    if polynomial.total_degree() != 10 or polynomial.degree(plane_ring.gen(1)) != 5:
        raise ValueError(f"unexpected plane bidegree in {path}")
    return {
        "path": path,
        "payload": payload,
        "prime": prime,
        "polynomial": polynomial,
    }


crt_records = tuple(plane_record(path) for path in input_paths)
validation_records = tuple(plane_record(path) for path in validation_paths)
all_records = crt_records+validation_records
primes = tuple(record["prime"] for record in crt_records)
if len(set(record["prime"] for record in all_records)) != len(all_records):
    raise ValueError("duplicate prime across CRT/validation artifacts")
modulus = prod(primes)
support = tuple(sorted({
    exponents
    for record in all_records
    for exponents in record["polynomial"].dict()
}))


def reconstruct(exponents):
    residues = [ZZ(record["polynomial"][exponents]) for record in crt_records]
    residue = ZZ(CRT_list(residues, list(primes)))
    return QQ(residue.rational_reconstruction(modulus))


coefficients = {exponents: reconstruct(exponents) for exponents in support}
qq_ring = PolynomialRing(QQ, names=("D", "Q"), order="degrevlex")
D, Q = qq_ring.gens()
qq_polynomial = sum(
    (coefficient*D**exponents[0]*Q**exponents[1]
     for exponents, coefficient in coefficients.items()),
    qq_ring.zero(),
)
assert qq_polynomial[(10, 0)] == 1


for record in all_records:
    finite = GF(record["prime"])
    reduced_ring = record["polynomial"].parent()
    reduced = reduced_ring.zero()
    for exponents, coefficient in qq_polynomial.dict().items():
        reduced += (
            finite(coefficient.numerator())/finite(coefficient.denominator())
        )*reduced_ring.monomial(*exponents)
    if reduced != record["polynomial"]:
        raise AssertionError(f"QQ plane equation misses {record['path']}")


series_path = Path(arguments.qq_series)
series_payload = json.loads(series_path.read_text())
if series_payload.get("schema") != "q80-cm24-qq-surface-series-v1":
    raise ValueError("unexpected exact-series schema")
order = ZZ(series_payload["order"])
centers = {"D": -QQ(1)/2, "Q": -QQ(9)/4}
series = {}
for name in ("D", "Q"):
    values = [QQ(value) for value in series_payload["series"][name]][:order]
    values[0] -= centers[name]
    series[name] = tuple(values)


def multiply(left, right):
    result = [QQ(0)]*order
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right[:order-i]):
            result[i+j] += left_value*right_value
    return tuple(result)


powers = {}
for name in ("D", "Q"):
    powers[(name, 0)] = (QQ(1),)+(QQ(0),)*(order-1)
    for exponent in range(1, 11):
        powers[(name, exponent)] = multiply(powers[(name, exponent-1)], series[name])
residual = [QQ(0)]*order
for exponents, coefficient in qq_polynomial.dict().items():
    term = multiply(powers[("D", exponents[0])], powers[("Q", exponents[1])])
    for index, value in enumerate(term):
        residual[index] += coefficient*value
if any(residual):
    raise AssertionError("reconstructed plane equation misses the exact QQ jet")


def metadata(record):
    path = record["path"]
    return {
        "path": str(path),
        "prime": int(record["prime"]),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


output = Path(arguments.output)
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps({
    "schema": "q80-cm24-qq-DQ-plane-v1",
    "status": "cross_prime_crt_reconstruction_with_withheld_and_exact_jet_checks",
    "scope": "candidate_global_plane_pending_rational_parametrization_or_direct_substitution",
    "slope": "8/87",
    "centered_coordinates": {"D": "d+1/2", "Q": "q+9/4"},
    "normalization": "coefficient(D^10)=1",
    "total_degree": int(qq_polynomial.total_degree()),
    "bidegree": [int(qq_polynomial.degree(D)), int(qq_polynomial.degree(Q))],
    "crt_modulus_bits": int(modulus.nbits()),
    "crt_inputs": [metadata(record) for record in crt_records],
    "withheld_inputs": [metadata(record) for record in validation_records],
    "exact_series": str(series_path),
    "exact_series_orders_checked": int(order),
    "polynomial": str(qq_polynomial),
    "terms": [
        {"exponents": list(exponents), "coefficient": str(coefficient)}
        for exponents, coefficient in sorted(qq_polynomial.dict().items())
    ],
}, indent=2, sort_keys=True)+"\n")
largest_bits = max(
    max(abs(value.numerator()).nbits(), value.denominator().nbits())
    for value in coefficients.values()
)
print(
    "Q80SURFACEPLANECRT|slope=8/87|coordinates=D,Q|"
    f"crt_primes={len(crt_records)}|withheld_primes={len(validation_records)}|"
    f"modulus_bits={modulus.nbits()}|terms={len(qq_polynomial.dict())}|"
    f"largest_numden_bits={largest_bits}|exact_orders={order}|output={output}|"
    "status=PASS_CROSS_PRIME_EXACT_PLANE_RECONSTRUCTION",
    flush=True,
)
