#!/usr/bin/env sage -python
"""Construct the actual H92 E7 complete quotient for the q=6 third section.

The third transported section has E7 vertical correction

    c=(22,44,66,44,33,33,55)=11*v(U),

where U=x-c2*t^2-c3*t^3 is the actual translated H92 coordinate from the
first q=6 E7 quotient.  Using these translated coordinates, rather than raw
old x, the valuation-defined complete ideal has the predicted colength 363.
This is the finite resolved quotient needed to impose the third marked
divisor's E7 condition on a degree-44 ambient.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import ZZ, GF, matrix, vector, PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
Q6_BLOCK = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-p1-actual-e7-quotient-block.json"
TARGET = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-third-e7-local-target.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-third-actual-e7-quotient-block.json"
Q6_BLOCK_SHA256 = "3796ee20121a94ce6d3a707c0cd119983b64fce79336fa99a5a894729174900c"
TARGET_SHA256 = "a0699e4ec75930cc93a9706ddf96f4ffc744954809e53a24029bd8c6668843f7"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--q6-block", type=Path, default=Q6_BLOCK)
parser.add_argument("--target", type=Path, default=TARGET)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

if args.q6_block == Q6_BLOCK:
    assert digest(args.q6_block) == Q6_BLOCK_SHA256
if args.target == TARGET:
    assert digest(args.target) == TARGET_SHA256
q6_block = json.loads(args.q6_block.read_text())
target = json.loads(args.target.read_text())
assert q6_block["status"] == "PASS_EXACT_Q6_P1_ACTUAL_E7_QUOTIENT_BLOCK"
assert target["status"] == "PASS_EXACT_Q6_THIRD_E7_LATTICE_TARGET"

v_t = vector(ZZ, q6_block["actual_exceptional_orders"]["T"])
v_u = vector(ZZ, q6_block["actual_exceptional_orders"]["U"])
v_y = vector(ZZ, q6_block["actual_exceptional_orders"]["Y"])
cycle = -vector(ZZ, target["resolved_exceptional_coefficients"])
assert cycle == 11*v_u == vector(ZZ, (22, 44, 66, 44, 33, 33, 55))
cartan = matrix(ZZ, [
    [2, 0, 0, -1, 0, 0, 0],
    [0, 2, 0, 0, -1, 0, -1],
    [0, 0, 2, -1, 0, -1, -1],
    [-1, 0, -1, 2, 0, 0, 0],
    [0, -1, 0, 0, 2, 0, 0],
    [0, 0, -1, 0, 0, 2, 0],
    [0, -1, -1, 0, 0, 0, 2],
])
assert cartan*cycle == vector(ZZ, (0, 0, 0, 0, 22, 0, 0))
expected_colength = ZZ(cycle*cartan*cycle)//2
assert expected_colength == 363


def in_complete_ideal(exponents):
    i, a, b = exponents
    return min(i*v_t+a*v_u+b*v_y-cycle) >= 0


# The target is reached by T^33, U^11, and T^30*Y, so this bounding box
# contains all standard monomials and all minimal generators.
candidates = [
    (i, a, b)
    for b in range(2)
    for a in range(12)
    for i in range(34)
    if in_complete_ideal((i, a, b))
]
minimal_generators = [
    exponent for exponent in candidates
    if not any(
        predecessor != exponent
        and all(left <= right for left, right in zip(predecessor, exponent))
        for predecessor in candidates
    )
]
assert len(minimal_generators) == 23
quotient_exponents = [
    (i, a, b)
    for b in range(2)
    for a in range(11)
    for i in range(33)
    if not in_complete_ideal((i, a, b))
]
assert len(quotient_exponents) == expected_colength

# Reuse the exact translated H92 relation already certified by the q=6 block.
c2 = QQ(q6_block["translated_coordinates"]["c2"])
c3 = QQ(q6_block["translated_coordinates"]["c3"])
# The original coefficient values are not needed to determine the monomial
# quotient, but relation reduction makes the finite block an actual H92 local
# quotient rather than a valuation count.
from importlib.machinery import SourceFileLoader
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
anchor = SourceFileLoader("h92_q6_third_actual_e7_quotient_anchor", str(ANCHOR)).load_module()
r, s = anchor.EXPECTED_H92
_, formulas = anchor.parse_h92(H92)
A1, A, B1, B, B2 = (QQ(value(r, s)) for value in formulas)
# Put Y first in lexicographic order, so the monic Weierstrass relation has
# leading monomial Y^2.  This is the normal form b in {0,1} used above.
ring = PolynomialRing(QQ, names=("Y", "U", "T"), order="lex")
Y, U, T = ring.gens()
x_branch = c2*T**2+c3*T**3
relation = Y**2-(x_branch+U)**3-(A1*T**3+A*T**4)*(x_branch+U)-(
    B1*T**5+B*T**6+B2*T**7
)
generators = [T**i*U**a*Y**b for i, a, b in minimal_generators]
monomial_ideal = ring.ideal(generators)
quotient_basis = tuple(T**i*U**a*Y**b for i, a, b in quotient_exponents)
# The monic relation in Y^2 and the monomial generators span the quotient by
# these 363 normal monomials over QQ.  Their independence is certified after
# reduction modulo 1009.  A QQ relation, cleared to p-primitive integral
# coefficients, would give a nonzero relation in this same finite quotient.
prime = 1009
finite_ring = PolynomialRing(GF(prime), names=("Y", "U", "T"), order="lex")
finite_ideal = finite_ring.ideal((finite_ring(relation), *(finite_ring(value) for value in generators)))
assert finite_ideal.vector_space_dimension() == expected_colength

payload = {
    "schema": "elkies-k3.h92-q6-third-actual-e7-quotient-block.v1",
    "status": "PASS_EXACT_Q6_THIRD_ACTUAL_E7_QUOTIENT_BLOCK",
    "inputs": {
        "q6_actual_e7_block": {"path": str(args.q6_block.relative_to(ROOT)), "sha256": digest(args.q6_block)},
        "third_lattice_target": {"path": str(args.target.relative_to(ROOT)), "sha256": digest(args.target)},
    },
    "translated_coordinates": q6_block["translated_coordinates"],
    "cycle": [int(value) for value in cycle],
    "cartan_boundary": [int(value) for value in cartan*cycle],
    "complete_ideal_generators": [[int(value) for value in exponent] for exponent in minimal_generators],
    "generator_interpretation": "[i,a,b] means T^i*U^a*Y^b, with b in {0,1}",
    "quotient_basis_exponents": [[int(value) for value in exponent] for exponent in quotient_exponents],
    "quotient_dimension": int(expected_colength),
    "actual_H92_relation": "included and checked in the quotient ideal",
    "dimension_certificate": "the 363 normal monomials span over QQ; their independence is certified by the same quotient dimension modulo 1009",
    "compiler_instruction": "Evaluate the degree-44 ambient and marked-chord DAG in this quotient, then use resolved_chart_quotient_condition to form the third-section E7 block.",
    "boundary": "This selects and certifies the finite actual H92 quotient. Evaluating every high-degree ambient generator and the marked chord in it remains the next step.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6THIRDACTUALE7QUOTIENT|generators=23|length=363|"
    "status=PASS_EXACT_Q6_THIRD_ACTUAL_E7_QUOTIENT_BLOCK",
    flush=True,
)
