#!/usr/bin/env sage -python
"""Lift the mod-131 degree-one q12/o5867 compiler branch to QQ.

The modular constructor supplies a P.O=0 section of the exact q8/o376 child
with polynomial degrees (4,6).  Its 13 coefficient equations have a rank-12
Jacobian in the twelve section coefficients.  Newton lift that regular branch,
rationally reconstruct it, and verify the Weierstrass equation literally.
No elimination or Groebner basis is used.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, Qp, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
Q8 = LOCAL / "q4o164-q8o376-smooth-rr-qq.json"
MODULAR = LOCAL / "q12o5867-degree1-compiler-branch-mod131.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--precision", type=int, default=900)
parser.add_argument(
    "--output",
    default="artifacts/local/elkies-k3/q12o5867-degree1-compiler-branch-qq.json",
)
args = parser.parse_args()
OUTPUT = Path(args.output)
if not OUTPUT.is_absolute():
    OUTPUT = ROOT / OUTPUT
precision = int(args.precision)
assert precision >= 100


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficient_bits(value):
    value = QQ(value)
    return max(abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())


started = time.monotonic()
q8 = json.loads(Q8.read_text())
modular = json.loads(MODULAR.read_text())
assert q8["status"] == "PASS_EXACT_QQ_Q8O376_4A1_RR_JACOBIAN_AND_P1229_ZERO"
assert modular["status"] == "PASS_MOD131_Q12O5867_DEGREE1_COMPILER_BRANCH_ON_Q8_CHILD"
prime = ZZ(modular["prime"])
assert prime == 131

RQ = PolynomialRing(QQ, "v")
vq = RQ.gen()
A_QQ = RQ([QQ(value) for value in q8["child"]["minimal_A_coefficients_low_to_high"]])
B_QQ = RQ([QQ(value) for value in q8["child"]["minimal_B_coefficients_low_to_high"]])


def reduce_qq(value, field):
    value = QQ(value)
    return field(value.numerator())/field(value.denominator())


F = GF(prime)
RF = PolynomialRing(F, "v")
vf = RF.gen()
A_F = RF([reduce_qq(value, F) for value in A_QQ])
B_F = RF([reduce_qq(value, F) for value in B_QQ])
seed_record = modular["q8_child_section_mod131"]
seed = vector(F, seed_record["x_coefficients_low_to_high"]+seed_record["y_coefficients_low_to_high"])
assert len(seed) == 12

K = Qp(prime, prec=precision, type="capped-rel")
RT = PolynomialRing(K, "v")
v = RT.gen()
A = RT([K(value) for value in A_QQ])
B = RT([K(value) for value in B_QQ])


def split(values, ring):
    return ring(list(values[:5])), ring(list(values[5:]))


def residual(values):
    X, Y = split(values, RT)
    equation = Y**2-X**3-A*X-B
    return vector(K, [equation[index] for index in range(13)])


def jacobian(values, ring, surface_A):
    X, Y = split(values, ring)
    dx = -3*X**2-surface_A
    dy = 2*Y
    zero = ring.base_ring().zero()
    return matrix(ring.base_ring(), [[
        dx[degree-shift] if 0 <= degree-shift <= dx.degree() else zero
        for shift in range(5)
    ]+[
        dy[degree-shift] if 0 <= degree-shift <= dy.degree() else zero
        for shift in range(7)
    ] for degree in range(13)])


X_F, Y_F = split(seed, RF)
assert Y_F**2 == X_F**3+A_F*X_F+B_F
J_F = jacobian(seed, RF, A_F)
rank = int(J_F.rank())
assert rank == 12
pivot_rows = list(map(int, J_F.transpose().pivots()))
assert len(pivot_rows) == 12
minor_determinant = int(matrix(F, [J_F.row(index) for index in pivot_rows]).det())
assert minor_determinant == 50


def minimum_valuation(values):
    nonzero = [value.valuation() for value in values if value]
    return min(nonzero) if nonzero else precision


values = vector(K, [K(value).add_bigoh(1) for value in seed])
known_precision = 1
iterations = []
while known_precision < precision:
    working_precision = min(2*known_precision, precision)
    values = vector(K, [K(value.lift()).add_bigoh(working_precision) for value in values])
    full = residual(values)
    chosen = vector(K, [full[index] for index in pivot_rows])
    J = jacobian(values, RT, A)
    square = matrix(K, [J.row(index) for index in pivot_rows])
    correction = square.solve_right(-chosen)
    values += correction
    iterations.append({
        "working_precision_p_adic_digits": working_precision,
        "minimum_full_residual_valuation_after": int(minimum_valuation(residual(values))),
        "minimum_correction_valuation": int(minimum_valuation(correction)),
    })
    known_precision = working_precision


def rational_reconstruct(value):
    if not value:
        return QQ.zero()
    usable_precision = min(precision-20, int(value.precision_absolute())-8)
    assert usable_precision >= 32
    modulus = prime**usable_precision
    residue = ZZ(value.lift()) % modulus
    return QQ(residue.rational_reconstruction(modulus))


reconstructed = [rational_reconstruct(value) for value in values]
assert [reduce_qq(value, F) for value in reconstructed] == list(seed)
X_QQ, Y_QQ = split(reconstructed, RQ)
assert Y_QQ**2 == X_QQ**3+A_QQ*X_QQ+B_QQ
assert X_QQ.degree() <= 4 and Y_QQ.degree() <= 6

payload = {
    "schema": "elkies-k3.h92-q12o5867-degree1-compiler-branch-qq.v1",
    "status": "PASS_EXACT_QQ_Q12O5867_DEGREE1_COMPILER_BRANCH_ON_Q8_CHILD",
    "compiler_branch": modular["compiler_branch"],
    "exact_lattice_to_equation_word": modular["exact_lattice_to_equation_word"],
    "section": {
        "x_coefficients_low_to_high": [str(value) for value in X_QQ.list()],
        "y_coefficients_low_to_high": [str(value) for value in Y_QQ.list()],
        "degrees_x_y": [int(X_QQ.degree()), int(Y_QQ.degree())],
        "P_dot_O": 0,
        "component_profile_at_good_reduction": modular["compiler_branch"]["current_4A1_component_pairings"],
        "maximum_rational_bits": max(map(coefficient_bits, reconstructed)),
        "exact_weierstrass_identity": True,
        "exact_reduction_to_mod131_seed": True,
    },
    "hensel": {
        "prime": int(prime),
        "working_precision_p_adic_digits": precision,
        "coefficient_equations": 13,
        "variables": 12,
        "mod131_jacobian_rank": rank,
        "selected_independent_equation_rows": pivot_rows,
        "selected_minor_determinant_mod131": minor_determinant,
        "iterations": iterations,
    },
    "method": {
        "large_Groebner_required": False,
        "elimination_required": False,
        "runtime_seconds": time.monotonic()-started,
    },
    "proof_boundary": (
        "The promoted degree-one q4/o164 parent branch is now an exact QQ P.O=0 section on the q8/o376 child, "
        "with literal Weierstrass substitution and exact reduction to the marked mod-131 seed. The other three "
        "q12/o5867 compiler sections and the q12 resolved RR pencil remain open."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (Q8, MODULAR)],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in (Q8, MODULAR)},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q12O5867DEG1QQ|h0_section_system=13x12|rank={}|minor={}|degrees={},{}|bits={}|"
    "status={}|output={}".format(
        rank, minor_determinant, X_QQ.degree(), Y_QQ.degree(),
        payload["section"]["maximum_rational_bits"], payload["status"], OUTPUT,
    ),
    flush=True,
)
