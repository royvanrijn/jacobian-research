#!/usr/bin/env sage -python
"""Lift q4/orbit323 simple-pole seeds in a resolved node chart.

Use Z=t-z, deg(X)<=6, deg(Y)<=9 and augment the 19 Weierstrass coefficient
equations by exact node equations at compact t=0,1,infinity.  These rows
restore full rank for the relevant p=59 branches despite cuspidal auxiliary
reduction.  Newton lifting and literal characteristic-zero verification use
only linear systems; no Groebner basis is used.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, Qp, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
COMPACT = LOCAL / "q4o208-compact-weierstrass-qq.json"
MODULAR = LOCAL / "q4o208-q4o323-horizontal-sums-mod59.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--precision", type=int, default=800)
parser.add_argument("--output", type=Path, default=LOCAL / "q4o208-q4o323-horizontal-resolved-qq.json")
args = parser.parse_args()
OUTPUT = args.output if args.output.is_absolute() else ROOT / args.output
PRECISION = int(args.precision)
PRIME = ZZ(59)
started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bits(value):
    value = QQ(value)
    return max(abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())


compact = json.loads(COMPACT.read_text())
modular = json.loads(MODULAR.read_text())
assert compact["status"] == "PASS_EXACT_QQ_Q4O208_COMPACT_WEIERSTRASS_NORMALIZATION"
assert modular["status"] == "PASS_MOD59_Q4O323_SIMPLE_POLE_SUM_SEEDS"

RQ = PolynomialRing(QQ, "t")
tq = RQ.gen()
A_QQ = RQ([QQ(value) for value in compact["compact_model"]["A_coefficients_low_to_high"]])
B_QQ = RQ([QQ(value) for value in compact["compact_model"]["B_coefficients_low_to_high"]])
node_zero_QQ = QQ(3)
node_one_QQ = -3*B_QQ(1)/(2*A_QQ(1))
node_infinity_QQ = -3*B_QQ[12]/(2*A_QQ[8])


def reduce_qq(value, field):
    value = QQ(value)
    return field(value.numerator())/field(value.denominator())


F = GF(PRIME)
RF = PolynomialRing(F, "t")
tf = RF.gen()
A_F = RF([reduce_qq(value, F) for value in A_QQ])
B_F = RF([reduce_qq(value, F) for value in B_QQ])
nodes_F = [reduce_qq(value, F) for value in (node_zero_QQ, node_one_QQ, node_infinity_QQ)]
KF = RF.fraction_field()


def finite_rational(record):
    return KF(RF(record["numerator_coefficients_low_to_high"]))/KF(
        RF(record["denominator_coefficients_low_to_high"])
    )


def seed_values(record):
    x_value = finite_rational(record["x"])
    y_value = finite_rational(record["y"])
    denominator = x_value.denominator().monic()
    factors = list(denominator.factor())
    assert len(factors) == 1 and factors[0][0].degree() == 1 and factors[0][1] == 2
    Z = factors[0][0].monic()
    z = -Z[0]
    X = RF(x_value*Z**2)
    Y = RF(y_value*Z**3)
    values = [z] + X.list() + [F.zero()]*(7-len(X.list()))
    values += Y.list() + [F.zero()]*(10-len(Y.list()))
    return vector(F, values)


# One y sign per z,X branch.
seeds = {}
for record in modular["horizontal_sum_seeds"]:
    seed = seed_values(record)
    seeds.setdefault(tuple(seed[:8]), (record, seed))
assert len(seeds) == 40

K = Qp(PRIME, prec=PRECISION, type="capped-rel")
RT = PolynomialRing(K, "t")
t = RT.gen()
A = RT([K(value) for value in A_QQ])
B = RT([K(value) for value in B_QQ])
nodes = [K(value) for value in (node_zero_QQ, node_one_QQ, node_infinity_QQ)]


def unpack(values, ring):
    z = values[0]
    Z = ring.gen()-z
    X = ring(list(values[1:8]))
    Y = ring(list(values[8:18]))
    return z, Z, X, Y


def residual(values):
    z, Z, X, Y = unpack(values, RT)
    equation = Y**2-X**3-A*X*Z**4-B*Z**6
    answer = [equation[degree] for degree in range(19)]
    for support, node in zip((K(0), K(1)), nodes[:2]):
        answer.extend((X(support)-node*(support-z)**2, Y(support)))
    answer.extend((X[6]-nodes[2], Y[9]))
    return vector(K, answer)


def jacobian(values, ring, surface_A, surface_B, node_values):
    z, Z, X, Y = unpack(values, ring)
    variable = ring.gen()
    derivatives = [4*surface_A*X*Z**3+6*surface_B*Z**5]
    derivatives.extend(-(3*X**2+surface_A*Z**4)*variable**degree for degree in range(7))
    derivatives.extend(2*Y*variable**degree for degree in range(10))
    zero = ring.base_ring().zero()
    rows = [[
        derivative[degree] if degree <= derivative.degree() else zero
        for derivative in derivatives
    ] for degree in range(19)]
    for support, node in zip((ring.base_ring()(0), ring.base_ring()(1)), node_values[:2]):
        rows.append([2*node*(support-z)]
                    + [support**degree for degree in range(7)] + [zero]*10)
        rows.append([zero]*8 + [support**degree for degree in range(10)])
    row = [zero]*18
    row[7] = ring.base_ring().one()
    rows.append(row)
    row = [zero]*18
    row[17] = ring.base_ring().one()
    rows.append(row)
    return matrix(ring.base_ring(), rows)


def minimum_valuation(values):
    nonzero = [value.valuation() for value in values if value]
    return min(nonzero) if nonzero else PRECISION


def reconstruct(value):
    digits = PRECISION-12
    modulus = PRIME**digits
    residue = ZZ(value.lift()) % modulus
    return QQ(residue.rational_reconstruction(modulus))


regular = []
qq_lifts = []
for branch_index, (source, seed) in enumerate(seeds.values()):
    J_F = jacobian(seed, RF, A_F, B_F, nodes_F)
    rank = int(J_F.rank())
    if rank != 18:
        continue
    pivot_rows = list(map(int, J_F.transpose().pivots()))
    determinant = int(matrix(F, [J_F.row(row) for row in pivot_rows]).det())
    values = vector(K, [K(value).add_bigoh(1) for value in seed])
    known_precision = 1
    iterations = []
    while known_precision < PRECISION:
        working_precision = min(2*known_precision, PRECISION)
        values = vector(K, [K(value.lift()).add_bigoh(working_precision) for value in values])
        full = residual(values)
        square = matrix(K, [jacobian(values, RT, A, B, nodes).row(row) for row in pivot_rows])
        correction = square.solve_right(-vector(K, [full[row] for row in pivot_rows]))
        values += correction
        iterations.append({
            "working_precision_p_adic_digits": working_precision,
            "minimum_augmented_residual_valuation_after": int(minimum_valuation(residual(values))),
        })
        known_precision = working_precision
    regular.append({
        "branch_index": branch_index,
        "source_one_node_index": int(source["one_node_index"]),
        "source_two_node_index": int(source["two_node_index"]),
        "seed_residues": [int(value) for value in seed],
        "pivot_rows": pivot_rows,
        "determinant_mod59": determinant,
        "iterations": iterations,
    })
    try:
        exact = [reconstruct(value) for value in values]
    except ArithmeticError:
        continue
    z_QQ, Z_QQ, X_QQ, Y_QQ = unpack(exact, RQ)
    if Y_QQ**2 != X_QQ**3+A_QQ*X_QQ*Z_QQ**4+B_QQ*Z_QQ**6:
        continue
    if not all(
        X_QQ(support) == node*(support-z_QQ)**2 and Y_QQ(support) == 0
        for support, node in ((QQ(0), node_zero_QQ), (QQ(1), node_one_QQ))
    ):
        continue
    if X_QQ[6] != node_infinity_QQ or Y_QQ[9] != 0:
        continue
    qq_lifts.append({
        "branch_index": branch_index,
        "z": str(z_QQ),
        "Z_coefficients_low_to_high": [str(value) for value in Z_QQ.list()],
        "X_coefficients_low_to_high": [str(value) for value in X_QQ.list()],
        "Y_coefficients_low_to_high": [str(value) for value in Y_QQ.list()],
        "maximum_rational_bits": max(map(bits, exact)),
        "exact_augmented_simple_pole_identity": True,
        "negative_section_also_certified": True,
    })

payload = {
    "schema": "elkies-k3.h3-q4o208-q4o323-horizontal-resolved-qq.v1",
    "status": (
        "PASS_EXACT_QQ_Q4O323_RESOLVED_SIMPLE_POLE_HORIZONTAL"
        if qq_lifts else "PASS_NO_QQ_Q4O323_AMONG_REGULAR_RESOLVED_MOD59_BRANCHES"
    ),
    "prime": int(PRIME),
    "precision_p_adic_digits": PRECISION,
    "search": {
        "distinct_mod59_x_z_seeds": len(seeds),
        "regular_augmented_branches": len(regular),
        "singular_augmented_branches": len(seeds)-len(regular),
        "exact_QQ_branches": len(qq_lifts),
    },
    "regular_branch_diagnostics": regular,
    "exact_QQ_horizontal_sections": qq_lifts,
    "method": {"large_Groebner_required": False, "runtime_seconds": time.monotonic()-started},
    "proof_boundary": (
        "All full-rank mod-59 branches in the six-node-equation simple-pole chart were "
        "Newton lifted. Listed QQ branches, if any, satisfy the Weierstrass and node equations "
        "literally. Algebraic non-QQ reconstruction and marked lattice matching remain separate gates."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (COMPACT, MODULAR)],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in (COMPACT, MODULAR)},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q4O323RESOLVED|seeds={}|regular={}|qq={}|status={}|output={}".format(
        len(seeds), len(regular), len(qq_lifts), payload["status"], OUTPUT,
    ), flush=True,
)
