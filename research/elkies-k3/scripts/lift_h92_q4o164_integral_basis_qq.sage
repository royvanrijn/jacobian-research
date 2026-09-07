#!/usr/bin/env sage -python
"""Lift a rank-nine q4/orbit164 integral-section basis from p=131 to QQ.

The naive polynomial-section Jacobian drops rank when a section passes through
a reducible-fibre node.  Following the repository's earlier resolved-Hensel
constructions, add the exact node-incidence equations selected by the modular
component profile, choose twelve independent rows, and Newton lift that
regular chart.  Rational reconstruction and literal substitution finish the
QQ certificate.  No Groebner basis is used.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, Qp, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q4o164-compact-weierstrass-qq.json"
MODULAR = LOCAL / "q4o164-integral-sections-mod131.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--precision", type=int, default=160)
parser.add_argument(
    "--output",
    default="artifacts/local/elkies-k3/q4o164-integral-basis-qq.json",
)
args = parser.parse_args()
OUTPUT = Path(args.output)
if not OUTPUT.is_absolute():
    OUTPUT = ROOT / OUTPUT

started = time.monotonic()
prime = ZZ(131)
precision = int(args.precision)
assert precision >= 60


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficient_bits(value):
    value = QQ(value)
    return max(abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())


compact = json.loads(MODEL.read_text())
modular = json.loads(MODULAR.read_text())
assert compact["status"] == "PASS_EXACT_QQ_Q4O164_COMPACT_WEIERSTRASS_NORMALIZATION"
assert modular["status"] in {
    "PASS_MOD131_Q4O164_INTEGRAL_SECTION_SUBGROUP_RANK9",
    "PASS_MOD131_Q4O164_PAIR_NODE_SECTION_SUBGROUP_RANK8",
}

RTQ = PolynomialRing(QQ, "t")
tq = RTQ.gen()
KQ = RTQ.fraction_field()
A_QQ = RTQ([QQ(value) for value in compact["compact_model"]["A_coefficients_low_to_high"]])
B_QQ = RTQ([QQ(value) for value in compact["compact_model"]["B_coefficients_low_to_high"]])
supports_QQ = [
    None if record["support"] == "infinity" else QQ(record["support"])
    for record in compact["compact_model"]["reducible_fibres"]
]
orders = [2, 2, 4, 4]
nodes_QQ = []
for support in supports_QQ:
    if support is None:
        nodes_QQ.append(QQ(3))
    else:
        a_value, b_value = A_QQ(support), B_QQ(support)
        node = -3 * b_value / (2 * a_value)
        assert a_value == -3 * node**2 and b_value == 2 * node**3
        nodes_QQ.append(node)


def reduce_mod_prime(value, field):
    value = QQ(value)
    return field(value.numerator()) / field(value.denominator())


F = GF(prime)
RTF = PolynomialRing(F, "t")
tf = RTF.gen()
A_F = RTF([reduce_mod_prime(value, F) for value in A_QQ])
B_F = RTF([reduce_mod_prime(value, F) for value in B_QQ])
supports_F = [None if value is None else reduce_mod_prime(value, F) for value in supports_QQ]
nodes_F = [reduce_mod_prime(value, F) for value in nodes_QQ]
c_F = reduce_mod_prime(compact["exact_coordinate_change"]["c"], F)
s_F = reduce_mod_prime(compact["exact_coordinate_change"]["s"], F)
assert c_F and s_F


def compact_seed(record):
    x_old = RTF(record["x_coefficients_low_to_high"])
    y_old = RTF(record["y_coefficients_low_to_high"])
    x_new = x_old(c_F * tf) / s_F**2
    y_new = y_old(c_F * tf) / s_F**3
    assert x_new.degree() <= 4 and y_new.degree() <= 6
    assert y_new**2 == x_new**3 + A_F * x_new + B_F
    return list(x_new) + [F.zero()] * (5 - len(x_new.list())) + list(y_new) + [F.zero()] * (7 - len(y_new.list()))


K = Qp(prime, prec=precision, type="capped-rel")
RT = PolynomialRing(K, "t")
t = RT.gen()
A = RT([K(value) for value in A_QQ])
B = RT([K(value) for value in B_QQ])
supports = [None if value is None else K(value) for value in supports_QQ]
nodes = [K(value) for value in nodes_QQ]


def polynomials(values, ring):
    return ring(list(values[:5])), ring(list(values[5:]))


def residual(values, profile):
    X, Y = polynomials(values, RT)
    equation = Y**2 - X**3 - A * X - B
    answers = [equation[index] for index in range(13)]
    for index, label in enumerate(profile):
        if not label:
            continue
        if supports[index] is None:
            answers.extend((X[4] - nodes[index], Y[6]))
        else:
            answers.extend((X(supports[index]) - nodes[index], Y(supports[index])))
    return vector(K, answers)


def jacobian(values, profile, ring, surface_A, support_values):
    X, Y = polynomials(values, ring)
    x_derivative = -3 * X**2 - surface_A
    y_derivative = 2 * Y
    zero = ring.base_ring().zero()
    rows = [[
        x_derivative[degree - shift]
        if 0 <= degree - shift <= x_derivative.degree() else zero
        for shift in range(5)
    ] + [
        y_derivative[degree - shift]
        if 0 <= degree - shift <= y_derivative.degree() else zero
        for shift in range(7)
    ] for degree in range(13)]
    for index, label in enumerate(profile):
        if not label:
            continue
        support = support_values[index]
        if support is None:
            rows.append([zero] * 4 + [ring.base_ring().one()] + [zero] * 7)
            rows.append([zero] * 11 + [ring.base_ring().one()])
        else:
            rows.append([support**degree for degree in range(5)] + [zero] * 7)
            rows.append([zero] * 5 + [support**degree for degree in range(7)])
    return matrix(ring.base_ring(), rows)


def minimum_valuation(values):
    nonzero = [value.valuation() for value in values if value]
    return min(nonzero) if nonzero else precision


def rational_reconstruct(value, digits):
    modulus = prime**digits
    residue = ZZ(value.lift()) % modulus
    return QQ(residue.rational_reconstruction(modulus))


def lift_one(index, record):
    profile = tuple(map(int, record["component_profile"]))
    seed = compact_seed(record)
    X_F, Y_F = polynomials(seed, RTF)
    J_F = jacobian(seed, profile, RTF, A_F, supports_F)
    rank = int(J_F.rank())
    assert rank == 12
    pivot_rows = list(map(int, J_F.transpose().pivots()))
    assert len(pivot_rows) == 12
    determinant = int(matrix(F, [J_F.row(row) for row in pivot_rows]).det())
    assert determinant

    values = vector(K, [K(value).add_bigoh(1) for value in seed])
    known_precision = 1
    iterations = []
    while known_precision < precision:
        working_precision = min(2 * known_precision, precision)
        values = vector(K, [K(value.lift()).add_bigoh(working_precision) for value in values])
        all_residual = residual(values, profile)
        chosen = vector(K, [all_residual[row] for row in pivot_rows])
        J = jacobian(values, profile, RT, A, supports)
        square = matrix(K, [J.row(row) for row in pivot_rows])
        correction = square.solve_right(-chosen)
        values += correction
        iterations.append({
            "working_precision_p_adic_digits": working_precision,
            "minimum_full_residual_valuation_after": int(minimum_valuation(residual(values, profile))),
            "minimum_correction_valuation": int(minimum_valuation(correction)),
        })
        known_precision = working_precision

    reconstruction_digits = precision - 10
    reconstructed = [rational_reconstruct(value, reconstruction_digits) for value in values]
    assert [reduce_mod_prime(value, F) for value in reconstructed] == seed
    X_QQ, Y_QQ = polynomials(reconstructed, RTQ)
    assert Y_QQ**2 == X_QQ**3 + A_QQ * X_QQ + B_QQ
    for fibre_index, label in enumerate(profile):
        if not label:
            continue
        support = supports_QQ[fibre_index]
        if support is None:
            assert X_QQ[4] == nodes_QQ[fibre_index] and Y_QQ[6] == 0
        else:
            assert X_QQ(support) == nodes_QQ[fibre_index] and Y_QQ(support) == 0
    return {
        "basis_index": index,
        "component_profile": list(profile),
        "mod131_augmented_jacobian_rank": rank,
        "selected_independent_equation_rows": pivot_rows,
        "selected_jacobian_determinant_mod131": determinant,
        "iterations": iterations,
        "x_coefficients_low_to_high": [str(value) for value in X_QQ.list()],
        "y_coefficients_low_to_high": [str(value) for value in Y_QQ.list()],
        "maximum_x_rational_bits": max(map(coefficient_bits, X_QQ)),
        "maximum_y_rational_bits": max(map(coefficient_bits, Y_QQ)),
        "exact_section_identity": True,
    }, (X_QQ, Y_QQ)


basis_records = modular["integral_subgroup"]["basis"]
expected_rank = int(
    modular["integral_subgroup"]["independence_certificate"]["signature_rank"]
)
assert len(basis_records) == expected_rank
lifted = []
coordinates = []
for index, record in enumerate(basis_records):
    result, point = lift_one(index, record)
    lifted.append(result)
    coordinates.append(point)
    print(
        f"Q4O164LIFT|section={index}|bits={result['maximum_x_rational_bits']}/"
        f"{result['maximum_y_rational_bits']}|elapsed={time.monotonic()-started:.3f}",
        flush=True,
    )

independence = modular["integral_subgroup"]["independence_certificate"]
ell = int(independence["quotient_prime"])
signature = matrix(GF(ell), independence["basis_signature_matrix"])
assert signature.ncols() == expected_rank and signature.rank() == expected_rank

status = (
    "PASS_EXACT_QQ_Q4O164_INTEGRAL_SECTION_BASIS_RANK9"
    if expected_rank == 9 else
    "PASS_EXACT_QQ_Q4O164_PAIR_NODE_SECTION_SUBGROUP_RANK8"
)

payload = {
    "schema": "elkies-k3.q4o164-integral-basis-qq.v2",
    "status": status,
    "prime": int(prime),
    "compact_model": str(MODEL.relative_to(ROOT)),
    "resolved_hensel": {
        "working_precision_p_adic_digits": precision,
        "variables": 12,
        "weierstrass_coefficient_equations": 13,
        "node_incidence_equations_per_nonidentity_support": 2,
        "sections": lifted,
    },
    "exact_independence_certificate": {
        "argument": (
            "The exact QQ sections reduce coefficientwise to the selected modular "
            "branches. Their images in the displayed product of smooth-fibre quotients "
            "E(F_131)/(ell) have the displayed full column rank, so no primitive integral relation can "
            "exist among the QQ sections."
        ),
        "quotient_prime": ell,
        "smooth_fibre_supports": independence["smooth_fibre_supports"],
        "signature_matrix": independence["basis_signature_matrix"],
        "signature_rank": int(signature.rank()),
    },
    "coefficient_growth": {
        "maximum_x_rational_bits": max(record["maximum_x_rational_bits"] for record in lifted),
        "maximum_y_rational_bits": max(record["maximum_y_rational_bits"] for record in lifted),
    },
    "method": {
        "large_Groebner_required": False,
        "runtime_seconds": time.monotonic() - started,
    },
    "proof_boundary": (
        "The selected modular polynomial branches were lifted in full-rank resolved node-incidence "
        "charts, rationally reconstructed, and verified by exact QQ substitution. Exact "
        "coefficientwise reduction plus the smooth-specialization quotient matrix proves "
        "their independence. The pair-node shell currently has rank eight; a ninth direction "
        "must be recovered from a zero- or one-node integral section before identifying the "
        "promoted q8/orbit376 horizontal."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (MODEL, MODULAR)],
        "sha256": {
            str(path.relative_to(ROOT)): sha256(path) for path in (MODEL, MODULAR)
        },
    },
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4O164LIFT|rank={}|sections={}|bits={}/{}|status={}|output={}".format(
        expected_rank, len(lifted),
        payload["coefficient_growth"]["maximum_x_rational_bits"],
        payload["coefficient_growth"]["maximum_y_rational_bits"],
        payload["status"], OUTPUT,
    ),
    flush=True,
)
