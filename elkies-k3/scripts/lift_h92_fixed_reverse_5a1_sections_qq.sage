#!/usr/bin/env sage
"""Lift regular q52 Abel-word section seeds from p=167 toward QQ.

The coefficient chart fixes Z monic and augments the Weierstrass identity by
the exact I2-node incidences recorded by the modular compiler.  The selected
square Jacobian is nonsingular mod 167.  A low-precision run is also useful as
a coefficient-growth probe; it writes a restartable residue checkpoint when
rational reconstruction is not yet complete.  No elimination or Groebner
basis is used.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, Qp, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "fixed-reverse-4a1-rr-qq.json"
POINTING = LOCAL / "fixed-reverse-4a1-pointing-qq.json"
SEED = LOCAL / "fixed-reverse-5a1-abel-word-seeds-mod167.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--precision", type=int, default=256)
parser.add_argument("--section-index", type=int, default=0)
parser.add_argument("--seed", type=Path, help="optional prior residue checkpoint")
parser.add_argument("--output", type=Path)
args = parser.parse_args()

prime = ZZ(167)
precision = int(args.precision)
section_index = int(args.section_index)
assert precision >= 32 and 0 <= section_index < 4
started = time.monotonic()


def read_json(path):
    return json.loads(path.read_text())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


model = read_json(MODEL)
pointing = read_json(POINTING)
seed_payload = read_json(SEED)
assert model["status"] == "PASS_EXACT_QQ_FIXED_REVERSE_4A1_RR_JACOBIAN"
assert pointing["status"] == "PASS_EXACT_QQ_FIXED_REVERSE_4A1_POINTING"
assert seed_payload["status"] == "PASS_MODP_FIXED_REVERSE_5A1_ABEL_WORD_REGULAR_SEEDS"
assert seed_payload["prime"] == prime
seed_record = seed_payload["sections"][section_index]
pole = int(seed_record["P_dot_O"])
variable_count = 6 * pole + 12
assert seed_record["resolved_jacobian_rank"] == variable_count

RTQ = PolynomialRing(QQ, "t")
tq = RTQ.gen()
A_QQ = RTQ([QQ(value) for value in model["child"]["minimal_A_coefficients_low_to_high"]])
B_QQ = RTQ([QQ(value) for value in model["child"]["minimal_B_coefficients_low_to_high"]])

# Two supports are explicit in the pointing certificate and the remaining
# quadratic splits exactly over QQ.  Match all four to their mod-167 residues.
supports_QQ = [
    QQ(record["child_I2_support"])
    for record in pointing["effective_horizontal_components"]
]
remaining_factor = RTQ([
    QQ(value) for value in
    pointing["remaining_vertical_components"][
        "child_I2_support_factor_coefficients_low_to_high"
    ]
])
supports_QQ.extend(remaining_factor.roots(QQ, multiplicities=False))
assert len(supports_QQ) == 4
nodes_QQ = [
    -3 * B_QQ(support) / (2 * A_QQ(support)) for support in supports_QQ
]


def reduce_mod_prime(value, field):
    value = QQ(value)
    return field(value.numerator()) / field(value.denominator())


F = GF(prime)
support_by_residue = {
    int(reduce_mod_prime(support, F)): (support, node)
    for support, node in zip(supports_QQ, nodes_QQ)
}
node_hits_QQ = []
for support_residue, node_residue in seed_record["I2_node_hits_support_and_node"]:
    support, node = support_by_residue[int(support_residue)]
    assert int(reduce_mod_prime(node, F)) == int(node_residue)
    node_hits_QQ.append((support, node))

K = Qp(prime, prec=precision, type="capped-rel")
RT = PolynomialRing(K, "t")
t = RT.gen()
A = RT([K(value) for value in A_QQ])
B = RT([K(value) for value in B_QQ])
node_hits = [(K(support), K(node)) for support, node in node_hits_QQ]

x_count = 2 * pole + 5
y_count = 3 * pole + 7
z_count = pole


def unpack(values, ring):
    X = ring(list(values[:x_count]))
    Y = ring(list(values[x_count:x_count + y_count]))
    Z = ring(list(values[x_count + y_count:]) + [ring.base_ring().one()])
    return X, Y, Z


def residual(values):
    X, Y, Z = unpack(values, RT)
    equation = Y ** 2 - X ** 3 - A * X * Z ** 4 - B * Z ** 6
    answer = [equation[index] for index in range(6 * pole + 13)]
    for support, node in node_hits:
        answer.extend((X(support) - node * Z(support) ** 2, Y(support)))
    return vector(K, answer)


def jacobian(values):
    X, Y, Z = unpack(values, RT)
    derivatives = []
    derivatives.extend(
        -(3 * X ** 2 + A * Z ** 4) * t ** degree
        for degree in range(x_count)
    )
    derivatives.extend(2 * Y * t ** degree for degree in range(y_count))
    derivatives.extend(
        -(4 * A * X * Z ** 3 + 6 * B * Z ** 5) * t ** degree
        for degree in range(z_count)
    )
    zero = K.zero()
    rows = [[
        derivative[degree] if degree <= derivative.degree() else zero
        for derivative in derivatives
    ] for degree in range(6 * pole + 13)]
    for support, node in node_hits:
        rows.append(
            [support ** degree for degree in range(x_count)]
            + [zero] * y_count
            + [
                -2 * node * Z(support) * support ** degree
                for degree in range(z_count)
            ]
        )
        rows.append(
            [zero] * x_count
            + [support ** degree for degree in range(y_count)]
            + [zero] * z_count
        )
    return matrix(K, rows)


def minimum_valuation(values):
    nonzero = [value.valuation() for value in values if value]
    return min(nonzero) if nonzero else precision


seed_values = (
    seed_record["X_coefficients_low_to_high"]
    + seed_record["Y_coefficients_low_to_high"]
    + seed_record["Z_coefficients_low_to_high"][:-1]
)
assert len(seed_values) == variable_count
known_precision = 1
if args.seed:
    prior = read_json(args.seed.resolve())
    assert prior["schema"] == "elkies-k3.fixed-reverse-5a1-section-hensel.v1"
    assert prior["prime"] == prime and prior["section_index"] == section_index
    seed_values = [ZZ(value) for value in prior["hensel"]["residues"]]
    known_precision = int(prior["hensel"]["precision"])
    assert precision > known_precision

values = vector(K, [K(value).add_bigoh(known_precision) for value in seed_values])
pivot_rows = list(map(int, seed_record["selected_independent_equation_rows"]))
assert len(pivot_rows) == variable_count
iterations = []
while known_precision < precision:
    working_precision = min(2 * known_precision, precision)
    values = vector(K, [K(value.lift()).add_bigoh(working_precision) for value in values])
    full_residual = residual(values)
    chosen = vector(K, [full_residual[row] for row in pivot_rows])
    square = matrix(K, [jacobian(values).row(row) for row in pivot_rows])
    correction = square.solve_right(-chosen)
    values += correction
    iterations.append({
        "precision": working_precision,
        "minimum_residual_valuation": int(minimum_valuation(residual(values))),
        "minimum_correction_valuation": int(minimum_valuation(correction)),
    })
    known_precision = working_precision
    print(
        "FIXEDREVERSE5A1HENSEL|section={}|precision={}|residual={}|seconds={:.3f}".format(
            section_index, known_precision,
            iterations[-1]["minimum_residual_valuation"],
            time.monotonic() - started,
        ),
        flush=True,
    )

modulus = prime ** precision
residues = [ZZ(value.lift()) % modulus for value in values]
reconstructed = []
for residue in residues:
    try:
        reconstructed.append(QQ(residue.rational_reconstruction(modulus)))
    except (ArithmeticError, ValueError):
        reconstructed.append(None)

exact = all(value is not None for value in reconstructed)
exact_record = None
if exact:
    X, Y, Z = unpack(reconstructed, RTQ)
    exact = Y ** 2 == X ** 3 + A_QQ * X * Z ** 4 + B_QQ * Z ** 6
    exact = exact and all(
        X(support) == node * Z(support) ** 2 and Y(support) == 0
        for support, node in node_hits_QQ
    )
    if exact:
        exact_record = {
            "X_coefficients_low_to_high": list(map(str, X.list())),
            "Y_coefficients_low_to_high": list(map(str, Y.list())),
            "Z_coefficients_low_to_high": list(map(str, Z.list())),
            "maximum_rational_bits": max(
                max(abs(value.numerator()).nbits(), value.denominator().nbits())
                for value in reconstructed
            ),
        }

payload = {
    "schema": "elkies-k3.fixed-reverse-5a1-section-hensel.v1",
    "status": (
        "PASS_EXACT_QQ_FIXED_REVERSE_5A1_SECTION"
        if exact else "PASS_REGULAR_PADIC_FIXED_REVERSE_5A1_SECTION_CHECKPOINT"
    ),
    "prime": int(prime),
    "section_index": section_index,
    "P_dot_O": pole,
    "resolved_node_hit_count": len(node_hits),
    "hensel": {
        "precision": precision,
        "residues": list(map(str, residues)),
        "iterations": iterations,
        "rational_reconstruction_count": sum(value is not None for value in reconstructed),
        "variable_count": variable_count,
    },
    "exact_section": exact_record,
    "method": {
        "resolved_full_rank_chart": True,
        "groebner_or_elimination": False,
        "runtime_seconds": time.monotonic() - started,
    },
    "inputs": {
        str(path.relative_to(ROOT)): sha256(path)
        for path in (MODEL, POINTING, SEED)
    },
}
default_output = LOCAL / "fixed-reverse-5a1-section{}-hensel-p167-prec{}.json".format(
    section_index, precision
)
output = args.output or default_output
output = output if output.is_absolute() else ROOT / output
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "FIXEDREVERSE5A1HENSEL|section={}|precision={}|reconstructed={}/{}|status={}|output={}".format(
        section_index, precision,
        payload["hensel"]["rational_reconstruction_count"], variable_count,
        payload["status"], output,
    ),
    flush=True,
)
