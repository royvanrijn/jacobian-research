#!/usr/bin/env sage -python
"""Lift q4/orbit1584 degree-one all-node sections from p=131 to QQ.

Use the resolved chart

    Z=t-z,  X=sum(x_i*t^i),  Y=L*sum(q_i*t^i),

where L is the product of the four finite reducible-fibre supports.  The
Weierstrass coefficient equations together with the four exact node equations
have a full-rank 11-variable Jacobian on each selected modular branch.  Newton
lifting, rational reconstruction, and literal QQ substitution certify the
sections.  No elimination or Groebner basis is used.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, Qp, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q4o1584-compact-weierstrass-qq.json"
MODULAR = LOCAL / "q4o1584-degree1-all-node-sections-mod131.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--precision", type=int, default=240)
parser.add_argument(
    "--output",
    default="artifacts/local/elkies-k3/q4o1584-degree1-all-node-sections-qq.json",
)
args = parser.parse_args()
OUTPUT = Path(args.output)
if not OUTPUT.is_absolute():
    OUTPUT = ROOT / OUTPUT

started = time.monotonic()
prime = ZZ(131)
precision = int(args.precision)
assert precision >= 80


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficient_bits(value):
    value = QQ(value)
    return max(abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())


compact = json.loads(MODEL.read_text())
modular = json.loads(MODULAR.read_text())
assert compact["status"] == "PASS_EXACT_QQ_Q4O1584_COMPACT_WEIERSTRASS_NORMALIZATION"
assert modular["status"] == "PASS_MOD131_Q4O1584_DEGREE1_ALL_NODE_SECTION_SEARCH"

RTQ = PolynomialRing(QQ, "t")
tq = RTQ.gen()
A_QQ = RTQ([QQ(value) for value in compact["compact_model"]["A_coefficients_low_to_high"]])
B_QQ = RTQ([QQ(value) for value in compact["compact_model"]["B_coefficients_low_to_high"]])
records = compact["compact_model"]["finite_reducible_fibres"]
supports_QQ = [QQ(record["support"]) for record in records]
nodes_QQ = []
for support in supports_QQ:
    node = -3 * B_QQ(support) / (2 * A_QQ(support))
    assert A_QQ(support) == -3 * node**2 and B_QQ(support) == 2 * node**3
    nodes_QQ.append(node)
L_QQ = RTQ.prod(tq - support for support in supports_QQ)


def reduce_mod_prime(value, field):
    value = QQ(value)
    return field(value.numerator()) / field(value.denominator())


F = GF(prime)
RTF = PolynomialRing(F, "t")
tf = RTF.gen()
A_F = RTF([reduce_mod_prime(value, F) for value in A_QQ])
B_F = RTF([reduce_mod_prime(value, F) for value in B_QQ])
supports_F = [reduce_mod_prime(value, F) for value in supports_QQ]
nodes_F = [reduce_mod_prime(value, F) for value in nodes_QQ]
L_F = RTF([reduce_mod_prime(value, F) for value in L_QQ])

K = Qp(prime, prec=precision, type="capped-rel")
RT = PolynomialRing(K, "t")
t = RT.gen()
A = RT([K(value) for value in A_QQ])
B = RT([K(value) for value in B_QQ])
supports = [K(value) for value in supports_QQ]
nodes = [K(value) for value in nodes_QQ]
L = RT([K(value) for value in L_QQ])


def unpack(values, ring):
    z = values[0]
    X = ring(list(values[1:7]))
    Q = ring(list(values[7:11]))
    Z = ring.gen() - z
    Y = ring(L.list() if ring is RT else (L_F.list() if ring is RTF else L_QQ.list())) * Q
    return z, Z, X, Q, Y


def residual(values):
    z, Z, X, Q, Y = unpack(values, RT)
    equation = Y**2 - X**3 - A * X * Z**4 - B * Z**6
    answers = [equation[index] for index in range(16)]
    answers.extend(
        X(support) - node * (support - z)**2
        for support, node in zip(supports, nodes)
    )
    return vector(K, answers)


def jacobian(values, ring, surface_A, surface_B, support_values, node_values, support_polynomial):
    z = values[0]
    X = ring(list(values[1:7]))
    Q = ring(list(values[7:11]))
    variable = ring.gen()
    Z = variable - z
    Y = support_polynomial * Q
    derivatives = [
        4 * surface_A * X * Z**3 + 6 * surface_B * Z**5,
    ]
    derivatives.extend(
        -(3 * X**2 + surface_A * Z**4) * variable**degree
        for degree in range(6)
    )
    derivatives.extend(
        2 * Y * support_polynomial * variable**degree
        for degree in range(4)
    )
    zero = ring.base_ring().zero()
    rows = [[
        derivative[degree] if degree <= derivative.degree() else zero
        for derivative in derivatives
    ] for degree in range(16)]
    for support, node in zip(support_values, node_values):
        rows.append(
            [2 * node * (support - z)]
            + [support**degree for degree in range(6)]
            + [zero] * 4
        )
    return matrix(ring.base_ring(), rows)


def minimum_valuation(values):
    nonzero = [value.valuation() for value in values if value]
    return min(nonzero) if nonzero else precision


def rational_reconstruct(value, digits):
    modulus = prime**digits
    residue = ZZ(value.lift()) % modulus
    return QQ(residue.rational_reconstruction(modulus))


# Keep one sign for each X,Z pair; negation gives the other section exactly.
branches = {}
for record in modular["search"]["sections"]:
    key = (tuple(record["Z_coefficients_low_to_high"]), tuple(record["X_coefficients_low_to_high"]))
    branches.setdefault(key, record)
branches = [branches[key] for key in sorted(branches)]
assert len(branches) == 4


def seed_values(record):
    z = F(record["z"])
    X = RTF(record["X_coefficients_low_to_high"])
    Q = RTF(record["Q_coefficients_low_to_high"])
    values = [z] + X.list() + [F.zero()] * (6 - len(X.list()))
    values += Q.list() + [F.zero()] * (4 - len(Q.list()))
    assert len(values) == 11
    _, Z, X_check, Q_check, Y = unpack(values, RTF)
    assert X_check == X and Q_check == Q
    assert Y**2 == X**3 + A_F * X * Z**4 + B_F * Z**6
    return vector(F, values)


def lift_one(index, record):
    seed = seed_values(record)
    J_F = jacobian(seed, RTF, A_F, B_F, supports_F, nodes_F, L_F)
    rank = int(J_F.rank())
    assert rank == 11
    pivot_rows = list(map(int, J_F.transpose().pivots()))
    assert len(pivot_rows) == 11
    determinant = int(matrix(F, [J_F.row(row) for row in pivot_rows]).det())
    assert determinant

    values = vector(K, [K(value).add_bigoh(1) for value in seed])
    known_precision = 1
    iterations = []
    while known_precision < precision:
        working_precision = min(2 * known_precision, precision)
        values = vector(K, [K(value.lift()).add_bigoh(working_precision) for value in values])
        full_residual = residual(values)
        chosen = vector(K, [full_residual[row] for row in pivot_rows])
        J = jacobian(values, RT, A, B, supports, nodes, L)
        square = matrix(K, [J.row(row) for row in pivot_rows])
        correction = square.solve_right(-chosen)
        values += correction
        iterations.append({
            "working_precision_p_adic_digits": working_precision,
            "minimum_full_residual_valuation_after": int(minimum_valuation(residual(values))),
            "minimum_correction_valuation": int(minimum_valuation(correction)),
        })
        known_precision = working_precision

    reconstruction_digits = precision - 12
    reconstructed = [rational_reconstruct(value, reconstruction_digits) for value in values]
    assert [reduce_mod_prime(value, F) for value in reconstructed] == list(seed)
    z, Z, X, Q, Y = unpack(reconstructed, RTQ)
    assert Y**2 == X**3 + A_QQ * X * Z**4 + B_QQ * Z**6
    assert all(
        X(support) == node * (support - z)**2 and Y(support) == 0
        for support, node in zip(supports_QQ, nodes_QQ)
    )
    assert z not in supports_QQ
    return {
        "branch_index": index,
        "mod131_augmented_jacobian_rank": rank,
        "selected_independent_equation_rows": pivot_rows,
        "selected_jacobian_determinant_mod131": determinant,
        "iterations": iterations,
        "z": str(z),
        "Z_coefficients_low_to_high": [str(value) for value in Z.list()],
        "X_coefficients_low_to_high": [str(value) for value in X.list()],
        "Q_coefficients_low_to_high": [str(value) for value in Q.list()],
        "Y_coefficients_low_to_high": [str(value) for value in Y.list()],
        "maximum_z_X_Q_rational_bits": max(map(coefficient_bits, reconstructed)),
        "maximum_Y_rational_bits": max(map(coefficient_bits, Y)),
        "exact_section_identity": True,
        "negative_section_also_certified": True,
    }


lifted = []
for index, record in enumerate(branches):
    result = lift_one(index, record)
    lifted.append(result)
    print(
        f"Q4O1584D1LIFT|branch={index}|z={result['z']}|bits="
        f"{result['maximum_z_X_Q_rational_bits']}/{result['maximum_Y_rational_bits']}|"
        f"elapsed={time.monotonic()-started:.3f}",
        flush=True,
    )

payload = {
    "schema": "elkies-k3.q4o1584-degree1-all-node-sections-qq.v1",
    "status": "PASS_EXACT_QQ_Q4O1584_FOUR_DEGREE1_ALL_NODE_SECTION_PAIRS",
    "prime": int(prime),
    "resolved_hensel": {
        "working_precision_p_adic_digits": precision,
        "variables": ["z", "six X coefficients", "four Q coefficients with Y=LQ"],
        "weierstrass_coefficient_equations": 16,
        "exact_node_incidence_equations": 4,
        "sign_pairs": len(lifted),
        "sections_including_negatives": 2 * len(lifted),
        "sections": lifted,
    },
    "coefficient_growth": {
        "maximum_z_X_Q_rational_bits": max(record["maximum_z_X_Q_rational_bits"] for record in lifted),
        "maximum_Y_rational_bits": max(record["maximum_Y_rational_bits"] for record in lifted),
    },
    "method": {
        "large_Groebner_required": False,
        "runtime_seconds": time.monotonic() - started,
    },
    "proof_boundary": (
        "Each regular modular branch was Newton lifted in a full-rank resolved chart, "
        "rationally reconstructed, reduced coefficientwise back to its seed, and verified "
        "by literal substitution over QQ. Lattice-class identification and transport to "
        "q4/orbit164 remain separate required checks."
    ),
    "next_required": (
        "Restrict the certified q4/orbit164 base function to these curves, select the two "
        "degree-one images, and match them to the target one-node lattice classes."
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
    "Q4O1584D1LIFT|pairs={}|sections={}|bits={}/{}|status={}|output={}".format(
        len(lifted), 2 * len(lifted),
        payload["coefficient_growth"]["maximum_z_X_Q_rational_bits"],
        payload["coefficient_growth"]["maximum_Y_rational_bits"],
        payload["status"], OUTPUT,
    ),
    flush=True,
)
