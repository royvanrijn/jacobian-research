#!/usr/bin/env sage -python
"""Lift regular zero-node q4/orbit164 modular sections to QQ.

Use the full-rank twelve-variable polynomial-section chart, Newton lift each
modular branch, rationally reconstruct, and verify literal substitution.  No
Groebner basis is used.  Independence and marking remain separate gates.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, Qp, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q4o164-compact-weierstrass-qq.json"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--modular", type=Path, default=LOCAL / "q4o164-zero-node-sections-mod23.json")
parser.add_argument("--output", type=Path, default=LOCAL / "q4o164-zero-node-sections-from23-qq.json")
parser.add_argument("--precision", type=int, default=180)
args = parser.parse_args()
MODULAR = args.modular if args.modular.is_absolute() else ROOT / args.modular
OUTPUT = args.output if args.output.is_absolute() else ROOT / args.output
PRECISION = int(args.precision)

started = time.monotonic()
model = json.loads(MODEL.read_text())
modular = json.loads(MODULAR.read_text())
assert model["status"] == "PASS_EXACT_QQ_Q4O164_COMPACT_WEIERSTRASS_NORMALIZATION"
assert modular["status"] == "PASS_MODP_Q4O164_ZERO_NODE_SECTION_SCAN"
PRIME = ZZ(modular["prime"])

RQ = PolynomialRing(QQ, "t")
A_QQ = RQ([QQ(value) for value in model["compact_model"]["A_coefficients_low_to_high"]])
B_QQ = RQ([QQ(value) for value in model["compact_model"]["B_coefficients_low_to_high"]])


def reduce_qq(value, field):
    value = QQ(value)
    return field(value.numerator()) / field(value.denominator())


F = GF(PRIME)
RF = PolynomialRing(F, "t")
A_F = RF([reduce_qq(value, F) for value in A_QQ])
B_F = RF([reduce_qq(value, F) for value in B_QQ])
K = Qp(PRIME, prec=PRECISION, type="capped-rel")
RT = PolynomialRing(K, "t")
A = RT([K(value) for value in A_QQ])
B = RT([K(value) for value in B_QQ])


def polynomials(values, ring):
    return ring(list(values[:5])), ring(list(values[5:]))


def residual(values):
    X, Y = polynomials(values, RT)
    equation = Y**2 - X**3 - A * X - B
    return vector(K, [equation[index] for index in range(13)])


def jacobian(values, ring, surface_A):
    X, Y = polynomials(values, ring)
    dx, dy = -3 * X**2 - surface_A, 2 * Y
    zero = ring.base_ring().zero()
    return matrix(ring.base_ring(), [[
        dx[degree-shift] if 0 <= degree-shift <= dx.degree() else zero
        for shift in range(5)
    ] + [
        dy[degree-shift] if 0 <= degree-shift <= dy.degree() else zero
        for shift in range(7)
    ] for degree in range(13)])


def minimum_valuation(values):
    nonzero = [value.valuation() for value in values if value]
    return min(nonzero) if nonzero else PRECISION


def reconstruct(value):
    modulus = PRIME ** (PRECISION - 12)
    return QQ((ZZ(value.lift()) % modulus).rational_reconstruction(modulus))


def bits(value):
    value = QQ(value)
    return max(abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())


lifted = []
for modular_index, row in enumerate(modular["sections"]):
    seed_values = row["x_coefficients_low_to_high"] + row["y_coefficients_low_to_high"]
    seed = vector(F, seed_values + [0] * (12-len(seed_values)))
    XF, YF = polynomials(seed, RF)
    assert YF**2 == XF**3 + A_F * XF + B_F
    JF = jacobian(seed, RF, A_F)
    assert JF.rank() == 12
    pivot_rows = list(map(int, JF.transpose().pivots()))
    determinant = int(matrix(F, [JF.row(index) for index in pivot_rows]).det())
    assert determinant
    values = vector(K, [K(value).add_bigoh(1) for value in seed])
    known = 1
    iterations = []
    while known < PRECISION:
        working = min(2 * known, PRECISION)
        values = vector(K, [K(value.lift()).add_bigoh(working) for value in values])
        full = residual(values)
        square = matrix(K, [jacobian(values, RT, A).row(index) for index in pivot_rows])
        correction = square.solve_right(-vector(K, [full[index] for index in pivot_rows]))
        values += correction
        iterations.append({
            "working_precision_p_adic_digits": working,
            "minimum_full_residual_valuation_after": int(minimum_valuation(residual(values))),
        })
        known = working
    try:
        reconstructed = [reconstruct(value) for value in values]
    except ArithmeticError:
        continue
    X, Y = polynomials(reconstructed, RQ)
    if Y**2 != X**3 + A_QQ * X + B_QQ:
        continue
    lifted.append({
        "modular_candidate_index": modular_index,
        "selected_independent_equation_rows": pivot_rows,
        "selected_jacobian_determinant_modp": determinant,
        "iterations": iterations,
        "x_coefficients_low_to_high": [str(value) for value in X.list()],
        "y_coefficients_low_to_high": [str(value) for value in Y.list()],
        "maximum_x_rational_bits": max(map(bits, X)),
        "maximum_y_rational_bits": max(map(bits, Y)),
        "exact_compact_weierstrass_identity": True,
    })

payload = {
    "schema": "elkies-k3.q4o164-zero-node-sections-qq.v1",
    "status": "PASS_EXACT_QQ_Q4O164_ZERO_NODE_SECTION_LIFTS",
    "prime": int(PRIME),
    "precision_p_adic_digits": PRECISION,
    "exact_rational_lifts": lifted,
    "method": {"large_Groebner_required": False, "runtime_seconds": time.monotonic()-started},
    "proof_boundary": (
        "Each retained branch is an exact QQ polynomial section. MW independence and "
        "marked-lattice identification remain separate gates."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (MODEL, MODULAR)],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (MODEL, MODULAR)
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4O164ZERONODELIFTS|prime={}|modular={}|exact={}|status={}|output={}".format(
        PRIME, len(modular["sections"]), len(lifted), payload["status"], OUTPUT,
    ),
    flush=True,
)
