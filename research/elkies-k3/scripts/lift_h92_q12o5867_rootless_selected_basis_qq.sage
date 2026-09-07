#!/usr/bin/env sage -python
"""Lift the selected 17 regular rootless mod-131 sections to QQ.

status: ACTIVE_PROOF
claim: exact QQ sections for the selected terminal short-vector basis
inputs: q12o5867-smooth-rr-qq.json, selected modular basis certificate
outputs: q12o5867-rootless-selected-basis-qq.json

Each 13-equation/12-variable branch is lifted independently by ordinary
Newton--Hensel iteration using a full-rank 12-row minor.  Adaptive rational
reconstruction is accepted only after literal QQ substitution and exact
reduction to the supplied seed.  No Groebner basis or elimination is used.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, Qp, ZZ, matrix, vector
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q12o5867-smooth-rr-qq.json"
SEEDS = LOCAL / "q12o5867-rootless-mod131-selected-basis.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--reconstruction-start", type=int, default=2048)
parser.add_argument("--maximum-precision", type=int, default=16384)
parser.add_argument("--limit", type=int, default=17)
parser.add_argument("--output", type=Path, default=LOCAL / "q12o5867-rootless-selected-basis-qq.json")
args = parser.parse_args()
OUTPUT = args.output if args.output.is_absolute() else ROOT / args.output
assert 1 <= args.limit <= 17
assert 64 <= args.reconstruction_start <= args.maximum_precision
started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficient_bits(value):
    value = QQ(value)
    return max(abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())


model = json.loads(MODEL.read_text())
seeds = json.loads(SEEDS.read_text())
assert model["status"] == "PASS_EXACT_QQ_Q12O5867_SMOOTH_RR_ROOTLESS_JACOBIAN"
assert seeds["status"] == "PASS_MOD131_Q12O5867_ROOTLESS_REGULAR_SHORT_BASIS"
prime = ZZ(seeds["prime"])
assert prime == 131

RQ = PolynomialRing(QQ, "u")
uq = RQ.gen()
A_QQ = RQ([QQ(value) for value in model["child"]["minimal_A_coefficients_low_to_high"]])
B_QQ = RQ([QQ(value) for value in model["child"]["minimal_B_coefficients_low_to_high"]])
F = GF(prime)
RF = PolynomialRing(F, "u")


def reduce_qq(value, field=F):
    value = QQ(value)
    return field(value.numerator())/field(value.denominator())


A_F = RF([reduce_qq(value) for value in A_QQ])
B_F = RF([reduce_qq(value) for value in B_QQ])
K = Qp(prime, prec=args.maximum_precision, type="capped-rel")
RT = PolynomialRing(K, "u")
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


def minimum_valuation(values, fallback):
    nonzero = [int(value.valuation()) for value in values if value]
    return min(nonzero) if nonzero else int(fallback)


def reconstruct_vector(values, usable_precision):
    modulus = prime**usable_precision
    answer = []
    for value in values:
        if not value:
            answer.append(QQ.zero())
            continue
        residue = ZZ(value.lift()) % modulus
        answer.append(QQ(residue.rational_reconstruction(modulus)))
    return answer


def lift_seed(seed_record):
    seed = vector(F, seed_record["x_coefficients_low_to_high"] + seed_record["y_coefficients_low_to_high"])
    X_F, Y_F = split(seed, RF)
    assert Y_F**2 == X_F**3+A_F*X_F+B_F
    J_F = jacobian(seed, RF, A_F)
    rank = int(J_F.rank())
    assert rank == 12
    pivot_rows = list(map(int, J_F.transpose().pivots()))
    minor = int(matrix(F, [J_F.row(index) for index in pivot_rows]).det())
    assert minor

    values = vector(K, [K(value).add_bigoh(1) for value in seed])
    precision = 1
    iterations = []
    attempts = []
    reconstructed = None
    while precision < args.maximum_precision:
        working = min(2*precision, args.maximum_precision)
        values = vector(K, [K(value.lift()).add_bigoh(working) for value in values])
        full = residual(values)
        chosen = vector(K, [full[index] for index in pivot_rows])
        square = matrix(K, [jacobian(values, RT, A).row(index) for index in pivot_rows])
        correction = square.solve_right(-chosen)
        values += correction
        after = residual(values)
        iterations.append({
            "working_precision_p_adic_digits": working,
            "minimum_full_residual_valuation_after": minimum_valuation(after, working),
            "minimum_correction_valuation": minimum_valuation(correction, working),
        })
        precision = working
        if working < args.reconstruction_start:
            continue
        usable = working-8
        attempt = {
            "working_precision_p_adic_digits": working,
            "usable_precision_p_adic_digits": usable,
            "rational_reconstruction_succeeded": False,
            "literal_weierstrass_identity": False,
        }
        try:
            candidate = reconstruct_vector(values, usable)
            attempt["rational_reconstruction_succeeded"] = True
            attempt["maximum_candidate_rational_bits"] = max(map(coefficient_bits, candidate))
            X_exact, Y_exact = split(candidate, RQ)
            exact_identity = Y_exact**2 == X_exact**3+A_QQ*X_exact+B_QQ
            exact_reduction = [reduce_qq(value) for value in candidate] == list(seed)
            attempt["literal_weierstrass_identity"] = bool(exact_identity)
            attempt["exact_reduction_to_mod131_seed"] = bool(exact_reduction)
            if exact_identity and exact_reduction:
                reconstructed = candidate
        except (ArithmeticError, ValueError, ZeroDivisionError) as error:
            attempt["failure"] = type(error).__name__
        attempts.append(attempt)
        if reconstructed is not None:
            break
    if reconstructed is None:
        raise ArithmeticError(
            f"basis section {seed_record['basis_index']} did not reconstruct by {args.maximum_precision} p-adic digits"
        )
    X_exact, Y_exact = split(reconstructed, RQ)
    assert X_exact.degree() <= 4 and Y_exact.degree() <= 6
    assert Y_exact**2 == X_exact**3+A_QQ*X_exact+B_QQ
    assert [reduce_qq(value) for value in reconstructed] == list(seed)
    return {
        "basis_index": seed_record["basis_index"],
        "shell_record_index": seed_record["shell_record_index"],
        "section": {
            "x_coefficients_low_to_high": [str(value) for value in X_exact.list()],
            "y_coefficients_low_to_high": [str(value) for value in Y_exact.list()],
            "degrees_x_y": [int(X_exact.degree()), int(Y_exact.degree())],
            "P_dot_O": 0,
            "maximum_rational_bits": max(map(coefficient_bits, reconstructed)),
            "exact_weierstrass_identity": True,
            "exact_reduction_to_mod131_seed": True,
        },
        "hensel": {
            "prime": int(prime),
            "coefficient_equations": 13,
            "variables": 12,
            "mod131_jacobian_rank": rank,
            "selected_independent_equation_rows": pivot_rows,
            "selected_minor_determinant_mod131": minor,
            "successful_working_precision_p_adic_digits": precision,
            "iterations": iterations,
            "rational_reconstruction_attempts": attempts,
        },
    }


lifted = []
for seed_record in seeds["selected_sections"][:args.limit]:
    section_started = time.monotonic()
    record = lift_seed(seed_record)
    record["runtime_seconds"] = time.monotonic()-section_started
    lifted.append(record)
    print(
        "Q12O5867ROOTLESSLIFT|basis={}|shell={}|minor={}|precision={}|bits={}|seconds={:.3f}|"
        "status=PASS_EXACT_QQ_SECTION".format(
            record["basis_index"], record["shell_record_index"],
            record["hensel"]["selected_minor_determinant_mod131"],
            record["hensel"]["successful_working_precision_p_adic_digits"],
            record["section"]["maximum_rational_bits"], record["runtime_seconds"],
        ), flush=True,
    )

complete = len(lifted) == 17
payload = {
    "schema": "elkies-k3.h92-q12o5867-rootless-selected-basis-qq.v1",
    "status": (
        "PASS_EXACT_QQ_Q12O5867_ROOTLESS_17_SELECTED_SECTIONS"
        if complete else "PARTIAL_EXACT_QQ_Q12O5867_ROOTLESS_SELECTED_SECTIONS"
    ),
    "surface": {
        "equation": "y^2=x^3+A(u)*x+B(u)",
        "rootless_model": str(MODEL.relative_to(ROOT)),
    },
    "lifted_section_count": len(lifted),
    "selected_mod131_height_gram": seeds["height_gram"],
    "selected_mod131_height_gram_determinant": seeds["height_gram_determinant"],
    "short_basis_to_pinned_basis": seeds["short_basis_to_pinned_basis"],
    "sections": lifted,
    "method": {
        "construction": "independent regular-branch Newton--Hensel lifts with adaptive rational reconstruction",
        "reconstruction_start_p_adic_digits": args.reconstruction_start,
        "maximum_precision_p_adic_digits": args.maximum_precision,
        "large_Groebner_required": False,
        "elimination_required": False,
        "sage_version": SAGE_VERSION,
        "runtime_seconds": time.monotonic()-started,
    },
    "proof_boundary": (
        "Every displayed section satisfies the exact terminal QQ equation and reduces to its "
        "selected regular mod-131 seed. Exact characteristic-zero pairings and saturation are "
        "certified separately."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (MODEL, SEEDS)],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in (MODEL, SEEDS)},
    },
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q12O5867ROOTLESSLIFT|sections={}|status={}|seconds={:.3f}|output={}".format(
        len(lifted), payload["status"], payload["method"]["runtime_seconds"], OUTPUT,
    ), flush=True,
)
