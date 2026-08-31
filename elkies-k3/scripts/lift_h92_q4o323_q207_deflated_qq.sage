#!/usr/bin/env sage-python
"""Lift the two third-order q207 seeds by two-step deflated Hensel.

The target-shape mod-61 sections 5887 and 5903 have coefficient-Jacobian
rank 69 in 70 projective variables.  The first deflation adjoins a normalized
kernel vector; the second adjoins a normalized kernel vector for that system.
Their ranks are 69 -> 139 -> 280, so ordinary Newton--Hensel applies to the
final 295-equation, 280-variable system.  Only polynomial coefficient
arithmetic and linear solves are used; there is no elimination or Groebner
basis.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import (
    GF, PolynomialRing, QQ, Qp, ZZ, block_matrix, matrix, vector, zero_matrix,
)


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
POINTING = LOCAL / "q4o323-component2-pointing-qq.json"
CANDIDATES = LOCAL / "q4o323-q207-four-section-words-mod61.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--candidate", type=int, choices=(5887, 5903))
parser.add_argument("--reconstruction-start", type=int, default=64)
parser.add_argument("--maximum-precision", type=int, default=4096)
parser.add_argument(
    "--output", type=Path,
    default=LOCAL / "q4o323-q207-deflated-horizontal-qq.json",
)
args = parser.parse_args()
candidate_indices = [args.candidate] if args.candidate is not None else [5887, 5903]
started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficient_bits(value):
    value = QQ(value)
    return max(abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())


pointing = json.loads(POINTING.read_text())
candidates = json.loads(CANDIDATES.read_text())
assert pointing["status"] == "PASS_EXACT_QQ_Q4O323_OLD_A11_COMPONENT2_POINTING"
assert candidates["prime"] == 61
prime = ZZ(61)

RQ = PolynomialRing(QQ, "u")
A_QQ = RQ(pointing["global_short_model"]["A_coefficients_low_to_high"])
B_QQ = RQ(pointing["global_short_model"]["B_coefficients_low_to_high"])
F = GF(prime)
RF = PolynomialRing(F, "u")
A_F = RF([F(QQ(c).numerator()) / F(QQ(c).denominator()) for c in A_QQ])
B_F = RF([F(QQ(c).numerator()) / F(QQ(c).denominator()) for c in B_QQ])

maximum_precision = int(args.maximum_precision)
reconstruction_start = int(args.reconstruction_start)
assert 16 <= reconstruction_start <= maximum_precision
KP = Qp(prime, prec=maximum_precision, type="capped-rel")
RP = PolynomialRing(KP, "u")
A_P = RP([KP(c) for c in A_QQ])
B_P = RP([KP(c) for c in B_QQ])


def shifted_matrix(base, px, py, pz):
    zero = base.zero()
    return matrix(base, 73, 70, lambda row, column: (
        px[row-column]
        if column < 25 and 0 <= row-column <= px.degree() else
        py[row-(column-25)]
        if 25 <= column < 60 and 0 <= row-(column-25) <= py.degree() else
        pz[row-(column-60)]
        if column >= 60 and 0 <= row-(column-60) <= pz.degree() else zero
    ))


def unpack_x(values, ring):
    X = ring(list(values[:25]))
    Y = ring(list(values[25:60]))
    Z = ring(list(values[60:70]) + [ring.base_ring().one()])
    return X, Y, Z


def jacobian_x(X, Y, Z, surface_A, surface_B):
    base = X.base_ring()
    return shifted_matrix(
        base,
        -3*X**2-surface_A*Z**4,
        2*Y,
        -4*surface_A*X*Z**3-6*surface_B*Z**5,
    )


def hessian_direction(X, Z, direction, ring, surface_A, surface_B):
    VX = ring(list(direction[:25]))
    VY = ring(list(direction[25:60]))
    VZ = ring(list(direction[60:70]))
    return shifted_matrix(
        X.base_ring(),
        -6*X*VX-4*surface_A*Z**3*VZ,
        2*VY,
        -4*surface_A*(VX*Z**3+3*X*Z**2*VZ)-30*surface_B*Z**4*VZ,
    )


def third_mixed_direction(X, Z, left, right, ring, surface_A, surface_B):
    LX = ring(list(left[:25]))
    LZ = ring(list(left[60:70]))
    RX = ring(list(right[:25]))
    RZ = ring(list(right[60:70]))
    return shifted_matrix(
        X.base_ring(),
        -6*RX*LX-12*surface_A*Z**2*RZ*LZ,
        ring.zero(),
        -4*surface_A*(
            3*LX*Z**2*RZ+3*RX*Z**2*LZ+6*X*Z*RZ*LZ
        )-120*surface_B*Z**3*RZ*LZ,
    )


def first_deflated_jacobian(x_values, lambda_values, ring, surface_A, surface_B,
                             lambda_anchor):
    X, Y, Z = unpack_x(x_values, ring)
    J = jacobian_x(X, Y, Z, surface_A, surface_B)
    HL = hessian_direction(X, Z, lambda_values, ring, surface_A, surface_B)
    base = X.base_ring()
    normalization = matrix(base, 1, 70, lambda unused, column: (
        base.one() if column == lambda_anchor else base.zero()
    ))
    return block_matrix(base, [
        [J, zero_matrix(base, 73, 70)],
        [HL, J],
        [zero_matrix(base, 1, 70), normalization],
    ])


def second_deflated_jacobian(values, ring, surface_A, surface_B,
                              lambda_anchor, mu_anchor):
    x_values = values[:70]
    lambda_values = values[70:140]
    mu_values = values[140:280]
    mu_x, mu_lambda = mu_values[:70], mu_values[70:140]
    X, unused_Y, Z = unpack_x(x_values, ring)
    A1 = first_deflated_jacobian(
        x_values, lambda_values, ring, surface_A, surface_B, lambda_anchor,
    )
    HV = hessian_direction(X, Z, mu_x, ring, surface_A, surface_B)
    HE = hessian_direction(X, Z, mu_lambda, ring, surface_A, surface_B)
    mixed = third_mixed_direction(
        X, Z, lambda_values, mu_x, ring, surface_A, surface_B,
    )
    base = X.base_ring()
    DA1 = block_matrix(base, [
        [HV, zero_matrix(base, 73, 70)],
        [mixed+HE, HV],
        [zero_matrix(base, 1, 70), zero_matrix(base, 1, 70)],
    ])
    normalization = matrix(base, 1, 140, lambda unused, column: (
        base.one() if column == mu_anchor else base.zero()
    ))
    return block_matrix(base, [
        [A1, zero_matrix(base, 147, 140)],
        [DA1, A1],
        [zero_matrix(base, 1, 140), normalization],
    ])


def residual(values, ring, surface_A, surface_B, lambda_anchor, mu_anchor):
    base = ring.base_ring()
    x_values = values[:70]
    lambda_values = vector(base, values[70:140])
    mu_values = vector(base, values[140:280])
    X, Y, Z = unpack_x(x_values, ring)
    equation = Y**2-X**3-surface_A*X*Z**4-surface_B*Z**6
    F_values = vector(base, [equation[index] for index in range(73)])
    J = jacobian_x(X, Y, Z, surface_A, surface_B)
    A1 = first_deflated_jacobian(
        x_values, lambda_values, ring, surface_A, surface_B, lambda_anchor,
    )
    return vector(base,
        list(F_values)
        + list(J*lambda_values)
        + [lambda_values[lambda_anchor]-base.one()]
        + list(A1*mu_values)
        + [mu_values[mu_anchor]-base.one()]
    )


def seed_for_candidate(candidate_index):
    record = candidates["search"]["candidates"][candidate_index]
    assert record["shape_Xnum_Xden_Ynum_Yden"] == [24, 20, 34, 30]
    X = RF(record["x"]["numerator_coefficients_low_to_high"])
    Y = RF(record["y"]["numerator_coefficients_low_to_high"])
    denominator_x = RF(record["x"]["denominator_coefficients_low_to_high"])
    denominator_y = RF(record["y"]["denominator_coefficients_low_to_high"])
    assert denominator_x.is_square()
    Z = denominator_x.sqrt()
    if Z**3 != denominator_y:
        Z = -Z
    assert Z.degree() == 10 and Z**3 == denominator_y
    assert Y**2 == X**3+A_F*X*Z**4+B_F*Z**6
    X_coefficients = X.list()
    Y_coefficients = Y.list()
    x_values = (
        X_coefficients+[F.zero()]*(25-len(X_coefficients))
        + Y_coefficients+[F.zero()]*(35-len(Y_coefficients))
        + Z.list()[:10]
    )
    X, Y, Z = unpack_x(x_values, RF)
    J = jacobian_x(X, Y, Z, A_F, B_F)
    assert J.rank() == 69
    lambda_values = list(J.right_kernel().basis()[0])
    lambda_anchor = next(index for index, value in enumerate(lambda_values) if value)
    scale = lambda_values[lambda_anchor]
    lambda_values = [value/scale for value in lambda_values]
    A1 = first_deflated_jacobian(
        x_values, lambda_values, RF, A_F, B_F, lambda_anchor,
    )
    assert A1.rank() == 139
    mu_values = list(A1.right_kernel().basis()[0])
    mu_anchor = next(index for index, value in enumerate(mu_values) if value)
    scale = mu_values[mu_anchor]
    mu_values = [value/scale for value in mu_values]
    values = x_values+lambda_values+mu_values
    B2 = second_deflated_jacobian(
        values, RF, A_F, B_F, lambda_anchor, mu_anchor,
    )
    assert B2.rank() == 280
    assert not residual(values, RF, A_F, B_F, lambda_anchor, mu_anchor)
    return record, values, lambda_anchor, mu_anchor, B2


def minimum_valuation(values, fallback):
    valuations = [int(value.valuation()) for value in values if value]
    return min(valuations) if valuations else int(fallback)


def reconstruct_x(values, usable_precision):
    modulus = prime**usable_precision
    answer = []
    for value in values[:70]:
        residue = ZZ.zero() if not value else ZZ(value.lift()) % modulus
        answer.append(QQ(residue.rational_reconstruction(modulus)))
    return answer


def reduce_qq(value):
    value = QQ(value)
    return F(value.numerator())/F(value.denominator())


def lift_candidate(candidate_index):
    record, seed, lambda_anchor, mu_anchor, B2_F = seed_for_candidate(candidate_index)
    pivot_rows = list(map(int, B2_F.transpose().pivots()))
    assert len(pivot_rows) == 280
    minor = matrix(F, [B2_F.row(index) for index in pivot_rows])
    assert minor.det()
    values = vector(KP, [KP(value).add_bigoh(1) for value in seed])
    known_precision = 1
    iterations, attempts = [], []
    exact = None
    while known_precision < maximum_precision:
        working_precision = min(2*known_precision, maximum_precision)
        values = vector(KP, [
            KP(value.lift()).add_bigoh(working_precision) for value in values
        ])
        full_residual = residual(
            values, RP, A_P, B_P, lambda_anchor, mu_anchor,
        )
        chosen_residual = vector(KP, [full_residual[index] for index in pivot_rows])
        B2 = second_deflated_jacobian(
            values, RP, A_P, B_P, lambda_anchor, mu_anchor,
        )
        square = matrix(KP, [B2.row(index) for index in pivot_rows])
        correction = square.solve_right(-chosen_residual)
        values += correction
        after = residual(values, RP, A_P, B_P, lambda_anchor, mu_anchor)
        iterations.append({
            "working_precision_p_adic_digits": working_precision,
            "minimum_residual_valuation": minimum_valuation(after, working_precision),
            "minimum_correction_valuation": minimum_valuation(correction, working_precision),
        })
        print(
            "Q4O323Q207DEFLIFT|candidate={}|precision={}|residual_val={}|correction_val={}".format(
                candidate_index, working_precision,
                iterations[-1]["minimum_residual_valuation"],
                iterations[-1]["minimum_correction_valuation"],
            ),
            flush=True,
        )
        known_precision = working_precision
        if working_precision < reconstruction_start:
            continue
        attempt = {
            "working_precision_p_adic_digits": working_precision,
            "rational_reconstruction_succeeded": False,
            "literal_weierstrass_identity": False,
        }
        try:
            candidate_values = reconstruct_x(values, working_precision-8)
            X, Y, Z = unpack_x(candidate_values, RQ)
            identity = Y**2 == X**3+A_QQ*X*Z**4+B_QQ*Z**6
            reduction = [reduce_qq(value) for value in candidate_values] == seed[:70]
            attempt.update({
                "rational_reconstruction_succeeded": True,
                "literal_weierstrass_identity": bool(identity),
                "exact_reduction_to_seed": bool(reduction),
                "maximum_candidate_rational_bits": max(map(coefficient_bits, candidate_values)),
            })
            if identity and reduction:
                exact = X, Y, Z, candidate_values
        except (ArithmeticError, ValueError, ZeroDivisionError) as error:
            attempt["failure"] = type(error).__name__
        attempts.append(attempt)
        if exact is not None:
            break
    if exact is None:
        return None, {
            "candidate_index": candidate_index,
            "status": "NO_RATIONAL_RECONSTRUCTION",
            "lambda_anchor": lambda_anchor,
            "mu_anchor": mu_anchor,
            "pivot_rows": pivot_rows,
            "minor_mod61": int(minor.det()),
            "iterations": iterations,
            "reconstruction_attempts": attempts,
        }
    X, Y, Z, exact_values = exact
    return (X, Y, Z), {
        "candidate_index": candidate_index,
        "status": "PASS_EXACT_QQ_SECTION",
        "lambda_anchor": lambda_anchor,
        "mu_anchor": mu_anchor,
        "pivot_rows": pivot_rows,
        "minor_mod61": int(minor.det()),
        "successful_precision_p_adic_digits": known_precision,
        "maximum_rational_bits": max(map(coefficient_bits, exact_values)),
        "iterations": iterations,
        "reconstruction_attempts": attempts,
    }


trials = []
exact = None
for candidate_index in candidate_indices:
    candidate_started = time.monotonic()
    section, trial = lift_candidate(candidate_index)
    trial["runtime_seconds"] = time.monotonic()-candidate_started
    trials.append(trial)
    if section is not None:
        exact = candidate_index, section, trial
        break

payload = {
    "schema": "elkies-k3.h92-q4o323-q207-deflated-horizontal-qq.v1",
    "status": (
        "PASS_EXACT_QQ_Q4O323_DEFLATED_CANDIDATE_SECTION"
        if exact is not None else "REJECTED_Q4O323_Q207_NO_RATIONAL_DEFLATED_LIFT"
    ),
    "prime": int(prime),
    "trials": trials,
    "method": {
        "coefficient_variables": 70,
        "first_deflated_variables": 140,
        "second_deflated_variables": 280,
        "second_deflated_equations": 295,
        "rank_sequence_mod61": [69, 139, 280],
        "large_Groebner_required": False,
        "elimination_required": False,
        "runtime_seconds": time.monotonic()-started,
    },
    "proof_boundary": (
        "An exact QQ section, if emitted, still requires the independent marked q207 "
        "class gate. Candidate 5887 is independently rejected as q207 by its mod-61 "
        "Abel trace and by its 6A1, rather than 5A1, smooth-chord child."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (POINTING, CANDIDATES)],
        "sha256": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in (POINTING, CANDIDATES)
        },
    },
}
if exact is not None:
    selected_index, (X, Y, Z), trial = exact
    payload["selected_candidate_index"] = selected_index
    payload["section"] = {
        "X_coefficients_low_to_high": [str(value) for value in X.list()],
        "Y_coefficients_low_to_high": [str(value) for value in Y.list()],
        "Z_coefficients_low_to_high": [str(value) for value in Z.list()],
        "degrees_X_Y_Z": [int(X.degree()), int(Y.degree()), int(Z.degree())],
        "exact_projective_weierstrass_identity": True,
        "P_dot_O": 10,
    }
output = args.output if args.output.is_absolute() else ROOT / args.output
output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q4O323Q207DEFLIFT|status={}|selected={}|runtime={:.3f}|output={}".format(
        payload["status"], payload.get("selected_candidate_index"),
        payload["method"]["runtime_seconds"], output,
    ),
    flush=True,
)
