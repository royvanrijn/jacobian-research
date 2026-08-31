#!/usr/bin/env sage-python
"""Lift regular q323 polynomial P.O=0 shell seeds to exact QQ sections.

Each seed has twelve coefficients (deg X<=4, deg Y<=6) and a rank-twelve
Jacobian for the thirteen Weierstrass coefficient equations.  Ordinary
Newton--Hensel lifting and rational reconstruction are therefore enough; no
elimination or Groebner basis is used.  Failure to reconstruct is recorded
per seed because most finite-field sections need not descend to QQ.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, Qp, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q4o323-component2-pointing-qq.json"
SHELL = LOCAL / "q4o323-p0-shell-mod61.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--indices", type=int, nargs="*")
parser.add_argument("--limit", type=int, default=10)
parser.add_argument("--reconstruction-start", type=int, default=64)
parser.add_argument("--maximum-precision", type=int, default=1024)
parser.add_argument(
    "--output", type=Path,
    default=LOCAL / "q4o323-regular-p0-shell-qq.json",
)
args = parser.parse_args()
assert 16 <= args.reconstruction_start <= args.maximum_precision
started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficient_bits(value):
    value = QQ(value)
    return max(abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())


model = json.loads(MODEL.read_text())
shell = json.loads(SHELL.read_text())
assert model["status"] == "PASS_EXACT_QQ_Q4O323_OLD_A11_COMPONENT2_POINTING"
assert shell["status"] == "PASS_MODP_Q4O323_COMPLETE_POLYNOMIAL_P0_SHELL"
prime = ZZ(shell["prime"])
assert prime == 61

RQ = PolynomialRing(QQ, "u")
A_QQ = RQ(model["global_short_model"]["A_coefficients_low_to_high"])
B_QQ = RQ(model["global_short_model"]["B_coefficients_low_to_high"])
F = GF(prime)
RF = PolynomialRing(F, "u")


def reduce_qq(value, field=F):
    value = QQ(value)
    return field(value.numerator())/field(value.denominator())


A_F = RF([reduce_qq(value) for value in A_QQ])
B_F = RF([reduce_qq(value) for value in B_QQ])
K = Qp(prime, prec=args.maximum_precision, type="capped-rel")
RP = PolynomialRing(K, "u")
A_P = RP([K(value) for value in A_QQ])
B_P = RP([K(value) for value in B_QQ])


def split(values, ring):
    return ring(list(values[:5])), ring(list(values[5:12]))


def residual(values):
    X, Y = split(values, RP)
    equation = Y**2-X**3-A_P*X-B_P
    return vector(K, [equation[index] for index in range(13)])


def jacobian(values, ring, surface_A):
    X, Y = split(values, ring)
    dx = -3*X**2-surface_A
    dy = 2*Y
    zero = ring.base_ring().zero()
    return matrix(ring.base_ring(), 13, 12, lambda row, column: (
        dx[row-column]
        if column < 5 and 0 <= row-column <= dx.degree() else
        dy[row-(column-5)]
        if column >= 5 and 0 <= row-(column-5) <= dy.degree() else zero
    ))


def minimum_valuation(values, fallback):
    nonzero = [int(value.valuation()) for value in values if value]
    return min(nonzero) if nonzero else int(fallback)


def reconstruct(values, precision):
    modulus = prime**precision
    answer = []
    for value in values:
        residue = ZZ.zero() if not value else ZZ(value.lift()) % modulus
        answer.append(QQ(residue.rational_reconstruction(modulus)))
    return answer


def padded_seed(record):
    x = list(map(F, record["x_coefficients_low_to_high"]))
    y = list(map(F, record["y_coefficients_low_to_high"]))
    return vector(F, x+[F.zero()]*(5-len(x))+y+[F.zero()]*(7-len(y)))


def lift_seed(shell_index, record):
    seed = padded_seed(record)
    X_F, Y_F = split(seed, RF)
    assert Y_F**2 == X_F**3+A_F*X_F+B_F
    J_F = jacobian(seed, RF, A_F)
    assert J_F.rank() == 12
    pivot_rows = list(map(int, J_F.transpose().pivots()))
    minor = matrix(F, [J_F.row(index) for index in pivot_rows]).det()
    assert minor

    values = vector(K, [K(value).add_bigoh(1) for value in seed])
    known = 1
    iterations = []
    attempts = []
    exact = None
    while known < args.maximum_precision:
        working = min(2*known, args.maximum_precision)
        values = vector(K, [K(value.lift()).add_bigoh(working) for value in values])
        full = residual(values)
        square = matrix(K, [
            jacobian(values, RP, A_P).row(index) for index in pivot_rows
        ])
        correction = square.solve_right(-vector(K, [full[index] for index in pivot_rows]))
        values += correction
        after = residual(values)
        iterations.append({
            "precision": working,
            "minimum_residual_valuation": minimum_valuation(after, working),
            "minimum_correction_valuation": minimum_valuation(correction, working),
        })
        known = working
        if working < args.reconstruction_start:
            continue
        attempt = {"precision": working, "rational_reconstruction_succeeded": False}
        try:
            candidate = reconstruct(values, working-8)
            X_QQ, Y_QQ = split(candidate, RQ)
            identity = Y_QQ**2 == X_QQ**3+A_QQ*X_QQ+B_QQ
            reduction = [reduce_qq(value) for value in candidate] == list(seed)
            attempt.update({
                "rational_reconstruction_succeeded": True,
                "literal_weierstrass_identity": bool(identity),
                "exact_seed_reduction": bool(reduction),
                "maximum_rational_bits": max(map(coefficient_bits, candidate)),
            })
            if identity and reduction:
                exact = X_QQ, Y_QQ, candidate
        except (ArithmeticError, ValueError, ZeroDivisionError) as error:
            attempt["failure"] = type(error).__name__
        attempts.append(attempt)
        if exact is not None:
            break

    base = {
        "shell_index": shell_index,
        "mod61_minor": int(minor),
        "iterations": iterations,
        "reconstruction_attempts": attempts,
    }
    if exact is None:
        base["status"] = "NO_QQ_RECONSTRUCTION"
        return base
    X_QQ, Y_QQ, candidate = exact
    base.update({
        "status": "PASS_EXACT_QQ_SECTION",
        "successful_precision": known,
        "section": {
            "x_coefficients_low_to_high": [str(value) for value in X_QQ.list()],
            "y_coefficients_low_to_high": [str(value) for value in Y_QQ.list()],
            "degrees_x_y": [int(X_QQ.degree()), int(Y_QQ.degree())],
            "P_dot_O": 0,
            "maximum_rational_bits": max(map(coefficient_bits, candidate)),
            "exact_weierstrass_identity": True,
            "exact_seed_reduction": True,
        },
    })
    return base


records = shell["shell"]["records"]
regular = [
    index for index, record in enumerate(records)
    if int(record["ordinary_coefficient_jacobian_rank"]) == 12
]
selected = args.indices if args.indices else regular[:args.limit]
if any(index not in set(regular) for index in selected):
    raise ValueError("every selected shell index must have Jacobian rank 12")

results = []
for shell_index in selected:
    section_started = time.monotonic()
    result = lift_seed(shell_index, records[shell_index])
    result["runtime_seconds"] = time.monotonic()-section_started
    results.append(result)
    print(
        "Q4O323P0LIFT|shell={}|status={}|precision={}|bits={}|runtime={:.3f}".format(
            shell_index, result["status"], result.get("successful_precision"),
            result.get("section", {}).get("maximum_rational_bits"),
            result["runtime_seconds"],
        ), flush=True,
    )

exact = [record for record in results if record["status"] == "PASS_EXACT_QQ_SECTION"]
payload = {
    "schema": "elkies-k3.h92-q4o323-regular-p0-shell-qq.v1",
    "status": "EXPERIMENTAL_Q4O323_REGULAR_P0_QQ_LIFTS",
    "prime": int(prime),
    "regular_shell_count": len(regular),
    "selected_shell_indices": selected,
    "exact_QQ_section_count": len(exact),
    "results": results,
    "method": {
        "coefficient_equations": 13,
        "variables": 12,
        "large_Groebner_required": False,
        "elimination_required": False,
        "runtime_seconds": time.monotonic()-started,
    },
    "proof_boundary": (
        "Each passing record is an exact QQ polynomial section. Its full q323 NS/MW "
        "name and whether the resulting subgroup contains q207 remain separate gates."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (MODEL, SHELL)],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in (MODEL, SHELL)},
    },
}
output = args.output if args.output.is_absolute() else ROOT/args.output
output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q4O323P0LIFT|selected={}|exact={}|runtime={:.3f}|output={}".format(
        len(results), len(exact), payload["method"]["runtime_seconds"], output,
    ), flush=True,
)
