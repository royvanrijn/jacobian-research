#!/usr/bin/env sage -python
"""Reconstruct the global D12 AJ section from degree-14 trace samples.

The input samples are produced by
``probe_h92_q24_a11_close_p24_quintic_modp.sage``.  Newton interpolation and
extended-Euclid Pade reconstruction recover the unique rational x-coordinate
whose denominator is a square and then the y-coordinate with the matching
cube denominator.  This is univariate finite-field arithmetic only.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, is_prime


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--samples",
    type=Path,
    default=LOCAL / "q24-close-p24-aj14-plus-700samples-mod100003.json",
)
parser.add_argument(
    "--output",
    type=Path,
    default=LOCAL / "q24-close-p24-aj14-plus-section-mod100003.json",
)
args = parser.parse_args()

SAMPLES = args.samples.resolve()
Q24 = LOCAL / "q24-d13-to-d12-component-valuation-qq.json"
for path in (SAMPLES, Q24):
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

payload_samples = json.loads(SAMPLES.read_text())
q24 = json.loads(Q24.read_text())
p = ZZ(payload_samples["prime"])
if not is_prime(p) or p in (2, 3):
    raise ArithmeticError("sample prime is not good")
F = GF(p)
R = PolynomialRing(F, "V")
V = R.gen()

raw = payload_samples["q24_stage"]["degree14_AJ_traces"]
branches = {row["branch"] for row in raw}
if len(branches) != 1:
    raise ArithmeticError("reconstruction input must contain one trace branch")
branch = next(iter(branches))
samples = [(F(row["tau"]), F(row["AJ_x"]), F(row["AJ_y"])) for row in raw]
if len({int(row[0]) for row in samples}) != len(samples):
    raise ArithmeticError("duplicate sample abscissa")


def interpolation_polynomial(values):
    interpolation = R.zero()
    modulus = R.one()
    for parameter, value in values:
        scale = modulus(parameter)
        if not scale:
            raise ArithmeticError("duplicate interpolation parameter")
        interpolation += ((value - interpolation(parameter)) / scale) * modulus
        modulus *= V - parameter
    interpolation %= modulus
    assert all(interpolation(parameter) == value for parameter, value in values)
    return interpolation, modulus


def pade_sequence(values):
    interpolation, modulus = interpolation_polynomial(values)
    r0, r1 = modulus, interpolation
    t0, t1 = R.zero(), R.one()
    candidates = []
    while r1:
        numerator, denominator = r1, t1
        common = numerator.gcd(denominator)
        numerator //= common
        denominator //= common
        if denominator:
            scale = denominator.leading_coefficient()
            numerator /= scale
            denominator /= scale
            if all(denominator(parameter) and numerator(parameter) == value * denominator(parameter) for parameter, value in values):
                candidates.append((numerator, denominator))
        quotient, r2 = r0.quo_rem(r1)
        t2 = t0 - quotient * t1
        r0, r1 = r1, r2
        t0, t1 = t1, t2
    return candidates


x_candidates = []
for numerator, denominator in pade_sequence([(t, x) for t, x, unused in samples]):
    if not denominator.is_square():
        continue
    Z = denominator.sqrt().monic()
    numerator /= denominator.leading_coefficient()
    denominator = Z**2
    pole = Z.degree()
    if numerator.degree() > 2 * pole + 4:
        continue
    x_candidates.append((numerator, denominator, Z))

if len(x_candidates) != 1:
    raise ArithmeticError(f"square-denominator x candidates: {len(x_candidates)}")
X, unused_x_denominator, Z = x_candidates[0]
pole_order = Z.degree()

y_candidates = []
for numerator, denominator in pade_sequence([(t, y) for t, unused, y in samples]):
    scale = denominator.leading_coefficient()
    numerator /= scale
    denominator /= scale
    if denominator != Z**3:
        continue
    if numerator.degree() > 3 * pole_order + 6:
        continue
    y_candidates.append((numerator, denominator))

if len(y_candidates) != 1:
    raise ArithmeticError(f"matching-cube y candidates: {len(y_candidates)}")
Y, unused_y_denominator = y_candidates[0]


def red(value):
    value = QQ(value)
    if value.denominator() % p == 0:
        raise ZeroDivisionError(f"bad denominator modulo {p}")
    return F(value.numerator()) / F(value.denominator())


A = R([red(value) for value in q24["child"]["minimal_A_coefficients_low_to_high"]])
B = R([red(value) for value in q24["child"]["minimal_B_coefficients_low_to_high"]])
assert Y**2 == X**3 + A * X * Z**4 + B * Z**6

payload = {
    "schema": "elkies-k3.h3-q24-close-p24-aj14-section-modp.v1",
    "status": "PASS_Q24_CLOSE_P24_AJ14_SECTION_RECONSTRUCTION_MODP",
    "prime": int(p),
    "branch": branch,
    "sample_count": len(samples),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (SAMPLES, Q24)],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (SAMPLES, Q24)
        },
    },
    "section": {
        "X_coefficients_low_to_high": [int(value) for value in X.list()],
        "Y_coefficients_low_to_high": [int(value) for value in Y.list()],
        "Z_coefficients_low_to_high": [int(value) for value in Z.list()],
        "degrees_X_Y_Z": [int(X.degree()), int(Y.degree()), int(Z.degree())],
        "P_dot_O": int(pole_order),
        "exact_modp_weierstrass_identity": True,
    },
    "method": "degree-14 L(15O) samples plus Newton/extended-Euclid Pade reconstruction",
    "proof_boundary": (
        "Exact over the declared finite field. The section still requires MW marking, "
        "multi-prime or characteristic-zero lifting, and literal QQ verification."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q24CLOSEAJSECTION|prime={}|branch={}|samples={}|degrees={},{},{}|PO={}|status={}".format(
        p, branch, len(samples), X.degree(), Y.degree(), Z.degree(), pole_order, payload["status"]
    ),
    flush=True,
)
print(f"OUTPUT|{args.output.resolve()}", flush=True)
