#!/usr/bin/env sage -python
"""Hensel-lift the pole-reduced H92 P2 section from one modular half.

The canonical three-neighbour divisor class has modular double with coordinate
degrees ``(184,180)``.  Its modular halving gives a section with the exact
pole profile forced by the H3 lattice:

    x = X/Z^2,  deg(X,Z)=(46,21),
    y = Y/Z^3,  deg(Y,Z)=(69,21).

For the H92 model in the base coordinate ``t``, the cleared Weierstrass
identity is a square polynomial system

    Y^2 - X^3 - a*X*Z^4 - b*Z^6 = 0.

This script checks its nonsingular modular Jacobian and uses p-adic Newton
iteration to lift the modular solution.  It deliberately reconstructs only
the low-pole half, rather than CRT-lifting the high-degree doubled section.
"""

from sage.all import GF, QQ, ZZ, PolynomialRing, Zp, matrix, vector

import argparse
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "artifacts/generated-results/h92-p2-half-mod-100003-v2.json"
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
parser.add_argument("--precision", type=int, default=32, help="p-adic digits")
parser.add_argument("--output", type=Path)
parser.add_argument("--seed", type=Path, help="previous p-adic residue output")
parser.add_argument("--rank-only", action="store_true")
args = parser.parse_args()

record = json.loads(args.input.read_text())
prime = ZZ(record["prime"])
finite = GF(prime)
finite_ring = PolynomialRing(finite, "t")
tf = finite_ring.gen()

xd = finite_ring(record["x"]["denominator"])
yd = finite_ring(record["y"]["denominator"])
assert xd.degree() == 42 and yd.degree() == 63 and xd.is_square()
zf = xd.sqrt()
if zf.leading_coefficient() != 1:
    zf = -zf
assert zf.leading_coefficient() == 1 and yd == zf**3
xf = finite_ring(record["x"]["numerator"])
yf = finite_ring(record["y"]["numerator"])
assert xf.degree() == 46 and yf.degree() == 69

anchor = SourceFileLoader("h92_p2_hensel_anchor", str(ANCHOR)).load_module()
r, s = anchor.EXPECTED_H92
_, formulas = anchor.parse_h92(H92)
A1, A, B1, B, B2 = tuple(QQ(value(r, s)) for value in formulas)


def polynomial_coefficients(value, degree):
    return [value[index] if index <= value.degree() else value.base_ring()(0) for index in range(degree + 1)]


def identity_and_jacobian(ring, z_values, x_values, y_values):
    """Return the identity coefficient vector and its Jacobian.

    ``z_values`` omits the fixed monic top coefficient of Z.
    """
    t = ring.gen()
    z = sum(ring(z_values[index]) * t**index for index in range(21)) + t**21
    x = sum(ring(x_values[index]) * t**index for index in range(47))
    y = sum(ring(y_values[index]) * t**index for index in range(70))
    a = ring(A1) * t**3 + ring(A) * t**4
    b = ring(B1) * t**5 + ring(B) * t**6 + ring(B2) * t**7
    identity = y**2 - x**3 - a*x*z**4 - b*z**6
    derivatives = []
    for index in range(21):
        derivatives.append(-4*a*x*z**3*t**index - 6*b*z**5*t**index)
    for index in range(47):
        derivatives.append((-3*x**2 - a*z**4)*t**index)
    for index in range(70):
        derivatives.append(2*y*t**index)
    return (
        vector(ring.base_ring(), polynomial_coefficients(identity, 138)),
        matrix(ring.base_ring(), [
            [polynomial_coefficients(derivative, 138)[row] for derivative in derivatives]
            for row in range(139)
        ]),
    )


z0 = polynomial_coefficients(zf, 20)
x0 = polynomial_coefficients(xf, 46)
y0 = polynomial_coefficients(yf, 69)
residue, jacobian = identity_and_jacobian(finite_ring, z0, x0, y0)
assert not residue
assert jacobian.rank() == 138
rows = tuple(jacobian.transpose().pivots())
assert len(rows) == 138
square = jacobian.matrix_from_rows(rows)
assert square.is_invertible()
print(
    f"H92P2HENSEL|prime={prime}|jacobian_rank={jacobian.rank()}|"
    f"rows={len(rows)}|status=PASS_NONSINGULAR",
    flush=True,
)

if args.rank_only:
    raise SystemExit(0)

padic = Zp(prime, prec=args.precision)
padic_ring = PolynomialRing(padic, "t")
# Coercing a finite-field element directly into ``Zp`` preserves only one
# p-adic digit.  Lift through its integer representative so Newton starts
# with full requested precision rather than a collection of O(p) entries.
coefficients = [padic(ZZ(value)) for value in z0 + x0 + y0]
seed_precision = 1
if args.seed:
    seed = json.loads(args.seed.read_text())
    if ZZ(seed["prime"]) != prime or "residues" not in seed:
        raise ValueError("seed has the wrong prime or no p-adic residues")
    seed_precision = int(seed["precision"])
    if seed_precision >= args.precision:
        raise ValueError("seed precision must be smaller than target precision")
    if len(seed["residues"]) != len(coefficients):
        raise ValueError("seed coefficient count is inconsistent")
    coefficients = [padic(ZZ(value)) for value in seed["residues"]]
    print(
        f"H92P2HENSEL|seed_precision={seed_precision}|target_precision={args.precision}",
        flush=True,
    )
target = args.precision


def valuation_floor(values):
    finite_values = [value.valuation() for value in values if value]
    return target if not finite_values else min(finite_values)


for iteration in range(1, 2 * target + 3):
    residual, derivative = identity_and_jacobian(
        padic_ring, coefficients[:21], coefficients[21:68], coefficients[68:]
    )
    valuation = valuation_floor(residual)
    print(
        f"H92P2HENSEL|iteration={iteration}|residual_valuation={valuation}",
        flush=True,
    )
    if valuation >= target:
        break
    correction = derivative.matrix_from_rows(rows).solve_right(
        -vector(padic, [residual[row] for row in rows])
    )
    coefficients = [value + delta for value, delta in zip(coefficients, correction, strict=True)]
else:
    raise ArithmeticError("p-adic Newton iteration did not reach target precision")

modulus = prime**target


def reconstruct(values):
    answer = []
    for value in values:
        try:
            answer.append(str(ZZ(value.lift()).rational_reconstruction(modulus)))
        except (ArithmeticError, ValueError):
            answer.append(None)
    return answer


payload = {
    "schema": "elkies-k3.h92-p2-hensel-lift.v1",
    "prime": int(prime),
    "precision": target,
    "seed_precision": int(seed_precision),
    "jacobian_rank": int(jacobian.rank()),
    "selected_rows": [int(row) for row in rows],
    "complete": False,
    "Z": reconstruct(coefficients[:21]) + ["1"],
    "X": reconstruct(coefficients[21:68]),
    "Y": reconstruct(coefficients[68:]),
    "residues": [str(ZZ(value.lift())) for value in coefficients],
}
payload["complete"] = not any(value is None for part in ("Z", "X", "Y") for value in payload[part])
if args.output:
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    f"H92P2HENSEL|precision={target}|complete={int(payload['complete'])}|"
    f"status=PASS",
    flush=True,
)
