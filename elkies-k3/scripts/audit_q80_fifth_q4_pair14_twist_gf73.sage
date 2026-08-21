#!/usr/bin/env sage
"""Audit the quadratic-twist unit in the pair14 fifth-q4 reconstruction.

The marked projection reconstructs a double cover over GF(73)(u).  Passing
from its factorization to the product of odd factors is valid only up to the
factorization unit.  Dropping a nonsquare unit changes the Jacobian by a
quadratic twist and can change the elliptic surface (in particular its Euler
number and reducible fibers).  This checker reconstructs the cover from the
pinned artifact, compares the monic odd-factor model with the unit-preserving
model, and prints exact minimal local Kodaira data.
"""

import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, FunctionField, GF, PolynomialRing, ZZ


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = (
    ROOT / "artifacts/generated-results/"
    "q80-fifth-q4-marked-projection-pair14-gf73.json"
)
artifact_bytes = ARTIFACT.read_bytes()
assert hashlib.sha256(artifact_bytes).hexdigest() == (
    "e46c9925c6870a6f9185f36994a5aef682382bba7a9bf8d2adc3d897420988fa"
)
data = json.loads(artifact_bytes)

finite = GF(73, impl="modn")
base = FunctionField(finite, "u")
u = base.gen()
cover_ring = PolynomialRing(base, "tau")
tau = cover_ring.gen()
cover = cover_ring(sum(
    base(coefficient)*tau**tau_degree*u**u_degree
    for tau_degree, u_degree, coefficient
    in data["integral_double_cover_terms_T_U_coefficient"]
))
factorization = cover.factor()
unit = base(factorization.unit())
odd_part = cover_ring(1)
factor_degrees_exponents = []
for factor, exponent in factorization:
    factor_degrees_exponents.append((int(factor.degree()), int(exponent)))
    if int(exponent) % 2:
        odd_part *= factor
odd_part = odd_part.monic()
corrected_quartic = cover_ring(unit*odd_part)


def jacobian(quartic):
    coefficients = list(quartic.list()) + [base(0)]*5
    e, d, c, b, a = coefficients[:5]
    invariant_i = 12*a*e - 3*b*d + c**2
    invariant_j = (
        72*a*c*e + 9*b*c*d - 27*a*d**2 - 27*b**2*e - 2*c**3
    )
    return EllipticCurve(base, [0, 0, 0, -27*invariant_i, -27*invariant_j])


def valuation_at(value, factor):
    return (
        int(value.numerator().valuation(factor))
        - int(value.denominator().valuation(factor))
    )


def infinity_valuation(value):
    return int(value.denominator().degree() - value.numerator().degree())


def minimal_orders(ord_a, ord_b, ord_delta):
    # x = pi^(2k) x', y = pi^(3k) y' sends the displayed orders to
    # (ord_a-4k, ord_b-6k, ord_delta-12k).  Choose the largest k preserving
    # integrality.  Python/Sage floor division has the required behavior for
    # negative valuations.
    k = min(ord_a//4, ord_b//6)
    return ord_a-4*k, ord_b-6*k, ord_delta-12*k


def kodaira(ord_a, ord_b, ord_delta):
    if ord_delta == 0:
        return "I0", 0, 0
    if ord_a == 0 or ord_b == 0:
        return f"I{ord_delta}", max(0, ord_delta-1), ord_delta
    if ord_delta == 2:
        return "II", 0, 2
    if ord_delta == 3:
        return "III", 1, 3
    if ord_delta == 4:
        return "IV", 2, 4
    if ord_delta == 6 and ord_a >= 2 and ord_b >= 3:
        return "I0*", 4, 6
    if ord_delta >= 7 and ord_a == 2 and ord_b == 3:
        n = ord_delta-6
        return f"I{n}*", ord_delta-2, ord_delta
    if ord_delta == 8:
        return "IV*", 6, 8
    if ord_delta == 9:
        return "III*", 7, 9
    if ord_delta == 10:
        return "II*", 8, 10
    raise ArithmeticError((ord_a, ord_b, ord_delta))


def signature(curve):
    delta = curve.discriminant()
    factors = set()
    for value in (curve.a4(), curve.a6(), delta):
        factors.update(factor for factor, _exponent in value.numerator().factor())
        factors.update(factor for factor, _exponent in value.denominator().factor())
    rows = []
    for factor in sorted(factors, key=lambda item: (item.degree(), str(item))):
        orders = minimal_orders(
            valuation_at(curve.a4(), factor),
            valuation_at(curve.a6(), factor),
            valuation_at(delta, factor),
        )
        if orders[2] > 0:
            fiber, rank, euler = kodaira(*orders)
            rows.append((str(factor), int(factor.degree()), orders, fiber, rank, euler))
    infinity_orders = minimal_orders(
        infinity_valuation(curve.a4()),
        infinity_valuation(curve.a6()),
        infinity_valuation(delta),
    )
    if infinity_orders[2] > 0:
        fiber, rank, euler = kodaira(*infinity_orders)
        rows.append(("infinity", 1, infinity_orders, fiber, rank, euler))
    root_rank = sum(degree*rank for _, degree, _, _, rank, _ in rows)
    euler_number = sum(degree*euler for _, degree, _, _, _, euler in rows)
    return tuple(rows), int(root_rank), int(euler_number)


monic_curve = jacobian(odd_part)
corrected_curve = jacobian(corrected_quartic)
monic_rows, monic_rank, monic_euler = signature(monic_curve)
corrected_rows, corrected_rank, corrected_euler = signature(corrected_curve)

assert tuple(factor_degrees_exponents) == tuple(
    tuple(row) for row in data["cover_factor_degrees_exponents"]
)
assert monic_euler == 48
assert corrected_euler % 12 == 0

print(
    "Q80FIFTHPAIR14TWIST|"
    f"factorization_unit={unit}|unit_num_factorization={unit.numerator().factor()}|"
    f"unit_den_factorization={unit.denominator().factor()}|"
    f"unit_square={int(unit.is_square())}|status=PASS_UNIT_AUDIT",
    flush=True,
)
for label, rows, root_rank, euler_number in (
    ("monic", monic_rows, monic_rank, monic_euler),
    ("unit_preserving", corrected_rows, corrected_rank, corrected_euler),
):
    for place, degree, orders, fiber, rank, euler in rows:
        print(
            "Q80FIFTHPAIR14FIBER|"
            f"model={label}|place={place}|degree={degree}|orders={orders}|"
            f"fiber={fiber}|root_rank={rank}|euler={euler}",
            flush=True,
        )
    print(
        "Q80FIFTHPAIR14SIGNATURE|"
        f"model={label}|geometric_root_rank={root_rank}|"
        f"euler_number={euler_number}|chi={ZZ(euler_number)//12}|"
        "status=PASS_EXACT_LOCAL_MINIMALIZATION",
        flush=True,
    )
