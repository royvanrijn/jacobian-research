#!/usr/bin/env sage
"""Padé-reconstruct Q80 marked coefficients over the CM compositum.

The input is the exact normalized formal branch at CM24.  We compose its
``h=P-P_CM`` series with the exact centered global parameter ``t`` and use a
strict training/holdout split.  Every accepted rational function is then an
exact candidate on the recovered coefficient line; a later global section
identity promotes the collection from local recognition to a certificate.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--jet",
    type=Path,
    default=ROOT / "artifacts/generated-results/"
    "q80-slope-8-87-qqcm-order28-marked-jet.json",
)
parser.add_argument(
    "--local-parameter",
    type=Path,
    default=ROOT / "artifacts/generated-results/"
    "q80-cm24-slope-8-87-qq-local-parameter.json",
)
parser.add_argument("--validation-orders", type=int, default=6)
parser.add_argument(
    "--output",
    type=Path,
    default=ROOT / "artifacts/generated-results/"
    "q80-slope-8-87-qqcm-marked-functions.json",
)
args = parser.parse_args()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


jet_payload = json.loads(args.jet.read_text())
parameter_payload = json.loads(args.local_parameter.read_text())
if jet_payload.get("kind") != "normalized_formal_jet" or jet_payload.get("prime") != "QQCM":
    raise ValueError("expected an exact QQCM marked jet")
if parameter_payload.get("schema") != "q80-cm24-formal-branch-parameter-v1":
    raise ValueError("unexpected local parameter schema")
order = ZZ(jet_payload["order"])
validation = ZZ(args.validation_orders)
fitting_order = order-validation
if validation < 2 or fitting_order < 4:
    raise ValueError("need at least two holdout orders and four fitting orders")

quadratic_six = QuadraticField(-6, "sqrt_minus_six")
s6_seed = quadratic_six.gen()
quadratic_three = QuadraticField(-3, "sqrt_minus_three")
s3_seed = quadratic_three.gen()
composite, embed_three, embed_six, _ = quadratic_three.composite_fields(
    quadratic_six, both_maps=True
)[0]
s3 = embed_three(s3_seed)
s6 = embed_six(s6_seed)

local_ring = PolynomialRing(composite, "t")
t = local_ring.gen()
local_field = local_ring.fraction_field()
series_ring = PowerSeriesRing(composite, "t", default_prec=order)
t_series = series_ring.gen()


def rational_function(record):
    numerator = local_ring(record["numerator"])
    denominator = local_ring(record["denominator"])
    if denominator == 0 or numerator.gcd(denominator) != 1:
        raise ArithmeticError("local parameter function is not reduced")
    return local_field(numerator/denominator)


qq_local_ring = PolynomialRing(QQ, "t")
qq_t = qq_local_ring.gen()
qq_series_ring = PowerSeriesRing(QQ, "t", default_prec=order)
qq_t_series = qq_series_ring.gen()
p_record = parameter_payload["functions"]["P"]
qq_p_numerator = qq_local_ring(p_record["numerator"])
qq_p_denominator = qq_local_ring(p_record["denominator"])
h_of_t = (
    qq_series_ring(qq_p_numerator(qq_t_series))
    / qq_series_ring(qq_p_denominator(qq_t_series))
).add_bigoh(order)
if h_of_t.valuation() != 1:
    raise ArithmeticError("centered P is not a uniformizer at CM24")
h_powers = [qq_series_ring.one()]
for _ in range(1, order):
    h_powers.append((h_powers[-1]*h_of_t).add_bigoh(order))

active_names = tuple(jet_payload["active_variables"])
rows = tuple(
    tuple(composite(value) for value in row)
    for row in jet_payload["coefficients"]
)
if len(rows) != order or any(len(row) != len(active_names) for row in rows):
    raise ValueError("marked jet has inconsistent dimensions")


def compose_h_series(column):
    coefficients = []
    for degree in range(order):
        coefficients.append(sum(
            (rows[index][column]*h_powers[index][degree] for index in range(order)),
            composite.zero(),
        ))
    return series_ring(coefficients).add_bigoh(order)


biquadratic_basis = (composite.one(), s3, s6, s3*s6)
biquadratic_matrix = Matrix(
    QQ,
    [list(value) for value in biquadratic_basis],
).transpose()


def reducer(prime):
    finite = GF(prime)
    r3 = min(finite(-3).sqrt(all=True), key=ZZ)
    r6 = min(finite(-6).sqrt(all=True), key=ZZ)
    images = (finite.one(), r3, r6, r3*r6)
    cache = {}

    def reduce_value(value):
        value = composite(value)
        if value in cache:
            return cache[value]
        coordinates = biquadratic_matrix.solve_right(vector(QQ, list(value)))
        result = sum(
            (finite(coordinate)*image for coordinate, image in zip(coordinates, images)),
            finite.zero(),
        )
        cache[value] = result
        return result

    return finite, reduce_value


def modular_pade_profile(value, prime):
    finite, reduce_value = reducer(prime)
    finite_ring = PolynomialRing(finite, "t")
    finite_t = finite_ring.gen()
    coefficients = [reduce_value(value[index]) for index in range(order)]
    interpolation = finite_ring(coefficients[:fitting_order])
    modulus = finite_t**fitting_order
    r0, r1 = modulus, interpolation
    d0, d1 = finite_ring.zero(), finite_ring.one()
    candidates = []
    while r1:
        if d1:
            common = r1.gcd(d1)
            numerator = r1//common
            denominator = d1//common
            scale = denominator.leading_coefficient()
            numerator /= scale
            denominator /= scale
            residual = []
            for degree in range(order):
                coefficient = sum(
                    (
                        denominator[index]*coefficients[degree-index]
                        for index in range(min(denominator.degree(), degree)+1)
                    ),
                    finite.zero(),
                )
                if degree <= numerator.degree():
                    coefficient -= numerator[degree]
                residual.append(coefficient)
            if not any(residual):
                candidates.append((
                    numerator.degree()+denominator.degree(),
                    numerator.degree(), denominator.degree(),
                ))
        quotient, r2 = r0.quo_rem(r1)
        d2 = d0-quotient*d1
        r0, r1 = r1, r2
        d0, d1 = d1, d2
    if not candidates:
        return None
    return min(candidates)


def pade_reconstruct(value):
    if not value:
        return local_ring.zero(), local_ring.one()
    if value.valuation() == 0 and not (value-value[0]).add_bigoh(order):
        return local_ring(value[0]), local_ring.one()
    profiles = tuple(modular_pade_profile(value, prime) for prime in (31, 73))
    if None in profiles or profiles[0][1:] != profiles[1][1:]:
        return None
    numerator_degree, denominator_degree = profiles[0][1:]
    unknown_count = numerator_degree+1+denominator_degree
    if unknown_count > fitting_order:
        return None
    rows = []
    right_hand_side = []
    for degree in range(unknown_count):
        rows.append(
            [composite.one() if degree == index else composite.zero()
             for index in range(numerator_degree+1)]
            + [
                -value[degree-index] if degree >= index else composite.zero()
                for index in range(denominator_degree)
            ]
        )
        right_hand_side.append(
            value[degree-denominator_degree]
            if degree >= denominator_degree else composite.zero()
        )
    solution = Matrix(composite, rows).solve_right(
        vector(composite, right_hand_side)
    )
    solution = list(solution)
    numerator = local_ring(solution[:numerator_degree+1])
    denominator = local_ring(
        solution[numerator_degree+1:] + [composite.one()]
    )
    residual = []
    for degree in range(order):
        coefficient = sum(
            (
                denominator[index]*value[degree-index]
                for index in range(min(denominator.degree(), degree)+1)
            ),
            composite.zero(),
        )
        if degree <= numerator.degree():
            coefficient -= numerator[degree]
        residual.append(coefficient)
    if any(residual):
        return None
    return numerator, denominator


surface_centers = {
    "D": QQ(-1)/2, "P": QQ(9)/4, "Q": QQ(-9)/4, "E": QQ(-27)/32,
}
functions = {}
unresolved = []
for column, name in enumerate(active_names):
    if name in surface_centers:
        centered = rational_function(parameter_payload["functions"][name])
        value = local_field(centered+surface_centers[name])
        numerator = local_ring(value.numerator())
        denominator = local_ring(value.denominator())
    else:
        reconstruction = pade_reconstruct(compose_h_series(column))
        if reconstruction is None:
            unresolved.append(name)
            continue
        numerator, denominator = reconstruction
    functions[name] = {
        "numerator": str(numerator),
        "denominator": str(denominator),
        "value": f"({numerator})/({denominator})",
        "degrees": [int(numerator.degree()), int(denominator.degree())],
    }
    print(
        f"Q80MARKEDPADE|coordinate={name}|degrees={numerator.degree()}/{denominator.degree()}|"
        f"fitting_orders={fitting_order}|holdout_orders={validation}|status=PASS",
        flush=True,
    )

output = {
    "schema": "q80-qqcm-marked-functions-v1",
    "status": "PASS_PARTIAL_EXACT_PADE" if unresolved else "PASS_EXACT_PADE_ALL_ACTIVE",
    "field": "QQ(sqrt(-3),sqrt(-6))",
    "parameter": "t with CM24 at t=0",
    "order": int(order),
    "fitting_orders": int(fitting_order),
    "holdout_orders": int(validation),
    "functions": functions,
    "unresolved": unresolved,
    "inputs": {
        "jet": {"path": str(args.jet.relative_to(ROOT)), "sha256": sha256(args.jet)},
        "local_parameter": {
            "path": str(args.local_parameter.relative_to(ROOT)),
            "sha256": sha256(args.local_parameter),
        },
    },
    "claim_boundary": (
        "Padé recognition with exact held-out Taylor orders; global marked-section "
        "identities are required before these functions are promoted."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
print(
    f"Q80MARKEDPADE|resolved={len(functions)}|unresolved={len(unresolved)}|"
    f"names={','.join(unresolved) if unresolved else '-'}|output={args.output}|"
    f"status={output['status']}",
    flush=True,
)
