#!/usr/bin/env sage -python
"""Recover modular samples of the rational H21 height-21/2 section.

At each requested prime this script specializes the pinned H21 entrance
plane cubic at small values of its new base.  It forces Sage's non-flex
conversion using the ancillary rational point (avoiding a generic flex
search), maps the canonical point ``(0,0)`` to the pinned short E7+E8 model,
and interpolates its x-coordinate.  The generic compact shape is degree
``(10,12)`` in the entrance base; the denominator is a scalar times
``u^4 Z_4(u)^2``.

This is a modular reconstruction stage.  It writes distinct window artifacts;
an exact characteristic-zero lift must be rationally reconstructed and checked
in the Weierstrass equation before the section or q=6 pencil is certified.
"""

from sage.all import *

import argparse
import hashlib
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def evaluate_function_mod(function, value, field):
    return field(function.numerator()(ZZ(value))) / field(
        function.denominator()(ZZ(value))
    )


def third_tangent_point(cubic, point):
    """Return the residual point on the tangent at a smooth cubic point."""

    ring = cubic.parent()
    field = ring.base_ring()
    variables = ring.gens()
    point = vector(field, point)
    assert cubic(*point) == 0
    gradient = vector(field, [cubic.derivative(variable)(*point) for variable in variables])
    if not gradient:
        raise ArithmeticError("singular cubic point")
    candidates = (
        vector(field, [gradient[1], -gradient[0], 0]),
        vector(field, [gradient[2], 0, -gradient[0]]),
        vector(field, [0, gradient[2], -gradient[1]]),
    )
    direction = next(
        (
            candidate
            for candidate in candidates
            if candidate and matrix(field, [point, candidate]).rank() == 2
        ),
        None,
    )
    if direction is None:
        raise ArithmeticError("failed to choose a tangent direction")
    assert gradient * direction == 0

    univariate = PolynomialRing(field, "lambda")
    parameter = univariate.gen()
    substituted = univariate(
        sum(
            coefficient
            * prod(
                (point[index] + parameter * direction[index]) ** exponent
                for index, exponent in enumerate(monomial)
            )
            for monomial, coefficient in cubic.dict().items()
        )
    )
    assert substituted[0] == substituted[1] == 0
    if substituted[3]:
        residual_parameter = -substituted[2] / substituted[3]
        residual = point + residual_parameter * direction
    elif substituted[2]:
        residual = direction
    else:
        raise ArithmeticError("tangent is a component of the cubic")
    assert cubic(*residual) == 0
    return residual


def forced_nonflex_curve(cubic, point):
    """Return the Weierstrass curve from Sage's non-flex construction.

    The public constructor first searches for all rational flexes.  That is
    harmless after numerical specialization but obscures the chosen origin
    modulo primes where a flex happens to be rational.  This is the exact
    non-flex branch of that constructor, with no generic elimination.
    """

    ring = cubic.parent()
    field = ring.base_ring()
    x, y, z = ring.gens()
    point = vector(field, point)
    point2 = third_tangent_point(cubic, point)
    point3 = third_tangent_point(cubic, point2)

    transform = matrix(field, [point, point2, point3]).transpose()
    transformed = transform.act_on_polynomial(cubic)
    first_substitution = [x * x, y * z, x * z]
    quartic = transformed(first_substitution) // (x**2 * z)
    leading_x = field(quartic.coefficient(x**3))
    leading_y = field(quartic.coefficient(y * y * z))
    if not leading_x or not leading_y:
        raise ArithmeticError("degenerate non-flex normalization")
    normalized = quartic(
        [-x, y / leading_y, leading_x * leading_y * z]
    ) / leading_x
    curve = EllipticCurve(normalized([x, y, 1]))
    assert curve(0, 0)
    return curve


def interpolate_x(samples, field):
    # In the converted standard base T=1/u, the height-21/2 section has
    # x=T^2*N_10/Z_4^2.  Hence its entrance-base degrees are (10,12).
    # Exceptional primes can cancel a common factor (p=101 gives (8,10));
    # those primes are deliberately skipped by the unique-kernel gate.
    numerator_degree = 10
    denominator_degree = 12
    rows = []
    for base, x_value in samples:
        rows.append(
            [base**index for index in range(numerator_degree + 1)]
            + [
                -x_value * base**index
                for index in range(denominator_degree + 1)
            ]
        )
    kernel = matrix(field, rows).right_kernel()
    if kernel.dimension() != 1:
        raise ArithmeticError(f"interpolation kernel has dimension {kernel.dimension()}")
    vector_value = kernel.basis()[0]
    polynomial_ring = PolynomialRing(field, "u")
    numerator = polynomial_ring(list(vector_value[: numerator_degree + 1]))
    denominator = polynomial_ring(list(vector_value[numerator_degree + 1 :]))
    if not numerator[0]:
        raise ArithmeticError("x numerator has zero constant coefficient")
    scale = numerator[0] ** -1
    numerator *= scale
    denominator *= scale
    if not all(
        denominator(base) and numerator(base) / denominator(base) == x_value
        for base, x_value in samples
    ):
        raise ArithmeticError("interpolated x-coordinate failed a sample")
    return numerator, denominator


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime-start", type=int, required=True)
parser.add_argument("--prime-count", type=int, default=12)
parser.add_argument("--target", choices=("h21", "h92"), default="h21")
parser.add_argument("--output", type=Path, required=True)
arguments = parser.parse_args()
if arguments.prime_count < 1 or arguments.prime_count > 36:
    raise SystemExit("--prime-count must lie in [1,36]")

# The anchor is pure Python despite its .sage suffix.
anchor = SourceFileLoader("h21_anchor", str(ANCHOR)).load_module()
entrance_cubic, _ = anchor.h21_entrance_cubic(*anchor.EXPECTED_H21)
point_x = anchor.EXPECTED_CUBIC_POINT_X
h21_coefficients = anchor.h21_coefficients(*anchor.EXPECTED_H21)
if arguments.target == "h21":
    short_coefficients = h21_coefficients
    base_multiplier = QQ(1)
else:
    h92_ring, h92_formulas = anchor.parse_h92(H92)
    r92, s92 = anchor.EXPECTED_H92
    short_coefficients = tuple(QQ(value(r92, s92)) for value in h92_formulas)
    base_scale, unused_twist = anchor.model_isomorphism(
        h21_coefficients, short_coefficients
    )
    # T_H21=mu*T_H92, hence u_H92=mu*u_H21.
    base_multiplier = base_scale

records = []
candidate = ZZ(arguments.prime_start)
while len(records) < arguments.prime_count:
    prime = candidate.next_prime()
    candidate = prime
    if prime in (2, 3):
        continue
    field = GF(prime)
    try:
        reduced_point_x = field(point_x)
        A1, A, B1, B, B2 = map(field, short_coefficients)
    except (ZeroDivisionError, TypeError, ValueError):
        continue
    if not A1 or not B1:
        continue

    projective_ring = PolynomialRing(field, names=("X", "T", "Z"))
    X, T, Z = projective_ring.gens()
    samples = []
    for integer_base in range(1, int(prime)):
        base = field(integer_base)
        try:
            cubic = projective_ring(
                sum(
                    evaluate_function_mod(coefficient, integer_base, field)
                    * X**x_degree * T**t_degree
                    * Z ** (3 - x_degree - t_degree)
                    for (x_degree, t_degree), coefficient
                    in entrance_cubic.dict().items()
                )
            )
            intermediate = forced_nonflex_curve(
                cubic, [reduced_point_x, 0, 1]
            )
            interpolation_base = field(base_multiplier) * base
            old_base = interpolation_base**-1
            target = EllipticCurve(
                field,
                [
                    0,
                    0,
                    0,
                    A1 * old_base**3 + A * old_base**4,
                    B1 * old_base**5 + B * old_base**6 + B2 * old_base**7,
                ],
            )
            image = intermediate.isomorphism_to(target)(intermediate(0, 0))
            samples.append((interpolation_base, image[0]))
        except (ArithmeticError, ZeroDivisionError, ValueError):
            continue
        if len(samples) >= 28:
            break
    if len(samples) < 23:
        continue

    try:
        numerator, denominator = interpolate_x(samples, field)
    except ArithmeticError as error:
        print(
            f"H21P1MOD|prime={prime}|status=SKIP_INTERPOLATION|error={error}",
            flush=True,
        )
        continue
    function_field = PolynomialRing(field, "u").fraction_field()
    u = function_field.gen()
    x_coordinate = function_field(numerator) / function_field(denominator)
    old_base = 1 / u
    square = (
        x_coordinate**3
        + (A1 * old_base**3 + A * old_base**4) * x_coordinate
        + B1 * old_base**5 + B * old_base**6 + B2 * old_base**7
    )
    if not square.is_square():
        raise ArithmeticError("interpolated x does not give a section modulo p")
    denominator_factorization = denominator.factor()
    records.append(
        {
            "prime": int(prime),
            "sample_count": len(samples),
            "x_numerator": [int(value) for value in numerator.list()],
            "x_denominator": [int(value) for value in denominator.list()],
            "denominator_factorization": str(denominator_factorization),
            "rhs_square": True,
        }
    )
    print(
        f"H21P1MOD|prime={prime}|samples={len(samples)}|"
        f"x_degree={numerator.degree()},{denominator.degree()}|rhs_square=1",
        flush=True,
    )

payload = {
    "schema": "elkies-k3.h21-p1-entrance-modular.v1",
    "status": "PASS_MODULAR_SECTION_WINDOWS",
    "inputs": {
        "anchor": str(ANCHOR.relative_to(ROOT)),
        "anchor_sha256": digest(ANCHOR),
    },
    "prime_start": arguments.prime_start,
    "prime_count": len(records),
    "target_model": arguments.target,
    "interpolation_base": (
        "u_H21" if arguments.target == "h21" else "u_H92=mu*u_H21"
    ),
    "records": records,
    "proof_boundary": (
        "These are exact finite-field interpolations of the ancillary "
        "non-flex construction.  They are reconstruction data, not yet a "
        "characteristic-zero section certificate."
    ),
}
arguments.output.parent.mkdir(parents=True, exist_ok=True)
arguments.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    f"H21P1MOD|primes={len(records)}|output={arguments.output}|"
    "status=PASS_MODULAR_SECTION_WINDOWS",
    flush=True,
)
