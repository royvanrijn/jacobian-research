#!/usr/bin/env sage -python
"""CRT-lift modular H21 section windows and verify over QQ when ready."""

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


def crt_lift(residues, primes):
    value = ZZ(residues[0])
    modulus = ZZ(primes[0])
    for residue, prime in zip(residues[1:], primes[1:]):
        value = CRT(value, ZZ(residue), modulus, ZZ(prime))
        modulus *= ZZ(prime)
    return value, modulus


def normalized_denominator_root(record):
    """Extract D=d4*u^4*Z4^2 with Z4(0)=1 from one modular record."""

    field = GF(record["prime"])
    polynomial_ring = PolynomialRing(field, "u")
    u_mod = polynomial_ring.gen()
    denominator = polynomial_ring(record["x_denominator"])
    if any(denominator[index] for index in range(4)):
        raise ArithmeticError("denominator is not divisible by u^4")
    d4 = denominator[4]
    if not d4:
        raise ArithmeticError("denominator has zero u^4 coefficient")
    quotient = denominator // (d4 * u_mod**4)
    if quotient * d4 * u_mod**4 != denominator or not quotient.is_square():
        raise ArithmeticError("denominator quotient is not a square")
    root = quotient.sqrt()
    if root[0] != 1:
        root = -root
    if root[0] != 1 or root.degree() != 4 or root**2 != quotient:
        raise ArithmeticError("failed to normalize quartic denominator root")
    return int(d4), [int(root[index]) for index in range(5)]


def third_tangent_point(cubic, point):
    """Residual point on the tangent, without a flex/factor search."""

    ring = cubic.parent()
    field = ring.base_ring()
    variables = ring.gens()
    point = vector(field, point)
    gradient = vector(
        field, [cubic.derivative(variable)(*point) for variable in variables]
    )
    candidates = (
        vector(field, [gradient[1], -gradient[0], 0]),
        vector(field, [gradient[2], 0, -gradient[0]]),
        vector(field, [0, gradient[2], -gradient[1]]),
    )
    direction = next(
        candidate
        for candidate in candidates
        if candidate and matrix(field, [point, candidate]).rank() == 2
    )
    parameter_ring = PolynomialRing(field, "lambda")
    parameter = parameter_ring.gen()
    substituted = parameter_ring(
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
        residual = point - substituted[2] / substituted[3] * direction
    elif substituted[2]:
        residual = direction
    else:
        raise ArithmeticError("tangent is a component of the cubic")
    assert cubic(*residual) == 0
    return residual


def forced_nonflex_curve(cubic, point):
    """Exact non-flex conversion with the ancillary point as origin marker."""

    ring = cubic.parent()
    field = ring.base_ring()
    x_projective, y_projective, z_projective = ring.gens()
    point = vector(field, point)
    point2 = third_tangent_point(cubic, point)
    point3 = third_tangent_point(cubic, point2)
    transform = matrix(field, [point, point2, point3]).transpose()
    transformed = transform.act_on_polynomial(cubic)
    quartic = transformed(
        [x_projective**2, y_projective * z_projective, x_projective * z_projective]
    ) // (x_projective**2 * z_projective)
    leading_x = field(quartic.coefficient(x_projective**3))
    leading_y = field(quartic.coefficient(y_projective**2 * z_projective))
    normalized = quartic(
        [-x_projective, y_projective / leading_y, leading_x * leading_y * z_projective]
    ) / leading_x
    curve = EllipticCurve(normalized([x_projective, y_projective, 1]))
    assert curve(0, 0)
    return curve


def exact_oriented_image(anchor, entrance_base, base_multiplier, coefficients):
    """Image of the ancillary marked point on one exact rational fiber."""

    entrance_cubic, _ = anchor.h21_entrance_cubic(*anchor.EXPECTED_H21)
    projective_ring = PolynomialRing(QQ, names=("X", "T", "Z"))
    X, T_projective, Z = projective_ring.gens()
    cubic = projective_ring(
        sum(
            QQ(coefficient.numerator()(entrance_base))
            / QQ(coefficient.denominator()(entrance_base))
            * X**x_degree
            * T_projective**t_degree
            * Z ** (3 - x_degree - t_degree)
            for (x_degree, t_degree), coefficient in entrance_cubic.dict().items()
        )
    )
    intermediate = forced_nonflex_curve(
        cubic, [anchor.EXPECTED_CUBIC_POINT_X, 0, 1]
    )
    target_base = base_multiplier * entrance_base
    old_base = 1 / target_base
    A1_value, A_value, B1_value, B_value, B2_value = coefficients
    target = EllipticCurve(
        QQ,
        [
            0,
            0,
            0,
            A1_value * old_base**3 + A_value * old_base**4,
            B1_value * old_base**5
            + B_value * old_base**6
            + B2_value * old_base**7,
        ],
    )
    image = intermediate.isomorphism_to(target)(intermediate(0, 0))
    return target_base, image


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("inputs", nargs="+", type=Path)
parser.add_argument("--output", required=True, type=Path)
arguments = parser.parse_args()

records = []
input_metadata = []
target_models = set()
for path in arguments.inputs:
    payload = json.loads(path.read_text())
    assert payload["status"] == "PASS_MODULAR_SECTION_WINDOWS"
    records.extend(payload["records"])
    target_models.add(payload.get("target_model", "h21"))
    input_metadata.append({"path": str(path), "sha256": digest(path)})
records.sort(key=lambda record: record["prime"])
assert len(target_models) == 1
target_model = target_models.pop()
primes = [ZZ(record["prime"]) for record in records]
assert len(primes) == len(set(primes))

modulus = prod(primes)

# Exploit the exact pole structure before attempting a lift.  Treating all
# thirteen denominator coefficients independently wastes nearly half of the
# CRT information and hides the square denominator forced by P.O=4.
structured_modular = [normalized_denominator_root(record) for record in records]
structured_residue_columns = [
    [record["x_numerator"][index] for record in records]
    for index in range(1, 11)
] + [
    [entry[0] for entry in structured_modular]
] + [
    [entry[1][index] for entry in structured_modular]
    for index in range(1, 5)
]
structured_residues = []
structured_coefficients = []
structured_reconstruction_complete = True
for column in structured_residue_columns:
    residue, row_modulus = crt_lift(column, primes)
    assert row_modulus == modulus
    structured_residues.append(residue)
    try:
        structured_coefficients.append(residue.rational_reconstruction(modulus))
    except ArithmeticError:
        structured_coefficients.append(None)
        structured_reconstruction_complete = False

coefficient_rows = []
residue_rows = []
reconstruction_complete = True
for key, length in (("x_numerator", 11), ("x_denominator", 13)):
    assert all(len(record[key]) == length for record in records)
    row = []
    residue_row = []
    for index in range(length):
        residue, row_modulus = crt_lift(
            [record[key][index] for record in records], primes
        )
        assert row_modulus == modulus
        residue_row.append(residue)
        try:
            row.append(residue.rational_reconstruction(modulus))
        except ArithmeticError:
            reconstruction_complete = False
            row.append(None)
    coefficient_rows.append(row)
    residue_rows.append(residue_row)

numerator_coefficients, denominator_coefficients = coefficient_rows
ring = PolynomialRing(QQ, "u")
u = ring.gen()
anchor = SourceFileLoader("h21_anchor_lift", str(ANCHOR)).load_module()
h21_coefficients = anchor.h21_coefficients(*anchor.EXPECTED_H21)
if target_model == "h21":
    A1, A, B1, B, B2 = h21_coefficients
    base_multiplier = QQ(1)
else:
    h92_ring, h92_formulas = anchor.parse_h92(H92)
    r92, s92 = anchor.EXPECTED_H92
    A1, A, B1, B, B2 = tuple(
        QQ(value(r92, s92)) for value in h92_formulas
    )
    base_multiplier, unused_twist = anchor.model_isomorphism(
        h21_coefficients, (A1, A, B1, B, B2)
    )


def exact_section(candidate_numerator, candidate_denominator):
    """Return (exact, rhs) for a proposed normalized x-coordinate."""

    numerator = ring(candidate_numerator)
    denominator = ring(candidate_denominator)
    if numerator[0] != 1 or numerator.degree() != 10 or denominator.degree() != 12:
        return False, None
    field = ring.fraction_field()
    u_field = field.gen()
    x = field(numerator) / field(denominator)
    T = 1 / u_field
    rhs = (
        x**3
        + (A1 * T**3 + A * T**4) * x
        + B1 * T**5 + B * T**6 + B2 * T**7
    )
    return rhs.is_square(), rhs


def structured_section(candidate):
    """Expand the 15 independent coefficients and test the exact section."""

    if len(candidate) != 15 or any(value is None for value in candidate):
        return False, None, None, None
    numerator = ring([QQ(1)] + list(candidate[:10]))
    d4 = QQ(candidate[10])
    denominator_root = ring([QQ(1)] + list(candidate[11:]))
    denominator = d4 * u**4 * denominator_root**2
    exact_candidate, candidate_rhs = exact_section(numerator, denominator)
    return exact_candidate, candidate_rhs, numerator, denominator


exact = False
rhs = None
reconstruction_method = "independent_rational_reconstruction"
structured_exact = False
structured_root = None
if structured_reconstruction_complete:
    structured_exact, rhs, numerator, denominator = structured_section(
        structured_coefficients
    )
    if structured_exact:
        exact = True
        reconstruction_complete = True
        reconstruction_method = "structured_independent_rational_reconstruction"
        numerator_coefficients = list(numerator)
        denominator_coefficients = list(denominator)
        structured_root = ring([QQ(1)] + list(structured_coefficients[11:]))
elif reconstruction_complete:
    exact, rhs = exact_section(numerator_coefficients, denominator_coefficients)

# Simultaneous rational reconstruction uses the fact that every coefficient
# is normalized by the same nonzero numerator constant.  The lattice consists
# of (n_1,...,n_23,d) with n_i == d*r_i (mod M).  A short vector recovers the
# common denominator and all numerators at once, usually with about half the
# CRT modulus required by coefficient-wise rational reconstruction.
if not exact:
    residues = structured_residues
    dimension = len(residues) + 1
    basis = matrix(ZZ, dimension, dimension)
    for index in range(dimension - 1):
        basis[index, index] = modulus
    for index, residue in enumerate(residues):
        basis[dimension - 1, index] = residue
    basis[dimension - 1, dimension - 1] = 1
    reduced = basis.LLL(delta=0.99)
    for vector_value in reduced.rows():
        common_denominator = ZZ(vector_value[-1])
        if not common_denominator:
            continue
        simultaneous = [
            QQ(value) / common_denominator for value in vector_value[:-1]
        ]
        candidate_exact, candidate_rhs, candidate_numerator, candidate_denominator = (
            structured_section(simultaneous)
        )
        if candidate_exact:
            numerator = candidate_numerator
            denominator = candidate_denominator
            numerator_coefficients = list(numerator)
            denominator_coefficients = list(denominator)
            structured_coefficients = simultaneous
            structured_root = ring([QQ(1)] + simultaneous[11:])
            reconstruction_complete = True
            structured_reconstruction_complete = True
            reconstruction_method = "structured_simultaneous_lll"
            exact = True
            rhs = candidate_rhs
            break

if exact:
    numerator = ring(numerator_coefficients)
    denominator = ring(denominator_coefficients)
    assert numerator[0] == 1

output = {
    "schema": "elkies-k3.h21-p1-entrance-lift.v1",
    "status": (
        f"PASS_EXACT_{target_model.upper()}_P1"
        if exact
        else "NEED_MORE_CRT_PRIMES"
    ),
    "inputs": input_metadata,
    "prime_count": len(primes),
    "primes": [int(prime) for prime in primes],
    "crt_modulus_bits": int(modulus.nbits()),
    "target_model": target_model,
    "structured_denominator": {
        "identity": "D(u)=d4*u^4*Z4(u)^2 with Z4(0)=1",
        "d4": (
            None
            if len(structured_coefficients) < 11
            or structured_coefficients[10] is None
            else str(structured_coefficients[10])
        ),
        "Z4_coefficients": [
            "1",
            *[
                None if value is None else str(value)
                for value in structured_coefficients[11:]
            ],
        ],
        "rational_reconstruction_complete": structured_reconstruction_complete,
    },
    "x_entrance_base": {
        "numerator_coefficients": [
            None if value is None else str(value)
            for value in numerator_coefficients
        ],
        "denominator_coefficients": [
            None if value is None else str(value)
            for value in denominator_coefficients
        ],
        "degrees": [10, 12],
    },
    "rational_reconstruction_complete": reconstruction_complete,
    "reconstruction_method": reconstruction_method,
    "exact_weierstrass_square": bool(exact),
}
if exact:
    y = rhs.sqrt()
    orientation_base, oriented_image = exact_oriented_image(
        anchor, QQ(2), base_multiplier, (A1, A, B1, B, B2)
    )
    function_field = ring.fraction_field()
    x_function = function_field(numerator) / function_field(denominator)
    if x_function(orientation_base) != oriented_image[0]:
        raise ArithmeticError("exact marked-fiber x incidence failed")
    if y(orientation_base) == -oriented_image[1]:
        y = -y
    if y(orientation_base) != oriented_image[1]:
        raise ArithmeticError("neither square-root sign matches the marked point")
    output["y_entrance_base"] = {
        "numerator_coefficients": [str(value) for value in y.numerator().list()],
        "denominator_coefficients": [str(value) for value in y.denominator().list()],
        "degrees": [int(y.numerator().degree()), int(y.denominator().degree())],
    }
    output["orientation_incidence"] = {
        "source": "ancillary non-flex point on the split A2+A6+E8 3-neighbor cubic",
        "h21_entrance_base": "2",
        "target_entrance_base": str(orientation_base),
        "target_point_x": str(oriented_image[0]),
        "target_point_y": str(oriented_image[1]),
        "selected_square_root_sign": 1,
        "exact_match": True,
    }
    output["point_identity"] = "y^2=x^3+(A1*T^3+A*T^4)x+B1*T^5+B*T^6+B2*T^7, T=1/u"

arguments.output.parent.mkdir(parents=True, exist_ok=True)
arguments.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
print(
    f"H21P1LIFT|primes={len(primes)}|modulus_bits={modulus.nbits()}|"
    f"exact_square={int(exact)}|status={output['status']}",
    flush=True,
)
