#!/usr/bin/env sage-python
"""Decode and certify joint resolved-depth13 NS0024 closed points.

The input is a zero-dimensional grevlex basis for a system exported by
``recover_lattice_foundry_ns0024_mw4_family_resolved_modp.sage``.  Unlike the
older P4 extractor, every surface and MW3 coordinate is read from the same
residue field.  The script scans all decoded points, replays the four section
identities, absolute component marking, and full Gram, and emits one compact
point accepted by the edge-1 adapter.

Factoring/selecting an irreducible eliminant is deliberately upstream: pass
its irreducible modulus here.  Thus the residue degree is arbitrary, not
hard-coded to two.
"""

import argparse
import json
from collections import Counter
from itertools import product
from pathlib import Path

from sage.all import (
    EllipticCurve,
    GF,
    PolynomialRing,
    QQ,
    ZZ,
    binomial,
    ceil,
    divisors,
    identity_matrix,
    matrix,
    sage_eval,
    sqrt,
    vector,
)


RESOLVED_PROFILES = ((6, 0, 0), (2, 1, 1), (4, 2, 0), (6, 4, 3))
RESOLVED_GRAM = (
    (-2, 1, 1, 1),
    (1, -2, 0, 2),
    (1, 0, -2, 2),
    (1, 2, 2, -2),
)
SOURCE_PO = (0, 0, 0, 1)
FIBRE_ORDERS = (7, 5, 4)


def fibre_correction(left, right, order):
    return QQ(min(left, right) * (order - max(left, right))) / order


SOURCE_HEIGHT_GRAM = matrix(
    QQ,
    4,
    4,
    lambda left, right: (
        2
        + SOURCE_PO[left]
        + SOURCE_PO[right]
        - RESOLVED_GRAM[left][right]
        - sum(
            fibre_correction(
                RESOLVED_PROFILES[left][place],
                RESOLVED_PROFILES[right][place],
                FIBRE_ORDERS[place],
            )
            for place in range(3)
        )
    ),
)
if SOURCE_HEIGHT_GRAM.det() != QQ(95) / 14:
    raise ArithmeticError("resolved source height Gram no longer has regulator 95/14")


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--system", type=Path, required=True, help="exported msolve input")
parser.add_argument("--gb", type=Path, required=True, help="msolve -g 2 output")
parser.add_argument("--prime", type=int, required=True)
parser.add_argument("--generator", default="z")
parser.add_argument("--modulus", required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()

prime = ZZ(args.prime)
if not prime.is_prime() or prime in (2, 3, 5, 7):
    raise SystemExit("require a good prime away from 2,3,5,7")
base = GF(prime)
modulus_ring = PolynomialRing(base, args.generator)
z0 = modulus_ring.gen()
modulus = modulus_ring(sage_eval(args.modulus, locals={args.generator: z0})).monic()
if modulus.degree() < 1 or not modulus.is_irreducible():
    raise SystemExit("--modulus must be irreducible over the prime field")
degree = modulus.degree()
extension = base if degree == 1 else GF(prime**degree, args.generator, modulus=modulus)
z = None if degree == 1 else extension.gen()

system_path = args.system.resolve()
gb_path = args.gb.resolve()
names = tuple(system_path.read_text().splitlines()[0].split(","))
required_names = {
    *(f"a{index}" for index in range(1, 7)),
    "r1", "ri", "c",
    *(f"p1x{index}" for index in range(1, 5)),
    *(f"p1y{index}" for index in range(1, 7)),
    "p2x2", "p2y2", "p2y3", "p2y4",
    "p3y3", "p3y4",
    *(f"p4x{index}" for index in range(1, 5)),
    *(f"p4y{index}" for index in range(1, 8)),
}
if not required_names.issubset(names):
    raise SystemExit("system is not a joint resolved-depth13 MW4 export")

ring = PolynomialRing(extension, names=names, order="degrevlex")
v = ring.gens_dict()
text = gb_path.read_text()
if "[" not in text or "]" not in text:
    raise SystemExit("GB file does not contain an msolve basis list")
body = text[text.index("[") + 1:text.rindex("]")]
polynomials = [
    ring(item.replace("^", "**"))
    for item in body.split(",\n")
    if item.strip()
]
ideal = ring.ideal(polynomials)
if ideal.dimension() != 0:
    raise SystemExit("joint ideal is not zero-dimensional")
quotient_dimension = ideal.vector_space_dimension()
print(
    "NS0024JOINTGB|p={}|degree={}|basis={}|quotient={}".format(
        prime, degree, len(polynomials), quotient_dimension
    ),
    flush=True,
)
solutions = ideal.variety(ring=extension)
print("NS0024JOINTVARIETY|solutions={}".format(len(solutions)), flush=True)

function_ring = PolynomialRing(extension, "t")
t = function_ring.gen()
function_field = function_ring.fraction_field()


def formal_center(a_jet, root, precision):
    center = [extension(root)]
    root_inverse = 1 / extension(root)
    for local_degree in range(1, precision):
        known = sum(
            center[left] * center[local_degree - left]
            for left in range(1, local_degree)
        )
        center.append((-a_jet[local_degree] / 3 - known) * root_inverse / 2)
    return center


def cube_coefficients(values, precision):
    local = PolynomialRing(extension, "s")
    s = local.gen()
    cube = 2 * local(values) ** 3
    return [extension(cube[index]) for index in range(precision)]


def finite_value(value, support):
    numerator, denominator = value.numerator(), value.denominator()
    return None if denominator(support) == 0 else numerator(support) / denominator(support)


def encoded(value):
    if degree == 1:
        return int(extension(value))
    coefficients = list(extension(value).polynomial())
    return [
        int(coefficients[index]) if index < len(coefficients) else 0
        for index in range(degree)
    ]


def frobenius_polynomial(polynomial, power):
    return function_ring(
        [extension(coefficient).frobenius(power) for coefficient in polynomial]
    )


def frobenius_function(value, power):
    return function_field(
        frobenius_polynomial(value.numerator(), power)
        / frobenius_polynomial(value.denominator(), power)
    )


def frobenius_point(point, power, curve):
    if point.is_zero():
        return point
    return curve(
        frobenius_function(point[0], power),
        frobenius_function(point[1], power),
    )


def polynomial_tuple_fixed(polynomials, power):
    return all(frobenius_polynomial(polynomial, power) == polynomial for polynomial in polynomials)


def point_tuple_fixed(points, power, curve):
    return all(frobenius_point(point, power, curve) == point for point in points)


def first_orbit_size(predicate):
    for orbit_size in divisors(degree):
        if predicate(int(orbit_size)):
            return int(orbit_size)
    raise ArithmeticError("finite-field tuple did not close under full Frobenius orbit")


def vectors_of_norm(target_norm):
    inverse = SOURCE_HEIGHT_GRAM.inverse()
    bounds = [
        int(ceil(sqrt(target_norm * inverse[index, index])))
        for index in range(4)
    ]
    for entries in product(*(range(-bound, bound + 1) for bound in bounds)):
        candidate = vector(ZZ, entries)
        if candidate * SOURCE_HEIGHT_GRAM * candidate == target_norm:
            yield candidate


def frobenius_mw4_audit(A, B, r1, ri, points, curve):
    """Separate the surface orbit from the relative action on marked MW4."""
    surface_orbit = first_orbit_size(
        lambda power: polynomial_tuple_fixed((A, B), power)
    )
    section_marking_orbit = first_orbit_size(
        lambda power: (
            polynomial_tuple_fixed((A, B), power)
            and point_tuple_fixed(points, power, curve)
        )
    )
    resolved_marking_orbit = first_orbit_size(
        lambda power: (
            polynomial_tuple_fixed((A, B), power)
            and extension(r1).frobenius(power) == r1
            and extension(ri).frobenius(power) == ri
            and point_tuple_fixed(points, power, curve)
        )
    )
    if section_marking_orbit % surface_orbit:
        raise ArithmeticError("section orbit is not divisible by the surface orbit")
    if resolved_marking_orbit % section_marking_orbit:
        raise ArithmeticError("resolved marking orbit is not divisible by the section orbit")

    relative_order = section_marking_orbit // surface_orbit
    columns = []
    failures = []
    for point_index, point in enumerate(points):
        conjugate = frobenius_point(point, surface_orbit, curve)
        matches = []
        for coordinates in vectors_of_norm(SOURCE_HEIGHT_GRAM[point_index, point_index]):
            candidate = sum(
                (ZZ(coordinates[index]) * points[index] for index in range(4)),
                curve(0),
            )
            if candidate == conjugate:
                matches.append(coordinates)
        if len(matches) != 1:
            failures.append(
                {
                    "section": point_index + 1,
                    "matching_lattice_vectors": [list(map(int, item)) for item in matches],
                }
            )
        else:
            columns.append(matches[0])

    audit = {
        "prime": int(prime),
        "residue_degree": int(degree),
        "surface_frobenius_orbit_size": surface_orbit,
        "section_marking_frobenius_orbit_size": section_marking_orbit,
        "resolved_oriented_marking_frobenius_orbit_size": resolved_marking_orbit,
        "relative_section_degree_over_surface_field": relative_order,
        "relative_orientation_degree_over_section_field": (
            resolved_marking_orbit // section_marking_orbit
        ),
        "source_height_gram": [list(map(str, row)) for row in SOURCE_HEIGHT_GRAM.rows()],
        "source_height_gram_determinant": str(SOURCE_HEIGHT_GRAM.det()),
        "action_closes_on_marked_mw4": not failures,
        "action_failures": failures,
        "interpretation_boundary": (
            "Only Frobenius powers fixing the unmarked surface define an action on this "
            "MW4 fibre. A surface orbit of size greater than one is inconclusive for "
            "descent of a characteristic-zero family over Q. Repeated good-prime fixed-rank "
            "data are an arithmetic-realizability gate, not a proof that NS=NS_Q."
        ),
    }
    if failures:
        audit["status"] = "WARN_RELATIVE_FROBENIUS_DOES_NOT_CLOSE_ON_MARKED_MW4"
        return audit

    action = matrix(ZZ, 4, 4)
    for column, coordinates in enumerate(columns):
        action.set_column(column, coordinates)
    identity = identity_matrix(ZZ, 4)
    isometry = action.transpose() * SOURCE_HEIGHT_GRAM * action == SOURCE_HEIGHT_GRAM
    order_check = action**relative_order == identity
    fixed_rank = 4 - (action - identity).rank()
    audit.update(
        {
            "relative_frobenius_action_columns": [
                list(map(int, action.column(index))) for index in range(4)
            ],
            "relative_frobenius_action_determinant": int(action.det()),
            "height_isometry": bool(isometry),
            "relative_order_check": bool(order_check),
            "relative_fixed_marked_mw4_rank": int(fixed_rank),
            "prime_field_fixed_marked_mw4_rank": (
                int(fixed_rank) if surface_orbit == 1 else None
            ),
        }
    )
    if not isometry or not order_check or abs(action.det()) != 1:
        audit["status"] = "FAIL_INVALID_RELATIVE_FROBENIUS_ACTION"
    elif surface_orbit != 1:
        audit["status"] = "INCONCLUSIVE_SURFACE_NOT_FIXED_OVER_PRIME_FIELD"
    elif fixed_rank == 4:
        audit["status"] = "PASS_FULL_MW4_FIXED_AT_THIS_REDUCTION"
    else:
        audit["status"] = "WARN_PRIME_FIELD_FIXED_MW_RANK_LT4"
    return audit


accepted = []
failure_counts = Counter()
for solution in solutions:
    value = lambda name: extension(solution[v[name]])
    r1, ri = value("r1"), value("ri")
    a = [extension(-3)] + [value(f"a{index}") for index in range(1, 7)]
    a += [None, -3 * ri**2]
    a[7] = -3 * r1**2 - sum(a[:7]) - a[8]
    A = function_ring(a)

    center0 = formal_center(a[:8], 1, 8)
    a_at_one = [
        sum(a[index] * binomial(index, jet) for index in range(jet, 9))
        for jet in range(6)
    ]
    center1 = formal_center(a_at_one, r1, 6)
    a_at_infinity = [a[8 - index] for index in range(5)]
    center_infinity = formal_center(a_at_infinity, ri, 5)
    b_bottom = cube_coefficients(center0, 7)
    b_at_one = cube_coefficients(center1, 5)
    b_top_reversed = cube_coefficients(center_infinity, 4)
    b = b_bottom + [None, None] + list(reversed(b_top_reversed))
    known_value = sum(b[index] for index in range(13) if index not in (7, 8))
    known_derivative = sum(
        index * b[index] for index in range(13) if index not in (7, 8)
    )
    value_rhs = b_at_one[0] - known_value
    derivative_rhs = b_at_one[1] - known_derivative
    b[8] = derivative_rhs - 7 * value_rhs
    b[7] = 8 * value_rhs - derivative_rhs
    B = function_ring(b)
    curve = EllipticCurve(function_field, [0, 0, 0, A, B])

    X1 = function_ring([1] + [value(f"p1x{index}") for index in range(1, 5)])
    Y1 = function_ring([0] + [value(f"p1y{index}") for index in range(1, 7)])
    x21 = center0[1]
    x22 = value("p2x2")
    x23 = r1 - 1 - ri - x21 - x22
    X2 = function_ring([1, x21, x22, x23, ri])
    y22, y23, y24 = (value(f"p2y{index}") for index in range(2, 5))
    Y2 = function_ring([0, 0, y22, y23, y24, -y22 - y23 - y24])

    known_at_one = 1 + center0[1] + center0[2]
    q3_value_rhs = r1 - known_at_one
    q3_derivative_rhs = center1[1] - center0[1] - 2 * center0[2]
    x34 = q3_derivative_rhs - 3 * q3_value_rhs
    x33 = 4 * q3_value_rhs - q3_derivative_rhs
    X3 = function_ring([1, center0[1], center0[2], x33, x34])
    y33, y34 = value("p3y3"), value("p3y4")
    Y3 = function_ring([0, 0, 0, y33, y34, -3 * y33 - 2 * y34, 2 * y33 + y34])

    c = value("c")
    h = t - c
    p4x = [value(f"p4x{index}") for index in range(1, 5)]
    x45 = r1 * (1 - c)**2 - c**2 - ri - sum(p4x)
    X4 = function_ring([c**2] + p4x + [x45, ri])
    p4y = [value(f"p4y{index}") for index in range(1, 8)]
    Y4 = function_ring([0] + p4y + [-sum(p4y)])

    identities = (
        Y1**2 == X1**3 + A * X1 + B,
        Y2**2 == X2**3 + A * X2 + B,
        Y3**2 == X3**3 + A * X3 + B,
        Y4**2 == X4**3 + A * X4 * h**4 + B * h**6,
    )
    if not all(identities):
        failure_counts["section_identity"] += 1
        continue
    points = [
        curve(function_field(X1), function_field(Y1)),
        curve(function_field(X2), function_field(Y2)),
        curve(function_field(X3), function_field(Y3)),
        curve(function_field(X4 / h**2), function_field(Y4 / h**3)),
    ]

    discriminant = function_ring(-16 * (4 * A**3 + 27 * B**2))
    orders = (
        discriminant.valuation(t),
        discriminant.valuation(t - 1),
        24 - discriminant.degree(),
    )
    if tuple(map(int, orders)) != (7, 5, 4):
        failure_counts["fibre_orders"] += 1
        continue
    residual = function_ring(discriminant // (t**7 * (t - 1)**5))
    if residual.degree() != 8 or residual.gcd(residual.derivative()).degree() != 0:
        failure_counts["residual_discriminant"] += 1
        continue

    def hits_node(point, fibre):
        if point.is_zero():
            return False
        if fibre < 2:
            x_value = finite_value(point[0], extension(fibre))
            y_value = finite_value(point[1], extension(fibre))
            return x_value == (r1 if fibre else 1) and y_value == 0
        x_num, x_den = point[0].numerator(), point[0].denominator()
        y_num, y_den = point[1].numerator(), point[1].denominator()
        x_excess = x_num.degree() - x_den.degree()
        y_excess = y_num.degree() - y_den.degree()
        x_value = (
            0 if x_excess < 4 else
            x_num.leading_coefficient() / x_den.leading_coefficient()
            if x_excess == 4 else None
        )
        y_value = (
            0 if y_excess < 6 else
            y_num.leading_coefficient() / y_den.leading_coefficient()
            if y_excess == 6 else None
        )
        return x_value == ri and y_value == 0

    def component_label(point, reference, order, fibre):
        labels = [
            multiplier for multiplier in range(order)
            if not hits_node(point - multiplier * reference, fibre)
        ]
        return labels[0] if len(labels) == 1 else -1

    def intersection(left, right):
        difference = left - right
        if difference.is_zero():
            return -2
        numerator, denominator = difference[0].numerator(), difference[0].denominator()
        intersection_degree = denominator.degree() + max(
            0, numerator.degree() - denominator.degree() - 4
        )
        if intersection_degree % 2:
            return None
        return intersection_degree // 2

    relative_profiles = tuple(
        tuple(
            component_label(point, points[3], order, fibre)
            for fibre, order in enumerate((7, 5, 4))
        )
        for point in points
    )
    p4_profile = RESOLVED_PROFILES[3]
    profiles = tuple(
        tuple(
            (relative_profiles[index][fibre] * p4_profile[fibre]) % order
            for fibre, order in enumerate((7, 5, 4))
        )
        for index in range(4)
    )
    gram = tuple(tuple(intersection(left, right) for right in points) for left in points)
    if profiles != RESOLVED_PROFILES:
        failure_counts["component_profile"] += 1
        continue
    if gram != RESOLVED_GRAM:
        failure_counts["gram"] += 1
        continue
    galois_audit = frobenius_mw4_audit(A, B, r1, ri, points, curve)
    accepted.append(
        (A, B, r1, ri, X1, Y1, X2, Y2, X3, Y3, X4, Y4, h, galois_audit)
    )

print(
    "NS0024JOINTMARKING|accepted={}|failures={}".format(
        len(accepted), dict(failure_counts)
    ),
    flush=True,
)
if not accepted:
    raise SystemExit("no decoded joint point has the exact resolved MW4 marking")

A, B, r1, ri, X1, Y1, X2, Y2, X3, Y3, X4, Y4, h, galois_audit = accepted[0]
payload = {
    "schema": "elkies-k3.lattice-foundry-ns0024-mw4-point-modp.v1",
    "status": "PASS_EXACT_MW4_MARKED_POINT_OVER_FINITE_FIELD",
    "prime": int(prime),
    "basis_marking": "resolved_component_depth_recommendation",
    "source": {
        "A_coefficients_low_to_high": [encoded(entry) for entry in A],
        "B_coefficients_low_to_high": [encoded(entry) for entry in B],
        "r1": encoded(r1),
        "ri": encoded(ri),
        "sections": {
            "P1": {
                "X_coefficients_low_to_high": [encoded(entry) for entry in X1],
                "Y_coefficients_low_to_high": [encoded(entry) for entry in Y1],
            },
            "P2": {
                "X_coefficients_low_to_high": [encoded(entry) for entry in X2],
                "Y_coefficients_low_to_high": [encoded(entry) for entry in Y2],
            },
            "P3": {
                "X_coefficients_low_to_high": [encoded(entry) for entry in X3],
                "Y_coefficients_low_to_high": [encoded(entry) for entry in Y3],
            },
        },
    },
    "P4": {
        "X_coefficients_low_to_high": [encoded(entry) for entry in X4],
        "Y_coefficients_low_to_high": [encoded(entry) for entry in Y4],
        "H_coefficients_low_to_high": [encoded(entry) for entry in h],
    },
    "component_profiles_I7_I5_I4": [list(row) for row in RESOLVED_PROFILES],
    "section_intersection_gram": [list(row) for row in RESOLVED_GRAM],
    "arithmetic_realizability": galois_audit,
    "residue_algebra_dimension": int(quotient_dimension),
    "decoded_solutions": len(solutions),
    "exact_marking_solutions": len(accepted),
    "inputs": {"system": str(system_path), "groebner_basis": str(gb_path)},
    "proof_boundary": (
        "This is an exact isolated finite-field marked point. It does not prove "
        "a positive-dimensional marked family, a characteristic-zero lift, NS=NS_Q, "
        "or generic Mordell-Weil rank 4 (and hence not rank 17) over Q(t)."
    ),
}
if degree > 1:
    payload["extension"] = {"generator": args.generator, "modulus": str(modulus)}
args.output.resolve().write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "NS0024JOINTPOINT|p={}|degree={}|output={}|status=PASS".format(
        prime, degree, args.output.resolve()
    ),
    flush=True,
)
