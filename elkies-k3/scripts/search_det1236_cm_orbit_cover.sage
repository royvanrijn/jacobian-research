#!/usr/bin/env sage
"""Search the discriminant -1236 cubic field for the missing cover divisor.

Let

    K = QQ[a]/(a^3-a^2+4*a+12),
    E = 618f1,
    G = (10,-29),
    A = (-3*a^2-9*a-20, 15*a^2-57*a-155).

The Galois-orbit sum of A is -16*G.  Hence the orbit sum of
P=m*G+n*A is (3*m-16*n)*G.  If the four branch points of a genus-three
cover E are r*G (r=+/-3) and the three conjugates of P, its L(6O)
representative has divisor

    (r*G) + Orbit(P) + 2*Q - 6*O,
    Q = -(r+3*m-16*n)/2 * G.

This script constructs that function by exact linear algebra and compares
its reductions with the unique local V4-compatible L(6O) classes at
p=11,13,17,19 found by explore_det1236_double_cover_local_gate.sage.
It is a theorem-directed diagnostic: a bounded coefficient search is not
an exclusion if no match is found.

Replay:
    sage elkies-k3/scripts/search_det1236_cm_orbit_cover.sage 30
"""

import json
import sys


BOUND = 20 if len(sys.argv) == 1 else ZZ(sys.argv[1])
assert BOUND >= 1

POLYNOMIAL_RING = PolynomialRing(QQ, "u")
u = POLYNOMIAL_RING.gen()
K = NumberField(u**3-u**2+4*u+12, "a")
a = K.gen()
E = EllipticCurve(K, [1, 0, 0, -185, 1401])
G = E(10, -29)
A = E(-3*a**2-9*a-20, 15*a**2-57*a-155)

LOCAL_CLASSES = {
    -3: {
        11: [1, 0, 10, 4, 9, 9],
        13: [1, 9, 10, 9, 8, 3],
        17: [1, 2, 13, 10, 4, 2],
        19: [1, 16, 13, 17, 9, 11],
    },
    3: {
        11: [1, 1, 1, 6, 2, 9],
        13: [1, 12, 3, 1, 5, 3],
        17: [1, 6, 4, 6, 13, 2],
        19: [1, 3, 6, 8, 10, 11],
    },
}


def evaluation_row(point, length=6):
    x_value, y_value = point[0], point[1]
    return vector(
        point.curve().base_ring(),
        [1, x_value, y_value, x_value**2, x_value*y_value, x_value**3][
            :length
        ],
    )


def derivative_row(point):
    x_value, y_value = point[0], point[1]
    denominator = 2*y_value+x_value
    numerator = 3*x_value**2-185-y_value
    if not denominator:
        return vector(QQ, [0, 0, 1, 0, x_value, 0])
    return vector(
        QQ,
        [
            0,
            denominator,
            numerator,
            2*x_value*denominator,
            y_value*denominator+x_value*numerator,
            3*x_value**2*denominator,
        ],
    )


def coordinate_rows(k_row):
    rows = [[QQ(0) for _ in k_row] for _ in range(3)]
    for column, value in enumerate(k_row):
        coordinates = list(K(value))
        coordinates += [QQ(0)]*(3-len(coordinates))
        for row in range(3):
            rows[row][column] = coordinates[row]
    return rows


def cover_function(branch_multiple, m_value, n_value):
    trace_multiple = 3*m_value-16*n_value
    total = branch_multiple+trace_multiple
    if total % 2:
        return None
    q_multiple = -total//2
    rational_curve = EllipticCurve(QQ, [1, 0, 0, -185, 1401])
    rational_generator = rational_curve(10, -29)
    rational_branch = branch_multiple*rational_generator
    q_point = q_multiple*rational_generator
    cubic_point = m_value*G+n_value*A
    if cubic_point.is_zero() or cubic_point[0] in QQ:
        return None

    if q_point.is_zero():
        row_length = 4
        rows = [list(evaluation_row(rational_branch, row_length))]
    else:
        row_length = 6
        rows = [
            list(evaluation_row(rational_branch)),
            list(evaluation_row(q_point)),
            list(derivative_row(q_point)),
        ]
    rows += coordinate_rows(evaluation_row(cubic_point, row_length))
    kernel = matrix(QQ, rows).right_kernel().basis()
    if len(kernel) != 1:
        return None
    coefficients = list(kernel[0])
    coefficients += [QQ(0)]*(6-len(coefficients))
    pivot = next(value for value in coefficients if value)
    coefficients = [value/pivot for value in coefficients]
    return coefficients, q_multiple, cubic_point


def reduce_projectively(coefficients, prime):
    field = GF(prime)
    try:
        reduced = [field(value) for value in coefficients]
    except (ZeroDivisionError, TypeError, ValueError):
        return None
    pivot = next((value for value in reduced if value), None)
    if pivot is None:
        return None
    return [int(value/pivot) for value in reduced]


# Verify the orbit-sum calculation independently in the splitting field.
splitting_field = K.defining_polynomial().splitting_field("b")
roots = K.defining_polynomial().change_ring(splitting_field).roots(
    multiplicities=False
)
split_curve = EllipticCurve(splitting_field, [1, 0, 0, -185, 1401])
orbit_sum = split_curve(0)
for root in roots:
    orbit_sum += split_curve(
        -3*root**2-9*root-20,
        15*root**2-57*root-155,
    )
assert orbit_sum == -16*split_curve(10, -29)

tested = 0
matches = []
best_agreement = []
best_count = -1
for branch_multiple in (-3, 3):
    for m_value in range(-BOUND, BOUND+1):
        for n_value in range(-BOUND, BOUND+1):
            if n_value == 0:
                continue
            result = cover_function(branch_multiple, m_value, n_value)
            if result is None:
                continue
            tested += 1
            coefficients, q_multiple, cubic_point = result
            reductions = {
                int(prime): reduce_projectively(coefficients, prime)
                for prime in LOCAL_CLASSES[branch_multiple]
            }
            agreeing_primes = [
                int(prime)
                for prime, reduction in reductions.items()
                if reduction == LOCAL_CLASSES[branch_multiple][prime]
            ]
            agreement_count = len(agreeing_primes)
            summary = {
                "branch_multiple": int(branch_multiple),
                "m": int(m_value),
                "n": int(n_value),
                "q_multiple": int(q_multiple),
                "coefficients": [str(value) for value in coefficients],
                "agreeing_primes": agreeing_primes,
                "reductions": reductions,
            }
            if agreement_count > best_count:
                best_count = agreement_count
                best_agreement = [summary]
            elif agreement_count == best_count:
                best_agreement.append(summary)
            if agreement_count == len(LOCAL_CLASSES[branch_multiple]):
                matches.append(summary)

payload = {
    "schema": "elkies-k3.det1236-cm-orbit-cover-search.v1",
    "status": (
        "CANDIDATE_MATCH_FOUND" if matches else "BOUNDED_SEARCH_NO_MATCH"
    ),
    "coefficient_bound": int(BOUND),
    "cubic_field_polynomial": str(K.defining_polynomial()),
    "cubic_field_discriminant": int(K.discriminant()),
    "trace_of_A": "-16*G",
    "tested_principal_divisors": int(tested),
    "local_filter_primes": [11, 13, 17, 19],
    "matches": matches,
    "best_agreement_count": int(best_count),
    "best_agreement": best_agreement[:20],
    "boundary": (
        "no-match is not an exclusion: the coefficient box is bounded and "
        "the cubic field point has not been certified as the Shimura CM "
        "branch orbit; a match is likewise only a reconstruction candidate"
    ),
}
print(json.dumps(payload, indent=2, sort_keys=True, default=int))
