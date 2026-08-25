#!/usr/bin/env sage -python
"""Find a full-rank integral-section subgroup on q4/orbit164 modulo 131.

For each pair among the finite I4, two finite I2 fibres, and the I4 at
infinity, impose the two corresponding node conditions on a degree-four x
polynomial.  The remaining three coefficients are searched exhaustively and
the cubic right side is tested for a polynomial square by coefficient
recursion.  Smooth-fibre specialization maps to E(k)/ell E(k) select a
genuinely independent rank-nine subset.  This is a finite-field construction
aid only.
"""

import argparse
import hashlib
import json
import time
from itertools import product
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, matrix


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q4o1584-physical-q4o164-rr-qq.json"
OUTPUT = LOCAL / "q4o164-integral-sections-mod131.json"
PRIME = 131

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--reuse-pool", action="store_true",
    help="reuse the 128 stored sections in the existing output and rerun only the rank selector",
)
args = parser.parse_args()

started = time.monotonic()


def log(stage, **fields):
    tail = "|".join(f"{key}={value}" for key, value in fields.items())
    print(
        f"Q4O164SECT131|stage={stage}|elapsed={time.monotonic()-started:.3f}"
        + (f"|{tail}" if tail else ""),
        flush=True,
    )


model = json.loads(MODEL.read_text())
assert model["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O164_2A3_2A1_RR_AND_JACOBIAN"

k = GF(PRIME)
R = PolynomialRing(k, "u")
u = R.gen()
K = R.fraction_field()
RX = PolynomialRing(k, "x")
x_variable = RX.gen()


def reduce_qq(value):
    value = QQ(value)
    return k(value.numerator()) / k(value.denominator())


A = R([reduce_qq(value) for value in model["child"]["minimal_A_coefficients_low_to_high"]])
B = R([reduce_qq(value) for value in model["child"]["minimal_B_coefficients_low_to_high"]])
E = EllipticCurve(K, [0, 0, 0, K(A), K(B)])
discriminant = -k(16) * (4 * A**3 + 27 * B**2)
assert discriminant.degree() == 20

finite_records = model["child"]["finite_reducible_fibres"]
finite_supports = [reduce_qq(item["support"]) for item in finite_records]
finite_orders = [int(item["delta_order"]) for item in finite_records]
assert sorted(finite_orders) == [2, 2, 4]
finite_nodes = []
for support in finite_supports:
    cubic = x_variable**3 + A(support) * x_variable + B(support)
    repeated = cubic.gcd(cubic.derivative())
    assert repeated.degree() == 1
    finite_nodes.append(-repeated[0] / repeated[1])

# At infinity use X=t^4*x, Y=t^6*y.  Since deg(A,B)=(8,12), the nodal
# cubic is determined by the leading coefficients.
A_infinity = A[8]
B_infinity = B[12]
cubic_infinity = x_variable**3 + A_infinity * x_variable + B_infinity
repeated_infinity = cubic_infinity.gcd(cubic_infinity.derivative())
assert repeated_infinity.degree() == 1
node_infinity = -repeated_infinity[0] / repeated_infinity[1]

supports = finite_supports + [None]
nodes = finite_nodes + [node_infinity]
orders = finite_orders + [4]
labels = [f"finite_{index}_{finite_records[index]['kodaira']}" for index in range(3)] + ["infinity_I4"]
log("LOAD", supports=",".join(labels), model_degrees="8/12/20")


def polynomial_square_roots(polynomial):
    polynomial = R(polynomial)
    if polynomial == 0:
        return (R.zero(),)
    shift = next(value for value in k if polynomial(value) != 0)
    shifted = polynomial(u + shift)
    constant = shifted[0]
    if not constant.is_square():
        return ()
    roots = []
    for first in constant.sqrt(all=True):
        coefficients = [first]
        for degree in range(1, 7):
            known = sum(
                (coefficients[left] * coefficients[degree-left]
                 for left in range(1, degree)),
                k.zero(),
            )
            coefficients.append((shifted[degree] - known) / (2 * first))
        candidate_shifted = R(coefficients)
        if candidate_shifted**2 == shifted:
            roots.append(candidate_shifted(u - shift))
    return tuple(roots)


def constraint_row(index):
    if supports[index] is None:
        return [k.zero()] * 4 + [k.one()]
    return [supports[index]**degree for degree in range(5)]


def search_pair(pair):
    constraint = matrix(k, [constraint_row(index) for index in pair])
    target = matrix(k, [[nodes[index] for index in pair]]).transpose()
    particular = constraint.solve_right(target).column(0)
    kernel = constraint.right_kernel().basis()
    assert len(kernel) == 3
    answers = []
    search_started = time.monotonic()
    for coefficients in product(k, repeat=3):
        vector = particular + sum(
            (coefficients[index] * kernel[index] for index in range(3)),
            particular.parent().zero(),
        )
        X = R(list(vector))
        for Y in polynomial_square_roots(X**3 + A * X + B):
            answers.append(E(K(X), K(Y)))
    log(
        "PAIR", pair=f"{pair[0]}{pair[1]}", tests=PRIME**3,
        sections=len(answers), seconds=f"{time.monotonic()-search_started:.3f}",
    )
    return answers


cached = json.loads(OUTPUT.read_text()) if args.reuse_pool and OUTPUT.exists() else None
if cached is not None:
    cached_rows = cached["integral_subgroup"]["all_sections"]
    assert len(cached_rows) == 128
    points = [
        E(K(R(row["x_coefficients_low_to_high"])), K(R(row["y_coefficients_low_to_high"])))
        for row in cached_rows
    ]
    pair_counts = cached["search"]["pairs"]
    log("REUSE", sections=len(points))
else:
    shells = {}
    for left in range(4):
        for right in range(left + 1, 4):
            shells[(left, right)] = search_pair((left, right))
    points = sorted(
        set(point for shell in shells.values() for point in shell),
        key=lambda point: (str(point[0]), str(point[1])),
    )
    pair_counts = {
        f"{left}{right}": len(shell)
        for (left, right), shell in shells.items()
    }
log("UNION", sections=len(points))


def hits_node(point, index):
    if point.is_zero():
        return False
    x_coordinate, y_coordinate = K(point[0]), K(point[1])
    if supports[index] is None:
        x_lead = x_coordinate.numerator()[4] if x_coordinate.denominator().degree() == 0 else None
        y_lead = y_coordinate.numerator()[6] if y_coordinate.denominator().degree() == 0 else None
        return x_lead == nodes[index] and y_lead == 0
    support = supports[index]
    return (
        x_coordinate.denominator()(support) != 0
        and y_coordinate.denominator()(support) != 0
        and x_coordinate(support) == nodes[index]
        and y_coordinate(support) == 0
    )


references = {}
for index, order in enumerate(orders):
    if order == 2:
        continue
    references[index] = next(
        (point for point in points if hits_node(2 * point, index)), None
    )


def component_label(point, index):
    order = orders[index]
    if order == 2:
        return int(hits_node(point, index))
    reference = references[index]
    if reference is None:
        # The searched subgroup can project only to {0,2} in this I4
        # component group.  Node specialization then distinguishes the two
        # labels without choosing an orientation.
        if any(hits_node(2 * candidate, index) for candidate in points):
            raise ArithmeticError("missing I4 reference despite an odd component label")
        return 2 if hits_node(point, index) else 0
    answers = [
        multiplier for multiplier in range(order)
        if not hits_node(point - multiplier * reference, index)
    ]
    assert len(answers) == 1
    return answers[0]


profiles = {point: tuple(component_label(point, index) for index in range(4)) for point in points}


def quotient_coordinates(curve, ell):
    """Return coordinates on E(k)/ell*E(k), by exhaustive finite-group cosets."""
    curve_points = curve.points()
    multiples = {ell * point for point in curve_points}
    unseen = set(curve_points)
    cosets = []
    representatives = []
    while unseen:
        representative = next(iter(unseen))
        coset = {representative + point for point in multiples}
        unseen -= coset
        representatives.append(representative)
        cosets.append(coset)
    if len(cosets) == 1:
        return None

    def coset_index(point):
        return next(index for index, coset in enumerate(cosets) if point in coset)

    zero_index = coset_index(curve(0))
    first_index = next(index for index in range(len(cosets)) if index != zero_index)
    first = representatives[first_index]
    first_span = {coset_index(multiplier * first): [multiplier] for multiplier in range(ell)}
    coordinates = first_span
    if len(cosets) == ell**2:
        second_index = next(index for index in range(len(cosets)) if index not in first_span)
        second = representatives[second_index]
        coordinates = {
            coset_index(left * first + right * second): [left, right]
            for left in range(ell) for right in range(ell)
        }
    assert len(coordinates) == len(cosets) in (ell, ell**2)
    return lambda point: coordinates[coset_index(point)]


# A relation among QQ lifts would reduce to a relation among these modular
# sections.  Select a genuinely independent rank-nine subset by evaluating at
# smooth fibres and accumulating maps to E(k)/ell E(k).  This avoids the
# invalid shortcut of reading intersections with O from denominators at bad
# fibres of the singular Weierstrass model.
basis_indices = None
specialization_record = None
best_rank = -1
best_indices = []
best_record = None
for ell in (2, 3, 5, 7, 11):
    signature_rows = []
    used_supports = []
    for support in k:
        a_value, b_value = A(support), B(support)
        if 4 * a_value**3 + 27 * b_value**2 == 0:
            continue
        curve = EllipticCurve(k, [0, 0, 0, a_value, b_value])
        coordinates = quotient_coordinates(curve, ell)
        if coordinates is None:
            continue
        point_coordinates = [
            coordinates(curve(point[0](support), point[1](support)))
            for point in points
        ]
        for coordinate_index in range(len(point_coordinates[0])):
            signature_rows.append([
                point_coordinates[point_index][coordinate_index]
                for point_index in range(len(points))
            ])
        used_supports.append(int(support))
        signature_matrix = matrix(GF(ell), signature_rows)
        current_rank = int(signature_matrix.rank())
        if current_rank > best_rank:
            best_rank = current_rank
            best_indices = list(map(int, signature_matrix.pivots()[:current_rank]))
            best_record = {
                "quotient_prime": ell,
                "smooth_fibre_supports": list(used_supports),
                "signature_rows": len(signature_rows),
                "signature_rank": current_rank,
                "basis_indices_in_union": best_indices,
                "basis_signature_matrix": [
                    [int(signature_matrix[row, column]) for column in best_indices]
                    for row in range(signature_matrix.nrows())
                ],
            }
        if signature_matrix.rank() >= 9:
            basis_indices = list(map(int, signature_matrix.pivots()[:9]))
            specialization_record = {
                "quotient_prime": ell,
                "smooth_fibre_supports": used_supports,
                "signature_rows": len(signature_rows),
                "signature_rank": int(signature_matrix.rank()),
                "basis_indices_in_union": basis_indices,
                "basis_signature_matrix": [
                    [int(signature_matrix[row, column]) for column in basis_indices]
                    for row in range(signature_matrix.nrows())
                ],
            }
            break
    if basis_indices is not None:
        break

basis_indices = basis_indices if basis_indices is not None else best_indices
specialization_record = specialization_record if specialization_record is not None else best_record
basis = [points[index] for index in basis_indices]
status = (
    "PASS_MOD131_Q4O164_INTEGRAL_SECTION_SUBGROUP_RANK9"
    if specialization_record["signature_rank"] == 9 else
    "PASS_MOD131_Q4O164_PAIR_NODE_SECTION_SUBGROUP_RANK8"
    if specialization_record["signature_rank"] == 8 else
    "INCOMPLETE_MOD131_Q4O164_INTEGRAL_SECTION_SUBGROUP"
)
log(
    "SPECIALIZATION", rank=(specialization_record or {}).get("signature_rank", 0),
    basis=len(basis), ell=(specialization_record or {}).get("quotient_prime"), status=status,
)


def point_record(point):
    return {
        "x_coefficients_low_to_high": [int(value) for value in point[0].numerator().list()],
        "y_coefficients_low_to_high": [int(value) for value in point[1].numerator().list()],
        "component_profile": list(profiles[point]),
    }


payload = {
    "schema": "elkies-k3.q4o164-integral-sections-mod131.v2",
    "status": status,
    "prime": PRIME,
    "fibres": [
        {"label": labels[index], "order": orders[index], "node": int(nodes[index])}
        for index in range(4)
    ],
    "search": {
        "tests_per_pair": PRIME**3,
        "pairs": pair_counts,
        "union_sections": len(points),
        "large_Groebner_required": False,
    },
    "integral_subgroup": {
        "independence_certificate": specialization_record,
        "basis": [point_record(point) for point in basis],
        "all_sections": [point_record(point) for point in points],
    },
    "method": {
        "proof_scope": "finite-field construction aid only",
        "runtime_seconds": time.monotonic() - started,
    },
    "inputs": {
        "paths": [str(MODEL.relative_to(ROOT))],
        "sha256": {str(MODEL.relative_to(ROOT)): hashlib.sha256(MODEL.read_bytes()).hexdigest()},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
log("DONE", output=OUTPUT, status=status)
