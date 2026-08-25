#!/usr/bin/env sage-python
"""Select the q8/orbit376 horizontal from an inherited-P1 Abel trace mod p.

The exact fourfold-height audit leaves eight embeddings of the equation-side
rank-eight subgroup in the C8-pointed marked MW9 lattice.  For each embedding,
express the marked difference between the inherited-P1 trace and the q8
horizontal using B0,...,B7 and the saturated C8-opposite point.  Reduce the
exact points modulo the selected good prime and use the certified P.O=4 degree
fingerprint to select the unique embedding.

Only exact elliptic-curve group law and finite lattice arithmetic are used; no
Groebner basis is used.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q4o164-compact-weierstrass-qq.json"
BASIS = LOCAL / "q4o164-integral-basis-qq.json"
AUDIT = LOCAL / "q4o164-integral-basis-height-gram-audit-qq.json"
C8 = LOCAL / "q4o164-c8-equation-marking-qq.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=131)
parser.add_argument("--trace", type=Path)
parser.add_argument(
    "--output", type=Path,
)
args = parser.parse_args()
PRIME = args.prime
TRACE = args.trace or LOCAL / f"q4o164-inherited-p1-abel-trace-section-mod{PRIME}.json"
TRACE = TRACE if TRACE.is_absolute() else ROOT / TRACE
OUTPUT = args.output or LOCAL / f"q4o164-q8o376-horizontal-from-abel-trace-mod{PRIME}.json"
OUTPUT = OUTPUT if OUTPUT.is_absolute() else ROOT / OUTPUT
INPUTS = (MODEL, BASIS, AUDIT, TRACE, C8)
started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


model = json.loads(MODEL.read_text())
basis = json.loads(BASIS.read_text())
audit = json.loads(AUDIT.read_text())
trace = json.loads(TRACE.read_text())
c8 = json.loads(C8.read_text())
assert model["status"] == "PASS_EXACT_QQ_Q4O164_COMPACT_WEIERSTRASS_NORMALIZATION"
assert basis["status"] == "PASS_EXACT_QQ_Q4O164_PAIR_NODE_SECTION_SUBGROUP_RANK8"
assert audit["status"] == "PASS_EXACT_QQ_Q4O164_FOURFOLD_HEIGHT_GRAM_AND_C8_MARKED_EMBEDDING_CENSUS"
assert trace["status"] == "PASS_MODP_Q4O164_INHERITED_P1_DEGREE7_ABEL_TRACE_SECTION"
assert trace["prime"] == PRIME
assert c8["status"] == "PASS_EXACT_QQ_Q4O164_C8_EQUATION_MARKING"

F = GF(PRIME)
R = PolynomialRing(F, "t")
t = R.gen()
K = R.fraction_field()


def reduce_qq(value):
    value = QQ(value)
    return F(value.numerator()) / F(value.denominator())


def polynomial(values):
    return R([reduce_qq(value) for value in values])


def rational_coordinate(record):
    return K(polynomial(record["numerator_coefficients_low_to_high"])) / K(
        polynomial(record["denominator_coefficients_low_to_high"])
    )


def substitute(function, scalar):
    return K(R(function.numerator()(scalar * t))) / K(R(function.denominator()(scalar * t)))


A = polynomial(model["compact_model"]["A_coefficients_low_to_high"])
B = polynomial(model["compact_model"]["B_coefficients_low_to_high"])
base_scale = reduce_qq(model["exact_coordinate_change"]["c"])
xy_scale = reduce_qq(model["exact_coordinate_change"]["s"])


def checked_point(x_coordinate, y_coordinate):
    """Use an affine pair after a literal equation check."""
    x_coordinate = K(x_coordinate)
    y_coordinate = K(y_coordinate)
    assert y_coordinate**2 == x_coordinate**3 + K(A) * x_coordinate + K(B)
    return (x_coordinate, y_coordinate)


def point_neg(point):
    return None if point is None else (point[0], -point[1])


def point_add(left, right):
    """Short-Weierstrass group law over K, with None as the zero point."""
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2:
        if y1 == -y2:
            return None
        slope = (3 * x1**2 + K(A)) / (2 * y1)
    else:
        slope = (y2 - y1) / (x2 - x1)
    x3 = slope**2 - x1 - x2
    y3 = slope * (x1 - x3) - y1
    result = (x3, y3)
    assert y3**2 == x3**3 + K(A) * x3 + K(B)
    return result


def point_mul(coefficient, point):
    coefficient = ZZ(coefficient)
    if coefficient < 0:
        return point_mul(-coefficient, point_neg(point))
    result = None
    addend = point
    while coefficient:
        if coefficient & 1:
            result = point_add(result, addend)
        addend = point_add(addend, addend)
        coefficient >>= 1
    return result


def old_to_compact(x_coordinate, y_coordinate):
    return checked_point(
        substitute(x_coordinate, base_scale) / xy_scale**2,
        substitute(y_coordinate, base_scale) / xy_scale**3,
    )


trace_record = trace["trace_section"]
trace_point = old_to_compact(
    K(R(trace_record["x"]["numerator_coefficients_low_to_high"]))
    / K(R(trace_record["x"]["denominator_coefficients_low_to_high"])),
    K(R(trace_record["y"]["numerator_coefficients_low_to_high"]))
    / K(R(trace_record["y"]["denominator_coefficients_low_to_high"])),
)
c8_record = c8["opposite_constant_support_section"]
c8_opposite = old_to_compact(
    rational_coordinate(c8_record["x"]), rational_coordinate(c8_record["y"]),
)
basis_points = [
    checked_point(
        K(polynomial(record["x_coefficients_low_to_high"])),
        K(polynomial(record["y_coefficients_low_to_high"])),
    )
    for record in basis["resolved_hensel"]["sections"]
]

# This saturated relation was found by fourfold pairings and is checked again
# after good reduction.  It lets rational B-coordinates be represented by an
# actual equation point without dividing on the elliptic curve.
c8_relation = vector(ZZ, [-2, -3, -4, 3, -2, 2, -1, -2])
relation_point = None
for coefficient, point in zip(c8_relation, basis_points):
    relation_point = point_add(relation_point, point_mul(coefficient, point))
assert point_mul(3, c8_opposite) == relation_point
c8_rational_basis_coordinates = vector(
    QQ, [QQ(-2) / 3, -1, QQ(-4) / 3, 1, QQ(-2) / 3, QQ(2) / 3, QQ(-1) / 3, QQ(-2) / 3]
)

marked = audit["marked_embedding_enumeration"]
q8_tail = vector(QQ, marked["q8_horizontal_marked_MW9_tail"])
trace_tail = vector(QQ, marked["inherited_P1_trace_marked_MW9_tail"])
full_difference = q8_tail - trace_tail


def coordinate_record(value):
    return {
        "numerator_coefficients_low_to_high": list(map(int, value.numerator().list())),
        "denominator_coefficients_low_to_high": list(map(int, value.denominator().list())),
        "degrees_numerator_denominator": [int(value.numerator().degree()), int(value.denominator().degree())],
    }


candidates = []
for embedding_record in marked["embeddings"]:
    if not embedding_record["compatible_with_first_seven_stored_profiles_up_to_fibre_symmetry"]:
        continue
    embedding = matrix(ZZ, embedding_record["rows_B0_through_B7_in_marked_MW9"])
    rational_word = embedding.transpose().solve_right(full_difference)
    saturated_words = []
    for c8_coefficient in range(-3, 4):
        integral_word = rational_word - c8_coefficient * c8_rational_basis_coordinates
        if all(value.denominator() == 1 for value in integral_word):
            saturated_words.append((c8_coefficient, vector(ZZ, integral_word)))
    assert saturated_words
    c8_coefficient, integral_word = min(
        saturated_words, key=lambda item: (abs(item[0]), item[0], tuple(item[1]))
    )
    candidate = point_add(trace_point, point_mul(c8_coefficient, c8_opposite))
    for coefficient, point in zip(integral_word, basis_points):
        candidate = point_add(candidate, point_mul(coefficient, point))
    assert candidate is not None
    x_coordinate, y_coordinate = candidate[0], candidate[1]
    degrees = (
        x_coordinate.numerator().degree(), x_coordinate.denominator().degree(),
        y_coordinate.numerator().degree(), y_coordinate.denominator().degree(),
    )
    pole_degree = max(degrees[1], degrees[0] - 4)
    candidates.append({
        "embedding_index": embedding_record["embedding_index"],
        "c8_opposite_coefficient": c8_coefficient,
        "integral_B0_through_B7_word": list(map(int, integral_word)),
        "all_equivalent_saturated_words": [
            {"c8_opposite_coefficient": coefficient, "integral_B_word": list(map(int, word))}
            for coefficient, word in saturated_words
        ],
        "coordinate_degrees_x_num_x_den_y_num_y_den": list(map(int, degrees)),
        "P_dot_O_from_compact_pole_degree": int(pole_degree // 2),
        "exact_modp_weierstrass_identity": True,
        "x": coordinate_record(x_coordinate),
        "y": coordinate_record(y_coordinate),
    })

selected = [
    record for record in candidates
    if record["coordinate_degrees_x_num_x_den_y_num_y_den"] == [12, 8, 18, 12]
    and record["P_dot_O_from_compact_pole_degree"] == 4
]
assert len(candidates) == 8
assert len(selected) == 1
selected = selected[0]
assert selected["c8_opposite_coefficient"] == -1
assert selected["integral_B0_through_B7_word"] == [-1, 2, 1, -3, -1, -2, 0, 1]

payload = {
    "schema": "elkies-k3.q4o164-q8o376-horizontal-from-abel-trace-modp.v1",
    "status": "PASS_EXACT_MODP_Q4O164_Q8O376_HORIZONTAL_FROM_ABEL_TRACE",
    "prime": PRIME,
    "selected_identity": (
        "H=T-C8opp-B0+2*B1+B2-3*B3-B4-2*B5+B7"
    ),
    "selected": selected,
    "embedding_candidates": candidates,
    "selection_gate": {
        "marked_q8_P_dot_O": 4,
        "required_compact_degrees_x_num_x_den_y_num_y_den": [12, 8, 18, 12],
        "candidate_count": len(candidates),
        "unique_degree_match": True,
    },
    "saturation_relation": {
        "identity": "3*C8opp=-2*B0-3*B1-4*B2+3*B3-2*B4+2*B5-B6-2*B7",
        "literal_modp_group_law_identity": True,
    },
    "method": {
        "large_Groebner_required": False,
        "exact_modp_group_law": True,
        "finite_marked_embedding_census": True,
        "runtime_seconds": time.monotonic() - started,
    },
    "proof_boundary": (
        f"The inherited-P1 Abel trace plus exact equation points produces a unique mod-{PRIME} "
        "section with the certified q8/orbit376 P.O=4 degree fingerprint. This selects the "
        "marked embedding at the pinned good prime and gives the complete modular horizontal. "
        "A characteristic-zero reconstruction of the Abel trace, followed by literal QQ(t) "
        "group law and the resolved q8 RR calculation, remains open."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in INPUTS},
    },
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4O164Q8ABEL|candidates={}|selected={}|degrees={}|word={}|status={}|output={}".format(
        len(candidates), selected["embedding_index"],
        selected["coordinate_degrees_x_num_x_den_y_num_y_den"],
        payload["selected_identity"], payload["status"], OUTPUT,
    ),
    flush=True,
)
