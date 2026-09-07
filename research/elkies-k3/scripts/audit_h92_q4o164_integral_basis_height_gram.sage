#!/usr/bin/env sage-python
"""Audit the q4/orbit164 equation basis by fourfold pole growth.

Raw specialization at a singular Weierstrass node does not distinguish all
resolved I_n components.  Clear the four component groups by replacing a
section P by 4P, recover its canonical height from the compact x-coordinate
pole degree, and polarize.  Then enumerate all isometric embeddings of the
resulting rank-eight height lattice in the certified marked MW9 lattice.

Only exact function-field group law and positive-definite lattice enumeration
are used; no Groebner basis is used.
"""

import argparse
import hashlib
import json
import time
from itertools import product
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
MODEL = LOCAL / "q4o164-compact-weierstrass-qq.json"
BASIS = LOCAL / "q4o164-integral-basis-qq.json"
MARKING = GENERATED / "elkies-k3-h3-q4o208-q4o1584-q4o164-old_a11_component_8-marking.json"
HANDOFF = GENERATED / "elkies-k3-h3-q4o1584-route-optimization-handoff.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--output",
    type=Path,
    default=LOCAL / "q4o164-integral-basis-height-gram-audit-qq.json",
)
args = parser.parse_args()
OUTPUT = args.output if args.output.is_absolute() else ROOT / args.output
INPUTS = (MODEL, BASIS, MARKING, HANDOFF)

started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


model = json.loads(MODEL.read_text())
basis = json.loads(BASIS.read_text())
marking = json.loads(MARKING.read_text())
handoff = json.loads(HANDOFF.read_text())
assert model["status"] == "PASS_EXACT_QQ_Q4O164_COMPACT_WEIERSTRASS_NORMALIZATION"
assert basis["status"] == "PASS_EXACT_QQ_Q4O164_PAIR_NODE_SECTION_SUBGROUP_RANK8"
assert marking["status"] == "PASS_EXACT_Q4O164_PHYSICAL_EFFECTIVE_ZERO_MARKING"
assert handoff["status"] == "PASS_EXACT_Q4O1584_Q4O164_Q8O376_Q12_ROOTLESS_OPTIONS_HANDOFF_PROMOTED"

RQ = PolynomialRing(QQ, "t")
KQ = RQ.fraction_field()
A = RQ([QQ(value) for value in model["compact_model"]["A_coefficients_low_to_high"]])
B = RQ([QQ(value) for value in model["compact_model"]["B_coefficients_low_to_high"]])
E = EllipticCurve(KQ, [0, 0, 0, KQ(A), KQ(B)])


def exact_point(record):
    return E(
        KQ(RQ([QQ(value) for value in record["x_coefficients_low_to_high"]])),
        KQ(RQ([QQ(value) for value in record["y_coefficients_low_to_high"]])),
    )


def fourfold_height(point):
    fourfold = 4 * point
    x_coordinate, y_coordinate = fourfold[0], fourfold[1]
    x_numerator = x_coordinate.numerator()
    x_denominator = x_coordinate.denominator()
    y_numerator = y_coordinate.numerator()
    y_denominator = y_coordinate.denominator()
    pole_degree = max(x_denominator.degree(), x_numerator.degree() - 4)
    assert pole_degree >= 0 and pole_degree % 2 == 0
    assert x_denominator.is_monic() and y_denominator.is_monic()
    assert x_denominator**3 == y_denominator**2
    expected_degrees = (
        4 + pole_degree,
        pole_degree,
        6 + 3 * pole_degree // 2,
        3 * pole_degree // 2,
    )
    actual_degrees = (
        x_numerator.degree(),
        x_denominator.degree(),
        y_numerator.degree(),
        y_denominator.degree(),
    )
    assert actual_degrees == expected_degrees
    return {
        "fourfold_compact_degrees_x_num_x_den_y_num_y_den": list(map(int, actual_degrees)),
        "fourfold_x_denominator_cubed_equals_y_denominator_squared": True,
        "cleared_pole_degree": int(pole_degree),
        "canonical_height": str(QQ(4 + pole_degree) / 16),
    }


points = [exact_point(record) for record in basis["resolved_hensel"]["sections"]]
diagonal_audits = [fourfold_height(point) for point in points]
diagonal = vector(QQ, [QQ(record["canonical_height"]) for record in diagonal_audits])

pair_audits = []
gram = matrix(QQ, 8, 8)
for i in range(8):
    gram[i, i] = diagonal[i]
    for j in range(i):
        audit = fourfold_height(points[i] + points[j])
        height_sum = QQ(audit["canonical_height"])
        pairing = (height_sum - diagonal[i] - diagonal[j]) / 2
        gram[i, j] = pairing
        gram[j, i] = pairing
        pair_audits.append({"indices": [j, i], "height_of_sum": str(height_sum), **audit})

scaled_gram = 4 * gram
assert scaled_gram.change_ring(ZZ).det() != 0
assert gram.det() == QQ(459) / 8

# The old profile for B7 treated a singular-node hit as I4 component 2.  Its
# literal group-law consequence N remains exact, but the fourfold audit gives
# height 13/4, so the coarse component label and the old height 3 are invalid.
N = 2 * points[0] + points[5] + points[7]
N_audit = fourfold_height(N)
assert QQ(N_audit["canonical_height"]) == QQ(13) / 4

# Enumerate isometric embeddings in the C8-pointed marked MW9 lattice used by
# the equation.  This frame is deliberately not the earlier pre-C8 frame.
frame_path = ROOT / marking["frame_output"]
frame_rows = [
    [ZZ(value) for value in line.split()]
    for line in frame_path.read_text().splitlines()
    if line and not line.startswith("#")
]
frame = matrix(ZZ, frame_rows)
assert frame.nrows() == frame.ncols() == 17
root_gram = frame[:8, :8]
coupling = frame[:8, 8:]
marked_height = frame[8:, 8:].change_ring(QQ) - (
    coupling.transpose() * root_gram.inverse() * coupling
)
assert marked_height.det() == QQ(237) / 16
marked_scaled = (4 * marked_height).change_ring(ZZ)
raw_columns = matrix(ZZ, pari(marked_scaled).qfminim(max(scaled_gram.diagonal()))[2]).columns()
shell = [vector(ZZ, column) for column in raw_columns] + [-vector(ZZ, column) for column in raw_columns]

candidates_by_row = [
    [values for values in shell if values * marked_scaled * values == scaled_gram[i, i]]
    for i in range(8)
]
embeddings = []


def extend_embedding(rows):
    i = len(rows)
    if i == 8:
        embeddings.append(matrix(ZZ, rows))
        return
    for candidate in candidates_by_row[i]:
        if all(candidate * marked_scaled * rows[j] == scaled_gram[i, j] for j in range(i)):
            extend_embedding(rows + [candidate])


extend_embedding([])
assert len(shell) == 136
assert len(embeddings) == 16

# Order the two A1 and two A3 factors as equation profiles do.  Chain
# orientations are allowed to reverse later, so only a deterministic endpoint
# is needed here.
unseen = set(range(8))
root_components = []
while unseen:
    first = min(unseen)
    unseen.remove(first)
    pending = [first]
    component = []
    while pending:
        index = pending.pop()
        component.append(index)
        for other in list(unseen):
            if root_gram[index, other]:
                unseen.remove(other)
                pending.append(other)
    root_components.append(tuple(sorted(component)))
root_components.sort(key=lambda component: (len(component), component))
assert tuple(map(len, root_components)) == (1, 1, 3, 3)


def chain_order(component):
    if len(component) == 1:
        return component
    endpoints = [
        index for index in component
        if sum(root_gram[index, other] != 0 for other in component if other != index) == 1
    ]
    order = [min(endpoints)]
    while len(order) < len(component):
        order.append(next(
            other for other in component
            if other not in order and root_gram[order[-1], other] != 0
        ))
    return tuple(order)


root_orders = [chain_order(component) for component in root_components]


def marked_profile(tail):
    z = coupling * vector(QQ, tail).column()
    return tuple(
        sum((position + 1) * z[index, 0] for position, index in enumerate(order))
        % (len(order) + 1)
        for order in root_orders
    )


stored_profiles = [
    tuple(map(ZZ, record["component_profile"]))
    for record in basis["resolved_hensel"]["sections"][:7]
]


def transform_profile(profile, symmetry):
    a, b, c, d = profile
    swap_i2, swap_i4, reverse_c, reverse_d = symmetry
    aa, bb = (b, a) if swap_i2 else (a, b)
    cc = (-c) % 4 if reverse_c else c
    dd = (-d) % 4 if reverse_d else d
    if swap_i4:
        cc, dd = dd, cc
    return aa, bb, cc, dd


component_compatible = []
q8_record = next(
    record for record in handoff["promoted_q323_free_pinned_route"]["compiler_frames"]
    if record["edge"] == "q8/o376"
)
q8_horizontal = vector(ZZ, q8_record["horizontal_section"])
inherited_p1 = vector(ZZ, marking["equation_explicit_curves_in_child"]["P1"])
ns = matrix(ZZ, [[0, 1], [1, 0]]).block_sum(-frame)
assert q8_horizontal * ns * q8_horizontal == -2
assert inherited_p1 * ns * inherited_p1 == -2
assert q8_horizontal[0] - q8_horizontal[1] == 4
assert inherited_p1[1] == 7
trace_tail = vector(ZZ, inherited_p1[-9:])
q8_tail = vector(ZZ, q8_horizontal[-9:])
named_tails = {
    index: vector(ZZ, marking["equation_explicit_curves_in_child"][f"old_A11_component_{index}"][-9:])
    for index in (1, 2, 3, 7)
}
named_correction = -8 * named_tails[1] + 5 * named_tails[2] - 2 * named_tails[3] - 7 * named_tails[7]
residual = q8_tail - trace_tail - named_correction
assert residual == vector(ZZ, [0, 0, 0, 0, -1, -1, 2, 0, 0])
assert q8_tail * marked_height * q8_tail == 11
assert trace_tail * marked_height * trace_tail == QQ(117) / 4
embedding_records = []
for index, embedding in enumerate(embeddings):
    profiles = [marked_profile(row) for row in embedding.rows()]
    compatible_symmetries = [
        symmetry
        for symmetry in product((False, True), repeat=4)
        if all(transform_profile(profiles[i], symmetry) == stored_profiles[i] for i in range(7))
    ]
    compatible = bool(compatible_symmetries)
    compatible_equation_profiles = [
        [transform_profile(profile, symmetry) for profile in profiles]
        for symmetry in compatible_symmetries
    ]
    if compatible:
        component_compatible.append(index)
    try:
        word = embedding.transpose().solve_right(residual)
        contains_residual = all(value.denominator() == 1 for value in word)
    except ValueError:
        word = None
        contains_residual = False
    if contains_residual:
        word = vector(ZZ, word)
    embedding_records.append({
        "embedding_index": index,
        "rows_B0_through_B7_in_marked_MW9": [list(map(int, row)) for row in embedding.rows()],
        "marked_component_profiles_B0_through_B7": [list(map(int, profile)) for profile in profiles],
        "compatible_with_first_seven_stored_profiles_up_to_fibre_symmetry": compatible,
        "compatible_fibre_symmetries_swapI2_swapI4_reverseI4a_reverseI4b": [
            list(map(bool, symmetry)) for symmetry in compatible_symmetries
        ],
        "compatible_equation_component_profiles_B0_through_B7": [
            [list(map(int, profile)) for profile in profile_list]
            for profile_list in compatible_equation_profiles
        ],
        "contains_q8_residual": contains_residual,
        "q8_residual_word_in_B_basis": None if word is None else list(map(int, word)),
    })

assert len(component_compatible) == 8
residual_embeddings = [
    record["embedding_index"]
    for record in embedding_records
    if record["compatible_with_first_seven_stored_profiles_up_to_fibre_symmetry"]
    and record["contains_q8_residual"]
]
assert len(residual_embeddings) == 8

payload = {
    "schema": "elkies-k3.q4o164-integral-basis-height-gram-audit-qq.v2",
    "status": "PASS_EXACT_QQ_Q4O164_FOURFOLD_HEIGHT_GRAM_AND_C8_MARKED_EMBEDDING_CENSUS",
    "basis_order": "resolved_hensel.sections B0,...,B7 in q4o164-integral-basis-qq.json",
    "fourfold_height_audit": {
        "component_group_exponent": 4,
        "formula": "height(P)=(4+max(deg den x(4P),deg num x(4P)-4))/16",
        "basis_sections": diagonal_audits,
        "pair_sums": pair_audits,
        "height_gram": [[str(value) for value in row] for row in gram.rows()],
        "four_times_height_gram": [list(map(int, row)) for row in scaled_gram.rows()],
        "height_gram_determinant": str(gram.det()),
    },
    "withdrawn_coarse_profile_consequence": {
        "literal_group_relation_retained": "N=2*B0+B5+B7",
        "old_B7_component_profile": [1, 0, 0, 2],
        "old_N_component_profile": [0, 0, 0, 2],
        "old_N_height": "3",
        "corrected_B7_equation_profile": [1, 0, 0, "1_or_3"],
        "corrected_N_equation_profile": [0, 0, 0, "1_or_3"],
        "corrected_N_fourfold_audit": N_audit,
    },
    "marked_embedding_enumeration": {
        "signed_vectors_in_shell": len(shell),
        "isometric_embeddings": len(embeddings),
        "component_compatible_embedding_indices": component_compatible,
        "component_compatible_embeddings_containing_q8_residual": residual_embeddings,
        "q8_residual_marked_MW9_tail": list(map(int, residual)),
        "q8_horizontal_marked_MW9_tail": list(map(int, q8_tail)),
        "inherited_P1_trace_marked_MW9_tail": list(map(int, trace_tail)),
        "named_trace_correction": "-8*C1+5*C2-2*C3-7*C7",
        "named_trace_correction_marked_MW9_tail": list(map(int, named_correction)),
        "q8_horizontal_height": "11",
        "inherited_P1_trace_height": "117/4",
        "residual_basis_words_are_embedding_dependent": True,
        "embeddings": embedding_records,
    },
    "method": {
        "large_Groebner_required": False,
        "exact_QQ_function_field_group_law": True,
        "finite_positive_definite_embedding_enumeration": True,
        "runtime_seconds": time.monotonic() - started,
    },
    "proof_boundary": (
        "Fourfold pole growth corrects the height Gram and invalidates the old coarse I4 "
        "component-2 label. In the correct C8-pointed frame, the first seven stored profiles "
        "leave eight marked embeddings and all eight contain the q8 residual, with different "
        "basis words. This exact QQ audit does not choose among them; a separate modular pole-"
        "degree comparison with the Abel trace supplies that construction-level selection."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS] + [str(frame_path.relative_to(ROOT))],
        "sha256": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in INPUTS + (frame_path,)
        },
    },
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4O164HEIGHTGRAM|det={}|embeddings={}|compatible={}|residual={}|status={}|output={}".format(
        gram.det(), len(embeddings), len(component_compatible), len(residual_embeddings),
        payload["status"], OUTPUT,
    ),
    flush=True,
)
