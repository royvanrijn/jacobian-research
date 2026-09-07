#!/usr/bin/env sage -python
"""Compile the exact q4 return by transposing the q4/o230 quartic.

The prescribed return fibre is exactly the original 2A5 fibre class, so its
H0 is the visible span <1,T>.  Viewing the existing bidegree-(4,4) equation as
a quartic in U gives the return genus-one curve without another RR solve.
Its two points at U=infinity have ordinates +/-Z(T): +Z is the original zero,
while -Z is the other component of the child I2 fibre and supplies the changed
zero.  No Groebner basis, polynomial gcd, or factorization is used.
"""

import hashlib
import json
from pathlib import Path

from sage.all import (
    IntegralLattice, PolynomialRing, QQ, ZZ, block_diagonal_matrix,
    identity_matrix, lcm, matrix, pari, vector, xgcd,
)


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
PARENT = LOCAL / "q24-a11-to-2a5-q8-resolved-rr-qq.json"
PARENT_MARKING = LOCAL / "q24-a11-to-2a5-q8-equation-marking-qq.json"
P230 = LOCAL / "q24-2a5-p230-scaled-x-qq.json"
Q4_WORD = LOCAL / "q24-2a5-q4o230-horizontal-word.json"
Q4_EQUATION = LOCAL / "q24-2a5-to-a1a4a5-q4o230-resolved-rr-qq.json"
Q4_MARKING = LOCAL / "q24-2a5-to-a1a4a5-q4o230-equation-marking-qq.json"
LOOPS = GENERATED / "elkies-k3-h3-a5a5-zero-changing-loop-search.json"
OUTPUT = LOCAL / "q24-a1a4a5-to-2a5-q4-return-resolved-rr-qq.json"
INPUTS = (PARENT, PARENT_MARKING, P230, Q4_WORD, Q4_EQUATION, Q4_MARKING, LOOPS)

parent = json.loads(PARENT.read_text())
parent_marking = json.loads(PARENT_MARKING.read_text())
p230 = json.loads(P230.read_text())
q4_word = json.loads(Q4_WORD.read_text())
q4_equation = json.loads(Q4_EQUATION.read_text())
q4_marking = json.loads(Q4_MARKING.read_text())
loops = json.loads(LOOPS.read_text())
assert parent["status"] == "PASS_EXACT_Q24_A11_Q8_2A5_RESOLVED_RR"
assert parent_marking["status"] == "PASS_EXACT_Q24_A11_Q8_2A5_EQUATION_MARKING"
assert p230["status"] == "PASS_EXACT_QQ_P230_SECTION_AND_H0"
assert q4_word["status"] == "PASS_EXACT_Q24_2A5_Q4O230_LOW_POLE_HORIZONTAL_WORD"
assert q4_equation["status"] == "PASS_EXACT_Q24_2A5_Q4O230_A1A4A5_RESOLVED_RR"
assert q4_marking["status"] == "PASS_EXACT_Q24_2A5_Q4O230_COMPONENT10_EQUATION_MARKING"
assert loops["status"] == "PASS_EXACT_ZERO_CHANGING_LOOP_SEARCH"

RB = PolynomialRing(QQ, "T")
T = RB.gen()
K = RB.fraction_field()
RV = PolynomialRing(K, "V")
V = RV.gen()


def poly_t(values):
    return RB([QQ(value) for value in values])


A_parent = poly_t(parent["child"]["minimal_A_coefficients_low_to_high"])
B_parent = poly_t(parent["child"]["minimal_B_coefficients_low_to_high"])
Z = poly_t(p230["P230"]["Z_coefficients_low_to_high"])

# Transpose the exact coefficient matrix q(T,U).  The return quartic variable
# V is the old new-base variable U; T is again the elliptic base.
by_t_degree = [
    [QQ(value) for value in coefficients]
    for coefficients in q4_equation["quartic"][
        "coefficients_by_old_degree_each_low_to_high_in_U"
    ]
]
assert len(by_t_degree) == 5 and all(len(row) <= 5 for row in by_t_degree)
by_u_degree = [
    RB([row[degree] if degree < len(row) else 0 for row in by_t_degree])
    for degree in range(5)
]
quartic_return = RV([K(value) for value in by_u_degree])
assert quartic_return.degree() == 4

e, d, c, b, a = quartic_return.list()
I = 12 * a * e - 3 * b * d + c**2
J = 72 * a * c * e + 9 * b * c * d - 27 * a * d**2 - 27 * b**2 * e - 2 * c**3
invariant_A = K(-27 * I)
invariant_B = K(-27 * J)
assert invariant_A == 1296 * K(A_parent)
assert invariant_B == 46656 * K(B_parent)

# At V=infinity use S=1/V and W'=W/V^2.  The reciprocal quartic has constant
# coefficient Z(T)^2.  Pointing it at -Z selects the changed zero; +Z is the
# original component-9 zero.
reciprocal = RV(list(reversed(quartic_return.list())))
q = reciprocal.list()
assert len(q) == 5 and q[0] == K(Z)**2
e0, d0, c0, b0, a0 = q
w_changed = -K(Z)
w_original = K(Z)
a1 = d0 / w_changed
a2 = c0 - d0**2 / (4 * w_changed**2)
a3 = 2 * w_changed * b0
a4 = -4 * w_changed**2 * a0
a6 = a2 * a4
b2 = a1**2 + 4 * a2
b4 = 2 * a4 + a1 * a3
b6 = a3**2 + 4 * a6
c4 = b2**2 - 24 * b4
c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
pointed_A = -c4 / 48
pointed_B = -c6 / 864
assert pointed_A == 16 * K(A_parent) and pointed_B == 64 * K(B_parent)

# The opposite infinity sign is the original component-9 zero, now a section
# on the changed-zero return model.
x_general = -a2
y_general = a1 * a2 - a3
x_original = K((x_general + b2 / 12) / 4)
y_original = K((y_general + (a1 * x_general + a3) / 2) / 8)
assert y_original**2 == x_original**3 + K(A_parent) * x_original + K(B_parent)

# -------------------------------------------------------------------------
# Equation-adapted returned marking and its exact isometry to the stored route
# chamber.  This repeats only the deterministic rank-17 lattice adaptation.
# -------------------------------------------------------------------------
record = next(
    row for row in loops["ranked_loops"]
    if row["first_edge_candidate_id"] == {"q": 4, "old_fibre_degree": 2, "orbit_index": 230}
    and row["explicit_zero_curve"] == "old_A11_component_10"
)
child_to_parent = matrix(ZZ, q4_marking["transport"][
    "component10_zero_child_to_component9_zero_parent_basis"
])
parent_to_child = matrix(ZZ, q4_marking["transport"][
    "component9_zero_parent_to_component10_zero_child_basis"
])
parent_frame_path = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frame.txt"
parent_frame = matrix(ZZ, [
    [ZZ(value) for value in line.split()]
    for line in parent_frame_path.read_text().splitlines()
    if line.strip() and not line.lstrip().startswith("#")
])
U2 = matrix(ZZ, ((0, 1), (1, 0)))
g_parent = block_diagonal_matrix(U2, -parent_frame)
g_explicit = child_to_parent * g_parent * child_to_parent.transpose()
return_fibre = vector(ZZ, record["return_fibre_in_explicit_child"])
old_fibre = vector(ZZ, [1, 0] + [0] * 17)
assert old_fibre * parent_to_child == return_fibre


def roots_and_data(gram):
    result = pari(gram).qfminim(2)
    count = ZZ(result[0])
    half = [vector(ZZ, column) for column in matrix(ZZ, result[2]).columns()]
    roots = tuple(half + [-root for root in half])
    basis = matrix(ZZ, [list(root) for root in roots]).row_module().basis_matrix()
    root_gram = basis * gram * basis.transpose()
    return roots, basis, (basis.rank(), count, abs(ZZ(root_gram.det())))


def deterministic_simple_roots(gram):
    roots, unused, data = roots_and_data(gram)
    positive = [root for root in roots if next(value for value in root if value) > 0]
    positive_set = {tuple(root) for root in positive}
    simple = [
        root for root in positive
        if not any(tuple(root-left) in positive_set for left in positive)
    ]
    result = matrix(ZZ, [list(root) for root in simple])
    assert result.nrows() == result.rank() == data[0]
    return result


def root_adaptation(child):
    unused, root_basis, data = roots_and_data(child)
    rank = data[0]
    smith, unused_left, right = root_basis.smith_form()
    assert tuple(abs(smith[index, index]) for index in range(rank)) == (1,) * rank
    simple = deterministic_simple_roots(child)
    completion = right.inverse()
    basis = simple.stack(completion[rank:])
    adapted = basis * child * basis.transpose()
    cartan = adapted[:rank, :rank]
    coupling = adapted[:rank, rank:]
    height = adapted[rank:, rank:] - coupling.transpose() * cartan.inverse() * coupling
    scale = lcm(value.denominator() for value in height.list())
    lll = matrix(ZZ, pari((scale * height).change_ring(ZZ)).qflllgram())
    basis = block_diagonal_matrix(identity_matrix(ZZ, rank), lll.transpose()) * basis
    adapted = basis * child * basis.transpose()
    assert abs(basis.det()) == 1
    return adapted, basis, data


def bezout_vector_for_pairing(ns, fibre):
    current = ZZ(0)
    answer = [ZZ(0)] * ns.nrows()
    for index, value in enumerate(ns * fibre):
        if not value:
            continue
        divisor, left, right = xgcd(current, ZZ(value))
        answer = [left * entry for entry in answer]
        answer[index] += right
        current = divisor
    assert abs(current) == 1
    if current == -1:
        answer = [-entry for entry in answer]
    return vector(ZZ, answer)


mate = bezout_vector_for_pairing(g_explicit, return_fibre)
mate -= ZZ(mate * g_explicit * mate // 2) * return_fibre
kernel = matrix(ZZ, [
    list(return_fibre * g_explicit), list(mate * g_explicit)
]).right_kernel_matrix()
returned_raw = -(kernel * g_explicit * kernel.transpose())
raw_transition = matrix(ZZ, [list(return_fibre), list(mate)] + list(kernel.rows()))
returned_frame, adaptation, root_data = root_adaptation(returned_raw)
assert tuple(root_data) == (10, 60, 36)
assert returned_frame == matrix(ZZ, record["returned_A5A5_frame"])
initial_transition = block_diagonal_matrix(identity_matrix(ZZ, 2), adaptation) * raw_transition
initial_inverse = initial_transition.inverse().change_ring(ZZ)
g_returned = block_diagonal_matrix(U2, -returned_frame)
assert initial_transition * g_explicit * initial_transition.transpose() == g_returned

changed_zero_explicit = vector(ZZ, initial_transition.row(1) - initial_transition.row(0))
i2_components = [
    vector(ZZ, values)
    for values in q4_marking["physical_fibres"]["I2_infinity"]["components_in_cycle_order"]
]
assert changed_zero_explicit == i2_components[1]
assert changed_zero_explicit * initial_inverse == vector(ZZ, [-1, 1] + [0] * 17)
original_zero_parent = vector(ZZ, [-1, 1] + [0] * 17)
original_zero_returned = original_zero_parent * parent_to_child * initial_inverse
assert original_zero_returned * g_returned * vector(ZZ, [1, 0] + [0] * 17) == 1

stored_transition = matrix(ZZ, record["return_transition"])
stored_to_equation_adapted = stored_transition * initial_inverse
equation_adapted_to_stored = stored_to_equation_adapted.inverse().change_ring(ZZ)
assert stored_to_equation_adapted * g_returned * stored_to_equation_adapted.transpose() == g_returned
assert abs(stored_to_equation_adapted.det()) == 1


def rational_record(value):
    value = K(value)
    numerator = RB(value.numerator())
    denominator = RB(value.denominator())
    if denominator.leading_coefficient() < 0:
        numerator, denominator = -numerator, -denominator
    return {
        "numerator_coefficients_low_to_high": [str(item) for item in numerator.list()],
        "denominator_coefficients_low_to_high": [str(item) for item in denominator.list()],
        "degrees_numerator_denominator": [int(numerator.degree()), int(denominator.degree())],
    }


def rows(value):
    return [[int(item) for item in row] for row in value.rows()]


def height_profile(polynomials):
    values = [value for polynomial in polynomials for value in polynomial.list()]
    return {
        "maximum_numerator_bits": int(max(abs(value.numerator()).nbits() for value in values)),
        "maximum_denominator_bits": int(max(value.denominator().nbits() for value in values)),
        "maximum_rational_bits": int(max(
            max(abs(value.numerator()).nbits(), value.denominator().nbits()) for value in values
        )),
    }


payload = {
    "schema": "elkies-k3.q24-a1a4a5-q4-return-2a5-resolved-rr-qq.v1",
    "status": "PASS_EXACT_Q24_A1A4A5_Q4_RETURN_2A5_CHANGED_ZERO",
    "edge": {
        "source": "A1+A4+A5/MW7 component10 zero",
        "q": 4,
        "target": "2A5/MW7 changed zero",
        "return_fibre_is_exact_original_2A5_fibre": True,
    },
    "resolved_RR": {
        "divisor": "original 2A5 fibre class",
        "h0": 2,
        "basis": ["1", "T"],
        "proof": "visible primitive base pencil on the exact bidegree-(4,4) model",
        "new_ambient_solve_required": False,
    },
    "quartic": {
        "base_variable": "T",
        "quartic_variable": "V (the previous U)",
        "coefficients_by_V_degree_each_low_to_high_in_T": [
            [str(item) for item in value.list()] for value in by_u_degree
        ],
        "bidegree": [4, 4],
        "exact_transpose_of_forward_quartic": True,
        "height_profile": height_profile(by_u_degree),
    },
    "pointed_infinity": {
        "reciprocal_coordinate": "S=1/V, W'=W/V^2",
        "original_zero_ordinate": rational_record(w_original),
        "changed_zero_ordinate": rational_record(w_changed),
        "changed_zero_curve": "nonidentity component of the forward child I2 at infinity",
        "selected_origin": "changed zero (-Z sign)",
        "exact_square_leading_coefficient": True,
    },
    "original_zero_on_changed_zero_model": {
        "x": rational_record(x_original),
        "y": rational_record(y_original),
        "exact_child_identity": True,
        "NS_coordinates": [int(item) for item in original_zero_returned],
    },
    "child": {
        "minimal_A_coefficients_low_to_high": [str(item) for item in A_parent.list()],
        "minimal_B_coefficients_low_to_high": [str(item) for item in B_parent.list()],
        "invariant_quartic_scaling": "A_invariant=1296*A_minimal, B_invariant=46656*B_minimal",
        "pointed_minimal_scaling": "x_pointed=4*x_minimal, y_pointed=8*y_minimal",
        "fibre_profile": "2I6+12I1",
        "root_lattice": "2A5",
        "root_rank": 10,
        "root_determinant": 36,
        "MW_rank_if_rho19": 7,
        "height_profile": height_profile((A_parent, B_parent)),
    },
    "transport": {
        "equation_adapted_return_to_component10_child_basis": rows(initial_transition),
        "component10_child_to_equation_adapted_return_basis": rows(initial_inverse),
        "stored_promoted_return_to_equation_adapted_return_basis": rows(stored_to_equation_adapted),
        "equation_adapted_return_to_stored_promoted_return_basis": rows(equation_adapted_to_stored),
        "equation_adapted_zero": "forward I2 nonidentity component",
        "stored_route_chamber_isometry_exact": True,
        "determinants": [
            int(initial_transition.det()), int(initial_inverse.det()),
            int(stored_to_equation_adapted.det()), int(equation_adapted_to_stored.det()),
        ],
        "Gram_transport_exact": True,
    },
    "verification": {
        "exact_H0_basis": True,
        "exact_quartic_transpose": True,
        "exact_binary_quartic_invariants": True,
        "exact_pointed_minimal_Jacobian": True,
        "exact_changed_zero_identification": True,
        "exact_route_chamber_isometry": True,
    },
    "large_Groebner_required": False,
    "next_required": "transport the promoted q6/orbit1315 fibre into this equation-adapted return marking",
    "proof_boundary": (
        "This proves the exact q4 return equation, H0, changed zero, minimized 2A5 Jacobian, "
        "and an exact integral isometry to the stored promoted-route chamber. The next q6 "
        "equation lift has not yet been compiled."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4RETURNQQ|H0=2|quartic=4|child=2I6+12I1|zero=opposite_I2_infinity|"
    "route_isometry_det={}|status={}|output={}".format(
        stored_to_equation_adapted.det(), payload["status"], OUTPUT
    ),
    flush=True,
)
