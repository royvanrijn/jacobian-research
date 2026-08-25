#!/usr/bin/env sage -python
"""Point the exact q4/o230 child at old_A11_component_10 over QQ.

At the parent I6 fibre where P230 meets old_A11_component_5, its Weierstrass
specialization is the node.  A variable line through that node has one nodal
residual point (component 5) and one smooth residual point (component 10).
Consequently the exact quartic slice is the square

    w_10^2,  w_10 = (N^2 - 3*X*Db^2)/Z^3.

Pointing the translated quartic at the smooth sign identifies component 10
as the child zero.  The physical I2/I5/I6 cycles and the bidirectional NS
transport are then replayed exactly.  No Groebner basis or factorization is
used.
"""

import hashlib
import json
import math
from pathlib import Path

from sage.all import (
    GF, PolynomialRing, QQ, ZZ, block_diagonal_matrix, identity_matrix,
    matrix, vector,
)


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
PARENT = LOCAL / "q24-a11-to-2a5-q8-resolved-rr-qq.json"
PARENT_MARKING = LOCAL / "q24-a11-to-2a5-q8-equation-marking-qq.json"
P230 = LOCAL / "q24-2a5-p230-scaled-x-qq.json"
Q4_WORD = LOCAL / "q24-2a5-q4o230-horizontal-word.json"
Q4_EQUATION = LOCAL / "q24-2a5-to-a1a4a5-q4o230-resolved-rr-qq.json"
SUPPORTS = LOCAL / "q24-2a5-q4o230-repeated-supports.json"
LOOPS = GENERATED / "elkies-k3-h3-a5a5-zero-changing-loop-search.json"
PARENT_FRAME = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frame.txt"
OUTPUT = LOCAL / "q24-2a5-to-a1a4a5-q4o230-equation-marking-qq.json"
INPUTS = (PARENT, PARENT_MARKING, P230, Q4_WORD, Q4_EQUATION, SUPPORTS, LOOPS, PARENT_FRAME)

parent = json.loads(PARENT.read_text())
parent_marking = json.loads(PARENT_MARKING.read_text())
p230 = json.loads(P230.read_text())
q4_word = json.loads(Q4_WORD.read_text())
q4_equation = json.loads(Q4_EQUATION.read_text())
supports = json.loads(SUPPORTS.read_text())
loops = json.loads(LOOPS.read_text())
assert parent["status"] == "PASS_EXACT_Q24_A11_Q8_2A5_RESOLVED_RR"
assert parent_marking["status"] == "PASS_EXACT_Q24_A11_Q8_2A5_EQUATION_MARKING"
assert p230["status"] == "PASS_EXACT_QQ_P230_SECTION_AND_H0"
assert q4_word["status"] == "PASS_EXACT_Q24_2A5_Q4O230_LOW_POLE_HORIZONTAL_WORD"
assert q4_equation["status"] == "PASS_EXACT_Q24_2A5_Q4O230_A1A4A5_RESOLVED_RR"
assert supports["status"] == "PASS_EXACT_Q24_2A5_Q4O230_REPEATED_SUPPORTS"
assert loops["status"] == "PASS_EXACT_ZERO_CHANGING_LOOP_SEARCH"

RT0 = PolynomialRing(QQ, "T")
T0 = RT0.gen()
RU = PolynomialRing(QQ, "U")
U = RU.gen()
KU = RU.fraction_field()
RT = PolynomialRing(KU, "T")
T = RT.gen()


def poly_t(values):
    return RT0([QQ(value) for value in values])


def poly_u(values):
    return RU([QQ(value) for value in values])


A_parent = poly_t(parent["child"]["minimal_A_coefficients_low_to_high"])
B_parent = poly_t(parent["child"]["minimal_B_coefficients_low_to_high"])
section = p230["P230"]
X = poly_t(section["X_coefficients_low_to_high"])
Y = poly_t(section["Y_coefficients_low_to_high"])
Z = poly_t(section["Z_coefficients_low_to_high"])
basis = p230["H0"]["constant_function_basis"]
AA0 = poly_t(basis[0]["AA_coefficients_low_to_high"])
BB0 = QQ(basis[0]["BB"])
AA1 = poly_t(basis[1]["AA_coefficients_low_to_high"])
BB1 = QQ(basis[1]["BB"])
assert AA0 == Z**2 and not BB0 and BB1 == 1

# Select the physical parent I6 at T=beta whose reduction is 89.  The pinned
# Tate orientation labels its identity component as old_A11_component_10 and
# the component met by P230 as old_A11_component_5.
Fp = GF(103)
i6_roots = []
for row in parent["child"]["discriminant_factorization"]:
    if int(row["multiplicity"]) != 6:
        continue
    factor = RT0(row["factor"])
    assert factor.degree() == 1
    root = QQ(-factor[0] / factor[1])
    i6_roots.append((int(Fp(root)), root))
assert sorted(value for value, unused in i6_roots) == [68, 89]
beta = next(root for value, root in i6_roots if value == 89)
x_beta, y_beta, z_beta = X(beta), Y(beta), Z(beta)
assert z_beta and not y_beta
assert A_parent(beta) * z_beta**4 == -3 * x_beta**2
assert B_parent(beta) * z_beta**6 == 2 * x_beta**3

N = KU(AA1(beta)) - KU(U) * KU(AA0(beta))
Db = KU(U * BB0 - BB1)
w10 = KU((N**2 - 3 * KU(x_beta) * Db**2) / KU(z_beta**3))
w5 = -w10

quartic = RT([
    KU(poly_u(values))
    for values in q4_equation["quartic"][
        "coefficients_by_old_degree_each_low_to_high_in_U"
    ]
])
assert quartic.degree() == 4
q = [KU(quartic.derivative(order)(KU(beta)) / math.factorial(order))
     for order in range(5)]
assert q[0] == w10**2 == w5**2

# The I5 support is exactly the tangent value on this nodal parent fibre.
support5 = QQ(supports["supports"]["I5"]["U"])
support6 = QQ(supports["supports"]["I6"]["U"])
assert w10(support5) == 0 and support5 != support6

# Point the translated quartic at (T-beta,w)=(0,w10).  The invariant 9/27
# scaling returns exactly the short Jacobian already compiled from I,J.
e, d, c, b, a = q
w0 = w10
a1 = d / w0
a2 = c - d**2 / (4 * w0**2)
a3 = 2 * w0 * b
a4 = -4 * w0**2 * a
a6 = a2 * a4
b2 = a1**2 + 4 * a2
b4 = 2 * a4 + a1 * a3
b6 = a3**2 + 4 * a6
c4 = b2**2 - 24 * b4
c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
pointed_A = -c4 / 48
pointed_B = -c6 / 864
A_child = KU(poly_u(q4_equation["child"]["minimal_A_coefficients_low_to_high"]))
B_child = KU(poly_u(q4_equation["child"]["minimal_B_coefficients_low_to_high"]))
assert 81 * pointed_A == A_child and 729 * pointed_B == B_child

# The opposite sign is old component 5, now an exact section with component
# 10 as origin.
x_general = -a2
y_general = a1 * a2 - a3
x5 = KU(9 * (x_general + b2 / 12))
y5 = KU(27 * (y_general + (a1 * x_general + a3) / 2))
assert y5**2 == x5**3 + A_child * x5 + B_child

# -------------------------------------------------------------------------
# Exact physical cycles and bidirectional NS transport.
# -------------------------------------------------------------------------
frame = matrix(ZZ, [
    [ZZ(value) for value in line.split()]
    for line in PARENT_FRAME.read_text().splitlines()
    if line.strip() and not line.lstrip().startswith("#")
])
g_parent = block_diagonal_matrix(matrix(ZZ, ((0, 1), (1, 0))), -frame)
parent_forward = matrix(ZZ, parent_marking["transport"][
    "equation_A11_to_component9_zero_2A5_basis"
])
parent_inverse = matrix(ZZ, parent_marking["transport"][
    "component9_zero_2A5_to_equation_A11_basis"
])
record = next(
    row for row in loops["ranked_loops"]
    if row["first_edge_candidate_id"] == {"q": 4, "old_fibre_degree": 2, "orbit_index": 230}
    and row["explicit_zero_curve"] == "old_A11_component_10"
)
direct_transition = matrix(ZZ, record["equation_A11_to_explicit_child_basis"])
child_to_parent = direct_transition * parent_inverse
parent_to_child = child_to_parent.inverse().change_ring(ZZ)
assert child_to_parent.base_ring() == ZZ and parent_to_child.base_ring() == ZZ
assert abs(child_to_parent.det()) == abs(parent_to_child.det()) == 1
assert child_to_parent * parent_to_child == identity_matrix(ZZ, 19)
g_child = child_to_parent * g_parent * child_to_parent.transpose()
assert g_child[:2, :2] == matrix(ZZ, ((0, 1), (1, 0)))
assert not g_child[:2, 2:] and not g_child[2:, :2]

old_fibre = vector(ZZ, [1, 0] + [0] * 17)
old_zero = vector(ZZ, [-1, 1] + [0] * 17)
new_fibre_parent = vector(ZZ, q4_word["q4_divisor_correction"]["fibre_NS_coordinates"])
old_components = {
    index: vector(ZZ, parent_marking["physical_2A5"]["child_coordinates"][
        f"old_A11_component_{index}"
    ])
    for index in range(11)
}
chains = parent_marking["physical_2A5"]["chains"]
old_affine = {
    index: old_fibre - sum((old_components[item] for item in chain), vector(ZZ, 19))
    for index, chain in enumerate(chains)
}
assert vector(ZZ, child_to_parent.row(0)) == new_fibre_parent
assert vector(ZZ, child_to_parent.row(1) - child_to_parent.row(0)) == old_components[10]
assert old_components[10] * parent_to_child == vector(ZZ, [-1, 1] + [0] * 17)


def child(curve):
    return vector(ZZ, curve * parent_to_child)


Fnew = vector(ZZ, [1, 0] + [0] * 17)
zero10 = child(old_components[10])
known_cycles = {
    "I2_infinity": [child(old_components[9])],
    "I5": [child(old_affine[0]), child(old_components[0]),
           child(old_components[3]), child(old_components[4])],
    "I6": [child(old_affine[1]), child(old_components[1]),
           child(old_components[2]), child(old_components[6]),
           child(old_components[7])],
}
cycles = {
    name: values + [Fnew - sum(values, vector(ZZ, 19))]
    for name, values in known_cycles.items()
}


def verify_cycle(values):
    size = len(values)
    assert sum(values, vector(ZZ, 19)) == Fnew
    assert all(value * g_child * value == -2 for value in values)
    for left in range(size):
        for right in range(left + 1, size):
            expected = 2 if size == 2 else int((right-left) in (1, size-1))
            assert values[left] * g_child * values[right] == expected


for values in cycles.values():
    verify_cycle(values)
identity_indices = {}
for name, values in cycles.items():
    hits = [index for index, value in enumerate(values) if zero10 * g_child * value == 1]
    assert len(hits) == 1
    identity_indices[name] = hits[0]
assert identity_indices == {"I2_infinity": 0, "I5": 0, "I6": 5}

root_components = (
    [cycles["I2_infinity"][1]]
    + cycles["I5"][1:]
    + cycles["I6"][:5]
)
root_matrix = matrix(ZZ, root_components)
root_cartan = -(root_matrix * g_child * root_matrix.transpose())
A1 = matrix(ZZ, [[2]])
A4 = matrix(ZZ, 4, 4, lambda i, j: 2 if i == j else -1 if abs(i-j) == 1 else 0)
A5 = matrix(ZZ, 5, 5, lambda i, j: 2 if i == j else -1 if abs(i-j) == 1 else 0)
assert root_cartan == block_diagonal_matrix(A1, A4, A5)
assert root_cartan.det() == 60

# Component 5 is the opposite quartic sign and the second degree-one old
# curve.  Its equation and NS class therefore use the same physical label.
component5_child = child(old_components[5])
assert component5_child * g_child * Fnew == 1
assert component5_child * g_child * component5_child == -2


def rational_record(value):
    value = KU(value)
    numerator = RU(value.numerator())
    denominator = RU(value.denominator())
    if denominator.leading_coefficient() < 0:
        numerator, denominator = -numerator, -denominator
    return {
        "numerator_coefficients_low_to_high": [str(item) for item in numerator.list()],
        "denominator_coefficients_low_to_high": [str(item) for item in denominator.list()],
        "degrees_numerator_denominator": [int(numerator.degree()), int(denominator.degree())],
    }


def rows(value):
    return [[int(item) for item in row] for row in value.rows()]


def vectors(values):
    return [[int(item) for item in value] for value in values]


payload = {
    "schema": "elkies-k3.q24-2a5-q4o230-equation-marking-qq.v1",
    "status": "PASS_EXACT_Q24_2A5_Q4O230_COMPONENT10_EQUATION_MARKING",
    "selected_zero": "old_A11_component_10",
    "nodal_parent_slice": {
        "parent_I6_base_value": str(beta),
        "base_value_mod_103": 89,
        "P230_specialization": "node on old_A11_component_5",
        "node_projective_X_Y_Z": [str(x_beta), str(y_beta), str(z_beta)],
        "old_A11_component_10_ordinate": rational_record(w10),
        "old_A11_component_5_ordinate": rational_record(w5),
        "exact_quartic_square_identity": True,
        "I5_tangent_support": str(support5),
    },
    "pointed_quartic": {
        "translated_old_base_coordinate": "T-beta",
        "selected_origin": "old_A11_component_10",
        "selected_origin_maps_to": "point at infinity on the pointed generalized Weierstrass model",
        "invariant_short_scaling": "x_child=9*x_short, y_child=27*y_short",
        "invariant_jacobian_identity": True,
    },
    "old_A11_component_5_on_component10_pointed_child": {
        "x": rational_record(x5),
        "y": rational_record(y5),
        "exact_child_identity": True,
        "NS_coordinates": [int(item) for item in component5_child],
    },
    "physical_fibres": {
        "I2_infinity": {
            "support": "infinity",
            "components_in_cycle_order": vectors(cycles["I2_infinity"]),
            "inherited_labels_before_missing_component": ["old_A11_component_9"],
            "identity_component_index": identity_indices["I2_infinity"],
        },
        "I5": {
            "support": str(support5),
            "components_in_cycle_order": vectors(cycles["I5"]),
            "inherited_labels_before_missing_component": [
                "first_I6_affine_component", "old_A11_component_0",
                "old_A11_component_3", "old_A11_component_4",
            ],
            "identity_component_index": identity_indices["I5"],
        },
        "I6": {
            "support": str(support6),
            "components_in_cycle_order": vectors(cycles["I6"]),
            "inherited_labels_before_missing_component": [
                "second_I6_affine_component", "old_A11_component_1",
                "old_A11_component_2", "old_A11_component_6",
                "old_A11_component_7",
            ],
            "identity_component_index": identity_indices["I6"],
        },
        "root_cartan": rows(root_cartan),
        "root_type": "A1+A4+A5",
        "root_rank": 10,
        "root_determinant": 60,
    },
    "transport": {
        "component10_zero_child_to_component9_zero_parent_basis": rows(child_to_parent),
        "component9_zero_parent_to_component10_zero_child_basis": rows(parent_to_child),
        "equation_A11_to_component10_zero_child_basis": rows(direct_transition),
        "forward_determinant": int(child_to_parent.det()),
        "inverse_determinant": int(parent_to_child.det()),
        "inverse_exact": True,
        "Gram_transport_exact": True,
    },
    "large_Groebner_required": False,
    "next_required": "compile the exact q4 return from this component10-zero child",
    "proof_boundary": (
        "This identifies old_A11_component_10 as an exact rational point of the q4 quartic, "
        "points the invariant Jacobian at that curve, attaches the physical I2/I5/I6 cycles, "
        "and verifies the full bidirectional unimodular NS transport. The q4 return equations "
        "are the next gate."
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
    "Q4O230MARKQQ|zero=old_A11_component_10|other=old_A11_component_5|"
    "fibres=I2+I5+I6|root=A1+A4+A5|det={}/{}|status={}|output={}".format(
        child_to_parent.det(), parent_to_child.det(), payload["status"], OUTPUT
    ),
    flush=True,
)
