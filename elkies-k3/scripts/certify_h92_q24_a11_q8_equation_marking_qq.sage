#!/usr/bin/env sage -python
"""Attach the exact orbit12 zero and physical 2A5 marking to the q8 child.

The exact q8 RR certificate produces a binary quartic but deliberately leaves
its two rational points over the old I12 fibre unordered.  The q8 horizontal
section specializes to the node of that fibre.  A line of slope ``m`` through
the node has the node and the affine point as its two residual intersections,
so the affine discriminant sign is ``m^2-3*x_node``.  In the fraction-free RR
normalization this gives

    w_affine = (N^2 - 3*X*Db^2)/Z^3,

where ``N=AA1-T*AA0`` and ``Db=T*BB0-BB1``.  The complementary sign is the
only other old-I12 curve of q8 degree one, namely old_A11_component_9, which
is the compiler-selected zero.

The script checks this identity exactly by a four-jet of the already-certified
quartic numerator, constructs the pointed quartic with component 9 at
infinity, and verifies that its short model is the invariant Jacobian used by
the equation lift.  It then attaches the ten degree-zero old-I12 components
as the two physical A5 chains and replays the full bidirectional unimodular NS
transport.  No Groebner basis, factorization, or nonlinear solve is used.
"""

import hashlib
import json
import math
from pathlib import Path

from sage.all import (
    PolynomialRing,
    PowerSeriesRing,
    QQ,
    ZZ,
    block_diagonal_matrix,
    identity_matrix,
    matrix,
    vector,
)


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
MODEL = LOCAL / "q24-d12-to-a11-orbit42-resolved-rr-qq.json"
DIFFERENCE = LOCAL / "q24-a11-q8-difference-section-qq.json"
RR = LOCAL / "q24-a11-to-2a5-q8-resolved-rr-qq.json"
LATTICE = GENERATED / "elkies-k3-h3-a11-q8-orbit12-lattice-certificate.json"
ZERO_FRAMES = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frames.json"
NEIGHBORS = LOCAL / "q24-a11-orbit64-q8-all.json"
OUTPUT = LOCAL / "q24-a11-to-2a5-q8-equation-marking-qq.json"
INPUTS = (MODEL, DIFFERENCE, RR, LATTICE, ZERO_FRAMES, NEIGHBORS)

for path in INPUTS:
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

model = json.loads(MODEL.read_text())
difference = json.loads(DIFFERENCE.read_text())
rr = json.loads(RR.read_text())
lattice = json.loads(LATTICE.read_text())
zero_frames = json.loads(ZERO_FRAMES.read_text())
neighbors = json.loads(NEIGHBORS.read_text())

assert model["status"] == "PASS_EXACT_Q24_D12_Q6_A11_COMPONENT_VALUATION_RR"
assert difference["status"] == "PASS_EXACT_Q24_A11_Q8_DIFFERENCE_SECTION_QQ"
assert rr["status"] == "PASS_EXACT_Q24_A11_Q8_2A5_RESOLVED_RR"
assert lattice["status"] == "PASS_EXACT_A11_Q8_EQUATION_COST_LATTICE_CERTIFICATE"
assert zero_frames["status"] == "PASS_EXACT_A11_Q8_ORBIT12_EXPLICIT_ZERO_FRAMES"
assert zero_frames["selected_zero_curve"] == "old_A11_component_9"

VQ = PolynomialRing(QQ, "V")
V = VQ.gen()
TQ = PolynomialRing(QQ, "T")
T = TQ.gen()
KT = TQ.fraction_field()


def polynomial(values):
    return VQ([QQ(value) for value in values])


section = difference["section"]
X = polynomial(section["X_coefficients_low_to_high"])
Y = polynomial(section["Y_coefficients_low_to_high"])
Z = polynomial(section["Z_coefficients_low_to_high"])
A = polynomial(model["child"]["minimal_A_coefficients_low_to_high"])
B = polynomial(model["child"]["minimal_B_coefficients_low_to_high"])

i12_factor_text = next(
    row["factor"]
    for row in model["child"]["discriminant_factorization"]
    if int(row["multiplicity"]) == 12
)
i12_factor = TQ(i12_factor_text)
if i12_factor.degree() != 1:
    raise ArithmeticError("old A11 I12 factor is not linear")
beta = QQ(-i12_factor[0] / i12_factor[1])
z_beta = QQ(Z(beta))
x_beta = QQ(X(beta))
y_beta = QQ(Y(beta))
if not z_beta or y_beta:
    raise ArithmeticError("q8 horizontal does not specialize to a finite I12 node")
if A(beta) * z_beta**4 != -3 * x_beta**2:
    raise ArithmeticError("I12 specialization misses the nodal derivative equation")
if B(beta) * z_beta**6 != 2 * x_beta**3:
    raise ArithmeticError("I12 specialization misses the nodal cubic equation")

pairs = []
for row in rr["resolved_RR"]["basis_pairs"]:
    pairs.append((
        polynomial(row["AA_coefficients_low_to_high"]),
        polynomial(row["BB_coefficients_low_to_high"]),
    ))
(AA0, BB0), (AA1, BB1) = pairs
N = KT(AA1(beta)) - KT(T) * KT(AA0(beta))
Db = KT(T) * KT(BB0(beta)) - KT(BB1(beta))
if not Db:
    raise ArithmeticError("q8 RR basis degenerates at the old I12 fibre")
w_affine = KT((N**2 - 3 * KT(x_beta) * Db**2) / KT(z_beta**3))
w_component9 = -w_affine


def taylor(poly, precision=5):
    """Return the exact V=beta+u jet in KT[[u]]."""

    return PS([
        KT(poly.derivative(order)(beta) / math.factorial(order))
        for order in range(precision)
    ])


# Recover only the four-jet of the exact quartic.  This is much smaller than
# expanding its multi-megabyte global coefficients again.
PS = PowerSeriesRing(KT, "u", default_prec=5)
Ns = taylor(AA1) - KT(T) * taylor(AA0)
Dbs = KT(T) * taylor(BB0) - taylor(BB1)
Xs, Ys, Zs, As = map(taylor, (X, Y, Z, A))
raw = (
    Ns**4
    - 6 * Xs * Ns**2 * Dbs**2
    - 8 * Ys * Ns * Dbs**3
    - 3 * Xs**2 * Dbs**4
    - 4 * As * Zs**4 * Dbs**4
)
quartic_jet = raw / Zs**6
q = [KT(quartic_jet[index]) for index in range(5)]
if q[0] != w_affine**2 or q[0] != w_component9**2:
    raise ArithmeticError("nodal-cubic sign does not lie on the exact q8 quartic")

# Point the quartic at (V,w)=(beta,w_component9).  These formulas send that
# point to infinity on a generalized Weierstrass equation.  The 9/27 scaling
# gives the exact invariant short Jacobian A=-27I, B=-27J.
e, d, c, b, a = q
w0 = w_component9
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
pointed_short_A = -c4 / 48
pointed_short_B = -c6 / 864

I = 12 * a * e - 3 * b * d + c**2
J = 72 * a * c * e + 9 * b * c * d - 27 * a * d**2 - 27 * b**2 * e - 2 * c**3
if 81 * pointed_short_A != -27 * I or 729 * pointed_short_B != -27 * J:
    raise ArithmeticError("component-9 pointed quartic misses the invariant Jacobian")

# The opposite ordinate is the old affine component.  Record its exact point
# on the same short model, with component 9 as the origin.
x_general = -a2
y_general = a1 * a2 - a3
x_affine = KT(9 * (x_general + b2 / 12))
y_affine = KT(27 * (y_general + (a1 * x_general + a3) / 2))
A_child = KT(-27 * I)
B_child = KT(-27 * J)
if y_affine**2 != x_affine**3 + A_child * x_affine + B_child:
    raise ArithmeticError("old affine component misses the pointed 2A5 Jacobian")

# -------------------------------------------------------------------------
# Physical 2A5 chains and the exact bidirectional NS transport.
# -------------------------------------------------------------------------
parent_frame_path = ROOT / neighbors["frame"]
parent_frame = matrix(ZZ, [
    [ZZ(value) for value in line.split()]
    for line in parent_frame_path.read_text().splitlines()
    if line.strip() and not line.lstrip().startswith("#")
])
g_parent = block_diagonal_matrix(matrix(ZZ, ((0, 1), (1, 0))), -parent_frame)
selected = zero_frames["selected"]
forward = matrix(ZZ, selected["equation_A11_to_explicit_zero_basis"])
inverse = matrix(ZZ, selected["explicit_zero_to_equation_A11_basis"])
child_frame = matrix(ZZ, selected["frame"])
g_child = block_diagonal_matrix(matrix(ZZ, ((0, 1), (1, 0))), -child_frame)
if forward * inverse != identity_matrix(ZZ, 19) or inverse * forward != identity_matrix(ZZ, 19):
    raise ArithmeticError("selected q8 NS transports are not inverse")
if abs(forward.det()) != 1 or abs(inverse.det()) != 1:
    raise ArithmeticError("selected q8 NS transport is not unimodular")
if forward * g_parent * forward.transpose() != g_child:
    raise ArithmeticError("selected q8 NS transport misses the Gram form")

fibre = vector(ZZ, [1, 0] + [0] * 17)
root_curves = {
    index: vector(ZZ, [0, 0] + [-ZZ(index == other) for other in range(17)])
    for index in range(11)
}
affine_curve = vector(ZZ, [1, 0] + [1] * 11 + [0] * 6)
physical = {**root_curves, 11: affine_curve}
degrees = list(map(int, lattice["selection"]["explicit_curve_degrees"]["physical_old_A11_fibre_components"]))
if degrees != [0] * 9 + [1, 0, 1]:
    raise ArithmeticError(f"unexpected old-I12 q8 degrees: {degrees}")

edges = []
for left in range(12):
    for right in range(left + 1, 12):
        if physical[left] * g_parent * physical[right] == 1:
            edges.append((left, right))
expected_edges = [
    (0, 3), (0, 11), (1, 2), (1, 11), (2, 6), (3, 4),
    (4, 5), (5, 10), (6, 7), (7, 8), (8, 9), (9, 10),
]
if edges != expected_edges:
    raise ArithmeticError(f"old-I12 physical cycle changed: {edges}")

chains = ((0, 3, 4, 5, 10), (1, 2, 6, 7, 8))
chain_matrix = matrix(ZZ, [list(physical[index]) for chain in chains for index in chain])
chain_cartan = -(chain_matrix * g_parent * chain_matrix.transpose())
A5 = matrix(ZZ, 5, 5, lambda i, j: 2 if i == j else -1 if abs(i - j) == 1 else 0)
if chain_cartan != block_diagonal_matrix(A5, A5):
    raise ArithmeticError("degree-zero old-I12 components do not form two A5 chains")

component9_child = root_curves[9] * inverse
affine_child = affine_curve * inverse
if component9_child != vector(ZZ, [-1, 1] + [0] * 17):
    raise ArithmeticError("compiler-selected component 9 is not the child zero")
if affine_child * g_child * fibre != 1 or affine_child * g_child * affine_child != -2:
    raise ArithmeticError("old affine curve is not a child section")


def rational_record(value):
    value = KT(value)
    numerator = TQ(value.numerator())
    denominator = TQ(value.denominator())
    common = numerator.gcd(denominator)
    numerator //= common
    denominator //= common
    if denominator.leading_coefficient() < 0:
        numerator = -numerator
        denominator = -denominator
    return {
        "numerator_coefficients_low_to_high": [str(value) for value in numerator.list()],
        "denominator_coefficients_low_to_high": [str(value) for value in denominator.list()],
        "degrees_numerator_denominator": [int(numerator.degree()), int(denominator.degree())],
    }


def rows(value):
    return [[int(item) for item in row] for row in value.rows()]


physical_child_coordinates = {
    ("old_A11_affine" if index == 11 else f"old_A11_component_{index}"):
        [int(value) for value in curve * inverse]
    for index, curve in physical.items()
}

payload = {
    "schema": "elkies-k3.h3-q24-a11-to-2a5-q8-equation-marking-qq.v1",
    "status": "PASS_EXACT_Q24_A11_Q8_2A5_EQUATION_MARKING",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
    "old_I12": {
        "base_value": str(beta),
        "horizontal_specialization": "finite node",
        "node_projective_X_Y_Z": [str(x_beta), str(y_beta), str(z_beta)],
        "physical_cycle_edges": [list(edge) for edge in edges],
        "q8_degrees_component0_through_component10_then_affine": degrees,
    },
    "pointed_quartic": {
        "parameter_convention": "N=AA1-T*AA0, Db=T*BB0-BB1",
        "affine_ordinate": rational_record(w_affine),
        "old_A11_component_9_ordinate": rational_record(w_component9),
        "affine_sign_derivation": "w=(Db^2/Z)*(m^2-3*x_node)=(N^2-3*X*Db^2)/Z^3",
        "quartic_four_jet_identity": True,
        "selected_origin": "old_A11_component_9",
        "selected_origin_maps_to": "point at infinity on the pointed generalized Weierstrass model",
        "invariant_short_scaling": "x_child=9*x_short, y_child=27*y_short",
        "invariant_jacobian_identity": True,
    },
    "old_A11_affine_section_on_component9_pointed_child": {
        "x": rational_record(x_affine),
        "y": rational_record(y_affine),
        "exact_child_identity": True,
        "NS_coordinates_in_selected_child_basis": [int(value) for value in affine_child],
    },
    "physical_2A5": {
        "chains": [list(chain) for chain in chains],
        "interpretation": "old A11 component indices, in physical adjacency order",
        "cartan": rows(chain_cartan),
        "root_rank": 10,
        "root_determinant": 36,
        "child_coordinates": physical_child_coordinates,
    },
    "transport": {
        "equation_A11_to_component9_zero_2A5_basis": rows(forward),
        "component9_zero_2A5_to_equation_A11_basis": rows(inverse),
        "forward_determinant": int(forward.det()),
        "inverse_determinant": int(inverse.det()),
        "inverse_exact": True,
        "Gram_transport_exact": True,
    },
    "resolved_RR_dimensions": {
        "ambient": 14,
        "collision_rank": 12,
        "h0": 2,
    },
    "large_Groebner_required": False,
    "next_required": "start the promoted q6 orbit1307 equation lift from this exact component9-zero 2A5 model",
    "proof_boundary": (
        "The exact q8 equation is now pointed at old_A11_component_9, its complementary old affine "
        "section is exact on the same short Jacobian, all ten physical root components are attached "
        "as two ordered A5 chains, and the full NS transport is exact in both directions. This promotes "
        "only A11--q8 orbit12-->2A5/MW7. The promoted q6/q4/q6 suffix remains a lattice route until "
        "its own equation lifts are executed."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A11Q8MARKQQ|zero=old_A11_component_9|other=old_A11_affine|"
    "chains=0-3-4-5-10,1-2-6-7-8|det={}/{}|ambient=14|rank=12|h0=2|status={}".format(
        forward.det(), inverse.det(), payload["status"]
    ),
    flush=True,
)
print(f"OUTPUT|{OUTPUT.resolve()}", flush=True)
