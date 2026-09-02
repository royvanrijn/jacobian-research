#!/usr/bin/env sage-python
"""Compile the genuine orbit-96 A7+D7 pencil and its MW Galois action.

status: ACTIVE_PROOF
claim: exact physical orbit-96 equation and 2+chi_-3 MW decomposition
inputs: elkies-k3-e6a1-rho19 dissection and genuine-q2 certificates
outputs: elkies-k3-e6a1-rho19-orbit96-rr-galois-v1.json

The earlier tangent audit left the simplified tangent slope in a fraction
field.  In that parent Sage's zero-argument ``discriminant()`` does not compute
the discriminant of the quadratic in the old x-coordinate.  Coercing the
slope back to the old-base polynomial ring gives the genuine binary quartic.
Its Jacobian has fibres I8+I3*+7I1, as required by the orbit-96 lattice split.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import (
    CartanMatrix,
    IntegralLattice,
    PolynomialRing,
    PowerSeriesRing,
    QQ,
    ZZ,
    block_diagonal_matrix,
    identity_matrix,
    lcm,
    matrix,
    vector,
)


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
GEN = ROOT / "artifacts/generated-results"
DISSECTION = GEN / "elkies-k3-e6a1-rho19-k3-dissection-v1.json"
Q2_SOURCE = GEN / "elkies-k3-e6a1-rho19-genuine-q2-neighbors-v1.json"
DEFAULT_OUTPUT = GEN / "elkies-k3-e6a1-rho19-orbit96-rr-galois-v1.json"

exec(
    compile(
        (HERE / "elliptic_neighbor_compiler.sage").read_text(),
        str(HERE / "elliptic_neighbor_compiler.sage"),
        "exec",
    ),
    globals(),
)
exec(
    compile(
        (HERE / "exact_neighbor_engine.sage").read_text(),
        str(HERE / "exact_neighbor_engine.sage"),
        "exec",
    ),
    globals(),
)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def entries(value):
    return [int(item) for item in vector(ZZ, value)]


def rational_rows(value):
    return [[str(item) for item in row] for row in matrix(QQ, value).rows()]


def polynomial_coefficients(value):
    value = value.parent()(value)
    return [str(value[index]) for index in range(value.degree() + 1)]


def rational_leading_term_at_zero(value, uniformizer):
    """Return the order and leading coefficient in a local parameter."""
    field = uniformizer.parent().fraction_field()
    value = field(value)
    numerator = uniformizer.parent()(value.numerator())
    denominator = uniformizer.parent()(value.denominator())
    numerator_order = numerator.valuation()
    denominator_order = denominator.valuation()
    numerator_lead = (numerator / uniformizer**numerator_order)(0)
    denominator_lead = (denominator / uniformizer**denominator_order)(0)
    return numerator_order - denominator_order, numerator_lead / denominator_lead


def permutation(source_to_chart):
    return matrix(
        ZZ,
        6,
        6,
        lambda source, chart: ZZ(source_to_chart[source] == chart),
    )


for source_path in (DISSECTION, Q2_SOURCE):
    if not source_path.exists():
        raise FileNotFoundError(source_path)
dissection = json.loads(DISSECTION.read_text())
q2_source = json.loads(Q2_SOURCE.read_text())
if dissection.get("status") != "PASS_EXACT_GENERIC_NS_T_AND_FOUR_SINGULAR_K3_BOUNDARIES":
    raise ArithmeticError("the E6+A1 dissection source is not exact")
if q2_source.get("status") != "PASS_EXACT_COMPLETE_GENUINE_Q2_CENSUS_AND_18_NEF_MW3_FRAMES":
    raise ArithmeticError("the genuine q=2 source is not exact")
target = q2_source["secondary_fibre_simple_compiler_target"]
if target["orbit"] != 96 or target["root_type"] != "A7+D7":
    raise ArithmeticError("the orbit-96 target changed")

ns = matrix(ZZ, dissection["generic_k3"]["ns_gram"])
divisor = vector(ZZ, target["divisor_in_ns_basis"])
expected_divisor = vector(
    ZZ, [0, 0, 2, 1, 2, 2, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0]
)
if divisor != expected_divisor or divisor * ns * divisor != 0:
    raise ArithmeticError("the orbit-96 divisor changed")


# Physical IV* marking.  The old abstract E6 order is Sage's Cartan order.
# The ordinary resolution order is leaf,+outer,-outer,+inner,-inner,central.
source_e6_cartan = CartanMatrix(["E", 6])
chart_cartan = matrix(
    ZZ,
    [
        [2, 0, 0, 0, 0, -1],
        [0, 2, 0, -1, 0, 0],
        [0, 0, 2, 0, -1, 0],
        [0, -1, 0, 2, 0, -1],
        [0, 0, -1, 0, 2, -1],
        [-1, 0, 0, -1, -1, 2],
    ],
)
source_to_chart = (1, 0, 3, 5, 4, 2)  # zero-based: (2,1,4,6,5,3)
change = permutation(source_to_chart)
if change * chart_cartan * change.transpose() != source_e6_cartan:
    raise ArithmeticError("the physical E6 graph attachment changed")

# P0 has (x/u^2,y/u^2)=(-1,c*lambda) and hence meets +outer.  This fixes the
# E6 arm involution and sends abstract root one to chart component +outer.
if source_to_chart[0] != 1:
    raise ArithmeticError("P0 no longer orients the physical E6 marking")
chart_u_orders = vector(ZZ, (2, 1, 1, 2, 2, 3))
chart_x_orders = vector(ZZ, (2, 2, 2, 3, 3, 4))
chart_y_orders = vector(ZZ, (3, 2, 2, 4, 4, 6))
if chart_cartan * chart_u_orders != vector(ZZ, (1, 0, 0, 0, 0, 0)):
    raise ArithmeticError("the IV* base-parameter valuations changed")
if chart_cartan * chart_x_orders != vector(ZZ, (0, 1, 1, 0, 0, 0)):
    raise ArithmeticError("the IV* x valuations changed")
if chart_cartan * chart_y_orders != vector(ZZ, (0, 0, 0, 0, 0, 1)):
    raise ArithmeticError("the IV* y valuations changed")

# For g=(y+y0+m(x-x0))/(x-x0)^2 the exact numerator orders are obtained by
# the two ordinary blowups.  On the minus outer arm both the constant and
# first exceptional coefficients cancel; on the minus inner arm exactly the
# first coefficient cancels.  These give the complete pole cycle below.
# Replay those cancellations over the function field of a generic component.
LOCAL_K_RING = PolynomialRing(QQ, "kappa")
kappa_polynomial = LOCAL_K_RING.gen()
LOCAL_K = LOCAL_K_RING.fraction_field()
GENERIC_RING = PolynomialRing(LOCAL_K, "U")
U = GENERIC_RING.gen()
GENERIC_FIELD = GENERIC_RING.fraction_field()
SERIES_RING = PowerSeriesRing(GENERIC_FIELD, "epsilon", default_prec=12)
epsilon = SERIES_RING.gen()
kappa = GENERIC_FIELD(kappa_polynomial)
local_D = 3 * kappa**2 - 4
local_c = 2 * kappa / local_D
local_lam = -(kappa**2 - 4) * local_D / 4
local_d = local_c * local_lam

# Second blowup, generic points of the two outer arms.
local_u = epsilon * U
local_a = local_lam - 3 * local_u
local_b = local_d**2 + local_lam * local_u - 2 * local_u**2
outer_square = (
    epsilon**2 * U
    + epsilon * U**2 * local_a
    + U**2 * local_b
).sqrt()
if outer_square[0] != -local_d * U:
    raise ArithmeticError("the IV* outer square-root orientation changed")
outer_numerator_orders = []
for branch_y in (-outer_square, outer_square):  # plus, minus
    branch_x = epsilon**2 * U
    branch_old_y = epsilon**2 * U * branch_y
    branch_numerator = (
        branch_old_y
        + local_d * local_u**2
        + local_u / (2 * local_c) * (branch_x + local_u**2)
    )
    outer_numerator_orders.append(int(branch_numerator.valuation()))
if outer_numerator_orders != [2, 4]:
    raise ArithmeticError("the IV* outer-arm tangent cancellations changed")

# Third blowup, generic points of the two inner arms.
local_u = epsilon**2 * U
branch_x = epsilon**3 * U**2
local_a = local_lam - 3 * local_u
local_b = local_d**2 + local_lam * local_u - 2 * local_u**2
inner_square = (epsilon * U**2 + epsilon * U * local_a + local_b).sqrt()
if inner_square[0] != -local_d:
    raise ArithmeticError("the IV* inner square-root orientation changed")
inner_numerator_orders = []
for branch_y in (-inner_square, inner_square):  # plus, minus
    branch_old_y = epsilon**4 * U**2 * branch_y
    branch_numerator = (
        branch_old_y
        + local_d * local_u**2
        + local_u / (2 * local_c) * (branch_x + local_u**2)
    )
    inner_numerator_orders.append(int(branch_numerator.valuation()))
if inner_numerator_orders != [4, 5]:
    raise ArithmeticError("the IV* inner-arm tangent cancellations changed")

denominator_orders = 2 * vector(ZZ, (2, 2, 2, 3, 3, 4))
numerator_orders = vector(ZZ, (3, 2, 4, 4, 5, 6))
large_chart_cycle = denominator_orders - numerator_orders
small_chart_cycle = vector(
    ZZ, [max(0, value) for value in (large_chart_cycle - chart_u_orders)]
)
if large_chart_cycle != vector(ZZ, (1, 2, 0, 2, 1, 2)):
    raise ArithmeticError("the unregularized P0 tangent pole cycle changed")
if small_chart_cycle != vector(ZZ, (0, 1, 0, 0, 0, 0)):
    raise ArithmeticError("the regularized P0 tangent pole cycle changed")
large_source_cycle = large_chart_cycle * change.transpose()
small_source_cycle = small_chart_cycle * change.transpose()
if large_source_cycle != vector(ZZ, (2, 1, 2, 2, 1, 0)):
    raise ArithmeticError("the first physical E6 divisor cycle changed")
if small_source_cycle != vector(ZZ, (1, 0, 0, 0, 0, 0)):
    raise ArithmeticError("the second physical E6 divisor cycle changed")
if vector(ZZ, divisor[2:8]) != large_source_cycle:
    raise ArithmeticError("orbit 96 missed the first physical E6 cycle")
if vector(ZZ, divisor[8:14]) != small_source_cycle:
    raise ArithmeticError("orbit 96 missed the second physical E6 cycle")
if vector(ZZ, divisor[14:17]) != vector(ZZ, (0, 1, 0)):
    raise ArithmeticError("orbit 96 lost its middle I4 component")


# Source equation and the resolved tangent pencil.
K_RING = PolynomialRing(QQ, "k")
k_polynomial = K_RING.gen()
K = K_RING.fraction_field()
Z_RING = PolynomialRing(K, "z")
z = Z_RING.gen()
KZ = Z_RING.fraction_field()
T_RING = PolynomialRing(KZ, "t")
t = T_RING.gen()
X_RING = PolynomialRing(T_RING, "x")
x = X_RING.gen()
k = KZ(k_polynomial)
D_parameter = 3 * k**2 - 4
c = 2 * k / D_parameter
lam = -(k**2 - 4) * D_parameter / 4
H = 1 - t**2
old_a = H**3 * (lam - 3 * H)
old_b = H**4 * (c**2 * lam**2 + lam * H - 2 * H**2)
x0 = -H**2
y0 = c * lam * H**2

# This coercion is essential.  Without it the ambient x-polynomial is over
# Frac(KZ[t]), and Sage's zero-argument discriminant follows the wrong parent.
tangent_slope = T_RING((3 * x0**2 + old_a) / (2 * y0))
expected_slope = T_RING(((-3 * k**2 / 4 + 1) / k) * (t - 1) * (t + 1))
if tangent_slope != expected_slope:
    raise ArithmeticError("the polynomial tangent slope changed")
regularizer = t + 1
pole_factor = x - x0

# Resolve the old I4 fibre at t=infinity.  In minimal coordinates
# X=s^4*x,Y=s^6*y, the two side charts have exponent one and the middle
# component has exponent two.
XI_ETA_RING = PolynomialRing(K, names=("xi", "eta"))
xi, eta = XI_ETA_RING.gens()
XI_ETA_FIELD = XI_ETA_RING.fraction_field()
S_RING = PolynomialRing(XI_ETA_FIELD, "s")
s = S_RING.gen()
FS = S_RING.fraction_field()
k_infinity = FS(k_polynomial)
D_infinity = 3 * k_infinity**2 - 4
c_infinity = 2 * k_infinity / D_infinity
lam_infinity = -(k_infinity**2 - 4) * D_infinity / 4
q_infinity = s**2 - 1
infinity_a = lam_infinity * s**2 * q_infinity**3 - 3 * q_infinity**4
infinity_b = (
    c_infinity**2 * lam_infinity**2 * s**4 * q_infinity**4
    + lam_infinity * s**2 * q_infinity**5
    - 2 * q_infinity**6
)
x0_infinity = -q_infinity**2
y0_infinity = c_infinity * lam_infinity * s**2 * q_infinity**2
i4_orders = {}
for label, exponent in (("side", 1), ("middle", 2)):
    x_chart = -1 + s**exponent * xi
    y_chart = s**exponent * eta
    equation_chart = y_chart**2 - x_chart**3 - infinity_a * x_chart - infinity_b
    numerator_chart = (
        y_chart
        + y0_infinity
        + q_infinity / (2 * c_infinity) * (x_chart - x0_infinity)
    )
    z_chart = FS(
        s * (1 + s) * numerator_chart / (x_chart - x0_infinity) ** 2
    )
    equation_order, equation_initial = rational_leading_term_at_zero(
        equation_chart, s
    )
    z_order, z_initial = rational_leading_term_at_zero(z_chart, s)
    i4_orders[label] = (int(equation_order), int(z_order))
    if not equation_initial or not z_initial:
        raise ArithmeticError("an I4 initial form vanished")
if i4_orders != {"side": (2, 0), "middle": (4, -1)}:
    raise ArithmeticError(f"the resolved orbit-96 I4 orders changed: {i4_orders}")

scaled_line = regularizer * (y0 + tangent_slope * pole_factor)
cleared = (
    (z * pole_factor**2 - scaled_line) ** 2
    - regularizer**2 * (x**3 + old_a * x + old_b)
)
if cleared.parent() != X_RING:
    raise ArithmeticError("the tangent elimination escaped the polynomial parent")
quadratic, remainder = cleared.quo_rem(pole_factor**2)
if remainder or quadratic.degree() != 2:
    raise ArithmeticError("the orbit-96 pencil did not leave a quadratic in old x")
manual_discriminant = quadratic[1] ** 2 - 4 * quadratic[2] * quadratic[0]
if quadratic.discriminant() != manual_discriminant:
    raise ArithmeticError("the polynomial-parent quadratic discriminant changed")
quartic, square_factor = squarefree_binary_quartic(
    T_RING(manual_discriminant), T_RING
)
if square_factor**2 != regularizer**2 or quartic.degree() != 4:
    raise ArithmeticError("the genuine orbit-96 quartic changed")
if quartic(1) != 4 or quartic(-1) != 0:
    raise ArithmeticError("the rational orbit-96 quartic points changed")

invariant_i, invariant_j = binary_quartic_invariants(quartic)
child_a = Z_RING(-27 * invariant_i)
child_b = Z_RING(-27 * invariant_j)
if child_a.degree() != 6 or child_b.degree() != 9:
    raise ArithmeticError("the orbit-96 Weierstrass degrees changed")
if child_a[0] != -27 or child_b[0] != 54:
    raise ArithmeticError("the I8 normalization at z=0 changed")

classification = classify_finite_short_weierstrass_fibres(
    Z_RING, child_a, child_b
)
finite_profile = sorted(
    (item["kodaira"], item["degree"])
    for item in classification["finite_fibres"]
)
if finite_profile != [("I1", 7), ("I8", 1)]:
    raise ArithmeticError(f"unexpected finite orbit-96 fibres: {finite_profile}")
if classification["infinity_boundary"]["normalized_orders"] != (2, 3, 9):
    raise ArithmeticError("the orbit-96 infinity fibre is no longer I3*")
if kodaira_data_from_short_orders(2, 3, 9) != (7, 9, 4, "I3*"):
    raise ArithmeticError("the infinity D7 Kodaira lookup changed")
if (
    classification["finite_root_rank"] + 7,
    classification["finite_euler_number"] + 9,
    classification["finite_root_determinant"] * 4,
) != (14, 24, 32):
    raise ArithmeticError("the global A7+D7 fibre data changed")
reduced_discriminant = Z_RING(4 * child_a**3 + 27 * child_b**2)
residual_i1 = Z_RING(reduced_discriminant / (reduced_discriminant[15] * z**8))
if residual_i1.degree() != 7 or residual_i1.gcd(z) != 1:
    raise ArithmeticError("the seven residual I1 fibres changed")


# Exact MW lattice and arithmetic Galois action.  Over QQ(k), every old IV*
# component and both old sections are rational.  The only nontrivial constant
# action is sqrt(-3)-conjugation on the I4 side components A3_1 <-> A3_3.
neighbor = primitive_hyperbolic_split(ns, divisor)
roots, root_basis, root_data = roots_and_data(neighbor["child_frame"])
if root_data != (14, 140, 32):
    raise ArithmeticError("the orbit-96 child root data changed")
root_rank = root_basis.rank()
saturated_roots = root_basis.row_module().saturation().basis_matrix()
smith, left, right = saturated_roots.smith_form()
if tuple(abs(ZZ(smith[i, i])) for i in range(root_rank)) != (1,) * root_rank:
    raise ArithmeticError("the orbit-96 root lattice lost primitivity")
completion = right.inverse()
adapted_basis = saturated_roots.stack(completion[root_rank:])
adapted_gram = adapted_basis * neighbor["child_frame"] * adapted_basis.transpose()
root_block = adapted_gram[:root_rank, :root_rank]
coupling = adapted_gram[:root_rank, root_rank:]
tail = adapted_gram[root_rank:, root_rank:]
height0 = tail - coupling.transpose() * root_block.inverse() * coupling
height_scale = lcm(value.denominator() for value in height0.list())
quotient_change = (height_scale * height0).change_ring(ZZ).LLL_gram().transpose()
adapted_basis = (
    block_diagonal_matrix(identity_matrix(ZZ, root_rank), quotient_change)
    * adapted_basis
)
adapted_gram = adapted_basis * neighbor["child_frame"] * adapted_basis.transpose()
root_block = adapted_gram[:root_rank, :root_rank]
coupling = adapted_gram[:root_rank, root_rank:]
tail = adapted_gram[root_rank:, root_rank:]
height = tail - coupling.transpose() * root_block.inverse() * coupling
if height != matrix(QQ, [[QQ(3) / 8, 0, 0], [0, 1, 0], [0, 0, 3]]):
    raise ArithmeticError("the orbit-96 MW height basis changed")

tail_frame = adapted_basis[root_rank:, :]
tail_parent = matrix(
    ZZ,
    [entries(lift_child_frame_vector(neighbor, row)) for row in tail_frame.rows()],
)
expected_tail_parent = matrix(
    ZZ,
    [
        [0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 1],
    ],
)
if tail_parent != expected_tail_parent:
    raise ArithmeticError("the orbit-96 MW quotient representatives changed")

galois = matrix(
    ZZ,
    19,
    19,
    lambda i, j: ZZ(
        (i == j and i not in (14, 16))
        or (i == 14 and j == 16)
        or (i == 16 and j == 14)
    ),
)
if galois * ns * galois.transpose() != ns or divisor * galois != divisor:
    raise ArithmeticError("sqrt(-3)-conjugation no longer preserves orbit 96")
adapted_inverse = adapted_basis.inverse()
galois_coordinates = []
for representative in tail_parent.rows():
    child_coordinates = transport_parent_vector_to_child(
        neighbor, vector(ZZ, representative) * galois
    )
    frame_coordinates = vector(ZZ, child_coordinates[2:]) * adapted_inverse
    galois_coordinates.append(frame_coordinates)
galois_coordinates = matrix(ZZ, galois_coordinates)
galois_on_mw = galois_coordinates[:, root_rank:]
if galois_on_mw != matrix(ZZ, ((1, 0, 0), (0, -1, 0), (0, 0, 1))):
    raise ArithmeticError("the orbit-96 MW Galois representation changed")
if galois_on_mw * height * galois_on_mw.transpose() != height:
    raise ArithmeticError("Galois stopped preserving the MW height lattice")

# Convert the three quotient directions into actual section divisor classes.
# E6a_6 is a rational degree-one curve and is chosen as the new zero.  Adding
# the indicated multiple of D makes each zero+quotient representative a root;
# nefness then forces it to be the effective section rather than its negative.
new_zero = vector(ZZ, [ZZ(index == 7) for index in range(19)])
if divisor * ns * new_zero != 1 or new_zero * ns * new_zero != -2:
    raise ArithmeticError("the selected orbit-96 zero is not a section")
section_classes = []
fibre_shifts = []
for representative in tail_parent.rows():
    base = new_zero + vector(ZZ, representative)
    shift = ZZ((-2 - base * ns * base) / 2)
    section = base + shift * divisor
    if section * ns * section != -2 or divisor * ns * section != 1:
        raise ArithmeticError("failed to reconstruct an orbit-96 generator section")
    section_classes.append(section)
    fibre_shifts.append(int(shift))
expected_sections = (
    vector(ZZ, [ZZ(index == 10) for index in range(19)]),
    vector(ZZ, [0, 0, 4, 2, 4, 4, 2, 1, 2, 0, -1, 0, 0, 0, 1, 2, 0, 4, 0]),
    vector(ZZ, [0, -1, 6, 3, 6, 6, 3, 1, 3, 0, -1, 0, 0, 0, 0, 3, 0, 6, 1]),
)
if tuple(section_classes) != expected_sections or fibre_shifts != [0, 2, 3]:
    raise ArithmeticError("the orbit-96 generator divisor classes changed")
if section_classes[0] * galois != section_classes[0]:
    raise ArithmeticError("the height-3/8 generator lost rationality")
if section_classes[2] * galois != section_classes[2]:
    raise ArithmeticError("the height-3 generator lost rationality")
if section_classes[1] * galois == section_classes[1]:
    raise ArithmeticError("the height-one generator became rational")


payload = {
    "schema": "elkies-k3.e6a1-rho19-orbit96-rr-galois.v1",
    "status": "PASS_EXACT_PHYSICAL_A7D7_WEIERSTRASS_AND_MW_2_PLUS_CHI_MINUS3",
    "inputs": {
        relative(DISSECTION): digest(DISSECTION),
        relative(Q2_SOURCE): digest(Q2_SOURCE),
    },
    "base_field": "QQ(k)(z)",
    "physical_divisor": {
        "orbit": 96,
        "class_in_old_ns_basis": entries(divisor),
        "decomposition": "2*P0+(2,1,2,2,1,0)_E6a+(1,0,0,0,0,0)_E6b+A3_2",
        "physical_E6_source_to_chart": [value + 1 for value in source_to_chart],
        "chart_order": ["leaf", "plus_outer", "minus_outer", "plus_inner", "minus_inner", "central"],
        "chart_coordinate_orders": {
            "u": entries(chart_u_orders),
            "x": entries(chart_x_orders),
            "y": entries(chart_y_orders),
        },
        "tangent_numerator_orders": entries(numerator_orders),
        "tangent_denominator_orders": entries(denominator_orders),
        "large_IVstar_pole_cycle": entries(large_chart_cycle),
        "small_IVstar_pole_cycle_after_t_plus_1": entries(small_chart_cycle),
        "I4_pole_cycle": [0, 1, 0],
        "h0": 2,
        "basis": ["1", "z"],
        "z": "(t+1)*(y+y0+m_tan*(x-x0))/(x-x0)^2",
        "m_tan": str(tangent_slope),
    },
    "elimination_parent_regression": {
        "required_parent": str(X_RING),
        "quadratic_discriminant_checked_as": "b^2-4*a*c",
        "historical_failure": (
            "Leaving m_tan in Frac(KZ[t]) changes the x-polynomial parent; "
            "zero-argument discriminant() then returned a different expression."
        ),
        "square_factor": str(square_factor),
        "quartic_coefficients_low_to_high": polynomial_coefficients(quartic),
        "rational_points": [["t=-1", "w=0"], ["t=1", "w=2"], ["t=1", "w=-2"]],
    },
    "weierstrass": {
        "equation": "y^2=x^3-27*I6(z)*x-27*J9(z)",
        "I6_coefficients_low_to_high": polynomial_coefficients(Z_RING(invariant_i)),
        "J9_coefficients_low_to_high": polynomial_coefficients(Z_RING(invariant_j)),
        "short_a_coefficients_low_to_high": polynomial_coefficients(child_a),
        "short_b_coefficients_low_to_high": polynomial_coefficients(child_b),
        "reduced_discriminant_factorization": str(reduced_discriminant.factor()),
        "residual_I1_coefficients_low_to_high": polynomial_coefficients(residual_i1.monic()),
        "fibres": "I8 at z=0, I3* at z=infinity, seven residual I1",
        "root_type": "A7+D7",
        "root_data": [14, 140, 32],
        "euler_number": 24,
    },
    "mordell_weil": {
        "geometric_rank": 3,
        "arithmetic_rank_over_QQ_k_z": 2,
        "torsion": "trivial",
        "height_gram": rational_rows(height),
        "regulator": str(height.det()),
        "new_zero_class_in_old_ns_basis": entries(new_zero),
        "generator_section_classes_in_old_ns_basis": [entries(item) for item in section_classes],
        "generator_heights": ["3/8", "1", "3"],
        "generator_fields": ["QQ(k)(z)", "QQ(k,sqrt(-3))(z)", "QQ(k)(z)"],
        "galois_action_in_height_basis": [[int(item) for item in row] for row in galois_on_mw.rows()],
        "rank_decomposition": "1 + chi_{-3} + 1",
        "field_proof": (
            "Over QQ(k) the split IV* components and old sections are rational; "
            "sqrt(-3)-conjugation only swaps A3_1 and A3_3.  Its exact action "
            "on the saturated MW quotient is diag(1,-1,1)."
        ),
    },
    "comparison_with_orbit103": {
        "orbit96": {
            "reducible_fibres": 2,
            "arithmetic_rank": 2,
            "nontrivial_character": "chi_{-3}",
        },
        "orbit103": {
            "reducible_fibres": 4,
            "arithmetic_rank": 2,
            "nontrivial_character": "chi_{-3}",
        },
        "conclusion": (
            "Orbit 96 is fibre-sparser but has neither a third rational MW "
            "direction nor a new quadratic character relative to orbit 103."
        ),
    },
    "proof_boundary": {
        "proved": (
            "Physical E6 attachment, full orbit-96 divisor, resolved H0 basis, "
            "binary quartic, A7+D7 Weierstrass equation, saturated MW height "
            "basis, and exact constant-field Galois action."
        ),
        "not_claimed": (
            "Polynomial Weierstrass coordinates for all three sections, "
            "specialized ranks at individual rational k,z values, or coefficient optimality."
        ),
    },
}

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()
output_path = arguments.output.resolve()
encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if not output_path.exists() or output_path.read_text() != encoded:
        raise SystemExit(f"stale artifact: {output_path}")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(encoded)

print(
    "E6A1O96|RR=2|quartic=4|fibres=I8+I3*+7I1|MW=3|"
    "rank_Qkz=2|galois=1+chi_-3+1|status=PASS_EXACT",
    flush=True,
)
print(f"OUTPUT|{output_path}", flush=True)
