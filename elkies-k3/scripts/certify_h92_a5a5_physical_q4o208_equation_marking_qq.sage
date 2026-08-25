#!/usr/bin/env sage -python
"""Point the physical q4/orbit208 quartic at C5 and attach its 3A3 marking.

The exact RR plane is spanned by ``g_i=a_i(T)+b_i*m``, where

    m = (y+Y_P1229)/(x-X_P1229),       U = g_1/g_0.

At each old split I6 fibre we use its resolved toric chart.  On oriented
component one over the root reducing to 89, both g_i vanish simply and their
leading-coefficient ratio is a Mobius function of the component parameter.
This is the exact new-base coordinate on old_A11_component_5.  The branch
ordinate

    W = (2*x+X_P1229-m^2)/((T-r_1)(T-r_2))

is linear in the same parameter.  Eliminating that parameter gives an exact
point on the quartic over QQ(U), which is selected as the new zero.

Only a six-jet at one old I6 root is used to point the quartic and recover the
opposite affine section.  Two further resolved leading terms attach the other
I4 supports.  The physical cycles and the complete determinant-one NS
transport are then replayed from the promoted lattice certificate.  No
Groebner basis and no large polynomial factorization are used.
"""

import hashlib
import json
import math
import time
from pathlib import Path

from sage.all import (
    EllipticCurve, GF, LaurentSeriesRing, PolynomialRing, QQ, ZZ, block_diagonal_matrix,
    identity_matrix, matrix, vector,
)


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
SURFACE = LOCAL / "q24-a11-to-2a5-q8-resolved-rr-qq.json"
P1229_PATH = LOCAL / "q24-2a5-p1229-scaled-x-qq.json"
RR_PATH = LOCAL / "q24-2a5-physical-q4o208-rr-qq.json"
MARKING = GENERATED / "elkies-k3-h3-a5a5-physical-component-chamber-marking.json"
ROUTE = GENERATED / "elkies-k3-h3-a5a5-physical-q4o208-to-pinned-r17-certificate.json"
OUTPUT = LOCAL / "q24-2a5-physical-q4o208-equation-marking-qq.json"
INPUTS = (SURFACE, P1229_PATH, RR_PATH, MARKING, ROUTE)

started = time.monotonic()


def log(stage, **fields):
    suffix = "|".join(f"{key}={value}" for key, value in fields.items())
    print(
        f"A5Q4O208MARKQQ|stage={stage}|elapsed={time.monotonic()-started:.3f}"
        + (f"|{suffix}" if suffix else ""),
        flush=True,
    )


surface = json.loads(SURFACE.read_text())
p1229 = json.loads(P1229_PATH.read_text())
rr = json.loads(RR_PATH.read_text())
marking = json.loads(MARKING.read_text())
route = json.loads(ROUTE.read_text())
assert surface["status"] == "PASS_EXACT_Q24_A11_Q8_2A5_RESOLVED_RR"
assert p1229["status"] == "PASS_EXACT_QQ_P1229_POLYNOMIAL_SECTION"
assert rr["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O208_3A3_RR_AND_JACOBIAN"
assert marking["status"] == "PASS_EXACT_A5A5_PHYSICAL_COMPONENT_CHAMBER_MARKING"
assert route["status"] == "PASS_EXACT_PHYSICAL_Q4O208_3A3_TO_PINNED_R17"
assert route["selection"]["zero"] == "old_A11_component_5"

RT = PolynomialRing(QQ, "T")
T = RT.gen()
RU = PolynomialRing(QQ, "U")
U = RU.gen()
KU = RU.fraction_field()


def poly_t(values):
    return RT([QQ(value) for value in values])


A = poly_t(surface["child"]["minimal_A_coefficients_low_to_high"])
B = poly_t(surface["child"]["minimal_B_coefficients_low_to_high"])
P = p1229["P1229"]
PX = poly_t(P["X_coefficients_low_to_high"])
PY = poly_t(P["Y_coefficients_low_to_high"])
assert PY**2 == PX**3 + A * PX + B

kernel = matrix(QQ, [[QQ(value) for value in row]
                     for row in rr["resolved_RR"]["kernel_basis"]])
assert kernel.nrows() == 2 and kernel.ncols() == 4 and kernel.rank() == 2
pairs = []
for row in kernel.rows():
    pairs.append((RT(row[0] + row[1] * T + row[2] * T**2), QQ(row[3])))
(a0, b0), (a1, b1) = pairs
denominator = KU(U) * KU(b0) - KU(b1)
assert denominator and b0

i6_roots = []
for record in surface["child"]["discriminant_factorization"]:
    if int(record["multiplicity"]) != 6:
        continue
    factor = RT(record["factor"])
    assert factor.degree() == 1
    i6_roots.append(QQ(-factor[0] / factor[1]))
assert len(i6_roots) == 2
Fp = GF(103)
assert sorted(int(Fp(root)) for root in i6_roots) == [68, 89]
roots_by_reduction = {int(Fp(root)): root for root in i6_roots}
log("LOAD", kernel_bits=rr["resolved_RR"]["maximum_kernel_rational_bits"])


# -------------------------------------------------------------------------
# Exact resolved I6 charts.
# -------------------------------------------------------------------------
r_ring = PolynomialRing(QQ, "r")
r_parameter = r_ring.gen()
Kr = r_ring.fraction_field()
PRECISION = 4
LS = LaurentSeriesRing(Kr, "s", default_prec=PRECISION)
s = LS.gen()


def shifted(poly, root):
    """Exact T=root+s jet in Kr((s))."""

    return LS([
        Kr(poly.derivative(order)(root) / math.factorial(order))
        for order in range(PRECISION)
    ])


def node_x(root):
    value = QQ(-3 * B(root) / (2 * A(root)))
    assert value**3 + A(root) * value + B(root) == 0
    assert 3 * value**2 + A(root) == 0
    return value


def oriented_rho(root, required_reduction):
    square = QQ(3 * node_x(root))
    if not square.is_square():
        raise ArithmeticError("split I6 tangent square is not rational")
    candidate = QQ(square.sqrt())
    for value in (candidate, -candidate):
        if Fp(value) == Fp(required_reduction):
            return value
    raise ArithmeticError("I6 tangent orientation does not lift the pinned mod-103 sign")


def component_values(root, rho0, component):
    """Return the two RR functions and raw branch ordinate on an I6 component."""

    Aloc = shifted(A, root)
    center = LS(Kr(node_x(root)))
    for unused in range(3):
        center = (center + (-Aloc / 3) / center) / 2

    aa = LS(Kr(r_parameter)) * s**component
    # The complementary toric coordinate has order 6-component.  For
    # components one and two this is beyond the precision needed for the
    # leading RR ratios, so its exact unit never enters this calculation.
    bb = LS.zero()
    yy = (aa + bb) / 2
    ww = (bb - aa) / 2
    rho = LS(Kr(rho0))
    for unused in range(3):
        rho -= (rho**3 - 3 * center * rho - ww) / (3 * rho**2 - 3 * center)
    xx = center + ww / rho
    chord = (yy + shifted(PY, root)) / (xx - shifted(PX, root))
    functions = [shifted(apoly, root) + LS(Kr(bvalue)) * chord
                 for apoly, bvalue in pairs]
    branch = 2 * xx + shifted(PX, root) - chord**2
    return functions, branch


def first_ratio(values):
    valuations = [int(value.valuation()) for value in values]
    if valuations[0] != valuations[1]:
        raise ArithmeticError(f"RR functions have unequal component orders {valuations}")
    order = valuations[0]
    return Kr(values[1][order] / values[0][order]), order


beta = roots_by_reduction[89]
gamma = roots_by_reduction[68]
rho_beta = oriented_rho(beta, 35)
rho_gamma = oriented_rho(gamma, 95)

# Component one at beta is old_A11_component_5.  Its new-base coordinate is
# genuinely variable, while its normalized branch ordinate is affine-linear.
c5_functions, c5_branch = component_values(beta, rho_beta, 1)
U_c5_of_r, c5_order = first_ratio(c5_functions)
if c5_order != 1 or int(c5_branch.valuation()) != 1:
    raise ArithmeticError("C5 does not have the expected simple resolved q4 slice")
other_root = gamma
W_c5_of_r = Kr(c5_branch[1] / (beta - other_root))

u_num = r_ring(U_c5_of_r.numerator())
u_den = r_ring(U_c5_of_r.denominator())
if (u_num.degree(), u_den.degree()) != (1, 1):
    raise ArithmeticError("C5 new-base map is not Mobius")
n0, n1 = QQ(u_num[0]), QQ(u_num[1])
d0, d1 = QQ(u_den[0]), QQ(u_den[1])
r_of_U = KU((KU(n0) - KU(U) * KU(d0)) / (KU(U) * KU(d1) - KU(n1)))


def evaluate_rational_at_u(value):
    value = Kr(value)
    numerator = r_ring(value.numerator())
    denominator_r = r_ring(value.denominator())
    return KU(numerator(r_of_U) / denominator_r(r_of_U))


W_c5 = evaluate_rational_at_u(W_c5_of_r)
if KU(U_c5_of_r.numerator()(r_of_U) / U_c5_of_r.denominator()(r_of_U)) != KU(U):
    raise ArithmeticError("Mobius inversion on C5 failed")
log(
    "C5_SLICE", U_degrees=f"{u_num.degree()}/{u_den.degree()}",
    W_degrees=f"{r_ring(W_c5_of_r.numerator()).degree()}/"
              f"{r_ring(W_c5_of_r.denominator()).degree()}",
)

# Component one at gamma is old_A11_component_7 in the pinned oriented I6
# chain (component two below is old_A11_component_6).  Its resolved slice
# labels the two exact signs over T=gamma: the opposite sign is the second
# I6 affine component.
c7_functions, c7_branch = component_values(gamma, rho_gamma, 1)
U_c7_of_r, c7_order = first_ratio(c7_functions)
if c7_order != 1 or int(c7_branch.valuation()) != 1:
    raise ArithmeticError("C7 does not have the expected simple resolved q4 slice")
W_c7_of_r = Kr(c7_branch[1] / (gamma - beta))
c7_u_num = r_ring(U_c7_of_r.numerator())
c7_u_den = r_ring(U_c7_of_r.denominator())
if (c7_u_num.degree(), c7_u_den.degree()) != (1, 1):
    raise ArithmeticError("C7 new-base map is not Mobius")
c7_n0, c7_n1 = QQ(c7_u_num[0]), QQ(c7_u_num[1])
c7_d0, c7_d1 = QQ(c7_u_den[0]), QQ(c7_u_den[1])
r_c7_of_U = KU(
    (KU(c7_n0) - KU(U) * KU(c7_d0))
    / (KU(U) * KU(c7_d1) - KU(c7_n1))
)
c7_w_num = r_ring(W_c7_of_r.numerator())
c7_w_den = r_ring(W_c7_of_r.denominator())
W_c7 = KU(c7_w_num(r_c7_of_U) / c7_w_den(r_c7_of_U))
if KU(U_c7_of_r.numerator()(r_c7_of_U) / U_c7_of_r.denominator()(r_c7_of_U)) != KU(U):
    raise ArithmeticError("Mobius inversion on C7 failed")
log(
    "C7_SLICE", U_degrees=f"{c7_u_num.degree()}/{c7_u_den.degree()}",
    W_degrees=f"{c7_w_num.degree()}/{c7_w_den.degree()}",
)

# Components two at beta and gamma are C4 and C6.  Their ratios are constant
# and attach the two non-special I4 supports without solving a global system.
first_functions, unused = component_values(beta, rho_beta, 2)
second_functions, unused = component_values(gamma, rho_gamma, 2)
U_first_r, first_order = first_ratio(first_functions)
U_second_r, second_order = first_ratio(second_functions)
if r_ring(U_first_r.numerator()).degree() or r_ring(U_first_r.denominator()).degree():
    raise ArithmeticError("C4 support still depends on its component parameter")
if r_ring(U_second_r.numerator()).degree() or r_ring(U_second_r.denominator()).degree():
    raise ArithmeticError("C6 support still depends on its component parameter")
U_first = QQ(U_first_r)
U_second = QQ(U_second_r)

# Along P1229 the chord has a pole, hence U tends to b1/b0.
U_special = QQ(b1 / b0)
supports = [U_special, U_first, U_second]
if len(set(supports)) != 3:
    raise ArithmeticError("the three physical I4 supports are not distinct")
repeated_support = RU(rr["child"]["finite_repeated_support"])
for support in supports:
    if repeated_support(support):
        raise ArithmeticError("physical I4 support misses the exact discriminant support")
log(
    "SUPPORTS", special_mod103=int(Fp(U_special)),
    first_mod103=int(Fp(U_first)), second_mod103=int(Fp(U_second)),
    component_orders=f"{first_order},{second_order}",
)


# -------------------------------------------------------------------------
# Six-jet of the quartic at beta and exact C5 pointing.
# -------------------------------------------------------------------------
JET_PRECISION = 7
JS = LaurentSeriesRing(KU, "z", default_prec=JET_PRECISION)
z = JS.gen()


def shifted_ku(poly, root):
    return JS([
        KU(poly.derivative(order)(root) / math.factorial(order))
        for order in range(JET_PRECISION)
    ])


a0s, a1s = shifted_ku(a0, beta), shifted_ku(a1, beta)
PXs, PYs, As = shifted_ku(PX, beta), shifted_ku(PY, beta), shifted_ku(A, beta)
chord = (a1s - KU(U) * a0s) / denominator
radicand = chord**4 - 6 * PXs * chord**2 - 8 * PYs * chord - 3 * PXs**2 - 4 * As
old_square = (JS(KU(beta)) + z - KU(beta)) * (JS(KU(beta)) + z - KU(gamma))
quartic_jet = radicand / old_square**2
q = [KU(quartic_jet[index]) for index in range(5)]
if q[0] != W_c5**2:
    raise ArithmeticError("C5 resolved ordinate misses the exact q4 quartic")

e, d, c, b, a = q
w0 = W_c5
a_1 = d / w0
a_2 = c - d**2 / (4 * w0**2)
a_3 = 2 * w0 * b
a_4 = -4 * w0**2 * a
a_6 = a_2 * a_4
b_2 = a_1**2 + 4 * a_2
b_4 = 2 * a_4 + a_1 * a_3
b_6 = a_3**2 + 4 * a_6
c_4 = b_2**2 - 24 * b_4
c_6 = -b_2**3 + 36 * b_2 * b_4 - 216 * b_6
pointed_A = -c_4 / 48
pointed_B = -c_6 / 864

I = 12 * a * e - 3 * b * d + c**2
J = 72 * a * c * e + 9 * b * c * d - 27 * a * d**2 - 27 * b**2 * e - 2 * c**3
if 81 * pointed_A != -27 * I or 729 * pointed_B != -27 * J:
    raise ArithmeticError("C5 pointed quartic misses its invariant Jacobian")

A_child = KU(RU([QQ(value) for value in rr["child"]["minimal_A_coefficients_low_to_high"]]))
B_child = KU(RU([QQ(value) for value in rr["child"]["minimal_B_coefficients_low_to_high"]]))
if 81 * pointed_A * denominator**8 != A_child:
    raise ArithmeticError("C5 pointed A invariant misses the global q4 Jacobian")
if 729 * pointed_B * denominator**12 != B_child:
    raise ArithmeticError("C5 pointed B invariant misses the global q4 Jacobian")

# The opposite quartic sign is the first old-I6 affine component.  Pointing
# at C5 turns it into the following exact section of the global short model.
x_general = -a_2
y_general = a_1 * a_2 - a_3
x_affine_raw = KU(9 * (x_general + b_2 / 12))
y_affine_raw = KU(27 * (y_general + (a_1 * x_general + a_3) / 2))
x_affine = KU(denominator**4 * x_affine_raw)
y_affine = KU(denominator**6 * y_affine_raw)
if y_affine**2 != x_affine**3 + A_child * x_affine + B_child:
    raise ArithmeticError("opposite old-I6 affine component misses the q4 child")

# Evaluate the same pointed quartic map at the resolved C7 slice.  This
# labels, over QQ(U), both constant-old-base sections above T=gamma.
z_gamma = KU(gamma - beta)
if W_c7**2 != sum((q[index] * z_gamma**index for index in range(5)), KU.zero()):
    raise ArithmeticError("C7 resolved ordinate misses the exact q4 quartic")


def pointed_coordinates(z_value, ordinate):
    x_value = (2 * w0 * (ordinate + w0) + d * z_value) / z_value**2
    y_value = (
        4 * w0**2 * (ordinate + w0) + 2 * w0 * d * z_value
        + (2 * w0 * c - d**2 / (2 * w0)) * z_value**2
    ) / z_value**3
    return (
        KU(denominator**4 * 9 * (x_value + b_2 / 12)),
        KU(denominator**6 * 27 * (y_value + (a_1 * x_value + a_3) / 2)),
    )


x_c7, y_c7 = pointed_coordinates(z_gamma, W_c7)
x_second_affine, y_second_affine = pointed_coordinates(z_gamma, -W_c7)
E_child = EllipticCurve(KU, [0, 0, 0, A_child, B_child])
C7_point = E_child(x_c7, y_c7)
second_affine_point = E_child(x_second_affine, y_second_affine)
first_affine_point = E_child(x_affine, y_affine)
if C7_point + second_affine_point != first_affine_point:
    raise ArithmeticError("C7 plus second affine misses the first affine section")
log(
    "POINT", quartic_identity=True, global_jacobian=True,
    opposite_section=True, gamma_pair_labeled=True,
)


# -------------------------------------------------------------------------
# Physical 3A3 cycles and bidirectional NS transport.
# -------------------------------------------------------------------------
source_frame_path = ROOT / marking["frame_output"]
source_frame = matrix(ZZ, [
    [ZZ(value) for value in line.split()]
    for line in source_frame_path.read_text().splitlines()
    if line.strip() and not line.lstrip().startswith("#")
])
U2 = matrix(ZZ, ((0, 1), (1, 0)))
g_parent = block_diagonal_matrix(U2, -source_frame)
child_to_parent = matrix(ZZ, route["transport"]["parent_to_effective_zero_child_basis"])
parent_to_child = matrix(ZZ, route["transport"]["effective_zero_child_to_parent_basis"])
if child_to_parent * parent_to_child != identity_matrix(ZZ, 19):
    raise ArithmeticError("q4 parent/child transports are not inverse")
if parent_to_child * child_to_parent != identity_matrix(ZZ, 19):
    raise ArithmeticError("q4 child/parent transports are not inverse")
if abs(child_to_parent.det()) != 1 or abs(parent_to_child.det()) != 1:
    raise ArithmeticError("q4 NS transport is not unimodular")
g_child = child_to_parent * g_parent * child_to_parent.transpose()
expected_child_gram = block_diagonal_matrix(U2, -matrix(ZZ, route["selection"]["child_frame"]))
if g_child != expected_child_gram:
    raise ArithmeticError("q4 transport misses the selected child Gram form")

curves = {
    name: vector(ZZ, value)
    for name, value in marking["equation_explicit_curves_in_child"].items()
}
assert curves["old_zero"] == curves["old_A11_component_9"]
F_parent = vector(ZZ, route["fibre"]["class_in_parent"])
F_child = vector(ZZ, [1, 0] + [0] * 17)
assert vector(ZZ, child_to_parent.row(0)) == F_parent


def child(curve):
    return vector(ZZ, curve * parent_to_child)


zero_c5 = child(curves["old_A11_component_5"])
if zero_c5 != vector(ZZ, [-1, 1] + [0] * 17):
    raise ArithmeticError("C5 is not the selected child zero")

cycle_names = {
    "special_I4": [
        "old_zero", "old_A11_component_10", "P1229", "old_A11_component_8",
    ],
    "first_old_I6_I4": [
        "old_A11_component_0", "old_A11_component_3", "old_A11_component_4",
    ],
    "second_old_I6_I4": [
        "old_A11_component_1", "old_A11_component_2", "old_A11_component_6",
    ],
}
cycles = {}
for name, names in cycle_names.items():
    values = [child(curves[label]) for label in names]
    if len(values) == 3:
        values.append(F_child - sum(values, vector(ZZ, 19)))
    cycles[name] = values


def verify_i4(values):
    if sum(values, vector(ZZ, 19)) != F_child:
        raise ArithmeticError("physical I4 components do not sum to the new fibre")
    gram = matrix(ZZ, [[left * g_child * right for right in values] for left in values])
    expected = matrix(ZZ, 4, 4, lambda i, j: -2 if i == j else int((j-i) % 4 in (1, 3)))
    if gram != expected:
        raise ArithmeticError("physical components do not form an I4 cycle")


for values in cycles.values():
    verify_i4(values)
identity_indices = {}
root_components = []
for name, values in cycles.items():
    hits = [index for index, value in enumerate(values) if zero_c5 * g_child * value == 1]
    if len(hits) != 1:
        raise ArithmeticError(f"C5 does not select one identity component on {name}")
    identity = hits[0]
    identity_indices[name] = identity
    root_components.extend(values[(identity + offset) % 4] for offset in (1, 2, 3))
root_matrix = matrix(ZZ, root_components)
root_cartan = -(root_matrix * g_child * root_matrix.transpose())
A3 = matrix(ZZ, 3, 3, lambda i, j: 2 if i == j else -1 if abs(i-j) == 1 else 0)
if root_cartan != block_diagonal_matrix(A3, A3, A3):
    raise ArithmeticError("physical roots do not form 3A3")
if root_cartan.det() != 64:
    raise ArithmeticError("physical 3A3 root determinant changed")
log("MARKING", fibres="3I4", roots="3A3", det=f"{child_to_parent.det()}/{parent_to_child.det()}")


def rational_record(value):
    value = KU(value)
    numerator = RU(value.numerator())
    denominator_u = RU(value.denominator())
    common = numerator.gcd(denominator_u)
    numerator //= common
    denominator_u //= common
    if denominator_u.leading_coefficient() < 0:
        numerator, denominator_u = -numerator, -denominator_u
    return {
        "numerator_coefficients_low_to_high": [str(item) for item in numerator.list()],
        "denominator_coefficients_low_to_high": [str(item) for item in denominator_u.list()],
        "degrees_numerator_denominator": [int(numerator.degree()), int(denominator_u.degree())],
    }


def rows(value):
    return [[int(item) for item in row] for row in value.rows()]


def vectors(values):
    return [[int(item) for item in value] for value in values]


support_by_cycle = {
    "special_I4": U_special,
    "first_old_I6_I4": U_first,
    "second_old_I6_I4": U_second,
}
payload = {
    "schema": "elkies-k3.q24-2a5-physical-q4o208-equation-marking-qq.v1",
    "status": "PASS_EXACT_QQ_PHYSICAL_Q4O208_C5_EQUATION_MARKING",
    "selected_zero": "old_A11_component_5",
    "resolved_C5_slice": {
        "parent_I6_base_value": str(beta),
        "base_value_mod_103": 89,
        "oriented_component": 1,
        "oriented_rho": str(rho_beta),
        "new_base_map_on_C5": rational_record(KU(U_c5_of_r.numerator()(r_of_U) /
                                                      U_c5_of_r.denominator()(r_of_U))),
        "new_base_map_degrees_in_component_parameter": [1, 1],
        "quartic_ordinate_on_C5": rational_record(W_c5),
        "raw_ordinate_degree_in_component_parameter": [
            int(r_ring(W_c5_of_r.numerator()).degree()),
            int(r_ring(W_c5_of_r.denominator()).degree()),
        ],
        "exact_quartic_square_identity": True,
    },
    "pointed_quartic": {
        "translated_old_base_coordinate": "T-beta",
        "selected_origin": "old_A11_component_5",
        "selected_origin_maps_to": "point at infinity on the pointed generalized Weierstrass model",
        "invariant_short_scaling": (
            "first x=9*x_short,y=27*y_short; then x_global=(U*b0-b1)^4*x, "
            "y_global=(U*b0-b1)^6*y"
        ),
        "six_jet_used": True,
        "invariant_jacobian_identity": True,
    },
    "first_I6_affine_component_on_C5_pointed_child": {
        "x": rational_record(x_affine),
        "y": rational_record(y_affine),
        "exact_child_identity": True,
        "NS_coordinates": [int(item) for item in child(curves["first_I6_affine_component"])],
    },
    "resolved_C7_slice": {
        "parent_I6_base_value": str(gamma),
        "base_value_mod_103": 68,
        "oriented_component": 1,
        "oriented_rho": str(rho_gamma),
        "new_base_map_degrees_in_component_parameter": [1, 1],
        "quartic_ordinate_on_C7": rational_record(W_c7),
        "exact_quartic_square_identity": True,
    },
    "old_A11_component_7_on_C5_pointed_child": {
        "x": rational_record(x_c7),
        "y": rational_record(y_c7),
        "exact_child_identity": True,
        "NS_coordinates": [int(item) for item in child(curves["old_A11_component_7"])],
    },
    "second_I6_affine_component_on_C5_pointed_child": {
        "x": rational_record(x_second_affine),
        "y": rational_record(y_second_affine),
        "exact_child_identity": True,
        "NS_coordinates": [int(item) for item in child(curves["second_I6_affine_component"])],
        "exact_group_relation": "C7 + second_I6_affine = first_I6_affine",
    },
    "physical_fibres": {
        name: {
            "support": str(support_by_cycle[name]),
            "support_mod_103": int(Fp(support_by_cycle[name])),
            "components_in_cycle_order": vectors(values),
            "inherited_labels_before_missing_component": cycle_names[name],
            "missing_component_appended": len(cycle_names[name]) == 3,
            "identity_component_index": identity_indices[name],
        }
        for name, values in cycles.items()
    },
    "root_lattice": {
        "cartan": rows(root_cartan),
        "type": "3A3",
        "rank": 9,
        "determinant": 64,
        "MW_rank_if_rho19": 8,
    },
    "transport": {
        "C5_zero_child_to_physical_2A5_parent_basis": rows(child_to_parent),
        "physical_2A5_parent_to_C5_zero_child_basis": rows(parent_to_child),
        "forward_determinant": int(child_to_parent.det()),
        "inverse_determinant": int(parent_to_child.det()),
        "inverse_exact": True,
        "Gram_transport_exact": True,
        "canonical_3A3_and_pinned_R17_transport": str(ROUTE.relative_to(ROOT)),
        "canonical_3A3_and_pinned_R17_transport_sha256": hashlib.sha256(ROUTE.read_bytes()).hexdigest(),
    },
    "resolved_RR_dimensions": {"ambient": 4, "condition_rank": 2, "h0": 2},
    "method": {
        "large_Groebner_required": False,
        "large_polynomial_factorization_required": False,
        "global_quartic_reparse_required": False,
        "runtime_seconds": float(time.monotonic() - started),
    },
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS] + [str(source_frame_path.relative_to(ROOT))],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS + (source_frame_path,)
        },
    },
    "next_required": "start the exact canonical 3A3-to-A3+2A2 suffix equation lift",
    "proof_boundary": (
        "This points the exact q4/orbit208 quartic at the effective C5 curve, attaches the "
        "opposite first-I6 affine section, identifies all three exact I4 supports and physical "
        "cycles, verifies the 3A3 root lattice, and replays the bidirectional unimodular NS "
        "transport to the already-certified canonical 3A3/pinned-R17 route. It promotes only "
        "the 2A5--q4/orbit208-->3A3 equation edge; later suffix equations remain separate gates."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A5Q4O208MARKQQ|zero=C5|other=first_I6_affine|fibres=3I4|root=3A3|"
    "ambient=4|rank=2|h0=2|det={}/{}|status={}|output={}".format(
        child_to_parent.det(), parent_to_child.det(), payload["status"], OUTPUT,
    ),
    flush=True,
)
