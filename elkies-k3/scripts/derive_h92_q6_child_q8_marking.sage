#!/usr/bin/env sage -python
"""Pin the selected H3 q=8 divisor to an exact degree-two child marking.

The selected D13/MW4 q=8 class has MW projection ``(0,-2,0)`` with respect
to the canonical q=6 zero supplied by the lattice neighbour frame.  The
explicit child equation instead uses the transported old zero.  This script
derives that zero translation before identifying the q=8 marked section.
Two exact child-Jacobian points are already available: the old E7_7 and
affine-E7 components, transported through the binary-quartic covariants.
Their differences from the transported child zero give the first two MW
directions; after the zero translation, the selected q=8 marking is the
first direction plus twice the second direction.

It establishes only the degree-two *generic-fibre* marking on the q=6 child.
Resolved II*/IV* coefficient modules are still required to construct the q=8
pencil, and no rootless bisection or rank statement is made here.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import (
    QQ, ZZ, EllipticCurve, PolynomialRing, block_diagonal_matrix,
    gcd, identity_matrix, matrix, pari, vector, xgcd,
)


ROOT = Path(__file__).resolve().parents[2]
FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
ZERO = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-zero-section.json"
COMPONENTS = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-e7-infinity-sections.json"
Q8_ORBITS = ROOT / "artifacts/generated-results/elkies-k3-h3-q6-q8-orbits.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-marking.json"

MW_LIFTS = matrix(ZZ, [
    [-5, -4, -3, 0, 0, 0, 0, 0, 0, 0, 0, -4, 1, 0, -4, 2, -2],
    [-10, -8, -6, 0, 0, 0, 0, 0, 0, 0, 0, -8, 4, 1, -8, 5, -4],
    [-5, -4, -3, 0, 0, 0, 0, 0, 0, 0, 0, -3, 2, 0, -4, 2, -2],
])
EXPECTED_HEIGHT = matrix(QQ, [
    [QQ(8) / 3, QQ(1) / 3, -1],
    [QQ(1) / 3, QQ(8) / 3, 1],
    [-1, 1, 46],
])
REFLECTIONS = (1, 2, 4, 3, 5, 4, 2, 6, 5, 4, 3, 1, 7, 6, 5, 4, 2, 3, 4, 5, 6, 7)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ])


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def rational_function(field, ring, numerator, denominator):
    return field(polynomial(ring, numerator)) / field(polynomial(ring, denominator))


def isotropic_mate(ns, fiber):
    current = ZZ(0)
    entries = [ZZ(0)] * ns.nrows()
    for index, value in enumerate(ns * fiber):
        if not value:
            continue
        divisor, left, right = xgcd(current, ZZ(value))
        entries = [left * entry for entry in entries]
        entries[index] += right
        current = divisor
    assert abs(current) == 1
    if current == -1:
        entries = [-entry for entry in entries]
    mate = vector(ZZ, entries)
    mate -= (mate * ns * mate // 2) * fiber
    assert mate * ns * mate == 0 and mate * ns * fiber == 1
    return mate


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--child", type=Path, default=CHILD)
parser.add_argument("--zero", type=Path, default=ZERO)
parser.add_argument("--components", type=Path, default=COMPONENTS)
parser.add_argument("--q8-orbits", type=Path, default=Q8_ORBITS)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
for name in ("child", "zero", "components", "q8_orbits", "output"):
    setattr(args, name, getattr(args, name).resolve())

child = json.loads(args.child.read_text())
zero_payload = json.loads(args.zero.read_text())
component_payload = json.loads(args.components.read_text())
q8 = json.loads(args.q8_orbits.read_text())
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert zero_payload["status"] == "PASS_EXACT_CHILD_ZERO_SECTION_TRANSPORT"
assert component_payload["status"] == "PASS_EXACT_CHILD_E7_INFINITY_TRANSPORT"
assert q8["status"] == "PASS_H3_Q6_CHILD_Q8_WEYL_CLASSIFICATION"

# First recover the exact standard-Jacobian points.  ``old_zero`` is the
# chosen zero for the q=6 child fibration, but is a finite point on the
# standard short Jacobian; the component points are likewise finite points.
base_ring = PolynomialRing(QQ, "T")
T = base_ring.gen()
base_field = base_ring.fraction_field()
coefficient_a = base_field(polynomial(
    base_ring, child["minimal_short_weierstrass"]["A_coefficients_low_to_high"]
))
coefficient_b = base_field(polynomial(
    base_ring, child["minimal_short_weierstrass"]["B_coefficients_low_to_high"]
))
curve = EllipticCurve(base_field, [0, 0, 0, coefficient_a, coefficient_b])
zero_data = zero_payload["section"]
old_zero = curve(
    rational_function(base_field, base_ring,
                      zero_data["x_numerator_coefficients_low_to_high"],
                      zero_data["x_denominator_coefficients_low_to_high"]),
    rational_function(base_field, base_ring,
                      zero_data["y_numerator_coefficients_low_to_high"],
                      zero_data["y_denominator_coefficients_low_to_high"]),
)
component_points = {}
for entry in component_payload["sections"]:
    component_points[entry["sign"]] = curve(
        rational_function(base_field, base_ring,
                          entry["x_numerator_coefficients_low_to_high"],
                          entry["x_denominator_coefficients_low_to_high"]),
        rational_function(base_field, base_ring,
                          entry["y_numerator_coefficients_low_to_high"],
                          entry["y_denominator_coefficients_low_to_high"]),
    )
assert set(component_points) == {"plus", "minus"}
affine = component_points[component_payload["source"]["affine_E7_sign"]]
e7_7 = component_points[component_payload["source"]["E7_7_sign"]]
assert all(point in curve for point in (old_zero, affine, e7_7))

# Relative to the child zero, these are the two marked MW directions.  On the
# standard Jacobian they are represented by the displayed differences.
first_direction = e7_7 - old_zero
second_direction = e7_7 - affine
selected_relative_section = 2 * first_direction + 2 * second_direction
assert not first_direction.is_zero()
assert not second_direction.is_zero()
assert not selected_relative_section.is_zero()

# The generic marking alone does not identify a resolved fibre component, but
# its exact specialization gives the indispensable starting jet.  The section
# is smooth at the II* Weierstrass fibre and runs into the IV* singular point
# with the depth-two jet that must be followed through the IV* blow-up charts.
def local_order(value, factor):
    return (
        value.numerator().valuation(factor)
        - value.denominator().valuation(factor)
    )


sx, sy = selected_relative_section.xy()


def monic_power_root(value, exponent):
    root = base_ring.one()
    factorization = value.factor()
    for irreducible, multiplicity in factorization:
        assert multiplicity % exponent == 0
        root *= irreducible.monic()**(multiplicity // exponent)
    return root.monic(), tuple(
        (int(irreducible.degree()), int(multiplicity // exponent))
        for irreducible, multiplicity in factorization
    )


# The two horizontal poles in O+S collide over a finite divisor on the old
# base.  It is essential data for a global pencil: raw generic chord functions
# have base jets there, even though their restriction to a generic old fibre
# is only the desired degree-two system.
x_denominator = base_ring(sx.denominator())
y_denominator = base_ring(sy.denominator())
collision_divisor, collision_factorization = monic_power_root(x_denominator, 2)
y_collision_divisor, y_collision_factorization = monic_power_root(y_denominator, 3)
assert collision_divisor == y_collision_divisor
assert x_denominator.degree() == 2 * collision_divisor.degree()
assert y_denominator.degree() == 3 * collision_divisor.degree()
assert collision_divisor.degree() == 46
discriminant = -16 * (4 * coefficient_a**3 + 27 * coefficient_b**2)
assert base_ring(discriminant.numerator()).gcd(collision_divisor).degree() == 0
# In the standard smooth Weierstrass chart at O, z=-x/y is the formal
# parameter.  Its divisor along the collision polynomial is reduced, so the
# two horizontal sections meet transversely at every closed collision point.
zero_parameter = -sx / sy
z_numerator = base_ring(zero_parameter.numerator())
z_denominator = base_ring(zero_parameter.denominator())
z_reduced_numerator, remainder = z_numerator.quo_rem(collision_divisor)
assert not remainder
assert z_reduced_numerator.gcd(collision_divisor).degree() == 0
assert z_denominator.gcd(collision_divisor).degree() == 0
for irreducible, multiplicity in collision_divisor.factor():
    assert multiplicity == 1 and local_order(zero_parameter, irreducible) == 1
collision_data = {
    "divisor_degree": int(collision_divisor.degree()),
    "x_denominator": "constant * collision_divisor^2",
    "y_denominator": "constant * collision_divisor^3",
    "irreducible_factor_degrees_and_multiplicities": [list(item) for item in collision_factorization],
    "is_coprime_to_weierstrass_discriminant": True,
    "zero_section_formal_parameter": "z=-x/y",
    "zero_parameter_order_at_every_collision": 1,
    "intersection_with_standard_zero": "transverse at every closed collision point",
    "interpretation": "At every closed point of this smooth-base divisor, S meets the standard zero section transversely.",
}

local_specialization = {}
for fibre in child["finite_fibres"]:
    kind = fibre["kodaira"]
    if kind not in ("II*", "IV*"):
        continue
    factor = base_ring(fibre["factor"])
    assert factor.degree() == 1
    base_point = -factor[0] / factor[1]
    x_order = local_order(sx, factor)
    y_order = local_order(sy, factor)
    local_specialization[kind] = {
        "factor": str(factor),
        "x_order": int(x_order),
        "y_order": int(y_order),
    }
    if kind in ("II*", "IV*"):
        assert (x_order, y_order) == (0, 0)
        assert (sx(base_point), sy(base_point)) != (0, 0)
        local_specialization[kind]["weierstrass_specialization"] = "smooth"
assert set(local_specialization) == {"II*", "IV*"}

# Recompute the coordinate bridge from the exact H3 NS frame.  This proves
# the orientation of the two component points in the same MW basis used by
# the q=8 orbit classifier, rather than identifying them from a height alone.
frame = load_gram(FRAME)
ns = block_diagonal_matrix(matrix(ZZ, ((0, 1), (1, 0))), -frame)
source_fiber = vector(ZZ, [1, 0] + [0] * 17)
source_zero = vector(ZZ, [-1, 1] + [0] * 17)
simple = tuple(
    vector(ZZ, [0, 0] + [ZZ(index == node) for index in range(17)])
    for node in range(15)
)
raw_q6_fiber = vector(ZZ, [3, 2] + [
    0, 0, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
])


def reflect(class_value, nodes):
    result = vector(ZZ, class_value)
    for node in nodes:
        curve_component = simple[node - 1]
        pairing = result * ns * curve_component
        assert pairing in ZZ
        result += pairing * curve_component
    return result


nef_q6_fiber = reflect(raw_q6_fiber, REFLECTIONS)
assert nef_q6_fiber * ns * nef_q6_fiber == 0
assert nef_q6_fiber * ns * source_fiber == 2
inverse_transport = lambda value: reflect(value, tuple(reversed(REFLECTIONS)))
raw_zero = inverse_transport(source_zero)
raw_e7_7 = inverse_transport(simple[6])
highest_e7 = (2, 2, 3, 4, 3, 2, 1)
affine_e7 = source_fiber - sum(
    (coefficient * simple[index] for index, coefficient in enumerate(highest_e7)),
    vector(ZZ, [0] * 19),
)
raw_affine = inverse_transport(affine_e7)
raw_mate = isotropic_mate(ns, raw_q6_fiber)
raw_orthogonal = matrix(
    ZZ, [list(raw_q6_fiber * ns), list(raw_mate * ns)]
).right_kernel_matrix()
raw_child = -(raw_orthogonal * ns * raw_orthogonal.transpose())
raw_roots = matrix(ZZ, pari(raw_child).qfminim(2)[2]).transpose().row_module().basis_matrix()
assert raw_roots.rank() == 14
root_curves = raw_roots * raw_orthogonal
root_gram = root_curves * ns * root_curves.transpose()
projection = identity_matrix(QQ, 19) - ns * root_curves.transpose() * root_gram.inverse() * root_curves


def shioda(section):
    horizontal = section - raw_zero - (section * ns * raw_zero + 2) * raw_q6_fiber
    assert horizontal * ns * raw_q6_fiber == horizontal * ns * raw_zero == 0
    return vector(QQ, horizontal) * projection


mw_sections = []
for lift in MW_LIFTS.rows():
    start = raw_zero + vector(ZZ, lift) * raw_orthogonal
    section = start + ((-2 - start * ns * start) // 2) * raw_q6_fiber
    assert section * ns * section == -2 and section * ns * raw_q6_fiber == 1
    mw_sections.append(section)
mw_shioda = [shioda(section) for section in mw_sections]
height = matrix(QQ, [[-left * ns * right for right in mw_shioda] for left in mw_shioda])
assert height == EXPECTED_HEIGHT


def mw_coordinates(value):
    rhs = vector(QQ, [-value * ns * basis for basis in mw_shioda])
    coordinates = rhs * height.inverse()
    assert all(entry in ZZ for entry in coordinates)
    return vector(ZZ, coordinates)


first_coordinates = mw_coordinates(shioda(raw_e7_7))
second_coordinates = mw_coordinates(shioda(raw_e7_7) - shioda(raw_affine))
assert first_coordinates == vector(ZZ, (-1, 0, 0))
assert second_coordinates == vector(ZZ, (0, -1, 0))

selected_hit = next(
    hit for hit in q8["q8"]["d13_mw4_hits"]
    if hit["mw_projection"] == [0, -2, 0]
)
assert q8["q8"]["old_fiber_degree"] == 2
selected_coordinates = vector(ZZ, selected_hit["mw_projection"])
assert selected_coordinates == vector(ZZ, (0, -2, 0))
canonical_zero = raw_mate - raw_q6_fiber
assert canonical_zero * ns * canonical_zero == -2
assert canonical_zero * ns * raw_q6_fiber == 1
canonical_zero_coordinates = mw_coordinates(shioda(canonical_zero))
assert canonical_zero_coordinates == vector(ZZ, (-1, 0, 0))
selected_relative_coordinates = selected_coordinates + 2 * canonical_zero_coordinates
assert selected_relative_coordinates == vector(ZZ, (-2, -2, 0))
assert selected_relative_coordinates == 2 * first_coordinates + 2 * second_coordinates
selected_q8_class = vector(ZZ, selected_hit["fiber_source_h3_ns"])
assert selected_q8_class * ns * raw_q6_fiber == 2
selected_q8_horizontal = selected_q8_class - 2 * raw_zero
selected_q8_horizontal -= (selected_q8_horizontal * ns * raw_zero) * raw_q6_fiber
assert selected_q8_horizontal * ns * raw_q6_fiber == 0
assert selected_q8_horizontal * ns * raw_zero == 0
selected_q8_projection = vector(QQ, selected_q8_horizontal) * projection
expected_q8_projection = sum(
    (coordinate * basis for coordinate, basis in zip(selected_relative_coordinates, mw_shioda)),
    vector(QQ, [0] * 19),
)
assert mw_coordinates(selected_q8_projection) == selected_relative_coordinates
assert selected_q8_projection - expected_q8_projection == -19 * raw_q6_fiber

# The corrected section is smooth at both additive Weierstrass fibres, hence
# meets the same identity components as the transported old zero.  Its exact
# NS class is therefore recovered from its Shioda projection by the unique
# fibre correction giving square -2.  Subtracting this horizontal divisor
# from the selected q8 fibre produces the concrete vertical target needed by
# future II*/IV* local-module compilers.
selected_section_base = raw_zero + expected_q8_projection
selected_section_fibre_correction = (-2 - selected_section_base * ns * selected_section_base) / 2
assert selected_section_fibre_correction in ZZ
selected_section_ns = selected_section_base + selected_section_fibre_correction * raw_q6_fiber
assert selected_section_ns * ns * selected_section_ns == -2
assert selected_section_ns * ns * raw_q6_fiber == 1
assert selected_section_ns * ns * root_curves.transpose() == raw_zero * ns * root_curves.transpose()
selected_q8_vertical = vector(QQ, selected_q8_class) - vector(QQ, raw_zero) - selected_section_ns
vertical_basis = matrix(QQ, [list(raw_q6_fiber)] + [list(row) for row in root_curves.rows()])
vertical_coordinates = vertical_basis.solve_left(selected_q8_vertical)
assert all(value in ZZ for value in vertical_coordinates)
assert selected_q8_vertical * ns * raw_q6_fiber == 0
assert vector(QQ, selected_q8_class) == vector(QQ, raw_zero) + selected_section_ns + selected_q8_vertical

# In standard-Jacobian coordinates the translated divisor O_child+S becomes
# O_standard + 2*first_direction + 2*second_direction.  The slope through -S
# has poles at the standard zero and S, hence gives the exact two-dimensional
# generic space L(O_standard+S) with basis (1,m).  This is only the generic
# RR ambient.
chord_ring = PolynomialRing(base_field, "m")
m = chord_ring.gen()
relation_ring = PolynomialRing(chord_ring, "x")
x = relation_ring.gen()
relation = (
    (relation_ring(-sy) + relation_ring(m) * (x - relation_ring(sx)))**2
    - x**3 - relation_ring(coefficient_a) * x - relation_ring(coefficient_b)
)
quadratic, remainder = relation.quo_rem(x - relation_ring(sx))
assert not remainder
assert quadratic.degree(x) == 2 and quadratic.leading_coefficient() == -1
quadratic = -quadratic


def point_data(point):
    x_value, y_value = point.xy()
    return {
        "x_numerator_coefficients_low_to_high": [str(value) for value in base_ring(x_value.numerator()).list()],
        "x_denominator_coefficients_low_to_high": [str(value) for value in base_ring(x_value.denominator()).list()],
        "y_numerator_coefficients_low_to_high": [str(value) for value in base_ring(y_value.numerator()).list()],
        "y_denominator_coefficients_low_to_high": [str(value) for value in base_ring(y_value.denominator()).list()],
        "coordinate_degrees": {
            "x": [int(x_value.numerator().degree()), int(x_value.denominator().degree())],
            "y": [int(y_value.numerator().degree()), int(y_value.denominator().degree())],
        },
    }


payload = {
    "schema": "elkies-k3.h92-q6-child-q8-marking.v1",
    "status": "PASS_EXACT_Q6_CHILD_Q8_MARKING",
    "inputs": {
        "frame": {"path": str(FRAME.relative_to(ROOT)), "sha256": digest(FRAME)},
        "child_jacobian": {"path": str(args.child.relative_to(ROOT)), "sha256": digest(args.child)},
        "child_zero": {"path": str(args.zero.relative_to(ROOT)), "sha256": digest(args.zero)},
        "child_e7_components": {"path": str(args.components.relative_to(ROOT)), "sha256": digest(args.components)},
        "q8_orbits": {"path": str(args.q8_orbits.relative_to(ROOT)), "sha256": digest(args.q8_orbits)},
    },
    "component_orientation": {
        "affine_E7_standard_jacobian_sign": component_payload["source"]["affine_E7_sign"],
        "E7_7_standard_jacobian_sign": component_payload["source"]["E7_7_sign"],
        "E7_7_minus_child_zero_MW_coordinates": list(map(int, first_coordinates)),
        "E7_7_minus_affine_E7_MW_coordinates": list(map(int, second_coordinates)),
        "canonical_q6_zero_minus_transported_old_zero_MW_coordinates": list(
            map(int, canonical_zero_coordinates)
        ),
        "height_gram": [[str(value) for value in row] for row in height.rows()],
    },
    "selected_q8": {
        "orbit_mw_projection_at_canonical_q6_zero": list(map(int, selected_coordinates)),
        "relative_child_section_MW_coordinates": list(map(int, selected_relative_coordinates)),
        "zero_translation": (
            "The q8 orbit coordinate is based at the canonical q6 zero; the "
            "explicit child uses the transported old zero, whose relative "
            "coordinate is (1,0,0)."
        ),
        "relative_child_section": "2*(E7_7-old_zero) + 2*(E7_7-affine_E7)",
        "relative_child_section_standard_jacobian_coordinates": point_data(selected_relative_section),
        "ns_horizontal_vertical_decomposition": {
            "identity_component_rule": "The corrected section is smooth at both additive Weierstrass fibres, hence has the same root-component intersections as the transported old zero.",
            "section_class_in_source_h3_ns": list(map(int, selected_section_ns)),
            "section_fibre_correction": int(selected_section_fibre_correction),
            "root_basis_in_source_h3_ns": [list(map(int, row)) for row in root_curves.rows()],
            "vertical_coordinates_in_fibre_then_root_basis": list(map(int, vertical_coordinates)),
            "identity": "q8_fibre=old_zero+corrected_marked_section+vertical",
        },
        "generic_divisor_after_translation": "O_standard + 2*(E7_7-old_zero) + 2*(E7_7-affine_E7)",
        "generic_rr_basis": ["1", "m"],
        "m": "(y+y(S))/(x-x(S)), S=2*(E7_7-old_zero)+2*(E7_7-affine_E7)",
        "quadratic_relation_for_x_over_QQ(T)(m)": str(quadratic),
        "zero_section_collision_divisor": collision_data,
        "additive_fibre_specialization": local_specialization,
    },
    "boundary": (
        "This pins the q8 generic-fibre marking and its exact section coordinate. "
        "The recorded additive and smooth-fibre collision jets still do not derive the II*/IV* local modules, a two-dimensional global q8 pencil, "
        "the D13 equation, a rootless bisection, an extension collision, or generic rank 18 or 19."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6CHILDQ8|canonical_mw=0,-2,0|relative_mw=-2,-2,0|"
    "generic_basis=2|zero_collisions=46:transverse|II*=smooth|IV*=smooth|"
    "status=PASS_EXACT_Q6_CHILD_Q8_MARKING",
    flush=True,
)
