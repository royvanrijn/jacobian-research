#!/usr/bin/env sage-python
"""Compile NS0024 q4/orbit1 from a certified finite-field MW4 family.

status: ACTIVE_COMPILER
claim: exact modular resolved RR plane, binary quartic, and target Jacobian

The input must use the schema
``elkies-k3.lattice-foundry-ns0024-mw4-family-modp.v1`` and bind its four
displayed sections to the exact minimum-pole source marking.  The compiler is
fail-closed: it replays the source fibre profile and section equations, checks
the P3 split-I5 orientation, constructs the complete 4-to-2 resolved toric
module for ``-C2-2*C3-C4``, and compiles the same kernel to the child
Jacobian.  Success requires the exact geometric fibre profile
``I1*+I5+I3+I2+7I1`` (root system ``D5+A4+A2+A1``).

This is a finite-field equation certificate.  It does not lift the family to
characteristic zero or prove a Picard number there.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PREPARATION = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-edge1-compiler-preparation.json"
BASIS = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-mw4-minimum-basis.json"

load(str(HERE / "elliptic_neighbor_compiler.sage"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def display_path(path):
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def jsonable(value):
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [jsonable(item) for item in value]
    try:
        return int(value) if value in ZZ else str(value)
    except (TypeError, ValueError):
        return str(value)


def polynomial_payload(poly):
    return [str(value) for value in poly.list()]


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()

input_path = args.input.resolve()
output_path = args.output.resolve()
model = json.loads(input_path.read_text())
preparation = json.loads(PREPARATION.read_text())
basis = json.loads(BASIS.read_text())
assert preparation["status"] == "PASS_EXACT_NS0024_EDGE1_COMPILER_PREPARATION"
assert basis["status"] == "PASS_EXACT_MINIMUM_POLE_FOUR_SECTION_BASIS"
if model.get("schema") != "elkies-k3.lattice-foundry-ns0024-mw4-family-modp.v1":
    raise ValueError("input is not an NS0024 MW4 modular family")
if model.get("status") != "PASS_EXACT_MODULAR_NS0024_MW4_FAMILY_MARKING":
    raise ValueError("input MW4 family marking is not certified")

prime = ZZ(model["prime"])
if not prime.is_prime() or prime in (2, 3, 5, 7):
    raise ValueError("edge-1 compiler requires a good prime outside 2,3,5,7")
prime_field = GF(prime)
extension_record = model.get("extension")
locals_map = {}
if extension_record is None:
    constant_field = prime_field
else:
    extension_generator_name = extension_record.get("generator")
    if (
        not isinstance(extension_generator_name, str)
        or not extension_generator_name.isidentifier()
        or extension_generator_name in ("t", "U")
    ):
        raise ValueError("finite-extension generator must be an identifier distinct from t,U")
    extension_polynomial_ring = PolynomialRing(prime_field, extension_generator_name)
    extension_indeterminate = extension_polynomial_ring.gen()
    extension_modulus = extension_polynomial_ring(
        sage_eval(
            str(extension_record.get("modulus")),
            locals={extension_generator_name: extension_indeterminate},
        )
    )
    if extension_modulus.degree() < 2 or not extension_modulus.is_irreducible():
        raise ValueError("finite-extension modulus must be irreducible of degree at least two")
    extension_modulus = extension_modulus.monic()
    constant_field = GF(
        prime ** extension_modulus.degree(),
        extension_generator_name,
        modulus=extension_modulus,
    )
    locals_map[extension_generator_name] = constant_field.gen()
parameter_names = tuple(model.get("parameters", ()))
if len(parameter_names) > 1:
    raise ValueError("edge-1 compiler currently supports at most one family parameter")
if parameter_names and parameter_names[0] in ("t", "U", *locals_map.keys()):
    raise ValueError("family parameter conflicts with an elliptic base coordinate")
if parameter_names:
    parameter_ring = PolynomialRing(constant_field, parameter_names[0])
    family_parameter = parameter_ring.gen()
    coefficient_field = parameter_ring.fraction_field()
    locals_map[parameter_names[0]] = coefficient_field(family_parameter)
else:
    coefficient_field = constant_field

old_ring = PolynomialRing(coefficient_field, "t")
t = old_ring.gen()
old_field = old_ring.fraction_field()
locals_map["t"] = old_field(t)


def parse_coefficient(value):
    if isinstance(value, int):
        return coefficient_field(value)
    return coefficient_field(sage_eval(str(value), locals=locals_map))


def parse_function(value):
    if isinstance(value, int):
        return old_field(value)
    return old_field(sage_eval(str(value), locals=locals_map))


surface = model["surface"]
A = old_ring([parse_coefficient(value) for value in surface["A_coefficients_low_to_high"]])
B = old_ring([parse_coefficient(value) for value in surface["B_coefficients_low_to_high"]])
if A.degree() > 8 or B.degree() > 12:
    raise ValueError("source model exceeds elliptic-K3 short-Weierstrass degree bounds")
discriminant = old_ring(-16 * (4 * A**3 + 27 * B**2))
if not discriminant:
    raise ValueError("source family is singular over its generic parameter field")

marking = model["marking"]
expected_basis_hash = digest(BASIS)
if marking.get("minimum_basis_sha256") != expected_basis_hash:
    raise ValueError("input marking does not bind the pinned NS0024 MW4 basis")
if marking.get("normalized_supports") != {"I7": "0", "I5": "1", "I4": "infinity"}:
    raise ValueError("input family uses a different reducible-fibre normalization")
expected_profiles = {
    item["name"]: item["components_I7_I5_I4"] for item in basis["basis"]
}
if marking.get("section_profiles_I7_I5_I4") != expected_profiles:
    raise ValueError("input section profiles do not match the pinned MW4 marking")
expected_intersection_gram = basis["section_intersection_gram"]
if marking.get("section_intersection_gram") != expected_intersection_gram:
    raise ValueError("input marking does not declare the pinned section intersection Gram matrix")

order_zero = int(discriminant.valuation(t))
order_one = int(discriminant.valuation(t - 1))
order_infinity = int(24 - discriminant.degree())
if (order_zero, order_one, order_infinity) != (7, 5, 4):
    raise ValueError(
        "source reducible orders are {}, expected (7,5,4)".format(
            (order_zero, order_one, order_infinity)
        )
    )
residual = old_ring(discriminant // (t**7 * (t - 1)**5))
if residual.degree() != 8 or residual(0) == 0 or residual(1) == 0:
    raise ValueError("source residual discriminant does not have eight separated degrees")
if residual.gcd(residual.derivative()).degree() != 0:
    raise ValueError("source residual discriminant is not squarefree")

section_records = model["sections"]
if set(section_records) != {"P1", "P2", "P3", "P4"}:
    raise ValueError("input must display exactly the pinned P1,P2,P3,P4 basis")
sections = {
    name: (parse_function(record["x"]), parse_function(record["y"]))
    for name, record in section_records.items()
}
for name, (x_value, y_value) in sections.items():
    if y_value**2 != x_value**3 + old_field(A) * x_value + old_field(B):
        raise ValueError("{} misses the source Weierstrass family".format(name))
    x_den = old_ring(x_value.denominator()).monic()
    y_den = old_ring(y_value.denominator()).monic()
    expected_pole = next(item["P_dot_O"] for item in basis["basis"] if item["name"] == name)
    if expected_pole == 0:
        if x_den.degree() or y_den.degree():
            raise ValueError("{} is not polynomial despite P.O=0".format(name))
    elif expected_pole == 1:
        if (x_den.degree(), y_den.degree()) != (2, 3) or x_den**3 != y_den**2:
            raise ValueError("{} does not have one common simple pole".format(name))
    else:
        raise ArithmeticError("unexpected pinned pole order")


def node_x(support):
    a_value = coefficient_field(A(support))
    b_value = coefficient_field(B(support))
    if not a_value:
        raise ValueError("multiplicative support has c4=0")
    value = coefficient_field(-3 * b_value / (2 * a_value))
    if value**3 + a_value * value + b_value or 3 * value**2 + a_value:
        raise ValueError("declared multiplicative support is not nodal")
    return value


def hits_finite_node(point, support):
    x_value, y_value = point
    if x_value.denominator()(support) == 0 or y_value.denominator()(support) == 0:
        raise ValueError("section pole collides with a normalized reducible support")
    return bool(x_value(support) == node_x(support) and y_value(support) == 0)


expected_node_hits = {
    "P1": (True, False),
    "P2": (True, True),
    "P3": (True, True),
    "P4": (True, True),
}
for name, point in sections.items():
    actual = (hits_finite_node(point, coefficient_field(0)), hits_finite_node(point, coefficient_field(1)))
    if actual != expected_node_hits[name]:
        raise ValueError("{} has the wrong I7/I5 node-incidence fingerprint".format(name))

# Replay the complete component and intersection marking rather than trusting
# the input labels.  P4 has component 1 in all three normalized fibres, hence
# it is a uniform local component-group reference.  A section Q has label j
# in I_n exactly when Q-j*P4 specializes to the identity component, i.e. no
# longer passes through the node in the displayed Weierstrass chart.
source_curve = EllipticCurve(old_field, [0, 0, 0, A, B])
section_points = {
    name: source_curve(old_field(x_value), old_field(y_value))
    for name, (x_value, y_value) in sections.items()
}


def finite_value(value, support):
    numerator, denominator = value.numerator(), value.denominator()
    return None if denominator(support) == 0 else numerator(support) / denominator(support)


def hits_normalized_node(point, fibre_index):
    if point.is_zero():
        return False
    if fibre_index < 2:
        support = coefficient_field(fibre_index)
        return (
            finite_value(point[0], support) == node_x(support)
            and finite_value(point[1], support) == 0
        )
    infinity_a = coefficient_field(A[8])
    infinity_b = coefficient_field(B[12])
    if not infinity_a:
        raise ValueError("normalized I4 fibre at infinity has c4=0")
    infinity_node = coefficient_field(-3 * infinity_b / (2 * infinity_a))
    if (
        infinity_node**3 + infinity_a * infinity_node + infinity_b
        or 3 * infinity_node**2 + infinity_a
    ):
        raise ValueError("normalized infinity support is not nodal")
    x_num, x_den = point[0].numerator(), point[0].denominator()
    y_num, y_den = point[1].numerator(), point[1].denominator()
    x_excess = x_num.degree() - x_den.degree()
    y_excess = y_num.degree() - y_den.degree()
    x_value = (
        coefficient_field.zero()
        if x_excess < 4
        else x_num.leading_coefficient() / x_den.leading_coefficient()
        if x_excess == 4
        else None
    )
    y_value = (
        coefficient_field.zero()
        if y_excess < 6
        else y_num.leading_coefficient() / y_den.leading_coefficient()
        if y_excess == 6
        else None
    )
    return x_value == infinity_node and y_value == 0


def component_label(point, reference, order, fibre_index):
    labels = [
        multiplier
        for multiplier in range(order)
        if not hits_normalized_node(point - multiplier * reference, fibre_index)
    ]
    return labels[0] if len(labels) == 1 else -1


p4_reference = section_points["P4"]
actual_profiles = {
    name: [
        component_label(point, p4_reference, order, fibre_index)
        for fibre_index, order in enumerate((7, 5, 4))
    ]
    for name, point in section_points.items()
}
if actual_profiles != expected_profiles:
    raise ValueError("displayed sections have the wrong complete component marking: {}".format(actual_profiles))


def section_intersection(left, right):
    difference = left - right
    if difference.is_zero():
        return -2
    numerator, denominator = difference[0].numerator(), difference[0].denominator()
    degree = denominator.degree() + max(
        0, numerator.degree() - denominator.degree() - 4
    )
    if degree % 2:
        raise ArithmeticError("section intersection degree is not even")
    return degree // 2


ordered_points = [section_points[name] for name in ("P1", "P2", "P3", "P4")]
actual_intersection_gram = [
    [section_intersection(left, right) for right in ordered_points]
    for left in ordered_points
]
if actual_intersection_gram != expected_intersection_gram:
    raise ValueError("displayed sections have the wrong intersection Gram matrix")

P3x, P3y = sections["P3"]
ambient = ("1", "t", "t^2", "m_P3")
chord_expansions = ((1, 0), (t, 0), (t**2, 0), (0, 1))
toric_block = split_multiplicative_toric_chord_condition(
    "NS0024 edge1 split-I5 C2+2C3+C4 quotient",
    ambient,
    chord_expansions,
    old_ring,
    coefficient_field(1),
    A,
    B,
    P3x,
    P3y,
    5,
    1,
    {2: 1, 3: 2, 4: 1},
    (
        "exact q4/orbit1 divisor D=O+P3+2F-C2-2C3-C4 from the pinned "
        "NS0024 source marking"
    ),
)
compiled = compile_resolved_conditions(
    ambient, (toric_block,), complete=True, coefficient_field=coefficient_field
)
if (
    compiled["ambient_dimension"], compiled["rank"], compiled["kernel_dimension"]
) != (4, 2, 2) or not compiled["h0_certified"]:
    raise ArithmeticError("NS0024 edge1 resolved RR dimension gate failed")
pencil_basis = compiled["kernel_basis"]

new_ring = PolynomialRing(coefficient_field, "U")
U_poly = new_ring.gen()
new_field = new_ring.fraction_field()
old_over_new_ring = PolynomialRing(new_field, "t")
tt = old_over_new_ring.gen()


def lift_polynomial(poly):
    return old_over_new_ring([new_field(value) for value in old_ring(poly).list()])


resolved_hop = compile_resolved_degree_two_chord_hop(
    compiled,
    pencil_basis,
    ((1, 0), (tt, 0), (tt**2, 0), (0, 1)),
    old_over_new_ring,
    new_field(U_poly),
    lift_polynomial(P3x.numerator()) / lift_polynomial(P3x.denominator()),
    lift_polynomial(P3y.numerator()) / lift_polynomial(P3y.denominator()),
    lift_polynomial(A),
    old_b=lift_polynomial(B),
)
conversion = resolved_hop["conversion"]
quartic = conversion["binary_quartic"]
if quartic.degree() != 4:
    raise ArithmeticError("edge1 resolved pencil did not compile to a binary quartic")
child_a = new_field(conversion["jacobian_a"])
child_b = new_field(conversion["jacobian_b"])
classification = classify_finite_short_weierstrass_fibres(new_ring, child_a, child_b)

finite_fibres = classification["finite_fibres"]
root_rank = ZZ(classification["finite_root_rank"])
root_determinant = ZZ(classification["finite_root_determinant"])
euler_number = ZZ(classification["finite_euler_number"])
root_symbols = []
for fibre in finite_fibres:
    if fibre["root_rank"]:
        root_symbols.extend([fibre["kodaira"]] * fibre["degree"])

infinity_record = classification["infinity_boundary"]
infinity_orders = tuple(map(ZZ, infinity_record["normalized_orders"]))
infinity_fibre = None
if infinity_orders[2] > 0:
    inf_rank, inf_euler, inf_det, inf_symbol = kodaira_data_from_short_orders(
        *infinity_orders
    )
    infinity_fibre = {
        "kodaira": inf_symbol,
        "orders": tuple(map(int, infinity_orders)),
        "root_rank": inf_rank,
        "euler_number": inf_euler,
        "root_determinant": inf_det,
    }
    root_rank += inf_rank
    root_determinant *= inf_det
    euler_number += inf_euler
    if inf_rank:
        root_symbols.append(inf_symbol)

if sorted(root_symbols) != sorted(("I1*", "I2", "I3", "I5")):
    raise ArithmeticError("edge1 child has wrong reducible fibres: {}".format(root_symbols))
if (root_rank, root_determinant, euler_number) != (12, 120, 24):
    raise ArithmeticError(
        "edge1 child target totals are {}, expected (12,120,24)".format(
            (root_rank, root_determinant, euler_number)
        )
    )

minimal = classification["finite_minimization"]
payload = {
    "schema": "elkies-k3.lattice-foundry-ns0024-edge1-compilation-modp.v1",
    "status": "PASS_EXACT_MODULAR_NS0024_EDGE1_RESOLVED_RR_AND_JACOBIAN",
    "prime": int(prime),
    "extension": extension_record,
    "parameters": list(parameter_names),
    "source": {
        "fibre_profile": "I7+I5+I4+8I1",
        "discriminant_orders_0_1_infinity": [order_zero, order_one, order_infinity],
        "residual_discriminant_degree": int(residual.degree()),
        "residual_discriminant_squarefree": True,
        "displayed_sections": sorted(sections),
        "section_profiles_I7_I5_I4": actual_profiles,
        "section_intersection_gram": actual_intersection_gram,
        "horizontal": "P3",
        "horizontal_component_I5": 1,
    },
    "edge": {
        "q": 4,
        "orbit_index": 1,
        "old_fibre_degree": 2,
        "divisor_identity": preparation["edge"]["divisor_identity"],
    },
    "resolved_RR": {
        "ambient_basis": list(ambient),
        "condition_matrix": [
            [str(value) for value in row] for row in compiled["condition_matrix"].rows()
        ],
        "condition_rank": int(compiled["rank"]),
        "kernel_dimension": int(compiled["kernel_dimension"]),
        "h0": 2,
        "kernel_basis": [
            [str(value) for value in row] for row in pencil_basis.rows()
        ],
        "toric_diagnostics": jsonable(toric_block["resolved_toric_diagnostics"]),
        "complete_resolved_chart_cover": True,
    },
    "quartic": {
        "degree": int(quartic.degree()),
        "coefficients_low_to_high": [str(value) for value in quartic.list()],
        "certified_square_factor": str(conversion["square_factor"]),
    },
    "child": {
        "raw_jacobian_A": str(child_a),
        "raw_jacobian_B": str(child_b),
        "finite_minimal_A_coefficients_low_to_high": polynomial_payload(minimal["minimal_a"]),
        "finite_minimal_B_coefficients_low_to_high": polynomial_payload(minimal["minimal_b"]),
        "finite_fibres": [
            {
                "factor": str(fibre["factor"]),
                "degree": fibre["degree"],
                "kodaira": fibre["kodaira"],
                "minimal_orders": list(fibre["minimal_orders"]),
                "root_rank": fibre["root_rank"],
                "root_determinant": fibre["root_determinant"],
                "euler_number": fibre["euler_number"],
            }
            for fibre in finite_fibres
        ],
        "infinity_fibre": infinity_fibre,
        "root_type": "A1+A2+A4+D5",
        "root_rank": int(root_rank),
        "root_determinant": int(root_determinant),
        "euler_number": int(euler_number),
        "lattice_predicted_mw_rank_for_rho_19": 5,
    },
    "inputs": {
        "paths": [
            display_path(input_path),
            display_path(PREPARATION),
            display_path(BASIS),
        ],
        "sha256": {
            display_path(input_path): digest(input_path),
            display_path(PREPARATION): digest(PREPARATION),
            display_path(BASIS): digest(BASIS),
        },
    },
    "proof_boundary": {
        "proved": (
            "For the supplied certified modular MW4 family, the exact q4/orbit1 "
            "resolved RR space is the displayed two-plane and its binary-quartic "
            "Jacobian has geometric fibre roots A1+A2+A4+D5."
        ),
        "not_proved": (
            "No characteristic-zero family, Picard-rank certificate, effective child "
            "zero, later route edge, or Mordell-Weil rank over characteristic zero is asserted."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/compile_lattice_foundry_ns0024_edge1_modp.sage "
        "--input {} --output {}"
    ).format(display_path(input_path), display_path(output_path)),
}

serialized = json.dumps(jsonable(payload), indent=2, sort_keys=True) + "\n"
if args.check:
    if output_path.read_text() != serialized:
        raise SystemExit("NS0024 edge-1 modular compilation artifact is stale")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)

print(
    "NS0024EDGE1|p={}|ambient=4|rank=2|h0=2|quartic=4|"
    "root=A1+A2+A4+D5|root_rank=12|root_det=120|euler=24|status=PASS".format(prime),
    flush=True,
)
