#!/usr/bin/env sage-python
"""Adapt an exact marked NS0024 point over GF(p^d) to the edge-1 compiler.

status: ACTIVE_COMPILER
claim: exact finite-extension source-marking replay and lossless compiler handoff

The recovery pipeline currently emits a compact P4 residue-algebra point plus
an oriented MW3 seed.  This adapter joins those two records, independently
replays the four section equations, fibre profile, component labels, and
intersection Gram matrix, and only then emits the certified source format
accepted by ``compile_lattice_foundry_ns0024_edge1_modp.sage``.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, ZZ, sage_eval


ROOT = Path(__file__).resolve().parents[2]
BASIS = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-mw4-minimum-basis.json"
EXPECTED_PROFILES = ((1, 0, 0), (2, 1, 3), (2, 1, 1), (1, 1, 1))
EXPECTED_INTERSECTIONS = (
    (-2, 1, 2, 1),
    (1, -2, 0, 1),
    (2, 0, -2, 1),
    (1, 1, 1, -2),
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def display_path(path):
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def parse_seed(path, prime_field):
    fields = {}
    for item in path.read_text().strip().split("|")[1:]:
        key, value = item.split("=", 1)
        fields[key] = value
    if fields.get("p") is None:
        raise ValueError("MW3 seed has no prime")

    def values(key):
        return [prime_field(int(value)) for value in fields[key].split(",")]

    return fields, values


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--point", type=Path, required=True)
parser.add_argument("--seed", type=Path, help="override the seed path recorded by the point")
parser.add_argument("--output", type=Path, required=True)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()

point_path = args.point.resolve()
output_path = args.output.resolve()
point = json.loads(point_path.read_text())
if point.get("schema") != "elkies-k3.lattice-foundry-ns0024-mw4-point-modp.v1":
    raise ValueError("input is not an NS0024 MW4 residue-algebra point")
if point.get("status") != "PASS_EXACT_MW4_MARKED_POINT_OVER_QUADRATIC_EXTENSION":
    raise ValueError("input residue-algebra point does not carry an exact marking certificate")

seed_path = (args.seed or Path(point["mw3_seed"])).resolve()
prime = ZZ(point["prime"])
if not prime.is_prime() or prime in (2, 3, 5, 7):
    raise ValueError("adapter requires a good prime outside 2,3,5,7")
prime_field = GF(prime)
seed_fields, seed_values = parse_seed(seed_path, prime_field)
if ZZ(seed_fields["p"]) != prime:
    raise ValueError("residue-algebra point and MW3 seed use different primes")

extension_record = point["extension"]
generator_name = extension_record["generator"]
extension_polynomial_ring = PolynomialRing(prime_field, generator_name)
extension_indeterminate = extension_polynomial_ring.gen()
modulus = extension_polynomial_ring(
    sage_eval(
        extension_record["modulus"],
        locals={generator_name: extension_indeterminate},
    )
).monic()
if modulus.degree() != 2 or not modulus.is_irreducible():
    raise ValueError("quadratic-point extension modulus must be irreducible of degree two")
constant_field = GF(prime ** modulus.degree(), generator_name, modulus=modulus)
generator = constant_field.gen()


def decode(entry):
    if not isinstance(entry, list) or len(entry) != modulus.degree():
        raise ValueError("extension coefficient has the wrong coordinate length")
    return sum(constant_field(value) * generator**index for index, value in enumerate(entry))


old_ring = PolynomialRing(constant_field, "t")
t = old_ring.gen()
old_field = old_ring.fraction_field()
A = old_ring(seed_values("A"))
B = old_ring(seed_values("B"))
curve = EllipticCurve(old_field, [0, 0, 0, A, B])


def seed_point(index):
    return curve(old_field(old_ring(seed_values(f"P{index}X"))), old_field(old_ring(seed_values(f"P{index}Y"))))


first_three = [seed_point(index) for index in (1, 2, 3)]
X4 = old_ring([decode(value) for value in point["P4"]["X_coefficients_low_to_high"]])
Y4 = old_ring([decode(value) for value in point["P4"]["Y_coefficients_low_to_high"]])
H4 = old_ring([decode(value) for value in point["P4"]["H_coefficients_low_to_high"]])
if H4.degree() != 1 or not H4.is_monic() or H4(0) == 0 or H4(1) == 0:
    raise ValueError("P4 pole divisor is not a single normalized point away from 0,1")
fourth = curve(old_field(X4 / H4**2), old_field(Y4 / H4**3))
points = first_three + [fourth]

discriminant = old_ring(-16 * (4 * A**3 + 27 * B**2))
orders = (
    discriminant.valuation(t),
    discriminant.valuation(t - 1),
    24 - discriminant.degree(),
)
if tuple(map(int, orders)) != (7, 5, 4):
    raise ValueError("MW3 seed has the wrong I7/I5/I4 discriminant orders")
residual = old_ring(discriminant // (t**7 * (t - 1)**5))
if residual.degree() != 8 or residual.gcd(residual.derivative()).degree() != 0:
    raise ValueError("MW3 seed residual discriminant is not eight separated I1 fibres")

r1 = constant_field(int(seed_fields["r1"]))
ri = constant_field(int(seed_fields["ri"]))


def finite_value(value, support):
    numerator, denominator = value.numerator(), value.denominator()
    return None if denominator(support) == 0 else numerator(support) / denominator(support)


def hits_node(point_value, fibre):
    if point_value.is_zero():
        return False
    if fibre < 2:
        x_value = finite_value(point_value[0], constant_field(fibre))
        y_value = finite_value(point_value[1], constant_field(fibre))
        node = r1 if fibre else constant_field.one()
        return x_value == node and y_value == 0
    x_num, x_den = point_value[0].numerator(), point_value[0].denominator()
    y_num, y_den = point_value[1].numerator(), point_value[1].denominator()
    x_excess = x_num.degree() - x_den.degree()
    y_excess = y_num.degree() - y_den.degree()
    x_value = (
        constant_field.zero()
        if x_excess < 4
        else x_num.leading_coefficient() / x_den.leading_coefficient()
        if x_excess == 4
        else None
    )
    y_value = (
        constant_field.zero()
        if y_excess < 6
        else y_num.leading_coefficient() / y_den.leading_coefficient()
        if y_excess == 6
        else None
    )
    return x_value == ri and y_value == 0


def component_label(point_value, reference, order, fibre):
    labels = [
        multiplier
        for multiplier in range(order)
        if not hits_node(point_value - multiplier * reference, fibre)
    ]
    return labels[0] if len(labels) == 1 else -1


def intersection(left, right):
    difference = left - right
    if difference.is_zero():
        return -2
    numerator, denominator = difference[0].numerator(), difference[0].denominator()
    degree = denominator.degree() + max(0, numerator.degree() - denominator.degree() - 4)
    if degree % 2:
        raise ArithmeticError("section intersection degree is not even")
    return degree // 2


profiles = tuple(
    tuple(component_label(point_value, fourth, order, fibre) for fibre, order in enumerate((7, 5, 4)))
    for point_value in points
)
intersections = tuple(tuple(intersection(left, right) for right in points) for left in points)
if profiles != EXPECTED_PROFILES:
    raise ValueError("joined point has the wrong component marking: {}".format(profiles))
if intersections != EXPECTED_INTERSECTIONS:
    raise ValueError("joined point has the wrong section intersection Gram matrix")
if tuple(tuple(item) for item in point["component_profiles_I7_I5_I4"]) != EXPECTED_PROFILES:
    raise ValueError("point metadata disagrees with the replayed component marking")
if tuple(tuple(item) for item in point["section_intersection_gram"]) != EXPECTED_INTERSECTIONS:
    raise ValueError("point metadata disagrees with the replayed intersection matrix")

basis = json.loads(BASIS.read_text())
if basis.get("status") != "PASS_EXACT_MINIMUM_POLE_FOUR_SECTION_BASIS":
    raise ValueError("pinned abstract MW4 basis is not certified")
profile_record = {
    item["name"]: item["components_I7_I5_I4"] for item in basis["basis"]
}

payload = {
    "schema": "elkies-k3.lattice-foundry-ns0024-mw4-family-modp.v1",
    "status": "PASS_EXACT_MODULAR_NS0024_MW4_FAMILY_MARKING",
    "prime": int(prime),
    "extension": {"generator": generator_name, "modulus": str(modulus)},
    "parameters": [],
    "surface": {
        "A_coefficients_low_to_high": [str(value) for value in A.list()],
        "B_coefficients_low_to_high": [str(value) for value in B.list()],
    },
    "sections": {
        "P1": {"x": str(first_three[0][0]), "y": str(first_three[0][1])},
        "P2": {"x": str(first_three[1][0]), "y": str(first_three[1][1])},
        "P3": {"x": str(first_three[2][0]), "y": str(first_three[2][1])},
        "P4": {"x": str(fourth[0]), "y": str(fourth[1])},
    },
    "marking": {
        "minimum_basis_sha256": digest(BASIS),
        "normalized_supports": {"I7": "0", "I5": "1", "I4": "infinity"},
        "section_profiles_I7_I5_I4": profile_record,
        "section_intersection_gram": [list(row) for row in intersections],
    },
    "inputs": {
        "paths": [display_path(point_path), display_path(seed_path), display_path(BASIS)],
        "sha256": {
            display_path(point_path): digest(point_path),
            display_path(seed_path): digest(seed_path),
            display_path(BASIS): digest(BASIS),
        },
    },
    "proof_boundary": {
        "proved": (
            "This exact GF(p^d) point has source profile I7+I5+I4+8I1 and the "
            "displayed sections realize the pinned NS0024 MW4 component and intersection marking."
        ),
        "not_proved": (
            "This isolated finite-field point does not prove a one-dimensional modular component, "
            "a characteristic-zero lift, or Picard rank 19."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/adapt_lattice_foundry_ns0024_mw4_point_for_edge1.sage "
        "--point {} --seed {} --output {}"
    ).format(display_path(point_path), display_path(seed_path), display_path(output_path)),
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if args.check:
    if output_path.read_text() != serialized:
        raise SystemExit("adapted NS0024 MW4 source artifact is stale")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)

print(
    "NS0024MW4ADAPTER|p={}|degree={}|profile=I7+I5+I4+8I1|"
    "sections=4|marking=PASS|status=PASS".format(prime, modulus.degree()),
    flush=True,
)
