#!/usr/bin/env sage -python
"""Decode a Q80 closure RUR and certify third-q12 horizontal orbits."""

import argparse
import ast
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--surface", type=Path, required=True)
parser.add_argument("--scheme", type=Path, required=True)
parser.add_argument("--solution", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
for name in ("surface", "scheme", "solution", "output"):
    setattr(args, name, getattr(args, name).resolve())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


surface = json.loads(args.surface.read_text())
scheme = json.loads(args.scheme.read_text())
if surface.get("schema") != "elkies-k3.q80-fixed-u-marked-third-q12-search.v1":
    raise ValueError("unexpected surface schema")
if scheme.get("status") != "PASS_EXACT_POLYNOMIAL_CLOSURE_PRODUCER_EXPORTED":
    raise ValueError("closure scheme is not certified")
parameter = surface["parameters"][0]
modular = parameter["modular"][0]
prime = int(modular["prime"])
if scheme["specialization"] != {"u": parameter["u"], "prime": prime}:
    raise ValueError("surface/scheme specialization mismatch")

solution_text = args.solution.read_text().strip()
solution = ast.literal_eval(solution_text[:-1] if solution_text.endswith(":") else solution_text)
if solution[0] != 0:
    raise ArithmeticError("msolve did not return a zero-dimensional RUR")
payload = solution[1]
accepted_variable_orders = (
    ["l", "x0", "x1", "x2", "x3", "sat"],
    ["a", "n0", "n1", "n2", "n3", "sat"],
)
if int(payload[0]) != prime or payload[3][:6] not in accepted_variable_orders:
    raise ArithmeticError("unexpected RUR characteristic or variable order")
parametrization = payload[5]
if parametrization[0] != 1:
    raise ArithmeticError("unexpected number of RUR blocks")
elimination_data, denominator_data, coordinate_data = parametrization[1]
if denominator_data != [0, [1]] or len(coordinate_data) != 6:
    raise ArithmeticError("unsupported RUR denominator or coordinate count")

base_finite = GF(prime)
elimination_ring = PolynomialRing(base_finite, "T")
elimination = elimination_ring(elimination_data[1])
squarefree = elimination.squarefree_part()
factorization = tuple(squarefree.factor())
if any(int(exponent) != 1 or factor.degree() not in (1, 2) for factor, exponent in factorization):
    raise ArithmeticError("producer currently requires square-free linear/quadratic RUR support")
coordinate_polynomials = [elimination_ring(block[0][1]) for block in coordinate_data]

# Decode every factor in one common quadratic field.  This makes pairwise
# sums between different quadratic factors available without composita.
quadratic_factor = next((factor for factor, unused in factorization if factor.degree() == 2), None)
if quadratic_factor is None:
    field = base_finite
    extension_degree = 1
    field_generator = None
else:
    field = GF(prime**2, "r", modulus=quadratic_factor)
    extension_degree = 2
    field_generator = field.gen()

base = PolynomialRing(field, "W")
W = base.gen()
function_field = base.fraction_field()


def reduce_rational(value):
    value = QQ(value)
    denominator = field(value.denominator())
    if not denominator:
        raise ZeroDivisionError("equation denominator vanishes modulo p")
    return field(value.numerator()) / denominator


equation = parameter["exact_equations"]["second_q4"]
A = base([reduce_rational(value) for value in equation["A_coefficients_low_to_high"]])
B = base([reduce_rational(value) for value in equation["B_coefficients_low_to_high"]])
curve = EllipticCurve(function_field, [0, 0, 0, function_field(A), function_field(B)])
delta = 4 * A**3 + 27 * B**2
star_factor = next(factor.monic() for factor, exponent in delta.factor() if int(exponent) == 7)
if star_factor.degree() != 1:
    raise ArithmeticError("finite I1* factor is not rational")
star_root = -star_factor[0] / star_factor[1]
cubic_ring = PolynomialRing(field, "Xnode")
Xnode = cubic_ring.gen()
node_cubic = Xnode**3 + A(star_root) * Xnode + B(star_root)
singular_roots = node_cubic.gcd(node_cubic.derivative()).roots(multiplicities=False)
if len(singular_roots) != 1:
    raise ArithmeticError("finite I1* cubic has no unique singular x")
singular_x = singular_roots[0]


def coordinates(value):
    values = list(field(value).list()) + [base_finite.zero()] * extension_degree
    return [int(values[index]) for index in range(extension_degree)]


def polynomial_record(value):
    return [coordinates(coefficient) for coefficient in value.list()]


def rational_record(value):
    return {
        "numerator_coefficients_low_to_high": polynomial_record(value.numerator()),
        "denominator_coefficients_low_to_high": polynomial_record(value.denominator()),
    }


def point_key(point):
    if point.is_zero():
        return ("O",)
    return (
        tuple(map(tuple, rational_record(point[0])["numerator_coefficients_low_to_high"])),
        tuple(map(tuple, rational_record(point[0])["denominator_coefficients_low_to_high"])),
        tuple(map(tuple, rational_record(point[1])["numerator_coefficients_low_to_high"])),
        tuple(map(tuple, rational_record(point[1])["denominator_coefficients_low_to_high"])),
    )


def po_from_x(point):
    if point.is_zero():
        return None
    x_value = point[0]
    twice_intersection = max(
        x_value.denominator().degree(), x_value.numerator().degree() - 4
    )
    if twice_intersection < 0 or twice_intersection % 2:
        raise ArithmeticError("invalid section pole divisor")
    return int(twice_intersection // 2)


decoded = []
for factor, unused in factorization:
    roots = base.change_ring(field)(factor).roots(multiplicities=False)
    if len(roots) != factor.degree():
        raise ArithmeticError("RUR factor did not split in the common quadratic field")
    for root in roots:
        values = [
            -sum(field(coefficient) * root**index for index, coefficient in enumerate(polynomial.list()))
            for polynomial in coordinate_polynomials
        ]
        l, x0, x1, x2, x3, saturation = values
        if not l or saturation * l != 1:
            raise ArithmeticError("decoded point violates sat*l=1")
        x_value = base([x0, x1, x2, x3, l**2])
        square = x_value**3 + A * x_value + B
        y_coefficients = [field.zero() for unused in range(7)]
        y_coefficients[6] = l**3
        for degree in range(11, 5, -1):
            index = degree - 6
            partial = sum(y_coefficients[j] * W**j for j in range(7))
            y_coefficients[index] = (square[degree] - (partial**2)[degree]) / (
                2 * y_coefficients[6]
            )
        y_value = sum(y_coefficients[index] * W**index for index in range(7))
        if y_value**2 != square:
            raise ArithmeticError("decoded closure section fails substitution")
        point = curve(function_field(x_value), function_field(y_value))
        decoded.append(
            {
                "factor": str(factor.monic()),
                "factor_degree": int(factor.degree()),
                "root": coordinates(root),
                "l_x0_x1_x2_x3_sat": [coordinates(value) for value in values],
                "point": point,
                "x_coefficients_low_to_high": polynomial_record(x_value),
                "y_coefficients_low_to_high": polynomial_record(y_value),
            }
        )
if len(decoded) != squarefree.degree():
    raise ArithmeticError("not every square-free RUR point was decoded")

# Pairwise closure is the prime-independent replacement for the p=19-specific
# Q+P4 decoder.  Up to overall sign, only P_i +/- P_j is needed.
hits_by_unsigned_key = {}
tested = 0
for left_index, left in enumerate(decoded):
    for right_index in range(left_index, len(decoded)):
        right = decoded[right_index]
        for relative_sign in (1, -1):
            tested += 1
            horizontal = left["point"] + relative_sign * right["point"]
            if horizontal.is_zero() or po_from_x(horizontal) != 2:
                continue
            x_horizontal, y_horizontal = horizontal[0], horizontal[1]
            if (
                x_horizontal.denominator().degree() != 4
                or y_horizontal.denominator().degree() != 6
                or x_horizontal.numerator().degree() - x_horizontal.denominator().degree() != 4
                or y_horizontal.numerator().degree() - y_horizontal.denominator().degree() != 6
            ):
                continue
            if not x_horizontal.denominator()(star_root) or x_horizontal(star_root) == singular_x:
                continue
            fourth = 4 * horizontal
            height = QQ(4 + 2 * po_from_x(fourth)) / 16
            eighth = 2 * fourth
            if QQ(4 + 2 * po_from_x(eighth)) / 64 != height:
                raise ArithmeticError("fourth/eighth height replay disagrees")
            if height != 8:
                continue
            key = min(point_key(horizontal), point_key(-horizontal))
            hits_by_unsigned_key.setdefault(
                key,
                {
                    "left_index_zero_based": left_index,
                    "right_index_zero_based": right_index,
                    "relative_sign": relative_sign,
                    "left_factor": left["factor"],
                    "right_factor": right["factor"],
                    "P_dot_O": 2,
                    "canonical_height": "8",
                    "finite_I1star_identity": True,
                    "infinity_I3star_identity": True,
                    "x": rational_record(x_horizontal),
                    "y": rational_record(y_horizontal),
                    "literal_curve_substitution": True,
                    "unsigned_key": repr(key),
                },
            )

hits = list(hits_by_unsigned_key.values())
if not hits:
    raise ArithmeticError("closure support produced no target-profile horizontal")

# Quotient by Frobenius as well as global section sign.
key_to_index = {hit["unsigned_key"]: index for index, hit in enumerate(hits)}
orbit_representatives = []
seen = set()
for index, hit in enumerate(hits):
    if index in seen:
        continue
    # Reconstruct only to apply coefficientwise p-th power.
    def from_record(record):
        numerator = base(
            [sum(field(c) * field_generator**j for j, c in enumerate(value)) if extension_degree == 2 else field(value[0]) for value in record["numerator_coefficients_low_to_high"]]
        )
        denominator = base(
            [sum(field(c) * field_generator**j for j, c in enumerate(value)) if extension_degree == 2 else field(value[0]) for value in record["denominator_coefficients_low_to_high"]]
        )
        return function_field(numerator / denominator)

    point = curve(from_record(hit["x"]), from_record(hit["y"]))
    frobenius_point = curve(
        function_field(
            base([coefficient**prime for coefficient in point[0].numerator().list()])
            / base([coefficient**prime for coefficient in point[0].denominator().list()])
        ),
        function_field(
            base([coefficient**prime for coefficient in point[1].numerator().list()])
            / base([coefficient**prime for coefficient in point[1].denominator().list()])
        ),
    )
    orbit_keys = {
        repr(min(point_key(point), point_key(-point))),
        repr(min(point_key(frobenius_point), point_key(-frobenius_point))),
    }
    orbit_indices = {key_to_index[key] for key in orbit_keys if key in key_to_index}
    seen.update(orbit_indices)
    representative = dict(hit)
    representative.pop("unsigned_key")
    representative["frobenius_orbit_size_up_to_sign"] = len(orbit_indices)
    orbit_representatives.append(representative)

output = {
    "schema": "elkies-k3.q80-third-q12-polynomial-closure-rur-modp.v1",
    "status": "PASS_EXACT_THIRD_Q12_HORIZONTAL_FROM_COMMON_CLOSURE_PRODUCER",
    "specialization": {
        "u": parameter["u"],
        "prime": prime,
        "common_extension_degree": extension_degree,
        "common_extension_modulus": str(field.modulus()) if extension_degree == 2 else None,
    },
    "rur": {
        "quotient_degree_with_multiplicity": int(payload[2]),
        "squarefree_support_degree": int(squarefree.degree()),
        "elimination_polynomial": str(elimination),
        "factorization": [[str(factor.monic()), int(exponent)] for factor, exponent in factorization],
        "decoded_support_points": len(decoded),
        "all_support_points_replayed": True,
    },
    "decoded_sections": [
        {key: value for key, value in record.items() if key != "point"}
        for record in decoded
    ],
    "pairwise_producer": {
        "pairs_with_relative_sign_tested": tested,
        "unsigned_target_hits": len(hits),
        "frobenius_orbits_up_to_sign": len(orbit_representatives),
        "candidates": orbit_representatives,
        "acceptance_profile": {
            "P_dot_O": 2,
            "canonical_height": "8",
            "finite_I1star_component": "identity",
            "infinity_I3star_component": "identity",
        },
    },
    "inputs": [
        {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
        for path in (args.surface, args.scheme, args.solution)
    ],
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            "all square-free closure RUR points decoded in one common quadratic field",
            "complete pairwise sum/difference traversal of the decoded support",
            "literal target height, pole order, and additive-fibre identity profile",
        ],
        "not_proved": [
            "a child pencil/Jacobian at this prime",
            "canonical PGL2/Weierstrass alignment with p=19",
            "a characteristic-zero horizontal",
        ],
    },
    "reproduce": (
        "sage -python elkies-k3/scripts/certify_q80_third_q12_polynomial_closure_rur_modp.sage "
        f"--surface {args.surface} --scheme {args.scheme} --solution {args.solution} --output {args.output}"
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
print(
    f"Q80THIRDQ12COMMONRUR|u={parameter['u']}|prime={prime}|support={len(decoded)}|"
    f"tested={tested}|unsigned_hits={len(hits)}|orbits={len(orbit_representatives)}|"
    "status=PASS_EXACT_THIRD_Q12_HORIZONTAL_FROM_COMMON_CLOSURE_PRODUCER",
    flush=True,
)
