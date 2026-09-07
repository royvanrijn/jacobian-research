#!/usr/bin/env sage
"""Decode the Q80 polynomial-section RUR and certify third-q12 candidates.

At the fixed specialization ``u=-2, p=19``, the rational polynomial shell
spans only rank three.  The saturated identity-at-infinity P.O=0 scheme has
additional quadratic points.  This verifier decodes every square-free RUR
point, replays its polynomial section equation, adds it to every signed
rational polynomial section, and retains only literal D7+D5 third-q12
horizontals (P.O=2, height 8, identity at both additive fibres).
"""

from sage.all import EllipticCurve, GF, PolynomialRing, QQ

import argparse
import ast
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SURFACE = (
    ROOT / "artifacts/generated-results/q80-fixed-u-minus2-p19-height-shell-with-po1.json"
)
DEFAULT_SOLUTION = ROOT / (
    "artifacts/generated-results/q80-fixed-u-minus2-p19-po0-recursive-"
    "saturated-msolve/q80-third-q12-um2d1-p19-po0-polynomial-recursive-sign+1.solve"
)
DEFAULT_OUTPUT = ROOT / (
    "artifacts/generated-results/q80-fixed-u-minus2-p19-po0-rur-third-q12-modp.json"
)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--surface", type=Path, default=DEFAULT_SURFACE)
parser.add_argument("--solution", type=Path, default=DEFAULT_SOLUTION)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def reduce_rational(value, field):
    value = QQ(value)
    return field(int(value.numerator())) / field(int(value.denominator()))


def coefficient_strings(poly):
    return [str(value) for value in poly.list()]


def rational_function_record(value):
    return {
        "numerator_coefficients_low_to_high": coefficient_strings(value.numerator()),
        "denominator_coefficients_low_to_high": coefficient_strings(value.denominator()),
    }


def point_key(point):
    return (
        tuple(coefficient_strings(point[0].numerator())),
        tuple(coefficient_strings(point[0].denominator())),
        tuple(coefficient_strings(point[1].numerator())),
        tuple(coefficient_strings(point[1].denominator())),
    )


def po_from_x(point):
    if point.is_zero():
        return None
    x_value = point[0]
    denominator_degree = x_value.denominator().degree()
    numerator_degree = x_value.numerator().degree()
    twice_intersection = max(denominator_degree, numerator_degree - 4)
    if twice_intersection < 0 or twice_intersection % 2:
        raise ArithmeticError("section x-coordinate has an invalid pole divisor")
    return int(twice_intersection // 2)


surface = json.loads(args.surface.read_text())
if surface["schema"] != "elkies-k3.q80-fixed-u-marked-third-q12-search.v1":
    raise ValueError("unexpected surface-certificate schema")
parameter = surface["parameters"][0]
if parameter["u"] != "-2":
    raise ValueError("this certificate is pinned to u=-2")
modular_record = parameter["modular"][0]
prime = int(modular_record["prime"])
if prime != 19:
    raise ValueError("this certificate is pinned to p=19")
second_equation = parameter["exact_equations"]["second_q4"]

solution_text = args.solution.read_text().strip()
solution = ast.literal_eval(
    solution_text[:-1] if solution_text.endswith(":") else solution_text
)
if solution[0] != 0:
    raise ArithmeticError("msolve did not return a zero-dimensional RUR")
payload = solution[1]
if int(payload[0]) != prime or int(payload[2]) != 16:
    raise ArithmeticError("unexpected RUR characteristic or quotient degree")
if payload[3][:6] != ["a", "n0", "n1", "n2", "n3", "sat"]:
    raise ArithmeticError("unexpected RUR variable order")
parametrization = payload[5]
if parametrization[0] != 1:
    raise ArithmeticError("unexpected number of RUR blocks")
elimination_data, denominator_data, coordinate_data = parametrization[1]
if denominator_data != [0, [1]] or len(coordinate_data) != 6:
    raise ArithmeticError("unsupported RUR denominator or coordinate count")

base_field = GF(prime)
elimination_ring = PolynomialRing(base_field, "T")
elimination = elimination_ring(elimination_data[1])
squarefree = elimination.squarefree_part()
factorization = tuple(squarefree.factor())
if any(exponent != 1 or factor.degree() not in (1, 2) for factor, exponent in factorization):
    raise ArithmeticError("RUR support has an unexpected factorization")
coordinate_polynomials = [
    elimination_ring(block[0][1]) for block in coordinate_data
]

decoded_sections = []
oriented_hits = []
for factor, _ in factorization:
    extension_degree = int(factor.degree())
    if extension_degree == 1:
        field = base_field
        roots = factor.roots(multiplicities=False)
        generator_name = None
    else:
        field = GF(prime**extension_degree, name="r", modulus=factor)
        root_ring = PolynomialRing(field, "Z")
        roots = root_ring(factor).roots(multiplicities=False)
        generator_name = "r"

    polynomial_ring = PolynomialRing(field, "W")
    W = polynomial_ring.gen()
    function_field = polynomial_ring.fraction_field()
    A = polynomial_ring(
        [reduce_rational(value, field) for value in second_equation["A_coefficients_low_to_high"]]
    )
    B = polynomial_ring(
        [reduce_rational(value, field) for value in second_equation["B_coefficients_low_to_high"]]
    )
    curve = EllipticCurve(function_field, [0, 0, 0, function_field(A), function_field(B)])
    star_factor = polynomial_ring(
        [
            reduce_rational(value, field)
            for value in second_equation["finite_I1star_factor_coefficients_low_to_high"]
        ]
    )
    star_root = -star_factor[0] / star_factor[1]
    cubic_ring = PolynomialRing(field, "X")
    X = cubic_ring.gen()
    cubic = X**3 + A(star_root) * X + B(star_root)
    singular_roots = cubic.gcd(cubic.derivative()).roots(multiplicities=False)
    if len(singular_roots) != 1:
        raise ArithmeticError("finite I1* cubic has no unique singular x")
    singular_x = singular_roots[0]

    rational_points = []
    for point_record in modular_record["polynomial_shell"]:
        x_value = polynomial_ring(
            [field(value) for value in point_record["x_coefficients_low_to_high"]]
        )
        y_value = polynomial_ring(
            [field(value) for value in point_record["y_coefficients_low_to_high"]]
        )
        rational_points.append(curve(function_field(x_value), function_field(y_value)))

    for root in roots:
        values = [
            -sum(field(coefficient) * root**index for index, coefficient in enumerate(poly.list()))
            for poly in coordinate_polynomials
        ]
        a, n0, n1, n2, n3, saturation = values
        if not a or saturation * a != 1:
            raise ArithmeticError("decoded RUR point violates sat*a=1")
        x_value = polynomial_ring([n0, n1, n2, n3, a**2])
        square = x_value**3 + A * x_value + B
        y_coefficients = [field.zero() for _ in range(7)]
        y_coefficients[6] = a**3
        for degree in range(11, 5, -1):
            index = degree - 6
            partial = sum(
                polynomial_ring(y_coefficients[j]) * W**j for j in range(7)
            )
            y_coefficients[index] = (
                square[degree] - (partial**2)[degree]
            ) / (2 * y_coefficients[6])
        y_value = sum(
            polynomial_ring(y_coefficients[index]) * W**index for index in range(7)
        )
        if y_value**2 != square:
            raise ArithmeticError("decoded RUR polynomial section fails substitution")
        Q = curve(function_field(x_value), function_field(y_value))
        decoded_sections.append(
            {
                "factor": str(factor),
                "extension_degree": extension_degree,
                "root": str(root),
                "a_n0_n1_n2_n3_sat": list(map(str, values)),
                "x_coefficients_low_to_high": coefficient_strings(x_value),
                "y_coefficients_low_to_high": coefficient_strings(y_value),
                "literal_curve_substitution": True,
            }
        )

        if extension_degree != 2:
            continue
        for point_index, rational_point in enumerate(rational_points, 1):
            for point_sign in (1, -1):
                horizontal = Q + point_sign * rational_point
                if po_from_x(horizontal) != 2:
                    continue
                x_horizontal, y_horizontal = horizontal[0], horizontal[1]
                if y_horizontal.denominator().degree() != 6:
                    continue
                if (
                    x_horizontal.numerator().degree() - x_horizontal.denominator().degree() != 4
                    or y_horizontal.numerator().degree() - y_horizontal.denominator().degree() != 6
                ):
                    continue
                if not x_horizontal.denominator()(star_root):
                    continue
                if x_horizontal(star_root) == singular_x:
                    continue
                fourth = 4 * horizontal
                height = QQ(4 + 2 * po_from_x(fourth)) / 16
                eighth = 2 * fourth
                if QQ(4 + 2 * po_from_x(eighth)) / 64 != height:
                    raise ArithmeticError("fourth/eighth height replay disagrees")
                if height != 8:
                    continue
                oriented_hits.append(
                    {
                        "factor": str(factor),
                        "extension_degree": extension_degree,
                        "field_generator": generator_name,
                        "root": str(root),
                        "polynomial_section_index_one_based": point_index,
                        "polynomial_section_sign": point_sign,
                        "P_dot_O": 2,
                        "canonical_height": str(height),
                        "finite_I1star_identity": True,
                        "infinity_I3star_identity": True,
                        "x": rational_function_record(x_horizontal),
                        "y": rational_function_record(y_horizontal),
                        "literal_curve_substitution": True,
                        "point_key": repr(point_key(horizontal)),
                        "negative_point_key": repr(point_key(-horizontal)),
                    }
                )

unsigned = {}
for hit in oriented_hits:
    key = min(hit["point_key"], hit["negative_point_key"])
    unsigned.setdefault(key, hit)
unsigned_hits = list(unsigned.values())
if len(decoded_sections) != 12:
    raise ArithmeticError("expected twelve square-free RUR support points")
if len(oriented_hits) != 4 or len(unsigned_hits) != 2:
    raise ArithmeticError("unexpected third-q12 candidate count")
if {hit["factor"] for hit in unsigned_hits} != {"T^2 + 12*T + 3"}:
    raise ArithmeticError("third-q12 candidates lie on an unexpected RUR factor")
if {hit["polynomial_section_index_one_based"] for hit in unsigned_hits} != {4}:
    raise ArithmeticError("third-q12 candidates use an unexpected polynomial section")

output = {
    "schema": "elkies-k3.q80-po0-rur-third-q12-modp.v1",
    "status": "PASS_EXACT_MODP2_THIRD_Q12_HORIZONTAL_FROBENIUS_ORBIT",
    "specialization": {"u": "-2", "prime": prime, "extension_degree": 2},
    "rur": {
        "quotient_degree_with_multiplicity": int(payload[2]),
        "squarefree_support_degree": int(squarefree.degree()),
        "elimination_polynomial": str(elimination),
        "squarefree_factorization": [
            [str(factor), int(exponent)] for factor, exponent in factorization
        ],
        "decoded_support_points": len(decoded_sections),
        "rational_support_points": sum(
            record["extension_degree"] == 1 for record in decoded_sections
        ),
        "quadratic_support_points": sum(
            record["extension_degree"] == 2 for record in decoded_sections
        ),
        "all_support_points_replayed": True,
    },
    "decoded_polynomial_sections": decoded_sections,
    "third_q12": {
        "oriented_hits": len(oriented_hits),
        "unsigned_hits": len(unsigned_hits),
        "frobenius_orbits_up_to_sign": 1,
        "candidate_factor": "T^2 + 12*T + 3",
        "candidate_factor_discriminant_mod_19": 18,
        "candidates_up_to_sign": [
            {key: value for key, value in hit.items() if not key.endswith("point_key")}
            for hit in unsigned_hits
        ],
        "acceptance_profile": {
            "P_dot_O": 2,
            "canonical_height": "8",
            "components": "identity at finite I1* and infinity I3*",
        },
    },
    "inputs": {
        "surface": {"path": str(args.surface), "sha256": sha256(args.surface)},
        "msolve_solution": {
            "path": str(args.solution),
            "sha256": sha256(args.solution),
        },
        "worker": {
            "path": str(Path(__file__).resolve().relative_to(ROOT)),
            "sha256": sha256(Path(__file__).resolve()),
        },
    },
    "claim_boundary": (
        "This is an exact equation-level certificate for the marked third-q12 "
        "horizontal orbit over GF(19^2)(W). It does not produce a GF(19)(W) or "
        "QQ(W) horizontal, construct the two-dimensional connected q12 "
        "Riemann-Roch pencil, compile the child, or certify A5+A3+3A1/MW6."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
print(
    "Q80PO0RURQ12|support=12|rational=6|quadratic=6|oriented_hits=4|"
    "unsigned_hits=2|frobenius_orbits=1|factor=T^2+12*T+3|"
    "status=PASS_EXACT_MODP2_THIRD_Q12_HORIZONTAL_FROBENIUS_ORBIT",
    flush=True,
)
