#!/usr/bin/env sage
"""Certify the first q=4 collision on the exact unmarked Q80 coefficient curve.

The reconstructed slope-8/87 coefficient curve is an exact rational map into
the Q80 ``E6+D5+A3`` chart.  This script makes one further, equation-level
check: its four coefficients lie identically on the degree-eight divisor
where the first q=4 Jacobian gains an I5 fibre.  Thus the explicit coordinate
``U=(x-T)/T^2`` gives a K3 fibration with generic reducible fibres
``I5*+I5`` and eight remaining nodal fibres over ``QQ(u)``.

This is deliberately unmarked.  It does not construct the old MW sections,
the later q12/q4/q6 pencils, a rootless equation, or a bisection.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *
from sage.misc.persist import load


ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--parameter",
    type=Path,
    default=ROOT / "artifacts/generated-results/"
    "q80-cm24-slope-8-87-qq-PDQE-parameter.json",
)
parser.add_argument(
    "--output",
    type=Path,
    default=ROOT / "artifacts/generated-results/"
    "q80-unmarked-first-q4-collision-qq.json",
)
parser.add_argument(
    "--witness-u",
    type=int,
    default=0,
    help="rational parameter value proving the required generic nonvanishings",
)
args = parser.parse_args()

load(str(ROOT / "elkies-k3/scripts/derive_q80_first_q4_pencil.sage"))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


payload = json.loads(args.parameter.read_text())
if payload.get("schema") != "q80-cm24-qq-PDQE-parameter-v1":
    raise ValueError("unexpected Q80 parameter schema")

parameter_ring = PolynomialRing(QQ, "u")
u = parameter_ring.gen()
parameter_field = parameter_ring.fraction_field()


def parameter_function(name):
    record = payload["original_functions"][name]
    numerator = parameter_ring(record["numerator"])
    denominator = parameter_ring(record["denominator"])
    if denominator == 0 or numerator.gcd(denominator) != 1:
        raise ArithmeticError(f"{name}(u) is not reduced")
    if [int(numerator.degree()), int(denominator.degree())] != record["degrees"]:
        raise ArithmeticError(f"{name}(u) has inconsistent recorded degrees")
    return parameter_field(numerator / denominator)


d_u, p_u, q_u, e_u = tuple(
    parameter_function(name) for name in ("d", "p", "q", "e")
)
specialization_homomorphism = parameters.hom(
    (d_u, p_u, q_u, e_u), parameter_field
)


def specialize_coefficient(value):
    """Evaluate an element of QQ(d,p,q,e) in QQ(u) without coercion guesses."""
    value = K(value)
    numerator = parameters(value.numerator())
    denominator = parameters(value.denominator())
    numerator_value = specialization_homomorphism(numerator)
    denominator_value = specialization_homomorphism(denominator)
    if denominator_value == 0:
        raise ZeroDivisionError("Q80 coefficient curve hits a chart denominator")
    return parameter_field(numerator_value / denominator_value)


collision_value = specialize_coefficient(rank19_collision_factor)
assert collision_value == 0
print("Q80UNMARKEDFIRSTQ4|collision_identity=0", flush=True)

def evaluate_at_witness(value, values):
    """Specialize QQ(d,p,q,e) to QQ at one exact parameter value."""
    value = K(value)
    map_to_qq = parameters.hom(values, QQ)
    numerator = map_to_qq(parameters(value.numerator()))
    denominator = map_to_qq(parameters(value.denominator()))
    if denominator == 0:
        raise ZeroDivisionError("witness lies on a Q80 chart denominator")
    return QQ(numerator / denominator)


witness_parameter = QQ(args.witness_u)
witness_values = tuple(
    QQ(parameter_function(name)(witness_parameter))
    for name in ("d", "p", "q", "e")
)
child_ring = PolynomialRing(QQ, "U")
U_u = child_ring.gen()
first_a = child_ring([evaluate_at_witness(value, witness_values) for value in jacobian_a.list()])
first_b = child_ring([evaluate_at_witness(value, witness_values) for value in jacobian_b.list()])
first_delta = 4 * first_a**3 + 27 * first_b**2
u0 = witness_values[0] - 1
linear = U_u - u0

# The ambient identity is Delta=(U-d+1)^4 R9.  The exact collision identity
# above gives the fifth factor over QQ(u).  This one rational specialization
# proves all remaining open nonvanishing conditions in the generic function
# field, without expanding the enormous QQ(u) child coefficients.
assert first_a.degree() == 6 and first_b.degree() == 9 and first_delta.degree() == 13
assert first_a(u0) != 0 and first_b(u0) != 0
assert first_delta % linear**5 == 0
assert first_delta % linear**6 != 0
residual = first_delta // linear**5
assert residual.degree() == 8
assert residual.gcd(residual.derivative()) == 1
assert residual(u0) != 0

# The degrees give the minimal infinity valuations (2,3,11), hence I5*.
infinity_valuations = (
    int(8 - first_a.degree()),
    int(12 - first_b.degree()),
    int(24 - first_delta.degree()),
)
assert infinity_valuations == (2, 3, 11)

output = {
    "schema": "q80-unmarked-first-q4-collision-qq-v1",
    "status": "PASS_EXACT_UNMARKED_FIRST_Q4_COLLISION",
    "inputs": {
        "parameter": {
            "path": str(args.parameter.relative_to(ROOT)),
            "sha256": sha256(args.parameter),
        },
        "pencil_derivation": {
            "path": "elkies-k3/scripts/derive_q80_first_q4_pencil.sage",
            "sha256": sha256(ROOT / "elkies-k3/scripts/derive_q80_first_q4_pencil.sage"),
        },
    },
    "coordinate": "U=(x-T)/T^2",
    "collision": {
        "ambient_divisor_degree": int(rank19_collision_factor.total_degree()),
        "ambient_divisor_terms": int(len(rank19_collision_factor.dict())),
        "exact_substitution": "0",
        "finite_fibre": {
            "base_point": "U=d(u)-1",
            "valuations_A_B_Delta": [int(0), int(0), int(5)],
            "Kodaira": "I5",
        },
        "remaining_finite_discriminant_degree": int(residual.degree()),
        "remaining_finite_discriminant_squarefree": True,
        "genericity_witness": {
            "u": str(witness_parameter),
            "d_p_q_e": [str(value) for value in witness_values],
            "checks": [
                "fifth but not sixth finite discriminant factor",
                "remaining degree-eight discriminant is squarefree",
                "finite A and B are units",
            ],
        },
        "infinity_valuations_A_B_Delta": list(infinity_valuations),
        "infinity_Kodaira": "I5*",
        "generic_reducible_fibres": ["I5*", "I5"],
        "generic_root_lattice": "D9+A4",
        "remaining_fibres": "8 I1",
    },
    "claim_boundary": {
        "proved": [
            "exact unmarked Q80 coefficient curve lies on the first-q4 collision divisor",
            "the first-q4 pencil is explicit over QQ(u)",
            "the resulting generic reducible fibre configuration is I5*+I5+8I1",
        ],
        "not_proved": [
            "global marked Mordell-Weil sections on the Q80 curve",
            "the second or later pencil over QQ(u)",
            "a rootless MW17 equation",
            "a bisection, extension collision, or generic rank-19 family",
        ],
    },
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")

print(
    "Q80UNMARKEDFIRSTQ4|collision=0|finite=I5|infinity=I5*|"
    "residual=8I1|ADE=D9+A4|status=PASS_EXACT_UNMARKED_FIRST_Q4_COLLISION",
    flush=True,
)
