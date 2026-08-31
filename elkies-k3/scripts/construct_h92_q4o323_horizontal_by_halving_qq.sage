#!/usr/bin/env sage
"""Construct the corrected q4/o323 horizontal by exact MW halving.

The marked lattice relation is

    2*T = P8^- + 2*P18^- + P33^- - 2*C7

modulo the q4/o208 trivial lattice.  All points on the right are already
exact over QQ(t).  We therefore factor only the univariate duplication
quartic and select the rational half; no Groebner basis or section shell is
used here.
"""

import hashlib
import json
import time
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
if not (ROOT / "MATH_STATUS.json").exists():
    ROOT = Path.cwd()
LOCAL = ROOT / "artifacts/local/elkies-k3"
COMPACT_PATH = LOCAL / "q4o208-compact-weierstrass-qq.json"
BRANCHES_PATH = LOCAL / "q4o208-q4o323-horizontal-resolved-qq.json"
MOD131_PATH = LOCAL / "q4o208-physical-q4o323-horizontal-mod131.json"
MARKING_PATH = LOCAL / "q24-2a5-physical-q4o208-equation-marking-qq.json"
OUTPUT = LOCAL / "q4o208-q4o323-horizontal-by-halving-qq.json"
INPUTS = (COMPACT_PATH, BRANCHES_PATH, MOD131_PATH, MARKING_PATH)

started = time.monotonic()
compact = json.loads(COMPACT_PATH.read_text())
branches = json.loads(BRANCHES_PATH.read_text())
mod131 = json.loads(MOD131_PATH.read_text())
marking = json.loads(MARKING_PATH.read_text())
assert compact["status"] == "PASS_EXACT_QQ_Q4O208_COMPACT_WEIERSTRASS_NORMALIZATION"
assert branches["status"] == "PASS_EXACT_QQ_Q4O323_RESOLVED_SIMPLE_POLE_HORIZONTAL"
assert mod131["status"] == "PASS_EXACT_Q4O323_POLYNOMIAL_SECTION_SUBGROUP_OBSTRUCTION"
assert marking["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O208_C5_EQUATION_MARKING"

R = PolynomialRing(QQ, "t")
t = R.gen()
K = R.fraction_field()
model = compact["compact_model"]
A = K(R([QQ(value) for value in model["A_coefficients_low_to_high"]]))
B = K(R([QQ(value) for value in model["B_coefficients_low_to_high"]]))
E = EllipticCurve(K, [0, 0, 0, A, B])


def polynomial(values):
    return R([QQ(value) for value in values])


def branch_point(index, sign):
    record = next(
        record for record in branches["exact_QQ_horizontal_sections"]
        if int(record["branch_index"]) == index
    )
    Z = polynomial(record["Z_coefficients_low_to_high"])
    X = polynomial(record["X_coefficients_low_to_high"])
    Y = polynomial(record["Y_coefficients_low_to_high"])
    return E(K(X/Z**2), K(sign*Y/Z**3))


c7_record = compact["transported_exact_section"]
assert c7_record["source_label"] == "old_A11_component_7_on_C5_pointed_child"
C7 = E(
    K(polynomial(c7_record["x_coefficients_low_to_high"])),
    K(polynomial(c7_record["y_coefficients_low_to_high"])),
)
P8 = branch_point(8, 1)
P18 = branch_point(18, 1)
P33 = branch_point(33, 1)
double_candidates = [
    {
        "orientation": 0,
        "word": "-P8 - 2*P18 - P33 - 2*C7",
        "point": -P8 - 2*P18 - P33 - 2*C7,
        "branch_coefficients": (-1, -2, -1),
    },
    {
        "orientation": 1,
        "word": "P8 + 2*P18 + P33 - 2*C7",
        "point": P8 + 2*P18 + P33 - 2*C7,
        "branch_coefficients": (1, 2, 1),
    },
]
assert all(not candidate["point"].is_zero() for candidate in double_candidates)

# If Q=(x,y), the x-coordinate of 2Q is
# (x^4-2*A*x^2-8*B*x+A^2)/(4*(x^3+A*x+B)).
S = PolynomialRing(K, "z")
z = S.gen()


def rational_halves(double):
    x2 = K(double[0])
    quartic = (
        z**4 - 4*x2*z**3 - 2*A*z**2
        - (4*A*x2 + 8*B)*z + A**2 - 4*B*x2
    )
    factorization = quartic.factor()
    halves = []
    for factor, exponent in factorization:
        if factor.degree() != 1:
            continue
        x = K(-factor[0]/factor[1])
        halves.extend(point for point in E.lift_x(x, all=True) if 2*point == double)
    return factorization, halves


def intersection_with_zero(point):
    x = K(point[0])
    numerator_degree = int(x.numerator().degree())
    denominator_degree = int(x.denominator().degree())
    assert denominator_degree % 2 == 0
    infinity_excess = max(0, numerator_degree-denominator_degree-4)
    assert infinity_excess % 2 == 0
    return denominator_degree//2 + infinity_excess//2


for candidate in double_candidates:
    factorization, halves = rational_halves(candidate["point"])
    candidate["factorization"] = factorization
    candidate["halves"] = halves
    candidate["P_dot_O_values"] = [intersection_with_zero(point) for point in halves]

expected_P_dot_O = int(mod131["target"]["P_dot_O"])
selected = [
    (candidate, point)
    for candidate in double_candidates for point in candidate["halves"]
    if intersection_with_zero(point) == expected_P_dot_O
]
if len(selected) != 1:
    raise ArithmeticError(f"expected one marked rational half, found {len(selected)}")
selected_candidate, target = selected[0]
double_target = selected_candidate["point"]
factorization = selected_candidate["factorization"]
halves = selected_candidate["halves"]

# Verify the doubled relation in the exact marked NS lattice.  Primitivity of
# the trivial lattice also rules out a torsion ambiguity between the two
# sections having this double.
selected_orientation = int(selected_candidate["orientation"])
branch_maps = mod131["exact_simple_pole_mapping"]["maps_by_global_component_orientation"][selected_orientation]


def branch_class(index, sign):
    record = next(
        record for record in branch_maps
        if int(record["stored_branch_index"]) == index
        and int(record["Y_sign_relative_to_stored"]) == sign
    )
    return vector(ZZ, record["NS_coordinates"])


target_class = vector(ZZ, mod131["target"]["NS_coordinates"])
c7_class = vector(ZZ, marking["old_A11_component_7_on_C5_pointed_child"]["NS_coordinates"])
coefficients = selected_candidate["branch_coefficients"]
rhs_class = (
    coefficients[0]*branch_class(8, 1)
    + coefficients[1]*branch_class(18, 1)
    + coefficients[2]*branch_class(33, 1)
    - 2*c7_class
)
fibre_class = vector(ZZ, [1, 0] + [0]*17)
fibres = marking["physical_fibres"]
first_cycle = [vector(ZZ, value) for value in fibres["first_old_I6_I4"]["components_in_cycle_order"]]
second_cycle = [vector(ZZ, value) for value in fibres["second_old_I6_I4"]["components_in_cycle_order"]]
special_cycle = [vector(ZZ, value) for value in fibres["special_I4"]["components_in_cycle_order"]]
trivial_difference = 2*target_class-rhs_class

zero_class = vector(ZZ, [-1, 1] + [0]*17)
root_components = []
for fibre in fibres.values():
    identity = int(fibre["identity_component_index"])
    root_components.extend(
        vector(ZZ, component)
        for index, component in enumerate(fibre["components_in_cycle_order"])
        if index != identity
    )
trivial_matrix = matrix(ZZ, [fibre_class, zero_class] + root_components)
trivial_smith_invariants = [int(value) for value in trivial_matrix.elementary_divisors()]
assert trivial_smith_invariants == [1]*11
trivial_coefficients = trivial_matrix.transpose().solve_right(trivial_difference)
assert all(value.denominator() == 1 for value in trivial_coefficients)
assert trivial_coefficients*trivial_matrix == trivial_difference


def rational_record(value):
    value = K(value)
    numerator = R(value.numerator())
    denominator = R(value.denominator())
    coefficients = list(numerator) + list(denominator)
    return {
        "numerator_coefficients_low_to_high": [str(value) for value in numerator.list()],
        "denominator_coefficients_low_to_high": [str(value) for value in denominator.list()],
        "degrees_numerator_denominator": [int(numerator.degree()), int(denominator.degree())],
        "maximum_rational_bit_length": max(
            max(abs(value.numerator()).nbits(), value.denominator().nbits())
            for value in coefficients
        ),
    }


payload = {
    "schema": "elkies-k3.h3-q4o208-q4o323-horizontal-by-halving-qq.v1",
    "status": "PASS_EXACT_QQ_Q4O323_HORIZONTAL_BY_MW_HALVING",
    "inputs": [
        {"path": str(path.relative_to(ROOT)), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        for path in INPUTS
    ],
    "method": {
        "selected_global_component_orientation": selected_orientation,
        "MW_relation_modulo_trivial_lattice": "2*T = " + selected_candidate["word"],
        "orientation_candidates": [
            {
                "global_component_orientation": int(candidate["orientation"]),
                "word": candidate["word"],
                "duplication_quartic_factor_degrees": [
                    int(factor.degree())
                    for factor, exponent in candidate["factorization"]
                    for unused in range(exponent)
                ],
                "rational_half_P_dot_O_values": [
                    int(value) for value in candidate["P_dot_O_values"]
                ],
            }
            for candidate in double_candidates
        ],
        "selected_P_dot_O": int(intersection_with_zero(target)),
        "large_Groebner_required": False,
        "duplication_quartic_factor_degrees": [
            int(factor.degree()) for factor, exponent in factorization for unused in range(exponent)
        ],
        "rational_halves": len(halves),
        "marked_NS_doubling_relation_verified": True,
        "trivial_lattice_relation_coefficients": [int(value) for value in trivial_coefficients],
        "trivial_lattice_smith_invariants": trivial_smith_invariants,
        "MW_torsion_ambiguity": "none: the trivial lattice is primitive in the marked NS lattice",
        "runtime_seconds": time.monotonic()-started,
    },
    "exact_QQ_horizontal": {
        "x": rational_record(target[0]),
        "y": rational_record(target[1]),
        "exact_compact_weierstrass_identity": bool(
            target[1]**2 == target[0]**3 + A*target[0] + B
        ),
        "exact_doubling_identity": bool(2*target == double_target),
    },
    "proof_boundary": (
        "This proves the exact QQ(t) section with the corrected q4/o323 marked NS class: "
        "the doubled lattice relation is exact and primitivity of the trivial lattice removes "
        "torsion ambiguity. The q4/o323 resolved RR pencil remains a separate gate."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q4O323HALF|status={}|degrees_x={}|degrees_y={}|bits_x={}|bits_y={}|runtime={:.3f}".format(
        payload["status"],
        payload["exact_QQ_horizontal"]["x"]["degrees_numerator_denominator"],
        payload["exact_QQ_horizontal"]["y"]["degrees_numerator_denominator"],
        payload["exact_QQ_horizontal"]["x"]["maximum_rational_bit_length"],
        payload["exact_QQ_horizontal"]["y"]["maximum_rational_bit_length"],
        payload["method"]["runtime_seconds"],
    ),
    flush=True,
)
