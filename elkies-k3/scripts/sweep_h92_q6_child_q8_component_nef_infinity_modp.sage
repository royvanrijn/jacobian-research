#!/usr/bin/env sage -python
"""Sweep signed infinity lattices for the component-nef q8 finite module.

This works in the standard q_regular frame; fibrewise translation to the
physical component-nef chord preserves the base cover and the infinity
valuation profile.  It allows negative required orders (allowed poles), unlike
the earlier exploratory probe.

For every signed order pair (a_order,b_order), it intersects the exact finite
module with the corresponding bounded infinity lattice.  Two-dimensional
kernels can optionally be screened by the branch degree of every ratio level.
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix


ROOT = Path(__file__).resolve().parents[2]
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
MARKING = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-marking.json"
NORMALIZER = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-q-pole-normalization-crt.json"
DEFAULT_FINITE = ROOT / "artifacts/local/elkies-k3/q8-qregular-finite-component-nef.json"
DEFAULT_OUTPUT = ROOT / "artifacts/local/elkies-k3/q8-component-nef-infinity-sweep.json"


def coefficient(field, value):
    value = QQ(value)
    denominator = field(ZZ(value.denominator()))
    if not denominator:
        raise ValueError("prime divides an input denominator")
    return field(ZZ(value.numerator())) / denominator


def polynomial(ring, finite, coefficients):
    return ring([coefficient(finite, value) for value in coefficients])


def monic_power_root(value, exponent):
    root = value.parent().one()
    for factor, multiplicity in value.factor():
        assert multiplicity % exponent == 0
        root *= factor.monic() ** (multiplicity // exponent)
    return root.monic()


def degree_or_minus_one(value):
    return -1 if not value else value.degree()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=43)
parser.add_argument("--finite", type=Path, default=DEFAULT_FINITE)
parser.add_argument("--max-s-degree", type=int, default=8)
parser.add_argument("--max-t-degree", type=int, default=8)
parser.add_argument("--a-order-min", type=int, default=-8)
parser.add_argument("--a-order-max", type=int, default=2)
parser.add_argument("--b-order-min", type=int, default=-8)
parser.add_argument("--b-order-max", type=int, default=2)
parser.add_argument("--screen-two-dimensional", action="store_true")
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
if not ZZ(args.prime).is_prime() or args.prime in (2, 3):
    raise ValueError("prime must be odd and different from 3")
if min(args.max_s_degree, args.max_t_degree) < 0:
    raise ValueError("degree bounds must be nonnegative")
args.finite = args.finite.resolve()
args.output = args.output.resolve()

child = json.loads(CHILD.read_text())
marking = json.loads(MARKING.read_text())
normalizer = json.loads(NORMALIZER.read_text())
finite_module = json.loads(args.finite.read_text())
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert marking["status"] == "PASS_EXACT_Q6_CHILD_Q8_MARKING"
assert normalizer["status"] == "PASS_EXACT_CRT_PRINCIPAL_PART_NORMALIZATION"
assert finite_module["status"] == "PASS_EXACT_Q_REGULAR_FINITE_MODULE"
assert finite_module["module"]["basis"] == [
    ["1", "lift(R/Nx)"],
    ["0", "f_II*^2*f_IV*^2"],
]

finite = GF(args.prime)
ring = PolynomialRing(finite, "T")
T = ring.gen()
field = ring.fraction_field()
section = marking["selected_q8"]["relative_child_section_standard_jacobian_coordinates"]
nx = polynomial(ring, finite, section["x_numerator_coefficients_low_to_high"])
dx = polynomial(ring, finite, section["x_denominator_coefficients_low_to_high"])
ny = polynomial(ring, finite, section["y_numerator_coefficients_low_to_high"])
dy = polynomial(ring, finite, section["y_denominator_coefficients_low_to_high"])
sx, sy = field(nx)/field(dx), field(ny)/field(dy)
h = monic_power_root(dx, 2)
R = polynomial(ring, finite, normalizer["normalizer"]["R_coefficients_low_to_high"])
assert (R*h*dy-ny) % nx == 0
ii = polynomial(ring, finite, PolynomialRing(QQ, "T")(next(
    item for item in child["finite_fibres"] if item["kodaira"] == "II*"
)["factor"]).list())
iv = polynomial(ring, finite, PolynomialRing(QQ, "T")(next(
    item for item in child["finite_fibres"] if item["kodaira"] == "IV*"
)["factor"]).list())
M = ii**2*iv**2
C1 = polynomial(
    ring, finite,
    finite_module["module"]["first_C_lift_coefficients_low_to_high"],
)
assert C1.degree() < 4 and (C1-R*nx.inverse_mod(M)) % M == 0

p0 = -sy/sx
alpha = -p0/field(h)-field(R)/field(nx)
beta = field(T**2)/field(h)
alpha_num, alpha_den = ring(alpha.numerator()), ring(alpha.denominator())
beta_num, beta_den = ring(beta.numerator()), ring(beta.denominator())
assert alpha_den(0) and beta_den(0)

labels = [("s", exponent) for exponent in range(args.max_s_degree+1)] + [
    ("t", exponent) for exponent in range(args.max_t_degree+1)
]


def module_pair(label):
    kind, exponent = label
    if kind == "s":
        return T**exponent, C1*T**exponent
    return ring.zero(), M*T**exponent


pairs = [module_pair(label) for label in labels]
a_values = [C*alpha_den+B*alpha_num for B, C in pairs]
b_values = [B*beta_num for B, C in pairs]


def rows_for(values, denominator_degree, required_order):
    cutoff = denominator_degree-required_order
    top = max((degree_or_minus_one(value) for value in values), default=-1)
    rows = [
        [value[degree] if degree <= value.degree() else finite.zero() for value in values]
        for degree in range(cutoff+1, top+1)
    ]
    return rows, cutoff, top


def rational_pair(row):
    s_poly = sum(
        row[index]*T**exponent
        for index, (kind, exponent) in enumerate(labels)
        if kind == "s"
    )
    t_poly = sum(
        row[index]*T**exponent
        for index, (kind, exponent) in enumerate(labels)
        if kind == "t"
    )
    B = s_poly
    C = C1*s_poly+M*t_poly
    a = field(C)-field(B)*p0/field(h)-field(B*R)/field(nx)
    b = field(B)/field(h)
    return a, b


Acurve = polynomial(ring, finite, child["minimal_short_weierstrass"]["A_coefficients_low_to_high"])
Bcurve = polynomial(ring, finite, child["minimal_short_weierstrass"]["B_coefficients_low_to_high"])


def branch_screen(kernel):
    a0, b0 = rational_pair(kernel.row(0))
    a1, b1 = rational_pair(kernel.row(1))
    x_ring = PolynomialRing(field, "x")
    x = x_ring.gen()
    histogram = {}
    good = []
    singular = []
    for level in finite:
        denominator = b1-field(level)*b0
        if not denominator:
            singular.append(int(level))
            continue
        m_value = -(a1-field(level)*a0)/denominator
        y = x_ring(m_value)*(x-x_ring(sx))-x_ring(sy)
        relation = y**2-x**3-x_ring(Acurve)*x-x_ring(Bcurve)
        quadratic, remainder = relation.quo_rem(x-x_ring(sx))
        assert not remainder and quadratic.degree() == 2
        discriminant = x_ring.base_ring()(quadratic[1]**2-4*quadratic[2]*quadratic[0])
        numerator = ring(discriminant.numerator())
        denominator_poly = ring(discriminant.denominator())
        odd_degree = sum(
            factor.degree()
            for value in (numerator, denominator_poly)
            for factor, multiplicity in value.squarefree_decomposition()
            if multiplicity % 2
        )
        degree = int(odd_degree+(denominator_poly.degree()-numerator.degree()) % 2)
        histogram[degree] = histogram.get(degree, 0)+1
        if degree == 4:
            good.append(int(level))
    return {
        "histogram": {str(key): value for key, value in sorted(histogram.items())},
        "genus_one_levels": good,
        "singular_levels": singular,
    }


records = []
screen_cache = {}
for a_order in range(args.a_order_min, args.a_order_max+1):
    a_rows, a_cutoff, _ = rows_for(a_values, alpha_den.degree(), a_order)
    for b_order in range(args.b_order_min, args.b_order_max+1):
        b_rows, b_cutoff, _ = rows_for(b_values, beta_den.degree(), b_order)
        condition = matrix(finite, a_rows+b_rows, ncols=len(labels))
        kernel = condition.right_kernel_matrix()
        record = {
            "a_order": a_order,
            "b_order": b_order,
            "rank": int(condition.rank()),
            "kernel_dimension": int(kernel.nrows()),
            "a_cutoff": int(a_cutoff),
            "b_cutoff": int(b_cutoff),
        }
        if args.screen_two_dimensional and kernel.nrows() == 2:
            canonical = kernel.echelon_form()
            key = tuple(tuple(int(value) for value in row) for row in canonical.rows())
            if key not in screen_cache:
                screen_cache[key] = branch_screen(canonical)
            record["branch_screen"] = screen_cache[key]
        records.append(record)

interesting = [
    record for record in records
    if record["kernel_dimension"] in (1, 2, 3)
]
payload = {
    "schema": "elkies-k3.h92-q6-child-q8-component-nef-infinity-sweep.v1",
    "status": "PASS_BOUNDED_SIGNED_INFINITY_SWEEP",
    "prime": args.prime,
    "ambient": {
        "s_degree": args.max_s_degree,
        "t_degree": args.max_t_degree,
        "dimension": len(labels),
    },
    "order_box": {
        "a": [args.a_order_min, args.a_order_max],
        "b": [args.b_order_min, args.b_order_max],
    },
    "records": records,
    "interesting": interesting,
    "boundary": (
        "This is a bounded modular intersection sweep.  Stable two-dimensional "
        "kernels with branch degree four are evidence for the correct infinity "
        "lattice, not a characteristic-zero line-bundle certificate."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")

print(
    "Q8COMPNEFINF|prime={}|ambient={}|interesting={}|two_dim={}|branch4={}|"
    "status=PASS_BOUNDED_SIGNED_INFINITY_SWEEP".format(
        args.prime,
        len(labels),
        ";".join(
            "a{}_b{}_k{}".format(
                record["a_order"], record["b_order"], record["kernel_dimension"]
            )
            for record in interesting
        ) or "none",
        sum(record["kernel_dimension"] == 2 for record in records),
        ";".join(
            "a{}_b{}_good{}_hist{}".format(
                record["a_order"],
                record["b_order"],
                len(record.get("branch_screen", {}).get("genus_one_levels", [])),
                ",".join(
                    "{}:{}".format(key, value)
                    for key, value in record.get("branch_screen", {}).get("histogram", {}).items()
                ) or "none",
            )
            for record in records
            if record["kernel_dimension"] == 2 and "branch_screen" in record
            and record["branch_screen"]["genus_one_levels"]
        ) or "none",
    ),
    flush=True,
)
