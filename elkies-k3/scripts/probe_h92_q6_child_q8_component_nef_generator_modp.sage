#!/usr/bin/env sage -python
"""Screen the corrected component-nef finite-generator ratio modulo p.

The component-nef finite module in q_regular coordinates has generators

    g1 = q_regular + C,          g2 = M=f_II^2*f_IV^2,

where q_regular=(m-p0)/h-R/Nx and C is the CRT lift of R/Nx modulo M.
For

    V = (T^d*g1 + T^e*g2) / g2

a constant level V=v gives

    m = p0 + h * ( M*(v-T^e)/T^d - C + R/Nx ).

The +R/Nx term is essential.  Older reconnaissance omitted it while claiming
to test q_regular.  Fibrewise translation by the old zero is an isomorphism
over the old base, so the branch divisor is the same for the physical
component-nef pencil.
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
MARKING = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-marking.json"
NORMALIZER = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-q-pole-normalization-crt.json"
DEFAULT_FINITE = ROOT / "artifacts/local/elkies-k3/q8-qregular-finite-component-nef.json"
DEFAULT_OUTPUT = ROOT / "artifacts/local/elkies-k3/q8-component-nef-generator-screen.json"


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


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=43)
parser.add_argument("--finite", type=Path, default=DEFAULT_FINITE)
parser.add_argument("--all-v", action="store_true")
parser.add_argument("--v", type=int, default=1)
parser.add_argument("--max-a-degree", type=int, default=0)
parser.add_argument("--max-b-degree", type=int, default=0)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
if not ZZ(args.prime).is_prime() or args.prime in (2, 3):
    raise ValueError("prime must be odd and different from 3")
if args.max_a_degree < 0 or args.max_b_degree < 0:
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
A = polynomial(ring, finite, child["minimal_short_weierstrass"]["A_coefficients_low_to_high"])
Bcurve = polynomial(ring, finite, child["minimal_short_weierstrass"]["B_coefficients_low_to_high"])
section = marking["selected_q8"]["relative_child_section_standard_jacobian_coordinates"]
nx = polynomial(ring, finite, section["x_numerator_coefficients_low_to_high"])
dx = polynomial(ring, finite, section["x_denominator_coefficients_low_to_high"])
ny = polynomial(ring, finite, section["y_numerator_coefficients_low_to_high"])
dy = polynomial(ring, finite, section["y_denominator_coefficients_low_to_high"])
sx, sy = field(nx)/field(dx), field(ny)/field(dy)
assert sy**2 == sx**3 + field(A)*sx + field(Bcurve)
h = monic_power_root(dx, 2)
ii = polynomial(ring, finite, PolynomialRing(QQ, "T")(next(
    item for item in child["finite_fibres"] if item["kodaira"] == "II*"
)["factor"]).list())
iv = polynomial(ring, finite, PolynomialRing(QQ, "T")(next(
    item for item in child["finite_fibres"] if item["kodaira"] == "IV*"
)["factor"]).list())
M = ii**2 * iv**2
assert M.degree() == 4 and nx.gcd(M).degree() == 0 and nx.gcd(h).degree() == 0

R = polynomial(
    ring, finite,
    normalizer["normalizer"]["R_coefficients_low_to_high"],
)
assert (R*h*dy-ny) % nx == 0
C = polynomial(
    ring, finite,
    finite_module["module"]["first_C_lift_coefficients_low_to_high"],
)
assert C.degree() < M.degree()
assert (C-R*nx.inverse_mod(M)) % M == 0
p0 = -sy/sx


def branch_at(level, a_degree, b_degree):
    multiplier = field(T**a_degree)
    translation = field(T**b_degree)
    q_regular_value = field(M)*(finite(level)-translation)/multiplier - field(C)
    m_value = p0 + field(h)*(q_regular_value + field(R)/field(nx))

    # Exact reconstruction of the stated finite-generator level.
    q_check = (m_value-p0)/field(h)-field(R)/field(nx)
    assert multiplier*(q_check+field(C))/field(M)+translation == finite(level)

    x_ring = PolynomialRing(field, "x")
    x = x_ring.gen()
    y = x_ring(m_value)*(x-x_ring(sx))-x_ring(sy)
    relation = y**2-x**3-x_ring(A)*x-x_ring(Bcurve)
    quadratic, remainder = relation.quo_rem(x-x_ring(sx))
    assert not remainder and quadratic.degree() == 2
    discriminant = x_ring.base_ring()(quadratic[1]**2-4*quadratic[2]*quadratic[0])
    numerator = ring(discriminant.numerator())
    denominator = ring(discriminant.denominator())
    odd_factors = []
    odd_degree = 0
    for side, value in (("num", numerator), ("den", denominator)):
        for factor, multiplicity in value.squarefree_decomposition():
            if multiplicity % 2:
                odd_degree += factor.degree()
                odd_factors.append([side, int(factor.degree()), int(multiplicity)])
    infinity = (denominator.degree()-numerator.degree()) % 2
    return {
        "level": int(level),
        "branch_degree": int(odd_degree+infinity),
        "finite_odd_degree": int(odd_degree),
        "infinity_branch": int(infinity),
        "odd_factor_degrees": odd_factors,
    }


levels = list(finite) if args.all_v else [finite(args.v)]
records = []
for a_degree in range(args.max_a_degree+1):
    for b_degree in range(args.max_b_degree+1):
        values = [branch_at(level, a_degree, b_degree) for level in levels]
        histogram = {}
        for value in values:
            degree = value["branch_degree"]
            histogram[degree] = histogram.get(degree, 0)+1
        records.append({
            "a_degree": a_degree,
            "b_degree": b_degree,
            "histogram": {str(key): value for key, value in sorted(histogram.items())},
            "genus_one_levels": [
                value["level"] for value in values if value["branch_degree"] == 4
            ],
            "levels": values,
        })

payload = {
    "schema": "elkies-k3.h92-q6-child-q8-component-nef-generator-screen.v1",
    "status": "PASS_CORRECTED_COMPONENT_NEF_GENERATOR_SCREEN",
    "prime": args.prime,
    "finite_module": str(args.finite),
    "ratio": "V=(T^d*(q_regular+C)+T^e*M)/M",
    "correction": "m=p+h*(M*(V-T^e)/T^d-C+R/Nx)",
    "records": records,
    "boundary": (
        "This is a necessary modular branch-degree screen.  Degree four is "
        "required for a genus-one level but does not by itself certify the "
        "global q8 line bundle or characteristic-zero pencil."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q8COMPNEFGEN|prime={}|levels={}|rows={}|screens={}|status=PASS_CORRECTED_COMPONENT_NEF_GENERATOR_SCREEN".format(
        args.prime,
        len(levels),
        ";".join(
            "d{}_e{}:{}:good={}".format(
                row["a_degree"],
                row["b_degree"],
                ",".join("{}:{}".format(key, value) for key, value in row["histogram"].items()),
                ",".join(map(str, row["genus_one_levels"])) or "none",
            )
            for row in records
        ),
        len(records),
    ),
    flush=True,
)
