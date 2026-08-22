#!/usr/bin/env sage -python
"""Evaluate one exact H92 q8 E7-node principal condition modulo a prime.

The six-node clearing atlas supplies an exact common-unit clearing for the
condition ``g*f/t^9 in R``.  For one chosen actual resolved chart, this script
maps a finite ambient into the infinite quotient ``R/(t^T)`` using local
standard normal forms.  It is a modular local regression, not a replacement
for the characteristic-zero condition matrix or for overlap gluing.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, sage_eval, singular


ROOT = Path(__file__).resolve().parents[2]
P1 = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
PULLBACKS = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-chart-pullbacks.json"
GLUING = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-actual-e7-gluing.json"
CLEARINGS = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e7-node-principal-clearings.json"
AMBIENT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-endpoint-rr-ambient.json"
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def invert_base(rational_u):
    numerator = rational_u.numerator()
    denominator = rational_u.denominator()
    t_ring = PolynomialRing(QQ, "t")
    t = t_ring.gen()
    field = t_ring.fraction_field()
    return field(
        t**(denominator.degree()-numerator.degree())
        * t_ring(list(reversed(numerator.list())))
        / t_ring(list(reversed(denominator.list())))
    )


def common_monomial_exponents(value):
    terms = list(value.dict())
    assert terms
    return tuple(min(exponent[index] for exponent in terms) for index in range(3))


def reduce_coefficient(value, finite):
    value = QQ(value)
    denominator = finite(value.denominator())
    if not denominator:
        raise ValueError("prime divides an input coefficient denominator")
    return finite(value.numerator()) / denominator


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--p1", type=Path, default=P1)
parser.add_argument("--pullbacks", type=Path, default=PULLBACKS)
parser.add_argument("--gluing", type=Path, default=GLUING)
parser.add_argument("--clearings", type=Path, default=CLEARINGS)
parser.add_argument("--ambient", type=Path, default=AMBIENT)
parser.add_argument("--chart", default="E7_4--E7_3")
parser.add_argument("--prime", type=int, default=43)
parser.add_argument("--output", type=Path)
args = parser.parse_args()
if args.prime <= 1:
    raise ValueError("prime must be greater than one")
if args.output is None:
    args.output = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-{}-principal-node-local-normal-form-mod-{}.json".format(
        args.chart.replace("--", "-"), args.prime
    )

p1 = json.loads(args.p1.read_text())
pullbacks = json.loads(args.pullbacks.read_text())
gluing = json.loads(args.gluing.read_text())
clearings = json.loads(args.clearings.read_text())
ambient = json.loads(args.ambient.read_text())
exec(compile(CORE.read_text(), str(CORE), "exec"))
assert p1["status"] == "PASS_EXACT_H92_P1"
assert pullbacks["status"] == "PASS_EXACT_H92_E7_CHART_PULLBACKS"
assert gluing["status"] == "PASS_EXACT_Q8_ACTUAL_E7_GLUING"
assert clearings["status"] == "PASS_EXACT_Q8_E7_NODE_PRINCIPAL_CLEARINGS"
assert ambient["status"] in {
    "PASS_EXACT_Q8_ENDPOINT_RR_AMBIENT",
    "PASS_EXACT_Q8_ENLARGED_ENDPOINT_RR_AMBIENT",
}
node = next((item for item in clearings["nodes"] if item["chart"] == args.chart), None)
if node is None:
    raise ValueError("unknown cleared E7 chart: {}".format(args.chart))
T, K = (int(clearings["common_parameters"][key]) for key in ("T", "K"))
assert K == max(int(entry["h_power"]) for entry in ambient["ambient_basis"])
assert T == 9 + max(
    int(entry["u_power"])-4*int(entry["h_power"])
    for entry in ambient["ambient_basis"]
)

u_ring = PolynomialRing(QQ, "u")
u_field = u_ring.fraction_field()
x_p = u_field(polynomial(u_ring, p1["x_entrance_base"]["numerator_coefficients"]))
x_p /= u_field(polynomial(u_ring, p1["x_entrance_base"]["denominator_coefficients"]))
y_p = u_field(polynomial(u_ring, p1["y_entrance_base"]["numerator_coefficients"]))
y_p /= u_field(polynomial(u_ring, p1["y_entrance_base"]["denominator_coefficients"]))
x_p_t = invert_base(x_p)
y_p_t = invert_base(y_p)
t_ring = x_p_t.parent()
t_formal = t_ring.gen()
r, s = x_p_t/t_formal**2, y_p_t/t_formal**3
assert r.valuation() == 0 and s.valuation() == 0
t_poly = PolynomialRing(QQ, "t")
r_num, r_den = t_poly(r.numerator()), t_poly(r.denominator())
s_num, s_den = t_poly(s.numerator()), t_poly(s.denominator())
h_reverse = t_poly(list(reversed(polynomial(
    u_ring, p1["structured_denominator"]["Z4_coefficients"]
).list())))

chart = next(item for item in pullbacks["charts"] if item["name"] == args.chart)
edge = next(item for item in gluing["actual_edge_chart_gluing"] if item["name"] == args.chart)
qq_ring = PolynomialRing(QQ, names=("Z", "U", "Y"))
Zq, Uq, Yq = qq_ring.gens()
surface_qq = qq_ring(sage_eval(
    chart["surface_equation"], locals={"Z": Zq, "U": Uq, "Y": Yq}
))
t_qq = qq_ring(sage_eval(chart["old_coordinate_pullback"]["t"], locals={"Z": Zq, "U": Uq, "Y": Yq}))
x_qq = qq_ring(sage_eval(chart["old_coordinate_pullback"]["x"], locals={"Z": Zq, "U": Uq, "Y": Yq}))
y_qq = qq_ring(sage_eval(chart["old_coordinate_pullback"]["y"], locals={"Z": Zq, "U": Uq, "Y": Yq}))
g_qq = qq_ring(sage_eval(edge["w_cartier_equation"], locals={"Z": Zq, "U": Uq, "Y": Yq}))
finite = GF(args.prime)
ring = PolynomialRing(finite, names=("Z", "U", "Y"))
Z, U, Y = ring.gens()
surface, t_value, x_value, y_value, g = tuple(
    ring(value) for value in (surface_qq, t_qq, x_qq, y_qq, g_qq)
)


def evaluate_t_polynomial(value):
    return ring(sum(
        reduce_coefficient(coefficient, finite)*t_value**degree
        for degree, coefficient in enumerate(value.list())
    ))


R_num, R_den, S_num, S_den, H_reverse = tuple(
    evaluate_t_polynomial(value)
    for value in (r_num, r_den, s_num, s_den, h_reverse)
)
numerator_x = x_value*R_den-t_value**2*R_num
numerator_y = y_value*S_den-t_value**3*S_num
mx, my = common_monomial_exponents(numerator_x), common_monomial_exponents(numerator_y)
m_exponents = tuple(my[index]-mx[index] for index in range(3))
assert all(exponent >= 0 for exponent in m_exponents)
mx_value, my_value = ring.monomial(*mx), ring.monomial(*my)
A, B = numerator_x // mx_value, numerator_y // my_value
assert A(0, 0, 0) and R_den(0, 0, 0) and S_den(0, 0, 0) and H_reverse(0, 0, 0)
t_exponents = common_monomial_exponents(t_value)
t_monomial = ring.monomial(*t_exponents)
assert (t_value // t_monomial)(0, 0, 0)


def truncate_t_power(value):
    """Discard only terms already in the actual local ideal (t^T)."""
    return ring({
        monomial: coefficient
        for monomial, coefficient in ring(value).dict().items()
        if not all(monomial[index] >= T*t_exponents[index] for index in range(3))
    })


def truncated_product(left, right):
    return truncate_t_power(ring(left)*ring(right))


def truncated_power(value, exponent):
    answer = ring.one()
    value = truncate_t_power(value)
    while exponent:
        if exponent & 1:
            answer = truncated_product(answer, value)
        value = truncated_product(value, value)
        exponent //= 2
    return answer


m_monomial = ring.monomial(*m_exponents)


def common_cleared_numerator(entry):
    a, b, i, k = tuple(int(entry[key]) for key in (
        "x_power", "m_power", "u_power", "h_power"
    ))
    t_exponent = T+4*k-i-9
    assert t_exponent >= 0 and 0 <= b <= 9
    answer = ring.one()
    for factor in (
        g, t_value**t_exponent, x_value**a, m_monomial**b,
        truncated_power(B, b), truncated_power(R_den, b),
        truncated_power(A, 9-b), truncated_power(S_den, 9-b),
        truncated_power(H_reverse, K-k),
    ):
        answer = truncated_product(answer, factor)
    return answer


singular.eval("ring q8node={},(Z,U,Y),ds;".format(args.prime))
singular.eval("poly surface={};".format(surface))
singular.eval("ideal principal=surface,({})^{};".format(t_monomial, T))
singular.eval("ideal standard_principal=std(principal);")

remainders = []
for entry in ambient["ambient_basis"]:
    numerator = common_cleared_numerator(entry)
    singular.eval("poly node_num={};".format(numerator))
    remainder = ring(singular("reduce(node_num,standard_principal)").sage())
    singular.eval("poly remainder={};".format(remainder))
    assert ring(singular("reduce(remainder,standard_principal)").sage()) == remainder
    remainders.append(remainder)

image = finite_ambient_image_condition(
    "actual {} principal local-normal-form image".format(args.chart),
    tuple(range(len(remainders))),
    lambda index: {tuple(monomial): coefficient for monomial, coefficient in remainders[index].dict().items()},
    lambda monomial: monomial,
    finite,
    "actual H92 resolved chart; local normal form modulo (surface,t^T)",
)
matrix = image["matrix"]
monomials = image["coordinate_keys"]
monomial_index = {monomial: index for index, monomial in enumerate(monomials)}
payload = {
    "schema": "elkies-k3.h92-q8-e7-node-principal-local-normal-form-modp.v1",
    "status": "EXPERIMENTAL_MODULAR_Q8_E7_NODE_LOCAL_NORMAL_FORM_BLOCK",
    "prime": args.prime,
    "inputs": {
        "checker_source": {"path": str(Path(__file__).relative_to(ROOT)), "sha256": digest(Path(__file__))},
        "p1": {"path": str(args.p1.relative_to(ROOT)), "sha256": digest(args.p1)},
        "actual_pullbacks": {"path": str(args.pullbacks.relative_to(ROOT)), "sha256": digest(args.pullbacks)},
        "q8_gluing": {"path": str(args.gluing.relative_to(ROOT)), "sha256": digest(args.gluing)},
        "node_clearings": {"path": str(args.clearings.relative_to(ROOT)), "sha256": digest(args.clearings)},
        "endpoint_ambient": {"path": str(args.ambient.relative_to(ROOT)), "sha256": digest(args.ambient)},
        "compiler_core": {"path": str(CORE.relative_to(ROOT)), "sha256": digest(CORE)},
    },
    "local_ring": {
        "chart": args.chart,
        "order": "Singular ds local degree order at (Z,U,Y)",
        "surface_equation": str(surface),
        "principal_ideal": "(surface,({})^{})".format(t_monomial, T),
        "t_is_monomial_times_unit": True,
    },
    "finite_ambient_image": {
        "ambient_dimension": len(remainders), "rows": len(monomials),
        "rank": int(matrix.rank()), "kernel_dimension": int(matrix.right_kernel().dimension()),
        "coordinate_monomials": ["Z^{}*U^{}*Y^{}".format(*monomial) for monomial in monomials],
        "sparse_columns": [
            [[monomial_index[tuple(monomial)], int(coefficient)] for monomial, coefficient in sorted(remainder.dict().items())]
            for remainder in remainders
        ],
    },
    "boundary": (
        "This is one modular chart image in the infinite local quotient. It "
        "does not certify a characteristic-zero q8 matrix, overlap gluing, "
        "a common q8 kernel, h0(D), a pencil, or a child equation."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q8E7NODELOCALNF|chart={}|prime={}|ambient={}|rows={}|rank={}|kernel={}|"
    "status=EXPERIMENTAL_MODULAR_Q8_E7_NODE_LOCAL_NORMAL_FORM_BLOCK".format(
        args.chart, args.prime, len(remainders), len(monomials), matrix.rank(),
        matrix.right_kernel().dimension(),
    ),
    flush=True,
)
