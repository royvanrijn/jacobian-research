#!/usr/bin/env sage -python
"""Probe the H92 source-q8 divisor in the collision-regular q,X frame.

This avoids the old endpoint ambient's independent h-denominator choices.

Coordinates:
    p = y(P1)/x(P1)
    q = (m-p)/h
    Xs = h^2*x

At the E8 singular chart h is a unit and
    q = u^-2 * J,  J in (X,Y).
Hence the q8 E8 module u^9*(u^2,X,Y) gives
    u^i*q^b: i >= 9+2b  (b>=1),
    u^i*Xs*q^b: i >= 13+2b.

At the generic point of E7_1:
    ord(t)=2, ord(m)=1, ord(p)=2, ord(h)=-8,
so ord(q)=9 and ord(Xs)=-14.  The q8 comparison therefore gives
    u^i*q^b:     i <= floor((9b-16)/2),
    u^i*Xs*q^b:  i <= floor((9b-30)/2).

Intersecting these exact endpoint bounds with
    1,q,...,q^9,Xs,...,Xs*q^7
leaves exactly 11 candidates:
    u^23 q^7
    u^25..u^28 q^8
    u^27..u^32 q^9
and no Xs terms.

The script then evaluates g*f/t^9 on each of the six actual resolved E7
edge-node charts.  It uses the same local-normal-form quotient as the
repository's existing node probe, but clears the q-frame directly.  Since the
largest remaining base pole is only five, the local quotient is (surface,t^5).

A stacked node kernel of dimension two is the expected q8-pencil signal.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix, sage_eval, singular


ROOT = Path.cwd()
P1 = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
PULLBACKS = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-chart-pullbacks.json"
GLUING = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-actual-e7-gluing.json"
DEFAULT_OUTPUT = ROOT / "artifacts/local/elkies-k3/q8-qframe-e7-node-probe.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def invert_base(rational_u):
    numerator, denominator = rational_u.numerator(), rational_u.denominator()
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
    denominator = finite(ZZ(value.denominator()))
    if not denominator:
        raise ValueError("prime divides an input coefficient denominator")
    return finite(ZZ(value.numerator()))/denominator


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=43)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--trace", action="store_true")
args = parser.parse_args()
if not ZZ(args.prime).is_prime() or args.prime in (2, 3):
    raise ValueError("prime must be an odd prime different from 3")
args.output = args.output.resolve()

for path in (P1, PULLBACKS, GLUING):
    if not path.exists():
        raise SystemExit("missing {}".format(path))

p1 = json.loads(P1.read_text())
pullbacks = json.loads(PULLBACKS.read_text())
gluing = json.loads(GLUING.read_text())
assert p1["status"] == "PASS_EXACT_H92_P1"
assert pullbacks["status"] == "PASS_EXACT_H92_E7_CHART_PULLBACKS"
assert gluing["status"] == "PASS_EXACT_Q8_ACTUAL_E7_GLUING"

# -------------------------------------------------------------------------
# Exact q-frame ambient derivation.
# -------------------------------------------------------------------------
# E7 component data already certified by the actual resolved H92 cover.
v_t = (2, 2, 4, 3, 1, 2, 3)
v_x = (2, 4, 6, 4, 2, 3, 5)
v_m = (1, 1, 3, 2, 0, 2, 2)
twist = (2, 5, 6, 4, 6, 3, 5)

# p=y(P1)/x(P1)=t*(unit), so ord_component(p)=ord_component(t).
# Therefore ord(m-p)>=min(ord(m),ord(t)); equality holds whenever unequal.
# h(1/t)=t^-4*(unit).
v_q_lower = tuple(min(v_m[j], v_t[j]) + 4*v_t[j] for j in range(7))
v_Xs = tuple(v_x[j] - 8*v_t[j] for j in range(7))

candidates = []
for b in range(10):
    e8_floor = 11 if b == 0 else 9+2*b
    e71_upper = (9*b-16)//2
    for i in range(e8_floor, e71_upper+1):
        candidates.append({"x_power": 0, "q_power": b, "u_power": i})
for b in range(8):
    e8_floor = 13+2*b
    e71_upper = (9*b-30)//2
    for i in range(e8_floor, e71_upper+1):
        candidates.append({"x_power": 1, "q_power": b, "u_power": i})

expected = (
    [{"x_power": 0, "q_power": 7, "u_power": 23}]
    + [{"x_power": 0, "q_power": 8, "u_power": i} for i in range(25, 29)]
    + [{"x_power": 0, "q_power": 9, "u_power": i} for i in range(27, 33)]
)
assert candidates == expected
assert len(candidates) == 11

# Every candidate satisfies every generic E7 component even using only the
# lower bound for ord(m-p); possible cancellation at E7_6 only helps.
generic_records = []
for entry in candidates:
    a, b, i = entry["x_power"], entry["q_power"], entry["u_power"]
    residuals = []
    for j in range(7):
        value = (
            twist[j] - (i+9)*v_t[j]
            + b*v_q_lower[j] + a*v_Xs[j]
        )
        residuals.append(int(value))
    assert min(residuals) >= 0
    generic_records.append({"candidate": entry, "residual_lower_bounds": residuals})

# The marked E7_5 smooth branch has m=(t/Z)*n with n=unit/W and g/t^6 a
# unit.  Since 1/h=t^4*(unit), the worst n^b term of u^i*q^b has coefficient
# t^(6+4b-i); all lower n powers are less singular in W.
for entry in candidates:
    assert 6+4*entry["q_power"]-entry["u_power"] >= 0

# -------------------------------------------------------------------------
# Common P1 t-expansions used on every actual E7 node chart.
# -------------------------------------------------------------------------
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
r = x_p_t/t_formal**2
s = y_p_t/t_formal**3
assert r.valuation() == 0 and s.valuation() == 0
r_num, r_den = t_ring(r.numerator()), t_ring(r.denominator())
s_num, s_den = t_ring(s.numerator()), t_ring(s.denominator())
assert r_num(0) and r_den(0) and s_num(0) and s_den(0)
h = polynomial(u_ring, p1["structured_denominator"]["Z4_coefficients"])
assert h.degree() == 4
h_reverse = t_ring(list(reversed(h.list())))
assert h_reverse(0)

finite = GF(args.prime)
qq_ring = PolynomialRing(QQ, names=("Z", "U", "Y"))
Zq, Uq, Yq = qq_ring.gens()
ring = PolynomialRing(finite, names=("Z", "U", "Y"))
Z, U, Y = ring.gens()

# For all 11 candidates, max(i+9-4b)=5.
T = max(entry["u_power"]+9-4*entry["q_power"] for entry in candidates)
assert T == 5
BMAX = 9

node_records = []
stacked_rows = []

for chart in sorted(pullbacks["charts"], key=lambda item: item["name"]):
    name = chart["name"]
    edge = next(item for item in gluing["actual_edge_chart_gluing"] if item["name"] == name)
    if args.trace:
        print("Q8QFRAME_NODE|prime={}|chart={}|stage=setup".format(args.prime, name), flush=True)

    surface_qq = qq_ring(sage_eval(
        chart["surface_equation"], locals={"Z": Zq, "U": Uq, "Y": Yq}
    ))
    t_qq = qq_ring(sage_eval(
        chart["old_coordinate_pullback"]["t"], locals={"Z": Zq, "U": Uq, "Y": Yq}
    ))
    x_qq = qq_ring(sage_eval(
        chart["old_coordinate_pullback"]["x"], locals={"Z": Zq, "U": Uq, "Y": Yq}
    ))
    y_qq = qq_ring(sage_eval(
        chart["old_coordinate_pullback"]["y"], locals={"Z": Zq, "U": Uq, "Y": Yq}
    ))
    g_qq = qq_ring(sage_eval(
        edge["w_cartier_equation"], locals={"Z": Zq, "U": Uq, "Y": Yq}
    ))
    surface, t_value, x_value, y_value, g = tuple(
        ring(value) for value in (surface_qq, t_qq, x_qq, y_qq, g_qq)
    )

    def eval_t_poly(value):
        return ring(sum(
            reduce_coefficient(coefficient, finite)*t_value**degree
            for degree, coefficient in enumerate(value.list())
        ))

    R_num, R_den, S_num, S_den, H_reverse = tuple(
        eval_t_poly(value) for value in (r_num, r_den, s_num, s_den, h_reverse)
    )
    assert all(value(0, 0, 0) for value in (R_num, R_den, S_num, S_den, H_reverse))

    # Existing exact chord factorization:
    # m = M_m * B*R_den/(A*S_den).
    numerator_x = x_value*R_den - t_value**2*R_num
    numerator_y = y_value*S_den - t_value**3*S_num
    mx = common_monomial_exponents(numerator_x)
    my = common_monomial_exponents(numerator_y)
    m_exponents = tuple(my[k]-mx[k] for k in range(3))
    assert all(e >= 0 for e in m_exponents)
    mx_value = ring.monomial(*mx)
    my_value = ring.monomial(*my)
    A = numerator_x // mx_value
    B = numerator_y // my_value
    assert A(0, 0, 0)

    m_monomial = ring.monomial(*m_exponents)

    # q=(m-p)/h with
    # p=t*S_num*R_den/(S_den*R_num), h=t^-4*H_reverse:
    #
    # q = t^4 * R_den *
    #     (m_monomial*B*R_num - t*S_num*A)
    #     /(A*S_den*R_num*H_reverse).
    C = m_monomial*B*R_num - t_value*S_num*A
    D = A*S_den*R_num*H_reverse
    assert D(0, 0, 0)

    t_exponents = common_monomial_exponents(t_value)
    t_monomial = ring.monomial(*t_exponents)
    assert (t_value // t_monomial)(0, 0, 0)

    def truncate(value):
        value = ring(value)
        return ring({
            monomial: coefficient
            for monomial, coefficient in value.dict().items()
            if not all(
                monomial[k] >= T*t_exponents[k] for k in range(3)
            )
        })

    def mul(left, right):
        return truncate(truncate(left)*truncate(right))

    def power(value, exponent):
        answer = ring.one()
        value = truncate(value)
        while exponent:
            if exponent & 1:
                answer = mul(answer, value)
            value = mul(value, value)
            exponent //= 2
        return answer

    numerators = []
    for entry in candidates:
        b, i = entry["q_power"], entry["u_power"]
        t_power = T + 4*b - i - 9
        assert t_power >= 0
        # Multiply g*f/t^9 by the common local unit D^9 and by t^T.
        # Regularity is equivalent to this numerator belonging to (t^T).
        value = ring.one()
        for factor in (
            g,
            power(t_value, t_power),
            power(R_den, b),
            power(C, b),
            power(D, BMAX-b),
        ):
            value = mul(value, factor)
        numerators.append(value)

    singular.eval("ring q8qframe={},(Z,U,Y),ds;".format(args.prime))
    singular.eval("poly surface={};".format(surface))
    singular.eval("ideal principal=surface,({})^{};".format(t_monomial, T))
    singular.eval("ideal standard_principal=std(principal);")

    remainders = []
    for index, numerator in enumerate(numerators):
        singular.eval("poly qnum={};".format(numerator))
        remainder = ring(singular("reduce(qnum,standard_principal)").sage())
        remainders.append(remainder)
        if args.trace:
            print(
                "Q8QFRAME_NODE|prime={}|chart={}|candidate={}|terms={}".format(
                    args.prime, name, index, len(remainder.dict())
                ),
                flush=True,
            )

    coordinates = sorted({
        tuple(monomial)
        for remainder in remainders
        for monomial in remainder.dict()
    })
    rows = []
    for monomial in coordinates:
        rows.append([
            remainder.dict().get(monomial, finite.zero())
            for remainder in remainders
        ])
    node_matrix = matrix(finite, rows, ncols=len(candidates))
    node_rank = int(node_matrix.rank())
    node_kernel = len(candidates)-node_rank
    stacked_rows.extend(rows)
    node_records.append({
        "chart": name,
        "rows": len(rows),
        "rank": node_rank,
        "kernel_dimension": node_kernel,
    })
    print(
        "Q8QFRAMENODE|prime={}|chart={}|rows={}|rank={}|kernel={}".format(
            args.prime, name, len(rows), node_rank, node_kernel
        ),
        flush=True,
    )

stacked = matrix(finite, stacked_rows, ncols=len(candidates))
stacked_rank = int(stacked.rank())
kernel = stacked.right_kernel_matrix()
kernel_dimension = int(kernel.nrows())

payload = {
    "schema": "elkies-k3.h92-q8-qframe-e7-node-probe.v1",
    "status": "EXPERIMENTAL_Q8_QFRAME_E7_NODE_PROBE",
    "prime": int(args.prime),
    "inputs": {
        "p1": {"path": str(P1.relative_to(ROOT)), "sha256": digest(P1)},
        "pullbacks": {"path": str(PULLBACKS.relative_to(ROOT)), "sha256": digest(PULLBACKS)},
        "gluing": {"path": str(GLUING.relative_to(ROOT)), "sha256": digest(GLUING)},
    },
    "q_frame": {
        "q": "(m-y(P1)/x(P1))/h",
        "Xs": "h^2*x",
        "generic_basis": "1,q,...,q^9,Xs,...,Xs*q^7",
        "E8_argument": "q=u^-2*J with J in (X,Y); q8 target u^9*(u^2,X,Y)",
        "E7_1_orders": {"t": 2, "q": 9, "Xs": -14, "twist": 2},
    },
    "ambient": {
        "dimension": len(candidates),
        "candidates": candidates,
        "generic_E7_residual_lower_bounds": generic_records,
        "marked_E7_5_worst_coefficient_orders": [
            6+4*entry["q_power"]-entry["u_power"] for entry in candidates
        ],
    },
    "node_local_quotient_t_power": T,
    "nodes": node_records,
    "stacked": {
        "rows": int(stacked.nrows()),
        "rank": stacked_rank,
        "kernel_dimension": kernel_dimension,
        "kernel_basis_rows": [
            [int(value) for value in row] for row in kernel.rows()
        ],
    },
    "boundary": (
        "This is a modular exact-chart E7-node probe in the q-frame. "
        "The E8 ideal argument, generic E7 inequalities and marked-E7 worst-term "
        "bound are encoded, but characteristic-zero lifting and a complete Cech "
        "overlap certificate remain separate."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")

print(
    "Q8QFRAME|prime={}|ambient=11|nodes={}|stacked_rows={}|rank={}|kernel={}|"
    "status=EXPERIMENTAL_Q8_QFRAME_E7_NODE_PROBE".format(
        args.prime, len(node_records), stacked.nrows(), stacked_rank, kernel_dimension
    ),
    flush=True,
)
