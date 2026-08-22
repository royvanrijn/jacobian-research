#!/usr/bin/env sage -python
"""Screen the first two-dimensional saturated q8-child coefficient slice.

The exact smooth module writes an element of the marked generic fibre space as

    f=a+b*m,       a=A/h^2, b=B/h,

with ``A*D+B*N=0 mod h^2``.  This script forms a bounded polynomial span of
such pairs, imposes the exact II*/IV* jet rows modulo a good prime, and, when
its kernel is two-dimensional, tests the ratio of its two deterministic kernel
basis elements.  At every constant level of that ratio, eliminating ``y``
leaves a quadratic in ``x``.  The odd squarefree degree of its discriminant is
the branch degree over the old base; degree four is the necessary genus-one
value.

This is a finite-field reconnaissance of one bounded coefficient slice.  It
does not establish the complete q8 divisor, a characteristic-zero pencil, or
a rootless equation.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix


ROOT = Path(__file__).resolve().parents[2]
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
MARKING = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-marking.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-saturated-pencil-mod-43.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficient(field, value):
    value = QQ(value)
    denominator = field(ZZ(value.denominator()))
    if not denominator:
        raise ValueError("the chosen prime divides an input coefficient denominator")
    return field(ZZ(value.numerator())) / denominator


def polynomial(ring, field, values):
    return ring([coefficient(field, value) for value in values])


def rational(function_field, ring, field, data, numerator, denominator):
    return function_field(polynomial(ring, field, data[numerator])) / function_field(
        polynomial(ring, field, data[denominator])
    )


def reduce_rational(value, modulus):
    ring = modulus.parent()
    numerator, denominator = ring(value.numerator()), ring(value.denominator())
    assert denominator.gcd(modulus).degree() == 0
    return (numerator * denominator.inverse_mod(modulus)).mod(modulus)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=43)
parser.add_argument("--max-b-degree", type=int, default=7)
parser.add_argument("--max-c-degree", type=int, default=-1)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
if not ZZ(args.prime).is_prime() or args.prime in (2, 3):
    raise ValueError("prime must be an odd prime different from 3")
if args.max_b_degree < 0 or args.max_c_degree < -1:
    raise ValueError("max-b-degree must be nonnegative and max-c-degree at least -1")

child = json.loads(CHILD.read_text())
marking = json.loads(MARKING.read_text())
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert marking["status"] == "PASS_EXACT_Q6_CHILD_Q8_MARKING"

finite = GF(args.prime)
ring = PolynomialRing(finite, "T")
T = ring.gen()
field = ring.fraction_field()
section = marking["selected_q8"]["relative_child_section_standard_jacobian_coordinates"]
x_s = rational(field, ring, finite, section, "x_numerator_coefficients_low_to_high", "x_denominator_coefficients_low_to_high")
y_s = rational(field, ring, finite, section, "y_numerator_coefficients_low_to_high", "y_denominator_coefficients_low_to_high")
A_curve = polynomial(ring, finite, child["minimal_short_weierstrass"]["A_coefficients_low_to_high"])
B_curve = polynomial(ring, finite, child["minimal_short_weierstrass"]["B_coefficients_low_to_high"])
assert y_s**2 == x_s**3 + field(A_curve)*x_s + field(B_curve)

h = ring.one()
for factor, multiplicity in ring(x_s.denominator()).factor():
    assert multiplicity % 2 == 0
    h *= factor.monic()**(multiplicity // 2)
h = h.monic()
assert h.degree() == 46
p = -y_s/x_s
N, p_denominator = ring(p.numerator()), ring(p.denominator())
D, remainder = p_denominator.quo_rem(h)
assert not remainder and D.gcd(h**2).degree() == 0
inverse_D = D.inverse_mod(h**2)

labels = [("B", degree) for degree in range(args.max_b_degree+1)]
labels += [("C", degree) for degree in range(args.max_c_degree+1)]
pairs = []
for kind, degree in labels:
    if kind == "B":
        B = T**degree
        A = (-B*N*inverse_D).mod(h**2)
    else:
        B = ring.zero()
        A = h**2*T**degree
    assert (A*D+B*N).mod(h**2) == 0
    pairs.append((A, B))


def additive_rows(kodaira):
    fibre = next(item for item in child["finite_fibres"] if item["kodaira"] == kodaira)
    source_ring = PolynomialRing(QQ, "T")
    factor = polynomial(ring, finite, source_ring(fibre["factor"]).list())
    point = -factor[0]/factor[1]
    u_ring = PolynomialRing(finite, "u")
    u = u_ring.gen()
    u_field = u_ring.fraction_field()

    def translate(value):
        return u_field(u_ring(value.numerator()(point+u))) / u_field(
            u_ring(value.denominator()(point+u))
        )

    a_values = [translate(field(A)/field(h**2)) for A, _ in pairs]
    b_values = [translate(field(B)/field(h)) for _, B in pairs]
    m_constant = translate(-y_s/x_s)
    if kodaira == "II*":
        return [
            [reduce_rational(a_values[index]+b_values[index]*m_constant, u**2)[jet]
             for index in range(len(labels))]
            for jet in range(2)
        ]
    b_u = u_ring(B_curve(point+u))
    unit_b = b_u // u**4
    if not unit_b(0).is_square():
        raise ValueError("the selected prime does not split the IV* branch constant")
    c = unit_b(0).sqrt()
    x_u, y_u = translate(x_s), translate(y_s)
    m_u = reduce_rational(m_constant+c/x_u(0)*u**2, u**3)
    x_coefficient = -y_u(0)/x_u(0)**2
    return [
        [reduce_rational(a_values[index]+b_values[index]*m_u, u**3)[jet]
         for index in range(len(labels))]
        for jet in range(3)
    ] + [[reduce_rational(b_values[index], u)[0]*x_coefficient for index in range(len(labels))]]


condition = matrix(finite, additive_rows("II*")+additive_rows("IV*"))
kernel = condition.right_kernel().basis_matrix()
if kernel.nrows() != 2:
    raise ValueError(
        "the declared bounded slice has kernel dimension {}, not two".format(kernel.nrows())
    )


def combine(row):
    return (
        sum((row[index]*pairs[index][0] for index in range(len(labels))), ring.zero()),
        sum((row[index]*pairs[index][1] for index in range(len(labels))), ring.zero()),
    )


(A0, B0), (A1, B1) = (combine(row) for row in kernel.rows())


def branch_degree(level):
    alpha = field(A0-level*A1)/field(h**2)
    beta = field(B0-level*B1)/field(h)
    if not beta:
        return None
    slope = -alpha/beta
    x_ring = PolynomialRing(field, "x")
    x = x_ring.gen()
    y = x_ring(slope)*(x-x_ring(x_s))-x_ring(y_s)
    relation = y**2-x**3-x_ring(A_curve)*x-x_ring(B_curve)
    quadratic, remainder = relation.quo_rem(x-x_ring(x_s))
    assert not remainder and quadratic.degree() == 2
    discriminant = x_ring.base_ring()(quadratic[1]**2-4*quadratic[2]*quadratic[0])
    numerator, denominator = ring(discriminant.numerator()), ring(discriminant.denominator())
    odd_degree = sum(
        factor.degree()
        for polynomial_value in (numerator, denominator)
        for factor, multiplicity in polynomial_value.squarefree_decomposition()
        if multiplicity % 2
    )
    infinity = (denominator.degree()-numerator.degree()) % 2
    return int(odd_degree+infinity)


levels = {int(level): branch_degree(finite(level)) for level in finite}
histogram = {}
for degree in levels.values():
    histogram[str(degree)] = histogram.get(str(degree), 0)+1
good_levels = sorted(level for level, degree in levels.items() if degree == 4)
payload = {
    "schema": "elkies-k3.h92-q6-child-q8-saturated-pencil-modp.v1",
    "status": "EXPERIMENTAL_BOUNDED_SATURATED_PENCIL_SCREEN",
    "inputs": {"child": digest(CHILD), "marking": digest(MARKING)},
    "prime": int(args.prime),
    "ansatz": {
        "coefficient_form": "a=A/h^2, b=B/h",
        "smooth_congruence": "A*D+B*N=0 mod h^2",
        "B_degree_at_most": args.max_b_degree,
        "C_degree_at_most": args.max_c_degree,
        "dimension": len(labels),
    },
    "additive_kernel": {
        "condition_rank": int(condition.rank()),
        "dimension": int(kernel.nrows()),
        "basis_labels": labels,
        "basis_rows": [[int(value) for value in row] for row in kernel.rows()],
    },
    "pencil_ratio": "(a0+b0*m)/(a1+b1*m)",
    "constant_level_branch_degree_histogram": histogram,
    "genus_one_levels": good_levels,
    "boundary": (
        "This is a one-prime screen of one bounded smooth-plus-additive coefficient slice. "
        "It does not prove that the slice is the q8 Riemann--Roch pencil, lift a "
        "surviving level to QQ, construct a D13/rootless equation, or give a bisection collision."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q6CHILDQ8SATPENCIL|prime={}|Bdeg={}|Cdeg={}|ambient={}|rank={}|kernel=2|"
    "branch_histogram={}|genus_one_levels={}|status=EXPERIMENTAL_BOUNDED_SATURATED_PENCIL_SCREEN".format(
        args.prime, args.max_b_degree, args.max_c_degree, len(labels), condition.rank(),
        ",".join("{}:{}".format(key, value) for key, value in sorted(histogram.items())),
        ",".join(map(str, good_levels)) or "none",
    ), flush=True,
)
