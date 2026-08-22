#!/usr/bin/env sage -python
"""Probe a bounded saturated global ansatz for the marked q=8 child chord.

This is explicitly reconnaissance, not a q=8 pencil certificate.  The
smooth collision module writes a prospective coefficient pair as

    a=A/h^2, b=B/h,  A*D+B*N=0 mod h^2,

where D=den(p)/h and N=num(p).  For a declared polynomial B bound and an
optional h^2*C correction, this script imposes the *exact* II* and IV*
quotient conditions on a+b*m after reduction at a declared good prime.  It
reports the common-kernel dimension.

The selected defaults use the first eight B monomials and no C correction.
They are a small diagnostic window only; their dimensions are not derived
from the q=8 divisor and therefore a two-dimensional kernel is not a proof.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix


ROOT = Path(__file__).resolve().parents[2]
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
MARKING = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-marking.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-saturated-ansatz-probe.json"


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


def rational(function_field, ring, coefficient_field, data, numerator, denominator):
    return function_field(polynomial(ring, coefficient_field, data[numerator])) / function_field(
        polynomial(ring, coefficient_field, data[denominator])
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

field = GF(args.prime)
T_ring = PolynomialRing(field, "T")
T = T_ring.gen()
T_field = T_ring.fraction_field()
section = marking["selected_q8"]["relative_child_section_standard_jacobian_coordinates"]
x_s = rational(T_field, T_ring, field, section, "x_numerator_coefficients_low_to_high", "x_denominator_coefficients_low_to_high")
y_s = rational(T_field, T_ring, field, section, "y_numerator_coefficients_low_to_high", "y_denominator_coefficients_low_to_high")
h = T_ring.one()
for factor, multiplicity in T_ring(x_s.denominator()).factor():
    assert multiplicity % 2 == 0
    h *= factor.monic()**(multiplicity // 2)
h = h.monic()
assert h.degree() == 46
p = -y_s / x_s
N, p_denominator = T_ring(p.numerator()), T_ring(p.denominator())
D, remainder = p_denominator.quo_rem(h)
assert not remainder and D.gcd(h**2).degree() == 0
inverse_D = D.inverse_mod(h**2)

# Each displayed label determines (A,B), with the smooth congruence built in.
labels = [("B", exponent) for exponent in range(args.max_b_degree + 1)]
labels += [("C", exponent) for exponent in range(args.max_c_degree + 1)]
pairs = []
for kind, exponent in labels:
    if kind == "B":
        B = T**exponent
        A = (-B * N * inverse_D).mod(h**2)
    else:
        B = T_ring.zero()
        A = h**2 * T**exponent
    assert (A * D + B * N).mod(h**2) == 0
    pairs.append((A, B))


def additive_rows(kodaira, order):
    fibre = next(item for item in child["finite_fibres"] if item["kodaira"] == kodaira)
    source_ring = PolynomialRing(QQ, "T")
    factor = polynomial(T_ring, field, source_ring(fibre["factor"]).list())
    A_curve = polynomial(T_ring, field, child["minimal_short_weierstrass"]["A_coefficients_low_to_high"])
    B_curve = polynomial(T_ring, field, child["minimal_short_weierstrass"]["B_coefficients_low_to_high"])

    def valuation(value):
        result = 0
        while value and value % factor == 0:
            value //= factor
            result += 1
        return result

    expected = {"II*": (4, 5, 10), "IV*": (3, 4, 8)}[kodaira]
    assert (valuation(A_curve), valuation(B_curve), valuation(4*A_curve**3 + 27*B_curve**2)) == expected
    point = -factor[0] / factor[1]
    u_ring = PolynomialRing(field, "u")
    u = u_ring.gen()
    u_field = u_ring.fraction_field()

    def translate(value):
        return u_field(u_ring(value.numerator()(point + u))) / u_field(u_ring(value.denominator()(point + u)))

    h_u = u_ring(h(point + u))
    a_values = []
    b_values = []
    for A, B in pairs:
        a_values.append(translate(T_field(A) / T_field(h**2)))
        b_values.append(translate(T_field(B) / T_field(h)))
    m_constant = translate(-y_s / x_s)

    if kodaira == "II*":
        modulus = u**2
        rows = []
        for jet in range(2):
            rows.append([
                reduce_rational(a_values[index] + b_values[index] * m_constant, modulus)[jet]
                for index in range(len(labels))
            ])
        return rows

    assert kodaira == "IV*" and order == 3
    b_u = u_ring(B_curve(point + u))
    unit_b = b_u // u**4
    if not unit_b(0).is_square():
        raise ValueError("the chosen prime does not split the IV* branch constant")
    c = unit_b(0).sqrt()
    assert c and c**2 == unit_b(0)
    x_u, y_u = translate(x_s), translate(y_s)
    m_u = reduce_rational(m_constant + c / x_u(0) * u**2, u**3)
    x_coefficient = -y_u(0) / x_u(0)**2
    rows = []
    for jet in range(3):
        rows.append([
            reduce_rational(a_values[index] + b_values[index] * m_u, u**3)[jet]
            for index in range(len(labels))
        ])
    rows.append([reduce_rational(b_values[index], u)[0] * x_coefficient for index in range(len(labels))])
    return rows


ii_rows = additive_rows("II*", 2)
iv_rows = additive_rows("IV*", 3)
condition = matrix(field, ii_rows + iv_rows)
kernel = condition.right_kernel()
payload = {
    "schema": "elkies-k3.h92-q6-child-q8-saturated-ansatz-probe.v1",
    "status": "EXPERIMENTAL_BOUNDED_SATURATED_ANSATZ",
    "inputs": {"child": digest(CHILD), "marking": digest(MARKING)},
    "prime": int(args.prime),
    "ansatz": {
        "coefficient_form": "a=A/h^2, b=B/h",
        "smooth_congruence": "A*D+B*N=0 mod h^2",
        "B_degree_at_most": args.max_b_degree,
        "C_degree_at_most": args.max_c_degree,
        "A_correction": "A may be shifted by h^2*C",
        "dimension": len(labels),
    },
    "additive_conditions": {"II*_rows": len(ii_rows), "IV*_rows": len(iv_rows), "rank": int(condition.rank())},
    "kernel_dimension": int(kernel.dimension()),
    "boundary": "The bounds are exploratory and do not give a complete q8 global coefficient cover or a pencil certificate.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("H92Q6CHILDQ8PROBE|prime={}|Bdeg={}|Cdeg={}|ambient={}|additive_rank={}|kernel={}|status=EXPERIMENTAL_BOUNDED_SATURATED_ANSATZ".format(args.prime, args.max_b_degree, args.max_c_degree, len(labels), condition.rank(), kernel.dimension()), flush=True)
