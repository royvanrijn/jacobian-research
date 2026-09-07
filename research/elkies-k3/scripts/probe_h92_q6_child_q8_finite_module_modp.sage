#!/usr/bin/env sage -python
"""Compute the complete finite q=8-child coefficient module modulo a prime.

This is the modular counterpart of the exact finite-saturation calculation.
It uses the smooth parametrization ``A=A0*B+h^2*C`` and imposes all two II*
and four oriented IV* quotient rows.  Unlike the earlier bounded B-window,
it calculates the entire A=F_p[T] submodule of pairs (B,C) satisfying those
finite conditions, then reports its Smith profile.  Infinity is not imposed.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix


ROOT = Path(__file__).resolve().parents[2]
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
MARKING = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-marking.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficient(field, value):
    value = QQ(value)
    denominator = field(ZZ(value.denominator()))
    if not denominator:
        raise ValueError("the chosen prime divides an input coefficient denominator")
    return field(ZZ(value.numerator())) / denominator


def polynomial(ring, field, coefficients):
    return ring([coefficient(field, value) for value in coefficients])


def rational(field, ring, coefficient_field, data, numerator, denominator):
    return field(polynomial(ring, coefficient_field, data[numerator])) / field(
        polynomial(ring, coefficient_field, data[denominator])
    )


def reduce_rational(value, modulus):
    ring = modulus.parent()
    numerator, denominator = ring(value.numerator()), ring(value.denominator())
    assert denominator.gcd(modulus).degree() == 0
    return (numerator * denominator.inverse_mod(modulus)).mod(modulus)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=43)
parser.add_argument("--output", type=Path, default=None)
args = parser.parse_args()
if not ZZ(args.prime).is_prime() or args.prime in (2, 3):
    raise ValueError("prime must be an odd prime different from 3")
if args.output is None:
    args.output = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-finite-module-mod-{}.json".format(args.prime)

child = json.loads(CHILD.read_text())
marking = json.loads(MARKING.read_text())
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert marking["status"] == "PASS_EXACT_Q6_CHILD_Q8_MARKING"

finite = GF(args.prime)
ring = PolynomialRing(finite, "T")
T = ring.gen()
field = ring.fraction_field()
section = marking["selected_q8"]["relative_child_section_standard_jacobian_coordinates"]
sx = rational(field, ring, finite, section, "x_numerator_coefficients_low_to_high", "x_denominator_coefficients_low_to_high")
sy = rational(field, ring, finite, section, "y_numerator_coefficients_low_to_high", "y_denominator_coefficients_low_to_high")
h = ring.one()
for factor, multiplicity in ring(sx.denominator()).factor():
    assert multiplicity % 2 == 0
    h *= factor.monic() ** (multiplicity // 2)
h = h.monic()
assert h.degree() == 46
p = -sy/sx
N, p_denominator = ring(p.numerator()), ring(p.denominator())
D, remainder = p_denominator.quo_rem(h)
assert not remainder and D.gcd(h**2).degree() == 0
A0 = (-N*D.inverse_mod(h**2)).mod(h**2)
assert (A0*D+N).mod(h**2) == 0


def local_condition(kodaira):
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

    def condition(B, C):
        A = A0*B+h**2*C
        a_value = translate(field(A)/field(h**2))
        b_value = translate(field(B)/field(h))
        m_value = translate(-sy/sx)
        if kodaira == "II*":
            modulus = u**2
            return tuple(reduce_rational(a_value+b_value*m_value, modulus)[jet]
                         for jet in range(2))

        B_curve = polynomial(ring, finite, child["minimal_short_weierstrass"]["B_coefficients_low_to_high"])
        b_curve = u_ring(B_curve(point+u))
        unit_b = b_curve // u**4
        c = unit_b(0).sqrt()
        if not c or c**2 != unit_b(0):
            raise ValueError("the chosen prime does not split the selected IV* branch")
        x_u, y_u = translate(sx), translate(sy)
        modulus = u**3
        m_u = reduce_rational(m_value+c/x_u(0)*u**2, modulus)
        x_coefficient = -y_u(0)/x_u(0)**2
        return tuple(
            [reduce_rational(a_value+b_value*m_u, modulus)[jet] for jet in range(3)]
            + [reduce_rational(b_value, u)[0]*x_coefficient]
        )

    return factor, condition


ii_factor, ii_condition = local_condition("II*")
iv_factor, iv_condition = local_condition("IV*")
modulus = ii_factor**2*iv_factor**3
assert modulus.degree() == 5 and ii_factor.gcd(iv_factor).degree() == 0
pairs = [(T**degree, ring.zero()) for degree in range(modulus.degree())]
pairs += [(ring.zero(), T**degree) for degree in range(modulus.degree())]
rows = []
for condition in (ii_condition, iv_condition):
    values = [condition(B, C) for B, C in pairs]
    rows.extend([[value[row] for value in values] for row in range(len(values[0]))])
condition_matrix = matrix(finite, rows)
assert condition_matrix.dimensions() == (6, 10) and condition_matrix.rank() == 6
kernel = condition_matrix.right_kernel().basis_matrix()
assert kernel.dimensions() == (4, 10)
generators = []
for vector in kernel.rows():
    generators.append((
        sum((vector[degree]*T**degree for degree in range(5)), ring.zero()),
        sum((vector[5+degree]*T**degree for degree in range(5)), ring.zero()),
    ))
generators += [(modulus, ring.zero()), (ring.zero(), modulus)]
generator_matrix = matrix(ring, 2, 6, lambda row, column: generators[column][row])
smith, smith_left, smith_right = generator_matrix.smith_form()
assert smith == smith_left*generator_matrix*smith_right and smith.rank() == 2
diagonal = tuple(ring(smith[index, index]) for index in range(2))
module_basis = smith_left.inverse()*smith[:, :2]
determinant = ring(module_basis.det())
assert determinant.degree() == 6 == sum(value.degree() for value in diagonal)

payload = {
    "schema": "elkies-k3.h92-q6-child-q8-finite-module-modp.v1",
    "status": "EXPERIMENTAL_FINITE_MODULE_MODULAR",
    "inputs": {"child": digest(CHILD), "marking": digest(MARKING)},
    "prime": int(args.prime),
    "coefficient_frame": "a=A/h^2, b=B/h, A=A0*B+h^2*C",
    "finite_conditions": {
        "CRT_modulus_degree": int(modulus.degree()),
        "rows": 6, "rank": int(condition_matrix.rank()),
        "residue_ambient_dimension": 10, "residue_kernel_dimension": 4,
    },
    "module": {
        "rank": 2,
        "smith_degrees": [int(value.degree()) for value in diagonal],
        "determinant_degree": int(determinant.degree()),
        "basis_column_B_C_degrees": [
            [int(ring(module_basis[row, column]).degree()) for row in range(2)]
            for column in range(2)
        ],
    },
    "boundary": "The finite module is complete at this prime, but no infinity condition or characteristic-zero q8 pencil is asserted.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q6CHILDQ8FINITEMODULE|prime={}|finite_rows=6|finite_rank=6|"
    "smith_degrees={}|determinant_degree=6|status=EXPERIMENTAL_FINITE_MODULE_MODULAR".format(
        args.prime, ",".join(str(value.degree()) for value in diagonal)
    ), flush=True,
)
