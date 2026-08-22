#!/usr/bin/env sage
"""Insert the reconstructed coefficient curve into the q=80 K3 chart.

This is deliberately an *unmarked* certificate.  It proves that the exact
rational functions reconstructed from the slope 8/87 formal branch give a
one-parameter family in the E6+D5+A3 ambient Weierstrass chart.  It does not
claim that the three local marked sections algebraize globally, nor that the
resulting coefficient line has already been identified with Elkies's
X(6,79)/<w_474> quotient.

The family is kept compact: the output records the four rational coefficient
functions by hash/reference and the universal formulas for A and B, instead
of expanding the same very large coefficients a second time.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


parser = argparse.ArgumentParser()
parser.add_argument(
    "--parameter",
    default="artifacts/generated-results/q80-cm24-slope-8-87-qq-PDQE-parameter.json",
)
parser.add_argument(
    "--output",
    default="artifacts/generated-results/q80-cm24-slope-8-87-qq-surface-family.json",
)
arguments = parser.parse_args()


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


parameter_path = Path(arguments.parameter)
payload = json.loads(parameter_path.read_text())
if payload.get("schema") != "q80-cm24-qq-PDQE-parameter-v1":
    raise ValueError("unexpected P,D,Q,E parameter schema")

parameter_ring = PolynomialRing(QQ, "u")
u = parameter_ring.gen()
parameter_field = parameter_ring.fraction_field()


def read_function(name):
    record = payload["original_functions"][name]
    numerator = parameter_ring(record["numerator"])
    denominator = parameter_ring(record["denominator"])
    if denominator == 0 or numerator.gcd(denominator).degree() != 0:
        raise ArithmeticError(f"{name}(u) is not a reduced rational function")
    if [int(numerator.degree()), int(denominator.degree())] != record["degrees"]:
        raise ArithmeticError(f"stored degree mismatch for {name}(u)")
    return parameter_field(numerator/denominator), numerator, denominator


function_data = {name: read_function(name) for name in ("d", "p", "q", "e")}
d_u, p_u, q_u, e_u = (function_data[name][0] for name in ("d", "p", "q", "e"))

# Derive the universal q=80 chart independently, rather than trusting copied
# formulas.  The first four jets at T=1 are the cubic binomial truncation of
# the multiplicative branch.  Hermite interpolation then gives b1,...,b4.
parameters = PolynomialRing(QQ, names=("d", "p", "q", "e"))
d, p, q, e = parameters.gens()
coefficient_field = parameters.fraction_field()
universal_T_ring = PolynomialRing(coefficient_field, "T")
T = universal_T_ring.gen()
jet_ring = PolynomialRing(coefficient_field, "s")
s = jet_ring.gen()

r = -3*d**2+3-p-q
A = T**2*(-3+p*T+q*T**2+r*T**3)
A_at_one = jet_ring(A(T=1+s))
assert A_at_one[0] == -3*d**2
correction = (A_at_one+3*d**2)/(-3*d**2)
branch = 2*d**3*(
    1+QQ(3)/2*correction+QQ(3)/8*correction**2-QQ(1)/16*correction**3
)
branch_jets = vector(coefficient_field, [branch[index] for index in range(4)])
jet_matrix = Matrix(
    coefficient_field,
    4,
    4,
    lambda row, column: jet_ring((1+s)**(4+column))[row],
)
fixed_jets = vector(
    coefficient_field,
    [jet_ring(2*(1+s)**3+e*(1+s)**8)[index] for index in range(4)],
)
b_coefficients = tuple(jet_matrix.solve_right(branch_jets-fixed_jets))
b1, b2, b3, b4 = b_coefficients
B = T**3*(2+b1*T+b2*T**2+b3*T**3+b4*T**4+e*T**5)

discriminant = 4*A**3+27*B**2
fixed_discriminant_factor = T**7*(T-1)**4
residual, remainder = discriminant.quo_rem(fixed_discriminant_factor)
assert remainder == 0 and residual.degree() == 5
zero_open = (discriminant//T**7)(T=0)
one_open = (discriminant//(T-1)**4)(T=1)
infinity_open = discriminant[16]
assert zero_open == 108*(p+b1)
assert infinity_open == 27*e**2
assert one_open != 0


def specialize_fraction(value):
    """Evaluate an element of QQ(d,p,q,e) in QQ(u)."""
    numerator = value.numerator()
    denominator = value.denominator()
    substitutions = {d: d_u, p: p_u, q: q_u, e: e_u}
    result_denominator = denominator.subs(substitutions)
    if result_denominator == 0:
        raise ZeroDivisionError("coefficient curve lies in a chart denominator")
    return parameter_field(numerator.subs(substitutions)/result_denominator)


specialized_open = {
    "T=0": specialize_fraction(zero_open),
    "T=1": specialize_fraction(one_open),
    "T=infinity": specialize_fraction(infinity_open),
}
if any(value == 0 for value in specialized_open.values()):
    raise ArithmeticError("the coefficient curve is contained in a Kodaira boundary")

# Locate the CM(-24) point intrinsically.  In centered coordinates it is
# (D,P,Q,E)=(0,0,0,0), equivalently
# (d,p,q,e)=(-1/2,9/4,-9/4,-27/32).  Cross multiplication avoids evaluating
# the enormous functions until the unique common linear factor is known.
cm24_values = {
    "d": -QQ(1)/2,
    "p": QQ(9)/4,
    "q": -QQ(9)/4,
    "e": -QQ(27)/32,
}
cm24_equations = []
for name, expected in cm24_values.items():
    _, numerator, denominator = function_data[name]
    cm24_equations.append(numerator-expected*denominator)
cm24_factor = cm24_equations[0]
for equation in cm24_equations[1:]:
    cm24_factor = cm24_factor.gcd(equation)
cm24_factor = cm24_factor.monic()
if cm24_factor.degree() != 1:
    raise ArithmeticError(
        f"expected a unique CM24 preimage, found gcd degree {cm24_factor.degree()}"
    )
cm24_parameter = -cm24_factor[0]/cm24_factor[1]
for name, expected in cm24_values.items():
    value, _, _ = function_data[name]
    assert value(cm24_parameter) == expected

# The universal identity plus a well-defined substitution is an exact proof
# of the discriminant shape after specialization.  Record independent facts
# showing that the substitution is a genuine curve and remains in the stated
# generic Kodaira chart.
function_degrees = {
    name: [int(part.degree()) for part in function_data[name][1:]]
    for name in ("d", "p", "q", "e")
}
assert any(numerator.degree() or denominator.degree() for numerator, denominator in (
    function_data[name][1:] for name in ("d", "p", "q", "e")
))

output = {
    "schema": "q80-qq-unmarked-surface-family-v1",
    "scope": "exact_unmarked_one_parameter_surface_family",
    "status": "PASS_EXACT_UNMARKED_SURFACE_FAMILY",
    "parameter": "u",
    "ambient_chart": {
        "equation": "y^2=x^3+A(u,T)*x+B(u,T)",
        "A": str(A),
        "B": str(B),
        "r": str(r),
        "b_coefficients": [str(value) for value in b_coefficients],
        "discriminant": "4*A^3+27*B^2=T^7*(T-1)^4*R5(T)",
        "residual_degree": int(5),
        "generic_fibers": ["I1* at T=0", "I4 at T=1", "IV* at T=infinity"],
        "generic_root_lattice": "D5+A3+E6",
    },
    "coefficient_functions": {
        "source": str(parameter_path),
        "source_sha256": sha256(parameter_path),
        "names": ["d", "p", "q", "e"],
        "degrees": function_degrees,
    },
    "cm24_anchor": {
        "parameter": str(cm24_parameter),
        "minimal_polynomial": str(cm24_factor),
        "values": {name: str(value) for name, value in cm24_values.items()},
    },
    "checks": {
        "universal_discriminant_remainder": "0",
        "universal_residual_degree": int(5),
        "specialized_open_coefficients_nonzero": {
            place: True for place in specialized_open
        },
        "cm24_specialization": "exact",
    },
    "claim_boundary": {
        "proved": [
            "exact rational coefficient functions",
            "exact inclusion in the q=80 E6+D5+A3 Weierstrass chart",
            "exact CM(-24) specialization",
        ],
        "not_proved": [
            "global marked Mordell-Weil sections",
            "generic Mordell-Weil rank three",
            "identification of the base or a cover with X(6,79)/<w_474>",
        ],
    },
}
output_path = Path(arguments.output)
output_path.parent.mkdir(parents=True, exist_ok=True)
output_path.write_text(json.dumps(output, indent=2, sort_keys=True)+"\n")

print(
    "Q80FAMILY|functions=d:{},p:{},q:{},e:{}|"
    "Delta=T^7*(T-1)^4*R5|R5_degree=5|"
    "fibers=I1*,I4,IV*|cm24_preimages=1|output={}|"
    "status=PASS_EXACT_UNMARKED_SURFACE_FAMILY".format(
        "/".join(map(str, function_degrees["d"])),
        "/".join(map(str, function_degrees["p"])),
        "/".join(map(str, function_degrees["q"])),
        "/".join(map(str, function_degrees["e"])),
        output_path,
    ),
    flush=True,
)
