#!/usr/bin/env sage -python
"""Modular branch screen for the exact component-nef fractional q8 pencil.

The exact module predicts the new base coordinate

    U = (m_component_nef-rho)/L,

where L=f_II*f_IV and rho interpolates the singular-point residues of the
physical chord.  After fibrewise translation by -P0 this becomes the standard
marked chord equation

    m_S = rho + L*U.

Translation preserves the degree-two cover of the old base, so we can compute
its branch squareclass without expanding the physical group-law pullback.
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
MARKING = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-marking.json"
DEFAULT_MODULE = ROOT / "artifacts/local/elkies-k3/q8-component-nef-fractional-module.json"
DEFAULT_OUTPUT = ROOT / "artifacts/local/elkies-k3/q8-component-nef-fractional-branch-modp.json"


def coefficient(finite, value):
    value = QQ(value)
    den = finite(ZZ(value.denominator()))
    if not den:
        raise ValueError("prime divides denominator")
    return finite(ZZ(value.numerator()))/den


def polynomial(ring, finite, coefficients):
    return ring([coefficient(finite, value) for value in coefficients])


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=43)
parser.add_argument("--module", type=Path, default=DEFAULT_MODULE)
parser.add_argument("--all-v", action="store_true")
parser.add_argument("--v", type=int, default=1)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
if not ZZ(args.prime).is_prime() or args.prime in (2,3):
    raise ValueError("prime must be odd and !=3")
args.module = args.module.resolve()
args.output = args.output.resolve()

child = json.loads(CHILD.read_text())
marking = json.loads(MARKING.read_text())
module = json.loads(args.module.read_text())
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert marking["status"] == "PASS_EXACT_Q6_CHILD_Q8_MARKING"
assert module["status"] == "PASS_EXACT_COMPONENT_NEF_FRACTIONAL_MODULE"

finite = GF(args.prime)
ring = PolynomialRing(finite, "T")
T = ring.gen()
field = ring.fraction_field()
A = polynomial(ring, finite, child["minimal_short_weierstrass"]["A_coefficients_low_to_high"])
Bcurve = polynomial(ring, finite, child["minimal_short_weierstrass"]["B_coefficients_low_to_high"])
sdata = marking["selected_q8"]["relative_child_section_standard_jacobian_coordinates"]
sx = field(polynomial(ring, finite, sdata["x_numerator_coefficients_low_to_high"])) / field(
    polynomial(ring, finite, sdata["x_denominator_coefficients_low_to_high"])
)
sy = field(polynomial(ring, finite, sdata["y_numerator_coefficients_low_to_high"])) / field(
    polynomial(ring, finite, sdata["y_denominator_coefficients_low_to_high"])
)
assert sy**2 == sx**3+field(A)*sx+field(Bcurve)

rho = polynomial(
    ring, finite,
    module["base_polynomials"]["rho_coefficients_low_to_high"],
)
ii = polynomial(ring, finite, PolynomialRing(QQ, "T")(
    next(item for item in child["finite_fibres"] if item["kodaira"]=="II*")["factor"]
).list()).monic()
iv = polynomial(ring, finite, PolynomialRing(QQ, "T")(
    next(item for item in child["finite_fibres"] if item["kodaira"]=="IV*")["factor"]
).list()).monic()
L = (ii*iv).monic()
assert L.degree()==2 and rho.degree()<2

xring = PolynomialRing(field, "x")
x = xring.gen()


def branch(level):
    m = field(rho + L*finite(level))
    y = xring(m)*(x-xring(sx))-xring(sy)
    relation = y**2-x**3-xring(A)*x-xring(Bcurve)
    quadratic, remainder = relation.quo_rem(x-xring(sx))
    assert not remainder and quadratic.degree()==2
    disc = xring.base_ring()(quadratic[1]**2-4*quadratic[2]*quadratic[0])
    num = ring(disc.numerator())
    den = ring(disc.denominator())
    odd_records=[]
    odd_degree=0
    for side, value in (("n",num),("d",den)):
        for factor, multiplicity in value.squarefree_decomposition():
            if multiplicity % 2:
                odd_degree += factor.degree()
                odd_records.append([side,int(factor.degree()),int(multiplicity)])
    infinity=(den.degree()-num.degree())%2
    return {
        "level": int(level),
        "branch_degree": int(odd_degree+infinity),
        "finite_odd_degree": int(odd_degree),
        "infinity_branch": int(infinity),
        "odd_factor_degrees": odd_records,
        "disc_num_degree": int(num.degree()),
        "disc_den_degree": int(den.degree()),
    }


levels = list(finite) if args.all_v else [finite(args.v)]
records=[branch(level) for level in levels]
hist={}
for record in records:
    d=record["branch_degree"]
    hist[d]=hist.get(d,0)+1
good=[record["level"] for record in records if record["branch_degree"]==4]

payload={
    "schema":"elkies-k3.h92-q6-child-q8-component-nef-fractional-branch-modp.v1",
    "status":"PASS_COMPONENT_NEF_FRACTIONAL_BRANCH_SCREEN",
    "prime":args.prime,
    "new_base":"U=(m_component_nef-rho)/L",
    "translated_equation":"m_S=rho+L*U",
    "histogram":{str(key):value for key,value in sorted(hist.items())},
    "genus_one_levels":good,
    "records":records,
    "boundary":"Necessary modular branch screen; characteristic-zero reconstruction is separate.",
}
args.output.parent.mkdir(parents=True,exist_ok=True)
args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(
    "Q8FRACBRANCH|prime={}|levels={}|hist={}|good={}|"
    "status=PASS_COMPONENT_NEF_FRACTIONAL_BRANCH_SCREEN".format(
        args.prime,len(levels),
        ",".join("{}:{}".format(k,v) for k,v in sorted(hist.items())),
        ",".join(map(str,good)) or "none",
    ),
    flush=True,
)
