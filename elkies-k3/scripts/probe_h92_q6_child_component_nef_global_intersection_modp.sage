#!/usr/bin/env sage -python
"""Exact-bound modular global intersection for the component-nef q8 divisor.

Finite local work gives the fractional submodule

    f = s(T)*(q_phys-rho)/M + t(T)/L,
    L=f_II*f_IV, M=L^2,

where rho is the exact transported q_regular residue modulo M.  Fibrewise
translation by -P0 sends q_phys to q_std, so the old-base double cover and the
infinity coefficient conditions may be computed in the standard q frame.

At smooth infinity

    q_std = alpha + beta*m_inf,
    alpha=-p/h-R/Nx, beta=T^2/h,

and a section is

    a + b*m_inf
    a = s*(alpha-rho)/M + t/L
    b = s*beta/M.

The component-nef vertical fibre coefficient is -2, so the exact local
line-bundle condition is ord_inf(a)>=2 and ord_inf(b)>=2.

The degree bounds are consequences, not search cutoffs:
  * ord(beta/M)=48, hence deg(s)<=46;
  * ord((alpha-rho)/M)=1 because deg(rho)=3;
    cancellation with t/L (order 2-deg(t)) gives deg(t)<=47.
Any coefficient above those bounds has a unique uncancellable leading term.

Thus the 95-dimensional polynomial ambient below is complete for this
fractional finite submodule.  If its infinity kernel has dimension two, the
two rows are explicit global sections of the declared divisor.  The optional
branch screen then checks the genus-one consequence.
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix


ROOT = Path(__file__).resolve().parents[2]
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
MARKING = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-marking.json"
NORMALIZER = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-q-pole-normalization-crt.json"
DEFAULT_MODULE = ROOT / "artifacts/local/elkies-k3/q8-component-nef-qreg-finite-module.json"
DEFAULT_OUTPUT = ROOT / "artifacts/local/elkies-k3/q8-component-nef-global-intersection-modp.json"


def coefficient(finite, value):
    value = QQ(value)
    den = finite(ZZ(value.denominator()))
    if not den:
        raise ValueError("prime divides an input denominator")
    return finite(ZZ(value.numerator())) / den


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
parser.add_argument("--module", type=Path, default=DEFAULT_MODULE)
parser.add_argument("--screen-levels", action="store_true")
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
if not ZZ(args.prime).is_prime() or args.prime in (2, 3):
    raise ValueError("prime must be odd and different from 3")
args.module = args.module.resolve()
args.output = args.output.resolve()

child = json.loads(CHILD.read_text())
marking = json.loads(MARKING.read_text())
normalizer = json.loads(NORMALIZER.read_text())
module = json.loads(args.module.read_text())
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert marking["status"] == "PASS_EXACT_Q6_CHILD_Q8_MARKING"
assert normalizer["status"] == "PASS_EXACT_CRT_PRINCIPAL_PART_NORMALIZATION"
assert module["status"] == "PASS_EXACT_COMPONENT_NEF_QREG_FINITE_MODULE"

finite = GF(args.prime)
ring = PolynomialRing(finite, "T")
T = ring.gen()
field = ring.fraction_field()

sdata = marking["selected_q8"]["relative_child_section_standard_jacobian_coordinates"]
nx = polynomial(ring, finite, sdata["x_numerator_coefficients_low_to_high"])
dx = polynomial(ring, finite, sdata["x_denominator_coefficients_low_to_high"])
ny = polynomial(ring, finite, sdata["y_numerator_coefficients_low_to_high"])
dy = polynomial(ring, finite, sdata["y_denominator_coefficients_low_to_high"])
sx, sy = field(nx) / field(dx), field(ny) / field(dy)
h = monic_power_root(dx, 2)
assert dx // h**2 in finite and dy // h**3 in finite

R = polynomial(
    ring, finite,
    normalizer["normalizer"]["R_coefficients_low_to_high"],
)
assert (R*h*dy-ny) % nx == 0
p = -sy/sx

rho = polynomial(
    ring, finite,
    module["module"]["rho_coefficients_low_to_high"],
)
ii = polynomial(ring, finite, PolynomialRing(QQ, "T")(
    next(item for item in child["finite_fibres"] if item["kodaira"] == "II*")["factor"]
).list()).monic()
iv = polynomial(ring, finite, PolynomialRing(QQ, "T")(
    next(item for item in child["finite_fibres"] if item["kodaira"] == "IV*")["factor"]
).list()).monic()
L = (ii*iv).monic()
M = L**2
assert L.degree() == 2 and M.degree() == 4 and rho.degree() == 3

alpha = -p/field(h) - field(R)/field(nx)
beta = field(T**2)/field(h)
abase = (alpha-field(rho))/field(M)
bbase = beta/field(M)
tbase = field.one()/field(L)

# Exact infinity orders before polynomial multipliers.
def infinity_order(value):
    numerator = ring(value.numerator())
    denominator = ring(value.denominator())
    return int(denominator.degree()-numerator.degree())

assert infinity_order(bbase) == 48
assert infinity_order(abase) == 1
assert infinity_order(tbase) == 2

max_s_degree = 46
max_t_degree = 47
required_order = 2

labels = [("s", exponent) for exponent in range(max_s_degree+1)]
labels += [("t", exponent) for exponent in range(max_t_degree+1)]
assert len(labels) == 95

# Put all columns into a common-denominator high-coefficient test separately
# for a and b.
a_columns = []
b_columns = []
for kind, exponent in labels:
    if kind == "s":
        a_columns.append(field(T**exponent)*abase)
        b_columns.append(field(T**exponent)*bbase)
    else:
        a_columns.append(field(T**exponent)*tbase)
        b_columns.append(field.zero())

def common_denominator(values):
    result = ring.one()
    for value in values:
        denominator = ring(value.denominator())
        result = result.lcm(denominator)
    return result.monic()

a_den = common_denominator(a_columns)
b_den = common_denominator([value for value in b_columns if value] or [field.one()])
a_nums = [
    ring(value*a_den) for value in a_columns
]
b_nums = [
    ring(value*b_den) if value else ring.zero()
    for value in b_columns
]
assert all(field(num)/field(a_den) == value for num, value in zip(a_nums, a_columns))
assert all(
    (not value and not num) or field(num)/field(b_den) == value
    for num, value in zip(b_nums, b_columns)
)

def high_rows(numerators, denominator, required):
    cutoff = denominator.degree()-required
    top = max((degree_or_minus_one(value) for value in numerators), default=-1)
    rows = [
        [
            value[degree] if value and degree <= value.degree() else finite.zero()
            for value in numerators
        ]
        for degree in range(cutoff+1, top+1)
    ]
    return rows, int(cutoff), int(top)

a_rows, a_cutoff, a_top = high_rows(a_nums, a_den, required_order)
b_rows, b_cutoff, b_top = high_rows(b_nums, b_den, required_order)
condition = matrix(finite, a_rows+b_rows, ncols=len(labels))
kernel = condition.right_kernel_matrix()

def st_from_row(row):
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
    return ring(s_poly), ring(t_poly)

def coefficient_pair(row):
    s_poly, t_poly = st_from_row(row)
    # Ignore the common 1/M for ratios.  Numerator section:
    # n=s*(q-rho)+t*L = A+B*m.
    Bcoef = field(s_poly)/field(h)
    Acoef = (
        field(t_poly*L)
        - field(s_poly)*p/field(h)
        - field(s_poly*R)/field(nx)
        - field(s_poly*rho)
    )
    return s_poly, t_poly, Acoef, Bcoef

kernel_records = []
for row in kernel.rows():
    s_poly, t_poly, Acoef, Bcoef = coefficient_pair(row)
    actual_a = (
        field(s_poly)*abase + field(t_poly)*tbase
    )
    actual_b = field(s_poly)*bbase
    assert infinity_order(actual_a) >= 2
    assert not actual_b or infinity_order(actual_b) >= 2
    kernel_records.append({
        "s_degree": -1 if not s_poly else int(s_poly.degree()),
        "t_degree": -1 if not t_poly else int(t_poly.degree()),
        "s_coefficients_low_to_high": [int(value) for value in s_poly.list()],
        "t_coefficients_low_to_high": [int(value) for value in t_poly.list()],
        "a_infinity_order": infinity_order(actual_a),
        "b_infinity_order": None if not actual_b else infinity_order(actual_b),
    })

branch = None
if args.screen_levels and kernel.nrows() == 2:
    _, _, A0, B0 = coefficient_pair(kernel.row(0))
    _, _, A1, B1 = coefficient_pair(kernel.row(1))
    Acurve = polynomial(
        ring, finite,
        child["minimal_short_weierstrass"]["A_coefficients_low_to_high"],
    )
    Bcurve = polynomial(
        ring, finite,
        child["minimal_short_weierstrass"]["B_coefficients_low_to_high"],
    )
    xring = PolynomialRing(field, "x")
    x = xring.gen()
    histogram = {}
    good = []
    singular = []
    records = []
    for level in finite:
        denominator = B1-field(level)*B0
        if not denominator:
            singular.append(int(level))
            continue
        mvalue = -(A1-field(level)*A0)/denominator
        y = xring(mvalue)*(x-xring(sx))-xring(sy)
        relation = y**2-x**3-xring(Acurve)*x-xring(Bcurve)
        quadratic, remainder = relation.quo_rem(x-xring(sx))
        assert not remainder and quadratic.degree() == 2
        discriminant = xring.base_ring()(
            quadratic[1]**2-4*quadratic[2]*quadratic[0]
        )
        num = ring(discriminant.numerator())
        den = ring(discriminant.denominator())
        odd_degree = 0
        odd_factors = []
        for side, value in (("n", num), ("d", den)):
            for factor, multiplicity in value.squarefree_decomposition():
                if multiplicity % 2:
                    odd_degree += factor.degree()
                    odd_factors.append(
                        [side, int(factor.degree()), int(multiplicity)]
                    )
        infinity_branch = (den.degree()-num.degree()) % 2
        branch_degree = int(odd_degree+infinity_branch)
        histogram[branch_degree] = histogram.get(branch_degree, 0)+1
        if branch_degree == 4:
            good.append(int(level))
        records.append({
            "level": int(level),
            "branch_degree": branch_degree,
            "finite_odd_degree": int(odd_degree),
            "infinity_branch": int(infinity_branch),
            "disc_num_degree": int(num.degree()),
            "disc_den_degree": int(den.degree()),
            "odd_factor_degrees": odd_factors,
        })
    branch = {
        "histogram": {str(key): value for key, value in sorted(histogram.items())},
        "genus_one_levels": good,
        "singular_levels": singular,
        "records": records,
    }

payload = {
    "schema": "elkies-k3.h92-q6-child-q8-component-nef-global-intersection-modp.v1",
    "status": "PASS_EXACT_BOUND_COMPONENT_NEF_GLOBAL_INTERSECTION_MODP",
    "prime": args.prime,
    "finite_fractional_module": {
        "form": "s*(q_phys-rho)/M+t/L",
        "L_degree": 2,
        "M_degree": 4,
        "rho_degree": 3,
    },
    "infinity": {
        "required_coefficient_order": 2,
        "base_orders": {
            "(alpha-rho)/M": infinity_order(abase),
            "beta/M": infinity_order(bbase),
            "1/L": infinity_order(tbase),
        },
        "derived_degree_bounds": {
            "s": max_s_degree,
            "t": max_t_degree,
            "proof": (
                "deg(s)>46 violates the m_inf coefficient uniquely; "
                "after deg(s)<=46, deg(t)>47 gives an uncancellable "
                "leading constant coefficient."
            ),
        },
        "a_rows": len(a_rows),
        "b_rows": len(b_rows),
        "a_cutoff": a_cutoff,
        "b_cutoff": b_cutoff,
    },
    "result": {
        "ambient_dimension": len(labels),
        "rank": int(condition.rank()),
        "kernel_dimension": int(kernel.nrows()),
        "kernel": kernel_records,
    },
    "branch_screen": branch,
    "boundary": (
        "This is a complete modular infinity intersection inside the exact "
        "fractional finite submodule.  A two-dimensional kernel whose levels "
        "have branch degree four is the expected equation-level q8 handoff; "
        "characteristic-zero reconstruction remains separate."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q8GLOBAL|prime={}|ambient=95|rank={}|kernel={}|degrees={}|hist={}|good={}|"
    "status=PASS_EXACT_BOUND_COMPONENT_NEF_GLOBAL_INTERSECTION_MODP".format(
        args.prime,
        condition.rank(),
        kernel.nrows(),
        ";".join(
            "s{}_t{}_a{}_b{}".format(
                item["s_degree"], item["t_degree"],
                item["a_infinity_order"],
                "inf" if item["b_infinity_order"] is None else item["b_infinity_order"],
            )
            for item in kernel_records
        ) or "none",
        "none" if branch is None else ",".join(
            "{}:{}".format(key, value)
            for key, value in branch["histogram"].items()
        ),
        "none" if branch is None or not branch["genus_one_levels"] else ",".join(
            map(str, branch["genus_one_levels"])
        ),
    ),
    flush=True,
)
