#!/usr/bin/env sage -python
"""Screen the first normalized finite-module generator as a q8 base coordinate.

In the globally base-regular frame q_regular=q-R/Nx, the first finite module
generator is f_IV*q_regular+C, where C is the degree-less-than-five lift of
f_IV*R/Nx modulo M=f_II^2*f_IV^3.  The ratio

    V=(f_IV*q_regular+C)/M

gives m=p+h*(M*V-C)/f_IV.  This script eliminates the old Weierstrass x over
F_p(T,V) and checks the necessary genus-one branch degree.  It is only a
modular screen for this one canonical finite-generator ratio, not a q8 pencil.
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
MARKING = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-marking.json"


def coefficient(field, value):
    value = QQ(value)
    denominator = field(ZZ(value.denominator()))
    if not denominator:
        raise ValueError("prime divides an input denominator")
    return field(ZZ(value.numerator())) / denominator


def polynomial(ring, field, coefficients):
    return ring([coefficient(field, value) for value in coefficients])


def monic_power_root(value, exponent):
    root = value.parent().one()
    for factor, multiplicity in value.factor():
        assert multiplicity % exponent == 0
        root *= factor.monic() ** (multiplicity // exponent)
    return root.monic()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=43)
parser.add_argument("--v", type=int, default=1)
parser.add_argument("--all-v", action="store_true", help="screen every constant level in F_p")
parser.add_argument(
    "--max-a-monomial-degree",
    type=int,
    help="also scan a=T^d for 0 <= d <= this bound in V=(a*g1)/M",
)
parser.add_argument(
    "--max-b-monomial-degree",
    type=int,
    help="scan b=T^e in V=(a*g1+b*M)/M (requires --max-a-monomial-degree)",
)
args = parser.parse_args()
if not ZZ(args.prime).is_prime() or args.prime in (2, 3):
    raise ValueError("prime must be odd and different from 3")

child = json.loads(CHILD.read_text())
marking = json.loads(MARKING.read_text())
finite = GF(args.prime)
ring = PolynomialRing(finite, "T")
T = ring.gen()
field = ring.fraction_field()
A = polynomial(ring, finite, child["minimal_short_weierstrass"]["A_coefficients_low_to_high"])
B = polynomial(ring, finite, child["minimal_short_weierstrass"]["B_coefficients_low_to_high"])
section = marking["selected_q8"]["relative_child_section_standard_jacobian_coordinates"]
nx = polynomial(ring, finite, section["x_numerator_coefficients_low_to_high"])
dx = polynomial(ring, finite, section["x_denominator_coefficients_low_to_high"])
ny = polynomial(ring, finite, section["y_numerator_coefficients_low_to_high"])
dy = polynomial(ring, finite, section["y_denominator_coefficients_low_to_high"])
sx, sy = field(nx)/field(dx), field(ny)/field(dy)
assert sy**2 == sx**3 + field(A)*sx + field(B)
h = monic_power_root(dx, 2)
ii = polynomial(ring, finite, PolynomialRing(QQ, "T")(next(
    item for item in child["finite_fibres"] if item["kodaira"] == "II*"
)["factor"]).list())
iv = polynomial(ring, finite, PolynomialRing(QQ, "T")(next(
    item for item in child["finite_fibres"] if item["kodaira"] == "IV*"
)["factor"]).list())
M = ii**2 * iv**3
assert nx.gcd(M).degree() == 0 and nx.gcd(h).degree() == 0
R = (ny * (h*dy).inverse_mod(nx)).mod(nx)
C = (iv * R * nx.inverse_mod(M)).mod(M)
assert (C - iv*R*nx.inverse_mod(M)) % M == 0

def screen(level, monomial_degree=0, translation_degree=None):
    """Return the exact branch data of the double cover at one V level."""

    # A level V=v is a double cover of the old base.  Its branch degree must
    # be four for the candidate level curve to have genus one.
    multiplier = field(T)**monomial_degree
    translation = field.zero() if translation_degree is None else field(T)**translation_degree
    m = -sy/sx + field(h) * (field(M)*(finite(level)-translation)/multiplier-field(C)) / field(iv)
    x_ring = PolynomialRing(field, "x")
    x = x_ring.gen()
    y = x_ring(m) * (x-x_ring(sx)) - x_ring(sy)
    relation = y**2-x**3-x_ring(A)*x-x_ring(B)
    quadratic, remainder = relation.quo_rem(x-x_ring(sx))
    assert remainder == 0 and quadratic.degree() == 2
    a, b, c = quadratic[2], quadratic[1], quadratic[0]
    discriminant = x_ring.base_ring()(b**2-4*a*c)
    numerator = ring(discriminant.numerator())
    denominator = ring(discriminant.denominator())
    odd_degree = sum(
        factor.degree()
        for polynomial_value in (numerator, denominator)
        for factor, multiplicity in polynomial_value.squarefree_decomposition()
        if multiplicity % 2
    )
    infinity = (denominator.degree()-numerator.degree()) % 2
    return odd_degree + infinity, odd_degree, infinity


if args.max_b_monomial_degree is not None and args.max_a_monomial_degree is None:
    raise ValueError("max b monomial degree requires max a monomial degree")
levels = range(args.prime) if args.all_v or args.max_a_monomial_degree is not None else (args.v,)


def summarize(monomial_degree, translation_degree=None):
    results = [(level, *screen(level, monomial_degree, translation_degree)) for level in levels]
    by_degree = {}
    for _, branch_degree, _, _ in results:
        by_degree[branch_degree] = by_degree.get(branch_degree, 0)+1
    return results, by_degree


if args.max_a_monomial_degree is not None:
    if args.max_a_monomial_degree < 0 or (
        args.max_b_monomial_degree is not None and args.max_b_monomial_degree < 0
    ):
        raise ValueError("monomial degrees must be nonnegative")
    rows = []
    for degree in range(args.max_a_monomial_degree+1):
        translations = (
            (None,) if args.max_b_monomial_degree is None
            else tuple(range(args.max_b_monomial_degree+1))
        )
        for translation_degree in translations:
            results, by_degree = summarize(degree, translation_degree)
            rows.append((degree, translation_degree, by_degree, [level for level, branch, _, _ in results if branch == 4]))
    print(
        "H92Q6CHILDQREGMONOMIALMODP|prime={}|levels={}|max_a_degree={}|"
        "rows={}|status={}".format(
            args.prime, len(levels), args.max_a_monomial_degree,
            ";".join(
                "d{}_e{}:{}:good={}".format(
                    degree,
                    "none" if translation_degree is None else translation_degree,
                    ",".join("{}:{}".format(branch, count) for branch, count in sorted(histogram.items())),
                    ",".join(str(level) for level in good_levels) or "none",
                )
                for degree, translation_degree, histogram, good_levels in rows
            ),
            "PASS_MONOMIAL_LEVEL_SCREEN",
        ),
        flush=True,
    )
elif args.all_v:
    results, by_degree = summarize(0)
    print(
        "H92Q6CHILDQREGGENMODP|prime={}|levels={}|C_degree={}|"
        "branch_degree_histogram={}|genus_one_levels={}|status={}".format(
            args.prime, len(results), C.degree(),
            ",".join("{}:{}".format(degree, count) for degree, count in sorted(by_degree.items())),
            ",".join(str(level) for level, degree, _, _ in results if degree == 4) or "none",
            "PASS_ALL_LEVELS_NECESSARY_SCREEN"
        ),
        flush=True,
    )
else:
    level = args.v
    branch_degree, odd_degree, infinity = screen(level)
    print(
        "H92Q6CHILDQREGGENMODP|prime={}|V={}|C_degree={}|branch_degree={}|"
        "finite_odd_degree={}|infinity_branch={}|status={}".format(
            args.prime, level, C.degree(), branch_degree, odd_degree, infinity,
            "PASS_NECESSARY_GENUS_ONE" if branch_degree == 4 else "REJECTED_BRANCH_DEGREE",
        ),
        flush=True,
    )
