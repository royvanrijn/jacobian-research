#!/usr/bin/env sage -python
"""Compute the actual q8 E7_4--E7_3 principal-node block modulo a prime.

This is deliberately not a rectangular jet calculation.  In the completed
resolved chart, q8 regularity of the common-cleared numerator is membership in
the principal ideal ``(t^T)`` (with ``T=17``), and ``R/(t^T)`` is infinite.
For a *fixed finite ambient*, however, its image in that infinite quotient is
a finite-dimensional vector space.  This script represents precisely that
finite image by local standard normal forms in the actual chart local ring.

The computation is modular and therefore a regression/probe, not a
characteristic-zero q8 condition certificate.  It nevertheless uses the
actual H92 blow-up equation, its actual q6 marked-module trivialization, and
the exact principal product ideal; it never substitutes ``(Z^51,Y^68)``.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, matrix, sage_eval, singular


ROOT = Path(__file__).resolve().parents[2]
P1 = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
PULLBACKS = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-chart-pullbacks.json"
AMBIENT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-endpoint-rr-ambient.json"
CLEARING = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e7-4-3-principal-node-clearing.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e7-4-3-principal-node-local-normal-form-mod-43.json"
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def invert_base(rational_u):
    """Rewrite a QQ(u) function in QQ(t), where t=1/u at the old E7 fibre."""
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


def parse_prime(value):
    value = int(value)
    if value <= 1:
        raise ValueError("prime must be greater than one")
    return value


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--p1", type=Path, default=P1)
parser.add_argument("--pullbacks", type=Path, default=PULLBACKS)
parser.add_argument("--ambient", type=Path, default=AMBIENT)
parser.add_argument("--clearing", type=Path, default=CLEARING)
parser.add_argument("--prime", type=parse_prime, default=43)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

p1 = json.loads(args.p1.read_text())
pullbacks = json.loads(args.pullbacks.read_text())
ambient = json.loads(args.ambient.read_text())
clearing = json.loads(args.clearing.read_text())
exec(compile(CORE.read_text(), str(CORE), "exec"))
assert p1["status"] == "PASS_EXACT_H92_P1"
assert pullbacks["status"] == "PASS_EXACT_H92_E7_CHART_PULLBACKS"
assert ambient["status"] in {
    "PASS_EXACT_Q8_ENDPOINT_RR_AMBIENT",
    "PASS_EXACT_Q8_ENLARGED_ENDPOINT_RR_AMBIENT",
}
assert clearing["status"] == "PASS_EXACT_Q8_E7_4_3_PRINCIPAL_NODE_CLEARING"
assert clearing["chart"]["name"] == "E7_4--E7_3"
T = int(clearing["common_clearing"]["T"])
K = max(int(entry["h_power"]) for entry in ambient["ambient_basis"])
# Enlarging every endpoint denominator power raises K but leaves
# max(i-4*k), and hence the actual principal exponent T, unchanged.
assert T == 9+max(
    int(entry["u_power"])-4*int(entry["h_power"])
    for entry in ambient["ambient_basis"]
)
assert T == 17

chart = next(item for item in pullbacks["charts"] if item["name"] == "E7_4--E7_3")
qq_ring = PolynomialRing(QQ, names=("Z", "U", "Y"))
Zq, Uq, Yq = qq_ring.gens()
surface_qq = qq_ring(sage_eval(
    chart["surface_equation"], locals={"Z": Zq, "U": Uq, "Y": Yq}
))
H_qq, remainder = (Yq**2-surface_qq).quo_rem(Uq)
assert not remainder and H_qq(0, 0, 0) == 1

u_ring = PolynomialRing(QQ, "u")
u_field = u_ring.fraction_field()
x_p = u_field(polynomial(u_ring, p1["x_entrance_base"]["numerator_coefficients"]))
x_p /= u_field(polynomial(u_ring, p1["x_entrance_base"]["denominator_coefficients"]))
y_p = u_field(polynomial(u_ring, p1["y_entrance_base"]["numerator_coefficients"]))
y_p /= u_field(polynomial(u_ring, p1["y_entrance_base"]["denominator_coefficients"]))
x_p_t = invert_base(x_p)
y_p_t = invert_base(y_p)
t_ring = PolynomialRing(QQ, "t")
# The local chord units are the entrance series with their prescribed
# t^2/t^3 orders removed.  The clearing certificate calls these r and s;
# using x(P1), y(P1) themselves here would insert two/three spurious powers
# of t into A_num and B_num and no longer evaluate that certificate.
t_formal = x_p_t.parent().gen()
r = x_p_t / t_formal**2
s = y_p_t / t_formal**3
assert r.valuation() == 0 and s.valuation() == 0
r_num, r_den = t_ring(r.numerator()), t_ring(r.denominator())
s_num, s_den = t_ring(s.numerator()), t_ring(s.denominator())
h_reverse = t_ring(list(reversed(polynomial(
    u_ring, p1["structured_denominator"]["Z4_coefficients"]
).list())))

finite = GF(args.prime)
ring = PolynomialRing(finite, names=("Z", "U", "Y"))
Z, U, Y = ring.gens()
surface = ring(surface_qq)
H = ring(H_qq)
t_value = Z**3*U**2
x_value = Z**4*U**3


def reduce_coefficient(value):
    value = QQ(value)
    denominator = finite(value.denominator())
    if not denominator:
        raise ValueError("prime divides an H92 coefficient denominator")
    return finite(value.numerator()) / denominator


def truncate_t_power(value):
    """Reduce modulo the actual principal monomial t^T=Z^(3T)U^(2T).

    This does *not* introduce a finite corner quotient.  It merely discards
    terms already zero in the actual principal ideal before costly local
    standard reduction, an exact operation in the quotient under study.
    """
    return ring({
        monomial: coefficient
        for monomial, coefficient in ring(value).dict().items()
        if not (monomial[0] >= 3*T and monomial[1] >= 2*T)
    })


def truncated_product(left, right):
    return truncate_t_power(ring(left)*ring(right))


def truncated_power(value, exponent):
    result = ring.one()
    value = truncate_t_power(value)
    exponent = int(exponent)
    while exponent:
        if exponent & 1:
            result = truncated_product(result, value)
        value = truncated_product(value, value)
        exponent //= 2
    return result


def evaluate_t_polynomial(value):
    """Evaluate and truncate a polynomial in t=Z^3 U^2 exactly modulo t^T."""
    return truncate_t_power(sum(
        reduce_coefficient(coefficient)*t_value**degree
        for degree, coefficient in enumerate(value.list())
        if degree < T
    ))


R_num, R_den, S_num, S_den, H_reverse = tuple(
    evaluate_t_polynomial(value)
    for value in (r_num, r_den, s_num, s_den, h_reverse)
)
A_num = truncate_t_power(R_den-Z**2*U*R_num)
B_num = truncate_t_power(H*S_den-Z**3*U*Y*S_num)
assert all(value(0, 0, 0) for value in (A_num, B_num, H, S_den, H_reverse))


def common_cleared_numerator(entry):
    """Return the exact local numerator from the certified principal frame."""
    a, b, i, k = tuple(int(entry[key]) for key in (
        "x_power", "m_power", "u_power", "h_power"
    ))
    t_exponent = T+4*k-i-9
    assert t_exponent >= 0
    factors = (
        Z**4*Y**6,
        t_value**t_exponent,
        x_value**a,
        (Z**2*U*Y)**b,
        truncated_power(B_num, b),
        truncated_power(R_den, b),
        truncated_power(A_num, 9-b),
        truncated_power(H, 9-b),
        truncated_power(S_den, 9-b),
        truncated_power(H_reverse, K-k),
    )
    answer = ring.one()
    for factor in factors:
        answer = truncated_product(answer, factor)
    return answer


# ``ds`` is Singular's local degree order.  Its standard normal form is in
# the local ring at (Z,U,Y), precisely the ring in which the chart-unit
# comparisons made by the clearing certificate are valid.
singular.eval("ring q8node={},(Z,U,Y),ds;".format(args.prime))
singular.eval("poly surface={};".format(surface))
singular.eval("ideal principal=surface,(Z3*U2)^{};".format(T))
singular.eval("ideal standard_principal=std(principal);")
assert int(singular.eval("size(standard_principal);")) >= 2

remainders = []
for entry in ambient["ambient_basis"]:
    numerator = common_cleared_numerator(entry)
    singular.eval("poly node_num={};".format(numerator))
    local_remainder = ring(singular("reduce(node_num,standard_principal)").sage())
    # A normal form is an exact representative of this ambient element in the
    # infinite local quotient.  The collection of such representatives has a
    # finite span because the ambient itself is finite.
    singular.eval("poly local_remainder={};".format(local_remainder))
    re_reduced = ring(singular("reduce(local_remainder,standard_principal)").sage())
    assert re_reduced == local_remainder
    remainders.append(local_remainder)

image_condition = finite_ambient_image_condition(
    "actual E7_4--E7_3 principal local-normal-form image",
    tuple(range(len(remainders))),
    lambda index: {tuple(monomial): coefficient for monomial, coefficient in remainders[index].dict().items()},
    lambda monomial: monomial,
    finite,
    "actual H92 resolved chart; local standard normal form modulo (surface,t^T)",
)
residue_matrix = image_condition["matrix"]
monomials = image_condition["coordinate_keys"]
monomial_index = {monomial: index for index, monomial in enumerate(monomials)}

payload = {
    "schema": "elkies-k3.h92-q8-e7-4-3-principal-node-local-normal-form-modp.v1",
    "status": "EXPERIMENTAL_MODULAR_Q8_E7_4_3_LOCAL_NORMAL_FORM_BLOCK",
    "prime": args.prime,
    "inputs": {
        "checker_source": {"path": str(Path(__file__).relative_to(ROOT)), "sha256": digest(Path(__file__))},
        "p1": {"path": str(args.p1.relative_to(ROOT)), "sha256": digest(args.p1)},
        "actual_pullbacks": {"path": str(args.pullbacks.relative_to(ROOT)), "sha256": digest(args.pullbacks)},
        "endpoint_ambient": {"path": str(args.ambient.relative_to(ROOT)), "sha256": digest(args.ambient)},
        "principal_clearing": {"path": str(args.clearing.relative_to(ROOT)), "sha256": digest(args.clearing)},
        "compiler_core": {"path": str(CORE.relative_to(ROOT)), "sha256": digest(CORE)},
    },
    "local_ring": {
        "chart": "E7_4--E7_3",
        "order": "Singular ds local degree order at (Z,U,Y)",
        "surface_equation": str(surface),
        "principal_ideal": "(surface,(Z^3*U^2)^%d)" % T,
        "completed_interpretation": "t^%d=Z^%d*Y^%d*unit" % (T, 3*T, 4*T),
    },
    "common_clearing": {
        "T": T,
        "K": K,
        "unit_multiplier": "A_num^9*H^9*s_den^9*h_reverse^{}".format(K),
        "formula": "the certified common-cleared numerator is reduced modulo (surface,t^T)",
        "exact_truncation": "during multiplication discard only monomials divisible by (Z^3*U^2)^T",
    },
    "good_reduction": {
        "all_input_coefficient_denominators_nonzero": True,
        "node_unit_residues": {
            "A_num": int(A_num(0, 0, 0)),
            "B_num": int(B_num(0, 0, 0)),
            "H": int(H(0, 0, 0)),
            "s_den": int(S_den(0, 0, 0)),
            "h_reverse": int(H_reverse(0, 0, 0)),
        },
        "argument": (
            "The common-clearing factors remain units at the chart origin after "
            "reduction. Thus a characteristic-zero local membership relation can "
            "be reduced at this prime after primitive normalization."
        ),
    },
    "finite_ambient_image": {
        "ambient_dimension": len(remainders),
        "coordinate_monomials": ["Z^{}*U^{}*Y^{}".format(*monomial) for monomial in monomials],
        "rows": len(monomials),
        "rank": int(residue_matrix.rank()),
        "kernel_dimension": int(residue_matrix.right_kernel().dimension()),
        "sparse_columns": [
            [
                [monomial_index[tuple(monomial)], int(coefficient)]
                for monomial, coefficient in sorted(remainder.dict().items(), key=lambda item: tuple(item[0]))
            ]
            for remainder in remainders
        ],
    },
    "boundary": (
        "The coordinate space is the finite image of this displayed ambient in "
        "the local quotient, not a finite presentation of R/(t^T). This is a "
        "single-prime resolved-chart computation and does not provide a "
        "characteristic-zero q8 matrix, overlap compatibility, h0, a pencil, "
        "or a child equation."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q8E743LOCALNF|prime={}|ambient={}|rows={}|rank={}|kernel={}|"
    "status=EXPERIMENTAL_MODULAR_Q8_E7_4_3_LOCAL_NORMAL_FORM_BLOCK".format(
        args.prime, len(remainders), len(monomials), residue_matrix.rank(),
        residue_matrix.right_kernel().dimension(),
    ),
    flush=True,
)
