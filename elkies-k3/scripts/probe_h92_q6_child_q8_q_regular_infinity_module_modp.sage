#!/usr/bin/env sage -python
"""Probe the q-regular finite module against the smooth infinity frame.

The exact finite q-regular module has polynomial generators

    (B,C)=(f_IV, C_IV), (0, M),  M=f_II^2*f_IV^3,

for a function ``C+B*q_regular``.  At the smooth fibre at infinity use the
minimal coordinates ``x=T^4 X, y=T^6 Y``.  If ``m_inf`` is the chord in
``(X,Y)``, then

    q_regular = alpha + beta*m_inf,
    alpha=-p/h-R/Nx, beta=T^2/h.

This gives a completely explicit *candidate* infinity filter on polynomial
combinations ``s*(f_IV,C_IV)+t*(0,M)``: both displayed coefficients of
``1,m_inf`` must have prescribed nonnegative u-orders.  It is a modular
linear-algebra probe of that natural frame; it does not yet certify that the
chosen orders are the complete q8 infinity sheaf.
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, matrix


ROOT = Path(__file__).resolve().parents[2]
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
MARKING = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-marking.json"
NORMALIZER = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-q-pole-normalization-crt.json"
FINITE = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-q-regular-finite-module-qq.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-q-regular-infinity-module-mod-43.json"


def coefficient(field, value):
    value = QQ(value)
    denominator = field(value.denominator())
    if not denominator:
        raise ValueError("prime divides an input denominator")
    return field(value.numerator()) / denominator


def polynomial(ring, field, coefficients):
    return ring([coefficient(field, value) for value in coefficients])


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
parser.add_argument("--finite", type=Path, default=FINITE)
parser.add_argument("--max-s-degree", type=int, default=42)
parser.add_argument("--max-t-degree", type=int, default=42)
parser.add_argument(
    "--a-infinity-order", type=int, default=1,
    help="required u-order of C+B*alpha in the minimal smooth frame",
)
parser.add_argument(
    "--b-infinity-order", type=int, default=1,
    help="required u-order of B*beta in the minimal smooth frame",
)
parser.add_argument(
    "--screen-branch-levels", action="store_true",
    help="when the candidate kernel has dimension two, screen every finite ratio level",
)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
if not args.prime > 3 or not QQ(args.prime).is_integer() or not args.max_s_degree >= 0 or not args.max_t_degree >= 0:
    raise ValueError("prime must exceed three and degree bounds must be nonnegative")
if not args.a_infinity_order >= 0 or not args.b_infinity_order >= 0:
    raise ValueError("infinity orders must be nonnegative")

child = json.loads(CHILD.read_text())
marking = json.loads(MARKING.read_text())
normalizer = json.loads(NORMALIZER.read_text())
finite_path = args.finite.resolve()
finite_module = json.loads(finite_path.read_text())
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert marking["status"] == "PASS_EXACT_Q6_CHILD_Q8_MARKING"
assert normalizer["status"] == "PASS_EXACT_CRT_PRINCIPAL_PART_NORMALIZATION"
assert finite_module["status"] == "PASS_EXACT_Q_REGULAR_FINITE_MODULE"

finite = GF(args.prime)
ring = PolynomialRing(finite, "T")
T = ring.gen()
field = ring.fraction_field()
section = marking["selected_q8"]["relative_child_section_standard_jacobian_coordinates"]
nx = polynomial(ring, finite, section["x_numerator_coefficients_low_to_high"])
dx = polynomial(ring, finite, section["x_denominator_coefficients_low_to_high"])
ny = polynomial(ring, finite, section["y_numerator_coefficients_low_to_high"])
dy = polynomial(ring, finite, section["y_denominator_coefficients_low_to_high"])
sx, sy = field(nx) / field(dx), field(ny) / field(dy)
h = monic_power_root(dx, 2)
assert dx // h**2 in finite and dy // h**3 in finite
normalizer_R = polynomial(ring, finite, normalizer["normalizer"]["R_coefficients_low_to_high"])
assert (normalizer_R*h*dy-ny) % nx == 0

ii = polynomial(ring, finite, PolynomialRing(QQ, "T")(next(
    item for item in child["finite_fibres"] if item["kodaira"] == "II*"
)["factor"]).list())
iv = polynomial(ring, finite, PolynomialRing(QQ, "T")(next(
    item for item in child["finite_fibres"] if item["kodaira"] == "IV*"
)["factor"]).list())
nef_ivstar = finite_module["module"]["basis"][0][0] == "1"
M = ii**2 * iv**(2 if nef_ivstar else 3)
C_iv = polynomial(ring, finite, finite_module["module"]["first_C_lift_coefficients_low_to_high"])
first_b = ring.one() if nef_ivstar else iv
assert M.degree() == (4 if nef_ivstar else 5) and iv.degree() == 1 and C_iv.degree() < M.degree()

# q_regular=(m-p)/h-R/Nx and m=T^2*m_inf in the minimal smooth model.
p = -sy / sx
alpha = -p / field(h) - field(normalizer_R) / field(nx)
beta = field(T**2) / field(h)
alpha_num, alpha_den = ring(alpha.numerator()), ring(alpha.denominator())
beta_num, beta_den = ring(beta.numerator()), ring(beta.denominator())
assert alpha_den(0) and beta_den(0)

labels = [("s", exponent) for exponent in range(args.max_s_degree + 1)] + [
    ("t", exponent) for exponent in range(args.max_t_degree + 1)
]

def column(label):
    kind, exponent = label
    if kind == "s":
        B, C = first_b*T**exponent, C_iv*T**exponent
    else:
        B, C = ring.zero(), M*T**exponent
    # C+B*alpha=(C*alpha_den+B*alpha_num)/alpha_den.
    a_numerator = C*alpha_den + B*alpha_num
    b_numerator = B*beta_num
    return a_numerator, b_numerator

columns = [column(label) for label in labels]

def high_coefficient_rows(values, denominator_degree, required_order):
    # A rational function N/D has u-order >= e iff deg(N)<=deg(D)-e.
    cutoff = denominator_degree - required_order
    top = max((degree_or_minus_one(value) for value in values), default=-1)
    return [
        [value[degree] if degree <= value.degree() else finite.zero() for value in values]
        for degree in range(cutoff + 1, top + 1)
    ], cutoff, top

a_rows, a_cutoff, a_top = high_coefficient_rows(
    [entry[0] for entry in columns], alpha_den.degree(), args.a_infinity_order
)
b_rows, b_cutoff, b_top = high_coefficient_rows(
    [entry[1] for entry in columns], beta_den.degree(), args.b_infinity_order
)
condition = matrix(finite, a_rows + b_rows, ncols=len(labels))
kernel = condition.right_kernel_matrix()


def pair_from_row(row):
    """Reconstruct C+B*q_regular from a module-coordinate row."""

    s = sum(row[index] * T**exponent for index, (kind, exponent) in enumerate(labels) if kind == "s")
    t = sum(row[index] * T**exponent for index, (kind, exponent) in enumerate(labels) if kind == "t")
    B = iv*s
    C = C_iv*s + M*t
    a = field(C) - field(B)*p/field(h) - field(B*normalizer_R)/field(nx)
    b = field(B)/field(h)
    return a, b


branch = None
if args.screen_branch_levels:
    if kernel.nrows() != 2:
        raise ValueError("branch-level screen requires a two-dimensional candidate kernel")
    A = polynomial(ring, finite, child["minimal_short_weierstrass"]["A_coefficients_low_to_high"])
    B_curve = polynomial(ring, finite, child["minimal_short_weierstrass"]["B_coefficients_low_to_high"])
    a0, b0 = pair_from_row(kernel.row(0))
    a1, b1 = pair_from_row(kernel.row(1))
    x_ring = PolynomialRing(field, "x")
    x = x_ring.gen()
    hist, good, singular_levels = {}, [], []
    for level in finite:
        denominator = b1-field(level)*b0
        if not denominator:
            singular_levels.append(int(level))
            continue
        m_value = -(a1-field(level)*a0)/denominator
        y = x_ring(m_value)*(x-x_ring(sx))-x_ring(sy)
        relation = y**2-x**3-x_ring(A)*x-x_ring(B_curve)
        quadratic, remainder = relation.quo_rem(x-x_ring(sx))
        assert not remainder and quadratic.degree() == 2
        discriminant = x_ring.base_ring()(quadratic[1]**2-4*quadratic[2]*quadratic[0])
        numerator, denominator = ring(discriminant.numerator()), ring(discriminant.denominator())
        odd_degree = sum(
            factor.degree()
            for value in (numerator, denominator)
            for factor, multiplicity in value.squarefree_decomposition()
            if multiplicity % 2
        )
        branch_degree = odd_degree + (denominator.degree()-numerator.degree()) % 2
        hist[branch_degree] = hist.get(branch_degree, 0) + 1
        if branch_degree == 4:
            good.append(int(level))
    branch = {
        "finite_levels": int(args.prime),
        "singular_levels": singular_levels,
        "branch_degree_histogram": {str(key): value for key, value in sorted(hist.items())},
        "genus_one_levels": good,
    }

payload = {
    "schema": "elkies-k3.h92-q6-child-q8-q-regular-infinity-module-modp.v1",
    "status": "EXPERIMENTAL_Q_REGULAR_INFINITY_MODULE_PROBE",
    "prime": int(args.prime),
    "candidate_frame": {
        "minimal_scaling": "x=T^4*X, y=T^6*Y, m=T^2*m_inf",
        "q_regular": "alpha+beta*m_inf",
        "alpha": "-p/h-R/Nx",
        "beta": "T^2/h",
        "alpha_degrees": [int(alpha_num.degree()), int(alpha_den.degree())],
        "beta_degrees": [int(beta_num.degree()), int(beta_den.degree())],
    },
    "finite_module": {
        "generators": (["(1,C_IV)", "(0,M=f_II^2*f_IV^2)"] if nef_ivstar else ["(f_IV,C_IV)", "(0,M=f_II^2*f_IV^3)"]),
        "C_IV_degree": int(C_iv.degree()),
        "M_degree": int(M.degree()),
    },
    "ambient": {
        "s_degree": int(args.max_s_degree),
        "t_degree": int(args.max_t_degree),
        "dimension": len(labels),
        "labels": [[kind, exponent] for kind, exponent in labels],
    },
    "infinity_filter": {
        "a_required_u_order": int(args.a_infinity_order),
        "b_required_u_order": int(args.b_infinity_order),
        "a_degree_cutoff": int(a_cutoff),
        "b_degree_cutoff": int(b_cutoff),
        "a_rows": len(a_rows),
        "b_rows": len(b_rows),
    },
    "result": {
        "rank": int(condition.rank()),
        "kernel_dimension": int(kernel.nrows()),
        "kernel_basis": [[int(value) for value in row] for row in kernel.rows()],
    },
    "branch_screen": branch,
    "boundary": (
        "This is a modular probe of a natural smooth-infinity frame for the exact "
        "finite q-regular module. The asserted infinity orders have not yet been "
        "identified with the complete q8 divisor sheaf, so it is not a global pencil, "
        "D13/rootless equation, bisection cover, collision, or rank certificate."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6CHILDQREGINFINITY|prime={}|ambient={}|a_rows={}|b_rows={}|rank={}|"
    "kernel={}|branch_histogram={}|genus_one_levels={}|status=EXPERIMENTAL_Q_REGULAR_INFINITY_MODULE_PROBE".format(
        args.prime, len(labels), len(a_rows), len(b_rows), condition.rank(), kernel.nrows(),
        "none" if branch is None else ",".join(
            "{}:{}".format(key, value) for key, value in branch["branch_degree_histogram"].items()
        ),
        "none" if branch is None or not branch["genus_one_levels"] else ",".join(map(str, branch["genus_one_levels"])),
    ),
    flush=True,
)
