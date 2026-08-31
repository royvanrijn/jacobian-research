#!/usr/bin/env sage -python
"""Identify the three remaining q12/o5867 compiler branches modulo p.

Invert each filtered polynomial P.O=0 section on the P1229-pointed q8 child
to a point (T(u),W(u)) on the resolved quartic and then to the corresponding
point (x(u),y(u)) on the old q4/o164 Weierstrass surface by the exact chord
formula.  At several ordinary old-base values, split T(u)-t0 over GF(p^6),
sum the conjugate q4 points, and compare literally with the exact marked Abel
word.  No Groebner basis, multivariate elimination, or QQ lift is used.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
Q8 = LOCAL / "q4o164-q8o376-smooth-rr-qq.json"
MODEL = LOCAL / "q4o164-compact-weierstrass-qq.json"
BASIS = LOCAL / "q4o164-integral-basis-qq.json"
C8 = LOCAL / "q4o164-c8-equation-marking-qq.json"
HORIZONTAL = LOCAL / "q4o164-q8o376-horizontal-crt-qq.json"
Q3 = LOCAL / "q12o5867-degree1-compiler-branch-qq.json"
INPUTS_BASE = (Q8, MODEL, BASIS, C8, HORIZONTAL, Q3)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=83)
parser.add_argument("--shell", type=Path)
parser.add_argument("--output", type=Path)
parser.add_argument("--fibres", type=int, default=4)
args = parser.parse_args()
prime = ZZ(args.prime)
if not prime.is_prime() or prime in (2, 3):
    raise ValueError("prime must be an odd prime other than 3")
SHELL = (
    args.shell.resolve() if args.shell else
    LOCAL / f"q12o5867-p0-shell-word-fingerprints-mod{prime}.json"
)
OUTPUT = (
    args.output.resolve() if args.output else
    LOCAL / f"q12o5867-abel-trace-named-seeds-mod{prime}.json"
)
started = time.monotonic()
F = GF(prime)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


q8 = json.loads(Q8.read_text())
model = json.loads(MODEL.read_text())
basis = json.loads(BASIS.read_text())
c8 = json.loads(C8.read_text())
horizontal = json.loads(HORIZONTAL.read_text())
q3_exact = json.loads(Q3.read_text())
shell = json.loads(SHELL.read_text())
assert q8["status"] == "PASS_EXACT_QQ_Q8O376_4A1_RR_JACOBIAN_AND_P1229_ZERO"
assert shell["prime"] == int(prime)
assert shell["status"].startswith("PASS_MODP_Q12O5867_P0_SHELL")

PT = PolynomialRing(F, "T")
T = PT.gen()
PK = PT.fraction_field()
U = PolynomialRing(F, "u")
u = U.gen()
KU = U.fraction_field()
TU = PolynomialRing(KU, "T")
Tvar = TU.gen()
TUP = PolynomialRing(U, "T")


def reduce_qq(value):
    value = QQ(value)
    assert value.denominator() % prime
    return F(value.numerator())/F(value.denominator())


def parent_poly(values):
    return PT([reduce_qq(value) for value in values])


def child_poly(values):
    return U([reduce_qq(value) for value in values])


def parent_function(record):
    return PK(parent_poly(record["numerator_coefficients_low_to_high"])) / PK(
        parent_poly(record["denominator_coefficients_low_to_high"])
    )


def child_function(record):
    return KU(child_poly(record["numerator_coefficients_low_to_high"])) / KU(
        child_poly(record["denominator_coefficients_low_to_high"])
    )


parent_A = parent_poly(model["compact_model"]["A_coefficients_low_to_high"])
parent_B = parent_poly(model["compact_model"]["B_coefficients_low_to_high"])
child_A = child_poly(q8["child"]["minimal_A_coefficients_low_to_high"])
child_B = child_poly(q8["child"]["minimal_B_coefficients_low_to_high"])


def checked_parent_point(x_coordinate, y_coordinate):
    point = PK(x_coordinate), PK(y_coordinate)
    assert point[1]**2 == point[0]**3+PK(parent_A)*point[0]+PK(parent_B)
    return point


def parent_neg(point):
    return None if point is None else (point[0], -point[1])


def parent_add(left, right):
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2:
        if y1 == -y2:
            return None
        slope = (3*x1**2+PK(parent_A))/(2*y1)
    else:
        slope = (y2-y1)/(x2-x1)
    x3 = slope**2-x1-x2
    return checked_parent_point(x3, slope*(x1-x3)-y1)


def parent_mul(coefficient, point):
    coefficient = ZZ(coefficient)
    if coefficient < 0:
        return parent_mul(-coefficient, parent_neg(point))
    answer = None
    addend = point
    while coefficient:
        if coefficient & 1:
            answer = parent_add(answer, addend)
        addend = parent_add(addend, addend)
        coefficient >>= 1
    return answer


parent_basis = [
    checked_parent_point(
        parent_poly(section["x_coefficients_low_to_high"]),
        parent_poly(section["y_coefficients_low_to_high"]),
    )
    for section in basis["resolved_hensel"]["sections"]
]
H = checked_parent_point(
    parent_function(horizontal["section"]["x"]),
    parent_function(horizontal["section"]["y"]),
)
c8_record = c8["opposite_constant_support_section"]
c8_old = parent_function(c8_record["x"]), parent_function(c8_record["y"])
base_scale = reduce_qq(model["exact_coordinate_change"]["c"])
xy_scale = reduce_qq(model["exact_coordinate_change"]["s"])


def substitute_scaled(function):
    return PK(PT(function.numerator()(base_scale*T))) / PK(
        PT(function.denominator()(base_scale*T))
    )


C8opposite = checked_parent_point(
    substitute_scaled(c8_old[0])/xy_scale**2,
    substitute_scaled(c8_old[1])/xy_scale**3,
)


def word_point(terms):
    answer = None
    for coefficient, point in terms:
        answer = parent_add(answer, parent_mul(coefficient, point))
    return answer


target_words = {
    "Q1": word_point(((1, H), (-1, parent_basis[3]))),
    "Q2": word_point((
        (1, H), (1, parent_basis[0]), (1, parent_basis[2]),
        (-1, parent_basis[3]), (1, parent_basis[4]),
    )),
    "Q4": word_point((
        (1, H), (-1, C8opposite), (-1, parent_basis[1]),
        (-1, parent_basis[2]), (1, parent_basis[3]),
        (1, parent_basis[5]), (-1, parent_basis[6]), (-1, parent_basis[7]),
    )),
    "Q3": word_point((
        (2, C8opposite), (1, parent_basis[1]), (2, parent_basis[2]),
        (-2, parent_basis[3]), (1, parent_basis[4]), (-2, parent_basis[5]),
        (1, parent_basis[6]), (1, parent_basis[7]),
    )),
}

# Reconstruct the normalized q8 quartic square factor.
XH = parent_poly(horizontal["section"]["x"]["numerator_coefficients_low_to_high"])
denominator_x = parent_poly(horizontal["section"]["x"]["denominator_coefficients_low_to_high"])
YH = parent_poly(horizontal["section"]["y"]["numerator_coefficients_low_to_high"])
denominator_y = parent_poly(horizontal["section"]["y"]["denominator_coefficients_low_to_high"])
Z = PT.one()
for factor, exponent in denominator_x.factor():
    assert int(exponent) % 2 == 0
    Z *= factor.monic()**(int(exponent)//2)
assert Z**2 == denominator_x and Z**3 == denominator_y
resolved_pairs = [
    (
        parent_poly(item["AA_coefficients_low_to_high"]),
        parent_poly(item["BB_coefficients_low_to_high"]),
    )
    for item in q8["resolved_RR"]["resolved_basis_pairs"]
]
(AA0, BB0), (AA1, BB1) = resolved_pairs


def lift_TU(poly):
    return TU([KU(value) for value in PT(poly).list()])


aa = lift_TU(AA0)+KU(u)*lift_TU(AA1)
bb = lift_TU(BB0)+KU(u)*lift_TU(BB1)
X_u, Y_u, Z_u, A_u = map(lift_TU, (XH, YH, Z, parent_A))
raw = (
    aa**4-6*X_u*aa**2*bb**2+8*Y_u*aa*bb**3
    -3*X_u**2*bb**4-4*A_u*bb**4*Z_u**4
)
after_collision, remainder = raw.quo_rem(Z_u**4)
assert not remainder
quartic = TUP([
    U([reduce_qq(value) for value in coefficient])
    for coefficient in q8["quartic"]["coefficients_in_old_T_low_to_high"]
])
square_quotient, remainder = TUP(after_collision).quo_rem(quartic)
assert not remainder
factorization = square_quotient.factor()
assert all(int(exponent) % 2 == 0 for unused, exponent in factorization)
unit = F(factorization.unit())
assert unit.is_square()
square_factor = TUP(unit.sqrt())
for factor, exponent in factorization:
    square_factor *= factor**(int(exponent)//2)
assert TUP(after_collision) == quartic*square_factor**2


def evaluate_bivariate(poly, old_t):
    answer = KU.zero()
    for coefficient in reversed(TUP(poly).list()):
        answer = answer*old_t+KU(U(coefficient))
    return answer


def evaluate_parent_poly(poly, old_t):
    answer = KU.zero()
    for coefficient in reversed(PT(poly).list()):
        answer = answer*old_t+KU(coefficient)
    return answer


def compose_parent(function, old_t):
    function = PK(function)
    return evaluate_parent_poly(function.numerator(), old_t) / evaluate_parent_poly(
        function.denominator(), old_t
    )


# Pointed-quartic inverse on the child.
pointing = q8["preferred_pointed_zero"]
q_origin = child_function(pointing["quartic_ordinate"])
p1229_support = reduce_qq(pointing["old_base_coordinate"])
quartic_coefficients = [KU(U(coefficient)) for coefficient in TUP(quartic).list()]
translated = []
from math import comb
for new_degree in range(5):
    translated.append(sum(
        quartic_coefficients[old_degree]*F(comb(old_degree, new_degree))
        * KU(p1229_support)**(old_degree-new_degree)
        for old_degree in range(new_degree, 5)
    ))
e, d, c, b, a = translated
assert e == q_origin**2
a1 = d/q_origin
a2 = c-d**2/(4*q_origin**2)
a3 = 2*q_origin*b
b2 = a1**2+4*a2


def invert_child_section(X, Y):
    X = KU(X)
    Y = KU(Y)
    x_general = X/9-b2/12
    y_general = Y/27-(a1*x_general+a3)/2
    local_t = 2*q_origin*(x_general+a2)/y_general
    old_t = KU(p1229_support)+local_t
    ordinate = (x_general*local_t**2-d*local_t)/(2*q_origin)-q_origin
    assert ordinate**2 == evaluate_bivariate(quartic, old_t)
    AA = evaluate_parent_poly(AA0, old_t)+KU(u)*evaluate_parent_poly(AA1, old_t)
    BB = evaluate_parent_poly(BB0, old_t)+KU(u)*evaluate_parent_poly(BB1, old_t)
    ZZ = evaluate_parent_poly(Z, old_t)
    slope = -AA/(BB*ZZ)
    square = evaluate_bivariate(square_factor, old_t)
    hx = compose_parent(H[0], old_t)
    hy = compose_parent(H[1], old_t)
    old_x = (ordinate*square/BB**2-hx+slope**2)/2
    old_y = slope*(old_x-hx)-hy
    assert old_y**2 == compose_parent(parent_A, old_t)*old_x + old_x**3 + compose_parent(parent_B, old_t)
    restrictions = [
        evaluate_parent_poly(AAi, old_t)+evaluate_parent_poly(BBi, old_t)*ZZ*slope
        for AAi, BBi in resolved_pairs
    ]
    assert KU(u) == -restrictions[0]/restrictions[1]
    return old_t, ordinate, old_x, old_y


def rational_degree(value):
    return max(value.numerator().degree(), value.denominator().degree())


def padded(poly, length):
    return [poly[i] if i <= poly.degree() else F.zero() for i in range(length)]


def jacobian_rank(X, Y):
    dx = -3*X**2-child_A
    dy = 2*Y
    columns = [padded(u**power*dx, 13) for power in range(5)]
    columns += [padded(u**power*dy, 13) for power in range(7)]
    return int(matrix(F, columns).transpose().rank())


# Prove the inverse formulas on the exact degree-one Q3 branch.
q3_section = q3_exact["section"]
Q3_child = (
    child_poly(q3_section["x_coefficients_low_to_high"]),
    child_poly(q3_section["y_coefficients_low_to_high"]),
)
q3_old_t, unused_q3_w, q3_old_x, q3_old_y = invert_child_section(*Q3_child)
assert rational_degree(q3_old_t) == 1
assert q3_old_x == compose_parent(target_words["Q3"][0], q3_old_t)
assert q3_old_y == compose_parent(target_words["Q3"][1], q3_old_t)

# A common degree-six extension splits every separable quadratic and cubic.
E = GF(prime**6, "z")
EU = PolynomialRing(E, "v")


def evaluate_over_extension(function, value):
    function = KU(function)

    def evaluate(poly):
        answer = E.zero()
        for coefficient in reversed(U(poly).list()):
            answer = answer*value+E(coefficient)
        return answer

    denominator = evaluate(function.denominator())
    if denominator == 0:
        raise ZeroDivisionError
    return evaluate(function.numerator())/denominator


def evaluate_parent_at(function, value):
    function = PK(function)
    numerator = E(PT(function.numerator())(value))
    denominator = E(PT(function.denominator())(value))
    if denominator == 0:
        raise ZeroDivisionError
    return numerator/denominator


def fibre_add(left, right, coefficient_A):
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2:
        if y1 == -y2:
            return None
        if not y1:
            return None
        slope = (3*x1**2+coefficient_A)/(2*y1)
    else:
        slope = (y2-y1)/(x2-x1)
    x3 = slope**2-x1-x2
    return x3, slope*(x1-x3)-y1


def trace_trials(old_t, old_x, old_y, target, expected_degree, needed):
    trials = []
    for integer_t0 in range(int(prime)):
        if len(trials) >= needed:
            break
        t0 = F(integer_t0)
        if (4*parent_A(t0)**3+27*parent_B(t0)**2) == 0:
            continue
        polynomial = U(old_t.numerator()-t0*old_t.denominator())
        if polynomial.degree() != expected_degree or not polynomial.is_squarefree():
            continue
        roots = EU([E(value) for value in polynomial.list()]).roots(multiplicities=False)
        if len(roots) != expected_degree:
            continue
        try:
            points = [
                (evaluate_over_extension(old_x, root), evaluate_over_extension(old_y, root))
                for root in roots
            ]
            target_value = (
                evaluate_parent_at(target[0], t0), evaluate_parent_at(target[1], t0)
            )
        except ZeroDivisionError:
            continue
        coefficient_A = E(parent_A(t0))
        coefficient_B = E(parent_B(t0))
        if any(y_value**2 != x_value**3+coefficient_A*x_value+coefficient_B for x_value, y_value in points):
            raise ArithmeticError("recovered conjugate missed the old q4 fibre")
        trace = None
        for point in points:
            trace = fibre_add(trace, point, coefficient_A)
        trials.append({
            "old_base_value": int(t0),
            "root_count": len(roots),
            "literal_target_identity": trace == target_value,
            "literal_negative_target_identity": (
                trace is not None
                and trace == (target_value[0], -target_value[1])
            ),
        })
    return trials


expected_degrees = {"Q1": 3, "Q2": 2, "Q4": 2}
branch_results = {}
for name in ("Q1", "Q2", "Q4"):
    candidates = shell["branches"][name]["candidates"]
    rows = []
    for candidate_index, candidate in enumerate(candidates):
        X = U(candidate["x_coefficients_low_to_high"])
        Y = U(candidate["y_coefficients_low_to_high"])
        old_t, unused_W, old_x, old_y = invert_child_section(X, Y)
        degree = rational_degree(old_t)
        assert degree == expected_degrees[name]
        trials = trace_trials(
            old_t, old_x, old_y, target_words[name], degree, args.fibres
        )
        passed = len(trials) == args.fibres and all(
            trial["literal_target_identity"] for trial in trials
        )
        rows.append({
            "candidate_index": candidate_index,
            "shell_index": candidate["shell_index"],
            "x_coefficients_low_to_high": candidate["x_coefficients_low_to_high"],
            "y_coefficients_low_to_high": candidate["y_coefficients_low_to_high"],
            "ordinary_coefficient_jacobian_rank": jacobian_rank(X, Y),
            "inverse_parent_degree": int(degree),
            "ordinary_fibre_trials": trials,
            "passes_all_Abel_trace_identities": passed,
        })
    survivors = [row for row in rows if row["passes_all_Abel_trace_identities"]]
    branch_results[name] = {
        "candidate_count": len(rows),
        "named_seed_count": len(survivors),
        "missing_named_match": len(survivors) == 0,
        "residual_ambiguity": len(survivors) > 1,
        "named_seeds": survivors,
        "all_candidates": rows,
    }

unique = all(branch_results[name]["named_seed_count"] == 1 for name in branch_results)
rank12 = unique and all(
    branch_results[name]["named_seeds"][0]["ordinary_coefficient_jacobian_rank"] == 12
    for name in branch_results
)
missing = any(branch_results[name]["named_seed_count"] == 0 for name in branch_results)
def display_path(path):
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


payload = {
    "schema": "elkies-k3.h92-q12o5867-abel-trace-named-seeds-modp.v1",
    "status": (
        "PASS_MODP_Q12O5867_THREE_NAMED_REGULAR_SEEDS_READY_FOR_QQ_HENSEL"
        if rank12 else (
            "PARTIAL_MODP_Q12O5867_ABEL_TRACE_MISSING_NAMED_BRANCH"
            if missing else
            "PASS_MODP_Q12O5867_ABEL_TRACE_WITH_RESIDUAL_AMBIGUITY"
        )
    ),
    "prime": int(prime),
    "exact_Q3_inverse_regression": {
        "pass": True,
        "inverse_parent_degree": 1,
        "literal_recovered_parent_word_identity": True,
    },
    "extension_degree_for_splitting": 6,
    "required_ordinary_fibres_per_candidate": args.fibres,
    "branches": branch_results,
    "method": {
        "large_Groebner_required": False,
        "elimination_required": False,
        "QQ_lift_attempted": False,
        "construction": "pointed-quartic inverse, exact chord inverse, finite-extension conjugate sum",
        "runtime_seconds": time.monotonic()-started,
    },
    "proof_boundary": (
        "Every displayed identity is exact over GF(p) or GF(p^6) for the displayed prime. The finite-fibre "
        "Abel comparisons name modular shell points; they do not constitute a QQ lift, "
        "a q12 resolved Riemann--Roch construction, or an endpoint equation certificate."
    ),
    "inputs": {
        "paths": [display_path(path) for path in INPUTS_BASE+(SHELL,)],
        "sha256": {display_path(path): sha256(path) for path in INPUTS_BASE+(SHELL,)},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q12O5867ABEL|prime={}|Q1={}/{}|Q2={}/{}|Q4={}/{}|rank12={}|runtime={:.3f}|"
    "status={}|output={}".format(
        prime, branch_results["Q1"]["named_seed_count"], branch_results["Q1"]["candidate_count"],
        branch_results["Q2"]["named_seed_count"], branch_results["Q2"]["candidate_count"],
        branch_results["Q4"]["named_seed_count"], branch_results["Q4"]["candidate_count"],
        int(rank12), payload["method"]["runtime_seconds"], payload["status"], OUTPUT,
    ), flush=True,
)
