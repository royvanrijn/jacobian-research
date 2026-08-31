#!/usr/bin/env sage-python
"""Filter q323 target-shape sections by their smooth-RR child over GF(61).

For every stored P.O=10 candidate of shape (24,20,34,30), compute the direct
D=O+P-4F congruence module (ambient dimension 22, modulus degree 20).  When
h0=2, compile its binary quartic and minimal Jacobian and retain precisely
the candidates whose child root rank is five.  This uses finite-field linear
algebra and bivariate factorization only; there is no elimination or
Groebner basis.
"""

import argparse
import hashlib
import json
import time
from collections import Counter
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
CANDIDATES = LOCAL / "q4o323-q207-four-section-words-mod61.json"
POINTING = LOCAL / "q4o323-component2-pointing-qq.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--start", type=int, default=0)
parser.add_argument("--stop", type=int)
parser.add_argument(
    "--output", type=Path,
    default=LOCAL / "q4o323-q207-smooth-rr-filter-mod61.json",
)
args = parser.parse_args()
started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


prime = ZZ(61)
F = GF(prime)
R = PolynomialRing(F, "u")
u = R.gen()
candidate_data = json.loads(CANDIDATES.read_text())
pointing = json.loads(POINTING.read_text())
assert candidate_data["status"] == "PASS_MOD61_Q4O323_Q207_FOUR_SECTION_CANDIDATES"
assert pointing["status"] == "PASS_EXACT_QQ_Q4O323_OLD_A11_COMPONENT2_POINTING"


def reduce_qq(value):
    value = QQ(value)
    return F(value.numerator()) / F(value.denominator())


A = R([reduce_qq(value) for value in pointing["global_short_model"]["A_coefficients_low_to_high"]])
B = R([reduce_qq(value) for value in pointing["global_short_model"]["B_coefficients_low_to_high"]])


def section_polynomials(record):
    X = R(record["x"]["numerator_coefficients_low_to_high"])
    Y = R(record["y"]["numerator_coefficients_low_to_high"])
    denominator_x = R(record["x"]["denominator_coefficients_low_to_high"])
    denominator_y = R(record["y"]["denominator_coefficients_low_to_high"])
    if not denominator_x.is_square():
        raise ArithmeticError("x denominator is not a square")
    Z = denominator_x.sqrt()
    if Z**3 != denominator_y:
        Z = -Z
    if Z**2 != denominator_x or Z**3 != denominator_y:
        raise ArithmeticError("incompatible projective denominators")
    if Y**2 != X**3+A*X*Z**4+B*Z**6:
        raise ArithmeticError("candidate misses the child Weierstrass equation")
    return X, Y, Z


def smooth_basis(X, Y, Z):
    ambient = (
        [(u**degree, R.zero()) for degree in range(17)]
        + [(R.zero(), u**degree) for degree in range(5)]
    )
    modulus = Z**2
    remainders = [(AA*X-BB*Y) % modulus for AA, BB in ambient]
    conditions = matrix(F, 20, 22, lambda row, column: remainders[column][row])
    kernel = conditions.right_kernel().basis_matrix()
    pairs = []
    for coefficients in kernel.rows():
        AA = sum((coefficients[index]*ambient[index][0] for index in range(22)), R.zero())
        BB = sum((coefficients[index]*ambient[index][1] for index in range(22)), R.zero())
        pairs.append((AA, BB))
    return conditions.rank(), pairs


def child_data(X, Y, Z, pairs):
    S = PolynomialRing(F, "w")
    w = S.gen()
    T = PolynomialRing(S, "u")

    def lift(poly):
        return T([S(value) for value in R(poly).list()])

    AA = lift(pairs[0][0]) + w*lift(pairs[1][0])
    BB = lift(pairs[0][1]) + w*lift(pairs[1][1])
    Xs, Ys, Zs, As = map(lift, (X, Y, Z, A))
    raw = (
        AA**4-6*Xs*AA**2*BB**2+8*Ys*AA*BB**3
        -3*Xs**2*BB**4-4*As*BB**4*Zs**4
    )
    after, remainder = raw.quo_rem(Zs**4)
    if remainder:
        raise ArithmeticError("collision factor Z^4 did not divide")

    # The odd squareclass is the binary quartic.  This is a small bivariate
    # factorization over GF(61), not a surface ideal or elimination problem.
    WU = PolynomialRing(S, "u")
    factorization = WU(after).factor()
    quartic = WU(factorization.unit())
    for factor, exponent in factorization:
        if int(exponent) % 2:
            quartic *= factor
    if quartic.degree() != 4:
        raise ArithmeticError("odd squareclass is not quartic")
    coefficients = quartic.list()+[S.zero()]*(5-len(quartic.list()))
    e, d, c, b, a = coefficients
    invariant_I = 12*a*e-3*b*d+c**2
    invariant_J = 72*a*c*e+9*b*c*d-27*a*d**2-27*b**2*e-2*c**3
    child_A = S(-27*invariant_I)
    child_B = S(-27*invariant_J)

    for factor, unused in child_A.gcd(child_B).factor():
        scale = min(child_A.valuation(factor)//4, child_B.valuation(factor)//6)
        if scale:
            child_A //= factor**(4*scale)
            child_B //= factor**(6*scale)
    delta = S(-16*(4*child_A**3+27*child_B**2))
    root_rank = 0
    fibre_orders = []
    for factor, exponent in delta.factor():
        order_delta = int(exponent)
        order_A = int(child_A.valuation(factor))
        order_B = int(child_B.valuation(factor))
        if order_A == 0 and order_B == 0:
            contribution = (order_delta-1)*factor.degree()
        elif order_delta in (1, 2):
            contribution = 0
        else:
            raise ArithmeticError("unclassified finite fibre")
        root_rank += contribution
        fibre_orders.append((int(factor.degree()), order_A, order_B, order_delta))
    infinity = (8-child_A.degree(), 12-child_B.degree(), 24-delta.degree())
    if infinity[2] > 1 and infinity[0] == infinity[1] == 0:
        root_rank += infinity[2]-1
    return {
        "root_rank": int(root_rank),
        "degrees_A_B_Delta": [int(child_A.degree()), int(child_B.degree()), int(delta.degree())],
        "finite_factor_degree_and_orders_A_B_Delta": fibre_orders,
        "infinity_orders_A_B_Delta": list(map(int, infinity)),
    }


records = candidate_data["search"]["candidates"]
stop = len(records) if args.stop is None else min(args.stop, len(records))
rows = []
counts = Counter()
for candidate_index in range(args.start, stop):
    record = records[candidate_index]
    if record["shape_Xnum_Xden_Ynum_Yden"] != [24, 20, 34, 30]:
        continue
    row = {"candidate_index": candidate_index}
    try:
        X, Y, Z = section_polynomials(record)
        rank, pairs = smooth_basis(X, Y, Z)
        row.update({"condition_rank": int(rank), "h0": len(pairs)})
        if rank != 20 or len(pairs) != 2:
            row["status"] = "REJECTED_RR_DIMENSION"
        else:
            row.update(child_data(X, Y, Z, pairs))
            row["status"] = "PASS_5A1_CHILD" if row["root_rank"] == 5 else "REJECTED_WRONG_ROOT_RANK"
    except (ArithmeticError, ValueError, ZeroDivisionError) as error:
        row.update({"status": "REJECTED_COMPILER_EXCEPTION", "reason": str(error)})
    counts[row["status"]] += 1
    rows.append(row)
    if len(rows) % 100 == 0:
        print(
            "Q4O323Q207RR61PROGRESS|tested={}|five_a1={}|runtime={:.3f}".format(
                len(rows), counts["PASS_5A1_CHILD"], time.monotonic()-started,
            ), flush=True,
        )

winners = [row for row in rows if row["status"] == "PASS_5A1_CHILD"]
payload = {
    "schema": "elkies-k3.h92-q4o323-q207-smooth-rr-filter-mod61.v1",
    "status": "PASS_MOD61_Q4O323_Q207_SMOOTH_RR_FILTER",
    "prime": 61,
    "range": [args.start, stop],
    "target_shape_records_tested": len(rows),
    "status_counts": dict(sorted(counts.items())),
    "five_A1_candidates": winners,
    "all_results": rows,
    "method": {
        "ambient_dimension": 22,
        "collision_condition_count": 20,
        "large_Groebner_required": False,
        "surface_elimination_required": False,
        "runtime_seconds": time.monotonic()-started,
    },
    "proof_boundary": (
        "A 5A1 result is an equation-side modular seed only. It still requires the "
        "independent marked q207 class gate, QQ lift, and exact QQ RR/Jacobian replay."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (CANDIDATES, POINTING)],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in (CANDIDATES, POINTING)},
    },
}
output = args.output if args.output.is_absolute() else ROOT/args.output
output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q4O323Q207RR61|tested={}|five_a1={}|counts={}|runtime={:.3f}|output={}".format(
        len(rows), len(winners), dict(sorted(counts.items())),
        payload["method"]["runtime_seconds"], output,
    ), flush=True,
)
