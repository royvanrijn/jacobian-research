#!/usr/bin/env sage-python
"""Select the reflected q207 seed by a degree-sixteen Abel trace mod 61.

Invert each q323 candidate to its degree-sixteen multisection on q4/o208.
At ordinary parent fibres, compute the sum of its sixteen points from the
unique relation in L(17 O).  Compare that sum with the exact q4/o208 MW word
of the transported q207 class.  This uses univariate quotient rings and
linear algebra only; it does not split the degree-sixteen polynomial and does
not use a Groebner basis.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
CANDIDATES_DEFAULT = LOCAL / "q4o323-q207-four-section-words-mod61.json"
POINTING = LOCAL / "q4o323-component2-pointing-qq.json"
Q323_RR = LOCAL / "q4o208-q4o323-a3-2a2-rr-qq.json"
Q208_MODEL = LOCAL / "q4o208-compact-weierstrass-qq.json"
Q208_BRANCHES = LOCAL / "q4o208-q4o323-horizontal-resolved-qq.json"
Q208_BRANCH_MARKING = LOCAL / "q4o208-q4o323-horizontal-marking-qq.json"
Q208_MARKING = LOCAL / "q24-2a5-physical-q4o208-equation-marking-qq.json"
Q208_PHYSICAL = GENERATED / "elkies-k3-h3-q4o208-physical-3a3-marking.json"
Q323_MODULAR = LOCAL / "q4o208-physical-q4o323-horizontal-mod131.json"
Q207_MARKING = LOCAL / "q4o323-reflected-fixed-suffix-component2-marking.json"
Q323_HORIZONTAL = LOCAL / "q4o208-q4o323-horizontal-by-halving-qq.json"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=61)
parser.add_argument("--fibres", type=int, default=3)
parser.add_argument("--max-candidates", type=int)
parser.add_argument("--candidate-index", type=int)
parser.add_argument("--candidates", type=Path, default=CANDIDATES_DEFAULT)
parser.add_argument("--chord-sign", type=int, choices=(-1, 1), default=1)
parser.add_argument(
    "--output", type=Path,
    default=LOCAL / "q4o323-q207-abel-trace-mod61.json",
)
args = parser.parse_args()
CANDIDATES = args.candidates.resolve()
INPUTS = (
    CANDIDATES, POINTING, Q323_RR, Q208_MODEL, Q208_BRANCHES,
    Q208_BRANCH_MARKING, Q208_MARKING, Q208_PHYSICAL, Q323_MODULAR,
    Q207_MARKING, Q323_HORIZONTAL,
)
if args.fibres < 1:
    raise ValueError("--fibres must be positive")
started = time.monotonic()
prime = ZZ(args.prime)
if not prime.is_prime() or prime in (2, 3):
    raise ValueError("--prime must be an odd prime other than 3")
F = GF(prime)
R = PolynomialRing(F, "u")
u = R.gen()
K = R.fraction_field()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def display_path(path):
    path = path.resolve()
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def reduce_qq(value):
    value = QQ(value)
    if value.denominator() % prime == 0:
        raise ZeroDivisionError("the selected prime divides a stored rational denominator")
    return F(value.numerator()) / F(value.denominator())


def polynomial(values):
    return R([reduce_qq(value) for value in values])


def rational_function(record):
    return K(polynomial(record["numerator_coefficients_low_to_high"])) / K(
        polynomial(record["denominator_coefficients_low_to_high"])
    )


def evaluate_polynomial(poly, value):
    answer = K.zero()
    for coefficient in reversed(list(poly)):
        answer = answer * value + K(coefficient)
    return answer


def reduce_mod_H(value, H):
    value = K(value)
    numerator = R(value.numerator())
    denominator = R(value.denominator())
    if denominator.gcd(H).degree() != 0:
        raise ZeroDivisionError("rational-function denominator meets the fibre")
    return (numerator * denominator.inverse_mod(H)) % H


def newton_power_sums(monic):
    degree = monic.degree()
    assert monic[degree] == 1
    sums = [F(degree)]
    for order in range(1, degree):
        total = F(order) * monic[degree - order]
        for index in range(1, order):
            total += monic[degree - index] * sums[order - index]
        sums.append(-total)
    return sums


candidate_data = json.loads(CANDIDATES.read_text())
pointing = json.loads(POINTING.read_text())
q323_rr = json.loads(Q323_RR.read_text())
q208_model = json.loads(Q208_MODEL.read_text())
q208_branches = json.loads(Q208_BRANCHES.read_text())
q208_branch_marking = json.loads(Q208_BRANCH_MARKING.read_text())
q208_marking = json.loads(Q208_MARKING.read_text())
q208_physical = json.loads(Q208_PHYSICAL.read_text())
q323_modular = json.loads(Q323_MODULAR.read_text())
q207_marking = json.loads(Q207_MARKING.read_text())
q323_horizontal = json.loads(Q323_HORIZONTAL.read_text())
assert candidate_data["status"] in (
    "PASS_MOD61_Q4O323_Q207_FOUR_SECTION_CANDIDATES",
    "PASS_MODP_Q4O323_Q207_UNIQUE_SHELL_ABEL_WORD",
    "PASS_MODP_Q4O323_Q207_SHELL_ABEL_WORDS_WITH_RESIDUAL_AMBIGUITY",
)
assert ZZ(candidate_data["prime"]) == prime
assert pointing["status"] == "PASS_EXACT_QQ_Q4O323_OLD_A11_COMPONENT2_POINTING"
assert q323_rr["status"] == "PASS_EXACT_QQ_Q4O323_A3_2A2_RR_AND_JACOBIAN"
assert q208_model["status"] == "PASS_EXACT_QQ_Q4O208_COMPACT_WEIERSTRASS_NORMALIZATION"
assert q208_branches["status"] == "PASS_EXACT_QQ_Q4O323_RESOLVED_SIMPLE_POLE_HORIZONTAL"
assert q208_branch_marking["status"] == "PASS_EXACT_QQ_Q4O323_LIFTED_SHELL_EXCLUDES_TARGET"
assert q208_marking["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O208_C5_EQUATION_MARKING"
assert q208_physical["status"] == "PASS_EXACT_Q4O208_PHYSICAL_3A3_MARKING"
assert q323_modular["status"] == "PASS_EXACT_Q4O323_POLYNOMIAL_SECTION_SUBGROUP_OBSTRUCTION"
assert q207_marking["status"] == "PASS_EXACT_Q4O323_REFLECTED_FIXED_SUFFIX_MARKING"
assert q323_horizontal["status"] == "PASS_EXACT_QQ_Q4O323_HORIZONTAL_BY_MW_HALVING"

# Parent q4/o208 model and the exact equation points used in the target word.
parent_A = polynomial(q208_model["compact_model"]["A_coefficients_low_to_high"])
parent_B = polynomial(q208_model["compact_model"]["B_coefficients_low_to_high"])
branch_records = {
    int(record["branch_index"]): record
    for record in q208_branches["exact_QQ_horizontal_sections"]
}
needed_branch_ids = (2, 6, 8, 16, 18)
assert all(index in branch_records for index in needed_branch_ids)


def evaluate_qq_coefficients(values, tau):
    """Evaluate over QQ before reduction, allowing coefficientwise p-poles to cancel."""

    tau = QQ(ZZ(tau))
    answer = QQ.zero()
    for coefficient in reversed(values):
        answer = answer * tau + QQ(coefficient)
    return answer


def branch_at(index, tau, curve):
    record = branch_records[index]
    Z_QQ = evaluate_qq_coefficients(record["Z_coefficients_low_to_high"], tau)
    if not Z_QQ:
        raise ZeroDivisionError("projective branch denominator vanished")
    X_QQ = evaluate_qq_coefficients(record["X_coefficients_low_to_high"], tau)
    Y_QQ = evaluate_qq_coefficients(record["Y_coefficients_low_to_high"], tau)
    return curve(reduce_qq(X_QQ / Z_QQ**2), reduce_qq(Y_QQ / Z_QQ**3))


c7_record = q208_model["transported_exact_section"]
target_coefficients = {2: -1, 6: 2, 8: 1, 16: -3, 18: 3}
c7_coefficient = 1
t_coefficient = 7

# Verify the exact marked-lattice origin of the parent MW word after applying
# the C5-child-to-physical-3A3 unimodular bridge.
matched_classes = {
    int(record["stored_branch_index"]): vector(ZZ, record["matched_NS_coordinates"])
    for record in q208_branch_marking["lattice_match"]["branch_matches"]
}
equation_to_physical = matrix(ZZ, q208_physical["C5_child_basis_in_physical_3A3"])
matched_classes = {
    index: section * equation_to_physical
    for index, section in matched_classes.items()
}
c7_class = vector(
    ZZ, q208_marking["old_A11_component_7_on_C5_pointed_child"]["NS_coordinates"]
) * equation_to_physical
t_class = vector(ZZ, q323_modular["target"]["NS_coordinates"]) * equation_to_physical
target_child = vector(
    ZZ,
    q207_marking["fixed_suffix_transport"]["q207_component_reduction"]
    ["equation_preflight"]["horizontal_section"],
)
child_in_source = matrix(ZZ, q207_marking["basis_in_source"])
target_source = target_child * child_in_source
target_word_class = sum(
    (coefficient * matched_classes[index] for index, coefficient in target_coefficients.items()),
    vector(ZZ, 19),
) + c7_coefficient * c7_class + t_coefficient * t_class
parent_fibre = vector(QQ, [1, 0] + [0] * 17)
parent_zero = vector(QQ, [-1, 1] + [0] * 17)
parent_roots = [
    vector(QQ, [0, 0] + [-int(index == other) for other in range(17)])
    for index in range(9)
]
parent_quotient_columns = matrix(
    QQ, [parent_fibre, parent_zero] + parent_roots
).right_kernel().basis_matrix().transpose()
target_signature = target_source * parent_quotient_columns
word_signature = target_word_class * parent_quotient_columns
assert target_signature == word_signature


def expected_target_at(tau):
    Atau, Btau = parent_A(tau), parent_B(tau)
    curve = EllipticCurve(F, [0, 0, 0, Atau, Btau])
    if not curve.discriminant():
        raise ZeroDivisionError("singular parent fibre")
    answer = curve(0)
    for index, coefficient in target_coefficients.items():
        answer += coefficient * branch_at(index, tau, curve)
    c7_x = reduce_qq(evaluate_qq_coefficients(c7_record["x_coefficients_low_to_high"], tau))
    c7_y = reduce_qq(evaluate_qq_coefficients(c7_record["y_coefficients_low_to_high"], tau))
    answer += c7_coefficient * curve(c7_x, c7_y)
    T = q323_horizontal["exact_QQ_horizontal"]
    t_x = reduce_qq(
        evaluate_qq_coefficients(T["x"]["numerator_coefficients_low_to_high"], tau)
        / evaluate_qq_coefficients(T["x"]["denominator_coefficients_low_to_high"], tau)
    )
    t_y = reduce_qq(
        evaluate_qq_coefficients(T["y"]["numerator_coefficients_low_to_high"], tau)
        / evaluate_qq_coefficients(T["y"]["denominator_coefficients_low_to_high"], tau)
    )
    answer += t_coefficient * curve(t_x, t_y)
    if answer.is_zero():
        raise ZeroDivisionError("target word specializes to zero")
    return curve, answer


# q323 child and inverse pointed-quartic data.
child_A = polynomial(pointing["global_short_model"]["A_coefficients_low_to_high"])
child_B = polynomial(pointing["global_short_model"]["B_coefficients_low_to_high"])
a1, a2, a3, unused_a4, unused_a6 = [
    rational_function(pointing["pointed_generalized_weierstrass"][name])
    for name in ("a1", "a2", "a3", "a4", "a6")
]
b2 = a1**2 + 4 * a2
w0 = K(polynomial(pointing["quartic_square_at_t0"]["L0_coefficients_low_to_high"]))
d_coefficient = a1 * w0
quartic = [
    K(polynomial(values))
    for values in q323_rr["quartic"]["coefficients_in_t_low_to_high"]
]
rr_pairs = [
    (polynomial(record["a_coefficients_low_to_high"]), reduce_qq(record["b"]))
    for record in pointing["quartic_map"]["rr_pairs"]
]
(aa0, bb0), (aa1, bb1) = rr_pairs
square_factor = [
    rational_function(record)
    for record in pointing["quartic_map"]["square_factor_coefficients_in_t_low_to_high"]
]
unit_root = rational_function(pointing["quartic_map"]["normalizing_unit_root"])

horizontal = q208_model["transported_exact_section"]
# The transported C7 section is not the q323 horizontal.  Load that horizontal
# from the exact halving artifact used by the compiler.
halving = json.loads(Q323_HORIZONTAL.read_text())
horizontal = halving["exact_QQ_horizontal"]
horizontal_x = rational_function(horizontal["x"])
horizontal_y = rational_function(horizontal["y"])
horizontal_Z = polynomial(horizontal["x"]["denominator_coefficients_low_to_high"])
assert horizontal_Z.is_square()
horizontal_Z = horizontal_Z.sqrt()


def invert_child_candidate(record):
    X = rational_function(record["x"])
    Y = rational_function(record["y"])
    assert Y**2 == X**3 + K(child_A) * X + K(child_B)
    x_general = X / 9 - b2 / 12
    y_general = Y / 27 - (a1 * x_general + a3) / 2
    if not y_general:
        raise ZeroDivisionError("pointed inverse has zero generalized ordinate")
    old_t = K(2 * w0 * (x_general + a2) / y_general)
    ordinate = K((x_general * old_t**2 - d_coefficient * old_t) / (2 * w0) - w0)
    assert ordinate**2 == sum(quartic[index] * old_t**index for index in range(5))
    return old_t, ordinate


def abel_trace_at(record, tau, chord_sign):
    curve, expected = expected_target_at(tau)
    Atau, Btau = parent_A(tau), parent_B(tau)
    old_t, ordinate = invert_child_candidate(record)
    H = R(old_t.numerator() - tau * old_t.denominator())
    if H.degree() != 16:
        raise ZeroDivisionError("degree-sixteen fibre dropped degree")
    H = H.monic()
    if H.gcd(H.derivative()).degree() != 0:
        raise ZeroDivisionError("degree-sixteen fibre is not etale")

    ordinate_A = reduce_mod_H(ordinate, H)
    aa = R(aa0(tau)) + u * R(aa1(tau))
    bb = R(bb0) + u * R(bb1)
    Ztau = horizontal_Z(tau)
    if not Ztau:
        raise ZeroDivisionError("q323 horizontal pole fibre")
    slope = reduce_mod_H(K(aa) / K(bb * Ztau), H)
    square_tau = reduce_mod_H(
        sum(square_factor[index] * K(tau)**index for index in range(len(square_factor))), H
    )
    unit_A = reduce_mod_H(unit_root, H)
    bb_A = reduce_mod_H(K(bb), H)
    if not bb_A.gcd(H).degree() == 0:
        raise ZeroDivisionError("RR slope denominator meets fibre")
    branch = (
        F(chord_sign) * ordinate_A * square_tau * unit_A
        * bb_A.inverse_mod(H)**2
    ) % H
    hx, hy = horizontal_x(tau), horizontal_y(tau)
    xA = ((slope**2 - F(hx) + branch) / 2) % H
    yA = (slope * (xA - F(hx)) + F(hy)) % H
    if (yA**2 - xA**3 - Atau * xA - Btau) % H:
        raise ArithmeticError("inverse chord misses the parent fibre")

    # L(17O) = <1,x,...,x^8,y,xy,...,x^7 y>.
    powers = [R.one()]
    for unused in range(8):
        powers.append((powers[-1] * xA) % H)
    columns = powers + [(yA * powers[index]) % H for index in range(8)]
    evaluation = matrix(F, 16, 17, lambda row, column: columns[column][row])
    kernel = evaluation.right_kernel_matrix()
    if kernel.nrows() != 1:
        raise ArithmeticError("L(17O) evaluation kernel is not one-dimensional")
    relation = kernel[0]
    SX = PolynomialRing(F, "X")
    XX = SX.gen()
    Afun = sum(relation[index] * XX**index for index in range(9))
    Bfun = sum(relation[9 + index] * XX**index for index in range(8))
    intersection = Afun**2 - (XX**3 + Atau * XX + Btau) * Bfun**2
    if intersection.degree() != 17:
        raise ArithmeticError("L(17O) intersection degree is not seventeen")
    root_sum = -intersection[16] / intersection[17]
    power_sums = newton_power_sums(H)
    divisor_x_sum = sum(xA[index] * power_sums[index] for index in range(16))
    residual_x = root_sum - divisor_x_sum
    if not Bfun(residual_x):
        raise ArithmeticError("residual point has zero B coefficient")
    residual_y = -Afun(residual_x) / Bfun(residual_x)
    trace = -curve(residual_x, residual_y)
    return trace, expected, H


all_candidates = candidate_data["search"].get(
    "candidates", candidate_data["search"].get("winners", [])
)
if args.candidate_index is not None:
    if not 0 <= args.candidate_index < len(all_candidates):
        raise IndexError("--candidate-index is outside the stored candidate list")
    indexed_candidates = [(args.candidate_index, all_candidates[args.candidate_index])]
elif args.max_candidates is not None:
    indexed_candidates = list(enumerate(all_candidates[:args.max_candidates]))
else:
    indexed_candidates = list(enumerate(all_candidates))
rows = []
survivors = []
for candidate_index, record in indexed_candidates:
    trials = []
    failures = []
    for integer_tau in range(int(prime)):
        if len(trials) >= args.fibres:
            break
        tau = F(integer_tau)
        try:
            trace, expected, H = abel_trace_at(record, tau, args.chord_sign)
        except (ArithmeticError, ZeroDivisionError, ValueError) as error:
            failures.append({"old_base_value": integer_tau, "reason": str(error)})
            continue
        matched = trace == expected
        trials.append({
            "old_base_value": integer_tau,
            "degree": int(H.degree()),
            "trace": [int(value) for value in trace.xy()] if not trace.is_zero() else None,
            "expected": [int(value) for value in expected.xy()],
            "literal_target_identity": bool(matched),
        })
        if not matched:
            break
    passed = len(trials) == args.fibres and all(row["literal_target_identity"] for row in trials)
    result = {
        "candidate_index": candidate_index,
        "shape_Xnum_Xden_Ynum_Yden": record["shape_Xnum_Xden_Ynum_Yden"],
        "records": record.get(
            "records", [{"shell_indices": record.get("shell_indices", [])}],
        ),
        "ordinary_fibre_trials": trials,
        "skipped_fibres_before_decision": failures,
        "passes_all_Abel_trace_identities": passed,
    }
    rows.append(result)
    if passed:
        survivors.append(result)
    if (candidate_index + 1) % 250 == 0:
        print(
            "Q4O323Q207ABELPROGRESS|tested={}|survivors={}|runtime={:.3f}".format(
                candidate_index + 1, len(survivors), time.monotonic() - started
            ),
            flush=True,
        )

complete = (
    args.candidate_index is None
    and (args.max_candidates is None or args.max_candidates >= len(all_candidates))
)
payload = {
    "schema": "elkies-k3.h92-q4o323-q207-abel-trace-modp.v1",
    "status": (
        "PASS_MODP_Q4O323_Q207_UNIQUE_ABEL_TRACE_SEED"
        if complete and len(survivors) == 1 else
        (
            "REJECTED_MODP_Q4O323_Q207_NO_ABEL_TRACE_SEED"
            if complete and not survivors else
            "PASS_MODP_Q4O323_Q207_ABEL_TRACE_WITH_RESIDUAL_AMBIGUITY"
        )
    ),
    "prime": int(prime),
    "chord_sign_relative_to_stored_quartic_ordinate": args.chord_sign,
    "required_ordinary_fibres_per_candidate": args.fibres,
    "target_parent_MW_word": {
        "word": "-P2+2*P6+P8-3*P16+3*P18+C7+7*T",
        "target_quotient_signature": [str(value) for value in target_signature],
        "word_quotient_signature": [str(value) for value in word_signature],
        "exact_marked_quotient_identity": True,
    },
    "search": {
        "available_candidates": len(all_candidates),
        "tested_candidates": len(rows),
        "complete_candidate_scan": complete,
        "survivor_count": len(survivors),
        "survivors": survivors,
        "all_results": rows,
    },
    "method": {
        "large_Groebner_required": False,
        "multisection_splitting_required": False,
        "elimination_required": False,
        "linear_algebra": "16x17 evaluation kernel in L(17O)",
        "runtime_seconds": time.monotonic() - started,
    },
    "proof_boundary": (
        "This is an exact good-prime construction fingerprint for the q207 equation-side seed. "
        "A surviving modular point still requires a characteristic-zero lift, exact resolved "
        "Riemann--Roch plane, quartic/Jacobian, fibres, and next marking."
    ),
    "inputs": {
        "paths": [display_path(path) for path in INPUTS],
        "sha256": {display_path(path): sha256(path) for path in INPUTS},
    },
}
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4O323Q207ABEL|prime={}|tested={}|survivors={}|fibres={}|sign={}|runtime={:.3f}|"
    "status={}|output={}".format(
        prime, len(rows), len(survivors), args.fibres, args.chord_sign,
        payload["method"]["runtime_seconds"], payload["status"], args.output,
    ),
    flush=True,
)
