#!/usr/bin/env sage-python
"""Abel-name a complete q323 polynomial P.O=0 shell modulo a good prime.

Each shell section is inverted to its degree-one-through-seven multisection on
q4/o208.  Its Abel trace is computed by a quotient-ring L((d+1)O) kernel and
matched against all 258 exact lattice P.O=0 classes expressed in a fixed
unimodular equation-side MW basis.  No multisection splitting, elimination,
or Groebner basis is used.
"""

import argparse
import hashlib
import itertools
import json
import time
from functools import lru_cache
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
LATTICE = LOCAL / "q4o323-p0-shell-anchor-domains-mod61.json"
POINTING = LOCAL / "q4o323-component2-pointing-qq.json"
Q323_RR = LOCAL / "q4o208-q4o323-a3-2a2-rr-qq.json"
Q208_MODEL = LOCAL / "q4o208-compact-weierstrass-qq.json"
Q208_BRANCHES = LOCAL / "q4o208-q4o323-horizontal-resolved-qq.json"
Q208_BRANCH_MARKING = LOCAL / "q4o208-q4o323-horizontal-marking-qq.json"
Q208_MARKING = LOCAL / "q24-2a5-physical-q4o208-equation-marking-qq.json"
Q208_PHYSICAL = GENERATED / "elkies-k3-h3-q4o208-physical-3a3-marking.json"
Q323_MODULAR_MARKING = LOCAL / "q4o208-physical-q4o323-horizontal-mod131.json"
Q207_MARKING = LOCAL / "q4o323-reflected-fixed-suffix-component2-marking.json"
Q323_HORIZONTAL = LOCAL / "q4o208-q4o323-horizontal-by-halving-qq.json"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=61)
parser.add_argument("--fibres", type=int, default=3)
parser.add_argument("--max-sections", type=int)
parser.add_argument("--shell", type=Path)
parser.add_argument("--output", type=Path)
args = parser.parse_args()
if args.fibres < 1:
    raise ValueError("--fibres must be positive")
started = time.monotonic()
prime = ZZ(args.prime)
if not prime.is_prime() or prime in (2, 3):
    raise ValueError("--prime must be an odd prime other than 3")
SHELL = args.shell or LOCAL / f"q4o323-p0-shell-mod{prime}.json"
OUTPUT = args.output or LOCAL / f"q4o323-p0-shell-abel-names-mod{prime}.json"
INPUTS = (
    SHELL, LATTICE, POINTING, Q323_RR, Q208_MODEL, Q208_BRANCHES,
    Q208_BRANCH_MARKING, Q208_MARKING, Q323_MODULAR_MARKING,
    Q208_PHYSICAL, Q207_MARKING, Q323_HORIZONTAL,
)
F = GF(prime)
R = PolynomialRing(F, "u")
u = R.gen()
K = R.fraction_field()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


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


def evaluate_qq_coefficients(values, tau):
    tau = QQ(ZZ(tau))
    answer = QQ.zero()
    for coefficient in reversed(values):
        answer = answer * tau + QQ(coefficient)
    return answer


def branch_at(record, tau, curve):
    Z = evaluate_qq_coefficients(record["Z_coefficients_low_to_high"], tau)
    if not Z:
        raise ZeroDivisionError("projective section denominator vanished")
    X = evaluate_qq_coefficients(record["X_coefficients_low_to_high"], tau)
    Y = evaluate_qq_coefficients(record["Y_coefficients_low_to_high"], tau)
    return curve(reduce_qq(X / Z**2), reduce_qq(Y / Z**3))


def reduce_mod_H(value, H):
    value = K(value)
    numerator = R(value.numerator())
    denominator = R(value.denominator())
    if denominator.gcd(H).degree() != 0:
        raise ZeroDivisionError("rational-function denominator meets the fibre")
    return (numerator * denominator.inverse_mod(H)) % H


def newton_power_sums(monic):
    degree = monic.degree()
    sums = [F(degree)]
    for order in range(1, degree):
        total = F(order) * monic[degree - order]
        for index in range(1, order):
            total += monic[degree - index] * sums[order - index]
        sums.append(-total)
    return sums


shell_data = json.loads(SHELL.read_text())
lattice_data = json.loads(LATTICE.read_text())
pointing = json.loads(POINTING.read_text())
q323_rr = json.loads(Q323_RR.read_text())
q208_model = json.loads(Q208_MODEL.read_text())
q208_branches = json.loads(Q208_BRANCHES.read_text())
branch_marking = json.loads(Q208_BRANCH_MARKING.read_text())
q208_marking = json.loads(Q208_MARKING.read_text())
q208_physical = json.loads(Q208_PHYSICAL.read_text())
q323_modular_marking = json.loads(Q323_MODULAR_MARKING.read_text())
q207_marking = json.loads(Q207_MARKING.read_text())
horizontal_data = json.loads(Q323_HORIZONTAL.read_text())
assert shell_data["status"] == "PASS_MODP_Q4O323_COMPLETE_POLYNOMIAL_P0_SHELL"
assert ZZ(shell_data["prime"]) == prime
assert lattice_data["status"] == "PASS_MOD61_Q4O323_REGULAR_P0_SHELL_ANCHOR_DOMAINS"
assert pointing["status"] == "PASS_EXACT_QQ_Q4O323_OLD_A11_COMPONENT2_POINTING"
assert q323_rr["status"] == "PASS_EXACT_QQ_Q4O323_A3_2A2_RR_AND_JACOBIAN"
assert q208_model["status"] == "PASS_EXACT_QQ_Q4O208_COMPACT_WEIERSTRASS_NORMALIZATION"
assert q208_branches["status"] == "PASS_EXACT_QQ_Q4O323_RESOLVED_SIMPLE_POLE_HORIZONTAL"
assert branch_marking["status"] == "PASS_EXACT_QQ_Q4O323_LIFTED_SHELL_EXCLUDES_TARGET"
assert q208_marking["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O208_C5_EQUATION_MARKING"
assert q208_physical["status"] == "PASS_EXACT_Q4O208_PHYSICAL_3A3_MARKING"
assert q323_modular_marking["status"] == "PASS_EXACT_Q4O323_POLYNOMIAL_SECTION_SUBGROUP_OBSTRUCTION"
assert q207_marking["status"] == "PASS_EXACT_Q4O323_REFLECTED_FIXED_SUFFIX_MARKING"
assert horizontal_data["status"] == "PASS_EXACT_QQ_Q4O323_HORIZONTAL_BY_MW_HALVING"

parent_A = polynomial(q208_model["compact_model"]["A_coefficients_low_to_high"])
parent_B = polynomial(q208_model["compact_model"]["B_coefficients_low_to_high"])
child_A = polynomial(pointing["global_short_model"]["A_coefficients_low_to_high"])
child_B = polynomial(pointing["global_short_model"]["B_coefficients_low_to_high"])
branch_records = {
    int(record["branch_index"]): record
    for record in q208_branches["exact_QQ_horizontal_sections"]
}
matched_classes = {
    int(record["stored_branch_index"]): vector(ZZ, record["matched_NS_coordinates"])
    for record in branch_marking["lattice_match"]["branch_matches"]
}

# Compute the actual q4/o208 MW quotient, rather than reading a coordinate
# tail from its non-root-block-adapted frame.
parent_fibre = vector(QQ, [1, 0] + [0] * 17)
parent_zero = vector(QQ, [-1, 1] + [0] * 17)
parent_roots = [
    vector(QQ, [0, 0] + [-int(index == other) for other in range(17)])
    for index in range(9)
]
parent_trivial = matrix(QQ, [parent_fibre, parent_zero] + parent_roots)
assert parent_trivial.rank() == 11
parent_quotient_columns = parent_trivial.right_kernel().basis_matrix().transpose()

# The simple-pole shell alone has index two in the exact parent MW quotient.
# Adding the exact C7 section and the halved q323 horizontal T supplies a
# determinant-one equation-side basis.
matched_classes["C7"] = vector(
    ZZ, q208_marking["old_A11_component_7_on_C5_pointed_child"]["NS_coordinates"]
)
matched_classes["T"] = vector(ZZ, q323_modular_marking["target"]["NS_coordinates"])
equation_to_physical = matrix(ZZ, q208_physical["C5_child_basis_in_physical_3A3"])
matched_classes = {
    label: section * equation_to_physical
    for label, section in matched_classes.items()
}
basis_labels = (2, 6, 8, 16, 18, 23, "C7", "T")
graph_solutions = branch_marking["lattice_match"]["all_graph_solution_branch_matches"]
for label in basis_labels[:6]:
    tails = {
        tuple(record["matched_NS_coordinates"])
        for solution in graph_solutions
        for record in solution["branch_matches"]
        if int(record["stored_branch_index"]) == label
    }
    assert len(tails) == 1
mw_basis_matrix = matrix(
    ZZ, [matched_classes[label] * parent_quotient_columns for label in basis_labels]
)
assert mw_basis_matrix.det() == -1
mw_basis_inverse = mw_basis_matrix.inverse()

child_in_source = matrix(ZZ, q207_marking["basis_in_source"])
lattice_sections = [
    vector(ZZ, record["NS_coordinates"])
    for record in lattice_data["lattice"]["sections"]
]
lattice_source_signatures = [
    section * child_in_source * parent_quotient_columns for section in lattice_sections
]
lattice_basis_words = []
for signature in lattice_source_signatures:
    coefficients = vector(ZZ, signature * mw_basis_inverse)
    assert coefficients * mw_basis_matrix == signature
    lattice_basis_words.append(coefficients)


@lru_cache(maxsize=None)
def expected_fingerprints(integer_tau):
    tau = F(integer_tau)
    Atau, Btau = parent_A(tau), parent_B(tau)
    curve = EllipticCurve(F, [0, 0, 0, Atau, Btau])
    if not curve.discriminant():
        raise ZeroDivisionError("singular parent fibre")
    points = [branch_at(branch_records[label], tau, curve) for label in basis_labels[:6]]
    c7_record = q208_model["transported_exact_section"]
    c7_x = reduce_qq(evaluate_qq_coefficients(c7_record["x_coefficients_low_to_high"], tau))
    c7_y = reduce_qq(evaluate_qq_coefficients(c7_record["y_coefficients_low_to_high"], tau))
    points.append(curve(c7_x, c7_y))
    horizontal_record = horizontal_data["exact_QQ_horizontal"]
    t_x = reduce_qq(
        evaluate_qq_coefficients(horizontal_record["x"]["numerator_coefficients_low_to_high"], tau)
        / evaluate_qq_coefficients(horizontal_record["x"]["denominator_coefficients_low_to_high"], tau)
    )
    t_y = reduce_qq(
        evaluate_qq_coefficients(horizontal_record["y"]["numerator_coefficients_low_to_high"], tau)
        / evaluate_qq_coefficients(horizontal_record["y"]["denominator_coefficients_low_to_high"], tau)
    )
    points.append(curve(t_x, t_y))
    rows = {}
    for lattice_index, coefficients in enumerate(lattice_basis_words):
        point = curve(0)
        for coefficient, basis_point in zip(coefficients, points):
            point += coefficient * basis_point
        key = None if point.is_zero() else tuple(map(int, point.xy()))
        rows.setdefault(key, []).append(lattice_index)
    return curve, rows


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
horizontal = horizontal_data["exact_QQ_horizontal"]
horizontal_x = rational_function(horizontal["x"])
horizontal_y = rational_function(horizontal["y"])
horizontal_Z = polynomial(horizontal["x"]["denominator_coefficients_low_to_high"])
assert horizontal_Z.is_square()
horizontal_Z = horizontal_Z.sqrt()


def invert_shell_section(record):
    X = K(R(record["x_coefficients_low_to_high"]))
    Y = K(R(record["y_coefficients_low_to_high"]))
    assert Y**2 == X**3 + K(child_A) * X + K(child_B)
    x_general = X / 9 - b2 / 12
    y_general = Y / 27 - (a1 * x_general + a3) / 2
    if not y_general:
        # This is the old-base infinity branch, hence vertical for q4/o208 and
        # zero in its generic-fibre MW quotient.
        return None, None, 0
    old_t = K(2 * w0 * (x_general + a2) / y_general)
    ordinate = K((x_general * old_t**2 - d_coefficient * old_t) / (2 * w0) - w0)
    assert ordinate**2 == sum(quartic[index] * old_t**index for index in range(5))
    degree = max(old_t.numerator().degree(), old_t.denominator().degree())
    assert 0 <= degree <= 7
    return old_t, ordinate, degree


def abel_trace_at(inverse_data, tau):
    old_t, ordinate, expected_degree = inverse_data
    curve, unused_expected = expected_fingerprints(int(tau))
    Atau, Btau = parent_A(tau), parent_B(tau)
    H = R(old_t.numerator() - tau * old_t.denominator())
    if H.degree() != expected_degree:
        raise ZeroDivisionError("inverse parent degree dropped")
    H = H.monic()
    if H.gcd(H.derivative()).degree() != 0:
        raise ZeroDivisionError("inverse parent fibre is not etale")
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
    bb_A = reduce_mod_H(K(bb), H)
    if bb_A.gcd(H).degree() != 0:
        raise ZeroDivisionError("RR slope denominator meets fibre")
    branch = (
        ordinate_A * square_tau * reduce_mod_H(unit_root, H)
        * bb_A.inverse_mod(H)**2
    ) % H
    hx, hy = horizontal_x(tau), horizontal_y(tau)
    xA = ((slope**2 - F(hx) + branch) / 2) % H
    yA = (slope * (xA - F(hx)) + F(hy)) % H
    if (yA**2 - xA**3 - Atau * xA - Btau) % H:
        raise ArithmeticError("inverse chord misses parent fibre")
    degree = H.degree()
    if degree == 1:
        return curve(xA[0], yA[0]), H

    pole_bound = degree + 1
    max_x_power = pole_bound // 2
    max_yx_power = (pole_bound - 3) // 2
    powers = [R.one()]
    for unused in range(max_x_power):
        powers.append((powers[-1] * xA) % H)
    columns = powers + [(yA * powers[index]) % H for index in range(max_yx_power + 1)]
    assert len(columns) == degree + 1
    evaluation = matrix(F, degree, degree + 1, lambda row, column: columns[column][row])
    kernel = evaluation.right_kernel_matrix()
    if kernel.nrows() != 1:
        raise ArithmeticError("Abel evaluation kernel is not one-dimensional")
    relation = kernel[0]
    SX = PolynomialRing(F, "X")
    XX = SX.gen()
    Afun = sum(relation[index] * XX**index for index in range(max_x_power + 1))
    offset = max_x_power + 1
    Bfun = sum(relation[offset + index] * XX**index for index in range(max_yx_power + 1))
    intersection = Afun**2 - (XX**3 + Atau * XX + Btau) * Bfun**2
    if intersection.degree() != degree + 1:
        raise ArithmeticError("Abel intersection has the wrong degree")
    root_sum = -intersection[degree] / intersection[degree + 1]
    power_sums = newton_power_sums(H)
    divisor_x_sum = sum(xA[index] * power_sums[index] for index in range(degree))
    residual_x = root_sum - divisor_x_sum
    if not Bfun(residual_x):
        raise ArithmeticError("residual point has zero B coefficient")
    residual_y = -Afun(residual_x) / Bfun(residual_x)
    return -curve(residual_x, residual_y), H


shell_records = shell_data["shell"]["records"]
if args.max_sections is not None:
    shell_records = shell_records[:args.max_sections]
results = []
named = []
zero_signature_indices = {
    index for index, signature in enumerate(lattice_source_signatures)
    if not any(signature)
}
for shell_index, record in enumerate(shell_records):
    inverse_data = invert_shell_section(record)
    possible = (
        set(zero_signature_indices)
        if inverse_data[2] == 0 else
        set(range(len(lattice_sections)))
    )
    trials = []
    skipped = []
    for integer_tau in (() if inverse_data[2] == 0 else range(int(prime))):
        if len(trials) >= args.fibres:
            break
        try:
            trace, H = abel_trace_at(inverse_data, F(integer_tau))
            unused_curve, expected = expected_fingerprints(integer_tau)
        except (ArithmeticError, ZeroDivisionError, ValueError) as error:
            skipped.append({"old_base_value": integer_tau, "reason": str(error)})
            continue
        key = None if trace.is_zero() else tuple(map(int, trace.xy()))
        matches_here = set(expected.get(key, ()))
        possible.intersection_update(matches_here)
        trials.append({
            "old_base_value": integer_tau,
            "inverse_degree": int(H.degree()),
            "trace": None if key is None else list(key),
            "matching_lattice_indices_after": sorted(possible),
        })
    row = {
        "shell_index": shell_index,
        "inverse_parent_degree": int(inverse_data[2]),
        "ordinary_fibre_trials": trials,
        "skipped_fibres": skipped,
        "matching_lattice_indices": sorted(possible),
    }
    results.append(row)
    if possible:
        named.append(row)

complete = args.max_sections is None or args.max_sections >= len(
    shell_data["shell"]["records"]
)
unique_names = [row for row in named if len(row["matching_lattice_indices"]) == 1]
payload = {
    "schema": "elkies-k3.h92-q4o323-p0-shell-abel-names-modp.v1",
    "status": (
        "PASS_MODP_Q4O323_P0_SHELL_ABEL_NAMES"
        if complete and named else
        "PARTIAL_MODP_Q4O323_P0_SHELL_ABEL_NAMES"
    ),
    "prime": int(prime),
    "required_ordinary_fibres_per_section": args.fibres,
    "equation_MW_basis": {
        "branch_labels": list(basis_labels),
        "quotient_signature_matrix": [list(map(int, row)) for row in mw_basis_matrix.rows()],
        "determinant": int(mw_basis_matrix.det()),
        "fixed_in_all_complete_graph_solutions": True,
    },
    "classification": {
        "available_shell_sections": len(shell_data["shell"]["records"]),
        "tested_shell_sections": len(results),
        "complete_shell_scan": complete,
        "named_shell_sections": len(named),
        "uniquely_named_shell_sections": len(unique_names),
        "ambiguous_named_shell_sections": len(named) - len(unique_names),
        "unmatched_shell_sections": len(results) - len(named),
        "records": results,
    },
    "method": {
        "large_Groebner_required": False,
        "multisection_splitting_required": False,
        "elimination_required": False,
        "maximum_inverse_parent_degree": 7,
        "construction": "unimodular parent MW words plus quotient-ring Abel traces",
        "runtime_seconds": time.monotonic() - started,
    },
    "proof_boundary": (
        "This Abel-names reductions of exact q323 lattice classes inside the complete mod-p "
        "polynomial shell. It is a modular marking certificate, not a QQ section lift or an "
        "outgoing q12 Riemann--Roch/Jacobian certificate."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in INPUTS},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4O323P0ABEL|prime={}|tested={}|named={}|unique={}|unmatched={}|fibres={}|runtime={:.3f}|"
    "status={}|output={}".format(
        prime, len(results), len(named), len(unique_names), len(results) - len(named), args.fibres,
        payload["method"]["runtime_seconds"], payload["status"], OUTPUT,
    ),
    flush=True,
)
