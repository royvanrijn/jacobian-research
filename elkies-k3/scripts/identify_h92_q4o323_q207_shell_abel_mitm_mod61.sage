#!/usr/bin/env sage-python
"""Find q207 as a four-point sum in a complete q323 shell modulo p.

Use the three-fibre Abel fingerprints of all rational polynomial P.O=0
sections.  A pair/pair meet-in-the-middle finds every four-point sum with the
exact transported q207 parent MW class.  Candidate sums are then checked on
the q323 function field for P.O=10 and inverse parent degree 16.  No splitting,
elimination, or Groebner basis is used.
"""

import argparse
import hashlib
import itertools
import json
import time
from collections import defaultdict
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
POINTING = LOCAL / "q4o323-component2-pointing-qq.json"
Q208_MODEL = LOCAL / "q4o208-compact-weierstrass-qq.json"
Q208_BRANCHES = LOCAL / "q4o208-q4o323-horizontal-resolved-qq.json"
Q323_HORIZONTAL = LOCAL / "q4o208-q4o323-horizontal-by-halving-qq.json"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=61)
parser.add_argument("--abel", type=Path)
parser.add_argument("--shell", type=Path)
parser.add_argument("--output", type=Path)
parser.add_argument(
    "--taus", type=int, nargs="+", default=(2, 3, 4),
    help="two or three ordinary parent-fibre values used for Abel fingerprints",
)
args = parser.parse_args()
started = time.monotonic()
prime = ZZ(args.prime)
if not prime.is_prime() or prime in (2, 3):
    raise ValueError("--prime must be an odd prime other than 3")
ABEL = args.abel or LOCAL / f"q4o323-p0-shell-abel-names-mod{prime}.json"
SHELL = args.shell or LOCAL / f"q4o323-p0-shell-mod{prime}.json"
OUTPUT = args.output or LOCAL / f"q4o323-q207-shell-abel-mitm-mod{prime}.json"
INPUTS = (ABEL, SHELL, POINTING, Q208_MODEL, Q208_BRANCHES, Q323_HORIZONTAL)
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


def evaluate_qq(values, tau):
    answer = QQ.zero()
    for coefficient in reversed(values):
        answer = answer * QQ(tau) + QQ(coefficient)
    return answer


def branch_at(record, tau, curve):
    Z = evaluate_qq(record["Z_coefficients_low_to_high"], tau)
    X = evaluate_qq(record["X_coefficients_low_to_high"], tau)
    Y = evaluate_qq(record["Y_coefficients_low_to_high"], tau)
    return curve(reduce_qq(X / Z**2), reduce_qq(Y / Z**3))


def point_key(point):
    return None if point.is_zero() else tuple(map(int, point.xy()))


abel = json.loads(ABEL.read_text())
shell = json.loads(SHELL.read_text())
pointing = json.loads(POINTING.read_text())
model = json.loads(Q208_MODEL.read_text())
branches = json.loads(Q208_BRANCHES.read_text())
horizontal = json.loads(Q323_HORIZONTAL.read_text())
assert abel["status"] == "PASS_MODP_Q4O323_P0_SHELL_ABEL_NAMES"
assert shell["status"] == "PASS_MODP_Q4O323_COMPLETE_POLYNOMIAL_P0_SHELL"
assert ZZ(abel["prime"]) == ZZ(shell["prime"]) == prime
assert pointing["status"] == "PASS_EXACT_QQ_Q4O323_OLD_A11_COMPONENT2_POINTING"
assert model["status"] == "PASS_EXACT_QQ_Q4O208_COMPACT_WEIERSTRASS_NORMALIZATION"
assert branches["status"] == "PASS_EXACT_QQ_Q4O323_RESOLVED_SIMPLE_POLE_HORIZONTAL"
assert horizontal["status"] == "PASS_EXACT_QQ_Q4O323_HORIZONTAL_BY_MW_HALVING"

taus = tuple(args.taus)
if len(taus) not in (2, 3) or len(set(taus)) != len(taus):
    raise ValueError("--taus must contain two or three distinct values")
parent_A = model["compact_model"]["A_coefficients_low_to_high"]
parent_B = model["compact_model"]["B_coefficients_low_to_high"]
curves = {}
for tau in taus:
    A = reduce_qq(evaluate_qq(parent_A, tau))
    B = reduce_qq(evaluate_qq(parent_B, tau))
    curves[tau] = EllipticCurve(F, [0, 0, 0, A, B])
    assert curves[tau].discriminant()

branch_records = {
    int(record["branch_index"]): record
    for record in branches["exact_QQ_horizontal_sections"]
}
basis_labels = (2, 6, 8, 16, 18, 23, "C7", "T")
target_coefficients = (-1, 2, 1, -3, 3, 0, 1, 7)


def basis_points_at(tau):
    curve = curves[tau]
    answer = [branch_at(branch_records[label], tau, curve) for label in basis_labels[:6]]
    c7 = model["transported_exact_section"]
    answer.append(curve(
        reduce_qq(evaluate_qq(c7["x_coefficients_low_to_high"], tau)),
        reduce_qq(evaluate_qq(c7["y_coefficients_low_to_high"], tau)),
    ))
    T = horizontal["exact_QQ_horizontal"]
    answer.append(curve(
        reduce_qq(
            evaluate_qq(T["x"]["numerator_coefficients_low_to_high"], tau)
            / evaluate_qq(T["x"]["denominator_coefficients_low_to_high"], tau)
        ),
        reduce_qq(
            evaluate_qq(T["y"]["numerator_coefficients_low_to_high"], tau)
            / evaluate_qq(T["y"]["denominator_coefficients_low_to_high"], tau)
        ),
    ))
    return answer


target_points = {}
for tau in taus:
    curve = curves[tau]
    target_points[tau] = sum(
        (coefficient * point for coefficient, point in zip(target_coefficients, basis_points_at(tau))),
        curve(0),
    )

# Retain the large common ordinary-fibre domain.  The omitted records have a
# ramified or dropped inverse fibre at one of 2,3,4 and can be revisited only
# if this complete common-domain meet-in-the-middle has no geometric winner.
fingerprints = {}
domains = {}
for record in abel["classification"]["records"]:
    trials = {
        int(trial["old_base_value"]): trial["trace"]
        for trial in record["ordinary_fibre_trials"]
    }
    if all(tau in trials for tau in taus):
        shell_index = int(record["shell_index"])
        fingerprints[shell_index] = tuple(
            curves[tau](trials[tau]) if trials[tau] is not None else curves[tau](0)
            for tau in taus
        )
        # The ten-fingerprint classifier attempts every value 2,...,11 for
        # every nonvertical shell section.  Later values are only retries
        # after failures, so restrict mixed-domain comparisons to this fixed,
        # uniformly audited universe.
        domains[shell_index] = frozenset(set(trials).intersection(range(2, 12)))
indices = sorted(fingerprints)
assert indices

pairs_by_fingerprint = defaultdict(list)
pair_count = 0
for left_position, left in enumerate(indices):
    for right in indices[left_position:]:
        key = tuple(
            point_key(fingerprints[left][slot] + fingerprints[right][slot])
            for slot in range(len(taus))
        )
        pairs_by_fingerprint[key].append((left, right))
        pair_count += 1

quadruples = set()
for left_position, left in enumerate(indices):
    for right in indices[left_position:]:
        complement = tuple(
            point_key(
                target_points[tau] - fingerprints[left][slot] - fingerprints[right][slot]
            )
            for slot, tau in enumerate(taus)
        )
        for other in pairs_by_fingerprint.get(complement, ()):
            quadruple = tuple(sorted((left, right) + other))
            if len(taus) == 2 and len(set.intersection(*(
                set(domains[index]) for index in quadruple
            ))) != 2:
                continue
            quadruples.add(quadruple)
quadruples = sorted(quadruples)

# Exact child function-field checks on the fingerprint hits.
child_A = polynomial(pointing["global_short_model"]["A_coefficients_low_to_high"])
child_B = polynomial(pointing["global_short_model"]["B_coefficients_low_to_high"])
child_curve = EllipticCurve(K, [0, 0, 0, K(child_A), K(child_B)])
shell_points = [
    child_curve(K(R(record["x_coefficients_low_to_high"])), K(R(record["y_coefficients_low_to_high"])))
    for record in shell["shell"]["records"]
]
a1, a2, a3 = [
    rational_function(pointing["pointed_generalized_weierstrass"][name])
    for name in ("a1", "a2", "a3")
]
b2 = a1**2 + 4 * a2
w0 = K(polynomial(pointing["quartic_square_at_t0"]["L0_coefficients_low_to_high"]))


def inverse_parent_base(point):
    if point.is_zero():
        return None
    X, Y = point.xy()
    x_general = K(X) / 9 - b2 / 12
    y_general = K(Y) / 27 - (a1 * x_general + a3) / 2
    if not y_general:
        return None
    return K(2 * w0 * (x_general + a2) / y_general)


def rational_degree(value):
    return None if value is None else max(value.numerator().degree(), value.denominator().degree())


def point_dot_zero(point):
    if point.is_zero():
        return -2
    X = K(point[0])
    numerator_degree = X.numerator().degree()
    denominator_degree = X.denominator().degree()
    infinity_excess = max(0, numerator_degree - denominator_degree - 4)
    if denominator_degree % 2 or infinity_excess % 2:
        raise ArithmeticError("nonintegral pole count")
    return denominator_degree // 2 + infinity_excess // 2


candidate_rows = []
seen_points = set()
for quadruple in quadruples:
    point = sum((shell_points[index] for index in quadruple), child_curve(0))
    if point.is_zero():
        continue
    X, Y = map(K, point.xy())
    key = (
        tuple(R(X.numerator()).list()), tuple(R(X.denominator()).list()),
        tuple(R(Y.numerator()).list()), tuple(R(Y.denominator()).list()),
    )
    if key in seen_points:
        continue
    seen_points.add(key)
    parent_base = inverse_parent_base(point)
    row = {
        "shell_indices": list(quadruple),
        "P_dot_O": int(point_dot_zero(point)),
        "inverse_parent_degree": int(rational_degree(parent_base)),
        "shape_Xnum_Xden_Ynum_Yden": [
            int(X.numerator().degree()), int(X.denominator().degree()),
            int(Y.numerator().degree()), int(Y.denominator().degree()),
        ],
        "x": {
            "numerator_coefficients_low_to_high": list(map(int, R(X.numerator()).list())),
            "denominator_coefficients_low_to_high": list(map(int, R(X.denominator()).list())),
        },
        "y": {
            "numerator_coefficients_low_to_high": list(map(int, R(Y.numerator()).list())),
            "denominator_coefficients_low_to_high": list(map(int, R(Y.denominator()).list())),
        },
    }
    row["passes_q207_degree_and_pole_gates"] = (
        row["P_dot_O"] == 10 and row["inverse_parent_degree"] == 16
    )
    candidate_rows.append(row)

winners = [row for row in candidate_rows if row["passes_q207_degree_and_pole_gates"]]
payload = {
    "schema": "elkies-k3.h92-q4o323-q207-shell-abel-mitm-modp.v1",
    "status": (
        "PASS_MODP_Q4O323_Q207_UNIQUE_SHELL_ABEL_WORD"
        if len(winners) == 1 else
        "PASS_MODP_Q4O323_Q207_SHELL_ABEL_WORDS_WITH_RESIDUAL_AMBIGUITY"
    ),
    "prime": int(prime),
    "ordinary_parent_fibres": list(taus),
    "target": {
        "equation_MW_basis": list(basis_labels),
        "coefficients": list(target_coefficients),
        "fibre_points": {str(tau): list(map(int, target_points[tau].xy())) for tau in taus},
        "required_P_dot_O": 10,
        "required_inverse_parent_degree": 16,
    },
    "search": {
        "common_domain_shell_sections": len(indices),
        "omitted_ramified_or_dropped_sections": len(shell["shell"]["records"]) - len(indices),
        "unordered_pairs": pair_count,
        "fingerprint_quadruples": len(quadruples),
        "distinct_child_points": len(candidate_rows),
        "degree_and_pole_winners": len(winners),
        "winners": winners,
        "all_candidate_points": candidate_rows,
    },
    "method": {
        "large_Groebner_required": False,
        "multisection_splitting_required": False,
        "elimination_required": False,
        "construction": "three-fibre Abel fingerprint pair/pair meet-in-the-middle",
        "runtime_seconds": time.monotonic() - started,
    },
    "proof_boundary": (
        "The winner set is exact over GF(p) on the displayed common ordinary-fibre domain. "
        "If unique, it is a modular q207 seed; characteristic-zero lifting and the q12 "
        "resolved Riemann--Roch/Jacobian certificate remain separate gates."
    ),
    "inputs": {
        "paths": [display_path(path) for path in INPUTS],
        "sha256": {display_path(path): sha256(path) for path in INPUTS},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4O323Q207MITM|prime={}|shell={}|pairs={}|quads={}|points={}|winners={}|runtime={:.3f}|"
    "status={}|output={}".format(
        prime, len(indices), pair_count, len(quadruples), len(candidate_rows), len(winners),
        payload["method"]["runtime_seconds"], payload["status"], OUTPUT,
    ),
    flush=True,
)
