#!/usr/bin/env sage -python
"""Lift the two new regular p=89 q12/o5867 replacement-word seeds to QQ.

The full-shell lattice classifier selects shells 206/class499 and 172/class511
for the replacement four-section word.  For each seed independently, Newton--Hensel lift the
13 coefficient equations in the twelve section coefficients, try rational
reconstruction at successively doubled precision, and accept only a literal
QQ Weierstrass identity with exact reduction to the supplied mod-89 seed.

No elimination or Groebner basis is used.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, Qp, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
Q8 = LOCAL / "q4o164-q8o376-smooth-rr-qq.json"
SHELL = LOCAL / "q12o5867-p0-shell-all-records-mod89.json"
CLASSIFIER = LOCAL / "q12o5867-p0-shell-lattice-classification-mod89.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--reconstruction-start", type=int, default=64)
parser.add_argument("--maximum-precision", type=int, default=4096)
parser.add_argument(
    "--class-index", action="append", type=int, default=[],
    help="lift the first regular uniquely classified seed for this physical class",
)
parser.add_argument(
    "--class-shell", action="append", default=[], metavar="CLASS:SHELL",
    help="lift an explicitly selected regular shell seed carrying the supplied modular class label",
)
parser.add_argument("--shell", type=Path, default=SHELL)
parser.add_argument("--classifier", type=Path, default=CLASSIFIER)
parser.add_argument(
    "--output",
    default="artifacts/local/elkies-k3/q12o5867-replacement-word-seeds-qq.json",
)
args = parser.parse_args()
SHELL = args.shell if args.shell.is_absolute() else ROOT / args.shell
CLASSIFIER = args.classifier if args.classifier.is_absolute() else ROOT / args.classifier
OUTPUT = Path(args.output)
if not OUTPUT.is_absolute():
    OUTPUT = ROOT / OUTPUT
reconstruction_start = int(args.reconstruction_start)
maximum_precision = int(args.maximum_precision)
assert 32 <= reconstruction_start <= maximum_precision


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficient_bits(value):
    value = QQ(value)
    return max(abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())


started = time.monotonic()
q8 = json.loads(Q8.read_text())
shell = json.loads(SHELL.read_text())
classifier = json.loads(CLASSIFIER.read_text())
assert q8["status"] == "PASS_EXACT_QQ_Q8O376_4A1_RR_JACOBIAN_AND_P1229_ZERO"
assert shell["prime"] == 89 and len(shell["all_records"]) > 0
assert classifier["prime"] == 89
prime = ZZ(shell["prime"])

RQ = PolynomialRing(QQ, "v")
vq = RQ.gen()
A_QQ = RQ([QQ(value) for value in q8["child"]["minimal_A_coefficients_low_to_high"]])
B_QQ = RQ([QQ(value) for value in q8["child"]["minimal_B_coefficients_low_to_high"]])


def reduce_qq(value, field):
    value = QQ(value)
    return field(value.numerator())/field(value.denominator())


F = GF(prime)
RF = PolynomialRing(F, "v")
A_F = RF([reduce_qq(value, F) for value in A_QQ])
B_F = RF([reduce_qq(value, F) for value in B_QQ])


def select_seed(shell_index, lattice_class_index):
    seed = dict(shell["all_records"][shell_index])
    classification = classifier["polynomial_shell"]["records"][shell_index]
    assert classification["shell_index"] == shell_index
    assert classification["ordinary_coefficient_jacobian_rank"] == 12
    assert lattice_class_index in classification["profile_compatible_lattice_class_indices"]
    seed["ordinary_fibre_trials"] = classification["ordinary_fibre_trials"]
    return seed


default_seed_specs = (
    {
        "key": "replacement_class499_shell206",
        "lattice_class_index": 499,
        "shell_index": 206,
        "seed": select_seed(206, 499),
        "selection_reason": "unique four-fibre Abel/profile class499 seed; ordinary coefficient Jacobian rank 12",
        "full_NS_name_certified": True,
    },
    {
        "key": "replacement_class511_shell172",
        "lattice_class_index": 511,
        "shell_index": 172,
        "seed": select_seed(172, 511),
        "selection_reason": "four-fibre Abel/profile class511 candidate selected from 511/546 ambiguity by exact Q2 smooth intersection; rank 12",
        "full_NS_name_certified": True,
    },
)

requested_class_shells = []
for value in args.class_shell:
    class_text, separator, shell_text = value.partition(":")
    if not separator:
        raise ValueError("--class-shell must have the form CLASS:SHELL")
    requested_class_shells.append((int(class_text), int(shell_text)))

if args.class_index or requested_class_shells:
    seed_specs = []
    requests = [(class_index, None) for class_index in args.class_index] + requested_class_shells
    for class_index, requested_shell_index in requests:
        matches = [
            row for row in classifier["polynomial_shell"]["records"]
            if row["ordinary_coefficient_jacobian_rank"] == 12
            and row["profile_compatible_lattice_class_indices"] == [class_index]
        ]
        if requested_shell_index is not None:
            matches = [row for row in matches if row["shell_index"] == requested_shell_index]
        if not matches:
            resolved_indices = {
                int(result["shell_index"])
                for result in classifier["polynomial_shell"].get(
                    "complete_pairwise_intersection_disambiguation", []
                )
                if result["resolved_uniquely"]
                and result["surviving_class_alternatives"] == [class_index]
            }
            matches = [
                classifier["polynomial_shell"]["records"][index]
                for index in sorted(resolved_indices)
                if classifier["polynomial_shell"]["records"][index][
                    "ordinary_coefficient_jacobian_rank"
                ] == 12
            ]
            if requested_shell_index is not None:
                matches = [
                    row for row in matches
                    if row["shell_index"] == requested_shell_index
                ]
        if not matches:
            raise ValueError(f"no regular uniquely classified seed for class {class_index}")
        selected = min(matches, key=lambda row: row["shell_index"])
        shell_index = int(selected["shell_index"])
        seed_specs.append({
            "key": f"physical_class{class_index}_shell{shell_index}",
            "lattice_class_index": class_index,
            "shell_index": shell_index,
            "seed": select_seed(shell_index, class_index),
            "selection_reason": (
                "first regular uniquely Abel/profile-classified modular seed in shell order; "
                "the supplied class index is a seed label, not a cross-prime NS-name proof"
            ),
            "full_NS_name_certified": False,
        })
    seed_specs = tuple(seed_specs)
else:
    seed_specs = default_seed_specs

K = Qp(prime, prec=maximum_precision, type="capped-rel")
RT = PolynomialRing(K, "v")
v = RT.gen()
A = RT([K(value) for value in A_QQ])
B = RT([K(value) for value in B_QQ])


def split(values, ring):
    return ring(list(values[:5])), ring(list(values[5:]))


def residual(values):
    X, Y = split(values, RT)
    equation = Y**2-X**3-A*X-B
    return vector(K, [equation[index] for index in range(13)])


def jacobian(values, ring, surface_A):
    X, Y = split(values, ring)
    dx = -3*X**2-surface_A
    dy = 2*Y
    zero = ring.base_ring().zero()
    return matrix(ring.base_ring(), [[
        dx[degree-shift] if 0 <= degree-shift <= dx.degree() else zero
        for shift in range(5)
    ]+[
        dy[degree-shift] if 0 <= degree-shift <= dy.degree() else zero
        for shift in range(7)
    ] for degree in range(13)])


def minimum_valuation(values, fallback):
    nonzero = [int(value.valuation()) for value in values if value]
    return min(nonzero) if nonzero else int(fallback)


def reconstruct_vector(values, usable_precision):
    modulus = prime**usable_precision
    answer = []
    for value in values:
        if not value:
            answer.append(QQ.zero())
            continue
        residue = ZZ(value.lift()) % modulus
        answer.append(QQ(residue.rational_reconstruction(modulus)))
    return answer


def lift_seed(spec):
    seed_record = spec["seed"]
    seed = vector(
        F,
        seed_record["x_coefficients_low_to_high"]
        + seed_record["y_coefficients_low_to_high"],
    )
    assert len(seed) == 12
    X_F, Y_F = split(seed, RF)
    assert Y_F**2 == X_F**3+A_F*X_F+B_F

    J_F = jacobian(seed, RF, A_F)
    rank = int(J_F.rank())
    assert rank == 12
    pivot_rows = list(map(int, J_F.transpose().pivots()))
    assert len(pivot_rows) == 12
    minor_determinant = int(matrix(F, [J_F.row(index) for index in pivot_rows]).det())
    assert minor_determinant != 0

    values = vector(K, [K(value).add_bigoh(1) for value in seed])
    known_precision = 1
    iterations = []
    reconstruction_attempts = []
    reconstructed = None
    X_exact = Y_exact = None

    while known_precision < maximum_precision:
        working_precision = min(2*known_precision, maximum_precision)
        values = vector(K, [
            K(value.lift()).add_bigoh(working_precision) for value in values
        ])
        full = residual(values)
        chosen = vector(K, [full[index] for index in pivot_rows])
        J = jacobian(values, RT, A)
        square = matrix(K, [J.row(index) for index in pivot_rows])
        correction = square.solve_right(-chosen)
        values += correction
        residual_after = residual(values)
        iterations.append({
            "working_precision_p_adic_digits": working_precision,
            "minimum_full_residual_valuation_after": minimum_valuation(
                residual_after, working_precision
            ),
            "minimum_correction_valuation": minimum_valuation(
                correction, working_precision
            ),
        })
        known_precision = working_precision

        if working_precision < reconstruction_start:
            continue
        usable_precision = working_precision-8
        attempt = {
            "working_precision_p_adic_digits": working_precision,
            "usable_precision_p_adic_digits": usable_precision,
            "rational_reconstruction_succeeded": False,
            "literal_weierstrass_identity": False,
        }
        try:
            candidate = reconstruct_vector(values, usable_precision)
            attempt["rational_reconstruction_succeeded"] = True
            attempt["maximum_candidate_rational_bits"] = max(
                map(coefficient_bits, candidate)
            )
            X_candidate, Y_candidate = split(candidate, RQ)
            exact_identity = Y_candidate**2 == X_candidate**3+A_QQ*X_candidate+B_QQ
            exact_reduction = [reduce_qq(value, F) for value in candidate] == list(seed)
            attempt["literal_weierstrass_identity"] = bool(exact_identity)
            attempt["exact_reduction_to_mod89_seed"] = bool(exact_reduction)
            if exact_identity and exact_reduction:
                reconstructed = candidate
                X_exact, Y_exact = X_candidate, Y_candidate
        except (ArithmeticError, ValueError, ZeroDivisionError) as error:
            attempt["failure"] = type(error).__name__
        reconstruction_attempts.append(attempt)
        if reconstructed is not None:
            break

    if reconstructed is None:
        raise ArithmeticError(
            "{} did not reconstruct by {} p-adic digits".format(
                spec["key"], maximum_precision
            )
        )
    assert X_exact.degree() <= 4 and Y_exact.degree() <= 6
    assert Y_exact**2 == X_exact**3+A_QQ*X_exact+B_QQ
    assert [reduce_qq(value, F) for value in reconstructed] == list(seed)

    return {
        "lattice_class_index": spec["lattice_class_index"],
        "full_NS_name_certified": spec["full_NS_name_certified"],
        "selection": {
            "shell_index": spec["shell_index"],
            "selection_reason": spec["selection_reason"],
            "inverse_parent_degree_mod89": int(seed_record["inverse_parent_degree"]),
            "equation_component_profile": seed_record["equation_component_profile"],
            "ordinary_fibre_trials": seed_record["ordinary_fibre_trials"],
        },
        "section": {
            "x_coefficients_low_to_high": [str(value) for value in X_exact.list()],
            "y_coefficients_low_to_high": [str(value) for value in Y_exact.list()],
            "degrees_x_y": [int(X_exact.degree()), int(Y_exact.degree())],
            "P_dot_O": 0,
            "maximum_rational_bits": max(map(coefficient_bits, reconstructed)),
            "exact_weierstrass_identity": True,
            "exact_reduction_to_mod89_seed": True,
            "mod89_x_coefficients_low_to_high": [int(value) for value in X_F.list()],
            "mod89_y_coefficients_low_to_high": [int(value) for value in Y_F.list()],
        },
        "hensel": {
            "prime": int(prime),
            "coefficient_equations": 13,
            "variables": 12,
            "mod89_jacobian_rank": rank,
            "selected_independent_equation_rows": pivot_rows,
            "selected_minor_determinant_mod89": minor_determinant,
            "successful_working_precision_p_adic_digits": known_precision,
            "iterations": iterations,
            "rational_reconstruction_attempts": reconstruction_attempts,
        },
    }


sections = {}
for spec in seed_specs:
    section_started = time.monotonic()
    record = lift_seed(spec)
    record["runtime_seconds"] = time.monotonic()-section_started
    sections[spec["key"]] = record
    print(
        "Q12O5867REPLACEMENTQQ|seed={}|rank={}|minor={}|precision={}|degrees={}|bits={}|"
        "status=PASS_EXACT_QQ_SECTION".format(
            spec["key"],
            record["hensel"]["mod89_jacobian_rank"],
            record["hensel"]["selected_minor_determinant_mod89"],
            record["hensel"]["successful_working_precision_p_adic_digits"],
            ",".join(map(str, record["section"]["degrees_x_y"])),
            record["section"]["maximum_rational_bits"],
        ),
        flush=True,
    )

payload = {
    "schema": (
        "elkies-k3.h92-q12o5867-selected-p0-sections-qq.v1"
        if args.class_index else
        "elkies-k3.h92-q12o5867-replacement-word-seeds-qq.v1"
    ),
    "status": (
        "PASS_EXACT_QQ_Q12O5867_SELECTED_P0_SECTIONS"
        if args.class_index else
        "PASS_EXACT_QQ_Q12O5867_REPLACEMENT_CLASS499_CLASS511_SECTIONS"
    ),
    "surface": {
        "model": "exact q8/o376 child",
        "equation": "y^2=x^3+A(v)*x+B(v)",
        "literal_QQ_coefficients_source": str(Q8.relative_to(ROOT)),
    },
    "sections": sections,
    "method": {
        "construction": "independent regular-branch Newton--Hensel lifts with adaptive rational reconstruction",
        "reconstruction_start_p_adic_digits": reconstruction_start,
        "maximum_precision_p_adic_digits": maximum_precision,
        "large_Groebner_required": False,
        "elimination_required": False,
        "runtime_seconds": time.monotonic()-started,
    },
    "proof_boundary": (
        "Both displayed degree-(4,6) polynomial sections satisfy the exact q8/o376 child "
        "Weierstrass equation over QQ and reduce literally to their selected mod-89 seeds. "
        "The class499 label is unique from Abel/profile/degree data; class511 is selected from "
        "the residual class511/class546 pair by the exact smooth intersection with Q2. "
        "The q12/o5867 resolved RR pencil and endpoint equation remain separate gates."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (Q8, SHELL, CLASSIFIER)],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in (Q8, SHELL, CLASSIFIER)},
    },
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q12O5867REPLACEMENTQQ|sections={}|status={}|output={}".format(
        len(sections), payload["status"], OUTPUT
    ),
    flush=True,
)
