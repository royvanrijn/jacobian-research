#!/usr/bin/env sage -python
"""Hensel-filter the polynomial-section pool and identify P1229 at p=103.

status: EXPERIMENT
claim: regular structured p-adic lift filter and coefficient-growth diagnostics
output: artifacts/local/elkies-k3/q24-2a5-p1229-hensel-p103.json

The zero-pole section has twelve coefficients (x0..x4,y0..y6).  Thirteen
coefficient identities define it, but the Jacobian has rank twelve.  This
script selects twelve independent rows modulo 103 and applies Newton--Hensel
directly to that square regular chart for all 130 modular candidates.  The
unused thirteenth equation filters out specialization-only points.  No
Groebner basis is used.
"""

import hashlib
import json
import time
import argparse
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, Qp, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
EXACT = LOCAL / "q24-a11-to-2a5-q8-resolved-rr-qq.json"
POOL = LOCAL / "q24-2a5-zero-pole-sections-p103.json"
MATCH = LOCAL / "q24-2a5-zero-pole-shell-match-p103.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--precision", type=int, default=160)
parser.add_argument(
    "--target-indices",
    help="comma-separated modular section indices; default probes the full pool",
)
parser.add_argument(
    "--checkpoint-digits",
    help="comma-separated reconstruction checkpoints below the working precision",
)
parser.add_argument(
    "--output",
    default="artifacts/local/elkies-k3/q24-2a5-p1229-hensel-p103.json",
)
args = parser.parse_args()
OUTPUT = Path(args.output)
if not OUTPUT.is_absolute():
    OUTPUT = ROOT / OUTPUT


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


started = time.monotonic()
exact = json.loads(EXACT.read_text())
pool = json.loads(POOL.read_text())
match = json.loads(MATCH.read_text())
assert exact["status"] == "PASS_EXACT_Q24_A11_Q8_2A5_RESOLVED_RR"
assert pool["status"] == "PASS_BOUNDED_MOD103_ZERO_POLE_SECTION_ENUMERATION"
assert match["status"] == "PASS_EXHAUSTIVE_MOD103_ZERO_POLE_SHELL_EMBEDDINGS_CANONICAL_MARKING"

p = ZZ(103)
p1229_candidate_indices = [
    int(index) for index in match["modular_pool"]["coarse_P1229_candidates"]
]
targeted = args.target_indices is not None
candidate_indices = (
    [int(index) for index in args.target_indices.split(",")]
    if targeted else list(range(len(pool["sections"])))
)
assert candidate_indices and len(set(candidate_indices)) == len(candidate_indices)
assert all(0 <= index < len(pool["sections"]) for index in candidate_indices)

A_QQ = [QQ(value) for value in exact["child"]["minimal_A_coefficients_low_to_high"]]
B_QQ = [QQ(value) for value in exact["child"]["minimal_B_coefficients_low_to_high"]]
assert all(value.denominator() % p for value in A_QQ + B_QQ)


def coefficient_height_bits(value):
    value = QQ(value)
    return max(abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())


surface_height = {
    "max_A_coefficient_bits": max(map(coefficient_height_bits, A_QQ)),
    "max_B_coefficient_bits": max(map(coefficient_height_bits, B_QQ)),
}

PRECISION = args.precision
assert PRECISION >= 40
CHECKPOINT_DIGITS = (
    tuple(int(value) for value in args.checkpoint_digits.split(","))
    if args.checkpoint_digits
    else ((20, 40, 80, 140) if PRECISION == 160 else (
        max(20, PRECISION // 8), PRECISION // 4, PRECISION // 2, PRECISION - 10
    ))
)
assert tuple(sorted(set(CHECKPOINT_DIGITS))) == CHECKPOINT_DIGITS
assert all(0 < digits < PRECISION for digits in CHECKPOINT_DIGITS)
K = Qp(p, prec=PRECISION, type="capped-rel")
RT = PolynomialRing(K, "T")
T = RT.gen()
A = RT([K(value) for value in A_QQ])
B = RT([K(value) for value in B_QQ])


def polynomials(values):
    # Sage vector slices are vectors, while the polynomial-ring constructor
    # expects a plain coefficient sequence here.
    return RT(list(values[:5])), RT(list(values[5:]))


def residual(values):
    X, Y = polynomials(values)
    R = Y**2 - X**3 - A * X - B
    return vector(K, [R[index] for index in range(13)])


def jacobian(values):
    X, Y = polynomials(values)
    x_derivative = -3 * X**2 - A
    y_derivative = 2 * Y
    rows = []
    for degree in range(13):
        rows.append([
            x_derivative[degree - shift] if 0 <= degree - shift <= x_derivative.degree() else K(0)
            for shift in range(5)
        ] + [
            y_derivative[degree - shift] if 0 <= degree - shift <= y_derivative.degree() else K(0)
            for shift in range(7)
        ])
    return matrix(K, rows)


F = GF(p)
RTF = PolynomialRing(F, "T")
TF = RTF.gen()
A_F = RTF([F(value.numerator()) / F(value.denominator()) for value in A_QQ])
B_F = RTF([F(value.numerator()) / F(value.denominator()) for value in B_QQ])


def minimum_valuation(values):
    nonzero = [value.valuation() for value in values if value]
    return min(nonzero) if nonzero else PRECISION


def residue_at_precision(value, digits):
    modulus = p**digits
    residue = ZZ(value.lift()) % modulus
    balanced = residue if residue <= modulus // 2 else residue - modulus
    return residue, balanced


def reconstruction_checkpoints(values):
    checkpoints = []
    previous_reconstructions = None
    for digits in CHECKPOINT_DIGITS:
        modulus = p**digits
        residues = []
        balanced_bits = []
        reconstructions = []
        for value in values:
            residue, balanced = residue_at_precision(value, digits)
            residues.append(str(residue))
            balanced_bits.append(abs(balanced).nbits())
            try:
                reconstructed = ZZ(residue).rational_reconstruction(modulus)
                reconstructions.append(str(reconstructed))
            except (ArithmeticError, ValueError):
                reconstructions.append(None)
        stable = 0
        if previous_reconstructions is not None:
            stable = sum(
                current is not None and current == previous
                for current, previous in zip(reconstructions, previous_reconstructions)
            )
        checkpoints.append({
            "p_adic_digits": digits,
            "modulus_bits": modulus.nbits(),
            "maximum_balanced_residue_bits": max(balanced_bits),
            "rational_reconstructions": reconstructions,
            "stable_from_previous_checkpoint": stable,
            "residues_0_to_p_power_minus_1": residues,
        })
        previous_reconstructions = reconstructions
    return checkpoints


def probe_candidate(section_index, include_checkpoints=False):
    section = pool["sections"][section_index]
    initial = (
        section["X_coefficients_low_to_high"]
        + section["Y_coefficients_low_to_high"]
    )
    assert len(initial) == 12
    X_F = RTF(initial[:5])
    Y_F = RTF(initial[5:])
    assert Y_F**2 == X_F**3 + A_F * X_F + B_F

    x_derivative_F = -3 * X_F**2 - A_F
    y_derivative_F = 2 * Y_F
    J_F = matrix(F, [[
        x_derivative_F[degree - shift]
        if 0 <= degree - shift <= x_derivative_F.degree() else F(0)
        for shift in range(5)
    ] + [
        y_derivative_F[degree - shift]
        if 0 <= degree - shift <= y_derivative_F.degree() else F(0)
        for shift in range(7)
    ] for degree in range(13)])
    rank = int(J_F.rank())
    if rank < 12:
        return {
            "section_index": section_index,
            "X_coefficients_low_to_high_mod103": initial[:5],
            "Y_coefficients_low_to_high_mod103": initial[5:],
            "mod103_jacobian_rank": rank,
            "full_thirteen_equation_hensel_lift": False,
            "regular_chart_available": False,
        }
    pivot_rows = list(map(int, J_F.transpose().pivots()))
    J_square_F = matrix(F, [J_F.row(index) for index in pivot_rows])
    determinant = int(J_square_F.det())
    assert determinant != 0

    values = vector(K, [K(value).add_bigoh(1) for value in initial])
    iterations = []
    known_precision = 1
    iteration = 0
    while known_precision < PRECISION:
        working_precision = min(2 * known_precision, PRECISION)
        values = vector(K, [
            K(value.lift()).add_bigoh(working_precision) for value in values
        ])
        all_residual = residual(values)
        square_residual = vector(K, [all_residual[index] for index in pivot_rows])
        square_before = int(minimum_valuation(square_residual))
        full_before = int(minimum_valuation(all_residual))
        full_jacobian = jacobian(values)
        J_square = matrix(K, [full_jacobian.row(index) for index in pivot_rows])
        # The determinant is a unit at the mod-p seed and remains a unit on
        # this Hensel ball; solve_right detects any contrary failure.  Avoid a
        # second cubic matrix operation at every high-precision doubling.
        correction = J_square.solve_right(-square_residual)
        values += correction
        iterations.append({
            "iteration": iteration,
            "working_precision_p_adic_digits": working_precision,
            "minimum_square_residual_valuation_before": square_before,
            "minimum_full_residual_valuation_before": full_before,
            "minimum_full_residual_valuation_after": int(
                minimum_valuation(residual(values))
            ),
            "minimum_correction_valuation": int(minimum_valuation(correction)),
        })
        known_precision = working_precision
        iteration += 1

    final_residual = residual(values)
    final_valuations = [
        int(value.valuation()) if value else PRECISION for value in final_residual
    ]
    full_lift = min(final_valuations) >= PRECISION - 8
    return {
        "section_index": section_index,
        "X_coefficients_low_to_high_mod103": initial[:5],
        "Y_coefficients_low_to_high_mod103": initial[5:],
        "mod103_jacobian_rank": rank,
        "selected_independent_equation_degrees": pivot_rows,
        "selected_jacobian_determinant_mod103": determinant,
        "iterations": iterations,
        "final_residual_valuations_by_degree": final_valuations,
        "final_minimum_full_residual_valuation": min(final_valuations),
        "full_thirteen_equation_hensel_lift": full_lift,
        "regular_chart_available": True,
        "checkpoints": (
            reconstruction_checkpoints(values)
            if full_lift and include_checkpoints else []
        ),
    }


candidate_probes = [
    probe_candidate(index, include_checkpoints=targeted)
    for index in candidate_indices
]
lifted = [probe for probe in candidate_probes if probe["full_thirteen_equation_hensel_lift"]]
p1229_lifted = [
    probe for probe in lifted if probe["section_index"] in p1229_candidate_indices
]
selected_probe = (
    probe_candidate(p1229_lifted[0]["section_index"], include_checkpoints=True)
    if len(p1229_lifted) == 1 else None
)
status = (
    (
        "PASS_TARGETED_REGULAR_MOD103_HENSEL_LIFTS"
        if len(lifted) == len(candidate_indices)
        else "INCONCLUSIVE_TARGETED_MOD103_HENSEL_LIFTS"
    )
    if targeted else (
        "PASS_REGULAR_MOD103_POOL_FILTER_UNIQUE_P1229_BRANCH"
        if len(p1229_lifted) == 1
        else "INCONCLUSIVE_MOD103_ZERO_POLE_HENSEL_BRANCH_COUNTS"
    )
)

payload = {
    "schema": "elkies-k3.q24-2a5-p1229-hensel-p103.v4",
    "status": status,
    "software": "SageMath 10.9 (conda-forge pinned repository environment)",
    "prime": int(p),
    "candidate_source": {
        "scope": "targeted" if targeted else "full_modular_pool",
        "probed_modular_polynomial_section_indices": candidate_indices,
        "shell_compatible_P1229_indices": p1229_candidate_indices,
        "previous_lexicographic_shell_choice": int(
            match["selected_P1229_mod103"]["section_index"]
        ),
    },
    "regular_system": {
        "variables": 12,
        "coefficient_equations": 13,
    },
    "hensel": {
        "working_precision_p_adic_digits": PRECISION,
        "rational_reconstruction_checkpoint_digits": list(CHECKPOINT_DIGITS),
        "candidate_probes": candidate_probes,
        "full_lift_indices": [probe["section_index"] for probe in lifted],
        "full_lift_count": len(lifted),
        "P1229_full_lift_indices": [
            probe["section_index"] for probe in p1229_lifted
        ],
        "unique_selected_probe": selected_probe,
    },
    "coefficient_growth_context": surface_height,
    "large_Groebner_required": False,
    "proof_boundary": (
        "The thirteenth coefficient identity is used as an obstruction after lifting "
        "regular twelve-equation charts for the requested modular candidates. "
        + (
            "This targeted run does not classify the rest of the modular pool. "
            if targeted else
            "Rank-deficient modular points are not classified, so the survivor count is "
            "not asserted to be the complete characteristic-zero zero-pole shell. "
        )
        + "A p-adic branch and rational-reconstruction diagnostics do not by themselves "
        "prove rationality or an exact QQ section identity."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (EXACT, POOL, MATCH)],
        "sha256": {
            str(path.relative_to(ROOT)): sha256(path) for path in (EXACT, POOL, MATCH)
        },
    },
    "elapsed_seconds": round(time.monotonic() - started, 6),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A5A5P1229HENSEL|"
    f"scope={'targeted' if targeted else 'pool'}|candidates={len(candidate_indices)}|full_lifts={len(lifted)}|"
    f"P1229_full_lifts={tuple(payload['hensel']['P1229_full_lift_indices'])}|"
    f"A_bits={surface_height['max_A_coefficient_bits']}|B_bits={surface_height['max_B_coefficient_bits']}|"
    f"status={payload['status']}|output={OUTPUT}",
    flush=True,
)
