#!/usr/bin/env sage -python
"""Hensel-lift a selected simple-pole section at p=103.

status: EXPERIMENT
claim: regular 24-variable projective lift and rational reconstruction audit
output: artifacts/local/elkies-k3/q24-2a5-p230-projective-hensel-p103.json

Write deg(X,Y,Z)=(8,12,2), normalize Z to be monic, and solve the 25
coefficient identities in 24 variables.  Optionally freeze a quadratic Z
candidate that reconstructed identically at the last two checkpoints of an
earlier run; this leaves a 22-variable exact-QQ specialization.  The modular
Jacobian has full column rank.  Newton--Hensel uses independent rows and checks
every coefficient identity at every iterate.  No Groebner basis is used.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, Qp, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
EXACT = LOCAL / "q24-a11-to-2a5-q8-resolved-rr-qq.json"
REFINEMENT = LOCAL / "q24-2a5-zero-pole-hensel-refinement-p103.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--section",
    choices=("p230", "q2-plus-q3", "p146"),
    default="p230",
)
parser.add_argument("--precision", type=int, default=2000)
parser.add_argument("--checkpoint-digits")
parser.add_argument(
    "--fixed-z-source",
    help=(
        "earlier output whose last two checkpoints have identical non-null "
        "rational reconstructions for both non-monic Z coefficients"
    ),
)
parser.add_argument(
    "--resume-source",
    help="earlier compatible output whose last checkpoint seeds this lift",
)
parser.add_argument(
    "--output",
    default="artifacts/local/elkies-k3/q24-2a5-p230-projective-hensel-p103.json",
)
args = parser.parse_args()
assert args.precision >= 40
output = Path(args.output)
if not output.is_absolute():
    output = ROOT / output
checkpoints_requested = (
    tuple(int(value) for value in args.checkpoint_digits.split(","))
    if args.checkpoint_digits else (
        max(20, args.precision // 8), args.precision // 4,
        args.precision // 2, args.precision - 10,
    )
)
assert tuple(sorted(set(checkpoints_requested))) == checkpoints_requested
assert all(0 < value < args.precision for value in checkpoints_requested)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


started = time.monotonic()
exact = json.loads(EXACT.read_text())
refinement = json.loads(REFINEMENT.read_text())
assert exact["status"] == "PASS_EXACT_Q24_A11_Q8_2A5_RESOLVED_RR"
assert refinement["status"] == "PASS_EXACT_SHELL_RELATION_REGULARLY_SELECTED_P230_MOD103_CONSTRUCTION"

p = ZZ(103)
precision = args.precision
A_QQ = [QQ(value) for value in exact["child"]["minimal_A_coefficients_low_to_high"]]
B_QQ = [QQ(value) for value in exact["child"]["minimal_B_coefficients_low_to_high"]]
if args.section == "p230":
    degrees_X_Y_Z = (8, 12, 2)
    seed_record = refinement["P230_mod103"]
elif args.section == "q2-plus-q3":
    degrees_X_Y_Z = (6, 9, 1)
    # This is the literal elliptic-curve sum of modular polynomial-section
    # indices 62 and 36 selected by the exact shell refinement.
    seed_record = {
        "construction": "section_62 + section_36 mod 103",
        "X_coefficients_low_to_high": [4, 99, 58, 28, 72, 29, 56],
        "Y_coefficients_low_to_high": [88, 43, 4, 38, 42, 71, 36, 32, 23, 9],
        "Z_coefficients_low_to_high": [41, 1],
    }
else:
    degrees_X_Y_Z = (6, 9, 1)
    # The exact shell relation for P146 is exact-index 2 + exact-index 33.
    # Every marking compatible with P1229 -> modular index 115 sends these
    # two classes to modular polynomial sections 62 and 107.  Their literal
    # elliptic-curve sum supplies this regular simple-pole seed.
    seed_record = {
        "construction": "section_62 + section_107 mod 103",
        "X_coefficients_low_to_high": [51, 81, 65, 33, 68, 52, 27],
        "Y_coefficients_low_to_high": [39, 13, 22, 37, 30, 94, 45, 102, 63, 22],
        "Z_coefficients_low_to_high": [23, 1],
    }
x_count = degrees_X_Y_Z[0] + 1
y_count = degrees_X_Y_Z[1] + 1
z_count = degrees_X_Y_Z[2]
y_start = x_count
z_start = x_count + y_count
initial_full = (
    seed_record["X_coefficients_low_to_high"]
    + seed_record["Y_coefficients_low_to_high"]
    + seed_record["Z_coefficients_low_to_high"][:z_count]
)
assert len(initial_full) == z_start + z_count
assert seed_record["Z_coefficients_low_to_high"][z_count] == 1

fixed_z_source = None if args.fixed_z_source is None else Path(args.fixed_z_source)
assert fixed_z_source is None or args.section == "p230"
if fixed_z_source is not None and not fixed_z_source.is_absolute():
    fixed_z_source = ROOT / fixed_z_source
fixed_Z_QQ = None
if fixed_z_source is not None:
    fixed_payload = json.loads(fixed_z_source.read_text())
    fixed_checkpoints = fixed_payload["hensel"]["checkpoints"]
    assert len(fixed_checkpoints) >= 2
    previous_reconstructions = fixed_checkpoints[-2]["rational_reconstructions"]
    current_reconstructions = fixed_checkpoints[-1]["rational_reconstructions"]
    assert len(previous_reconstructions) == len(current_reconstructions) == len(initial_full)
    fixed_z_strings = current_reconstructions[z_start:z_start + z_count]
    assert all(value is not None for value in fixed_z_strings)
    assert previous_reconstructions[z_start:z_start + z_count] == fixed_z_strings
    fixed_Z_QQ = [QQ(value) for value in fixed_z_strings] + [QQ(1)]

initial = initial_full if fixed_Z_QQ is None else initial_full[:z_start]

resume_source = None if args.resume_source is None else Path(args.resume_source)
if resume_source is not None and not resume_source.is_absolute():
    resume_source = ROOT / resume_source
resume_digits = None
resume_residues = None
if resume_source is not None:
    resume_payload = json.loads(resume_source.read_text())
    assert resume_payload.get("section", "p230") == args.section
    assert tuple(resume_payload["projective_chart"]["degrees_X_Y_Z"]) == degrees_X_Y_Z
    resume_checkpoint = resume_payload["hensel"]["checkpoints"][-1]
    resume_digits = int(resume_checkpoint["p_adic_digits"])
    assert 0 < resume_digits < precision
    resume_residues = [ZZ(value) for value in
                       resume_checkpoint["residues_0_to_p_power_minus_1"]]
    assert len(resume_residues) == len(initial)
    assert all(int(value % p) == seed for value, seed in zip(resume_residues, initial))

K = Qp(p, prec=precision, type="capped-rel")
RT = PolynomialRing(K, "T")
T = RT.gen()
A = RT([K(value) for value in A_QQ])
B = RT([K(value) for value in B_QQ])


def polynomials(values):
    Z = (
        RT(list(values[z_start:z_start + z_count]) + [K(1)])
        if fixed_Z_QQ is None
        else RT([K(value) for value in fixed_Z_QQ])
    )
    return (
        RT(list(values[:x_count])),
        RT(list(values[y_start:z_start])),
        Z,
    )


def residual(values):
    X, Y, Z = polynomials(values)
    identity = Y**2 - X**3 - A * X * Z**4 - B * Z**6
    equation_count = max(
        2 * degrees_X_Y_Z[1],
        3 * degrees_X_Y_Z[0],
        A.degree() + degrees_X_Y_Z[0] + 4 * degrees_X_Y_Z[2],
        B.degree() + 6 * degrees_X_Y_Z[2],
    ) + 1
    return vector(K, [identity[index] for index in range(equation_count)])


def jacobian(values):
    X, Y, Z = polynomials(values)
    columns = []
    x_derivative = -3 * X**2 - A * Z**4
    y_derivative = 2 * Y
    z_derivative = -4 * A * X * Z**3 - 6 * B * Z**5
    columns.extend(x_derivative * T**shift for shift in range(x_count))
    columns.extend(y_derivative * T**shift for shift in range(y_count))
    if fixed_Z_QQ is None:
        columns.extend(z_derivative * T**shift for shift in range(z_count))
    return matrix(
        K, len(residual(values)), len(initial),
        lambda row, column: columns[column][row],
    )


F = GF(p)
RTF = PolynomialRing(F, "T")
TF = RTF.gen()
A_F = RTF([F(value.numerator()) / F(value.denominator()) for value in A_QQ])
B_F = RTF([F(value.numerator()) / F(value.denominator()) for value in B_QQ])
X_F = RTF(initial_full[:x_count])
Y_F = RTF(initial_full[y_start:z_start])
Z_F = RTF(list(initial_full[z_start:z_start + z_count]) + [1])
if fixed_Z_QQ is not None:
    fixed_Z_F = RTF([
        F(value.numerator()) / F(value.denominator()) for value in fixed_Z_QQ
    ])
    assert fixed_Z_F == Z_F
assert Y_F**2 == X_F**3 + A_F * X_F * Z_F**4 + B_F * Z_F**6
columns_F = []
columns_F.extend(
    (-3 * X_F**2 - A_F * Z_F**4) * TF**shift for shift in range(x_count)
)
columns_F.extend(2 * Y_F * TF**shift for shift in range(y_count))
if fixed_Z_QQ is None:
    columns_F.extend(
        (-4 * A_F * X_F * Z_F**3 - 6 * B_F * Z_F**5) * TF**shift
        for shift in range(z_count)
    )
equation_count = len(residual(vector(K, [K(value) for value in initial])))
J_F = matrix(
    F, equation_count, len(initial),
    lambda row, column: columns_F[column][row],
)
assert J_F.rank() == len(initial), (args.section, J_F.rank(), len(initial))
pivot_rows = tuple(map(int, J_F.transpose().pivots()))
J_square_F = matrix(F, [J_F.row(index) for index in pivot_rows])
determinant = int(J_square_F.det())
if args.section == "p230" and fixed_Z_QQ is None:
    assert determinant == 71


def minimum_valuation(values):
    valuations = [value.valuation() for value in values if value]
    return min(valuations) if valuations else precision


values = vector(K, [
    K(value).add_bigoh(1 if resume_digits is None else resume_digits)
    for value in (initial if resume_residues is None else resume_residues)
])
iterations = []
known_precision = 1 if resume_digits is None else resume_digits
iteration = 0
while known_precision < precision:
    working_precision = min(2 * known_precision, precision)
    values = vector(K, [
        K(value.lift()).add_bigoh(working_precision) for value in values
    ])
    all_residual = residual(values)
    square_residual = vector(K, [all_residual[index] for index in pivot_rows])
    before = int(minimum_valuation(square_residual))
    full_jacobian = jacobian(values)
    square_jacobian = matrix(K, [full_jacobian.row(index) for index in pivot_rows])
    # Full-rank reduction mod p certifies a unit determinant throughout this
    # Hensel ball.  Do not recompute a redundant high-precision determinant.
    correction = square_jacobian.solve_right(-square_residual)
    values += correction
    iterations.append({
        "iteration": iteration,
        "working_precision_p_adic_digits": working_precision,
        "minimum_square_residual_valuation_before": before,
        "minimum_correction_valuation": int(minimum_valuation(correction)),
        "minimum_full_residual_valuation_after": int(minimum_valuation(residual(values))),
    })
    known_precision = working_precision
    iteration += 1
final_residual_valuations = [
    int(value.valuation()) if value else precision for value in residual(values)
]
assert min(final_residual_valuations) >= precision - 8

# Canonical 6 -> 2 resolved RR plane for P230.  Use AA constant/linear
# coefficients as free pivots and solve for AA_2..AA_4 and constant BB.
rr_reconstruction_values = vector(K, [])
if args.section == "p230":
    X_lift, Y_lift, Z_lift = polynomials(values)
    collision_modulus = Z_lift**2
    ambient_pairs = (
        [(T**degree, RT.zero()) for degree in range(5)]
        + [(RT.zero(), RT.one())]
    )
    collision_remainders = [
        RT((AA * X_lift - BB * Y_lift) % collision_modulus)
        for AA, BB in ambient_pairs
    ]
    collision_matrix = matrix(K, [
        [remainder[degree] for remainder in collision_remainders]
        for degree in range(4)
    ])
    unknown_columns = (2, 3, 4, 5)
    square_collision = matrix(K, [
        [collision_matrix[row, column] for column in unknown_columns]
        for row in range(4)
    ])
    assert square_collision.det().valuation() == 0
    rr_rows = []
    for free_column in (0, 1):
        solution = square_collision.solve_right(-collision_matrix.column(free_column))
        row = [K(0)] * 6
        row[free_column] = K(1)
        for column, value in zip(unknown_columns, solution):
            row[column] = value
        assert collision_matrix * vector(K, row) == 0
        rr_rows.append(row)
    rr_reconstruction_values = vector(K, rr_rows[0][2:] + rr_rows[1][2:])
    AA0_lift = RT(rr_rows[0][:5])
    BB0_lift = rr_rows[0][5]
    AA1_lift = RT(rr_rows[1][:5])
    BB1_lift = rr_rows[1][5]
    rr_determinant_lift = AA0_lift * BB1_lift - AA1_lift * BB0_lift
    rr_determinant_ratio, rr_determinant_remainder = rr_determinant_lift.quo_rem(
        collision_modulus
    )
    assert rr_determinant_ratio.degree() == 0 and rr_determinant_ratio[0]
    assert int(minimum_valuation(rr_determinant_remainder.list())) >= precision - 8


def reconstruct_at_digits(value, digits):
    modulus = p**digits
    residue = ZZ(value.lift()) % modulus
    try:
        reconstruction = ZZ(residue).rational_reconstruction(modulus)
    except (ArithmeticError, ValueError):
        reconstruction = None
    return residue, reconstruction


checkpoint_records = []
previous = None
previous_rr = None
exact_reconstruction = None
exact_rr_reconstruction = None
for digits in checkpoints_requested:
    modulus = p**digits
    residues = []
    reconstructions = []
    for value in values:
        residue, reconstruction = reconstruct_at_digits(value, digits)
        residues.append(str(residue))
        reconstructions.append(reconstruction)
    stable = 0 if previous is None else sum(
        current is not None and current == old
        for current, old in zip(reconstructions, previous)
    )
    exact_identity = False
    if all(value is not None for value in reconstructions):
        X_QQ = PolynomialRing(QQ, "T")(reconstructions[:x_count])
        Y_QQ = PolynomialRing(QQ, "T")(reconstructions[y_start:z_start])
        Z_QQ = PolynomialRing(QQ, "T")(
            list(reconstructions[z_start:z_start + z_count]) + [1]
            if fixed_Z_QQ is None else fixed_Z_QQ
        )
        ring_QQ = X_QQ.parent()
        exact_identity = (
            Y_QQ**2
            == X_QQ**3 + ring_QQ(A_QQ) * X_QQ * Z_QQ**4 + ring_QQ(B_QQ) * Z_QQ**6
        )
        if exact_identity:
            exact_reconstruction = reconstructions
    rr_residues = []
    rr_reconstructions = []
    for value in rr_reconstruction_values:
        residue, reconstruction = reconstruct_at_digits(value, digits)
        rr_residues.append(str(residue))
        rr_reconstructions.append(reconstruction)
    rr_stable = 0 if previous_rr is None else sum(
        current is not None and current == old
        for current, old in zip(rr_reconstructions, previous_rr)
    )
    rr_exact_square_determinant = False
    rr_recovered_Z = None
    if args.section == "p230" and all(
        value is not None for value in rr_reconstructions
    ):
        ring_QQ = PolynomialRing(QQ, "T")
        rr0 = [QQ(1), QQ(0)] + rr_reconstructions[:4]
        rr1 = [QQ(0), QQ(1)] + rr_reconstructions[4:]
        AA0_QQ, BB0_QQ = ring_QQ(rr0[:5]), rr0[5]
        AA1_QQ, BB1_QQ = ring_QQ(rr1[:5]), rr1[5]
        determinant_QQ = AA0_QQ * BB1_QQ - AA1_QQ * BB0_QQ
        if determinant_QQ.degree() == 4:
            monic = determinant_QQ / determinant_QQ.leading_coefficient()
            rr_exact_square_determinant = monic.is_square()
            if rr_exact_square_determinant:
                root = ring_QQ(monic.sqrt()).monic()
                if root.degree() == 2:
                    rr_recovered_Z = root
                    exact_rr_reconstruction = rr_reconstructions
    checkpoint_records.append({
        "p_adic_digits": digits,
        "modulus_bits": int(modulus.nbits()),
        "non_null_rational_reconstructions": sum(value is not None for value in reconstructions),
        "stable_from_previous_checkpoint": stable,
        "exact_Weierstrass_identity": exact_identity,
        "rational_reconstructions": [
            None if value is None else str(value) for value in reconstructions
        ],
        "residues_0_to_p_power_minus_1": residues,
        "resolved_RR": {
            "non_null_rational_reconstructions": sum(
                value is not None for value in rr_reconstructions
            ),
            "stable_from_previous_checkpoint": rr_stable,
            "determinant_is_scalar_times_square_quadratic": rr_exact_square_determinant,
            "recovered_monic_Z_coefficients_low_to_high": (
                None if rr_recovered_Z is None
                else [str(value) for value in rr_recovered_Z.list()]
            ),
            "rational_reconstructions_AA2_AA3_AA4_BB_for_two_rows": [
                None if value is None else str(value) for value in rr_reconstructions
            ],
            "residues_0_to_p_power_minus_1": rr_residues,
        },
    })
    previous = reconstructions
    previous_rr = rr_reconstructions

section_label = {
    "p230": "P230",
    "q2-plus-q3": "Q2_PLUS_Q3",
    "p146": "P146",
}[args.section]
status = (
    f"PASS_EXACT_QQ_{section_label}_PROJECTIVE_SECTION"
    if exact_reconstruction is not None
    else f"PASS_REGULAR_P103_{section_label}_PROJECTIVE_HENSEL_DIAGNOSTIC"
)
payload = {
    "schema": "elkies-k3.q24-2a5-projective-hensel-p103.v2",
    "status": status,
    "software": "SageMath 10.9 (conda-forge pinned repository environment)",
    "prime": int(p),
    "section": args.section,
    "projective_chart": {
        "degrees_X_Y_Z": list(degrees_X_Y_Z),
        "normalization": "Z monic",
        "variables": len(initial),
        "coefficient_equations": equation_count,
        "mod103_jacobian_rank": int(J_F.rank()),
        "selected_independent_equation_degrees": list(pivot_rows),
        "selected_jacobian_determinant_mod103": determinant,
        "seed_X_Y_Z": seed_record,
        "fixed_Z_QQ_candidate": (
            None if fixed_Z_QQ is None else [str(value) for value in fixed_Z_QQ]
        ),
    },
    "resolved_RR_chart": None if args.section != "p230" else {
        "ambient_dimension": 6,
        "AA_degree_bound": 4,
        "BB_degree_bound": 0,
        "collision_modulus": "Z^2",
        "collision_rows": 4,
        "kernel_dimension": 2,
        "canonical_free_columns": ["AA_0", "AA_1"],
        "determinant_identity_over_Qp": "AA0*BB1-AA1*BB0 = unit*Z^2",
        "exact_rational_reconstruction_found": exact_rr_reconstruction is not None,
    },
    "hensel": {
        "working_precision_p_adic_digits": precision,
        "iterations": iterations,
        "final_residual_valuations_by_degree": final_residual_valuations,
        "checkpoints": checkpoint_records,
    },
    "exact_QQ_section": (
        None if exact_reconstruction is None else {
            "X_coefficients_low_to_high": [
                str(value) for value in exact_reconstruction[:x_count]
            ],
            "Y_coefficients_low_to_high": [
                str(value) for value in exact_reconstruction[y_start:z_start]
            ],
            "Z_coefficients_low_to_high": [
                *(
                    [
                        str(value) for value in
                        exact_reconstruction[z_start:z_start + z_count]
                    ]
                    if fixed_Z_QQ is None
                    else [str(value) for value in fixed_Z_QQ[:z_count]]
                ),
                "1",
            ],
            "exact_Weierstrass_identity": True,
        }
    ),
    "large_Groebner_required": False,
    "proof_boundary": (
        "Full residual valuation proves a regular p-adic branch in the normalized "
        "projective chart. Only a displayed exact QQ reconstruction with literal "
        "Weierstrass substitution promotes this experiment to an exact section."
    ),
    "inputs": {
        "paths": [
            str(path.relative_to(ROOT)) for path in (EXACT, REFINEMENT)
        ] + ([] if fixed_z_source is None else [str(fixed_z_source)])
        + ([] if resume_source is None else [str(resume_source)]),
        "sha256": {
            str(path.relative_to(ROOT)): sha256(path) for path in (EXACT, REFINEMENT)
        } | (
            {} if fixed_z_source is None else {str(fixed_z_source): sha256(fixed_z_source)}
        ) | (
            {} if resume_source is None else {str(resume_source): sha256(resume_source)}
        ),
    },
    "elapsed_seconds": round(time.monotonic() - started, 6),
}
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A5A5P230HENSEL|rank={}|det={}|precision={}|iterations={}|exact={}|"
    "status={}|output={}".format(
        J_F.rank(), determinant, precision, len(iterations),
        exact_reconstruction is not None, status, output
    ),
    flush=True,
)
