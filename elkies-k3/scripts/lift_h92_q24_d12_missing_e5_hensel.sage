#!/usr/bin/env sage -python
"""Lift a normalized low-pole section by structured coefficient Hensel.

For D12 E5 there are 36 normalized section coefficients; for the A11 bridge
there are 60.  In either case the coefficients of

    Y^2 - X^3 - A*X*Z^4 - B*Z^6

give an overdetermined isolated system with one extra equation.  Starting
from the certified section modulo 100003, select a full-rank square subsystem
once and apply
quadratically convergent multivariable Hensel/Newton lifting.  The unused
equation is checked at every precision.  This avoids both a Groebner basis
and the tens of thousands of independent primes that ordinary CRT would need
for the very large rational normalization of this D12 model.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, Zmod, lcm, matrix


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--target",
    choices=(
        "d12-e5",
        "a11-bridge",
        "a11-bridge-marked",
        "a11-q8-residual",
        "a11-q8-difference",
        "a11-q8-zero",
        "a11-q8-trace-s5",
        "a11-q8-trace-s7",
        "a11-q8-trace-s17",
    ),
    default="d12-e5",
)
parser.add_argument("--max-exponent", type=int, default=131072)
parser.add_argument("--reconstruct-min-exponent", type=int, default=1024)
parser.add_argument("--seed", type=Path)
parser.add_argument("--checkpoint", type=Path)
parser.add_argument("--output", type=Path)
args = parser.parse_args()
if args.max_exponent < 1:
    raise SystemExit("--max-exponent must be positive")

CONFIG = {
    "d12-e5": {
        "seed": LOCAL / "q24-d12-missing-e5-section-mod100003.json",
        "model": LOCAL / "q24-d13-to-d12-component-valuation-qq.json",
        "checkpoint": LOCAL / "q24-d12-missing-e5-hensel-p100003.json",
        "output": LOCAL / "q24-d12-missing-e5-section-qq.json",
        "seed_status": "PASS_Q24_D12_MISSING_E5_SECTION_RECONSTRUCTION_MODP",
        "checkpoint_schema": "elkies-k3.h3-q24-d12-missing-e5-hensel.v1",
        "checkpoint_status": "PARTIAL_Q24_D12_MISSING_E5_HENSEL",
        "output_schema": "elkies-k3.h3-q24-d12-missing-e5-section-qq.v1",
        "output_status": "PASS_EXACT_Q24_D12_MISSING_E5_SECTION_QQ",
        "variable": "V",
        "label": "Q24E5",
        "section_metadata": {
            "P_dot_O": 4,
            "D12_MW_Abel_Jacobi": [0, 0, 0, 0, 1],
        },
        "proof_boundary": (
            "The displayed rational polynomials satisfy the exact characteristic-zero D12 Weierstrass identity. "
            "Their MW marking is inherited from the exact shell word and pinned-good-prime trace orientation."
        ),
    },
    "a11-bridge": {
        "seed": LOCAL / "q24-a11-bridge-m-section-mod100003.json",
        "model": LOCAL / "q24-d12-to-a11-orbit42-resolved-rr-qq.json",
        "checkpoint": LOCAL / "q24-a11-bridge-m-hensel-p100003.json",
        "output": LOCAL / "q24-a11-bridge-m-section-qq.json",
        "seed_status": "PASS_Q24_A11_BRIDGE_M_SECTION_RECONSTRUCTION_MODP",
        "checkpoint_schema": "elkies-k3.h3-q24-a11-bridge-m-hensel.v1",
        "checkpoint_status": "PARTIAL_Q24_A11_BRIDGE_M_HENSEL",
        "output_schema": "elkies-k3.h3-q24-a11-bridge-m-section-qq.v1",
        "output_status": "PASS_EXACT_Q24_A11_BRIDGE_M_SECTION_QQ",
        "variable": "T",
        "label": "A11BRIDGEM",
        "section_metadata": {
            "equation_P_dot_O": 8,
            "pinned_lattice_P_dot_O": 5,
            "equation_MW_Abel_Jacobi": [0, 2, -1, 1, 0, 1],
            "pinned_lattice_MW_Abel_Jacobi": [1, 0, 0, 0, 0, 1],
            "equation_zero_translation_MW": [1, -2, 1, -1, 0, 0],
            "physical_orientation": "C10",
        },
        "proof_boundary": (
            "The displayed rational polynomials satisfy the exact characteristic-zero A11 Weierstrass identity. "
            "The pinned lattice marking is inherited from the exact degree-13 bridge word, the C10 zero translation, "
            "and the unique pinned-good-prime orientation."
        ),
    },
    "a11-bridge-marked": {
        "seed": LOCAL / "q24-a11-bridge-m-section-marked-mod100003.json",
        "model": LOCAL / "q24-d12-to-a11-orbit42-resolved-rr-qq.json",
        "checkpoint": LOCAL / "q24-a11-bridge-m-marked-hensel-p100003.json",
        "output": LOCAL / "q24-a11-bridge-m-section-marked-qq.json",
        "seed_status": "PASS_Q24_A11_BRIDGE_M_SECTION_MARKED_RECONSTRUCTION_MODP",
        "checkpoint_schema": "elkies-k3.h3-q24-a11-bridge-m-marked-hensel.v1",
        "checkpoint_status": "PARTIAL_Q24_A11_BRIDGE_M_MARKED_HENSEL",
        "output_schema": "elkies-k3.h3-q24-a11-bridge-m-section-marked-qq.v1",
        "output_status": "PASS_EXACT_Q24_A11_BRIDGE_M_SECTION_MARKED_QQ",
        "variable": "T",
        "label": "A11BRIDGEMMARKED",
        "section_metadata": {
            "equation_P_dot_O": 8,
            "equation_I12_component_depth_up_to_negation": 1,
            "pinned_lattice_P_dot_O": 5,
            "equation_MW_Abel_Jacobi": [0, 2, -1, 1, 0, 1],
            "pinned_lattice_MW_Abel_Jacobi": [1, 0, 0, 0, 0, 1],
            "equation_zero_translation_MW": [1, -2, 1, -1, 0, 0],
            "physical_orientation": "C10",
        },
        "proof_boundary": (
            "The displayed rational polynomials satisfy the exact characteristic-zero A11 Weierstrass identity. "
            "The pinned lattice marking is inherited from the corrected modular trace word and the independent "
            "O_pinned-M and M-S7 profile gates at the pinned good prime."
        ),
    },
    "a11-q8-residual": {
        "seed": LOCAL / "q24-a11-q8-horizontal-points-mod100003.json",
        "seed_section_path": ("residual_P12_minus_M", "section"),
        "model": LOCAL / "q24-d12-to-a11-orbit42-resolved-rr-qq.json",
        "checkpoint": LOCAL / "q24-a11-q8-residual-hensel-p100003.json",
        "output": LOCAL / "q24-a11-q8-residual-section-qq.json",
        "seed_status": "PASS_Q24_A11_Q8_HORIZONTAL_POINTS_RECONSTRUCTION_MODP",
        "checkpoint_schema": "elkies-k3.h3-q24-a11-q8-residual-hensel.v1",
        "checkpoint_status": "PARTIAL_Q24_A11_Q8_RESIDUAL_HENSEL",
        "output_schema": "elkies-k3.h3-q24-a11-q8-residual-section-qq.v1",
        "output_status": "PASS_EXACT_Q24_A11_Q8_RESIDUAL_SECTION_QQ",
        "variable": "T",
        "label": "A11Q8RESIDUAL",
        "section_metadata": {
            "equation_P_dot_O": 0,
            "equation_I12_component_depth_up_to_negation": 3,
            "pinned_A11_MW_Abel_Jacobi": [-1, 0, -1, 0, 0, 0],
            "relation": "P12-M",
        },
        "proof_boundary": (
            "The displayed polynomial section satisfies the exact characteristic-zero A11 Weierstrass identity. "
            "Its marking is inherited from the exact marked M section, the corrected trace word, and the pinned-good-prime profile gates."
        ),
    },
    "a11-q8-difference": {
        "seed": LOCAL / "q24-a11-q8-horizontal-points-mod100003.json",
        "seed_section_path": ("marked_difference_P12_minus_Opinned", "section"),
        "model": LOCAL / "q24-d12-to-a11-orbit42-resolved-rr-qq.json",
        "checkpoint": LOCAL / "q24-a11-q8-difference-hensel-p100003.json",
        "output": LOCAL / "q24-a11-q8-difference-section-qq.json",
        "seed_status": "PASS_Q24_A11_Q8_HORIZONTAL_POINTS_RECONSTRUCTION_MODP",
        "checkpoint_schema": "elkies-k3.h3-q24-a11-q8-difference-hensel.v1",
        "checkpoint_status": "PARTIAL_Q24_A11_Q8_DIFFERENCE_HENSEL",
        "output_schema": "elkies-k3.h3-q24-a11-q8-difference-section-qq.v1",
        "output_status": "PASS_EXACT_Q24_A11_Q8_DIFFERENCE_SECTION_QQ",
        "variable": "T",
        "label": "A11Q8DIFFERENCE",
        "section_metadata": {
            "equation_P_dot_O": 6,
            "equation_I12_component_depth_up_to_negation": 6,
            "relation": "P12-O_pinned",
        },
        "proof_boundary": (
            "The displayed section satisfies the exact characteristic-zero A11 Weierstrass identity. "
            "Its P12-O_pinned marking is inherited from the exact pinned zero, the q8 trace word, "
            "and the unique central-I12 pinned-good-prime profile gate."
        ),
    },
    "a11-q8-zero": {
        "seed": LOCAL / "q24-a11-q8-horizontal-points-mod100003.json",
        "seed_section_path": ("pinned_zero", "section"),
        "model": LOCAL / "q24-d12-to-a11-orbit42-resolved-rr-qq.json",
        "checkpoint": LOCAL / "q24-a11-q8-zero-hensel-p100003.json",
        "output": LOCAL / "q24-a11-q8-pinned-zero-section-qq.json",
        "seed_status": "PASS_Q24_A11_Q8_HORIZONTAL_POINTS_RECONSTRUCTION_MODP",
        "checkpoint_schema": "elkies-k3.h3-q24-a11-q8-zero-hensel.v1",
        "checkpoint_status": "PARTIAL_Q24_A11_Q8_ZERO_HENSEL",
        "output_schema": "elkies-k3.h3-q24-a11-q8-pinned-zero-section-qq.v1",
        "output_status": "PASS_EXACT_Q24_A11_Q8_PINNED_ZERO_SECTION_QQ",
        "variable": "T",
        "label": "A11Q8ZERO",
        "section_metadata": {
            "equation_P_dot_O": 3,
            "equation_I12_component_depth_up_to_negation": 2,
            "equation_MW_Abel_Jacobi": [-1, 2, -1, 1, 0, 0],
            "relation": "O_pinned against C10 equation zero",
        },
        "proof_boundary": (
            "The displayed rational section satisfies the exact characteristic-zero A11 Weierstrass identity. "
            "Its pinned-zero marking is inherited from the exact zero-translation word and the pinned-good-prime Qminus/S7 pairing."
        ),
    },
    **{
        f"a11-q8-trace-{name.lower()}": {
            "seed": LOCAL / "q24-a11-pinned-zero-section-mod100003.json",
            "seed_section_path": ("selected_trace_sections", name),
            "model": LOCAL / "q24-d12-to-a11-orbit42-resolved-rr-qq.json",
            "checkpoint": LOCAL / f"q24-a11-q8-trace-{name.lower()}-hensel-p100003.json",
            "output": LOCAL / f"q24-a11-q8-trace-{name.lower()}-section-qq.json",
            "seed_status": "PASS_Q24_A11_PINNED_ZERO_SECTION_RECONSTRUCTION_MODP",
            "checkpoint_schema": f"elkies-k3.h3-q24-a11-q8-trace-{name.lower()}-hensel.v1",
            "checkpoint_status": f"PARTIAL_Q24_A11_Q8_TRACE_{name.upper()}_HENSEL",
            "output_schema": f"elkies-k3.h3-q24-a11-q8-trace-{name.lower()}-section-qq.v1",
            "output_status": f"PASS_EXACT_Q24_A11_Q8_TRACE_{name.upper()}_SECTION_QQ",
            "variable": "T",
            "label": f"A11Q8TRACE{name.upper()}",
            "section_metadata": {
                "trace_name": name,
                "equation_profile_from_seed": "retained in the hashed modular seed",
            },
            "proof_boundary": (
                f"The displayed {name} trace section satisfies the exact characteristic-zero A11 Weierstrass identity. "
                "Its signed MW marking is inherited from the pinned-good-prime trace-profile and Qminus pairing gates."
            ),
        }
        for name in ("S5", "S7", "S17")
    },
}[args.target]

SEED = (args.seed or CONFIG["seed"]).resolve()
MODEL_PATH = CONFIG["model"].resolve()
checkpoint = (args.checkpoint or CONFIG["checkpoint"]).resolve()
output = (args.output or CONFIG["output"]).resolve()
for path in (SEED, MODEL_PATH):
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

seed = json.loads(SEED.read_text())
model_artifact = json.loads(MODEL_PATH.read_text())
assert seed["status"] == CONFIG["seed_status"]
p = ZZ(seed["prime"])
A_QQ = [QQ(value) for value in model_artifact["child"]["minimal_A_coefficients_low_to_high"]]
B_QQ = [QQ(value) for value in model_artifact["child"]["minimal_B_coefficients_low_to_high"]]
A_denominator = lcm([ZZ(value.denominator()) for value in A_QQ])
B_denominator = lcm([ZZ(value.denominator()) for value in B_QQ])
A_numerators = [ZZ(value * A_denominator) for value in A_QQ]
B_numerators = [ZZ(value * B_denominator) for value in B_QQ]
seed_section = seed
for key in CONFIG.get("seed_section_path", ("section",)):
    seed_section = seed_section[key]
degrees = [int(value) for value in seed_section["degrees_X_Y_Z"]]
assert degrees == [
    len(seed_section["X_coefficients_low_to_high"]) - 1,
    len(seed_section["Y_coefficients_low_to_high"]) - 1,
    len(seed_section["Z_coefficients_low_to_high"]) - 1,
]
assert seed_section["Z_coefficients_low_to_high"][-1] == 1

NX, NY, NZ = degrees[0] + 1, degrees[1] + 1, degrees[2]
unknown_count = NX + NY + NZ
equation_count = max(
    2 * degrees[1],
    3 * degrees[0],
    len(A_QQ) - 1 + degrees[0] + 4 * degrees[2],
    len(B_QQ) - 1 + 6 * degrees[2],
) + 1
assert equation_count == unknown_count + 1
seed_values = [
    ZZ(value)
    for key in (
        "X_coefficients_low_to_high",
        "Y_coefficients_low_to_high",
    )
    for value in seed_section[key]
] + [ZZ(value) for value in seed_section["Z_coefficients_low_to_high"][:NZ]]
assert len(seed_values) == unknown_count
input_hashes = {
    str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
    for path in (SEED, MODEL_PATH)
}


model_reduction_cache = {}


def reduced_model(modulus):
    """Reduce A and B with one common-denominator inverse apiece."""

    modulus = ZZ(modulus)
    if modulus in model_reduction_cache:
        return model_reduction_cache[modulus]
    if A_denominator.gcd(modulus) != 1 or B_denominator.gcd(modulus) != 1:
        raise ZeroDivisionError(f"model denominator is not a unit modulo {modulus}")
    inverse_A = (A_denominator % modulus).inverse_mod(modulus)
    inverse_B = (B_denominator % modulus).inverse_mod(modulus)
    answer = (
        [(value % modulus) * inverse_A % modulus for value in A_numerators],
        [(value % modulus) * inverse_B % modulus for value in B_numerators],
    )
    # Only the current and target precision are useful.  Bounding this cache
    # avoids retaining every large p-power coefficient vector in a long run.
    model_reduction_cache[modulus] = answer
    if len(model_reduction_cache) > 2:
        del model_reduction_cache[next(iter(model_reduction_cache))]
    return answer


def mod_ring(modulus):
    """Construct Z/modulus without asking Sage to primality-test p^e."""

    modulus = ZZ(modulus)
    S = Zmod(modulus, is_field=(modulus == p))
    if modulus != p:
        # PolynomialRing otherwise calls is_field(), whose generic
        # implementation attempts a pointless primality proof for the known
        # composite integer p^e.  Seed the cached mathematical answer.
        S.is_field.set_cache(False)
    return S


def polynomials(values, modulus):
    S = mod_ring(modulus)
    R = PolynomialRing(S, "V")
    V = R.gen()
    X = R([S(value) for value in values[:NX]])
    Y = R([S(value) for value in values[NX:NX + NY]])
    Z = R([S(value) for value in values[NX + NY:]] + [S.one()])
    reduced_A, reduced_B = reduced_model(modulus)
    A = R([S(value) for value in reduced_A])
    B = R([S(value) for value in reduced_B])
    return R, V, X, Y, Z, A, B


def coefficient_vector(poly, length=equation_count):
    return [poly[index] for index in range(length)]


def residual(values, modulus):
    unused_R, unused_V, X, Y, Z, A, B = polynomials(values, modulus)
    return coefficient_vector(Y**2 - X**3 - A * X * Z**4 - B * Z**6)


def jacobian(values, modulus):
    unused_R, V, X, Y, Z, A, B = polynomials(values, modulus)
    columns = []
    for index in range(NX):
        columns.append(coefficient_vector((-3 * X**2 - A * Z**4) * V**index))
    for index in range(NY):
        columns.append(coefficient_vector(2 * Y * V**index))
    for index in range(NZ):
        columns.append(
            coefficient_vector((-4 * A * X * Z**3 - 6 * B * Z**5) * V**index)
        )
    return matrix(mod_ring(modulus), equation_count, unknown_count, lambda i, j: columns[j][i])


Jp = jacobian(seed_values, p)
if Jp.rank() != unknown_count:
    raise ArithmeticError(f"seed Jacobian rank {Jp.rank()}, expected {unknown_count}")
pivot_rows = list(Jp.transpose().pivots())
assert len(pivot_rows) == unknown_count
assert matrix(GF(p), Jp[pivot_rows, :]).det() != 0


def structured_remainder_matrix(values, modulus):
    """Linearized equations after eliminating delta-Y by division by 2Y."""

    unused_R, V, X, Y, Z, A, B = polynomials(values, modulus)
    divisor = 2 * Y
    if not divisor.leading_coefficient().is_unit():
        raise ArithmeticError("2Y has nonunit leading coefficient")
    x_factor = -3 * X**2 - A * Z**4
    z_factor = -4 * A * X * Z**3 - 6 * B * Z**5
    columns = []
    for factor, count in ((x_factor, NX), (z_factor, NZ)):
        for index in range(count):
            unused_quotient, remainder = (factor * V**index).quo_rem(divisor)
            columns.append([remainder[row] for row in range(NY - 1)])
    return matrix(
        mod_ring(modulus),
        NY - 1,
        NX + NZ,
        lambda row, column: columns[column][row],
    )


structured_p = structured_remainder_matrix(seed_values, p)
if structured_p.rank() != NX + NZ:
    raise ArithmeticError(
        f"structured seed Jacobian rank {structured_p.rank()}, expected {NX + NZ}"
    )
structured_pivot_rows = list(structured_p.transpose().pivots())
assert len(structured_pivot_rows) == NX + NZ


def structured_correction(values, modulus, rhs_coefficients):
    """Solve the Newton correction using a small X/Z remainder system."""

    R, V, X, Y, Z, A, B = polynomials(values, modulus)
    divisor = 2 * Y
    x_factor = -3 * X**2 - A * Z**4
    z_factor = -4 * A * X * Z**3 - 6 * B * Z**5
    rhs = R(rhs_coefficients)
    unused_quotient, rhs_remainder = rhs.quo_rem(divisor)
    small = structured_remainder_matrix(values, modulus)
    square = small[structured_pivot_rows, :]
    right = matrix(
        mod_ring(modulus),
        NX + NZ,
        1,
        [rhs_remainder[index] for index in structured_pivot_rows],
    )
    xz = square.solve_right(right).column(0)
    delta_x = R(list(xz[:NX]))
    delta_z = R(list(xz[NX:]))
    delta_y, remainder = (
        rhs - x_factor * delta_x - z_factor * delta_z
    ).quo_rem(divisor)
    if remainder:
        raise ArithmeticError("structured Newton correction failed unused remainder equation")
    if delta_y.degree() >= NY:
        raise ArithmeticError("structured Newton delta-Y exceeds its degree bound")
    answer = list(xz[:NX]) + [delta_y[index] for index in range(NY)] + list(xz[NX:])
    if len(answer) != unknown_count:
        raise AssertionError("structured correction length mismatch")
    return answer

if checkpoint.exists():
    saved = json.loads(checkpoint.read_text())
    if saved.get("schema") != CONFIG["checkpoint_schema"]:
        raise ArithmeticError("checkpoint schema mismatch")
    if saved.get("input_sha256") != input_hashes or ZZ(saved.get("prime")) != p:
        raise ArithmeticError("checkpoint prerequisites changed")
    exponent = int(saved["exponent"])
    values = [ZZ(value) for value in saved["coefficient_residues"]]
    if saved["pivot_equation_rows"] != pivot_rows:
        raise ArithmeticError("checkpoint pivot equations changed")
    if (
        "structured_pivot_remainder_rows" in saved
        and saved["structured_pivot_remainder_rows"] != structured_pivot_rows
    ):
        raise ArithmeticError("checkpoint structured pivot equations changed")
else:
    exponent = 1
    values = [value % p for value in seed_values]


def write_checkpoint():
    modulus = p**exponent
    payload = {
        "schema": CONFIG["checkpoint_schema"],
        "status": CONFIG["checkpoint_status"],
        "target": args.target,
        "prime": int(p),
        "exponent": exponent,
        "modulus_bits": int(modulus.nbits()),
        "pivot_equation_rows": pivot_rows,
        "structured_pivot_remainder_rows": structured_pivot_rows,
        "structured_unknown_count": NX + NZ,
        "coefficient_order": (
            f"X[0..{NX - 1}],Y[0..{NY - 1}],Z[0..{NZ - 1}] with Z[{NZ}]=1"
        ),
        "coefficient_residues": [str(value % modulus) for value in values],
        "input_sha256": input_hashes,
        "large_Groebner_required": False,
    }
    checkpoint.parent.mkdir(parents=True, exist_ok=True)
    checkpoint.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def reconstruct_scalar(value, modulus):
    try:
        return QQ((ZZ(value) % modulus).rational_reconstruction(modulus))
    except (ArithmeticError, ValueError, ZeroDivisionError):
        return None


def reduced_q_scalar(value, modulus):
    value = QQ(value)
    numerator = ZZ(value.numerator()) % modulus
    denominator = ZZ(value.denominator()) % modulus
    if denominator.gcd(modulus) != 1:
        raise ZeroDivisionError("rational denominator is not a modular unit")
    return numerator * denominator.inverse_mod(modulus) % modulus


def exact_candidate():
    modulus = p**exponent
    # Only X and the non-leading coefficients of monic Z need rational
    # reconstruction.  Once they are exact, Y is the unique polynomial square
    # root of the Weierstrass right-hand side with seed-compatible sign.  This
    # avoids expensive, redundant reconstruction of all NY coefficients of Y.
    reconstructed_X = [reconstruct_scalar(values[index], modulus) for index in range(NX)]
    reconstructed_Z = [
        reconstruct_scalar(values[NX + NY + index], modulus)
        for index in range(NZ)
    ]
    print(
        f"{CONFIG['label']}RECONSTRUCT|exponent={exponent}"
        f"|X={sum(value is not None for value in reconstructed_X)}/{NX}"
        f"|Z={sum(value is not None for value in reconstructed_Z)}/{NZ}"
    )
    if any(value is None for value in reconstructed_X + reconstructed_Z):
        return None
    RQ = PolynomialRing(QQ, CONFIG["variable"])
    X = RQ(reconstructed_X)
    Z = RQ(reconstructed_Z + [QQ.one()])
    A = RQ(A_QQ)
    B = RQ(B_QQ)
    right = X**3 + A * X * Z**4 + B * Z**6
    degree = right.degree()
    if degree != 2 * (NY - 1) or not QQ(right[degree]).is_square():
        return None
    leading = QQ(right[degree]).sqrt()
    # Choose the square-root sign compatible with the certified seed.
    seed_leading = seed_values[NX + NY - 1] % p
    if reduced_q_scalar(leading, p) != seed_leading:
        leading = -leading
    coefficients = [QQ.zero()] * NY
    coefficients[-1] = leading
    for target_degree in range(degree - 1, NY - 2, -1):
        unknown_index = target_degree - (NY - 1)
        known = sum(
            coefficients[i] * coefficients[target_degree - i]
            for i in range(NY)
            if 0 <= target_degree - i < NY
            and i != unknown_index
            and target_degree - i != unknown_index
        )
        coefficients[unknown_index] = (QQ(right[target_degree]) - known) / (2 * leading)
    Y = RQ(coefficients)
    if Y**2 != right:
        return None
    return X, Y, Z


if any(int(value) for value in residual(values, p)):
    raise ArithmeticError("seed misses the normalized section equations")

answer = None
while exponent < args.max_exponent:
    next_exponent = min(2 * exponent, args.max_exponent)
    current_modulus = p**exponent
    quotient_modulus = p**(next_exponent - exponent)
    target_modulus = current_modulus * quotient_modulus

    target_residual = [ZZ(value) for value in residual(values, target_modulus)]
    if any(value % current_modulus for value in target_residual):
        raise ArithmeticError("current point lost its certified Hensel precision")
    rhs_all = [(-value // current_modulus) % quotient_modulus for value in target_residual]
    correction = structured_correction(values, quotient_modulus, rhs_all)
    values = [
        (values[index] + current_modulus * ZZ(correction[index])) % target_modulus
        for index in range(unknown_count)
    ]
    exponent = next_exponent
    if any(int(value) for value in residual(values, target_modulus)):
        raise ArithmeticError("unused equation or Newton update failed at raised precision")
    write_checkpoint()
    print(
        f"{CONFIG['label']}HENSEL|exponent={exponent}|bits={target_modulus.nbits()}|"
        f"rank={unknown_count}|all_equations={equation_count}|status=PASS_STEP",
        flush=True,
    )
    if exponent >= args.reconstruct_min_exponent:
        answer = exact_candidate()
        if answer is not None:
            break

if answer is None and exponent >= args.reconstruct_min_exponent:
    answer = exact_candidate()

if answer is None:
    print(
        f"{CONFIG['label']}HENSEL_RESULT|exponent={exponent}|bits={(p**exponent).nbits()}|"
        f"status=NEED_MORE_{args.target.upper().replace('-', '_')}_HENSEL_PRECISION",
        flush=True,
    )
    print(f"CHECKPOINT|{checkpoint}", flush=True)
    raise SystemExit(0)

X, Y, Z = answer
max_bits = max(
    max(abs(ZZ(value.numerator())).nbits(), abs(ZZ(value.denominator())).nbits())
    for poly in (X, Y, Z)
    for value in poly
)
payload = {
    "schema": CONFIG["output_schema"],
    "status": CONFIG["output_status"],
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (SEED, MODEL_PATH)],
        "sha256": input_hashes,
        "hensel_checkpoint": str(checkpoint.relative_to(ROOT)),
        "prime": int(p),
        "exponent": exponent,
        "modulus_bits": int((p**exponent).nbits()),
    },
    "section": {
        "X_coefficients_low_to_high": [str(value) for value in X.list()],
        "Y_coefficients_low_to_high": [str(value) for value in Y.list()],
        "Z_coefficients_low_to_high": [str(value) for value in Z.list()],
        "degrees_X_Y_Z": [int(X.degree()), int(Y.degree()), int(Z.degree())],
        **CONFIG["section_metadata"],
        "max_coefficient_bits": int(max_bits),
        "exact_QQ_weierstrass_identity": True,
    },
    "method": (
        f"quadratically convergent {unknown_count}-coefficient Hensel/Newton lift "
        f"with {NX + NZ}-variable polynomial-remainder corrections"
    ),
    "large_Groebner_required": False,
    "proof_boundary": CONFIG["proof_boundary"],
}
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "{}QQ|degrees={},{},{}|max_bits={}|hensel_exponent={}|status={}".format(
        CONFIG["label"], X.degree(), Y.degree(), Z.degree(), max_bits, exponent, payload["status"]
    ),
    flush=True,
)
print(f"OUTPUT|{output.resolve()}", flush=True)
