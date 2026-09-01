#!/usr/bin/env sage -python
"""Hensel-lift the exact q12 cubic discriminant factors at p=19."""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, Matrix, PolynomialRing, Qp, infinity


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
DEFAULT_INPUT = RESULTS / "q80-third-q12-exact-pencil-p19-adic-precision64.json"
DEFAULT_OUTPUT = RESULTS / "q80-third-q12-discriminant-factors-p19-adic-precision5.json"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
parser.add_argument("--digits", type=int, default=5)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()
args.input = args.input.resolve()
args.output = args.output.resolve()
if args.digits < 2:
    raise ValueError("at least two Hensel digits are required")


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


payload = json.loads(args.input.read_text())
if payload.get("status") != "PASS_EXACT_THIRD_Q12_PENCIL_REDUCTION_MOD_19_POWER":
    raise ValueError("high-precision exact-pencil reduction is not certified")
if payload["specialization"]["precision_digits"] < args.digits + 2:
    raise ValueError("input does not have enough p-adic precision")

prime = 19
working_precision = args.digits + 3
omega_square_residue = payload["quadratic_field"]["omega_square_mod_19"]

# Finite seed and its fixed nine-by-nine Hensel Jacobian inverse.
prime_field = GF(prime)
finite_modulus_ring = PolynomialRing(prime_field, "z_finite")
z_finite = finite_modulus_ring.gen()
finite = GF(
    prime**2,
    "omega_bar",
    modulus=z_finite**2 - prime_field(omega_square_residue),
)
omega_bar = finite.gen()
finite_u_ring = PolynomialRing(finite, "U")
U_bar = finite_u_ring.gen()
finite_u_field = finite_u_ring.fraction_field()
finite_w_ring = PolynomialRing(finite_u_field, "W")
W_bar = finite_w_ring.gen()
finite_x_ring = PolynomialRing(finite_w_ring, "old_x")
finite_coefficients = [finite_w_ring.zero() for unused in range(4)]
for v_degree, w_degree, x_degree, coordinates in payload["pencil"][
    "terms_V_W_old_x_coefficient_1_omega"
]:
    finite_coefficients[x_degree] += (
        (finite(coordinates[0]) + finite(coordinates[1]) * omega_bar)
        * finite_u_field(U_bar) ** v_degree
        * W_bar**w_degree
    )
finite_cubic = finite_x_ring(
    [value / finite_coefficients[3] for value in finite_coefficients]
)
finite_b, finite_c, finite_d = finite_cubic[2], finite_cubic[1], finite_cubic[0]
finite_discriminant = (
    finite_b**2 * finite_c**2
    - 4 * finite_c**3
    - 4 * finite_b**3 * finite_d
    - 27 * finite_d**2
    + 18 * finite_b * finite_c * finite_d
)
factor_by_exponent = {
    int(exponent): factor.monic() for factor, exponent in finite_discriminant.factor()
}
if set(factor_by_exponent) != {1, 2, 3}:
    raise ArithmeticError("finite discriminant does not have exponents 3,2,1")
L_finite = factor_by_exponent[3]
Q_finite = factor_by_exponent[2]
D_finite = factor_by_exponent[1]
if [value.degree() for value in (L_finite, Q_finite, D_finite)] != [1, 4, 4]:
    raise ArithmeticError("finite discriminant factors do not have degrees 1,4,4")

finite_jacobian_columns = (
    [3 * L_finite**2 * Q_finite**2 * D_finite]
    + [2 * L_finite**3 * Q_finite * D_finite * W_bar**index for index in range(4)]
    + [L_finite**3 * Q_finite**2 * W_bar**index for index in range(4)]
)
finite_jacobian = Matrix(
    finite_u_field,
    9,
    9,
    lambda row, column: finite_jacobian_columns[column][row],
)
if finite_jacobian.rank() != 9:
    raise ArithmeticError("discriminant-factor Hensel Jacobian is singular")
finite_jacobian_inverse = finite_jacobian.inverse()

# Repeated root modulo the squarefree conductor L*Q.  Recover separately
# modulo L and Q because their product quotient has zero divisors.
finite_x_over_u = PolynomialRing(finite_u_field, "x_repeat")
w_at_L = -L_finite[0] / L_finite[1]
cubic_at_L = finite_x_over_u(
    [coefficient(w_at_L) for coefficient in finite_cubic.list()]
)
gcd_at_L = cubic_at_L.gcd(cubic_at_L.derivative()).monic()
if gcd_at_L.degree() != 1:
    raise ArithmeticError("finite cubic has no unique repeated root modulo L")
root_at_L = -gcd_at_L[0]
mod_Q = finite_w_ring.quotient(Q_finite, "W_mod_Q")
x_mod_Q = PolynomialRing(mod_Q, "x_repeat_Q")
cubic_mod_Q = x_mod_Q([mod_Q(coefficient) for coefficient in finite_cubic.list()])
gcd_at_Q = cubic_mod_Q.gcd(cubic_mod_Q.derivative()).monic()
if gcd_at_Q.degree() != 1:
    raise ArithmeticError("finite cubic has no unique repeated root modulo Q")
root_at_Q = -gcd_at_Q[0]
finite_conductor = L_finite * Q_finite
finite_repeated_root = finite_w_ring(
    root_at_L
    + L_finite * ((root_at_Q - mod_Q(root_at_L)) / mod_Q(L_finite)).lift()
) % finite_conductor
finite_second_derivative_inverse = (
    finite_w_ring(finite_cubic.derivative().derivative()(finite_repeated_root))
    % finite_conductor
).inverse_mod(finite_conductor)
finite_integral_A = (finite_cubic[2] + finite_repeated_root) % finite_conductor
finite_integral_B = (
    -finite_cubic[2] * finite_repeated_root - 2 * finite_repeated_root**2
) % finite_conductor

# The same exact source over the unramified p-adic quadratic field.
padic = Qp(prime, prec=working_precision, type="capped-rel")
padic_modulus_ring = PolynomialRing(padic, "z_padic")
z_padic = padic_modulus_ring.gen()
padic_field = padic.extension(
    z_padic**2 - padic(payload["quadratic_field"]["omega_square_modulus"]),
    names="omega",
)
omega = padic_field.gen()
padic_u_ring = PolynomialRing(padic_field, "U")
U = padic_u_ring.gen()
padic_u_field = padic_u_ring.fraction_field()
padic_w_ring = PolynomialRing(padic_u_field, "W")
W = padic_w_ring.gen()
padic_x_ring = PolynomialRing(padic_w_ring, "old_x")
padic_coefficients = [padic_w_ring.zero() for unused in range(4)]
for v_degree, w_degree, x_degree, coordinates in payload["pencil"][
    "terms_V_W_old_x_coefficient_1_omega"
]:
    padic_coefficients[x_degree] += (
        (padic_field(coordinates[0]) + padic_field(coordinates[1]) * omega)
        * padic_u_field(U) ** v_degree
        * W**w_degree
    )
padic_cubic = padic_x_ring(
    [value / padic_coefficients[3] for value in padic_coefficients]
)
padic_b, padic_c, padic_d = padic_cubic[2], padic_cubic[1], padic_cubic[0]
padic_discriminant = (
    padic_b**2 * padic_c**2
    - 4 * padic_c**3
    - 4 * padic_b**3 * padic_d
    - 27 * padic_d**2
    + 18 * padic_b * padic_c * padic_d
)
target = padic_discriminant / padic_discriminant.leading_coefficient()


def lift_finite_constant(value):
    coordinates = list(finite(value))
    coordinates += [prime_field.zero()] * (2 - len(coordinates))
    return padic_field(int(coordinates[0])) + padic_field(int(coordinates[1])) * omega


def lift_finite_rational(value):
    value = finite_u_field(value)
    numerator = padic_u_ring(
        [lift_finite_constant(coefficient) for coefficient in value.numerator().list()]
    )
    denominator = padic_u_ring(
        [lift_finite_constant(coefficient) for coefficient in value.denominator().list()]
    )
    return padic_u_field(numerator) / padic_u_field(denominator)


def lift_finite_polynomial(value):
    return padic_w_ring([lift_finite_rational(coefficient) for coefficient in value.list()])


L, Q, D = map(lift_finite_polynomial, (L_finite, Q_finite, D_finite))


def polynomial_valuation(polynomial):
    values = [padic_field(value).valuation() for value in polynomial.list() if value]
    return min(values) if values else infinity


def rational_valuation(value):
    value = padic_u_field(value)
    return polynomial_valuation(value.numerator()) - polynomial_valuation(value.denominator())


def w_polynomial_valuation(polynomial):
    values = [rational_valuation(value) for value in polynomial.list() if value]
    return min(values) if values else infinity


def residue_constant_scaled(value, shift):
    coordinates = padic_field(value).polynomial().list()
    coordinates += [padic.zero()] * (2 - len(coordinates))
    return finite(
        [prime_field((coordinates[index] / prime**shift).residue()) for index in range(2)]
    )


def residue_rational_scaled(value, shift):
    value = padic_u_field(value)
    numerator = finite_u_ring(
        [residue_constant_scaled(coefficient, shift) for coefficient in value.numerator().list()]
    )
    denominator = finite_u_ring(
        [residue_constant_scaled(coefficient, 0) for coefficient in value.denominator().list()]
    )
    return finite_u_field(numerator) / finite_u_field(denominator)


hensel_sample_points = [
    padic_field(a_value) + padic_field(b_value) * omega
    for b_value in range(prime)
    for a_value in range(prime)
]


def compress_rational_hensel(value, numerator_degree, denominator_degree):
    """Compact one intermediate digit lift; final certification is below."""
    value = padic_u_field(value)
    unknown_count = numerator_degree + 1 + denominator_degree
    rows = []
    right_hand_side = []
    for point in hensel_sample_points:
        denominator_value = value.denominator()(point)
        if not denominator_value or denominator_value.valuation() != 0:
            continue
        image = value.numerator()(point) / denominator_value
        rows.append(
            [point**index for index in range(numerator_degree + 1)]
            + [-image * point**index for index in range(denominator_degree)]
        )
        right_hand_side.append(image * point**denominator_degree)
        if len(rows) == unknown_count:
            break
    if len(rows) != unknown_count:
        raise ArithmeticError("insufficient intermediate Hensel compression points")
    solution = Matrix(padic_field, rows).solve_right(
        Matrix(padic_field, unknown_count, 1, right_hand_side)
    ).column(0)
    numerator = padic_u_ring(
        [padic_field(solution[index]) for index in range(numerator_degree + 1)]
    )
    denominator = padic_u_ring(
        [
            padic_field(solution[index])
            for index in range(numerator_degree + 1, unknown_count)
        ]
        + [padic_field.one()]
    )
    return padic_u_field(numerator) / padic_u_field(denominator)


def compress_w_hensel(value, finite_value):
    return padic_w_ring(
        [
            compress_rational_hensel(
                value[index],
                int(finite_u_field(finite_value[index]).numerator().degree()),
                int(finite_u_field(finite_value[index]).denominator().degree()),
            )
            for index in range(finite_value.degree() + 1)
        ]
    )


valuation_history = []
while True:
    residual = L**3 * Q**2 * D - target
    valuation = int(w_polynomial_valuation(residual))
    valuation_history.append(valuation)
    print(f"Q80THIRDQ12PADICFACTORS_PROGRESS|stage=factors|valuation={valuation}", flush=True)
    if valuation >= args.digits:
        break
    right_hand_side = Matrix(
        finite_u_field,
        9,
        1,
        [-residue_rational_scaled(residual[index], valuation) for index in range(9)],
    )
    correction = finite_jacobian_inverse * right_hand_side
    scale = padic(prime) ** valuation
    lifted = [scale * lift_finite_rational(correction[index, 0]) for index in range(9)]
    L = padic_w_ring([L[0] + lifted[0], padic_u_field.one()])
    Q = padic_w_ring(
        [Q[index] + lifted[1 + index] for index in range(4)] + [padic_u_field.one()]
    )
    D = padic_w_ring(
        [D[index] + lifted[5 + index] for index in range(4)] + [padic_u_field.one()]
    )
    L = compress_w_hensel(L, L_finite)
    Q = compress_w_hensel(Q, Q_finite)
    D = compress_w_hensel(D, D_finite)

# Lift the repeated root of the cubic derivative modulo the lifted conductor.
conductor = L * Q
repeated_root = lift_finite_polynomial(finite_repeated_root)
root_valuation_history = []
while True:
    root_residual = padic_w_ring(padic_cubic.derivative()(repeated_root)) % conductor
    root_valuation = int(w_polynomial_valuation(root_residual))
    root_valuation_history.append(root_valuation)
    print(f"Q80THIRDQ12PADICFACTORS_PROGRESS|stage=root|valuation={root_valuation}", flush=True)
    if root_valuation >= args.digits:
        break
    finite_residual = finite_w_ring(
        [
            residue_rational_scaled(coefficient, root_valuation)
            for coefficient in root_residual.list()
        ]
    )
    finite_correction = (
        -finite_residual * finite_second_derivative_inverse
    ) % finite_conductor
    repeated_root = (
        repeated_root
        + padic(prime) ** root_valuation * lift_finite_polynomial(finite_correction)
    ) % conductor
    repeated_root = compress_w_hensel(repeated_root, finite_repeated_root)

integral_A = (padic_b + repeated_root) % conductor
integral_B = (-padic_b * repeated_root - 2 * repeated_root**2) % conductor

# Digit corrections over a generic fraction field accumulate artificial
# denominator products.  The mod-19 solution pins the true numerator and
# denominator degrees.  Recover that compact support by p-adic rational
# interpolation and verify on held-out residue classes.
sample_points = [
    padic_field(a_value) + padic_field(b_value) * omega
    for b_value in range(prime)
    for a_value in range(prime)
]


def compress_rational(value, numerator_degree, denominator_degree):
    value = padic_u_field(value)
    unknown_count = numerator_degree + 1 + denominator_degree
    rows = []
    right_hand_side = []
    used = []
    for point in sample_points:
        denominator_value = value.denominator()(point)
        if not denominator_value or denominator_value.valuation() != 0:
            continue
        image = value.numerator()(point) / denominator_value
        row = [point**index for index in range(numerator_degree + 1)]
        row += [-image * point**index for index in range(denominator_degree)]
        rows.append(row)
        right_hand_side.append(image * point**denominator_degree)
        used.append((point, image))
        if len(rows) == unknown_count:
            break
    if len(rows) != unknown_count:
        raise ArithmeticError("insufficient unit-denominator interpolation points")
    solution = Matrix(padic_field, rows).solve_right(Matrix(padic_field, unknown_count, 1, right_hand_side)).column(0)
    numerator = padic_u_ring(
        [padic_field(solution[index]) for index in range(numerator_degree + 1)]
    )
    denominator = padic_u_ring(
        [
            padic_field(solution[index])
            for index in range(numerator_degree + 1, unknown_count)
        ]
        + [padic_field.one()]
    )
    result = padic_u_field(numerator) / padic_u_field(denominator)
    held_out = 0
    for point in sample_points:
        old_denominator = value.denominator()(point)
        new_denominator = result.denominator()(point)
        if (
            not old_denominator
            or not new_denominator
            or old_denominator.valuation() != 0
            or new_denominator.valuation() != 0
            or any(point == used_point for used_point, unused in used)
        ):
            continue
        difference = (
            value.numerator()(point) / old_denominator
            - result.numerator()(point) / new_denominator
        )
        if difference and difference.valuation() < args.digits:
            raise ArithmeticError("held-out rational interpolation replay failed")
        held_out += 1
        if held_out == 3:
            break
    if held_out != 3:
        raise ArithmeticError("insufficient held-out interpolation points")
    return result


def compress_w_polynomial(value, finite_value):
    if value.degree() != finite_value.degree():
        raise ArithmeticError("p-adic/finite W degrees disagree during compression")
    result = []
    for index in range(value.degree() + 1):
        finite_coefficient = finite_u_field(finite_value[index])
        result.append(
            compress_rational(
                value[index],
                int(finite_coefficient.numerator().degree()),
                int(finite_coefficient.denominator().degree()),
            )
        )
    return padic_w_ring(result)


L = compress_w_polynomial(L, L_finite)
print("Q80THIRDQ12PADICFACTORS_PROGRESS|stage=compress|object=L", flush=True)
Q = compress_w_polynomial(Q, Q_finite)
print("Q80THIRDQ12PADICFACTORS_PROGRESS|stage=compress|object=Q", flush=True)
D = compress_w_polynomial(D, D_finite)
print("Q80THIRDQ12PADICFACTORS_PROGRESS|stage=compress|object=D", flush=True)
integral_A = compress_w_polynomial(integral_A, finite_integral_A)
print("Q80THIRDQ12PADICFACTORS_PROGRESS|stage=compress|object=A", flush=True)
integral_B = compress_w_polynomial(integral_B, finite_integral_B)
print("Q80THIRDQ12PADICFACTORS_PROGRESS|stage=compress|object=B", flush=True)
conductor = L * Q
compressed_factor_residual_valuation = int(
    w_polynomial_valuation(L**3 * Q**2 * D - target)
)
if compressed_factor_residual_valuation < args.digits:
    raise ArithmeticError("compressed discriminant factors lose p-adic precision")
if w_polynomial_valuation(integral_B % L) < args.digits:
    raise ArithmeticError("integral-basis constant numerator is not divisible by L")

output_modulus = prime**args.digits


def constant_coordinates(value):
    coordinates = padic_field(value).polynomial().list()
    coordinates += [padic.zero()] * (2 - len(coordinates))
    return [int(coordinates[index].lift()) % output_modulus for index in range(2)]


def rational_record(value):
    value = padic_u_field(value)
    numerator = padic_u_ring(value.numerator())
    denominator = padic_u_ring(value.denominator())
    leading = denominator.leading_coefficient()
    numerator /= leading
    denominator /= leading
    return {
        "numerator_coefficients_low_to_high_U_1_omega": [
            constant_coordinates(coefficient) for coefficient in numerator.list()
        ],
        "denominator_coefficients_low_to_high_U_1_omega": [
            constant_coordinates(coefficient) for coefficient in denominator.list()
        ],
        "degrees_numerator_denominator": [int(numerator.degree()), int(denominator.degree())],
    }


def factor_record(value):
    return {
        "degree_W": int(value.degree()),
        "coefficients_low_to_high_W": [rational_record(coefficient) for coefficient in value.list()],
    }


output = {
    "schema": "elkies-k3.q80-third-q12-discriminant-factors-p19-adic.v1",
    "status": "PASS_EXACT_THIRD_Q12_DISCRIMINANT_FACTOR_HENSEL_LIFT_P19",
    "specialization": {
        "u": "-2",
        "prime": prime,
        "digits": args.digits,
        "modulus": output_modulus,
    },
    "factorization": {
        "identity": "monic_discriminant=L^3*Q^2*D mod 19^digits",
        "degrees_W_L_Q_D": [1, 4, 4],
        "L": factor_record(L),
        "Q": factor_record(Q),
        "D": factor_record(D),
        "valuation_history": valuation_history,
        "final_residual_valuation": valuation_history[-1],
        "compressed_residual_valuation": compressed_factor_residual_valuation,
    },
    "integral_basis_candidate": {
        "candidate_basis": ["1", "z", "e"],
        "candidate_e_formula": "(z^2+A*z+B)/(L*Q)",
        "A": factor_record(integral_A),
        "B": factor_record(integral_B),
        "degrees_W_A_B": [int(integral_A.degree()), int(integral_B.degree())],
        "repeated_root_valuation_history": root_valuation_history,
        "B_divisible_by_L_to_requested_precision": True,
        "generic_integrality_check": "pending trace/second-symmetric/determinant divisibility worker",
        "compression": "mod-19 numerator/denominator degree bounds with three held-out p-adic evaluations",
    },
    "hensel": {
        "unknown_coefficients": 9,
        "finite_jacobian_rank": 9,
        "pivot_coefficient_degrees": list(range(9)),
        "algorithm": "fixed inverse of the mod-19 Jacobian with digit corrections",
    },
    "input": {"path": str(args.input.relative_to(ROOT)), "sha256": sha256(args.input)},
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            "the exact p-adic cubic discriminant has a unique lift of the mod-19 L^3 Q^2 D factorization",
            "the nine-variable factor system has nonsingular mod-19 Jacobian",
            "literal factor identity through the displayed p-adic precision",
            "a lifted repeated root modulo LQ and the uniquely normalized candidate cubic integral-basis element",
        ],
        "not_proved": [
            "an exact characteristic-zero discriminant factorization",
            "generic integrality of the displayed candidate basis element",
            "a p-adic Jacobian or maps",
            "rational reconstruction of L, Q, or D",
        ],
    },
    "reproduce": (
        "sage -python elkies-k3/scripts/lift_q80_third_q12_discriminant_factors_p19_adic.sage "
        f"--digits {args.digits} --output {args.output.relative_to(ROOT)}"
    ),
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
if args.check:
    if not args.output.exists() or args.output.read_text() != serialized:
        raise SystemExit(f"p-adic discriminant-factor artifact is stale: {args.output}")
else:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(serialized)
print(
    f"Q80THIRDQ12PADICFACTORS|p=19|digits={args.digits}|history={valuation_history}|"
    f"root_history={root_valuation_history}|candidate_basis=1,z,e|"
    "degrees=1,4,4|rank=9|status=PASS_EXACT_THIRD_Q12_DISCRIMINANT_FACTOR_HENSEL_LIFT_P19"
)
