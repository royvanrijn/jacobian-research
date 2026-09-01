#!/usr/bin/env sage -python
"""Lift the two p=19 third-q12 closure operands exactly over QQ."""

import argparse
import hashlib
import json
import time
from pathlib import Path

from fpylll import IntegerMatrix, LLL
from sage.all import (
    GF,
    QQ,
    ZZ,
    Zmod,
    PolynomialRing,
    gcd,
    matrix,
    vector,
)


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--surface", type=Path,
    default=RESULTS / "q80-fixed-u-minus2-p19-height-shell-with-po1.json",
)
parser.add_argument(
    "--horizontal", type=Path,
    default=RESULTS / "q80-third-q12-um2-p19-common-producer-horizontal.json",
)
parser.add_argument(
    "--modular-horizontal", type=Path, action="append", default=[],
    help="additional unique-orbit common-producer artifact to replay; repeatable",
)
parser.add_argument("--precision", type=int, default=32768)
parser.add_argument(
    "--quadratic-only", action="store_true",
    help="certify the exact trace-zero quadratic operand without lifting its rational companion",
)
parser.add_argument(
    "--residue-only", action="store_true",
    help="write the high-precision rational-branch Hensel residues without LLL reconstruction",
)
parser.add_argument(
    "--biquadratic-operands", action="store_true",
    help="reconstruct and certify both trace-zero quadratic closure operands",
)
parser.add_argument(
    "--output", type=Path,
    default=RESULTS / "q80-third-q12-um2-closure-operands-p19-hensel-qq.json",
)
args = parser.parse_args()
for name in ("surface", "horizontal", "output"):
    setattr(args, name, getattr(args, name).resolve())
if args.precision < 2:
    raise ValueError("p-adic precision must be at least two")
if not args.modular_horizontal:
    args.modular_horizontal = [
        RESULTS / f"q80-third-q12-um2-p{prime}-common-producer-horizontal.json"
        for prime in (19, 61, 67, 83, 89, 103, 131)
    ]
args.modular_horizontal = [path.resolve() for path in args.modular_horizontal]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


started = time.monotonic()
surface = json.loads(args.surface.read_text())
horizontal = json.loads(args.horizontal.read_text())
if surface.get("schema") != "elkies-k3.q80-fixed-u-marked-third-q12-search.v1":
    raise ValueError("unexpected surface schema")
if horizontal.get("status") != "PASS_EXACT_THIRD_Q12_HORIZONTAL_FROM_COMMON_CLOSURE_PRODUCER":
    raise ValueError("p=19 horizontal is not the certified common-producer artifact")
if horizontal["specialization"] != {
    "common_extension_degree": 2,
    "common_extension_modulus": "x^2 + 4*x + 1",
    "prime": 19,
    "u": "-2",
}:
    raise ValueError("unexpected p=19 horizontal specialization")

parameter = surface["parameters"][0]
equation = parameter["exact_equations"]["second_q4"]
W_ring = PolynomialRing(QQ, "W")
A = W_ring(equation["A_coefficients_low_to_high"])
B = W_ring(equation["B_coefficients_low_to_high"])

# Rebuild the exact closure equations used by the prime-independent producer.
closure_ring = PolynomialRing(QQ, names=("l", "x0", "x1", "x2", "x3"), order="degrevlex")
l, x0, x1, x2, x3 = closure_ring.gens()
closure_fraction = closure_ring.fraction_field()
section_ring = PolynomialRing(closure_fraction, "W_section")
W_section = section_ring.gen()


def lift_polynomial(value):
    return section_ring([closure_fraction(coefficient) for coefficient in value.list()])


X = l**2 * W_section**4 + x3 * W_section**3 + x2 * W_section**2 + x1 * W_section + x0
square = X**3 + lift_polynomial(A) * X + lift_polynomial(B)
y_coefficients = [closure_fraction.zero() for unused in range(7)]
y_coefficients[6] = closure_fraction(l**3)
for degree in range(11, 5, -1):
    index = degree - 6
    partial = sum(y_coefficients[j] * W_section**j for j in range(7))
    y_coefficients[index] = (
        square[degree] - (partial**2)[degree]
    ) / (2 * y_coefficients[6])
Y = sum(y_coefficients[j] * W_section**j for j in range(7))
identity = Y**2 - square
closure_equations = [closure_ring(identity[index].numerator()) for index in range(6)]
if any(identity[index] for index in range(6, 13)):
    raise ArithmeticError("exact top-down closure recursion did not close")

# The selected quadratic operand has trace(l)=0 and rational x coefficients.
# Every closure equation is even in l, so q=l^2 gives a rational five-variable
# system rather than a quadratic p-adic computation.
trace_zero_ring = PolynomialRing(QQ, names=("q", "x0", "x1", "x2", "x3"), order="degrevlex")
q, *trace_zero_x = trace_zero_ring.gens()
trace_zero_equations = []
for polynomial in closure_equations:
    reduced = trace_zero_ring.zero()
    for exponents, coefficient in polynomial.dict().items():
        l_exponent, *x_exponents = exponents
        if l_exponent % 2:
            raise ArithmeticError("closure equation is not even on the trace-zero branch")
        term = coefficient * q**(l_exponent // 2)
        for variable, exponent in zip(trace_zero_x, x_exponents):
            term *= variable**exponent
        reduced += term
    trace_zero_equations.append(reduced)

candidate = horizontal["pairwise_producer"]["candidates"][0]
if horizontal["pairwise_producer"]["frobenius_orbits_up_to_sign"] != 1:
    raise ValueError("p=19 horizontal no longer has one target orbit")
left = horizontal["decoded_sections"][candidate["left_index_zero_based"]]
right = horizontal["decoded_sections"][candidate["right_index_zero_based"]]
left_values = left["l_x0_x1_x2_x3_sat"]
right_values = right["l_x0_x1_x2_x3_sat"]
if any(value[1] for value in left_values):
    raise ArithmeticError("selected rational operand is not rational mod 19")
if any(right_values[index][1] for index in range(1, 5)):
    raise ArithmeticError("selected trace-zero operand has nonrational x coefficients mod 19")
quadratic_seed = [
    ZZ(right_values[0][1]) * ZZ(right_values[5][1]).inverse_mod(19) % 19,
    *[ZZ(right_values[index][0]) for index in range(1, 5)],
]
rational_seed = [ZZ(left_values[index][0]) for index in range(5)]


def independent_rows(polynomials, ring, seed, prime=19):
    finite_ring = ring.change_ring(GF(prime))
    finite_polynomials = [finite_ring(polynomial) for polynomial in polynomials]
    jacobian = matrix(
        GF(prime),
        [
            [polynomial.derivative(variable)(*seed) for variable in finite_ring.gens()]
            for polynomial in finite_polynomials
        ],
    )
    if jacobian.rank() != len(seed):
        raise ArithmeticError("selected closure branch is not Jacobian-regular mod 19")
    row_indices = tuple(jacobian.transpose().pivots())
    if len(row_indices) != len(seed):
        raise ArithmeticError("could not select a square Hensel subsystem")
    return row_indices, int(jacobian.rank())


quadratic_rows, quadratic_rank = independent_rows(
    trace_zero_equations, trace_zero_ring, quadratic_seed
)
rational_rows, rational_rank = independent_rows(
    closure_equations, closure_ring, rational_seed
)


def hensel_lift(polynomials, ring, row_indices, seed, prime, precision, reconstruct=True):
    selected = [polynomials[index] for index in row_indices]
    values = [ZZ(value) % prime for value in seed]
    exponent = 1
    iterations = []
    while exponent < precision:
        next_exponent = min(2 * exponent, precision)
        modulus = ZZ(prime) ** next_exponent
        residue_ring = Zmod(modulus)
        modular_ring = ring.change_ring(residue_ring)
        modular_polynomials = [modular_ring(polynomial) for polynomial in selected]
        modular_values = [residue_ring(value) for value in values]
        jacobian = matrix(
            residue_ring,
            [
                [polynomial.derivative(variable)(*modular_values) for variable in modular_ring.gens()]
                for polynomial in modular_polynomials
            ],
        )
        residual = vector(
            residue_ring,
            [polynomial(*modular_values) for polynomial in modular_polynomials],
        )
        correction = jacobian.solve_right(residual)
        lifted = vector(residue_ring, modular_values) - correction
        if any(polynomial(*lifted) for polynomial in modular_polynomials):
            raise ArithmeticError("Newton-Hensel correction failed literal replay")
        values = [ZZ(value) for value in lifted]
        exponent = next_exponent
        iterations.append({
            "precision_exponent": exponent,
            "modulus_bits": int(modulus.nbits()),
        })
    modulus = ZZ(prime) ** precision
    if not reconstruct:
        return values, modulus, iterations
    reconstructed = []
    for value in values:
        reconstructed.append(Zmod(modulus)(value).rational_reconstruction())
    return reconstructed, modulus, iterations


quadratic_values, hensel_modulus, quadratic_iterations = hensel_lift(
    trace_zero_equations,
    trace_zero_ring,
    quadratic_rows,
    quadratic_seed,
    19,
    args.precision,
)
if any(polynomial(*quadratic_values) for polynomial in trace_zero_equations):
    raise ArithmeticError("reconstructed trace-zero operand fails exact closure equations")
q_value, *quadratic_x = map(QQ, quadratic_values)
if not q_value or q_value.is_square():
    raise ArithmeticError("reconstructed quadratic operand has the wrong field shape")


def reduce_rational_early(value, prime):
    value = QQ(value)
    if value.denominator() % prime == 0:
        raise ZeroDivisionError(f"reconstructed denominator vanishes at p={prime}")
    return int(value.numerator() * ZZ(value.denominator()).inverse_mod(prime) % prime)


if args.quadratic_only:
    modular_replays = []
    for path in args.modular_horizontal:
        payload = json.loads(path.read_text())
        if payload.get("status") != "PASS_EXACT_THIRD_Q12_HORIZONTAL_FROM_COMMON_CLOSURE_PRODUCER":
            raise ValueError(f"uncertified modular horizontal: {path}")
        if payload["pairwise_producer"]["frobenius_orbits_up_to_sign"] != 1:
            raise ValueError(f"modular horizontal is not unique-orbit: {path}")
        prime = int(payload["specialization"]["prime"])
        q_residue = reduce_rational_early(q_value, prime)
        finite = GF(prime)
        if not q_residue or finite(q_residue).is_square():
            raise ArithmeticError(f"exact quadratic field is not inert at p={prime}")
        modulus = PolynomialRing(finite, "x")(payload["specialization"]["common_extension_modulus"])
        local_discriminant = modulus[1]**2 - 4 * modulus[0]
        ratio = finite(q_residue) / local_discriminant
        if not ratio.is_square():
            raise ArithmeticError(f"exact/local quadratic fields disagree at p={prime}")
        operand_replayed = False
        if prime == 19:
            modular_candidate = payload["pairwise_producer"]["candidates"][0]
            modular_right = payload["decoded_sections"][modular_candidate["right_index_zero_based"]]
            coordinates = modular_right["l_x0_x1_x2_x3_sat"]
            observed = [
                int(ZZ(coordinates[0][1]) * ZZ(coordinates[5][1]).inverse_mod(prime) % prime),
                *[int(coordinates[index][0]) for index in range(1, 5)],
            ]
            expected = [reduce_rational_early(value, prime) for value in quadratic_values]
            operand_replayed = observed == expected
            if not operand_replayed:
                raise ArithmeticError("exact quadratic operand does not replay at p=19")
        modular_replays.append({
            "prime": prime,
            "path": str(path.relative_to(ROOT)),
            "sha256": sha256(path),
            "quadratic_field_inert_and_isomorphic": True,
            "quadratic_operand_replayed": operand_replayed,
        })

    def early_record(value):
        value = QQ(value)
        return {"numerator": str(value.numerator()), "denominator": str(value.denominator())}

    def early_bits(value):
        value = QQ(value)
        return max(abs(value.numerator()).nbits(), value.denominator().nbits())

    output = {
        "schema": "elkies-k3.q80-third-q12-quadratic-closure-operand-p19-hensel-qq.v1",
        "status": "PASS_EXACT_QQ_THIRD_Q12_QUADRATIC_CLOSURE_OPERAND_P19_HENSEL",
        "specialization": {"u": "-2"},
        "quadratic_field": {
            "generator": "l",
            "polynomial": f"z^2 - ({q_value})",
            "l_squared": early_record(q_value),
            "nonsquare_over_QQ": True,
        },
        "quadratic_operand": {
            "l": "l",
            "x_coefficients_low_to_high": [
                *[early_record(value) for value in quadratic_x],
                early_record(q_value),
            ],
            "sat": f"l/({q_value})",
            "trace_l": "0",
            "exact_closure_substitution": True,
        },
        "hensel": {
            "prime": 19,
            "precision_exponent": args.precision,
            "modulus_bits": int(hensel_modulus.nbits()),
            "jacobian_rank": quadratic_rank,
            "rows_zero_based": list(quadratic_rows),
            "iterations": quadratic_iterations,
        },
        "coefficient_height": {
            "max_bits": max(map(early_bits, quadratic_values)),
        },
        "checks": {
            "exact_six_equation_quadratic_branch": True,
            "literal_operand_replay_primes": [19],
            "compatible_inert_field_primes": [record["prime"] for record in modular_replays],
            "modular_field_comparisons": modular_replays,
        },
        "inputs": {
            "surface": {"path": str(args.surface.relative_to(ROOT)), "sha256": sha256(args.surface)},
            "p19_horizontal": {
                "path": str(args.horizontal.relative_to(ROOT)),
                "sha256": sha256(args.horizontal),
            },
        },
        "worker": {
            "path": str(Path(__file__).resolve().relative_to(ROOT)),
            "sha256": sha256(Path(__file__).resolve()),
        },
        "claim_boundary": {
            "proved": [
                "exact quadratic number field generated by the trace-zero polynomial-section coefficient l",
                "exact quadratic closure operand at u=-2",
                "literal exact substitution into all six characteristic-zero closure equations",
                "literal operand reduction replay at p=19",
                "abstract inert residue-field compatibility at every listed reconstruction prime",
            ],
            "not_proved": [
                "the exact rational companion closure operand",
                "the exact sum as a target horizontal section over the quadratic field",
                "a characteristic-zero child equation, maps, marking, or Mordell--Weil rank",
            ],
        },
        "runtime_seconds": time.monotonic() - started,
        "reproduce": (
            "sage -python elkies-k3/scripts/lift_q80_third_q12_closure_operands_p19_qq.sage "
            f"--precision {args.precision} --quadratic-only --output {args.output}"
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        f"Q80THIRDQ12QUADRATICHENSEL|precision={args.precision}|bits={hensel_modulus.nbits()}|"
        f"coefficient_bits={output['coefficient_height']['max_bits']}|"
        f"field_primes={','.join(map(str, output['checks']['compatible_inert_field_primes']))}|"
        "status=PASS_EXACT_QQ_THIRD_Q12_QUADRATIC_CLOSURE_OPERAND_P19_HENSEL",
        flush=True,
    )
    raise SystemExit(0)

rational_residues, rational_modulus, rational_iterations = hensel_lift(
    closure_equations,
    closure_ring,
    rational_rows,
    rational_seed,
    19,
    args.precision,
    reconstruct=False,
)
if rational_modulus != hensel_modulus:
    raise ArithmeticError("Hensel moduli disagree")

if args.residue_only:
    output = {
        "schema": "elkies-k3.q80-third-q12-closure-operands-p19-hensel-residues.v1",
        "status": "PARTIAL_EXACT_THIRD_Q12_RATIONAL_OPERAND_HENSEL_RESIDUES",
        "specialization": {"u": "-2"},
        "prime": 19,
        "precision_exponent": args.precision,
        "modulus": str(rational_modulus),
        "modulus_bits": int(rational_modulus.nbits()),
        "variable_order": [str(variable) for variable in closure_ring.gens()],
        "residues": [str(value) for value in rational_residues],
        "quadratic_operand_exact": {
            "l_squared": str(q_value),
            "x_coefficients_low_to_high_without_leading": [str(value) for value in quadratic_x],
        },
        "hensel": {
            "rational_branch_jacobian_rank": rational_rank,
            "rational_branch_rows_zero_based": list(rational_rows),
            "rational_iterations": rational_iterations,
            "literal_selected_equation_replay_modulo_p_power": True,
        },
        "inputs": {
            "surface": {"path": str(args.surface.relative_to(ROOT)), "sha256": sha256(args.surface)},
            "p19_horizontal": {
                "path": str(args.horizontal.relative_to(ROOT)),
                "sha256": sha256(args.horizontal),
            },
        },
        "worker": {
            "path": str(Path(__file__).resolve().relative_to(ROOT)),
            "sha256": sha256(Path(__file__).resolve()),
        },
        "claim_boundary": {
            "proved": [
                "unique Jacobian-regular p-adic lift of the selected p=19 rational closure operand to the stated precision",
                "literal replay of the selected five closure equations modulo the full p-power",
            ],
            "not_proved": [
                "rational reconstruction of the companion operand",
                "exact substitution over QQ",
                "the characteristic-zero target horizontal or child",
            ],
        },
        "runtime_seconds": time.monotonic() - started,
        "reproduce": (
            "sage -python elkies-k3/scripts/lift_q80_third_q12_closure_operands_p19_qq.sage "
            f"--precision {args.precision} --residue-only --output {args.output}"
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        f"Q80THIRDQ12HENSELRESIDUES|precision={args.precision}|bits={rational_modulus.nbits()}|"
        "status=PARTIAL_EXACT_THIRD_Q12_RATIONAL_OPERAND_HENSEL_RESIDUES",
        flush=True,
    )
    raise SystemExit(0)

if args.biquadratic_operands:
    second_values = [
        Zmod(rational_modulus)(rational_residues[0]**2).rational_reconstruction(),
        *[
            Zmod(rational_modulus)(value).rational_reconstruction()
            for value in rational_residues[1:]
        ],
    ]
    if any(polynomial(*second_values) for polynomial in trace_zero_equations):
        raise ArithmeticError("reconstructed second quadratic operand fails exact closure equations")
    q2_value, *second_x = map(QQ, second_values)
    if not q2_value or q2_value.is_square() or (q2_value / q_value).is_square():
        raise ArithmeticError("the two reconstructed quadratic square classes are not independent")

    def bq_record(value):
        value = QQ(value)
        return {"numerator": str(value.numerator()), "denominator": str(value.denominator())}

    def bq_bits(value):
        value = QQ(value)
        return max(abs(value.numerator()).nbits(), value.denominator().nbits())

    local_behavior = []
    literal_replay = False
    for path in args.modular_horizontal:
        payload = json.loads(path.read_text())
        prime = int(payload["specialization"]["prime"])
        q1_residue = reduce_rational_early(q_value, prime)
        q2_residue = reduce_rational_early(q2_value, prime)
        finite = GF(prime)
        symbols = [
            -1 if finite(value).is_square() is False else (0 if not value else 1)
            for value in (q1_residue, q2_residue, q1_residue * q2_residue % prime)
        ]
        record = {
            "prime": prime,
            "q1_q2_q1q2_quadratic_characters": symbols,
            "path": str(path.relative_to(ROOT)),
            "sha256": sha256(path),
        }
        if prime == 19:
            modular_candidate = payload["pairwise_producer"]["candidates"][0]
            modular_first = payload["decoded_sections"][modular_candidate["right_index_zero_based"]]
            modular_second = payload["decoded_sections"][modular_candidate["left_index_zero_based"]]
            first_coordinates = modular_first["l_x0_x1_x2_x3_sat"]
            second_coordinates = modular_second["l_x0_x1_x2_x3_sat"]
            observed_first = [
                int(ZZ(first_coordinates[0][1]) * ZZ(first_coordinates[5][1]).inverse_mod(prime) % prime),
                *[int(first_coordinates[index][0]) for index in range(1, 5)],
            ]
            observed_second = [
                int(second_coordinates[0][0] ** 2 % prime),
                *[int(second_coordinates[index][0]) for index in range(1, 5)],
            ]
            expected_first = [reduce_rational_early(value, prime) for value in quadratic_values]
            expected_second = [reduce_rational_early(value, prime) for value in second_values]
            literal_replay = observed_first == expected_first and observed_second == expected_second
            if not literal_replay:
                raise ArithmeticError("exact biquadratic operands do not replay at p=19")
            record["both_operands_replayed"] = True
        local_behavior.append(record)

    output = {
        "schema": "elkies-k3.q80-third-q12-biquadratic-closure-operands-p19-hensel-qq.v1",
        "status": "PASS_EXACT_QQ_THIRD_Q12_BIQUADRATIC_CLOSURE_OPERANDS_P19_HENSEL",
        "specialization": {"u": "-2"},
        "biquadratic_field": {
            "generators": ["a", "b"],
            "relations": [f"a^2=({q_value})", f"b^2=({q2_value})"],
            "q1": bq_record(q_value),
            "q2": bq_record(q2_value),
            "q1_nonsquare": True,
            "q2_nonsquare": True,
            "q2_over_q1_nonsquare": True,
            "degree": 4,
        },
        "operands": [
            {
                "generator": "a",
                "l_squared": bq_record(q_value),
                "x_coefficients_low_to_high": [
                    *[bq_record(value) for value in quadratic_x], bq_record(q_value)
                ],
                "sat": f"a/({q_value})",
                "exact_six_equation_closure_substitution": True,
            },
            {
                "generator": "b",
                "l_squared": bq_record(q2_value),
                "x_coefficients_low_to_high": [
                    *[bq_record(value) for value in second_x], bq_record(q2_value)
                ],
                "sat": f"b/({q2_value})",
                "exact_six_equation_closure_substitution": True,
            },
        ],
        "hensel": {
            "prime": 19,
            "precision_exponent": args.precision,
            "modulus_bits": int(rational_modulus.nbits()),
            "first_branch_rows_zero_based": list(quadratic_rows),
            "second_branch_rows_zero_based": list(rational_rows),
            "first_iterations": quadratic_iterations,
            "second_iterations": rational_iterations,
        },
        "coefficient_height": {
            "first_max_bits": max(map(bq_bits, quadratic_values)),
            "second_max_bits": max(map(bq_bits, second_values)),
        },
        "checks": {
            "exact_six_equation_first_branch": True,
            "exact_six_equation_second_branch": True,
            "literal_p19_both_operand_replay": literal_replay,
            "local_splitting_behavior": local_behavior,
        },
        "inputs": {
            "surface": {"path": str(args.surface.relative_to(ROOT)), "sha256": sha256(args.surface)},
            "p19_horizontal": {
                "path": str(args.horizontal.relative_to(ROOT)),
                "sha256": sha256(args.horizontal),
            },
        },
        "worker": {
            "path": str(Path(__file__).resolve().relative_to(ROOT)),
            "sha256": sha256(Path(__file__).resolve()),
        },
        "claim_boundary": {
            "proved": [
                "two exact independent quadratic closure operands at u=-2",
                "literal characteristic-zero substitution of both operands into all six closure equations",
                "literal reduction of both operands to the selected p=19 pair",
                "exact local splitting table for the two quadratic square classes",
            ],
            "not_proved": [
                "the exact sum as the target horizontal over the biquadratic field",
                "a characteristic-zero resolved pencil, child equation, maps, marking, or Mordell--Weil rank",
            ],
        },
        "runtime_seconds": time.monotonic() - started,
        "reproduce": (
            "sage -python elkies-k3/scripts/lift_q80_third_q12_closure_operands_p19_qq.sage "
            f"--precision {args.precision} --biquadratic-operands --output {args.output}"
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        f"Q80THIRDQ12BIQUADRATIC|precision={args.precision}|bits={rational_modulus.nbits()}|"
        f"first_bits={output['coefficient_height']['first_max_bits']}|"
        f"second_bits={output['coefficient_height']['second_max_bits']}|"
        "status=PASS_EXACT_QQ_THIRD_Q12_BIQUADRATIC_CLOSURE_OPERANDS_P19_HENSEL",
        flush=True,
    )
    raise SystemExit(0)


def simultaneous_projective_reconstructions(residues, modulus):
    dimension = len(residues)
    basis_rows = [[ZZ(0) for unused in range(dimension + 1)] for unused in range(dimension + 1)]
    for index in range(dimension):
        basis_rows[index][index] = modulus
    basis_rows[dimension] = list(map(ZZ, residues)) + [ZZ(1)]
    basis = IntegerMatrix.from_matrix([[int(value) for value in row] for row in basis_rows])
    LLL.reduction(basis, delta=0.99)
    reduced_rows = [
        vector(ZZ, [ZZ(basis[row, column]) for column in range(dimension + 1)])
        for row in range(dimension + 1)
    ]
    candidates = []
    for row in sorted(reduced_rows, key=lambda value: value * value):
        scale = ZZ(row[-1])
        if not scale or gcd(scale, modulus) != 1:
            continue
        if any((row[index] - scale * residues[index]) % modulus for index in range(dimension)):
            continue
        values = [QQ(row[index]) / scale for index in range(dimension)]
        candidates.append({
            "values": values,
            "projective_scale": scale,
            "primitive_vector_max_bits": max(abs(ZZ(value)).nbits() for value in row),
            "lll_norm_squared": row * row,
        })
    return candidates


rational_projective_candidates = simultaneous_projective_reconstructions(
    rational_residues, rational_modulus
)
valid_rational_candidates = [
    record
    for record in rational_projective_candidates
    if not any(polynomial(*record["values"]) for polynomial in closure_equations)
]
if len(valid_rational_candidates) != 1:
    raise ArithmeticError(
        "projective reconstruction did not produce one exact rational operand "
        f"({len(valid_rational_candidates)} exact among "
        f"{len(rational_projective_candidates)} modular candidates)"
    )
rational_projective = valid_rational_candidates[0]
rational_values = rational_projective["values"]

rational_l, *rational_x = map(QQ, rational_values)
if not q_value or q_value.is_square() or not rational_l:
    raise ArithmeticError("reconstructed operand has the wrong rational/quadratic shape")


def reduce_rational(value, prime):
    value = QQ(value)
    if value.denominator() % prime == 0:
        raise ZeroDivisionError(f"reconstructed denominator vanishes at p={prime}")
    return int(value.numerator() * ZZ(value.denominator()).inverse_mod(prime) % prime)


modular_replays = []
for path in args.modular_horizontal:
    payload = json.loads(path.read_text())
    if payload.get("status") != "PASS_EXACT_THIRD_Q12_HORIZONTAL_FROM_COMMON_CLOSURE_PRODUCER":
        raise ValueError(f"uncertified modular horizontal: {path}")
    if payload["pairwise_producer"]["frobenius_orbits_up_to_sign"] != 1:
        raise ValueError(f"modular horizontal is not unique-orbit: {path}")
    prime = int(payload["specialization"]["prime"])
    modular_candidate = payload["pairwise_producer"]["candidates"][0]
    modular_left = payload["decoded_sections"][modular_candidate["left_index_zero_based"]]
    modular_right = payload["decoded_sections"][modular_candidate["right_index_zero_based"]]
    left_coordinates = modular_left["l_x0_x1_x2_x3_sat"]
    right_coordinates = modular_right["l_x0_x1_x2_x3_sat"]
    if any(value[1] for value in left_coordinates + right_coordinates[1:5]):
        raise ArithmeticError(f"operand rational-coordinate shape changed at p={prime}")
    observed_q = (
        ZZ(right_coordinates[0][1])
        * ZZ(right_coordinates[5][1]).inverse_mod(prime)
        % prime
    )
    expected_quadratic = [
        reduce_rational(q_value, prime),
        *[reduce_rational(value, prime) for value in quadratic_x],
    ]
    observed_quadratic = [
        int(observed_q),
        *[int(right_coordinates[index][0]) for index in range(1, 5)],
    ]
    expected_rational = [reduce_rational(value, prime) for value in rational_values]
    observed_rational = [int(left_coordinates[index][0]) for index in range(5)]
    if expected_quadratic != observed_quadratic or expected_rational != observed_rational:
        raise ArithmeticError(f"exact closure operands do not replay at p={prime}")
    modular_replays.append({
        "prime": prime,
        "path": str(path.relative_to(ROOT)),
        "sha256": sha256(path),
        "quadratic_operand_replayed": True,
        "rational_operand_replayed": True,
    })


def rational_record(value):
    value = QQ(value)
    return {"numerator": str(value.numerator()), "denominator": str(value.denominator())}


def bit_size(value):
    value = QQ(value)
    return max(abs(value.numerator()).nbits(), value.denominator().nbits())


output = {
    "schema": "elkies-k3.q80-third-q12-closure-operands-p19-hensel-qq.v1",
    "status": "PASS_EXACT_QQ_THIRD_Q12_QUADRATIC_CLOSURE_OPERANDS_P19_HENSEL",
    "specialization": {"u": "-2"},
    "quadratic_field": {
        "generator": "l",
        "polynomial": f"z^2 - ({q_value})",
        "l_squared": rational_record(q_value),
        "nonsquare_over_QQ": True,
    },
    "quadratic_operand": {
        "l": "l",
        "x_coefficients_low_to_high": [
            *[rational_record(value) for value in quadratic_x],
            rational_record(q_value),
        ],
        "sat": f"l/({q_value})",
        "trace_l": "0",
        "exact_closure_substitution": True,
    },
    "rational_operand": {
        "l": rational_record(rational_l),
        "x_coefficients_low_to_high": [
            *[rational_record(value) for value in rational_x],
            rational_record(rational_l**2),
        ],
        "sat": rational_record(1 / rational_l),
        "exact_closure_substitution": True,
    },
    "hensel": {
        "prime": 19,
        "precision_exponent": args.precision,
        "modulus_bits": int(hensel_modulus.nbits()),
        "quadratic_branch_jacobian_rank": quadratic_rank,
        "quadratic_branch_rows_zero_based": list(quadratic_rows),
        "rational_branch_jacobian_rank": rational_rank,
        "rational_branch_rows_zero_based": list(rational_rows),
        "quadratic_iterations": quadratic_iterations,
        "rational_iterations": rational_iterations,
        "rational_projective_reconstruction": {
            "candidate_count": len(rational_projective_candidates),
            "exact_candidate_count": len(valid_rational_candidates),
            "projective_scale": str(rational_projective["projective_scale"]),
            "primitive_vector_max_bits": rational_projective["primitive_vector_max_bits"],
            "lll_norm_squared": str(rational_projective["lll_norm_squared"]),
        },
    },
    "coefficient_height": {
        "quadratic_max_bits": max(map(bit_size, quadratic_values)),
        "rational_max_bits": max(map(bit_size, rational_values)),
    },
    "checks": {
        "exact_six_equation_quadratic_branch": True,
        "exact_six_equation_rational_branch": True,
        "modular_operand_replay_primes": [record["prime"] for record in modular_replays],
        "modular_operand_replays": modular_replays,
    },
    "inputs": {
        "surface": {"path": str(args.surface.relative_to(ROOT)), "sha256": sha256(args.surface)},
        "p19_horizontal": {
            "path": str(args.horizontal.relative_to(ROOT)),
            "sha256": sha256(args.horizontal),
        },
    },
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            "exact quadratic number field generated by the trace-zero polynomial-section coefficient l",
            "exact rational and quadratic closure operands at u=-2",
            "literal exact substitution into all six characteristic-zero closure equations",
            "literal reduction replay at every listed unique-orbit reconstruction prime",
        ],
        "not_proved": [
            "the exact sum as a target horizontal section over the quadratic field",
            "a characteristic-zero resolved pencil, child Jacobian, or birational maps",
            "the characteristic-zero A5+A3+3A1 marking or Mordell--Weil rank",
        ],
    },
    "runtime_seconds": time.monotonic() - started,
    "reproduce": (
        "sage -python elkies-k3/scripts/lift_q80_third_q12_closure_operands_p19_qq.sage "
        f"--precision {args.precision} --output {args.output}"
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
print(
    f"Q80THIRDQ12CLOSUREHENSEL|precision={args.precision}|bits={hensel_modulus.nbits()}|"
    f"quadratic_bits={output['coefficient_height']['quadratic_max_bits']}|"
    f"rational_bits={output['coefficient_height']['rational_max_bits']}|"
    f"replay_primes={','.join(map(str, output['checks']['modular_operand_replay_primes']))}|"
    "status=PASS_EXACT_QQ_THIRD_Q12_QUADRATIC_CLOSURE_OPERANDS_P19_HENSEL",
    flush=True,
)
